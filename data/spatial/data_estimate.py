#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
data_estimate.py -- exploratory analysis and first-pass moments for the spatial
teachers model (`julia/spatial_model/spatial_continuous.jl`).

Consumes the tables `fetch_data.py` left in `raw/` and produces the spatial
moments that discipline the model's L = 2 locations.  Following
`notes/spatial_calibration_notes.md`, a "location" is a GROUP OF SCHOOL
DISTRICTS INSIDE A COMMUTING ZONE, so every headline number here is a
*within-CZ* gap, enrollment-weighted across CZs.  Locations 1 and 2 are ordered
so that 2 is the advantaged group (suburb / high-revenue / high-SES), matching
the model's `κ = [0.75, 0.9]`, `B = [0.0, 0.1]`.

What it does
    1. inventory  -- report which raw tables exist and which are still missing
    2. panel      -- district x year: enrollment, FTE teachers, F-33 revenue and
                     salaries, SAIPE poverty, CWIFT deflator, SEDA learning
                     rates, county -> commuting zone
    3. partition  -- split each CZ in two, under three interchangeable schemes
    4. moments    -- within-CZ gaps in the objects the model needs
                     (κ_2/κ_1, Q_2/Q_1, H̃_T/M, M_2/M_1, t_l, revenue shares),
                     plus the share of each gap's variance that is within-CZ.
                     Location totals are SUMMED and per-pupil objects are ratios
                     of totals; an average of district ratios would answer a
                     question the model does not ask.
    5. regressions-- σ (class-size congestion), σ_ν (mobility), the κ gradient,
                     and the local-revenue -> teacher-pay pass-through
    6. anchors    -- literature values for the parameters no public district
                     table can identify (λ, ψ, m, ρ_z), carried as placeholders
    7. emit       -- one JSON + one markdown report, each estimate tagged with
                     its source and a one-line note on how it disciplines the model

Every block degrades gracefully: a missing input becomes a `missing-input`
estimate carrying the source that would supply it, so the report doubles as a
to-do list for the empirical section.

Run against the full download of 2026-08-27: CCD, F-33, SAIPE, SEDA 6.0, CWIFT
and Dorn's crosswalk are all exercised against the real tables.  ACS is the one
source still absent -- it needs a Census key -- so the effective tax rate t_l is
the only moment here with no data behind it.  Three things the raw files force,
and the reason each resolver looks the way it does:

  * SEDA estimate files are LONG IN SUBGROUP.  Race, gender and
    economic-disadvantage rows sit under the same district, and rows with
    `gap == 1` are differences between two subgroups, not levels.  Only the
    all-students level rows are read.
  * CWIFT stacks district, county and state universes in one table with `LEAID`
    blank on the last two, and is published for the 2013-14 universe onward, so
    it is merged on the nearest release within CWIFT_YEAR_TOLERANCE years and
    is simply absent from the earlier cross-sections.
  * F-33 gives some agencies an id with a letter in it (`20D0001`), which is a
    different agency from `0200001`; district keys are normalised but never
    stripped of letters, and every merge is de-duplicated first.

Usage
    uv run data_estimate.py                        # base year 2018, locale split
    uv run data_estimate.py --scheme all           # all three CZ partitions
    uv run data_estimate.py --base-year 2010 --scheme revenue
    uv run data_estimate.py --years 1994 2000 2010 2018 --dump-panel
    uv run data_estimate.py --geo cbsa             # metro instead of commuting zone

Needs pandas, numpy, pyarrow.
Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import re
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

LOG = logging.getLogger("data_estimate")
HERE = Path(__file__).resolve().parent


# =============================================================================
# Constants: sample definition and the model objects each moment maps onto
# =============================================================================

# CCD agency types kept: 1 = regular local district, 2 = component of a
# supervisory union.  Charters (7), regional service agencies (4-6) and
# state/federally operated agencies (8-9) are dropped -- none of them is the
# taxing, hiring unit the model's location is.
REGULAR_AGENCY_TYPES = (1.0, 2.0)

# NCES urban-centric locale: 11-13 city, 21-23 suburb, 31-33 town, 41-43 rural.
LOCALE_GROUPS = {1: "city", 2: "suburb", 3: "town", 4: "rural"}

# F-33 teacher-salary components; summed to get the teaching wage bill.
TEACHER_SALARY_PARTS = ("salaries_teachers_regular_prog", "salaries_teachers_sped",
                        "salaries_teachers_vocational", "salaries_teachers_other_ed")

# Urban Institute / F-33 sentinel codes for missing, suppressed, N/A.
SENTINELS = (-1, -2, -3, -9, -99)

# Test-score -> log-earnings conversion used to put SEDA learning rates in the
# model's units (log Q enters log h one-for-one, and log h is log earnings).
# 1 SD of test score ~ 0.13 log points of adult earnings; Chetty-Friedman-Rockoff
# (2014b) get ~0.12 from teacher VA, Hanushek's survey brackets 0.10-0.20.
# EVERY estimate that uses it is flagged, because it is an assumption, not data.
SD_SCORE_TO_LOG_EARNINGS = 0.13

# SEDA publishes every estimate on four scales.  `cs` (cohort-standardised) is
# the one whose units are SDs of the national grade-cohort distribution, which
# is what SD_SCORE_TO_LOG_EARNINGS converts; `gcs`/`gys` are grade equivalents
# (a mean of ~5.5 at the grade-5.5 center) and would need a different one.
# `ol` is the unshrunken estimate: the empirical-Bayes version pulls small
# districts toward the national mean, which is exactly the wrong thing to do
# before an enrollment-weighted comparison of two groups of districts.
SEDA_SCALE = "cs"
SEDA_ESTIMATOR = "ol"
# SEDA fits the learning rate over grades 3-8, so a district's rate -- SDs of
# achievement gained per grade -- cumulates over five grades into the level
# difference that Q_l is.  The regressions below work in that cumulative unit;
# the reported gap stays in SEDA's published per-grade unit.
SEDA_GRADE_SPAN = 5.0

# CWIFT starts with the 2013-14 district universe, so an exact-year merge would
# leave every earlier cross-section undeflated.  The index moves slowly, so the
# nearest release within this many years is used -- and no further, which is why
# the deflated salary gap simply does not exist before 2013.
CWIFT_YEAR_TOLERANCE = 2

# Benchmarks the model is not fitted to but should reproduce (see notes §6).
LIT_SIGMA_NU = 0.20          # Eckert-Kleineberg spatial-taste dispersion
LIT_EXPOSURE_RATE = 0.04     # Chetty-Hendren convergence per year of exposure
LIT_LOCAL_REV_SHARE = 0.44   # NCES, FY2020-21 national local share


# =============================================================================
# The estimate registry -- one record per number that leaves this script
# =============================================================================

@dataclass
class Estimate:
    key: str                       # short identifier, stable across runs
    label: str                     # human description of the quantity
    model_object: str              # the parameter/moment in spatial_continuous.jl
    value: float | None
    units: str
    source: str                    # dataset, years, unit of observation
    note: str                      # what we learn / how it disciplines the model
    status: str = "ok"             # ok | missing-input | assumption | literature
    se: float | None = None
    n: int | None = None
    extra: dict[str, Any] = field(default_factory=dict)


class Registry:
    """Accumulates `Estimate`s and writes the JSON and markdown reports."""

    def __init__(self) -> None:
        self.items: list[Estimate] = []

    def add(self, est: Estimate) -> Estimate:
        self.items.append(est)
        v = "n/a" if est.value is None else f"{est.value:.4g}"
        LOG.info("  [%-12s] %-28s = %s", est.status, est.key, v)
        return est

    def missing(self, key: str, label: str, model_object: str, source: str,
                note: str, units: str = "") -> Estimate:
        return self.add(Estimate(key, label, model_object, None, units, source,
                                 note, status="missing-input"))

    def to_dict(self) -> list[dict[str, Any]]:
        return [asdict(e) for e in self.items]


# =============================================================================
# Raw-table access and column resolution
# =============================================================================

class RawStore:
    """Reads whatever `fetch_data.py` managed to write, by stem name.

    Column selection happens before the read, so a 50 MB CCD panel costs only
    the handful of columns actually wanted.
    """

    def __init__(self, raw_dir: Path) -> None:
        self.dir = raw_dir
        self.paths: dict[str, Path] = {}
        for ext in (".parquet", ".csv"):
            for p in sorted(raw_dir.glob(f"*{ext}")):
                self.paths.setdefault(p.stem, p)

    def stems(self, pattern: str) -> list[str]:
        rx = re.compile(pattern, re.I)
        return sorted(s for s in self.paths if rx.search(s))

    def columns(self, stem: str) -> list[str]:
        path = self.paths[stem]
        if path.suffix == ".parquet":
            import pyarrow.parquet as pq
            return list(pq.ParquetFile(path).schema_arrow.names)
        return list(pd.read_csv(path, nrows=0).columns)

    def load(self, stem: str, keep: Sequence[str] | None = None
             ) -> pd.DataFrame | None:
        """Load `stem`, keeping columns whose names fullmatch any regex in `keep`."""
        if stem not in self.paths:
            LOG.warning("  missing raw table: %s", stem)
            return None
        path, cols = self.paths[stem], None
        if keep:
            avail = self.columns(stem)
            cols = [c for c in avail if any(re.fullmatch(k, c, re.I) for k in keep)]
            if not cols:
                LOG.warning("  %s: none of the requested columns present", stem)
                return None
        df = (pd.read_parquet(path, columns=cols) if path.suffix == ".parquet"
              else pd.read_csv(path, usecols=cols, low_memory=False))
        LOG.info("  loaded %-28s %7d rows x %d cols", stem, len(df), df.shape[1])
        return df

    def inventory(self, expected: Iterable[str]) -> dict[str, dict[str, Any]]:
        """One row per expected source.  SEDA ships ~100 tables; listing each of
        them would bury the three that are missing, so a prefix matching more
        than one file is reported as a count plus the tables actually read."""
        out: dict[str, dict[str, Any]] = {}
        for stem in expected:
            hits = self.stems(f"^{stem}")
            if not hits:
                out[stem] = {"status": "missing"}
            elif len(hits) == 1:
                h = hits[0]
                out[h] = {"status": "ok", "path": self.paths[h].name,
                          "n_cols": len(self.columns(h))}
            else:
                out[f"{stem}*"] = {"status": "ok", "n_tables": len(hits),
                                   "tables": hits[:4] + (["..."] if len(hits) > 4
                                                         else [])}
        return out


def pick(df: pd.DataFrame, *patterns: str) -> str | None:
    """First column whose name fullmatches one of the patterns, in priority order."""
    for pat in patterns:
        for c in df.columns:
            if re.fullmatch(pat, c, re.I):
                return c
    return None


def id7(s: pd.Series) -> pd.Series:
    """NCES LEAID as a 7-character string, whatever dtype it arrived in.

    CWIFT stores it as a float (`100005.0`), SEDA as an int, the CCD as text; a
    bare `astype(str).zfill(7)` silently produces three different keys and the
    merge quietly returns nothing.  Letters are KEPT -- F-33 gives Alaska's
    supervisory agencies ids like `20D0001`, and stripping the `D` would collide
    them onto district `0200001` and double-count it -- while the Urban `-1`
    sentinel and a stringified NaN are dropped.
    """
    x = (s.astype(str).str.strip().str.upper()
         .str.replace(r"\.0+$", "", regex=True))
    return x.where(x.str.fullmatch(r"(?=.*\d)[0-9A-Z]{1,7}")).str.zfill(7)


def num(s: pd.Series, positive: bool = False) -> pd.Series:
    """Coerce to float, blank out the sentinel codes, optionally require > 0."""
    x = pd.to_numeric(s, errors="coerce")
    x = x.mask(x.isin(SENTINELS))
    return x.mask(x <= 0) if positive else x


def rowsum(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    """Sum the columns that exist; NaN where every component is missing."""
    have = [c for c in cols if c in df.columns]
    if not have:
        return pd.Series(np.nan, index=df.index)
    block = df[have].apply(num)
    return block.sum(axis=1, min_count=1)


# =============================================================================
# Weighted statistics and a small fixed-effects OLS
# =============================================================================

def wmean(x: pd.Series, w: pd.Series) -> float:
    m = x.notna() & w.notna() & (w > 0)
    return float(np.average(x[m], weights=w[m])) if m.any() else np.nan


def wquantile(x: pd.Series, w: pd.Series, q: float) -> float:
    m = x.notna() & w.notna() & (w > 0)
    if not m.any():
        return np.nan
    xs, ws = x[m].to_numpy(float), w[m].to_numpy(float)
    o = np.argsort(xs)
    xs, ws = xs[o], ws[o]
    c = np.cumsum(ws) - 0.5 * ws
    return float(np.interp(q, c / ws.sum(), xs))


def _wdemean(M: np.ndarray, g: pd.Series, w: np.ndarray) -> np.ndarray:
    """Weighted within-transform: subtract each column's w-weighted group mean."""
    idx = pd.RangeIndex(len(w))
    gs = pd.Series(np.asarray(g), index=idx)
    den = pd.Series(w, index=idx).groupby(gs).transform("sum").to_numpy()
    out = np.empty_like(M, dtype=float)
    for j in range(M.shape[1]):
        numj = pd.Series(M[:, j] * w, index=idx).groupby(gs).transform("sum").to_numpy()
        out[:, j] = M[:, j] - numj / den
    return out


def ols_fe(df: pd.DataFrame, y: str, x: Sequence[str], fe: str | None = None,
           w: str | None = None, cluster: str | None = None) -> dict[str, Any] | None:
    """WLS with optionally absorbed fixed effects and cluster-robust SEs.

    Returns None when too little survives listwise deletion.  Written out rather
    than pulled from statsmodels to keep this script's dependency set equal to
    `fetch_data.py`'s.
    """
    cols = list(dict.fromkeys([y, *x] + [c for c in (fe, w, cluster) if c]))
    d = df[[c for c in cols if c in df.columns]].copy()
    if len(d.columns) < len(cols):
        return None
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 20 + len(x):
        return None

    Y = d[[y]].to_numpy(float)
    X = d[list(x)].to_numpy(float)
    wt = d[w].to_numpy(float) if w else np.ones(len(d))
    names = list(x)
    if fe:
        Y, X = _wdemean(Y, d[fe], wt), _wdemean(X, d[fe], wt)
    else:
        X, names = np.column_stack([np.ones(len(d)), X]), ["const", *x]

    sw = np.sqrt(wt)[:, None]
    b, *_ = np.linalg.lstsq(X * sw, Y * sw, rcond=None)
    u = (Y - X @ b).ravel()
    bread = np.linalg.pinv((X * wt[:, None]).T @ X)
    S = X * (wt * u)[:, None]
    if cluster:
        grp = d[cluster].to_numpy()
        Sc = pd.DataFrame(S).groupby(grp).sum().to_numpy()
        G, n, k = len(Sc), len(d), X.shape[1]
        meat = Sc.T @ Sc * (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
    else:
        G, meat = None, S.T @ S
    V = bread @ meat @ bread
    tss = float(np.sum(wt * (Y.ravel() - np.average(Y.ravel(), weights=wt)) ** 2))
    return {"coef": dict(zip(names, b.ravel())),
            "se": dict(zip(names, np.sqrt(np.diag(V)))),
            "n": len(d), "n_clusters": G,
            "r2_within": 1 - float(np.sum(wt * u ** 2)) / tss if tss > 0 else np.nan}


def precision(b: float, se: float | None) -> str:
    """A one-clause verdict on whether a slope is distinguishable from zero.

    Several of these regressions come back with the wrong sign and a standard
    error twice the coefficient. That is a finding -- the within-CZ cross-section
    does not identify the parameter -- and it belongs in the note rather than in
    a reader's arithmetic.
    """
    if se is None or not np.isfinite(se) or se <= 0:
        return "no standard error"
    t = b / se
    if abs(t) < 1.96:
        return (f"NOT distinguishable from zero (t = {t:+.2f}), so this is an "
                f"upper bound on how much the cross-section knows")
    return f"significant (t = {t:+.2f})"


def fmt_reg(r: dict[str, Any] | None, key: str) -> str:
    if not r or key not in r["coef"]:
        return "n/a"
    return f"{r['coef'][key]:+.4f} ({r['se'][key]:.4f})"


# =============================================================================
# Builders -- one per raw source, each returning a district-keyed frame
# =============================================================================

def build_ccd(store: RawStore, years: Sequence[int]) -> pd.DataFrame | None:
    """CCD directory: the sample frame, enrollment, FTE teachers, locale, county."""
    df = store.load("urban_ccd_directory",
                    ["year", "leaid", "fips", "lea_name", "county_code",
                     "urban_centric_locale", "enrollment", "teachers_total_fte",
                     "agency_type", "agency_level", "number_of_schools", "cbsa"])
    if df is None:
        return None
    df = df[df.year.isin(years)].copy()
    df = df[df.agency_type.isin(REGULAR_AGENCY_TYPES)]
    df["leaid"] = id7(df.leaid)
    df["county_fips"] = (df.county_code.astype(str).str.extract(r"(\d+)")[0]
                         .str.zfill(5))
    df["pupils"] = num(df.enrollment, positive=True)
    df["teachers"] = num(df.teachers_total_fte, positive=True)
    df["str_ratio"] = df.pupils / df.teachers
    # Locale groups only exist in the 2-digit urban-centric coding (2006 on);
    # earlier vintages use a different 1-8 scheme, left as NaN on purpose.
    loc = num(df.urban_centric_locale)
    df["locale_grp"] = (loc // 10).where(loc >= 11).map(LOCALE_GROUPS)
    return df.drop(columns=["enrollment", "teachers_total_fte", "county_code"])


def build_finance(store: RawStore, years: Sequence[int]) -> pd.DataFrame | None:
    """F-33: revenue by source, current expenditure, the teacher wage bill."""
    df = store.load("urban_ccd_finance",
                    ["year", "leaid", "rev_total", "rev_local_total",
                     "rev_state_total", "rev_fed_total", "rev_local_prop_tax",
                     "exp_total", "exp_current_elsec_total",
                     "exp_current_instruction_total", "salaries_total",
                     "salaries_instruction", "benefits_employee_total",
                     "enrollment_fall_responsible", *TEACHER_SALARY_PARTS])
    if df is None:
        return None
    df = df[df.year.isin(years)].copy()
    df["leaid"] = id7(df.leaid)
    df["wagebill_teachers"] = rowsum(df, TEACHER_SALARY_PARTS)
    # Pre-2000 F-33 vintages report only the instruction aggregate; fall back to
    # it so the early cross-sections are not silently empty.
    if "salaries_instruction" in df:
        df["wagebill_teachers"] = df.wagebill_teachers.fillna(
            num(df.salaries_instruction))
    for c in ("rev_total", "rev_local_total", "rev_state_total", "rev_fed_total",
              "rev_local_prop_tax", "exp_current_elsec_total",
              "exp_current_instruction_total", "salaries_total",
              "benefits_employee_total", "enrollment_fall_responsible"):
        if c in df:
            df[c] = num(df[c], positive=c.startswith(("rev_total", "exp_", "enrol")))
    return df.drop(columns=[c for c in ("salaries_instruction",
                                        *TEACHER_SALARY_PARTS) if c in df])


def build_saipe(store: RawStore, years: Sequence[int]) -> pd.DataFrame | None:
    """SAIPE: child poverty and the district population base."""
    df = store.load("urban_saipe",
                    ["year", "leaid", "est_population_total", "est_population_5_17",
                     "est_population_5_17_poverty"])
    if df is None:
        return None
    df = df[df.year.isin(years)].copy()
    df["leaid"] = id7(df.leaid)
    df["pop_5_17"] = num(df.est_population_5_17, positive=True)
    df["pop_5_17_poverty"] = num(df.est_population_5_17_poverty)
    df["pov_rate_5_17"] = df.pop_5_17_poverty / df.pop_5_17
    df["pop_total"] = num(df.est_population_total, positive=True)
    return df[["year", "leaid", "pop_total", "pop_5_17", "pop_5_17_poverty",
               "pov_rate_5_17"]]


def build_cwift(store: RawStore) -> pd.DataFrame | None:
    """CWIFT: the comparable-wage deflator that strips local price levels out of κ.

    The EDGE table stacks three universes in one file -- one block of district
    rows, one of county rows and one of state rows per release -- with `LEAID`
    blank on the last two, so the district index has to be selected explicitly
    (`LEA_CWIFTEST`, not the `CNTY_`/`ST_` columns sitting beside it).  The
    school year of each district block lives in `source_file`
    (`EDGE_ACS_CWIFT2018_LEA1819` -> 2018 on the CCD's fall-year convention),
    which is the year the index actually describes; `cwift_year` is the ACS
    vintage and runs one to two years ahead of it.
    """
    df = store.load("edge_cwift")
    if df is None:
        return None
    lea = pick(df, r"leaid", r"lea_id", r"ncesid")
    est = pick(df, r"lea_cwiftest", r"lea_cwift", r"cwiftest", r".*lea.*cwift.*est.*",
               r".*cwift.*est.*", r"cwift")
    if not lea or not est:
        LOG.warning("  edge_cwift: could not resolve leaid/index columns from %s",
                    list(df.columns)[:12])
        return None
    out = pd.DataFrame({"leaid": id7(df[lea]),
                        "cwift": num(df[est], positive=True),
                        "year": _cwift_year(df)})
    out = out.dropna(subset=["leaid", "cwift"])
    if out.empty:
        LOG.warning("  edge_cwift: %s resolved but every district row is blank", est)
        return None
    LOG.info("  cwift: index=%s, %d district rows, years %s", est, len(out),
             sorted(int(y) for y in out.year.dropna().unique())
             if out.year.notna().any() else "unknown")
    if out.year.notna().any():
        return (out.dropna(subset=["year"])
                .groupby(["leaid", "year"], as_index=False).cwift.mean())
    return out.groupby("leaid", as_index=False).cwift.mean()


def _cwift_year(df: pd.DataFrame) -> pd.Series:
    """Fall year of the district universe each CWIFT block describes."""
    src = pick(df, r"source_file", r"file", r"source")
    if src:
        sy = df[src].astype(str).str.extract(r"LEA(\d{2})\d{2}", flags=re.I)[0]
        yr = pd.to_numeric(sy, errors="coerce")
        yr = yr.where(yr.isna(), yr + np.where(yr < 80, 2000, 1900))
        if yr.notna().any():
            return yr
    col = pick(df, r"cwift_year", r"year")
    return num(df[col]) if col else pd.Series(np.nan, index=df.index)


def _seda_pick_stem(store: RawStore, pattern: str, exclude: str | None = None
                    ) -> str | None:
    """The pooled district file on the configured scale, if it is there.

    SEDA ships the same estimates at five aggregations x four scales, so the
    choice has to be made on the file NAME rather than by taking whichever
    string sorts first.
    """
    cands = [s for s in store.stems(pattern)
             if exclude is None or not re.search(exclude, s)]
    if not cands:
        return None
    def rank(stem: str) -> tuple:
        parts = stem.split("_")
        return (SEDA_SCALE not in parts and any(
                    sc in parts for sc in ("cs", "gcs", "ys", "gys")),
                "pool" not in parts, "long" in parts, len(stem))
    return sorted(cands, key=rank)[0]


def _seda_all_students(df: pd.DataFrame, stem: str) -> pd.DataFrame:
    """Keep the all-students LEVEL rows.

    Every SEDA estimate file is long in subgroup: race, gender and
    economic-disadvantage rows sit under the same district, and rows with
    `gap == 1` are DIFFERENCES between two subgroups, not levels.  Averaging
    across them -- which is what happens if the long dimension is ignored --
    mixes a white-black test-score gap into a district mean and shrinks it
    toward zero.
    """
    n0 = len(df)
    for col, keep in (("subcat", "all"), ("subgroup", "all")):
        if col in df.columns:
            df = df[df[col].astype(str).str.lower().eq(keep)]
    if "gap" in df.columns:
        df = df[num(df["gap"]).fillna(0).eq(0)]
    if len(df) < n0:
        LOG.info("  seda: %s -> %d of %d rows are all-students levels",
                 stem, len(df), n0)
    return df


def build_seda(store: RawStore) -> pd.DataFrame | None:
    """SEDA district files: mean score, LEARNING RATE (the Q_l analogue), and SES.

    Two files are read.  The estimate file carries achievement, and learning
    rates -- grade-to-grade growth -- are the right Q analogue because they net
    out the SES level that score *levels* mostly measure.  SES itself lives only
    in the separate district COVARIATE file, so `seda_ses` comes from there.
    """
    stem = _seda_pick_stem(store, r"^seda_.*(geodist|dist)", exclude=r"cov|crosswalk")
    if stem is None:
        LOG.warning("  no SEDA district estimate table in raw/")
        return None
    df = store.load(stem)
    if df is None:
        return None
    lea = pick(df, r"sedalea", r"leaid", r"leaidC", r"sedaleaid")
    if not lea:
        LOG.warning("  %s: no district id column found", stem)
        return None
    df = _seda_all_students(df, stem)
    sc, ev = SEDA_SCALE, SEDA_ESTIMATOR
    # 6.0 names the learning rate `<scale>_mn_lrn_<estimator>`; 4.x/5.x called it
    # `mn_grd`.  Both are tried before giving up, newest first.
    grd = pick(df, rf"{sc}_mn_lrn_{ev}", rf"{sc}_mn_lrn_eb", r".*_mn_lrn_(ol|eb)",
               rf"{sc}_mn_grd_all", r".*mn_grd.*all.*")
    lvl = pick(df, rf"{sc}_mn_avg_{ev}", rf"{sc}_mn_avg_eb", r".*_mn_avg_(ol|eb)",
               rf"{sc}_mn_all", r".*mn_avg.*")
    wgt = pick(df, r"tot_asmts", r"totgyb_all", r"n_all", r"totenrl")
    if not grd:
        LOG.warning("  %s: no learning-rate column; Q_2/Q_1 and σ will be blank "
                    "(saw %s)", stem, [c for c in df.columns if "mn_" in c][:8])
    LOG.info("  seda: %s -> id=%s learn=%s level=%s weight=%s", stem, lea, grd,
             lvl, wgt)
    out = pd.DataFrame({"leaid": id7(df[lea])})
    out["learn_rate"] = num(df[grd]) if grd else np.nan
    out["score_mean"] = num(df[lvl]) if lvl else np.nan
    out["seda_w"] = num(df[wgt], positive=True) if wgt else 1.0
    out = out.dropna(subset=["leaid"])
    # Collapse any subject / grade / year dimension to one row per district.
    # Vectorised weighted mean: sum(w*x)/sum(w over the rows where x is present).
    vars_ = ["learn_rate", "score_mean"]
    w = out.seda_w.fillna(1.0)
    num_ = pd.DataFrame({v: out[v] * w for v in vars_})
    den_ = pd.DataFrame({v: w.where(out[v].notna()) for v in vars_})
    key = out.leaid
    agg = (num_.groupby(key).sum(min_count=1) /
           den_.groupby(key).sum(min_count=1)).reset_index()
    cov = build_seda_covariates(store)
    if cov is not None:
        agg = agg.merge(cov, on="leaid", how="outer")
    return agg


def build_seda_covariates(store: RawStore) -> pd.DataFrame | None:
    """SEDA district covariates: the composite SES index and median income.

    `sesavgall` is the ACS-based composite the `ses` partition is meant to split
    on, and `lninc50avgall` is log median household income -- the closest public
    counterpart to the income gap the model's sorting produces.
    """
    stem = _seda_pick_stem(store, r"^seda_cov_.*(geodist|admindist)")
    if stem is None:
        LOG.warning("  no SEDA district covariate table in raw/; the ses split "
                    "will fall back to SAIPE poverty")
        return None
    df = store.load(stem)
    if df is None:
        return None
    lea = pick(df, r"sedalea", r"leaid", r"sedaleaid")
    ses = pick(df, r"sesavgall", r"sesall", r"ses", r".*ses.*all.*")
    inc = pick(df, r"lninc50avgall", r"lninc50all", r".*lninc.*all.*")
    if not lea or not ses:
        LOG.warning("  %s: could not resolve leaid/SES from %s", stem,
                    list(df.columns)[:12])
        return None
    df = _seda_all_students(df, stem)
    out = pd.DataFrame({"leaid": id7(df[lea]), "seda_ses": num(df[ses])})
    if inc:
        out["seda_lninc"] = num(df[inc])
    LOG.info("  seda cov: %s -> ses=%s income=%s", stem, ses, inc)
    return out.dropna(subset=["leaid"]).groupby("leaid", as_index=False).mean()


def build_czones(store: RawStore) -> pd.DataFrame | None:
    """Dorn's county -> 1990 commuting-zone crosswalk: the labor market boundary."""
    df = store.load("czone_crosswalk")
    if df is None:
        return None
    cty = pick(df, r"cty_fips", r"county_fips", r".*cty.*", r".*county.*", r"fips")
    cz = pick(df, r"czone", r"cz", r".*czone.*", r".*cz\d*")
    if not cty or not cz:
        LOG.warning("  czone_crosswalk: unresolved columns %s", list(df.columns)[:10])
        return None
    out = pd.DataFrame({
        "county_fips": df[cty].astype(str).str.extract(r"(\d+)")[0].str.zfill(5),
        "czone": num(df[cz])})
    return out.dropna().drop_duplicates("county_fips")


def build_acs(store: RawStore, years: Sequence[int]) -> pd.DataFrame | None:
    """ACS district tables: the income base for the effective local tax rate t_l."""
    df = store.load("acs_school_districts")
    if df is None:
        return None
    # The API names the geography column after the layer, so the stacked table
    # carries three of them -- "school district (unified)", "(elementary)",
    # "(secondary)" -- each filled only on its own rows.  Taking the first would
    # key the unified districts and silently drop the other two layers.
    sd_cols = [c for c in df.columns
               if re.fullmatch(r"school district.*|sd_.*code|district", c, re.I)]
    st = pick(df, r"state")
    if not sd_cols or not st:
        LOG.warning("  acs: could not build leaid from %s", list(df.columns)[:12])
        return None
    code = df[sd_cols[0]].astype(str)
    for c in sd_cols[1:]:
        code = code.where(df[c].isna(), df[c].astype(str))
    out = pd.DataFrame({
        "leaid": (df[st].astype(str).str.zfill(2)
                  + code.str.extract(r"(\d+)")[0].str.zfill(5)),
        "year": num(df["year"]) if "year" in df else np.nan})
    for c in ("median_hh_income", "per_capita_income", "median_home_value",
              "median_property_tax", "population"):
        if c in df:
            out[c] = num(df[c], positive=True)
    out = out[out.year.isin(years)] if out.year.notna().any() else out
    return out.drop_duplicates(["leaid", "year"])


# =============================================================================
# Assembly and the derived district-level variables
# =============================================================================

# Derived ratios carry impossible tails -- F-33 wage bills matched to a CCD
# teacher count from a different reporting universe give $2M per FTE, class sizes
# of 1000, and per-pupil spending of $860k.  These are reporting errors, not
# small districts, so they are blanked rather than winsorised.  The cut is on
# WITHIN-YEAR quantiles because the panel is nominal and spans 1994-2018.
SCREEN_VARS = ("str_ratio", "teachers_pp", "salary_per_teacher", "salary_real",
               "exp_pp", "rev_pp_local", "rev_pp_state", "rev_pp_total",
               "proptax_pp", "t_eff")


def screen_outliers(p: pd.DataFrame, tail: float = 0.005) -> pd.DataFrame:
    """Blank each derived ratio outside its within-year [tail, 1-tail] range."""
    for c in SCREEN_VARS:
        if c not in p:
            continue
        lo = p.groupby("year")[c].transform(lambda x: x.quantile(tail))
        hi = p.groupby("year")[c].transform(lambda x: x.quantile(1 - tail))
        bad = (p[c] < lo) | (p[c] > hi)
        if bad.any():
            LOG.info("  screened %5d/%d values of %s", int(bad.sum()),
                     int(p[c].notna().sum()), c)
        p[c] = p[c].mask(bad)
    return p


def dedup(df: pd.DataFrame, keys: Sequence[str], name: str) -> pd.DataFrame:
    """Guarantee one row per key before a left merge.

    A right-hand table with a repeated key silently multiplies panel rows, and
    the location totals the moments are built from would then double-count the
    affected districts.  Duplicates are dropped rather than tolerated, loudly.
    """
    d = df.dropna(subset=list(keys))
    dup = int(d.duplicated(list(keys)).sum())
    if dup:
        LOG.warning("  %s: %d duplicate %s rows dropped before the merge",
                    name, dup, "/".join(keys))
        d = d.drop_duplicates(list(keys))
    return d


def merge_nearest_year(panel: pd.DataFrame, other: pd.DataFrame, value: str,
                       tolerance: int) -> pd.DataFrame:
    """Merge a district x year table on the nearest vintage within `tolerance`.

    CWIFT is published for a handful of school years, none of which has to be a
    year the panel is built for.  Merging on the exact year would silently drop
    the deflator everywhere; extrapolating it back two decades would invent one.
    The nearest release within `tolerance` years is used, the year actually used
    is carried in `<value>_year`, and everything further away stays missing.
    """
    if "year" not in other.columns or other.year.isna().all():
        return panel.merge(other, on="leaid", how="left")
    # merge_asof refuses to join an int year to a float one, and the CCD's is
    # int while a parsed vintage is float, so both are cast before the merge.
    as_float = lambda d: pd.to_numeric(d.year, errors="coerce").astype(float)
    left = (panel[["leaid", "year"]].reset_index()
            .assign(year=as_float).dropna(subset=["year"]).sort_values("year"))
    right = (other.dropna(subset=["leaid", "year", value])
             .assign(year=as_float).sort_values("year"))
    right[f"{value}_year"] = right.year
    m = pd.merge_asof(left, right, on="year", by="leaid", direction="nearest",
                      tolerance=tolerance).set_index("index")
    panel[value] = m[value]
    panel[f"{value}_year"] = m[f"{value}_year"]
    hit = panel.groupby("year")[value].apply(lambda x: x.notna().mean())
    LOG.info("  %s merged within %d years: coverage %s", value, tolerance,
             {int(k): round(float(v), 2) for k, v in hit.items()})
    return panel


def assemble(store: RawStore, years: Sequence[int], geo: str = "auto"
             ) -> pd.DataFrame:
    """Merge every available source onto the CCD frame and derive model objects.

    `geo` picks the labor-market unit the two locations sit inside: "czone" is
    what the calibration notes argue for, "cbsa" is the CCD's own metro code and
    the fallback when Dorn's crosswalk has not been downloaded yet.  A CBSA is a
    defensible labor market -- it is just not the one the notes settled on, and
    it drops non-metro districts entirely, so the choice is recorded in the
    report rather than made silently.
    """
    ccd = build_ccd(store, years)
    if ccd is None:
        raise SystemExit("urban_ccd_directory is required; run fetch_data.py first")
    panel = ccd

    for name, builder, on in (
            ("finance", build_finance, ["year", "leaid"]),
            ("saipe", build_saipe, ["year", "leaid"])):
        part = builder(store, years)
        if part is not None:
            panel = panel.merge(dedup(part, on, name), on=on, how="left",
                                suffixes=("", f"_{name}"))

    cw = build_cwift(store)
    if cw is not None:
        panel = merge_nearest_year(panel, cw, "cwift", CWIFT_YEAR_TOLERANCE)

    seda = build_seda(store)
    if seda is not None:      # SEDA pools 2009-2019: one cross-section per district
        panel = panel.merge(seda, on="leaid", how="left")

    cz = build_czones(store) if geo in ("auto", "czone") else None
    if cz is not None:
        panel = panel.merge(cz, on="county_fips", how="left")
    geo_unit = "czone (Dorn 1990 commuting zones)"
    if cz is None or panel.czone.isna().all():
        if geo == "czone":
            raise SystemExit("no commuting-zone crosswalk in raw/; rerun "
                             "fetch_data.py --sources czones, or pass --geo cbsa")
        panel["czone"] = panel.cbsa
        if geo == "cbsa":
            geo_unit = "cbsa (requested; non-metro districts dropped)"
            LOG.info("  labor market: the CCD's CBSA code, as requested")
        else:
            geo_unit = "cbsa (FALLBACK -- Dorn crosswalk absent; non-metro dropped)"
            LOG.warning("  no CZ crosswalk: falling back to the CCD's CBSA code")

    acs = build_acs(store, years)
    if acs is not None:
        keys = ["leaid", "year"] if acs.year.notna().any() else ["leaid"]
        panel = panel.merge(acs.drop(columns=[] if "year" in keys else ["year"]),
                            on=keys, how="left")

    # ---- derived: everything the moments below read ------------------------
    p = panel
    if "enrollment_fall_responsible" in p:
        p["pupils"] = p.pupils.fillna(p.enrollment_fall_responsible)
    pupils = p.pupils
    if "wagebill_teachers" in p:
        p["salary_per_teacher"] = p.wagebill_teachers / p.teachers
    if {"cwift", "salary_per_teacher"} <= set(p.columns):
        # κ net of the local price of labor -- CWIFT is what makes the salary gap
        # a teaching-wage gap rather than a cost-of-living gap.
        p["salary_real"] = p.salary_per_teacher / p.cwift
    for src in ("total", "local", "state", "fed"):
        col = f"rev_{src}_total" if src != "total" else "rev_total"
        if col in p:
            p[f"rev_pp_{src}"] = p[col] / pupils
    if {"rev_total", "rev_local_total"} <= set(p.columns):
        for src in ("local", "state", "fed"):
            p[f"share_{src}"] = p[f"rev_{src}_total"] / p.rev_total
    if "exp_current_elsec_total" in p:
        p["exp_pp"] = p.exp_current_elsec_total / pupils
    if "rev_local_prop_tax" in p:
        p["proptax_pp"] = p.rev_local_prop_tax / pupils
    p["teachers_pp"] = p.teachers / pupils      # H̃_T/M, up to the β weighting
    # t_l as the model writes it: local revenue over the local income base.  ACS
    # per-capita income x population is the clean version; the SAIPE fallback
    # assumes 2.5 people per household, coarse enough that t_eff is an order of
    # magnitude rather than a target.
    if {"per_capita_income", "population"} <= set(p.columns):
        p["income_base"] = p.per_capita_income * p.population
    elif {"median_hh_income", "pop_total"} <= set(p.columns):
        p["income_base"] = p.median_hh_income * p.pop_total / 2.5
    if {"income_base", "rev_local_total"} <= set(p.columns):
        p["t_eff"] = p.rev_local_total / p.income_base
    p = screen_outliers(p)
    for c in ("pupils", "teachers", "teachers_pp", "str_ratio", "exp_pp",
              "salary_per_teacher", "salary_real", "rev_pp_local", "rev_pp_state",
              "rev_pp_total", "proptax_pp"):
        if c in p:
            p[f"log_{c}"] = np.log(p[c].where(p[c] > 0))
    LOG.info("panel: %d district-years, %d districts, %d labor markets",
             len(p), p.leaid.nunique(), p.czone.nunique() if "czone" in p else 0)
    p.attrs["geo_unit"] = geo_unit      # set last: a merge would drop it
    return p


# =============================================================================
# The two-location partition inside each commuting zone
# =============================================================================

def partition(panel: pd.DataFrame, scheme: str, min_pupils: int = 500,
              min_group_share: float = 0.10) -> pd.DataFrame:
    """Assign every district a location 1 or 2 inside its CZ.

    Location 2 is the ADVANTAGED group throughout, so a positive log gap always
    means "the good location has more of it".

        locale   1 = city (11-13), 2 = suburb (21-23); town/rural dropped
        revenue  1/2 = below/above the CZ enrollment-weighted median of local
                 revenue per pupil
        ses      1/2 = below/above the CZ median SEDA SES (falls back to the
                 inverse of the SAIPE child-poverty rate)

    A CZ is kept only if both groups clear `min_group_share` of its enrollment;
    that is what makes a two-location representative metro a sensible object.
    """
    d = panel.copy()
    if "czone" not in d or d.czone.isna().all():
        raise SystemExit("no labor-market unit merged; rerun fetch_data.py "
                         "--sources czones")
    d = d[d.czone.notna()]
    d = d[d.pupils.fillna(0) >= min_pupils]

    if scheme == "locale":
        d = d[d.locale_grp.isin(["city", "suburb"])]
        d["loc_id"] = np.where(d.locale_grp.eq("suburb"), 2, 1)
    elif scheme in ("revenue", "ses"):
        var = "rev_pp_local" if scheme == "revenue" else "seda_ses"
        if var not in d or d[var].isna().all():
            if scheme == "ses" and "pov_rate_5_17" in d:
                d["seda_ses"] = -d.pov_rate_5_17     # poverty as an inverse SES proxy
                var = "seda_ses"
            else:
                raise SystemExit(f"scheme {scheme!r} needs column {var!r}")
        cuts = (d.groupby(["year", "czone"])
                .apply(lambda g: wquantile(g[var], g.pupils, 0.5),
                       include_groups=False)
                .rename("cut").reset_index())
        d = d.merge(cuts, on=["year", "czone"])
        d = d[d[var].notna()]
        d["loc_id"] = np.where(d[var] > d.cut, 2, 1)
    else:
        raise SystemExit(f"unknown scheme {scheme!r}")

    # Keep only CZ-years with a genuine two-sided split.
    tot = d.groupby(["year", "czone"]).pupils.transform("sum")
    grp = d.groupby(["year", "czone", "loc_id"]).pupils.transform("sum")
    d["group_share"] = grp / tot
    ok = (d.groupby(["year", "czone"])
          .agg(n_loc=("loc_id", "nunique"), minshare=("group_share", "min")))
    ok = ok[(ok.n_loc == 2) & (ok.minshare >= min_group_share)].index
    d = d[pd.MultiIndex.from_frame(d[["year", "czone"]]).isin(ok)]
    LOG.info("partition %-8s: %d district-years, %d CZ-years",
             scheme, len(d), d.groupby(["year", "czone"]).ngroups)
    return d


# =============================================================================
# Within-CZ gaps -- the model's spatial moments
# =============================================================================

@dataclass(frozen=True)
class Gap:
    """One within-CZ moment, and how to aggregate a group of districts into it.

    The distinction that matters: the model's M_l and H̃_T,l are LOCATION TOTALS,
    and its per-pupil objects are ratios of totals, not averages of district
    ratios.  Taking a district-weighted mean of enrollment answers "how big is
    the typical district", which is not a question the model asks.

        total   sum the column over the group            -> M_l, H̃_T,l
        ratio   sum(num)/sum(den) over the group         -> per-pupil, per-FTE
        mean    weighted mean of a district-level index  -> scores, deflated pay
    """
    key: str
    label: str
    model_object: str
    stat: str                       # total | ratio | mean
    num: str = ""                   # summed for total/ratio
    den: str = ""                   # denominator for ratio
    col: str = ""                   # averaged for mean
    weight: str = "pupils"          # weight for mean
    contrast: str = "log"           # log ratio (2/1) or diff in levels (2-1)

    def needs(self) -> list[str]:
        return [c for c in (self.num, self.den, self.col,
                            self.weight if self.stat == "mean" else "") if c]


GAPS = [
    Gap("pupils", "enrollment, group total", "M_2/M_1", "total", num="pupils"),
    Gap("teachers", "FTE teachers, group total", "H̃_T,2/H̃_T,1 (before the β weighting)",
        "total", num="teachers"),
    Gap("teachers_pp", "FTE teachers per pupil", "H̃_T,l / M_l",
        "ratio", num="teachers", den="pupils"),
    Gap("str_ratio", "student-teacher ratio", "class size N(h)",
        "ratio", num="pupils", den="teachers"),
    Gap("salary_per_teacher", "teacher salary per FTE, nominal", "κ_2/κ_1 (undeflated)",
        "ratio", num="wagebill_teachers", den="teachers"),
    Gap("salary_real", "teacher salary per FTE, CWIFT-deflated", "κ_2/κ_1",
        "mean", col="salary_real", weight="teachers"),
    Gap("exp_pp", "current expenditure per pupil", "W_l / M_l",
        "ratio", num="exp_current_elsec_total", den="pupils"),
    Gap("rev_pp_local", "local revenue per pupil", "t_l I_l / M_l",
        "ratio", num="rev_local_total", den="pupils"),
    Gap("rev_pp_state", "state revenue per pupil", "G_l (notes §5.2)",
        "ratio", num="rev_state_total", den="pupils"),
    Gap("rev_pp_total", "total revenue per pupil", "W_l / M_l",
        "ratio", num="rev_total", den="pupils"),
    Gap("proptax_pp", "property-tax revenue per pupil", "t_l",
        "ratio", num="rev_local_prop_tax", den="pupils"),
    Gap("t_eff", "local revenue / local income base", "t_l",
        "ratio", num="rev_local_total", den="income_base", contrast="diff"),
    Gap("share_local", "local share of revenue", "the pure local budget of §10",
        "ratio", num="rev_local_total", den="rev_total", contrast="diff"),
    Gap("share_state", "state share of revenue", "G_l",
        "ratio", num="rev_state_total", den="rev_total", contrast="diff"),
    Gap("pov_rate_5_17", "child poverty rate 5-17", "ability/income sorting",
        "ratio", num="pop_5_17_poverty", den="pop_5_17", contrast="diff"),
    Gap("learn_rate", "SEDA learning rate", "Q_2/Q_1",
        "mean", col="learn_rate", contrast="diff"),
    Gap("score_mean", "SEDA mean score", "Q level (composition-laden)",
        "mean", col="score_mean", contrast="diff"),
    Gap("seda_ses", "SEDA composite SES", "ability/income sorting",
        "mean", col="seda_ses", contrast="diff"),
    Gap("seda_lninc", "log median household income", "I_2/I_1 (the income base)",
        "mean", col="seda_lninc", contrast="diff"),
]


# Which `fetch_data.py` source supplies a column, so a blocked moment can say
# what to rerun instead of just disappearing from the report.
COLUMN_SOURCE = {
    "pupils": "CCD directory (--sources urban)",
    "teachers": "CCD directory (--sources urban)",
    "wagebill_teachers": "F-33 (--urban-topics ccd_finance)",
    "rev_local_total": "F-33 (--urban-topics ccd_finance)",
    "rev_state_total": "F-33 (--urban-topics ccd_finance)",
    "rev_total": "F-33 (--urban-topics ccd_finance)",
    "rev_local_prop_tax": "F-33 (--urban-topics ccd_finance)",
    "exp_current_elsec_total": "F-33 (--urban-topics ccd_finance)",
    "salary_real": "NCES EDGE CWIFT (--sources edge)",
    "income_base": "Census ACS (--sources acs, needs a free API key)",
    "pop_5_17": "SAIPE (--urban-topics saipe)",
    "pop_5_17_poverty": "SAIPE (--urban-topics saipe)",
    "learn_rate": "SEDA 6.0 (--sources seda)",
    "score_mean": "SEDA 6.0 (--sources seda)",
    "seda_ses": "SEDA 6.0 district covariates (--sources seda)",
    "seda_lninc": "SEDA 6.0 district covariates (--sources seda)",
}


def _group_stat(gl: pd.DataFrame, g: Gap) -> float:
    """Collapse one location's districts to the single number the model has."""
    if g.stat == "total":
        return float(gl[g.num].sum(min_count=1))
    if g.stat == "ratio":
        m = gl[g.num].notna() & gl[g.den].notna() & (gl[g.den] > 0)
        if not m.any():
            return np.nan
        den = gl.loc[m, g.den].sum()
        return float(gl.loc[m, g.num].sum() / den) if den > 0 else np.nan
    return wmean(gl[g.col], gl[g.weight])


def czone_gaps(panel: pd.DataFrame) -> pd.DataFrame:
    """One row per CZ-year: each moment for location 1, location 2, and the gap."""
    have = [g for g in GAPS if set(g.needs()) <= set(panel.columns)]
    skipped = [g.key for g in GAPS if g not in have]
    if skipped:
        LOG.info("  gaps skipped for want of inputs: %s", ", ".join(skipped))
    rows: list[dict[str, Any]] = []
    for (yr, cz), grp in panel.groupby(["year", "czone"]):
        rec: dict[str, Any] = {"year": yr, "czone": cz,
                               "pupils": grp.pupils.sum(), "n_districts": len(grp)}
        parts = {int(l): gl for l, gl in grp.groupby("loc_id")}
        for l, gl in parts.items():
            rec[f"n_districts_{l}"] = len(gl)
        for g in have:
            lo, hi = (_group_stat(parts[l], g) if l in parts else np.nan
                      for l in (1, 2))
            rec[f"{g.key}_1"], rec[f"{g.key}_2"] = lo, hi
            if g.contrast == "log":
                rec[f"gap_{g.key}"] = (np.log(hi / lo) if lo > 0 and hi > 0
                                       else np.nan)
            else:
                rec[f"gap_{g.key}"] = hi - lo
        rows.append(rec)
    return pd.DataFrame(rows)


def summarize_gaps(gaps: pd.DataFrame, reg: Registry, scheme: str, year: int,
                   panel: pd.DataFrame, report_missing: bool = True) -> None:
    """Turn each within-CZ gap distribution into an `Estimate`."""
    g = gaps[gaps.year == year]
    if g.empty:
        LOG.warning("no CZ gaps for %d", year)
        return
    src = (f"CCD directory + F-33 + SAIPE + CWIFT + SEDA, {year}; "
           f"district x CZ ({scheme} split), enrollment-weighted across "
           f"{len(g)} CZs")
    for gv in GAPS:
        col = f"gap_{gv.key}"
        if col not in g or g[col].notna().sum() < 5:
            if report_missing:
                lack = [c for c in gv.needs() if c not in panel.columns]
                reg.missing(f"gap_{gv.key}", f"within-CZ gap, {gv.label}",
                            gv.model_object,
                            "; ".join(dict.fromkeys(COLUMN_SOURCE.get(c, c)
                                                    for c in lack)) or
                            "present but too few CZs with both locations",
                            GAP_NOTES.get(gv.key, "not computable from the tables "
                                                  "currently in raw/"))
            continue
        val = wmean(g[col], g.pupils)
        reg.add(Estimate(
            key=f"gap_{gv.key}_{scheme}",
            label=f"within-CZ gap, {gv.label} (location 2 - location 1)",
            model_object=gv.model_object,
            value=val,
            units="log points" if gv.contrast == "log" else "levels",
            source=src,
            note=GAP_NOTES.get(gv.key, "within-CZ dispersion the two-location "
                                       "model has to reproduce"),
            n=int(g[col].notna().sum()),
            extra={"sd": float(g[col].std()),
                   "p25": wquantile(g[col], g.pupils, 0.25),
                   "p50": wquantile(g[col], g.pupils, 0.50),
                   "p75": wquantile(g[col], g.pupils, 0.75),
                   "ratio_2_over_1": (float(np.exp(val))
                                      if gv.contrast == "log" else None),
                   "level_1": wmean(g[f"{gv.key}_1"], g.pupils),
                   "level_2": wmean(g[f"{gv.key}_2"], g.pupils),
                   "within_cz_var_share": within_share(panel, gv.key,
                                                       gv.contrast, year)}))


GAP_NOTES = {
    "pupils": "M_2/M_1, the enrollment split. B_l is the residual backed out to "
              "match it, so this is a target, never a check.",
    "teachers": "The teacher stock by location. Together with the enrollment gap "
                "it says whether the good location's advantage is more teachers "
                "or better ones -- the model's answer is 'better', through β.",
    "teachers_pp": "This is H̃_T,l/M_l up to the β-weighting of teacher human "
                   "capital, so it maps into Q_l = (2H̃_T,l/M_l)^σ directly. If "
                   "it is near zero while the score gap is not, the quality gap "
                   "is composition, not class size.",
    "str_ratio": "Class size. Prop. 1 has better teachers take LARGER classes, "
                 "which is not what districts do -- check the mapping (notes "
                 "§5.3) before targeting this.",
    "salary_real": "The direct target for κ_2/κ_1; CWIFT-deflated, so it is a "
                   "teaching-wage gap rather than a local price level.",
    "salary_per_teacher": "Compare with the deflated version: the wedge between "
                          "the two is how much of the raw salary gap is cost of "
                          "living rather than pay.",
    "exp_pp": "The resource gap the budget constraint has to deliver.",
    "rev_pp_local": "The local half of the budget; with the state half it pins "
                    "how far pure local finance is from the data.",
    "rev_pp_state": "Motivates the exogenous transfer G_l of notes §5.2 -- state "
                    "aid is compensatory, so expect this gap to be NEGATIVE even "
                    "where the local gap is positive.",
    "rev_pp_total": "Local plus state plus federal. If this gap has the opposite "
                    "sign to the local one, redistribution is doing more than "
                    "local finance, and a pure local budget gets the sign of the "
                    "resource gap wrong.",
    "proptax_pp": "The property-tax base per pupil, the closest thing to t_l I_l "
                  "in the data.",
    "t_eff": "The effective local tax rate. Keep it endogenous from the budget "
             "and use this as an untargeted check.",
    "share_local": "If this is far below 1 the pure local-finance budget in §10 "
                   "is counterfactual -- see G_l.",
    "share_state": "The size of the transfer the model is currently missing.",
    "pov_rate_5_17": "The observable counterpart of sorting on z: how much richer "
                     "the good location's families are.",
    "learn_rate": "The closest public analogue to Q_2/Q_1; growth rather than "
                  "levels, so it nets out the SES composition the model is "
                  "supposed to generate endogenously.",
    "score_mean": "NOT a Q target -- it mixes school quality with who enrolls. "
                  "Useful as the untargeted check that the model's sorting "
                  "reproduces a level gap it did not fit.",
    "seda_lninc": "The income base I_l, in logs. This is what t_l is levied on, "
                  "so together with the local-revenue gap it says whether the "
                  "advantaged location raises more because it is richer or "
                  "because it taxes itself harder.",
    "seda_ses": "Same as the poverty gap, on SEDA's composite. Reported in SEDA "
                "index units, whose cross-district SD is about 0.9 -- divide by "
                "that before setting it beside the model's mean-ability gap "
                "(`idx_work_signed`, ≈ +0.08 at the baseline), which is in units "
                "of z. The sign is the part that compares directly.",
}


def within_share(panel: pd.DataFrame, var: str, kind: str,
                 year: int | None = None) -> float | None:
    """Share of the cross-district variance that is WITHIN commuting zones.

    The single most important descriptive number here: if most dispersion in
    teacher pay and school quality is within CZ, then a within-CZ two-location
    model is where the action is, and Option 3 (across metros) is not.

    Computed on ONE cross-section.  Pooling 1994-2018 would count two decades of
    nominal drift in every dollar variable as within-CZ variance -- each CZ
    appears in every year -- and inflate the share for exactly the moments the
    number is supposed to adjudicate.
    """
    if var not in panel or "czone" not in panel:
        return None
    if year is not None and "year" in panel:
        panel = panel[panel.year == year]
    x = panel[var]
    x = np.log(x.where(x > 0)) if kind == "log" else x
    d = pd.DataFrame({"x": x, "cz": panel.czone, "w": panel.pupils}).dropna()
    if len(d) < 100:
        return None
    tot = np.average((d.x - np.average(d.x, weights=d.w)) ** 2, weights=d.w)
    dm = d.x - (d.w * d.x).groupby(d.cz).transform("sum") / d.w.groupby(d.cz).transform("sum")
    return float(np.average(dm ** 2, weights=d.w) / tot) if tot > 0 else None


# =============================================================================
# Regressions -- the parameters that need a slope, not a gap
# =============================================================================

def run_regressions(panel: pd.DataFrame, reg: Registry, year: int) -> dict[str, Any]:
    """Within-CZ slopes for σ, σ_ν, the κ gradient, and budget pass-through."""
    d = panel[panel.year == year].copy()
    d["cz"] = d.czone
    out: dict[str, Any] = {}
    # SEDA publishes the learning rate as SDs gained per grade; the model's Q_l
    # is a level, so the regressions below work with the gain cumulated over the
    # grades SEDA fits (3-8).  The gap moment stays in the published unit.
    if "learn_rate" in d:
        d["learn_gain"] = d.learn_rate * SEDA_GRADE_SPAN

    def have(*cols: str) -> bool:
        """A column full of NaN is a missing input, not a present one."""
        return all(c in d.columns and d[c].notna().any() for c in cols)

    def why(*cols: str) -> str:
        """Say which of it is a download and which is a coverage year.

        CWIFT is fetched but starts in 2013-14 and SEDA pools 2009-2019, so a
        moment can be blocked in 2010 by data that is sitting in raw/. Naming a
        `fetch_data.py` flag in that case would send the reader off to re-run a
        download that would change nothing.
        """
        empty = [c for c in cols if c in panel.columns
                 and panel[c].notna().any() and not d[c].notna().any()]
        if not empty:
            return ""
        return (f" Note that {', '.join(empty)} IS in raw/ -- it just does not "
                f"cover {year}; try a base year the source reaches.")

    def thin(key: str, label: str, model_object: str, source: str) -> None:
        """`ols_fe` returned None: the inputs are there but do not overlap."""
        reg.missing(key, label, model_object, source,
                    f"The inputs are in raw/ but too few districts have all of "
                    f"them in {year} for a within-CZ regression. Widen --years, "
                    f"or check the merge coverage reported above.")

    # --- σ: the congestion elasticity in Q_l = (2 H̃_T,l / M_l)^σ -------------
    # log Q = σ log(teachers per pupil) + const, and log Q is in log-h units, so
    # the SEDA learning rate must first be converted with SD_SCORE_TO_LOG_EARNINGS.
    if not have("learn_rate", "log_teachers_pp"):
        reg.missing("sigma_congestion", "class-size / congestion elasticity σ",
                    "σ in Q_l = (2 H̃_T,l / M_l)^σ",
                    "SEDA 6.0 learning rates (--sources seda) x CCD teachers per pupil",
                    "σ is a free 0.5 in the aspatial draft and 0.25 in the code. In "
                    "the spatial model it governs whether move-to-opportunity is "
                    "self-defeating in equilibrium, so it should be estimated." + why("learn_rate"))
    else:
        r = ols_fe(d, "learn_gain", ["log_teachers_pp"], fe="cz", w="pupils",
                   cluster="cz")
        out["sigma"] = r
        if r is None:
            thin("sigma_congestion", "class-size / congestion elasticity σ",
                 "σ in Q_l = (2 H̃_T,l / M_l)^σ",
                 "SEDA 6.0 learning rates x CCD teachers per pupil")
        else:
            b, se = r["coef"]["log_teachers_pp"], r["se"]["log_teachers_pp"]
            reg.add(Estimate(
                key="sigma_congestion", label="class-size / congestion elasticity σ",
                model_object="σ in Q_l = (2 H̃_T,l / M_l)^σ",
                value=b * SD_SCORE_TO_LOG_EARNINGS,
                se=abs(se * SD_SCORE_TO_LOG_EARNINGS), n=r["n"],
                units="elasticity of log Q w.r.t. teachers per pupil",
                source=f"SEDA learning rates x CCD teachers per pupil, {year}, "
                       f"district level, CZ fixed effects, {r['n_clusters']} CZ clusters",
                status="assumption",
                note=f"The slope is {precision(b, se)}. "
                     f"σ is a free 0.25 in the code. This converts the raw slope "
                     f"({b:+.3f} SD of achievement gained over grades 3-8 per log "
                     f"point of teachers per pupil) into log-h units "
                     f"with 1 SD ≈ {SD_SCORE_TO_LOG_EARNINGS} log earnings -- the "
                     f"conversion is an assumption, and OLS here is a correlation, "
                     f"not a class-size experiment. Use as a prior, then discipline "
                     f"with STAR / Angrist-Lavy.",
                extra={"raw_slope": b, "raw_se": se,
                       "conversion": SD_SCORE_TO_LOG_EARNINGS,
                       "grade_span": SEDA_GRADE_SPAN,
                       "r2_within": r["r2_within"]}))

    # --- σ_ν: the Gumbel taste scale from the location-share elasticity -------
    # log(M_2/M_1) = (ΔB + ΔV)/σ_ν, so regressing the log enrollment share on the
    # quality differential gives 1/σ_ν once quality is in consumption units.
    if not have("learn_rate", "pupils"):
        reg.missing("sigma_nu", "location-taste dispersion σ_ν",
                    "σν (Gumbel scale on the location choice)",
                    "SEDA 6.0 learning rates (--sources seda) x CCD enrollment shares",
                    f"Currently σν = {LIT_SIGMA_NU} in the code, imported from "
                    f"Eckert-Kleineberg rather than estimated. The within-CZ share "
                    f"regression is the way to earn it -- and it needs an instrument.")
    else:
        d["log_share_cz"] = np.log(d.pupils / d.groupby("cz").pupils.transform("sum"))
        d["q_util"] = d.learn_gain * SD_SCORE_TO_LOG_EARNINGS
        r = ols_fe(d, "log_share_cz", ["q_util"], fe="cz", w="pupils", cluster="cz")
        out["sigma_nu"] = r
        if r is None or abs(r["coef"]["q_util"]) <= 1e-8:
            thin("sigma_nu", "location-taste dispersion σ_ν",
                 "σν (Gumbel scale on the location choice)",
                 "CCD enrollment shares x SEDA learning rates")
        elif r:
            slope = r["coef"]["q_util"]
            sign_clause = ("; a non-positive slope has no σ_ν behind it, hence "
                           "the blank value" if slope <= 0 else "")
            reg.add(Estimate(
                key="sigma_nu", label="location-taste dispersion σ_ν",
                model_object="σν (Gumbel scale on the location choice)",
                value=1.0 / slope if slope > 0 else None,
                n=r["n"], units="consumption-equivalent utils",
                source=f"CCD enrollment shares x SEDA learning rates, {year}, "
                       f"district within CZ, CZ fixed effects",
                status="assumption",
                note=f"Slope {slope:+.3f} = 1/σ_ν under the logit, and it is "
                     f"{precision(slope, r['se']['q_util'])}{sign_clause}. "
                     f"OLS is biased "
                     f"here -- quality is jointly determined with the share "
                     f"through congestion -- so treat it as an upper bound on "
                     f"1/σ_ν and keep Eckert-Kleineberg's {LIT_SIGMA_NU} as the "
                     f"prior. Needs an instrument (boundary discontinuities, "
                     f"state-aid formula kinks).",
                extra={"slope": slope, "se": r["se"]["q_util"],
                       "literature_prior": LIT_SIGMA_NU}))

    # --- the κ gradient: do better districts pay teachers more? --------------
    if not have("log_salary_real", "learn_rate"):
        reg.missing("kappa_quality_gradient",
                    "log real teacher salary on district learning rate, within CZ",
                    "κ_l alongside Q_l",
                    "NCES EDGE CWIFT (--sources edge) + SEDA 6.0 (--sources seda)",
                    "Whether κ_l and Q_l co-move or offset decides how much of the "
                    "sorting the paper attributes to the teacher spillover is "
                    "really just the wage." + why("salary_real", "learn_rate"))
    else:
        xs = [c for c in ("learn_gain", "log_str_ratio", "pov_rate_5_17")
              if c in d.columns]
        r = ols_fe(d, "log_salary_real", xs, fe="cz", w="teachers", cluster="cz")
        out["kappa_gradient"] = r
        if r is None:
            thin("kappa_quality_gradient",
                 "log real teacher salary on district learning rate, within CZ",
                 "κ_l alongside Q_l",
                 "CWIFT-deflated F-33 salaries x SEDA learning rates")
        else:
            b = r["coef"].get("learn_gain")
            reg.add(Estimate(
                key="kappa_quality_gradient",
                label="log real teacher salary on district learning rate, within CZ",
                model_object="κ_l alongside Q_l",
                value=b, se=r["se"].get("learn_gain"), n=r["n"],
                units="log points per SD of cumulative achievement growth",
                source=f"F-33 teacher salaries / CCD FTE, CWIFT-deflated, x SEDA, "
                       f"{year}, CZ fixed effects",
                note=f"The gradient is {precision(b, r['se'].get('learn_gain', float('nan')))}. "
                     "Sign tells you whether κ_l and Q_l co-move (agglomeration of "
                     "advantage) or offset (compensating differentials). The model "
                     "sets both exogenously; if they co-move strongly in the data, "
                     "κ_2 > κ_1 is doing the sorting work the paper attributes to "
                     "the spillover, and the two need separating.",
                extra={"controls": xs, "coefs": r["coef"], "ses": r["se"]}))

    # --- budget pass-through: does local revenue become teacher pay? ---------
    if not have("log_salary_per_teacher", "log_rev_pp_local"):
        reg.missing("revenue_to_salary_passthrough",
                    "elasticity of teacher salary to local revenue per pupil",
                    "the balanced budget t_l I_l = W_l",
                    "F-33 salaries and local revenue (--urban-topics ccd_finance)",
                    "Without it there is no evidence on how much of a local "
                    "dollar reaches teacher pay, which is the whole content of "
                    "the model's budget constraint.")
    else:
        r = ols_fe(d, "log_salary_per_teacher", ["log_rev_pp_local"], fe="cz",
                   w="teachers", cluster="cz")
        out["passthrough"] = r
        if r is None:
            thin("revenue_to_salary_passthrough",
                 "elasticity of teacher salary to local revenue per pupil",
                 "the balanced budget t_l I_l = W_l", "F-33")
        else:
            b = r["coef"]["log_rev_pp_local"]
            reg.add(Estimate(
                key="revenue_to_salary_passthrough",
                label="elasticity of teacher salary to local revenue per pupil",
                model_object="the balanced budget t_l I_l = W_l",
                value=b, se=r["se"]["log_rev_pp_local"], n=r["n"], units="elasticity",
                source=f"F-33 {year}, district within CZ, CZ fixed effects",
                note=f"The elasticity is {precision(b, r['se']['log_rev_pp_local'])}. "
                     "The model spends the whole local budget on teachers, i.e. an "
                     "elasticity of 1. Anything well below 1 says the budget also "
                     "buys class-size reduction, facilities and non-teaching staff, "
                     "and that the mapping from t_l to κ_l needs a wedge.",
                extra={"r2_within": r["r2_within"]}))

    # --- the state-aid share, which motivates G_l -----------------------------
    if not have("share_local"):
        reg.missing("local_revenue_share",
                    "local share of district revenue (enrollment-weighted mean)",
                    "the pure local budget of §10 vs. the G_l extension",
                    "F-33 revenue by source (--urban-topics ccd_finance)",
                    "The size of the state transfer the model currently omits.")
    else:
        reg.add(Estimate(
            key="local_revenue_share",
            label="local share of district revenue (enrollment-weighted mean)",
            model_object="the pure local budget of §10 vs. the G_l extension",
            value=wmean(d.share_local, d.pupils), n=int(d.share_local.notna().sum()),
            units="share",
            source=f"F-33 {year}, all sample districts, enrollment-weighted",
            note=f"Benchmark: NCES reports {LIT_LOCAL_REV_SHARE:.2f} nationally. "
                 f"A pure local budget describes no year in the sample, which is "
                 f"the case for adding the exogenous transfer G_l (notes §5.2) "
                 f"and makes finance equalization the natural counterfactual.",
            extra={f"share_{s}": (wmean(d[f"share_{s}"], d.pupils)
                                  if f"share_{s}" in d else None)
                   for s in ("state", "fed")}
            | {"nces_national": LIT_LOCAL_REV_SHARE}))
    return out


# =============================================================================
# Trends and literature anchors
# =============================================================================

def trend_table(gaps: pd.DataFrame, reg: Registry, scheme: str) -> pd.DataFrame:
    """Gaps by year -- the spatial counterpart of the 1970/1990/2010 steady states."""
    cols = [c for c in gaps.columns if c.startswith("gap_")]
    rows = []
    for yr, g in gaps.groupby("year"):
        rec = {"year": int(yr), "n_czones": len(g)}
        rec.update({c: wmean(g[c], g.pupils) for c in cols})
        rows.append(rec)
    tr = pd.DataFrame(rows).sort_values("year")
    for c in ("gap_salary_real", "gap_rev_pp_local", "gap_exp_pp", "gap_pupils",
              "gap_share_local"):
        if c in tr and tr[c].notna().sum() >= 2:
            first, last = tr[tr[c].notna()].iloc[0], tr[tr[c].notna()].iloc[-1]
            reg.add(Estimate(
                key=f"trend_{c}_{scheme}",
                label=f"change in the within-CZ {c[4:]} gap, "
                      f"{int(first.year)}->{int(last.year)}",
                model_object="comparative statics across steady states",
                value=float(last[c] - first[c]), units="log points (or levels)",
                source=f"F-33 / CCD panel, {scheme} split, CZ-years",
                note="Whether the spatial gaps widened tells you whether a single "
                     "calibrated steady state is enough or the spatial block needs "
                     "its own time series. See notes §7: F-33 is thin before 1990.",
                extra={"first_year": int(first.year), "last_year": int(last.year),
                       "first": float(first[c]), "last": float(last[c])}))
    return tr


def add_literature_anchors(reg: Registry) -> None:
    """Parameters no district-level public table can identify: carry the source."""
    reg.add(Estimate(
        key="lambda_wtp_level", label="average WTP for school quality",
        model_object="λ (altruism strength)", value=None, units="$/month, 1990",
        source="Bayer, Ferreira & McMillan (2007), NBER 13236, Tables 5-7",
        status="literature",
        note="~$33/month per SD of test score with boundary fixed effects, ~$17 "
             "once neighborhood sociodemographics are controlled -- 1-2% of "
             "monthly housing user cost. This pins the LEVEL of λ; ψ is pinned "
             "by how WTP varies with parental income (BFM Table 7). Target them "
             "separately or they fight for the same moment (notes §4)."))
    reg.add(Estimate(
        key="psi_wtp_gradient", label="income gradient in WTP for school quality",
        model_object="ψ (warm-glow curvature)", value=None, units="elasticity",
        source="Bayer, Ferreira & McMillan (2007) Table 7; ACS district "
               "composition by income",
        status="literature",
        note="The baseline sets ψ = 0 (linear in h'), which is the strongest "
             "sorting the power kernel allows. The gradient in BFM says how much "
             "of that is real. This is the single most valuable number still "
             "missing from the empirical section."))
    reg.add(Estimate(
        key="moving_cost_m", label="school-quality capitalization at district borders",
        model_object="mcost / τmove (the goods cost of moving)",
        value=None, units="share of consumption, annualized",
        source="Black (1999, QJE); Bayer, Ferreira & McMillan (2007)",
        status="literature",
        note="The boundary-discontinuity house-price premium is the observed "
             "price of moving between locations, so `mcost` stops being a free "
             "normalization (currently MCOST_BENCH ≈ 0.0294, reverse-engineered "
             "from a 0.20 utility cost) and becomes a measured object."))
    reg.add(Estimate(
        key="rho_z_sigma_xi", label="parent-child ability persistence",
        model_object="ρz, σξ", value=None, units="AR(1) coefficient",
        source="CNLSY (NLSY79 mothers x Child/Young Adult PIAT and AFQT)",
        status="literature",
        note="Must be the parent-child COGNITIVE correlation, not the income IGE: "
             "z is ability. Currently ρz = 0.9, σξ = 0.20 by assumption. "
             "Cross-check against the stationary AFQT dispersion already matched "
             "in the aspatial calibration."))
    reg.add(Estimate(
        key="sigma_nu_prior", label="spatial taste dispersion, literature prior",
        model_object="σν", value=LIT_SIGMA_NU, units="consumption-equivalent utils",
        source="Eckert & Kleineberg, 'Saving the American Dream?', spatial GE "
               "with education policy",
        status="literature",
        note="Equals the code's current σν = 0.20 exactly. Reassuring, but it "
             "means σν is presently an import rather than an estimate -- the "
             "within-CZ share regression above is the way to earn it."))
    reg.add(Estimate(
        key="exposure_effect", label="neighborhood exposure effect on adult income rank",
        model_object="untargeted validation of the Q_l mechanism",
        value=LIT_EXPOSURE_RATE, units="convergence per year of childhood exposure",
        source="Chetty & Hendren (2018), NBER 23001; Opportunity Atlas",
        status="literature",
        note="The single best untargeted test available. Simulate a child moved "
             "from location 1 to location 2 at each age and read off the implied "
             "convergence rate; ~4%/year without having targeted it would "
             "validate the whole Q_l channel (notes §6)."))


DATA_GAPS = [
    {"what": "Residence-workplace split for teachers",
     "why": "Stage 2 gives one Gumbel draw for where you live, pay tax and teach. "
            "Within a metro that is the binding assumption: city teachers commonly "
            "live in the suburbs. Two independent logits plus a commuting cost buy "
            "two independent moment sets (family sorting AND teacher flows).",
     "source": "SASS/NTPS restricted-use; state administrative files (WI via "
               "Biasi-Fu-Stromme, NC, TX); Bates-Dinerstein-Johnston-Sorkin on "
               "commute as the dominant teacher preference",
     "blocks": "notes §5.1"},
    {"what": "Teacher ability composition by location",
     "why": "The paper's distinctive claim is about WHO teaches, but every moment "
            "above is a quantity or a dollar. Nothing here measures teacher "
            "ability, so β and γ are untouched by the data.",
     "source": "NLSY79/97 restricted geocode files (county/MSA of residence) "
               "crossed with the AFQT-by-occupation moments already in hand",
     "blocks": "β, γ, the H̃_T aggregator"},
    {"what": "An instrument for the quality-share relationship",
     "why": "σ_ν and σ are both estimated off an equilibrium relationship between "
            "school quality and enrollment that the model says is simultaneous.",
     "source": "district-boundary discontinuities; state-aid formula kinks "
               "(Lafortune-Rothstein-Schanzenbach); court-ordered finance reforms",
     "blocks": "σ, σν"},
    {"what": "A pre-1990 cross-section",
     "why": "The aspatial calibration has 1970/1990/2010. F-33 is thin before 1990 "
            "and SEDA starts in 2009, so the spatial block realistically calibrates "
            "to one date.",
     "source": "Census of Governments school district finances (published back to "
               "1957); Digest of Education Statistics",
     "blocks": "notes §7"},
    {"what": "Tract -> school district crosswalk for the Opportunity Atlas",
     "why": "OI is tract/county/CZ, never district. The crosswalk is what turns the "
            "exposure-effect validation from a hand-wave into a number.",
     "source": "NCES EDGE district boundaries, population-weighted onto tracts",
     "blocks": "the untargeted validation in notes §6"},
    {"what": "Class size vs. teacher quality within districts",
     "why": "Prop. 1 implies N(h) = h^(β/σ) Q^(-1/σ): better teachers get BIGGER "
            "classes. If districts do the opposite, the model's class sizes cannot "
            "be read as observed student-teacher ratios.",
     "source": "CRDC school-level teacher experience/certification "
               "(--urban-topics crdc_teachers); state admin data",
     "blocks": "notes §5.3, and whether gap_str_ratio is a legitimate target"},
]


# =============================================================================
# Reporting
# =============================================================================

def write_reports(reg: Registry, meta: dict[str, Any], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {**meta, "estimates": reg.to_dict(), "data_gaps": DATA_GAPS}
    (outdir / "spatial_moments.json").write_text(
        json.dumps(payload, indent=2, default=str))

    L: list[str] = []
    add = L.append
    add("# Spatial moments for `spatial_continuous.jl`\n")
    add(f"Generated {meta['generated']} by `data_estimate.py`. "
        f"Base year **{meta['base_year']}**, partition scheme "
        f"**{meta['scheme']}**, labor market **{meta['geo_unit']}**; a location "
        f"is a group of school districts inside one labor market, and location 2 "
        f"is the advantaged group.\n")

    add("\n## Inputs\n")
    add("| table | status | columns |")
    add("|---|---|---|")
    for k, v in sorted(meta["inputs"].items()):
        detail = (v.get("n_cols") if "n_cols" in v else
                  f"{v['n_tables']} tables" if "n_tables" in v else "")
        add(f"| `{k}` | {v['status']} | {detail} |")

    add(f"\n`within-CZ var share` is the fraction of the {meta['base_year']} "
        "cross-district variance in that variable that survives commuting-zone "
        "fixed effects. It is the "
        "descriptive case for the whole design: where it is high, a within-CZ "
        "two-location model is where the variation is, and a cross-metro model "
        "(Option 3 in the calibration notes) would be explaining the smaller half.\n")
    add(f"\nPanel: {meta['panel']['district_years']} district-years, "
        f"{meta['panel']['districts']} districts, {meta['panel']['czones']} CZs, "
        f"{meta['panel']['cz_years_partitioned']} CZ-years surviving the split.\n")

    for status, title in (("ok", "Estimates from the data"),
                          ("assumption", "Estimates that rest on a stated assumption"),
                          ("literature", "Literature anchors (not estimated here)"),
                          ("missing-input", "Blocked on missing data")):
        items = [e for e in reg.items if e.status == status]
        if not items:
            continue
        add(f"\n## {title}\n")
        add("| key | model object | value | ratio 2/1 | within-CZ var share | n |")
        add("|---|---|---|---|---|---|")
        for e in items:
            v = "—" if e.value is None else f"{e.value:.4g}"
            if e.se:
                v += f" ({e.se:.3g})"
            r = e.extra.get("ratio_2_over_1")
            w = e.extra.get("within_cz_var_share")
            add(f"| `{e.key}` | {e.model_object} | {v} | "
                f"{'—' if r is None else f'{r:.3f}'} | "
                f"{'—' if w is None else f'{w:.2f}'} | {e.n or '—'} |")
        add("")
        for e in items:
            add(f"- **`{e.key}`** — {e.label}. {e.note}  \n  *Source:* {e.source}")

    add("\n## What the empirical section still needs\n")
    for g in DATA_GAPS:
        add(f"- **{g['what']}** ({g['blocks']}). {g['why']}  \n"
            f"  *Source:* {g['source']}")

    (outdir / "spatial_moments.md").write_text("\n".join(L) + "\n")
    LOG.info("wrote %s and %s", outdir / "spatial_moments.json",
             outdir / "spatial_moments.md")


# =============================================================================
# CLI
# =============================================================================

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indir", default=str(HERE),
                   help="root holding raw/ (default: this script's directory)")
    p.add_argument("--outdir", default=None,
                   help="where the report goes (default: <indir>/estimates)")
    p.add_argument("--years", nargs="+", type=int,
                   default=[1994, 2000, 2006, 2010, 2014, 2018],
                   help="cross-sections to build the panel from")
    p.add_argument("--base-year", type=int, default=2018,
                   help="the year the headline moments are computed for")
    p.add_argument("--geo", default="auto", choices=["auto", "czone", "cbsa"],
                   help="labor market the two locations sit inside (default: "
                        "commuting zone, falling back to CBSA if the crosswalk "
                        "has not been downloaded)")
    p.add_argument("--scheme", default="locale",
                   choices=["locale", "revenue", "ses", "all"],
                   help="how to split each CZ into two locations")
    p.add_argument("--min-pupils", type=int, default=500,
                   help="drop districts below this enrollment")
    p.add_argument("--min-group-share", type=float, default=0.10,
                   help="each location must hold at least this share of CZ enrollment")
    p.add_argument("--dump-panel", action="store_true",
                   help="also write the assembled district panel to estimates/")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")

    indir = Path(args.indir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve() if args.outdir else indir / "estimates"
    store = RawStore(indir / "raw")

    LOG.info("== inputs ==")
    expected = ["urban_ccd_directory", "urban_ccd_finance", "urban_ccd_enrollment",
                "urban_saipe", "urban_edfacts_assessments", "urban_edfacts_grad_rates",
                "seda_", "edge_cwift", "edge_geocode_lea", "czone_crosswalk",
                "acs_school_districts"]
    inputs = store.inventory(expected)
    for k, v in sorted(inputs.items()):
        LOG.info("  %-30s %s", k, v["status"])

    years = sorted(set(args.years) | {args.base_year})
    LOG.info("== assembling panel ==")
    panel = assemble(store, years, args.geo)
    LOG.info("labor-market unit: %s", panel.attrs["geo_unit"])

    reg = Registry()
    schemes = ["locale", "revenue", "ses"] if args.scheme == "all" else [args.scheme]
    gaps_all, trends, cz_years = {}, {}, 0
    for scheme in schemes:
        LOG.info("== moments: %s split ==", scheme)
        try:
            part = partition(panel, scheme, args.min_pupils, args.min_group_share)
        except SystemExit as exc:
            LOG.warning("  skipping %s: %s", scheme, exc)
            reg.missing(f"gaps_{scheme}", f"within-CZ gaps under the {scheme} split",
                        "all spatial moments", "see raw/ inventory above", str(exc))
            continue
        g = czone_gaps(part)
        gaps_all[scheme] = g
        cz_years = max(cz_years, g[g.year == args.base_year].shape[0])
        summarize_gaps(g, reg, scheme, args.base_year, part,
                       report_missing=(scheme == schemes[0]))
        trends[scheme] = trend_table(g, reg, scheme)
        if scheme == schemes[0]:
            LOG.info("== regressions ==")
            run_regressions(part, reg, args.base_year)

    LOG.info("== literature anchors ==")
    add_literature_anchors(reg)

    meta = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "base_year": args.base_year, "scheme": args.scheme,
        "geo_unit": panel.attrs["geo_unit"],
        "years": years, "inputs": inputs,
        "sample": {"agency_types": list(REGULAR_AGENCY_TYPES),
                   "min_pupils": args.min_pupils,
                   "min_group_share": args.min_group_share,
                   "ratio_trim_tail": 0.005},
        "panel": {"district_years": len(panel), "districts": int(panel.leaid.nunique()),
                  "czones": int(panel.czone.nunique()) if "czone" in panel else 0,
                  "cz_years_partitioned": int(cz_years)},
        "assumptions": {"sd_score_to_log_earnings": SD_SCORE_TO_LOG_EARNINGS},
    }
    write_reports(reg, meta, outdir)

    outdir.mkdir(parents=True, exist_ok=True)
    for scheme, g in gaps_all.items():
        g.to_csv(outdir / f"czone_gaps_{scheme}.csv", index=False)
    for scheme, t in trends.items():
        t.to_csv(outdir / f"gap_trends_{scheme}.csv", index=False)
    if args.dump_panel:
        panel.to_parquet(outdir / "district_panel.parquet", index=False)
    LOG.info("done -- %d estimates, %d blocked on missing data",
             len(reg.items), sum(e.status == "missing-input" for e in reg.items))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
