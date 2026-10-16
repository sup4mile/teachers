#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
district_settings.py -- the audited district moments under three "settings", the
way a commuting zone (CZ) is cut into two locations (calibration item
"revenue-based setting", 2026-10-02).

The locale calibration targets are the audited, bootstrapped versions of the
tracked moments (`district_audit.py`, `district_bootstrap.py`):

    gap_pupils       enrollment gap, log(M2/M1)                       (tracked definition)
    gap_salary       audit variant 28: CWIFT-deflated salary per K-12 FTE, zero-coded
                     F-33 wage bills filled with 0.881 x instruction salaries, NYC as
                     one district, pre-K teacher FTE out of the denominator
    gap_teachers_pp  audit FTE variant 11: (teachers - pre-K teachers) /
                     (pupils - pre-K enrollment)
    gap_income       `pm_ftyr_all`: log pooled median earnings of full-time year-round
                     workers (ACS B20005 bins summed over a location's districts);
                     `pm_earn_all` (B20001, all workers) is carried beside it

This script computes that set, with the same corrections, for each partition
scheme of `data_estimate.partition` -- `locale` (city / suburb), `revenue`
(below / above the CZ's pupil-weighted median of local revenue per pupil) and `ses`
(SEDA SES) -- at the same sample rules (base year 2018, >= 500 pupils per district,
each location >= 10% of the CZ's partitioned pupils).  Location 2 is the advantaged
group throughout.  Every moment gets a whole-CZ bootstrap (B = 1999, seed 20260929,
one multiplicity matrix per CZ set, shared by all statistics): value, SE, 95%
percentile interval, number of CZs, and the full correlation matrix of the draws.

Also computed:
  * validation moments (child-poverty gap, SEDA growth and score gaps, SEDA SES and
    log income, revenue per pupil by source, current expenditure, teacher payroll
    per K-12 pupil, local revenue share) under every scheme;
  * coverage of each partition, and the "common CZs" (in all three partitions) on
    which the three settings are compared like for like (the pairwise overlaps with
    locale are in the JSON);
  * the state-aid offset under the revenue split: a dollar decomposition of the
    revenue gaps and within-CZ district regressions of state (and other) revenue on
    local revenue per pupil;
  * the pre-registered hold-out split of the locale CZs (pairs of CZs ordered by
    partitioned pupils, one drawn from `default_rng(20261002)` per pair): half A is
    reported; half B goes ONLY to `sealed/holdout_B.json` and is never printed;
  * a sensitivity in which NYC is one fully observed district (F-33 revenue, SAIPE,
    SEDA and ACS records filed under the Chancellor's id), not part of the targets.

The panel, partitions, `czone_gaps`, screen, bootstrap and audit corrections are
`data_estimate.py`, `district_bootstrap.py` and `district_audit.py`'s own functions;
those scripts are unchanged.  The locale full-sample targets must reproduce
0.0086 (0.0092), -0.0177 (0.0131), 0.3764 (0.1154) and 0.1122 (0.0168); the script
stops if they do not.

Usage
    python3 district_settings.py                  # B = 1999, seed 20260929
    python3 district_settings.py --B 199 --outdir /tmp/settings_test    # quick check
Writes estimates/settings/district_settings.{json,md},
district_settings_targets.toml, district_settings_draws_<set>.csv and
sealed/holdout_B.json (about a minute on a compute node).

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_estimate as de  # noqa: E402
import district_audit as da  # noqa: E402
import district_bootstrap as db  # noqa: E402

LOG = logging.getLogger("district_settings")
YEAR = db.YEAR
SCHEMES = ("locale", "revenue", "ses")
MIN_PUPILS = 500
MIN_GROUP_SHARE = 0.10
HOLDOUT_SEED = 20261002
NYC = da.NYC_CHANCELLOR
SCREEN_TAIL = 0.005

# The four locale targets as reported in district_audit.md / district_bootstrap.md,
# (value, bootstrap SE) to the digits those reports print.
REPORTED = {"gap_salary": (0.0086, 0.0092), "gap_teachers_pp": (-0.0177, 0.0131),
            "gap_pupils": (0.3764, 0.1154), "gap_income": (0.1122, 0.0168)}
TARGET_KEYS = ["gap_pupils", "gap_salary", "gap_teachers_pp", "gap_income"]
TOML_KEYS = TARGET_KEYS + ["poverty_gap", "seda_growth_gap", "score_gap"]


# --------------------------------------------------------------------------
# Moment definitions
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Moment:
    key: str        # name in the outputs: the model's key where there is one
    col: str        # column of the CZ table ("" for the level statistic)
    role: str       # target | validation | diagnostic | memo | usd
    kind: str       # log (log ratio 2/1) | diff (levels, 2 - 1) | level | usd
    label: str


MOMENTS: list[Moment] = [
    Moment("gap_pupils", "gap_pupils", "target", "log",
           "enrollment, log(M2/M1) (CCD enrollment incl. pre-K)"),
    Moment("gap_salary", "gap_salary_audited", "target", "log",
           "CWIFT-deflated teacher salary per K-12 FTE, audit variant 28"),
    Moment("gap_teachers_pp", "gap_teachers_k12", "target", "log",
           "K-12 FTE teachers per K-12 pupil, audit FTE variant 11"),
    Moment("gap_income", "gap_pm_ftyr_all", "target", "log",
           "log pooled median earnings, full-time year-round workers (`pm_ftyr_all`)"),
    Moment("pm_earn_all", "gap_pm_earn_all", "target", "log",
           "log pooled median earnings, all workers with earnings (`pm_earn_all`)"),
    Moment("poverty_gap", "gap_pov_rate_5_17", "validation", "diff",
           "child poverty rate 5-17 (SAIPE), difference in levels"),
    Moment("seda_growth_gap", "gap_learn_rate", "validation", "diff",
           "SEDA learning rate (`learn_rate`), difference"),
    Moment("score_gap", "gap_score_mean", "validation", "diff",
           "SEDA mean score (`score_mean`), difference"),
    Moment("seda_ses", "gap_seda_ses", "validation", "diff",
           "SEDA composite SES, difference"),
    Moment("seda_lninc", "gap_seda_lninc", "validation", "diff",
           "log median household income (SEDA), difference"),
    Moment("rev_pp_local", "gap_rev_pp_local", "diagnostic", "log",
           "local revenue per pupil"),
    Moment("rev_pp_state", "gap_rev_pp_state", "diagnostic", "log",
           "state revenue per pupil"),
    Moment("rev_pp_fed", "gap_rev_pp_fed", "diagnostic", "log",
           "federal revenue per pupil"),
    Moment("rev_pp_total", "gap_rev_pp_total", "diagnostic", "log",
           "total revenue per pupil"),
    Moment("exp_pp", "gap_exp_pp", "diagnostic", "log",
           "current expenditure per pupil"),
    Moment("payroll_pp_k12", "gap_payroll_k12", "diagnostic", "log",
           "audited teacher wage bill per K-12 pupil"),
    Moment("share_local", "gap_share_local", "diagnostic", "diff",
           "local share of revenue, difference in levels"),
    Moment("share_state", "gap_share_state", "diagnostic", "diff",
           "state share of revenue, difference in levels"),
    Moment("local_revenue_share", "", "diagnostic", "level",
           "local share of revenue, LEVEL (pupil-weighted over the partition's districts)"),
    Moment("salary_tracked", "gap_salary_real", "memo", "log",
           "tracked `salary_real` gap (zero-coded wage bills kept, NYC absent, pre-K in FTE)"),
    Moment("teachers_pp_tracked", "gap_teachers_pp", "memo", "log",
           "tracked FTE-per-pupil gap (pre-K teachers and pupils included)"),
]
USD_PARTS = {"rev_pp_local": "local revenue", "rev_pp_state": "state revenue",
             "rev_pp_fed": "federal revenue", "rev_pp_total": "total revenue",
             "exp_pp": "current expenditure", "payroll_k12": "teacher payroll (K-12 pupils)"}
MOMENTS += [Moment(f"usd_{k}", f"dlev_{k}", "usd", "usd", f"{v}, $ per pupil, 2 minus 1")
            for k, v in USD_PARTS.items()]
BY_KEY = {m.key: m for m in MOMENTS}


def settings_gaps() -> list["de.Gap"]:
    """`data_estimate.Gap`s computed in every CZ table: the tracked ones the moments
    use (or replace), plus the audited definitions."""
    keep = ("pupils", "teachers_pp", "salary_real", "exp_pp", "rev_pp_local",
            "rev_pp_state", "rev_pp_total", "share_local", "share_state",
            "pov_rate_5_17", "learn_rate", "score_mean", "seda_ses", "seda_lninc")
    gaps = [g for g in de.GAPS if g.key in keep]
    gaps += [
        de.Gap("salary_audited", "audited salary per K-12 FTE", "κ_2/κ_1", "mean",
               col="salary_audited", weight="teachers", contrast="log"),
        de.Gap("teachers_k12", "K-12 FTE per K-12 pupil", "H̃_T,l / M_l", "ratio",
               num="teachers_k12", den="pupils_k12", contrast="log"),
        de.Gap("payroll_k12", "audited wage bill per K-12 pupil", "W_l / M_l", "ratio",
               num="wb_hybrid_pos", den="pupils_k12", contrast="log"),
        de.Gap("rev_pp_fed", "federal revenue per pupil", "federal transfer", "ratio",
               num="rev_fed_total", den="pupils", contrast="log"),
    ]
    return gaps


# --------------------------------------------------------------------------
# The audited panel (district_audit.py's corrections, on the base-year panel)
# --------------------------------------------------------------------------

def build_audited_panel(indir: Path) -> dict[str, Any]:
    """Base panel + raw F-33/CCD fields + ACS earnings columns, NYC collapsed to
    one district, and the audited columns.

    Everything follows `district_audit.main` step for step: the hybrid share is the
    median of teacher salaries / instruction salaries among the non-zero districts
    of the locale partition (0.881); the collapsed panel carries the audited salary
    and the wage bill, while the screen's quantiles are the panel-wide ones.
    """
    P, f, d = da.load_inputs(indir)
    store = de.RawStore(indir / "raw")
    cw_all = de.build_cwift(store)
    Q = da.add_raw_fields(P, f, d)
    Q = db.load_income_columns(Q, indir / "raw")
    z = de.partition(Q, "locale", MIN_PUPILS, MIN_GROUP_SHARE)
    nz = z[z.wagebill_teachers > 0]
    hyb = float((nz.wagebill_teachers / nz.salaries_instruction).median())
    A = da.derive(da.collapse_nyc(Q, f, d, cw_all), hybrid_share=hyb)
    add_audited_columns(A)
    LOG.info("panel: %d regular districts (NYC collapsed from %d geographic districts), "
             "hybrid share %.4f", len(A), int(Q.lea_name.astype(str).str.upper()
                                              .str.startswith(da.NYC_PREFIX).sum()), hyb)
    return dict(A=A, Q=Q, f=f, d=d, store=store, cw_all=cw_all, hybrid_share=hyb)


def add_audited_columns(A: pd.DataFrame) -> None:
    """In place: the columns the audited moments read."""
    prek_t = A.teachers_prek_fte.fillna(0)
    tkn = (A.teachers - prek_t).where(A.teachers > prek_t)
    # variant 28: hybrid wage bill / K-12 FTE, deflated, screened panel-wide at 0.5/99.5
    A["salary_audited"] = da.screen((A.wb_hybrid.where(A.wb_hybrid > 0) / tkn) / A.cwift,
                                    SCREEN_TAIL)
    # FTE variant 11
    A["teachers_k12"] = tkn
    pk12 = A.pupils - A.prek_enr
    A["pupils_k12"] = pk12.where(pk12 > 0)
    # teacher payroll per K-12 pupil: audited wage bill, unscreened ratio of totals
    A["wb_hybrid_pos"] = A.wb_hybrid.where(A.wb_hybrid > 0)
    # district-level versions of the per-pupil variables for the offset regressions;
    # data_estimate screens its own per-pupil ratios at the same within-year tails
    A["payroll_pp_k12"] = da.screen(A.wb_hybrid_pos / A.pupils_k12, SCREEN_TAIL)
    A["rev_pp_fed_s"] = da.screen(A.rev_pp_fed, SCREEN_TAIL)
    for c in ("payroll_pp_k12", "rev_pp_fed_s"):
        A[f"log_{c}"] = np.log(A[c].where(A[c] > 0))


def fill_nyc(A: pd.DataFrame, store: "de.RawStore", f: pd.DataFrame, indir: Path
             ) -> pd.DataFrame:
    """Sensitivity panel: the collapsed NYC row with every record the Chancellor's
    Office id carries -- F-33 revenue and expenditure, SAIPE, SEDA, ACS earnings --
    instead of the missing values the 32 geographic districts have.  Per-pupil
    ratios are screened with the cuts the panel's own [0.5%, 99.5%] screen used."""
    X = A.copy()
    i = X.index[X.leaid == NYC]
    if len(i) != 1:
        raise SystemExit("expected exactly one collapsed NYC row")
    i = i[0]
    fr = f[f.leaid == NYC].iloc[0]
    for c in ("rev_total", "rev_fed_total", "rev_state_total", "rev_local_total",
              "rev_local_prop_tax", "exp_current_elsec_total",
              "exp_current_instruction_total"):
        X.loc[i, c] = de.num(pd.Series([fr[c]]), positive=c.startswith(("rev_total", "exp_"))).iloc[0]
    pupils = X.loc[i, "pupils"]
    # unscreened ratios of the others fix the screen cuts exactly as `assemble` did
    for col, num in (("rev_pp_total", "rev_total"), ("rev_pp_local", "rev_local_total"),
                     ("rev_pp_state", "rev_state_total"), ("exp_pp", "exp_current_elsec_total")):
        raw = A[num] / A.pupils
        lo, hi = raw.quantile(SCREEN_TAIL), raw.quantile(1 - SCREEN_TAIL)
        v = X.loc[i, num] / pupils
        X.loc[i, col] = v if lo <= v <= hi else np.nan
        X.loc[i, f"log_{col}"] = np.log(v) if v > 0 and lo <= v <= hi else np.nan
    X.loc[i, "rev_pp_fed"] = X.loc[i, "rev_fed_total"] / pupils
    for s in ("local", "state", "fed"):
        X.loc[i, f"share_{s}"] = X.loc[i, f"rev_{s}_total"] / X.loc[i, "rev_total"]
    sa = de.build_saipe(store, [YEAR])
    sa = sa[sa.leaid == NYC]
    for c in ("pop_total", "pop_5_17", "pop_5_17_poverty", "pov_rate_5_17"):
        if len(sa) and c in X:
            X.loc[i, c] = sa.iloc[0][c]
    sd = de.build_seda(store)
    sd = sd[sd.leaid == NYC] if sd is not None else pd.DataFrame()
    for c in ("learn_rate", "score_mean", "seda_ses", "seda_lninc"):
        if len(sd) and c in sd and c in X:
            X.loc[i, c] = sd.iloc[0][c]
    inc = db.load_income_columns(pd.DataFrame({"leaid": [NYC]}), indir / "raw")
    for c in inc.columns:
        if c != "leaid" and c in X:
            X.loc[i, c] = inc.iloc[0][c]
    add_audited_columns(X)
    return X


# --------------------------------------------------------------------------
# CZ tables, bootstrap
# --------------------------------------------------------------------------

def scheme_tables(A: pd.DataFrame, scheme: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Partitioned districts and the CZ x moment table of one scheme."""
    part, g = db.cz_gap_table(A, settings_gaps(), scheme, MIN_PUPILS, MIN_GROUP_SHARE)
    g, _ = db.pooled_median_gaps(part, g)
    for k in USD_PARTS:
        g[f"dlev_{k}"] = g[f"{k}_2"] - g[f"{k}_1"]
    LOG.info("%-8s %d CZs, %d districts, %.2f m pupils", scheme, len(g), len(part),
             part.pupils.sum() / 1e6)
    return part, g


def summarise(stat: dict[str, Any]) -> dict[str, Any]:
    d = stat["draws"]
    d = d[np.isfinite(d)]
    lo, hi = np.percentile(d, [2.5, 97.5])
    return dict(value=float(stat["est"]), se=float(np.std(d, ddof=1)),
                ci95_lo=float(lo), ci95_hi=float(hi), n_cz=int(stat["n"]))


def level_share_stat(part: pd.DataFrame, g: pd.DataFrame, counts: np.ndarray
                     ) -> dict[str, Any]:
    """Pupil-weighted mean of the districts' local revenue share, whole-CZ bootstrap
    (the `local_revenue_share` of district_bootstrap.py, for any partition)."""
    d = part[(part.year == YEAR) & part.share_local.notna() & (part.pupils > 0)]
    tab = (d.assign(num=d.pupils * d.share_local).groupby("czone")
           .agg(num=("num", "sum"), den=("pupils", "sum")))
    tab = tab.reindex(g.czone.to_numpy())
    den = tab.den.fillna(0).to_numpy(float)
    num = tab.num.fillna(0).to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        draws = (counts @ num) / (counts @ den)
    return dict(est=float(num.sum() / den.sum()), draws=draws, n=int((den > 0).sum()),
                n_districts=int(len(d)))


def moment_set(g: pd.DataFrame, part: pd.DataFrame, B: int, seed: int,
               moments: list[Moment] = MOMENTS) -> dict[str, Any]:
    """Every moment on the CZs of `g`: value, bootstrap SE, 95% interval, n CZs, and
    the draws.  One (B x n_CZ) multiplicity matrix is shared by all statistics."""
    counts = db.draw_counts(len(g), B, seed)
    w = g.pupils.to_numpy(float)
    stats: dict[str, Any] = {}
    draws: dict[str, np.ndarray] = {}
    for m in moments:
        if m.kind == "level":
            s = level_share_stat(part, g, counts)
        elif m.col in g.columns and np.isfinite(g[m.col].to_numpy(float)).any():
            s = db.boot_wmean(g[m.col].to_numpy(float), w, counts)
        else:
            continue
        stats[m.key] = dict(summarise(s), role=m.role, kind=m.kind, label=m.label)
        if "n_districts" in s:
            stats[m.key]["n_districts"] = s["n_districts"]
        draws[m.key] = s["draws"]
    D = pd.DataFrame(draws)
    return dict(stats=stats, draws=D, corr=D.corr(), n_cz_table=len(g),
                pupils=float(g.pupils.sum()))


def restrict(g: pd.DataFrame, part: pd.DataFrame, czs: set[float]
             ) -> tuple[pd.DataFrame, pd.DataFrame]:
    return (g[g.czone.isin(czs)].reset_index(drop=True), part[part.czone.isin(czs)])


def corr_matrix(D: pd.DataFrame, keys: list[str]) -> list[list[float]]:
    return D[keys].corr().to_numpy(float).tolist()


# --------------------------------------------------------------------------
# Coverage
# --------------------------------------------------------------------------

def coverage(part: pd.DataFrame, g: pd.DataFrame, A: pd.DataFrame) -> dict[str, Any]:
    """Size of a partition and the share of its pupils with each input."""
    pu = float(part.pupils.sum())
    share = lambda mask: float(part.pupils[mask].sum() / pu)
    nd = part.groupby("czone").size()
    nyc = part.leaid == NYC
    grp = part.groupby("loc_id").pupils.sum() / pu
    bins_ok = part[[c for c in part.columns if c.startswith("bin_ftyr_all_")]].notna().all(axis=1)
    sal_gap = g.gap_salary_audited.notna()
    wn = g.pupils.to_numpy(float) / g.pupils.sum()
    return dict(
        n_cz=int(len(g)), n_districts=int(len(part)), pupils=pu,
        n_eff_cz=float(1.0 / (wn ** 2).sum()), top5_weight_share=float(np.sort(wn)[::-1][:5].sum()),
        pupil_share_of_regular_districts=float(pu / A.pupils.sum()),
        districts_per_cz_median=float(nd.median()), districts_per_cz_mean=float(nd.mean()),
        share_cz_le3_districts=float((nd <= 3).mean()),
        pupil_share_loc1=float(grp.get(1, np.nan)), pupil_share_loc2=float(grp.get(2, np.nan)),
        nyc_pupils_in_partition=float(part.pupils[nyc].sum()),
        valid_salary_audited=share(part.salary_audited.notna()),
        valid_salary_tracked=share(part.salary_real > 0),
        valid_k12_fte=share(part.teachers_k12.notna() & part.pupils_k12.notna()),
        valid_earnings_bins=share(bins_ok),
        valid_poverty=share(part.pov_rate_5_17.notna()),
        valid_seda_growth=share(part.learn_rate.notna()),
        valid_seda_score=share(part.score_mean.notna()),
        valid_rev_local=share(part.rev_local_total.notna() & part.rev_total.notna()),
        n_cz_salary_gap=int(sal_gap.sum()),
        pupil_share_cz_salary_gap=float(g.pupils[sal_gap].sum() / g.pupils.sum()))


# --------------------------------------------------------------------------
# State-aid offsets
# --------------------------------------------------------------------------

REG_SPECS: list[tuple[str, str, str, str]] = [
    # (name, y, x, kind)
    ("state on local, $ per pupil", "rev_pp_state", "rev_pp_local", "levels"),
    ("state on local, logs", "log_rev_pp_state", "log_rev_pp_local", "logs"),
    ("federal on local, $ per pupil", "rev_pp_fed_s", "rev_pp_local", "levels"),
    ("total revenue on local, $ per pupil", "rev_pp_total", "rev_pp_local", "levels"),
    ("total revenue on local, logs", "log_rev_pp_total", "log_rev_pp_local", "logs"),
    ("current expenditure on local, $ per pupil", "exp_pp", "rev_pp_local", "levels"),
    ("current expenditure on local, logs", "log_exp_pp", "log_rev_pp_local", "logs"),
    ("teacher payroll per K-12 pupil on local, $", "payroll_pp_k12", "rev_pp_local", "levels"),
    ("teacher payroll per K-12 pupil on local, logs", "log_payroll_pp_k12", "log_rev_pp_local", "logs"),
]


def offset_regressions(part: pd.DataFrame, sample: str) -> list[dict[str, Any]]:
    """Within-CZ district regressions (CZ fixed effects, pupil weights, CZ-clustered
    SEs, `data_estimate.ols_fe`) of each revenue / spending aggregate on local
    revenue per pupil.  The slope in levels is dollars per dollar of local revenue."""
    d = part[part.year == YEAR].copy()
    d["cz"] = d.czone
    rows = []
    for name, y, x, kind in REG_SPECS:
        r = de.ols_fe(d, y, [x], fe="cz", w="pupils", cluster="cz")
        if r is None:
            continue
        rows.append(dict(sample=sample, spec=name, kind=kind, y=y, x=x,
                         coef=float(r["coef"][x]), se=float(r["se"][x]), n=int(r["n"]),
                         n_cz=int(r["n_clusters"]), r2_within=float(r["r2_within"])))
    return rows


def decomposition(g: pd.DataFrame, part: pd.DataFrame, B: int, seed: int) -> dict[str, Any]:
    """Dollar decomposition of a split's gaps: CZ-weighted mean of (level in
    location 2 - level in location 1) for each revenue / spending aggregate, and
    dollars per dollar of the local gap (bootstrap on one set of draws)."""
    ms = [BY_KEY[f"usd_{k}"] for k in USD_PARTS] + [BY_KEY[k] for k in (
        "rev_pp_local", "rev_pp_state", "rev_pp_fed", "rev_pp_total", "exp_pp",
        "payroll_pp_k12", "share_local", "share_state")]
    s = moment_set(g, part, B, seed, ms)
    D = s["draws"]
    out: dict[str, Any] = dict(stats=s["stats"], per_dollar={}, n_cz=s["n_cz_table"])
    base = s["stats"]["usd_rev_pp_local"]["value"]
    out["local_gap_t"] = float(base / s["stats"]["usd_rev_pp_local"]["se"])
    for k in USD_PARTS:
        if k == "rev_pp_local":
            continue
        num = s["stats"][f"usd_{k}"]["value"]
        ratio = D[f"usd_{k}"] / D["usd_rev_pp_local"]
        r = ratio[np.isfinite(ratio)]
        out["per_dollar"][k] = dict(value=float(num / base), se=float(np.std(r, ddof=1)),
                                    ci95_lo=float(np.percentile(r, 2.5)),
                                    ci95_hi=float(np.percentile(r, 97.5)))
    return out


# --------------------------------------------------------------------------
# Hold-out split
# --------------------------------------------------------------------------

def holdout_split(g: pd.DataFrame, seed: int = HOLDOUT_SEED) -> dict[str, Any]:
    """Pre-registered split of the locale CZs.  Sort by partitioned pupils
    (descending, ties by ascending CZ id), pair consecutive CZs (1-2, 3-4, ...; an unpaired last
    CZ goes to A); for each pair one `rng.integers(2)` from `default_rng(seed)`, in pair
    order: 0 sends the first CZ of the pair to hold-out B, 1 sends the second; the
    other goes to A."""
    cz = g.czone.to_numpy(float)
    order = np.lexsort((cz, -g.pupils.to_numpy(float)))
    ordered = cz[order]
    rng = np.random.default_rng(seed)
    a, b, pairs = [], [], []
    for k in range(len(ordered) // 2):
        first, second = float(ordered[2 * k]), float(ordered[2 * k + 1])
        r = int(rng.integers(2))
        to_b, to_a = (first, second) if r == 0 else (second, first)
        b.append(to_b)
        a.append(to_a)
        pairs.append(dict(pair=k + 1, first=first, second=second, draw=r, to_B=to_b))
    unpaired = float(ordered[-1]) if len(ordered) % 2 else None
    if unpaired is not None:
        a.append(unpaired)
    pu = np.sort(g.pupils.to_numpy(float))
    return dict(A=sorted(a), B=sorted(b), pairs=pairs, seed=seed, n=len(ordered),
                unpaired=unpaired, n_ties=int((np.diff(pu) == 0).sum()))


# --------------------------------------------------------------------------
# Sensitivity panels
# --------------------------------------------------------------------------

def k12_revenue_panel(A: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity panel: the revenue split variable is local revenue per K-12 pupil
    (pre-K enrollment out of the denominator, as in the FTE target) instead of per CCD
    pupil, screened at the same within-year tails.  Only the partition changes: the
    moments read totals, not this column."""
    X = A.copy()
    X["rev_pp_local"] = da.screen(X.rev_local_total / X.pupils_k12, SCREEN_TAIL)
    X["log_rev_pp_local"] = np.log(X.rev_pp_local.where(X.rev_pp_local > 0))
    return X


def reclassified(part0: pd.DataFrame, part1: pd.DataFrame) -> dict[str, float]:
    """How many pupils sit in a different location when the split variable changes."""
    m = (part0[["leaid", "loc_id", "pupils"]]
         .merge(part1[["leaid", "loc_id"]], on="leaid", how="inner", suffixes=("", "_v")))
    return dict(pupil_share_changed=float(m.pupils[m.loc_id != m.loc_id_v].sum() / m.pupils.sum()),
                districts_changed=int((m.loc_id != m.loc_id_v).sum()), districts_in_both=int(len(m)),
                pupil_share_of_baseline_in_both=float(m.pupils.sum() / part0.pupils.sum()))


# --------------------------------------------------------------------------
# Reproduction check
# --------------------------------------------------------------------------

# moment -> statistic of district_bootstrap.json that it must equal at the locale scheme
STORED_BOOT = {"gap_income": "pm_ftyr_all", "pm_earn_all": "pm_earn_all",
               "poverty_gap": "pov_rate_5_17", "seda_growth_gap": "learn_rate",
               "score_gap": "score_mean", "seda_ses": "seda_ses", "seda_lninc": "seda_lninc",
               "rev_pp_local": "rev_pp_local", "rev_pp_state": "rev_pp_state",
               "rev_pp_total": "rev_pp_total", "exp_pp": "exp_pp",
               "share_local": "share_local", "share_state": "share_state",
               "local_revenue_share": "local_revenue_share",
               "salary_tracked": "salary_real", "teachers_pp_tracked": "teachers_pp"}


# moment -> `gap_<key>_<scheme>` of spatial_moments.json (the tracked definitions)
TRACKED_MAP = {"gap_pupils": "pupils", "teachers_pp_tracked": "teachers_pp",
               "salary_tracked": "salary_real", "exp_pp": "exp_pp", "rev_pp_local": "rev_pp_local",
               "rev_pp_state": "rev_pp_state", "rev_pp_total": "rev_pp_total",
               "share_local": "share_local", "share_state": "share_state",
               "poverty_gap": "pov_rate_5_17", "seda_growth_gap": "learn_rate",
               "score_gap": "score_mean", "seda_ses": "seda_ses", "seda_lninc": "seda_lninc"}


def tracked_check(res: dict[str, dict[str, Any]], est_dir: Path) -> dict[str, Any]:
    """The tracked-definition moments of each scheme against spatial_moments.json
    (`data_estimate.py --scheme all`), value and number of CZs."""
    path = est_dir / "spatial_moments.json"
    if not path.exists():
        return {}
    tr = {e["key"]: e for e in json.loads(path.read_text())["estimates"]}
    out = {}
    for sch in SCHEMES:
        diffs, same_n = [], 0
        for k, tk in TRACKED_MAP.items():
            e = tr.get(f"gap_{tk}_{sch}")
            if e is None or k not in res[sch]["stats"]:
                continue
            diffs.append(abs(res[sch]["stats"][k]["value"] - e["value"]))
            same_n += int(res[sch]["stats"][k]["n_cz"] == e["n"])
        out[sch] = dict(n_compared=len(diffs), same_n=same_n,
                        max_abs_diff=float(max(diffs)) if diffs else None)
    return out


def reproduction(res: dict[str, dict[str, Any]], est_dir: Path, check_se: bool) -> dict[str, Any]:
    """The locale full-sample moments against (a) the digits the audit and bootstrap
    reports print for the four targets and (b) the full-precision values stored in
    district_audit.json / district_bootstrap.json for every moment that exists there;
    and the tracked-definition moments of all three schemes against spatial_moments.json."""
    stats = res["locale"]["stats"]
    stored: dict[str, tuple[float, float, int]] = {}
    ap, bp = est_dir / "district_audit.json", est_dir / "district_bootstrap.json"
    if ap.exists():
        a = json.loads(ap.read_text())
        for k, r in (("gap_salary", da._find(a["salary"], "hybrid + NYC collapsed + FTE")),
                     ("gap_teachers_pp", da._find(a["fte"], "K-12 consistent")),
                     ("gap_pupils", next(x for x in a["other"] if x["group"] == "pupils"))):
            stored[k] = (r["value"], r["se"], r["n"])
    if bp.exists():
        b = json.loads(bp.read_text())["statistics"]
        for k, bk in STORED_BOOT.items():
            stored[k] = (b[bk]["value"], b[bk]["boot_se"], b[bk]["n"])
    rows = []
    for k in list(REPORTED) + [k for k in STORED_BOOT if k not in REPORTED]:
        s = stats[k]
        row = dict(moment=k, value=s["value"], se=s["se"], n_cz=s["n_cz"])
        if k in REPORTED:
            row.update(reported_value=REPORTED[k][0], reported_se=REPORTED[k][1],
                       ok_value=round(s["value"], 4) == REPORTED[k][0],
                       ok_se=(round(s["se"], 4) == REPORTED[k][1]) if check_se else None)
        if k in stored:
            row.update(stored_value=stored[k][0], stored_se=stored[k][1], stored_n=stored[k][2],
                       diff_value=s["value"] - stored[k][0], diff_se=s["se"] - stored[k][1])
        rows.append(row)
    ok = all(r["ok_value"] and (r["ok_se"] in (True, None)) for r in rows if "ok_value" in r)
    diffs = [abs(r["diff_value"]) for r in rows if "diff_value" in r]
    sdiffs = [abs(r["diff_se"]) for r in rows if "diff_se" in r]
    return dict(rows=rows, ok=bool(ok), se_checked=bool(check_se),
                max_abs_diff_value=float(max(diffs)) if diffs else None,
                max_abs_diff_se=float(max(sdiffs)) if sdiffs else None,
                n_compared_stored=len(diffs), tracked_schemes=tracked_check(res, est_dir))


# --------------------------------------------------------------------------
# Serialisation
# --------------------------------------------------------------------------

def clean(o: Any) -> Any:
    """JSON-safe copy: numpy scalars to Python, NaN / inf to null."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (pd.Timestamp, dt.datetime)):
        return o.isoformat()
    return o


def toml_val(v: Any) -> str:
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if np.isnan(v):
            return "nan"
        return ("inf" if v > 0 else "-inf") if np.isinf(v) else repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    v = list(v)
    if v and isinstance(v[0], (list, tuple, np.ndarray)):
        return "[\n" + "".join(f"  {toml_val(r)},\n" for r in v) + "]"
    return "[" + ", ".join(toml_val(x) for x in v) + "]"


SET_LABELS = {
    "locale": "locale split (city vs suburb), the 188 CZs of the tracked moments",
    "revenue": "revenue split (below / above the CZ's pupil-weighted median of local revenue per pupil)",
    "ses": "ses split (below / above the CZ's pupil-weighted median SEDA SES)",
    "locale_common": "locale split on the CZs that are in all three partitions",
    "revenue_common": "revenue split on the CZs that are in all three partitions",
    "ses_common": "ses split on the CZs that are in all three partitions",
    "locale_holdout_A": "locale split, pre-registered hold-out half A (the sealed half B is not here)",
}
SET_ORDER = list(SET_LABELS)


def write_toml(path: Path, res: dict[str, dict[str, Any]], meta: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    a("# district_settings_targets.toml -- audited district moments by setting, for the Julia calibration.")
    a(f"# Generated {meta['generated']} by data/spatial/district_settings.py; see district_settings.md.")
    a("# value / se / n_cz per moment: weighted mean over CZs of the within-CZ gap (location 2 minus location 1),")
    a("# whole-CZ bootstrap SE (B = %d, seed %d), number of CZs with the statistic. corr = bootstrap correlation" %
      (meta["B"], meta["seed"]))
    a("# of the four targets, in corr_keys order (draws shared across moments within a table).")
    a("")
    a("[meta]")
    for k in ("generated", "script", "year", "B", "seed", "holdout_seed", "min_pupils",
              "min_group_share", "n_common_cz", "hybrid_share"):
        a(f"{k} = {toml_val(meta[k])}")
    a("")
    for name in SET_ORDER:
        s = res[name]
        a(f"[{name}]")
        a(f"description = {toml_val(SET_LABELS[name])}")
        a(f"n_cz = {toml_val(s['n_cz_table'])}")
        a(f"corr_keys = {toml_val(TARGET_KEYS)}")
        a(f"corr = {toml_val(corr_matrix(s['draws'], TARGET_KEYS))}")
        a("")
        for k in TOML_KEYS:
            m = s["stats"].get(k)
            a(f"[{name}.{k}]")
            a(f"value = {toml_val(m['value'] if m else float('nan'))}")
            a(f"se = {toml_val(m['se'] if m else float('nan'))}")
            a(f"n_cz = {toml_val(m['n_cz'] if m else 0)}")
            a("")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def _f(v: Any, d: int = 4) -> str:
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def _cell(m: dict[str, Any] | None, d: int | None = None) -> str:
    if not m:
        return ""
    d = (0 if m["kind"] == "usd" else 4) if d is None else d
    return f"{m['value']:.{d}f} ({m['se']:.{d}f})"


def _pct(v: Any, d: int = 1) -> str:
    return "" if v is None else f"{100 * v:.{d}f}%"


ROLE_ORDER = ["target", "validation", "diagnostic", "memo", "usd"]
ROLE_NAME = {"target": "target", "validation": "validation", "diagnostic": "diagnostic",
             "memo": "memo (tracked definition)", "usd": "$ per pupil"}


def moment_table(a: Callable[[str], None], stats: dict[str, Any]) -> None:
    a("| moment | what | role | value | boot SE | 95% interval | n CZs |")
    a("|---|---|---|---|---|---|---|")
    for role in ROLE_ORDER:
        for k, m in stats.items():
            if m["role"] != role:
                continue
            d = 0 if m["kind"] == "usd" else 4
            a(f"| `{k}` | {m['label']} | {ROLE_NAME[role]} | {_f(m['value'], d)} | {_f(m['se'], d)} | "
              f"[{_f(m['ci95_lo'], d)}, {_f(m['ci95_hi'], d)}] | {m['n_cz']} |")


def corr_table(a: Callable[[str], None], mat: list[list[float]], keys: list[str]) -> None:
    a("| | " + " | ".join(f"`{k}`" for k in keys) + " |")
    a("|---|" + "---|" * len(keys))
    for k, row in zip(keys, mat):
        a(f"| `{k}` | " + " | ".join(f"{v:.2f}".replace("-0.00", "0.00") for v in row) + " |")


def tracked_sentence(rp: dict[str, Any]) -> str:
    t = rp.get("tracked_schemes") or {}
    if not t:
        return ""
    return ("The tracked-definition moments of all three partitions equal `spatial_moments.json` "
            "(" + "; ".join(f"{k}: {v['n_compared']} moments, max |diff| {v['max_abs_diff']:.1e}, "
                             f"same n for {v['same_n']}" for k, v in t.items()) + ").")


def write_summary(a: Callable[[str], None], R: dict[str, Any]) -> None:
    S = R["sets"]
    loc, rev, ses = S["locale"], S["revenue"], S["ses"]
    lc, rc, sc = S["locale_common"], S["revenue_common"], S["ses_common"]
    rp, c = R["reproduction"], R["common"]
    mv = lambda s, k: s["moments"][k]["value"]
    a("## Summary\n")
    a(f"1. *Reproduction.* The locale full-sample targets match the adopted values: gap_salary "
      f"{_cell(loc['moments']['gap_salary'])}, gap_teachers_pp {_cell(loc['moments']['gap_teachers_pp'])}, "
      f"gap_pupils {_cell(loc['moments']['gap_pupils'])}, gap_income {_cell(loc['moments']['gap_income'])}. "
      f"{rp['n_compared_stored']} moments are also compared with the full-precision values stored in "
      f"`district_audit.json` / `district_bootstrap.json`: largest absolute difference "
      f"{rp['max_abs_diff_value']:.1e} in the value and {rp['max_abs_diff_se']:.1e} in the bootstrap SE"
      + ("" if rp["se_checked"] else " (a short run: B differs from the reported runs)")
      + ". " + tracked_sentence(rp))
    a(f"2. *Coverage.* The locale split needs a city and a suburb in the same CZ ({loc['n_cz']} CZs, "
      f"{loc['n_districts']:,} districts, {loc['pupils'] / 1e6:.1f} m pupils); the revenue and ses splits "
      f"need only two groups of districts ({rev['n_cz']} CZs, {rev['n_districts']:,} districts, "
      f"{rev['pupils'] / 1e6:.1f} m pupils; {ses['n_cz']} CZs, {ses['n_districts']:,} districts, "
      f"{ses['pupils'] / 1e6:.1f} m). {c['n']} CZs are in all three. They hold "
      f"{_pct(c['pupil_share']['locale'])} of the locale pupils but {_pct(c['pupil_share']['revenue'])} "
      f"and {_pct(c['pupil_share']['ses'])} of the revenue and ses pupils, so the full-sample revenue and "
      f"ses averages are weighted toward many smaller CZs (effective number of CZs "
      f"{loc['coverage']['n_eff_cz']:.0f} locale, {rev['coverage']['n_eff_cz']:.0f} revenue, "
      f"{ses['coverage']['n_eff_cz']:.0f} ses). New York City (32 geographic districts, "
      f"{loc['coverage']['nyc_pupils_in_partition']:,.0f} pupils) is in the locale partition and in neither "
      f"of the others: its split variables (F-33 revenue, SEDA SES) are filed under the Chancellor's "
      f"Office id and are missing for the geographic districts.")
    a(f"3. *The four targets, same {c['n']} CZs* (value, bootstrap SE):\n")
    a("| target | locale | revenue | ses |")
    a("|---|---|---|---|")
    for k in TARGET_KEYS:
        a(f"| `{k}` | {_cell(lc['moments'][k])} | {_cell(rc['moments'][k])} | {_cell(sc['moments'][k])} |")
    a("")
    a(f"   The pupil-weighted median split equalises enrollment within a CZ, so the enrollment gap is about "
      f"zero under revenue ({_f(mv(rc, 'gap_pupils'), 3)}) and ses ({_f(mv(sc, 'gap_pupils'), 3)}) against "
      f"{_f(mv(lc, 'gap_pupils'), 3)} for locale: it no longer measures a size difference. The earnings gap "
      f"is {_f(mv(lc, 'gap_income'), 3)} (locale), {_f(mv(rc, 'gap_income'), 3)} (revenue), "
      f"{_f(mv(sc, 'gap_income'), 3)} (ses); the poverty gap {_f(mv(lc, 'poverty_gap'), 3)}, "
      f"{_f(mv(rc, 'poverty_gap'), 3)}, {_f(mv(sc, 'poverty_gap'), 3)}; the SEDA growth gap "
      f"{_f(mv(lc, 'seda_growth_gap'), 4)}, {_f(mv(rc, 'seda_growth_gap'), 4)}, "
      f"{_f(mv(sc, 'seda_growth_gap'), 4)}.")
    da_ = R["state_aid"]["decomposition"]["revenue_common"]
    sd = da_["stats"]
    pdl = da_["per_dollar"]
    reg0 = next(r for r in R["state_aid"]["regressions"]
                if r["sample"] == "revenue" and r["spec"].startswith("state on local, $"))
    a(f"4. *State-aid offsets (revenue split, common CZs).* The local revenue gap is "
      f"${sd['usd_rev_pp_local']['value']:,.0f} per pupil (log {_f(sd['rev_pp_local']['value'], 3)}); state aid "
      f"offsets {_f(-pdl['rev_pp_state']['value'], 2)} of each local dollar (SE {_f(pdl['rev_pp_state']['se'], 2)}) and "
      f"federal aid {_f(-pdl['rev_pp_fed']['value'], 2)}, leaving a total-revenue gap of "
      f"${sd['usd_rev_pp_total']['value']:,.0f} (log {_f(sd['rev_pp_total']['value'], 3)}), a current-expenditure "
      f"gap of ${sd['usd_exp_pp']['value']:,.0f} and a teacher-payroll gap of "
      f"${sd['usd_payroll_k12']['value']:,.0f} per K-12 pupil (log {_f(sd['payroll_pp_k12']['value'], 3)}). The "
      f"within-CZ district regression gives a state-on-local slope of {_f(reg0['coef'], 2)} "
      f"(SE {_f(reg0['se'], 2)}, {reg0['n']:,} districts, {reg0['n_cz']} CZs).")
    nyc = R["sensitivity"]["nyc_filled"]
    k12 = R["sensitivity"]["k12_split"]
    a(f"5. *NYC and pre-K.* With NYC filled in as one fully observed district, the revenue-split enrollment gap "
      f"moves from {_f(rev['moments']['gap_pupils']['value'], 3)} to "
      f"{_f(nyc['revenue']['moments']['gap_pupils']['value'], 3)} (SE {_f(nyc['revenue']['moments']['gap_pupils']['se'], 3)}) "
      f"and the locale earnings gap from {_f(loc['moments']['gap_income']['value'], 3)} to "
      f"{_f(nyc['locale']['moments']['gap_income']['value'], 3)} (section 9). Splitting on local revenue per "
      f"K-12 rather than per CCD pupil moves {_pct(k12['reclass']['pupil_share_changed'])} of pupils to the other "
      f"location and changes the revenue-split targets by at most "
      f"{max(abs(k12['moments'][k]['value'] - rev['moments'][k]['value']) for k in TARGET_KEYS):.4f} (section 9).")
    ho = R["holdout"]
    a(f"6. *Hold-out.* The pre-registered split puts {ho['n_A']} of the {ho['n']} locale CZs in half A (values in "
      f"section 8) and {ho['n_B']} in the sealed half B; B's results are in `sealed/holdout_B.json` and are not "
      f"reported here.\n")


def write_md(path: Path, R: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    m = R["meta"]
    S = R["sets"]
    a("# District moments under three settings: locale, revenue and SES splits\n")
    a(f"Generated {m['generated']} by `district_settings.py` (year {m['year']}, whole-CZ bootstrap "
      f"B = {m['B']}, seed {m['seed']}). For each way of cutting a commuting zone in two it computes the "
      "audited moment set used for the locale calibration targets (`district_audit.md` variant 28 and FTE "
      "variant 11, `district_bootstrap.md` pooled earnings medians) with the same corrections, the validation "
      "moments, and a whole-CZ bootstrap of every one of them. The scripts it builds on are unchanged.\n")
    write_summary(a, R)

    a("## 1. Method\n")
    a("**Settings.** `data_estimate.partition` at its tracked sample rules (districts with >= 500 pupils in a CZ, "
      "each location >= 10% of the CZ's partitioned pupils). Location 2 is the advantaged group: suburb "
      "(locale), above the CZ's pupil-weighted median of local revenue per pupil (revenue), above the median "
      "SEDA SES (ses). The revenue and ses cut points are computed on the districts that have the split "
      "variable, after the panel's [0.5%, 99.5%] screen of `rev_pp_local`.\n")
    a("**Moments.** Each is the pupil-weighted mean across CZs (weights = the CZ's partitioned pupils) of a "
      "within-CZ gap, location 2 minus location 1; a CZ missing the statistic drops out of the mean (and from "
      "that statistic's `n CZs`). Log gaps are `log(x_2/x_1)` of location aggregates, level gaps are `x_2 - x_1`.\n")
    a("| moment | definition |")
    a("|---|---|")
    a("| `gap_pupils` | log ratio of location totals of CCD enrollment (pre-K included), the tracked definition |")
    a(f"| `gap_salary` | audit variant 28: teacher wage bill per K-12 FTE, divided by CWIFT, teacher-FTE-weighted mean "
      f"over a location's districts, [0.5%, 99.5%] screen panel-wide. Zero-coded F-33 wage bills (IL, AK, NM and "
      f"others) are filled with {R['inputs']['hybrid_share']:.3f} x instruction salaries (the median teacher / "
      f"instruction ratio of the locale sample, the same constant in every setting), the 32 NYC geographic "
      f"districts are one district with the Chancellor's Office F-33 record, and the FTE denominator excludes "
      f"pre-K teacher FTE |")
    a("| `gap_teachers_pp` | audit FTE variant 11: log ratio of location totals of (teachers - pre-K teacher FTE) / "
      "(pupils - pre-K enrollment) |")
    a("| `gap_income` | `pm_ftyr_all`: log of the location median of full-time year-round earnings of workers 16+, "
      "read off ACS 5-year 2014-18 B20005 bin counts summed over the location's districts, linear interpolation "
      "inside bins; `pm_earn_all` is the same from B20001 (all workers with earnings) |")
    a("| `poverty_gap` | difference in the SAIPE child poverty rate 5-17 (ratio of location totals), levels |")
    a("| `seda_growth_gap`, `score_gap` | SEDA 6.0 learning rate and mean score, enrollment-weighted district "
      "means, difference in levels (SEDA's published units) |")
    a("| `seda_ses`, `seda_lninc` | SEDA composite SES and log median household income, difference |")
    a("| `rev_pp_*`, `exp_pp` | log ratio of location totals of local, state, federal, total revenue and of "
      "current expenditure per CCD pupil; `rev_pp_fed` is F-33 federal revenue (total = local + state + federal "
      "holds exactly) |")
    a("| `payroll_pp_k12` | log ratio of location totals of the audited (hybrid) teacher wage bill per K-12 pupil, "
      "unscreened like the other per-pupil ratios, districts with a positive wage bill |")
    a("| `share_local`, `share_state` | difference in the location ratio of local (state) to total revenue |")
    a("| `local_revenue_share` | LEVEL: pupil-weighted mean of the districts' local revenue share over the "
      "partition's districts, as in `district_bootstrap.py` |")
    a("| `salary_tracked`, `teachers_pp_tracked` | memo: the tracked definitions the corrections replace |")
    a("| `usd_*` | memo: location 2 minus location 1 in dollars per pupil, which add up across sources |")
    a("")
    a(f"**Bootstrap.** Whole CZs resampled with replacement, B = {m['B']}, seed {m['seed']}, "
      "`numpy.random.default_rng`; one (B x n_CZ) multiplicity matrix per CZ set shared by all of its "
      "statistics, so the correlations are coherent. SE = SD of the draws; the interval is the 2.5-97.5 "
      "percentile. As in `district_bootstrap.py` the panel, screen, partition and CZ filters are held fixed "
      "and only CZ membership of the average is resampled; a restricted CZ set (common CZs, hold-out half) is "
      "resampled over its own CZs.\n")

    a("## 2. Reproduction check (locale, full sample)\n")
    rp = R["reproduction"]
    a("The four adopted targets must reproduce the digits printed in `district_audit.md` / "
      "`district_bootstrap.md`; every other moment that exists in `district_audit.json` / "
      "`district_bootstrap.json` is compared at full precision.\n")
    a("| moment | reported (value, SE) | this script | n CZs | value - stored | SE - stored |")
    a("|---|---|---|---|---|---|")
    for r in rp["rows"]:
        rep = f"{r['reported_value']:.4f} ({r['reported_se']:.4f})" if "reported_value" in r else ""
        a(f"| `{r['moment']}` | {rep} | {r['value']:.4f} ({r['se']:.4f}) | {r['n_cz']} | "
          f"{'' if 'diff_value' not in r else format(r['diff_value'], '+.1e')} | "
          f"{'' if 'diff_se' not in r else format(r['diff_se'], '+.1e')} |")
    a(f"\nAll four targets match to the reported digits: **{'yes' if rp['ok'] else 'NO'}** "
      f"(standard errors checked: {'yes' if rp['se_checked'] else 'no, B differs from the reported runs'}). "
      + tracked_sentence(rp) + "\n")

    a("## 3. Coverage\n")
    cov = {k: S[k]["coverage"] for k in SET_ORDER}
    names = ["locale", "revenue", "ses", "locale_common", "revenue_common", "ses_common"]
    a("| | " + " | ".join(f"{n.replace('_common', ', common CZs')}" for n in names) + " |")
    a("|---|" + "---|" * len(names))
    rows = [("CZs", "n_cz", "{:,}"), ("districts", "n_districts", "{:,}"),
            ("pupils (millions)", "pupils", "{:.2f}"),
            ("share of all regular-district pupils", "pupil_share_of_regular_districts", "%"),
            ("districts per CZ, median", "districts_per_cz_median", "{:.0f}"),
            ("CZs with <= 3 districts", "share_cz_le3_districts", "%"),
            ("effective number of CZs (1 / sum of squared weights)", "n_eff_cz", "{:.0f}"),
            ("weight on the 5 largest CZs", "top5_weight_share", "%"),
            ("pupils in location 2", "pupil_share_loc2", "%"),
            ("NYC pupils in the partition", "nyc_pupils_in_partition", "{:,.0f}"),
            ("pupils with a valid audited salary (variant 28)", "valid_salary_audited", "%"),
            ("pupils with a valid tracked salary (positive, screened)", "valid_salary_tracked", "%"),
            ("CZs with an audited salary gap", "n_cz_salary_gap", "{:,}"),
            ("... their share of the partition's pupils", "pupil_share_cz_salary_gap", "%"),
            ("pupils with valid K-12 FTE and K-12 pupils", "valid_k12_fte", "%"),
            ("pupils with ACS full-time year-round earnings bins", "valid_earnings_bins", "%"),
            ("pupils with a SAIPE poverty rate", "valid_poverty", "%"),
            ("pupils with a SEDA learning rate", "valid_seda_growth", "%"),
            ("pupils with a SEDA mean score", "valid_seda_score", "%"),
            ("pupils with F-33 local and total revenue", "valid_rev_local", "%")]
    for lab, key, fmt in rows:
        cells = []
        for n in names:
            v = cov[n].get(key)
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("")
            elif fmt == "%":
                cells.append(f"{100 * v:.1f}%")
            elif key == "pupils":
                cells.append(f"{v / 1e6:.2f}")
            else:
                cells.append(fmt.format(v))
        a(f"| {lab} | " + " | ".join(cells) + " |")
    ex = R["common"]["excluded"]
    a(f"\nThe {R['common']['n']} common CZs are those in all three partitions. " + " ".join(
        f"The locale partition has CZ {r['czone']} ({r['largest_district']}, {_pct(r['locale_weight'], 2)} of locale "
        f"weight) which is not in the {sch} partition: {_pct(r['pupil_share_with_split_variable'], 0)} of its "
        f"pupils in districts of >= 500 pupils have the split variable." for sch in ("revenue", "ses")
        for r in ex[sch]) + f" {R['common']['n_revenue_only']} revenue-only and {R['common']['n_ses_only']} "
      f"ses-only CZs are outside the locale partition (they lack a city group or a suburb group with >= 10% of "
      f"the CZ's partitioned pupils). "
      "Hybrid share (median teacher / instruction salaries, non-zero districts, NYC excluded) by partition: "
      + ", ".join(f"{k} {v:.3f}" for k, v in R["inputs"]["hybrid_share_by_scheme"].items())
      + f"; the {R['inputs']['hybrid_share']:.3f} of the audit is used everywhere.\n")

    a("## 4. Audited moments by setting (full samples)\n")
    for sch in SCHEMES:
        s = S[sch]
        a(f"### {sch}\n")
        a(f"{SET_LABELS[sch]}: {s['n_cz']} CZs, {s['n_districts']:,} districts, {s['pupils'] / 1e6:.2f} m pupils.\n")
        moment_table(a, s["moments"])
        a("")

    a("## 5. The same CZs\n")
    a(f"All three settings restricted to the {R['common']['n']} CZs in every partition (bootstrap over those "
      f"CZs, CZ weights unchanged). `full` columns repeat section 4.\n")
    keys = [k for k in S["locale"]["moments"] if S["locale"]["moments"][k]["role"] in
            ("target", "validation", "diagnostic", "memo")]
    a("| moment | locale, full | locale, common | revenue, full | revenue, common | ses, full | ses, common |")
    a("|---|---|---|---|---|---|---|")
    for k in keys:
        cells = [_cell(S[n]["moments"].get(k)) for n in ("locale", "locale_common", "revenue",
                                                          "revenue_common", "ses", "ses_common")]
        a(f"| `{k}` | " + " | ".join(cells) + " |")
    a("")

    a("## 6. Bootstrap correlation of the four targets\n")
    for n in ("locale", "revenue", "ses", "locale_common", "revenue_common", "ses_common"):
        a(f"**{n}**\n")
        corr_table(a, S[n]["corr_targets"], TARGET_KEYS)
        a("")
    a("The correlation among the targets and the validation moments (every statistic, draw by draw) is in "
      "`district_settings.json` (`corr`) and the draws in `district_settings_draws_<set>.csv`.\n")

    a("## 7. State-aid offsets under the revenue split\n")
    sa = R["state_aid"]
    a("### 7a. Decomposition of the gaps\n")
    a("Within-CZ gap between location 2 and location 1 (above- and below-median local revenue per pupil for the "
      "revenue split; suburb and city, high and low SES for the others) of each component, in logs (ratio of "
      "location totals) and in dollars per pupil (location 2 minus 1, CZ-weighted); the dollar gaps add up "
      "(total = local + state + federal). `per $ of local gap` is the component's dollar gap over the local "
      "dollar gap, so a negative number is the offset (n.m. where the local gap is within two SEs of zero).\n")
    a("| set | component | log gap | boot SE | $ per pupil | boot SE | per $ of local gap | boot SE |")
    a("|---|---|---|---|---|---|---|---|")
    comps = [("local revenue", "rev_pp_local", "rev_pp_local"), ("state revenue", "rev_pp_state", "rev_pp_state"),
             ("federal revenue", "rev_pp_fed", "rev_pp_fed"), ("total revenue", "rev_pp_total", "rev_pp_total"),
             ("current expenditure", "exp_pp", "exp_pp"), ("teacher payroll, K-12 pupils", "payroll_pp_k12", "payroll_k12")]
    for sname in ("revenue", "revenue_common", "locale", "ses"):
        dec = sa["decomposition"][sname]
        for lab, mk, uk in comps:
            lg, us = dec["stats"][mk], dec["stats"][f"usd_{uk}"]
            pd_ = dec["per_dollar"].get(uk)
            nm = abs(dec["local_gap_t"]) < 2
            a(f"| {sname} | {lab} | {_f(lg['value'], 3)} | {_f(lg['se'], 3)} | {us['value']:,.0f} | {us['se']:,.0f} | "
              f"{'1' if pd_ is None else 'n.m.' if nm else _f(pd_['value'], 3)} | "
              f"{'' if pd_ is None or nm else _f(pd_['se'], 3)} |")
    a("")
    a("### 7b. Within-CZ district regressions\n")
    a("Each aggregate per pupil on local revenue per pupil, district level, CZ fixed effects, pupil weights, "
      "CZ-clustered SEs (`data_estimate.ols_fe`); per-pupil ratios screened at the panel's [0.5%, 99.5%] "
      "within-year tails (state revenue and local revenue as `data_estimate` screens them; federal revenue and "
      "payroll screened here the same way). The slope in dollars is dollars of the left-hand side per extra "
      "dollar of local revenue per pupil: the state-on-local rows are the offset rate.\n")
    a("| sample | regression | slope | SE | districts | CZs | within R2 |")
    a("|---|---|---|---|---|---|---|")
    for r in sa["regressions"]:
        a(f"| {r['sample']} | {r['spec']} | {r['coef']:.3f} | {r['se']:.3f} | {r['n']:,} | {r['n_cz']} | "
          f"{r['r2_within']:.2f} |")
    a("")

    a("## 8. Pre-registered hold-out (locale)\n")
    ho = R["holdout"]
    a(f"Rule (fixed before any hold-out number was computed): sort the {ho['n']} locale CZs by partitioned pupils "
      f"(descending, ties by ascending CZ id; {ho['n_ties']} tied pairs), pair consecutive CZs (1-2, 3-4, ...; an unpaired last CZ goes to half A), "
      f"and for each pair take one `rng.integers(2)` from `numpy.random.default_rng({ho['seed']})`, in pair "
      "order: 0 sends the first CZ of the pair to half B, 1 sends the second; the other CZ goes to half A. "
      f"Half A has {ho['n_A']} CZs and half B {ho['n_B']}. **Half B's results are written only to "
      "`sealed/holdout_B.json` and are not printed or summarised in any log or report.** The weights of "
      "the two halves are public, so B's point estimates are implied by the full-sample and A values; "
      "the seal is procedural.\n")
    a("**Half A**\n")
    moment_table(a, S["locale_holdout_A"]["moments"])
    a("")
    a("**Bootstrap correlation of the four targets, half A**\n")
    corr_table(a, S["locale_holdout_A"]["corr_targets"], TARGET_KEYS)
    a("")
    a("**CZ assignment** (CZ id: half), in pair order of partitioned pupils:\n")
    a("```")
    items = [f"{int(p['first'])}:{'B' if p['to_B'] == p['first'] else 'A'}"
             f" {int(p['second'])}:{'B' if p['to_B'] == p['second'] else 'A'}" for p in ho["pairs"]]
    for i in range(0, len(items), 6):
        a("  ".join(items[i:i + 6]))
    if ho["n"] % 2:
        a(f"{int(ho['unpaired'])}:A")
    a("```\n")

    a("## 9. Sensitivities (not part of the targets)\n")
    nyc = R["sensitivity"]["nyc_filled"]
    a("### 9a. NYC as one fully observed district\n")
    a("The headline panel collapses NYC's 32 geographic districts into one row but keeps, for the revenue, "
      "SAIPE, SEDA and ACS moments, the geographic districts' missing values, because those records are filed "
      "under the Chancellor's Office id (`3620580`, agency type 3): F-33 revenue and expenditure, SAIPE, SEDA "
      "(`New York City Public Schools`) and the ACS school-district tables all have NYC in one record. This "
      "panel attaches all of them to the collapsed NYC row (per-pupil ratios screened with the panel's own "
      "cuts). The enrollment, FTE and salary moments are unchanged under locale; the others move as shown.\n")
    nv = R["sensitivity"]["nyc_attached_values"]
    a(f"NYC's attached values: {nv['pupils']:,.0f} CCD pupils; local revenue ${nv['rev_pp_local']:,.0f}, state "
      f"${nv['rev_pp_state']:,.0f} and total revenue ${nv['rev_pp_total']:,.0f} per pupil; SAIPE {nv['pop_5_17']:,.0f} "
      f"children 5-17, {_pct(nv['pov_rate_5_17'])} poor; SEDA SES {nv['seda_ses']:.3f}, growth {nv['learn_rate']:.4f}, "
      f"score {nv['score_mean']:.4f}; median earnings (B20002) ${nv['earn_all']:,.0f}. It sits in location "
      + ", ".join(f"{nyc[s]['nyc_loc_id']} ({s})" for s in SCHEMES) + ".\n")
    a("| setting | moment | headline | NYC filled | difference |")
    a("|---|---|---|---|---|")
    for sch in SCHEMES:
        for k in TOML_KEYS + ["pm_earn_all", "rev_pp_local", "rev_pp_state", "rev_pp_total"]:
            h, x = S[sch]["moments"].get(k), nyc[sch]["moments"].get(k)
            if h and x:
                a(f"| {sch} | `{k}` | {_cell(h)} | {_cell(x)} | {x['value'] - h['value']:+.4f} |")
    a("\nNYC pupils in the partition with NYC filled: "
      + ", ".join(f"{s} {nyc[s]['coverage']['nyc_pupils_in_partition']:,.0f}" for s in SCHEMES) + ".\n")
    k12 = R["sensitivity"]["k12_split"]
    a("### 9b. Revenue split on local revenue per K-12 pupil\n")
    a("The audited FTE and payroll moments use K-12 pupils, but the revenue split divides local revenue by "
      "CCD enrollment, which includes pre-K pupils. This variant splits on local revenue per K-12 pupil "
      "(`pupils - pre-K enrollment`), screened at the same tails; moments unchanged. "
      f"{_pct(k12['reclass']['pupil_share_changed'])} of the pupils in districts that are in both partitions "
      f"({k12['reclass']['districts_changed']} of {k12['reclass']['districts_in_both']:,} districts) switch "
      f"location; the new partition has {k12['n_cz']} CZs and {k12['n_districts']:,} districts "
      f"(baseline {S['revenue']['n_cz']} and {S['revenue']['n_districts']:,}).\n")
    a("| moment | per CCD pupil (headline) | per K-12 pupil | difference |")
    a("|---|---|---|---|")
    for k in TOML_KEYS + ["pm_earn_all", "rev_pp_local", "rev_pp_state", "rev_pp_total"]:
        h, x = S["revenue"]["moments"].get(k), k12["moments"].get(k)
        if h and x:
            a(f"| `{k}` | {_cell(h)} | {_cell(x)} | {x['value'] - h['value']:+.4f} |")
    a("")

    a("## 10. Judgment calls and caveats\n")
    for t in R["caveats"]:
        a(f"- {t}")
    a("\n## Files\n")
    a("- `district_settings.json`: everything above, machine readable (moments with SE, interval and n; the full "
      "bootstrap correlation matrix of every set; coverage; regressions; hold-out assignment and half A).")
    a("- `district_settings_targets.toml`: value / se / n_cz of the four targets and three validation moments, and "
      "the 4 x 4 target correlation, for each setting, the common-CZ versions and half A.")
    a("- `district_settings_draws_<set>.csv`: the bootstrap draws (one row per draw, one column per moment) for "
      "`locale`, `revenue`, `ses`, `locale_common`, `revenue_common`, `ses_common`, `locale_holdout_A`.")
    a("- `district_settings-<jobid>.out`: the SLURM log of the run; it carries no hold-out numbers.")
    a("- `sealed/holdout_B.json`: hold-out half B, written sealed. Nothing in the logs or in this report reads "
      "it back; open it only when the researchers decide to.")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Sealed hold-out half
# --------------------------------------------------------------------------

def write_sealed(path: Path, g: pd.DataFrame, part: pd.DataFrame, B: int, seed: int,
                 ho: dict[str, Any]) -> None:
    """Compute the full moment set on hold-out half B and write it to `path`.
    Nothing computed here is logged or returned."""
    s = moment_set(g, part, B, seed)
    payload = dict(
        sealed=True,
        note=("Hold-out half B of the pre-registered locale split (district_settings.py). Sealed until "
              "the researchers open it: no log, report or summary reads or prints these numbers."),
        created=dt.datetime.now().isoformat(timespec="seconds"),
        holdout_seed=ho["seed"], B=B, bootstrap_seed=seed, n_cz=int(len(g)),
        cz_ids=[int(c) for c in ho["B"]],
        moments=s["stats"], corr_targets=corr_matrix(s["draws"], TARGET_KEYS),
        corr=s["corr"].round(6).to_dict())
    path.write_text(json.dumps(clean(payload), indent=2))
    LOG.info("wrote %s (contents not printed)", path.name)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def excluded_cz_info(A: pd.DataFrame, g_loc: pd.DataFrame, czs: set[float], var: str
                     ) -> list[dict[str, Any]]:
    """Locale CZs missing from another partition: largest district, share of locale
    weight, and the share of the CZ's pupils (districts >= 500) with the split variable."""
    rows = []
    for c in sorted(czs):
        d = A[(A.czone == c) & (A.pupils >= MIN_PUPILS)]
        rows.append(dict(
            czone=int(c), districts=int(len(d)),
            largest_district=str(d.sort_values("pupils", ascending=False).lea_name.iloc[0]) if len(d) else "",
            locale_weight=float(g_loc.pupils[g_loc.czone == c].sum() / g_loc.pupils.sum()),
            pupil_share_with_split_variable=float(d.pupils[d[var].notna()].sum() / d.pupils.sum())
            if len(d) else float("nan")))
    return rows


def build_caveats(R: dict[str, Any]) -> list[str]:
    S, nyc, k12 = R["sets"], R["sensitivity"]["nyc_filled"], R["sensitivity"]["k12_split"]
    rev, ses, loc = S["revenue"], S["ses"], S["locale"]
    d = lambda sch, k: nyc[sch]["moments"][k]["value"] - S[sch]["moments"][k]["value"]
    hb = R["inputs"]["hybrid_share_by_scheme"]
    cv = {s: S[s]["coverage"] for s in SCHEMES}
    ex = R["common"]["excluded"]
    return [
        "**NYC.** F-33 revenue and expenditure, SAIPE, SEDA and the ACS school-district tables file New York "
        "City under the Chancellor's Office id (`3620580`), which the frame filter drops; the CCD frame has its "
        f"32 geographic districts ({loc['coverage']['nyc_pupils_in_partition']:,.0f} pupils), whose values for those "
        "sources are missing. Locale: NYC enters `gap_pupils`, `gap_teachers_pp` and, through the audit's collapse "
        "to one row with the Chancellor's wage bill and CWIFT, `gap_salary`, but not `gap_income`, the poverty and "
        "SEDA gaps or the revenue moments, so the New York CZ's city location in those is the CZ's other cities. "
        "Revenue and ses: the split variable (local revenue per pupil, SEDA SES) is missing for NYC with or without "
        "the collapse, so NYC is outside both partitions, the collapse changes nothing there, and the New York CZ "
        "is split over its other districts. Section 9a attaches NYC's own records: it moves the locale "
        f"`gap_income` by {d('locale', 'gap_income'):+.4f} and `poverty_gap` by {d('locale', 'poverty_gap'):+.4f}, and "
        "puts NYC into the revenue and ses partitions, moving the revenue-split `gap_pupils` by "
        f"{d('revenue', 'gap_pupils'):+.4f} (to {nyc['revenue']['moments']['gap_pupils']['value']:.4f}, SE "
        f"{nyc['revenue']['moments']['gap_pupils']['se']:.4f}) and `gap_teachers_pp` by "
        f"{d('revenue', 'gap_teachers_pp'):+.4f}. The headline numbers follow the instruction (tracked sample rules, "
        "the audit's corrections) and are not NYC-complete.",
        "**Pre-K corrections.** They are district-level data corrections (pre-K FTE out of the salary denominator, "
        "K-12 FTE over K-12 pupils) and apply unchanged inside each location under every split. The partition "
        "itself, the enrollment gap and the per-pupil revenue and expenditure moments still divide by CCD "
        "enrollment, which includes pre-K pupils, so the revenue split ranks districts by local revenue per "
        "CCD pupil while the FTE and payroll moments use K-12 pupils. Section 9b splits on revenue per K-12 pupil "
        f"instead: {_pct(k12['reclass']['pupil_share_changed'])} of pupils change location and the four revenue "
        "targets move by "
        + ("less than" if all(abs(k12['moments'][k]['value'] - rev['moments'][k]['value']) < rev['moments'][k]['se']
                              for k in TARGET_KEYS) else "up to more than")
        + " their bootstrap SEs "
        f"(largest change {max(abs(k12['moments'][k]['value'] - rev['moments'][k]['value']) for k in TARGET_KEYS):.4f}).",
        f"**Hybrid fill and NYC collapse.** One constant ({R['inputs']['hybrid_share']:.3f}, the locale sample's "
        "median ratio of teacher to instruction salaries) fills the zero-coded wage bills in every setting, as the "
        f"target is defined; the sample medians are {hb['locale']:.3f} (locale), {hb['revenue']:.3f} (revenue) and "
        f"{hb['ses']:.3f} (ses), NYC excluded. Without the fill and the collapse a tracked salary exists for "
        f"{_pct(cv['locale']['valid_salary_tracked'])} (locale), {_pct(cv['revenue']['valid_salary_tracked'])} (revenue) "
        f"and {_pct(cv['ses']['valid_salary_tracked'])} (ses) of the partitioned pupils; with them "
        f"{_pct(cv['locale']['valid_salary_audited'])}, {_pct(cv['revenue']['valid_salary_audited'])} and "
        f"{_pct(cv['ses']['valid_salary_audited'])}. The locale-sample constant is not re-estimated per setting.",
        f"**Common CZs.** The comparison set is the {R['common']['n']} CZs in all three partitions, so that one set "
        "serves all three settings. Pairwise with locale the common sets would have "
        f"{R['common']['pairwise_with_locale']['revenue']} (revenue) and {R['common']['pairwise_with_locale']['ses']} "
        "(ses) CZs, because two locale CZs are each missing from one other partition: "
        + "; ".join(f"CZ {r['czone']} ({r['largest_district']}) from the {sch} partition, {_pct(r['locale_weight'], 2)} "
                    "of the locale weight" for sch in ("revenue", "ses") for r in ex[sch])
        + f". Their weights are negligible: on the pairwise sets (`common.pairwise` in the JSON) the seven "
        f"TOML moments and `pm_earn_all` differ from the three-way ones by at most "
        f"{max(v[sch]['max_abs_diff_vs_three_way'] for v in R['common']['pairwise'].values() for sch in ('locale', 'revenue', 'ses') if sch in v):.4f}. "
        "The CZ weights are each setting's own partitioned pupils, so the common-CZ revenue and ses averages "
        "weight a CZ by how many pupils its revenue (ses) partition covers, which is more than its city + suburb pupils.",
        f"**Weights and effective sample.** Weights are partitioned pupils, so a few large metros dominate: the "
        f"effective number of CZs is {loc['coverage']['n_eff_cz']:.0f} (locale), {rev['coverage']['n_eff_cz']:.0f} "
        f"(revenue) and {ses['coverage']['n_eff_cz']:.0f} (ses), and {_pct(rev['coverage']['share_cz_le3_districts'], 0)} "
        "of the revenue CZs have three or fewer partitioned districts, where the median split is a comparison of one "
        "district with the rest. The revenue and ses locations are not geographic and include town and rural "
        "districts: location 1/2 means below/above the CZ median, wherever the districts lie.",
        "**Mechanical gaps.** Under the revenue split the local-revenue gap, and under the ses split the SES gap, "
        "are produced by the split and are not findings; the pupil-weighted median cut also equalises enrollment, "
        "so `gap_pupils` is about zero by construction there. They are reported as diagnostics, not targets.",
        "**Bootstrap scope.** As in `district_bootstrap.py`, only CZ membership of the average is resampled; the "
        "district panel, the [0.5%, 99.5%] screens, the revenue / SES cut points and the CZ filters are held "
        "fixed. A CZ missing a statistic drops out of that statistic's draws, so n varies by draw. The sets use "
        "separate multiplicity matrices (same seed, different n), so correlations are comparable within a set, "
        "not across sets.",
        "**Offsets are associations.** The state-on-local slopes are within-CZ district correlations; state aid "
        "formulas respond to local wealth by design, so they are the offset rate the model needs, not the causal "
        "effect of exogenous local revenue.",
        "**Seal.** B's weights and CZ ids are public, so its point estimates are implied by the full-sample and A "
        "numbers; the seal only keeps B out of the logs and reports.",
    ]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(HERE))
    ap.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    ap.add_argument("--B", type=int, default=db.B_DEFAULT)
    ap.add_argument("--seed", type=int, default=db.SEED)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates" / "settings"
    (outdir / "sealed").mkdir(parents=True, exist_ok=True)
    B, seed = args.B, args.seed
    full_run = B == db.B_DEFAULT and seed == db.SEED

    # ---- panel, partitions, full-sample moment sets ------------------------
    inp = build_audited_panel(indir)
    A = inp["A"]
    cz_tabs: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}     # name -> (CZ table, districts)
    res: dict[str, dict[str, Any]] = {}
    for s in SCHEMES:
        part, g = scheme_tables(A, s)
        cz_tabs[s] = (g, part)
        res[s] = moment_set(g, part, B, seed)

    # ---- reproduction: stop if the locale targets do not come back ----------
    rep = reproduction(res, indir / "estimates", check_se=full_run)
    for r in rep["rows"]:
        if "reported_value" in r:
            LOG.info("reproduction %-16s %.4f (%.4f)  reported %.4f (%.4f)  %s", r["moment"], r["value"],
                     r["se"], r["reported_value"], r["reported_se"],
                     ("ok" if r["ok_value"] and r["ok_se"] in (True, None) else "MISMATCH")
                     + ("" if full_run else " (value only: B differs from the reported runs)"))
    LOG.info("reproduction: %d moments vs stored files, max |diff| value %.1e, SE %.1e",
             rep["n_compared_stored"], rep["max_abs_diff_value"] or 0.0, rep["max_abs_diff_se"] or 0.0)
    for sch, v in rep["tracked_schemes"].items():
        LOG.info("reproduction vs spatial_moments.json, %-7s: %d moments, max |diff| %.1e, same n for %d",
                 sch, v["n_compared"], v["max_abs_diff"] or 0.0, v["same_n"])
    if not rep["ok"]:
        raise SystemExit("the locale targets do not reproduce the reported values; stopping "
                         "(see the reproduction lines above)")

    # ---- common CZs ---------------------------------------------------------
    czsets = {s: set(cz_tabs[s][0].czone) for s in SCHEMES}
    common = set.intersection(*czsets.values())
    for s in SCHEMES:
        g, part = cz_tabs[s]
        cz_tabs[f"{s}_common"] = restrict(g, part, common)
        res[f"{s}_common"] = moment_set(*cz_tabs[f"{s}_common"], B, seed)
    LOG.info("common CZs (in all three partitions): %d", len(common))
    pairwise: dict[str, Any] = {}          # the literal reading: revenue (ses) and locale on their own overlap
    for other in ("revenue", "ses"):
        cs = czsets["locale"] & czsets[other]
        pairwise[other] = dict(n_cz=len(cs), cz_not_in_three_way=sorted(int(c) for c in cs - common))
        for sch in ("locale", other):
            r = moment_set(*restrict(*cz_tabs[sch], cs), B, seed)
            pairwise[other][sch] = dict(
                moments={k: r["stats"][k] for k in TOML_KEYS + ["pm_earn_all"]},
                max_abs_diff_vs_three_way=float(max(
                    abs(r["stats"][k]["value"] - res[f"{sch}_common"]["stats"][k]["value"])
                    for k in TOML_KEYS + ["pm_earn_all"])))

    # ---- hold-out: A reported, B sealed -------------------------------------
    ho = holdout_split(cz_tabs["locale"][0])
    cz_tabs["locale_holdout_A"] = restrict(*cz_tabs["locale"], set(ho["A"]))
    res["locale_holdout_A"] = moment_set(*cz_tabs["locale_holdout_A"], B, seed)
    write_sealed(outdir / "sealed" / "holdout_B.json", *restrict(*cz_tabs["locale"], set(ho["B"])),
                 B, seed, ho)

    # ---- coverage, hybrid share by scheme ------------------------------------
    cov = {n: coverage(cz_tabs[n][1], cz_tabs[n][0], A) for n in SET_ORDER}
    hyb_by = {}
    for s in SCHEMES:
        p = cz_tabs[s][1]
        nz = p[(p.wagebill_teachers > 0) & (p.leaid != NYC)]
        hyb_by[s] = float((nz.wagebill_teachers / nz.salaries_instruction).median())

    # ---- state-aid offsets ---------------------------------------------------
    decs = {n: decomposition(cz_tabs[n][0], cz_tabs[n][1], B, seed)
            for n in ("revenue", "revenue_common", "locale", "ses")}
    regs: list[dict[str, Any]] = []
    for n in ("revenue", "revenue_common", "locale", "ses"):
        regs += offset_regressions(cz_tabs[n][1], n)

    # ---- sensitivities --------------------------------------------------------
    LOG.info("sensitivity: NYC as one fully observed district")
    X = fill_nyc(A, inp["store"], inp["f"], indir)
    nrow = X[X.leaid == NYC].iloc[0]
    nyc_attached = {k: float(nrow[k]) for k in (
        "pupils", "teachers", "rev_pp_local", "rev_pp_state", "rev_pp_total", "exp_pp", "share_local",
        "pop_5_17", "pov_rate_5_17", "learn_rate", "score_mean", "seda_ses", "seda_lninc", "earn_all")}
    nyc_res = {}
    for s in SCHEMES:
        part, g = scheme_tables(X, s)
        in_part = part.leaid == NYC
        nyc_res[s] = dict(moments=moment_set(g, part, B, seed)["stats"], n_cz=int(len(g)),
                          coverage=coverage(part, g, X),
                          nyc_loc_id=int(part.loc_id[in_part].iloc[0]) if in_part.any() else None)
    LOG.info("sensitivity: revenue split on local revenue per K-12 pupil")
    partK, gK = scheme_tables(k12_revenue_panel(A), "revenue")
    k12_res = dict(moments=moment_set(gK, partK, B, seed)["stats"], n_cz=int(len(gK)),
                   n_districts=int(len(partK)), reclass=reclassified(cz_tabs["revenue"][1], partK))

    # ---- assemble ------------------------------------------------------------
    nyc_rows = inp["Q"][inp["Q"].lea_name.astype(str).str.upper().str.startswith(da.NYC_PREFIX)]
    R: dict[str, Any] = dict(
        meta=dict(generated=dt.datetime.now().isoformat(timespec="seconds"),
                  script="data/spatial/district_settings.py", year=YEAR, B=B, seed=seed,
                  holdout_seed=HOLDOUT_SEED, min_pupils=MIN_PUPILS, min_group_share=MIN_GROUP_SHARE,
                  n_common_cz=len(common), hybrid_share=inp["hybrid_share"],
                  geo_unit=A.attrs.get("geo_unit"), full_run=full_run,
                  resampling="whole CZ with replacement; one multiplicity matrix per CZ set"),
        inputs=dict(hybrid_share=inp["hybrid_share"], hybrid_share_by_scheme=hyb_by,
                    nyc=dict(chancellor_id=NYC, n_geographic_districts=int(len(nyc_rows)),
                             pupils=float(nyc_rows.pupils.sum())),
                    sets=SET_LABELS),
        reproduction=rep, sets={}, state_aid=dict(decomposition=decs, regressions=regs),
        sensitivity=dict(nyc_filled=nyc_res, nyc_attached_values=nyc_attached, k12_split=k12_res))
    for n in SET_ORDER:
        r = res[n]
        R["sets"][n] = dict(label=SET_LABELS[n], n_cz=r["n_cz_table"], n_districts=int(len(cz_tabs[n][1])),
                            pupils=r["pupils"], coverage=cov[n], moments=r["stats"],
                            corr=r["corr"].round(6).to_dict(),
                            corr_targets=corr_matrix(r["draws"], TARGET_KEYS))
    g_loc = cz_tabs["locale"][0]
    R["common"] = dict(
        cz=sorted(int(c) for c in common), n=len(common),
        excluded=dict(revenue=excluded_cz_info(A, g_loc, czsets["locale"] - czsets["revenue"], "rev_pp_local"),
                      ses=excluded_cz_info(A, g_loc, czsets["locale"] - czsets["ses"], "seda_ses")),
        pupil_share={s: float(cz_tabs[f"{s}_common"][1].pupils.sum() / cz_tabs[s][1].pupils.sum())
                     for s in SCHEMES},
        locale_only_vs_revenue=sorted(czsets["locale"] - czsets["revenue"]),
        locale_only_vs_ses=sorted(czsets["locale"] - czsets["ses"]),
        n_revenue_only=len(czsets["revenue"] - czsets["locale"]),
        n_ses_only=len(czsets["ses"] - czsets["locale"]),
        pairwise_with_locale=dict(revenue=len(czsets["locale"] & czsets["revenue"]),
                                  ses=len(czsets["locale"] & czsets["ses"])),
        pairwise=pairwise)
    R["holdout"] = dict(
        seed=HOLDOUT_SEED, n=ho["n"], n_A=len(ho["A"]), n_B=len(ho["B"]), unpaired=ho["unpaired"],
        n_ties=ho["n_ties"],
        rule=holdout_split.__doc__.split("\n\n")[0].replace("\n    ", " "),
        assignment={str(int(c)): ("B" if c in set(ho["B"]) else "A")
                    for c in cz_tabs["locale"][0].czone},
        pairs=ho["pairs"], A_cz=[int(c) for c in ho["A"]], B_cz=[int(c) for c in ho["B"]])
    R["caveats"] = build_caveats(R)

    # ---- write ---------------------------------------------------------------
    (outdir / "district_settings.json").write_text(json.dumps(clean(R), indent=2))
    write_toml(outdir / "district_settings_targets.toml", res, R["meta"])
    write_md(outdir / "district_settings.md", R)
    for n in SET_ORDER:
        res[n]["draws"].round(10).to_csv(outdir / f"district_settings_draws_{n}.csv", index_label="draw")
    LOG.info("wrote district_settings.{json,md}, district_settings_targets.toml, %d draw files, "
             "sealed/holdout_B.json in %s", len(SET_ORDER), outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
