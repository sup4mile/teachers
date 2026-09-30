#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
district_audit.py -- audit of the pay and FTE-count moments (calibration item T3).

Rebuilds the 2018 locale-scheme salary gap (`gap_salary_real_locale`) and FTE-per-
pupil gap (`gap_teachers_pp_locale`) with `data_estimate.py`'s own functions, then
(1) documents how each is built and how many districts / pupils survive each step,
(2) checks whether numerator and denominator describe the same staff and the same
    districts (F-33 vs CCD), and
(3) recomputes each gap under reasonable alternatives, with whole-CZ bootstrap SEs
    (same machinery and seed as `district_bootstrap.py`).

Every variant changes the statistic, not the partition: a district that a variant
drops is blanked in that statistic (as the [0.5%, 99.5%] screen does), so the CZ
set, CZ weights (partitioned pupils) and group shares stay those of the baseline
unless the variant itself changes the sample filters (min pupils, group share).

Needs, beyond the `fetch_data.py` tables, two small Urban API pulls that it fetches
itself when absent: pre-K enrollment by district (raw/audit_urban_ccd_enrollment_prek_2018)
and the CCD directories of 2017 and 2019 (raw/audit_urban_ccd_directory_*), the latter
only for the F-33/CCD year-alignment check.

Usage
    python3 district_audit.py                 # B = 1999, seed 20260929 (about 4 minutes)
Writes estimates/district_audit.md (and district_audit.json).

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_estimate as de  # noqa: E402
import district_bootstrap as db  # noqa: E402

LOG = logging.getLogger("district_audit")
YEAR = db.YEAR
NYC_CHANCELLOR = "3620580"          # F-33 reports all of New York City under this id
NYC_PREFIX = "NEW YORK CITY GEOGRAPHIC DISTRICT"


# --------------------------------------------------------------------------
# Data: the baseline panel plus the raw F-33 / CCD fields the audit needs
# --------------------------------------------------------------------------

def load_inputs(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(panel as data_estimate builds it, raw F-33 table, raw CCD directory)."""
    P = db.build_panel(indir)
    raw = indir / "raw"
    f = pd.read_parquet(raw / "urban_ccd_finance.parquet")
    f["leaid"] = de.id7(f.leaid)
    f = f[f.year == YEAR].drop_duplicates("leaid")
    d = pd.read_parquet(raw / "urban_ccd_directory.parquet")
    d["leaid"] = de.id7(d.leaid)
    d = d[d.year == YEAR].drop_duplicates("leaid")
    return P, f, d


def add_raw_fields(P: pd.DataFrame, f: pd.DataFrame, d: pd.DataFrame) -> pd.DataFrame:
    """Merge the extra F-33 and CCD fields onto the panel (left join on leaid)."""
    fx = pd.DataFrame({"leaid": f.leaid})
    for c in ("salaries_instruction", "salaries_total", "benefits_employee_instruction",
              "benefits_employee_total", "enrollment_fall_school"):
        fx[c] = de.num(f[c]).to_numpy()
    fx["wb_regular"] = de.num(f["salaries_teachers_regular_prog"]).to_numpy()
    dx = pd.DataFrame({"leaid": d.leaid})
    for c in ("teachers_prek_fte", "instructional_aides_fte", "coordinators_fte",
              "staff_total_fte", "spec_ed_students"):
        dx[c] = de.num(d[c]).to_numpy()
    dx["teachers_ccd_raw"] = d["teachers_total_fte"].to_numpy(float)
    pk_path = Path(__file__).resolve().parent / "raw" / "audit_urban_ccd_enrollment_prek_2018.parquet"
    if not pk_path.exists():
        fetch_prek_enrollment(pk_path)
    pk = pd.read_parquet(pk_path)
    pk["leaid"] = de.id7(pk.leaid)
    dx = dx.merge(pk.drop_duplicates("leaid")[["leaid", "enrollment"]]
                  .rename(columns={"enrollment": "prek_enr"}), on="leaid", how="left")
    dx["prek_enr"] = de.num(dx.prek_enr).fillna(0)     # no pre-K row = no pre-K pupils
    attrs = dict(P.attrs)
    fx = fx[[c for c in fx.columns if c == "leaid" or c not in P.columns]]
    dx = dx[[c for c in dx.columns if c == "leaid" or c not in P.columns]]
    Q = P.merge(fx, on="leaid", how="left").merge(dx, on="leaid", how="left")
    Q.attrs.update(attrs)
    return Q


def fetch_prek_enrollment(path: Path) -> None:
    """Pre-K enrollment by district, Urban API CCD enrollment, grade-pk, all races and
    sexes (the directory's `enrollment` total includes pre-K pupils)."""
    import requests
    s = requests.Session()
    s.headers.update({"User-Agent": "teachers-spatial-calibration/1.0 "
                                    "(academic research; jpryan7@wisc.edu)"})
    url = ("https://educationdata.urban.org/api/v1/school-districts/ccd/enrollment/"
           f"{YEAR}/grade-pk/?race=99&sex=99")
    frames = []
    while url:
        r = s.get(url, timeout=180)
        r.raise_for_status()
        j = r.json()
        frames.append(pd.DataFrame(j["results"]))
        url = j.get("next")
    pd.concat(frames, ignore_index=True).to_parquet(path, index=False)


def fetch_directory(year: int, path: Path) -> None:
    """CCD directory for another year (Urban API), only the columns the year-alignment
    check needs; used to test that F-33 '2018' matches CCD 2018, not 2017 or 2019."""
    import requests
    s = requests.Session()
    s.headers.update({"User-Agent": "teachers-spatial-calibration/1.0 "
                                    "(academic research; jpryan7@wisc.edu)"})
    url = f"https://educationdata.urban.org/api/v1/school-districts/ccd/directory/{year}/"
    frames = []
    while url:
        r = s.get(url, timeout=180)
        r.raise_for_status()
        j = r.json()
        frames.append(pd.DataFrame(j["results"]))
        url = j.get("next")
    df = pd.concat(frames, ignore_index=True)
    keep = [c for c in ("year", "leaid", "lea_name", "fips", "county_code",
                        "urban_centric_locale", "enrollment", "teachers_total_fte",
                        "agency_type", "instructional_aides_fte", "staff_total_fte",
                        "number_of_schools", "agency_level") if c in df.columns]
    df[keep].to_parquet(path, index=False)


def screen(x: pd.Series, tail: float) -> pd.Series:
    """Blank values outside the panel-wide [tail, 1-tail] quantiles (as
    `data_estimate.screen_outliers`, which works on the whole panel)."""
    if tail <= 0:
        return x
    lo, hi = x.quantile(tail), x.quantile(1 - tail)
    return x.mask((x < lo) | (x > hi))


def collapse_nyc(Q: pd.DataFrame, f: pd.DataFrame, d: pd.DataFrame,
                 cw_all: pd.DataFrame | None) -> pd.DataFrame:
    """Replace the 32 CCD 'New York City Geographic District' rows (which carry
    -2 in F-33) by one NYC row: CCD pupils / staff summed over the geographic
    districts, F-33 fields from the Chancellor's Office record (`3620580`)."""
    is_nyc = Q.lea_name.astype(str).str.upper().str.startswith(NYC_PREFIX)
    nyc = Q[is_nyc]
    if nyc.empty or NYC_CHANCELLOR not in set(f.leaid):
        return Q
    fr = f[f.leaid == NYC_CHANCELLOR].iloc[0]
    row = nyc.sort_values("pupils", ascending=False).iloc[[0]].copy()
    row["leaid"] = NYC_CHANCELLOR
    row["lea_name"] = "NEW YORK CITY (32 geographic districts collapsed)"
    for c in ("pupils", "teachers", "teachers_prek_fte", "instructional_aides_fte",
              "coordinators_fte", "staff_total_fte", "teachers_ccd_raw", "prek_enr"):
        row[c] = nyc[c].sum(min_count=1)
    parts = ["salaries_teachers_regular_prog", "salaries_teachers_sped",
             "salaries_teachers_vocational", "salaries_teachers_other_ed"]
    row["wagebill_teachers"] = de.num(pd.Series([fr[c] for c in parts])).sum(min_count=1)
    for c in ("salaries_instruction", "salaries_total", "benefits_employee_instruction",
              "benefits_employee_total", "enrollment_fall_school"):
        row[c] = de.num(pd.Series([fr[c]])).iloc[0]
    row["enrollment_fall_responsible"] = float(fr["enrollment_fall_responsible"])
    cw = build_cwift_2018(cw_all)
    row["cwift"] = cw.get(NYC_CHANCELLOR, np.nan) if cw is not None else np.nan
    out = pd.concat([Q[~is_nyc], row], ignore_index=True)
    out.attrs.update(Q.attrs)
    return out


def build_cwift_2018(cw: pd.DataFrame | None) -> dict[str, float] | None:
    if cw is None:
        return None
    c = cw[cw.year == YEAR]
    return dict(zip(c.leaid, c.cwift))


def derive(Q: pd.DataFrame, tail: float = 0.005, hybrid_share: float | None = None
           ) -> pd.DataFrame:
    """Derived columns for the variants, from raw fields (not screened yet)."""
    Q = Q.copy()
    wb = Q.wagebill_teachers
    Q["wb_hybrid"] = wb.where(wb > 0, hybrid_share * Q.salaries_instruction) \
        if hybrid_share is not None else wb
    Q["one"] = 1.0
    return Q


# --------------------------------------------------------------------------
# Variant machinery
# --------------------------------------------------------------------------

class Runner:
    """Computes a CZ-level gap under a variant and its whole-CZ bootstrap SE."""

    def __init__(self, B: int, seed: int) -> None:
        self.B, self.seed = B, seed
        self._counts: dict[int, np.ndarray] = {}
        self.rows: list[dict[str, Any]] = []
        self.common: set[float] | None = None

    def counts(self, n: int) -> np.ndarray:
        if n not in self._counts:
            self._counts[n] = db.draw_counts(n, self.B, self.seed)
        return self._counts[n]

    def stat(self, g: pd.DataFrame, col: str, weight: str = "pupils") -> dict[str, Any]:
        """Weighted mean over CZs, its n and bootstrap SE."""
        x = g[col].to_numpy(float)
        w = (np.ones(len(g)) if weight == "unweighted" else g[weight].to_numpy(float))
        s = db.boot_wmean(x, w, self.counts(len(g)))
        d = s["draws"][np.isfinite(s["draws"])]
        return dict(value=s["est"], n=s["n"], se=float(np.std(d, ddof=1)))

    def run(self, group: str, label: str, panel: pd.DataFrame, gap: "de.Gap",
            note: str = "", scheme: str = "locale", min_pupils: int = 500,
            min_group_share: float = 0.10, weight: str = "pupils",
            restrict: Callable[[pd.DataFrame], pd.Series] | None = None
            ) -> dict[str, Any]:
        """One variant.  `value`/`n`/`se` use the CZs the variant itself yields;
        `value_c`/`n_c`/`se_c` restrict to the CZs of the baseline statistic
        (`self.common`), to separate a coverage change from a data change."""
        _, g = db.cz_gap_table(panel, [gap], scheme, min_pupils, min_group_share)
        col = f"gap_{gap.key}"
        if restrict is not None:
            g = g[restrict(g)].reset_index(drop=True)
        r = self.stat(g, col, weight)
        if self.common is not None:
            gc = g[g.czone.isin(self.common)].reset_index(drop=True)
            rc = self.stat(gc, col, weight)
            r.update(value_c=rc["value"], n_c=rc["n"], se_c=rc["se"])
        r.update(group=group, label=label, note=note)
        self.rows.append(r)
        return r


def mean_gap(key: str, col: str, weight: str = "teachers") -> "de.Gap":
    return de.Gap(key, key, "", "mean", col=col, weight=weight, contrast="log")


def ratio_gap(key: str, num: str, den: str) -> "de.Gap":
    return de.Gap(key, key, "", "ratio", num=num, den=den, contrast="log")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(HERE))
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--B", type=int, default=db.B_DEFAULT)
    ap.add_argument("--seed", type=int, default=db.SEED)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates"
    outdir.mkdir(parents=True, exist_ok=True)

    P, f, d = load_inputs(indir)
    store = de.RawStore(indir / "raw")
    cw_all = de.build_cwift(store)
    Q = add_raw_fields(P, f, d)
    R = Runner(args.B, args.seed)
    out: dict[str, Any] = {}

    # ---- baseline CZ table (also checks the reproduction) -------------------
    part, g0 = db.cz_gap_table(Q, list(de.GAPS))
    base = {k: de.wmean(g0[f"gap_{k}"], g0.pupils)
            for k in ("pupils", "salary_real", "teachers_pp", "salary_per_teacher",
                      "seda_lninc")}
    out["baseline"] = base
    LOG.info("baseline: %s", {k: round(v, 5) for k, v in base.items()})

    # ---- (1) attrition table -------------------------------------------------
    n_all_ccd, n_types12 = len(d), int(d.agency_type.isin(de.REGULAR_AGENCY_TYPES).sum())
    p_all = Q
    pu = lambda x: float(x.pupils.sum())
    z = part
    part_nf = de.partition(Q, "locale", 500, 0.0)
    attr = [
        ("CCD directory 2018, all agencies", n_all_ccd, float(d.enrollment.clip(lower=0).sum())),
        ("regular districts (agency type 1-2)", len(p_all), pu(p_all)),
        ("... with a commuting zone, >= 500 pupils, city or suburb locale, in a CZ "
         "that has both", int(len(part_nf)), pu(part_nf)),
        ("... after dropping CZs where a group has < 10% of pupils "
         "(the 188 CZs of the tracked moment)", int(len(z)), pu(z)),
        ("... with valid CCD FTE teachers (> 0)", int(z.teachers.notna().sum()),
         pu(z[z.teachers.notna()])),
        ("... with an F-33 teacher wage bill (non-missing)", int(z.wagebill_teachers.notna().sum()),
         pu(z[z.wagebill_teachers.notna()])),
        ("... and wage bill > 0 (F-33 zero-coded otherwise)",
         int((z.wagebill_teachers > 0).sum()), pu(z[z.wagebill_teachers > 0])),
        ("... and a CWIFT index", int(((z.wagebill_teachers > 0) & z.cwift.notna()).sum()),
         pu(z[(z.wagebill_teachers > 0) & z.cwift.notna()])),
        ("... surviving the [0.5%, 99.5%] screen = districts with `salary_real` > 0",
         int((z.salary_real > 0).sum()), pu(z[z.salary_real > 0])),
        ("`salary_real` non-missing including zeros (what enters the tracked gap)",
         int(z.salary_real.notna().sum()), pu(z[z.salary_real.notna()])),
    ]
    out["attrition"] = attr
    ncz_sal = int(g0.gap_salary_real.notna().sum())
    lost = g0[g0.gap_salary_real.isna()]
    nm = part.sort_values("pupils", ascending=False).groupby("czone").lea_name.first()
    out["lost_cz"] = [dict(czone=float(r.czone), name=str(nm.get(r.czone, "")),
                           pupils=float(r.pupils), weight=float(r.pupils / g0.pupils.sum()),
                           pupils_valid_salary_share=float(
                               part[(part.czone == r.czone) & (part.salary_real > 0)].pupils.sum()
                               / r.pupils))
                      for _, r in lost.sort_values("pupils", ascending=False).iterrows()]
    out["cz_salary"] = dict(n_cz=len(g0), n_with_gap=ncz_sal, n_lost=len(lost),
                            pupil_share_lost=float(lost.pupils.sum() / g0.pupils.sum()))

    # districts by group: coverage, city vs suburb
    cov = []
    for lg, gp in z.groupby("locale_grp"):
        cov.append(dict(group=lg, districts=len(gp), pupils=pu(gp),
                        teachers_valid=float(gp.teachers.notna().mean()),
                        wb_missing=float(gp.wagebill_teachers.isna().mean()),
                        wb_zero=float((gp.wagebill_teachers == 0).mean()),
                        wb_zero_pupils=float(gp[gp.wagebill_teachers == 0].pupils.sum() / gp.pupils.sum()),
                        cwift_missing=float(gp.cwift.isna().mean()),
                        salary_real_valid_pos=float((gp.salary_real > 0).mean())))
    out["coverage_by_group"] = cov

    # zero-coded states (all regular districts)
    Q["st"] = Q.leaid.str[:2]
    zc = (Q.assign(zero=Q.wagebill_teachers.eq(0))
          .groupby("st").agg(n=("leaid", "size"), zero=("zero", "sum"), pupils=("pupils", "sum")))
    zc["zero_share"] = zc.zero / zc.n
    out["zero_coded_states"] = zc[zc.zero_share >= 0.05].sort_values("zero_share", ascending=False) \
        .reset_index().to_dict("records")
    out["zero_all"] = dict(n=int(Q.wagebill_teachers.eq(0).sum()),
                           share=float(Q.wagebill_teachers.eq(0).mean()))
    raw_spt = Q.wagebill_teachers / Q.teachers
    out["screen_cuts"] = dict(
        nominal_lo=float(raw_spt.quantile(0.005)), nominal_hi=float(raw_spt.quantile(0.995)),
        real_lo=float((raw_spt / Q.cwift).quantile(0.005)),
        real_hi=float((raw_spt / Q.cwift).quantile(0.995)),
        share_zero_all=float((Q.wagebill_teachers == 0).mean()),
        n_screened_real=int(((raw_spt / Q.cwift).notna() & Q.salary_real.isna()).sum()))

    # NYC: the F-33 id mismatch
    is_nyc = Q.lea_name.astype(str).str.upper().str.startswith(NYC_PREFIX)
    out["nyc"] = dict(
        n_geo_districts=int(is_nyc.sum()), pupils=float(Q[is_nyc].pupils.sum()),
        teachers=float(Q[is_nyc].teachers.sum()),
        f33_wagebill_nan=int(Q[is_nyc].wagebill_teachers.isna().sum()),
        chancellor_wagebill=float(de.rowsum(f[f.leaid == NYC_CHANCELLOR],
                                            list(de.TEACHER_SALARY_PARTS)).iloc[0])
        if NYC_CHANCELLOR in set(f.leaid) else None,
        chancellor_enrollment=float(f[f.leaid == NYC_CHANCELLOR].enrollment_fall_responsible.iloc[0])
        if NYC_CHANCELLOR in set(f.leaid) else None,
        czone=float(Q[is_nyc].czone.mode().iloc[0]) if is_nyc.any() else None)

    # ---- (2) universe checks --------------------------------------------------
    unv: dict[str, Any] = {}
    nz = z[z.wagebill_teachers > 0]
    ratio_t_i = (nz.wagebill_teachers / nz.salaries_instruction).replace([np.inf, -np.inf], np.nan)
    unv["teacher_over_instruction_salary"] = ratio_t_i.describe(
        percentiles=[.01, .05, .5, .95, .99]).to_dict()
    unv["n_teacher_salary_above_instruction"] = int((ratio_t_i > 1.001).sum())
    r_enr = (z.pupils / z.enrollment_fall_responsible).replace([np.inf, -np.inf], np.nan)
    unv["ccd_over_f33_enrollment"] = dict(
        exact=float((r_enr == 1).mean()), within2pct=float(((r_enr - 1).abs() < 0.02).mean()),
        within10pct=float(((r_enr - 1).abs() < 0.10).mean()),
        n_off_by_over_10pct=int(((r_enr - 1).abs() >= 0.10).sum()))
    # year alignment of F-33 and CCD
    yr = {}
    for y, path in ((2017, "audit_urban_ccd_directory_2017.parquet"),
                    (2018, "urban_ccd_directory.parquet"),
                    (2019, "audit_urban_ccd_directory_2019.parquet")):
        pth = indir / "raw" / path
        if not pth.exists():
            LOG.info("fetching CCD directory %d for the year-alignment check", y)
            fetch_directory(y, pth)
        dd = pd.read_parquet(pth)
        dd["leaid"] = de.id7(dd.leaid)
        dd = dd[dd.agency_type.isin([1, 2, 1.0, 2.0])]
        m = dd.merge(f[["leaid", "enrollment_fall_responsible"]], on="leaid")
        m = m[(m.enrollment > 0) & (m.enrollment_fall_responsible > 0)]
        yr[y] = dict(n=len(m), exact=float((m.enrollment == m.enrollment_fall_responsible).mean()),
                     within2pct=float(((m.enrollment / m.enrollment_fall_responsible - 1).abs() < 0.02).mean()))
    unv["f33_enrollment_vs_ccd_year"] = yr
    # pre-K teachers and aides: do they load on pay per FTE within CZ?
    zz = nz.copy()
    zz["prek_share"] = zz.teachers_prek_fte / zz.teachers
    zz["aide_ratio"] = zz.instructional_aides_fte / zz.teachers
    zz["ln_spt"] = np.log(zz.wagebill_teachers / zz.teachers)
    zz["cz"] = zz.czone
    reg = de.ols_fe(zz[(zz.ln_spt > 0)], "ln_spt", ["prek_share", "aide_ratio"], fe="cz",
                    w="pupils", cluster="cz")
    unv["pay_on_prek_and_aides"] = reg
    # CCD teacher components sum
    unv["ccd_prek_teachers_share_mean"] = float((z.teachers_prek_fte / z.teachers).mean())
    out["universe"] = unv

    # ---- (3) sensitivity: salary ---------------------------------------------
    hyb = float((nz.wagebill_teachers / nz.salaries_instruction).median())
    out["hybrid_share"] = hyb
    A = derive(Q, hybrid_share=hyb)
    # unscreened building blocks on the whole panel
    A["spt_raw"] = A.wagebill_teachers / A.teachers
    A["real_raw"] = A.spt_raw / A.cwift
    A["wb_real"] = A.wagebill_teachers / A.cwift
    zero = A.wagebill_teachers.eq(0)
    A["spt_nz"] = A.spt_raw.where(~zero)
    A["real_nz"] = A.real_raw.where(~zero)
    A["real_nz_s"] = screen(A.real_nz, 0.005)          # trim among non-zero districts
    for tag, t in (("1", 0.01), ("25", 0.025), ("5", 0.05)):
        A[f"real_nz_s{tag}"] = screen(A.real_nz, t)
    A["real_noscreen"] = A.real_raw
    A["real_screen_only"] = A.salary_real                          # baseline column
    A["sal_real_nozero"] = A.salary_real.where(A.wagebill_teachers > 0)
    A["real_nz_range"] = A.real_nz.where((A.spt_nz >= 25_000) & (A.spt_nz <= 125_000))
    A["real_nz_lo25"] = A.real_nz.where(A.spt_nz >= 25_000)
    enr_ok = (A.pupils / A.enrollment_fall_responsible - 1).abs() <= 0.10
    A["real_enr_ok"] = A.real_nz.where(enr_ok)
    share_ti = A.wagebill_teachers / A.salaries_instruction
    A["real_share_ok"] = A.real_nz.where((share_ti >= 0.5) & (share_ti <= 1.0005))
    raw_str = A.pupils / A.teachers
    A["real_str_ok"] = A.real_nz.where((raw_str >= 6) & (raw_str <= 30))
    # alternative numerators / denominators, each screened at 0.5/99.5 across the panel
    def real_of(num: pd.Series, den: pd.Series, tail: float = 0.005) -> pd.Series:
        return screen((num / den) / A.cwift, tail)
    A["real_instr"] = real_of(A.salaries_instruction.where(A.salaries_instruction > 0), A.teachers)
    A["real_instr_aides"] = real_of(A.salaries_instruction.where(A.salaries_instruction > 0),
                                    A.teachers + A.instructional_aides_fte.fillna(0))
    A["real_total_staff"] = real_of(A.salaries_total.where(A.salaries_total > 0), A.staff_total_fte)
    A["real_regular"] = real_of(A.wb_regular.where(A.wb_regular > 0), A.teachers)
    comp = (A.wagebill_teachers + A.benefits_employee_instruction).where(A.wagebill_teachers > 0)
    A["real_comp"] = real_of(comp, A.teachers)
    A["real_hybrid"] = real_of(A.wb_hybrid.where(A.wb_hybrid > 0), A.teachers)
    A["nom_hybrid"] = screen(A.wb_hybrid.where(A.wb_hybrid > 0) / A.teachers, 0.005)
    A["nom_instr"] = screen(A.salaries_instruction.where(A.salaries_instruction > 0) / A.teachers, 0.005)
    tk = (A.teachers - A.teachers_prek_fte.fillna(0)).where(A.teachers > A.teachers_prek_fte.fillna(0))
    A["real_exprek"] = screen((A.wagebill_teachers.where(A.wagebill_teachers > 0) / tk) / A.cwift, 0.005)
    # deflator alternatives: other CWIFT releases
    for yy in (2015, 2019):
        cwy = dict(zip(*(cw_all[cw_all.year == yy][["leaid", "cwift"]].to_numpy().T)))
        A[f"cwift_{yy}"] = A.leaid.map(cwy)
        A[f"real_cw{yy}"] = screen(A.spt_nz / A[f"cwift_{yy}"], 0.005)

    S = "salary gap"
    R.common = set(g0[g0.gap_salary_real.notna()].czone)
    add = R.run
    add(S, "baseline: teacher-weighted mean of screened district `salary_real` "
           "(zeros kept; [0.5,99.5] screen)", A, mean_gap("s", "salary_real"))
    # aggregation / deflation
    add(S, "pupil-weighted (not teacher-weighted) mean within location", A,
        mean_gap("s", "salary_real", "pupils"))
    add(S, "unweighted district mean within location", A, mean_gap("s", "salary_real", "one"))
    add(S, "ratio of totals, deflated: sum(wage bill/CWIFT) / sum(FTE)", A,
        ratio_gap("s", "wb_real", "teachers"))
    add(S, "undeflated, same aggregation (`salary_per_teacher`, screened, teacher-weighted mean)",
        A, mean_gap("s", "salary_per_teacher"))
    add(S, "undeflated ratio of totals (tracked `gap_salary_per_teacher_locale`)", A,
        ratio_gap("s", "wagebill_teachers", "teachers"))
    # zeros
    add(S, "baseline with zero-coded wage bills set to missing", A,
        mean_gap("s", "sal_real_nozero"))
    add(S, "zeros missing, no trimming", A, mean_gap("s", "real_nz"))
    add(S, "zeros missing, trim [0.5%, 99.5%] recomputed on non-zero districts", A,
        mean_gap("s", "real_nz_s"))
    add(S, "zeros missing, trim [1%, 99%]", A, mean_gap("s", "real_nz_s1"))
    add(S, "zeros missing, trim [2.5%, 97.5%]", A, mean_gap("s", "real_nz_s25"))
    add(S, "zeros missing, trim [5%, 95%]", A, mean_gap("s", "real_nz_s5"))
    add(S, "zeros missing, drop nominal salary/FTE < $25,000 or > $125,000", A,
        mean_gap("s", "real_nz_range"))
    add(S, "zeros missing, drop nominal salary/FTE < $25,000 only", A,
        mean_gap("s", "real_nz_lo25"))
    add(S, "zeros missing, drop districts whose CCD and F-33 enrollment differ by > 10%", A,
        mean_gap("s", "real_enr_ok"))
    add(S, "zeros missing, drop teacher pay outside 50-100% of instruction salaries", A,
        mean_gap("s", "real_share_ok"))
    add(S, "zeros missing, drop students-per-teacher outside [6, 30]", A,
        mean_gap("s", "real_str_ok"))
    add(S, "zeros missing, FTE denominator excludes pre-K teacher FTE", A,
        mean_gap("s", "real_exprek"))
    add(S, "zeros missing, CWIFT 2015 release", A, mean_gap("s", "real_cw2015"))
    add(S, "zeros missing, CWIFT 2019 release", A, mean_gap("s", "real_cw2019"))
    # numerators / denominators
    add(S, "instruction salaries (F-33 `salaries_instruction`) per FTE teacher, deflated", A,
        mean_gap("s", "real_instr"))
    add(S, "instruction salaries per (FTE teachers + instructional aides), deflated", A,
        mean_gap("s", "real_instr_aides"))
    add(S, "total salaries per total staff FTE, deflated", A, mean_gap("s", "real_total_staff"))
    add(S, "regular-program teacher salaries only per FTE teacher, deflated", A,
        mean_gap("s", "real_regular"))
    add(S, "teacher salaries + instruction benefits per FTE teacher, deflated", A,
        mean_gap("s", "real_comp"))
    add(S, f"hybrid: teacher salaries, else {hyb:.2f} x instruction salaries where zero-coded "
           "(recovers IL, AK, NM)", A, mean_gap("s", "real_hybrid"))
    # NYC collapse
    Qn = collapse_nyc(Q, f, d, cw_all)
    An = derive(Qn, hybrid_share=hyb)
    An["spt_raw"] = An.wagebill_teachers / An.teachers
    An["real_raw"] = An.spt_raw / An.cwift
    znn = An.wagebill_teachers.eq(0)
    An["real_nz"] = An.real_raw.where(~znn)
    An["real_nz_s"] = screen(An.real_nz, 0.005)
    An["real_hybrid_raw"] = (An.wb_hybrid.where(An.wb_hybrid > 0) / An.teachers) / An.cwift
    An["real_hybrid"] = screen(An.real_hybrid_raw, 0.005)
    An["real_instr"] = screen((An.salaries_instruction.where(An.salaries_instruction > 0)
                               / An.teachers) / An.cwift, 0.005)
    add(S, "zeros missing + NYC collapsed to one district ([0.5,99.5] on non-zero)", An,
        mean_gap("s", "real_nz_s"))
    add(S, "hybrid + NYC collapsed", An, mean_gap("s", "real_hybrid"))
    tkn = (An.teachers - An.teachers_prek_fte.fillna(0)).where(An.teachers > An.teachers_prek_fte.fillna(0))
    An["real_hybrid_exprek"] = screen((An.wb_hybrid.where(An.wb_hybrid > 0) / tkn) / An.cwift, 0.005)
    add(S, "hybrid + NYC collapsed + FTE denominator excludes pre-K teacher FTE", An,
        mean_gap("s", "real_hybrid_exprek"))
    add(S, "instruction salaries per FTE teacher + NYC collapsed", An, mean_gap("s", "real_instr"))
    # sample rules (baseline statistic)
    for mp in (0, 250, 1000, 2000):
        add(S, f"min district size {mp} pupils (baseline 500)", A, mean_gap("s", "salary_real"),
            min_pupils=mp)
    for mg in (0.05, 0.20, 0.30):
        add(S, f"min group share {mg:.0%} of CZ pupils (baseline 10%)", A,
            mean_gap("s", "salary_real"), min_group_share=mg)
    # CZ aggregation
    add(S, "unweighted average across CZs", A, mean_gap("s", "salary_real"), weight="unweighted")
    for k in (1, 5, 10):
        big = set(g0.sort_values("pupils", ascending=False).czone.head(k))
        add(S, f"drop the {k} largest CZs by pupils", A, mean_gap("s", "salary_real"),
            restrict=lambda g, big=big: ~g.czone.isin(big))
    # partitions
    for sch in ("revenue", "ses"):
        add(S, f"partition: {sch} split (tracked: {'0.0154' if sch == 'revenue' else '0.0188'})",
            A, mean_gap("s", "salary_real"), scheme=sch)
    salary_rows = list(R.rows)
    R.rows = []

    # ---- (3b) sensitivity: FTE per pupil ------------------------------------
    T = "FTE per pupil gap"
    R.common = set(g0.czone)
    B_ = Q.copy()
    B_["one"] = 1.0
    B_["raw_str"] = B_.pupils / B_.teachers
    B_["tpp"] = B_.teachers / B_.pupils
    B_["teachers_s"] = B_.teachers.where(B_.teachers_pp.notna())                       # baseline screen
    for tag, t in (("1", 0.01), ("25", 0.025), ("5", 0.05)):
        B_[f"tpp_s{tag}"] = screen(B_.tpp, t)
        B_[f"teachers_s{tag}"] = B_.teachers.where(B_[f"tpp_s{tag}"].notna())
    B_["teachers_str"] = B_.teachers.where((B_.raw_str >= 6) & (B_.raw_str <= 30))
    B_["teachers_exprek"] = (B_.teachers - B_.teachers_prek_fte.fillna(0)).where(B_.teachers > 0)
    B_["instr_staff"] = (B_.teachers + B_.instructional_aides_fte).where(B_.teachers > 0)
    B_["teachers_f33ok"] = B_.teachers.where((B_.pupils / B_.enrollment_fall_responsible - 1).abs() <= 0.10)
    B_["teachers_salok"] = B_.teachers.where(B_.salary_real.notna() & (B_.wagebill_teachers > 0))
    B_["teachers_nzst"] = B_.teachers.where(~B_.leaid.str[:2].isin(["17", "02", "35"]))
    tk12 = B_.teachers - B_.teachers_prek_fte.fillna(0)
    B_["teachers_k12"] = tk12.where(tk12 > 0)
    pk12 = B_.pupils - B_.prek_enr
    B_["pupils_k12"] = pk12.where(pk12 > 0)
    B_["staff"] = B_.staff_total_fte
    B_["pupils_f33"] = B_.enrollment_fall_responsible.where(B_.enrollment_fall_responsible > 0)
    add(T, "baseline: ratio of location totals, sum(teachers)/sum(pupils) (screen has no effect)",
        B_, ratio_gap("t", "teachers", "pupils"))
    add(T, "apply the [0.5%, 99.5%] screen of `teachers_pp` (drop flagged districts)", B_,
        ratio_gap("t", "teachers_s", "pupils"))
    add(T, "screen [1%, 99%]", B_, ratio_gap("t", "teachers_s1", "pupils"))
    add(T, "screen [2.5%, 97.5%]", B_, ratio_gap("t", "teachers_s25", "pupils"))
    add(T, "screen [5%, 95%]", B_, ratio_gap("t", "teachers_s5", "pupils"))
    add(T, "drop students-per-teacher outside [6, 30]", B_, ratio_gap("t", "teachers_str", "pupils"))
    add(T, "mean of district ratios, teacher-weighted", B_, mean_gap("t", "tpp", "teachers"))
    add(T, "mean of district ratios, unweighted", B_, mean_gap("t", "tpp", "one"))
    add(T, "denominator: F-33 fall enrollment (responsible) instead of CCD enrollment", B_,
        ratio_gap("t", "teachers", "pupils_f33"))
    add(T, "drop districts whose CCD and F-33 enrollment differ by > 10%", B_,
        ratio_gap("t", "teachers_f33ok", "pupils"))
    add(T, "numerator excludes pre-K teacher FTE, denominator unchanged (pre-K pupils stay in: "
           "mismatched)", B_, ratio_gap("t", "teachers_exprek", "pupils"))
    add(T, "K-12 consistent: (teachers - pre-K teacher FTE) / (pupils - pre-K enrollment)", B_,
        ratio_gap("t", "teachers_k12", "pupils_k12"))
    add(T, "numerator: teachers + instructional aides", B_, ratio_gap("t", "instr_staff", "pupils"))
    add(T, "numerator: total staff FTE", B_, ratio_gap("t", "staff", "pupils"))
    add(T, "restrict to districts in the salary sample (valid, positive `salary_real`)", B_,
        ratio_gap("t", "teachers_salok", "pupils"))
    add(T, "drop IL, AK and NM (zero-coded F-33 states)", B_, ratio_gap("t", "teachers_nzst", "pupils"))
    for mp in (0, 250, 1000, 2000):
        add(T, f"min district size {mp} pupils (baseline 500)", B_,
            ratio_gap("t", "teachers", "pupils"), min_pupils=mp)
    for mg in (0.05, 0.20, 0.30):
        add(T, f"min group share {mg:.0%} of CZ pupils (baseline 10%)", B_,
            ratio_gap("t", "teachers", "pupils"), min_group_share=mg)
    add(T, "unweighted average across CZs", B_, ratio_gap("t", "teachers", "pupils"),
        weight="unweighted")
    for k in (1, 5, 10):
        big = set(g0.sort_values("pupils", ascending=False).czone.head(k))
        add(T, f"drop the {k} largest CZs by pupils", B_, ratio_gap("t", "teachers", "pupils"),
            restrict=lambda g, big=big: ~g.czone.isin(big))
    add(T, "restrict to the CZs that have a salary gap", B_, ratio_gap("t", "teachers", "pupils"),
        restrict=lambda g: g.czone.isin(g0[g0.gap_salary_real.notna()].czone))
    for sch in ("revenue", "ses"):
        add(T, f"partition: {sch} split (tracked: "
               f"{'0.0366' if sch == 'revenue' else '-0.0036'})", B_,
            ratio_gap("t", "teachers", "pupils"), scheme=sch)
    fte_rows = list(R.rows)
    R.rows = []

    # ---- (3c) other targets: sample-rule sensitivity ------------------------
    C = Q.copy()
    C["one"] = 1.0
    others = []
    for key, gap in (("pupils", de.Gap("p", "p", "", "total", num="pupils", contrast="log")),
                     ("seda_lninc", de.Gap("i", "i", "", "mean", col="seda_lninc", contrast="diff"))):
        R.rows = []
        R.common = set(g0[g0["gap_pupils" if key == "pupils" else "gap_seda_lninc"].notna()].czone)
        R.run(key, "baseline", C, gap)
        for mp in (0, 250, 1000, 2000):
            R.run(key, f"min district size {mp} pupils (baseline 500)", C, gap, min_pupils=mp)
        for mg in (0.05, 0.20, 0.30):
            R.run(key, f"min group share {mg:.0%} (baseline 10%)", C, gap, min_group_share=mg)
        R.run(key, "unweighted average across CZs", C, gap, weight="unweighted")
        for k in (1, 5, 10):
            big = set(g0.sort_values("pupils", ascending=False).czone.head(k))
            R.run(key, f"drop the {k} largest CZs by pupils", C, gap,
                  restrict=lambda g, big=big: ~g.czone.isin(big))
        R.run(key, "restrict to the CZs that have a salary gap", C, gap,
              restrict=lambda g: g.czone.isin(g0[g0.gap_salary_real.notna()].czone))
        for sch in ("revenue", "ses"):
            R.run(key, f"partition: {sch} split", C, gap, scheme=sch)
        others += R.rows
    R.rows = []

    # influence of the largest CZs on the four targets
    wts = g0.pupils / g0.pupils.sum()
    big5 = g0.assign(w=wts).sort_values("w", ascending=False).head(5)
    names = (part.sort_values("pupils", ascending=False).groupby("czone").lea_name.first())
    out["top_cz"] = [dict(czone=float(r.czone), name=str(names.get(r.czone, "")),
                          weight=float(r.w),
                          gaps={k: (None if not np.isfinite(r[f"gap_{k}"]) else float(r[f"gap_{k}"]))
                                for k in ("pupils", "salary_real", "teachers_pp", "seda_lninc")})
                     for _, r in big5.iterrows()]

    out.update(salary=salary_rows, fte=fte_rows, other=others,
               meta=dict(generated=dt.datetime.now().isoformat(timespec="seconds"),
                         B=args.B, seed=args.seed, year=YEAR))
    (outdir / "district_audit.json").write_text(json.dumps(out, indent=2, default=str))
    write_md(outdir / "district_audit.md", out)
    LOG.info("wrote district_audit.{md,json}")
    return 0


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def _f(v: Any, d: int = 4) -> str:
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def var_table(rows: list[dict[str, Any]], base_value: float, hdr: str = "variant") -> list[str]:
    has_c = any("value_c" in r for r in rows)
    L = [f"| # | {hdr} | gap | boot SE | n CZs | diff vs baseline |"
         + (" gap on baseline CZs | boot SE | n | diff |" if has_c else ""),
         "|---|---|---|---|---|---|" + ("---|---|---|---|" if has_c else "")]
    b_c = rows[0].get("value_c")
    for i, r in enumerate(rows):
        dv = "" if i == 0 else "%+.4f" % (r["value"] - base_value)
        line = f"| {i} | {r['label']} | {_f(r['value'])} | {_f(r['se'])} | {r['n']} | {dv} |"
        if has_c:
            dc = "" if i == 0 or "value_c" not in r else "%+.4f" % (r["value_c"] - b_c)
            line += (f" {_f(r.get('value_c'))} | {_f(r.get('se_c'))} | {r.get('n_c', '')} | {dc} |")
        L.append(line)
    return L


def _find(rows: list[dict[str, Any]], prefix: str) -> dict[str, Any]:
    for r in rows:
        if r["label"].startswith(prefix):
            return r
    raise KeyError(prefix)


def write_summary(L: list[str], o: dict[str, Any]) -> None:
    a = L.append
    S, T = o["salary"], o["fte"]
    b, bt = S[0], T[0]
    cz = o["cz_salary"]
    und = _find(S, "undeflated, same aggregation")
    unr = _find(S, "undeflated ratio of totals")
    pre = _find(S, "zeros missing, FTE denominator excludes pre-K")
    ins = _find(S, "instruction salaries (F-33")
    ben = _find(S, "teacher salaries + instruction benefits")
    reg = _find(S, "regular-program teacher salaries only")
    hyb = _find(S, "hybrid + NYC collapsed")
    hybp = _find(S, "hybrid + NYC collapsed + FTE")
    nyc = _find(S, "zeros missing + NYC collapsed")
    notrim = _find(S, "zeros missing, no trimming")
    uw = _find(S, "unweighted average across CZs")
    d10 = _find(S, "drop the 10 largest")
    dist_uw = _find(S, "unweighted district mean")
    k12 = _find(T, "K-12 consistent")
    tuw = _find(T, "unweighted average across CZs")
    t5 = _find(T, "drop the 5 largest")
    tsal = _find(T, "restrict to the CZs that have a salary gap")
    taid = _find(T, "numerator: teachers + instructional aides")
    a("## Summary\n")
    a("**Salary gap** (`gap_salary_real_locale` = "
      f"{b['value']:.4f}, whole-CZ bootstrap SE {b['se']:.4f}, i.e. {b['value'] / b['se']:.1f} SEs "
      f"from zero; {b['n']} CZs).\n")
    a(f"1. *Coverage.* It averages over {100 * (1 - cz['pupil_share_lost']):.1f}% of the "
      f"partitioned pupils: {cz['n_lost']} CZs drop out, among them New York (F-33 files all of "
      "NYC under the Chancellor's Office, which the frame filter removes; the 32 geographic "
      "districts are `-2` in F-33), Chicago and the other Illinois CZs, Albuquerque and Alaska "
      "(F-33 reports zero in all four teacher-salary fields for those states; 8.3% of all regular "
      "districts). The other three targets use all 188 CZs.")
    pref = ("pupil-weighted", "zeros missing, trim", "zeros missing, drop nominal",
            "zeros missing, CWIFT", "zeros missing, FTE denominator", "teacher salaries + instruction",
            "hybrid", "min district size", "min group share", "drop the", "baseline with zero",
            "instruction salaries per (FTE")
    cand = [(r.get("value_c", r["value"]), r["label"]) for r in S if r["label"].startswith(pref)]
    lo_, hi_ = min(cand), max(cand)
    a(f"2. *Sign and size.* The positive gap comes entirely from the CWIFT deflator: undeflated "
      f"it is {und['value']:.4f} (teacher-weighted mean) and {unr.get('value_c', np.nan):.4f} "
      "(ratio of totals, same 176 CZs). Across the trimming, deflator-vintage, denominator, "
      "sample-rule and CZ-dropping variants (same 176 CZs) it ranges from "
      f"{lo_[0]:.4f} ({lo_[1]}) to {hi_[0]:.4f} ({hi_[1]}); it is {dist_uw['value']:.4f} "
      f"with unweighted district means and {uw['value']:.4f} with unweighted CZ averaging. "
      f"Treat it as a small gap, imprecisely measured (bootstrap SE {b['se']:.4f}).")
    a(f"3. *Trimming.* Because > 5% of the panel is exactly zero, the lower [0.5%, 99.5%] cut is 0 "
      "and only the upper tail is trimmed. The upper cut "
      f"(${o['screen_cuts']['real_hi'] / 1000:.0f}k deflated) removes genuine "
      "high-pay districts (Long Island, Westchester, Greenwich, Bucks County). On the "
      f"common 176 CZs this changes the gap by {notrim['value_c'] - b['value_c']:+.4f}; but "
      "without the trim one more CZ (New York, weight 5.3%) enters with a spurious city group "
      "(White Plains and other Westchester cities, no NYC) and a gap of -0.23, moving the "
      f"aggregate to {notrim['value']:.4f}. The screen is thus what keeps that CZ out.")
    e_ = o["universe"]["ccd_over_f33_enrollment"]
    a(f"4. *Universe.* F-33 '2018' and CCD 2018 are the same year and the same enrollment "
      f"({100 * e_['exact']:.0f}% exactly equal, {100 * e_['within10pct']:.1f}% within 10%), "
      "and the four teacher fields are a strict subset of instruction salaries. "
      "One mismatch matters: CCD FTE includes pre-K teachers, and the within-CZ regression "
      "of log pay per FTE on the pre-K share gives "
      f"{o['universe']['pay_on_prek_and_aides']['coef']['prek_share']:+.2f} "
      f"(SE {o['universe']['pay_on_prek_and_aides']['se']['prek_share']:.2f}), i.e. pre-K pay "
      f"is mostly missing from the numerator. Dividing by non-pre-K FTE moves the gap to "
      f"{pre['value']:.4f}.")
    a(f"5. *Other numerators.* Adding instruction benefits: {ben['value']:.4f}; instruction "
      f"salaries per teacher: {ins['value']:.4f} ({ins['n']} CZs; {ins['value_c']:.4f} on the "
      f"common CZs); regular-program teacher salaries only (which drops special-education, "
      f"vocational and other teachers from the numerator but not the denominator): "
      f"{reg['value']:.4f}.")
    a(f"6. *Full-coverage version.* Filling zero-coded wage bills with 0.88 x instruction "
      f"salaries and collapsing NYC into one district gives {hyb['value']:.4f} (SE "
      f"{hyb['se']:.4f}) on all {hyb['n']} CZs; also dividing by non-pre-K FTE gives "
      f"{hybp['value']:.4f} (SE {hybp['se']:.4f}). Collapsing NYC alone (zeros missing): "
      f"{nyc['value']:.4f} on {nyc['n']} CZs.\n")
    a("**FTE per pupil** (`gap_teachers_pp_locale` = "
      f"{bt['value']:.4f}, bootstrap SE {bt['se']:.4f}, {bt['n']} CZs).\n")
    a("1. The screen has no effect on it: the gap is a ratio of location totals of the "
      "unscreened `teachers` and `pupils`. Applying the screen changes it by "
      f"{_find(T, 'apply the [0.5%')['value'] - bt['value']:+.4f}.")
    a(f"2. CCD FTE includes pre-K teachers and CCD enrollment includes pre-K pupils; removing "
      f"both gives {k12['value']:.4f} (SE {k12['se']:.4f}). Removing only the teachers "
      f"(mismatched) would give {_find(T, 'numerator excludes pre-K')['value']:.4f}.")
    a(f"3. The sign is carried by the largest CZs. The unweighted CZ average is "
      f"{tuw['value']:.4f} (SE {tuw['se']:.4f}), dropping the five largest CZs gives "
      f"{t5['value']:.4f}, and restricting to the CZs with a salary gap gives "
      f"{tsal['value']:.4f}. Adding aides to the numerator gives {taid['value']:.4f}.")
    a("4. Numerator and denominator are the same CCD district universe (FTE teachers over "
      "students, both with pre-K), so the gap is not affected by the F-33 issues above.\n")
    a("**Sample rules matter for the enrollment gap** (section 5): the pupils gap moves from "
      "0.376 to 0.155 when each location must hold >= 20% of CZ pupils, and to -0.015 with "
      "unweighted CZ averaging.\n")


def write_md(path: Path, o: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    m = o["meta"]
    a("# District moments: audit of pay and FTE counts (T3)\n")
    a(f"Generated {m['generated']} by `district_audit.py` (year {m['year']}, `locale` "
      f"partition, whole-CZ bootstrap B = {m['B']}, seed {m['seed']}). Every number here "
      "is recomputed with `data_estimate.py`'s own panel, partition and `czone_gaps`; the "
      "baseline rows reproduce the tracked moments exactly.\n")

    write_summary(L, o)
    a("## 1. How the salary gap is built\n")
    a("`gap_salary_real_locale` (2018) = pupil-weighted mean across commuting zones of\n\n"
      "    log( salary_real[suburb] / salary_real[city] ),\n\n"
      "where, within a location (the CZ's suburb or city districts), the statistic is the "
      "**teacher-FTE-weighted mean of district `salary_real`** (`stat = 'mean'`, "
      "`weight = 'teachers'`), and the CZ weight is the CZ's partitioned pupils. The steps:\n")
    a("1. **Sample frame:** CCD directory 2018 (Urban Institute API `school-districts/ccd/"
      "directory/2018`), regular local districts only (agency type 1, or 2 = supervisory-union "
      "component). `pupils` = CCD `enrollment` (> 0; F-33 responsible enrollment if missing); "
      "`teachers` = CCD `teachers_total_fte` (> 0), i.e. FTE teachers of all levels "
      "(pre-K, K, elementary, secondary, ungraded), not aides, coordinators or counselors.")
    a("2. **Numerator (F-33 2018, Urban API `school-districts/ccd/finance/2018`):** "
      "`wagebill_teachers` = row sum (NaN-skipping, sentinel codes -1/-2/-3/-9/-99 set to "
      "missing) of the four teacher-salary fields `salaries_teachers_regular_prog`, "
      "`salaries_teachers_sped`, `salaries_teachers_vocational`, `salaries_teachers_other_ed`. "
      "These are **teacher salaries only**: no benefits, no aides, substitutes, "
      "coordinators or other instruction salaries (those sit in `salaries_instruction`, "
      "not used; the fall-back to it exists only when all four fields are missing). "
      "`salary_per_teacher = wagebill_teachers / teachers`.")
    a("3. **Deflator (NCES EDGE CWIFT):** the district-level `LEA_CWIFTEST` of the release "
      "labelled 2018 (school year 2018-19), merged on the nearest release within 2 years "
      "(exact for 2018: coverage 96% of districts; on this machine the other LEA vintages are 2013, "
      "2015, 2019, 2021, 2022). `salary_real = salary_per_teacher / cwift`. CWIFT is a "
      "labor-cost index built from wages of comparable college-educated non-teachers, so "
      "the deflated gap is teacher pay relative to local comparable-worker pay levels, not "
      "cost-of-living-adjusted pay.")
    a("4. **Screen:** `screen_outliers` blanks `salary_real` (and `salary_per_teacher`, "
      "`teachers_pp`, `str_ratio`, ...) outside the panel-wide (all regular districts, one "
      "year) [0.5%, 99.5%] quantiles. **It acts on the derived column only**: `gap_teachers_pp` "
      "and `gap_salary_per_teacher` are ratios of location totals of the unscreened "
      "`teachers`, `pupils`, `wagebill_teachers`, so the screen never touches them "
      "(only the `mean`-type `salary_real` gap is screened).")
    a("5. **Partition:** districts with a CZ (Dorn 1990 crosswalk via CCD county), >= 500 "
      "pupils, NCES locale 11-13 (city) or 21-23 (suburb); a CZ is kept if both groups hold "
      ">= 10% of its partitioned pupils (188 CZs). A CZ has a salary gap only if both "
      "locations contain at least one district with a non-missing `salary_real`; that gives 176.")
    a("6. **Gap and average:** log of the ratio of the suburb to the city teacher-weighted "
      "mean; averaged over CZs with weights = CZ partitioned pupils (all partitioned "
      "districts, including those without a salary).\n")

    a("### Sample attrition\n")
    a("| step | districts | pupils (millions) |")
    a("|---|---|---|")
    for lab, n, p in o["attrition"]:
        a(f"| {lab} | {n:,} | {p / 1e6:.2f} |")
    cz = o["cz_salary"]
    a(f"\n{cz['n_lost']} of {cz['n_cz']} CZs have no salary gap, holding "
      f"**{100 * cz['pupil_share_lost']:.1f}% of the partitioned pupils** (the other targets "
      "use all of them).\n")
    a("| group | districts | pupils (m) | FTE valid | F-33 wage bill missing | wage bill == 0 | "
      "... share of pupils | CWIFT missing | `salary_real` > 0 |")
    a("|---|---|---|---|---|---|---|---|---|")
    for c in o["coverage_by_group"]:
        a(f"| {c['group']} | {c['districts']:,} | {c['pupils'] / 1e6:.2f} | "
          f"{100 * c['teachers_valid']:.1f}% | {100 * c['wb_missing']:.1f}% | "
          f"{100 * c['wb_zero']:.1f}% | {100 * c['wb_zero_pupils']:.1f}% | "
          f"{100 * c['cwift_missing']:.1f}% | {100 * c['salary_real_valid_pos']:.1f}% |")

    a("\nCZs without a salary gap (largest first):\n")
    a("| CZ | largest district | pupils | share of partition | share of CZ pupils in districts with `salary_real` > 0 |")
    a("|---|---|---|---|---|")
    for r in o["lost_cz"][:8]:
        a(f"| {r['czone']:.0f} | {r['name']} | {r['pupils']:,.0f} | {100 * r['weight']:.1f}% | "
          f"{100 * r['pupils_valid_salary_share']:.0f}% |")
    a("\n## 2. Do numerator and denominator cover the same universe?\n")
    zs = o["zero_coded_states"]
    a("**Findings that change what the salary gap measures:**\n")
    a(f"1. **Zero-coded F-33 teacher salaries.** {o['zero_all']['n']:,} of the regular districts "
      f"({100 * o['zero_all']['share']:.1f}%) report 0 in all four teacher-salary fields while "
      "instruction and total salaries are positive: it is missing data coded as zero, not "
      "a sentinel. It is concentrated in whole states: "
      + "; ".join(f"state {r['st']}: {int(r['zero'])}/{int(r['n'])} districts "
                  f"({100 * r['zero_share']:.0f}%)" for r in zs)
      + " (02 = Alaska, 35 = New Mexico, 17 = Illinois, 23 = Maine, 50 = Vermont, 33 = New "
      "Hampshire). Because > 5% of the panel is exactly zero, the lower screen quantile "
      f"is {o['screen_cuts']['nominal_lo']:.0f}: **the screen has no lower bite**; "
      f"the upper cuts are ${o['screen_cuts']['nominal_hi']:,.0f} nominal and "
      f"${o['screen_cuts']['real_hi']:,.0f} deflated. Zeros are kept as salaries of 0 in "
      "the location averages. They drop a CZ only when a whole group is zero (log of 0). "
      "In CZs where some districts are zero and others positive (mostly OH, CT, MA, ME/VT/NH "
      "and a few IL/AK ones) the zeros lower that location's mean (see the "
      "'zeros missing' rows in section 3).")
    n = o["nyc"]
    a(f"2. **New York City has no salary data in the panel.** F-33 reports the whole city "
      f"under the NYC Chancellor's Office (`{NYC_CHANCELLOR}`, agency type 3, dropped by the "
      f"frame filter; F-33 enrollment {n['chancellor_enrollment']:,.0f}, teacher wage bill "
      f"${(n['chancellor_wagebill'] or 0) / 1e9:.2f} billion), while the CCD frame splits it "
      f"into {n['n_geo_districts']} geographic districts ({n['pupils']:,.0f} pupils, "
      f"{n['teachers']:,.0f} FTE) whose F-33 records are all `-2` (not reported). "
      f"{n['f33_wagebill_nan']} of {n['n_geo_districts']} have a missing wage bill, so the "
      f"largest city group in the country has no salary and CZ {n['czone']:.0f} (New York) "
      "drops from the salary gap while staying in the pupils, FTE-per-pupil and income "
      "targets.")
    a(f"3. **Combined effect:** the salary gap averages over {cz['n_with_gap']} CZs holding "
      f"{100 * (1 - cz['pupil_share_lost']):.1f}% of the partitioned pupils. The excluded CZs "
      "include the second and third largest (New York, Chicago) and Albuquerque. The salary "
      "target and the other three targets are therefore averages over different populations "
      "(see 'restrict to the CZs that have a salary gap' rows for the other targets).\n")
    u = o["universe"]
    t = u["teacher_over_instruction_salary"]
    a("**Other checks.**\n")
    a(f"- *Teacher salaries are a strict subset of instruction salaries.* Among non-zero "
      f"districts the ratio (four teacher fields / `salaries_instruction`) has median "
      f"{t['50%']:.2f}, 5th-95th percentile [{t['5%']:.2f}, {t['95%']:.2f}], and exceeds 1.001 "
      f"in {u['n_teacher_salary_above_instruction']} districts. The remainder (about 12%) is "
      "instruction pay for aides, substitutes and others, which the CCD teacher FTE does not "
      "count: using `salaries_instruction` over teachers alone (section 3) overstates pay per "
      "teacher by a factor that varies with aide intensity.")
    e = u["ccd_over_f33_enrollment"]
    a(f"- *Same districts, same year.* CCD 2018 enrollment equals F-33 `enrollment_fall_"
      f"responsible` exactly in {100 * e['exact']:.1f}% of sample districts and within 10% in "
      f"{100 * e['within10pct']:.1f}% ({e['n_off_by_over_10pct']} differ by more). F-33 "
      "'2018' matches CCD 2018 (exact match with the 2017 and 2019 CCD in only "
      + ", ".join(f"{100 * v['exact']:.1f}% ({y})" for y, v in u['f33_enrollment_vs_ccd_year'].items() if y != 2018)
      + "), so there is no one-year offset between numerator and denominator.")
    rg = u["pay_on_prek_and_aides"]
    if rg:
        a(f"- *Pre-K teachers and aides.* CCD FTE includes pre-K teachers (mean {100 * u['ccd_prek_teachers_share_mean']:.1f}% "
          "of teacher FTE), whose pay may sit outside the four F-33 fields. Within-CZ "
          f"regression of log salary per FTE on the pre-K FTE share and the aide/teacher ratio "
          f"(non-zero districts, pupil-weighted, CZ-clustered, n = {rg['n']}): pre-K share "
          f"{rg['coef']['prek_share']:+.3f} (SE {rg['se']['prek_share']:.3f}); aide ratio "
          f"{rg['coef']['aide_ratio']:+.3f} (SE {rg['se']['aide_ratio']:.3f}). A coefficient near "
          "-1 on the pre-K share would indicate that pre-K pay is missing from the numerator.")
    a("- *Timing:* F-33 wage bills are fiscal-year flows; CCD FTE are fall counts; both "
      "refer to the school year labelled 2018 in the Urban API.")
    a("- *Teacher quality composition:* pay per FTE mixes the wage schedule and teacher "
      "experience/education; the target is `Delta log[kappa_l E(h_T^gamma | T, l)]` per "
      "the calibration note, so this is intended.\n")

    b = o["baseline"]
    a("## 3. Sensitivity of the salary gap\n")
    a(f"Baseline {b['salary_real']:.5f} on {o['cz_salary']['n_with_gap']} CZs. "
      "Whole-CZ bootstrap SEs, same seed and B as `district_bootstrap.json`; a variant that "
      "blanks a district keeps the CZ weights and partition of the baseline unless the "
      "variant changes the sample rule.\n")
    a("`gap on baseline CZs` restricts the same variant to the CZs that have a baseline "
      "gap (176 for salary), so that a change in coverage is separated from a change in "
      "the data; the bootstrap resamples the CZs of each column.\n")
    a(f"The hybrid rows use {o['hybrid_share']:.3f} = median of (teacher salaries / "
      "instruction salaries) among non-zero sample districts to fill zero-coded wage bills.\n")
    L.extend(var_table(o["salary"], o["salary"][0]["value"]))

    a("\n## 4. Sensitivity of the FTE-per-pupil gap\n")
    a(f"Baseline {b['teachers_pp']:.5f} on {o['fte'][0]['n']} CZs (equals minus the "
      "student-teacher-ratio gap).\n")
    L.extend(var_table(o["fte"], o["fte"][0]["value"]))

    a("\n## 5. Other targets under the same sample rules\n")
    for key in ("pupils", "seda_lninc"):
        rows = [r for r in o["other"] if r["group"] == key]
        a(f"\n**`{key}`**\n")
        L.extend(var_table(rows, rows[0]["value"]))

    a("\n## 6. Influence of the largest CZs\n")
    a("| CZ | largest district in the CZ | weight | gap pupils | gap salary | gap FTE/pupil | gap income |")
    a("|---|---|---|---|---|---|---|")
    for r in o["top_cz"]:
        gp = r["gaps"]
        a(f"| {r['czone']:.0f} | {r['name']} | {100 * r['weight']:.1f}% | {_f(gp['pupils'], 3)} | "
          f"{_f(gp['salary_real'], 3)} | {_f(gp['teachers_pp'], 3)} | {_f(gp['seda_lninc'], 3)} |")
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
