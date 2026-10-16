#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
crdc_teacher_allocation.py -- SEALED, exploratory: how teachers of different
experience and certification are allocated across the two locations of a commuting
zone (CZ), from the school-level CRDC counts of the Urban Institute Education Data
API (`crdc_teachers_staff.py`, school year 2017-18).

PRE-REGISTERED OUTCOMES (fixed before any gap was computed).  Locale partition
(suburb NCES locale 21-23 minus city 11-13; `data_estimate.partition`, 500-pupil and
10%-group-share cuts), the within-CZ gap in a location ratio, location 2 minus
location 1, averaged over CZs with the CZ weights of `district_bootstrap.py`
(partitioned pupils), whole-CZ bootstrap (`district_bootstrap.draw_counts`, B = 1999,
seed 20260929):

  PRIMARY (i)   new-teacher share    = sum(first-year + second-year teacher FTE) /
                                       sum(teacher FTE)
  PRIMARY (ii)  certified share      = sum(certified teacher FTE) / sum(teacher FTE)

each a ratio of totals over the schools of the districts of a location (the
convention of `data_estimate.Gap`, stat = "ratio").  A school enters (i) if teacher
FTE > 0 and the first- and second-year counts are non-negative, and (ii) if teacher
FTE > 0 and the certified count is non-negative; negative values are the API's codes
for missing / not applicable / suppressed.  Schools without a district id (`leaid`)
cannot be linked and are dropped.

Data note (school-level file, no location split): `teachers_fte_crdc` is integer
valued while the certified, uncertified, first- and second-year counts carry decimals,
and certified + uncertified exceeds the total by less than one FTE in about half of
the schools; the total looks truncated to an integer.  The primary definitions follow
the pre-registration (denominator = teacher FTE); the secondary variants, fixed after
this school-level check and before any gap value was viewed, are `new_alt` / `cert_alt` (denominator = certified plus
uncertified FTE, in the same decimal units as the numerators) and `new_trim` /
`cert_trim` (drop schools whose numerator exceeds teacher FTE by more than one FTE,
i.e. by more than the integer truncation).

SEALED: the gap values are written ONLY to
    estimates/settings/sealed/crdc_teacher_allocation.json
Nothing in this script prints, logs or writes a gap, a level or a standard error
anywhere else; the log reports coverage (schools, districts, commuting zones) only.

Usage
    python3 crdc_teacher_allocation.py                  # reads raw/crdc_teachers_staff_2017.parquet
    python3 crdc_teacher_allocation.py --selftest       # synthetic planted gaps only

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import data_estimate as de  # noqa: E402
import district_bootstrap as db  # noqa: E402

LOG = logging.getLogger("crdc_teacher_allocation")

SEED = db.SEED
B_DEFAULT = db.B_DEFAULT
YEAR_DEFAULT = 2017
SCHEME = "locale"

# outcome key -> (label, primary?)
OUTCOMES: dict[str, tuple[str, bool]] = {
    "new": ("share of first- and second-year teachers in teacher FTE", True),
    "cert": ("share of certified teachers in teacher FTE", True),
    "new_trim": ("new-teacher share, schools with first+second-year FTE <= teacher FTE + 1", False),
    "cert_trim": ("certified share, schools with certified FTE <= teacher FTE + 1", False),
    "new_alt": ("(first- + second-year FTE) / (certified + uncertified FTE)", False),
    "cert_alt": ("certified FTE / (certified + uncertified FTE)", False),
}


# --------------------------------------------------------------------------
# Schools -> districts
# --------------------------------------------------------------------------

def school_outcomes(raw: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Per outcome: school-level numerator and denominator, NaN where the school does not
    enter (negative API codes are missing)."""
    t, c, u = raw.teachers_fte_crdc, raw.teachers_certified_fte, raw.teachers_uncertified_fte
    f1, f2 = raw.teachers_first_year_fte, raw.teachers_second_year_fte
    linked = raw.leaid.astype(str).str.fullmatch(r"\d{7}")
    new_ok = linked & (t > 0) & (f1 >= 0) & (f2 >= 0)
    cert_ok = linked & (t > 0) & (c >= 0)
    alt_ok = linked & (c >= 0) & (u >= 0) & ((c + u) > 0)
    newalt_ok = alt_ok & (f1 >= 0) & (f2 >= 0)
    spec = {
        "new": (new_ok, f1 + f2, t),
        "cert": (cert_ok, c, t),
        "new_trim": (new_ok & ((f1 + f2) <= t + 1), f1 + f2, t),
        "cert_trim": (cert_ok & (c <= t + 1), c, t),
        "new_alt": (newalt_ok, f1 + f2, c + u),
        "cert_alt": (alt_ok, c, c + u),
    }
    out = {}
    for k, (ok, num, den) in spec.items():
        out[k] = pd.DataFrame({"leaid": raw.leaid.astype(str), "ok": ok.to_numpy(),
                               "num": num.where(ok), "den": den.where(ok)})
    return out


def district_sums(schools: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per leaid: sums of the numerator and denominator of every outcome."""
    parts = []
    for k, df in schools.items():
        d = df[df.ok].groupby("leaid").agg(**{f"crdc_{k}_num": ("num", "sum"),
                                              f"crdc_{k}_den": ("den", "sum"),
                                              f"crdc_{k}_n": ("num", "size")})
        parts.append(d)
    return pd.concat(parts, axis=1).reset_index().rename(columns={"index": "leaid"})


def gap_defs() -> list["de.Gap"]:
    G = []
    for k, (label, _) in OUTCOMES.items():
        G.append(de.Gap(f"crdc_{k}", label, "teacher allocation (sealed)", "ratio",
                        num=f"crdc_{k}_num", den=f"crdc_{k}_den", contrast="diff"))
    return G


# --------------------------------------------------------------------------
# Gaps and bootstrap (values leave this function only through the returned dict)
# --------------------------------------------------------------------------

def compute(panel: pd.DataFrame, B: int, seed: int) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    part, g = db.cz_gap_table(panel, gap_defs(), SCHEME)
    w = g.pupils.to_numpy(float)
    counts = db.draw_counts(len(g), B, seed)
    res: dict[str, Any] = {}
    for k, (label, primary) in OUTCOMES.items():
        col = f"gap_crdc_{k}"
        if col not in g or g[col].notna().sum() < 5:
            res[k] = {"label": label, "primary": primary, "available": False}
            continue
        st = db.boot_wmean(g[col].to_numpy(float), w, counts)
        s = db.summarize(st)
        m = g[col].notna().to_numpy() & (w > 0)
        res[k] = {"label": label, "primary": primary, "available": True,
                  "gap_suburb_minus_city": s["value"], "boot_se": s["boot_se"],
                  "ci95_lo": s["ci95_lo"], "ci95_hi": s["ci95_hi"],
                  "n_cz": s["n"], "mean_n_cz_per_draw": s["mean_n_per_draw"],
                  "n_eff_cz": s["n_eff_cz"],
                  "level_city": float(np.average(g[f"crdc_{k}_1"].to_numpy(float)[m], weights=w[m])),
                  "level_suburb": float(np.average(g[f"crdc_{k}_2"].to_numpy(float)[m], weights=w[m]))}
    return res, part, g


def coverage(raw: pd.DataFrame, schools: dict[str, pd.DataFrame], panel: pd.DataFrame,
             part: pd.DataFrame, g: pd.DataFrame) -> dict[str, Any]:
    """Counts only: how much of the CRDC file reaches the locale partition."""
    d = part[part.year == db.YEAR]
    out: dict[str, Any] = {
        "schools_in_file": int(len(raw)),
        "schools_with_leaid": int(raw.leaid.astype(str).str.fullmatch(r"\d{7}").sum()),
        "districts_in_file": int(raw.leaid[raw.leaid.astype(str).str.fullmatch(r"\d{7}")].nunique()),
        "panel_districts": int(panel.leaid.nunique()),
        "partition_districts": int(len(d)), "partition_czs": int(d.czone.nunique()),
        "partition_districts_by_location": {("city" if k == 1 else "suburb"): int(v)
                                            for k, v in d.groupby("loc_id").size().items()},
        "by_outcome": {}}
    for k in OUTCOMES:
        den = f"crdc_{k}_den"
        has = d[den].notna() & (d[den] > 0) if den in d else pd.Series(False, index=d.index)
        both = d[has].groupby("czone").loc_id.nunique().eq(2)
        out["by_outcome"][k] = {
            "schools_entering": int(schools[k].ok.sum()),
            "districts_entering_file": int(schools[k][schools[k].ok].leaid.nunique()),
            "partition_schools_entering": int(d.loc[has, f"crdc_{k}_n"].sum()) if has.any() else 0,
            "partition_schools_entering_by_location": {
                ("city" if kk == 1 else "suburb"): int(v)
                for kk, v in d[has].groupby("loc_id")[f"crdc_{k}_n"].sum().items()},
            "partition_districts_with_data": int(has.sum()),
            "share_of_partition_districts_with_data": float(has.mean()),
            "partition_districts_with_data_by_location": {
                ("city" if kk == 1 else "suburb"): int(v) for kk, v in d[has].groupby("loc_id").size().items()},
            "czs_with_both_locations": int(both.sum()),
            "czs_in_bootstrap_table": int(g[f"gap_crdc_{k}"].notna().sum()) if f"gap_crdc_{k}" in g else 0}
    return out


def finite_check(res: dict[str, Any]) -> bool:
    ok = True
    for r in res.values():
        if r.get("available"):
            ok &= all(np.isfinite(r[x]) for x in ("gap_suburb_minus_city", "boot_se", "ci95_lo",
                                                 "ci95_hi", "level_city", "level_suburb"))
    return bool(ok)


# --------------------------------------------------------------------------
# Self-test on synthetic data (planted gaps; nothing here comes from the CRDC file)
# --------------------------------------------------------------------------

def selftest(indir: Path, B: int = 199) -> int:
    panel = db.build_panel(indir)
    rng = np.random.default_rng(7)
    sub = panel.locale_grp.eq("suburb").to_numpy()
    den = panel.teachers.where(panel.teachers > 0)
    plant = {"new": 0.04, "cert": -0.02, "new_trim": 0.04, "cert_trim": -0.02, "new_alt": 0.04,
             "cert_alt": -0.02}
    base = {"new": 0.10, "cert": 0.95, "new_trim": 0.10, "cert_trim": 0.95, "new_alt": 0.10,
            "cert_alt": 0.95}
    for noise, tag in ((0.0, "noise-free"), (0.05, "noisy")):
        p = panel.copy()
        for k in OUTCOMES:
            sh = base[k] + plant[k] * sub + noise * rng.standard_normal(len(p))
            p[f"crdc_{k}_den"] = den
            p[f"crdc_{k}_num"] = den * sh
        res, _, g = compute(p, B, SEED)
        print(f"[selftest {tag}]")
        for k, r in res.items():
            print(f"  {k:10s} planted {plant[k]:+.4f}  recovered {r['gap_suburb_minus_city']:+.4f} "
                  f"(boot SE {r['boot_se']:.4f}, n CZ {r['n_cz']})")
        # bootstrap SE against the linearized cluster-robust SE of the same weighted mean
        col = "gap_crdc_new"
        st = db.boot_wmean(g[col].to_numpy(float), g.pupils.to_numpy(float),
                           db.draw_counts(len(g), B, SEED))
        print(f"  boot SD of draws {np.std(st['draws'], ddof=1):.5f}  linearized SE {st['se_lin']:.5f}")
    return 0


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(HERE))
    ap.add_argument("--year", type=int, default=YEAR_DEFAULT, help="CRDC year (2017 = 2017-18)")
    ap.add_argument("--B", type=int, default=B_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--outdir", default=None,
                    help="default: <indir>/estimates/settings/sealed")
    ap.add_argument("--selftest", action="store_true",
                    help="run on synthetic school data with planted gaps (prints synthetic numbers only)")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    indir = Path(args.indir).resolve()
    if args.selftest:
        return selftest(indir)
    outdir = (Path(args.outdir).resolve() if args.outdir
              else indir / "estimates" / "settings" / "sealed")
    outdir.mkdir(parents=True, exist_ok=True)

    path = indir / "raw" / f"crdc_teachers_staff_{args.year}.parquet"
    if not path.exists():
        LOG.error("no %s: run crdc_teachers_staff.py first", path.name)
        return 1
    raw = pd.read_parquet(path)
    schools = school_outcomes(raw)
    dsum = district_sums(schools)
    panel = db.build_panel(indir)
    panel = panel.merge(dsum, on="leaid", how="left")
    res, part, g = compute(panel, args.B, args.seed)
    cov = coverage(raw, schools, panel, part, g)
    LOG.info("CRDC %d: %d schools in the file, %d with a district id (%d districts)", args.year,
             cov["schools_in_file"], cov["schools_with_leaid"], cov["districts_in_file"])
    LOG.info("locale partition: %d districts in %d CZs (%s)", cov["partition_districts"],
             cov["partition_czs"], cov["partition_districts_by_location"])
    for k, c in cov["by_outcome"].items():
        LOG.info("  outcome %-9s schools entering %6d (%6d in partition districts) | partition "
                 "districts with data %5d (%.1f%%; %s) | CZs with both locations %3d", k,
                 c["schools_entering"], c["partition_schools_entering"],
                 c["partition_districts_with_data"], 100 * c["share_of_partition_districts_with_data"],
                 c["partition_districts_with_data_by_location"], c["czs_with_both_locations"])
    meta = {
        "sealed": True,
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "data/spatial/crdc_teacher_allocation.py",
        "data": f"Urban Institute Education Data API schools/crdc/teachers-staff/{args.year}/ "
                f"(CRDC, school year {args.year}-{str(args.year + 1)[-2:]}), raw/crdc_teachers_staff_{args.year}.parquet",
        "year": args.year, "partition": SCHEME,
        "geo_unit": panel.attrs.get("geo_unit"),
        "gap": "location 2 (suburb, NCES locale 21-23) minus location 1 (city, 11-13) of the location "
               "ratio of totals, enrollment-weighted mean across CZs (weights = partitioned pupils)",
        "preregistered_primary": ["new", "cert"],
        "definitions": {
            "new": "sum(first-year + second-year teacher FTE) / sum(teacher FTE); schools with "
                   "teacher FTE > 0 and non-negative first- and second-year counts",
            "cert": "sum(certified teacher FTE) / sum(teacher FTE); schools with teacher FTE > 0 and a "
                    "non-negative certified count",
            "new_trim": "new, dropping schools with first+second-year FTE > teacher FTE + 1",
            "cert_trim": "cert, dropping schools with certified FTE > teacher FTE + 1",
            "new_alt": "sum(first-year + second-year FTE) / sum(certified + uncertified FTE); schools "
                       "with non-negative counts and a positive certified + uncertified sum",
            "cert_alt": "sum(certified FTE) / sum(certified + uncertified FTE); schools with both "
                        "non-negative and a positive sum"},
        "data_note": "teachers_fte_crdc is integer valued (truncated) while certified/uncertified/"
                     "first-/second-year FTE carry decimals; certified + uncertified exceeds the "
                     "total by < 1 FTE in about half of the schools",
        "missing_codes": "negative values (-1, -2, -3, ...) are masked; schools without leaid dropped",
        "B": args.B, "seed": args.seed,
        "rng": "numpy default_rng (PCG64) via district_bootstrap.draw_counts; one (B x n_CZ) "
               "multiplicity matrix; resampling unit = whole CZ",
        "min_pupils": 500, "min_group_share": 0.10,
        "coverage": cov,
    }
    payload = {"meta": meta, "outcomes": res}

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        return o

    out_path = outdir / "crdc_teacher_allocation.json"
    out_path.write_text(json.dumps(clean(payload), indent=1))
    LOG.info("wrote %s (%d outcomes, all finite: %s); values are sealed and not logged",
             out_path.name, sum(bool(r.get("available")) for r in res.values()), finite_check(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
