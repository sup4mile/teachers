#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
housing_costs.py -- suburb-minus-city gaps in housing costs from the ACS 5-year
2014-18 school-district tables (`acs_housing.py`), as discipline for a recurring
residence cost in the spatial model.

The model's residence cost is p_2 - p_1 = r_p * Cbar, with Cbar mean consumption per
adult: a fixed (non-homothetic) goods cost of living in location 2.  This script
measures the monthly housing-cost gap between the two locations of a commuting zone
(CZ) and converts it to model units, 12 x (dollar gap per month) / (mean annual
expenditures per consumer unit, BLS Consumer Expenditure Survey 2018).

Panel, partitions, CZ weights and bootstrap are exactly those of
`district_bootstrap.py`: the district panel of `data_estimate.assemble`, the ACS
tables merged on `leaid` as the earnings tables are (`load_income_columns`: the
"Remainder of <state>" rows are dropped, annotation codes are masked), the
`partition` of `data_estimate.py` (locale = suburb minus city; revenue and ses =
high minus low half of the CZ), `czone_gaps` (location statistic, then
location-2-minus-location-1 gap within the CZ) and the enrollment-weighted mean across
CZs; whole-CZ bootstrap with `district_bootstrap.draw_counts` (B = 1999, seed
20260929, one multiplicity matrix per partition shared by all statistics).  Same
seed and CZ order as `district_bootstrap.py`, so the locale draws are the same ones
(the report's reproduction check compares the published `pm_ftyr_all` with the
one rebuilt here).

Statistics, for each partition
  (a) pooled-location median monthly housing cost: B25104 bin counts summed over
      the districts of each location, median by linear interpolation inside the
      bin, "No cash rent" units excluded (they have no cost).  Gap in logs and in
      dollars.  (The report gives the mean log error of the interpolated district
      medians against the published B25105.)
  (b) occupied-unit-weighted mean of district medians (weights B25003_001) of
      B25105 (monthly housing costs), B25064 (gross rent), B25088 (selected
      monthly owner costs: all, with, without a mortgage) and B25077 (home value):
      the gap in the average log median (`ln_`), the log ratio of average levels
      (`lvl_`), and the gap in average dollars (`usd_`).  Tenure-specific weights
      (renters for rent, owners for owner costs and value) are a `_tw` variant.
  (c) owner share: owner-occupied over occupied units (B25003), pooled over the
      districts of a location, gap in levels.
  Sensitivity: pooled-location MEAN monthly cost from the bin midpoints (top bin
  valued at $3,500 and $4,000).
  Scale: the dollar gaps are also expressed as 12 x gap / CE 2018 mean annual
  expenditures per consumer unit ($61,224.13; BLS Table R-1) = r_p, and as a share
  of the CZ-pooled median full-time year-round (FTYR) earnings built from the
  B20005 bins of `acs_earnings.py` (both locations pooled; ratio within the CZ,
  averaged with the CZ weights, and ratio of the two averages).

Writes (under estimates/settings/):
    housing_costs.json, housing_costs.md, housing_costs_draws.csv (locale draws).

Usage
    python3 housing_costs.py                  # B = 1999, seed 20260929
    python3 housing_costs.py --B 199          # quick check

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

LOG = logging.getLogger("housing_costs")

YEAR = db.YEAR
SEED = db.SEED
B_DEFAULT = db.B_DEFAULT
SCHEMES = ("locale", "revenue", "ses")
SCHEME_DESC = {"locale": "locale (suburb 21-23 minus city 11-13)",
               "revenue": "revenue (above minus below the CZ median local revenue per pupil)",
               "ses": "ses (above minus below the CZ median SEDA SES)"}

# BLS Consumer Expenditure Survey 2018, Table R-1 "All consumer units: annual
# detailed expenditure means" (calendar year 2018, dollars of 2018).
CE_TOTAL_EXP = 61224.13          # average annual expenditures per consumer unit
CE_HOUSING = 20090.59            # of which housing
CE_SHELTER = 11747.43            # of which shelter (owned + rented dwellings + other lodging)
CE_PENSIONS = 7295.55            # personal insurance and pensions (saving-like items in the total)
CE_UNITS_K = 131439              # consumer units, thousands
CE_PERSONS, CE_EARNERS, CE_INCOME = 2.5, 1.3, 78635
CE_URL = "https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf"
CE_REPORT_URL = "https://www.bls.gov/opub/reports/consumer-expenditures/2018/"

# ACS district medians: key -> (column, label, tenure of the natural weight, bottom code)
SERIES: dict[str, tuple[str, str, str, float]] = {
    "h_all": ("B25105_001E", "median monthly housing costs, occupied units (B25105)", "occ", 99),
    "h_rent": ("B25064_001E", "median gross rent (B25064)", "rent", 99),
    "h_own": ("B25088_001E", "median selected monthly owner costs, all owner units (B25088_001)",
              "own", 99),
    "h_own_mort": ("B25088_002E", "median selected monthly owner costs, units with a mortgage "
                   "(B25088_002)", "own", 99),
    "h_own_nomort": ("B25088_003E", "median selected monthly owner costs, units without a "
                     "mortgage (B25088_003)", "own", 99),
    "h_value": ("B25077_001E", "median value of owner-occupied units (B25077)", "own", 9999),
}
TENURE_WEIGHT = {"occ": "acsh_occ", "own": "acsh_own", "rent": "acsh_rent"}
NO_TW = ("h_all", "h_own_mort", "h_own_nomort")     # no tenure-weighted variant

# B25104 monthly housing cost bins: 15 cost intervals (_002 .. _016) and "No cash rent"
H_EDGES = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500,
                    3000, np.inf])
H_MIDS = np.array([50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1250, 1750, 2250, 2750])
TOP_BIN_VALUES = (3500, 4000)


# --------------------------------------------------------------------------
# Merge of the ACS housing tables onto the panel (as load_income_columns does)
# --------------------------------------------------------------------------

def load_housing_columns(panel: pd.DataFrame, raw_dir: Path) -> pd.DataFrame:
    """Merge the ACS housing columns onto the district panel on `leaid`.

    As for the earnings tables: 'Remainder of <state>' rows (code 99999) are not
    districts and are dropped, ids are de-duplicated, and Census annotation codes
    (negative) are masked to NaN.  Bottom-coded medians (<= $99 a month, <= $9,999
    of value) are dropped as the earnings code drops the bottom-coded 2,499; top
    codes (4,001 / 3,501 / 1,501 / 2,000,001) are kept as lower bounds.
    """
    path = raw_dir / f"acs_housing_school_districts_{YEAR}.parquet"
    if not path.exists():
        raise SystemExit(f"no {path.name}: run acs_housing.py first")
    a = pd.read_parquet(path)
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid")
    out = pd.DataFrame({"leaid": a.leaid})
    for k, (col, _, _, bottom) in SERIES.items():
        x = a[col].where(a[col] > bottom)               # drops negative annotation codes
        out[k] = x
        out[f"ln_{k}"] = np.log(x)
    out["acsh_occ"] = a["B25003_001E"].where(a["B25003_001E"] > 0)
    out["acsh_own"] = a["B25003_002E"].where(a["B25003_002E"] >= 0)    # zero owners is a count
    out["acsh_rent"] = a["B25003_003E"].where(a["B25003_003E"] >= 0)
    bpath = raw_dir / f"acs_housing_bins_school_districts_{YEAR}.parquet"
    if not bpath.exists():
        raise SystemExit(f"no {bpath.name}: run acs_housing.py first")
    bn = pd.read_parquet(bpath)
    bn = bn[~bn.leaid.str.endswith("99999")].drop_duplicates("leaid")
    bi = bn.set_index("leaid")
    newcols = {f"bin_h_{j:02d}": bi[f"B25104_{j + 2:03d}E"].where(
        bi[f"B25104_{j + 2:03d}E"] >= 0).reindex(out.leaid).to_numpy() for j in range(15)}
    newcols["bin_h_nocash"] = bi["B25104_017E"].where(
        bi["B25104_017E"] >= 0).reindex(out.leaid).to_numpy()
    out = pd.concat([out, pd.DataFrame(newcols, index=out.index)], axis=1)
    keep = [c for c in out.columns if c not in panel.columns or c == "leaid"]
    attrs = dict(panel.attrs)
    merged = panel.merge(out[keep], on="leaid", how="left")
    merged.attrs.update(attrs)
    return merged


# --------------------------------------------------------------------------
# Gap definitions (data_estimate.Gap, computed by czone_gaps) and pooled medians
# --------------------------------------------------------------------------

def gap_defs() -> list["de.Gap"]:
    G: list[de.Gap] = []
    for k, (col, label, tenure, _) in SERIES.items():
        variants = [("", TENURE_WEIGHT["occ"], "occupied-unit weights")]
        if k not in NO_TW and tenure != "occ":
            variants.append(("_tw", TENURE_WEIGHT[tenure],
                             "renter weights" if tenure == "rent" else "owner weights"))
        for tag, wcol, wlab in variants:
            G.append(de.Gap(f"ln_{k}{tag}", f"log {label}, average of district log medians "
                            f"({wlab})", "residence cost (log gap)", "mean", col=f"ln_{k}",
                            weight=wcol, contrast="diff"))
            G.append(de.Gap(f"lvl_{k}{tag}", f"log ratio of average {label} ({wlab})",
                            "residence cost (log gap, level aggregation)", "mean", col=k,
                            weight=wcol, contrast="log"))
            G.append(de.Gap(f"usd_{k}{tag}", f"dollar gap in the average {label} ({wlab})",
                            "residence cost (dollars)", "mean", col=k, weight=wcol,
                            contrast="diff"))
    G.append(de.Gap("own_share", "owner-occupied share of occupied units (B25003), pooled over "
                    "the districts of a location", "tenure", "ratio", num="acsh_own",
                    den="acsh_occ", contrast="diff"))
    return G


def grouped_median(counts: np.ndarray, edges: np.ndarray) -> float:
    """Median of a distribution known through counts in `edges` intervals, by linear
    interpolation inside the bin that holds the 50th percentile (NaN if that is the
    open-ended top bin or there are no observations)."""
    c = np.asarray(counts, float)
    tot = c.sum()
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    cum = np.cumsum(c)
    j = int(np.searchsorted(cum, tot / 2.0))
    if j >= len(c) - 1 and c[-1] > 0 and cum[-2] < tot / 2.0:
        return np.nan
    below = cum[j - 1] if j > 0 else 0.0
    lo, hi = edges[j], edges[j + 1]
    return float(lo + (hi - lo) * (tot / 2.0 - below) / c[j])


def grouped_mean(counts: np.ndarray, top: float) -> float:
    """Mean from the 15 monthly-cost bins: bin midpoints, the open top bin at `top`."""
    c = np.asarray(counts, float)
    tot = c.sum()
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    mids = np.append(H_MIDS, top)
    return float((c * mids).sum() / tot)


def _by_location(d: pd.DataFrame, cols: list[str], fn) -> pd.DataFrame:
    """fn(row of summed bin counts) for every (czone, loc_id): a CZ x {1, 2} frame."""
    sums = d.groupby(["czone", "loc_id"])[cols].sum(min_count=1)
    vals = pd.Series({idx: fn(row.to_numpy()) for idx, row in sums.iterrows()})
    vals.index = pd.MultiIndex.from_tuples(vals.index, names=["czone", "loc_id"])
    return vals.unstack("loc_id")


def add_pooled(part: pd.DataFrame, g: pd.DataFrame) -> pd.DataFrame:
    """Add to the CZ table: pooled-location median / mean monthly housing cost (B25104
    bins; levels `pm_h_{1,2}`, `pmean_h<top>_{1,2}` and the gaps) and FTYR earnings
    medians (B20005 bins) by location and for the two locations pooled."""
    d = part[part.year == YEAR]
    hcols = [f"bin_h_{j:02d}" for j in range(15)]
    if all(c in d.columns for c in hcols):
        wide = _by_location(d, hcols, lambda r: grouped_median(r, H_EDGES))
        for loc in (1, 2):
            g[f"pm_h_{loc}"] = g.czone.map(wide[loc]).to_numpy()
        g["gap_pm_h_log"] = np.log(g.pm_h_2 / g.pm_h_1)
        g["gap_pm_h_usd"] = g.pm_h_2 - g.pm_h_1
        for top in TOP_BIN_VALUES:
            wm = _by_location(d, hcols, lambda r, t=top: grouped_mean(r, t))
            for loc in (1, 2):
                g[f"pmean_h{top}_{loc}"] = g.czone.map(wm[loc]).to_numpy()
            g[f"gap_pmean_h{top}_log"] = np.log(g[f"pmean_h{top}_2"] / g[f"pmean_h{top}_1"])
            g[f"gap_pmean_h{top}_usd"] = g[f"pmean_h{top}_2"] - g[f"pmean_h{top}_1"]
    ecols = [f"bin_ftyr_all_{j:02d}" for j in range(20)]
    if all(c in d.columns for c in ecols):
        we = _by_location(d, ecols, lambda r: db.grouped_median(r))
        for loc in (1, 2):
            g[f"ftyr_{loc}"] = g.czone.map(we[loc]).to_numpy()
        sums = d.groupby("czone")[ecols].sum(min_count=1)
        pooled = pd.Series({cz: db.grouped_median(row.to_numpy()) for cz, row in sums.iterrows()})
        g["ftyr_12"] = g.czone.map(pooled).to_numpy()
        g["gap_pm_ftyr_all"] = np.log(g.ftyr_2 / g.ftyr_1)       # = district_bootstrap pm_ftyr_all
        if "gap_pm_h_usd" in g:
            g["share_ftyr12_pm_h"] = 12.0 * g.gap_pm_h_usd / g.ftyr_12
    return g


# --------------------------------------------------------------------------
# Bootstrap statistics
# --------------------------------------------------------------------------

def boot_ratio_of_means(xn: np.ndarray, xd: np.ndarray, w: np.ndarray, counts: np.ndarray,
                        scale: float = 1.0) -> dict[str, Any]:
    """scale * (weighted mean of xn) / (weighted mean of xd) over the CZs that have both,
    with its whole-CZ bootstrap draws (the weights cancel in the draw ratio)."""
    m = np.isfinite(xn) & np.isfinite(xd) & np.isfinite(w) & (w > 0)
    wn, wd = np.where(m, w * xn, 0.0), np.where(m, w * xd, 0.0)
    est = scale * wn.sum() / wd.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        draws = scale * (counts @ wn) / (counts @ wd)
    d = draws[np.isfinite(draws)]
    lo, hi = np.percentile(d, [2.5, 97.5])
    return dict(value=float(est), n=int(m.sum()), boot_se=float(np.std(d, ddof=1)),
                ci95_lo=float(lo), ci95_hi=float(hi), n_draws_finite=int(len(d)),
                draws=draws)


def stat_from_column(g: pd.DataFrame, col: str, w: np.ndarray, counts: np.ndarray,
                     lev1: str | None, lev2: str | None) -> dict[str, Any]:
    """Pupil-weighted mean over CZs of a CZ-level column, with bootstrap SE and interval;
    `level_1` / `level_2` are the same-sample weighted means of the location levels."""
    x = g[col].to_numpy(float)
    st = db.boot_wmean(x, w, counts)
    s = db.summarize(st)
    keep = ("value", "n", "boot_se", "ci95_lo", "ci95_hi", "n_draws_finite",
            "mean_n_per_draw", "placeholder_se", "n_eff_cz", "top5_weight_share")
    out = {k: s[k] for k in keep}
    out["draws"] = st["draws"]
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    for tag, lev in (("level_1", lev1), ("level_2", lev2)):
        if lev is not None and lev in g:
            y = g[lev].to_numpy(float)
            mm = m & np.isfinite(y)
            out[tag] = float(np.average(y[mm], weights=w[mm])) if mm.any() else None
    return out


def coverage(part: pd.DataFrame) -> dict[str, Any]:
    d = part[part.year == YEAR]
    pu = d.pupils.fillna(0)
    has_h = d["h_all"].notna() & d["bin_h_00"].notna()
    has_e = d["bin_ftyr_all_00"].notna() if "bin_ftyr_all_00" in d else has_h
    return {"n_districts": int(len(d)),
            "share_districts_with_housing_tables": float(has_h.mean()),
            "share_pupils_with_housing_tables": float(pu[has_h].sum() / pu.sum()),
            "share_districts_with_ftyr_bins": float(has_e.mean()),
            "n_cz": int(d.czone.nunique()),
            "n_districts_by_location": {int(k): int(v) for k, v in
                                        d.groupby("loc_id").size().items()},
            "top_coded_districts": {k: int((d[k] >= {"h_value": 2000001, "h_rent": 3501,
                                                     "h_own_nomort": 1501}.get(k, 4001)).sum())
                                    for k in SERIES}}


def stat_specs() -> list[dict[str, Any]]:
    """The statistics reported for every partition: key, CZ-table column, kind, label,
    level columns."""
    specs: list[dict[str, Any]] = []

    def add(key, col, kind, label, l1=None, l2=None, family=""):
        specs.append(dict(key=key, col=col, kind=kind, label=label, l1=l1, l2=l2,
                          family=family))

    add("pm_h_log", "gap_pm_h_log", "log", "(a) log pooled-location median monthly housing "
        "cost (B25104 bins), suburb - city", family="a")
    add("pm_h_usd", "gap_pm_h_usd", "usd", "(a) pooled-location median monthly housing cost "
        "(B25104 bins), dollar gap", "pm_h_1", "pm_h_2", family="a")
    for top in TOP_BIN_VALUES:
        add(f"pmean_h{top}_log", f"gap_pmean_h{top}_log", "log", "sensitivity: log pooled-location "
            f"MEAN monthly housing cost, bin midpoints, top bin at ${top:,}", family="a_mean")
        add(f"pmean_h{top}_usd", f"gap_pmean_h{top}_usd", "usd", "sensitivity: pooled-location "
            f"MEAN monthly housing cost (bin midpoints, top bin at ${top:,}), dollar gap",
            f"pmean_h{top}_1", f"pmean_h{top}_2", family="a_mean")
    for gp in gap_defs():
        if gp.key == "own_share":
            continue
        kind = ("log" if gp.key.startswith(("ln_", "lvl_")) else
                "usd_value" if gp.key.startswith("usd_h_value") else "usd")
        fam = "b_tw" if gp.key.endswith("_tw") else "b"
        add(gp.key, f"gap_{gp.key}", kind, gp.label,
            f"{gp.key}_1" if kind != "log" or gp.key.startswith("lvl_") else None,
            f"{gp.key}_2" if kind != "log" or gp.key.startswith("lvl_") else None, family=fam)
    add("own_share", "gap_own_share", "share", "(c) owner-occupied share of occupied units, "
        "suburb - city", "own_share_1", "own_share_2", family="c")
    add("pm_ftyr_all", "gap_pm_ftyr_all", "log", "reference: log pooled-location median FTYR "
        "earnings (B20005 bins) -- the `pm_ftyr_all` of district_bootstrap.py", "ftyr_1",
        "ftyr_2", family="ref")
    return specs


def interpolation_check(panel: pd.DataFrame) -> dict[str, Any]:
    """Median interpolated from a district's own B25104 bins against its published B25105
    median (districts with both): mean and SD of the log error."""
    d = panel[panel.year == YEAR]
    hcols = [f"bin_h_{j:02d}" for j in range(15)]
    ok = d[hcols].notna().all(axis=1) & d["h_all"].notna()
    gm = np.array([grouped_median(r, H_EDGES) for r in d.loc[ok, hcols].to_numpy()])
    err = np.log(gm) - np.log(d.loc[ok, "h_all"].to_numpy(float))
    err = err[np.isfinite(err)]
    return {"n": int(err.size), "mean_log_error": float(err.mean()),
            "sd_log_error": float(err.std(ddof=1)),
            "p01": float(np.percentile(err, 1)), "p99": float(np.percentile(err, 99))}


def run_scheme(panel: pd.DataFrame, scheme: str, B: int, seed: int) -> dict[str, Any]:
    LOG.info("== %s partition ==", scheme)
    part, g = db.cz_gap_table(panel, gap_defs(), scheme)
    g = add_pooled(part, g)
    w = g.pupils.to_numpy(float)
    counts = db.draw_counts(len(g), B, seed)
    LOG.info("  %d CZs, %d districts, B = %d", len(g), len(part), B)
    stats: dict[str, dict[str, Any]] = {}
    for sp in stat_specs():
        if sp["col"] not in g.columns or g[sp["col"]].notna().sum() < 5:
            continue
        st = stat_from_column(g, sp["col"], w, counts, sp["l1"], sp["l2"])
        st.update(kind=sp["kind"], label=sp["label"], family=sp["family"])
        if sp["kind"] == "usd":                 # model units: r_p = 12 * gap / CE expenditures
            f = 12.0 / CE_TOTAL_EXP
            st["rp"] = dict(value=st["value"] * f, boot_se=st["boot_se"] * f,
                            ci95_lo=st["ci95_lo"] * f, ci95_hi=st["ci95_hi"] * f)
            st["annual_dollars"] = 12.0 * st["value"]
        stats[sp["key"]] = st
    # share of the CZ-pooled FTYR earnings median
    ratios: dict[str, dict[str, Any]] = {}
    if "share_ftyr12_pm_h" in g:
        s1 = stat_from_column(g, "share_ftyr12_pm_h", w, counts, None, None)
        s1["label"] = ("12 x pooled-median dollar gap (a) / CZ-pooled median FTYR earnings, "
                       "computed within each CZ and averaged with the CZ weights")
        ratios["pm_h_over_ftyr12_mean_of_ratios"] = s1
        s2 = boot_ratio_of_means(g.gap_pm_h_usd.to_numpy(float), g.ftyr_12.to_numpy(float), w,
                                 counts, scale=12.0)
        s2["label"] = ("12 x (average pooled-median dollar gap) / (average CZ-pooled median "
                       "FTYR earnings): ratio of the two CZ-weighted averages")
        ratios["pm_h_over_ftyr12_ratio_of_means"] = s2
        for lab, col in (("city", "ftyr_1"), ("suburb", "ftyr_2")):
            r = boot_ratio_of_means(g.gap_pm_h_usd.to_numpy(float), g[col].to_numpy(float), w,
                                    counts, scale=12.0)
            r["label"] = (f"12 x (average pooled-median dollar gap) / (average {lab} median "
                          "FTYR earnings)")
            ratios[f"pm_h_over_ftyr_{lab}_ratio_of_means"] = r
        for lev in ("ftyr_1", "ftyr_2", "ftyr_12"):
            x = g[lev].to_numpy(float)
            m = np.isfinite(x) & (w > 0)
            ratios[f"level_{lev}"] = dict(value=float(np.average(x[m], weights=w[m])),
                                          n=int(m.sum()), label=f"CZ-weighted mean of {lev}")
        # CE-based r_p for the other dollar-gap measures also against FTYR earnings
    return dict(scheme=scheme, n_cz=int(len(g)), cz_weights_total_pupils=float(w.sum()),
                coverage=coverage(part), statistics=stats, ftyr_ratios=ratios, cz_table=g,
                counts_shape=list(counts.shape))


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def clean(o: Any) -> Any:
    """JSON-safe copy: NaN/inf -> None, numpy -> python, draws dropped."""
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items() if k not in ("draws", "cz_table")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


def f(v: float | None, d: int = 3) -> str:
    return "" if v is None or not np.isfinite(v) else f"{v:.{d}f}"


def fse(v: float | None, se: float | None, d: int = 3) -> str:
    return "" if v is None or not np.isfinite(v) else f"{v:,.{d}f} ({se:,.{d}f})"


def fd(v: float | None, d: int = 0) -> str:
    """Dollars with a thousands separator."""
    return "" if v is None or not np.isfinite(v) else f"{v:,.{d}f}"


SHORT = {
    "pm_h_usd": "pooled-location median, B25104 bins (a)",
    "usd_h_all": "average of district medians, B25105, occupied-unit weights (b)",
    "usd_h_rent": "average gross rent, B25064, occupied-unit weights (b)",
    "usd_h_rent_tw": "average gross rent, B25064, renter weights (b)",
    "usd_h_own": "average owner costs, B25088_001, occupied-unit weights (b)",
    "usd_h_own_tw": "average owner costs, B25088_001, owner weights (b)",
    "usd_h_own_mort": "average owner costs with a mortgage, B25088_002 (b)",
    "usd_h_own_nomort": "average owner costs without a mortgage, B25088_003 (b)",
    "pmean_h3500_usd": "pooled-location MEAN from bin midpoints, top bin $3,500",
    "pmean_h4000_usd": "pooled-location MEAN from bin midpoints, top bin $4,000",
}


def write_md(path: Path, meta: dict[str, Any], res: dict[str, dict[str, Any]],
             repro: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    loc, rev, ses = res["locale"], res.get("revenue"), res.get("ses")
    S = lambda r, k: r["statistics"].get(k) if r else None
    pm, pmu = S(loc, "pm_h_log"), S(loc, "pm_h_usd")
    parts = (("locale", loc), ("revenue", rev), ("ses", ses))

    a("# Housing costs by school district: suburb-minus-city gaps and the residence cost r_p\n")
    a(f"Generated {meta['generated']} by `housing_costs.py`. ACS 5-year 2014-18 school-district "
      "tables (`acs_housing.py`, dollars of 2018), merged to the district panel and the "
      f"CZ partitions of `district_bootstrap.py`; whole-CZ bootstrap, B = {meta['B']}, seed "
      f"{meta['seed']}.\n")

    a("## Summary\n")
    if pm and pmu:
        a(f"- **Locale partition (suburb minus city), pooled-location median monthly housing cost "
          f"(B25104 bins).** ${fd(pmu['level_1'])} in the city and ${fd(pmu['level_2'])} in "
          f"the suburb (CZ-weighted means of the CZ-level medians): a gap of **${fd(pmu['value'])} "
          f"a month (SE {fd(pmu['boot_se'])})**, **{f(pm['value'])} log points (SE "
          f"{f(pm['boot_se'])})**, over {pmu['n']} CZs.")
        rp = pmu["rp"]
        rl, ra = meta["rp_range_locale_core_measures"], meta["rp_range_core_all_partitions"]
        a(f"- **In model units.** r_p = 12 x gap / CE-2018 mean annual expenditures per consumer unit "
          f"(${CE_TOTAL_EXP:,.2f}, BLS) = **{f(rp['value'], 4)} (SE {f(rp['boot_se'], 4)})**, i.e. "
          f"{100 * rp['value']:.1f}% of Cbar. **Range across the all-tenure dollar measures, locale partition: "
          f"{rl[0]:.3f} ({rl[2]}) to {rl[1]:.3f} ({rl[3]})**; the pooled median sits inside. "
          f"Owner subsets (with / without a mortgage) are lower, see the scale table. "
          f"Pooled median under the revenue and ses partitions: "
          f"{f(S(rev, 'pm_h_usd')['rp']['value'], 4) if rev else ''} and "
          f"{f(S(ses, 'pm_h_usd')['rp']['value'], 4) if ses else ''}; across the same measures and "
          f"the three partitions {ra[0]:.3f} to {ra[1]:.3f}.")
    fr = loc["ftyr_ratios"]
    if "pm_h_over_ftyr12_mean_of_ratios" in fr:
        r1, r2 = fr["pm_h_over_ftyr12_mean_of_ratios"], fr["pm_h_over_ftyr12_ratio_of_means"]
        a(f"- **Relative to earnings.** 12 x the pooled-median dollar gap is "
          f"**{100 * r2['value']:.1f}% (SE {100 * r2['boot_se']:.1f}) of the CZ-pooled median "
          f"full-time year-round (FTYR) earnings** (${fd(fr['level_ftyr_12']['value'])}; ratio of "
          f"CZ-weighted averages), {100 * r1['value']:.1f}% (SE {100 * r1['boot_se']:.1f}) as the "
          "average of the within-CZ ratios.")
    o = S(loc, "own_share")
    ha, hr, ho, hv = (S(loc, "ln_h_all"), S(loc, "ln_h_rent"), S(loc, "ln_h_own"),
                      S(loc, "ln_h_value"))
    if all((o, ha, hr, ho, hv)):
        a(f"- **District medians, locale (log gaps, SE).** monthly housing costs (B25105) "
          f"{fse(ha['value'], ha['boot_se'])}; gross rent (B25064) {fse(hr['value'], hr['boot_se'])}; "
          f"owner costs (B25088) {fse(ho['value'], ho['boot_se'])}; home value (B25077) "
          f"{fse(hv['value'], hv['boot_se'])}. **Owner share** {f(o['level_1'])} in the city, "
          f"{f(o['level_2'])} in the suburb: gap {fse(o['value'], o['boot_se'])}.")
    a("- **Upper bound.** The whole level gap in housing spending is not a fixed cost: it mixes "
      "(i) a price per unit of housing, (ii) more or better housing (rooms, lot size, quality, "
      "tenure mix) in the suburb, which scales with income, and (iii) for owners, mortgage "
      "principal, which is saving. Treating the entire gap as the fixed, non-homothetic cost "
      "r_p x Cbar therefore gives an **upper bound** on that component; the homothetic part "
      "(richer households buy more housing wherever they live) belongs in the budget share, not "
      "in r_p.\n")

    ic = meta["interpolation_check"]
    a("## Data and method\n")
    a("- **Tables** (ACS 5-year 2014-18, school districts, three layers): B25105 median monthly "
      "housing costs; B25064 median gross rent; B25088 median selected monthly owner costs (all, "
      "with, without a mortgage); B25077 median home value; B25003 tenure (occupied, owner, renter "
      "units); B25104 distribution of monthly housing costs (15 cost bins from 'less than $100' to "
      "'$3,000 or more', plus 'No cash rent', which is excluded from the cost distribution). "
      "Census annotation codes (-666666666 etc.) are masked as missing; bottom-coded medians "
      "(<= $99, <= $9,999 of value) are dropped and top codes (4,001; 3,501 rent; 1,501 owner costs "
      "without a mortgage; 2,000,001 value) kept as lower bounds, as the earnings code does.")
    a("- **Merge and partition** exactly as the earnings tables in `district_bootstrap.py`: "
      "`leaid`, 'Remainder of <state>' rows dropped, `data_estimate.partition` (min 500 pupils, "
      "min 10% group share), `czone_gaps`, CZ weight = partitioned pupils. Within a location the "
      "median is pooled from summed bin counts (a), or the district medians are averaged with "
      "occupied-unit weights (b), or tenure counts are summed (c). Districts without a table are "
      "skipped.")
    a(f"- **Bin interpolation check.** The median interpolated from a district's B25104 bins "
      f"reproduces its published B25105 median with a mean log error of {f(ic['mean_log_error'], 4)} "
      f"(SD {f(ic['sd_log_error'], 4)}, {ic['n']} districts).")
    a("- **Bootstrap:** whole CZs resampled with replacement (`district_bootstrap.draw_counts`), "
      f"B = {meta['B']}, seed {meta['seed']}; one multiplicity matrix per partition shared by every "
      "statistic; SE = SD of the draws, interval = 2.5-97.5 percentiles. The locale draws are "
      "the ones `district_bootstrap.py` uses (same n and seed). Reproduction: rebuilt "
      f"`pm_ftyr_all` = {f(repro['value_here'], 5)} (SE {f(repro['boot_se_here'], 5)}) vs the "
      f"published {f(repro['value_published'], 5)} (SE {f(repro['boot_se_published'], 5)}).")
    a("- **Model units.** r_p = 12 x (monthly dollar gap) / (mean annual expenditures per consumer "
      f"unit, CE 2018) with the denominator **${CE_TOTAL_EXP:,.2f}**: BLS Consumer Expenditure "
      f"Survey, Table R-1 'All consumer units: annual detailed expenditure means', calendar year "
      f"2018 ({CE_UNITS_K:,} thousand consumer units, {CE_PERSONS} persons and {CE_EARNERS} earners "
      f"per unit, income before taxes ${CE_INCOME:,}; housing ${CE_HOUSING:,.2f}, of which shelter "
      f"${CE_SHELTER:,.2f}). Source: <{CE_URL}>; narrative report <{CE_REPORT_URL}>.\n")

    if pmu:
        alt = 12.0 * pmu["value"] / (CE_TOTAL_EXP - CE_PENSIONS)
        a(f"- **Denominator sensitivity.** CE total expenditures include personal insurance and "
          f"pensions (${CE_PENSIONS:,.2f}, {100 * CE_PENSIONS / CE_TOTAL_EXP:.1f}% of the total), "
          f"saving-like items; dropping them from the denominator (${CE_TOTAL_EXP - CE_PENSIONS:,.2f}) "
          f"raises every r_p by {100 * CE_PENSIONS / (CE_TOTAL_EXP - CE_PENSIONS):.1f}% (pooled-median, "
          f"locale: {alt:.4f}). For scale, 12 x the pooled-median gap is "
          f"{100 * 12 * pmu['value'] / CE_SHELTER:.1f}% of CE shelter spending (${CE_SHELTER:,.2f}).")

    def cover(r):
        c = r["coverage"]
        return (f"{r['n_cz']} CZs, {c['n_districts']} districts; ACS housing tables matched for "
                f"{100 * c['share_districts_with_housing_tables']:.1f}% of districts "
                f"({100 * c['share_pupils_with_housing_tables']:.1f}% of pupils)")

    for tag, r in parts:
        if r:
            a(f"- `{tag}` ({SCHEME_DESC[tag]}): {cover(r)}.")
    a("")

    a("## (a) Pooled-location median monthly housing cost (B25104)\n")
    a("Location levels are CZ-weighted means of the CZ-level pooled medians (dollars of 2018 per "
      "month); the gap is the CZ-weighted mean of the within-CZ gap (location 2 minus location 1; "
      "city/low = 1).\n")
    a("| partition | loc 1 $ | loc 2 $ | gap $ (SE) | gap log (SE) | 95% interval (log) | r_p = 12 x gap / "
      "CE (SE) | 12 x gap / FTYR earnings % (SE) | n CZs |")
    a("|---|---|---|---|---|---|---|---|---|")
    for tag, r in parts:
        if not r:
            continue
        u, lg = S(r, "pm_h_usd"), S(r, "pm_h_log")
        rr = r["ftyr_ratios"].get("pm_h_over_ftyr12_ratio_of_means")
        a(f"| {tag} | {fd(u['level_1'])} | {fd(u['level_2'])} | {fse(u['value'], u['boot_se'], 0)} | "
          f"{fse(lg['value'], lg['boot_se'])} | [{f(lg['ci95_lo'])}, {f(lg['ci95_hi'])}] | "
          f"{fse(u['rp']['value'], u['rp']['boot_se'], 4)} | "
          f"{fse(100 * rr['value'], 100 * rr['boot_se'], 1) if rr else ''} | {u['n']} |")
    a("")
    a("Sensitivity, pooled-location MEAN monthly housing cost from bin midpoints (top bin '$3,000 or "
      "more' valued at $3,500 / $4,000). The CE denominator is a mean, but the median gap is the "
      "better proxy for a cost every household pays; the mean gap also contains the upper tail "
      "of larger suburban homes.\n")
    a("| partition | top bin | gap $ (SE) | gap log (SE) | r_p (SE) |")
    a("|---|---|---|---|---|")
    for tag, r in parts:
        if not r:
            continue
        for top in TOP_BIN_VALUES:
            u, lg = S(r, f"pmean_h{top}_usd"), S(r, f"pmean_h{top}_log")
            if u:
                a(f"| {tag} | ${top:,} | {fse(u['value'], u['boot_se'], 0)} | "
                  f"{fse(lg['value'], lg['boot_se'])} | {fse(u['rp']['value'], u['rp']['boot_se'], 4)} |")
    a("")

    a("## (b) Occupied-unit-weighted means of district medians\n")
    a("Weights B25003_001 (occupied units) within a location; CZ weights = pupils. `ln` = gap in the "
      "average log median (the earnings convention), `lvl` = log ratio of the average level, `$` = "
      "gap in the average level; `r_p` as above (not defined for the home-value stock). Weights "
      "'tenure' (`_tw` in the JSON) weight rent by renter units and owner costs and value by owner "
      "units.\n")
    a("| series | partition | weights | ln gap (SE) | lvl gap (SE) | $ gap (SE) | loc 1 $ | loc 2 $ | "
      "r_p (SE) | n CZs |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for k, (col, lab, tenure, _) in SERIES.items():
        for tag in ("", "_tw"):
            for tag2, r in parts:
                if not r:
                    continue
                ln, lv, us = S(r, f"ln_{k}{tag}"), S(r, f"lvl_{k}{tag}"), S(r, f"usd_{k}{tag}")
                if not ln or not us:
                    continue
                rp = us.get("rp")
                a(f"| `{col[:-1]}` | {tag2} | {'tenure' if tag else 'occupied'} | "
                  f"{fse(ln['value'], ln['boot_se'])} | {fse(lv['value'], lv['boot_se'])} | "
                  f"{fse(us['value'], us['boot_se'], 0)} | {fd(us.get('level_1'))} | "
                  f"{fd(us.get('level_2'))} | "
                  f"{fse(rp['value'], rp['boot_se'], 4) if rp else ''} | {ln['n']} |")
    a("")

    a("## (c) Owner share\n")
    a("Owner-occupied over occupied units (B25003), pooled over the districts of a location.\n")
    a("| partition | loc 1 | loc 2 | gap (SE) | 95% interval | n CZs |")
    a("|---|---|---|---|---|---|")
    for tag, r in parts:
        if not r:
            continue
        o = S(r, "own_share")
        a(f"| {tag} | {f(o['level_1'])} | {f(o['level_2'])} | {fse(o['value'], o['boot_se'])} | "
          f"[{f(o['ci95_lo'])}, {f(o['ci95_hi'])}] | {o['n']} |")
    a("")

    a("## Scale: r_p from each dollar measure (locale partition)\n")
    a("| measure | $ gap per month (SE) | 12 x gap per year | r_p (SE) | 95% interval for r_p |")
    a("|---|---|---|---|---|")
    for k in meta["rp_measures"]:
        s_ = S(loc, k)
        if s_ and "rp" in s_:
            a(f"| {SHORT.get(k, k)} | {fse(s_['value'], s_['boot_se'], 0)} | "
              f"{fd(s_['annual_dollars'])} | {fse(s_['rp']['value'], s_['rp']['boot_se'], 4)} | "
              f"[{f(s_['rp']['ci95_lo'], 4)}, {f(s_['rp']['ci95_hi'], 4)}] |")
    a("")
    if fr:
        a("Pooled-median dollar gap against the earnings scale (locale; B20005 FTYR bins), CZ-weighted "
          "means with whole-CZ bootstrap SEs:\n")
        a("| denominator | 12 x gap / earnings % (SE) | 95% interval % | n CZs |")
        a("|---|---|---|---|")
        for k, v in fr.items():
            if k.startswith("pm_h_over") and "boot_se" in v:
                a(f"| {v['label']} | {fse(100 * v['value'], 100 * v['boot_se'], 1)} | "
                  f"[{100 * v['ci95_lo']:.1f}, {100 * v['ci95_hi']:.1f}] | {v['n']} |")
        a("")
        a("Earnings levels (CZ-weighted means, dollars of 2018): median FTYR earnings city "
          f"${fd(fr['level_ftyr_1']['value'])}, suburb ${fd(fr['level_ftyr_2']['value'])}, the two "
          f"pooled ${fd(fr['level_ftyr_12']['value'])}.\n")

    a("## Caveats\n")
    a("- **Upper bound, not a point estimate of the fixed cost.** The gap in total housing "
      "*spending* bundles the price of a unit of housing, the quantity and quality of housing "
      "(suburban homes have more rooms and land), tenure (owner costs include mortgage principal, "
      "which is saving, and property taxes, part of which pay for the local schools) "
      "and the income gradient of housing demand. Only the first is the model's p_2 - p_1; "
      "the homothetic part of the second is in the preferences. Hedonic or repeat-sales price "
      "indices would give a smaller number; district medians cannot.")
    a("- **Units.** B25105/B25064/B25088 are per housing unit (household); the CE denominator is per "
      f"consumer unit ({CE_PERSONS} persons, {CE_EARNERS} earners), the model's Cbar is per adult. The "
      "ratio is a unit-free share if housing cost and consumption scale together with household "
      "size, which is the assumption used here.")
    a("- **Medians against a mean.** The numerator is a median (or an average of medians) and the "
      "denominator a mean; the share of expenditures would be larger against a median "
      "expenditure. The mean-based sensitivity in (a) brackets the numerator side.")
    a("- **Timing and geography.** ACS 2014-18 (5-year pooled, dollars of 2018) merged to the 2018 CCD "
      "panel by district id. Locations are groups of whole districts inside a 1990 commuting zone, "
      "so a 'city' district mixes in suburban neighborhoods and vice versa; the gaps are within-CZ "
      "averages weighted by pupils, so they are dominated by large CZs (see `n_eff_cz` in the JSON).")
    a("- **Revenue and ses partitions** compare the advantaged with the disadvantaged half of a CZ, "
      "not the suburb with the city; their housing gaps are price differences between rich and "
      "poor school districts.")
    a("- **B25077 (home value)** is a stock; it is shown for reference and not converted to r_p.")
    path.write_text("\n".join(L) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(HERE))
    ap.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    ap.add_argument("--B", type=int, default=B_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--schemes", nargs="+", default=list(SCHEMES), choices=SCHEMES)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates" / "settings"
    outdir.mkdir(parents=True, exist_ok=True)

    panel = db.build_panel(indir)
    panel = db.load_income_columns(panel, indir / "raw")
    panel = load_housing_columns(panel, indir / "raw")
    LOG.info("panel: %d districts; ACS housing medians matched: %d", len(panel),
             int(panel["h_all"].notna().sum()))
    icheck = interpolation_check(panel)
    LOG.info("bin interpolation vs published B25105: mean log error %.4f (SD %.4f, n = %d)",
             icheck["mean_log_error"], icheck["sd_log_error"], icheck["n"])

    res: dict[str, dict[str, Any]] = {}
    for sch in args.schemes:
        res[sch] = run_scheme(panel, sch, args.B, args.seed)

    # reproduction check against district_bootstrap.json (same CZ table, same draws)
    repro: dict[str, Any] = {"value_here": np.nan, "boot_se_here": np.nan,
                             "value_published": np.nan, "boot_se_published": np.nan}
    pub = indir / "estimates" / "district_bootstrap.json"
    if "locale" in res and pub.exists() and "pm_ftyr_all" in res["locale"]["statistics"]:
        p = json.loads(pub.read_text())["statistics"]["pm_ftyr_all"]
        h = res["locale"]["statistics"]["pm_ftyr_all"]
        repro = {"value_here": h["value"], "boot_se_here": h["boot_se"],
                 "value_published": p["value"], "boot_se_published": p["boot_se"],
                 "n_here": h["n"], "n_published": p["n"],
                 "max_abs_diff_value": abs(h["value"] - p["value"]),
                 "max_abs_diff_se": abs(h["boot_se"] - p["boot_se"]),
                 "same_B_and_seed": bool(args.B == db.B_DEFAULT and args.seed == db.SEED)}
        LOG.info("reproduction of pm_ftyr_all: here %.6f (SE %.6f), published %.6f (SE %.6f)",
                 h["value"], h["boot_se"], p["value"], p["boot_se"])

    # r_p ranges over the dollar measures of monthly housing cost (not the home-value stock):
    # `core` = every all-tenure measure of the cost a typical household pays (pooled medians and
    # means, averages of district medians, rent, owner costs); `all` adds the owner subsets.
    rp_core = ["pm_h_usd", "usd_h_all", "usd_h_rent", "usd_h_rent_tw", "usd_h_own",
               "usd_h_own_tw", "pmean_h3500_usd", "pmean_h4000_usd"]
    rp_measures = rp_core + ["usd_h_own_mort", "usd_h_own_nomort"]

    def rp_vals(r, keys):
        return [(r["statistics"][k]["rp"]["value"], k) for k in keys
                if k in r["statistics"] and "rp" in r["statistics"][k]]

    core_loc, all_loc = rp_vals(res["locale"], rp_core), rp_vals(res["locale"], rp_measures)
    lo, hi = min(core_loc), max(core_loc)
    all_vals = [v for r in res.values() for v, _ in rp_vals(r, rp_core)]
    meta = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "data/spatial/housing_costs.py", "year": YEAR, "B": args.B, "seed": args.seed,
        "rng": "numpy default_rng (PCG64) via district_bootstrap.draw_counts; one (B x n_CZ) "
               "multiplicity matrix per partition shared by every statistic; resampling unit = whole CZ",
        "geo_unit": panel.attrs.get("geo_unit"), "min_pupils": 500, "min_group_share": 0.10,
        "cz_weights": "partitioned pupils of the CZ",
        "tables": {"B25105_001": "median monthly housing costs", "B25064_001": "median gross rent",
                   "B25088_001/002/003": "median selected monthly owner costs (all / mortgage / none)",
                   "B25077_001": "median value", "B25003_001/002/003": "occupied / owner / renter units",
                   "B25104_001-017": "monthly housing cost bins (counts)"},
        "ce_2018": {"average_annual_expenditures_per_consumer_unit": CE_TOTAL_EXP,
                    "housing": CE_HOUSING, "shelter": CE_SHELTER,
                    "personal_insurance_and_pensions": CE_PENSIONS,
                    "consumer_units_thousands": CE_UNITS_K, "persons_per_unit": CE_PERSONS,
                    "earners_per_unit": CE_EARNERS, "income_before_taxes": CE_INCOME,
                    "source_table": CE_URL, "source_report": CE_REPORT_URL,
                    "source": "BLS Consumer Expenditure Surveys, Table R-1 'All consumer units: "
                              "Annual detailed expenditure means', 2018"},
        "rp_definition": "r_p = 12 x (dollar gap in monthly housing cost, suburb - city) / "
                         "(mean annual expenditures per consumer unit, CE 2018)",
        "upper_bound_note": "The whole level gap in housing spending includes quantity/quality "
                            "and income-scaled (homothetic) differences and mortgage principal; "
                            "treating it all as the fixed (non-homothetic) cost is an upper bound.",
        "rp_measures": rp_measures, "rp_core_measures": rp_core,
        "rp_range_locale_core_measures": [lo[0], hi[0], SHORT[lo[1]], SHORT[hi[1]]],
        "rp_range_locale_incl_owner_subsets": [min(all_loc)[0], max(all_loc)[0],
                                               SHORT[min(all_loc)[1]], SHORT[max(all_loc)[1]]],
        "rp_range_core_all_partitions": [min(all_vals), max(all_vals)],
        "interpolation_check": icheck,
        "reproduction_vs_district_bootstrap": repro,
    }
    payload = clean({"meta": meta, "partitions": {k: {kk: vv for kk, vv in v.items()
                                                      if kk != "cz_table"} for k, v in res.items()}})
    (outdir / "housing_costs.json").write_text(json.dumps(payload, indent=1))
    write_md(outdir / "housing_costs.md", meta, res, repro)

    if "locale" in res:       # locale draws (same draws as district_bootstrap.py), main statistics
        D = {k: res["locale"]["statistics"][k]["draws"] for k in
             ("pm_h_log", "pm_h_usd", "usd_h_all", "usd_h_rent", "usd_h_own", "ln_h_all",
              "own_share", "pm_ftyr_all") if k in res["locale"]["statistics"]}
        pd.DataFrame(D).round(8).to_csv(outdir / "housing_costs_draws.csv", index_label="draw")
    LOG.info("wrote housing_costs.{json,md} and housing_costs_draws.csv to %s", outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
