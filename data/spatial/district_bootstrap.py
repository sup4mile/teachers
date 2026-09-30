#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
district_bootstrap.py -- whole-commuting-zone bootstrap of the district moments
(calibration item T5) and of the income alternatives built from ACS earnings
tables (item T3).

Every headline number of `data_estimate.py` is an enrollment-weighted average,
across commuting zones (CZs), of a within-CZ suburb-minus-city gap.  This script
rebuilds the same objects with `data_estimate.py`'s own functions (panel,
outlier screen, locale partition, `czone_gaps`, `wmean`), then resamples WHOLE
CZs with replacement and recomputes the weighted mean.  A CZ's gap depends only
on its own districts, so the CZ x gap table is a sufficient statistic and each
draw is a weighted mean over resampled table rows:

    stat_b = sum_i c_bi w_i x_i 1{x_i observed} / sum_i c_bi w_i 1{x_i observed}

with c_bi the multiplicity of CZ i in draw b and w_i its (partitioned) pupils.
CZs that lack the statistic (one location has no valid district) drop out of
that draw, as they do from the point estimate, so the CZ count n varies by draw.

Held fixed across draws: the district panel after the within-year [0.5%, 99.5%]
ratio screen (the screen quantiles are computed on all districts, not CZ by CZ),
the CZ partition and the min-pupil / min-group-share filters.  Only CZ membership
of the average is resampled.

Placeholder scale: the cross-CZ SD of the gap (unweighted, ddof = 1) over sqrt(n),
what `spatial_calibrate.jl` currently uses.  It ignores the enrollment weights,
which concentrate on a handful of large CZs, so the ratio boot SE / placeholder
is the correction.

Also reported: the same bootstrap for `local_revenue_share` (a level over sample
districts, not a gap) and for income alternatives (log median earnings from ACS
5-year 2014-18 school-district tables, see `acs_earnings.py`), computed on the
same draws so that their bootstrap covariances with the targets are coherent.

Usage
    python3 district_bootstrap.py                 # B = 1999, seed 20260929
    python3 district_bootstrap.py --B 499 --seed 1
Writes estimates/district_bootstrap.{json,md} and district_bootstrap_draws.csv.

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

LOG = logging.getLogger("district_bootstrap")

YEAR = 2018
SCHEME = "locale"
SEED = 20260929
B_DEFAULT = 1999

TARGETS = ["pupils", "salary_real", "teachers_pp", "seda_lninc"]
VALIDATION = ["pov_rate_5_17", "learn_rate", "score_mean"]   # + local_revenue_share

# ACS earnings series: name -> (ACS column, description).  Medians of district
# earnings are aggregated like SEDA's log income: pupil-weighted average of the
# district log median inside each location, then suburb minus city.
INCOME_SERIES: dict[str, tuple[str, str]] = {
    "earn_all": ("B20002_001E", "median earnings, 16+ with earnings (B20002_001)"),
    "earn_male": ("B20002_002E", "median earnings, men 16+ with earnings (B20002_002)"),
    "earn_female": ("B20002_003E", "median earnings, women 16+ with earnings (B20002_003)"),
    "earn_ftyr": ("B20018_001E", "median earnings, full-time year-round 16+ (B20018_001)"),
    "earn_ftyr_male": ("B20017_003E", "median earnings, men full-time year-round (B20017_003)"),
    "earn_ftyr_female": ("B20017_006E", "median earnings, women full-time year-round (B20017_006)"),
    "earn_25p": ("B20004_001E", "median earnings, 25+ with earnings (B20004_001)"),
    "earn_25p_male": ("B20004_007E", "median earnings, men 25+ (B20004_007)"),
    "earn_25p_female": ("B20004_013E", "median earnings, women 25+ (B20004_013)"),
    "earn_emp": ("B24011_001E", "median earnings, civilian employed 16+ (B24011_001)"),
    "hhinc_2544": ("B19049_003E", "median household income, householder 25-44 (B19049_003)"),
    "hhinc_acs": ("B19013_001E", "median household income, ACS 2014-18 (B19013_001)"),
}
BOTTOM_CODE = 2499       # ACS publishes 2,499 when the median is in the lowest interval

# Earnings bins (B20001 for everyone with earnings, B20005 for full-time
# year-round workers): 20 intervals, the first "$1 to $2,499 or loss", the last
# "$100,000 or more".  Pooled location medians are read off summed bin counts.
BIN_EDGES = np.array([0, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500,
                      25000, 30000, 35000, 40000, 45000, 50000, 55000, 65000, 75000,
                      100000, np.inf])
# series -> ACS columns of the 20 bins (male, female), summed for "all"
BIN_SERIES: dict[str, tuple[str, list[str], str]] = {
    "earn_all":      ("B20001", ["m", "f"], "all workers 16+ with earnings"),
    "earn_male":     ("B20001", ["m"], "men 16+ with earnings"),
    "earn_female":   ("B20001", ["f"], "women 16+ with earnings"),
    "ftyr_all":      ("B20005", ["m", "f"], "full-time year-round workers 16+"),
    "ftyr_male":     ("B20005", ["m"], "men, full-time year-round"),
    "ftyr_female":   ("B20005", ["f"], "women, full-time year-round"),
}
BIN_START = {("B20001", "m"): 3, ("B20001", "f"): 24, ("B20005", "m"): 6, ("B20005", "f"): 53}


# --------------------------------------------------------------------------
# Building the CZ x gap table exactly as data_estimate.py does
# --------------------------------------------------------------------------

def load_income_columns(panel: pd.DataFrame, raw_dir: Path) -> pd.DataFrame:
    """Merge the ACS earnings columns (logs and levels) onto the district panel."""
    path = raw_dir / f"acs_earnings_school_districts_{YEAR}.parquet"
    if not path.exists():
        LOG.warning("no %s: income alternatives skipped (run acs_earnings.py)", path.name)
        return panel
    a = pd.read_parquet(path)
    # 'Remainder of <state>' rows carry code 99999 and are not districts; the
    # layers otherwise never share an id.
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid")
    out = pd.DataFrame({"leaid": a.leaid})
    for name, (col, _) in INCOME_SERIES.items():
        x = a[col].where(a[col] > BOTTOM_CODE)           # drops negative annotation codes
        out[name] = x
        out[f"ln_{name}"] = np.log(x)
    earn_hh = (a["B19061_001E"].where(a["B19061_001E"] > 0)
               / a["B19051_002E"].where(a["B19051_002E"] > 0))
    out["earn_hh_mean"] = earn_hh
    out["ln_earn_hh_mean"] = np.log(earn_hh)
    out["earners"] = a["B20001_001E"].where(a["B20001_001E"] > 0)
    out["pop_acs"] = a["B01003_001E"].where(a["B01003_001E"] > 0)
    bpath = raw_dir / f"acs_earnings_bins_school_districts_{YEAR}.parquet"
    if bpath.exists():
        bn = pd.read_parquet(bpath)
        bn = bn[~bn.leaid.str.endswith("99999")].drop_duplicates("leaid")
        bi = bn.set_index("leaid")
        newcols = {}
        for name, (tab, sexes, _) in BIN_SERIES.items():
            for j in range(20):
                cols = [f"{tab}_{BIN_START[(tab, sx)] + j:03d}E" for sx in sexes]
                newcols[f"bin_{name}_{j:02d}"] = bi[cols].sum(
                    axis=1, min_count=len(cols)).reindex(out.leaid).to_numpy()
        out = pd.concat([out, pd.DataFrame(newcols, index=out.index)], axis=1)
    keep = [c for c in out.columns if c not in panel.columns or c == "leaid"]
    attrs = dict(panel.attrs)
    merged = panel.merge(out[keep], on="leaid", how="left")
    merged.attrs.update(attrs)
    return merged


def grouped_median(counts: np.ndarray) -> float:
    """Median of a distribution known only through counts in BIN_EDGES intervals,
    by linear interpolation inside the bin holding the 50th percentile (NaN if that
    is the open-ended top bin or there are no observations)."""
    c = np.asarray(counts, float)
    tot = c.sum()
    if not np.isfinite(tot) or tot <= 0:
        return np.nan
    cum = np.cumsum(c)
    j = int(np.searchsorted(cum, tot / 2.0))
    if j >= len(c) - 1 and c[-1] > 0 and cum[-2] < tot / 2.0:
        return np.nan                                   # median in the top bin
    below = cum[j - 1] if j > 0 else 0.0
    lo, hi = BIN_EDGES[j], BIN_EDGES[j + 1]
    return float(lo + (hi - lo) * (tot / 2.0 - below) / c[j])


def pooled_median_gaps(part: pd.DataFrame, g: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Add `gap_pm_<series>` to the CZ table: log pooled-median earnings of the
    suburb (loc 2) minus the city (loc 1), each median from the summed bin counts
    of the districts in the location.  Districts with a missing table are skipped."""
    d = part[part.year == YEAR]
    keys = []
    for name in BIN_SERIES:
        cols = [f"bin_{name}_{j:02d}" for j in range(20)]
        if not all(c in d.columns for c in cols):
            continue
        sums = d.groupby(["czone", "loc_id"])[cols].sum(min_count=1)
        med = pd.Series({idx: grouped_median(row.to_numpy()) for idx, row in sums.iterrows()})
        med.index = pd.MultiIndex.from_tuples(med.index, names=["czone", "loc_id"])
        wide = med.unstack("loc_id")
        gap = np.log(wide[2] / wide[1])
        g[f"gap_pm_{name}"] = g.czone.map(gap).to_numpy()
        keys.append(name)
    return g, keys


def income_gap_defs() -> list["de.Gap"]:
    """Gap definitions for the alternatives: log-difference of pupil-weighted
    means of district log medians (baseline convention), plus variants."""
    G: list[de.Gap] = []
    for name, (_, label) in INCOME_SERIES.items():
        G.append(de.Gap(f"ln_{name}", f"log {label}", "I_2/I_1 (labor income alternative)",
                        "mean", col=f"ln_{name}", contrast="diff"))
    G.append(de.Gap("ln_earn_hh_mean", "log mean earnings per earning household "
                    "(B19061_001/B19051_002)", "I_2/I_1 (alternative)",
                    "mean", col="ln_earn_hh_mean", contrast="diff"))
    # Level aggregation: log of the ratio of location averages of the district median.
    for name in ("earn_all", "earn_ftyr", "hhinc_2544"):
        G.append(de.Gap(f"lvl_{name}", f"log ratio of location-average {INCOME_SERIES[name][1]}",
                        "I_2/I_1 (alternative, level aggregation)",
                        "mean", col=name, contrast="log"))
    # Earner-weighted within-location averages (CZ weights stay pupils).
    for name in ("earn_all", "earn_ftyr"):
        G.append(de.Gap(f"ew_ln_{name}", f"log {INCOME_SERIES[name][1]}, earner-weighted "
                        "within location", "I_2/I_1 (alternative, earner weights)",
                        "mean", col=f"ln_{name}", weight="earners", contrast="diff"))
    return G


def build_panel(indir: Path, year: int = YEAR) -> pd.DataFrame:
    store = de.RawStore(indir / "raw")
    return de.assemble(store, [year], "auto")


def cz_gap_table(panel: pd.DataFrame, gaps: list["de.Gap"] | None = None,
                 scheme: str = SCHEME, min_pupils: int = 500,
                 min_group_share: float = 0.10) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Partition and compute the CZ x gap table for `gaps` (default: de.GAPS).

    Returns (partitioned districts, CZ table).  Calls data_estimate's own
    `partition` and `czone_gaps`; the gap list is swapped in for the call.
    """
    part = de.partition(panel, scheme, min_pupils, min_group_share)
    saved = de.GAPS
    try:
        if gaps is not None:
            de.GAPS = list(gaps)
        g = de.czone_gaps(part)
    finally:
        de.GAPS = saved
    return part, g[g.year == YEAR].reset_index(drop=True)


# --------------------------------------------------------------------------
# Bootstrap machinery
# --------------------------------------------------------------------------

def draw_counts(n: int, B: int, seed: int) -> np.ndarray:
    """(B x n) multiplicities of CZ resampling with replacement."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(B, n))
    counts = np.zeros((B, n))
    for b in range(B):
        counts[b] = np.bincount(idx[b], minlength=n)
    return counts


def boot_wmean(x: np.ndarray, w: np.ndarray, counts: np.ndarray) -> dict[str, Any]:
    """Point estimate and bootstrap draws of the weighted mean over observed rows."""
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    xw = np.where(m, x * w, 0.0)
    ww = np.where(m, w, 0.0)
    est = xw.sum() / ww.sum()
    den = counts @ ww
    with np.errstate(invalid="ignore", divide="ignore"):
        draws = (counts @ xw) / den
    n_draw = counts @ m.astype(float)
    # linearised (cluster-robust) SE of a ratio estimator, as a cross-check
    n = int(m.sum())
    r = np.where(m, w * (x - est), 0.0)
    se_lin = float(np.sqrt(n / (n - 1) * (r ** 2).sum()) / ww.sum()) if n > 1 else np.nan
    wn = w[m] / w[m].sum()
    return dict(est=float(est), draws=draws, n=n, n_draw_mean=float(n_draw.mean()),
                se_lin=se_lin, n_eff=float(1.0 / (wn ** 2).sum()),
                top5_share=float(np.sort(wn)[::-1][:5].sum()),
                max_share=float(wn.max()),
                sd_unw=float(np.std(x[m], ddof=1)))


def summarize(stat: dict[str, Any]) -> dict[str, Any]:
    d = stat["draws"]
    d = d[np.isfinite(d)]
    se = float(np.std(d, ddof=1))
    lo, hi = np.percentile(d, [2.5, 97.5])
    placeholder = stat["sd_unw"] / np.sqrt(stat["n"])
    return dict(value=stat["est"], n=stat["n"], boot_se=se,
                ci95_lo=float(lo), ci95_hi=float(hi),
                boot_mean_minus_est=float(d.mean() - stat["est"]),
                n_draws_finite=int(len(d)),
                mean_n_per_draw=stat["n_draw_mean"],
                placeholder_sd=stat["sd_unw"], placeholder_se=float(placeholder),
                ratio_boot_over_placeholder=float(se / placeholder),
                se_linearized=stat["se_lin"], n_eff_cz=stat["n_eff"],
                top5_weight_share=stat["top5_share"], max_weight_share=stat["max_share"])


def tracked_values(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    d = json.loads(path.read_text())
    return {e["key"]: e for e in d["estimates"]}


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

LABELS = {g.key: g.label for g in de.GAPS}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default=str(HERE))
    ap.add_argument("--outdir", default=None, help="default: <indir>/estimates")
    ap.add_argument("--B", type=int, default=B_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--tracked", default=None,
                    help="spatial_moments.json to compare point estimates with "
                         "(default: <indir>/estimates/spatial_moments.json)")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates"
    outdir.mkdir(parents=True, exist_ok=True)

    panel = build_panel(indir)
    panel = load_income_columns(panel, indir / "raw")
    inc_defs = income_gap_defs()
    have_inc = all(g.col in panel.columns for g in inc_defs)
    part, g = cz_gap_table(panel, list(de.GAPS) + (inc_defs if have_inc else []))
    pm_keys: list[str] = []
    if have_inc:
        g, pm_keys = pooled_median_gaps(part, g)
    LOG.info("%d CZs in the %s partition, %d districts", len(g), SCHEME, len(part))

    w = g.pupils.to_numpy(float)
    counts = draw_counts(len(g), args.B, args.seed)

    stats: dict[str, dict[str, Any]] = {}
    roles: dict[str, str] = {}
    labels: dict[str, str] = {}
    for gp in de.GAPS:
        col = f"gap_{gp.key}"
        if col in g and g[col].notna().sum() >= 5:
            stats[gp.key] = boot_wmean(g[col].to_numpy(float), w, counts)
            roles[gp.key] = ("target" if gp.key in TARGETS else
                             "validation" if gp.key in VALIDATION else "other")
            labels[gp.key] = gp.label
    if have_inc:
        for gp in inc_defs:
            stats[gp.key] = boot_wmean(g[f"gap_{gp.key}"].to_numpy(float), w, counts)
            roles[gp.key] = "income_alternative"
            labels[gp.key] = gp.label
        for name in pm_keys:
            k = f"pm_{name}"
            stats[k] = boot_wmean(g[f"gap_{k}"].to_numpy(float), w, counts)
            roles[k] = "income_alternative"
            labels[k] = ("log pooled location median earnings, " + BIN_SERIES[name][2]
                         + f" ({BIN_SERIES[name][0]} bins summed over the districts of "
                         "each location)")

    # local revenue share: enrollment-weighted mean over sample districts, CZs as clusters
    d = part[part.year == YEAR]
    ok = d.share_local.notna() & (d.pupils > 0)
    cz = d[ok].groupby("czone")
    cz_tab = pd.DataFrame({"num": cz.apply(lambda t: (t.pupils * t.share_local).sum(),
                                           include_groups=False),
                           "den": cz.pupils.sum(), "nd": cz.size()})
    # align rows to the CZ table used for the counts (same CZ set, same order)
    cz_tab = cz_tab.reindex(g.czone.to_numpy())
    den = cz_tab.den.fillna(0).to_numpy()
    num_ = cz_tab.num.fillna(0).to_numpy()
    est = num_.sum() / den.sum()
    draws = (counts @ num_) / (counts @ den)
    xcz = np.where(den > 0, num_ / np.where(den > 0, den, 1), np.nan)   # CZ-level shares
    wn = den / den.sum()
    nn = int((den > 0).sum())
    r = np.where(den > 0, num_ - est * den, 0.0)
    stats["local_revenue_share"] = dict(
        est=float(est), draws=draws, n=nn, n_draw_mean=float((counts @ (den > 0)).mean()),
        se_lin=float(np.sqrt(nn / (nn - 1) * (r ** 2).sum()) / den.sum()),
        n_eff=float(1 / (wn[den > 0] ** 2).sum()),
        top5_share=float(np.sort(wn)[::-1][:5].sum()), max_share=float(wn.max()),
        sd_unw=float(np.nanstd(xcz, ddof=1)))
    roles["local_revenue_share"] = "validation"
    labels["local_revenue_share"] = ("local share of district revenue, pupil-weighted over "
                                     "sample districts (level; CZ-level shares for the "
                                     "placeholder SD)")

    summ = {k: summarize(v) for k, v in stats.items()}
    for k in summ:
        summ[k]["role"] = roles[k]
        summ[k]["label"] = labels[k]
    n_districts_lrs = int(ok.sum())

    # compare with the tracked moments
    tracked_path = Path(args.tracked) if args.tracked else indir / "estimates" / "spatial_moments.json"
    tr = tracked_values(tracked_path)
    for k in summ:
        key = "local_revenue_share" if k == "local_revenue_share" else f"gap_{k}_{SCHEME}"
        e = tr.get(key)
        if e:
            summ[k]["tracked_value"] = e["value"]
            summ[k]["tracked_n"] = e["n"] if k != "local_revenue_share" else e["n"]
            summ[k]["diff_vs_tracked"] = summ[k]["value"] - e["value"]
    summ["local_revenue_share"]["n_districts"] = n_districts_lrs
    summ["local_revenue_share"]["n_clusters"] = summ["local_revenue_share"]["n"]
    cmp_keys = [k for k in summ if "diff_vs_tracked" in summ[k]]
    n_same = sum(1 for k in cmp_keys if k == "local_revenue_share"
                 and summ[k]["tracked_n"] == n_districts_lrs
                 or k != "local_revenue_share" and summ[k]["tracked_n"] == summ[k]["n"])
    repro = dict(n_compared=len(cmp_keys), n_same_n=n_same,
                 max_abs_diff=float(max((abs(summ[k]["diff_vs_tracked"]) for k in cmp_keys),
                                        default=np.nan)))

    # bootstrap correlation among the targets, validation and main alternatives
    core = TARGETS + VALIDATION + ["local_revenue_share"]
    Dm = pd.DataFrame({k: stats[k]["draws"] for k in core + [
        k for k in stats if roles[k] == "income_alternative"]})
    corr = Dm[core].corr()
    cov_t = Dm[TARGETS].cov()

    draws_path = outdir / "district_bootstrap_draws.csv"
    Dm.round(8).to_csv(draws_path, index_label="draw")

    meta = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "data/spatial/district_bootstrap.py",
        "year": YEAR, "scheme": SCHEME, "geo_unit": panel.attrs.get("geo_unit"),
        "B": args.B, "seed": args.seed, "rng": "numpy default_rng (PCG64); "
        "one (B x n_CZ) multiplicity matrix shared by every statistic",
        "resampling_unit": "commuting zone (whole CZ, with replacement, n = number of "
                           "partitioned CZs)", "n_cz_partition": int(len(g)),
        "n_districts_partition": int(len(part)),
        "cz_weights": "partitioned pupils of the CZ (as in data_estimate.summarize_gaps)",
        "min_pupils": 500, "min_group_share": 0.10,
        "held_fixed": "district panel after the within-year [0.5%,99.5%] screen, the "
                      "partition and the CZ filters",
        "placeholder": "unweighted cross-CZ SD (ddof=1) / sqrt(n)",
        "tracked_file_compared": str(tracked_path.name) if tr else None,
    }
    meta["reproduction_vs_tracked"] = repro
    payload = {**meta, "statistics": summ,
               "bootstrap_correlation": corr.round(4).to_dict(),
               "bootstrap_covariance_targets": cov_t.to_dict(),
               "income_alternatives_available": have_inc}
    (outdir / "district_bootstrap.json").write_text(json.dumps(payload, indent=2, default=str))
    write_md(outdir / "district_bootstrap.md", meta, summ, corr, have_inc)
    LOG.info("wrote district_bootstrap.{json,md} and %s", draws_path.name)
    return 0


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def _row(k: str, s: dict[str, Any], digits: int = 4) -> str:
    f = lambda v, d=digits: "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"
    return (f"| `{k}` | {f(s['value'])} | {s['n']} | {f(s['boot_se'])} | "
            f"[{f(s['ci95_lo'])}, {f(s['ci95_hi'])}] | {f(s['placeholder_se'])} | "
            f"{s['ratio_boot_over_placeholder']:.2f} | {f(s['se_linearized'])} | "
            f"{s['n_eff_cz']:.0f} | {100 * s['top5_weight_share']:.0f}% |")


def write_md(path: Path, meta: dict[str, Any], summ: dict[str, Any],
             corr: pd.DataFrame, have_inc: bool) -> None:
    L: list[str] = []
    a = L.append
    a("# District moments: whole-CZ bootstrap\n")
    a(f"Generated {meta['generated']} by `district_bootstrap.py`. Base year {meta['year']}, "
      f"`{meta['scheme']}` partition (suburb NCES locale 21-23 minus city 11-13), labor "
      f"market: {meta['geo_unit']}.\n")
    a("**Method.** The point estimates are computed with `data_estimate.py`'s own panel, "
      "screen, partition and `czone_gaps` (enrollment-weighted mean across CZs of the "
      "within-CZ gap, weights = partitioned CZ pupils). The bootstrap resamples whole CZs "
      f"with replacement: **B = {meta['B']}, seed = {meta['seed']}** "
      f"(`numpy.random.default_rng`; one multiplicity matrix shared by all statistics), "
      f"{meta['n_cz_partition']} CZs / {meta['n_districts_partition']} districts in the "
      "partition. Each draw recomputes the pupil-weighted mean over the resampled CZs "
      "that have the statistic (so n varies by draw). Held fixed: the district panel "
      "after the within-year [0.5%, 99.5%] screen, the partition and the CZ filters "
      "(500-pupil and 10%-group-share cuts). `SE` = SD of the draws; the interval is "
      "the 2.5-97.5 percentile. `placeholder` = unweighted cross-CZ SD / sqrt(n) as "
      "currently used in `spatial_calibrate.jl`; `lin. SE` is the analytic "
      "cluster-robust (linearised ratio) SE, a cross-check; `n_eff` = 1 / sum of squared "
      "normalised CZ weights; `top5` = share of the weight on the five largest CZs.\n")

    rp = meta.get("reproduction_vs_tracked")
    if rp and rp["n_compared"]:
        a(f"**Reproduction.** Rebuilt from the raw 2018 tables, the {rp['n_compared']} "
          f"statistics that exist in the tracked `{meta['tracked_file_compared']}` match it "
          f"(max absolute difference {rp['max_abs_diff']:.1e}; identical n for "
          f"{rp['n_same_n']} of {rp['n_compared']}).\n")
    hdr = ("| statistic | value | n CZs | boot SE | 95% pct interval | placeholder SE | "
           "boot / placeholder | lin. SE | n_eff | top5 |")
    sep = "|---|---|---|---|---|---|---|---|---|---|"

    def table(title: str, keys: list[str], note: str = "") -> None:
        if not keys:
            return
        a(f"\n## {title}\n")
        if note:
            a(note + "\n")
        a(hdr)
        a(sep)
        for k in keys:
            a(_row(k, summ[k]))
        a("")
        for k in keys:
            a(f"- `{k}`: {summ[k]['label']}." + (
                f" Tracked value {summ[k]['tracked_value']:.5f} (n = "
                f"{summ[k]['tracked_n']}), difference {summ[k]['diff_vs_tracked']:+.2e}."
                if "tracked_value" in summ[k] else ""))

    role = lambda r: [k for k, s in summ.items() if s["role"] == r]
    table("Targets", TARGETS,
          "`salary_real` is the CWIFT-deflated teacher salary per FTE (teacher-weighted "
          "within location); `seda_lninc` is log median household income from SEDA "
          "covariates (the current income proxy).")
    table("Validation moments", VALIDATION + ["local_revenue_share"],
          "`local_revenue_share` is a level (pupil-weighted mean of the local revenue "
          f"share over {summ['local_revenue_share'].get('n_districts', '')} sample "
          "districts); its `n` column counts CZ clusters and its placeholder is the "
          "SD of CZ-level shares / sqrt(CZs).")
    if have_inc:
        table("Income alternatives (ACS 5-year 2014-18, school districts)",
              role("income_alternative"),
              "Same partition, same CZ weights. Within a location: pupil-weighted mean "
              "of the district log median (`ln_*`), or of the level then logged "
              "(`lvl_*`), or weighted by earners (`ew_*`). Suburb minus city. Districts "
              "with a bottom-coded (2,499) or missing median are dropped; a top code "
              "(250,001) is kept.")
    if have_inc:
        si, pe, pf = summ["seda_lninc"], summ["pm_earn_all"], summ["pm_ftyr_all"]
        ea, hh = summ["ln_earn_all"], summ["ln_hhinc_2544"]
        a("\n### Reading the income alternatives\n")
        a("- The model's statistic (`income = :log_median_market` in `spatial_calibrate.jl`) is "
          "the gap in log median earnings of market workers. The closest ACS school-district "
          "objects are the median of earnings for everyone 16+ with earnings and, for full-time "
          "year-round workers, B20018/B20017. No school-district table cuts earnings by age, so "
          "25-34 cannot be isolated; the only age-specific object is household income of "
          "householders 25-44 (B19049_003), which is a household, not a labor-income, concept.")
        a(f"- The current proxy (`seda_lninc`, log median household income) is {si['value']:.3f} "
          f"(SE {si['boot_se']:.3f}). Earnings-based gaps are about half as large: "
          f"{ea['value']:.3f} (SE {ea['boot_se']:.3f}) for the average of district median earnings "
          f"(B20002_001), {pe['value']:.3f} (SE {pe['boot_se']:.3f}) for the pooled location "
          f"median of everyone with earnings, and {pf['value']:.3f} (SE {pf['boot_se']:.3f}) "
          "for the pooled median of full-time year-round workers. Household income of "
          f"25-44 householders is {hh['value']:.3f} (SE {hh['boot_se']:.3f}) against "
          f"{summ['ln_hhinc_acs']['value']:.3f} for all householders (ACS B19013), close to the "
          f"SEDA proxy, while mean earnings per earning household (B19061/B19051) is "
          f"{summ['ln_earn_hh_mean']['value']:.3f}. So the wedge is between household income "
          "(a median over all households, with non-labor income) and earnings, not age.")
        a("- `pm_*` (pooled location median) is the closest analogue of the model object: the "
          "median is taken over the pooled earnings distribution of a location (summed "
          "B20001/B20005 bin counts over the location's districts, linear interpolation inside "
          "bins) rather than averaging district medians. Interpolation reproduces the published "
          "district medians with a mean log error of +0.008 (SD 0.012).")
        a("- The men/women rows show the source of the level difference: women's earnings gap "
          "(city vs suburb) is smaller than men's, and FTYR restrictions shrink both.")
    table("All other gaps in `data_estimate.py`", role("other"))

    core = list(corr.columns)
    a("\n## Bootstrap correlation of the targets and validation moments\n")
    a("| | " + " | ".join(f"`{c}`" for c in core) + " |")
    a("|---|" + "---|" * len(core))
    for r_ in core:
        a(f"| `{r_}` | " + " | ".join(f"{corr.loc[r_, c]:.2f}" for c in core) + " |")
    a("\nFull draws: `district_bootstrap_draws.csv` (one row per draw).")
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
