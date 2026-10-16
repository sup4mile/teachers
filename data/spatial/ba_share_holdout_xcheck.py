#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ba_share_holdout_xcheck.py -- independent recomputation of the SEALED hold-out B1 (BA-or-more share of adults 25+,
suburb minus city, ACS 2014-18 by school district) written by `ba_share_holdout.py seal`.

The point is a second implementation: nothing here imports `ba_share_holdout.py`, `data_estimate.Gap`,
`district_bootstrap.boot_wmean` or `district_settings.summarise`.  What is shared is only the INPUT panel of the
four district targets (`district_settings.build_audited_panel`, the audited district panel); the partition itself
(city / suburb, >= 500 pupils, each location >= 10% of the CZ) is re-derived below from that panel and compared
with the one the calibration stack builds.  Everything specific to B1 is written again from the raw tables:

    * the ACS records (new pull `raw/acs_b15003_school_districts_2018.parquet` for BA-or-more; the OLD file
      `raw/acs_school_districts.parquet`, `ba_holders` / `pop_25plus`, for the `ba_only` variant);
    * the overlap rules (adults once from unified + elementary records; all layers; unified + secondary; unified
      only), the Detroit refile, New York City, the weights;
    * the ratio-of-totals location shares, the CZ gaps, the CZ-weighted means, the pooled version, the log ratio;
    * the whole-CZ bootstrap (multiplicity matrix from numpy's generator with seed 20260929, B draws; SE = std with
      ddof 1; 95% percentile interval), on the full sample and on the two pre-registered halves (read from
      `estimates/settings/district_settings.json`, not recomputed);
    * the bootstrap correlations with the four district targets, whose draws are recomputed here from the CZ table
      (full sample: and compared with the stored `district_settings_draws_locale.csv`; half A: and compared with
      the stored `district_settings_draws_locale_holdout_A.csv`; half B: recomputed only).

It reads the sealed file, compares every statistic, correlation, bootstrap draw and per-CZ row with the
recomputation, and prints ONLY the block names, PASS/FAIL, counts and the maximum absolute differences.  It never
prints a value, a standard error, a share or a correlation of B1 (or of any variant); a difference is printed, a
value never is.  Run it after `ba_share_holdout.py seal`; it does not write anything.

Usage (from data/spatial)
    python3 ba_share_holdout_xcheck.py [--sealed PATH] [--indir DIR]
On Slurm the pyarrow caveat of `ba_share_holdout.py` applies (ARROW_DEFAULT_MEMORY_POOL=system).

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

YEAR = 2018
MIN_PUPILS = 500
MIN_SHARE = 0.10
NYC = "3620580"
DETROIT_SRC, DETROIT_DST = "2612000", "2601103"
TARGET_COLS = {"gap_pupils": "gap_pupils", "gap_salary": "gap_salary_audited",
               "gap_teachers_pp": "gap_teachers_k12", "gap_income": "gap_pm_ftyr_all"}
REPORTED = {"gap_pupils": (0.3764, 0.1154), "gap_salary": (0.0086, 0.0092),
            "gap_teachers_pp": (-0.0177, 0.0131), "gap_income": (0.1122, 0.0168)}
STATS = ("ba_share_gap", "ba_share_gap_pp", "ba_log_ratio", "ba_share_city", "ba_share_suburb",
         "pooled_gap", "pooled_city", "pooled_suburb")
DRAWN = ("ba_share_gap", "ba_log_ratio", "pooled_gap")
SAMPLES = ("full", "A", "B")

# (key, layers, measure, nyc adults kept, Detroit refiled, NY CZ kept, CZ weights): the ten specifications of the sealed
# file, written out again here from their documented meaning (not read from the file)
SPECS = [
    ("headline", "ue", "plus", True, True, True, "pupils"),
    ("all_layers", "all", "plus", True, True, True, "pupils"),
    ("unified_secondary", "us", "plus", True, True, True, "pupils"),
    ("unified_only", "u", "plus", True, True, True, "pupils"),
    ("ba_only", "ue", "only", True, True, True, "pupils"),
    ("nyc_dropped", "ue", "plus", False, True, True, "pupils"),
    ("drop_ny_cz", "ue", "plus", True, True, False, "pupils"),
    ("no_detroit_crosswalk", "ue", "plus", True, False, True, "pupils"),
    ("adult_weights", "ue", "plus", True, True, True, "adults"),
    ("equal_weights", "ue", "plus", True, True, True, "equal"),
]


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def new_records(path: Path) -> pd.DataFrame:
    """One row per leaid: layer, adults, BA+ (lines 22-25), exactly-BA (line 22); empty/invalid records are NaN."""
    r = pd.read_parquet(path)
    r = r[~r.leaid.astype(str).str.endswith("99999")].drop_duplicates("leaid")
    lines = {j: pd.to_numeric(r[f"B15003_{j:03d}E"], errors="coerce") for j in (1, 22, 23, 24, 25)}
    plus = lines[22] + lines[23] + lines[24] + lines[25]
    good = (lines[1] >= 1) & (plus <= lines[1])
    for j in (22, 23, 24, 25):
        good &= lines[j] >= 0
    out = pd.DataFrame({"leaid": r.leaid.astype(str).to_numpy(),
                        "layer": r.sd_layer.str.replace("school district (", "", regex=False).str.rstrip(")").to_numpy(),
                        "adults": lines[1].where(good).to_numpy(), "plus": plus.where(good).to_numpy(),
                        "only": lines[22].where(good).to_numpy()})
    return out


def old_records(path: Path) -> pd.DataFrame:
    """The existing file (fetch_data.py): `ba_holders` over `pop_25plus`, leaid from state + district code."""
    o = pd.read_parquet(path)
    layer = o.sd_layer.str.replace("school district (", "", regex=False).str.rstrip(")")
    code = np.where(layer == "unified", o["school district (unified)"],
                    np.where(layer == "elementary", o["school district (elementary)"], o["school district (secondary)"]))
    leaid = o.state.astype(str).str.zfill(2) + pd.Series(code, index=o.index).astype(str).str.zfill(5)
    pop = pd.to_numeric(o.pop_25plus, errors="coerce")
    ba = pd.to_numeric(o.ba_holders, errors="coerce")
    good = (pop >= 1) & (ba >= 0) & (ba <= pop)
    out = pd.DataFrame({"leaid": leaid.to_numpy(), "layer": layer.to_numpy(), "adults": pop.where(good).to_numpy(),
                        "only": ba.where(good).to_numpy()})
    out = out[~out.leaid.str.endswith("99999")].drop_duplicates("leaid")
    return out.reset_index(drop=True)


def independent_partition(A: pd.DataFrame) -> pd.DataFrame:
    """City (locale 11-13) / suburb (21-23) districts of >= 500 pupils in a CZ, kept when both locations hold >= 10% of
    the CZ's partitioned pupils: written again from the audited panel."""
    d = A[A.czone.notna() & (A.pupils.fillna(0) >= MIN_PUPILS) & A.locale_grp.isin(["city", "suburb"])].copy()
    if "year" in d:
        d = d[d.year == YEAR]
    d["loc"] = np.where(d.locale_grp.eq("suburb"), 2, 1)
    tot = d.groupby("czone").pupils.transform("sum")
    grp = d.groupby(["czone", "loc"]).pupils.transform("sum")
    d["share"] = grp / tot
    ok = d.groupby("czone").agg(n=("loc", "nunique"), mn=("share", "min"))
    keep = ok.index[(ok.n == 2) & (ok.mn >= MIN_SHARE)]
    return d[d.czone.isin(keep)][["leaid", "czone", "loc", "pupils"]].reset_index(drop=True)


# --------------------------------------------------------------------------
# Recomputation
# --------------------------------------------------------------------------

def multiplicities(n: int, B: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(B, n))
    return np.stack([np.bincount(row, minlength=n) for row in idx]).astype(float)


def wmean_boot(x: np.ndarray, w: np.ndarray, counts: np.ndarray) -> dict[str, Any]:
    use = np.isfinite(x) & np.isfinite(w) & (w > 0)
    num = np.where(use, x * w, 0.0)
    den = np.where(use, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        draws = (counts @ num) / (counts @ den)
    return summary(float(num.sum() / den.sum()), draws, int(use.sum()))


def summary(est: float, draws: np.ndarray, n: int) -> dict[str, Any]:
    f = draws[np.isfinite(draws)]
    lo, hi = np.percentile(f, [2.5, 97.5])
    return dict(value=est, se=float(np.std(f, ddof=1)), ci95_lo=float(lo), ci95_hi=float(hi), n_cz=n, draws=draws)


def spec_cz_table(part: pd.DataFrame, new: pd.DataFrame, old: pd.DataFrame, spec: tuple) -> pd.DataFrame:
    """Per CZ (ascending czone): pupils, and for each location the adults and the BA totals of the used records."""
    key, layers, measure, keep_nyc, detroit, keep_ny_cz, wt = spec
    rec = (old if measure == "only" and key == "ba_only" else new).copy()
    if measure == "only":
        rec["ba"] = rec["only"]
    else:
        rec["ba"] = rec["plus"]
    if detroit:
        rec["leaid"] = rec.leaid.where(rec.leaid != DETROIT_SRC, DETROIT_DST)
    d = part.merge(rec[["leaid", "layer", "adults", "ba"]], on="leaid", how="left")
    valid = d.adults.notna() & d.ba.notna()
    lay = d.layer
    if layers == "ue":
        use = lay.isin(["unified", "elementary"])
    elif layers == "all":
        use = lay.notna()
    elif layers == "u":
        use = lay.eq("unified")
    else:                                     # unified + secondary; elementary only where the cell has no valid secondary record
        sec_cells = set(map(tuple, d.loc[lay.eq("secondary") & valid, ["czone", "loc"]].to_numpy().tolist()))
        in_sec = pd.Series([(c, l) in sec_cells for c, l in zip(d.czone.tolist(), d["loc"].tolist())], index=d.index)
        use = lay.eq("unified") | lay.eq("secondary") | (lay.eq("elementary") & ~in_sec)
    if not keep_nyc:
        use &= d.leaid.ne(NYC)
    d["use"] = valid & use
    d["a"] = d.adults.where(d.use)
    d["b"] = d.ba.where(d.use)
    cz = d.groupby("czone").agg(pupils=("pupils", "sum")).sort_index()
    for loc in (1, 2):
        s = d[d["loc"] == loc].groupby("czone").agg(a=("a", lambda v: v.sum(min_count=1)), b=("b", lambda v: v.sum(min_count=1)))
        cz[f"adults_{loc}"] = s.a.reindex(cz.index)
        cz[f"ba_{loc}"] = s.b.reindex(cz.index)
    for loc in (1, 2):                        # share of the location's pupils whose district has a used record
        dl = d[d["loc"] == loc]
        cz[f"cov_{loc}"] = (dl.pupils.where(dl.use, 0.0).groupby(dl.czone).sum() / dl.pupils.groupby(dl.czone).sum()).reindex(cz.index)
    cz["share_1"], cz["share_2"] = cz.ba_1 / cz.adults_1, cz.ba_2 / cz.adults_2
    cz["gap"] = cz.share_2 - cz.share_1
    ok = (cz.share_1 > 0) & (cz.share_2 > 0)
    cz["logr"] = np.where(ok, np.log(cz.share_2.where(ok) / cz.share_1.where(ok)), np.nan)
    if not keep_ny_cz:
        ny = d.czone[d.leaid == NYC]
        if len(ny):
            cz.loc[ny.iloc[0], ["gap", "logr", "share_1", "share_2"]] = np.nan
    if wt == "pupils":
        cz["w"] = cz.pupils
    elif wt == "adults":
        cz["w"] = cz.adults_1 + cz.adults_2
    else:
        cz["w"] = 1.0
    return cz.reset_index()


def sample_stats(cz: pd.DataFrame, counts: np.ndarray) -> dict[str, dict[str, Any]]:
    """The eight sealed statistics (with draws) for one set of CZs and one multiplicity matrix."""
    w = cz.w.to_numpy(float)
    gap = cz.gap.to_numpy(float)
    have = np.isfinite(gap)
    out = {"ba_share_gap": wmean_boot(gap, w, counts),
           "ba_log_ratio": wmean_boot(cz.logr.to_numpy(float), w, counts),
           "ba_share_city": wmean_boot(np.where(have, cz.share_1.to_numpy(float), np.nan), w, counts),
           "ba_share_suburb": wmean_boot(np.where(have, cz.share_2.to_numpy(float), np.nan), w, counts)}
    g = out["ba_share_gap"]
    out["ba_share_gap_pp"] = dict(value=100 * g["value"], se=100 * g["se"], ci95_lo=100 * g["ci95_lo"],
                                  ci95_hi=100 * g["ci95_hi"], n_cz=g["n_cz"])
    n1, d1 = np.where(have, cz.ba_1, 0.0), np.where(have, cz.adults_1, 0.0)
    n2, d2 = np.where(have, cz.ba_2, 0.0), np.where(have, cz.adults_2, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        dr1, dr2 = (counts @ n1) / (counts @ d1), (counts @ n2) / (counts @ d2)
    e1, e2 = n1.sum() / d1.sum(), n2.sum() / d2.sum()
    n = int(have.sum())
    out["pooled_city"], out["pooled_suburb"] = summary(float(e1), dr1, n), summary(float(e2), dr2, n)
    out["pooled_gap"] = summary(float(e2 - e1), dr2 - dr1, n)
    return out


def target_draws(g: pd.DataFrame, counts: np.ndarray) -> pd.DataFrame:
    """Draws of the four district targets on the CZ rows of `g` (pupil-weighted means of the CZ gaps)."""
    w = g.pupils.to_numpy(float)
    return pd.DataFrame({k: wmean_boot(g[c].to_numpy(float), w, counts)["draws"] for k, c in TARGET_COLS.items()})


# --------------------------------------------------------------------------
# Comparison (prints differences, never values)
# --------------------------------------------------------------------------

class Report:
    def __init__(self) -> None:
        self.rows: list[tuple[str, float, bool]] = []

    def add(self, name: str, diff: float, tol: float) -> None:
        ok = bool(np.isfinite(diff) and diff <= tol)
        self.rows.append((name, float(diff), ok))
        print(f"  {name:<74s} max|diff| {diff:9.2e}   {'PASS' if ok else 'FAIL'}")

    def flag(self, name: str, ok: bool, note: str = "") -> None:
        self.rows.append((name, 0.0 if ok else float("inf"), bool(ok)))
        print(f"  {name:<74s} {'PASS' if ok else 'FAIL'}{('   ' + note) if note else ''}")


def num(x: Any) -> float:
    return float("nan") if x is None else float(x)


def stat_diff(mine: dict[str, Any], sealed: dict[str, Any]) -> float:
    d = [abs(mine[k] - num(sealed[k])) for k in ("value", "se", "ci95_lo", "ci95_hi")]
    d.append(float(abs(mine["n_cz"] - sealed["n_cz"])))
    return float(np.nanmax(d)) if np.isfinite(d).all() else float("inf")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sealed", default=str(HERE / "estimates" / "settings" / "sealed" / "ba_share_holdout_B1.json"))
    ap.add_argument("--indir", default=str(HERE))
    args = ap.parse_args(argv)
    indir = Path(args.indir)
    S = json.loads(Path(args.sealed).read_text())
    B, seed = int(S["meta"]["B"]), int(S["meta"]["bootstrap_seed"])
    rep = Report()
    print(f"xcheck of {Path(args.sealed).name}: B = {B}, bootstrap seed {seed}")

    # ---- inputs shared with the district targets: the audited panel, then the partition re-derived ----
    import district_settings as ds                     # the input panel and the stored split only
    inp = ds.build_audited_panel(indir)
    A = inp["A"]
    part_stack, g = ds.scheme_tables(A, "locale")
    part = independent_partition(A)
    same = (set(part.leaid) == set(part_stack.leaid) and len(part) == len(part_stack)
            and part.set_index("leaid").loc[part_stack.leaid, "loc"].to_numpy().tolist() == part_stack.loc_id.tolist())
    rep.flag("partition re-derived from the audited panel = the calibration stack's", same,
             f"{part.czone.nunique()} CZs, {len(part)} districts")
    new = new_records(indir / "raw" / "acs_b15003_school_districts_2018.parquet")
    old = old_records(indir / "raw" / "acs_school_districts.parquet")
    st = json.loads((indir / "estimates" / "settings" / "district_settings.json").read_text())["holdout"]
    cz_ids = np.sort(part.czone.unique())
    halves = {"full": cz_ids, "A": np.array(sorted(st["A_cz"]), float), "B": np.array(sorted(st["B_cz"]), float)}
    rep.flag("hold-out halves from district_settings.json partition the CZs",
             set(halves["A"]) | set(halves["B"]) == set(cz_ids) and not set(halves["A"]) & set(halves["B"]),
             f"{len(halves['A'])} / {len(halves['B'])}")
    rep.flag("sealed file's CZ lists = the stored halves", sorted(S["samples"]["A"]["cz_ids"]) == sorted(int(c) for c in halves["A"])
             and sorted(S["samples"]["B"]["cz_ids"]) == sorted(int(c) for c in halves["B"])
             and sorted(S["samples"]["full"]["cz_ids"]) == sorted(int(c) for c in cz_ids))
    g = g.set_index("czone").sort_index()
    rep.flag("CZ weights of the targets = partitioned pupils", np.allclose(g.pupils.reindex(cz_ids).to_numpy(float),
                                                                           part.groupby("czone").pupils.sum().reindex(cz_ids).to_numpy(float)))

    # ---- multiplicity matrices and the four public targets ----
    counts = {s: multiplicities(len(halves[s]), B, seed) for s in SAMPLES}
    tdr, tstats = {}, {}
    for s in SAMPLES:
        gs = g.reindex(halves[s]).reset_index()
        tdr[s] = target_draws(gs, counts[s])
    for k, c in TARGET_COLS.items():
        r = wmean_boot(g[c].to_numpy(float), g.pupils.to_numpy(float), counts["full"])
        tstats[k] = r
    full_run = B == 1999 and seed == 20260929
    ok = all(round(tstats[k]["value"], 4) == REPORTED[k][0] and (not full_run or round(tstats[k]["se"], 4) == REPORTED[k][1])
             for k in REPORTED)
    rep.flag("the four public targets reproduced by this code (value" + (" and SE)" if full_run else ")"), ok)
    stored = indir / "estimates" / "settings" / "district_settings_draws_locale.csv"
    if stored.exists() and full_run:
        c = pd.read_csv(stored)
        rep.add("target draws (full sample) vs district_settings_draws_locale.csv",
                max(float(np.nanmax(np.abs(c[k].to_numpy(float) - tdr["full"][k].to_numpy(float)))) for k in TARGET_COLS), 1e-9)
    stored_a = indir / "estimates" / "settings" / "district_settings_draws_locale_holdout_A.csv"
    if stored_a.exists() and full_run:
        c = pd.read_csv(stored_a)
        rep.add("target draws (half A) vs district_settings_draws_locale_holdout_A.csv",
                max(float(np.nanmax(np.abs(c[k].to_numpy(float) - tdr["A"][k].to_numpy(float)))) for k in TARGET_COLS), 1e-9)

    # ---- every specification, every sample ----
    mine: dict[str, dict[str, dict[str, dict[str, Any]]]] = {}
    tables: dict[str, pd.DataFrame] = {}
    for spec in SPECS:
        tab = spec_cz_table(part, new, old, spec)
        tables[spec[0]] = tab
        mine[spec[0]] = {}
        for s in SAMPLES:
            sub = tab[tab.czone.isin(halves[s])].reset_index(drop=True)
            mine[spec[0]][s] = sample_stats(sub, counts[s])

    def sealed_stats(key: str, s: str) -> dict[str, Any]:
        return S["samples"][s]["stats"] if key == "headline" else S["variants"][key]["samples"][s]["stats"]

    d_head = max(stat_diff(mine["headline"][s][k], sealed_stats("headline", s)[k]) for s in SAMPLES for k in STATS)
    rep.add(f"headline: {len(STATS)} statistics x {len(SAMPLES)} samples (value, se, ci95, n_cz)", d_head, 1e-9)
    d_hd = max(stat_diff(mine["headline"][s]["ba_share_gap"], S["headline"][s]) for s in SAMPLES)
    rep.add("headline block: ba_share_gap for full / A / B", d_hd, 1e-9)
    d_pp = max(stat_diff(mine["headline"][s]["ba_share_gap_pp"], S["headline"]["pp"][s]) for s in SAMPLES)
    rep.add("headline block: percentage points", d_pp, 1e-7)
    for spec in SPECS[1:]:
        k = spec[0]
        d = max(stat_diff(mine[k][s][t], sealed_stats(k, s)[t]) for s in SAMPLES for t in STATS)
        rep.add(f"variant {k}: {len(STATS)} statistics x {len(SAMPLES)} samples", d, 1e-9)
    rep.flag("n_cz and n_cz_with_gap of every sample, headline and variants",
             all(S["samples"][s]["n_cz"] == len(halves[s]) and S["samples"][s]["n_cz_with_gap"] == mine["headline"][s]["ba_share_gap"]["n_cz"]
                 for s in SAMPLES)
             and all(S["variants"][sp[0]]["samples"][s]["n_cz_with_gap"] == mine[sp[0]][s]["ba_share_gap"]["n_cz"]
                     for sp in SPECS[1:] for s in SAMPLES))

    # ---- correlations ----
    keys = list(TARGET_COLS) + list(DRAWN)
    d_corr, d_with = 0.0, 0.0
    for s in SAMPLES:
        D = pd.concat([tdr[s].reset_index(drop=True),
                       pd.DataFrame({k: mine["headline"][s][k]["draws"] for k in DRAWN})], axis=1)[keys]
        C = D.corr().to_numpy(float)
        sealed_m = np.array(S["samples"][s]["corr"]["matrix"], dtype=float)
        d_corr = max(d_corr, float(np.nanmax(np.abs(C - sealed_m))))
        for t in TARGET_COLS:
            d_with = max(d_with, abs(float(D[["ba_share_gap", t]].corr().iloc[0, 1]) - S["samples"][s]["corr"]["ba_share_gap_with"][t]))
    rep.add("bootstrap correlation matrices (4 targets + 3 B1 statistics), 3 samples", d_corr, 2e-6)
    rep.add("correlation of ba_share_gap with the four targets, 3 samples", d_with, 1e-9)
    d_var = 0.0
    for spec in SPECS[1:]:
        for s in SAMPLES:
            D = pd.concat([tdr[s].reset_index(drop=True), pd.DataFrame({"x": mine[spec[0]][s]["ba_share_gap"]["draws"]})], axis=1)
            for t in TARGET_COLS:
                d_var = max(d_var, abs(float(D[["x", t]].corr().iloc[0, 1]) - S["variants"][spec[0]]["samples"][s]["corr_with_targets"][t]))
    rep.add("variants: correlation of the gap with the four targets", d_var, 1e-9)

    # ---- draws and per-CZ table ----
    d_dr = max(float(np.nanmax(np.abs(np.asarray(S["draws_full"][k], dtype=float) - mine["headline"]["full"][k]["draws"]))) for k in DRAWN)
    rep.add("stored bootstrap draws of the headline (3 statistics)", d_dr, 1e-9)
    tab = tables["headline"].set_index("czone")
    rows = pd.DataFrame(S["cz_table"]).set_index("czone")
    rows.index = rows.index.astype(float)
    pairs = [("pupils", "pupils"), ("share_city", "share_1"), ("share_suburb", "share_2"), ("gap", "gap"), ("log_ratio", "logr"),
             ("ba_city", "ba_1"), ("ba_suburb", "ba_2"), ("adults_city", "adults_1"), ("adults_suburb", "adults_2"),
             ("pupil_cov_city", "cov_1"), ("pupil_cov_suburb", "cov_2")]
    d_cz = max(float(np.nanmax(np.abs(rows[a].astype(float).reindex(tab.index).to_numpy() - tab[b].to_numpy(float)))) for a, b in pairs)
    nan_same = all((rows[a].astype(float).reindex(tab.index).isna().to_numpy() == tab[b].isna().to_numpy()).all() for a, b in pairs)
    rep.add("per-CZ table of the headline (11 columns, 188 rows)", d_cz, 1e-9)
    rep.flag("per-CZ table: the same CZs have a gap", nan_same)
    half_ok = all((rows.half.reindex(tab.index) == np.where(tab.index.isin(halves["A"]), "A", "B")))
    rep.flag("per-CZ table: half assignment", half_ok)

    # ---- provenance ----
    code = S["meta"].get("code_sha256", {}).get("ba_share_holdout.py")
    rep.flag("the sealed file was written by the ba_share_holdout.py now on disk", code == sha256_file(HERE / "ba_share_holdout.py"))
    inputs = S["meta"].get("inputs_sha256", {})
    rep.flag("the input files are those the sealed file was built from",
             all(sha256_file(indir / "raw" / n) == h for n, h in inputs.items()) and len(inputs) == 3)

    bad = [r for r in rep.rows if not r[2]]
    tol_max = max((r[1] for r in rep.rows if np.isfinite(r[1])), default=0.0)
    print(f"xcheck: {len(rep.rows) - len(bad)} of {len(rep.rows)} blocks PASS; largest difference {tol_max:.2e}"
          + ("" if not bad else "; FAIL: " + ", ".join(r[0] for r in bad)))
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
