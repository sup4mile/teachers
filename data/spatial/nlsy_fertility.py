#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_fertility.py -- probability of ever having a (biological) child, from the
NLSY79 (born 1957-64), for the "child shock" of the spatial teachers model.

In the model each adult (both genders) has, with exogenous probability p, one
child at the transition to adulthood.  p is calibrated to the share of adults who
EVER have a child (completed fertility): the NLSY79 NUMKIDyy item, "number of
children ever born", read at each respondent's last valid interview, restricted to
respondents last observed at age >= 45.  Same cohort as the project's city/suburb
move-rate target (T3) and mother-child ability correlation.

    cd data/spatial && python3 nlsy_fertility.py              # B = 500 bootstrap
    python3 nlsy_fertility.py --extract                       # re-cut the CSV in the zip
    python3 nlsy_fertility.py --quick                         # B = 50
    python3 nlsy_fertility.py --extract --nrows 200 --scratch /tmp/x   # debug path

Extract step: variables are selected from the NLS `.sdf` index by question name
(`fetch_data.read_nlsy_sdf` / `select_nlsy_vars`) and cut from the full-cohort CSV
inside raw/_downloads/nlsy79_all_1979-2022.zip into raw/nlsy79_fertility.parquet
(+ metadata/nlsy79_fertility_vars.csv).  Reused when present unless --extract.
Estimate step: weighted shares (last-round sampling weight / 100), household-
cluster bootstrap (HHID clusters, weights carried).  Writes
estimates/nlsy_fertility.json and .md (no microdata).

Written with help from Claude Code.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fetch_data as fd  # noqa: E402
import nlsy_ability as na  # noqa: E402
import nlsy_extensions as ne  # noqa: E402

LOG = logging.getLogger("nlsy_fertility")
ZIP = HERE / "raw" / "_downloads" / "nlsy79_all_1979-2022.zip"
SEED = 20261018
COMPLETE_AGE = 45

BLOCKS = {
    "id": r"CASEID|HHID|SAMPLE_ID|SAMPLE_SEX",
    "birth": r"Q1-3_A~Y",
    "weight": r"SAMPWEIGHT",
    "age": r"AGEATINT",
    "kids": r"NUMKID(\d\d)?",
    "geo": r"SMSARES",
}

# Published cross-checks (looked up 2026-10-18; see the .md for URLs).
CROSSCHECK = {
    "nsfg_2015_19_women_40_49_ever_birth": 84.3,
    "nsfg_2015_19_men_40_49_ever_fathered": 76.5,
    "nsfg_2011_15_men_40_44_ever_fathered": 80.4,
    "cps_2018_women_40_44_childless_pct": 15,
    "cps_2016_women_40_44_mothers_pct": 86,
}


# =============================================================================
# Extract
# =============================================================================

def extract(raw: Path, meta: Path, nrows: int | None) -> pd.DataFrame:
    LOG.info("extract: reading %s", ZIP)
    with zipfile.ZipFile(ZIP) as zf:
        members = {Path(m).suffix.lower(): m for m in zf.namelist()}
        index = fd.read_nlsy_sdf(zf.read(members[".sdf"]).decode("latin-1"))
        sel = fd.select_nlsy_vars(index, BLOCKS)
        LOG.info("  %d variables selected of %d", len(sel), len(index))
        t0 = time.time()
        with zf.open(members[".csv"]) as fh:
            df = pd.read_csv(fh, usecols=sel["rnum"].tolist(), nrows=nrows)
        LOG.info("  read %d rows in %.0fs", len(df), time.time() - t0)
    df = df[sel["rnum"].tolist()]
    raw.mkdir(parents=True, exist_ok=True)
    meta.mkdir(parents=True, exist_ok=True)
    df.to_parquet(raw / "nlsy79_fertility.parquet", index=False)
    sel.to_csv(meta / "nlsy79_fertility_vars.csv", index=False)
    return df


class Panel:
    """Raw table + dictionary; `get` replaces NLS negative codes by NaN."""

    def __init__(self, raw: Path, meta: Path) -> None:
        self.df = pd.read_parquet(raw / "nlsy79_fertility.parquet")
        self.meta = pd.read_csv(meta / "nlsy79_fertility_vars.csv", dtype=str)
        self.map = {(r.qname, str(r.year)): r.rnum for r in self.meta.itertuples()}

    def get(self, q: str, y) -> np.ndarray:
        x = self.df[self.map[(q, str(y))]].to_numpy(dtype=float)
        return np.where(x < 0, np.nan, x)

    def has(self, q: str, y) -> bool:
        return (q, str(y)) in self.map


# =============================================================================
# Respondent-level construction
# =============================================================================

def build(P: Panel) -> dict:
    n = len(P.df)
    sample = P.get("SAMPLE_ID", 1979)
    pop = np.isin(np.nan_to_num(sample, nan=-1).astype(int), na.POP_SAMPLES)
    female = P.get("SAMPLE_SEX", 1979) == 2
    hh = P.get("HHID", 1979)
    w79 = P.get("SAMPWEIGHT", 1979) / 100.0
    age79 = P.get("AGEATINT", 1979)
    byear = P.get("Q1-3_A~Y", 1979)

    # rounds with a NUMKIDyy item
    rounds = sorted(int(y) for (q, y) in P.map if q.startswith("NUMKID") and y.isdigit())
    K = np.full((n, len(rounds)), np.nan)
    A = np.full((n, len(rounds)), np.nan)
    W = np.full((n, len(rounds)), np.nan)
    for j, Y in enumerate(rounds):
        yy = f"{Y % 100:02d}"
        K[:, j] = P.get(f"NUMKID{yy}", Y)
        A[:, j] = P.get("AGEATINT", Y)
        W[:, j] = P.get("SAMPWEIGHT", Y) / 100.0
    kx = P.get("NUMKID", "XRND") if P.has("NUMKID", "XRND") else np.full(n, np.nan)

    valid = np.isfinite(K) & np.isfinite(A)            # valid kids count AND age
    has_any = valid.any(axis=1)
    last = np.where(has_any, valid.shape[1] - 1 - np.argmax(valid[:, ::-1], axis=1), 0)
    r = np.arange(n)
    last_age = np.where(has_any, A[r, last], np.nan)
    last_w = np.where(has_any, W[r, last], np.nan)
    last_year = np.where(has_any, np.array(rounds)[last], 0)
    last_k = np.where(has_any, K[r, last], np.nan)

    allk = np.column_stack([K, kx])
    maxk = np.where(np.isfinite(allk).any(axis=1), np.nanmax(np.where(np.isfinite(allk), allk, -1), axis=1), np.nan)

    # status at the LAST interview at ages 40-44 (cumulative max up to that round)
    in4044 = valid & (A >= 40) & (A < 45)
    h4 = in4044.any(axis=1)
    l4 = np.where(h4, in4044.shape[1] - 1 - np.argmax(in4044[:, ::-1], axis=1), 0)
    cum = np.fmax.accumulate(np.where(np.isfinite(K), K, -1), axis=1)
    k4 = np.where(h4, cum[r, l4], np.nan)
    w4 = np.where(h4, W[r, l4], np.nan)

    # origin sample of the T3 move-rate target (build_trans79)
    smsa = P.get("SMSARES", 1979)
    t = ne.smsa_type79(smsa, 1979)                      # 0 central city, 1 suburb, 2 not in SMSA
    orig = (pop & (age79 >= ne.ORIGIN_AGE_79[0]) & (age79 <= ne.ORIGIN_AGE_79[1])
            & np.isfinite(w79) & (w79 > 0) & np.isfinite(smsa) & np.isin(t, (0, 1)))

    # diagnostics on NUMKID80 (info-sheet item) and monotonicity
    k80 = K[:, rounds.index(1980)] if 1980 in rounds else np.full(n, np.nan)
    later = np.nanmax(np.where(np.isfinite(K[:, 1:]), K[:, 1:], -1), axis=1) if K.shape[1] > 1 else np.full(n, -1)
    diag = {
        "n_numkid80_valid": int(np.isfinite(k80).sum()),
        "n_numkid80_pos": int((k80 >= 1).sum()),
        "n_numkid80_pos_but_later_all_zero": int(((k80 >= 1) & (later == 0)).sum()),
        "n_last_below_max": int((has_any & (last_k < maxk)).sum()),
        "n_xrnd_valid": int(np.isfinite(kx).sum()),
        "n_xrnd_differs_from_max_rounds": int((np.isfinite(kx) & has_any & (kx != np.nanmax(np.where(np.isfinite(K), K, -1), axis=1))).sum()),
        "rounds": rounds,
        "birth_year_min_max": [float(np.nanmin(byear[pop])), float(np.nanmax(byear[pop]))],
    }
    return dict(n=n, pop=pop, female=female, hh=hh, w79=w79, has_any=has_any, last_age=last_age,
                last_w=last_w, last_year=last_year, maxk=maxk, h4=h4, k4=k4, w4=w4, orig=orig, diag=diag)


# =============================================================================
# Specs: (mask, weight, outcome) -> share; bootstrap via cluster sums
# =============================================================================

def specs(D: dict) -> dict:
    pop, f = D["pop"], D["female"]
    ever = (D["maxk"] >= 1).astype(float)
    ok_last = pop & D["has_any"] & np.isfinite(D["last_w"]) & (D["last_w"] > 0)
    out = {}

    def add(name, mask, w, y, grp=("pooled", "women", "men")):
        for g in grp:
            gm = {"pooled": np.ones_like(f), "women": f, "men": ~f}[g].astype(bool)
            m = mask & gm & np.isfinite(w) & np.isfinite(y)
            out[(name, g)] = (m, np.where(m, w, 0.0), np.where(m, y, 0.0))

    c45 = lambda a: ok_last & (D["last_age"] >= a)
    add("main", c45(COMPLETE_AGE), D["last_w"], ever)
    add("a_w79", pop & D["has_any"] & (D["last_age"] >= COMPLETE_AGE) & np.isfinite(D["w79"]) & (D["w79"] > 0),
        D["w79"], ever)
    add("b_age40", c45(40), D["last_w"], ever)
    add("c_unweighted", pop & D["has_any"] & (D["last_age"] >= COMPLETE_AGE), np.ones(D["n"]), ever)
    add("d_t3_origin", c45(COMPLETE_AGE) & D["orig"], D["last_w"], ever)
    ok4 = pop & D["h4"] & np.isfinite(D["w4"]) & (D["w4"] > 0)
    add("e_women_40_44", ok4, D["w4"], (D["k4"] >= 1).astype(float), grp=("women",))
    add("f_mean_kids", c45(COMPLETE_AGE), D["last_w"], D["maxk"])
    # share of parents among mothers/fathers is derived: mean kids | any child
    add("f_mean_kids_parents", c45(COMPLETE_AGE) & (D["maxk"] >= 1), D["last_w"], D["maxk"])
    return out


def point(sp: dict) -> dict:
    res = {}
    for k, (m, w, y) in sp.items():
        den = w.sum()
        res[k] = {"est": float((w * y).sum() / den) if den > 0 else np.nan,
                  "n": int(m.sum()), "sum_w": float(den)}
    return res


def bootstrap(sp: dict, hh: np.ndarray, pop: np.ndarray, B: int) -> dict:
    """Resample HHID clusters (of all POP_SAMPLES respondents) with replacement."""
    ids = np.where(np.isfinite(hh), hh, -np.arange(len(hh)) - 1)    # missing HHID: own cluster
    cl = pd.factorize(ids[pop])[0]
    nc = cl.max() + 1
    rng = np.random.default_rng(SEED)
    keys = list(sp)
    num = np.zeros((nc, len(keys)))
    den = np.zeros((nc, len(keys)))
    for j, k in enumerate(keys):
        _, w, y = sp[k]
        num[:, j] = np.bincount(cl, weights=(w * y)[pop], minlength=nc)
        den[:, j] = np.bincount(cl, weights=w[pop], minlength=nc)
    draws = np.empty((B, len(keys)))
    for b in range(B):
        c = np.bincount(rng.integers(0, nc, nc), minlength=nc).astype(float)
        draws[b] = (c @ num) / np.where(c @ den > 0, c @ den, np.nan)
    return {k: float(np.nanstd(draws[:, j], ddof=1)) for j, k in enumerate(keys)}


# =============================================================================
# Report
# =============================================================================

def pct(x): return f"{100 * x:.1f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", action="store_true", help="re-cut the CSV inside the zip")
    ap.add_argument("--B", type=int, default=500)
    ap.add_argument("--quick", action="store_true", help="B = 50")
    ap.add_argument("--nrows", type=int, default=None, help="debug: read only the first N CSV rows")
    ap.add_argument("--scratch", type=Path, default=None,
                    help="debug: write raw/metadata/estimates under this directory instead of data/spatial")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    B = 50 if a.quick else a.B
    root = a.scratch or HERE
    raw, meta, est = root / "raw", root / "metadata", root / "estimates"
    est.mkdir(parents=True, exist_ok=True)

    if a.extract or a.nrows or not (raw / "nlsy79_fertility.parquet").exists():
        extract(raw, meta, a.nrows)
    P = Panel(raw, meta)
    LOG.info("loaded %d respondents x %d columns", *P.df.shape)
    D = build(P)
    sp = specs(D)
    res = point(sp)
    LOG.info("bootstrap B=%d", B)
    se = bootstrap(sp, D["hh"], D["pop"], B)

    # ---- coverage ----
    pop, f = D["pop"], D["female"]
    def cov(mask):
        w = D["w79"]
        return {"n": int(mask.sum()), "n_women": int((mask & f).sum()), "n_men": int((mask & ~f).sum())}
    la = D["last_age"]
    coverage = {
        "POP_SAMPLES respondents": cov(pop),
        "no valid NUMKID observation": cov(pop & ~D["has_any"]),
        "last valid obs at age < 40": cov(pop & D["has_any"] & (la < 40)),
        "last valid obs at age 40-44": cov(pop & D["has_any"] & (la >= 40) & (la < 45)),
        "last valid obs at age >= 45 (main sample before weight screen)": cov(pop & D["has_any"] & (la >= 45)),
        "main sample (age >= 45, positive last-round weight)": cov(sp[("main", "pooled")][0]),
    }
    w79 = D["w79"]
    okw = pop & np.isfinite(w79) & (w79 > 0)
    cov_share = float(w79[okw & D["has_any"] & (la >= 45)].sum() / w79[okw].sum())
    # share with kids among those NOT in the main sample but observed earlier (descriptive)
    notcomp = pop & D["has_any"] & (la < 45)
    ever = (D["maxk"] >= 1)
    short = float(np.average(ever[notcomp & okw], weights=w79[notcomp & okw])) if (notcomp & okw).any() else np.nan

    # ---- json ----
    def row(name):
        return {g: {"est": res[(name, g)]["est"], "se": se.get((name, g)), "n": res[(name, g)]["n"]}
                for g in ("pooled", "women", "men") if (name, g) in res}
    names = ["main", "a_w79", "b_age40", "c_unweighted", "d_t3_origin", "e_women_40_44",
             "f_mean_kids", "f_mean_kids_parents"]
    out = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "definition": "share with >=1 biological child ever born (max over NUMKIDyy rounds and NUMKID XRND), "
                      "last valid obs at age >= 45, last-round SAMPWEIGHT/100, POP_SAMPLES, both genders pooled",
        "B": B, "seed": SEED, "nrows_debug": a.nrows,
        "main_p": res[("main", "pooled")]["est"], "main_p_se": se[("main", "pooled")],
        "variants": {nm: row(nm) for nm in names},
        "coverage": coverage,
        "coverage_weighted_share_w79_with_complete_obs": cov_share,
        "share_ever_child_among_observed_but_younger_than_45_at_last_obs_w79": short,
        "diagnostics": D["diag"],
        "crosscheck": CROSSCHECK,
    }
    (est / "nlsy_fertility.json").write_text(json.dumps(out, indent=2, default=float))

    # ---- markdown ----
    def cell(nm, g, share=True):
        if (nm, g) not in res: return "--"
        r = res[(nm, g)]
        e = pct(r["est"]) if share else f"{r['est']:.2f}"
        s = se.get((nm, g))
        s = (f" ({100 * s:.1f})" if share else f" ({s:.2f})") if s is not None and np.isfinite(s) else ""
        return f"{e}{s} [n={r['n']}]"
    lab = {"main": "MAIN: last obs age >= 45, last-round weight",
           "a_w79": "(a) 1979 base-year weight",
           "b_age40": "(b) age cutoff >= 40",
           "c_unweighted": "(c) unweighted",
           "d_t3_origin": "(d) T3 origin sample (MSA youths 14-16 in 1979)",
           "e_women_40_44": "(e) women, status at an interview at ages 40-44"}
    m = res[("main", "pooled")]
    md = []
    md.append("# NLSY79: probability of ever having a child\n")
    md.append(f"Generated by `data/spatial/nlsy_fertility.py` (B = {B} household-cluster bootstrap draws, seed {SEED}). "
              "No microdata.\n")
    md.append("## Definition and sample\n")
    md.append("- Outcome: ever had a biological child = max over rounds of NUMKIDyy (\"number of children ever born\", "
              "1980/1982-2022) and the cross-round NUMKID is >= 1. NUMKID80 is the 1980 info-sheet item; it is included in the max "
              f"(diagnostic: {D['diag']['n_numkid80_pos']} respondents report >= 1 child in 1980, "
              f"{D['diag']['n_numkid80_pos_but_later_all_zero']} of them report 0 in every later round).")
    md.append("- Sample: `nlsy_ability.POP_SAMPLES` respondents (cross-sectional + retained supplemental samples; military and "
              f"dropped supplements excluded), born 19{int(D['diag']['birth_year_min_max'][0])}-19{int(D['diag']['birth_year_min_max'][1])} (Q1-3_A~Y, two-digit years); "
              "completed-fertility sample = last valid (NUMKID, AGEATINT) observation at age >= 45; both genders pooled; "
              "weight = SAMPWEIGHT/100 of that last round (zero/missing dropped).")
    md.append("- Mapping to the model: p = share of adults (both genders) who ever have a child; model child arrives at the "
              "transition to adulthood with probability p.\n")
    md.append("## Main estimate\n")
    md.append(f"**p = {pct(m['est'])}% (SE {100 * se[('main', 'pooled')]:.1f}), n = {m['n']}.** "
              f"Women {cell('main', 'women')}, men {cell('main', 'men')} (percent; SE in parentheses).\n")
    md.append("## Variants (percent with >= 1 child; bootstrap SE in parentheses)\n")
    md.append("| Variant | Pooled | Women | Men |\n|---|---|---|---|")
    for nm, l in lab.items():
        md.append(f"| {l} | {cell(nm, 'pooled')} | {cell(nm, 'women')} | {cell(nm, 'men')} |")
    md.append("")
    md.append("(f) Descriptive, main sample (last-round weight): mean number of children ever born "
              f"(max over rounds) pooled {cell('f_mean_kids', 'pooled', False)}, women {cell('f_mean_kids', 'women', False)}, "
              f"men {cell('f_mean_kids', 'men', False)}; mean among those with any child: pooled "
              f"{cell('f_mean_kids_parents', 'pooled', False)}, women {cell('f_mean_kids_parents', 'women', False)}, "
              f"men {cell('f_mean_kids_parents', 'men', False)}. (SEs for the means are bootstrap SEs too.)\n")
    md.append("## Coverage and attrition\n")
    md.append("| Group | n | women | men |\n|---|---|---|---|")
    for k, v in coverage.items():
        md.append(f"| {k} | {v['n']} | {v['n_women']} | {v['n_men']} |")
    md.append("")
    d_att = res[("main", "pooled")]["est"] - res[("a_w79", "pooled")]["est"]
    md.append(f"Attrition check: main (last-round weights) {pct(m['est'])}% vs (a) 1979 base-year weights "
              f"{pct(res[('a_w79', 'pooled')]['est'])}%: difference {100 * d_att:+.1f} pp. "
              f"Respondents with a complete (age >= 45) observation hold {100 * cov_share:.1f}% of the 1979 weight among "
              f"POP_SAMPLES respondents with a positive 1979 weight. Among those observed only before age 45, "
              f"the 1979-weighted share with a child (as last seen) is {pct(short)}% (a descriptive contrast only: these respondents are younger when last seen, so their status is not completed fertility).\n")
    md.append("## Male under-reporting of fathered children\n")
    mw, mm = res[("main", "women")]["est"], res[("main", "men")]["est"]
    md.append(f"Women {pct(mw)}% vs men {pct(mm)}% ever had a child (gap {100 * (mw - mm):+.1f} pp). Surveys of men are known to "
              "miss fathered children (births outside marriage, non-resident fathers, children of earlier relationships); "
              "the NSFG shows the same pattern (men 76.5% vs women 84.3% at 40-49). The model's p applies to both genders, "
              "so the women's rate is an upper bound for adults and the men's rate a lower bound; the pooled rate "
              "is the default.\n")
    md.append("## Cross-check against published figures\n")
    md.append("Looked up by web search on 2026-10-18 (secondary sources as noted; the Census Bureau's own June 2018 table was not "
              "retrievable, so the CPS figure is from published analyses of it).\n")
    md.append(f"- Women 40-44, CPS June Fertility Supplement: 15% childless in 2018 (NCFMR/BGSU Family Profile FP-20-04, "
              "Guzzo & Schweizer), i.e. about 85% ever had a child; Pew (2018) reports 86% of women 40-44 are mothers in 2016 "
              "(80% in 2006). NLSY79 variant (e): "
              f"{pct(res[('e_women_40_44', 'women')]['est'])}% (SE {100 * se[('e_women_40_44', 'women')]:.1f}) -- compare directly. "
              "Pew's 86% is from search-result text; the page itself could not be fetched.")
    md.append("- NSFG 2015-2019, ages 40-49 (NCHS NHSR No. 179): 84.3% of women have given birth and 76.5% of men "
              f"have fathered a child (men 80.4% at ages 40-44 in NSFG 2011-2015, NHSR No. 113). NLSY79 main: women {pct(mw)}%, men {pct(mm)}%.")
    md.append("- Period note: the NLSY79 cohort (born 1957-64) was 40-44 in 1997-2008; NCFMR reports CPS childlessness at 40-44 "
              "of 19% in 2000 (10% in 1980), so a completed-fertility rate near 81% mothers is the cohort-relevant CPS benchmark; "
              "the 2016-18 values (85-86%) are for later cohorts.")
    md.append("- Sources: https://www.bgsu.edu/ncfmr/resources/data/family-profiles/guzzo-schweizer-number-children-women-40-44-1980-2018-fp-20-04.html ; "
              "https://www.pewresearch.org/social-trends/2018/01/18/theyre-waiting-longer-but-u-s-women-today-more-likely-to-have-children-than-a-decade-ago/ ; "
              "https://www.cdc.gov/nchs/data/nhsr/nhsr179.pdf ; https://www.cdc.gov/nchs/data/nhsr/nhsr113.pdf\n")
    md.append("## Diagnostics\n")
    md.append(f"- Respondents whose last valid NUMKID is below their maximum over rounds (reporting inconsistency): {D['diag']['n_last_below_max']}.")
    md.append(f"- Cross-round NUMKID valid for {D['diag']['n_xrnd_valid']}; differs from max of round items for {D['diag']['n_xrnd_differs_from_max_rounds']}.")
    md.append(f"- Rounds used: {D['diag']['rounds']}.\n")
    md.append("*Written with help from Claude Code.*\n")
    (est / "nlsy_fertility.md").write_text("\n".join(md))
    LOG.info("main p = %.4f (SE %.4f), n=%d", m["est"], se[("main", "pooled")], m["n"])
    LOG.info("wrote %s", est / "nlsy_fertility.md")


if __name__ == "__main__":
    main()
