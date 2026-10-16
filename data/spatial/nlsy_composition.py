#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_composition.py -- NLSY79: residential composition (central city vs suburb of an
SMSA) by the respondent's OWN ability and OWN income.

`nlsy_extensions.py` (T3) tabulates city/suburb transitions by PARENTAL-income
tercile.  A model in which people sort across locations by their own type needs the
same objects by OWN type.  This script uses the T3 machinery (`build_trans79`,
`Trans`, `stats_from`, `tercile_of`, `collapse_rows`; same samples, weights,
definitions, household-cluster bootstrap with B = 500 and the seed of the NLSY79
main table, seed + 2 = 20260931) to produce

  1. the published parental-income tercile table, rebuilt as a check: city to
     suburb %, suburb to city %, move rate %, origin city share %;
  2. the same table by tercile of the respondent's own standardized AFQT
     (`afqt_c` of `nlsy_ability.py`: equal-weight mean of the four ASVAB IRT
     z-scores, each standardized, re-standardized; 1979 weights, block-A sample;
     terciles cut over the origin-age population with a valid score, as the
     income terciles are cut over the origin-age population with valid income);
  3. the adult (25-34) suburb share, among person-rounds in an SMSA with known
     central-city status, by tercile of own AFQT (cut over the whole block-A
     sample) and by tercile of own adult earnings (mean over the person's 25-34
     rounds of wage-and-salary income, CPI-U deflated to 2010 dollars, zeros
     included; terciles over all respondents with at least one earnings
     observation, weighted by their mean round weight).  The sample is every
     NLSY79 population-sample respondent observed at 25-34 (no origin
     requirement); the 2x2-sample version (respondents with a known 1979 origin
     at 14-16) is a sensitivity, as are earnings terciles over positive earners,
     over full-year full-time wage workers (the ACS-style sample of
     `nlsy_ability.py`) and within sex.

Top-minus-bottom differences carry bootstrap SEs.  Every draw resamples NLSY79
households (HHID) with replacement, the same multiplicities for every table
(common random numbers), re-standardizes AFQT with the draw weights and recomputes
every tercile cut.  Terciles are cut with `tercile_of` (weighted quantiles).

    cd data/spatial && python3 nlsy_composition.py              # B = 500
    python3 nlsy_composition.py --quick                         # B = 50

Reads raw/nlsy79_ability.parquet (via `nlsy_ability.Cohort`) and
estimates/nlsy_extensions.json (published table, for the check).  Writes
estimates/settings/nlsy_composition.json and .md (no microdata).

Written with help from Claude Code.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402
import nlsy_extensions as ne  # noqa: E402

LOG = logging.getLogger("nlsy_composition")

SEED = 20260929
SEED_OFFSET = 2           # nlsy_extensions.t3 bootstraps the NLSY79 main table with seed + 2
B_DEFAULT = 500
TRANS_KEYS = ("city_to_suburb", "suburb_to_city", "move_rate", "origin_share_city")


# =============================================================================
# Data: persons, AFQT, earnings, adult rows
# =============================================================================

def afqt_composite(D: na.N79, w: np.ndarray) -> np.ndarray:
    """`afqt_c` of nlsy_ability.pop_afqt without the factor fit (which the composite does
    not use): each of the four IRT z-scores standardized in the block-A sample with
    weights `w`, their mean re-standardized."""
    m = D.asamp
    Z, wm = D.z[m], w[m]
    Zs = np.column_stack([na.standardize(Z[:, k], wm) for k in range(Z.shape[1])])
    out = np.full(len(w), np.nan)
    out[m] = na.standardize(Zs.mean(axis=1), wm)
    return out


def adult_rows(c: na.Cohort, orig: np.ndarray | None) -> pd.DataFrame:
    """Person-rounds at 25-34 in an SMSA (or not) with a known classification, in the
    layout of `build_trans79`: (p, a, w, rw) with a = 0 city, 1 suburb, 2 non-SMSA and
    rw = round weight / the person's number of rows.  `orig` (if given) restricts to
    persons whose 1979 origin at 14-16 is known (orig >= 0): then the rows are exactly
    those of `build_trans79`."""
    pop = np.isin(c.get("SAMPLE_ID", 1979).astype(int), na.POP_SAMPLES)
    P, A, W = [], [], []
    for Y in sorted(int(y) for y in c.by_q["SMSARES"] if str(y).isdigit()):
        if not (c.has("AGEATINT", Y) and c.has("SAMPWEIGHT", Y)):
            continue
        age = c.get("AGEATINT", Y)
        m = c.get("SMSARES", Y)
        w = c.get("SAMPWEIGHT", Y) / 100.0
        t = ne.smsa_type79(m, Y)
        ok = ((age >= ne.ADULT_AGE[0]) & (age <= ne.ADULT_AGE[1]) & np.isfinite(m)
              & np.isfinite(w) & (w > 0) & (t >= 0) & pop)
        if orig is not None:
            ok &= orig >= 0
        i = np.flatnonzero(ok)
        P.append(i)
        A.append(t[i])
        W.append(w[i])
    p, a, w = (np.concatenate(v) for v in (P, A, W))
    return ne.collapse_rows(p, a, w)


def earnings_person(c: na.Cohort, D: na.N79, pop: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Own adult earnings, per person: mean over the 25-34 person-rounds of real (2010
    dollars, CPI-U of the reference year) wage-and-salary income, and the person's weight
    (mean of his or her round weights over those rounds, the `collapse_rows` convention).

    `all`   every round with an income report (zeros included)
    `pos`   rounds with positive income only
    `fyft`  rounds in the ACS-style full-year full-time wage sample (`PY79.fy`)
    Returns {name: (earnings[n], weight[n])}, NaN / 0 for persons with no such round."""
    py = D.py
    n = D.n
    wR = np.full(len(py.p), np.nan)
    for Y in np.unique(py.Y):
        sel = py.Y == Y
        wR[sel] = c.get("SAMPWEIGHT", int(Y))[py.p[sel]] / 100.0
    real = py.inc * ne.CPI2010 / np.array([ne.CPI[int(y) - 1] for y in py.Y])
    base = np.isfinite(real) & np.isfinite(wR) & (wR > 0) & pop[py.p]
    out = {}
    for name, sel in (("all", base), ("pos", base & (real > 0)), ("fyft", base & py.fy)):
        g = pd.DataFrame({"p": py.p[sel], "e": real[sel], "w": wR[sel]}).groupby("p")
        e, ew = np.full(n, np.nan), np.zeros(n)
        mean_e, mean_w = g.e.mean(), g.w.mean()
        e[mean_e.index.to_numpy()] = mean_e.to_numpy()
        ew[mean_w.index.to_numpy()] = mean_w.to_numpy()
        out[name] = (e, ew)
    return out


class Frame:
    """Everything that does not change across bootstrap draws."""

    def __init__(self, c: na.Cohort) -> None:
        self.c = c
        D = self.D = na.build_n79(c)
        self.n = D.n
        self.pop = np.isin(c.get("SAMPLE_ID", 1979).astype(int), na.POP_SAMPLES)
        self.female = D.female.astype(bool)
        age0 = c.get("AGEATINT", 1979)
        self.w0 = c.get("SAMPWEIGHT", 1979) / 100.0
        self.origin_age = ((age0 >= ne.ORIGIN_AGE_79[0]) & (age0 <= ne.ORIGIN_AGE_79[1])
                           & self.pop & np.isfinite(self.w0) & (self.w0 > 0))
        self.T = ne.build_trans79(c)                       # 2x2 rows, parental income
        rows_all = adult_rows(c, None)
        rows_chk = adult_rows(c, self.T.orig)
        a = rows_chk.sort_values(["p", "a", "w"]).reset_index(drop=True)
        b = pd.DataFrame({"p": self.T.p, "a": self.T.a, "w": self.T.rows.w.to_numpy(),
                          "rw": self.T.rw}).sort_values(["p", "a", "w"]).reset_index(drop=True)
        self.rows_match_t3 = bool(np.array_equal(a.p, b.p) and np.array_equal(a.a, b.a)
                                  and np.allclose(a.rw, b.rw))
        self.rows = {"all": rows_all, "2x2": rows_chk}
        self.earn = earnings_person(c, D, self.pop)
        self.cl_idx = self.T.cl_idx
        self.ncl = self.T.ncl


# =============================================================================
# Tercile assignment (recomputed in every draw) and statistics
# =============================================================================

def groupings(fr: Frame, mult: np.ndarray) -> dict[str, tuple[np.ndarray, Any]]:
    """Tercile code per person (0, 1, 2; -1 = no value) and cut points, for every
    grouping variable, under person multiplicities `mult`."""
    D, T = fr.D, fr.T
    G: dict[str, tuple[np.ndarray, Any]] = {}
    G["income"] = ne.tercile_of(T.income, T.origw * mult * T.inc_pop)
    afqt = afqt_composite(D, D.w * mult)
    G["afqt_origin"] = ne.tercile_of(afqt, fr.w0 * mult * (fr.origin_age & np.isfinite(afqt)))
    G["afqt_full"] = ne.tercile_of(afqt, D.w * mult * D.asamp)
    for k, (e, ew) in fr.earn.items():
        G[f"earn_{k}"] = ne.tercile_of(e, ew * mult)
    e, ew = fr.earn["all"]
    t = np.full(fr.n, -1)
    cuts = []
    for sx in (False, True):
        s = fr.female == sx
        ts, cs = ne.tercile_of(np.where(s, e, np.nan), np.where(s, ew * mult, 0.0))
        t = np.where(s, ts, t)
        cuts.append(cs)
    G["earn_all_sex"] = (t, tuple(cuts))
    return G


def trans_stats(T: ne.Trans, terc: np.ndarray, mult: np.ndarray) -> dict[str, float]:
    """Origin x adult counts by tercile (as `Trans.matrices`) -> the four table statistics
    by tercile and their top-minus-bottom differences."""
    mp = mult[T.p]
    o = T.orig[T.p]
    wt = T.rw * mp
    t = terc[T.p]
    okt = (o >= 0) & (t >= 0)
    M = np.zeros((3, 3, 3))
    np.add.at(M, (t[okt], o[okt], T.a[okt]), wt[okt])
    out: dict[str, float] = {}
    st = [ne.stats_from(M[k]) for k in range(3)]
    for key in TRANS_KEYS:
        for k in range(3):
            out[f"{key}|{k + 1}"] = st[k][key]
        out[f"{key}|diff"] = st[2][key] - st[0][key]
    return out


def adult_stats(rows: dict[str, np.ndarray], terc: np.ndarray, mult: np.ndarray
                ) -> dict[str, float]:
    """Suburb share of the person-rounds with a known city/suburb status, by tercile of the
    person, top minus bottom, and overall."""
    p, a, rw = rows["p"], rows["a"], rows["rw"]
    wt = rw * mult[p]
    in2 = a <= 1
    out: dict[str, float] = {}
    with np.errstate(invalid="ignore", divide="ignore"):
        out["overall"] = float((wt[in2] * (a[in2] == 1)).sum() / wt[in2].sum())
        t = terc[p]
        ok = in2 & (t >= 0)
        num = np.bincount(t[ok], weights=wt[ok] * (a[ok] == 1), minlength=3)
        den = np.bincount(t[ok], weights=wt[ok], minlength=3)
        sh = num / den
    for k in range(3):
        out[f"share|{k + 1}"] = float(sh[k])
    out["share|diff"] = float(sh[2] - sh[0])
    out["share|diff_21"] = float(sh[1] - sh[0])
    out["share|diff_32"] = float(sh[2] - sh[1])
    return out


def rows_arrays(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {"p": df.p.to_numpy(), "a": df.a.to_numpy(), "rw": df.rw.to_numpy()}


# (name, grouping, rows) of the adult-share tables
ADULT_TABLES = [
    ("afqt", "afqt_full", "all"),
    ("earn_all", "earn_all", "all"),
    ("earn_pos", "earn_pos", "all"),
    ("earn_fyft", "earn_fyft", "all"),
    ("earn_all_sex", "earn_all_sex", "all"),
    ("afqt_2x2sample", "afqt_full", "2x2"),
    ("earn_all_2x2sample", "earn_all", "2x2"),
]
TRANS_TABLES = [("income", "income"), ("afqt", "afqt_origin")]


def one_draw(fr: Frame, mult: np.ndarray, rows: dict[str, dict[str, np.ndarray]]
             ) -> tuple[dict[str, float], dict[str, Any]]:
    G = groupings(fr, mult)
    out: dict[str, float] = {}
    for name, g in TRANS_TABLES:
        for k, v in trans_stats(fr.T, G[g][0], mult).items():
            out[f"trans|{name}|{k}"] = v
        cuts = G[g][1]
        out[f"trans|{name}|cut|1"], out[f"trans|{name}|cut|2"] = float(cuts[0]), float(cuts[1])
    for name, g, rs in ADULT_TABLES:
        for k, v in adult_stats(rows[rs], G[g][0], mult).items():
            out[f"adult|{name}|{k}"] = v
    for g in ("afqt_full", "earn_all", "earn_pos", "earn_fyft"):
        out[f"cut|{g}|1"], out[f"cut|{g}|2"] = float(G[g][1][0]), float(G[g][1][1])
    return out, G


# =============================================================================
# Bootstrap
# =============================================================================

def run(fr: Frame, B: int, seed: int) -> dict[str, Any]:
    rows = {k: rows_arrays(v) for k, v in fr.rows.items()}
    point, G0 = one_draw(fr, np.ones(fr.n), rows)
    rng = np.random.default_rng(seed)
    draws: list[dict[str, float]] = []
    t0 = time.time()
    for b in range(B):
        cm = np.bincount(rng.integers(0, fr.ncl, fr.ncl), minlength=fr.ncl).astype(float)
        draws.append(one_draw(fr, cm[fr.cl_idx], rows)[0])     # same stream as t3_estimate
        if (b + 1) % 100 == 0 or b + 1 == B:
            LOG.info("  bootstrap %d/%d (%.1fs)", b + 1, B, time.time() - t0)
    keys = list(point)
    arr = np.array([[d[k] for k in keys] for d in draws])
    se = np.nanstd(arr, axis=0, ddof=1)
    lo, hi = np.nanpercentile(arr, [2.5, 97.5], axis=0)
    res = {k: {"value": float(point[k]), "se": float(se[i]), "ci95_lo": float(lo[i]),
               "ci95_hi": float(hi[i])} for i, k in enumerate(keys)}
    counts = sample_counts(fr, G0, rows)
    return {"stats": res, "counts": counts, "B": B, "seed": seed}


def sample_counts(fr: Frame, G0: dict[str, tuple[np.ndarray, Any]],
                  rows: dict[str, dict[str, np.ndarray]]) -> dict[str, Any]:
    """Person-years and persons per tercile (point-estimate classification)."""
    T = fr.T
    out: dict[str, Any] = {"trans": {}, "adult": {}}
    for name, g in TRANS_TABLES:
        t = G0[g][0][T.p]
        o = T.orig[T.p]
        d = {}
        for k in range(3):
            sel = (t == k) & (o >= 0) & (o <= 1) & (T.a <= 1)
            d[str(k + 1)] = {"person_years": int(sel.sum()), "persons": int(np.unique(T.p[sel]).size)}
        out["trans"][name] = d
    for name, g, rs in ADULT_TABLES:
        r = rows[rs]
        t = G0[g][0][r["p"]]
        d = {}
        for k in range(3):
            sel = (t == k) & (r["a"] <= 1)
            d[str(k + 1)] = {"person_years": int(sel.sum()),
                             "persons": int(np.unique(r["p"][sel]).size)}
        sel = r["a"] <= 1
        d["all"] = {"person_years": int(sel.sum()), "persons": int(np.unique(r["p"][sel]).size)}
        out["adult"][name] = d
    return out


# =============================================================================
# Checks
# =============================================================================

def reproduction(fr: Frame, res: dict[str, Any], outdir: Path, indir: Path, seed: int,
                 B: int) -> dict[str, Any]:
    """Published parental-income table vs (i) `nlsy_extensions.t3_estimate` re-run with the
    published seed and (ii) this script's generic loop."""
    pub_path = indir / "estimates" / "nlsy_extensions.json"
    chk: dict[str, Any] = {"published_file": str(pub_path.relative_to(indir)) if pub_path.exists() else None}
    rr = ne.t3_estimate(fr.T, B, seed)                       # the code that produced the table
    chk["rerun_of_t3_estimate"] = {t: {"stats": rr["tercile"][t], "se": rr["tercile_se"][t]}
                                   for t in ("1", "2", "3")}
    chk["rerun_cuts"] = rr["tercile_cuts"]
    diffs = []
    for k, t in enumerate(("1", "2", "3")):
        for key in TRANS_KEYS:
            mine = res["stats"][f"trans|income|{key}|{k + 1}"]
            diffs.append(abs(mine["value"] - rr["tercile"][t][key]))
            diffs.append(abs(mine["se"] - rr["tercile_se"][t][key]))
    chk["max_abs_diff_generic_vs_t3_estimate"] = float(max(diffs))
    if pub_path.exists():
        pub = json.loads(pub_path.read_text())["T3"]["nlsy79_main"]
        d2 = []
        for t in ("1", "2", "3"):
            for key in TRANS_KEYS:
                d2.append(abs(pub["tercile"][t][key] - rr["tercile"][t][key]))
                d2.append(abs(pub["tercile_se"][t][key] - rr["tercile_se"][t][key]))
        chk["max_abs_diff_rerun_vs_published"] = float(max(d2))
        chk["published"] = {"tercile": pub["tercile"], "tercile_se": pub["tercile_se"],
                            "tercile_cuts": pub["tercile_cuts"], "tercile_diff": pub["tercile_diff"],
                            "B": json.loads(pub_path.read_text())["meta"]["B"],
                            "seed_used_for_nlsy79_main": json.loads(pub_path.read_text())["meta"]["seed"] + SEED_OFFSET}
    return chk


def afqt_checks(fr: Frame) -> dict[str, Any]:
    D = fr.D
    ref, *_ = na.pop_afqt(D.z, D.asamp, D.w)
    mine = afqt_composite(D, D.w)
    both = np.isfinite(ref) & np.isfinite(mine)
    e, ew = fr.earn["all"]
    m = np.isfinite(mine) & np.isfinite(e)
    tf = ne.tercile_of(mine, D.w * D.asamp)[0]
    to = ne.tercile_of(mine, fr.w0 * (fr.origin_age & np.isfinite(mine)))[0]
    s = fr.origin_age & np.isfinite(mine)
    return {"afqt_composite_equals_pop_afqt_max_abs_diff": float(np.max(np.abs(ref[both] - mine[both]))),
            "n_afqt_sample": int(D.asamp.sum()),
            "n_origin_age_with_afqt": int(s.sum()),
            "tercile_agreement_origin_vs_full_cut": float(np.mean(to[s] == tf[s])),
            "weighted_corr_afqt_earnings_all": float(na.wcorr(mine[m], e[m], ew[m])),
            "n_corr": int(m.sum()),
            "n_persons_with_earnings": {k: int(np.isfinite(v[0]).sum()) for k, v in fr.earn.items()},
            "adult_rows_equal_t3_rows_when_origin_required": fr.rows_match_t3}


# =============================================================================
# Report
# =============================================================================

def pse(r: dict[str, float], scale: float = 100.0, d: int = 1) -> str:
    return f"{scale * r['value']:.{d}f} ({scale * r['se']:.{d}f})"


def write_md(path: Path, meta: dict[str, Any], res: dict[str, Any], chk: dict[str, Any],
             ac: dict[str, Any]) -> None:
    S, C = res["stats"], res["counts"]
    L: list[str] = []
    a = L.append
    a("# NLSY79: residential composition by own ability and own income\n")
    a(f"Generated {meta['generated']} by `nlsy_composition.py`. NLSY79 (`raw/nlsy79_ability`), the "
      "transition machinery of `nlsy_extensions.py` (T3): same samples, weights and definitions. "
      f"Household-cluster bootstrap, B = {meta['B']}, seed {meta['seed']} (= 20260929 + 2, the seed "
      "`nlsy_extensions.t3` uses for the NLSY79 main table); every table uses the same household "
      "multiplicities in every draw, re-standardizes AFQT with the draw weights and recomputes all "
      "tercile cut points. SEs are bootstrap SDs. *City* = central city of an SMSA, *suburb* = in "
      "an SMSA, not in the central city; person-rounds with 'central city not known' are "
      "excluded.\n")

    def trans_table(name: str, title: str) -> None:
        rows = []
        for k in ("1", "2", "3"):
            rows.append([f"T{k}", pse(S[f"trans|{name}|city_to_suburb|{k}"]),
                         pse(S[f"trans|{name}|suburb_to_city|{k}"]),
                         pse(S[f"trans|{name}|move_rate|{k}"]),
                         pse(S[f"trans|{name}|origin_share_city|{k}"]),
                         f"{C['trans'][name][k]['person_years']} / {C['trans'][name][k]['persons']}"])
        a(title)
        a(ne.tab(rows, ["tercile", "city to suburb %", "suburb to city %", "move rate %",
                        "origin city share %", "person-years / persons"]))
        d = [S[f"trans|{name}|{key}|diff"] for key in TRANS_KEYS]
        a("Top minus bottom tercile (pp): " + "; ".join(
            f"{lab} {pse(r)}" for lab, r in zip(("city to suburb", "suburb to city", "move rate",
                                                 "origin city share"), d)) + ".\n")

    a("## Summary\n")
    dA = S["adult|afqt|share|diff"]
    dE = S["adult|earn_all|share|diff"]
    a(f"- **Adult suburb share (25-34, SMSA with known status), overall {pse(S['adult|afqt|overall'])}%.** "
      f"By own AFQT tercile: {pse(S['adult|afqt|share|1'])}, {pse(S['adult|afqt|share|2'])}, "
      f"{pse(S['adult|afqt|share|3'])}; top minus bottom **{pse(dA)} pp** "
      f"(middle minus bottom {pse(S['adult|afqt|share|diff_21'])}, top minus middle "
      f"{pse(S['adult|afqt|share|diff_32'])}). By own-earnings tercile: "
      f"{pse(S['adult|earn_all|share|1'])}, {pse(S['adult|earn_all|share|2'])}, "
      f"{pse(S['adult|earn_all|share|3'])}; top minus bottom **{pse(dE)} pp**.")
    dm, dc, ds = (S["trans|afqt|move_rate|diff"], S["trans|afqt|city_to_suburb|diff"],
                  S["trans|afqt|suburb_to_city|diff"])
    a(f"- **Transitions by own AFQT tercile** (origin 14-16, adult 25-34): top minus bottom, move "
      f"rate {pse(dm)} pp, city to suburb {pse(dc)} pp, suburb to city {pse(ds)} pp; origin "
      f"city share {pse(S['trans|afqt|origin_share_city|diff'])} pp. The same four differences by "
      f"parental-income tercile (the published table; origin city share from the same draws): "
      f"{pse(S['trans|income|move_rate|diff'])}, {pse(S['trans|income|city_to_suburb|diff'])}, "
      f"{pse(S['trans|income|suburb_to_city|diff'])} and "
      f"{pse(S['trans|income|origin_share_city|diff'])} pp.")
    d_pub = chk.get("max_abs_diff_rerun_vs_published")
    a("- Reproduction of the published parental-income table (statistics and SEs): maximum absolute "
      f"difference {'n/a' if d_pub is None else f'{d_pub:.1e}'} (rerun of `t3_estimate` with the "
      f"published seed vs the saved table), {chk['max_abs_diff_generic_vs_t3_estimate']:.1e} "
      "(this script's generic code vs the rerun).\n")

    a("## 1. Reproduction: transitions by parental-income tercile (NLSY79 main table)\n")
    cuts = [S["trans|income|cut|1"]["value"], S["trans|income|cut|2"]["value"]]
    a(f"Terciles of 1978 family income (`TNFI_TRUNC`) over respondents aged 14-16 in 1979 with valid "
      f"income; cut points {cuts[0]:,.0f} and {cuts[1]:,.0f} (nominal dollars). Rebuilt here / "
      "published (`estimates/nlsy_extensions.md`, section B):\n")
    pub = chk.get("published")
    rows = []
    for k in ("1", "2", "3"):
        mine = [pse(S[f"trans|income|{key}|{k}"]) for key in TRANS_KEYS[:3]]
        mine.append(f"{100 * S[f'trans|income|origin_share_city|{k}']['value']:.1f}")
        if pub:
            ps = pub["tercile"][k]
            pe = pub["tercile_se"][k]
            pb = [f"{100 * ps[key]:.1f} ({100 * pe[key]:.1f})" for key in TRANS_KEYS[:3]]
            pb.append(f"{100 * ps['origin_share_city']:.1f}")
        else:
            pb = [""] * 4
        rows.append([f"T{k}"] + [f"{m} / {p}" for m, p in zip(mine, pb)] +
                    [f"{C['trans']['income'][k]['person_years']} / {C['trans']['income'][k]['persons']}"])
    a(ne.tab(rows, ["tercile", "city to suburb %", "suburb to city %", "move rate %",
                    "origin city share %", "person-years / persons"]))
    if pub:
        d = pub["tercile_diff"]
        a("Top minus bottom, published: move rate {:.1f} ({:.1f}) pp, city to suburb {:.1f} ({:.1f}) "
          "pp, suburb to city {:.1f} ({:.1f}) pp; rebuilt: move rate {}, city to suburb {}, suburb "
          "to city {}.\n".format(
              100 * d["move_rate_T3_minus_T1"][0], 100 * d["move_rate_T3_minus_T1"][1],
              100 * d["city_to_suburb_T3_minus_T1"][0], 100 * d["city_to_suburb_T3_minus_T1"][1],
              100 * d["suburb_to_city_T3_minus_T1"][0], 100 * d["suburb_to_city_T3_minus_T1"][1],
              pse(S["trans|income|move_rate|diff"]), pse(S["trans|income|city_to_suburb|diff"]),
              pse(S["trans|income|suburb_to_city|diff"])))

    a("## 2. Transitions by tercile of own AFQT\n")
    c1, c2 = S["trans|afqt|cut|1"]["value"], S["trans|afqt|cut|2"]["value"]
    a(f"Own AFQT = `afqt_c` of `nlsy_ability.py` (equal-weight mean of the four standardized ASVAB IRT "
      f"z-scores, re-standardized; 1979 weights; block-A sample, N = {ac['n_afqt_sample']}); terciles "
      f"cut over respondents aged 14-16 in 1979 with a valid score (N = "
      f"{ac['n_origin_age_with_afqt']}), cut points {c1:.3f} and {c2:.3f} SD. Same origin "
      "(1979 residence at 14-16), adult rounds (25-34) and weights as the table above.\n")
    trans_table("afqt", "")

    a("## 3. Adult (25-34) suburb share by own AFQT and by own adult earnings\n")
    a("Share of person-rounds in the suburb among person-rounds in an SMSA with known central-city "
      "status; every NLSY79 population-sample respondent observed at 25-34 (no origin "
      "requirement); weights = round `SAMPWEIGHT` / the person's number of classified rounds at 25-34 "
      "(city, suburb or non-SMSA; as in T3). "
      "AFQT terciles: cut over the whole block-A sample (cut points "
      f"{S['cut|afqt_full|1']['value']:.3f} and {S['cut|afqt_full|2']['value']:.3f} SD; the origin-age cut "
      f"classifies {100 * ac['tercile_agreement_origin_vs_full_cut']:.1f}% of 14-16-year-olds identically). "
      "Earnings: mean over the person's 25-34 rounds of wage-and-salary income (Q13-5; CPI-U of the "
      "reference year, 2010 dollars; zeros included), terciles over respondents with at least one "
      f"earnings observation (N = {ac['n_persons_with_earnings']['all']}), weighted by their mean round "
      f"weight; cut points ${S['cut|earn_all|1']['value']:,.0f} and ${S['cut|earn_all|2']['value']:,.0f}. "
      f"The weighted correlation of own AFQT and own mean earnings is {ac['weighted_corr_afqt_earnings_all']:.2f} "
      f"(N = {ac['n_corr']}).\n")

    def adult_row(name: str, label: str) -> list[str]:
        r = [label]
        for k in ("1", "2", "3"):
            r.append(pse(S[f"adult|{name}|share|{k}"]))
        r.append(pse(S[f"adult|{name}|share|diff"]))
        r.append(pse(S[f"adult|{name}|share|diff_21"]))
        r.append(pse(S[f"adult|{name}|share|diff_32"]))
        r.append(pse(S[f"adult|{name}|overall"]))
        cc = C["adult"][name]
        r.append(" / ".join(f"{cc[k]['person_years']}-{cc[k]['persons']}" for k in ("1", "2", "3")))
        return r

    head = ["grouping", "T1 %", "T2 %", "T3 %", "T3 - T1 (pp)", "T2 - T1 (pp)", "T3 - T2 (pp)",
            "overall %", "person-years-persons by tercile"]
    a("Main tables:\n")
    a(ne.tab([adult_row("afqt", "own AFQT"), adult_row("earn_all", "own adult earnings")], head))
    a("Sensitivities:\n")
    sens = [("earn_pos", "earnings, positive earners (positive rounds only)"),
            ("earn_fyft", "earnings, full-year full-time wage workers (ACS-style rounds)"),
            ("earn_all_sex", "earnings, terciles within sex"),
            ("afqt_2x2sample", "own AFQT, 2x2 sample (origin at 14-16 known)"),
            ("earn_all_2x2sample", "own earnings, 2x2 sample")]
    a(ne.tab([adult_row(n, lab) for n, lab in sens], head))
    a("Earnings tercile cut points (2010 dollars): "
      f"all rounds {S['cut|earn_all|1']['value']:,.0f} / {S['cut|earn_all|2']['value']:,.0f}; "
      f"positive {S['cut|earn_pos|1']['value']:,.0f} / {S['cut|earn_pos|2']['value']:,.0f}; "
      f"FYFT {S['cut|earn_fyft|1']['value']:,.0f} / {S['cut|earn_fyft|2']['value']:,.0f}.\n")

    a("## Notes and caveats\n")
    a("- **Same machinery, different conditioning.** Sections 1-2 condition on a known SMSA status at 14-16 "
      "(in 1979 29% of 14-16-year-olds are outside any SMSA and 21% are in an SMSA whose central-city "
      "status is not known) and on being classified at 25-34; section 3 does not require an origin.")
    a("- **Origin at 14-16 is the 1979 interview residence.** Cross-SMSA movers cannot be separated "
      "(see the T3 caveats in `nlsy_extensions.md`); AFQT was measured in 1980 when respondents were "
      "15-23, after the 1979 residence, and its terciles are age-normed scores.")
    a("- **'Central city not known'** is 33-41% of SMSA residents at 25-34 through 1996 and is dropped; "
      "if unknown status is related to ability or income the shares are selected.")
    a("- **Own earnings** are wage-and-salary income reports (top-coded in the revised variable) at 25-34 "
      "for men and women: women's zeros enter the bottom tercile unless the within-sex or "
      "positive-earner variants are used. Terciles use the person's mean over rounds, which averages "
      "out transitory income but mixes survey years.")
    a("- **SEs** cluster on the NLSY79 household (siblings share a draw); AFQT standardization and all cut "
      "points are recomputed in every draw. Differences are computed draw by draw, so their SEs "
      "include the covariance between terciles.")
    path.write_text("\n".join(L) + "\n")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bootstrap", type=int, default=B_DEFAULT)
    p.add_argument("--seed", type=int, default=SEED,
                   help="base seed; the bootstrap uses seed + 2 as nlsy_extensions.t3 does for NLSY79")
    p.add_argument("--quick", action="store_true", help="B = 50")
    p.add_argument("--indir", default=str(HERE))
    p.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    B = 50 if args.quick else args.bootstrap
    seed = args.seed + SEED_OFFSET
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates" / "settings"
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    c = na.Cohort("nlsy79", indir)
    fr = Frame(c)
    LOG.info("persons %d, households %d; 2x2 person-rounds %d; all adult person-rounds %d",
             fr.n, fr.ncl, len(fr.T.p), len(fr.rows["all"]))
    res = run(fr, B, seed)
    chk = reproduction(fr, res, outdir, indir, seed, B)
    ac = afqt_checks(fr)
    LOG.info("reproduction: rerun vs published max abs diff %s; generic vs rerun %.2e",
             chk.get("max_abs_diff_rerun_vs_published"), chk["max_abs_diff_generic_vs_t3_estimate"])
    meta = {"generated": dt.datetime.now().isoformat(timespec="seconds"),
            "script": "data/spatial/nlsy_composition.py", "B": B, "seed": seed, "base_seed": args.seed,
            "bootstrap": "household (HHID) clusters, multinomial multiplicities, same draws for "
                         "every table; AFQT re-standardized and all terciles recomputed in each draw",
            "tercile_conventions": {
                "income": "NLSY79 TNFI_TRUNC 1978 family income, over respondents aged 14-16 in 1979 "
                          "with valid income (as nlsy_extensions)",
                "afqt_origin": "afqt_c over respondents aged 14-16 in 1979 with a valid score",
                "afqt_full": "afqt_c over the whole block-A sample (1979 weights)",
                "earn_all": "mean over 25-34 rounds of real wage-and-salary income (2010 dollars, "
                            "zeros included); terciles over persons with >= 1 observation, weighted "
                            "by mean round weight",
                "earn_pos": "same, rounds with positive income only",
                "earn_fyft": "same, rounds of the ACS-style full-year full-time wage sample",
                "earn_all_sex": "earn_all terciles cut within sex"},
            "samples": {"adult_rows_all": int(len(fr.rows["all"])),
                        "adult_rows_2x2": int(len(fr.rows["2x2"]))}}
    payload = {"meta": meta, "reproduction": chk, "afqt_and_sample_checks": ac,
               "counts": res["counts"], "statistics": res["stats"]}
    (outdir / "nlsy_composition.json").write_text(json.dumps(payload, indent=1, default=float))
    write_md(outdir / "nlsy_composition.md", meta, res, chk, ac)
    LOG.info("done in %.1fs -> %s", time.time() - t0, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
