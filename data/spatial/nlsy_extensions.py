#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_extensions.py -- two NLSY extensions to the ability block of the spatial
calibration (`spatial_calibration.md`, items T6a and T3).

  T6a  CFR replication in the CNLSY.  Chetty-Friedman-Rockoff (2014b) find about
       12% higher earnings at 28 per student SD of test scores (single-year
       scores in grades 3-8, earnings with zeros as a share of the mean,
       controls).  The NLSY79 slope of log hourly wages on LATENT ability is 0.213
       (ACS-style sample).  This script regresses CNLSY children's earnings at
       ~28 (Young Adult survey) on OBSERVED standardized PIAT/PPVT scores at ages
       8-14 and walks from CFR's definitions to the ACS-style ones one step at a
       time, ending with the disattenuated slope per latent SD.
  T3   Residential-transition target.  Share of children schooled in the central
       city (suburb) of an MSA who live in the suburb (central city) as adults:
       NLSY97 (residence at 12-16 vs 25-34), by parental-income tercile; NLSY79
       (1979 vs 25-34, ages 14-16 in 1979) as a check.

    cd data/spatial && python3 nlsy_extensions.py              # B = 500 bootstrap
    python3 nlsy_extensions.py --quick                         # B = 50
    python3 nlsy_extensions.py --only t3                       # or t6a

Reuses the data access, weighting, norming and factor-model code of
`nlsy_ability.py` (same age-cell rank-normed test scores, same mother weights,
same ACS weeks/hours/income cutoffs and CPI deflator).  Inputs are the tables
`fetch_data.py --sources nlsy` writes to raw/ (this script needs the variables
added to NLSY_BLOCKS for it: geography, family income, YA job history).
Writes estimates/nlsy_extensions.json and .md (no microdata).

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
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402

LOG = logging.getLogger("nlsy_extensions")
CPI = na.CPI_U
CPI2010 = na.CPI_U[2010]

# ---- T6a constants ----------------------------------------------------------
SCORE_AGE_M = (96, 179)         # test ages 8-14 (months)
TARGET_AGE = 28                 # earnings age; interviews at 27-29, closest to 28
TARGET_WINDOW = (27, 29)
ACS_WINDOW = (25, 34)
MIN_TESTS = 2                   # a child-round composite needs >= 2 of the 4 tests
FY_MONTHS = 11                  # >= 11 of 12 reference-year months with a job ~ >= 48 weeks
MIN_HGC = na.MIN_HGC
FT_HOURS = na.FT_HOURS
RACE_LABELS = {1: "Hispanic", 2: "Black", 3: "Other"}

# ---- T3 constants -----------------------------------------------------------
ORIGIN_AGE_97 = (12, 16)
ORIGIN_AGE_79 = (14, 16)
ADULT_AGE = (25, 34)


# =============================================================================
# Small statistics helpers
# =============================================================================

def wls_cluster(y: np.ndarray, X: np.ndarray, w: np.ndarray, cl: np.ndarray,
                j: int = 1) -> tuple[float, float, int]:
    """Weighted OLS; CR1 cluster-robust SE of coefficient j (probability weights).
    Returns (coef, se, number of clusters)."""
    n, k = X.shape
    Xw = X * w[:, None]
    Ai = np.linalg.pinv(Xw.T @ X)
    b = Ai @ (Xw.T @ y)
    u = X * (w * (y - X @ b))[:, None]
    _, inv = np.unique(cl, return_inverse=True)
    G = int(inv.max()) + 1
    Sg = pd.DataFrame(u).groupby(inv).sum().to_numpy()
    V = Ai @ (Sg.T @ Sg) @ Ai * (G / max(G - 1, 1)) * ((n - 1) / max(n - k, 1))
    return float(b[j]), float(np.sqrt(max(V[j, j], 0.0))), G


def dummies(codes: np.ndarray) -> np.ndarray:
    lv = np.unique(codes)[1:]
    return (codes[:, None] == lv[None, :]).astype(float)


def fmt(v: float | None, se: float | None = None, d: int = 3) -> str:
    if v is None or not np.isfinite(v):
        return "n/a"
    return f"{v:.{d}f}" if se is None else f"{v:.{d}f} ({se:.{d}f})"


def get_or_nan(c: na.Cohort, q: str, y: Any, n: int) -> np.ndarray:
    v = c.try_get(q, y)
    return np.full(n, np.nan) if v is None else v


# =============================================================================
# T6a: CFR replication in the CNLSY
# =============================================================================

def ya_panel(ccy: na.Cohort) -> pd.DataFrame:
    """One row per (child, Young Adult round) at ages 25-34.

    Earnings are wage-and-salary income in the calendar year before the interview
    (Q15-5, top-coded from 2006).  Annual weeks and hours are NOT released for the
    CNLSY YA sample, so ACS-style weeks/hours are proxied: months of the reference
    year covered by at least one job (job-history start/stop dates, current jobs
    run to the interview month) and usual weekly hours over all current jobs at the
    interview date (TOTHOURS).  Real income uses CPI-U of the reference year."""
    n = ccy.df.shape[0]
    frames = []
    for yr in range(1994, 2021, 2):
        s, yy = str(yr), f"{yr % 100:02d}"
        an = f"AGEINT{yy}" if ccy.has(f"AGEINT{yy}", s) else f"AGEINT{yr}"
        if not ccy.has(an, s):
            continue
        age = ccy.get(an, s)
        inc = ccy.get("Q15-5", s) if ccy.has("Q15-5", s) else ccy.try_get("Q15-5-TOP", s)
        if inc is None:
            continue
        keep = np.isfinite(age) & (age >= ACS_WINDOW[0]) & (age <= ACS_WINDOW[1])
        if not keep.any():
            continue
        # interview month/year (the round label can differ from the calendar year)
        im = iy = None
        for qm, qy in (("Q-1C_M", "Q-1C_Y"), ("SYMBOL!CURDATE~M", "SYMBOL!CURDATE~Y"),
                       ("CURDATE~M", "CURDATE~Y")):
            if ccy.has(qm, s):
                im, iy = ccy.get(qm, s), ccy.get(qy, s)
                break
        im = np.full(n, np.nan) if im is None else im
        iy = np.full(n, np.nan) if iy is None else iy
        iy = np.where(iy < 100, iy + np.where(iy > 50, 1900, 2000), iy)     # 2-digit years
        iy = np.where(np.isfinite(iy), iy, yr)                              # fallback: round label
        im = np.where(np.isfinite(im), im, 6.0)
        wq = f"YA{yy}WEIGHT_REVISED" if ccy.has(f"YA{yy}WEIGHT_REVISED", s) else f"YA{yy}WEIGHT"
        wt = np.where(get_or_nan(ccy, wq, s, n) > 0, get_or_nan(ccy, wq, s, n), np.nan) / 100.0
        hgc = get_or_nan(ccy, f"HGC{yr}", s, n)
        hgc = np.where(hgc >= 95, np.nan, hgc)             # ungraded
        # HGC{yr} is years of schooling through 2012 and a 1-14 category scale from 2014
        # (2 = some high school): put both on "at least some high school" = 1
        hgc_hs = (hgc >= MIN_HGC) if yr <= 2012 else (hgc >= 2)
        hgc_hs = np.where(np.isfinite(hgc), hgc_hs.astype(float), np.nan)
        tot = get_or_nan(ccy, f"TOTHOURS{yy}", s, n)
        # reference-year months covered by any job
        ref0 = 12 * (iy - 1)
        cov = np.zeros((n, 12), bool)
        unk = np.zeros(n, bool)                 # a job overlaps the year with an unknown stop
        anyjob = np.zeros(n, bool)
        for jj in range(1, 12):
            j = f"{jj:02d}"
            if not ccy.has(f"JOB-HISTORY_START-DATE.{j}~M", s):
                continue
            sm, sy = get_or_nan(ccy, f"JOB-HISTORY_START-DATE.{j}~M", s, n), get_or_nan(ccy, f"JOB-HISTORY_START-DATE.{j}~Y", s, n)
            pm, py = get_or_nan(ccy, f"JOB-HISTORY_STOP-DATE.{j}~M", s, n), get_or_nan(ccy, f"JOB-HISTORY_STOP-DATE.{j}~Y", s, n)
            cur = get_or_nan(ccy, f"JOB-HISTORY_CURRFLAG.{j}", s, n) == 1
            start = 12 * sy + sm
            stop = np.where(cur, 12 * iy + im, 12 * py + pm)
            have_start = np.isfinite(start)
            anyjob |= have_start
            unk |= have_start & ~np.isfinite(stop) & (start <= ref0 + 12)
            ok = have_start & np.isfinite(stop)
            for m in range(12):
                idx = ref0 + m + 1
                cov[:, m] |= ok & (start <= idx) & (stop >= idx)
        months = cov.sum(axis=1).astype(float)
        months[unk & (months < 12)] = np.nan
        idx = np.flatnonzero(keep)
        frames.append(pd.DataFrame({
            "child": idx, "yr": yr, "iy": iy[idx].astype(int), "age": age[idx].astype(int),
            "inc": inc[idx], "hgc_hs": hgc_hs[idx], "tot": tot[idx], "months": months[idx],
            "anyjob": anyjob[idx], "ya_w": wt[idx]}))
    E = pd.concat(frames, ignore_index=True)
    E = E[np.isfinite(E.inc) & E.iy.between(1990, 2023)].reset_index(drop=True)
    E["ref_year"] = E.iy - 1
    E["cpi_ratio"] = CPI2010 / E.ref_year.map(CPI)
    E["inc_real"] = E.inc * E.cpi_ratio
    E["inc_cap"] = np.minimum(E.inc_real, 100000.0)     # CFR cap ($100,000)
    E["floor_ok"] = E.inc >= na.nominal_income_floor(E.iy.to_numpy())
    E["fy"] = E.months >= FY_MONTHS
    E["ft"] = E.tot >= FT_HOURS
    E["hgc_ok"] = ~(E.hgc_hs == 0)                  # missing schooling retained, as in NLSY97
    E["acs"] = E.fy & E.ft & E.floor_ok & E.hgc_ok
    with np.errstate(divide="ignore", invalid="ignore"):
        E["lw"] = np.where(E.acs, np.log(E.inc / (52.0 * E.months / 12.0 * E.tot)), np.nan)
        E["learn"] = np.where(E.inc > 0, np.log(E.inc), np.nan)
    return E


def pick_age28(E: pd.DataFrame) -> pd.DataFrame:
    """Per child, the interview at 27-29 closest to 28 (ties -> the older age:
    the reference-year age is one lower than the interview age)."""
    S = E[E.age.between(*TARGET_WINDOW)].copy()
    S["d"] = (S.age - TARGET_AGE).abs()
    S["o"] = -S.age
    return S.sort_values(["child", "d", "o"]).drop_duplicates("child").drop(columns=["d", "o"])


def child_scores(K: na.Kids) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """(single-year child-round scores, child-level round-average scores, 4-test
    matrix of test-level round means) at ages 8-14.  z = rank-normal within
    3-month age cells (`nlsy_ability.build_kids`); imputed COMP rows excluded."""
    L = K.long
    L = L[(L.age >= SCORE_AGE_M[0]) & (L.age <= SCORE_AGE_M[1]) & ~L.imp & L.z_rank.notna()]
    wide = L.pivot_table(index=["child", "round"], columns="test", values="z_rank")
    wide = wide.reindex(columns=[0, 1, 2, 3])
    age_m = L.groupby(["child", "round"]).age.first()
    Sy = pd.DataFrame({
        "ntests": wide.notna().sum(axis=1),
        "comp": wide.mean(axis=1),
        "math": wide[0],
        "read": wide[[1, 2]].mean(axis=1),
        "age_m": age_m}).reset_index()
    Sy["z0"], Sy["z1"], Sy["z2"], Sy["z3"] = (wide[t].to_numpy() for t in range(4))
    Sy = Sy[Sy.ntests >= MIN_TESTS].reset_index(drop=True)
    Sy["age_y"] = (Sy.age_m // 12).astype(int)
    # previous assessment round (any age 5-14 in the norming window): lagged math,
    # reading and composite, for CFR's cubic-in-prior-scores control
    L0 = K.long[~K.long.imp & K.long.z_rank.notna()]
    W0 = L0.pivot_table(index=["child", "round"], columns="test", values="z_rank").reindex(columns=[0, 1, 2, 3])
    lg = pd.DataFrame({"ntests": W0.notna().sum(axis=1), "lmath": W0[0],
                       "lread": W0[[1, 2]].mean(axis=1), "lcomp": W0.mean(axis=1)}).reset_index()
    lg = lg[lg.ntests >= MIN_TESTS].sort_values(["child", "round"])
    g = lg.groupby("child")
    lg["lag_math"], lg["lag_read"], lg["lag_comp"] = g.lmath.shift(1), g.lread.shift(1), g.lcomp.shift(1)
    lg["lag_gap"] = lg["round"] - g["round"].shift(1)
    Sy = Sy.merge(lg[["child", "round", "lag_math", "lag_read", "lag_comp", "lag_gap"]],
                  on=["child", "round"], how="left")
    tm = L.groupby(["child", "test"]).z_rank.mean().unstack().reindex(columns=[0, 1, 2, 3])
    tm = tm[tm.notna().sum(axis=1) >= MIN_TESTS]
    nr = L.groupby("child")["round"].nunique()
    Sa = pd.DataFrame({"comp": tm.mean(axis=1), "nrounds": nr.reindex(tm.index)}).reset_index()
    T = np.full((K.n, 4), np.nan)
    T[tm.index.to_numpy()] = tm.to_numpy()
    return Sy, Sa, T


def parental_income(c79: na.Cohort, D: na.N79, K: na.Kids) -> np.ndarray:
    """Mother's total net family income (2010 dollars), averaged over survey years
    whose income year falls when the child is 8-14 (NaN if none)."""
    years = sorted(int(y) for y in c79.by_q["TNFI_TRUNC"] if str(y).isdigit())
    inc = np.column_stack([c79.get("TNFI_TRUNC", Y) * CPI2010 / CPI[Y - 1] for Y in years])
    incyr = np.array(years) - 1
    cy = K.cyrb
    ok = (K.mrow >= 0) & np.isfinite(cy)
    out = np.full(K.n, np.nan)
    M = inc[np.maximum(K.mrow, 0)]
    for k in np.flatnonzero(ok):
        sel = (incyr >= cy[k] + 8) & (incyr <= cy[k] + 14)
        v = M[k, sel]
        v = v[np.isfinite(v)]
        if v.size:
            out[k] = v.mean()
    return out


def child_covariates(c79: na.Cohort, ccy: na.Cohort, D: na.N79, K: na.Kids) -> pd.DataFrame:
    hgc_m = c79.get("HGC_EVER", "XRND")
    by = np.where(D.birth_cell >= 0, D.birth_cell // 10, np.nan).astype(float)
    mr = np.maximum(K.mrow, 0)
    valid = K.mrow >= 0
    C = pd.DataFrame({
        "child": np.arange(K.n),
        "mrow": K.mrow,
        "hh": np.where(valid, D.hh_idx[mr], -1),
        "wm": np.where(valid & D.msamp[mr], D.w[mr], np.nan),
        "female": K.female.astype(float),
        "race": ccy.get("CRACE"),
        "m_hgc": np.where(valid, hgc_m[mr], np.nan),
        "m_age_birth": np.where(valid, K.cyrb - by[mr], np.nan),
        "par_inc": parental_income(c79, D, K)})
    C["ln_par_inc"] = np.log(np.maximum(C.par_inc, 1000.0))
    return C


def cubic(x: np.ndarray) -> list[np.ndarray]:
    return [x, x ** 2, x ** 3]


def rung_design(C: pd.DataFrame, obs: pd.DataFrame, score: np.ndarray, weighted: bool,
                ctrl: bool, pooled_age: bool, single: bool, lag: str | None = None
                ) -> np.ndarray:
    """Design matrix (const, standardized score, [controls], age/year dummies,
    [cubics in the previous round's scores])."""
    X = [np.ones(len(obs)), score]
    if ctrl:
        cc = C.loc[obs.child, ["female", "race", "m_hgc", "m_age_birth", "ln_par_inc"]].to_numpy()
        miss = ~np.isfinite(cc[:, 4])
        cc[miss, 4] = 0.0
        X += [cc[:, 0], (cc[:, 1] == 2).astype(float), (cc[:, 1] == 3).astype(float),
              cc[:, 2], cc[:, 3], cc[:, 4], miss.astype(float)]
    else:
        X += [C.loc[obs.child, "female"].to_numpy()]
    dm = [dummies(obs.iy.to_numpy())]
    if pooled_age:
        dm.append(dummies(obs.age.to_numpy()))
    else:
        dm.append(dummies(obs.age.to_numpy()))
    if single:
        dm.append(dummies(obs.age_y.to_numpy()))
    if lag == "subj":            # CFR: cubics in prior-year math and English
        dm.append(np.column_stack(cubic(obs.lag_math.to_numpy()) + cubic(obs.lag_read.to_numpy())))
    elif lag == "comp":
        dm.append(np.column_stack(cubic(obs.lag_comp.to_numpy())))
    return np.column_stack(X + dm)


def std_w(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    return (x - na.wmean(x, w)) / np.sqrt(na.wvar(x, w))


def design_reliability(X: np.ndarray, w: np.ndarray) -> tuple[float, np.ndarray]:
    """One-factor ULS on the 4 test scores; reliability of the mean of the available
    tests as a measure of the common factor (`comp_reliability_avail`)."""
    ok = np.isfinite(w) & (w > 0)
    R = na.wcorr_cross(X[ok], X[ok], w[ok])
    fit = na.fit_factor_uls(R, [0] * 4)
    return na.comp_reliability_avail(fit.lam, R, X[ok], w[ok]), fit.lam


def t6a(c79: na.Cohort, ccy: na.Cohort, D: na.N79, K: na.Kids, out: dict, B: int = 200,
         seed: int = 1) -> None:
    LOG.info("== T6a: CFR replication in the CNLSY ==")
    E = ya_panel(ccy)
    E28 = pick_age28(E)
    Sy, Sa, Tavg = child_scores(K)
    C = child_covariates(c79, ccy, D, K)
    C = C.set_index("child", drop=False)
    elig = C.index[np.isfinite(C.wm) & (C.hh >= 0)]
    LOG.info("  YA person-rounds at 25-34: %d (children %d); at 27-29 closest to 28: %d children",
             len(E), E.child.nunique(), len(E28))
    out["T6a.sample"] = {
        "ya_person_rounds_25_34": int(len(E)), "ya_children_25_34": int(E.child.nunique()),
        "age28_children": int(len(E28)),
        "age28_age_counts": {int(k): int(v) for k, v in E28.age.value_counts().sort_index().items()},
        "age28_interview_years": {int(k): int(v) for k, v in E28.iy.value_counts().sort_index().items()},
        "single_year_score_obs": int(len(Sy)), "children_with_scores_8_14": int(Sa.shape[0]),
        "mean_rounds_per_child": float(Sa.nrounds.mean()),
        "months_covered_dist_age28": {int(k): int(v) for k, v in
                                      E28.months.dropna().astype(int).value_counts().sort_index().items()},
    }

    def prep(Ebase: pd.DataFrame, ycol: str, rows_mask: np.ndarray, pooled: bool):
        """Merge outcome rows with scores and covariates; returns dict of frames."""
        Eb = Ebase[rows_mask & Ebase.child.isin(elig)].copy()
        Eb = Eb[np.isfinite(Eb[ycol])]
        res = {}
        # single-year: one row per (child-round score) x (earnings row)
        M = Eb.merge(Sy, on="child", how="inner")
        M["wm"] = C.loc[M.child, "wm"].to_numpy()
        M["hh"] = C.loc[M.child, "hh"].to_numpy()
        res["single"] = M
        A = Eb.merge(Sa[["child", "comp"]], on="child", how="inner")
        A["wm"] = C.loc[A.child, "wm"].to_numpy()
        A["hh"] = C.loc[A.child, "hh"].to_numpy()
        res["avg"] = A
        return res

    def weights(frame: pd.DataFrame, weighted: bool) -> np.ndarray:
        """Family weight: mother's 1979 weight split across her observations so that
        each CHILD carries an equal share (children of a mother: nc; child's
        earnings rows x score rows: nobs)."""
        if not weighted:
            return np.ones(len(frame))
        key = frame.child.to_numpy()
        nobs = pd.Series(key).map(pd.Series(key).value_counts()).to_numpy()
        mkey = frame.mrow.to_numpy() if "mrow" in frame else C.loc[frame.child, "mrow"].to_numpy()
        nkids = frame.groupby(mkey).child.transform("nunique").to_numpy()
        return frame.wm.to_numpy() / nkids / nobs

    def one_reg(frame: pd.DataFrame, ycol: str, score: str, weighted: bool, ctrl: bool,
                kind: str, transform: str, pooled: bool, lag: str | None = None,
                restrict_lag: bool = False) -> dict:
        m = np.isfinite(frame[score]) & np.isfinite(frame.wm)
        if (lag or restrict_lag) and kind == "single":
            cols = {"comp": ["lag_comp"]}.get(lag, ["lag_math", "lag_read"])
            for cc in cols:
                m &= np.isfinite(frame[cc])
        f = frame[m].copy()
        w = weights(f, weighted)
        y = f[ycol].to_numpy()
        if transform == "mean":
            y = y / na.wmean(y, w)
        s = std_w(f[score].to_numpy(), w)
        X = rung_design(C, f, s, weighted, ctrl, pooled, single=(kind == "single"), lag=lag)
        b, se, G = wls_cluster(y, X, w, f.hh.to_numpy())
        r = {"coef": b, "se": se, "n_obs": int(len(f)), "n_children": int(f.child.nunique()),
             "n_mothers": int(G), "sd_score_raw": float(np.sqrt(na.wvar(f[score].to_numpy(), w)))}
        if kind == "single" and lag:
            r["share_lag_gap_2"] = float((f.lag_gap == 2).mean())
        return r

    SPECS = {"single": [("base", False, None, False), ("ctrl", True, None, False),
                        ("base_lagsample", False, None, True), ("lag", False, "subj", False),
                        ("ctrl_lag", True, "subj", False)],
             "avg": [("base", False, None, False), ("ctrl", True, None, False)]}

    # ---- reliabilities on the CFR-definition sample ---------------------------
    base = prep(E28, "inc_real", np.ones(len(E28), bool), False)
    rel = {}
    for weighted in (True, False):
        fs = base["single"]
        fs = fs[np.isfinite(fs.wm)]
        ws = weights(fs, weighted)
        Rs, lam_s = design_reliability(fs[["z0", "z1", "z2", "z3"]].to_numpy(), ws)
        fa = base["avg"]
        fa = fa[np.isfinite(fa.wm)]
        wa = weights(fa, weighted)
        Ra, lam_a = design_reliability(Tavg[fa.child.to_numpy()], wa)
        rel["weighted" if weighted else "unweighted"] = {
            "single": Rs, "avg": Ra, "lam_single": lam_s.tolist(), "lam_avg": lam_a.tolist()}
        LOG.info("  reliability (%s): single-year %.3f, round-average %.3f",
                 "weighted" if weighted else "unweighted", Rs, Ra)
    out["T6a.reliability"] = rel
    out["T6a.reliability_saved_child_composite_T2a"] = 0.892

    # ---- the ladder --------------------------------------------------------------
    allrows = np.ones(len(E28), bool)
    pos = (E28.inc > 0).to_numpy()
    acs28 = E28.acs.to_numpy()
    E25 = E                                                     # ages 25-34 pooled
    acs25 = E25.acs.to_numpy()
    rungs = [
        ("1", "CFR: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28", E28, "inc_cap", allrows, "mean", False),
        ("2", "Drop zeros (earnings > 0), capped earnings / mean", E28, "inc_cap", pos, "mean", False),
        ("3", "Log earnings, earnings > 0", E28, "learn", pos, "none", False),
        ("4", "Log earnings, ACS-style sample (full-year, full-time, >= $1,000 in 2010$, hgc >= 9)",
         E28, "learn", acs28, "none", False),
        ("5", "Log hourly wage, ACS-style sample", E28, "lw", acs28, "none", False),
        ("6", "Log hourly wage, ACS-style sample, ages 25-34 pooled", E25, "lw", acs25, "none", True),
    ]
    ladder: dict[str, Any] = {}
    for weighted in (True, False):
        wl = "weighted" if weighted else "unweighted"
        for rid, label, Eb, ycol, mask, tr, pooled in rungs:
            fr = prep(Eb, ycol, mask, pooled)
            for kind in ("single", "avg"):
                score = "comp"
                for spec, ctrl, lag, restrict in SPECS[kind]:
                    r = one_reg(fr[kind], ycol, score, weighted, ctrl, kind, tr, pooled, lag, restrict)
                    R = rel[wl][kind]
                    if lag is None:                     # disattenuation is defined without a lag control
                        r["disatt_coef"] = r["coef"] / np.sqrt(R)
                        r["disatt_se"] = r["se"] / np.sqrt(R)
                    ladder[f"{wl}|{rid}|{kind}|{spec}"] = {
                        **r, "rung": rid, "label": label, "y": ycol, "weighting": wl,
                        "design": kind, "controls": ctrl, "lag": lag}
        LOG.info("  ladder (%s) done", wl)
    out["T6a.ladder"] = ladder

    # ---- household-cluster bootstrap of the disattenuated slopes (weighted, base spec) ---
    # the analytic SEs above treat the reliability as known; here it is re-estimated in each draw
    rng = np.random.default_rng(seed)
    hhu, _ = np.unique(C.hh[elig], return_inverse=True)
    hh_of = {h: i for i, h in enumerate(np.unique(C.hh[elig]))}
    nh = len(hh_of)
    boot_jobs = {}
    for rid, ycol, Eb, mask, pooled in (("5", "lw", E28, acs28, False), ("6", "lw", E25, acs25, True)):
        fr = prep(Eb, ycol, mask, pooled)
        for kind in ("single", "avg"):
            f = fr[kind]
            f = f[np.isfinite(f.comp) & np.isfinite(f.wm)].copy()
            w0 = weights(f, True)
            X = rung_design(C, f, np.zeros(len(f)), True, False, pooled, single=(kind == "single"))
            fb = base[kind][np.isfinite(base[kind].wm)]
            Xt = fb[["z0", "z1", "z2", "z3"]].to_numpy() if kind == "single" else Tavg[fb.child.to_numpy()]
            boot_jobs[(rid, kind)] = (f.comp.to_numpy(), f[ycol].to_numpy(), X, w0,
                                      np.array([hh_of[h] for h in f.hh]), Xt, weights(fb, True),
                                      np.array([hh_of[h] for h in fb.hh]))
    draws = {k: [] for k in boot_jobs}
    for b in range(B):
        cm = np.bincount(rng.integers(0, nh, nh), minlength=nh).astype(float)
        for k, (sc, y, X, w0, hi, Xt, wb, hib) in boot_jobs.items():
            w = w0 * cm[hi]
            X = X.copy()
            X[:, 1] = std_w(sc, w)
            coef = na.wls(y, X, w)[1]
            try:
                R, _ = design_reliability(Xt, wb * cm[hib])
                draws[k].append(coef / np.sqrt(R))
            except Exception:
                draws[k].append(np.nan)
    out["T6a.bootstrap_disattenuated"] = {
        f"{rid}|{kind}": {"point": ladder[f"weighted|{rid}|{kind}|base"]["disatt_coef"],
                          "se_analytic": ladder[f"weighted|{rid}|{kind}|base"]["disatt_se"],
                          "se_boot": float(np.nanstd(v, ddof=1)),
                          "ci95": [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))],
                          "B": B}
        for (rid, kind), v in draws.items()}
    LOG.info("  disattenuated-slope bootstrap done (B=%d)", B)

    # ---- subject-specific single-year scores under the CFR definition -------------
    subj = {}
    fr = prep(E28, "inc_cap", allrows, False)["single"]
    for weighted in (True, False):
        for sc in ("math", "read", "comp"):
            for spec, ctrl, lag, restrict in SPECS["single"]:
                r = one_reg(fr, "inc_cap", sc, weighted, ctrl, "single", "mean", False, lag, restrict)
                subj[f"{'weighted' if weighted else 'unweighted'}|{sc}|{spec}"] = r
    out["T6a.subject_single_year_cfr"] = subj

    # ---- sensitivities of rung 1 and rung 5 (single-year, avg; weighted, ctrl) ------
    sens: dict[str, Any] = {}
    # (a) wider earnings-age window 26-30 (closest to 28)
    E2630 = E[E.age.between(26, 30)].copy()
    E2630["d"] = (E2630.age - 28).abs()
    E2630["o"] = -E2630.age
    E2630 = E2630.sort_values(["child", "d", "o"]).drop_duplicates("child")
    # (b) stricter full-year: 12 covered months
    E28s = E28.copy()
    E28s["acs"] = E28s.acs & (E28s.months >= 12)
    E28s["lw"] = np.where(E28s.acs, E28s.lw, np.nan)
    # (c) earnings-year age exactly 28 at interview
    E28x = E28[E28.age == 28]
    for lab, Eb, ycol, mask, tr in (
            ("rung1_uncapped", E28, "inc_real", allrows, "mean"),
            ("rung1_age26_30", E2630, "inc_cap", np.ones(len(E2630), bool), "mean"),
            ("rung1_interview_age28_only", E28x, "inc_cap", np.ones(len(E28x), bool), "mean"),
            ("rung5_age26_30", E2630, "lw", E2630.acs.to_numpy(), "none"),
            ("rung5_12_months", E28s, "lw", E28s.acs.to_numpy(), "none"),
            ("rung5_interview_age28_only", E28x, "lw", E28x.acs.to_numpy(), "none")):
        fr = prep(Eb, ycol, mask, False)
        for kind in ("single", "avg"):
            for spec, ctrl, lag, restrict in ((("base", False, None, False), ("ctrl", True, None, False),
                                               ("lag", False, "subj", False), ("ctrl_lag", True, "subj", False))
                                              if kind == "single" else (("base", False, None, False),
                                                                        ("ctrl", True, None, False))):
                r = one_reg(fr[kind], ycol, "comp", True, ctrl, kind, tr, False, lag, restrict)
                if lag is None:
                    R = rel["weighted"][kind]
                    r["disatt_coef"] = r["coef"] / np.sqrt(R)
                    r["disatt_se"] = r["se"] / np.sqrt(R)
                sens[f"{lab}|{kind}|{spec}"] = r
    # (d) lag control as a cubic in the previous round's COMPOSITE (rather than math and reading)
    fr = prep(E28, "inc_cap", allrows, False)["single"]
    for ctrl in (False, True):
        sens[f"rung1_lag_composite_cubic|single|{'ctrl_lag' if ctrl else 'lag'}"] = one_reg(
            fr, "inc_cap", "comp", True, ctrl, "single", "mean", False, "comp", False)
    out["T6a.sensitivities"] = sens

    # ---- diagnostics of the constructed wage sample ----------------------------------
    d = E28
    out["T6a.wage_sample_age28"] = {
        "n": int(len(d)),
        "share_zero_earnings": float((d.inc == 0).mean()),
        "share_positive": float((d.inc > 0).mean()),
        "share_fy_months_ge_11": float(d.fy.mean()),
        "share_ft_hours_ge_30_of_known": float(d.ft[np.isfinite(d.tot)].mean()),
        "share_hours_known": float(np.isfinite(d.tot).mean()),
        "share_acs_style": float(d.acs.mean()),
        "n_acs_style": int(d.acs.sum()),
        "mean_earnings_2010usd": float(d.inc_real.mean()),
        "mean_capped_earnings_2010usd": float(d.inc_cap.mean()),
        "share_uncapped_above_100k": float((d.inc_real > 100000).mean()),
        "corr_lw_learn_acs": float(np.corrcoef(d.lw[d.acs], d.learn[d.acs])[0, 1]) if d.acs.sum() > 5 else None,
        "weighted_share_inc0_months_ge_1": float(((d.inc == 0) & (d.months >= 1)).mean()),
        "share_inc_positive_months_0": float(((d.inc > 0) & (d.months == 0)).mean()),
    }


# =============================================================================
# T3: residential transitions
# =============================================================================

def tercile_of(x: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, tuple[float, float]]:
    q1, q2 = na.wquantile(x, w, 1 / 3), na.wquantile(x, w, 2 / 3)
    t = np.where(~np.isfinite(x), -1, np.where(x <= q1, 0, np.where(x <= q2, 1, 2)))
    return t.astype(int), (q1, q2)


class Trans:
    """Origin (city/suburb/non-MSA) x adult residence transition data.

    orig[p]      0 city, 1 suburb, 2 not in MSA/SMSA, -1 excluded (unknown, missing, age)
    rows: person p, adult type a (0/1/2), row weight rw (already divided by the person's
          number of adult rows, so each person carries the mean of his or her weights)
    """

    def __init__(self, name: str, n: int, cluster: np.ndarray, orig: np.ndarray, origw: np.ndarray,
                 income: np.ndarray, inc_pop: np.ndarray, rows: pd.DataFrame,
                 note: str = "") -> None:
        self.name, self.n, self.cluster, self.orig = name, n, cluster, orig
        self.origw, self.income, self.inc_pop, self.rows, self.note = origw, income, inc_pop, rows, note
        self.p = rows.p.to_numpy()
        self.a = rows.a.to_numpy()
        self.rw = rows.rw.to_numpy()
        _, self.cl_idx = np.unique(cluster, return_inverse=True)
        self.ncl = int(self.cl_idx.max()) + 1

    def matrices(self, mult_p: np.ndarray, with_terciles: bool = True) -> dict[str, np.ndarray]:
        """Weighted counts M[o, a] (3x3) overall and by income tercile (3x3x3)."""
        mp = mult_p[self.p]
        o = self.orig[self.p]
        wt = self.rw * mp
        ok = o >= 0
        res = {}
        M = np.zeros((3, 3))
        np.add.at(M, (o[ok], self.a[ok]), wt[ok])
        res["all"] = M
        if with_terciles:
            # cut points from the whole origin population (weights x multiplicity)
            pw = self.origw * mult_p * self.inc_pop
            terc, cuts = tercile_of(self.income, pw)
            res["cuts"] = np.array(cuts)
            T = np.zeros((3, 3, 3))
            t = terc[self.p]
            okt = ok & (t >= 0)
            np.add.at(T, (t[okt], o[okt], self.a[okt]), wt[okt])
            res["terc"] = T
        return res


def stats_from(M: np.ndarray) -> dict[str, float]:
    """2x2 (city, suburb) statistics from a 3x3 count matrix."""
    m = M[:2, :2]
    tot = m.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        return {
            "city_to_suburb": m[0, 1] / m[0].sum(), "suburb_to_city": m[1, 0] / m[1].sum(),
            "city_stay": m[0, 0] / m[0].sum(), "suburb_stay": m[1, 1] / m[1].sum(),
            "move_rate": (m[0, 1] + m[1, 0]) / tot,
            "origin_share_city": m[0].sum() / tot,
            "adult_share_city": m[:, 0].sum() / tot,
            "city_to_nonMSA": M[0, 2] / M[0].sum(), "suburb_to_nonMSA": M[1, 2] / M[1].sum(),
            "nonMSA_to_city": M[2, 0] / M[2].sum(), "nonMSA_to_suburb": M[2, 1] / M[2].sum(),
            "weighted_n_2x2": tot,
        }


def t3_estimate(T: Trans, B: int, seed: int) -> dict[str, Any]:
    ones = np.ones(T.n)
    pt = T.matrices(ones)
    res: dict[str, Any] = {"name": T.name, "note": T.note}
    res["n_persons_origin_city"] = int(((T.orig == 0) & np.isin(np.arange(T.n), T.p)).sum())
    res["n_persons_origin_suburb"] = int(((T.orig == 1) & np.isin(np.arange(T.n), T.p)).sum())
    inn = np.isin(np.arange(T.n), T.p)
    res["n_persons_2x2"] = int(((T.orig >= 0) & (T.orig <= 1) & inn & _has_msa_row(T)).sum())
    res["n_person_years_2x2"] = int(((T.orig[T.p] >= 0) & (T.orig[T.p] <= 1) & (T.a <= 1)).sum())
    res["matrix_3x3_rowpct"] = (pt["all"] / pt["all"].sum(axis=1, keepdims=True)).tolist()
    res["stats"] = stats_from(pt["all"])
    res["tercile_cuts"] = pt["cuts"].tolist()
    res["tercile"] = {}
    for k in range(3):
        st = stats_from(pt["terc"][k])
        o = T.orig[T.p]
        t = tercile_of(T.income, T.origw * T.inc_pop)[0][T.p]
        sel = (t == k) & (o >= 0) & (o <= 1) & (T.a <= 1)
        st["n_person_years_2x2"] = int(sel.sum())
        st["n_persons_2x2"] = int(np.unique(T.p[sel]).size)
        res["tercile"][str(k + 1)] = st
    # cluster bootstrap over persons' clusters (households for NLSY79, persons for NLSY97)
    rng = np.random.default_rng(seed)
    keys = ["city_to_suburb", "suburb_to_city", "move_rate", "city_stay", "suburb_stay",
            "origin_share_city", "adult_share_city"]
    dr = {k: [] for k in keys}
    drt = {str(k + 1): {kk: [] for kk in keys} for k in range(3)}
    cutd = []
    diff_t = {"move_rate_T3_minus_T1": [], "city_to_suburb_T3_minus_T1": [], "suburb_to_city_T3_minus_T1": []}
    for b in range(B):
        cm = np.bincount(rng.integers(0, T.ncl, T.ncl), minlength=T.ncl).astype(float)
        mp = cm[T.cl_idx]
        mb = T.matrices(mp)
        st = stats_from(mb["all"])
        for k in keys:
            dr[k].append(st[k])
        cutd.append(mb["cuts"])
        sts = [stats_from(mb["terc"][k]) for k in range(3)]
        for k in range(3):
            for kk in keys:
                drt[str(k + 1)][kk].append(sts[k][kk])
        diff_t["move_rate_T3_minus_T1"].append(sts[2]["move_rate"] - sts[0]["move_rate"])
        diff_t["city_to_suburb_T3_minus_T1"].append(sts[2]["city_to_suburb"] - sts[0]["city_to_suburb"])
        diff_t["suburb_to_city_T3_minus_T1"].append(sts[2]["suburb_to_city"] - sts[0]["suburb_to_city"])
    res["boot_B"] = B
    res["stats_se"] = {k: float(np.nanstd(v, ddof=1)) for k, v in dr.items()}
    res["tercile_se"] = {t: {k: float(np.nanstd(v, ddof=1)) for k, v in d.items()} for t, d in drt.items()}
    res["tercile_diff"] = {}
    pt_t = [stats_from(pt["terc"][k]) for k in range(3)]
    res["tercile_diff"]["move_rate_T3_minus_T1"] = (
        pt_t[2]["move_rate"] - pt_t[0]["move_rate"], float(np.nanstd(diff_t["move_rate_T3_minus_T1"], ddof=1)))
    res["tercile_diff"]["city_to_suburb_T3_minus_T1"] = (
        pt_t[2]["city_to_suburb"] - pt_t[0]["city_to_suburb"],
        float(np.nanstd(diff_t["city_to_suburb_T3_minus_T1"], ddof=1)))
    res["tercile_diff"]["suburb_to_city_T3_minus_T1"] = (
        pt_t[2]["suburb_to_city"] - pt_t[0]["suburb_to_city"],
        float(np.nanstd(diff_t["suburb_to_city_T3_minus_T1"], ddof=1)))
    res["tercile_cuts_se"] = np.nanstd(np.array(cutd), axis=0, ddof=1).tolist()
    return res


def _has_msa_row(T: Trans) -> np.ndarray:
    m = np.zeros(T.n, bool)
    m[T.p[T.a <= 1]] = True
    return m


def collapse_rows(p: np.ndarray, a: np.ndarray, w: np.ndarray, extra: dict | None = None) -> pd.DataFrame:
    """Row weight = round weight / (person's number of adult rows in the sample)."""
    df = pd.DataFrame({"p": p, "a": a, "w": w})
    if extra:
        for k, v in extra.items():
            df[k] = v
    df["rw"] = df.w / df.groupby("p").p.transform("size")
    return df


# ---- NLSY97 ----------------------------------------------------------------------

def msa_type97(x: np.ndarray) -> np.ndarray:
    """CV_MSA: 1 not in MSA, 2 in MSA not in central city (suburb), 3 central city,
    4 in MSA unknown, 5 not in country.  -> 0 city, 1 suburb, 2 non-MSA, -1 unclassifiable."""
    return np.select([x == 3, x == 2, x == 1], [0, 1, 2], default=-1)


def build_trans97(c: na.Cohort, origin: str = "r1", adult: str = "pooled",
                  years: tuple[int, int] | None = None, same_region: bool = False) -> Trans:
    n = c.df.shape[0]
    pid = c.get("PUBID", 1997)
    age0 = c.get("CV_AGE_INT_DATE", 1997)
    w0 = c.get("SAMPLING_WEIGHT_CC", 1997) / 100.0
    inc = c.get("CV_INCOME_GROSS_YR", 1997)
    reg0 = c.get("CV_CENSUS_REGION", 1997)
    if origin == "r1":
        o_raw = c.get("CV_MSA", 1997)
        reg_o = reg0
        okage = (age0 >= ORIGIN_AGE_97[0]) & (age0 <= ORIGIN_AGE_97[1])
    elif origin == "age12_parent":
        o_raw = c.get("CV_MSA_AGE_12", 1997)
        reg_o = c.get("CV_CENSUS_REGION_AGE_12", 1997)
        okage = np.ones(n, bool)
    elif origin == "age12_youth":
        o_raw = c.get("CV_MSA_AGE_12_YCHR", "XRND") if c.has("CV_MSA_AGE_12_YCHR", "XRND") else \
            c.get("CV_MSA_AGE_12_YCHR", 1997)
        reg_o = c.get("CV_CENSUS_REGION_AGE_12_YCHR", "XRND") if c.has("CV_CENSUS_REGION_AGE_12_YCHR", "XRND") \
            else c.get("CV_CENSUS_REGION_AGE_12_YCHR", 1997)
        okage = np.ones(n, bool)
    else:
        raise ValueError(origin)
    orig = np.where(np.isfinite(o_raw) & okage & np.isfinite(w0) & (w0 > 0), msa_type97(o_raw), -1)
    # income terciles are cut over everyone aged 12-16 in round 1 with valid income
    inc_pop = (np.isfinite(inc) & np.isfinite(w0) & (w0 > 0)
               & (age0 >= ORIGIN_AGE_97[0]) & (age0 <= ORIGIN_AGE_97[1])).astype(float)
    P, A, W = [], [], []
    yrs = sorted(int(y) for y in c.by_q["CV_MSA"] if str(y).isdigit())
    for Y in yrs:
        if years is not None and not (years[0] <= Y <= years[1]):
            continue
        if not (c.has("CV_AGE_INT_DATE", Y) and c.has("SAMPLING_WEIGHT_CC", Y)):
            continue
        age = c.get("CV_AGE_INT_DATE", Y)
        m = c.get("CV_MSA", Y)
        w = c.get("SAMPLING_WEIGHT_CC", Y) / 100.0
        ok = (age >= ADULT_AGE[0]) & (age <= ADULT_AGE[1]) & np.isfinite(m) & np.isfinite(w) & (w > 0)
        if adult == "first25_27":          # young-adult window: rounds at 25-27 only
            ok &= (age <= 27)
        elif adult == "near30":            # rounds at 28-32
            ok &= (age >= 28) & (age <= 32)
        ok &= (orig >= 0)
        if same_region:
            ok &= (c.get("CV_CENSUS_REGION", Y) == reg_o)
        t = msa_type97(m)
        ok &= t >= 0
        i = np.flatnonzero(ok)
        P.append(i)
        A.append(t[i])
        W.append(w[i])
    p, a, w = (np.concatenate(v) for v in (P, A, W))
    df = collapse_rows(p, a, w)
    return Trans(f"NLSY97 {origin}/{adult}", n, pid, orig, w0, inc, inc_pop, df)


# ---- NLSY79 ----------------------------------------------------------------------

def smsa_type79(x: np.ndarray, Y: int) -> np.ndarray:
    """SMSARES: <= 1998: 0 not in SMSA, 1 SMSA not central city, 2 central city not
    known, 3 central city.  >= 2000: 1 not, 2 suburb, 3 central city, 4 unknown."""
    if Y <= 1998:
        return np.select([x == 3, x == 1, x == 0], [0, 1, 2], default=-1)
    return np.select([x == 3, x == 2, x == 1], [0, 1, 2], default=-1)


def build_trans79(c: na.Cohort, years: tuple[int, int] | None = None) -> Trans:
    n = c.df.shape[0]
    sample = c.get("SAMPLE_ID", 1979).astype(int)
    pop = np.isin(sample, na.POP_SAMPLES)
    hh = c.get("HHID", 1979)
    age0 = c.get("AGEATINT", 1979)
    w0 = c.get("SAMPWEIGHT", 1979) / 100.0
    o_raw = c.get("SMSARES", 1979)
    inc = c.get("TNFI_TRUNC", 1979)
    okage = (age0 >= ORIGIN_AGE_79[0]) & (age0 <= ORIGIN_AGE_79[1]) & pop
    orig = np.where(np.isfinite(o_raw) & okage & np.isfinite(w0) & (w0 > 0), smsa_type79(o_raw, 1979), -1)
    inc_pop = (np.isfinite(inc) & okage & np.isfinite(w0) & (w0 > 0)).astype(float)
    P, A, W = [], [], []
    for Y in sorted(int(y) for y in c.by_q["SMSARES"] if str(y).isdigit()):
        if not (c.has("AGEATINT", Y) and c.has("SAMPWEIGHT", Y)):
            continue
        if years is not None and not (years[0] <= Y <= years[1]):
            continue
        age = c.get("AGEATINT", Y)
        m = c.get("SMSARES", Y)
        w = c.get("SAMPWEIGHT", Y) / 100.0
        t = smsa_type79(m, Y)
        ok = ((age >= ADULT_AGE[0]) & (age <= ADULT_AGE[1]) & np.isfinite(m) & np.isfinite(w) & (w > 0)
              & (t >= 0) & (orig >= 0) & pop)
        i = np.flatnonzero(ok)
        P.append(i)
        A.append(t[i])
        W.append(w[i])
    p, a, w = (np.concatenate(v) for v in (P, A, W))
    df = collapse_rows(p, a, w)
    return Trans("NLSY79 1979 (ages 14-16) vs 25-34", n, hh, orig, w0, inc, inc_pop, df,
                 note="SMSA 'central city not known' excluded from the 2x2")


def origin_coverage(c: na.Cohort, kind: str) -> dict[str, Any]:
    """How many eligible children fall in each origin category (before adult matching)."""
    n = c.df.shape[0]
    if kind == "97":
        age0 = c.get("CV_AGE_INT_DATE", 1997)
        w0 = c.get("SAMPLING_WEIGHT_CC", 1997) / 100.0
        m = c.get("CV_MSA", 1997)
        ok = (age0 >= ORIGIN_AGE_97[0]) & (age0 <= ORIGIN_AGE_97[1]) & np.isfinite(w0) & (w0 > 0)
        lab = {1: "not in MSA", 2: "MSA, suburb", 3: "MSA, central city", 4: "MSA, unknown", 5: "not in country"}
    else:
        pop = np.isin(c.get("SAMPLE_ID", 1979).astype(int), na.POP_SAMPLES)
        age0 = c.get("AGEATINT", 1979)
        w0 = c.get("SAMPWEIGHT", 1979) / 100.0
        m = c.get("SMSARES", 1979)
        ok = (age0 >= ORIGIN_AGE_79[0]) & (age0 <= ORIGIN_AGE_79[1]) & pop & np.isfinite(w0) & (w0 > 0)
        lab = {0: "not in SMSA", 1: "SMSA, suburb", 2: "SMSA, central city not known", 3: "SMSA, central city"}
    out = {"n_eligible": int(ok.sum()), "weighted_shares": {}, "n": {}}
    tot = w0[ok].sum()
    for k, v in lab.items():
        s = ok & (m == k)
        out["weighted_shares"][v] = float(w0[s].sum() / tot)
        out["n"][v] = int(s.sum())
    return out


def switch_pairs(c: na.Cohort, kind: str) -> list[dict[str, Any]]:
    """Round-to-round switching between central city and suburb among persons who are
    classified in an MSA at both consecutive rounds and are >= 18 at the first.  A
    jump at a change of MSA standard (NLSY97: 2003->2004 from 1990 to 2000 standards,
    2011->2013 to 2010 standards) measures mechanical reclassification."""
    rows = []
    if kind == "97":
        yrs = sorted(int(y) for y in c.by_q["CV_MSA"] if str(y).isdigit())
        typ = lambda Y: msa_type97(c.get("CV_MSA", Y))
        age = lambda Y: c.get("CV_AGE_INT_DATE", Y)
        wt = lambda Y: c.get("SAMPLING_WEIGHT_CC", Y) / 100.0
        pop = np.ones(c.df.shape[0], bool)
    else:
        yrs = sorted(int(y) for y in c.by_q["SMSARES"] if str(y).isdigit() and int(y) <= 1998)
        typ = lambda Y: smsa_type79(c.get("SMSARES", Y), Y)
        age = lambda Y: c.get("AGEATINT", Y)
        wt = lambda Y: c.get("SAMPWEIGHT", Y) / 100.0
        pop = np.isin(c.get("SAMPLE_ID", 1979).astype(int), na.POP_SAMPLES)
    for Y1, Y2 in zip(yrs[:-1], yrs[1:]):
        t1, t2, a1, w2 = typ(Y1), typ(Y2), age(Y1), wt(Y2)
        ok = pop & (t1 >= 0) & (t1 <= 1) & (t2 >= 0) & (t2 <= 1) & (a1 >= 18) & np.isfinite(w2) & (w2 > 0)
        if ok.sum() < 30:
            continue
        sw = (t1 != t2)[ok].astype(float)
        w = w2[ok]
        pr = na.wmean(sw, w)
        neff = w.sum() ** 2 / (w ** 2).sum()
        rows.append({"from": Y1, "to": Y2, "gap_years": Y2 - Y1, "n": int(ok.sum()), "switch_rate": pr,
                     "se_approx": float(np.sqrt(pr * (1 - pr) / neff)),
                     "share_city_at_to": na.wmean((t2[ok] == 0).astype(float), w)})
    return rows


def unknown_share79(c: na.Cohort) -> list[dict[str, Any]]:
    """Weighted share of SMSA residents at ages 25-34 whose central-city status is unknown."""
    pop = np.isin(c.get("SAMPLE_ID", 1979).astype(int), na.POP_SAMPLES)
    rows = []
    for Y in sorted(int(y) for y in c.by_q["SMSARES"] if str(y).isdigit() and int(y) <= 1998):
        age, w, m = c.get("AGEATINT", Y), c.get("SAMPWEIGHT", Y) / 100.0, c.get("SMSARES", Y)
        t = smsa_type79(m, Y)
        ok = pop & (age >= 25) & (age <= 34) & np.isfinite(w) & (w > 0) & np.isfinite(m) & (t != 2)
        if ok.sum() < 30:
            continue
        rows.append({"year": Y, "n": int(ok.sum()),
                     "unknown_share_of_smsa_residents": na.wmean((t[ok] == -1).astype(float), w[ok])})
    return rows


def t3(c97: na.Cohort, c79: na.Cohort, B: int, seed: int, out: dict) -> None:
    LOG.info("== T3: residential transitions ==")
    res: dict[str, Any] = {}
    LOG.info("  NLSY97 main")
    main97 = build_trans97(c97)
    res["nlsy97_main"] = t3_estimate(main97, B, seed)
    res["nlsy97_origin_coverage"] = origin_coverage(c97, "97")
    variants = {
        "nlsy97_origin_age12_parent": dict(origin="age12_parent"),
        "nlsy97_origin_age12_youth": dict(origin="age12_youth"),
        "nlsy97_same_region": dict(same_region=True),
        "nlsy97_adult_first25_27": dict(adult="first25_27"),
        "nlsy97_adult_near30": dict(adult="near30"),
        "nlsy97_adult_years_2004_2011": dict(years=(2004, 2011)),
        "nlsy97_adult_years_2013_2019": dict(years=(2013, 2019)),
    }
    for k, kw in variants.items():
        try:
            LOG.info("  %s", k)
            res[k] = t3_estimate(build_trans97(c97, **kw), max(B // 2, 20), seed + 1)
        except Exception as exc:                                  # a sparse variant must not sink the run
            LOG.warning("  %s failed: %s", k, exc)
            res[k] = {"error": str(exc)}
    LOG.info("  NLSY79")
    m79 = build_trans79(c79)
    res["nlsy79_main"] = t3_estimate(m79, B, seed + 2)
    res["nlsy79_origin_coverage"] = origin_coverage(c79, "79")
    for k, yy in (("nlsy79_adult_years_le_1996", (1979, 1996)), ("nlsy79_adult_years_1987_1990", (1987, 1990))):
        try:
            LOG.info("  %s", k)
            res[k] = t3_estimate(build_trans79(c79, years=yy), max(B // 2, 20), seed + 3)
        except Exception as exc:
            LOG.warning("  %s failed: %s", k, exc)
            res[k] = {"error": str(exc)}
    res["nlsy97_switch_pairs"] = switch_pairs(c97, "97")
    res["nlsy79_switch_pairs"] = switch_pairs(c79, "79")
    res["nlsy79_unknown_share_adult"] = unknown_share79(c79)
    out["T3"] = res


# =============================================================================
# Report
# =============================================================================

def tab(rows: list[list[str]], head: list[str]) -> str:
    s = "| " + " | ".join(head) + " |\n|" + "|".join("---" for _ in head) + "|\n"
    return s + "\n".join("| " + " | ".join(r) + " |" for r in rows) + "\n"


def write_md(out: dict, meta: dict, path: Path) -> None:
    L: list[str] = []
    L.append("# NLSY extensions: CFR replication in the CNLSY (T6a) and residential transitions (T3)\n")
    L.append(f"Generated {meta['generated']} by `nlsy_extensions.py`. T3: cluster bootstrap, B = {meta['B']}. "
             "T6a: analytic CR1 SEs clustered by the NLSY79 mother's household, plus a household bootstrap "
             f"(B = {max(meta['B'] // 2, 20)}) for the disattenuated slopes. Seed {meta['seed']}. "
             "Data: `raw/nlsy79_ability`, `raw/nlscya_ability`, `raw/nlsy97_ability` from `fetch_data.py --sources nlsy`.\n")
    if "T6a.ladder" in out:
        L += t6a_md(out)
    if "T3" in out:
        L += t3_md(out)
    path.write_text("\n".join(L))


def t6a_md(out: dict) -> list[str]:
    L = []
    lad, rel = out["T6a.ladder"], out["T6a.reliability"]
    smp = out["T6a.sample"]
    L.append("## A. T6a: CFR replication in the CNLSY\n")
    L.append("Children of NLSY79 women, Young Adult survey, earnings (wages and salary, prior calendar year) at the "
             "interview at 27-29 closest to 28 (ties to the older age), regressed on **observed** PIAT/PPVT "
             "composites at ages 8-14. Scores are rank-normed within 3-month age cells (as in `nlsy_ability.py`); "
             "the composite is the mean of the available tests (at least 2 of MATH, RECOG, COMP, PPVT), then "
             "standardized to weighted SD 1 in each regression sample, so a coefficient is per observed score SD. "
             "*Single-year*: one observation per child-round score (all rounds at ages 8-14, so the child's earnings "
             "row repeats), the CFR design. *Round average*: the child's average over rounds (no lag defined). "
             "Standard errors cluster on the mother's household. Weights: mother's 1979 weight split so each child "
             "counts equally (the `nlsy_ability.py` dyad convention); unweighted below.\n")
    L.append("### Which CFR / NLSY79 number each column is comparable to\n")
    L.append("CFR (2014b, Appendix Table 3): earnings at 28 are individual W-2 wages, zeros included (33.1% have none at "
             "28), capped at $100,000 and expressed relative to the mean ($21,622). Unconditionally, one student SD "
             "of the score goes with $7,709, **36% of the mean** (0.305 in logs, = log 1.357). Conditioning on "
             "cubics in prior-year math and English scores and teacher fixed effects gives $2,585, **12% of the "
             "mean** (13.9% math, 10.1% English); the model's chi = log 1.12 = 0.113 is that conditional number.\n")
    L.append("- **Unconditional columns** (`base`, `controls`, round average; Table 1) have no lagged-score control. "
             "Compare them with CFR's unconditional **36%** (steps 1-2, whose definitions are CFR's) and, for the "
             "log-wage steps 5-6 disattenuated, with the NLSY79 latent slope **0.213** "
             "(controls: age, year, sex; ACS sample; no lagged score). The conditional CFR 12% is not the "
             "comparison for any of these. The `controls` columns add parental characteristics, which CFR's "
             "unconditional number does not have; they fall between CFR's unconditional and conditional numbers.\n"
             "- **Lag-conditional columns** (Table 2) add a cubic in the child's PIAT math and reading scores "
             "from the previous assessment round (two years earlier when the child was tested in consecutive "
             "biennial rounds), the analogue of CFR's cubics in prior-year scores (there is no teacher fixed "
             "effect in the NLSY). Compare these with CFR's conditional **12%** (steps 1-2 only, which share CFR's "
             "definitions), and with chi = 0.113 in logs (steps 3-6 are logs, not earnings/mean, so the "
             "comparison to 12% is loose). The reference for the lag sample is `base (lag sample)`, the "
             "unconditional coefficient on the same rows, so the effect of the control is separated from "
             "the effect of the sample change. Disattenuation is not applied to lag-conditional coefficients: "
             "the lag absorbs part of the measurement error and the raw reliability no longer applies.\n")
    L.append("**Controls.** *Base*: earnings-year and age-at-earnings dummies, female (and age-at-test dummies for "
             "single-year), the specification of the NLSY79 slope. *Controls*: adds race (Hispanic, Black, other), "
             "mother's highest grade, mother's age at birth and log mother's family income (mean over the child's "
             "ages 8-14, 2010$, with a missing indicator), CFR-style parental characteristics; mother's AFQT is "
             "deliberately excluded.\n")
    L.append("**Sample.** " + f"{smp['ya_children_25_34']} children have a YA interview at 25-34 "
             f"({smp['ya_person_rounds_25_34']} person-rounds); {smp['age28_children']} have one at 27-29. "
             f"{smp['children_with_scores_8_14']} children have a composite at ages 8-14 "
             f"(mean {smp['mean_rounds_per_child']:.2f} rounds); {smp['single_year_score_obs']} child-round scores. "
             f"Interview ages at the chosen round: {smp['age28_age_counts']}.\n")
    L.append("**Reliability** (one-factor ULS on the 4 tests in the CFR-definition sample; reliability of the "
             "available-test mean as a measure of the common factor, `comp_reliability_avail`): "
             + f"single-year composite {rel['weighted']['single']:.3f}, round-average composite "
             f"{rel['weighted']['avg']:.3f} (weighted); {rel['unweighted']['single']:.3f} and "
             f"{rel['unweighted']['avg']:.3f} (unweighted). `nlsy_ability.py` reports 0.892 for the child "
             "composite over ages 5-14. The disattenuated slope is the observed slope divided by the "
             "square root of the reliability.\n")

    def cell(k: str, dis: bool = False) -> str:
        r = lad.get(k)
        if r is None:
            return "n/a"
        if dis and "disatt_coef" not in r:
            return "n/a"
        b, se = (r["disatt_coef"], r["disatt_se"]) if dis else (r["coef"], r["se"])
        return f"{b:.3f} ({se:.3f})"

    names = (("1", "CFR definition: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28"),
             ("2", "Drop zeros (earnings > 0), capped earnings / mean"), ("3", "Log earnings, earnings > 0"),
             ("4", "Log earnings, ACS-style sample"), ("5", "Log hourly wage, ACS-style sample"),
             ("6", "Log hourly wage, ACS-style, ages 25-34 pooled"))
    for wl in ("weighted", "unweighted"):
        L.append(f"### Table 1 ({wl}): no lagged-score control\n")
        head = ["step", "definition", "single-year, base", "single-year, controls", "round avg, base",
                "round avg, controls", "N obs / children (single; avg)", "compare with"]
        rows = []
        for rid, name in names:
            comp = {"1": "CFR unconditional 36%", "2": "(no CFR analogue)", "3": "(no CFR analogue)",
                    "4": "(no analogue)", "5": "NLSY79 0.213 after disattenuation", "6": "NLSY79 0.213 after disattenuation"}[rid]
            rows.append([rid, name] + [cell(f"{wl}|{rid}|{kind}|{c}") for kind in ("single", "avg") for c in ("base", "ctrl")]
                        + [f"{lad[f'{wl}|{rid}|single|base']['n_obs']} / {lad[f'{wl}|{rid}|single|base']['n_children']}; "
                           f"{lad[f'{wl}|{rid}|avg|base']['n_obs']} / {lad[f'{wl}|{rid}|avg|base']['n_children']}", comp])
        for rid, name in (("5", "5, per latent SD (observed / sqrt R)"), ("6", "6, per latent SD (observed / sqrt R)")):
            rows.append([f"{rid}d", name] + [cell(f"{wl}|{rid}|{kind}|{c}", True) for kind in ("single", "avg")
                                             for c in ("base", "ctrl")] + ["", "NLSY79 0.213 (SE 0.008)"])
        L.append(tab(rows, head))
        L.append(f"### Table 2 ({wl}): with cubics in the previous round's math and reading scores (single-year design)\n")
        head = ["step", "definition", "base, lag sample", "lag control", "controls + lag control",
                "N obs / children (lag sample)", "compare with"]
        rows = []
        for rid, name in names:
            comp = {"1": "CFR conditional 12%", "2": "(no CFR analogue)", "3": "chi = 0.113 (loose)",
                    "4": "chi = 0.113 (loose)", "5": "chi = 0.113 (loose)", "6": "chi = 0.113 (loose)"}[rid]
            k0 = f"{wl}|{rid}|single|"
            rows.append([rid, name, cell(k0 + "base_lagsample"), cell(k0 + "lag"), cell(k0 + "ctrl_lag"),
                         f"{lad[k0 + 'lag']['n_obs']} / {lad[k0 + 'lag']['n_children']}", comp])
        L.append(tab(rows, head))
    bd = out.get("T6a.bootstrap_disattenuated")
    if bd:
        L.append("### Disattenuated log-wage slopes with the reliability re-estimated (weighted, base spec)\n")
        rows = []
        for k, r in bd.items():
            rid, kind = k.split("|")
            rows.append([f"{rid}d", kind, f"{r['point']:.3f}", f"{r['se_analytic']:.3f}", f"{r['se_boot']:.3f}",
                         f"[{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]"])
        L.append(tab(rows, ["step", "design", "per latent SD", "SE (analytic, R fixed)", "SE (household bootstrap, R re-estimated)",
                            "95% percentile CI"]) + f"B = {next(iter(bd.values()))['B']}. NLSY79: 0.213 (SE 0.008).\n")
    L.append("Steps 1-3 and 6 use the whole age-28 (or 25-34) earnings sample the definition allows; step 4 restricts "
             "to the ACS-style sample without changing the outcome; step 5 divides by constructed hours.\n")
    L.append("### Subject-specific single-year scores under CFR's definition (step 1)\n")
    sub = out["T6a.subject_single_year_cfr"]
    rows = []
    for wl in ("weighted", "unweighted"):
        for sc, nm in (("comp", "composite"), ("math", "PIAT math"), ("read", "reading (RECOG, COMP)")):
            rows.append([wl, nm] + [f"{sub[f'{wl}|{sc}|{c}']['coef']:.3f} ({sub[f'{wl}|{sc}|{c}']['se']:.3f})"
                                    for c in ("base", "ctrl", "base_lagsample", "lag", "ctrl_lag")])
    L.append(tab(rows, ["weights", "score", "base", "controls", "base, lag sample", "lag control", "controls + lag"]))
    L.append("### Sensitivities (weighted; specification as labelled: base, controls, lag control, controls + lag control)\n")
    se = out["T6a.sensitivities"]
    rows = []
    for k, r in se.items():
        lab, kind, spec = k.split("|")
        dis = f"{r['disatt_coef']:.3f} ({r['disatt_se']:.3f})" if "disatt_coef" in r else "n/a"
        rows.append([lab, kind, spec, f"{r['coef']:.3f} ({r['se']:.3f})", dis, f"{r['n_obs']} / {r['n_children']}"])
    L.append(tab(rows, ["variant", "design", "spec", "observed", "per latent SD", "N obs / children"]))
    d = out["T6a.wage_sample_age28"]
    L.append("### Constructed wage sample at ~28\n")
    L.append(f"N = {d['n']}; zero earnings {d['share_zero_earnings']:.3f} (CFR: 0.331); job in >= 11 of 12 reference-year "
             f"months {d['share_fy_months_ge_11']:.3f}; usual hours known at the interview {d['share_hours_known']:.3f} "
             f"(of those, >= 30: {d['share_ft_hours_ge_30_of_known']:.3f}); ACS-style sample "
             f"{d['n_acs_style']} ({d['share_acs_style']:.3f}); mean earnings ${d['mean_earnings_2010usd']:.0f} "
             f"(${d['mean_capped_earnings_2010usd']:.0f} capped at $100,000; {d['share_uncapped_above_100k']:.3f} above the cap), "
             f"2010$; zero earnings with a job in the year {d['weighted_share_inc0_months_ge_1']:.3f}; positive earnings "
             f"with no job in the history {d['share_inc_positive_months_0']:.3f} (unweighted shares).\n")
    L.append("### Caveats for T6a\n")
    L.append("- **Annual weeks and hours are not released for the CNLSY Young Adults**, unlike NLSY79 "
             "(`WKSWK-PCY`, `HRSWK-PCY`) and NLSY97 (`CVC_*_YR_ALL`). ACS-style full-year is proxied by a job in at "
             "least 11 of the 12 reference-year months (job-history start/stop dates; about 48 weeks), full-time by "
             "usual weekly hours at all jobs current on the interview date (`TOTHOURS`, >= 30), and hourly wage "
             "by income / (52 x months/12 x usual hours). Hours at the interview can differ from hours in the "
             "reference year; the sample is therefore restricted to workers still employed at the interview.\n"
             "- Earnings are wages and salary (as in NLSY79), top-coded from 2006; CFR use W-2 wages from tax "
             "records. Zeros are reported zeros (`Q15-5` universe is all YAs).\n"
             "- Test ages 8-14 are PIAT/PPVT raw scores normed within age cells, not grade-level state tests. "
             "CNLSY assessments are biennial, so the 'previous round' is two years earlier (see "
             "`share_lag_gap_2` in the JSON), not the prior grade; there is no teacher fixed effect.\n"
             "- The CNLSY is children of women aged 14-22 in 1979 (born mostly 1975-1992), a different population "
             "from CFR's NYC 1989-2009 cohorts. Mother-family weights assume the NLSY79 weights carry over to children.\n")
    return L


def t3_md(out: dict) -> list[str]:
    L = ["## B. T3: residential transitions between the city and the suburb of an MSA\n"]
    T = out["T3"]
    m97, m79 = T["nlsy97_main"], T["nlsy79_main"]
    s97, e97, s79, e79 = m97["stats"], m97["stats_se"], m79["stats"], m79["stats_se"]
    L.append("**Summary.** Among children in an MSA at both ages, the share living in the other type of location "
             f"as adults (the model's move rate) is **{100 * s97['move_rate']:.1f}% ({100 * e97['move_rate']:.1f}) in "
             f"the NLSY97** (city to suburb {100 * s97['city_to_suburb']:.1f}%, suburb to city "
             f"{100 * s97['suburb_to_city']:.1f}%) and **{100 * s79['move_rate']:.1f}% ({100 * e79['move_rate']:.1f}) in "
             f"the NLSY79** (city to suburb {100 * s79['city_to_suburb']:.1f}%, suburb to city "
             f"{100 * s79['suburb_to_city']:.1f}%). The two are not comparable measurements of one number: the NLSY97 "
             "rate contains a mechanical reclassification from the change of MSA standard between origin and adult "
             "rounds (the switch rate roughly doubles at the 2003-04 break; about 6-8 pp), while the NLSY79 rate "
             "drops the 30-40% of SMSA residents with unknown central-city status (33-41% at ages 25-34 through 1996). Read the NLSY79 figure as the "
             "cleaner point and the NLSY97 figure as an upper bound. By parental-income tercile the move rate is "
             "nearly flat in both cohorts (below); the composition of moves by origin is not.\n")
    L.append("Origin: residence at 12-16, the round-1 (1997) NLSY97 interview (`CV_MSA`, respondents aged 12-16 at "
             "interview) or the 1979 NLSY79 interview (`SMSARES`, aged 14-16). Adult: interview rounds at 25-34; "
             "each person carries the mean of his or her round weights over the rounds observed at 25-34 in an "
             "MSA (the transition share is the weighted share of person-rounds). *City* = MSA central city, "
             "*suburb* = in MSA, not in central city (labelled 'not central city' in the codebook). The "
             "2x2 conditions on being in an MSA at both ages and on a known central-city status. Weights: "
             "NLSY97 round `SAMPLING_WEIGHT_CC`; NLSY79 round `SAMPWEIGHT`. SEs are cluster bootstrap "
             "(persons in NLSY97, households in NLSY79); income tercile cuts are recomputed in each draw. "
             "Parental income terciles: NLSY97 `CV_INCOME_GROSS_YR` (1996 gross household income, topcoded at "
             "2%), NLSY79 `TNFI_TRUNC` (1978 family income), over all respondents of the origin age with valid "
             "income.\n")

    def mat(r: dict) -> str:
        m = np.array(r["matrix_3x3_rowpct"])
        st, se = r["stats"], r["stats_se"]
        rows = [["origin city", f"{100 * (1 - st['city_to_suburb']):.1f}", f"{100 * st['city_to_suburb']:.1f} ({100 * se['city_to_suburb']:.1f})"],
                ["origin suburb", f"{100 * st['suburb_to_city']:.1f} ({100 * se['suburb_to_city']:.1f})", f"{100 * (1 - st['suburb_to_city']):.1f}"]]
        return tab(rows, ["row %", "adult city", "adult suburb"])

    for key, title in (("nlsy97_main", "NLSY97 (main)"), ("nlsy79_main", "NLSY79 (check)")):
        r = T[key]
        st, se = r["stats"], r["stats_se"]
        L.append(f"### {title}\n")
        L.append(mat(r))
        L.append(f"Person-years in the 2x2: {r['n_person_years_2x2']} ({r['n_persons_2x2']} persons). "
                 f"Move rate (share of origin-city plus origin-suburb children living in the other type as adults): "
                 f"**{100 * st['move_rate']:.1f}% ({100 * se['move_rate']:.1f})**. Origin share city "
                 f"{100 * st['origin_share_city']:.1f}% ({100 * se['origin_share_city']:.1f}); adult share city "
                 f"{100 * st['adult_share_city']:.1f}% ({100 * se['adult_share_city']:.1f}).\n")
        cov = T["nlsy97_origin_coverage" if key == "nlsy97_main" else "nlsy79_origin_coverage"]
        L.append("Origin coverage (all eligible respondents, weighted share; N): " + "; ".join(
            f"{k} {100 * v:.1f}% ({cov['n'][k]})" for k, v in cov["weighted_shares"].items()) + ".\n")
        rows = []
        for t in ("1", "2", "3"):
            s, e = r["tercile"][t], r["tercile_se"][t]
            rows.append([f"T{t}", f"{100 * s['city_to_suburb']:.1f} ({100 * e['city_to_suburb']:.1f})",
                         f"{100 * s['suburb_to_city']:.1f} ({100 * e['suburb_to_city']:.1f})",
                         f"{100 * s['move_rate']:.1f} ({100 * e['move_rate']:.1f})",
                         f"{100 * s['origin_share_city']:.1f}", f"{s['n_person_years_2x2']} / {s['n_persons_2x2']}"])
        cuts = r["tercile_cuts"]
        L.append(f"By parental-income tercile (cut points {cuts[0]:,.0f} and {cuts[1]:,.0f} in nominal dollars):\n")
        L.append(tab(rows, ["tercile", "city to suburb %", "suburb to city %", "move rate %", "origin city share %",
                            "person-years / persons"]))
        d = r["tercile_diff"]
        L.append("Top minus bottom tercile: move rate {:.1f} ({:.1f}) pp, city to suburb {:.1f} ({:.1f}) pp, "
                 "suburb to city {:.1f} ({:.1f}) pp.\n".format(
                     100 * d["move_rate_T3_minus_T1"][0], 100 * d["move_rate_T3_minus_T1"][1],
                     100 * d["city_to_suburb_T3_minus_T1"][0], 100 * d["city_to_suburb_T3_minus_T1"][1],
                     100 * d["suburb_to_city_T3_minus_T1"][0], 100 * d["suburb_to_city_T3_minus_T1"][1]))
        M3 = np.array(r["matrix_3x3_rowpct"])
        L.append("Full 3x3 with non-MSA (row %, rows origin city / suburb / non-MSA; columns adult city / suburb / non-MSA): "
                 + "; ".join("[" + ", ".join(f"{100 * x:.1f}" for x in row) + "]" for row in M3) + ".\n")

    L.append("### Variants (move rate; city to suburb; suburb to city; in %)\n")
    rows = []
    for k, r in T.items():
        if not isinstance(r, dict) or not (k.startswith("nlsy97_") or k.startswith("nlsy79_")):
            continue
        if k in ("nlsy97_main", "nlsy79_main") or "coverage" in k or "switch" in k or "unknown" in k:
            continue
        if "error" in r:
            rows.append([k, "error", "", "", ""])
            continue
        st, se = r["stats"], r["stats_se"]
        rows.append([k,
                     f"{100 * st['move_rate']:.1f} ({100 * se['move_rate']:.1f})",
                     f"{100 * st['city_to_suburb']:.1f} ({100 * se['city_to_suburb']:.1f})",
                     f"{100 * st['suburb_to_city']:.1f} ({100 * se['suburb_to_city']:.1f})",
                     f"{r['n_person_years_2x2']} / {r['n_persons_2x2']}"])
    L.append(tab(rows, ["variant", "move rate", "city to suburb", "suburb to city", "person-years / persons"]))
    L.append("### Definitional breaks: round-to-round switching between city and suburb\n")
    L.append("Persons classified as city or suburb at two consecutive rounds (>= 18 at the first), weighted by the second "
             "round's weight. **NLSY97 changes MSA standard between the 2003 and 2004 rounds (to the 2000 standards; the "
             "earlier rounds carry no label in the codebook, presumably the 1990 standards) and between 2011 and 2013 "
             "(2010 standards)**; the switch rate jumps at both, and the city share jumps by 6 pp at the first. The "
             "1997 origin is on the earlier standard and every adult round at 25-34 (2005 onward) is on a later one, "
             "so the NLSY97 transition rate includes this reclassification (the excess switching at the 2003-04 break "
             "is about 6-8 pp against the 9-12% annual switching either side).\n")
    rows = [[f"{r['from']}-{r['to']}", str(r["gap_years"]), str(r["n"]), f"{100 * r['switch_rate']:.1f} ({100 * r['se_approx']:.1f})",
             f"{100 * r['share_city_at_to']:.1f}"] for r in T["nlsy97_switch_pairs"]]
    L.append("NLSY97:\n")
    L.append(tab(rows, ["rounds", "gap (yrs)", "N", "switch rate % (approx SE)", "share city at second round %"]))
    rows = [[f"{r['from']}-{r['to']}", str(r["gap_years"]), str(r["n"]), f"{100 * r['switch_rate']:.1f} ({100 * r['se_approx']:.1f})",
             f"{100 * r['share_city_at_to']:.1f}"] for r in T["nlsy79_switch_pairs"]]
    L.append("NLSY79 (SMSA definitions: no spike in switching, but the 1998 round has no 'central city not known' "
             "category at all, so the classified sample changes there; the falling city share is the cohort moving to "
             "the suburbs with age):\n")
    L.append(tab(rows, ["rounds", "gap (yrs)", "N", "switch rate % (approx SE)", "share city at second round %"]))
    rows = [[str(r["year"]), str(r["n"]), f"{100 * r['unknown_share_of_smsa_residents']:.1f}"]
            for r in T["nlsy79_unknown_share_adult"]]
    L.append("NLSY79: share of SMSA residents aged 25-34 whose central-city status is 'not known' (excluded from the 2x2):\n")
    L.append(tab(rows, ["year", "N", "unknown share %"]))
    L.append("### Caveats for T3\n")
    L.append("- **Cross-MSA movers cannot be separated** in the public files (no county/CBSA code). A child who leaves a "
             "metro area and lands in another metro's central city or suburb counts as a city/suburb transition (or a "
             "stay). The model's move rate is within one commuting zone, so a cross-MSA move that changes type inflates "
             "the rate and one that keeps type deflates it (a cross-MSA 'stayer' counts as no move). The same-region "
             "variant drops cross-region movers only (a region change requires leaving the MSA); same-region "
             "cross-MSA moves remain.\n"
             "- **Central-city definitions change over time.** The NLSY97 origin (1997) is on the pre-2004 (unlabelled, "
             "presumably 1990) standard; adult rounds use 2000 standards through 2011 and 2010 standards for 2013-2019 "
             "(splits by adult year above). New CBSA "
             "definitions reclassify places between 'central city', 'suburb' and 'not in MSA', mechanically moving "
             "people across categories without moves.\n"
             "- NLSY79: 'SMSA, central city not known' is 29% of 1979 SMSA residents aged 14-16 and 33-41% of SMSA "
             "residents at 25-34 through 1996 (none in 1998, when the classification became complete); those are "
             "excluded from the 2x2 (selection into unclassifiable SMSAs), and 1979 age-14-16 residence is the "
             "interview residence (at home), not residence at 14.\n"
             "- Transition shares depend on the ages: residence at 25-34 is the average over person-rounds; parental "
             "income is a single noisy year, missing for a share of respondents, and terciles use the whole "
             "origin-age population, so the MSA sample within a tercile is not balanced.\n")
    return L


# =============================================================================
# Main
# =============================================================================

def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260929)
    p.add_argument("--quick", action="store_true", help="B = 50")
    p.add_argument("--only", choices=["t6a", "t3"], default=None)
    p.add_argument("--md-only", action="store_true", help="rewrite the .md from the saved .json")
    p.add_argument("--outdir", default=str(HERE / "estimates"))
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    B = 50 if args.quick else args.bootstrap
    t0 = time.time()
    out: dict[str, Any] = {}
    if args.md_only:
        outdir = Path(args.outdir).expanduser().resolve()
        saved = json.loads((outdir / "nlsy_extensions.json").read_text())
        meta = saved.pop("meta")
        write_md(saved, meta, outdir / "nlsy_extensions.md")
        return 0
    if args.only != "t3":
        c79, ccy = na.Cohort("nlsy79", HERE), na.Cohort("nlscya", HERE)
        D = na.build_n79(c79)
        D.kids = na.build_kids(ccy, D)
        t6a(c79, ccy, D, D.kids, out, B=max(B // 2, 20), seed=args.seed)
    if args.only != "t6a":
        c97 = na.Cohort("nlsy97", HERE)
        c79 = na.Cohort("nlsy79", HERE)
        t3(c97, c79, B, args.seed, out)
    meta = {"generated": dt.datetime.now().isoformat(timespec="seconds"),
            "script": "data/spatial/nlsy_extensions.py", "B": B, "seed": args.seed}
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "nlsy_extensions.json").write_text(json.dumps({"meta": meta, **out}, indent=1, default=float))
    write_md(out, meta, outdir / "nlsy_extensions.md")
    LOG.info("done in %.1fs -> %s", time.time() - t0, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
