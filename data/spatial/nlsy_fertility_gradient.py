#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_fertility_gradient.py -- how the probability of EVER having a child varies with
ability (AFQT) and with own earnings / family income, NLSY79.

Calibrates the ability-dependent child shock of the spatial teachers model:
    P(child | z) = logistic(a + b * x),  x = latent log ability in population-SD units,
with the population mean fixed at 0.802 (estimates/nlsy_fertility.json).

    cd data/spatial && python3 nlsy_fertility_gradient.py            # B = 200
    python3 nlsy_fertility_gradient.py --quick                       # B = 20

Nothing is re-extracted: uses raw/nlsy79_fertility.parquet (nlsy_fertility.py) and
raw/nlsy79_ability.parquet (nlsy_ability.py).  Both pipelines are imported, not
re-implemented: `afqt_c` and `R_comp_pop` come from nlsy_ability.block_a, and the
completed-fertility sample / outcome from nlsy_fertility.build + specs("main").

Sample: completed fertility (last valid NUMKID obs at age >= 45), POP_SAMPLES,
last-round SAMPWEIGHT/100; merged to the ability panel on CASEID.
Inference: HHID-cluster bootstrap (multiplicity weights); afqt_c and R_comp_pop are
held fixed across draws (their own sampling error, SE(R_comp_pop) = 0.002, is ignored).
Writes estimates/nlsy_fertility_gradient.{json,md} (no microdata).

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

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402
import nlsy_fertility as nf  # noqa: E402

LOG = logging.getLogger("nlsy_fertility_gradient")
SEED = 20261019
GROUPS = ("pooled", "women", "men")
EARN_AGE = (30, 44)          # ages at interview for the earnings / income windows
MIN_YEARS = 3                # min valid person-years for permanent earnings / income
MIN_YEARS_WAGE = 2           # min valid wage person-years for the permanent hourly wage
MIN_HOURS = 500              # annual hours for a wage-rate observation
BLOGWAGE = 0.213             # log earnings per latent SD (block C target)
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(80)
GH_W = GH_W / GH_W.sum()     # E over x ~ N(0,1)


# =============================================================================
# Weighted helpers
# =============================================================================

def wols(y, X, w):
    XtW = X.T * w
    return np.linalg.solve(XtW @ X, XtW @ y)


def wlogit(y, X, w, it=40):
    b = np.zeros(X.shape[1])
    b[0] = np.log(np.average(y, weights=w) / (1 - np.average(y, weights=w)))
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ b))
        g = X.T @ (w * (y - p))
        H = (X.T * (w * p * (1 - p))) @ X
        step = np.linalg.solve(H, g)
        b = b + step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def wrank(v, w):
    """Weighted mid-rank percentile in (0,1)."""
    o = np.argsort(v, kind="stable")
    cw = np.cumsum(w[o])
    r = np.empty(len(v))
    r[o] = (cw - 0.5 * w[o]) / cw[-1]
    return r


def wquint_groups(v, w):
    """Quintile labels 0..4 by weighted cumulative share."""
    return np.minimum((wrank(v, w) * 5).astype(int), 4)


def eiv_two(y, e, x, w, R):
    """y ~ 1 + e + x_lat, x observed with classical error (reliability R).
    Returns (coef on e, coef on STANDARDISED latent x), standardisation = latent SD
    sqrt(R * var_obs(x))."""
    X = np.column_stack([np.ones(len(y)), e, x])
    M = (X.T * w) @ X / w.sum()
    mu = (X.T @ w) / w.sum()
    C = M - np.outer(mu, mu)                                # centred 2nd moments
    Cx = C[1:, 1:].copy()
    Cx[1, 1] *= R                                           # var(x_lat) = R var(x_obs)
    cxy = ((X[:, 1:] * (w * (y - np.average(y, weights=w)))[:, None]).sum(0) / w.sum())
    g = np.linalg.solve(Cx, cxy)
    return float(g[0]), float(g[1] * np.sqrt(Cx[1, 1]))


def deconv_logit(p, s_lat):
    """(a, b) with E[logistic(a + b x)] = p and E[logistic(a+bx) x] = s_lat, x ~ N(0,1)
    (s_lat = covariance of y with the unit-variance latent ability = latent LPM slope)."""
    v = np.array([np.log(p / (1 - p)), s_lat / (p * (1 - p))])
    for _ in range(100):
        q = 1 / (1 + np.exp(-(v[0] + v[1] * GH_X)))
        dq = q * (1 - q)
        F = np.array([(GH_W * q).sum() - p, (GH_W * q * GH_X).sum() - s_lat])
        J = np.array([[(GH_W * dq).sum(), (GH_W * dq * GH_X).sum()],
                      [(GH_W * dq * GH_X).sum(), (GH_W * dq * GH_X ** 2).sum()]])
        step = np.linalg.solve(J, F)
        v = v - step
        if np.max(np.abs(step)) < 1e-12:
            return float(v[0]), float(v[1])
    return np.nan, np.nan


def deconv_logit_a(b, p):
    """Intercept with E[logistic(a + b x)] = p, x ~ N(0,1)."""
    a = np.log(p / (1 - p))
    for _ in range(100):
        q = 1 / (1 + np.exp(-(a + b * GH_X)))
        step = ((GH_W * q).sum() - p) / (GH_W * q * (1 - q)).sum()
        a -= step
        if abs(step) < 1e-13:
            break
    return float(a)


# =============================================================================
# Data
# =============================================================================

def earnings_panels(c: na.Cohort, age_lo: int, age_hi: int):
    """Per-respondent permanent measures from person-years aged lo..hi at interview.
    Real = 2010 dollars via the pipeline's CPI_U; income reported for year Y-1."""
    n = c.df.shape[0]
    years = sorted(int(y) for y in c.by_q["AGEATINT"] if y.isdigit())
    E, F, LW = [np.full((n, len(years)), np.nan) for _ in range(3)]
    ages = np.full((n, len(years)), np.nan)
    for j, Y in enumerate(years):
        age = c.get("AGEATINT", Y)
        ok = (age >= age_lo) & (age <= age_hi)
        if not ok.any() or (Y - 1) not in na.CPI_U:
            continue
        defl = na.CPI_U[2010] / na.CPI_U[Y - 1]
        inc = na._income79(c, Y) * defl
        hrs = c.get("HRSWK-PCY", Y)
        wks = c.get("WKSWK-PCY", Y)
        tnfi = c.get("TNFI_TRUNC", Y) * defl if c.has("TNFI_TRUNC", Y) else np.full(n, np.nan)
        floor = na.MIN_INCOME_2010
        with np.errstate(divide="ignore", invalid="ignore"):
            lw = np.where((inc >= floor) & (hrs >= MIN_HOURS), np.log(inc / hrs), np.nan)
        E[:, j] = np.where(ok, inc, np.nan)
        F[:, j] = np.where(ok, tnfi, np.nan)
        LW[:, j] = np.where(ok, lw, np.nan)
        ages[:, j] = np.where(ok, age, np.nan)
    nE, nF, nW = (np.isfinite(A).sum(1) for A in (E, F, LW))
    with np.errstate(all="ignore"):
        mE = np.nanmean(E, axis=1)
        mF = np.nanmean(F, axis=1)
        mW = np.nanmean(LW, axis=1)
        logE = np.where((nE >= MIN_YEARS) & (mE >= na.MIN_INCOME_2010), np.log(mE), np.nan)
        logF = np.where((nF >= MIN_YEARS) & (mF > 0), np.log(mF), np.nan)
        logW = np.where(nW >= MIN_YEARS_WAGE, mW, np.nan)
    info = dict(years=years, nE=nE, nF=nF, nW=nW, mE=mE, mF=mF,
                zeroE=(nE >= MIN_YEARS) & (mE < na.MIN_INCOME_2010),
                age_min=np.nanmin(ages), age_max=np.nanmax(ages))
    return logE, logF, logW, info


def load():
    c = na.Cohort("nlsy79", HERE)
    D = na.build_n79(c)
    out: dict = {}
    afqt, Rc = na.block_a(D, D.w, out, None)
    P = nf.Panel(HERE / "raw", HERE / "metadata")
    F = nf.build(P)
    sp = nf.specs(F)
    m, wlast, ever = sp[("main", "pooled")]
    cid_f = P.get("CASEID", 1979)
    # merge: ability row for each fertility row, by CASEID
    pos = pd.Series(np.arange(len(D.caseid)), index=D.caseid.astype(int))
    idx = pos.reindex(cid_f.astype(int)).to_numpy()
    matched = np.isfinite(idx)
    idx = np.where(matched, idx, 0).astype(int)
    take = lambda a: np.where(matched, a[idx], np.nan)
    logE, logF, logW, einfo = earnings_panels(c, *EARN_AGE)
    data = dict(
        n=len(cid_f), main=m, matched=matched, w=np.where(m, wlast, 0.0), y=np.where(m, ever, 0.0),
        female=F["female"], hh=F["hh"], pop=F["pop"],
        x=take(afqt), logE=take(logE), logF=take(logF), logW=take(logW),
        Rc=Rc, afqt_pop=afqt, asamp=D.asamp, einfo=einfo, p_main=float(np.average(ever[m], weights=wlast[m])),
        eiv_ok=True,
    )
    for k in ("mE", "mF"):
        data[k] = take(einfo[k])
    data["nE"] = take(einfo["nE"])
    return data


# =============================================================================
# Statistics given a weight vector (point estimate or bootstrap multiplicity)
# =============================================================================

def stats(d: dict, mult: np.ndarray, detail: bool = False) -> dict:
    out: dict = {}
    R = d["Rc"]
    sR = np.sqrt(R)
    f = d["female"]
    gm = {"pooled": np.ones(d["n"], bool), "women": f, "men": ~f}
    w0 = d["w"] * mult
    y = d["y"]

    # ---------------- ability ----------------
    base = (w0 > 0) & np.isfinite(d["x"])
    for g in GROUPS:
        m = base & gm[g]
        w, yy, x = w0[m], y[m], d["x"][m]
        out[f"A|{g}|n"] = int(m.sum()) if detail else np.nan
        out[f"A|{g}|mean_y"] = np.average(yy, weights=w)
        out[f"A|{g}|sd_x"] = np.sqrt(np.average((x - np.average(x, weights=w)) ** 2, weights=w))
        b = wols(yy, np.column_stack([np.ones(len(x)), x]), w)
        out[f"A|{g}|lpm_obs"] = b[1]
        out[f"A|{g}|lpm_lat"] = b[1] / sR
        bl = wlogit(yy, np.column_stack([np.ones(len(x)), x]), w)
        out[f"A|{g}|logit_obs"] = bl[1]
        out[f"A|{g}|logit_lat_approx"] = bl[1] / sR
        a_d, b_d = deconv_logit(out[f"A|{g}|mean_y"], b[1] / sR)
        out[f"A|{g}|logit_lat_deconv"] = b_d
        out[f"A|{g}|logit_a_deconv"] = a_d
        # quadratic LPM
        bq = wols(yy, np.column_stack([np.ones(len(x)), x, x ** 2]), w)
        out[f"A|{g}|quad_x"] = bq[1]
        out[f"A|{g}|quad_x2"] = bq[2]
        # quintiles of x
        q = wquint_groups(x, w)
        for k in range(5):
            out[f"A|{g}|q{k + 1}"] = np.average(yy[q == k], weights=w[q == k])
            out[f"A|{g}|qx{k + 1}"] = np.average(x[q == k], weights=w[q == k])
    # pooled with a female dummy
    m = base
    X = np.column_stack([np.ones(m.sum()), d["x"][m], f[m].astype(float)])
    b = wols(y[m], X, w0[m])
    out["A|pooled_fe|lpm_obs"] = b[1]
    out["A|pooled_fe|lpm_lat"] = b[1] / sR
    out["A|pooled_fe|logit_lat_approx"] = wlogit(y[m], X, w0[m])[1] / sR

    # ---------------- income / earnings ----------------
    for meas, key, lab in (("E", "logE", "own earnings"), ("W", "logW", "own wage"),
                           ("F", "logF", "family income")):
        v = d[key]
        base = (w0 > 0) & np.isfinite(d["x"]) & np.isfinite(v)
        for g in GROUPS:
            m = base & gm[g]
            w, yy, x, vv = w0[m], y[m], d["x"][m], v[m]
            one = np.ones(len(x))
            p = f"{meas}|{g}|"
            out[p + "n"] = int(m.sum()) if detail else np.nan
            out[p + "mean_y"] = np.average(yy, weights=w)
            out[p + "sd_v"] = np.sqrt(np.average((vv - np.average(vv, weights=w)) ** 2, weights=w))
            out[p + "corr_vx"] = (np.average((vv - np.average(vv, weights=w)) * (x - np.average(x, weights=w)), weights=w)
                                  / (out[p + "sd_v"] * np.sqrt(np.average((x - np.average(x, weights=w)) ** 2, weights=w))))
            # log-level slope
            out[p + "lpm_log"] = wols(yy, np.column_stack([one, vv]), w)[1]
            bc = wols(yy, np.column_stack([one, vv, x]), w)
            out[p + "lpm_log_ctl_obs"] = bc[1]
            out[p + "lpm_afqt_obs_ctl"] = bc[2]
            ce, cx = eiv_two(yy, vv, x, w, R)
            out[p + "lpm_log_ctl_lat"] = ce
            out[p + "lpm_afqt_lat_ctl"] = cx
            # rank (within the group; pooled uses the pooled rank), per 10 percentiles
            r = wrank(vv, w)
            out[p + "lpm_rank10"] = wols(yy, np.column_stack([one, r]), w)[1] / 10.0
            bcr = wols(yy, np.column_stack([one, r, x]), w)
            out[p + "lpm_rank10_ctl_obs"] = bcr[1] / 10.0
            out[p + "lpm_afqt_obs_ctl_rank"] = bcr[2]
            if g == "pooled":
                # within-gender rank
                rw = np.empty(len(vv))
                for gg in (True, False):
                    k = f[m] == gg
                    rw[k] = wrank(vv[k], w[k])
                out[p + "lpm_rank10_wg"] = wols(yy, np.column_stack([one, rw]), w)[1] / 10.0
                out[p + "lpm_rank10_wg_ctl_obs"] = wols(yy, np.column_stack([one, rw, x]), w)[1] / 10.0
                Xf = np.column_stack([one, vv, f[m].astype(float)])
                out[p + "lpm_log_fe"] = wols(yy, Xf, w)[1]
                out[p + "lpm_log_fe_ctl"] = wols(yy, np.column_stack([Xf, x]), w)[1]
                # within-gender rank quintiles for the pooled sample
                qq = np.minimum((rw * 5).astype(int), 4)
                for k in range(5):
                    out[p + f"qwg{k + 1}"] = np.average(yy[qq == k], weights=w[qq == k])
            q = np.minimum((r * 5).astype(int), 4)
            for k in range(5):
                out[p + f"q{k + 1}"] = np.average(yy[q == k], weights=w[q == k])
    return out


# =============================================================================
# Bootstrap
# =============================================================================

def bootstrap(d: dict, B: int) -> dict:
    ids = np.where(np.isfinite(d["hh"]), d["hh"], -np.arange(d["n"]) - 1)
    pop = d["pop"]
    cl = pd.factorize(ids[pop])[0]
    nc = cl.max() + 1
    rng = np.random.default_rng(SEED)
    draws = []
    t0 = time.time()
    for b in range(B):
        c = np.bincount(rng.integers(0, nc, nc), minlength=nc).astype(float)
        mult = np.zeros(d["n"])
        mult[pop] = c[cl]
        draws.append(stats(d, mult))
        if b % 25 == 0:
            LOG.info("  draw %d/%d (%.0fs)", b, B, time.time() - t0)
    keys = draws[0].keys()
    arr = {k: np.array([dr[k] for dr in draws], float) for k in keys}
    return arr


# =============================================================================
# Report
# =============================================================================

def pc(x, d=1): return f"{100 * x:.{d}f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--outdir", type=Path, default=HERE / "estimates")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    B = 20 if a.quick else a.B
    d = load()
    D_ = d
    LOG.info("loaded; Rc=%.4f", d["Rc"])

    # ---- merge / sample accounting ----
    m = d["main"]
    f = d["female"]
    w = d["w"]
    cov = {
        "fertility main sample": {"n": int(m.sum()), "women": int((m & f).sum()), "men": int((m & ~f).sum())},
        "matched to ability panel by CASEID": {"n": int((m & d["matched"]).sum())},
        "with valid afqt_c (4 valid ASVAB IRT z)": {k: int((m & np.isfinite(d["x"]) & g).sum())
                                                     for k, g in (("n", np.ones(d["n"], bool)), ("women", f), ("men", ~f))},
        "afqt_c and own earnings (>=3 yrs, age 30-44)": {k: int((m & np.isfinite(d["x"]) & np.isfinite(d["logE"]) & g).sum())
                                                          for k, g in (("n", np.ones(d["n"], bool)), ("women", f), ("men", ~f))},
        "afqt_c and permanent hourly wage": {k: int((m & np.isfinite(d["x"]) & np.isfinite(d["logW"]) & g).sum())
                                              for k, g in (("n", np.ones(d["n"], bool)), ("women", f), ("men", ~f))},
        "afqt_c and family income": {k: int((m & np.isfinite(d["x"]) & np.isfinite(d["logF"]) & g).sum())
                                      for k, g in (("n", np.ones(d["n"], bool)), ("women", f), ("men", ~f))},
    }
    ok = m & np.isfinite(d["x"])
    sel = {
        "weighted share of main sample with afqt_c": float(w[ok].sum() / w[m].sum()),
        "ever-child share, main sample": float(np.average(d["y"][m], weights=w[m])),
        "ever-child share, afqt sample": float(np.average(d["y"][ok], weights=w[ok])),
        "ever-child share, main sample WITHOUT afqt": float(np.average(d["y"][m & ~np.isfinite(d["x"])], weights=w[m & ~np.isfinite(d["x"])])),
        "n main sample WITHOUT afqt": int((m & ~np.isfinite(d["x"])).sum()),
        "merge rate (matched CASEID / main sample)": float((m & d["matched"]).sum() / m.sum()),
        "mean afqt_c main/afqt sample (weighted)": float(np.average(d["x"][ok], weights=w[ok])),
    }
    ei = d["einfo"]
    inm = m & np.isfinite(d["x"])
    earn_info = {
        "age window": list(EARN_AGE),
        "min valid person-years": MIN_YEARS,
        "median # valid earnings person-years (afqt sample)": float(np.median(ei["nE"][inm])),
        "share of afqt sample with >= 3 earnings years": float((ei["nE"][inm] >= MIN_YEARS).mean()),
        "share with >= 3 yrs but mean real earnings < $1000 (dropped from log)": float(ei["zeroE"][inm].mean()),
        "share with >= 3 family-income years": float((ei["nF"][inm] >= MIN_YEARS).mean()),
        "share with >= 2 wage-rate years": float((ei["nW"][inm] >= MIN_YEARS_WAGE).mean()),
        "dollars": "2010 dollars, BLS CPI-U annual average of year Y-1 (nlsy_ability.CPI_U), pipeline floor $1,000",
    }

    LOG.info("point estimates")
    pt = stats(d, np.ones(d["n"]), detail=True)
    LOG.info("bootstrap B=%d", B)
    draws = bootstrap(d, B)
    se = {k: float(np.nanstd(v, ddof=1)) for k, v in draws.items()}
    res = {k: {"est": float(v), "se": se[k]} for k, v in pt.items()}

    # ---- model mapping (pooled) ----
    p_bar = pt["A|pooled|mean_y"]
    p_tgt = d["p_main"]
    lat = pt["A|pooled|lpm_lat"]
    b_lpm = lat / (p_bar * (1 - p_bar))
    b_lpm_draws = draws["A|pooled|lpm_lat"] / (draws["A|pooled|mean_y"] * (1 - draws["A|pooled|mean_y"]))
    mapping = {
        "p_target": p_tgt, "p_afqt_sample": p_bar,
        "b_from_lpm_over_pq": {"est": float(b_lpm), "se": float(np.std(b_lpm_draws, ddof=1))},
        "b_logit_latent_approx": res["A|pooled|logit_lat_approx"],
        "b_logit_latent_deconv": res["A|pooled|logit_lat_deconv"],
        "a_deconv_at_N01": res["A|pooled|logit_a_deconv"],
        "by_group_b_deconv": {g: res[f"A|{g}|logit_lat_deconv"] for g in GROUPS},
        "earn_implied_ability_slope": {},
    }
    # intercept that holds the population mean at p_tgt for the deconvolved b
    b_d = pt["A|pooled|logit_lat_deconv"]
    a_fix = deconv_logit_a(b_d, p_tgt)
    mapping["a_for_mean_p_target_at_b_deconv"] = float(a_fix)
    for g in GROUPS:
        for var, lab in (("E", "own earnings"), ("W", "own wage")):
            for kk in ("lpm_log", "lpm_log_ctl_lat"):
                e = pt[f"{var}|{g}|{kk}"]
                mapping["earn_implied_ability_slope"][f"{var}|{g}|{kk}"] = {
                    "dP_dlogv": e, "se_dP_dlogv": se[f"{var}|{g}|{kk}"],
                    "implied_per_latent_sd": e * BLOGWAGE, "se": se[f"{var}|{g}|{kk}"] * BLOGWAGE}

    out = dict(generated=dt.datetime.now().isoformat(timespec="seconds"), script="data/spatial/nlsy_fertility_gradient.py",
               B=B, seed=SEED, R_comp_pop=d["Rc"], coverage=cov, selection=sel, earnings_construction=earn_info,
               estimates=res, mapping=mapping,
               note="afqt_c and R_comp_pop held fixed across bootstrap draws; HHID-cluster bootstrap")
    a.outdir.mkdir(parents=True, exist_ok=True)
    (a.outdir / "nlsy_fertility_gradient.json").write_text(json.dumps(out, indent=2, default=float))
    write_md(a.outdir / "nlsy_fertility_gradient.md", d, out, pt, se, draws)
    LOG.info("done")


def write_md(path: Path, d: dict, out: dict, pt: dict, se: dict, draws: dict) -> None:
    R = d["Rc"]
    B = out["B"]
    mp = out["mapping"]
    cov, sel, ei = out["coverage"], out["selection"], out["earnings_construction"]
    def c(k, dig=3, scale=1.0):
        return f"{scale * pt[k]:.{dig}f} ({scale * se[k]:.{dig}f})"
    L = []
    L.append("# NLSY79: gradient of ever having a child in ability and earnings\n")
    L.append(f"Generated by `data/spatial/nlsy_fertility_gradient.py` (B = {B} HHID-cluster bootstrap draws, seed {out['seed']}). "
             "No microdata. SEs in parentheses.\n")
    L.append("## Definitions\n")
    L.append("- Outcome y = 1 if ever had a biological child (`nlsy_fertility.py`): completed-fertility sample = last valid NUMKID observation at age >= 45, "
             "POP_SAMPLES, last-round SAMPWEIGHT/100; both genders pooled and by gender. Mean p = 0.802 (SE 0.005).")
    L.append(f"- x = `afqt_c` from `nlsy_ability.block_a`: equal-weight mean of the four standardised ASVAB IRT z-scores, re-standardised "
             f"(1979 weights) in the Block A population sample (mean 0, SD 1 there). Reliability R_comp_pop = {R:.4f}. Latent slope = observed slope / sqrt(R) "
             "(classical measurement error; same correction as block C). Within the fertility sample SD(x) is "
             f"{pt['A|pooled|sd_x']:.3f} (pooled).")
    L.append("- Logit: weighted logit of y on x; 'latent (approx)' = b_obs / sqrt(R), a first-order approximation (exact only for a linear probability model). "
             "'Latent (deconvolved)' solves E[logistic(a + b x)] = p and E[logistic(a + b x) x] = latent LPM slope for x ~ N(0,1), i.e. it matches the "
             "population-average derivative of the logistic to the measurement-error-corrected LPM slope (no small-b approximation).")
    L.append(f"- Own permanent earnings: mean over person-years at interview ages {EARN_AGE[0]}-{EARN_AGE[1]} (>= {MIN_YEARS} valid years) of wage-and-salary income "
             "in the previous calendar year (`nlsy_ability._income79`: revised-topcode Q13-5), deflated to 2010 dollars with the pipeline's CPI-U (year Y-1); "
             "zeros (non-work years) are kept in the mean; log of the mean, dropped if the mean is < $1,000. Wage-rate variant: mean over years with income >= $1,000 and "
             f">= {MIN_HOURS} annual hours of log(income / hours) (>= {MIN_YEARS_WAGE} years). Family income: same window, TNFI_TRUNC (total net family income, "
             "topcoded), same deflator, log of the mean over >= 3 years. Ranks are weighted mid-rank percentiles within the estimation sample of each row "
             "(pooled rank for the pooled row; 'within-gender' rank is computed separately for women and men in the pooled sample). "
             "Rank slopes are per 10 percentile points.")
    L.append("- Survey interviews are annual to 1994 and biennial after, so a 30-44 window holds up to 15 annual or about 10 mixed person-years; "
             f"median # valid earnings years in the sample = {ei['median # valid earnings person-years (afqt sample)']:.0f}.")
    L.append("- Income in the 30-44 window is largely observed AFTER the first birth for most parents, and family income is jointly determined with marriage and children. "
             "Women's own earnings also respond to children (reverse causality). Own AFQT (ASVAB administered in 1980, ages 15-23) is predetermined and is the cleaner gradient.\n")
    L.append("## Merge and sample\n")
    L.append("| Step | n | women | men |\n|---|---|---|---|")
    for k, v in cov.items():
        L.append(f"| {k} | {v['n']} | {v.get('women', '')} | {v.get('men', '')} |")
    L.append("")
    L.append(f"Merge rate on CASEID: {100 * sel['merge rate (matched CASEID / main sample)']:.1f}% of the fertility main sample. "
             f"Share of the (weighted) main sample with a valid AFQT composite: {100 * sel['weighted share of main sample with afqt_c']:.1f}%. "
             f"Ever-child share: main sample {pc(sel['ever-child share, main sample'])}%, AFQT sample {pc(sel['ever-child share, afqt sample'])}%, "
             f"main sample without AFQT (n = {sel['n main sample WITHOUT afqt']}) {pc(sel['ever-child share, main sample WITHOUT afqt'])}%.\n")

    # ---- ability ----
    L.append("## 1. Ability gradient (y on afqt_c)\n")
    L.append("Slopes are the change in P(ever a child) per 1 SD of AFQT (percentage points for the LPM; log-odds for the logit). "
             "'Latent' divides by sqrt(R_comp_pop).\n")
    L.append("| Group | n | mean y (%) | LPM obs (pp/SD) | LPM latent (pp/SD) | logit obs | logit latent (approx) | logit latent (deconvolved) | intercept a (deconv) |\n|---|---|---|---|---|---|---|---|---|")
    for g in GROUPS:
        p = f"A|{g}|"
        L.append(f"| {g} | {pt[p + 'n']} | {pc(pt[p + 'mean_y'])} | {c(p + 'lpm_obs', 2, 100)} | {c(p + 'lpm_lat', 2, 100)} | "
                 f"{c(p + 'logit_obs', 3)} | {c(p + 'logit_lat_approx', 3)} | {c(p + 'logit_lat_deconv', 3)} | {pt[p + 'logit_a_deconv']:.3f} |")
    L.append(f"| pooled + female dummy | {pt['A|pooled|n']} | | {c('A|pooled_fe|lpm_obs', 2, 100)} | {c('A|pooled_fe|lpm_lat', 2, 100)} | | {c('A|pooled_fe|logit_lat_approx', 3)} | | |")
    L.append("")
    L.append("### Shares with a child by AFQT quintile (percent; weighted; quintiles of afqt_c within the group)\n")
    L.append("| Group | Q1 (low) | Q2 | Q3 | Q4 | Q5 (high) | Q5 - Q1 (pp) | mean x in Q1 / Q5 |\n|---|---|---|---|---|---|---|---|")
    for g in GROUPS:
        p = f"A|{g}|"
        dq = draws[p + "q5"] - draws[p + "q1"]
        L.append(f"| {g} | " + " | ".join(f"{pc(pt[p + f'q{k}'])} ({pc(se[p + f'q{k}'])})" for k in range(1, 6))
                 + f" | {100 * (pt[p + 'q5'] - pt[p + 'q1']):.1f} ({100 * np.std(dq, ddof=1):.1f}) | {pt[p + 'qx1']:.2f} / {pt[p + 'qx5']:.2f} |")
    L.append("")
    L.append("### Quadratic term (LPM of y on x and x^2, observed x)\n")
    L.append("| Group | x (pp/SD) | x^2 (pp/SD^2) | t(x^2) |\n|---|---|---|---|")
    for g in GROUPS:
        p = f"A|{g}|"
        L.append(f"| {g} | {c(p + 'quad_x', 2, 100)} | {c(p + 'quad_x2', 2, 100)} | {pt[p + 'quad_x2'] / se[p + 'quad_x2']:.2f} |")
    L.append("")

    # ---- earnings ----
    names = {"E": "Own permanent earnings (mean real wage/salary income, ages 30-44)",
             "W": "Own permanent hourly wage (mean log, ages 30-44, >= 500 hrs)",
             "F": "Total net family income (mean real TNFI, ages 30-44)"}
    L.append("## 2. Earnings and income gradients\n")
    L.append("All rows use the sample with a valid afqt_c AND the income measure, so the with- and without-AFQT columns share a sample. "
             "'AFQT ctl (obs)' adds observed afqt_c; 'AFQT ctl (latent)' corrects the two-regressor normal equations for AFQT measurement error "
             "(reliability R_comp_pop) and reports the earnings coefficient plus the AFQT coefficient per latent SD. Slopes are in pp of probability.\n")
    for meas in ("E", "W", "F"):
        L.append(f"### {names[meas]}\n")
        L.append("| Group | n | mean y (%) | corr(v, afqt) | per log pt (pp) | per log pt, AFQT ctl obs | per log pt, AFQT ctl latent | AFQT coef, latent ctl (pp/SD) | per 10 pctile | per 10 pctile, AFQT ctl obs |\n|---|---|---|---|---|---|---|---|---|---|")
        for g in GROUPS:
            p = f"{meas}|{g}|"
            L.append(f"| {g} | {pt[p + 'n']} | {pc(pt[p + 'mean_y'])} | {pt[p + 'corr_vx']:.2f} | {c(p + 'lpm_log', 2, 100)} | {c(p + 'lpm_log_ctl_obs', 2, 100)} | "
                     f"{c(p + 'lpm_log_ctl_lat', 2, 100)} | {c(p + 'lpm_afqt_lat_ctl', 2, 100)} | {c(p + 'lpm_rank10', 2, 100)} | {c(p + 'lpm_rank10_ctl_obs', 2, 100)} |")
        p = f"{meas}|pooled|"
        L.append(f"| pooled, within-gender rank | {pt[p + 'n']} | | | | | | | {c(p + 'lpm_rank10_wg', 2, 100)} | {c(p + 'lpm_rank10_wg_ctl_obs', 2, 100)} |")
        L.append(f"| pooled + female dummy | {pt[p + 'n']} | | | {c(p + 'lpm_log_fe', 2, 100)} | {c(p + 'lpm_log_fe_ctl', 2, 100)} (obs ctl) | | | | |")
        L.append("")
        L.append("Shares with a child by quintile of the measure (percent; weighted; Q1 = lowest):\n")
        L.append("| Group | Q1 | Q2 | Q3 | Q4 | Q5 | Q5 - Q1 (pp) |\n|---|---|---|---|---|---|---|")
        for g in GROUPS:
            p = f"{meas}|{g}|"
            dq = draws[p + "q5"] - draws[p + "q1"]
            L.append(f"| {g} | " + " | ".join(f"{pc(pt[p + f'q{k}'])} ({pc(se[p + f'q{k}'])})" for k in range(1, 6))
                     + f" | {100 * (pt[p + 'q5'] - pt[p + 'q1']):.1f} ({100 * np.std(dq, ddof=1):.1f}) |")
        p = f"{meas}|pooled|"
        dq = draws[p + "qwg5"] - draws[p + "qwg1"]
        L.append(f"| pooled, within-gender quintiles | " + " | ".join(f"{pc(pt[p + f'qwg{k}'])} ({pc(se[p + f'qwg{k}'])})" for k in range(1, 6))
                 + f" | {100 * (pt[p + 'qwg5'] - pt[p + 'qwg1']):.1f} ({100 * np.std(dq, ddof=1):.1f}) |")
        L.append("")

    # ---- mapping ----
    L.append("## 3. Mapping to the model\n")
    L.append("Model: P(child | x) = logistic(a + b x), x ~ N(0,1) latent log ability, mean P fixed at p = 0.802.\n")
    bl = mp["b_from_lpm_over_pq"]
    L.append(f"- Pooled latent LPM slope = {c('A|pooled|lpm_lat', 2, 100)} pp per latent SD; sample mean p = {pc(mp['p_afqt_sample'], 2)}%. "
             f"Small-slope map b = slope / (p(1-p)) = **{bl['est']:.3f}** ({bl['se']:.3f}).")
    L.append(f"- Direct weighted logit, latent (approx = b_obs / sqrt(R)): **{c('A|pooled|logit_lat_approx', 3)}**.")
    L.append(f"- Deconvolved logit slope (matches mean p and latent LPM slope exactly for x ~ N(0,1)): **{c('A|pooled|logit_lat_deconv', 3)}**; "
             f"intercept that keeps the population mean at p = 0.802: a = {mp['a_for_mean_p_target_at_b_deconv']:.3f}.")
    L.append("- By group (deconvolved b): " + "; ".join(f"{g} {c(f'A|{g}|logit_lat_deconv', 3)}" for g in GROUPS) + ".")
    L.append("")
    L.append(f"Cross-check from the earnings gradient: (dP/d log earnings) x {BLOGWAGE} (log earnings per latent SD of AFQT, block C target) "
             "= the ability slope that would arise if the entire earnings gradient were driven by ability (pp per latent SD). Compare with the direct latent LPM slope above.\n")
    L.append("| Measure | Group | dP/dlog (pp) | x 0.213 (pp per latent SD) | with AFQT control: dP/dlog (pp) | x 0.213 |\n|---|---|---|---|---|---|")
    for var, nm in (("E", "own earnings"), ("W", "own wage")):
        for g in GROUPS:
            e0 = mp["earn_implied_ability_slope"][f"{var}|{g}|lpm_log"]
            e1 = mp["earn_implied_ability_slope"][f"{var}|{g}|lpm_log_ctl_lat"]
            L.append(f"| {nm} | {g} | {100 * e0['dP_dlogv']:.2f} ({100 * e0['se_dP_dlogv']:.2f}) | {100 * e0['implied_per_latent_sd']:.2f} ({100 * e0['se']:.2f}) | "
                     f"{100 * e1['dP_dlogv']:.2f} ({100 * e1['se_dP_dlogv']:.2f}) | {100 * e1['implied_per_latent_sd']:.2f} ({100 * e1['se']:.2f}) |")
    L.append("")
    L.append("## Caveats\n")
    L.append("- Men under-report fathered children (non-resident fathers, earlier relationships), which plausibly biases men's rates down, "
             "and more so at low AFQT / low earnings, so men's ability slope is likely biased toward positive values (and the earnings slope upward).")
    L.append("- afqt_c is held fixed in the bootstrap; SEs do not include uncertainty in R_comp_pop (SE 0.002) or in the factor fit.")
    L.append("- Earnings and family income are measured after most births and respond to them (especially women's earnings); the earnings gradients are descriptive, not causal.")
    L.append("- The latent logit corrections assume classical measurement error and a normal latent ability distribution.")
    L.append("\n*Written with help from Claude Code.*\n")
    path.write_text("\n".join(L))


if __name__ == "__main__":
    main()
