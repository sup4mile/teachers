#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_ability.py -- the DATA side of the first-stage indirect inference for the
ability block (item T2a of `spatial_calibration.md`, Table 2).

The model (`julia/spatial_model/spatial_continuous.jl`) has log ability
`log z' = rho_z log z + xi` with stationary SD `s_z`, and log wages load on
log z with coefficient `b = 1/(1-eta)`.  T2b will simulate mother-child and
sibling panels from the model, add a score measurement equation and recompute
the auxiliary statistics below with the SAME samples, ages, weights and
composites.  So every statistic here has a crisp definition, restated in the
"Definitions for T2b" section of the report.

    cd data/spatial && .venv/bin/python nlsy_ability.py          # B = 200
    .venv/bin/python nlsy_ability.py --quick                     # B = 20
    .venv/bin/python nlsy_ability.py --selftest                  # unit checks
    .venv/bin/python nlsy_ability.py --bootstrap 500 --seed 1 -v

Writes `estimates/nlsy_ability.json` and `estimates/nlsy_ability.md` (no
microdata).  Needs numpy, pandas, pyarrow only (no scipy / statsmodels).

Inputs (from `fetch_data.py`): raw/nlsy79_ability.parquet, raw/nlscya_ability.parquet,
raw/nlsy97_ability.parquet, with metadata/<cohort>_ability_vars.csv as the
variable dictionary.  Every variable is resolved by (qname, year); every
negative value is an NLS missing code and becomes NaN.

Blocks
  A  NLSY79 population measurement model: one-factor ULS on the four age-normed
     ASVAB IRT z-scores, composite reliability R_comp_pop, KR-20 benchmarks,
     the observed score `afqt_c`, and an own-norming sensitivity.
  B  Mother-child measurement system (NLSY79 women x CNLSY children): two-factor
     ULS with a factor correlation `rho_pc_latent` (the primary moment for
     rho_z), composite reliabilities, two-year child stabilities, nine
     sensitivities.
  S  Siblings: cross-sibling cross-test correlations, `rho_sib_latent`,
     implied rho = rho_sib / rho_pc.
  C  Log-wage slope on latent AFQT at ACS-matched ages (25-34, full-year
     full-time): `wage_slope_latent = slope_obs / sqrt(R_comp_pop)`, the
     primary moment for b * s_z; implied s_z and the gap to chi; 90/10;
     participation checks; sensitivities.
  V  NLSY97 vintage check on block C (CAT-ASVAB, born 1980-84 = the ACS
     2009-13 ages 25-34 population).
  E  Earnings checks: permanent-earnings correlations, mother-child and
     sibling, with the implied rho.  Checks, not targets.

Child scores: PIAT comprehension is IMPUTED for weak readers.  NLS gives children with
RECOG < 19 no comprehension test and copies the recognition score into COMP (~15% of
in-window assessments, mostly ages 5-7).  Those rows stay in the COMP norming reference
but are excluded from child-level COMP means and from COMP two-year stabilities.  Child
age is CSAGE (child supplement), falling back to MSAGE where missing.  The child
composite reliability accounts for which tests each child has.

Inference: cluster bootstrap over NLSY79 households (`HHID`; CNLSY children
inherit their mother's), implemented as multinomial household multiplicities
applied to the weights -- exactly equivalent to resampling households with
replacement and treating each draw of a household as a distinct copy (pairs
form only within a copy).  Every estimate is re-run on every draw except the
child and NLSY97 age-norming, which are held fixed.  NLSY97: cluster = PUBID
(the extract has no household id, so SEs ignore NLSY97 sibling clustering).
SE = bootstrap SD, CI = 95% percentile interval.

Written with help from Claude Code.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from statistics import NormalDist
from typing import Any, Sequence

import numpy as np
import pandas as pd

LOG = logging.getLogger("nlsy_ability")
HERE = Path(__file__).resolve().parent

# =============================================================================
# Constants
# =============================================================================

ETA = 0.0807                    # goods-investment elasticity (Table 1)
CHI = float(np.log(1.12))       # earnings-achievement bridge, 0.11333
AGE_LO, AGE_HI = 25, 34         # ACS-matched age window
MIN_HGC = 9                     # "at least one year of high school"
FT_HOURS = 35                   # full-time: usual hours per week
FY_WEEKS = 50                   # full-year: weeks per year
HP_HOURS = 15 * 52              # note's 15-hour home-production cutoff, per year
TRIM = (0.01, 0.99)             # weighted log-wage trimming within year

TESTS79 = ["AR", "WK", "PC", "MK"]
IRT_Q = {"AR": "ASVAB-ARITHREASON-IRT-ZSCORE", "WK": "ASVAB-WORDKNOW-IRT-ZSCORE",
         "PC": "ASVAB-PARACOMP-IRT-ZSCORE", "MK": "ASVAB-MATHKNOW-IRT-ZSCORE"}
SEC_Q = {"AR": "ASVAB-34", "WK": "ASVAB-35", "PC": "ASVAB-36", "MK": "ASVAB-40"}
ITEM_Q = {"AR": ("ASVAB-ARITHMETIC-REASONING", 30), "WK": ("ASVAB-WORD-KNOWLEDGE", 35),
          "PC": ("ASVAB-PARAGRAPH-COMPREHENSION", 15),
          "MK": ("ASVAB-MATHEMATICS-KNOWLEDGE", 25)}
CH_TESTS = ["MATH", "RECOG", "COMP", "PPVT"]
CH_ROUNDS = list(range(1986, 2017, 2))
CH_WINDOW = (60, 179)           # child age in months at assessment: ages 5-14
CH_WINDOW_OLD = (120, 179)      # ages 10-14 only (sensitivity f)
MOTHER_SAMPLES = (5, 6, 7, 8, 13, 14)       # women in the NLSY79
POP_SAMPLES = (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14)
XS_SAMPLES = (1, 2, 3, 4, 5, 6, 7, 8)
GROUPS_PC = [0, 0, 0, 0, 1, 1, 1, 1]
DROP_DOMAIN = ((0, 3), (1, 2), (5, 6))      # AR-MK, WK-PC, RECOG-COMP
IND8 = TESTS79 + CH_TESTS

_NORM = NormalDist()
_PPF = np.frompyfunc(_NORM.inv_cdf, 1, 1)


# =============================================================================
# The estimate registry -- one record per number that leaves this script
# =============================================================================

@dataclass
class Estimate:
    key: str
    label: str
    model_object: str
    value: float | None
    se: float | None
    ci_lo: float | None
    ci_hi: float | None
    n: int | None
    units: str
    source: str
    note: str
    extra: dict[str, Any] = field(default_factory=dict)


class Registry:
    def __init__(self) -> None:
        self.items: list[Estimate] = []

    def add(self, est: Estimate) -> Estimate:
        self.items.append(est)
        return est

    def to_dict(self) -> list[dict[str, Any]]:
        return [asdict(e) for e in self.items]

    def get(self, key: str) -> Estimate:
        return next(e for e in self.items if e.key == key)


# =============================================================================
# Cohort access: (qname, year) -> column, negatives -> NaN
# =============================================================================

class Cohort:
    """One raw table plus its variable dictionary.

    `get(qname, year)` returns a float array with EVERY negative value (the NLS
    missing codes -1..-7) replaced by NaN.  R-numbers never appear in code.
    """

    def __init__(self, name: str, root: Path) -> None:
        self.name = name
        self.df = pd.read_parquet(root / "raw" / f"{name}_ability.parquet")
        self.meta = pd.read_csv(root / "metadata" / f"{name}_ability_vars.csv",
                                dtype=str)
        self.map: dict[tuple[str, str], str] = {}
        self.by_q: dict[str, list[str]] = {}
        for r in self.meta.itertuples():
            self.map[(r.qname, str(r.year))] = r.rnum
            self.by_q.setdefault(r.qname, []).append(str(r.year))
        LOG.info("  loaded %-8s %6d rows x %d cols", name, *self.df.shape)

    def has(self, q: str, y: str | int | None = None) -> bool:
        if y is None:
            return len(self.by_q.get(q, [])) == 1
        return (q, str(y)) in self.map

    def get(self, q: str, y: str | int | None = None, scale: float = 1.0,
            shift: float = 0.0) -> np.ndarray:
        if y is None:
            ys = self.by_q[q]
            if len(ys) != 1:
                raise KeyError(f"{q}: year needed, has {ys}")
            y = ys[0]
        rn = self.map[(q, str(y))]
        x = self.df[rn].to_numpy(dtype=float)
        x = np.where(x < 0, np.nan, x)
        if LOG.isEnabledFor(logging.DEBUG):
            v = x[np.isfinite(x)]
            LOG.debug("    %-8s %-32s %-6s n=%6d min=%g max=%g", self.name, q, y,
                      v.size, v.min() if v.size else np.nan,
                      v.max() if v.size else np.nan)
        return x / scale + shift

    def try_get(self, q: str, y: str | int | None = None, **kw) -> np.ndarray | None:
        try:
            return self.get(q, y, **kw)
        except KeyError:
            return None


# =============================================================================
# Weighted statistics, rank-normalisation, ULS factor model
# =============================================================================

def wmean(x: np.ndarray, w: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    return float((x[m] * w[m]).sum() / w[m].sum()) if m.any() else np.nan


def wvar(x: np.ndarray, w: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    mu = (x[m] * w[m]).sum() / w[m].sum()
    return float((w[m] * (x[m] - mu) ** 2).sum() / w[m].sum())


def wcov(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    ww = w[m] / w[m].sum()
    return float((ww * (x[m] - (ww * x[m]).sum()) * (y[m] - (ww * y[m]).sum())).sum())


def wcorr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Weighted correlation on the rows where both are present (the SAME rows
    for both variances)."""
    m = np.isfinite(x) & np.isfinite(y)
    return _wcorr_m(x[m], y[m], w[m])


def _wcorr_m(x, y, w) -> float:
    m = np.isfinite(w) & (w > 0)
    if m.sum() < 3:
        return np.nan
    x, y, w = x[m], y[m], w[m] / w[m].sum()
    xc, yc = x - (w * x).sum(), y - (w * y).sum()
    d = np.sqrt((w * xc ** 2).sum() * (w * yc ** 2).sum())
    return float((w * xc * yc).sum() / d) if d > 0 else np.nan


def wquantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    xs, ws = x[m], w[m]
    o = np.argsort(xs, kind="stable")
    xs, ws = xs[o], ws[o]
    c = np.cumsum(ws) - 0.5 * ws
    return float(np.interp(q, c / ws.sum(), xs))


def wcorr_cross(X: np.ndarray, Y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """K x L matrix of weighted, pairwise-complete correlations between the
    columns of X and of Y.  Each (k, l) entry uses only rows where BOTH x_k and
    y_l are present, and the same rows for the two means and variances."""
    w = np.where(np.isfinite(w), w, 0.0)[:, None]
    MX, MY = np.isfinite(X), np.isfinite(Y)
    X0, Y0 = np.where(MX, X, 0.0), np.where(MY, Y, 0.0)
    W = (w * MX).T @ MY.astype(float)
    Sx, Sy = (w * X0).T @ MY, (w * MX).T @ Y0
    Sxx, Syy = (w * X0 ** 2).T @ MY, (w * MX).T @ Y0 ** 2
    Sxy = (w * X0).T @ Y0
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = Sx / W, Sy / W
        vx, vy = Sxx / W - mx ** 2, Syy / W - my ** 2
        return (Sxy / W - mx * my) / np.sqrt(vx * vy)


def phi_inv(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    return np.asarray(_PPF(p), dtype=float) if p.size else p.copy()


def rank_normal(x: np.ndarray, w: np.ndarray, cell: np.ndarray,
                min_ref: int = 30) -> np.ndarray:
    """Weighted within-cell rank -> N(0,1).

    Reference = rows in the cell with valid x and w > 0.  EVERY row with valid x
    (reference or not) gets p = (W_below + W_tied/2) / W_total against its
    cell's reference, then z = Phi^{-1}(p).  Non-reference rows beyond the
    reference range are clipped to p in [1e-4, 1 - 1e-4].
    """
    x = np.asarray(x, dtype=float)
    out = np.full(x.shape, np.nan)
    valid = np.isfinite(x)
    ref = valid & np.isfinite(w) & (w > 0)
    small = 0
    for c in np.unique(cell[valid]):
        inc = valid & (cell == c)
        r = ref & (cell == c)
        if r.sum() < min_ref:
            small += 1
        if r.sum() == 0:
            continue
        o = np.argsort(x[r], kind="stable")
        xs = x[r][o]
        cw = np.concatenate([[0.0], np.cumsum(w[r][o])])
        lo = np.searchsorted(xs, x[inc], side="left")
        hi = np.searchsorted(xs, x[inc], side="right")
        p = 0.5 * (cw[lo] + cw[hi]) / cw[-1]
        out[inc] = phi_inv(np.clip(p, 1e-4, 1 - 1e-4))
    if small:
        LOG.debug("rank_normal: %d cell(s) with < %d reference rows", small, min_ref)
    return out


def linear_norm(x: np.ndarray, w: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Weighted (x - mean_cell) / sd_cell (sensitivity to rank-normalisation)."""
    out = np.full(np.shape(x), np.nan)
    ok = np.isfinite(x)
    for c in np.unique(cell[ok]):
        inc = ok & (cell == c)
        r = inc & np.isfinite(w) & (w > 0)
        if r.sum() < 3:
            continue
        mu, sd = wmean(x[r], w[r]), np.sqrt(wvar(x[r], w[r]))
        if sd > 0:
            out[inc] = (x[inc] - mu) / sd
    return out


@dataclass
class FactorFit:
    lam: np.ndarray
    rho: float | None
    comm: np.ndarray
    rmsr: float
    converged: bool
    heywood: bool
    n_iter: int


def fit_factor_uls(R: np.ndarray, groups: Sequence[int],
                   drop_pairs: Sequence[tuple[int, int]] = ()) -> FactorFit:
    """ULS fit of a one- or two-factor model to the OFF-DIAGONAL of a correlation
    matrix.  Standardised indicators, unit factor variances.  Sigma_ij =
    lam_i lam_j inside a group and lam_i lam_j rho across the two groups.
    Levenberg-Marquardt with the analytic Jacobian.  Heywood cases (|lam| > 1)
    are flagged, not clipped."""
    R = np.asarray(R, dtype=float)
    g = np.asarray(groups)
    p = len(g)
    two = len(set(g.tolist())) > 1
    drop = {tuple(sorted(d)) for d in drop_pairs}
    pr = [(i, j) for i in range(p) for j in range(i + 1, p)
          if (i, j) not in drop and np.isfinite(R[i, j])]
    I = np.array([a for a, _ in pr])
    J = np.array([b for _, b in pr])
    cross = (g[I] != g[J]).astype(float)
    m = len(I)
    ar = np.arange(m)

    # start values: triads inside groups, rho from the mean cross correlation
    lam0 = np.full(p, 0.7)
    for i in range(p):
        vals = []
        for j in range(p):
            for k in range(j + 1, p):
                if i in (j, k) or g[j] != g[i] or g[k] != g[i]:
                    continue
                if np.isfinite(R[i, j] * R[i, k] * R[j, k]) and R[j, k] > 0.02:
                    vals.append(R[i, j] * R[i, k] / R[j, k])
        if vals:
            lam0[i] = np.sqrt(np.clip(np.mean(vals), 0.01, 0.9))
    th = lam0.copy()
    if two:
        cr = R[np.ix_(g == g.min(), g != g.min())]
        cr = cr[np.isfinite(cr)]
        rho0 = cr.mean() / (lam0[g == g.min()].mean() * lam0[g != g.min()].mean())
        th = np.append(th, np.clip(rho0, -0.95, 0.95))

    def resid(t):
        r = np.where(cross > 0, t[p] if two else 1.0, 1.0)
        return R[I, J] - t[I] * t[J] * r

    def jac(t):
        r = np.where(cross > 0, t[p] if two else 1.0, 1.0)
        Jm = np.zeros((m, len(t)))
        Jm[ar, I] = -t[J] * r
        Jm[ar, J] = -t[I] * r
        if two:
            Jm[:, p] = -t[I] * t[J] * cross
        return Jm

    r0 = resid(th)
    cost = float(r0 @ r0)
    mu, ok, it = 1e-3, False, 0
    for it in range(1, 501):
        Jm = jac(th)
        A = Jm.T @ Jm
        gr = Jm.T @ r0
        step = np.linalg.solve(A + mu * np.diag(np.diag(A) + 1e-12), -gr)
        r1 = resid(th + step)
        c1 = float(r1 @ r1)
        if c1 < cost:
            th, r0, dc, cost = th + step, r1, cost - c1, c1
            mu = max(mu / 3, 1e-12)
            if np.abs(step).max() < 1e-10 or dc < 1e-16:
                ok = True
                break
        else:
            mu *= 4
            if mu > 1e12:
                break
    # converged only if the final gradient J'r vanishes, whichever rule stopped the loop
    ok = bool(np.abs(jac(th).T @ r0).max() < 1e-7)
    lam = th[:p]
    rho = float(th[p]) if two else None
    return FactorFit(lam=lam, rho=rho, comm=lam ** 2, rmsr=float(np.sqrt(cost / m)),
                     converged=ok, heywood=bool(np.any(np.abs(lam) > 1)), n_iter=it)


def comp_reliability(lam: np.ndarray, S: np.ndarray) -> tuple[float, float]:
    """Equal-weight composite of standardised indicators.
    R_comp = (sum lam)^2 / (1' S 1) with S the SAMPLE correlation matrix (share
    of the actual composite's variance that is common factor);
    omega = (sum lam)^2 / ((sum lam)^2 + sum(1 - lam^2))."""
    num = float(lam.sum() ** 2)
    S = np.where(np.isfinite(S), S, 0.0)
    return num / float(S.sum()), num / (num + float((1 - lam ** 2).sum()))


def comp_reliability_avail(lam: np.ndarray, S: np.ndarray, X: np.ndarray,
                           w: np.ndarray) -> float:
    """Reliability of the composite actually used: the mean of the standardised
    tests each child has.  R = sum_i w_i (mean_{k in K_i} lam_k)^2 /
    sum_i w_i (1' S_{K_i} 1 / |K_i|^2), K_i = child i's available tests."""
    M = np.isfinite(X) & (w > 0)[:, None]
    keep = M.any(axis=1)
    pats, inv = np.unique(M[keep], axis=0, return_inverse=True)
    wk = w[keep]
    Wp = np.bincount(inv.ravel(), weights=wk, minlength=len(pats))
    num = den = 0.0
    S0 = np.where(np.isfinite(S), S, 0.0)
    for pat, wp in zip(pats, Wp):
        K = np.flatnonzero(pat)
        num += wp * lam[K].mean() ** 2
        den += wp * S0[np.ix_(K, K)].sum() / len(K) ** 2
    return num / den


def wls(y: np.ndarray, X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted least squares via the Gram matrix; lstsq handles the all-zero
    dummy columns a bootstrap draw can produce (minimum-norm solution)."""
    Xw = X * w[:, None]
    b, *_ = np.linalg.lstsq(Xw.T @ X, Xw.T @ y, rcond=None)
    return b


def dummies(codes: np.ndarray, levels: Sequence[int], drop_first: bool = True
            ) -> np.ndarray:
    lv = list(levels)[1:] if drop_first else list(levels)
    return (codes[:, None] == np.array(lv)[None, :]).astype(float)


def trim_mask(x: np.ndarray, w: np.ndarray, grp: np.ndarray,
              lohi: tuple[float, float] = TRIM) -> np.ndarray:
    """Keep rows inside the weighted [lo, hi] quantiles of x within each group."""
    keep = np.zeros(len(x), dtype=bool)
    for gval in np.unique(grp):
        i = np.flatnonzero(grp == gval)
        lo, hi = wquantile(x[i], w[i], lohi[0]), wquantile(x[i], w[i], lohi[1])
        keep[i] = (x[i] >= lo) & (x[i] <= hi)
    return keep


# =============================================================================
# Data preparation
# =============================================================================

@dataclass
class PY79:
    """NLSY79 person-years at ages 25-34 (all survey rounds)."""
    p: np.ndarray            # row in NLSY79 person table
    Y: np.ndarray            # survey year
    age: np.ndarray
    hgc: np.ndarray
    esr: np.ndarray
    wks: np.ndarray
    hrs: np.ndarray
    inc: np.ndarray
    hrp1: np.ndarray         # dollars per hour, current/most recent job
    lw: np.ndarray = None    # log hourly wage from income / hours
    fy: np.ndarray = None    # full-year full-time wage sample (before trimming)
    elig: np.ndarray = None  # C5 sample: hgc >= 9, not armed forces
    fyft: np.ndarray = None  # C5 outcome 1{FYFT wage worker}, NaN if unknown
    hp: np.ndarray = None    # C5 outcome 1{hours < 15*52}, NaN if unknown
    age_d: np.ndarray = None
    yr_d: np.ndarray = None


@dataclass
class N79:
    n: int
    caseid: np.ndarray
    hh_idx: np.ndarray       # household index 0..H-1
    nhh: int
    sample: np.ndarray
    female: np.ndarray
    w: np.ndarray            # 1979 SAMPWEIGHT / 100
    birth_cell: np.ndarray   # birth year x 4-month group
    z: np.ndarray            # n x 4 IRT z-scores (AR, WK, PC, MK)
    sec: np.ndarray          # n x 4 section standard scores (own-norming)
    items: dict[str, np.ndarray]
    asamp: np.ndarray        # block A sample mask
    msamp: np.ndarray        # block B mother sample mask
    py: PY79 = None
    kids: "Kids" = None
    ya: "YA" = None
    sets: dict[str, "DyadSet"] = field(default_factory=dict)
    stab: dict[str, Any] = field(default_factory=dict)
    years: list[int] = field(default_factory=list)


def _income79(c: Cohort, Y: int) -> np.ndarray:
    """Wage and salary income past calendar year: revised topcode where it exists
    (1982-2000), else the 1979-81 variable, else the 2002+ truncated one."""
    for q in ("Q13-5_TRUNC_REVISED", "Q13-5", "Q13-5_TRUNC"):
        if c.has(q, Y):
            return c.get(q, Y)
    return np.full(c.df.shape[0], np.nan)


def build_n79(c: Cohort) -> N79:
    n = c.df.shape[0]
    caseid = c.get("CASEID", 1979)
    hh = c.get("HHID", 1979)
    _, hh_idx = np.unique(hh, return_inverse=True)
    sample = c.get("SAMPLE_ID", 1979).astype(int)
    female = (c.get("SAMPLE_SEX", 1979) == 2)
    w = c.get("SAMPWEIGHT", 1979, scale=100.0)
    by = np.where(np.isfinite(c.get("Q1-3_A~Y", 1981)), c.get("Q1-3_A~Y", 1981),
                  c.get("Q1-3_A~Y", 1979))
    bm = np.where(np.isfinite(c.get("Q1-3_A~M", 1981)), c.get("Q1-3_A~M", 1981),
                  c.get("Q1-3_A~M", 1979))
    cell = np.where(np.isfinite(by) & np.isfinite(bm),
                    by * 10 + (np.nan_to_num(bm, nan=1) - 1) // 4, -1).astype(int)
    z = np.column_stack([c.get(IRT_Q[t], "XRND", scale=100.0, shift=-5.0)
                         for t in TESTS79])
    sec = np.column_stack([c.get(SEC_Q[t], 1981) for t in TESTS79])
    items = {}
    for t, (stem, k) in ITEM_Q.items():
        items[t] = np.column_stack([c.get(f"{stem}-{i}", "XRND") for i in range(1, k + 1)])
    zok = np.isfinite(z).all(axis=1)
    asamp = np.isin(sample, POP_SAMPLES) & zok & np.isfinite(w)
    msamp = np.isin(sample, MOTHER_SAMPLES) & zok & np.isfinite(w)
    D = N79(n, caseid, hh_idx, int(hh_idx.max() + 1), sample, female, w, cell, z,
            sec, items, asamp, msamp)
    D.py = build_py79(c, D)
    return D


def build_py79(c: Cohort, D: N79) -> PY79:
    cols: dict[str, list[np.ndarray]] = {k: [] for k in
        ("p", "Y", "age", "hgc", "esr", "wks", "hrs", "inc", "hrp1")}
    hgc_ever = c.get("HGC_EVER", "XRND")
    years = sorted(int(y) for y in c.by_q["AGEATINT"] if y.isdigit())
    for Y in years:
        age = c.get("AGEATINT", Y)
        keep = np.flatnonzero((age >= AGE_LO) & (age <= AGE_HI))
        if keep.size == 0:
            continue
        hg = c.try_get(f"HGCREV{str(Y)[-2:]}", Y)
        hg = hgc_ever if hg is None else np.where(np.isfinite(hg), hg, hgc_ever)
        esr = c.try_get("ESR_COL", Y)
        esr = np.full(D.n, np.nan) if esr is None else esr
        hrp = c.try_get("HRP1", Y)
        hrp = np.full(D.n, np.nan) if hrp is None else hrp / 100.0
        vals = dict(p=keep, Y=np.full(keep.size, Y), age=age[keep], hgc=hg[keep],
                    esr=esr[keep], wks=c.get("WKSWK-PCY", Y)[keep],
                    hrs=c.get("HRSWK-PCY", Y)[keep], inc=_income79(c, Y)[keep],
                    hrp1=hrp[keep])
        for k, v in vals.items():
            cols[k].append(v)
    a = {k: np.concatenate(v) for k, v in cols.items()}
    py = PY79(**a)
    with np.errstate(divide="ignore", invalid="ignore"):
        hpw = py.hrs / py.wks
        py.fy = ((py.inc > 0) & (py.wks >= FY_WEEKS) & (hpw >= FT_HOURS)
                 & (py.hrs > 0) & ~(py.esr == 4) & (py.hgc >= MIN_HGC))
        py.lw = np.where(py.fy, np.log(py.inc / py.hrs), np.nan)
        py.elig = (py.hgc >= MIN_HGC) & ~(py.esr == 4)
        # outcome known-zero if a defining condition demonstrably fails
        fails = (py.inc == 0) | (py.wks < FY_WEEKS) | (hpw < FT_HOURS)
        known = np.isfinite(py.inc) & np.isfinite(py.wks) & np.isfinite(py.hrs)
        py.fyft = np.where(py.fy, 1.0, np.where(fails | (known & ~py.fy), 0.0, np.nan))
        py.fyft = np.where(py.fy | fails | known, py.fyft, np.nan)
        py.hp = np.where(np.isfinite(py.hrs), (py.hrs < HP_HOURS).astype(float), np.nan)
    D.years = sorted(set(py.Y.tolist()))
    py.age_d = dummies(py.age.astype(int), range(AGE_LO, AGE_HI + 1))
    py.yr_d = dummies(py.Y.astype(int), D.years)
    return py


# ---- CNLSY children ---------------------------------------------------------

@dataclass
class Kids:
    n: int
    cpubid: np.ndarray
    mrow: np.ndarray         # row of the mother in the NLSY79 table (-1 if none)
    female: np.ndarray
    cyrb: np.ndarray
    bthordr: np.ndarray
    long: pd.DataFrame       # one row per (child, round, test) assessment
    scores: dict[str, np.ndarray]   # variant -> n x 4 child-level scores
    counts: dict[str, np.ndarray]
    cw: dict[str, np.ndarray]       # variant -> mean child weight
    norm_info: dict[str, Any] = field(default_factory=dict)
    nrounds: np.ndarray = None      # distinct assessment rounds in the main window


@dataclass
class YA:
    krow: np.ndarray
    Y: np.ndarray
    age: np.ndarray
    inc: np.ndarray
    w: np.ndarray            # YA round weight (fallback: mother's weight)
    age_d: np.ndarray = None
    yr_d: np.ndarray = None


def build_kids(cy: Cohort, D: N79) -> Kids:
    n = cy.df.shape[0]
    cpub = cy.get("CPUBID")
    mp = cy.get("MPUBID")
    mrow = pd.Index(D.caseid).get_indexer(mp)
    female = cy.get("CSEX") == 2
    cyrb, bth = cy.get("CYRB"), cy.get("BTHORDR")
    frames = []
    n_fallback = n_valid = 0
    for yr in CH_ROUNDS:
        s = str(yr)
        if not cy.has(f"MATH{yr}", s):
            continue
        # age in months: child supplement (CSAGE), MSAGE only where CSAGE is missing
        age_m = cy.get(f"MSAGE{yr}", s)
        age_c = cy.try_get(f"CSAGE{yr}", s)
        age = age_m if age_c is None else np.where(np.isfinite(age_c), age_c, age_m)
        recog = cy.get(f"RECOG{yr}", s)
        wt = cy.try_get(f"CSAMWT{yr}_REV", s)
        wt0 = cy.get(f"CSAMWT{yr}", s)
        wt = wt0 if wt is None else wt
        wt = np.where(wt > 0, wt, np.nan) / 100.0
        for ti, t in enumerate(CH_TESTS):
            raw = cy.get(f"{t}{yr}", s)
            ok = np.isfinite(raw) & np.isfinite(age)
            idx = np.flatnonzero(ok)
            # NLS gives children with RECOG < 19 no comprehension test and copies
            # the recognition score into COMP: flag those rows as imputed
            imp = (t == "COMP") & (recog[idx] < 19)
            n_valid += idx.size
            if age_c is not None:
                n_fallback += int((~np.isfinite(age_c[idx])).sum())
            frames.append(pd.DataFrame({"child": idx, "round": yr, "test": ti,
                                        "raw": raw[idx], "age": age[idx],
                                        "w": wt[idx], "imp": imp}))
    long = pd.concat(frames, ignore_index=True)
    long = long[(long.age >= CH_WINDOW[0]) & (long.age <= CH_WINDOW[1])].reset_index(drop=True)
    n_imp = int(long.imp.sum())
    LOG.info("  age from MSAGE where CSAGE missing: %d of %d assessments; imputed COMP "
             "(RECOG < 19): %d of %d in-window assessments (%.1f%%)", n_fallback, n_valid,
             n_imp, len(long), 100 * n_imp / len(long))
    long["cell"] = (long.age // 3).astype(int)
    long["z_rank"] = np.nan
    long["z_lin"] = np.nan
    small = 0
    for ti in range(4):
        m = (long.test == ti).to_numpy()
        x, w, cell = long.raw.to_numpy()[m], long.w.to_numpy()[m], long.cell.to_numpy()[m]
        long.loc[m, "z_rank"] = rank_normal(x, w, cell)
        long.loc[m, "z_lin"] = linear_norm(x, w, cell)
        for cc in np.unique(cell):
            small += int(((cell == cc) & np.isfinite(w) & (w > 0)).sum() < 30)
    LOG.info("  child assessments in window: %d; test x 3-month cells with < 30 "
             "reference rows: %d", len(long), small)
    K = Kids(n, cpub, mrow, female, cyrb, bth, long, {}, {}, {})
    K.nrounds = np.bincount(long.drop_duplicates(["child", "round"]).child, minlength=n)
    K.norm_info = {"n_assess": int(len(long)), "cells_lt30": small, "n_imputed_comp": n_imp,
                   "n_age_fallback_msage": n_fallback, "n_comp": int((long.test == 2).sum()), "n_assess_all_ages": n_valid}
    for var, zc, lo in (("main", "z_rank", CH_WINDOW[0]), ("w10", "z_rank", CH_WINDOW_OLD[0]),
                        ("linear", "z_lin", CH_WINDOW[0])):
        sub = long[(long.age >= lo) & ~long.imp]     # imputed COMP rows never enter child means
        S = np.full((n, 4), np.nan)
        C = np.zeros((n, 4))
        for ti in range(4):
            s = sub[(sub.test == ti) & sub[zc].notna()]
            sm = np.bincount(s.child, weights=s[zc], minlength=n)
            ct = np.bincount(s.child, minlength=n)
            S[ct > 0, ti] = sm[ct > 0] / ct[ct > 0]
            C[:, ti] = ct
        rd = sub.drop_duplicates(["child", "round"])
        rd = rd[rd.w.notna()]
        cwt = np.full(n, np.nan)
        cnt = np.bincount(rd.child, minlength=n)
        sm = np.bincount(rd.child, weights=rd.w, minlength=n)
        cwt[cnt > 0] = sm[cnt > 0] / cnt[cnt > 0]
        K.scores[var], K.counts[var], K.cw[var] = S, C, cwt
    return K


def build_stability(K: Kids) -> dict[str, tuple[np.ndarray, ...]]:
    """Two-year stability pairs: normed z at round y and y+2, same child, both in
    the window.  Weight = round-y child weight."""
    out = {}
    L = K.long[K.long.z_rank.notna() & ~K.long.imp]
    for ti, t in enumerate(CH_TESTS):
        a = L[L.test == ti][["child", "round", "z_rank", "w"]]
        b = a[["child", "round", "z_rank"]].copy()
        b["round"] -= 2
        mg = a.merge(b, on=["child", "round"], suffixes=("1", "2"))
        mg = mg[mg.w.notna()]
        out[t] = (mg.child.to_numpy(), mg.z_rank1.to_numpy(), mg.z_rank2.to_numpy(),
                  mg.w.to_numpy())
    return out


def build_ya(cy: Cohort, D: N79, K: Kids) -> YA:
    rows = []
    ages_names = {}
    for yr in range(1994, 2021, 2):
        yy = str(yr)[-2:]
        s = str(yr)
        an = f"AGEINT{yy}" if cy.has(f"AGEINT{yy}", s) else f"AGEINT{yr}"
        if not cy.has(an, s):
            continue
        age = cy.get(an, s)
        inc = cy.get("Q15-5", s) if cy.has("Q15-5", s) else cy.try_get("Q15-5-TOP", s)
        if inc is None:
            continue
        wq = f"YA{yy}WEIGHT_REVISED" if cy.has(f"YA{yy}WEIGHT_REVISED", s) else f"YA{yy}WEIGHT"
        wt = cy.get(wq, s) if cy.has(wq, s) else np.full(K.n, np.nan)
        wt = np.where(wt > 0, wt, np.nan) / 100.0
        idx = np.flatnonzero((age >= AGE_LO) & (age <= AGE_HI) & (inc > 0) & (K.mrow >= 0))
        rows.append((idx, np.full(idx.size, yr), age[idx], inc[idx], wt[idx]))
    krow, Y, age, inc, wt = (np.concatenate([r[i] for r in rows]) for i in range(5))
    wt = np.where(np.isfinite(wt), wt, D.w[K.mrow[krow]])       # fallback: mother's weight
    ya = YA(krow, Y, age, inc, wt)
    ya.age_d = dummies(age.astype(int), range(AGE_LO, AGE_HI + 1))
    ya.yr_d = dummies(Y.astype(int), sorted(set(Y.tolist())))
    return ya


# ---- dyads -------------------------------------------------------------------

@dataclass
class DyadSet:
    name: str
    krow: np.ndarray
    mrow: np.ndarray
    xm: np.ndarray           # n x 4 mother IRT z
    xc: np.ndarray           # n x 4 child scores
    w0: np.ndarray           # base weight (before bootstrap multiplicity)
    pI: np.ndarray = None    # ordered sibling pairs (indices into the dyad arrays)
    pJ: np.ndarray = None
    pw0: np.ndarray = None
    pm: np.ndarray = None    # mother row of each pair
    n_mothers: int = 0
    n_children: int = 0
    n_families: int = 0
    n_pairs: int = 0


def build_dyadset(D: N79, name: str, var: str = "main", mother_mask: np.ndarray | None = None,
                  w_mode: str = "mother", kid_mask: np.ndarray | None = None,
                  first_only: bool = False) -> DyadSet:
    K = D.kids
    S = K.scores[var]
    mm = D.msamp if mother_mask is None else mother_mask
    ok = (K.mrow >= 0) & np.isfinite(S).any(axis=1)
    ok &= np.where(K.mrow >= 0, mm[np.maximum(K.mrow, 0)], False)
    if kid_mask is not None:
        ok &= kid_mask
    if w_mode == "child":
        ok &= np.isfinite(K.cw[var])
    idx = np.flatnonzero(ok)
    if first_only:
        o = np.lexsort((idx, K.bthordr[idx], K.mrow[idx]))
        idx = idx[o]
        first = np.concatenate([[True], K.mrow[idx][1:] != K.mrow[idx][:-1]])
        idx = idx[first]
    mrow = K.mrow[idx]
    nc = np.bincount(mrow, minlength=D.n)[mrow]
    if w_mode == "mother":
        w0 = D.w[mrow] / nc
    elif w_mode == "mother_only":
        w0 = D.w[mrow].copy()
    elif w_mode == "unit":
        w0 = np.ones(len(idx))
    elif w_mode == "child":
        w0 = K.cw[var][idx]
    else:
        raise ValueError(w_mode)
    ds = DyadSet(name, idx, mrow, D.z[mrow], S[idx], w0)
    ds.n_mothers, ds.n_children = int(np.unique(mrow).size), int(idx.size)
    # ordered sibling pairs within a mother; each family's total weight = mother's weight
    o = np.argsort(mrow, kind="stable")
    ms = mrow[o]
    starts = np.flatnonzero(np.concatenate([[True], ms[1:] != ms[:-1]]))
    ends = np.append(starts[1:], len(ms))
    I, J, PW, PM = [], [], [], []
    for s, e in zip(starts, ends):
        k = e - s
        if k < 2:
            continue
        g = o[s:e]
        a, b = np.triu_indices(k, 1)
        npair = a.size
        wf = D.w[ms[s]] / (2 * npair)
        I += [g[a], g[b]]
        J += [g[b], g[a]]
        PW += [np.full(2 * npair, wf)]
        PM += [np.full(2 * npair, ms[s])]
    if I:
        ds.pI, ds.pJ = np.concatenate(I), np.concatenate(J)
        ds.pw0, ds.pm = np.concatenate(PW), np.concatenate(PM)
        ds.n_pairs = int(ds.pI.size // 2)
        ds.n_families = int(np.unique(ds.pm).size)
    else:
        ds.pI = ds.pJ = np.zeros(0, int)
        ds.pw0 = np.zeros(0)
        ds.pm = np.zeros(0, int)
    return ds


def build_sets(D: N79) -> None:
    K = D.kids
    xs = np.isin(D.sample, (5, 6, 7, 8)) & D.msamp
    D.sets["main"] = build_dyadset(D, "main")
    D.sets["a"] = build_dyadset(D, "a", mother_mask=xs)
    D.sets["b"] = build_dyadset(D, "b", w_mode="unit")
    D.sets["c"] = build_dyadset(D, "c", w_mode="child")
    D.sets["d"] = build_dyadset(D, "d", kid_mask=(K.cyrb >= 1981) & (K.cyrb <= 2000))
    D.sets["e"] = build_dyadset(D, "e", w_mode="mother_only", first_only=True)
    D.sets["f"] = build_dyadset(D, "f", var="w10")
    D.sets["g"] = build_dyadset(D, "g", var="linear")
    D.sets["h"] = D.sets["main"]
    D.sets["i"] = D.sets["main"]
    D.sets["j"] = D.sets["main"]
    D.stab = build_stability(K)


# ---- NLSY97 ------------------------------------------------------------------

@dataclass
class N97:
    n: int
    female: np.ndarray
    w: np.ndarray
    theta: np.ndarray        # n x 4 raw CAT-ASVAB thetas
    postvar: np.ndarray
    zn: np.ndarray           # n x 4 age-normed (birth year x quarter)
    p: np.ndarray            # person-year frame (already filtered)
    Y: np.ndarray
    age: np.ndarray
    lw: np.ndarray
    age_d: np.ndarray = None
    yr_d: np.ndarray = None
    sample_type: np.ndarray = None
    n_rows_pos_neg: dict = field(default_factory=dict)


def build_n97(c: Cohort) -> N97:
    n = c.df.shape[0]
    female = c.get("KEY!SEX", 1997) == 2
    w = c.get("SAMPLING_WEIGHT_CC", 1997, scale=100.0)
    by, bm = c.get("KEY!BDATE_Y", 1997), c.get("KEY!BDATE_M", 1997)
    cell = (by * 10 + (bm - 1) // 3).astype(int)
    th = np.full((n, 4), np.nan)
    pv = np.full((n, 4), np.nan)
    info = {}
    for k, t in enumerate(TESTS79):
        pos = c.get(f"ASVAB_{t}_ABILITY_EST_POS", "XRND")
        neg = c.get(f"ASVAB_{t}_ABILITY_EST_NEG", "XRND")
        info[t] = dict(both=int((np.isfinite(pos) & np.isfinite(neg)).sum()),
                       neither=int((~np.isfinite(pos) & ~np.isfinite(neg)).sum()))
        th[:, k] = np.where(np.isfinite(pos), pos / 1000.0,
                            np.where(np.isfinite(neg), -neg / 1000.0, np.nan))
        pv[:, k] = c.get(f"ASVAB_{t}_POST_VARIANCE", "XRND", scale=1000.0)
    zn = np.column_stack([rank_normal(th[:, k], w, cell) for k in range(4)])
    p_, Y_, a_, l_ = [], [], [], []
    sexes = female
    for Y in sorted(int(y) for y in c.by_q["YINC-1700"]):
        age = c.get("CV_AGE_INT_DATE", Y)
        inc = c.get("YINC-1700", Y)
        t = Y - 1
        tt = f"{t % 100:02d}"
        wk = c.try_get(f"CVC_WKSWK_YR_ALL.{tt}", "XRND")
        hr = c.try_get(f"CVC_HOURS_WK_YR_ALL.{tt}", "XRND")
        if wk is None or hr is None:
            continue
        hg = c.try_get("CV_HGC_EVER_EDT", Y)
        hg = np.full(n, np.nan) if hg is None else np.where(hg >= 95, np.nan, hg)
        with np.errstate(divide="ignore", invalid="ignore"):
            ok = ((age >= AGE_LO) & (age <= AGE_HI) & (inc > 0) & (wk >= FY_WEEKS)
                  & (hr > 0) & (hr / wk >= FT_HOURS) & ~(hg < MIN_HGC))
        i = np.flatnonzero(ok)
        p_.append(i)
        Y_.append(np.full(i.size, Y))
        a_.append(age[i])
        l_.append(np.log(inc[i] / hr[i]))
    p, Yv, av, lw = (np.concatenate(v) for v in (p_, Y_, a_, l_))
    R = N97(n, female, w, th, pv, zn, p, Yv, av, lw)
    R.age_d = dummies(av.astype(int), range(AGE_LO, AGE_HI + 1))
    R.yr_d = dummies(Yv.astype(int), sorted(set(Yv.tolist())))
    R.sample_type = c.get("CV_SAMPLE_TYPE", 1997)
    R.n_rows_pos_neg = info
    return R


# =============================================================================
# Estimators (each returns flat dict of scalars; re-run on every bootstrap draw)
# =============================================================================

def standardize(x: np.ndarray, w: np.ndarray) -> np.ndarray:
    return (x - wmean(x, w)) / np.sqrt(wvar(x, w))


def composite(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Equal-weight mean of the available standardised columns."""
    Zs = np.column_stack([standardize(X[:, k], w) for k in range(X.shape[1])])
    with np.errstate(all="ignore"):
        return np.nanmean(Zs, axis=1)


def _put(out: dict, prefix: str, names: Sequence[str], vals: Sequence[float]) -> None:
    for nme, v in zip(names, vals):
        out[f"{prefix}{nme}"] = float(v)


def pop_afqt(Zfull: np.ndarray, m: np.ndarray, w: np.ndarray
             ) -> tuple[np.ndarray, float, float, FactorFit, np.ndarray]:
    """One-factor fit on the four z's in sample `m`, R_comp, omega, and afqt_c =
    equal-weight mean of the four z's each standardised (weighted) in `m`, then
    re-standardised in `m`.  Everything uses weights `w`."""
    Z = Zfull[m]
    wm = w[m]
    R = wcorr_cross(Z, Z, wm)
    fit = fit_factor_uls(R, [0] * 4)
    rc, om = comp_reliability(fit.lam, R)
    Zs = np.column_stack([standardize(Z[:, k], wm) for k in range(4)])
    afqt = np.full(len(w), np.nan)
    afqt[m] = standardize(Zs.mean(axis=1), wm)
    return afqt, rc, om, fit, R


def block_a(D: N79, w: np.ndarray, out: dict, info: dict | None) -> tuple[np.ndarray, float]:
    m = D.asamp
    wa = w * m
    afqt, rc, om, fit, _ = pop_afqt(D.z, m, w)
    out["conv.A"] = float(fit.converged)
    _put(out, "A.lam_", TESTS79, fit.lam)
    _put(out, "A.comm_", TESTS79, fit.comm)
    out.update({"A.rmsr": fit.rmsr, "R_comp_pop": rc, "A.omega": om,
                "A.heywood": float(fit.heywood)})
    # A2: KR-20 on number-right scores, persons with a valid IRT z for the subtest
    for k, t in enumerate(TESTS79):
        mk = np.isin(D.sample, POP_SAMPLES) & np.isfinite(D.z[:, k]) & np.isfinite(w)
        it = np.nan_to_num(D.items[t][mk], nan=0.0)
        ww = w[mk]
        ps = (it * ww[:, None]).sum(0) / ww.sum()
        vs = wvar(it.sum(1), ww)
        kk = it.shape[1]
        out[f"A.alpha_{t}"] = kk / (kk - 1) * (1 - (ps * (1 - ps)).sum() / vs)
    if info is not None:
        info["n.A"] = int(m.sum())
        info["n.A_alpha"] = {t: int((np.isin(D.sample, POP_SAMPLES) & np.isfinite(D.z[:, k])).sum())
                             for k, t in enumerate(TESTS79)}
        info["fit.A"] = fit
    return afqt, rc


def own_normed(D: N79, w: np.ndarray) -> np.ndarray:
    """A4: section standard scores rank-normed within birth-year x 4-month cells,
    NLSY79 1979 weight, reference = every respondent with a score."""
    return np.column_stack([rank_normal(D.sec[:, k], w, D.birth_cell) for k in range(4)])


def block_a4(D: N79, w: np.ndarray, Zo: np.ndarray, afqt: np.ndarray, out: dict,
             info: dict | None) -> tuple[np.ndarray, float]:
    # evaluated on the Block A sample (intersection with the own-normed sample)
    mo = D.asamp & np.isfinite(Zo).all(axis=1)
    afqt_o, rc, _, fit, _ = pop_afqt(Zo, mo, w)
    out["conv.A4"] = float(fit.converged)
    out["A4.R_comp_own"] = rc
    _put(out, "A4.lam_", TESTS79, fit.lam)
    out["A4.corr_afqt"] = wcorr(afqt, afqt_o, w)
    if info is not None:
        info["n.A4"] = int(mo.sum())
    return afqt_o, rc


# ---- block C -----------------------------------------------------------------

def _design(pd_: PY79 | N97, rows: np.ndarray, extra: list[np.ndarray]) -> np.ndarray:
    return np.column_stack([np.ones(rows.sum())] + [e[rows] for e in extra]
                           + [pd_.age_d[rows], pd_.yr_d[rows]])


def slope_set(py: PY79, y: np.ndarray, w: np.ndarray, afqt: np.ndarray,
              female: np.ndarray, base: np.ndarray, sqrtR: float, prefix: str,
              out: dict, info: dict | None, extras: bool = False) -> None:
    """Trim within year, then weighted OLS of y on afqt (+ female) + age + year
    dummies: pooled, male, female.  `base` = person-year mask before trimming;
    `w` = person-year weight (0 excludes)."""
    rows = base & (w > 0) & np.isfinite(afqt) & np.isfinite(y)
    keep = np.zeros(len(y), bool)
    r = np.flatnonzero(rows)
    keep[r] = trim_mask(y[r], w[r], py.Y[r])
    rows = keep
    for lab, sel, use_f in (("pooled", rows, True), ("male", rows & ~female, False),
                            ("female", rows & female, False)):
        ex = [afqt] + ([female.astype(float)] if use_f else [])
        X = _design(py, sel, ex)
        b = wls(y[sel], X, w[sel])[1]
        out[f"{prefix}wage_slope_obs_{lab}"] = b
        out[f"{prefix}wage_slope_latent_{lab}"] = b / sqrtR
        if info is not None:
            info[f"n.{prefix}py_{lab}"] = int(sel.sum())
            info[f"n.{prefix}pers_{lab}"] = int(np.unique(py.p[sel]).size)
    if extras:
        # 90/10 of year-demeaned log wage within each group
        for lab, sel in (("pooled", rows), ("male", rows & ~female), ("female", rows & female)):
            yy = y[sel].copy()
            ww = w[sel]
            for Yv in np.unique(py.Y[sel]):
                k = py.Y[sel] == Yv
                yy[k] -= wmean(yy[k], ww[k])
            if lab == "pooled":
                out["C.sd_afqt_wage_sample"] = float(np.sqrt(wvar(afqt[sel], ww)))
            out[f"C.p9010_{lab}"] = float(np.exp(wquantile(yy, ww, 0.9) - wquantile(yy, ww, 0.1)))


def block_c(D: N79, w: np.ndarray, afqt: np.ndarray, Rc: float, afqt_o: np.ndarray,
            Rc_o: float, out: dict, info: dict | None) -> None:
    py = D.py
    female = D.female[py.p]
    a_py = afqt[py.p]
    wpy = w[py.p]
    slope_set(py, py.lw, wpy, a_py, female, py.fy, np.sqrt(Rc), "", out, info, extras=True)
    out["implied_s_z"] = out["wage_slope_latent_pooled"] * (1 - ETA)
    out["gap_to_chi"] = out["wage_slope_latent_pooled"] - CHI
    # sensitivities (pooled slope only)
    tmp: dict = {}
    xs = np.isin(D.sample, XS_SAMPLES)
    wxs = w * xs
    afx, rcx, _, fx, _ = pop_afqt(D.z, D.asamp & xs, w)
    out["conv.C.xs"] = float(fx.converged)
    slope_set(py, py.lw, wxs[py.p], afx[py.p], female, py.fy, np.sqrt(rcx), "sens.xs.", tmp, None)
    out["sens.C.xs"] = tmp["sens.xs.wage_slope_latent_pooled"]
    out["sens.C.xs_R"] = rcx
    tmp = {}
    # fully unweighted: block A (standardisation, R_comp_pop) re-run with unit weights
    afu, rcu, _, fu, _ = pop_afqt(D.z, D.asamp, D.mult_p)
    out["conv.C.unw"] = float(fu.converged)
    slope_set(py, py.lw, D.mult_p[py.p], afu[py.p], female, py.fy, np.sqrt(rcu), "sens.unw.", tmp, None)
    out["sens.C.unw"] = tmp["sens.unw.wage_slope_latent_pooled"]
    tmp = {}
    slope_set(py, py.lw, wpy, afqt_o[py.p], female, py.fy, np.sqrt(Rc_o), "sens.own.", tmp, None)
    out["sens.C.own"] = tmp["sens.own.wage_slope_latent_pooled"]
    out["sens.C.own_obs"] = tmp["sens.own.wage_slope_obs_pooled"]
    tmp = {}
    lw_h = np.where(py.fy & (py.hrp1 > 0), np.log(np.where(py.hrp1 > 0, py.hrp1, 1.0)), np.nan)
    slope_set(py, lw_h, wpy, a_py, female, py.fy & np.isfinite(lw_h), np.sqrt(Rc), "sens.hrp.", tmp, None)
    out["sens.C.hrp"] = tmp["sens.hrp.wage_slope_latent_pooled"]
    # C5: participation, all Block A person-years at 25-34 (hgc >= 9, not armed forces)
    base = py.elig & (wpy > 0) & np.isfinite(a_py)
    for oname, yv in (("fyft", py.fyft), ("hp", py.hp)):
        for lab, sel in (("pooled", base), ("male", base & ~female), ("female", base & female)):
            s = sel & np.isfinite(yv)
            ex = [a_py] + ([female.astype(float)] if lab == "pooled" else [])
            X = _design(py, s, ex)
            out[f"C5.{oname}_{lab}"] = wls(yv[s], X, wpy[s])[1]
            if info is not None and oname == "fyft":
                info[f"n.C5_{lab}"] = int(s.sum())


# ---- block B -----------------------------------------------------------------

def est_B(ds: DyadSet, mult_p: np.ndarray, out: dict, prefix: str, full: bool,
          drop: Sequence[tuple[int, int]] = (), xm: np.ndarray | None = None,
          info: dict | None = None) -> FactorFit:
    w = ds.w0 * mult_p[ds.mrow]
    X = np.hstack([ds.xm if xm is None else xm, ds.xc])
    R = wcorr_cross(X, X, w)
    fit = fit_factor_uls(R, GROUPS_PC, drop)
    lam = fit.lam
    RM, oM = comp_reliability(lam[:4], R[:4, :4])
    RC4, oC = comp_reliability(lam[4:], R[4:, 4:])           # all four tests present
    RC = comp_reliability_avail(lam[4:], R[4:, 4:], X[:, 4:], w)
    out[f"{prefix}rho_pc_latent"] = fit.rho
    out[f"{prefix}R_comp_mother"] = RM
    out[f"{prefix}R_comp_child"] = RC
    out[f"{prefix}R_comp_child_all4"] = RC4
    out[f"conv.B.{prefix or 'main.'}"] = float(fit.converged)
    if full:
        _put(out, "B.lam_", IND8, lam)
        _put(out, "B.comm_", IND8, fit.comm)
        out.update({"B.rmsr": fit.rmsr, "B.omega_mother": oM, "B.omega_child": oC,
                    "B.heywood": float(fit.heywood)})
        cm, cc = composite(X[:, :4], w), composite(X[:, 4:], w)
        r = wcorr(cm, cc, w)
        out["r_pc_composite"] = r
        out["B3.check_rho"] = r / np.sqrt(RM * RC)
        # residual structure: mean within-domain vs across-domain cross-correlations
        if info is not None:
            info["fit.B"] = fit
            info["R.B"] = R
    return fit


def est_S(ds: DyadSet, mult_p: np.ndarray, lam_c: np.ndarray, rho_pc: float, out: dict,
          info: dict | None) -> None:
    pw = ds.pw0 * mult_p[ds.pm]
    w = ds.w0 * mult_p[ds.mrow]
    Xc = ds.xc
    C = wcorr_cross(Xc[ds.pI], Xc[ds.pJ], pw)
    ll = np.outer(lam_c, lam_c)
    off = ~np.eye(4, dtype=bool)
    rs_off = np.nansum(C[off] * ll[off]) / np.sum(ll[off] ** 2)
    rs_all = np.nansum(C * ll) / np.sum(ll ** 2)
    out["rho_sib_latent"] = rs_off
    out["rho_sib_latent_incl"] = rs_all
    for k in range(4):
        for l in range(4):
            out[f"S.C_{IND8[4 + k]}_{IND8[4 + l]}"] = C[k, l]
    comp = composite(Xc, w)
    out["S.sib_comp_corr"] = wcorr(comp[ds.pI], comp[ds.pJ], pw)
    out["implied_rho_from_sib"] = rs_off / rho_pc


def est_stability(D: N79, mult_p: np.ndarray, out: dict) -> None:
    K = D.kids
    km = np.where(K.mrow >= 0, mult_p[np.maximum(K.mrow, 0)], 0.0)
    for t, (ch, z1, z2, w) in D.stab.items():
        out[f"B4.stab_{t}"] = wcorr(z1, z2, w * km[ch])


def block_b(D: N79, mult_p: np.ndarray, Zo: np.ndarray, out: dict, info: dict | None) -> None:
    S = D.sets
    fit = est_B(S["main"], mult_p, out, "", True, info=info)
    est_S(S["main"], mult_p, fit.lam[4:], fit.rho, out, info)
    est_stability(D, mult_p, out)
    sens = {"a": {}, "b": {}, "c": {}, "d": {}, "e": {}, "f": {}, "g": {},
            "h": dict(drop=DROP_DOMAIN), "i": dict(xm=Zo[S["main"].mrow]),
            "j": dict(drop=((5, 6),))}
    for k, kw in sens.items():
        est_B(S[k], mult_p, out, f"sens.B.{k}.", False, **kw)
    if info is not None:
        for k, ds in S.items():
            info[f"n.B_{k}"] = (ds.n_mothers, ds.n_children)
        s = S["main"]
        info["n.B"] = (s.n_mothers, s.n_children)
        info["n.S"] = (s.n_families, s.n_pairs)
        K = D.kids
        mn = K.counts["main"][s.krow]
        info["B.mean_assess"] = float(K.nrounds[s.krow].mean())
        info["B.share_ppvt"] = float((mn[:, 3] > 0).mean())
        info["B.mean_scores"] = mn.mean(0).tolist()
        info["B.share_female_child"] = float(K.female[s.krow].mean())
        info["n.kids_all"] = int(K.n)
        info["n.assess_all"] = int(K.norm_info["n_assess"])
        info["n.comp_all"] = int(K.norm_info["n_comp"])
        info["n.imputed_comp"] = int(K.norm_info["n_imputed_comp"])
        info["n.age_fallback_msage"] = int(K.norm_info["n_age_fallback_msage"])


# ---- block E -----------------------------------------------------------------

def _resid_perm(y: np.ndarray, pid: np.ndarray, Y: np.ndarray, age_d: np.ndarray,
                yr_d: np.ndarray, w: np.ndarray, npers: int) -> tuple[np.ndarray, np.ndarray]:
    """Trim log income 1/99 within year, residualise on age and year dummies
    (weighted), and return each person's mean residual and number of years."""
    rows = np.flatnonzero((w > 0) & np.isfinite(y))
    keep = rows[trim_mask(y[rows], w[rows], Y[rows])]
    X = np.column_stack([np.ones(keep.size), age_d[keep], yr_d[keep]])
    b = wls(y[keep], X, w[keep])
    res = y[keep] - X @ b
    cnt = np.bincount(pid[keep], minlength=npers)
    sm = np.bincount(pid[keep], weights=res, minlength=npers)
    perm = np.full(npers, np.nan)
    perm[cnt > 0] = sm[cnt > 0] / cnt[cnt > 0]
    return perm, cnt


def block_e(D: N79, w: np.ndarray, mult_p: np.ndarray, out: dict, info: dict | None) -> None:
    py, K, ya = D.py, D.kids, D.ya
    # mothers
    mrows = D.msamp[py.p] & (py.inc > 0)
    ym = np.where(mrows, np.log(np.where(py.inc > 0, py.inc, 1.0)), np.nan)
    wm = w[py.p] * mrows
    perm_m, cnt_m = _resid_perm(ym, py.p, py.Y, py.age_d, py.yr_d, wm, D.n)
    # young adults, by sex
    ky = ya.krow
    mo = K.mrow[ky]
    okya = D.msamp[mo]
    yy = np.where(okya, np.log(ya.inc), np.nan)
    wy = ya.w * mult_p[mo] * okya
    perm_k = np.full(K.n, np.nan)
    cnt_k = np.zeros(K.n, int)
    for sx in (False, True):
        s = K.female[ky] == sx
        pk, ck = _resid_perm(np.where(s, yy, np.nan), ky, ya.Y, ya.age_d, ya.yr_d, wy * s, K.n)
        perm_k = np.where(ck > 0, pk, perm_k)
        cnt_k += ck
    for tag, minyrs in (("", 1), ("_2yr", 2)):
        kid = (cnt_k >= minyrs) & (K.mrow >= 0)
        kid &= np.where(K.mrow >= 0, cnt_m[np.maximum(K.mrow, 0)] >= minyrs, False)
        ki = np.flatnonzero(kid)
        mm = K.mrow[ki]
        nk = np.bincount(mm, minlength=D.n)[mm]
        wd = w[mm] / nk
        xm_, yk_ = perm_m[mm], perm_k[ki]
        out[f"E1.corr{tag}"] = wcorr(xm_, yk_, wd)
        cv = wcov(xm_, yk_, wd)
        out[f"E1.cov{tag}"] = cv
        out[f"E1.slope{tag}"] = cv / wvar(xm_, wd)
        if tag == "":
            for lab, sel in (("male", ~K.female[ki]), ("female", K.female[ki])):
                out[f"E1.corr_{lab}"] = wcorr(xm_[sel], yk_[sel], wd[sel])
                out[f"E1.slope_{lab}"] = wcov(xm_[sel], yk_[sel], wd[sel]) / wvar(xm_[sel], wd[sel])
            if info is not None:
                info["n.E1"] = (int(np.unique(mm).size), int(ki.size))
                info["n.E1_male"] = int((~K.female[ki]).sum())
                info["n.E1_female"] = int(K.female[ki].sum())
                info["n.E_mothers"] = int((cnt_m > 0).sum())
                info["n.E_py_mothers"] = int(mrows.sum())
                info["n.E_py_ya"] = int(okya.sum())
    # E2 siblings: all YA children of sampled mothers with permanent earnings
    have = (cnt_k > 0) & (K.mrow >= 0)
    hi = np.flatnonzero(have)
    hm = K.mrow[hi]
    o = np.argsort(hm, kind="stable")
    hi, hm = hi[o], hm[o]
    starts = np.flatnonzero(np.concatenate([[True], hm[1:] != hm[:-1]])) if hi.size else []
    ends = np.append(starts[1:], hi.size) if hi.size else []
    A, Bk, PWt = [], [], []
    for s_, e_ in zip(starts, ends):
        k = e_ - s_
        if k < 2:
            continue
        a, b = np.triu_indices(k, 1)
        g = hi[s_:e_]
        A += [g[a], g[b]]
        Bk += [g[b], g[a]]
        PWt += [np.full(2 * a.size, w[hm[s_]] / (2 * a.size))]
    if A:
        A, Bk, PWt = np.concatenate(A), np.concatenate(Bk), np.concatenate(PWt)
        out["E2.corr_sib"] = wcorr(perm_k[A], perm_k[Bk], PWt)
        if info is not None:
            info["n.E2"] = (int(np.unique(K.mrow[A]).size), int(A.size // 2))
    else:
        out["E2.corr_sib"] = np.nan
    out["E3.implied_rho"] = out["E2.corr_sib"] / out["E1.corr"]


# ---- block V -----------------------------------------------------------------

def run97(R: N97, mult: np.ndarray, info: dict | None) -> dict[str, float]:
    out: dict[str, float] = {}
    w = R.w * mult
    ok = np.isfinite(R.zn).all(axis=1)
    m = ok & (w > 0)
    Z = R.zn[m]
    Rm = wcorr_cross(Z, Z, w[m])
    fit = fit_factor_uls(Rm, [0] * 4)
    rc, om = comp_reliability(fit.lam, Rm)
    _put(out, "V.lam_", TESTS79, fit.lam)
    out["R_comp_97"] = rc
    out["V.omega"] = om
    out["V.heywood"] = float(fit.heywood)
    out["conv.V"] = float(fit.converged)
    for k, t in enumerate(TESTS79):
        okk = np.isfinite(R.theta[:, k]) & (w > 0)
        vt, pvm = wvar(R.theta[okk, k], w[okk]), wmean(R.postvar[okk, k], w[okk])
        out[f"V.postvar_rel_{t}"] = vt / (vt + pvm)      # Bayesian thetas: var(theta_hat)/(var + E[post_var])
    comp = np.full(R.n, np.nan)
    comp[m] = Z.mean(1)
    af = np.full(R.n, np.nan)
    af[m] = standardize(comp[m], w[m])
    female = R.female[R.p]
    wpy = w[R.p]
    tmp: dict = {}
    slope_set(R, R.lw, wpy, af[R.p], female, np.ones(len(R.p), bool), np.sqrt(rc), "V.", tmp, info)
    for lab in ("pooled", "male", "female"):
        out[f"wage_slope_obs97_{lab}"] = tmp[f"V.wage_slope_obs_{lab}"]
        out[f"wage_slope_latent97_{lab}"] = tmp[f"V.wage_slope_latent_{lab}"]
    if info is not None:
        info["n.V_persons"] = int(m.sum())
        info["n.V_persons_all"] = int(R.n)
        info["fit.V"] = fit
    return out


# =============================================================================
# Orchestration
# =============================================================================

def run79(D: N79, mult_hh: np.ndarray | None, info: dict | None = None) -> dict[str, float]:
    mult_p = np.ones(D.n) if mult_hh is None else mult_hh[D.hh_idx].astype(float)
    D.mult_p = mult_p
    w = D.w * mult_p
    out: dict[str, float] = {}
    afqt, Rc = block_a(D, w, out, info)
    Zo = own_normed(D, w)
    afqt_o, Rc_o = block_a4(D, w, Zo, afqt, out, info)
    block_c(D, w, afqt, Rc, afqt_o, Rc_o, out, info)
    block_b(D, mult_p, Zo, out, info)
    block_e(D, w, mult_p, out, info)
    return out


def run_all(D: N79, R97: N97, mh: np.ndarray | None, m97: np.ndarray | None,
            info: dict | None = None) -> dict[str, float]:
    out = run79(D, mh, info)
    out.update(run97(R97, np.ones(R97.n) if m97 is None else m97, info))
    return out


def bootstrap(D: N79, R97: N97, B: int, seed: int) -> tuple[dict[str, np.ndarray], list[float]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, float]] = []
    times: list[float] = []
    for b in range(B):
        t0 = time.time()
        mh = rng.multinomial(D.nhh, np.full(D.nhh, 1.0 / D.nhh))
        m97 = rng.multinomial(R97.n, np.full(R97.n, 1.0 / R97.n))
        try:
            rows.append(run_all(D, R97, mh, m97))
        except Exception as exc:                       # a failed draw is a NaN row
            LOG.warning("  draw %d failed: %s", b, exc)
            rows.append({})
        times.append(time.time() - t0)
        if (b + 1) % 10 == 0 or b + 1 == B:
            LOG.info("  bootstrap %d/%d  (%.2fs/draw)", b + 1, B, float(np.mean(times)))
    keys = sorted({k for r in rows for k in r})
    return {k: np.array([r.get(k, np.nan) for r in rows]) for k in keys}, times


# =============================================================================
# Records
# =============================================================================

class Boot:
    def __init__(self, point: dict[str, float], draws: dict[str, np.ndarray]) -> None:
        self.point, self.draws = point, draws

    def val(self, k: str) -> float | None:
        v = self.point.get(k)
        return None if v is None or not np.isfinite(v) else float(v)

    def se(self, k: str) -> float | None:
        d = self.draws.get(k)
        if d is None or np.isfinite(d).sum() < 2:
            return None
        return float(np.nanstd(d, ddof=1))

    def ci(self, k: str) -> tuple[float | None, float | None]:
        d = self.draws.get(k)
        if d is None or np.isfinite(d).sum() < 2:
            return None, None
        return (float(np.nanpercentile(d, 2.5)), float(np.nanpercentile(d, 97.5)))

    def pair(self, k: str) -> dict[str, float | None]:
        lo, hi = self.ci(k)
        return {"value": self.val(k), "se": self.se(k), "ci_lo": lo, "ci_hi": hi}


SRC79 = "NLSY79 (ASVAB 1980-81; interviews 1979-2000), weights = 1979 SAMPWEIGHT"
SRCB = "NLSY79 women x CNLSY children (PIAT/PPVT 1986-2014)"
SRCC = "NLSY79, ages 25-34, full-year full-time wage workers"
SRC97 = "NLSY97 CAT-ASVAB; rounds with ages 25-34"
OBJ_RHO = "ρz (plus the Q_l channel and assortative mating)"
OBJ_B = "b·s_z, b = 1/(1−η)"
OBJ_U = "var(u) in the score measurement equation"
OBJ_CHK = "check"


def make_records(reg: Registry, bt: Boot, info: dict, args: argparse.Namespace) -> None:
    def add(key, label, obj, stat, units, source, note, n=None, extra=None, sk=None):
        lo, hi = bt.ci(stat)
        reg.add(Estimate(key, label, obj, bt.val(stat), bt.se(stat), lo, hi, n, units,
                         source, note, extra or {}))

    def ex(prefix: str, names: Sequence[str]) -> dict:
        return {n: bt.pair(f"{prefix}{n}") for n in names}

    nB, nS, nA = info["n.B"], info["n.S"], info["n.A"]
    add("rho_pc_latent", "Latent mother-child ability correlation (two-factor ULS)", OBJ_RHO,
        "rho_pc_latent", "correlation", SRCB,
        "Primary moment for rho_z. Mother factor: AR, WK, PC, MK (IRT z). Child factor: PIAT math, "
        "recognition, comprehension, PPVT, age-normed.", nB[1],
        {"lam": ex("B.lam_", IND8), "comm": ex("B.comm_", IND8), "rmsr": bt.pair("B.rmsr"),
         "n_mothers": nB[0], "n_children": nB[1],
         "heywood": bool(info["fit.B"].heywood), "converged": bool(info["fit.B"].converged)})
    add("R_comp_mother", "Reliability of the equal-weight mother composite (4 IRT z)", OBJ_U,
        "R_comp_mother", "share", SRCB, "(Σλ)²/(1'S1), S = sample correlations, dyad sample.",
        nB[0], {"omega": bt.pair("B.omega_mother")})
    add("R_comp_child", "Reliability of the child composite actually used (mean of available tests)", OBJ_U,
        "R_comp_child", "share", SRCB, "(Σλ)²/(1'S1), S = sample correlations of the child-level scores.",
        nB[1], {"omega": bt.pair("B.omega_child"), "R_comp_child_all4": bt.pair("R_comp_child_all4")})
    add("r_pc_composite", "Weighted correlation of the two equal-weight composites", OBJ_CHK,
        "r_pc_composite", "correlation", SRCB, "Observed mother-child score correlation.", nB[1],
        {"check_r_over_sqrt_RMRC": bt.pair("B3.check_rho")})
    for t in CH_TESTS:
        add(f"B4.stab_{t}", f"Two-year stability of age-normed {t}", OBJ_CHK, f"B4.stab_{t}",
            "correlation", "CNLSY, rounds y and y+2, ages 5-14",
            "Lower bound on single-assessment reliability (not used).", None)
    add("R_comp_pop", "Reliability of the 4-subtest AFQT composite, NLSY79 population", OBJ_U,
        "R_comp_pop", "share", SRC79, "One-factor ULS on the four IRT z-scores; equal-weight composite.",
        nA, {"lam": ex("A.lam_", TESTS79), "comm": ex("A.comm_", TESTS79),
             "omega": bt.pair("A.omega"), "rmsr": bt.pair("A.rmsr"),
             "kr20": ex("A.alpha_", TESTS79),
             "heywood": bool(info["fit.A"].heywood)})
    for lab in ("pooled", "male", "female"):
        add(f"wage_slope_obs_{lab}", f"Log-wage slope on observed afqt_c ({lab})", OBJ_CHK,
            f"wage_slope_obs_{lab}", "log points per SD of afqt_c", SRCC,
            "Weighted OLS with age and year dummies (+ female when pooled); trimmed 1/99 within year.",
            info[f"n.py_{lab}"], {"n_persons": info[f"n.pers_{lab}"]})
        add(f"wage_slope_latent_{lab}", f"Log-wage slope per SD of latent ability ({lab})", OBJ_B,
            f"wage_slope_latent_{lab}", "log points per SD of f", SRCC,
            "slope_obs / sqrt(R_comp_pop)." + (" Primary moment for b·s_z." if lab == "pooled" else ""),
            info[f"n.py_{lab}"])
    add("implied_s_z", "Implied stationary SD of log ability, s_z = slope_latent (1−η)", "s_z",
        "implied_s_z", "log points", SRCC, f"η = {ETA}.", info["n.py_pooled"])
    add("gap_to_chi", "Pooled latent slope minus χ = log 1.12", OBJ_B, "gap_to_chi", "log points",
        SRCC, f"χ = {CHI:.5f}. Positive: scores predict wages more than the CFR bridge.",
        info["n.py_pooled"])
    for lab in ("pooled", "male", "female"):
        add(f"C.p9010_{lab}", f"Weighted 90/10 of year-demeaned hourly wage ({lab})", OBJ_CHK,
            f"C.p9010_{lab}", "ratio", SRCC, "Trimmed sample; ACS 2009-13 pooled is 3.74 (different vintage).",
            info[f"n.py_{lab}"])
    for oname, desc in (("fyft", "1{full-year full-time wage worker}"), ("hp", "1{annual hours < 780}")):
        for lab in ("pooled", "male", "female"):
            add(f"C5.{oname}_{lab}", f"LPM slope of {desc} on afqt_c ({lab})", OBJ_CHK,
                f"C5.{oname}_{lab}", "probability per SD", SRC79 + ", ages 25-34, hgc >= 9",
                "Model makes work independent of z; a nonzero slope is a caveat.", info[f"n.C5_{lab}"])
    add("rho_sib_latent", "Latent sibling correlation (off-diagonal cross-test)", OBJ_CHK,
        "rho_sib_latent", "correlation", SRCB,
        "Σ_{k≠l} C_kl λ_k λ_l / Σ_{k≠l} (λ_k λ_l)² with child loadings from B2.", nS[1],
        {"including_k_eq_l": bt.pair("rho_sib_latent_incl"), "n_families": nS[0], "n_pairs": nS[1]})
    add("S.sib_comp_corr", "Raw sibling correlation of the child composite", OBJ_CHK, "S.sib_comp_corr",
        "correlation", SRCB, "Ordered sibling pairs.", nS[1])
    add("implied_rho_from_sib", "Implied rho = rho_sib_latent / rho_pc_latent", OBJ_CHK,
        "implied_rho_from_sib", "ratio", SRCB,
        "AR(1): siblings correlate ρ², mother-child ρ. Expected ≥ 1 (shared family inputs).", nS[1])
    for lab in ("pooled", "male", "female"):
        add(f"wage_slope_latent97_{lab}", f"NLSY97 latent log-wage slope ({lab})", OBJ_B,
            f"wage_slope_latent97_{lab}", "log points per SD of f", SRC97,
            "Same definition as block C; slope_obs / sqrt(R_comp_97).", info[f"n.V.py_{lab}"],
            {"obs": bt.pair(f"wage_slope_obs97_{lab}")})
    add("R_comp_97", "Reliability of the NLSY97 4-subtest composite", OBJ_U, "R_comp_97", "share",
        SRC97, "One-factor ULS on age-normed CAT-ASVAB thetas.", info["n.V_persons"],
        {"lam": ex("V.lam_", TESTS79), "postvar_reliability": ex("V.postvar_rel_", TESTS79),
         "omega": bt.pair("V.omega")})
    add("E1.corr", "Mother-child correlation of permanent log earnings", OBJ_CHK, "E1.corr",
        "correlation", "NLSY79 mothers x CNLSY young adults, ages 25-34", "Check, not a target.",
        info["n.E1"][1],
        {"cov": bt.pair("E1.cov"), "slope": bt.pair("E1.slope"),
         "male_children": bt.pair("E1.corr_male"), "female_children": bt.pair("E1.corr_female"),
         "slope_male": bt.pair("E1.slope_male"), "slope_female": bt.pair("E1.slope_female"),
         "at_least_2_years": {"corr": bt.pair("E1.corr_2yr"), "slope": bt.pair("E1.slope_2yr"),
                              "cov": bt.pair("E1.cov_2yr")}})
    add("E2.corr_sib", "Sibling correlation of permanent log earnings", OBJ_CHK, "E2.corr_sib",
        "correlation", "CNLSY young adults, ages 25-34", "Check, not a target.", info["n.E2"][1])
    add("E3.implied_rho", "Earnings-implied rho = sibling / mother-child correlation", OBJ_CHK,
        "E3.implied_rho", "ratio", "NLSY79 x CNLSY", "Check.", info["n.E1"][1])
    add("A4.corr_afqt", "Correlation of afqt_c (NLS age-normed) with own-normed afqt_c", OBJ_CHK,
        "A4.corr_afqt", "correlation", SRC79, "Sensitivity of the observed score to age-norming.",
        info["n.A4"], {"R_comp_own": bt.pair("A4.R_comp_own")})

    # ---- sensitivities -------------------------------------------------------
    sens_rows = [("main", "Main dyad sample", ""), ("a", "Cross-sectional mothers only (SAMPLE_ID 5-8)", "a"),
                 ("b", "Unweighted", "b"), ("c", "Child weights (mean CSAMWT_REV), no 1/n_c", "c"),
                 ("d", "Children born 1981-2000", "d"), ("e", "One child per mother (first-born)", "e"),
                 ("f", "Child scores from ages 10-14 only", "f"), ("g", "Linear age-norming", "g"),
                 ("h", "Drop within-domain pairs AR-MK, WK-PC, RECOG-COMP", "h"),
                 ("i", "Mothers' own-normed section scores", "i"),
                 ("j", "Drop only the RECOG-COMP pair", "j")]
    for k, lab, _ in sens_rows:
        pre = "" if k == "main" else f"sens.B.{k}."
        nn = info[f"n.B_{k}"] if k != "main" else info["n.B"]
        for st, sname in (("rho_pc_latent", "rho_pc_latent"), ("R_comp_mother", "R_comp_mother"),
                          ("R_comp_child", "R_comp_child"), ("R_comp_child_all4", "R_comp_child_all4")):
            lo, hi = bt.ci(pre + st)
            reg_s.add(Estimate(f"sens.B.{k}.{sname}", f"{lab}: {sname}", OBJ_RHO if st == "rho_pc_latent" else OBJ_U,
                               bt.val(pre + st), bt.se(pre + st), lo, hi, nn[1], "", SRCB,
                               lab, {"n_mothers": nn[0], "n_children": nn[1]}))
    for st, lab in (("xs", "Cross-sectional sample only (SAMPLE_ID 1-8)"), ("unw", "Fully unweighted (incl. standardization and R_comp_pop)"),
                    ("own", "Own-normed afqt (A4), Block A sample"), ("hrp", "Hourly wage from HRP1")):
        k = f"sens.C.{st}"
        lo, hi = bt.ci(k)
        reg_s.add(Estimate(f"sens.C.{st}.wage_slope_latent_pooled", f"{lab}: pooled latent wage slope",
                           OBJ_B, bt.val(k), bt.se(k), lo, hi, None, "log points per SD of f", SRCC, lab, {}))


reg_s = Registry()


# =============================================================================
# Report
# =============================================================================

def fmt(v: float | None, se: float | None = None, d: int = 3) -> str:
    if v is None:
        return "n/a"
    s = f"{v:.{d}f}"
    return s + (f" ({se:.{d}f})" if se is not None else "")


def write_reports(reg: Registry, meta: dict[str, Any], outdir: Path, bt: Boot) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {**meta, "estimates": reg.to_dict(), "sensitivities": reg_s.to_dict()}
    (outdir / "nlsy_ability.json").write_text(json.dumps(payload, indent=2, default=str))
    E = {e.key: e for e in reg.items}
    S = {e.key: e for e in reg_s.items}
    info = meta["samples"]
    L: list[str] = []
    add = L.append

    def row(k: str, d: int = 3) -> str:
        e = E[k]
        ci = "" if e.ci_lo is None else f"[{e.ci_lo:.{d}f}, {e.ci_hi:.{d}f}]"
        return f"| `{k}` | {e.label} | {fmt(e.value, e.se, d)} | {ci} | {e.n if e.n else '—'} |"

    hdr = "| key | quantity | estimate (SE) | 95% CI | n |\n|---|---|---|---|---|"
    add("# NLSY ability moments for the first-stage indirect inference (T2a)\n")
    add(f"Generated {meta['generated']} by `nlsy_ability.py`; cluster bootstrap "
        f"B = {meta['assumptions']['B']} (seed {meta['assumptions']['seed']}), SE = bootstrap SD, "
        f"CI = 95% percentile. Data side of T2b: every statistic is to be recomputed on simulated "
        f"model panels as defined in the last section.\n")
    add("## What feeds T2b\n")
    add(hdr)
    for k in ("rho_pc_latent", "wage_slope_latent_pooled", "R_comp_mother", "R_comp_child", "R_comp_pop",
              "implied_s_z", "gap_to_chi"):
        add(row(k))
    add(f"\nη = {ETA}, χ = log 1.12 = {CHI:.5f}. `rho_pc_latent` is the target for ρz; "
        "`wage_slope_latent_pooled` is the target for b·s_z; the R_comp's set var(u) in the score "
        "measurement equation (mother composite, child composite, population composite).\n")

    add("## A. NLSY79 population measurement model\n")
    e = E["R_comp_pop"]
    add(f"Sample: SAMPLE_ID 1-8, 10, 11, 13, 14 (drops poor-white supplement 9/12 and military 15-20), "
        f"all four IRT z valid, N = {info['n.A']}. Weight = 1979 SAMPWEIGHT. IRT z-scores are already "
        f"age-normed by NLS (four-month groups within birth year).\n")
    add("| subtest | loading | SE | communality | KR-20 α (number-right) |\n|---|---|---|---|---|")
    for t in TESTS79:
        l, c, a = e.extra["lam"][t], e.extra["comm"][t], e.extra["kr20"][t]
        add(f"| {t} | {l['value']:.3f} | {l['se']:.3f} | {c['value']:.3f} | {a['value']:.3f} |")
    add(f"\nR_comp_pop = {fmt(e.value, e.se)}; omega = {fmt(e.extra['omega']['value'], e.extra['omega']['se'])}; "
        f"RMSR = {e.extra['rmsr']['value']:.4f}; Heywood: {e.extra['heywood']}. "
        "KR-20 is for the number-right score, not the IRT score: a benchmark for the communalities. "
        f"Own-norming sensitivity (section standard scores ranked within birth-year x 4-month cells): "
        f"corr(afqt_c, afqt_c_own) = {fmt(E['A4.corr_afqt'].value, E['A4.corr_afqt'].se)}, "
        f"R_comp_own = {fmt(E['A4.corr_afqt'].extra['R_comp_own']['value'])}.\n")

    add("## B. Mother-child measurement system\n")
    e = E["rho_pc_latent"]
    add(f"Mothers: SAMPLE_ID 5-8, 13, 14 with four valid IRT z; children: CNLSY, PIAT math / recognition / "
        f"comprehension and PPVT-R at ages 5-14 (60-179 months), rank-normed within 3-month age cells over all "
        f"CNLSY assessments (reference weight = round child weight), averaged within child. Dyad weight = "
        f"mother's 1979 weight / number of her children in the sample. N mothers = {e.extra['n_mothers']}, "
        f"N children = {e.extra['n_children']}; mean assessment rounds per child = {info['B.mean_assess']:.2f}; "
        f"share with PPVT = {info['B.share_ppvt']:.2f}. Child age is `CSAGE` (child supplement), with `MSAGE` only "
        f"where CSAGE is missing ({info['n.age_fallback_msage']} assessments). **Imputed comprehension:** NLS gives "
        f"children with RECOG < 19 no comprehension test and copies the recognition score into COMP "
        f"({info['n.imputed_comp']} in-window assessments, {100 * info['n.imputed_comp'] / info['n.comp_all']:.1f}% of in-window COMP assessments). "
        "Those rows stay in the COMP norming reference (so real scores are ranked against the whole age cell) but are "
        "excluded from child-level COMP means and from the COMP two-year stabilities. `R_comp_child` is the reliability of "
        "the mean of the available tests.\n")
    add(hdr)
    for k in ("rho_pc_latent", "R_comp_mother", "R_comp_child", "r_pc_composite"):
        add(row(k))
    add("\n| indicator | loading | SE | communality |\n|---|---|---|---|")
    for t in IND8:
        add(f"| {t} | {e.extra['lam'][t]['value']:.3f} | {e.extra['lam'][t]['se']:.3f} | "
            f"{e.extra['comm'][t]['value']:.3f} |")
    add(f"\nRMSR = {e.extra['rmsr']['value']:.4f}; Heywood: {e.extra['heywood']}; omega (mother, child) = "
        f"{fmt(E['R_comp_mother'].extra['omega']['value'])}, {fmt(E['R_comp_child'].extra['omega']['value'])}; "
        f"check r/sqrt(R_M R_C) = {fmt(E['r_pc_composite'].extra['check_r_over_sqrt_RMRC']['value'])} "
        f"(SE {fmt(E['r_pc_composite'].extra['check_r_over_sqrt_RMRC']['se'])}) vs the ULS rho above.\n")
    add("Two-year child stabilities (not used):\n")
    add("| test | corr |\n|---|---|")
    for t in CH_TESTS:
        add(f"| {t} | {fmt(E[f'B4.stab_{t}'].value, E[f'B4.stab_{t}'].se)} |")
    add("\n### Sensitivities\n")
    add("| variant | rho_pc_latent | R_comp_mother | R_comp_child | N mothers | N children |\n|---|---|---|---|---|---|")
    labs = {"main": "Main"}
    for k in ("main", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"):
        if k == "main":
            r = [E["rho_pc_latent"], E["R_comp_mother"], E["R_comp_child"]]
            nn = (info["n.B"][0], info["n.B"][1])
            nm = "Main"
        else:
            r = [S[f"sens.B.{k}.{s}"] for s in ("rho_pc_latent", "R_comp_mother", "R_comp_child")]
            nn = info[f"n.B_{k}"]
            nm = r[0].note
        add(f"| {nm} | " + " | ".join(fmt(x.value, x.se) for x in r) + f" | {nn[0]} | {nn[1]} |")

    add("\n## S. Siblings\n")
    e = E["rho_sib_latent"]
    add(f"Families with ≥ 2 children in the dyad sample: {e.extra['n_families']}; sibling pairs: {e.extra['n_pairs']}. "
        "Ordered pairs, each family's weight = mother's weight.\n")
    add(hdr)
    for k in ("rho_sib_latent", "S.sib_comp_corr", "implied_rho_from_sib"):
        add(row(k))
    add(f"\nIncluding k = l terms: rho_sib_latent = {fmt(e.extra['including_k_eq_l']['value'], e.extra['including_k_eq_l']['se'])}.\n")
    add("Cross-sibling cross-test correlations C[k, l] (child k, sibling l):\n")
    add("| | " + " | ".join(CH_TESTS) + " |\n|---|" + "---|" * 4)
    for k in CH_TESTS:
        add(f"| {k} | " + " | ".join(f"{bt.val(f'S.C_{k}_{l}'):.3f}" for l in CH_TESTS) + " |")

    add("\n## C. Log-wage slope on latent AFQT (NLSY79)\n")
    add(f"Person-years of the block A sample, survey rounds with ages {AGE_LO}-{AGE_HI}; income > 0, weeks ≥ {FY_WEEKS}, "
        f"hours/week ≥ {FT_HOURS}, not armed forces (ESR ≠ 4), highest grade ≥ {MIN_HGC}; log(income / annual hours); "
        "weighted 1/99 trim within year; weight = 1979 SAMPWEIGHT; controls: age and year dummies (+ female).\n")
    add("| slope | pooled | male | female |\n|---|---|---|---|")
    for nm, pre in (("observed, per SD of afqt_c", "wage_slope_obs_"), ("latent, per SD of f", "wage_slope_latent_")):
        add(f"| {nm} | " + " | ".join(fmt(E[f'{pre}{l}'].value, E[f'{pre}{l}'].se) for l in ("pooled", "male", "female")) + " |")
    add("| N person-years | " + " | ".join(str(E[f'wage_slope_obs_{l}'].n) for l in ("pooled", "male", "female")) + " |")
    add("| N persons | " + " | ".join(str(E[f'wage_slope_obs_{l}'].extra['n_persons']) for l in ("pooled", "male", "female")) + " |")
    add("| 90/10 (year-demeaned) | " + " | ".join(fmt(E[f'C.p9010_{l}'].value, E[f'C.p9010_{l}'].se, 2) for l in ("pooled", "male", "female")) + " |")
    add(f"\nImplied s_z = {fmt(E['implied_s_z'].value, E['implied_s_z'].se)}; gap to χ = "
        f"{fmt(E['gap_to_chi'].value, E['gap_to_chi'].se)}.\n")
    add("Participation (LPM slope per SD of afqt_c; all block A person-years, hgc ≥ 9, not armed forces, no labor-supply filter):\n")
    add("| outcome | pooled | male | female |\n|---|---|---|---|")
    for o, nm in (("fyft", "1{FYFT wage worker}"), ("hp", "1{hours < 780}")):
        add(f"| {nm} | " + " | ".join(fmt(E[f'C5.{o}_{l}'].value, E[f'C5.{o}_{l}'].se) for l in ("pooled", "male", "female")) + " |")
    add("\nSensitivities (pooled latent slope):\n")
    add("| variant | slope (SE) |\n|---|---|")
    add(f"| Main | {fmt(E['wage_slope_latent_pooled'].value, E['wage_slope_latent_pooled'].se)} |")
    for st in ("xs", "unw", "own", "hrp"):
        s = S[f"sens.C.{st}.wage_slope_latent_pooled"]
        add(f"| {s.note} | {fmt(s.value, s.se)} |")

    add("\n## V. NLSY97 vintage check\n")
    e = E["R_comp_97"]
    add(f"CAT-ASVAB thetas age-normed within birth-year x quarter (1997 weight, fixed across draws), one-factor ULS: "
        f"R_comp_97 = {fmt(e.value, e.se)} (N = {e.n}). Wage sample as block C with `YINC-1700` and CVC hours/weeks of "
        "calendar year Y−1; both NLSY97 samples, weighted by 1997 SAMPLING_WEIGHT_CC; bootstrap by person.\n")
    add("| | pooled | male | female |\n|---|---|---|---|")
    add("| latent slope | " + " | ".join(fmt(E[f'wage_slope_latent97_{l}'].value, E[f'wage_slope_latent97_{l}'].se) for l in ("pooled", "male", "female")) + " |")
    add("| observed slope | " + " | ".join(fmt(E[f'wage_slope_latent97_{l}'].extra['obs']['value'], E[f'wage_slope_latent97_{l}'].extra['obs']['se']) for l in ("pooled", "male", "female")) + " |")
    add("| N person-years | " + " | ".join(str(E[f'wage_slope_latent97_{l}'].n) for l in ("pooled", "male", "female")) + " |")
    add("\n| subtest | loading | posterior-variance reliability |\n|---|---|---|")
    for t in TESTS79:
        add(f"| {t} | {e.extra['lam'][t]['value']:.3f} | {e.extra['postvar_reliability'][t]['value']:.3f} |")

    add("\n## E. Earnings checks (not targets)\n")
    e = E["E1.corr"]
    add("Permanent earnings = mean over ages 25-34 of log wage-and-salary income residualized on age and year "
        "dummies within generation (and, for young adults, sex), after 1/99 trimming within year.\n")
    add("| statistic | estimate (SE) | n |\n|---|---|---|")
    add(f"| E1 mother-child correlation | {fmt(e.value, e.se)} | {e.n} |")
    add(f"| E1 covariance | {fmt(e.extra['cov']['value'], e.extra['cov']['se'])} | |")
    add(f"| E1 slope (child on mother) | {fmt(e.extra['slope']['value'], e.extra['slope']['se'])} | |")
    add(f"| E1 correlation, male / female children | {fmt(e.extra['male_children']['value'], e.extra['male_children']['se'])} / "
        f"{fmt(e.extra['female_children']['value'], e.extra['female_children']['se'])} | {info['n.E1_male']} / {info['n.E1_female']} |")
    add(f"| E1 ≥ 2 years each: correlation, slope | {fmt(e.extra['at_least_2_years']['corr']['value'], e.extra['at_least_2_years']['corr']['se'])}, "
        f"{fmt(e.extra['at_least_2_years']['slope']['value'], e.extra['at_least_2_years']['slope']['se'])} | |")
    add(f"| E2 sibling correlation | {fmt(E['E2.corr_sib'].value, E['E2.corr_sib'].se)} | {E['E2.corr_sib'].n} |")
    add(f"| E3 implied rho (E2 / E1) | {fmt(E['E3.implied_rho'].value, E['E3.implied_rho'].se)} | |")

    add("\n## Definitions for T2b\n")
    sd_w = bt.val("C.sd_afqt_wage_sample")
    add("In every case the model panel must use the same sample, ages, weights and composites. 'Score signal' means "
        "`θ0 + θ1 (log z + c log Q_l)` without the noise u; var(u) comes from the R_comp's.\n")
    add("- **`rho_pc_latent`** (ρz). Data: mothers' scores are the 1980 ASVAB taken at ages 15-23 (four IRT z-scores, "
        "age-normed by NLS within birth cohort); children's scores are PIAT math, recognition, comprehension and PPVT-R "
        "at ages 5-14, each age-normed within 3-month cells and averaged over the child's assessment rounds "
        f"(about {info['B.mean_assess']:.1f} rounds per child, so the child factor is a multi-round average; imputed "
        "comprehension scores excluded). Estimate: two-factor ULS on the 8x8 weighted pairwise-complete correlation matrix, "
        "weights = mother's 1979 weight / n_c. Model counterpart: the correlation between the latent score signals of a "
        "parent and a child, one child per parent, parents weighted by population mass. Includes assortative mating and "
        "the Q_l channel because the data factor does too.")
    add(f"- **`wage_slope_latent_pooled`** (b·s_z). `afqt_c` = equal-weight mean of the four standardized IRT z-scores, "
        "re-standardized (weighted) in the Block A POPULATION sample (all persons with four valid scores, any labor-market "
        f"state; NOT the wage sample). Pooled person-years at ages {AGE_LO}-{AGE_HI} at interview, with income, weeks and "
        f"hours for the previous calendar year; filters income > 0, weeks ≥ {FY_WEEKS}, hours/weeks ≥ {FT_HOURS}, "
        f"ESR ≠ 4, highest grade ≥ {MIN_HGC}; weighted 1st/99th percentile trim of log hourly wage within survey year; "
        "weighted OLS on afqt_c with female, age and year dummies; slope divided by sqrt(R_comp_pop). Model counterpart: "
        "the same regression on the noisy simulated composite standardized in the simulated population, divided by the "
        "square root of its reliability, or equivalently the slope on the noise-free signal standardized in the population. "
        f"Standardizing in the wage sample instead would move the slope by about {abs(1 / sd_w - 1) * 100:.0f}% "
        f"(SD of afqt_c in the full-year full-time sample is {sd_w:.2f}). `implied_s_z` = slope × (1 − η).")
    add("- **`R_comp_mother`, `R_comp_pop`**: (Σλ)² / (1'S1) for the equal-weight composite of the four standardized IRT "
        "z-scores, S the sample correlation matrix, in the dyad sample and the population sample respectively. "
        "**`R_comp_child`**: reliability of the composite actually used, the mean of each child's available standardized "
        "tests, R = Σ_i w_i (mean_{k∈K_i} λ_k)² / Σ_i w_i (1'S_{K_i}1 / |K_i|²) over dyad children (`R_comp_child_all4` "
        "in the JSON assumes all four tests present). The model composite's noise variance is set so that the simulated "
        "composite's reliability matches these.")
    add("- **`rho_sib_latent`**: Σ_{k≠l} C_kl λ_k λ_l / Σ_{k≠l} (λ_k λ_l)² over ordered sibling pairs (mother-level "
        "weights); a check, expected ≥ ρ² but sibling-shared inputs are outside the model.")
    add("- **Earnings checks (E1-E3)**: permanent earnings = mean of age-and-year-residualized log wage-and-salary income "
        "over ages 25-34 per person; correlations weighted with mother's weight / (number of children with earnings).")
    add("- **Age norming**: data scores are age-normed to N(0,1) before averaging. In the model, scores are already "
        "age-free.")
    (outdir / "nlsy_ability.md").write_text("\n".join(L) + "\n")
    LOG.info("wrote %s and %s", outdir / "nlsy_ability.json", outdir / "nlsy_ability.md")


# =============================================================================
# Self-test
# =============================================================================

def selftest() -> int:
    rng = np.random.default_rng(1)
    fails = 0

    def check(name: str, ok: bool, detail: str) -> None:
        nonlocal fails
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
        fails += (not ok)

    # (1) two-factor ULS recovers a known structure
    lm = np.array([.8, .85, .7, .75])
    lc = np.array([.7, .8, .75, .6])
    rho = 0.45
    N = 5000
    fm = rng.standard_normal(N)
    fc = rho * fm + np.sqrt(1 - rho ** 2) * rng.standard_normal(N)
    XM = lm * fm[:, None] + np.sqrt(1 - lm ** 2) * rng.standard_normal((N, 4))
    XC = lc * fc[:, None] + np.sqrt(1 - lc ** 2) * rng.standard_normal((N, 4))
    X = np.hstack([XM, XC])
    X[rng.random(X.shape) < 0.15] = np.nan
    X[:, :4] = XM
    w = rng.uniform(0.2, 3.0, N)
    R = wcorr_cross(X, X, w)
    fit = fit_factor_uls(R, GROUPS_PC)
    err = np.abs(fit.lam - np.concatenate([lm, lc])).max()
    check("fit_factor_uls", abs(fit.rho - rho) < 0.03 and err < 0.03 and fit.converged,
          f"rho={fit.rho:.3f} (true .45), max|lam err|={err:.3f}, rmsr={fit.rmsr:.4f}")
    # (2) rank_normal
    n = 20000
    cell = rng.integers(0, 5, n)
    x = rng.standard_normal(n) * (1 + cell) + cell
    ww = rng.uniform(0.5, 2, n)
    z = rank_normal(x, ww, cell)
    mu, sd = wmean(z, ww), np.sqrt(wvar(z, ww))
    check("rank_normal", abs(mu) < 0.02 and abs(sd - 1) < 0.03, f"mean={mu:.4f}, sd={sd:.4f}")
    # (3) disattenuation
    n = 400000
    f = rng.standard_normal(n)
    e = rng.standard_normal(n) * 0.5           # var(e)=.25 -> R = 1/1.25 = .8
    xo = f + e
    y = 0.2 * f + 0.3 * rng.standard_normal(n)
    xs = (xo - xo.mean()) / xo.std()
    Rt = 1 / 1.25
    b = wls(y, np.column_stack([np.ones(n), xs]), np.ones(n))[1]
    check("disattenuation", abs(b / np.sqrt(Rt) - 0.2) < 0.01, f"b_obs={b:.4f}, latent={b / np.sqrt(Rt):.4f} (true .2)")
    # composite_reliability sanity: equal loadings .8 x 4 -> R_comp = 16*.64/(4+12*.64)
    S = np.full((4, 4), .64) + np.eye(4) * .36
    rc, _ = comp_reliability(np.full(4, .8), S)
    check("comp_reliability", abs(rc - 10.24 / 11.68) < 1e-9, f"{rc:.4f}")
    # phi_inv agrees with NormalDist
    check("phi_inv", abs(phi_inv(np.array([0.975]))[0] - 1.959964) < 1e-5, "Phi^-1(.975)")
    print("selftest:", "ALL PASS" if not fails else f"{fails} FAILED")
    return 1 if fails else 0


# =============================================================================
# CLI
# =============================================================================

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bootstrap", type=int, default=200, help="bootstrap draws B (default 200)")
    p.add_argument("--seed", type=int, default=20260928)
    p.add_argument("--outdir", default=str(HERE / "estimates"))
    p.add_argument("--quick", action="store_true", help="B = 20")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    if args.selftest:
        return selftest()
    B = 20 if args.quick else args.bootstrap
    t0 = time.time()
    LOG.info("== loading ==")
    c79, ccy, c97 = (Cohort(n, HERE) for n in ("nlsy79", "nlscya", "nlsy97"))
    inputs = {f"{n}_ability": {"rows": int(c.df.shape[0]), "cols": int(c.df.shape[1])}
              for n, c in (("nlsy79", c79), ("nlscya", ccy), ("nlsy97", c97))}
    LOG.info("== building ==")
    D = build_n79(c79)
    D.kids = build_kids(ccy, D)
    D.ya = build_ya(ccy, D, D.kids)
    build_sets(D)
    R97 = build_n97(c97)
    LOG.info("== point estimates ==")
    info: dict[str, Any] = {}
    point = run_all(D, R97, None, None, info)
    LOG.info("point estimates done in %.1fs", time.time() - t0)
    draws, times = ({}, [])
    if B > 0:
        LOG.info("== bootstrap B=%d ==", B)
        draws, times = bootstrap(D, R97, B, args.seed)
    bt = Boot(point, draws)
    reg = Registry()
    make_records(reg, bt, info, args)
    samples = {k: v for k, v in info.items() if not k.startswith(("fit.", "R."))}
    nfail = int(sum(1 for r in range(B) if not np.isfinite(draws.get("rho_pc_latent", np.full(B, np.nan))[r])))
    meta = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "script": "data/spatial/nlsy_ability.py",
        "inputs": inputs, "samples": samples,
        "assumptions": {
            "eta": ETA, "chi": CHI, "B": B, "seed": args.seed, "failed_draws": nfail,
            "age_window": [AGE_LO, AGE_HI], "min_hgc": MIN_HGC, "fy_weeks": FY_WEEKS,
            "ft_hours_per_week": FT_HOURS, "trim": list(TRIM),
            "weights": {"nlsy79": "1979 SAMPWEIGHT/100", "cnlsy_dyad": "mother's weight / n_c",
                        "nlsy97": "1997 SAMPLING_WEIGHT_CC/100"},
            "clusters": {"nlsy79_cnlsy": "HHID (multinomial household multiplicities)",
                         "nlsy97": "PUBID (no household id in extract; sibling clustering ignored)"},
            "norming": {"nlsy79_irt": "NLS four-month age groups, as delivered",
                        "children": f"rank-normal within 3-month age cells, window {CH_WINDOW} months, fixed across draws",
                        "nlsy97": "rank-normal within birth year x quarter, fixed across draws"},
            "child_window_months": list(CH_WINDOW),
            "bootstrap_flags": {k: (float(np.nanmean(draws[k])) if k in draws else None) for k in
                                ("A.heywood", "B.heywood", "V.heywood")},
            "nonconverged_fits": {"point_estimate": sum(1 for k, v in point.items()
                                                        if k.startswith("conv.") and v < 1),
                                  "n_fits_per_draw": sum(1 for k in point if k.startswith("conv.")),
                                  "draws_with_any": int(sum(
                                      any(draws[k][b] < 1 for k in draws if k.startswith("conv."))
                                      for b in range(B))) if B else 0,
                                  "total_over_draws": int(sum(
                                      np.nansum(draws[k] < 1) for k in draws if k.startswith("conv.")))},
            "mother_sample_ids": list(MOTHER_SAMPLES), "population_sample_ids": list(POP_SAMPLES),
        },
    }
    write_reports(reg, meta, Path(args.outdir).expanduser().resolve(), bt)
    LOG.info("done in %.1fs (%d estimates, %d sensitivities, %d failed draws)", time.time() - t0,
             len(reg.items), len(reg_s.items), nfail)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
