#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_ba_slopes.py -- NLSY79: the latent wage-AFQT slope holding a BA indicator fixed, the same slope
for K-12 teachers and for non-teacher BA workers (and the teacher pay curvature gamma it implies), and
the logit shift that takes the NLSY79 BA share to the ACS non-teacher BA share (task D3b of the
2026-10b round, DESIGN_v2.md section M3).

Why.  M3 lets the young choose a BA before the occupational shocks.  The ability scale s_z can then no
longer be read off the unconditional NLSY79 wage-AFQT slope 0.2127 (SE 0.0077) of `nlsy_ability.py`
(part of that slope is the BA premium and the selection into college); it is re-estimated jointly with
the college parameters against the slope holding the BA status fixed.  The model's ratio of teachers'
to non-teachers' log-wage slopes on log z, gamma (1 - eta) / (1 - gamma eta), checks the teacher wage
curvature gamma.  This script supplies those data moments.

Statistics (all NLSY79, ages 25-34, SAME wage sample, wage definition, weights, latent-AFQT
standardization and household bootstrap as `nlsy_ability.py` block C; the script first reproduces the
published slope 0.2127):

  1. the latent wage-AFQT slope controlling for a BA indicator, pooled (with a female dummy) and by sex;
     with the BA indicator interacted (slopes within BA and within non-BA); the BA coefficient itself.
  2. the same slope (the BA-controlled one of part 1) among K-12 teachers (the `nlsy_teachers.py`
     definition: Census 1980 codes, all K-12 codes, a teacher in the majority of the employed 25-34
     rounds) and among non-teacher BA workers (for whom the BA control is constant); their ratio, and the
     gamma that solves ratio = gamma (1 - eta) / (1 - gamma eta) at eta = 0.080.  Beside it: the teachers'
     slope without the BA control, among BA teachers only, and against all non-teachers.
  3. the logit level shift that moves the NLSY79 BA share (primary sample of `education_occupation.py`)
     to the ACS 2009-13 share of the SAME population, and the implied shares by tercile of the AFQT:
     all persons (teachers included) to the all-persons ACS share 0.3391, and the non-teachers of the
     sample to the ACS non-teacher share 0.3231.  The mixed version (all persons to the non-teacher
     share, 0.440) is kept as a flagged alternative: the NLSY79 all-persons share contains the teachers,
     about nine in ten of whom have a BA.
     P(BA | AFQT tercile) and the BA-non-BA latent gap themselves are NOT recomputed for use: they are
     those of `education_occupation.json` (D3); the script only reproduces them as a check.

Definitions.
  Wage sample: person-years at 25-34 with the ACS cutoffs (>= 48 weeks, >= 30 hours/week, wage income
    >= $1,000 in 2010 dollars), highest grade >= 9, not armed forces; log hourly wage = log(income / annual
    hours); no trimming (`nlsy_ability.PY79.fy`, `.lw`).  Weight: the person's 1979 SAMPWEIGHT x household
    multiplicity, repeated on every person-year (block C).  Controls: age and survey-year dummies and,
    in pooled samples, a female dummy (`nlsy_ability._design`).
  Ability: `afqt_c` (standardized equal-weight composite of the four age-normed ASVAB IRT z-scores,
    block-A population sample, 1979 weights); R = R_comp_pop (0.916; its sampling error is re-estimated
    in every draw).  A slope per SD of the latent factor is the slope per afqt_c SD corrected for
    measurement error; two corrections are reported side by side (plus a check).  The CONVENTION is the
    headline value of every slope in the targets file (it is the pipeline of the 0.2127 target, so the
    BA-conditional slope and the unconditional one are on one definition); EIV is carried as the alternative:
      convention  slope_obs / sqrt(R)   -- `nlsy_ability.py`'s latent slope (0.2127 for the pooled
                  unconditional regression);
      EIV         the classical errors-in-variables correction (Fuller 1987) of the regression itself:
                  afqt_c = f + u with var(u) = 1 - R, var(f) = R; the coefficients solve (M - S) b = m,
                  M, m the weighted moments of the regression and S = (1 - R) E[c c'] the error
                  covariance of the regressors (c = 1 for the AFQT column, BA for AFQT x BA, 0 else);
                  the slope per SD of f is b_AFQT sqrt(R).  With a single mismeasured regressor this is
                  the OLS slope divided by the conditional reliability lam = 1 - (1 - R) / var(afqt_c |
                  controls) times sqrt(R).  For an unconditional regression on the full block-A
                  sample (var(afqt_c) = 1) lam = R and it is the convention; in subsamples (the wage
                  sample, BA holders, teachers: var(afqt_c) 0.5-0.9) and with a control that is
                  correlated with AFQT (BA) the reliability is LOWER than R and the convention
                  understates the latent slope.  The model's counterpart (`ability_moments`) is a
                  regression on the true latent score; both versions estimate that slope, EIV allowing
                  for the lower reliability.  Pairing: the convention BA-conditional slope goes with the
                  convention 0.2127; the EIV one with the EIV unconditional 0.2152 (a PREREG amendment).
                  In a draw whose conditional reliability is below LAM_MIN = 0.5 the EIV statistic is set
                  to missing (it happens only for the 24 male teachers).
      split-half IV (a check, not a third target): one half of the four tests instrumented by the other
                  half (two mixed splits; the math-verbal split is shown but not averaged), times the
                  loading of the composite on the standardized one-factor latent variable; it assumes
                  independent errors across halves instead of an error variance.
  BA: highest grade completed >= 16 at some round at 25-34 (`nlsy_teachers.COLLEGE_HGC`; the project's
    NLSY79 extract has no degree variable, although the survey asks for the highest degree received),
    as `education_occupation.py` (D3).  Sensitivity: >= 16 at the round itself.
  Teachers: `nlsy_teachers.py` person-level classification, spec `c80-emp-all-maj`; a teacher's wage
    person-years are all of his or her wage-sample rounds.  The teachers' slope carries the BA indicator
    (key TC_*), so the credential mix of the teachers (12% without a BA) is held fixed as in part 1; the
    version without it (T_*) is a reference.  Non-teacher BA workers: persons of that spec's non-teacher
    group with a BA.  Sensitivities: core K-12 codes, at-least-half and ever rules, any 25-34 occupation,
    1970 codes, round-level K-12 occupation, and the controls.

Inference.  Household-cluster bootstrap, B = 500, seed 20260929, the draws of `nlsy_teachers.run_bootstrap`
(cluster multiplicities from `rng.integers`; common random numbers with D3 and the teacher gaps); in every
draw AFQT is re-standardized, R is re-estimated and every regression and the logit shift are re-run.  SE =
bootstrap SD; CI = 95% percentile interval.  Reproduction check: the stream of `nlsy_ability.bootstrap`
(multinomial draws, seed 20260928, B = 200, interleaved NLSY97 draw, which is NOT read) is replayed for
the three unconditional slopes and gives the published SEs.

    cd data/spatial && uv run python nlsy_ba_slopes.py              # B = 500
    uv run python nlsy_ba_slopes.py --quick --outdir /tmp/x         # B = 20 (needs an explicit --outdir)
    uv run python nlsy_ba_slopes.py --md-only                       # rewrite .md and .toml from the .json
    uv run python nlsy_ba_slopes.py --selftest                      # unit checks, no data needed

Inputs: raw/nlsy79_ability.parquet (via `nlsy_ability.Cohort`), estimates/nlsy_ability.json (published
slopes), estimates/settings/education_occupation.json (ACS non-teacher BA share, D3 tercile shares).
NLSY97 and the sealed hold-out are not used.  Writes estimates/settings/nlsy_ba_slopes.{json,md} and
nlsy_ba_slopes_targets.toml (no microdata).

Written with help from Claude Code.
"""
from __future__ import annotations

import os

# the linear algebra here is tiny; BLAS/OpenMP pools sized to a shared 72-core node can fail at start-up
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import argparse
import datetime as dt
import json
import logging
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402
import nlsy_extensions as ne  # noqa: E402  (tercile_of)
import nlsy_teachers as nt  # noqa: E402

LOG = logging.getLogger("nlsy_ba_slopes")

SEED = 20260929                 # the draws of nlsy_teachers.run_bootstrap and education_occupation.py
B_DEFAULT = 500
ETA = 0.080                     # goods-investment elasticity (calibration_log.md Table 1; nlsy_ability.ETA = 0.0807)
GAMMA_BASELINE = 0.6194         # baseline fit runs/exact-polish-base (calibration_log.md, exactly identified fit)
COLLEGE_HGC = nt.COLLEGE_HGC    # 16
MAIN_SPEC = "c80-emp-all-maj"
PUB_SEED, PUB_B = 20260928, 200  # nlsy_ability.json: the stream behind the published SEs
GROUPS = ("pooled", "male", "female")
# split-half instrumental-variable check of the disattenuation: the composite of two tests (A) instrumented by the composite of the other two (B)
# and vice versa; tests in the order of nlsy_ability.TESTS79 = AR, WK, PC, MK; two mixed splits (each half has a math and a verbal test) and the math-verbal split
IV_SPLITS = {"awpm": ((0, 1), (2, 3), "AR+WK vs PC+MK"), "apwm": ((0, 2), (1, 3), "AR+PC vs WK+MK"),
             "amwp": ((0, 3), (1, 2), "AR+MK (math) vs WK+PC (verbal)")}
IV_MIXED = ("awpm", "apwm")     # the average `mean4` uses the mixed splits; the math-verbal split is reported only
OUTDIR_DEFAULT = HERE / "estimates" / "settings"
STAT_KEYS = ("value", "se", "ci95_lo", "ci95_hi", "n_draws")
LAM_MIN = 0.5                   # EIV statistics are set to missing when the conditional reliability of the score is below this
                                # (the correction would more than double the slope and is not stable: only the 24 male teachers in a few draws)


# =============================================================================
# Estimators
# =============================================================================

def gamma_of_ratio(r: float, eta: float = ETA) -> float:
    """gamma that solves r = gamma (1 - eta) / (1 - gamma eta)."""
    return float(r / (1.0 - eta + r * eta))


def ratio_of_gamma(g: float, eta: float = ETA) -> float:
    return float(g * (1.0 - eta) / (1.0 - g * eta))


def eiv_fit(y: np.ndarray, X: np.ndarray, w: np.ndarray, ucols: dict[int, np.ndarray], s2u: float
            ) -> tuple[np.ndarray, np.ndarray, float]:
    """Weighted OLS and the classical errors-in-variables corrected coefficients.

    X holds the regressors as observed.  `ucols` maps a column j to the vector c_j such that the error in
    that column is c_ij u_i with var(u) = s2u (1 for a score column, BA_i for score x BA).  The corrected
    coefficients solve (M - S) b = m with M = X'WX / sum(w), m = X'Wy / sum(w), S = s2u E[c c'].  Returns
    (OLS coefficients, corrected coefficients, conditional reliability of the first listed column:
    1 - s2u / var(x_j | the other columns)).  Rank-deficient designs (an empty dummy in a draw) get the
    minimum-norm solution, as `nlsy_ability.wls`.
    """
    W = w / w.sum()
    XW = X * W[:, None]
    M = XW.T @ X
    m = XW.T @ y
    b_ols = np.linalg.lstsq(M, m, rcond=None)[0]
    js = sorted(ucols)
    C = np.column_stack([ucols[j] for j in js])
    G = (C * W[:, None]).T @ C
    S = np.zeros_like(M)
    S[np.ix_(js, js)] = s2u * G
    b_eiv = np.linalg.lstsq(M - S, m, rcond=None)[0]
    j = js[0]
    rest = [i for i in range(X.shape[1]) if i != j]
    mzx = M[rest, j]
    vc = M[j, j] - mzx @ np.linalg.lstsq(M[np.ix_(rest, rest)], mzx, rcond=None)[0]
    return b_ols, b_eiv, float(1.0 - s2u / vc)


def partial_iv(y: np.ndarray, A: np.ndarray, B: np.ndarray, Zc: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    """Weighted IV slopes of y on A instrumented by B and on B instrumented by A, with exogenous controls Zc
    (partialled out): cov(y, B) / cov(A, B) and cov(y, A) / cov(A, B) of the Zc-residuals."""
    W = w / w.sum()
    ZW = Zc * W[:, None]
    V = np.column_stack([y, A, B])
    coef = np.linalg.lstsq(ZW.T @ Zc, ZW.T @ V, rcond=None)[0]
    ry, rA, rB = (V - Zc @ coef).T
    cab = (W * rA * rB).sum()
    return float((W * rB * ry).sum() / cab), float((W * rA * ry).sum() / cab)


def wlogit(x: np.ndarray, y: np.ndarray, w: np.ndarray, iters: int = 30) -> np.ndarray:
    """Weighted logistic regression of y on (1, x) by Newton steps; (a, b).  Same as education_occupation.wlogit."""
    X = np.column_stack([np.ones_like(x), x])
    b = np.zeros(2)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-(X @ b)))
        g = X.T @ (w * (y - p))
        H = X.T @ ((w * p * (1 - p))[:, None] * X)
        step = np.linalg.solve(H, g)
        b += step
        if np.abs(step).max() < 1e-10:
            break
    return b


def logit_shift(eta: np.ndarray, w: np.ndarray, target: float) -> float:
    """delta with sum(w Lambda(eta + delta)) / sum(w) = target (bisection, as education_occupation)."""
    lo, hi = -30.0, 30.0
    ww = w / w.sum()
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (ww / (1.0 + np.exp(-(eta + mid)))).sum() < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def lam_fn(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


# =============================================================================
# Data frame: everything that does not change across draws
# =============================================================================

@dataclass
class Reg:
    """One regression of the log wage on the AFQT (and controls) on a fixed set of person-years."""
    key: str
    rows: np.ndarray                      # bool over person-years, before the draw's w > 0 and finite-AFQT filters
    kind: str = "slope"                   # slope | ba (BA dummy) | inter (BA dummy and AFQT x BA)
    female: bool = False                  # female dummy (pooled samples)
    ba: Optional[np.ndarray] = None       # BA indicator over person-years (kinds ba, inter)
    deflate: bool = False                 # CPI-deflated log wage
    controls: str = "ay"                  # ay = age + year dummies; a = age dummies; none
    desc: str = ""


class Data:
    def __init__(self, fr: nt.Frame) -> None:
        self.fr, D, py = fr, fr.D, fr.py
        self.D, self.py = D, py
        self.female_py = D.female[py.p]
        hg = fr.hg_max
        self.ba_person = np.isfinite(hg) & (hg >= COLLEGE_HGC)
        self.ba_py = self.ba_person[py.p]
        self.ba_round_py = np.isfinite(py.hgc) & (py.hgc >= COLLEGE_HGC)
        self.base = py.fy & np.isfinite(py.lw)
        cpi = np.array([na.CPI_U[int(y) - 1] for y in py.Y]) / na.CPI_U[2010]
        self.lw = py.lw
        self.lw_defl = py.lw - np.log(cpi)
        self.specs = {s.name: s for s in fr.specs}
        self.teacher_person = self.specs[MAIN_SPEC].teacher        # person-level K-12 teacher flag of the main spec (BA-share populations)
        self.regs: list[Reg] = []
        self.ratios: dict[str, tuple[str, str, str]] = {}
        self._build()

    # -- regression registry ---------------------------------------------------------------
    def _gmask(self, g: str) -> np.ndarray:
        return {"pooled": np.ones_like(self.female_py), "male": ~self.female_py, "female": self.female_py}[g]

    def _add(self, **kw: Any) -> None:
        self.regs.append(Reg(**kw))

    def _build(self) -> None:
        b, ba, bar = self.base, self.ba_py, self.ba_round_py
        py, fr = self.py, self.fr
        # ---- 1. unconditional, BA control, interacted, separate samples
        for g in GROUPS:
            gm, fem = self._gmask(g), g == "pooled"
            self._add(key=f"u_{g}", rows=b & gm, female=fem, desc=f"no BA control ({g})")
            self._add(key=f"ba_{g}", rows=b & gm, female=fem, kind="ba", ba=ba, desc=f"BA control ({g})")
            self._add(key=f"bar_{g}", rows=b & gm, female=fem, kind="ba", ba=bar,
                      desc=f"BA control, BA = highest grade >= 16 at the round ({g})")
            self._add(key=f"int_{g}", rows=b & gm, female=fem, kind="inter", ba=ba,
                      desc=f"BA interacted with AFQT ({g})")
            self._add(key=f"sepba_{g}", rows=b & gm & ba, female=fem, desc=f"BA holders only ({g})")
            self._add(key=f"sepnba_{g}", rows=b & gm & ~ba, female=fem, desc=f"non-BA only ({g})")
        # ---- 2. teachers and non-teacher BA workers.  The teachers' slope is the SAME slope as in part 1 (key TC_*:
        # a BA indicator in the regression), so that the 12% of teachers without a BA (lower AFQT, lower pay) do not
        # put the credential composition into the slope the closed form compares with non-teacher BA workers.  T_*
        # (no BA control) is kept as a reference; for BA workers (NTBA_*, TBA_*) the BA indicator is constant.
        sp = self.specs[MAIN_SPEC]
        T, NT = sp.teacher[py.p] & b, sp.nonteacher[py.p] & b
        fm = self.female_py
        self._add(key="TC_main", rows=T, female=True, kind="ba", ba=ba, desc="K-12 teachers (main spec), pooled, BA control")
        self._add(key="TC_women", rows=T & fm, kind="ba", ba=ba, desc="K-12 teachers, women, BA control")
        self._add(key="TC_men", rows=T & ~fm, kind="ba", ba=ba, desc="K-12 teachers, men, BA control")
        self._add(key="T_main", rows=T, female=True, desc="K-12 teachers (main spec), pooled, no BA control")
        self._add(key="T_women", rows=T & fm, desc="K-12 teachers, women, no BA control")
        self._add(key="T_men", rows=T & ~fm, desc="K-12 teachers, men, no BA control")
        self._add(key="TBA_main", rows=T & ba, female=True, desc="K-12 teachers with a BA, pooled")
        self._add(key="NT_main", rows=NT, female=True, desc="all non-teachers, pooled, no BA control")
        self._add(key="NTC_main", rows=NT, female=True, kind="ba", ba=ba, desc="all non-teachers, pooled, BA control")
        self._add(key="NTBA_main", rows=NT & ba, female=True, desc="non-teacher BA workers, pooled")
        self._add(key="NTBA_women", rows=NT & ba & fm, desc="non-teacher BA workers, women")
        self._add(key="NTBA_men", rows=NT & ba & ~fm, desc="non-teacher BA workers, men")
        self._add(key="NTnBA_main", rows=NT & ~ba, female=True, desc="non-teachers without a BA, pooled")
        self.ratios["main"] = ("TC_main", "NTBA_main", "all K-12 teachers (BA control) / non-teacher BA workers (main)")
        self.ratios["nobac"] = ("T_main", "NTBA_main", "all K-12 teachers, no BA control / non-teacher BA workers (first version)")
        self.ratios["tba"] = ("TBA_main", "NTBA_main", "K-12 teachers with a BA / non-teacher BA workers")
        self.ratios["ntc"] = ("TC_main", "NTC_main", "all K-12 teachers / all non-teachers (BA control in both)")
        self.ratios["nt"] = ("T_main", "NT_main", "all K-12 teachers / all non-teachers (no BA control)")
        self.ratios["women"] = ("TC_women", "NTBA_women", "women: teachers (BA control) / non-teacher BA workers")
        self.ratios["men"] = ("TC_men", "NTBA_men", "men: teachers (BA control) / non-teacher BA workers")
        for tag, name, lab in (("core", "c80-emp-core-maj", "core K-12 codes (elementary, secondary, special education)"),
                               ("half", "c80-emp-all-half", "teacher in at least half of the rounds"),
                               ("ever", "c80-emp-all-ever", "ever a teacher at 25-34"),
                               ("any", "c80-any-all-maj", "any 25-34 occupation (not only employed rounds)"),
                               ("c70", "c70-emp-all-maj", "1970 Census codes")):
            s2 = self.specs[name]
            self._add(key=f"TC_{tag}", rows=s2.teacher[py.p] & b, female=True, kind="ba", ba=ba, desc=f"teachers (BA control), {lab}")
            self._add(key=f"NTBA_{tag}", rows=s2.nonteacher[py.p] & b & ba, female=True,
                      desc=f"non-teacher BA workers, {lab}")
            self.ratios[tag] = (f"TC_{tag}", f"NTBA_{tag}", f"teachers (BA control) / non-teacher BA workers, {lab}")
        k12r = fr.k12["c80"]
        self._add(key="TC_round", rows=b & k12r, female=True, kind="ba", ba=ba, desc="person-years in a K-12 occupation at the round (BA control)")
        self._add(key="NTBA_round", rows=b & fr.valid["c80"] & ~k12r & ba, female=True,
                  desc="BA person-years in another occupation at the round")
        self.ratios["round"] = ("TC_round", "NTBA_round", "round-level K-12 occupation (BA control) / other BA person-years")
        for tag, ctrl, lab in (("deflage", "a", "CPI-deflated log wage, age dummies and female (no year dummies)"),
                               ("deflnone", "none", "CPI-deflated log wage, female dummy only")):
            self._add(key=f"TC_{tag}", rows=T, female=True, kind="ba", ba=ba, deflate=True, controls=ctrl, desc=f"teachers (BA control), {lab}")
            self._add(key=f"NTBA_{tag}", rows=NT & ba, female=True, deflate=True, controls=ctrl,
                      desc=f"non-teacher BA workers, {lab}")
            self.ratios[tag] = (f"TC_{tag}", f"NTBA_{tag}", f"teachers (BA control) / non-teacher BA workers, {lab}")

    # -- design ----------------------------------------------------------------------------
    def design(self, rg: Reg, rows: np.ndarray, a: np.ndarray) -> tuple[np.ndarray, dict[int, np.ndarray]]:
        """X = [1, afqt, (BA, AFQT x BA), (female), age dummies, year dummies] and the error multipliers."""
        n = int(rows.sum())
        cols = [np.ones(n), a]
        uc = {1: np.ones(n)}
        if rg.kind in ("ba", "inter"):
            bv = rg.ba[rows].astype(float)
            cols.append(bv)
            if rg.kind == "inter":
                cols.append(a * bv)
                uc[len(cols) - 1] = bv
        if rg.female:
            cols.append(self.female_py[rows].astype(float))
        X = np.column_stack(cols)
        blocks = [X]
        if rg.controls in ("ay", "a"):
            blocks.append(self.py.age_d[rows])
        if rg.controls == "ay":
            blocks.append(self.py.yr_d[rows])
        return np.column_stack(blocks), uc


def fit_reg(dat: Data, rg: Reg, a_py: np.ndarray, w_py: np.ndarray, s2u: float, sqR: float,
            counts: Optional[dict[str, Any]] = None) -> dict[str, float]:
    """Statistics of one regression for one draw (keys `{rg.key}|{stat}`)."""
    rows = rg.rows & (w_py > 0) & np.isfinite(a_py)
    k = rg.key
    if counts is not None:
        counts[k] = {"person_years": int(rows.sum()), "persons": int(np.unique(dat.py.p[rows]).size), "desc": rg.desc}
    if rows.sum() < 30:
        return {}
    y = (dat.lw_defl if rg.deflate else dat.lw)[rows]
    X, uc = dat.design(rg, rows, a_py[rows])
    bo, be, lam = eiv_fit(y, X, w_py[rows], uc, s2u)
    if rg.kind != "inter" and lam < LAM_MIN:
        be = np.full_like(be, np.nan)
    out: dict[str, float] = {}
    if rg.kind == "inter":
        for tag, jj in (("nonba", [1]), ("ba", [1, 3])):
            out[f"{k}|{tag}_obs"] = float(bo[jj].sum())
            out[f"{k}|{tag}_naive"] = float(bo[jj].sum() / sqR)
            out[f"{k}|{tag}_eiv"] = float(be[jj].sum() * sqR)
        out[f"{k}|diff_naive"] = out[f"{k}|ba_naive"] - out[f"{k}|nonba_naive"]
        out[f"{k}|diff_eiv"] = out[f"{k}|ba_eiv"] - out[f"{k}|nonba_eiv"]
        out[f"{k}|ba_coef_naive"], out[f"{k}|ba_coef_eiv"] = float(bo[2]), float(be[2])
    else:
        out[f"{k}|obs"], out[f"{k}|naive"], out[f"{k}|eiv"] = float(bo[1]), float(bo[1] / sqR), float(be[1] * sqR)
        out[f"{k}|lam"] = lam
        if rg.kind == "ba":
            out[f"{k}|ba_coef_naive"], out[f"{k}|ba_coef_eiv"] = float(bo[2]), float(be[2])
    return out


# =============================================================================
# One draw
# =============================================================================

def iv_block(dat: Data, w: np.ndarray, lam_tests: np.ndarray, w_py: np.ndarray) -> dict[str, float]:
    """Split-half IV slopes per SD of the latent factor: A <- B scaled by the loading of composite A on the
    standardized factor (mean test loading of the one-factor ULS fit over the composite's own SD), B <- A by
    that of B.  Pooled wage sample; controls as the main regressions (female dummy, age and year dummies and,
    for `ba`, the BA indicator).  Keys iv_{split}_{ctrl}|AB, |BA, |mean and iv_{ctrl}|mean4 (the two mixed splits, both directions)."""
    D, py = dat.D, dat.py
    m = D.asamp
    Zs = np.column_stack([na.standardize(D.z[m][:, k], w[m]) for k in range(4)])
    out: dict[str, float] = {}
    rows = dat.base & (w_py > 0)
    comp: dict[str, tuple[np.ndarray, float]] = {}
    for tag, (ia, ib, _) in IV_SPLITS.items():
        for nm, idx in (("A", list(ia)), ("B", list(ib))):
            raw = Zs[:, idx].mean(axis=1)
            sd = float(np.sqrt(na.wvar(raw, w[m])))
            full = np.full(D.n, np.nan)
            full[m] = (raw - na.wmean(raw, w[m])) / sd
            comp[f"{tag}{nm}"] = (full[py.p], float(lam_tests[idx].mean() / sd))
    for ctrl in ("none", "ba"):
        vals = []
        for tag in IV_SPLITS:
            (A, lA), (B, lB) = comp[f"{tag}A"], comp[f"{tag}B"]
            r = rows & np.isfinite(A) & np.isfinite(B)
            cols = [np.ones(int(r.sum()))] + ([dat.ba_py[r].astype(float)] if ctrl == "ba" else []) + [dat.female_py[r].astype(float), py.age_d[r], py.yr_d[r]]
            b_ab, b_ba = partial_iv(dat.lw[r], A[r], B[r], np.column_stack(cols), w_py[r])
            out[f"iv_{tag}_{ctrl}|AB"], out[f"iv_{tag}_{ctrl}|BA"] = b_ab * lA, b_ba * lB
            out[f"iv_{tag}_{ctrl}|mean"] = 0.5 * (b_ab * lA + b_ba * lB)
            if tag in IV_MIXED:
                vals += [b_ab * lA, b_ba * lB]
        out[f"iv_{ctrl}|mean4"] = float(np.mean(vals))
    return out


def ba_pop_block(prefix: str, afqt: np.ndarray, w: np.ndarray, mask: np.ndarray, ba: np.ndarray, t: np.ndarray,
                 targets: dict[str, float]) -> dict[str, float]:
    """The BA-share statistics of the persons in `mask`, the AFQT terciles `t` being fixed by the caller (so that
    every population is cut at the same observed-score values): the weighted BA share, the BA share (terc{q}) and the
    population share (pi{q}) of each tercile, the individual-level logit of BA on the AFQT, and for every target tag
    the level shifts that take the weighted BA share to the target -- `shift_cell` (one delta on the logits of the
    tercile shares: sum_k pi_k Lambda(logit s_k + delta) = target) and `shift_indiv` (one delta on the individual
    logit: D3's 'level shift of the same gradient') -- with the implied tercile shares of both."""
    out: dict[str, float] = {}
    wm_, tm_, xm_, bm_ = w[mask], t[mask], afqt[mask], ba[mask]
    out[f"{prefix}|share"] = float(wm_[bm_].sum() / wm_.sum())
    s_k, pi_k = [], []
    for q in range(3):
        mq = tm_ == q
        s_k.append(float(wm_[mq & bm_].sum() / wm_[mq].sum()))
        pi_k.append(float(wm_[mq].sum() / wm_.sum()))
        out[f"{prefix}|terc{q + 1}"] = s_k[-1]
        out[f"{prefix}|pi{q + 1}"] = pi_k[-1]
    out[f"{prefix}|terc_diff"] = s_k[2] - s_k[0]
    a_, b_ = wlogit(xm_, bm_.astype(float), wm_)
    out[f"{prefix}|logit_a"], out[f"{prefix}|logit_b"] = float(a_), float(b_)
    eta_i = a_ + b_ * xm_
    eta_k = np.log(np.array(s_k) / (1.0 - np.array(s_k)))
    p0 = lam_fn(eta_i)                                  # the logit fit itself, no shift
    for q in range(3):
        mq = tm_ == q
        out[f"{prefix}|fit_indiv_terc{q + 1}"] = float((wm_[mq] * p0[mq]).sum() / wm_[mq].sum())
    pi_arr = np.array(pi_k)
    for tag, tgt in targets.items():
        d_i = logit_shift(eta_i, wm_, tgt)
        d_c = logit_shift(eta_k, pi_arr, tgt)
        out[f"{prefix}|shift_indiv_{tag}"], out[f"{prefix}|shift_cell_{tag}"] = d_i, d_c
        p_i = lam_fn(eta_i + d_i)
        p_c = lam_fn(eta_k + d_c)
        # d(delta)/d(target) = 1 / E[p (1 - p)]: carries the ACS sampling error of the target into the shift
        out[f"{prefix}|dshift_indiv_{tag}"] = float(1.0 / ((wm_ / wm_.sum()) * p_i * (1 - p_i)).sum())
        out[f"{prefix}|dshift_cell_{tag}"] = float(1.0 / (pi_arr * p_c * (1 - p_c)).sum())
        imp_i = []
        for q in range(3):
            mq = tm_ == q
            imp_i.append(float((wm_[mq] * p_i[mq]).sum() / wm_[mq].sum()))
            out[f"{prefix}|imp_indiv_{tag}_terc{q + 1}"] = imp_i[-1]
            out[f"{prefix}|imp_cell_{tag}_terc{q + 1}"] = float(p_c[q])
        out[f"{prefix}|imp_indiv_{tag}_diff"] = imp_i[2] - imp_i[0]
        out[f"{prefix}|imp_cell_{tag}_diff"] = float(p_c[2] - p_c[0])
        out[f"{prefix}|imp_indiv_{tag}_share"] = float((wm_ * p_i).sum() / wm_.sum())
    return out


def ba_block(dat: Data, afqt: np.ndarray, w: np.ndarray, targets_all: dict[str, float],
             targets_nt: dict[str, float]) -> dict[str, float]:
    """BA share by AFQT tercile (D3 definitions) of two populations of D3's primary sample: ALL persons (`ba|...`; teachers
    included) and the NON-TEACHERS among them (`banT|...`; everyone not classified as a K-12 teacher in the main spec,
    persons never employed at 25-34 included, as the ACS non-teacher universe includes home production).  The terciles
    are those of all persons (the same cut points for both).  Each population is shifted to its own targets."""
    fr, D = dat.fr, dat.D
    hg = fr.hg_max
    okg = np.isfinite(hg)
    ba = okg & (hg >= COLLEGE_HGC)
    m = fr.observed & D.asamp & okg & (hg >= na.MIN_HGC) & np.isfinite(afqt) & (w > 0)   # D3's primary sample
    t, cuts = ne.tercile_of(np.where(m, afqt, np.nan), w * m)
    out: dict[str, float] = {"ba|cut1": float(cuts[0]), "ba|cut2": float(cuts[1])}
    out.update(ba_pop_block("ba", afqt, w, m, ba, t, targets_all))
    teach = dat.teacher_person
    out.update(ba_pop_block("banT", afqt, w, m & ~teach, ba, t, targets_nt))
    out["ba|teacher_wshare"] = float(w[m & teach].sum() / w[m].sum())
    out["ba|teacher_ba_share"] = float(w[m & teach & ba].sum() / w[m & teach].sum())
    return out


def evaluate(dat: Data, mult_hh: Optional[np.ndarray], targets: dict[str, dict[str, float]],
             counts: Optional[dict[str, Any]] = None) -> dict[str, float]:
    """All statistics for household multiplicities `mult_hh` (None = the sample itself).  `targets` = {"all": the
    ACS BA shares that the all-persons NLSY79 shares are shifted to, "nt": those for the non-teachers}."""
    D, py = dat.D, dat.py
    mult_p = np.ones(D.n) if mult_hh is None else mult_hh[D.hh_idx].astype(float)
    w = D.w * mult_p
    afqt, Rc, _, fit, _ = na.pop_afqt(D.z, D.asamp, w)
    sqR, s2u = float(np.sqrt(Rc)), float(1.0 - Rc)
    a_py, w_py = afqt[py.p], w[py.p]
    out: dict[str, float] = {"R_comp_pop": float(Rc)}
    for rg in dat.regs:
        out.update(fit_reg(dat, rg, a_py, w_py, s2u, sqR, counts))
    nan = float("nan")
    for tag, (num, den, _) in dat.ratios.items():
        for sc in ("naive", "eiv"):
            a, b = out.get(f"{num}|{sc}", nan), out.get(f"{den}|{sc}", nan)
            r = a / b if np.isfinite(a) and np.isfinite(b) and b > 0 else nan
            out[f"r_{tag}|{sc}"] = r
            out[f"gamma_{tag}|{sc}"] = gamma_of_ratio(r) if np.isfinite(r) else nan
    for g in GROUPS:
        for sc in ("naive", "eiv"):
            a, b = out.get(f"ba_{g}|{sc}", nan), out.get(f"u_{g}|{sc}", nan)
            out[f"ratio_ba_u_{g}|{sc}"] = a / b if np.isfinite(a) and np.isfinite(b) and b > 0 else nan
    out.update(iv_block(dat, w, np.asarray(fit.lam), w_py))
    out.update(ba_block(dat, afqt, w, targets["all"], targets["nt"]))
    return out


def run_bootstrap(dat: Data, B: int, seed: int, targets: dict[str, dict[str, float]]) -> dict[str, Any]:
    counts: dict[str, Any] = {}
    point = evaluate(dat, None, targets, counts)
    keys = list(point)
    rng = np.random.default_rng(seed)
    draws = np.full((B, len(keys)), np.nan)
    nhh = dat.D.nhh
    t0 = time.time()
    for b in range(B):
        cm = np.bincount(rng.integers(0, nhh, nhh), minlength=nhh).astype(float)     # the draws of nlsy_teachers.run_bootstrap
        try:
            d = evaluate(dat, cm, targets)
            draws[b] = [d.get(k, np.nan) for k in keys]
        except Exception as exc:                                                     # noqa: BLE001 -- a failed draw is a NaN row
            LOG.warning("  draw %d failed: %s", b, exc)
        if (b + 1) % 25 == 0 or b + 1 == B:
            LOG.info("  bootstrap %d/%d (%.0fs, %.2fs/draw)", b + 1, B, time.time() - t0, (time.time() - t0) / (b + 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        se = np.nanstd(draws, axis=0, ddof=1) if B > 1 else np.full(len(keys), np.nan)
        lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0) if B > 1 else (np.full(len(keys), np.nan),) * 2
    n_ok = np.sum(np.isfinite(draws), axis=0)
    stats = {k: {"value": float(point[k]), "se": float(se[i]), "ci95_lo": float(lo[i]), "ci95_hi": float(hi[i]),
                 "n_draws": int(n_ok[i])} for i, k in enumerate(keys)}
    return {"stats": stats, "counts": counts, "B": B, "seed": seed, "draws": draws, "keys": keys,
            "failed_draws": int(np.sum(~np.isfinite(draws[:, :1]).any(axis=1))) if B else 0}


def delta_method_gamma(res: dict[str, Any], num: str, den: str, sc: str = "eiv") -> dict[str, float]:
    """Delta-method SE of ratio = a / b and gamma(ratio) from the bootstrap covariance of the two slopes."""
    keys = res["keys"]
    ia, ib = keys.index(f"{num}|{sc}"), keys.index(f"{den}|{sc}")
    dr = res["draws"][:, [ia, ib]]
    dr = dr[np.isfinite(dr).all(axis=1)]
    a, b = res["stats"][f"{num}|{sc}"]["value"], res["stats"][f"{den}|{sc}"]["value"]
    r = a / b
    if dr.shape[0] < 10 or not np.isfinite(r):
        nan = float("nan")
        return {"ratio": r, "se_ratio": nan, "gamma": nan, "se_gamma": nan, "dgamma_dratio": nan, "corr_slopes": nan, "n_draws": int(dr.shape[0])}
    cov = np.cov(dr.T, ddof=1)
    grad = np.array([1.0 / b, -a / b ** 2])
    se_r = float(np.sqrt(grad @ cov @ grad))
    dg = (1.0 - ETA) / (1.0 - ETA + r * ETA) ** 2
    return {"ratio": r, "se_ratio": se_r, "gamma": gamma_of_ratio(r), "se_gamma": float(dg * se_r),
            "dgamma_dratio": float(dg), "corr_slopes": float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])),
            "n_draws": int(dr.shape[0])}


# =============================================================================
# Reproduction of the published slope and its SE
# =============================================================================

def replay_published(dat: Data, B: int, seed: int, n97: int) -> dict[str, Any]:
    """Replay the bootstrap stream of nlsy_ability.bootstrap (multinomial household multiplicities for the
    NLSY79 and, interleaved, for the NLSY97, whose draw is generated but not used) for the unconditional
    slopes of block C, using `nlsy_ability.slope_set` itself."""
    D, py = dat.D, dat.py
    rng = np.random.default_rng(seed)
    female = D.female[py.p]
    vals: dict[str, list[float]] = {g: [] for g in GROUPS}
    t0 = time.time()
    for b in range(B):
        mh = rng.multinomial(D.nhh, np.full(D.nhh, 1.0 / D.nhh))
        rng.multinomial(n97, np.full(n97, 1.0 / n97))
        w = D.w * mh[D.hh_idx].astype(float)
        afqt, Rc = na.pop_afqt(D.z, D.asamp, w)[:2]
        out: dict[str, float] = {}
        na.slope_set(py, py.lw, w[py.p], afqt[py.p], female, py.fy, np.sqrt(Rc), "", out, None)
        for g in GROUPS:
            vals[g].append(out[f"wage_slope_latent_{g}"])
        if (b + 1) % 50 == 0:
            LOG.info("  replay %d/%d (%.0fs)", b + 1, B, time.time() - t0)
    return {g: {"se": float(np.std(v, ddof=1)), "mean": float(np.mean(v)), "B": B, "seed": seed} for g, v in vals.items()}


def read_published() -> dict[str, Any]:
    j = json.loads((HERE / "estimates" / "nlsy_ability.json").read_text())
    est = {e["key"]: e for e in j["estimates"]}
    return {"inputs": j["inputs"], "assumptions": {k: j["assumptions"][k] for k in ("B", "seed", "eta")},
            "slopes": {g: {"value": est[f"wage_slope_latent_{g}"]["value"], "se": est[f"wage_slope_latent_{g}"]["se"],
                           "n": est[f"wage_slope_latent_{g}"]["n"]} for g in GROUPS},
            "R_comp_pop": est["R_comp_pop"]["value"], "implied_s_z": est["implied_s_z"]["value"]}


def read_d3() -> dict[str, Any]:
    j = json.loads((OUTDIR_DEFAULT / "education_occupation.json").read_text())
    s = j["nlsy"]["stats"]
    T = j["targets"]
    return {"ba_share_nonteacher": T["ba_share_nonteacher"]["value"],
            "ba_share_nonteacher_se": T["ba_share_nonteacher"]["se"],
            "ba_share_all": j["nlsy"]["meta"]["share_target"]["pooled"],
            "ba_share_all_se": T["ba_share_all"]["se"],
            "acs_teacher_share": T["teacher_ba_share"]["n_universe"] / T["ba_share_all"]["n_universe"],
            "teacher_gap_vs_ba_nonteachers": T["teacher_z_gap_nonteacher_ba"]["value"], "teacher_gap_within_ba": T["teacher_z_gap_within_ba"]["value"],
            "stats": {k: s[k] for k in s if k.startswith(("terc|all_hgc9|pooled", "logit|all_hgc9|pooled", "all_hgc9|pooled|ba_share",
                                                           "all_hgc9|pooled|gap", "terc|all_hgc9|cut"))},
            "generated": j["meta"]["generated"]}


# =============================================================================
# Report
# =============================================================================

def finite(x: Any) -> bool:
    return x is not None and bool(np.isfinite(x))


def clean(obj: Any) -> Any:
    """JSON-safe copy: non-finite floats become null."""
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


class Rep:
    """Formatting helpers over the saved payload."""

    def __init__(self, P: dict[str, Any]) -> None:
        self.P, self.S = P, P["statistics"]

    def r(self, key: str) -> Optional[dict[str, float]]:
        r = self.S.get(key)
        return r if r is not None and finite(r.get("value")) else None

    def v(self, key: str, d: int = 3) -> str:
        r = self.r(key)
        return "n/a" if r is None else f"{r['value']:.{d}f}"

    def se(self, key: str, d: int = 3) -> str:
        r = self.r(key)
        return "n/a" if r is None or not finite(r.get("se")) else f"{r['se']:.{d}f}"

    def c(self, key: str, d: int = 3) -> str:
        r = self.r(key)
        return "n/a" if r is None else f"{r['value']:.{d}f} ({r['se']:.{d}f})"

    def ci(self, key: str, d: int = 2) -> str:
        r = self.r(key)
        return "n/a" if r is None else f"[{r['ci95_lo']:.{d}f}, {r['ci95_hi']:.{d}f}]"

    def val(self, key: str) -> float:
        r = self.r(key)
        return float("nan") if r is None else r["value"]

    def cnt(self, key: str) -> dict[str, Any]:
        return self.P["counts"].get(key, {"person_years": 0, "persons": 0})


def _tab(header: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" if i == 0 else "---:" for i in range(len(header))) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


SENS_LABEL = {"core": "core K-12 codes", "half": "teacher in at least half of the rounds", "ever": "ever a teacher", "any": "any 25-34 occupation",
              "c70": "1970 codes", "round": "round-level K-12 occupation", "deflage": "deflated wage, no year dummies", "deflnone": "deflated wage, female dummy only"}


def write_md(P: dict[str, Any], path: Path) -> None:
    R = Rep(P)
    M = P["meta"]
    pub, rep, d3 = P["published"], P["replay"], P["d3"]
    vers = [R.val("ba_pooled|naive"), R.val("ba_pooled|eiv"), R.val("iv_ba|mean4")]
    spread = max(vers) - min(vers)                      # convention, EIV and split-half IV versions of the BA-conditional slope
    L: list[str] = []
    ap = L.append
    g3 = list(GROUPS)
    glab = {"pooled": "Pooled", "male": "Men", "female": "Women"}
    dmm, dme = P["delta_method"]["main"]["naive"], P["delta_method"]["main"]["eiv"]
    rb = ratio_of_gamma(GAMMA_BASELINE)

    ap("# NLSY79: wage-ability slopes holding the BA status fixed, teachers' slope ratio, and the BA-share logit shift")
    ap("")
    ap(f"Generated {M['generated']} by `nlsy_ba_slopes.py` (task D3b of the 2026-10b round, DESIGN_v2.md section M3; second version, after the review of the first). NLSY79, ages 25-34, "
       f"the wage sample, wage definition, 1979-weight convention, latent-AFQT standardization (`afqt_c`, R_comp_pop = {R.v('R_comp_pop', 3)} "
       f"({R.se('R_comp_pop', 4)})) and controls (age and survey-year dummies; a female dummy in pooled samples) of `nlsy_ability.py` block C. "
       f"Household-cluster bootstrap, B = {M['B']}, seed {M['seed']} (the draws of `nlsy_teachers.py` and D3); every draw re-standardizes AFQT, "
       f"re-estimates R and re-runs every regression and shift. SE = bootstrap SD in parentheses; 95% percentile CI in brackets. "
       f"All slopes are log hourly-wage points per SD of latent AFQT. NLSY97 and the sealed hold-out are not used (see the sealed-set note under Definitions and caveats).")
    ap("")
    ap("Two disattenuation versions are reported side by side, with a split-half IV check (section 2). **Convention** (the `value` of every slope in the targets file): observed slope / sqrt(R), "
       "the project's latent slope; it gives the published 0.2127, so the BA-conditional slope and the unconditional target are on one definition. **EIV** (carried as the alternative): "
       "the classical errors-in-variables correction of the regression itself (error variance 1 - R in afqt_c units); it equals the convention for an unconditional regression on the full block-A sample, "
       f"but in subsamples and with a BA control the conditional reliability of the score is lower than R, so it is {100 * (R.val('ba_pooled|eiv') / R.val('ba_pooled|naive') - 1):.0f}% higher "
       "for the BA-controlled slope. Both estimate the slope on the true latent score, which is what the model's `wage_slope` is; the choice between them is a PREREG matter (V12, the unconditional slope, is on the convention).")
    ap("")

    # ------------------------------------------------------------------ summary
    ap("## Summary")
    ap("")
    u, b_, i_ = "u_pooled", "ba_pooled", "int_pooled"
    ap(f"- **Reproduction.** The pooled unconditional latent slope is {R.v(u + '|naive', 4)} against the published {pub['slopes']['pooled']['value']:.4f} "
       f"(SE {R.se(u + '|naive', 4)} here, {pub['slopes']['pooled']['se']:.4f} published; replaying the published bootstrap stream, {rep['pooled']['se']:.4f}). "
       f"Men {R.v('u_male|naive', 4)}, women {R.v('u_female|naive', 4)} reproduce too (section 1).")
    ap(f"- **Slope holding BA fixed (the M3 target `wage_slope_ba`; convention).** Adding a BA indicator (highest grade >= 16) lowers the pooled latent slope from {R.v(u + '|naive', 4)} to "
       f"**{R.c(b_ + '|naive', 4)}** on the project's convention (observed slope / sqrt(R), the pipeline of the 0.2127 target): a fall of "
       f"{100 * (1 - R.val('ratio_ba_u_pooled|naive')):.0f}% (SE {100 * float(R.r('ratio_ba_u_pooled|naive')['se']):.0f} points). "
       f"Men {R.c('ba_male|naive', 4)}, women {R.c('ba_female|naive', 4)}. The EIV version, which lets the reliability of the score fall to {R.v(b_ + '|lam', 3)} with the BA control, is "
       f"{R.c(b_ + '|eiv', 4)} (unconditional EIV {R.c(u + '|eiv', 4)}); a split-half IV check that assumes no error variance gives {R.c('iv_ba|mean4', 4)}; "
       f"the three versions span {spread:.3f} ({'less' if spread < float(R.r(b_ + '|naive')['se']) else 'more'} than one SE). **Pairing:** {R.v(b_ + '|naive', 4)} goes with the unconditional 0.2127 (both convention); "
       f"if EIV is preferred, use {R.v(b_ + '|eiv', 4)} together with {R.v(u + '|eiv', 4)} (and amend PREREG V12 before any M3 result). The BA coefficient at a given latent ability is "
       f"{R.c(b_ + '|ba_coef_eiv')} (EIV; {R.c(b_ + '|ba_coef_naive')} uncorrected). First-stage reading s_z = (1 - eta) x slope = "
       f"{(1 - ETA) * R.val(b_ + '|naive'):.3f} (convention) or {(1 - ETA) * R.val(b_ + '|eiv'):.3f} (EIV) against {pub['implied_s_z']:.4f} now; the joint estimation with the college margin re-estimates it.")
    zd = R.val(i_ + '|diff_naive') / float(R.r(i_ + '|diff_naive')['se'])
    ap(f"- **BA interacted.** Within non-BA {R.c(i_ + '|nonba_naive')} / {R.c(i_ + '|nonba_eiv')}; within BA {R.c(i_ + '|ba_naive')} / {R.c(i_ + '|ba_eiv')} (convention / EIV); "
       f"BA minus non-BA {R.c(i_ + '|diff_naive')} / {R.c(i_ + '|diff_eiv')} (z = {zd:.1f} for the convention): the slope is about the same with and without a BA, so the fall from the unconditional slope is the BA level, not a flatter BA gradient.")
    tm, nb = R.cnt("TC_main"), R.cnt("NTBA_main")
    sens = [t for t in SENS_LABEL if R.r(f"gamma_{t}|naive") is not None]
    lo_t = min(sens, key=lambda t: R.val(f"gamma_{t}|naive"))
    hi_t = max(sens, key=lambda t: R.val(f"gamma_{t}|naive"))
    nob_sd = abs(R.val("r_main|naive") - rb) / float(R.r("r_main|naive")["se"])
    ap(f"- **Teachers against non-teacher BA workers (the BA-controlled slope).** With the BA indicator in the regression (the slope of `wage_slope_ba`), K-12 teachers ({tm['persons']} persons, {tm['person_years']} wage person-years) "
       f"have a latent slope of {R.c('TC_main|naive')} (convention; EIV {R.c('TC_main|eiv')}) and non-teacher BA workers ({nb['persons']:,} persons, {nb['person_years']:,} person-years) {R.c('NTBA_main|naive')} (EIV {R.c('NTBA_main|eiv')}). "
       f"Ratio **{R.c('r_main|naive', 2)}** (EIV {R.c('r_main|eiv', 2)}), implied **gamma = {R.c('gamma_main|naive', 2)}** (EIV {R.c('gamma_main|eiv', 2)}) at eta = {ETA:.3f} "
       f"(bootstrap CI {R.ci('gamma_main|naive')}; delta-method SE {dmm['se_gamma']:.2f}). The baseline gamma = {GAMMA_BASELINE} implies a ratio of {rb:.2f}, {nob_sd:.1f} SE from the estimate: "
       f"the ratio does not reject the baseline and is too imprecise to be a target. The first version of this file compared the teachers' slope WITHOUT the BA control ({R.c('T_main|naive')}; ratio {R.c('r_nobac|naive', 2)}, gamma {R.c('gamma_nobac|naive', 2)}): "
       f"the {100 * (1 - R.cnt('TBA_main')['persons'] / max(1, tm['persons'])):.0f}% of teachers without a BA score lower and, on the point estimate, are paid less, so that slope contains the credential composition, which the closed form has no margin for; "
       f"the BA control removes it. The answer also moves with the teacher definition: gamma runs from {R.v('gamma_' + lo_t + '|naive', 2)} ({SENS_LABEL[lo_t]}) to {R.v('gamma_' + hi_t + '|naive', 2)} ({SENS_LABEL[hi_t]}); "
       f"BA teachers only {R.c('gamma_tba|naive', 2)}; against all non-teachers (BA control in both) {R.c('gamma_ntc|naive', 2)}.")
    ap(f"- **BA-share shift, like for like.** The NLSY79 share of ALL persons of D3's primary sample ({R.c('ba|share', 3)}; it contains the K-12 teachers, {100 * R.val('ba|teacher_ba_share'):.0f}% of whom have a BA) moves to the ACS all-persons share "
       f"{d3['ba_share_all']:.4f} with a common logit shift of the D3 tercile shares of **{R.c('ba|shift_cell_all')}**; implied P(BA | AFQT tercile) = {R.c('ba|imp_cell_all_terc1', 3)}, {R.c('ba|imp_cell_all_terc2', 3)}, {R.c('ba|imp_cell_all_terc3', 3)} "
       f"(NLSY79: {R.v('ba|terc1')}, {R.v('ba|terc2')}, {R.v('ba|terc3')}); D3's individual-logit shift to the same share, {R.c('ba|shift_indiv_all')}, is reproduced. The same exercise for the NON-TEACHERS (NLSY79 share {R.c('banT|share', 3)}, "
       f"ACS non-teacher share {d3['ba_share_nonteacher']:.4f}) gives **{R.c('banT|shift_cell_nt')}** and {R.c('banT|imp_cell_nt_terc1', 3)}, {R.c('banT|imp_cell_nt_terc2', 3)}, {R.c('banT|imp_cell_nt_terc3', 3)}. "
       f"Shifting the all-persons shares to the non-teacher share ({R.c('ba|shift_cell_mixed')}; shares {R.v('ba|imp_cell_mixed_terc1')}, {R.v('ba|imp_cell_mixed_terc2')}, {R.v('ba|imp_cell_mixed_terc3')}), "
       f"the version the task text asked for and the first version of this file, mixes two populations and understates the shift by {R.val('ba|shift_cell_all') - R.val('ba|shift_cell_mixed'):.2f} log points: it is kept only as a flagged alternative.")
    ap("")

    # ------------------------------------------------------------------ 1
    ap("## 1. Reproduction of the published slope")
    ap("")
    rows = [["Published (`nlsy_ability.json`, B = 200, seed 20260928)"] + [f"{pub['slopes'][g]['value']:.5f} ({pub['slopes'][g]['se']:.4f})" for g in g3],
            ["This script, point estimate (convention)"] + [R.v(f"u_{g}|naive", 5) for g in g3],
            ["  absolute difference to the published value"] + [f"{abs(R.val(f'u_{g}|naive') - pub['slopes'][g]['value']):.1e}" for g in g3],
            ["SE, replay of the published bootstrap stream (B = 200)"] + [f"{rep[g]['se']:.5f}" for g in g3],
            [f"SE, draws of D3 and `nlsy_teachers.py` (B = {M['B']}, seed {M['seed']})"] + [R.se(f"u_{g}|naive", 5) for g in g3],
            ["EIV version (point and SE)"] + [R.c(f"u_{g}|eiv", 4) for g in g3],
            ["Conditional reliability of afqt_c in the wage sample (R_comp_pop = " + R.v("R_comp_pop", 3) + ")"] + [R.v(f"u_{g}|lam", 3) for g in g3],
            ["Person-years / persons"] + [f"{R.cnt('u_' + g)['person_years']:,} / {R.cnt('u_' + g)['persons']:,}" for g in g3]]
    ap(_tab([""] + [glab[g] for g in g3], rows))
    ap("")
    ap(f"The published stream is replayed exactly (multinomial household draws with the interleaved NLSY97 draw, which is generated for the stream and not read; "
       f"`nlsy_ability.slope_set` itself computes the slopes): maximum absolute difference to the published SEs "
       f"{max(abs(rep[g]['se'] - pub['slopes'][g]['se']) for g in g3):.1e}. D3 is reproduced the same way: BA share by tercile "
       f"{d3['check']['max_abs_diff_terc_shares']:.1e}, logit coefficients {d3['check']['max_abs_diff_logit']:.1e}, shift to the all-persons ACS share {d3['check']['diff_shift_all']:.1e} "
       f"(maximum absolute differences to `education_occupation.json`).")
    ap("")

    # ------------------------------------------------------------------ 2
    ap("## 2. Slopes holding the BA status fixed")
    ap("")
    ap(f"BA = highest grade completed >= 16 at some round at 25-34 (a proxy: the project's NLSY79 extract has no degree variable); {M['ba_share_wage']:.3f} of the weighted wage-sample person-years. "
       f"Pooled regressions include a female dummy; by-sex regressions do not. 'Interacted' = one regression with BA, AFQT x BA and the common controls; "
       f"'separate samples' = the regression run within each BA group with its own controls.")
    ap("")
    spec_rows = [("No BA control (the existing slope)", lambda g, sc: R.c(f"u_{g}|{sc}")),
                 ("**BA control**", lambda g, sc: R.c(f"ba_{g}|{sc}")),
                 ("BA control, BA = highest grade >= 16 at the round itself", lambda g, sc: R.c(f"bar_{g}|{sc}")),
                 ("Interacted: slope within non-BA", lambda g, sc: R.c(f"int_{g}|nonba_{sc}")),
                 ("Interacted: slope within BA", lambda g, sc: R.c(f"int_{g}|ba_{sc}")),
                 ("Interacted: BA minus non-BA", lambda g, sc: R.c(f"int_{g}|diff_{sc}")),
                 ("Separate samples: non-BA", lambda g, sc: R.c(f"sepnba_{g}|{sc}")),
                 ("Separate samples: BA", lambda g, sc: R.c(f"sepba_{g}|{sc}")),
                 ("BA control / no control", lambda g, sc: R.c(f"ratio_ba_u_{g}|{sc}", 3))]
    for sc, title in (("naive", "Convention (observed slope / sqrt(R)): the target values"), ("eiv", "EIV (classical errors-in-variables correction): the alternative")):
        ap(f"**{title}**")
        ap("")
        ap(_tab([""] + [glab[g] for g in g3], [[lab] + [fn(g, sc) for g in g3] for lab, fn in spec_rows]))
        ap("")
    ap("**BA coefficient and conditional reliabilities**")
    ap("")
    rows = [["BA coefficient in the BA-control regression, uncorrected"] + [R.c(f"ba_{g}|ba_coef_naive") for g in g3],
            ["BA coefficient in the BA-control regression, EIV (premium at given latent ability)"] + [R.c(f"ba_{g}|ba_coef_eiv") for g in g3],
            ["Conditional reliability of afqt_c: BA-control regression"] + [R.v(f"ba_{g}|lam") for g in g3],
            ["Conditional reliability: BA holders only"] + [R.v(f"sepba_{g}|lam") for g in g3],
            ["Conditional reliability: non-BA only"] + [R.v(f"sepnba_{g}|lam") for g in g3],
            ["Person-years / persons, BA holders"] + [f"{R.cnt('sepba_' + g)['person_years']:,} / {R.cnt('sepba_' + g)['persons']:,}" for g in g3],
            ["Person-years / persons, non-BA"] + [f"{R.cnt('sepnba_' + g)['person_years']:,} / {R.cnt('sepnba_' + g)['persons']:,}" for g in g3]]
    ap(_tab([""] + [glab[g] for g in g3], rows))
    ap("")
    ap(f"The reliability of afqt_c is {R.v('R_comp_pop', 3)} in the block-A sample; in the wage sample it is {R.v('u_pooled|lam', 3)}, with a BA control {R.v('ba_pooled|lam', 3)}, and "
       f"within a BA group {R.v('sepba_pooled|lam', 3)} (weighted SD of afqt_c: {M['sd_afqt_wage']:.2f} in the wage sample, {M['sd_afqt_ba']:.2f} among BA holders, {M['sd_afqt_nba']:.2f} among non-BA). "
       f"Dividing by sqrt(R) therefore under-corrects by "
       f"{100 * (R.val('ba_pooled|eiv') / R.val('ba_pooled|naive') - 1):.0f}% with the BA control and by {100 * (R.val('sepba_pooled|eiv') / R.val('sepba_pooled|naive') - 1):.0f}% within BA, "
       f"if the error variance of the score is 1 - R in every group. "
       f"The BA coefficient falls by {R.val('ba_pooled|ba_coef_naive') - R.val('ba_pooled|ba_coef_eiv'):.3f} under EIV because the noisy AFQT leaves part of ability to the BA dummy.")
    ap("")

    ap("**Check of the disattenuation: split-half instrumental variables**")
    ap("")
    rows = [["Convention (observed / sqrt(R))", R.c("u_pooled|naive", 4), R.c("ba_pooled|naive", 4)],
            ["EIV (one-factor reliability)", R.c("u_pooled|eiv", 4), R.c("ba_pooled|eiv", 4)]]
    for tg, (_, _, lab) in IV_SPLITS.items():
        rows.append([f"IV, {lab}: A instrumented by B / B by A", f"{R.c(f'iv_{tg}_none|AB', 4)} / {R.c(f'iv_{tg}_none|BA', 4)}", f"{R.c(f'iv_{tg}_ba|AB', 4)} / {R.c(f'iv_{tg}_ba|BA', 4)}"])
    rows.append(["**IV, mean of the two mixed splits (4 estimates)**", R.c("iv_none|mean4", 4), R.c("iv_ba|mean4", 4)])
    ap(_tab(["Pooled, per SD of latent AFQT", "No BA control", "BA control"], rows))
    ap("")
    ap("Each half of the four ASVAB tests (AR, WK, PC, MK) is a composite; the regressor is one half instrumented by the other half and rescaled by the one-factor loading of its composite "
       "(so the factor is the project's, but no error variance is assumed and no group-specific reliability is needed). The two mixed splits give nearly the same slope in both directions; "
       "the math-verbal split does not (math and verbal composites carry wage-relevant skills beyond the common factor), so it is shown but not averaged. "
       f"The mixed-split IV puts the BA-controlled slope at {R.v('iv_ba|mean4', 3)} ({R.se('iv_ba|mean4', 3)}): {100 * (R.val('iv_ba|mean4') / R.val('ba_pooled|naive') - 1):+.0f}% against the convention and "
       f"{100 * (R.val('iv_ba|mean4') / R.val('ba_pooled|eiv') - 1):+.0f}% against EIV; the three versions span {spread:.3f}, {'less' if spread < float(R.r('ba_pooled|naive')['se']) else 'more'} than one SE ({R.se('ba_pooled|naive', 3)}). "
       f"The unconditional IV slope is {R.v('iv_none|mean4', 3)}, {100 * (R.val('iv_none|mean4') / R.val('u_pooled|naive') - 1):+.0f}% against the published {R.v('u_pooled|naive', 3)}: "
       "the IV agrees with EIV that conditioning on BA raises the needed correction (IV / convention rises from "
       f"{R.val('iv_none|mean4') / R.val('u_pooled|naive'):.3f} to {R.val('iv_ba|mean4') / R.val('ba_pooled|naive'):.3f}, EIV / convention from {R.val('u_pooled|eiv') / R.val('u_pooled|naive'):.3f} to {R.val('ba_pooled|eiv') / R.val('ba_pooled|naive'):.3f}) "
       "but puts the level lower than the one-factor reliability does.")
    ap("")

    # ------------------------------------------------------------------ 3
    ap("## 3. Teachers and non-teacher BA workers: slopes, ratio and implied gamma")
    ap("")
    ap(f"Teachers: the `nlsy_teachers.py` person-level classification (Census 1980 codes, all K-12 codes, a K-12 teacher in the majority of the employed 25-34 rounds, "
       f"spec `{MAIN_SPEC}`; {M['n_teacher_persons']} persons, {M['n_teacher_women']} women and {M['n_teacher_men']} men), all of whose wage-sample person-years enter. "
       f"Non-teacher BA workers: the non-teachers of that spec with a BA. The teachers' slope is the SAME slope as in section 2, i.e. from a regression with a BA indicator (the non-teacher BA workers all have a BA, so there the control is constant): "
       f"{100 * (1 - R.cnt('TBA_main')['persons'] / max(1, tm['persons'])):.0f}% of the teachers have no BA; they score lower (D3: the teachers' latent AFQT gap is {d3['teacher_gap_vs_ba_nonteachers']:+.3f} SD against non-teacher BA workers and {d3['teacher_gap_within_ba']:+.3f} SD for BA teachers only) and, on the point estimate, are paid less "
       f"(BA coefficient among teachers {R.c('TC_main|ba_coef_naive', 3)}, uncorrected), so a slope without the control would put the credential composition into the teachers' slope, "
       f"for which the closed form has no margin. The gamma formula is the model's: ratio = gamma (1 - eta) / (1 - gamma eta), so gamma = ratio / (1 - eta + eta ratio), eta = {ETA:.3f} (0.0807 changes gamma by less than 0.001).")
    ap("")
    grp = [("K-12 teachers, BA control (main)", "TC_main"), ("  women", "TC_women"), ("  men", "TC_men"),
           ("K-12 teachers, no BA control (first version)", "T_main"), ("  women", "T_women"), ("  men", "T_men"), ("K-12 teachers with a BA", "TBA_main"),
           ("All non-teachers, no BA control", "NT_main"), ("All non-teachers, BA control", "NTC_main"),
           ("Non-teacher BA workers", "NTBA_main"), ("  women", "NTBA_women"), ("  men", "NTBA_men"), ("Non-teachers without a BA", "NTnBA_main")]
    rows = [[lab, f"{R.cnt(k)['persons']:,}", f"{R.cnt(k)['person_years']:,}", R.c(f"{k}|obs"), R.c(f"{k}|naive"), R.c(f"{k}|eiv"), R.v(f"{k}|lam")] for lab, k in grp]
    ap(_tab(["Group (wage sample)", "Persons", "Person-years", "Observed slope per afqt_c SD", "Latent, convention", "Latent, EIV", "Cond. reliability"], rows))
    ap("")
    ap(f"Pooled groups include a female dummy, women and men do not; all regressions have age and year dummies. The teacher slopes rest on {tm['persons']} persons. "
       f"Among teachers the BA indicator's coefficient is {R.c('TC_main|ba_coef_eiv', 3)} (EIV; {R.c('TC_main|ba_coef_naive', 3)} uncorrected), identified from the {R.cnt('T_main')['persons'] - R.cnt('TBA_main')['persons']} teachers without a BA; "
       f"for the {R.cnt('TC_men')['persons']} male teachers, of whom very few lack a BA, the coefficient and so the men's row are unreliable.")
    ap("")
    ap("**Ratio of slopes and implied gamma** (convention = the values; EIV = the alternative)")
    ap("")
    order = ["main", "nobac", "tba", "ntc", "nt", "women", "men", "core", "half", "ever", "any", "c70", "round", "deflage", "deflnone"]
    rows = []
    for tg in order:
        num, den, desc = P["ratios"][tg]
        cn = R.cnt(num)
        rows.append([("**" + desc + "**") if tg == "main" else desc, f"{cn['persons']:,} / {cn['person_years']:,}", R.c(f"r_{tg}|naive", 2), R.c(f"gamma_{tg}|naive", 2),
                     R.ci(f"gamma_{tg}|naive"), R.c(f"r_{tg}|eiv", 2), R.c(f"gamma_{tg}|eiv", 2)])
    ap(_tab(["Comparison (numerator = teachers)", "Teacher persons / person-years", "Ratio, convention", "gamma, convention", "gamma 95% CI", "Ratio, EIV", "gamma, EIV"], rows))
    ap("")
    ap(f"Main specification, delta method (bootstrap covariance of the two slopes; correlation {dmm['corr_slopes']:.2f} convention, {dme['corr_slopes']:.2f} EIV): convention ratio {dmm['ratio']:.3f} (SE {dmm['se_ratio']:.3f}), "
       f"gamma {dmm['gamma']:.3f} (SE {dmm['se_gamma']:.3f}; d gamma / d ratio = {dmm['dgamma_dratio']:.3f}), against the bootstrap SDs {R.se('r_main|naive')} and {R.se('gamma_main|naive')}; "
       f"EIV ratio {dme['ratio']:.3f} (SE {dme['se_ratio']:.3f}), gamma {dme['gamma']:.3f} (SE {dme['se_gamma']:.3f}), against {R.se('r_main|eiv')} and {R.se('gamma_main|eiv')}. "
       f"Model reference: gamma = {GAMMA_BASELINE} (baseline fit) implies ratio {rb:.3f}; gamma = 1 implies 1.")
    ap("")
    ap("Why the rows differ: the teacher slope is estimated from 100-400 persons, so the ratio has an SE of 0.2-0.4 (0.8 for men) and the rows are not independent evidence; they change who counts "
       "as a teacher (core K-12 codes drop the pre-K and 'other teacher' codes, whose members have lower AFQT and pay; 'ever' and the round-level row add persons who teach in a minority of rounds), "
       "the comparison group (all non-teachers carry the BA premium in their slope unless the BA control is in both regressions) or the credential control (the first-version row, without it, "
       "puts the teachers without a BA into the teachers' slope). The main interval contains both the baseline gamma and gamma = 1.")
    ap("")

    # ------------------------------------------------------------------ 4
    ap("## 4. The BA share: logit level shift and implied tercile shares")
    ap("")
    ap(f"Primary NLSY79 sample of D3 (persons observed at 25-34, highest grade >= 9, N = {M['n_ba_sample']:,}); BA share {R.c('ba|share', 4)}; terciles of the observed `afqt_c` cut over the weighted sample "
       f"(cut points {R.v('ba|cut1')}, {R.v('ba|cut2')} SDs), as in D3, for both populations below. **A shift is only meaningful against the ACS share of the same population.** The NLSY79 share above is that of ALL persons, "
       f"{100 * R.val('ba|teacher_wshare'):.1f}% of whom (weighted) are K-12 teachers (main spec) with a BA share of {100 * R.val('ba|teacher_ba_share'):.0f}%; the non-teachers' NLSY79 share is {R.c('banT|share', 4)}. Two like-for-like targets: "
       f"(1) all persons to the ACS all-persons share {d3['ba_share_all']:.4f} (`ba_share_all`, D3's own target; the model's `ba_share` and `ba_share_by_z` are over all adults), "
       f"(2) the non-teachers (everyone not a teacher, persons never employed included, as the ACS non-teacher universe includes home production) to the ACS non-teacher share {d3['ba_share_nonteacher']:.4f} "
       f"({d3['ba_share_nonteacher_se']:.4f}, T2 convention, `education_targets.toml` `ba_share_nonteacher`, the M3 target). The third version, all persons to the non-teacher share, is NOT like for like and is shown only because the task text and the first version of this file used it. "
       f"**Cell shift**: one delta with sum_k pi_k Lambda(logit s_k + delta) = target, s_k the D3 tercile shares of the population and pi_k its tercile weights. **Individual logit**: P(BA | afqt) = Lambda(a + delta + b afqt) "
       f"with (a, b) = ({R.v('ba|logit_a', 3)}, {R.v('ba|logit_b', 3)}) for all persons ({R.v('banT|logit_a', 3)}, {R.v('banT|logit_b', 3)} for non-teachers) fitted by weighted logit, delta set so the weighted mean is the target (D3's 'level shift of the same gradient'); "
       f"the tercile shares are its weighted averages within the NLSY79 terciles.")
    ap("")
    rows = [["All persons, as observed (D3, not recomputed for use)", R.c("ba|share", 4), R.c("ba|terc1", 4), R.c("ba|terc2", 4), R.c("ba|terc3", 4), R.c("ba|terc_diff", 3), "-"],
            [f"**All persons, cell shift to the ACS all-persons share {d3['ba_share_all']:.4f} (like for like; primary)**", f"{d3['ba_share_all']:.4f}", R.c("ba|imp_cell_all_terc1", 4), R.c("ba|imp_cell_all_terc2", 4),
             R.c("ba|imp_cell_all_terc3", 4), R.c("ba|imp_cell_all_diff", 3), R.c("ba|shift_cell_all")],
            [f"All persons, individual logit shifted to {d3['ba_share_all']:.4f} (D3's method; D3: {d3['stats']['logit|all_hgc9|pooled|shift']['value']:.4f})", R.v("ba|imp_indiv_all_share", 4), R.c("ba|imp_indiv_all_terc1", 4),
             R.c("ba|imp_indiv_all_terc2", 4), R.c("ba|imp_indiv_all_terc3", 4), R.c("ba|imp_indiv_all_diff", 3), R.c("ba|shift_indiv_all")],
            ["All persons, individual logit fit, no shift", R.v("ba|share", 4), R.c("ba|fit_indiv_terc1", 4), R.c("ba|fit_indiv_terc2", 4), R.c("ba|fit_indiv_terc3", 4),
             f"{R.val('ba|fit_indiv_terc3') - R.val('ba|fit_indiv_terc1'):.3f}", "0"],
            ["Non-teachers, as observed (same terciles)", R.c("banT|share", 4), R.c("banT|terc1", 4), R.c("banT|terc2", 4), R.c("banT|terc3", 4), R.c("banT|terc_diff", 3), "-"],
            [f"**Non-teachers, cell shift to the ACS non-teacher share {d3['ba_share_nonteacher']:.4f} (like for like)**", f"{d3['ba_share_nonteacher']:.4f}", R.c("banT|imp_cell_nt_terc1", 4), R.c("banT|imp_cell_nt_terc2", 4),
             R.c("banT|imp_cell_nt_terc3", 4), R.c("banT|imp_cell_nt_diff", 3), R.c("banT|shift_cell_nt")],
            [f"Non-teachers, individual logit shifted to {d3['ba_share_nonteacher']:.4f}", R.v("banT|imp_indiv_nt_share", 4), R.c("banT|imp_indiv_nt_terc1", 4), R.c("banT|imp_indiv_nt_terc2", 4),
             R.c("banT|imp_indiv_nt_terc3", 4), R.c("banT|imp_indiv_nt_diff", 3), R.c("banT|shift_indiv_nt")],
            [f"*Mixed, NOT like for like: all persons shifted to the non-teacher share {d3['ba_share_nonteacher']:.4f} (task text; first version)*", f"{d3['ba_share_nonteacher']:.4f}", R.c("ba|imp_cell_mixed_terc1", 4),
             R.c("ba|imp_cell_mixed_terc2", 4), R.c("ba|imp_cell_mixed_terc3", 4), R.c("ba|imp_cell_mixed_diff", 3), R.c("ba|shift_cell_mixed")],
            [f"*Mixed, individual logit (D3's method) to {d3['ba_share_nonteacher']:.4f}*", R.v("ba|imp_indiv_mixed_share", 4), R.c("ba|imp_indiv_mixed_terc1", 4), R.c("ba|imp_indiv_mixed_terc2", 4),
             R.c("ba|imp_indiv_mixed_terc3", 4), R.c("ba|imp_indiv_mixed_diff", 3), R.c("ba|shift_indiv_mixed")]]
    ap(_tab(["", "BA share", "Tercile 1 (low AFQT)", "Tercile 2", "Tercile 3 (high)", "T3 - T1", "Logit shift delta"], rows))
    ap("")
    dca = R.val("ba|dshift_cell_all") * d3["ba_share_all_se"]
    dcn = R.val("banT|dshift_cell_nt") * d3["ba_share_nonteacher_se"]
    ap(f"The two like-for-like cell shifts agree ({R.val('ba|shift_cell_all'):.3f} for all persons, {R.val('banT|shift_cell_nt'):.3f} for non-teachers), as they should: both bring a cohort whose BA share is lower than the ACS one to the ACS level, "
       f"whereas the mixed version asks the NLSY79 all-persons distribution, which includes the high-BA teachers, to reach a level that excludes them. The cell shift keeps the empirical tercile pattern (its middle tercile is the D3 one shifted); "
       f"the individual logit imposes a log-odds line in afqt_c, which fits the middle tercile {R.val('ba|fit_indiv_terc2') - R.val('ba|terc2'):+.3f} off and so shifts it differently. The SEs are NLSY79 sampling errors; the ACS targets' own SEs "
       f"({d3['ba_share_all_se']:.4f} all persons, {d3['ba_share_nonteacher_se']:.4f} non-teachers) add {dca:.4f} and {dcn:.4f} to delta, negligible. The shares are by tercile of the OBSERVED score (reliability {R.v('R_comp_pop', 3)}): they understate the latent gradient, "
       f"so the model counterpart should add classical noise with variance (1 - R) / R = {(1 - R.val('R_comp_pop')) / R.val('R_comp_pop'):.3f} times the variance of the latent score before cutting terciles. "
       f"The NLSY79 teacher share of persons ({100 * R.val('ba|teacher_wshare'):.1f}%) is lower than the ACS one ({100 * d3['acs_teacher_share']:.1f}% of the Table 2 universe), which is why the non-teacher and all-persons levels are closer in NLSY79 "
       f"({R.val('ba|share') - R.val('banT|share'):.3f}) than in the ACS ({d3['ba_share_all'] - d3['ba_share_nonteacher']:.3f}).")
    ap("")

    # ------------------------------------------------------------------ definitions and caveats
    ap("## Definitions and caveats")
    ap("")
    ap("- **Wage sample, weights, controls**: as `nlsy_ability.py` block C (>= 48 weeks, >= 30 hours/week, wage income >= $1,000 in 2010 dollars, highest grade >= 9, not armed forces; "
       "log of annual income over annual hours, no trimming; person's 1979 SAMPWEIGHT on each person-year; age and year dummies). The CPI-deflated rows of the ratio table replace the year dummies by deflation with the prior-year CPI-U.")
    ap("- **Convention and EIV**: the convention is OLS / sqrt(R), the project's latent slope (0.2127 for the pooled unconditional regression). EIV: afqt_c = f + u, var(u) = 1 - R, var(f) = R, u independent of the other regressors; the coefficients solve (M - S) b = m with S = (1 - R) E[c c'] "
       "(c = 1 for the AFQT column, BA for AFQT x BA); the slope per SD of f is b_AFQT sqrt(R). The EIV assumptions are strong: the error variance is taken equal across groups, independent of BA, "
       "and AFQT was measured at 15-23, when schooling already raises scores; a one-factor ULS model supplies R. The targets file carries the convention as the value and EIV as `alt_eiv`; the two are not to be mixed within a pair "
       "(BA-conditional slope and unconditional slope).")
    dropped = [(k[:-4], int(P["statistics"][k.replace("|eiv", "|naive")]["n_draws"] - P["statistics"][k]["n_draws"])) for k in P["statistics"]
               if k.endswith("|eiv") and k.replace("|eiv", "|naive") in P["statistics"]
               and P["statistics"][k.replace("|eiv", "|naive")]["n_draws"] > P["statistics"][k]["n_draws"]]
    if dropped:
        Bn = M["B"]
        drop_txt = ", ".join("`%s` (%d of %d draws)" % (k, n, Bn) for k, n in dropped)
        ap(f"- **Unstable EIV draws**: an EIV statistic is set to missing in a draw whose conditional reliability of the score is below {LAM_MIN} (the correction would more than double the slope). "
           f"This affects only {drop_txt}; their SEs use the remaining draws.")
    ap("- **Split-half IV**: one half of the tests (a composite of two standardized IRT z-scores, standardized in block A) instrumented by the other half, controls partialled out, times the loading of the composite on the standardized one-factor "
       "latent variable (mean ULS loading of its two tests over the composite's SD); re-estimated in every draw. It assumes independent errors across halves, not an error variance.")
    ap("- **BA** is highest grade completed >= 16 at some round at 25-34, a proxy for a bachelor's degree: some respondents with 16 years report no bachelor's degree and others with 17-18 years do, and graduate degrees are not separated. "
       "NLSY79 itself records degrees (questions 1979-2020; the created variable HIGHEST_DEGREE_EVER was added in the 2022, round 30, release: NLSY79 Codebook Supplement, Appendix 8), "
       "but the project's extract has no degree variable; adding it needs a new NLS Investigator extract. The model's d = 1 is a BA, and its teachers are mostly BA holders (87% in the ACS).")
    ap(f"- **Teacher sample**: {tm['persons']} persons, {tm['person_years']} person-years in the wage sample ({R.cnt('TC_women')['persons']} women, {R.cnt('TC_men')['persons']} men), so teacher slopes carry SEs of "
       f"{float(R.r('TC_main|naive')['se']):.2f} pooled and {float(R.r('TC_men|naive')['se']):.2f} for men; the classification itself is held fixed in the bootstrap. Teachers' wage-sample person-years "
       f"include rounds in other occupations ({100 * M['teacher_rows_other_occ_share']:.0f}%). The round-level row uses person-years in a K-12 occupation instead.")
    ap("- **Cohort**: NLSY79 persons born 1957-64 observed 1982-98, against the ACS 2009-13 cohorts of the model's BA share. The NLSY97 BA-conditional slope (the closer cohort; DESIGN_v2.md lists it among the unsealed M3 inputs) is outside this NLSY79 task and was not computed.")
    ap("- **Not included in the SEs**: the ACS targets' SEs (shown above for the shifts), eta (0.080 against 0.0807), the choice of BA and teacher definitions (the sensitivity rows) and the sampling error of the assumed EIV structure.")
    ap(f"- **By sex**: women's slopes are larger than men's (BA-conditional {R.v('ba_female|naive')} against {R.v('ba_male|naive')}, convention; unconditional {R.v('u_female|naive')} against {R.v('u_male|naive')}); the model has one slope for both sexes.")
    ap("- **Sealed set**: no sealed value was read or used and no NLSY97 data were read. One protocol deviation, from the first run of this task (2026-10-05): the sealed file `estimates/settings/sealed/nlsy97_holdout.json` was loaded once and only its "
       "top-level key names were printed (no value, no statistic). PREREG.md has the sealed set 'opened only after selection', so this is an access before selection; it is to be recorded in the round's LOG.md. "
       "Nothing in this file or in the targets file depends on it.")
    ap("")
    ap("## Not done")
    ap("")
    ap(f"- The project's extract has no degree variable, so BA is a schooling-years proxy; no separate treatment of graduate degrees or of the {100 * (1 - R.cnt('TBA_main')['persons'] / max(1, tm['persons'])):.0f}% of teachers without a BA beyond the BA control and the BA-teacher rows.")
    ap("- The NLSY97 counterparts (the cohort closer to the ACS 2009-13 population) were not computed: the BA-conditional slope is a separate unsealed input, and the teachers' wage-AFQT slope of the NLSY97 belongs to the sealed hold-out (D4), which is not touched here.")
    ap("- gamma is not pinned down: the teacher sample (about 150 persons) gives a ratio SE of 0.3 and the point estimate moves with the teacher definition; no larger NLSY79 teacher sample exists.")
    ap("- Uncertainty in the EIV structure (R, equal error variance across groups, independence from BA) is not in the SEs; the convention and EIV columns bracket the effect of that choice.")
    ap("- No saved bootstrap draws (only SEs and percentile CIs), so the covariances of these moments with the D3 and teacher-gap moments (same draws) are not provided.")
    ap("")
    path.write_text("\n".join(L) + "\n")


def toml_str(x: str) -> str:
    return json.dumps(x, ensure_ascii=True)


def write_toml(P: dict[str, Any], path: Path) -> None:
    R = Rep(P)
    M, d3 = P["meta"], P["d3"]
    L: list[str] = []
    ap = L.append
    ap("# nlsy_ba_slopes_targets.toml -- NLSY79 inputs for the joint estimation of the ability scale s_z with the college decision (DESIGN_v2.md section M3, task D3b)")
    ap(f"# and for checking the teacher wage curvature gamma. Generated {M['generated']} by data/spatial/nlsy_ba_slopes.py; see nlsy_ba_slopes.md (definitions, SEs, caveats).")
    ap("# role: target (M3), target_by_gender, validation. value, se: point estimate and household-bootstrap standard error (B = %d, seed %d)." % (M["B"], M["seed"]))
    ap("# model_moment: the name of the model-side moment (M3 spec) that the data moment is the counterpart of (\"none\" where the model has no counterpart yet).")
    ap("# convention: 'obs/sqrt(R)' = observed slope / sqrt(R_comp_pop), the project's latent slope (the pipeline of the 0.2127 wage_slope target); it is the `value` of every slope in this file.")
    ap("# 'eiv' = classical errors-in-variables correction of the regression itself (error variance 1 - R_comp_pop in afqt_c units), carried in the alt_eiv* fields. A BA-conditional and an")
    ap("# unconditional slope must be on one convention: wage_slope_ba (obs/sqrt(R)) goes with the wage_slope 0.2127 (0.0077); the EIV pair is wage_slope_ba alt_eiv with alt_unconditional_eiv")
    ap("# (it needs an amendment of the pre-registered V12 before any M3 result). alt_*: alternative definitions or versions. Nothing here is fitted; no value was chosen to match a model output.")
    ap("")
    ap("[source]")
    ap('nlsy79 = "NLSY79 block A, ages 25-34 wage sample of nlsy_ability.py block C (>= 48 weeks, >= 30 hours/week, wage income >= $1,000 in 2010 dollars, highest grade >= 9), 1979 weights, household bootstrap"')
    ap('acs = "ACS 2009-13 BA shares of all persons and of non-teachers (education_targets.toml ba_share_all and ba_share_nonteacher, T2 convention for the latter)"')
    ap("")
    ap("[meta]")
    ap(f"generated = {toml_str(M['generated'])}")
    ap('script = "data/spatial/nlsy_ba_slopes.py"')
    ap(f"nlsy_B = {M['B']}")
    ap(f"nlsy_seed = {M['seed']}")
    ap(f"eta = {ETA}")
    ap(f"gamma_baseline = {GAMMA_BASELINE}")
    ap('value_convention = "obs/sqrt(R)"')
    ap('ba_proxy = "highest grade completed >= 16 at some round at 25-34 (the extract has no degree variable)"')
    ap('sealed_set = "not used; NLSY97 not read (see nlsy_ba_slopes.md, Definitions and caveats)"')
    ap("")

    def num(x: float) -> str:
        return repr(float(x))

    def fmt(v: Any) -> Optional[str]:
        if isinstance(v, (bool, np.bool_)):
            return "true" if v else "false"
        if isinstance(v, (float, np.floating)):
            return num(v) if np.isfinite(v) else None
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        if isinstance(v, str):
            return toml_str(v)
        return str(v)

    def emit(name: str, head: dict[str, Any], extra: Optional[dict[str, Any]] = None) -> None:
        ap(f"[moment.{name}]")
        for k, v in {**head, **(extra or {})}.items():
            s = fmt(v)
            if s is not None:
                ap(f"{k} = {s}")
        ap("")

    def moment(name: str, role: str, key: str, units: str, convention: str, model_moment: str, definition: str, model: str,
               extra: Optional[dict[str, Any]] = None, scale: Optional[float] = None, scale_rule: Optional[str] = None,
               n_key: Optional[str] = None, counts: Optional[dict[str, int]] = None) -> None:
        r = R.r(key)
        if r is None:
            return
        head: dict[str, Any] = {"role": role, "value": r["value"], "se": r["se"]}
        if scale is not None:
            head["scale"] = scale
            head["scale_rule"] = scale_rule or ""
        c = counts if counts is not None else (R.cnt(n_key) if n_key is not None else None)
        if c is not None:
            head["n_person_years"], head["n_persons"] = c["person_years"], c["persons"]
        head.update({"units": units, "convention": convention, "source": "nlsy79", "model_moment": model_moment,
                     "definition": definition, "model": model})
        emit(name, head, extra)

    slope_units = "log hourly-wage points per SD of latent AFQT"
    CONV = "obs/sqrt(R)"
    versions = [R.val("ba_pooled|naive"), R.val("ba_pooled|eiv"), R.val("iv_ba|mean4")]
    wage_slope_ba_scale = float(np.sqrt(R.r("ba_pooled|naive")["se"] ** 2 + (0.5 * (max(versions) - min(versions))) ** 2))
    se_rule = "bootstrap SE, as the existing wage_slope target (0.2127, SE 0.0077)"
    mod_ba = ("weighted OLS slope of log wage on the standardized adult score with a BA (d = 1) dummy{fem}, over market workers (teachers included, home producers excluded); "
              "wage_slope_ba in spatial_moments.jl (spec M3_SPEC_v2.md)")
    moment("wage_slope_ba", "target", "ba_pooled|naive", slope_units, CONV, "wage_slope_ba",
           "latent wage-AFQT slope in a regression of the log hourly wage on latent AFQT, a BA indicator (highest grade >= 16 at 25-34), a female dummy and age and survey-year dummies; "
           "NLSY79 wage sample, 1979 weights; latent slope = observed slope / sqrt(R_comp_pop), the pipeline of the 0.2127 wage_slope target",
           mod_ba.format(fem=" and a female dummy"),
           {"alt_eiv": R.val("ba_pooled|eiv"), "alt_eiv_se": R.r("ba_pooled|eiv")["se"],
            "alt_split_half_iv": R.val("iv_ba|mean4"), "alt_split_half_iv_se": R.r("iv_ba|mean4")["se"],
            "alt_unconditional_convention": R.val("u_pooled|naive"), "alt_unconditional_eiv": R.val("u_pooled|eiv"), "alt_unconditional_eiv_se": R.r("u_pooled|eiv")["se"],
            "alt_ba_at_round": R.val("bar_pooled|naive"), "alt_ba_at_round_eiv": R.val("bar_pooled|eiv"),
            "ba_coefficient_eiv": R.val("ba_pooled|ba_coef_eiv"), "ba_coefficient_uncorrected": R.val("ba_pooled|ba_coef_naive"),
            "conditional_reliability": R.val("ba_pooled|lam"),
            "pairing": "same convention as the wage_slope target %.4f (0.0077): %.4f pairs with %.4f. If EIV is preferred, switch BOTH to the EIV values (alt_eiv with alt_unconditional_eiv); that needs an amendment of the pre-registered V12 before any M3 result"
                       % (R.val("u_pooled|naive"), R.val("ba_pooled|naive"), R.val("u_pooled|naive")),
            "note": "convention value; with the BA control the reliability of afqt_c given the controls falls to the conditional_reliability, so EIV and the split-half IV are %.0f%% and %.0f%% higher; the three versions span %.4f (SE %.4f)"
                    % (100 * (R.val("ba_pooled|eiv") / R.val("ba_pooled|naive") - 1), 100 * (R.val("iv_ba|mean4") / R.val("ba_pooled|naive") - 1),
                       max(versions) - min(versions), R.r("ba_pooled|naive")["se"])},
           scale=wage_slope_ba_scale, scale_rule="bootstrap SE in quadrature with half the spread of the three disattenuation versions (convention, split-half IV, EIV); the existing wage_slope target uses the bare SE (0.0077)",
           n_key="ba_pooled")
    for g, nm in (("male", "m"), ("female", "f")):
        moment(f"wage_slope_ba_{nm}", "target_by_gender", f"ba_{g}|naive", slope_units, CONV, "wage_slope_ba",
               f"as wage_slope_ba, {'men' if g == 'male' else 'women'} only (no female dummy)", mod_ba.format(fem="") + f", {'men' if g == 'male' else 'women'} only",
               {"alt_eiv": R.val(f"ba_{g}|eiv"), "alt_eiv_se": R.r(f"ba_{g}|eiv")["se"], "alt_unconditional_convention": R.val(f"u_{g}|naive"),
                "alt_unconditional_eiv": R.val(f"u_{g}|eiv"), "conditional_reliability": R.val(f"ba_{g}|lam")},
               scale=R.r(f"ba_{g}|naive")["se"], scale_rule=se_rule, n_key=f"ba_{g}")
    reg_n = {"n_person_years_regression": R.cnt("int_pooled")["person_years"], "n_persons_regression": R.cnt("int_pooled")["persons"]}
    moment("wage_slope_within_nonba", "validation", "int_pooled|nonba_naive", slope_units, CONV, "none",
           "slope within non-BA workers in the regression with BA, AFQT x BA, a female dummy and the common controls; n_* are the non-BA workers' counts, n_*_regression those of the whole interacted regression",
           "the same slope over d = 0 market workers (no model counterpart yet)",
           {"alt_eiv": R.val("int_pooled|nonba_eiv"), "alt_eiv_se": R.r("int_pooled|nonba_eiv")["se"],
            "alt_separate_samples": R.val("sepnba_pooled|naive"), "alt_separate_samples_eiv": R.val("sepnba_pooled|eiv"), **reg_n},
           counts=R.cnt("sepnba_pooled"))
    moment("wage_slope_within_ba", "validation", "int_pooled|ba_naive", slope_units, CONV, "none",
           "slope within BA workers in the same interacted regression; n_* are the BA workers' counts, n_*_regression those of the whole interacted regression",
           "the same slope over d = 1 market workers (teachers included; no model counterpart yet)",
           {"alt_eiv": R.val("int_pooled|ba_eiv"), "alt_eiv_se": R.r("int_pooled|ba_eiv")["se"],
            "alt_separate_samples": R.val("sepba_pooled|naive"), "alt_separate_samples_eiv": R.val("sepba_pooled|eiv"),
            "diff_ba_minus_nonba": R.val("int_pooled|diff_naive"), "diff_se": R.r("int_pooled|diff_naive")["se"],
            "diff_ba_minus_nonba_eiv": R.val("int_pooled|diff_eiv"), "diff_se_eiv": R.r("int_pooled|diff_eiv")["se"], **reg_n},
           counts=R.cnt("sepba_pooled"))
    moment("ba_coefficient_given_afqt", "validation", "ba_pooled|ba_coef_eiv", "log hourly-wage points", "eiv", "none",
           "coefficient of the BA indicator in the wage_slope_ba regression (premium at given LATENT AFQT, so errors-in-variables corrected: the uncorrected coefficient mixes in the ability the noisy score leaves out)",
           "BA minus non-BA log wage at the same z over market workers (credential premium rho_c plus the return to the extra schooling time), ability held fixed; no model counterpart yet",
           {"alt_uncorrected": R.val("ba_pooled|ba_coef_naive"), "alt_uncorrected_se": R.r("ba_pooled|ba_coef_naive")["se"]}, n_key="ba_pooled")

    # ---- teachers against non-teacher BA workers: the BA-controlled slope is the headline; the first-version (no control) and BA-teacher rows sit beside it
    rat_defs = (("main", "teacher_slope_ratio", "gamma_implied",
                 "slope of the log wage on latent AFQT among all K-12 teachers, in a regression with a BA indicator (the slope of wage_slope_ba), over the same slope among non-teacher BA workers"),
                ("tba", "teacher_slope_ratio_ba_teachers", "gamma_implied_ba_teachers",
                 "slope of the log wage on latent AFQT among K-12 teachers with a BA over the same slope among non-teacher BA workers"),
                ("nobac", "teacher_slope_ratio_no_ba_control", "gamma_implied_no_ba_control",
                 "slope of the log wage on latent AFQT among all K-12 teachers WITHOUT a BA control (it contains the credential composition of the teachers) over the same slope among non-teacher BA workers (the first version of this file)"))
    sens_tags = ("core", "half", "ever", "any", "c70", "round", "deflage", "deflnone")
    for tg, nm, gnm, desc in rat_defs:
        num_k, den_k, dsc = P["ratios"][tg]
        dm_n, dm_e = P["delta_method"][tg]["naive"], P["delta_method"][tg]["eiv"]
        ex = {"alt_eiv": R.val(f"r_{tg}|eiv"), "alt_eiv_se": R.r(f"r_{tg}|eiv")["se"],
              "ci95_lo": R.r(f"r_{tg}|naive")["ci95_lo"], "ci95_hi": R.r(f"r_{tg}|naive")["ci95_hi"],
              "teacher_slope": R.val(f"{num_k}|naive"), "teacher_slope_se": R.r(f"{num_k}|naive")["se"],
              "comparison_slope": R.val(f"{den_k}|naive"), "comparison_slope_se": R.r(f"{den_k}|naive")["se"],
              "teacher_slope_eiv": R.val(f"{num_k}|eiv"), "teacher_slope_eiv_se": R.r(f"{num_k}|eiv")["se"],
              "comparison_slope_eiv": R.val(f"{den_k}|eiv"), "comparison_slope_eiv_se": R.r(f"{den_k}|eiv")["se"],
              "n_teacher_persons": R.cnt(num_k)["persons"], "n_teacher_person_years": R.cnt(num_k)["person_years"],
              "n_comparison_persons": R.cnt(den_k)["persons"], "n_comparison_person_years": R.cnt(den_k)["person_years"],
              "ba_control_in_teacher_slope": tg == "main", "delta_method_se": dm_n["se_ratio"], "delta_method_se_eiv": dm_e["se_ratio"]}
        if tg == "main":
            ex.update({"alt_no_ba_control": R.val("r_nobac|naive"), "alt_no_ba_control_se": R.r("r_nobac|naive")["se"],
                       "alt_ba_teachers_only": R.val("r_tba|naive"), "alt_ba_teachers_only_se": R.r("r_tba|naive")["se"],
                       "alt_vs_all_nonteachers_ba_control": R.val("r_ntc|naive"), "alt_vs_all_nonteachers_ba_control_se": R.r("r_ntc|naive")["se"],
                       "baseline_implied_ratio": ratio_of_gamma(GAMMA_BASELINE)})
        ex["note"] = "small teacher sample (147 persons): a weak check, not a target"
        moment(nm, "validation", f"r_{tg}|naive", "ratio", CONV, "teacher_slope_ratio",
               f"{desc} (nlsy_teachers.py spec {MAIN_SPEC}, person-level majority rule); {dsc}",
               "ratio of teachers' to non-teacher BA workers' log-wage slopes on log z in the model (not yet a model moment); closed form gamma (1 - eta) / (1 - gamma eta), which has no credential margin",
               ex)
        exg = {"alt_eiv": R.val(f"gamma_{tg}|eiv"), "alt_eiv_se": R.r(f"gamma_{tg}|eiv")["se"],
               "ci95_lo": R.r(f"gamma_{tg}|naive")["ci95_lo"], "ci95_hi": R.r(f"gamma_{tg}|naive")["ci95_hi"],
               "eta": ETA, "baseline_gamma": GAMMA_BASELINE, "baseline_implied_ratio": ratio_of_gamma(GAMMA_BASELINE),
               "delta_method_se": dm_n["se_gamma"], "delta_method_se_eiv": dm_e["se_gamma"]}
        if tg == "main":
            exg.update({"range_over_definitions_lo": min(R.val(f"gamma_{t}|naive") for t in sens_tags),
                        "range_over_definitions_hi": max(R.val(f"gamma_{t}|naive") for t in sens_tags),
                        "alt_no_ba_control": R.val("gamma_nobac|naive"), "alt_no_ba_control_se": R.r("gamma_nobac|naive")["se"],
                        "alt_ba_teachers_only": R.val("gamma_tba|naive"), "alt_ba_teachers_only_se": R.r("gamma_tba|naive")["se"]})
        moment(gnm, "validation", f"gamma_{tg}|naive", "gamma", CONV, "none",
               f"teacher wage curvature solving ratio = gamma (1 - eta) / (1 - gamma eta) at eta = {ETA} for the slope ratio of {nm}",
               "gamma (internal parameter, Table 3), checked against, not targeted", exg)

    # ---- BA-share logit shift: LIKE FOR LIKE.  Primary: all persons (teachers included) to the ACS all-persons share; alternatives: non-teachers to the ACS
    # non-teacher share; and the MIXED version (all persons to the non-teacher share), kept as a flagged alternative.
    n_all = f"all persons of the D3 primary sample (persons observed at 25-34, highest grade >= 9; teachers included; NLSY79 share {R.val('ba|share'):.4f})"
    mixed_note = ("NOT like for like: the NLSY79 all-persons share contains the K-12 teachers (%.0f%% with a BA), the ACS non-teacher share does not; this is the version of the first release of this file "
                  "and of the task text, and it understates the shift by %.2f" % (100 * R.val("ba|teacher_ba_share"), R.val("ba|shift_cell_all") - R.val("ba|shift_cell_mixed")))
    shift_extra = {"alt_individual_logit": R.val("ba|shift_indiv_all"), "alt_individual_logit_se": R.r("ba|shift_indiv_all")["se"],
                   "alt_nonteachers": R.val("banT|shift_cell_nt"), "alt_nonteachers_se": R.r("banT|shift_cell_nt")["se"],
                   "alt_nonteachers_individual_logit": R.val("banT|shift_indiv_nt"), "alt_nonteachers_individual_logit_se": R.r("banT|shift_indiv_nt")["se"],
                   "alt_mixed_population": R.val("ba|shift_cell_mixed"), "alt_mixed_population_se": R.r("ba|shift_cell_mixed")["se"],
                   "alt_mixed_population_individual_logit": R.val("ba|shift_indiv_mixed"),
                   "mixed_population_note": mixed_note,
                   "target_share": d3["ba_share_all"], "target_share_se": d3["ba_share_all_se"], "nlsy79_share": R.val("ba|share"),
                   "nlsy79_share_nonteachers": R.val("banT|share"), "target_share_nonteachers": d3["ba_share_nonteacher"],
                   "nlsy79_teacher_weight_share": R.val("ba|teacher_wshare"), "nlsy79_teacher_ba_share": R.val("ba|teacher_ba_share"),
                   "shift_se_from_acs_target": R.val("ba|dshift_cell_all") * d3["ba_share_all_se"]}
    moment("ba_logit_shift", "validation", "ba|shift_cell_all", "log-odds", "cell shift of the weighted D3 tercile shares", "ba_share_by_z",
           f"common logit shift of the D3 tercile shares of P(BA | AFQT tercile) of {n_all} that moves the weighted NLSY79 BA share to the ACS 2009-13 all-persons BA share "
           f"{d3['ba_share_all']:.4f} (education_targets.toml ba_share_all): like for like (all persons to all persons)",
           "cohort adjustment of the BA level before comparing P(d = 1 | tercile of z) over ALL adults (teachers and non-teachers, home production included) with the data; the model's ba_share and ba_share_by_z are over all adults",
           shift_extra)
    for q in (1, 2, 3):
        ex = {"nlsy79_unshifted": R.val(f"ba|terc{q}"), "alt_individual_logit": R.val(f"ba|imp_indiv_all_terc{q}"), "alt_individual_logit_se": R.r(f"ba|imp_indiv_all_terc{q}")["se"],
              "alt_nonteachers": R.val(f"banT|imp_cell_nt_terc{q}"), "alt_nonteachers_se": R.r(f"banT|imp_cell_nt_terc{q}")["se"],
              "alt_mixed_population": R.val(f"ba|imp_cell_mixed_terc{q}"), "alt_mixed_population_se": R.r(f"ba|imp_cell_mixed_terc{q}")["se"],
              "mixed_population_note": "NOT like for like (all persons shifted to the non-teacher share); see ba_logit_shift"}
        moment(f"p_ba_afqt_tercile_{q}_at_acs_share", "validation", f"ba|imp_cell_all_terc{q}", "share", "weighted, cell-shifted to the ACS share", "ba_share_by_z",
               f"share with 16+ years of schooling in tercile {q} of own AFQT (observed score; cut points {R.val('ba|cut1'):.3f}, {R.val('ba|cut2'):.3f}) among ALL persons of the D3 primary sample, the D3 NLSY79 share shifted by "
               f"ba_logit_shift to the ACS all-persons BA share {d3['ba_share_all']:.4f} (like for like)",
               "P(d = 1 | tercile of the observed score x = sqrt(R) theta + sqrt(1 - R) u, theta = standardized log z, R = 0.916; terciles of x over ALL adults, teachers and non-teachers) at the model's own BA share of all adults",
               ex)
    moment("p_ba_afqt_top_minus_bottom_at_acs_share", "validation", "ba|imp_cell_all_diff", "share difference", "weighted, cell-shifted to the ACS share", "ba_share_by_z",
           "top minus bottom AFQT tercile of the shifted shares (all persons, like for like)", "as p_ba_afqt_tercile_q_at_acs_share",
           {"nlsy79_unshifted": R.val("ba|terc_diff"), "alt_individual_logit": R.val("ba|imp_indiv_all_diff"),
            "alt_nonteachers": R.val("banT|imp_cell_nt_diff"), "alt_nonteachers_se": R.r("banT|imp_cell_nt_diff")["se"],
            "alt_mixed_population": R.val("ba|imp_cell_mixed_diff"), "alt_mixed_population_se": R.r("ba|imp_cell_mixed_diff")["se"],
            "mixed_population_note": "NOT like for like (all persons shifted to the non-teacher share); see ba_logit_shift"})
    path.write_text("\n".join(L))


# =============================================================================
# Self-test
# =============================================================================

def selftest() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'ok' if cond else 'FAIL'}] {name}")

    # gamma formula
    check("gamma <-> ratio are inverse", all(abs(ratio_of_gamma(gamma_of_ratio(r)) - r) < 1e-12 for r in (0.0, 0.3, 0.756, 1.0, 1.4)))
    check("gamma = 1 <=> ratio = 1; ratio(0.6194) = 0.5996", abs(gamma_of_ratio(1.0) - 1.0) < 1e-12 and abs(ratio_of_gamma(0.6194) - 0.5996) < 5e-4)
    # logit shift
    rng = np.random.default_rng(7)
    eta = rng.normal(-1.0, 1.2, 5000)
    w = rng.uniform(0.5, 2.0, 5000)
    d = logit_shift(eta, w, 0.40)
    check("logit_shift hits the target", abs((w * lam_fn(eta + d)).sum() / w.sum() - 0.40) < 1e-10)
    check("logit_shift is 0 at the current share", abs(logit_shift(eta, w, float((w * lam_fn(eta)).sum() / w.sum()))) < 1e-8)
    # weighted logit recovers (a, b)
    x = rng.normal(0, 1, 200000)
    yb = (rng.uniform(size=x.size) < lam_fn(-1.0 + 1.5 * x)).astype(float)
    a_, b_ = wlogit(x, yb, np.ones(x.size))
    check("weighted logit recovers (a, b)", abs(a_ + 1.0) < 0.03 and abs(b_ - 1.5) < 0.03)
    # errors-in-variables: latent f with var R, error u with var 1 - R, BA correlated with f
    n, Rl = 600000, 0.80
    f = rng.normal(0, np.sqrt(Rl), n)
    u = rng.normal(0, np.sqrt(1 - Rl), n)
    a = f + u                                                   # var(afqt) = 1
    ba = (1.2 * f + rng.normal(0, 0.6, n) > 0.4).astype(float)
    b0, bba, bint = 0.25, 0.20, -0.05                            # per latent SD (SD of f = sqrt(Rl))
    sdf = np.sqrt(Rl)
    y = b0 * f / sdf + bba * ba + bint * (f / sdf) * ba + rng.normal(0, 0.3, n)
    ww = rng.uniform(0.5, 1.5, n)
    one = np.ones(n)
    bo, be, lam = eiv_fit(y, np.column_stack([one, a]), ww, {1: one}, 1 - Rl)
    check("EIV, no control: reproduces OLS / sqrt(R) for var(afqt) = 1 (the project's convention)",
          abs(be[1] * sdf - bo[1] / sdf) < 0.004 and abs(lam - Rl) < 0.004)
    lam_direct = 1 - (1 - Rl) / np.var(a - np.polyval(np.polyfit(ba, a, 1), ba))
    y2 = b0 * f / sdf + bba * ba + rng.normal(0, 0.3, n)           # no interaction: BA-control version
    bo, be, lam = eiv_fit(y2, np.column_stack([one, a, ba]), ww, {1: one}, 1 - Rl)
    check("EIV with a BA control recovers the latent slope and BA coefficient; OLS / sqrt(R) and OLS do not",
          abs(be[1] * sdf - b0) < 0.012 and abs(be[2] - bba) < 0.012 and abs(bo[1] / sdf - b0) > 0.012 and abs(bo[2] - bba) > 0.012)
    check("conditional reliability matches 1 - var(u) / var(afqt | BA)", abs(lam - lam_direct) < 0.01 and lam < Rl)
    bo, be, lam = eiv_fit(y, np.column_stack([one, a, ba, a * ba]), ww, {1: one, 3: ba}, 1 - Rl)
    check("EIV interacted: slopes within non-BA and BA recovered",
          abs(be[1] * sdf - b0) < 0.02 and abs((be[1] + be[3]) * sdf - (b0 + bint)) < 0.02)
    # rank deficiency: an all-zero dummy does not break the solver
    z = np.zeros(n)
    bo, be, lam = eiv_fit(y2, np.column_stack([one, a, ba, z]), ww, {1: one}, 1 - Rl)
    check("rank-deficient design (empty dummy) solves", np.isfinite(be).all() and abs(be[1] * sdf - b0) < 0.012 and be[3] == 0.0)
    # BA-share populations: shifts hit their own targets; a population without the high-BA teachers has the lower share
    m_ = 30000
    xs = rng.normal(0, 1, m_)
    wv = rng.uniform(0.5, 2.0, m_)
    teach = (rng.uniform(size=m_) < 0.06) & (xs > -0.6)
    bav = (rng.uniform(size=m_) < lam_fn(-1.2 + 1.8 * xs)) | (teach & (rng.uniform(size=m_) < 0.9))
    tq = np.digitize(xs, np.quantile(xs, [1 / 3, 2 / 3]))
    o_all = ba_pop_block("p", xs, wv, np.ones(m_, bool), bav, tq, {"a": 0.40, "b": 0.30})
    o_nt = ba_pop_block("n", xs, wv, ~teach, bav, tq, {"a": 0.40})
    check("ba_pop_block: the cell and the individual-logit shift both hit the target share",
          abs(sum(o_all[f"p|pi{q}"] * o_all[f"p|imp_cell_a_terc{q}"] for q in (1, 2, 3)) - 0.40) < 1e-9 and abs(o_all["p|imp_indiv_a_share"] - 0.40) < 1e-9
          and abs(sum(o_nt[f"n|pi{q}"] * o_nt[f"n|imp_cell_a_terc{q}"] for q in (1, 2, 3)) - 0.40) < 1e-9)
    check("ba_pop_block: tercile weights sum to one, a lower target needs a smaller shift, and teachers (high BA) lift the all-persons share",
          abs(sum(o_all[f"p|pi{q}"] for q in (1, 2, 3)) - 1.0) < 1e-12 and o_all["p|shift_cell_b"] < o_all["p|shift_cell_a"] and o_nt["n|share"] < o_all["p|share"])
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


# =============================================================================
# Orchestration
# =============================================================================

def wmean_b(x: np.ndarray, w: np.ndarray) -> float:
    return float((x * w).sum() / w.sum())


def point_meta(dat: Data) -> dict[str, Any]:
    """Counts and sample facts at the point estimate (unweighted counts, weighted shares and SDs)."""
    fr, D, py = dat.fr, dat.D, dat.py
    afqt, Rc = na.pop_afqt(D.z, D.asamp, D.w)[:2]
    a_py, w_py = afqt[py.p], D.w[py.p]
    rows = dat.base & np.isfinite(a_py) & (w_py > 0)
    sp = dat.specs[MAIN_SPEC]
    tp = rows & sp.teacher[py.p]
    hg = fr.hg_max
    m = fr.observed & D.asamp & np.isfinite(hg) & (hg >= na.MIN_HGC) & np.isfinite(afqt) & (D.w > 0)

    def sd(mask: np.ndarray) -> float:
        return float(np.sqrt(na.wvar(a_py[mask], w_py[mask])))

    return {"ba_share_wage": wmean_b(dat.ba_py[rows].astype(float), w_py[rows]),
            "sd_afqt_wage": sd(rows), "sd_afqt_ba": sd(rows & dat.ba_py), "sd_afqt_nba": sd(rows & ~dat.ba_py),
            "teacher_rows_other_occ_share": float(1.0 - (tp & fr.k12["c80"]).sum() / tp.sum()),
            "n_teacher_persons": int(np.unique(py.p[tp]).size),
            "n_teacher_women": int(np.unique(py.p[tp & dat.female_py]).size),
            "n_teacher_men": int(np.unique(py.p[tp & ~dat.female_py]).size),
            "n_ba_sample": int(m.sum()), "n_wage_rows": int(rows.sum()), "n_wage_persons": int(np.unique(py.p[rows]).size),
            "R_comp_pop_point": float(Rc)}


def d3_check(stats: dict[str, dict[str, float]], d3: dict[str, Any]) -> dict[str, float]:
    s = d3["stats"]
    terc = max(abs(stats[f"ba|terc{q}"]["value"] - s[f"terc|all_hgc9|pooled|ba_share|{q}"]["value"]) for q in (1, 2, 3))
    lg = max(abs(stats["ba|logit_a"]["value"] - s["logit|all_hgc9|pooled|intercept"]["value"]),
             abs(stats["ba|logit_b"]["value"] - s["logit|all_hgc9|pooled|slope_obs"]["value"]))
    return {"max_abs_diff_terc_shares": float(terc), "max_abs_diff_logit": float(lg),
            "diff_shift_all": float(abs(stats["ba|shift_indiv_all"]["value"] - s["logit|all_hgc9|pooled|shift"]["value"])),
            "max_abs_diff_cuts": float(max(abs(stats["ba|cut1"]["value"] - s["terc|all_hgc9|cut1"]["value"]),
                                           abs(stats["ba|cut2"]["value"] - s["terc|all_hgc9|cut2"]["value"]))),
            "diff_ba_share": float(abs(stats["ba|share"]["value"] - s["all_hgc9|pooled|ba_share"]["value"]))}


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bootstrap", type=int, default=B_DEFAULT)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--quick", action="store_true", help="B = 20 (and a 20-draw replay)")
    p.add_argument("--no-replay", action="store_true", help="skip the replay of the published bootstrap stream")
    p.add_argument("--md-only", action="store_true", help="rewrite the .md and .toml from the saved .json")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--outdir", default=str(OUTDIR_DEFAULT))
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    if args.selftest:
        return selftest()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if args.md_only:
        P = json.loads((outdir / "nlsy_ba_slopes.json").read_text())
        write_md(P, outdir / "nlsy_ba_slopes.md")
        write_toml(P, outdir / "nlsy_ba_slopes_targets.toml")
        return 0
    B = 20 if args.quick else args.bootstrap
    if args.quick and outdir == OUTDIR_DEFAULT.resolve():
        p.error("--quick results are not valid inference; pass an explicit --outdir so the published files are not overwritten")
    t0 = time.time()
    c = na.Cohort("nlsy79", HERE)
    fr = nt.Frame(c)
    dat = Data(fr)
    pub, d3 = read_published(), read_d3()
    # like-for-like targets: all persons (teachers included) to the ACS all-persons share, the non-teachers to the ACS non-teacher share;
    # "mixed" = all persons to the non-teacher share (the version the task text asked for; NOT like-for-like, kept as a flagged alternative)
    targets = {"all": {"all": d3["ba_share_all"], "mixed": d3["ba_share_nonteacher"]}, "nt": {"nt": d3["ba_share_nonteacher"]}}
    LOG.info("== bootstrap B=%d ==", B)
    res = run_bootstrap(dat, B, args.seed, targets)
    S = res["stats"]
    LOG.info("pooled unconditional latent slope %.5f (published %.5f), BA control %.4f (convention) / %.4f (EIV)",
             S["u_pooled|naive"]["value"], pub["slopes"]["pooled"]["value"], S["ba_pooled|naive"]["value"], S["ba_pooled|eiv"]["value"])
    if args.no_replay:
        replay = {g: {"se": float("nan"), "mean": float("nan"), "B": 0, "seed": PUB_SEED} for g in GROUPS}
    else:
        LOG.info("== replay of the published bootstrap stream (B=%d, seed %d) ==", 20 if args.quick else PUB_B, PUB_SEED)
        replay = replay_published(dat, 20 if args.quick else PUB_B, PUB_SEED, pub["inputs"]["nlsy97_ability"]["rows"])
    meta = point_meta(dat)
    meta.update({"generated": dt.datetime.now().isoformat(timespec="seconds"), "script": "data/spatial/nlsy_ba_slopes.py",
                 "B": B, "seed": args.seed, "failed_draws": res["failed_draws"], "eta": ETA, "gamma_baseline": GAMMA_BASELINE,
                 "main_spec": MAIN_SPEC, "college_hgc": COLLEGE_HGC, "age_window": [na.AGE_LO, na.AGE_HI],
                 "years": sorted(set(dat.py.Y.tolist())), "quick": bool(args.quick)})
    delta = {tg: {sc: delta_method_gamma(res, num, den, sc) for sc in ("naive", "eiv")} for tg, (num, den, _) in dat.ratios.items()
             if f"{num}|eiv" in S and f"{den}|eiv" in S}
    P = clean({
        "meta": meta, "published": pub, "replay": replay,
        "d3": {**{k: d3[k] for k in ("ba_share_nonteacher", "ba_share_nonteacher_se", "ba_share_all", "ba_share_all_se", "acs_teacher_share", "teacher_gap_vs_ba_nonteachers", "teacher_gap_within_ba", "stats", "generated")},
               "check": d3_check(S, d3)},
        "specs": {rg.key: rg.desc for rg in dat.regs},
        "ratios": {tg: list(v) for tg, v in dat.ratios.items()},
        "counts": res["counts"], "statistics": S, "delta_method": delta,
        "definitions": {"statistic_keys": {
            "{reg}|obs": "coefficient of afqt_c (log points per afqt_c SD)",
            "{reg}|naive": "convention: obs / sqrt(R_comp_pop), per SD of latent AFQT",
            "{reg}|eiv": "classical errors-in-variables corrected slope per SD of latent AFQT",
            "{reg}|lam": "conditional reliability of afqt_c given the regression's controls",
            "{reg}|ba_coef_naive / ba_coef_eiv": "coefficient of the BA indicator (OLS, EIV-corrected)",
            "int_*|nonba_*, ba_*, diff_*": "slopes within non-BA, within BA and BA minus non-BA in the interacted regression",
            "r_{tag}|{naive,eiv}, gamma_{tag}|{naive,eiv}": "teacher / non-teacher-BA slope ratio and gamma = ratio / (1 - eta + eta ratio); tag main = teachers with a BA control, nobac = without, tba = BA teachers only, ntc/nt = against all non-teachers",
            "ba|...": "ALL persons of D3's primary sample: BA share, tercile cuts and shares (D3 definitions), logit (a, b), level shifts to the ACS targets 'all' (all-persons share 0.3391; like for like) and "
                      "'mixed' (non-teacher share 0.3231; NOT like for like): shift_cell_*, shift_indiv_*, imp_cell_*_terc{q}, imp_indiv_*_terc{q}, fit_indiv_terc{q} (unshifted logit fit), dshift_* (d delta / d target); "
                      "teacher_wshare, teacher_ba_share (the main-spec teachers among them)",
            "banT|...": "the NON-TEACHERS of the same sample (everyone not a main-spec teacher), cut at the all-persons terciles: share, terc{q}, pi{q} (tercile weights), logit (a, b), level shifts to the ACS non-teacher share ('nt')"},
            "regressions": {rg.key: {"kind": rg.kind, "female_dummy": rg.female, "deflated": rg.deflate, "controls": rg.controls} for rg in dat.regs}}})
    (outdir / "nlsy_ba_slopes.json").write_text(json.dumps(P, indent=1, default=float, allow_nan=False))
    write_md(P, outdir / "nlsy_ba_slopes.md")
    write_toml(P, outdir / "nlsy_ba_slopes_targets.toml")
    LOG.info("done in %.1fs -> %s", time.time() - t0, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
