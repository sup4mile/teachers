#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
education_occupation.py -- evidence for a joint education (BA) and occupational choice: the data side of
DESIGN.md section M3 (task D3 of the 2026-10b model experiments).

M3 lets the young choose college (d = 1, a BA, which teaching requires) before the occupational shocks are
drawn: a log credential premium rho_c on non-teaching wages, a utility cost c_C of college and its gradient
c_Cz in log z, and phi re-fit to non-teachers' schooling.  This script measures what disciplines them and what
validates them, with standard errors, and writes a report, a JSON and a TOML for the Julia calibration.

  (a) ACS PUMS 2009-13 (the 2013 5-year file cached by `acs_occupations.py`; its sample builder, HHJK rules,
      crosswalk and K-12 codes are reused unchanged).  Ages 25-34, the Table 2 universe with home production,
      the wage sample for wages, years of schooling as in the Mincer workbook (capped 9-17), by gender and
      pooled: BA+ (SCHL >= 21) shares of everyone, of non-teachers and of K-12 teachers; graduate-degree share of
      teachers; K-12 teachers by Census code; mean years by BA status; the log-wage premium of BA+ over non-BA
      among non-teacher market workers (raw, age-controlled, within occupation); teachers' mean log wage minus
      all non-teacher market workers' (the existing targets +0.0557 / +0.0860: cell means weighted by the Table 2
      counts, reproduced) and minus non-teacher BA+ workers' (same convention), for all K-12 teachers, for BA+
      teachers (the like-for-like counterpart of a model whose teachers all have a BA) and for the core codes
      2310-2330; BA+ minus non-BA K-12 teachers; graduate-degree minus BA-only teachers; the share of full-time
      workers outside the Table 2 wage sample by weeks worked and the premia when the 40-47 (and 27-39) week
      brackets are admitted.  Standard errors: the Census successive-difference formula on the 80 replicate
      weights (`--fetch-repw` downloads them once, 34 MB), every statistic recomputed per replicate.
  (b) NCES Digest Table 211.20 (average base salary of full-time public school teachers by highest degree and
      years of experience; 2017-18 NTPS and the other editions): the master's over bachelor's premium at low
      experience, delta-method SEs.
  (c) NLSY79 (the `nlsy_teachers.py` Frame: latent AFQT, 1979 weights, household bootstrap B = 500 with the
      draws of that script): the latent AFQT gap of BA+ (16+ years) over non-BA at 25-34, all and employed;
      P(BA+ | AFQT tercile) on the observed scale, level-shifted to the ACS BA share, and on the latent scale
      (an errors-in-variables logit with the composite reliability); the existing teacher gaps (+0.614, -0.334)
      recomputed and compared with `nlsy_teachers.json`.  NLSY97 is not touched (reserved for the sealed
      hold-out, DESIGN.md D4).
  (d) Tuition: NCES Tables 330.10 (published prices), 331.30 and the College Board (net prices), 303.70
      (enrollment weights), total cost of a BA, and its conversion to the model's units: against CE mean annual
      expenditures per ADULT ($61,224 per consumer unit over 1.9 adults; per consumer unit and per earner are
      reported too), annual (the unit of a recurring cost) and total over working-life consumption, undiscounted
      and at 3% (the unit of a one-off cost), plus a unit-free aggregate-flow check from the CE college tuition line.

Two aggregations are reported side by side (see section 1 of the code): P, plain PWGTP-weighted population
means, and T2, the convention of the existing targets (cell means combined with the unweighted Table 2 counts).

    cd data/spatial
    python3 education_occupation.py --fetch-repw         # once: ACS replicate weights (about 4 minutes)
    python3 education_occupation.py                      # everything, about 10 minutes on a compute node
    python3 education_occupation.py --quick              # NLSY B = 20, 8 ACS replicates (SEs not valid)
    python3 education_occupation.py --md-only            # rewrite the .md and .toml from the saved .json
    python3 education_occupation.py --selftest           # unit checks, no data needed

Inputs: raw/acs_pums/acs5_2013/ and acs5_2013_repw/, estimates/acs_occupations_crosswalk.csv and
acs_occupations.json (reproduction check), raw/nlsy79_ability.parquet, estimates/settings/nlsy_teachers.json,
NCES pages, the College Board report and the BLS CE table (downloaded once into raw/education_occupation/).  The
Census key (CENSUS_API_KEY) is only used by --fetch-repw, through `acs_occupations.Api`, and is never printed or
written.
Writes estimates/settings/education_occupation.{md,json} and education_targets.toml.

Written with help from Claude Code.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acs_occupations as ao  # noqa: E402

LOG = logging.getLogger("education_occupation")

REPW_DIR = ao.RAW / "acs5_2013_repw"      # raw/acs_pums/acs5_2013_repw/<state>.parquet
N_REP = 80


# =============================================================================
# 0. ACS replicate weights (Census API), for design-based standard errors
# =============================================================================

def _fetch_repw_state(st: int) -> tuple[int, int]:
    """Download PWGTP1..PWGTP80 for ages 25-34 of one state of the 2013 5-year PUMS (two calls of
    40: the API allows 50 variables per call) and cache them as int32 in REPW_DIR/<state>.parquet."""
    f = REPW_DIR / f"{st:02d}.parquet"
    if f.exists():
        return st, len(pd.read_parquet(f, columns=["SPORDER"]))
    api = ao.Api()                                   # one session (and key scrubber) per task
    parts = []
    for lo in (1, 41):
        want = ["SERIALNO", "SPORDER", "PWGTP"] + [f"PWGTP{k}" for k in range(lo, lo + 40)]
        j = api.get_json(ao.API.format(year=2013, product="acs5"),
                         {"get": ",".join(want), "for": f"state:{st:02d}", "AGEP": f"{ao.AGE_MIN}:{ao.AGE_MAX}"})
        if not j:
            raise RuntimeError(f"state {st:02d}: empty response")
        df = pd.DataFrame(j[1:], columns=j[0])
        df = df.loc[:, ~df.columns.duplicated()]
        parts.append(df[want])
    a, b = parts
    if not (len(a) == len(b) and (a.SERIALNO.values == b.SERIALNO.values).all()
            and (a.SPORDER.values == b.SPORDER.values).all()):
        b = a[["SERIALNO", "SPORDER"]].merge(b, on=["SERIALNO", "SPORDER"], how="left")
        if len(b) != len(a):
            raise RuntimeError(f"state {st:02d}: the two replicate-weight calls do not line up")
    out = pd.concat([a, b.drop(columns=["SERIALNO", "SPORDER", "PWGTP"])], axis=1)
    rep = [f"PWGTP{k}" for k in range(1, N_REP + 1)]
    for c in ["PWGTP"] + rep:
        out[c] = pd.to_numeric(out[c], errors="raise").astype("int32")
    out["SPORDER"] = pd.to_numeric(out["SPORDER"]).astype("int16")
    REPW_DIR.mkdir(parents=True, exist_ok=True)
    out.to_parquet(f)
    return st, len(out)


def fetch_repw(workers: int = 6, states: list[int] | None = None) -> None:
    REPW_DIR.mkdir(parents=True, exist_ok=True)
    todo = [s for s in (states or ao.STATES) if not (REPW_DIR / f"{s:02d}.parquet").exists()]
    LOG.info("replicate weights: %d of %d states to fetch (%d workers)", len(todo), len(states or ao.STATES), workers)
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_fetch_repw_state, s): s for s in todo}
        for fu in as_completed(futs):
            s = futs[fu]
            try:
                _, n = fu.result()
            except Exception as e:                      # noqa: BLE001 -- report and keep going
                LOG.error("state %02d failed: %s", s, ao.Api().scrub(f"{type(e).__name__}: {e}")[:300])
                continue
            done += 1
            LOG.info("state %02d: %d persons  (%d/%d, %.0fs)", s, n, done, len(todo), time.time() - t0)
    missing = [s for s in (states or ao.STATES) if not (REPW_DIR / f"{s:02d}.parquet").exists()]
    if missing:
        raise SystemExit(f"replicate weights missing for states {missing}; rerun")


# =============================================================================
# 1. ACS PUMS 2009-13: the Table 2 sample, cell sums and statistics
# =============================================================================
#
# Everything below is a function of one weight vector `pw` (PWGTP or a replicate weight), so a
# standard error is the Census successive-difference formula  V = (4/80) sum_r (theta_r - theta)^2
# applied to the same code.  The person sample is `acs_occupations.build_sample(..., rules="workbook")`
# on the 2013 5-year PUMS (data years 2009-13), ages 25-34, at least grade 9 (SCHL >= 12), armed
# forces and the unemployed with 15+ usual hours dropped, i.e. the Table 2 / HHJK universe.
#
# Each person enters one or two (person, cell) rows with unweighted weight u: the reported
# occupation (u = w_mkt: 1, or 0.5 for 15-29 usual hours) and home production (u = w_home), so a
# person's u sum to 1.  A cell is (gender, group) with the 21 groups of `acs_occupations.ALL_GROUPS`
# (19 market occupations, home production, K-12 teachers).  Table 2's counts are the sums of u.
#
# Two ways to aggregate cell means into one number are reported throughout:
#   P   population: every row weighted by PWGTP x u (the plain weighted mean of the population in
#       question);
#   T2  the Table 2 convention of the existing targets (e.g. teacher_premium_m/f = +0.0557/+0.0860):
#       cell means weighted by PWGTP x u, combined across cells with the Table 2 counts (sums of
#       unweighted u), because the model's occupation shares match those counts.  For a mean over
#       a subset (BA holders, say) the weights are the counts of that subset in each cell.
# The two differ because PWGTP-weighted shares differ from unweighted ones (K-12 teachers are 1.94%
# of men unweighted and 1.76% PWGTP-weighted).

K12, HOME = ao.K12, ao.HOME
GROUPS = list(ao.ALL_GROUPS)                   # 19 market groups, home production, K-12 teachers
GI = {g: i for i, g in enumerate(GROUPS)}
NG = len(GROUPS)                               # 21
G_HOME, G_K12 = GI[HOME], GI[K12]
N_MKT = G_HOME                                 # market groups are 0..18
NC = 2 * NG                                    # cells: gender (0 men, 1 women) x group
SCHL_BA, SCHL_GRAD = 21, 22                    # bachelor's degree; master's, professional, doctorate
NAGE = ao.AGE_MAX - ao.AGE_MIN + 1             # ages 25..34
GENDERS = ("men", "women", "pooled")
K12_CODE_LABEL = {"2300": "pre-K and kindergarten teachers", "2310": "elementary and middle school teachers",
                  "2320": "secondary school teachers", "2330": "special education teachers",
                  "2340": "other teachers and instructors"}
CORE_CODES = ("2310", "2320", "2330")          # the K-12 teachers of licensure: elementary and middle, secondary, special education
# weeks-worked brackets of the 2009-13 PUMS (WKW) -> weeks used for an hourly wage; the Table 2 rule admits brackets 1 and 2 only
WKW_WEEKS = {1: 51.0, 2: 49.0}
WKW_LABEL = {1: "50-52 weeks", 2: "48-49 weeks", 3: "40-47 weeks", 4: "27-39 weeks", 5: "14-26 weeks", 6: "13 or fewer weeks"}
WEEKS_RULES = {"wks40": {3: 43.5}, "wks27": {3: 43.5, 4: 33.0}}   # sensitivity: also admit 40-47 (midpoint 43.5) / and 27-39 (33.0) weeks


def cell_mask(gender: str, pop: str) -> np.ndarray:
    """Boolean mask over the NC cells: gender in {men, women, pooled}; pop in {all, teacher, nonteacher,
    market_nonteacher}."""
    c = np.arange(NC)
    g, gi = c // NG, c % NG
    m = np.ones(NC, bool)
    if gender == "men":
        m &= g == 0
    elif gender == "women":
        m &= g == 1
    if pop == "teacher":
        m &= gi == G_K12
    elif pop == "nonteacher":
        m &= gi != G_K12
    elif pop == "market_nonteacher":
        m &= gi < N_MKT
    return m


def load_acs(need_repw: bool = True) -> "AcsData":
    """The Table 2 sample of the 2009-13 5-year PUMS, with the replicate weights attached."""
    persons = ao.load_person_file("acs5_2009_13")
    # SPORDER and the weeks-worked bracket WKW are not kept by load_person_file (which keeps only the 48+ week
    # brackets, as `weeks`): read them from the same files in the same order.
    files = sorted((ao.RAW / "acs5_2013").glob("*.parquet"))
    keys = pd.concat([pd.read_parquet(f, columns=["SERIALNO", "SPORDER", "WKW"]) for f in files], ignore_index=True)
    if len(keys) != len(persons) or not (keys.SERIALNO.values == persons.serialno.values).all():
        raise RuntimeError("person file and key columns do not line up")
    persons["sporder"] = pd.to_numeric(keys.SPORDER).astype("int16").values
    persons["wkw"] = pd.to_numeric(keys.WKW, errors="coerce").values
    cw = pd.read_csv(ao.EST / "acs_occupations_crosswalk.csv", dtype={"vintage": str, "occ": str})
    d = ao.build_sample(persons, cw, rules="workbook").reset_index(drop=True)
    rep = None
    if need_repw:
        fr = []
        for st in ao.STATES:
            f = REPW_DIR / f"{st:02d}.parquet"
            if not f.exists():
                raise SystemExit(f"{f} missing: run `python3 education_occupation.py --fetch-repw` first")
            fr.append(pd.read_parquet(f).assign(st=f"{st:02d}"))
        R = pd.concat(fr, ignore_index=True)
        m = d[["serialno", "sporder", "st", "pwgtp"]].merge(
            R, left_on=["serialno", "sporder", "st"], right_on=["SERIALNO", "SPORDER", "st"], how="left")
        if len(m) != len(d) or m.PWGTP.isna().any() or not (m.PWGTP.values == m.pwgtp.values).all():
            raise RuntimeError("replicate weights do not match the person file")
        rep = np.ascontiguousarray(m[[f"PWGTP{k}" for k in range(1, N_REP + 1)]].to_numpy(np.float64).T)
        del R, m
    return AcsData(d, rep)


class AcsData:
    """Person arrays, the (person, cell) stack and the wage sample of the Table 2 sample."""

    def __init__(self, d: pd.DataFrame, rep: np.ndarray | None) -> None:
        self.d, self.rep = d, rep
        n = self.n = len(d)
        self.pw = d.pwgtp.to_numpy(float)
        self.g = (d.gender.to_numpy() == "F").astype(np.int64)
        schl = d.schl.to_numpy(float)
        self.ba = (schl >= SCHL_BA).astype(float)
        self.grad = (schl >= SCHL_GRAD).astype(float)
        self.bao = (schl == SCHL_BA).astype(float)
        self.yrs = d.yrs.to_numpy(float)
        self.age_i = d.age.to_numpy(int) - ao.AGE_MIN
        wmk, whm = d.w_mkt.to_numpy(float), d.w_home.to_numpy(float)
        gi = d.group.map(GI).to_numpy(float)                      # NaN where there is no market occupation
        # ---- the (person, cell) stack
        im, ih = np.flatnonzero(wmk > 0), np.flatnonzero(whm > 0)
        if np.isnan(gi[im]).any():
            raise RuntimeError("a market row has no group")
        self.s_person = np.concatenate([im, ih])
        self.s_u = np.concatenate([wmk[im], whm[ih]])
        s_gi = np.concatenate([gi[im].astype(np.int64), np.full(len(ih), G_HOME)])
        self.s_cell = self.g[self.s_person] * NG + s_gi
        if abs(np.bincount(self.s_person, weights=self.s_u, minlength=n) - 1).max() > 1e-12:
            raise RuntimeError("a person's cell weights do not sum to one")
        sp = self.s_person
        self.svars = {"one": np.ones(len(sp)), "ba": self.ba[sp], "nba": 1 - self.ba[sp], "grad": self.grad[sp],
                      "yrs": self.yrs[sp], "yba": (self.yrs * self.ba)[sp], "ynba": (self.yrs * (1 - self.ba))[sp]}
        # Table 2 counts (sums of the unweighted u) and their BA / non-BA parts, by cell
        self.cnt = {k: np.bincount(self.s_cell, weights=self.s_u * v, minlength=NC) for k, v in self.svars.items()}
        # K-12 teacher rows by Census occupation code (small arrays, reused in every replicate)
        k12 = np.flatnonzero((self.s_cell % NG) == G_K12)
        self.k12_rows = k12
        self.k12_code = d.occ.to_numpy(object)[sp[k12]]
        self.k12_g = self.g[sp[k12]]
        # ---- the wage sample: employed 30+ hours, 48+ weeks, income floor, with a market group
        ws = d.wage_ok.to_numpy(bool) & d.group.notna().to_numpy()
        wi = self.wi = np.flatnonzero(ws)
        self.w_gi = gi[wi].astype(np.int64)
        self.w_g = self.g[wi]
        self.w_cell = self.w_g * NG + self.w_gi
        self.w_lw = np.log(d.wage99.to_numpy(float)[wi])
        self.w_age = self.age_i[wi]
        self.w_ba, self.w_grad, self.w_bao = self.ba[wi], self.grad[wi], self.bao[wi]
        self.w_teacher = self.w_gi == G_K12
        # K-12 teachers by Census code in the wage sample: the core codes (2310-2330), and the code index 0..4 for fixed effects
        w_occ = d.occ.to_numpy(object)[wi]
        self.w_core = np.isin(w_occ, CORE_CODES) & self.w_teacher
        self.w_code = np.where(self.w_teacher, pd.Series(w_occ).map({c: i for i, c in enumerate(K12_CODE_LABEL)}).fillna(-1).to_numpy(), -1).astype(np.int64)
        if self.w_teacher.any() and (self.w_code[self.w_teacher] < 0).any():
            raise RuntimeError("a K-12 teacher record has a Census code outside 2300-2340")
        # full-time workers with a market group, for the weeks-rule exclusion shares (needs the WKW bracket)
        self.has_weeks = "wkw" in d.columns
        if self.has_weeks:
            self.x_full = (d.full.to_numpy(bool) & d.group.notna().to_numpy())
            self.x_teacher = self.x_full & (d.group.to_numpy(object) == K12)
            self.x_wkw = d.wkw.to_numpy(float)
            self.x_wage_ok = d.wage_ok.to_numpy(bool)
        # (gender, market group) fixed effects of non-teachers: 0..37; teachers -1
        self.w_fe = np.where(self.w_gi < N_MKT, self.w_g * N_MKT + self.w_gi, -1)
        self.n_fe = 2 * N_MKT
        self.wage_count = np.bincount(self.w_cell, minlength=NC).astype(float)        # wage-sample records per cell
        self.wage_count_ba = np.bincount(self.w_cell, weights=self.w_ba, minlength=NC)
        self.kt = {"men": float(np.sum(self.w_teacher & (self.w_g == 0))),
                   "women": float(np.sum(self.w_teacher & (self.w_g == 1)))}

    def pweights(self, pw: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(P-weights of the stack rows, P-weights of the wage sample) for person weights `pw`."""
        return pw[self.s_person] * self.s_u, pw[self.wi]


def _ratio_stat(S: dict, C: dict, m: np.ndarray, num: str, den: str, conv: str) -> float:
    """Mean of `num` over `den` across the cells in `m` under convention P or T2 (see above)."""
    if conv == "P":
        d_ = S[den][m].sum()
        return float(S[num][m].sum() / d_) if d_ > 0 else np.nan
    ok = m & (S[den] > 0)
    w = C[den][ok]
    return float((w * S[num][ok] / S[den][ok]).sum() / w.sum()) if w.sum() > 0 else np.nan


def _wls_coef(X: np.ndarray, W: np.ndarray, L: np.ndarray, col: int = 1) -> float:
    """Coefficient `col` of the weighted regression of y on the cell design X (cells x p) from the cell
    sums W = sum w and L = sum w y:  beta = (X'WX)^-1 X'L.  Cells with no weight and columns with no
    variation among the occupied cells are dropped."""
    ok = W > 0
    X = X[ok]
    keep = np.abs(X).sum(axis=0) > 0
    new_col = int(np.cumsum(keep)[col] - 1)
    X = X[:, keep]
    A = X.T @ (W[ok, None] * X)
    b = X.T @ L[ok]
    try:
        return float(np.linalg.solve(A, b)[new_col])
    except np.linalg.LinAlgError:
        return float(np.linalg.lstsq(A, b, rcond=None)[0][new_col])


def age_flag_coef(age: np.ndarray, flag: np.ndarray, g: np.ndarray, lw: np.ndarray, w: np.ndarray,
                  pooled: bool, fe: np.ndarray | None = None, n_fe: int = 0) -> float:
    """Coefficient on `flag` in the weighted regression of log wage on flag, single-year age dummies,
    [female if `pooled`] and [fixed-effect dummies `fe` in 0..n_fe-1, which absorb gender], computed on the
    cells (age, flag, [female], [fe])."""
    nfe = n_fe if fe is not None else 1
    code = age + NAGE * (flag.astype(np.int64) + 2 * (g if (pooled and fe is None) else 0))
    if fe is not None:
        code = code + NAGE * 4 * fe
    K = NAGE * 4 * nfe
    W = np.bincount(code, weights=w, minlength=K)
    L = np.bincount(code, weights=w * lw, minlength=K)
    cc = np.arange(K)
    a_, f_, g_, k_ = cc % NAGE, (cc // NAGE) % 2, (cc // (NAGE * 2)) % 2, cc // (NAGE * 4)
    cols = [np.ones(K), f_.astype(float)] + [(a_ == a).astype(float) for a in range(1, NAGE)]
    if pooled and fe is None:
        cols.append(g_.astype(float))
    if fe is not None:
        cols += [(k_ == k).astype(float) for k in range(1, nfe)]
    return _wls_coef(np.column_stack(cols), W, L, 1)


def acs_stats(A: AcsData, pw: np.ndarray, wage_only: bool = False) -> dict[str, float]:
    """All ACS statistics for person weights `pw` (point estimate: PWGTP).  `wage_only` skips the composition,
    schooling, K-12 code and weeks-rule statistics, which do not depend on the wage sample (the weeks-rule
    sensitivity recomputes only the wage statistics)."""
    w_s, w_w = A.pweights(pw)
    S = {k: np.bincount(A.s_cell, weights=w_s * v, minlength=NC) for k, v in A.svars.items()}
    C = A.cnt
    out: dict[str, float] = {}
    if not wage_only:
        # ---- composition and schooling in the Table 2 universe (home production included)
        for gk in GENDERS:
            mall, mnt, mt, mm = (cell_mask(gk, p) for p in ("all", "nonteacher", "teacher", "market_nonteacher"))
            for conv in ("P", "T2"):
                for pop, m in (("all", mall), ("nonteacher", mnt)):
                    out[f"ba_share|{gk}|{pop}|{conv}"] = _ratio_stat(S, C, m, "ba", "one", conv)
                    out[f"years|{gk}|{pop}|{conv}"] = _ratio_stat(S, C, m, "yrs", "one", conv)
                    out[f"years_ba|{gk}|{pop}|{conv}"] = _ratio_stat(S, C, m, "yba", "ba", conv)
                    out[f"years_nba|{gk}|{pop}|{conv}"] = _ratio_stat(S, C, m, "ynba", "nba", conv)
                out[f"ba_share|{gk}|market_nonteacher|{conv}"] = _ratio_stat(S, C, mm, "ba", "one", conv)
            out[f"teacher_share|{gk}|P"] = float(S["one"][mt].sum() / S["one"][mall].sum())
            out[f"teacher_share|{gk}|T2"] = float(C["one"][mt].sum() / C["one"][mall].sum())
            out[f"ba_share|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "ba", "one", "P")
            out[f"grad_share|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "grad", "one", "P")
            out[f"grad_share_of_ba|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "grad", "ba", "P")
            out[f"years|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "yrs", "one", "P")
            out[f"years_ba|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "yba", "ba", "P")
            out[f"years_nba|{gk}|teacher|P"] = _ratio_stat(S, C, mt, "ynba", "nba", "P")
            # teachers' schooling over non-teachers' (the s_ratio of the calibration: 1.176 from the Mincer workbook)
            out[f"s_ratio|{gk}|P"] = out[f"years|{gk}|teacher|P"] / out[f"years|{gk}|nonteacher|P"]
            out[f"s_ratio|{gk}|T2"] = out[f"years|{gk}|teacher|P"] / out[f"years|{gk}|nonteacher|T2"]
        # ---- K-12 teachers by Census occupation code (population weights)
        wk, bk = w_s[A.k12_rows], A.svars["ba"][A.k12_rows]
        for gk in GENDERS:
            gsel = np.ones(len(wk), bool) if gk == "pooled" else (A.k12_g == (0 if gk == "men" else 1))
            tot = wk[gsel].sum()
            for code in list(K12_CODE_LABEL) + ["core"]:
                sel = gsel & (np.isin(A.k12_code, CORE_CODES) if code == "core" else (A.k12_code == code))
                den = wk[sel].sum()
                out[f"k12_code_share|{gk}|{code}"] = float(den / tot)
                out[f"k12_code_ba_share|{gk}|{code}"] = float((wk[sel] * bk[sel]).sum() / den) if den > 0 else np.nan
        # ---- who is outside the wage sample: full-time workers (employed, 30+ usual hours) by weeks-worked bracket
        if A.has_weeks:
            xp = {"teacher": A.x_teacher, "teacher_ba": A.x_teacher & (A.ba > 0), "nonteacher": A.x_full & ~A.x_teacher,
                  "nonteacher_ba": A.x_full & ~A.x_teacher & (A.ba > 0)}
            for gk in GENDERS:
                gsel = np.ones(A.n, bool) if gk == "pooled" else (A.g == (0 if gk == "men" else 1))
                for pop, m0 in xp.items():
                    m = m0 & gsel
                    den = pw[m].sum()
                    out[f"excl|{gk}|{pop}|in_wage_sample"] = float(pw[m & A.x_wage_ok].sum() / den)
                    for b in WKW_LABEL:
                        out[f"excl|{gk}|{pop}|wkw{b}"] = float(pw[m & (A.x_wkw == b)].sum() / den)
                    out[f"excl|{gk}|{pop}|wkw_missing"] = float(pw[m & np.isnan(A.x_wkw)].sum() / den)
                    out[f"excl|{gk}|{pop}|wk48_outside"] = float(pw[m & (A.x_wkw <= 2) & ~A.x_wage_ok].sum() / den)
    # ---- log hourly wages (wage sample)
    isT = A.w_teacher
    Wk = {"all": w_w, "ba": w_w * A.w_ba, "nba": w_w * (1 - A.w_ba)}
    Wc = {k: np.bincount(A.w_cell, weights=v, minlength=NC) for k, v in Wk.items()}
    Lc = {k: np.bincount(A.w_cell, weights=v * A.w_lw, minlength=NC) for k, v in Wk.items()}
    # BA+ minus non-BA among non-teacher market workers
    for gk in GENDERS:
        mm = cell_mask(gk, "market_nonteacher")
        pooled = gk == "pooled"
        gsel = np.ones(len(isT), bool) if pooled else (A.w_g == (0 if gk == "men" else 1))
        for conv in ("P", "T2"):
            for k in ("ba", "nba"):
                if conv == "P":
                    out[f"mean_lw_{k}|{gk}|nonteacher|P"] = float(Lc[k][mm].sum() / Wc[k][mm].sum())
                else:
                    ok = mm & (Wc[k] > 0)
                    out[f"mean_lw_{k}|{gk}|nonteacher|T2"] = float((C[k][ok] * Lc[k][ok] / Wc[k][ok]).sum() / C[k][ok].sum())
            out[f"ba_premium|{gk}|raw|{conv}"] = out[f"mean_lw_ba|{gk}|nonteacher|{conv}"] - out[f"mean_lw_nba|{gk}|nonteacher|{conv}"]
        sel = gsel & ~isT
        out[f"ba_premium|{gk}|age|P"] = age_flag_coef(A.w_age[sel], A.w_ba[sel], A.w_g[sel], A.w_lw[sel], w_w[sel], pooled)
        out[f"ba_premium|{gk}|age_occupation|P"] = age_flag_coef(
            A.w_age[sel], A.w_ba[sel], A.w_g[sel], A.w_lw[sel], w_w[sel], pooled, fe=A.w_fe[sel], n_fe=A.n_fe)
    # ---- K-12 teachers minus the comparison workers, by gender.  Teacher groups: "all" K-12 teachers (suffix "": the existing
    # targets), those with a BA or more ("_ba": every teacher has a BA in a model with a licensure rule) and the core codes
    # 2310-2330 ("_core"); comparison workers: all non-teacher market workers or those with a BA or more.
    tsub = {"": np.ones(len(isT), bool), "_ba": A.w_ba > 0, "_core": A.w_core}
    comps = (("all_nonteachers", "all"), ("ba_nonteachers", "ba"))
    comp_mean: dict[tuple[str, str, str], float] = {}
    for gk in ("men", "women"):
        mm = cell_mask(gk, "market_nonteacher")
        for comp, k in comps:
            ok = mm & (Wc[k] > 0)
            cnt_t2 = C["one"] if k == "all" else C["ba"]                 # Table 2 counts of the comparison workers
            cnt_ws = A.wage_count if k == "all" else A.wage_count_ba      # wage-sample records
            cnt_pw = S["one"] if k == "all" else S["ba"]                 # PWGTP-weighted counts of the same workers
            cm = Lc[k][ok] / Wc[k][ok]                                   # cell mean log wages
            for tag, cnt in (("T2", cnt_t2), ("wsw", cnt_ws), ("PW", cnt_pw)):
                comp_mean[gk, comp, tag] = float((cnt[ok] * cm).sum() / cnt[ok].sum())
            comp_mean[gk, comp, "P"] = float(Lc[k][mm].sum() / Wc[k][mm].sum())
    wts = {k: A.kt[k] / (A.kt["men"] + A.kt["women"]) for k in A.kt}
    out["teacher_wage_weight|men"], out["teacher_wage_weight|women"] = wts["men"], wts["women"]
    for sfx, tm_mask in tsub.items():
        for gk in ("men", "women"):
            gidx = 0 if gk == "men" else 1
            tsel = isT & tm_mask & (A.w_g == gidx)
            t_mean = float((w_w[tsel] * A.w_lw[tsel]).sum() / w_w[tsel].sum())
            out[f"mean_lw_teacher{sfx}|{gk}"] = t_mean
            for comp, k in comps:
                for tag in ("T2", "wsw", "PW", "P"):
                    out[f"teacher_premium{sfx}|{gk}|{comp}|{tag}"] = t_mean - comp_mean[gk, comp, tag]
                if sfx == "":
                    out[f"mean_lw_{comp}|{gk}|T2"] = comp_mean[gk, comp, "T2"]
                sel = (A.w_g == gidx) & ((isT & tm_mask) | (~isT & ((A.w_ba > 0) if comp == "ba_nonteachers" else True)))
                out[f"teacher_premium{sfx}|{gk}|{comp}|age"] = age_flag_coef(
                    A.w_age[sel], (isT & tm_mask)[sel].astype(float), A.w_g[sel], A.w_lw[sel], w_w[sel], False)
        tsel = isT & tm_mask
        out[f"mean_lw_teacher{sfx}|pooled"] = float((w_w[tsel] * A.w_lw[tsel]).sum() / w_w[tsel].sum())
        for comp, k in comps:
            for tag in ("T2", "wsw", "PW", "P", "age"):
                out[f"teacher_premium{sfx}|pooled|{comp}|{tag}"] = sum(
                    wts[g_] * out[f"teacher_premium{sfx}|{g_}|{comp}|{tag}"] for g_ in ("men", "women"))
            mm = cell_mask("pooled", "market_nonteacher")
            out[f"teacher_premium{sfx}|pooled|{comp}|P_nogenderadj"] = out[f"mean_lw_teacher{sfx}|pooled"] - float(Lc[k][mm].sum() / Wc[k][mm].sum())
    # ---- teachers with a graduate degree minus teachers with exactly a bachelor's degree
    for gk in GENDERS:
        pooled = gk == "pooled"
        sel = isT & (A.w_ba > 0) & (np.ones(len(isT), bool) if pooled else (A.w_g == (0 if gk == "men" else 1)))
        wg, wb = w_w[sel] * A.w_grad[sel], w_w[sel] * A.w_bao[sel]
        out[f"grad_teacher_premium|{gk}|raw|P"] = float((wg * A.w_lw[sel]).sum() / wg.sum() - (wb * A.w_lw[sel]).sum() / wb.sum())
        out[f"grad_teacher_premium|{gk}|age|P"] = age_flag_coef(
            A.w_age[sel], A.w_grad[sel], A.w_g[sel], A.w_lw[sel], w_w[sel], pooled)
    # ---- K-12 teachers with a BA or more minus those without (pay penalty of the non-BA teaching track): raw, age, age and Census code
    for gk in GENDERS:
        pooled = gk == "pooled"
        sel = isT & (np.ones(len(isT), bool) if pooled else (A.w_g == (0 if gk == "men" else 1)))
        wb_, wn_ = w_w[sel] * A.w_ba[sel], w_w[sel] * (1 - A.w_ba[sel])
        out[f"teacher_ba_gap|{gk}|raw|P"] = float((wb_ * A.w_lw[sel]).sum() / wb_.sum() - (wn_ * A.w_lw[sel]).sum() / wn_.sum())
        out[f"teacher_ba_gap|{gk}|age|P"] = age_flag_coef(A.w_age[sel], A.w_ba[sel], A.w_g[sel], A.w_lw[sel], w_w[sel], pooled)
        nk = len(K12_CODE_LABEL)
        fe = (A.w_g[sel] * nk + A.w_code[sel]) if pooled else A.w_code[sel]          # code x gender absorbs sex when pooled
        out[f"teacher_ba_gap|{gk}|age_code|P"] = age_flag_coef(
            A.w_age[sel], A.w_ba[sel], A.w_g[sel], A.w_lw[sel], w_w[sel], pooled, fe=fe, n_fe=2 * nk if pooled else nk)
    return out


def acs_counts(A: AcsData) -> dict[str, Any]:
    """Unweighted record counts behind the statistics."""
    out: dict[str, Any] = {"persons": {"men": int((A.g == 0).sum()), "women": int((A.g == 1).sum()), "pooled": int(A.n)}}
    for gk in GENDERS:
        gm = np.ones(A.n, bool) if gk == "pooled" else (A.g == (0 if gk == "men" else 1))
        # K-12 teachers as person-equivalents (u sums) and BA holders
        mt = cell_mask(gk, "teacher")
        out.setdefault("teachers_u", {})[gk] = float(A.cnt["one"][mt].sum())
        out.setdefault("nonteachers_u", {})[gk] = float(A.cnt["one"][cell_mask(gk, "nonteacher")].sum())
        out.setdefault("persons_ba", {})[gk] = int((gm & (A.ba > 0)).sum())
    wg = lambda gk: np.ones(len(A.wi), bool) if gk == "pooled" else (A.w_g == (0 if gk == "men" else 1))   # noqa: E731
    out["wage_sample"] = {gk: {"all": int(wg(gk).sum()), "teachers": int((wg(gk) & A.w_teacher).sum()),
                               "nonteachers": int((wg(gk) & ~A.w_teacher).sum()),
                               "nonteachers_ba": int((wg(gk) & ~A.w_teacher & (A.w_ba > 0)).sum()),
                               "teachers_ba": int((wg(gk) & A.w_teacher & (A.w_ba > 0)).sum()),
                               "teachers_nonba": int((wg(gk) & A.w_teacher & (A.w_ba == 0)).sum()),
                               "teachers_core": int((wg(gk) & A.w_core).sum()),
                               "teachers_grad": int((wg(gk) & A.w_teacher & (A.w_grad > 0)).sum()),
                               "teachers_ba_only": int((wg(gk) & A.w_teacher & (A.w_bao > 0)).sum())} for gk in GENDERS}
    if A.has_weeks:
        out["full_time_workers"] = {gk: {"teachers": int((A.x_teacher & (np.ones(A.n, bool) if gk == "pooled" else (A.g == (0 if gk == "men" else 1)))).sum()),
                                         "nonteachers": int(((A.x_full & ~A.x_teacher) & (np.ones(A.n, bool) if gk == "pooled" else (A.g == (0 if gk == "men" else 1)))).sum())}
                                    for gk in GENDERS}
    return out


def alt_wage_sample(d: pd.DataFrame, extra_weeks: dict[int, float]) -> pd.DataFrame:
    """Copy of the person frame whose hourly-wage sample also admits the weeks-worked brackets in `extra_weeks`
    ({WKW bracket: weeks used}), under the workbook rules of `acs_occupations.build_sample` (income floor on ADJINC income,
    wage = reference-period income deflated by the survey-year CPI-U to 1999 dollars over usual hours x weeks).  The
    recomputation is checked against the original `wage99` / `wage_ok` for the 48+ week brackets, where they must agree."""
    wk = d.wkw.map({**WKW_WEEKS, **extra_weeks})
    raw = d.wagp.fillna(0) + d.semp.fillna(0)
    cpi_y = d.data_year.map(ao.CPI_U)
    with np.errstate(divide="ignore", invalid="ignore"):
        wage99 = raw / cpi_y * ao.CPI_U[1999] / (d.wkhp * wk)
    ok = d.full & wk.notna() & (d.inc2010 >= ao.INCOME_FLOOR_2010) & (d.wkhp > 0) & (wage99 > 0)
    base = d.wkw.isin(list(WKW_WEEKS))
    if not (d.wage_ok[base].to_numpy() == ok[base].to_numpy()).all() or not np.allclose(
            d.wage99[base & d.wage_ok].to_numpy(float), wage99[base & ok].to_numpy(float), rtol=1e-12, atol=0):
        raise RuntimeError("the recomputed 48+ week wage sample differs from acs_occupations.build_sample")
    out = d.copy()
    out["wage99"], out["wage_ok"] = wage99, ok
    return out


def replicate_se(A: AcsData, point: dict[str, float], **kw: Any) -> dict[str, float]:
    """Census successive-difference replicate-weight standard errors, V = (4/80) sum_r (theta_r - theta)^2.
    Keyword arguments go to `acs_stats` (e.g. wage_only=True)."""
    if A.rep is None:
        return {k: np.nan for k in point}
    keys = list(point)
    acc = np.zeros(len(keys))
    t0 = time.time()
    for r in range(N_REP):
        st = acs_stats(A, A.rep[r], **kw)
        v = np.array([st[k] for k in keys])
        acc += (v - np.array([point[k] for k in keys])) ** 2
        if (r + 1) % 20 == 0:
            LOG.info("  replicate %d/%d (%.0fs)", r + 1, N_REP, time.time() - t0)
    se = np.sqrt(4.0 / N_REP * acc)
    return dict(zip(keys, se.tolist()))


# =============================================================================
# 2. NLSY79: latent AFQT of BA holders, P(BA | AFQT tercile), the teacher gaps (confirmation)
# =============================================================================
#
# Machinery of `nlsy_teachers.py` (and so of `nlsy_ability.py`): block-A population sample, ability =
# `afqt_c` (standardized, age-normed composite of the four ASVAB IRT z-scores, 1979 weights), a gap in
# means of `afqt_c` divided by sqrt(R_comp_pop) is a gap in SDs of the latent factor, one observation
# per person, household-cluster bootstrap with the draws of `nlsy_teachers.run_bootstrap` (B = 500,
# seed 20260929), AFQT standardization and R_comp_pop re-estimated in every draw.  BA+ = highest
# grade completed >= 16 at some round at ages 25-34 (`nlsy_teachers.COLLEGE_HGC`; the extract has no
# degree variable).

NLSY_SEED = 20260929
NLSY_B = 500
COLLEGE_HGC = 16
BA_SAMPLES = {
    "all_hgc9": "all persons observed at 25-34, highest grade >= 9 (the ACS Table 2 universe): primary",
    "all": "all persons observed at 25-34, no grade rule",
    "emp_hgc9": "persons with an employed round at 25-34, highest grade >= 9",
    "emp": "persons with an employed round at 25-34 (the 'employed' sample of nlsy_teachers.py)",
    "emp_nonteacher_hgc9": "employed non-K-12-teachers, highest grade >= 9",
    "emp_nonteacher": "employed non-K-12-teachers (comparison group A of nlsy_teachers.py)",
}
TEACHER_SPEC = "c80-emp-all-maj"      # 1980 Census codes, all K-12 codes, majority of rounds, employed


def wlogit(x: np.ndarray, y: np.ndarray, w: np.ndarray, iters: int = 30) -> np.ndarray:
    """Weighted logistic regression of y on (1, x) by Newton steps; returns (a, b)."""
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


def share_adjusted_gap(x: np.ndarray, w: np.ndarray, a: float, b: float, target: float) -> tuple[float, float, float]:
    """Gap E[x | BA] - E[x | no BA] (in the units of x) if P(BA | x) = logit(a + delta + b x) with delta set so that
    the BA share is `target`: the logistic gradient is kept and only its level is shifted.  Returns (gap, delta, share)."""
    lo, hi = -30.0, 30.0
    sig = lambda e: 1.0 / (1.0 + np.exp(-e))                                             # noqa: E731
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (w * sig(a + mid + b * x)).sum() / w.sum() < target:
            lo = mid
        else:
            hi = mid
    delta = 0.5 * (lo + hi)
    p = sig(a + delta + b * x)
    m1 = (w * p * x).sum() / (w * p).sum()
    m0 = (w * (1 - p) * x).sum() / (w * (1 - p)).sum()
    return float(m1 - m0), float(delta), float((w * p).sum() / w.sum())


def gh_nodes(n: int = 41) -> tuple[np.ndarray, np.ndarray]:
    """Nodes and weights of Gauss-Hermite quadrature for E f(Z), Z ~ N(0, 1) (the weights sum to one)."""
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return x, w / np.sqrt(2.0 * np.pi)


def eiv_logit(x: np.ndarray, y: np.ndarray, w: np.ndarray, R: float, n_nodes: int = 41,
              iters: int = 100, tol: float = 1e-10) -> tuple[float, float]:
    """Maximum likelihood (a, b) of P(y = 1 | theta) = logistic(a + b theta) for a standard-normal latent theta that is observed
    only through the standardized score x = sqrt(R) theta + sqrt(1 - R) u (u standard normal, R the reliability), i.e. the
    errors-in-variables logit: P(y = 1 | x) = E logistic(a + b theta) with theta | x ~ N(sqrt(R) x, 1 - R), integrated by
    Gauss-Hermite quadrature; weights `w`; Fisher scoring with step halving, started at the naive logit with b / sqrt(R)."""
    t, om = gh_nodes(n_nodes)
    TH = np.sqrt(R) * x[:, None] + np.sqrt(1.0 - R) * t[None, :]
    wn = w / w.sum()

    def parts(beta: np.ndarray):
        L = 1.0 / (1.0 + np.exp(-(beta[0] + beta[1] * TH)))
        p = np.clip(L @ om, 1e-12, 1.0 - 1e-12)
        dL = L * (1.0 - L)
        ga, gb = dL @ om, (dL * TH) @ om                      # dp/da, dp/db
        ll = float((wn * (y * np.log(p) + (1.0 - y) * np.log(1.0 - p))).sum())
        r = wn * (y - p) / (p * (1.0 - p))
        v = wn / (p * (1.0 - p))
        score = np.array([(r * ga).sum(), (r * gb).sum()])
        info = np.array([[(v * ga * ga).sum(), (v * ga * gb).sum()], [(v * ga * gb).sum(), (v * gb * gb).sum()]])
        return ll, score, info

    a0, b0 = wlogit(x, y, w)
    beta = np.array([a0, b0 / np.sqrt(R)])
    ll, score, info = parts(beta)
    for _ in range(iters):
        step = np.linalg.solve(info, score)
        lam = 1.0
        while lam > 1e-8:
            new = beta + lam * step
            ll2, s2, i2 = parts(new)
            if ll2 >= ll - 1e-14:
                break
            lam *= 0.5
        else:
            break
        beta, ll, score, info = new, ll2, s2, i2
        if np.abs(lam * step).max() < tol:
            break
    return float(beta[0]), float(beta[1])


_TH_GRID = np.linspace(-8.0, 8.0, 3201)
_TH_PDF = np.exp(-0.5 * _TH_GRID ** 2)
_TH_PDF /= _TH_PDF.sum()
_TERC_CUT = NormalDist().inv_cdf(2.0 / 3.0)               # terciles of a standard normal: -0.4307, +0.4307


def latent_logit_stats(a: float, b: float, target: float | None = None) -> dict[str, float]:
    """For P(BA | theta) = logistic(a + shift + b theta) with theta ~ N(0, 1): the BA share, P(BA | tercile of theta) (q = 1, 2, 3)
    and the gap E[theta | BA] - E[theta | no BA].  With `target`, the intercept is shifted first so that the share equals it
    (the logistic gradient is kept, only the level changes).  Returns also the shift."""
    sig = lambda e: 1.0 / (1.0 + np.exp(-e))                                             # noqa: E731
    shift = 0.0
    if target is not None:
        lo, hi = -30.0, 30.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if (_TH_PDF * sig(a + mid + b * _TH_GRID)).sum() < target:
                lo = mid
            else:
                hi = mid
        shift = 0.5 * (lo + hi)
    p = sig(a + shift + b * _TH_GRID)
    out = {"share": float((_TH_PDF * p).sum()), "shift": float(shift),
           "gap": float((_TH_PDF * p * _TH_GRID).sum() / (_TH_PDF * p).sum() - (_TH_PDF * (1 - p) * _TH_GRID).sum() / (_TH_PDF * (1 - p)).sum())}
    for q, sel in enumerate((_TH_GRID <= -_TERC_CUT, (_TH_GRID > -_TERC_CUT) & (_TH_GRID <= _TERC_CUT), _TH_GRID > _TERC_CUT)):
        out[f"terc{q + 1}"] = float((_TH_PDF[sel] * p[sel]).sum() / _TH_PDF[sel].sum())
    return out


def evaluate_ba(fr, nt, ne, na, mult_hh: np.ndarray | None, share_target: dict[str, float]) -> dict[str, float]:
    """BA statistics for household multiplicities `mult_hh` (None = the sample itself)."""
    D = fr.D
    mult_p = np.ones(D.n) if mult_hh is None else mult_hh[D.hh_idx].astype(float)
    w = D.w * mult_p
    afqt, Rc, *_ = na.pop_afqt(D.z, D.asamp, w)
    sR = float(np.sqrt(Rc))
    hg = fr.hg_max
    okg = np.isfinite(hg)
    ba = okg & (hg >= COLLEGE_HGC)
    g9 = okg & (hg >= na.MIN_HGC)
    sp = next(s for s in fr.specs if s.name == TEACHER_SPEC)
    base = {"all": fr.observed & D.asamp & okg, "emp": sp.member & okg, "emp_nonteacher": sp.nonteacher & okg}
    samples = {k: v for k, v in base.items()}
    for k in ("all", "emp", "emp_nonteacher"):
        samples[f"{k}_hgc9"] = base[k] & g9
    out: dict[str, float] = {"R_comp_pop": float(Rc)}
    for sk, sm in samples.items():
        for gk, gm in fr.masks.items():
            m = sm & gm & np.isfinite(afqt) & (w > 0)
            m1, m0 = m & ba, m & ~ba
            k = f"{sk}|{gk}"
            y1, y0 = nt.wm(afqt, w, m1), nt.wm(afqt, w, m0)
            out[f"{k}|ba_share"] = nt.ratio(w[m1].sum(), w[m].sum())
            out[f"{k}|mean_ba"] = y1 / sR
            out[f"{k}|mean_nba"] = y0 / sR
            out[f"{k}|gap"] = (y1 - y0) / sR
            out[f"{k}|gap_obs"] = y1 - y0
    # ---- P(BA | tercile of AFQT): terciles cut over the (weighted) sample itself, pooled over gender
    tdict: dict[str, np.ndarray] = {}
    for sk in ("all_hgc9", "all"):
        m = samples[sk] & np.isfinite(afqt) & (w > 0)
        t, cuts = ne.tercile_of(np.where(m, afqt, np.nan), w * m)
        tdict[sk] = t
        out[f"terc|{sk}|cut1"], out[f"terc|{sk}|cut2"] = float(cuts[0]), float(cuts[1])
        for gk, gm in fr.masks.items():
            mg = m & gm
            sh = []
            for q in range(3):
                mq = mg & (t == q)
                s_ = nt.ratio(w[mq & ba].sum(), w[mq].sum())
                sh.append(s_)
                out[f"terc|{sk}|{gk}|ba_share|{q + 1}"] = s_
                out[f"terc|{sk}|{gk}|pop_share|{q + 1}"] = nt.ratio(w[mq].sum(), w[mg].sum())
            out[f"terc|{sk}|{gk}|ba_share|diff"] = sh[2] - sh[0]
            out[f"terc|{sk}|{gk}|ba_share|diff_21"] = sh[1] - sh[0]
            out[f"terc|{sk}|{gk}|ba_share|diff_32"] = sh[2] - sh[1]
    # ---- logit gradient and the gap at the ACS BA share (a level shift of the logit)
    for sk in ("all_hgc9",):
        for gk, gm in fr.masks.items():
            m = samples[sk] & gm & np.isfinite(afqt) & (w > 0)
            a_, b_ = wlogit(afqt[m], ba[m].astype(float), w[m])
            out[f"logit|{sk}|{gk}|slope_obs"] = float(b_)
            out[f"logit|{sk}|{gk}|slope_latent_approx"] = float(b_ / sR)   # E[latent | score] = sqrt(R) x score: the slope on the latent factor is b / sqrt(R) to first order
            out[f"logit|{sk}|{gk}|intercept"] = float(a_)
            gap_obs, delta, sh_ = share_adjusted_gap(afqt[m], w[m], a_, b_, share_target[gk])
            out[f"logit|{sk}|{gk}|gap_at_acs_share"] = gap_obs / sR
            out[f"logit|{sk}|{gk}|shift"] = delta
            # (1) observed-scale terciles: the fitted logit, and the same logit with its level shifted to the ACS BA share
            xs, ws, tm = afqt[m], w[m], tdict[sk][m]
            pf, pa = 1.0 / (1.0 + np.exp(-(a_ + b_ * xs))), 1.0 / (1.0 + np.exp(-(a_ + delta + b_ * xs)))
            for q in range(3):
                mq = tm == q
                out[f"terc_fit|{sk}|{gk}|ba_share|{q + 1}"] = float((ws[mq] * pf[mq]).sum() / ws[mq].sum())
                out[f"terc_acs|{sk}|{gk}|ba_share|{q + 1}"] = float((ws[mq] * pa[mq]).sum() / ws[mq].sum())
            out[f"terc_acs|{sk}|{gk}|ba_share|diff"] = out[f"terc_acs|{sk}|{gk}|ba_share|3"] - out[f"terc_acs|{sk}|{gk}|ba_share|1"]
            # (2) latent scale: errors-in-variables logit (the score has reliability R_comp_pop), its latent terciles and gap,
            #     at the NLSY BA share and at the ACS BA share
            aL, bL = eiv_logit(xs, ba[m].astype(float), ws, Rc)
            lt, la = latent_logit_stats(aL, bL), latent_logit_stats(aL, bL, share_target[gk])
            out[f"logit|{sk}|{gk}|slope_latent_eiv"] = bL
            out[f"logit|{sk}|{gk}|intercept_latent_eiv"] = aL
            out[f"logit|{sk}|{gk}|gap_latent_eiv"] = lt["gap"]
            out[f"logit|{sk}|{gk}|gap_latent_eiv_at_acs_share"] = la["gap"]
            for q in range(3):
                out[f"terc_latent|{sk}|{gk}|ba_share|{q + 1}"] = lt[f"terc{q + 1}"]
                out[f"terc_latent_acs|{sk}|{gk}|ba_share|{q + 1}"] = la[f"terc{q + 1}"]
            out[f"terc_latent|{sk}|{gk}|ba_share|diff"] = lt["terc3"] - lt["terc1"]
            out[f"terc_latent_acs|{sk}|{gk}|ba_share|diff"] = la["terc3"] - la["terc1"]
    # ---- the existing teacher gaps, recomputed (confirmation of nlsy_teachers.json)
    ev = nt.evaluate(fr, mult_hh)
    for gk in ("pooled", "men", "women"):
        for cl in ("nt", "col", "cc"):
            for stt in ("gap", "mean_t", "mean_c"):
                out[f"teacher|{gk}|{cl}|{stt}"] = ev[f"{TEACHER_SPEC}|{gk}|{cl}|{stt}"]
        out[f"teacher|{gk}|col_share_t"] = ev[f"{TEACHER_SPEC}|{gk}|col_share_t"]
        out[f"teacher|{gk}|col_share_nt"] = ev[f"{TEACHER_SPEC}|{gk}|col_share_nt"]
    return out


def run_nlsy(B: int, seed: int, share_target: dict[str, float]) -> dict[str, Any]:
    """Point estimates and household-bootstrap SEs of the NLSY79 statistics."""
    import nlsy_ability as na
    import nlsy_extensions as ne
    import nlsy_teachers as nt
    c = na.Cohort("nlsy79", HERE)
    fr = nt.Frame(c)
    point = evaluate_ba(fr, nt, ne, na, None, share_target)
    keys = list(point)
    rng = np.random.default_rng(seed)
    nhh = fr.D.nhh
    draws = np.full((B, len(keys)), np.nan)
    t0 = time.time()
    for b in range(B):
        cm = np.bincount(rng.integers(0, nhh, nhh), minlength=nhh).astype(float)    # the draws of nlsy_teachers.run_bootstrap
        try:
            d = evaluate_ba(fr, nt, ne, na, cm, share_target)
            draws[b] = [d.get(k, np.nan) for k in keys]
        except Exception as exc:                                                    # noqa: BLE001 -- a failed draw is a NaN row
            LOG.warning("  draw %d failed: %s", b, exc)
        if (b + 1) % 50 == 0 or b + 1 == B:
            LOG.info("  NLSY bootstrap %d/%d (%.0fs)", b + 1, B, time.time() - t0)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        se = np.nanstd(draws, axis=0, ddof=1) if B > 1 else np.full(len(keys), np.nan)
        lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0) if B > 1 else (np.full(len(keys), np.nan),) * 2
    stats = {k: {"value": float(point[k]), "se": float(se[i]), "ci95_lo": float(lo[i]), "ci95_hi": float(hi[i]),
                 "n_draws": int(np.isfinite(draws[:, i]).sum())} for i, k in enumerate(keys)}
    # unweighted counts behind the point estimates
    D = fr.D
    hg = fr.hg_max
    okg = np.isfinite(hg)
    ba = okg & (hg >= COLLEGE_HGC)
    sp = next(s for s in fr.specs if s.name == TEACHER_SPEC)
    base = {"all": fr.observed & D.asamp & okg, "emp": sp.member & okg, "emp_nonteacher": sp.nonteacher & okg}
    samples = dict(base)
    for k in list(base):
        samples[f"{k}_hgc9"] = base[k] & (hg >= na.MIN_HGC)
    counts = {sk: {gk: {"persons": int((sm & gm).sum()), "ba": int((sm & gm & ba).sum()),
                        "non_ba": int((sm & gm & ~ba).sum())} for gk, gm in fr.masks.items()}
              for sk, sm in samples.items()}
    counts["k12_teachers"] = {gk: int((sp.teacher & D.asamp & (D.w > 0) & gm).sum()) for gk, gm in fr.masks.items()}
    meta = {"B": B, "seed": seed, "failed_draws": int(np.sum(~np.isfinite(draws).any(axis=1))),
            "n_block_a": int(D.asamp.sum()), "n_observed_25_34": int((fr.observed & D.asamp).sum()),
            "n_missing_highest_grade": int((fr.observed & D.asamp & ~okg).sum()),
            "age_window": [na.AGE_LO, na.AGE_HI], "years": sorted(set(fr.py.Y.tolist())),
            "share_target": share_target}
    return {"stats": stats, "counts": counts, "meta": meta}


# =============================================================================
# 3. NCES Digest tables: teacher salary schedules (211.20), tuition (330.10), net price (331.30)
# =============================================================================
#
# The tables are read from the NCES HTML pages themselves (cached under raw/education_occupation/),
# not from a summary of them.

EDU_RAW = HERE / "raw" / "education_occupation"
NCES = "https://nces.ed.gov/programs/digest/"
NCES_TABLES = {
    "211.20": ("d22/tables/dt22_211.20.asp", "Digest of Education Statistics 2022, Table 211.20"),
    "211.20_d20": ("d20/tables/dt20_211.20.asp", "Digest of Education Statistics 2020, Table 211.20"),
    "330.10": ("d18/tables/dt18_330.10.asp", "Digest of Education Statistics 2018, Table 330.10"),
    "331.30": ("d19/tables/dt19_331.30.asp", "Digest of Education Statistics 2019, Table 331.30"),
    "303.70": ("d18/tables/dt18_303.70.asp", "Digest of Education Statistics 2018, Table 303.70"),
}


def nces_html(key: str) -> str:
    """The NCES page of table `key`, downloaded once into EDU_RAW."""
    path, _ = NCES_TABLES[key]
    f = EDU_RAW / (Path(path).name.replace(".asp", ".html"))
    if not f.exists():
        import requests
        EDU_RAW.mkdir(parents=True, exist_ok=True)
        r = requests.get(NCES + path, timeout=120, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        f.write_bytes(r.content)
    return f.read_text(encoding="iso-8859-1")


def table_rows(html: str) -> list[list[str]]:
    """Every <tr> of an NCES page as a list of cell texts (standard library parser)."""
    import re
    from html.parser import HTMLParser

    class P(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.rows: list[list[str]] = []
            self.cur: list[str] | None = None
            self.cell: list[str] | None = None

        def handle_starttag(self, tag, attrs):
            if tag == "tr":
                self.cur = []
            elif tag in ("td", "th"):
                self.cell = []

        def handle_endtag(self, tag):
            if tag in ("td", "th") and self.cell is not None and self.cur is not None:
                self.cur.append(re.sub(r"\s+", " ", "".join(self.cell)).strip())
                self.cell = None
            elif tag == "tr" and self.cur is not None:
                self.rows.append(self.cur)
                self.cur = None

        def handle_data(self, data):
            if self.cell is not None:
                self.cell.append(data)

    p = P()
    p.feed(html)
    return p.rows


def num(s: str) -> float:
    """'$12,615' -> 12615.0; '(365)' -> 365.0; '‡', '†', '—', '' -> NaN."""
    t = s.replace("$", "").replace(",", "").replace("(", "").replace(")", "").strip()
    try:
        return float(t)
    except ValueError:
        return float("nan")


def table_notes(html: str, start: str = "NOTE:") -> str:
    import html as H
    import re
    t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", html, flags=re.S)
    t = H.unescape(re.sub(r"<[^>]+>", "\n", t))
    t = re.sub(r"\n\s*\n+", "\n", re.sub(r"[ \t\r\f]+", " ", t))
    i = t.find(start)
    return re.split(r"\n20\d\d Tables and Figures", t[i:i + 3000])[0].strip() if i >= 0 else ""


# ---- 211.20: average base salary of full-time public school teachers by degree and experience
SAL_DEGREES = ("all", "ba", "ma", "eds", "phd")
SAL_ROWS = ["1 year or less", "2 years", "3 years", "4 years", "5 years", "6 to 9 years", "10 to 14 years",
            "15 to 19 years", "20 to 24 years", "25 to 29 years", "30 to 34 years", "35 years or more", "Total"]
LOW_EXPERIENCE = {"1-3 years": ["1 year or less", "2 years", "3 years"],
                  "1-5 years": ["1 year or less", "2 years", "3 years", "4 years", "5 years"],
                  "1-9 years": ["1 year or less", "2 years", "3 years", "4 years", "5 years", "6 to 9 years"]}


def parse_211_20(html: str) -> dict[str, dict[str, dict[str, Any]]]:
    """{school year: {experience row: {'n': (v, se), 'cur': {degree: (v, se)}, 'const': {...}}}}.
    Each statistic occupies three cells (value, blank, '(se)'); 11 statistics follow the row label: the number of
    teachers, base salary in current dollars for all teachers / bachelor's / master's / education specialist /
    doctor's, and the same five in constant dollars."""
    import re
    out: dict[str, dict[str, dict[str, Any]]] = {}
    year = None
    for r in table_rows(html):
        if len(r) != 34:
            continue
        if re.fullmatch(r"\d{4}-\d{2,4}", r[0]) and all(x == "" for x in r[1:]):
            year = r[0]
            out[year] = {}
            continue
        if year is None or r[0] not in SAL_ROWS:
            continue
        v = [num(x) for x in r[1::3]]
        se = [num(x) for x in r[3::3]]
        pairs = list(zip(v, se))
        out[year][r[0]] = {"n": pairs[0], "cur": dict(zip(SAL_DEGREES, pairs[1:6])),
                           "const": dict(zip(SAL_DEGREES, pairs[6:11]))}
    return out


def salary_premia(tab: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Master's over bachelor's base salary by experience row: log ratio, ratio, dollar gap, SE (delta method,
    the two means treated as independent), and the averages over the low-experience rows weighted by the number
    of teachers in the row (all degrees)."""
    rows: dict[str, Any] = {}
    for lab in SAL_ROWS:
        if lab not in tab:
            continue
        (b, sb), (m, sm) = tab[lab]["cur"]["ba"], tab[lab]["cur"]["ma"]
        if not (np.isfinite(b) and np.isfinite(m)):
            continue
        rows[lab] = {"ba": b, "ma": m, "se_ba": sb, "se_ma": sm, "n_teachers": tab[lab]["n"][0],
                     "log": float(np.log(m / b)), "se_log": float(np.hypot(sm / m, sb / b)),
                     "ratio_minus_1": float(m / b - 1), "dollars": float(m - b),
                     "se_dollars": float(np.hypot(sm, sb))}
    agg = {}
    for name, labs in LOW_EXPERIENCE.items():
        labs = [x for x in labs if x in rows]
        if not labs:
            continue
        n = np.array([rows[x]["n_teachers"] for x in labs])
        w = n / n.sum()
        agg[name] = {"rows": labs,
                     "log": float(sum(w_ * rows[x]["log"] for w_, x in zip(w, labs))),
                     "se_log": float(np.sqrt(sum((w_ * rows[x]["se_log"]) ** 2 for w_, x in zip(w, labs)))),
                     "dollars": float(sum(w_ * rows[x]["dollars"] for w_, x in zip(w, labs))),
                     "ba": float(sum(w_ * rows[x]["ba"] for w_, x in zip(w, labs))),
                     "ma": float(sum(w_ * rows[x]["ma"] for w_, x in zip(w, labs)))}
    return {"rows": rows, "aggregates": agg}


def nces_salary() -> dict[str, Any]:
    html = nces_html("211.20")
    tab = parse_211_20(html)
    out: dict[str, Any] = {"source": NCES + NCES_TABLES["211.20"][0], "title": NCES_TABLES["211.20"][1],
                           "notes": table_notes(html), "years": {}}
    for y in ("2011-12", "2015-16", "2017-18", "2020-21"):
        if y in tab:
            out["years"][y] = salary_premia(tab[y])
            out["years"][y]["total_row"] = {"ba": tab[y]["Total"]["cur"]["ba"][0], "ma": tab[y]["Total"]["cur"]["ma"][0],
                                            "n_teachers": tab[y]["Total"]["n"][0]}
    # the 2017-18 rows of the 2020 edition must equal those of the 2022 edition
    try:
        t20 = parse_211_20(nces_html("211.20_d20"))
        diffs = []
        for lab in SAL_ROWS:
            a, b = tab["2017-18"][lab], t20["2017-18"][lab]
            for dg in SAL_DEGREES:
                for i in (0, 1):
                    if not (np.isnan(a["cur"][dg][i]) and np.isnan(b["cur"][dg][i])):
                        diffs.append(abs(a["cur"][dg][i] - b["cur"][dg][i]))
        out["edition_check_2017_18_vs_digest_2020"] = {"n_compared": len(diffs), "max_abs_diff": float(max(diffs))}
    except Exception as exc:                                              # noqa: BLE001
        out["edition_check_2017_18_vs_digest_2020"] = {"error": str(exc)}
    return out


# ---- 330.10 / 331.30 / 303.70: tuition, net price, enrollment by control
def parse_330_10(html: str) -> dict[str, dict[str, dict[str, float]]]:
    """{control: {year: {measure_level: dollars}}}, current dollars, from the 25-cell rows: label, constant-dollar
    blocks (total; tuition and required fees; dormitory; board; each for all / 4-year / 2-year institutions), then the
    same four blocks in current dollars."""
    import re
    out: dict[str, dict[str, dict[str, float]]] = {}
    ctl = None
    for r in table_rows(html):
        if len(r) != 25:
            continue
        if all(x == "" for x in r[1:]) and r[0]:
            t = r[0].replace("- ", "")
            ctl = ("all" if t.startswith("All") else "public" if t.startswith("Public")
                   else "private" if t.startswith("Private") else "nonprofit" if t.startswith("Nonprofit")
                   else "forprofit" if t.startswith("For") else None)
            continue
        if ctl and re.fullmatch(r"\d{4}-\d{2}", r[0]):
            v = [num(x) for x in r[13:25]]
            d = {}
            for blk, name in enumerate(("total", "tuition_fees", "dorm", "board")):
                for j, lvl in enumerate(("all", "4yr", "2yr")):
                    d[f"{name}_{lvl}"] = v[3 * blk + j]
            out.setdefault(ctl, {})[r[0]] = d
    return out


def parse_331_30(html: str) -> dict[str, Any]:
    """Current-dollar rows of Table 331.30 for 2017-18 (the last four of sixteen columns: all institutions, public,
    private nonprofit, private for-profit): {level: {measure: {income level: [4 values]}}}."""
    out: dict[str, Any] = {}
    cur = False
    lvl = meas = None
    for r in table_rows(html):
        if len(r) == 2 and r[0] == "" and "dollars" in r[1]:
            cur = r[1].startswith("Current")
            continue
        if not cur:
            continue
        if r[0].startswith(("4-year", "2-year")) and all(x == "" for x in r[1:]):
            lvl = "4yr" if r[0].startswith("4") else "2yr"
            continue
        if r[0].startswith(("Grant and scholarship", "Net price")) and all(x == "" for x in r[1:]):
            meas = "grant" if r[0].startswith("Grant") else "net_price"
            continue
        if lvl and meas and len(r) == 17 and (r[0].startswith("All income") or r[0].startswith("$")):
            out.setdefault(lvl, {}).setdefault(meas, {})[r[0]] = [num(x) for x in r[13:17]]
    return out


def parse_303_70(html: str) -> dict[str, float]:
    """Fall 2017 undergraduate enrollment at 4-year institutions: total, public, private nonprofit, private for-profit."""
    sec = None
    for r in table_rows(html):
        if len(r) == 14 and all(x == "" for x in r[1:]) and r[0]:
            sec = r[0]
        if sec == "4-year institutions" and len(r) == 14 and r[0] == "2017":
            v = [num(x) for x in r[1:]]
            return {"total": v[0], "public": v[9], "private": v[10], "nonprofit": v[11], "forprofit": v[12]}
    raise RuntimeError("Table 303.70: 4-year 2017 row not found")


# College Board, Trends in College Pricing 2018 (Ma, Pender and Welch), Figures 9 and 10 (printed pages 18 and 19):
# full-time undergraduates, 2017-18, dollars of 2018, enrollment-weighted; net = published price minus grant aid
# from all sources and federal education tax benefits.  https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf
COLLEGE_BOARD_2018 = {
    "url": "https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf",
    "citation": "Ma, Pender and Welch (2018), Trends in College Pricing 2018, College Board; Figure 9 (public four-year, "
                "in-state, p. 18) and Figure 10 (private nonprofit four-year, p. 19), 2017-18, dollars of 2018",
    "public_4yr_in_state": {"published_tuition_fees": 10270, "published_tfrb": 21400, "net_tuition_fees": 3780,
                            "net_tfrb": 14910, "grant_aid_tax_benefits": 6490},
    "private_nonprofit_4yr": {"published_tuition_fees": 35720, "published_tfrb": 48380, "net_tuition_fees": 14500,
                              "net_tfrb": 27160, "grant_aid_tax_benefits": 21220},
}

def verify_college_board() -> dict[str, Any]:
    """Re-read the College Board figures from the cached report (`pdftotext -layout`) and compare them with the constants:
    in each of the two figure blocks the five price rows have nine columns (1998-99 ... 2018-19); 2017-18 is the eighth."""
    import re
    import shutil
    import subprocess
    pdf = EDU_RAW / "trends-college-pricing-2018-full-report.pdf"
    if not pdf.exists():
        try:
            import requests
            EDU_RAW.mkdir(parents=True, exist_ok=True)
            r = requests.get(COLLEGE_BOARD_2018["url"], timeout=180, headers={"User-Agent": "Mozilla/5.0"})
            r.raise_for_status()
            pdf.write_bytes(r.content)
        except Exception as exc:                                              # noqa: BLE001
            return {"verified": None, "note": f"report not cached and not downloadable: {type(exc).__name__}"}
    if shutil.which("pdftotext") is None:
        return {"verified": None, "note": "pdftotext not available"}
    txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True, check=True).stdout.splitlines()
    rows = {"published_tuition_fees": "Published Tuition and Fees", "published_tfrb": "and Room and Board (TFRB)", "net_tuition_fees": "Net Tuition and Fees",
            "net_tfrb": "Net TFRB", "grant_aid_tax_benefits": "per Student"}
    found: dict[str, dict[str, int]] = {}
    head = "Average Published and Net Prices and Average Grant Aid and Tax Benefits per Student in 2018 Dollars"
    for i, ln in enumerate(txt):
        if not ln.strip().startswith(head):
            continue
        title = " ".join(x.strip() for x in txt[i:i + 3])
        sector = ("public_4yr_in_state" if "Public Four-Year" in title else "private_nonprofit_4yr" if "Nonprofit Four-Year" in title else None)
        if sector is None:
            continue
        got: dict[str, int] = {}
        for ln2 in txt[i:i + 14]:
            for key, lab in rows.items():
                if ln2.strip().startswith(lab):
                    nums = [int(x.replace(",", "")) for x in re.findall(r"\$([\d,]+)", ln2)]
                    if len(nums) == 9 and key not in got:
                        got[key] = nums[7]
        found[sector] = got
    checks = [(sec, k, v, found.get(sec, {}).get(k)) for sec in ("public_4yr_in_state", "private_nonprofit_4yr")
              for k, v in COLLEGE_BOARD_2018[sec].items()]
    ok = all(g == v for _, _, v, g in checks)
    return {"verified": bool(ok), "n_checked": len(checks), "mismatches": [(a, b, v, g) for a, b, v, g in checks if g != v]}


# BLS Consumer Expenditure Survey 2018, Table R-1 (as housing_costs.py): mean annual expenditures per consumer unit (CU) and the
# average CU: 2.5 people, 0.6 children under 18, 1.3 earners (one decimal as published).  A CU is a household-like unit; the
# model's agent is an adult, so a cost paid by one person (tuition) is measured against expenditures per adult (or per earner),
# not per CU.  Housing costs, which scale with the size of the CU, are unit-free ratios (housing_costs.md); tuition is not.
CE_URL = "https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf"
CE_TOTAL_EXP = 61224.13
CE_TOTAL_EXP_SE = 775.65
CE_PERSONS, CE_CHILDREN, CE_EARNERS = 2.5, 0.6, 1.3
CE_ADULTS = CE_PERSONS - CE_CHILDREN                    # 1.9 adults (18 and over); the one-decimal inputs leave a range of about 1.8-2.0
CE_COLLEGE_TUITION, CE_COLLEGE_TUITION_SE = 798.09, 55.80     # mean annual college tuition per CU, 2018 ("College tuition" under Education)
CE_UNITS = {"per_cu": 1.0, "per_adult": CE_ADULTS, "per_earner": CE_EARNERS}   # persons (adults, earners) per CU: the divisor of CE_TOTAL_EXP
WORKING_YEARS = 40
DISCOUNT = 0.03
BA_YEARS = 4


def annuity(r: float, n: int = WORKING_YEARS) -> float:
    """Present value of one dollar a year for n years, paid at the end of each year."""
    return float(n) if r == 0 else float((1 - (1 + r) ** (-n)) / r)


def compound_to_entry(r: float, years: int = BA_YEARS) -> float:
    """Value at the start of working life of one dollar paid at the start of each of the `years` college years."""
    return float(years) if r == 0 else float((1 + r) * ((1 + r) ** years - 1) / r)


def model_units(annual: float, years: int = BA_YEARS, ce: float = CE_TOTAL_EXP, r: float = DISCOUNT) -> dict[str, float]:
    """A cost of `annual` dollars a year for `years` years in the model's units of mean annual expenditures `ce`.
    (i) annual / ce: right for a cost paid EVERY year of adult life, like the residence cost r_p;
    (ii) total / working-life consumption: right for a one-off cost, because the model's single adult period stands for
    the whole working life (40 years of consumption ce): undiscounted, with the 40 years discounted at r, and with the
    college payments also compounded at r to the start of working life."""
    tot = annual * years
    return {"annual": annual, "total": tot,
            "i_annual_over_ce": annual / ce,
            "ii_total_over_40ce": tot / (WORKING_YEARS * ce),
            "ii_total_over_pv40ce_r": tot / (annuity(r) * ce),
            "ii_compounded_total_over_pv40ce_r": annual * compound_to_entry(r, years) / (annuity(r) * ce)}


def model_units_by_unit(annual: float, years: int = BA_YEARS, r: float = DISCOUNT) -> dict[str, dict[str, float]]:
    """`model_units` against expenditures per CU, per adult (CE / 1.9) and per earner (CE / 1.3)."""
    return {u: model_units(annual, years, CE_TOTAL_EXP / n, r) for u, n in CE_UNITS.items()}


def verify_ce_table() -> dict[str, Any]:
    """Re-read the constants above from the cached BLS Table R-1 PDF (`pdftotext -layout`; downloaded once if missing) and compare."""
    import re
    import shutil
    import subprocess
    pdf = EDU_RAW / "cu-all-detail-2018.pdf"
    if not pdf.exists():
        try:
            import requests
            EDU_RAW.mkdir(parents=True, exist_ok=True)
            r = requests.get(CE_URL, timeout=180, headers={"User-Agent": "Mozilla/5.0"})
            r.raise_for_status()
            pdf.write_bytes(r.content)
        except Exception as exc:                                              # noqa: BLE001
            return {"verified": None, "note": f"table not cached and not downloadable: {type(exc).__name__}"}
    if shutil.which("pdftotext") is None:
        return {"verified": None, "note": "pdftotext not available"}
    txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True, check=True).stdout.splitlines()

    def last_number(line: str) -> float | None:
        tail = re.split(r"\.{3,}", line)[-1]
        nums = re.findall(r"\$?(\d[\d,]*\.\d+|\d[\d,]*|\.\d+)", tail)
        return float(nums[-1].replace(",", "")) if nums else None

    def find(label: str, start: int = 0) -> int:
        return next(i for i in range(start, len(txt)) if txt[i].strip().startswith(label))

    got: dict[str, float | None] = {}
    got["persons"] = last_number(txt[find("People [I]")])
    got["children_under_18"] = last_number(txt[find("Children under 18 [I]")])
    got["earners"] = last_number(txt[find("Earners [I]")])
    i = find("Average annual expenditures")
    got["total_expenditures"] = next((last_number(t) for t in txt[i + 1:i + 4] if "$" in t), None)
    got["total_expenditures_se"] = next((last_number(t) for t in txt[i + 1:i + 6] if t.strip().startswith("SE")), None)
    j = find("College tuition [I]")
    got["college_tuition"] = next((last_number(t) for t in txt[j + 1:j + 6] if t.strip().startswith("Mean")), None)
    got["college_tuition_se"] = next((last_number(t) for t in txt[j + 1:j + 6] if t.strip().startswith("SE")), None)
    want = {"persons": CE_PERSONS, "children_under_18": CE_CHILDREN, "earners": CE_EARNERS, "total_expenditures": CE_TOTAL_EXP,
            "total_expenditures_se": CE_TOTAL_EXP_SE, "college_tuition": CE_COLLEGE_TUITION, "college_tuition_se": CE_COLLEGE_TUITION_SE}
    bad = {k: (want[k], got.get(k)) for k in want if got.get(k) is None or abs(got[k] - want[k]) > 1e-9}
    return {"verified": not bad, "n_checked": len(want), "read": got, "mismatches": bad, "url": CE_URL}


def tuition_block(ba_share: float) -> dict[str, Any]:
    h330, h331, h303 = nces_html("330.10"), nces_html("331.30"), nces_html("303.70")
    t330 = parse_330_10(h330)
    y = "2017-18"
    pub, allc, nonp = t330["public"][y], t330["all"][y], t330["nonprofit"][y]
    t331 = parse_331_30(h331)
    e303 = parse_303_70(h303)
    cb = COLLEGE_BOARD_2018
    w_pub = e303["public"] / (e303["public"] + e303["nonprofit"])
    nt_pub = cb["public_4yr_in_state"]["net_tuition_fees"]
    nt_np = cb["private_nonprofit_4yr"]["net_tuition_fees"]
    measures: dict[str, dict[str, Any]] = {
        "sticker_tf_public_4yr_in_state": {"label": "published tuition and required fees, public 4-year (in-state)",
                                           "annual": pub["tuition_fees_4yr"], "source": "NCES Digest 2018 Table 330.10, 2017-18"},
        "sticker_tf_all_4yr": {"label": "published tuition and required fees, all 4-year institutions",
                               "annual": allc["tuition_fees_4yr"], "source": "NCES Digest 2018 Table 330.10, 2017-18"},
        "sticker_tf_private_nonprofit_4yr": {"label": "published tuition and required fees, private nonprofit 4-year",
                                             "annual": nonp["tuition_fees_4yr"], "source": "NCES Digest 2018 Table 330.10, 2017-18"},
        "sticker_tfrb_public_4yr_in_state": {"label": "published tuition, fees, room and board, public 4-year (in-state)",
                                             "annual": pub["total_4yr"], "source": "NCES Digest 2018 Table 330.10, 2017-18"},
        "sticker_tfrb_all_4yr": {"label": "published tuition, fees, room and board, all 4-year institutions",
                                 "annual": allc["total_4yr"], "source": "NCES Digest 2018 Table 330.10, 2017-18"},
        "net_tf_public_4yr_in_state": {"label": "net tuition and fees, public 4-year in-state (after grants and tax benefits)",
                                       "annual": float(nt_pub), "source": "College Board Trends in College Pricing 2018, Figure 9, 2017-18 (2018 dollars)"},
        "net_tf_private_nonprofit_4yr": {"label": "net tuition and fees, private nonprofit 4-year",
                                         "annual": float(nt_np), "source": "College Board Trends in College Pricing 2018, Figure 10, 2017-18 (2018 dollars)"},
        "net_tf_blend_4yr": {"label": "net tuition and fees, public and private nonprofit 4-year averaged with fall-2017 undergraduate "
                                      f"enrollment weights ({w_pub:.3f} public)",
                             "annual": float(w_pub * nt_pub + (1 - w_pub) * nt_np),
                             "source": "College Board (prices) and NCES Digest 2018 Table 303.70 (weights; for-profit 4-year, 6% of 4-year undergraduates, excluded)"},
        "net_price_ipeds_public_4yr": {"label": "net price (cost of attendance minus grants), public 4-year, Title IV aid recipients",
                                       "annual": t331["4yr"]["net_price"]["All income levels"][1],
                                       "source": "NCES Digest 2019 Table 331.30, 2017-18; includes room, board, books and other living costs"},
        "net_price_ipeds_all_4yr": {"label": "net price (cost of attendance minus grants), all 4-year, Title IV aid recipients",
                                    "annual": t331["4yr"]["net_price"]["All income levels"][0],
                                    "source": "NCES Digest 2019 Table 331.30, 2017-18; includes room, board, books and other living costs"},
    }
    for k, m in measures.items():
        by = model_units_by_unit(m["annual"])
        by5 = model_units_by_unit(m["annual"], years=5)
        m["model_units"] = by["per_cu"]                       # per consumer unit (the first release's denominator; understates a per-person cost)
        m["model_units_per_adult"], m["model_units_per_earner"] = by["per_adult"], by["per_earner"]
        m["model_units_5yr"], m["model_units_5yr_per_adult"] = by5["per_cu"], by5["per_adult"]
    # aggregate-flow cross-check, unit-free: household college tuition over household expenditures in the CE, divided by the
    # BA+ share, is the per-college-goer cost relative to annual consumption in a stationary economy (aggregate flow = share x cost)
    share = CE_COLLEGE_TUITION / CE_TOTAL_EXP
    share_se = share * float(np.hypot(CE_COLLEGE_TUITION_SE / CE_COLLEGE_TUITION, CE_TOTAL_EXP_SE / CE_TOTAL_EXP))     # ignores the covariance
    flow = {"college_tuition_per_cu": CE_COLLEGE_TUITION, "college_tuition_se": CE_COLLEGE_TUITION_SE, "total_expenditures": CE_TOTAL_EXP,
            "share_of_expenditures": share, "share_se_approx": share_se, "ba_share": ba_share,
            "implied_Fc": share / ba_share, "implied_Fc_se_approx": share_se / ba_share,
            "note": "all postsecondary tuition (two-year, four-year, graduate and professional, and non-completers) paid by households in 2018, per BA+ "
                    "holder; the SE ignores the covariance of the numerator and denominator and the sampling error of the BA share"}
    return {"ce_total_expenditures": CE_TOTAL_EXP, "ce_total_expenditures_se": CE_TOTAL_EXP_SE, "ce_persons_per_cu": CE_PERSONS,
            "ce_children_per_cu": CE_CHILDREN, "ce_adults_per_cu": CE_ADULTS, "ce_earners_per_cu": CE_EARNERS, "ce_url": CE_URL,
            "ce_table_check": verify_ce_table(), "ce_flow_check": flow,
            "denominators": {u: {"persons_per_cu": n, "expenditures": CE_TOTAL_EXP / n} for u, n in CE_UNITS.items()},
            "working_years": WORKING_YEARS, "discount_rate": DISCOUNT, "years_of_college": BA_YEARS,
            "annuity_factor_r": annuity(DISCOUNT), "compounding_factor_r": compound_to_entry(DISCOUNT),
            "nces_330_10_2017_18": {k: v[y] for k, v in t330.items() if y in v},
            "nces_331_30_2017_18_4yr": t331.get("4yr"), "nces_303_70_fall_2017_4yr_undergraduates": e303,
            "college_board_2018": cb, "college_board_check": verify_college_board(), "weight_public": w_pub, "measures": measures,
            "notes_330_10": table_notes(h330), "notes_331_30": table_notes(h331)}


# =============================================================================
# 4. Reproduction of the published ACS block, orchestration
# =============================================================================

DEFINITIONS = {
    "acs_sample": "ACS 2013 5-year PUMS (data years 2009-13), ages 25-34, SCHL >= 12, armed forces and unemployed with 15+ usual hours dropped, "
                  "group quarters in (acs_occupations.build_sample, rules='workbook'); each person has total weight 1 split between the reported occupation "
                  "(u = 1, or 0.5 for 15-29 usual hours) and home production",
    "ba": "ACS: SCHL >= 21 (bachelor's, master's, professional, doctorate); graduate: SCHL >= 22; BA only: SCHL == 21. NLSY79: highest grade completed >= 16 at some round at 25-34",
    "k12_teachers": "Census occupation codes 2300-2340 (pre-K and kindergarten, elementary and middle, secondary, special education, other teachers and instructors); "
                    "core = 2310, 2320, 2330",
    "years_of_schooling": "HHJK highgrade from SCHL (grade 9 = 9, diploma/GED 12, some college 13, associate 14, bachelor 16, graduate 19) bounded to [9, 17]",
    "wage_sample": "employed, 30+ usual hours, 48+ weeks, wage and business income >= $1,000 in 2010 dollars; log hourly wage in 1999 dollars (workbook rules)",
    "conventions": {"P": "population weights: every (person, cell) row weighted by PWGTP x u",
                    "T2": "occupation x gender cell means weighted by PWGTP x u, combined with the unweighted Table 2 counts (sums of u) of the population in question; "
                          "SEs hold the counts fixed",
                    "wsw": "cell means combined with wage-sample record counts", "PW": "cell means combined with PWGTP-weighted counts",
                    "age": "PWGTP-weighted regression with single-year age dummies (and a female dummy when pooled); age_occupation adds gender x occupation fixed effects",
                    "pooled teacher premium": "average of the men's and women's premia with the K-12 wage-sample weights 0.2444 / 0.7556; P_nogenderadj is the unadjusted pooled gap"},
    "acs_stat_keys": "family|group|population|convention; group in men, women, pooled; family in ba_share, years, years_ba, years_nba, s_ratio, grad_share, grad_share_of_ba, "
                     "teacher_share, k12_code_share, k12_code_ba_share, mean_lw_* (mean_lw_teacher, _teacher_ba, _teacher_core), ba_premium, grad_teacher_premium, "
                     "teacher_premium (all K-12 teachers), teacher_premium_ba (BA+ teachers), teacher_premium_core (codes 2310-2330) with population all_nonteachers or ba_nonteachers "
                     "and convention T2, wsw, PW, P, age, P_nogenderadj (pooled); teacher_ba_gap|group|{raw,age,age_code}|P (BA+ minus non-BA K-12 teachers); "
                     "excl|group|{teacher,teacher_ba,nonteacher,nonteacher_ba}|{in_wage_sample,wkw1..wkw6,wkw_missing,wk48_outside} (PWGTP shares of full-time workers by weeks-worked bracket)",
    "weeks_rules": "acs.weeks_sensitivity.{wks40,wks27}: the wage statistics recomputed with the wage sample also admitting WKW bracket 3 (40-47 weeks, weeks = 43.5) / brackets 3 and 4 (27-39 weeks, weeks = 33), "
                   "for teachers and non-teachers alike (workbook wage formula; the 48+ week records are unchanged and checked against build_sample)",
    "nlsy_stat_keys": "sample|group|stat with stat in ba_share, mean_ba, mean_nba, gap (latent SD), gap_obs; terc|sample|group|ba_share|q (observed terciles); terc_fit (the observed-score logit's fitted shares), "
                      "terc_acs (that logit with its level shifted to the ACS BA share), terc_latent / terc_latent_acs (terciles of the latent factor from the errors-in-variables logit, at the NLSY79 / ACS share); "
                      "logit|sample|group|{slope_obs, slope_latent_approx (= slope_obs / sqrt(R)), slope_latent_eiv, intercept_latent_eiv, gap_latent_eiv, ...}; "
                      "teacher|group|{nt,col,cc}|{gap,mean_t,mean_c} (the published teacher gaps)",
    "tuition_units": "model_units: against CE expenditures per consumer unit (the first release); model_units_per_adult: per adult (CE / 1.9, central); model_units_per_earner: per earner (CE / 1.3); each with "
                     "i_annual_over_ce, ii_total_over_40ce (undiscounted, central), ii_total_over_pv40ce_r (3%), ii_compounded_total_over_pv40ce_r; ce_flow_check is the unit-free aggregate-flow check",
    "nlsy_samples": BA_SAMPLES,
}


def acs_reproduction(A: AcsData) -> dict[str, Any]:
    """The sample built here against `estimates/acs_occupations.json` (the published 2009-13 block) and against
    the teacher-premium targets coded in `spatial_calibrate.jl`."""
    J = json.loads((ao.EST / "acs_occupations.json").read_text())["samples"]["acs5_2009_13"]
    w_s, w_w = A.pweights(A.pw)
    S = {k: np.bincount(A.s_cell, weights=w_s * v, minlength=NC) for k, v in A.svars.items()}
    Wc = np.bincount(A.w_cell, weights=w_w, minlength=NC)
    Lc = np.bincount(A.w_cell, weights=w_w * A.w_lw, minlength=NC)
    d_cnt, d_wc, d_lw, d_yr = [], [], [], []
    for gi_, g in enumerate(("M", "F")):
        for grp in GROUPS:
            c = gi_ * NG + GI[grp]
            d_cnt.append(abs(A.cnt["one"][c] - J["counts"][g][grp]))
            d_yr.append(abs(S["yrs"][c] / S["one"][c] - J["years_schooling"][g][grp]))
            if grp != HOME:
                d_wc.append(abs(A.wage_count[c] - J["wage_sample_counts"][g][grp]))
                d_lw.append(abs(Lc[c] / Wc[c] - J["mean_log_hourly_wage_1999usd"][g][grp]))
    st = acs_stats(A, A.pw)
    sc = J["schooling_k12_vs_nonteachers"]
    return {"max_abs_diff_table2_counts": float(max(d_cnt)), "max_abs_diff_wage_sample_counts": float(max(d_wc)),
            "max_abs_diff_mean_log_wage_by_cell": float(max(d_lw)), "max_abs_diff_mean_years_by_cell": float(max(d_yr)),
            "teacher_premium_men_T2": st["teacher_premium|men|all_nonteachers|T2"], "teacher_premium_men_target": 0.0557,
            "teacher_premium_women_T2": st["teacher_premium|women|all_nonteachers|T2"], "teacher_premium_women_target": 0.0860,
            "teacher_premium_men_wsw": st["teacher_premium|men|all_nonteachers|wsw"],
            "teacher_premium_women_wsw": st["teacher_premium|women|all_nonteachers|wsw"],
            "teacher_share_men_T2": st["teacher_share|men|T2"], "teacher_share_women_T2": st["teacher_share|women|T2"],
            "published_teacher_share_men": J["shares"]["M"][K12], "published_teacher_share_women": J["shares"]["F"][K12],
            "years_k12_men": st["years|men|teacher|P"], "published_years_k12_men": sc["men"]["k12_teachers"],
            "years_nonteachers_men": st["years|men|nonteacher|P"], "published_years_nonteachers_men": sc["men"]["all_nonteachers_incl_home"],
            "years_k12_women": st["years|women|teacher|P"], "published_years_k12_women": sc["women"]["k12_teachers"],
            "years_nonteachers_women": st["years|women|nonteacher|P"], "published_years_nonteachers_women": sc["women"]["all_nonteachers_incl_home"],
            "years_nonteachers_pooled": st["years|pooled|nonteacher|P"], "published_years_nonteachers_pooled": sc["pooled"]["all_nonteachers_incl_home"],
            "workbook_years_k12": 16.04, "workbook_years_nonteachers": 13.64}


def run_acs(n_rep: int = N_REP) -> dict[str, Any]:
    t0 = time.time()
    A = load_acs(need_repw=True)
    LOG.info("ACS sample: %d persons, %d in the wage sample (%.0fs)", A.n, len(A.wi), time.time() - t0)
    point = acs_stats(A, A.pw)
    if n_rep != N_REP:
        A.rep = A.rep[:n_rep]
    se = replicate_se(A, point) if n_rep == N_REP else _quick_se(A, point, n_rep)
    stats = {k: {"value": float(v), "se": float(se[k])} for k, v in point.items()}
    # sensitivity of the wage statistics to the weeks rule: also admit the 40-47 (and 27-39) week brackets
    sens: dict[str, Any] = {}
    for tag, extra in WEEKS_RULES.items():
        A2 = AcsData(alt_wage_sample(A.d, extra), A.rep)
        pt = acs_stats(A2, A2.pw, wage_only=True)
        se2 = replicate_se(A2, pt, wage_only=True) if n_rep == N_REP else _quick_se(A2, pt, n_rep, wage_only=True)
        sens[tag] = {"extra_weeks": {WKW_LABEL[b]: wk for b, wk in extra.items()}, "n_wage_sample": int(len(A2.wi)),
                     "counts": acs_counts(A2)["wage_sample"], "stats": {k: {"value": float(v), "se": float(se2[k])} for k, v in pt.items()}}
        LOG.info("weeks-rule sensitivity %s: %d records in the wage sample (%.0fs)", tag, len(A2.wi), time.time() - t0)
        del A2
    return {"stats": stats, "counts": acs_counts(A), "reproduction": acs_reproduction(A), "weeks_sensitivity": sens,
            "n_persons": A.n, "n_wage_sample": int(len(A.wi)), "n_replicates": n_rep}


def _quick_se(A: AcsData, point: dict[str, float], n_rep: int, **kw: Any) -> dict[str, float]:
    """Replicate SE with fewer than 80 replicates (development only)."""
    keys = list(point)
    acc = np.zeros(len(keys))
    for r in range(n_rep):
        st = acs_stats(A, A.rep[r], **kw)
        acc += (np.array([st[k] for k in keys]) - np.array([point[k] for k in keys])) ** 2
    return dict(zip(keys, np.sqrt(4.0 / n_rep * acc).tolist()))


# =============================================================================
# 5. Moments for the model (education_targets.toml)
# =============================================================================

# Existing calibration constants (spatial_calibrate.jl): Mincer workbook schooling and the teacher-premium targets
EXISTING = {"years_nonteach": 13.6419, "years_teach": 16.0437, "horizon": 25.0}
SOURCES = {
    "acs": "ACS 2009-13 5-year PUMS (data years 2009-13; Census API https://api.census.gov/data/2013/acs/acs5/pums), ages 25-34, Table 2 universe (acs_occupations.py sample, HHJK rules); SEs from the 80 Census replicate weights",
    "nlsy79": "NLSY79 block A, ages 25-34 (nlsy_teachers.py / nlsy_ability.py machinery), household bootstrap B = 500, seed 20260929",
    "mincer": "data/Mincer_YrsSchool.xlsx, sheet years_of_school, 2013 column (the constants of spatial_calibrate.jl SCHOOLING); PUMS cross-checks from acs",
    "nces_211_20": "NCES Digest of Education Statistics 2022, Table 211.20, https://nces.ed.gov/programs/digest/d22/tables/dt22_211.20.asp (2017-18 NTPS), current dollars; delta-method SEs",
    "tuition": "College Board, Trends in College Pricing 2018 (Ma, Pender and Welch), Figures 9 and 10, printed pp. 18-19, https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf (net tuition and fees); "
               "NCES Digest 2018 Table 330.10 https://nces.ed.gov/programs/digest/d18/tables/dt18_330.10.asp and Table 303.70 https://nces.ed.gov/programs/digest/d18/tables/dt18_303.70.asp; "
               "NCES Digest 2019 Table 331.30 https://nces.ed.gov/programs/digest/d19/tables/dt19_331.30.asp; "
               "BLS Consumer Expenditure Survey 2018 Table R-1 https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf (expenditures, people, children and earners per consumer unit, college tuition)",
}
GEND = (("pooled", ""), ("men", "_m"), ("women", "_f"))


def build_targets(P: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """The moments of DESIGN.md M3 as {key: entry}, each with role, value, se, definition, the model counterpart (`model`) and the
    name of the model-side moment (`model_moment`)."""
    S, N = P["acs"]["stats"], P["nlsy"]["stats"]
    cnt = P["acs"]["counts"]
    sal = P["nces_salary"]["years"]["2017-18"]
    tu = P["tuition"]
    WS = P["acs"].get("weeks_sensitivity", {})
    T: dict[str, dict[str, Any]] = {}

    def A(k: str) -> tuple[float, float]:
        return S[k]["value"], S[k]["se"]

    def wk(tag: str, k: str) -> float | None:
        return WS[tag]["stats"][k]["value"] if tag in WS and k in WS[tag]["stats"] else None

    lab_g = {"pooled": "men and women pooled", "men": "men", "women": "women"}
    ba_choice = ("Alternative forms of ONE moment (choice_group ba_premium): use one of them in a criterion. ba_premium (raw) is the total return to a BA, including the "
                 "BA holders' sorting into better-paid occupations, which the model's non-teaching occupation choice (independent of d) cannot generate, so rho_c would "
                 "absorb it; ba_premium_within_occupation is the premium within occupation, which is the model's premium by construction, but it leaves the model's "
                 "total return to college below the data's. Which one is primary is the M3 designer's choice.")
    for g, suf in GEND:
        lab = lab_g[g]
        role = "target" if g == "pooled" else "target_by_gender"
        v, se = A(f"ba_share|{g}|nonteacher|T2")
        T[f"ba_share_nonteacher{suf}"] = {
            "role": role, "value": v, "se": se, "scale": 0.05 * v, "scale_rule": "5% of the value (modelling tolerance; the sampling SE is far smaller)",
            "n_universe": int(round(cnt["nonteachers_u"][g])), "n_universe_note": "Table 2 person-equivalents of non-teachers (home production included)",
            "units": "share", "convention": "T2", "source": "acs", "model_moment": "ba_share_nonteach",
            "definition": f"share with a bachelor's degree or more (SCHL >= 21) among non-teachers aged 25-34, {lab}; non-teachers = every "
                          "non-K-12-teaching occupation plus home production (Table 2 universe); cell shares weighted by PWGTP, combined with the Table 2 counts",
            "model": "share of adults with d = 1 among non-teachers (all non-teaching occupations and home production)" + ("" if g == "pooled" else f", {lab}"),
            "alt_population_weights": A(f"ba_share|{g}|nonteacher|P")[0], "alt_market_workers_only_T2": A(f"ba_share|{g}|market_nonteacher|T2")[0]}
        v, se = A(f"ba_premium|{g}|raw|T2")
        T[f"ba_premium{suf}"] = {
            "role": role, "choice_group": "ba_premium", "value": v, "se": se, "scale": 0.05, "scale_rule": "0.05 (modelling tolerance, as the teacher premia)",
            "n_wage_sample": cnt["wage_sample"][g]["nonteachers"], "units": "log points", "convention": "T2", "source": "acs", "model_moment": "ba_premium",
            "definition": f"mean log hourly wage of non-teacher market workers with a BA or more minus those without, {lab}, full-time full-year wage sample; "
                          "occupation x gender cell means (PWGTP) weighted by the Table 2 counts of each education group (raw: no occupation control)",
            "model": "mean log wage of d = 1 minus d = 0 among non-teacher market workers (home production out), cells weighted by the model's shares; the model's "
                     "non-teaching occupation choice does not depend on d, so this equals its within-occupation premium",
            "tradeoff": ba_choice, "alt_population_weights": A(f"ba_premium|{g}|raw|P")[0], "alt_age_controlled": A(f"ba_premium|{g}|age|P")[0],
            "alt_age_and_occupation_controlled": A(f"ba_premium|{g}|age_occupation|P")[0], "alt_weeks_40plus": wk("wks40", f"ba_premium|{g}|raw|T2")}
        v, se = A(f"ba_premium|{g}|age_occupation|P")
        T[f"ba_premium_within_occupation{suf}"] = {
            "role": role, "choice_group": "ba_premium", "value": v, "se": se, "scale": 0.05, "scale_rule": "0.05 (the existing practice for premia)",
            "n_wage_sample": cnt["wage_sample"][g]["nonteachers"], "units": "log points", "convention": "P (regression)", "source": "acs", "model_moment": "ba_premium",
            "definition": f"coefficient on BA+ in a PWGTP-weighted regression of log hourly wage on BA+, single-year age dummies and occupation x gender fixed effects, "
                          f"non-teacher market workers, {lab}: the BA premium within occupation",
            "model": "ba_premium when occupation cells are held fixed: the model has one credential premium rho_c in every non-teaching occupation and no BA-specific "
                     "occupational sorting, so this is its counterpart of rho_c plus selection on z",
            "tradeoff": ba_choice, "alt_weeks_40plus": wk("wks40", f"ba_premium|{g}|age_occupation|P")}
        key = f"all_hgc9|{g}|gap"
        T[f"afqt_gap_ba{suf}"] = {
            "role": role, "value": N[key]["value"], "se": N[key]["se"], "scale": 0.05,
            "scale_rule": "0.05 (the bootstrap SE is smaller; AFQT was taken at 15-23 and partly reflects schooling, and the cohort is 1957-64)",
            "n_persons": P["nlsy"]["counts"]["all_hgc9"][g]["persons"], "units": "latent AFQT SDs", "convention": "weighted means, latent SD = observed / sqrt(R_comp_pop)",
            "source": "nlsy79", "model_moment": "ba_z_gap_sd",
            "definition": f"mean latent AFQT of persons with 16+ years of schooling minus those with less (highest grade >= 9), {lab}, observed at 25-34",
            "model": "mean log z of d = 1 adults minus d = 0 adults, in SDs of log z (as teacher_z_gap_market_sd; stationary s_z = 0.1955), all adults, home production included",
            "alt_all_grades": N[f"all|{g}|gap"]["value"], "alt_employed": N[f"emp_hgc9|{g}|gap"]["value"],
            "alt_employed_nonteachers": N[f"emp_nonteacher_hgc9|{g}|gap"]["value"],
            "alt_at_acs_ba_share": N[f"logit|all_hgc9|{g}|gap_at_acs_share"]["value"]}
    v, se = A("ba_share|pooled|all|P")
    T["ba_share_all"] = {
        "role": "validation", "value": v, "se": se, "n_universe": cnt["persons"]["pooled"], "units": "share", "convention": "P", "source": "acs",
        "definition": "share with a bachelor's degree or more (SCHL >= 21) among all persons aged 25-34 of the Table 2 universe (home production included), PWGTP weights",
        "model": "share of adults with d = 1 (teachers and non-teachers)", "model_moment": "ba_share",
        "alt_men": A("ba_share|men|all|P")[0], "alt_women": A("ba_share|women|all|P")[0]}
    v, se = A("ba_premium|pooled|age|P")
    T["ba_premium_age_adj"] = {
        "role": "target_alt", "choice_group": "ba_premium", "value": v, "se": se, "scale": 0.05, "scale_rule": "0.05", "n_wage_sample": cnt["wage_sample"]["pooled"]["nonteachers"],
        "units": "log points", "convention": "P (regression)", "source": "acs", "model_moment": "ba_premium",
        "definition": "coefficient on BA+ in a PWGTP-weighted regression of log hourly wage on BA+, single-year age dummies and a female dummy, non-teacher market workers (no occupation control)",
        "model": "as ba_premium, with the model's gender composition entering through the female dummy", "tradeoff": ba_choice}
    v = EXISTING["years_nonteach"] / EXISTING["horizon"]
    T["s_O_mean"] = {
        "role": "target", "value": v, "scale": 0.005, "scale_rule": "0.005 (0.125 years of schooling): about three times the gap between the workbook and the PUMS T2 mean (0.0015), "
        "far above the PUMS sampling SE; phi is re-fit to hit the value, so the scale matters only if the moment enters a criterion",
        "units": "share of 25 years", "convention": "Mincer workbook, Table 2 counts (existing target, spatial_calibrate.jl SCHOOLING)", "source": "mincer", "model_moment": "s_O_mean",
        "definition": "mean years of schooling (capped 9-17) of non-teachers, home production included, over 25: 13.6419 / 25",
        "model": "mean time share s of non-teachers; re-fit by phi (calibrate_phi)",
        "note": "the workbook value cannot be reproduced exactly from PUMS (RMSE 0.059 years by cell); the PUMS SE of the T2 mean is in se_pums_T2",
        "se_pums_T2": A("years|pooled|nonteacher|T2")[1] / EXISTING["horizon"],
        "alt_pums_T2": A("years|pooled|nonteacher|T2")[0] / EXISTING["horizon"], "alt_pums_population_weights": A("years|pooled|nonteacher|P")[0] / EXISTING["horizon"]}
    # ---- validation
    T["s_ratio"] = {
        "role": "validation", "value": EXISTING["years_teach"] / EXISTING["years_nonteach"], "units": "ratio", "se_pums_T2": A("s_ratio|pooled|T2")[1],
        "convention": "Mincer workbook (existing); PUMS in alt", "source": "mincer", "model_moment": "s_ratio",
        "definition": "mean years of schooling of K-12 teachers over non-teachers' (16.0437 / 13.6419)", "model": "s_T / s_O (mean time shares)",
        "alt_pums_T2": A("s_ratio|pooled|T2")[0], "alt_pums_population_weights": A("s_ratio|pooled|P")[0]}
    T["s_ba_nonteacher"] = {
        "role": "validation", "value": A("years_ba|pooled|nonteacher|T2")[0] / EXISTING["horizon"], "se": A("years_ba|pooled|nonteacher|T2")[1] / EXISTING["horizon"],
        "units": "share of 25 years", "convention": "T2", "source": "acs", "model_moment": "none (mean s of d = 1 non-teachers)",
        "definition": "mean years of schooling (capped 9-17) of BA+ non-teachers over 25 (a BA is 16, a graduate degree 17 after the cap)",
        "model": "mean s of d = 1 non-teachers (s >= s_BA = 16/25 = 0.64 binds for many)"}
    T["s_nonba_nonteacher"] = {
        "role": "validation", "value": A("years_nba|pooled|nonteacher|T2")[0] / EXISTING["horizon"], "se": A("years_nba|pooled|nonteacher|T2")[1] / EXISTING["horizon"],
        "units": "share of 25 years", "convention": "T2", "source": "acs", "model_moment": "none (mean s of d = 0 non-teachers)",
        "definition": "mean years of schooling (capped 9-17) of non-BA non-teachers over 25", "model": "mean s of d = 0 non-teachers (s unconstrained)"}
    T["s_T_mean_ba"] = {
        "role": "validation", "value": A("years_ba|pooled|teacher|P")[0] / EXISTING["horizon"], "se": A("years_ba|pooled|teacher|P")[1] / EXISTING["horizon"],
        "units": "share of 25 years", "convention": "P", "source": "acs", "model_moment": "s_T_mean_ba",
        "definition": "mean years of schooling (capped 9-17) of K-12 teachers with a BA or more over 25: a BA is 16, a graduate degree 17 after the cap "
                      f"({100 * A('grad_share_of_ba|pooled|teacher|P')[0]:.0f}% of BA+ teachers hold one)",
        "model": "mean s of d = 1 teachers", "alt_men": A("years_ba|men|teacher|P")[0] / EXISTING["horizon"], "alt_women": A("years_ba|women|teacher|P")[0] / EXISTING["horizon"]}
    v, se = A("ba_share|pooled|teacher|P")
    T["teacher_ba_share"] = {
        "role": "validation", "value": v, "se": se, "n_universe": int(round(cnt["teachers_u"]["pooled"])), "n_universe_note": "Table 2 person-equivalents of K-12 teachers",
        "units": "share", "convention": "P", "source": "acs", "model_moment": "ba_share_teach",
        "definition": "BA+ share of K-12 teachers (Census codes 2300-2340: pre-K and kindergarten, elementary and middle, secondary, special education, "
                      "other teachers and instructors)",
        "model": "1 if the licensure rule (s >= 16/25) is imposed on every teacher; the data share is lower because of the pre-K/K and 'other instructors' codes "
                 "(a non-BA teaching track would reproduce it)",
        "alt_core_codes_2310_2330": A("k12_code_ba_share|pooled|core")[0], "alt_prek_kindergarten_2300": A("k12_code_ba_share|pooled|2300")[0],
        "alt_other_instructors_2340": A("k12_code_ba_share|pooled|2340")[0], "core_share_of_teachers": A("k12_code_share|pooled|core")[0]}
    v, se = A("grad_share|pooled|teacher|P")
    T["teacher_graduate_share"] = {
        "role": "validation", "value": v, "se": se, "units": "share", "convention": "P", "source": "acs", "model_moment": "none",
        "definition": "share of K-12 teachers with a master's, professional or doctoral degree (SCHL >= 22), ages 25-34",
        "model": "no graduate margin in the model (outside M3)", "alt_share_of_ba_teachers": A("grad_share_of_ba|pooled|teacher|P")[0]}
    # teacher premia: all K-12 teachers (the existing targets) and the BA+ and core-code teacher groups, against all and against BA+ non-teachers
    tdesc = {"": "all K-12 teachers (codes 2300-2340, including the pre-K/K and other instructors, 12.6% of whom have no BA)",
             "_ba": "K-12 teachers with a BA or more", "_core": "K-12 teachers of the core codes 2310-2330 (96% BA+)"}
    cdesc = {"all_nonteachers": "all non-teacher market workers (cell means weighted by Table 2 counts)",
             "ba_nonteachers": "non-teacher market workers with a BA or more (cell means weighted by the Table 2 counts of the BA+ workers)"}
    tp = [("teacher_premium_all", "", "all_nonteachers", "teacher_premium_m / teacher_premium_f (the existing targets +0.0557 / +0.0860)",
           "teacher_premium_m/f", "all teachers against all non-teachers; the existing targets"),
          ("teacher_premium_ba", "", "ba_nonteachers", "teacher premium relative to d = 1 non-teacher market workers when the model has a non-BA teaching track "
           "(all teachers against d = 1 non-teachers)", "teacher_premium_ba_m/f", "all teachers against BA+ non-teachers"),
          ("teacher_premium_ba_teachers_all", "_ba", "all_nonteachers", "d = 1 teachers against all non-teachers", "teacher_premium_m/f when every teacher has d = 1",
           "BA+ teachers against all non-teachers"),
          ("teacher_premium_within_ba", "_ba", "ba_nonteachers", "d = 1 teachers against d = 1 non-teacher market workers: the counterpart of 'teachers versus d = 1 "
           "non-teachers' when every teacher is forced to d = 1 (the wage analogue of NLSY comparison C, teacher_z_gap_within_ba)", "teacher_premium_ba_m/f when every teacher has d = 1",
           "BA+ teachers against BA+ non-teachers"),
          ("teacher_premium_core_all", "_core", "all_nonteachers", "teachers of the licensure codes against all non-teachers", "none (robustness)",
           "core-code teachers against all non-teachers"),
          ("teacher_premium_core_ba", "_core", "ba_nonteachers", "teachers of the licensure codes against d = 1 non-teachers", "none (robustness)",
           "core-code teachers against BA+ non-teachers")]
    for key, sfx, comp, model, mm, short in tp:
        fam = f"teacher_premium{sfx}"
        for g, suf in GEND:
            v, se = A(f"{fam}|{g}|{comp}|T2")
            T[f"{key}{suf}"] = {
                "role": "validation", "value": v, "se": se, "units": "log points", "convention": "T2", "source": "acs", "model_moment": mm,
                "definition": f"mean log hourly wage of {tdesc[sfx]} minus that of {cdesc[comp]}, {lab_g[g]}"
                              + ("; pooled = average of the two with the K-12 wage-sample weights 0.2444/0.7556" if g == "pooled" else ""),
                "model": model, "comparison": short, "alt_population_weights": A(f"{fam}|{g}|{comp}|P")[0],
                "alt_wage_sample_weights": A(f"{fam}|{g}|{comp}|wsw")[0], "alt_age_controlled": A(f"{fam}|{g}|{comp}|age")[0],
                "alt_weeks_40plus": wk("wks40", f"{fam}|{g}|{comp}|T2"), "alt_weeks_27plus": wk("wks27", f"{fam}|{g}|{comp}|T2")}
    for g, suf in GEND:
        v, se = A(f"teacher_ba_gap|{g}|raw|P")
        T[f"teacher_ba_wage_gap{suf}"] = {
            "role": "validation", "value": v, "se": se, "units": "log points", "convention": "P", "source": "acs", "model_moment": "none (pay penalty of a non-BA teaching track)",
            "definition": f"mean log hourly wage of K-12 teachers with a BA or more minus those without (pre-K/K and other instructors, mostly), {lab_g[g]}, wage sample",
            "model": "mean log wage of d = 1 teachers minus d = 0 teachers when the model has a non-BA teaching track (its pay penalty, rho_T0, in the model's terms)",
            "alt_age_controlled": A(f"teacher_ba_gap|{g}|age|P")[0], "alt_age_and_census_code_controlled": A(f"teacher_ba_gap|{g}|age_code|P")[0]}
    v, se = A("grad_teacher_premium|pooled|age|P")
    T["grad_teacher_premium"] = {
        "role": "validation", "value": v, "se": se, "units": "log points", "convention": "P (regression)", "source": "acs", "model_moment": "none (compare with ma_premium_implied)",
        "definition": "log hourly wage of K-12 teachers with a graduate degree minus teachers with exactly a bachelor's degree, controlling for single-year age (and sex)",
        "model": "none (no graduate margin); compare with the salary-schedule premium", "alt_raw": A("grad_teacher_premium|pooled|raw|P")[0],
        "alt_men_age": A("grad_teacher_premium|men|age|P")[0], "alt_women_age": A("grad_teacher_premium|women|age|P")[0]}
    ag = sal["aggregates"]["1-5 years"]
    T["ma_ba_salary_premium"] = {
        "role": "validation", "value": ag["log"], "se": ag["se_log"], "units": "log points", "convention": "delta method, experience rows weighted by the number of teachers",
        "source": "nces_211_20", "model_moment": "ma_premium_implied",
        "definition": "log of the average base salary of full-time public school teachers with a master's degree over that with a bachelor's degree, 1-5 years of teaching experience",
        "model": "theta_s x log(17.5/16) (the schedule premium implied by the credential-lane elasticity); no graduate margin otherwise", "alt_first_year": sal["rows"]["1 year or less"]["log"],
        "alt_1_3_years": sal["aggregates"]["1-3 years"]["log"], "alt_1_9_years": sal["aggregates"]["1-9 years"]["log"]}
    for key, cl, lab, mm in (("teacher_z_gap_employed_nonteachers", "nt", "employed non-teachers (comparison A)", "teacher_z_gap_market_sd"),
                             ("teacher_z_gap_nonteacher_ba", "col", "non-teacher persons with 16+ years (comparison B)", "teacher_z_gap_ba_sd"),
                             ("teacher_z_gap_within_ba", "cc", "non-teachers with 16+ years, teachers with 16+ years only (comparison C)", "teacher_z_gap_within_ba_sd")):
        r = N[f"teacher|pooled|{cl}|gap"]
        T[key] = {"role": "validation", "value": r["value"], "se": r["se"], "units": "latent AFQT SDs", "convention": "weighted means",
                  "source": "nlsy79", "note": "reproduces nlsy_teachers.json (same draws)", "model_moment": mm,
                  "definition": f"mean latent AFQT of K-12 teachers minus that of {lab}, employed sample, 1980 Census codes, majority-of-rounds rule",
                  "model": "teacher_z_gap_market_sd (mean log z of teachers minus non-teacher workers over s_z)" if cl == "nt" else "teacher z gap against / within BA+ (needs d)"}
    for q in (1, 2, 3):
        r = N[f"terc|all_hgc9|pooled|ba_share|{q}"]
        T[f"p_ba_afqt_tercile_{q}"] = {
            "role": "validation", "value": r["value"], "se": r["se"], "units": "share", "convention": "weighted", "source": "nlsy79", "model_moment": "ba_share_by_z",
            "definition": f"share with 16+ years of schooling in tercile {q} of own OBSERVED AFQT (terciles of the weighted sample of persons observed at 25-34, highest grade >= 9; "
                          f"cut points {N['terc|all_hgc9|cut1']['value']:.3f}, {N['terc|all_hgc9|cut2']['value']:.3f} observed SDs), at the NLSY79 BA share "
                          f"({N['all_hgc9|pooled|ba_share']['value']:.3f}, against {S['ba_share|pooled|all|P']['value']:.3f} in the ACS)",
            "model": "P(d = 1 | tercile) of a NOISY score x = sqrt(R) theta + sqrt(1 - R) u with theta = standardized log z, u standard normal, R = "
                     f"{N['R_comp_pop']['value']:.3f} (terciles of x over adults); for a noise-free model use alt_latent_scale (terciles of log z itself)",
            "alt_at_acs_share": N[f"terc_acs|all_hgc9|pooled|ba_share|{q}"]["value"], "alt_latent_scale": N[f"terc_latent|all_hgc9|pooled|ba_share|{q}"]["value"],
            "alt_latent_scale_at_acs_share": N[f"terc_latent_acs|all_hgc9|pooled|ba_share|{q}"]["value"],
            "alt_note": "alt_at_acs_share: the observed-score logit with its level shifted to the ACS BA share; alt_latent_scale: terciles of the latent factor from the "
                        "errors-in-variables logit (assumes a normal latent factor and a logistic P(BA | theta)); *_at_acs_share: level shifted to the ACS BA share"}
    r = N["terc|all_hgc9|pooled|ba_share|diff"]
    T["p_ba_afqt_top_minus_bottom"] = {"role": "validation", "value": r["value"], "se": r["se"], "units": "share difference", "convention": "weighted",
                                      "source": "nlsy79", "model_moment": "ba_share_by_z", "definition": "top minus bottom AFQT tercile, share with 16+ years of schooling",
                                      "model": "as p_ba_afqt_tercile_q", "alt_at_acs_share": N["terc_acs|all_hgc9|pooled|ba_share|diff"]["value"],
                                      "alt_latent_scale": N["terc_latent|all_hgc9|pooled|ba_share|diff"]["value"],
                                      "alt_latent_scale_at_acs_share": N["terc_latent_acs|all_hgc9|pooled|ba_share|diff"]["value"]}
    r = N["logit|all_hgc9|pooled|slope_obs"]
    T["afqt_logit_slope"] = {
        "role": "validation", "value": r["value"], "se": r["se"], "units": "log-odds per observed AFQT SD", "convention": "weighted logit", "source": "nlsy79",
        "model_moment": "ba_share_by_z (slope form)", "definition": "slope of the logit of P(16+ years of schooling) on observed AFQT (standardized)",
        "model": "slope of logit P(d = 1) on the NOISY score x = sqrt(R) theta + sqrt(1 - R) u (theta = log z standardized, u standard normal, R = "
                 f"{N['R_comp_pop']['value']:.3f}); the observed score is attenuated, so a noise-free model needs afqt_logit_slope_latent instead",
        "alt_latent_first_order": N["logit|all_hgc9|pooled|slope_latent_approx"]["value"], "alt_latent_eiv": N["logit|all_hgc9|pooled|slope_latent_eiv"]["value"]}
    r = N["logit|all_hgc9|pooled|slope_latent_eiv"]
    T["afqt_logit_slope_latent"] = {
        "role": "validation", "value": r["value"], "se": r["se"], "units": "log-odds per latent AFQT SD", "convention": "errors-in-variables logit", "source": "nlsy79",
        "model_moment": "ba_share_by_z (slope form)",
        "definition": "slope of the logit of P(16+ years of schooling) on the latent AFQT factor (SD units), by maximum likelihood with a normal latent factor and the composite "
                      f"reliability R = {N['R_comp_pop']['value']:.3f} (Gauss-Hermite integration over the factor given the score); the first-order correction slope_obs / sqrt(R) "
                      "understates it",
        "model": "d logit P(d = 1 | z) / d (log z / s_z): the slope of the logit of the model's BA probability on standardized log z (weighted over adults)",
        "alt_first_order": N["logit|all_hgc9|pooled|slope_latent_approx"]["value"], "alt_men": N["logit|all_hgc9|men|slope_latent_eiv"]["value"],
        "alt_women": N["logit|all_hgc9|women|slope_latent_eiv"]["value"], "alt_implied_latent_gap": N["logit|all_hgc9|pooled|gap_latent_eiv"]["value"]}
    # ---- cost
    m = tu["measures"]
    cen = m["net_tf_blend_4yr"]
    T["tuition_Fc"] = {
        "role": "cost", "value": cen["model_units_per_adult"]["ii_total_over_40ce"], "se": float("nan"),
        "units": f"units of Cbar = mean annual expenditures per ADULT (${tu['denominators']['per_adult']['expenditures']:,.0f} = CE ${tu['ce_total_expenditures']:,.0f} / {tu['ce_adults_per_cu']:.1f} adults per consumer unit)",
        "convention": "(ii-a) total BA cost over 40 years of consumption per adult, undiscounted", "source": "tuition", "model_moment": "Fc",
        "definition": "central: net tuition and fees (after grants and education tax benefits), public and private nonprofit 4-year weighted by undergraduate enrollment, "
                      "4 years, over 40 years of annual consumption per adult: the undiscounted, per-person counterpart of a one-off goods cost paid from adult consumption in a "
                      "model without an interest rate whose other education parameters are calibrated on annual aggregate flows",
        "model": "Fc, goods cost of a BA in Node.m (and Xi), common to all college goers",
        "annual_dollars": cen["annual"], "total_dollars": cen["model_units_per_adult"]["total"],
        "alt_annuitized_3pct_ii_b": cen["model_units_per_adult"]["ii_total_over_pv40ce_r"],
        "alt_compounded_3pct_ii_c": cen["model_units_per_adult"]["ii_compounded_total_over_pv40ce_r"],
        "alt_aggregate_flow_check": tu["ce_flow_check"]["implied_Fc"], "alt_aggregate_flow_check_se": tu["ce_flow_check"]["implied_Fc_se_approx"],
        "alt_per_earner_ii_a": cen["model_units_per_earner"]["ii_total_over_40ce"], "alt_per_earner_ii_b": cen["model_units_per_earner"]["ii_total_over_pv40ce_r"],
        "alt_per_earner_i": cen["model_units_per_earner"]["i_annual_over_ce"], "alt_per_earner_ii_c": cen["model_units_per_earner"]["ii_compounded_total_over_pv40ce_r"],
        "alt_per_cu_i": cen["model_units"]["i_annual_over_ce"], "alt_per_cu_ii_a": cen["model_units"]["ii_total_over_40ce"], "alt_per_cu_ii_b": cen["model_units"]["ii_total_over_pv40ce_r"],
        "alt_per_cu_ii_c": cen["model_units"]["ii_compounded_total_over_pv40ce_r"],
        "per_cu_note": "per consumer unit (the first release of this file) understates a per-person cost by the number of adults per unit, 1.9",
        "not_for_one_off_costs_i_annual_over_ce_per_adult": cen["model_units_per_adult"]["i_annual_over_ce"],
        "range_low": m["net_tf_public_4yr_in_state"]["model_units_per_adult"]["ii_total_over_40ce"],
        "range_high": m["sticker_tf_all_4yr"]["model_units_per_adult"]["ii_total_over_pv40ce_r"],
        "range_note": "per adult; low = net tuition and fees at public 4-year, 4 years, undiscounted; high = published tuition and fees, all 4-year, discounted at 3%"}
    return T


# =============================================================================
# 6. Writers: JSON, TOML, markdown
# =============================================================================

def clean(obj: Any) -> Any:
    """JSON-safe copy: numpy scalars to Python, NaN/inf to None."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _toml_val(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_val(x) for x in v) + "]"
    return json.dumps(str(v), ensure_ascii=True)


def write_toml(targets: dict[str, dict[str, Any]], path: Path, meta: dict[str, Any]) -> None:
    L = [
        "# education_targets.toml -- moments for the joint education and occupational choice (DESIGN.md section M3), for the Julia calibration.",
        f"# Generated {meta['generated']} by data/spatial/education_occupation.py; see education_occupation.md (definitions, SEs, caveats).",
        "# role: target (the M3 targets), target_by_gender / target_alt (other forms of a target), validation, cost.",
        "# model_moment: the name of the model-side moment (M3 spec) that the data moment is the counterpart of (\"none\" where the model has no counterpart yet).",
        "# choice_group: moments with the same group are alternative forms of ONE moment; a criterion should contain only one of them (the choice is the designer's).",
        "# value, se: point estimate and standard error (ACS: Census replicate weights; NLSY79: household bootstrap B = 500, seed 20260929; NCES: delta method).",
        "# scale: suggested criterion scale of a target (the sampling SE is much smaller than the other errors). alt_*: alternative definitions, weights or samples",
        "# (alt_weeks_40plus / alt_weeks_27plus: the wage sample also admits 40-47 / 27-47 weeks worked; alt_at_acs_share: level-shifted to the ACS BA share).",
        "# convention T2 = cell means weighted by PWGTP, combined with the Table 2 (unweighted) counts, as the existing teacher-premium targets;",
        "# P = population weights. Nothing here is fitted; no value was chosen to match a model output.",
        "",
        "[source]",
    ] + [f"{k} = {_toml_val(v)}" for k, v in SOURCES.items()] + [
        "",
        "[meta]",
        f"generated = {_toml_val(meta['generated'])}",
        'script = "data/spatial/education_occupation.py"',
        f"acs_persons = {int(meta['acs_persons'])}",
        f"nlsy_B = {int(meta['nlsy_B'])}",
        f"nlsy_seed = {int(meta['nlsy_seed'])}",
        "",
        "[units]",
        "# dollars to model units for goods costs: C-bar is mean annual expenditures PER ADULT (CE 2018 Table R-1, per consumer unit / adults per unit); a recurring annual cost is",
        "# annual dollars / C-bar, a one-off cost is total dollars / (working_years x C-bar) undiscounted (or / (annuity factor x C-bar) at the discount rate)",
        f"ce_expenditures_per_consumer_unit = {CE_TOTAL_EXP!r}",
        f"people_per_consumer_unit = {CE_PERSONS!r}",
        f"adults_per_consumer_unit = {CE_ADULTS!r}",
        f"earners_per_consumer_unit = {CE_EARNERS!r}",
        f"cbar_dollars_per_adult = {CE_TOTAL_EXP / CE_ADULTS!r}",
        f"cbar_dollars_per_earner = {CE_TOTAL_EXP / CE_EARNERS!r}",
        f"working_years = {WORKING_YEARS}",
        f"discount_rate = {DISCOUNT!r}",
        f"annuity_factor = {annuity(DISCOUNT)!r}",
        "",
    ]
    order = ["target", "target_by_gender", "target_alt", "validation", "cost"]
    for role in order:
        for k, e in targets.items():
            if e["role"] != role:
                continue
            L.append(f"[moment.{k}]")
            for f, v in e.items():
                if v is None or (isinstance(v, float) and not np.isfinite(v)):
                    continue
                L.append(f"{f} = {_toml_val(v)}")
            L.append("")
    path.write_text("\n".join(L), encoding="utf-8")


def _vs(r: dict[str, float], nd: int = 3, sign: bool = False) -> str:
    """'value (se)' with the SE to one more decimal (at most four); no SE for fixed counts (SE = 0)."""
    if r is None or r.get("value") is None or not np.isfinite(r["value"]):
        return "n/a"
    v = f"{r['value']:+.{nd}f}" if sign else f"{r['value']:.{nd}f}"
    se = r.get("se")
    if se is None or not np.isfinite(se) or se == 0:
        return v
    d = min(nd + 1, 4)
    return f"{v} (<{10.0 ** -d:.{d}f})" if se < 0.5 * 10.0 ** -d else f"{v} ({se:.{d}f})"


def _md(headers: list[str], rows: list[list[str]], left: tuple[int, ...] = ()) -> str:
    """Markdown table; the first column and the columns listed in `left` are left-aligned, the rest right-aligned."""
    align = ["---" if (i == 0 or i in left) else "---:" for i in range(len(headers))]
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join(align) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def write_md(P: dict[str, Any], path: Path) -> None:
    S, N = P["acs"]["stats"], P["nlsy"]["stats"]
    NC_ = P["nlsy"]["counts"]
    cnt = P["acs"]["counts"]
    rep = P["acs"]["reproduction"]
    WS = P["acs"].get("weeks_sensitivity", {})
    sal = P["nces_salary"]
    tu = P["tuition"]
    T = P["targets"]
    meta = P["meta"]
    L: list[str] = []
    ap = L.append

    def s(k, nd=3, sign=False):
        return _vs(S[k], nd, sign)

    def sv(k):
        return S[k]["value"]

    def n_(k, nd=3, sign=False):
        return _vs(N[k], nd, sign)

    def nv(k):
        return N[k]["value"]

    def cols(tmpl: str, nd=3, sign=False, gs=("men", "women", "pooled")):
        return [s(tmpl.format(g=g), nd, sign) for g in gs]

    def wsc(tag, k, nd=3, sign=False):
        return _vs(WS[tag]["stats"][k], nd, sign)

    def wsv(tag, k):
        return WS[tag]["stats"][k]["value"]

    mu = tu["measures"]
    cen = mu["net_tf_blend_4yr"]
    adults = tu["ce_adults_per_cu"]
    ex_ad = tu["denominators"]["per_adult"]["expenditures"]
    flow = tu["ce_flow_check"]
    fc = T["tuition_Fc"]

    ap("# Education and occupation: evidence for a joint BA and occupational choice")
    ap("")
    if meta.get("quick"):
        ap("> **Development run (`--quick`): few replicates and bootstrap draws, the standard errors are not valid.**")
        ap("")
    ap(f"Generated {meta['generated']} by `education_occupation.py` (DESIGN.md section M3, data task D3). Sources are listed below; "
       "SEs are Census successive-difference replicate weights (80 replicates) for the ACS, the household bootstrap "
       f"(B = {meta['nlsy_B']}, seed {meta['nlsy_seed']}, the draws of `nlsy_teachers.py`) for the NLSY79, and the delta method for NCES. "
       "Machine-readable results: `education_occupation.json`; moments for the Julia code: `education_targets.toml`.")
    ap("")
    ap("**Sources.**")
    ap("")
    ap("- ACS 2009-13 5-year PUMS (data years 2009-13), Census API `https://api.census.gov/data/2013/acs/acs5/pums`, cached in `raw/acs_pums/`; sample, HHJK rules and K-12 codes of "
       "`acs_occupations.py`; the 80 replicate weights PWGTP1-PWGTP80 (`--fetch-repw`).")
    ap("- NCES, *Digest of Education Statistics* 2022, Table 211.20 (average base salary of full-time public school teachers by highest degree and years of experience), "
       "<https://nces.ed.gov/programs/digest/d22/tables/dt22_211.20.asp> (2020 edition, `d20/tables/dt20_211.20.asp`, for the check of the 2017-18 rows).")
    ap("- NCES, *Digest* 2018, Table 330.10 (average undergraduate tuition and fees and room and board), <https://nces.ed.gov/programs/digest/d18/tables/dt18_330.10.asp>; "
       "Table 303.70 (undergraduate enrollment), <https://nces.ed.gov/programs/digest/d18/tables/dt18_303.70.asp>; *Digest* 2019, Table 331.30 (net price), "
       "<https://nces.ed.gov/programs/digest/d19/tables/dt19_331.30.asp>.")
    ap("- College Board, *Trends in College Pricing 2018* (Ma, Pender and Welch), Figure 9 (public four-year, in-state; printed p. 18) and Figure 10 (private nonprofit four-year; p. 19), "
       "<https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf>.")
    ap(f"- BLS, Consumer Expenditure Survey 2018, Table R-1 (mean annual expenditures per consumer unit, people, children and earners per unit, college tuition), <{tu['ce_url']}>.")
    ap("- NLSY79 through the `nlsy_teachers.py` / `nlsy_ability.py` machinery (latent AFQT, 1979 weights); the published teacher gaps are in `nlsy_teachers.json`.")
    ap("")
    ap("## Summary")
    ap("")
    ap(f"- **BA shares (ACS, ages 25-34, Table 2 universe, home production included).** BA or more (SCHL >= 21): {s('ba_share|pooled|all|P')} of everyone "
       f"(men {s('ba_share|men|all|P')}, women {s('ba_share|women|all|P')}); {s('ba_share|pooled|nonteacher|T2')} of non-teachers (T2 convention; men "
       f"{s('ba_share|men|nonteacher|T2')}, women {s('ba_share|women|nonteacher|T2')}); {s('ba_share|pooled|teacher|P')} of the K-12 teachers of the project's "
       f"definition but {s('k12_code_ba_share|pooled|core')} of core K-12 teachers (elementary, secondary, special education; pre-K/kindergarten "
       f"{s('k12_code_ba_share|pooled|2300')}, other instructors {s('k12_code_ba_share|pooled|2340')}). Graduate degree: {s('grad_share|pooled|teacher|P')} of teachers.")
    ap(f"- **Schooling by BA status.** Mean years (capped 9-17): BA+ non-teachers {s('years_ba|pooled|nonteacher|T2', 2)}, non-BA non-teachers "
       f"{s('years_nba|pooled|nonteacher|T2', 2)}; all non-teachers {s('years|pooled|nonteacher|T2', 2)} (workbook 13.64), teachers {s('years|pooled|teacher|P', 2)} "
       f"(BA+ teachers {s('years_ba|pooled|teacher|P', 2)}).")
    ap(f"- **BA wage premium** among non-teacher market workers (log hourly wage, BA+ minus non-BA): {s('ba_premium|pooled|raw|T2')} raw "
       f"(men {s('ba_premium|men|raw|T2')}, women {s('ba_premium|women|raw|T2')}), {s('ba_premium|pooled|age|P')} controlling for age and sex, "
       f"{s('ba_premium|pooled|age_occupation|P')} also within occupation. The model has no BA-specific occupational sorting, so the raw (total return) and the within-occupation "
       "premium are alternative targets with different meanings; the choice is the M3 designer's (section 1.2).")
    ap(f"- **Teachers' pay (section 1.3).** All K-12 teachers (codes 2300-2340) earn {s('teacher_premium|men|all_nonteachers|T2', 4, True)} (men) and "
       f"{s('teacher_premium|women|all_nonteachers|T2', 4, True)} (women) against all non-teacher market workers, the existing targets +0.0557 / +0.0860 reproduced, and "
       f"{s('teacher_premium|men|ba_nonteachers|T2', 3, True)} and {s('teacher_premium|women|ba_nonteachers|T2', 3, True)} against non-teacher BA+ workers. "
       f"The {100 * (1 - sv('ba_share|pooled|teacher|P')):.1f}% of them without a BA (pre-K/kindergarten, other instructors) earn {s('teacher_ba_gap|pooled|raw|P', 2)} less than those with one, so for a model whose "
       f"teachers all have a BA the like-for-like rows are BA+ teachers: {s('teacher_premium_ba|men|all_nonteachers|T2', 3, True)} and {s('teacher_premium_ba|women|all_nonteachers|T2', 3, True)} against all non-teachers, "
       f"{s('teacher_premium_ba|men|ba_nonteachers|T2', 3, True)} and {s('teacher_premium_ba|women|ba_nonteachers|T2', 3, True)} against BA+ non-teachers (core codes {s('teacher_premium_core|men|ba_nonteachers|T2', 3, True)} "
       f"and {s('teacher_premium_core|women|ba_nonteachers|T2', 3, True)}). "
       f"Teachers with a graduate degree earn {s('grad_teacher_premium|pooled|age|P', 3, True)} more than BA-only teachers (age-controlled); the NTPS 2017-18 salary schedule gives "
       f"{_vs({'value': sal['years']['2017-18']['aggregates']['1-5 years']['log'], 'se': sal['years']['2017-18']['aggregates']['1-5 years']['se_log']}, 3, True)} "
       "(master's over bachelor's, 1-5 years of experience).")
    if WS:
        ap(f"- **Weeks rule (section 1.4).** The Table 2 wage sample (48+ weeks worked) leaves out {100 * (1 - sv('excl|pooled|teacher|in_wage_sample')):.1f}% of full-time K-12 teachers against "
           f"{100 * (1 - sv('excl|pooled|nonteacher|in_wage_sample')):.1f}% of other full-time workers, mostly teachers reporting 40-47 weeks ({100 * sv('excl|pooled|teacher|wkw3'):.1f}%), whose hourly wage is higher "
           f"for the same salary. Admitting those records raises the teacher premium against all non-teachers from {sv('teacher_premium|men|all_nonteachers|T2'):+.3f} / {sv('teacher_premium|women|all_nonteachers|T2'):+.3f} "
           f"(men / women) to {wsv('wks40', 'teacher_premium|men|all_nonteachers|T2'):+.3f} / {wsv('wks40', 'teacher_premium|women|all_nonteachers|T2'):+.3f}, and the BA+-teacher premium against BA+ non-teachers from "
           f"{sv('teacher_premium_ba|men|ba_nonteachers|T2'):+.3f} / {sv('teacher_premium_ba|women|ba_nonteachers|T2'):+.3f} to {wsv('wks40', 'teacher_premium_ba|men|ba_nonteachers|T2'):+.3f} / "
           f"{wsv('wks40', 'teacher_premium_ba|women|ba_nonteachers|T2'):+.3f}; the existing teacher-premium targets, which use the rule, are probably 0.02-0.04 lower than a sample that keeps those teachers would give.")
    ap(f"- **Ability (NLSY79).** Latent AFQT of persons with 16+ years of schooling minus the rest: {n_('all_hgc9|pooled|gap', 3, True)} SD "
       f"(employed {n_('emp_hgc9|pooled|gap', 3, True)}; men {n_('all_hgc9|men|gap', 2, True)}, women {n_('all_hgc9|women|gap', 2, True)}). P(BA | observed AFQT tercile) = "
       f"{n_('terc|all_hgc9|pooled|ba_share|1', 3)}, {n_('terc|all_hgc9|pooled|ba_share|2', 3)}, {n_('terc|all_hgc9|pooled|ba_share|3', 3)} at the NLSY79 BA share ({nv('all_hgc9|pooled|ba_share'):.3f}); "
       f"{nv('terc_acs|all_hgc9|pooled|ba_share|1'):.3f}, {nv('terc_acs|all_hgc9|pooled|ba_share|2'):.3f}, {nv('terc_acs|all_hgc9|pooled|ba_share|3'):.3f} at the ACS share ({sv('ba_share|pooled|all|P'):.3f}); "
       f"the latent logit slope is {n_('logit|all_hgc9|pooled|slope_latent_eiv', 2)} per latent SD (observed scale {n_('logit|all_hgc9|pooled|slope_obs', 2)}). The existing teacher gaps are "
       f"confirmed to six digits: {n_('teacher|pooled|nt|gap', 3, True)} against employed non-teachers and {n_('teacher|pooled|col|gap', 3, True)} against non-teacher college graduates.")
    ap(f"- **Cost of a BA (2017-18).** Published tuition and fees ${mu['sticker_tf_public_4yr_in_state']['annual']:,.0f} a year at public 4-year institutions (in-state) and "
       f"${mu['sticker_tf_all_4yr']['annual']:,.0f} at all 4-year institutions; net of grants and tax benefits ${mu['net_tf_public_4yr_in_state']['annual']:,.0f} (public) to "
       f"${mu['net_tf_private_nonprofit_4yr']['annual']:,.0f} (private nonprofit). A one-off cost paid by one person is measured against consumption *per adult* (${ex_ad:,.0f}: the CE's "
       f"${tu['ce_total_expenditures']:,.0f} per consumer unit over {adults:.1f} adults), over the working life. Central F_c = {fc['value']:.4f} of C-bar (net tuition and fees, "
       f"undiscounted); {fc['alt_annuitized_3pct_ii_b']:.4f} if annuitized at 3%; an aggregate-flow check gives {fc['alt_aggregate_flow_check']:.4f}. The first release's per-consumer-unit "
       f"figures ({fc['alt_per_cu_ii_a']:.3f} undiscounted, {fc['alt_per_cu_ii_b']:.3f} at 3%) understate F_c by the {adults:.1f} adults per unit (section 4).")
    ap("")
    ap(f"*Changes from the first release (2026-10-05 01:14).* (1) tuition is measured against expenditures per adult, not per consumer unit (a factor 1.9), and the undiscounted figure is central (F_c 0.0189 -> {fc['value']:.4f}); "
       "(2) the teacher comparisons are extended to BA+ teachers and the core codes, which are the like-for-like counterparts of a model whose teachers all have a BA; "
       "(3) the wage-sample weeks rule and its effect on the teacher premia are quantified; (4) the BA-premium choice (raw against within occupation) is left to the M3 designer; "
       "(5) the ability gradient is reported on the observed and the latent scale and at the ACS BA share; (6) the TOML carries `model_moment` names, `choice_group`s, scales for every target, "
       "correct universe sizes and URLs for every source. No ACS, NLSY or NCES statistic of the first release changed; only the tuition ratios did.")
    ap("")
    # ------------------------------------------------------------------ 1. ACS
    ap("## 1. ACS PUMS 2009-13: education, occupation and wages at 25-34")
    ap("")
    ap(f"**Sample.** `acs_occupations.build_sample(rules=\"workbook\")` on the 2013 5-year PUMS (data years 2009-13): ages 25-34, SCHL >= 12 (at least grade 9), armed forces and unemployed "
       f"persons with 15+ usual hours dropped, 50 states and DC, group quarters in: {meta['acs_persons']:,} persons ({cnt['persons']['men']:,} men, {cnt['persons']['women']:,} women), "
       f"{P['acs']['n_wage_sample']:,} in the wage sample. Home production is out of the labor force or under 15 hours; 15-29 hours are split 0.5/0.5 between home production and the "
       "occupation, so every person has total weight 1 (Table 2 counts are sums of these unweighted weights). **K-12 teachers** are Census codes 2300-2340. **BA+** is SCHL >= 21 "
       "(bachelor's, master's, professional or doctoral degree); **graduate degree** SCHL >= 22. **Years of schooling** are the workbook's mapping (grade 9 = 9, high school or GED 12, "
       "some college 13, associate 14, bachelor 16, graduate 19) capped at 9-17, so a graduate degree counts 17. **Wages:** the wage sample of Table 2 (employed, 30+ hours, 48+ weeks, "
       "income >= $1,000 in 2010 dollars), log hourly wage in 1999 dollars (annual wage and business income over usual hours times weeks), PWGTP-weighted.")
    ap("")
    wsm = cnt["wage_sample"]
    own = {kk: wsm["men"][kk] / (wsm["men"][kk] + wsm["women"][kk]) for kk in ("teachers_ba", "teachers_core")}
    dpool = max(abs((own[kk] - sv("teacher_wage_weight|men")) * (sv(f"{fam}|men|{comp}|T2") - sv(f"{fam}|women|{comp}|T2")))
                for fam, kk in (("teacher_premium_ba", "teachers_ba"), ("teacher_premium_core", "teachers_core")) for comp in ("all_nonteachers", "ba_nonteachers"))
    ap("**Two aggregation conventions.** *P*: every row weighted by PWGTP (the plain population mean). *T2*: the convention of the existing premium targets "
       "(`teacher_premium_m/f`): occupation x gender cell means (PWGTP-weighted) combined with the *unweighted* Table 2 counts of the population in question "
       "(for a BA subset, the counts of BA holders in each cell), because the model's occupation shares match those counts. They differ by up to 0.01 (K-12 teachers are "
       "1.94% of men by counts and 1.76% PWGTP-weighted). The pooled teacher premia average the two gender premia with the K-12 wage-sample weights 0.2444 / 0.7556 "
       f"(the `teacher_p9010` convention; used for every teacher group: the men's share of the BA+ and core-code teacher records, {own['teachers_ba']:.3f} and {own['teachers_core']:.3f}, would change the pooled figures by at most {dpool:.4f}), "
       "which removes the gender-mix effect (the unadjusted pooled gap is "
       f"{s('teacher_premium|pooled|all_nonteachers|P_nogenderadj', 3, True)} against all non-teachers and {s('teacher_premium|pooled|ba_nonteachers|P_nogenderadj', 3, True)} against BA+). "
       "SEs hold the Table 2 counts fixed.")
    ap("")
    ap("### 1.1 BA shares and years of schooling")
    ap("")
    rows = [
        ["BA+ share, all (P)"] + cols("ba_share|{g}|all|P"),
        ["BA+ share, non-teachers (T2)"] + cols("ba_share|{g}|nonteacher|T2"),
        ["BA+ share, non-teachers (P)"] + cols("ba_share|{g}|nonteacher|P"),
        ["BA+ share, non-teacher market workers, home production out (T2)"] + cols("ba_share|{g}|market_nonteacher|T2"),
        ["BA+ share, K-12 teachers (P)"] + cols("ba_share|{g}|teacher|P"),
        ["Graduate degree, share of K-12 teachers"] + cols("grad_share|{g}|teacher|P"),
        ["Graduate degree, share of BA+ K-12 teachers"] + cols("grad_share_of_ba|{g}|teacher|P"),
        ["Teachers' share of all 25-34 (T2 counts)"] + cols("teacher_share|{g}|T2", 4),
        ["Mean years, non-teachers (T2)"] + cols("years|{g}|nonteacher|T2", 2),
        ["Mean years, non-teachers (P)"] + cols("years|{g}|nonteacher|P", 2),
        ["Mean years, K-12 teachers (P)"] + cols("years|{g}|teacher|P", 2),
        ["Mean years, BA+ non-teachers (T2)"] + cols("years_ba|{g}|nonteacher|T2", 2),
        ["Mean years, non-BA non-teachers (T2)"] + cols("years_nba|{g}|nonteacher|T2", 2),
        ["Mean years, BA+ teachers (P)"] + cols("years_ba|{g}|teacher|P", 2),
        ["Mean years, non-BA teachers (P)"] + cols("years_nba|{g}|teacher|P", 2),
        ["Teachers' years over non-teachers' (T2)"] + cols("s_ratio|{g}|T2", 3),
    ]
    ap(_md(["", "Men", "Women", "Pooled"], rows))
    ap("")
    ap("Dividing years by 25 gives the model's time shares: non-BA non-teachers "
       f"{S['years_nba|pooled|nonteacher|T2']['value'] / 25:.3f}, BA+ non-teachers {S['years_ba|pooled|nonteacher|T2']['value'] / 25:.3f}, BA+ teachers {S['years_ba|pooled|teacher|P']['value'] / 25:.3f} "
       f"(a BA is 16/25 = 0.64; graduate degrees, 17 years after the cap, are "
       f"{100 * (S['years_ba|pooled|nonteacher|T2']['value'] - 16):.0f}% of BA+ non-teachers and {100 * S['grad_share_of_ba|pooled|teacher|P']['value']:.0f}% of BA+ teachers).")
    ap("")
    ap("### 1.2 The BA wage premium, and which premium to match")
    ap("")
    rows = [
        ["BA+ minus non-BA, non-teacher market workers: raw (T2)"] + cols("ba_premium|{g}|raw|T2"),
        ["  raw (P)"] + cols("ba_premium|{g}|raw|P"),
        ["  controlling for single-year age (pooled: and sex) (P)"] + cols("ba_premium|{g}|age|P"),
        ["  also occupation x gender fixed effects (P)"] + cols("ba_premium|{g}|age_occupation|P"),
    ]
    ap(_md(["", "Men", "Women", "Pooled"], rows))
    ap("")
    ap(f"Holding occupation x gender fixed lowers the age-controlled premium from {sv('ba_premium|pooled|age|P'):.3f} to {sv('ba_premium|pooled|age_occupation|P'):.3f}: about "
       f"{100 * (1 - sv('ba_premium|pooled|age_occupation|P') / sv('ba_premium|pooled|age|P')):.0f}% of it is BA holders' sorting into better-paid occupations. "
       "**Which premium to match.** The model has one credential premium rho_c in every non-teaching occupation, and its non-teaching occupation choice does not depend on d, so the model's "
       "raw d = 1 minus d = 0 premium among non-teaching market workers equals its within-occupation premium (plus selection on z) by construction. Matching the data's *raw* premium "
       f"({sv('ba_premium|pooled|raw|T2'):.3f}) assigns the between-occupation sorting (about {sv('ba_premium|pooled|raw|T2') - sv('ba_premium|pooled|age_occupation|P'):.2f} log points) to rho_c: it is the "
       "total return to the college decision, which is what college-going responds to, but it overstates the wage structure within occupations. Matching the *within-occupation* premium "
       f"({sv('ba_premium|pooled|age_occupation|P'):.3f}) fits that wage structure but leaves the model's total return to college below the data's, which could understate the incentive to attend. "
       "Both are in the TOML with the same `choice_group`; a criterion should contain only one, and which one is primary is the M3 designer's choice.")
    ap("")
    ap("### 1.3 Teachers' pay: which teachers, which comparison")
    ap("")
    ap("The K-12 group of the existing targets (Census codes 2300-2340) includes pre-K/kindergarten teachers and other teachers and instructors, "
       f"{100 * (1 - sv('ba_share|pooled|teacher|P')):.1f}% of the group in all lack a BA (section 1.5). A model in which teaching requires a BA has only d = 1 teachers, so the premium "
       "of the model's teachers has several data counterparts (T2 convention, log hourly wage; the first two rows are the existing figures):")
    ap("")
    grp = [("all K-12 teachers (2300-2340)", "teacher_premium", "all non-teacher workers"),
           ("all K-12 teachers (2300-2340)", "teacher_premium", "BA+ non-teacher workers"),
           ("BA+ K-12 teachers", "teacher_premium_ba", "all non-teacher workers"),
           ("BA+ K-12 teachers", "teacher_premium_ba", "BA+ non-teacher workers"),
           ("core-code teachers (2310-2330)", "teacher_premium_core", "all non-teacher workers"),
           ("core-code teachers (2310-2330)", "teacher_premium_core", "BA+ non-teacher workers")]
    rows = []
    for i, (tl, fam, cl) in enumerate(grp):
        comp = "all_nonteachers" if i % 2 == 0 else "ba_nonteachers"
        rows.append([tl, cl] + [s(f"{fam}|{g}|{comp}|T2", 4, True) for g in ("men", "women", "pooled")]
                    + [s(f"{fam}|{g}|{comp}|P", 3, True) for g in ("men", "women", "pooled")] + [s(f"{fam}|pooled|{comp}|age", 3, True)])
    ap(_md(["Teachers", "Compared with", "Men (T2)", "Women (T2)", "Pooled (T2)", "Men (P)", "Women (P)", "Pooled (P)", "Pooled, age-controlled (P)"], rows, left=(1,)))
    ap("")
    ap(f"Wage-sample record-count weights for the first row give {s('teacher_premium|men|all_nonteachers|wsw', 4, True)} / {s('teacher_premium|women|all_nonteachers|wsw', 4, True)} "
       "(the calibration log's check, +0.042 / +0.070). Which row is the counterpart of 'teachers versus d = 1 non-teachers' depends on the model: "
       "(a) if every teacher is forced to d = 1, the model's teachers are the BA+ teachers, and the counterpart is **BA+ teachers against BA+ non-teachers** (fourth row), the wage "
       "analogue of NLSY comparison C (teachers with 16+ years against non-teachers with 16+ years, section 3), with BA+ teachers against all non-teachers (third row) as the counterpart of the "
       "existing teacher_premium targets; core-code teachers (last two rows) are a robustness variant that drops the non-BA codes; (b) if the model also has a non-BA teaching track, "
       "all K-12 teachers against d = 1 non-teachers (second row, the figure of the first release) is the counterpart of the all-teacher comparison and the fourth row that of "
       f"the d = 1 teachers. Moving from the second to the fourth row raises teachers' relative pay by {sv('teacher_premium_ba|pooled|ba_nonteachers|T2') - sv('teacher_premium|pooled|ba_nonteachers|T2'):.3f} "
       f"(pooled), {100 * (sv('teacher_premium_ba|pooled|ba_nonteachers|T2') - sv('teacher_premium|pooled|ba_nonteachers|T2')) / abs(sv('teacher_premium|pooled|ba_nonteachers|T2')):.0f}% of the shortfall against "
       f"BA+ workers (women {sv('teacher_premium_ba|women|ba_nonteachers|T2') - sv('teacher_premium|women|ba_nonteachers|T2'):.3f}, men {sv('teacher_premium_ba|men|ba_nonteachers|T2') - sv('teacher_premium|men|ba_nonteachers|T2'):.3f}): "
       "that part is the low pay of the non-BA teachers, and it matters most for women, who are most of the non-BA teachers.")
    ap("")
    rows = [
        ["BA+ minus non-BA K-12 teachers: raw (P)"] + cols("teacher_ba_gap|{g}|raw|P", 3, True),
        ["  controlling for single-year age (pooled: and sex)"] + cols("teacher_ba_gap|{g}|age|P", 3, True),
        ["  also Census code (2300-2340; pooled: code x gender)"] + cols("teacher_ba_gap|{g}|age_code|P", 3, True),
        ["Teachers with a graduate degree minus BA-only teachers: raw"] + cols("grad_teacher_premium|{g}|raw|P", 3, True),
        ["  controlling for single-year age (pooled: and sex)"] + cols("grad_teacher_premium|{g}|age|P", 3, True),
    ]
    ap(_md(["", "Men", "Women", "Pooled"], rows))
    ap("")
    cw_ = cnt["wage_sample"]
    ap(f"Records in the wage sample: non-teachers {cw_['pooled']['nonteachers']:,} ({cw_['pooled']['nonteachers_ba']:,} BA+), K-12 teachers {cw_['pooled']['teachers']:,} "
       f"({cw_['pooled']['teachers_ba']:,} BA+, {cw_['pooled']['teachers_nonba']:,} without a BA, {cw_['pooled']['teachers_core']:,} in the core codes; {cw_['pooled']['teachers_grad']:,} with a graduate degree, "
       f"{cw_['pooled']['teachers_ba_only']:,} with exactly a BA). Mean log wages (T2): teachers {s('mean_lw_teacher|men', 3)} (men) and {s('mean_lw_teacher|women', 3)} (women); BA+ teachers "
       f"{s('mean_lw_teacher_ba|men', 3)} and {s('mean_lw_teacher_ba|women', 3)}; all non-teachers {s('mean_lw_all_nonteachers|men|T2', 3)} and {s('mean_lw_all_nonteachers|women|T2', 3)}; "
       f"BA+ non-teachers {s('mean_lw_ba_nonteachers|men|T2', 3)} and {s('mean_lw_ba_nonteachers|women|T2', 3)}. The BA-minus-non-BA teacher gap is the data counterpart of the pay penalty "
       "of a non-BA teaching track (rho_T0 in the model's terms); the graduate-degree rows are the ACS counterpart of the salary-schedule premium of section 2 (different objects: hourly wages at ages 25-34, "
       "including professional and doctoral degrees and district composition, against base salary at equal experience).")
    ap("")
    ap("### 1.4 Weeks worked and the wage sample")
    ap("")
    if WS:
        cf = cnt["full_time_workers"]["pooled"]
        ap("The Table 2 wage sample keeps full-time workers (30+ usual hours) who report 48 or more weeks worked in the past year (bracket 48-49 or 50-52 weeks) and at least $1,000 of income. "
           "The weeks bracket is selective for teachers (plausibly because many are on school-year contracts): "
           f"of the {cf['teachers']:,} full-time K-12 teachers in the sample, {100 * (1 - sv('excl|pooled|teacher|in_wage_sample')):.1f}% are outside the wage sample, against "
           f"{100 * (1 - sv('excl|pooled|nonteacher|in_wage_sample')):.1f}% of the {cf['nonteachers']:,} other full-time workers (PWGTP shares of full-time workers):")
        ap("")
        rows = []
        for lab, k in (("in the wage sample", "in_wage_sample"), ("48+ weeks but outside (income floor)", "wk48_outside"), ("40-47 weeks (WKW 3)", "wkw3"),
                       ("27-39 weeks (WKW 4)", "wkw4"), ("14-26 weeks (WKW 5)", "wkw5"), ("13 or fewer weeks (WKW 6)", "wkw6"), ("weeks missing", "wkw_missing")):
            rows.append([lab] + [s(f"excl|pooled|{pop}|{k}", 3) for pop in ("teacher", "teacher_ba", "nonteacher", "nonteacher_ba")])
        ap(_md(["Full-time workers", "K-12 teachers", "BA+ K-12 teachers", "Non-teachers", "BA+ non-teachers"], rows))
        ap("")
        ap("Teachers who report fewer weeks (40-47, a school-year contract, often with the summer unpaid or not counted) have a higher implied hourly wage for the same annual salary, so dropping them "
           "removes relatively well-paid teachers from the hourly sample and the teacher premia of the Table 2 rule are probably biased down. Sensitivity: the wage sample also admits the 40-47 week bracket "
           "(weeks set to the midpoint 43.5), for teachers and non-teachers alike, or also the 27-39 bracket (weeks 33); SEs are recomputed with the replicate weights.")
        ap("")
        rows = []
        for lab, fam, comp in (("All K-12 teachers minus all non-teachers (existing targets)", "teacher_premium", "all_nonteachers"),
                               ("All K-12 teachers minus BA+ non-teachers", "teacher_premium", "ba_nonteachers"),
                               ("BA+ teachers minus all non-teachers", "teacher_premium_ba", "all_nonteachers"),
                               ("BA+ teachers minus BA+ non-teachers", "teacher_premium_ba", "ba_nonteachers")):
            for g in ("men", "women", "pooled"):
                k = f"{fam}|{g}|{comp}|T2"
                rows.append([f"{lab}: {g}", s(k, 4, True), wsc("wks40", k, 4, True), wsc("wks27", k, 4, True)])
        for g in ("men", "women", "pooled"):
            k = f"ba_premium|{g}|raw|T2"
            rows.append([f"BA+ minus non-BA non-teachers, raw: {g}", s(k, 4), wsc("wks40", k, 4), wsc("wks27", k, 4)])
        ap(_md(["Premium (log points, T2)", "48+ weeks (Table 2 rule)", "also 40-47 weeks", "also 27-47 weeks"], rows))
        ap("")
        n40, n27 = WS["wks40"]["n_wage_sample"], WS["wks27"]["n_wage_sample"]
        d40 = [WS["wks40"]["stats"][f"{fam}|{g}|{comp}|T2"]["value"] - sv(f"{fam}|{g}|{comp}|T2") for fam in ("teacher_premium", "teacher_premium_ba")
               for comp in ("all_nonteachers", "ba_nonteachers") for g in ("men", "women", "pooled")]
        d27 = [WS["wks27"]["stats"][f"{fam}|{g}|{comp}|T2"]["value"] - sv(f"{fam}|{g}|{comp}|T2") for fam in ("teacher_premium", "teacher_premium_ba")
               for comp in ("all_nonteachers", "ba_nonteachers") for g in ("men", "women", "pooled")]
        ap(f"The wage sample grows from {P['acs']['n_wage_sample']:,} to {n40:,} and {n27:,} records. Admitting the 40-47 week records raises every teacher premium, by {min(d40):.3f} to {max(d40):.3f} "
           f"({min(d27):.3f} to {max(d27):.3f} with the 27-39 bracket as well): well above the sampling SEs (0.003-0.005) and within the 0.05 modelling scale; the BA premium is unchanged "
           "(raw) or moves by less than 0.003. The sign is the one argued above (the 48+ week rule drops teachers whose implied hourly wage is high), and the teacher premium targets, "
           "which use the rule, are probably 0.02-0.04 lower than a sample that keeps those teachers would give (the hours of a 40-47 week worker and the weeks within the bracket are measured with error, "
           "so the size is indicative). The rule is inherited from the workbook (HHJK) and the existing "
           "targets use it, so the main tables keep it; both variants are in the JSON and as `alt_weeks_40plus` / `alt_weeks_27plus` in the TOML.")
        ap("")
    ap("### 1.5 K-12 teachers by Census occupation code")
    ap("")
    rows = []
    for code, lab in list(K12_CODE_LABEL.items()) + [("core", "core K-12 (2310 + 2320 + 2330)")]:
        lab = f"{code} {lab}" if code != "core" else lab
        rows.append([lab, s(f"k12_code_share|pooled|{code}"), s(f"k12_code_ba_share|men|{code}"), s(f"k12_code_ba_share|women|{code}"), s(f"k12_code_ba_share|pooled|{code}")])
    ap(_md(["Code", "Share of K-12 teachers (pooled)", "BA+ men", "BA+ women", "BA+ pooled"], rows))
    ap("")
    ap(f"Licensure binds for the core K-12 codes ({100 * S['k12_code_ba_share|pooled|core']['value']:.0f}% BA+). The project's K-12 group also contains pre-K/kindergarten teachers and 'other teachers and instructors' "
       f"({100 * (1 - S['k12_code_share|pooled|core']['value']):.0f}% of teachers by PWGTP, {100 * S['k12_code_ba_share|pooled|2300']['value']:.0f}% and {100 * S['k12_code_ba_share|pooled|2340']['value']:.0f}% BA+), "
       f"so only {100 * S['ba_share|pooled|teacher|P']['value']:.0f}% of the group has a BA. A model that imposes s >= 16/25 on every teacher is matched to teaching shares that include these codes.")
    ap("")
    ap("### 1.6 Reproduction of the published block")
    ap("")
    ap(f"The sample reproduces `acs_occupations.json` exactly: maximum absolute differences are {rep['max_abs_diff_table2_counts']:.1f} persons in the 42 Table 2 counts, "
       f"{rep['max_abs_diff_wage_sample_counts']:.1f} wage-sample records, {rep['max_abs_diff_mean_log_wage_by_cell']:.1e} in cell mean log wages and "
       f"{rep['max_abs_diff_mean_years_by_cell']:.1e} years in cell mean schooling. The T2 teacher premia are {rep['teacher_premium_men_T2']:.4f} / {rep['teacher_premium_women_T2']:.4f} "
       f"(coded targets {rep['teacher_premium_men_target']:.4f} / {rep['teacher_premium_women_target']:.4f}); with wage-sample weights {rep['teacher_premium_men_wsw']:.4f} / {rep['teacher_premium_women_wsw']:.4f}. "
       f"Years of schooling: K-12 teachers {rep['years_k12_men']:.3f} (men) and {rep['years_k12_women']:.3f} (women), non-teachers {rep['years_nonteachers_pooled']:.3f} pooled, "
       f"against 16.04 and 13.64 in the Mincer workbook (the calibration's constants 16.0437 and 13.6419 come from the workbook, which cannot be matched exactly).")
    ap("")
    # ------------------------------------------------------------------ 2. NCES salary
    ap("## 2. Teacher salary schedules (NCES Digest Table 211.20)")
    ap("")
    y = sal["years"]["2017-18"]
    ap(f"Source: {sal['title']} (<{sal['source']}>), average base salary (current dollars) of full-time public school teachers by highest degree and years of teaching experience; "
       "the 2017-18 rows are the National Teacher and Principal Survey (NTPS) and are identical in the 2020 and 2022 editions "
       f"(max absolute difference {sal['edition_check_2017_18_vs_digest_2020'].get('max_abs_diff', float('nan')):.0f} over {sal['edition_check_2017_18_vs_digest_2020'].get('n_compared', 0)} cells). "
       "Base salary excludes supplemental contracts and nonschool income; experience counts all school years, full or part time. Standard errors are the table's; the SE of the log ratio treats "
       "the two means as independent (the two groups are disjoint teachers).")
    ap("")
    rows = []
    for lab in ["1 year or less", "2 years", "3 years", "4 years", "5 years", "6 to 9 years"]:
        r = y["rows"][lab]
        rows.append([lab, f"{r['ba']:,.0f} ({r['se_ba']:,.0f})", f"{r['ma']:,.0f} ({r['se_ma']:,.0f})", f"{r['dollars']:,.0f}", f"{r['log']:.3f} ({r['se_log']:.3f})"])
    for lab in ("1-3 years", "1-5 years", "1-9 years"):
        a = y["aggregates"][lab]
        rows.append([f"**{lab}** (weighted by teachers)", f"{a['ba']:,.0f}", f"{a['ma']:,.0f}", f"{a['dollars']:,.0f}", f"**{a['log']:.3f}** ({a['se_log']:.3f})"])
    tr = y["total_row"]
    rows.append(["All experience (composition differs)", f"{tr['ba']:,.0f}", f"{tr['ma']:,.0f}", f"{tr['ma'] - tr['ba']:,.0f}", f"{np.log(tr['ma'] / tr['ba']):.3f}"])
    ap(_md(["2017-18, experience", "Bachelor's", "Master's", "Master's - bachelor's ($)", "log(M/B) (SE)"], rows))
    ap("")
    rows = []
    survey = {"2011-12": "2011-12 (SASS)", "2015-16": "2015-16 (NTPS)", "2017-18": "2017-18 (NTPS)", "2020-21": "2020-21 (NTPS)"}
    for yr in ("2011-12", "2015-16", "2017-18", "2020-21"):
        if yr in sal["years"]:
            a1, a2 = sal["years"][yr]["aggregates"]["1-5 years"], sal["years"][yr]["rows"]["1 year or less"]
            rows.append([survey[yr], f"{a1['log']:.3f} ({a1['se_log']:.3f})", f"{a1['dollars']:,.0f}", f"{a2['log']:.3f} ({a2['se_log']:.3f})"])
    ap(_md(["Survey year (survey)", "log(M/B), 1-5 years (SE)", "dollars, 1-5 years", "log(M/B), first year (SE)"], rows))
    ap("")
    lo_hi = [sal["years"][yr]["aggregates"]["1-5 years"]["log"] for yr in sal["years"]]
    ap(f"The master's premium at the same experience is {y['rows']['1 year or less']['log']:.2f} in the first year and {y['aggregates']['1-5 years']['log']:.2f} over years 1-5, and the 1-5 year figure is stable across "
       f"surveys ({min(lo_hi):.2f}-{max(lo_hi):.2f}); the 2011-12 row is from the Schools and Staffing Survey (SASS), the others from the NTPS. It is below the ACS hourly-wage premium of teachers aged 25-34 with a graduate degree over BA-only teachers "
       f"({s('grad_teacher_premium|pooled|age|P')} age-controlled, {s('grad_teacher_premium|pooled|raw|P')} raw), which also contains professional and doctoral degrees, district composition and hours; "
       "at a given age the graduate teachers also have fewer years of experience, which works the other way.")
    ap("")
    # ------------------------------------------------------------------ 3. NLSY79
    ap("## 3. NLSY79: ability, BA status and teachers")
    ap("")
    m3 = P["nlsy"]["meta"]
    R_ = nv("R_comp_pop")
    ap(f"Block A (population sample, N = {m3['n_block_a']:,}), persons observed at 25-34 (survey years 1982-1998; {m3['n_observed_25_34']:,} persons). Ability is `afqt_c`; a gap in means divided by "
       f"sqrt(R_comp_pop) = {np.sqrt(R_):.3f} is a gap in latent SDs. BA+ = highest grade completed >= 16 at some round at 25-34 (the extract has no degree variable; "
       "the convention of `nlsy_teachers.py`). 1979 weights, one observation per person, household-cluster bootstrap with the draws of `nlsy_teachers.py`. "
       f"Failed draws: {m3['failed_draws']}.")
    ap("")
    rows = []
    for sk in ("all_hgc9", "all", "emp_hgc9", "emp", "emp_nonteacher_hgc9", "emp_nonteacher"):
        rows.append([BA_SAMPLES[sk], f"{NC_[sk]['pooled']['persons']:,}", n_(f"{sk}|pooled|ba_share", 3), n_(f"{sk}|pooled|gap", 3, True), n_(f"{sk}|men|gap", 3, True), n_(f"{sk}|women|gap", 3, True)])
    ap(_md(["Sample", "Persons", "BA+ share", "Gap BA+ minus non-BA (SD), pooled", "men", "women"], rows))
    ap("")
    ap(f"The primary sample (first row) is the NLSY counterpart of the ACS Table 2 universe. The gap is large (about 1.2 latent SDs) and insensitive to the sample; at the ACS BA share "
       f"({m3['share_target']['pooled']:.3f} against {n_('all_hgc9|pooled|ba_share', 3)} here) a logit level shift of the same gradient gives {n_('logit|all_hgc9|pooled|gap_at_acs_share', 3, True)}, "
       "so the cohort difference in BA attainment does not move it.")
    ap("")
    ap(f"**BA status by AFQT tercile and the logit slope.** The observed score has reliability R = {R_:.3f}, so statistics of its terciles and of the logit of BA on the observed score are attenuated "
       "relative to the latent ability that the model's z represents. Three forms are reported. *Observed* (the data as they are): terciles of the weighted sample, pooled over gender "
       f"(cut points {nv('terc|all_hgc9|cut1'):.3f} and {nv('terc|all_hgc9|cut2'):.3f} observed SDs), at the NLSY79 BA share. *At the ACS share*: the fitted observed-score logit with its level shifted so that the BA share is the ACS's "
       f"({m3['share_target']['pooled']:.3f}; the NLSY79 cohort has fewer BAs). *Latent*: the logit of BA on the latent factor, estimated by maximum likelihood with a normal factor and reliability R "
       "(Gauss-Hermite integration of the logistic over the factor given the score; assumes a logistic P(BA | factor) and a normal factor), with P(BA | tercile of the factor) at both shares.")
    ap("")
    rows = []
    for q, lab in ((1, "T1 (lowest)"), (2, "T2"), (3, "T3 (highest)")):
        rows.append([lab] + [n_(f"terc|all_hgc9|{g}|ba_share|{q}", 3) for g in ("pooled", "men", "women")] + [n_(f"terc_fit|all_hgc9|pooled|ba_share|{q}", 3)]
                    + [n_(f"terc_acs|all_hgc9|pooled|ba_share|{q}", 3), n_(f"terc_latent|all_hgc9|pooled|ba_share|{q}", 3), n_(f"terc_latent_acs|all_hgc9|pooled|ba_share|{q}", 3)])
    rows.append(["T3 - T1"] + [n_(f"terc|all_hgc9|{g}|ba_share|diff", 3, True) for g in ("pooled", "men", "women")] + [""]
                + [n_("terc_acs|all_hgc9|pooled|ba_share|diff", 3, True), n_("terc_latent|all_hgc9|pooled|ba_share|diff", 3, True), n_("terc_latent_acs|all_hgc9|pooled|ba_share|diff", 3, True)])
    ap(_md(["P(BA+ | tercile)", "Observed, pooled", "men", "women", "Logit fit to the observed, pooled", "Observed logit, at ACS share", "Latent terciles", "Latent, at ACS share"], rows))
    ap("")
    rows = [["Slope per observed AFQT SD (logit, observed score)"] + [n_(f"logit|all_hgc9|{g}|slope_obs", 2) for g in ("pooled", "men", "women")],
            ["Slope per latent SD, first-order correction (observed / sqrt(R))"] + [n_(f"logit|all_hgc9|{g}|slope_latent_approx", 2) for g in ("pooled", "men", "women")],
            ["Slope per latent SD, errors-in-variables logit (MLE)"] + [n_(f"logit|all_hgc9|{g}|slope_latent_eiv", 2) for g in ("pooled", "men", "women")]]
    ap(_md(["Logit of BA+ on ability", "Pooled", "Men", "Women"], rows))
    ap("")
    ap(f"The first-order correction understates the latent slope ({nv('logit|all_hgc9|pooled|slope_latent_approx'):.2f} against {nv('logit|all_hgc9|pooled|slope_latent_eiv'):.2f}, "
       f"{100 * (nv('logit|all_hgc9|pooled|slope_latent_eiv') / nv('logit|all_hgc9|pooled|slope_latent_approx') - 1):.0f}% more) because measurement error also flattens the "
       "non-linearity of the logistic (the attenuation of a logit slope exceeds the linear factor sqrt(R)). The latent logit is a parametric deconvolution: it implies a latent gap of "
       f"{nv('logit|all_hgc9|pooled|gap_latent_eiv'):.3f} against {nv('all_hgc9|pooled|gap'):.3f} measured directly (men {nv('logit|all_hgc9|men|gap_latent_eiv'):.3f} against {nv('all_hgc9|men|gap'):.3f}, "
       f"women {nv('logit|all_hgc9|women|gap_latent_eiv'):.3f} against {nv('all_hgc9|women|gap'):.3f}), and the observed-score logit fits the observed terciles to within "
       f"{max(abs(nv(f'terc_fit|all_hgc9|pooled|ba_share|{q}') - nv(f'terc|all_hgc9|pooled|ba_share|{q}')) for q in (1, 2, 3)):.3f}. "
       "**Model counterparts.** A model without score noise should be compared with the *latent* terciles and slope. To compare with the observed-scale statistics instead, add noise to the model's standardized "
       f"log z (theta): x = sqrt(R) theta + sqrt(1 - R) u with u standard normal and R = {R_:.3f}, and form terciles of x (over adults) and the logit slope of P(d = 1) on x. "
       "The first columns are at the NLSY79 BA share; the model's BA share will be near the ACS's, so the 'at ACS share' columns (a level shift of the same logit) are the relevant levels.")
    ap("")
    ap("**Teacher gaps recomputed (confirmation).** Same Frame, spec `c80-emp-all-maj` (employed, 1980 Census codes, all K-12 codes, majority of rounds), same bootstrap draws:")
    ap("")
    pub = P["nlsy"].get("published_teacher_gaps", {})
    rows = []
    for cl, lab in (("nt", "A. teachers minus employed non-teachers"), ("col", "B. teachers minus non-teacher persons with 16+ years"), ("cc", "C. teachers with 16+ years minus non-teachers with 16+ years")):
        rows.append([lab, n_(f"teacher|pooled|{cl}|gap", 3, True), f"{pub[f'pooled|{cl}']['value']:+.3f} ({pub[f'pooled|{cl}']['se']:.4f})" if pub else "n/a",
                     n_(f"teacher|men|{cl}|gap", 3, True), n_(f"teacher|women|{cl}|gap", 3, True)])
    ap(_md(["Comparison (latent AFQT SDs)", "Here, pooled", "Published `nlsy_teachers.json`, pooled", "men", "women"], rows))
    ap("")
    ap(f"Teachers have 16+ years of schooling in {n_('teacher|pooled|col_share_t', 3)} of cases, employed non-teachers in {n_('teacher|pooled|col_share_nt', 3)}. "
       f"Point estimates and SEs equal the published ones (max abs. difference {P['nlsy']['confirmation']['max_abs_diff_value']:.1e} in values, {P['nlsy']['confirmation']['max_abs_diff_se']:.1e} in SEs).")
    ap("")
    # ------------------------------------------------------------------ 4. Cost
    ap("## 4. The cost of a BA, in model units")
    ap("")
    ap(f"Sources: published prices from NCES Digest 2018 Table 330.10 (<https://nces.ed.gov/programs/digest/d18/tables/dt18_330.10.asp>; average undergraduate tuition and required fees, dormitory and board charges of full-time students, 2017-18, "
       "FTE-weighted, in-state rates at public institutions; current = constant dollars in the base year); net prices from the College Board (*Trends in College Pricing 2018*, Figures 9 and 10, printed pp. 18-19: published "
       "price minus grant aid from all sources and federal education tax benefits, 2017-18, 2018 dollars; <https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf>) and NCES Digest 2019 Table 331.30 "
       "(<https://nces.ed.gov/programs/digest/d19/tables/dt19_331.30.asp>; cost of attendance minus grants for first-time full-time "
       "students awarded Title IV aid, which includes room, board, books and other expenses and excludes the 17% with no Title IV aid). Enrollment weights for the blend: fall 2017 undergraduates at "
       f"4-year institutions, Table 303.70 (<https://nces.ed.gov/programs/digest/d18/tables/dt18_303.70.asp>; {100 * tu['weight_public']:.1f}% public among public and private nonprofit; for-profits, 6%, excluded). "
       f"The College Board's published public in-state price "
       f"(${tu['college_board_2018']['public_4yr_in_state']['published_tuition_fees']:,}, 2018 dollars) is above NCES's (${mu['sticker_tf_public_4yr_in_state']['annual']:,.0f}): different coverage and enrollment weights; the net "
       f"prices are the College Board's own (re-read from the report at run time: {tu['college_board_check'].get('n_checked', 0)} figures checked, "
       f"{'all equal' if tu['college_board_check'].get('verified') else 'NOT verified'}). A BA costs {tu['years_of_college']} years of the annual amount (5 years in the JSON).")
    ap("")
    ap(f"**Units: two choices.** The model is a two-period OLG whose adult period stands for the whole working life; its flows were converted as annual dollars over CE mean annual expenditures "
       f"(${tu['ce_total_expenditures']:,.2f} per consumer unit, CE 2018 Table R-1, <{tu['ce_url']}>, as in `housing_costs.md`).")
    ap("")
    r_b = tu["working_years"] / tu["annuity_factor_r"]
    r_c = r_b * tu["compounding_factor_r"] / tu["years_of_college"]
    ap(f"1. *Per consumer unit or per person.* The CE unit has {tu['ce_persons_per_cu']:.1f} people, {adults:.1f} adults ({tu['ce_persons_per_cu']:.1f} minus {tu['ce_children_per_cu']:.1f} children under 18; the one-decimal "
       f"inputs leave a range of about 1.8-2.0) and {tu['ce_earners_per_cu']:.1f} earners. A housing cost is per household and scales with the household, so its ratio to expenditures per unit is unit-free "
       "(`housing_costs.md`, *Units*). Tuition is paid by one student and does not scale with household size. The model's agent is one adult (every person of the Table 2 universe, home production included, "
       f"is in the model), so the consistent denominator is expenditures **per adult**, ${tu['ce_total_expenditures']:,.0f} / {adults:.1f} = ${ex_ad:,.0f}. Per earner "
       f"(${tu['denominators']['per_earner']['expenditures']:,.0f}) is the alternative if only workers count; it gives {tu['ce_earners_per_cu'] / adults:.2f} times the per-adult F_c. The first release divided "
       f"by the per-unit figure, which understates F_c by the factor {adults:.1f} for a one-adult model household.")
    ap(f"2. *Annual or one-off.* A *recurring* annual cost such as the residence cost r_p enters as annual dollars over annual expenditures (column i). A *one-off* goods cost paid from adult consumption "
       f"(the education cost e and the moving cost m are paid this way, `spatial_teachers_model_eps0.md` section 2, *Timing*) is spread over the {tu['working_years']} working years: the consistent "
       f"counterpart is **(ii) total cost over working-life consumption**. Undiscounted, total / ({tu['working_years']} x annual consumption) (ii-a), is the counterpart in a model without an interest rate and "
       "with the other education parameters calibrated on annual aggregate flows (eta from education spending over labor income; the schooling time share s from years over 25): in a stationary economy "
       "an annual aggregate flow ratio equals the undiscounted lifetime ratio per person. Discounting the 40 years at "
       f"{100 * tu['discount_rate']:.0f}% (ii-b, total / ({tu['annuity_factor_r']:.3f} x annual consumption), the annuitized debt service as a share of annual consumption) or also compounding the four college payments to the start of "
       f"working life (ii-c, factor {tu['compounding_factor_r']:.3f} instead of 4) adds the cost of financing and gives {r_b:.2f} and {r_c:.2f} times (ii-a); it is the right figure only if the model's period "
       f"is meant to carry the interest on tuition. Column (i) charges the cost every year of adult life and overstates a four-year cost by {tu['working_years']}/{tu['years_of_college']} = {tu['working_years'] // tu['years_of_college']} undiscounted "
       f"({tu['annuity_factor_r']:.1f}/{tu['years_of_college']} = {tu['annuity_factor_r'] / tu['years_of_college']:.1f} at 3%). A single-year reading of the adult period would reverse this; it is not the model's.")
    ap("")
    rows = []
    for k, mm in mu.items():
        u = mm["model_units_per_adult"]
        rows.append([mm["label"], f"{u['annual']:,.0f}", f"{u['total']:,.0f}", f"{u['i_annual_over_ce']:.4f}", f"{u['ii_total_over_40ce']:.4f}", f"{u['ii_total_over_pv40ce_r']:.4f}", f"{u['ii_compounded_total_over_pv40ce_r']:.4f}"])
    ap(f"Per adult (denominator ${ex_ad:,.0f}):")
    ap("")
    ap(_md(["Measure (2017-18)", "Annual $", f"{tu['years_of_college']}-year total $", "(i) annual / CE", "(ii-a) total / 40 CE", "(ii-b) at 3%", "(ii-c) compounded, at 3%"], rows))
    ap("")
    rows = []
    for lab, un in (("per consumer unit (the first release)", "model_units"), ("per adult (central)", "model_units_per_adult"), ("per earner", "model_units_per_earner")):
        u = cen[un]
        rows.append([lab, f"{u['i_annual_over_ce']:.4f}", f"{u['ii_total_over_40ce']:.4f}", f"{u['ii_total_over_pv40ce_r']:.4f}", f"{u['ii_compounded_total_over_pv40ce_r']:.4f}"])
    ap("Denominator sensitivity for the central measure (net tuition and fees, enrollment-weighted):")
    ap("")
    ap(_md(["Denominator", "(i) annual / CE", "(ii-a) total / 40 CE", "(ii-b) at 3%", "(ii-c) compounded, at 3%"], rows))
    ap("")
    rows = []
    for k in ("net_tf_blend_4yr", "net_tf_public_4yr_in_state", "sticker_tf_all_4yr"):
        mm = mu[k]
        row = [mm["label"] + (" [central]" if k == "net_tf_blend_4yr" else "")]
        for un in ("model_units", "model_units_per_adult", "model_units_per_earner"):
            row += [f"{mm[un]['ii_total_over_40ce']:.4f}", f"{mm[un]['ii_total_over_pv40ce_r']:.4f}"]
        rows.append(row)
    ap("The same for the central measure and the two ends of the range, (ii-a) / (ii-b):")
    ap("")
    ap(_md(["Measure", "CU (ii-a)", "CU (ii-b)", "Adult (ii-a)", "Adult (ii-b)", "Earner (ii-a)", "Earner (ii-b)"], rows))
    ap("")
    ap(f"**Aggregate-flow cross-check (unit-free).** The CE records ${flow['college_tuition_per_cu']:,.2f} of college tuition per consumer unit in 2018 (SE {flow['college_tuition_se']:.2f}), "
       f"{100 * flow['share_of_expenditures']:.2f}% of expenditures (SE {100 * flow['share_se_approx']:.2f} points, ignoring the covariance). In a stationary economy this share equals P(d = 1) x F_c, "
       f"so with the ACS BA+ share {flow['ba_share']:.3f} it implies F_c = {flow['implied_Fc']:.4f} (SE about {flow['implied_Fc_se_approx']:.4f}) of consumption, with no per-unit or per-adult conversion and "
       "no discounting. It is larger than the central BA-only figure because it counts all postsecondary tuition households pay (graduate and professional schools, two-year colleges, non-completers) "
       f"over the BA+ share, and close to the annuitized per-adult figure ({fc['alt_annuitized_3pct_ii_b']:.4f}); it supports a range of about {fc['value']:.2f}-{fc['alt_aggregate_flow_check']:.2f} per adult.")
    ap("")
    ap(f"**Recommended.** The education goods cost is tuition and fees net of grants, not room and board (living costs are paid whether or not one attends; the model has no young-age consumption). "
       f"Central F_c = **{fc['value']:.4f}** of C-bar per adult (net tuition and fees, enrollment-weighted public/private nonprofit, four years, (ii-a), undiscounted: ${cen['model_units_per_adult']['total']:,.0f} over 40 x ${ex_ad:,.0f}); "
       f"{fc['alt_annuitized_3pct_ii_b']:.4f} annuitized at 3% (ii-b), {fc['alt_per_earner_ii_a']:.4f} / {fc['alt_per_earner_ii_b']:.4f} per earner; range {fc['range_low']:.4f} (public 4-year net, undiscounted) to "
       f"{fc['range_high']:.4f} (published tuition, all 4-year, 3%). The choice between undiscounted and annuitized is a judgment about whether the model's period carries the interest on tuition, and the level "
       "is small next to the 0.48 log-point BA wage premium: tuition is a minor part of the cost of a BA next to the time cost (the s margin) and any utility cost c_C. Note that eta is calibrated on total education "
       "spending, 6.0% of GDP in 2018 (NCES Digest 2021 Table 605.20, government and private expenditures on educational institutions at all levels, postsecondary included, "
       "<https://nces.ed.gov/programs/digest/d21/tables/dt21_605.20.asp>): if eta is left unchanged, adding F_c counts the college part of that spending twice (F_c = 0 with eta unchanged is the alternative).")
    ap("")
    # ------------------------------------------------------------------ 5. targets
    ap("## 5. Moments for the model")
    ap("")
    ap("`education_targets.toml` carries each moment with its definition, SE, a suggested criterion scale for targets, the model counterpart and the name of the model-side moment (`model_moment`). Roles follow DESIGN.md M3: *target* "
       "(BA share of non-teachers, BA wage premium in its two forms, NLSY79 AFQT gap, s_O mean), alternative forms by gender or with controls, and *validation* (schooling ratio, teacher premia for each teacher group and comparison, teacher z gaps, "
       "BA status by AFQT tercile, the latent and observed logit slopes). `ba_premium` and `ba_premium_within_occupation` (and the age-controlled form) are alternatives of one moment (`choice_group`).")
    ap("")
    rows = []
    for k, e in T.items():
        if e["role"] in ("target", "target_by_gender", "target_alt", "validation", "cost"):
            se = e.get("se")
            val = f"{e['value']:.4f}" + (f" ({se:.4f})" if se is not None and np.isfinite(se) and se >= 5e-5 else (" (<0.0001)" if se is not None and np.isfinite(se) else ""))
            c = e.get("convention", "")
            agg = ("T2" if c.startswith("T2") else "P" if c.startswith("P") else "weighted means" if c.startswith("weighted") else
                   "workbook" if c.startswith("Mincer") else "delta method" if c.startswith("delta") else "(ii-a)" if c.startswith("(ii") else
                   "EIV logit" if c.startswith("errors") else c)
            rows.append([f"`{k}`", e["role"], val, e.get("units", "").split(" (")[0].split(" =")[0], agg, e.get("model_moment", "")])
    ap(_md(["Key", "Role", "Value (SE)", "Units", "Aggregation", "Model moment"], rows, left=(1, 3, 4, 5)))
    ap("")
    # ------------------------------------------------------------------ caveats
    ap("## Caveats")
    ap("")
    ap("- **Sampling error is not the error that matters.** ACS replicate-weight SEs are 0.001-0.005 on shares and premia; the conventions (T2 against P, comparison group, wage-sample rules) move the "
       "numbers by 0.005-0.01 (the weeks rule by more, section 1.4), and the existing targets use a 0.05 modelling scale for premia. The suggested scales in the TOML follow that practice. SEs treat the Table 2 counts as fixed.")
    ap("- **Degrees against years.** ACS BA+ is a degree (SCHL >= 21); NLSY79 BA+ is 16+ completed years (no degree variable in the extract). The NLSY79 cohort (born 1957-64, observed 1982-98) had a BA share "
       f"of {n_('all_hgc9|pooled|ba_share', 3)} against {s('ba_share|pooled|all|P')} in the ACS 2009-13; the gap is not sensitive to this (section 3), the tercile shares are (hence the 'at ACS share' columns).")
    ap("- **AFQT timing and measurement.** The ASVAB was taken at ages 15-23 and partly reflects schooling already attained, so the BA gap overstates the ability difference at the schooling decision "
       "(the caveat of `nlsy_teachers.md`); latent SDs assume classical error and the same error variance in every group; terciles of the observed score flatten the latent gradient. The latent logit "
       "slope and terciles are model-based (a logistic P(BA | factor), a normal factor and the composite reliability R) and not direct data moments.")
    ap("- **Teachers.** The K-12 definition (codes 2300-2340) includes pre-K/kindergarten and other instructors; 'teaching requires a BA' holds for the core codes only. "
       "The existing teacher premia and the teaching-share targets use the whole group; the BA+ and core-code rows of section 1.3 are the like-for-like counterparts of a model whose teachers all have d = 1. "
       f"The +0.614 / -0.334 NLSY79 gaps rest on {NC_['k12_teachers']['pooled']} teachers ({NC_['k12_teachers']['men']} men, {NC_['k12_teachers']['women']} women).")
    ap("- **Wages.** Hourly wages are annual wage and business income over usual hours x weeks (weeks bracketed to 49 or 51). The 48+ week rule is selective for teachers (section 1.4): it leaves out "
       + (f"{100 * (1 - sv('excl|pooled|teacher|in_wage_sample')):.1f}% of full-time K-12 teachers against {100 * (1 - sv('excl|pooled|nonteacher|in_wage_sample')):.1f}% of other full-time workers, " if WS else "")
       + "mostly teachers who report 40-47 weeks and whose implied hourly wage is higher for the same salary, so the teacher premia are probably biased down; the sensitivity is reported, and the rule is inherited from "
       "the workbook and the existing targets. The salary-schedule premium (base salary, same experience) and the ACS hourly premium (same age) are different objects.")
    ap("- **NCES salary table.** Experience-specific SEs are the table's (survey design included); the log-ratio SEs assume independence; the aggregates weight rows by all-teacher counts because counts "
       "by degree are not published.")
    ap("- **Costs.** Published prices are averages of sticker prices (in-state at public institutions); net prices average over students including those paying nothing; the IPEDS net price covers Title IV "
       f"recipients only and includes living costs. The CE denominator is per consumer unit ({tu['ce_persons_per_cu']:.1f} people, {adults:.1f} adults, {tu['ce_earners_per_cu']:.1f} earners; one-decimal inputs, so adults per unit are known to about "
       f"+/-5%); the central F_c divides by adults, the per-earner and per-unit figures are in section 4. 2017-18 academic-year and 2018 calendar-year dollars "
       "differ by about 1-2% and are not adjusted. The 3% rate and 40 working years are the task's assumptions. BA duration is the nominal four years; completers take longer and non-completers also pay. "
       "As for the other goods costs, mapping the CE's expenditures per adult to the model's C-bar is a unit convention, not a calibration of C-bar.")
    ap("- **Untouched sample.** The NLSY97 counterparts are reserved for the sealed hold-out (DESIGN.md D4) and were not computed here.")
    ap("")
    ap("Written with help from Claude Code.")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


# =============================================================================
# 7. Confirmation of the published teacher gaps, self-test, main
# =============================================================================

def confirm_teacher_gaps(nl: dict[str, Any], path: Path) -> None:
    """Compare the recomputed teacher gaps with `nlsy_teachers.json` (values and bootstrap SEs)."""
    pub = json.loads(path.read_text())["statistics"]
    pg, dv, ds = {}, [], []
    for g in ("pooled", "men", "women"):
        for cl in ("nt", "col", "cc"):
            p_ = pub[f"{TEACHER_SPEC}|{g}|{cl}|gap"]
            r = nl["stats"][f"teacher|{g}|{cl}|gap"]
            pg[f"{g}|{cl}"] = {"value": p_["value"], "se": p_["se"]}
            dv.append(abs(r["value"] - p_["value"]))
            ds.append(abs(r["se"] - p_["se"]))
    nl["published_teacher_gaps"] = pg
    nl["confirmation"] = {"max_abs_diff_value": float(max(dv)), "max_abs_diff_se": float(max(ds)), "file": str(path.name)}


def selftest() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'ok' if cond else 'FAIL'}] {name}")

    rng = np.random.default_rng(1)
    # units
    check("annuity and compounding factors", abs(annuity(0.03) - 23.114772) < 1e-5 and annuity(0.0) == 40.0
          and abs(compound_to_entry(0.03) - 4.30913581) < 1e-7 and compound_to_entry(0.0) == 4.0)
    u = model_units(10000.0, 4, 50000.0, 0.03)
    check("model units identities", abs(u["i_annual_over_ce"] - 0.2) < 1e-12 and abs(u["ii_total_over_40ce"] - 0.02) < 1e-12
          and abs(u["ii_total_over_pv40ce_r"] - 40000 / (annuity(0.03) * 50000)) < 1e-12
          and abs(u["ii_total_over_40ce"] / u["i_annual_over_ce"] - 0.1) < 1e-12)
    check("table cell parser", num("$12,615") == 12615.0 and num("(365)") == 365.0 and np.isnan(num("‡")) and np.isnan(num("")))
    # logit and the share adjustment
    x = rng.standard_normal(200000)
    p = 1 / (1 + np.exp(-(-1.0 + 1.5 * x)))
    y = (rng.random(len(x)) < p).astype(float)
    a, b = wlogit(x, y, np.ones(len(x)))
    check("weighted logit recovers (a, b)", abs(a + 1.0) < 0.03 and abs(b - 1.5) < 0.03)
    w = np.ones(len(x))
    gap, delta, sh = share_adjusted_gap(x, w, -1.0, 1.5, float(p.mean()))
    m1, m0 = (p * x).sum() / p.sum(), ((1 - p) * x).sum() / (1 - p).sum()
    check("share adjustment is the identity at the sample share", abs(delta) < 1e-6 and abs(gap - (m1 - m0)) < 1e-6 and abs(sh - p.mean()) < 1e-9)
    gap2, d2, sh2 = share_adjusted_gap(x, w, -1.0, 1.5, 0.5)
    check("share adjustment hits the target share", abs(sh2 - 0.5) < 1e-9 and d2 > 0)
    # cell regression against person-level weighted least squares
    n = 4000
    age = rng.integers(0, NAGE, n)
    flag = (rng.random(n) < 0.4).astype(float)
    g = rng.integers(0, 2, n)
    w_ = rng.integers(5, 100, n).astype(float)
    lw = 1.0 + 0.5 * flag + 0.02 * age - 0.1 * g + rng.normal(0, 0.3, n)
    X = np.column_stack([np.ones(n), flag] + [(age == a_).astype(float) for a_ in range(1, NAGE)] + [g.astype(float)])
    beta = np.linalg.lstsq(X * np.sqrt(w_)[:, None], lw * np.sqrt(w_), rcond=None)[0]
    check("cell-level WLS equals person-level WLS (pooled)", abs(age_flag_coef(age, flag, g, lw, w_, True) - beta[1]) < 1e-9)
    sel = g == 0
    Xm = np.column_stack([np.ones(sel.sum()), flag[sel]] + [(age[sel] == a_).astype(float) for a_ in range(1, NAGE)])
    bm = np.linalg.lstsq(Xm * np.sqrt(w_[sel])[:, None], lw[sel] * np.sqrt(w_[sel]), rcond=None)[0]
    check("cell-level WLS equals person-level WLS (one gender)", abs(age_flag_coef(age[sel], flag[sel], g[sel], lw[sel], w_[sel], False) - bm[1]) < 1e-9)
    fe = rng.integers(0, 5, n)
    Xf = np.column_stack([np.ones(n), flag] + [(age == a_).astype(float) for a_ in range(1, NAGE)] + [(fe == k).astype(float) for k in range(1, 5)])
    bf = np.linalg.lstsq(Xf * np.sqrt(w_)[:, None], lw * np.sqrt(w_), rcond=None)[0]
    check("cell-level WLS with fixed effects", abs(age_flag_coef(age, flag, g, lw, w_, True, fe=fe, n_fe=5) - bf[1]) < 1e-9)
    # a synthetic sample through AcsData / acs_stats
    m = 6000
    names = [x_ for x_ in GROUPS if x_ != HOME]
    status = rng.choice(["full", "part", "home"], m, p=[0.6, 0.1, 0.3])
    w_mkt = np.where(status == "full", 1.0, np.where(status == "part", 0.5, 0.0))
    w_home = 1.0 - w_mkt
    grp = np.where(w_mkt > 0, rng.choice(names, m), None)
    grp = np.where((grp != None) & (rng.random(m) < 0.15), K12, grp)             # noqa: E711
    schl = np.where(grp == K12, rng.choice([21, 22, 20], m, p=[0.5, 0.4, 0.1]), rng.choice([15, 18, 20, 21, 22], m, p=[0.3, 0.2, 0.15, 0.25, 0.1]))
    full = status == "full"
    wkw = np.where(full, rng.choice([1, 2, 3, 4, 5, 6], m, p=[0.55, 0.15, 0.12, 0.08, 0.06, 0.04]), np.nan)
    wkhp = np.where(full, rng.integers(30, 60, m), np.where(status == "part", rng.integers(15, 30, m), 0)).astype(float)
    wagp = np.where(full, rng.uniform(300, 90000, m), 0.0)                       # some below the $1,000 floor
    data_year = rng.choice([2009, 2010, 2011, 2012, 2013], m)
    weeks = pd.Series(wkw).map(WKW_WEEKS)
    cpi_y = pd.Series(data_year).map(ao.CPI_U)
    with np.errstate(divide="ignore", invalid="ignore"):
        wage99_ = wagp / cpi_y * ao.CPI_U[1999] / (wkhp * weeks)
    wage_ok = full & weeks.notna().to_numpy() & (wagp >= ao.INCOME_FLOOR_2010) & (wkhp > 0) & (wage99_.to_numpy() > 0)
    d = pd.DataFrame({"pwgtp": rng.integers(10, 200, m).astype(float), "gender": rng.choice(["M", "F"], m), "schl": schl.astype(float),
                      "age": rng.integers(25, 35, m), "w_mkt": w_mkt, "w_home": w_home, "group": grp, "wage_ok": wage_ok,
                      "occ": np.where(grp == K12, rng.choice(list(K12_CODE_LABEL), m), None),
                      "full": full, "wkw": wkw, "wkhp": wkhp, "wagp": wagp, "semp": np.zeros(m), "data_year": data_year, "inc2010": wagp})
    d["yrs"] = ao.schl_to_highgrade(d.schl).clip(9, 17)
    d["wage99"] = wage99_.to_numpy()                                              # the build_sample formula (checked by alt_wage_sample)
    A = AcsData(d, None)
    st = acs_stats(A, A.pw)
    ba = (d.schl >= 21).astype(float)
    check("P: BA share of everyone", abs(st["ba_share|pooled|all|P"] - (d.pwgtp * ba).sum() / d.pwgtp.sum()) < 1e-12)
    tw = d.pwgtp * d.w_mkt * (d.group == K12)
    check("P: BA share of teachers (weights PWGTP x w_mkt)", abs(st["ba_share|pooled|teacher|P"] - (tw * ba).sum() / tw.sum()) < 1e-12)
    nt_w = d.pwgtp * (1 - d.w_mkt * (d.group == K12))
    check("P: BA share of non-teachers (home production half of part-timers included)", abs(st["ba_share|pooled|nonteacher|P"] - (nt_w * ba).sum() / nt_w.sum()) < 1e-12)
    # T2 teacher premium by hand for men
    ws = d[d.wage_ok & d.group.notna() & (d.gender == "M")].copy()
    ws["lw"] = np.log(ws.wage99)
    cm = pd.Series({k: np.average(t.lw, weights=t.pwgtp) for k, t in ws.groupby("group")})
    cnt = d[(d.gender == "M") & (d.w_mkt > 0)].groupby("group").w_mkt.sum()
    comp = float((cnt[cm.index.drop(K12)] * cm.drop(K12)).sum() / cnt[cm.index.drop(K12)].sum())
    check("T2 teacher premium (men) by hand", abs(st["teacher_premium|men|all_nonteachers|T2"] - (cm[K12] - comp)) < 1e-12)
    # BA+ comparison by hand
    wsb = ws[ws.schl >= 21]
    cmb = pd.Series({k: np.average(t.lw, weights=t.pwgtp) for k, t in wsb.groupby("group")})
    cntb = d[(d.gender == "M") & (d.w_mkt > 0) & (d.schl >= 21)].groupby("group").w_mkt.sum()
    idx = [i for i in cmb.index if i != K12]
    compb = float((cntb[idx] * cmb[idx]).sum() / cntb[idx].sum())
    check("T2 teacher premium against BA+ workers (men) by hand", abs(st["teacher_premium|men|ba_nonteachers|T2"] - (cm[K12] - compb)) < 1e-12)
    check("pooled teacher premium is the K-12 wage-sample-weighted average", abs(
        st["teacher_premium|pooled|all_nonteachers|T2"] - (st["teacher_wage_weight|men"] * st["teacher_premium|men|all_nonteachers|T2"]
                                                           + st["teacher_wage_weight|women"] * st["teacher_premium|women|all_nonteachers|T2"])) < 1e-12)
    check("a person's cell weights sum to one", abs(np.bincount(A.s_person, weights=A.s_u, minlength=A.n) - 1).max() < 1e-12)
    # BA+ teachers and core-code teachers, and BA+ minus non-BA teachers, by hand (men)
    tbm = ws[(ws.group == K12) & (ws.schl >= 21)]
    check("T2 BA+ teachers against BA+ workers (men) by hand",
          abs(st["teacher_premium_ba|men|ba_nonteachers|T2"] - (np.average(tbm.lw, weights=tbm.pwgtp) - compb)) < 1e-12)
    check("T2 BA+ teachers against all workers (men) by hand",
          abs(st["teacher_premium_ba|men|all_nonteachers|T2"] - (np.average(tbm.lw, weights=tbm.pwgtp) - comp)) < 1e-12)
    tcm = ws[(ws.group == K12) & ws.occ.isin(CORE_CODES)]
    check("T2 core-code teachers against BA+ workers (men) by hand",
          abs(st["teacher_premium_core|men|ba_nonteachers|T2"] - (np.average(tcm.lw, weights=tcm.pwgtp) - compb)) < 1e-12)
    tn = ws[(ws.group == K12) & (ws.schl < 21)]
    check("BA+ minus non-BA K-12 teachers, raw (men) by hand",
          abs(st["teacher_ba_gap|men|raw|P"] - (np.average(tbm.lw, weights=tbm.pwgtp) - np.average(tn.lw, weights=tn.pwgtp))) < 1e-12)
    check("P teacher premium against BA+ workers is unchanged by the generalisation (all teachers, men)",
          abs(st["teacher_premium|men|ba_nonteachers|P"] - (np.average(ws[ws.group == K12].lw, weights=ws[ws.group == K12].pwgtp)
                                                            - np.average(wsb[wsb.group != K12].lw, weights=wsb[wsb.group != K12].pwgtp))) < 1e-12)
    # who is outside the wage sample: full-time workers by weeks bracket, by hand (pooled teachers)
    ft = d[d.full & d.group.notna()]
    tt = ft[ft.group == K12]
    den = tt.pwgtp.sum()
    check("exclusion shares of full-time teachers by hand",
          abs(st["excl|pooled|teacher|in_wage_sample"] - tt.pwgtp[tt.wage_ok].sum() / den) < 1e-12
          and abs(st["excl|pooled|teacher|wkw3"] - tt.pwgtp[tt.wkw == 3].sum() / den) < 1e-12
          and abs(sum(st[f"excl|pooled|teacher|wkw{b}"] for b in WKW_LABEL) + st["excl|pooled|teacher|wkw_missing"] - 1) < 1e-12
          and abs(st["excl|pooled|teacher|wk48_outside"] - tt.pwgtp[(tt.wkw <= 2) & ~tt.wage_ok].sum() / den) < 1e-12)
    # the weeks-rule sample: more records, only from the added brackets, the 48+ week records unchanged
    d2 = alt_wage_sample(d, {3: 43.5})
    added = d2.wage_ok & ~d.wage_ok
    check("weeks-rule sample adds exactly the 40-47 week full-time records above the income floor",
          bool((d2.wkw[added] == 3).all()) and bool((d2.wage_ok[d.wage_ok]).all())
          and int(added.sum()) == int((d.full & (d.wkw == 3) & (d.wagp >= 1000) & (d.wkhp > 0)).sum())
          and np.allclose(d2.wage99[added], d.wagp[added] / d.data_year[added].map(ao.CPI_U) * ao.CPI_U[1999] / (d.wkhp[added] * 43.5)))
    A2 = AcsData(d2, None)
    st2 = acs_stats(A2, A2.pw, wage_only=True)
    check("wage-only statistics are a subset of the full ones", set(st2) <= set(st) and "ba_share|pooled|all|P" not in st2)
    # errors-in-variables logit and the latent logit statistics
    nn, Rr, a0, b0 = 150000, 0.9, -1.0, 2.0
    th = rng.standard_normal(nn)
    xx = np.sqrt(Rr) * th + np.sqrt(1 - Rr) * rng.standard_normal(nn)
    yy = (rng.random(nn) < 1 / (1 + np.exp(-(a0 + b0 * th)))).astype(float)
    ae, be = eiv_logit(xx, yy, np.ones(nn), Rr)
    an, bn = wlogit(xx, yy, np.ones(nn))
    check("errors-in-variables logit recovers (a, b); the naive slope over sqrt(R) understates it", abs(ae - a0) < 0.06 and abs(be - b0) < 0.09 and bn / np.sqrt(Rr) < be - 0.02)
    ls = latent_logit_stats(-0.5, 0.0)
    check("latent logit statistics at b = 0: no gap, equal terciles", abs(ls["gap"]) < 1e-9 and max(abs(ls[f"terc{q}"] - ls["share"]) for q in (1, 2, 3)) < 1e-9)
    ls2 = latent_logit_stats(-1.0, 2.0, 0.34)
    check("latent logit statistics: level shift hits the target share, terciles increase, gap positive",
          abs(ls2["share"] - 0.34) < 1e-9 and ls2["terc1"] < ls2["terc2"] < ls2["terc3"] and ls2["gap"] > 0 and abs(sum(ls2[f"terc{q}"] for q in (1, 2, 3)) / 3 - 0.34) < 1e-3)
    # tuition: per-adult and per-earner denominators scale the per-CU ratios
    mu = model_units_by_unit(6700.0)
    check("tuition ratios per adult are 1.9 times those per consumer unit; per earner 1.3",
          all(abs(mu["per_adult"][k] / mu["per_cu"][k] - CE_ADULTS) < 1e-12 and abs(mu["per_earner"][k] / mu["per_cu"][k] - CE_EARNERS) < 1e-12
              for k in ("i_annual_over_ce", "ii_total_over_40ce", "ii_total_over_pv40ce_r", "ii_compounded_total_over_pv40ce_r")))
    # replicate weights that equal PWGTP give a zero SE
    A.rep = np.tile(A.pw, (N_REP, 1))
    se = replicate_se(A, st)
    check("replicate SE is zero when all replicates equal the point weights", max(se.values()) < 1e-12)
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch-repw", action="store_true", help="download the ACS replicate weights (cached in raw/acs_pums/acs5_2013_repw) and stop")
    ap.add_argument("--states", type=int, nargs="*", default=None, help="with --fetch-repw: FIPS codes (default all)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--B", type=int, default=NLSY_B, help="NLSY79 household bootstrap draws")
    ap.add_argument("--seed", type=int, default=NLSY_SEED)
    ap.add_argument("--quick", action="store_true", help="NLSY B = 20 and 8 ACS replicates: development only, the SEs are not valid")
    ap.add_argument("--md-only", action="store_true", help="rebuild the .md and .toml from the saved .json")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--outdir", default=None, help="default: estimates/settings")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    if a.selftest:
        return selftest()
    if a.fetch_repw:
        fetch_repw(a.workers, a.states)
        return 0
    outdir = Path(a.outdir).resolve() if a.outdir else ao.EST / "settings"
    outdir.mkdir(parents=True, exist_ok=True)
    jpath = outdir / "education_occupation.json"
    if a.md_only:
        P = json.loads(jpath.read_text())
    else:
        t0 = time.time()
        acs = run_acs(8 if a.quick else N_REP)
        share = {g: acs["stats"][f"ba_share|{g}|all|P"]["value"] for g in ("pooled", "men", "women")}
        nl = run_nlsy(20 if a.quick else a.B, a.seed, share)
        confirm_teacher_gaps(nl, ao.EST / "settings" / "nlsy_teachers.json")
        P = {"meta": {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"), "script": "data/spatial/education_occupation.py",
                      "acs_persons": acs["n_persons"], "nlsy_B": nl["meta"]["B"], "nlsy_seed": nl["meta"]["seed"], "quick": bool(a.quick),
                      "acs_replicates": acs["n_replicates"], "definitions": DEFINITIONS},
             "acs": acs, "nces_salary": nces_salary(), "nlsy": nl, "tuition": tuition_block(share["pooled"])}
        LOG.info("all estimates done in %.0fs", time.time() - t0)
    P["targets"] = build_targets(P)
    P = clean(P)
    jpath.write_text(json.dumps(P, indent=1, allow_nan=False))
    write_toml(P["targets"], outdir / "education_targets.toml", P["meta"])
    write_md(P, outdir / "education_occupation.md")
    LOG.info("wrote %s, education_targets.toml and education_occupation.md", jpath)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
