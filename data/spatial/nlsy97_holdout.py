#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy97_holdout.py -- SEALED hold-out sample: NLSY97 counterparts of the NLSY79 residential-composition,
teacher-selection and schooling moments, computed once, before the next model specification is chosen.

The spatial teachers model is being specified against NLSY79 moments (`nlsy_composition.py`,
`nlsy_extensions.py`, `nlsy_teachers.py`).  This script computes the same objects in the NLSY97, which
no specification has seen, and writes them to a sealed file that is opened only AFTER the next
specification is chosen:

  1. origin central-city share (MSA central city vs suburb at the 1997 interview, ages 12-16) by tercile
     of the respondent's OWN ability, and top minus bottom;
  2. city->suburb and suburb->city transition rates (origin -> adult, ages 25-34) by own-ability tercile;
  3. move rate by own-ability tercile;
  4. K-12 teachers at ages 25-34: gap in mean LATENT ability against employed non-teachers (A) and
     against non-teacher BA+ workers (B), also BA+ teachers against BA+ non-teachers (C), latent SD units;
  5. BA+ share by ability tercile, and the BA+ vs non-BA gap in latent ability.  Headline: the universe of
     the NLSY79 objects of `education_occupation.py` (16+ years of schooling, highest grade >= 9, terciles cut
     over that sample itself: `edu9_cutm|c16|share`, `edu9|c16|pooled|gap`); BA+ by degree, and the
     whole-sample tercile cuts, are alternatives.

(The origin central-city share by PARENTAL-income tercile is already known and is not recomputed.)
Added: teaching shares at 25-34 (a check on the occupation coding), the adult (25-34) suburb share by own-
ability tercile (`nlsy_composition.py` section 3), and the variants listed in the methods file.

Review round.  An independent review of the first sealed file led to: design-based SEs (`se_psu`), households
from reciprocated links only, the item-5 headline in the NLSY79 universe (`edu9_cutm`, `edu9`), the CBSA-vs-MSA
universe departure documented with an urban-only sensitivity (`adult_urban`), and armed-forces rounds
identified from the roster (`build_pr97`: a military roster employer held at the interview and no civilian one)
and dropped from every sample, the `ESR_COL` = 4 analogue.  Section 10 of the methods file lists the changes
and every place where a count that bears on a sealed statistic was printed or listed.

Ability.  The NLSY79 `afqt_c` is the equal-weight mean of the four age-normed ASVAB IRT z-scores (AR, WK,
PC, MK), each standardized, re-standardized, with 1979 weights.  The NLSY97 counterpart is the same
composite of the four CAT-ASVAB ability estimates, age-normed by `nlsy_ability.build_n97` (weighted rank-
normal within birth-year x quarter cells, 1997 weights); a difference of group means is divided by
sqrt(R_comp_97), the one-factor ULS reliability of the composite (re-estimated in every draw), to give
latent-SD units, exactly as `nlsy_teachers.py` does with R_comp_pop.

Inference and the seal.  Two bootstraps, B = 500 each, on the same statistics.  (i) The household-cluster
bootstrap (`se`, `ci95_*`), seed 20260929 (the seed `nlsy_extensions.t3` uses for the NLSY97 main table):
every draw resamples round-1 households with replacement, the SAME multiplicities for every statistic
(common random numbers as in `nlsy_composition.py`), re-standardizes the ability composite with the draw
weights, re-estimates R_comp_97 and recomputes every tercile cut.  It is the NLSY97 counterpart of the
NLSY79 SEs.  NLSY97 households are the connected components of the RECIPROCATED round-1 roster links
`YOUTH_SIBID01-04.01` (the public file has no household id).  (ii) The design-based bootstrap (`se_psu`,
`ci95_*_psu`), the one that respects the sample design: the NLSY97 is an area sample, the public file carries
the variance stratum `VSTRAT` (117) and the variance PSU `VPSU` (two per stratum) "to correct for clustering in
the sample design" (NLS), and city/suburb residence is strongly clustered by PSU, so household-cluster SEs
understate the sampling error of the residential statistics.  Each draw is a stratified Rao-Wu (1988) rescaling
bootstrap: in every stratum one of the two PSUs is drawn (n_h - 1 = 1 draw) and its weights are multiplied by
n_h/(n_h - 1) = 2, the other PSU's weights by 0 (households nest in PSUs, so households need no separate
resampling); seed 20260929 + 97; the same multipliers for every statistic.  SEs are bootstrap SDs.

SEALING RULES.  No computed moment value or standard error is printed, logged or returned: values are
written ONLY with json.dump to estimates/settings/sealed/nlsy97_holdout.json.  The terminal shows sample
sizes, missing-data counts and PASS/FAIL results of checks that involve no sealed value (the published
NLSY97 parental-income table and the published R_comp_97 reproduced; the generic functions reproduce the
published NLSY79 numbers).  Do not open the sealed file before the specification is chosen.

    cd data/spatial && python3 nlsy97_holdout.py --check     # sample sizes + checks, no sealed values
    python3 nlsy97_holdout.py --seal                         # B = 500 -> sealed json + methods file
    python3 nlsy97_holdout.py --seal --quick                 # B = 50 (the same files; rerun with --force for B = 500)
    python3 nlsy97_holdout.py --selftest                     # synthetic unit tests, no data
    python3 nlsy97_holdout.py --verify-seal                  # structure of the sealed file, no values

`--seal` prints exactly one line, `sealed: <path> (<k> moments)`, and refuses to overwrite an existing sealed
file without `--force`; `nlsy97_holdout_xcheck.py` recomputes a subset of the sealed statistics independently and
prints only PASS/FAIL.  Reads raw/nlsy97_ability.parquet and
metadata/nlsy97_ability_vars.csv (via `nlsy_ability.Cohort`; the `household`, `jobs`, `design` and `military`
blocks of `fetch_data.NLSY_BLOCKS['nlsy97']` are needed), estimates/nlsy_ability.json and
estimates/nlsy_extensions.json (published numbers, for the checks), raw/nlsy79_ability.parquet (equivalence
checks).  Writes estimates/settings/sealed/nlsy97_holdout.json and
estimates/settings/nlsy97_holdout_methods.md (no values).

Written with help from Claude Code.
"""
from __future__ import annotations

import os

# tiny linear algebra: cap BLAS/OpenMP pools before numpy loads (shared nodes can fail at thread start-up)
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import argparse
import contextlib
import datetime as dt
import hashlib
import io
import json
import logging
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402
import nlsy_composition as nc  # noqa: E402
import nlsy_extensions as ne  # noqa: E402
import nlsy_teachers as nt  # noqa: E402

LOG = logging.getLogger("nlsy97_holdout")

SEED = 20260929                 # nlsy_extensions.t3 bootstraps the NLSY97 main table with this seed
PSU_SEED_OFFSET = 97            # the Rao-Wu PSU draws use SEED + 97 (the household draws use SEED itself)
B_DEFAULT = 500
AGE_LO, AGE_HI = ne.ADULT_AGE   # 25-34
BA_DEGREE = 4                   # CV_HIGHEST_DEGREE_EVER_EDT: 4 = bachelor's, 5 master's, 6 PhD, 7 professional
COLLEGE_HGC = nt.COLLEGE_HGC    # 16 years, the NLSY79 'college graduate' of nlsy_teachers.py
MIN_HGC = na.MIN_HGC            # 9
MAX_LOOP = 15                   # roster positions of the job variables (YEMP_*.01 - .15)
OCC_VALID = (10, 9750)          # 2002 Census 4-digit civilian occupations (codebook bands 10-9750)
OCC_MIL = (9800, 9840)          # 'military specific occupations' (codes at a non-military employer; counted, not used:
                                # armed-forces rounds are identified from the roster, see `build_pr97`)
URBAN = 1                       # CV_URBAN-RURAL: 0 rural, 1 urban, 2 unknown
TYPES = nt.TYPES                # prek, elem, sec, sped, nec
K12_CODES = nt.TEACHER_CODES["c02"]    # 2300 pre-K/K, 2310 elementary, 2320 secondary, 2330 special ed, 2340 other
GROUPS = nt.GROUPS              # pooled, men, women
SPEC_VINT = "c02"
COLLEGE_DEFS = ("ba", "c16")    # BA+ by degree (the ACS-style definition); 16+ years of schooling (the NLSY79 definition: headline)
TRANS_VARIANTS = {              # the T3 variants of nlsy_extensions.t3 plus 'latest' (see build_trans97_adult)
    "age12_parent": dict(origin="age12_parent"),
    "age12_youth": dict(origin="age12_youth"),
    "same_region": dict(same_region=True),
    "adult_first25_27": dict(adult="first25_27"),
    "adult_near30": dict(adult="near30"),
    "adult_2004_2011": dict(years=(2004, 2011)),
    "adult_2013_2019": dict(years=(2013, 2019)),
}
ADULT_PREFIX = {"all": "adult", "2x2": "adult_2x2", "urban": "adult_urban"}      # adult-share rows -> key prefix
EXTRA_TRANS_VARIANTS = {        # built here (`build_trans97_adult`), not by `ne.build_trans97`
    "adult_latest25": dict(adult="each person's latest classified round at age >= 25"),
    "adult_urban": dict(adult="rounds at 25-34 with CV_URBAN-RURAL = 1 (urban) only: a partial sensitivity to the "
                              "change of universe from MSA (origin) to CBSA (adult rounds); see the methods file"),
}


# =============================================================================
# Small helpers
# =============================================================================

def occ_valid(x: np.ndarray) -> np.ndarray:
    return np.isfinite(x) & (x >= OCC_VALID[0]) & (x <= OCC_VALID[1])


def occ_military(x: np.ndarray) -> np.ndarray:
    return np.isfinite(x) & (x >= OCC_MIL[0]) & (x <= OCC_MIL[1])


def household_clusters(pid: np.ndarray, sib: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    """Households = connected components of the RECIPROCATED round-1 sibling links.  `pid[i]` is respondent
    i's PUBID and `sib[i, k]` the PUBID of the k-th other NLSY97 respondent in i's household (NaN if none).
    A link i -> j is used only if j -> i is listed as well: 11 of the 5,038 directed links of the public file
    are not reciprocated, and two of those join respondents in different design PSUs (a household lives in
    one PSU), so using them would merge two households that are not one and break the nesting of households
    in PSUs.  Returns (cluster index 0..H-1 per respondent, counts: directed links between two different
    respondents in the sample, links pointing outside the sample, self-links, reciprocated directed links,
    unreciprocated links that were ignored)."""
    n = len(pid)
    row_of = {int(p): i for i, p in enumerate(pid)}
    parent = np.arange(n)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    links: set[tuple[int, int]] = set()
    missing = selfl = 0
    for i in range(n):
        for s in sib[i]:
            if not np.isfinite(s):
                continue
            j = row_of.get(int(s))
            if j is None:
                missing += 1
                continue
            if j == i:
                selfl += 1                      # a respondent listed as his or her own household member
            else:
                links.add((i, j))
    recip = {(i, j) for (i, j) in links if (j, i) in links}
    for i, j in recip:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    roots = np.array([find(i) for i in range(n)])
    _, first, inv = np.unique(roots, return_index=True, return_inverse=True)
    rank = np.empty(len(first), int)
    rank[np.argsort(first)] = np.arange(len(first))        # canonical numbering: households in order of their lowest respondent index
    return rank[inv].astype(int), {"links": len(links), "links_not_in_sample": missing, "self_links_ignored": selfl,
                                   "reciprocated_links": len(recip), "unreciprocated_links_ignored": len(links) - len(recip)}


def design_psus(vstrat: np.ndarray, vpsu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Design PSUs of the NLSY97 public file: the pair (variance stratum `VSTRAT`, variance PSU `VPSU`).
    Returns (PSU index 0..P-1 of every respondent, array (H, 2) with the two PSU indices of each stratum,
    ordered by VPSU = 1, 2).  The public file has exactly two PSUs in each of 117 strata (NLS codebook); a
    different structure raises, because `rao_wu_multipliers` is the n_h = 2 case."""
    if not (np.isfinite(vstrat).all() and np.isfinite(vpsu).all()):
        raise ValueError("VSTRAT / VPSU missing for some respondent")
    s_codes, s_idx = np.unique(vstrat, return_inverse=True)
    H = len(s_codes)
    pairs = np.full((H, 2), -1)
    psu = np.full(len(vstrat), -1)
    nxt = 0
    for h in range(H):
        for k, v in enumerate((1, 2)):
            sel = (s_idx == h) & (vpsu == v)
            if not sel.any():
                raise ValueError(f"stratum {int(s_codes[h])} has no PSU {v}: the n_h = 2 design is assumed")
            psu[sel] = nxt
            pairs[h, k] = nxt
            nxt += 1
        if ((s_idx == h) & ~np.isin(vpsu, (1, 2))).any():
            raise ValueError(f"stratum {int(s_codes[h])}: VPSU outside {{1, 2}}")
    return psu, pairs


def rao_wu_multipliers(pairs: np.ndarray, n_psu: int, rng: np.random.Generator) -> np.ndarray:
    """One draw of the Rao-Wu (1988) rescaling bootstrap for a stratified design with n_h = 2 PSUs per stratum:
    in every stratum h draw m_h = n_h - 1 = 1 PSU with equal probability and multiply the weights of its
    respondents by n_h/(n_h - 1) = 2 (the general multiplier is n_h/(n_h - 1) times the number of times the
    PSU is drawn, because m_h = n_h - 1); the other PSU's respondents get 0.  For any weighted total the
    variance over the draws is sum_h (t_h1 - t_h2)^2, the standard 2-PSU-per-stratum variance estimator.
    Returns the multiplier of every PSU (length `n_psu`); one `rng.integers(0, 2, H)` call per draw."""
    pick = rng.integers(0, 2, pairs.shape[0])
    m = np.zeros(n_psu)
    m[pairs[np.arange(pairs.shape[0]), pick]] = 2.0
    return m


def terciles(af: np.ndarray, w: np.ndarray, origin_age: np.ndarray, ok: np.ndarray
             ) -> tuple[tuple[np.ndarray, tuple[float, float]], tuple[np.ndarray, tuple[float, float]]]:
    """(afqt_origin, afqt_full) of nlsy_composition.groupings: terciles cut over respondents of the origin
    age with a valid score, and over the whole sample with a valid score; weights `w` (draw weights)."""
    t_o = ne.tercile_of(af, w * (origin_age & np.isfinite(af)))
    t_f = ne.tercile_of(af, w * ok)
    return t_o, t_f


def origin_all_stats(orig: np.ndarray, w0: np.ndarray, mult: np.ndarray, terc: np.ndarray) -> dict[str, float]:
    """Origin central-city share by tercile among ALL respondents with a known city/suburb origin (no
    requirement of an adult observation, unlike the 2x2 statistic `origin_share_city`)."""
    wt = w0 * mult
    ok = (orig >= 0) & (orig <= 1) & (terc >= 0) & (wt > 0)
    num = np.bincount(terc[ok], weights=wt[ok] * (orig[ok] == 0), minlength=3)
    den = np.bincount(terc[ok], weights=wt[ok], minlength=3)
    with np.errstate(invalid="ignore", divide="ignore"):
        sh = num / den
    out = {f"share|{k + 1}": float(sh[k]) for k in range(3)}
    out["share|diff"] = float(sh[2] - sh[0])
    return out


# =============================================================================
# Data: NLSY97 person table, person-rounds at 25-34, transitions
# =============================================================================

@dataclass
class PR97:
    """NLSY97 person-rounds at ages 25-34 (interview rounds 2004-2019)."""
    p: np.ndarray        # row in the person table
    Y: np.ndarray        # survey year of the round
    age: np.ndarray
    wr: np.ndarray       # round weight (SAMPLING_WEIGHT_CC / 100)
    occ: np.ndarray      # 2002 Census code of the main job (CV_MAINJOB_FLG roster position), NaN if none
    cur: np.ndarray      # 1 if the main job is held at the interview date (YEMP_CURFLAG), 0 if not, NaN
    deg: np.ndarray      # CV_HIGHEST_DEGREE_EVER_EDT at the round
    hgc: np.ndarray      # CV_HGC_EVER_EDT at the round (95 = ungraded -> NaN)
    doi: np.ndarray      # CV_DOI_EMPLOYED (rounds 2015+), a check on `cur` and on `armed`
    milany: np.ndarray   # bool: some roster employer of the round is a military employer (YEMP_MILFLAG = 1), current or not
    main_mil: np.ndarray  # bool: the main job (CV_MAINJOB_FLG) is a military employer (a check: never true in the public file)
    milc: np.ndarray     # bool: a roster employer with YEMP_MILFLAG = 1 is held at the interview (YEMP_CURFLAG = 1)
    civc: np.ndarray     # bool: a roster employer that is not military is held at the interview
    armed: np.ndarray    # bool: serves in the armed forces at the interview and holds no civilian job (milc & ~civc)


def build_pr97(c: na.Cohort) -> PR97:
    """Person-rounds at 25-34.  `armed` is the analogue of the NLSY79 `ESR_COL` = 4 (armed forces at the interview
    date): the respondent holds a military roster job (`YEMP_MILFLAG.NN` = 1 and `YEMP_CURFLAG.NN` = 1) and no
    civilian job (no roster employer with `YEMP_CURFLAG.NN` = 1 and `YEMP_MILFLAG.NN` != 1).  A reservist or guard
    member who also holds a civilian job is a civilian worker, as in the CPS and in NLS's own `CV_DOI_EMPLOYED`
    (2 = 'military service (but no job)', 2015+ rounds, which this rule reproduces).  `CV_MAINJOB_FLG` never
    points to a military employer, so active-duty rounds with no civilian job have no main job."""
    n = c.df.shape[0]
    cols: dict[str, list[np.ndarray]] = {k: [] for k in ("p", "Y", "age", "wr", "occ", "cur", "deg", "hgc", "doi",
                                                         "milany", "main_mil", "milc", "civc", "armed")}
    for Y in sorted(int(y) for y in c.by_q["CV_AGE_INT_DATE"] if str(y).isdigit()):
        age = c.get("CV_AGE_INT_DATE", Y)
        keep = np.flatnonzero(np.isfinite(age) & (age >= AGE_LO) & (age <= AGE_HI))
        if keep.size == 0:
            continue
        main = ne.get_or_nan(c, "CV_MAINJOB_FLG", Y, n)
        occ, cur = np.full(n, np.nan), np.full(n, np.nan)
        milc, civc = np.zeros(n, bool), np.zeros(n, bool)
        milany, main_mil = np.zeros(n, bool), np.zeros(n, bool)
        for loop in range(1, MAX_LOOP + 1):
            mf = c.try_get(f"YEMP_MILFLAG.{loop:02d}", Y)
            is_mil = (mf == 1) if mf is not None else np.zeros(n, bool)
            milany |= is_mil
            main_mil |= is_mil & (main == loop)
            f = c.try_get(f"YEMP_CURFLAG.{loop:02d}", Y)
            if f is not None:
                milc |= is_mil & (f == 1)
                civc |= ~is_mil & (f == 1)
            m = main == loop
            if not m.any():
                continue
            o = c.try_get(f"YEMP_OCCODE-2002.{loop:02d}", Y)
            if o is not None:
                occ[m] = o[m]
            if f is not None:
                cur[m] = f[m]
        hgc = ne.get_or_nan(c, "CV_HGC_EVER_EDT", Y, n)
        hgc = np.where(hgc >= 95, np.nan, hgc)
        vals = dict(p=keep, Y=np.full(keep.size, Y), age=age[keep], wr=c.get("SAMPLING_WEIGHT_CC", Y)[keep] / 100.0,
                    occ=occ[keep], cur=cur[keep], deg=ne.get_or_nan(c, "CV_HIGHEST_DEGREE_EVER_EDT", Y, n)[keep],
                    hgc=hgc[keep], doi=ne.get_or_nan(c, "CV_DOI_EMPLOYED", Y, n)[keep],
                    milany=milany[keep], main_mil=main_mil[keep], milc=milc[keep], civc=civc[keep],
                    armed=(milc & ~civc)[keep])
        for k, v in vals.items():
            cols[k].append(v)
    return PR97(**{k: np.concatenate(v) for k, v in cols.items()})


def build_trans97_adult(c: na.Cohort, base: ne.Trans, ages: tuple[int, int] = (AGE_LO, AGE_HI),
                        pick: str = "all", require_origin: bool = True, urban_only: bool = False) -> ne.Trans:
    """`ne.build_trans97` (origin = round 1, adult = pooled rounds at 25-34) with the adult rows generalized:
    `ages` is the adult age window and `pick` is 'all' (every classified round, as in T3) or 'latest' (only
    each person's latest classified round).  The filters are those of `ne.build_trans97`; `base` supplies
    the origin arrays; `require_origin=False` keeps persons with no classified origin (the adult-residence
    table of `nlsy_composition.adult_rows(c, None)`); `urban_only=True` keeps rounds with `CV_URBAN-RURAL` =
    1 only (the adult rounds are classified by CBSA, metropolitan plus micropolitan, the origin by the 1990
    MSA standard: a partial sensitivity to that change of universe)."""
    n = c.df.shape[0]
    orig = base.orig
    P, A, W, YY = [], [], [], []
    for Y in sorted(int(y) for y in c.by_q["CV_MSA"] if str(y).isdigit()):
        if not (c.has("CV_AGE_INT_DATE", Y) and c.has("SAMPLING_WEIGHT_CC", Y)):
            continue
        age = c.get("CV_AGE_INT_DATE", Y)
        m = c.get("CV_MSA", Y)
        w = c.get("SAMPLING_WEIGHT_CC", Y) / 100.0
        ok = (age >= ages[0]) & (age <= ages[1]) & np.isfinite(m) & np.isfinite(w) & (w > 0)
        if require_origin:
            ok &= orig >= 0
        if urban_only:
            ok &= ne.get_or_nan(c, "CV_URBAN-RURAL", Y, n) == URBAN
        t = ne.msa_type97(m)
        ok &= t >= 0
        i = np.flatnonzero(ok)
        P.append(i)
        A.append(t[i])
        W.append(w[i])
        YY.append(np.full(i.size, Y))
    p, a, w, y = (np.concatenate(v) for v in (P, A, W, YY))
    if pick == "latest":
        o = np.lexsort((y, p))
        p, a, w, y = p[o], a[o], w[o], y[o]
        last = np.r_[p[1:] != p[:-1], True]
        p, a, w = p[last], a[last], w[last]
    elif pick != "all":
        raise ValueError(pick)
    df = ne.collapse_rows(p, a, w)
    return ne.Trans(f"NLSY97 r1/{pick}/{ages[0]}-{ages[1]}", n, base.cluster, orig, base.origw, base.income,
                    base.inc_pop, df)


@dataclass
class Spec97:
    name: str
    samp: str            # "emp" | "any" | "pop"
    defn: str            # key of nt.DEFS
    rule: str            # "maj" | "half" | "ever"
    hgc9: bool = False
    wmode: str = "w97"   # "w97" | "cs" | "unit" | "round"
    member: np.ndarray | None = None
    teacher: np.ndarray | None = None
    nonteacher: np.ndarray | None = None

    def describe(self) -> dict[str, Any]:
        return {"vintage": SPEC_VINT, "sample": self.samp, "definition": self.defn, "rule": self.rule,
                "hgc_ge_9": self.hgc9, "weights": self.wmode}


def make_specs97() -> list[Spec97]:
    """The NLSY79 grid of `nlsy_teachers.make_specs` (1980 -> 2002 Census codes; no 1970 codes; the
    cross-sectional-sample check `xs` becomes `cs`)."""
    specs: list[Spec97] = []

    def add(samp: str, defn: str, rule: str, hgc9: bool = False, wmode: str = "w97") -> None:
        name = "-".join([SPEC_VINT, samp, defn, rule] + (["hgc9"] if hgc9 else [])
                        + ([] if wmode == "w97" else [wmode]))
        specs.append(Spec97(name, samp, defn, rule, hgc9, wmode))

    for samp in ("emp", "any"):
        for defn in ("all", "core"):
            for rule in ("maj", "half", "ever"):
                add(samp, defn, rule)
        for defn in ("prek", "elem", "sec", "sped", "nec"):
            add(samp, defn, "maj")
        add(samp, "all", "maj", hgc9=True)
        for wmode in ("cs", "unit", "round"):
            add(samp, "all", "maj", wmode=wmode)
    for rule in ("maj", "half", "ever"):
        add("pop", "all", rule)
    return specs


class Frame97:
    """Everything that does not change across bootstrap draws."""

    def __init__(self, c: na.Cohort) -> None:
        self.c = c
        n = self.n = c.df.shape[0]
        R = self.R97 = na.build_n97(c)              # female, w (1997 weight), zn (age-normed z) -- as block V
        self.w0 = R.w
        self.zn = R.zn
        self.zok = np.isfinite(R.zn).all(axis=1) & np.isfinite(R.w)
        self.female = R.female.astype(bool)
        self.cs = c.get("CV_SAMPLE_TYPE", 1997) == 1
        age0 = c.get("CV_AGE_INT_DATE", 1997)
        self.origin_age = ((age0 >= ne.ORIGIN_AGE_97[0]) & (age0 <= ne.ORIGIN_AGE_97[1])
                           & np.isfinite(self.w0) & (self.w0 > 0))
        self.pct = c.get("ASVAB_MATH_VERBAL_SCORE_PCT", "XRND") / 1000.0    # NLS's own AFQT-like percentile
        # households (reciprocated round-1 roster links) and the sample design (variance stratum x PSU)
        pid = c.get("PUBID", 1997)
        sib = np.column_stack([c.get(f"YOUTH_SIBID0{k}.01", 1997) for k in range(1, 5)])
        self.hh_idx, self.hh_info = household_clusters(pid, sib)
        self.ncl = int(self.hh_idx.max()) + 1
        self.vstrat, self.vpsu = c.get("VSTRAT"), c.get("VPSU")
        self.psu_idx, self.psu_pairs = design_psus(self.vstrat, self.vpsu)
        self.n_psu = int(self.psu_idx.max()) + 1
        self.n_strata = int(self.psu_pairs.shape[0])
        span = pd.Series(self.psu_idx).groupby(self.hh_idx).nunique()
        self.hh_info["households_spanning_psu"] = int((span > 1).sum())
        if self.hh_info["households_spanning_psu"]:
            raise RuntimeError("a household spans more than one design PSU: households no longer nest in PSUs")
        # residential transitions (T3 machinery of nlsy_extensions): person clusters in `T0`, used only
        # for the published-table check; the bootstrap here uses household / PSU multiplicities
        self.T0 = ne.build_trans97(c)
        self.T = self.T0
        self.orig = self.T.orig
        self.TV: dict[str, ne.Trans] = {}
        for name, kw in TRANS_VARIANTS.items():
            self.TV[name] = ne.build_trans97(c, **kw)
        self.TV["adult_latest25"] = build_trans97_adult(c, self.T0, (AGE_LO, 200), "latest")
        self.TV["adult_urban"] = build_trans97_adult(c, self.T0, (AGE_LO, AGE_HI), "all", True, True)
        self.rows_match_t3 = self._rows_match(build_trans97_adult(c, self.T0, (AGE_LO, AGE_HI), "all"), self.T0)
        # adult (25-34) residence rows for the adult suburb share by own ability (nlsy_composition section 3):
        # every classified round, no origin requirement ("all"), the 2x2-sample version ("2x2"), and the
        # urban-only sensitivity ("urban", no origin requirement)
        self.adult_rows = {"all": nc.rows_arrays(build_trans97_adult(c, self.T0, (AGE_LO, AGE_HI), "all", False).rows),
                           "2x2": nc.rows_arrays(self.T0.rows),
                           "urban": nc.rows_arrays(build_trans97_adult(c, self.T0, (AGE_LO, AGE_HI), "all", False, True).rows)}
        # person-rounds at 25-34: occupation of the main job, armed-forces status, schooling
        self.pr = pr = build_pr97(c)
        self.mil_occ = occ_military(pr.occ)           # military-specific occupation codes (counted only)
        self.valid = occ_valid(pr.occ)
        self.typ = {t: self.valid & np.isin(pr.occ, K12_CODES[t]) for t in TYPES}
        self.k12 = np.any([self.typ[t] for t in TYPES], axis=0)
        self.core = np.any([self.typ[t] for t in nt.DEFS["core"]], axis=0)
        self.base = ~pr.armed                         # armed-forces rounds are dropped (ESR_COL = 4 analogue)
        self.emp = pr.cur == 1
        self.cnt: dict[str, dict[str, Any]] = {}
        for samp in ("any", "emp"):
            m = self.base & self.valid & (self.emp if samp == "emp" else True)
            nr = np.bincount(pr.p[m], minlength=n)
            ntyp = {t: np.bincount(pr.p[m & self.typ[t]], minlength=n) for t in TYPES}
            okw = m & np.isfinite(pr.wr) & (pr.wr > 0)
            sw = np.bincount(pr.p[okw], weights=pr.wr[okw], minlength=n)
            cw = np.bincount(pr.p[okw], minlength=n)
            self.cnt[samp] = {"mask": m, "nr": nr, "nt": ntyp, "wround": np.where(cw > 0, sw / np.maximum(cw, 1), 0.0)}
        self.observed = np.bincount(pr.p[self.base], minlength=n) > 0        # observed at 25-34, armed-forces rounds dropped
        self.observed_any = np.bincount(pr.p, minlength=n) > 0                # observed at 25-34 in any round
        dm, hm = np.full(n, np.nan), np.full(n, np.nan)
        np.fmax.at(dm, pr.p, pr.deg)                  # highest over all rounds at 25-34 (as `nlsy_teachers.Frame.hg_max`)
        np.fmax.at(hm, pr.p, pr.hgc)
        self.deg_max, self.hgc_max = dm, hm
        self.colleges = {"ba": np.isfinite(dm) & (dm >= BA_DEGREE), "c16": np.isfinite(hm) & (hm >= COLLEGE_HGC)}
        self.edu_known = {"ba": np.isfinite(dm), "c16": np.isfinite(hm)}
        self.hgc9 = np.isfinite(hm) & (hm >= MIN_HGC)
        # schooling samples of item 5, as the NLSY79 `education_occupation.evaluate_ba` sample `all`: persons observed
        # at 25-34 (armed-forces rounds dropped) with the schooling variable known
        self.edu_member = {d: self.observed & self.edu_known[d] for d in COLLEGE_DEFS}
        self.masks = {"pooled": np.ones(n, bool), "men": ~self.female, "women": self.female}
        self.specs = make_specs97()
        for sp in self.specs:
            self._build(sp)

    @staticmethod
    def _rows_match(a: ne.Trans, b: ne.Trans) -> bool:
        x = pd.DataFrame({"p": a.p, "a": a.a, "rw": a.rw}).sort_values(["p", "a", "rw"]).reset_index(drop=True)
        y = pd.DataFrame({"p": b.p, "a": b.a, "rw": b.rw}).sort_values(["p", "a", "rw"]).reset_index(drop=True)
        return bool(len(x) == len(y) and np.array_equal(x.p, y.p) and np.array_equal(x.a, y.a)
                    and np.allclose(x.rw, y.rw))

    def _build(self, sp: Spec97) -> None:
        R = self.cnt["any" if sp.samp == "pop" else sp.samp]
        nr = R["nr"]
        k = sum(R["nt"][t] for t in nt.DEFS[sp.defn])
        k_all = sum(R["nt"][t] for t in nt.DEFS["all"])
        member = self.observed if sp.samp == "pop" else nr > 0
        if sp.hgc9:
            member = member & self.hgc9
        if sp.wmode == "cs":
            member = member & self.cs
        teacher_all = nt.classify(nr, k_all, sp.rule) & member
        sp.member = member
        sp.teacher = nt.classify(nr, k, sp.rule) & member
        sp.nonteacher = member & ~teacher_all          # persons in other K-12 types are in neither group


# =============================================================================
# Statistics (re-run on every bootstrap draw)
# =============================================================================

def gap_stats(prefix: str, a: np.ndarray, sR: float, wv: np.ndarray, member: np.ndarray, teacher: np.ndarray,
              nonteacher: np.ndarray, masks: dict[str, np.ndarray], colleges: dict[str, np.ndarray]
              ) -> dict[str, float]:
    """Teacher gaps in latent AFQT for one specification: the inner loop of `nlsy_teachers.evaluate`, with
    the college-graduate flag generalized to several definitions.  Keys `prefix|group|comparison|stat`:
    comparison nt = A (all non-teachers); col_<d> = B (non-teacher college graduates, all teachers);
    cc_<d> = C (college-graduate teachers vs non-teacher college graduates); stat = mean_t, mean_c (latent
    SD from the sample mean 0), gap, gap_obs (afqt_c units); plus `share`, `col_share_t_<d>`, `col_share_nt_<d>`."""
    ok = member & np.isfinite(a) & (wv > 0)
    t, nonT = teacher & ok, nonteacher & ok
    out: dict[str, float] = {}

    def put(k: str, y1: float, y0: float) -> None:
        out[f"{k}|mean_t"] = y1 / sR
        out[f"{k}|mean_c"] = y0 / sR
        out[f"{k}|gap"] = (y1 - y0) / sR
        out[f"{k}|gap_obs"] = y1 - y0

    for g, gm in masks.items():
        k = f"{prefix}|{g}"
        m1, m0 = t & gm, nonT & gm
        mt, mc = nt.wm(a, wv, m1), nt.wm(a, wv, m0)
        out[f"{k}|share"] = nt.ratio(wv[m1].sum(), wv[m1 | m0].sum())
        put(f"{k}|nt", mt, mc)
        for d, col in colleges.items():
            m1c, m0c = m1 & col, m0 & col
            mtc, mcc = nt.wm(a, wv, m1c), nt.wm(a, wv, m0c)
            out[f"{k}|col_share_t_{d}"] = nt.ratio(wv[m1c].sum(), wv[m1].sum())
            out[f"{k}|col_share_nt_{d}"] = nt.ratio(wv[m0c].sum(), wv[m0].sum())
            put(f"{k}|col_{d}", mt, mcc)
            put(f"{k}|cc_{d}", mtc, mcc)
    return out


def edu_stats(prefix: str, af: np.ndarray, sR: float, w: np.ndarray, terc: np.ndarray, member: np.ndarray,
              flag: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, float]:
    """BA+ (or `flag`) share by ability tercile and the BA+ vs non-BA gap in latent ability, over `member`.
    Keys: `prefix|share|{1,2,3}`, `diff` (3-1), `diff_21`, `diff_32`; `prefix|overall|group`; and
    `prefix|group|{mean_t,mean_c,gap,gap_obs}` (flag vs not, latent SD)."""
    ok = member & np.isfinite(af) & (w > 0)
    out: dict[str, float] = {}
    mt = ok & (terc >= 0)
    num = np.bincount(terc[mt], weights=w[mt] * flag[mt], minlength=3)
    den = np.bincount(terc[mt], weights=w[mt], minlength=3)
    with np.errstate(invalid="ignore", divide="ignore"):
        sh = num / den
    for k in range(3):
        out[f"{prefix}|share|{k + 1}"] = float(sh[k])
    out[f"{prefix}|share|diff"] = float(sh[2] - sh[0])
    out[f"{prefix}|share|diff_21"] = float(sh[1] - sh[0])
    out[f"{prefix}|share|diff_32"] = float(sh[2] - sh[1])
    for g, gm in masks.items():
        m = ok & gm
        out[f"{prefix}|overall|{g}"] = nt.ratio(w[m & flag].sum(), w[m].sum())
        y1, y0 = nt.wm(af, w, m & flag), nt.wm(af, w, m & ~flag)
        out[f"{prefix}|{g}|mean_t"] = y1 / sR
        out[f"{prefix}|{g}|mean_c"] = y0 / sR
        out[f"{prefix}|{g}|gap"] = (y1 - y0) / sR
        out[f"{prefix}|{g}|gap_obs"] = y1 - y0
    return out


def terc_share_stats(prefix: str, af: np.ndarray, w: np.ndarray, terc: np.ndarray, member: np.ndarray,
                     flag: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, float]:
    """P(flag | ability tercile) as `education_occupation.evaluate_ba` computes `terc|<sample>|<group>|ba_share|q`:
    the tercile codes `terc` are cut once over the pooled sample (`member`, finite ability, positive weight) and
    the weighted share of `flag` within tercile q is computed for each group.  Keys `prefix|share|{1,2,3,diff,
    diff_21,diff_32}` for the pooled group and `prefix|<group>|share|...` for men and women."""
    ok = member & np.isfinite(af) & (w > 0)
    out: dict[str, float] = {}
    for g, gm in masks.items():
        mg = ok & gm
        sh = []
        for q in range(3):
            mq = mg & (terc == q)
            sh.append(nt.ratio(w[mq & flag].sum(), w[mq].sum()))
        pre = prefix if g == "pooled" else f"{prefix}|{g}"
        for q in range(3):
            out[f"{pre}|share|{q + 1}"] = sh[q]
        out[f"{pre}|share|diff"] = sh[2] - sh[0]
        out[f"{pre}|share|diff_21"] = sh[1] - sh[0]
        out[f"{pre}|share|diff_32"] = sh[2] - sh[1]
    return out


def share_stats_generic(*, n: int, p: np.ndarray, w_person: np.ndarray, female: np.ndarray, base: np.ndarray,
                        hgc: np.ndarray, emp: np.ndarray, k12: np.ndarray, core: np.ndarray,
                        typs: dict[str, np.ndarray], nr: np.ndarray, nk: np.ndarray, asamp: np.ndarray,
                        masks: dict[str, np.ndarray]) -> dict[str, float]:
    """Teaching shares at 25-34, `nlsy_teachers.share_stats`: person-round shares over the ACS-like
    denominator (highest grade >= 9, not armed forces; `acs`) or without the grade restriction (`all`), weights
    = the person's baseline weight; and person shares (persons classified as K-12 teachers in the employed
    sample over persons observed at 25-34 with highest grade >= 9)."""
    wpy = w_person[p]
    fpy = female[p]
    colr = np.isfinite(hgc) & (hgc >= COLLEGE_HGC)
    gmr = {"pooled": np.ones(len(fpy), bool), "men": ~fpy, "women": fpy}
    out: dict[str, float] = {}
    dens = {"acs": base & np.isfinite(hgc) & (hgc >= MIN_HGC), "all": base}
    for dl, den0 in dens.items():
        for g, gm in gmr.items():
            d = den0 & gm & (wpy > 0)
            tot = wpy[d].sum()
            num = {"emp_k12": emp & k12, "any_k12": k12, "emp_core": emp & core}
            num.update({f"emp_{t}": emp & typs[t] for t in TYPES})
            for nl, nm in num.items():
                out[f"share|{dl}|{nl}|{g}"] = nt.ratio(wpy[d & nm].sum(), tot)
            out[f"share|{dl}|emp_rate|{g}"] = nt.ratio(wpy[d & emp].sum(), tot)
            out[f"share|{dl}|col|{g}"] = nt.ratio(wpy[d & colr].sum(), tot)
            out[f"share|{dl}|emp_k12_given_emp|{g}"] = nt.ratio(wpy[d & emp & k12].sum(), wpy[d & emp].sum())
            out[f"share|{dl}|k12_given_col|{g}"] = nt.ratio(wpy[d & colr & emp & k12].sum(), wpy[d & colr].sum())
            out[f"share|{dl}|k12_given_noncol|{g}"] = nt.ratio(wpy[d & ~colr & emp & k12].sum(), wpy[d & ~colr].sum())
    obs_acs = np.bincount(p[dens["acs"]], minlength=n) > 0
    for rule in ("maj", "half", "ever"):
        tch = nt.classify(nr, nk, rule) & obs_acs & asamp
        for g, gm in masks.items():
            den = obs_acs & asamp & gm
            out[f"pshare|{rule}|{g}"] = nt.ratio(w_person[tch & gm].sum(), w_person[den].sum())
    return out


def evaluate97(fr: Frame97, mult_p: np.ndarray) -> dict[str, float]:
    """All statistics for person multiplicities `mult_p` (ones = the sample itself)."""
    w = fr.w0 * mult_p
    ok = fr.zok & (w > 0)
    af, Rc, *_ = na.pop_afqt(fr.zn, ok, w)
    sR = float(np.sqrt(Rc))
    out: dict[str, float] = {"R_comp_97": float(Rc)}
    # the two standardizations behind the weight sensitivities (cross-sectional sample; unit weights)
    afc, Rcs, *_ = na.pop_afqt(fr.zn, ok & fr.cs, w)
    afu, Ru, *_ = na.pop_afqt(fr.zn, fr.zok, mult_p)
    out["R_comp_97_cs"], out["R_comp_97_unit"] = float(Rcs), float(Ru)
    (terc_o, cuts_o), (terc_f, cuts_f) = terciles(af, w, fr.origin_age, ok)
    out["cut|afqt_origin|1"], out["cut|afqt_origin|2"] = float(cuts_o[0]), float(cuts_o[1])
    out["cut|afqt_full|1"], out["cut|afqt_full|2"] = float(cuts_f[0]), float(cuts_f[1])
    # ---- items 1-3: transitions by own-ability tercile (main sample, T3 NLSY97 definitions) ----------------
    for k, v in nc.trans_stats(fr.T, terc_o, mult_p).items():
        out[f"trans|afqt|{k}"] = v
    for k, v in origin_all_stats(fr.orig, fr.w0, mult_p, terc_o).items():
        out[f"origin_all|afqt|{k}"] = v
    for name, T in fr.TV.items():
        for k, v in nc.trans_stats(T, terc_o, mult_p).items():
            out[f"transv|{name}|afqt|{k}"] = v
    # adult (25-34) suburb share among person-rounds in a CBSA with known status, by tercile (whole-sample cuts)
    for name, rows in fr.adult_rows.items():
        for k, v in nc.adult_stats(rows, terc_f, mult_p).items():
            out[f"{ADULT_PREFIX[name]}|afqt|{k}"] = v
    # tercile assignment by the NLS percentile (`ASVAB_MATH_VERBAL_SCORE_PCT`) instead of the composite
    terc_p, cuts_p = ne.tercile_of(fr.pct, w * (fr.origin_age & np.isfinite(fr.pct)))
    for k, v in nc.trans_stats(fr.T, terc_p, mult_p).items():
        out[f"transv|nls_pct|afqt|{k}"] = v
    # ---- item 4: teachers --------------------------------------------------------------------------------
    A = {"w97": (af, np.sqrt(Rc)), "round": (af, np.sqrt(Rc)), "cs": (afc, np.sqrt(Rcs)), "unit": (afu, np.sqrt(Ru))}
    for sp in fr.specs:
        a, s_ = A[sp.wmode]
        if sp.wmode == "unit":
            wv = mult_p
        elif sp.wmode == "round":
            wv = fr.cnt[sp.samp]["wround"] * mult_p
        else:
            wv = w
        out.update(gap_stats(sp.name, a, float(s_), wv, sp.member, sp.teacher, sp.nonteacher, fr.masks, fr.colleges))
    # ---- item 5: BA+ share by tercile and the BA+ vs non-BA gap --------------------------------------------
    #   edu     whole-sample tercile cuts (`afqt_full`), all persons observed at 25-34 with the schooling variable
    #   edu9    the same, highest grade >= 9 only (gap: the NLSY79 `all_hgc9|<group>|gap` universe, `c16`)
    #   edu_cutm   terciles cut over the analysis sample itself (no grade rule)
    #   edu9_cutm  grade >= 9 and terciles cut over that sample: the universe of the NLSY79 `terc|all_hgc9|...|ba_share`
    for d in COLLEGE_DEFS:
        mem = fr.edu_member[d]                                      # observed at 25-34 (armed-forces rounds dropped), schooling known
        out.update(edu_stats(f"edu|{d}", af, sR, w, terc_f, mem, fr.colleges[d], fr.masks))
        out.update(edu_stats(f"edu9|{d}", af, sR, w, terc_f, mem & fr.hgc9, fr.colleges[d], fr.masks))
        terc_m, _ = ne.tercile_of(af, w * (mem & ok))             # terciles cut over the analysis sample itself
        for k, v in edu_stats(f"edu_cutm|{d}", af, sR, w, terc_m, mem, fr.colleges[d], {"pooled": fr.masks["pooled"]}).items():
            if "|share|" in k:
                out[k] = v
        mem9 = mem & fr.hgc9
        terc_9, cuts_9 = ne.tercile_of(af, w * (mem9 & ok))
        out.update(terc_share_stats(f"edu9_cutm|{d}", af, w, terc_9, mem9, fr.colleges[d], fr.masks))
        out[f"edu9_cutm|{d}|cut|1"], out[f"edu9_cutm|{d}|cut|2"] = float(cuts_9[0]), float(cuts_9[1])
    # ---- teaching shares at 25-34 (a check on the occupation coding) -----------------------------------------
    nk = sum(fr.cnt["emp"]["nt"][t] for t in TYPES)
    out.update(share_stats_generic(n=fr.n, p=fr.pr.p, w_person=w, female=fr.female, base=fr.base & ok[fr.pr.p], hgc=fr.pr.hgc,
                                   emp=fr.emp, k12=fr.k12, core=fr.core, typs=fr.typ, nr=fr.cnt["emp"]["nr"], nk=nk,
                                   asamp=ok, masks=fr.masks))
    return out


def sample_counts(fr: Frame97) -> dict[str, Any]:
    """Unweighted sample sizes (point-estimate classification): these, and only these, may be printed."""
    w = fr.w0
    ok = fr.zok & (w > 0)
    af, Rc, *_ = na.pop_afqt(fr.zn, ok, w)
    (terc_o, _), (terc_f, _) = terciles(af, w, fr.origin_age, ok)
    T = fr.T
    o = T.orig[T.p]
    inn = (o >= 0) & (o <= 1) & (T.a <= 1)
    cnt: dict[str, Any] = {"respondents": int(fr.n), "households": int(fr.ncl), "household_sizes":
                           {int(k): int(v) for k, v in zip(*np.unique(np.bincount(fr.hh_idx), return_counts=True))},
                           "link_checks": fr.hh_info,
                           "valid_weight_1997": int((np.isfinite(fr.w0) & (fr.w0 > 0)).sum()),
                           "four_casasvab_scores_and_weight": int(ok.sum()),
                           "nls_percentile_available": int(np.isfinite(fr.pct).sum()),
                           "origin_age_12_16": int(fr.origin_age.sum()),
                           "origin_age_with_afqt": int((fr.origin_age & ok).sum()),
                           "origin_city_or_suburb": int(((fr.orig == 0) | (fr.orig == 1)).sum()),
                           "origin_classified_in_origin_age": int((fr.orig >= 0).sum()),
                           "origin_with_afqt_city_or_suburb": int((((fr.orig == 0) | (fr.orig == 1)) & ok).sum()),
                           "trans_2x2": {"persons": int(np.unique(T.p[inn]).size), "person_years": int(inn.sum())},
                           "trans_2x2_by_afqt_tercile": {}}
    for k in range(3):
        sel = inn & (terc_o[T.p] == k)
        cnt["trans_2x2_by_afqt_tercile"][str(k + 1)] = {"persons": int(np.unique(T.p[sel]).size),
                                                          "person_years": int(sel.sum())}
    cnt["trans_2x2_with_afqt"] = int(np.unique(T.p[inn & (terc_o[T.p] >= 0)]).size)
    ra = fr.adult_rows["all"]
    cnt["adult_rows_city_or_suburb"] = {"persons": int(np.unique(ra["p"][ra["a"] <= 1]).size), "person_years": int((ra["a"] <= 1).sum())}
    cnt["origin_afqt_tercile_sizes"] = {str(k + 1): int(((terc_o == k) & fr.origin_age).sum()) for k in range(3)}
    # person-rounds at 25-34
    pr = fr.pr
    cnt["person_rounds_25_34"] = int(len(pr.p))
    cnt["persons_observed_25_34"] = int(fr.observed.sum())
    cnt["persons_observed_25_34_any_round"] = int(fr.observed_any.sum())
    cnt["person_rounds_by_year"] = {int(y): int((pr.Y == y).sum()) for y in np.unique(pr.Y)}
    cnt["rounds_with_main_job_occupation"] = int(np.isfinite(pr.occ).sum())
    cnt["rounds_valid_civilian_occupation"] = int((fr.base & fr.valid).sum())
    cnt["rounds_military_occupation_code_not_counted"] = int(fr.mil_occ.sum())
    cnt["rounds_employed_valid_occupation"] = int((fr.base & fr.valid & fr.emp).sum())
    cnt["rounds_missing_current_flag_with_occupation"] = int((np.isfinite(pr.occ) & ~np.isfinite(pr.cur)).sum())
    cnt["rounds_k12_any"] = int((fr.base & fr.k12).sum())
    cnt["rounds_k12_employed"] = int((fr.base & fr.k12 & fr.emp).sum())
    cnt["persons_with_occupation"] = {s: int((fr.cnt[s]["nr"] > 0).sum()) for s in ("any", "emp")}
    cnt["persons_with_schooling"] = {"degree": int(fr.edu_known["ba"].sum()), "hgc": int(fr.edu_known["c16"].sum()),
                                     "degree_and_afqt": int((fr.edu_known["ba"] & ok).sum()),
                                     "hgc_and_afqt": int((fr.edu_known["c16"] & ok).sum()),
                                     "observed_without_degree": int((fr.observed_any & ~fr.edu_known["ba"]).sum()),
                                     "observed_without_hgc": int((fr.observed_any & ~fr.edu_known["c16"]).sum())}
    cnt["schooling_definitions_discordant_persons"] = int((fr.edu_known["ba"] & fr.edu_known["c16"]
                                                           & (fr.colleges["ba"] != fr.colleges["c16"])).sum())
    cnt["gender"] = {"men_observed": int((fr.observed & ~fr.female).sum()), "women_observed": int((fr.observed & fr.female).sum())}
    # armed forces at the interview (the ESR_COL = 4 analogue), design and residence universe: counts only
    doi_ok = np.isfinite(pr.doi) & (pr.doi <= 2)
    doi_mil = doi_ok & (pr.doi == 2)
    cnt["armed_forces"] = {
        "rounds_with_a_military_employer_on_the_roster": int(pr.milany.sum()),
        "rounds_main_job_is_a_military_employer": int(pr.main_mil.sum()),
        "rounds_with_a_current_military_job": int(pr.milc.sum()),
        "rounds_current_military_and_current_civilian_job": int((pr.milc & pr.civc).sum()),
        "rounds_armed_dropped": int(pr.armed.sum()),
        "persons_with_an_armed_round": int(np.unique(pr.p[pr.armed]).size),
        "armed_rounds_with_a_main_job_on_record": int((pr.armed & np.isfinite(pr.occ)).sum()),
        "persons_observed_only_in_armed_rounds": int((fr.observed_any & ~fr.observed).sum()),
        "cv_doi_employed_check": {"rounds_compared": int(doi_ok.sum()),
                                  "armed_and_doi_military": int((pr.armed & doi_mil).sum()),
                                  "armed_not_doi_military": int((pr.armed & doi_ok & ~doi_mil).sum()),
                                  "doi_military_not_armed": int((~pr.armed & doi_mil).sum())}}
    psu_sizes = np.bincount(fr.psu_idx, minlength=fr.n_psu)
    cnt["design"] = {"strata": fr.n_strata, "psus": fr.n_psu,
                     "respondents_per_psu": {"min": int(psu_sizes.min()), "median": float(np.median(psu_sizes)),
                                             "max": int(psu_sizes.max())},
                     "households_spanning_psu": fr.hh_info["households_spanning_psu"]}
    cnt["residence_universe"] = residence_counts(fr.c)
    # teacher/comparison group sizes (persons, unweighted, with a valid score)
    sc: dict[str, Any] = {}
    for sp in fr.specs:
        okp = sp.member & ok
        t, c0 = sp.teacher & okp, sp.nonteacher & okp
        sc[sp.name] = {g: {"teachers": int((t & gm).sum()), "non_teachers": int((c0 & gm).sum()),
                           **{f"non_teacher_{d}": int((c0 & gm & fr.colleges[d]).sum()) for d in COLLEGE_DEFS},
                           **{f"teacher_{d}": int((t & gm & fr.colleges[d]).sum()) for d in COLLEGE_DEFS}}
                       for g, gm in fr.masks.items()}
    cnt["teacher_groups"] = sc
    ed = {}
    for d in COLLEGE_DEFS:
        mem = fr.edu_member[d] & ok
        mem9 = mem & fr.hgc9
        t9, _ = ne.tercile_of(af, w * mem9)
        ed[d] = {"members": int(mem.sum()), "by_afqt_tercile_full": {str(k + 1): int((mem & (terc_f == k)).sum()) for k in range(3)},
                 "members_grade_ge_9": int(mem9.sum()),
                 "grade_ge_9_by_afqt_tercile_cut_over_that_sample": {str(k + 1): int((mem9 & (t9 == k)).sum()) for k in range(3)},
                 "members_grade_below_9_or_unknown": int((mem & ~fr.hgc9).sum())}
    cnt["schooling_samples"] = ed
    return cnt


def residence_counts(c: na.Cohort) -> dict[str, Any]:
    """Unweighted counts of the residence codes behind the city/suburb classification: the origin (1997 `CV_MSA`,
    ages 12-16, positive weight; 1990-standard MSAs) and the adult person-rounds at 25-34 (`CV_MSA` = CBSA,
    metropolitan plus micropolitan, 2000 standards 2004-2011 and 2010 standards from 2013), and `CV_URBAN-RURAL`
    among the adult city/suburb rows."""
    n = c.df.shape[0]
    lab = {1: "not_in_msa_or_cbsa", 2: "suburb", 3: "central_city", 4: "in_msa_or_cbsa_status_unknown", 5: "not_in_country"}
    org = {v: 0 for v in lab.values()}
    adu = {v: 0 for v in lab.values()}
    urb = {"rural": 0, "urban": 0, "unknown": 0, "missing": 0}
    for Y in sorted(int(y) for y in c.by_q["CV_MSA"] if str(y).isdigit()):
        if not (c.has("CV_AGE_INT_DATE", Y) and c.has("SAMPLING_WEIGHT_CC", Y)):
            continue
        age, m, w = c.get("CV_AGE_INT_DATE", Y), c.get("CV_MSA", Y), c.get("SAMPLING_WEIGHT_CC", Y)
        u = ne.get_or_nan(c, "CV_URBAN-RURAL", Y, n)
        base = np.isfinite(m) & np.isfinite(w) & (w > 0) & np.isfinite(age)
        if Y == 1997:
            sel = base & (age >= ne.ORIGIN_AGE_97[0]) & (age <= ne.ORIGIN_AGE_97[1])
            for code, name in lab.items():
                org[name] += int((sel & (m == code)).sum())
        sel = base & (age >= AGE_LO) & (age <= AGE_HI)
        for code, name in lab.items():
            adu[name] += int((sel & (m == code)).sum())
        s2 = sel & ((m == 2) | (m == 3))
        urb["rural"] += int((s2 & (u == 0)).sum())
        urb["urban"] += int((s2 & (u == 1)).sum())
        urb["unknown"] += int((s2 & (u == 2)).sum())
        urb["missing"] += int((s2 & ~np.isfinite(u)).sum())
    return {"origin_1997_ages_12_16": org, "adult_person_rounds_25_34": adu, "adult_city_suburb_rows_by_urban_rural": urb}


# =============================================================================
# Bootstrap
# =============================================================================

def summarize_draws(draws: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Bootstrap SD (ddof = 1), 2.5 and 97.5 percentiles and the number of finite draws, column by column."""
    B = draws.shape[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        se = np.nanstd(draws, axis=0, ddof=1) if B > 1 else np.full(draws.shape[1], np.nan)
        lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    return se, lo, hi, np.sum(np.isfinite(draws), axis=0)


def bootstrap97(fr: Frame97, B: int, seed: int) -> dict[str, Any]:
    """Point estimates and two bootstraps of every statistic: household clusters (`se`, `ci95_*`; the NLSY79
    analogue; draws `rng(seed)`: `rng.integers(0, H, H)` per draw) and design-based Rao-Wu PSU draws (`se_psu`,
    `ci95_*_psu`; draws `rng(seed + PSU_SEED_OFFSET)`: `rng.integers(0, 2, strata)` per draw).  Within a method
    every statistic uses the same draws."""
    point = evaluate97(fr, np.ones(fr.n))
    keys = list(point)
    rng_h = np.random.default_rng(seed)
    rng_p = np.random.default_rng(seed + PSU_SEED_OFFSET)
    dr_h = np.full((B, len(keys)), np.nan)
    dr_p = np.full((B, len(keys)), np.nan)
    for b in range(B):
        cm = np.bincount(rng_h.integers(0, fr.ncl, fr.ncl), minlength=fr.ncl).astype(float)
        mp = rao_wu_multipliers(fr.psu_pairs, fr.n_psu, rng_p)
        for dr, mult in ((dr_h, cm[fr.hh_idx]), (dr_p, mp[fr.psu_idx])):
            try:
                d = evaluate97(fr, mult)
                dr[b] = [d.get(k, np.nan) for k in keys]
            except Exception as exc:                              # a failed draw is a NaN row; never echo values
                LOG.warning("draw %d failed (%s)", b, type(exc).__name__)
    se, lo, hi, n_ok = summarize_draws(dr_h)
    se_p, lo_p, hi_p, n_ok_p = summarize_draws(dr_p)
    stats = {k: {"value": float(point[k]), "se": float(se[i]), "ci95_lo": float(lo[i]), "ci95_hi": float(hi[i]),
                 "n_draws": int(n_ok[i]), "se_psu": float(se_p[i]), "ci95_lo_psu": float(lo_p[i]),
                 "ci95_hi_psu": float(hi_p[i]), "n_draws_psu": int(n_ok_p[i])} for i, k in enumerate(keys)}
    return {"stats": stats, "B": B, "seed": seed, "failed_draws": int(np.sum(~np.isfinite(dr_h).any(axis=1))),
            "failed_draws_psu": int(np.sum(~np.isfinite(dr_p).any(axis=1)))}


# =============================================================================
# Public checks (no sealed value is printed)
# =============================================================================

def ok_str(flag: bool) -> str:
    return "PASS" if flag else "FAIL"


def check_published_97(c: na.Cohort, fr: Frame97, B: int) -> bool:
    """Published NLSY97 numbers rebuilt with this script's data: R_comp_97 (nlsy_ability.json) and the
    parental-income tercile table of T3 (nlsy_extensions.json, person clusters, seed 20260929).  Both are
    unsealed published results."""
    good = True
    pub = json.loads((HERE / "estimates" / "nlsy_ability.json").read_text())
    r_pub = next(e["value"] for e in pub["estimates"] if e["key"] == "R_comp_97")
    w = fr.w0
    ok = fr.zok & (w > 0)
    _, Rc, *_ = na.pop_afqt(fr.zn, ok, w)
    d = abs(Rc - r_pub)
    print(f"  R_comp_97 rebuilt vs published nlsy_ability.json: abs diff {d:.1e}  [{ok_str(d < 1e-10)}]; N = {int(ok.sum())}")
    good &= d < 1e-10
    ext = json.loads((HERE / "estimates" / "nlsy_extensions.json").read_text())["T3"]["nlsy97_main"]
    res = ne.t3_estimate(fr.T0, B, SEED)
    diffs = []
    for t in ("1", "2", "3"):
        for key in nc.TRANS_KEYS:
            diffs.append(abs(res["tercile"][t][key] - ext["tercile"][t][key]))
            if B == 500:
                diffs.append(abs(res["tercile_se"][t][key] - ext["tercile_se"][t][key]))
    for key in ("move_rate", "city_to_suburb", "suburb_to_city", "origin_share_city"):
        diffs.append(abs(res["stats"][key] - ext["stats"][key]))
    d = max(diffs)
    print(f"  published NLSY97 parental-income table (T3, persons as clusters, B = {B}) rebuilt: max abs diff {d:.1e}  "
          f"[{ok_str(d < 1e-9)}]; 2x2 persons {res['n_persons_2x2']} (published {ext['n_persons_2x2']}), "
          f"person-years {res['n_person_years_2x2']} (published {ext['n_person_years_2x2']})")
    good &= d < 1e-9 and res["n_persons_2x2"] == ext["n_persons_2x2"] and res["n_person_years_2x2"] == ext["n_person_years_2x2"]
    print(f"  adult rows of build_trans97_adult(25-34, all) equal those of ne.build_trans97: [{ok_str(fr.rows_match_t3)}]")
    good &= fr.rows_match_t3
    return good


def check_nlsy79_equivalence(B: int = 0) -> bool:
    """The generic functions applied to NLSY79 inputs reproduce `nlsy_teachers.evaluate` / `share_stats` (all
    teacher gaps, shares and teaching shares, at the sample and at one bootstrap draw's multiplicities), the
    published `nlsy_teachers.json`, and `nlsy_composition.groupings` (tercile codes and cuts).  NLSY79
    results are unsealed."""
    good = True
    c79 = na.Cohort("nlsy79", HERE)
    fr = nt.Frame(c79)
    D = fr.D
    rng = np.random.default_rng(SEED)
    mh = np.bincount(rng.integers(0, D.nhh, D.nhh), minlength=D.nhh).astype(float)
    for lab, mult_hh in (("sample", None), ("one bootstrap draw", mh)):
        mult_p = np.ones(D.n) if mult_hh is None else mult_hh[D.hh_idx].astype(float)
        ref = nt.evaluate(fr, mult_hh)
        w79 = D.w * mult_p
        afqt, Rc, *_ = na.pop_afqt(D.z, D.asamp, w79)
        afx, Rx, *_ = na.pop_afqt(D.z, D.asamp & fr.xs, w79)
        afu, Ru, *_ = na.pop_afqt(D.z, D.asamp, mult_p)
        A = {"w79": (afqt, np.sqrt(Rc)), "round": (afqt, np.sqrt(Rc)), "xs": (afx, np.sqrt(Rx)), "unit": (afu, np.sqrt(Ru))}
        mine: dict[str, float] = {}
        for sp in fr.specs:
            a, s_ = A[sp.wmode]
            wv = mult_p if sp.wmode == "unit" else (fr.R[(sp.vint, sp.samp)]["wround"] * mult_p if sp.wmode == "round" else w79)
            r = gap_stats(sp.name, a, float(s_), wv, sp.member, sp.teacher, sp.nonteacher, fr.masks, {"c16": fr.college})
            for k, v in r.items():                                   # map the new key names to the nlsy_teachers ones
                k2 = k.replace("|col_c16|", "|col|").replace("|cc_c16|", "|cc|").replace("|col_share_t_c16", "|col_share_t") \
                      .replace("|col_share_nt_c16", "|col_share_nt")
                mine[k2] = v
        diffs = [abs(mine[k] - ref[k]) for k in mine if np.isfinite(ref.get(k, np.nan)) or np.isfinite(mine[k])]
        miss = [k for k in mine if k not in ref]
        bad = [k for k in mine if k in ref and (np.isfinite(mine[k]) != np.isfinite(ref[k]))]
        d = max(diffs) if diffs else np.nan
        print(f"  NLSY79 {lab}: gap_stats vs nlsy_teachers.evaluate over {len(mine)} statistics: max abs diff {d:.1e}, "
              f"keys missing in reference {len(miss)}, finiteness mismatches {len(bad)}  [{ok_str(d < 1e-12 and not miss and not bad)}]")
        good &= bool(d < 1e-12 and not miss and not bad)
        # teaching shares
        Rm = fr.R[("c80", "emp")]
        nk = sum(Rm["nt"][t] for t in nt.TYPES)
        core = np.any([fr.typ["c80"][t] for t in nt.DEFS["core"]], axis=0)
        sh = share_stats_generic(n=D.n, p=fr.py.p, w_person=w79, female=fr.female, base=fr.base, hgc=fr.py.hgc,
                                 emp=fr.py.esr == 1, k12=fr.k12["c80"], core=core, typs=fr.typ["c80"], nr=Rm["nr"], nk=nk,
                                 asamp=D.asamp, masks=fr.masks)
        d2 = [abs(v - ref[k]) for k, v in sh.items() if k in ref and (np.isfinite(v) or np.isfinite(ref[k]))]
        n_cmp = len([k for k in sh if k in ref])
        print(f"  NLSY79 {lab}: share_stats_generic vs nlsy_teachers.share_stats over {n_cmp} statistics: max abs diff "
              f"{max(d2):.1e}  [{ok_str(max(d2) < 1e-12 and n_cmp == len(sh))}]")
        good &= bool(max(d2) < 1e-12 and n_cmp == len(sh))
    # published nlsy_teachers.json (point estimates at the sample)
    pubj = json.loads((HERE / "estimates" / "settings" / "nlsy_teachers.json").read_text())["statistics"]
    ref = nt.evaluate(fr, None)
    d3 = max(abs(ref[k] - pubj[k]["value"]) for k in pubj if pubj[k]["value"] is not None and k in ref and np.isfinite(ref[k]))
    print(f"  NLSY79 nlsy_teachers.evaluate vs published nlsy_teachers.json (point estimates): max abs diff {d3:.1e}  [{ok_str(d3 < 1e-9)}]")
    good &= bool(d3 < 1e-9)
    # terciles: generic helper vs nlsy_composition.groupings
    cf = nc.Frame(c79)
    cm = np.bincount(np.random.default_rng(SEED + 2).integers(0, cf.ncl, cf.ncl), minlength=cf.ncl).astype(float)
    for lab, m_p in (("sample", np.ones(cf.n)), ("one draw", cm[cf.cl_idx])):
        G = nc.groupings(cf, m_p)
        w = cf.D.w * m_p
        afq = nc.afqt_composite(cf.D, w)
        (to, co), (tf, cf_) = terciles(afq, cf.w0 * m_p, cf.origin_age, cf.D.asamp & (w > 0))
        eq = (np.array_equal(to, G["afqt_origin"][0]) and np.array_equal(tf, G["afqt_full"][0])
              and np.allclose(co, G["afqt_origin"][1]) and np.allclose(cf_, G["afqt_full"][1]))
        print(f"  NLSY79 {lab}: terciles() reproduces nlsy_composition.groupings (codes and cuts): [{ok_str(eq)}]")
        good &= bool(eq)
    # the transition table by own-AFQT tercile against the published nlsy_composition.json
    pubc = json.loads((HERE / "estimates" / "settings" / "nlsy_composition.json").read_text())["statistics"]
    G = nc.groupings(cf, np.ones(cf.n))
    st = nc.trans_stats(cf.T, G["afqt_origin"][0], np.ones(cf.n))
    d4 = max(abs(v - pubc[f"trans|afqt|{k}"]["value"]) for k, v in st.items())
    print(f"  NLSY79 trans_stats by own-AFQT tercile vs published nlsy_composition.json: max abs diff {d4:.1e}  [{ok_str(d4 < 1e-9)}]")
    good &= bool(d4 < 1e-9)
    # item 5: the generic share-by-tercile and gap functions against the published education_occupation.json (sample all_hgc9)
    ed = json.loads((HERE / "estimates" / "settings" / "education_occupation.json").read_text())["nlsy"]["stats"]
    w79 = fr.D.w
    afqt, Rc, *_ = na.pop_afqt(fr.D.z, fr.D.asamp, w79)
    sR = float(np.sqrt(Rc))
    hg = fr.hg_max
    okg = np.isfinite(hg)
    ba = okg & (hg >= nt.COLLEGE_HGC)
    samp = fr.observed & fr.D.asamp & okg & (hg >= na.MIN_HGC)
    m = samp & np.isfinite(afqt) & (w79 > 0)
    t, cuts = ne.tercile_of(np.where(m, afqt, np.nan), w79 * m)
    mine5: dict[str, float] = {}
    for k, v in terc_share_stats("x", afqt, w79, t, samp, ba, fr.masks).items():
        parts = k.split("|")
        g, rest = ("pooled", parts[2:]) if parts[1] == "share" else (parts[1], parts[3:])
        mine5[f"terc|all_hgc9|{g}|ba_share|{rest[0]}"] = v
    for k, v in edu_stats("x", afqt, sR, w79, t, samp, ba, fr.masks).items():
        parts = k.split("|")
        if parts[1] == "share" or parts[1] == "overall":
            continue                                    # the pooled-tercile shares are compared above; `overall` below
        g, stt = parts[1], {"mean_t": "mean_ba", "mean_c": "mean_nba"}.get(parts[2], parts[2])
        mine5[f"all_hgc9|{g}|{stt}"] = v
    for g, gm in fr.masks.items():
        mg = m & gm
        mine5[f"all_hgc9|{g}|ba_share"] = nt.ratio(w79[mg & ba].sum(), w79[mg].sum())
    mine5["terc|all_hgc9|cut1"], mine5["terc|all_hgc9|cut2"] = float(cuts[0]), float(cuts[1])
    miss5 = [k for k in mine5 if k not in ed or ed[k]["value"] is None]
    d5 = max(abs(v - ed[k]["value"]) for k, v in mine5.items() if k in ed and ed[k]["value"] is not None)
    print(f"  NLSY79 sample all_hgc9: terc_share_stats / edu_stats vs published education_occupation.json over {len(mine5) - len(miss5)} statistics: "
          f"max abs diff {d5:.1e}, keys missing in reference {len(miss5)}  [{ok_str(d5 < 1e-9 and not miss5)}]")
    good &= bool(d5 < 1e-9 and not miss5)
    return good


def is_bounded_key(k: str) -> bool:
    """Statistics that are shares or rates, hence in [0, 1] (everything except differences, gaps, means, cuts)."""
    last = k.split("|")[-1]
    if k.startswith(("share|", "pshare|")) or "|overall|" in k or "col_share" in k or k.endswith("|share"):
        return True
    if k.startswith(("adult|", "adult_2x2|", "adult_urban|")) and (k.endswith("|overall") or "|share|" in k and last in ("1", "2", "3")):
        return True
    return (last in ("1", "2", "3") and "|cut|" not in k
            and k.startswith(("trans|", "transv|", "origin_all|", "edu|", "edu9|", "edu_cutm|", "edu9_cutm|")))


def check_structure(fr: Frame97) -> bool:
    """Structural properties of the sealed statistics (finiteness, bounds, ordering of cuts); prints PASS/FAIL
    only.  Substantive properties (signs, magnitudes) are deliberately not checked."""
    good = True
    pt = evaluate97(fr, np.ones(fr.n))
    rng = np.random.default_rng(SEED)
    cm = np.bincount(rng.integers(0, fr.ncl, fr.ncl), minlength=fr.ncl).astype(float)
    dr = evaluate97(fr, cm[fr.hh_idx])
    for lab, d in (("sample", pt), ("one draw", dr)):
        keys = list(d)
        bounded = [k for k in keys if is_bounded_key(k)]
        out_of_range = [k for k in bounded if np.isfinite(d[k]) and not (-1e-12 <= d[k] <= 1 + 1e-12)]
        nan_keys = [k for k in keys if not np.isfinite(d[k])]
        n_nan = len(nan_keys)
        if lab == "sample" and nan_keys:
            byblk: dict[str, int] = {}
            for k in nan_keys:
                byblk[k.split("|")[0] + "|" + k.split("|")[1]] = byblk.get(k.split("|")[0] + "|" + k.split("|")[1], 0) + 1
            print("    non-finite statistics at the sample, by spec|group (cells without observations): " + str(byblk))
        cuts_ok = (d["cut|afqt_origin|1"] < d["cut|afqt_origin|2"]) and (d["cut|afqt_full|1"] < d["cut|afqt_full|2"])
        r_ok = all(0 < d[k] < 1 for k in ("R_comp_97", "R_comp_97_cs", "R_comp_97_unit"))
        ok = not out_of_range and cuts_ok and r_ok
        print(f"  {lab}: {len(keys)} statistics, {n_nan} not finite, bounded shares in [0, 1]: {not out_of_range}, "
              f"ordered tercile cuts: {cuts_ok}, reliabilities in (0, 1): {r_ok}  [{ok_str(ok)}]")
        good &= ok
    # tercile masses (weights) are about equal
    w = fr.w0
    ok_ = fr.zok & (w > 0)
    af, *_ = na.pop_afqt(fr.zn, ok_, w)
    (to, _), (tf, _) = terciles(af, w, fr.origin_age, ok_)
    for lab, t, pop in (("origin-age population", to, fr.origin_age & ok_), ("whole sample", tf, ok_)):
        mass = np.array([w[(t == k) & pop].sum() for k in range(3)]) / w[pop].sum()
        flag = bool(np.all(np.abs(mass - 1 / 3) < 0.02))
        print(f"  terciles over the {lab} carry weight shares within 0.02 of one third: [{ok_str(flag)}]")
        good &= flag
    # CURFLAG against CV_DOI_EMPLOYED in the rounds where the latter exists (counts only)
    pr = fr.pr
    both = np.isfinite(pr.doi) & np.isfinite(pr.cur) & (pr.doi <= 1)
    disc = int((both & (pr.doi != pr.cur)).sum())
    print(f"  employment at the interview: main-job YEMP_CURFLAG vs CV_DOI_EMPLOYED in rounds 2015+: {int(both.sum())} "
          f"person-rounds compared, {disc} discordant")
    return good


def check_design(fr: Frame97, B: int = 500) -> bool:
    """The design variables and the PSU bootstrap, on published/unsealed quantities only: 117 strata x 2 PSUs (the NLS
    codebook), households nested in PSUs, the Rao-Wu multipliers have mean one, and (information, no ability statistic is
    involved) the ratio of the design-based SE to the household-cluster SE for the published parental-income T3 statistics."""
    good = fr.n_strata == 117 and fr.n_psu == 234 and fr.hh_info["households_spanning_psu"] == 0
    rng = np.random.default_rng(SEED + PSU_SEED_OFFSET)
    mm = np.mean([rao_wu_multipliers(fr.psu_pairs, fr.n_psu, rng) for _ in range(2000)], axis=0)
    mean_ok = bool(np.abs(mm - 1.0).max() < 0.1)
    good &= mean_ok
    print(f"  design: {fr.n_strata} variance strata x 2 PSUs = {fr.n_psu} PSUs (NLS: 117 strata, two PSUs each); households spanning a PSU: "
          f"{fr.hh_info['households_spanning_psu']}; Rao-Wu multipliers have mean 1 +- 0.1 for every PSU over 2000 draws: {mean_ok}  [{ok_str(good)}]")
    T = fr.T0
    keys = ("city_to_suburb", "suburb_to_city", "move_rate", "origin_share_city")

    def vec(mult: np.ndarray) -> list[float]:
        pt = T.matrices(mult)
        sts = [ne.stats_from(pt["terc"][k]) for k in range(3)]
        out: list[float] = []
        for key in keys:
            out += [sts[k][key] for k in range(3)] + [sts[2][key] - sts[0][key]]
        out += [ne.stats_from(pt["all"])[k] for k in keys]
        return out

    rng_h = np.random.default_rng(SEED)
    rng_p = np.random.default_rng(SEED + PSU_SEED_OFFSET)
    dh, dp = [], []
    for _ in range(B):
        cm = np.bincount(rng_h.integers(0, fr.ncl, fr.ncl), minlength=fr.ncl).astype(float)
        dh.append(vec(cm[fr.hh_idx]))
        dp.append(vec(rao_wu_multipliers(fr.psu_pairs, fr.n_psu, rng_p)[fr.psu_idx]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        sd_p = np.nanstd(dp, axis=0, ddof=1)
        ratio = sd_p / np.nanstd(dh, axis=0, ddof=1)
    # the Rao-Wu SE of a ratio against the textbook linearization SE for a two-PSU-per-stratum design (the whole-sample origin
    # city share of the parental-income T3 sample: a ratio of person totals); Var = sum_h (Z_h1 - Z_h2)^2 with
    # Z_psu = sum of z_i = (y_i - r x_i) / sum x over the PSU's respondents
    o = T.orig[T.p]
    inn = (o >= 0) & (o <= 1) & (T.a <= 1)
    x_i = np.bincount(T.p[inn], weights=T.rw[inn], minlength=fr.n)
    y_i = np.bincount(T.p[inn], weights=T.rw[inn] * (o[inn] == 0), minlength=fr.n)
    z_i = (y_i - (y_i.sum() / x_i.sum()) * x_i) / x_i.sum()
    Z = np.bincount(fr.psu_idx, weights=z_i, minlength=fr.n_psu)
    se_lin = float(np.sqrt(((Z[fr.psu_pairs[:, 0]] - Z[fr.psu_pairs[:, 1]]) ** 2).sum()))
    lin_ratio = float(sd_p[-1] / se_lin)
    lin_ok = abs(lin_ratio - 1.0) < 0.10
    good &= bool(lin_ok)
    parts = []
    for i, key in enumerate(keys):
        r = ratio[4 * i:4 * i + 4]
        parts.append(f"{key} T1/T2/T3/T3-T1 {r[0]:.2f}/{r[1]:.2f}/{r[2]:.2f}/{r[3]:.2f}")
    allr = ratio[16:]
    print(f"  design effect on the published parental-income T3 table (B = {B}): SE(PSU) / SE(household), " + "; ".join(parts)
          + "; whole sample " + "/".join(f"{x:.2f}" for x in allr) + " (city_to_suburb, suburb_to_city, move_rate, origin_share_city)")
    print(f"  Rao-Wu bootstrap SE of the whole-sample origin city share (parental-income T3 sample) over the textbook linearization SE "
          f"for a ratio in a two-PSU-per-stratum design: {lin_ratio:.2f} (within 10%; the Monte Carlo error of a bootstrap SD with B = {B} is about 3%)  [{ok_str(lin_ok)}]")
    return bool(good)


def run_checks(c: na.Cohort, fr: Frame97, with79: bool = True) -> tuple[bool, list[str]]:
    """Every check that involves no sealed value, with its printed lines captured (so `--seal` can keep the
    terminal to one line and record the results in the methods file)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        good = check_published_97(c, fr, 500)
        good &= check_design(fr)
        good &= check_structure(fr)
        if with79:
            good &= check_nlsy79_equivalence()
    return bool(good), buf.getvalue().splitlines()


def public_counts(cnt: dict[str, Any]) -> dict[str, Any]:
    """The counts that may leave the sealed file: sample sizes and missing-data counts.  Counts of BA+ (or
    16+ years) persons inside the teacher and comparison groups are dropped: their ratios to the group
    sizes are unweighted versions of sealed statistics (`col_share_t_*`, `col_share_nt_*`, the BA+ share),
    so they stay in the sealed file's `counts` only."""
    out = dict(cnt)
    out["teacher_groups"] = {sp: {g: {"teachers": v["teachers"], "non_teachers": v["non_teachers"]} for g, v in gv.items()}
                             for sp, gv in cnt["teacher_groups"].items()}
    return out


def print_counts(cnt: dict[str, Any]) -> None:
    print("Sample sizes and missing-data counts (unweighted persons / person-rounds):")
    for k in ("respondents", "valid_weight_1997", "four_casasvab_scores_and_weight", "nls_percentile_available",
              "origin_age_12_16", "origin_age_with_afqt", "origin_classified_in_origin_age", "origin_city_or_suburb",
              "origin_with_afqt_city_or_suburb"):
        print(f"  {k}: {cnt[k]}")
    print(f"  households: {cnt['households']} (sizes {cnt['household_sizes']}); links {cnt['link_checks']}")
    d = cnt["design"]
    print(f"  design: {d['strata']} variance strata, {d['psus']} variance PSUs; respondents per PSU {d['respondents_per_psu']}; "
          f"households spanning more than one PSU: {d['households_spanning_psu']}")
    ru = cnt["residence_universe"]
    print(f"  residence codes (unweighted): origin 1997 ages 12-16 {ru['origin_1997_ages_12_16']}; adult person-rounds 25-34 "
          f"{ru['adult_person_rounds_25_34']}; adult city/suburb rows by CV_URBAN-RURAL {ru['adult_city_suburb_rows_by_urban_rural']}")
    print(f"  transition 2x2 sample: {cnt['trans_2x2']['persons']} persons, {cnt['trans_2x2']['person_years']} person-years; "
          f"with a valid score {cnt['trans_2x2_with_afqt']} persons")
    for t, v in cnt["trans_2x2_by_afqt_tercile"].items():
        print(f"    AFQT tercile {t}: {v['persons']} persons, {v['person_years']} person-years")
    print(f"  origin-age respondents by own-AFQT tercile: {cnt['origin_afqt_tercile_sizes']}")
    print(f"  adult residence rows (25-34, city or suburb, no origin requirement): {cnt['adult_rows_city_or_suburb']['persons']} persons, "
          f"{cnt['adult_rows_city_or_suburb']['person_years']} person-years")
    print(f"  person-rounds at 25-34: {cnt['person_rounds_25_34']} (persons observed {cnt['persons_observed_25_34']}; "
          f"men {cnt['gender']['men_observed']}, women {cnt['gender']['women_observed']}); by year {cnt['person_rounds_by_year']}")
    for k in ("rounds_with_main_job_occupation", "rounds_valid_civilian_occupation", "rounds_military_occupation_code_not_counted",
              "rounds_employed_valid_occupation", "rounds_missing_current_flag_with_occupation", "rounds_k12_any",
              "rounds_k12_employed"):
        print(f"  {k}: {cnt[k]}")
    print(f"  armed forces: {cnt['armed_forces']}")
    print(f"  persons with a main-job occupation: {cnt['persons_with_occupation']}")
    print(f"  schooling: {cnt['persons_with_schooling']}; persons whose BA+ (degree) and 16+ years differ: "
          f"{cnt['schooling_definitions_discordant_persons']}")
    for d, v in cnt["schooling_samples"].items():
        print(f"  schooling sample {d}: {v['members']} persons with a valid score; by own-AFQT tercile (whole-sample cuts) {v['by_afqt_tercile_full']}; "
              f"highest grade >= 9: {v['members_grade_ge_9']} (by tercile cut over that sample {v['grade_ge_9_by_afqt_tercile_cut_over_that_sample']}); "
              f"grade below 9 or unknown: {v['members_grade_below_9_or_unknown']}")
    for nm in ("c02-emp-all-maj", "c02-any-all-maj", "c02-emp-all-ever", "c02-pop-all-maj"):
        g = cnt["teacher_groups"][nm]
        print(f"  {nm}: " + "; ".join(f"{k}: {v['teachers']} teachers vs {v['non_teachers']} non-teachers" for k, v in g.items()))


# =============================================================================
# Methods file (definitions, samples, N, departures: no values)
# =============================================================================

def write_methods(path: Path, cnt: dict[str, Any], meta: dict[str, Any], specs: list[Spec97],
                  check_lines: list[str] | None = None) -> None:
    L: list[str] = []
    a = L.append
    T = cnt["trans_2x2"]
    ru, af, dz, lk = cnt["residence_universe"], cnt["armed_forces"], cnt["design"], cnt["link_checks"]
    org, adu, urb = ru["origin_1997_ages_12_16"], ru["adult_person_rounds_25_34"], ru["adult_city_suburb_rows_by_urban_rural"]
    org_n, adu_n = sum(org.values()), sum(adu.values())
    chk = af["cv_doi_employed_check"]
    disc = chk["armed_not_doi_military"] + chk["doi_military_not_armed"]
    a("# NLSY97 hold-out: methods (no values)\n")
    a(f"Generated {meta['generated']} by `nlsy97_holdout.py` (re-sealed after an independent review: section 10). **This file contains definitions, "
      "samples and sample sizes only.** The NLSY97 counterparts of the NLSY79 residential-composition, teacher-selection and schooling moments are in "
      "`estimates/settings/sealed/nlsy97_holdout.json` (git-ignored until unsealed), written with `json.dump` and neither printed nor logged when it was "
      "produced; since then it has been opened only by programs that print structure and PASS/FAIL (`--verify-seal`, `nlsy97_holdout_xcheck.py`). "
      "It is to be opened by a person only after the next model specification is chosen. Section 10 lists every place where a count that bears on a sealed "
      "statistic was printed or listed.\n")
    a("## 1. What is sealed\n")
    a("All objects are NLSY97 analogues of NLSY79 moments already in the calibration log. The last column gives the key in the sealed json; the third the "
      "NLSY79 key (file in `estimates/settings/`) that each headline statistic matches.\n")
    a("| # | Object | NLSY79 counterpart (file, key) | Key in the sealed json |\n|---|---|---|---|")
    a("| 1 | Origin central-city share by tercile of own ability, T3 - T1 | `nlsy_composition.json`, `trans\\|afqt\\|origin_share_city\\|{1,2,3,diff}` | `trans\\|afqt\\|origin_share_city\\|{1,2,3,diff}` |")
    a("| 2 | City to suburb and suburb to city rates by own-ability tercile | `nlsy_composition.json`, `trans\\|afqt\\|{city_to_suburb,suburb_to_city}\\|{1,2,3,diff}` | same keys |")
    a("| 3 | Move rate by own-ability tercile | `nlsy_composition.json`, `trans\\|afqt\\|move_rate\\|{1,2,3,diff}` | same key |")
    a("| 4 | K-12 teachers at 25-34: gap in mean latent ability against employed non-teachers (A), against non-teacher college workers (B), BA+ teachers against BA+ non-teachers (C) | `nlsy_teachers.json`, `c80-emp-all-maj\\|{pooled,men,women}\\|{nt,col,cc}\\|gap` (college = 16+ years) | `c02-emp-all-maj\\|{pooled,men,women}\\|{nt,col_c16,cc_c16}\\|gap`; by degree `col_ba`, `cc_ba` |")
    a("| 5a | **Headline.** Share with 16+ years of schooling by own-ability tercile, highest grade >= 9, terciles cut over that sample | `education_occupation.json` (`nlsy/stats`), `terc\\|all_hgc9\\|{pooled,men,women}\\|ba_share\\|{1,2,3,diff}` | `edu9_cutm\\|c16\\|share\\|{1,2,3,diff}` (pooled), `edu9_cutm\\|c16\\|{men,women}\\|share\\|…` |")
    a("| 5b | **Headline.** Latent-ability gap, 16+ years vs less, highest grade >= 9 | `education_occupation.json`, `all_hgc9\\|{pooled,men,women}\\|gap` | `edu9\\|c16\\|{pooled,men,women}\\|gap` |")
    a("| 5c | The same with BA+ by degree (the ACS-style definition) | none (the NLSY79 extract has no degree variable) | `edu9_cutm\\|ba\\|…`, `edu9\\|ba\\|…` |")
    a("| 5d | Whole-sample tercile cuts, no grade rule (the previous headline: the `afqt_full` convention of `nlsy_composition.py`) | none: `education_occupation.py` cuts the terciles over its analysis sample | `edu\\|{ba,c16}\\|share\\|…`, `edu\\|{ba,c16}\\|{g}\\|gap`, `edu_cutm\\|…`, `edu9\\|{ba,c16}\\|share\\|…` |")
    a("| extra | Adult (25-34) suburb share by own-ability tercile | `nlsy_composition.json`, `adult\\|afqt\\|share\\|{1,2,3,diff}` | `adult\\|afqt\\|share\\|{1,2,3,diff}`, `adult_2x2\\|afqt\\|…`, `adult_urban\\|afqt\\|…` |")
    a("\nThe origin central-city share by parental-income tercile is already known (`nlsy_extensions.md`) and is not recomputed. Teaching shares at 25-34 "
      "(`share|…`, `pshare|…`), the adult suburb share and many variants are included as checks. Every statistic carries two standard errors (section 4).\n")
    a("## 2. Data\n")
    a("`raw/nlsy97_ability.parquet`, cut from the cached public file `raw/_downloads/nlsy97_all_1997-2023.zip` (8,984 respondents, rounds 1997-2023) by "
      "`fetch_data.py --sources nlsy --nlsy-cohorts nlsy97`. **Extract change** (in two additive steps). `NLSY_BLOCKS['nlsy97']` gained `household` "
      "(`YOUTH_SIBID01-04.01`, the round-1 roster IDs of the other NLSY97 respondents in the youth's household) and `jobs` (`CV_MAINJOB_FLG`, "
      "`YEMP_OCCODE-2002.02-14`, `YEMP_CURFLAG.01-15`, `CV_DOI_EMPLOYED`; `YEMP_OCCODE-2002.01` was already cut by the earlier `work` block): 410 more variables, "
      "911 in all. In the review round it gained `design` (`VSTRAT` R14897.00 and `VPSU` R14898.00, 1997) and `military` (`YEMP_MILFLAG.01-15`, 196 variables over rounds "
      "1998-2023): 198 more, 1,109 in all. Every earlier block is cut first and in the same order, so each new file is a strict superset: every earlier column and every earlier "
      "row of `metadata/nlsy97_ability_vars.csv` is bit-for-bit unchanged (verified for the last re-cut against the 911-variable file and against the original 501-variable file; "
      "the re-cut was written to a scratch root and installed only after the check; the previous files are kept in `raw/_nlsy97_extract_backup_20261005/` (501 variables) and "
      "`raw/_nlsy97_extract_backup_20261005b/` (911 variables)). The README's count of NLSY97 variables (340) is out of date and is not edited. No NLSY79 output was touched.\n")
    a("## 3. Definitions\n")
    a("**Ability (own AFQT).** The NLSY79 `afqt_c` is the equal-weight mean of the four age-normed ASVAB IRT z-scores (AR, WK, PC, MK), each standardized and the "
      "mean re-standardized, 1979 weights. The NLSY97 counterpart uses the four CAT-ASVAB ability estimates (`ASVAB_{AR,WK,PC,MK}_ABILITY_EST_POS/NEG`, 1997-98), "
      "age-normed exactly as `nlsy_ability.build_n97` does for block V (weighted rank-normal within birth-year x quarter cells, round-1 `SAMPLING_WEIGHT_CC`), "
      "then the equal-weight mean of the four standardized z-scores, re-standardized with the draw weights (`nlsy_ability.pop_afqt`). "
      f"Sample: respondents with all four scores and a positive weight (N = {cnt['four_casasvab_scores_and_weight']}, the figure behind the published R_comp_97). "
      "**Latent units.** A difference of group means is divided by sqrt(R_comp_97), the one-factor ULS reliability of the composite, re-estimated in every draw "
      "(the NLSY79 gaps use sqrt(R_comp_pop)). Terciles: weighted quantiles (`nlsy_extensions.tercile_of`).\n")
    a(f"**Items 1-3 (residential transitions).** The T3 NLSY97 machinery of `nlsy_extensions.py` (`build_trans97`, `stats_from`) through `nlsy_composition.trans_stats`. "
      f"Origin: `CV_MSA` at the 1997 interview for respondents aged 12-16 (central city = city, in MSA but not central city = suburb; not in MSA and 'not known' are "
      f"outside the 2x2). Adult: every round at ages 25-34 (`CV_AGE_INT_DATE`; rounds 2004-2019) with a classified `CV_MSA`, which from 2004 codes residence in a CBSA "
      f"(metropolitan plus micropolitan; see the departures table); each person carries the mean of his or her "
      f"round weights. The 2x2 sample is persons classified city or suburb at both dates ({T['persons']} persons, {T['person_years']} person-years; the same as the published NLSY97 table). "
      f"Terciles of own ability are cut over respondents aged 12-16 at the 1997 interview with a valid score (the analogue of `afqt_origin`), weights `w0` x multiplicity; "
      f"statistics by tercile: city to suburb, suburb to city, move rate, origin city share (all within the 2x2 sample), and top minus bottom. Variants (`transv|…`): origin at age 12 "
      f"(parent and youth reports), same-region movers only, adult windows 25-27, 28-32, rounds 2004-2011, 2013-2019, the *latest* classified round at age 25 or older "
      f"(any age up to 43; the 'latest available' reading of the request), adult rounds with `CV_URBAN-RURAL` = urban only (`adult_urban`, a partial sensitivity to the change of "
      f"universe: urban areas include micropolitan cores, so it does not isolate metropolitan residence), the tercile assignment by NLS's own percentile "
      f"`ASVAB_MATH_VERBAL_SCORE_PCT`, and `origin_all|…` (origin city share among everyone with a known city/suburb origin, no adult requirement).\n")
    a("**Extra: adult suburb share.** Share of person-rounds in the suburb among person-rounds at 25-34 in a CBSA with a known central-city status, by tercile of own ability cut over the "
      "whole sample with a valid score (`afqt_full`), `nlsy_composition.adult_stats`; sample = every respondent observed at 25-34 (`adult`), the 2x2 transition sample (`adult_2x2`) or "
      "urban rounds only (`adult_urban`).\n")
    a("**Item 4 (teachers).** Occupation of the *main job* in each round at 25-34: `CV_MAINJOB_FLG` is the roster position of the current or most recent employer at the "
      "interview, `YEMP_OCCODE-2002.NN` the 2002 Census 4-digit code of roster job NN. K-12 teachers are codes 2300 (pre-K and kindergarten), 2310 (elementary and middle), "
      "2320 (secondary), 2330 (special education) and 2340 (other teachers and instructors), the project's ACS definition (`nlsy_teachers.TEACHER_CODES['c02']`; the labels "
      "are those of the 2002 Census listing in NLSY79 Codebook Supplement Attachment 3, `att300.pdf`, https://www.nlsinfo.org/sites/default/files/attachments/121217/att300.pdf; "
      "postsecondary teachers (2200) and teacher assistants (2540) are not counted). A round counts if it is not an armed-forces round (below) and the code is a civilian occupation "
      f"(10-9750); the {cnt['rounds_military_occupation_code_not_counted']} rounds with a 'military specific' code (9800-9840) are at non-military employers, are not counted as an occupation "
      "and stay in the denominators. Samples as in `nlsy_teachers.py`: `emp` = rounds in which the main job is held at the interview date "
      "(`YEMP_CURFLAG` of the main job = 1; the analogue of `ESR_COL` = 1), `any` = every round with a valid occupation (current or most recent job), `pop` = every respondent "
      "observed at 25-34, a person with no occupation counting as a non-teacher. A person is a K-12 teacher if more than half of his or her counted rounds are K-12 "
      "('maj'; 'half' counts ties as teachers, 'ever' = at least one round). Comparison groups A, B, C as in the NLSY79 report (B and C use college, below). Weights: round-1 weight "
      "(`SAMPLING_WEIGHT_CC` 1997), one observation per person. The specification grid is that of `nlsy_teachers.make_specs` (all/core definitions x three rules, by type, "
      "highest grade >= 9, cross-sectional sample only (`CV_SAMPLE_TYPE` = 1), unweighted, mean round weights, for `emp` and `any`; three rules for `pop`).\n")
    a(f"**Armed forces (the `ESR_COL` = 4 analogue).** `CV_MAINJOB_FLG` never points to a military employer ({af['rounds_main_job_is_a_military_employer']} of "
      f"{af['rounds_with_a_military_employer_on_the_roster']} person-rounds with a military employer on the roster), so active-duty rounds have no main job and leave `emp` and `any` by construction; "
      "the earlier occupation-code rule (9800-9840) did not identify them (it dropped rounds at non-military employers) and the rounds stayed in the `pop` sample, the teaching-share denominators and the "
      "schooling samples, which the NLSY79 drops through `ESR_COL` = 4. A round is now an armed-forces round if a military roster employer (`YEMP_MILFLAG` = 1) is held at the interview "
      "(`YEMP_CURFLAG` = 1) and no civilian employer is. "
      f"{af['rounds_armed_dropped']} person-rounds ({af['persons_with_an_armed_round']} persons) are dropped from every sample; {af['persons_observed_only_in_armed_rounds']} persons are observed only in such rounds. "
      f"A reserve or guard member who also holds a civilian job ({af['rounds_current_military_and_current_civilian_job']} rounds) is a civilian worker, as in the CPS and in NLS's own "
      f"`CV_DOI_EMPLOYED` (2 = 'military service (but no job)', 2015+ rounds). The rule and `CV_DOI_EMPLOYED` agree in {chk['rounds_compared'] - disc} of {chk['rounds_compared']} "
      f"person-rounds with a valid `CV_DOI_EMPLOYED` ({chk['armed_and_doi_military']} flagged by both).\n")
    a("**Item 5 (schooling).** Two college definitions: `c16` = 16 or more years of schooling (`CV_HGC_EVER_EDT`, maximum over the person's rounds at 25-34, as `nlsy_teachers.py` uses the highest `HGCREV`: "
      "the NLSY79 definition) and `ba` = bachelor's degree or higher (`CV_HIGHEST_DEGREE_EVER_EDT`, maximum over the rounds at 25-34, codes 4-7; the NLSY79 extract has no degree variable). "
      "Sample: persons observed at 25-34 (an armed-forces-only person excluded, as the NLSY79 `observed`) with the schooling variable known. **Headline universe = that of the NLSY79 objects of "
      "`education_occupation.py`** (`evaluate_ba`, sample `all_hgc9`, the ACS Table 2 universe): highest grade >= 9, terciles of ability cut over that sample itself with the draw weights "
      "(`edu9_cutm|c16|share|…`, pooled; men and women are shares within the pooled terciles), and the gap in mean latent ability between 16+ years and the rest of the same sample (`edu9|c16|{pooled,men,women}|gap`). "
      "Alternatives: the same with `ba` (the ACS-style definition); the earlier headline universe (`edu|…`: no grade rule, terciles cut over the whole sample with a valid score, the `afqt_full` convention of "
      f"`nlsy_composition.py`; {cnt['schooling_samples']['c16']['members_grade_below_9_or_unknown']} persons with a valid score and a highest grade below 9 (`c16`; for `ba`, below 9 or unknown: {cnt['schooling_samples']['ba']['members_grade_below_9_or_unknown']}) are in `edu` and out of `edu9`, which is where the two "
      "universes differ); `edu_cutm` (no grade rule, terciles over the analysis sample). Gaps are in latent SD units.\n")
    a("**Teaching shares (check).** Employed K-12 person-rounds over all person-rounds at 25-34 with highest grade >= 9 (home production included, armed-forces rounds excluded), "
      "`nlsy_teachers.share_stats`; compare with the ACS 2009-13 shares of the calibration (NLSY97 respondents are 25-34 in 2004-2019).\n")
    a("## 4. Inference\n")
    a(f"Two bootstraps, B = {meta['B']} each, on the same statistics; both re-standardize the composite, re-estimate R_comp_97 and recompute all tercile cuts in every draw, hold the age-norming of the raw "
      "CAT-ASVAB estimates fixed (as `nlsy_ability.py` block V), use the same multipliers for every statistic (common random numbers, so correlations across statistics are valid), and store the SD of the draws as the SE "
      "and the 2.5 and 97.5 percentiles as a 95% interval.\n")
    a(f"1. **Household clusters** (`se`, `ci95_lo`, `ci95_hi`, `n_draws`; seed {meta['seed']}): round-1 households resampled with replacement. NLSY97 households are the connected components of the "
      f"*reciprocated* round-1 roster links `YOUTH_SIBID01-04.01` ({cnt['households']} households for {cnt['respondents']} respondents; sizes {cnt['household_sizes']}; {lk['links']} directed links between two respondents, "
      f"{lk['reciprocated_links']} reciprocated and used, {lk['unreciprocated_links_ignored']} not reciprocated and ignored, none pointing outside the sample; {lk['self_links_ignored']} respondent listed as his or her own household member, ignored). NLS documents 6,819 unique households, 1,862 of them with more than one "
      f"respondent (https://www.nlsinfo.org/content/cohorts/nlsy97/intro-to-the-sample/sample-design-screening-process); the {sum(v for k, v in cnt['household_sizes'].items() if k > 1)} multi-respondent components here match that "
      f"figure, so the difference of {6819 - cnt['households']} households is in single-respondent households (two respondents who sit in a multi-respondent component here would be alone in NLS's count; "
      f"the public file does not resolve which). Two of the 11 unreciprocated links join respondents in different design PSUs, so using them "
      f"(the earlier version did) merged two households that are not one and broke the nesting of households in PSUs; no household spans a PSU now ({dz['households_spanning_psu']}). This is the analogue of the NLSY79 `HHID` "
      "cluster bootstrap and is kept for comparability with the NLSY79 SEs.\n")
    a(f"2. **Design-based PSU bootstrap** (`se_psu`, `ci95_lo_psu`, `ci95_hi_psu`, `n_draws_psu`; seed {meta['seed']} + {PSU_SEED_OFFSET}). The NLSY97 is an area sample, and the public file carries "
      f"`VSTRAT` and `VPSU` (NLS: 'variable for use with variance PSU to correct for clustering in the sample design'; {dz['strata']} strata, two PSUs each, {dz['psus']} PSUs; respondents per PSU "
      f"{dz['respondents_per_psu']['min']} to {dz['respondents_per_psu']['max']}, median {dz['respondents_per_psu']['median']:g}). City/suburb residence is strongly clustered by PSU, so household-cluster SEs "
      "understate the sampling error of the residential statistics (section 7 reports the ratio of the two SEs for the published parental-income table, on which no ability statistic is involved). Each draw is a stratified Rao-Wu (1988) rescaling bootstrap: in every stratum one of the two PSUs is drawn "
      "(n_h - 1 = 1 draw) and its respondents' weights are multiplied by n_h/(n_h - 1) = 2, the other PSU's by 0. Households nest in PSUs, so they need no separate resampling; for a weighted total this reproduces the "
      "standard two-PSU-per-stratum variance estimator. **Use `se_psu` for any test of the NLSY97 against the model**; `se` is for comparison with the NLSY79 SEs, which the NLSY79 design (household clusters) allows.\n")
    a("## 5. Sample sizes (unweighted)\n")
    rows = [["respondents", str(cnt["respondents"])],
            ["households (reciprocated links) / variance strata x PSUs", f"{cnt['households']} / {dz['strata']} x 2 = {dz['psus']}"],
            ["all four CAT-ASVAB scores and positive weight", str(cnt["four_casasvab_scores_and_weight"])],
            ["aged 12-16 at the 1997 interview", str(cnt["origin_age_12_16"])],
            ["... with a valid score", str(cnt["origin_age_with_afqt"])],
            ["classified city/suburb/non-MSA at the origin (ages 12-16)", str(cnt["origin_classified_in_origin_age"])],
            ["2x2 transition sample: persons / person-years", f"{T['persons']} / {T['person_years']}"],
            ["... of which with a valid score", str(cnt["trans_2x2_with_afqt"])],
            ["adult (25-34) rows in a city or suburb, no origin requirement: persons / person-years", f"{cnt['adult_rows_city_or_suburb']['persons']} / {cnt['adult_rows_city_or_suburb']['person_years']}"],
            ["person-rounds at 25-34 / persons observed in any round", f"{cnt['person_rounds_25_34']} / {cnt['persons_observed_25_34_any_round']}"],
            ["armed-forces rounds dropped / persons observed after dropping them (the `pop` teacher sample)", f"{af['rounds_armed_dropped']} / {cnt['persons_observed_25_34']}"],
            ["rounds with a valid civilian main-job occupation / employed at the interview", f"{cnt['rounds_valid_civilian_occupation']} / {cnt['rounds_employed_valid_occupation']}"],
            ["K-12 rounds: any / employed", f"{cnt['rounds_k12_any']} / {cnt['rounds_k12_employed']}"],
            ["rounds with a military-specific occupation code (not counted as an occupation, not dropped)", str(cnt["rounds_military_occupation_code_not_counted"])],
            ["persons with a known degree / grade at 25-34", f"{cnt['persons_with_schooling']['degree']} / {cnt['persons_with_schooling']['hgc']}"],
            ["persons whose BA+ (degree) and 16+ years classifications differ", str(cnt["schooling_definitions_discordant_persons"])]]
    a(ne.tab(rows, ["sample", "N"]))
    a("Transition sample by tercile of own ability (persons / person-years): " + "; ".join(
        f"T{k}: {v['persons']} / {v['person_years']}" for k, v in cnt["trans_2x2_by_afqt_tercile"].items()) + ".\n")
    a("Residence codes (unweighted counts of `CV_MSA`, persons at the origin and person-rounds at 25-34; the origin is the 1997 interview at ages 12-16 on the 1990 MSA standard, the adult rounds are on CBSA standards; "
      "the NLS codebook lists the whole-sample frequencies of the same codes). Coverage: not in an MSA / CBSA is "
      f"{100 * org['not_in_msa_or_cbsa'] / org_n:.1f}% of the origin sample and {100 * adu['not_in_msa_or_cbsa'] / adu_n:.1f}% of the adult person-rounds, status unknown "
      f"{100 * org['in_msa_or_cbsa_status_unknown'] / org_n:.1f}% and {100 * adu['in_msa_or_cbsa_status_unknown'] / adu_n:.1f}%.\n")
    lab_res = {"not_in_msa_or_cbsa": "1 not in MSA (origin) / not in CBSA (adult)", "suburb": "2 in MSA/CBSA, not in central city (suburb)",
               "central_city": "3 in MSA/CBSA, in central city", "in_msa_or_cbsa_status_unknown": "4 in MSA/CBSA, central-city status not known",
               "not_in_country": "5 not in country"}
    rows = [[lab_res[k], str(org[k]), str(adu[k])] for k in org]
    a(ne.tab(rows, ["CV_MSA code", "origin N", "adult person-rounds N"]))
    a(f"Among the {urb['rural'] + urb['urban'] + urb['unknown'] + urb['missing']} adult city/suburb person-rounds, `CV_URBAN-RURAL` is urban in {urb['urban']}, rural in {urb['rural']} and unknown/missing in "
      f"{urb['unknown'] + urb['missing']}.\n")
    rows = []
    for nm, lab in (("c02-emp-all-maj", "employed, all K-12, majority (main)"), ("c02-any-all-maj", "any occupation, all K-12, majority"),
                    ("c02-emp-all-ever", "employed, all K-12, ever"), ("c02-pop-all-maj", "all observed, all K-12, majority")):
        for g, v in cnt["teacher_groups"][nm].items():
            rows.append([lab if g == "pooled" else "", g, str(v["teachers"]), str(v["non_teachers"])])
    a("Teacher and comparison groups (persons with a valid score; group A):\n")
    a(ne.tab(rows, ["specification", "group", "teachers", "non-teachers"]))
    a("The sizes of the BA+ comparison groups (B: non-teacher BA+ workers; C: BA+ teachers and BA+ non-teachers), under both college definitions and for "
      "every specification, are in the sealed file's `counts`, not here: their ratios to the group sizes above are unweighted versions of sealed "
      "statistics. (The ratio of the group sizes in this table is itself an unweighted version of the sealed teaching share, a check statistic; see section 10.)\n")
    a("Schooling samples (persons with a valid score): " + "; ".join(
        f"{d}: {v['members']} persons (by whole-sample ability tercile {', '.join(str(x) for x in v['by_afqt_tercile_full'].values())}); highest grade >= 9: {v['members_grade_ge_9']} "
        f"(by tercile cut over that sample {', '.join(str(x) for x in v['grade_ge_9_by_afqt_tercile_cut_over_that_sample'].values())}); highest grade below 9 or unknown: {v['members_grade_below_9_or_unknown']}"
        for d, v in cnt["schooling_samples"].items()) + ".\n")
    a("## 6. Departures from the NLSY79 definitions\n")
    dep = [
        ["Ability test", "Paper ASVAB 1980 at 15-23, NLS age-normed IRT z-scores; 1979 weights", "CAT-ASVAB 1997-98 at 12-18, ability estimates age-normed here (birth-year x quarter rank-normal, 1997 weights, held fixed in the bootstrap); NLS built its own percentile with custom weights for the 7,093 with all scores, not public (`ASVAB_MATH_VERBAL_SCORE_PCT` is a tercile-assignment variant)"],
        ["Timing of the score", "Measured in 1980, after the 1979 origin residence", "Measured within a year of the 1997 origin interview, before any move"],
        ["Population", "NLSY79 population samples (SAMPLE_ID 1-8, 10, 11, 13, 14)", "Cross-sectional sample and oversample together, round-1 weights (`cs` variant: cross-sectional sample only)"],
        ["Origin residence", "`SMSARES` 1979, ages 14-16; 'central city not known' is 29% of the SMSA residents aged 14-16 (20.6% of all 14-16-year-olds; `nlsy_extensions.json`, T3 `nlsy79_origin_coverage`)", "`CV_MSA` 1997, ages 12-16 at the interview; 'not known' is about 1% of the origin sample (1.3% of MSA residents; `nlsy97_origin_coverage`); origin on the pre-2004 (unlabelled) MSA standard, metropolitan areas only"],
        ["Adult residence", "`SMSARES` rounds 1982-98 at 25-34, a metropolitan (SMSA) universe at both dates; 33-41% of SMSA residents 'not known', dropped",
         f"`CV_MSA` rounds 2004-2019 at 25-34 codes residence in a **CBSA, metropolitan plus micropolitan** (2000 standards to 2011, 2010 standards from 2013; the codebook: 'MSA-equivalent (CBSA)'), whereas the 1997 origin is a metropolitan MSA. "
         f"Not being in a CBSA is {100 * adu['not_in_msa_or_cbsa'] / adu_n:.1f}% of the adult person-rounds against {100 * org['not_in_msa_or_cbsa'] / org_n:.1f}% not in an MSA at the origin (19.6% weighted, `nlsy97_origin_coverage`). "
         "Adults in micropolitan areas therefore enter the 2x2, their principal cities counted as 'central city'; the NLSY79 adult universe is metropolitan only. A move from a metropolitan to a micropolitan area can count as a move or as a stay, "
         "so the move rates and the adult suburb share have no signed bias from this change of universe, and neither do their by-tercile gradients. Mechanical reclassification at the changes of MSA standard (T3 caveat) raises the *level* of the move rates; "
         "cross-MSA moves cannot be separated. `adult_urban` (urban rounds only) is a partial sensitivity"],
        ["Cluster and design", "`HHID` household clusters", f"Household bootstrap on the connected components of the reciprocated round-1 roster links `YOUTH_SIBID01-04.01` ({cnt['households']} households; siblings living together in 1997 only; NLS documents 6,819), `se`; and the design-based Rao-Wu bootstrap over {dz['strata']} variance strata x 2 PSUs, `se_psu`. The NLSY79 public file has no PSU or stratum variable (its only design-related variables are `HHID` and `SAMPLE_ID`), so the NLSY79 SEs can only be household-clustered, and its residential SEs may be understated in the same way"],
        ["Occupation", "CPS job, 1980 Census 3-digit codes (`CPSOCC80`), 1982-98", "Main job = roster job `CV_MAINJOB_FLG` (current or most recent employer), 2002 Census 4-digit codes (`YEMP_OCCODE-2002`), 2004-2019; K-12 codes 2300-2340; the main job can be a job held since the last interview that has ended"],
        ["Employed at the interview", "`ESR_COL` = 1", "`YEMP_CURFLAG` of the main job = 1 (checked against `CV_DOI_EMPLOYED` in 2015+ rounds)"],
        ["Armed forces", "Rounds with `ESR_COL` = 4 (armed forces at the interview) dropped from every sample", f"Rounds in which a military roster employer (`YEMP_MILFLAG` = 1) is held at the interview and no civilian employer is ({af['rounds_armed_dropped']} rounds) dropped from every sample; `CV_MAINJOB_FLG` never names a military employer, so the rule uses the roster; reserve or guard members with a civilian job stay civilian ({af['rounds_current_military_and_current_civilian_job']} rounds); agrees with `CV_DOI_EMPLOYED` = 2 in {chk['rounds_compared'] - disc} of {chk['rounds_compared']} 2015+ rounds"],
        ["College", "16+ years of schooling (highest `HGCREV` at 25-34); in `education_occupation.py` the BA+ objects use highest grade >= 9 and terciles cut over the sample", "Headline: 16+ years (`CV_HGC_EVER_EDT`), grade >= 9, terciles cut over the sample (`edu9_cutm|c16`, `edu9|c16`); BA+ by degree reported alongside (`ba`); the whole-sample cuts without a grade rule are alternatives (`edu`)"],
        ["Cohort and period", "Born 1957-64, 25-34 in 1982-98", "Born 1980-84, 25-34 in 2004-2019, a period that includes the ACS 2009-13 window of the occupational and teaching-share targets"],
        ["Round weights", "1979 weights for person-level statistics, round weights for transitions", "Same: round-1 weights for person-level statistics, round weights for transitions"],
    ]
    a(ne.tab(dep, ["Item", "NLSY79 (nlsy_composition / nlsy_teachers / education_occupation)", "NLSY97 (this script)"]))
    a("## 7. Checks run before sealing (no sealed value involved)\n")
    a("`python3 nlsy97_holdout.py --check` prints the sample sizes and PASS/FAIL for: the published `R_comp_97` and the published NLSY97 parental-income table "
      "(`nlsy_extensions.json`, persons as clusters, seed 20260929) rebuilt from this extract; the adult rows of `build_trans97_adult` equal those of `ne.build_trans97`; "
      "the design variables (117 strata x 2 PSUs, households nested in PSUs, Rao-Wu multipliers with mean one); "
      "the generic gap, teaching-share and schooling functions reproduce `nlsy_teachers.evaluate` / `share_stats` and the published `education_occupation.json` BA-by-tercile and gap objects on the NLSY79 data (at the sample and, for the teacher functions, at one bootstrap draw) and the "
      "published `nlsy_teachers.json`; `terciles()` reproduces `nlsy_composition.groupings`; the transition statistics reproduce the published `nlsy_composition.json`; and "
      "structural properties of the NLSY97 statistics (finite, shares in [0, 1], ordered cuts, tercile masses). `--selftest` runs synthetic unit tests (including the Rao-Wu variance identity, the household "
      "components and the armed-forces rule), and `--verify-seal` "
      "prints the structure (not the values) of the sealed file. `--seal` runs the same checks first and refuses to seal if one fails; it prints one line, "
      "`sealed: <path> (<k> moments)`, and refuses to overwrite an existing sealed file without `--force`. **Scope of the cross-check.** The checks above validate the generic code on "
      "NLSY79 data. `nlsy97_holdout_xcheck.py` recomputes 66 of the sealed statistics (transitions by tercile, the teacher gaps A, B and C of the main specification for pooled, men and women, the BA+ / 16+ share by tercile and "
      "the BA+ vs non-BA gap in both universes, R_comp_97) in separately written pandas code at the sample (against the sealed point estimates), at the first draw of each bootstrap (against `evaluate97`) and over both whole "
      "bootstraps (against the sealed SEs and 95% limits), printing only PASS/FAIL and the maximum absolute difference. In this version it also derives the household clusters (breadth-first search on the reciprocated links), "
      "the PSU index and the armed-forces rounds independently from the raw columns. It does not check what the raw variables *mean* (the codebook readings: CBSA vs MSA, `CV_MAINJOB_FLG`, `YEMP_CURFLAG`, "
      "`CV_HGC_EVER_EDT`), which the earlier version, reusing `Frame97`'s household index and reading the variables the same way, could not either; those readings are documented in this file and were reviewed against the codebook.\n")
    if check_lines:
        a("Result at sealing (published and NLSY79 quantities only; the NLSY97 statistics enter only through counts and PASS/FAIL):\n")
        a("```\n" + "\n".join(check_lines) + "\n```\n")
    a("## 8. The sealed file\n")
    a("`estimates/settings/sealed/nlsy97_holdout.json`: `meta`, `definitions`, `counts` (the sample sizes above), `statistics` "
      "(key -> `value`, `se`, `ci95_lo`, `ci95_hi`, `n_draws` for the household bootstrap and `se_psu`, `ci95_lo_psu`, `ci95_hi_psu`, `n_draws_psu` for the design-based bootstrap) and `headline` (the five requested objects, with values and both SEs). "
      "Keys are `<block>|…|<stat>`: `trans|afqt|<stat>|{1,2,3,diff}`, "
      "`transv|<variant>|afqt|…`, `origin_all|afqt|share|…`, `<spec>|<group>|<comparison>|{mean_t,mean_c,gap,gap_obs}` with spec = "
      "`c02-<sample>-<definition>-<rule>[-hgc9|-cs|-unit|-round]`, `edu|{ba,c16}|…`, `edu9|…`, `edu_cutm|…`, `edu9_cutm|{ba,c16}|[men|women|]share|…`, `edu9_cutm|…|cut|…`, `adult|…`, `adult_2x2|…`, `adult_urban|…`, `share|{acs,all}|…`, `pshare|…`, `cut|…`, `R_comp_97*`. Latent-SD units divide by sqrt(R_comp_97); "
      "`gap_obs` is in standardized-composite units.\n")
    a("## 9. Caveats\n")
    for ln in [
        "The NLSY97 occupation is the main job of the round, which is not the NLSY79 CPS item: a person with two jobs is classified by the roster position NLS flags as main, and an "
        "'any' round can be a job that ended before the interview.",
        "Respondents' later rounds (after 34) are not used for teachers or schooling; the youngest respondents contribute fewer rounds at 25-34 than the oldest.",
        "Residential transitions: the adult rounds are classified by CBSA (metropolitan plus micropolitan) and the origin by metropolitan MSA, so the 2x2 universe differs between the two dates and the NLSY79 adult universe is "
        "metropolitan only (departures table). The move rates and the adult suburb share can be biased in either direction by this, and so can their by-tercile gradients; only the mechanical reclassification at the MSA-standard "
        "changes has a known sign (it raises move-rate levels). Cross-MSA moves cannot be separated. The origin-to-adult comparison is more favorable to the NLSY79 rate in cleanliness and to the NLSY97 rate in sample size.",
        "`se` (household clusters) understates the sampling error of clustered outcomes, above all residence; the design-based `se_psu` is the SE to use against the model. `se_psu` is conditional on NLS's variance PSUs, which are "
        "constructed for variance estimation (two per stratum) and are not the 147 sampled PSUs; both SEs treat the age-norming of the raw CAT-ASVAB estimates as known.",
        "Respondents without all four CAT-ASVAB scores (about a fifth) are outside every statistic that uses ability, and no correction is made for them or for attrition "
        "(round-1 weights, as `SAMPWEIGHT` 1979 in the NLSY79 analyses).",
        "The NLSY79 public file has no PSU or stratum variable (its only design-related variables are `HHID` and `SAMPLE_ID`), so the NLSY79 SEs of the calibration targets can only be household-clustered; "
        "for residential statistics they may be understated in the same way as the household-cluster SEs here.",
        "With half of the PSUs given weight 0 in every design draw, thin cells (men in pre-K or special-education teaching) are empty in more draws than in the household bootstrap: see `n_draws_psu`.",
        "This file reports no value; the checks above use only published or unsealed NLSY79 and NLSY97 results.",
    ]:
        a(f"- {ln}")
    a("\n## 10. Changes after the independent review, and disclosure log\n")
    a("The first sealed file was reviewed independently; this version fixes what the review found. Changes: (1) the design-based bootstrap (`se_psu`) was added to every statistic, and the households were rebuilt from reciprocated "
      "links; (2) the item-5 headline moved to the NLSY79 universe (`edu9_cutm|c16`, `edu9|c16`; the earlier whole-sample-cut statistics are kept as `edu`); (3) the residence-universe departure (CBSA vs MSA) is stated and an "
      "urban-only sensitivity added (`adult_urban`); (4) armed-forces rounds are identified from the roster, not from occupation codes (the earlier text described the 9800-9840 rule as the `ESR_COL` = 4 analogue, which it was not), "
      "and dropped from every sample; (5) wording: the NLSY79 'central city not known' share, the contents of the `jobs` block, the household count and the scope of the cross-check. The previous sealed file was not read by a "
      "person; a regression check of this version against it printed only counts (statistics compared, statistics whose point estimate differs, by block) and PASS/FAIL for the blocks the changes cannot affect "
      "(transitions, adult shares, tercile cuts, R_comp_97, the employed-sample teacher statistics), and the previous file was deleted afterwards. Nothing in this round printed a sealed value.\n")
    a("Disclosure log (counts that bear on sealed statistics): (a) the first `--check` run, before the code was changed, printed the numbers of BA+ and 16+ persons inside the teacher and comparison groups; their ratios to the group sizes "
      "approximate sealed `col_share_t_*`, `col_share_nt_*` (composition of the teacher and comparison groups: part of item 4, not a requested moment); that log was deleted, the printed and listed counts no longer contain them, "
      "and they are in the sealed file's `counts` only. (b) Section 5 lists the numbers of teachers and non-teachers (group A), whose ratio is an unweighted version of the sealed teaching `share` (a check statistic, "
      "not one of the five objects); they are the N behind the gaps and are kept. (c) Section 5 lists schooling sample sizes by ability tercile (the denominators of the sealed BA+ shares, not their numerators), the "
      "number of persons with a highest grade below 9 and the number of persons whose BA+ (degree) and 16+ years classifications differ. (d) The residence-code counts in section 5 are unweighted frequencies of the codebook categories "
      "(the NLS codebook lists the same frequencies for the whole sample), not statistics by ability; only the coverage shares (not in an MSA / CBSA, status unknown) are given as percentages, not the city/suburb split. "
      "(e) While the review's points were diagnosed, the terminal showed unweighted counts: the `CV_MSA` codes by round at ages 25-34 and `CV_URBAN-RURAL` among the city/suburb rows (so the overall adult city/suburb split by round is known, "
      "an unweighted and ability-free version of the overall level of the extra statistics `adult|afqt|overall`; the by-tercile statistics are not), the armed-forces round counts with their cross-tabulation against `CV_DOI_EMPLOYED`, "
      "the household and PSU structure, and the NLS codebook entries of the variables involved. None is by ability.\n")
    path.write_text("\n".join(L) + "\n")


# =============================================================================
# Self-test (synthetic)
# =============================================================================

def selftest() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'ok' if cond else 'FAIL'}] {name}")

    # households: 0-1-2 chained, 3-4 pair, 5 alone; link 6->7 not reciprocated (ignored: 6 and 7 are separate households)
    pid = np.array([10, 11, 12, 13, 14, 15, 16, 17], float)
    sib = np.full((8, 4), np.nan)
    sib[0, 0], sib[1, 0], sib[1, 1], sib[2, 0] = 11, 10, 12, 11
    sib[3, 0], sib[4, 0] = 14, 13
    sib[6, 0] = 17
    idx, info = household_clusters(pid, sib)
    check("household components (reciprocated links only)", idx[0] == idx[1] == idx[2] and idx[3] == idx[4] and len(set(idx.tolist())) == 5
          and idx[6] != idx[7] and idx[5] not in (idx[0], idx[3], idx[6], idx[7]))
    check("households are numbered by their lowest respondent index", idx.tolist() == [0, 0, 0, 1, 1, 2, 3, 4])
    check("link counts", info["links"] == 7 and info["reciprocated_links"] == 6 and info["unreciprocated_links_ignored"] == 1
          and info["links_not_in_sample"] == 0)
    sib[5, 0] = 9999
    check("sibling id outside the sample is counted, not linked", household_clusters(pid, sib)[1]["links_not_in_sample"] == 1)
    sib[5, 0] = 15
    idx2, info2 = household_clusters(pid, sib)
    check("a self-link is counted and ignored", info2["self_links_ignored"] == 1 and info2["links"] == 7 and idx2.tolist() == idx.tolist())
    sib[5, 0] = np.nan
    # design PSUs: strata 1 and 2, two PSUs each
    vs = np.array([1, 1, 1, 1, 2, 2, 2], float)
    vp = np.array([1, 2, 1, 2, 1, 2, 2], float)
    psu, pairs = design_psus(vs, vp)
    check("design PSU index", psu.tolist() == [0, 1, 0, 1, 2, 3, 3] and pairs.tolist() == [[0, 1], [2, 3]])
    try:
        design_psus(np.array([1.0, 1.0, 2.0]), np.array([1.0, 2.0, 1.0]))
        check("a stratum with one PSU is refused", False)
    except ValueError:
        check("a stratum with one PSU is refused", True)

    class _PickRng:                                        # a stand-in generator that returns a preset pick per stratum
        def __init__(self, pick: np.ndarray) -> None:
            self.pick = pick

        def integers(self, lo: int, hi: int, size: int) -> np.ndarray:
            return self.pick

    # Rao-Wu: one PSU per stratum with weight 2, the other 0; over ALL 2^H draws the mean multiplier of every PSU is 1 and the
    # variance of a weighted total is sum_h (t_h1 - t_h2)^2, the standard two-PSU-per-stratum variance estimator
    t_psu = np.array([3.0, 1.0, 5.0, 4.5])                # PSU totals of some variable (strata (0, 1) and (2, 3))
    est, mults = [], []
    for k0 in (0, 1):
        for k1 in (0, 1):
            mm = rao_wu_multipliers(pairs, 4, _PickRng(np.array([k0, k1])))
            mults.append(mm)
            est.append(float((mm * t_psu).sum()))
    check("Rao-Wu multipliers are 2 for one PSU per stratum and 0 for the other", all(sorted(m_[:2].tolist()) == [0.0, 2.0]
          and sorted(m_[2:].tolist()) == [0.0, 2.0] for m_ in mults))
    check("Rao-Wu multipliers have mean one over all draws", np.allclose(np.mean(mults, axis=0), 1.0))
    check("Rao-Wu total is unbiased with the textbook 2-PSU variance",
          np.isclose(np.mean(est), t_psu.sum()) and np.isclose(np.mean((np.array(est) - t_psu.sum()) ** 2),
                                                               (3.0 - 1.0) ** 2 + (5.0 - 4.5) ** 2))
    mm = rao_wu_multipliers(pairs, 4, np.random.default_rng(3))
    check("Rao-Wu draw with a real generator", mm.sum() == 4.0 and set(mm.tolist()) == {0.0, 2.0})
    # occupation validity
    x = np.array([np.nan, 0, 10, 2310, 9750, 9800, 9840, 9950])
    check("valid civilian codes", occ_valid(x).tolist() == [False, False, True, True, True, False, False, False])
    check("military codes", occ_military(x).tolist() == [False, False, False, False, False, True, True, False])
    check("K-12 codes are the project's ACS definition", [c for t in TYPES for c in K12_CODES[t]] == [2300, 2310, 2320, 2330, 2340])
    # edu_stats on a toy sample: terciles 0,1,2 with 4 persons each, flag share 0.25, 0.5, 0.75; weights unit
    af = np.repeat([-1.0, 0.0, 1.0], 4)
    terc = np.repeat([0, 1, 2], 4)
    flag = np.array([1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0], bool)
    masks = {"pooled": np.ones(12, bool), "men": np.arange(12) % 2 == 0, "women": np.arange(12) % 2 == 1}
    r = edu_stats("e", af, 0.5, np.ones(12), terc, np.ones(12, bool), flag, masks)
    check("share by tercile", np.allclose([r["e|share|1"], r["e|share|2"], r["e|share|3"]], [0.25, 0.5, 0.75])
          and np.isclose(r["e|share|diff"], 0.5) and np.isclose(r["e|share|diff_21"], 0.25) and np.isclose(r["e|share|diff_32"], 0.25))
    mt, mc = af[flag].mean(), af[~flag].mean()
    check("gap in latent units", np.isclose(r["e|pooled|gap"], (mt - mc) / 0.5) and np.isclose(r["e|pooled|gap_obs"], mt - mc)
          and np.isclose(r["e|overall|pooled"], 0.5))
    # weights and multiplicities enter linearly: doubling one person's weight
    w = np.ones(12)
    w[0] = 3.0
    r2 = edu_stats("e", af, 0.5, w, terc, np.ones(12, bool), flag, masks)
    check("weighted share", np.isclose(r2["e|share|1"], 3 / 6))
    # gap_stats: teachers 0..3 (afqt 1.0, 1.0, 0.0, 0.0), non-teachers 4..9 (afqt 0, 0, 0, 0, -1, -1); persons 0, 2, 4, 6 are 'college'
    a = np.array([1.0, 1.0, 0.0, 0.0, 0, 0, 0, 0, -1, -1])
    tch = np.arange(10) < 4
    col = np.array([1, 0, 1, 0, 1, 0, 1, 0, 0, 0], bool)
    mk = {"pooled": np.ones(10, bool)}
    g = gap_stats("s", a, 1.0, np.ones(10), np.ones(10, bool), tch, ~tch, mk, {"ba": col})
    check("A gap: teachers 0.5 vs non-teachers -1/3", np.isclose(g["s|pooled|nt|gap"], 0.5 - (-2 / 6)))
    check("B gap: teachers 0.5 vs non-teacher college 0", np.isclose(g["s|pooled|col_ba|gap"], 0.5) and np.isclose(g["s|pooled|col_ba|mean_c"], 0.0))
    check("C gap: college teachers 0.5 vs non-teacher college 0", np.isclose(g["s|pooled|cc_ba|gap"], 0.5))
    check("shares", np.isclose(g["s|pooled|share"], 0.4) and np.isclose(g["s|pooled|col_share_t_ba"], 0.5)
          and np.isclose(g["s|pooled|col_share_nt_ba"], 2 / 6))
    check("latent scaling", np.isclose(gap_stats("s", a, 0.5, np.ones(10), np.ones(10, bool), tch, ~tch, mk, {"ba": col})["s|pooled|nt|gap"],
                                      2 * g["s|pooled|nt|gap"]))
    # classify (reused from nlsy_teachers): majority rule
    check("majority rule", nt.classify(np.array([0, 2, 2, 3]), np.array([0, 1, 2, 2]), "maj").tolist() == [False, False, True, True])
    # origin_all_stats
    orig = np.array([0, 0, 1, 1, 1, 2, -1, 0])
    tc = np.array([0, 0, 0, 1, 1, 1, 1, 2])
    o = origin_all_stats(orig, np.ones(8), np.ones(8), tc)
    check("origin city share among known origins", np.isclose(o["share|1"], 2 / 3) and np.isclose(o["share|2"], 0.0)
          and np.isclose(o["share|3"], 1.0) and np.isclose(o["share|diff"], 1 / 3))
    # terciles(): cut over the origin-age population vs the whole sample
    rng = np.random.default_rng(1)
    z = rng.normal(size=300)
    org = np.arange(300) < 150
    (to, co), (tf, cf) = terciles(z, np.ones(300), org, np.ones(300, bool))
    check("terciles are balanced in their own population", abs(np.mean(to[org] == 0) - 1 / 3) < 0.03 and abs(np.mean(tf == 0) - 1 / 3) < 0.03
          and co[0] < co[1] and cf[0] < cf[1])
    # specs
    names = [s.name for s in make_specs97()]
    check("spec names unique", len(names) == len(set(names)) and "c02-emp-all-maj" in names)
    # terc_share_stats: terciles cut once over the pooled sample, shares within tercile for each group
    af9 = np.repeat([-1.0, 0.0, 1.0], 4)
    tc9 = np.repeat([0, 1, 2], 4)
    fl9 = np.array([1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0], bool)
    mk9 = {"pooled": np.ones(12, bool), "men": np.arange(12) % 2 == 0, "women": np.arange(12) % 2 == 1}
    r9 = terc_share_stats("p", af9, np.ones(12), tc9, np.ones(12, bool), fl9, mk9)
    check("terc_share_stats pooled", np.allclose([r9["p|share|1"], r9["p|share|2"], r9["p|share|3"]], [0.25, 0.5, 0.75])
          and np.isclose(r9["p|share|diff"], 0.5))
    men1 = np.array([0, 2])                                # men are the even rows; tercile 1 holds rows 0-3, flag only on row 0
    check("terc_share_stats by group", np.isclose(r9["p|men|share|1"], 1 / 2) and np.isclose(r9["p|women|share|1"], 0.0)
          and np.isclose(r9["p|men|share|3"], 1.0) and men1.size == 2)
    # bounded keys: shares are, cut points are not
    check("bounded keys", is_bounded_key("edu9_cutm|c16|share|2") and is_bounded_key("edu9_cutm|c16|men|share|1")
          and not is_bounded_key("edu9_cutm|c16|cut|1") and is_bounded_key("adult_urban|afqt|share|3")
          and not is_bounded_key("edu9_cutm|c16|share|diff"))
    # armed-forces rule on a synthetic cohort: military job held and no civilian job held
    class _Fake:
        def __init__(self, data: dict, n: int, by_q: dict) -> None:
            self.data, self.by_q = data, by_q
            self.df = np.zeros((n, 1))

        def try_get(self, q: str, y: Any = None) -> np.ndarray | None:
            v = self.data.get((q, str(y)))
            return None if v is None else np.where(v < 0, np.nan, np.asarray(v, float))

        def get(self, q: str, y: Any = None) -> np.ndarray:
            v = self.try_get(q, y)
            if v is None:
                raise KeyError(q)
            return v

    nan = np.nan
    Y = "2005"
    dta = {("CV_AGE_INT_DATE", Y): np.full(6, 27.0), ("SAMPLING_WEIGHT_CC", Y): np.full(6, 100.0),
           ("CV_MAINJOB_FLG", Y): np.array([1, nan, 2, 2, 2, nan]),
           # roster job 1
           ("YEMP_MILFLAG.01", Y): np.array([0, 1, 1, 1, 1, nan]), ("YEMP_CURFLAG.01", Y): np.array([1, 1, 1, 0, 1, nan]),
           ("YEMP_OCCODE-2002.01", Y): np.array([2310, nan, nan, nan, nan, nan]),
           # roster job 2 (the civilian main job of persons 2, 3, 4)
           ("YEMP_MILFLAG.02", Y): np.array([nan, nan, 0, 0, 0, nan]), ("YEMP_CURFLAG.02", Y): np.array([nan, nan, 1, 0, 0, nan]),
           ("YEMP_OCCODE-2002.02", Y): np.array([nan, nan, 2320, 4000, 4000, nan])}
    fake = _Fake(dta, 6, {"CV_AGE_INT_DATE": [Y]})
    pr = build_pr97(fake)                                # type: ignore[arg-type]
    check("armed rule: civilian worker, active duty, reservist with a civilian job, ended military job, active duty with an "
          "ended civilian job, no job", pr.armed.tolist() == [False, True, False, False, True, False]
          and pr.milc.tolist() == [False, True, True, False, True, False] and pr.civc.tolist() == [True, False, True, False, False, False]
          and pr.main_mil.sum() == 0)
    check("main-job occupation and current flag follow CV_MAINJOB_FLG", pr.occ[0] == 2310 and pr.cur[0] == 1 and pr.occ[2] == 2320
          and pr.cur[2] == 1 and pr.occ[3] == 4000 and pr.cur[3] == 0 and np.isnan(pr.occ[1]))
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


# =============================================================================
# CLI
# =============================================================================

HEADLINE_FIELDS = ("value", "se", "ci95_lo", "ci95_hi", "se_psu", "ci95_lo_psu", "ci95_hi_psu")


def build_headline(S: dict[str, dict[str, float]]) -> dict[str, Any]:
    """The requested objects, with values.  `se` is the household-cluster bootstrap SE (the NLSY79 analogue);
    `se_psu` the design-based (Rao-Wu PSU) SE, the one that respects the sample design."""
    def pick(keys: dict[str, str]) -> dict[str, Any]:
        return {lab: {"key": k, **{f: S[k][f] for f in HEADLINE_FIELDS}} for lab, k in keys.items() if k in S}

    h: dict[str, Any] = {"se_note": "se = household-cluster bootstrap SE (comparable to the NLSY79 SEs); se_psu = design-based "
                                    "stratified Rao-Wu PSU bootstrap SE (VSTRAT x VPSU), the one that respects the NLSY97 sample design"}
    for item, stat in (("1_origin_city_share_by_afqt_tercile", "origin_share_city"), ("2a_city_to_suburb_by_afqt_tercile", "city_to_suburb"),
                       ("2b_suburb_to_city_by_afqt_tercile", "suburb_to_city"), ("3_move_rate_by_afqt_tercile", "move_rate")):
        h[item] = pick({"T1": f"trans|afqt|{stat}|1", "T2": f"trans|afqt|{stat}|2", "T3": f"trans|afqt|{stat}|3",
                        "T3_minus_T1": f"trans|afqt|{stat}|diff"})
    spec = "c02-emp-all-maj"
    t: dict[str, Any] = {}
    for g in GROUPS:
        for cl, lab in (("nt", "A_vs_employed_non_teachers"), ("col_ba", "B_vs_non_teacher_BA"), ("cc_ba", "C_BA_teachers_vs_BA_non_teachers"),
                        ("col_c16", "B_vs_non_teacher_16yrs"), ("cc_c16", "C_16yrs_teachers_vs_16yrs_non_teachers")):
            t[f"{g}|{lab}"] = pick({"gap_latent_sd": f"{spec}|{g}|{cl}|gap"})["gap_latent_sd"] if f"{spec}|{g}|{cl}|gap" in S else None
    h["4_teacher_latent_afqt_gap__" + spec] = t
    # item 5.  Main: the universe of the NLSY79 objects of education_occupation.evaluate_ba (sample `all_hgc9`: highest grade
    # >= 9, 16+ years of schooling, terciles cut over that sample).  Alternatives: BA+ by degree (the ACS-style definition),
    # and tercile cuts over the whole sample with a valid score without the grade rule (`afqt_full` of nlsy_composition)
    e: dict[str, Any] = {}
    for name, d, shp, gp, nlsy79 in (
            ("main_c16_grade_ge_9_terciles_cut_over_sample", "c16", "edu9_cutm", "edu9",
             "terc|all_hgc9|pooled|ba_share|{1,2,3,diff} and all_hgc9|pooled|gap (education_occupation.py)"),
            ("alt_ba_degree_grade_ge_9_terciles_cut_over_sample", "ba", "edu9_cutm", "edu9",
             "no NLSY79 counterpart (no degree variable); ACS-style BA+; same universe as the main"),
            ("alt_c16_whole_sample_terciles", "c16", "edu", "edu", "terc|all|pooled|ba_share|... without the grade rule, whole-sample cuts"),
            ("alt_ba_degree_whole_sample_terciles", "ba", "edu", "edu", "no NLSY79 counterpart")):
        e[name] = {"nlsy79_counterpart": nlsy79,
                   **pick({"T1": f"{shp}|{d}|share|1", "T2": f"{shp}|{d}|share|2", "T3": f"{shp}|{d}|share|3",
                           "T3_minus_T1": f"{shp}|{d}|share|diff", "gap_latent_sd": f"{gp}|{d}|pooled|gap",
                           "gap_latent_sd_men": f"{gp}|{d}|men|gap", "gap_latent_sd_women": f"{gp}|{d}|women|gap"})}
    h["5_BA_share_by_afqt_tercile_and_BA_vs_nonBA_gap"] = e
    h["6_extra_adult_suburb_share_by_afqt_tercile"] = pick({"T1": "adult|afqt|share|1", "T2": "adult|afqt|share|2", "T3": "adult|afqt|share|3",
                                                            "T3_minus_T1": "adult|afqt|share|diff", "overall": "adult|afqt|overall"})
    return h


def run_seal(fr: Frame97, B: int, seed: int, outdir: Path, quick: bool, check_lines: list[str] | None = None) -> tuple[Path, int]:
    res = bootstrap97(fr, B, seed)
    cnt = sample_counts(fr)
    S = res["stats"]
    meta = {"generated": dt.datetime.now().isoformat(timespec="seconds"), "script": "data/spatial/nlsy97_holdout.py",
            "sealed": True, "B": B, "seed": seed, "quick": bool(quick), "failed_draws": res["failed_draws"],
            "failed_draws_psu": res["failed_draws_psu"],
            "bootstrap": "household clusters: NLSY97 round-1 households (connected components of the reciprocated links "
                         "YOUTH_SIBID01-04.01), multinomial multiplicities, same draws for every statistic; composite re-standardized, "
                         "R_comp_97 re-estimated and all terciles recomputed in each draw; age-norming of the CAT-ASVAB estimates held fixed",
            "bootstrap_psu": f"design-based: stratified Rao-Wu rescaling bootstrap over the {fr.n_strata} variance strata x 2 variance PSUs "
                             f"(VSTRAT, VPSU): one PSU drawn per stratum with weight x 2, the other x 0; seed {seed} + {PSU_SEED_OFFSET}; "
                             "same draws for every statistic; everything else as above",
            "weights": "round-1 SAMPLING_WEIGHT_CC for person-level statistics, round weights for transitions",
            "samples": {"respondents": cnt["respondents"], "households": cnt["households"], "strata": fr.n_strata, "psus": fr.n_psu},
            "extract": {"file": "raw/nlsy97_ability.parquet", "respondents": int(fr.n), "variables": int(fr.c.df.shape[1]),
                        "blocks_added_for_holdout": ["household", "jobs", "design", "military"],
                        "sha256": hashlib.sha256((HERE / "raw" / "nlsy97_ability.parquet").read_bytes()).hexdigest()},
            "code_sha256": {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest()
                            for f in ("nlsy97_holdout.py", "nlsy_ability.py", "nlsy_composition.py", "nlsy_extensions.py", "nlsy_teachers.py")}}
    defs = {"college_defs": {"ba": "highest degree at 25-34 >= bachelor's (CV_HIGHEST_DEGREE_EVER_EDT >= 4)",
                             "c16": "highest grade at 25-34 >= 16 (CV_HGC_EVER_EDT)"},
            "comparisons": {"nt": "A: teachers vs all non-teachers", "col_<d>": "B: all teachers vs non-teachers with <d>",
                            "cc_<d>": "C: teachers with <d> vs non-teachers with <d>"},
            "units": "gap/mean_t/mean_c in latent SD (afqt_c difference / sqrt(R_comp_97)); gap_obs in afqt_c units; shares and rates in [0, 1]",
            "standard_errors": {"se, ci95_lo, ci95_hi, n_draws": "household-cluster bootstrap (NLSY79 analogue)",
                                "se_psu, ci95_lo_psu, ci95_hi_psu, n_draws_psu": "design-based Rao-Wu PSU bootstrap (use for tests against the model)"},
            "specs": {sp.name: sp.describe() for sp in fr.specs},
            "tercile_variable": "equal-weight composite of the four age-normed CAT-ASVAB estimates (nlsy_ability.pop_afqt)",
            "armed_forces": "round dropped if a military roster job (YEMP_MILFLAG = 1) is held at the interview (YEMP_CURFLAG = 1) and no "
                            "civilian job is held (the ESR_COL = 4 analogue; agrees with CV_DOI_EMPLOYED = 2 in 2015+ rounds)",
            "item5_keys": {"edu|<d>|...": "tercile cuts over the whole sample with a valid score (afqt_full); members observed at 25-34 with the "
                                          "schooling variable known",
                           "edu9|<d>|...": "the same, highest grade >= 9 only; gap keys edu9|<d>|{pooled,men,women}|gap = NLSY79 all_hgc9|<group>|gap",
                           "edu_cutm|<d>|share|q": "terciles cut over the analysis sample itself, no grade rule",
                           "edu9_cutm|<d>|[men|women|]share|q": "grade >= 9 and terciles cut over that sample = NLSY79 terc|all_hgc9|<group>|ba_share|q",
                           "edu9_cutm|<d>|cut|q": "tercile cut points (afqt_c units)"},
            "trans_variants": {**{k: {kk: list(v) if isinstance(v, tuple) else v for kk, v in kw.items()} for k, kw in TRANS_VARIANTS.items()},
                               **EXTRA_TRANS_VARIANTS}}
    payload = {"sealed": True, "meta": meta, "definitions": defs, "counts": cnt, "statistics": S, "headline": build_headline(S)}
    payload = nt.clean(payload)
    path = outdir / "sealed" / "nlsy97_holdout.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1, default=float, allow_nan=False)
    os.replace(tmp, path)
    write_methods(outdir / "nlsy97_holdout_methods.md", public_counts(cnt), meta, fr.specs, check_lines)
    return path, len(S)


def verify_seal(path: Path) -> int:
    """Structure of a sealed file, without printing any value: keys, number of statistics, how many have no value
    (empty cells), B, seed, failed draws and the smallest number of finite draws behind a statistic (household
    and PSU bootstraps)."""
    d = json.loads(path.read_text())
    S = d["statistics"]
    n_null = sum(1 for v in S.values() if v["value"] is None)
    n_nose = sum(1 for v in S.values() if v["value"] is not None and v["se"] is None)
    n_nose_p = sum(1 for v in S.values() if v["value"] is not None and v.get("se_psu") is None)
    dmin = min(v["n_draws"] for v in S.values() if v["value"] is not None)
    dmin_p = min(v.get("n_draws_psu", 0) for v in S.values() if v["value"] is not None)
    print(f"sealed file {path.name} ({path.stat().st_size} bytes): top-level keys {list(d)}; {len(S)} statistics, {n_null} without a value, "
          f"{n_nose} with a value but no household SE, {n_nose_p} with a value but no PSU SE; B = {d['meta']['B']}, seed {d['meta']['seed']}, "
          f"failed draws {d['meta']['failed_draws']} (household) / {d['meta'].get('failed_draws_psu')} (PSU); "
          f"fewest finite draws behind a statistic {dmin} (household) / {dmin_p} (PSU); headline blocks {list(d['headline'])}; "
          f"counts keys {len(d['counts'])}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seal", action="store_true", help="compute and write the sealed json (prints one line)")
    p.add_argument("--check", action="store_true", help="sample sizes and checks that involve no sealed value (default)")
    p.add_argument("--selftest", action="store_true", help="synthetic unit tests")
    p.add_argument("--verify-seal", action="store_true", help="print the structure (no values) of the sealed file in --outdir")
    p.add_argument("--bootstrap", type=int, default=B_DEFAULT)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--quick", action="store_true", help="B = 50")
    p.add_argument("--force", action="store_true", help="--seal: overwrite an existing sealed file")
    p.add_argument("--no-nlsy79", action="store_true", help="--check: skip the NLSY79 equivalence checks (--seal always runs them)")
    p.add_argument("--indir", default=str(HERE))
    p.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    args = p.parse_args(argv)
    if args.selftest:
        return selftest()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates" / "settings"
    B = 50 if args.quick else args.bootstrap
    if args.verify_seal:
        return verify_seal(outdir / "sealed" / "nlsy97_holdout.json")
    if args.seal:
        target = outdir / "sealed" / "nlsy97_holdout.json"
        if target.exists() and not args.force:
            print(f"refusing to overwrite the sealed file {target} (use --force)", file=sys.stderr)
            return 2
        logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")
        warnings.simplefilter("ignore")
        c = na.Cohort("nlsy97", indir)
        fr = Frame97(c)
        good, lines = run_checks(c, fr)
        if not good:
            print("checks failed (see `--check`); nothing sealed", file=sys.stderr)
            return 1
        path, k = run_seal(fr, B, args.seed, outdir, args.quick, lines)
        print(f"sealed: {path} ({k} moments)")
        return 0
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    t0 = time.time()
    c = na.Cohort("nlsy97", indir)
    fr = Frame97(c)
    print_counts(public_counts(sample_counts(fr)))
    print("Checks that involve no sealed value:")
    good, lines = run_checks(c, fr, with79=not args.no_nlsy79)
    print("\n".join(lines))
    print(f"all checks: {ok_str(good)}  ({time.time() - t0:.0f}s)")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
