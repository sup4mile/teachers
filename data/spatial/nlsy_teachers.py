#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy_teachers.py -- NLSY79: how the latent ability (AFQT) of K-12 teachers compares with that of
non-teachers at ages 25-34.

The spatial teachers model lets people sort into K-12 teaching on log persistent ability z (and
on a taste draw), and currently has teachers NEGATIVELY selected on z.  The data counterpart is
the gap in mean LATENT AFQT between K-12 teachers and non-teachers, in SDs of latent ability
(the model's SDs of z).  This script computes it, with the gap against non-teacher college
graduates (teachers are almost all college graduates) and the teaching share at 25-34 (men,
women) as a sanity check against the ACS.

Ability.  `afqt_c` of `nlsy_ability.py` (`pop_afqt`): the equal-weight mean of the four age-normed
ASVAB IRT z-scores (AR, WK, PC, MK), each standardized and the mean re-standardized with 1979
`SAMPWEIGHT` in the block-A population sample (SAMPLE_ID 1-8, 10, 11, 13, 14 with four valid
scores, N = 9033).  A difference of group means of `afqt_c` is divided by sqrt(R_comp_pop)
(R_comp_pop = 0.916, one-factor ULS reliability of the composite, re-estimated in every bootstrap
draw) to give the difference in SDs of the latent factor, exactly as the wage slope per latent SD
of `nlsy_ability.py` (block C) is the observed slope over sqrt(R_comp_pop).

Occupation.  The occupation of the CPS job (current job, or most recent job when not working at
the interview), a Census 3-digit code coded by the NLS: `CPSOCC80` (1980 Census codes) in every
round with ages 25-34 (survey years 1982-1998; the cohort is 14-22 in 1979, so no later round
has anyone aged 25-34); `CPSOCC70` (1970 codes, 1979-1993) and `OCCALL-EMP.01` (job #1, 1970
codes through 2000) give a 1970-vintage classification used as a sensitivity.  The extract
`raw/nlsy79_ability.parquet` already holds these (block `work` of `fetch_data.NLSY_BLOCKS`).
Each variable's vintage is read off its codebook title and, unless `--no-codebook`, checked
against the code bands of the codebook (`raw/_downloads/nlsy79_all_1979-2022.zip`, member
`.cdb`, streamed, nothing extracted).  K-12 teacher codes (labels checked against NLSY79
Codebook Supplement Attachment 3, see `ATT3_FILES`):

    1970   142 elementary, 143 pre-K and kindergarten, 144 secondary, 145 teachers n.e.c.
           (102-140 college and university teachers, 141 adult education: not K-12)
    1980   155 pre-K and kindergarten, 156 elementary, 157 secondary, 158 special education,
           159 teachers n.e.c.  (113-154 postsecondary teachers, 387 teachers' aides: not K-12)
    2000   230, 231, 232, 233, 234 (same order); 2002 4-digit 2300-2340 (rounds 2004+, not used
           here because nobody is 25-34 then)

"All K-12" (the main definition) is the project's ACS definition (Census 2300-2340: pre-K and
kindergarten, elementary and middle, secondary, special education, other teachers and
instructors).  "Core" drops pre-K/kindergarten and n.e.c.; pre-K/K, elementary, secondary,
special education and n.e.c. are also reported separately.

Samples and rule.  Person-rounds at AGEATINT 25-34 in the block-A sample, not in the armed
forces (ESR_COL 4), with a valid occupation code (1980 codes 3-889, 1970 codes 1-984; 000, 990,
995, 999 are not occupations).  Sample `emp`: rounds with ESR_COL 1 (employed at the interview);
`any`: every round with a valid occupation (including the most recent job of the unemployed and
of persons out of the labor force); `pop` (sensitivity): every respondent observed at 25-34, a
person with no occupation counting as a non-teacher.  A person is a K-12 teacher if more than
half of his or her counted rounds are K-12 teaching ("majority"; ties are non-teachers; "half"
counts ties as teachers; "ever" = at least one round).  Comparison groups: (A) all non-teachers of
the same sample; (B) non-teachers with 16+ years of schooling (highest HGCREV over the person's
25-34 rounds), against ALL teachers; (C) the within-college comparison, teachers with 16+ years
against non-teachers with 16+ years.  Weight = 1979 SAMPWEIGHT (as in `nlsy_ability.py` block C),
one observation per person.

Inference.  Household-cluster bootstrap (HHID), B = 500, seed 20260929, the draw rule of
`nlsy_extensions.t3_estimate` (cluster multiplicities from `rng.integers`); every draw
re-standardizes AFQT and re-estimates R_comp_pop with the draw weights.  SE = bootstrap SD;
CI = 95% percentile interval.

    cd data/spatial && python3 nlsy_teachers.py                 # B = 500 (about 1 min on a compute node)
    python3 nlsy_teachers.py --quick                            # B = 50
    python3 nlsy_teachers.py --no-codebook                      # skip the 165 MB codebook scan
    python3 nlsy_teachers.py --md-only                          # rewrite the .md from the .json
    python3 nlsy_teachers.py --selftest                         # unit checks, no data needed

Inputs: raw/nlsy79_ability.parquet and metadata/nlsy79_ability_vars.csv (via `nlsy_ability.Cohort`),
raw/_downloads/nlsy79_all_1979-2022.zip (codebook), raw/_downloads/nlsy79_attachment3/*.txt
(`pdftotext -layout` of the three Attachment 3 PDFs, nlsinfo.org), estimates/nlsy_ability.json
(s_z), estimates/acs_occupations.json and raw/acs_pums/acs1_20{09..13}/ (ACS benchmark shares).
Writes estimates/settings/nlsy_teachers.json and .md (no microdata).

Written with help from Claude Code.
"""
from __future__ import annotations

import os

# the linear algebra here is tiny; on a shared 72-core node BLAS/OpenMP thread pools sized to the machine
# can fail at start-up ("std::system_error: Resource temporarily unavailable"), so cap them before numpy loads
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import argparse
import datetime as dt
import glob
import io
import json
import logging
import re
import shutil
import subprocess
import sys
import time
import warnings
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nlsy_ability as na  # noqa: E402
import nlsy_extensions as ne  # noqa: E402  (table helper)

LOG = logging.getLogger("nlsy_teachers")

SEED = 20260929
B_DEFAULT = 500
COLLEGE_HGC = 16                       # "college graduate": 16+ years of schooling
LAST_CPS70_YEAR = 1993                 # CPSOCC70 covers 1979-1993; OCCALL-EMP.01 (1970 codes) after
GROUPS = ("pooled", "men", "women")
COMPS = ("nt", "col", "cc")            # A all non-teachers; B non-teacher college graduates; C college vs college

# =============================================================================
# Occupation vintages and K-12 teacher codes
# =============================================================================

VINTAGE_LABEL = {"c70": "1970 Census, 3-digit", "c80": "1980 Census, 3-digit",
                 "c00": "2000 Census, 3-digit", "c02": "2000/2002 Census, 4-digit"}
# substantive codes; outside these: 000 none, 990 same as present job, 995 did not work (1970),
# 999 not reported (1980).  590 (current armed forces) stays in range but ESR_COL 4 rounds are dropped
VALID_RANGE = {"c70": (1, 984), "c80": (3, 889)}
TYPES = ("prek", "elem", "sec", "sped", "nec")
TYPE_LABEL = {"prek": "pre-K and kindergarten", "elem": "elementary", "sec": "secondary",
              "sped": "special education", "nec": "teachers n.e.c. / other"}
TEACHER_CODES: dict[str, dict[str, tuple[int, ...]]] = {
    "c70": {"prek": (143,), "elem": (142,), "sec": (144,), "sped": (), "nec": (145,)},
    "c80": {"prek": (155,), "elem": (156,), "sec": (157,), "sped": (158,), "nec": (159,)},
    "c00": {"prek": (230,), "elem": (231,), "sec": (232,), "sped": (233,), "nec": (234,)},
    "c02": {"prek": (2300,), "elem": (2310,), "sec": (2320,), "sped": (2330,), "nec": (2340,)},
}
# neighbours that are NOT K-12 teachers: (code, why not)
NEIGHBOURS: dict[str, tuple[tuple[int, str], ...]] = {
    "c70": ((102, "college and university teachers (102-140)"), (140, "college and university teachers, "
            "subject not specified"), (141, "adult education teachers")),
    "c80": ((113, "postsecondary teachers (113-154)"), (154, "postsecondary teachers, subject not "
            "specified")),
    "c00": ((220, "postsecondary teachers"), (254, "teacher assistants")),
    "c02": ((2200, "postsecondary teachers"), (2540, "teacher assistants")),
}
DEFS = {"all": TYPES, "core": ("elem", "sec", "sped"), "prek": ("prek",), "elem": ("elem",),
        "sec": ("sec",), "sped": ("sped",), "nec": ("nec",)}
DEF_LABEL = {"all": "all K-12 codes", "core": "core (elementary, secondary, special education)",
             "prek": "pre-K/kindergarten only", "elem": "elementary only", "sec": "secondary only",
             "sped": "special education only", "nec": "teachers n.e.c. only"}
RULE_LABEL = {"maj": "majority of rounds", "half": "at least half of rounds",
              "ever": "ever at 25-34"}
SAMPLE_LABEL = {"emp": "employed at the interview", "any": "any 25-34 occupation",
                "pop": "all respondents observed at 25-34"}

ATT3_BASE = "https://www.nlsinfo.org/sites/default/files/attachments/121217/"
ATT3_FILES = {"c70": "Attachment-3-1970_Industry_and_Occupational_Codes",
              "c80": "Attachment-3-1980-Census-Codes", "c00": "att300", "c02": "att300"}


def vintage_from_title(title: str) -> str | None:
    """Census vintage from an NLSY79 variable title, e.g. '(70 CENSUS 3 DIGIT)',
    '(CENSUS 3 DIGIT, 70 CODES)', '(CENSUS 4 DIGIT, 00 CODES)'."""
    t = title.upper()
    if "80 CENSUS" in t:
        return "c80"
    if "70 CENSUS" in t or "70 CODES" in t:
        return "c70"
    if "4 DIGIT, 00 CODES" in t:
        return "c02"
    if "3 DIGIT, 00 CODES" in t:
        return "c00"
    return None


def vintage_from_bands(text: str) -> str | None:
    """Census vintage from the first band labels of a codebook entry."""
    t = re.sub(r"\s+", "", text.upper())
    if "PROFESSIONAL,TECHNICALANDKINDRED" in t:
        return "c70"
    if "MANAGERIALANDPROFESSIONALSPECIALITIES" in t:
        return "c80"
    if "EXECUTIVE,ADMINISTRATIVEANDMANAGERIAL" in t:
        return "c02"
    if "BUSINESSANDFINANCIALOPERATIONS" in t:
        return "c00"
    return None


def valid_code(x: np.ndarray, vint: str) -> np.ndarray:
    lo, hi = VALID_RANGE[vint]
    return np.isfinite(x) & (x >= lo) & (x <= hi)


def classify(nr: np.ndarray, k: np.ndarray, rule: str) -> np.ndarray:
    """K-12 teacher indicator from the number of counted rounds `nr` and of K-12 rounds `k`."""
    has = nr > 0
    if rule == "maj":
        t = 2 * k > nr
    elif rule == "half":
        t = 2 * k >= nr
    elif rule == "ever":
        t = k >= 1
    else:
        raise ValueError(rule)
    return has & t


# =============================================================================
# Codebook, Attachment 3 and ACS benchmark
# =============================================================================

def read_codebook(zip_path: Path, qnames: set[str]) -> dict[tuple[str, str], dict[str, Any]]:
    """Stream the `.cdb` member of the NLS zip and keep the entries of `qnames`:
    {(qname, survey year): refnum, title, n_valid, bands (first labels), codes}."""
    hdr = re.compile(r"^([A-Z]\d{5}\.\d\d)\s+\[([^\]]+)\]\s+Survey Year:\s*(\S+)")
    code_line = re.compile(r"^\s*(\d+)\s+(-?\d+)\s+([A-Za-z].*?)\s*$")
    band_line = re.compile(r"^\s*(\d+)\s+(\d+)(?: TO (\d+))?:\s*(.*?)\s*$")
    out: dict[tuple[str, str], dict[str, Any]] = {}

    def keep(block: list[str]) -> None:
        m = hdr.match(block[0]) if block else None
        if not m or m.group(2) not in qnames:
            return
        ref, q, yr = m.groups()
        title, bands, codes, n_valid = [], [], {}, None
        in_title = False
        after_rule = False
        for ln in block[1:]:
            if ln.strip().startswith("PRIMARY VARIABLE"):
                in_title = True
                continue
            if in_title:
                if ln.strip():
                    title.append(ln.strip())
                elif title:
                    in_title = False
                continue
            if ln.strip().startswith("-------"):
                after_rule = True
                continue
            if after_rule and n_valid is None and ln.strip().isdigit():
                n_valid = int(ln.strip())
                continue
            if after_rule:
                continue
            b = band_line.match(ln)
            if b:
                bands.append(b.group(4))
                if b.group(3) is None:                       # '1: 1  EMPLOYED' (later rounds)
                    codes[int(b.group(2))] = re.sub(r"^\d+\s+", "", b.group(4)).strip()
                continue
            c = code_line.match(ln)
            if c:
                codes[int(c.group(2))] = c.group(3)
        out[(q, yr)] = {"refnum": ref, "title": " ".join(title), "n_valid": n_valid,
                        "bands": bands[:3], "codes": codes}

    with zipfile.ZipFile(zip_path) as zf:
        member = next(m for m in zf.namelist() if m.lower().endswith(".cdb"))
        with zf.open(member) as fh:
            cur: list[str] = []
            for raw in io.TextIOWrapper(fh, encoding="latin-1"):
                line = raw.rstrip("\n")
                if line.startswith("-" * 80):
                    keep(cur)
                    cur = []
                else:
                    cur.append(line)
            keep(cur)
    return out


def _find_label(lines: list[str], code: int) -> str | None:
    """Label of `code` in a pdftotext listing; the label must name teachers."""
    pat = re.compile(rf"^\s*(?:[A-Z]\s*\(\s*)?{code}\s*\)?\s+(?:\d{{1,3}}\s{{2,}})?(.*?)\s*$")
    for ln in lines:
        m = pat.match(ln)
        if m:
            lab = re.sub(r"\s{2,}\d{2}-\d{4}\s*$", "", m.group(1)).strip()
            if "teacher" in lab.lower():
                return lab
    return None


def attachment3_labels(att_dir: Path) -> dict[str, dict[int, str | None]] | None:
    """Official labels of the K-12 teacher codes (and of the neighbouring non-K-12 codes) in
    NLSY79 Attachment 3, from `pdftotext -layout` listings in `att_dir`."""
    res: dict[str, dict[int, str | None]] = {}
    for v, stem in ATT3_FILES.items():
        txt, pdf = att_dir / f"{stem}.txt", att_dir / f"{stem}.pdf"
        if not txt.exists() and pdf.exists() and shutil.which("pdftotext"):
            subprocess.run(["pdftotext", "-layout", str(pdf), str(txt)], check=True)
        if not txt.exists():
            LOG.warning("  Attachment 3 listing for %s not found (%s): labels unverified", v, txt)
            continue
        lines = txt.read_text(errors="ignore").splitlines()
        codes = {c for ts in TEACHER_CODES[v].values() for c in ts} | {c for c, _ in NEIGHBOURS[v]}
        res[v] = {c: _find_label(lines, c) for c in sorted(codes)}
    return res or None


def acs_benchmark(root: Path) -> dict[str, Any] | None:
    """ACS 2009-13 (pooled 1-year PUMS in raw/acs_pums), ages 25-34, SCHL >= 12 (grade 9+), not in
    the armed forces: K-12 teaching share (OCCP 2300-2340, employed) and its split by
    education, PWGTP-weighted persons, simple rule (no hours-based home-production weights)."""
    files = [f for y in range(2009, 2014)
             for f in sorted(glob.glob(str(root / "raw" / "acs_pums" / f"acs1_{y}" / "*.parquet")))]
    if not files:
        LOG.warning("  raw/acs_pums not found: ACS benchmark skipped")
        return None
    cols = ["AGEP", "SEX", "SCHL", "ESR", "OCCP", "PWGTP"]
    A = pd.concat([pd.read_parquet(f, columns=cols) for f in files], ignore_index=True)
    for c in cols:
        A[c] = pd.to_numeric(A[c], errors="coerce")
    A = A[(A.AGEP >= 25) & (A.AGEP <= 34) & (A.SCHL >= 12) & ~A.ESR.isin([4, 5])]
    A = A.assign(emp=A.ESR.isin([1, 2]), ba=A.SCHL >= 21)
    A["k12"] = A.emp & A.OCCP.between(2300, 2340)
    out: dict[str, Any] = {"n_persons": int(len(A)), "files": len(files),
                           "rule": "ages 25-34, SCHL >= 12, ESR not 4/5; K-12 = OCCP 2300-2340 and "
                                   "ESR 1/2; PWGTP weights; BA+ = SCHL >= 21"}
    for sx, lab in ((1, "men"), (2, "women")):
        d = A[A.SEX == sx]
        w = d.PWGTP.to_numpy(float)
        k, ba, emp = (d[c].to_numpy(bool) for c in ("k12", "ba", "emp"))
        out[lab] = {"n": int(len(d)), "k12_share": float((w * k).sum() / w.sum()),
                    "emp_rate": float((w * emp).sum() / w.sum()),
                    "ba_share": float((w * ba).sum() / w.sum()),
                    "k12_given_ba": float((w * k * ba).sum() / (w * ba).sum()),
                    "k12_given_noba": float((w * k * ~ba).sum() / (w * ~ba).sum())}
    return out


# =============================================================================
# Data: person-rounds, person-level classification
# =============================================================================

@dataclass
class Spec:
    name: str
    vint: str            # "c80" | "c70"
    samp: str            # "emp" | "any" | "pop"
    defn: str            # key of DEFS
    rule: str            # "maj" | "half" | "ever"
    hgc9: bool = False   # restrict to persons with highest grade >= 9 at 25-34
    wmode: str = "w79"   # "w79" | "xs" | "unit" | "round"
    member: np.ndarray | None = None
    teacher: np.ndarray | None = None
    nonteacher: np.ndarray | None = None

    def describe(self) -> dict[str, Any]:
        return {"vintage": self.vint, "sample": self.samp, "definition": self.defn, "rule": self.rule,
                "hgc_ge_9": self.hgc9, "weights": self.wmode}


def make_specs() -> list[Spec]:
    specs: list[Spec] = []

    def add(vint, samp, defn, rule, hgc9=False, wmode="w79"):
        name = "-".join([vint, samp, defn, rule] + (["hgc9"] if hgc9 else [])
                        + ([] if wmode == "w79" else [wmode]))
        specs.append(Spec(name, vint, samp, defn, rule, hgc9, wmode))

    for samp in ("emp", "any"):
        for defn in ("all", "core"):
            for rule in ("maj", "half", "ever"):
                add("c80", samp, defn, rule)
        for defn in ("prek", "elem", "sec", "sped", "nec"):
            add("c80", samp, defn, "maj")
        add("c70", samp, "all", "maj")
        for defn in ("core", "prek", "elem", "sec", "nec"):
            add("c70", samp, defn, "maj")
        add("c80", samp, "all", "maj", hgc9=True)
        for wm in ("xs", "unit", "round"):
            add("c80", samp, "all", "maj", wmode=wm)
    for rule in ("maj", "half", "ever"):
        add("c80", "pop", "all", rule)
    return specs


class Frame:
    """Everything that does not change across bootstrap draws."""

    def __init__(self, c: na.Cohort) -> None:
        self.c = c
        D = self.D = na.build_n79(c)
        py = self.py = D.py
        n = self.n = D.n
        npr = len(py.p)
        self.female = D.female.astype(bool)
        self.masks = {"pooled": np.ones(n, bool), "men": ~self.female, "women": self.female}
        self.xs = np.isin(D.sample, na.XS_SAMPLES)
        self.pop_py = D.asamp[py.p]
        # ---- occupation codes by person-round --------------------------------------------
        occ = {k: np.full(npr, np.nan) for k in ("cps70", "all70", "cps80")}
        wr = np.full(npr, np.nan)
        for Y in sorted(set(py.Y.tolist())):
            s = py.Y == Y
            wr[s] = c.get("SAMPWEIGHT", Y)[py.p[s]] / 100.0
            for key, q in (("cps70", "CPSOCC70"), ("all70", "OCCALL-EMP.01"), ("cps80", "CPSOCC80")):
                v = c.try_get(q, Y)
                if v is not None:
                    occ[key][s] = v[py.p[s]]
        self.round_weight = wr
        self.occ80 = occ["cps80"]
        self.occ70 = np.where(py.Y <= LAST_CPS70_YEAR, occ["cps70"], occ["all70"])
        self.valid = {"c80": valid_code(self.occ80, "c80"), "c70": valid_code(self.occ70, "c70")}
        self.occ = {"c80": self.occ80, "c70": self.occ70}
        self.typ = {v: {t: self.valid[v] & np.isin(self.occ[v], TEACHER_CODES[v][t]) for t in TYPES}
                    for v in ("c80", "c70")}
        self.k12 = {v: np.any([self.typ[v][t] for t in TYPES], axis=0) for v in ("c80", "c70")}
        # ---- rounds that count, person-level counts -------------------------------------
        self.base = self.pop_py & (py.esr != 4)
        self.R: dict[tuple[str, str], dict[str, Any]] = {}
        for v in ("c80", "c70"):
            for s in ("any", "emp"):
                m = self.base & self.valid[v] & ((py.esr == 1) if s == "emp" else True)
                nr = np.bincount(py.p[m], minlength=n)
                nt = {t: np.bincount(py.p[m & self.typ[v][t]], minlength=n) for t in TYPES}
                okw = m & np.isfinite(wr) & (wr > 0)
                sw = np.bincount(py.p[okw], weights=wr[okw], minlength=n)
                cw = np.bincount(py.p[okw], minlength=n)
                self.R[(v, s)] = {"mask": m, "nr": nr, "nt": nt,
                                  "wround": np.where(cw > 0, sw / np.maximum(cw, 1), 0.0)}
        self.observed = np.bincount(py.p[self.base], minlength=n) > 0
        hg = np.full(n, np.nan)
        np.fmax.at(hg, py.p[self.pop_py], py.hgc[self.pop_py])
        self.hg_max = hg
        self.college = np.isfinite(hg) & (hg >= COLLEGE_HGC)
        self.hgc9 = np.isfinite(hg) & (hg >= na.MIN_HGC)
        self.specs = make_specs()
        for sp in self.specs:
            self._build(sp)
        LOG.info("  persons %d (block A %d), households %d; person-rounds at 25-34 %d (block A %d)",
                 n, int(D.asamp.sum()), D.nhh, npr, int(self.pop_py.sum()))

    def _build(self, sp: Spec) -> None:
        R = self.R[(sp.vint, "any" if sp.samp == "pop" else sp.samp)]
        nr = R["nr"]
        k = sum(R["nt"][t] for t in DEFS[sp.defn])
        k_all = sum(R["nt"][t] for t in DEFS["all"])
        member = (self.observed if sp.samp == "pop" else nr > 0) & self.D.asamp
        if sp.hgc9:
            member = member & self.hgc9
        if sp.wmode == "xs":
            member = member & self.xs
        teacher_all = classify(nr, k_all, sp.rule) & member
        sp.member = member
        sp.teacher = classify(nr, k, sp.rule) & member
        sp.nonteacher = member & ~teacher_all          # persons in other K-12 types are in neither group


# =============================================================================
# Estimates (re-run on every bootstrap draw)
# =============================================================================

def wm(a: np.ndarray, w: np.ndarray, m: np.ndarray) -> float:
    s = w[m].sum()
    return float((a[m] * w[m]).sum() / s) if s > 0 else np.nan


def ratio(num: float, den: float) -> float:
    return float(num / den) if den > 0 else np.nan


def evaluate(fr: Frame, mult_hh: np.ndarray | None, info: dict[str, Any] | None = None
             ) -> dict[str, float]:
    """All statistics for household multiplicities `mult_hh` (None = the sample itself).
    `info`, if given, receives unweighted counts.  Keys: `spec|group|comparison|stat` with stat in
    mean_t, mean_c (latent SDs from the block-A mean), gap (latent SDs), gap_obs (afqt_c units)."""
    D = fr.D
    mult_p = np.ones(D.n) if mult_hh is None else mult_hh[D.hh_idx].astype(float)
    w79 = D.w * mult_p
    out: dict[str, float] = {}
    afqt, Rc, *_ = na.pop_afqt(D.z, D.asamp, w79)
    afx, Rx, *_ = na.pop_afqt(D.z, D.asamp & fr.xs, w79)
    afu, Ru, *_ = na.pop_afqt(D.z, D.asamp, mult_p)
    out.update({"R_comp_pop": Rc, "R_comp_pop_xs": Rx, "R_comp_pop_unit": Ru})
    A = {"w79": (afqt, np.sqrt(Rc)), "round": (afqt, np.sqrt(Rc)), "xs": (afx, np.sqrt(Rx)),
         "unit": (afu, np.sqrt(Ru))}
    cnt: dict[str, Any] = {}
    for sp in fr.specs:
        a, sR = A[sp.wmode]
        if sp.wmode == "unit":
            wv = mult_p
        elif sp.wmode == "round":
            wv = fr.R[(sp.vint, sp.samp)]["wround"] * mult_p
        else:
            wv = w79
        ok = sp.member & np.isfinite(a) & (wv > 0)
        t, nt = sp.teacher & ok, sp.nonteacher & ok
        for g, gm in fr.masks.items():
            k = f"{sp.name}|{g}"
            m1, m0 = t & gm, nt & gm
            m1c, m0c = m1 & fr.college, m0 & fr.college
            mt, mtc, mc, mcc = (wm(a, wv, m) for m in (m1, m1c, m0, m0c))
            out[f"{k}|share"] = ratio(wv[m1].sum(), wv[m1 | m0].sum())
            out[f"{k}|col_share_t"] = ratio(wv[m1c].sum(), wv[m1].sum())
            out[f"{k}|col_share_nt"] = ratio(wv[m0c].sum(), wv[m0].sum())
            cnt.setdefault(sp.name, {})[g] = {}
            # A: all teachers vs all non-teachers; B: all teachers vs non-teacher college graduates;
            # C: college-graduate teachers vs non-teacher college graduates
            for cl, (x1, y1, x0, y0) in {"nt": (m1, mt, m0, mc), "col": (m1, mt, m0c, mcc),
                                         "cc": (m1c, mtc, m0c, mcc)}.items():
                out[f"{k}|{cl}|mean_t"] = y1 / sR
                out[f"{k}|{cl}|mean_c"] = y0 / sR
                out[f"{k}|{cl}|gap"] = (y1 - y0) / sR
                out[f"{k}|{cl}|gap_obs"] = y1 - y0
                cnt[sp.name][g][cl] = {"n_t": int(x1.sum()), "n_c": int(x0.sum())}
    out.update(share_stats(fr, w79))
    if info is not None:
        info["counts"] = cnt
    return out


def share_stats(fr: Frame, w79: np.ndarray) -> dict[str, float]:
    """Teaching shares at 25-34.  Round shares: person-rounds of the block-A sample with the
    ACS-like denominator (highest grade >= 9 at the round, not armed forces; `acs`) or without the
    grade restriction (`all`), weights = the person's 1979 weight.  Person shares: persons
    classified as K-12 teachers (employed sample, 1980 codes, all K-12 codes) over persons observed
    at 25-34 (highest grade >= 9)."""
    py = fr.py
    wpy = w79[py.p]
    fpy = fr.female[py.p]
    emp = py.esr == 1
    k12 = fr.k12["c80"]
    core = np.any([fr.typ["c80"][t] for t in DEFS["core"]], axis=0)
    colr = np.isfinite(py.hgc) & (py.hgc >= COLLEGE_HGC)
    gmr = {"pooled": np.ones(len(fpy), bool), "men": ~fpy, "women": fpy}
    out: dict[str, float] = {}
    dens = {"acs": fr.base & np.isfinite(py.hgc) & (py.hgc >= na.MIN_HGC), "all": fr.base}
    for dl, den0 in dens.items():
        for g, gm in gmr.items():
            d = den0 & gm & (wpy > 0)
            tot = wpy[d].sum()
            num = {"emp_k12": emp & k12, "any_k12": k12, "emp_core": emp & core,
                   "emp_k12_c70": emp & fr.k12["c70"]}
            num.update({f"emp_{t}": emp & fr.typ["c80"][t] for t in TYPES})
            for nl, nm in num.items():
                out[f"share|{dl}|{nl}|{g}"] = ratio(wpy[d & nm].sum(), tot)
            out[f"share|{dl}|emp_rate|{g}"] = ratio(wpy[d & emp].sum(), tot)
            out[f"share|{dl}|col|{g}"] = ratio(wpy[d & colr].sum(), tot)
            out[f"share|{dl}|emp_k12_given_emp|{g}"] = ratio(wpy[d & emp & k12].sum(), wpy[d & emp].sum())
            out[f"share|{dl}|k12_given_col|{g}"] = ratio(wpy[d & colr & emp & k12].sum(), wpy[d & colr].sum())
            out[f"share|{dl}|k12_given_noncol|{g}"] = ratio(wpy[d & ~colr & emp & k12].sum(), wpy[d & ~colr].sum())
    obs_acs = np.bincount(py.p[dens["acs"]], minlength=fr.n) > 0
    Rm = fr.R[("c80", "emp")]
    k_all = sum(Rm["nt"][t] for t in TYPES)
    for rule in ("maj", "half", "ever"):
        tch = classify(Rm["nr"], k_all, rule) & obs_acs & fr.D.asamp
        for g, gm in fr.masks.items():
            den = obs_acs & fr.D.asamp & gm
            out[f"pshare|{rule}|{g}"] = ratio(w79[tch & gm].sum(), w79[den].sum())
    return out


def run_bootstrap(fr: Frame, B: int, seed: int) -> dict[str, Any]:
    info: dict[str, Any] = {}
    point = evaluate(fr, None, info)
    keys = list(point)
    rng = np.random.default_rng(seed)
    draws = np.full((B, len(keys)), np.nan)
    nhh = fr.D.nhh
    t0 = time.time()
    for b in range(B):
        cm = np.bincount(rng.integers(0, nhh, nhh), minlength=nhh).astype(float)
        try:
            d = evaluate(fr, cm)
            draws[b] = [d.get(k, np.nan) for k in keys]
        except Exception as exc:                                  # a failed draw is a NaN row
            LOG.warning("  draw %d failed: %s", b, exc)
        if (b + 1) % 50 == 0 or b + 1 == B:
            LOG.info("  bootstrap %d/%d (%.1fs, %.2fs/draw)", b + 1, B, time.time() - t0,
                     (time.time() - t0) / (b + 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        se = np.nanstd(draws, axis=0, ddof=1) if B > 1 else np.full(len(keys), np.nan)
        lo, hi = np.nanpercentile(draws, [2.5, 97.5], axis=0)
    n_ok = np.sum(np.isfinite(draws), axis=0)
    stats = {k: {"value": float(point[k]), "se": float(se[i]), "ci95_lo": float(lo[i]),
                 "ci95_hi": float(hi[i]), "n_draws": int(n_ok[i])} for i, k in enumerate(keys)}
    return {"stats": stats, "counts": info["counts"], "B": B, "seed": seed,
            "failed_draws": int(np.sum(~np.isfinite(draws).any(axis=1)))}


# =============================================================================
# Descriptive tables that do not need the bootstrap
# =============================================================================

def round_map(fr: Frame, meta: pd.DataFrame, cdb: dict | None) -> list[dict[str, Any]]:
    """One row per (occupation variable, survey year): reference number, vintage from the title
    (and from the codebook bands), valid counts, and, for the person-rounds at 25-34 (block A, not
    armed forces), the valid and K-12 counts of that variable itself."""
    py, c = fr.py, fr.c
    rows = []
    for r in meta[meta.qname.isin(["CPSOCC70", "CPSOCC80", "OCCALL-EMP.01"])].itertuples():
        Y = int(r.year)
        x = c.get(r.qname, Y)
        vint = vintage_from_title(r.title)
        s = (py.Y == Y) & fr.base
        row: dict[str, Any] = {"year": Y, "variable": r.qname, "rnum": r.rnum, "refnum": r.refnum,
                               "title": r.title, "vintage": vint, "n_valid_extract": int(np.isfinite(x).sum()),
                               "n_person_rounds_25_34": int(s.sum())}
        if cdb is not None and (r.qname, str(Y)) in cdb:
            e = cdb[(r.qname, str(Y))]
            row.update({"codebook_refnum": e["refnum"], "codebook_title": e["title"],
                        "n_valid_codebook": e["n_valid"], "vintage_from_bands": vintage_from_bands(" ".join(e["bands"]))})
        if s.any() and vint in VALID_RANGE:
            xs_ = x[py.p[s]]
            ok = valid_code(xs_, vint)
            codes = [cd for t in TYPES for cd in TEACHER_CODES[vint][t]]
            row["valid_25_34"] = int(ok.sum())
            row["k12_25_34"] = int((ok & np.isin(xs_, codes)).sum())
            if r.qname == "CPSOCC80":
                row["use"] = "main (1980 codes)"
            elif r.qname == "CPSOCC70" or (r.qname == "OCCALL-EMP.01" and Y > LAST_CPS70_YEAR):
                row["use"] = "1970-code sensitivity"
        rows.append(row)
    return sorted(rows, key=lambda d: (d["year"], d["variable"]))


def code_counts(fr: Frame) -> dict[str, Any]:
    """Person-rounds (unweighted, block A, ages 25-34, not armed forces) by K-12 teacher code, and
    the 1970-1980 concordance on the same person-rounds."""
    py = fr.py
    res: dict[str, Any] = {"by_code": {}}
    for v in ("c80", "c70"):
        res["by_code"][v] = {}
        for t in TYPES:
            for code in TEACHER_CODES[v][t]:
                m = fr.base & fr.valid[v] & (fr.occ[v] == code)
                res["by_code"][v][str(code)] = {"type": t, "person_rounds": int(m.sum()),
                                                "employed_person_rounds": int((m & (py.esr == 1)).sum())}
    both = fr.base & fr.valid["c80"] & fr.valid["c70"]
    rows70, cols80 = [141, 142, 143, 144, 145], [155, 156, 157, 158, 159]
    mat = {}
    for r in rows70 + ["other"]:
        mr = both & ((fr.occ["c70"] == r) if r != "other" else ~np.isin(fr.occ["c70"], rows70))
        mat[str(r)] = {str(cc): int((mr & ((fr.occ["c80"] == cc) if cc != "other"
                                           else ~np.isin(fr.occ["c80"], cols80))).sum())
                       for cc in cols80 + ["other"]}
    k70, k80 = fr.k12["c70"] & both, fr.k12["c80"] & both
    res["concordance"] = {
        "matrix_rows_1970_cols_1980": mat, "n_person_rounds_both_valid": int(both.sum()),
        "k12_both": int((k70 & k80).sum()), "k12_only_1970": int((k70 & ~k80).sum()),
        "k12_only_1980": int((~k70 & k80).sum()), "neither": int((both & ~k70 & ~k80).sum()),
        "years": sorted(set(py.Y[both].tolist()))}
    return res


# =============================================================================
# Report
# =============================================================================

def finite(x: Any) -> bool:
    return x is not None and bool(np.isfinite(x))


def clean(obj: Any) -> Any:
    """JSON-safe copy: non-finite floats (cells with no persons) become null."""
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def cell(S: dict[str, dict[str, float]], key: str, d: int = 3, ci: bool = False) -> str:
    r = S.get(key)
    if r is None or not finite(r["value"]):
        return "n/a"
    s = f"{r['value']:.{d}f} ({r['se']:.{d}f})"
    if ci:
        s += f" [{r['ci95_lo']:.{d}f}, {r['ci95_hi']:.{d}f}]"
    return s


def pct(S: dict[str, dict[str, float]], key: str, d: int = 2) -> str:
    r = S.get(key)
    if r is None or not finite(r["value"]):
        return "n/a"
    return f"{100 * r['value']:.{d}f} ({100 * r['se']:.{d}f})"


def year_ranges(ys: list[int]) -> str:
    ys = sorted(ys)
    out, i = [], 0
    while i < len(ys):
        j = i
        while j + 1 < len(ys) and ys[j + 1] - ys[j] <= 2:
            j += 1
        out.append(str(ys[i]) if i == j else f"{ys[i]}-{ys[j]}")
        i = j + 1
    return ", ".join(out)


def val(S: dict[str, dict[str, float]], key: str) -> float:
    return S[key]["value"]


def spec_label(sp: dict[str, Any]) -> str:
    parts = [DEF_LABEL[sp["definition"]], RULE_LABEL[sp["rule"]]]
    if sp["vintage"] == "c70":
        parts.append("1970 codes")
    if sp["hgc_ge_9"]:
        parts.append("highest grade >= 9")
    if sp["weights"] != "w79":
        parts.append({"xs": "cross-sectional sample (SAMPLE_ID 1-8)", "unit": "unweighted",
                      "round": "mean round weights"}[sp["weights"]])
    return "; ".join(parts)


def gap_table(S: dict, C: dict, comp: str, rows_spec: list[tuple[str, str]]) -> str:
    rows = []
    for spec, lab in rows_spec:
        for g in GROUPS:
            n = C[spec][g][comp]
            k = f"{spec}|{g}|{comp}|"
            rows.append([lab if g == "pooled" else "", g, str(n["n_t"]), str(n["n_c"]), cell(S, k + "mean_t"),
                         cell(S, k + "mean_c"), cell(S, k + "gap", 3, True)])
    return ne.tab(rows, ["sample", "group", "teachers", "comparison persons", "mean latent AFQT, teachers",
                         "mean latent AFQT, comparison", "gap (teachers - comparison)"])


def write_md(payload: dict[str, Any], path: Path) -> None:
    meta, S, C, SP = payload["meta"], payload["statistics"], payload["counts"], payload["specs"]
    cc = payload["code_counts"]["concordance"]
    ac = payload["acs_benchmark"]
    L: list[str] = []
    a = L.append
    sz = meta.get("s_z")
    E, N = "c80-emp-all-maj", "c80-any-all-maj"            # main specifications
    yrs = meta["years"]

    def g_(spec: str, g: str, comp: str, ci: bool = False) -> str:
        return cell(S, f"{spec}|{g}|{comp}|gap", 3, ci)

    def n_(spec: str, g: str, comp: str = "nt", which: str = "n_t") -> int:
        return C[spec][g][comp][which]

    closed: dict[str, float] = {}
    cfs: dict[str, float] = {}
    if ac:
        for g in ("men", "women"):
            cfs[g] = (ac[g]["ba_share"] * val(S, f"share|acs|k12_given_col|{g}")
                      + (1 - ac[g]["ba_share"]) * val(S, f"share|acs|k12_given_noncol|{g}"))
            closed[g] = (cfs[g] - val(S, f"share|acs|emp_k12|{g}")) / (ac[g]["k12_share"] - val(S, f"share|acs|emp_k12|{g}"))
    ratio_m = val(S, "share|acs|emp_k12|men") / meta["acs_published_pwgtp"]["M"]
    ratio_w = val(S, "share|acs|emp_k12|women") / meta["acs_published_pwgtp"]["F"]

    a("# NLSY79: latent ability (AFQT) of K-12 teachers and non-teachers at ages 25-34\n")
    a(f"Generated {meta['generated']} by `nlsy_teachers.py`. NLSY79 population sample (block A of "
      f"`nlsy_ability.py`, N = {meta['n_block_a']}), person-rounds at ages 25-34 (survey years {yrs[0]}-{yrs[-1]}). "
      "Ability = `afqt_c` of `nlsy_ability.py` (standardized, age-normed composite of the four ASVAB IRT z-scores, "
      "1979 weights); a gap is a difference of group means of `afqt_c` divided by sqrt(R_comp_pop), i.e. in SDs of "
      f"the latent factor (R_comp_pop = {meta['R_comp_pop']:.3f}). Weights: 1979 `SAMPWEIGHT`, one observation per "
      f"person. Household-cluster bootstrap, B = {meta['B']}, seed {meta['seed']}; SE = bootstrap SD, in "
      "parentheses; 95% percentile CI in brackets; AFQT standardization and R_comp_pop are re-estimated in every "
      "draw. Comparison A: all non-teachers; B: non-teachers with 16+ years of schooling (against all teachers); "
      "C: teachers with 16+ years against non-teachers with 16+ years.\n")

    # ---- summary -----------------------------------------------------------------------
    a("## Summary\n")
    a(f"- **Relative to all non-teachers, K-12 teachers are positively selected on latent AFQT: "
      f"{g_(E, 'pooled', 'nt')} SD** (employed; teacher = K-12 teacher in the majority of 25-34 rounds; "
      f"{n_(E, 'pooled')} teachers, {n_(E, 'pooled', 'nt', 'n_c')} non-teachers). Men {g_(E, 'men', 'nt')} "
      f"({n_(E, 'men')} teachers), women {g_(E, 'women', 'nt')} ({n_(E, 'women')} teachers). "
      f"All respondents with any 25-34 occupation: {g_(N, 'pooled', 'nt')} ({n_(N, 'pooled')} teachers). "
      f"The gap runs through schooling: {100 * val(S, E + '|pooled|col_share_t'):.0f}% of these teachers and "
      f"{100 * val(S, E + '|pooled|col_share_nt'):.0f}% of employed non-teachers have 16+ years.")
    a(f"- **Relative to non-teacher college graduates, teachers are negatively selected: {g_(E, 'pooled', 'col')} SD** "
      f"(B; employed; {n_(E, 'pooled', 'col', 'n_c')} comparison persons); men {g_(E, 'men', 'col')}, women "
      f"{g_(E, 'women', 'col')}. Any 25-34 occupation: {g_(N, 'pooled', 'col')}. Within college graduates only (C): "
      f"{g_(E, 'pooled', 'cc')} ({n_(E, 'pooled', 'cc')} teachers with 16+ years against {n_(E, 'pooled', 'cc', 'n_c')}).")
    if sz:
        a(f"- In the model's log-z units a gap of g latent SDs is g x s_z log points, s_z = {sz:.4f} "
          f"(`implied_s_z` of `nlsy_ability.py`; its sampling error is not propagated): comparison A "
          f"{val(S, E + '|pooled|nt|gap'):+.3f} SD = {val(S, E + '|pooled|nt|gap') * sz:+.3f}; B "
          f"{val(S, E + '|pooled|col|gap'):+.3f} SD = {val(S, E + '|pooled|col|gap') * sz:+.3f}; C "
          f"{val(S, E + '|pooled|cc|gap'):+.3f} SD = {val(S, E + '|pooled|cc|gap') * sz:+.3f}.")
    a(f"- **Ever a K-12 teacher at 25-34** (employed sample, {n_('c80-emp-all-ever', 'pooled')} persons, "
      f"{n_('c80-emp-all-ever', 'men')} men and {n_('c80-emp-all-ever', 'women')} women): A "
      f"{g_('c80-emp-all-ever', 'pooled', 'nt')}, B {g_('c80-emp-all-ever', 'pooled', 'col')}, C "
      f"{g_('c80-emp-all-ever', 'pooled', 'cc')}.")
    a(f"- **Teaching shares at 25-34** (employed K-12 teachers over all person-rounds with highest grade >= 9, not "
      f"armed forces, home production included): men {pct(S, 'share|acs|emp_k12|men')}%, women "
      f"{pct(S, 'share|acs|emp_k12|women')}%, against ACS 2009-13 {100 * meta['acs_published']['M']:.2f}% and "
      f"{100 * meta['acs_published']['F']:.2f}% (project headline, unweighted persons; PWGTP-weighted "
      f"{100 * meta['acs_published_pwgtp']['M']:.2f}% and {100 * meta['acs_published_pwgtp']['F']:.2f}%): "
      f"{100 * ratio_m:.0f}% and {100 * ratio_w:.0f}% of the weighted ACS shares."
      + (f" Giving the NLSY rates the ACS share of bachelor's holders closes {100 * closed['women']:.0f}% of the gap for "
         f"women and {100 * closed['men']:.0f}% for men (section 5).\n" if closed else "\n"))

    # ---- section 1: variables and vintages ----------------------------------------------
    a("## 1. Occupation variables, code vintages and K-12 teacher codes\n")
    a("The extract `raw/nlsy79_ability.parquet` already contains the CPS-job occupation variables (block `work` of "
      "`fetch_data.NLSY_BLOCKS`), so nothing was cut from the full-cohort CSV. Ages 25-34 occur only in survey "
      f"years {yrs[0]}-{yrs[-1]} (born 1957-64). Per round, with reference numbers and, for block-A person-rounds at "
      "ages 25-34 not in the armed forces, the number with a valid occupation code and with a K-12 teacher code "
      "(`80` = `CPSOCC80`, 1980 Census codes, the main classification; `70` = `CPSOCC70` through 1993 and "
      "`OCCALL-EMP.01`, job #1, 1970 codes, from 1994, the sensitivity):\n")
    by = {}
    for r in payload["round_map"]:
        by.setdefault(r["year"], {})[r["variable"]] = r
    rows = []
    for y in yrs:
        r80 = by[y]["CPSOCC80"]
        r70 = by[y].get("CPSOCC70") if y <= LAST_CPS70_YEAR else by[y].get("OCCALL-EMP.01")
        rows.append([str(y), str(r80["n_person_rounds_25_34"]), r80["refnum"], r70["refnum"] if r70 else "-",
                     str(r80.get("valid_25_34", "")), str(r80.get("k12_25_34", "")),
                     str(r70.get("valid_25_34", "")) if r70 else "", str(r70.get("k12_25_34", "")) if r70 else ""])
    a(ne.tab(rows, ["survey year", "person-rounds", "CPSOCC80", "1970-code variable", "valid (80)", "K-12 (80)",
                    "valid (70)", "K-12 (70)"]))
    oc = [r["n_valid_extract"] for r in payload["round_map"] if r["variable"] == "OCCALL-EMP.01" and 1982 <= r["year"] < 1994]
    occ_cov = f"{min(oc)}-{max(oc)} respondents a round in 1982-92, against {max(by[1994]['OCCALL-EMP.01']['n_valid_extract'], 0)} in 1994"
    grp: dict[tuple[str, str], list[int]] = {}
    for r in payload["round_map"]:
        if r["year"] not in yrs:
            grp.setdefault((r["variable"], r["vintage"]), []).append(r["year"])
    a("Other rounds (nobody aged 25-34; vintage from the codebook title): " + "; ".join(
        f"`{v}` {year_ranges(ys)} ({VINTAGE_LABEL.get(vt, vt)})" for (v, vt), ys in sorted(grp.items(), key=lambda kv: min(kv[1]))) + ". "
      f"`OCCALL-EMP.01` before 1994 covers only a first employer other than the current job ({occ_cov}) "
      "and is not used. `CPSOCC80` for 1993 is the occupation of job #1 (the "
      "codebook notes that 14 first jobs are not the CPS job and carry 1970 codes; `CPSOCC70` for 1993 is the CPS item). From "
      "1994 the CPS occupation exists for respondents with an employer for >= 10 hours a week and >= 9 weeks "
      "since the last interview.\n")
    cbk = [r for r in payload["round_map"] if r.get("codebook_refnum")]
    if cbk:
        agree = sum(r.get("vintage_from_bands") == r["vintage"] for r in cbk)
        a(f"Codebook check ({len(cbk)} variable-rounds): the vintage in the variable title agrees with the vintage of "
          "the code bands in the `.cdb` entry (`PROFESSIONAL, TECHNICAL AND KINDRED` = 1970, `MANAGERIAL AND "
          "PROFESSIONAL SPECIALITIES` = 1980, `Business and Financial Operations` second group = 2000 3-digit, "
          f"`Executive, Administrative and Managerial` = 2002 4-digit) for {agree}; the extract's reference number "
          f"equals the codebook's for {sum(r['codebook_refnum'] == r['refnum'] for r in cbk)} and its number of valid "
          f"values equals the codebook's for {sum(r['n_valid_codebook'] == r['n_valid_extract'] for r in cbk)}; "
          f"`ESR_COL` reads 1 = employed, 4 = in active forces in {meta['esr_rounds_confirmed']} of "
          f"{meta['esr_rounds']} rounds. The `.cdb` lists only code bands, not individual codes: the individual "
          "teacher codes below are checked against Attachment 3 and against the two codings of the same jobs.\n")
    rows = []
    for v in ("c70", "c80", "c00", "c02"):
        lab = payload["code_lists"].get(v, {})
        for t in TYPES:
            for code in TEACHER_CODES[v][t]:
                e = lab.get(str(code), {})
                cn = payload["code_counts"]["by_code"].get(v, {}).get(str(code))
                rows.append([VINTAGE_LABEL[v], str(code), TYPE_LABEL[t], e.get("label") or "unverified",
                             str(cn["person_rounds"]) if cn else "-"])
    a("K-12 teacher codes with the official labels of NLSY79 Codebook Supplement Attachment 3 (person-rounds = "
      "block A, ages 25-34, valid occupation; the 2000 and 2002 vintages do not occur at these ages):\n")
    a(ne.tab(rows, ["vintage", "code", "type", "Attachment 3 label", "person-rounds 25-34"]))
    nb = "; ".join(f"{VINTAGE_LABEL[v]}: " + ", ".join(
        f"{c} {payload['code_lists'][v].get(str(c), {}).get('label') or '?'}" for c, _ in NEIGHBOURS[v])
        for v in ("c70", "c80", "c00"))
    a(f"Not counted (neighbouring codes): {nb}.\n")
    tot = cc["k12_both"] + cc["k12_only_1970"] + cc["k12_only_1980"]
    a(f"Concordance of the two codings on the same person-rounds ({cc['n_person_rounds_both_valid']} person-rounds with "
      f"a valid code in both, years {cc['years'][0]}-{cc['years'][-1]}; rows = 1970 code, columns = 1980 code). K-12 under "
      f"both {cc['k12_both']}, only 1970 {cc['k12_only_1970']}, only 1980 {cc['k12_only_1980']}: the codings agree on "
      f"{100 * cc['k12_both'] / tot:.0f}% of the person-rounds that either calls K-12 teaching.\n")
    cols = ["155", "156", "157", "158", "159", "other"]
    rows = [[r] + [str(cc["matrix_rows_1970_cols_1980"][r][k]) for k in cols]
            for r in ("141", "142", "143", "144", "145", "other")]
    a(ne.tab(rows, ["1970 \\ 1980"] + cols))
    a("1970 elementary (142) maps to 1980 elementary (156), pre-K/kindergarten (143) to 155, secondary (144) to 157, "
      "and teachers n.e.c. (145) and adult education (141) to 159; special education (158) has no 1970 code and "
      "comes mostly from 142 and 145.\n")

    # ---- section 2: main tables ---------------------------------------------------------
    a("## 2. Gaps in latent AFQT: teachers minus non-teachers\n")
    a("Main specification: 1980 Census codes, all K-12 teacher codes, a teacher is a K-12 teacher in the majority "
      "of his or her 25-34 rounds with a valid occupation. Latent SDs; mean latent AFQT is measured from the "
      "block-A population mean (0). Counts are unweighted persons.\n")
    rows_spec = [(E, "employed"), (N, "any 25-34 occupation"), ("c80-pop-all-maj", "all respondents observed at 25-34")]
    for comp, title in (("nt", "A. Teachers against all non-teachers"),
                        ("col", "B. All teachers against non-teacher college graduates (16+ years)"),
                        ("cc", "C. College-graduate teachers against non-teacher college graduates")):
        a(f"**{title}**\n")
        a(gap_table(S, C, comp, rows_spec))

    # ---- section 3: sensitivities ---------------------------------------------------------
    a("## 3. Sensitivities\n")
    heads = ["variant", "teachers (pooled)", "A pooled", "A men", "A women", "B pooled", "B men", "B women", "C pooled"]
    for samp, title in (("emp", "Employed sample"), ("any", "Any 25-34 occupation")):
        a(f"**{title}**: gaps in latent SDs; A = against all non-teachers, B = against non-teacher college "
          "graduates, C = college-graduate teachers against non-teacher college graduates.\n")
        names = [f"c80-{samp}-all-maj", f"c80-{samp}-all-half", f"c80-{samp}-all-ever", f"c80-{samp}-core-maj",
                 f"c80-{samp}-core-ever", f"c70-{samp}-all-maj", f"c80-{samp}-all-maj-hgc9",
                 f"c80-{samp}-all-maj-xs", f"c80-{samp}-all-maj-unit", f"c80-{samp}-all-maj-round"]
        rows = []
        for nm in names:
            rows.append([spec_label(SP[nm]) + (" (main)" if nm == f"c80-{samp}-all-maj" else ""),
                         str(n_(nm, "pooled")), g_(nm, "pooled", "nt"), g_(nm, "men", "nt"), g_(nm, "women", "nt"),
                         g_(nm, "pooled", "col"), g_(nm, "men", "col"), g_(nm, "women", "col"), g_(nm, "pooled", "cc")])
        a(ne.tab(rows, heads))
    a("**No occupation required** (every respondent observed at 25-34; a person with no reported occupation is a "
      "non-teacher):\n")
    rows = []
    for rule in ("maj", "half", "ever"):
        nm = f"c80-pop-all-{rule}"
        rows.append([RULE_LABEL[rule], str(n_(nm, "pooled")), g_(nm, "pooled", "nt"), g_(nm, "men", "nt"),
                     g_(nm, "women", "nt"), g_(nm, "pooled", "col"), g_(nm, "men", "col"), g_(nm, "women", "col"),
                     g_(nm, "pooled", "cc")])
    a(ne.tab(rows, heads))
    a(f"Ties: {meta['n_tie_persons']} employed-sample persons have exactly half of their rounds in K-12 teaching "
      "(non-teachers under the majority rule, teachers under 'at least half').\n")

    # ---- section 4: by type -------------------------------------------------------------
    a("## 4. By type of teacher\n")
    a("Employed sample, majority rule within the type (a person is a type-T teacher if more than half of his or "
      "her counted rounds are in type-T codes); the comparison groups are those of the all-K-12 definition, so "
      "persons who teach mostly but not mainly in one type are in neither group. 1980 codes (main):\n")
    heads4 = ["type", "teachers", "of which women", "A pooled", "A women", "B pooled", "B women", "C pooled"]

    def type_rows(vint: str) -> list[list[str]]:
        rows = []
        for t in list(TYPES) + ["core"]:
            if t != "core" and not TEACHER_CODES[vint][t]:
                continue
            nm = f"{vint}-emp-{t}-maj"
            codes = ("+".join(str(cd) for tt in DEFS["core"] for cd in TEACHER_CODES[vint][tt]) if t == "core"
                     else ", ".join(str(cd) for cd in TEACHER_CODES[vint][t]))
            lab = ("core (elementary + secondary" + (" + special education" if vint == "c80" else "") + ")") \
                if t == "core" else TYPE_LABEL[t]
            rows.append([f"{lab} ({codes})", str(n_(nm, "pooled")), str(n_(nm, "women")), g_(nm, "pooled", "nt"),
                         g_(nm, "women", "nt"), g_(nm, "pooled", "col"), g_(nm, "women", "col"), g_(nm, "pooled", "cc")])
        return rows

    a(ne.tab(type_rows("c80"), heads4))
    a("1970 codes (`CPSOCC70` through 1993, `OCCALL-EMP.01` from 1994; no separate special education code):\n")
    a(ne.tab(type_rows("c70"), heads4))

    # ---- section 5: shares --------------------------------------------------------------
    a("## 5. Teaching shares at 25-34 (sanity check against the ACS)\n")
    a("Round shares: weighted person-rounds of the block-A sample in an employed K-12 teacher occupation (1980 "
      "codes) over all person-rounds, employed or not (home production included); the denominator is restricted "
      "to highest grade >= 9 and not armed forces, as in the ACS sample (`acs`), or only to not armed forces "
      "(`no grade rule`). Person shares: persons classified as teachers (employed sample) over persons observed at "
      "25-34 with highest grade >= 9. Percent.\n")
    rows = []
    for g in ("men", "women"):
        gk = {"men": "M", "women": "F"}[g]
        rows.append([g, pct(S, f"share|acs|emp_k12|{g}"), pct(S, f"share|all|emp_k12|{g}"),
                     pct(S, f"share|acs|any_k12|{g}"), pct(S, f"share|acs|emp_core|{g}"),
                     pct(S, f"pshare|maj|{g}"), pct(S, f"pshare|ever|{g}"), pct(S, f"share|acs|emp_k12_given_emp|{g}"),
                     f"{100 * meta['acs_published'][gk]:.2f} / {100 * meta['acs_published_pwgtp'][gk]:.2f}"
                     + (f" / {100 * ac[g]['k12_share']:.2f}" if ac else "")])
    a(ne.tab(rows, ["", "NLSY: employed K-12 / all (acs)", "NLSY: no grade rule", "NLSY: K-12 occupation, any status",
                    "NLSY: core only", "NLSY persons, majority", "NLSY persons, ever", "NLSY: K-12 / employed",
                    "ACS 2009-13: project unweighted / project PWGTP / PWGTP simple rule (here)"]))
    a("By code (employed K-12 / all, acs denominator, percent): " + "; ".join(
        f"{g}: " + ", ".join(f"{TYPE_LABEL[t]} {pct(S, f'share|acs|emp_{t}|{g}')}" for t in TYPES)
        for g in ("men", "women")) + ".\n")
    if ac:
        a("Split by schooling (NLSY: 16+ years at the round; ACS 2009-13 pooled 1-year PUMS in `raw/acs_pums`: "
          "bachelor's or more, simple rule, PWGTP weights; percent). 'NLSY rates at the ACS BA+ share' applies the "
          "NLSY K-12 rates by schooling to the ACS share of bachelor's holders; 'gap closed' is the fraction of the "
          "distance from the NLSY to the ACS simple-rule share that this closes:\n")
        rows = []
        for g in ("men", "women"):
            cf = cfs[g]
            rows.append([g, pct(S, f"share|acs|col|{g}", 1), pct(S, f"share|acs|k12_given_col|{g}", 1),
                         pct(S, f"share|acs|k12_given_noncol|{g}", 2), f"{100 * ac[g]['ba_share']:.1f}",
                         f"{100 * ac[g]['k12_given_ba']:.1f}", f"{100 * ac[g]['k12_given_noba']:.2f}",
                         f"{100 * cf:.2f}", f"{100 * ac[g]['k12_share']:.2f}", f"{100 * closed[g]:.0f}%"])
        a(ne.tab(rows, ["", "NLSY: 16+ share", "NLSY: K-12 given 16+", "NLSY: K-12 given < 16", "ACS: BA+ share",
                        "ACS: K-12 given BA+", "ACS: K-12 given < BA", "NLSY rates at the ACS BA+ share",
                        "ACS simple-rule share", "gap closed"]))
        a(f"For women the education mix closes about {100 * closed['women']:.0f}% of the gap (the K-12 rate among "
          f"college graduates, {pct(S, 'share|acs|k12_given_col|women', 1)}% against {100 * ac['women']['k12_given_ba']:.1f}%, "
          "is close to the ACS); for men it closes about "
          f"{100 * closed['men']:.0f}%, and the rate among college graduates itself is lower "
          f"({pct(S, 'share|acs|k12_given_col|men', 1)}% against {100 * ac['men']['k12_given_ba']:.1f}%). The NLSY79 shares are "
          "not evidence of undercounted teachers: the occupation codes recover the ACS rate for college-educated "
          "women. What remains (the male rate among college graduates and the non-college rates) is not explained "
          "here; it may be a genuine difference between the 1957-64 cohort in 1982-98 and the 2009-13 cohorts, or a "
          "difference between the 1980 and 2002/2010 Census codings of 'other' teachers.\n")

    # ---- caveats ------------------------------------------------------------------------
    a("## Notes and caveats\n")
    men_n = n_(E, "men")
    ever_n, maj_n = n_("c80-emp-all-ever", "pooled"), n_(E, "pooled")
    for ln in [
        "**What the three comparisons mean.** Comparison A is positive because teaching requires a college degree "
        f"({100 * val(S, E + '|pooled|col_share_t'):.0f}% of majority teachers against "
        f"{100 * val(S, E + '|pooled|col_share_nt'):.0f}% of employed non-teachers have 16+ years) and college "
        "graduates score far higher on the AFQT than everyone else. The model has no schooling margin, so its "
        "non-teachers include people who would be non-college workers in the data and `z` carries whatever the "
        "AFQT measures, schooling-related or not; A is the literal counterpart of 'teachers minus non-teachers'; "
        "whether B or C is the better counterpart depends on how much of the schooling difference the model's z is "
        "meant to carry. B and C isolate selection among college graduates, where teachers are below other graduates.",
        "**Timing of the test.** The ASVAB was taken in 1980 at ages 15-23 and is age-normed, not grade-normed, so "
        "it partly reflects schooling attained by the test date (respondents already in college score higher than "
        "their latent ability alone implies). A compares groups with very different final schooling and may "
        "overstate the ability difference; B and C compare groups with (almost) the same final schooling and are less exposed.",
        "**Latent scale.** Dividing a difference of means of `afqt_c` by sqrt(R_comp_pop) assumes classical "
        "measurement error with the same error variance in every group (the convention of the wage slope in "
        "`nlsy_ability.py`).",
        f"**Persistent teachers.** The majority rule keeps persistent teachers: {ever_n} persons teach in at least one "
        f"employed round at 25-34 but only {maj_n} in a majority of them (teaching spells are short for many). "
        f"'Ever a teacher' gives a smaller A ({g_('c80-emp-all-ever', 'pooled', 'nt')} against {g_(E, 'pooled', 'nt')}) "
        f"and a more negative B ({g_('c80-emp-all-ever', 'pooled', 'col')} against {g_(E, 'pooled', 'col')}), with C "
        f"about unchanged ({g_('c80-emp-all-ever', 'pooled', 'cc')} against {g_(E, 'pooled', 'cc')}); the person-round "
        "teaching share in section 5 is the cross-sectional analogue of the model's equilibrium teachers.",
        f"**Men.** Only {men_n} men are majority teachers in the employed sample ({n_('c80-emp-all-ever', 'men')} are "
        "ever teachers), so the male gaps are imprecise; by type the cells are smaller still (section 4).",
        "**Occupation measure.** The CPS job is the current job or, for the non-employed, the most recent job, so "
        "the 'any 25-34 occupation' sample mixes current and past jobs. Occupations are 3-digit Census codes assigned "
        f"by the NLS from verbatim answers; the 1970 and 1980 codings agree on {100 * cc['k12_both'] / tot:.0f}% of "
        "K-12 teacher person-rounds (section 1). Teacher aides, postsecondary and adult education teachers are not "
        "counted; private-school teachers are (occupation-based, like the ACS). 'Teachers n.e.c.' and pre-K/"
        "kindergarten teachers are in the main definition because they are in the ACS definition (2300-2340) that "
        "the model's teaching shares use; sections 3 and 4 show the core and the single types.",
        "**Weights and attrition.** 1979 weights, no correction for later attrition; a person enters the person-level "
        "classification with whatever rounds he or she was observed at 25-34 (1 to 10). SEs cluster on the NLSY79 "
        "household.",
        "**Teaching shares.** NLSY79 person-rounds cover 1982-98 and the 1957-64 birth cohort; the ACS benchmark is "
        "2009-13 and the project's own shares use unweighted persons with hours-based home-production weights, "
        "while the NLSY shares are weighted person-rounds employed at the interview.",
    ]:
        a(f"- {ln}")
    path.write_text("\n".join(L) + "\n")


# =============================================================================
# Self-test and CLI
# =============================================================================

def selftest() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  [{'ok' if cond else 'FAIL'}] {name}")

    nr = np.array([0, 1, 2, 2, 3, 4, 4])
    k = np.array([0, 1, 1, 2, 2, 2, 3])
    check("majority: ties are non-teachers, no rounds is not a teacher",
          classify(nr, k, "maj").tolist() == [False, True, False, True, True, False, True])
    check("half: ties are teachers, no rounds is not a teacher",
          classify(nr, k, "half").tolist() == [False, True, True, True, True, True, True])
    check("ever", classify(nr, k, "ever").tolist() == [False, True, True, True, True, True, True])
    check("vintage titles", vintage_from_title("OCCUPATION AT CURRENT JOB/MOST RECENT JOB (80 CENSUS 3 DIGIT) CPS ITEM") == "c80"
          and vintage_from_title("OCCUPATION (CENSUS 3 DIGIT, 70 CODES) (ALL) JOB #01") == "c70"
          and vintage_from_title("OCCUPATION (70 CENSUS 3 DIGIT) CPS") == "c70"
          and vintage_from_title("OCCUPATION (CENSUS 3 DIGIT, 00 CODES) (ALL)") == "c00"
          and vintage_from_title("OCCUPATION (CENSUS 4 DIGIT, 00 CODES) (ALL) JOB #01") == "c02")
    check("vintage bands", vintage_from_bands("001-195 PROFESSIONAL,TECHNICAL AND KINDRED") == "c70"
          and vintage_from_bands("003-199 MANAGERIAL AND PROFESSIONAL SPECIALITIES") == "c80"
          and vintage_from_bands("10 TO 430: Executive, Administrative and Managerial Occupations") == "c02"
          and vintage_from_bands("Management Business and Financial Operations Computer and Mathematical") == "c00")
    x = np.array([0, 1, 3, 889, 890, 984, 995, 999, np.nan])
    check("valid 1980 codes", valid_code(x, "c80").tolist() == [False, False, True, True, False, False, False, False, False])
    check("valid 1970 codes", valid_code(x, "c70").tolist() == [False, True, True, True, True, True, False, False, False])
    check("code sets", TEACHER_CODES["c70"]["elem"] == (142,) and TEACHER_CODES["c70"]["prek"] == (143,)
          and TEACHER_CODES["c70"]["sec"] == (144,) and all(155 <= cd <= 159 for t in TYPES for cd in TEACHER_CODES["c80"][t])
          and sum(len(TEACHER_CODES["c70"][t]) + len(TEACHER_CODES["c80"][t]) for t in TYPES) == 9)
    lines = ["142        71                             Elementary school teachers",
             "155 Teachers, prekindergarten and kindergarten (231)", "N (156) Teachers, elementary school (232)",
             "  2310     Elementary and Middle School Teachers                             25-2020",
             "144        99                             Something else"]
    check("label parser", _find_label(lines, 142) == "Elementary school teachers"
          and _find_label(lines, 155) == "Teachers, prekindergarten and kindergarten (231)"
          and _find_label(lines, 156) == "Teachers, elementary school (232)"
          and _find_label(lines, 2310) == "Elementary and Middle School Teachers"
          and _find_label(lines, 144) is None)
    names = [s.name for s in make_specs()]
    check("spec names unique", len(names) == len(set(names)))
    print("selftest", "passed" if ok else "FAILED")
    return 0 if ok else 1


def read_json(path: Path) -> Any:
    return json.loads(path.read_text())


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bootstrap", type=int, default=B_DEFAULT)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--quick", action="store_true", help="B = 50")
    p.add_argument("--no-codebook", action="store_true", help="skip the scan of the NLS codebook (.cdb)")
    p.add_argument("--md-only", action="store_true", help="rewrite the .md from the saved .json")
    p.add_argument("--selftest", action="store_true")
    p.add_argument("--indir", default=str(HERE))
    p.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    if args.selftest:
        return selftest()
    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "estimates" / "settings"
    outdir.mkdir(parents=True, exist_ok=True)
    if args.md_only:
        write_md(read_json(outdir / "nlsy_teachers.json"), outdir / "nlsy_teachers.md")
        return 0
    B = 50 if args.quick else args.bootstrap
    t0 = time.time()

    c = na.Cohort("nlsy79", indir)
    fr = Frame(c)
    # ---- codebook: vintages, ESR_COL codes ----------------------------------------------
    cdb = None
    esr_rounds = esr_ok = 0
    zip_path = indir / "raw" / "_downloads" / "nlsy79_all_1979-2022.zip"
    if not args.no_codebook and zip_path.exists():
        LOG.info("scanning the codebook in %s", zip_path.name)
        cdb = read_codebook(zip_path, {"CPSOCC70", "CPSOCC80", "OCCALL-EMP.01", "ESR_COL"})
        esr = {y: e["codes"] for (q, y), e in cdb.items() if q == "ESR_COL"}
        esr_rounds = len(esr)
        esr_ok = sum("EMPLOYED" in cd.get(1, "").upper() and "UN" not in cd.get(1, "").upper()
                     and "FORCES" in cd.get(4, "").upper() for cd in esr.values())
        LOG.info("  ESR_COL labels (1 = employed, 4 = in active forces) confirmed in %d of %d rounds", esr_ok, esr_rounds)
    rmap = round_map(fr, c.meta[c.meta.qname.isin(["CPSOCC70", "CPSOCC80", "OCCALL-EMP.01"])], cdb)
    mism = [(r["variable"], r["year"]) for r in rmap if r.get("vintage_from_bands") not in (None, r["vintage"])]
    if mism:
        LOG.warning("  vintage from title and from codebook bands disagree for: %s", mism)
    att3 = attachment3_labels(indir / "raw" / "_downloads" / "nlsy79_attachment3")
    code_lists: dict[str, dict[str, dict[str, Any]]] = {v: {} for v in TEACHER_CODES}
    for v in TEACHER_CODES:
        for t in TYPES:
            for code in TEACHER_CODES[v][t]:
                code_lists[v][str(code)] = {"label": (att3 or {}).get(v, {}).get(code), "type": t}
        for code, why in NEIGHBOURS[v]:
            code_lists[v][str(code)] = {"label": (att3 or {}).get(v, {}).get(code), "type": None, "not_k12": why}
    cc = code_counts(fr)
    LOG.info("  concordance: K-12 under both codings %d, only 1970 %d, only 1980 %d",
             cc["concordance"]["k12_both"], cc["concordance"]["k12_only_1970"], cc["concordance"]["k12_only_1980"])

    # ---- estimates ------------------------------------------------------------------------
    res = run_bootstrap(fr, B, args.seed)
    S = res["stats"]
    nan = float("nan")
    acs_pub = acs_pub_w = {"M": nan, "F": nan}
    try:
        aj = read_json(indir / "estimates" / "acs_occupations.json")
        key = "Kindergarten - Secondary Teachers"
        acs_pub = {g: aj["samples"]["acs5_2009_13"]["shares"][g][key] for g in ("M", "F")}
        acs_pub_w = {g: aj["samples"]["acs5_2009_13"]["shares_pwgtp_weighted"][g][key] for g in ("M", "F")}
    except Exception as exc:                                    # pragma: no cover
        LOG.warning("  ACS shares not read: %s", exc)
    sz = None
    try:
        sz = next(e["value"] for e in read_json(indir / "estimates" / "nlsy_ability.json")["estimates"]
                  if e["key"] == "implied_s_z")
    except Exception as exc:                                    # pragma: no cover
        LOG.warning("  s_z not read: %s", exc)
    Rm = fr.R[("c80", "emp")]
    k_all = sum(Rm["nt"][t] for t in TYPES)
    meta = {"generated": dt.datetime.now().isoformat(timespec="seconds"),
            "script": "data/spatial/nlsy_teachers.py", "B": B, "seed": args.seed,
            "bootstrap": "NLSY79 households (HHID), cluster multiplicities from rng.integers as in "
                         "nlsy_extensions.t3_estimate; AFQT standardization and R_comp_pop re-estimated in every draw",
            "age_window": [na.AGE_LO, na.AGE_HI], "years": sorted(set(fr.py.Y.tolist())),
            "weights": "1979 SAMPWEIGHT / 100, one observation per person (nlsy_ability.py block C convention); "
                       "mean round weights as a sensitivity",
            "n_block_a": int(fr.D.asamp.sum()), "n_persons_observed_25_34": int((fr.observed & fr.D.asamp).sum()),
            "n_person_rounds_25_34_block_a": int(fr.pop_py.sum()), "R_comp_pop": S["R_comp_pop"]["value"],
            "R_comp_pop_se": S["R_comp_pop"]["se"], "s_z": sz, "acs_published": acs_pub,
            "acs_published_pwgtp": acs_pub_w, "failed_draws": res["failed_draws"],
            "esr_rounds": esr_rounds, "esr_rounds_confirmed": esr_ok,
            "n_persons_with_occupation": {f"{v}_{s}": int((R["nr"] > 0).sum()) for (v, s), R in fr.R.items()},
            "n_tie_persons": int(np.sum((Rm["nr"] > 0) & (2 * k_all == Rm["nr"]) & (k_all > 0)))}
    headline = {}
    for samp, spec in (("employed", "c80-emp-all-maj"), ("any_occupation", "c80-any-all-maj"),
                       ("employed_ever", "c80-emp-all-ever")):
        headline[samp] = {g: {cl: {**{k: S[f"{spec}|{g}|{cl}|gap"][k] for k in ("value", "se", "ci95_lo", "ci95_hi")},
                                   **res["counts"][spec][g][cl]} for cl in COMPS} for g in GROUPS}
    payload = {"meta": meta, "headline_gap_latent_sd": headline,
               "definitions": {"vintages": VINTAGE_LABEL, "valid_code_range": VALID_RANGE,
                               "teacher_codes": {v: {t: list(cs) for t, cs in d.items()} for v, d in TEACHER_CODES.items()},
                               "definitions": {k: list(v) for k, v in DEFS.items()}, "rules": RULE_LABEL,
                               "samples": SAMPLE_LABEL,
                               "comparisons": {"nt": "A: teachers vs all non-teachers",
                                               "col": "B: all teachers vs non-teachers with 16+ years",
                                               "cc": "C: teachers with 16+ years vs non-teachers with 16+ years"},
                               "attachment3_sources": {v: ATT3_BASE + f + ".pdf" for v, f in ATT3_FILES.items()},
                               "statistic_keys": {
                                   "spec|group|comparison|stat": "spec = a key of `specs` (vintage-sample-definition-rule[-hgc9][-weights]); "
                                       "group = pooled/men/women; comparison = nt/col/cc; stat = mean_t, mean_c (mean latent AFQT of "
                                       "teachers / comparison group, latent SDs from the block-A mean), gap (= mean_t - mean_c, latent SDs), "
                                       "gap_obs (the same in afqt_c units; gap = gap_obs / sqrt(R_comp_pop))",
                                   "spec|group|share": "weighted share of teachers among the persons of the spec classified as teacher or comparison",
                                   "spec|group|col_share_t / col_share_nt": "weighted share with 16+ years among teachers / among non-teachers",
                                   "share|{acs,all}|{measure}|group": "teaching shares over person-rounds at 25-34 (acs = highest grade >= 9 "
                                       "and not armed forces; all = not armed forces); measure emp_k12 = employed K-12 teachers over all "
                                       "person-rounds, any_k12 = K-12 occupation whatever the status, emp_core, emp_k12_c70 (1970 codes), "
                                       "emp_{prek,elem,sec,sped,nec}, emp_rate, col (share with 16+ years), emp_k12_given_emp, "
                                       "k12_given_col, k12_given_noncol",
                                   "pshare|rule|group": "persons classified as K-12 teachers (employed sample, 1980 codes) over persons "
                                       "observed at 25-34 with highest grade >= 9"}},
               "round_map": rmap, "code_lists": code_lists, "code_counts": cc,
               "specs": {sp.name: sp.describe() for sp in fr.specs},
               "counts": res["counts"], "statistics": S, "acs_benchmark": acs_benchmark(indir)}
    payload = clean(payload)
    (outdir / "nlsy_teachers.json").write_text(json.dumps(payload, indent=1, default=float, allow_nan=False))
    write_md(payload, outdir / "nlsy_teachers.md")
    LOG.info("done in %.1fs -> %s", time.time() - t0, outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
