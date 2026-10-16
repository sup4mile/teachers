#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
teacher_residence.py -- evidence to discipline a teacher RESIDENCE / WORKPLACE split in the spatial
model (DESIGN.md section M2, calibration_log.md section 5 item 2).

In the bundled model a teacher lives where she teaches.  The M2 extension lets a teacher choose a
residence r and a workplace w, pays a commuting cost when r != w (a goods cost k_rw = r_k * Cbar in
M2_SPEC.md; a proportional time cost, income x (1 - tau_k), in M2_SPEC_v2.md), and adds three data
objects: the gap in teachers by RESIDENCE against the gap by WORKPLACE (`net_commute_gap`), the
residence-based teacher share of employed residents (`gap_teach_res_emp`), and the commuting cost
(r_k, tau_k).  This script builds them.

Naming.  The model's `gap_teach_res_emp` and `net_commute_gap` are about K-12 teachers.  The ACS
school-district table has only the whole education, training and library group (EDU, SOC 25-0000), so
the district gaps are EDU PROXIES, named `gap_edu_res_emp` and `net_commute_gap_edu`; the K-12
versions (`k12.gap_teach_res_emp`, `k12.net_commute_gap`) add a PUMS composition term (the share of EDU
workers who are public K-12 teachers, suburban against city PUMAs of the same commuting zone) to them,
with a joint CZ bootstrap and a specification range.  Residents are counted once (the elementary and
secondary ACS layers overlap; `variant_panel`), which matters for counts but not for ratios.

(a) Residence-based teachers.  ACS 5-year 2014-18 (vintage 2018) school-district table C24010 (Sex by
    Occupation for the Civilian Employed Population 16 Years and Over): total employed residents and
    "Educational instruction, and library occupations" (the 2018 vintage's label for SOC 25-0000,
    Education, Training, and Library; men + women).  The three school-district layers (unified,
    elementary, secondary) are pulled state by state as `acs_earnings.py` / `acs_housing.py` do and
    merged to the district panel on `leaid` as `district_bootstrap.load_income_columns` merges the
    B20005 earnings tables ("Remainder of <state>" rows dropped, one row per leaid, annotation codes
    masked).  That rule is right for ratios and medians but not for COUNTS: a resident of an
    elementary district that lies inside a secondary district appears in both records and would be
    counted twice when both districts sit in the same location.  The headline therefore counts
    residents once (unified + elementary records; the secondary LEAs keep their FTE and pupils but
    contribute no residents; the unified and elementary records tile the states); the earnings rule,
    residents from unified + secondary records, unified districts only, ... are variants.
(b) Workplace-based teachers.  K-12 FTE teachers (`teachers_k12`: pre-K teachers removed, the audited
    variant behind the `gap_teachers_pp` target) and K-12 pupils (`pupils_k12`) of the audited panel
    of `district_settings.py`, New York City as one district.
(c) Locale-partition moments on the 188 CZs of the calibration, whole-CZ bootstrap with the SAME
    draws as `district_settings.py` (B = 1999, seed 20260929, one multiplicity matrix shared by every
    statistic, so correlations across moments are valid):
        gap_edu_res_emp      Delta log (EDU residents / employed residents)      (EDU proxy)
        gap_edu_res          Delta log EDU residents
        gap_fte              Delta log K-12 FTE teachers
        net_commute_gap_edu  gap_edu_res - gap_fte = Delta log (EDU residents per FTE)   (EDU proxy)
        EDU residents per FTE by location (levels), EDU residents per pupil gap
    on the districts that have an ACS record AND an audited FTE and pupil count, with variants (the
    earnings rule, residents from the secondary layer, unified districts only, New York City dropped,
    New York CZ left out, each statistic on its own sample) and the ACS sampling error of C24010
    (margins of error).  The four adopted targets (0.3764, 0.0086, -0.0177, 0.1122) and their SEs are
    reproduced as a check and the script stops if they are not.
(d) Composition of SOC 25-0000.  ACS PUMS 2018 1-year (Census API, one state at a time, the 80
    replicate weights included for the successive-difference SEs): the share of K-12 teachers
    (elementary/middle, secondary, special education) within education, training and library
    occupations, nationally, across PUMAs (dispersion net of sampling noise) and by principal-city
    versus suburban PUMAs, classified from the Census 2010 PUMA equivalency files (2010 population of
    each place within each PUMA) and the OMB September-2018 list of principal cities of metropolitan
    areas (cuts varied).  The composition term, computed CZ by CZ, converts the EDU gaps of (c) into the
    K-12 counterparts of the model moments (`k12_block`): joint CZ bootstrap of the district gap and
    the composition term on the draws of (c), plus the survey (replicate-weight) and ACS sampling
    variances, and the range over specifications (class cuts, sectors, group quarters, overlap rule).
(e) Commuting.  Mean one-way commute (JWMNP) of K-12 teachers and other workers 25-64 by sex, the share
    commuting more than 30 minutes (and 30 minutes or more: 12% to 14% of commuters report exactly
    30), the share working outside the residence POWPUMA (Census PUMA -> POWPUMA conversion of the 2018
    PUMS code lists), the extra commute of working elsewhere from a regression with residence-PUMA
    fixed effects (crossing a POWPUMA; working in a POWPUMA of the other class, city or suburb, with
    residence and workplace on the SAME POWPUMA geography), teachers by residence against by workplace
    on one geography, and a check of the PUMS against the New York City IBO count of teachers living in
    the city.  Literature (Boyd et al. 2005, Reininger 2012, NYC IBO, Santelli and Grissom 2024, USDOT
    value of travel time, NCES school year length) is carried with its source.
(f) Model units.  The extra commuting implied by working in the other location costs a worker an annual
    time cost.  r_k (goods cost, M2_SPEC.md) = that cost / Cbar, in two units: per CE consumer unit
    (Cbar = CE 2018 mean annual expenditures per consumer unit, $61,224.13, as for r_p in
    `housing_costs.py`, the unit the task fixed) and per ADULT (the model's agents are single adults:
    Cbar / 1.9 adults per consumer unit = $32,223; r_p is a household cost over household spending
    and unit-free, a per-worker cost is not).  tau_k (proportional time cost, M2_SPEC_v2.md) = the
    time cost over the worker's annual earnings = (share of the own wage at which time is valued) x
    extra hours / annual hours, for teachers (185 days) and for other workers (250 days).  Central
    values, ranges and one-at-a-time sensitivities.

Usage
    set -a; source .env; set +a                  # CENSUS_API_KEY (or let the script read .env)
    python3 teacher_residence.py fetch           # C24010 tables, PUMS 2018 1-year (about 6 minutes),
                                                 # PUMA / OMB files; cached under raw/
    python3 teacher_residence.py fetch --only acs        # or: pums, geo
    python3 teacher_residence.py build           # B = 1999, seed 20260929; about 3-4 minutes
    python3 teacher_residence.py build --B 199 --outdir <dir>   # quick check
    python3 teacher_residence.py                 # fetch (cached downloads are skipped) and build
    python3 teacher_residence.py selftest        # unit checks of the helpers (no data, no network)
Writes estimates/settings/teacher_residence.{md,json}, teacher_residence_targets.toml and
teacher_residence_draws.csv.  Run on a compute node (the build peaks at about 11 GB of memory).

The Census API key is read from the environment or from data/spatial/.env and is never printed or
written; error strings are scrubbed of it.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acs_earnings as ae  # noqa: E402  (school-district pull machinery; key scrubbing lives there)
import acs_occupations as ao  # noqa: E402  (PUMS API wrapper)
import data_estimate as de  # noqa: E402
import district_bootstrap as db  # noqa: E402
import district_settings as ds  # noqa: E402
import housing_costs as hc  # noqa: E402  (CE 2018 denominator of r_p)

LOG = logging.getLogger("teacher_residence")

YEAR = 2018
ACS_MOE_Z = 1.645                      # ACS margins of error are at the 90% level
RAW = HERE / "raw"
PUMS_DIR = ao.RAW / "acs1_2018_workers"          # raw/acs_pums/acs1_2018_workers/<ST>.parquet
GEO_DIR = RAW / "puma_geo"
ACS_C24010 = RAW / f"acs_c24010_school_districts_{YEAR}.parquet"
ACS_C24010_ST = RAW / f"acs_c24010_states_{YEAR}.parquet"
SEED = db.SEED
B_DEFAULT = db.B_DEFAULT

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------

# ACS 5-year 2014-18 table C24010 (Sex by Occupation for the Civilian Employed Population
# 16 Years and Over).  Vintages through 2017 label C24010_014 / _050 (the SOC 25-0000 group) "Education,
# training, and library occupations"; the 2018 vintage and later (SOC 2018 renamed the major group)
# label them "Educational instruction, and library occupations" (api.census.gov variables.json).
C24010 = {
    "C24010_001": "total, civilian employed 16+",
    "C24010_002": "male, civilian employed 16+",
    "C24010_038": "female, civilian employed 16+",
    "C24010_014": "male, Educational instruction, and library occupations",
    "C24010_050": "female, Educational instruction, and library occupations",
}

CENSUS = "https://www2.census.gov"
PUMSEQ_URL = CENSUS + "/geo/docs/reference/puma/PUMSEQ10_{st}.txt"
OMB_BASE = CENSUS + "/programs-surveys/metro-micro/geographies/reference-files/2018/delineation-files/"
OMB_FILES = {"list1": "list1_Sep_2018.xls", "list2": "list2_Sep_2018.xls"}
PUMS_CODELIST_URL = CENSUS + "/programs-surveys/acs/tech_docs/pums/code_lists/ACSPUMS2018CodeLists.xls"
UA = "Mozilla/5.0 (research; teachers-spatial-calibration; jpryan7@wisc.edu)"

# PUMS 2018 1-year person variables.  `core` is requested together with the first 28 replicate
# weights; the other 52 follow in two more calls keyed on (SERIALNO, SPORDER) (the API takes at
# most 50 variables per call).
PUMS_CORE = ["SERIALNO", "SPORDER", "PUMA", "PWGTP", "AGEP", "SEX", "ESR", "COW", "OCCP", "SOCP",
             "INDP", "POWSP", "POWPUMA", "JWMNP", "JWTR", "SCHL", "WAGP", "WKHP", "WKW", "ADJINC"]
N_REP = 80
REP = [f"PWGTP{i}" for i in range(1, N_REP + 1)]
PUMS_CHUNKS = [PUMS_CORE + REP[:28], ["SERIALNO", "SPORDER"] + REP[28:76],
               ["SERIALNO", "SPORDER"] + REP[76:]]
STR_COLS = ("SERIALNO", "SOCP")

# 2018-based Census occupation codes in SOC 25-0000 (education, training and library), PUMS
# OCCP labels "EDU-...": SOCP in parentheses.
EDU_CODES = {
    2205: "postsecondary teachers (25-1000)",
    2300: "preschool and kindergarten teachers (25-2010)",
    2310: "elementary and middle school teachers (25-2020)",
    2320: "secondary school teachers (25-2030)",
    2330: "special education teachers (25-2050)",
    2350: "tutors (25-3041)",
    2360: "other teachers and instructors (25-30XX)",
    2400: "archivists, curators and museum technicians (25-4010)",
    2435: "librarians and media collections specialists (25-4022)",
    2440: "library technicians (25-4031)",
    2545: "teaching assistants (25-9040)",
    2555: "other educational instruction and library workers (25-90XX)",
}
K12_CODES = (2310, 2320, 2330)           # elementary/middle, secondary, special education
K12_PLUS = (2300,) + K12_CODES           # ... plus preschool and kindergarten teachers

# BLS CPI-U, all items, U.S. city average, annual averages (1982-84 = 100)
CPI_U = {2015: 237.017, 2018: 251.107}

# BLS Consumer Expenditure Survey 2018, Table 1400 (size of consumer unit), all consumer units: 2.5 people,
# 0.6 children under 18, 0.4 adults 65 and older, 1.3 earners, average annual expenditures $61,224
# (published to one decimal, so adults = 2.5 - 0.6 = 1.9 carries a rounding range of +-0.1).
CE_PERSONS, CE_CHILDREN = 2.5, 0.6
CE_ADULTS = round(CE_PERSONS - CE_CHILDREN, 1)
CE_ADULTS_RANGE = (round(CE_ADULTS - 0.1, 1), round(CE_ADULTS + 0.1, 1))
CE_URL = ("https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error/"
          "cu-size-2018.xlsx")

# commuting days a year: teachers (NCES 179 school days, 2011-12 SASS, plus about a week of contract days;
# assumption) and other workers (50 weeks of 5 days, less holidays and days at home; assumption)
NASHVILLE_EXTRA_MIN = 13.8             # Santelli and Grissom (2024, WP p. 15): out-of-county teachers' commute, minutes longer
DAYS_TEACHER = dict(low=180.0, central=185.0, high=190.0)
DAYS_OTHER = dict(low=240.0, central=250.0, high=260.0)

# USDOT, "The Value of Travel Time Savings: Departmental Guidance for Conducting Economic
# Evaluations, Revision 2 (2016 Update)", Sept. 27, 2016: Table 1 (50% of hourly income, local
# personal travel), Table 2 (plausible range 35%-60%), Table 3 ($27.20 an hour: 2015 median
# household income $56,516 / 2,080), Table 4 ($13.60 local personal), Table 5 ($9.50-$16.30).
VTT_2015 = {"central": 13.60, "low": 9.50, "high": 16.30, "income_hour": 27.20,
            "share_central": 0.50, "share_low": 0.35, "share_high": 0.60}
USDOT_URL = ("https://www.transportation.gov/sites/dot.gov/files/docs/"
             "2016%20Revised%20Value%20of%20Travel%20Time%20Guidance.pdf")


def _key() -> str:
    return ao._read_key()


# --------------------------------------------------------------------------
# Fetch: ACS C24010 by school district
# --------------------------------------------------------------------------

def fetch_acs(year: int = YEAR, workers: int = 8) -> Path:
    """C24010 for the unified, elementary and secondary school-district layers, one state per
    request, with margins of error; raw/acs_c24010_school_districts_<year>.parquet."""
    key = _key()
    sess = ae.make_session()
    base = f"https://api.census.gov/data/{year}/acs/acs5"
    r = sess.get(f"{base}?get=NAME&for=state:*&key={key}", timeout=60)
    r.raise_for_status()
    states = [row[-1] for row in r.json()[1:]]
    getvars = [f"{k}{s}" for k in C24010 for s in ("E", "M")]
    out = ae.pull(sess, base, states, key, ",".join(["NAME"] + getvars), workers, "C24010")
    if out.empty:
        raise SystemExit("C24010: nothing returned")
    out["year"] = year
    for c in out.columns:
        if c[:1] == "C" and c[-1] in "EM":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    ACS_C24010.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(ACS_C24010, index=False)
    LOG.info("wrote %s (%d rows x %d cols); layers: %s", ACS_C24010.name, len(out), out.shape[1],
             out.sd_layer.value_counts().to_dict())
    # state totals of the same table: the check that unified + elementary records tile each state
    r = sess.get(f"{base}?get=NAME,{','.join(k + 'E' for k in C24010)}&for=state:*&key={key}", timeout=120)
    r.raise_for_status()
    j = r.json()
    st = pd.DataFrame(j[1:], columns=j[0])
    for c in st.columns:
        if c[:1] == "C" and c[-1] == "E":
            st[c] = pd.to_numeric(st[c], errors="coerce")
    st.to_parquet(ACS_C24010_ST, index=False)
    LOG.info("wrote %s (%d states)", ACS_C24010_ST.name, len(st))
    return ACS_C24010


# --------------------------------------------------------------------------
# Fetch: ACS PUMS 2018 1-year, employed persons 16+, with replicate weights
# --------------------------------------------------------------------------

def _pums_state(st: int) -> tuple[int, int, str]:
    """One state: the three chunks merged on (SERIALNO, SPORDER) and cached as parquet.
    Returns (state, rows, message)."""
    f = PUMS_DIR / f"{st:02d}.parquet"
    if f.exists():
        have = set(pd.read_parquet(f).columns)
        if set(PUMS_CORE) | set(REP) <= have:
            return st, -1, "cached"
    api = ao.Api()
    url = ao.API.format(year=YEAR, product="acs1")
    frames = []
    for i, chunk in enumerate(PUMS_CHUNKS):
        j = api.get_json(url, {"get": ",".join(chunk), "for": f"state:{st:02d}",
                               "ESR": "1:2", "AGEP": "16:99"})
        if not j:
            return st, 0, f"chunk {i}: empty response"
        d = pd.DataFrame(j[1:], columns=j[0])
        d = d.loc[:, ~d.columns.duplicated()]                 # predicates come back twice
        frames.append(d[[c for c in chunk if c in d.columns]])   # drop the predicate / state columns
    d = frames[0]
    for e in frames[1:]:
        if len(e) != len(d):
            return st, 0, f"chunk row counts differ: {len(d)} vs {len(e)}"
        d = d.merge(e, on=["SERIALNO", "SPORDER"], how="inner", validate="one_to_one")
    if len(d) != len(frames[0]):
        return st, 0, "merge on (SERIALNO, SPORDER) lost rows"
    d["ST"] = f"{st:02d}"
    for c in d.columns:
        if c not in STR_COLS and c != "ST":
            d[c] = pd.to_numeric(d[c], errors="coerce")
    for c in REP + ["PWGTP"]:
        d[c] = d[c].astype(np.int32)
    PUMS_DIR.mkdir(parents=True, exist_ok=True)
    d.to_parquet(f, index=False)
    return st, len(d), "ok"


def fetch_pums(workers: int = 3) -> Path:
    """Employed (ESR 1-2) persons 16+ of the ACS 2018 1-year PUMS, all states and DC, cached
    state by state under raw/acs_pums/acs1_2018_workers/."""
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for st, n, msg in pool.map(_pums_state, ao.STATES):
            done += 1
            LOG.info("PUMS 2018 state %02d: %s%s  [%d/%d, %.0fs]", st,
                     msg, "" if n < 0 else f" ({n} rows)", done, len(ao.STATES), time.time() - t0)
            if msg not in ("ok", "cached"):
                raise SystemExit(f"PUMS state {st:02d} failed: {msg}")
    return PUMS_DIR


# --------------------------------------------------------------------------
# Fetch: PUMA equivalency files, OMB delineation lists, PUMS code lists
# --------------------------------------------------------------------------

def _download(url: str, dest: Path, tries: int = 4) -> None:
    import requests
    if dest.exists() and dest.stat().st_size > 0:
        return
    last = ""
    for k in range(tries):
        try:
            r = requests.get(url, headers={"User-Agent": UA}, timeout=180)
            if r.status_code == 200 and len(r.content) > 0:
                dest.write_bytes(r.content)
                return
            last = f"http {r.status_code}"
        except Exception as exc:                                  # noqa: BLE001
            last = f"{type(exc).__name__}: {exc}"
        time.sleep(3 * (k + 1))
    raise SystemExit(f"could not download {url}: {last}")


def fetch_geo() -> Path:
    GEO_DIR.mkdir(parents=True, exist_ok=True)
    for st in ao.STATES:
        _download(PUMSEQ_URL.format(st=f"{st:02d}"), GEO_DIR / f"PUMSEQ10_{st:02d}.txt")
    for k, name in OMB_FILES.items():
        _download(OMB_BASE + name, GEO_DIR / name)
    _download(PUMS_CODELIST_URL, GEO_DIR / "ACSPUMS2018CodeLists.xls")
    LOG.info("geography files in %s", GEO_DIR)
    return GEO_DIR


# --------------------------------------------------------------------------
# Districts (a)-(c): residence-based EDU workers against workplace-based K-12 FTE
# --------------------------------------------------------------------------

MATCH = ("m_edu", "m_emp", "m_fte", "m_pup")       # the matched sample's columns
VARS = ("m_vedu", "m_vemp")                         # sampling variances of the two resident counts (ACS MOEs)

# Gaps computed by `data_estimate.czone_gaps` on the matched columns (location totals, ratios of
# location totals, then log(location 2 / location 1) within the CZ), exactly as the audited moments.
GAPS_TR = [
    de.Gap("tr_teach_emp", "EDU residents per employed resident", "n_teach_res / residents at work",
           "ratio", num="m_edu", den="m_emp", contrast="log"),
    de.Gap("tr_edu_res", "EDU residents", "n_teach_res", "total", num="m_edu"),
    de.Gap("tr_emp_res", "employed residents", "residents at work", "total", num="m_emp"),
    de.Gap("tr_fte", "K-12 FTE teachers (workplace)", "n_teach_work", "total", num="m_fte"),
    de.Gap("tr_pup", "K-12 pupils", "M", "total", num="m_pup"),
    de.Gap("tr_edu_fte", "EDU residents per K-12 FTE", "n_teach_res / n_teach_work", "ratio",
           num="m_edu", den="m_fte", contrast="log"),
    de.Gap("tr_edu_pup", "EDU residents per K-12 pupil", "n_teach_res / M", "ratio",
           num="m_edu", den="m_pup", contrast="log"),
    de.Gap("tr_fte_pup", "K-12 FTE per K-12 pupil (matched)", "n_teach_work / M", "ratio",
           num="m_fte", den="m_pup", contrast="log"),
    # each statistic on its own maximal sample (not the matched one): variant `unmatched`
    de.Gap("tr_edu_all", "EDU residents, every district with a C24010 record", "n_teach_res",
           "total", num="a_edu"),
    de.Gap("tr_fte_all", "K-12 FTE, every district with an audited FTE", "n_teach_work", "total",
           num="a_fte"),
    # location totals of the ACS sampling variances (used for the ACS sampling SE of the log gaps)
    de.Gap("tr_vedu", "ACS sampling variance of the EDU count", "-", "total", num="m_vedu"),
    de.Gap("tr_vemp", "ACS sampling variance of the employed count", "-", "total", num="m_vemp"),
]

# The statistics reported: key -> (CZ-table column, kind, label).  The district tables only have the whole
# education, training and library group (EDU, SOC 25-0000), so these are EDU PROXIES of the model's
# `gap_teach_res_emp` and `net_commute_gap` (K-12 teachers), which `k12_block` builds with a PUMS
# composition term; the proxies carry their own names so that they cannot be mistaken for the model's.
STATS_TR: list[tuple[str, str, str, str]] = [
    ("gap_edu_res_emp", "gap_tr_teach_emp", "log",
     "EDU proxy for gap_teach_res_emp (NOT the K-12 model counterpart): Delta log (EDU residents / employed residents)"),
    ("gap_edu_res", "gap_tr_edu_res", "log", "Delta log EDU residents (workers by residence)"),
    ("gap_fte", "gap_tr_fte", "log", "Delta log K-12 FTE teachers (workers by workplace), matched sample"),
    ("net_commute_gap_edu", "gap_tr_edu_fte", "log",
     "EDU proxy for net_commute_gap (NOT the K-12 model counterpart): gap_edu_res - gap_fte = Delta log (EDU residents per K-12 FTE)"),
    ("gap_edu_per_pupil", "gap_tr_edu_pup", "log", "Delta log (EDU residents per K-12 pupil)"),
    ("gap_emp_res", "gap_tr_emp_res", "log", "Delta log employed residents"),
    ("gap_pupils_k12", "gap_tr_pup", "log", "Delta log K-12 pupils, matched sample"),
    ("gap_fte_per_pupil", "gap_tr_fte_pup", "log",
     "Delta log (K-12 FTE per K-12 pupil), matched sample (the gap_teachers_pp definition)"),
]
LEVELS_TR: list[tuple[str, str, str]] = [
    ("edu_per_fte", "tr_edu_fte", "EDU residents per K-12 FTE teacher"),
    ("edu_share_emp", "tr_teach_emp", "EDU residents / employed residents"),
    ("edu_per_pupil", "tr_edu_pup", "EDU residents per K-12 pupil"),
    ("fte_per_pupil", "tr_fte_pup", "K-12 FTE per K-12 pupil"),
]
CORR_TARGETS = ["gap_pupils", "gap_salary", "gap_teachers_pp", "gap_income"]


def load_c24010() -> pd.DataFrame:
    """C24010 records, one row per leaid, with the rule of `district_bootstrap.load_income_columns`:
    'Remainder of <state>' rows (code 99999) are not districts and are dropped, one row per leaid
    (the three layers otherwise never share an id), annotation codes (negative) are masked."""
    if not ACS_C24010.exists():
        raise SystemExit(f"no {ACS_C24010.name}: run `teacher_residence.py fetch --only acs`")
    a = pd.read_parquet(ACS_C24010)
    n_all = len(a)
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid")

    def pos(col: str, lo: float) -> np.ndarray:
        x = a[col].astype(float)
        return x.where(x >= lo).to_numpy()

    out = pd.DataFrame({"leaid": a.leaid.to_numpy(),
                        "sd_layer": a.sd_layer.str.extract(r"\((\w+)\)")[0].to_numpy()})
    out["emp_res"] = pos("C24010_001E", 1)             # a district with no employed resident has no ratio
    out["emp_res_m"], out["emp_res_f"] = pos("C24010_002E", 0), pos("C24010_038E", 0)
    out["edu_res_m"], out["edu_res_f"] = pos("C24010_014E", 0), pos("C24010_050E", 0)
    out["edu_res"] = out.edu_res_m + out.edu_res_f
    # sampling variance of the two counts from the published 90% margins of error (SE = MOE / 1.645; the
    # men's and women's EDU counts are independent estimates, so the variances add; no annotation codes
    # occur in this table, but negative values would be masked to zero variance)
    moe = lambda col: np.nan_to_num(a[col].astype(float).where(a[col].astype(float) >= 0).to_numpy(), nan=0.0)   # noqa: E731
    out["edu_var"] = (moe("C24010_014M") ** 2 + moe("C24010_050M") ** 2) / ACS_MOE_Z ** 2
    out["emp_var"] = moe("C24010_001M") ** 2 / ACS_MOE_Z ** 2
    LOG.info("C24010: %d records, %d after dropping 'Remainder' rows; layers %s", n_all, len(out),
             out.sd_layer.value_counts().to_dict())
    return out


def tiling_check() -> dict[str, Any]:
    """Do the unified and elementary records of C24010 tile the states?  Their employed residents against the
    state totals of the same table, and with the secondary layer added (the overlap of the earnings rule)."""
    if not ACS_C24010_ST.exists():
        LOG.warning("no %s: tiling check skipped (re-run `fetch --only acs`)", ACS_C24010_ST.name)
        return {}
    a = pd.read_parquet(ACS_C24010)
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid")
    a["layer"] = a.sd_layer.str.extract(r"\((\w+)\)")[0]
    tot = a.groupby(["state", "layer"]).C24010_001E.sum().unstack("layer").fillna(0.0)
    st = pd.read_parquet(ACS_C24010_ST).set_index("state")
    tab = tot.join(st[["NAME", "C24010_001E"]].rename(columns={"C24010_001E": "state_total"}))
    tab["ue"] = (tab.unified + tab.elementary) / tab.state_total
    tab["all"] = (tab.unified + tab.elementary + tab.secondary) / tab.state_total
    nat_ue = float((tot.unified.sum() + tot.elementary.sum()) / st.C24010_001E.sum())
    nat_all = float((tot.unified.sum() + tot.elementary.sum() + tot.secondary.sum()) / st.C24010_001E.sum())
    top = tab.sort_values("all", ascending=False).head(6)
    edu_state = float((st.C24010_014E + st.C24010_050E).sum())
    return dict(national_unified_plus_elementary=nat_ue, national_all_layers=nat_all, edu_state_sum=edu_state,
                emp_state_sum=float(st.C24010_001E.sum()),
                top_states=[dict(state=str(r.NAME), all_layers=float(r["all"]), unified_plus_elementary=float(r.ue))
                            for _, r in top.iterrows()])


def add_c24010(A: pd.DataFrame, c: pd.DataFrame) -> pd.DataFrame:
    """ACS columns onto the audited panel (after the NYC collapse: the collapsed row carries the
    Chancellor's id 3620580, which is the id the ACS files New York City under)."""
    n0 = len(A)
    X = A.merge(c, on="leaid", how="left", validate="many_to_one")
    if len(X) != n0:
        raise SystemExit("merge of C24010 changed the number of districts")
    X["a_edu"] = X.edu_res.where(X.emp_res.notna())
    ok = X.edu_res.notna() & X.emp_res.notna() & X.teachers_k12.notna() & X.pupils_k12.notna() \
        & (X.pupils_k12 > 0) & (X.teachers_k12 > 0)
    X["a_fte"] = X.teachers_k12.where(X.pupils_k12.notna() & (X.pupils_k12 > 0))
    X["ok_match"] = ok
    X["m_edu"], X["m_emp"] = X.edu_res.where(ok), X.emp_res.where(ok)
    X["m_fte"], X["m_pup"] = X.teachers_k12.where(ok), X.pupils_k12.where(ok)
    X["m_vedu"], X["m_vemp"] = X.edu_var.where(ok), X.emp_var.where(ok)
    X.attrs.update(A.attrs)
    return X


def czone_table(part: pd.DataFrame, gaps: list["de.Gap"]) -> pd.DataFrame:
    """`czone_gaps` for `gaps` on an already partitioned panel (the swap `cz_gap_table` does)."""
    saved = de.GAPS
    try:
        de.GAPS = list(gaps)
        g = de.czone_gaps(part)
    finally:
        de.GAPS = saved
    return g[g.year == YEAR].reset_index(drop=True)


def variant_panel(part: pd.DataFrame, kind: str) -> pd.DataFrame:
    """The partitioned districts with the matched columns re-defined for a variant:
        earnings_rule  every LEA of the audited panel takes the record of its own layer (the rule of the
                       B20005 earnings moments; right for medians and ratios, but residents of an
                       elementary district inside a secondary district are counted twice when both
                       sit in one location)
        no_overlap     HEADLINE: residents (EDU, employed) from unified and elementary records only, so
                       that nobody is counted twice (the unified and elementary layers tile the states);
                       FTE and pupils from every LEA (the secondary LEAs keep their FTE and pupils)
        no_overlap_sec residents from unified and secondary records instead: the elementary LEAs'
                       residents are dropped in the CZ locations that hold a matched secondary LEA
                       (FTE and pupils from every LEA)
        unified        unified districts only: elementary and secondary LEAs leave every sum
        nyc_dropped    New York City dropped from the resident AND the FTE / pupil counts (the earnings
                       moments drop NYC from the ACS statistics)
    The columns changed: MATCH (the matched sample), VARS (ACS sampling variances) and `a_edu` (EDU
    residents on every district with a record, for the `unmatched` variant)."""
    p = part.copy()
    cols = list(MATCH) + list(VARS)
    if kind == "earnings_rule":
        return p
    if kind == "unified":
        p.loc[p.sd_layer != "unified", cols] = np.nan
    elif kind == "no_overlap":
        sec = p.sd_layer == "secondary"
        have = p.teachers_k12.notna() & p.pupils_k12.notna() & (p.pupils_k12 > 0) & (p.teachers_k12 > 0)
        p.loc[sec & have, "m_fte"] = p.loc[sec & have, "teachers_k12"]
        p.loc[sec & have, "m_pup"] = p.loc[sec & have, "pupils_k12"]
        p.loc[sec & have, ["m_edu", "m_emp", "m_vedu", "m_vemp"]] = 0.0
        p.loc[sec & ~have, cols] = np.nan
        p.loc[sec & p.a_edu.notna(), "a_edu"] = 0.0
    elif kind == "no_overlap_sec":
        sec = (p.sd_layer == "secondary") & p.ok_match
        cells = set(map(tuple, p.loc[sec, ["czone", "loc_id"]].drop_duplicates().to_numpy()))
        in_cell = np.array([(cz, lo) in cells for cz, lo in zip(p.czone, p.loc_id)], dtype=bool)
        el = (p.sd_layer == "elementary") & in_cell
        p.loc[el & p.ok_match, ["m_edu", "m_emp", "m_vedu", "m_vemp"]] = 0.0
        p.loc[el & p.a_edu.notna(), "a_edu"] = 0.0
    elif kind == "nyc_dropped":
        p.loc[p.leaid == ds.NYC, cols] = np.nan
    else:
        raise ValueError(kind)
    return p


def boot_stat(g: pd.DataFrame, col: str, w: np.ndarray, counts: np.ndarray) -> dict[str, Any]:
    s = db.boot_wmean(g[col].to_numpy(float), w, counts)
    out = ds.summarise(s)
    out["draws"] = s["draws"]
    return out


def level_stat(g: pd.DataFrame, key: str, gapcol: str, w: np.ndarray, counts: np.ndarray
               ) -> dict[str, Any]:
    """CZ-weighted mean level of location 1 (city) and 2 (suburb) on the CZs that have the gap."""
    out: dict[str, Any] = {}
    have = g[gapcol].notna().to_numpy()
    for loc in (1, 2):
        x = g[f"{key}_{loc}"].where(have).to_numpy(float)
        out[f"loc{loc}"] = boot_stat(g.assign(_x=x), "_x", w, counts)
    return out


def acs_sampling_se(g: pd.DataFrame, w: np.ndarray) -> dict[str, float]:
    """Sampling SE of the pupil-weighted CZ means of the EDU log gaps from the C24010 margins of error:
    each district's EDU (employed) count has variance MOE^2 / 1.645^2 (men's and women's EDU counts added),
    a location's count the sum over its districts (the Census root-sum-of-squares rule: independence across
    districts), the log gap of a CZ the delta method (FTE are exact counts), and the CZ means independent
    across CZs.  EDU residents are inside employed residents, so cov(EDU, employed) is approximated by the
    variance of EDU in the ratio.  Approximate; the whole-CZ bootstrap already contains most of this noise,
    so it is reported next to the bootstrap SE rather than inside it."""
    col = lambda c: g[c].to_numpy(float)   # noqa: E731
    E1, E2, V1, V2 = col("tr_edu_res_1"), col("tr_edu_res_2"), col("tr_vedu_1"), col("tr_vedu_2")
    N1, N2, W1, W2 = col("tr_emp_res_1"), col("tr_emp_res_2"), col("tr_vemp_1"), col("tr_vemp_2")
    with np.errstate(invalid="ignore", divide="ignore"):
        var_edu = V1 / E1 ** 2 + V2 / E2 ** 2
        one = lambda V, E, Wv, N: np.maximum(V / E ** 2 + Wv / N ** 2 - 2.0 * V / (E * N), 0.0)   # noqa: E731
        var_ratio = one(V1, E1, W1, N1) + one(V2, E2, W2, N2)

    def se_of(var: np.ndarray, gapcol: str) -> float:
        m = np.isfinite(var) & g[gapcol].notna().to_numpy() & (w > 0)
        wn = w[m] / w[m].sum()
        return float(np.sqrt(np.sum(wn ** 2 * var[m])))

    return dict(gap_edu_res=se_of(var_edu, "gap_tr_edu_res"), net_commute_gap_edu=se_of(var_edu, "gap_tr_edu_fte"),
                gap_edu_res_emp=se_of(var_ratio, "gap_tr_teach_emp"))


def district_block(B: int, seed: int) -> dict[str, Any]:
    """Sections (a)-(c): the audited panel, C24010 merged on leaid, the locale CZ table and the
    whole-CZ bootstrap on the draws of `district_settings.py`.  The headline counts residents once
    (`variant_panel(..., 'no_overlap')`); the earnings rule is the variant `earnings_rule`."""
    inp = ds.build_audited_panel(HERE)
    A = add_c24010(inp["A"], load_c24010())
    # --- the audited locale partition and its standard moments (reproduction check) -----------
    gaps = ds.settings_gaps() + GAPS_TR
    part, g_er = db.cz_gap_table(A, gaps, "locale", ds.MIN_PUPILS, ds.MIN_GROUP_SHARE)
    g_er, _ = db.pooled_median_gaps(part, g_er)
    for k in ds.USD_PARTS:
        g_er[f"dlev_{k}"] = g_er[f"{k}_2"] - g_er[f"{k}_1"]
    std = ds.moment_set(g_er, part, B, seed)
    full_run = B == db.B_DEFAULT and seed == db.SEED
    repro = []
    for k in ds.TARGET_KEYS:
        s = std["stats"][k]
        row = dict(moment=k, value=s["value"], se=s["se"], n_cz=s["n_cz"],
                   reported_value=ds.REPORTED[k][0], reported_se=ds.REPORTED[k][1],
                   ok_value=round(s["value"], 4) == ds.REPORTED[k][0],
                   ok_se=(round(s["se"], 4) == ds.REPORTED[k][1]) if full_run else None)
        repro.append(row)
        LOG.info("reproduction %-16s %.4f (%.4f) reported %.4f (%.4f) %s", k, s["value"], s["se"],
                 row["reported_value"], row["reported_se"],
                 "ok" if row["ok_value"] and row["ok_se"] in (True, None) else "MISMATCH")
    if not all(r["ok_value"] and r["ok_se"] in (True, None) for r in repro):
        raise SystemExit("the four locale targets do not reproduce the reported values; stopping")
    counts = db.draw_counts(len(g_er), B, seed)
    w = g_er.pupils.to_numpy(float)
    if not np.allclose(counts, db.draw_counts(len(g_er), B, seed)):
        raise SystemExit("bootstrap draws are not deterministic")

    # --- the headline CZ table: residents counted once ------------------------------------------------
    part_h = variant_panel(part, "no_overlap")
    g_h = czone_table(part_h, GAPS_TR)
    if not g_h.czone.equals(g_er.czone):
        raise SystemExit("the headline CZ table differs from the audited one")
    g = g_er.copy()
    trcols = [c for c in g_h.columns if c.startswith(("gap_tr_", "tr_"))]
    g[trcols] = g_h[trcols].to_numpy()

    # --- the statistics, headline specification --------------------------------------------------------
    def stats_for(gt: pd.DataFrame) -> dict[str, dict[str, Any]]:
        out = {}
        for key, col, kind, label in STATS_TR:
            if col in gt and gt[col].notna().sum() >= 5:
                out[key] = dict(boot_stat(gt, col, w, counts), kind=kind, label=label)
        # identity: net_commute_gap_edu = gap_edu_res - gap_fte, CZ by CZ
        if {"gap_edu_res", "gap_fte", "net_commute_gap_edu"} <= set(out):
            d = (gt.gap_tr_edu_res - gt.gap_tr_fte - gt.gap_tr_edu_fte).abs().max()
            out["net_commute_gap_edu"]["identity_max_abs_diff"] = float(d)
        return out

    main_stats = stats_for(g)
    levels = {key: level_stat(g, gk, f"gap_{gk}", w, counts) for key, gk, _ in LEVELS_TR}
    D_std = std["draws"][CORR_TARGETS]
    D_new = pd.DataFrame({k: v["draws"] for k, v in main_stats.items()})
    D = pd.concat([D_std.reset_index(drop=True), D_new.reset_index(drop=True)], axis=1)
    corr = D.corr()
    draws_main = {k: v["draws"] for k, v in main_stats.items()}
    main_stats = {k: {kk: vv for kk, vv in v.items() if kk != "draws"} for k, v in main_stats.items()}
    for lv in levels.values():
        for s in lv.values():
            s.pop("draws", None)
    acs_se = acs_sampling_se(g, w)
    for k, v in acs_se.items():
        main_stats[k]["se_acs_sampling"] = v
        main_stats[k]["se_total"] = float(np.hypot(main_stats[k]["se"], v))
        LOG.info("ACS sampling SE of %s: %.4f (bootstrap SE %.4f)", k, v, main_stats[k]["se"])

    # --- variants: overlap rule, layers, NYC, unmatched samples ----------------------------------------
    variants: dict[str, Any] = {}
    panels = {"earnings_rule": variant_panel(part, "earnings_rule"), "no_overlap_sec": variant_panel(part, "no_overlap_sec"),
              "unified": variant_panel(part, "unified"), "nyc_dropped": variant_panel(part_h, "nyc_dropped")}
    for kind, pv in panels.items():
        gv = g_er if kind == "earnings_rule" else czone_table(pv, GAPS_TR)
        if not gv.czone.equals(g.czone):
            raise SystemExit("variant CZ table differs from the main one")
        variants[kind] = stats_for(gv)
        variants[kind]["_n_matched_districts"] = int(pv.m_edu.notna().sum())
        variants[kind]["_gap_cols"] = {c: gv[c].to_numpy(float) for c in ("gap_tr_edu_fte", "gap_tr_teach_emp", "gap_tr_edu_res")}
        LOG.info("variant %-14s net_commute_gap_edu %.4f (%.4f), gap_edu_res %.4f, gap_fte %.4f, gap_edu_res_emp %.4f", kind,
                 variants[kind]["net_commute_gap_edu"]["value"], variants[kind]["net_commute_gap_edu"]["se"],
                 variants[kind]["gap_edu_res"]["value"], variants[kind]["gap_fte"]["value"],
                 variants[kind]["gap_edu_res_emp"]["value"])
    # leave the New York CZ out (NYC is one district of 0.96 m pupils with the Chancellor's record)
    ny_cz = float(part.czone[part.leaid == ds.NYC].iloc[0]) if (part.leaid == ds.NYC).any() else None
    if ny_cz is not None:
        gn = g.copy()
        gcols = [c for c in gn.columns if c.startswith("gap_tr_")]
        gn.loc[gn.czone == ny_cz, gcols] = np.nan
        variants["drop_ny_cz"] = stats_for(gn)
        variants["drop_ny_cz"]["_gap_cols"] = {c: gn[c].to_numpy(float) for c in ("gap_tr_edu_fte", "gap_tr_teach_emp", "gap_tr_edu_res")}
    # each statistic on its own maximal sample (residents counted once)
    gu = g.copy()
    gu["gap_net_unmatched"] = gu.gap_tr_edu_all - gu.gap_tr_fte_all
    variants["unmatched"] = {
        "gap_edu_res": dict(boot_stat(gu, "gap_tr_edu_all", w, counts),
                            label="Delta log EDU residents, every district with a C24010 record"),
        "gap_fte": dict(boot_stat(gu, "gap_tr_fte_all", w, counts),
                        label="Delta log K-12 FTE, every district with an audited FTE"),
        "net_commute_gap_edu": dict(boot_stat(gu, "gap_net_unmatched", w, counts),
                                    label="gap_edu_res - gap_fte on the two maximal samples")}
    for v in variants.values():
        for s in v.values():
            if isinstance(s, dict):
                s.pop("draws", None)

    # --- coverage and the elementary / secondary overlap ---------------------------------------
    cover = coverage_tr(part, g_er, A, part_h)
    cover["tiling"] = tiling_check()
    cover["ccd_fte_k12_all_regular_districts"] = float(A.teachers_k12.sum())
    cover["ccd_pupils_k12_all_regular_districts"] = float(A.pupils_k12.sum())
    return dict(A=A, part=part, g=g, std=std, counts=counts, w=w, counts_shape=list(counts.shape), reproduction=repro,
                stats=main_stats, draws_main=draws_main, levels=levels, corr=corr, draws=D, variants=variants,
                coverage=cover, hybrid_share=inp["hybrid_share"], full_run=full_run,
                geo_unit=A.attrs.get("geo_unit"), acs_se=acs_se)


def coverage_tr(part: pd.DataFrame, g: pd.DataFrame, A: pd.DataFrame, part_h: pd.DataFrame | None = None
                ) -> dict[str, Any]:
    """Size of the matched sample and the elementary / secondary overlap of the ACS records (`part` carries
    the earnings-rule matched columns, `part_h` the headline's)."""
    pu = float(part.pupils.sum())
    layer = part.sd_layer.fillna("none")
    ok = part.ok_match
    out: dict[str, Any] = dict(
        n_cz=int(len(g)), n_districts=int(len(part)), pupils=pu,
        n_matched=int(ok.sum()), share_districts_matched=float(ok.mean()),
        share_pupils_matched=float(part.pupils[ok].sum() / pu),
        share_pupils_with_acs=float(part.pupils[part.edu_res.notna()].sum() / pu),
        share_pupils_with_fte=float(part.pupils[part.teachers_k12.notna() & part.pupils_k12.notna()].sum() / pu),
        n_by_layer={k: int(v) for k, v in layer.value_counts().items()},
        pupil_share_by_layer={k: float(part.pupils[layer == k].sum() / pu) for k in layer.unique()},
        employed_res_matched=float(part.m_emp.sum()), edu_res_matched=float(part.m_edu.sum()),
        fte_matched=float(part.m_fte.sum()), pupils_k12_matched=float(part.m_pup.sum()))
    mm = part[ok]
    tot = float(mm.m_emp.sum())
    out["employed_res_share_by_layer"] = {k: float(mm.m_emp[mm.sd_layer == k].sum() / tot)
                                          for k in ("unified", "elementary", "secondary")}
    ue = float(mm.m_emp[mm.sd_layer.isin(["unified", "elementary"])].sum())
    out["overlap_factor_employed_residents"] = float(tot / ue)       # >1: residents counted twice
    byl = mm.groupby(["loc_id", "sd_layer"]).m_emp.sum().unstack("sd_layer").fillna(0.0)
    out["overlap_factor_by_location"] = {int(lo): float(r.sum() / (r.get("unified", 0.0) + r.get("elementary", 0.0)))
                                         for lo, r in byl.iterrows()}
    ee = mm.groupby(["czone", "loc_id", "sd_layer"]).m_emp.sum().unstack("sd_layer").fillna(0.0)
    ee = ee[ee.get("secondary", 0) > 0]                              # CZ locations that hold a secondary LEA
    out["secondary_vs_elementary_employed_residents"] = dict(
        n_cz_locations=int(len(ee)), elementary=float(ee.elementary.sum()), secondary=float(ee.secondary.sum()),
        by_location={int(lo): dict(elementary=float(v.elementary.sum()), secondary=float(v.secondary.sum()))
                     for lo, v in ee.groupby(level="loc_id")})
    sec_cz = set(mm.czone[mm.sd_layer == "secondary"])
    wcz = g.set_index("czone").pupils
    out["n_cz_with_secondary_lea"] = int(len(sec_cz))
    out["pupil_weight_cz_with_secondary_lea"] = float(wcz[wcz.index.isin(sec_cz)].sum() / wcz.sum())
    out["n_by_location"] = {int(k): int(v) for k, v in part.groupby("loc_id").size().items()}
    out["n_matched_by_location"] = {int(k): int(v) for k, v in mm.groupby("loc_id").size().items()}
    ny = part[part.leaid == ds.NYC]
    if len(ny):
        out["nyc"] = dict(czone=float(ny.czone.iloc[0]), loc_id=int(ny.loc_id.iloc[0]),
                          pupils=float(ny.pupils.iloc[0]), matched=bool(ny.ok_match.iloc[0]),
                          edu_res=float(ny.edu_res.iloc[0]) if ny.edu_res.notna().iloc[0] else None,
                          emp_res=float(ny.emp_res.iloc[0]) if ny.emp_res.notna().iloc[0] else None,
                          fte=float(ny.teachers_k12.iloc[0]) if ny.teachers_k12.notna().iloc[0] else None)
    nm = part[~ok]
    out["unmatched_reasons"] = dict(
        no_acs_record=int(nm.edu_res.isna().sum()),
        no_fte_or_pupils=int((nm.teachers_k12.isna() | nm.pupils_k12.isna()).sum()),
        pupils_of_unmatched=float(nm.pupils.sum()))
    w = g.pupils.to_numpy(float)
    wn = w / w.sum()
    out["n_eff_cz"] = float(1.0 / (wn ** 2).sum())
    out["top5_weight_share"] = float(np.sort(wn)[::-1][:5].sum())
    if part_h is not None:                                           # the headline (residents counted once)
        okh = part_h.m_fte.notna() & part_h.m_pup.notna()
        sec0 = (part_h.sd_layer == "secondary") & part_h.m_edu.eq(0)
        out["single_count"] = dict(
            n_matched=int(okh.sum()), share_pupils_matched=float(part_h.pupils[okh].sum() / pu),
            n_secondary_without_residents=int(sec0.sum()),
            employed_res=float(part_h.m_emp.sum()), edu_res=float(part_h.m_edu.sum()),
            fte=float(part_h.m_fte.sum()), pupils_k12=float(part_h.m_pup.sum()),
            n_matched_by_location={int(k): int(v) for k, v in part_h[okh].groupby("loc_id").size().items()})
    return out


# --------------------------------------------------------------------------
# PUMA geography: principal-city / suburban classification and the PUMA -> POWPUMA conversion
# --------------------------------------------------------------------------

PUMSEQ_COLS = [(0, 3), (3, 5), (13, 18), (18, 21), (29, 34), (42, 47), (55, 61), (61, 70),
               (70, 79), (79, 179)]
PUMSEQ_NAMES = ["sl", "st", "puma", "county", "cousub", "place", "tract", "pop10", "hu10", "name"]


def read_pumseq(st: int) -> pd.DataFrame:
    """One state's 2010 PUMA equivalency file (PUMSEQ10_<ST>.txt, fixed width, layout in
    2010_PUMA_Equivalency_Format_Layout.pdf): summary levels 795 (PUMA), 796 (PUMA-county),
    797 (county subdivision, 20 MCD states), 798 (place or remainder) and 799 (tract), each with
    the 2010 population."""
    path = GEO_DIR / f"PUMSEQ10_{st:02d}.txt"
    if not path.exists():
        raise SystemExit(f"no {path.name}: run `teacher_residence.py fetch --only geo`")
    d = pd.read_fwf(path, colspecs=PUMSEQ_COLS, names=PUMSEQ_NAMES, header=None, dtype=str,
                    encoding="latin-1", keep_default_na=False)
    d["pop10"] = pd.to_numeric(d.pop10, errors="coerce")
    return d


def read_omb() -> tuple[pd.DataFrame, pd.DataFrame]:
    """OMB September 2018 delineation: principal cities of metropolitan statistical areas
    (list 2: state and place FIPS) and the counties of metropolitan statistical areas (list 1)."""
    l2 = pd.read_excel(GEO_DIR / OMB_FILES["list2"], header=2, dtype=str)
    l2 = l2[l2["FIPS State Code"].notna() & l2["FIPS Place Code"].notna()]
    l2 = l2[l2["Metropolitan/Micropolitan Statistical Area"].str.startswith("Metropolitan")]
    l2 = l2[l2["FIPS State Code"].str.zfill(2) != "72"]                   # Puerto Rico is not in the PUMS used
    pcs = pd.DataFrame({"st": l2["FIPS State Code"].str.zfill(2).to_numpy(),
                        "place": l2["FIPS Place Code"].str.zfill(5).to_numpy(),
                        "pc_name": l2["Principal City Name"].to_numpy(),
                        "cbsa": l2["CBSA Title"].to_numpy()})
    l1 = pd.read_excel(GEO_DIR / OMB_FILES["list1"], header=2, dtype=str)
    l1 = l1[l1["FIPS State Code"].notna() & l1["FIPS County Code"].notna()]
    l1 = l1[l1["Metropolitan/Micropolitan Statistical Area"].str.startswith("Metropolitan")]
    cty = pd.DataFrame({"st": l1["FIPS State Code"].str.zfill(2).to_numpy(),
                        "county": l1["FIPS County Code"].str.zfill(3).to_numpy(),
                        "status": l1["Central/Outlying County"].to_numpy()})
    return pcs.drop_duplicates(["st", "place"]), cty.drop_duplicates(["st", "county"])


def read_puma_powpuma() -> pd.DataFrame:
    """Census PUMA -> POWPUMA conversion of the 2018 PUMS code lists (sheet
    `PUMA_POWPUMA_MIGPUMA`): for each residence PUMA the POWSP and POWPUMA code of the place-of-work
    area that contains it (POWPUMAs are county based and aggregate one or more PUMAs)."""
    x = pd.read_excel(GEO_DIR / "ACSPUMS2018CodeLists.xls", sheet_name="PUMA_POWPUMA_MIGPUMA",
                      header=3, dtype=str)
    x.columns = ["st", "puma", "powsp", "powpuma"]
    x = x.dropna(subset=["st", "puma", "powsp", "powpuma"])
    x["powpuma"] = x.powpuma.str.replace("*", "", regex=False)       # footnote flags (**, ***)
    x = x[x.st.str.fullmatch(r"\d+") & x.powpuma.str.fullmatch(r"\d+")]
    for c in x.columns:
        x[c] = x[c].astype(int)
    return x.drop_duplicates(["st", "puma"]).reset_index(drop=True)


def set_classes(geo: dict[str, Any], pc_cut: float = 0.5, metro_cut: float = 0.5) -> dict[str, Any]:
    """Assign `cls` to the PUMAs and POWPUMAs of `geo` (in place) from their principal-city and
    metropolitan shares: 'city' if pc_share >= pc_cut; 'suburb' if metro_share >= metro_cut and
    pc_share < pc_cut; else 'other'."""
    G, P = geo["puma"], geo["powpuma"]
    for df in (G, P):
        df["cls"] = np.where(df.pc_share >= pc_cut, "city",
                             np.where(df.metro_share >= metro_cut, "suburb", "other"))
    metro = G.cls != "other"
    geo["info"].update(
        pc_cut=pc_cut, metro_cut=metro_cut,
        pop_by_class={k: float(v) for k, v in G.groupby("cls")["pop"].sum().items()},
        n_puma_by_class={k: int(v) for k, v in G.groupby("cls").size().items()},
        n_powpuma_by_class={k: int(v) for k, v in P.groupby("cls").size().items()},
        pc_share_of_metro_pop=float(G.pop_pc[metro].sum() / G["pop"][metro].sum()))
    return geo


def build_puma_geo(pc_cut: float = 0.5, metro_cut: float = 0.5) -> dict[str, Any]:
    """Principal-city / suburban / other classification of the 2010 PUMAs and of the POWPUMAs.

    For each PUMA: pc_share = 2010 population living in places that are principal cities of a
    metropolitan statistical area (OMB September 2018, matched on state and place FIPS) over
    the PUMA's 2010 population; metro_share = 2010 population in counties of metropolitan
    statistical areas (`set_classes` turns them into classes).  A POWPUMA takes the
    population-weighted shares of the PUMAs it aggregates."""
    eq = pd.concat([read_pumseq(st) for st in ao.STATES], ignore_index=True)
    pcs, mcty = read_omb()
    to_int = lambda df: df.assign(st=df.st.astype(int), puma=df.puma.astype(int))     # noqa: E731
    pup = to_int(eq[eq.sl == "795"]).groupby(["st", "puma"]).pop10.sum().rename("pop").reset_index()
    pl = eq[eq.sl == "798"].merge(pcs, on=["st", "place"], how="inner")
    pc = to_int(pl).groupby(["st", "puma"]).pop10.sum().rename("pop_pc").reset_index()
    ct = eq[eq.sl == "796"].merge(mcty, on=["st", "county"], how="inner")
    me = to_int(ct).groupby(["st", "puma"]).pop10.sum().rename("pop_metro").reset_index()
    G = pup.merge(pc, on=["st", "puma"], how="left").merge(me, on=["st", "puma"], how="left")
    G[["pop_pc", "pop_metro"]] = G[["pop_pc", "pop_metro"]].fillna(0.0)
    names = to_int(eq[eq.sl == "795"])
    G = G.merge(names[["st", "puma", "name"]].drop_duplicates(["st", "puma"]), on=["st", "puma"], how="left")
    G = G.merge(read_puma_powpuma(), on=["st", "puma"], how="left")
    # the commuting zone (Dorn 1990, as in the district panel) that holds most of each PUMA's 2010 population
    cz = de.build_czones(de.RawStore(RAW))
    if cz is not None:
        c796 = eq[eq.sl == "796"].copy()
        c796["county_fips"] = c796.st.str.zfill(2) + c796.county.str.zfill(3)
        c796 = to_int(c796.merge(cz, on="county_fips", how="left"))
        tot = c796.groupby(["st", "puma"]).pop10.sum().rename("pop_all").reset_index()
        best = (c796.dropna(subset=["czone"]).groupby(["st", "puma", "czone"]).pop10.sum().reset_index()
                .sort_values("pop10", ascending=False).drop_duplicates(["st", "puma"]).merge(tot, on=["st", "puma"]))
        best["cz_share"] = best.pop10 / best.pop_all
        G = G.merge(best[["st", "puma", "czone", "cz_share"]], on=["st", "puma"], how="left")
    else:
        G["czone"], G["cz_share"] = np.nan, np.nan
    P = G.groupby(["powsp", "powpuma"])[["pop", "pop_pc", "pop_metro"]].sum().reset_index()
    for df in (G, P):
        df["pc_share"] = df.pop_pc / df["pop"]
        df["metro_share"] = df.pop_metro / df["pop"]
    found = pl[["st", "place"]].drop_duplicates()
    unmatched = pcs.merge(found, on=["st", "place"], how="left", indicator=True)
    unmatched = unmatched[unmatched._merge == "left_only"]
    info = dict(n_puma=int(len(G)), n_powpuma=int(len(P)), n_pc_omb=int(len(pcs)), n_pc_found=int(len(found)),
                n_pc_unmatched=int(len(unmatched)),
                unmatched_examples=[f"{r.pc_name} (state {r.st}, place {r.place}; {r.cbsa})" for r in unmatched.head(12).itertuples()],
                pop_total=float(G["pop"].sum()), n_puma_without_powpuma=int(G.powpuma.isna().sum()))
    return set_classes(dict(puma=G, powpuma=P, info=info), pc_cut, metro_cut)


# --------------------------------------------------------------------------
# PUMS (d)-(e): replicate-weight engine
# --------------------------------------------------------------------------

class Rep:
    """ACS PUMS person weight and its 80 successive-difference replicate weights.  Census
    variance estimator (ACS PUMS accuracy documentation): SE(x) = sqrt(4/80 sum_r (x_r - x)^2),
    x_r the statistic computed with replicate weight r.  `W` is (n, 80), Fortran order so that a
    replicate column is contiguous."""

    def __init__(self, w0: np.ndarray, W: np.ndarray) -> None:
        self.w0 = np.asarray(w0, float)
        self.W = np.asfortranarray(W)

    def sub(self, mask: np.ndarray) -> "Rep":
        return Rep(self.w0[mask], self.W[mask])

    def run(self, fn: Callable[[np.ndarray], Any]) -> tuple[np.ndarray, np.ndarray]:
        """fn(weights) -> scalar or array; returns (estimate, replicates (80, ...))."""
        est = np.asarray(fn(self.w0), float)
        reps = np.empty((N_REP,) + est.shape)
        for r in range(N_REP):
            reps[r] = fn(self.W[:, r].astype(float))
        return est, reps


def sdr_se(est: np.ndarray, reps: np.ndarray) -> np.ndarray:
    return np.sqrt(4.0 / N_REP * ((reps - est) ** 2).sum(axis=0))


class Ratios:
    """Ratio-of-weighted-sums statistics sum(w * num) / sum(w * den), many at once, with the
    replicate-weight standard error; derived statistics (differences, log ratios) are taken on the
    replicates, so their standard errors respect the covariance."""

    def __init__(self, rep: Rep, specs: dict[str, tuple[np.ndarray, np.ndarray]]) -> None:
        self.names = list(specs)
        self.ix = {k: i for i, k in enumerate(self.names)}
        NUM = np.column_stack([np.asarray(specs[k][0], float) for k in self.names])
        DEN = np.column_stack([np.asarray(specs[k][1], float) for k in self.names])

        def fn(w: np.ndarray) -> np.ndarray:
            with np.errstate(invalid="ignore", divide="ignore"):
                return (w @ NUM) / (w @ DEN)

        self.est, self.reps = rep.run(fn)

    def get(self, k: str) -> dict[str, float]:
        i = self.ix[k]
        return dict(value=float(self.est[i]), se=float(sdr_se(self.est[i], self.reps[:, i])))

    def derive(self, f: Callable[[Callable[[str], Any]], Any]) -> dict[str, float]:
        """f(get) with get(name) the statistic; evaluated at the estimate and at every replicate."""
        e = f(lambda k: self.est[self.ix[k]])
        r = np.array([f(lambda k, j=j: self.reps[j, self.ix[k]]) for j in range(N_REP)])
        return dict(value=float(e), se=float(sdr_se(np.asarray(e), r)))


def wquantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    ok = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x, w = x[ok], w[ok]
    if len(x) == 0:
        return float("nan")
    o = np.argsort(x, kind="stable")
    x, w = x[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    return float(np.interp(q * w.sum(), c, x))


def fe_coefs(y: np.ndarray, X: np.ndarray, g: np.ndarray, w: np.ndarray, G: int) -> np.ndarray:
    """Weighted least squares of y on X with group fixed effects (within transformation); the
    coefficients on the columns of X."""
    sw = np.bincount(g, weights=w, minlength=G)
    sw = np.where(sw > 0, sw, 1.0)

    def dm(v: np.ndarray) -> np.ndarray:
        return v - (np.bincount(g, weights=w * v, minlength=G) / sw)[g]

    Xd = np.column_stack([dm(X[:, j]) for j in range(X.shape[1])])
    yd = dm(y)
    A = (Xd * w[:, None]).T @ Xd
    b = (Xd * w[:, None]).T @ yd
    return np.linalg.solve(A, b)


# --------------------------------------------------------------------------
# PUMS: person file, geography merged
# --------------------------------------------------------------------------

def attach_classes(d: pd.DataFrame, geo: dict[str, Any]) -> None:
    """Class of the residence PUMA (`res_cls`), of the POWPUMA that contains it (`res_pow_cls`) and of the
    place-of-work POWPUMA (`wk_cls`) of every person (in place; needs `pid`)."""
    G, P = geo["puma"], geo["powpuma"]
    gi = G.assign(pid=G.st * 100000 + G.puma).set_index("pid")
    pk = P.assign(key=P.powsp * 100000 + P.powpuma).set_index("key")
    for c in ("cls", "pc_share", "powsp", "powpuma"):
        d[f"res_{c}"] = d.pid.map(gi[c])
    d["res_pow_key"] = d.res_powsp * 100000 + d.res_powpuma
    d["res_pow_cls"] = d.res_pow_key.map(pk.cls)
    d["wk_key"] = d.POWSP * 100000 + d.POWPUMA
    d["wk_cls"] = d.wk_key.map(pk.cls)


def load_pums(geo: dict[str, Any]) -> tuple[pd.DataFrame, Rep]:
    """Employed (civilian) persons 16+ of the ACS 2018 1-year PUMS in households (group quarters
    dropped), with the PUMA classification of residence and of the place-of-work area."""
    files = sorted(PUMS_DIR.glob("*.parquet"))
    if len(files) < len(ao.STATES):
        raise SystemExit(f"only {len(files)} PUMS state files in {PUMS_DIR}; run `fetch --only pums`")
    d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    n_all = len(d)
    # pipeline check against the published national means (ACS 2018 1-year: 27.1 minutes of travel time to work
    # for workers 16+ who do not work at home, 5.3% working at home)
    mt = ((d.JWMNP > 0) & (d.JWTR != 11)).to_numpy()
    check = dict(mean_travel_time_all_persons=float(np.average(d.JWMNP.to_numpy()[mt], weights=d.PWGTP.to_numpy()[mt])),
                 work_at_home_share_all_persons=float(d.PWGTP[d.JWTR == 11].sum() / d.PWGTP.sum()),
                 published_mean_travel_time=27.1, published_work_at_home_share=0.053)
    d = d[~d.SERIALNO.str.contains("GQ")].reset_index(drop=True)
    LOG.info("PUMS 2018: %d employed persons 16+ in %d files, %d in households (%.2f m weighted)",
             n_all, len(files), len(d), d.PWGTP.sum() / 1e6)
    d["st"] = d.ST.astype(int)
    d["pid"] = d.st * 100000 + d.PUMA
    attach_classes(d, geo)
    occ = d.OCCP.to_numpy()
    d["edu"] = np.isin(occ, list(EDU_CODES))
    d["k12"] = np.isin(occ, K12_CODES)
    d["k12p"] = np.isin(occ, K12_PLUS)
    d["pub"] = d.COW.isin([3, 4])                       # local or state government employee
    d["female"] = d.SEX.eq(2)
    d["home"] = d.JWTR.eq(11)
    d["comm"] = d.JWMNP.gt(0) & ~d.home                 # commutes to a place of work away from home
    d["pow_ok"] = d.comm & d.POWSP.between(1, 56) & d.POWPUMA.gt(1) & d.res_powpuma.notna()
    d["out_pow"] = d.pow_ok & ((d.POWSP != d.res_powsp) | (d.POWPUMA != d.res_powpuma))
    rep = Rep(d.PWGTP.to_numpy(float), d[REP].to_numpy(np.int32))
    d = d.drop(columns=REP)
    d.attrs["national_check"] = check
    return d, rep


def stat_row(r: dict[str, float], d: int = 4) -> dict[str, float]:
    return {k: (round(v, d + 4) if isinstance(v, float) else v) for k, v in r.items()}


# --------------------------------------------------------------------------
# (d) composition of SOC 25-0000
# --------------------------------------------------------------------------

CLASSES = ("city", "suburb", "other")


def share_specs(prefix: str, sub: pd.DataFrame, mask: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Shares within the EDU workers selected by `mask`: each EDU occupation, K-12 teachers
    (elementary/middle, secondary, special education), with preschool/kindergarten added, and the
    public-sector (local or state government employer) K-12 teachers."""
    m = mask.astype(float)
    occ = sub.OCCP.to_numpy()
    k12, k12p, pub = sub.k12.to_numpy(), sub.k12p.to_numpy(), sub.pub.to_numpy()
    sp = {f"{prefix}_occ_{c}": ((occ == c) * m, m) for c in EDU_CODES}
    sp[f"{prefix}_k12"] = (k12 * m, m)
    sp[f"{prefix}_k12p"] = (k12p * m, m)
    sp[f"{prefix}_pub_k12"] = ((k12 & pub) * m, m)
    sp[f"{prefix}_pub_k12p"] = ((k12p & pub) * m, m)
    sp[f"{prefix}_pub_k12_sch"] = ((k12 & pub & (sub.INDP.to_numpy() == 7860)) * m, m)
    return sp


def puma_dispersion(sub: pd.DataFrame, rep: Rep, flag: np.ndarray, min_n: int = 1) -> dict[str, Any]:
    """Dispersion across PUMAs of the share of EDU workers who are `flag`: weighted (by the PUMA's
    EDU workers) and unweighted standard deviation, percentiles, the sampling variance of each
    PUMA's share from its replicates and the SD net of it."""
    pid = sub.pid.to_numpy()
    ids, inv = np.unique(pid, return_inverse=True)
    P = len(ids)
    n_samp = np.bincount(inv, minlength=P)
    fl = flag.astype(float)

    def shares(w: np.ndarray) -> np.ndarray:
        E = np.bincount(inv, weights=w, minlength=P)
        K = np.bincount(inv, weights=w * fl, minlength=P)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.concatenate([np.where(E > 0, K / E, np.nan), E])

    est, reps = rep.run(shares)
    s0, E0 = est[:P], est[P:]
    var_p = 4.0 / N_REP * np.nansum((reps[:, :P] - s0) ** 2, axis=0)       # sampling variance of each PUMA share
    keep = (n_samp >= min_n) & np.isfinite(s0) & (E0 > 0)
    vbar = float(np.sum(E0[keep] * var_p[keep]) / np.sum(E0[keep]))
    vbar_u = float(np.mean(var_p[keep]))

    def disp(vec: np.ndarray, vb: float, vbu: float) -> np.ndarray:
        s, E = vec[:P], vec[P:]
        k = keep & np.isfinite(s) & (E > 0)
        mu = np.sum(E[k] * s[k]) / np.sum(E[k])
        v = np.sum(E[k] * (s[k] - mu) ** 2) / np.sum(E[k])
        mu_u = np.mean(s[k])
        v_u = np.mean((s[k] - mu_u) ** 2)
        q = [wquantile(s[k], E[k], p) for p in (0.1, 0.25, 0.5, 0.75, 0.9)]
        return np.array([mu, np.sqrt(v), np.sqrt(max(v - vb, 0.0)), mu_u, np.sqrt(v_u),
                         np.sqrt(max(v_u - vbu, 0.0))] + q)

    d0 = disp(est, vbar, vbar_u)
    dr = np.array([disp(reps[j], vbar, vbar_u) for j in range(N_REP)])
    se = sdr_se(d0, dr)
    names = ["mean_w", "sd_w", "sd_w_net_of_sampling", "mean_u", "sd_u", "sd_u_net_of_sampling",
             "p10", "p25", "p50", "p75", "p90"]
    out = {k: dict(value=float(v), se=float(s)) for k, v, s in zip(names, d0, se)}
    out.update(n_puma=int(keep.sum()), mean_sampling_var_w=vbar, mean_sampling_sd_w=float(np.sqrt(vbar)),
               mean_edu_sample_per_puma=float(n_samp[keep].mean()),
               min_edu_sample_per_puma=int(n_samp[keep].min()))
    return out


def composition_block(d: pd.DataFrame, rep: Rep) -> dict[str, Any]:
    """Section (d): composition of the education, training and library group of the PUMS."""
    e = d.edu.to_numpy()
    sub, rs = d[e].reset_index(drop=True), rep.sub(e)
    n_edu = len(sub)
    cls = sub.res_cls.to_numpy()
    masks = {"all": np.ones(n_edu, bool), "city": cls == "city", "suburb": cls == "suburb",
             "other": cls == "other", "metro": (cls == "city") | (cls == "suburb")}
    specs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for k, m in masks.items():
        specs.update(share_specs(k, sub, m))
    # the EDU workers' own distribution across PUMA classes (weighted shares)
    for k in CLASSES:
        specs[f"edu_in_{k}"] = ((cls == k).astype(float), np.ones(n_edu))
    R = Ratios(rs, specs)
    out: dict[str, Any] = dict(n_sample=int(n_edu), weighted=float(rs.w0.sum()))
    est, _ = rs.run(lambda w: np.array([w.sum(), w @ sub.k12.to_numpy(), w @ (sub.k12 & sub.pub).to_numpy(),
                                        w @ sub.k12p.to_numpy()]))
    out["counts"] = dict(edu=float(est[0]), k12=float(est[1]), pub_k12=float(est[2]), k12p=float(est[3]))
    out["share_by_occ"] = {str(c): dict(label=lab, **R.get(f"all_occ_{c}")) for c, lab in EDU_CODES.items()}
    # the SOCP codes found with each EDU occupation code (one SOC group per occupation code)
    out["socp_by_occ"] = {str(c): sorted(set(sub.SOCP[sub.OCCP == c])) for c in EDU_CODES}
    for k in ("k12", "k12p", "pub_k12", "pub_k12p", "pub_k12_sch"):
        out[f"share_{k}"] = R.get(f"all_{k}")
    out["by_class"] = {}
    for c in ("city", "suburb", "other", "metro"):
        row = {k: R.get(f"{c}_{k}") for k in ("k12", "k12p", "pub_k12", "pub_k12p")}
        row["occ"] = {str(o): R.get(f"{c}_occ_{o}") for o in EDU_CODES}
        if c != "metro":
            row["edu_workers_share"] = R.get(f"edu_in_{c}")
        out["by_class"][c] = row
    # suburb minus city: share differences and log ratios
    out["split"] = {}
    for k in ("k12", "k12p", "pub_k12", "pub_k12p"):
        out["split"][k] = dict(
            diff=R.derive(lambda g, k=k: g(f"suburb_{k}") - g(f"city_{k}")),
            log_ratio=R.derive(lambda g, k=k: np.log(g(f"suburb_{k}")) - np.log(g(f"city_{k}"))))
    out["split"]["occ"] = {str(o): dict(
        diff=R.derive(lambda g, o=o: g(f"suburb_occ_{o}") - g(f"city_occ_{o}")),
        log_ratio=R.derive(lambda g, o=o: np.log(g(f"suburb_occ_{o}")) - np.log(g(f"city_occ_{o}"))))
        for o in EDU_CODES}
    out["split"]["postsecondary_etc_share_city"] = R.get("city_occ_2205")
    # dispersion across PUMAs
    out["dispersion"] = {
        "k12": puma_dispersion(sub, rs, sub.k12.to_numpy()),
        "pub_k12": puma_dispersion(sub, rs, (sub.k12 & sub.pub).to_numpy()),
    }
    out["_ratios"] = R
    return out


def cz_composition(d: pd.DataFrame, rep: Rep, geo: dict[str, Any], cz_w: pd.Series) -> dict[str, Any]:
    """The composition term inside the CZs of the district partition: for each CZ the share of EDU workers
    living in suburban PUMAs who are (public) K-12 teachers over the same share in city PUMAs, log ratio, and
    the pupil-weighted mean over the CZs with both classes (the weights of the district moments).  PUMAs are
    assigned to the CZ that holds most of their 2010 population."""
    G = geo["puma"]
    gi = G.assign(pid=G.st * 100000 + G.puma).set_index("pid")
    e = d.edu.to_numpy()
    sub, rs = d[e].reset_index(drop=True), rep.sub(e)
    n = len(cz_w)
    pos = pd.Series(np.arange(n), index=cz_w.index.to_numpy(float))
    ci = sub.pid.map(gi.czone).map(pos)
    cls = sub.res_cls.to_numpy()
    ok = ci.notna().to_numpy() & ((cls == "city") | (cls == "suburb"))
    code = (ci.to_numpy()[ok].astype(int) * 2 + (cls[ok] == "suburb")).astype(int)
    wz = cz_w.to_numpy(float)
    out: dict[str, Any] = dict(n_edu_sample=int(ok.sum()), n_edu_sample_all=int(len(sub)),
                               n_cz_with_pumas=int(len(np.unique(code // 2))))
    for lab, flag in (("pub_k12", (sub.k12 & sub.pub).to_numpy()), ("k12", sub.k12.to_numpy())):
        fl = flag[ok].astype(float)

        def fn(w: np.ndarray, fl=fl) -> np.ndarray:
            wv = w[ok]
            E = np.bincount(code, weights=wv, minlength=2 * n)
            K = np.bincount(code, weights=wv * fl, minlength=2 * n)
            with np.errstate(invalid="ignore", divide="ignore"):
                sh = np.where(E > 0, K / E, np.nan)
            c, sb = sh[0::2], sh[1::2]
            good = np.isfinite(c) & np.isfinite(sb) & (c > 0) & (sb > 0)
            wt = wz[good]
            return np.array([np.sum(wt * np.log(sb[good] / c[good])) / np.sum(wt), good.sum(), wt.sum() / wz.sum()])

        est, reps = rs.run(fn)
        se = sdr_se(est, reps)
        out[lab] = dict(value=float(est[0]), se=float(se[0]), n_cz=int(est[1]), weight_share=float(est[2]))
    return out


def cz_terms(d: pd.DataFrame, rep: Rep, geo: dict[str, Any], cz_w: pd.Series, num: np.ndarray,
             den: np.ndarray | None = None) -> dict[str, Any]:
    """Composition terms CZ BY CZ: for each CZ of the district partition the log of the share of the `den`
    workers (default: EDU workers) living in suburban PUMAs who are `num`, over the same share in city PUMAs
    (PUMAs go to the CZ that holds most of their 2010 population; classes as attached to `d`).  Returns the
    per-CZ point estimates `est` (NaN where a class is missing or a share is zero), the replicate estimates
    `reps` (N_REP x n_cz) and the pupil-weighted mean over the CZs that have the term (`value`, replicate-weight
    `se_survey`, `n_cz`, `weight_share`); `cz_composition` is the same mean for two fixed `num`."""
    G = geo["puma"]
    gi = G.assign(pid=G.st * 100000 + G.puma).set_index("pid")
    den = d.edu.to_numpy() if den is None else den
    n = len(cz_w)
    pos = pd.Series(np.arange(n), index=cz_w.index.to_numpy(float))
    ci = d.pid.map(gi.czone).map(pos)
    cls = d.res_cls.to_numpy()
    ok = den & ci.notna().to_numpy() & ((cls == "city") | (cls == "suburb"))
    code = (ci.to_numpy()[ok].astype(int) * 2 + (cls[ok] == "suburb")).astype(int)
    fl = np.asarray(num)[ok].astype(float)
    rs = rep.sub(ok)
    wz = cz_w.to_numpy(float)

    def fn(w: np.ndarray) -> np.ndarray:
        E = np.bincount(code, weights=w, minlength=2 * n)
        K = np.bincount(code, weights=w * fl, minlength=2 * n)
        with np.errstate(invalid="ignore", divide="ignore"):
            sh = np.where(E > 0, K / E, np.nan)
            c, sb = sh[0::2], sh[1::2]
            good = np.isfinite(c) & np.isfinite(sb) & (c > 0) & (sb > 0)
            return np.where(good, np.log(sb / c), np.nan)

    est, reps = rs.run(fn)

    def mean(t: np.ndarray) -> float:
        m = np.isfinite(t)
        return float(np.sum(wz[m] * t[m]) / np.sum(wz[m]))

    m0 = mean(est)
    means = np.array([mean(r) for r in reps])
    have = np.isfinite(est)
    return dict(est=est, reps=reps, value=m0, se_survey=float(sdr_se(np.asarray(m0), means)),
                n_cz=int(have.sum()), weight_share=float(wz[have].sum() / wz.sum()), n_sample=int(ok.sum()))


# --------------------------------------------------------------------------
# (e) commuting
# --------------------------------------------------------------------------

def commute_groups(d: pd.DataFrame) -> dict[str, np.ndarray]:
    age = d.AGEP.between(25, 64).to_numpy()
    k12, pub = d.k12.to_numpy(), d.pub.to_numpy()
    f = d.female.to_numpy()
    cls = d.res_cls.to_numpy()
    metro = (cls == "city") | (cls == "suburb")
    bach = d.SCHL.ge(21).to_numpy()              # bachelor's degree or higher
    G: dict[str, np.ndarray] = {}
    base = {"T": age & k12, "Tpub": age & k12 & pub, "O": age & ~k12,
            "Obach": age & ~k12 & bach, "Otherpub": age & ~k12 & pub}
    for k, m in base.items():
        G[k] = m
        G[f"{k}_f"] = m & f
        G[f"{k}_m"] = m & ~f
        G[f"{k}@metro"] = m & metro
        G[f"{k}_f@metro"] = m & f & metro
        G[f"{k}_m@metro"] = m & ~f & metro
        G[f"{k}@city"] = m & (cls == "city")
        G[f"{k}@suburb"] = m & (cls == "suburb")
    return G


def commute_specs(d: pd.DataFrame, groups: dict[str, np.ndarray]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    jw = d.JWMNP.to_numpy(float)
    comm = d.comm.to_numpy()
    pow_ok = d.pow_ok.to_numpy()
    out = d.out_pow.to_numpy()
    home = d.home.to_numpy()
    sp: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for k, m in groups.items():
        mc, mp = (m & comm).astype(float), (m & pow_ok).astype(float)
        sp[f"{k}|mean"] = (jw * mc, mc)
        sp[f"{k}|gt30"] = ((jw > 30) * mc, mc)
        sp[f"{k}|ge30"] = ((jw >= 30) * mc, mc)      # 30 minutes or more (reported times heap at 30)
        sp[f"{k}|eq30"] = ((jw == 30) * mc, mc)
        sp[f"{k}|home"] = ((home & m).astype(float), m.astype(float))
        sp[f"{k}|out"] = ((out & m).astype(float), mp)
    return sp


def _metro_class(v: np.ndarray) -> np.ndarray:
    return (v == "city") | (v == "suburb")


def fe_extra_minutes(d: pd.DataFrame, rep: Rep, mask: np.ndarray, kind: str = "pow", se: bool = True
                     ) -> dict[str, Any]:
    """Extra one-way minutes of working outside the residence area: weighted regression of JWMNP on
    the crossing indicator, sex, age and age squared with residence-PUMA fixed effects (the
    comparison is between workers who live in the same PUMA).
        kind 'pow'        crossing = place of work in another POWPUMA (or state) than the residence's
        kind 'type'       crossing = the place-of-work POWPUMA is of the other class (city / suburb) than
                          the POWPUMA that contains the residence PUMA (`res_pow_cls` against `wk_cls`:
                          residence and workplace on the SAME geography, both ends city or suburb)
        kind 'type_mixed' the earlier, superseded definition, kept to show the artifact: class of the
                          residence PUMA (`res_cls`) against class of the place-of-work POWPUMA (two
                          geographies: 19% of metropolitan teachers live in a PUMA whose class differs
                          from its POWPUMA's, so 24% were coded as crossing against 9% on one geography)
    `se=False` skips the replicate-weight standard error (point estimate only)."""
    if kind == "pow":
        keep = mask & d.pow_ok.to_numpy()
        x = d.out_pow.to_numpy()[keep].astype(float)
    else:
        rc = (d.res_pow_cls if kind == "type" else d.res_cls).to_numpy()
        wc = d.wk_cls.to_numpy()
        keep = mask & d.pow_ok.to_numpy() & _metro_class(rc) & _metro_class(wc)
        x = (rc[keep] != wc[keep]).astype(float)
    s = d[keep]
    ids, g = np.unique(s.pid.to_numpy(), return_inverse=True)
    age = s.AGEP.to_numpy(float)
    cols = [x]
    fem = s.female.to_numpy(float)
    if 0.0 < fem.mean() < 1.0:                       # a sex-specific sample has no sex control
        cols.append(fem)
    X = np.column_stack(cols + [age - 45.0, (age - 45.0) ** 2 / 100.0])
    y = s.JWMNP.to_numpy(float)
    wt = rep.w0[keep]
    if se:
        rs = rep.sub(keep)
        est, reps = rs.run(lambda w: fe_coefs(y, X, g, w, len(ids)))
        se_, reps0 = float(sdr_se(est, reps)[0]), reps[:, 0].copy()
    else:
        est, se_, reps0 = fe_coefs(y, X, g, wt, len(ids)), float("nan"), None
    return dict(coef=float(est[0]), se=se_, n=int(keep.sum()),
                n_cross=int((x > 0).sum()), share_cross=float(wt[x > 0].sum() / wt.sum()),
                n_puma=int(len(ids)), kind=kind, _reps=reps0)


def commute_block(d: pd.DataFrame, rep: Rep) -> dict[str, Any]:
    groups = commute_groups(d)
    keep_any = np.zeros(len(d), bool)
    for m in groups.values():
        keep_any |= m
    sub, rs = d[keep_any].reset_index(drop=True), rep.sub(keep_any)
    g_sub = {k: m[keep_any] for k, m in groups.items()}
    R = Ratios(rs, commute_specs(sub, g_sub))
    out: dict[str, Any] = dict(n_sample={k: int(m.sum()) for k, m in g_sub.items()}, stats={})
    for k in g_sub:
        out["stats"][k] = {s: R.get(f"{k}|{s}") for s in ("mean", "gt30", "ge30", "eq30", "home", "out")}
    out["diff_T_minus_O"] = {}
    for suf in ("", "_f", "_m", "@metro", "@city", "@suburb"):
        out["diff_T_minus_O"][suf or "all"] = {s: R.derive(lambda g, s=s, suf=suf: g(f"T{suf}|{s}") - g(f"O{suf}|{s}"))
                                               for s in ("mean", "gt30", "ge30", "home", "out")}
    out["diff_T_minus_Obach"] = {s: R.derive(lambda g, s=s: g(f"T|{s}") - g(f"Obach|{s}"))
                                 for s in ("mean", "gt30", "ge30", "home", "out")}
    # extra commute of crossing, residence-PUMA fixed effects
    fe: dict[str, Any] = {}
    for lab, m in (("T", groups["T"]), ("T@metro", groups["T@metro"]), ("Tpub@metro", groups["Tpub@metro"]),
                   ("T_f@metro", groups["T_f@metro"]), ("T_m@metro", groups["T_m@metro"]),
                   ("T@suburb", groups["T@suburb"]), ("T@city", groups["T@city"]),
                   ("O", groups["O"]), ("O@metro", groups["O@metro"]),
                   ("Obach@metro", groups["Obach@metro"])):
        fe[lab] = fe_extra_minutes(d, rep, m, "pow")
        LOG.info("extra minutes of working outside the residence POWPUMA, %-12s %.2f (%.2f), n = %d",
                 lab, fe[lab]["coef"], fe[lab]["se"], fe[lab]["n"])
    # working in a POWPUMA of the other class (city / suburb): residence and workplace both on the POWPUMA
    # geography, so the sample is every teacher (other worker) 25-64 whose residence and workplace POWPUMAs
    # are of class city or suburb
    for lab, m in (("T@metro", groups["T"]), ("O@metro", groups["O"])):
        fe[f"{lab}:type"] = fe_extra_minutes(d, rep, m, "type")
        LOG.info("extra minutes of working in the other class of POWPUMA, %-9s %.2f (%.2f), n = %d",
                 lab, fe[f"{lab}:type"]["coef"], fe[f"{lab}:type"]["se"], fe[f"{lab}:type"]["n"])
    # robustness: only residence PUMAs whose own class equals the class of their POWPUMA
    same = (d.res_cls == d.res_pow_cls).to_numpy()
    fe["T@metro:type_same_class"] = fe_extra_minutes(d, rep, groups["T"] & same, "type")
    # the superseded definition (class of the residence PUMA against the class of the work POWPUMA), for the record
    fe["T@metro:type_mixed"] = fe_extra_minutes(d, rep, groups["T@metro"], "type_mixed")
    for k in ("T@metro:type_same_class", "T@metro:type_mixed"):
        LOG.info("extra minutes, %-24s %.2f (%.2f), n = %d, share crossing %.3f", k, fe[k]["coef"], fe[k]["se"],
                 fe[k]["n"], fe[k]["share_cross"])
    out["fe"] = fe
    # raw difference of mean commute, outside minus inside the residence POWPUMA, teachers in metro areas
    for lab in ("T@metro", "O@metro"):
        m = groups[lab] & d.pow_ok.to_numpy()
        jw = d.JWMNP.to_numpy(float)
        o = d.out_pow.to_numpy()
        r2 = rep.sub(m)

        def fn(w: np.ndarray, jw=jw[m], o=o[m]) -> np.ndarray:
            return np.array([w @ (jw * o) / (w @ o), w @ (jw * ~o) / (w @ ~o)])

        est, reps = r2.run(fn)
        diff = est[0] - est[1]
        dse = float(sdr_se(np.asarray(diff), reps[:, 0] - reps[:, 1]))
        out.setdefault("raw_outside_minus_inside", {})[lab] = dict(
            outside=float(est[0]), inside=float(est[1]), diff=float(diff), se=dse, n=int(m.sum()))
    out["_ratios"] = R
    return out


def flows_block(d: pd.DataFrame, rep: Rep) -> dict[str, Any]:
    """Residence against workplace of teachers on one geography: both by the class (city / suburb /
    other) of the POWPUMA (the finest area on which the PUMS reports the place of work), the PUMS
    counterpart of the district-based net_commute_gap, and the same with EDU workers by residence
    (so that the difference is the composition term)."""
    ok = (d.POWSP.between(1, 56) & d.POWPUMA.gt(1) & d.res_powpuma.notna()).to_numpy()
    sub, rs = d[ok].reset_index(drop=True), rep.sub(ok)
    rc, wc = sub.res_pow_cls.to_numpy(), sub.wk_cls.to_numpy()
    sel = {"k12": sub.k12.to_numpy(), "pub_k12": (sub.k12 & sub.pub).to_numpy(), "edu": sub.edu.to_numpy()}
    specs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    one = np.ones(len(sub))
    for c in CLASSES:
        specs[f"all|res|{c}"] = ((rc == c).astype(float), one)               # all employed residents of the class
    for k, m in sel.items():
        for c in CLASSES:
            specs[f"{k}|res|{c}"] = (((rc == c) & m).astype(float), one)    # shares of all workers
            specs[f"{k}|work|{c}"] = (((wc == c) & m).astype(float), one)
    # cross-tabulation of residence class by workplace class for the K-12 teachers
    for a in CLASSES:
        for b in CLASSES:
            specs[f"k12|x|{a}|{b}"] = (((rc == a) & (wc == b) & sel["k12"]).astype(float), one)
            specs[f"pub_k12|x|{a}|{b}"] = (((rc == a) & (wc == b) & sel["pub_k12"]).astype(float), one)
    R = Ratios(rs, specs)
    out: dict[str, Any] = {"n_sample": int(ok.sum())}
    for k in sel:
        out[k] = {"res": {c: R.get(f"{k}|res|{c}") for c in CLASSES},
                  "work": {c: R.get(f"{k}|work|{c}") for c in CLASSES}}
        out[k]["gap_res"] = R.derive(lambda g, k=k: np.log(g(f"{k}|res|suburb")) - np.log(g(f"{k}|res|city")))
        out[k]["gap_work"] = R.derive(lambda g, k=k: np.log(g(f"{k}|work|suburb")) - np.log(g(f"{k}|work|city")))
        # the group's share of the employed residents of the class: the counterpart of gap_teach_res_emp
        out[k]["gap_res_emp"] = R.derive(
            lambda g, k=k: (np.log(g(f"{k}|res|suburb")) - np.log(g("all|res|suburb")))
            - (np.log(g(f"{k}|res|city")) - np.log(g("all|res|city"))))
    for k in ("k12", "pub_k12"):
        out[k]["net_commute_gap"] = R.derive(
            lambda g, k=k: (np.log(g(f"{k}|res|suburb")) - np.log(g(f"{k}|res|city")))
            - (np.log(g(f"{k}|work|suburb")) - np.log(g(f"{k}|work|city"))))
        out[k]["share_cross_class"] = R.derive(
            lambda g, k=k: 1.0 - sum(g(f"{k}|x|{c}|{c}") for c in CLASSES) / sum(
                g(f"{k}|x|{a}|{b}") for a in CLASSES for b in CLASSES))
        out[k]["share_metro_cross"] = R.derive(
            lambda g, k=k: (g(f"{k}|x|city|suburb") + g(f"{k}|x|suburb|city")) / sum(
                g(f"{k}|x|{a}|{b}") for a in ("city", "suburb") for b in ("city", "suburb")))
        out[k]["x"] = {f"{a}->{b}": R.get(f"{k}|x|{a}|{b}") for a in CLASSES for b in CLASSES}
        # of the teachers who work in a class of POWPUMA, the share who live in each class
        out[k]["live_given_work"] = {
            f"work_{b}": {f"live_{a}": R.derive(lambda g, k=k, a=a, b=b: g(f"{k}|x|{a}|{b}") / sum(
                g(f"{k}|x|{c}|{b}") for c in CLASSES)) for a in CLASSES} for b in ("city", "suburb")}
    # EDU workers by residence against public K-12 teachers by workplace (the district measure's PUMS analogue)
    out["edu_res_vs_pubk12_work"] = dict(
        gap_edu_res=R.derive(lambda g: np.log(g("edu|res|suburb")) - np.log(g("edu|res|city"))),
        gap_pubk12_work=R.derive(lambda g: np.log(g("pub_k12|work|suburb")) - np.log(g("pub_k12|work|city"))),
        net_commute_gap=R.derive(lambda g: (np.log(g("edu|res|suburb")) - np.log(g("edu|res|city")))
                                 - (np.log(g("pub_k12|work|suburb")) - np.log(g("pub_k12|work|city")))),
        composition=R.derive(lambda g: (np.log(g("edu|res|suburb")) - np.log(g("edu|res|city")))
                             - (np.log(g("pub_k12|res|suburb")) - np.log(g("pub_k12|res|city")))))
    return out


CUTS = ((0.25, 0.5), (0.5, 0.5), (0.75, 0.5), (0.9, 0.5), (0.5, 0.9))


def class_sensitivity(d: pd.DataFrame, rep: Rep, geo: dict[str, Any], cz_w: pd.Series,
                      cuts: tuple[tuple[float, float], ...] = CUTS) -> list[dict[str, Any]]:
    """The composition term and the class-of-area extra commute (POWPUMA classes at both ends) under other
    cuts of the PUMA classification (principal-city share >= pc_cut is 'city'; metropolitan share >=
    metro_cut and below the cut is 'suburb').  `d` and `geo` are left with the original classification."""
    main = (geo["info"]["pc_cut"], geo["info"]["metro_cut"])
    rows: list[dict[str, Any]] = []
    age = d.AGEP.between(25, 64).to_numpy()
    for pc_cut, metro_cut in cuts:
        set_classes(geo, pc_cut, metro_cut)
        attach_classes(d, geo)
        e = d.edu.to_numpy()
        sub, rs = d[e].reset_index(drop=True), rep.sub(e)
        cls = sub.res_cls.to_numpy()
        sp: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for c in ("city", "suburb"):
            m = (cls == c).astype(float)
            sp[f"{c}_k12"] = (sub.k12.to_numpy() * m, m)
            sp[f"{c}_pub"] = ((sub.k12 & sub.pub).to_numpy() * m, m)
        R = Ratios(rs, sp)
        # extra minutes of working in a POWPUMA of the other class: residence and workplace on the POWPUMA geography
        fe = fe_extra_minutes(d, rep, age & d.k12.to_numpy(), "type")
        feo = fe_extra_minutes(d, rep, age & ~d.k12.to_numpy(), "type", se=False)        # other workers: point estimate
        wcz = cz_composition(d, rep, geo, cz_w)
        tz = cz_terms(d, rep, geo, cz_w, (d.k12 & d.pub).to_numpy())                       # the same term CZ by CZ
        if abs(tz["value"] - wcz["pub_k12"]["value"]) > 1e-10 or abs(tz["se_survey"] - wcz["pub_k12"]["se"]) > 1e-10:
            raise SystemExit("cz_terms does not reproduce cz_composition")
        rows.append(dict(
            pc_cut=pc_cut, metro_cut=metro_cut, main=(pc_cut, metro_cut) == main,
            n_puma=geo["info"]["n_puma_by_class"], pop_by_class=geo["info"]["pop_by_class"],
            share_edu_workers_city=float(np.mean(cls == "city")),
            s_pub_city=R.get("city_pub"), s_pub_suburb=R.get("suburb_pub"),
            delta_log_s_pub=R.derive(lambda g: np.log(g("suburb_pub")) - np.log(g("city_pub"))),
            delta_log_s_k12=R.derive(lambda g: np.log(g("suburb_k12")) - np.log(g("city_k12"))),
            within_cz=wcz, _terms=tz["est"],
            extra_minutes_class=dict(coef=fe["coef"], se=fe["se"], n=fe["n"], share_cross=fe["share_cross"]),
            extra_minutes_class_other=dict(coef=feo["coef"], n=feo["n"], share_cross=feo["share_cross"])))
        LOG.info("classes pc >= %.2f, metro >= %.2f: Delta log s (public K-12) %.3f, extra minutes %.2f (others %.2f)",
                 pc_cut, metro_cut, rows[-1]["delta_log_s_pub"]["value"], fe["coef"], feo["coef"])
    set_classes(geo, *main)
    attach_classes(d, geo)
    return rows


# --------------------------------------------------------------------------
# Wages of teachers, and (f) the commuting cost r_k
# --------------------------------------------------------------------------

def wage_block(d: pd.DataFrame, rep: Rep) -> dict[str, Any]:
    """Median hourly wage and salary earnings (2018 dollars) of full-time (usual hours >= 35),
    full-year (48-52 weeks) workers 25-64: WAGP * ADJINC / (WKHP * weeks), weeks 51 for the 50-52
    bracket and 48.5 for 48-49; and their mean annual hours worked WKHP * weeks (the hours over which
    annual earnings are earned, for tau_k)."""
    age = d.AGEP.between(25, 64)
    fty = (d.WKW.isin([1, 2]) & d.WKHP.ge(35) & age).to_numpy()
    ftyr = (fty & d.WAGP.gt(0).to_numpy())
    weeks = np.where(d.WKW.to_numpy() == 1, 51.0, 48.5)
    hours = d.WKHP.to_numpy(float) * weeks
    h = (d.WAGP * d.ADJINC).to_numpy() / hours
    out: dict[str, Any] = {}
    k12, pub, bach = d.k12.to_numpy(), d.pub.to_numpy(), d.SCHL.ge(21).to_numpy()
    for lab, m, mh in (("k12", ftyr & k12, fty & k12), ("k12_public", ftyr & k12 & pub, fty & k12 & pub),
                       ("other", ftyr & ~k12, fty & ~k12), ("other_bachelor", ftyr & ~k12 & bach, fty & ~k12 & bach)):
        x, rs = h[m], rep.sub(m)
        order = np.argsort(x, kind="stable")
        xs = x[order]

        def fn(w: np.ndarray, xs=xs, order=order) -> float:
            ws = w[order]
            c = np.cumsum(ws) - 0.5 * ws
            return float(np.interp(0.5 * ws.sum(), c, xs))

        est, reps = rs.run(fn)
        lx = np.log(x)
        gm, gr = rs.run(lambda w, lx=lx: float(np.exp(w @ lx / w.sum())))
        hh, rh = rep.sub(mh), hours[mh]
        ha, hr = hh.run(lambda w, rh=rh: float(w @ rh / w.sum()))
        out[lab] = dict(median_hourly=float(est), se_median=float(sdr_se(est, reps)),
                        geometric_mean_hourly=float(gm), se_geometric_mean=float(sdr_se(gm, gr)),
                        share_at_median=float(rs.w0[np.isclose(x, est)].sum() / rs.w0.sum()), n=int(m.sum()),
                        annual_hours_mean=float(ha), se_annual_hours_mean=float(sdr_se(ha, hr)), n_hours=int(mh.sum()))
    return out


def rk_value(delta_min: float, days: float, vtt: float, cbar: float) -> float:
    """r_k = 2 trips a day * extra one-way hours * days a year * value of an hour / Cbar."""
    return 2.0 * delta_min / 60.0 * days * vtt / cbar


def tau_value(delta_min: float, days: float, vot_share: float, hours_year: float) -> float:
    """tau_k = (value of the extra commuting time) / (annual earnings)
             = vot_share * w * (2 * extra one-way hours * days) / (w * annual hours)
    with the value of an hour of travel a share `vot_share` of the worker's own hourly wage w and annual
    earnings w * annual hours: the share of income lost (M2_SPEC_v2: income x (1 - tau_k) when r != w)."""
    return vot_share * 2.0 * delta_min / 60.0 * days / hours_year


def rk_block(commute: dict[str, Any], wages: dict[str, Any], sens_cls: list[dict[str, Any]]) -> dict[str, Any]:
    """Section (f).  The extra one-way commute of working in the other location is estimated
    two ways (weighted regression of the commute on the crossing indicator with residence-PUMA fixed
    effects): crossing into another POWPUMA (county based) and working in a POWPUMA of the other class,
    city or suburb, residence and workplace on the same POWPUMA geography.  Both measure crossings of a
    county line, so they may OVERSTATE the increment of working in the other location where city and suburb
    are in one county; the Nashville out-of-county figure of Santelli and Grissom (13.8 minutes) anchors the
    low end.  Central = the average of the two estimates.

    The annual time cost is 2 x extra hours x days x value of an hour.  r_k divides it by mean consumption:
    per CE consumer unit (Cbar = $61,224.13: the denominator of r_p and the unit the task fixed) and per
    ADULT (Cbar / 1.9 adults = $32,223: the model's agents are single adults and a commuting cost is a cost
    per worker, so per-adult consumption is the matching denominator).  tau_k divides it by annual earnings
    (value of an hour = a share of the worker's own wage), for teachers and for other workers."""
    cbar = hc.CE_TOTAL_EXP
    cbar_ad = cbar / CE_ADULTS
    cpi = CPI_U[2018] / CPI_U[2015]
    vtt = {k: VTT_2015[k] * cpi for k in ("central", "low", "high")}
    share = dict(central=VTT_2015["share_central"], low=VTT_2015["share_low"], high=VTT_2015["share_high"])
    fe = commute["fe"]
    a, b = fe["T@metro"], fe["T@metro:type"]
    reps = 0.5 * (a["_reps"] + b["_reps"])
    cen = 0.5 * (a["coef"] + b["coef"])
    se = float(sdr_se(np.asarray(cen), reps))
    # range: the smallest and largest teacher estimate over the specifications (crossing a POWPUMA, the class
    # of area under each PUMA class cut, Nashville); the 95% intervals of single estimates are narrower
    cand = [a["coef"]] + [r_["extra_minutes_class"]["coef"] for r_ in sens_cls] + [NASHVILLE_EXTRA_MIN]
    delta = dict(central=cen, low=min(cand), high=max(cand),
                 ci_low=min(a["coef"] - 1.96 * a["se"], b["coef"] - 1.96 * b["se"]),
                 ci_high=max(a["coef"] + 1.96 * a["se"], b["coef"] + 1.96 * b["se"]))
    ao, bo = fe["O@metro"], fe["O@metro:type"]
    reps_o = 0.5 * (ao["_reps"] + bo["_reps"])
    cen_o = 0.5 * (ao["coef"] + bo["coef"])
    cand_o = [ao["coef"]] + [r_["extra_minutes_class_other"]["coef"] for r_ in sens_cls]
    delta_o = dict(central=cen_o, low=min(cand_o), high=max(cand_o), se=float(sdr_se(np.asarray(cen_o), reps_o)))
    hrs = dict(teachers=wages["k12"]["annual_hours_mean"], other_workers=wages["other"]["annual_hours_mean"])
    days = dict(DAYS_TEACHER)
    wage_vtt = 0.5 * wages["k12"]["median_hourly"]
    out: dict[str, Any] = dict(
        cbar_cu=cbar, adults_per_cu=CE_ADULTS, adults_per_cu_range=list(CE_ADULTS_RANGE), cbar_adult=cbar_ad,
        ce_persons=CE_PERSONS, ce_children_under_18=CE_CHILDREN, ce_source=CE_URL,
        cpi_ratio_2018_over_2015=cpi, vtt_2018=vtt, vtt_2015=VTT_2015, days=days, days_other=dict(DAYS_OTHER),
        delta_min=delta, delta_se=se, delta_min_other=delta_o, annual_hours=hrs, delta_estimates={
            "pow_crossing_teachers_metro": dict(coef=a["coef"], se=a["se"], n=a["n"]),
            "class_crossing_teachers_metro": dict(coef=b["coef"], se=b["se"], n=b["n"]),
            "pow_crossing_all_other_workers_metro": dict(coef=ao["coef"], se=ao["se"], n=ao["n"]),
            "class_crossing_all_other_workers_metro": dict(coef=bo["coef"], se=bo["se"], n=bo["n"]),
            "santelli_grissom_2024_out_of_county": dict(coef=NASHVILLE_EXTRA_MIN, se=None, n=None)},
        teacher_median_hourly_wage=wages["k12"]["median_hourly"], vtt_half_teacher_wage=wage_vtt)
    out["central_per_cu"] = rk_value(delta["central"], days["central"], vtt["central"], cbar)
    out["central_per_adult"] = rk_value(delta["central"], days["central"], vtt["central"], cbar_ad)
    out["central_se_from_delta_per_cu"] = out["central_per_cu"] * se / delta["central"]
    out["central_se_from_delta_per_adult"] = out["central_per_adult"] * se / delta["central"]
    out["low_per_cu"] = rk_value(delta["low"], days["low"], vtt["low"], cbar)
    out["high_per_cu"] = rk_value(delta["high"], days["high"], vtt["high"], cbar)
    out["low_per_adult"] = rk_value(delta["low"], days["low"], vtt["low"], cbar_ad)
    out["high_per_adult"] = rk_value(delta["high"], days["high"], vtt["high"], cbar_ad)
    out["annual_time_cost_central_usd"] = out["central_per_cu"] * cbar
    out["annual_hours_central"] = 2.0 * delta["central"] / 60.0 * days["central"]
    # tau_k: time cost over the worker's annual earnings
    tau: dict[str, Any] = {}
    for grp, dl, dy, hr in (("teachers", delta, DAYS_TEACHER, hrs["teachers"]),
                            ("other_workers", delta_o, DAYS_OTHER, hrs["other_workers"])):
        t = dict(extra_one_way_minutes=[dl["low"], dl["central"], dl["high"]],
                 days=[dy["low"], dy["central"], dy["high"]], annual_hours=hr,
                 value_of_time_share_of_own_wage=[share["low"], share["central"], share["high"]],
                 central=tau_value(dl["central"], dy["central"], share["central"], hr),
                 low=tau_value(dl["low"], dy["low"], share["low"], hr),
                 high=tau_value(dl["high"], dy["high"], share["high"], hr),
                 central_with_2080_hours=tau_value(dl["central"], dy["central"], share["central"], 2080.0),
                 se_from_extra_minutes=tau_value(dl["central"], dy["central"], share["central"], hr) * (se if grp == "teachers" else dl["se"]) / dl["central"])
        tau[grp] = t
    out["tau_k"] = tau
    adj = {1.8: CE_ADULTS_RANGE[0], 2.0: CE_ADULTS_RANGE[1]}
    sens = {}
    for name, (dd, dy, vv, sh) in {
            "extra minutes at its low end (smallest estimate)": (delta["low"], days["central"], vtt["central"], share["central"]),
            "extra minutes at its high end (largest estimate)": (delta["high"], days["central"], vtt["central"], share["central"]),
            "value of time at 35% of hourly household income": (delta["central"], days["central"], vtt["low"], share["low"]),
            "value of time at 60% of hourly household income": (delta["central"], days["central"], vtt["high"], share["high"]),
            "value of time at 50% of the teacher median hourly wage": (delta["central"], days["central"], wage_vtt, share["central"]),
            "days 180": (delta["central"], 180.0, vtt["central"], share["central"]),
            "days 190": (delta["central"], 190.0, vtt["central"], share["central"]),
            "days 250 (a full-year worker)": (delta["central"], 250.0, vtt["central"], share["central"]),
            "extra minutes: crossing a POWPUMA (teachers, metro)": (a["coef"], days["central"], vtt["central"], share["central"]),
            "extra minutes: crossing a POWPUMA of the other class, city / suburb (teachers)": (b["coef"], days["central"], vtt["central"], share["central"]),
            "extra minutes: Santelli-Grissom out-of-county, 13.8": (NASHVILLE_EXTRA_MIN, days["central"], vtt["central"], share["central"])}.items():
        sens[name] = dict(delta_min=dd, days=dy, vtt=vv, r_k_per_cu=rk_value(dd, dy, vv, cbar),
                          r_k_per_adult=rk_value(dd, dy, vv, cbar_ad), tau_k=tau_value(dd, dy, sh, hrs["teachers"]))
    for r_ in sens_cls:                              # the class-of-area estimate under other PUMA class cuts
        if not r_["main"]:
            dd = 0.5 * (a["coef"] + r_["extra_minutes_class"]["coef"])
            sens[f"central with the class-of-area minutes at PUMA class cut pc >= {r_['pc_cut']}, metro >= {r_['metro_cut']}"] = dict(
                delta_min=dd, days=days["central"], vtt=vtt["central"], r_k_per_cu=rk_value(dd, days["central"], vtt["central"], cbar),
                r_k_per_adult=rk_value(dd, days["central"], vtt["central"], cbar_ad),
                tau_k=tau_value(dd, days["central"], share["central"], hrs["teachers"]))
    for na, av in adj.items():                       # the rounding of the published CE persons / children counts
        sens[f"adults per consumer unit {na} (CE counts are rounded to one decimal)"] = dict(
            delta_min=delta["central"], days=days["central"], vtt=vtt["central"], r_k_per_cu=out["central_per_cu"],
            r_k_per_adult=rk_value(delta["central"], days["central"], vtt["central"], cbar / av),
            tau_k=out["tau_k"]["teachers"]["central"])
    out["sensitivity"] = sens
    # vehicle operating costs of the extra miles (not in r_k or tau_k): size only, with assumed speed and cost per mile
    mph, cpm = 35.0, (0.10, 0.30)
    miles = 2.0 * delta["central"] / 60.0 * mph * days["central"]
    out["vehicle_cost_note"] = dict(mph=mph, miles_per_year=miles, usd_per_mile=list(cpm), usd=[miles * c for c in cpm],
                                    per_cu=[miles * c / cbar for c in cpm], per_adult=[miles * c / cbar_ad for c in cpm])
    return out


def nyc_check(d: pd.DataFrame, rep: Rep, geo: dict[str, Any]) -> dict[str, Any]:
    """Validation of the PUMS against the NYC Independent Budget Office count (67% of New York City teachers
    live in the five boroughs, 2011-12): commuters whose place of work is in one of the five NYC POWPUMAs
    (the POWPUMAs that contain the 55 `NYC-` PUMAs), and the share living in those PUMAs."""
    G = geo["puma"]
    nyc = G[(G.st == 36) & G.name.str.startswith("NYC")]
    pid = set(nyc.st * 100000 + nyc.puma)
    pw = set(nyc.powsp.astype(int) * 100000 + nyc.powpuma.astype(int))
    base = ((d.JWMNP > 0) & d.POWSP.eq(36) & d.wk_key.isin(pw)).to_numpy()
    sub, rs = d[base].reset_index(drop=True), rep.sub(base)
    live = sub.pid.isin(pid).to_numpy().astype(float)
    same = (sub.res_pow_key == sub.wk_key).to_numpy().astype(float)
    out: dict[str, Any] = {"n_nyc_pumas": int(len(nyc)), "n_nyc_powpumas": len(pw)}
    specs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for lab, g in (("k12", sub.k12.to_numpy()), ("pub_k12", (sub.k12 & sub.pub).to_numpy()),
                   ("all_workers", np.ones(len(sub), bool))):
        gm = g.astype(float)
        specs[f"{lab}|live"] = (live * gm, gm)
        specs[f"{lab}|same_borough"] = (same * gm, gm)
        out[f"n_{lab}"] = int(g.sum())
        out[f"workers_{lab}"] = float(rs.w0[g].sum())
    R = Ratios(rs, specs)
    for lab in ("k12", "pub_k12", "all_workers"):
        out[f"{lab}_live_in_nyc"] = R.get(f"{lab}|live")
        out[f"{lab}_live_in_workplace_powpuma"] = R.get(f"{lab}|same_borough")
    return out


def cmd_fetch(only: str | None) -> None:
    if only in (None, "acs"):
        fetch_acs()
    if only in (None, "geo"):
        fetch_geo()
    if only in (None, "pums"):
        fetch_pums()


def selftest() -> None:
    """Unit checks of the helpers on synthetic data."""
    rng = np.random.default_rng(1)
    n, G = 400, 12
    g = rng.integers(0, G, n)
    x, z = rng.normal(size=n), rng.normal(size=n)
    w = rng.integers(1, 20, n).astype(float)
    y = 2.0 * x - 1.0 * z + rng.normal(size=G)[g] + rng.normal(size=n)
    X = np.column_stack([x, z])
    b = fe_coefs(y, X, g, w, G)
    sw = np.sqrt(w)
    b2 = np.linalg.lstsq(np.column_stack([X, np.eye(G)[g]]) * sw[:, None], y * sw, rcond=None)[0][:2]
    assert np.allclose(b, b2), (b, b2)                       # within transformation = dummy-variable WLS
    W = rng.integers(1, 20, (n, N_REP))
    w0 = W.mean(axis=1)
    num = (rng.random(n) < 0.3).astype(float)
    R = Ratios(Rep(w0, W), {"a": (num, np.ones(n)), "b": (num * z, np.ones(n))})
    est = (w0 @ num) / w0.sum()
    reps = np.array([(W[:, r] @ num) / W[:, r].sum() for r in range(N_REP)])
    assert np.isclose(R.get("a")["value"], est)
    assert np.isclose(R.get("a")["se"], np.sqrt(4.0 / N_REP * ((reps - est) ** 2).sum()))    # Census SDR formula
    d_ = R.derive(lambda gt: gt("a") - gt("a"))
    assert abs(d_["value"]) < 1e-15 and abs(d_["se"]) < 1e-15
    assert np.isclose(wquantile(np.array([1.0, 2.0, 3.0]), np.ones(3), 0.5), 2.0)
    assert np.isclose(rk_value(15, 180, 20.0, 50000.0), 0.036)                              # 2 x 0.25 h x 180 d x $20 / 50,000
    assert np.isclose(tau_value(15, 180, 0.5, 2000.0), 0.5 * 2 * 0.25 * 180 / 2000)         # share x hours lost / hours worked
    assert np.isclose(CE_ADULTS, 1.9) and np.isclose(hc.CE_TOTAL_EXP / CE_ADULTS, 32223.23, atol=0.01)
    # --- variant_panel: residents counted once -------------------------------------------------------------
    toy = pd.DataFrame({
        "czone": [1.0] * 5, "loc_id": [1, 2, 2, 2, 1], "leaid": list("abcde"),
        "sd_layer": ["unified", "elementary", "elementary", "secondary", "secondary"],
        "teachers_k12": [10.0, 4.0, 6.0, 12.0, 5.0], "pupils_k12": [100.0, 40.0, 60.0, 90.0, 50.0],
        "m_edu": [8.0, 3.0, 5.0, 8.0, 4.0], "m_emp": [80.0, 30.0, 50.0, 80.0, 40.0], "m_fte": [10.0, 4.0, 6.0, 12.0, 5.0],
        "m_pup": [100.0, 40.0, 60.0, 90.0, 50.0], "m_vedu": [1.0] * 5, "m_vemp": [2.0] * 5,
        "a_edu": [8.0, 3.0, 5.0, 8.0, 4.0], "ok_match": True})
    one = variant_panel(toy, "no_overlap")
    assert one.m_edu.sum() == 16.0 and one.m_fte.sum() == 37.0 and one.m_pup.sum() == 340.0      # residents once, FTE / pupils from every LEA
    assert one.loc[one.sd_layer == "secondary", "m_edu"].eq(0).all() and one.a_edu.sum() == 16.0
    sec = variant_panel(toy, "no_overlap_sec")          # elementary residents dropped where the CZ location holds a secondary LEA
    assert sec.m_edu.tolist() == [8.0, 0.0, 0.0, 8.0, 4.0] and sec.m_fte.sum() == 37.0
    assert variant_panel(toy, "earnings_rule").m_edu.sum() == 28.0
    assert variant_panel(toy, "unified").m_edu.notna().sum() == 1
    # --- cz_terms: composition CZ by CZ ---------------------------------------------------------------------
    geo_t = {"puma": pd.DataFrame({"st": [1, 1, 1, 1], "puma": [1, 2, 3, 4], "czone": [10.0, 10.0, 20.0, 20.0]})}
    dd = pd.DataFrame({"pid": [100001, 100002, 100003, 100004] * 2, "res_cls": ["city", "suburb"] * 4,
                       "edu": True, "k12": [True, True, False, True, True, False, True, True]})
    dd["pid"] = [100001, 100002, 100003, 100004, 100001, 100002, 100003, 100004]
    dd["res_cls"] = ["city", "suburb", "city", "suburb"] * 2
    wts = np.array([2.0, 1.0, 3.0, 1.0, 2.0, 3.0, 1.0, 1.0])
    ct = cz_terms(dd, Rep(wts, np.tile(wts[:, None], (1, N_REP))), geo_t, pd.Series([5.0, 3.0], index=[10.0, 20.0]),
                  dd.k12.to_numpy())
    # CZ 10: city share = (2 + 2) / (2 + 2) = 1.0 (rows 0, 4); suburb = (1 + 0) / (1 + 3) = 0.25 (rows 1, 5)
    # CZ 20: city = (0 + 1) / (3 + 1) = 0.25 (rows 2, 6); suburb = (1 + 1) / (1 + 1) = 1.0 (rows 3, 7)
    assert np.allclose(ct["est"], [np.log(0.25), np.log(4.0)]), ct["est"]
    assert np.isclose(ct["value"], (5 * np.log(0.25) + 3 * np.log(4.0)) / 8.0) and ct["n_cz"] == 2
    a_, b_ = dict(value=0.1, se=0.03), dict(value=0.2, se=0.04)
    assert np.isclose(comb(a_, b_)["se"], 0.05) and np.isclose(comb(a_, b_, -1.0)["value"], -0.1)
    toml = ds.toml_val([[1.0, 0.5], [0.5, 1.0]])
    assert toml.startswith("[\n") and "0.5" in toml
    assert load_c24010.__doc__ and pcs(dict(value=0.123, se=0.004)) == "12.3% (0.4)"
    print("selftest ok")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    f = sub.add_parser("fetch", help="download and cache the raw inputs")
    f.add_argument("--only", choices=["acs", "pums", "geo"])
    sub.add_parser("selftest", help="unit checks of the helpers")
    b = sub.add_parser("build", help="compute the moments and write the reports")
    b.add_argument("--B", type=int, default=B_DEFAULT)
    b.add_argument("--seed", type=int, default=SEED)
    b.add_argument("--outdir", default=None, help="default: <here>/estimates/settings")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    if args.cmd == "fetch":
        cmd_fetch(args.only)
    elif args.cmd == "build":
        cmd_build(args)
    elif args.cmd == "selftest":
        selftest()
    else:
        cmd_fetch(None)
        cmd_build(argparse.Namespace(B=B_DEFAULT, seed=SEED, outdir=None))
    return 0


# --------------------------------------------------------------------------
# Literature
# --------------------------------------------------------------------------

LITERATURE: list[dict[str, Any]] = [
    dict(key="boyd2005",
         ref="Boyd, Lankford, Loeb and Wyckoff (2005), 'The Draw of Home: How Teachers' Preferences for Proximity "
             "Disadvantage Urban Schools', Journal of Policy Analysis and Management 24(1): 113-132 (NBER WP 9953, "
             "Sept. 2003)",
         url="https://www.nber.org/papers/w9953",
         read="NBER working paper PDF, text and Tables 2-5",
         facts=["New York State, teachers entering public school teaching 1999-2002: 61% (60.8%) took their first job "
                "within 15 miles of their hometown and 85% within 40 miles (Table 2; 90% within 40 miles in New York City)",
                "urban districts import teachers from the suburbs: 36% of urban teachers have a suburban hometown, "
                "6% of suburban teachers an urban hometown; the ratio of suburban-hometown teachers taking a city "
                "job to city-hometown teachers taking a suburban job runs from 2.6 (Buffalo) to 7.5 (New York City) "
                "(p. 7)"]),
    dict(key="reininger2012",
         ref="Reininger (2012), 'Hometown Disadvantage? It Depends on Where You're From', Educational Evaluation and "
             "Policy Analysis 34(2): 127-145",
         url="https://eric.ed.gov/?id=EJ963419",
         read="not read directly: the figure is as summarised at https://www.deansforimpact.org/about/news-and-blog/"
              "2017/03/06/13-miles-the-inherent-localism-of-teaching",
         facts=["NELS:88: the median distance between the high school a teacher attended in 10th grade and her "
                "residence eight years later is 13 miles, about a quarter of the distance moved by other professionals; "
                "4 miles for teachers who grew up in dense cities"]),
    dict(key="ibo_nyc",
         ref="New York City Independent Budget Office, 'Where Do NYC's Teachers and Principals Live Compared With Where "
             "They Work?'",
         url="https://www.ibo.nyc.ny.us/iboreports/printnycbtn19.pdf",
         read="the PDF could not be opened from this environment (the old IBO address now redirects to the IBO home page, "
              "ibo.nyc.gov); figures and the date (Aug. 15, 2013) as given in the search summary of the report",
         facts=["school year 2011-12: 67% of 70,328 New York City teachers live within the five boroughs; 81% of Staten "
                "Island teachers live in the borough where they work, 57% of Brooklyn teachers; teachers in Queens and "
                "the Bronx are the most likely to live outside the five boroughs"]),
    dict(key="santelli_grissom",
         ref="Santelli and Grissom (2024), 'A Bad Commute: Travel Time to Work Predicts Teacher Turnover and Other "
             "Workplace Outcomes', AERA Open 10; working paper: EdWorkingPaper No. 22-691 (Dec. 2022), 'A Bad Commute: "
             "Does Travel Time to Work Predict Teacher and Leader Turnover and Other Workplace Outcomes?'",
         url="https://edworkingpapers.com/ai22-691",
         read="EdWorkingPaper No. 22-691 (Dec. 2022), PDF text, pp. 10, 15, 16, 17, 22 (printed page numbers); the 0.8-1.0 "
              "percentage points per 5 minutes is as stated in the published AERA Open abstract (via search; the published "
              "article was not opened)",
         facts=["Metropolitan Nashville Public Schools (countywide), 2008-09 to 2019-20, home and school addresses with HERE "
                "travel times: mean one-way commute 21.8 minutes for teachers (p. 10); 19% of teachers travel 30 minutes or more "
                "(p. 15; compare the PUMS share of 30 minutes or more, section 4, not the share above 30); teachers who live out "
                "of the county commute 13.8 minutes longer than those who live in it (p. 15)",
                "working paper: teachers who commute 45 minutes or more each way are 10 percentage points more likely to transfer "
                "schools (overall rate 12%) and 3 percentage points more likely to exit the district (overall rate 7%) than "
                "teachers in the same school who commute 0-5 minutes (pp. 16-17); predicted transfer probability 8% against 18%, "
                "exit 8% against 11% (p. 22). Published version: each 5 minutes of one-way commute raises the transfer probability "
                "by 0.8-1.0 percentage points over most of the distribution"]),
    dict(key="usdot2016",
         ref="U.S. Department of Transportation (2016), 'The Value of Travel Time Savings: Departmental Guidance for "
             "Conducting Economic Evaluations, Revision 2 (2016 Update)', memorandum of Sept. 27, 2016",
         url=USDOT_URL,
         read="read from a web.archive.org copy (transportation.gov refuses scripted requests)",
         facts=["local personal travel (commuting included) is valued at 50% of hourly median household income "
                "($56,516 / 2,080 = $27.20 in 2015), $13.60 an hour (Tables 1, 3, 4; text p. 11); plausible range "
                "35%-60%, $9.50-$16.30 (Tables 2 and 5, pp. 14 and 19); business travel 100% of gross compensation, "
                "$25.40"]),
    dict(key="nces_days",
         ref="NCES, Digest of Education Statistics 2017, Table 203.90, 'Average daily attendance (ADA) as a percentage of "
             "total enrollment, school day length, and school year length in public schools, by school level and state: "
             "2007-08 and 2011-12' (source: SASS Public School Data File, 2007-08 and 2011-12)",
         url="https://nces.ed.gov/programs/digest/d17/tables/dt17_203.90.asp",
         read="web page (United States row)",
         facts=["United States average of 179 days in the school year (SE 0.1) and 1,203 hours, 2011-12 (the table reports "
                "school year length for 2011-12 only; 2007-08 has attendance and hours per day, 6.6, only)"]),
    dict(key="nces_sass",
         ref="NCES SASS / NTPS tabulations",
         url="https://nces.ed.gov/surveys/ntps/",
         read="searched; nothing found",
         facts=["no NCES SASS / NTPS table of the share of public school teachers living in the district where they "
                "teach was found: the surveys report teacher mobility between schools and districts, not residence"]),
]


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------

def comb(a: dict[str, float], b: dict[str, float], sign: float = 1.0) -> dict[str, float]:
    """a + sign * b with the standard errors added in quadrature (the district bootstrap and the PUMS
    replicate weights come from independent samples)."""
    return dict(value=a["value"] + sign * b["value"], se=float(np.hypot(a["se"], b["se"])))


def composition_terms(d: pd.DataFrame, rep: Rep, geo: dict[str, Any], cz_w: pd.Series) -> dict[str, Any]:
    """The composition terms CZ by CZ (`cz_terms`) at the main PUMA class cut:
        pub_k12            public K-12 teachers / EDU workers (headline: the FTE are public-school teachers)
        k12                all-sector K-12 teachers / EDU workers
        pub_k12p           public K-12 plus preschool and kindergarten teachers / EDU workers (the FTE count
                           kindergarten teachers, which OCCP 2310-2330 leave out, and drop pre-K)
        direct_pub_k12_emp public K-12 teachers / all employed workers: the model's gap_teach_res_emp measured
                           directly in the PUMS (suburban against city PUMAs of the same CZ)
        direct_k12_emp     all-sector K-12 teachers / all employed workers"""
    pubk = (d.k12 & d.pub).to_numpy()
    allw = np.ones(len(d), bool)
    out = {"pub_k12": cz_terms(d, rep, geo, cz_w, pubk),
           "k12": cz_terms(d, rep, geo, cz_w, d.k12.to_numpy()),
           "pub_k12p": cz_terms(d, rep, geo, cz_w, (d.k12p & d.pub).to_numpy()),
           "direct_pub_k12_emp": cz_terms(d, rep, geo, cz_w, pubk, allw),
           "direct_k12_emp": cz_terms(d, rep, geo, cz_w, d.k12.to_numpy(), allw)}
    for k, v in out.items():
        LOG.info("within-CZ term %-20s %.4f (survey SE %.4f), %d CZs, %.1f%% of the weight", k, v["value"],
                 v["se_survey"], v["n_cz"], 100 * v["weight_share"])
    return out


def edu_with_group_quarters(geo: dict[str, Any], cz_w: pd.Series) -> dict[str, Any]:
    """The within-CZ composition term with GROUP QUARTERS included.  C24010 counts every resident, the
    household PUMS sample used elsewhere drops group quarters; EDU workers in group quarters (student teaching
    and research assistants in dormitories) sit mostly in college towns.  Only the EDU workers are read."""
    files = sorted(PUMS_DIR.glob("*.parquet"))
    cols = ["SERIALNO", "PUMA", "ST", "PWGTP", "OCCP", "COW", "POWSP", "POWPUMA"] + REP
    parts = []
    for f in files:
        x = pd.read_parquet(f, columns=cols)
        parts.append(x[x.OCCP.isin(list(EDU_CODES))])
    d = pd.concat(parts, ignore_index=True)
    d["gq"] = d.SERIALNO.str.contains("GQ")
    d["st"] = d.ST.astype(int)
    d["pid"] = d.st * 100000 + d.PUMA
    attach_classes(d, geo)
    d["edu"] = True
    d["k12"] = d.OCCP.isin(K12_CODES)
    d["pub"] = d.COW.isin([3, 4])
    rep = Rep(d.PWGTP.to_numpy(float), d[REP].to_numpy(np.int32))
    d = d.drop(columns=REP)
    t = cz_terms(d, rep, geo, cz_w, (d.k12 & d.pub).to_numpy())
    w = d.PWGTP.to_numpy(float)
    cls = d.res_cls.to_numpy()
    gq = d.gq.to_numpy()
    t["gq_share_of_edu_workers"] = float(w[gq].sum() / w.sum())
    t["gq_share_by_class"] = {c: float(w[gq & (cls == c)].sum() / w[cls == c].sum()) for c in CLASSES}
    t["n_edu_records"], t["n_gq_records"] = int(len(d)), int(gq.sum())
    LOG.info("EDU workers with group quarters: %.2f%% of them in group quarters (city %.2f%%, suburb %.2f%%); within-CZ term %.4f",
             100 * t["gq_share_of_edu_workers"], 100 * t["gq_share_by_class"]["city"], 100 * t["gq_share_by_class"]["suburb"],
             t["value"])
    return t


def k12_block(dist: dict[str, Any], comp: dict[str, Any], flows: dict[str, Any], terms: dict[str, Any],
              gq: dict[str, Any], sens_cls: list[dict[str, Any]]) -> dict[str, Any]:
    """The K-12 counterparts of the model moments (M2_SPEC_v2: `gap_teach_res_emp`, `net_commute_gap`).

    The ACS district table has the whole education group (EDU), of which public K-12 teachers are the share s
    (PUMS), so Delta log (public K-12 residents) = Delta log (EDU residents) + Delta log s, with Delta log s the
    log of the suburban over the city share inside each CZ of the partition (`cz_terms`):
        gap_teach_res_emp  = gap_edu_res_emp  + Delta log s       (K-12 teachers / employed residents)
        net_commute_gap    = net_commute_gap_edu + Delta log s    (K-12 residents per K-12 FTE)
        gap_k12_res        = gap_edu_res + Delta log s            (K-12 residents)
    HEADLINE: residents counted once, s = public K-12 / EDU workers on the PUMA classes at the main cut, the
    mean of the term over the CZs that have both classes added to the district gap of all 188 CZs.

    Standard error: the district gap and the composition term are bootstrapped TOGETHER over whole CZs on the
    draws of `district_settings.py` (their covariance is in the draws), and the survey (replicate-weight)
    variance of the composition term and the ACS sampling variance of C24010 are added.  The CZ bootstrap
    already holds part of that sampling noise (each CZ's PUMS and ACS estimates are noisy and the bootstrap
    resamples them), so the total is CONSERVATIVE; `se_boot` alone is the lower figure.
    Range: the smallest and largest value over the specifications (overlap rule, sample, class cuts, sectors,
    group quarters, geography of the composition)."""
    counts, w, g = dist["counts"], dist["w"], dist["g"]
    st, V, acs = dist["stats"], dist["variants"], dist["acs_se"]
    wm = lambda x: db.boot_wmean(np.asarray(x, float), w, counts)   # noqa: E731
    tz = terms["pub_k12"]
    comp_b = wm(tz["est"])
    base = dict(net=g.gap_tr_edu_fte.to_numpy(float), emp=g.gap_tr_teach_emp.to_numpy(float),
                res=g.gap_tr_edu_res.to_numpy(float))
    D_main = comp_b["est"]
    if abs(D_main - tz["value"]) > 1e-10:
        raise SystemExit("the CZ-weighted composition term differs between cz_terms and the bootstrap helper")

    def headline(base_key: str, edu_key: str, acs_key: str, label: str) -> dict[str, Any]:
        b = wm(base[base_key])
        draws = b["draws"] + comp_b["draws"]
        se_boot = float(np.std(draws[np.isfinite(draws)], ddof=1))
        se_surv, se_acs = tz["se_survey"], acs[acs_key]
        se = float(np.sqrt(se_boot ** 2 + se_surv ** 2 + se_acs ** 2))
        both = wm(np.where(np.isfinite(tz["est"]), base[base_key] + tz["est"], np.nan))
        v = b["est"] + D_main
        return dict(label=label, value=float(v), se=se, se_boot=se_boot, se_survey=float(se_surv), se_acs=float(se_acs),
                    se_quadrature_of_district_and_survey=float(np.hypot(st[edu_key]["se"], se_surv)),
                    ci95_lo=float(v - 1.96 * se), ci95_hi=float(v + 1.96 * se), n_cz=int(b["n"]),
                    n_cz_composition=int(comp_b["n"]), weight_share_composition=tz["weight_share"],
                    corr_district_vs_composition=float(np.corrcoef(b["draws"][np.isfinite(draws)],
                                                                   comp_b["draws"][np.isfinite(draws)])[0, 1]),
                    district_gap_edu=float(b["est"]), composition_term=float(D_main),
                    per_cz_sum_common_czs=dict(value=float(both["est"]), se_boot=float(np.std(both["draws"][np.isfinite(both["draws"])], ddof=1)),
                                               n_cz=int(both["n"])),
                    _draws=draws)

    H = {"gap_teach_res_emp": headline("emp", "gap_edu_res_emp", "gap_edu_res_emp",
                                       "public K-12 teacher residents / employed residents (gap_teach_res_emp)"),
         "net_commute_gap": headline("net", "net_commute_gap_edu", "net_commute_gap_edu",
                                     "public K-12 teacher residents per K-12 FTE (net_commute_gap)"),
         "gap_k12_res": headline("res", "gap_edu_res", "gap_edu_res", "public K-12 teacher residents (gap_edu_res + composition)")}

    # ---- specifications --------------------------------------------------------------------------------
    spec: list[dict[str, Any]] = []

    def add(key: str, label: str, group: str, net: float | None, emp: float | None, res: float | None, n_cz: int | None = None,
            in_range: bool = True) -> None:
        spec.append(dict(key=key, label=label, group=group, net_commute_gap=net, gap_teach_res_emp=emp, gap_k12_res=res, n_cz=n_cz,
                         in_range=in_range))

    Gh, Eh, Rh = (float(st[k]["value"]) for k in ("net_commute_gap_edu", "gap_edu_res_emp", "gap_edu_res"))
    add("headline", "headline: residents counted once, public K-12 / EDU within CZs, PUMA classes 0.5 / 0.5", "headline",
        Gh + D_main, Eh + D_main, Rh + D_main, st["net_commute_gap_edu"]["n_cz"])
    for key, lab in (("earnings_rule", "earnings rule: every LEA takes its own layer's record (residents counted twice)"),
                     ("no_overlap_sec", "residents once from unified + secondary records"),
                     ("unified", "unified districts only"),
                     ("nyc_dropped", "New York City dropped from both counts (not in the range: the New York CZ keeps its weight "
                                     "with almost no city location left, so its gap is dominated by what remains; `drop_ny_cz` is the clean version)"),
                     ("drop_ny_cz", "New York CZ left out")):
        if key in V:
            add(key, lab, "district sample or overlap rule", float(V[key]["net_commute_gap_edu"]["value"]) + D_main,
                float(V[key]["gap_edu_res_emp"]["value"]) + D_main, float(V[key]["gap_edu_res"]["value"]) + D_main,
                V[key]["net_commute_gap_edu"]["n_cz"], in_range=key != "nyc_dropped")
    for key, lab, dv in (("all_sector", "composition: all-sector K-12 teachers / EDU", terms["k12"]["value"]),
                         ("incl_prek", "composition: public K-12 plus preschool / kindergarten (OCCP 2300) / EDU", terms["pub_k12p"]["value"]),
                         ("national_classes", "composition: national city / suburb PUMA split (not CZ by CZ)",
                          comp["split"]["pub_k12"]["log_ratio"]["value"]),
                         ("powpuma_classes", "composition: POWPUMA classes (county based, national)",
                          -flows["edu_res_vs_pubk12_work"]["composition"]["value"]),
                         ("group_quarters", "composition: group quarters included (as in C24010)", gq["value"])):
        add(key, lab, "composition term", Gh + dv, Eh + dv, Rh + dv)
    for r_ in sens_cls:
        if not r_["main"]:
            dv = r_["within_cz"]["pub_k12"]["value"]
            add(f"cut_{r_['pc_cut']}_{r_['metro_cut']}", f"composition: PUMA class cut pc >= {r_['pc_cut']}, metro >= {r_['metro_cut']}",
                "composition term", Gh + dv, Eh + dv, Rh + dv, r_["within_cz"]["pub_k12"]["n_cz"])
    add("per_cz_sum", "district gap plus composition added CZ by CZ on the CZs that have both", "composition term",
        H["net_commute_gap"]["per_cz_sum_common_czs"]["value"], H["gap_teach_res_emp"]["per_cz_sum_common_czs"]["value"],
        H["gap_k12_res"]["per_cz_sum_common_czs"]["value"], H["net_commute_gap"]["per_cz_sum_common_czs"]["n_cz"])
    for k in ("gap_teach_res_emp", "net_commute_gap", "gap_k12_res"):
        vals = [x[k] for x in spec if x[k] is not None and x["in_range"]]
        H[k]["spec_lo"], H[k]["spec_hi"] = float(min(vals)), float(max(vals))
        H[k]["spec_n"] = len(vals)
    # cross-checks that do not use the district gap
    direct = dict(
        gap_teach_res_emp_within_cz_pums=dict(value=terms["direct_pub_k12_emp"]["value"], se_survey=terms["direct_pub_k12_emp"]["se_survey"],
                                             n_cz=terms["direct_pub_k12_emp"]["n_cz"], weight_share=terms["direct_pub_k12_emp"]["weight_share"],
                                             label="PUMS: public K-12 teachers / employed residents, suburban against city PUMAs of the same CZ, CZ-weighted"),
        gap_teach_res_emp_all_sector_within_cz_pums=dict(value=terms["direct_k12_emp"]["value"], se_survey=terms["direct_k12_emp"]["se_survey"],
                                                         n_cz=terms["direct_k12_emp"]["n_cz"]),
        gap_teach_res_emp_powpuma_classes_national=flows["pub_k12"]["gap_res_emp"],
        net_commute_gap_powpuma_classes_national=flows["pub_k12"]["net_commute_gap"])
    # bootstrap of the direct within-CZ measure on the shared CZ draws (CZ-level uncertainty, no district part)
    dwm = wm(terms["direct_pub_k12_emp"]["est"])
    direct["gap_teach_res_emp_within_cz_pums"]["se_cz_bootstrap"] = float(np.std(dwm["draws"][np.isfinite(dwm["draws"])], ddof=1))
    direct["gap_teach_res_emp_within_cz_pums"]["se"] = float(np.hypot(direct["gap_teach_res_emp_within_cz_pums"]["se_cz_bootstrap"],
                                                                    terms["direct_pub_k12_emp"]["se_survey"]))
    return dict(headline=H, specs=spec, direct=direct, composition=dict(
        pub_k12=dict(value=tz["value"], se_survey=tz["se_survey"], se_cz_bootstrap=float(np.std(comp_b["draws"][np.isfinite(comp_b["draws"])], ddof=1)),
                     n_cz=tz["n_cz"], weight_share=tz["weight_share"], n_sample=tz["n_sample"]),
        k12=dict(value=terms["k12"]["value"], se_survey=terms["k12"]["se_survey"], n_cz=terms["k12"]["n_cz"]),
        pub_k12p=dict(value=terms["pub_k12p"]["value"], se_survey=terms["pub_k12p"]["se_survey"], n_cz=terms["pub_k12p"]["n_cz"]),
        group_quarters=dict(value=gq["value"], se_survey=gq["se_survey"], n_cz=gq["n_cz"], gq_share_of_edu_workers=gq["gq_share_of_edu_workers"],
                            gq_share_by_class=gq["gq_share_by_class"], n_gq_records=gq["n_gq_records"], n_edu_records=gq["n_edu_records"])),
        _draws={k: v["_draws"] for k, v in H.items()})


def strip(o: Any) -> Any:
    """Drop private keys, arrays and frames from a result tree (for JSON)."""
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if not (isinstance(k, str) and k.startswith("_")) and k != "draws"}
    if isinstance(o, (list, tuple)):
        return [strip(v) for v in o]
    if isinstance(o, pd.DataFrame):
        return o.round(6).to_dict()
    return o


CORR_KEYS_TOML = CORR_TARGETS + ["k12_gap_teach_res_emp", "k12_net_commute_gap", "gap_edu_res", "gap_fte"]
CORR_KEYS_MD = CORR_TARGETS + ["k12_gap_teach_res_emp", "k12_net_commute_gap", "k12_gap_k12_res", "gap_edu_res_emp",
                               "gap_edu_res", "gap_fte", "net_commute_gap_edu"]


def write_toml(path: Path, R: dict[str, Any]) -> None:
    D, m, K, RK = R["districts"], R["meta"], R["k12"], R["rk"]
    tv = ds.toml_val
    L: list[str] = []
    a = L.append
    a("# teacher_residence_targets.toml -- residence-based against workplace-based teachers, for the M2 split.")
    a(f"# Generated {m['generated']} by data/spatial/teacher_residence.py; see teacher_residence.md.")
    a("#")
    a("# WHICH VALUE MAPS TO WHICH MODEL MOMENT (M2_SPEC_v2.md; all gaps are log(suburb) - log(city), location 2 = suburb)")
    a("#   gap_teach_res_emp  (model: dlog teachers residing / market workers residing)  <-  [k12.gap_teach_res_emp]")
    a("#   net_commute_gap    (model: dlog n_teach_res - dlog n_teach_work)               <-  [k12.net_commute_gap]")
    a("#   dlog n_teach_res   (K-12 teacher residents)                                    <-  [k12.gap_k12_res]")
    a("#   tau_k  (proportional commuting time cost, income x (1 - tau_k); spec v2)       <-  [tau_k.teachers]; [tau_k.other_workers] for non-teachers")
    a("#   r_k    (commuting goods cost in units of mean consumption; spec v1)            <-  [r_k.per_adult] (the model's agents are adults); [r_k.per_cu] is per CE consumer unit")
    a("# The [k12.*] values are K-12 (public-school) teachers: the district gap of the whole education, training and library group")
    a("# (EDU, SOC 25-0000) plus a PUMS composition term (approximate; see the specification range spec_lo / spec_hi).")
    a("# The [locale.*] tables are the EDU district gaps themselves: EDU PROXIES (gap_edu_res_emp, net_commute_gap_edu) that include")
    a("# postsecondary teachers, aides, librarians and private-school staff. They are NOT the model counterparts and have a")
    a("# different (negative) sign for the net gap; do not target them against a K-12 teacher object.")
    a("# value / se / n_cz: pupil-weighted mean over CZs of the within-CZ gap, whole-CZ bootstrap SE (B = %d, seed %d; the draws" %
      (m["B"], m["seed"]))
    a("# of district_settings.py, so the correlations are valid). Residents are counted once (unified + elementary ACS records).")
    a("")
    a("[meta]")
    for k in ("generated", "script", "year", "B", "seed"):
        a(f"{k} = {tv(m[k])}")
    a(f"n_cz = {tv(D['coverage']['n_cz'])}")
    a("")

    # ---------------------------------------------------------------- K-12 counterparts of the model moments
    a("[k12]")
    a(f"description = {tv('K-12 (public-school) teachers: district EDU gap (residents counted once) plus the PUMS composition term, joint CZ bootstrap on the district_settings draws; approximate')}")
    a(f"n_cz = {tv(D['coverage']['n_cz'])}")
    a(f"corr_keys = {tv(CORR_KEYS_TOML)}")
    a(f"corr = {tv(D['corr_matrix_toml'])}")
    a("# corr: bootstrap correlation (district and composition draws together) with the four adopted targets")
    a("")
    status = {"gap_teach_res_emp": "inactive target / validation candidate: model gap_teach_res_emp",
              "net_commute_gap": "inactive target / validation candidate: model net_commute_gap",
              "gap_k12_res": "memo: K-12 teacher residents gap (model: dlog n_teach_res)"}
    for key in ("gap_teach_res_emp", "net_commute_gap", "gap_k12_res"):
        h = K["headline"][key]
        a(f"[k12.{key}]")
        a(f"label = {tv(h['label'])}")
        a(f"status = {tv(status[key])}")
        a(f"value = {tv(h['value'])}")
        a(f"se = {tv(h['se'])}")
        a("# se = sqrt(se_boot^2 + se_survey^2 + se_acs^2): conservative (the CZ bootstrap already holds part of the survey noise)")
        a(f"se_boot = {tv(h['se_boot'])}")
        a(f"se_survey = {tv(h['se_survey'])}")
        a(f"se_acs = {tv(h['se_acs'])}")
        a(f"ci95 = {tv([h['ci95_lo'], h['ci95_hi']])}")
        a(f"spec_lo = {tv(h['spec_lo'])}")
        a(f"spec_hi = {tv(h['spec_hi'])}")
        a(f"n_spec = {tv(h['spec_n'])}")
        a(f"n_cz = {tv(h['n_cz'])}")
        a(f"n_cz_composition = {tv(h['n_cz_composition'])}")
        a(f"weight_share_composition = {tv(h['weight_share_composition'])}")
        a(f"corr_district_vs_composition = {tv(h['corr_district_vs_composition'])}")
        a(f"district_gap_edu = {tv(h['district_gap_edu'])}")
        a(f"composition_term = {tv(h['composition_term'])}")
        a("")
    cp = K["composition"]["pub_k12"]
    a("[k12.composition]")
    a("# Delta log s, s = public K-12 teachers / EDU workers, suburban against city PUMAs of the same CZ (ACS 2018 1-year PUMS)")
    a(f"delta_log_s = {tv(cp['value'])}")
    a(f"se_survey = {tv(cp['se_survey'])}")
    a(f"se_cz_bootstrap = {tv(cp['se_cz_bootstrap'])}")
    a(f"n_cz = {tv(cp['n_cz'])}")
    a(f"weight_share = {tv(cp['weight_share'])}")
    a("")
    sp = K["specs"]
    a("[k12.specifications]")
    a("# the value of each K-12 moment under every specification (order of `names`; nan = not defined for that moment)")
    a(f"names = {tv([x['key'] for x in sp])}")
    for key in ("gap_teach_res_emp", "net_commute_gap", "gap_k12_res"):
        a(f"{key} = {tv([float('nan') if x[key] is None else x[key] for x in sp])}")
    a(f"in_range = {tv([bool(x['in_range']) for x in sp])}")
    a(f"labels = {tv([x['label'] for x in sp])}")
    a("")
    dr = K["direct"]
    a("[k12.cross_checks]")
    a("# computed without the district table")
    d1 = dr["gap_teach_res_emp_within_cz_pums"]
    a(f"gap_teach_res_emp_within_cz_pums = {tv(d1['value'])}")
    a(f"gap_teach_res_emp_within_cz_pums_se = {tv(d1['se'])}")
    a(f"gap_teach_res_emp_powpuma_classes_national = {tv(dr['gap_teach_res_emp_powpuma_classes_national']['value'])}")
    a(f"gap_teach_res_emp_powpuma_classes_national_se = {tv(dr['gap_teach_res_emp_powpuma_classes_national']['se'])}")
    a(f"net_commute_gap_powpuma_classes_national = {tv(dr['net_commute_gap_powpuma_classes_national']['value'])}")
    a(f"net_commute_gap_powpuma_classes_national_se = {tv(dr['net_commute_gap_powpuma_classes_national']['se'])}")
    a("")

    # ---------------------------------------------------------------- commuting cost: tau_k and r_k
    a("[tau_k]")
    a("# proportional time cost of working in the other location: income x (1 - tau_k) when r != w (M2_SPEC_v2.md)")
    a("# tau_k = value-of-time share of the own wage x 2 x extra one-way hours x days / annual hours worked")
    a(f"formula = {tv('tau_k = vot_share * 2 * extra_one_way_minutes / 60 * days_per_year / annual_hours')}")
    for grp, lab in (("teachers", "K-12 teachers (185 days)"), ("other_workers", "other workers (250 days)")):
        t = RK["tau_k"][grp]
        a(f"[tau_k.{grp}]")
        a(f"description = {tv(lab)}")
        a(f"central = {tv(t['central'])}")
        a(f"low = {tv(t['low'])}")
        a(f"high = {tv(t['high'])}")
        a(f"central_se_from_extra_minutes = {tv(t['se_from_extra_minutes'])}")
        a(f"central_with_2080_hours = {tv(t['central_with_2080_hours'])}")
        a(f"extra_one_way_minutes = {tv(t['extra_one_way_minutes'])}")
        a(f"days_per_year = {tv(t['days'])}")
        a(f"annual_hours = {tv(t['annual_hours'])}")
        a(f"value_of_time_share_of_own_wage = {tv(t['value_of_time_share_of_own_wage'])}")
        a("")
    a("[r_k]")
    a("# goods cost of the extra commuting in units of mean consumption Cbar (M2_SPEC.md): annual time cost / Cbar")
    a("# the model's agents are single adults and the cost is per worker: use [r_k.per_adult]; [r_k.per_cu] divides by CE")
    a("# consumption per consumer unit (2.5 persons, 1.9 adults), the unit of r_p, and is NOT the per-adult number")
    a(f"formula = {tv('r_k = 2 * extra_one_way_minutes / 60 * days_per_year * value_of_hour / cbar')}")
    a(f"extra_one_way_minutes = {tv([RK['delta_min']['low'], RK['delta_min']['central'], RK['delta_min']['high']])}")
    a(f"days_per_year = {tv([RK['days']['low'], RK['days']['central'], RK['days']['high']])}")
    a(f"value_of_hour_2018_usd = {tv([RK['vtt_2018']['low'], RK['vtt_2018']['central'], RK['vtt_2018']['high']])}")
    a(f"annual_time_cost_central_usd = {tv(RK['annual_time_cost_central_usd'])}")
    a("")
    a("[r_k.per_adult]")
    a(f"cbar_usd = {tv(RK['cbar_adult'])}")
    a(f"adults_per_consumer_unit = {tv(RK['adults_per_cu'])}")
    a(f"central = {tv(RK['central_per_adult'])}")
    a(f"central_se_from_extra_minutes = {tv(RK['central_se_from_delta_per_adult'])}")
    a(f"low = {tv(RK['low_per_adult'])}")
    a(f"high = {tv(RK['high_per_adult'])}")
    a("")
    a("[r_k.per_cu]")
    a(f"cbar_usd = {tv(RK['cbar_cu'])}")
    a(f"central = {tv(RK['central_per_cu'])}")
    a(f"central_se_from_extra_minutes = {tv(RK['central_se_from_delta_per_cu'])}")
    a(f"low = {tv(RK['low_per_cu'])}")
    a(f"high = {tv(RK['high_per_cu'])}")
    a("")

    # ---------------------------------------------------------------- EDU proxies (district gaps)
    keys = D["corr_keys"]
    a("[locale]")
    a(f"description = {tv('EDU (SOC 25-0000) PROXIES, not the model counterparts: locale split (city vs suburb), the 188 CZs of the tracked moments; residents counted once (unified + elementary ACS records); matched sample of districts with a C24010 record and an audited K-12 FTE and pupil count')}")
    a(f"n_cz = {tv(D['coverage']['n_cz'])}")
    a(f"corr_keys = {tv(keys)}")
    a(f"corr = {tv(D['corr_matrix'])}")
    a("")
    for key in ("gap_edu_res_emp", "net_commute_gap_edu", "gap_edu_res", "gap_fte", "gap_edu_per_pupil", "gap_emp_res"):
        s = D["stats"][key]
        a(f"[locale.{key}]")
        a(f"label = {tv(s['label'])}")
        a(f"value = {tv(s['value'])}")
        a(f"se = {tv(s['se'])}")
        if "se_acs_sampling" in s:
            a(f"se_acs_sampling = {tv(s['se_acs_sampling'])}")
        a(f"ci95 = {tv([s['ci95_lo'], s['ci95_hi']])}")
        a(f"n_cz = {tv(s['n_cz'])}")
        a("")
    for sub, lab in (("edu_per_fte", "EDU residents per K-12 FTE"), ("edu_share_emp", "EDU residents per employed resident")):
        lv = D["levels"][sub]
        a(f"[locale.level_{sub}]")
        a(f"label = {tv(lab + ', CZ-weighted mean of the location level (city = loc1, suburb = loc2)')}")
        a(f"loc1 = {tv(lv['loc1']['value'])}")
        a(f"loc1_se = {tv(lv['loc1']['se'])}")
        a(f"loc2 = {tv(lv['loc2']['value'])}")
        a(f"loc2_se = {tv(lv['loc2']['se'])}")
        a("")
    for var, lab in (("earnings_rule", "earnings rule: every LEA takes its own layer's record (residents of overlapping elementary and secondary districts counted twice)"),
                     ("no_overlap_sec", "residents once from unified + secondary records (elementary residents dropped where a secondary LEA is present)"),
                     ("unified", "unified districts only"),
                     ("nyc_dropped", "New York City dropped from both the resident and the FTE / pupil counts"),
                     ("drop_ny_cz", "New York CZ left out")):
        if var not in D["variants"]:
            continue
        a(f"[locale_{var}]")
        a(f"description = {tv(lab)}")
        a("")
        for key in ("gap_edu_res_emp", "net_commute_gap_edu", "gap_edu_res"):
            s = D["variants"][var][key]
            a(f"[locale_{var}.{key}]")
            a(f"value = {tv(s['value'])}")
            a(f"se = {tv(s['se'])}")
            a(f"n_cz = {tv(s['n_cz'])}")
            a("")
    C, F = R["composition"], R["flows"]
    a("[composition]")
    a("# ACS 2018 1-year PUMS, employed 16+ in households; replicate-weight SEs. s = share of EDU workers who are K-12 teachers")
    a(f"s_k12 = {tv(C['share_k12']['value'])}")
    a(f"s_k12_se = {tv(C['share_k12']['se'])}")
    a(f"s_public_k12 = {tv(C['share_pub_k12']['value'])}")
    a(f"s_public_k12_se = {tv(C['share_pub_k12']['se'])}")
    for c in ("city", "suburb"):
        a(f"s_public_k12_{c} = {tv(C['by_class'][c]['pub_k12']['value'])}")
        a(f"s_public_k12_{c}_se = {tv(C['by_class'][c]['pub_k12']['se'])}")
    lr = C["split"]["pub_k12"]["log_ratio"]
    a(f"delta_log_s_public_k12_national_classes = {tv(lr['value'])}")
    a(f"delta_log_s_public_k12_national_classes_se = {tv(lr['se'])}")
    a(f"delta_log_s_public_k12_within_cz = {tv(cp['value'])}")
    a(f"delta_log_s_public_k12_within_cz_se = {tv(cp['se_survey'])}")
    a(f"delta_log_s_public_k12_within_cz_n_cz = {tv(cp['n_cz'])}")
    a(f"delta_log_s_public_k12_within_cz_weight_share = {tv(cp['weight_share'])}")
    dd = C["dispersion"]["pub_k12"]["sd_w_net_of_sampling"]
    a(f"sd_s_public_k12_across_pumas = {tv(dd['value'])}")
    a(f"sd_s_public_k12_across_pumas_se = {tv(dd['se'])}")
    a("")
    a("[pums_k12_flows]")
    a("# K-12 teachers by residence and by workplace on the same POWPUMA classes (national, 2018 1-year PUMS)")
    for k in ("k12", "pub_k12"):
        a(f"[pums_k12_flows.{k}]")
        a(f"net_commute_gap = {tv(F[k]['net_commute_gap']['value'])}")
        a(f"net_commute_gap_se = {tv(F[k]['net_commute_gap']['se'])}")
        a(f"gap_res_emp = {tv(F[k]['gap_res_emp']['value'])}")
        a(f"gap_res_emp_se = {tv(F[k]['gap_res_emp']['se'])}")
        a(f"share_work_in_other_class = {tv(F[k]['share_cross_class']['value'])}")
        a(f"share_work_in_other_class_se = {tv(F[k]['share_cross_class']['se'])}")
        a(f"share_metro_city_suburb_cross = {tv(F[k]['share_metro_cross']['value'])}")
        a(f"share_metro_city_suburb_cross_se = {tv(F[k]['share_metro_cross']['se'])}")
        a("")
    Kc = R["commute"]["stats"]
    a("[commute]")
    a("# K-12 teachers (OCCP 2310, 2320, 2330) and other employed workers 25-64, one-way commute of commuters")
    for lab, g in (("teachers", "T"), ("other_workers", "O")):
        a(f"[commute.{lab}]")
        a(f"mean_minutes = {tv(Kc[g]['mean']['value'])}")
        a(f"mean_minutes_se = {tv(Kc[g]['mean']['se'])}")
        a(f"share_over_30 = {tv(Kc[g]['gt30']['value'])}")
        a(f"share_over_30_se = {tv(Kc[g]['gt30']['se'])}")
        a(f"share_30_or_more = {tv(Kc[g]['ge30']['value'])}")
        a(f"share_30_or_more_se = {tv(Kc[g]['ge30']['se'])}")
        a(f"share_exactly_30 = {tv(Kc[g]['eq30']['value'])}")
        a(f"share_outside_residence_powpuma = {tv(Kc[g]['out']['value'])}")
        a(f"share_outside_residence_powpuma_se = {tv(Kc[g]['out']['se'])}")
        a("")
    path.write_text("\n".join(L) + "\n")


def fx(v: Any, d: int = 3) -> str:
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def fs(s: dict[str, Any], d: int = 3) -> str:
    return f"{s['value']:.{d}f} ({s['se']:.{d}f})"


def fci(s: dict[str, Any], d: int = 3) -> str:
    return f"[{s['ci95_lo']:.{d}f}, {s['ci95_hi']:.{d}f}]"


def pc(v: float, d: int = 1) -> str:
    return f"{100 * v:.{d}f}%"


def pcs(s: dict[str, Any], d: int = 1) -> str:
    return f"{100 * s['value']:.{d}f}% ({100 * s['se']:.{d}f})"


def write_md(path: Path, R: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    m, D, C, K, F = R["meta"], R["districts"], R["composition"], R["commute"], R["flows"]
    W, RK, KB, GEO = R["wages"], R["rk"], R["k12"], R["geo"]
    S, cov = D["stats"], D["coverage"]
    lv, V = D["levels"], D["variants"]
    nc, gt, ge, gf = S["net_commute_gap_edu"], S["gap_edu_res_emp"], S["gap_edu_res"], S["gap_fte"]
    H = KB["headline"]
    h_t, h_n, h_r = H["gap_teach_res_emp"], H["net_commute_gap"], H["gap_k12_res"]
    cp = KB["composition"]["pub_k12"]
    sh = C["by_class"]
    kt, ko = K["stats"]["T"], K["stats"]["O"]
    fe = K["fe"]
    TK = RK["tau_k"]
    repro_txt = ", ".join("%s %.4f" % (r["moment"], r["value"]) for r in D["reproduction"])
    nm1, nm2 = cov["single_count"]["n_matched_by_location"][1], cov["single_count"]["n_matched_by_location"][2]
    cuts = [x for x in KB["specs"] if x["key"].startswith("cut_")]
    ow = W["k12"]
    er = V["earnings_rule"]
    vmax = max((abs(V[k]["net_commute_gap_edu"]["value"] - nc["value"]), k) for k in ("earnings_rule", "no_overlap_sec", "unified", "nyc_dropped", "drop_ny_cz", "unmatched"))
    dr = KB["direct"]
    d1 = dr["gap_teach_res_emp_within_cz_pums"]
    ci = lambda h: f"[{h['ci95_lo']:.2f}, {h['ci95_hi']:.2f}]"   # noqa: E731
    rg = lambda h: f"[{h['spec_lo']:.2f}, {h['spec_hi']:.2f}]"   # noqa: E731
    fk = lambda h, d=3: f"{h['value']:.{d}f} ({h['se']:.{d}f})"   # noqa: E731

    a("# Teacher residence versus workplace: evidence for the M2 split\n")
    a(f"Generated {m['generated']} by `teacher_residence.py` (DESIGN.md section M2, M2_SPEC_v2.md). Residence-based teachers: ACS "
      f"5-year 2014-18 school-district table C24010; workplace-based teachers: the audited K-12 FTE of "
      f"`district_settings.py`; locale partition, {cov['n_cz']} CZs, whole-CZ bootstrap B = {m['B']}, seed "
      f"{m['seed']} (the draws of `district_settings.py`). Composition and commuting: ACS 2018 1-year PUMS "
      f"(replicate-weight standard errors) with the Census 2010 PUMA equivalency files and the OMB September-2018 "
      f"principal-city list. EDU is the occupation group 'Educational instruction, and library occupations' (SOC 25-0000, education, training "
      f"and library). Sign convention: suburb minus city, in logs. **Names.** The model's `gap_teach_res_emp` and `net_commute_gap` are about "
      f"K-12 teachers; the district table can only give the EDU group, so the district gaps are EDU proxies (`gap_edu_res_emp`, "
      f"`net_commute_gap_edu`) and the K-12 counterparts are built in section 3 (`[k12]` of the TOML). Residents are counted once.\n")

    a("## Summary\n")
    a(f"- **K-12 counterparts of the model moments ({h_n['n_cz']} CZs; approximate).** The ACS district table holds the whole EDU group, of which "
      f"public K-12 teachers are the share s (PUMS). Adding the log ratio of s in suburban to city PUMAs of the same CZ "
      f"(Delta log s = {cp['value']:.3f}, survey SE {cp['se_survey']:.3f}, {cp['n_cz']} CZs with both classes, {pc(cp['weight_share'], 0)} of the weight) "
      f"to the EDU district gaps gives **`gap_teach_res_emp` = {fk(h_t)}** (public K-12 teachers per employed resident; range over "
      f"specifications {rg(h_t)}) and **`net_commute_gap` = {fk(h_n)}** (K-12 residents per K-12 FTE; range {rg(h_n)}); the K-12 residents' gap is "
      f"{fk(h_r)} against a K-12 FTE gap of {fs(gf)}. The SE adds the CZ-bootstrap variance of the district gap and the composition term "
      f"(together, correlation {h_n['corr_district_vs_composition']:.2f}), the PUMS survey variance and the ACS sampling variance; it is conservative. "
      f"Cross-checks without the district table: public K-12 teachers per employed resident, suburban against city PUMAs of the same CZ, "
      f"{fs(d1)}; national residence against workplace of public K-12 teachers on POWPUMA classes, {fs(F['pub_k12']['net_commute_gap'])} "
      f"(by residence {fs(F['pub_k12']['gap_res'])}, by workplace {fs(F['pub_k12']['gap_work'])}). Reading: more public K-12 teachers per resident worker live in "
      f"the suburb (positive in every specification, {rg(h_t)}); the net commuting gap is small and positive in the headline (more teachers live in the suburb than "
      f"work there: a net flow of teachers from suburban homes to city schools) but its sign is not secure across specifications {rg(h_n)}.")
    a(f"- **EDU proxies (not the model counterparts).** `gap_edu_res_emp` (EDU residents / employed residents) {fs(gt)}; `gap_edu_res` {fs(ge)}; "
      f"`gap_fte` (K-12 FTE) {fs(gf)}; **`net_commute_gap_edu` = `gap_edu_res` - `gap_fte` = {fs(nc)}** (95% interval {fci(nc)}), the log gap of "
      f"EDU residents per K-12 FTE: {fx(lv['edu_per_fte']['loc1']['value'], 2)} in the city, {fx(lv['edu_per_fte']['loc2']['value'], 2)} in the "
      f"suburb. Its negative sign is composition: public K-12 teachers are {pcs(sh['city']['pub_k12'])} of EDU workers in principal-city PUMAs and "
      f"{pcs(sh['suburb']['pub_k12'])} in suburban PUMAs (postsecondary teachers {pcs(sh['city']['occ']['2205'])} and {pcs(sh['suburb']['occ']['2205'])}). "
      f"The four adopted targets reproduce ({repro_txt}); the bootstrap correlations of the new moments with them are in section 2 "
      f"(`net_commute_gap_edu` with `gap_teachers_pp` {D['corr_matrix'][CORR_KEYS_MD.index('net_commute_gap_edu')][2]:.2f}, `gap_salary` "
      f"{D['corr_matrix'][CORR_KEYS_MD.index('net_commute_gap_edu')][1]:.2f}).")
    a(f"- **Residents are counted once.** The earnings rule (each LEA takes its own layer's record) counts residents of an elementary district inside "
      f"a secondary district twice when both sit in one location (employed residents overcounted by {100 * (cov['overlap_factor_employed_residents'] - 1):.1f}% "
      f"in the matched sample, more in the suburbs). The headline uses the unified and elementary records, which tile the states, and keeps the "
      f"secondary LEAs' FTE and pupils; the earnings rule would give `net_commute_gap_edu` {fs(er['net_commute_gap_edu'])} and `gap_edu_res` "
      f"{fs(er['gap_edu_res'])} (bias {er['net_commute_gap_edu']['value'] - nc['value']:+.3f}), residents from unified + secondary records "
      f"{fs(V['no_overlap_sec']['net_commute_gap_edu'])}; the ratio `gap_edu_res_emp` is insensitive ({fs(er['gap_edu_res_emp'])} under the earnings rule).")
    a(f"- **Commuting (PUMS, commuters 25-64).** K-12 teachers commute {fs(kt['mean'], 1)} minutes one way against {fs(ko['mean'], 1)} "
      f"for other workers; {pcs(kt['gt30'])} commute more than 30 minutes against {pcs(ko['gt30'])} ({pcs(kt['ge30'])} against {pcs(ko['ge30'])} "
      f"for 30 minutes or more: reported times heap at exactly 30, {pcs(kt['eq30'])} and {pcs(ko['eq30'])}); {pcs(kt['out'])} work outside the "
      f"POWPUMA they live in against {pcs(ko['out'])}, and {pcs(F['pub_k12']['share_metro_cross'])} of public K-12 teachers who live and work in a "
      f"city or suburban POWPUMA cross between the two classes. Working in another POWPUMA adds **{fe['T@metro']['coef']:.1f} minutes** (SE "
      f"{fe['T@metro']['se']:.1f}) one way for teachers in metropolitan areas (residence-PUMA fixed effects) and working in a POWPUMA of the other class, "
      f"city or suburb, **{fe['T@metro:type']['coef']:.1f}** ({fe['T@metro:type']['se']:.1f}; residence and workplace classified on the same POWPUMA "
      f"geography). Both measure crossings of a county line, so they may overstate the increment within a county; Nashville teachers who live out of the "
      f"county commute {NASHVILLE_EXTRA_MIN} minutes longer (Santelli and Grissom 2024).")
    tt, to = TK["teachers"], TK["other_workers"]
    a(f"- **Model units.** Extra commute {RK['delta_min']['central']:.1f} minutes one way (range {RK['delta_min']['low']:.1f}-{RK['delta_min']['high']:.1f}), "
      f"{RK['days']['central']:.0f} days, value of an hour ${RK['vtt_2018']['central']:.2f} (USDOT: 50% of hourly household income; range 35%-60%): an annual "
      f"time cost of ${RK['annual_time_cost_central_usd']:,.0f} ({RK['annual_hours_central']:.0f} hours). **r_k per adult = {RK['central_per_adult']:.4f} "
      f"(range {RK['low_per_adult']:.4f} to {RK['high_per_adult']:.4f})** in units of mean consumption per adult (CE 2018 ${RK['cbar_cu']:,.2f} per consumer unit "
      f"/ {RK['adults_per_cu']:.1f} adults = ${RK['cbar_adult']:,.0f}); per CE consumer unit (the unit of r_p) it is {RK['central_per_cu']:.4f} "
      f"({RK['low_per_cu']:.4f} to {RK['high_per_cu']:.4f}). **tau_k** (share of income lost, M2_SPEC_v2) is **{tt['central']:.3f}** for teachers "
      f"(range {tt['low']:.3f}-{tt['high']:.3f}; 185 days) and {to['central']:.3f} ({to['low']:.3f}-{to['high']:.3f}) for other workers (250 days, "
      f"{to['extra_one_way_minutes'][1]:.1f} extra minutes). Time cost only (section 5).\n")

    a("## 1. Data, definitions and units\n")
    a("- **Residents (a).** ACS 5-year 2014-18 table C24010, 'Sex by Occupation for the Civilian Employed Population 16 Years and Over' "
      "(`teacher_residence.py fetch`, school-district layers unified, elementary and secondary, one state per request as in "
      "`acs_earnings.py`): C24010_001 total, C24010_014 (men) + C24010_050 (women) 'Educational instruction, and library occupations' "
      "(the 2018 vintage's label for SOC 25-0000, Education, Training, and Library; vintages through 2017 call it 'Education, training, and library "
      "occupations'). Annotation codes (negative) are masked (none occur); 'Remainder of <state>' rows (code 99999) are dropped; one record per "
      f"`leaid`, as `district_bootstrap.load_income_columns` does for B20005. {cov['n_by_layer'].get('unified', 0):,} unified, "
      f"{cov['n_by_layer'].get('elementary', 0)} elementary and {cov['n_by_layer'].get('secondary', 0)} secondary LEAs of the partition have their own "
      f"layer's record; {cov['n_by_layer'].get('none', 0)} have none. **Residents are counted once** (headline): the unified and elementary records "
      "tile each state (section 2) and the secondary layer overlays them, so the secondary LEAs contribute their FTE and pupils but no residents; "
      "the earnings rule of the B20005 moments (every LEA its own record) is the variant `earnings_rule`.")
    a("- **Workplace (b).** `teachers_k12` (CCD teacher FTE minus pre-K teacher FTE) and `pupils_k12` (pupils minus pre-K enrollment) "
      "of the audited panel of `district_settings.py`: the FTE variant 11 behind `gap_teachers_pp`; New York City is one district "
      "(the collapsed row carries the Chancellor's id 3620580, under which the ACS files New York City). Regular districts only "
      "(CCD agency types 1-2), so charter LEAs are not in the FTE.")
    a("- **Matched sample.** Every district moment uses the districts with a C24010 EDU count, a C24010 employed total, an audited K-12 FTE "
      f"and K-12 pupils (and, for the secondary LEAs of the headline, an FTE and pupils): {cov['single_count']['n_matched']:,} of {cov['n_districts']:,} partition districts, "
      f"{pc(cov['single_count']['share_pupils_matched'])} of partition pupils, so that `net_commute_gap_edu` is exactly the log gap of EDU residents per FTE on one set of districts. "
      "Each statistic on its own maximal sample is the `unmatched` variant in section 2.")
    a("- **Moments.** As the audited targets: within each CZ, location totals (or ratios of totals) summed over the districts of the "
      "location, then log(suburb / city); pupil-weighted mean over CZs (weights = the CZ's partitioned pupils); whole-CZ bootstrap with "
      f"`district_bootstrap.draw_counts` (B = {m['B']}, seed {m['seed']}), SE = SD of the draws, interval = 2.5-97.5 percentiles. "
      "The draws are the ones `district_settings.py` uses, so the correlations in section 2 can weight the new moments with the four "
      "existing targets.")
    a("- **PUMS (d, e).** ACS 2018 1-year PUMS through the Census API, one state at a time (`acs_occupations.py`'s `Api`): employed "
      "civilians (ESR 1-2) 16+ in households, with OCCP, SOCP, INDP, COW, PUMA, POWPUMA, POWSP, JWMNP, JWTR, AGEP, SEX, SCHL, wage "
      f"variables and the 80 replicate weights; {R['pums']['n_households']:,} records, {R['pums']['weighted_millions']:.1f} million "
      "weighted. Standard errors are the Census successive-difference replication, SE = sqrt(4/80 sum_r (x_r - x)^2), applied to every "
      "statistic including differences and regressions. K-12 teachers are OCCP 2310, 2320 and 2330 (elementary and middle, secondary and "
      "special-education teachers; SOC 25-2020, 25-2030, 25-2050); adding preschool and kindergarten teachers (2300) is a sensitivity (the "
      "PUMS cannot separate kindergarten from preschool, and the FTE count kindergarten teachers). 'Public' is a local or state government employer (COW 3-4). Pipeline check: "
      f"the national mean travel time to work of workers 16+ who do not work at home is {R['pums']['national_check']['mean_travel_time_all_persons']:.2f} "
      f"minutes (published: 27.1) and {100 * R['pums']['national_check']['work_at_home_share_all_persons']:.1f}% work at home (published: 5.3%), "
      "all persons, ACS 2018 1-year (Census Bureau, commuting guidance 'Mean travel time to work', "
      "www2.census.gov/programs-surveys/commuting/guidance/acs-1yr/Mean-travel-time.pdf).")
    a(f"- **PUMA classes.** 2010 PUMAs from the Census PUMA equivalency files (`PUMSEQ10_<ST>.txt`: 2010 population of each place in each "
      f"PUMA) and the OMB September-2018 list of principal cities of *metropolitan* statistical areas (list 2; the NCES city locale is "
      f"territory in a principal city of an urbanized area). pc_share = population in principal-city places / PUMA population; "
      f"**city** if pc_share >= {GEO['pc_cut']:.2f}; **suburb** if at least {GEO['metro_cut']:.0%} of the population lives in "
      f"counties of metropolitan areas (list 1) and pc_share is below the cut; **other** otherwise. {GEO['n_pc_found']} of the "
      f"{GEO['n_pc_omb']} principal cities of metropolitan areas outside Puerto Rico are found in the files "
      f"({GEO['n_pc_unmatched']} unmatched: {'; '.join(GEO['unmatched_examples']) or 'none'}). {GEO['n_puma_by_class'].get('city', 0)} city, "
      f"{GEO['n_puma_by_class'].get('suburb', 0)} suburban and {GEO['n_puma_by_class'].get('other', 0)} other PUMAs hold "
      f"{GEO['pop_by_class']['city'] / 1e6:.0f}, {GEO['pop_by_class']['suburb'] / 1e6:.0f} and {GEO['pop_by_class']['other'] / 1e6:.0f} million "
      f"people (2010). POWPUMAs (county based; the PUMA to POWPUMA conversion of the 2018 PUMS code lists) take the population-weighted "
      f"shares of the PUMAs they aggregate: {GEO['n_powpuma_by_class'].get('city', 0)} city, {GEO['n_powpuma_by_class'].get('suburb', 0)} "
      f"suburban and {GEO['n_powpuma_by_class'].get('other', 0)} other. Spot checks: the PUMAs of New York City, Chicago, Milwaukee, "
      "Hartford and Boston are all city; those of Nassau, Suffolk and Westchester (NY), Lake (IL) and Oakland (MI) counties are suburb "
      "(White Plains, 35% principal city, is suburb); DuPage County's six PUMAs are suburb and a seventh, which holds Naperville, a "
      "principal city, is city.")
    a("- **Units.** r_k is a goods cost per period in units of mean consumption: the annual time cost of the extra commute divided by "
      f"CE 2018 mean annual expenditures, ${RK['cbar_cu']:,.2f} per consumer unit (BLS Table R-1, the denominator of `housing_costs.py`; "
      f"{CE_PERSONS} persons, {CE_CHILDREN} children under 18, so {RK['adults_per_cu']:.1f} adults, BLS Table 1400) or, because the model's agents "
      f"are adults and the cost is per worker, per adult, ${RK['cbar_adult']:,.0f}. tau_k is the same cost over the worker's annual earnings. The value "
      f"of an hour is the USDOT local-personal rate (50% of hourly median household income, ${VTT_2015['central']:.2f} in 2015 dollars) uprated by the CPI-U ratio "
      f"{CPI_U[2018]:.3f} / {CPI_U[2015]:.3f} = {RK['cpi_ratio_2018_over_2015']:.4f} to 2018 dollars.\n")

    a("## 2. Districts: residents against workplace (a)-(c)\n")
    sc = cov["single_count"]
    a(f"Locale partition (city = NCES locale 11-13, suburb = 21-23; districts of >= 500 pupils; each location >= 10% of the CZ's "
      f"partitioned pupils): {cov['n_cz']} CZs, {cov['n_districts']:,} districts, {cov['pupils'] / 1e6:.2f} million pupils "
      f"(effective number of CZs {cov['n_eff_cz']:.0f}, the five largest carry {pc(cov['top5_weight_share'], 0)} of the weight). "
      f"Headline matched sample (residents counted once): {sc['n_matched']:,} districts ({nm1} city, {nm2} suburban), {pc(sc['share_pupils_matched'])} of the pupils; "
      f"the {cov['n_districts'] - sc['n_matched']} others have no ACS record ({cov['unmatched_reasons']['no_acs_record']} under the earnings rule) "
      f"or no audited FTE or pupil count ({cov['unmatched_reasons']['no_fte_or_pupils']}). In the headline sample: "
      f"{sc['edu_res'] / 1e6:.2f} million EDU residents, {sc['employed_res'] / 1e6:.1f} million employed residents, {sc['fte'] / 1e6:.2f} million K-12 FTE teachers "
      f"({pc(sc['fte'] / cov['ccd_fte_k12_all_regular_districts'], 0)} of the {cov['ccd_fte_k12_all_regular_districts'] / 1e6:.2f} million "
      f"audited K-12 FTE of all regular districts; {sc['n_secondary_without_residents']} secondary LEAs contribute FTE and pupils only).\n")
    a("### EDU moments (suburb minus city, log gaps)\n")
    a("| moment | what | value | boot SE | ACS sampling SE | 95% interval | n CZs |")
    a("|---|---|---|---|---|---|---|")
    for k, s in S.items():
        acs = f"{s['se_acs_sampling']:.4f}" if "se_acs_sampling" in s else ""
        a(f"| `{k}` | {s['label']} | {s['value']:.4f} | {s['se']:.4f} | {acs} | {fci(s, 4)} | {s['n_cz']} |")
    a(f"\n`net_commute_gap_edu` equals `gap_edu_res` - `gap_fte` CZ by CZ (largest absolute difference {nc['identity_max_abs_diff']:.1e}), "
      f"so it is also the log gap of EDU residents per K-12 FTE. `gap_fte` = `gap_pupils_k12` + `gap_fte_per_pupil` "
      f"(the matched-sample pupils gap is {S['gap_pupils_k12']['value']:.4f} against the tracked `gap_pupils` {D['gap_pupils']['value']:.4f}), and "
      f"the matched `gap_fte_per_pupil` {S['gap_fte_per_pupil']['value']:.4f} is the adopted `gap_teachers_pp` "
      f"({D['gap_teachers_pp']['value']:.4f}) on the districts that also have an ACS record. The ACS sampling SE is computed from the published "
      f"margins of error of C24010 (90% MOE / 1.645, men's and women's EDU counts added, location counts by the root-sum-of-squares rule, independence across "
      f"districts and CZs, delta method for the log gaps; FTE are exact counts): small next to the bootstrap SE, which already contains most of it.\n")
    a("### Levels (CZ-weighted means of the location values, on the CZs that have the gap)\n")
    a("| level | city | suburb |")
    a("|---|---|---|")
    for k, _, lab in LEVELS_TR:
        d_ = 2 if k == "edu_per_fte" else 4
        a(f"| {lab} | {fs(lv[k]['loc1'], d_)} | {fs(lv[k]['loc2'], d_)} |")
    a("\nEDU residents per FTE is about 3-4: the group includes postsecondary teachers, preschool teachers, teaching assistants, "
      "librarians and private-school staff that are not in the public K-12 FTE; only the gaps are comparable with the model.\n")
    a("### Bootstrap correlation with the four adopted targets\n")
    keys = D["corr_keys"]
    a("| | " + " | ".join(f"`{k}`" for k in keys) + " |")
    a("|---|" + "---|" * len(keys))
    for i, k in enumerate(keys):
        a(f"| `{k}` | " + " | ".join(f"{v:.2f}".replace("-0.00", "0.00") for v in D["corr_matrix"][i]) + " |")
    a("\n`k12_*` are the K-12 counterparts of section 3 (district gap and composition term resampled together). `gap_edu_res` and `gap_fte` correlate 0.9 or more "
      "with `gap_pupils`: both are size gaps. The net gaps are nearly uncorrelated with the size gap and correlate with `gap_teachers_pp` and "
      "`gap_salary`, the information that diagonal weights discard. The full draws are in `teacher_residence_draws.csv`.\n")
    a("### Overlap of the ACS layers and other variants\n")
    a(f"Layers in the earnings-rule matched sample (share of employed residents): unified {pc(cov['employed_res_share_by_layer']['unified'])}, "
      f"elementary {pc(cov['employed_res_share_by_layer']['elementary'])}, secondary {pc(cov['employed_res_share_by_layer']['secondary'])}; "
      f"pupils: unified {pc(cov['pupil_share_by_layer']['unified'])}, elementary {pc(cov['pupil_share_by_layer']['elementary'])}, "
      f"secondary {pc(cov['pupil_share_by_layer']['secondary'])}. {cov['n_cz_with_secondary_lea']} CZs, carrying "
      f"{pc(cov['pupil_weight_cz_with_secondary_lea'], 0)} of the weight, have a secondary LEA. Summing every layer's record counts "
      f"{cov['overlap_factor_employed_residents']:.3f} times the residents of the unified and elementary records. In the CZ locations that hold a secondary LEA "
      f"the elementary and secondary layers hold nearly the same number of employed residents (city {cov['secondary_vs_elementary_employed_residents']['by_location'][1]['elementary'] / 1e6:.2f} against "
      f"{cov['secondary_vs_elementary_employed_residents']['by_location'][1]['secondary'] / 1e6:.2f} million, suburb "
      f"{cov['secondary_vs_elementary_employed_residents']['by_location'][2]['elementary'] / 1e6:.2f} against "
      f"{cov['secondary_vs_elementary_employed_residents']['by_location'][2]['secondary'] / 1e6:.2f}; the elementary layer is "
      f"{100 * (1 - cov['secondary_vs_elementary_employed_residents']['by_location'][1]['elementary'] / cov['secondary_vs_elementary_employed_residents']['by_location'][1]['secondary']):.1f}% and "
      f"{100 * (1 - cov['secondary_vs_elementary_employed_residents']['by_location'][2]['elementary'] / cov['secondary_vs_elementary_employed_residents']['by_location'][2]['secondary']):.1f}% smaller), "
      "so counting either layer once leaves at most a small coverage gap (elementary districts below the 500-pupil screen or outside the partition).\n")
    tl = cov.get("tiling") or {}
    if tl:
        a(f"The unified and elementary records tile the states: their employed residents are {100 * tl['national_unified_plus_elementary']:.2f}% of the "
          f"state totals of the same table, and the secondary layer adds {100 * (tl['national_all_layers'] - tl['national_unified_plus_elementary']):.1f}% "
          "nationally (" + ", ".join(f"{t['state']} {100 * (t['all_layers'] - 1):.0f}%" for t in tl["top_states"][:4]) + "). The state sums "
          f"of EDU residents ({tl['edu_state_sum'] / 1e6:.3f} million) match the PUMS count ({C['counts']['edu'] / 1e6:.3f} million).\n")
    a("| variant | `gap_edu_res_emp` | `gap_edu_res` | `gap_fte` | `net_commute_gap_edu` | n CZs |")
    a("|---|---|---|---|---|---|")
    a(f"| **headline: residents counted once (unified + elementary records); FTE and pupils from every LEA** | {fs(gt)} | {fs(ge)} | {fs(gf)} | {fs(nc)} | {nc['n_cz']} |")
    names = {"no_overlap_sec": "residents once from unified + secondary records (elementary residents dropped where a secondary LEA is present)",
             "earnings_rule": "earnings rule: every LEA takes its own layer's record (overlapping residents counted twice)",
             "unified": "unified districts only (elementary and secondary LEAs leave every sum)",
             "nyc_dropped": "New York City dropped from both the resident and the FTE / pupil counts",
             "drop_ny_cz": "New York CZ left out",
             "unmatched": "each statistic on its own maximal sample"}
    for k, lab in names.items():
        if k not in V:
            continue
        v = V[k]
        cells = [fs(v[x]) if x in v else "" for x in ("gap_edu_res_emp", "gap_edu_res", "gap_fte", "net_commute_gap_edu")]
        a(f"| {lab} | " + " | ".join(cells) + f" | {v['net_commute_gap_edu']['n_cz']} |")
    a(f"\nThe earnings rule suits ratios and medians, not counts: it adds {100 * (cov['overlap_factor_employed_residents'] - 1):.1f}% to the employed residents, more in the suburbs "
      f"(secondary records add {pc(cov['overlap_factor_by_location'][1] - 1)} to the city's and {pc(cov['overlap_factor_by_location'][2] - 1)} to the suburb's employed residents in the matched sample), so it "
      f"inflates the suburban count: it moves `net_commute_gap_edu` by {er['net_commute_gap_edu']['value'] - nc['value']:+.3f} ({(er['net_commute_gap_edu']['value'] - nc['value']) / nc['se']:.1f} SE) and "
      f"`gap_edu_res` by {er['gap_edu_res']['value'] - ge['value']:+.3f}; the two single-count constructions agree within "
      f"{abs(V['no_overlap_sec']['net_commute_gap_edu']['value'] - nc['value']):.3f}. The ratio `gap_edu_res_emp` moves by {er['gap_edu_res_emp']['value'] - gt['value']:+.4f} (numerator "
      f"and denominator are overcounted alike). The largest change in `net_commute_gap_edu` across the variants is {vmax[0]:.3f} ({vmax[0] / nc['se']:.1f} SE; `{vmax[1]}`). "
      f"The unified-only sample drops {nc['n_cz'] - V['unified']['net_commute_gap_edu']['n_cz']} CZs; leaving New York out or dropping it from both counts changes the net gap by "
      f"{abs(V['drop_ny_cz']['net_commute_gap_edu']['value'] - nc['value']):.3f} and {abs(V['nyc_dropped']['net_commute_gap_edu']['value'] - nc['value']):.3f}.\n")
    a("### Reproduction check\n")
    a("| target | adopted (value, SE) | this script (value, SE) | n CZs |")
    a("|---|---|---|---|")
    for r in D["reproduction"]:
        a(f"| `{r['moment']}` | {r['reported_value']:.4f} ({r['reported_se']:.4f}) | {r['value']:.4f} ({r['se']:.4f}) | {r['n_cz']} |")
    a(f"\nAll four match to the printed digits: **{'yes' if all(r['ok_value'] and r['ok_se'] in (True, None) for r in D['reproduction']) else 'NO'}** "
      f"(standard errors checked: {'yes' if D['full_run'] else 'no, B differs from the reported runs'}).\n")

    a("## 3. Composition of SOC 25-0000 (d) and the K-12 counterparts\n")
    cn = C["counts"]
    a(f"ACS 2018 1-year PUMS, employed 16+ in households: {cn['edu'] / 1e6:.2f} million workers ({C['n_sample']:,} records) in education, "
      f"training and library occupations, {cn['k12'] / 1e6:.2f} million of them K-12 teachers (OCCP 2310, 2320, 2330), "
      f"{cn['pub_k12'] / 1e6:.2f} million in the public sector.\n")
    a("| occupation (OCCP; SOCP in the records) | workers (millions) | share of EDU (SE) |")
    a("|---|---|---|")
    for c, v in C["share_by_occ"].items():
        a(f"| {c} {v['label'].split(' (')[0]} (SOCP {', '.join(C['socp_by_occ'][c])}) | {v['value'] * cn['edu'] / 1e6:.3f} | {pcs(v, 2)} |")
    a(f"| **K-12 teachers (2310 + 2320 + 2330)** | {cn['k12'] / 1e6:.3f} | **{pcs(C['share_k12'], 2)}** |")
    a(f"| K-12 plus preschool and kindergarten (+ 2300) | {cn['k12p'] / 1e6:.3f} | {pcs(C['share_k12p'], 2)} |")
    a(f"| public-sector K-12 (COW 3-4) | {cn['pub_k12'] / 1e6:.3f} | **{pcs(C['share_pub_k12'], 2)}** |")
    a(f"| public K-12 in elementary and secondary schools (INDP 7860) | | {pcs(C['share_pub_k12_sch'], 2)} |")
    a("\nAs a check on the count, the PUMS public K-12 teachers (headcount) are "
      f"{cn['pub_k12'] / 1e6:.2f} million against {cov['ccd_fte_k12_all_regular_districts'] / 1e6:.2f} million audited FTE of regular districts.\n")
    a("### By class of residence PUMA\n")
    a("Share of EDU workers (by residence PUMA) in each occupation group; SE in parentheses.\n")
    a("| | city | suburb | other | metropolitan | suburb - city | log ratio suburb / city |")
    a("|---|---|---|---|---|---|---|")
    a(f"| EDU workers' distribution over classes | {pcs(sh['city']['edu_workers_share'])} | {pcs(sh['suburb']['edu_workers_share'])} | "
      f"{pcs(sh['other']['edu_workers_share'])} | | | |")
    for k, lab in (("k12", "K-12 teachers"), ("k12p", "K-12 + preschool / kindergarten"), ("pub_k12", "public K-12 teachers"),
                   ("pub_k12p", "public K-12 + preschool / kindergarten")):
        sp = C["split"][k]
        a(f"| {lab} | {pcs(sh['city'][k])} | {pcs(sh['suburb'][k])} | {pcs(sh['other'][k])} | {pcs(sh['metro'][k])} | "
          f"{100 * sp['diff']['value']:.1f} pp ({100 * sp['diff']['se']:.1f}) | {fs(sp['log_ratio'])} |")
    for o, lab in EDU_CODES.items():
        sp = C["split"]["occ"][str(o)]
        a(f"| {o} {lab.split(' (')[0]} | {pcs(sh['city']['occ'][str(o)])} | {pcs(sh['suburb']['occ'][str(o)])} | "
          f"{pcs(sh['other']['occ'][str(o)])} | | {100 * sp['diff']['value']:.1f} pp ({100 * sp['diff']['se']:.1f}) | {fs(sp['log_ratio'])} |")
    a("\nThe city's EDU workforce is tilted toward postsecondary teachers (colleges and universities are in principal cities) and away "
      "from elementary and secondary teachers, so residents' EDU gaps understate the suburb's teacher residents.\n")
    a("**Is the PUMA classification defensible?** It uses only public Census files (the 2010 population of every place within every PUMA and the OMB "
      "list of principal cities), uses the principal-city criterion that underlies the NCES city locale and passes the spot checks of section 1. It is not the "
      "NCES locale partition of the districts: a PUMA holds 100,000 people or more, so PUMAs near the cut mix city and suburb, and PUMAs in "
      "outlying counties of a metropolitan area count as suburb whether or not they are in the urbanized area. The composition result does not "
      f"depend on it in sign: the log ratio of the public K-12 share is {min(r_['delta_log_s_pub']['value'] for r_ in R['class_sensitivity']):.2f}-"
      f"{max(r_['delta_log_s_pub']['value'] for r_ in R['class_sensitivity']):.2f} on national classes and "
      f"{min(r_['within_cz']['pub_k12']['value'] for r_ in R['class_sensitivity']):.2f}-"
      f"{max(r_['within_cz']['pub_k12']['value'] for r_ in R['class_sensitivity']):.2f} inside CZs across the five cuts of the table below.\n")
    a("### Dispersion across PUMAs\n")
    dp = C["dispersion"]
    a(f"Share of a PUMA's EDU workers who are K-12 teachers, {dp['k12']['n_puma']:,} PUMAs with on average "
      f"{dp['k12']['mean_edu_sample_per_puma']:.0f} sampled EDU workers (minimum {dp['k12']['min_edu_sample_per_puma']}), so a PUMA's share "
      f"has a sampling SD of about {dp['k12']['mean_sampling_sd_w']:.2f}. 'Net of sampling' subtracts the average replicate-weight "
      "sampling variance of the PUMA shares from the variance across PUMAs (weights: the PUMA's EDU workers).\n")
    a("| | K-12 teachers | public K-12 teachers |")
    a("|---|---|---|")
    rows = [("mean (EDU-weighted)", "mean_w"), ("SD across PUMAs (EDU-weighted)", "sd_w"),
            ("SD net of sampling noise (EDU-weighted)", "sd_w_net_of_sampling"), ("SD across PUMAs (unweighted)", "sd_u"),
            ("SD net of sampling noise (unweighted)", "sd_u_net_of_sampling"), ("10th percentile", "p10"), ("25th percentile", "p25"),
            ("median", "p50"), ("75th percentile", "p75"), ("90th percentile", "p90")]
    for lab, k in rows:
        a(f"| {lab} | {fs(dp['k12'][k])} | {fs(dp['pub_k12'][k])} |")
    a(f"\nThe composition of EDU workers differs across areas by about {dp['k12']['sd_w_net_of_sampling']['value']:.1f} (one SD, net of sampling "
      f"error) around a mean share of {dp['k12']['mean_w']['value']:.2f}, of the same order as the {100 * C['split']['k12']['diff']['value']:.0f} point city / suburb "
      "difference in the K-12 share, so one area's EDU count is a noisy proxy for its teachers; averages over the CZs of the partition are "
      "less noisy.\n")
    a("### The K-12 counterparts of the model moments (approximate)\n")
    a(f"Inside each CZ of the partition the PUMAs are assigned to the CZ that holds most of their 2010 population (Dorn's crosswalk, as in the district panel); "
      f"the share s of EDU workers who are public K-12 teachers is computed for the city and the suburban PUMAs of the CZ (class cut as above), "
      f"log(suburb / city) per CZ, and averaged with the CZ pupil weights over the {cp['n_cz']} CZs with PUMAs of both classes and a positive share in both "
      f"({pc(cp['weight_share'], 1)} of the weight; {cp['n_sample']:,} sampled EDU workers). Since EDU residents = public K-12 residents / s, "
      "Delta log (public K-12 residents) = Delta log (EDU residents) + Delta log s: the K-12 moments are the headline EDU moments of section 2 plus this term "
      "(the district gaps of all CZs plus the mean of the term over the CZs that have it).\n")
    a(f"Standard errors. The district gap and the composition term are resampled TOGETHER over whole CZs with the draws of `district_settings.py` "
      f"(the two are correlated: {h_n['corr_district_vs_composition']:.2f} for the net gap; the CZ-bootstrap SE of the composition term is {cp['se_cz_bootstrap']:.3f} against its "
      f"replicate-weight survey SE {cp['se_survey']:.3f}); the survey variance of the composition term and the ACS sampling variance of C24010 are added. "
      "The CZ bootstrap already contains part of that sampling noise (each CZ's PUMS and ACS estimates are noisy and the bootstrap resamples them), so the "
      "total SE is conservative; `se_boot` alone is the lower figure. Adding the district SE and the survey SE in quadrature instead would ignore the "
      "CZ variance of the composition term and its covariance with the district gap.\n")
    a("| moment | what | value | SE | bootstrap | survey | ACS | 95% interval | range over specifications | n CZs |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for k, h in H.items():
        a(f"| `{k}` | {h['label']} | **{h['value']:.4f}** | {h['se']:.4f} | {h['se_boot']:.4f} | {h['se_survey']:.4f} | {h['se_acs']:.4f} | {ci(h)} | "
          f"{rg(h)} ({h['spec_n']} specifications) | {h['n_cz']} |")
    a("\nValue of each K-12 moment under the specifications (district gap of the variant plus the composition term of the specification):\n")
    a("| specification | `net_commute_gap` | `gap_teach_res_emp` | `gap_k12_res` |")
    a("|---|---|---|---|")
    for x in KB["specs"]:
        cells = ["" if x[k] is None else f"{x[k]:.3f}" for k in ("net_commute_gap", "gap_teach_res_emp", "gap_k12_res")]
        a(f"| {'**' + x['label'] + '**' if x['key'] == 'headline' else x['label']} | " + " | ".join(cells) + " |")
    gq = KB["composition"]["group_quarters"]
    a(f"\nCross-checks that do not use the district table: **`gap_teach_res_emp` measured directly in the PUMS** (public K-12 teachers / all employed "
      f"residents, suburban against city PUMAs of the same CZ, CZ-weighted, {d1['n_cz']} CZs) is {fs(d1)} (survey and CZ-bootstrap variance added); "
      f"on the national POWPUMA classes it is {fs(dr['gap_teach_res_emp_powpuma_classes_national'])}; the national net commuting gap of public K-12 teachers on POWPUMA "
      f"classes (section 4) is {fs(dr['net_commute_gap_powpuma_classes_national'])}. "
      f"The composition term with group quarters included (C24010 counts them, the household PUMS sample does not): {gq['value']:.3f} against {cp['value']:.3f}; "
      f"{100 * gq['gq_share_of_edu_workers']:.1f}% of EDU workers live in group quarters ({100 * gq['gq_share_by_class']['city']:.1f}% in city PUMAs, "
      f"{100 * gq['gq_share_by_class']['suburb']:.1f}% in suburban ones).\n")
    a("The composition term and the extra commute of working in a POWPUMA of the other class under other cuts of the PUMA classes (pc_share is the share of a PUMA's 2010 "
      "population living in principal-city places; the metropolitan cut is the share in metropolitan counties):\n")
    a("| city if pc_share >= | suburb if metro share >= | PUMAs city / suburb | Delta log s, national classes | Delta log s, inside the CZs (n CZs) | K-12 `net_commute_gap` | extra minutes, other class of POWPUMA | share crossing |")
    a("|---|---|---|---|---|---|---|---|")
    cutrows = {x["key"]: x for x in KB["specs"]}
    for r_ in R["class_sensitivity"]:
        if r_["main"]:
            net_v = h_n["value"]
        else:
            net_v = cutrows[f"cut_{r_['pc_cut']}_{r_['metro_cut']}"]["net_commute_gap"]
        a(f"| {r_['pc_cut']:.2f}{' (main)' if r_['main'] else ''} | {r_['metro_cut']:.2f} | {r_['n_puma'].get('city', 0)} / {r_['n_puma'].get('suburb', 0)} | "
          f"{fs(r_['delta_log_s_pub'])} | {fs(r_['within_cz']['pub_k12'])} ({r_['within_cz']['pub_k12']['n_cz']}) | {net_v:.3f} | "
          f"{r_['extra_minutes_class']['coef']:.2f} ({r_['extra_minutes_class']['se']:.2f}) | {100 * r_['extra_minutes_class']['share_cross']:.1f}% |")
    a("")

    a("## 4. Commuting (e)\n")
    a("ACS 2018 1-year PUMS, employed 25-64 in households. Commuters are workers with a positive one-way travel time (JWMNP > 0) "
      "who do not work at home; 'outside' is a place of work in another POWPUMA (or state) than the POWPUMA that contains the worker's "
      "residence PUMA (Census PUMA to POWPUMA conversion of the 2018 PUMS code lists), among commuters with a US place of work. Teachers: "
      "OCCP 2310, 2320, 2330. Reported times heap at multiples of five: "
      f"{pcs(kt['eq30'])} of teachers and {pcs(ko['eq30'])} of other workers report exactly 30 minutes, so 'more than 30' and '30 or more' differ by that much; the "
      "literature figure (Santelli and Grissom, 19% of Nashville teachers) is a '30 or more' share.\n")
    a("| group | n | mean one-way minutes | share > 30 min | share >= 30 min | share working at home | share working outside the residence POWPUMA |")
    a("|---|---|---|---|---|---|---|")
    glabs = [("T", "K-12 teachers"), ("T_f", "  women"), ("T_m", "  men"), ("Tpub", "  public-sector"),
             ("O", "other workers"), ("O_f", "  women"), ("O_m", "  men"), ("Obach", "  with a bachelor's degree or more"),
             ("Otherpub", "  public-sector (non-teachers)"), ("T@metro", "K-12 teachers, metropolitan PUMAs"),
             ("O@metro", "other workers, metropolitan PUMAs"), ("T@city", "K-12 teachers living in city PUMAs"),
             ("T@suburb", "K-12 teachers living in suburban PUMAs"), ("O@city", "other workers, city PUMAs"),
             ("O@suburb", "other workers, suburban PUMAs")]
    for k, lab in glabs:
        s = K["stats"][k]
        a(f"| {lab} | {K['n_sample'][k]:,} | {fs(s['mean'], 2)} | {pcs(s['gt30'])} | {pcs(s['ge30'])} | {pcs(s['home'])} | {pcs(s['out'])} |")
    a("\nTeachers minus other workers (SE):\n")
    a("| sample | mean minutes | share > 30 min | share >= 30 min | share at home | share outside POWPUMA |")
    a("|---|---|---|---|---|---|")
    dl = {"all": "all", "_f": "women", "_m": "men", "@metro": "metropolitan PUMAs", "@city": "living in city PUMAs",
          "@suburb": "living in suburban PUMAs"}
    for k, v in K["diff_T_minus_O"].items():
        a(f"| {dl[k]} | {fs(v['mean'], 2)} | {100 * v['gt30']['value']:.1f} pp ({100 * v['gt30']['se']:.1f}) | {100 * v['ge30']['value']:.1f} pp ({100 * v['ge30']['se']:.1f}) | "
          f"{100 * v['home']['value']:.1f} pp ({100 * v['home']['se']:.1f}) | {100 * v['out']['value']:.1f} pp ({100 * v['out']['se']:.1f}) |")
    v = K["diff_T_minus_Obach"]
    a(f"| vs workers with a bachelor's degree or more | {fs(v['mean'], 2)} | {100 * v['gt30']['value']:.1f} pp ({100 * v['gt30']['se']:.1f}) | {100 * v['ge30']['value']:.1f} pp ({100 * v['ge30']['se']:.1f}) | "
      f"{100 * v['home']['value']:.1f} pp ({100 * v['home']['se']:.1f}) | {100 * v['out']['value']:.1f} pp ({100 * v['out']['se']:.1f}) |")
    a("\n### Extra commute from working elsewhere\n")
    a("Weighted regression of one-way minutes on the crossing indicator, sex, age and age squared with residence-PUMA fixed effects "
      "(workers who live in the same PUMA are compared); replicate-weight SEs. Two definitions of crossing: another POWPUMA (or state) "
      "than the one that contains the residence PUMA, and a POWPUMA of the other class (city / suburb) than the one that contains the residence PUMA "
      "(both ends on the POWPUMA geography; both ends city or suburb).\n")
    a("| sample | crossing | extra minutes (SE) | n | n crossing | share crossing |")
    a("|---|---|---|---|---|---|")
    flab = {"T": "K-12 teachers", "T@metro": "K-12 teachers, metropolitan PUMAs", "Tpub@metro": "public K-12 teachers, metropolitan",
            "T_f@metro": "women teachers, metropolitan", "T_m@metro": "men teachers, metropolitan",
            "T@suburb": "teachers living in suburban PUMAs", "T@city": "teachers living in city PUMAs",
            "O": "other workers", "O@metro": "other workers, metropolitan", "Obach@metro": "other workers with a bachelor's degree, metropolitan",
            "T@metro:type": "K-12 teachers living and working in a city or suburban POWPUMA",
            "O@metro:type": "other workers living and working in a city or suburban POWPUMA",
            "T@metro:type_same_class": "K-12 teachers, as above, residence PUMAs whose own class equals their POWPUMA's",
            "T@metro:type_mixed": "K-12 teachers, SUPERSEDED: class of the residence PUMA against the class of the work POWPUMA"}
    for k, lab in flab.items():
        e = fe[k]
        crossing = {"pow": "other POWPUMA or state", "type": "POWPUMA of the other class (city / suburb)",
                    "type_mixed": "mixed geographies (see below)"}[e["kind"]]
        a(f"| {lab} | {crossing} | {e['coef']:.2f} ({e['se']:.2f}) | {e['n']:,} | {e['n_cross']:,} | {100 * e['share_cross']:.1f}% |")
    rw = K["raw_outside_minus_inside"]
    a(f"\nWithout controls, K-12 teachers in metropolitan PUMAs commute {rw['T@metro']['outside']:.1f} minutes when they work outside the "
      f"residence POWPUMA and {rw['T@metro']['inside']:.1f} when they work inside: {rw['T@metro']['diff']:.1f} minutes (SE {rw['T@metro']['se']:.2f}); "
      f"other workers {rw['O@metro']['diff']:.1f}. **Bracketing.** Classifying the residence by its PUMA and the workplace by its POWPUMA (two "
      f"geographies; the SUPERSEDED row of the table) codes {100 * fe['T@metro:type_mixed']['share_cross']:.0f}% of metropolitan teachers as crossing against "
      f"{100 * fe['T@metro:type']['share_cross']:.0f}% on one geography, most of them working in their own residence POWPUMA, and biases the "
      f"estimate down ({fe['T@metro:type_mixed']['coef']:.1f} against {fe['T@metro:type']['coef']:.1f} minutes); it is not used. "
      "Both definitions used now are crossings of a county line (POWPUMAs are whole counties or groups of counties): the PUMS cannot see a move between "
      "a central city and its suburbs inside one county, and a county-line crossing is on average longer than a district-to-district move inside a county, "
      "so both can OVERSTATE the increment of working in the other location. Santelli and Grissom's out-of-county figure for Nashville "
      f"({NASHVILLE_EXTRA_MIN} minutes) anchors the low end of the range.\n")
    a("### Residence against workplace of K-12 teachers on one geography (national)\n")
    a("Both residence and workplace are classified by the POWPUMA (the finest area on which the PUMS reports the workplace), all employed "
      "persons 16+ with a US POWPUMA; shares of all workers, so the log gaps are those of counts.\n")
    a("| | residence: city | suburb | other | workplace: city | suburb | other | gap by residence | gap by workplace | net commuting gap | gap, residents per employed resident | works in the other class | city <-> suburb |")
    a("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for k, lab in (("k12", "K-12 teachers"), ("pub_k12", "public K-12 teachers")):
        f_ = F[k]
        a(f"| {lab} | " + " | ".join(f"{100 * f_['res'][c]['value']:.2f}%" for c in CLASSES) + " | "
          + " | ".join(f"{100 * f_['work'][c]['value']:.2f}%" for c in CLASSES) + f" | {fs(f_['gap_res'])} | {fs(f_['gap_work'])} | "
          f"**{fs(f_['net_commute_gap'])}** | {fs(f_['gap_res_emp'])} | {pcs(f_['share_cross_class'])} | {pcs(f_['share_metro_cross'])} |")
    cmp_ = F["edu_res_vs_pubk12_work"]
    a(f"\nEDU workers by residence against public K-12 teachers by workplace (the PUMS version of the district measure): gap by residence "
      f"{fs(cmp_['gap_edu_res'])}, by workplace {fs(cmp_['gap_pubk12_work'])}, net {fs(cmp_['net_commute_gap'])}; of which composition "
      f"(EDU against public K-12 residents) {fs(cmp_['composition'])}. The sign agrees with the district measure and the size is smaller; "
      "the classes are national and county-based.\n")
    a("Where the teachers who work in a class of area live (public K-12 teachers; share of those working in the class):\n")
    a("| work in | live in city | live in suburb | live in other |")
    a("|---|---|---|---|")
    lg = F["pub_k12"]["live_given_work"]
    for b_ in ("city", "suburb"):
        a(f"| {b_} POWPUMAs | " + " | ".join(pcs(lg[f"work_{b_}"][f"live_{a_}"]) for a_ in CLASSES) + " |")
    ny = R["nyc_check"]
    a(f"\n**Check against the IBO count for New York City.** Public K-12 teachers (COW 3-4) with a workplace in the five NYC POWPUMAs "
      f"({ny['workers_pub_k12'] / 1e3:.0f} thousand weighted, {ny['n_pub_k12']:,} records): **{pcs(ny['pub_k12_live_in_nyc'])} live in the NYC PUMAs** "
      f"({pcs(ny['k12_live_in_nyc'])} for all K-12 teachers, {pcs(ny['all_workers_live_in_nyc'])} for all commuters) against 67% of the "
      f"70,328 teachers in the IBO count for 2011-12; {pcs(ny['pub_k12_live_in_workplace_powpuma'])} live in the borough group (POWPUMA) where "
      "they work. The PUMS counts persons who report a teaching occupation (aides and substitutes included) in 2018, the IBO the "
      "Department of Education's teachers in schools, so the two differ by a few points.\n")
    a("### Literature\n")
    a("| source | what it shows |")
    a("|---|---|")
    for lit in R["literature"]:
        a(f"| [{lit['ref']}]({lit['url']}) | " + "; ".join(lit["facts"]) + f" (read: {lit['read']}) |")
    a("")

    a("## 5. The commuting cost: r_k and tau_k (f)\n")
    a("Annual time cost = 2 x (extra one-way minutes / 60) x days a year x (value of an hour). "
      "r_k = annual time cost / Cbar (M2_SPEC.md, a goods cost in units of mean consumption); tau_k = annual time cost / annual earnings "
      "(M2_SPEC_v2.md: income x (1 - tau_k) when r != w; the value of an hour is a share of the worker's own wage, so "
      "tau_k = share x 2 x extra hours x days / annual hours).\n")
    a("| input | low | central | high | source |")
    a("|---|---|---|---|---|")
    a(f"| extra one-way minutes of working in the other location, teachers | {RK['delta_min']['low']:.1f} | {RK['delta_min']['central']:.1f} | "
      f"{RK['delta_min']['high']:.1f} | central = average of the two main teacher estimates in section 4 (crossing a POWPUMA {fe['T@metro']['coef']:.1f}, crossing to a POWPUMA "
      f"of the other class {fe['T@metro:type']['coef']:.1f}; SE {RK['delta_se']:.2f}); low / high = the smallest and largest teacher estimate over the specifications (POWPUMA crossing, "
      "the POWPUMA-class crossing under five PUMA class cuts, Nashville); the 95% interval of the central is "
      f"{RK['delta_min']['ci_low']:.1f}-{RK['delta_min']['ci_high']:.1f} |")
    a(f"| extra one-way minutes, other workers | {RK['delta_min_other']['low']:.1f} | {RK['delta_min_other']['central']:.1f} | "
      f"{RK['delta_min_other']['high']:.1f} | the same two estimators on other workers 25-64 ({fe['O@metro']['coef']:.1f} and {fe['O@metro:type']['coef']:.1f}) and the class cuts |")
    a(f"| days a year, teachers | {RK['days']['low']:.0f} | {RK['days']['central']:.0f} | {RK['days']['high']:.0f} | NCES: 179 school days "
      "(Digest 2017, Table 203.90, 2011-12), plus about a week of contract days (assumption) |")
    a(f"| days a year, other workers | {RK['days_other']['low']:.0f} | {RK['days_other']['central']:.0f} | {RK['days_other']['high']:.0f} | 50 weeks of 5 days less holidays and days at home (assumption) |")
    a(f"| value of an hour (2018 dollars) | {RK['vtt_2018']['low']:.2f} | {RK['vtt_2018']['central']:.2f} | {RK['vtt_2018']['high']:.2f} | "
      "USDOT 2016: 35% / 50% / 60% of hourly median household income ($27.20 in 2015), x CPI-U 1.0594 |")
    a(f"| annual hours worked (full-time year-round workers 25-64) | | teachers {RK['annual_hours']['teachers']:,.0f}, others {RK['annual_hours']['other_workers']:,.0f} | | "
      f"PUMS: usual weekly hours x weeks worked (51 for 50-52 weeks, 48.5 for 48-49), mean |")
    a(f"| Cbar per consumer unit | | {RK['cbar_cu']:,.2f} | | CE 2018 (BLS Table R-1; as r_p) |")
    a(f"| Cbar per adult | {RK['cbar_cu'] / RK['adults_per_cu_range'][0]:,.0f} | {RK['cbar_adult']:,.0f} | {RK['cbar_cu'] / RK['adults_per_cu_range'][1]:,.0f} | "
      f"CE 2018 Table 1400: {CE_PERSONS} persons less {CE_CHILDREN} children under 18 = {RK['adults_per_cu']:.1f} adults per consumer unit (one decimal: 1.8-2.0) |")
    a(f"| **r_k per adult** | **{RK['low_per_adult']:.4f}** | **{RK['central_per_adult']:.4f}** | **{RK['high_per_adult']:.4f}** | all-low and all-high combinations; central SE from the extra minutes alone {RK['central_se_from_delta_per_adult']:.4f} |")
    a(f"| r_k per consumer unit | {RK['low_per_cu']:.4f} | {RK['central_per_cu']:.4f} | {RK['high_per_cu']:.4f} | the unit of r_p; central SE from the extra minutes alone {RK['central_se_from_delta_per_cu']:.4f} |")
    a(f"| **tau_k, teachers** | **{tt['low']:.4f}** | **{tt['central']:.4f}** | **{tt['high']:.4f}** | value of time 35% / 50% / 60% of the own wage; central SE from the extra minutes alone {tt['se_from_extra_minutes']:.4f}; with 2,080 hours {tt['central_with_2080_hours']:.4f} |")
    a(f"| **tau_k, other workers** | **{to['low']:.4f}** | **{to['central']:.4f}** | **{to['high']:.4f}** | same shares; central SE {to['se_from_extra_minutes']:.4f}; with 2,080 hours {to['central_with_2080_hours']:.4f} |\n")
    a(f"Central: {RK['annual_hours_central']:.0f} hours a year ({2 * RK['delta_min']['central']:.0f} extra minutes a day) x "
      f"${RK['vtt_2018']['central']:.2f} = ${RK['annual_time_cost_central_usd']:,.0f} a year: {100 * RK['central_per_cu']:.1f}% of Cbar per consumer unit, "
      f"{100 * RK['central_per_adult']:.1f}% of Cbar per adult. For scale, "
      f"full-time year-round K-12 teachers 25-64 earn a median ${ow['median_hourly']:.2f} an hour on reported hours (a heap: "
      f"{pc(ow['share_at_median'], 0)} of them are exactly at it, so its SE is not informative; geometric mean "
      f"${ow['geometric_mean_hourly']:.2f}, SE {ow['se_geometric_mean']:.2f}; other full-time year-round workers' median "
      f"${W['other']['median_hourly']:.2f}), so half of the teacher median is ${0.5 * ow['median_hourly']:.2f}. The model's agents are adults, so the per-adult r_k is "
      f"the one to use with the per-adult consumption scale (DESIGN_v2 'Units'); r_p (housing, `housing_costs.md`, 0.028) is a household cost over household spending and unit-free, "
      "a per-worker cost is not. M2_SPEC_v2 applies one tau_k to every worker: the other-worker value fits the 90%+ of the population that are not teachers; "
      "the teacher value fits teachers (fewer commuting days a year).\n")
    a("One-at-a-time sensitivity (other inputs at their central values; tau_k for teachers):\n")
    a("| change | extra minutes | days | $ an hour | r_k per adult | r_k per consumer unit | tau_k |")
    a("|---|---|---|---|---|---|---|")
    for k, s_ in RK["sensitivity"].items():
        a(f"| {k} | {s_['delta_min']:.1f} | {s_['days']:.0f} | {s_['vtt']:.2f} | {s_['r_k_per_adult']:.4f} | {s_['r_k_per_cu']:.4f} | {s_['tau_k']:.4f} |")
    a("")

    a("## 6. Caveats\n")
    cav = [
        "**The K-12 moments are approximate.** The ACS table has no K-12 detail at the school-district level; the group (EDU) also holds "
        "preschool and postsecondary teachers, aides, librarians and private-school staff, and the K-12 share differs by area (section 3). "
        "The K-12 versions add a PUMS composition term (a different survey, PUMA classes instead of the NCES locale partition, 2018 only, CZs assigned by population share) "
        f"and are reported with a specification range ({rg(h_t)} for `gap_teach_res_emp`, {rg(h_n)} for `net_commute_gap`): the sign of the net gap is not secure, "
        "the level of `gap_teach_res_emp` is. The EDU proxies must not be targeted against a K-12 teacher object. Charter-school teachers are in the residents but "
        "not in the FTE (charter LEAs are not regular districts). The PUMS K-12 count includes aides who report a teaching occupation, "
        "excludes kindergarten teachers (OCCP 2300 mixes them with preschool teachers; the FTE count kindergarten and drop pre-K; the `incl_prek` specification brackets it), "
        "and the household PUMS sample excludes group quarters while C24010 includes them (the `group_quarters` specification adds them).",
        "**Timing.** C24010 is ACS 2014-18 (occupation coded on the 2010 list for 2014-17, the 2018 list for 2018); FTE and pupils are fall "
        "2018; PUMS is 2018 only. School-district boundaries are those of the ACS 2018 vintage.",
        "**Who is in the locations.** Only districts of the partition (>= 500 pupils, city or suburb locale) enter either count, so teachers who "
        "live in town or rural districts, in small districts or outside the CZ and work in the partition count in the FTE and not in the "
        "residents (and the reverse). The net gap is a net of all such flows, not just city-suburb commuting.",
        "**Overlap and NYC.** Residents are counted once (unified + elementary records); the earnings rule that counts overlapping residents twice is a variant "
        "(it biases the net gap upward by the amount in section 2). New York City enters with the Chancellor's record, matching the FTE, pupil and salary "
        "targets; the `nyc_dropped` variant drops it from both the resident and the FTE / pupil counts (the earnings moments drop NYC from the ACS-based statistics only).",
        "**Residents are persons, FTE are full-time equivalents.** Levels of EDU residents per FTE are not comparable with the model's "
        "teacher counts; only the gaps are.",
        "**POWPUMA geography.** POWPUMAs are whole counties or groups of counties, so the PUMS cannot see commuting between a central city and "
        "its suburbs inside one county, and both extra-commute estimates are crossings of a county line, which may overstate the increment of working in the other "
        "location (the Nashville out-of-county figure anchors the low end). PUMAs near the class cut are "
        "classified by a threshold (section 3 shows the cuts); the classes are not the NCES locale partition.",
        "**r_k and tau_k are time costs only.** Vehicle operating costs of the extra miles are not included. For scale (assumptions: 35 mph, 10-30 cents a mile of variable "
        f"cost; the IRS 2018 standard business rate is 54.5 cents a mile, 25 of them depreciation, Notice 2018-03): {RK['vehicle_cost_note']['miles_per_year']:,.0f} extra miles a year would add "
        f"${RK['vehicle_cost_note']['usd'][0]:,.0f}-${RK['vehicle_cost_note']['usd'][1]:,.0f}, i.e. {RK['vehicle_cost_note']['per_cu'][0]:.3f}-{RK['vehicle_cost_note']['per_cu'][1]:.3f} "
        f"per consumer unit and {RK['vehicle_cost_note']['per_adult'][0]:.3f}-{RK['vehicle_cost_note']['per_adult'][1]:.3f} per adult to r_k. "
        "The value of an hour is the USDOT guidance for personal travel (50% of household income), not an estimate for teachers; days a year are assumptions; "
        "the extra minutes are estimated on teachers and on other workers separately.",
        "**Units.** r_p is a cost per household over spending per household and is unit-free; a commuting cost is per worker, so it is divided by CE spending per ADULT "
        f"({RK['adults_per_cu']:.1f} adults per consumer unit from the rounded published counts: 2.5 persons less 0.6 children; the per-adult scale is uncertain by about 5%). "
        "The per-consumer-unit r_k (the unit the task fixed) understates the per-adult value by that factor (about 1.9). tau_k is unit-free.",
        "**Standard errors.** District: whole-CZ bootstrap with the audited panel, screens and partition held fixed; the ACS sampling error of C24010 (margins of error) is "
        "reported separately (small). PUMS: Census replicate weights; the PUMA classes, OMB list and POWPUMA conversion are treated as known. The K-12 moments resample "
        "the district gap and the composition term together over CZs and add the survey and ACS variances: conservative.",
        "**Literature.** Boyd et al. (2005) and Reininger (2012) are about the distance between a teacher's hometown and her school, not "
        "residence and workplace; the IBO figure was not verified against the PDF; no national NCES table of teachers' residence relative to "
        "the district of employment was found in SASS / NTPS.",
    ]
    for t in cav:
        a(f"- {t}")
    a("\n## Files\n")
    a("- `teacher_residence.json`: everything above, machine readable.")
    a("- `teacher_residence_targets.toml`: the K-12 counterparts of the model moments with their SEs and specification ranges (`[k12.*]`), tau_k and r_k "
      "(`[tau_k.*]`, `[r_k.*]`), the EDU district gaps (`[locale.*]`: proxies), their bootstrap correlation with the four adopted targets, the overlap variants, "
      "the composition terms and the commuting statistics.")
    a("- `teacher_residence_draws.csv`: the bootstrap draws (one row per draw) of the four adopted targets, the EDU moments and the K-12 moments (`k12_*`).")
    a("- Raw inputs under `raw/` (gitignored): `acs_c24010_school_districts_2018.parquet`, `acs_pums/acs1_2018_workers/<ST>.parquet`, "
      "`puma_geo/`.")
    path.write_text("\n".join(L) + "\n")


def cmd_build(args: argparse.Namespace) -> None:
    B, seed = args.B, args.seed
    outdir = Path(args.outdir).resolve() if args.outdir else HERE / "estimates" / "settings"
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    dist = district_block(B, seed)
    LOG.info("districts done in %.0fs", time.time() - t0)
    geo = build_puma_geo()
    d, rep = load_pums(geo)
    comp = composition_block(d, rep)
    comm = commute_block(d, rep)
    flows = flows_block(d, rep)
    wages = wage_block(d, rep)
    nyc = nyc_check(d, rep, geo)
    cz_w = dist["g"].set_index("czone").pupils
    sens_cls = class_sensitivity(d, rep, geo, cz_w)
    rk = rk_block(comm, wages, sens_cls)
    within = cz_composition(d, rep, geo, cz_w)
    terms = composition_terms(d, rep, geo, cz_w)
    if abs(terms["pub_k12"]["value"] - within["pub_k12"]["value"]) > 1e-10:
        raise SystemExit("composition_terms does not reproduce cz_composition")
    gq = edu_with_group_quarters(geo, cz_w)
    k12 = k12_block(dist, comp, flows, terms, gq, sens_cls)
    LOG.info("PUMS blocks done in %.0fs", time.time() - t0)

    # bootstrap draws and correlations: the four adopted targets, the EDU moments and the K-12 moments
    Dr = dist["draws"].copy()
    for k, v in k12["_draws"].items():
        Dr[f"k12_{k}"] = v
    corr = Dr.corr()
    std = dist["std"]["stats"]
    districts = dict(
        stats=dist["stats"], levels=dist["levels"], variants=dist["variants"], coverage=dist["coverage"],
        reproduction=dist["reproduction"], full_run=dist["full_run"], hybrid_share=dist["hybrid_share"],
        corr_keys=CORR_KEYS_MD, corr_matrix=corr.loc[CORR_KEYS_MD, CORR_KEYS_MD].to_numpy(float).tolist(),
        corr_keys_toml=CORR_KEYS_TOML, corr_matrix_toml=corr.loc[CORR_KEYS_TOML, CORR_KEYS_TOML].to_numpy(float).tolist(),
        corr_all=corr.round(6).to_dict(), gap_pupils=std["gap_pupils"], gap_teachers_pp=std["gap_teachers_pp"],
        adopted_targets={k: std[k] for k in CORR_TARGETS}, acs_sampling_se=dist["acs_se"])
    R = dict(
        meta=dict(generated=dt.datetime.now().isoformat(timespec="seconds"), script="data/spatial/teacher_residence.py",
                  year=YEAR, B=B, seed=seed, geo_unit=dist["geo_unit"],
                  resampling="whole CZ with replacement; the multiplicity matrix of district_settings.py (n = 188, B, seed)",
                  c24010=C24010, pums_vintage="ACS 2018 1-year PUMS", replicate_weights=N_REP),
        districts=districts, composition=comp, commute=comm, flows=flows, wages=wages, rk=rk, k12=k12,
        class_sensitivity=sens_cls, nyc_check=nyc, composition_within_cz=within,
        geo=geo["info"], literature=LITERATURE,
        pums=dict(n_households=int(len(d)), weighted_millions=float(rep.w0.sum() / 1e6),
                  national_check=d.attrs.get("national_check")))
    clean = ds.clean(strip(R))
    (outdir / "teacher_residence.json").write_text(json.dumps(clean, indent=1))
    Rr = strip(R)
    write_md(outdir / "teacher_residence.md", Rr)
    write_toml(outdir / "teacher_residence_targets.toml", Rr)
    Dr.round(8).to_csv(outdir / "teacher_residence_draws.csv", index_label="draw")
    LOG.info("wrote teacher_residence.{md,json}, teacher_residence_targets.toml and teacher_residence_draws.csv to %s "
             "(%.0fs)", outdir, time.time() - t0)


if __name__ == "__main__":
    raise SystemExit(main())
