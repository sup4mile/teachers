#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acs_occupations.py -- the ACS occupational block of the spatial teachers model
(items T4 "mean wages by occupation" and T7 "2018 occupational block" of
`spatial_calibration.md`), rebuilt from Census API PUMS microdata.

The 2009-13 block in `data/LaborMarketData/wages_occ_shares_v2.xlsx` (built
elsewhere from IPUMS extracts following Hsieh, Hurst, Jones and Klenow 2019,
"HHJK") gives the shares, hourly-wage 90/10 ratios and wage-sample weights that
Table 2 targets.  This script

  1. downloads ACS PUMS microdata for ages 25-34 from the Census API, one state
     at a time, and caches it under raw/acs_pums/ (fetch step);
  2. builds a crosswalk from each Census occupation-code vintage (2002-, 2010-
     and 2018-based; the 2012+ PUMS collapses some 2010 codes) to IPUMS
     `occ1990`, and from `occ1990` to the workbook's 21 groups (HHJK's `occ_broad`
     with K-12 teachers carved out by Census code) (crosswalk section);
  3. applies the HHJK sample and wage definitions (`build_sample`);
  4. computes counts, shares, cell 90/10s, wage-sample counts, the pooled 90/10,
     mean wages / log wages and schooling by group x gender (`compute_moments`);
  5. validates the 2009-13 reproduction against the workbook and writes
     estimates/acs_occupations.json, estimates/acs_occupations.md and the
     crosswalk table estimates/acs_occupations_crosswalk.csv (build step).

    python acs_occupations.py                    # everything (downloads are cached)
    python acs_occupations.py fetch              # downloads only (~10 datasets x 51 states)
    python acs_occupations.py fetch --samples acs1_2016_19
    python acs_occupations.py build              # crosswalk + moments + report (~6 minutes)
    python acs_occupations.py --selftest         # unit checks on the helpers (no network)
    python acs_occupations.py --check-bls-table  # re-read the BLS appendix PDF (needs pdftotext)

The Census API key is read from the environment (CENSUS_API_KEY); if it is not
set, `data/spatial/.env` is parsed.  The key is never printed or written.
Needs pandas, pyarrow, numpy, requests, openpyxl, xlrd (crosswalk .xls files).

Samples (see SAMPLES): `acs5_2009_13` is the 2013 5-year PUMS (data years
2009-13), validated against the workbook with the "workbook" income rules;
`acs1_2009_13` pools the five 1-year PUMS files (the same person records, with
1-year weights and survey-year dollars) under the T7 income rules, the
like-for-like comparator for `acs1_2016_19`, the pooled 2016-19 1-year files.

What the workbook's `readme` leaves out and this script had to infer (all
verified against the workbook, see the report): SCHL >= 12; unemployed persons
with under 15 usual hours stay in as home producers; K-12 teachers are Census
codes 2300-2340; the wage-sample floor is ADJINC-adjusted income >= $1,000;
cell 90/10s are unweighted; hourly wages are reference-period income deflated by
the survey-year CPI-U; and eleven occupation-code departures from the BLS/IPUMS
mapping (`BLS_DEPARTURES_2000`, `OCC1990_OVERRIDES_2010`).

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import os
import re
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import requests
except ImportError:  # pragma: no cover
    print("This script needs `requests`.", file=sys.stderr)
    raise

LOG = logging.getLogger("acs_occupations")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent                       # repository root
RAW = HERE / "raw" / "acs_pums"                 # cache of API downloads
XW_DIR = RAW / "crosswalks"                     # cache of crosswalk source files
META_DIR = RAW / "meta"                         # variables.json snapshots
EST = HERE / "estimates"
WORKBOOK = ROOT / "data" / "LaborMarketData" / "wages_occ_shares_v2.xlsx"
MINCER = ROOT / "data" / "Mincer_YrsSchool.xlsx"

API = "https://api.census.gov/data/{year}/acs/{product}/pums"
API_VARS = API + "/variables.json"

# 50 states + DC (the workbook is from IPUMS ACS, which excludes Puerto Rico)
STATES = [1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22,
          23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
          40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56]
AGE_MIN, AGE_MAX = 25, 34

# Variables requested from every dataset (those it lacks are skipped).  The
# occupation variable differs by vintage: OCCP (1-year files; 2002 codes in
# 2009, 2010 codes 2010-17, 2018 codes 2018+); OCCP02/OCCP10/OCCP12 in the 2013
# 5-year file (data years 2009 / 2010-11 / 2012-13); weeks worked is the bracket
# WKW through 2018 and exact WKWN from 2019; TYPE (TYPEHUGQ from 2019) flags
# group quarters.
WANT_VARS = ["SERIALNO", "SPORDER", "AGEP", "SEX", "SCHL", "RAC1P", "HISP", "ESR",
             "MIL", "WKW", "WKWN", "WKHP", "WAGP", "SEMP", "PERNP", "ADJINC",
             "PWGTP", "COW", "TYPE", "TYPEHUGQ", "OCCP", "OCCP02", "OCCP10",
             "OCCP12", "SOCP", "SOCP00", "SOCP10", "SOCP12"]

# `rules` selects the income definitions (see build_sample):
#   "workbook": reproduces wages_occ_shares_v2.xlsx (floor on ADJINC-adjusted income,
#               hourly wage from reference-period income deflated by the survey-year CPI-U);
#   "t7":       the instructed rule -- ADJINC puts incomes in survey-year dollars, then
#               CPI-U deflates to 1999 dollars (floor: $1,000 in 2010 dollars).
SAMPLES = {
    "acs5_2009_13": dict(label="ACS 2009-13 (5-year PUMS)", product="acs5", years=[2013],
                         rules="workbook"),
    "acs1_2009_13": dict(label="ACS 2009-13 (pooled 1-year PUMS)", product="acs1",
                         years=[2009, 2010, 2011, 2012, 2013], rules="t7"),
    "acs1_2016_19": dict(label="ACS 2016-19 (pooled 1-year PUMS)", product="acs1",
                         years=[2016, 2017, 2018, 2019], rules="t7"),
}


# =============================================================================
# Census API access (key handling, retries, caching)
# =============================================================================

def _read_key() -> str:
    """CENSUS_API_KEY from the environment, else from data/spatial/.env."""
    key = os.environ.get("CENSUS_API_KEY", "").strip()
    if not key:
        env = HERE / ".env"
        if env.exists():
            for line in env.read_text().splitlines():
                m = re.match(r"\s*(?:export\s+)?CENSUS_API_KEY\s*=\s*(.*)", line)
                if m:
                    key = m.group(1).strip().strip("'\"")
    if not key:
        sys.exit("CENSUS_API_KEY not set (environment or data/spatial/.env)")
    return key


class Api:
    """Thin requests wrapper: the key goes in `params`, and every message that
    could contain it is scrubbed."""

    def __init__(self) -> None:
        self.key = _read_key()
        self.s = requests.Session()
        self.s.headers["User-Agent"] = "acs_occupations.py (research; spatial teachers model)"

    def scrub(self, text: str) -> str:
        return str(text).replace(self.key, "***")

    def get_json(self, url: str, params: dict, tries: int = 6, timeout: int = 600):
        p = dict(params)
        p["key"] = self.key
        for k in range(tries):
            try:
                r = self.s.get(url, params=p, timeout=timeout)
                if r.status_code == 204 or not r.content.strip():
                    return []
                r.raise_for_status()
                return r.json()
            except Exception as e:  # noqa: BLE001 -- retry on anything transient
                msg = self.scrub(f"{type(e).__name__}: {e}")[:160]
                LOG.warning("retry %d/%d for %s (%s)", k + 1, tries, url.split("/data/")[-1], msg)
                if k == tries - 1:
                    raise RuntimeError(f"giving up on {url.split('/data/')[-1]}: {msg}") from None
                time.sleep(min(120, 5 * 2 ** k))


def dataset_variables(api: Api, product: str, year: int) -> dict:
    """variables.json for one PUMS dataset, cached (it also carries the code
    labels of the occupation variables)."""
    META_DIR.mkdir(parents=True, exist_ok=True)
    f = META_DIR / f"{product}_{year}_variables.json"
    if f.exists():
        return json.loads(f.read_text())["variables"]
    r = api.s.get(API_VARS.format(year=year, product=product), timeout=300)
    r.raise_for_status()
    f.write_text(r.text)
    return r.json()["variables"]


def fetch_dataset(api: Api, product: str, year: int) -> Path:
    """Download ages 25-34 for every state of one PUMS dataset into
    raw/acs_pums/<product>_<year>/<ST>.parquet (skipping states already cached
    with all the requested columns).  Returns the directory."""
    out = RAW / f"{product}_{year}"
    out.mkdir(parents=True, exist_ok=True)
    avail = dataset_variables(api, product, year)
    want = [v for v in WANT_VARS if v in avail]
    LOG.info("%s %d: %d variables (%s ...)", product, year, len(want), ",".join(want[:6]))
    for st in STATES:
        f = out / f"{st:02d}.parquet"
        if f.exists():
            have = set(pd.read_parquet(f).columns)
            if set(want) <= have:
                continue
        j = api.get_json(API.format(year=year, product=product),
                         {"get": ",".join(want), "for": f"state:{st:02d}",
                          "AGEP": f"{AGE_MIN}:{AGE_MAX}"})
        if not j:
            LOG.warning("%s %d state %02d: empty response", product, year, st)
            continue
        df = pd.DataFrame(j[1:], columns=j[0])
        df = df.loc[:, ~df.columns.duplicated()]        # AGEP comes back twice
        df["ST"] = f"{st:02d}"
        df.to_parquet(f)
        LOG.info("%s %d state %02d: %d rows", product, year, st, len(df))
    return out


def cmd_fetch(samples: list[str]) -> None:
    api = Api()
    for name in samples:
        spec = SAMPLES[name]
        for y in spec["years"]:
            fetch_dataset(api, spec["product"], y)


# =============================================================================
# Occupation groups (HHJK `occ_broad`, K-12 teachers carved out)
# =============================================================================

MARKET_GROUPS = [
    "Executives, Administrative, and Managerial",
    "Management Related",
    "Architects, Engineers, Math, and Computer Science",
    "Natural and Social Scientists, Recreation, Religious, Arts, Athletes",
    "Doctors and Lawyers",
    "Nurses, Therapists, and Other Health Service",
    "Teachers, Postsecondary",
    "Teachers, Non-Postsecondary and Librarians",
    "Health and Science Technicians",
    "Sales, All",
    "Administrative Support, Clerks, Record",
    "Fire, Police, and Guards",
    "Food, Cleaning, and Personal Services and Private Household",
    "Farm, Related Agrigulture, Logging, and Extraction",
    "Mechanics and Construction",
    "Precision Manufacturing",
    "Manufacturing Operators",
    "Fabricators, Inspectors, and Material Handlers",
    "Vehicle Operators",
]
HOME = "Home Production (out of labor force, <15 hrs per week)"
K12 = "Kindergarten - Secondary Teachers"
NONTEACH = MARKET_GROUPS + [HOME]              # the 20 non-teaching groups
ALL_GROUPS = NONTEACH + [K12]                  # the 21 model groups
OTHER_TEACH = MARKET_GROUPS[7]                 # group 8, net of K-12 teachers

# HHJK `occ_broad` (create_2010_IPUMS_extract.do, E3): occ1990 ranges -> group 1..19.
# Codes in no range fall through to 0 ("Home") in the Stata code; we count them.
OCC_BROAD_RANGES = {
    1: [(3, 22)], 2: [(23, 37), (200, 200)], 3: [(43, 68), (867, 867)],
    4: [(69, 83), (166, 177), (183, 199)], 5: [(84, 89), (178, 179)],
    6: [(95, 106), (445, 447)], 7: [(113, 154)], 8: [(155, 165)],
    9: [(203, 235)], 10: [(243, 290)], 11: [(303, 391)], 12: [(413, 427)],
    13: [(403, 408), (433, 444), (448, 469)],
    14: [(473, 499), (613, 617), (868, 868)],
    15: [(503, 599), (866, 866), (869, 869)], 16: [(628, 688)],
    17: [(694, 779)], 18: [(783, 799), (689, 693), (875, 890), (874, 874)],
    19: [(803, 815), (823, 834), (843, 865)],
}
# K-12 teachers are identified by Census occupation code, not by an occ1990 range:
# the codes 2300-2340 (pre-K and kindergarten, elementary and middle, secondary,
# special-education teachers, and "other teachers and instructors"; in 2018 codes
# 2300-2330 plus 2350 tutors and 2360 other teachers and instructors, the two
# pieces of the old 2340).  Teacher assistants (2540) and "other education,
# training and library workers" (2550), which map to occ1990 159 like 2340, stay
# in "Teachers, Non-Postsecondary and Librarians" (occ1990 155-165) with
# counselors, librarians and archivists.  This is the split that reproduces the
# workbook's K-12 and remainder counts exactly, for both genders and all code
# vintages (occ1990 155-158 alone misses code 2340, 2,224 men and 3,396 women;
# 155-159 adds teacher assistants and code 2550, 1,353 men).
K12_CODES = {"02": {"2300", "2310", "2320", "2330", "2340"},
             "10": {"2300", "2310", "2320", "2330", "2340"},
             "18": {"2300", "2310", "2320", "2330", "2350", "2360"}}


def occ1990_to_group(j: int) -> str | None:
    """HHJK `occ_broad` group of an occ1990 code (None if it is in no range)."""
    for g, rngs in OCC_BROAD_RANGES.items():
        for lo, hi in rngs:
            if lo <= j <= hi:
                return MARKET_GROUPS[g - 1]
    return None


# =============================================================================
# Crosswalk: Census occupation code (2002 / 2010 / 2018 vintages) -> occ1990
# =============================================================================
# IPUMS builds `occ1990` from the Census Bureau's double-coded crosswalks by
# modal assignment ("a plurality of records").  IPUMS does not publish the
# resulting OCC->OCC1990 table for the 2010+ codes, so the chain is rebuilt:
#
#   2002 code k --(BLS/IPUMS category of Census 2000 code k/10; BLS_2000_TO_OCC1990)--> j
#   2010 code c --(Census 2002->2010 crosswalk: direct matches + the ACS
#                  conversion rates for split codes, weighted by 2002 counts)--> k
#   2018 code a --(Census 2010->2018 change list, weighted by 2010 counts)--> c
#
# Each stage yields a distribution; the modal occ1990 (the IPUMS convention) is
# used.  Where that distribution straddles two groups the assignment is flagged
# (`group_mass` gives the share of the modal group) and reported.

XW_BASE = "https://www2.census.gov/programs-surveys/demo/guidance/industry-occupation/"
XW_URLS = {
    "OCCBLS_paper.pdf": "https://usa.ipums.org/usa/resources/chapter4/OCCBLS_paper.pdf",
    "census_2010_from_2002.xls": XW_BASE + "2010-occ-codes-with-crosswalk-from-2002-2011.xls",
    "acs_2006_2010_conversion.xlsx": XW_BASE + "2006-2010-acs-pums-occupation-conversion-rates.xlsx",
    "census_2018_crosswalk.xlsx": XW_BASE + "2018-occupation-code-list-and-crosswalk.xlsx",
}


def xw_file(name: str) -> Path:
    """Crosswalk source file, downloaded once into raw/acs_pums/crosswalks/."""
    XW_DIR.mkdir(parents=True, exist_ok=True)
    f = XW_DIR / name
    if not f.exists():
        r = requests.get(XW_URLS[name], timeout=180, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        f.write_bytes(r.content)
        LOG.info("downloaded crosswalk source %s", name)
    return f


def _norm_code(x) -> str | None:
    """Occupation code as a zero-padded 4-digit string (None if not a code)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    t = str(x).strip()
    if t.endswith(".0"):
        t = t[:-2]
    return t.zfill(4) if re.fullmatch(r"\d{1,4}", t) else None


# Census 2000 (three-digit) occupation code -> IPUMS `occ1990`, from Appendix A of
# Meyer and Osborne (2005), "Proposed Category System for 1960-2000 Census
# Occupations" (BLS Working Paper 383, linked from IPUMS USA's documentation of
# OCC1990: usa.ipums.org/usa/resources/chapter4/OCCBLS_paper.pdf).  The table was
# read off the PDF (`parse_bls_appendix_a`, which re-derives it with `pdftotext`)
# and checked: every Census 2000 code appears exactly once.  Two IPUMS departures
# from the BLS paper are applied (IPUMS OCC1990 documentation): parking
# enforcement workers (384) -> 423, and "helper-production workers" (895) stay 874.
# ACS 2002-vintage codes (2005-09, and 2009 in the pooled files) are the Census
# 2000 codes times ten.
BLS_2000_TO_OCC1990 = {
    1: 4, 2: 22, 3: 3, 4: 13, 5: 13, 6: 13, 10: 22, 11: 22, 12: 7, 13: 8, 14: 22, 15: 33, 16:
    373, 20: 475, 21: 473, 22: 22, 23: 14, 30: 22, 31: 17, 32: 19, 33: 21, 34: 17, 35: 15, 36:
    21, 40: 16, 41: 18, 42: 21, 43: 22, 50: 34, 51: 28, 52: 29, 53: 33, 54: 375, 56: 36, 60: 22,
    62: 27, 70: 65, 71: 26, 72: 21, 73: 37, 80: 23, 81: 254, 82: 25, 83: 25, 84: 25, 85: 25, 86:
    24, 90: 36, 91: 25, 93: 23, 94: 25, 95: 25, 100: 64, 101: 229, 102: 229, 104: 64, 106: 64,
    110: 64, 111: 64, 120: 66, 121: 68, 122: 65, 123: 67, 124: 68, 130: 43, 131: 218, 132: 44,
    133: 59, 134: 59, 135: 48, 136: 53, 140: 55, 141: 55, 142: 59, 143: 56, 144: 59, 145: 45,
    146: 57, 150: 59, 151: 59, 152: 47, 153: 59, 154: 217, 155: 214, 156: 218, 160: 77, 161: 78,
    164: 79, 165: 83, 170: 69, 171: 74, 172: 73, 174: 75, 176: 76, 180: 166, 181: 166, 182: 167,
    183: 168, 184: 173, 186: 169, 190: 223, 191: 223, 192: 224, 193: 225, 194: 235, 196: 214,
    200: 163, 201: 174, 202: 465, 204: 176, 205: 176, 206: 176, 210: 178, 211: 179, 214: 234,
    215: 234, 220: 154, 230: 155, 231: 156, 232: 157, 233: 158, 234: 159, 240: 165, 243: 164,
    244: 329, 254: 159, 255: 159, 260: 188, 263: 185, 270: 187, 271: 187, 272: 199, 274: 193,
    275: 186, 276: 194, 280: 198, 281: 195, 282: 13, 283: 195, 284: 184, 285: 183, 286: 194,
    290: 228, 291: 189, 292: 195, 296: 228, 300: 89, 301: 85, 303: 97, 304: 87, 305: 96, 306:
    84, 311: 106, 312: 88, 313: 95, 314: 104, 315: 99, 316: 103, 320: 105, 321: 105, 322: 98,
    323: 104, 324: 105, 325: 86, 326: 89, 330: 203, 331: 204, 332: 206, 340: 208, 341: 678, 350:
    207, 351: 205, 352: 677, 353: 208, 354: 208, 360: 447, 361: 99, 362: 103, 363: 469, 364:
    445, 365: 446, 370: 423, 371: 418, 372: 417, 373: 415, 374: 417, 375: 417, 380: 423, 382:
    418, 383: 423, 384: 423, 385: 418, 386: 418, 390: 427, 391: 418, 392: 426, 394: 425, 395:
    427, 400: 436, 401: 436, 402: 436, 403: 444, 404: 434, 405: 439, 406: 443, 411: 435, 412:
    443, 413: 443, 414: 444, 415: 469, 416: 444, 420: 448, 421: 485, 422: 453, 423: 405, 424:
    455, 425: 486, 430: 22, 432: 456, 434: 479, 435: 487, 440: 459, 441: 773, 442: 462, 443:
    459, 446: 469, 450: 457, 451: 458, 452: 458, 453: 464, 454: 461, 455: 463, 460: 468, 461:
    447, 462: 175, 464: 468, 465: 469, 470: 243, 471: 243, 472: 276, 474: 274, 475: 274, 476:
    275, 480: 256, 481: 253, 482: 255, 483: 318, 484: 274, 485: 274, 490: 283, 492: 254, 493:
    258, 494: 274, 495: 277, 496: 274, 500: 303, 501: 348, 502: 348, 503: 349, 510: 378, 511:
    344, 512: 337, 513: 276, 514: 338, 515: 365, 516: 383, 520: 336, 521: 326, 522: 389, 523:
    316, 524: 376, 525: 377, 526: 335, 530: 317, 531: 316, 532: 329, 533: 376, 534: 316, 535:
    326, 536: 328, 540: 319, 541: 318, 542: 336, 550: 364, 551: 357, 552: 359, 553: 366, 554:
    354, 555: 355, 556: 354, 560: 373, 561: 364, 562: 365, 563: 368, 570: 313, 580: 308, 581:
    385, 582: 315, 583: 315, 584: 375, 585: 356, 586: 379, 590: 347, 591: 384, 592: 386, 593:
    389, 600: 496, 601: 489, 602: 475, 604: 488, 605: 479, 610: 498, 611: 498, 612: 496, 613:
    496, 620: 558, 621: 643, 622: 563, 623: 567, 624: 563, 625: 588, 626: 869, 630: 594, 631:
    599, 632: 844, 633: 573, 635: 575, 636: 589, 640: 593, 642: 579, 643: 583, 644: 585, 646:
    584, 650: 597, 651: 595, 652: 596, 653: 597, 660: 866, 666: 35, 670: 543, 671: 599, 672:
    593, 673: 869, 674: 889, 675: 889, 676: 599, 680: 614, 682: 598, 683: 615, 684: 616, 691:
    617, 692: 614, 693: 869, 694: 617, 700: 503, 701: 525, 702: 527, 703: 533, 704: 577, 705:
    533, 710: 523, 711: 533, 712: 523, 713: 575, 714: 508, 715: 514, 716: 514, 720: 505, 721:
    507, 722: 516, 724: 509, 726: 516, 730: 539, 731: 534, 732: 526, 733: 518, 734: 549, 735:
    519, 736: 544, 741: 577, 742: 527, 743: 535, 751: 804, 752: 199, 754: 536, 755: 549, 756:
    549, 760: 577, 761: 865, 762: 549, 770: 628, 771: 785, 772: 785, 773: 785, 774: 597, 775:
    785, 780: 687, 781: 686, 783: 763, 784: 688, 785: 769, 790: 233, 792: 755, 793: 713, 794:
    707, 795: 706, 796: 708, 800: 709, 801: 703, 802: 703, 803: 637, 804: 766, 806: 645, 810:
    719, 812: 684, 813: 634, 814: 783, 815: 724, 816: 646, 820: 723, 821: 644, 822: 726, 823:
    679, 824: 734, 825: 736, 826: 736, 830: 748, 831: 747, 832: 744, 833: 669, 834: 745, 835:
    666, 836: 749, 840: 743, 841: 739, 842: 738, 843: 755, 844: 645, 845: 668, 846: 749, 850:
    657, 851: 658, 852: 645, 853: 727, 854: 729, 855: 733, 860: 695, 861: 696, 862: 694, 863:
    699, 864: 757, 865: 756, 871: 769, 872: 755, 873: 766, 874: 799, 875: 535, 876: 678, 880:
    754, 881: 759, 883: 774, 884: 779, 885: 753, 886: 764, 890: 779, 891: 649, 892: 675, 893:
    765, 894: 779, 895: 874, 896: 779, 900: 803, 903: 226, 904: 227, 911: 809, 912: 808, 913:
    804, 914: 809, 915: 809, 920: 824, 923: 825, 924: 823, 926: 824, 930: 829, 931: 829, 933:
    829, 934: 834, 935: 813, 936: 885, 941: 463, 942: 883, 950: 876, 951: 848, 952: 853, 956:
    848, 960: 804, 961: 887, 962: 889, 963: 878, 964: 888, 965: 859, 972: 875, 973: 859, 974:
    876, 975: 454, 980: 905, 981: 905, 982: 905, 983: 905,
}


# Departures from the BLS table found by matching the workbook.  With the sample
# rules fixed, 13 of the 21 groups reproduce the workbook's counts exactly for both
# genders on the BLS table alone.  Six departures at the Census-2000-code level
# (`BLS_DEPARTURES_2000`) and five 2010 codes (`OCC1990_OVERRIDES_2010`) close the
# other eight.  They were found by asking which small set of occupation codes
# reproduces each group's residual in FOUR dimensions at once (person counts and
# wage-sample counts, men and women); the search returns a unique solution with 11
# flips, each one semantically sensible, and afterwards all 21 groups match the
# workbook exactly in all four dimensions.  The destination GROUP is pinned down by
# the sums; the occ1990 code inside it is a placeholder.  Each departure, with the
# persons it moves (men / women, unweighted), reading "code -> group (BLS group)":
#   341  health diagnosing and treating practitioner support technicians (ACS 2002
#        code 3410, 2010 code 3420) -> Nurses, therapists, other health (Precision):
#        1,282 / 5,436
#   73   other business operations specialists (2002 code 0730; 2010 0740, 0425)
#        -> Executives (Management related): 996.5 / 1,461
#   455  transportation attendants (2002 code 4550; 2010 9050 flight attendants,
#        9415) -> Health and science technicians (Food and personal services): 176 / 422.5
#   824, 826  printing machine operators and job printers (2002 8240, 8260; 2010 8255,
#        8256) -> Precision manufacturing (Manufacturing operators): 1,145 / 312.
#        825 (prepress, 8250) stays
#   441  motion picture projectionists (2002 and 2010 code 4410) -> Food, cleaning and
#        personal services (Manufacturing operators): 70.5 / 9.5
#   2010 0726 fundraisers -> Executives (Sales; the Census crosswalk sends 25% of
#        2002 code 4960 "other sales workers" there): 173.5 / 678
#   2010 0735 market research analysts and marketing specialists -> Executives
#        (Natural and social scientists): 1,042 / 1,727.5
#   2010 2015, 2016 probation officers, social and human service assistants ->
#        Natural and social scientists (Food and personal services): 629 / 1,285
#   2010 4465 morticians, undertakers and funeral directors -> Food, cleaning and
#        personal services (Executives): 154 / 81.5
BLS_DEPARTURES_2000 = {341: 446, 73: 22, 455: 235, 824: 684, 826: 684, 441: 462}   # census-2000 code -> occ1990
OCC1990_OVERRIDES_2010 = {"0726": 13, "0735": 13, "2015": 174, "2016": 174, "4465": 469}   # 2010 codes; later vintages inherit
for _c, _j in BLS_DEPARTURES_2000.items():
    BLS_2000_TO_OCC1990[_c] = _j


def p1990_by_2002() -> dict[str, pd.Series]:
    """P(occ1990 | 2002 code): degenerate, from the table above."""
    return {f"{c * 10:04d}": pd.Series({float(j): 1.0}) for c, j in BLS_2000_TO_OCC1990.items()}


def load_2002_to_2010() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(aligned, rates): Census direct 2002->2010 code pairs (rows of the
    side-by-side crosswalk that carry both codes) and the ACS conversion rates
    for the 28 split/merged 2002 codes."""
    d = pd.read_excel(xw_file("census_2010_from_2002.xls"), sheet_name="2002to2010xwalk",
                      header=None, dtype=str).iloc[4:, :6]
    d.columns = ["soc02", "occ02", "t02", "soc10", "occ10", "t10"]
    d["occ02"] = d.occ02.map(_norm_code)
    d["occ10"] = d.occ10.map(_norm_code)
    aligned = d.dropna(subset=["occ02", "occ10"])[["occ02", "occ10"]].reset_index(drop=True)
    t = pd.read_excel(xw_file("acs_2006_2010_conversion.xlsx"), sheet_name="2. Total Conversion Rate",
                      header=None, dtype=str)
    i = t.index[t[0] == "Occupation"][0]
    r = t.iloc[i + 3:, [0, 3, 6]].copy()
    r.columns = ["occ02", "occ10", "rate"]
    r["occ02"] = r.occ02.map(_norm_code)
    r["occ10"] = r.occ10.map(_norm_code)
    r = r.dropna()
    r["rate"] = r.rate.astype(float)
    return aligned, r.reset_index(drop=True)


def p2002_by_2010(aligned: pd.DataFrame, rates: pd.DataFrame, n02: pd.Series) -> dict[str, pd.Series]:
    """P(2002 code | 2010 code) for every 2010 code in the Census crosswalk.
    A 2010 code inherits from (i) its direct 2002 partner if that code was not
    split, and (ii) every split 2002 code that sends it a share `rate`; sources
    are weighted by their 2002 employment counts n02 (P(k|c) ~ n_k * rate)."""
    split = set(rates.occ02)
    w: dict[str, dict[str, float]] = {}
    for k, c in aligned.itertuples(index=False):
        if k not in split:
            w.setdefault(c, {})[k] = 1.0
    for k, c, rate in rates.itertuples(index=False):
        w.setdefault(c, {})[k] = w.get(c, {}).get(k, 0.0) + rate
    out = {}
    for c, src in w.items():
        s = pd.Series({k: rate * float(n02.get(k, 0.0)) for k, rate in src.items()})
        if s.sum() <= 0:                        # no 2002 data: weight by rate only
            s = pd.Series(src)
        out[c] = s / s.sum()
    return out


def load_2010_to_2018() -> tuple[pd.DataFrame, set[str]]:
    """Census '2010 -> 2018 occupation code changes': (pairs [occ10, occ18] of
    2010 codes and the 2018 codes they were split into or renamed to, set of
    2010 codes that changed)."""
    d = pd.read_excel(xw_file("census_2018_crosswalk.xlsx"), sheet_name="Occ Code Changes",
                      header=None, dtype=str).iloc[3:, :3]
    d.columns = ["occ10", "t10", "occ18"]
    d["occ10"] = d.occ10.map(_norm_code).ffill()
    d["occ18"] = d.occ18.map(_norm_code)
    d = d.dropna(subset=["occ10", "occ18"])[["occ10", "occ18"]].reset_index(drop=True)
    return d, set(d.occ10)


def p2010_by_2018(pairs: pd.DataFrame, changed: set[str], n10: pd.Series, n18: pd.Series,
                  codes18: list[str]) -> dict[str, pd.Series]:
    """P(2010 code | 2018 code).  A 2018 code with the same number as a 2010 code
    that did not change comes from it; codes created by a split or rename come
    from the 2010 codes listing them, weighted by n10 * share (the share of the
    source's 2018 employment that each target holds)."""
    tot18 = pairs.merge(n18.rename("n18"), left_on="occ18", right_index=True, how="left") \
                 .fillna({"n18": 0.0}).groupby("occ10").n18.sum()
    out = {}
    for a in codes18:
        cand = {}
        if a not in changed and a in n10.index:
            cand[a] = float(n10[a]) + 1e-9
        for c in pairs.loc[pairs.occ18 == a, "occ10"]:
            share = float(n18.get(a, 0.0)) / tot18[c] if tot18.get(c, 0) > 0 else 1.0
            cand[c] = cand.get(c, 0.0) + float(n10.get(c, 0.0)) * share + 1e-9
        if a not in changed and a not in cand:
            cand[a] = 1.0                       # unchanged code, no 2010 count
        if not cand:
            continue
        s = pd.Series(cand)
        out[a] = s / s.sum()
    return out


def _occ_dict(api: Api, product: str, year: int, var: str) -> dict[str, str]:
    """Code -> label of an occupation variable from the dataset's variables.json."""
    v = dataset_variables(api, product, year)[var]["values"]["item"]
    return {_norm_code(k): lab for k, lab in v.items() if _norm_code(k)}


def count_codes(persons: pd.DataFrame) -> dict[str, pd.Series]:
    """Employed 25-34 records per occupation code and vintage (weights for the
    crosswalk's many-to-many links)."""
    e = persons[persons.esr.isin([1, 2]) & persons.occ.notna()]
    return {v: g.groupby("occ").size().astype(float) for v, g in e.groupby("vintage")}


def build_crosswalks(api: Api, counts: dict[str, pd.Series]) -> pd.DataFrame:
    """Crosswalk table [vintage, occ, label, occ1990, p_occ1990, group, group_alt,
    group_mass, ambiguous] for every valid occupation code of the three vintages.

    `group` is the group of the modal occ1990 (IPUMS convention); `group_alt` is the
    group with the largest probability mass; they differ only where the propagated
    distribution straddles groups.  `ambiguous` flags codes whose modal group holds
    less than 75% of the mass."""
    valid = valid_codes(api)
    labels = {"02": _occ_dict(api, "acs5", 2013, "OCCP02"),
              "10": {**_occ_dict(api, "acs1", 2017, "OCCP"), **_occ_dict(api, "acs5", 2013, "OCCP10")},
              "18": _occ_dict(api, "acs1", 2018, "OCCP")}
    p02 = p1990_by_2002()
    aligned, rates = load_2002_to_2010()
    p2002_2010 = p2002_by_2010(aligned, rates, counts.get("02", pd.Series(dtype=float)))
    p10: dict[str, pd.Series] = {}
    for c, ps in p2002_2010.items():
        acc = pd.Series(dtype=float)
        for k, pk in ps.items():
            if k in p02:
                acc = acc.add(p02[k] * pk, fill_value=0.0)
        if len(acc):
            p10[c] = acc / acc.sum()
    for c, j in OCC1990_OVERRIDES_2010.items():
        p10[c] = pd.Series({float(j): 1.0})
    pairs, changed = load_2010_to_2018()
    p10_18 = p2010_by_2018(pairs, changed, counts.get("10", pd.Series(dtype=float)),
                           counts.get("18", pd.Series(dtype=float)), sorted(valid["18"]))
    p18: dict[str, pd.Series] = {}
    for a, ps in p10_18.items():
        acc = pd.Series(dtype=float)
        for c, pc in ps.items():
            if c in p10:
                acc = acc.add(p10[c] * pc, fill_value=0.0)
        if len(acc):
            p18[a] = acc / acc.sum()
    dist = {"02": p02, "10": p10, "18": p18}
    rows = []
    for vint in ("02", "10", "18"):
        for code in sorted(valid[vint]):
            lab = labels[vint].get(code, "")
            if code.startswith("98"):                    # military
                rows.append((vint, code, lab, 905, 1.0, None, None, 1.0, False, "{}")); continue
            if code == "9920":                           # unemployed, never worked
                rows.append((vint, code, lab, 991, 1.0, None, None, 1.0, False, "{}")); continue
            if code < "0010":                            # not an occupation ("N/A")
                continue
            ps = dist[vint].get(code)
            if ps is None or ps.empty:
                near = sorted(dist[vint], key=lambda k: abs(int(k) - int(code)))
                if not near:
                    continue
                LOG.warning("no crosswalk for %s code %s (%s): using nearest code %s", vint, code, lab[:40], near[0])
                ps = dist[vint][near[0]]
            j = int(ps.idxmax())
            gm = ps.groupby(lambda x: occ1990_to_group(int(x)) or "unclassified").sum()
            g_alt = gm.idxmax()
            g_mod = occ1990_to_group(j) or "unclassified"
            if code in K12_CODES[vint]:                  # K-12 teachers by Census code
                g_mod = g_alt = K12
                gm = pd.Series({K12: 1.0})
            rows.append((vint, code, lab, j, float(ps.max()), g_mod, g_alt, float(gm.get(g_mod, 0.0)),
                         bool(gm.get(g_mod, 0.0) < 0.75), json.dumps({k: round(float(v), 4) for k, v in gm.items()})))
    cw = pd.DataFrame(rows, columns=["vintage", "occ", "label", "occ1990", "p_occ1990", "group",
                                     "group_alt", "group_mass", "ambiguous", "group_dist"])
    return cw


# =============================================================================
# Sample construction (HHJK definitions)
# =============================================================================

# BLS annual-average CPI-U (all urban consumers, 1982-84=100).  1999 puts wages
# in HHJK's units; 2010 defines the $1,000 income floor; the survey year is the
# dollar year of the ACS income variables after ADJINC (2013 for the 5-year file).
CPI_U = {1999: 166.6, 2009: 214.537, 2010: 218.056, 2011: 224.939, 2012: 229.594,
         2013: 232.957, 2014: 236.736, 2015: 237.017, 2016: 240.007, 2017: 245.120,
         2018: 251.107, 2019: 255.657}
INCOME_FLOOR_2010 = 1000.0
MIN_SCHL = 12            # ACS SCHL 12 = grade 9: "at least one year of high school"
MIN_HOURS_FULL, MIN_HOURS_PART = 30, 15
MIN_WEEKS = 48

# years of schooling as in HHJK's `highgrade` (create_*_IPUMS_extract.do, step C3),
# from the ACS SCHL categories through IPUMS `educ`; the Mincer workbook then
# bounds it to [9, 17] ("min 9, max 17 years of school").
def schl_to_highgrade(schl: pd.Series) -> pd.Series:
    """ACS SCHL -> HHJK highgrade (grade 9 = 9 ... 12th grade / diploma / GED = 12,
    some college = 13, associate = 14, bachelor = 16, graduate = 19)."""
    m = {12: 9, 13: 10, 14: 11, 15: 12, 16: 12, 17: 12, 18: 13, 19: 13, 20: 14,
         21: 16, 22: 19, 23: 19, 24: 19}
    out = schl.map(m)
    out[schl <= 11] = 8                       # below grade 9 (bounded to 9 below)
    out[schl <= 1] = 0
    return out.astype(float)


def _to_num(d: pd.DataFrame, cols) -> None:
    for c in cols:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")


def occ_vintage(product: str, ds_year: int, data_year: pd.Series) -> pd.Series:
    """Census occupation-code vintage ('02', '10', '18') of each record."""
    if product == "acs5":
        return data_year.map(lambda y: "02" if y <= 2009 else "10")
    return pd.Series("02" if ds_year <= 2009 else ("10" if ds_year <= 2017 else "18"),
                     index=data_year.index)


def load_person_file(name: str, api: Api | None = None) -> pd.DataFrame:
    """Harmonised person file for one sample: age 25-34, all states, with the
    columns sex, age, schl, race, esr, mil, wkhp, weeks (49/51 or NaN), wagp, semp,
    adjinc, pwgtp, occ (4-digit string), vintage, data_year, ds_year, gq."""
    spec = SAMPLES[name]
    frames = []
    for y in spec["years"]:
        files = sorted((RAW / f"{spec['product']}_{y}").glob("*.parquet"))
        if len(files) < len(STATES):
            sys.exit(f"{name}: only {len(files)} state files for {spec['product']} {y}; run `fetch`")
        d = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
        _to_num(d, ["AGEP", "SEX", "SCHL", "RAC1P", "ESR", "MIL", "WKW", "WKWN", "WKHP",
                    "WAGP", "SEMP", "ADJINC", "PWGTP", "TYPE", "TYPEHUGQ", "COW"])
        d["ds_year"] = y
        d["data_year"] = d.SERIALNO.str[:4].astype(int) if spec["product"] == "acs5" else y
        vint = occ_vintage(spec["product"], y, d.data_year)
        if spec["product"] == "acs5":
            occ = d.OCCP02.where(d.data_year <= 2009, d.OCCP10.where(d.data_year <= 2011, d.OCCP12))
        else:
            occ = d.OCCP
        d["occ"] = occ.map(_norm_code)
        d["vintage"] = vint
        gq = d["TYPE"] if "TYPE" in d else d.get("TYPEHUGQ", pd.Series(1, index=d.index))
        d["gq"] = gq.isin([2, 3])
        if "WKW" in d and d.WKW.notna().any():
            wk = d.WKW.map({1: 51.0, 2: 49.0})             # 50-52 / 48-49 weeks
        else:
            wn = d.WKWN
            wk = pd.Series(np.where(wn >= 50, 51.0, np.where(wn >= 48, 49.0, np.nan)), index=d.index)
        d["weeks"] = wk
        keep = ["SERIALNO", "ST", "ds_year", "data_year", "vintage", "occ", "gq", "weeks", "AGEP", "SEX",
                "SCHL", "RAC1P", "ESR", "MIL", "WKHP", "WAGP", "SEMP", "ADJINC", "PWGTP"]
        frames.append(d[keep])
    d = pd.concat(frames, ignore_index=True)
    d.columns = [c.lower() if c not in ("ds_year", "data_year", "vintage", "occ", "gq", "weeks") else c
                 for c in d.columns]
    d = d.rename(columns={"agep": "age"})
    return d


def valid_codes(api: Api) -> dict[str, set[str]]:
    """Valid occupation codes of each vintage (from the variables.json dictionaries)."""
    d02 = _occ_dict(api, "acs5", 2013, "OCCP02")
    d10 = set(_occ_dict(api, "acs5", 2013, "OCCP10")) | set(_occ_dict(api, "acs5", 2013, "OCCP12")) \
        | set(_occ_dict(api, "acs1", 2017, "OCCP"))
    d18 = _occ_dict(api, "acs1", 2018, "OCCP")
    return {"02": set(d02), "10": d10, "18": set(d18)}


def build_sample(d: pd.DataFrame, cw: pd.DataFrame, *, keep_gq: bool = True,
                 min_schl: int | None = MIN_SCHL, races: tuple | None = None,
                 drop_mil_last_occ: bool = False, unemp_hours_rule: bool = True,
                 group_col: str = "group", rules: str = "t7") -> pd.DataFrame:
    """Apply the HHJK sample and definitions.  `cw` is the crosswalk table with
    columns [vintage, occ, occ1990, group] (see build_crosswalks).

    Adds: group (market group of the reported occupation, K-12 carved out),
    w_mkt / w_home (1, or 0.5 for 15-29 hour workers), full, part, home,
    inc2010, wage99 (hourly, 1999 dollars) and wage_ok (the hourly-wage sample:
    employed >= 30 hours/week, >= 48 weeks, wage + business income >= $1,000 in
    2010 dollars), and years of schooling bounded to [9, 17]."""
    d = d[(d.age >= AGE_MIN) & (d.age <= AGE_MAX)].copy()
    if not keep_gq:
        d = d[~d.gq]
    if min_schl is not None:
        d = d[d.schl >= min_schl]
    if races is not None:
        d = d[d.rac1p.isin(races)]
    d = d.merge(cw[["vintage", "occ", "occ1990", group_col]].rename(columns={group_col: "group"}),
                on=["vintage", "occ"], how="left")
    mil = d.esr.isin([4, 5]) | ((d.occ1990 == 905) & drop_mil_last_occ)
    # The workbook's "drop unemployed" removes the unemployed who usually worked
    # 15+ hours/week; unemployed persons with under 15 usual hours (mostly no work
    # in the past year) stay in the sample as home producers.  This is the rule
    # that reproduces the workbook's Home Production counts exactly.
    unemp = (d.esr == 3) & (d.wkhp.fillna(0) >= MIN_HOURS_PART if unemp_hours_rule else True)
    d = d[~mil & ~unemp].copy()
    emp = d.esr.isin([1, 2])
    miss = emp & d.occ1990.isna()                       # employed, occupation missing
    LOG.info("dropping %d employed with no valid occupation code", int(miss.sum()))
    d = d[~miss].copy()
    emp = d.esr.isin([1, 2])
    d["wkhp"] = d.wkhp.fillna(0)
    d["full"] = emp & (d.wkhp >= MIN_HOURS_FULL)
    d["part"] = emp & (d.wkhp < MIN_HOURS_FULL) & (d.wkhp >= MIN_HOURS_PART)
    d["home"] = ~d.full & ~d.part
    d["w_mkt"] = np.where(d.full, 1.0, np.where(d.part, 0.5, 0.0))
    d["w_home"] = np.where(d.home, 1.0, np.where(d.part, 0.5, 0.0))
    grp = d.group.copy()
    d["group"] = grp.where(d.w_mkt > 0)
    unclassified = (d.w_mkt > 0) & d.group.isna()
    if unclassified.any():                               # HHJK: occ_broad = 0 ("Home")
        LOG.info("%d employed records with occ1990 in no HHJK range -> Home", int(unclassified.sum()))
        d.loc[unclassified, ["w_home"]] = d.loc[unclassified, "w_home"] + d.loc[unclassified, "w_mkt"]
        d.loc[unclassified, "w_mkt"] = 0.0
    # income and hourly wage
    adj = np.where(d.adjinc > 10, d.adjinc / 1e6, d.adjinc)          # the API gives 1.085467 or 1085467
    raw = d.wagp.fillna(0) + d.semp.fillna(0)                        # reference-period dollars
    hours_weeks = d.wkhp * d.weeks
    if rules == "workbook":
        # Reverse-engineered from the workbook (exact on 28 wage-sample counts and within
        # 0.2% on the 90/10 of the 14 cells whose composition is exact): the $1,000
        # floor is applied to ADJINC-adjusted income (2013 dollars in the 5-year file),
        # and hourly wages are reference-period income deflated by the CPI-U of the
        # survey year (the SERIALNO year) to 1999 dollars.
        floor_ok = raw * adj >= INCOME_FLOOR_2010
        cpi_y = d.data_year.map(CPI_U)
        d["inc2010"] = raw * adj
        d["wage99"] = raw / cpi_y * CPI_U[1999] / hours_weeks
    else:
        # ADJINC puts income in survey-year dollars; the floor is $1,000 in 2010
        # dollars and wages are in 1999 dollars (BLS annual CPI-U).
        cpi_y = d.ds_year.map(CPI_U)
        inc = raw * adj
        d["inc2010"] = inc * CPI_U[2010] / cpi_y
        floor_ok = d.inc2010 >= INCOME_FLOOR_2010
        d["wage99"] = inc / hours_weeks * CPI_U[1999] / cpi_y
    d["wage_ok"] = d.full & d.weeks.notna() & floor_ok & (d.wkhp > 0) & (d.wage99 > 0)
    d["yrs"] = schl_to_highgrade(d.schl).clip(9, 17)
    d["gender"] = np.where(d.sex == 1, "M", "F")
    return d


# =============================================================================
# Moments
# =============================================================================

def wquantile(x, w, q: float, method: str = "stata") -> float:
    """Weighted quantile.  `stata` (default): smallest x whose cumulative weight
    reaches q of the total, averaging with the next value when it is hit exactly
    (Stata's `summarize, detail` rule for frequency weights); `lower`: the same
    without averaging; `interp`: linear interpolation on cumulative weights."""
    x = np.asarray(x, float)
    w = np.asarray(w, float)
    o = np.argsort(x, kind="stable")
    x, w = x[o], w[o]
    cw = np.cumsum(w)
    t = q * cw[-1]
    i = int(np.searchsorted(cw, t - 1e-9 * cw[-1], side="left"))
    i = min(i, len(x) - 1)
    if method == "stata" and np.isclose(cw[i], t, rtol=1e-12, atol=1e-9) and i + 1 < len(x):
        return 0.5 * (x[i] + x[i + 1])
    if method == "interp":
        return float(np.interp(t, cw - 0.5 * w, x))
    return float(x[i])


def ratio_9010(x, w, method: str = "stata") -> float:
    return wquantile(x, w, 0.90, method) / wquantile(x, w, 0.10, method)


def pooled_9010(ratios: dict, weights: dict) -> float:
    """Weighted average of cell 90/10 ratios (weights normalised over the cells)."""
    keys = [k for k in ratios if k in weights and np.isfinite(ratios[k]) and weights[k] > 0]
    tw = sum(weights[k] for k in keys)
    return float(sum(ratios[k] * weights[k] for k in keys) / tw) if tw > 0 else float("nan")


def compute_moments(d: pd.DataFrame, pct_method: str = "stata") -> dict:
    """All block moments for one built sample (see `build_sample`)."""
    out: dict = {}
    for g_ in ("M", "F"):
        out[g_] = {}
    # ---- counts (person-adjusted, unweighted and PWGTP-weighted)
    for gnd in ("M", "F"):
        sub = d[d.gender == gnd]
        mk = sub[sub.w_mkt > 0]
        cnt = mk.groupby("group").w_mkt.sum()
        cntw = (mk.w_mkt * mk.pwgtp).groupby(mk.group).sum()
        home = sub.w_home.sum()
        homew = (sub.w_home * sub.pwgtp).sum()
        c, cw_ = {}, {}
        for grp in NONTEACH:
            if grp == HOME:
                c[grp], cw_[grp] = float(home), float(homew)
            else:
                c[grp], cw_[grp] = float(cnt.get(grp, 0.0)), float(cntw.get(grp, 0.0))
        c[K12], cw_[K12] = float(cnt.get(K12, 0.0)), float(cntw.get(K12, 0.0))
        tot = sum(c[g] for g in NONTEACH)
        totw = sum(cw_[g] for g in NONTEACH)
        # K-12 teachers in the workbook are excluded from the "Total" and the
        # 20 non-teaching shares; their share is over everyone (incl. home).
        out[gnd]["count"] = {**c, "Total": tot}
        out[gnd]["count_w"] = {**cw_, "Total": totw}
        dv = lambda a, b: a / b if b > 0 else float("nan")          # noqa: E731
        out[gnd]["share"] = {**{g: dv(c[g], tot) for g in NONTEACH}, K12: dv(c[K12], tot + c[K12])}
        out[gnd]["share_w"] = {**{g: dv(cw_[g], totw) for g in NONTEACH}, K12: dv(cw_[K12], totw + cw_[K12])}
        out[gnd]["n_persons"] = int(len(sub))
    # ---- wage sample (hourly wages, 1999 dollars).  The workbook's cell 90/10s are
    # unweighted (record-level; PWGTP-weighted ratios are reported too), its pooled
    # 90/10 averages the cell ratios with wage-sample record counts.
    ws = d[d.wage_ok & d.group.notna()]
    for gnd in ("M", "F"):
        sub = ws[ws.gender == gnd]
        keys = ("wage_count", "wage_count_w", "p90_p10", "p90_p10_pwgtp", "p10_wage99", "p90_wage99",
                "mean_wage99", "mean_wage99_unweighted", "mean_logwage99", "mean_logwage99_unweighted")
        o = {k: {} for k in keys}
        for grp in MARKET_GROUPS + [K12]:
            x = sub[sub.group == grp]
            o["wage_count"][grp] = float(len(x))
            o["wage_count_w"][grp] = float(x.pwgtp.sum())
            if len(x) >= 2:
                v, w1, wp = x.wage99.values, np.ones(len(x)), x.pwgtp.values
                o["p90_p10"][grp] = ratio_9010(v, w1, pct_method)
                o["p90_p10_pwgtp"][grp] = ratio_9010(v, wp, pct_method)
                o["p10_wage99"][grp] = wquantile(v, w1, 0.10, pct_method)
                o["p90_wage99"][grp] = wquantile(v, w1, 0.90, pct_method)
                o["mean_wage99"][grp] = float(np.average(v, weights=wp))
                o["mean_logwage99"][grp] = float(np.average(np.log(v), weights=wp))
                o["mean_wage99_unweighted"][grp] = float(v.mean())
                o["mean_logwage99_unweighted"][grp] = float(np.log(v).mean())
            else:
                for k in keys[2:]:
                    o[k][grp] = float("nan")
        out[gnd].update(o)
    # ---- pooled 90/10s: cell ratios averaged with wage-sample counts (non-teachers)
    def _pool(key_r: str, key_c: str) -> float:
        ratios = {(g_, gr): out[g_][key_r][gr] for g_ in ("M", "F") for gr in MARKET_GROUPS}
        wts = {(g_, gr): out[g_][key_c][gr] for g_ in ("M", "F") for gr in MARKET_GROUPS}
        return pooled_9010(ratios, wts)
    pooled = {"nonteacher": _pool("p90_p10", "wage_count"),
              "nonteacher_pwgtp": _pool("p90_p10_pwgtp", "wage_count_w")}
    kr = {g_: out[g_]["p90_p10"][K12] for g_ in ("M", "F")}
    kc = {g_: out[g_]["wage_count"][K12] for g_ in ("M", "F")}
    pooled["k12"] = pooled_9010(kr, kc)                      # count-weighted men and women
    pooled["k12_men"], pooled["k12_women"] = kr["M"], kr["F"]
    pooled["k12_male_wage_weight"] = kc["M"] / (kc["M"] + kc["F"]) if kc["M"] + kc["F"] > 0 else float("nan")
    krp = {g_: out[g_]["p90_p10_pwgtp"][K12] for g_ in ("M", "F")}
    kcp = {g_: out[g_]["wage_count_w"][K12] for g_ in ("M", "F")}
    pooled["k12_pwgtp"] = pooled_9010(krp, kcp)
    pooled["k12_men_pwgtp"], pooled["k12_women_pwgtp"] = krp["M"], krp["F"]
    x = ws[ws.group == K12]
    pooled["k12_both_genders_one_sample"] = ratio_9010(x.wage99.values, np.ones(len(x)), pct_method)
    out["pooled"] = pooled
    # ---- schooling: PWGTP x person-adjusted mean of bounded years, by group x gender
    for gnd in ("M", "F"):
        sub = d[d.gender == gnd]
        ys, yu = {}, {}
        mk = sub[sub.w_mkt > 0]
        for grp in ALL_GROUPS:
            if grp == HOME:
                w, y = sub.w_home * sub.pwgtp, sub.yrs
                wu = sub.w_home
            else:
                x = mk[mk.group == grp]
                w, y, wu = x.w_mkt * x.pwgtp, x.yrs, x.w_mkt
                if len(x) == 0:
                    ys[grp] = yu[grp] = float("nan")
                    continue
            ys[grp] = _wmean(y, w)
            yu[grp] = _wmean(y, wu)
        out[gnd]["years_school"] = ys
        out[gnd]["years_school_unweighted"] = yu
    out["schooling_comparison"] = schooling_comparison(d)
    return out


def _wmean(y, w) -> float:
    return float(np.average(y, weights=w)) if len(y) and np.sum(w) > 0 else float("nan")


def schooling_comparison(d: pd.DataFrame) -> dict:
    """Weighted (PWGTP x person adjustment) mean years of schooling of K-12
    teachers vs all non-teachers (home production included) and vs market
    non-teachers, by gender and pooled."""
    res = {}
    for label, sub in (("men", d[d.gender == "M"]), ("women", d[d.gender == "F"]), ("pooled", d)):
        t = sub[sub.group == K12]
        nt_mkt = sub[(sub.w_mkt > 0) & (sub.group != K12)]
        t_y = _wmean(t.yrs, t.w_mkt * t.pwgtp)
        mkt_y = _wmean(nt_mkt.yrs, nt_mkt.w_mkt * nt_mkt.pwgtp)
        home_w = sub.w_home * sub.pwgtp
        den = np.sum(nt_mkt.w_mkt * nt_mkt.pwgtp) + np.sum(home_w)
        all_y = (np.sum(nt_mkt.yrs * nt_mkt.w_mkt * nt_mkt.pwgtp) + np.sum(sub.yrs * home_w)) / den if den > 0 else float("nan")
        res[label] = {"k12_teachers": t_y, "all_nonteachers_incl_home": float(all_y),
                      "market_nonteachers": mkt_y, "diff_vs_all": t_y - float(all_y),
                      "diff_vs_market": t_y - mkt_y,
                      "ratio_vs_all": t_y / float(all_y), "ratio_vs_market": t_y / mkt_y}
    return res


# =============================================================================
# The workbook (validation targets)
# =============================================================================

def _find_col(ws, header_text: str, row: int, start: int = 1) -> int:
    for c in range(start, ws.max_column + 1):
        if str(ws.cell(row=row, column=c).value or "").strip() == header_text:
            return c
    raise KeyError(header_text)


def _first_block(ws, col_name: int, col_val: int, names: list[str]) -> dict[str, float]:
    """Values in column `col_val` for the first occurrence of each name in
    column `col_name` (the sheets stack a counts block on a shares block)."""
    seen: dict[str, float] = {}
    for r in range(1, ws.max_row + 1):
        nm = ws.cell(row=r, column=col_name).value
        if nm in names and nm not in seen:
            v = ws.cell(row=r, column=col_val).value
            seen[nm] = float(v) if isinstance(v, (int, float)) else float("nan")
    return seen


def read_workbook() -> dict:
    """ACS 2009-2013 columns of `wages_occ_shares_v2.xlsx` (and the 2013 schooling
    column of `Mincer_YrsSchool.xlsx`), by gender ('M', 'F') and group.  Names
    follow the `moments_shares` sheet."""
    import openpyxl
    wb = openpyxl.load_workbook(WORKBOOK, data_only=True, read_only=False)
    out: dict = {"M": {}, "F": {}}
    names = ALL_GROUPS + ["Total"]
    for sheet, key, block in (("moments_shares", "count", 0), ("occ_gender_weights", "wage_count", 0),
                              ("90_10_hr_wages_weighted", "p90_p10", 0)):
        ws = wb[sheet]
        hdr = next(r for r in range(1, 8) if str(ws.cell(row=r, column=1).value or "").startswith("Male"))
        cm = _find_col(ws, "ACS 2009-2013", hdr)
        cf = _find_col(ws, "ACS 2009-2013", hdr, cm + 1)
        out["M"][key] = _first_block(ws, 1, cm, names)
        out["F"][key] = _first_block(ws, cf - 6, cf, names)
        if sheet == "moments_shares":
            # shares block: the SECOND occurrence of each name
            for gnd, cn, cv in (("M", 1, cm), ("F", cf - 6, cf)):
                sh, seen = {}, set()
                for r in range(1, ws.max_row + 1):
                    nm = ws.cell(row=r, column=cn).value
                    if nm in names:
                        if nm in seen:
                            sh.setdefault(nm, float(ws.cell(row=r, column=cv).value))
                        seen.add(nm)
                out[gnd]["share"] = sh
    ws = wb["90_10_hr_wages_weighted"]
    cm = _find_col(ws, "ACS 2009-2013", 1)
    pooled = {}
    for r in range(1, ws.max_row + 1):
        nm = ws.cell(row=r, column=1).value
        if nm == "Weighted non-teaching":
            pooled["nonteacher"] = float(ws.cell(row=r, column=7).value)
        if nm == K12 and isinstance(ws.cell(row=r, column=1).value, str) and r > 40:
            pooled["k12"] = float(ws.cell(row=r, column=7).value)
    out["pooled"] = pooled
    # schooling, 2013 column (years), from the Mincer workbook
    wm = openpyxl.load_workbook(MINCER, data_only=True)["years_of_school"]
    sch = {"M": {}, "F": {}}
    cur = None
    for r in range(1, wm.max_row + 1):
        v = wm.cell(row=r, column=1).value
        if v in ("Male", "Female"):
            cur = "M" if v == "Male" else "F"
        elif cur and isinstance(v, str):
            key = ("Administrative Support, Clerks, Record" if v.startswith("Administrative Support")
                   else HOME if v.startswith("Home Production") else v)
            val = wm.cell(row=r, column=7).value           # 2013 years (col G)
            if key in ALL_GROUPS and isinstance(val, (int, float)):
                sch[cur][key] = float(val)
    out["M"]["years_school"], out["F"]["years_school"] = sch["M"], sch["F"]
    return out


# =============================================================================
# Provenance check of the embedded BLS table
# =============================================================================

def check_bls_table(verbose: bool = True) -> dict:
    """Re-read Appendix A of the BLS paper (`pdftotext -bbox`) and compare with
    `BLS_2000_TO_OCC1990`: every Census 2000 code must be listed exactly once, and
    assigning each cell of the 2000-codes column to the nearest standard-code row
    (by vertical position) must reproduce the embedded table except for cells the
    layout makes ambiguous (multi-line rows).  Returns counts; needs `pdftotext`."""
    import subprocess
    pdf = xw_file("OCCBLS_paper.pdf")
    html = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"], capture_output=True, text=True,
                          check=True).stdout
    pages = re.split(r"<page ", html)[1:]
    wre = re.compile(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">(.*?)</word>')

    def col(x0):                      # column of the table by x position (pt)
        return "title" if x0 < 250 else "code" if x0 < 283 else "other" if x0 < 495 else "c00"

    anchors, cells = [], []
    for pno in range(16, 28):         # Appendix A is on PDF pages 17-28
        ws = [(float(a), float(b), t) for a, b, _, _, t in wre.findall(pages[pno]) if float(b) < 730]
        if pno == 16:
            y0 = min(b for a, b, t in ws if t == "Legislators")
            ws = [w for w in ws if w[1] >= y0 - 1]
        else:
            hy = max(b for a, b, t in ws if t == "codes" and a > 490 and b < 130)
            ws = [w for w in ws if w[1] > hy + 2]
        anchors += [(pno * 1000 + b, int(t)) for a, b, t in ws if col(a) == "code" and re.fullmatch(r"\d+", t)]
        lines: dict[float, list] = {}
        for a, b, t in ws:
            if col(a) == "c00":
                lines.setdefault(round(b, 1), []).append((a, t))
        cur, ys = [], []
        for y in sorted(lines):
            txt = " ".join(t for a, t in sorted(lines[y]))
            cur.append(txt)
            ys.append(y)
            if not txt.endswith(";"):
                cells.append((pno * 1000 + float(np.mean(ys)), " ".join(cur)))
                cur, ys = [], []
    anchors.sort()
    codes_seen, mism = [], []
    for y, txt in cells:
        if not re.fullmatch(r"[\d; ]+", txt):
            continue
        cs = [int(x) for x in re.findall(r"\d+", txt)]
        codes_seen += cs
        near = min(anchors, key=lambda a: abs(a[0] - y - 3))[1]
        for c in cs:
            if c in BLS_2000_TO_OCC1990 and BLS_2000_TO_OCC1990[c] != near and c not in BLS_DEPARTURES_2000:
                mism.append((c, BLS_2000_TO_OCC1990[c], near))
    dup = sorted({c for c in codes_seen if codes_seen.count(c) > 1})
    missing = sorted(set(BLS_2000_TO_OCC1990) - set(codes_seen))
    out = {"codes_in_pdf": len(set(codes_seen)), "embedded": len(BLS_2000_TO_OCC1990),
           "duplicates": dup, "in_table_not_in_pdf": missing,
           "nearest_row_disagreements": len(mism)}
    if verbose:
        print(json.dumps(out, indent=1))
        print("nearest-row disagreements (code, embedded, nearest):", mism[:40])
    return out


# =============================================================================
# Build: all samples -> validation -> JSON and Markdown report
# =============================================================================

SHORT = {
    "Executives, Administrative, and Managerial": "Executives/managers",
    "Management Related": "Management related",
    "Architects, Engineers, Math, and Computer Science": "Architects, engineers, math, CS",
    "Natural and Social Scientists, Recreation, Religious, Arts, Athletes": "Nat./social scientists, arts, religious",
    "Doctors and Lawyers": "Doctors and lawyers",
    "Nurses, Therapists, and Other Health Service": "Nurses, therapists, health service",
    "Teachers, Postsecondary": "Teachers, postsecondary",
    "Teachers, Non-Postsecondary and Librarians": "Other teachers, counselors, librarians",
    "Health and Science Technicians": "Health and science technicians",
    "Sales, All": "Sales",
    "Administrative Support, Clerks, Record": "Administrative support",
    "Fire, Police, and Guards": "Fire, police, guards",
    "Food, Cleaning, and Personal Services and Private Household": "Food, cleaning, personal services",
    "Farm, Related Agrigulture, Logging, and Extraction": "Farm, logging, extraction",
    "Mechanics and Construction": "Mechanics and construction",
    "Precision Manufacturing": "Precision manufacturing",
    "Manufacturing Operators": "Manufacturing operators",
    "Fabricators, Inspectors, and Material Handlers": "Fabricators, inspectors, handlers",
    "Vehicle Operators": "Vehicle operators",
    HOME: "Home production",
    K12: "K-12 teachers",
}


def _clean(o):
    """JSON-safe copy (numpy scalars -> float/int, NaN -> None)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def _group_block(m: dict, key: str, groups=ALL_GROUPS) -> dict:
    return {g_: {gr: m[g_][key].get(gr) for gr in groups if gr in m[g_][key]} for g_ in ("M", "F")}


def sample_summary(m: dict, d: pd.DataFrame) -> dict:
    """JSON block for one sample (all 21 groups by gender)."""
    out = {
        "n_persons": {g_: m[g_]["n_persons"] for g_ in ("M", "F")},
        "counts": _group_block(m, "count", ALL_GROUPS + ["Total"]),
        "counts_pwgtp_weighted": _group_block(m, "count_w", ALL_GROUPS + ["Total"]),
        "shares": _group_block(m, "share"),
        "shares_pwgtp_weighted": _group_block(m, "share_w"),
        "wage_sample_counts": _group_block(m, "wage_count", MARKET_GROUPS + [K12]),
        "wage_sample_counts_pwgtp_weighted": _group_block(m, "wage_count_w", MARKET_GROUPS + [K12]),
        "p90_p10": _group_block(m, "p90_p10", MARKET_GROUPS + [K12]),
        "p90_p10_pwgtp_weighted": _group_block(m, "p90_p10_pwgtp", MARKET_GROUPS + [K12]),
        "p10_hourly_wage_1999usd": _group_block(m, "p10_wage99", MARKET_GROUPS + [K12]),
        "p90_hourly_wage_1999usd": _group_block(m, "p90_wage99", MARKET_GROUPS + [K12]),
        "pooled_p90_p10": m["pooled"],
        "mean_hourly_wage_1999usd": _group_block(m, "mean_wage99", MARKET_GROUPS + [K12]),
        "mean_log_hourly_wage_1999usd": _group_block(m, "mean_logwage99", MARKET_GROUPS + [K12]),
        "mean_hourly_wage_1999usd_unweighted": _group_block(m, "mean_wage99_unweighted", MARKET_GROUPS + [K12]),
        "mean_log_hourly_wage_1999usd_unweighted": _group_block(m, "mean_logwage99_unweighted", MARKET_GROUPS + [K12]),
        "years_schooling": _group_block(m, "years_school"),
        "years_schooling_unweighted": _group_block(m, "years_school_unweighted"),
        "schooling_k12_vs_nonteachers": m["schooling_comparison"],
    }
    # group shares by survey year (continuity of the code vintages)
    by_year = {}
    for y, dy in d.groupby("data_year"):
        mm = {}
        for gnd in ("M", "F"):
            sub = dy[dy.gender == gnd]
            mk = sub[sub.w_mkt > 0]
            cnt = {gr: float(mk.loc[mk.group == gr, "w_mkt"].sum()) for gr in MARKET_GROUPS + [K12]}
            cnt[HOME] = float(sub.w_home.sum())
            tot = sum(cnt[g] for g in NONTEACH)
            mm[gnd] = {**{g: cnt[g] / tot for g in NONTEACH}, K12: cnt[K12] / (tot + cnt[K12])}
        by_year[int(y)] = mm
    out["shares_by_data_year"] = by_year
    ws = d[d.wage_ok & d.group.notna()]
    py = {}
    for y, dy in ws.groupby("data_year"):
        ratios, wts = {}, {}
        for gnd in ("M", "F"):
            for gr in MARKET_GROUPS:
                x = dy[(dy.group == gr) & (dy.gender == gnd)]
                if len(x) > 1:
                    ratios[(gnd, gr)] = ratio_9010(x.wage99.values, np.ones(len(x)))
                    wts[(gnd, gr)] = len(x)
        py[int(y)] = pooled_9010(ratios, wts)
    out["pooled_p90_p10_nonteacher_by_data_year"] = py
    return out


def validation_tables(m: dict, w: dict) -> dict:
    """Workbook vs reproduction (acs5_2009_13), all in one dict."""
    v: dict = {"workbook": {}, "difference": {}}
    for key, wkey, groups in (("counts", "count", ALL_GROUPS + ["Total"]), ("shares", "share", ALL_GROUPS),
                              ("wage_sample_counts", "wage_count", MARKET_GROUPS + [K12]),
                              ("p90_p10", "p90_p10", MARKET_GROUPS + [K12]),
                              ("years_schooling", "years_school", ALL_GROUPS)):
        mkey = {"counts": "count", "shares": "share", "wage_sample_counts": "wage_count",
                "p90_p10": "p90_p10", "years_schooling": "years_school"}[key]
        v["workbook"][key] = {g_: {gr: w[g_][wkey].get(gr) for gr in groups} for g_ in ("M", "F")}
        v["difference"][key] = {g_: {gr: (m[g_][mkey][gr] - w[g_][wkey][gr]) if gr in m[g_][mkey] and
                                     w[g_][wkey].get(gr) is not None else None for gr in groups} for g_ in ("M", "F")}
    v["workbook"]["pooled_p90_p10"] = w["pooled"]
    v["difference"]["pooled_p90_p10"] = {"nonteacher": m["pooled"]["nonteacher"] - w["pooled"]["nonteacher"],
                                         "k12": m["pooled"]["k12"] - w["pooled"]["k12"]}
    return v


def _fmt(x, nd=2, pct=False, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "."
    if pct:
        return f"{100 * x:+.{nd}f}" if sign else f"{100 * x:.{nd}f}"
    return f"{x:+,.{nd}f}" if sign else f"{x:,.{nd}f}"


def _md(headers: list[str], rows: list[list[str]], align: str | None = None) -> str:
    a = align or ("l" + "r" * (len(headers) - 1))
    sep = "|" + "|".join(":---" if c == "l" else "---:" for c in a) + "|"
    out = ["| " + " | ".join(headers) + " |", sep]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def render_report(P: dict, moments: dict, w: dict, cw: pd.DataFrame) -> str:
    S = P["samples"]
    s5, s1, s19 = S["acs5_2009_13"], S["acs1_2009_13"], S["acs1_2016_19"]
    m5, m1, m19 = moments["acs5_2009_13"], moments["acs1_2009_13"], moments["acs1_2016_19"]
    val = s5["validation"]
    sh = SHORT
    L: list[str] = []
    ap = L.append
    mx = lambda dct: max(abs(v) for g_ in dct.values() for v in g_.values() if v is not None)   # noqa: E731

    ap("# ACS occupational block: 2009-13 reproduction and 2016-19 update\n")
    ap("Generated by `data/spatial/acs_occupations.py` from Census API PUMS microdata "
       f"({P['meta']['generated']}). Machine-readable version: `acs_occupations.json`; the crosswalk table is "
       "`acs_occupations_crosswalk.csv`. Items T4 (mean wages by occupation) and T7 (2018 occupational block) of "
       "`spatial_calibration.md`.\n")

    # ---------------- headline
    ap("## Headline\n")
    pn, pk = m5["pooled"]["nonteacher"], m5["pooled"]["k12"]
    le = np.array([np.log(m5[g_]["p90_p10"][gr] / w[g_]["p90_p10"][gr]) for g_ in "MF" for gr in MARKET_GROUPS + [K12]])
    ap(f"- **The 2009-13 block is reproduced exactly where it is a count.** All 21 group counts, the 20-group total, the "
       f"K-12 counts and all 38 wage-sample counts equal the workbook's, for both genders (largest absolute difference "
       f"{mx(val['difference']['counts']):.1f} persons and {mx(val['difference']['wage_sample_counts']):.1f} wage-sample records). "
       f"The 40 cell 90/10s have an RMSE log error of {np.sqrt((le**2).mean()):.4f} and a maximum of {np.abs(le).max():.4f}; "
       f"they agree to a median of {100*np.median([abs(v/val['workbook']['p90_p10'][g_][gr]) for g_ in ('M','F') for gr,v in val['difference']['p90_p10'][g_].items() if v is not None]):.2f}% "
       f"and at most {100*max(abs(v/val['workbook']['p90_p10'][g_][gr]) for g_ in ('M','F') for gr,v in val['difference']['p90_p10'][g_].items() if v is not None):.2f}%. "
       f"Pooled non-teacher 90/10: {pn:.4f} against the workbook's {w['pooled']['nonteacher']:.4f}; K-12 teachers {pk:.4f} against "
       f"{w['pooled']['k12']:.4f} (men {m5['pooled']['k12_men']:.4f} against {w['M']['p90_p10'][K12]:.4f}, women "
       f"{m5['pooled']['k12_women']:.4f} against {w['F']['p90_p10'][K12]:.4f}).")
    ap("- The reproduction needed several rules that the workbook `readme` does not state (section 1) and eleven occupation-code "
       "departures from the published BLS/IPUMS mapping (section 2). Each was pinned down by matching counts in four dimensions "
       "at once, not by fitting 90/10s, but they are inferences about how the workbook was built and should be read as such.")
    a_, b_ = s19["shares"], s5["shares"]
    ap(f"- **2016-19** (pooled ACS 1-year PUMS, n = {s19['n_persons']['M']:,} men and {s19['n_persons']['F']:,} women): "
       f"K-12 teaching share {100*a_['M'][K12]:.2f}% of men and {100*a_['F'][K12]:.2f}% of women (2009-13: "
       f"{100*b_['M'][K12]:.2f}% and {100*b_['F'][K12]:.2f}%); home production {100*a_['M'][HOME]:.1f}% of men and "
       f"{100*a_['F'][HOME]:.1f}% of women (2009-13: {100*b_['M'][HOME]:.1f}% and {100*b_['F'][HOME]:.1f}%); pooled non-teacher 90/10 "
       f"{m19['pooled']['nonteacher']:.3f} (same rules in 2009-13: {m1['pooled']['nonteacher']:.3f}); K-12 teacher 90/10 "
       f"{m19['pooled']['k12']:.3f} (men {m19['pooled']['k12_men']:.3f}, women {m19['pooled']['k12_women']:.3f}).\n")

    # ---------------- 1. definitions
    ap("## 1. Data and definitions\n")
    ap("**Data.** Census API, `.../data/<year>/acs/acs5/pums` (2013 5-year file, data years 2009-13) and `.../acs1/pums` "
       "(1-year 2009-13 and 2016-19), one state at a time, ages 25-34, 50 states and DC, group quarters included. The 5-year file "
       "carries the same person records as the five 1-year files (1,784,437 persons aged 25-34 in both), so the pooled 1-year "
       "sample differs from it only in weights and income dollars; it is kept as the like-for-like comparator for 2016-19. "
       "Occupation variables: `OCCP02` (2009), `OCCP10` (2010-11), `OCCP12` (2012-13) in the 5-year file; `OCCP` in the 1-year "
       "files (2002 codes in 2009, 2010 codes 2010-17, 2018 codes 2018-19). Weeks worked: bracket `WKW` (48-49 and 50-52 weeks "
       "become 49 and 51, as in HHJK's code), exact `WKWN` in 2019.\n")
    ap("**Sample and definitions** (each verified against the workbook unless marked):\n")
    ap("- Ages 25-34 with `SCHL` >= 12 (at least grade 9). Armed forces (`ESR` 4, 5) are dropped. The unemployed (`ESR` 3) are "
       "dropped **only if they usually worked 15+ hours a week; unemployed persons under 15 hours (mostly no work in the past year) "
       "stay in the sample as home producers**: this rule (30,429 men and 32,013 women) is what reproduces the workbook's Home "
       "Production counts exactly. Persons out of the labor force whose last job was military stay in.")
    ap("- Home production: out of the labor force, or employed under 15 usual hours (weight 1); employed 15-29 hours: 0.5 in home "
       "production and 0.5 in the reported occupation. Counts are **unweighted persons** with that 0.5 rule.")
    ap("- 21 groups: HHJK's `occ_broad` (19 market occupations, from occ1990 ranges in `create_2010_IPUMS_extract.do`), home "
       "production, and K-12 teachers carved out of \"Teachers, non-postsecondary and librarians\" **by Census occupation code** "
       "(2300-2340: pre-K, elementary and middle, secondary, special education and \"other teachers and instructors\"). The "
       "remainder (teacher assistants, other education workers, counselors, librarians, archivists) stays in the market group. "
       "Using occ1990 155-158 or 155-159 instead misses the workbook by 2,224 or 1,353 men.")
    ap("- Wage sample: employed >= 30 usual hours/week, 48+ weeks, wage-and-salary plus business/farm income above the floor. "
       "*Workbook rule* (2009-13 5-year file): income x ADJINC >= $1,000, hourly wage = reference-period income deflated by the "
       "CPI-U of the survey year, in 1999 dollars. This reproduces all 38 wage-sample counts exactly; a strict $1,000 in 2010 "
       "dollars floor (= $1,068 in the file's 2013 dollars) drops 35 records and is shown as a sensitivity. *T7 rule* (1-year "
       "files): ADJINC to survey-year dollars, then CPI-U (BLS annual averages, 1999 = 166.6, 2010 = 218.056) to 2010 dollars for the "
       "floor and 1999 dollars for wages.")
    ap("- Cell 90/10s are **unweighted** within occupation x gender (weighting by `PWGTP` moves them by up to 6%, and the unweighted "
       "ratios match the workbook to 0.2%), with Stata's weighted-percentile rule; the pooled non-teacher 90/10 averages cell ratios "
       "with wage-sample record counts (the workbook's `occ_gender_weights`). PWGTP-weighted versions are in the JSON.")
    ap("- Means (wages, log wages, schooling) are `PWGTP`-weighted. Schooling is HHJK's `highgrade` (grade 9 = 9 ... 12th grade, "
       "diploma or GED = 12, some college 13, associate 14, bachelor 16, graduate 19) bounded to [9, 17], as in "
       "`Mincer_YrsSchool.xlsx`.\n")

    # ---------------- 2. crosswalk
    ap("## 2. Crosswalk to occ1990 and the 21 groups\n")
    n = P["meta"]["crosswalk"]["codes"]
    ap(f"IPUMS builds `occ1990` from Census double-coded crosswalks by modal assignment and does not publish an OCC-to-OCC1990 "
       f"table for the 2010+ codes, so the chain is rebuilt ({n['02']} codes of the 2002 vintage, {n['10']} of the 2010/2012 vintage, "
       f"{n['18']} of the 2018 vintage; table in `acs_occupations_crosswalk.csv`):\n")
    ap("1. **Census 2000 code -> occ1990:** Appendix A of Meyer and Osborne (2005), BLS Working Paper 383, the table IPUMS documents "
       "as the basis of OCC1990 (`usa.ipums.org/usa/resources/chapter4/OCCBLS_paper.pdf`), read from the PDF and embedded in the "
       "script (`BLS_2000_TO_OCC1990`; `--check-bls-table` re-reads the PDF: 511 codes, none twice) with IPUMS's two documented "
       "departures. ACS 2002-vintage codes are the 2000 codes times ten.")
    ap("2. **2010 code -> 2002 code:** the Census 2002-to-2010 code crosswalk (direct matches) plus the Census ACS conversion "
       "rates for the 28 split codes, weighted by 2002-code counts. The 2012+ PUMS lists are subsets of the 2010 list.")
    ap("3. **2018 code -> 2010 code:** the Census 2010-to-2018 change list (splits and renamings), weighted by counts; unchanged "
       "codes map to themselves. One 2018 code with no 2010 counterpart (6540, solar installers) borrows the nearest code.\n")
    ap("**Where the chain departs from the workbook.** On the BLS table alone 13 of 21 groups match exactly. The rest differ by "
       "occupation codes that IPUMS classified differently from BLS (or from the crosswalk chain). Eleven code-level departures "
       "close all eight remaining groups, exactly in person counts and wage-sample counts for both genders (four numbers per "
       "group). They were found by searching for small sets of codes whose counts sum exactly to a group's surplus or deficit "
       "(single codes or short lists for most; a minimal-flip integer program with six flips for the last four groups). "
       "Uniqueness is not claimed, but each flip is a plausible reclassification and none was chosen to fit a 90/10:\n")
    dep = [
        ("341 (2002 code 3410, 2010 code 3420)", "health diagnosing and treating practitioner support technicians", "Nurses, therapists, health service", "Precision manufacturing", "1,282 / 5,436"),
        ("73 (2002 0730; 2010 0740, 0425)", "other business operations specialists; emergency management directors", "Executives", "Management related", "996.5 / 1,461"),
        ("455 (2002 4550; 2010 9050, 9415)", "transportation attendants incl. flight attendants", "Health and science technicians", "Food, cleaning, personal services", "176 / 422.5"),
        ("824, 826 (2002 8240, 8260; 2010 8255, 8256)", "printing machine operators, job printers, print binding", "Precision manufacturing", "Manufacturing operators", "1,145 / 312"),
        ("441 (2002 and 2010 4410)", "motion picture projectionists", "Food, cleaning, personal services", "Manufacturing operators", "70.5 / 9.5"),
        ("2010 0726", "fundraisers", "Executives", "Sales", "173.5 / 678"),
        ("2010 0735", "market research analysts and marketing specialists", "Executives", "Nat./social scientists", "1,042 / 1,727.5"),
        ("2010 2015, 2016", "probation officers; social and human service assistants", "Nat./social scientists", "Food, cleaning, personal services", "629 / 1,285"),
        ("2010 4465", "morticians, undertakers, funeral directors", "Food, cleaning, personal services", "Executives", "154 / 81.5"),
    ]
    ap(_md(["Census code", "Occupation", "Workbook group", "BLS chain group", "Men / women moved"],
           [list(r) for r in dep], "lllll"))
    ap("\nThese apply to every vintage (2018 codes inherit them through the 2010-code step). The occ1990 code inside the destination "
       "group is a placeholder (`BLS_DEPARTURES_2000`, `OCC1990_OVERRIDES_2010` in the script); only the group matters for the block. "
       f"Codes whose propagated distribution puts under 75% of its mass in one group are flagged `ambiguous` in the crosswalk "
       f"table ({P['meta']['crosswalk']['ambiguous_codes']} codes, {100*P['meta']['crosswalk']['share_of_2016_19_employed_in_ambiguous_codes']:.1f}% of 2016-19 employed persons).\n")

    # ---------------- 3. validation
    ap("## 3. Validation against the workbook (ACS 2009-13, 5-year PUMS)\n")
    ap("### 3.1 Counts by group and gender (`moments_shares`, \"ACS 2009-2013\")\n")
    rows = []
    for g in ALL_GROUPS:
        rows.append([sh[g], _fmt(w["M"]["count"][g], 1), _fmt(m5["M"]["count"][g], 1), _fmt(w["F"]["count"][g], 1),
                     _fmt(m5["F"]["count"][g], 1), _fmt(max(abs(val["difference"]["counts"][x][g]) for x in "MF"), 1)])
    rows.append(["**Total, 20 non-teaching groups**", _fmt(w["M"]["count"]["Total"], 1), _fmt(m5["M"]["count"]["Total"], 1),
                 _fmt(w["F"]["count"]["Total"], 1), _fmt(m5["F"]["count"]["Total"], 1),
                 _fmt(max(abs(val["difference"]["counts"][x]["Total"]) for x in "MF"), 1)])
    ap(_md(["Group", "Men workbook", "Men here", "Women workbook", "Women here", "Max abs diff"], rows))
    ap("\nShares (over the 20 non-teaching groups; K-12 over everyone) are therefore identical: "
       f"max absolute difference {100*mx(val['difference']['shares']):.6f} percentage points.\n")
    ap("### 3.2 Wage-sample counts (`occ_gender_weights`)\n")
    rows = []
    for g in MARKET_GROUPS + [K12]:
        rows.append([sh[g], _fmt(w["M"]["wage_count"][g], 0), _fmt(m5["M"]["wage_count"][g], 0), _fmt(w["F"]["wage_count"][g], 0),
                     _fmt(m5["F"]["wage_count"][g], 0)])
    ap(_md(["Group", "Men workbook", "Men here", "Women workbook", "Women here"], rows))
    ap("\nAll 40 counts are equal.\n")
    ap("### 3.3 Within-cell 90/10 of hourly wages (`90_10_hr_wages_weighted`)\n")
    rows = []
    for g in MARKET_GROUPS + [K12]:
        rows.append([sh[g], _fmt(w["M"]["p90_p10"][g], 3), _fmt(m5["M"]["p90_p10"][g], 3),
                     _fmt(m5["M"]["p90_p10"][g] / w["M"]["p90_p10"][g] - 1, 2, pct=True, sign=True),
                     _fmt(w["F"]["p90_p10"][g], 3), _fmt(m5["F"]["p90_p10"][g], 3),
                     _fmt(m5["F"]["p90_p10"][g] / w["F"]["p90_p10"][g] - 1, 2, pct=True, sign=True)])
    ap(_md(["Group", "Men workbook", "Men here", "diff %", "Women workbook", "Women here", "diff %"], rows))
    ap(f"\nRMSE of the log ratio over the 40 cells {np.sqrt((le**2).mean()):.4f}, maximum {np.abs(le).max():.4f} (doctors and lawyers, men); "
       f"{int((np.abs(le) < 1e-3).sum())} of 40 cells are within 0.1%. The remaining differences come from how hourly wages are built "
       "(the workbook's weeks-worked convention or deflator rounding, not identified); they are as large in cells whose composition "
       "is exact, so they do not reflect crosswalk error.\n")
    ap("### 3.4 Pooled 90/10s\n")
    mm = m5["pooled"]
    ap(_md(["Moment", "Workbook", "Here", "Difference"],
           [["Weighted non-teacher (19 occupations x gender, wage-sample counts)", f"{w['pooled']['nonteacher']:.4f}", f"{mm['nonteacher']:.4f}", f"{mm['nonteacher']-w['pooled']['nonteacher']:+.4f}"],
            ["K-12 teachers (24% men, 76% women)", f"{w['pooled']['k12']:.4f}", f"{mm['k12']:.4f}", f"{mm['k12']-w['pooled']['k12']:+.4f}"],
            ["K-12 teachers, men", f"{w['M']['p90_p10'][K12]:.4f}", f"{mm['k12_men']:.4f}", f"{mm['k12_men']-w['M']['p90_p10'][K12]:+.4f}"],
            ["K-12 teachers, women", f"{w['F']['p90_p10'][K12]:.4f}", f"{mm['k12_women']:.4f}", f"{mm['k12_women']-w['F']['p90_p10'][K12]:+.4f}"]]))
    ap("")
    ap("### 3.5 Sensitivity of the pooled 90/10s\n")
    st = P["meta"]["sensitivity"]["acs5_strict_2010_dollar_floor"]
    ap(_md(["Variant", "Non-teacher", "K-12 (pooled)", "K-12 men", "K-12 women"],
           [["Workbook rules, unweighted cells (baseline)", f"{mm['nonteacher']:.4f}", f"{mm['k12']:.4f}", f"{mm['k12_men']:.4f}", f"{mm['k12_women']:.4f}"],
            ["Cell ratios weighted by PWGTP", f"{mm['nonteacher_pwgtp']:.4f}", f"{mm['k12_pwgtp']:.4f}", f"{mm['k12_men_pwgtp']:.4f}", f"{mm['k12_women_pwgtp']:.4f}"],
            ["Strict $1,000 in 2010 dollars floor", f"{st['pooled_p90_p10']['nonteacher']:.4f}", f"{st['pooled_p90_p10']['k12']:.4f}", f"{st['pooled_p90_p10']['k12_men']:.4f}", f"{st['pooled_p90_p10']['k12_women']:.4f}"],
            ["Pooled 1-year files, T7 rules (same persons)", f"{m1['pooled']['nonteacher']:.4f}", f"{m1['pooled']['k12']:.4f}", f"{m1['pooled']['k12_men']:.4f}", f"{m1['pooled']['k12_women']:.4f}"]]))
    ap("\nThe strict floor drops "
       f"{sum(m5[g_]['wage_count'][gr] for g_ in 'MF' for gr in MARKET_GROUPS) - sum(st['wage_sample_total'].values()):.0f} of "
       f"{sum(m5[g_]['wage_count'][gr] for g_ in 'MF' for gr in MARKET_GROUPS):,.0f} non-teacher wage-sample records. "
       "The unweighted-versus-PWGTP choice moves the pooled non-teacher ratio by about 0.03 and is the only material choice; "
       "the workbook's numbers are the unweighted ones.\n")
    ap("### 3.6 Schooling: `Mincer_YrsSchool.xlsx`, sheet `years_of_school`, 2013 column\n")
    rows = []
    dd = []
    for g in ALL_GROUPS:
        r = [sh[g]]
        for x in "MF":
            a2, b2 = m5[x]["years_school"][g], w[x]["years_school"][g]
            r += [_fmt(b2, 3), _fmt(a2, 3), _fmt(a2 - b2, 3, sign=True)]
            dd.append(a2 - b2)
        rows.append(r)
    ap(_md(["Group", "Men workbook", "Men here", "diff", "Women workbook", "Women here", "diff"], rows))
    dd = np.array(dd)
    ap(f"\nThe mean years of schooling are close but not identical (RMSE {np.sqrt((dd**2).mean()):.3f} years, mean difference "
       f"{dd.mean():+.3f}, largest {np.abs(dd).max():.2f} for men's home production, which the workbook puts at "
       f"{w['M']['years_school'][HOME]:.2f} against {m5['M']['years_school'][HOME]:.2f} here). The `Mincer` workbook's sample and "
       "education recode are not documented; the counts workbook could be matched exactly, this one could not, and the gap is "
       "unexplained (it is small next to the 2-3 year gap between K-12 teachers and non-teachers).\n")

    # ---------------- 4. T4
    ap("## 4. T4: mean wages and schooling by occupation and gender\n")
    for nm, sm, lab in (("acs5_2009_13", s5, "ACS 2009-13 (5-year file, workbook rules)"), ("acs1_2016_19", s19, "ACS 2016-19 (T7 rules)")):
        ap(f"### 4.{'1' if nm.startswith('acs5') else '2'} Mean hourly wage and mean log hourly wage, 1999 dollars: {lab}\n")
        rows = []
        for g in MARKET_GROUPS + [K12]:
            r = [sh[g]]
            for x in "MF":
                r += [_fmt(sm["wage_sample_counts"][x][g], 0), _fmt(sm["mean_hourly_wage_1999usd"][x][g], 2), _fmt(sm["mean_log_hourly_wage_1999usd"][x][g], 3)]
            rows.append(r)
        ap(_md(["Group", "Men N", "Mean wage", "Mean log wage", "Women N", "Mean wage", "Mean log wage"], rows))
        ap("\n`PWGTP`-weighted on the wage sample; unweighted means are in the JSON. Hourly wages in 1999 dollars.\n")
    ap("### 4.3 Mean years of schooling (bounded to [9, 17])\n")
    rows = []
    for g in ALL_GROUPS:
        rows.append([sh[g]] + [_fmt(x[y]["years_schooling"][gg][g], 2) for gg in "MF" for x, y in ((S, "acs5_2009_13"), (S, "acs1_2016_19"))])
    ap(_md(["Group", "Men 2009-13", "Men 2016-19", "Women 2009-13", "Women 2016-19"], rows))
    ap("\n### 4.4 K-12 teachers against non-teachers\n")
    rows = []
    for nm, lab in (("acs5_2009_13", "2009-13"), ("acs1_2016_19", "2016-19")):
        c = S[nm]["schooling_k12_vs_nonteachers"]
        for x, xl in (("men", "Men"), ("women", "Women"), ("pooled", "Both")):
            r = c[x]
            rows.append([lab, xl, _fmt(r["k12_teachers"], 2), _fmt(r["all_nonteachers_incl_home"], 2), _fmt(r["market_nonteachers"], 2),
                         _fmt(r["ratio_vs_all"], 3), _fmt(r["ratio_vs_market"], 3)])
    ap(_md(["Sample", "Gender", "K-12 teachers", "All non-teachers (incl. home)", "Market non-teachers", "Ratio to all", "Ratio to market"], rows))
    ap("\nMean years, `PWGTP` x 0.5-rule weights. Dividing years by 25 gives the model's schooling shares; the ratio column is the "
       "`known miss` of the calibration note (model 0.91, data about 1.2).\n")

    # ---------------- 5. T7
    ap("## 5. T7: the 2016-19 block\n")
    ap("### 5.1 Shares and counts\n")
    rows = []
    for g in ALL_GROUPS:
        r = [sh[g]]
        for x in "MF":
            r += [_fmt(s5["shares"][x][g], 2, pct=True), _fmt(s19["shares"][x][g], 2, pct=True),
                  _fmt(s19["shares"][x][g] - s5["shares"][x][g], 2, pct=True, sign=True)]
        rows.append(r)
    ap(_md(["Group (share, %)", "Men 09-13", "Men 16-19", "change (pp)", "Women 09-13", "Women 16-19", "change (pp)"], rows))
    ap("\nShares over the 20 non-teaching groups; the K-12 row is over everyone including home production. "
       "The 2009-13 columns are the workbook's.\n")
    ch = P["changes"]["shares_2016_19_minus_2009_13"]
    big = sorted(((abs(v), g_, gr, v) for g_ in "MF" for gr, v in ch[g_].items()), reverse=True)[:8]
    ap("Largest share changes (percentage points): " + "; ".join(f"{'men' if g_=='M' else 'women'} {sh[gr]} {100*v:+.2f}" for _, g_, gr, v in big) + ".\n")
    ap("Counts and wage-sample counts, 2016-19 (unweighted with the 0.5 rule; weighted = PWGTP, millions of persons):\n")
    rows = []
    for g in ALL_GROUPS:
        rows.append([sh[g]] + [_fmt(s19["counts"][x][g], 1) for x in "MF"] + [_fmt(s19["counts_pwgtp_weighted"][x][g] / 1e6, 3) for x in "MF"] +
                    ([_fmt(s19["wage_sample_counts"]["M"][g], 0), _fmt(s19["wage_sample_counts"]["F"][g], 0)] if g not in (HOME,) else [".", "."]))
    rows.append(["**Total (20 non-teaching groups)**"] + [_fmt(s19["counts"][x]["Total"], 1) for x in "MF"] +
                [_fmt(s19["counts_pwgtp_weighted"][x]["Total"] / 1e6, 3) for x in "MF"] +
                [_fmt(sum(s19["wage_sample_counts"]["M"][g] for g in MARKET_GROUPS), 0), _fmt(sum(s19["wage_sample_counts"]["F"][g] for g in MARKET_GROUPS), 0)])
    ap(_md(["Group", "Men count", "Women count", "Men weighted (millions)", "Women weighted (millions)", "Men wage sample", "Women wage sample"], rows))
    ap("")
    ap("### 5.2 Within-cell 90/10s and the pooled ratios\n")
    rows = []
    for g in MARKET_GROUPS + [K12]:
        r = [sh[g]]
        for x in "MF":
            r += [_fmt(s5["validation"]["workbook"]["p90_p10"][x][g], 3), _fmt(s1["p90_p10"][x][g], 3), _fmt(s19["p90_p10"][x][g], 3)]
        rows.append(r)
    ap(_md(["Group", "Men wb 09-13", "Men 09-13 same rules", "Men 16-19", "Women wb 09-13", "Women 09-13 same rules", "Women 16-19"], rows))
    ap("\n\"Same rules\" is the pooled 2009-13 1-year sample under the T7 income rules, the like-for-like comparator for 2016-19; "
       "it differs from the workbook column by the income-rule change only.\n")
    p19, p1 = m19["pooled"], m1["pooled"]
    ap(_md(["Pooled 90/10", "Workbook 2009-13", "2009-13 same rules", "2016-19", "2016-19 PWGTP-weighted"],
           [["Weighted non-teacher", f"{w['pooled']['nonteacher']:.4f}", f"{p1['nonteacher']:.4f}", f"{p19['nonteacher']:.4f}", f"{p19['nonteacher_pwgtp']:.4f}"],
            ["K-12 teachers", f"{w['pooled']['k12']:.4f}", f"{p1['k12']:.4f}", f"{p19['k12']:.4f}", f"{p19['k12_pwgtp']:.4f}"],
            ["K-12 teachers, men", f"{w['M']['p90_p10'][K12]:.4f}", f"{p1['k12_men']:.4f}", f"{p19['k12_men']:.4f}", f"{p19['k12_men_pwgtp']:.4f}"],
            ["K-12 teachers, women", f"{w['F']['p90_p10'][K12]:.4f}", f"{p1['k12_women']:.4f}", f"{p19['k12_women']:.4f}", f"{p19['k12_women_pwgtp']:.4f}"],
            ["Men's share of K-12 wage weight", f"{w['M']['wage_count'][K12]/(w['M']['wage_count'][K12]+w['F']['wage_count'][K12]):.4f}", f"{p1['k12_male_wage_weight']:.4f}", f"{p19['k12_male_wage_weight']:.4f}", "."]]))
    ap("")
    ap("### 5.3 Continuity across the occupation-code change (2017 to 2018)\n")
    by = s19["shares_by_data_year"]
    jumps = []
    for x in "MF":
        for g in ALL_GROUPS:
            v = [by[str(y)][x][g] if str(y) in by else by[y][x][g] for y in (2016, 2017, 2018, 2019)]
            jumps.append((abs(v[2] - v[1]), x, g, v))
    jumps.sort(reverse=True)
    typ = np.median([abs(j[3][1] - j[3][0]) for j in jumps] + [abs(j[3][3] - j[3][2]) for j in jumps])
    rows = [[("men" if x == "M" else "women"), sh[g]] + [_fmt(vv, 2, pct=True) for vv in v] + [_fmt(v[2] - v[1], 2, pct=True, sign=True)] for _, x, g, v in jumps[:8]]
    ap(_md(["Gender", "Group", "2016", "2017 (2010 codes)", "2018 (2018 codes)", "2019", "2017 to 2018 (pp)"], rows))
    ap(f"\nThe eight largest 2017-to-2018 moves in group shares (percentage points; the median within-vintage year-to-year move is "
       f"{100*typ:.3f} pp). A break at 2018 here would point to the 2018-to-2010 code step; the K-12 share moves "
       f"{100*(by['2018']['F'][K12]-by['2017']['F'][K12]) if '2018' in by else 100*(by[2018]['F'][K12]-by[2017]['F'][K12]):+.2f} pp for women "
       f"and {100*(by['2018']['M'][K12]-by['2017']['M'][K12]) if '2018' in by else 100*(by[2018]['M'][K12]-by[2017]['M'][K12]):+.2f} pp for men.\n")
    pyr = s19["pooled_p90_p10_nonteacher_by_data_year"]
    pyr = {int(k): v for k, v in pyr.items()}
    ap("Pooled non-teacher 90/10 by year: " + ", ".join(f"{y} {pyr[y]:.3f}" for y in sorted(pyr)) +
       " (2016-17 use 2010 codes, 2018-19 the 2018 codes): no break at the code change. One reclassification is real and is not "
       "undone: the 2018 codes move graduate teaching assistants out of postsecondary teachers (2200) into \"teaching assistants\" "
       "(2545, which stays with teacher assistants in the other-teachers group), so postsecondary teachers lose "
       f"{100*(by['2017']['M'][MARKET_GROUPS[6]]-by['2018']['M'][MARKET_GROUPS[6]]) if '2018' in by else 100*(by[2017]['M'][MARKET_GROUPS[6]]-by[2018]['M'][MARKET_GROUPS[6]]):.2f} pp of men "
       f"and {100*(by['2017']['F'][MARKET_GROUPS[6]]-by['2018']['F'][MARKET_GROUPS[6]]) if '2018' in by else 100*(by[2017]['F'][MARKET_GROUPS[6]]-by[2018]['F'][MARKET_GROUPS[6]]):.2f} pp of women from 2017 to 2018.\n")
    ap("### 5.4 Reading for the calibration\n")
    ap(f"- Teaching margin targets (Table 3): K-12 shares {100*a_['M'][K12]:.2f}% (men) and {100*a_['F'][K12]:.2f}% (women) in 2016-19 against "
       f"{100*b_['M'][K12]:.2f}% and {100*b_['F'][K12]:.2f}% in 2009-13; teacher 90/10 {p19['k12']:.3f} against the workbook's {w['pooled']['k12']:.3f} "
       f"({p1['k12']:.3f} on the same rules).")
    ap(f"- Home production falls for men ({100*b_['M'][HOME]:.1f}% to {100*a_['M'][HOME]:.1f}%) and women ({100*b_['F'][HOME]:.1f}% to {100*a_['F'][HOME]:.1f}%), "
       "the recession effect that caveat 3 of Table 2 flags.")
    ap(f"- The pooled non-teacher 90/10, the $\\sigma_\\epsilon$ target, is {p19['nonteacher']:.3f} in 2016-19 against {p1['nonteacher']:.3f} on the same "
       f"2009-13 rules (workbook {w['pooled']['nonteacher']:.3f}): wage dispersion keeps rising. Refitting $\\sigma_\\epsilon$ to 3.83 instead of 3.74 is "
       "a T7 step for the Julia side; this report gives only the moments.\n")

    # ---------------- 6. caveats
    ap("## 6. Caveats\n")
    ap("1. The workbook was not built from these files. Everything above is matched to it in what it reports (counts, wage-sample "
       "counts, 90/10s, schooling), and the rules that make it match are inferences (the unemployed-under-15-hours rule, the K-12 "
       "definition by Census code, the ADJINC floor, unweighted 90/10s, eleven code departures). If the 2016-19 block is used as a "
       "target, note that the departures were fitted on 2009-13 and are assumed to hold for later years and for the 2018 codes.")
    ap("2. The 2018-code step (2018 to 2010 to 2002 to 2000 codes) rests on Census change lists with count-weighted many-to-many links; "
       "section 5.3 checks it for a break in group shares.")
    ap("3. Two income conventions coexist: the workbook rule (used for the 5-year reproduction) and the instructed T7 rule (2016-19 and the "
       "pooled 2009-13 comparator). They differ by 35 wage-sample records of 515,000 and by 0.001 in the pooled 90/10 in 2009-13.")
    ap("4. 90/10s of thin cells (women in mechanics, fire and police, vehicle operators, precision manufacturing; a few hundred to a few "
       "thousand records) are noisy; the 2016-19 sample is about 10% smaller than the 2009-13 one.")
    ap("5. Group quarters are in the sample (as in the workbook); institutional inmates enter as out of the labor force or, if "
       "employed, in the occupation reported.")
    ap("6. Standard errors are not computed. Counts are unweighted person counts, PWGTP-weighted counts are in the JSON.\n")
    return "\n".join(L)


def cmd_build(out_json: Path | None = None, out_md: Path | None = None) -> None:
    api = Api()
    names = list(SAMPLES)
    LOG.info("loading person files ...")
    persons = {n: load_person_file(n) for n in names}
    counts = count_codes(pd.concat(persons.values(), ignore_index=True))
    LOG.info("building crosswalks ...")
    cw = build_crosswalks(api, counts)
    w = read_workbook()
    built, moments, summaries = {}, {}, {}
    for n in names:
        LOG.info("sample %s", n)
        d = build_sample(persons[n], cw, rules=SAMPLES[n]["rules"])
        built[n] = d
        moments[n] = compute_moments(d)
        summaries[n] = sample_summary(moments[n], d)
    # sensitivity of the 2009-13 reproduction to the income rules and the pooled-1-year weights
    sens = {}
    d_strict = build_sample(persons["acs5_2009_13"], cw, rules="t7")
    m_strict = compute_moments(d_strict)
    sens["acs5_strict_2010_dollar_floor"] = {
        "wage_sample_total": {g_: sum(m_strict[g_]["wage_count"][gr] for gr in MARKET_GROUPS) for g_ in ("M", "F")},
        "pooled_p90_p10": m_strict["pooled"]}
    summaries["acs5_2009_13"]["validation"] = validation_tables(moments["acs5_2009_13"], w)
    # 2016-19 vs 2009-13 changes
    a, b = summaries["acs1_2016_19"], summaries["acs5_2009_13"]
    a1 = summaries["acs1_2009_13"]
    changes = {
        "shares_2016_19_minus_2009_13": {g_: {gr: a["shares"][g_][gr] - b["shares"][g_][gr] for gr in ALL_GROUPS}
                                         for g_ in ("M", "F")},
        "p90_p10_2016_19_minus_2009_13_same_rules": {g_: {gr: a["p90_p10"][g_][gr] - a1["p90_p10"][g_][gr]
                                                          for gr in MARKET_GROUPS + [K12]} for g_ in ("M", "F")},
        "pooled_p90_p10_2016_19_minus_2009_13_same_rules": {
            k: a["pooled_p90_p10"][k] - a1["pooled_p90_p10"][k] for k in ("nonteacher", "k12", "k12_men", "k12_women")},
    }
    amb = cw[cw.ambiguous]
    e = persons["acs1_2016_19"]
    e = e[e.esr.isin([1, 2]) & e.occ.notna() & (e.age.between(AGE_MIN, AGE_MAX)) & (e.schl >= MIN_SCHL)]
    ambshare = float(e.merge(amb[["vintage", "occ"]], on=["vintage", "occ"]).shape[0] / len(e))
    meta = {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "script": "data/spatial/acs_occupations.py",
        "groups": {"market": MARKET_GROUPS, "home_production": HOME, "k12_teachers": K12,
                   "non_teaching_20": NONTEACH, "model_21": ALL_GROUPS},
        "samples": {n: {k: v for k, v in SAMPLES[n].items()} for n in names},
        "definitions": {
            "population": "ages 25-34, ACS PUMS via the Census API (50 states + DC, group quarters included), "
                          "SCHL >= 12 (at least grade 9), armed forces (ESR 4, 5) dropped, unemployed (ESR 3) "
                          "with 15+ usual hours dropped (unemployed with <15 hours stay as home production)",
            "home_production": "not in the labor force, or employed under 15 usual hours/week (weight 1); "
                               "employed 15-29 hours: weight 0.5 in home production and 0.5 in the occupation",
            "counts": "unweighted person counts with the 0.5 rule (the workbook's convention); "
                      "*_pwgtp_weighted use PWGTP x the 0.5 rule",
            "shares": "over the 20 non-teaching groups (19 market occupations + home production); "
                      "K-12 share is over the 20 groups plus K-12 teachers",
            "wage_sample": "employed >= 30 usual hours/week, 48+ weeks, wage+business income floor "
                           "(workbook rule: ADJINC-adjusted >= $1,000; t7 rule: $1,000 in 2010 dollars via CPI-U)",
            "hourly_wage": "income / (usual hours x 49 or 51 weeks), 1999 dollars",
            "p90_p10": "unweighted within occupation x gender (the workbook's convention); PWGTP-weighted alternative given",
            "pooled_p90_p10": "cell ratios averaged with wage-sample record counts over the 19 market occupations x gender "
                              "(K-12 teachers: men and women)",
            "mean_wage": "PWGTP-weighted mean of hourly wage (1999 dollars) and of its log on the wage sample; unweighted alternative given",
            "years_schooling": "HHJK highgrade from SCHL, bounded to [9, 17]; PWGTP x 0.5-rule weighted (all persons for home production)",
            "k12_teachers": "Census occupation codes 2300-2340 (2018 codes 2300-2330, 2350, 2360)",
            "floor_note": "T7 rule = ADJINC to survey-year dollars, then CPI-U to 2010 dollars for the floor and to 1999 dollars for wages",
        },
        "cpi_u_annual_average": CPI_U,
        "crosswalk": {
            "codes": {v: int((cw.vintage == v).sum()) for v in ("02", "10", "18")},
            "bls_departures_2000_code": BLS_DEPARTURES_2000,
            "occ1990_overrides_2010_code": OCC1990_OVERRIDES_2010,
            "k12_codes": {k: sorted(v) for k, v in K12_CODES.items()},
            "ambiguous_codes": int(cw.ambiguous.sum()),
            "share_of_2016_19_employed_in_ambiguous_codes": ambshare,
        },
        "sensitivity": sens,
    }
    payload = {"meta": meta, "samples": summaries, "changes": changes}
    EST.mkdir(parents=True, exist_ok=True)
    out_json = out_json or EST / "acs_occupations.json"
    out_json.write_text(json.dumps(_clean(payload), indent=1))
    LOG.info("wrote %s", out_json)
    cw.drop(columns=["group_dist"]).to_csv(EST / "acs_occupations_crosswalk.csv", index=False)
    out_md = out_md or EST / "acs_occupations.md"
    out_md.write_text(render_report(payload, moments, w, cw))
    LOG.info("wrote %s", out_md)


# =============================================================================
# Self-test
# =============================================================================

def selftest() -> None:
    """Unit checks of the helpers (no network)."""
    # weighted quantiles: Stata's rule averages when the cumulative weight hits the target exactly
    assert wquantile([1, 2, 3, 4], [1, 1, 1, 1], 0.5) == 2.5
    assert wquantile([1, 2, 3, 4], [1, 1, 1, 1], 0.5, "lower") == 2
    assert wquantile([1, 2, 3, 4], [1, 1, 1, 1], 0.9) == 4
    assert wquantile([4, 3, 2, 1], [1, 1, 1, 97], 0.1) == 1
    assert abs(ratio_9010(np.arange(1, 101), np.ones(100)) - (90.5 / 10.5)) < 1e-12
    assert abs(pooled_9010({"a": 2.0, "b": 4.0}, {"a": 1, "b": 3}) - 3.5) < 1e-12
    # schooling: grade 9 = 9, HS diploma or GED 12, some college 13, associate 14, bachelor 16, graduate 19 -> 17
    hg = schl_to_highgrade(pd.Series([1, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24]))
    assert hg.tolist() == [0, 8, 9, 12, 12, 12, 13, 13, 14, 16, 19, 19]
    assert hg.clip(9, 17).tolist() == [9, 9, 9, 12, 12, 12, 13, 13, 14, 16, 17, 17]
    # groups
    assert occ1990_to_group(22) == MARKET_GROUPS[0] and occ1990_to_group(23) == MARKET_GROUPS[1]
    assert occ1990_to_group(155) == OTHER_TEACH and occ1990_to_group(229) == MARKET_GROUPS[8]
    assert occ1990_to_group(999) is None and occ1990_to_group(905) is None
    assert len(MARKET_GROUPS) == 19 and len(NONTEACH) == 20 and len(ALL_GROUPS) == 21
    # codes and the BLS table
    assert _norm_code("10") == "0010" and _norm_code(2310.0) == "2310" and _norm_code("N.A.") is None and _norm_code(None) is None
    assert BLS_2000_TO_OCC1990[102] == 229 and BLS_2000_TO_OCC1990[231] == 156 and BLS_2000_TO_OCC1990[341] == 446
    assert len(BLS_2000_TO_OCC1990) == 509 and all(isinstance(v, int) for v in BLS_2000_TO_OCC1990.values())
    assert all(occ1990_to_group(j) for c, j in BLS_2000_TO_OCC1990.items() if c not in (0, 975, 980, 981, 982, 983, 990, 995, 999)
               ), "every civilian Census 2000 code maps into an HHJK group"
    assert "2340" in K12_CODES["10"] and "2350" in K12_CODES["18"] and "2540" not in K12_CODES["10"]
    # vintages
    assert occ_vintage("acs1", 2009, pd.Series([2009]))[0] == "02" and occ_vintage("acs1", 2017, pd.Series([2017]))[0] == "10"
    assert occ_vintage("acs1", 2018, pd.Series([2018]))[0] == "18"
    v5 = occ_vintage("acs5", 2013, pd.Series([2009, 2010, 2013]))
    assert v5.tolist() == ["02", "10", "10"]
    # build_sample on a toy frame: 15-29 hour workers split 0.5/0.5, unemployed rule, income floor
    cw = pd.DataFrame({"vintage": ["10", "10"], "occ": ["2310", "4020"], "occ1990": [156, 436],
                       "group": [K12, MARKET_GROUPS[12]], "group_alt": [K12, MARKET_GROUPS[12]]})
    base = dict(ds_year=2013, data_year=2013, vintage="10", gq=False, weeks=51.0, age=30, sex=2, schl=16, rac1p=1, mil=4,
                adjinc=1.0, pwgtp=10, semp=0, wagp=51000.0)
    rows = [dict(base, esr=1, wkhp=40, occ="2310"),               # full-time teacher, in the wage sample
            dict(base, esr=1, wkhp=20, occ="4020"),               # part-time cook: 0.5 / 0.5
            dict(base, esr=6, wkhp=0, occ=None),                  # not in the labor force: home
            dict(base, esr=3, wkhp=0, occ=None),                  # unemployed, no work last year: kept as home
            dict(base, esr=3, wkhp=40, occ="4020"),               # unemployed, 40 hours: dropped
            dict(base, esr=4, wkhp=40, occ=None),                 # armed forces: dropped
            dict(base, esr=1, wkhp=40, occ="2310", schl=11)]      # below grade 9: dropped
    d = build_sample(pd.DataFrame(rows), cw, rules="workbook")
    assert len(d) == 4 and abs(d.w_home.sum() - 2.5) < 1e-12 and abs(d.w_mkt.sum() - 1.5) < 1e-12
    assert d.wage_ok.sum() == 1 and abs(d.wage99[d.wage_ok].iloc[0] - 51000 / CPI_U[2013] * CPI_U[1999] / (40 * 51)) < 1e-9
    m = compute_moments(d.assign(gender="F"))
    assert m["F"]["count"][K12] == 1.0 and m["F"]["count"][HOME] == 2.5
    print("selftest ok")


# =============================================================================
# main
# =============================================================================

def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", nargs="?", default="all", choices=["all", "fetch", "build"])
    ap.add_argument("--check-bls-table", action="store_true",
                    help="re-read the BLS appendix PDF and check the embedded table (needs pdftotext)")
    ap.add_argument("--samples", nargs="+", default=list(SAMPLES), choices=list(SAMPLES),
                    help="samples to download (the build step always uses all three)")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if a.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    if a.selftest:
        selftest()
        return
    if a.check_bls_table:
        check_bls_table()
        return
    if a.step in ("all", "fetch"):
        cmd_fetch(a.samples)
    if a.step in ("all", "build"):
        cmd_build()


if __name__ == "__main__":
    main()
