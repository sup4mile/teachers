#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_data.py -- download and package the public data for the spatial teachers model.

Acquisition only.  Every table is downloaded, converted to a DataFrame, and
written beside this script.  No merging, harmonising, or moment construction --
that belongs downstream (see `notes/spatial_calibration_notes.md` for the moment
map these data feed).  The unit those notes assume -- and the reason this pulls
district-level tables plus a commuting-zone crosswalk -- is a school district
within a commuting zone.

Sources
    urban    Urban Institute Education Data API -- CCD directory/enrollment,
             F-33 finance, SAIPE child poverty, EDFacts assessments/grad rates
    seda     Stanford Education Data Archive -- district/county/metro/CZ mean
             scores and learning rates, covariates, school-district crosswalk
    edge     NCES EDGE -- CWIFT teacher wage index, district geocode/locale
    czones   Dorn's 1990 county -> 1990 commuting-zone crosswalk
    acs      Census ACS school-district tables (optional; needs a free API key)
    nlsy     NLSY79, NLSY79 Child/YA and NLSY97 test batteries for the ability
             block (optional; ~1.1 GB of zips, no key)

Output layout (under this script's directory unless --outdir says otherwise)
    raw/             one table per source-endpoint, as downloaded
    raw/_downloads/  the original csv/zip files
    metadata/        API endpoint and variable lists, SEDA docs, manifest.json

Usage
    python fetch_data.py                      # everything, full panel
    python fetch_data.py --dry-run            # print every URL, download nothing
    python fetch_data.py --sources urban --urban-topics ccd_finance
    python fetch_data.py --years 1990 2000 2010 2018    # benchmark years only
    python fetch_data.py --sources nlsy                 # the ability block only

Needs pandas, requests, pyarrow (+ openpyxl/xlrd for zipped spreadsheets).
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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:  # pragma: no cover
    print("This script needs `requests`:  pip install requests", file=sys.stderr)
    raise

LOG = logging.getLogger("fetch_data")
HERE = Path(__file__).resolve().parent   # data is stored beside the script
OUT_FMT = "parquet"                      # set from --format in main()


# =============================================================================
# Source addresses
# =============================================================================

URBAN_BASE = "https://educationdata.urban.org/api/v1"
URBAN_CSV_BASE = "https://educationdata.urban.org/csv"

# Stanford Digital Repository druids for SEDA releases.  The file list is read
# from the SDR at run time, so a new release only needs a new druid here.
SEDA_DRUIDS = {"6.0": "xh833nn4025", "5.0": "cs829jn7849",
               "4.1": "db586ns4974", "3.0": "mk782vn8293"}
SDR_PURL_JSON = "https://purl.stanford.edu/{druid}.json"
SDR_FILE = "https://stacks.stanford.edu/file/druid:{druid}/{filename}"

NCES_EDGE_DATA = "https://nces.ed.gov/programs/edge/data"
DORN_CZ_CROSSWALK = "https://www.ddorn.net/data/cw_cty_czone.zip"

# NCES has published CWIFT under several naming conventions; try each in turn.
CWIFT_YEARS = tuple(range(2013, 2024))
CWIFT_PATTERNS = ("EDGE_ACS_CWIFT{year}.zip", "EDGE_ACS_CWIFT{year}_LEA.zip",
                  "CWIFT{year}.zip")
# EDGE district geocode files are labelled by school year: 1314, 1415, ... 2324.
GEOCODE_URL = NCES_EDGE_DATA + "/EDGE_GEOCODE_PUBLICLEA_{sy}.zip"
GEOCODE_SCHOOL_YEARS = tuple(f"{y % 100:02d}{(y + 1) % 100:02d}"
                             for y in range(2013, 2024))

# ACS 5-year variables for school districts: the local income/property base.
ACS_SD_VARS = {
    "B19013_001E": "median_hh_income",
    "B19301_001E": "per_capita_income",
    "B25077_001E": "median_home_value",
    "B25103_001E": "median_property_tax",
    "B01003_001E": "population",
    "B15003_022E": "ba_holders",
    "B15003_001E": "pop_25plus",
    "B17001_002E": "pop_below_poverty",
}

# NLSY full-cohort files.  The access page links the current release; the
# pinned names below are the fallback if it cannot be read.
NLS_ACCESS_PAGE = "https://www.nlsinfo.org/accessing-data-cohorts"
NLS_COHORT_BASE = "https://nlsinfo.org/cohort-data"
NLSY_COHORT_FILES = {"nlsy79": "nlsy79_all_1979-2022",
                     "nlscya": "nlscya_all_1979-2020",
                     "nlsy97": "nlsy97_all_1997-2023"}
NLSY_TAGSET_EXT = {"nlsy79": ".NLSY79", "nlscya": ".CHILDYA", "nlsy97": ".NLSY97"}

# Variables for the ability measurement model (calibration item T2), as
# {block: regex over NLS question names}; a pattern must match the whole name.
# NLSY79 gives the parents' test battery, the Child/YA file the children's (with
# MPUBID = the mother's CASEID), and NLSY97 a second, CAT-administered battery
# with posterior variances.  Occupation, wage and schooling come along because
# test scores alone cannot separate common ability from comparative advantage.
NLSY_BLOCKS: dict[str, dict[str, str]] = {
    "nlsy79": {
        "id": r"CASEID|SAMPLE_ID|SAMPLE_RACE|SAMPLE_SEX|Q1-3_A~[MY]",
        "weights": r"SAMPWEIGHT|SAMPWEIGHT_ASVAB",
        # 1980 ASVAB, published 1981: section raw/scale/standard scores and
        # standard errors (ASVAB-1..43), AFQT (1980, 1989, 2006 norms), and the
        # IRT z-scores
        "asvab": r"ASVAB-\d+|AFQT-\d|ASVAB-[A-Z-]+-IRT-ZSCORE(_PERCENTILE)?",
        # item responses for the four AFQT sections, for reliability
        "asvab_items": r"ASVAB-(ARITHMETIC-REASONING|WORD-KNOWLEDGE"
                       r"|PARAGRAPH-COMPREHENSION|MATHEMATICS-KNOWLEDGE)-\d+",
        "age": r"AGEATINT",
        "family": r"HGC-(MOTHER|FATHER)|FAMOCC-(19|26)",
        "schooling": r"HGCREV\d\d|HGC_EVER",
        "work": r"CPSOCC(70|80)|OCCALL-EMP\.01|CPSHRP|HRP1",
    },
    "nlscya": {
        "id": r"CPUBID|MPUBID|CRACE|CSEX|CMOB|CYRB|BTHORDR",
        "weights": r"CSAMWT\d{4}(_REV)?",
        "assess_age": r"MSAGE\d{4}",          # age in months at assessment
        # PIAT math, reading recognition and comprehension, PPVT and digit span:
        # raw, percentile (P) and standard (Z) scores, 1986-2014
        "scores": r"(MATH|RECOG|COMP|PPVT)[PZ]?\d{4}|DIGITZ?\d{4}",
        "family": r"HGCREV\d{4}",             # the mother's schooling, by round
        "schooling": r"HGC\d{4}",             # the young adult's own schooling
    },
    "nlsy97": {
        "id": r"PUBID|KEY!SEX|KEY!BDATE_[MY]|KEY!RACE_ETHNICITY|CV_SAMPLE_TYPE",
        "weights": r"SAMPLING_WEIGHT_CC(_\d{4})?|SAMPLING_PANEL_WEIGHT(_\d{4})?",
        # 1997-98 CAT-ASVAB: ability estimates, posterior variances and items
        # completed per subtest; the AFQT-like percentile; the test date
        "asvab": r"ASVAB_[A-Z]{2}_(ABILITY_EST_(POS|NEG)|ITEMS_COMPLETE|POST_VARIANCE)"
                 r"|ASVAB_MATH_VERBAL_SCORE_PCT|ASVAB_TEST_DATE~[MY]|ASVAB!SPECIAL",
        "piat": r"CV_PIAT_(PERCENTILE|STANDARD)_(SCORE|UPD)|PIAT_RAW_SCORE_REVISED",
        "family": r"CV_HGC_(BIO|RES)_(MOM|DAD)",
        "schooling": r"CV_HGC_EVER_EDT|CV_HIGHEST_DEGREE_EVER_EDT",
        "work": r"YEMP_OCCODE-2002\.01|CV_HRLY_PAY\.01",
    },
}


# ---------------------------------------------------------------------------
# Urban Institute endpoints.  `years` is the documented coverage, intersected
# with --years at run time; `path` is the {level}/{source}/{topic} fragment and
# `subpath` is appended after the year (CCD enrollment needs a grade).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UrbanEndpoint:
    key: str
    path: str
    years: tuple[int, ...]
    subpath: str = ""
    filters: dict[str, str] = field(default_factory=dict)
    bulk_hint: str = ""      # file under /csv/ if a bulk CSV exists
    allow_bulk: bool = True  # False where the bulk file carries every breakout
    note: str = ""


def _yrs(*specs: Any) -> tuple[int, ...]:
    """_yrs(1986, (1994, 2018)) -> (1986, 1994, 1995, ..., 2018)"""
    out: list[int] = []
    for s in specs:
        if isinstance(s, tuple):
            out.extend(range(s[0], s[1] + 1))
        else:
            out.append(int(s))
    return tuple(sorted(set(out)))


URBAN_ENDPOINTS: dict[str, UrbanEndpoint] = {
    "ccd_directory": UrbanEndpoint(
        key="ccd_directory", path="school-districts/ccd/directory",
        years=_yrs((1986, 2022)), bulk_hint="ccd/districts_ccd_directory.csv",
        note="enrollment, teachers_total_fte, urban_centric_locale, county_code, cbsa"),
    "ccd_enrollment": UrbanEndpoint(
        key="ccd_enrollment", path="school-districts/ccd/enrollment",
        years=_yrs((1986, 2022)), subpath="grade-99",  # 99 = total across grades
        allow_bulk=False,  # the bulk file carries every grade x race x sex cell
        note="total enrollment by district-year"),
    "ccd_finance": UrbanEndpoint(
        key="ccd_finance", path="school-districts/ccd/finance",
        years=_yrs(1991, (1994, 2018)), bulk_hint="ccd/districts_ccd_finance.csv",
        note="F-33: revenue by source, expenditure, salaries"),
    "saipe": UrbanEndpoint(
        key="saipe", path="school-districts/saipe",
        years=_yrs(1995, 1997, (1999, 2021)), bulk_hint="saipe/districts_saipe.csv",
        note="child poverty 5-17; proxy for the local income base"),
    "edfacts_assessments": UrbanEndpoint(
        key="edfacts_assessments", path="school-districts/edfacts/assessments",
        years=_yrs((2009, 2018), 2020), subpath="grade-99",  # 99 = all grades
        note="proficiency rates; coarse but long-running quality measure"),
    "edfacts_grad_rates": UrbanEndpoint(
        key="edfacts_grad_rates", path="school-districts/edfacts/grad-rates",
        years=_yrs((2010, 2019)), note="adjusted cohort graduation rate"),
    # Opt-in only (school-level and heavy):
    "crdc_teachers": UrbanEndpoint(
        key="crdc_teachers", path="schools/crdc/teachers-staff",
        years=_yrs(2011, 2013, 2015, 2017),
        note="teacher experience/certification"),
}

DEFAULT_URBAN_TOPICS = ["ccd_directory", "ccd_enrollment", "ccd_finance", "saipe",
                        "edfacts_assessments", "edfacts_grad_rates"]
DEFAULT_SOURCES = ["urban", "seda", "edge", "czones"]
ALL_SOURCES = DEFAULT_SOURCES + ["acs", "nlsy"]


# =============================================================================
# Small utilities: HTTP, file I/O, and the fetch manifest
# =============================================================================

def make_session(retries: int = 5, backoff: float = 1.0, timeout: int = 120
                 ) -> requests.Session:
    """A requests session with retry/backoff for flaky federal servers."""
    s = requests.Session()
    retry = Retry(total=retries, connect=retries, read=retries, status=retries,
                  backoff_factor=backoff, status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset(["GET", "HEAD"]), raise_on_status=False)
    adapter = HTTPAdapter(max_retries=retry, pool_maxsize=16)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers.update({"User-Agent": "teachers-spatial-calibration/1.0 "
                                    "(academic research; contact: jpryan7@wisc.edu)"})
    s.request_timeout = timeout  # type: ignore[attr-defined]
    return s


def http_get(session: requests.Session, url: str, timeout: int | None = None,
             stream: bool = False) -> requests.Response:
    """GET with the session's default timeout unless one is given."""
    return session.get(url, stream=stream,
                       timeout=timeout or getattr(session, "request_timeout", 120))


def _ext(path: Path, new: str) -> Path:
    """Like Path.with_suffix, but safe for stems containing dots (e.g. `..._6.0`)."""
    return path.parent / (path.name + new)


def _stringify_mixed_object_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Stringify object columns that mix Python types across concatenated
    vintages (e.g. EDGE's STFIP: zero-padded str in one year, int in the
    next) -- pyarrow can infer one Arrow type per column but not a union,
    so left alone these raise "Expected bytes, got a 'int' object" deep in
    to_parquet. Every other column is returned untouched.
    """
    mixed = [c for c in df.columns if df[c].dtype == object
             and df[c].dropna().map(type).nunique() > 1]
    if not mixed:
        return df
    df = df.copy()
    for col in mixed:
        df[col] = df[col].map(lambda x: str(x) if pd.notna(x) else x)
    return df


def write_table(df: pd.DataFrame, path_no_ext: Path, fmt: str | None = None) -> None:
    """Package one downloaded table as parquet and/or CSV.

    Parquet is the default: the original CSV/zip is already kept verbatim under
    raw/_downloads, so a second plain-text copy of a multi-GB panel is waste.
    """
    fmt = fmt or OUT_FMT
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    if fmt in ("parquet", "both"):
        df = _stringify_mixed_object_columns(df)
        try:
            df.to_parquet(_ext(path_no_ext, ".parquet"), index=False)
        except Exception as exc:  # pyarrow missing, or mixed dtypes in a column
            LOG.warning("parquet write failed for %s (%s); writing CSV instead",
                        path_no_ext.name, exc)
            fmt = "csv"
    if fmt in ("csv", "both"):
        df.to_csv(_ext(path_no_ext, ".csv"), index=False)
    LOG.info("wrote %s  (%d rows x %d cols)", path_no_ext.name, len(df), df.shape[1])


def sniff_sep(head: bytes) -> str:
    """Guess a delimiter from a header line.

    Necessary because NCES ships tab-delimited files under a `.txt` extension
    (CWIFT is one), and reading those with a comma either collapses every field
    into one column or throws outright.
    """
    line = head.split(b"\n", 1)[0].decode("latin-1")
    return max(("\t", ",", "|", ";"), key=line.count)


def read_any(path: Path) -> pd.DataFrame:
    """Read csv / txt / tsv / dta / xlsx / parquet by extension."""
    suf = path.suffix.lower()
    if suf in (".csv", ".txt", ".tsv"):
        with open(path, "rb") as fh:
            sep = sniff_sep(fh.read(4096))
        return pd.read_csv(path, sep=sep, low_memory=False, encoding="latin-1")
    if suf == ".parquet":
        return pd.read_parquet(path)
    if suf == ".dta":
        return pd.read_stata(path, convert_categoricals=False)
    if suf in (".xls", ".xlsx"):
        return pd.read_excel(path)
    raise ValueError(f"don't know how to read {path.name}")


def download_zip_tables(session: requests.Session, url: str, dest_dir: Path
                        ) -> dict[str, pd.DataFrame]:
    """GET a zip, keep the archive, and return every tabular member as a frame.

    Two archive quirks are handled here.  macOS AppleDouble members
    (`__MACOSX/`, `._name`) are skipped -- several of these files were zipped on
    a Mac and the resource forks otherwise parse as one-row "tables" that shadow
    the real file.  And where several members share a stem, the richest format
    wins: an EDGE geocode archive ships the same table as .xlsx and as a
    header-less pipe-delimited .TXT, and only the spreadsheet names its columns.
    """
    resp = http_get(session, url, timeout=300)
    resp.raise_for_status()
    dest_dir.mkdir(parents=True, exist_ok=True)
    (dest_dir / Path(url).name).write_bytes(resp.content)

    rank = {".dta": 0, ".xlsx": 1, ".xls": 1, ".csv": 2, ".tsv": 2, ".txt": 3}
    out: dict[str, pd.DataFrame] = {}
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        members = [m for m in zf.namelist()
                   if not m.startswith("__MACOSX/")
                   and not Path(m).name.startswith("._")
                   and Path(m).suffix.lower() in rank]
        for member in sorted(members, key=lambda m: rank[Path(m).suffix.lower()]):
            stem = Path(member).stem
            if stem in out:
                continue
            data = zf.read(member)
            low = member.lower()
            try:
                if low.endswith((".xls", ".xlsx")):
                    df = pd.read_excel(io.BytesIO(data))
                elif low.endswith(".dta"):
                    df = pd.read_stata(io.BytesIO(data), convert_categoricals=False)
                else:
                    df = pd.read_csv(io.BytesIO(data), sep=sniff_sep(data[:4096]),
                                     low_memory=False, encoding="latin-1")
            except ImportError as exc:   # a missing reader silently degrades data
                LOG.warning("    cannot read %s (%s); see pyproject.toml deps", member, exc)
                continue
            except Exception as exc:
                LOG.debug("    skipping %s (%s)", member, exc)
                continue
            out[stem] = df
    return out


class Manifest:
    """Records what was attempted, what succeeded, and which URL it came from."""

    def __init__(self, path: Path, refreshing: Sequence[str] = ()):
        self.path = path
        self.refreshing = set(refreshing)
        self.entries: list[dict[str, Any]] = []

    def _carried_over(self) -> list[dict[str, Any]]:
        """Entries from earlier runs for sources this run did not touch.

        Re-pulling one source overwrites only that source's tables, so the
        manifest has to survive the same way -- a `--sources acs` rerun must not
        erase the provenance of the CCD panel sitting next to it.
        """
        if not self.path.exists():
            return []
        try:
            prior = json.loads(self.path.read_text())
        except Exception as exc:
            LOG.warning("manifest: could not read %s (%s); starting fresh",
                        self.path, exc)
            return []
        return [e for e in prior if e.get("source") not in self.refreshing]

    def add(self, source: str, name: str, url: str, status: str,
            rows: int | None = None, detail: str = "") -> None:
        self.entries.append({"source": source, "name": name, "url": url,
                             "status": status, "rows": rows, "detail": detail,
                             "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%S")})

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        kept = self._carried_over()
        self.path.write_text(json.dumps(kept + self.entries, indent=2))
        ok = sum(1 for e in self.entries if e["status"] == "ok")
        LOG.info("manifest: %d/%d fetches ok (%d carried over) -> %s",
                 ok, len(self.entries), len(kept), self.path)


# =============================================================================
# Urban Institute Education Data Portal
# =============================================================================

def urban_url(ep: UrbanEndpoint, year: int) -> str:
    """Build the API URL for one endpoint-year."""
    parts = [URBAN_BASE, ep.path, str(year)] + ([ep.subpath] if ep.subpath else [])
    url = "/".join(parts) + "/"
    if ep.filters:
        url += "?" + "&".join(f"{k}={v}" for k, v in ep.filters.items())
    return url


def urban_get_paged(session: requests.Session, url: str,
                    max_pages: int = 10_000) -> pd.DataFrame:
    """Follow the API's `next` links and concatenate `results`.

    The portal caps a page at 10,000 records, so a district-year pull is a
    handful of requests; a school-level pull is a few hundred.
    """
    frames: list[pd.DataFrame] = []
    nxt: str | None = url
    n_pages = 0
    while nxt and n_pages < max_pages:
        resp = http_get(session, nxt)
        if resp.status_code == 404:
            raise FileNotFoundError(f"404 {nxt}")
        resp.raise_for_status()
        payload = resp.json()
        results = payload.get("results", payload if isinstance(payload, list) else [])
        if results:
            frames.append(pd.DataFrame(results))
        nxt = payload.get("next") if isinstance(payload, dict) else None
        n_pages += 1
        if n_pages % 10 == 0:
            LOG.debug("  ... %d pages, %d rows so far", n_pages, sum(map(len, frames)))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def fetch_urban_metadata(session: requests.Session, meta_dir: Path,
                         manifest: Manifest, dry_run: bool) -> dict[str, pd.DataFrame]:
    """Pull the portal's own metadata: endpoint list, variable lists, bulk downloads.

    Worth doing every run -- it is the authoritative record of which variable
    names exist in which vintage, and it is what lets the bulk-CSV path resolve
    file names instead of hard-coding them.
    """
    out: dict[str, pd.DataFrame] = {}
    for name in ("api-endpoints", "api-downloads", "api-variables",
                 "api-endpoint-varlist"):
        url = f"{URBAN_BASE}/{name}/"
        if dry_run:
            print(f"[dry-run] GET {url}")
            manifest.add("urban", name, url, "dry-run")
            continue
        try:
            df = urban_get_paged(session, url)
            if not df.empty:
                write_table(df, meta_dir / f"urban_{name.replace('-', '_')}", fmt="csv")
                out[name] = df
            manifest.add("urban", name, url, "ok", len(df))
        except Exception as exc:
            LOG.warning("metadata %s failed: %s", name, exc)
            manifest.add("urban", name, url, "failed", detail=str(exc))
    return out


def fetch_urban_bulk(session: requests.Session, ep: UrbanEndpoint,
                     downloads: pd.DataFrame | None, raw_dir: Path,
                     manifest: Manifest, dry_run: bool) -> pd.DataFrame | None:
    """Try the bulk CSV for an endpoint (all years in one file) instead of the API.

    Much faster for a full annual panel.  Returns None if no bulk file resolves,
    in which case the caller falls back to year-by-year API calls.
    """
    candidates: list[str] = []
    # The api-downloads metadata table lists the bulk file for each endpoint.
    if downloads is not None and not downloads.empty:
        cols = {c.lower(): c for c in downloads.columns}
        ep_col = cols.get("endpoint_url") or cols.get("endpoint")
        dir_col = cols.get("file_dir") or cols.get("directory")
        name_col = cols.get("file_name") or cols.get("filename")
        if ep_col and dir_col and name_col:
            hit = downloads[downloads[ep_col].astype(str)
                            .str.contains(ep.path, regex=False, na=False)]
            candidates += [f"{URBAN_CSV_BASE}/{r[dir_col]}/{r[name_col]}"
                           for _, r in hit.iterrows()]
    if ep.bulk_hint:
        candidates.append(f"{URBAN_CSV_BASE}/{ep.bulk_hint}")

    seen: set[str] = set()
    for url in candidates:
        if url in seen:
            continue
        seen.add(url)
        if dry_run:
            print(f"[dry-run] GET (bulk) {url}")
            manifest.add("urban", f"{ep.key}:bulk", url, "dry-run")
            return None
        try:
            LOG.info("  bulk CSV: %s", url)
            resp = http_get(session, url, stream=True)
            if resp.status_code != 200:
                continue
            local = raw_dir / "_downloads" / Path(url).name
            local.parent.mkdir(parents=True, exist_ok=True)
            with open(local, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
            mb = local.stat().st_size / 1e6
            LOG.info("  downloaded %.0f MB", mb)
            if mb > 1500:
                LOG.warning("  %s is %.0f MB; if this exhausts memory, rerun with "
                            "--no-bulk to page the API instead", local.name, mb)
            df = pd.read_csv(local, low_memory=False, encoding="latin-1")
            manifest.add("urban", f"{ep.key}:bulk", url, "ok", len(df))
            return df
        except Exception as exc:
            LOG.debug("  bulk attempt failed (%s): %s", url, exc)
            manifest.add("urban", f"{ep.key}:bulk", url, "failed", detail=str(exc))
    return None


def fetch_urban(session: requests.Session, topics: Sequence[str],
                years: Sequence[int] | None, raw_dir: Path, meta_dir: Path,
                manifest: Manifest, dry_run: bool, use_bulk: bool) -> None:
    """Download each requested Urban endpoint and write one table per endpoint."""
    downloads = fetch_urban_metadata(session, meta_dir, manifest, dry_run).get("api-downloads")

    for topic in topics:
        ep = URBAN_ENDPOINTS.get(topic)
        if ep is None:
            LOG.warning("unknown urban topic %r (choices: %s)", topic,
                        ", ".join(URBAN_ENDPOINTS))
            continue
        wanted = [y for y in ep.years if years is None or y in set(years)]
        if not wanted:
            LOG.info("%s: no requested year overlaps its coverage %s-%s",
                     ep.key, min(ep.years), max(ep.years))
            continue
        LOG.info("== %s: %d years (%d-%d) ==", ep.key, len(wanted),
                 min(wanted), max(wanted))

        # One bulk file beats hundreds of API pages, but only covers all years.
        if use_bulk and ep.allow_bulk and years is None:
            bulk = fetch_urban_bulk(session, ep, downloads, raw_dir, manifest, dry_run)
            if bulk is not None and not bulk.empty:
                write_table(bulk, raw_dir / f"urban_{ep.key}")
                continue

        frames: list[pd.DataFrame] = []
        for year in wanted:
            url = urban_url(ep, year)
            if dry_run:
                print(f"[dry-run] GET {url}")
                manifest.add("urban", f"{ep.key}:{year}", url, "dry-run")
                continue
            try:
                df = urban_get_paged(session, url)
                if df.empty:
                    LOG.info("  %s %s: empty", ep.key, year)
                    manifest.add("urban", f"{ep.key}:{year}", url, "empty", 0)
                    continue
                if "year" not in df.columns:
                    df["year"] = year
                frames.append(df)
                LOG.info("  %s %s: %d rows", ep.key, year, len(df))
                manifest.add("urban", f"{ep.key}:{year}", url, "ok", len(df))
            except FileNotFoundError:
                LOG.info("  %s %s: not available (404)", ep.key, year)
                manifest.add("urban", f"{ep.key}:{year}", url, "404")
            except Exception as exc:
                LOG.warning("  %s %s failed: %s", ep.key, year, exc)
                manifest.add("urban", f"{ep.key}:{year}", url, "failed", detail=str(exc))
        if frames:
            write_table(pd.concat(frames, ignore_index=True), raw_dir / f"urban_{ep.key}")


# =============================================================================
# SEDA, via the Stanford Digital Repository
# =============================================================================

def sdr_list_files(session: requests.Session, druid: str) -> list[str]:
    """Ask the SDR for a druid's file list, so file names need not be hard-coded."""
    url = SDR_PURL_JSON.format(druid=druid)
    try:
        resp = http_get(session, url)
        resp.raise_for_status()
        meta = resp.json()
    except Exception as exc:
        LOG.warning("SDR listing failed for %s (%s); falling back to name patterns",
                    druid, exc)
        return []

    names: list[str] = []

    def walk(node: Any) -> None:  # the PURL JSON nests files a few levels deep
        if isinstance(node, dict):
            fn = node.get("filename") or node.get("label")
            if isinstance(fn, str) and re.search(r"\.(csv|dta|zip|pdf|xlsx?)$", fn, re.I):
                names.append(fn)
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)

    walk(meta)
    return sorted(set(names))


def seda_default_filenames(version: str) -> list[str]:
    """Fallback names if the SDR listing is unavailable.  SEDA's convention since
    v3.0 is seda_{level}_{form}_{scale}_{version}.csv."""
    names = [f"seda_{lv}_{form}_{scale}_{version}.csv"
             for lv in ("geodist", "county", "metro", "commzone")
             for form in ("annual", "long", "pool")
             for scale in ("cs", "gcs", "gys")]
    names += [f"seda_cov_{lv}_{form}_{version}.csv"
              for lv in ("geodist", "county", "metro")
              for form in ("long", "poolyr", "pool")]
    return names + [f"seda_crosswalk_{version}.csv"]


def fetch_seda(session: requests.Session, version: str, druid: str | None,
               raw_dir: Path, meta_dir: Path, manifest: Manifest,
               dry_run: bool, levels: Sequence[str]) -> None:
    """Download SEDA mean scores and learning rates at the requested geographies,
    plus the covariate files, the crosswalk, and the technical documentation."""
    druid = druid or SEDA_DRUIDS.get(version)
    if not druid:
        LOG.error("no SDR druid known for SEDA version %r; pass --seda-druid", version)
        return
    LOG.info("== SEDA %s (druid %s) ==", version, druid)

    listed = [] if dry_run else sdr_list_files(session, druid)
    if listed:
        (meta_dir / "seda_file_list.txt").parent.mkdir(parents=True, exist_ok=True)
        (meta_dir / "seda_file_list.txt").write_text("\n".join(listed))
        LOG.info("  SDR lists %d files", len(listed))
        csvs = [f for f in listed if f.lower().endswith(".csv")]
        wanted = sorted(set(
            [f for f in csvs if any(f"_{lv}_" in f.lower()
                                    or f.lower().startswith(f"seda_{lv}")
                                    for lv in levels)]
            + [f for f in csvs if "cov" in f.lower() or "crosswalk" in f.lower()]
        ))
        docs = [f for f in listed if f.lower().endswith(".pdf")]
    else:
        wanted = [f for f in seda_default_filenames(version)
                  if any(lv in f.lower() for lv in levels)
                  or "cov" in f.lower() or "crosswalk" in f.lower()]
        docs = [f"SEDA_documentation_{version}.pdf"]

    for fn in wanted + docs:
        url = SDR_FILE.format(druid=druid, filename=fn)
        if dry_run:
            print(f"[dry-run] GET {url}")
            manifest.add("seda", fn, url, "dry-run")
            continue
        local = raw_dir / "_downloads" / fn
        local.parent.mkdir(parents=True, exist_ok=True)
        try:
            resp = http_get(session, url, timeout=300, stream=True)
            if resp.status_code != 200:
                LOG.debug("  %s -> HTTP %s", fn, resp.status_code)
                manifest.add("seda", fn, url, f"http{resp.status_code}")
                continue
            with open(local, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
            if fn.lower().endswith(".pdf"):     # codebooks live with the metadata
                (meta_dir / fn).parent.mkdir(parents=True, exist_ok=True)
                (meta_dir / fn).write_bytes(local.read_bytes())
                manifest.add("seda", fn, url, "ok")
                continue
            df = read_any(local)
            stem = "seda_" + re.sub(r"[^a-z0-9]+", "_", Path(fn).stem.lower()).strip("_")
            write_table(df, raw_dir / stem.replace("seda_seda_", "seda_"))
            manifest.add("seda", fn, url, "ok", len(df))
            LOG.info("  %s: %d rows", fn, len(df))
        except Exception as exc:
            LOG.warning("  SEDA %s failed: %s", fn, exc)
            manifest.add("seda", fn, url, "failed", detail=str(exc))


# =============================================================================
# NCES EDGE and the commuting-zone crosswalk
# =============================================================================

def fetch_edge(session: requests.Session, raw_dir: Path, manifest: Manifest,
               dry_run: bool, cwift_years: Sequence[int],
               geocode_years: Sequence[str]) -> None:
    """CWIFT (the teacher-wage deflator) and the EDGE district geocode/locale file."""
    LOG.info("== NCES EDGE ==")

    # ---- CWIFT: one zip per year, name pattern varies by vintage ----------
    cwift_frames: list[pd.DataFrame] = []
    for year in cwift_years:
        got = False
        for pat in CWIFT_PATTERNS:
            url = f"{NCES_EDGE_DATA}/{pat.format(year=year)}"
            if dry_run:
                print(f"[dry-run] GET {url}")
                manifest.add("edge", f"cwift:{year}", url, "dry-run")
                got = True
                break
            try:
                tables = download_zip_tables(session, url, raw_dir / "_downloads")
            except Exception as exc:
                LOG.debug("  cwift %s (%s): %s", year, pat, exc)
                continue
            for name, df in tables.items():
                cwift_frames.append(df.assign(cwift_year=year, source_file=name))
            if tables:
                LOG.info("  CWIFT %s: %s", year, ", ".join(tables))
                manifest.add("edge", f"cwift:{year}", url, "ok",
                             sum(map(len, tables.values())))
                got = True
                break
        if not got:
            manifest.add("edge", f"cwift:{year}",
                         f"{NCES_EDGE_DATA}/{CWIFT_PATTERNS[0].format(year=year)}",
                         "not-found")
    if cwift_frames:
        write_table(pd.concat(cwift_frames, ignore_index=True), raw_dir / "edge_cwift")

    # ---- district geocode / locale: one zip per school year ---------------
    geo_frames: list[pd.DataFrame] = []
    for sy in geocode_years:
        url = GEOCODE_URL.format(sy=sy)
        if dry_run:
            print(f"[dry-run] GET {url}")
            manifest.add("edge", f"geocode:{sy}", url, "dry-run")
            continue
        try:
            tables = download_zip_tables(session, url, raw_dir / "_downloads")
        except Exception as exc:
            LOG.debug("  geocode %s: %s", sy, exc)
            manifest.add("edge", f"geocode:{sy}", url, "not-found", detail=str(exc))
            continue
        for df in tables.values():
            geo_frames.append(df.assign(school_year=sy))
        if tables:
            LOG.info("  EDGE geocode %s: %s", sy, ", ".join(tables))
            manifest.add("edge", f"geocode:{sy}", url, "ok",
                         sum(map(len, tables.values())))
    if geo_frames:
        geo = pd.concat(geo_frames, ignore_index=True)
        # SURVYEAR is only present in a couple of vintages, and formatted
        # differently in each ("2015-2016" vs. 2016); it duplicates
        # school_year, so drop it rather than let the mixed str/int/NaN
        # column break parquet's type inference.
        geo = geo.drop(columns=["SURVYEAR"], errors="ignore")
        write_table(geo, raw_dir / "edge_geocode_lea")


def fetch_czones(session: requests.Session, raw_dir: Path, manifest: Manifest,
                 dry_run: bool) -> None:
    """Dorn's 1990 county -> 1990 commuting-zone crosswalk (741 CZs)."""
    LOG.info("== commuting zone crosswalk ==")
    url = DORN_CZ_CROSSWALK
    if dry_run:
        print(f"[dry-run] GET {url}")
        manifest.add("czones", "cw_cty_czone", url, "dry-run")
        return
    try:
        tables = download_zip_tables(session, url, raw_dir / "_downloads")
        if not tables:
            raise RuntimeError("no tabular member in the archive")
        # The archive also ships a README; the crosswalk is the biggest table.
        name, df = max(tables.items(), key=lambda kv: len(kv[1]))
        write_table(df, raw_dir / "czone_crosswalk")
        manifest.add("czones", name, url, "ok", len(df))
        LOG.info("  %s: %d rows, columns %s", name, len(df), df.columns.tolist())
    except Exception as exc:
        LOG.warning("  CZ crosswalk failed: %s", exc)
        manifest.add("czones", "cw_cty_czone", url, "failed", detail=str(exc))


# =============================================================================
# Census ACS (optional; requires a free API key)
# =============================================================================

ACS_SD_LAYERS = ("school district (unified)", "school district (elementary)",
                 "school district (secondary)")
ACS_WORKERS = 8          # states fetched at once; the adapter pools 16
_ACS_UNKNOWN_VAR = re.compile(r"unknown variable '([^']+)'")


def acs_url(year: int, layer: str, state: str, getvars: Sequence[str],
            api_key: str | None = None) -> str:
    """One ACS 5-year district query.  `api_key` is left off for the manifest."""
    url = (f"https://api.census.gov/data/{year}/acs/acs5"
           f"?get={','.join(['NAME'] + list(getvars))}"
           f"&for={layer.replace(' ', '%20')}:*&in=state:{state}")
    return url + (f"&key={api_key}" if api_key else "")


def acs_states(session: requests.Session, year: int,
               api_key: str | None) -> list[str]:
    """State FIPS the vintage publishes.  The district layers refuse
    `in=state:*` before 2022, so every year is fetched state by state."""
    url = f"https://api.census.gov/data/{year}/acs/acs5?get=NAME&for=state:*"
    resp = http_get(session, url + (f"&key={api_key}" if api_key else ""),
                    timeout=60)
    resp.raise_for_status()
    return [row[-1] for row in resp.json()[1:]]


def acs_supported_vars(session: requests.Session, year: int, state: str,
                       api_key: str | None) -> list[str]:
    """Drop the variables this vintage does not publish.

    B25103 (median property tax) starts in 2010 and B15003 (detailed
    attainment) in 2012, so the full `get=` is a 400 for the earliest years.
    The API names one unknown variable per error, so ask until it stops.
    """
    getvars = list(ACS_SD_VARS)
    for _ in range(len(ACS_SD_VARS)):
        resp = http_get(session, acs_url(year, ACS_SD_LAYERS[0], state,
                                         getvars, api_key), timeout=180)
        if resp.status_code in (200, 204):
            break
        missing = _ACS_UNKNOWN_VAR.search(resp.text)
        if not missing:
            break             # some other failure; the caller records it per state
        LOG.info("  ACS %s: %s not published this vintage", year,
                 missing.group(1))
        getvars = [v for v in getvars if v != missing.group(1)]
    return getvars


def fetch_acs_school_districts(session: requests.Session, years: Sequence[int],
                               api_key: str | None, raw_dir: Path,
                               manifest: Manifest, dry_run: bool) -> None:
    """ACS 5-year income, home value and property tax by school district.

    Key from https://api.census.gov/data/key_signup.html; the three district
    layers (unified / elementary / secondary) are stacked into one table.
    """
    LOG.info("== Census ACS (school districts) ==")
    frames: list[pd.DataFrame] = []
    for year in years:
        if dry_run:
            for layer in ACS_SD_LAYERS:
                url = acs_url(year, layer, "*", ACS_SD_VARS)
                print(f"[dry-run] GET {url}")
                manifest.add("acs", f"{year}:{layer}", url, "dry-run")
            continue
        try:
            states = acs_states(session, year, api_key)
        except Exception as exc:
            LOG.warning("  ACS %s: state list failed: %s", year, exc)
            manifest.add("acs", f"{year}:states",
                         f"https://api.census.gov/data/{year}/acs/acs5?get=NAME"
                         "&for=state:*", "failed", detail=str(exc))
            continue
        getvars = acs_supported_vars(session, year, states[0], api_key)
        for layer in ACS_SD_LAYERS:
            public = acs_url(year, layer, "*", getvars)
            # One state per request (the wildcard is refused before 2022) and
            # ~6 s each, so 14 years x 3 layers x 52 states is hours serially.
            def one_state(state: str) -> tuple[str, list[list[str]] | None, str]:
                try:
                    resp = http_get(session, acs_url(year, layer, state, getvars,
                                                     api_key), timeout=180)
                    if resp.status_code == 204:   # no district of this layer here
                        return state, [], ""
                    if resp.status_code != 200:
                        return state, None, f"{state}:http{resp.status_code}"
                    return state, resp.json(), ""
                except Exception as exc:
                    return state, None, f"{state}:{exc}"

            with ThreadPoolExecutor(max_workers=ACS_WORKERS) as pool:
                results = list(pool.map(one_state, states))
            rows_out: list[list[str]] = []
            header: list[str] | None = None
            bad: list[str] = []
            for _state, rows, err in results:
                if err:
                    bad.append(err)
                elif rows:
                    header = header or rows[0]
                    rows_out.extend(rows[1:])
            if header is None:
                LOG.warning("  ACS %s %s: nothing returned (%s)", year, layer,
                            "; ".join(bad[:3]) or "no districts")
                manifest.add("acs", f"{year}:{layer}", public, "failed",
                             detail="; ".join(bad[:10]))
                continue
            df = pd.DataFrame(rows_out, columns=header).rename(columns=ACS_SD_VARS)
            frames.append(df.assign(year=year, sd_layer=layer))
            LOG.info("  ACS %s %s: %d districts%s", year, layer, len(df),
                     f" ({len(bad)} states failed)" if bad else "")
            manifest.add("acs", f"{year}:{layer}", public,
                         "ok" if not bad else "partial", len(df),
                         detail="; ".join(bad[:10]))
    if frames:
        write_table(pd.concat(frames, ignore_index=True), raw_dir / "acs_school_districts")


# =============================================================================
# NLSY test batteries (the ability block; optional, ~1.1 GB of zips)
# =============================================================================

def nlsy_cohort_urls(session: requests.Session, cohorts: Sequence[str],
                     dry_run: bool) -> dict[str, str]:
    """Current full-cohort zip for each cohort, read off the NLS access page.

    The file name carries the last survey round (`nlsy79_all_1979-2022.zip`),
    so it changes with every release; the pinned names are only a fallback.
    """
    urls = {c: f"{NLS_COHORT_BASE}/{NLSY_COHORT_FILES[c]}.zip" for c in cohorts}
    if dry_run:
        return urls
    try:
        resp = http_get(session, NLS_ACCESS_PAGE, timeout=60)
        resp.raise_for_status()
        for c in cohorts:
            hits = sorted(set(re.findall(rf"{c}_all_\d{{4}}-\d{{4}}\.zip", resp.text)))
            if hits:
                urls[c] = f"{NLS_COHORT_BASE}/{hits[-1]}"
    except Exception as exc:
        LOG.warning("NLS access page unreadable (%s); using pinned file names", exc)
    return urls


def download_file(session: requests.Session, url: str, dest: Path) -> Path:
    """Stream a large file to disk, reusing a local copy of the same size.

    The NLSY zips run to ~450 MB each, so they are never held in memory, and a
    rerun that only changes the variable selection does not fetch them again.
    """
    head = session.head(url, allow_redirects=True, timeout=60)
    size = int(head.headers.get("content-length", -1))
    if dest.exists() and dest.stat().st_size == size:
        LOG.info("  reusing %s (%.0f MB)", dest.name, size / 2**20)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.parent / (dest.name + ".part")
    with http_get(session, url, timeout=600, stream=True) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 22):
                fh.write(chunk)
    tmp.replace(dest)
    LOG.info("  downloaded %s (%.0f MB)", dest.name, dest.stat().st_size / 2**20)
    return dest


def read_nlsy_sdf(text: str) -> pd.DataFrame:
    """Parse an NLS `.sdf` variable index: one fixed-width line per variable.

    Column offsets are read off the header line rather than hard-coded.  The
    reference number `R06150.00` is the CSV column `R0615000`.
    """
    lines = text.splitlines()
    hdr_at = next(i for i, ln in enumerate(lines) if ln.startswith("Number"))
    hdr = lines[hdr_at]
    c_year, c_title, c_qname = (hdr.index(k) for k in
                                ("Year", "Variable Description", "Question Name"))
    rows = []
    for ln in lines[hdr_at + 2:]:          # skip the dashed rule under the header
        if not ln.strip():
            continue
        refnum = ln[:c_year].strip()
        rows.append((refnum.replace(".", ""), refnum, ln[c_year:c_title].strip(),
                     ln[c_qname:].strip(), ln[c_title:c_qname].strip()))
    return pd.DataFrame(rows, columns=["rnum", "refnum", "year", "qname", "title"])


def select_nlsy_vars(index: pd.DataFrame, blocks: dict[str, str]) -> pd.DataFrame:
    """Rows of the variable index whose question name matches a block's regex.

    Question names are stable across releases (MATH1996, ASVAB-3, KEY!SEX);
    reference numbers are too, but naming the variables this way picks up new
    survey rounds without editing a list of R-numbers.
    """
    picked = []
    for block, pattern in blocks.items():
        hit = index[index["qname"].str.fullmatch(pattern)]
        if hit.empty:
            LOG.warning("  block %r matched no variables (%s)", block, pattern)
        picked.append(hit.assign(block=block))
    out = pd.concat(picked, ignore_index=True)
    return out.drop_duplicates("rnum").reset_index(drop=True)


def fetch_nlsy(session: requests.Session, cohorts: Sequence[str], raw_dir: Path,
               meta_dir: Path, manifest: Manifest, dry_run: bool) -> None:
    """Test scores, ages, family links and outcomes for the ability measurement
    model, cut from the public full-cohort NLSY files.

    NLS Investigator needs a login, but the NLS also publishes each cohort
    whole: one zip with a CSV of every variable (columns are reference numbers),
    a codebook, and a variable index.  The zip is kept verbatim; the variables
    in NLSY_BLOCKS are cut from the CSV into `raw/<cohort>_ability`, and each
    cohort's selection is written to `metadata/` twice -- as a dictionary
    (reference number, survey year, question name, title, block) and as an
    Investigator tagset, which can be uploaded to browse the codebook entries.
    Values are left as released: negative codes are NLS missing-value flags.
    """
    LOG.info("== NLSY (ability block) ==")
    urls = nlsy_cohort_urls(session, cohorts, dry_run)
    for cohort in cohorts:
        url = urls[cohort]
        if dry_run:
            print(f"[dry-run] GET {url}")
            manifest.add("nlsy", cohort, url, "dry-run")
            continue
        try:
            local = download_file(session, url, raw_dir / "_downloads" / Path(url).name)
            with zipfile.ZipFile(local) as zf:
                members = {Path(m).suffix.lower(): m for m in zf.namelist()}
                index = read_nlsy_sdf(zf.read(members[".sdf"]).decode("latin-1"))
                sel = select_nlsy_vars(index, NLSY_BLOCKS[cohort])
                with zf.open(members[".csv"]) as fh:
                    df = pd.read_csv(fh, usecols=sel["rnum"].tolist())
            df = df[sel["rnum"].tolist()]           # index order, not CSV order
            write_table(df, raw_dir / f"{cohort}_ability")
            sel.to_csv(meta_dir / f"{cohort}_ability_vars.csv", index=False)
            (meta_dir / f"{cohort}_ability{NLSY_TAGSET_EXT[cohort]}").write_text(
                "\n".join(sel["rnum"]) + "\n")
            LOG.info("  %s: %d respondents x %d variables (%s)", cohort, len(df),
                     df.shape[1], ", ".join(f"{b} {n}" for b, n in
                                             sel["block"].value_counts(sort=False).items()))
            manifest.add("nlsy", cohort, url, "ok", len(df),
                         detail=f"{df.shape[1]} variables of {len(index)}")
        except Exception as exc:
            LOG.warning("  NLSY %s failed: %s", cohort, exc)
            manifest.add("nlsy", cohort, url, "failed", detail=str(exc))


# =============================================================================
# CLI
# =============================================================================

def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=str(HERE),
                   help="output root (default: this script's directory)")
    p.add_argument("--sources", nargs="+", default=DEFAULT_SOURCES, choices=ALL_SOURCES,
                   help=f"which sources to fetch (default: {' '.join(DEFAULT_SOURCES)})")
    p.add_argument("--urban-topics", nargs="+", default=DEFAULT_URBAN_TOPICS,
                   choices=sorted(URBAN_ENDPOINTS),
                   help="which Urban Institute endpoints to pull")
    p.add_argument("--years", nargs="+", type=int, default=None,
                   help="restrict to these years (default: the full documented panel)")
    p.add_argument("--seda-version", default="6.0",
                   choices=sorted(SEDA_DRUIDS) + ["custom"], help="SEDA release")
    p.add_argument("--seda-druid", default=None,
                   help="override the Stanford Digital Repository druid for SEDA")
    p.add_argument("--seda-levels", nargs="+",
                   default=["geodist", "commzone", "county", "metro"],
                   help="SEDA geographies to download")
    p.add_argument("--census-key", default=os.environ.get("CENSUS_API_KEY"),
                   help="Census API key for --sources acs (or set CENSUS_API_KEY)")
    p.add_argument("--nlsy-cohorts", nargs="+", default=list(NLSY_COHORT_FILES),
                   choices=list(NLSY_COHORT_FILES),
                   help="NLSY cohorts for --sources nlsy (default: all three)")
    p.add_argument("--no-bulk", action="store_true",
                   help="always page the API instead of trying the bulk CSVs")
    p.add_argument("--format", choices=["parquet", "csv", "both"], default="parquet",
                   help="format for the packaged tables (default: parquet)")
    p.add_argument("--dry-run", action="store_true",
                   help="print every URL that would be requested and exit")
    p.add_argument("--timeout", type=int, default=180, help="per-request timeout (s)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    global OUT_FMT
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    OUT_FMT = args.format

    out = Path(args.outdir).expanduser().resolve()
    raw_dir, meta_dir = out / "raw", out / "metadata"
    for d in (raw_dir, raw_dir / "_downloads", meta_dir):
        d.mkdir(parents=True, exist_ok=True)
    LOG.info("output root: %s", out)

    manifest = Manifest(meta_dir / "manifest.json", refreshing=args.sources)
    session = make_session(timeout=args.timeout)

    if "urban" in args.sources:
        fetch_urban(session, args.urban_topics, args.years, raw_dir, meta_dir,
                    manifest, args.dry_run, use_bulk=not args.no_bulk)
    if "seda" in args.sources:
        fetch_seda(session, args.seda_version, args.seda_druid, raw_dir, meta_dir,
                   manifest, args.dry_run, args.seda_levels)
    if "edge" in args.sources:
        fetch_edge(session, raw_dir, manifest, args.dry_run,
                   CWIFT_YEARS, GEOCODE_SCHOOL_YEARS)
    if "czones" in args.sources:
        fetch_czones(session, raw_dir, manifest, args.dry_run)
    if "acs" in args.sources:
        fetch_acs_school_districts(session, args.years or list(range(2009, 2023)),
                                   args.census_key, raw_dir, manifest, args.dry_run)
    if "nlsy" in args.sources:
        fetch_nlsy(session, args.nlsy_cohorts, raw_dir, meta_dir, manifest, args.dry_run)
    if not args.dry_run:          # a dry run must not rewrite the real record
        manifest.save()

    if args.dry_run:
        LOG.info("dry run complete -- nothing downloaded")
    else:
        LOG.info("done. data written to %s", raw_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
