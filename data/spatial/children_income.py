#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
children_income.py -- a within-commuting-zone counterpart of the NLSY79 parental-income sorting moment V6
(DESIGN_v2.md task D5): the share of the CZ's school-age children in each national family-income tercile
who live in the city, and top minus bottom, from the NCES EDGE ACS-ED school-district tabulations.

V6 is the NLSY79 central-city share of MSA children by tercile of parental income (45.8 / 36.0 / 16.3 %,
top minus bottom -0.295, SE 0.037): a 1979 statistic pooled across SMSAs, and silent on the 20-40% of
SMSA residents with unknown central-city status.  Here the same object is built inside each commuting
zone of the calibration, on ACS 2014-18, so that it is free of between-MSA composition.

Source (documented in detail in children_income.md, section 1)
    NCES EDGE, ACS-ED Children's Tabulation, ACS 5-year 2014-2018 ("2014-18" vintage, school year 2017-18),
    Profile Tables, school-district level (all 11,910 district records: 10,037 unified, 1,437 elementary,
    436 secondary), downloaded with the EDGE Table Viewer (https://nces.ed.gov/programs/edge/TableViewer/
    acsProfile/2018, "All Districts", Download Text File).  Population groups ("iterations"):
        102  Relevant Children                          primary: every child who lives in the district and
                                                        whose assigned grade the district serves, whatever
                                                        the school type (NLSY-comparable: all children)
        104  Relevant Children - Enrolled (public and private)
        105  Relevant Children - Enrolled Public        children enrolled in public school
    Tables (variables named as in the layout file, ACS-ED_2014-2018_RecordLayouts.xlsx):
        CDP03.2  INCOME AND BENEFITS: Families (CDP03_34est) by FAMILY INCOME in ten brackets (CDP03_35est-
                 CDP03_44est: < $10k, 10-15, 15-25, 25-35, 35-50, 50-75, 75-100, 100-150, 150-200, $200k+),
                 median (CDP03_45est) and mean (CDP03_46est) family income
        CDP05.1  SEX AND AGE: total population of the iteration = the number of children (CDP05_1est)
        CDP02.1  HOUSEHOLDS BY TYPE: average number of children in families with children (CDP02_9est)
        CDP03.4  poverty: related children 5-17 below the poverty level, percent (CDP03_67pct)
    The unit of the income table is the FAMILY (household) in which the iteration's children live, not the
    child.  The ACS-ED has no child-by-family-income cross-tabulation.  Children are obtained by scaling
    each district's brackets by its children per family (CDP05_1est / CDP03_34est): exact for the number
    of children in a district, an assumption (equal family size across income brackets within a district)
    for the distribution.  A family with children relevant to two overlapping districts (an elementary and
    a secondary district) appears in both records, children are counted once; the family-weighted variant
    is carried as a sensitivity.  Census API B17024 (age by income-to-poverty ratio, children 6-17) is a
    cross-check in a different income concept and on residence rather than relevance.

Definitions (children_income.md section 1)
    The calibration's panel: `district_settings.build_audited_panel` (New York City one district, carrying
    the Chancellor's id 3620580 under which the ACS-ED files it), locale partition (city = NCES locale 11-13,
    suburb = 21-23; districts of >= 500 pupils; each location >= 10% of the CZ's partitioned pupils: 188 CZs),
    CZ weights = partitioned pupils, whole-CZ bootstrap B = 1999, seed 20260929 with the multiplicity matrix of
    `district_settings.py` (so the bootstrap correlations with the four district targets are valid).  The
    four targets (0.3764, 0.0086, -0.0177, 0.1122) are reproduced and the script stops if they are not.
    National terciles: the cut points of the national distribution of children by family income, from the
    sum of all ACS-ED district records, interpolated linearly inside the bracket that holds the cut (the
    bracket is split in the same proportions in every district).  Within a CZ and location, the children in
    tercile t are the bracket counts times the share of each bracket in t; the city share of tercile t is
    T_city,t / (T_city,t + T_suburb,t); the within-CZ analogue of V6 is the share of tercile 3 minus the share
    of tercile 1, averaged over the CZs with the CZ weights.  Also: the composition of each location by
    bracket and tercile; the log suburb-minus-city gap in the median (pooled bracket counts, linear
    interpolation), in the mean and in the average district median family income of children; the
    suburb-minus-city gap in the child poverty rate; the variants listed in children_income.md section 4.

Hold-out.  The pre-registered split of the locale CZs (`district_settings.holdout_split`, seed 20261002) is
reproduced and checked against `district_settings.json`.  Half A is reported; half B goes ONLY to
estimates/settings/sealed/children_income_holdout_B.json and is never printed, logged or returned.

Usage (from data/spatial; on a compute node for `build`, see the project notes)
    uv run python children_income.py fetch                 # ACS-ED exports + Census B17024 (a few minutes)
    uv run python children_income.py fetch --only acsed    # or: --only b17024
    uv run python children_income.py build                 # B = 1999, seed 20260929
    uv run python children_income.py build --B 199 --outdir <dir>     # quick check
    uv run python children_income.py selftest              # unit checks of the helpers (no data, no network)
Writes estimates/settings/children_income.{md,json}, children_income_targets.toml,
children_income_draws.csv, children_income_draws_holdout_A.csv and sealed/children_income_holdout_B.json.
Raw downloads are cached under raw/children_income/ (gitignored).  The Census API key (CENSUS_API_KEY, from
the environment or data/spatial/.env) is used only for the B17024 cross-check and is never printed or
written; error strings are scrubbed of it.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import math
import re
import sys
import time
import zipfile
from dataclasses import dataclass
from html import unescape
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

LOG = logging.getLogger("children_income")
de = da = db = ds = None          # the calibration stack, imported by `_stack()` (not needed by `fetch`)


def _stack() -> None:
    global de, da, db, ds
    if ds is None:
        import data_estimate as _de
        import district_audit as _da
        import district_bootstrap as _db
        import district_settings as _ds
        de, da, db, ds = _de, _da, _db, _ds

RAW = HERE / "raw" / "children_income"
YEAR = 2018                                   # ACS 5-year vintage 2014-18, dollars of 2018
SEED = 20260929                               # = district_bootstrap.SEED
B_DEFAULT = 1999                              # = district_bootstrap.B_DEFAULT
UA = "teachers-spatial-calibration/1.0 (academic research; jpryan7@wisc.edu)"

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------

TV = "https://nces.ed.gov/programs/edge/TableViewer/"
TV_PAGE = TV + "acsProfile/2018"
TV_DATASET = "acsprofile"
ITERS = {
    "102": "Relevant Children",
    "104": "Relevant Children - Enrolled (public and private)",
    "105": "Relevant Children - Enrolled Public",
}
PRIMARY_IT = "102"
ACSED_TABLES = ("CDP03.2", "CDP05.1", "CDP02.1", "CDP03.4")
GEO_ALL = "USSchoolDistrictAll"
LAYER_BY_PREFIX = {"95000US": "elementary", "96000US": "secondary", "97000US": "unified"}
NYC_LEAID = "3620580"
# ACS-ED / Census files the Detroit school district under its pre-2016 id (Detroit City School District,
# 2612000); the CCD frame carries it as Detroit Public Schools Community District (2601103), 49,931 pupils
# in the city location of CZ 11600.  Source id -> panel id, applied to both ACS-ED and B17024.
DETROIT = ("2612000", "2601103")

# family-income brackets of CDP03.2 (dollars of 2018), CDP03_35est ... CDP03_44est
INC_EDGES = np.array([0, 10_000, 15_000, 25_000, 35_000, 50_000, 75_000, 100_000, 150_000, 200_000, np.inf])
INC_LABELS = ["< $10k", "$10-15k", "$15-25k", "$25-35k", "$35-50k", "$50-75k", "$75-100k", "$100-150k",
              "$150-200k", "$200k+"]
NB = len(INC_LABELS)
INC_COLS = [f"CDP03_{35 + k}est" for k in range(NB)]

# B17024 (age by ratio of income to poverty level), children 6 to 11 and 12 to 17, 12 bins each
IPR_EDGES = np.array([0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 1.85, 2.0, 3.0, 4.0, 5.0, np.inf])
IPR_LABELS = ["<0.50", "0.50-0.74", "0.75-0.99", "1.00-1.24", "1.25-1.49", "1.50-1.74", "1.75-1.84",
              "1.85-1.99", "2.00-2.99", "3.00-3.99", "4.00-4.99", "5.00+"]
B17_6_11 = [f"B17024_{j:03d}E" for j in range(16, 28)]
B17_12_17 = [f"B17024_{j:03d}E" for j in range(29, 41)]
B17_VARS = ["B17024_001E", "B17024_015E", "B17024_028E"] + B17_6_11 + B17_12_17

# NLSY79 V6 (calibration_log.md section 5; nlsy_extensions.md T3; settings/nlsy_composition.md)
NLSY_V6 = dict(city_share=[0.458, 0.360, 0.163], gradient=-0.295, gradient_se=0.037)


def acsed_zip(table: str, it: str) -> Path:
    return RAW / f"acsed_{table}_{it}_districts.zip"


def us_json(table: str, it: str) -> Path:
    return RAW / f"acsed_US_{table}_{it}.json"


B17024_PATH = RAW / f"acs_b17024_school_districts_{YEAR}.parquet"


# --------------------------------------------------------------------------
# Small utilities
# --------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fmt(v: Any, d: int = 3) -> str:
    return "" if v is None or (isinstance(v, (float, np.floating)) and not np.isfinite(v)) else f"{v:,.{d}f}"


def fse(v: Any, se: Any, d: int = 3) -> str:
    if v is None or se is None or not np.isfinite(v) or not np.isfinite(se):
        return fmt(v, d)
    return f"{v:,.{d}f} ({se:,.{d}f})"


def pc(v: Any, d: int = 1) -> str:
    return "" if v is None or not np.isfinite(v) else f"{100 * v:.{d}f}"


def pcs(v: Any, se: Any, d: int = 1) -> str:
    if v is None or se is None or not np.isfinite(v) or not np.isfinite(se):
        return pc(v, d)
    return f"{100 * v:.{d}f} ({100 * se:.{d}f})"


# --------------------------------------------------------------------------
# fetch: ACS-ED exports (NCES EDGE Table Viewer) and the Census B17024 cross-check
# --------------------------------------------------------------------------

def _session():
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    s = requests.Session()
    retry = Retry(total=5, connect=5, read=5, status=5, backoff_factor=1.5,
                  status_forcelist=(429, 500, 502, 503, 504), allowed_methods=None,
                  raise_on_status=False)
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=4))
    s.headers.update({"User-Agent": UA})
    return s


def _tv_post(sess, endpoint: str, body: str) -> Any:
    """POST to the Table Viewer's JSON API (the page's own AJAX calls; the body is the page's
    non-strict JSON string)."""
    r = sess.post(TV + endpoint, data=body.encode("utf-8"),
                  headers={"Content-Type": "application/json; charset=utf-8"}, timeout=600)
    r.raise_for_status()
    return r.json()


def export_table(sess, table: str, it: str, geo: str, dest: Path, tries: int = 4) -> None:
    """One table, one population group, one geography: the Table Viewer builds a zip holding the data
    file (`<table>_<it>_<geo>_*.txt`, pipe-delimited), the ReadMe and the record layout workbook."""
    attrs = f"{table}---{it}-{geo}-{geo}"
    body = "{iYear : %d, sDataset : '%s', sTableAttrSets : '%s'}" % (YEAR, TV_DATASET, attrs)
    last = ""
    for k in range(tries):
        try:
            j = _tv_post(sess, "TableExportToCSV", body)
            fn = j.get("filename") or ""
            if not fn.endswith(".zip"):
                last = f"export refused: {fn[:80]}"
                time.sleep(5 * (k + 1))
                continue
            url = "https://nces.ed.gov/" + (j["path"] or "").replace("\\", "/").strip("/") + "/" + fn
            r = sess.get(url, timeout=600)
            r.raise_for_status()
            if r.content[:2] != b"PK":
                last = "not a zip"
                continue
            dest.write_bytes(r.content)
            return
        except Exception as exc:                                         # noqa: BLE001
            last = f"{type(exc).__name__}: {exc}"[:160]
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"could not export {attrs}: {last}")


def fetch_acsed(refresh: bool = False) -> dict[str, Any]:
    RAW.mkdir(parents=True, exist_ok=True)
    sess = _session()
    manifest: dict[str, Any] = dict(source=TV_PAGE, dataset=TV_DATASET, year_label="ACS 2014-2018",
                                    fetched=dt.datetime.now().isoformat(timespec="seconds"), files={})
    for it in ITERS:
        for tb in ACSED_TABLES:
            dest = acsed_zip(tb, it)
            if dest.exists() and not refresh:
                LOG.info("cached   %s", dest.name)
            else:
                export_table(sess, tb, it, GEO_ALL, dest)
                LOG.info("fetched  %s (%s bytes)", dest.name, f"{dest.stat().st_size:,}")
                time.sleep(1.0)
            with zipfile.ZipFile(dest) as z:
                data = [n for n in z.namelist() if n.endswith(".txt")]
                nrows = sum(1 for _ in z.open(data[0])) - 1
                manifest["files"][dest.name] = dict(bytes=dest.stat().st_size, sha256=_sha256(dest),
                                                    data_file=data[0], rows=nrows,
                                                    attrs=f"{tb}---{it}-{GEO_ALL}-{GEO_ALL}")
        # the national row of the income table, for the cross-check (the "US" export is empty; the
        # preview call returns the national table)
        dj = us_json("CDP03.2", it)
        if refresh or not dj.exists():
            j = _tv_post(sess, "TableGet",
                         "{iYear : %d, sDataset : '%s', sTableId : 'CDP03.2', sUniverse : null, sGradeAge : null, "
                         "sEnrollmentType : '%s', sGeoType : 'US', sGeoId : 'US'}" % (YEAR, TV_DATASET, it))
            dj.write_text(json.dumps(j))
            time.sleep(0.5)
        manifest["files"][dj.name] = dict(bytes=dj.stat().st_size, sha256=_sha256(dj))
    # documentation shipped inside every export
    z0 = zipfile.ZipFile(acsed_zip("CDP03.2", PRIMARY_IT))
    for n in z0.namelist():
        if n.startswith("ReadMe") or n.endswith("RecordLayouts.xlsx"):
            tgt = RAW / f"acsed_{n}"
            if not tgt.exists() or refresh:
                tgt.write_bytes(z0.read(n))
            manifest["files"][tgt.name] = dict(bytes=tgt.stat().st_size, sha256=_sha256(tgt))
    (RAW / "manifest_acsed.json").write_text(json.dumps(manifest, indent=1))
    return manifest


def fetch_b17024(workers: int = 8, refresh: bool = False) -> Path:
    """Census API B17024 (age by ratio of income to poverty level) for the three school-district layers,
    one state per request as `acs_earnings.py` does; raw/children_income/acs_b17024_school_districts_<year>."""
    if B17024_PATH.exists() and not refresh:
        LOG.info("cached   %s", B17024_PATH.name)
        return B17024_PATH
    import acs_earnings as ae
    import acs_occupations as ao
    key = ao._read_key()
    sess = ae.make_session()
    base = f"https://api.census.gov/data/{YEAR}/acs/acs5"
    # the state list needs no key; keeping it out of the URL keeps it out of any HTTPError traceback
    r = sess.get(f"{base}?get=NAME&for=state:*", timeout=60)
    r.raise_for_status()
    states = [row[-1] for row in r.json()[1:]]
    out = ae.pull(sess, base, states, key, ",".join(["NAME"] + B17_VARS), workers, "B17024")
    if out.empty:
        raise SystemExit("B17024: nothing returned")
    out["year"] = YEAR
    for c in out.columns:
        if c[:1] == "B" and c[-1] == "E":
            out[c] = pd.to_numeric(out[c], errors="coerce")
    RAW.mkdir(parents=True, exist_ok=True)
    out.to_parquet(B17024_PATH, index=False)
    LOG.info("wrote %s (%d rows); layers: %s", B17024_PATH.name, len(out), out.sd_layer.value_counts().to_dict())
    return B17024_PATH


def cmd_fetch(args: argparse.Namespace) -> None:
    if args.only in (None, "acsed"):
        m = fetch_acsed(args.refresh)
        LOG.info("ACS-ED: %d files cached under %s", len(m["files"]), RAW)
    if args.only in (None, "b17024"):
        fetch_b17024(args.workers, args.refresh)


# --------------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------------

def _read_zip_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"no {path.name}: run `children_income.py fetch`")
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith(".txt")][0]
        return pd.read_csv(z.open(name), sep="|", dtype=str)


def _num(s: pd.Series) -> pd.Series:
    """Counts and estimates as published: thousands separators removed; '-', '**', '(X)', '250,000+' (a median in
    the open top interval) and other annotations become NaN."""
    return pd.to_numeric(s.astype(str).str.replace(",", "", regex=False), errors="coerce")


def load_acsed(it: str) -> pd.DataFrame:
    """One row per ACS-ED school-district record for population group `it`: the family-income brackets
    (families), children, children per family, published medians, child poverty."""
    inc = _read_zip_table(acsed_zip("CDP03.2", it))
    kid = _read_zip_table(acsed_zip("CDP05.1", it))
    hh = _read_zip_table(acsed_zip("CDP02.1", it))
    pv = _read_zip_table(acsed_zip("CDP03.4", it))
    if inc.LEAID.duplicated().any():
        raise SystemExit("duplicate LEAID in the ACS-ED income table")
    out = pd.DataFrame({"leaid": inc.LEAID.to_numpy(), "name": inc.Geography.to_numpy(),
                        "layer": inc.GeoId.str[:7].map(LAYER_BY_PREFIX).to_numpy()})
    out["fam"] = _num(inc.CDP03_34est).to_numpy()
    out["hh"] = _num(inc.CDP03_5est).to_numpy()
    for k, c in enumerate(INC_COLS):
        out[f"f{k}"] = _num(inc[c]).to_numpy()
    out["med_fam"] = _num(inc.CDP03_45est).to_numpy()
    out["mean_fam"] = _num(inc.CDP03_46est).to_numpy()
    for tab, cols in ((kid, {"kids": "CDP05_1est"}), (hh, {"kpf": "CDP02_9est"}),
                      (pv, {"pov517": "CDP03_67pct", "pov_u18": "CDP03_65pct"})):
        t = tab.drop_duplicates("LEAID").set_index("LEAID")
        for new, old in cols.items():
            out[new] = _num(t[old]).reindex(out.leaid).to_numpy()
    if out.layer.isna().any():
        raise SystemExit("unexpected GeoId prefix in the ACS-ED export")
    return out


def load_b17024() -> pd.DataFrame:
    """Children 6-17 by income-to-poverty bin, one row per leaid (the rule of `district_bootstrap.
    load_income_columns`: 'Remainder of <state>' rows dropped, one row per leaid, annotation codes masked)."""
    if not B17024_PATH.exists():
        raise SystemExit(f"no {B17024_PATH.name}: run `children_income.py fetch --only b17024`")
    a = pd.read_parquet(B17024_PATH)
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid")
    out = pd.DataFrame({"leaid": a.leaid.to_numpy(),
                        "layer": a.sd_layer.str.extract(r"\((\w+)\)")[0].to_numpy()})
    for j in range(12):
        x = a[B17_6_11[j]].astype(float) + a[B17_12_17[j]].astype(float)
        bad = (a[B17_6_11[j]] < 0) | (a[B17_12_17[j]] < 0)
        out[f"i{j}"] = x.where(~bad).to_numpy()
    cols = [f"i{j}" for j in range(12)]
    out.loc[out[cols].isna().any(axis=1), cols] = np.nan             # a district has all twelve bins or none
    out["kids6_17"] = out[cols].sum(axis=1, min_count=12)
    return out


def parse_us_families(path: Path) -> dict[str, Any]:
    """National 'Families' rows of the CDP03.2 preview (family-income brackets, total, median, mean)."""
    j = json.loads(path.read_text())
    if not j.get("html"):                       # the national row is not defined for iteration 102
        return {}
    rows = re.findall(r"<tr>(.*?)</tr>", j["html"], flags=re.S)
    section, fam, total, median, mean = None, [], None, None, None
    for r in rows:
        cells = [unescape(re.sub(r"<[^>]+>", "", c)).replace("\xa0", "").strip()
                 for c in re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", r, flags=re.S)]
        ind = re.search(r'indent(\d)', r)
        if not cells or not ind:
            continue
        lab, est = cells[0], cells[1] if len(cells) > 1 else ""
        if ind.group(1) == "1":
            section = lab
            if lab == "Families":
                total = float(est.replace(",", ""))
        elif section == "Families" and ind.group(1) == "2":
            if lab.startswith("Median family income"):
                median = float(est.replace(",", ""))
            elif lab.startswith("Mean family income"):
                mean = float(est.replace(",", ""))
            else:
                fam.append(float(est.replace(",", "")))
        elif section == "Families" and lab.startswith(("Median family income", "Mean family income")):
            v = float(est.replace(",", ""))
            if lab.startswith("Median"):
                median = v
            else:
                mean = v
    return dict(total=total, brackets=fam, median=median, mean=mean)


# --------------------------------------------------------------------------
# Binned income distributions: cuts, tercile weights, medians
# --------------------------------------------------------------------------

def cuts_from_counts(N: np.ndarray, edges: np.ndarray, ps: tuple[float, ...] = (1 / 3, 2 / 3)) -> np.ndarray:
    """Quantiles `ps` of a distribution known through counts in the intervals of `edges`, by linear
    interpolation inside the interval that holds the quantile (uniform density within an interval)."""
    N = np.asarray(N, float)
    cum = np.cumsum(N) / N.sum()
    out = []
    for p in ps:
        j = int(min(np.searchsorted(cum, p), len(N) - 1))
        below = cum[j - 1] if j > 0 else 0.0
        lo, hi = edges[j], edges[j + 1]
        if not np.isfinite(hi) or cum[j] <= below:
            out.append(float(lo))
        else:
            out.append(float(lo + (hi - lo) * (p - below) / (cum[j] - below)))
    return np.array(out)


def tercile_weights(cuts: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """(n_bins x 3): the share of each interval that lies in tercile t (cuts = the two interior cut
    points).  An open top interval goes entirely to the tercile that holds its lower edge."""
    nb = len(edges) - 1
    c = np.concatenate([[-np.inf], np.asarray(cuts, float), [np.inf]])
    W = np.zeros((nb, 3))
    for k in range(nb):
        lo, hi = edges[k], edges[k + 1]
        for t in range(3):
            if not np.isfinite(hi):
                W[k, t] = 1.0 if c[t] <= lo < c[t + 1] else 0.0
            else:
                W[k, t] = max(0.0, min(hi, c[t + 1]) - max(lo, c[t])) / (hi - lo)
    return W


def bracket_weights(N: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, dict[str, Any]]:
    """Terciles without interpolation: the two interior interval edges whose cumulative national share is
    closest to 1/3 and 2/3 (a 0/1 matrix) and the national shares of the three groups that result."""
    N = np.asarray(N, float)
    share = N / N.sum()
    cum = np.cumsum(share)[:-1]                                # cumulative share at the upper edge of interval k
    j1 = int(np.argmin(np.abs(cum - 1 / 3)))
    j2 = int(np.argmin(np.abs(cum - 2 / 3)))
    if j2 <= j1:
        j2 = min(j1 + 1, len(cum) - 1)
    W = np.zeros((len(N), 3))
    for k in range(len(N)):
        W[k, 0 if k <= j1 else (1 if k <= j2 else 2)] = 1.0
    return W, dict(cut_edges=[float(edges[j1 + 1]), float(edges[j2 + 1])],
                   group_shares=[float(share[: j1 + 1].sum()), float(share[j1 + 1: j2 + 1].sum()),
                                 float(share[j2 + 1:].sum())])


def grouped_median(counts: np.ndarray, edges: np.ndarray) -> float:
    """Median of a distribution known through counts in the intervals of `edges`, by linear interpolation
    inside the interval that holds the 50th percentile (NaN if that is the open top interval or if there is
    no observation: the rule of `district_bootstrap.grouped_median`)."""
    c = np.asarray(counts, float)
    tot = c.sum()
    if not np.isfinite(tot) or tot <= 0 or np.isnan(c).any():
        return float("nan")
    cum = np.cumsum(c)
    j = int(min(np.searchsorted(cum, tot / 2.0), len(c) - 1))
    if not np.isfinite(edges[j + 1]):
        return float("nan")
    below = cum[j - 1] if j > 0 else 0.0
    return float(edges[j] + (edges[j + 1] - edges[j]) * (tot / 2.0 - below) / c[j])


@dataclass
class National:
    """The national distribution of one population group and weighting: interval totals, shares, cut points and
    the tercile weight matrices (interpolated and bracket-based)."""
    key: str
    edges: np.ndarray
    labels: list[str]
    N: np.ndarray
    n_records: int
    cuts: np.ndarray
    W: np.ndarray
    W_bracket: np.ndarray
    bracket_info: dict[str, Any]

    @classmethod
    def build(cls, key: str, N: np.ndarray, edges: np.ndarray, labels: list[str], n_records: int) -> "National":
        cuts = cuts_from_counts(N, edges)
        Wb, info = bracket_weights(N, edges)
        return cls(key, edges, labels, np.asarray(N, float), n_records, cuts, tercile_weights(cuts, edges), Wb, info)

    def summary(self) -> dict[str, Any]:
        sh = self.N / self.N.sum()
        return dict(total=float(self.N.sum()), n_records=self.n_records, bin_shares=sh.tolist(),
                    cut_points=self.cuts.tolist(), bracket_cut_edges=self.bracket_info["cut_edges"],
                    bracket_group_shares=self.bracket_info["group_shares"],
                    tercile_shares_check=(sh @ self.W).tolist())


# --------------------------------------------------------------------------
# District-level inputs of one specification
# --------------------------------------------------------------------------

def district_counts(P: pd.DataFrame, it: str, weight: str) -> pd.DataFrame:
    """The ten family-income interval counts of every district of `P` (columns `a<it>_*` of the merged panel):
    `child` = families times the district's children per family (so that the interval counts add up to the
    district's children), `family` = families as published.  A district counts only if all ten intervals are
    published and it has families (and, for `child`, children); the others are NaN.  Also: the total weight
    `wt`, the published mean and median family income, and the child poverty rate (5-17, percent)."""
    p = f"a{it}_"
    F = P[[p + f"f{k}" for k in range(NB)]].to_numpy(float)
    fam = P[p + "fam"].to_numpy(float)
    kids = P[p + "kids"].to_numpy(float)
    ok = np.isfinite(F).all(axis=1) & np.isfinite(fam) & (fam > 0)
    if weight == "child":
        ok &= np.isfinite(kids) & (kids > 0)
        s = np.where(ok, kids / np.where(fam > 0, fam, np.nan), np.nan)
    elif weight == "family":
        s = np.ones(len(P))
    else:
        raise ValueError(weight)
    N = F * s[:, None]
    N[~ok] = np.nan
    out = pd.DataFrame(N, columns=[f"n{k}" for k in range(NB)], index=P.index)
    out["wt"] = np.where(ok, np.nansum(np.where(np.isfinite(N), N, 0.0), axis=1), np.nan)
    out["kids"] = np.where(ok, kids, np.nan)
    out["mean_fam"] = np.where(ok, P[p + "mean_fam"].to_numpy(float), np.nan)
    out["med_fam"] = np.where(ok, P[p + "med_fam"].to_numpy(float), np.nan)
    out["pov517"] = np.where(ok & np.isfinite(kids) & (kids > 0), P[p + "pov517"].to_numpy(float), np.nan)
    out["ok"] = ok
    return out


def cz_sums(df: pd.DataFrame, czone: np.ndarray, loc_id: np.ndarray, cols: list[str], czones: np.ndarray
            ) -> np.ndarray:
    """(n_cz, 2, len(cols)): column sums over the valid rows of each (CZ, location); NaN where a (CZ, location)
    has no valid row, which removes the CZ from every statistic that needs both locations."""
    d = df[cols].copy()
    d["_cz"], d["_loc"] = czone, loc_id
    sums = d.groupby(["_cz", "_loc"])[cols].sum(min_count=1)
    out = np.full((len(czones), 2, len(cols)), np.nan)
    pos = {float(c): i for i, c in enumerate(czones)}
    for (c, l), row in sums.iterrows():
        i = pos.get(float(c))
        if i is not None and int(l) in (1, 2):
            out[i, int(l) - 1] = row.to_numpy(float)
    return out


def cz_pairs(df: pd.DataFrame, czone: np.ndarray, loc_id: np.ndarray, num: np.ndarray, den: np.ndarray,
             czones: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(n_cz, 2) sums of `num` and of `den` over the rows where both are finite and den > 0."""
    ok = np.isfinite(num) & np.isfinite(den) & (den > 0)
    t = pd.DataFrame({"num": np.where(ok, num, np.nan), "den": np.where(ok, den, np.nan)}, index=df.index)
    a = cz_sums(t, czone, loc_id, ["num", "den"], czones)
    return a[:, :, 0], a[:, :, 1]


# --------------------------------------------------------------------------
# Bootstrap: the multiplicity matrix of district_settings.py
# --------------------------------------------------------------------------

class Boot:
    """The whole-CZ multiplicity matrix (B x n_CZ) of `district_bootstrap.draw_counts`, shared by every
    statistic of a sample, and the CZ-weighted mean / ratio statistics on it."""

    def __init__(self, n_cz: int, B: int, seed: int, w: np.ndarray, counts: np.ndarray | None = None):
        self.counts = db.draw_counts(n_cz, B, seed) if counts is None else counts
        self.w = np.asarray(w, float)
        self.B, self.seed, self.n_cz = B, seed, n_cz

    def with_weights(self, w: np.ndarray) -> "Boot":
        return Boot(self.n_cz, self.B, self.seed, w, self.counts)

    def wmean(self, x: np.ndarray) -> dict[str, Any]:
        """CZ-weighted mean of x over the CZs where it is observed (value, bootstrap SE, 95% interval, n)."""
        x = np.asarray(x, float)
        if not (np.isfinite(x) & (self.w > 0)).any():
            return dict(value=float("nan"), se=float("nan"), ci95_lo=float("nan"), ci95_hi=float("nan"), n_cz=0,
                        draws=np.full(self.B, np.nan))
        s = db.boot_wmean(x, self.w, self.counts)
        out = ds.summarise(s)
        out["draws"] = s["draws"]
        return out

    def ratio(self, num: np.ndarray, den: np.ndarray) -> dict[str, Any]:
        """sum(num) / sum(den) over the CZs that have both, with draws (the weights are the sizes of num/den)."""
        num, den = np.asarray(num, float), np.asarray(den, float)
        m = np.isfinite(num) & np.isfinite(den) & (den > 0)
        if not m.any():
            return dict(value=float("nan"), se=float("nan"), ci95_lo=float("nan"), ci95_hi=float("nan"), n_cz=0,
                        draws=np.full(self.B, np.nan))
        n_, d_ = np.where(m, num, 0.0), np.where(m, den, 0.0)
        est = n_.sum() / d_.sum()
        with np.errstate(invalid="ignore", divide="ignore"):
            draws = (self.counts @ n_) / (self.counts @ d_)
        s = dict(est=est, draws=draws, n=int(m.sum()))
        out = ds.summarise(s)
        out["draws"] = draws
        return out


# --------------------------------------------------------------------------
# Statistics of one sample and specification
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Spec:
    key: str
    label: str
    it: str = PRIMARY_IT
    weight: str = "child"          # child | family
    layers: str = "all"            # all | unified
    nyc: bool = True               # False: New York City without its ACS-ED record
    detroit: bool = True           # False: Detroit without its record (no id crosswalk)
    ny_cz: bool = True             # False: the New York CZ leaves the sample
    cuts: str = "national"         # national | bracket | cz
    cz_weights: str = "pupils"     # pupils | equal
    universe: str = "partition"    # partition (city vs suburb) | rest (city vs every other district of the CZ)


SPECS: list[Spec] = [
    Spec("primary", "primary: relevant children (102), children by scaling, national terciles, all layers"),
    Spec("public", "relevant children enrolled in public school (105)", it="105"),
    Spec("enrolled", "relevant children enrolled in public or private school (104)", it="104"),
    Spec("family", "families as published (no scaling to children), group 102", weight="family"),
    Spec("unified", "unified districts only", layers="unified"),
    Spec("nyc_no_record", "New York City without its ACS-ED record", nyc=False),
    Spec("detroit_no_record", "Detroit without its ACS-ED record (no id crosswalk)", detroit=False),
    Spec("drop_ny_cz", "New York CZ left out", ny_cz=False),
    Spec("bracket_cuts", "terciles at the bracket edges closest to 1/3 and 2/3 (no interpolation)", cuts="bracket"),
    Spec("within_cz_cuts", "terciles of each CZ's own distribution (cuts in every CZ)", cuts="cz"),
    Spec("equal_cz_weights", "every CZ weighted equally", cz_weights="equal"),
    Spec("city_vs_rest", "city against every other district of the CZ (suburb, town and rural), CZs where the city holds "
         "at least 10% of the pupils", universe="rest"),
]
PRIMARY = SPECS[0]
HEADLINE = ("city_share_t1", "city_share_t2", "city_share_t3", "city_share_gradient", "city_share_children",
            "pooled_city_share_gradient", "gap_median_family_income", "gap_mean_family_income")


def apply_spec_mask(P: pd.DataFrame, D: pd.DataFrame, spec: Spec) -> pd.DataFrame:
    """Rows excluded by the specification become invalid (NaN counts, `ok` False)."""
    D = D.copy()
    drop = np.zeros(len(P), bool)
    if spec.layers == "unified":
        drop |= (P["a_layer"] != "unified").to_numpy()
    if not spec.nyc:
        drop |= (P["leaid"] == NYC_LEAID).to_numpy()
    if not spec.detroit:
        drop |= (P["leaid"] == DETROIT[1]).to_numpy()
    if drop.any():
        cols = [c for c in D.columns if c != "ok"]
        D.loc[drop, cols] = np.nan
        D.loc[drop, "ok"] = False
    return D


def core_stats(cnt: np.ndarray, W: np.ndarray, boot: Boot) -> dict[str, dict[str, Any]]:
    """The tercile statistics of an (n_cz, 2, n_bins) array of interval counts (location 1 = city, 2 = suburb)
    and tercile weights W ((n_bins, 3), or (n_cz, n_bins, 3) for CZ-specific cuts)."""
    n = cnt.shape[0]
    Wc = np.broadcast_to(W, (n,) + W.shape[-2:]) if W.ndim == 2 else W
    with np.errstate(invalid="ignore", divide="ignore"):
        T = np.einsum("clk,ckt->clt", cnt, Wc)                       # children in tercile t, by CZ and location
        tot = T.sum(axis=2)
        city = T[:, 0, :] / (T[:, 0, :] + T[:, 1, :])
    out: dict[str, dict[str, Any]] = {}
    for t in range(3):
        out[f"city_share_t{t + 1}"] = boot.wmean(city[:, t])
    out["city_share_gradient"] = boot.wmean(city[:, 2] - city[:, 0])
    out["city_share_mid_minus_bottom"] = boot.wmean(city[:, 1] - city[:, 0])
    out["city_share_top_minus_mid"] = boot.wmean(city[:, 2] - city[:, 1])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["city_share_children"] = boot.wmean(tot[:, 0] / (tot[:, 0] + tot[:, 1]))
    pooled = [boot.ratio(T[:, 0, t], T[:, 0, t] + T[:, 1, t]) for t in range(3)]
    for t in range(3):
        out[f"pooled_city_share_t{t + 1}"] = pooled[t]
    # pooled gradient: difference of the two pooled ratios on the same draws
    m = np.isfinite(city).all(axis=1)
    est = pooled[2]["value"] - pooled[0]["value"]
    dr = pooled[2]["draws"] - pooled[0]["draws"]
    s = ds.summarise(dict(est=est, draws=dr, n=int(m.sum())))
    s["draws"] = dr
    out["pooled_city_share_gradient"] = s
    with np.errstate(invalid="ignore", divide="ignore"):
        comp = T / tot[:, :, None]
    for t in range(3):
        out[f"comp_city_t{t + 1}"] = boot.wmean(comp[:, 0, t])
        out[f"comp_suburb_t{t + 1}"] = boot.wmean(comp[:, 1, t])
    return out


def level_stats(D: pd.DataFrame, P: pd.DataFrame, cnt: np.ndarray, edges: np.ndarray, czones: np.ndarray,
                boot: Boot) -> dict[str, dict[str, Any]]:
    """Gaps (suburb minus city) in the median, the mean and the average published median of family income, and the
    child poverty rate; the location levels behind them."""
    out: dict[str, dict[str, Any]] = {}
    n = cnt.shape[0]
    med = np.array([[grouped_median(cnt[c, l], edges) for l in (0, 1)] for c in range(n)])
    out["gap_median_family_income"] = boot.wmean(np.log(med[:, 1] / med[:, 0]))
    for l, nm in ((0, "city"), (1, "suburb")):
        both = np.where(np.isfinite(med).all(axis=1), med[:, l], np.nan)
        out[f"median_family_income_{nm}"] = boot.wmean(both)
    cz, loc = P["czone"].to_numpy(), P["loc_id"].to_numpy()
    wt = D["wt"].to_numpy(float)
    # exact aggregate mean: sum over districts of wt * published mean / sum of wt
    nm_, dn_ = cz_pairs(D, cz, loc, wt * D["mean_fam"].to_numpy(float), wt, czones)
    mean = nm_ / dn_
    out["gap_mean_family_income"] = boot.wmean(np.log(mean[:, 1] / mean[:, 0]))
    for l, nm in ((0, "city"), (1, "suburb")):
        out[f"mean_family_income_{nm}"] = boot.wmean(np.where(np.isfinite(mean).all(axis=1), mean[:, l], np.nan))
    # average of district log median family incomes (the `ln_` convention of the housing scripts)
    medd = D["med_fam"].to_numpy(float)
    nl, dl = cz_pairs(D, cz, loc, wt * np.log(np.where(medd > 0, medd, np.nan)), wt, czones)
    lnm = nl / dl
    out["gap_lnmedian_published"] = boot.wmean(lnm[:, 1] - lnm[:, 0])
    # child poverty (related children 5-17 below the poverty level, percent), weighted by children
    kids = D["kids"].to_numpy(float)
    npv, dpv = cz_pairs(D, cz, loc, kids * D["pov517"].to_numpy(float) / 100.0, kids, czones)
    pov = npv / dpv
    out["gap_child_poverty"] = boot.wmean(pov[:, 1] - pov[:, 0])
    for l, nm in ((0, "city"), (1, "suburb")):
        out[f"child_poverty_{nm}"] = boot.wmean(np.where(np.isfinite(pov).all(axis=1), pov[:, l], np.nan))
    return out


def distribution_stats(cnt: np.ndarray, boot: Boot, labels: list[str]) -> dict[str, dict[str, Any]]:
    """Interval shares within each location (CZ-weighted means), their suburb-minus-city difference, and the city
    share of the children in each interval."""
    out: dict[str, dict[str, Any]] = {}
    n, _, nb = cnt.shape
    tot = cnt.sum(axis=2, keepdims=True)
    sh = cnt / tot
    for k in range(nb):
        out[f"bin{k}_city"] = boot.wmean(sh[:, 0, k])
        out[f"bin{k}_suburb"] = boot.wmean(sh[:, 1, k])
        out[f"bin{k}_diff"] = boot.wmean(sh[:, 1, k] - sh[:, 0, k])
        out[f"bin{k}_city_share"] = boot.wmean(cnt[:, 0, k] / (cnt[:, 0, k] + cnt[:, 1, k]))
    return out


def spec_inputs(part: pd.DataFrame, g: pd.DataFrame, spec: Spec) -> dict[str, Any]:
    """District counts and the CZ x location x interval array for `spec` on the partition `part`."""
    czones = g.czone.to_numpy(float)
    D = apply_spec_mask(part, district_counts(part, spec.it, spec.weight), spec)
    cnt = cz_sums(D, part["czone"].to_numpy(), part["loc_id"].to_numpy(), [f"n{k}" for k in range(NB)], czones)
    ny_idx = None
    ny_rows = part.loc[part["leaid"] == NYC_LEAID, "czone"]
    if len(ny_rows):
        hit = np.where(czones == float(ny_rows.iloc[0]))[0]
        ny_idx = int(hit[0]) if len(hit) else None
    if not spec.ny_cz and ny_idx is not None:
        cnt[ny_idx] = np.nan
    return dict(D=D, cnt=cnt, czones=czones, ny_idx=ny_idx)


def run_spec(part: pd.DataFrame, g: pd.DataFrame, boot: Boot, nat: dict[str, National], spec: Spec,
             full: bool = False) -> dict[str, Any]:
    """Statistics of one specification: the tercile block always; the level, median / mean / poverty and interval
    blocks when `full` (the primary specification)."""
    inp = spec_inputs(part, g, spec)
    cnt, D, czones = inp["cnt"], inp["D"], inp["czones"]
    nt = nat[f"{spec.it}:{spec.weight}"]
    bt = boot if spec.cz_weights == "pupils" else boot.with_weights(np.ones(len(g)))
    if spec.cuts == "national":
        W = nt.W
    elif spec.cuts == "bracket":
        W = nt.W_bracket
    else:                                                          # within-CZ terciles
        W = np.zeros((len(g), NB, 3))
        for c in range(len(g)):
            tot = cnt[c].sum(axis=0)
            W[c] = tercile_weights(cuts_from_counts(tot, INC_EDGES), INC_EDGES) if np.isfinite(tot).all() else np.nan
    stats = core_stats(cnt, W, bt)
    med = np.array([[grouped_median(cnt[c, l], INC_EDGES) for l in (0, 1)] for c in range(len(g))])
    stats["gap_median_family_income"] = bt.wmean(np.log(med[:, 1] / med[:, 0]))
    if full:
        stats.update(level_stats(D, part, cnt, INC_EDGES, inp["czones"], bt))
        stats.update(distribution_stats(cnt, bt, INC_LABELS))
    else:
        nm_, dn_ = cz_pairs(D, part["czone"].to_numpy(), part["loc_id"].to_numpy(),
                            D["wt"].to_numpy(float) * D["mean_fam"].to_numpy(float), D["wt"].to_numpy(float), czones)
        mean_pair = nm_ / dn_
        stats["gap_mean_family_income"] = bt.wmean(np.log(mean_pair[:, 1] / mean_pair[:, 0]))
        kids = D["kids"].to_numpy(float)
        npv, dpv = cz_pairs(D, part["czone"].to_numpy(), part["loc_id"].to_numpy(),
                            kids * D["pov517"].to_numpy(float) / 100.0, kids, czones)
        stats["gap_child_poverty"] = bt.wmean(npv[:, 1] / dpv[:, 1] - npv[:, 0] / dpv[:, 0])
    return dict(stats=stats, cnt=cnt, D=D, n_districts=int(D["ok"].sum()),
                children=float(np.nansum(D["kids"].to_numpy(float))) if spec.weight == "child" else None)


# --------------------------------------------------------------------------
# Cross-check: Census B17024, children 6-17 by income-to-poverty ratio
# --------------------------------------------------------------------------

def run_b17024(part: pd.DataFrame, g: pd.DataFrame, boot: Boot, nat: National) -> dict[str, Any]:
    """Same statistics on the Census API table B17024 (children 6-17 by ratio of income to the poverty level,
    twelve bins), unified and elementary district records only (they tile the states, so every child is counted
    once; the secondary layer overlays them), national terciles of the IPR distribution of those children."""
    czones = g.czone.to_numpy(float)
    cols = [f"b_i{j}" for j in range(12)]
    N = part[cols].to_numpy(float)
    ok = np.isfinite(N).all(axis=1) & part["b_layer"].isin(["unified", "elementary"]).to_numpy()
    N = np.where(ok[:, None], N, np.nan)
    D = pd.DataFrame(N, columns=[f"n{j}" for j in range(12)], index=part.index)
    cnt = cz_sums(D, part["czone"].to_numpy(), part["loc_id"].to_numpy(), list(D.columns), czones)
    stats = core_stats(cnt, nat.W, boot)
    med = np.array([[grouped_median(cnt[c, l], IPR_EDGES) for l in (0, 1)] for c in range(len(g))])
    stats["gap_median_ipr"] = boot.wmean(np.log(med[:, 1] / med[:, 0]))
    tot = cnt.sum(axis=2)
    pov = cnt[:, :, :3].sum(axis=2) / tot                               # IPR below 1.00
    stats["gap_poverty_ipr"] = boot.wmean(pov[:, 1] - pov[:, 0])
    for l, nm in ((0, "city"), (1, "suburb")):
        stats[f"poverty_ipr_{nm}"] = boot.wmean(np.where(np.isfinite(pov).all(axis=1), pov[:, l], np.nan))
    return dict(stats=stats, cnt=cnt, n_districts=int(ok.sum()),
                children=float(np.nansum(N)))


# --------------------------------------------------------------------------
# National distributions, merge, checks
# --------------------------------------------------------------------------

def build_nationals(acs: dict[str, pd.DataFrame], b17: pd.DataFrame | None) -> dict[str, National]:
    """National interval totals of every group and weighting from the sum of all district records."""
    nat: dict[str, National] = {}
    for it, a in acs.items():
        P = a.add_prefix(f"a{it}_")
        for weight in ("child", "family"):
            D = district_counts(P, it, weight)
            N = D[[f"n{k}" for k in range(NB)]].sum().to_numpy()
            nat[f"{it}:{weight}"] = National.build(f"{it}:{weight}", N, INC_EDGES, INC_LABELS, int(D["ok"].sum()))
    if b17 is not None:
        cols = [f"i{j}" for j in range(12)]
        ok = b17.layer.isin(["unified", "elementary"]) & b17[cols].notna().all(axis=1)
        nat["b17024"] = National.build("b17024", b17.loc[ok, cols].sum().to_numpy(), IPR_EDGES, IPR_LABELS,
                                       int(ok.sum()))
    return nat


def merge_ext(A: pd.DataFrame, acs: dict[str, pd.DataFrame], b17: pd.DataFrame | None) -> pd.DataFrame:
    """ACS-ED columns (`a<it>_*`, layer `a_layer`) and the B17024 bins (`b_i*`, `b_layer`) onto the audited panel
    by `leaid`, after the NYC collapse: the collapsed row carries the Chancellor's id 3620580, under which both
    sources file New York City."""
    n0 = len(A)
    X = A.copy()
    src, dst = DETROIT
    fix = dst in set(A.leaid)

    def refile(d: pd.DataFrame) -> pd.DataFrame:
        d = d.copy()
        if fix:
            d["leaid"] = d["leaid"].where(d["leaid"] != src, dst)
        return d

    for it, a in acs.items():
        b = refile(a).drop(columns=["name", "layer"]).add_prefix(f"a{it}_").rename(columns={f"a{it}_leaid": "leaid"})
        X = X.merge(b, on="leaid", how="left", validate="many_to_one")
    X = X.merge(refile(acs[PRIMARY_IT])[["leaid", "layer"]].rename(columns={"layer": "a_layer"}), on="leaid",
                how="left", validate="many_to_one")
    if b17 is not None:
        b = refile(b17).rename(columns={"layer": "b_layer", **{f"i{j}": f"b_i{j}" for j in range(12)},
                                        "kids6_17": "b_kids"})
        X = X.merge(b, on="leaid", how="left", validate="many_to_one")
    if len(X) != n0:
        raise SystemExit("merge of the ACS-ED / B17024 tables changed the number of districts")
    X.attrs.update(A.attrs)
    return X


def us_row_check(acs: dict[str, pd.DataFrame], nat: dict[str, National]) -> dict[str, Any]:
    """National 'Families' rows of the EDGE preview against the sum of the district records (group 102 has no
    national row: for the nation, 'relevant' children are all children)."""
    out: dict[str, Any] = {}
    for it in ITERS:
        p = us_json("CDP03.2", it)
        us = parse_us_families(p) if p.exists() else {}
        if not us or not us.get("brackets"):
            out[it] = dict(available=False)
            continue
        us_sh = np.array(us["brackets"]) / np.sum(us["brackets"])
        a = acs[it]
        F = a[[f"f{k}" for k in range(NB)]].sum().to_numpy()
        d_sh = F / F.sum()
        out[it] = dict(available=True, us_families=us["total"], districts_families=float(a.fam.sum()),
                       ratio_districts_to_us=float(a.fam.sum() / us["total"]),
                       max_abs_share_diff=float(np.abs(us_sh - d_sh).max()), us_shares=us_sh.tolist(),
                       district_shares=d_sh.tolist(), us_median=us.get("median"), us_mean=us.get("mean"),
                       district_children=float(a.kids.sum()))
    return out


def interpolation_check(a: pd.DataFrame) -> dict[str, Any]:
    """Median family income interpolated from a district's own ten intervals against its published median
    (districts with both; the open top interval gives no median)."""
    F = a[[f"f{k}" for k in range(NB)]].to_numpy(float)
    ok = np.isfinite(F).all(axis=1) & a.med_fam.notna().to_numpy() & (a.fam > 0).to_numpy()
    gm = np.array([grouped_median(r, INC_EDGES) for r in F[ok]])
    pub = a.med_fam.to_numpy(float)[ok]
    err = np.log(gm) - np.log(pub)
    wt = a.fam.to_numpy(float)[ok]
    m = np.isfinite(err)
    return dict(n=int(m.sum()), mean_log_error=float(err[m].mean()),
                weighted_mean_log_error=float(np.average(err[m], weights=wt[m])),
                sd_log_error=float(err[m].std(ddof=1)), p01=float(np.percentile(err[m], 1)),
                p99=float(np.percentile(err[m], 99)),
                n_published_median_missing=int(a.med_fam.isna().sum()))


def rounding_check(a: pd.DataFrame, part: pd.DataFrame, it: str) -> dict[str, Any]:
    """Counts are rounded to 5 and 1-4 appear as 4: the share of such cells and their mass."""
    F = a[[f"f{k}" for k in range(NB)]].to_numpy(float)
    F = F[np.isfinite(F).all(axis=1)]
    P = part[[f"a{it}_f{k}" for k in range(NB)]].to_numpy(float)
    P = P[np.isfinite(P).all(axis=1)]
    return dict(share_cells_equal_4=float((F == 4).mean()), share_cells_1_3=float(((F >= 1) & (F <= 3)).mean()),
                share_cells_multiple_of_5=float((F % 5 == 0).mean()), mass_of_4_cells=float(F[F == 4].sum() / F.sum()),
                partition_share_cells_equal_4=float((P == 4).mean()),
                partition_mass_of_4_cells=float(P[P == 4].sum() / P.sum()))


def kids_per_family_check(acs: dict[str, pd.DataFrame]) -> dict[str, Any]:
    """The scaling factor (children / families of the group's table) against the published average number of
    children in families with children (CDP02_9est): national means, and the correlation across districts."""
    a = acs[PRIMARY_IT].copy()
    a["r"] = a.kids / a.fam
    d = a.dropna(subset=["r", "kpf"])
    big = d[d.fam >= 500]
    return dict(national_mean={it: float(acs[it].kids.sum() / acs[it].fam.sum()) for it in ITERS},
                corr_all=float(d[["r", "kpf"]].corr().iloc[0, 1]), n_all=int(len(d)),
                corr_ge500=float(big[["r", "kpf"]].corr().iloc[0, 1]), n_ge500=int(len(big)),
                mean_r_ge500=float(big.r.mean()), mean_kpf_ge500=float(big.kpf.mean()))


def coverage(part: pd.DataFrame, g: pd.DataFrame, nat: dict[str, National], acs: dict[str, pd.DataFrame]
             ) -> dict[str, Any]:
    """Size of the partition and of the matched ACS-ED sample."""
    pu = float(part.pupils.sum())
    out: dict[str, Any] = dict(n_cz=int(len(g)), n_districts=int(len(part)), pupils=pu)
    wn = g.pupils.to_numpy(float) / g.pupils.sum()
    out["n_eff_cz"] = float(1.0 / (wn ** 2).sum())
    out["top5_weight_share"] = float(np.sort(wn)[::-1][:5].sum())
    for it in ITERS:
        D = district_counts(part, it, "child")
        ok = D["ok"].to_numpy()
        byloc = {}
        for l, nm in ((1, "city"), (2, "suburb")):
            m = (part.loc_id == l).to_numpy()
            byloc[nm] = dict(n_districts=int(m.sum()), n_matched=int((m & ok).sum()),
                             pupil_share_matched=float(part.pupils[m & ok].sum() / part.pupils[m].sum()),
                             pupils=float(part.pupils[m].sum()),
                             children=float(np.nansum(D["kids"].to_numpy(float)[m & ok])))
        out[f"matched_{it}"] = dict(n_districts=int(ok.sum()), pupil_share=float(part.pupils[ok].sum() / pu),
                                    by_location=byloc,
                                    layers={k: int(v) for k, v in
                                            part.loc[ok, "a_layer"].value_counts().items()})
    unm = part[district_counts(part, PRIMARY_IT, "child")["ok"].to_numpy() == False]      # noqa: E712
    out["unmatched"] = dict(n=int(len(unm)), pupils=float(unm.pupils.sum()),
                            largest=[dict(name=str(r.get("lea_name", "")), pupils=float(r.pupils), loc=int(r.loc_id))
                                     for _, r in unm.sort_values("pupils", ascending=False).head(5).iterrows()])
    nyc = part[part.leaid == NYC_LEAID]
    out["nyc_in_partition"] = bool(len(nyc))
    out["nyc_pupils"] = float(nyc.pupils.sum()) if len(nyc) else 0.0
    out["nyc_weight_share_of_cz"] = None
    return out


# --------------------------------------------------------------------------
# One sample: every specification, the cross-check
# --------------------------------------------------------------------------

def run_sample(g: pd.DataFrame, part: pd.DataFrame, B: int, seed: int, nat: dict[str, National],
               part_alt: pd.DataFrame | None = None) -> dict[str, Any]:
    """All specifications on the CZs of `g` (the full sample or a hold-out half), on the whole-CZ multiplicity
    matrix of that sample (n CZs, B, seed).  `part_alt` is the partition of the `city_vs_rest` variant.  Nothing
    is logged here."""
    boot = Boot(len(g), B, seed, g.pupils.to_numpy(float))
    specs = {}
    for sp in SPECS:
        if sp.universe == "rest" and part_alt is None:
            continue
        pp = part_alt if sp.universe == "rest" else part
        specs[sp.key] = run_spec(pp, g, boot, nat, sp, full=(sp.key == "primary"))
    b = run_b17024(part, g, boot, nat["b17024"]) if "b17024" in nat and "b_i0" in part.columns else None
    stats = specs["primary"]["stats"]
    # the gradient relative to the overall city share (reading aid, not a target)
    gr, ch = stats["city_share_gradient"], stats["city_share_children"]
    rel = gr["draws"] / ch["draws"]
    s = ds.summarise(dict(est=gr["value"] / ch["value"], draws=rel, n=gr["n_cz"]))
    s["draws"] = rel
    stats["city_share_gradient_relative"] = s
    return dict(boot=boot, specs=specs, b17024=b)


def strip_draws(o: Any) -> Any:
    """Copy without the bootstrap draws and the array payloads (cnt, D) of the run dictionaries."""
    if isinstance(o, dict):
        return {k: strip_draws(v) for k, v in o.items() if k not in ("draws", "cnt", "D", "boot")}
    if isinstance(o, (list, tuple)):
        return [strip_draws(v) for v in o]
    return o


def draws_frame(stats: dict[str, dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame({k: v["draws"] for k, v in stats.items() if "draws" in v})


# --------------------------------------------------------------------------
# Child poverty cross-check against the SAIPE-based `poverty_gap` of district_settings.py
# --------------------------------------------------------------------------

def poverty_check(part: pd.DataFrame, g: pd.DataFrame, boot: Boot) -> dict[str, Any]:
    """Suburb-minus-city gap in the child poverty rate (5-17) from three sources on the same districts: ACS-ED
    (related children 5-17, weighted by children), Census B17024 (children 6-17 below IPR 1.00), and SAIPE
    (`pop_5_17_poverty / pop_5_17`, the `poverty_gap` of district_settings.py).  The matched sample is the
    districts that have all three (New York City, whose SAIPE record is filed under another id, drops out)."""
    czones = g.czone.to_numpy(float)
    cz, loc = part["czone"].to_numpy(), part["loc_id"].to_numpy()
    D = district_counts(part, PRIMARY_IT, "child")
    kids = D["kids"].to_numpy(float)
    acs_num = kids * D["pov517"].to_numpy(float) / 100.0
    cols = [f"b_i{j}" for j in range(12)]
    Nb = part[cols].to_numpy(float)
    okb = np.isfinite(Nb).all(axis=1) & part["b_layer"].isin(["unified", "elementary"]).to_numpy()
    b_den = np.where(okb, Nb.sum(axis=1), np.nan)
    b_num = np.where(okb, Nb[:, :3].sum(axis=1), np.nan)
    s_den = part["pop_5_17"].to_numpy(float)
    s_num = part["pop_5_17_poverty"].to_numpy(float)
    have = (np.isfinite(acs_num) & np.isfinite(kids) & (kids > 0) & np.isfinite(b_den) & (b_den > 0)
            & np.isfinite(s_den) & (s_den > 0) & np.isfinite(s_num))
    out: dict[str, Any] = dict(n_districts_matched=int(have.sum()),
                               pupil_share_matched=float(part.pupils[have].sum() / part.pupils.sum()))
    gaps = {}
    for key, num, den, lab in (("acsed", acs_num, kids, "ACS-ED, related children 5-17 (CDP03_67pct) x children"),
                               ("b17024", b_num, b_den, "Census B17024, children 6-17 with IPR below 1.00"),
                               ("saipe", s_num, s_den, "SAIPE, children 5-17 in poverty / children 5-17")):
        nm, dn = cz_pairs(D, cz, loc, np.where(have, num, np.nan), np.where(have, den, np.nan), czones)
        rate = nm / dn
        st = boot.wmean(rate[:, 1] - rate[:, 0])
        st["label"] = lab
        st["city"] = boot.wmean(np.where(np.isfinite(rate).all(axis=1), rate[:, 0], np.nan))["value"]
        st["suburb"] = boot.wmean(np.where(np.isfinite(rate).all(axis=1), rate[:, 1], np.nan))["value"]
        gaps[key] = st
    out["gaps"] = gaps
    return out


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

CORR_TARGETS = ["gap_pupils", "gap_salary", "gap_teachers_pp", "gap_income"]
CORR_NEW = ["city_share_gradient", "city_share_t1", "city_share_t2", "city_share_t3", "gap_median_family_income",
            "gap_mean_family_income", "gap_child_poverty", "comp_city_t1", "comp_suburb_t3"]


def reproduction_check(std: dict[str, Any], full_run: bool) -> list[dict[str, Any]]:
    rows = []
    for k in ds.TARGET_KEYS:
        s = std["stats"][k]
        row = dict(moment=k, value=s["value"], se=s["se"], n_cz=s["n_cz"], reported_value=ds.REPORTED[k][0],
                   reported_se=ds.REPORTED[k][1], ok_value=round(s["value"], 4) == ds.REPORTED[k][0],
                   ok_se=(round(s["se"], 4) == ds.REPORTED[k][1]) if full_run else None)
        rows.append(row)
        LOG.info("reproduction %-16s %.4f (%.4f) reported %.4f (%.4f) %s", k, s["value"], s["se"],
                 row["reported_value"], row["reported_se"],
                 "ok" if row["ok_value"] and row["ok_se"] in (True, None) else "MISMATCH")
    if not all(r["ok_value"] and r["ok_se"] in (True, None) for r in rows):
        raise SystemExit("the four locale targets do not reproduce the reported values; stopping")
    return rows


def check_holdout(ho: dict[str, Any], indir: Path) -> dict[str, Any]:
    """The split drawn here against the assignment stored by district_settings.py (compared, not printed)."""
    path = indir / "estimates" / "settings" / "district_settings.json"
    if not path.exists():
        LOG.warning("no %s: hold-out assignment not compared", path.name)
        return dict(compared=False)
    st = json.loads(path.read_text())["holdout"]
    same = (sorted(int(c) for c in st["A_cz"]) == sorted(int(c) for c in ho["A"])
            and sorted(int(c) for c in st["B_cz"]) == sorted(int(c) for c in ho["B"]) and st["seed"] == ho["seed"])
    if not same:
        raise SystemExit("the hold-out split differs from the one stored in district_settings.json; stopping")
    LOG.info("hold-out split identical to district_settings.json (seed %d, %d / %d CZs)", ho["seed"],
             len(ho["A"]), len(ho["B"]))
    return dict(compared=True, identical=True, seed=ho["seed"])


def write_sealed_b(path: Path, gB: pd.DataFrame, partB: pd.DataFrame, B: int, seed: int, ho: dict[str, Any],
                   nat: dict[str, National], partB_alt: pd.DataFrame | None = None) -> None:
    """Every half-B statistic goes to `path` and nowhere else: nothing computed here is logged or returned."""
    run = run_sample(gB, partB, B, seed, nat, partB_alt)
    std = ds.moment_set(gB, partB, B, seed)
    D = pd.concat([std["draws"][CORR_TARGETS].reset_index(drop=True),
                   draws_frame({k: run["specs"]["primary"]["stats"][k] for k in CORR_NEW
                                if k in run["specs"]["primary"]["stats"]})], axis=1)
    payload = dict(
        sealed=True,
        note=("Hold-out half B of the pre-registered locale split (district_settings.py), D5 children by family income "
              "(children_income.py). Sealed until the researchers open it: no log, report or summary reads or prints "
              "these numbers."),
        created=dt.datetime.now().isoformat(timespec="seconds"), holdout_seed=ho["seed"], B=B, bootstrap_seed=seed,
        n_cz=int(len(gB)), cz_ids=[int(c) for c in ho["B"]],
        national=dict(primary=nat[f"{PRIMARY_IT}:child"].summary()),
        moments={k: strip_draws(r["stats"]) for k, r in run["specs"].items()},
        b17024=strip_draws(run["b17024"]["stats"]) if run["b17024"] else None,
        corr_keys=CORR_TARGETS + [k for k in CORR_NEW if k in D.columns],
        corr=D.corr().round(6).to_numpy().tolist())
    path.write_text(json.dumps(ds.clean(payload), indent=1))
    LOG.info("wrote %s (contents not printed)", path.name)


def cmd_build(args: argparse.Namespace) -> int:
    _stack()
    B, seed = args.B, args.seed
    indir = HERE
    outdir = Path(args.outdir).resolve() if args.outdir else HERE / "estimates" / "settings"
    (outdir / "sealed").mkdir(parents=True, exist_ok=True)
    full_run = B == B_DEFAULT and seed == SEED
    if (B_DEFAULT, SEED) != (db.B_DEFAULT, db.SEED):
        raise SystemExit("B_DEFAULT / SEED differ from district_bootstrap.py")
    t0 = time.time()

    acs = {it: load_acsed(it) for it in ITERS}
    b17 = load_b17024() if B17024_PATH.exists() else None
    if b17 is None:
        LOG.warning("no B17024 file: the cross-check is skipped")
    nat = build_nationals(acs, b17)
    LOG.info("national distributions: %s", {k: v.n_records for k, v in nat.items()})

    inp = ds.build_audited_panel(indir)
    A = merge_ext(inp["A"], acs, b17)
    part, g = ds.scheme_tables(A, "locale")
    part = part[part.year == YEAR]
    # the `city_vs_rest` variant: town and rural districts join the suburb (partition rules unchanged), same CZs
    A_rest = A.copy()
    A_rest["locale_grp"] = A_rest["locale_grp"].replace({"town": "suburb", "rural": "suburb"})
    part_alt = de.partition(A_rest, "locale", ds.MIN_PUPILS, ds.MIN_GROUP_SHARE)
    part_alt = part_alt[(part_alt.year == YEAR) & part_alt.czone.isin(set(g.czone))]

    # --- reproduction: stop if the four locale targets do not come back -----------------------------
    std = ds.moment_set(g, part, B, seed)
    repro = reproduction_check(std, full_run)

    # --- hold-out split: the pre-registered one -------------------------------------------------------
    ho = ds.holdout_split(g)
    ho_check = check_holdout(ho, indir)
    gA, partA = ds.restrict(g, part, set(ho["A"]))
    gB, partB = ds.restrict(g, part, set(ho["B"]))
    partA_alt = part_alt[part_alt.czone.isin(set(ho["A"]))]
    partB_alt = part_alt[part_alt.czone.isin(set(ho["B"]))]

    # --- full sample and half A --------------------------------------------------------------------
    full = run_sample(g, part, B, seed, nat, part_alt)
    pov = poverty_check(part, g, full["boot"]) if "b_i0" in part.columns else None
    stdA = ds.moment_set(gA, partA, B, seed)
    runA = run_sample(gA, partA, B, seed, nat, partA_alt)
    LOG.info("full sample and half A done in %.0fs", time.time() - t0)
    for lab, r in (("full", full), ("half A", runA)):
        st = r["specs"]["primary"]["stats"]
        LOG.info("%s: city share by tercile %.3f / %.3f / %.3f, top - bottom %.4f (%.4f), n_cz %d", lab,
                 st["city_share_t1"]["value"], st["city_share_t2"]["value"], st["city_share_t3"]["value"],
                 st["city_share_gradient"]["value"], st["city_share_gradient"]["se"], st["city_share_gradient"]["n_cz"])

    # --- correlations with the four district targets --------------------------------------------------
    def corr_block(std_: dict[str, Any], run: dict[str, Any]) -> dict[str, Any]:
        P = run["specs"]["primary"]["stats"]
        keys = CORR_TARGETS + [k for k in CORR_NEW if k in P]
        D = pd.concat([std_["draws"][CORR_TARGETS].reset_index(drop=True),
                       draws_frame({k: P[k] for k in CORR_NEW if k in P})], axis=1)
        C = D.corr()
        return dict(keys=keys, matrix=C.loc[keys, keys].to_numpy(float).tolist(), draws=D,
                    targets={k: dict(value=std_["stats"][k]["value"], se=std_["stats"][k]["se"])
                             for k in CORR_TARGETS})

    corr_full, corr_A = corr_block(std, full), corr_block(stdA, runA)
    # the multiplicities are those of district_settings.py: the same draws for two of its statistics
    dchk = {}
    for k, col in (("gap_pupils", "gap_pupils"), ("gap_income", "gap_pm_ftyr_all")):
        mine = full["boot"].wmean(g[col].to_numpy(float))["draws"]
        dchk[k] = float(np.nanmax(np.abs(mine - std["draws"][k].to_numpy(float))))
    if max(dchk.values()) > 1e-9:
        raise SystemExit(f"bootstrap multiplicities differ from district_settings.py: {dchk}")
    LOG.info("bootstrap draws identical to district_settings.py (max abs difference %.1e)", max(dchk.values()))

    # --- sealed half B: written, not read back --------------------------------------------------------
    write_sealed_b(outdir / "sealed" / "children_income_holdout_B.json", gB, partB, B, seed, ho, nat, partB_alt)

    # --- checks and coverage ------------------------------------------------------------------------
    checks = dict(
        us_row=us_row_check(acs, nat), interpolation=interpolation_check(acs[PRIMARY_IT]),
        rounding=rounding_check(acs[PRIMARY_IT], part, PRIMARY_IT), coverage=coverage(part, g, nat, acs),
        reproduction=repro, holdout=ho_check, poverty=pov, draws_identical=dchk,
        kids_per_family=kids_per_family_check(acs))
    n_cz_ny = None
    ny = part.loc[part.leaid == NYC_LEAID, "czone"]
    if len(ny):
        n_cz_ny = float(g.pupils[g.czone == float(ny.iloc[0])].sum() / g.pupils.sum())
    checks["coverage"]["nyc_weight_share_of_cz"] = n_cz_ny
    checks["rest_partition"] = dict(
        n_cz=int(part_alt.czone.nunique()), n_districts=int(len(part_alt)), pupils=float(part_alt.pupils.sum()),
        pupils_city=float(part_alt.pupils[part_alt.loc_id == 1].sum()),
        pupils_other=float(part_alt.pupils[part_alt.loc_id == 2].sum()),
        districts_by_locale=A_rest.loc[A_rest.leaid.isin(part_alt.leaid), "locale_grp"].value_counts().to_dict(),
        pupils_town_rural=float(part_alt.pupils[part_alt.leaid.isin(A.leaid[A.locale_grp.isin(["town", "rural"])])].sum()))
    checks["acsed_records"] = dict(total=int(len(acs[PRIMARY_IT])),
                                   layers={k: int(v) for k, v in acs[PRIMARY_IT].layer.value_counts().items()})
    checks["saipe_poverty_gap"] = dict(value=std["stats"]["poverty_gap"]["value"], se=std["stats"]["poverty_gap"]["se"],
                                       n_cz=std["stats"]["poverty_gap"]["n_cz"])
    R = dict(
        meta=dict(generated=dt.datetime.now().isoformat(timespec="seconds"), script="data/spatial/children_income.py",
                  year=YEAR, B=B, seed=seed, holdout_seed=ds.HOLDOUT_SEED, n_cz=int(len(g)), full_run=full_run,
                  geo_unit=A.attrs.get("geo_unit"), primary_iteration=PRIMARY_IT, iterations=ITERS,
                  resampling="whole CZ with replacement; the multiplicity matrix of district_settings.py (n, B, seed)",
                  hybrid_share=inp["hybrid_share"]),
        source=dict(page=TV_PAGE, dataset="ACS-ED Children's Tabulation, ACS 5-year 2014-2018, Profile Tables",
                    tables=list(ACSED_TABLES), income_columns=INC_COLS, income_labels=INC_LABELS,
                    layout_file="ACS-ED_2014-2018_RecordLayouts.xlsx (children_income/ in raw/)",
                    geography=GEO_ALL, b17024_vars=B17_VARS),
        national={k: v.summary() for k, v in nat.items()},
        nlsy79_v6=NLSY_V6,
        locale=dict(n_cz=int(len(g)), specs=strip_draws({k: r["stats"] for k, r in full["specs"].items()}),
                    spec_labels={s.key: s.label for s in SPECS}, b17024=strip_draws(full["b17024"]["stats"])
                    if full["b17024"] else None,
                    spec_sizes={k: dict(n_districts=r["n_districts"], children=r["children"])
                                for k, r in full["specs"].items()},
                    corr_keys=corr_full["keys"], corr=corr_full["matrix"], targets=corr_full["targets"]),
        holdout_A=dict(n_cz=int(len(gA)), specs=strip_draws({k: r["stats"] for k, r in runA["specs"].items()}),
                       b17024=strip_draws(runA["b17024"]["stats"]) if runA["b17024"] else None,
                       corr_keys=corr_A["keys"], corr=corr_A["matrix"], targets=corr_A["targets"],
                       rule=ds.holdout_split.__doc__.split("\n\n")[0].replace("\n    ", " ").replace(
                           "Pre-registered split of the locale CZs.  ", "").strip(), seed=ho["seed"],
                       n_B=int(len(gB)), sealed_file="sealed/children_income_holdout_B.json"),
        checks=checks)
    R["caveats"] = build_caveats(R)

    clean = ds.clean(strip_draws(R))
    (outdir / "children_income.json").write_text(json.dumps(clean, indent=1))
    write_toml(outdir / "children_income_targets.toml", R)
    write_md(outdir / "children_income.md", R)
    corr_full["draws"].round(8).to_csv(outdir / "children_income_draws.csv", index_label="draw")
    corr_A["draws"].round(8).to_csv(outdir / "children_income_draws_holdout_A.csv", index_label="draw")
    LOG.info("wrote children_income.{md,json}, children_income_targets.toml, children_income_draws.csv, "
             "children_income_draws_holdout_A.csv and sealed/children_income_holdout_B.json to %s (%.0fs)", outdir,
             time.time() - t0)
    return 0


# --------------------------------------------------------------------------
# Caveats
# --------------------------------------------------------------------------

def build_caveats(R: dict[str, Any]) -> list[str]:
    L = R["locale"]["specs"]
    P = L["primary"]
    chk = R["checks"]
    cov = chk["coverage"]
    g = lambda k, key="city_share_gradient": L[k][key]["value"]
    pov = chk.get("poverty")
    out = [
        "**The income table counts families, not children.** The ACS-ED has no child-by-family-income cross-tabulation: "
        "CDP03.2 gives the families (households) in which the group's children live, by family income. Children are "
        "obtained by scaling every district's ten counts by its children per family (CDP05_1est / CDP03_34est), which "
        "reproduces the number of children in each district exactly and the distribution under the assumption that the "
        "number of children per family does not vary with family income inside a district (it does vary across "
        "districts, which the scaling captures). A family with children in two overlapping districts (an elementary "
        "and a secondary district) is counted in both records, a child in one; the unscaled variant (families as "
        f"published) gives {g('family'):.3f} against {g('primary'):.3f}.",
        "**Relevant children, not pupils.** Group 102 is every child who lives in the district and whose grade it "
        "serves, whatever the school type; it is the population the NLSY79 and the model's children correspond to. "
        "Group 105 (enrolled in public school) drops private-school and unenrolled children, who are richer and more "
        f"urban on average: the gradient is {g('public'):.3f} (SE {L['public']['city_share_gradient']['se']:.3f}). "
        "Group 102 is the primary for that reason; the public-school and enrolled groups are reported beside it. The "
        "CZ weights are the partitioned public pupils of the CCD, as for every district target, while the children "
        f"counted here are residents: {cov['matched_102']['by_location']['city']['children'] / 1e6:.1f} "
        f"million relevant children in the matched city districts against "
        f"{cov['matched_102']['by_location']['city']['pupils'] / 1e6:.1f} million pupils (charter and other "
        "independent LEAs are not in the pupil count; their children are counted in their district of residence).",
        "**Not like-for-like with the NLSY79 statistic.** V6 is the 1979 central-city share of SMSA children aged 14-16 "
        "with known central-city status (20-40% of SMSA residents have none), by tercile of 1978 family income of all "
        "14-16-year-olds, pooled across SMSAs. Here: 2014-18, city = the NCES locale of the child's district, "
        "terciles of the national distribution of children by family income in 2018 dollars, averaged over CZs. The "
        "levels differ for that reason (the city holds "
        f"{100 * P['city_share_children']['value']:.1f}% of the partition's children, 30.2% of the NLSY79's MSA children); "
        "the comparison of gradients is suggestive, not an estimate of a change over time. The children also differ in "
        "age: relevant children are of school age (about 5-17, and 18-19 without a diploma; pre-school children only "
        "where a district serves them), the NLSY79's are 14-16.",
        "**City and suburb are district-level.** A district is in one location, whatever its internal geography: large "
        "central-city districts include suburban neighborhoods and the reverse, so the city share by income is a "
        "district-composition statistic. Only districts of at least 500 pupils in CZs with both locations are used "
        f"(the effective number of CZs is {cov['n_eff_cz']:.0f}, the five largest carry "
        f"{100 * cov['top5_weight_share']:.0f}% of the weight, New York "
        f"{100 * (cov.get('nyc_weight_share_of_cz') or 0):.1f}% alone). Town and rural districts (NCES locale 31-43) are "
        "outside the partition, whereas the NLSY79's 'suburb' is everything in the SMSA outside the central city; the "
        "`city_vs_rest` variant puts them in the second location (against the first, the city) and gives "
        + (f"{g('city_vs_rest'):.3f} on {L['city_vs_rest']['city_share_gradient']['n_cz']} CZs "
           f"(the city then holds {100 * L['city_vs_rest']['city_share_children']['value']:.1f}% of the children)"
           if "city_vs_rest" in L else "no estimate") + ".",
        "**Matching.** "
        f"{cov['matched_102']['n_districts']:,} of the {cov['n_districts']:,} partition districts "
        f"({100 * cov['matched_102']['pupil_share']:.2f}% of pupils) have an ACS-ED record. The "
        f"{cov['unmatched']['n']} that do not ({cov['unmatched']['pupils'] / 1e3:.0f} thousand pupils) are mostly "
        "state-wide, county-office, regional vocational, virtual and charter LEAs with no residential boundary, plus a "
        "few the id match misses (e.g. North Colonie CSD, NY, 5,844 pupils: its ACS-ED record appears to carry the id "
        "3600971 against the CCD's 3621000; not crosswalked). Detroit's record is filed under the pre-2016 id 2612000; it is matched to "
        "2601103 (Detroit Public Schools Community District, 49,931 pupils) by an explicit crosswalk: without it the "
        f"gradient is {g('detroit_no_record'):.4f} (see section 4). New York City's record (filed under the "
        "Chancellor's id 3620580, which the audited panel carries for its collapsed row) is used, unlike the "
        f"earnings, poverty and SEDA targets, whose New York CZ has no city record; without it the gradient is "
        f"{g('nyc_no_record'):.3f}, the largest of the variants.",
        "**Interpolation.** The national cut points fall inside the $35-50k bracket (first) and at the edge of the "
        "$75-100k bracket (second, $100,063 for group 102), so the first tercile boundary assumes a uniform density "
        "inside $35-50k; the variant with terciles at the bracket edges closest to 1/3 and 2/3 ($50k and $100k, "
        f"national group shares {', '.join(f'{100 * s:.1f}' for s in R['national']['102:child']['bracket_group_shares'])}%) "
        f"gives {g('bracket_cuts'):.3f}. Medians interpolated from a district's own ten brackets exceed its published "
        f"median by {100 * chk['interpolation']['mean_log_error']:.1f}% on average (SD "
        f"{100 * chk['interpolation']['sd_log_error']:.1f}%), which shifts the level of the pooled-median gap little "
        "because both locations are interpolated.",
        "**Rounding.** ACS-ED counts are rounded to the nearest 5, and the values 1-3 never occur while 4 occurs in "
        f"{100 * chk['rounding']['share_cells_equal_4']:.1f}% of cells (so a count of 1-4 is evidently published as 4); "
        f"these cells hold {100 * chk['rounding']['mass_of_4_cells']:.2f}% of all families nationally and "
        f"{100 * chk['rounding']['partition_mass_of_4_cells']:.3f}% in the partition: negligible for pooled counts.",
        "**Standard errors.** The whole-CZ bootstrap resamples CZs, so it carries the sampling error of the ACS "
        "estimates (it enters the CZ-level statistics) and the dispersion of the statistic across CZs, but not the "
        "uncertainty of the national cut points (a fraction of a percent of income) or of the scaling assumption.",
    ]
    if pov:
        gp = pov["gaps"]
        out.append(
            "**Child poverty: two ACS sources agree, SAIPE shows a smaller gap.** On the districts that have all "
            f"three, the suburb-minus-city gap in the child poverty rate is {gp['acsed']['value']:.3f} (ACS-ED), "
            f"{gp['b17024']['value']:.3f} (B17024) and {gp['saipe']['value']:.3f} (SAIPE, the estimator behind the "
            f"scored `poverty_gap` = V3): the ACS gap is {gp['acsed']['value'] / gp['saipe']['value']:.2f} times the SAIPE "
            "gap. SAIPE district estimates are partly model-based (they may smooth differences between districts of one "
            "CZ; not tested here). If the direct ACS estimates are the better measure of the within-CZ gap, V3 "
            "understates the sorting by that factor; V3 is not changed here.")
    out.append(
        "**Not done.** (1) No child-weighted income distribution: the ACS-ED cross-tabulation does not exist (children "
        "are scaled, see above), and a PUMS version (children by exact family income, PUMAs classified as principal-city "
        "or suburban as in `teacher_residence.py`) was not attempted: PUMAs do not nest in school districts, so it "
        "would be a different city/suburb definition. (2) Income detail is limited to ten brackets (the "
        "Census detailed table B19101 has 16, but for all families, not families with children). (3) The national "
        "row of the EDGE preview is not defined for group 102 (for the nation 'relevant' children are all children), "
        "so the national cross-check of the district sum covers groups 104 and 105 only. (4) Census B17024 is a "
        "cross-check in another income concept (income relative to the poverty line, family size adjusted), for "
        "children 6-17 by residence (unified and elementary records, each child once), not a substitute. "
        "(5) Nothing here is a model-implied moment: the model's moment is the share of the young schooled in "
        "location 1 by tercile of the parent's pre-tax income across all adults (`parent_income_composition`), the "
        "data counterpart of which is `city_share_t1..3` with a single CZ's terciles replaced by national ones "
        "(the `within_cz_cuts` variant uses each CZ's own).")
    out.append(
        "**Seal.** The CZ ids and CZ weights of half B are public (`district_settings.json`), so its point estimates "
        "are in principle implied by the full-sample and half-A numbers; the seal only keeps them out of the logs and "
        "reports. Nothing in this script's log, this report, the JSON or the TOML reads or prints half B.")
    return out


# --------------------------------------------------------------------------
# Output: TOML
# --------------------------------------------------------------------------

TOML_LABELS = {
    "city_share_t1": "share of the CZ's children in national family-income tercile 1 (bottom) who live in the city",
    "city_share_t2": "same, tercile 2",
    "city_share_t3": "same, tercile 3 (top)",
    "city_share_gradient": "TARGET (within-CZ analogue of V6, NLSY79 -0.295 (0.037)): top minus bottom tercile, city share",
    "city_share_gradient_relative": "top minus bottom divided by the overall city share of children (reading aid)",
    "city_share_mid_minus_bottom": "middle minus bottom tercile, city share",
    "city_share_top_minus_mid": "top minus middle tercile, city share",
    "city_share_children": "share of the CZ's children (all terciles) who live in the city",
    "pooled_city_share_t1": "tercile 1 city share, pooled over CZs (ratio of sums over CZs; the NLSY79 statistic's structure)",
    "pooled_city_share_t2": "same, tercile 2",
    "pooled_city_share_t3": "same, tercile 3",
    "pooled_city_share_gradient": "pooled top minus bottom",
    "comp_city_t1": "share of the city's children in national tercile 1",
    "comp_city_t2": "same, tercile 2",
    "comp_city_t3": "same, tercile 3",
    "comp_suburb_t1": "share of the suburb's children in national tercile 1",
    "comp_suburb_t2": "same, tercile 2",
    "comp_suburb_t3": "same, tercile 3",
    "gap_median_family_income": "log(suburb / city) of the pooled median family income of children (ten brackets, linear interpolation)",
    "median_family_income_city": "pooled median family income of children, city, dollars (CZ-weighted mean)",
    "median_family_income_suburb": "same, suburb",
    "gap_mean_family_income": "log(suburb / city) of the mean family income of children (exact aggregate of the published district means)",
    "mean_family_income_city": "mean family income of children, city, dollars (CZ-weighted mean)",
    "mean_family_income_suburb": "same, suburb",
    "gap_lnmedian_published": "suburb minus city, children-weighted average of district log median family income (published medians)",
    "gap_child_poverty": "suburb minus city, poverty rate of related children 5-17 (ACS-ED CDP03_67pct, weighted by children), levels",
    "child_poverty_city": "child poverty rate, city",
    "child_poverty_suburb": "child poverty rate, suburb",
}


def write_toml(path: Path, R: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    m = R["meta"]
    P = R["locale"]["specs"]["primary"]
    nat = R["national"][f"{PRIMARY_IT}:child"]
    a("# children_income_targets.toml -- within-CZ parental-income composition of school-age children (task D5).")
    a(f"# Generated {m['generated']} by data/spatial/children_income.py; see children_income.md for definitions and caveats.")
    a("# NCES EDGE ACS-ED Children's Tabulation, ACS 2014-18, group 102 (relevant children), family income of the")
    a("# children's families in ten brackets scaled to children; national terciles of children by family income;")
    a("# locale partition (city = NCES locale 11-13, suburb = 21-23), the 188 CZs of the calibration.")
    a("# value / se / ci95 / n_cz: pupil-weighted mean over CZs of the within-CZ statistic (location 1 = city, 2 = suburb),")
    a(f"# whole-CZ bootstrap SE (B = {m['B']}, seed {m['seed']}; the multiplicity matrix of district_settings.py, so the")
    a("# correlations with the four district targets are valid). Full-sample values only: half A is in the report and")
    a("# the JSON, half B is sealed.")
    a("")
    a("[meta]")
    for k in ("generated", "script", "year", "B", "seed", "n_cz"):
        a(f"{k} = {ds.toml_val(m[k])}")
    a(f"primary_iteration = {ds.toml_val(m['primary_iteration'])}")
    a("")
    a("[locale]")
    a('description = "locale split (city vs suburb), the 188 CZs of the tracked moments; relevant children (ACS-ED group 102)"')
    a(f"n_cz = {ds.toml_val(R['locale']['n_cz'])}")
    a(f"corr_keys = {ds.toml_val(R['locale']['corr_keys'])}")
    a(f"corr = {ds.toml_val(R['locale']['corr'])}")
    a("")
    for k, lab in TOML_LABELS.items():
        s = P.get(k)
        if not s:
            continue
        a(f"[locale.{k}]")
        a(f"label = {ds.toml_val(lab)}")
        a(f"value = {ds.toml_val(s['value'])}")
        a(f"se = {ds.toml_val(s['se'])}")
        a(f"ci95 = {ds.toml_val([s['ci95_lo'], s['ci95_hi']])}")
        a(f"n_cz = {ds.toml_val(s['n_cz'])}")
        a("")
    a("[locale.brackets]")
    a('description = "ten family-income brackets (dollars of 2018): shares of each location\'s children, and the city share of the children in each bracket"')
    a(f"labels = {ds.toml_val(INC_LABELS)}")
    for nm, key in (("city_share_of_children_in_bracket", "city_share"), ("location_share_city", "city"),
                    ("location_share_suburb", "suburb"), ("location_share_diff_suburb_minus_city", "diff")):
        v = [P[f"bin{k}_{key}"]["value"] for k in range(NB)]
        se = [P[f"bin{k}_{key}"]["se"] for k in range(NB)]
        a(f"{nm} = {ds.toml_val(v)}")
        a(f"{nm}_se = {ds.toml_val(se)}")
    a("")
    a("[national]")
    a('description = "national distribution of relevant children by family income: sum of the ACS-ED district records (child-scaled counts)"')
    a(f"n_records = {ds.toml_val(nat['n_records'])}")
    a(f"children = {ds.toml_val(nat['total'])}")
    a(f"bin_labels = {ds.toml_val(INC_LABELS)}")
    a(f"bin_shares = {ds.toml_val(nat['bin_shares'])}")
    a(f"tercile_cut_points_dollars = {ds.toml_val(nat['cut_points'])}")
    a(f"bracket_cut_edges_dollars = {ds.toml_val(nat['bracket_cut_edges'])}")
    a(f"bracket_group_shares = {ds.toml_val(nat['bracket_group_shares'])}")
    a("")
    a("[nlsy79_v6]")
    a('description = "reference: NLSY79 central-city share of SMSA children 14-16 by tercile of parental income (nlsy_extensions.md T3; SE of the gradient from nlsy_composition.md)"')
    a(f"city_share = {ds.toml_val(NLSY_V6['city_share'])}")
    a(f"gradient = {ds.toml_val(NLSY_V6['gradient'])}")
    a(f"gradient_se = {ds.toml_val(NLSY_V6['gradient_se'])}")
    a("")
    for key, st in R["locale"]["specs"].items():
        if key == "primary":
            continue
        a(f"[variants.{key}]")
        a(f"label = {ds.toml_val(R['locale']['spec_labels'][key])}")
        for t in (1, 2, 3):
            a(f"city_share_t{t} = {ds.toml_val(st[f'city_share_t{t}']['value'])}")
        a(f"city_share_gradient = {ds.toml_val(st['city_share_gradient']['value'])}")
        a(f"city_share_gradient_se = {ds.toml_val(st['city_share_gradient']['se'])}")
        a(f"gap_median_family_income = {ds.toml_val(st['gap_median_family_income']['value'])}")
        a(f"n_cz = {ds.toml_val(st['city_share_gradient']['n_cz'])}")
        a("")
    if R["locale"].get("b17024"):
        b = R["locale"]["b17024"]
        a("[variants.b17024]")
        a('label = "Census B17024: children 6-17 by income-to-poverty ratio, unified + elementary records, national IPR terciles (cut points 1.65 and 3.74)"')
        for t in (1, 2, 3):
            a(f"city_share_t{t} = {ds.toml_val(b[f'city_share_t{t}']['value'])}")
        a(f"city_share_gradient = {ds.toml_val(b['city_share_gradient']['value'])}")
        a(f"city_share_gradient_se = {ds.toml_val(b['city_share_gradient']['se'])}")
        a(f"n_cz = {ds.toml_val(b['city_share_gradient']['n_cz'])}")
        a("")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Output: report
# --------------------------------------------------------------------------

def write_md(path: Path, R: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    m, Lc, HA, chk = R["meta"], R["locale"], R["holdout_A"], R["checks"]
    S, SA = Lc["specs"], HA["specs"]
    P, PA = S["primary"], SA["primary"]
    nat = R["national"][f"{PRIMARY_IT}:child"]
    natf = R["national"][f"{PRIMARY_IT}:family"]
    cov = chk["coverage"]
    cuts = nat["cut_points"]
    gr = P["city_share_gradient"]
    nl = NLSY_V6
    z_nlsy = (gr["value"] - nl["gradient"]) / math.hypot(gr["se"], nl["gradient_se"])
    rel = P["city_share_gradient_relative"]
    nl_rel = nl["gradient"] / 0.302
    keys = Lc["corr_keys"]
    C = np.array(Lc["corr"], float)
    ci = {k: i for i, k in enumerate(keys)}

    def cor(new: str, tgt: str) -> float:
        return float(C[ci[new], ci[tgt]])

    def f3(s: dict[str, Any]) -> str:
        return fse(s["value"], s["se"], 3)

    def f4(s: dict[str, Any]) -> str:
        return fse(s["value"], s["se"], 4)

    def pp(s: dict[str, Any], d: int = 1) -> str:
        return pcs(s["value"], s["se"], d)

    def ival(s: dict[str, Any], d: int = 3) -> str:
        return f"[{s['ci95_lo']:.{d}f}, {s['ci95_hi']:.{d}f}]"

    others = {k: S[k]["city_share_gradient"]["value"] for k in S if k != "primary"}
    lo_k, hi_k = min(others, key=others.get), max(others, key=others.get)
    far = sorted(others, key=lambda k: -abs(others[k] - gr["value"]))[:2]
    b = Lc.get("b17024")
    pov = chk.get("poverty")
    sd = "district_settings.json"
    repro_str = ", ".join("%s %.4f (%.4f)" % (r["moment"], r["value"], r["se"]) for r in chk["reproduction"])

    a("# Within-CZ parental-income composition of school-age children (task D5, DESIGN_v2.md)")
    a("")
    a(f"Generated {m['generated']} by `children_income.py`. NCES EDGE ACS-ED Children's Tabulation, ACS 5-year 2014-18 "
      f"(Profile Tables CDP03.2, CDP05.1, CDP02.1, CDP03.4; group 102 'Relevant Children' is the primary, 104 and 105 are "
      f"variants), merged to the audited district panel of `district_settings.py` and cut by its locale partition (city = "
      f"NCES locale 11-13, suburb = 21-23; districts of at least 500 pupils; each location at least 10% of the CZ's "
      f"partitioned pupils; {cov['n_cz']} CZs, {cov['n_districts']:,} districts, {cov['pupils'] / 1e6:.2f} million pupils). "
      f"Statistics are within-CZ and are averaged over CZs with the partitioned pupils as weights; whole-CZ bootstrap, "
      f"B = {m['B']}, seed {m['seed']}, the multiplicity matrix of `district_settings.py` (so the bootstrap correlations "
      f"with the four district targets are valid). Hold-out half B ({HA['n_B']} CZs) is in "
      f"`{HA['sealed_file']}` and is not read, printed or summarised anywhere here; half A ({HA['n_cz']} CZs) is "
      f"reported in section 7. Dollars are 2018 dollars.")
    a("")
    a("## Summary")
    a("")
    a(f"- **Within-CZ analogue of V6.** The share of a CZ's relevant children who live in the city is "
      f"**{pp(P['city_share_t1'])}% in the bottom, {pp(P['city_share_t2'])}% in the middle and {pp(P['city_share_t3'])}% in "
      f"the top national family-income tercile** (cut points ${cuts[0]:,.0f} and ${cuts[1]:,.0f}). **Top minus bottom = "
      f"{gr['value']:.3f} (SE {gr['se']:.3f})**, 95% interval {ival(gr)}, on {gr['n_cz']} CZs. The NLSY79 V6 is 45.8 / 36.0 / "
      f"16.3% and {nl['gradient']:.3f} ({nl['gradient_se']:.3f}): the within-CZ gradient of 2014-18 is "
      f"{100 * gr['value'] / nl['gradient']:.0f}% as steep, {abs(z_nlsy):.1f} combined SEs flatter. Pooling the CZs (ratio "
      f"of sums, the NLSY79 statistic's structure) gives {pp(P['pooled_city_share_t1'])} / {pp(P['pooled_city_share_t2'])} / "
      f"{pp(P['pooled_city_share_t3'])}% and {f3(P['pooled_city_share_gradient'])}.")
    a(f"- **Levels differ, so the gradient should be read against the level.** The city holds "
      f"{pp(P['city_share_children'])}% of the partition's children (NLSY79: 30.2% of SMSA children with known status), "
      f"so the gradient relative to the overall city share is **{f3(rel)}** here against {nl_rel:.3f} for the NLSY79. "
      f"The model's location 1 holds 0.407 of the population (experiments-2026-10b LOG.md, I0), close to the "
      f"{100 * cov['matched_102']['by_location']['city']['pupils'] / cov['pupils']:.0f}% of pupils and "
      f"{100 * P['city_share_children']['value']:.0f}% of children here, not to the NLSY79's 30%.")
    a(f"- **The two locations.** Composition by national tercile (bottom / middle / top): city "
      f"{pp(P['comp_city_t1'])} / {pp(P['comp_city_t2'])} / {pp(P['comp_city_t3'])}%, suburb {pp(P['comp_suburb_t1'])} / "
      f"{pp(P['comp_suburb_t2'])} / {pp(P['comp_suburb_t3'])}%. Median family income of children (pooled brackets) "
      f"${P['median_family_income_city']['value']:,.0f} in the city, ${P['median_family_income_suburb']['value']:,.0f} in the "
      f"suburb: **log gap {f3(P['gap_median_family_income'])}** (mean family income {f3(P['gap_mean_family_income'])}; average "
      f"of district medians {f3(P['gap_lnmedian_published'])}). Child poverty (related children 5-17) "
      f"{pp(P['child_poverty_city'])}% against {pp(P['child_poverty_suburb'])}%, gap {f3(P['gap_child_poverty'])}.")
    a(f"- **Robustness.** Across the {len(others)} variants of section 4 the gradient runs from {others[lo_k]:.3f} "
      f"(`{lo_k}`) to {others[hi_k]:.3f} (`{hi_k}`); the largest moves from {gr['value']:.3f} are `{far[0]}` "
      f"({others[far[0]]:.3f}) and `{far[1]}` ({others[far[1]]:.3f}). Terciles at the "
      f"bracket edges (${nat['bracket_cut_edges'][0]:,.0f} and ${nat['bracket_cut_edges'][1]:,.0f}) give "
      f"{others['bracket_cuts']:.3f}; each CZ's own terciles {others['within_cz_cuts']:.3f}; families as published "
      f"{others['family']:.3f}."
      + (f" The Census B17024 cross-check (children 6-17 by income-to-poverty ratio, residence-based, a different income "
         f"concept) gives {f3(b['city_share_gradient'])}." if b else ""))
    a(f"- **Correlation with the four district targets** (bootstrap): the gradient with `gap_income` "
      f"{cor('city_share_gradient', 'gap_income'):+.2f}, `gap_teachers_pp` {cor('city_share_gradient', 'gap_teachers_pp'):+.2f}, "
      f"`gap_salary` {cor('city_share_gradient', 'gap_salary'):+.2f}, `gap_pupils` {cor('city_share_gradient', 'gap_pupils'):+.2f}. "
      f"The tercile shares themselves move with the size of the city (`city_share_t1` with `gap_pupils` "
      f"{cor('city_share_t1', 'gap_pupils'):+.2f}); the difference does not.")
    a(f"- **Half A** ({HA['n_cz']} CZs): {pp(PA['city_share_t1'])} / {pp(PA['city_share_t2'])} / {pp(PA['city_share_t3'])}%, "
      f"top minus bottom {f3(PA['city_share_gradient'])} (full sample {f3(gr)}); log median gap "
      f"{f3(PA['gap_median_family_income'])} (full {f3(P['gap_median_family_income'])}). Half B is sealed.")
    a(f"- **For the calibration.** `city_share_gradient` {gr['value']:.3f} ({gr['se']:.3f}) is the within-CZ counterpart of "
      f"the scored V6 ({nl['gradient']:.3f}); it is a different object, not a re-estimate of the same number (section 1, "
      f"caveats). For orientation the baseline model's V6 is +0.011 (calibration log, section 4), "
      f"{(0.011 - gr['value']) / gr['se']:.0f} SEs from it. The model's own counterpart is `parent_income_composition` "
      f"(terciles of parental pre-tax income over all adults of one CZ); the `within_cz_cuts` variant is the closest "
      f"data analogue of that single-CZ definition.")
    a("")
    a("## 1. Data, definitions and units")
    a("")
    a(f"- **Source.** NCES EDGE, ACS-ED Children's Tabulation, 'ACS 2014-2018 Profile' (school year 2017-18), Table Viewer "
      f"<{TV_PAGE}>: 'All Districts' export, one table and one population group per file ('Download Text File'; "
      f"pipe-delimited, one row per district record). {chk['acsed_records']['total']:,} district records: "
      f"{chk['acsed_records']['layers'].get('unified', 0):,} unified (GeoId prefix 97000US), "
      f"{chk['acsed_records']['layers'].get('elementary', 0):,} elementary (95000US), "
      f"{chk['acsed_records']['layers'].get('secondary', 0):,} secondary (96000US), keyed by the 7-digit NCES LEAID. "
      f"Cached under `raw/children_income/` with a manifest (`manifest_acsed.json`: URLs, bytes, SHA-256) together with "
      f"the ReadMe and the record-layout workbook (`ACS-ED_2014-2018_RecordLayouts.xlsx`) that ship in every export "
      f"(the workbook's variable descriptions say 'in 2017 inflation-adjusted dollars', the Table Viewer and the ACS "
      f"2014-18 documentation say 2018 dollars; terciles and log gaps do not depend on it).")
    a("- **Population groups** ('iterations'): 102 *Relevant Children* (every child 0-17, or 18-19 without a high school "
      "diploma, who lives in the district and whose assigned grade is within the district's grade range, each child "
      "relevant to one district only, so overlapping elementary and secondary districts do not double count children); "
      "104 *Relevant Children - Enrolled (public and private)*; 105 *Relevant Children - Enrolled Public*. **Primary: 102**, the "
      "population the NLSY79 statistic and the model's children correspond to (every child, whatever the school "
      "type); 105 is the 'children enrolled in public school' population of the task statement and is carried as a "
      "variant. The ACS-ED has no Children's detailed tables after the 2008-12 vintage, only the Profile Tables used here.")
    a("- **Variables** (record-layout names). CDP03.2 *Income and benefits*: `CDP03_34est` families; `CDP03_35est`-"
      "`CDP03_44est` families by **family income** in ten brackets (< $10k, $10-15k, $15-25k, $25-35k, $35-50k, $50-75k, "
      "$75-100k, $100-150k, $150-200k, $200k+); `CDP03_45est` / `CDP03_46est` median / mean family income. CDP05.1 *Sex and "
      "age*: `CDP05_1est` total population of the group = the number of children. CDP02.1: `CDP02_9est` average number "
      "of children in families with children (a check on the scaling). CDP03.4: `CDP03_67pct` related children 5-17 below the "
      "poverty level, percent. The household-income brackets (`CDP03_6`-`15`) are not used.")
    a("- **Unit and scaling.** The income table counts the *families in which the group's children live*; it is not a "
      "child-by-income tabulation. Each district's ten counts are scaled by its children per family (`CDP05_1est` / "
      f"`CDP03_34est`; national mean {chk['kids_per_family']['national_mean'][PRIMARY_IT]:.2f}; correlation with the "
      f"published `CDP02_9est` {chk['kids_per_family']['corr_ge500']:.3f} across the {chk['kids_per_family']['n_ge500']:,} "
      f"districts of at least 500 families, {chk['kids_per_family']['corr_all']:.2f} over all {chk['kids_per_family']['n_all']:,} "
      "where the rounding of tiny counts dominates), so that the district's "
      "counts add up to its children; this assumes constant family size across brackets inside a district. A family with "
      "children in two overlapping districts is in both records. The unscaled version is the `family` variant.")
    a(f"- **National terciles.** The national distribution of relevant children by family income is the sum of the "
      f"{nat['n_records']:,} district records (all valid): {nat['total'] / 1e6:.2f} million children. The cut points at "
      f"one third and two thirds are interpolated linearly inside the bracket that holds them (uniform density inside "
      f"a bracket, the same split in every district): **${cuts[0]:,.0f}** (inside $35-50k) and **${cuts[1]:,.0f}** "
      f"(inside $75-100k). Families as published: ${natf['cut_points'][0]:,.0f} and ${natf['cut_points'][1]:,.0f}. "
      f"The terciles are of the nation, as for V6 (NLSY79 cut points are over all 14-16-year-olds with income).")
    a("- **Statistics.** For CZ c, location l (1 = city, 2 = suburb) and tercile t, `T[c,l,t]` is the sum over the "
      "districts of the location of each bracket's count times the share of the bracket in t. "
      "`city_share_t` = T[c,1,t] / (T[c,1,t] + T[c,2,t]); `city_share_gradient` = city_share_3 - city_share_1 (the "
      "within-CZ analogue of V6); both averaged over CZs with the CZ weights (partitioned pupils). `pooled_*` is the "
      "ratio of sums over CZs. `comp_*` is the tercile composition of a location. Medians come from the pooled "
      "bracket counts of a location by linear interpolation (NaN if the median is in the open $200k+ bracket); the mean "
      "gap is the exact aggregate of the published district means (weights = children); child poverty is weighted by children.")
    a(f"- **Panel and sample.** `district_settings.build_audited_panel` (New York City one district; its collapsed row "
      f"carries the Chancellor's id 3620580, under which the ACS-ED files it), locale partition, base year 2018: "
      f"{cov['n_districts']:,} districts. {cov['matched_102']['n_districts']:,} have an ACS-ED record "
      f"({100 * cov['matched_102']['pupil_share']:.2f}% of pupils; city {100 * cov['matched_102']['by_location']['city']['pupil_share_matched']:.2f}%, "
      f"suburb {100 * cov['matched_102']['by_location']['suburb']['pupil_share_matched']:.2f}%): layers "
      f"{cov['matched_102']['layers']}. The four district targets reproduce ({repro_str}); the bootstrap "
      f"draws of `gap_pupils` and `gap_income` rebuilt here equal those of `district_settings.py` (max abs difference "
      f"{max(chk['draws_identical'].values()):.0e}); the hold-out split is identical to the one stored in `{sd}` (seed {HA['seed']}).")
    a("")
    a("## 2. City share by national family-income tercile")
    a("")
    a("City share of the CZ's relevant children (group 102, children by scaling), CZ-weighted mean of the within-CZ "
      "ratios, whole-CZ bootstrap SE in parentheses; 'pooled' is the ratio of sums over CZs. NLSY79: central-city share "
      "of SMSA children 14-16 by tercile of 1978 family income (cut points $14,000 and $23,000, nominal).")
    a("")
    a("| tercile | family income cut | share of the nation's children | city share (SE) | 95% interval | pooled over CZs (SE) | NLSY79 V6 |")
    a("|---|---|---|---|---|---|---|")
    labs = [f"below ${cuts[0]:,.0f}", f"${cuts[0]:,.0f} to ${cuts[1]:,.0f}", f"${cuts[1]:,.0f} and above"]
    for t in range(3):
        s, pl = P[f"city_share_t{t + 1}"], P[f"pooled_city_share_t{t + 1}"]
        a(f"| T{t + 1}{' (bottom)' if t == 0 else ' (top)' if t == 2 else ''} | {labs[t]} | 33.3% | {pp(s)}% | "
          f"[{100 * s['ci95_lo']:.1f}, {100 * s['ci95_hi']:.1f}] | {pp(pl)}% | {100 * nl['city_share'][t]:.1f}% |")
    a(f"| **top minus bottom** | | | **{f3(gr)}** | {ival(gr)} | {f3(P['pooled_city_share_gradient'])} | {nl['gradient']:.3f} ({nl['gradient_se']:.3f}) |")
    a(f"| middle minus bottom | | | {f3(P['city_share_mid_minus_bottom'])} | {ival(P['city_share_mid_minus_bottom'])} | | |")
    a(f"| top minus middle | | | {f3(P['city_share_top_minus_mid'])} | {ival(P['city_share_top_minus_mid'])} | | |")
    a(f"| all children | | | {pp(P['city_share_children'])}% | [{100 * P['city_share_children']['ci95_lo']:.1f}, {100 * P['city_share_children']['ci95_hi']:.1f}] | | 30.2% (1.5) |")
    a(f"| gradient / overall city share | | | {f3(rel)} | {ival(rel)} | | {nl_rel:.3f} |")
    a("")
    a(f"The within-CZ gradient is {gr['value']:.3f}; the NLSY79 V6 is {nl['gradient']:.3f} ({nl['gradient_se']:.3f}); the "
      f"difference is {gr['value'] - nl['gradient']:+.3f} ({z_nlsy:+.1f} combined SEs; the two samples are independent). "
      f"n = {gr['n_cz']} of {m['n_cz']} CZs (a CZ drops out only if a location has no district with an ACS-ED record).")
    a("")
    a("### Composition of the two locations by national tercile")
    a("")
    a("| tercile | city: share of its children (SE) | suburb (SE) | suburb minus city, pp |")
    a("|---|---|---|---|")
    for t in range(3):
        c, s_ = P[f"comp_city_t{t + 1}"], P[f"comp_suburb_t{t + 1}"]
        a(f"| T{t + 1} | {pp(c)}% | {pp(s_)}% | {100 * (s_['value'] - c['value']):+.1f} |")
    a("")
    a("## 3. Children by family income in the two locations")
    a("")
    a("Shares of each location's children in the ten brackets (CZ-weighted means of the within-CZ shares), suburb minus "
      "city, and the city share of the children in each bracket (CZ-weighted mean of the within-CZ ratios).")
    a("")
    a("| family income | city, % (SE) | suburb, % (SE) | suburb minus city, pp (SE) | city share of the bracket's children, % (SE) | nation, % |")
    a("|---|---|---|---|---|---|")
    for k in range(NB):
        a(f"| {INC_LABELS[k]} | {pp(P[f'bin{k}_city'])} | {pp(P[f'bin{k}_suburb'])} | "
          f"{pcs(P[f'bin{k}_diff']['value'], P[f'bin{k}_diff']['se'], 1)} | {pp(P[f'bin{k}_city_share'])} | "
          f"{100 * nat['bin_shares'][k]:.1f} |")
    a("")
    a("| statistic | city | suburb | suburb minus city (SE) | 95% interval |")
    a("|---|---|---|---|---|")
    a(f"| pooled median family income of children, $ | {P['median_family_income_city']['value']:,.0f} | "
      f"{P['median_family_income_suburb']['value']:,.0f} | log gap **{f4(P['gap_median_family_income'])}** | {ival(P['gap_median_family_income'])} |")
    a(f"| mean family income of children, $ | {P['mean_family_income_city']['value']:,.0f} | {P['mean_family_income_suburb']['value']:,.0f} | "
      f"log gap {f4(P['gap_mean_family_income'])} | {ival(P['gap_mean_family_income'])} |")
    a(f"| average district log median family income (published medians) | | | {f4(P['gap_lnmedian_published'])} | {ival(P['gap_lnmedian_published'])} |")
    a(f"| child poverty rate, related children 5-17 | {pp(P['child_poverty_city'])}% | {pp(P['child_poverty_suburb'])}% | "
      f"{f4(P['gap_child_poverty'])} (levels) | {ival(P['gap_child_poverty'])} |")
    a("")
    a(f"For scale: the audited earnings gap `gap_income` (log pooled median earnings of full-time year-round workers, "
      f"suburb minus city) is {R['locale']['targets']['gap_income']['value']:.4f} "
      f"({R['locale']['targets']['gap_income']['se']:.4f}); the family income of *children* differs by "
      f"{P['gap_median_family_income']['value']:.2f} in the median and {P['gap_mean_family_income']['value']:.2f} in the mean "
      f"because family income aggregates two earners, transfers and capital income, and children are concentrated in "
      f"high-income suburban families.")
    a("")
    a("## 4. Sensitivity")
    a("")
    a("Full sample and half A, same national terciles unless stated. T1-T3: city share by tercile; the last columns are the "
      "log suburb-minus-city gap in median family income and the child poverty gap.")
    a("")
    a("| variant | what | T1, % | T2, % | T3, % | top minus bottom (SE) | half A | median gap (SE) | poverty gap | n CZs |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for key, st in S.items():
        sa = SA[key]["city_share_gradient"]
        a(f"| `{key}` | {Lc['spec_labels'][key]} | {pc(st['city_share_t1']['value'])} | {pc(st['city_share_t2']['value'])} | "
          f"{pc(st['city_share_t3']['value'])} | {f3(st['city_share_gradient'])} | {sa['value']:.3f} | "
          f"{f3(st['gap_median_family_income'])} | {st['gap_child_poverty']['value']:.3f} | {st['city_share_gradient']['n_cz']} |")
    if b:
        ba = HA["b17024"]
        a(f"| `b17024` | Census B17024, children 6-17 by income-to-poverty ratio, unified + elementary records, national IPR terciles "
          f"(cut points {R['national']['b17024']['cut_points'][0]:.2f} and {R['national']['b17024']['cut_points'][1]:.2f}) | "
          f"{pc(b['city_share_t1']['value'])} | {pc(b['city_share_t2']['value'])} | {pc(b['city_share_t3']['value'])} | "
          f"{f3(b['city_share_gradient'])} | {ba['city_share_gradient']['value']:.3f} | {f3(b['gap_median_ipr'])} (median IPR) | "
          f"{b['gap_poverty_ipr']['value']:.3f} | {b['city_share_gradient']['n_cz']} |")
    a("")
    a("The `equal_cz_weights` variant weights every CZ equally (levels rise because the city share of children is larger in "
      "small CZs); `within_cz_cuts` uses each CZ's own terciles, so T1-T3 are within-CZ terciles and the gradient is free of "
      "between-CZ income differences. The gradient steepens with CZ size: "
      f"{S['equal_cz_weights']['city_share_gradient']['value']:.3f} with equal weights, {gr['value']:.3f} with pupil "
      f"weights, {P['pooled_city_share_gradient']['value']:.3f} pooled over CZs (the pooled ratio is dominated by the "
      "largest CZs).")
    a("")
    if pov:
        gp = pov["gaps"]
        a("## 5. Child poverty: three sources on the same districts")
        a("")
        a(f"Suburb minus city gap in the poverty rate of children, districts with an ACS-ED record, a B17024 record and a "
          f"SAIPE record ({pov['n_districts_matched']:,} districts, {100 * pov['pupil_share_matched']:.1f}% of pupils; "
          f"New York City, whose SAIPE record is filed under another id, is not in it), CZ-weighted, same draws.")
        a("")
        a("| source | city | suburb | gap (SE) |")
        a("|---|---|---|---|")
        for k in ("acsed", "b17024", "saipe"):
            s = gp[k]
            a(f"| {s['label']} | {100 * s['city']:.1f}% | {100 * s['suburb']:.1f}% | {f4(s)} |")
        a("")
        a(f"The two ACS sources agree; SAIPE, the estimator behind the scored `poverty_gap` (V3), shows a gap "
          f"{100 * (1 - gp['saipe']['value'] / gp['acsed']['value']):.0f}% smaller on the same districts (the ACS gap is "
          f"{gp['acsed']['value'] / gp['saipe']['value']:.2f} times the SAIPE gap). SAIPE district estimates are partly "
          f"model-based; whether that explains the difference was not tested here. On the full ACS-ED sample the gap is "
          f"{f4(P['gap_child_poverty'])} (with New York City's record) and the audited SAIPE-based `poverty_gap` is "
          f"{chk['saipe_poverty_gap']['value']:.4f} ({chk['saipe_poverty_gap']['se']:.4f}).")
        a("")
    a("## 6. Bootstrap correlation with the four district targets")
    a("")
    a("Correlation of the whole-CZ bootstrap draws of the new statistics with the draws of the four audited locale "
      "targets (the same multiplicity matrix, B = 1999, seed 20260929).")
    a("")
    a("| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |")
    a("|---|---|---|---|---|")
    for k in [k for k in keys if k not in CORR_TARGETS]:
        a(f"| `{k}` | " + " | ".join(f"{cor(k, t):+.2f}" for t in CORR_TARGETS) + " |")
    a("")
    a("Targets (value, SE): " + "; ".join(f"`{k}` {v['value']:.4f} ({v['se']:.4f})" for k, v in Lc["targets"].items()) + ".")
    a("")
    a("## 7. Hold-out half A")
    a("")
    a(f"Pre-registered split of the locale CZs (`district_settings.holdout_split`: {HA['rule']}; seed {HA['seed']}); half A "
      f"has {HA['n_cz']} CZs, half B {HA['n_B']}. Bootstrap on half A's own multiplicity matrix (n = {HA['n_cz']}, B = {m['B']}, "
      f"seed {m['seed']}). Half B was computed by the same code and written only to `{HA['sealed_file']}`.")
    a("")
    a("| statistic | full sample (SE) | half A (SE) |")
    a("|---|---|---|")
    for k, lab in (("city_share_t1", "city share, tercile 1, %"), ("city_share_t2", "tercile 2, %"),
                   ("city_share_t3", "tercile 3, %"), ("city_share_children", "all children, %")):
        a(f"| {lab} | {pp(P[k])} | {pp(PA[k])} |")
    for k, lab in (("city_share_gradient", "**top minus bottom**"), ("pooled_city_share_gradient", "top minus bottom, pooled over CZs"),
                   ("city_share_gradient_relative", "gradient / overall city share"),
                   ("gap_median_family_income", "log gap, median family income of children"),
                   ("gap_mean_family_income", "log gap, mean family income"), ("gap_child_poverty", "child poverty gap")):
        a(f"| {lab} | {f3(P[k])} | {f3(PA[k])} |")
    a("")
    tA = HA["targets"]
    a("Half-A district targets (value, SE): " + "; ".join(f"`{k}` {v['value']:.4f} ({v['se']:.4f})" for k, v in tA.items()) + ".")
    a("")
    a("## 8. Checks")
    a("")
    us = chk["us_row"]
    for it in ("104", "105"):
        if us.get(it, {}).get("available"):
            u = us[it]
            a(f"- **National row, group {it}.** The EDGE preview's national table has {u['us_families']:,.0f} families; the sum of the "
              f"district records has {u['districts_families']:,.0f} ({100 * (u['ratio_districts_to_us'] - 1):+.2f}%: families with "
              f"children in two overlapping districts are in both records); bracket shares differ by at most "
              f"{100 * u['max_abs_share_diff']:.2f} pp; national median / mean family income ${u['us_median']:,.0f} / ${u['us_mean']:,.0f}.")
    a("- **National row, group 102.** Not defined by the EDGE viewer: for the nation 'relevant' children are all children.")
    ic = chk["interpolation"]
    a(f"- **Interpolated against published medians** ({ic['n']:,} districts): mean log error {ic['mean_log_error']:+.4f} "
      f"(children-weighted {ic['weighted_mean_log_error']:+.4f}), SD {ic['sd_log_error']:.3f}, 1st / 99th percentile "
      f"{ic['p01']:+.3f} / {ic['p99']:+.3f}; {ic['n_published_median_missing']} districts have no published median "
      f"(open top interval or suppressed).")
    rc = chk["rounding"]
    a(f"- **Rounding.** {100 * rc['share_cells_multiple_of_5']:.1f}% of bracket counts are multiples of 5, none equals 1-3, "
      f"{100 * rc['share_cells_equal_4']:.1f}% equal 4 (a count of 1-4 is evidently published as 4: {100 * rc['mass_of_4_cells']:.2f}% of national "
      f"family mass; {100 * rc['partition_mass_of_4_cells']:.3f}% in the partition).")
    a(f"- **Unmatched partition districts.** {cov['unmatched']['n']} ({cov['unmatched']['pupils']:,.0f} pupils, "
      f"{100 * cov['unmatched']['pupils'] / cov['pupils']:.2f}%); largest: "
      + "; ".join(f"{u['name']} ({u['pupils']:,.0f})" for u in cov['unmatched']['largest'][:4]) + ".")
    a(f"- **Children against pupils.** Matched city districts: {cov['matched_102']['by_location']['city']['children'] / 1e6:.2f} million "
      f"relevant children, {cov['matched_102']['by_location']['city']['pupils'] / 1e6:.2f} million CCD pupils; suburb "
      f"{cov['matched_102']['by_location']['suburb']['children'] / 1e6:.2f} / {cov['matched_102']['by_location']['suburb']['pupils'] / 1e6:.2f} million.")
    a("")
    a("## 9. Caveats and what could not be done")
    a("")
    for c in R["caveats"]:
        a(f"- {c}")
    a("")
    a("## Files")
    a("")
    a("- `children_income.py` (fetch / build / selftest); `children_income.json` (all statistics, full sample and half A); "
      "`children_income_targets.toml` (full-sample values and SEs); `children_income_draws.csv` and "
      "`children_income_draws_holdout_A.csv` (bootstrap draws with the four district targets); "
      f"`{HA['sealed_file']}` (half B, sealed).")
    a("- Raw: `raw/children_income/` (ACS-ED exports for groups 102, 104, 105, the ReadMe and record layouts, "
      "`manifest_acsed.json`, B17024 parquet).")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Self-test and entry point
# --------------------------------------------------------------------------

def selftest() -> None:
    """Unit checks of the helpers on synthetic data (no data files, no network)."""
    _stack()
    e = np.array([0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, np.inf])
    uni = np.ones(10)
    c = cuts_from_counts(uni, e)
    assert np.allclose(c, [10 / 3, 20 / 3]), c
    W = tercile_weights(c, e)
    assert np.allclose(W.sum(axis=1), 1.0)
    assert np.allclose(uni @ W / uni.sum(), 1 / 3)                       # each tercile holds a third
    skew = np.array([5.0, 4, 3, 2, 1, 1, 1, 1, 1, 1])
    cs = cuts_from_counts(skew, e)
    Ws = tercile_weights(cs, e)
    assert np.allclose(skew @ Ws / skew.sum(), 1 / 3), skew @ Ws / skew.sum()
    Wb, info = bracket_weights(skew, e)
    assert set(np.unique(Wb)) <= {0.0, 1.0} and np.allclose(Wb.sum(axis=1), 1.0)
    assert abs(sum(info["group_shares"]) - 1.0) < 1e-12
    assert grouped_median(np.array([0.0, 10, 0, 0]), np.array([0.0, 1, 2, 3, np.inf])) == 1.5
    assert np.isnan(grouped_median(np.array([0.0, 0, 0, 5]), np.array([0.0, 1, 2, 3, np.inf])))
    assert np.isnan(grouped_median(np.zeros(4), np.array([0.0, 1, 2, 3, np.inf])))
    # an open top bracket goes to the tercile that holds its lower edge
    assert tercile_weights(np.array([2.5, 6.0]), e)[9, 2] == 1.0
    # city holds only poor children, suburb only rich: shares 1 / ? / 0, gradient -1
    n = 4
    cnt = np.zeros((n, 2, 10))
    cnt[:, 0, 0], cnt[:, 1, 9] = 100.0, 100.0
    boot = Boot(n, 50, 1, np.array([1.0, 2.0, 3.0, 4.0]))
    st = core_stats(cnt, tercile_weights(c, e), boot)
    assert abs(st["city_share_t1"]["value"] - 1.0) < 1e-12 and abs(st["city_share_t3"]["value"]) < 1e-12
    assert abs(st["city_share_gradient"]["value"] + 1.0) < 1e-12
    assert abs(st["city_share_children"]["value"] - 0.5) < 1e-12
    assert abs(st["comp_city_t1"]["value"] - 1.0) < 1e-12 and abs(st["comp_suburb_t3"]["value"] - 1.0) < 1e-12
    # weighted mean equals the weighted mean; bootstrap SE positive for a varying statistic
    x = np.array([0.1, 0.2, 0.3, 0.4])
    s = boot.wmean(x)
    assert abs(s["value"] - np.average(x, weights=boot.w)) < 1e-12 and s["se"] > 0
    # district scaling: counts add up to the district's children
    P = pd.DataFrame({"a102_fam": [100.0, 50.0], "a102_kids": [180.0, 60.0],
                      **{f"a102_f{k}": [10.0, 5.0] for k in range(NB)}, "a102_mean_fam": [1.0, 2.0],
                      "a102_med_fam": [1.0, 2.0], "a102_pov517": [10.0, 20.0]})
    D = district_counts(P, "102", "child")
    assert np.allclose(D[[f"n{k}" for k in range(NB)]].sum(axis=1), [180.0, 60.0])
    assert np.allclose(district_counts(P, "102", "family")["wt"], [100.0, 50.0])
    # national row parser on a synthetic preview
    html = ('<table><tr><td class="label, indent1">Families</td><td>1,000</td></tr>'
            + "".join(f'<tr><td class="label, indent2">b{k}</td><td>{100 * (k + 1):,}</td></tr>' for k in range(NB))
            + '<tr><td class="label, indent2">Median family income (dollars)</td><td>70,000</td></tr>'
            '<tr><td class="label, indent2">Mean family income (dollars)</td><td>90,000</td></tr></table>')
    tmp = RAW / "_selftest_us.json"
    RAW.mkdir(parents=True, exist_ok=True)
    tmp.write_text(json.dumps({"html": html}))
    try:
        r = parse_us_families(tmp)
    finally:
        tmp.unlink()
    assert r["total"] == 1000 and len(r["brackets"]) == NB and r["median"] == 70000 and r["mean"] == 90000
    print("selftest passed")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    f = sub.add_parser("fetch", help="download and cache the ACS-ED exports and the Census B17024 table")
    f.add_argument("--only", choices=["acsed", "b17024"])
    f.add_argument("--refresh", action="store_true", help="download again even if cached")
    f.add_argument("--workers", type=int, default=8)
    sub.add_parser("selftest", help="unit checks of the helpers")
    bd = sub.add_parser("build", help="compute the statistics and write the reports")
    bd.add_argument("--B", type=int, default=B_DEFAULT)
    bd.add_argument("--seed", type=int, default=SEED)
    bd.add_argument("--outdir", default=None, help="default: <here>/estimates/settings")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    if args.cmd == "fetch":
        cmd_fetch(args)
    elif args.cmd == "selftest":
        selftest()
    elif args.cmd == "build":
        return cmd_build(args)
    else:
        ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
