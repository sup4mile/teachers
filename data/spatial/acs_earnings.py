#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acs_earnings.py -- ACS 5-year earnings tables by school district, for the
income statistic of the spatial calibration (item T3).

The model's income moment is the suburb-minus-city gap in LABOR income of
25-34-year-olds.  `fetch_data.py --sources acs` only pulls median household
income and per-capita income, which mix transfers, capital income and household
composition.  This script adds the earnings tables, for the three school-district
layers (unified, elementary, secondary), one state per request as `fetch_data.py`
does (the district layers refuse `in=state:*` before the 2022 vintage).

Tables (ACS 5-year, default 2018 = 2014-2018, dollars of 2018):
    B20002  median earnings, population 16+ with earnings: total / male / female
    B20017  median earnings by sex x work experience (full-time year-round)
    B20018  median earnings, full-time year-round workers 16+
    B20004  median earnings, population 25+, by sex and by education (total)
    B24011  median earnings, civilian employed 16+
    B19049  median household income by age of householder (25-44 is _003)
    B19013  median household income (all households)
    B19061 / B19051  aggregate earnings / households with earnings
    B20001  earners 16+ (count; a population weight)
    B01003  total population

The full earnings distributions (bins) come from two whole-group pulls,
`group(B20001)` (sex x earnings bins, population 16+ with earnings) and
`group(B20005)` (sex x work experience x earnings bins), so that a LOCATION's
median can be computed from pooled bin counts rather than from an average of
district medians.

Estimates come with the matching margins of error for the medians.  Census
annotation codes (-666666666 no sample, -222222222 median in the open-ended top
interval, ...) are kept exactly as published; they are negative and must be
masked before any estimate.

Output (under raw/, which is gitignored):
    raw/acs_earnings_school_districts_<year>.parquet    one row per district
        with `leaid` (state FIPS + 5-digit district code), `sd_layer`, `year`.
    raw/acs_earnings_bins_school_districts_<year>.parquet   B20001 and B20005
        estimate columns (counts by earnings bin), same keys.

Usage
    set -a; source .env; set +a          # CENSUS_API_KEY
    python3 acs_earnings.py              # vintage 2018
    python3 acs_earnings.py --year 2018 --workers 8

The key is read from the environment and is never printed or written: error
strings are scrubbed of it before they reach the log.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

LOG = logging.getLogger("acs_earnings")
HERE = Path(__file__).resolve().parent

LAYERS = ("school district (unified)", "school district (elementary)",
          "school district (secondary)")

MEDIANS = {                      # estimate + margin of error requested
    "B20002_001": "median earnings, 16+ with earnings, total",
    "B20002_002": "median earnings, 16+ with earnings, male",
    "B20002_003": "median earnings, 16+ with earnings, female",
    "B20017_001": "median earnings, 16+ with earnings, total (by work experience)",
    "B20017_002": "median earnings, male, total",
    "B20017_003": "median earnings, male, full-time year-round",
    "B20017_005": "median earnings, female, total",
    "B20017_006": "median earnings, female, full-time year-round",
    "B20018_001": "median earnings, full-time year-round 16+",
    "B20004_001": "median earnings, 25+ with earnings, total",
    "B24011_001": "median earnings, civilian employed 16+",
    "B19049_001": "median household income, all householders",
    "B19049_003": "median household income, householder 25-44",
    "B19013_001": "median household income",
}
COUNTS_AND_OTHER = {             # estimate only
    "B20004_002": "median earnings 25+, less than high school",
    "B20004_003": "median earnings 25+, high school graduate",
    "B20004_004": "median earnings 25+, some college / associate",
    "B20004_005": "median earnings 25+, bachelor's",
    "B20004_006": "median earnings 25+, graduate / professional",
    "B20004_007": "median earnings 25+, male total",
    "B20004_013": "median earnings 25+, female total",
    "B19049_002": "median household income, householder under 25",
    "B19049_004": "median household income, householder 45-64",
    "B19049_005": "median household income, householder 65+",
    "B19061_001": "aggregate household earnings",
    "B19051_001": "households",
    "B19051_002": "households with earnings",
    "B20001_001": "earners, population 16+ with earnings",
    "B01003_001": "total population",
}


def get_vars() -> list[str]:
    v = [f"{k}E" for k in MEDIANS] + [f"{k}M" for k in MEDIANS]
    v += [f"{k}E" for k in COUNTS_AND_OTHER]
    return v


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=5, connect=5, read=5, status=5, backoff_factor=1.0,
                  status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset(["GET"]), raise_on_status=False)
    ad = HTTPAdapter(max_retries=retry, pool_maxsize=16)
    s.mount("https://", ad)
    s.headers.update({"User-Agent": "teachers-spatial-calibration/1.0 "
                                    "(academic research; jpryan7@wisc.edu)"})
    return s


def pull(sess: requests.Session, base: str, states: list[str], key: str,
         get_clause: str, workers: int, label: str) -> pd.DataFrame:
    """One `get=` clause for all three district layers, state by state."""
    def scrub(x: object) -> str:
        return str(x).replace(key, "<key>")

    frames: list[pd.DataFrame] = []
    for layer in LAYERS:
        def one(state: str):
            url = (f"{base}?get={get_clause}"
                   f"&for={layer.replace(' ', '%20')}:*&in=state:{state}&key={key}")
            try:
                resp = sess.get(url, timeout=300)
                if resp.status_code == 204:
                    return state, [], ""
                if resp.status_code != 200:
                    return state, None, f"{state}:http{resp.status_code}"
                return state, resp.json(), ""
            except Exception as exc:               # noqa: BLE001
                return state, None, f"{state}:{scrub(exc)}"

        t0 = time.time()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(one, states))
        rows, header, bad = [], None, []
        for st, data, err in results:
            if err:
                bad.append(err)
            elif data:
                header = header or data[0]
                rows.extend(data[1:])
        if header is None:
            LOG.warning("%s %s: nothing returned (%s)", label, layer, "; ".join(bad[:3]))
            continue
        df = pd.DataFrame(rows, columns=header)
        df["leaid"] = df["state"].str.zfill(2) + df[layer].str.zfill(5)
        df = df.drop(columns=[layer])
        df["sd_layer"] = layer
        frames.append(df)
        LOG.info("%s %s: %d districts in %.0fs%s", label, layer, len(df), time.time() - t0,
                 f" ({len(bad)} states failed: {'; '.join(bad[:5])})" if bad else "")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=2018, help="ACS 5-year vintage")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--outdir", default=str(HERE / "raw"))
    ap.add_argument("--no-bins", action="store_true",
                    help="skip the B20001 / B20005 earnings-bin pulls")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        print("CENSUS_API_KEY is not set (source data/spatial/.env)", file=sys.stderr)
        return 2

    sess = make_session()
    base = f"https://api.census.gov/data/{args.year}/acs/acs5"
    r = sess.get(f"{base}?get=NAME&for=state:*&key={key}", timeout=60)
    r.raise_for_status()
    states = [row[-1] for row in r.json()[1:]]
    getvars = get_vars()
    LOG.info("year %d: %d states, %d variables", args.year, len(states), len(getvars))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def finish(df: pd.DataFrame, keep_prefix: str | None) -> pd.DataFrame:
        df["year"] = args.year
        for c in df.columns:
            if c[:1] == "B" and c[-1] in "EM":
                df[c] = pd.to_numeric(df[c], errors="coerce")
        if keep_prefix:      # whole-group pulls also return M, EA, MA columns: keep E only
            drop = [c for c in df.columns
                    if c[:1] == "B" and not (c.endswith("E") and not c.endswith("EA"))]
            df = df.drop(columns=drop)
        return df

    out = finish(pull(sess, base, states, key, ",".join(["NAME"] + getvars),
                      args.workers, "medians"), None)
    if out.empty:
        return 1
    path = outdir / f"acs_earnings_school_districts_{args.year}.parquet"
    out.to_parquet(path, index=False)
    LOG.info("wrote %s (%d rows x %d cols)", path.name, len(out), out.shape[1])

    if not args.no_bins:
        parts = []
        for grp in ("B20001", "B20005"):
            df = finish(pull(sess, base, states, key, f"NAME,group({grp})",
                             args.workers, grp), grp)
            if df.empty:
                LOG.warning("%s: no rows", grp)
                continue
            parts.append(df.set_index(["leaid", "sd_layer"]).drop(
                columns=["NAME", "GEO_ID", "state", "year"], errors="ignore"))
        if parts:
            bins = pd.concat(parts, axis=1).reset_index()
            bins["year"] = args.year
            bpath = outdir / f"acs_earnings_bins_school_districts_{args.year}.parquet"
            bins.to_parquet(bpath, index=False)
            LOG.info("wrote %s (%d rows x %d cols)", bpath.name, len(bins), bins.shape[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
