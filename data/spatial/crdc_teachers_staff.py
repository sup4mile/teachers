#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
crdc_teachers_staff.py -- school-level CRDC teacher and staff counts from the Urban
Institute Education Data API, for the exploratory teacher-allocation check
(`crdc_teacher_allocation.py`).

Endpoint:  https://educationdata.urban.org/api/v1/schools/crdc/teachers-staff/{year}/
(Civil Rights Data Collection, one row per school; `year` is the first calendar
year of the school year, so 2017 = 2017-18 and 2015 = 2015-16).  Each record has
the school's NCES id (`ncessch`), its district (`leaid`) and state (`fips`) and

    teachers_fte_crdc         full-time-equivalent (FTE) teachers
    teachers_certified_fte    FTE certified teachers
    teachers_uncertified_fte  FTE uncertified teachers
    teachers_first_year_fte   FTE first-year teachers
    teachers_second_year_fte  FTE second-year teachers
    teachers_current_sy / teachers_previous_sy / teachers_absent_fte
    counselors_fte, psychologists_fte, social_workers_fte, nurses_fte,
    security_guard_fte, law_enforcement_fte / _ind, teachers_current_male / female

The API returns the numbers as strings; negative values are Urban/CRDC codes for
missing, not applicable or suppressed (-1, -2, -3, -9, ...).  They are kept exactly
as published in the raw file and masked by the analysis.

The portal pages at 10,000 records; the script follows the `next` links one request
at a time with a short pause (the API publishes no rate limit), retries with backoff
on 429/5xx, checks the row count against the API's own `count`, and writes
    raw/crdc_teachers_staff_<year>.parquet
(raw/ is gitignored).  No key is needed.

Usage
    python3 crdc_teachers_staff.py                 # 2017 (school year 2017-18)
    python3 crdc_teachers_staff.py --year 2015
    python3 crdc_teachers_staff.py --probe         # only report whether the year is served

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

LOG = logging.getLogger("crdc_teachers_staff")
HERE = Path(__file__).resolve().parent
BASE = "https://educationdata.urban.org/api/v1/schools/crdc/teachers-staff"
ID_COLS = ("crdc_id", "ncessch", "leaid")


def make_session() -> requests.Session:
    s = requests.Session()
    retry = Retry(total=6, connect=6, read=6, status=6, backoff_factor=2.0,
                  status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset(["GET"]), raise_on_status=False)
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=2))
    s.headers.update({"User-Agent": "teachers-spatial-calibration/1.0 "
                                    "(academic research; jpryan7@wisc.edu)"})
    return s


def pull_year(sess: requests.Session, year: int, pause: float, timeout: int = 300
              ) -> tuple[pd.DataFrame, int]:
    """Page through one year; returns (rows, the API's `count`)."""
    url: str | None = f"{BASE}/{year}/"
    frames: list[pd.DataFrame] = []
    count = -1
    page = 0
    while url:
        resp = sess.get(url, timeout=timeout)
        if resp.status_code == 404:
            raise FileNotFoundError(f"404 {url}")
        resp.raise_for_status()
        payload = resp.json()
        count = int(payload.get("count", count))
        if payload.get("results"):
            frames.append(pd.DataFrame(payload["results"]))
        url = payload.get("next")
        page += 1
        LOG.info("  page %d: %d rows so far (API count %d)", page,
                 sum(map(len, frames)), count)
        if url:
            time.sleep(pause)
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), count


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=2017, help="2017 = school year 2017-18")
    ap.add_argument("--pause", type=float, default=1.0, help="seconds between page requests")
    ap.add_argument("--outdir", default=str(HERE / "raw"))
    ap.add_argument("--probe", action="store_true", help="only check that the year is served")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    sess = make_session()
    if args.probe:
        r = sess.get(f"{BASE}/{args.year}/", timeout=120)
        ok = r.status_code == 200
        LOG.info("%s: HTTP %d%s", args.year, r.status_code,
                 f", {r.json()['count']} schools" if ok else "")
        return 0 if ok else 1
    try:
        df, count = pull_year(sess, args.year, args.pause)
    except FileNotFoundError as exc:
        LOG.error("not served: %s", exc)
        return 1
    if df.empty:
        LOG.error("no rows returned")
        return 1
    for c in df.columns:
        if c in ID_COLS:
            df[c] = df[c].astype(str)
        elif c not in ("year", "fips"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df["year"] = args.year
    LOG.info("rows %d (API count %d); schools %d unique; districts %d; states %d",
             len(df), count, df.ncessch.nunique(), df.leaid.nunique(), df.fips.nunique())
    if len(df) != count:
        LOG.warning("row count differs from the API count")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"crdc_teachers_staff_{args.year}.parquet"
    df.to_parquet(path, index=False)
    LOG.info("wrote %s (%d rows x %d cols)", path.name, *df.shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
