#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acs_housing.py -- ACS 5-year housing-cost tables by school district, to discipline
a recurring residence cost (the goods cost of living in location 2, the suburb) in
the spatial calibration.

The model's residence cost is p_2 - p_1 = r_p * Cbar, a fixed (non-homothetic) goods
cost of living in the advantaged location.  `housing_costs.py` turns the
suburb-minus-city gaps in these tables into dollars and into units of Cbar.  This
script only acquires the tables, for the three school-district layers (unified,
elementary, secondary), one state per request as `acs_earnings.py` does (the
district layers refuse `in=state:*` before the 2022 vintage).  The pull machinery
(`pull`, `make_session`, `LAYERS`) is imported from `acs_earnings.py`, which scrubs
the API key from every error string.

Tables (ACS 5-year, default 2018 = 2014-2018, dollars of 2018):
    B25105_001  median monthly housing costs, occupied units (rent + owner costs)
    B25064_001  median gross rent (renter-occupied units paying cash rent)
    B25088_001  median selected monthly owner costs, all owner units with owner costs
    B25088_002  ... housing units with a mortgage
    B25088_003  ... housing units without a mortgage
    B25077_001  median value of owner-occupied units
    B25003_001  occupied housing units (a weight)
    B25003_002  owner-occupied units
    B25003_003  renter-occupied units
    B25104_001-017  distribution of MONTHLY HOUSING COSTS (counts by bin): less than
                $100, $100-199, ..., $900-999, $1,000-1,499, $1,500-1,999,
                $2,000-2,499, $2,500-2,999, $3,000 or more, No cash rent.  The bins
                let a LOCATION's median be read off pooled counts instead of from an
                average of district medians (as the earnings bins do in
                `district_bootstrap.py`).

Medians come with their margins of error.  Census annotation codes (-666666666 no
sample, -333333333 median in the open-ended interval, -222222222, -555555555,
-888888888, -999999999) are kept exactly as published in the raw files; they are
negative and are masked as missing before any estimate (`housing_costs.py` does).

Output (under raw/, which is gitignored):
    raw/acs_housing_school_districts_<year>.parquet    one row per district with
        `leaid` (state FIPS + 5-digit district code), `sd_layer`, `year`: medians,
        margins of error and the B25003 counts.
    raw/acs_housing_bins_school_districts_<year>.parquet   B25104 bin counts, same
        keys.

Usage
    set -a; source .env; set +a          # CENSUS_API_KEY
    python3 acs_housing.py               # vintage 2018
    python3 acs_housing.py --year 2018 --workers 8

The key is read from the environment and is never printed or written: error
strings are scrubbed of it before they reach the log.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import acs_earnings as ae  # noqa: E402  (pull machinery; key scrubbing lives there)

LOG = logging.getLogger("acs_housing")

MEDIANS = {                      # estimate + margin of error requested
    "B25105_001": "median monthly housing costs, occupied units",
    "B25064_001": "median gross rent",
    "B25088_001": "median selected monthly owner costs, all owner units",
    "B25088_002": "median selected monthly owner costs, units with a mortgage",
    "B25088_003": "median selected monthly owner costs, units without a mortgage",
    "B25077_001": "median value, owner-occupied units",
}
COUNTS = {                       # estimate only
    "B25003_001": "occupied housing units",
    "B25003_002": "owner-occupied housing units",
    "B25003_003": "renter-occupied housing units",
}
N_BINS = 17                      # B25104_001 (total) ... _017 (no cash rent)
BIN_LABELS = ["total", "<100", "100-199", "200-299", "300-399", "400-499", "500-599",
              "600-699", "700-799", "800-899", "900-999", "1000-1499", "1500-1999",
              "2000-2499", "2500-2999", "3000+", "no cash rent"]


def bin_vars() -> list[str]:
    return [f"B25104_{j:03d}E" for j in range(1, N_BINS + 1)]


def median_vars() -> list[str]:
    return [f"{k}E" for k in MEDIANS] + [f"{k}M" for k in MEDIANS] + [f"{k}E" for k in COUNTS]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=2018, help="ACS 5-year vintage")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--outdir", default=str(HERE / "raw"))
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s",
                        datefmt="%H:%M:%S")
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        print("CENSUS_API_KEY is not set (source data/spatial/.env)", file=sys.stderr)
        return 2

    sess = ae.make_session()
    base = f"https://api.census.gov/data/{args.year}/acs/acs5"
    r = sess.get(f"{base}?get=NAME&for=state:*&key={key}", timeout=60)
    r.raise_for_status()
    states = [row[-1] for row in r.json()[1:]]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    LOG.info("year %d: %d states, %d median/count variables, %d bin variables", args.year,
             len(states), len(median_vars()), N_BINS)

    def finish(df: pd.DataFrame) -> pd.DataFrame:
        df["year"] = args.year
        for c in df.columns:
            if c[:1] == "B" and c[-1] in "EM":
                df[c] = pd.to_numeric(df[c], errors="coerce")
        return df

    # one request per state and layer carries the medians, counts and bins together
    # (33 variables including NAME, below the API limit of 50); they are split into two files below.
    allvars = median_vars() + bin_vars()
    out = finish(ae.pull(sess, base, states, key, ",".join(["NAME"] + allvars),
                         args.workers, "housing"))
    if out.empty:
        LOG.error("nothing returned")
        return 1
    keys = ["leaid", "sd_layer"]
    med_cols = ["NAME", "state"] + keys + ["year"] + median_vars()
    med = out[[c for c in med_cols if c in out.columns]]
    path = outdir / f"acs_housing_school_districts_{args.year}.parquet"
    med.to_parquet(path, index=False)
    LOG.info("wrote %s (%d rows x %d cols)", path.name, len(med), med.shape[1])

    bins = out[keys + ["year"] + bin_vars()]
    bpath = outdir / f"acs_housing_bins_school_districts_{args.year}.parquet"
    bins.to_parquet(bpath, index=False)
    LOG.info("wrote %s (%d rows x %d cols)", bpath.name, len(bins), bins.shape[1])

    # completeness check against the earnings pull (same school-district universe)
    epath = outdir / f"acs_earnings_school_districts_{args.year}.parquet"
    if epath.exists():
        e = pd.read_parquet(epath, columns=["leaid", "sd_layer"])
        same = set(map(tuple, e.to_numpy())) == set(map(tuple, out[keys].to_numpy()))
        LOG.info("same (leaid, layer) universe as the earnings pull: %s (%d vs %d rows)",
                 same, len(e), len(out))
    # share of districts with a published (non-annotation) estimate, per series
    for k in MEDIANS:
        x = out[f"{k}E"]
        LOG.info("  %-11s valid %5.1f%%  annotation codes %5.1f%%  (%s)", k,
                 100 * float((x > 0).mean()), 100 * float((x < 0).mean()), MEDIANS[k])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
