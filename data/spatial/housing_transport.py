#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
housing_transport.py -- evidence for income-driven demand for housing and an endogenous
housing market in the spatial model (DESIGN.md M1/M1b, task D2).  Five pieces:

 (a) the suburb-minus-city gap in the monthly housing-cost DISTRIBUTION (10th, 25th, 50th,
     75th, 90th percentiles) from the ACS 5-year 2014-18 school-district bin tables: all
     units (B25104, cached by `acs_housing.py`), renters (B25063 gross rent) and owners
     (B25087 selected monthly owner costs, with and without a mortgage, pooled on a common
     grid).  Same panel, locale partition (188 CZs), CZ weights and whole-CZ bootstrap draws
     (B = 1999, seed 20260929) as `housing_costs.py` / `district_settings.py`; the all-units
     median gap must reproduce housing_costs.md (0.114, SE 0.017).  Question: does the
     suburb lack cheap units (a larger gap at the bottom)?  Also who occupies the cheapest
     units (renters, owners with a mortgage, owners without): the bottom of the pooled
     distribution is mostly outright owners, so the renter profile is the validation moment
     for the lower tail of R H.
 (b) the housing expenditure share by income: BLS Consumer Expenditure Survey 2018, quintiles
     and deciles of income before taxes (Tables 1101 and 1110): shares of average annual
     expenditures for shelter, utilities/fuels/public services and housing; Stone-Geary fit
     s(Y) = h/Y + aH (1 - h/Y) (Y = total expenditures) by nonlinear least squares, h as a
     share of mean expenditures ($61,224); tercile shares for the model's
     `housing_share_by_income`; cross-check with ACS national housing cost burdens by
     household income (B25074 renters, B25095 owners, B25106).  Household size rises with
     income across the CE groups (1.7 to 3.2 persons per consumer unit), so the fit is
     repeated on EQUIVALIZED expenditure Y / n^theta (square root of persons, OECD-modified,
     per adult, per capita; the central choice is the square root scale) and the CU-level fit
     (no size adjustment) is kept as the other end; CE Table 1400 (size of consumer unit) gives
     a one-person cross-check.  The report states which aggregate the model's `housing_share`
     must compute (a ratio of aggregates, like the CE shares).
 (c) housing supply elasticities: Saiz (2010, QJE) MSA estimates (author's data, 269 MSAs,
     Table VI = the 95 with population over 500,000) and Baum-Snow and Han (2024, JPE 132(6))
     tract-level predicted elasticities (authors' September 2023 release): means by distance
     to the CBD, a city/suburb split matched to the pupil shares of the locale partition,
     developed fraction; a recommended long-run eps_H for the city and the suburb.  ONE horizon
     convention is applied to every anchor (the long-run stock elasticity of Baum-Snow and
     Duranton's slide 76 at 3.5% depreciation: x1.52 for Saiz's 30-year, x3.34 for the
     10-year estimates); Baum-Snow and Han's Online Appendix D (Saiz's method on 2000-10 data)
     is used for a period check.
 (d) transport: ACS 2014-18 school-district B25044 (tenure by vehicles available), B25046
     (aggregate vehicles) and B08201 (household size by vehicles): zero-vehicle share and
     vehicles per household by location (locale partition + bootstrap); within-CZ district
     regressions on a suburb indicator with and without log median household income and
     within HOUSEHOLD-SIZE classes (1, 2, 3, 4+ persons: suburban households are larger);
     AAA "Your Driving Costs" 2018 and the CE 2018 transportation lines converted into a
     fixed goods cost in units of mean consumption ($61,224).
 (e) the income gradient of willingness to pay in BFM (2007, JPE 115(4), Table 8): the
     elasticity of WTP with respect to household income at the means and its SE, and the
     elasticity that also counts the college-degree term (and the children term) through the
     1990 Bay Area association with household income, measured in the Census 1990 PUMS 5%
     sample (six counties), for bachelor's and for associate-or-more definitions of college.

Everything downloaded is cached under raw/housing_transport/ (gitignored); `fetch` fills it
(idempotent), `build` computes.  The Census key is read from the environment
(CENSUS_API_KEY, `set -a; source .env; set +a`) and never printed or written; error strings
are scrubbed of it (the pull machinery is `acs_earnings.py`'s).

Writes (under estimates/settings/):
    housing_transport.md, housing_transport.json, housing_transport_targets.toml, and
    housing_transport_draws.csv (the bootstrap draws of the main gap statistics, on the same
    CZ multiplicities as `housing_costs_draws.csv` and `district_settings_draws_locale.csv`,
    so correlations with the district targets can be formed)

Usage (from data/spatial; on a compute node, see the project notes)
    set -a; source .env; set +a
    uv run python housing_transport.py fetch            # ACS tables, BLS CE (incl. Table 1400), papers, slides, PUMS (a minute)
    uv run python housing_transport.py build            # B = 1999, seed 20260929 (about 30 s)
    uv run python housing_transport.py build --B 199 --outdir <dir>    # quick check
    uv run python housing_transport.py build --parts a b c d e         # subsets (d computes the CE tables it needs)

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import math
import os
import re
import sys
import zipfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

LOG = logging.getLogger("housing_transport")

RAW = HERE / "raw" / "housing_transport"
YEAR = 2018                                   # ACS 5-year vintage 2014-18, dollars of 2018
SEED = 20260929
B_DEFAULT = 1999
UA_ACADEMIC = "teachers-spatial-calibration/1.0 (academic research; jpryan7@wisc.edu)"
UA_BROWSER = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
              "Chrome/120.0 Safari/537.36")

# Consumer Expenditure Survey 2018 (BLS): mean annual expenditures per consumer unit.
CE_TOTAL_EXP = 61224.13                       # Table R-1 (housing_costs.py)
CE_BASE = "https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error"
CE_FILES = {
    "ce_quintiles_xlsx": f"{CE_BASE}/cu-income-quintiles-before-taxes-2018.xlsx",
    "ce_quintiles_pdf": f"{CE_BASE}/cu-income-quintiles-before-taxes-2018.pdf",
    "ce_deciles_xlsx": f"{CE_BASE}/cu-income-deciles-before-taxes-2018.xlsx",
    "ce_deciles_pdf": f"{CE_BASE}/cu-income-deciles-before-taxes-2018.pdf",
    "ce_size_xlsx": f"{CE_BASE}/cu-size-2018.xlsx",                    # Table 1400, size of consumer unit
    "ce_size_pdf": f"{CE_BASE}/cu-size-2018.pdf",
}
# Other resources: name -> (url, file name under RAW, minimum plausible size in bytes, magic)
RESOURCES: dict[str, tuple[str, str, int, bytes]] = {
    "ce_quintiles_xlsx": (CE_FILES["ce_quintiles_xlsx"], "cu-income-quintiles-before-taxes-2018.xlsx",
                          10_000, b"PK"),
    "ce_quintiles_pdf": (CE_FILES["ce_quintiles_pdf"], "cu-income-quintiles-before-taxes-2018.pdf",
                         10_000, b"%PDF"),
    "ce_deciles_xlsx": (CE_FILES["ce_deciles_xlsx"], "cu-income-deciles-before-taxes-2018.xlsx",
                        10_000, b"PK"),
    "ce_deciles_pdf": (CE_FILES["ce_deciles_pdf"], "cu-income-deciles-before-taxes-2018.pdf",
                       10_000, b"%PDF"),
    "ce_size_xlsx": (CE_FILES["ce_size_xlsx"], "cu-size-2018.xlsx", 10_000, b"PK"),
    "ce_size_pdf": (CE_FILES["ce_size_pdf"], "cu-size-2018.pdf", 10_000, b"%PDF"),
    # N. Baum-Snow and G. Duranton, 'Housing Supply and Housing Affordability' (Handbook of Regional
    # and Urban Economics vol. 6 chapter; slides of the Urban Economics Association summer school, 2025)
    "bd_slides": ("https://urbaneconomics.org/workshops/summerschool2025/slides/baumsnow.pdf",
                  "BaumSnowDuranton_slides_2025.pdf", 500_000, b"%PDF"),
    "aaa_2018": ("https://exchange.aaa.com/wp-content/uploads/2018/09/"
                 "18-0090_2018-Your-Driving-Costs-Brochure_FNL-Lo-5-2.pdf",
                 "AAA_Your_Driving_Costs_2018.pdf", 100_000, b"%PDF"),
    "saiz_data": ("https://drive.google.com/uc?export=download&id=1RnTT3f6w2LHjH7iPKnXqlTTkeHeMQ8ea",
                  "saiz_elasticity.dta", 10_000, b""),
    "saiz_qje": ("https://joeornstein.github.io/pols-4641/readings/"
                 "Saiz%20-%202010%20-%20The%20Geographic%20Determinants%20of%20Housing%20Supply.pdf",
                 "Saiz2010_QJE.pdf", 100_000, b"%PDF"),
    "bsh_paper": ("https://drive.google.com/uc?export=download&id=1L4iUwSdUjqy1VMsvalJ6SgAWeYDkVwU2",
                  "BaumSnowHan_2023_accepted.pdf", 100_000, b"%PDF"),
    "bsh_estimates": ("https://drive.google.com/uc?export=download&id=1Qn3PZ0u8PO2IbyyPETCrca3rGmPyTVL7",
                      "BaumSnowHan_elasticity_estimates.zip", 1_000_000, b"PK"),
    "bfm_jpe": ("https://matthewturner.org/ec2410/readings/Bayer_Ferreira_Mcmillan_JPE_2007.pdf",
                "BFM_JPE_2007.pdf", 100_000, b"%PDF"),
    "pums90_ca": ("https://www2.census.gov/census_1990/1990_PUMS_A/PUMSAXCA.zip",
                  "PUMSAXCA.zip", 10_000_000, b"PK"),
    "pums90_dd": ("https://www2.census.gov/census_1990/1990_PUMS_A/DOCUMENT/PUMSUSDD.doc",
                  "PUMSUSDD.doc", 10_000, b""),
}

# ACS tables pulled at the school-district level (three layers, state by state)
DISTRICT_PULLS = {            # output file -> (variables, label)
    "acs_tenure_bins_school_districts": (
        [f"B25063_{j:03d}E" for j in range(1, 28)], "B25063 gross rent"),
    "acs_owner_bins_school_districts": (
        [f"B25087_{j:03d}E" for j in range(1, 40)], "B25087 owner costs"),
    "acs_vehicles_school_districts": (
        [f"B25044_{j:03d}E" for j in range(1, 16)] + [f"B25046_{j:03d}E" for j in range(1, 4)]
        + [f"B08201_{j:03d}E" for j in range(1, 31)], "B25044/B25046/B08201 vehicles"),
}
NATIONAL_TABLES = ("B25074", "B25095", "B25106", "B19001", "B19025", "B25003", "B25044",
                   "B25046", "B08201")


# --------------------------------------------------------------------------
# Small utilities
# --------------------------------------------------------------------------

def clean(o: Any) -> Any:
    """JSON-safe copy: NaN/inf -> None, numpy -> python, arrays -> lists, `draws` dropped."""
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items() if k not in ("draws", "cz_table")}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return [clean(v) for v in o.tolist()]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def fmt(v: Any, d: int = 3) -> str:
    return "" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:,.{d}f}"


def fse(v: Any, se: Any, d: int = 3) -> str:
    if v is None or se is None or not np.isfinite(v) or not np.isfinite(se):
        return fmt(v, d)
    return f"{v:,.{d}f} ({se:,.{d}f})"


def pct(v: Any, d: int = 1) -> str:
    return "" if v is None or not np.isfinite(v) else f"{100 * v:.{d}f}"


def boot_se_of(draws: np.ndarray) -> tuple[float, float, float]:
    """(SD, 2.5th, 97.5th percentile) of the finite bootstrap draws."""
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    lo, hi = np.percentile(d, [2.5, 97.5])
    return float(np.std(d, ddof=1)), float(lo), float(hi)


# --------------------------------------------------------------------------
# fetch: caches every external input under raw/housing_transport/
# --------------------------------------------------------------------------

def _scrub(x: object, key: str | None) -> str:
    s = str(x)
    return s.replace(key, "<key>") if key else s


def download(url: str, dest: Path, minsize: int, magic: bytes, refresh: bool = False) -> bool:
    """Download `url` to `dest` (skipped if present).  Tries the academic user agent, then a
    browser one (BLS, AAA and Google Drive each accept a different one).  Checks the size and
    the first bytes so that an HTML error page is never cached as data."""
    import requests
    if dest.exists() and dest.stat().st_size >= minsize and not refresh:
        LOG.info("cached   %s (%s bytes)", dest.name, f"{dest.stat().st_size:,}")
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    for ua in (UA_ACADEMIC, UA_BROWSER):
        try:
            r = requests.get(url, headers={"User-Agent": ua}, timeout=600, allow_redirects=True)
        except Exception as exc:                                     # noqa: BLE001
            LOG.warning("  %s: %s", dest.name, exc)
            continue
        html = r.content[:20].lstrip().lower().startswith((b"<!doctype", b"<html"))
        ok = (r.status_code == 200 and len(r.content) >= minsize and not html
              and r.content[:len(magic)] == magic)
        if ok:
            dest.write_bytes(r.content)
            LOG.info("fetched  %s (%s bytes)", dest.name, f"{len(r.content):,}")
            return True
        LOG.warning("  %s: http %s, %s bytes, starts %r", dest.name, r.status_code,
                    f"{len(r.content):,}", r.content[:8])
    LOG.error("could not fetch %s from %s", dest.name, url.split("?")[0])
    return False


def fetch_district_tables(year: int, workers: int, refresh: bool) -> None:
    import acs_earnings as ae
    todo = {k: v for k, v in DISTRICT_PULLS.items()
            if refresh or not (RAW / f"{k}_{year}.parquet").exists()}
    if not todo:
        LOG.info("cached   school-district ACS tables")
        return
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        raise SystemExit("CENSUS_API_KEY is not set (set -a; source .env; set +a)")
    sess = ae.make_session()
    base = f"https://api.census.gov/data/{year}/acs/acs5"
    r = sess.get(f"{base}?get=NAME&for=state:*&key={key}", timeout=60)
    r.raise_for_status()
    states = [row[-1] for row in r.json()[1:]]
    for name, (vars_, label) in todo.items():
        # B25044 + B25046 + B08201 is 48 variables + NAME: below the API limit of 50
        out = ae.pull(sess, base, states, key, ",".join(["NAME"] + vars_), workers, label)
        if out.empty:
            LOG.error("%s: nothing returned", label)
            continue
        out["year"] = year
        for c in out.columns:
            if c[:1] == "B" and c[-1] in "EM":
                out[c] = pd.to_numeric(out[c], errors="coerce")
        keep = ["leaid", "sd_layer", "year"] + vars_
        path = RAW / f"{name}_{year}.parquet"
        out[keep].to_parquet(path, index=False)
        LOG.info("wrote    %s (%d rows x %d cols)", path.name, len(out), len(keep))


def fetch_national_tables(refresh: bool) -> None:
    import acs_earnings as ae
    missing = [(prod, tbl) for prod in ("acs5", "acs1") for tbl in NATIONAL_TABLES
               if refresh or not (RAW / "acs_national" / f"{prod}_{YEAR}_{tbl}.json").exists()]
    if not missing:
        LOG.info("cached   national ACS tables")
        return
    key = os.environ.get("CENSUS_API_KEY")
    if not key:
        raise SystemExit("CENSUS_API_KEY is not set (set -a; source .env; set +a)")
    sess = ae.make_session()
    (RAW / "acs_national").mkdir(parents=True, exist_ok=True)
    for prod, label in (("acs5", "5-year 2014-18"), ("acs1", "1-year 2018")):
        for tbl in NATIONAL_TABLES:
            path = RAW / "acs_national" / f"{prod}_{YEAR}_{tbl}.json"
            if path.exists() and not refresh:
                continue
            url = f"https://api.census.gov/data/{YEAR}/acs/{prod}?get=NAME,group({tbl})&for=us:1&key={key}"
            try:
                r = sess.get(url, timeout=120)
                r.raise_for_status()
                data = r.json()
            except Exception as exc:                                   # noqa: BLE001
                LOG.warning("national %s %s: %s", prod, tbl, _scrub(exc, key))
                continue
            head, row = data[0], data[1]
            rec = {h: v for h, v in zip(head, row)}
            path.write_text(json.dumps(rec))
        LOG.info("national ACS %s tables cached", label)


def cmd_fetch(args: argparse.Namespace) -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    bad = []
    for name, (url, fn, minsize, magic) in RESOURCES.items():
        if not download(url, RAW / fn, minsize, magic, args.refresh):
            bad.append(name)
    fetch_district_tables(YEAR, args.workers, args.refresh)
    fetch_national_tables(args.refresh)
    # the Baum-Snow-Han archive holds the tract file; unpack it once
    z = RAW / "BaumSnowHan_elasticity_estimates.zip"
    out = RAW / "BaumSnowHan"
    if z.exists() and not (out / "gammas_hat_all.dta").exists():
        with zipfile.ZipFile(z) as zf:
            zf.extractall(out)
        LOG.info("unpacked %s", z.name)
    if bad:
        LOG.error("missing: %s", ", ".join(bad))
        return 1
    return 0


# --------------------------------------------------------------------------
# Panel, partition and the CZ weights / bootstrap draws (as housing_costs.py)
# --------------------------------------------------------------------------

def clean_acs_counts(a: pd.DataFrame, cols: list[str], total: str) -> pd.DataFrame:
    """One row per district id: the count columns with Census annotation codes (negative)
    masked; a district whose table total is invalid has the whole table set to NaN, so a
    pooled bin is never a sum over a different set of districts than the others."""
    a = a[~a.leaid.str.endswith("99999")].drop_duplicates("leaid").copy()
    x = a[cols].where(a[cols] >= 0)
    bad = ~(a[total] >= 0) | x.isna().any(axis=1)
    x.loc[bad, :] = np.nan
    x.insert(0, "leaid", a.leaid.to_numpy())
    return x


def load_panel(indir: Path) -> pd.DataFrame:
    """District panel of data_estimate.assemble with the ACS earnings columns and the B25104
    housing-cost bins (as housing_costs.py), plus the tenure-specific cost bins (B25063,
    B25087) and the vehicle tables (B25044, B25046, B08201) of `fetch`."""
    import district_bootstrap as db
    import housing_costs as hc
    panel = db.build_panel(indir)
    panel = db.load_income_columns(panel, indir / "raw")
    panel = hc.load_housing_columns(panel, indir / "raw")
    attrs = dict(panel.attrs)
    # B25003 occupied units (weights) already merged by hc (acsh_occ, acsh_own, acsh_rent)
    for fn, total in (("acs_tenure_bins_school_districts", "B25063_001E"),
                      ("acs_owner_bins_school_districts", "B25087_001E"),
                      ("acs_vehicles_school_districts", "B25044_001E")):
        path = RAW / f"{fn}_{YEAR}.parquet"
        if not path.exists():
            raise SystemExit(f"no {path.name}: run `housing_transport.py fetch` first")
        a = pd.read_parquet(path)
        cols = [c for c in a.columns if c[:1] == "B" and c.endswith("E")]
        x = clean_acs_counts(a, cols, total)
        panel = panel.merge(x, on="leaid", how="left", suffixes=("", f"_{fn}"))
    panel.attrs.update(attrs)
    return panel


def cz_table(panel: pd.DataFrame, scheme: str = "locale") -> tuple[pd.DataFrame, pd.DataFrame]:
    """(partitioned districts, CZ table) as `district_bootstrap.cz_gap_table` with no gap
    columns: CZ rows in the order (sorted CZ code) that fixes the bootstrap multiplicities."""
    import district_bootstrap as db
    part, g = db.cz_gap_table(panel, [], scheme)
    part = part[part.year == YEAR]
    return part, g


class Boot:
    """The whole-CZ multiplicity matrix (B x n_CZ) of `district_bootstrap.draw_counts`,
    shared by every statistic of this script, and the weighted-mean statistic."""

    def __init__(self, n_cz: int, B: int, seed: int, weights: np.ndarray):
        import district_bootstrap as db
        self.counts = db.draw_counts(n_cz, B, seed)
        self.w = np.asarray(weights, float)
        self.B, self.seed, self.n_cz = B, seed, n_cz

    def wmean(self, x: np.ndarray, l1: np.ndarray | None = None,
              l2: np.ndarray | None = None) -> dict[str, Any]:
        """CZ-weighted mean of x over the CZs where it is observed, with bootstrap SE, 95%
        percentile interval and draws; `level_1/2` are same-sample weighted means of l1/l2."""
        import district_bootstrap as db
        x = np.asarray(x, float)
        st = db.boot_wmean(x, self.w, self.counts)
        se, lo, hi = boot_se_of(st["draws"])
        out = dict(value=float(st["est"]), boot_se=se, ci95_lo=lo, ci95_hi=hi, n=int(st["n"]),
                   draws=st["draws"])
        m = np.isfinite(x) & np.isfinite(self.w) & (self.w > 0)
        for tag, lev in (("level_1", l1), ("level_2", l2)):
            if lev is not None:
                y = np.asarray(lev, float)
                mm = m & np.isfinite(y)
                out[tag] = float(np.average(y[mm], weights=self.w[mm])) if mm.any() else None
        return out

    def ratio_of_means(self, num: np.ndarray, den: np.ndarray) -> dict[str, Any]:
        """(weighted mean of num) / (weighted mean of den) over CZs with both, with draws."""
        m = (np.isfinite(num) & np.isfinite(den) & np.isfinite(self.w) & (self.w > 0))
        wn, wd = np.where(m, self.w * num, 0.0), np.where(m, self.w * den, 0.0)
        est = wn.sum() / wd.sum()
        with np.errstate(invalid="ignore", divide="ignore"):
            draws = (self.counts @ wn) / (self.counts @ wd)
        se, lo, hi = boot_se_of(draws)
        return dict(value=float(est), boot_se=se, ci95_lo=lo, ci95_hi=hi, n=int(m.sum()),
                    draws=draws)

    def stat_of_draws(self, est: float, draws: np.ndarray, n: int = 0) -> dict[str, Any]:
        se, lo, hi = boot_se_of(draws)
        return dict(value=float(est), boot_se=se, ci95_lo=lo, ci95_hi=hi, n=n, draws=draws)


# --------------------------------------------------------------------------
# (a) The monthly housing-cost distribution by location
# --------------------------------------------------------------------------

PCTS = (0.10, 0.25, 0.50, 0.75, 0.90)
PCT_KEYS = tuple(f"p{int(100 * p)}" for p in PCTS)
INF = np.inf

# B25104 (all occupied units): 15 bins, cost = gross rent or selected monthly owner cost
EDGES_ALL = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500,
                      3000, INF], float)
COLS_ALL = [f"bin_h_{j:02d}" for j in range(15)]
# B25063 (renter-occupied units paying cash rent), gross rent: _003 .. _026
EDGES_RENT = np.array([0, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                       800, 900, 1000, 1250, 1500, 2000, 2500, 3000, 3500, INF], float)
COLS_RENT = [f"B25063_{j:03d}E" for j in range(3, 27)]
# B25087 (owner-occupied units), selected monthly owner costs: with a mortgage _003 .. _019,
# without a mortgage _021 .. _039
EDGES_OWN_M = np.array([0, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 2500,
                        3000, 3500, 4000, INF], float)
COLS_OWN_M = [f"B25087_{j:03d}E" for j in range(3, 20)]
EDGES_OWN_N = np.array([0, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000,
                        1100, 1200, 1300, 1400, 1500, INF], float)
COLS_OWN_N = [f"B25087_{j:03d}E" for j in range(21, 40)]
# common grid of the two owner groups: the union of their finite edges
EDGES_OWN = np.array(sorted(set(EDGES_OWN_M[:-1]) | set(EDGES_OWN_N[:-1])) + [INF], float)


def split_matrix(edges: np.ndarray, union: np.ndarray, open_alloc: np.ndarray | None) -> np.ndarray:
    """(n_union x n_bins) matrix M with union_counts = M @ counts: each finite bin's count is
    split over the union bins it contains in proportion to width (uniform within a bin); the
    open top bin's count is split over the union bins at or above its lower edge in the
    proportions `open_alloc` (a vector over those union bins; the identity if only one)."""
    nb, nu = len(edges) - 1, len(union) - 1
    M = np.zeros((nu, nb))
    for j in range(nb - 1):
        lo, hi = edges[j], edges[j + 1]
        for k in range(nu):
            a, b = max(lo, union[k]), min(hi, union[k + 1])
            if b > a:
                M[k, j] = (b - a) / (hi - lo)
    top = [k for k in range(nu) if union[k] >= edges[nb - 1]]
    if open_alloc is None:
        if len(top) != 1:
            raise ValueError("open bin spans several union bins: an allocation is needed")
        M[top[0], nb - 1] = 1.0
    else:
        M[top, nb - 1] = open_alloc
    return M


def pct_from_bins(counts: np.ndarray, edges: np.ndarray, ps: Iterable[float],
                  top: str = "pareto") -> tuple[np.ndarray, np.ndarray]:
    """Percentiles of a distribution known through counts in `edges` intervals (last edge inf),
    by linear interpolation inside the bin holding the percentile (the rule of
    `housing_costs.grouped_median`).  If the percentile lies in the open top bin: `top` =
    'pareto' extrapolates a Pareto tail through the survival function at the lower edge of the
    top bin and of the nearest non-empty bin below it (shape alpha >= 1), 'lower' returns the bin's
    lower edge, 'nan' NaN.
    Returns (values, flags) with flags = 1 where the top bin was used."""
    ps = list(ps)
    val, flag = np.full(len(ps), np.nan), np.zeros(len(ps))
    c = np.asarray(counts, float)
    tot = np.nansum(c)
    if not np.isfinite(tot) or tot <= 0 or np.isnan(c).any():
        return val, flag
    cum = np.cumsum(c)
    k = len(c)
    for i, p in enumerate(ps):
        target = p * tot
        j = min(int(np.searchsorted(cum, target)), k - 1)
        below = cum[j - 1] if j > 0 else 0.0
        if np.isinf(edges[j + 1]):
            flag[i] = 1
            if top == "lower":
                val[i] = edges[j]
            elif top == "pareto" and c[j] > 0:
                m = j - 1
                while m > 0 and c[m] <= 0:
                    m -= 1
                if m >= 1 and c[m] > 0:
                    s_top, s_m = c[j] / tot, c[m:].sum() / tot
                    alpha = math.log(s_m / s_top) / math.log(edges[j] / edges[m])
                    alpha = max(alpha, 1.0)
                    val[i] = edges[j] * (s_top / (1.0 - p)) ** (1.0 / alpha)
            continue
        val[i] = edges[j] + (edges[j + 1] - edges[j]) * (target - below) / c[j]
    return val, flag


def cdf_at(counts: np.ndarray, edges: np.ndarray, x: float) -> float:
    """Share of the distribution below x (linear interpolation inside bins; x inside the open
    top bin is treated as at its lower edge, i.e. the share below the lower edge)."""
    c = np.asarray(counts, float)
    tot = np.nansum(c)
    if not np.isfinite(tot) or tot <= 0 or np.isnan(c).any() or not np.isfinite(x):
        return np.nan
    j = int(np.searchsorted(edges, x, side="right")) - 1
    j = min(max(j, 0), len(c) - 1)
    below = c[:j].sum()
    if np.isinf(edges[j + 1]):
        return float(below / tot)
    frac = (x - edges[j]) / (edges[j + 1] - edges[j])
    return float((below + frac * c[j]) / tot)


def pooled_by_location(part: pd.DataFrame, cols: list[str], czones: np.ndarray) -> np.ndarray:
    """(n_cz, 2, k) bin counts summed over the districts of each location of each CZ; NaN
    where no district of the location has the table."""
    sums = part.groupby(["czone", "loc_id"])[cols].sum(min_count=1)
    out = np.full((len(czones), 2, len(cols)), np.nan)
    pos = {cz: i for i, cz in enumerate(czones)}
    for (cz, loc), row in sums.iterrows():
        if cz in pos:
            out[pos[cz], int(loc) - 1] = row.to_numpy(float)
    return out


def distribution_blocks(part: pd.DataFrame, czones: np.ndarray) -> dict[str, dict[str, Any]]:
    """Pooled counts by CZ x location x bin for each housing-cost distribution, on its own
    grid: all units (B25104), renters (B25063), owners (B25087, common grid), and the two
    owner subgroups."""
    blocks: dict[str, dict[str, Any]] = {}
    blocks["all"] = dict(edges=EDGES_ALL, label="all occupied units (B25104: gross rent or selected "
                         "monthly owner costs)", counts=pooled_by_location(part, COLS_ALL, czones))
    blocks["rent"] = dict(edges=EDGES_RENT, label="renter-occupied units paying cash rent (B25063 "
                          "gross rent)", counts=pooled_by_location(part, COLS_RENT, czones))
    cm = pooled_by_location(part, COLS_OWN_M, czones)
    cn = pooled_by_location(part, COLS_OWN_N, czones)
    blocks["own_mort"] = dict(edges=EDGES_OWN_M, label="owner-occupied units with a mortgage (B25087 "
                              "selected monthly owner costs)", counts=cm)
    blocks["own_nomort"] = dict(edges=EDGES_OWN_N, label="owner-occupied units without a mortgage "
                                "(B25087)", counts=cn)
    # owners pooled on the union grid.  The no-mortgage open bin ($1,500 or more) is spread over
    # the union bins at or above $1,500 like the with-mortgage owners of the partitioned districts
    # (all 188 CZs and both locations pooled) there.
    nu = len(EDGES_OWN) - 1
    top_idx = [k for k in range(nu) if EDGES_OWN[k] >= 1500]
    Mm_full = split_matrix(EDGES_OWN_M, EDGES_OWN, None)
    nat_m = np.nansum(part[COLS_OWN_M].to_numpy(float), axis=0)
    nat_union = Mm_full @ nat_m
    nat = nat_union[top_idx]
    alloc = nat / nat.sum()
    Mn = split_matrix(EDGES_OWN_N, EDGES_OWN, alloc)
    own = np.einsum("uj,clj->clu", Mm_full, cm) + np.einsum("uj,clj->clu", Mn, cn)
    blocks["own"] = dict(edges=EDGES_OWN, label="all owner-occupied units (B25087, with and without a "
                         "mortgage pooled on a common grid)", counts=own,
                         note_alloc=[float(a) for a in alloc], alloc_edges=[float(EDGES_OWN[k]) for k in top_idx])
    # all units on a finer grid: renters (B25063) + owners (B25087) on the union of their grids, the open
    # renter bin ($3,500 or more) spread over the union bins above it like the with-mortgage owners of the
    # partitioned districts.
    # A cross-check of the coarse B25104 bins (top bin $3,000 or more) in the upper tail.
    U = np.array(sorted(set(EDGES_RENT[:-1]) | set(EDGES_OWN[:-1])) + [INF], float)
    nat_u = split_matrix(EDGES_OWN_M, U, None) @ nat_m
    top_r = [k for k in range(len(U) - 1) if U[k] >= 3500]
    top_n = [k for k in range(len(U) - 1) if U[k] >= 1500]
    Mr = split_matrix(EDGES_RENT, U, nat_u[top_r] / nat_u[top_r].sum())
    Mm_u = split_matrix(EDGES_OWN_M, U, None)
    Mn_u = split_matrix(EDGES_OWN_N, U, nat_u[top_n] / nat_u[top_n].sum())
    fine = (np.einsum("uj,clj->clu", Mr, blocks["rent"]["counts"]) + np.einsum("uj,clj->clu", Mm_u, cm)
            + np.einsum("uj,clj->clu", Mn_u, cn))
    blocks["all_fine"] = dict(edges=U, label="all occupied units on a finer grid: renters (B25063) plus owners (B25087) "
                              "pooled (cross-check of the B25104 bins)", counts=fine)
    return blocks


def percentile_tables(blocks: dict[str, dict[str, Any]], top: str) -> dict[str, dict[str, np.ndarray]]:
    """Per distribution: P (n_cz, 2, 5) percentiles and F (same) top-bin flags."""
    out: dict[str, dict[str, np.ndarray]] = {}
    for tab, b in blocks.items():
        cnt = b["counts"]
        n, _, _ = cnt.shape
        P, F = np.full((n, 2, len(PCTS)), np.nan), np.zeros((n, 2, len(PCTS)))
        for i in range(n):
            for loc in range(2):
                P[i, loc], F[i, loc] = pct_from_bins(cnt[i, loc], b["edges"], PCTS, top)
        out[tab] = dict(P=P, F=F)
    return out


TENURE_GROUPS = (("rent", "renters (B25063 gross rent, cash renters)"),
                 ("own_mort", "owners with a mortgage (B25087)"),
                 ("own_nomort", "owners without a mortgage (B25087)"))


def tail_composition(blocks: dict[str, dict[str, Any]], czones: np.ndarray, boot: Boot) -> dict[str, Any]:
    """Who occupies the cheapest units.  For each CZ the cut is the 10th / 25th / 50th percentile of the
    all-units cost distribution (B25104) of the two locations pooled; the entry is, for each location,
    the share of the units below the cut that are rented (cash rent), owned with a mortgage or owned
    without a mortgage (each group's count below the cut by linear interpolation inside its own bins,
    B25063 / B25087).  Also the occupancy shares of the three groups.  CZ-weighted means with the shared
    whole-CZ bootstrap (gap = suburb minus city)."""
    tot = {k: blocks[k]["counts"].sum(axis=2) for k, _ in TENURE_GROUPS}          # (n_cz, 2), NaN-propagating
    all_tot = sum(tot.values())
    pooled_all = np.nansum(blocks["all"]["counts"], axis=1)
    out: dict[str, Any] = dict(groups=[k for k, _ in TENURE_GROUPS], labels={k: v for k, v in TENURE_GROUPS},
                               occupancy={}, below={}, cut_usd={})
    for k, _ in TENURE_GROUPS:
        s1, s2 = tot[k][:, 0] / all_tot[:, 0], tot[k][:, 1] / all_tot[:, 1]
        out["occupancy"][k] = boot.wmean(s2 - s1, s1, s2)
    n = len(czones)
    for q in (0.10, 0.25, 0.50):
        comp = {k: np.full((n, 2), np.nan) for k, _ in TENURE_GROUPS}
        cuts = np.full(n, np.nan)
        for ci in range(n):
            cut, _ = pct_from_bins(pooled_all[ci], blocks["all"]["edges"], [q], "pareto")
            cuts[ci] = cut[0]
            if not np.isfinite(cut[0]):
                continue
            for loc in range(2):
                below = {k: tot[k][ci, loc] * cdf_at(blocks[k]["counts"][ci, loc], blocks[k]["edges"], cut[0])
                         for k, _ in TENURE_GROUPS}
                s = sum(below.values())
                for k in below:
                    comp[k][ci, loc] = below[k] / s
        key = f"q{int(100 * q)}"
        out["below"][key] = {k: boot.wmean(comp[k][:, 1] - comp[k][:, 0], comp[k][:, 0], comp[k][:, 1])
                             for k, _ in TENURE_GROUPS}
        out["cut_usd"][key] = boot.wmean(cuts)
    return out


def part_a(part: pd.DataFrame, g: pd.DataFrame, boot: Boot, indir: Path) -> dict[str, Any]:
    """(a) percentile gaps of the housing-cost distribution, all units / renters / owners."""
    import housing_costs as hc
    czones = g.czone.to_numpy()
    blocks = distribution_blocks(part, czones)
    T = {top: percentile_tables(blocks, top) for top in ("pareto", "nan", "lower")}
    res: dict[str, Any] = {"percentiles_tabulated": list(PCT_KEYS), "tables": {}}

    for tab, b in blocks.items():
        P, F = T["pareto"][tab]["P"], T["pareto"][tab]["F"]
        d: dict[str, Any] = dict(label=b["label"], edges=[float(e) if np.isfinite(e) else None
                                                           for e in b["edges"]], percentiles={},
                                 diff_vs_median={}, sensitivity={}, cheap_share={})
        log_gap = np.log(P[:, 1, :] / P[:, 0, :])
        for i, key in enumerate(PCT_KEYS):
            st = boot.wmean(log_gap[:, i], P[:, 0, i], P[:, 1, i])
            su = boot.wmean(P[:, 1, i] - P[:, 0, i])
            d["percentiles"][key] = dict(
                log_gap=st, usd_gap=su, level_1=st["level_1"], level_2=st["level_2"],
                n_cells_top_bin=[int(F[:, 0, i].sum()), int(F[:, 1, i].sum())])
        for i, key in enumerate(PCT_KEYS):
            if key == "p50":
                continue
            d["diff_vs_median"][key] = boot.wmean(log_gap[:, i] - log_gap[:, 2])
        # top-bin handling and a common sample without any extrapolated percentile
        for top in ("nan", "lower"):
            Pt = T[top][tab]["P"]
            lg = np.log(Pt[:, 1, :] / Pt[:, 0, :])
            d["sensitivity"][f"top_{top}"] = {k: boot.wmean(lg[:, i]) for i, k in enumerate(PCT_KEYS)}
        ok = (F.sum(axis=(1, 2)) == 0)
        lg = np.where(ok[:, None], log_gap, np.nan)
        d["sensitivity"]["common_no_extrapolation"] = {k: boot.wmean(lg[:, i])
                                                        for i, k in enumerate(PCT_KEYS)}
        d["sensitivity"]["common_no_extrapolation_n_cz"] = int(ok.sum())
        # share of each location's units below the CZ-pooled 10th / 25th / 50th percentile
        pooled = np.nansum(b["counts"], axis=1)
        for q in (0.10, 0.25, 0.50):
            s1, s2 = np.full(len(czones), np.nan), np.full(len(czones), np.nan)
            for ci in range(len(czones)):
                cut, _ = pct_from_bins(pooled[ci], b["edges"], [q], "pareto")
                s1[ci] = cdf_at(b["counts"][ci, 0], b["edges"], cut[0])
                s2[ci] = cdf_at(b["counts"][ci, 1], b["edges"], cut[0])
            st = boot.wmean(s2 - s1, s1, s2)
            d["cheap_share"][f"q{int(100 * q)}"] = st
        res["tables"][tab] = d
    res["tail_composition"] = tail_composition(blocks, czones, boot)

    # reproduction of housing_costs.md and of the earnings target of district_settings.py
    gh = hc.add_pooled(part, g.copy())
    mine = res["tables"]["all"]["percentiles"]["p50"]["log_gap"]
    ref_json = indir / "estimates" / "settings" / "housing_costs.json"
    rep: dict[str, Any] = {"value_here": mine["value"], "boot_se_here": mine["boot_se"],
                           "max_abs_cz_diff_vs_housing_costs": float(np.nanmax(np.abs(
                               gh["gap_pm_h_log"].to_numpy(float)
                               - np.log(T["nan"]["all"]["P"][:, 1, 2] / T["nan"]["all"]["P"][:, 0, 2]))))}
    if ref_json.exists():
        pub = json.loads(ref_json.read_text())["partitions"]["locale"]["statistics"]["pm_h_log"]
        rep.update(value_published=pub["value"], boot_se_published=pub["boot_se"],
                   abs_diff_value=abs(mine["value"] - pub["value"]),
                   abs_diff_se=abs(mine["boot_se"] - pub["boot_se"]), n_published=pub["n"])
    if "gap_pm_ftyr_all" in gh:
        st = boot.wmean(gh["gap_pm_ftyr_all"].to_numpy(float))
        rep.update(earnings_gap_here=st["value"], earnings_gap_se_here=st["boot_se"],
                   earnings_gap_reported=0.1122, earnings_gap_se_reported=0.0168, earnings_stat=st)
    res["reproduction"] = rep
    res["bin_check"] = dict(
        total_units_B25104_bins=float(np.nansum(part[COLS_ALL].to_numpy(float))),
        renters_B25063_bins=float(np.nansum(part[COLS_RENT].to_numpy(float))),
        owners_B25087_total=float(np.nansum(part["B25087_001E"].to_numpy(float))),
        owner_open_bin_allocation=dict(edges=blocks["own"]["alloc_edges"],
                                       shares=blocks["own"]["note_alloc"],
                                       basis="with-mortgage owners (B25087) of the partitioned districts, all 188 CZs "
                                             "and both locations pooled"))
    res["n_cz"] = int(len(g))
    return res


# --------------------------------------------------------------------------
# (b) Housing expenditure shares by income: BLS CE 2018 and the Stone-Geary fit
# --------------------------------------------------------------------------

CE_STATS = ("Mean", "Share", "SE", "CV(%)")


def parse_ce(path: Path) -> dict[str, Any]:
    """Parse a BLS CE 'mean-item-share-average-standard-error' workbook (Tables 1101 / 1110).

    Each item is a label row followed by Mean / Share / SE / CV(%) rows; characteristics
    (people, earners, vehicles, tenure) are single rows of values.  Only the FIRST occurrence
    of an item label is kept (labels such as 'Housing' recur under other headings later in the
    table).  Returns column names, consumer units (thousands), lower limits, characteristics
    and items as lists over the columns (all consumer units first)."""
    import openpyxl
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wb = openpyxl.load_workbook(path, data_only=True)
    rows = [list(r) for r in wb.worksheets[0].iter_rows(values_only=True)]
    title = str(rows[0][0])
    cols = [str(c).replace("\n", " ") for c in rows[2][1:] if c is not None]
    k = len(cols)
    num = lambda r: [(float(v) if isinstance(v, (int, float)) else None) for v in r[1:1 + k]]
    out: dict[str, Any] = dict(title=title, columns=cols, n_cu_thousands=num(rows[3]),
                               lower_limit=num(rows[4]), chars={}, items={})
    cur: str | None = None
    for r in rows[5:]:
        lab = r[0]
        if lab is None or str(lab).strip() == "":
            continue
        lab = str(lab).strip()
        vals = r[1:1 + k]
        if lab in CE_STATS:
            if cur is not None and lab not in out["items"][cur]:
                out["items"][cur][lab] = num(r)
        elif all(v is None for v in vals):
            cur = lab
            out["items"].setdefault(cur, {})
        elif cur is None or lab in ("Vehicles", "People", "Earners", "Homeowner", "Renter",
                                    "Age of reference person", "Children under 18", "College",
                                    "At least one vehicle owned or leased", "With mortgage",
                                    "Without mortgage", "Adults 65 and older"):
            out["chars"].setdefault(lab, num(r))
    # 'Income before taxes' / 'Income after taxes' appear as items with Mean/SE/CV
    return out


def ce_item(ce: dict[str, Any], name: str, stat: str = "Mean") -> np.ndarray:
    x = ce["items"][name][stat]
    return np.array([np.nan if v is None else v for v in x], float)


HOUSING_DEFS = ("shelter", "utilities", "shelter_utilities", "shelter_utilities_ex_phone", "housing")
DEF_LABELS = {"shelter": "shelter", "utilities": "utilities", "shelter_utilities": "shelter + utilities",
              "shelter_utilities_ex_phone": "shelter + utilities ex. telephone", "housing": "housing"}


def housing_defs(ce: dict[str, Any]) -> dict[str, dict[str, np.ndarray]]:
    """Mean and SE of each housing definition by column: shelter, utilities, shelter +
    utilities (SE assuming independent sampling error of the two items: an upper bound),
    shelter + utilities without telephone services (shelter, natural gas, electricity, fuel oil
    and other fuels, water and other public services: the items an ACS gross rent or owner cost
    covers) and the CE 'Housing' group."""
    sh, ut, ho = (ce_item(ce, "Shelter"), ce_item(ce, "Utilities, fuels, and public services"),
                  ce_item(ce, "Housing"))
    se_sh, se_ut, se_ho = (ce_item(ce, "Shelter", "SE"),
                           ce_item(ce, "Utilities, fuels, and public services", "SE"),
                           ce_item(ce, "Housing", "SE"))
    parts = ["Natural gas", "Electricity", "Fuel oil and other fuels", "Water and other public services"]
    ex = sum(ce_item(ce, nm) for nm in parts)
    se_ex = np.sqrt(sum(ce_item(ce, nm, "SE") ** 2 for nm in parts))
    return {"shelter": dict(mean=sh, se=se_sh),
            "utilities": dict(mean=ut, se=se_ut),
            "shelter_utilities": dict(mean=sh + ut, se=np.sqrt(se_sh ** 2 + se_ut ** 2)),
            "shelter_utilities_ex_phone": dict(mean=sh + ex, se=np.sqrt(se_sh ** 2 + se_ex ** 2)),
            "housing": dict(mean=ho, se=se_ho)}


def ols_ab(Y: np.ndarray, s: np.ndarray) -> tuple[float, float]:
    """OLS of the share on 1/Y: s = a + b / Y.  With a = aH and b = (1 - aH) h this is the
    Stone-Geary share s(Y) = h/Y + aH (1 - h/Y) written linearly, so its least-squares solution
    is the nonlinear least-squares fit in (h, aH)."""
    X = np.column_stack([np.ones_like(Y), 1.0 / Y])
    coef, *_ = np.linalg.lstsq(X, s, rcond=None)
    return float(coef[0]), float(coef[1])


def sg_share(Y: np.ndarray, h: float, aH: float) -> np.ndarray:
    return h / Y + aH * (1.0 - h / Y)


def nls_sg(Y: np.ndarray, s: np.ndarray, h0: float, a0: float, iters: int = 200) -> tuple[float, float]:
    """Levenberg-Marquardt for min sum (s - sg_share(Y; h, aH))^2, in (log h, aH), from (h0, a0).
    Used to confirm that the generic nonlinear fit reproduces the linear-in-(a, b) solution."""
    th = np.array([math.log(h0), a0])
    lam = 1e-3

    def resid(t):
        return s - sg_share(Y, math.exp(t[0]), t[1])

    r = resid(th)
    for _ in range(iters):
        h, a = math.exp(th[0]), th[1]
        J = np.column_stack([-(1 - a) * h / Y, -(1 - h / Y)])           # d resid / d(log h, a)
        step = np.linalg.solve(J.T @ J + lam * np.diag(np.diag(J.T @ J)), -J.T @ r)
        tn = th + step
        rn = resid(tn)
        if rn @ rn < r @ r:
            th, r, lam = tn, rn, lam / 3
            if np.max(np.abs(step)) < 1e-13:
                break
        else:
            lam *= 4
    return math.exp(th[0]), float(th[1])


def fit_stone_geary(Y: np.ndarray, M: np.ndarray, se_Y: np.ndarray, se_M: np.ndarray,
                    ybar: float, rng: np.random.Generator, ndraw: int = 4000,
                    rho: float = 0.0) -> dict[str, Any]:
    """Stone-Geary fit of the group shares s_q = M_q / Y_q (M = item mean, Y = total expenditures
    per consumer unit) on Y by nonlinear least squares, equal weight on the groups.

    Reports aH, h (dollars) and k = h / ybar (subsistence housing expenditure as a share of mean
    expenditures), the expenditure elasticity of housing at the mean, the fit and standard
    errors from (i) the residual variance (delta method, n - 2 degrees of freedom) and (ii) a
    parametric bootstrap of the CE sampling error: item mean and total drawn normal around
    their published values with the published SEs and correlation `rho`, independently across
    groups."""
    s = M / Y
    a, b = ols_ab(Y, s)
    h = b / (1.0 - a)
    ha, aa = nls_sg(Y, s, max(h, 1.0), a)
    fit = sg_share(Y, h, a)
    res = s - fit
    n = len(Y)
    X = np.column_stack([np.ones_like(Y), 1.0 / Y])
    dof = n - 2
    s2 = float(res @ res) / dof if dof > 0 else np.nan
    cov_ab = s2 * np.linalg.inv(X.T @ X)
    # delta method for (aH, k = h / ybar), h = b / (1 - a)
    J = np.array([[1.0, 0.0], [b / (1 - a) ** 2 / ybar, 1.0 / ((1 - a) * ybar)]])
    cov_ak = J @ cov_ab @ J.T
    s_mean = a + (1 - a) * h / ybar
    out: dict[str, Any] = dict(
        alphaH=a, h_usd=h, h_share_of_mean_exp=h / ybar, mean_exp_used=ybar,
        s_at_mean_exp=s_mean, elasticity_at_mean=a / s_mean,
        nls_check_max_abs_diff=float(max(abs(ha - h) / h, abs(aa - a))),
        rmse_pp=100 * math.sqrt(float(res @ res) / n),
        max_abs_resid_pp=100 * float(np.max(np.abs(res))),
        r2=1.0 - float(res @ res) / float(((s - s.mean()) ** 2).sum()),
        n_groups=n, dof=dof,
        se_resid=dict(alphaH=float(math.sqrt(cov_ak[0, 0])), h_share=float(math.sqrt(cov_ak[1, 1])),
                      corr=float(cov_ak[0, 1] / math.sqrt(cov_ak[0, 0] * cov_ak[1, 1]))),
        shares=s.tolist(), fitted=fit.tolist(), resid_pp=(100 * res).tolist(),
        Y=Y.tolist(), subsistence_share_by_group=(h / Y).tolist())
    # parametric bootstrap
    ab = np.empty(ndraw), np.empty(ndraw), np.empty(ndraw)
    for d in range(ndraw):
        z1, z2 = rng.standard_normal(n), rng.standard_normal(n)
        Md = M + se_M * z1
        Yd = Y + se_Y * (rho * z1 + math.sqrt(1 - rho ** 2) * z2)
        ad, bd = ols_ab(Yd, Md / Yd)
        ab[0][d], ab[1][d] = ad, bd / (1 - ad) / ybar
        ab[2][d] = ad / (ad + (1 - ad) * bd / (1 - ad) / ybar)
    out["se_ce"] = dict(alphaH=float(np.std(ab[0], ddof=1)), h_share=float(np.std(ab[1], ddof=1)),
                        elasticity=float(np.std(ab[2], ddof=1)),
                        alphaH_ci95=[float(v) for v in np.percentile(ab[0], [2.5, 97.5])],
                        h_share_ci95=[float(v) for v in np.percentile(ab[1], [2.5, 97.5])],
                        rho=rho, draws=ndraw)
    return out


def tercile_weights(n_groups: int) -> np.ndarray:
    """(3 x n_groups) fractions of each equal-sized income group (decile or quintile) that fall
    in the bottom, middle and top tercile of the income distribution."""
    edges = np.linspace(0, 1, n_groups + 1)
    W = np.zeros((3, n_groups))
    for t in range(3):
        lo, hi = t / 3.0, (t + 1) / 3.0
        for g in range(n_groups):
            W[t, g] = max(0.0, min(hi, edges[g + 1]) - max(lo, edges[g])) / (edges[g + 1] - edges[g])
    return W


def ce_block(path_q: Path, path_d: Path, rng: np.random.Generator) -> dict[str, Any]:
    """CE quintile and decile tables, shares, Stone-Geary fits, tercile aggregates."""
    res: dict[str, Any] = {}
    for tag, path in (("quintiles", path_q), ("deciles", path_d)):
        ce = parse_ce(path)
        cols = ce["columns"]
        ng = len(cols) - 1
        n_cu = np.array(ce["n_cu_thousands"][1:], float)
        Y = ce_item(ce, "Average annual expenditures")[1:]
        seY = ce_item(ce, "Average annual expenditures", "SE")[1:]
        Yall = float(ce_item(ce, "Average annual expenditures")[0])
        inc = ce_item(ce, "Income before taxes")[1:] if "Income before taxes" in ce["items"] else np.full(ng, np.nan)
        defs = housing_defs(ce)
        pens = ce_item(ce, "Personal insurance and pensions")[1:]
        block: dict[str, Any] = dict(
            title=ce["title"], columns=cols[1:], n_cu_thousands=n_cu.tolist(),
            lower_limit=ce["lower_limit"][1:], mean_income_before_taxes=inc.tolist(),
            mean_expenditures=Y.tolist(), se_expenditures=seY.tolist(),
            all_cu=dict(mean_expenditures=Yall, n_cu_thousands=ce["n_cu_thousands"][0]),
            chars={k: v for k, v in ce["chars"].items()}, defs={}, fits={})
        def share_se(Mv: np.ndarray, sM: np.ndarray, Yv: np.ndarray, sY: np.ndarray, rho: float) -> np.ndarray:
            """Delta-method SE of M / Y with correlation rho between the sampling errors of M and Y."""
            cm, cy = sM / Mv, sY / Yv
            return (Mv / Yv) * np.sqrt(np.maximum(cm ** 2 + cy ** 2 - 2 * rho * cm * cy, 0.0))

        for name, d in defs.items():
            M = d["mean"][1:]
            pub = {"shelter": "Shelter", "utilities": "Utilities, fuels, and public services", "housing": "Housing"}.get(name)
            pubs = None if pub is None else [None if v is None else v / 100 for v in ce["items"][pub]["Share"][1:]]
            block["defs"][name] = dict(
                mean=M.tolist(), se=d["se"][1:].tolist(), share=(M / Y).tolist(),
                share_se_rho0=share_se(M, d["se"][1:], Y, seY, 0.0).tolist(),
                share_se_rho05=share_se(M, d["se"][1:], Y, seY, 0.5).tolist(),
                all_cu_share_se_rho0=float(share_se(d["mean"][:1], d["se"][:1], np.array([Yall]),
                                                    np.array([ce_item(ce, "Average annual expenditures", "SE")[0]]), 0.0)[0]),
                all_cu_share_se_rho05=float(share_se(d["mean"][:1], d["se"][:1], np.array([Yall]),
                                                     np.array([ce_item(ce, "Average annual expenditures", "SE")[0]]), 0.5)[0]),
                max_abs_diff_vs_published_pp=None if pubs is None else float(max(
                    abs(100 * (M[i] / Y[i]) - 100 * pubs[i]) for i in range(ng) if pubs[i] is not None)),
                share_published=pubs,
                all_cu_mean=float(d["mean"][0]), all_cu_share=float(d["mean"][0] / Yall),
                mean_of_group_shares=float(np.average(M / Y, weights=n_cu)))
        # fits (equal weight on the income groups); k = h / Ybar with Ybar = CE mean expenditures
        for name, d in defs.items():
            M, seM = d["mean"][1:], d["se"][1:]
            for rho_tag, rho in (("rho0", 0.0), ("rho05", 0.5)):
                f = fit_stone_geary(Y, M, seY, seM, CE_TOTAL_EXP, rng, rho=rho)
                if rho_tag == "rho0":
                    block["fits"][name] = f
                else:
                    block["fits"][name]["se_ce_rho05"] = f["se_ce"]
        # robustness fits: drop the lowest income group; deflate by expenditures net of personal
        # insurance and pensions (saving-like items)
        block["fits_dropbottom"], block["fits_net_of_pensions"] = {}, {}
        Yc = Y - pens
        Yc_all = CE_TOTAL_EXP - float(ce_item(ce, "Personal insurance and pensions")[0])
        for name, d in defs.items():
            M, seM = d["mean"][1:], d["se"][1:]
            block["fits_dropbottom"][name] = fit_stone_geary(Y[1:], M[1:], seY[1:], seM[1:],
                                                             CE_TOTAL_EXP, rng, ndraw=1000)
            block["fits_net_of_pensions"][name] = fit_stone_geary(
                Yc, M, np.sqrt(seY ** 2), seM, Yc_all, rng, ndraw=1000)
        block["net_pensions_mean_exp"] = Yc_all
        # terciles of income before taxes from the groups, with fractional weights
        W = tercile_weights(ng)
        block["terciles"] = {}
        for name, d in defs.items():
            M = d["mean"][1:]
            agg = [float((W[t] * n_cu * M).sum() / (W[t] * n_cu * Y).sum()) for t in range(3)]
            mor = [float((W[t] * n_cu * M / Y).sum() / (W[t] * n_cu).sum()) for t in range(3)]
            sM = d["se"][1:]
            se_t = {}
            for rho_tag, rho in (("rho0", 0.0), ("rho05", 0.5)):
                dr = np.empty((2000, 3))
                for k_ in range(2000):
                    z1, z2 = rng.standard_normal(ng), rng.standard_normal(ng)
                    Md = M + sM * z1
                    Yd = Y + seY * (rho * z1 + math.sqrt(1 - rho ** 2) * z2)
                    dr[k_] = [(W[t] * n_cu * Md).sum() / (W[t] * n_cu * Yd).sum() for t in range(3)]
                se_t[rho_tag] = np.std(dr, axis=0, ddof=1).tolist()
            block["terciles"][name] = dict(ratio_of_means=agg, mean_of_group_shares=mor, se_rho0=se_t["rho0"], se_rho05=se_t["rho05"])
        block["terciles"]["mean_expenditures"] = [float((W[t] * n_cu * Y).sum() / (W[t] * n_cu).sum())
                                                  for t in range(3)]
        block["terciles"]["mean_income_before_taxes"] = [
            float((W[t] * n_cu * inc).sum() / (W[t] * n_cu).sum()) for t in range(3)] if np.isfinite(inc).all() else None
        # vehicles, tenure and transport lines per group
        block["transport"] = {name: ce_item(ce, name)[1:].tolist() for name in
                              ("Transportation", "Vehicle purchases (net outlay)",
                               "Gasoline, other fuels, and motor oil", "Other vehicle expenses",
                               "Vehicle finance charges", "Maintenance and repairs",
                               "Vehicle insurance", "Vehicle rental, leases, licenses, and other charges",
                               "Public and other transportation")}
        block["transport_all"] = {name: float(ce_item(ce, name)[0]) for name in block["transport"]}
        block["personal_insurance_pensions"] = pens.tolist()
        res[tag] = block
    return res


# --------------------------------------------------------------------------
# Household size and the income gradient: Stone-Geary fits on equivalized expenditure
# --------------------------------------------------------------------------
#
# The CE counts consumer units (CUs).  Household size rises with income across the CE groups
# (1.7 to 3.2 persons per CU from the lowest to the highest quintile), so the fall of the housing
# share with expenditure mixes nonhomotheticity with economies of scale in shared housing.  The
# model's household is one adult.  With an equivalence scale e = n^theta (n persons) the share of a
# household is s = h1 / y + aH (1 - h1 / y) with y = Y / e the equivalized expenditure and h1 the
# subsistence housing expenditure per equivalent adult.  h1 / (mean expenditure per equivalent adult)
# is the subsistence share of a household of ANY size when housing and consumption scale together
# with household size (the assumption of housing_costs.md), i.e. the model's hmin in units of Cbar.
# theta = 0 is the CU-level fit (no adjustment), theta = 1 per capita (full sharing of nothing).
# Scales: square root of persons (the OECD income-distribution convention, equivalence elasticity
# 0.5, central), OECD-modified (1 + 0.5 (adults - 1) + 0.3 children; elasticity about 0.53 for typical
# households; CE 'children' are under 18 where the OECD's are under 14), per adult (children count
# zero: the model's literal 'one household per adult') and per capita.  Source for the scales:
# Foerster and Mira d'Ercole (OECD), 'The OECD approach to measuring income distribution and
# poverty', Table 2 (elasticities: per capita 1, OECD-modified 0.53, square root 0.50, household 0).

EQ_SCALES = ("none", "sqrt", "oecd_modified", "per_adult", "per_capita")
EQ_LABELS = {"none": "none (consumer unit)", "sqrt": "square root of persons",
             "oecd_modified": "OECD-modified", "per_adult": "per adult (children count zero)",
             "per_capita": "per capita"}
CENTRAL_SCALE = "sqrt"


def equivalence_scales(ce: dict[str, Any]) -> dict[str, Any]:
    """Persons, children under 18, adults (persons minus children) and earners per CU of each income
    group (CE characteristics rows, one decimal) and the equivalence scales built from them."""
    people = np.array(ce["chars"]["People"][1:], float)
    kids = np.array(ce["chars"]["Children under 18"][1:], float)
    adults = people - kids
    earners = np.array(ce["chars"]["Earners"][1:], float)
    scales = {"none": np.ones_like(people), "sqrt": np.sqrt(people),
              "oecd_modified": 1.0 + 0.5 * (adults - 1.0) + 0.3 * kids,
              "per_adult": adults, "per_capita": people}
    return dict(people=people, children=kids, adults=adults, earners=earners, scales=scales)


def theta_profile(Y: np.ndarray, M: np.ndarray, people: np.ndarray, n_cu: np.ndarray) -> dict[str, Any]:
    """Stone-Geary fit of the shares M / Y on equivalized expenditure Y / n^theta over a grid of theta
    (persons as the size measure): aH, h over mean expenditure per equivalent adult and the RMSE.
    The theta with the smallest RMSE on a 0.01 grid is reported (two degrees of freedom: suggestive)."""
    s = M / Y

    def one(th: float) -> dict[str, float]:
        e = people ** th
        ybar = CE_TOTAL_EXP / float((n_cu * e).sum() / n_cu.sum())
        a, b = ols_ab(Y / e, s)
        h = b / (1.0 - a)
        res = s - sg_share(Y / e, h, a)
        return dict(theta=float(th), alphaH=a, h_share=h / ybar, rmse_pp=100 * math.sqrt(float(res @ res) / len(s)))

    rows = [one(th) for th in np.round(np.arange(0.0, 1.0001, 0.05), 2)]
    best = min((one(th) for th in np.round(np.arange(0.0, 1.0001, 0.01), 2)), key=lambda r: r["rmse_pp"])
    return dict(rows=rows, best=best)


def ce_equivalized(path_q: Path, path_d: Path, rng: np.random.Generator) -> dict[str, Any]:
    """Stone-Geary fits of every housing definition on equivalized expenditure for each scale but
    'none' (the CU-level fits are in `ce_block`), for the quintile and decile tables, with the lowest
    group dropped as a robustness check.  y = Y / e; h is in dollars per equivalent adult and h / ybar
    uses ybar = CE mean expenditures per CU / mean scale (aggregate expenditure per equivalent adult)."""
    res: dict[str, Any] = {}
    for tag, path in (("quintiles", path_q), ("deciles", path_d)):
        ce = parse_ce(path)
        ng = len(ce["columns"]) - 1
        n_cu = np.array(ce["n_cu_thousands"][1:], float)
        Y = ce_item(ce, "Average annual expenditures")[1:]
        seY = ce_item(ce, "Average annual expenditures", "SE")[1:]
        defs = housing_defs(ce)
        hh = equivalence_scales(ce)
        blk: dict[str, Any] = dict(
            people=hh["people"].tolist(), children=hh["children"].tolist(), adults=hh["adults"].tolist(),
            earners=hh["earners"].tolist(), scales={}, fits={}, fits_dropbottom={})
        for sname in EQ_SCALES:
            e = hh["scales"][sname]
            ebar = float((n_cu * e).sum() / n_cu.sum())
            blk["scales"][sname] = dict(values=e.tolist(), mean_scale=ebar, mean_exp_per_equivalent_adult=CE_TOTAL_EXP / ebar,
                                        mean_of_group_equivalized_exp=float((n_cu * Y / e).sum() / n_cu.sum()))
            if sname == "none":
                continue
            ybar = CE_TOTAL_EXP / ebar
            blk["fits"][sname], blk["fits_dropbottom"][sname] = {}, {}
            for name, d in defs.items():
                M, seM = d["mean"][1:], d["se"][1:]
                blk["fits"][sname][name] = fit_stone_geary(Y / e, M / e, seY / e, seM / e, ybar, rng)
                blk["fits_dropbottom"][sname][name] = fit_stone_geary(
                    (Y / e)[1:], (M / e)[1:], (seY / e)[1:], (seM / e)[1:], ybar, rng, ndraw=1000)
        blk["theta_profile"] = theta_profile(Y, defs["shelter_utilities"]["mean"][1:], hh["people"], n_cu)
        res[tag] = blk
    return res


def ce_size_table(path_size: Path) -> dict[str, Any]:
    """CE Table 1400 (size of consumer unit): all CUs, one person, two or more people.  Expenditures,
    shelter + utilities shares, persons, vehicles per CU and the vehicle spending lines."""
    ce = parse_ce(path_size)
    cols = ce["columns"]
    item = lambda nm, st="Mean": ce_item(ce, nm, st)
    Y, seY = item("Average annual expenditures"), item("Average annual expenditures", "SE")
    sh, seSh = item("Shelter"), item("Shelter", "SE")
    ut, seUt = item("Utilities, fuels, and public services"), item("Utilities, fuels, and public services", "SE")
    inc = item("Income before taxes")
    people = np.array(ce["chars"]["People"], float)
    veh = np.array(ce["chars"]["Vehicles"], float)
    out: dict[str, Any] = dict(title=ce["title"], columns=cols, groups={})
    for j, nm in enumerate(cols):
        M = sh[j] + ut[j]
        seM = math.sqrt(seSh[j] ** 2 + seUt[j] ** 2)
        cm, cy = seM / M, seY[j] / Y[j]
        out["groups"][nm] = dict(
            people=float(people[j]), vehicles=float(veh[j]), expenditures=float(Y[j]), se_expenditures=float(seY[j]),
            income_before_taxes=float(inc[j]), shelter=float(sh[j]), utilities=float(ut[j]), shelter_utilities=float(M),
            share=float(M / Y[j]), share_se_rho0=float(M / Y[j] * math.sqrt(cm ** 2 + cy ** 2)),
            private_vehicle_spending=float(item("Vehicle purchases (net outlay)")[j] + item("Gasoline, other fuels, and motor oil")[j]
                                           + item("Other vehicle expenses")[j]))
    return out


def size_check(size: dict[str, Any], fit_cu: dict[str, Any], fit_eq: dict[str, Any]) -> dict[str, Any]:
    """Compare the shelter + utilities shares of one-person and two-or-more-person CUs (CE Table 1400)
    with the Stone-Geary curves of the income quintiles (CU-level fit and the central equivalized fit),
    and give the subsistence share implied for each group by the central marginal share aH."""
    out: dict[str, Any] = {}
    a_eq, h_eq = fit_eq["alphaH"], fit_eq["h_usd"]
    for nm, grp in size["groups"].items():
        if nm == size["columns"][0]:
            continue
        n = grp["people"]
        Y = grp["expenditures"]
        s = grp["share"]
        pred_cu = float(sg_share(np.array([Y]), fit_cu["h_usd"], fit_cu["alphaH"])[0])
        e = math.sqrt(n)
        pred_eq = float(sg_share(np.array([Y / e]), h_eq, a_eq)[0])
        out[nm] = dict(people=n, expenditures=Y, share=s, share_se=grp["share_se_rho0"], predicted_cu_level=pred_cu,
                       predicted_sqrt_scale=pred_eq, excess_over_cu_level_pp=100 * (s - pred_cu),
                       excess_over_sqrt_scale_pp=100 * (s - pred_eq),
                       implied_subsistence_share_given_aH_sqrt=(s - a_eq) / (1.0 - a_eq),
                       implied_subsistence_usd=(s - a_eq) / (1.0 - a_eq) * Y)
    return out


# ACS national housing-cost burdens by household income (B25074 renters, B25095 owners, B25106)
BURDEN_EDGES = [0.0, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, np.inf]
INCOME_BINS_7 = [(0, 10000), (10000, 20000), (20000, 35000), (35000, 50000), (50000, 75000),
                 (75000, 100000), (100000, np.inf)]


def load_national(prod: str, tbl: str, suffix: str = "E") -> dict[str, float]:
    """Estimates (suffix 'E') or margins of error ('M', 90%) of a national ACS table, keyed by
    variable name with the suffix removed."""
    path = RAW / "acs_national" / f"{prod}_{YEAR}_{tbl}.json"
    if not path.exists():
        raise SystemExit(f"no {path.name}: run `housing_transport.py fetch` first")
    d = json.loads(path.read_text())
    return {k[:-1] + "E": float(v) for k, v in d.items() if re.fullmatch(r"B\d+_\d+" + suffix, k) and v is not None}


def acs_mean_income_by_bin(prod: str) -> list[float]:
    """Mean household income of each of the 7 income bins, from the B19001 distribution (bin
    midpoints) and, for the open $100,000+ bin, from aggregate household income B19025_001."""
    b = load_national(prod, "B19001")
    agg = load_national(prod, "B19025")["B19025_001E"]
    edges = [0, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 60000, 75000, 100000,
             125000, 150000, 200000]
    cnt = [b[f"B19001_{j:03d}E"] for j in range(2, 18)]
    mids = [(edges[i] + edges[i + 1]) / 2 for i in range(15)]
    lo = sum(c * m for c, m in zip(cnt[:12], mids[:12]))
    n_top = sum(cnt[12:])
    top_mean = (agg - lo) / n_top
    grp = [(0, 1), (1, 3), (3, 6), (6, 9), (9, 11), (11, 12)]
    out = []
    for a_, b_ in grp:
        n = sum(cnt[a_:b_])
        out.append(sum(c * m for c, m in zip(cnt[a_:b_], mids[a_:b_])) / n)
    out.append(top_mean)
    return out


def acs_burden(prod: str, lt20: float = 0.12, top50: float = 0.70) -> dict[str, Any]:
    """Distribution of housing cost / household income by income bin for all households
    (renters: gross rent as a percentage of income, B25074; owners: selected monthly owner
    costs, B25095), households with a computed burden.  Mean burden uses class midpoints, with
    `lt20` for the class below 20% and `top50` for the open class of 50% and over; the median
    is read by linear interpolation inside classes (flagged when it lies in an open class)."""
    r, o = load_national(prod, "B25074"), load_national(prod, "B25095")
    inc = acs_mean_income_by_bin(prod)
    rows = []
    for i in range(7):
        cr = [r[f"B25074_{2 + 9 * i + 1 + c:03d}E"] for c in range(7)]
        if i < 6:
            co = [o[f"B25095_{2 + 9 * i + 1 + c:03d}E"] for c in range(7)]
        else:
            co = [o[f"B25095_{2 + 9 * 6 + 1 + c:03d}E"] + o[f"B25095_{2 + 9 * 7 + 1 + c:03d}E"]
                  for c in range(7)]
        cnt = np.array(cr, float) + np.array(co, float)
        tot = cnt.sum()
        mids = np.array([lt20, 0.225, 0.275, 0.325, 0.375, 0.45, top50])
        cum = np.cumsum(cnt) / tot
        j = int(np.searchsorted(cum, 0.5))
        lo_e, hi_e = BURDEN_EDGES[j], BURDEN_EDGES[j + 1]
        below = cum[j - 1] if j > 0 else 0.0
        med_open = (j == 0) or (j == 6)
        med = lo_e + (min(hi_e, 1.0) - lo_e) * (0.5 - below) / (cum[j] - below) if j not in (6,) else np.nan
        rows.append(dict(
            income_bin=[INCOME_BINS_7[i][0], None if np.isinf(INCOME_BINS_7[i][1]) else INCOME_BINS_7[i][1]],
            mean_income=inc[i], households=float(tot), households_renters=float(sum(cr)),
            households_owners=float(sum(co)), share_30plus=float(cnt[3:].sum() / tot),
            share_50plus=float(cnt[6] / tot), share_below_20=float(cnt[0] / tot),
            mean_burden=float((cnt * mids).sum() / tot), median_burden=float(med),
            median_in_open_class=bool(med_open)))
    # B25106: coarse income bins, share with costs of 30 percent or more of income
    b6 = load_national(prod, "B25106")
    m6 = load_national(prod, "B25106", "M")
    coarse = []
    tot_n = tot_h = 0.0
    var_n = var_h = 0.0
    for lab, ro, rr in (("Less than $20,000", 3, 25), ("$20,000-34,999", 7, 29), ("$35,000-49,999", 11, 33),
                        ("$50,000-74,999", 15, 37), ("$75,000 or more", 19, 41)):
        n_all = sum(b6[f"B25106_{base + 0:03d}E"] for base in (ro, rr))
        heavy = sum(b6[f"B25106_{base + 3:03d}E"] for base in (ro, rr))
        moe_n = math.sqrt(sum(m6[f"B25106_{base + 0:03d}E"] ** 2 for base in (ro, rr)))
        moe_h = math.sqrt(sum(m6[f"B25106_{base + 3:03d}E"] ** 2 for base in (ro, rr)))
        p_ = heavy / n_all
        moe_p = math.sqrt(max(moe_h ** 2 - p_ ** 2 * moe_n ** 2, 0.0)) / n_all     # Census proportion formula
        coarse.append(dict(income_bin=lab, households=float(n_all), share_30plus=float(p_), share_30plus_moe90=float(moe_p)))
        tot_n, tot_h, var_n, var_h = tot_n + n_all, tot_h + heavy, var_n + moe_n ** 2, var_h + moe_h ** 2
    p_all = tot_h / tot_n
    all_b6 = dict(share_30plus=float(p_all), moe90=float(math.sqrt(max(var_h - p_all ** 2 * var_n, 0.0)) / tot_n),
                  households=float(tot_n))
    ten = load_national(prod, "B25003")
    total_hh = b6["B25106_001E"]
    out = dict(product=prod, rows=rows, b25106=coarse, b25106_all=all_b6, assumptions=dict(lt20=lt20, top50=top50),
               share_30plus_all=float(sum(rw["share_30plus"] * rw["households"] for rw in rows)
                                      / sum(rw["households"] for rw in rows)),
               households_total_b25106=float(total_hh),
               mean_household_income=float(load_national(prod, "B19025")["B19025_001E"]
                                           / load_national(prod, "B25003")["B25003_001E"]))
    return out


def acs_fit(burden: dict[str, Any]) -> dict[str, Any]:
    """Stone-Geary form on income: mean burden s(y) = a + b / y by weighted least squares on the
    income bins (weights = households), a = marginal share of income spent on housing, k =
    h / mean income with h = b / (1 - a)."""
    rows = burden["rows"]
    y = np.array([r["mean_income"] for r in rows])
    s = np.array([r["mean_burden"] for r in rows])
    w = np.array([r["households"] for r in rows])
    X = np.column_stack([np.ones_like(y), 1 / y])
    sw = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * sw, s * sw.ravel(), rcond=None)
    a, b = float(coef[0]), float(coef[1])
    ybar = float(np.average(y, weights=w))
    res = s - X @ coef
    return dict(alpha=a, h_usd=b / (1 - a), h_share_of_mean_income=b / (1 - a) / ybar, mean_income=ybar,
                rmse_pp=100 * float(math.sqrt(np.average(res ** 2, weights=w))),
                fitted=(X @ coef).tolist(), shares=s.tolist())


def acs_block() -> dict[str, Any]:
    out: dict[str, Any] = {}
    for prod in ("acs1", "acs5"):
        base = acs_burden(prod)
        out[prod] = dict(burden=base, fit=acs_fit(base), sensitivity={})
        for lt20, top50 in ((0.10, 0.60), (0.14, 0.70), (0.16, 0.80)):
            bs = acs_burden(prod, lt20, top50)
            out[prod]["sensitivity"][f"lt20_{lt20}_top50_{top50}"] = acs_fit(bs)
    return out


# --------------------------------------------------------------------------
# (c) Housing supply elasticities: Saiz (2010) and Baum-Snow and Han (2024)
# --------------------------------------------------------------------------

BSH_TYPES = {"units": "housing units", "space": "floorspace", "devl": "land development"}

# Horizon convention.  In a steady state of a stock with annual depreciation delta, a response measured
# over T years is a fraction 1 - (1 - delta)^T of the long-run response (Baum-Snow and Duranton, slide 76
# of 90: eps_LR = eps_T / (1 - (1 - delta)^T); with delta = 3.5% the long run is 3.33 times a 10-year and
# 1.5 times a 30-year elasticity).  Saiz (2010) uses 1970-2000 changes (30 years); Baum-Snow and Han (2024)
# and the Saiz-method replication in their Online Appendix D use 2000-2010 changes (10 years).  The SAME
# convention is applied to every anchor of `eps_recommendation`.
SAIZ_YEARS, BSH_YEARS, DELTA_CENTRAL = 30, 10, 0.035
DELTAS = (0.015, 0.02, 0.035)


def lr_factor(years: float, delta: float) -> float:
    """Long-run elasticity / T-year elasticity = 1 / (1 - (1 - delta)^T)."""
    return 1.0 / (1.0 - (1.0 - delta) ** years)


# Numbers read from the cached accepted manuscript of Baum-Snow and Han (September 5, 2023), checked against
# the PDF text by `verify_phrases` (the pages are those of the accepted-manuscript file: Online Appendix D is
# on appendix pp. 17-18 = PDF pp. 71-72, Table A8 on appendix p. 26 = PDF p. 80, the metro-region means of
# Figure 4 in the main text on p. 33 = PDF p. 34).
BSH_APPENDIX_D = dict(
    source="Baum-Snow and Han (2024), accepted manuscript of September 5, 2023: Online Appendix D, 'Comparison With "
           "Estimates in Saiz (2010)' (appendix pp. 17-18) and Table A8 'Metro Level Supply Regressions' (appendix p. 26)",
    saiz_1970_2000_mean=2.54, saiz_1970_2000_sd=1.43, saiz_n_metros=237,
    saiz_method_2000_2010_mean=1.25,
    table_a8_means_1970_2000=[1.54, 2.19, 2.56, 2.54],
    table_a8_means_2000_2010=[0.803, 0.910, 0.803, 1.270, 0.949, 1.253, 0.855, 1.300, 1.011],
    bsh_metro_region_means_units=[0.38, 0.41], bsh_metro_region_sd=0.11,
    text="'Table A8 replicates Saiz' estimation for the 2000-2010 period using the same set of metros, specifications and "
         "variables. The average implied metro supply elasticity declines to 1.25 for 2000-2010, indicating that many "
         "metropolitan areas have become increasingly supply constrained over time.' Saiz' 1970-2000 estimates have a mean "
         "of 2.54 (s.d. 1.43) for the 237 metro regions of the BSH data, against 0.38-0.41 (s.d. 0.11) for BSH's metro-region "
         "aggregates of their tract elasticities.")
BSH_PHRASES = ["declines to 1.25 for 2000-2010", "with a mean of 2.54 and a s.d. of 1.43", "with means of 0.38-0.41",
               "Table A8: Metro Level Supply Regressions", "mean unit supply elasticity across metro regions of 0.38"]
BD_PHRASES = ["Long-run elasticity is 3.33 times a 10-year elasticity and 1.5 times a 30-year elasticity",
              "Falling supply elasticities in US supply-constrained markets only",
              "Population-weighted average = 1.6", "Handbook of Regional and Urban Economics Volume 6",
              "Gorback & Keys (2024): Zip code level 2010-2020.", "Find even smaller supply elasticities"]


def verify_phrases(pdf: Path, phrases: list[str]) -> dict[str, Any]:
    """Check that each phrase occurs in the cached PDF's text (pdftotext; whitespace and line breaks
    normalised), so that a constant typed from a paper is tied to the cached file."""
    import shutil
    import subprocess
    if not pdf.exists() or shutil.which("pdftotext") is None:
        return dict(checked=False, reason="pdftotext or the cached file is missing")
    txt = " ".join(subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True).stdout.split())
    found = {p: (" ".join(p.split()) in txt) for p in phrases}
    return dict(checked=True, found=found, all_found=all(found.values()))


def wavg(x: Any, w: Any) -> float:
    """Weighted mean over the entries where x and w are finite and w > 0."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    return float(np.average(x[m], weights=w[m])) if m.any() else float("nan")


def wq(x: np.ndarray, w: np.ndarray, q: float) -> float:
    m = np.isfinite(x) & (w > 0)
    xs, ws = x[m], w[m]
    o = np.argsort(xs)
    xs, ws = xs[o], ws[o]
    c = (np.cumsum(ws) - 0.5 * ws) / ws.sum()
    return float(np.interp(q, c, xs))


def part_c(city_share: float) -> dict[str, Any]:
    """Saiz (2010) MSA elasticities and Baum-Snow and Han (2024) tract elasticities, with the
    tract file read as the authors released it (September 2023): predicted 2000-2010
    elasticities of housing units ('units'), floorspace ('space') and land development
    ('devl') by tract, from linear ('a') and quadratic-in-developed-fraction ('b') models,
    estimated by IV and by the finite mixture model (FMM-IV); `pctdis` is the fraction of the
    way from the region's CBD to its edge, `occhu0` the tract's occupied units in 2000, `wgt` the
    tract's weight (1 / number of regions it appears in)."""
    out: dict[str, Any] = {}
    # ---- Saiz (2010): author's file, 269 MSAs; Table VI lists the 95 with population > 500,000
    sz = pd.read_stata(RAW / "saiz_elasticity.dta")
    big = sz[sz.population > 500000]

    def sstats(x: pd.DataFrame) -> dict[str, Any]:
        e, p = x.elasticity.to_numpy(float), x.population.to_numpy(float)
        return dict(n=int(len(x)), mean=float(e.mean()), pop_weighted_mean=float(np.average(e, weights=p)),
                    median=float(np.median(e)), p10=float(np.percentile(e, 10)),
                    p25=float(np.percentile(e, 25)), p75=float(np.percentile(e, 75)),
                    p90=float(np.percentile(e, 90)), min=float(e.min()), max=float(e.max()))

    sel = {"San Francisco": "San Francisco, CA", "Oakland": "Oakland, CA", "San Jose": "San Jose, CA",
           "New York": "New York, NY", "Los Angeles": "Los Angeles-Long Beach, CA", "Chicago": "Chicago, IL",
           "Boston": "Boston", "Seattle": "Seattle", "Washington": "Washington, DC", "Atlanta": "Atlanta, GA",
           "Dallas": "Dallas, TX", "Houston": "Houston, TX", "Phoenix": "Phoenix-Mesa, AZ",
           "Minneapolis": "Minneapolis-St. Paul", "Indianapolis": "Indianapolis, IN"}
    named = {}
    for k, pat in sel.items():
        m = sz[sz.msaname.str.contains(pat, regex=False)]
        if len(m):
            named[k] = float(m.elasticity.iloc[0])
    out["saiz"] = dict(all_msas=sstats(sz), msas_over_500k=sstats(big), named=named,
                       source="Saiz (2010, QJE 125(3): 1253-1296), author's MSA file (269 MSAs; "
                              "elasticity, population in 2000); Table VI = the 95 MSAs with population "
                              "over 500,000 (QJE pp. 1282-1284)",
                       paper_text="population-weighted average 1.75, unweighted 2.5 (p. 1281)")
    # ---- Baum-Snow and Han (2024)
    path = RAW / "BaumSnowHan" / "gammas_hat_all.dta"
    cols = ["region_all", "pctdis", "wgt", "occhu0", "Trct_FrcDev_01", "lsf_1_flat_plains",
            "primary_sample"] + [f"gamma01a_{t}_FMM" for t in BSH_TYPES] \
        + [f"gamma01b_{t}_FMM" for t in BSH_TYPES] + [f"gamma01a_{t}_IV" for t in BSH_TYPES] \
        + [f"gamma01b_{t}_IV" for t in BSH_TYPES] + [f"sd_gamma01b_{t}_FMM" for t in ("units", "space")]
    d = pd.read_stata(path, columns=cols)
    d["w_tract"] = d.wgt.astype(float)
    d["w_hu"] = d.wgt.astype(float) * d.occhu0.astype(float)
    sp = {"FMM_quad": "gamma01b_{t}_FMM", "FMM_lin": "gamma01a_{t}_FMM", "IV_quad": "gamma01b_{t}_IV",
          "IV_lin": "gamma01a_{t}_IV"}
    bsh: dict[str, Any] = dict(
        n_rows=int(len(d)), n_regions=int(d.region_all.nunique()),
        n_tracts_weighted=float(d.w_tract.sum()), n_primary_sample=int(d.primary_sample.sum()))
    # overall: tract-weighted (paper's 'mean across urban neighborhoods') and housing-unit weighted
    overall: dict[str, Any] = {}
    for spec, pat in sp.items():
        for t in BSH_TYPES:
            x = d[pat.format(t=t)].to_numpy(float)
            overall[f"{spec}_{t}"] = dict(
                mean_tract=wavg(x, d.w_tract), mean_housing_units=wavg(x, d.w_hu),
                p25_tract=wq(x, d.w_tract.to_numpy(), 0.25), p75_tract=wq(x, d.w_tract.to_numpy(), 0.75))
    bsh["overall"] = overall
    # profile by distance to the CBD (fraction of the way to the region edge)
    bins = [0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.0001]
    d["dist_bin"] = pd.cut(d.pctdis, bins, right=False)
    prof = []
    for iv, grp in d.groupby("dist_bin", observed=True):
        rec = dict(bin=f"{iv.left:.2f}-{min(iv.right, 1.0):.2f}", share_housing_units=float(grp.w_hu.sum() / d.w_hu.sum()),
                   mean_frac_developed=wavg(grp.Trct_FrcDev_01, grp.w_hu))
        for spec in ("FMM_quad", "FMM_lin", "IV_lin"):
            for t in ("units", "space"):
                rec[f"{spec}_{t}"] = wavg(grp[sp[spec].format(t=t)], grp.w_hu)
        prof.append(rec)
    bsh["by_cbd_distance"] = prof
    # developed fraction
    dbins = [0, 0.10, 0.30, 0.60, 1.0001]
    d["dev_bin"] = pd.cut(d.Trct_FrcDev_01, dbins, right=False)
    dev = []
    for iv, grp in d.groupby("dev_bin", observed=True):
        rec = dict(bin=f"{iv.left:.2f}-{min(iv.right, 1.0):.2f}", share_housing_units=float(grp.w_hu.sum() / d.w_hu.sum()),
                   mean_cbd_distance=wavg(grp.pctdis, grp.w_hu))
        for spec in ("FMM_quad", "FMM_lin", "IV_lin"):
            for t in ("units", "space"):
                rec[f"{spec}_{t}"] = wavg(grp[sp[spec].format(t=t)], grp.w_hu)
        dev.append(rec)
    bsh["by_developed_fraction"] = dev
    # city / suburb: within each region the tracts closest to the CBD holding `share` of its
    # occupied units are the city, the rest the suburb
    d = d.sort_values(["region_all", "pctdis"]).reset_index(drop=True)
    d["hu"] = d.occhu0.astype(float) * d.wgt.astype(float)
    tot = d.groupby("region_all").hu.transform("sum")
    cum = d.groupby("region_all").hu.cumsum() - 0.5 * d.hu            # cumulative share at the tract midpoint
    d["cum_share"] = cum / tot
    splits: dict[str, Any] = {}
    for share in sorted({0.25, 0.33, round(city_share, 3), 0.50}):
        d["city"] = d.cum_share <= share
        rec: dict[str, Any] = dict(city_share_of_housing_units=share)
        for spec in sp:
            for t in BSH_TYPES:
                x = d[sp[spec].format(t=t)]
                c, s_ = wavg(x[d.city], d.hu[d.city]), wavg(x[~d.city], d.hu[~d.city])
                rec[f"{spec}_{t}"] = dict(city=c, suburb=s_, ratio=s_ / c if c > 0 else np.nan)
        # within-region ratio of suburb to city (region-level means of each side, housing-unit
        # weighted); the median over regions is robust to regions with a city mean near zero
        for spec in ("FMM_quad", "FMM_lin"):
            for t in ("units", "space"):
                col = sp[spec].format(t=t)
                dd = d[d[col].notna()]
                gm = dd.assign(vh=dd[col] * dd.hu).groupby(["region_all", "city"])[["vh", "hu"]].sum()
                gm["m"] = gm.vh / gm.hu
                wide = gm["m"].unstack("city")
                hu_r = gm["hu"].groupby("region_all").sum()
                ok = wide[True].notna() & wide[False].notna() & (wide[True] > 0.02)
                ratio = (wide[False] / wide[True])[ok]
                rec[f"{spec}_{t}_within_region_ratio"] = dict(
                    median=float(np.median(ratio)), weighted_median=wq(ratio.to_numpy(), hu_r[ok].to_numpy(), 0.5),
                    n_regions=int(ok.sum()))
        rec["mean_frac_developed"] = dict(city=wavg(d.Trct_FrcDev_01[d.city], d.hu[d.city]),
                                          suburb=wavg(d.Trct_FrcDev_01[~d.city], d.hu[~d.city]))
        splits[f"city_share_{share}"] = rec
    bsh["city_suburb"] = splits
    bsh["city_share_used"] = round(city_share, 3)
    out["baum_snow_han"] = bsh
    # long-run scaling of a T-year elasticity under annual depreciation delta (Baum-Snow and Duranton,
    # slide 76: eps_LR = eps_T / (1 - (1 - delta)^T); 3.33 for T = 10 and 1.5 for T = 30 at delta = 3.5%)
    out["longrun_factor"] = {f"delta_{int(1000 * dl)}": lr_factor(BSH_YEARS, dl) for dl in DELTAS}
    out["longrun_factor_30y"] = {f"delta_{int(1000 * dl)}": lr_factor(SAIZ_YEARS, dl) for dl in DELTAS}
    out["appendix_d"] = BSH_APPENDIX_D
    out["verification"] = dict(
        bsh_appendix=verify_phrases(RAW / "BaumSnowHan_2023_accepted.pdf", BSH_PHRASES),
        baum_snow_duranton_slides=verify_phrases(RAW / "BaumSnowDuranton_slides_2025.pdf", BD_PHRASES))
    return out


# --------------------------------------------------------------------------
# (d) Transport: vehicles by location, a within-CZ regression, and the cost of a car
# --------------------------------------------------------------------------

# AAA, "Your Driving Costs", 2018 edition: the composite "Average" vehicle (nine categories
# weighted by sales; brochure p. 9) for a new vehicle driven 15,000 miles a year.  The
# ownership ("fixed") costs do not depend on miles; the operating costs are per mile.
AAA_2018 = dict(
    source="AAA, Your Driving Costs, 2018 edition, composite 'Average' column, brochure p. 9 "
           "(https://exchange.aaa.com/wp-content/uploads/2018/09/18-0090_2018-Your-Driving-Costs-"
           "Brochure_FNL-Lo-5-2.pdf)",
    insurance=1189.0, license_registration_taxes=738.0, depreciation_15k=3289.0, finance=744.0,
    ownership_per_year=5960.0, fuel_cents_per_mile=11.05, maintenance_cents_per_mile=8.21,
    operating_cents_per_mile=19.26, total_per_year_10k=7531.0, total_per_year_15k=8849.0,
    total_per_year_20k=10213.0, cents_per_mile_15k=58.99,
    note="new vehicle; depreciation over 5 years / 75,000 miles; finance on a 5-year loan with 10% "
         "down; full-coverage insurance; 55/45 city/highway fuel economy at $2.523 per gallon")


def verify_aaa(pdf: Path) -> dict[str, Any]:
    """Check that the AAA constants above appear in the cached brochure text (pdftotext)."""
    import shutil
    import subprocess
    if not pdf.exists() or shutil.which("pdftotext") is None:
        return dict(checked=False, reason="pdftotext or the cached brochure is missing")
    txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True).stdout
    need = ["$5,960", "$8,849", "$7,531", "$10,213", "$1,189", "$738", "$3,289", "$744",
            "11.05 cents", "8.21 cents", "19.26 cents", "58.99 cents"]
    return dict(checked=True, found={n: (n in txt) for n in need}, all_found=all(n in txt for n in need))


def wls_fe_boot(y: np.ndarray, X: np.ndarray, w: np.ndarray, cz: np.ndarray, boot: Boot,
                cz_index: dict[Any, int]) -> dict[str, Any]:
    """Weighted least squares of y on X with CZ fixed effects (weighted within-transform), the
    cluster-robust SE by CZ, and a whole-CZ bootstrap SE: a CZ drawn m times carries weight m x
    w for all its districts, and the fixed effect absorbs the CZ means, so each draw solves
    (sum_c m_c S_c)^-1 (sum_c m_c s_c) with the per-CZ cross-products S_c, s_c."""
    k = X.shape[1]
    ci = np.array([cz_index[c] for c in cz])
    ncz = boot.n_cz
    sw = np.bincount(ci, weights=w, minlength=ncz)
    Xm = np.column_stack([np.bincount(ci, weights=w * X[:, j], minlength=ncz) for j in range(k)]) / sw[:, None]
    ym = np.bincount(ci, weights=w * y, minlength=ncz) / sw
    Xt, yt = X - Xm[ci], y - ym[ci]
    S = np.zeros((ncz, k, k))
    s = np.zeros((ncz, k))
    for a in range(k):
        s[:, a] = np.bincount(ci, weights=w * Xt[:, a] * yt, minlength=ncz)
        for b in range(a, k):
            S[:, a, b] = S[:, b, a] = np.bincount(ci, weights=w * Xt[:, a] * Xt[:, b], minlength=ncz)
    A = S.sum(axis=0)
    beta = np.linalg.solve(A, s.sum(axis=0))
    u = yt - Xt @ beta
    # cluster-robust (CZ) sandwich
    sc = np.zeros((ncz, k))
    for a in range(k):
        sc[:, a] = np.bincount(ci, weights=w * Xt[:, a] * u, minlength=ncz)
    present = np.bincount(ci, minlength=ncz) > 0
    G, n = int(present.sum()), len(y)
    meat = sc.T @ sc * (G / max(G - 1, 1)) * ((n - 1) / max(n - ncz - k, 1))
    Ainv = np.linalg.inv(A)
    V = Ainv @ meat @ Ainv
    # bootstrap
    cnt = boot.counts * present[None, :]
    num = np.einsum("bc,ck->bk", cnt, s)
    den = np.einsum("bc,cjk->bjk", cnt, S)
    draws = np.array([np.linalg.solve(den[b], num[b]) for b in range(len(cnt))])
    tss = float(np.sum(w * (yt - np.average(yt, weights=w)) ** 2))
    return dict(beta=beta.tolist(), se_cluster=np.sqrt(np.diag(V)).tolist(),
                se_boot=np.std(draws, axis=0, ddof=1).tolist(), n=int(n), n_cz=G,
                r2_within=1.0 - float(np.sum(w * u ** 2)) / tss, draws=draws)


# Household-size classes of B08201 (household size by vehicles available): the first variable of each size
# class is its total (1-person 007, 2-person 013, 3-person 019, 4-or-more-person 025), followed by the
# households with 0, 1, 2, 3 and 4 or more vehicles.
SIZE_BASE = {1: 7, 2: 13, 3: 19, 4: 25}
SIZE_LABELS = {1: "1 person", 2: "2 persons", 3: "3 persons", 4: "4 or more persons"}
SIZE_KS = (1, 2, 3, 4)


def b08201_col(k: int, j: int) -> str:
    """B08201 variable of size class k: j = 0 the class total, j = 1..5 the households with 0, 1, 2, 3 and
    4 or more vehicles."""
    return f"B08201_{SIZE_BASE[k] + j:03d}E"


def size_pools(part: pd.DataFrame, czones: np.ndarray) -> dict[str, Any]:
    """CZ x location pools of the B08201 cells.  Vehicles per household within a size class use the
    class midpoints 0, 1, 2, 3 and, for '4 or more', a mean m calibrated by location so that the
    class-based total equals the aggregate B25046_001 / B25044_001 pooled over all CZs (m = 4.37 in both
    locations).  Returns households N and zero-vehicle households N0 (n_cz x 2 x 4), vehicles per household
    V and zero-vehicle shares Z by class, and the class-based and B25046-based totals."""
    cols = [b08201_col(k, j) for k in SIZE_KS for j in range(6)] + ["B25046_001E", "B25044_001E"]
    P = pooled_by_location(part, cols, czones)
    ix = {c: i for i, c in enumerate(cols)}
    c_ = lambda name: P[:, :, ix[name]]
    A = np.nansum(c_("B25046_001E"), axis=0)
    low, n4 = np.zeros(2), np.zeros(2)
    for k in SIZE_KS:
        for j, mult in ((2, 1), (3, 2), (4, 3)):
            low += mult * np.nansum(c_(b08201_col(k, j)), axis=0)
        n4 += np.nansum(c_(b08201_col(k, 5)), axis=0)
    m = (A - low) / n4
    N = np.stack([c_(b08201_col(k, 0)) for k in SIZE_KS], axis=2)
    N0 = np.stack([c_(b08201_col(k, 1)) for k in SIZE_KS], axis=2)
    V = np.zeros_like(N)
    for ci, k in enumerate(SIZE_KS):
        v = sum(mult * c_(b08201_col(k, j)) for j, mult in ((2, 1), (3, 2), (4, 3))) + c_(b08201_col(k, 5)) * m[None, :]
        V[:, :, ci] = v / N[:, :, ci]
    return dict(N=N, N0=N0, V=V, Z=N0 / N, m=m, total_class=(V * N).sum(axis=2) / N.sum(axis=2),
                total_b25046=c_("B25046_001E") / c_("B25044_001E"))


def add_size_columns(d: pd.DataFrame, m: np.ndarray) -> None:
    """District-level household-size columns (in place): households n1..n4 by size class, vehicles per
    household v1..v4 and zero-vehicle shares z1..z4 within each class, class shares sh1..sh4."""
    mloc = np.where(d.loc_id == 1, m[0], m[1])
    hh = d["B25044_001E"]
    for k in SIZE_KS:
        nk = d[b08201_col(k, 0)]
        vk = sum(mult * d[b08201_col(k, j)] for j, mult in ((2, 1), (3, 2), (4, 3))) + mloc * d[b08201_col(k, 5)]
        d[f"n{k}"] = nk
        d[f"v{k}"] = vk / nk
        d[f"z{k}"] = d[b08201_col(k, 1)] / nk
        d[f"sh{k}"] = nk / hh


def size_class_analysis(pools: dict[str, Any], d: pd.DataFrame, g: pd.DataFrame, boot: Boot,
                        cz_index: dict[Any, int]) -> dict[str, Any]:
    """Vehicles per household by household-size class, the size-standardized suburb-city gap and the
    within-class district regressions with CZ fixed effects and log median household income.

    Standardization: in each CZ the class gaps (suburb minus city, within the class) are averaged with the
    weights of the CZ's pooled (city + suburb) household-size distribution (also with the city's and the
    suburb's own weights); the district regressions give the income-adjusted class gaps (coefficient on the
    suburb indicator, weights = households of the class) and are averaged with the national shares of the
    classes among the partitioned districts' households.  Bootstrap: whole CZs, shared draws; the
    standardized regression gap uses the draws of the class regressions."""
    N, V, Z, N0 = pools["N"], pools["V"], pools["Z"], pools["N0"]
    out: dict[str, Any] = dict(classes=[SIZE_LABELS[k] for k in SIZE_KS], four_plus_vehicles_city=float(pools["m"][0]),
                               four_plus_vehicles_suburb=float(pools["m"][1]))
    diff = pools["total_class"] - pools["total_b25046"]
    out["check_vs_b25046"] = dict(
        max_abs_diff=float(np.nanmax(np.abs(diff))), mean_diff=float(np.nanmean(diff)),
        raw_gap_class_based=boot.wmean(pools["total_class"][:, 1] - pools["total_class"][:, 0])["value"],
        raw_gap_b25046=boot.wmean(pools["total_b25046"][:, 1] - pools["total_b25046"][:, 0])["value"])
    share = N / N.sum(axis=2, keepdims=True)
    out["composition"] = {SIZE_LABELS[k]: boot.wmean(share[:, 1, ci] - share[:, 0, ci], share[:, 0, ci], share[:, 1, ci])
                          for ci, k in enumerate(SIZE_KS)}
    out["by_class"] = {SIZE_LABELS[k]: dict(
        vehicles=boot.wmean(V[:, 1, ci] - V[:, 0, ci], V[:, 0, ci], V[:, 1, ci]),
        zero_share=boot.wmean(Z[:, 1, ci] - Z[:, 0, ci], Z[:, 0, ci], Z[:, 1, ci])) for ci, k in enumerate(SIZE_KS)}
    pooled = N[:, 0, :] + N[:, 1, :]
    for tag, wts in (("pooled", pooled), ("city_weights", N[:, 0, :]), ("suburb_weights", N[:, 1, :])):
        pi = wts / wts.sum(axis=1, keepdims=True)
        out[f"standardized_gap_{tag}"] = boot.wmean((pi * (V[:, 1, :] - V[:, 0, :])).sum(axis=1))
        out[f"standardized_zero_share_gap_{tag}"] = boot.wmean((pi * (Z[:, 1, :] - Z[:, 0, :])).sum(axis=1))
    # within-class district regressions (CZ fixed effects, households of the class as weights)
    nat = np.array([float(d[f"n{k}"].sum()) for k in SIZE_KS])
    nat = nat / nat.sum()
    out["national_class_shares"] = nat.tolist()
    regs: dict[str, Any] = {}
    draws: dict[str, np.ndarray] = {}

    def run(y: str, xs: list[str], w: str) -> dict[str, Any]:
        dd = d[["czone"] + xs + [y, w]].replace([np.inf, -np.inf], np.nan).dropna()
        dd = dd[dd[w] > 0]
        return wls_fe_boot(dd[y].to_numpy(float), dd[xs].to_numpy(float), dd[w].to_numpy(float), dd.czone.to_numpy(), boot, cz_index)

    for ci, k in enumerate(SIZE_KS):
        for outcome, yv in (("vehicles", f"v{k}"), ("zero_share", f"z{k}")):
            for spec, xs in (("no_controls", ["suburb"]), ("log_income", ["suburb", "ln_inc"])):
                r = run(yv, xs, f"n{k}")
                key = f"{outcome}|{SIZE_LABELS[k]}|{spec}"
                regs[key] = dict(outcome=outcome, size_class=SIZE_LABELS[k], spec=spec, regressors=xs,
                                 **{kk: vv for kk, vv in r.items() if kk != "draws"})
                draws[key] = r["draws"][:, 0]
    out["regressions"] = regs

    def combo(outcome: str, spec: str) -> dict[str, Any]:
        est = float(sum(nat[ci] * regs[f"{outcome}|{SIZE_LABELS[k]}|{spec}"]["beta"][0] for ci, k in enumerate(SIZE_KS)))
        dr = sum(nat[ci] * draws[f"{outcome}|{SIZE_LABELS[k]}|{spec}"] for ci, k in enumerate(SIZE_KS))
        se, lo, hi = boot_se_of(dr)
        return dict(value=est, boot_se=se, ci95_lo=lo, ci95_hi=hi, draws=dr)

    for outcome in ("vehicles", "zero_share"):
        out[f"standardized_regression_{outcome}"] = dict(
            income_adjusted=combo(outcome, "log_income"), raw=combo(outcome, "no_controls"))
    return out


def part_d(panel: pd.DataFrame, part: pd.DataFrame, g: pd.DataFrame, boot: Boot,
           ce: dict[str, Any], indir: Path, size_table: dict[str, Any] | None = None) -> dict[str, Any]:
    """(d) vehicles by location (pooled over districts, whole-CZ bootstrap) and the within-CZ
    district regressions; conversion into a fixed goods cost in units of mean consumption."""
    czones = g.czone.to_numpy()
    cz_index = {c: i for i, c in enumerate(czones)}

    def pooled(cols: list[str]) -> np.ndarray:
        return pooled_by_location(part, cols, czones)[:, :, :]

    # pooled sums by CZ x location
    cols = ["B25044_001E", "B25044_002E", "B25044_003E", "B25044_009E", "B25044_010E", "B25046_001E",
            "B25046_002E", "B25046_003E", "B08201_001E", "B08201_002E", "B08201_007E", "B08201_008E",
            "B25044_004E", "B25044_005E", "B25044_006E", "B25044_007E", "B25044_008E",
            "B25044_011E", "B25044_012E", "B25044_013E", "B25044_014E", "B25044_015E"]
    P = pooled(cols)
    ix = {c: i for i, c in enumerate(cols)}
    c_ = lambda name: P[:, :, ix[name]]
    stats: dict[str, Any] = {}

    def stat(key: str, label: str, num: np.ndarray, den: np.ndarray, kind: str) -> None:
        r = num / den                                      # (n_cz, 2)
        if kind == "log":
            x = np.log(r[:, 1] / r[:, 0])
        else:
            x = r[:, 1] - r[:, 0]
        st = boot.wmean(x, r[:, 0], r[:, 1])
        stats[key] = dict(label=label, kind=kind, **st)

    z_all = c_("B25044_003E") + c_("B25044_010E")
    stat("zero_vehicle_share", "households with no vehicle available, share of occupied units (B25044)",
         z_all, c_("B25044_001E"), "diff")
    stat("vehicles_per_hh", "vehicles per occupied unit: aggregate vehicles (B25046) / occupied units (B25044)",
         c_("B25046_001E"), c_("B25044_001E"), "diff")
    stat("log_vehicles_per_hh", "log ratio, vehicles per occupied unit (suburb / city)",
         c_("B25046_001E"), c_("B25044_001E"), "log")
    stat("zero_vehicle_share_owners", "owner-occupied units with no vehicle (B25044)",
         c_("B25044_003E"), c_("B25044_002E"), "diff")
    stat("zero_vehicle_share_renters", "renter-occupied units with no vehicle (B25044)",
         c_("B25044_010E"), c_("B25044_009E"), "diff")
    stat("vehicles_per_owner_hh", "vehicles per owner-occupied unit (B25046_002 / B25044_002)",
         c_("B25046_002E"), c_("B25044_002E"), "diff")
    stat("vehicles_per_renter_hh", "vehicles per renter-occupied unit (B25046_003 / B25044_009)",
         c_("B25046_003E"), c_("B25044_009E"), "diff")
    stat("renter_share", "renter share of occupied units (B25044_009 / B25044_001)",
         c_("B25044_009E"), c_("B25044_001E"), "diff")
    stat("two_plus_vehicles", "units with two or more vehicles (B25044), share",
         c_("B25044_005E") + c_("B25044_006E") + c_("B25044_007E") + c_("B25044_008E")
         + c_("B25044_012E") + c_("B25044_013E") + c_("B25044_014E") + c_("B25044_015E"),
         c_("B25044_001E"), "diff")
    stat("zero_vehicle_share_b08201", "households with no vehicle (B08201_002 / B08201_001)",
         c_("B08201_002E"), c_("B08201_001E"), "diff")
    stat("zero_vehicle_share_one_person", "one-person households with no vehicle (B08201_008 / B08201_007)",
         c_("B08201_008E"), c_("B08201_007E"), "diff")

    # ---- household-size classes (B08201): pools by CZ x location, and district-level columns
    pools = size_pools(part, czones)

    # ---- within-CZ district regressions
    d = part.copy()
    d["suburb"] = (d.loc_id == 2).astype(float)
    d["hh"] = d["B25044_001E"]
    d["zero_share"] = (d["B25044_003E"] + d["B25044_010E"]) / d["B25044_001E"]
    d["veh_per_hh"] = d["B25046_001E"] / d["B25044_001E"]
    d["renter_share"] = d["B25044_009E"] / d["B25044_001E"]
    d["ln_inc"] = np.log(d["median_hh_income"].where(d["median_hh_income"] > 0))
    d["ln_inc_sq"] = (d["ln_inc"] - d["ln_inc"].mean()) ** 2
    add_size_columns(d, pools["m"])
    reg: dict[str, Any] = {}
    for yname in ("zero_share", "veh_per_hh"):
        for wname, wcol in (("households", "hh"), ("pupils", "pupils")):
            for spec, xs in (("no_controls", ["suburb"]), ("log_income", ["suburb", "ln_inc"]),
                             ("log_income_quadratic", ["suburb", "ln_inc", "ln_inc_sq"]),
                             ("log_income_renter_share", ["suburb", "ln_inc", "renter_share"]),
                             ("log_income_size_shares", ["suburb", "ln_inc", "sh1", "sh3", "sh4"])):
                # the sample of the first four specifications is the one of the first version of this report
                # (the size-share columns are required only by the last)
                need = ["czone", "suburb", "ln_inc", "ln_inc_sq", "renter_share"] + (["sh1", "sh3", "sh4"] if "sh1" in xs else [])
                dd = d[need + [yname, wcol]].replace([np.inf, -np.inf], np.nan).dropna()
                dd = dd[(dd[wcol] > 0)]
                r = wls_fe_boot(dd[yname].to_numpy(float), dd[xs].to_numpy(float), dd[wcol].to_numpy(float),
                                dd.czone.to_numpy(), boot, cz_index)
                reg[f"{yname}|{wname}|{spec}"] = dict(
                    y=yname, weights=wname, spec=spec, regressors=xs, **{k: v for k, v in r.items() if k != "draws"})
    # unweighted district means of the outcomes for reference (levels by location)
    lev = {}
    for yname in ("zero_share", "veh_per_hh"):
        for loc, lab in ((1, "city"), (2, "suburb")):
            m = d.loc[d.loc_id == loc]
            lev[f"{yname}_{lab}_household_weighted"] = float(np.average(
                m[yname].dropna(), weights=m.loc[m[yname].notna(), "hh"]))
    out: dict[str, Any] = dict(stats=stats, regressions=reg, district_means=lev,
                               n_districts=int(len(d)), n_cz=int(len(g)))
    out["household_size"] = size_class_analysis(pools, d, g, boot, cz_index)
    out["ce_size_table"] = size_table

    # ---- conversion into a fixed goods cost in units of mean consumption
    cet = ce["quintiles"]
    veh = cet["chars"].get("Vehicles")
    tr = cet["transport"]
    allt = cet["transport_all"]
    veh_all = veh[0]
    private_all = allt["Vehicle purchases (net outlay)"] + allt["Gasoline, other fuels, and motor oil"] \
        + allt["Other vehicle expenses"]
    fixed_like = allt["Vehicle purchases (net outlay)"] + allt["Vehicle finance charges"] \
        + allt["Vehicle insurance"] + allt["Vehicle rental, leases, licenses, and other charges"]
    per_veh = {
        "ce_private_vehicle_spending_per_vehicle": private_all / veh_all,
        "ce_ownership_like_per_vehicle": fixed_like / veh_all,
        "aaa_ownership_per_vehicle_year": AAA_2018["ownership_per_year"],
        "aaa_total_per_vehicle_year_15k_miles": AAA_2018["total_per_year_15k"],
        "aaa_total_per_vehicle_year_10k_miles": AAA_2018["total_per_year_10k"]}
    out["cost_per_vehicle_year"] = per_veh
    # the CE publishes vehicles per consumer unit to one decimal: +-0.05 on 1.9 vehicles is +-2.6% on every CE
    # cost per vehicle (and the vehicle-count elasticity below is computed from rounded counts)
    out["ce_vehicles_rounding"] = dict(published_vehicles_per_cu=float(veh_all), half_unit=0.05,
                                       relative_bound=0.05 / float(veh_all),
                                       per_vehicle_private_lo=private_all / (veh_all + 0.05),
                                       per_vehicle_private_hi=private_all / (veh_all - 0.05))
    out["ce_per_quintile_vehicle_spending"] = dict(
        vehicles=veh[1:],
        private_vehicle_spending=[tr["Vehicle purchases (net outlay)"][q] + tr["Gasoline, other fuels, and motor oil"][q]
                                  + tr["Other vehicle expenses"][q] for q in range(5)],
        per_vehicle=[(tr["Vehicle purchases (net outlay)"][q] + tr["Gasoline, other fuels, and motor oil"][q]
                      + tr["Other vehicle expenses"][q]) / veh[1 + q] for q in range(5)],
        transportation=tr["Transportation"], share_transportation=[tr["Transportation"][q] / cet["mean_expenditures"][q]
                                                                    for q in range(5)],
        at_least_one_vehicle_pct=cet["chars"].get("At least one vehicle owned or leased"),
        mean_expenditures=cet["mean_expenditures"])
    # income (expenditure) elasticity of the vehicle count and of vehicle spending across the CE groups
    for tag in ("quintiles", "deciles"):
        b = ce[tag]
        v = b["chars"].get("Vehicles")
        Yg = np.array(b["mean_expenditures"])
        if v is not None:
            vg = np.array(v[1:], float)
            sl = np.polyfit(np.log(Yg), np.log(vg), 1)[0]
            tt = b["transport"]
            pv = np.array(tt["Vehicle purchases (net outlay)"]) + np.array(tt["Gasoline, other fuels, and motor oil"]) \
                + np.array(tt["Other vehicle expenses"])
            sl2 = np.polyfit(np.log(Yg), np.log(pv), 1)[0]
            out[f"ce_expenditure_elasticity_{tag}"] = dict(vehicles=float(sl), private_vehicle_spending=float(sl2))
    return out


def car_requirement(d: dict[str, Any]) -> dict[str, Any]:
    """Convert the suburb-city difference in vehicles into a fixed goods cost r_car in units of mean
    consumption (CE 2018 mean annual expenditures per consumer unit, $61,224.13):
    r_car = (extra vehicles per household) x (annual cost of a vehicle) / mean expenditures.

    Extra vehicles (per household), from the unadjusted to the composition-adjusted:
      raw pooled gap (B25046 / B25044);
      district regression with CZ fixed effects and log median household income (the first version's
        central value, NOT adjusted for household size: suburban households are larger);
      size-standardized gap: the within-size-class gaps (1, 2, 3, 4+ persons, B08201) averaged with the
        CZ's pooled household-size distribution;
      size-standardized and income-adjusted: within-class district regressions with log income, averaged
        with the national class shares (CENTRAL);
      one-person households (the model's single adult), raw and income-adjusted;
      the district regression with log income and the district shares of 1-, 3- and 4+-person households
        as controls (LOW: it over-controls, since the size mix is itself a location outcome).
    Cost per vehicle: CE private-vehicle spending per vehicle (purchases net outlay + fuel + other vehicle
    expenses over vehicles per consumer unit; published to one decimal), its ownership-like part, AAA
    ownership costs, AAA total cost at 10,000 and 15,000 miles (new vehicle).  A one-person household's
    ratio uses the one-person consumer unit's mean expenditures (CE Table 1400) as denominator: the ratio of
    a per-household cost to per-household consumption is unit-free only if both scale with household size."""
    reg = d["regressions"]
    hs = d["household_size"]
    one = SIZE_LABELS[1]
    cost = d["cost_per_vehicle_year"]
    st = d["stats"]
    c_ = lambda r, k=0: dict(value=r["beta"][k], se=r["se_boot"][k])
    ext = {
        "vehicles_raw_pooled": dict(value=st["vehicles_per_hh"]["value"], se=st["vehicles_per_hh"]["boot_se"]),
        "vehicles_income_adjusted": c_(reg["veh_per_hh|households|log_income"]),
        "vehicles_size_standardized": dict(value=hs["standardized_gap_pooled"]["value"], se=hs["standardized_gap_pooled"]["boot_se"]),
        "vehicles_size_standardized_income_adjusted": dict(
            value=hs["standardized_regression_vehicles"]["income_adjusted"]["value"],
            se=hs["standardized_regression_vehicles"]["income_adjusted"]["boot_se"]),
        "vehicles_one_person_raw": dict(value=hs["by_class"][one]["vehicles"]["value"], se=hs["by_class"][one]["vehicles"]["boot_se"]),
        "vehicles_one_person_income_adjusted": c_(hs["regressions"][f"vehicles|{one}|log_income"]),
        "vehicles_size_controls_regression": c_(reg["veh_per_hh|households|log_income_size_shares"]),
        "vehicles_raw_regression": c_(reg["veh_per_hh|households|no_controls"]),
        "no_vehicle_share_raw_pooled": dict(value=-st["zero_vehicle_share"]["value"], se=st["zero_vehicle_share"]["boot_se"]),
        "no_vehicle_share_income_adjusted": dict(value=-reg["zero_share|households|log_income"]["beta"][0],
                                                 se=reg["zero_share|households|log_income"]["se_boot"][0]),
        "no_vehicle_share_size_standardized_income_adjusted": dict(
            value=-hs["standardized_regression_zero_share"]["income_adjusted"]["value"],
            se=hs["standardized_regression_zero_share"]["income_adjusted"]["boot_se"]),
    }
    labels = {
        "vehicles_raw_pooled": "extra vehicles, raw (not adjusted)",
        "vehicles_income_adjusted": "extra vehicles, income-adjusted (not size-adjusted: first version)",
        "vehicles_size_standardized": "extra vehicles, size-standardized",
        "vehicles_size_standardized_income_adjusted": "extra vehicles, size-standardized and income-adjusted",
        "vehicles_one_person_raw": "extra vehicles, one-person households, raw",
        "vehicles_one_person_income_adjusted": "extra vehicles, one-person households, income-adjusted",
        "vehicles_size_controls_regression": "extra vehicles, income and household-size shares as controls",
        "vehicles_raw_regression": "extra vehicles, regression without controls",
        "no_vehicle_share_raw_pooled": "fewer zero-vehicle households, raw",
        "no_vehicle_share_income_adjusted": "fewer zero-vehicle households, income-adjusted (not size-adjusted)",
        "no_vehicle_share_size_standardized_income_adjusted": "fewer zero-vehicle households, size-standardized and income-adjusted"}
    size = d.get("ce_size_table")
    one_exp = size["groups"]["One person"]["expenditures"] if size else float("nan")
    rows = []
    for dk in ("vehicles_raw_pooled", "vehicles_income_adjusted", "vehicles_size_standardized", "vehicles_size_standardized_income_adjusted",
               "vehicles_size_controls_regression", "no_vehicle_share_income_adjusted",
               "no_vehicle_share_size_standardized_income_adjusted"):
        for ck, cv in cost.items():
            e = ext[dk]
            rows.append(dict(extra=dk, extra_label=labels[dk], cost=ck, extra_value=e["value"], extra_se=e["se"], cost_usd=cv,
                             denominator=CE_TOTAL_EXP, denominator_label="mean expenditures per CU",
                             r_car=e["value"] * cv / CE_TOTAL_EXP, r_car_se=e["se"] * cv / CE_TOTAL_EXP, annual_usd=e["value"] * cv))
    for dk in ("vehicles_one_person_raw", "vehicles_one_person_income_adjusted"):
        for ck in ("ce_private_vehicle_spending_per_vehicle", "ce_ownership_like_per_vehicle", "aaa_total_per_vehicle_year_15k_miles"):
            e, cv = ext[dk], cost[ck]
            rows.append(dict(extra=dk, extra_label=labels[dk], cost=ck, extra_value=e["value"], extra_se=e["se"], cost_usd=cv,
                             denominator=one_exp, denominator_label="mean expenditures of one-person CUs",
                             r_car=e["value"] * cv / one_exp, r_car_se=e["se"] * cv / one_exp, annual_usd=e["value"] * cv))
    cp = cost["ce_private_vehicle_spending_per_vehicle"]
    aaa = cost["aaa_total_per_vehicle_year_15k_miles"]
    cen = ext["vehicles_size_standardized_income_adjusted"]
    prev = ext["vehicles_income_adjusted"]
    return dict(
        extra_vehicles=ext, extra_labels=labels, table=rows,
        central=dict(r_car=cen["value"] * cp / CE_TOTAL_EXP, se=cen["se"] * cp / CE_TOTAL_EXP,
                     rule="size-standardized, income-adjusted extra vehicles x CE private-vehicle spending per vehicle / mean "
                          "expenditures per consumer unit"),
        low=dict(r_car=ext["vehicles_size_controls_regression"]["value"] * cp / CE_TOTAL_EXP,
                 rule="suburb coefficient with log income and household-size shares as controls x CE private-vehicle spending per vehicle"),
        high=dict(r_car=ext["vehicles_size_standardized"]["value"] * aaa / CE_TOTAL_EXP,
                  rule="size-standardized extra vehicles (no income adjustment) x AAA total cost per vehicle (15,000 miles)"),
        single_adult=dict(r_car=ext["vehicles_one_person_income_adjusted"]["value"] * cp / one_exp,
                          se=ext["vehicles_one_person_income_adjusted"]["se"] * cp / one_exp,
                          denominator=one_exp, rule="income-adjusted extra vehicles of one-person households x CE private-vehicle "
                                                    "spending per vehicle / mean expenditures of one-person consumer units (CE Table 1400)"),
        first_version=dict(r_car=prev["value"] * cp / CE_TOTAL_EXP,
                           rule="income-adjusted extra vehicles, no household-size adjustment (the first version's central value)"),
        unadjusted_upper=dict(r_car=ext["vehicles_raw_pooled"]["value"] * aaa / CE_TOTAL_EXP,
                              rule="raw pooled extra vehicles x AAA total cost per vehicle (15,000 miles): no adjustment (the first version's upper end)"),
        extensive_margin_central=ext["no_vehicle_share_size_standardized_income_adjusted"]["value"] * cp / CE_TOTAL_EXP,
        ce_total_exp=CE_TOTAL_EXP)


# --------------------------------------------------------------------------
# (e) BFM (2007) WTP-income gradient and the 1990 Bay Area college / income association
# --------------------------------------------------------------------------

BFM = dict(
    source="Bayer, Ferreira and McMillan (2007), JPE 115(4): 588-638, Table 8 (p. 625), preferred "
           "specification (boundary fixed effects, neighborhood sociodemographics), 'Average Test Score "
           "+1 SD' column; household income and college row are interaction terms of the "
           "heterogeneous-choice model",
    mean_mwtp=19.69, mean_mwtp_se=7.41, per_10k_income=1.38, per_10k_income_se=0.33,
    college=13.03, college_se=3.57, children=7.41, children_se=3.58,
    mean_income=54103.0, sd_income=50719.0, share_college_heads=0.44,
    income_source="BFM (2004), NBER w10871, Table 1 (same 242,100 households), via literature_targets.md "
                  "section 1d")
BAY_PMSAS = {"5775": "Oakland PMSA (Alameda, Contra Costa)", "7360": "San Francisco PMSA (San Francisco, "
             "Marin, San Mateo)", "7400": "San Jose PMSA (Santa Clara)"}


def read_pums90_bay(path: Path) -> pd.DataFrame:
    """Households of the six Bay Area counties in the Census 1990 PUMS 5% sample (California
    file PUMSAXCA): housing weight, household income in 1989 (RHHINC), tenure, and the
    householder's educational attainment (YEARSCH of the RELAT1 = 00 person record).  Fixed-width
    layout from the PUMS data dictionary (PUMSUSDD): H record MSAPMSA cols 20-23, HOUSWGT
    29-32, PERSONS 33-34, GQINST 35, TENURE 44, RHHINC 141-147, R18UNDR 162; P record RELAT1
    9-10, AGE 15-16, YEARSCH 51-52.  Group quarters and vacant units are dropped."""
    rows = []
    cur: dict[str, Any] | None = None
    with zipfile.ZipFile(path) as z:
        name = z.infolist()[0].filename
        with z.open(name) as f:
            for raw in f:
                line = raw.decode("latin-1")
                t = line[0]
                if t == "H":
                    if line[19:23] in BAY_PMSAS and line[34] == "0" and int(line[32:34]) >= 1:
                        cur = dict(pmsa=line[19:23], w=int(line[28:32]), persons=int(line[32:34]),
                                   tenure=int(line[43]), hhinc=int(line[140:147]), under18=int(line[161]),
                                   yearsch=-1, age=-1)
                        rows.append(cur)
                    else:
                        cur = None
                elif cur is not None and cur["yearsch"] < 0 and line[8:10] == "00":
                    cur["yearsch"] = int(line[50:52])
                    cur["age"] = int(line[14:16])
    d = pd.DataFrame(rows)
    d = d[d.yearsch >= 0]
    d["college"] = (d.yearsch >= 14).astype(float)               # bachelor's degree or more
    d["some_college_plus"] = (d.yearsch >= 11).astype(float)
    d["assoc_plus"] = (d.yearsch >= 12).astype(float)            # associate degree (occupational or academic) or more
    d["child"] = (d.under18 > 0).astype(float)                   # R18UNDR: a person under 18 in the household
    return d


def part_e(path: Path, rng: np.random.Generator) -> dict[str, Any]:
    """(e) elasticity of BFM's marginal WTP for school quality with respect to household income at
    the means, with and without the college-degree term folded in through the 1990 Bay Area
    association of the householder's college degree with household income (PUMS 5%)."""
    b = BFM
    scale = b["mean_income"] / b["mean_mwtp"]
    by, se_by = b["per_10k_income"] / 1e4, b["per_10k_income_se"] / 1e4
    eps = by * scale
    out: dict[str, Any] = dict(bfm=b, elasticity_income=dict(
        value=eps, se_conditional_on_means=se_by * scale,
        se_with_mean_mwtp_independent=math.sqrt((se_by * scale) ** 2 + (eps * b["mean_mwtp_se"] / b["mean_mwtp"]) ** 2),
        formula="(1.38 / 10,000) x 54,103 / 19.69; SE = 0.33 / 10,000 x 54,103 / 19.69 holding the means fixed",
        implied_share_elasticity=eps - 1.0))
    d = read_pums90_bay(path)
    w, y, c = d.w.to_numpy(float), d.hhinc.to_numpy(float), d.college.to_numpy(float)

    def moments(w, y, c):
        ybar, cbar = np.average(y, weights=w), np.average(c, weights=w)
        vy = np.average((y - ybar) ** 2, weights=w)
        vc = np.average((c - cbar) ** 2, weights=w)
        cov = np.average((y - ybar) * (c - cbar), weights=w)
        return ybar, cbar, math.sqrt(vy), cov / vy, cov / math.sqrt(vy * vc), cbar, vc

    ybar, cbar, sdy, gamma, corr, _, vc = moments(w, y, c)
    ca, ch = d.assoc_plus.to_numpy(float), d.child.to_numpy(float)
    gam_boot, gam_boot_a, gam_boot_ch = [], [], []
    n = len(d)
    for _ in range(300):
        i = rng.integers(0, n, n)                                   # the same resamples for every slope
        gam_boot.append(moments(w[i], y[i], c[i])[3])
        gam_boot_a.append(moments(w[i], y[i], ca[i])[3])
        gam_boot_ch.append(moments(w[i], y[i], ch[i])[3])
    gam_se = float(np.std(gam_boot, ddof=1))
    gam_se_a, gam_se_ch = float(np.std(gam_boot_a, ddof=1)), float(np.std(gam_boot_ch, ddof=1))
    _, cbar_a, _, gamma_a, corr_a, _, _ = moments(w, y, ca)
    _, cbar_ch, _, gamma_ch, corr_ch, _, _ = moments(w, y, ch)
    # local slope of E[college | income] at the BFM mean income: weighted logit of college on log income
    pos = y > 0
    x = np.log(y[pos])
    beta = np.array([0.0, 0.0])
    Xl = np.column_stack([np.ones(pos.sum()), x])
    for _ in range(50):
        p = 1 / (1 + np.exp(-(Xl @ beta)))
        W = w[pos] * p * (1 - p)
        beta = beta + np.linalg.solve((Xl * W[:, None]).T @ Xl, Xl.T @ (w[pos] * (c[pos] - p)))
    p0 = 1 / (1 + np.exp(-(beta[0] + beta[1] * math.log(b["mean_income"]))))
    slope_logit = float(p0 * (1 - p0) * beta[1] / b["mean_income"])
    # local linear slope in a +- $20,000 window around the BFM mean
    m = (np.abs(y - b["mean_income"]) <= 20000)
    Xw = np.column_stack([np.ones(m.sum()), y[m] - b["mean_income"]])
    coef = np.linalg.solve((Xw * w[m][:, None]).T @ Xw, Xw.T @ (w[m] * c[m]))
    slope_local = float(coef[1])
    # college / non-college mean incomes
    inc1, inc0 = np.average(y[c == 1], weights=w[c == 1]), np.average(y[c == 0], weights=w[c == 0])
    out["pums90"] = dict(
        source="Census 1990 PUMS 5% sample, California file (https://www2.census.gov/census_1990/1990_PUMS_A/"
               "PUMSAXCA.zip), households in the Oakland (5775), San Francisco (7360) and San Jose (7400) "
               "PMSAs = the six counties of BFM (Alameda, Contra Costa, Marin, San Francisco, San Mateo, "
               "Santa Clara; 1990 OMB definitions); non-group-quarters occupied units; weight HOUSWGT; income = RHHINC "
               "(total household income in 1989 dollars); college = householder with a bachelor's degree "
               "or more (YEARSCH >= 14)",
        n_households=int(n), weighted_households=float(w.sum()), mean_income=float(ybar),
        sd_income=float(sdy), share_college_head=float(cbar), corr_college_income=float(corr),
        gamma_per_dollar=float(gamma), gamma_per_10k=float(gamma * 1e4), gamma_se_boot=gam_se,
        gamma_per_cv_income=float(gamma * ybar),
        slope_logit_at_bfm_mean_per_10k=slope_logit * 1e4, slope_local_linear_per_10k=slope_local * 1e4,
        mean_income_college=float(inc1), mean_income_noncollege=float(inc0),
        p_college_at_mean=float(p0), bfm_mean_income=b["mean_income"], bfm_sd_income=b["sd_income"],
        bfm_share_college_heads=b["share_college_heads"], by_pmsa={k: dict(n=int((d.pmsa == k).sum()),
                                                     share_college=float(np.average(d.college[d.pmsa == k], weights=d.w[d.pmsa == k])),
                                                     mean_income=float(np.average(d.hhinc[d.pmsa == k], weights=d.w[d.pmsa == k])))
                                               for k in BAY_PMSAS})
    # folded-in elasticities: income slope + college coefficient x slope of college on income
    def folded(g_per_dollar: float, g_se: float) -> dict[str, float]:
        tot = by + b["college"] * g_per_dollar
        se = math.sqrt(se_by ** 2 + (b["college"] * g_se) ** 2 + (g_per_dollar * b["college_se"]) ** 2)
        return dict(slope_per_10k=tot * 1e4, elasticity=tot * scale, se=se * scale,
                    college_term_per_10k=b["college"] * g_per_dollar * 1e4,
                    college_term_elasticity=b["college"] * g_per_dollar * scale)
    g_bfm = corr * math.sqrt(b["share_college_heads"] * (1 - b["share_college_heads"])) / b["sd_income"]
    out["pums90"].update(
        share_assoc_plus_head=float(cbar_a), corr_assoc_plus_income=float(corr_a), gamma_assoc_plus_per_10k=float(gamma_a * 1e4),
        gamma_assoc_plus_se_boot=gam_se_a, share_child_present=float(cbar_ch), corr_child_income=float(corr_ch),
        gamma_child_per_10k=float(gamma_ch * 1e4), gamma_child_se_boot=gam_se_ch)

    def folded_terms(terms: list[tuple[float, float, float, float]]) -> dict[str, float]:
        """Income slope plus sum of coef x slope-of-characteristic-on-income terms; terms = (Table 8 coefficient, its SE,
        slope per dollar, its bootstrap SE); independent sampling errors."""
        tot = by + sum(c0 * g0 for c0, _, g0, _ in terms)
        var = se_by ** 2 + sum((c0 * gs) ** 2 + (g0 * cs) ** 2 for c0, cs, g0, gs in terms)
        return dict(slope_per_10k=tot * 1e4, elasticity=tot * scale, se=math.sqrt(var) * scale,
                    terms_per_10k=[c0 * g0 * 1e4 for c0, _, g0, _ in terms], terms_elasticity=[c0 * g0 * scale for c0, _, g0, _ in terms])

    col_b = (b["college"], b["college_se"], gamma, gam_se)
    col_a = (b["college"], b["college_se"], gamma_a, gam_se_a)
    kid = (b["children"], b["children_se"], gamma_ch, gam_se_ch)
    out["elasticity_variants"] = dict(
        bachelors=folded_terms([col_b]), associate_plus=folded_terms([col_a]),
        bachelors_children=folded_terms([col_b, kid]), associate_plus_children=folded_terms([col_a, kid]),
        definitions="college = householder with a bachelor's degree or more (YEARSCH >= 14) or an associate degree or more "
                    "(YEARSCH >= 12); BFM's text does not define 'college degree or more' (their Table 1 has 44% of the heads "
                    "of household of a census block group with a college degree or more, block-group average); children = a "
                    "person under 18 in the household (R18UNDR); slopes are the linear projections Cov(x, income) / Var(income)")
    out["elasticity_with_college"] = dict(
        linear_projection=folded(gamma, gam_se),
        local_linear=folded(slope_local, gam_se),
        local_logit=folded(slope_logit, gam_se),
        bfm_moments_pums_correlation=folded(g_bfm, gam_se),
        assumptions="independent sampling errors of the Table 8 coefficients (their covariance is not "
                    "published; income and college are positively correlated so it is likely negative, "
                    "which would lower the SE of the sum); the other Table 8 interactions (children, race, "
                    "work status, age) are held fixed; the college dummy changes by gamma per dollar of "
                    "income, gamma = Cov(college, income) / Var(income) in the PUMS Bay Area households "
                    "(linear projection), or the local slope of E[college | income] at $54,103 (logit, local "
                    "linear window +-$20,000)")
    # WTP as a share of monthly income at income multiples of the mean, linear Table 8 gradient
    rows = []
    for mult in (0.5, 1.0, 2.0):
        inc = mult * b["mean_income"]
        wtp = b["mean_mwtp"] + by * (inc - b["mean_income"])
        rows.append(dict(income_multiple=mult, income=inc, wtp=wtp, wtp_share_monthly_income=wtp / (inc / 12)))
    out["wtp_share_by_income"] = rows
    return out


# --------------------------------------------------------------------------
# Recommended model inputs
# --------------------------------------------------------------------------

SPECS = ("FMM_quad", "FMM_lin", "IV_quad", "IV_lin")


def implied_delta(ratio: float) -> float:
    """Annual depreciation rate at which the 10-year / 30-year ratio of observed elasticities,
    (1 - (1 - delta)^10) / (1 - (1 - delta)^30), equals `ratio` (bisection; the ratio rises from 1/3 at
    delta -> 0 to 1 at delta -> 1)."""
    lo, hi = 1e-6, 0.999
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        r = lr_factor(SAIZ_YEARS, mid) / lr_factor(BSH_YEARS, mid)
        lo, hi = (mid, hi) if r < ratio else (lo, mid)
    return 0.5 * (lo + hi)


def eps_recommendation(c: dict[str, Any], city_share: float) -> dict[str, Any]:
    """Long-run housing supply elasticities for the city and the suburb, ONE horizon convention.

    Every anchor is a long-run elasticity of the housing stock (units), obtained from the published
    T-year estimate with eps_LR = eps_T / (1 - (1 - delta)^T), delta = 3.5% a year (Baum-Snow and
    Duranton, slide 76): Saiz (2010) is a 30-year estimate (x1.52), Baum-Snow and Han (2024) and the
    Saiz-method replication of their Online Appendix D are 10-year estimates (x3.34).
      central metro anchor = Saiz's population-weighted mean over the 95 MSAs of Table VI, converted;
      low  = Baum-Snow and Han's tract elasticities (housing-unit-weighted mean, FMM linear), converted;
      high = Saiz's unweighted mean over the 269 MSAs, converted.
    Period check: the horizon alone predicts that Saiz's method on 2000-10 data gives 0.46 of its
    1970-2000 value; Online Appendix D finds 1.25 / 2.54 = 0.49, so no extra period adjustment is made.
    Within-metro ratio suburb / city: Baum-Snow and Han's predicted tract elasticities on the city /
    suburb split matched to the pupil share of the locale partition (central: mean of the four
    estimators for units; range: every variant, i.e. also units / floorspace, the other city shares
    and the within-region medians).  With city share s the metro elasticity is s eps_c + (1 - s) eps_s
    (equal price changes in the two places, perfect substitutes, Baum-Snow and Duranton slide 77), so
    eps_c = anchor / (s + (1 - s) ratio) and eps_s = ratio x eps_c.  The values are constructed from
    the literature, not a direct estimate for the model's two locations."""
    sz, bsh, ad = c["saiz"], c["baum_snow_han"], c["appendix_d"]
    s = city_share
    f30, f10 = lr_factor(SAIZ_YEARS, DELTA_CENTRAL), lr_factor(BSH_YEARS, DELTA_CENTRAL)
    shares = bsh["city_suburb"]
    split = shares[f"city_share_{bsh['city_share_used']}"]
    ru = [split[f"{sp}_units"]["ratio"] for sp in SPECS]
    rs = [split[f"{sp}_space"]["ratio"] for sp in SPECS]
    rho = float(np.mean(ru))
    variants: list[dict[str, Any]] = []
    for rec in shares.values():
        sh = rec["city_share_of_housing_units"]
        for sp in SPECS:
            for t in ("units", "space"):
                variants.append(dict(city_share=sh, kind="pooled means", spec=sp, quantity=t, ratio=rec[f"{sp}_{t}"]["ratio"]))
        for sp in ("FMM_lin", "FMM_quad"):
            for t in ("units", "space"):
                w_ = rec.get(f"{sp}_{t}_within_region_ratio")
                if w_ is not None:
                    variants.append(dict(city_share=sh, kind="within-region median", spec=sp, quantity=t, ratio=w_["median"]))
                    variants.append(dict(city_share=sh, kind="within-region housing-weighted median", spec=sp, quantity=t,
                                         ratio=w_["weighted_median"]))
    vr = np.array([v["ratio"] for v in variants], float)
    vr_ok = np.isfinite(vr)
    rho_lo, rho_hi = float(vr[vr_ok].min()), float(vr[vr_ok].max())
    v_lo, v_hi = variants[int(np.nanargmin(vr))], variants[int(np.nanargmax(vr))]

    sz95, sz269 = sz["msas_over_500k"], sz["all_msas"]
    bsh_hw = bsh["overall"]["FMM_lin_units"]["mean_housing_units"]
    a8 = ad["table_a8_means_2000_2010"]
    ratio_obs = ad["saiz_method_2000_2010_mean"] / ad["saiz_1970_2000_mean"]
    ratio_pred = (1.0 - (1.0 - DELTA_CENTRAL) ** BSH_YEARS) / (1.0 - (1.0 - DELTA_CENTRAL) ** SAIZ_YEARS)
    anchors = {
        "saiz_pop95": dict(label="Saiz (2010), population-weighted mean over the 95 MSAs of Table VI", years=SAIZ_YEARS,
                           published=sz95["pop_weighted_mean"], long_run=sz95["pop_weighted_mean"] * f30, role="central"),
        "saiz_pop269": dict(label="Saiz (2010), population-weighted mean over all 269 MSAs", years=SAIZ_YEARS,
                            published=sz269["pop_weighted_mean"], long_run=sz269["pop_weighted_mean"] * f30, role=""),
        "saiz_mean269": dict(label="Saiz (2010), unweighted mean over all 269 MSAs", years=SAIZ_YEARS,
                             published=sz269["mean"], long_run=sz269["mean"] * f30, role="high"),
        "saiz_p25_95": dict(label="Saiz (2010), 25th percentile of the 95 MSAs (dispersion across metros)", years=SAIZ_YEARS,
                            published=sz95["p25"], long_run=sz95["p25"] * f30, role="dispersion"),
        "saiz_p75_95": dict(label="Saiz (2010), 75th percentile of the 95 MSAs (dispersion across metros)", years=SAIZ_YEARS,
                            published=sz95["p75"], long_run=sz95["p75"] * f30, role="dispersion"),
        "saiz_method_2000_10": dict(label="Saiz's method on 2000-10 data (Baum-Snow and Han, Online Appendix D), mean over metros",
                                    years=BSH_YEARS, published=ad["saiz_method_2000_2010_mean"],
                                    published_range=[min(a8), max(a8)], long_run=ad["saiz_method_2000_2010_mean"] * f10,
                                    long_run_range=[min(a8) * f10, max(a8) * f10], role="period check"),
        "saiz_pop95_period_refreshed": dict(label="Saiz's population-weighted mean x (1.25 / 2.54), 10-year equivalent", years=BSH_YEARS,
                                            published=sz95["pop_weighted_mean"] * ratio_obs,
                                            long_run=sz95["pop_weighted_mean"] * ratio_obs * f10, role="period check"),
        "bsh_metro_regions": dict(label="Baum-Snow and Han, metro-region aggregates of the tract elasticities (units, FMM-IV, mean)",
                                  years=BSH_YEARS, published_range=list(ad["bsh_metro_region_means_units"]),
                                  long_run_range=[x * f10 for x in ad["bsh_metro_region_means_units"]], role="reference"),
        "bsh_tracts": dict(label="Baum-Snow and Han, tract elasticities, housing-unit-weighted mean (units, FMM linear)", years=BSH_YEARS,
                           published=bsh_hw, long_run=bsh_hw * f10, role="low"),
    }
    low, central, high = anchors["bsh_tracts"]["long_run"], anchors["saiz_pop95"]["long_run"], anchors["saiz_mean269"]["long_run"]

    def split_eps(a_: float, r_: float) -> tuple[float, float]:
        ec = a_ / (s + (1 - s) * r_)
        return ec, r_ * ec

    cen_c, cen_s = split_eps(central, rho)
    grid = [(a_, r_) for a_ in (low, central, high) for r_ in (rho_lo, rho_hi)]
    cs = [split_eps(a_, r_) for a_, r_ in grid]
    unconv_c, unconv_s = split_eps(sz95["pop_weighted_mean"], rho)
    delta_rows = []
    for dl in DELTAS:
        a_ = sz95["pop_weighted_mean"] * lr_factor(SAIZ_YEARS, dl)
        ec, es = split_eps(a_, rho)
        delta_rows.append(dict(delta=dl, factor_30y=lr_factor(SAIZ_YEARS, dl), factor_10y=lr_factor(BSH_YEARS, dl), metro=a_, city=ec,
                               suburb=es, bsh_tracts_metro=bsh_hw * lr_factor(BSH_YEARS, dl)))
    return dict(
        city_share=s, ratio_suburb_city=rho, ratio_range=[rho_lo, rho_hi],
        ratio_range_matched_split=[float(min(ru + rs)), float(max(ru + rs))], ratio_variants_n=int(vr_ok.sum()),
        ratio_min_variant=v_lo, ratio_max_variant=v_hi,
        ratio_units_mean_by_city_share={k: float(np.mean([rec[f"{sp}_units"]["ratio"] for sp in SPECS])) for k, rec in shares.items()},
        anchors=anchors, anchor_central=central, anchor_low=low, anchor_high=high,
        longrun_factor_30y=f30, longrun_factor_10y=f10, delta=DELTA_CENTRAL,
        longrun_factor_delta35=f10,
        eps_city=dict(value=cen_c, lo=min(x[0] for x in cs), hi=max(x[0] for x in cs)),
        eps_suburb=dict(value=cen_s, lo=min(x[1] for x in cs), hi=max(x[1] for x in cs)),
        eps_metro=dict(value=central, lo=low, hi=high),
        no_horizon_conversion=dict(city=unconv_c, suburb=unconv_s, metro=sz95["pop_weighted_mean"],
                                   note="Saiz's 30-year estimate read as already long-run (the first version of this report)"),
        delta_sensitivity=delta_rows,
        period_check=dict(observed_ratio=ratio_obs, predicted_ratio=ratio_pred, delta=DELTA_CENTRAL,
                          predicted_2000_10_mean=ad["saiz_1970_2000_mean"] * ratio_pred,
                          observed_2000_10_mean=ad["saiz_method_2000_2010_mean"], implied_delta=implied_delta(ratio_obs),
                          table_a8_means_2000_2010=a8, table_a8_means_1970_2000=ad["table_a8_means_1970_2000"]),
        short_run_10yr_units=dict(city=split["FMM_lin_units"]["city"], suburb=split["FMM_lin_units"]["suburb"]),
        short_run_10yr_space=dict(city=split["FMM_lin_space"]["city"], suburb=split["FMM_lin_space"]["suburb"]),
        grid_for_sensitivity=[0.5, 1.0, 2.0, 4.0, "elastic (1e6)"])


def model_definitions(share_aggregate: float, share_mean_of_ratios: float) -> str:
    """What the model must compute to be comparable with the housing inputs (text for the report and the TOML)."""
    return (
        "Preferences (1 - aH) log c + aH log(H - hmin) over non-housing consumption c and housing H; Y' is the household's outlay on c and "
        "housing (income after taxes, education goods, moving, residence-charge and other fixed goods costs), C = Y' - R hmin the "
        "supernumerary expenditure. Then H = hmin + aH C / R, R H = R hmin + aH C, c = (1 - aH) C and Y' = C + R hmin "
        "(DESIGN_v2: the line of DESIGN.md M1b 'Y' = (1 + aH) C + R hmin' double counts aH C; it holds only for preferences "
        "log c + aH log(H - hmin), whose marginal housing share is aH / (1 + aH)). aH is the MARGINAL budget share of housing, the "
        "quantity estimated by the Stone-Geary fit s(y) = aH + (1 - aH) h / y; hmin is h over mean expenditure per equivalent adult "
        "(units of Cbar at R = 1). `housing_share` must be computed as a RATIO OF AGGREGATES, sum_i R H_i / sum_i Y'_i over the adults "
        "(and `housing_share_by_income` as the same ratio within terciles of pre-tax income): the CE targets are aggregate shelter + "
        "utilities over aggregate expenditures. A mean over adults of R H_i / Y'_i (M1B_SPEC) is a different moment: because the share "
        f"falls with expenditure it is larger (the consumer-unit-weighted mean of the five CE quintile shares, {share_mean_of_ratios:.3f}, is a lower bound "
        f"for a household-level mean of ratios, against the aggregate {share_aggregate:.3f}). Units: dollar amounts per household (r_car, the housing-cost gap) are "
        "divided by the mean expenditure per consumer unit, $61,224, as in housing_costs.md, which is unit-free if per-household costs and consumption scale together "
        "with household size; hmin is already a ratio of two amounts in the same unit (h per equivalent adult over mean expenditure per equivalent adult), so "
        "the dollar h per equivalent adult must not be divided by $61,224 (a household of mean size has h = hmin x $61,224).")


def recommendations(res: dict[str, Any]) -> dict[str, Any]:
    rec: dict[str, Any] = {}
    if "a" in res:
        t = res["a"]["tables"]
        g50 = t["all"]["percentiles"]["p50"]["log_gap"]
        rec["gap_housing_cost"] = dict(value=g50["value"], se=g50["boot_se"], n_cz=g50["n"],
                                       basis="log gap of the pooled median monthly housing cost, suburb - "
                                             "city (B25104 bins): reproduces housing_costs.md (0.114, SE 0.017)")
        rec["lower_tail_validation"] = dict(
            distribution="renters, gross rent (B25063)",
            log_gap={k: (t["rent"]["percentiles"][k]["log_gap"]["value"], t["rent"]["percentiles"][k]["log_gap"]["boot_se"]) for k in PCT_KEYS},
            usd_gap={k: (t["rent"]["percentiles"][k]["usd_gap"]["value"], t["rent"]["percentiles"][k]["usd_gap"]["boot_se"]) for k in PCT_KEYS},
            cheap_share_q25=(t["rent"]["cheap_share"]["q25"]["value"], t["rent"]["cheap_share"]["q25"]["boot_se"]),
            basis="gross rent is the price of the unit; the bottom of the all-units and owner distributions is mostly outright owners "
                  "whose cash costs (taxes, insurance, utilities) do not measure the housing consumed")
    if "b" in res:
        q, d = res["b"]["ce"]["quintiles"], res["b"]["ce"]["deciles"]
        eq = res["b"]["equivalized"]
        fq = eq["quintiles"]["fits"][CENTRAL_SCALE]["shelter_utilities"]
        fd = eq["deciles"]["fits"][CENTRAL_SCALE]["shelter_utilities"]
        f_cu, fd_cu = q["fits"]["shelter_utilities"], d["fits"]["shelter_utilities"]
        allfits: list[dict[str, Any]] = []
        defs4 = ("shelter", "shelter_utilities", "shelter_utilities_ex_phone", "housing")
        for tag, blk in (("quintiles", q), ("deciles", d)):
            for tg in ("fits", "fits_dropbottom"):
                allfits += [blk[tg][nm] for nm in defs4]
                for sname in eq[tag]["fits"]:
                    allfits += [eq[tag][tg][sname][nm] for nm in defs4]
        by_scale = {}
        for sname in EQ_SCALES:
            f = f_cu if sname == "none" else eq["quintiles"]["fits"][sname]["shelter_utilities"]
            by_scale[sname] = dict(alphaH=f["alphaH"], hmin=f["h_share_of_mean_exp"], h_usd=f["h_usd"], rmse_pp=f["rmse_pp"],
                                   mean_exp_per_unit=f["mean_exp_used"])
        rec["alphaH"] = dict(value=fq["alphaH"], se=fq["se_ce"]["alphaH"], se_resid=fq["se_resid"]["alphaH"],
                             lo=min(x["alphaH"] for x in allfits), hi=max(x["alphaH"] for x in allfits),
                             decile_fit=fd["alphaH"], cu_level=f_cu["alphaH"], by_scale=by_scale,
                             basis="Stone-Geary NLS on CE 2018 quintile means of shelter + utilities / total expenditures, "
                                   "expenditure per equivalent adult (square-root scale)")
        rec["hmin"] = dict(value=fq["h_share_of_mean_exp"], se=fq["se_ce"]["h_share"], se_resid=fq["se_resid"]["h_share"],
                           lo=min(x["h_share_of_mean_exp"] for x in allfits), hi=max(x["h_share_of_mean_exp"] for x in allfits),
                           decile_fit=fd["h_share_of_mean_exp"], cu_level=f_cu["h_share_of_mean_exp"], h_usd=fq["h_usd"],
                           h_usd_mean_size_cu=fq["h_share_of_mean_exp"] * CE_TOTAL_EXP, mean_exp_per_equivalent_adult=fq["mean_exp_used"],
                           by_scale=by_scale,
                           basis="subsistence housing expenditure per equivalent adult over mean expenditure per equivalent adult "
                                 "(= that of a household of any size if housing and consumption scale together with size): the "
                                 "model's hmin in units of Cbar at R = 1")
        du, qu = d["terciles"]["shelter_utilities"], q["terciles"]["shelter_utilities"]
        size = res["b"]["size_table"]["groups"]
        rec["housing_share"] = dict(
            value=q["defs"]["shelter_utilities"]["all_cu_share"], se=q["defs"]["shelter_utilities"]["all_cu_share_se_rho0"],
            mean_of_group_shares=q["defs"]["shelter_utilities"]["mean_of_group_shares"],
            mean_of_group_shares_deciles=d["defs"]["shelter_utilities"]["mean_of_group_shares"],
            by_income_tercile=du["ratio_of_means"], by_income_tercile_se=du["se_rho0"],
            by_income_tercile_mean_of_group_shares=du["mean_of_group_shares"],
            by_income_tercile_quintiles=qu["ratio_of_means"], by_income_tercile_quintiles_se=qu["se_rho0"],
            one_person_cu=size["One person"]["share"], one_person_cu_se=size["One person"]["share_se_rho0"],
            two_plus_cu=size["Two or more people"]["share"],
            definition="ratio of aggregates: aggregate shelter + utilities / aggregate total expenditures, CE 2018, consumer-unit "
                       "weights (sum R H / sum Y' over adults in the model)",
            basis="shelter + utilities / total expenditures, CE 2018; terciles of income before taxes built from the CE deciles "
                  "(fractional weights; the quintile-based values are in the JSON)")
        rec["model_definitions"] = model_definitions(rec["housing_share"]["value"], rec["housing_share"]["mean_of_group_shares"])
    if "c" in res:
        rec["eps"] = eps_recommendation(res["c"], res["c"]["_city_share"])
    if "d" in res:
        cr = res["d"]["car_requirement"]
        rec["r_car"] = dict(value=cr["central"]["r_car"], se=cr["central"]["se"], lo=cr["low"]["r_car"], hi=cr["high"]["r_car"],
                            extensive_margin=cr["extensive_margin_central"], single_adult=cr["single_adult"]["r_car"],
                            single_adult_se=cr["single_adult"]["se"], first_version=cr["first_version"]["r_car"],
                            unadjusted_upper=cr["unadjusted_upper"]["r_car"], basis=cr["central"]["rule"])
    if "e" in res:
        e = res["e"]
        el = e["elasticity_with_college"]
        va = e["elasticity_variants"]
        slope_measures = [v["elasticity"] for k, v in el.items() if isinstance(v, dict)]
        variant_vals = [v["elasticity"] for k, v in va.items() if isinstance(v, dict)]
        rec["wtp_income_elasticity"] = dict(
            value=e["elasticity_income"]["value"], se=e["elasticity_income"]["se_conditional_on_means"],
            se_with_mean_mwtp=e["elasticity_income"]["se_with_mean_mwtp_independent"],
            with_college=el["linear_projection"]["elasticity"], with_college_se=el["linear_projection"]["se"],
            with_college_lo=min(slope_measures + variant_vals), with_college_hi=max(slope_measures + variant_vals),
            with_college_slope_measures_lo=min(slope_measures), with_college_slope_measures_hi=max(slope_measures),
            with_college_and_children=va["bachelors_children"]["elasticity"], with_college_and_children_se=va["bachelors_children"]["se"],
            associate_plus=va["associate_plus"]["elasticity"], associate_plus_se=va["associate_plus"]["se"])
    return rec


# --------------------------------------------------------------------------
# Markdown report
# --------------------------------------------------------------------------

def table(headers: list[str], rows: list[list[Any]]) -> list[str]:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    out.append("")
    return out


DIST_NAMES = {"all": "all occupied units (B25104)", "rent": "renters, gross rent (B25063)",
              "own": "owners, selected monthly owner costs (B25087)",
              "own_mort": "owners with a mortgage (B25087)", "own_nomort": "owners without a mortgage (B25087)"}


def md_summary(res: dict[str, Any]) -> list[str]:
    L: list[str] = ["## Summary\n"]
    rec = res.get("recommended", {})
    if "a" in res:
        t = res["a"]["tables"]
        tc = res["a"]["tail_composition"]
        P = lambda tab, k: t[tab]["percentiles"][k]["log_gap"]
        D = lambda tab, k: t[tab]["diff_vs_median"][k]
        cs = lambda tab, q: t[tab]["cheap_share"][q]
        bl = tc["below"]["q10"]
        L.append("- **(a) Housing-cost distribution.** The median gap reproduces `housing_costs.md`: "
                 f"**{fse(P('all', 'p50')['value'], P('all', 'p50')['boot_se'], 4)}** log points. All units, log gap "
                 "(SE) at the 10th / 25th / 50th / 75th / 90th percentile: "
                 + " / ".join(fse(P('all', k)['value'], P('all', k)['boot_se']) for k in PCT_KEYS)
                 + f"; p10 minus median {fse(D('all', 'p10')['value'], D('all', 'p10')['boot_se'])}: **pooled over tenure the suburb does not "
                 "lack cheap units**. The pooled distribution mixes three groups with very different cost levels, and its bottom is mostly "
                 f"outright owners, not renters: of the units below the CZ-pooled 10th percentile {pct(bl['own_nomort']['level_1'], 0)}% (city) and "
                 f"{pct(bl['own_nomort']['level_2'], 0)}% (suburb) are owned without a mortgage, {pct(bl['rent']['level_1'], 0)}% and "
                 f"{pct(bl['rent']['level_2'], 0)}% rented; their cash costs (taxes, insurance, utilities) are not the price of the housing "
                 "consumed, so the all-units and owner profiles do not measure R H at the bottom. **Gross rent is the validation moment for "
                 "the lower tail of R_r H_i**: renters' log gap "
                 + " / ".join(fse(P('rent', k)['value'], P('rent', k)['boot_se']) for k in PCT_KEYS)
                 + f" (p10 minus median {fse(D('rent', 'p10')['value'], D('rent', 'p10')['boot_se'])}; dollar gap "
                 f"${fmt(t['rent']['percentiles']['p10']['usd_gap']['value'], 0)} at p10, "
                 f"${fmt(t['rent']['percentiles']['p90']['usd_gap']['value'], 0)} at p90): the suburb lacks cheap rentals, not a uniform "
                 f"fixed charge; rental units below the CZ-pooled 25th percentile {pct(cs('rent', 'q25')['level_1'])}% city vs "
                 f"{pct(cs('rent', 'q25')['level_2'])}% suburb (gap {fse(100 * cs('rent', 'q25')['value'], 100 * cs('rent', 'q25')['boot_se'], 1)} pp). "
                 f"The all-units 90th percentile is partly extrapolated ({sum(t['all']['percentiles']['p90']['n_cells_top_bin'])} of "
                 f"{2 * res['a']['n_cz']} CZ-location cells fall in the open $3,000+ bin): "
                 f"{fmt(t['all']['sensitivity']['top_lower']['p90']['value'], 3)}-{fmt(t['all']['sensitivity']['top_nan']['p90']['value'], 3)} "
                 "across treatments.")
    if "b" in res:
        q = res["b"]["ce"]["quintiles"]
        df = q["defs"]
        eq = res["b"]["equivalized"]["quintiles"]
        f = eq["fits"][CENTRAL_SCALE]["shelter_utilities"]
        f0 = q["fits"]["shelter_utilities"]
        fd = res["b"]["equivalized"]["deciles"]["fits"][CENTRAL_SCALE]["shelter_utilities"]
        acs = res["b"]["acs"]["acs1"]
        hm, al = rec["hmin"], rec["alphaH"]
        sh = lambda nm: " / ".join(pct(x) for x in df[nm]["share"])
        hs = rec["housing_share"]
        tp = res["b"]["equivalized"]["quintiles"]["theta_profile"]["best"]
        tpd = res["b"]["equivalized"]["deciles"]["theta_profile"]["best"]
        L.append("- **(b) Housing share and Stone-Geary fit.** CE 2018 shares of average annual expenditures by quintile of "
                 f"income before taxes (lowest to highest): shelter {sh('shelter')}%; utilities, fuels and public services "
                 f"{sh('utilities')}%; **shelter + utilities {sh('shelter_utilities')}%** (all consumer units "
                 f"{pct(df['shelter_utilities']['all_cu_share'])}%); housing {sh('housing')}% ({pct(df['housing']['all_cu_share'])}%). "
                 "Household size rises with income across the CE groups (people per consumer unit "
                 f"{' / '.join(fmt(x, 1) for x in eq['people'])}), so the fall of the share with expenditure partly reflects economies of "
                 "scale in shared housing, not only nonhomotheticity; the model's household is one adult. **Central fit, expenditure per "
                 f"equivalent adult (square-root scale): aH = {fmt(f['alphaH'])} (SE {fmt(f['se_ce']['alphaH'])}), hmin = "
                 f"{fmt(f['h_share_of_mean_exp'], 4)} (SE {fmt(f['se_ce']['h_share'], 4)}) of mean expenditure per equivalent adult "
                 f"(${fmt(f['mean_exp_used'], 0)}; h = ${fmt(f['h_usd'], 0)} per equivalent adult, ${fmt(hm['h_usd_mean_size_cu'], 0)} for a "
                 f"household of mean size)**, RMSE {fmt(f['rmse_pp'], 2)} pp (deciles {fmt(fd['alphaH'])} and {fmt(fd['h_share_of_mean_exp'], 4)}). "
                 "Other scales (aH; hmin): "
                 + "; ".join(f"{EQ_LABELS[s_].split(' (')[0]} {fmt(v['alphaH'])}, {fmt(v['hmin'], 4)}" for s_, v in al["by_scale"].items())
                 + f" (none = the consumer-unit fit of the first version). The curve fits best at an equivalence "
                 f"elasticity of {fmt(tp['theta'], 2)} (quintiles) / {fmt(tpd['theta'], 2)} (deciles) (two degrees of freedom: suggestive), in line with the square root. "
                 f"Envelope over scales, definitions, quintiles / deciles and dropping the lowest group: aH {fmt(al['lo'])}-{fmt(al['hi'])}, "
                 f"hmin {fmt(hm['lo'], 3)}-{fmt(hm['hi'], 3)}. **The model's `housing_share` must be a ratio of aggregates** (sum R H / sum Y'): "
                 f"**{fse(hs['value'], hs['se'], 3)}** in the data, terciles of income {' / '.join(fmt(x, 3) for x in hs['by_income_tercile'])} "
                 f"(deciles); a mean of household-level ratios is larger (the mean of the five quintile ratios, {fmt(hs['mean_of_group_shares'], 3)}, is a lower bound). "
                 f"ACS cross-check (national, 2018 1-year, household income, no size adjustment): {pct(acs['burden']['share_30plus_all'])}% of households "
                 f"pay 30% or more of income; the same form fitted to mean burdens gives aH = {fmt(acs['fit']['alpha'])} "
                 f"({fmt(min(v['alpha'] for v in acs['sensitivity'].values()))}-{fmt(max(v['alpha'] for v in acs['sensitivity'].values()))} over the open-class assumptions).")
    if "c" in res:
        e = rec["eps"]
        sz = res["c"]["saiz"]
        pc = e["period_check"]
        L.append("- **(c) Supply elasticities.** One horizon convention for every anchor: the long-run stock elasticity "
                 f"eps_T / (1 - (1 - delta)^T) at delta = 3.5% (Baum-Snow and Duranton, slide 76), a factor {fmt(e['longrun_factor_30y'], 2)} for Saiz's "
                 f"30-year and {fmt(e['longrun_factor_10y'], 2)} for the 10-year estimates. Saiz (2010): population-weighted mean "
                 f"{fmt(sz['msas_over_500k']['pop_weighted_mean'], 2)} over the 95 MSAs of Table VI ({fmt(sz['all_msas']['mean'], 2)} unweighted over 269), "
                 f"long run **{fmt(e['anchor_central'], 2)}** ({fmt(e['anchor_low'], 2)} from Baum-Snow and Han's 10-year tract elasticities, "
                 f"{fmt(e['anchor_high'], 2)} from the unweighted mean). Period check (Baum-Snow and Han, Online Appendix D): Saiz's method on 2000-10 data gives a mean "
                 f"of {fmt(pc['observed_2000_10_mean'], 2)} against 2.54 for 1970-2000, a ratio {fmt(pc['observed_ratio'], 2)}; the horizon alone "
                 f"predicts {fmt(pc['predicted_ratio'], 2)} at delta = 3.5% ({fmt(pc['implied_delta'] * 100, 1)}% would match), so the lower 2000-10 estimates are "
                 "mostly the shorter horizon and no further period adjustment is made. Suburban tracts are "
                 f"{fmt(e['ratio_suburb_city'], 2)} times as elastic as central-city tracts "
                 f"(all variants {fmt(e['ratio_range'][0], 2)}-{fmt(e['ratio_range'][1], 2)}: units / floorspace, four estimators, city shares 0.25-0.5, "
                 f"within-region medians). **Recommended long-run eps_H: city {fmt(e['eps_city']['value'], 2)} "
                 f"({fmt(e['eps_city']['lo'], 2)}-{fmt(e['eps_city']['hi'], 2)}), suburb {fmt(e['eps_suburb']['value'], 2)} "
                 f"({fmt(e['eps_suburb']['lo'], 2)}-{fmt(e['eps_suburb']['hi'], 2)})**, metro average "
                 f"{fmt(e['eps_metro']['value'], 2)} ({fmt(e['eps_metro']['lo'], 2)}-{fmt(e['eps_metro']['hi'], 2)}); "
                 f"constructed from the literature, not estimated. Without the horizon conversion (Saiz read as long run, the first version) "
                 f"{fmt(e['no_horizon_conversion']['city'], 2)} / {fmt(e['no_horizon_conversion']['suburb'], 2)}.")
    if "d" in res:
        st = res["d"]["stats"]
        cr = res["d"]["car_requirement"]
        rg = res["d"]["regressions"]
        hsz = res["d"]["household_size"]
        one = SIZE_LABELS[1]
        g_ = lambda k: st[k]
        ev = cr["extra_vehicles"]
        L.append("- **(d) Transport.** Zero-vehicle households: "
                 f"{pct(g_('zero_vehicle_share')['level_1'])}% city vs {pct(g_('zero_vehicle_share')['level_2'])}% suburb "
                 f"(gap {fse(100 * g_('zero_vehicle_share')['value'], 100 * g_('zero_vehicle_share')['boot_se'], 1)} pp); "
                 f"vehicles per household {fmt(g_('vehicles_per_hh')['level_1'], 2)} vs {fmt(g_('vehicles_per_hh')['level_2'], 2)} "
                 f"(gap {fse(g_('vehicles_per_hh')['value'], g_('vehicles_per_hh')['boot_se'])}). Within CZ, controlling for log "
                 "median household income, the suburb coefficient is "
                 f"{fse(rg['veh_per_hh|households|log_income']['beta'][0], rg['veh_per_hh|households|log_income']['se_boot'][0])} "
                 "vehicles, **but suburban households are larger** (one-person households "
                 f"{pct(hsz['composition'][one]['level_1'], 0)}% of city and {pct(hsz['composition'][one]['level_2'], 0)}% of suburban households; 4+ persons "
                 f"{pct(hsz['composition'][SIZE_LABELS[4]]['level_1'], 0)}% vs {pct(hsz['composition'][SIZE_LABELS[4]]['level_2'], 0)}%), and the model's household is one adult. "
                 "Within household-size classes the income-adjusted gaps are "
                 + " / ".join(fmt(hsz['regressions'][f'vehicles|{SIZE_LABELS[k]}|log_income']['beta'][0], 3) for k in SIZE_KS)
                 + f" vehicles (1 / 2 / 3 / 4+ persons); **size-standardized and income-adjusted: {fse(ev['vehicles_size_standardized_income_adjusted']['value'], ev['vehicles_size_standardized_income_adjusted']['se'])}**; "
                 f"one-person households {fse(ev['vehicles_one_person_income_adjusted']['value'], ev['vehicles_one_person_income_adjusted']['se'])}; with household-size shares "
                 f"as controls {fse(ev['vehicles_size_controls_regression']['value'], ev['vehicles_size_controls_regression']['se'])}. "
                 f"**Car requirement as a fixed goods cost: r_car = {fmt(cr['central']['r_car'], 4)} (SE {fmt(cr['central']['se'], 4)}) "
                 f"of mean consumption** (${fmt(cr['central']['r_car'] * CE_TOTAL_EXP, 0)} a year; {fmt(cr['low']['r_car'], 4)} "
                 f"to {fmt(cr['high']['r_car'], 4)}); the first version's {fmt(cr['first_version']['r_car'], 4)} did not adjust for household size "
                 f"({fmt(100 * (cr['central']['r_car'] / cr['first_version']['r_car'] - 1), 0)}%); a single adult reading (one-person households, "
                 f"one-person consumption) gives {fmt(cr['single_adult']['r_car'], 4)}. Vehicles are a normal good (vehicle count elasticity "
                 f"{fmt(res['d']['ce_expenditure_elasticity_quintiles']['vehicles'], 2)} with respect to expenditures across CE income groups).")
    if "e" in res:
        e = res["e"]
        w = rec["wtp_income_elasticity"]
        pm = e["pums90"]
        L.append("- **(e) BFM income gradient.** Elasticity of WTP for school quality with respect to household income at the "
                 f"means: **{fmt(w['value'], 3)} (SE {fmt(w['se'], 3)}; {fmt(w['se_with_mean_mwtp'], 3)} if the SE of the mean "
                 "MWTP is also counted)**. Folding in the college-degree term through the 1990 Bay Area association "
                 f"(PUMS: {pct(pm['share_college_head'])}% of householders with a bachelor's degree or more, {pct(pm['share_assoc_plus_head'])}% with an associate degree or more "
                 f"(BFM's block groups: 44%), correlation with household income {fmt(pm['corr_college_income'], 2)}, slope {fmt(pm['gamma_per_10k'], 4)} per $10,000): "
                 f"**{fmt(w['with_college'], 3)} (SE {fmt(w['with_college_se'], 3)})**; associate-or-more {fmt(w['associate_plus'], 3)}; adding the children term "
                 f"{fmt(w['with_college_and_children'], 3)}; {fmt(w['with_college_lo'], 2)}-{fmt(w['with_college_hi'], 2)} across all variants.")
    L.append("")
    return L


def md_recommended(res: dict[str, Any]) -> list[str]:
    rec = res["recommended"]
    rows = []
    if "gap_housing_cost" in rec:
        r = rec["gap_housing_cost"]
        rows.append(["`gap_housing_cost` (target for a_H)", fse(r["value"], r["se"], 4), "", r["basis"]])
    if "lower_tail_validation" in rec:
        v = rec["lower_tail_validation"]
        rows.append(["lower-tail validation moments for R_r H_i: renters' log gap at p10 / p25 / p50 / p75 / p90",
                     " / ".join(fse(v["log_gap"][k][0], v["log_gap"][k][1]) for k in PCT_KEYS), "",
                     "gross rent, B25063; dollar gap " + " / ".join(fmt(v["usd_gap"][k][0], 0) for k in PCT_KEYS)
                     + "; " + v["basis"]])
    if "alphaH" in rec:
        r, h, s_ = rec["alphaH"], rec["hmin"], rec["housing_share"]
        rows.append(["`alphaH` (marginal budget share of housing)", fse(r["value"], r["se"], 3), f"{fmt(r['lo'], 3)} - {fmt(r['hi'], 3)}",
                     r["basis"] + f"; no size adjustment (consumer-unit level) {fmt(r['cu_level'], 3)}"])
        rows.append(["`hmin` (units of Cbar at R = 1)", fse(h["value"], h["se"], 4), f"{fmt(h['lo'], 3)} - {fmt(h['hi'], 3)}",
                     h["basis"] + f" (h = ${fmt(h['h_usd'], 0)} per equivalent adult of ${fmt(h['mean_exp_per_equivalent_adult'], 0)}; "
                     f"${fmt(h['h_usd_mean_size_cu'], 0)} for a household of mean size); no size adjustment {fmt(h['cu_level'], 4)}"])
        rows.append(["`housing_share` (sum R H / sum Y'), ratio of aggregates", fse(s_["value"], s_["se"], 3),
                     "household-level mean of ratios larger; lower bound " + fmt(s_["mean_of_group_shares"], 3),
                     s_["definition"] + f"; one-person consumer units {fmt(s_['one_person_cu'], 3)}, two or more persons {fmt(s_['two_plus_cu'], 3)} (CE Table 1400)"])
        rows.append(["`housing_share_by_income`, terciles of income before taxes (ratio of aggregates within tercile)",
                     " / ".join(fse(x, se, 3) for x, se in zip(s_["by_income_tercile"], s_["by_income_tercile_se"])), "",
                     "CE deciles with fractional weights (quintile-based: " + " / ".join(fmt(x, 3) for x in s_["by_income_tercile_quintiles"])
                     + "); SE: CE sampling error, upper bound; consumer-unit income terciles, so household-size differences are included"])
    if "eps" in rec:
        e = rec["eps"]
        rows.append(["`epsH` city", fmt(e["eps_city"]["value"], 2), f"{fmt(e['eps_city']['lo'], 2)} - {fmt(e['eps_city']['hi'], 2)}",
                     f"long run (delta = 3.5%): Saiz's population-weighted mean over 95 MSAs x {fmt(e['longrun_factor_30y'], 2)} = {fmt(e['anchor_central'], 2)}, split by the "
                     f"Baum-Snow-Han suburb/city ratio; without the horizon conversion {fmt(e['no_horizon_conversion']['city'], 2)}"])
        rows.append(["`epsH` suburb", fmt(e["eps_suburb"]["value"], 2), f"{fmt(e['eps_suburb']['lo'], 2)} - {fmt(e['eps_suburb']['hi'], 2)}",
                     f"ratio suburb/city {fmt(e['ratio_suburb_city'], 2)} (all variants {fmt(e['ratio_range'][0], 2)}-{fmt(e['ratio_range'][1], 2)}); "
                     f"sensitivity grid 0.5, 1, 2, 4, elastic; without the horizon conversion {fmt(e['no_horizon_conversion']['suburb'], 2)}"])
    if "r_car" in rec:
        r = rec["r_car"]
        rows.append(["`r_car` (fixed goods cost of the suburban car requirement, units of Cbar)", fse(r["value"], r["se"], 4),
                     f"{fmt(r['lo'], 4)} - {fmt(r['hi'], 4)}",
                     r["basis"] + f"; single-adult reading {fmt(r['single_adult'], 4)}; extensive margin only {fmt(r['extensive_margin'], 4)}; "
                     f"first version (no household-size adjustment) {fmt(r['first_version'], 4)}"])
    if "wtp_income_elasticity" in rec:
        w = rec["wtp_income_elasticity"]
        rows.append(["`wtp_income_elasticity`", fse(w["value"], w["se"], 3), "",
                     "BFM Table 8: 1.38 per $10,000 at mean income $54,103 and mean WTP $19.69"])
        rows.append(["`wtp_income_elasticity`, college term folded in", fse(w["with_college"], w["with_college_se"], 3),
                     f"{fmt(w['with_college_lo'], 3)} - {fmt(w['with_college_hi'], 3)}",
                     "plus 13.03 x slope of the householder college share on income (1990 Bay Area PUMS); associate-or-more "
                     f"{fmt(w['associate_plus'], 3)}, with the children term {fmt(w['with_college_and_children'], 3)}"])
    L = ["## Recommended inputs for the model\n"] + table(["Input", "Value (SE)", "Range", "Basis"], rows)
    L += ["How the pieces map into M1 / M1b (DESIGN_v2.md): the supply shifter a_H is fitted to `gap_housing_cost` (the all-units median, as in "
          "`housing_costs.md`); the **renter** percentile profile and renters' units-below-the-cutoff shares of (a) are the validation moments for the lower "
          "tail of R_r H_i (at the bottom the all-units and owner profiles are dominated by outright owners' cash costs, which do not measure housing "
          "consumed); `alphaH`, `hmin` and the tercile `housing_share` come from (b), all for shelter + utilities and, for the first two, in units of "
          "expenditure per equivalent adult (the tercile shares are validation of the fitted curve; they are consumer-unit income terciles and so include "
          "household-size differences); `epsH` is a long-run elasticity for the stationary equilibrium, with the grid of (c) for sensitivity; `r_car` is a "
          "separate suburban fixed cost (ACS housing costs do not contain vehicles, so it does not double count the housing gap), standardized for household "
          "size because suburban households are larger; the WTP elasticity targets the model's `wtp_income_elasticity`.\n"]
    if "model_definitions" in rec:
        L += ["**Model definitions (what the model must compute to be comparable with these inputs).** " + rec["model_definitions"] + "\n"]
    return L


def md_a(res: dict[str, Any]) -> list[str]:
    a = res["a"]
    L = ["## (a) Housing-cost distribution by location\n"]
    L.append("Pooled monthly housing cost (dollars of 2018) by location: bin counts summed over the districts of a "
             "location within each CZ, percentiles by linear interpolation inside the bin (the rule of "
             "`housing_costs.grouped_median`), gap = log(suburb) - log(city) per CZ, then the CZ-weighted mean. "
             "All units: B25104 (renters' gross rent and owners' selected monthly owner costs in 15 bins, top bin "
             "$3,000 or more; units with no cash rent excluded). Renters: B25063 gross rent (24 bins, top $3,500 or "
             "more). Owners: B25087 selected monthly owner costs, with a mortgage (17 bins, top $4,000 or more) and "
             "without (19 bins, top $1,500 or more), pooled on the union of the two grids; the $1,500-or-more no-mortgage "
             "count is spread over the union bins at or above $1,500 like the with-mortgage owners of the partitioned districts (all CZs pooled) "
             f"(shares {', '.join(fmt(s_, 3) for s_ in a['bin_check']['owner_open_bin_allocation']['shares'])} over "
             f"{', '.join(fmt(x, 0) for x in a['bin_check']['owner_open_bin_allocation']['edges'])}). When a percentile "
             "falls in an open top bin it is extrapolated with a Pareto tail through the survival function at the lower edge of the top bin and "
             "of the nearest non-empty bin below it (shape at least 1); the number of CZ-location cells that needed this is in the last "
             "column and the alternatives are in the sensitivity table. Bin check: B25104 bins sum to "
             f"{fmt(a['bin_check']['total_units_B25104_bins'] / 1e6, 2)} million units = renters (B25063) "
             f"{fmt(a['bin_check']['renters_B25063_bins'] / 1e6, 2)} + owners (B25087) "
             f"{fmt(a['bin_check']['owners_B25087_total'] / 1e6, 2)} million.\n")
    rp = a["reproduction"]
    L.append(f"**Reproduction.** The all-units median gap is {fse(rp['value_here'], rp['boot_se_here'], 5)} here and "
             f"{fse(rp.get('value_published'), rp.get('boot_se_published'), 5)} in `housing_costs.json` "
             f"(max CZ-level difference {rp['max_abs_cz_diff_vs_housing_costs']:.1e}); the earnings target `gap_income` "
             f"of `district_settings.py` is {fse(rp.get('earnings_gap_here'), rp.get('earnings_gap_se_here'), 5)} on the same "
             "draws (reported 0.1122, 0.0168).\n")
    for tab in ("all", "rent", "own"):
        d = a["tables"][tab]
        L.append(f"### {DIST_NAMES[tab]}\n")
        rows = []
        for k in PCT_KEYS:
            v = d["percentiles"][k]
            lg = v["log_gap"]
            df = d["diff_vs_median"].get(k)
            rows.append([k[1:], fmt(v["level_1"], 0), fmt(v["level_2"], 0), fse(v["usd_gap"]["value"], v["usd_gap"]["boot_se"], 0),
                         fse(lg["value"], lg["boot_se"]), f"[{fmt(lg['ci95_lo'])}, {fmt(lg['ci95_hi'])}]",
                         "" if df is None else fse(df["value"], df["boot_se"]),
                         f"{v['n_cells_top_bin'][0]} / {v['n_cells_top_bin'][1]}"])
        L += table(["percentile", "city $", "suburb $", "gap $ (SE)", "log gap (SE)", "95% interval", "log gap minus median gap (SE)",
                    "top-bin cells (city / suburb)"], rows)
    d1, d2 = a["tables"]["own_mort"], a["tables"]["own_nomort"]
    L.append("### Owners by mortgage status (log gap, SE)\n")
    L += table(["percentile"] + [DIST_NAMES[t] for t in ("own_mort", "own_nomort")],
               [[k[1:]] + [fse(a["tables"][t]["percentiles"][k]["log_gap"]["value"], a["tables"][t]["percentiles"][k]["log_gap"]["boot_se"])
                           for t in ("own_mort", "own_nomort")] for k in PCT_KEYS])
    L.append("### Units below the CZ-pooled percentile\n")
    L.append("For each CZ the cut is the 10th / 25th / 50th percentile of the two locations pooled; the entry is the share of "
             "each location's units below it (equal shares of 10 / 25 / 50% would mean the same distribution in both "
             "locations). Gap = suburb minus city, in percentage points.\n")
    rows = []
    for tab in ("all", "rent", "own"):
        for q in ("q10", "q25", "q50"):
            v = a["tables"][tab]["cheap_share"][q]
            rows.append([DIST_NAMES[tab], q[1:], pct(v["level_1"]), pct(v["level_2"]), fse(100 * v["value"], 100 * v["boot_se"], 1)])
    L += table(["distribution", "cut (percentile)", "city %", "suburb %", "gap pp (SE)"], rows)
    L.append("### Sensitivity to the top-bin treatment (log gap)\n")
    rows = []
    for tab in ("all", "rent", "own"):
        s_ = a["tables"][tab]["sensitivity"]
        main = a["tables"][tab]["percentiles"]
        rows.append([DIST_NAMES[tab], "Pareto tail (main)"] + [fmt(main[k]["log_gap"]["value"]) for k in PCT_KEYS] + [str(main["p90"]["log_gap"]["n"])])
        rows.append(["", "CZs with an undetermined percentile dropped"] + [fmt(s_["top_nan"][k]["value"]) for k in PCT_KEYS]
                    + [f"{s_['top_nan']['p90']['n']} (p90)"])
        rows.append(["", "percentile set at the bin's lower edge"] + [fmt(s_["top_lower"][k]["value"]) for k in PCT_KEYS]
                    + [str(s_["top_lower"]["p90"]["n"])])
        rows.append(["", f"common sample, no extrapolated cell ({s_['common_no_extrapolation_n_cz']} CZs)"]
                    + [fmt(s_["common_no_extrapolation"][k]["value"]) for k in PCT_KEYS] + [str(s_["common_no_extrapolation_n_cz"])])
        if tab == "all" and "all_fine" in a["tables"]:
            fp = a["tables"]["all_fine"]["percentiles"]
            rows.append(["", "finer grid: renters (B25063) + owners (B25087) pooled, Pareto tail above $4,000"]
                        + [fmt(fp[k]["log_gap"]["value"]) for k in PCT_KEYS] + [str(fp["p90"]["log_gap"]["n"])])
    L += table(["distribution", "treatment", "p10", "p25", "p50", "p75", "p90", "CZs"], rows)
    tc = a["tail_composition"]
    L.append("### Who occupies the cheapest units\n")
    L.append("The pooled distribution mixes three groups with very different cost levels, and their mix differs by location. For each CZ the cut is "
             "the 10th / 25th / 50th percentile of the all-units distribution of the two locations pooled; the entries are the shares of the units "
             "below the cut that are rented (cash rent, B25063), owned with a mortgage or owned without a mortgage (B25087), each group's count below "
             "the cut taken by linear interpolation inside its own bins; the occupancy rows are shares of occupied units (cash renters plus owners). "
             "CZ-weighted means, whole-CZ bootstrap, gap = suburb minus city in percentage points.\n")
    rows = []
    for k, lab in TENURE_GROUPS:
        v = tc["occupancy"][k]
        rows.append([lab, "all occupied units", pct(v["level_1"]), pct(v["level_2"]), fse(100 * v["value"], 100 * v["boot_se"], 1)])
    for q in ("q10", "q25", "q50"):
        for k, lab in TENURE_GROUPS:
            v = tc["below"][q][k]
            rows.append([lab, f"units below the pooled {q[1:]}th percentile (about ${fmt(tc['cut_usd'][q]['value'], 0)})", pct(v["level_1"]),
                         pct(v["level_2"]), fse(100 * v["value"], 100 * v["boot_se"], 1)])
    L += table(["group", "population", "city %", "suburb %", "gap pp (SE)"], rows)
    L.append("Monthly cost by group, city / suburb in dollars (CZ-weighted means of the percentiles):\n")
    rows = []
    for tab, lab in (("rent", "renters (B25063)"), ("own_mort", "owners with a mortgage (B25087)"), ("own_nomort", "owners without a mortgage (B25087)")):
        dd = a["tables"][tab]["percentiles"]
        rows.append([lab] + [f"{fmt(dd[k]['level_1'], 0)} / {fmt(dd[k]['level_2'], 0)}" for k in PCT_KEYS])
    L += table(["group", "p10", "p25", "p50", "p75", "p90"], rows)
    own_nm, rent = tc["below"]["q10"]["own_nomort"], tc["below"]["q10"]["rent"]
    own_nm25, rent25 = tc["below"]["q25"]["own_nomort"], tc["below"]["q25"]["rent"]
    occ_r, occ_m, occ_n = tc["occupancy"]["rent"], tc["occupancy"]["own_mort"], tc["occupancy"]["own_nomort"]
    rp = a["tables"]["rent"]["percentiles"]
    L.append("**Reading.** For the pooled distribution (all units) the log gap is about the same at the 10th, 50th and 90th percentile and lower "
             "at the 25th: no bottom-heavy gap. But the pooled gap is not a clean measure of the price of housing at the bottom: the cheapest units "
             "are mostly owned outright, not rented. Of the units below the CZ-pooled 10th percentile, "
             f"{pct(own_nm['level_1'], 0)}% (city) and {pct(own_nm['level_2'], 0)}% (suburb) are owned without a mortgage and "
             f"{pct(rent['level_1'], 0)}% and {pct(rent['level_2'], 0)}% are rented (below the 25th percentile: "
             f"{pct(own_nm25['level_1'], 0)}% / {pct(own_nm25['level_2'], 0)}% and {pct(rent25['level_1'], 0)}% / {pct(rent25['level_2'], 0)}%); "
             "mortgage owners are a few percent of the bottom. The suburb has more outright owners "
             f"({pct(occ_n['level_2'], 0)}% of occupied units against {pct(occ_n['level_1'], 0)}% in the city) and more mortgage owners "
             f"({pct(occ_m['level_2'], 0)}% against {pct(occ_m['level_1'], 0)}%) but far fewer renters ({pct(occ_r['level_2'], 0)}% against "
             f"{pct(occ_r['level_1'], 0)}%), so a gap between pooled percentiles combines within-group price differences with composition. The "
             "cash cost of an outright owner is property tax, insurance and utilities, not the price of the housing consumed, so the bottom of the "
             "all-units and owner distributions does not measure R_r H_i. **Gross rent does** (a rent is the price of the unit): the renter log "
             f"gap falls from {fmt(rp['p10']['log_gap']['value'], 3)} at the 10th to {fmt(rp['p90']['log_gap']['value'], 3)} at the 90th percentile "
             f"and in dollars from ${fmt(rp['p10']['usd_gap']['value'], 0)} to ${fmt(rp['p90']['usd_gap']['value'], 0)}, so the suburb lacks cheap "
             "rentals rather than facing a uniform fixed charge. **The renter (gross rent) profile and renters' units-below-the-cutoff shares are the "
             "validation moments for the lower tail of R_r H_i**. The all-units median gap remains the target for the supply shifter (it is the "
             "`housing_costs.md` target), but it mixes the tenure groups too (units owned without a mortgage are "
             f"{pct(tc['below']['q50']['own_nomort']['level_1'], 0)}% of the city's and {pct(tc['below']['q50']['own_nomort']['level_2'], 0)}% of the suburb's units below the "
             f"pooled median): the like-for-like renter median gap is {fse(rp['p50']['log_gap']['value'], rp['p50']['log_gap']['boot_se'])}, about half of "
             f"{fmt(a['tables']['all']['percentiles']['p50']['log_gap']['value'], 3)}. The all-units and owner profiles are descriptive only.\n")
    return L


def by_theta(eqb: dict[str, Any], theta: float) -> dict[str, Any]:
    """Row of the theta profile closest to `theta`."""
    return min(eqb["theta_profile"]["rows"], key=lambda r_: abs(r_["theta"] - theta))


def md_b(res: dict[str, Any]) -> list[str]:
    b = res["b"]
    ce = b["ce"]
    q, d = ce["quintiles"], ce["deciles"]
    L = ["## (b) Housing expenditure shares by income and the Stone-Geary fit\n"]
    L.append("Source: BLS Consumer Expenditure Survey, calendar year 2018, **Table 1101** (quintiles of income before taxes) and "
             "**Table 1110** (deciles), 'Annual expenditure means, shares, standard errors, and coefficients of variation' "
             "(https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error/"
             "cu-income-quintiles-before-taxes-2018.xlsx and ...-deciles-before-taxes-2018.xlsx; PDFs cached beside them). "
             "Shares are item means over mean total expenditures of the group (they reproduce the published one-decimal shares). "
             "Consumer units (CUs), not adults: see the household-size subsection; dollars of 2018.\n")
    L.append("### Quintiles of income before taxes\n")
    rows = []
    for i, nm in enumerate(q["columns"]):
        df = q["defs"]
        rows.append([nm, fmt(q["n_cu_thousands"][i], 0), fmt(q["mean_income_before_taxes"][i], 0), fmt(q["mean_expenditures"][i], 0),
                     fmt(df["shelter"]["mean"][i], 0), pct(df["shelter"]["share"][i]), fmt(df["utilities"]["mean"][i], 0),
                     pct(df["utilities"]["share"][i]), fmt(df["shelter_utilities"]["mean"][i], 0),
                     f"{pct(df['shelter_utilities']['share'][i])} ({pct(df['shelter_utilities']['share_se_rho0'][i])})",
                     pct(df["shelter_utilities_ex_phone"]["share"][i]), fmt(df["housing"]["mean"][i], 0), pct(df["housing"]["share"][i])])
    df = q["defs"]
    rows.append(["All consumer units", fmt(q["all_cu"]["n_cu_thousands"], 0), "78,635", fmt(q["all_cu"]["mean_expenditures"], 0),
                 fmt(df["shelter"]["all_cu_mean"], 0), pct(df["shelter"]["all_cu_share"]), fmt(df["utilities"]["all_cu_mean"], 0),
                 pct(df["utilities"]["all_cu_share"]), fmt(df["shelter_utilities"]["all_cu_mean"], 0), pct(df["shelter_utilities"]["all_cu_share"]),
                 pct(df["shelter_utilities_ex_phone"]["all_cu_share"]), fmt(df["housing"]["all_cu_mean"], 0), pct(df["housing"]["all_cu_share"])])
    rows[-1][9] = f"{pct(df['shelter_utilities']['all_cu_share'])} ({pct(df['shelter_utilities']['all_cu_share_se_rho0'])})"
    L += table(["group", "CUs (000)", "income before taxes $", "expenditures $", "shelter $", "shelter %", "utilities $", "utilities %",
                "shelter + utilities $", "shelter + utilities % (SE)", "shelter + utilities ex. telephone %", "housing $", "housing %"], rows)
    mx = max(v["max_abs_diff_vs_published_pp"] for blk in (q, d) for v in blk["defs"].values() if v.get("max_abs_diff_vs_published_pp") is not None)
    L.append("The SE of a share (percentage points, in parentheses) is the delta-method SE of item mean over total with the published SEs of the "
             "CE means and independent sampling errors (an upper bound; with a correlation of 0.5 between item and total the SE of the "
             f"shelter + utilities share is {pct(df['shelter_utilities']['share_se_rho05'][0])} pp for the lowest quintile). The computed shares "
             f"reproduce the published one-decimal shares to within {fmt(mx, 2)} pp. The shelter + utilities SE treats the two items as independent. "
             "Mean of the five group shares "
             f"(a mean-of-ratios counterpart weighted by consumer units): shelter + utilities "
             f"{pct(df['shelter_utilities']['mean_of_group_shares'])}%.\n")
    L.append("CE 'Utilities, fuels, and public services' includes telephone services (about a third of it); 'shelter + utilities ex. telephone' "
             "(shelter, natural gas, electricity, fuel oil and other fuels, water and other public services) covers what an ACS gross rent or "
             "owner cost covers and is the closest match to the ACS-based `gap_housing_cost` target; the requested shelter + utilities "
             "remains the central definition.\n")
    L.append("### Deciles of income before taxes (shares of expenditures, %)\n")
    dd = d["defs"]
    rows = [[nm, fmt(d["mean_income_before_taxes"][i], 0), fmt(d["mean_expenditures"][i], 0), pct(dd["shelter"]["share"][i]),
             pct(dd["utilities"]["share"][i]), pct(dd["shelter_utilities"]["share"][i]), pct(dd["shelter_utilities_ex_phone"]["share"][i]),
             pct(dd["housing"]["share"][i])] for i, nm in enumerate(d["columns"])]
    L += table(["group", "income before taxes $", "expenditures $", "shelter", "utilities", "shelter + utilities", "shelter + utilities ex. telephone", "housing"], rows)
    L.append("### Stone-Geary fit on consumer-unit expenditure (no household-size adjustment): s(Y) = h/Y + aH (1 - h/Y), Y = total expenditures of the group\n")
    L.append("Nonlinear least squares in (h, aH) on the group shares, equal weight per group. The form is linear in (aH, "
             "(1 - aH) h) after writing s = aH + (1 - aH) h / Y, so the NLS solution equals that regression; a Levenberg-Marquardt "
             "fit from a different start reproduces it (largest relative difference "
             f"{max(f['nls_check_max_abs_diff'] for blk in (q, d) for f in blk['fits'].values()):.0e}). h/mean = h over CE mean "
             f"expenditures per consumer unit (${CE_TOTAL_EXP:,.2f}). SE (CE) = parametric bootstrap of the CE sampling error "
             "(item mean and total drawn normal with the published SEs, independent across groups and of each other: "
             "conservative); SE (resid) = delta method from the residual variance (n - 2 degrees of freedom). "
             "Elasticity = aH / s(mean): the elasticity of housing expenditure with respect to total expenditures at the "
             "mean.\n")
    rows = []
    for tag, blk in (("quintiles", q), ("deciles", d)):
        for nm in HOUSING_DEFS:
            f = blk["fits"][nm]
            rows.append([tag, DEF_LABELS[nm], fse(f["alphaH"], f["se_ce"]["alphaH"]), fmt(f["se_resid"]["alphaH"]), fmt(f["h_usd"], 0),
                         fse(f["h_share_of_mean_exp"], f["se_ce"]["h_share"], 4), fmt(f["se_resid"]["h_share"], 4),
                         pct(f["s_at_mean_exp"]), fmt(f["elasticity_at_mean"], 2), fmt(f["rmse_pp"], 2), fmt(f["max_abs_resid_pp"], 2), fmt(f["r2"], 3)])
    L += table(["groups", "definition", "aH (SE, CE)", "SE resid", "h $", "h / mean exp. (SE, CE)", "SE resid", "s at mean %",
                "elasticity at mean", "RMSE pp", "max abs resid pp", "R2"], rows)
    f = q["fits"]["shelter_utilities"]
    L.append("Fit of the central definition (shelter + utilities, quintiles):\n")
    L += table(["group", "expenditures $", "data share %", "fitted share %", "residual pp", "subsistence h/Y %"],
               [[nm, fmt(f["Y"][i], 0), pct(f["shares"][i], 2), pct(f["fitted"][i], 2), fmt(f["resid_pp"][i], 2), pct(f["subsistence_share_by_group"][i], 1)]
                for i, nm in enumerate(q["columns"])])
    L.append("Robustness (aH; h / mean expenditures):\n")
    rows = []
    for nm in ("shelter", "shelter_utilities", "shelter_utilities_ex_phone", "housing"):
        rows.append([DEF_LABELS[nm]] + [f"{fmt(blk[tag][nm]['alphaH'])}; {fmt(blk[tag][nm]['h_share_of_mean_exp'], 4)}"
                                                for blk in (q, d) for tag in ("fits", "fits_dropbottom", "fits_net_of_pensions")])
    L += table(["definition", "quintiles", "quintiles, lowest dropped", "quintiles, net of pensions", "deciles", "deciles, lowest dropped", "deciles, net of pensions"], rows)
    L.append(f"'Net of pensions' deflates by expenditures excluding personal insurance and pensions (mean ${fmt(q['net_pensions_mean_exp'], 0)}), "
             "a closer match to consumption (these entries are in units of that smaller mean). Dropping the lowest group removes "
             "consumer units whose expenditures exceed income (students, retirees, business losses).\n")
    L.append("### Household size and the income gradient: fits on equivalized expenditure\n")
    eqq, eqd = b["equivalized"]["quintiles"], b["equivalized"]["deciles"]
    sc = eqq["scales"]
    L.append("The CE counts consumer units (CUs), and household size rises with income: people per CU "
             f"{' / '.join(fmt(x, 1) for x in eqq['people'])} and earners {' / '.join(fmt(x, 1) for x in eqq['earners'])} from the lowest to the "
             "highest quintile. The fall of the housing share with expenditure therefore mixes nonhomotheticity with economies of scale in shared "
             "housing, and the model's household is one adult. With an equivalence scale e = n^theta the Stone-Geary share of a household is "
             "s = aH + (1 - aH) h1 / y with y = Y / e (equivalized expenditure) and h1 the subsistence housing expenditure per equivalent adult. "
             "h1 over the mean expenditure per equivalent adult (CE mean expenditures per CU $61,224.13 over the CU-weighted mean scale) is the "
             "subsistence share of a household of ANY size if housing and consumption scale together with household size, the unit-free assumption "
             "of `housing_costs.md`; it is the model's hmin in units of Cbar. Scales (theta = 0 no adjustment, 1 per capita; Foerster and Mira "
             "d'Ercole (OECD), Table 2): square root of persons (the OECD income-distribution convention, elasticity 0.5; central); OECD-modified "
             "(1 + 0.5 (adults - 1) + 0.3 children, elasticity about 0.53; adults = persons minus children under 18, the CE's child age limit "
             "instead of the OECD's 14); per adult (children count zero: the model's literal 'one household per adult', DESIGN_v2 mean expenditure "
             "about $32,000); per capita. People per CU are published to one decimal and treated as exact. Shelter + utilities, equal weight on the groups, "
             "SE (CE) as above.\n")
    rows = []
    for i, nm in enumerate(q["columns"]):
        rows.append([nm, fmt(eqq["people"][i], 1), fmt(eqq["children"][i], 1), fmt(eqq["adults"][i], 1), fmt(eqq["earners"][i], 1)]
                    + [fmt(q["mean_expenditures"][i] / sc[s_]["values"][i], 0) for s_ in ("sqrt", "oecd_modified", "per_adult", "per_capita")])
    rows.append(["mean scale (CU-weighted) / mean expenditure per equivalent adult $", "", "", "", ""]
                + [f"{fmt(sc[s_]['mean_scale'], 3)} / {fmt(sc[s_]['mean_exp_per_equivalent_adult'], 0)}" for s_ in ("sqrt", "oecd_modified", "per_adult", "per_capita")])
    L += table(["group", "people", "children under 18", "adults", "earners", "expenditure per sqrt-equivalent $", "per OECD-modified $",
                "per adult $", "per capita $"], rows)
    L.append("Central definition (shelter + utilities), Stone-Geary fit by scale:\n")
    rows = []
    for tag, blk, eqb in (("quintiles", q, eqq), ("deciles", d, eqd)):
        for s_ in EQ_SCALES:
            f = blk["fits"]["shelter_utilities"] if s_ == "none" else eqb["fits"][s_]["shelter_utilities"]
            rows.append([tag, EQ_LABELS[s_] + (" (central)" if s_ == CENTRAL_SCALE else ""), fse(f["alphaH"], f["se_ce"]["alphaH"]), fmt(f["h_usd"], 0),
                         fmt(f["mean_exp_used"], 0), fse(f["h_share_of_mean_exp"], f["se_ce"]["h_share"], 4), pct(f["s_at_mean_exp"]),
                         fmt(f["elasticity_at_mean"], 2), fmt(f["rmse_pp"], 2), fmt(f["r2"], 3)])
    L += table(["groups", "scale", "aH (SE, CE)", "h $ per equivalent adult", "mean expenditure per equivalent adult $", "h / mean (SE, CE)",
                "s at mean %", "elasticity at mean", "RMSE pp", "R2"], rows)
    L.append("All housing definitions by scale (aH; h / mean expenditure per equivalent adult):\n")
    rows = []
    for tag, blk, eqb in (("quintiles", q, eqq), ("deciles", d, eqd)):
        for nm in ("shelter", "shelter_utilities", "shelter_utilities_ex_phone", "housing"):
            rows.append([tag, DEF_LABELS[nm]] + [f"{fmt((blk['fits'][nm] if s_ == 'none' else eqb['fits'][s_][nm])['alphaH'])}; "
                                                f"{fmt((blk['fits'][nm] if s_ == 'none' else eqb['fits'][s_][nm])['h_share_of_mean_exp'], 4)}"
                                                for s_ in EQ_SCALES])
    L += table(["groups", "definition"] + [EQ_LABELS[s_].split(" (")[0] for s_ in EQ_SCALES], rows)
    L.append("Fit by equivalence elasticity theta (persons as the size measure; shelter + utilities; aH; h / mean expenditure per equivalent adult; RMSE in pp):\n")
    rows = []
    for tag, eqb in (("quintiles", eqq), ("deciles", eqd)):
        tp = eqb["theta_profile"]
        pick = {r_["theta"]: r_ for r_ in tp["rows"]}
        rows.append([tag] + [f"{fmt(pick[th]['alphaH'])}; {fmt(pick[th]['h_share'], 3)}; {fmt(pick[th]['rmse_pp'], 2)}" for th in (0.0, 0.25, 0.5, 0.75, 1.0)]
                    + [f"theta = {fmt(tp['best']['theta'], 2)}: {fmt(tp['best']['alphaH'])}; {fmt(tp['best']['h_share'], 3)}; {fmt(tp['best']['rmse_pp'], 2)}"])
    L += table(["groups", "theta = 0", "0.25", "0.5", "0.75", "1", "best fit on a 0.01 grid"], rows)
    f = eqq["fits"][CENTRAL_SCALE]["shelter_utilities"]
    L.append("Fit of the central equivalized specification (shelter + utilities, quintiles, square-root scale):\n")
    L += table(["group", "equivalized expenditure $", "data share %", "fitted share %", "residual pp", "subsistence h/y %"],
               [[nm, fmt(f["Y"][i], 0), pct(f["shares"][i], 2), pct(f["fitted"][i], 2), fmt(f["resid_pp"][i], 2), pct(f["subsistence_share_by_group"][i], 1)]
                for i, nm in enumerate(q["columns"])])
    hm, al = res["recommended"]["hmin"], res["recommended"]["alphaH"]
    tpq, tpd = eqq["theta_profile"]["best"], eqd["theta_profile"]["best"]
    L.append(f"**Central choice: the square-root scale, aH = {fmt(al['value'])} and hmin = {fmt(hm['value'], 4)}.** (i) It is the OECD's "
             "income-distribution convention (OECD-modified, elasticity about 0.53, gives "
             f"{fmt(al['by_scale']['oecd_modified']['alphaH'])} and {fmt(al['by_scale']['oecd_modified']['hmin'], 3)}); (ii) the Stone-Geary curve fits the group means best at an "
             f"equivalence elasticity of {fmt(tpq['theta'], 2)} for the quintiles (RMSE {fmt(tpq['rmse_pp'], 2)} pp, against "
             f"{fmt(by_theta(eqq, 0.0)['rmse_pp'], 2)} pp with no adjustment and {fmt(by_theta(eqq, 1.0)['rmse_pp'], 2)} pp per capita) and "
             f"{fmt(tpd['theta'], 2)} for the deciles (two degrees of freedom: suggestive, not an estimate of the scale); (iii) no adjustment and per capita bracket it, and "
             "they are the ends of the range together with the other definitions. In dollars a household of mean size has "
             f"h = ${fmt(hm['h_usd_mean_size_cu'], 0)} a year (${fmt(hm['h_usd_mean_size_cu'] / 12, 0)} a month; ${fmt(hm['h_usd'], 0)} per equivalent adult "
             f"on mean expenditure per equivalent adult of ${fmt(hm['mean_exp_per_equivalent_adult'], 0)}), against ${fmt(q['fits']['shelter_utilities']['h_usd'], 0)} with no adjustment. "
             "The marginal share aH does not depend on the unit of the model, only on the adjustment for size: it falls from "
             f"{fmt(al['by_scale']['none']['alphaH'])} (no adjustment) to {fmt(al['by_scale']['per_capita']['alphaH'])} (per capita). Per adult, the model's literal unit, "
             f"gives {fmt(al['by_scale']['per_adult']['alphaH'])} and {fmt(al['by_scale']['per_adult']['hmin'], 3)}. The 'mean expenditure per equivalent adult' denominator is aggregate "
             "expenditure over aggregate equivalent adults; the CU-weighted mean of the groups' equivalized expenditures would give a slightly larger "
             f"hmin ({fmt(hm['h_usd'] / eqq['scales'][CENTRAL_SCALE]['mean_of_group_equivalized_exp'], 3)} at the square-root scale).\n")
    sck = b["size_check"]
    o, t2 = sck["One person"], sck["Two or more people"]
    L.append("**Cross-check with CE Table 1400 (size of consumer unit).** One-person CUs (expenditures "
             f"${fmt(o['expenditures'], 0)}) spend {pct(o['share'])}% (SE {fmt(100 * o['share_se'], 1)} pp) of expenditures on shelter + utilities, CUs of two or more people "
             f"(${fmt(t2['expenditures'], 0)}, {fmt(t2['people'], 1)} persons) {pct(t2['share'])}%. The CU-level income curve predicts {pct(o['predicted_cu_level'])}% and "
             f"{pct(t2['predicted_cu_level'])}%, the square-root curve {pct(o['predicted_sqrt_scale'])}% and {pct(t2['predicted_sqrt_scale'])}%: one-person CUs spend "
             f"{fmt(o['excess_over_sqrt_scale_pp'], 1)} pp more than the equivalized curve implies, i.e. economies of scale in housing are stronger than the square root "
             "(a second person adds little to subsistence housing). With the central aH the one-person share implies a subsistence share of "
             f"{fmt(o['implied_subsistence_share_given_aH_sqrt'], 3)} of their expenditure (${fmt(o['implied_subsistence_usd'], 0)} a year) against "
             f"{fmt(t2['implied_subsistence_share_given_aH_sqrt'], 3)} (${fmt(t2['implied_subsistence_usd'], 0)}) for CUs of two or more. Read literally for a model of single adults "
             f"the subsistence share would be nearer {fmt(o['implied_subsistence_share_given_aH_sqrt'], 2)} than {fmt(hm['value'], 3)}, and the aggregate share nearer "
             f"{fmt(o['share'], 3)} than {fmt(res['recommended']['housing_share']['value'], 3)}. This is a one-point cross-check (one-person CUs are older and more often "
             "retired), not a fit; it is reported, not included in the range.\n")
    L.append("### Terciles of income before taxes (for `housing_share_by_income`)\n")
    L.append("Terciles are built from the equal-sized groups with fractional weights (from the deciles: tercile 1 = the lowest three deciles plus "
             "one third of the fourth, and so on; from the quintiles: tercile 1 = the lowest quintile plus two thirds of the second). The "
             "decile-based terciles are closer to true tercile cuts and are the recommended `housing_share_by_income` values. Ratio of means "
             "(aggregate dollars over aggregate expenditures, the moment the model must compute) in %, with the parametric-bootstrap "
             "SE of the CE sampling error (independent errors, an upper bound) in parentheses; the last three columns are the consumer-unit-weighted "
             "means of the groups' ratios (a mean-of-ratios counterpart). These are terciles of consumer-unit income before taxes, so they include "
             "household-size differences (people per CU rise from 1.7 to 3.2 across the quintiles), unlike the equivalized fit above:\n")
    rows = []
    for tag, blk in (("quintiles", q), ("deciles", d)):
        t = blk["terciles"]
        for nm in HOUSING_DEFS:
            rows.append([tag, DEF_LABELS[nm]] + [f"{pct(x)} ({pct(se)})" for x, se in zip(t[nm]["ratio_of_means"], t[nm]["se_rho0"])]
                        + [pct(x) for x in t[nm]["mean_of_group_shares"]])
        rows.append([tag, "mean expenditures $"] + [fmt(x, 0) for x in t["mean_expenditures"]] + ["", "", ""])
    L += table(["groups", "definition", "tercile 1 (ratio of means)", "tercile 2", "tercile 3", "T1 (mean of group shares)", "T2", "T3"], rows)
    acs = b["acs"]
    L.append("### Cross-check: ACS national housing cost burdens by household income\n")
    L.append("ACS 2018 1-year, United States (tables B25074 renters' gross rent and B25095 owners' selected monthly owner costs as a "
             "percentage of household income, households with a computed burden, B25106 for the coarse 30%-or-more shares; mean "
             "income by bin from B19001 midpoints and, for $100,000 and over, from aggregate income B19025). The 5-year 2014-18 file "
             "is in the JSON and agrees to about 0.01. Mean burden uses class midpoints, with 12% for 'less than 20%' and 70% for "
             "'50% or more'; the median is interpolated inside classes (shown when it lies in a bounded class). ACS costs include "
             "utilities and, for owners, mortgage principal, taxes and insurance (not maintenance); CE shelter excludes principal "
             "and CE utilities include telephone service. Burdens are household-level ratios to household income with no adjustment for household "
             "size, so the comparison is with the consumer-unit-level (no adjustment) CE fit.\n")
    bl = acs["acs1"]["burden"]
    rows = [[(f"${r['income_bin'][0]:,.0f} and over" if r["income_bin"][1] is None else f"${r['income_bin'][0]:,.0f} - ${r['income_bin'][1]:,.0f}"), fmt(r["households"] / 1e6, 2),
             fmt(r["mean_income"], 0), pct(r["share_30plus"]), pct(r["share_50plus"]), pct(r["mean_burden"]),
             "" if r["median_in_open_class"] or not np.isfinite(r["median_burden"]) else pct(r["median_burden"])] for r in bl["rows"]]
    L += table(["household income", "households (M)", "mean income $", "burden >= 30% (%)", "burden >= 50% (%)", "mean burden (%)", "median burden (%)"], rows)
    L.append("B25106 (coarse bins): share paying 30% or more of income: " + "; ".join(f"{c['income_bin']} {pct(c['share_30plus'])}% (90% margin of error "
             f"+-{pct(c['share_30plus_moe90'], 2)} pp)" for c in bl["b25106"])
             + f". All households with a computed burden: {pct(bl['share_30plus_all'])}% (B25106: {pct(bl['b25106_all']['share_30plus'])}%, margin of error "
             f"+-{pct(bl['b25106_all']['moe90'], 2)} pp): national ACS sampling error is negligible here, the uncertainty is the class values for the open "
             "classes.\n")
    fa = acs["acs1"]["fit"]
    qi = q["mean_income_before_taxes"]
    sh_u = df["shelter_utilities"]["mean"]
    ce_burden = [sh_u[i] / qi[i] for i in range(5)]
    L.append(f"Stone-Geary form on income, s(y) = a + b / y, fitted to the mean burdens (weights = households): "
             f"a = {fmt(fa['alpha'])}, h = ${fmt(fa['h_usd'], 0)} = {fmt(fa['h_share_of_mean_income'], 3)} of mean household income "
             f"(${fmt(fa['mean_income'], 0)}), RMSE {fmt(fa['rmse_pp'], 1)} pp (the lowest bins are fitted poorly); a ranges "
             f"{fmt(min(v['alpha'] for v in acs['acs1']['sensitivity'].values()))}-{fmt(max(v['alpha'] for v in acs['acs1']['sensitivity'].values()))} and "
             f"h/mean income {fmt(min(v['h_share_of_mean_income'] for v in acs['acs1']['sensitivity'].values()), 3)}-"
             f"{fmt(max(v['h_share_of_mean_income'] for v in acs['acs1']['sensitivity'].values()), 3)} when the class values are "
             "varied (less than 20%: 10-16%; 50% or more: 60-80%). In dollars h is "
             f"{fmt(fa['h_usd'] / CE_TOTAL_EXP, 3)} of CE mean expenditures, against {fmt(q['fits']['shelter_utilities']['h_share_of_mean_exp'], 3)} "
             "from the CE fit. CE shelter + utilities over income before taxes, by quintile: "
             + " / ".join(pct(x) for x in ce_burden) + "% (mean incomes " + " / ".join(fmt(x, 0) for x in qi)
             + "), against ACS mean burdens of " + " / ".join(pct(r["mean_burden"]) for r in bl["rows"][2:]) + "% at mean incomes "
             + " / ".join(fmt(r["mean_income"], 0) for r in bl["rows"][2:]) + " for the bins from $20,000 up.\n")
    return L


def md_c(res: dict[str, Any]) -> list[str]:
    c = res["c"]
    e = res["recommended"]["eps"]
    sz, bsh = c["saiz"], c["baum_snow_han"]
    L = ["## (c) Housing supply elasticities\n"]
    L.append("### Saiz (2010)\n")
    L.append("A. Saiz, 'The Geographic Determinants of Housing Supply', *Quarterly Journal of Economics* 125(3): 1253-1296. "
             "Price elasticity of the supply of housing units by metropolitan area, estimated from 1970-2000 changes (a 30-year horizon) with "
             "instrumented demand shocks interacted with the share of land unavailable for development (geography) and the "
             "Wharton regulation index (Table V column 2; estimates in Table VI, pp. 1282-1284, for the 95 MSAs with population over "
             "500,000 in 2000). Text (p. 1281): population-weighted average 1.75, unweighted 2.5. The author's file "
             "(269 MSAs; downloaded from the link on urbaneconomics.mit.edu) reproduces both and Table VI:\n")
    rows = []
    for lab, k in (("95 MSAs over 500,000 (Table VI)", "msas_over_500k"), ("all 269 MSAs", "all_msas")):
        s_ = sz[k]
        rows.append([lab, s_["n"], fmt(s_["mean"], 2), fmt(s_["pop_weighted_mean"], 2), fmt(s_["median"], 2),
                     f"{fmt(s_['p25'], 2)} - {fmt(s_['p75'], 2)}", f"{fmt(s_['p10'], 2)} - {fmt(s_['p90'], 2)}", f"{fmt(s_['min'], 2)} - {fmt(s_['max'], 2)}"])
    L += table(["sample", "MSAs", "mean", "population-weighted mean", "median", "p25 - p75", "p10 - p90", "min - max"], rows)
    L.append("Selected MSAs (Table VI): " + ", ".join(f"{k} {fmt(v, 2)}" for k, v in sz["named"].items()) + ". Land-constrained coastal "
             "metros are below 1; large interior metros are above 2. These are metro-wide averages of housing units, with no "
             "separation of the city from the suburbs.\n")
    L.append("### Baum-Snow and Han (2024)\n")
    L.append("N. Baum-Snow and L. Han, 'The Microgeography of Housing Supply', *Journal of Political Economy* 132(6): 1897-1946. "
             "Tract-level supply elasticities for 2000-2010 changes (a ten-year horizon; the authors attribute most of the gap to Saiz's "
             "much larger 1970-2000 elasticities to the later, shorter period, introduction p. 4), identified from Bartik-type labor demand shocks to commuting destinations; the tract file below is "
             "the authors' predicted-elasticity release of September 2023 (50,410 tracts in 306 metro regions; "
             f"{bsh['n_rows']:,} tract-region rows, `wgt` = 1 / number of regions of a tract). Text (accepted "
             "manuscript, September 5, 2023): in the estimation sample the average elasticities are 0.35 (housing units), 0.42 (floorspace) and 0.09 "
             "(land development); predicted for all 50,410 urban tracts, unit and floorspace elasticities average 0.29 and 0.51 (quartiles "
             "0.14-0.44 and 0.33-0.70), with more variation within than between metro regions (introduction, p. 2); at the CBD the implied "
             "unit elasticity is 0.11, rising to 0.59 near halfway to the region edge (Table 5 column 1; section 5.2, p. 25); supply "
             "elasticities increase with CBD distance, flatten in the suburbs (about 30 to 90 percent of the way out) and rise again at the "
             "fringe because the developed fraction falls with distance while regulation is tighter in the suburbs than in central cities "
             "(Figure 3 discussion, p. 30). The file reproduces the text (FMM-IV linear specification, tract-weighted): "
             f"units {fmt(bsh['overall']['FMM_lin_units']['mean_tract'], 3)} (quartiles {fmt(bsh['overall']['FMM_lin_units']['p25_tract'], 2)}-"
             f"{fmt(bsh['overall']['FMM_lin_units']['p75_tract'], 2)}), floorspace {fmt(bsh['overall']['FMM_lin_space']['mean_tract'], 3)} "
             f"({fmt(bsh['overall']['FMM_lin_space']['p25_tract'], 2)}-{fmt(bsh['overall']['FMM_lin_space']['p75_tract'], 2)}). Aggregated to the "
             "metro-region level the FMM-IV unit elasticities average "
             f"{fmt(c['appendix_d']['bsh_metro_region_means_units'][0], 2)}-{fmt(c['appendix_d']['bsh_metro_region_means_units'][1], 2)} "
             f"(s.d. {fmt(c['appendix_d']['bsh_metro_region_sd'], 2)}; main text p. 33 and Online Appendix D), against 2.54 for Saiz's 1970-2000 estimates for the same "
             "metros; the authors attribute the gap to the later and shorter period, to aggregation (their aggregation weights tracts with large initial stocks, "
             "which respond less, more heavily) and to differences in the price and quantity measures.\n")
    rows = []
    for spec in ("FMM_lin", "FMM_quad", "IV_lin", "IV_quad"):
        rows.append([spec.replace("_", " "), *[fmt(bsh["overall"][f"{spec}_{t}"]["mean_tract"], 3) for t in ("units", "space", "devl")],
                     *[fmt(bsh["overall"][f"{spec}_{t}"]["mean_housing_units"], 3) for t in ("units", "space")]])
    L += table(["estimator (2001 developed fraction)", "units", "floorspace", "land development", "units, housing-unit weighted", "floorspace, housing-unit weighted"], rows)
    L.append("**By distance to the CBD** (fraction of the way from the CBD to the region edge; housing-unit weighted):\n")
    rows = [[r["bin"], pct(r["share_housing_units"]), fmt(r["mean_frac_developed"], 2), fmt(r["FMM_lin_units"], 2), fmt(r["FMM_lin_space"], 2),
             fmt(r["FMM_quad_units"], 2), fmt(r["FMM_quad_space"], 2), fmt(r["IV_lin_units"], 2), fmt(r["IV_lin_space"], 2)] for r in bsh["by_cbd_distance"]]
    L += table(["distance", "housing units %", "developed fraction", "FMM lin units", "FMM lin space", "FMM quad units", "FMM quad space",
                "IV lin units", "IV lin space"], rows)
    L.append("**By developed fraction of the tract in 2001:**\n")
    rows = [[r["bin"], pct(r["share_housing_units"]), fmt(r["mean_cbd_distance"], 2), fmt(r["FMM_lin_units"], 2), fmt(r["FMM_lin_space"], 2),
             fmt(r["FMM_quad_units"], 2), fmt(r["FMM_quad_space"], 2)] for r in bsh["by_developed_fraction"]]
    L += table(["developed fraction", "housing units %", "mean CBD distance", "FMM lin units", "FMM lin space", "FMM quad units", "FMM quad space"], rows)
    s_used = bsh["city_share_used"]
    sp_ = bsh["city_suburb"][f"city_share_{s_used}"]
    L.append(f"**Central city versus suburb.** Within each region the tracts closest to the CBD holding {pct(s_used)}% of occupied units "
             "(the city share of partitioned pupils in the locale partition) are the 'city' and the rest the 'suburb'; means are "
             "housing-unit weighted pooled over regions. The classification is by distance, not by the principal-city boundary, so it is "
             "a proxy for the locale split.\n")
    rows = []
    for spec in ("FMM_lin", "FMM_quad", "IV_lin", "IV_quad"):
        for t in ("units", "space"):
            v = sp_[f"{spec}_{t}"]
            w_ = sp_.get(f"{spec}_{t}_within_region_ratio")
            rows.append([spec.replace("_", " "), t, fmt(v["city"], 3), fmt(v["suburb"], 3), fmt(v["ratio"], 2),
                         "" if w_ is None else f"{fmt(w_['median'], 2)} / {fmt(w_['weighted_median'], 2)}"])
    L += table(["estimator", "quantity", "city (10-year)", "suburb (10-year)", "suburb / city", "within-region ratio, median / housing-weighted median"], rows)
    L.append("Mean developed fraction in 2001: city " + fmt(sp_["mean_frac_developed"]["city"], 2) + ", suburb "
             + fmt(sp_["mean_frac_developed"]["suburb"], 2) + ". Ratio of suburban to central-city elasticity (FMM linear, units) at other city shares: "
             + ", ".join(f"{k.replace('city_share_', '')}: {fmt(v['FMM_lin_units']['ratio'], 2)}" for k, v in bsh["city_suburb"].items()) + ".\n")
    L.append("### Recommended long-run eps_H\n")
    pc = e["period_check"]
    ad = c["appendix_d"]
    L.append("**Convention.** The model's housing market is a steady state, so the relevant elasticity is that of the housing stock after full adjustment. For a stock "
             "with annual depreciation delta, a response measured over T years is the fraction 1 - (1 - delta)^T of the long-run response: "
             "eps_LR = eps_T / (1 - (1 - delta)^T) (N. Baum-Snow and G. Duranton, 'Housing Supply and Housing Affordability', chapter for the Handbook of Regional and "
             "Urban Economics vol. 6, slides of the Urban Economics Association summer school 2025, slide 76 of 90). At delta = 3.5% the long run is "
             f"{fmt(e['longrun_factor_10y'], 2)} times a 10-year and {fmt(e['longrun_factor_30y'], 2)} times a 30-year elasticity "
             f"({fmt(c['longrun_factor']['delta_20'], 2)} and {fmt(c['longrun_factor_30y']['delta_20'], 2)} at 2%; {fmt(c['longrun_factor']['delta_15'], 2)} and "
             f"{fmt(c['longrun_factor_30y']['delta_15'], 2)} at 1.5%). **One convention is applied to every anchor**: Saiz (2010) is a 30-year estimate (1970-2000 changes), "
             "Baum-Snow and Han (2024) and the Saiz-method replication of their Online Appendix D are 10-year estimates (2000-10 changes). The slide adds that "
             "supply elasticities have fallen 'in US supply-constrained markets only' when Baum-Snow and Han are compared with Saiz. The slide and the "
             f"Appendix D statements quoted below were checked against the text of the cached files (slides: {res['c']['verification']['baum_snow_duranton_slides'].get('all_found')}; "
             f"manuscript: {res['c']['verification']['bsh_appendix'].get('all_found')}).\n")
    L.append("**Period check (Online Appendix D, pp. 17-18, and Table A8, p. 26, of the accepted manuscript).** Baum-Snow and Han replicate Saiz's estimation "
             "for 2000-2010 'using the same set of metros, specifications and variables': 'The average implied metro supply elasticity declines to "
             f"{fmt(ad['saiz_method_2000_2010_mean'], 2)} for 2000-2010, indicating that many metropolitan areas have become increasingly supply constrained over time', "
             f"against a mean of {fmt(ad['saiz_1970_2000_mean'], 2)} (s.d. {fmt(ad['saiz_1970_2000_sd'], 2)}) for Saiz's 1970-2000 estimates (Table A8 reports means of "
             f"{fmt(min(ad['table_a8_means_2000_2010']), 2)}-{fmt(max(ad['table_a8_means_2000_2010']), 2)} across its 2000-10 specifications and "
             f"{fmt(min(ad['table_a8_means_1970_2000']), 2)}-{fmt(max(ad['table_a8_means_1970_2000']), 2)} for the 1970-2000 ones). The ratio is {fmt(pc['observed_ratio'], 3)}. "
             f"The horizon alone predicts (1 - (1 - delta)^10) / (1 - (1 - delta)^30) = {fmt(pc['predicted_ratio'], 3)} at delta = 3.5%, i.e. a 2000-10 mean of "
             f"{fmt(pc['predicted_2000_10_mean'], 2)}, and a depreciation rate of {fmt(100 * pc['implied_delta'], 1)}% would match the ratio exactly. So the lower 2000-10 "
             "estimates are mostly the shorter horizon, and no additional period adjustment is made to the central anchor; the 2000-10 vintage converted with "
             f"the same factor gives {fmt(e['anchors']['saiz_pop95_period_refreshed']['long_run'], 2)} for the population-weighted 95-MSA mean (Saiz's mean x "
             f"{fmt(pc['observed_ratio'], 2)} x {fmt(e['longrun_factor_10y'], 2)}), close to the central anchor. Baum-Snow and Duranton's reading, that the fall is "
             "concentrated in supply-constrained markets, and the post-2010 zip-code evidence of Gorback and Keys (2024, 2010-20; slide 74: 'Find even smaller supply elasticities'), are "
             "reasons for a lower value in constrained coastal CZs, which the model's CZ-average locations do not distinguish; they are covered by the low anchor, and the "
             "unconverted 10-year values (0.30-0.76) are not long-run elasticities but are covered by the grid's 0.5.\n")
    rows = []
    for k, a_ in e["anchors"].items():
        pub = fmt(a_["published"], 2) if "published" in a_ else "-".join(fmt(x, 2) for x in a_["published_range"])
        if "published" in a_ and "published_range" in a_:
            pub += f" ({'-'.join(fmt(x, 2) for x in a_['published_range'])} across Table A8 specifications)"
        lr_ = fmt(a_["long_run"], 2) if "long_run" in a_ else "-".join(fmt(x, 2) for x in a_["long_run_range"])
        rows.append([a_["label"], a_["years"], pub, lr_, a_["role"]])
    L += table(["metro anchor (housing units)", "horizon (years)", "as published", f"long run (delta = 3.5%)", "role"], rows)
    L.append("Construction: eps_H(city) and eps_H(suburb) follow from the metro anchor and the within-metro ratio suburb / city of Baum-Snow and Han's predicted tract "
             f"elasticities, pooled housing-unit-weighted means on the city / suburb split matched to the pupil share of the locale partition (city share s = {fmt(e['city_share'], 3)}); "
             "with equal price changes in the two places the metro elasticity is s eps_c + (1 - s) eps_s (perfect substitutes, slide 77), so eps_c = anchor / (s + (1 - s) ratio) "
             f"and eps_s = ratio x eps_c. The central ratio {fmt(e['ratio_suburb_city'], 2)} is the mean of the four estimators for units; the range "
             f"{fmt(e['ratio_range'][0], 2)}-{fmt(e['ratio_range'][1], 2)} covers all {e['ratio_variants_n']} variants: units and floorspace, the four estimators, city shares "
             "0.25-0.5, and the within-region medians (matched split only: "
             f"{fmt(e['ratio_range_matched_split'][0], 2)}-{fmt(e['ratio_range_matched_split'][1], 2)}; the extremes are the {e['ratio_min_variant']['kind']} of "
             f"{e['ratio_min_variant']['spec'].replace('_', ' ')} {e['ratio_min_variant']['quantity']} at s = {e['ratio_min_variant']['city_share']} "
             f"and the {e['ratio_max_variant']['kind']} of {e['ratio_max_variant']['spec'].replace('_', ' ')} {e['ratio_max_variant']['quantity']} at s = {e['ratio_max_variant']['city_share']}). "
             "Mean ratio for units by city share: " + ", ".join(f"{k.replace('city_share_', '')}: {fmt(v, 2)}" for k, v in e["ratio_units_mean_by_city_share"].items()) + ".\n")
    L += table(["", "central", "low - high", "ten-year Baum-Snow-Han value (units; floorspace)"],
               [["eps_H city", fmt(e["eps_city"]["value"], 2), f"{fmt(e['eps_city']['lo'], 2)} - {fmt(e['eps_city']['hi'], 2)}",
                 f"{fmt(e['short_run_10yr_units']['city'], 2)}; {fmt(e['short_run_10yr_space']['city'], 2)}"],
                ["eps_H suburb", fmt(e["eps_suburb"]["value"], 2), f"{fmt(e['eps_suburb']['lo'], 2)} - {fmt(e['eps_suburb']['hi'], 2)}",
                 f"{fmt(e['short_run_10yr_units']['suburb'], 2)}; {fmt(e['short_run_10yr_space']['suburb'], 2)}"],
                ["metro average", fmt(e["eps_metro"]["value"], 2), f"{fmt(e['eps_metro']['lo'], 2)} - {fmt(e['eps_metro']['hi'], 2)}", ""],
                ["memo: no horizon conversion (Saiz's 30-year estimate read as long run: the first version of this report)",
                 f"{fmt(e['no_horizon_conversion']['city'], 2)} city / {fmt(e['no_horizon_conversion']['suburb'], 2)} suburb", "", ""]])
    L.append("Low - high is the extreme over the three metro anchors (low: Baum-Snow and Han tract elasticities, housing-unit weighted, x 3.34; high: Saiz's unweighted "
             "mean over 269 MSAs x 1.52) and the lowest and highest suburb / city ratio. **Why the central anchor is Saiz's population-weighted mean:** the model's CZs are "
             "pupil-weighted and mostly large metros, which the population weights approximate; it is the long-run object under the one convention; and the "
             "period check shows that the 2000-10 literature does not call for a lower long-run value once horizons are aligned. Baum-Snow and Han's tract "
             "elasticities, a within-metro identification closer to the model's city-versus-suburb margin, give a much lower value (low anchor; "
             f"{fmt(e['anchor_low'], 2)} for the metro, {fmt(e['eps_city']['lo'], 2)} for the city at the extremes), and they are a floor for sensitivity runs rather than "
             "the central case because the tract estimates are within-region elasticities identified from relative changes and aggregate to metro values "
             f"({fmt(ad['bsh_metro_region_means_units'][0], 2)}-{fmt(ad['bsh_metro_region_means_units'][1], 2)} x 3.34 = {fmt(e['anchors']['bsh_metro_regions']['long_run_range'][0], 2)}-"
             f"{fmt(e['anchors']['bsh_metro_regions']['long_run_range'][1], 2)}) well below Saiz's. Sensitivity to the depreciation rate (central anchor with each delta, the "
             "same ratio):\n")
    rows = [[fmt(100 * r_["delta"], 1) + "%", fmt(r_["factor_30y"], 2), fmt(r_["factor_10y"], 2), fmt(r_["metro"], 2), fmt(r_["city"], 2), fmt(r_["suburb"], 2),
             fmt(r_["bsh_tracts_metro"], 2)] for r_ in e["delta_sensitivity"]]
    L += table(["depreciation delta", "factor, 30 years", "factor, 10 years", "central metro anchor", "eps_H city", "eps_H suburb", "Baum-Snow-Han tract anchor"], rows)
    L.append("Suggested sensitivity grid for the model: eps_H in {0.5, 1, 2, 4} in both places and the elastic limit (which nests the fixed-cost model). "
             f"Floorspace elasticities are about {fmt(c['baum_snow_han']['overall']['FMM_lin_space']['mean_tract'] / c['baum_snow_han']['overall']['FMM_lin_units']['mean_tract'], 1)} "
             "times the unit ones in the Baum-Snow-Han data; if the model's H is read as floorspace or quality-adjusted housing services "
             "the values scale up by that factor. The recommended values are constructed from the literature (a metro anchor split by a within-metro ratio), "
             "not estimated for the model's two locations.\n")
    return L


def md_d(res: dict[str, Any]) -> list[str]:
    d = res["d"]
    st = d["stats"]
    cr = d["car_requirement"]
    L = ["## (d) Transport: the suburban car requirement\n"]
    L.append("ACS 5-year 2014-18 school-district tables, same panel, locale partition and draws as (a): B25044 tenure by vehicles "
             "available (occupied units), B25046 aggregate vehicles available (vehicles per household = B25046_001 / B25044_001), "
             "B08201 household size by vehicles available (a second zero-vehicle measure, one-person households and the household-size classes used below). Pooled over "
             "the districts of a location within each CZ, suburb minus city, CZ-weighted mean, whole-CZ bootstrap.\n")
    rows = []
    for k, v in st.items():
        rows.append([v["label"], fmt(v["level_1"], 3), fmt(v["level_2"], 3), fse(v["value"], v["boot_se"], 4),
                     f"[{fmt(v['ci95_lo'], 3)}, {fmt(v['ci95_hi'], 3)}]", v["n"]])
    L += table(["statistic", "city", "suburb", "gap (SE)", "95% interval", "CZs"], rows)
    L.append("### Within-CZ district regressions\n")
    L.append("District outcome on a suburb indicator, CZ fixed effects (weighted within-transformation), weights = occupied units "
             "(households) or pupils, log median household income of the district (`median_hh_income`, ACS B19013) as the income "
             "control, and (a further check) the district's renter share, which is a choice correlated with location, or the district's shares of "
             "1-, 3- and 4-or-more-person households (the 2-person share is the omitted category). SE (cluster) = "
             "CZ-clustered; SE (boot) = whole-CZ bootstrap with the same draws as everywhere else. 'log income quadratic' adds the squared "
             "deviation of log income from its mean.\n")
    rows = []
    names = {"zero_share": "zero-vehicle share", "veh_per_hh": "vehicles per household"}
    xnames = {"ln_inc": "log income", "ln_inc_sq": "log income squared (centered)", "renter_share": "renter share",
              "sh1": "share of 1-person households", "sh3": "share of 3-person households", "sh4": "share of 4+-person households"}
    for key, r in d["regressions"].items():
        b_ = r["beta"]
        xs = r["regressors"]
        other = "; ".join(f"{xnames[x]} {fmt(b_[i], 4)} ({fmt(r['se_cluster'][i], 4)}; {fmt(r['se_boot'][i], 4)})" for i, x in enumerate(xs) if i > 0)
        rows.append([names[r["y"]], r["weights"], r["spec"].replace("_", " "),
                     f"{fmt(b_[0], 4)} ({fmt(r['se_cluster'][0], 4)}; {fmt(r['se_boot'][0], 4)})", other,
                     r["n"], r["n_cz"], fmt(r["r2_within"], 3)])
    L += table(["outcome", "weights", "controls", "suburb (SE cluster; SE boot)", "other coefficients (SE cluster; SE boot)", "districts", "CZs", "within R2"], rows)
    rv = d["regressions"]
    L.append(f"Income explains {pct(1 - rv['veh_per_hh|households|log_income']['beta'][0] / rv['veh_per_hh|households|no_controls']['beta'][0], 0)}% "
             f"of the within-CZ suburb-city difference in vehicles per household and "
             f"{pct(1 - rv['zero_share|households|log_income']['beta'][0] / rv['zero_share|households|no_controls']['beta'][0], 0)}% of "
             "the difference in the zero-vehicle share; adding the renter share removes most of the rest, but tenure is itself a "
             "location-related choice (within tenure, the pooled gaps are +0.25 vehicles for owners and +0.24 for renters). "
             f"Adding the district household-size shares cuts the suburb coefficient on vehicles per household from "
             f"{fmt(rv['veh_per_hh|households|log_income']['beta'][0], 3)} to {fmt(rv['veh_per_hh|households|log_income_size_shares']['beta'][0], 3)}: "
             "most of what remains after income is household size.\n")
    hs = d["household_size"]
    L.append("### Household size: the vehicles gap holding household size fixed\n")
    L.append("Suburban households are larger and larger households own more vehicles, while the model's household is one adult, so the vehicles gap per household "
             "is separated from household composition with B08201 (household size by vehicles available), the school-district table of 1-, 2-, 3- and "
             "4-or-more-person households by vehicles (0, 1, 2, 3, 4+). Vehicles per household within a size class use the midpoints 0, 1, 2, 3 and, for '4 or more', "
             f"a mean calibrated by location to the aggregate B25046 / B25044 ({fmt(hs['four_plus_vehicles_city'], 2)} in the city, "
             f"{fmt(hs['four_plus_vehicles_suburb'], 2)} in the suburb); the class-based total reproduces the B25046-based one (largest difference over CZ x location cells "
             f"{fmt(hs['check_vs_b25046']['max_abs_diff'], 3)}; pooled raw gap {fmt(hs['check_vs_b25046']['raw_gap_class_based'], 3)} against "
             f"{fmt(hs['check_vs_b25046']['raw_gap_b25046'], 3)}). The 4-or-more-person class is open-ended, so larger suburban families inside it are not removed. Pooled "
             "over the districts of a location within each CZ, CZ-weighted, whole-CZ bootstrap:\n")
    rows = []
    for k in SIZE_KS:
        lab = SIZE_LABELS[k]
        cmp_, bc = hs["composition"][lab], hs["by_class"][lab]
        rv1 = hs["regressions"][f"vehicles|{lab}|log_income"]
        rz1 = hs["regressions"][f"zero_share|{lab}|log_income"]
        rows.append([lab, pct(cmp_["level_1"]), pct(cmp_["level_2"]), fse(100 * cmp_["value"], 100 * cmp_["boot_se"], 1),
                     fmt(bc["vehicles"]["level_1"], 2), fmt(bc["vehicles"]["level_2"], 2), fse(bc["vehicles"]["value"], bc["vehicles"]["boot_se"]),
                     fse(rv1["beta"][0], rv1["se_boot"][0]), fse(100 * bc["zero_share"]["value"], 100 * bc["zero_share"]["boot_se"], 1),
                     fse(100 * rz1["beta"][0], 100 * rz1["se_boot"][0], 1)])
    L += table(["household size", "share of households, city %", "suburb %", "gap pp (SE)", "vehicles per household, city", "suburb",
                "gap (SE)", "gap, income-adjusted (SE)", "zero-vehicle share gap pp (SE)", "zero-vehicle gap, income-adjusted pp (SE)"], rows)
    rows = [
        ["raw gap, all households (B25046 / B25044)", fse(d["stats"]["vehicles_per_hh"]["value"], d["stats"]["vehicles_per_hh"]["boot_se"], 3), ""],
        ["size-standardized gap, CZ-pooled size distribution (within-class gaps averaged with each CZ's pooled household-size shares)",
         fse(hs["standardized_gap_pooled"]["value"], hs["standardized_gap_pooled"]["boot_se"], 3),
         fse(100 * hs["standardized_zero_share_gap_pooled"]["value"], 100 * hs["standardized_zero_share_gap_pooled"]["boot_se"], 1)],
        ["  with the city's own household-size shares", fse(hs["standardized_gap_city_weights"]["value"], hs["standardized_gap_city_weights"]["boot_se"], 3), ""],
        ["  with the suburb's own household-size shares", fse(hs["standardized_gap_suburb_weights"]["value"], hs["standardized_gap_suburb_weights"]["boot_se"], 3), ""],
        ["size-standardized regression gap, no income control (district regressions within class, national class shares " + " / ".join(pct(x, 0) for x in hs["national_class_shares"]) + "%)",
         fse(hs["standardized_regression_vehicles"]["raw"]["value"], hs["standardized_regression_vehicles"]["raw"]["boot_se"], 3),
         fse(100 * hs["standardized_regression_zero_share"]["raw"]["value"], 100 * hs["standardized_regression_zero_share"]["raw"]["boot_se"], 1)],
        ["**size-standardized and income-adjusted regression gap** (log median household income control)",
         fse(hs["standardized_regression_vehicles"]["income_adjusted"]["value"], hs["standardized_regression_vehicles"]["income_adjusted"]["boot_se"], 3),
         fse(100 * hs["standardized_regression_zero_share"]["income_adjusted"]["value"], 100 * hs["standardized_regression_zero_share"]["income_adjusted"]["boot_se"], 1)],
    ]
    L += table(["standardized suburb-minus-city gap", "vehicles per household (SE)", "zero-vehicle share, pp (SE)"], rows)
    L.append("### Cost of a vehicle and the fixed goods cost\n")
    a_ = AAA_2018
    L.append(f"AAA, *Your Driving Costs*, 2018 edition (composite average vehicle, nine categories weighted by sales; new vehicle; "
             f"brochure p. 9; every number below was found in the cached brochure text: {res['d']['aaa_verification'].get('all_found')}): ownership ('fixed') costs "
             f"${fmt(a_['ownership_per_year'], 0)} a year = insurance ${fmt(a_['insurance'], 0)} + license, registration and taxes "
             f"${fmt(a_['license_registration_taxes'], 0)} + depreciation ${fmt(a_['depreciation_15k'], 0)} + finance ${fmt(a_['finance'], 0)}; "
             f"operating costs {fmt(a_['operating_cents_per_mile'], 2)} cents a mile (fuel {fmt(a_['fuel_cents_per_mile'], 2)}, maintenance, repair and "
             f"tires {fmt(a_['maintenance_cents_per_mile'], 2)}); total per year ${fmt(a_['total_per_year_10k'], 0)} at 10,000 miles, "
             f"**${fmt(a_['total_per_year_15k'], 0)} at 15,000** ({fmt(a_['cents_per_mile_15k'], 2)} cents a mile), ${fmt(a_['total_per_year_20k'], 0)} at 20,000. "
             "These are costs of a new vehicle held five years, so they overstate what the typical household pays.\n")
    cq = d["ce_per_quintile_vehicle_spending"]
    rd = d["ce_vehicles_rounding"]
    L.append("CE 2018 (quintiles of income before taxes) transportation lines and vehicles per consumer unit; private-vehicle spending = "
             "vehicle purchases (net outlay) + gasoline, other fuels and motor oil + other vehicle expenses (finance charges, maintenance and "
             "repairs, insurance, rentals, leases, licenses):\n")
    L += table(["group", "vehicles per CU", "with a vehicle %", "transportation $", "transportation % of expenditures", "private-vehicle spending $", "per vehicle $"],
               [[res["b"]["ce"]["quintiles"]["columns"][i], fmt(cq["vehicles"][i], 1), fmt(cq["at_least_one_vehicle_pct"][i + 1], 0), fmt(cq["transportation"][i], 0),
                 pct(cq["share_transportation"][i]), fmt(cq["private_vehicle_spending"][i], 0), fmt(cq["per_vehicle"][i], 0)] for i in range(5)])
    L.append(f"Across the five groups the elasticity of the vehicle count with respect to total expenditures is "
             f"{fmt(d['ce_expenditure_elasticity_quintiles']['vehicles'], 2)} (deciles {fmt(d['ce_expenditure_elasticity_deciles']['vehicles'], 2)}) and that of private-vehicle spending "
             f"{fmt(d['ce_expenditure_elasticity_quintiles']['private_vehicle_spending'], 2)}: spending on vehicles rises one for one with "
             "total expenditures, so a car is only partly a fixed cost. Cost per vehicle in the CE (all consumer units): private-vehicle "
             f"spending ${fmt(d['cost_per_vehicle_year']['ce_private_vehicle_spending_per_vehicle'], 0)}; ownership-like items (purchases, finance, insurance, "
             f"leases and licenses) ${fmt(d['cost_per_vehicle_year']['ce_ownership_like_per_vehicle'], 0)}. **The CE publishes vehicles per consumer unit to one decimal** "
             f"({fmt(rd['published_vehicles_per_cu'], 1)} overall; 1.0-2.7 by quintile), so every CE cost per vehicle carries a rounding error of up to "
             f"+-{pct(rd['relative_bound'])}% (all-CU private-vehicle spending ${fmt(rd['per_vehicle_private_lo'], 0)}-${fmt(rd['per_vehicle_private_hi'], 0)}), and the "
             "vehicle-count elasticity above is computed from the rounded counts.\n")
    cr = d["car_requirement"]
    one_exp = cr["single_adult"]["denominator"]
    L.append("**Conversion.** r_car = (extra vehicles per household) x (annual cost of a vehicle) / (mean expenditures per consumer unit, "
             f"${CE_TOTAL_EXP:,.2f}), a fixed goods cost of living in the suburb in units of mean consumption. The ratio of a per-household cost to per-household "
             "consumption is unit-free if both scale with household size (`housing_costs.md`), so the extra vehicles are standardized for household size; "
             f"for one-person households the denominator is the mean expenditure of one-person consumer units (${fmt(one_exp, 0)}, CE Table 1400). Extra vehicles per household, "
             "from the unadjusted to the composition-adjusted: raw pooled "
             f"{fse(cr['extra_vehicles']['vehicles_raw_pooled']['value'], cr['extra_vehicles']['vehicles_raw_pooled']['se'])}; income-adjusted regression, not adjusted for "
             f"household size (the first version) {fse(cr['extra_vehicles']['vehicles_income_adjusted']['value'], cr['extra_vehicles']['vehicles_income_adjusted']['se'])}; "
             f"size-standardized {fse(cr['extra_vehicles']['vehicles_size_standardized']['value'], cr['extra_vehicles']['vehicles_size_standardized']['se'])}; "
             "**size-standardized and income-adjusted "
             f"{fse(cr['extra_vehicles']['vehicles_size_standardized_income_adjusted']['value'], cr['extra_vehicles']['vehicles_size_standardized_income_adjusted']['se'])} (central)**; "
             f"one-person households {fse(cr['extra_vehicles']['vehicles_one_person_raw']['value'], cr['extra_vehicles']['vehicles_one_person_raw']['se'])} raw and "
             f"{fse(cr['extra_vehicles']['vehicles_one_person_income_adjusted']['value'], cr['extra_vehicles']['vehicles_one_person_income_adjusted']['se'])} income-adjusted; "
             "income and household-size shares as controls "
             f"{fse(cr['extra_vehicles']['vehicles_size_controls_regression']['value'], cr['extra_vehicles']['vehicles_size_controls_regression']['se'])} (an over-control, "
             "since the size mix is itself a location outcome: the low end). On the extensive margin the fall in the zero-vehicle share is "
             f"{fse(cr['extra_vehicles']['no_vehicle_share_raw_pooled']['value'], cr['extra_vehicles']['no_vehicle_share_raw_pooled']['se'])} raw, "
             f"{fse(cr['extra_vehicles']['no_vehicle_share_income_adjusted']['value'], cr['extra_vehicles']['no_vehicle_share_income_adjusted']['se'])} income-adjusted and "
             f"{fse(cr['extra_vehicles']['no_vehicle_share_size_standardized_income_adjusted']['value'], cr['extra_vehicles']['no_vehicle_share_size_standardized_income_adjusted']['se'])} "
             "size-standardized and income-adjusted, times the cost of one vehicle.\n")
    names_c = {"ce_private_vehicle_spending_per_vehicle": "CE private-vehicle spending per vehicle",
               "ce_ownership_like_per_vehicle": "CE ownership-like spending per vehicle",
               "aaa_ownership_per_vehicle_year": "AAA ownership (fixed) costs",
               "aaa_total_per_vehicle_year_15k_miles": "AAA total, 15,000 miles",
               "aaa_total_per_vehicle_year_10k_miles": "AAA total, 10,000 miles"}
    shown = ("ce_private_vehicle_spending_per_vehicle", "ce_ownership_like_per_vehicle", "aaa_total_per_vehicle_year_15k_miles")
    rows = [[r["extra_label"], fmt(r["extra_value"], 3), names_c[r["cost"]], fmt(r["cost_usd"], 0), fmt(r["annual_usd"], 0),
             fmt(r["denominator"], 0), fse(r["r_car"], r["r_car_se"], 4)] for r in cr["table"] if r["cost"] in shown]
    L += table(["extra", "value", "cost per vehicle-year", "$", "annual $ per household", "denominator $", "r_car (SE from the extra-vehicle SE only)"], rows)
    L.append("(The JSON has all five cost measures, including AAA ownership costs and the AAA total at 10,000 miles.)\n")
    L.append(f"**Central: r_car = {fmt(cr['central']['r_car'], 4)} (SE {fmt(cr['central']['se'], 4)})** = size-standardized, income-adjusted extra vehicles x CE private-vehicle spending per "
             f"vehicle (${fmt(cr['central']['r_car'] * CE_TOTAL_EXP, 0)} a year), against {fmt(cr['first_version']['r_car'], 4)} in the first version, which adjusted for income but not for "
             f"household size (a {fmt(-100 * (cr['central']['r_car'] / cr['first_version']['r_car'] - 1), 0)}% reduction). Range {fmt(cr['low']['r_car'], 4)} (income and household-size shares as controls, CE private-vehicle "
             f"cost) to {fmt(cr['high']['r_car'], 4)} (size-standardized, no income adjustment, AAA total cost at 15,000 miles; the unadjusted raw gap with the AAA cost gives "
             f"{fmt(cr['unadjusted_upper']['r_car'], 4)}, the first version's upper end). **Single-adult reading** (one-person households, income-adjusted, divided by one-person "
             f"consumer units' expenditures): {fse(cr['single_adult']['r_car'], cr['single_adult']['se'], 4)}; it is higher than the central value because vehicles scale less "
             "than proportionally with household expenditure (economies of scale in car use). Extensive margin only (size-standardized, income-adjusted): "
             f"{fmt(cr['extensive_margin_central'], 4)}. For scale, the housing-cost gap of `housing_costs.md` is r_p = 0.0276 (0.0044), so the car requirement adds "
             f"about {fmt(100 * cr['central']['r_car'] / 0.0276, 0)}% to the residence charge if the whole housing gap were also a fixed cost (the same household-size caveat applies "
             "to the housing-cost gap: suburban households are larger and buy more housing).\n")
    return L


def md_e(res: dict[str, Any]) -> list[str]:
    e = res["e"]
    b = e["bfm"]
    pm = e["pums90"]
    ei = e["elasticity_income"]
    L = ["## (e) Income gradient of willingness to pay (BFM 2007)\n"]
    L.append(f"{b['source']}. Mean MWTP ${b['mean_mwtp']} (SE {b['mean_mwtp_se']}) per month per school-average SD of test scores (monthly user cost of "
             f"housing, 1990 dollars); household income (+$10,000) {b['per_10k_income']} (SE {b['per_10k_income_se']}); college degree or more vs some "
             f"college or less {b['college']} (SE {b['college_se']}); children under 18 vs none {b['children']} (SE {b['children_se']}). Mean household income "
             f"${b['mean_income']:,.0f} (SD ${b['sd_income']:,.0f}; {b['income_source']}). The education variable is the household head's.\n")
    L.append(f"**Elasticity at the means.** (1.38 / 10,000) x 54,103 / 19.69 = **{fmt(ei['value'], 4)}**; SE {fmt(ei['se_conditional_on_means'], 4)} holding the means "
             f"fixed (0.33 / 10,000 x 54,103 / 19.69), {fmt(ei['se_with_mean_mwtp_independent'], 3)} if the SE of the mean MWTP (7.41) is added "
             "independently (Table 8 gives no covariance; with demeaned interaction characteristics the mean and the slopes are approximately "
             "uncorrelated). The elasticity of the WTP *share* of income is " + fmt(ei["implied_share_elasticity"], 3) + ": school quality is a necessity "
             "in BFM's data. In the linear Table 8 gradient, WTP as a share of monthly income is "
             + ", ".join(f"{fmt(100 * r['wtp_share_monthly_income'], 2)}% at {fmt(r['income_multiple'], 1)} times the mean income" for r in e["wtp_share_by_income"]) + ".\n")
    L.append("### Folding in the college-degree term\n")
    L.append("Table 8 holds the other characteristics fixed. If households with more income are also more likely to be college educated, the "
             "total gradient of WTP with respect to income is 1.38 per $10,000 plus 13.03 times the slope of the college probability on income. "
             "That slope is measured in the Census 1990 PUMS 5% sample for the six Bay Area counties (" + pm["source"].split(", households in")[0].replace("Census 1990 PUMS 5% sample, ", "")
             + "; households in the Oakland, San Francisco and San Jose PMSAs; householder with a bachelor's degree or more; total household income in 1989).\n")
    L += table(["statistic", "PUMS 1990 Bay Area", "BFM sample"],
               [["households (records / weighted)", f"{pm['n_households']:,} / {pm['weighted_households']:,.0f}", "242,100"],
                ["mean household income $", fmt(pm["mean_income"], 0), fmt(b["mean_income"], 0)],
                ["SD of household income $", fmt(pm["sd_income"], 0), fmt(b["sd_income"], 0)],
                ["householders with a bachelor's degree or more", pct(pm["share_college_head"]) + "%", "44% (block-group average)"],
                ["mean income, college / non-college $", f"{fmt(pm['mean_income_college'], 0)} / {fmt(pm['mean_income_noncollege'], 0)}", ""],
                ["correlation of college with income", fmt(pm["corr_college_income"], 3), ""],
                ["slope of college on income, per $10,000 (linear projection)", f"{fmt(pm['gamma_per_10k'], 4)} (bootstrap SE {fmt(10000 * pm['gamma_se_boot'], 5)})", ""],
                ["local slope at $54,103: logit / local linear, per $10,000", f"{fmt(pm['slope_logit_at_bfm_mean_per_10k'], 4)} / {fmt(pm['slope_local_linear_per_10k'], 4)}", ""]])
    el = e["elasticity_with_college"]
    rows = [[lab, fmt(el[k]["slope_per_10k"] - el[k]["college_term_per_10k"] + el[k]["college_term_per_10k"], 3), fmt(el[k]["college_term_per_10k"], 3),
             fmt(el[k]["elasticity"], 3), fmt(el[k]["se"], 3), fmt(el[k]["college_term_elasticity"], 3)]
            for lab, k in (("linear projection (Cov / Var), PUMS", "linear_projection"), ("local logit at the BFM mean, PUMS", "local_logit"),
                           ("local linear, +-$20,000 window, PUMS", "local_linear"),
                           ("linear projection with BFM's moments (44% college, SD $50,719) and the PUMS correlation", "bfm_moments_pums_correlation"))]
    L += table(["slope of college on income", "total WTP gradient per $10,000", "of which college term", "elasticity", "SE", "college part of the elasticity"], rows)
    L.append(f"Assumptions: independent sampling errors of the two Table 8 coefficients "
             "(income and college are positively correlated, so their covariance is probably negative and the SE conservative); the children, "
             "race, work-status and age terms are held fixed; the PUMS Bay Area sample differs from BFM's (mean income $50,775 vs $54,103, 35.5% "
             "vs 44% college heads), so the last row rescales the slope to BFM's reported moments. **Approximate elasticity with the college term: "
             f"{fmt(el['linear_projection']['elasticity'], 3)} (SE {fmt(el['linear_projection']['se'], 3)}), {fmt(min(v['elasticity'] for v in el.values() if isinstance(v, dict)), 2)}"
             f"-{fmt(max(v['elasticity'] for v in el.values() if isinstance(v, dict)), 2)} across the four slope measures.**\n")
    va = e["elasticity_variants"]
    L.append("### Education definition and the children term\n")
    L.append("Table 8 does not define 'college degree or more' and BFM's Table 1 has 44% of the heads of household of a census block group with a college degree or more "
             f"(block-group average); the PUMS share of Bay Area householders with a bachelor's degree or more is {pct(pm['share_college_head'])}% and with an associate "
             f"degree or more (YEARSCH >= 12) {pct(pm['share_assoc_plus_head'])}%, so both are shown. Table 8 also has a children-under-18 coefficient "
             f"({b['children']}, SE {b['children_se']}); the 1990 PUMS share of households with a person under 18 is {pct(pm['share_child_present'])}% and its linear slope on "
             f"household income {fmt(pm['gamma_child_per_10k'], 4)} per $10,000 (bootstrap SE {fmt(10000 * pm['gamma_child_se_boot'], 5)}). Slopes are linear projections Cov / Var; "
             "SEs combine the independent sampling errors of the Table 8 coefficients and the bootstrap SEs of the slopes (300 resamples of the PUMS households).\n")
    rows = []
    for lab, k, share, slope, slope_se in (
            ("bachelor's degree or more", "bachelors", pm["share_college_head"], pm["gamma_per_10k"], pm["gamma_se_boot"]),
            ("associate degree or more", "associate_plus", pm["share_assoc_plus_head"], pm["gamma_assoc_plus_per_10k"], pm["gamma_assoc_plus_se_boot"]),
            ("bachelor's degree or more + children", "bachelors_children", pm["share_college_head"], pm["gamma_per_10k"], pm["gamma_se_boot"]),
            ("associate degree or more + children", "associate_plus_children", pm["share_assoc_plus_head"], pm["gamma_assoc_plus_per_10k"], pm["gamma_assoc_plus_se_boot"])):
        v = va[k]
        rows.append([lab, pct(share) + "%", f"{fmt(slope, 4)} ({fmt(10000 * slope_se, 5)})", fmt(v["slope_per_10k"], 3),
                     " + ".join(fmt(x, 3) for x in v["terms_elasticity"]), fmt(v["elasticity"], 3), fmt(v["se"], 3)])
    L += table(["college definition (and children term)", "PUMS share", "slope of college on income per $10,000 (SE)", "total WTP gradient per $10,000",
                "college (+ children) part of the elasticity", "elasticity", "SE"], rows)
    L.append(f"The associate-or-more definition reproduces BFM's 44% (and the slope, 0.0357 against 0.0358 per $10,000) and leaves the elasticity unchanged at "
             f"{fmt(va['associate_plus']['elasticity'], 3)}; adding the children term raises it by {fmt(va['bachelors_children']['elasticity'] - va['bachelors']['elasticity'], 3)} to "
             f"{fmt(va['bachelors_children']['elasticity'], 3)}. The range across the slope measures and these variants is "
             f"{fmt(res['recommended']['wtp_income_elasticity']['with_college_lo'], 3)}-{fmt(res['recommended']['wtp_income_elasticity']['with_college_hi'], 3)}.\n")
    return L


def md_caveats(res: dict[str, Any]) -> list[str]:
    L = ["## Caveats\n"]
    L += ["- **(a) Cost, not price.** ACS monthly housing costs are gross rent (rent plus utilities) and owners' selected monthly costs "
          "(mortgage payments including principal, taxes, insurance, utilities; no maintenance); a percentile gap mixes the price of housing "
          "services, quantity and quality, and, for owners, saving. Outright owners' cash costs (taxes, insurance, utilities) are not the price of the "
          "housing consumed, and they fill the bottom of the pooled distribution, so only the renter (gross rent) profile is recommended for the lower tail. "
          "Units with no cash rent are excluded. Percentiles come from binned counts "
          "(linear interpolation inside bins; the median reproduces the published district medians to a mean log error of 0.002, "
          "`housing_costs.md`). Elementary and secondary district layers overlap geographically with unified ones in a few places, as in "
          "`housing_costs.py`. The tenure split is within the same CZ partition, so the renter / owner gaps include the tenure-specific mix of "
          "structures, and none of the percentile gaps is adjusted for household size (suburban households are larger and buy more housing).",
          "- **(b) Units, household size and income.** The CE counts consumer units (2.5 persons, 1.3 earners), not adults, and expenditures, not income; "
          "its income quintiles mix permanent and transitory income (expenditures of the lowest quintile are 2.3 times its income). Household size "
          "rises with income (1.7 to 3.2 persons per CU), which the equivalence scales address by assumption: the scale is a convention (the OECD "
          "table of elasticities runs from 0 to 1), the central square root is supported by the fit of the curve (two degrees of freedom) and not "
          "estimated, people and children per CU are published to one decimal and applied to the group means (the within-group dispersion of household size is "
          "ignored), 'children' are under 18, and the CE Table 1400 cross-check suggests that "
          "housing economies of scale exceed the square root (one-person households spend more than the equivalized curve implies), which would raise the "
          "subsistence share of single adults above the central value. The ratio h / (mean expenditure per equivalent adult) is the model's hmin only "
          "under the assumption that housing and consumption scale together with household size (`housing_costs.md`). Total expenditures include personal "
          "insurance and pensions (11.9%), which the model treats as saving; the net-of-pensions fit is reported at the consumer-unit level only. CE shelter "
          "excludes mortgage principal; utilities include telephone service. The Stone-Geary shares ignore the price of housing (relative prices differ across "
          "areas), and two parameters (three with the scale) are fitted to five or ten group means, so the residual SEs have 3 (8) degrees of freedom. The SE (CE) "
          "treats item and total as independent draws (conservative). The tercile shares are consumer-unit income terciles, so they include household-size "
          "differences. ACS burdens are relative to household income with no size adjustment, include utilities and principal, and the mean burden needs class "
          "values for the open classes.",
          "- **(c) Elasticities.** The recommended values are constructed (a metro anchor split by a within-metro ratio), not estimated for the model's two locations. "
          "The horizon conversion assumes full stock adjustment with constant depreciation (3.5% a year, from the Baum-Snow and Duranton slide; 2% would give "
          "factors of 2.20 and 5.47) and equal price changes in city and suburb; it is applied to Saiz's 30-year estimates (1970-2000, metro-wide, housing units, no "
          "separation of city and suburb) and to Baum-Snow and Han's 2000-2010 tract predictions (identified from within-region variation in access to jobs) alike. The "
          "Online Appendix D numbers are read from the accepted manuscript of September 2023 (the text refers to 237 metros, the Table A8 notes to predicted values for "
          "the full Saiz sample of 269 metros; the 1.25 is the text's figure, Table A8 reports 0.80-1.30 across specifications) and may differ in the published JPE version. "
          "The city / suburb split of the tracts is by distance to the CBD at the pupil share of the city, not by principal-city boundaries. The model's H is a composite "
          "of quantity and quality: floorspace elasticities are higher than unit elasticities.",
          "- **(d) Vehicles.** B25044, B25046 and B08201 count vehicles available to occupied housing units (households, not adults). The household-size classes use "
          "midpoints for 0-3 vehicles and a calibrated mean for '4 or more'; the 4-or-more-person class is open-ended; the within-class regressions use district medians of "
          "household income, so they control for income imperfectly (within-district dispersion), and the standardized gap uses national class shares; the renter share absorbs "
          "most of the remaining difference but tenure is itself a location choice, and the regression with household-size shares as controls over-controls (the size mix is a "
          "location outcome). A car is partly a normal good (vehicle count elasticity about 0.65 in the CE, spending elasticity about 1), so the fixed-cost reading applies to "
          "the first vehicle; the CE vehicles per consumer unit are published to one decimal (+-2.6% on the cost per vehicle, and the elasticity is computed from rounded counts); "
          "the AAA figures are for new vehicles driven 10,000-20,000 miles. The conversion divides household costs by expenditures per consumer unit as `housing_costs.py` does; "
          "the one-person variant divides by one-person consumer units' expenditures.",
          "- **(e) BFM.** The WTP gradient is a coefficient of a heterogeneous-choice model with 135 interactions, conditional on the other "
          "household characteristics; Table 8 gives no covariances; mean income comes from the working paper version. The 1990 association "
          "uses the public PUMS 5% sample for the six counties, not BFM's restricted long-form data, and the householder's education; BFM do not define "
          "'college degree or more' (bachelor's and associate-or-more are both shown) and the children term is the household-level presence of a person under 18."]
    L.append("")
    return L


def md_sources() -> list[str]:
    L = ["## Sources\n"]
    L += ["- ACS 5-year 2014-18 (vintage 2018) school-district tables B25104 (via `acs_housing.py`), B25063, B25087, B25044, B25046, B08201 "
          "(Census API, `housing_transport.py fetch`); national B25074, B25095, B25106, B19001, B19025, B25003 (acs5 and acs1 2018).",
          "- BLS Consumer Expenditure Survey 2018: Table 1101 (quintiles of income before taxes), Table 1110 (deciles) and Table 1400 (size of consumer unit) at "
          "https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error.htm; Table R-1 mean expenditures "
          "($61,224.13) via `housing_costs.py`.",
          "- M. F. Foerster and M. Mira d'Ercole (OECD), 'The OECD approach to measuring income distribution and poverty: strengths, limits and statistical issues', Table 2 "
          "(equivalence scales and elasticities: per capita 1, 'old' OECD 0.73, OECD-modified 0.53, square root 0.50, household 0), "
          "https://cipe.umd.edu/conferences/oecdumd/conf_papers/Papers/The_OECD_Approach_to_Measuring_Income_Distribution_and_Poverty.pdf; B. Buhmann, L. Rainwater, "
          "G. Schmaus and T. Smeeding (1988), Review of Income and Wealth 34(2): 115-142 (sensitivity of inequality and poverty to the equivalence elasticity, LIS).",
          "- AAA, Your Driving Costs, 2018 edition: https://exchange.aaa.com/wp-content/uploads/2018/09/18-0090_2018-Your-Driving-Costs-Brochure_FNL-Lo-5-2.pdf",
          "- A. Saiz (2010), QJE 125(3): 1253-1296, https://doi.org/10.1162/qjec.2010.125.3.1253; MSA elasticities file via https://urbaneconomics.mit.edu/research/data.",
          "- N. Baum-Snow and L. Han (2024), JPE 132(6): 1897-1946, https://doi.org/10.1086/728110; accepted version (September 5, 2023; introduction p. 4, Figure 4 "
          "discussion p. 33, Online Appendix D pp. 17-18 and Table A8 p. 26) and tract elasticity estimates (gammas_hat_all.dta, 'Explanation of Housing Supply Elasticity "
          "Estimates', September 2023) from the authors' page https://sites.google.com/site/baumsnow/research; replication package https://doi.org/10.7910/DVN/LFGAAW (not downloaded).",
          "- N. Baum-Snow and G. Duranton, 'Housing Supply and Housing Affordability', chapter for the Handbook of Regional and Urban Economics vol. 6; slides of the Urban Economics "
          "Association summer school, July 2025, https://urbaneconomics.org/workshops/summerschool2025/slides/baumsnow.pdf (slide 73 Saiz, slide 74 Baum-Snow and Han, slide 76 "
          "long-run conversion, slide 77 aggregation).",
          "- P. Bayer, F. Ferreira and R. McMillan (2007), JPE 115(4): 588-638, Table 8 (p. 625), https://doi.org/10.1086/522381; Table 1 (44% of block-group heads with a "
          "college degree or more); mean household income from NBER w10871 Table 1 (`literature_targets.md`).",
          "- Census 1990 PUMS 5% sample (PUMS A), California: https://www2.census.gov/census_1990/1990_PUMS_A/PUMSAXCA.zip, data dictionary "
          "PUMSUSDD.doc in the same directory's DOCUMENT folder (YEARSCH 12-13 associate, 14 bachelor's, R18UNDR).",
          "- Internal: `housing_costs.md` (r_p, median gap), `district_settings.py` (draws), `literature_targets.md` (BFM), `calibration_log.md`, "
          "`runs/experiments-2026-10b/DESIGN.md` and `DESIGN_v2.md` (M1, M1b), `M1B_SPEC.md`."]
    L.append("")
    return L


# --------------------------------------------------------------------------
# TOML targets
# --------------------------------------------------------------------------

def tv(v: Any) -> str:
    """TOML value."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(tv(x) for x in v) + "]"
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "nan"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return repr(float(v))


def toml_section(L: list[str], name: str, items: dict[str, Any]) -> None:
    L.append(f"[{name}]")
    for k, v in items.items():
        if v is not None:
            L.append(f"{k} = {tv(v)}")
    L.append("")


def toml_text(res: dict[str, Any], meta: dict[str, Any]) -> str:
    rec = res["recommended"]
    L = ["# housing_transport_targets.toml -- evidence for income-driven demand and an endogenous housing market (DESIGN.md M1, M1b).",
         f"# Generated {meta['generated']} by data/spatial/housing_transport.py; see housing_transport.md for definitions, sources and caveats.",
         "# value / se: estimate and bootstrap (whole-CZ, B = %d, seed %d) or sampling SE; lo / hi: plausible range used for sensitivity runs." % (meta["B"], meta["seed"]),
         "# Log gaps are suburb minus city (locale partition, 188 CZs, CZ weights = partitioned pupils); h-type quantities are in units of mean",
         "# consumption (CE 2018 mean expenditures per consumer unit $61,224.13; for hmin the mean expenditure per EQUIVALENT ADULT, see housing_demand.hmin), the model's Cbar.", ""]
    toml_section(L, "meta", dict(generated=meta["generated"], script="data/spatial/housing_transport.py", year=YEAR, B=meta["B"], seed=meta["seed"],
                                 n_cz=meta["n_cz"], ce_mean_expenditures=CE_TOTAL_EXP))
    if "a" in res:
        a = res["a"]
        g50 = a["tables"]["all"]["percentiles"]["p50"]["log_gap"]
        toml_section(L, "gap_housing_cost", dict(value=g50["value"], se=g50["boot_se"], n_cz=g50["n"],
                                                 note="log gap of the pooled median monthly housing cost (B25104), the target for a_H; = housing_costs.md"))
        toml_section(L, "validation_lower_tail", dict(
            use="housing_cost_profile.renters and cheap_units.renters",
            do_not_use="housing_cost_profile.all_units and housing_cost_profile.owners below the median: the bottom of those distributions is mostly outright owners whose cash costs are not the price of the housing consumed",
            note="gross rent is the price of the unit; see tail_composition"))
        for tab, nm in (("all", "all_units"), ("rent", "renters"), ("own", "owners")):
            d = a["tables"][tab]
            for k in PCT_KEYS:
                v = d["percentiles"][k]
                items = dict(value=v["log_gap"]["value"], se=v["log_gap"]["boot_se"], n_cz=v["log_gap"]["n"], city_usd=v["level_1"], suburb_usd=v["level_2"],
                             usd_gap=v["usd_gap"]["value"], usd_gap_se=v["usd_gap"]["boot_se"])
                if k in d["diff_vs_median"]:
                    items.update(minus_median_gap=d["diff_vs_median"][k]["value"], minus_median_gap_se=d["diff_vs_median"][k]["boot_se"])
                toml_section(L, f"housing_cost_profile.{nm}.{k}", items)
        for tab, nm in (("all", "all_units"), ("rent", "renters"), ("own", "owners")):
            for q in ("q10", "q25", "q50"):
                v = a["tables"][tab]["cheap_share"][q]
                toml_section(L, f"cheap_units.{nm}.{q}", dict(city=v["level_1"], suburb=v["level_2"], gap=v["value"], se=v["boot_se"]))
        tc = a["tail_composition"]
        for k, lab in TENURE_GROUPS:
            v = tc["occupancy"][k]
            toml_section(L, f"tail_composition.occupancy.{k}", dict(city=v["level_1"], suburb=v["level_2"], gap=v["value"], se=v["boot_se"], label=lab))
        for q in ("q10", "q25", "q50"):
            for k, lab in TENURE_GROUPS:
                v = tc["below"][q][k]
                toml_section(L, f"tail_composition.below_{q}.{k}", dict(city=v["level_1"], suburb=v["level_2"], gap=v["value"], se=v["boot_se"],
                                                                          note="share of the units below the CZ-pooled all-units percentile that belong to the group"))
    if "alphaH" in rec:
        r, h, s_ = rec["alphaH"], rec["hmin"], rec["housing_share"]
        toml_section(L, "housing_demand.alphaH", dict(value=r["value"], se=r["se"], se_resid=r["se_resid"], lo=r["lo"], hi=r["hi"], decile_fit=r["decile_fit"],
                                                      cu_level=r["cu_level"], scale=CENTRAL_SCALE, note=r["basis"]))
        toml_section(L, "housing_demand.hmin", dict(value=h["value"], se=h["se"], se_resid=h["se_resid"], lo=h["lo"], hi=h["hi"], decile_fit=h["decile_fit"],
                                                    cu_level=h["cu_level"], scale=CENTRAL_SCALE, h_usd_per_equivalent_adult=h["h_usd"],
                                                    h_usd_mean_size_cu=h["h_usd_mean_size_cu"], mean_exp_per_equivalent_adult=h["mean_exp_per_equivalent_adult"],
                                                    note=h["basis"]))
        for sname, v in r["by_scale"].items():
            toml_section(L, f"housing_demand.by_scale.{sname}", dict(alphaH=v["alphaH"], hmin=v["hmin"], h_usd=v["h_usd"], rmse_pp=v["rmse_pp"],
                                                                    mean_exp_per_unit=v["mean_exp_per_unit"],
                                                                    note="quintiles, shelter + utilities; none = consumer-unit level (the first version's central values)"))
        toml_section(L, "housing_demand.housing_share", dict(
            aggregate=s_["value"], aggregate_se=s_["se"], mean_of_group_shares=s_["mean_of_group_shares"],
            mean_of_group_shares_deciles=s_["mean_of_group_shares_deciles"],
            by_income_tercile=s_["by_income_tercile"], by_income_tercile_se=s_["by_income_tercile_se"],
            by_income_tercile_mean_of_group_shares=s_["by_income_tercile_mean_of_group_shares"],
            by_income_tercile_from_quintiles=s_["by_income_tercile_quintiles"], by_income_tercile_from_quintiles_se=s_["by_income_tercile_quintiles_se"],
            one_person_cu=s_["one_person_cu"], two_plus_cu=s_["two_plus_cu"],
            model_moment="sum_i R H_i / sum_i Y'_i over adults (ratio of aggregates); by tercile of pre-tax income: the same ratio within the tercile",
            definition=s_["definition"], note=s_["basis"]))
        toml_section(L, "housing_demand.model_definitions", dict(note=rec["model_definitions"]))
        ce = res["b"]["ce"]
        for tag in ("quintiles", "deciles"):
            for nm in HOUSING_DEFS:
                f = ce[tag]["fits"][nm]
                toml_section(L, f"housing_demand.fits.{tag}.{nm}", dict(alphaH=f["alphaH"], alphaH_se=f["se_ce"]["alphaH"], hmin=f["h_share_of_mean_exp"],
                                                                        hmin_se=f["se_ce"]["h_share"], h_usd=f["h_usd"], share_at_mean=f["s_at_mean_exp"],
                                                                        elasticity_at_mean=f["elasticity_at_mean"], rmse_pp=f["rmse_pp"], r2=f["r2"],
                                                                        note="consumer-unit level, no household-size adjustment"))
        eqb = res["b"]["equivalized"]
        for tag in ("quintiles", "deciles"):
            for sname in eqb[tag]["fits"]:
                for nm in HOUSING_DEFS:
                    f = eqb[tag]["fits"][sname][nm]
                    toml_section(L, f"housing_demand.fits_equivalized.{sname}.{tag}.{nm}", dict(
                        alphaH=f["alphaH"], alphaH_se=f["se_ce"]["alphaH"], hmin=f["h_share_of_mean_exp"], hmin_se=f["se_ce"]["h_share"],
                        h_usd=f["h_usd"], mean_exp_per_equivalent_adult=f["mean_exp_used"], share_at_mean=f["s_at_mean_exp"],
                        elasticity_at_mean=f["elasticity_at_mean"], rmse_pp=f["rmse_pp"], r2=f["r2"]))
            tp = eqb[tag]["theta_profile"]
            toml_section(L, f"housing_demand.theta_profile.{tag}", dict(
                theta=[r_["theta"] for r_ in tp["rows"]], alphaH=[r_["alphaH"] for r_ in tp["rows"]], hmin=[r_["h_share"] for r_ in tp["rows"]],
                rmse_pp=[r_["rmse_pp"] for r_ in tp["rows"]], best_theta=tp["best"]["theta"], best_alphaH=tp["best"]["alphaH"], best_hmin=tp["best"]["h_share"],
                best_rmse_pp=tp["best"]["rmse_pp"], note="shelter + utilities, persons as the size measure; best on a 0.01 grid (two degrees of freedom: suggestive)"))
        for nm in HOUSING_DEFS:
            toml_section(L, f"housing_demand.shares.{nm}", dict(by_quintile=ce["quintiles"]["defs"][nm]["share"], se_by_quintile=ce["quintiles"]["defs"][nm]["share_se_rho0"],
                                                                aggregate=ce["quintiles"]["defs"][nm]["all_cu_share"], aggregate_se=ce["quintiles"]["defs"][nm]["all_cu_share_se_rho0"],
                                                                by_tercile=ce["quintiles"]["terciles"][nm]["ratio_of_means"], se_by_tercile=ce["quintiles"]["terciles"][nm]["se_rho0"],
                                                                by_tercile_from_deciles=ce["deciles"]["terciles"][nm]["ratio_of_means"],
                                                                se_by_tercile_from_deciles=ce["deciles"]["terciles"][nm]["se_rho0"]))
        for nm_, v in res["b"]["size_check"].items():
            toml_section(L, "housing_demand.size_check." + ("one_person" if nm_ == "One person" else "two_plus"), dict(
                people=v["people"], expenditures=v["expenditures"], share=v["share"], share_se=v["share_se"], predicted_cu_level=v["predicted_cu_level"],
                predicted_sqrt_scale=v["predicted_sqrt_scale"], excess_over_sqrt_scale_pp=v["excess_over_sqrt_scale_pp"],
                implied_subsistence_share_given_aH=v["implied_subsistence_share_given_aH_sqrt"], implied_subsistence_usd=v["implied_subsistence_usd"],
                note="CE Table 1400; one-point cross-check, not a fit"))
        a1 = res["b"]["acs"]["acs1"]
        toml_section(L, "housing_demand.acs_check", dict(alpha_income=a1["fit"]["alpha"], h_share_of_mean_income=a1["fit"]["h_share_of_mean_income"],
                                                         share_cost_burdened_30plus=a1["burden"]["share_30plus_all"],
                                                         mean_burden_by_income_bin=[r["mean_burden"] for r in a1["burden"]["rows"]],
                                                         note="ACS 2018 1-year, Stone-Geary form on household income (no size adjustment), mean burdens by 7 income bins"))
    if "eps" in rec:
        e = rec["eps"]
        toml_section(L, "housing_supply.eps_city", dict(value=e["eps_city"]["value"], lo=e["eps_city"]["lo"], hi=e["eps_city"]["hi"],
                                                        no_horizon_conversion=e["no_horizon_conversion"]["city"],
                                                        note="long-run elasticity of housing supply with respect to rent, central city (one horizon convention, delta = 3.5%)"))
        toml_section(L, "housing_supply.eps_suburb", dict(value=e["eps_suburb"]["value"], lo=e["eps_suburb"]["lo"], hi=e["eps_suburb"]["hi"],
                                                          no_horizon_conversion=e["no_horizon_conversion"]["suburb"],
                                                          note="long-run elasticity of housing supply with respect to rent, suburb (one horizon convention, delta = 3.5%)"))
        toml_section(L, "housing_supply.derivation", dict(
            ratio_suburb_city=e["ratio_suburb_city"], ratio_lo=e["ratio_range"][0], ratio_hi=e["ratio_range"][1],
            ratio_lo_matched_split=e["ratio_range_matched_split"][0], ratio_hi_matched_split=e["ratio_range_matched_split"][1],
            city_share=e["city_share"], metro_anchor=e["anchor_central"], metro_anchor_lo=e["anchor_low"], metro_anchor_hi=e["anchor_high"],
            depreciation=e["delta"], longrun_factor_30y=e["longrun_factor_30y"], longrun_factor_10y=e["longrun_factor_10y"],
            ten_year_units_city=e["short_run_10yr_units"]["city"], ten_year_units_suburb=e["short_run_10yr_units"]["suburb"],
            ten_year_space_city=e["short_run_10yr_space"]["city"], ten_year_space_suburb=e["short_run_10yr_space"]["suburb"],
            note="eps_c = anchor / (s + (1 - s) ratio), eps_s = ratio x eps_c; anchor = Saiz population-weighted mean over 95 MSAs x 30-year factor (low: Baum-Snow and Han tracts x 10-year factor; high: Saiz unweighted mean over 269 MSAs x 30-year factor)"))
        for k, a_ in e["anchors"].items():
            items = dict(years=a_["years"], role=a_["role"], label=a_["label"])
            if "published" in a_:
                items["published"] = a_["published"]
            if "published_range" in a_:
                items["published_range"] = a_["published_range"]
            if "long_run" in a_:
                items["long_run"] = a_["long_run"]
            if "long_run_range" in a_:
                items["long_run_range"] = a_["long_run_range"]
            toml_section(L, f"housing_supply.anchors.{k}", items)
        pc = e["period_check"]
        toml_section(L, "housing_supply.period_check", dict(
            saiz_1970_2000_mean=res["c"]["appendix_d"]["saiz_1970_2000_mean"], saiz_method_2000_2010_mean=res["c"]["appendix_d"]["saiz_method_2000_2010_mean"],
            observed_ratio=pc["observed_ratio"], predicted_ratio_delta35=pc["predicted_ratio"], implied_delta=pc["implied_delta"],
            table_a8_means_2000_2010=pc["table_a8_means_2000_2010"], table_a8_means_1970_2000=pc["table_a8_means_1970_2000"],
            note="Baum-Snow and Han, Online Appendix D (pp. 17-18) and Table A8 (p. 26): Saiz's method on 2000-10 data"))
        ds = e["delta_sensitivity"]
        toml_section(L, "housing_supply.delta_sensitivity", dict(
            delta=[r_["delta"] for r_ in ds], factor_30y=[r_["factor_30y"] for r_ in ds], factor_10y=[r_["factor_10y"] for r_ in ds],
            metro_anchor=[r_["metro"] for r_ in ds], eps_city=[r_["city"] for r_ in ds], eps_suburb=[r_["suburb"] for r_ in ds],
            bsh_tracts_metro=[r_["bsh_tracts_metro"] for r_ in ds]))
        sz = res["c"]["saiz"]
        toml_section(L, "housing_supply.saiz", dict(pop_weighted_mean_95=sz["msas_over_500k"]["pop_weighted_mean"], mean_95=sz["msas_over_500k"]["mean"],
                                                    p25_95=sz["msas_over_500k"]["p25"], p75_95=sz["msas_over_500k"]["p75"], pop_weighted_mean_269=sz["all_msas"]["pop_weighted_mean"],
                                                    mean_269=sz["all_msas"]["mean"], horizon_years=SAIZ_YEARS))
    if "r_car" in rec:
        r = rec["r_car"]
        d = res["d"]
        cr = d["car_requirement"]
        toml_section(L, "transport.r_car", dict(value=r["value"], se=r["se"], lo=r["lo"], hi=r["hi"], extensive_margin=r["extensive_margin"],
                                                single_adult=r["single_adult"], single_adult_se=r["single_adult_se"], first_version_no_size_adjustment=r["first_version"],
                                                unadjusted_upper=r["unadjusted_upper"], annual_usd=r["value"] * CE_TOTAL_EXP, note=r["basis"]))
        for k, v in cr["extra_vehicles"].items():
            toml_section(L, f"transport.extra_vehicles.{k}", dict(value=v["value"], se=v["se"], label=cr["extra_labels"][k]))
        hs = d["household_size"]
        toml_section(L, "transport.household_size", dict(
            classes=hs["classes"], national_class_shares=hs["national_class_shares"],
            share_city=[hs["composition"][c_]["level_1"] for c_ in hs["classes"]], share_suburb=[hs["composition"][c_]["level_2"] for c_ in hs["classes"]],
            vehicles_city=[hs["by_class"][c_]["vehicles"]["level_1"] for c_ in hs["classes"]], vehicles_suburb=[hs["by_class"][c_]["vehicles"]["level_2"] for c_ in hs["classes"]],
            vehicles_gap=[hs["by_class"][c_]["vehicles"]["value"] for c_ in hs["classes"]], vehicles_gap_se=[hs["by_class"][c_]["vehicles"]["boot_se"] for c_ in hs["classes"]],
            vehicles_gap_income_adjusted=[hs["regressions"][f"vehicles|{c_}|log_income"]["beta"][0] for c_ in hs["classes"]],
            vehicles_gap_income_adjusted_se=[hs["regressions"][f"vehicles|{c_}|log_income"]["se_boot"][0] for c_ in hs["classes"]],
            four_plus_vehicles_city=hs["four_plus_vehicles_city"], four_plus_vehicles_suburb=hs["four_plus_vehicles_suburb"],
            standardized_gap=hs["standardized_gap_pooled"]["value"], standardized_gap_se=hs["standardized_gap_pooled"]["boot_se"],
            standardized_gap_income_adjusted=hs["standardized_regression_vehicles"]["income_adjusted"]["value"],
            standardized_gap_income_adjusted_se=hs["standardized_regression_vehicles"]["income_adjusted"]["boot_se"],
            note="B08201 household-size classes; suburb minus city; standardized = within-class gaps averaged with household-size shares"))
        for k in ("zero_vehicle_share", "vehicles_per_hh", "zero_vehicle_share_owners", "zero_vehicle_share_renters", "vehicles_per_owner_hh", "vehicles_per_renter_hh"):
            v = d["stats"][k]
            toml_section(L, f"transport.locations.{k}", dict(city=v["level_1"], suburb=v["level_2"], gap=v["value"], se=v["boot_se"], n_cz=v["n"]))
        for k in ("veh_per_hh|households|no_controls", "veh_per_hh|households|log_income", "veh_per_hh|households|log_income_size_shares",
                  "zero_share|households|no_controls", "zero_share|households|log_income", "zero_share|households|log_income_size_shares"):
            v = d["regressions"][k]
            toml_section(L, "transport.regressions." + k.replace("|", "_"), dict(suburb=v["beta"][0], se_boot=v["se_boot"][0], se_cluster=v["se_cluster"][0], n=v["n"]))
        toml_section(L, "transport.cost_per_vehicle_year", dict(d["cost_per_vehicle_year"], ce_vehicles_rounding_bound=d["ce_vehicles_rounding"]["relative_bound"]))
    if "wtp_income_elasticity" in rec:
        w = rec["wtp_income_elasticity"]
        toml_section(L, "wtp_income.elasticity", dict(value=w["value"], se=w["se"], se_with_mean_mwtp=w["se_with_mean_mwtp"],
                                                      note="BFM (2007) Table 8: 1.38 (0.33) per $10,000 at mean income $54,103 and mean MWTP $19.69"))
        toml_section(L, "wtp_income.elasticity_with_college", dict(value=w["with_college"], se=w["with_college_se"], lo=w["with_college_lo"], hi=w["with_college_hi"],
                                                                  note="adds 13.03 x slope of the householder college share (bachelor's or more) on income, 1990 Bay Area PUMS; lo / hi over the slope measures and the variants"))
        va = res["e"]["elasticity_variants"]
        for k in ("bachelors", "associate_plus", "bachelors_children", "associate_plus_children"):
            toml_section(L, f"wtp_income.variants.{k}", dict(value=va[k]["elasticity"], se=va[k]["se"], slope_per_10k=va[k]["slope_per_10k"]))
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

def cmd_build(args: argparse.Namespace) -> int:
    indir = Path(args.indir).resolve()
    default_out = indir / "estimates" / "settings"
    if args.outdir is None and (args.B != B_DEFAULT or args.seed != SEED):
        raise SystemExit("non-default B or seed: pass --outdir so the tracked outputs are not overwritten")
    outdir = Path(args.outdir).resolve() if args.outdir else default_out
    outdir.mkdir(parents=True, exist_ok=True)
    parts = set(args.parts)
    t0 = dt.datetime.now()
    panel = load_panel(indir)
    part, g = cz_table(panel, "locale")
    pup = part.groupby("loc_id").pupils.sum()
    city_share = float(pup[1] / (pup[1] + pup[2]))
    LOG.info("locale partition: %d CZs, %d districts; city share of partitioned pupils %.4f", len(g), len(part), city_share)
    boot = Boot(len(g), args.B, args.seed, g.pupils.to_numpy(float))
    rng = np.random.default_rng(args.seed)
    res: dict[str, Any] = {}
    draws: dict[str, np.ndarray] = {}
    if "a" in parts:
        LOG.info("(a) housing-cost distributions")
        res["a"] = part_a(part, g, boot, indir)
        for tab, nm in (("all", "all"), ("rent", "rent"), ("own", "own")):
            for k in PCT_KEYS:
                draws[f"loggap_{nm}_{k}"] = res["a"]["tables"][tab]["percentiles"][k]["log_gap"]["draws"]
        es = res["a"]["reproduction"].get("earnings_stat")
        if es is not None:
            draws["gap_income"] = es["draws"]
    ce = None
    if parts & {"b", "d"}:
        LOG.info("(b) CE tables and Stone-Geary fits")
        ce = ce_block(RAW / "cu-income-quintiles-before-taxes-2018.xlsx", RAW / "cu-income-deciles-before-taxes-2018.xlsx", rng)
        # the equivalized fits use their own random stream, so the CU-level numbers and the later draws of `rng` are unchanged
        eq = ce_equivalized(RAW / "cu-income-quintiles-before-taxes-2018.xlsx", RAW / "cu-income-deciles-before-taxes-2018.xlsx",
                            np.random.default_rng([args.seed, 11]))
        size = ce_size_table(RAW / "cu-size-2018.xlsx")
        res["b"] = dict(ce=ce, acs=acs_block(), equivalized=eq, size_table=size,
                        size_check=size_check(size, ce["quintiles"]["fits"]["shelter_utilities"],
                                              eq["quintiles"]["fits"][CENTRAL_SCALE]["shelter_utilities"]))
    if "c" in parts:
        LOG.info("(c) supply elasticities")
        res["c"] = part_c(city_share)
        res["c"]["_city_share"] = city_share
    if "d" in parts:
        LOG.info("(d) transport")
        dd = part_d(panel, part, g, boot, ce, indir, size_table=res["b"]["size_table"])
        dd["car_requirement"] = car_requirement(dd)
        dd["aaa_verification"] = verify_aaa(RAW / "AAA_Your_Driving_Costs_2018.pdf")
        dd["aaa_2018"] = AAA_2018
        res["d"] = dd
        for k in ("zero_vehicle_share", "vehicles_per_hh"):
            draws[f"gap_{k}"] = dd["stats"][k]["draws"]
        hs_ = dd["household_size"]
        draws["gap_vehicles_size_standardized"] = hs_["standardized_gap_pooled"]["draws"]
        draws["gap_vehicles_one_person"] = hs_["by_class"][SIZE_LABELS[1]]["vehicles"]["draws"]
        draws["gap_vehicles_size_standardized_income_adjusted"] = hs_["standardized_regression_vehicles"]["income_adjusted"]["draws"]
    if "e" in parts:
        LOG.info("(e) BFM gradient and the 1990 Bay Area college / income association")
        res["e"] = part_e(RAW / "PUMSAXCA.zip", rng)
    res["recommended"] = recommendations(res)
    meta = dict(generated=t0.isoformat(timespec="seconds"), script="data/spatial/housing_transport.py", year=YEAR, B=args.B, seed=args.seed,
                n_cz=int(len(g)), n_districts=int(len(part)), city_share_of_pupils=city_share, parts=sorted(parts),
                rng="numpy default_rng (PCG64) via district_bootstrap.draw_counts; one (B x n_CZ) multiplicity matrix shared by every statistic; "
                    "resampling unit = whole CZ; CZ weights = partitioned pupils",
                ce_mean_expenditures=CE_TOTAL_EXP)
    payload = clean({"meta": meta, **res})
    (outdir / "housing_transport.json").write_text(json.dumps(payload, indent=1))
    md = write_md(res, meta)
    (outdir / "housing_transport.md").write_text(md)
    (outdir / "housing_transport_targets.toml").write_text(toml_text(res, meta))
    if draws:
        pd.DataFrame(draws).round(8).to_csv(outdir / "housing_transport_draws.csv", index_label="draw")
    LOG.info("wrote housing_transport.{md,json}, housing_transport_targets.toml and housing_transport_draws.csv to %s", outdir)
    return 0


def write_md(res: dict[str, Any], meta: dict[str, Any]) -> str:
    L = ["# Housing market and transport evidence for income-driven demand (task D2, DESIGN.md M1 / M1b)\n"]
    L.append(f"Generated {meta['generated']} by `housing_transport.py`. ACS 5-year 2014-18 school-district tables (dollars of 2018), merged to the "
             "district panel and the locale partition of `district_bootstrap.py` (city = NCES locale 11-13, suburb = 21-23; base year 2018, at least "
             f"500 pupils per district, each location at least 10% of the CZ's partitioned pupils; {meta['n_cz']} CZs, {meta['n_districts']} districts). "
             "Location statistics pool the districts of a location within a CZ; gaps are suburb minus city per CZ and are averaged over CZs with the "
             f"partitioned pupils as weights; whole-CZ bootstrap, B = {meta['B']}, seed {meta['seed']}: the draws of `housing_costs.py` and "
             "`district_settings.py`, shared by every statistic here. CE = BLS Consumer Expenditure Survey, calendar year 2018 (mean annual "
             f"expenditures per consumer unit ${CE_TOTAL_EXP:,.2f}, the unit of the model's Cbar).\n")
    rec = res.get("recommended", {})
    if {"alphaH", "eps", "r_car"} <= set(rec):
        a_, h_, e_, r_ = rec["alphaH"], rec["hmin"], rec["eps"], rec["r_car"]
        L.append("**Revision after review (2026-10-05).** Compared with the first version of this report: (b) the central `alphaH` / `hmin` come from the fit on "
                 f"expenditure per equivalent adult (square-root scale), {fmt(a_['value'], 3)} / {fmt(h_['value'], 4)} against {fmt(a_['cu_level'], 3)} / "
                 f"{fmt(h_['cu_level'], 4)} at the consumer-unit level, with the ranges widened to {fmt(a_['lo'], 2)}-{fmt(a_['hi'], 2)} / {fmt(h_['lo'], 2)}-{fmt(h_['hi'], 2)}; "
                 "`housing_share` is defined as a ratio of aggregates and the recommended terciles come from the deciles; (c) eps_H applies one horizon convention to every "
                 f"anchor: city {fmt(e_['eps_city']['value'], 2)} / suburb {fmt(e_['eps_suburb']['value'], 2)} (ranges {fmt(e_['eps_city']['lo'], 2)}-{fmt(e_['eps_city']['hi'], 2)} / "
                 f"{fmt(e_['eps_suburb']['lo'], 2)}-{fmt(e_['eps_suburb']['hi'], 2)}) against {fmt(e_['no_horizon_conversion']['city'], 2)} / "
                 f"{fmt(e_['no_horizon_conversion']['suburb'], 2)}; (d) r_car is standardized for household size, {fmt(r_['value'], 4)} ({fmt(r_['lo'], 4)}-{fmt(r_['hi'], 4)}) against "
                 f"{fmt(r_['first_version'], 4)}; (a) the validation moments for the lower tail are the renters' (the bottom of the all-units distribution is mostly outright owners); "
                 "(e) associate-or-more and children-term variants are added. The partition, the bootstrap draws and every other statistic are unchanged.\n")
    if "recommended" in res and res["recommended"]:
        L += md_summary(res)
        L += md_recommended(res)
    if "a" in res:
        L += md_a(res)
    if "b" in res:
        L += md_b(res)
    if "c" in res:
        L += md_c(res)
    if "d" in res:
        L += md_d(res)
    if "e" in res:
        L += md_e(res)
    L += md_caveats(res)
    L += md_sources()
    return "\n".join(L) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pf = sub.add_parser("fetch", help="download and cache every external input under raw/housing_transport/")
    pf.add_argument("--refresh", action="store_true")
    pf.add_argument("--workers", type=int, default=8)
    pb = sub.add_parser("build", help="compute parts (a)-(e) and write the report, JSON and TOML")
    pb.add_argument("--indir", default=str(HERE))
    pb.add_argument("--outdir", default=None, help="default: <indir>/estimates/settings")
    pb.add_argument("--B", type=int, default=B_DEFAULT)
    pb.add_argument("--seed", type=int, default=SEED)
    pb.add_argument("--parts", nargs="+", default=list("abcde"), choices=list("abcde"))
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    return cmd_fetch(args) if args.cmd == "fetch" else cmd_build(args)


if __name__ == "__main__":
    raise SystemExit(main())
