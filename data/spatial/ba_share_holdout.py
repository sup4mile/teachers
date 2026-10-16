#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ba_share_holdout.py -- SEALED hold-out moment B1 of the 2026-10-15 round: the within-commuting-zone gap in the
share of adults 25+ with a bachelor's degree or more, suburb minus city, ACS 5-year 2014-18 by school district.

B1 is the data counterpart of the model's `ba_share_gap` (BA share of adults by residence, suburb minus city; M3
only).  It is pre-registered as "Sealed holdout B1" in julia/spatial_model/calibration/runs/experiments-2026-10d/
PREREG.md, was never used before, and is computed here before any result of the round is selected.  The value
goes ONLY to estimates/settings/sealed/ba_share_holdout_B1.json; it is opened after the specification is chosen
(PREREG A7) and no specification is re-selected on it.

What `ba_holders` measures.  `fetch_data.ACS_SD_VARS` maps B15003_022E to `ba_holders` and B15003_001E to
`pop_25plus`.  B15003 is "Educational attainment for the population 25 years and over", one line per HIGHEST level
completed; line 22 is "Bachelor's degree", and the master's (23), professional-school (24) and doctorate (25)
degrees are separate lines.  `ba_holders / pop_25plus` of `raw/acs_school_districts.parquet` is therefore the
share whose highest degree is a bachelor's, NOT the share with a bachelor's or more.  The pre-registered concept
(and the model's college graduate) is lines 22-25, so `fetch` pulls B15003 lines 1-25 for the same vintage and the
same three school-district layers (and the state totals, for a tiling check); the exactly-bachelor's share
computed from the same lines is the variant `ba_only`.  The pull is checked line by line against the existing
parquet (same records, same `pop_25plus`, same `ba_holders`) and against the Census metadata labels.

Definition (headline, fixed before any value was computed; methods file section 3)
    Panel and partition: `district_settings.build_audited_panel` (New York City one district under the
    Chancellor's id 3620580, under which the ACS files it) and `scheme_tables(A, "locale")`: city = NCES locale
    11-13, suburb = 21-23, districts of >= 500 pupils, each location >= 10% of the CZ's partitioned pupils: the
    188 CZs of the four district targets.  CZ weights = partitioned pupils.
    Per CZ and location l: s_l = (sum of BA+ over the location's districts) / (sum of adults 25+), BA+ = lines
    22-25 of B15003, a ratio of location totals as for every per-pupil object of `data_estimate.Gap`; the CZ
    gap is s_2 - s_1 (suburb minus city).  B1 = the CZ-weighted mean of the CZ gaps (a CZ with no valid record in
    one location drops out), reported as a fraction and in percentage points; also the log ratio
    log(s_2/s_1), the CZ-weighted city and suburb shares, and the version pooled over CZs (adult-weighted:
    the ratio of sums over all CZs, not differenced within a CZ).
    Overlapping layers: an adult who lives in an elementary district inside a secondary district appears in both
    records.  The headline counts adults once, from the unified and elementary records, which tile each state
    (the rule of `teacher_residence.py` for employed residents and of `children_income.py` for children); the
    secondary LEAs keep their pupils in the CZ weights and contribute no adults.  The rule of the earnings
    targets (`district_bootstrap.load_income_columns`: every LEA takes the record of its own layer) is the
    variant `all_layers`; the others are in SPECS.
    Matching: ACS record to LEA on the 7-digit leaid; the "Remainder of <state>" rows (code 99999) are dropped;
    Detroit's record, filed under the pre-2016 id 2612000, is matched to 2601103 as in `children_income.py`.
    Inference: whole-CZ bootstrap, B = 1999, seed 20260929, the multiplicity matrix of `district_settings.py`
    (`district_bootstrap.draw_counts(188, B, seed)`), so the draws are the draws of the four district targets and
    the bootstrap correlations with them are valid.  The four targets (0.3764, 0.0086, -0.0177, 0.1122) and
    their SEs are reproduced, and the draws compared with `district_settings_draws_locale.csv`; the script
    stops if either fails.
    Hold-out halves: the pre-registered split of the locale CZs (`district_settings.holdout_split`, seed
    20261002) is reproduced and checked against `district_settings.json`; B1 is reported on half A and half B
    (each with the multiplicity matrix of its own n CZs, as `district_settings.py` does).

SEALING RULES.  No value, standard error, share or correlation of B1 or of any variant is printed, logged or
returned: they are written ONLY with json.dump to estimates/settings/sealed/ba_share_holdout_B1.json.  The
terminal, the logs and the methods file show sample sizes, coverage, the reproduction of the four public targets
and PASS/FAIL results of checks that involve no sealed value.  `seal` prints the public coverage tables and the
check results, then the counts of the leak scan of the methods file, `sealed: <path> (<k> statistics)` and the
SHA-256 of the file; it refuses to overwrite an existing sealed file without `--force`, and deletes the methods
file (and exits 1) if the leak scan finds a sealed number in it.  `leak-scan FILE...` runs the same scan on other
texts (the run logs).  The scan looks for high-precision renderings of the sealed numbers (fractions at 5-6
decimals, percentage points at 4: none may appear) and compares the count of coarse renderings (4 decimals, 2-3 for
percentage points: these collide with public numbers by chance) with the counts for random decoys of the sealed
numbers.  `ba_share_holdout_xcheck.py` recomputes the sealed statistics from the raw tables in separate code and
prints only PASS/FAIL and maximum absolute differences.  Do not open the sealed file before the specification is
chosen.

Usage (from data/spatial; every step but `selftest` also runs on a compute node: see below)
    python3 ba_share_holdout.py fetch               # B15003 by school district + state totals (Census API) -> raw/
    python3 ba_share_holdout.py check               # coverage and checks that involve no sealed value
    python3 ba_share_holdout.py seal                # B = 1999 -> sealed json + methods file
    python3 ba_share_holdout.py verify-seal         # structure of the sealed file, no values
    python3 ba_share_holdout.py leak-scan FILE ...  # scan texts (logs) for renderings of the sealed numbers
    python3 ba_share_holdout.py selftest            # unit checks of the helpers (no data, no network)
    python3 ba_share_holdout_xcheck.py              # independent recomputation, PASS/FAIL only
Reads raw/acs_school_districts.parquet, raw/acs_b15003_*.parquet (from `fetch`) and the tables of the
calibration stack (`district_settings.build_audited_panel`); reads estimates/settings/district_settings.json and
district_settings_draws_locale.csv (unsealed) for the checks.  Writes estimates/settings/sealed/
ba_share_holdout_B1.json and estimates/settings/ba_share_holdout_methods.md (no values).  The Census API key
(CENSUS_API_KEY, from the environment or data/spatial/.env) is used by `fetch` only and is never printed or
written; error strings are scrubbed of it.
On Slurm (srun, econ-grad): `--mem` is enforced as an address-space limit that pyarrow's per-thread jemalloc
arenas exceed (segfaults, `ArrowMemoryError: malloc failed`); run with ARROW_DEFAULT_MEMORY_POOL=system and
OMP_NUM_THREADS equal to the allocated CPUs.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import logging
import os
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

LOG = logging.getLogger("ba_share_holdout")
de = da = db = ds = None          # the calibration stack, imported by `_stack()` (not needed by `fetch`)


def _stack() -> None:
    global de, da, db, ds
    if ds is None:
        import data_estimate as _de
        import district_audit as _da
        import district_bootstrap as _db
        import district_settings as _ds
        de, da, db, ds = _de, _da, _db, _ds


RAW = HERE / "raw"
EST = HERE / "estimates" / "settings"
YEAR = 2018                                   # ACS 5-year vintage 2014-18
SEED = 20260929                               # = district_bootstrap.SEED
B_DEFAULT = 1999                              # = district_bootstrap.B_DEFAULT
HOLDOUT_SEED = 20261002                       # = district_settings.HOLDOUT_SEED
OLD_FILE = RAW / "acs_school_districts.parquet"         # fetch_data.py: ba_holders = B15003_022E, pop_25plus = B15003_001E
B15003_PATH = RAW / f"acs_b15003_school_districts_{YEAR}.parquet"
B15003_ST_PATH = RAW / f"acs_b15003_states_{YEAR}.parquet"
B15003_LABELS = RAW / f"acs_b15003_labels_{YEAR}.json"
SEALED_NAME = "ba_share_holdout_B1.json"
METHODS_NAME = "ba_share_holdout_methods.md"

N_LINES = 25                                  # B15003_001 (total) ... _025 (doctorate)
B15003_VARS = [f"B15003_{j:03d}E" for j in range(1, N_LINES + 1)]
# lines the construction relies on: the label the Census metadata must carry (last level of "Estimate!!Total:!!...")
LINE_LABEL = {22: "Bachelor's degree", 23: "Master's degree", 24: "Professional school degree",
              25: "Doctorate degree"}
NYC_LEAID = "3620580"                         # the Chancellor's id: F-33, SAIPE, SEDA and the ACS file New York City under it
# Detroit's ACS record is filed under the pre-2016 id; the CCD frame carries Detroit Public Schools Community
# District (2601103).  Source id -> panel id, as children_income.DETROIT.
DETROIT = ("2612000", "2601103")
TARGETS = ["gap_pupils", "gap_salary", "gap_teachers_pp", "gap_income"]     # order of district_settings.TARGET_KEYS
SAMPLES = ("full", "A", "B")
STAT_KEYS = ("ba_share_gap", "ba_share_gap_pp", "ba_log_ratio", "ba_share_city", "ba_share_suburb",
             "pooled_gap", "pooled_city", "pooled_suburb")
DRAW_KEYS = ("ba_share_gap", "ba_log_ratio", "pooled_gap")
DRAW_STAT = {k: k for k in DRAW_KEYS}         # the statistic each stored draw vector belongs to
COVERAGE_MIN = 0.90                           # sanity floor: share of CZs / of the CZ weight with a valid headline gap


# --------------------------------------------------------------------------
# Fetch: ACS B15003 by school district (Census API)
# --------------------------------------------------------------------------

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


def _norm_label(lab: str) -> str:
    """Last level of a Census variable label ('Estimate!!Total:!!Bachelor's degree' -> "Bachelor's degree")."""
    return lab.split("!!")[-1].strip().rstrip(":").strip()


def check_labels(labels: dict[str, str]) -> None:
    """The Census metadata labels of the lines the construction relies on (stop if the lines moved)."""
    if _norm_label(labels.get("B15003_001E", "")) != "Total":
        raise SystemExit(f"B15003_001E is not the total: {labels.get('B15003_001E')!r}")
    for j, lab in LINE_LABEL.items():
        got = labels.get(f"B15003_{j:03d}E", "")
        if _norm_label(got) != lab:
            raise SystemExit(f"B15003_{j:03d}E is {got!r}, expected {lab!r}")


def _api_json(sess: Any, url: str, key: str, tries: int = 4) -> Any:
    """GET a Census API url with the key attached as a parameter; every message that could contain the key is
    scrubbed before it leaves this function."""
    last = ""
    for k in range(tries):
        try:
            r = sess.get(url, params={"key": key}, timeout=300)
            if r.status_code == 200:
                return r.json()
            last = f"http {r.status_code}"
        except Exception as exc:                                  # noqa: BLE001
            last = f"{type(exc).__name__}: {str(exc).replace(key, '<key>')[:160]}"
        time.sleep(2 * (k + 1))
    raise SystemExit(f"Census API request failed ({url.split('/data/')[-1].split('?')[0]}): {last}")


def fetch_b15003(workers: int = 8, refresh: bool = False) -> Path:
    """B15003 lines 1-25 for the unified, elementary and secondary school-district layers, one state per request
    as `acs_earnings.py` does, the state totals of the same lines, and the Census labels of the lines;
    raw/acs_b15003_{school_districts,states,labels}_<year>."""
    if B15003_PATH.exists() and B15003_ST_PATH.exists() and B15003_LABELS.exists() and not refresh:
        LOG.info("cached   %s", B15003_PATH.name)
        return B15003_PATH
    import acs_earnings as ae
    key = _read_key()
    sess = ae.make_session()
    base = f"https://api.census.gov/data/{YEAR}/acs/acs5"
    var = _api_json(sess, f"{base}/groups/B15003.json", key)["variables"]
    labels = {v: var[v]["label"] for v in B15003_VARS}
    check_labels(labels)
    LOG.info("labels: B15003_001E total; lines 022-025 = %s", [_norm_label(labels[f"B15003_{j:03d}E"])
                                                              for j in range(22, 26)])
    states = [row[-1] for row in _api_json(sess, f"{base}?get=NAME&for=state:*", key)[1:]]
    out = ae.pull(sess, base, states, key, ",".join(["NAME"] + B15003_VARS), workers, "B15003")
    if out.empty:
        raise SystemExit("B15003: nothing returned")
    out["year"] = YEAR
    for c in B15003_VARS:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    RAW.mkdir(parents=True, exist_ok=True)
    out.to_parquet(B15003_PATH, index=False)
    LOG.info("wrote %s (%d rows); layers: %s", B15003_PATH.name, len(out), out.sd_layer.value_counts().to_dict())
    # state totals of the same lines: the check that unified + elementary records tile each state
    j = _api_json(sess, f"{base}?get=NAME,{','.join(B15003_VARS)}&for=state:*", key)
    st = pd.DataFrame(j[1:], columns=j[0])
    for c in B15003_VARS:
        st[c] = pd.to_numeric(st[c], errors="coerce")
    st.to_parquet(B15003_ST_PATH, index=False)
    B15003_LABELS.write_text(json.dumps(dict(vintage=f"ACS 5-year {YEAR - 4}-{YEAR}", table="B15003",
                                             fetched=dt.datetime.now().isoformat(timespec="seconds"),
                                             labels=labels), indent=1))
    LOG.info("wrote %s (%d states) and %s", B15003_ST_PATH.name, len(st), B15003_LABELS.name)
    return B15003_PATH


def cmd_fetch(args: argparse.Namespace) -> int:
    fetch_b15003(args.workers, args.refresh)
    return 0


# --------------------------------------------------------------------------
# ACS records: one row per leaid
# --------------------------------------------------------------------------

def layer_of(sd_layer: pd.Series) -> pd.Series:
    """'school district (unified)' -> 'unified'."""
    return sd_layer.str.extract(r"\((\w+)\)")[0]


def records_from_raw(a: pd.DataFrame) -> pd.DataFrame:
    """B15003 records of the three layers (columns NAME, B15003_001E..025E, leaid, sd_layer) to one row per leaid:
        'Remainder of <state>' rows (code 99999) are not districts and are dropped, one row per leaid (the
        three layers otherwise never share an id), annotation codes (negative) are masked;
        adults   B15003_001E (a record needs at least one adult)
        ba_plus  lines 22-25 (bachelor's, master's, professional, doctorate: highest degree a bachelor's or more)
        ba_only  line 22 (highest degree exactly a bachelor's: the old `ba_holders`)
    A record is valid if it has at least one adult, the adult count and the four degree lines are published and
    BA+ does not exceed the adults; an invalid record has NaN in all three columns (it leaves every sum).  In the
    2018 pull the invalid records are exactly the empty districts (0 adults 25+: "School District Not Defined",
    unorganized territories, a few elementary and secondary districts with no residents)."""
    n_all = len(a)
    a = a[~a.leaid.str.endswith("99999")]
    n_remainder = n_all - len(a)
    n_dup = int(a.leaid.duplicated().sum())
    a = a.drop_duplicates("leaid")

    def cell(j: int) -> pd.Series:
        return a[f"B15003_{j:03d}E"].astype(float)

    adults = cell(1).where(cell(1) >= 1)
    top = pd.concat([cell(j).where(cell(j) >= 0) for j in (22, 23, 24, 25)], axis=1)
    ba_plus = top.sum(axis=1, min_count=4)
    ba_only = cell(22).where(cell(22) >= 0)
    ok = adults.notna() & ba_plus.notna() & (ba_plus <= adults)
    out = pd.DataFrame({"leaid": a.leaid.to_numpy(), "name": a.NAME.to_numpy(),
                        "layer": layer_of(a.sd_layer).to_numpy(), "adults": adults.where(ok).to_numpy(),
                        "ba_plus": ba_plus.where(ok).to_numpy(), "ba_only": ba_only.where(ok).to_numpy()})
    out.attrs.update(n_records=int(n_all), n_remainder_dropped=int(n_remainder), n_duplicates_dropped=n_dup,
                     n_invalid=int((~ok).sum()), n_zero_adults=int((cell(1) < 1).sum()))
    return out


def load_b15003(path: Path = B15003_PATH) -> pd.DataFrame:
    if not path.exists():
        raise SystemExit(f"no {path.name}: run `ba_share_holdout.py fetch`")
    return records_from_raw(pd.read_parquet(path))


def load_old_file(path: Path = OLD_FILE) -> pd.DataFrame:
    """`raw/acs_school_districts.parquet` (fetch_data.py) with its leaid, layer, `pop_25plus` and `ba_holders`."""
    o = pd.read_parquet(path)
    lay = layer_of(o.sd_layer)
    code = pd.Series(np.where(lay == "unified", o["school district (unified)"],
                              np.where(lay == "elementary", o["school district (elementary)"],
                                       o["school district (secondary)"])), index=o.index)
    return pd.DataFrame({"leaid": (o.state.astype(str).str.zfill(2) + code.astype(str).str.zfill(5)).to_numpy(),
                         "sd_layer": o.sd_layer.to_numpy(), "name": o.NAME.to_numpy(),
                         "pop_25plus": pd.to_numeric(o.pop_25plus, errors="coerce").to_numpy(),
                         "ba_holders": pd.to_numeric(o.ba_holders, errors="coerce").to_numpy()})


def pull_vs_old(raw_new: pd.DataFrame, old: pd.DataFrame) -> dict[str, Any]:
    """The new pull against the old file: same records, line 1 = `pop_25plus`, line 22 = `ba_holders`, and the
    internal consistency of B15003 (lines 2-25 add up to line 1; no annotation codes)."""
    m = raw_new.merge(old, on=["leaid", "sd_layer"], how="outer", indicator=True, suffixes=("", "_old"))
    both = m[m["_merge"] == "both"]
    lines = [f"B15003_{j:03d}E" for j in range(1, N_LINES + 1)]
    inner = raw_new[[f"B15003_{j:03d}E" for j in range(2, N_LINES + 1)]].sum(axis=1, min_count=N_LINES - 1)
    return dict(n_new=int(len(raw_new)), n_old=int(len(old)), n_both=int(len(both)),
                n_only_new=int((m["_merge"] == "left_only").sum()), n_only_old=int((m["_merge"] == "right_only").sum()),
                n_equal_pop=int((both.B15003_001E == both.pop_25plus).sum()),
                n_equal_ba=int((both.B15003_022E == both.ba_holders).sum()),
                n_lines_sum_to_total=int((inner == raw_new.B15003_001E).sum()),
                n_with_negative_code=int((raw_new[lines] < 0).any(axis=1).sum()),
                n_with_missing=int(raw_new[lines].isna().any(axis=1).sum()))


def tiling_check(raw_new: pd.DataFrame, states: pd.DataFrame) -> dict[str, Any]:
    """Do the unified and elementary records tile the states?  Their adults 25+ against the state totals of the
    same line, and with the secondary layer added (the overlap of the earnings rule)."""
    a = raw_new[~raw_new.leaid.str.endswith("99999")].drop_duplicates("leaid").copy()
    a["layer"] = layer_of(a.sd_layer)
    tot = a.groupby(["state", "layer"]).B15003_001E.sum().unstack("layer").fillna(0.0)
    st = states.set_index("state")
    tab = tot.join(st[["NAME", "B15003_001E"]].rename(columns={"B15003_001E": "state_total"}))
    tab["ue"] = (tab.unified + tab.elementary) / tab.state_total
    tab["all"] = (tab.unified + tab.elementary + tab.secondary) / tab.state_total
    nat = float(st.B15003_001E.sum())
    top = tab.sort_values("all", ascending=False).head(5)
    return dict(n_states=int(len(tab)), national_adults_state_totals=nat,
                national_unified_plus_elementary=float((tot.unified.sum() + tot.elementary.sum()) / nat),
                national_all_layers=float((tot.unified.sum() + tot.elementary.sum() + tot.secondary.sum()) / nat),
                n_states_ue_within_1pct=int(((tab.ue - 1.0).abs() <= 0.01).sum()),
                n_states_ue_within_5pct=int(((tab.ue - 1.0).abs() <= 0.05).sum()),
                min_ue=float(tab.ue.min()), max_ue=float(tab.ue.max()),
                top_states_all_layers=[dict(state=str(r.NAME), all_layers=float(r["all"]),
                                            unified_plus_elementary=float(r.ue)) for _, r in top.iterrows()])


# --------------------------------------------------------------------------
# Specifications
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Spec:
    key: str
    label: str
    layers: str = "ue"             # ue unified + elementary records | all every record | us unified + secondary | u unified only
    measure: str = "ba_plus"       # ba_plus lines 22-25 | ba_only line 22
    nyc: bool = True               # False: New York City's adults leave the sums (its pupils stay in the CZ weights)
    detroit: bool = True           # False: Detroit's record is not refiled to 2601103 (no id crosswalk)
    ny_cz: bool = True             # False: the New York CZ leaves the CZ mean
    cz_weights: str = "pupils"     # pupils (partitioned pupils: the four targets) | adults | equal


SPECS: list[Spec] = [
    Spec("headline", "HEADLINE: BA or more (B15003 lines 22-25), adults counted once (unified + elementary records), "
                     "Detroit refiled to 2601103, CZ weights = partitioned pupils"),
    Spec("all_layers", "every LEA takes the record of its own layer (the rule of the earnings targets): adults in an "
                       "elementary district inside a secondary district are counted twice when both sit in one location",
         layers="all"),
    Spec("unified_secondary", "unified + secondary records: the elementary LEAs' adults are dropped in the CZ locations "
                              "that hold a valid secondary record", layers="us"),
    Spec("unified_only", "unified districts only: elementary and secondary LEAs leave the adult sums", layers="u"),
    Spec("ba_only", "bachelor's degree as the HIGHEST degree (line 22 = the old `ba_holders`, over `pop_25plus`), "
                    "adults counted once", measure="ba_only"),
    Spec("nyc_dropped", "New York City's adults leave the sums (its pupils stay in the CZ weights)", nyc=False),
    Spec("drop_ny_cz", "the New York CZ leaves the CZ mean", ny_cz=False),
    Spec("no_detroit_crosswalk", "Detroit's record is not refiled (Detroit has no adults)", detroit=False),
    Spec("adult_weights", "CZ weights = adults 25+ with a valid record in the CZ's partitioned districts",
         cz_weights="adults"),
    Spec("equal_weights", "every CZ weighted equally", cz_weights="equal"),
]
SPEC_BY_KEY = {s.key: s for s in SPECS}


def gaps_b1() -> list["de.Gap"]:
    """`data_estimate.Gap`s of the B1 tables: ratios of location totals (BA+ over adults 25+), the CZ gap in levels
    (suburb minus city) and in logs, and the location totals."""
    return [
        de.Gap("b1_share", "BA-or-more share of adults 25+, suburb minus city", "ba_share_gap",
               "ratio", num="m_ba", den="m_adults", contrast="diff"),
        de.Gap("b1_logratio", "BA-or-more share of adults 25+, log(suburb / city)", "log ba_share ratio",
               "ratio", num="m_ba", den="m_adults", contrast="log"),
        de.Gap("b1_ba", "adults 25+ with the degree level, location total", "-", "total", num="m_ba"),
        de.Gap("b1_adults", "adults 25+, location total", "-", "total", num="m_adults"),
    ]


def add_b15003(A: pd.DataFrame, c: pd.DataFrame) -> pd.DataFrame:
    """ACS columns `b_layer`, `b_adults`, `b_ba` (lines 22-25), `b_ba_only` and `b_name` onto the audited panel by
    `leaid`, after the NYC collapse (the collapsed row carries the Chancellor's id 3620580, under which the ACS files
    New York City).  Detroit's record is refiled from 2612000 to 2601103."""
    n0 = len(A)
    src, dst = DETROIT
    cc = c.copy()
    ids = set(A.leaid)
    if dst in ids and src not in ids:
        cc["leaid"] = cc["leaid"].where(cc["leaid"] != src, dst)
    if cc.leaid.duplicated().any():
        raise SystemExit("the refiled B15003 records share a leaid")
    cc = cc.rename(columns={"name": "b_name", "layer": "b_layer", "adults": "b_adults", "ba_plus": "b_ba",
                            "ba_only": "b_ba_only"})
    X = A.merge(cc, on="leaid", how="left", validate="many_to_one")
    if len(X) != n0:
        raise SystemExit("merge of the B15003 records changed the number of districts")
    X.attrs.update(A.attrs)
    return X


def spec_panel(part: pd.DataFrame, spec: Spec) -> pd.DataFrame:
    """The partitioned districts with `m_adults` and `m_ba` (the columns the Gaps sum) set for a specification: NaN
    where the district's record is invalid, absent, or excluded by the specification."""
    p = part.copy()
    ba = p["b_ba_only" if spec.measure == "ba_only" else "b_ba"]
    ad = p["b_adults"]
    layer = p["b_layer"]
    ok = ad.notna() & ba.notna()
    if spec.layers == "ue":
        keep = layer.isin(["unified", "elementary"])
    elif spec.layers == "all":
        keep = layer.notna()
    elif spec.layers == "u":
        keep = layer.eq("unified")
    elif spec.layers == "us":
        sec = layer.eq("secondary") & ok
        cells = set(zip(p.czone[sec].to_numpy(), p.loc_id[sec].to_numpy()))
        in_cell = np.fromiter(((cz, lo) in cells for cz, lo in zip(p.czone.to_numpy(), p.loc_id.to_numpy())),
                              dtype=bool, count=len(p))
        keep = layer.eq("unified") | sec | (layer.eq("elementary") & ~in_cell)
    else:
        raise ValueError(spec.layers)
    if not spec.nyc:
        keep = keep & p.leaid.ne(NYC_LEAID)
    if not spec.detroit:
        keep = keep & p.leaid.ne(DETROIT[1])
    use = ok & keep
    p["m_adults"], p["m_ba"] = ad.where(use), ba.where(use)
    return p


def cz_table(part: pd.DataFrame, gaps: list["de.Gap"]) -> pd.DataFrame:
    """`czone_gaps` for `gaps` on an already partitioned panel (the swap `district_bootstrap.cz_gap_table` does)."""
    saved = de.GAPS
    try:
        de.GAPS = list(gaps)
        t = de.czone_gaps(part)
    finally:
        de.GAPS = saved
    return t[t.year == YEAR].reset_index(drop=True)


def spec_table(part: pd.DataFrame, g: pd.DataFrame, spec: Spec) -> pd.DataFrame:
    """The B1 columns of one specification on the rows (CZs) of the audited table `g`: czone, the audited CZ weight
    `pupils`, and `b1_*` / `gap_b1_*` (location values and gaps)."""
    t = cz_table(spec_panel(part, spec), gaps_b1())
    if not t.czone.equals(g.czone):
        raise SystemExit("the B1 CZ table does not have the CZs of the audited locale table")
    if not np.allclose(t.pupils.to_numpy(float), g.pupils.to_numpy(float)):
        raise SystemExit("the B1 CZ weights differ from the audited ones")
    gz = g[["czone", "pupils"]].copy()
    for c in t.columns:
        if c.startswith(("b1_", "gap_b1_")):
            gz[c] = t[c].to_numpy(float)
    if not spec.ny_cz:
        ny = part.czone[part.leaid == NYC_LEAID]
        if len(ny):
            cols = [c for c in gz.columns if c.startswith(("b1_", "gap_b1_"))]
            gz.loc[gz.czone == float(ny.iloc[0]), cols] = np.nan
    return gz


def cz_weights(gz: pd.DataFrame, spec: Spec) -> np.ndarray:
    if spec.cz_weights == "pupils":
        return gz.pupils.to_numpy(float)
    if spec.cz_weights == "adults":
        return (gz.b1_adults_1 + gz.b1_adults_2).to_numpy(float)
    if spec.cz_weights == "equal":
        return np.ones(len(gz))
    raise ValueError(spec.cz_weights)


# --------------------------------------------------------------------------
# Statistics (every value they return goes to the sealed file and nowhere else)
# --------------------------------------------------------------------------

def _nan_summary() -> dict[str, Any]:
    return dict(value=float("nan"), se=float("nan"), ci95_lo=float("nan"), ci95_hi=float("nan"), n_cz=0)


def scale_summary(s: dict[str, Any], k: float) -> dict[str, Any]:
    return dict(value=s["value"] * k, se=s["se"] * k, ci95_lo=s["ci95_lo"] * k, ci95_hi=s["ci95_hi"] * k,
                n_cz=s["n_cz"])


def wmean_stat(x: np.ndarray, w: np.ndarray, counts: np.ndarray) -> tuple[dict[str, Any], np.ndarray]:
    """Weighted mean over the CZs where `x` is observed: summary (value, bootstrap SE, 95% percentile interval, n)
    and the B draws, on the multiplicity matrix `counts` (the CZs with a NaN statistic drop out)."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    if not (np.isfinite(x) & np.isfinite(w) & (w > 0)).any():
        return _nan_summary(), np.full(counts.shape[0], np.nan)
    s = db.boot_wmean(x, w, counts)
    return ds.summarise(s), s["draws"]


def pooled_stats(gt: pd.DataFrame, counts: np.ndarray) -> tuple[dict[str, dict[str, Any]], np.ndarray]:
    """Pooled over the CZs that have the gap: sum of BA+ over sum of adults in each location (adult-weighted), and
    the pooled suburb minus city gap, with the same bootstrap draws."""
    have = np.isfinite(gt["gap_b1_share"].to_numpy(float))
    if not have.any():
        return ({k: _nan_summary() for k in ("pooled_gap", "pooled_city", "pooled_suburb")},
                np.full(counts.shape[0], np.nan))

    def one(ba: np.ndarray, ad: np.ndarray) -> tuple[float, np.ndarray]:
        n_, d_ = np.where(have, ba, 0.0), np.where(have, ad, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return float(n_.sum() / d_.sum()), (counts @ n_) / (counts @ d_)

    col = lambda c: gt[c].to_numpy(float)   # noqa: E731
    e1, d1 = one(col("b1_ba_1"), col("b1_adults_1"))
    e2, d2 = one(col("b1_ba_2"), col("b1_adults_2"))
    n = int(have.sum())
    out = {"pooled_city": ds.summarise(dict(est=e1, draws=d1, n=n)),
           "pooled_suburb": ds.summarise(dict(est=e2, draws=d2, n=n)),
           "pooled_gap": ds.summarise(dict(est=e2 - e1, draws=d2 - d1, n=n))}
    return out, d2 - d1


def sample_stats(gt: pd.DataFrame, w: np.ndarray, counts: np.ndarray
                 ) -> tuple[dict[str, dict[str, Any]], dict[str, np.ndarray]]:
    """All B1 statistics of a set of CZs (rows of a `spec_table`) with weights `w` and the multiplicity matrix
    `counts` (B x n): the level gap (suburb minus city; B1), the same in percentage points, the log ratio, the
    CZ-weighted city and suburb shares over the CZs that have the gap, and the pooled versions.  Also the draws of
    the statistics in DRAW_KEYS (for correlations)."""
    gap = gt["gap_b1_share"].to_numpy(float)
    have = np.isfinite(gap)
    stats: dict[str, dict[str, Any]] = {}
    draws: dict[str, np.ndarray] = {}

    def put(key: str, x: np.ndarray) -> None:
        stats[key], d = wmean_stat(x, w, counts)
        if key in DRAW_KEYS:
            draws[key] = d

    put("ba_share_gap", gap)
    stats["ba_share_gap_pp"] = scale_summary(stats["ba_share_gap"], 100.0)
    put("ba_log_ratio", gt["gap_b1_logratio"].to_numpy(float))
    put("ba_share_city", np.where(have, gt["b1_share_1"].to_numpy(float), np.nan))
    put("ba_share_suburb", np.where(have, gt["b1_share_2"].to_numpy(float), np.nan))
    pooled, d_pooled = pooled_stats(gt, counts)
    stats.update(pooled)
    draws["pooled_gap"] = d_pooled
    return {k: stats[k] for k in STAT_KEYS}, draws


def corr_block(draws: dict[str, np.ndarray], tdraws: pd.DataFrame) -> dict[str, Any]:
    """Bootstrap correlations among the draws of the B1 statistics and of the four district targets (same
    multiplicity matrix)."""
    D = pd.concat([tdraws[TARGETS].reset_index(drop=True), pd.DataFrame(draws)], axis=1)
    C = D.corr()
    return dict(keys=list(D.columns), matrix=C.to_numpy(float).round(6).tolist(),
                ba_share_gap_with={t: float(C.loc["ba_share_gap", t]) for t in TARGETS},
                ba_share_gap_with_pooled=float(C.loc["ba_share_gap", "pooled_gap"]),
                ba_share_gap_with_logratio=float(C.loc["ba_share_gap", "ba_log_ratio"]))


def sample_masks(gz: pd.DataFrame, ho: dict[str, Any]) -> dict[str, np.ndarray]:
    cz = gz.czone.to_numpy(float)
    return {"full": np.ones(len(gz), dtype=bool), "A": np.isin(cz, np.asarray(ho["A"], float)),
            "B": np.isin(cz, np.asarray(ho["B"], float))}


def run_all(R: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Every statistic of every specification on the full sample and the two hold-out halves (in memory only)."""
    res: dict[str, dict[str, Any]] = {}
    for spec in SPECS:
        gz = R["tables"][spec.key]
        w = cz_weights(gz, spec)
        res[spec.key] = {}
        for s in SAMPLES:
            m = R["masks"][s]
            gt = gz[m].reset_index(drop=True)
            stats, dr = sample_stats(gt, w[m], R["count_mats"][s])
            res[spec.key][s] = dict(n_cz=int(m.sum()), n_cz_gap=int(stats["ba_share_gap"]["n_cz"]), stats=stats,
                                    draws=dr, corr=corr_block(dr, R["std_draws"][s]))
    return res


# --------------------------------------------------------------------------
# Checks that involve no sealed value
# --------------------------------------------------------------------------

def ok_str(flag: bool) -> str:
    return "PASS" if flag else "FAIL"


def reproduction_check(std: dict[str, Any], full_run: bool) -> list[dict[str, Any]]:
    rows = []
    for k in ds.TARGET_KEYS:
        s = std["stats"][k]
        row = dict(moment=k, value=s["value"], se=s["se"], n_cz=s["n_cz"], reported_value=ds.REPORTED[k][0],
                   reported_se=ds.REPORTED[k][1], ok_value=bool(round(s["value"], 4) == ds.REPORTED[k][0]),
                   ok_se=bool(round(s["se"], 4) == ds.REPORTED[k][1]) if full_run else None)
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
        return dict(compared=False, identical=None)
    st = json.loads(path.read_text())["holdout"]
    same = (sorted(int(c) for c in st["A_cz"]) == sorted(int(c) for c in ho["A"])
            and sorted(int(c) for c in st["B_cz"]) == sorted(int(c) for c in ho["B"]) and st["seed"] == ho["seed"])
    if not same:
        raise SystemExit("the hold-out split differs from the one stored in district_settings.json; stopping")
    LOG.info("hold-out split identical to district_settings.json (seed %d, %d / %d CZs)", ho["seed"],
             len(ho["A"]), len(ho["B"]))
    return dict(compared=True, identical=True, seed=int(ho["seed"]), n_A=len(ho["A"]), n_B=len(ho["B"]))


def draws_identity(g: pd.DataFrame, std: dict[str, Any], counts: np.ndarray, indir: Path, full_run: bool
                   ) -> dict[str, Any]:
    """The multiplicities used here are those of district_settings.py: the draws of each of the four targets
    recomputed from the CZ table with `counts` equal the draws of `moment_set`, and (full run) the stored
    `district_settings_draws_locale.csv`."""
    w = g.pupils.to_numpy(float)
    mine = {k: float(np.nanmax(np.abs(db.boot_wmean(g[ds.BY_KEY[k].col].to_numpy(float), w, counts)["draws"]
                                      - std["draws"][k].to_numpy(float)))) for k in TARGETS}
    out: dict[str, Any] = dict(recomputed_max_abs_diff=mine, csv_max_abs_diff=None)
    csv = indir / "estimates" / "settings" / "district_settings_draws_locale.csv"
    if full_run and csv.exists():
        c = pd.read_csv(csv)
        out["csv_max_abs_diff"] = {k: float(np.nanmax(np.abs(c[k].to_numpy(float) - std["draws"][k].to_numpy(float))))
                                   for k in TARGETS}
        out["csv_rows"] = int(len(c))
    ok = max(mine.values()) <= 1e-12 and (out["csv_max_abs_diff"] is None or max(out["csv_max_abs_diff"].values()) <= 1e-9)
    out["ok"] = bool(ok)
    return out


def machinery_check(part: pd.DataFrame, g: pd.DataFrame, std: dict[str, Any], counts: np.ndarray) -> dict[str, Any]:
    """The two code paths B1 uses, applied to PUBLIC quantities: `wmean_stat` on the income gap reproduces the
    income target (value and SE), and `cz_table` with a `Gap` reproduces the enrollment gap column of the audited
    table."""
    s, _ = wmean_stat(g[ds.BY_KEY["gap_income"].col].to_numpy(float), g.pupils.to_numpy(float), counts)
    t = cz_table(part, [de.Gap("pupils", "enrollment", "M_2/M_1", "total", num="pupils")])
    d_gap = float(np.nanmax(np.abs(t.gap_pupils.to_numpy(float) - g.gap_pupils.to_numpy(float))))
    d_val = abs(s["value"] - std["stats"]["gap_income"]["value"])
    d_se = abs(s["se"] - std["stats"]["gap_income"]["se"])
    return dict(max_abs_diff_wmean=float(max(d_val, d_se)), max_abs_diff_cz_table=d_gap,
                ok=bool(max(d_val, d_se) <= 1e-12 and d_gap <= 1e-12))


# --------------------------------------------------------------------------
# Coverage (public counts)
# --------------------------------------------------------------------------

def cz_pupil_coverage(part: pd.DataFrame) -> pd.DataFrame:
    """Share of each CZ location's partitioned pupils whose district carries a valid headline ACS record (columns 1 =
    city, 2 = suburb; index czone).  A public count of pupils: it says how much of the location the adult share
    covers, not what the share is.  (Pupils of a secondary LEA count as unmatched although the same children's
    adults are counted through the elementary districts, so this understates the coverage of adults.)"""
    hp = spec_panel(part, SPEC_BY_KEY["headline"])
    hp["matched"] = hp.pupils.where(hp.m_adults.notna(), 0.0)
    t = hp.groupby(["czone", "loc_id"]).agg(pu=("pupils", "sum"), mu=("matched", "sum"))
    return (t.mu / t.pu).unstack("loc_id").sort_index()


def coverage(part: pd.DataFrame, g: pd.DataFrame, tables: dict[str, pd.DataFrame], ho: dict[str, Any]) -> dict[str, Any]:
    """Sizes of the partition and of the matched ACS sample under every specification.  Counts, pupils and adult
    totals only: nothing here involves a degree count."""
    pu = float(part.pupils.sum())
    layer = part.b_layer.fillna("none")
    out: dict[str, Any] = dict(
        n_cz=int(len(g)), n_districts=int(len(part)), pupils=pu,
        n_by_location={int(k): int(v) for k, v in part.groupby("loc_id").size().items()},
        pupils_by_location={int(k): float(v) for k, v in part.groupby("loc_id").pupils.sum().items()},
        n_by_layer={k: int(v) for k, v in layer.value_counts().items()},
        pupil_share_by_layer={k: float(part.pupils[layer == k].sum() / pu) for k in layer.unique()},
        n_A=int(len(ho["A"])), n_B=int(len(ho["B"])))
    w = g.pupils.to_numpy(float)
    wn = w / w.sum()
    out["n_eff_cz"] = float(1.0 / (wn ** 2).sum())
    out["top5_weight_share"] = float(np.sort(wn)[::-1][:5].sum())
    specs: dict[str, Any] = {}
    for s in SPECS:
        p = spec_panel(part, s)
        v = p.m_adults.notna()
        gz = tables[s.key]
        have = np.isfinite(gz.gap_b1_share.to_numpy(float))
        specs[s.key] = dict(
            n_valid=int(v.sum()), share_districts_valid=float(v.mean()), share_pupils_valid=float(p.pupils[v].sum() / pu),
            share_pupils_valid_by_location={int(lo): float(p.pupils[v & (p.loc_id == lo)].sum() / p.pupils[p.loc_id == lo].sum())
                                            for lo in (1, 2)},
            n_cz_gap=int(have.sum()), weight_share_cz_gap=float(gz.pupils.to_numpy(float)[have].sum() / gz.pupils.sum()),
            n_cz_gap_A=int((have & np.isin(gz.czone.to_numpy(float), np.asarray(ho["A"], float))).sum()),
            n_cz_gap_B=int((have & np.isin(gz.czone.to_numpy(float), np.asarray(ho["B"], float))).sum()))
    out["specs"] = specs
    cc = cz_pupil_coverage(part)
    flat = cc.to_numpy(float).ravel()
    wcz = g.set_index("czone").pupils
    below = lambda th: cc.index[(cc < th).any(axis=1)]                      # noqa: E731
    out["cz_location_coverage"] = dict(
        n=int(len(flat)), min=float(flat.min()), p5=float(np.quantile(flat, 0.05)), median=float(np.median(flat)),
        n_below_50=int((flat < 0.5).sum()), n_below_80=int((flat < 0.8).sum()),
        weight_share_cz_below_80=float(wcz.reindex(below(0.8)).sum() / wcz.sum()))
    head, allp = spec_panel(part, SPEC_BY_KEY["headline"]), spec_panel(part, SPEC_BY_KEY["all_layers"])
    low = {(float(c), int(l)) for c in cc.index for l in cc.columns if cc.loc[c, l] < 0.8}
    in_low = np.fromiter(((float(c), int(l)) in low for c, l in zip(head.czone.to_numpy(), head.loc_id.to_numpy())),
                         dtype=bool, count=len(head))
    um = head.m_adults.isna().to_numpy() & in_low
    sec_um = um & head.b_layer.eq("secondary").to_numpy() & head.b_adults.notna().to_numpy()
    out["cz_location_coverage"]["secondary_share_of_unmatched_pupils_below_80"] = (
        float(head.pupils.to_numpy(float)[sec_um].sum() / head.pupils.to_numpy(float)[um].sum()) if um.any() else float("nan"))
    out["adults"] = dict(
        headline=float(head.m_adults.sum()), all_layers=float(allp.m_adults.sum()),
        overlap_factor=float(allp.m_adults.sum() / head.m_adults.sum()),
        overlap_factor_by_location={int(lo): float(allp.m_adults[allp.loc_id == lo].sum() / head.m_adults[head.loc_id == lo].sum())
                                    for lo in (1, 2)})
    sec = part.b_layer.eq("secondary") & part.b_adults.notna()
    out["secondary_excluded"] = dict(n=int(sec.sum()), pupils=float(part.pupils[sec].sum()),
                                     pupil_share=float(part.pupils[sec].sum() / pu))
    # districts that leave the headline sample, and why
    bad = head[head.m_adults.isna()].sort_values("pupils", ascending=False)
    reason = np.where(bad.b_layer.isna(), "no ACS record",
                      np.where(bad.b_layer.eq("secondary") & bad.b_adults.notna(),
                               "secondary record (adults counted from the elementary records)", "invalid record"))
    out["unmatched"] = dict(
        n=int(len(bad)), pupils=float(bad.pupils.sum()), pupil_share=float(bad.pupils.sum() / pu),
        by_reason={k: dict(n=int((reason == k).sum()), pupils=float(bad.pupils.to_numpy(float)[reason == k].sum()))
                   for k in sorted(set(reason))},
        top=[dict(leaid=str(r.leaid), name=str(r.lea_name), pupils=float(r.pupils), locale=str(r.locale_grp), reason=str(rs))
             for r, rs in zip(bad.head(10).itertuples(), reason[:10])])
    nr = bad[bad.b_layer.isna()]
    out["unmatched"]["top_no_record"] = [dict(leaid=str(r.leaid), name=str(r.lea_name), pupils=float(r.pupils), locale=str(r.locale_grp))
                                         for r in nr.head(10).itertuples()]
    for nm, lid in (("nyc", NYC_LEAID), ("detroit", DETROIT[1])):
        r = part[part.leaid == lid]
        if len(r):
            r = r.iloc[0]
            cz_pu = float(part.pupils[part.czone == r.czone].sum())
            out[nm] = dict(in_partition=True, leaid=lid, location=int(r.loc_id), pupils=float(r.pupils),
                           pupil_share_of_cz=float(r.pupils / cz_pu), pupil_share_of_weights=float(cz_pu / pu),
                           record_layer=None if pd.isna(r.b_layer) else str(r.b_layer),
                           has_valid_record=bool(pd.notna(r.b_adults)),
                           in_headline=bool(pd.notna(head.m_adults[head.leaid == lid].iloc[0])))
        else:
            out[nm] = dict(in_partition=False, leaid=lid)
    return out


# --------------------------------------------------------------------------
# prepare: everything that involves no sealed statistic
# --------------------------------------------------------------------------

def prepare(B: int, seed: int, indir: Path = HERE) -> dict[str, Any]:
    """Panel, partition, hold-out split, multiplicity matrices, the district targets' draws and every public check;
    the B1 tables of every specification (data, not statistics).  Nothing sealed is computed here."""
    _stack()
    t0 = time.time()
    full_run = B == B_DEFAULT and seed == SEED
    if (B_DEFAULT, SEED, HOLDOUT_SEED) != (db.B_DEFAULT, db.SEED, ds.HOLDOUT_SEED):
        raise SystemExit("B_DEFAULT / SEED / HOLDOUT_SEED differ from district_bootstrap.py / district_settings.py")
    checks: list[tuple[str, bool, str]] = []

    def chk(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))
        LOG.info("check %-58s %s%s", name, ok_str(bool(ok)), f"  {detail}" if detail else "")

    # --- ACS records ------------------------------------------------------------------------------
    raw_new = pd.read_parquet(B15003_PATH)
    states = pd.read_parquet(B15003_ST_PATH)
    labels = json.loads(B15003_LABELS.read_text())
    check_labels(labels["labels"])
    chk("Census labels of lines 1 and 22-25", True, "total; bachelor's; master's; professional school; doctorate")
    rec = records_from_raw(raw_new)
    pull = pull_vs_old(raw_new, load_old_file())
    chk("new pull = old file (same records, pop_25plus, ba_holders)",
        pull["n_both"] == pull["n_new"] == pull["n_old"] and pull["n_equal_pop"] == pull["n_both"]
        and pull["n_equal_ba"] == pull["n_both"], f"{pull['n_equal_ba']:,} of {pull['n_both']:,} records")
    chk("B15003 lines 2-25 add up to line 1; no annotation codes",
        pull["n_lines_sum_to_total"] == pull["n_new"] and pull["n_with_negative_code"] == 0 and pull["n_with_missing"] == 0)
    chk("one record per leaid after dropping the 'Remainder' rows", rec.attrs["n_duplicates_dropped"] == 0,
        f"{rec.attrs['n_remainder_dropped']} remainder rows dropped, {len(rec):,} records")
    chk("every record that is not used is an empty district (0 adults 25+)",
        rec.attrs["n_invalid"] == rec.attrs["n_zero_adults"],
        f"{rec.attrs['n_zero_adults']} records with no adults, {rec.attrs['n_invalid']} not used")
    tiling = tiling_check(raw_new, states)
    records = dict(n_records=rec.attrs["n_records"], n_remainder_dropped=rec.attrs["n_remainder_dropped"],
                   n_unique=int(len(rec)), n_invalid=rec.attrs["n_invalid"],
                   n_by_layer={k: int(v) for k, v in rec.layer.value_counts().items()},
                   adults_by_layer={k: float(v) for k, v in rec.groupby("layer").adults.sum().items()})

    # --- panel, partition, the four targets ---------------------------------------------------------
    inp = ds.build_audited_panel(indir)
    A = add_b15003(inp["A"], rec)
    part, g = ds.scheme_tables(A, "locale")
    part = part[part.year == YEAR].reset_index(drop=True)
    std = ds.moment_set(g, part, B, seed)
    repro = reproduction_check(std, full_run)
    chk("the four locale targets reproduce the reported values", True,
        "value" + (" and SE" if full_run else " only (B differs from the reported runs)"))

    # --- hold-out split and multiplicity matrices ------------------------------------------------------
    ho = ds.holdout_split(g)
    ho_check = check_holdout(ho, indir)
    chk("hold-out split = district_settings.json", ho_check.get("identical") is True,
        f"{len(ho['A'])} / {len(ho['B'])} CZs, seed {ho['seed']}")
    masks_g = {"full": np.ones(len(g), dtype=bool), "A": np.isin(g.czone.to_numpy(float), np.asarray(ho["A"], float)),
               "B": np.isin(g.czone.to_numpy(float), np.asarray(ho["B"], float))}
    count_mats, std_draws = {}, {}
    for s in SAMPLES:
        n_s = int(masks_g[s].sum())
        count_mats[s] = db.draw_counts(n_s, B, seed)
        if s == "full":
            std_draws[s] = std["draws"]
        else:
            gH, partH = ds.restrict(g, part, set(ho[s]))
            std_draws[s] = ds.moment_set(gH, partH, B, seed)["draws"]
            if len(gH) != n_s:
                raise SystemExit("hold-out half has a different number of CZs")
    dchk = draws_identity(g, std, count_mats["full"], indir, full_run)
    chk("bootstrap draws = district_settings.py (recomputed; stored csv)", dchk["ok"],
        f"max abs diff {max(dchk['recomputed_max_abs_diff'].values()):.1e}"
        + ("" if dchk["csv_max_abs_diff"] is None else f"; csv {max(dchk['csv_max_abs_diff'].values()):.1e}"))
    mach = machinery_check(part, g, std, count_mats["full"])
    chk("statistic machinery reproduces the income and enrollment targets", mach["ok"],
        f"max abs diff {max(mach['max_abs_diff_wmean'], mach['max_abs_diff_cz_table']):.1e}")

    # --- B1 tables of every specification --------------------------------------------------------------
    tables = {s.key: spec_table(part, g, s) for s in SPECS}
    chk("B1 tables: the 188 CZs and weights of the audited table, every specification",
        all(t.czone.equals(g.czone) for t in tables.values()), f"{len(g)} CZs x {len(SPECS)} specifications")
    cov = coverage(part, g, tables, ho)
    hs = cov["specs"]["headline"]
    chk("headline: valid gap in >= 90% of the CZs and of the CZ weight",
        hs["n_cz_gap"] >= COVERAGE_MIN * len(g) and hs["weight_share_cz_gap"] >= COVERAGE_MIN,
        f"{hs['n_cz_gap']} of {len(g)} CZs, {100 * hs['weight_share_cz_gap']:.1f}% of the weight")
    chk("New York City and Detroit have a record in the headline sample",
        bool(cov["nyc"].get("in_headline")) and bool(cov["detroit"].get("in_headline")))
    for s in SAMPLES:
        if cov["specs"]["headline"].get(f"n_cz_gap_{s}", 1) == 0:
            chk(f"headline: half {s} has CZs with a valid gap", False)
    R = dict(B=B, seed=seed, full_run=full_run, A=A, part=part, g=g, std=std, repro=repro, ho=ho, ho_check=ho_check,
             cz_cov=cz_pupil_coverage(part),
             masks=masks_g, count_mats=count_mats, std_draws=std_draws, draw_checks=dchk, machinery=mach,
             tables=tables, coverage=cov, pull=pull, tiling=tiling, records=records, labels=labels,
             checks=checks, checks_ok=all(c[1] for c in checks), seconds=time.time() - t0,
             hybrid_share=inp["hybrid_share"], geo_unit=A.attrs.get("geo_unit"))
    return R


# --------------------------------------------------------------------------
# The sealed file
# --------------------------------------------------------------------------

DEFINITION = dict(
    object=("B1 (PREREG, experiments-2026-10d): within-CZ gap in the share of adults 25+ with a bachelor's degree or "
            "more, suburb minus city; model counterpart `ba_share_gap` (M3)"),
    ba_or_more="B15003_022E + 023E + 024E + 025E (bachelor's, master's, professional school, doctorate)",
    adults="B15003_001E (population 25 years and over)",
    ba_only="B15003_022E / B15003_001E (= `ba_holders` / `pop_25plus` of raw/acs_school_districts.parquet)",
    ba_holders_note=("PREREG names `ba_holders` / `pop_25plus` of raw/acs_school_districts.parquet as the input of B1. `ba_holders` is "
                     "B15003_022E = adults whose HIGHEST degree is a bachelor's (graduate and professional degrees are lines 23-25), "
                     "so the headline uses lines 22-25 of the same table (pulled again, checked equal to the old file on lines 1 and 22); "
                     "the literal PREREG input is the variant `ba_only`."),
    location_share="s_l = sum of BA+ over the location's districts / sum of adults 25+ (a ratio of location totals)",
    ba_share_gap="CZ-weighted mean (weights = partitioned pupils) of s_suburb - s_city; a fraction",
    ba_share_gap_pp="100 x ba_share_gap",
    ba_log_ratio="CZ-weighted mean of log(s_suburb / s_city)",
    ba_share_city="CZ-weighted mean of s_city over the CZs that have the gap",
    ba_share_suburb="CZ-weighted mean of s_suburb over the CZs that have the gap",
    pooled_gap=("sum of BA+ over sum of adults, all suburb districts of the CZs that have the gap, minus the same for "
                "the city districts (adult-weighted, not differenced within a CZ)"),
    pooled_city="the city term of pooled_gap", pooled_suburb="the suburb term of pooled_gap",
    partition=("district_settings locale partition: city = NCES locale 11-13, suburb = 21-23, districts of >= 500 "
               "pupils, each location >= 10% of the CZ's partitioned pupils"),
    bootstrap=("whole-CZ resampling with replacement, B draws, multiplicity matrix district_bootstrap.draw_counts(n_cz, B, "
               "seed); se = std (ddof 1) of the draws, ci95 = 2.5 / 97.5 percentiles; CZs with a NaN statistic drop out"),
    halves="A and B: the pre-registered split of district_settings.holdout_split (seed 20261002); each half bootstrapped on its own n CZs",
    layers=dict(ue="unified + elementary records (headline: adults counted once)", all="every record (earnings rule)",
                us="unified + secondary; elementary dropped in the CZ locations that hold a valid secondary record",
                u="unified records only"),
    units="shares and gaps are fractions in [0, 1] (ba_share_gap_pp in percentage points); ba_log_ratio in logs")


def code_hashes() -> dict[str, str]:
    return {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("ba_share_holdout.py", "ba_share_holdout_xcheck.py") if (HERE / f).exists()}


def input_hashes() -> dict[str, str]:
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (B15003_PATH, B15003_ST_PATH, OLD_FILE)}


def cz_rows(R: dict[str, Any]) -> list[dict[str, Any]]:
    """Per-CZ table of the headline (sealed)."""
    gz = R["tables"]["headline"]
    half = np.where(R["masks"]["A"], "A", "B")
    cc = R["cz_cov"].reindex(gz.czone.to_numpy(float))
    rows = []
    for i, r in enumerate(gz.itertuples()):
        rows.append(dict(czone=int(r.czone), half=str(half[i]), pupils=float(r.pupils), share_city=r.b1_share_1,
                         share_suburb=r.b1_share_2, gap=r.gap_b1_share, log_ratio=r.gap_b1_logratio,
                         ba_city=r.b1_ba_1, ba_suburb=r.b1_ba_2, adults_city=r.b1_adults_1, adults_suburb=r.b1_adults_2,
                         pupil_cov_city=float(cc[1].iloc[i]), pupil_cov_suburb=float(cc[2].iloc[i])))
    return rows


def build_payload(R: dict[str, Any], res: dict[str, dict[str, Any]]) -> dict[str, Any]:
    h = res["headline"]
    cz = R["tables"]["headline"].czone.to_numpy(float)
    payload = dict(
        sealed=True,
        note=("SEALED hold-out B1 of the 2026-10-15 round (ACS 2014-18 school districts, BA-or-more share of adults 25+, "
              "suburb minus city, audited locale partition). Opened only after the specification is chosen (PREREG A7): "
              "no log, report or summary reads or prints these numbers."),
        meta=dict(created=dt.datetime.now().isoformat(timespec="seconds"), script="data/spatial/ba_share_holdout.py",
                  B=R["B"], bootstrap_seed=R["seed"], holdout_seed=int(R["ho"]["seed"]), year=YEAR,
                  vintage=f"ACS 5-year {YEAR - 4}-{YEAR}", table="B15003", n_cz=int(len(R["g"])), full_run=R["full_run"],
                  geo_unit=R["geo_unit"], hybrid_share=R["hybrid_share"], code_sha256=code_hashes(),
                  inputs_sha256=input_hashes()),
        definition=DEFINITION,
        specs={s.key: dict(asdict(s), label=s.label) for s in SPECS},
        headline=dict(key="ba_share_gap", spec="headline",
                      **{s: h[s]["stats"]["ba_share_gap"] for s in SAMPLES},
                      pp={s: h[s]["stats"]["ba_share_gap_pp"] for s in SAMPLES}),
        samples={s: dict(n_cz=h[s]["n_cz"], n_cz_with_gap=h[s]["n_cz_gap"],
                         cz_ids=[int(c) for c in cz[R["masks"][s]]], stats=h[s]["stats"], corr=h[s]["corr"])
                 for s in SAMPLES},
        variants={sp.key: dict(label=sp.label, spec=asdict(sp),
                               samples={s: dict(n_cz=res[sp.key][s]["n_cz"], n_cz_with_gap=res[sp.key][s]["n_cz_gap"],
                                                stats=res[sp.key][s]["stats"],
                                                corr_with_targets=res[sp.key][s]["corr"]["ba_share_gap_with"])
                                        for s in SAMPLES})
                  for sp in SPECS[1:]},
        draws_full={k: [round(float(x), 10) for x in h["full"]["draws"][k]] for k in DRAW_KEYS},
        cz_table=cz_rows(R))
    return ds.clean(payload)


def n_statistics(payload: dict[str, Any]) -> int:
    """Number of (sample, statistic) cells in the sealed file: headline + variants."""
    n = sum(len(v["stats"]) for v in payload["samples"].values())
    return n + sum(len(s["stats"]) for v in payload["variants"].values() for s in v["samples"].values())


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_sealed(path: Path, payload: dict[str, Any]) -> None:
    """json.dump to a temporary file, then rename: the only place where a statistic is written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=1, allow_nan=False)
    os.replace(tmp, path)


def _all_stats(payload: dict[str, Any]) -> list[dict[str, Any]]:
    out = [s for v in payload["samples"].values() for s in v["stats"].values()]
    out += [s for v in payload["variants"].values() for sm in v["samples"].values() for s in sm["stats"].values()]
    return out


def verify_seal(path: Path) -> int:
    """Structure of the sealed file, without printing any value: keys, number of statistics, whether they are
    finite and consistent, B, seed, sizes of the samples, SHA-256."""
    d = json.loads(path.read_text())
    S = _all_stats(d)
    fin = lambda x: x is not None and np.isfinite(x)   # noqa: E731
    n_null = sum(1 for s in S if not fin(s["value"]))
    h = d["samples"]
    full = h["full"]["stats"]
    checks = {
        "every statistic of the headline and variants has a value": n_null == 0,
        "every statistic has a positive SE and a 95% interval":
            all(fin(s["se"]) and s["se"] > 0 and fin(s["ci95_lo"]) and fin(s["ci95_hi"]) for s in S),
        "intervals contain the point estimates": all(s["ci95_lo"] <= s["value"] <= s["ci95_hi"] for s in S),
        "shares in [0, 1]": all(0.0 <= v["stats"][k]["value"] <= 1.0 for v in h.values()
                                for k in ("ba_share_city", "ba_share_suburb", "pooled_city", "pooled_suburb")),
        "pp = 100 x fraction": all(abs(v["stats"]["ba_share_gap_pp"]["value"] - 100 * v["stats"]["ba_share_gap"]["value"]) < 1e-9
                                   for v in h.values()),
        "suburb share - city share = gap (same CZs, same weights)":
            all(abs(v["stats"]["ba_share_suburb"]["value"] - v["stats"]["ba_share_city"]["value"]
                    - v["stats"]["ba_share_gap"]["value"]) < 1e-10 for v in h.values()),
        "pooled suburb - city = pooled gap": all(abs(v["stats"]["pooled_suburb"]["value"] - v["stats"]["pooled_city"]["value"]
                                                    - v["stats"]["pooled_gap"]["value"]) < 1e-10 for v in h.values()),
        "halves partition the CZs": (h["A"]["n_cz"] + h["B"]["n_cz"] == h["full"]["n_cz"]
                                     and not set(h["A"]["cz_ids"]) & set(h["B"]["cz_ids"])
                                     and sorted(h["A"]["cz_ids"] + h["B"]["cz_ids"]) == sorted(h["full"]["cz_ids"])),
        "correlations finite and in [-1, 1]": all(fin(x) and -1.0 <= x <= 1.0 for v in h.values()
                                                  for row in v["corr"]["matrix"] for x in row),
        "four target correlations present in every sample": all(set(v["corr"]["ba_share_gap_with"]) == set(TARGETS)
                                                                for v in h.values()),
        "draws: B finite draws per statistic, the mean of the draws is within 5 SE of the estimate":
            all(len(d["draws_full"][k]) == d["meta"]["B"]
                and np.isfinite(np.asarray(d["draws_full"][k], dtype=float)).all()
                and abs(np.mean(np.asarray(d["draws_full"][k], dtype=float)) - full[DRAW_STAT[k]]["value"])
                <= 5 * full[DRAW_STAT[k]]["se"] for k in DRAW_KEYS),
        "draws: the SE of the draws is the reported SE": all(
            abs(np.std(np.asarray(d["draws_full"][k], dtype=float), ddof=1) - full[DRAW_STAT[k]]["se"])
            <= 1e-6 * max(1.0, abs(full[DRAW_STAT[k]]["se"])) + 1e-9 for k in DRAW_KEYS),
        "cz table: one row per CZ, finite gap": len(d["cz_table"]) == h["full"]["n_cz"]
                                                and sum(1 for r in d["cz_table"] if fin(r["gap"])) == h["full"]["n_cz_with_gap"],
    }
    print(f"sealed file {path.name} ({path.stat().st_size:,} bytes): top-level keys {list(d)}")
    print(f"  {len(S)} statistics (sample x statistic cells), {len(d['variants'])} variants + headline, "
          f"B = {d['meta']['B']}, bootstrap seed {d['meta']['bootstrap_seed']}, hold-out seed {d['meta']['holdout_seed']}")
    print(f"  CZs: full {h['full']['n_cz']} ({h['full']['n_cz_with_gap']} with a gap), A {h['A']['n_cz']} "
          f"({h['A']['n_cz_with_gap']}), B {h['B']['n_cz']} ({h['B']['n_cz_with_gap']}); {len(d['cz_table'])} rows in cz_table; "
          f"{len(full)} statistics per sample")
    ok_all = True
    for name, ok in checks.items():
        ok_all &= bool(ok)
        print(f"  {name}: {ok_str(bool(ok))}")
    print(f"  sha256 {sha256_file(path)}")
    return 0 if ok_all else 1


def _floats(o: Any):
    if isinstance(o, dict):
        for v in o.values():
            yield from _floats(v)
    elif isinstance(o, (list, tuple)):
        for v in o:
            yield from _floats(v)
    elif isinstance(o, float):
        yield o


DECIMALS = re.compile(r"(?<![\d.])\d+\.\d+(?!\d)")        # a maximal decimal number in a text
FINE = dict(dec_frac=(5, 6), dec_pp=(4,))                    # high-precision renderings: a hit is never a coincidence
COARSE = dict(dec_frac=(4,), dec_pp=(2, 3))                  # the renderings a careless print would use
N_DECOYS = 200


def decimal_strings(text: str) -> set[str]:
    """Every maximal decimal number of a text, as a string (not preceded or followed by a digit or a point)."""
    return set(DECIMALS.findall(text))


def sealed_floats(payload: dict[str, Any]) -> np.ndarray:
    """Every float of the headline, samples and variants blocks with 0 < |x| < 100."""
    v = [x for blk in (payload["headline"], payload["samples"], payload["variants"]) for x in _floats(blk)]
    return np.array([x for x in v if x != 0.0 and abs(x) < 100.0], dtype=float)


def render(vals: np.ndarray, dec_frac: tuple, dec_pp: tuple) -> set[str]:
    """Strings under which sealed numbers could appear: |x| at `dec_frac` decimals, and 100 |x| at `dec_pp` decimals."""
    out: set[str] = set()
    for v in vals:
        for d in dec_frac:
            out.add(f"{abs(v):.{d}f}")
        for d in dec_pp:
            out.add(f"{abs(100.0 * v):.{d}f}")
    return out


def leak_tokens(payload: dict[str, Any], **kind: Any) -> set[str]:
    return render(sealed_floats(payload), **(kind or FINE))


def leak_scan(text: str, payload: dict[str, Any], **kind: Any) -> int:
    """Number of sealed-number renderings found in `text`.  With the default (FINE: fractions at 5 and 6 decimals,
    percentage points at 4) a hit cannot be a coincidence with a public number -- the public numbers of the methods
    file and the logs have at most four decimals as fractions -- and the scan must find none."""
    return len(leak_tokens(payload, **kind) & decimal_strings(text))


def coarse_scan(text: str, payload: dict[str, Any], n_decoys: int = N_DECOYS, seed: int = 20261015) -> dict[str, Any]:
    """The renderings of a careless print (fractions at 4 decimals, percentage points at 2 and 3) do collide with the
    public numbers of a text by chance (shares of pupils, tiling ratios, the four reported targets).  Their count is
    therefore compared with the count for `n_decoys` random decoys of the sealed numbers (each multiplied by an
    independent factor in [0.5, 2]): the real count must not exceed the largest decoy count."""
    vals = sealed_floats(payload)
    nums = decimal_strings(text)
    real = len(render(vals, **COARSE) & nums)
    rng = np.random.default_rng(seed)
    dec = [len(render(vals * np.exp(rng.uniform(-0.7, 0.7, size=len(vals))), **COARSE) & nums) for _ in range(n_decoys)]
    return dict(real=real, decoy_mean=float(np.mean(dec)), decoy_max=int(max(dec)), ok=bool(real <= max(dec)))


def scan_report(name: str, text: str, payload: dict[str, Any]) -> bool:
    """One line per scanned text: counts and PASS/FAIL only."""
    fine = leak_scan(text, payload)
    co = coarse_scan(text, payload)
    ok = fine == 0 and co["ok"]
    print(f"leak scan of {name}: {len(leak_tokens(payload)):,} high-precision renderings of the sealed numbers, {fine} found; "
          f"coarse renderings {co['real']} found vs {co['decoy_mean']:.1f} on average (max {co['decoy_max']}) in {N_DECOYS} random decoys  "
          f"[{ok_str(ok)}]")
    return ok


# --------------------------------------------------------------------------
# The methods note (public counts only)
# --------------------------------------------------------------------------

def fn(x: Any) -> str:
    return f"{int(round(float(x))):,}"


def fp(x: float, d: int = 1) -> str:
    return f"{100 * float(x):.{d}f}%"


def write_methods(path: Path, R: dict[str, Any], k_stats: int, sha: str | None = None) -> None:
    """The methods note: definitions, samples, coverage and PASS/FAIL checks.  Reads only `prepare` output (public
    counts, no statistic of B1), never the results of `run_all`."""
    cov, pull, til, rec = R["coverage"], R["pull"], R["tiling"], R["records"]
    hs = cov["specs"]["headline"]
    L: list[str] = []
    a = L.append
    gen = dt.datetime.now().isoformat(timespec="seconds")
    a("# B1 hold-out (BA-or-more share gap, suburb minus city): methods (no values)\n")
    a(f"Generated {gen} by `ba_share_holdout.py`. **This file contains definitions, samples, coverage and PASS/FAIL checks only; it reports no "
      "value of B1 or of any variant.** The statistics are in `estimates/settings/sealed/ba_share_holdout_B1.json` (git-ignored), written "
      "with `json.dump` and neither printed nor logged. It is to be opened by a person only after the model specification of the 2026-10-15 round "
      "is chosen (`julia/spatial_model/calibration/runs/experiments-2026-10d/PREREG.md`, A7); no specification is re-selected on it.\n")
    a("## 1. What is sealed\n")
    a("**B1**, pre-registered under \"Sealed holdout\" in the 2026-10-15 PREREG: the within-CZ gap in the share of adults 25+ with a bachelor's "
      "degree or more, suburb minus city, ACS 5-year 2014-18 by school district, on the audited locale partition and CZ weights of "
      "`district_settings.py` (the CZs of the four district targets). The model counterpart is the BA share of adults by residence "
      "(`ba_share_by_residence` / `ba_share_gap`, M3 only). The sealed file holds, for the full sample and for the two pre-registered "
      "hold-out halves, the level gap (fraction and percentage points), the log ratio, the CZ-weighted city and suburb shares, the version pooled "
      "over CZs, whole-CZ bootstrap SEs and 95% intervals, the bootstrap correlations with the four district targets "
      "(`gap_income`, `gap_pupils`, `gap_salary`, `gap_teachers_pp`), the bootstrap draws of the headline, a per-CZ table, and the same "
      f"statistics for {len(SPECS) - 1} variants ({k_stats} sample x statistic cells in all).\n")
    a("## 2. What `ba_holders` measures (and why the data were pulled again)\n")
    a("The PREREG names `raw/acs_school_districts.parquet` (`ba_holders`, `pop_25plus`) as the input. `fetch_data.ACS_SD_VARS` maps "
      "`B15003_022E` to `ba_holders` and `B15003_001E` to `pop_25plus`. B15003 is *Educational attainment for the population 25 years and over*: "
      "one line per **highest** level completed. The Census metadata labels (read from the API at fetch time and stored in "
      "`raw/acs_b15003_labels_2018.json`) are: line 1 Total; line 22 **Bachelor's degree**; line 23 Master's degree; line 24 Professional school "
      "degree; line 25 Doctorate degree. `ba_holders` is therefore the number of adults whose highest degree is *exactly* a bachelor's: it "
      "**excludes** holders of graduate and professional degrees, so `ba_holders / pop_25plus` is *not* the \"BA-or-more\" share the PREREG defines. "
      "The headline therefore uses lines 22-25 of the same table, pulled for the same vintage (ACS 2018 5-year) and the same three school-district "
      "layers (`python3 ba_share_holdout.py fetch`, Census API, state by state). The exactly-bachelor's share is the variant `ba_only`.\n")
    a(f"Check of the new pull against the old file ({fn(pull['n_new'])} records each): {fn(pull['n_equal_pop'])} of {fn(pull['n_both'])} records have "
      f"`B15003_001E` = `pop_25plus` and {fn(pull['n_equal_ba'])} have `B15003_022E` = `ba_holders` (records only in one file: "
      f"{pull['n_only_new']} / {pull['n_only_old']}); lines 2-25 add up to line 1 in {fn(pull['n_lines_sum_to_total'])} records; "
      f"{pull['n_with_negative_code']} records carry an annotation code. **Departure from the PREREG wording:** the data file named there "
      "does not measure the concept it names; the concept (BA or more) is kept and the data are replaced by the full table.\n")
    a("## 3. Definitions\n")
    a("**Records.** One row per `leaid` (state FIPS + 5-digit district code) and layer (unified, elementary, secondary). The \"Remainder of <state>\" "
      f"rows (code 99999; {rec['n_remainder_dropped']} rows) are not districts and are dropped (as `district_bootstrap.load_income_columns`); "
      f"{fn(rec['n_unique'])} records remain, none with a repeated id. `adults` = `B15003_001E`; `BA+` = `B15003_022E + 023E + 024E + 025E`; a record is "
      "valid if the adult count and the four degree lines are published (no annotation code), there is at least one adult, and BA+ does not exceed the adults. "
      f"No record carries an annotation code and lines 2-25 add up to line 1 in every record; the {rec['n_invalid']} records that are not used are districts with no "
      "adults 25+ (\"School District Not Defined\", unorganized territories, a few elementary and secondary districts without residents).\n")
    a("**Panel and partition.** `district_settings.build_audited_panel` (New York City collapsed to one district under the Chancellor's id 3620580, "
      "under which the ACS files it) and `scheme_tables(A, \"locale\")`: city = NCES locale 11-13, suburb = 21-23, districts with at least 500 pupils, "
      f"each location holding at least 10% of the CZ's partitioned pupils: {cov['n_cz']} CZs, {fn(cov['n_districts'])} districts "
      f"({fn(cov['n_by_location'][1])} city, {fn(cov['n_by_location'][2])} suburb), {cov['pupils'] / 1e6:.2f} million pupils. The B15003 records are "
      "merged on `leaid` after the NYC collapse. CZ weights are the partitioned pupils, as for the four targets.\n")
    a("**Location shares and the CZ gap.** For CZ c and location l (1 = city, 2 = suburb), s_{c,l} = (sum of BA+ over the location's districts with a valid "
      "record) / (sum of adults 25+ over the same districts): a ratio of location totals, as every per-pupil object of `data_estimate.Gap` (not a mean of "
      "district shares). The CZ gap is s_{c,2} - s_{c,1}; a CZ with no valid record in one location has no gap and drops out of every statistic. "
      "**B1 = the CZ-weighted mean of the CZ gaps** (weights: partitioned pupils), a fraction; the sealed file also has it in percentage points, the "
      "CZ-weighted mean of log(s_{c,2} / s_{c,1}), the CZ-weighted means of the city and suburb shares over the CZs that have the gap (so that their "
      "difference is the gap), and the version pooled over CZs (sum of BA+ over sum of adults, per location, adult-weighted, then suburb minus city).\n")
    a("**Overlapping layers.** An adult who lives in an elementary district inside a secondary district appears in both records, so summing every layer "
      "double-counts adults (the unified and elementary records tile each state; the secondary layer overlays the elementary one). "
      f"The headline counts adults **once, from the unified and elementary records** (the rule of `teacher_residence.py` for employed residents and of "
      f"`children_income.py` for children); the {fn(cov['secondary_excluded']['n'])} secondary LEAs of the partition ({fp(cov['secondary_excluded']['pupil_share'], 2)} of "
      "its pupils) keep their pupils in the CZ weights and contribute no adults. This rule was fixed on those precedents before any share was "
      "computed. The rule of the four earnings targets (`load_income_columns`: every LEA takes the record of its own layer) is the variant `all_layers`; "
      "the two other layer rules (`unified_secondary`, `unified_only`) are in section 4.\n")
    a("**Matching.** ACS record to LEA on the 7-digit `leaid`. Detroit's record is filed under the pre-2016 id 2612000 (Detroit City School District) while "
      "the CCD frame carries its successor Detroit Public Schools Community District (2601103; the same boundary): it is refiled to 2601103, as in "
      "`children_income.py` (variant `no_detroit_crosswalk` leaves it out). No other crosswalk is applied: the other unmatched LEAs are listed in section 6 "
      "(non-geographic LEAs such as county offices, regional vocational districts and virtual schools have no ACS record).\n")
    a("## 4. Specifications (all sealed)\n")
    a("| key | description |\n|---|---|")
    for s in SPECS:
        a(f"| `{s.key}` | {s.label} |")
    a("\nEvery specification is evaluated on the full sample and on each hold-out half. The headline is `headline`; the variants are reported for robustness and "
      "are not used to choose among specifications.\n")
    a("## 5. Inference and hold-out\n")
    a(f"Whole-CZ bootstrap, B = {R['B']}, seed {R['seed']}, the multiplicity matrix of `district_settings.py` "
      f"(`district_bootstrap.draw_counts({cov['n_cz']}, B, seed)`): SE = standard deviation (ddof 1) of the draws, 95% interval = 2.5 and 97.5 percentiles, CZs with an "
      "undefined statistic drop out of the draw. The draws are therefore the same draws as those of the four district targets, and the bootstrap "
      "correlations with them are valid. The pre-registered split of the locale CZs (`district_settings.holdout_split`, seed 20261002; pairs of CZs "
      f"sorted by pupils, one to each half) is reproduced here ({cov['n_A']} CZs in A, {cov['n_B']} in B) and compared with the assignment stored in "
      "`district_settings.json`; each half is bootstrapped on its own number of CZs with the same seed, as `district_settings.py` does. "
      "Only the ACS sampling variance of the district counts is not added (section 9).\n")
    a("## 6. Sample and coverage\n")
    a(f"Partition: {cov['n_cz']} CZs, {fn(cov['n_districts'])} districts, {cov['pupils'] / 1e6:.2f} million partitioned pupils (CZ weights: effective number of CZs "
      f"{cov['n_eff_cz']:.1f}, five largest CZs {fp(cov['top5_weight_share'])} of the weight). ACS record layers of the partitioned districts: "
      + ", ".join(f"{k} {fn(v)} ({fp(cov['pupil_share_by_layer'][k])} of pupils)" for k, v in sorted(cov["n_by_layer"].items())) + ".\n")
    a("| specification | districts with a valid record | share of districts | share of partition pupils | city / suburb share of pupils | CZs with a gap (of "
      f"{cov['n_cz']}) | share of the CZ weight | in half A / B |\n|---|---|---|---|---|---|---|---|")
    for s in SPECS:
        v = cov["specs"][s.key]
        a(f"| `{s.key}` | {fn(v['n_valid'])} | {fp(v['share_districts_valid'])} | {fp(v['share_pupils_valid'])} | "
          f"{fp(v['share_pupils_valid_by_location'][1])} / {fp(v['share_pupils_valid_by_location'][2])} | {v['n_cz_gap']} | "
          f"{fp(v['weight_share_cz_gap'])} | {v['n_cz_gap_A']} / {v['n_cz_gap_B']} |")
    un = cov["unmatched"]
    a(f"\nDistricts of the partition outside the headline sample: {fn(un['n'])} ({fp(un['pupil_share'], 2)} of partition pupils). By reason: "
      + "; ".join(f"{k} {v['n']} ({fn(v['pupils'])} pupils)" for k, v in un["by_reason"].items())
      + ". Largest by pupils (public CCD names):\n")
    a("| LEA | leaid | pupils | locale | reason |\n|---|---|---|---|---|")
    for t in un["top"]:
        a(f"| {t['name']} | {t['leaid']} | {fn(t['pupils'])} | {t['locale']} | {t['reason']} |")
    a("\nLargest districts with no ACS record (non-geographic LEAs -- county offices, regional vocational and virtual districts, statewide "
      "turnaround districts -- and districts whose CCD and Census ids differ; none is refiled except Detroit):\n")
    a("| LEA | leaid | pupils | locale |\n|---|---|---|---|")
    for t in un["top_no_record"]:
        a(f"| {t['name']} | {t['leaid']} | {fn(t['pupils'])} | {t['locale']} |")
    cl = cov["cz_location_coverage"]
    a(f"\nCoverage by CZ location (share of the location's partitioned pupils whose district carries a valid headline record; {cl['n']} CZ locations): "
      f"minimum {fp(cl['min'])}, fifth percentile {fp(cl['p5'])}, median {fp(cl['median'])}; {cl['n_below_50']} below 50% and {cl['n_below_80']} below 80% "
      f"(their CZs carry {fp(cl['weight_share_cz_below_80'])} of the CZ weight). In those locations {fp(cl['secondary_share_of_unmatched_pupils_below_80'])} of the unmatched pupils are in "
      "secondary LEAs (whose pupils count as unmatched although the same households' adults are counted through the elementary districts) and the rest in LEAs "
      "without a record; the per-CZ coverage is in the sealed `cz_table` (`pupil_cov_city`, `pupil_cov_suburb`).")
    ad = cov["adults"]
    a(f"\nOverlap: the adults of all valid records of the partition exceed those counted once by the factor {ad['overlap_factor']:.4f} "
      f"(city {ad['overlap_factor_by_location'][1]:.4f}, suburb {ad['overlap_factor_by_location'][2]:.4f}); the headline counts {fn(ad['headline'])} adults, "
      f"the `all_layers` variant {fn(ad['all_layers'])}.\n")
    a("Tiling of the states by the unified and elementary records (adults 25+, `B15003_001E`): national "
      f"sum of unified + elementary records over the state totals of the same line = {til['national_unified_plus_elementary']:.4f}; with the secondary layer added "
      f"{til['national_all_layers']:.4f}; {til['n_states_ue_within_1pct']} of {til['n_states']} states within 1% and {til['n_states_ue_within_5pct']} within 5% "
      f"(smallest {til['min_ue']:.2f}, largest {til['max_ue']:.2f}). States with the largest double-count when all layers are summed: "
      + "; ".join(f"{t['state']} {t['all_layers']:.2f}" for t in til["top_states_all_layers"]) + ".\n")
    for nm, lab in (("nyc", "New York City"), ("detroit", "Detroit")):
        v = cov[nm]
        if v.get("in_partition"):
            a(f"{lab} (`{v['leaid']}`): in the partition (location {v['location']}, {fn(v['pupils'])} pupils, {fp(v['pupil_share_of_cz'])} of its CZ's partitioned "
              f"pupils, its CZ carries {fp(v['pupil_share_of_weights'], 2)} of the CZ weight); ACS record layer {v['record_layer']}; in the headline sample: {v['in_headline']}.")
    a("")
    a("## 7. Checks run before sealing (no sealed value involved)\n")
    a("`seal` runs every check first and refuses to seal if one fails; `check` prints the same lines.\n")
    a("```")
    for name, ok, detail in R["checks"]:
        a(f"  {name}: {ok_str(ok)}" + (f"  [{detail}]" if detail else ""))
    a("```\n")
    a("Reproduction of the four district targets on the same draws (value (bootstrap SE); the digits of `district_audit.md` / `district_bootstrap.md`):\n")
    a("| target | computed here | reported | n CZs |\n|---|---|---|---|")
    for r in R["repro"]:
        a(f"| `{r['moment']}` | {r['value']:.4f} ({r['se']:.4f}) | {r['reported_value']:.4f} ({r['reported_se']:.4f}) | {r['n_cz']} |")
    dc = R["draw_checks"]
    a(f"\nThe {R['B']} x {cov['n_cz']} multiplicity matrix reproduces the draws of the four targets: maximum absolute difference between the draws recomputed from "
      f"the CZ table and those of `moment_set` {max(dc['recomputed_max_abs_diff'].values()):.1e}"
      + ("" if dc["csv_max_abs_diff"] is None else f", against the stored `district_settings_draws_locale.csv` ({dc['csv_rows']} draws) "
         f"{max(dc['csv_max_abs_diff'].values()):.1e} (the file is rounded to 10 decimals)") + ". "
      f"The two code paths that B1 uses are validated on public quantities: the weighted-mean/bootstrap routine reproduces the income target "
      f"(difference {R['machinery']['max_abs_diff_wmean']:.1e}) and the CZ-table routine the enrollment gap column ({R['machinery']['max_abs_diff_cz_table']:.1e}). "
      "`ba_share_holdout_xcheck.py` recomputes every sealed statistic in separately written code (including the `ba_only` variant from the old parquet) and prints only "
      "PASS/FAIL and the maximum absolute differences.\n")
    a("## 8. The sealed file\n")
    a("`estimates/settings/sealed/ba_share_holdout_B1.json`: `meta` (B, seeds, code and input SHA-256), `definition`, `specs`, `headline` "
      "(`ba_share_gap` and `pp` for `full`, `A`, `B`: value, se, ci95_lo, ci95_hi, n_cz), `samples` (`full`, `A`, `B`: `n_cz`, `cz_ids`, `stats` with keys "
      f"{', '.join('`' + k + '`' for k in STAT_KEYS)}, and `corr`: the bootstrap correlation matrix of the draws of "
      f"{', '.join('`' + k + '`' for k in TARGETS + list(DRAW_KEYS))}), `variants` (the same `stats` per sample and specification, plus the correlation of the "
      "gap with the four targets), `draws_full` (the B draws of the headline gap, log ratio and pooled gap) and `cz_table` (per CZ: weight, location shares, "
      "BA+ and adult totals, gap, pupil coverage of the two locations). `python3 ba_share_holdout.py verify-seal` prints the structure and consistency checks (finite, shares in [0, 1], "
      "pp = 100 x fraction, city + gap = suburb, halves partition the CZs), never a value.\n")
    a("## 9. Caveats\n")
    a("- The ACS 5-year school-district tabulation is a sample (2014-18, about 5% of households over five years): the district counts have sampling error "
      "(margins of error are not pulled). The SE reported is the whole-CZ bootstrap SE, as for the four district targets; it does not add the ACS sampling "
      "variance of the district counts, which is small relative to the variation across CZs but not zero for small districts.")
    a("- The location of an adult is the locale of the school district he lives in, not of his household's neighbourhood; adults without children are included "
      "(B15003 is the whole population 25+, including group quarters). B1 is the education sorting of *residents*, not of parents.")
    a("- The headline counts adults from the unified and elementary records. A CZ location made up only of secondary-layer LEAs has no adults under the headline "
      "(the counts of such CZs are in section 6; `unified_secondary` and `all_layers` bracket the rule). Dual-system states (elementary + secondary districts) "
      "are where the overlap matters.")
    a("- The partition is on 2017-18 CCD districts of at least 500 pupils; the ACS boundaries are those of 2017. Non-geographic LEAs and LEAs whose "
      "id differs between the CCD and the ACS (apart from Detroit) have no record and leave the sample (section 6).")
    a("- CZ weights are the partitioned pupils (the four targets' weights), not adults. The `adult_weights` and `equal_weights` variants show the sensitivity.")
    a("- `ba_holders` of the existing parquet is the exactly-bachelor's count (section 2). `fetch_data.py` writes `ba_holders` and `pop_25plus`; no other script of "
      "the repository reads them (checked by `grep`), so no earlier estimate used them as a college share.")
    a("- This file reports no value; the checks above use only public quantities (the four reported district targets, counts, and Census totals of adults).\n")
    a("## 10. Disclosure log\n")
    a("What was looked at before sealing, and what the runs printed. (a) Structure inspection of `raw/acs_school_districts.parquet`, run twice: column names, dtypes, layer counts and "
      "its first three rows (Choctaw County, Claiborne County and Clarksdale Municipal, Mississippi: median income, home value, population, `ba_holders` and `pop_25plus`, "
      "i.e. the exactly-bachelor's counts of three small districts); no aggregate of any BA share was computed from it by location, CZ or sample. (b) The Census metadata labels of B15003 lines 1 and 22-25 (section 2). (c) `fetch` printed record counts only. (d) The comparison of the new "
      "pull with the old file printed counts of equal records. (e) The partition was matched to the ACS records and the unmatched LEAs listed by pupils "
      "(names and pupils: public inputs, section 6); the adult count (`B15003_001E`, a total, not a degree count) of the 2612000 record was printed to confirm that it is "
      "Detroit's before the refile was fixed; no degree count entered that decision. (f) `check` and `seal` print sample sizes, "
      "coverage shares of pupils, the reproduction of the four public targets and PASS/FAIL lines; nothing computed from a degree count is printed, logged or written "
      "outside the sealed file. (g) The sealed file was verified by `verify-seal` and the independent `ba_share_holdout_xcheck.py` (structure, PASS/FAIL, maximum absolute "
      "differences). (h) `seal` scans this file, and `leak-scan` the run logs, for renderings of the sealed numbers: the high-precision ones (fractions at 5-6 decimals, "
      "percentage points at 4) must not appear, and the coarse ones (fractions at 4 decimals, percentage points at 2-3, which collide with public numbers by chance) must "
      "not appear more often than for random decoys of the sealed numbers; the counts are in the run record. This file is generated from `prepare` output, never from a statistic.")
    if sha:
        a(f"\nSHA-256 of the sealed file at the time of writing: `{sha}`.")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------

def print_public(R: dict[str, Any]) -> None:
    cov = R["coverage"]
    print(f"partition: {cov['n_cz']} CZs, {cov['n_districts']:,} districts ({cov['n_by_location'][1]:,} city, {cov['n_by_location'][2]:,} suburb), "
          f"{cov['pupils'] / 1e6:.2f} m pupils; hold-out halves {cov['n_A']} / {cov['n_B']} CZs")
    print("ACS record layers of the partitioned districts: " + ", ".join(f"{k} {v:,}" for k, v in sorted(cov["n_by_layer"].items())))
    print(f"{'specification':<22}{'valid districts':>16}{'pupils matched':>16}{'CZs with gap':>14}{'CZ weight':>11}{'A':>5}{'B':>5}")
    for s in SPECS:
        v = cov["specs"][s.key]
        print(f"{s.key:<22}{v['n_valid']:>16,}{fp(v['share_pupils_valid'], 2):>16}{v['n_cz_gap']:>14}{fp(v['weight_share_cz_gap']):>11}"
              f"{v['n_cz_gap_A']:>5}{v['n_cz_gap_B']:>5}")
    un = cov["unmatched"]
    print(f"outside the headline sample: {un['n']:,} districts, {fp(un['pupil_share'], 2)} of pupils; "
          + "; ".join(f"{k}: {v['n']}" for k, v in un["by_reason"].items()))
    cl = cov["cz_location_coverage"]
    print(f"pupil coverage by CZ location (headline): min {fp(cl['min'])}, p5 {fp(cl['p5'])}, median {fp(cl['median'])}; "
          f"{cl['n_below_50']} of {cl['n']} below 50%, {cl['n_below_80']} below 80%")
    ad = cov["adults"]
    print(f"adults: headline {ad['headline']:,.0f}, all layers {ad['all_layers']:,.0f} (overlap factor {ad['overlap_factor']:.4f})")
    til = R["tiling"]
    print(f"tiling: unified + elementary / state totals {til['national_unified_plus_elementary']:.4f}, with secondary {til['national_all_layers']:.4f}; "
          f"{til['n_states_ue_within_1pct']} of {til['n_states']} states within 1%")
    print("largest districts outside the headline sample: " + "; ".join(f"{t['name'][:34]} ({t['pupils']:,.0f})" for t in un["top"][:5]))
    for nm in ("nyc", "detroit"):
        v = cov[nm]
        if v.get("in_partition"):
            print(f"{nm}: location {v['location']}, record layer {v['record_layer']}, in headline {v['in_headline']}")
    bad = [c for c in R["checks"] if not c[1]]
    print(f"checks: {len(R['checks']) - len(bad)} of {len(R['checks'])} PASS" + ("" if not bad else "; FAIL: " + ", ".join(c[0] for c in bad)))


def cmd_check(args: argparse.Namespace) -> int:
    R = prepare(args.B, args.seed)
    print_public(R)
    print(f"all checks: {ok_str(R['checks_ok'])}  ({R['seconds']:.0f}s)")
    return 0 if R["checks_ok"] else 1


def cmd_seal(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).resolve() if args.outdir else EST
    target = outdir / "sealed" / SEALED_NAME
    if target.exists() and not args.force:
        print(f"refusing to overwrite the sealed file {target} (use --force)", file=sys.stderr)
        return 2
    R = prepare(args.B, args.seed)
    print_public(R)
    if not R["checks_ok"]:
        print("checks failed (see above); nothing sealed", file=sys.stderr)
        return 1
    res = run_all(R)
    payload = build_payload(R, res)
    write_sealed(target, payload)
    k = n_statistics(payload)
    sha = sha256_file(target)
    md = outdir / METHODS_NAME
    write_methods(md, R, k, sha)
    if not scan_report(md.name, md.read_text(), payload):
        md.unlink()
        print("a sealed number appears in the methods file: it was deleted; the sealed file is kept (re-run with --force after fixing)",
              file=sys.stderr)
        return 1
    print(f"sealed: {target} ({k} statistics)")
    print(f"sha256: {sha}")
    return 0


def cmd_leak_scan(args: argparse.Namespace) -> int:
    """Scan text files (run logs, methods file) for renderings of the sealed numbers; counts and PASS/FAIL only."""
    outdir = Path(args.outdir).resolve() if args.outdir else EST
    payload = json.loads((outdir / "sealed" / SEALED_NAME).read_text())
    ok = True
    for f in args.files:
        ok &= scan_report(Path(f).name, Path(f).read_text(errors="replace"), payload)
    return 0 if ok else 1


def cmd_verify(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).resolve() if args.outdir else EST
    return verify_seal(outdir / "sealed" / SEALED_NAME)


# --------------------------------------------------------------------------
# Self-test (synthetic data only)
# --------------------------------------------------------------------------

def selftest() -> int:
    import tempfile
    _stack()
    bad: list[str] = []

    def check(name: str, cond: bool) -> None:
        print(f"  {name}: {ok_str(bool(cond))}")
        if not cond:
            bad.append(name)

    # labels
    good = {"B15003_001E": "Estimate!!Total:", "B15003_022E": "Estimate!!Total:!!Bachelor's degree",
            "B15003_023E": "Estimate!!Total:!!Master's degree", "B15003_024E": "Estimate!!Total:!!Professional school degree",
            "B15003_025E": "Estimate!!Total:!!Doctorate degree"}
    check("label check accepts the Census labels", _ok(lambda: check_labels(good)))
    check("label check rejects a shifted line", not _ok(lambda: check_labels({**good, "B15003_022E": "Estimate!!Total:!!Master's degree"})))

    # records: remainder rows, annotation codes, BA+ sum, invalid
    n = 6
    raw = pd.DataFrame({"NAME": [f"d{i}" for i in range(n)],
                        "leaid": ["0100001", "0100002", "0100003", "0199999", "0199999", "0100004"],
                        "sd_layer": ["school district (unified)", "school district (elementary)", "school district (secondary)",
                                     "school district (unified)", "school district (elementary)", "school district (unified)"]})
    for j in range(1, N_LINES + 1):
        raw[f"B15003_{j:03d}E"] = 0.0
    raw["B15003_001E"] = [100.0, 200.0, 50.0, 10.0, 10.0, 100.0]
    raw["B15003_022E"] = [10.0, 40.0, 5.0, 1.0, 1.0, 10.0]
    raw["B15003_023E"] = [5.0, 10.0, 5.0, 0.0, 0.0, -666666666.0]       # annotation code on the last record
    raw["B15003_024E"] = [1.0, 2.0, 0.0, 0.0, 0.0, 1.0]
    raw["B15003_025E"] = [1.0, 3.0, 0.0, 0.0, 0.0, 1.0]
    rec = records_from_raw(raw)
    check("remainder rows (99999) dropped", len(rec) == 4 and rec.attrs["n_remainder_dropped"] == 2)
    r = rec.set_index("leaid")
    check("BA+ = lines 22-25, BA only = line 22", r.ba_plus["0100001"] == 17 and r.ba_only["0100001"] == 10
          and r.ba_plus["0100002"] == 55 and r.adults["0100002"] == 200)
    check("a record with an annotation code in lines 22-25 is invalid", np.isnan(r.ba_plus["0100004"]) and np.isnan(r.adults["0100004"])
          and rec.attrs["n_invalid"] == 1)

    # merge with the Detroit refile
    A = pd.DataFrame({"leaid": ["0100001", "2601103", "3620580"], "x": [1, 2, 3]})
    c = pd.DataFrame({"leaid": ["0100001", "2612000", "3620580"], "name": ["a", "detroit", "nyc"],
                      "layer": ["unified", "unified", "unified"], "adults": [100.0, 400.0, 900.0],
                      "ba_plus": [20.0, 100.0, 300.0], "ba_only": [10.0, 60.0, 200.0]})
    X = add_b15003(A, c)
    check("Detroit's record is refiled to 2601103", float(X.b_adults[X.leaid == "2601103"].iloc[0]) == 400.0
          and X.b_name[X.leaid == "2601103"].iloc[0] == "detroit")

    # a toy partition: two CZs; hand-computed answers
    #   CZ 1: city  = unified a (100 adults, 40 BA+, 25 BA only) + secondary b (100, 50, 30) ;
    #         suburb = elementary c (200, 120, 100) + unified d (100, 20, 10)
    #   CZ 2: city  = unified e (50, 5, 5);  suburb = unified f (50, 10, 10)
    rows = [("a", 1, 1.0, 1000, "unified", 100, 40, 25), ("b", 1, 1.0, 400, "secondary", 100, 50, 30),
            ("c", 1, 2.0, 600, "elementary", 200, 120, 100), ("d", 1, 2.0, 400, "unified", 100, 20, 10),
            ("e", 2, 1.0, 500, "unified", 50, 5, 5), ("f", 2, 2.0, 500, "unified", 50, 10, 10)]
    part = pd.DataFrame(rows, columns=["leaid", "czone", "loc_id", "pupils", "b_layer", "b_adults", "b_ba", "b_ba_only"])
    part["year"] = YEAR
    part["lea_name"] = part.leaid
    g = pd.DataFrame({"czone": [1.0, 2.0], "pupils": [2400.0, 1000.0]})
    part["czone"] = part.czone.astype(float)
    # spec: headline -> adults once (b is secondary and drops out)
    gz = spec_table(part, g, SPEC_BY_KEY["headline"])
    s1, s2 = 40 / 100, (120 + 20) / 300           # CZ 1 city, suburb
    t1, t2 = 5 / 50, 10 / 50                       # CZ 2
    check("headline: location shares are ratios of location totals, secondary records out",
          np.allclose(gz.b1_share_1, [s1, t1]) and np.allclose(gz.b1_share_2, [s2, t2]))
    check("headline: gap = suburb minus city", np.allclose(gz.gap_b1_share, [s2 - s1, t2 - t1]))
    check("headline: log ratio", np.allclose(gz.gap_b1_logratio, [np.log(s2 / s1), np.log(t2 / t1)]))
    gz_all = spec_table(part, g, SPEC_BY_KEY["all_layers"])
    check("all_layers: every record counted", np.isclose(gz_all.b1_share_1[0], (40 + 50) / 200)
          and np.isclose(gz_all.b1_share_2[0], 140 / 300))
    gz_us = spec_table(part, g, SPEC_BY_KEY["unified_secondary"])
    check("unified_secondary: the elementary record is dropped where a secondary record exists (CZ 1 city has b; the suburb has none)",
          np.isclose(gz_us.b1_share_1[0], (40 + 50) / 200) and np.isclose(gz_us.b1_share_2[0], 140 / 300))
    gz_u = spec_table(part, g, SPEC_BY_KEY["unified_only"])
    check("unified_only: elementary and secondary records leave the sums", np.isclose(gz_u.b1_share_1[0], 40 / 100)
          and np.isclose(gz_u.b1_share_2[0], 20 / 100))
    gz_bo = spec_table(part, g, SPEC_BY_KEY["ba_only"])
    check("ba_only: line 22 only", np.isclose(gz_bo.b1_share_1[0], 25 / 100) and np.isclose(gz_bo.b1_share_2[0], (100 + 10) / 300))
    gz_nz = spec_table(part, g, SPEC_BY_KEY["drop_ny_cz"])
    check("drop_ny_cz without an NYC row leaves the table unchanged", np.allclose(gz_nz.gap_b1_share, gz.gap_b1_share))
    part2 = part.copy()
    part2.loc[part2.leaid == "e", "leaid"] = NYC_LEAID
    gz_nyc = spec_table(part2, g, SPEC_BY_KEY["drop_ny_cz"])
    check("drop_ny_cz removes the NYC CZ", np.isnan(gz_nyc.gap_b1_share[1]) and np.isfinite(gz_nyc.gap_b1_share[0]))
    gz_nd = spec_table(part2, g, SPEC_BY_KEY["nyc_dropped"])
    check("nyc_dropped removes NYC's adults (the CZ loses its city location)", np.isnan(gz_nd.gap_b1_share[1]))

    # statistics: weighted mean, pooled, bootstrap draws against hand computation
    B_ = 40
    counts = db.draw_counts(2, B_, 7)
    st, dr = sample_stats(gz, g.pupils.to_numpy(float), counts)
    w = g.pupils.to_numpy(float)
    gap = np.array([s2 - s1, t2 - t1])
    est = float((w * gap).sum() / w.sum())
    draws = (counts @ (w * gap)) / (counts @ w)
    check("B1 = pupil-weighted mean of the CZ gaps", np.isclose(st["ba_share_gap"]["value"], est))
    check("B1 draws = multiplicity matrix x weights", np.allclose(dr["ba_share_gap"], draws))
    check("B1 SE = std (ddof 1) of the draws", np.isclose(st["ba_share_gap"]["se"], np.std(draws[np.isfinite(draws)], ddof=1)))
    check("pp = 100 x fraction", np.isclose(st["ba_share_gap_pp"]["value"], 100 * est) and np.isclose(st["ba_share_gap_pp"]["se"], 100 * st["ba_share_gap"]["se"]))
    check("city and suburb means: suburb - city = gap", np.isclose(st["ba_share_suburb"]["value"] - st["ba_share_city"]["value"], est))
    pc = (40 + 5) / (100 + 50)
    ps = (140 + 10) / (300 + 50)
    check("pooled shares are ratios of sums over CZs", np.isclose(st["pooled_city"]["value"], pc) and np.isclose(st["pooled_suburb"]["value"], ps)
          and np.isclose(st["pooled_gap"]["value"], ps - pc))
    # a CZ without a gap drops out, with the multiplicities kept
    gz_nan = gz.copy()
    gz_nan.loc[1, ["b1_share_1", "b1_share_2", "gap_b1_share", "gap_b1_logratio"]] = np.nan
    st2, _ = sample_stats(gz_nan, w, counts)
    check("a CZ without a gap drops out of the mean", np.isclose(st2["ba_share_gap"]["value"], s2 - s1) and st2["ba_share_gap"]["n_cz"] == 1)
    # corr
    tdraws = pd.DataFrame({k: np.random.default_rng(i).normal(size=B_) for i, k in enumerate(TARGETS)})
    cb = corr_block(dr, tdraws)
    check("corr block has the four target correlations", set(cb["ba_share_gap_with"]) == set(TARGETS))

    # sealing: no overwrite, structure, leak scan
    with tempfile.TemporaryDirectory() as td:
        payload = dict(sealed=True, headline=dict(full=st["ba_share_gap"]), samples={"full": dict(stats=st)}, variants={})
        p = Path(td) / "x.json"
        write_sealed(p, ds.clean(payload))
        check("sealed json written without NaN", p.exists() and "NaN" not in p.read_text())
    toks_payload = ds.clean(dict(headline=dict(full=dict(value=0.123456, se=0.0123456)), samples={}, variants={}))
    check("leak scan finds a sealed number", leak_scan("the gap is 0.12346 here", toks_payload) >= 1)
    check("leak scan ignores unrelated text and low-precision renderings",
          leak_scan("0.5 and 12 and 0.99999 and 0.1235 and 12.35 pp", toks_payload) == 0)
    check("leak scan finds a percentage-point rendering", leak_scan("12.3456 pp", toks_payload) >= 1)
    check("leak scan does not match inside a longer number", leak_scan("0.123461 and 10.12346", toks_payload) == 0)
    check("coarse scan: a 4-decimal print is found and exceeds the decoys",
          not coarse_scan("the gap is 0.1235 and the se 0.0123", toks_payload)["ok"])
    check("coarse scan: unrelated numbers do not exceed the decoys", coarse_scan("shares 96.5% and 28.55 m", toks_payload)["ok"])
    check("verify-able file scan: leak-scan reads the sealed json layout",
          len(sealed_floats(toks_payload)) == 2)

    print(f"selftest: {len(bad)} failure(s)")
    return 1 if bad else 0


def _ok(fn_: Any) -> bool:
    try:
        fn_()
        return True
    except SystemExit:
        return False


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    f = sub.add_parser("fetch", help="download and cache B15003 by school district, the state totals and the labels")
    f.add_argument("--refresh", action="store_true", help="download again even if cached")
    f.add_argument("--workers", type=int, default=8)
    for name, hlp in (("check", "coverage and every check that involves no sealed value"),
                      ("seal", "compute, write the sealed json and the methods file (prints one line and the SHA-256)")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("--B", type=int, default=B_DEFAULT)
        p.add_argument("--seed", type=int, default=SEED)
        p.add_argument("--outdir", default=None, help="default: estimates/settings")
        if name == "seal":
            p.add_argument("--force", action="store_true", help="overwrite an existing sealed file")
    v = sub.add_parser("verify-seal", help="structure of the sealed file, no values")
    v.add_argument("--outdir", default=None)
    ls = sub.add_parser("leak-scan", help="scan text files for renderings of the sealed numbers (counts and PASS/FAIL only)")
    ls.add_argument("files", nargs="+")
    ls.add_argument("--outdir", default=None)
    sub.add_parser("selftest", help="unit checks on synthetic data")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    for n in ("data_estimate", "district_audit"):
        logging.getLogger(n).setLevel(logging.WARNING)
    if args.cmd == "fetch":
        return cmd_fetch(args)
    if args.cmd == "check":
        return cmd_check(args)
    if args.cmd == "seal":
        return cmd_seal(args)
    if args.cmd == "verify-seal":
        _stack()
        return cmd_verify(args)
    if args.cmd == "leak-scan":
        return cmd_leak_scan(args)
    if args.cmd == "selftest":
        return selftest()
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
