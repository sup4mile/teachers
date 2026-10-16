#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
commuting_all.py -- residence against workplace of ALL workers and of NON-TEACHER workers (task D1b).

`teacher_residence.py` tabulates, for K-12 teachers, the class (city / suburb / other) of the POWPUMA of
residence against the class of the POWPUMA of work (its `flows` block; section "Residence against
workplace of K-12 teachers on one geography" of teacher_residence.md).  This script builds the same
table for (i) all employed persons 16+ and (ii) all non-teacher workers (everyone except OCCP
2310-2330), on the same POWPUMA classes, with the same replicate-weight standard errors.

It IMPORTS teacher_residence.py (principal-city / suburban classification of the PUMAs and POWPUMAs,
the PUMS loader, the successive-difference-replicate engine `Ratios`) and changes nothing in it.  The
K-12 teacher rows are recomputed by the same code and must reproduce `flows.k12` / `flows.pub_k12` of
estimates/settings/teacher_residence.json (the script stops if they do not).

Sample: ACS 2018 1-year PUMS, employed persons 16+ in households with a place-of-work POWPUMA in the
United States (POWSP 1-56, POWPUMA > 1; the sample of `flows_block`; workers who work at home carry
their own area as the place of work).

Statistics for each group g (shares of the group's workers; PWGTP; 80 replicate weights, SE =
sqrt(4/80 sum (x_r - x)^2); derived statistics are evaluated on every replicate):
    res / work   shares living / working in a city, suburb, other POWPUMA
    x[a->b]      the 3x3 residence x workplace matrix (joint shares)
    p[b|a]       P(work in b | live in a);  q[a|b] = P(live in a | work in b)
    cross        P(work in suburb | live in city), P(work in city | live in suburb), the share crossing
                 between city and suburb among those who live AND work in a city or suburb POWPUMA
                 (`share_metro_cross` of teacher_residence.py), the share working in another class
                 (`share_cross_class`)
    block        among workers who live and work in city or suburb: the city share of jobs, the city
                 share of residents, the net inflow of workers into the city
    marginal     the city share of residents among workers living in city or suburb, the city share of
                 jobs among workers working in city or suburb (the margins `gap_res` / `gap_work` use)
    gaps         log(suburb/city) by residence and by workplace and their difference (net commuting gap)
    ... for all workers, for metropolitan-area residents only, for workers away from home, and under
    the other class cuts of `teacher_residence.CUTS`.

Usage
    python commuting_all.py                  # build (about 2 minutes, a compute node; peaks near 8 GB)
    python commuting_all.py selftest         # unit check of the statistics registry (no data)
    python commuting_all.py render           # rewrite the .md from the existing .json (no data)
    python commuting_all.py build --outdir <dir>
Writes estimates/settings/commuting_all.{md,json}.

Written with help from Claude Code.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import district_settings as ds  # noqa: E402  (ds.clean: JSON-safe copy)
import teacher_residence as tr  # noqa: E402

LOG = logging.getLogger("commuting_all")

CLS = tr.CLASSES                         # ('city', 'suburb', 'other')
CITY, SUB = "city", "suburb"
METRO_ROWS = (CITY, SUB)
REF_JSON = HERE / "estimates" / "settings" / "teacher_residence.json"

GROUP_LABEL = {
    "all": "all workers",
    "nonteacher": "non-teacher workers",
    "k12": "K-12 teachers (reference)",
    "pub_k12": "public K-12 teachers (reference)",
}
GROUP_NOTE = {
    "all": "employed persons 16+ with a US POWPUMA of residence and of work",
    "nonteacher": "all of the above except OCCP 2310-2330 (elementary/middle, secondary and special-education "
                  "teachers; OCCP 2300 preschool/kindergarten teachers stay in)",
    "k12": "OCCP 2310-2330 (`k12` of teacher_residence.py), recomputed here as a check on the code",
    "pub_k12": "OCCP 2310-2330 employed by a local or state government (COW 3-4), recomputed as a check",
}
UNIV_LABEL = {
    "all": "all workers with a US POWPUMA of residence and of work",
    "metro_pow": "residents of metropolitan areas: residence POWPUMA class city or suburb (rows city and suburb of the matrix)",
    "metro_puma": "residents of metropolitan PUMAs: residence PUMA class city or suburb (the `metro` flag of the "
                  "commute block of teacher_residence.py), cells on the POWPUMA classes",
    "away": "workers who work away from home (JWTR != 11)",
}


# --------------------------------------------------------------------------
# The statistics registry: every statistic is a function of the cell shares
# --------------------------------------------------------------------------

Getter = Callable[[str], Any]


def cell_key(key: str, a: str, b: str) -> str:
    return f"{key}|x|{a}|{b}"


def registry(rows: tuple[str, ...] = CLS) -> dict[str, Callable[[Getter, str], Any]]:
    """name -> f(get, key).  `get(name)` returns a cell share (the residence class a, the workplace class b) of the
    group `key`; `rows` are the residence classes in the universe (all three, or the metropolitan city and suburb).
    Every statistic is a ratio of sums of cells, so it does not depend on what the cells are shares of."""
    def x(g: Getter, k: str, a: str, b: str) -> Any:
        return g(cell_key(k, a, b))

    def tot(g: Getter, k: str) -> Any:
        return sum(x(g, k, a, b) for a in rows for b in CLS)

    def res_(g: Getter, k: str, a: str) -> Any:            # workers living in a
        return sum(x(g, k, a, b) for b in CLS)

    def work_(g: Getter, k: str, b: str) -> Any:           # workers working in b (from the universe's rows)
        return sum(x(g, k, a, b) for a in rows)

    def blk(g: Getter, k: str) -> Any:                     # live AND work in city or suburb
        return sum(x(g, k, a, b) for a in METRO_ROWS for b in METRO_ROWS)

    F: dict[str, Callable[[Getter, str], Any]] = {}
    for a in rows:
        F[f"res.{a}"] = lambda g, k, a=a: res_(g, k, a) / tot(g, k)
    for b in CLS:
        F[f"work.{b}"] = lambda g, k, b=b: work_(g, k, b) / tot(g, k)
    for a in rows:
        for b in CLS:
            F[f"x.{a}.{b}"] = lambda g, k, a=a, b=b: x(g, k, a, b) / tot(g, k)
            F[f"p.{a}.{b}"] = lambda g, k, a=a, b=b: x(g, k, a, b) / res_(g, k, a)        # P(work in b | live in a)
            F[f"q.{b}.{a}"] = lambda g, k, a=a, b=b: x(g, k, a, b) / work_(g, k, b)       # P(live in a | work in b)
    F["p_work_suburb_given_live_city"] = lambda g, k: x(g, k, CITY, SUB) / res_(g, k, CITY)
    F["p_work_city_given_live_suburb"] = lambda g, k: x(g, k, SUB, CITY) / res_(g, k, SUB)
    F["share_cross_city_suburb"] = lambda g, k: (x(g, k, CITY, SUB) + x(g, k, SUB, CITY)) / blk(g, k)
    F["share_cross_class"] = lambda g, k: 1.0 - sum(x(g, k, a, a) for a in rows) / tot(g, k)
    F["block.share_of_group"] = lambda g, k: blk(g, k) / tot(g, k)
    F["block.city_share_jobs"] = lambda g, k: (x(g, k, CITY, CITY) + x(g, k, SUB, CITY)) / blk(g, k)
    F["block.city_share_residents"] = lambda g, k: (x(g, k, CITY, CITY) + x(g, k, CITY, SUB)) / blk(g, k)
    F["block.net_inflow_to_city"] = lambda g, k: (x(g, k, SUB, CITY) - x(g, k, CITY, SUB)) / blk(g, k)
    F["marginal.city_share_residents"] = lambda g, k: res_(g, k, CITY) / (res_(g, k, CITY) + res_(g, k, SUB))
    F["marginal.city_share_jobs"] = lambda g, k: work_(g, k, CITY) / (work_(g, k, CITY) + work_(g, k, SUB))
    F["gap_res"] = lambda g, k: np.log(res_(g, k, SUB)) - np.log(res_(g, k, CITY))
    F["gap_work"] = lambda g, k: np.log(work_(g, k, SUB)) - np.log(work_(g, k, CITY))
    F["net_commute_gap"] = lambda g, k: F["gap_res"](g, k) - F["gap_work"](g, k)
    return F


# the statistics compared across groups and across class cuts
HEADLINE = ["p_work_suburb_given_live_city", "p_work_city_given_live_suburb", "share_cross_city_suburb",
            "share_cross_class", "block.city_share_jobs", "block.city_share_residents", "block.net_inflow_to_city",
            "marginal.city_share_residents", "marginal.city_share_jobs", "gap_res", "gap_work", "net_commute_gap"]


def flow_block(R: tr.Ratios, key: str, rows: tuple[str, ...] = CLS) -> dict[str, Any]:
    """All statistics of one group in one universe, with replicate-weight standard errors."""
    F = registry(rows)

    def D(name: str) -> dict[str, float]:
        return R.derive(lambda g: F[name](g, key))

    out: dict[str, Any] = dict(
        res={a: D(f"res.{a}") for a in rows},
        work={b: D(f"work.{b}") for b in CLS},
        x={f"{a}->{b}": D(f"x.{a}.{b}") for a in rows for b in CLS},
        work_given_live={a: {b: D(f"p.{a}.{b}") for b in CLS} for a in rows},
        live_given_work={b: {a: D(f"q.{b}.{a}") for a in rows} for b in CLS})
    for name in ("p_work_suburb_given_live_city", "p_work_city_given_live_suburb", "share_cross_city_suburb",
                 "share_cross_class", "gap_res", "gap_work", "net_commute_gap"):
        out[name] = D(name)
    out["block"] = {n: D(f"block.{n}") for n in ("share_of_group", "city_share_jobs", "city_share_residents",
                                                  "net_inflow_to_city")}
    out["marginal"] = {n: D(f"marginal.{n}") for n in ("city_share_residents", "city_share_jobs")}
    return out


def diff_block(R: tr.Ratios, k1: str, k2: str, rows: tuple[str, ...] = CLS) -> dict[str, Any]:
    """Statistic(k1) - statistic(k2), standard errors from the replicates (covariance respected)."""
    F = registry(rows)
    return {name: R.derive(lambda g, name=name: F[name](g, k1) - F[name](g, k2)) for name in HEADLINE}


def cell_specs(rc: np.ndarray, wc: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Cell shares of each masked population: weighted share of its persons who live in class a and work in b."""
    specs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for key, m in masks.items():
        den = m.astype(float)
        for a in CLS:
            ma = m & (rc == a)
            for b in CLS:
                specs[cell_key(key, a, b)] = ((ma & (wc == b)).astype(float), den)
    return specs


def selftest() -> None:
    """Hand-computed check of the registry on a known 3x3 flow matrix (rows: live city / suburb / other, columns:
    work city / suburb / other)."""
    N = np.array([[50, 10, 2], [6, 120, 4], [1, 3, 30]])
    pairs = [(a, b, N[i, j]) for i, a in enumerate(CLS) for j, b in enumerate(CLS)]
    rc = np.concatenate([[a] * n for a, b, n in pairs])
    wc = np.concatenate([[b] * n for a, b, n in pairs])
    n = len(rc)
    rng = np.random.default_rng(3)
    W = rng.integers(1, 5, (n, tr.N_REP))
    R = tr.Ratios(tr.Rep(np.ones(n), W), cell_specs(rc, wc, {"t": np.ones(n, bool)}))
    # the point estimate uses the unit weights, so the cells are the counts
    F = registry(CLS)
    g = lambda name: R.est[R.ix[name]]                                          # noqa: E731
    T = N.sum()
    assert np.isclose(F["p_work_suburb_given_live_city"](g, "t"), 10 / 62)
    assert np.isclose(F["p_work_city_given_live_suburb"](g, "t"), 6 / 130)
    assert np.isclose(F["share_cross_city_suburb"](g, "t"), 16 / 186)
    assert np.isclose(F["share_cross_class"](g, "t"), 1 - 200 / T)
    assert np.isclose(F["block.city_share_jobs"](g, "t"), 56 / 186)
    assert np.isclose(F["block.city_share_residents"](g, "t"), 60 / 186)
    assert np.isclose(F["block.net_inflow_to_city"](g, "t"), -4 / 186)
    assert np.isclose(F["marginal.city_share_residents"](g, "t"), 62 / 192)
    assert np.isclose(F["marginal.city_share_jobs"](g, "t"), 57 / 190)
    assert np.isclose(F["gap_res"](g, "t"), np.log(130 / 62)) and np.isclose(F["gap_work"](g, "t"), np.log(133 / 57))
    assert np.isclose(F["p.city.other"](g, "t"), 2 / 62) and np.isclose(F["q.suburb.city"](g, "t"), 10 / 133)
    assert np.isclose(F["res.other"](g, "t"), 34 / T) and np.isclose(F["work.other"](g, "t"), 36 / T)
    Fm = registry(METRO_ROWS)                                                   # metropolitan residents only
    assert np.isclose(Fm["res.city"](g, "t"), 62 / 192) and np.isclose(Fm["work.other"](g, "t"), 6 / 192)
    assert np.isclose(Fm["p_work_suburb_given_live_city"](g, "t"), 10 / 62)       # invariant to the restriction
    assert np.isclose(Fm["share_cross_class"](g, "t"), 1 - 170 / 192)
    d = diff_block(R, "t", "t")
    assert all(abs(v["value"]) < 1e-15 and abs(v["se"]) < 1e-15 for v in d.values())
    blk = flow_block(R, "t", CLS)
    assert blk["x"]["city->suburb"]["se"] > 0 and abs(sum(v["value"] for v in blk["res"].values()) - 1) < 1e-12
    print("selftest ok")


# --------------------------------------------------------------------------
# Sample
# --------------------------------------------------------------------------

def load_sample() -> tuple[pd.DataFrame, tr.Rep, dict[str, Any], np.ndarray, dict[str, Any]]:
    """The `flows_block` sample of teacher_residence.py: employed persons 16+ in households with a US POWPUMA of work,
    classified on the POWPUMA geography at both ends.  Returns (persons, replicate weights, geo, selection mask over the
    persons, diagnostics)."""
    geo = tr.build_puma_geo()
    d, rep = tr.load_pums(geo)
    ok = (d.POWSP.between(1, 56) & d.POWPUMA.gt(1) & d.res_powpuma.notna()).to_numpy()      # as in flows_block
    cl = (d.res_pow_cls.isin(CLS) & d.wk_cls.isin(CLS)).to_numpy()
    sel = ok & cl
    diag: dict[str, Any] = dict(
        n_employed_in_households=int(len(d)), weighted_employed_millions=float(rep.w0.sum() / 1e6),
        n_us_powpuma=int(ok.sum()), weighted_us_powpuma_millions=float(rep.w0[ok].sum() / 1e6),
        n_unclassified_dropped=int((ok & ~cl).sum()), n_sample=int(sel.sum()),
        weighted_sample_millions=float(rep.w0[sel].sum() / 1e6),
        n_not_in_sample_no_workplace=int((~ok).sum()),
        weighted_not_in_sample_millions=float(rep.w0[~ok].sum() / 1e6))
    sub = d[sel].reset_index(drop=True)
    diag["work_at_home_share_of_sample"] = float(sub.PWGTP[sub.home].sum() / sub.PWGTP.sum())
    h = sub[sub.home]
    other_pow = ((h.POWSP != h.res_powsp) | (h.POWPUMA != h.res_powpuma))
    diag["n_work_at_home"] = int(len(h))
    diag["n_work_at_home_other_powpuma"] = int(other_pow.sum())
    diag["work_at_home_other_powpuma_share"] = float(h.PWGTP[other_pow].sum() / h.PWGTP.sum())
    k12 = sub.k12.to_numpy()
    assert np.array_equal(k12, sub.OCCP.between(2310, 2330).to_numpy()), "OCCP 2310-2330 differs from k12"
    diag["k12_share_of_sample"] = float(sub.PWGTP[k12].sum() / sub.PWGTP.sum())
    diag["occp_codes_in_2310_2330"] = sorted(int(c) for c in np.unique(sub.OCCP[sub.OCCP.between(2310, 2330)]))
    diag["national_check"] = d.attrs.get("national_check")
    return d, rep, geo, sel, diag


def group_masks(sub: pd.DataFrame) -> dict[str, np.ndarray]:
    k12 = sub.k12.to_numpy()
    return {"all": np.ones(len(sub), bool), "nonteacher": ~k12, "k12": k12, "pub_k12": k12 & sub.pub.to_numpy()}


def universe_masks(sub: pd.DataFrame) -> dict[str, np.ndarray]:
    rc = sub.res_cls.to_numpy()
    return {"all": np.ones(len(sub), bool), "away": ~sub.home.to_numpy(),
            "metro_puma": (rc == CITY) | (rc == SUB)}


# --------------------------------------------------------------------------
# Class-cut sensitivity
# --------------------------------------------------------------------------

def class_cuts(d: pd.DataFrame, rep: tr.Rep, geo: dict[str, Any], sel: np.ndarray,
               cuts: tuple[tuple[float, float], ...] = tr.CUTS) -> list[dict[str, Any]]:
    """The headline statistics of every group under the other cuts of the PUMA / POWPUMA classification
    (`teacher_residence.set_classes`); `d` and `geo` are restored to the original classification."""
    main = (geo["info"]["pc_cut"], geo["info"]["metro_cut"])
    sub0 = d[sel]
    gm = group_masks(sub0)
    rs = rep.sub(sel)
    out: list[dict[str, Any]] = []
    for pc_cut, metro_cut in cuts:
        tr.set_classes(geo, pc_cut, metro_cut)
        tr.attach_classes(d, geo)
        rc, wc = d.res_pow_cls.to_numpy()[sel], d.wk_cls.to_numpy()[sel]
        assert np.isin(rc, CLS).all() and np.isin(wc, CLS).all()
        specs = cell_specs(rc, wc, {f"{g}:all": m for g, m in gm.items()})
        for g, m in gm.items():
            specs[f"{g}:all|pop"] = (m.astype(float), np.ones(len(m)))
        R = tr.Ratios(rs, specs)
        F = registry(CLS)
        row: dict[str, Any] = dict(
            pc_cut=pc_cut, metro_cut=metro_cut, main=(pc_cut, metro_cut) == main,
            n_powpuma_by_class=dict(geo["info"]["n_powpuma_by_class"]), n_puma_by_class=dict(geo["info"]["n_puma_by_class"]),
            groups={})
        for g in gm:
            key = f"{g}:all"
            st = {name: R.derive(lambda gt, name=name: F[name](gt, key)) for name in HEADLINE}
            st["res"] = {a: R.derive(lambda gt, a=a: F[f"res.{a}"](gt, key)) for a in CLS}
            st["work"] = {b: R.derive(lambda gt, b=b: F[f"work.{b}"](gt, key)) for b in CLS}
            row["groups"][g] = st
        out.append(row)
        LOG.info("cut pc >= %.2f, metro >= %.2f: P(work suburb | live city) all %.4f, non-teachers %.4f, K-12 %.4f",
                 pc_cut, metro_cut, row["groups"]["all"]["p_work_suburb_given_live_city"]["value"],
                 row["groups"]["nonteacher"]["p_work_suburb_given_live_city"]["value"],
                 row["groups"]["k12"]["p_work_suburb_given_live_city"]["value"])
    tr.set_classes(geo, *main)
    tr.attach_classes(d, geo)
    return out


# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------

def reproduction_check(groups: dict[str, Any], n_sample_flows: int) -> dict[str, Any]:
    """The K-12 rows computed here against `flows.k12` / `flows.pub_k12` of teacher_residence.json."""
    if not REF_JSON.exists():
        return dict(checked=False, note=f"{REF_JSON.name} not found")
    J = json.loads(REF_JSON.read_text())["flows"]
    out: dict[str, Any] = dict(checked=True, reference=str(REF_JSON.relative_to(HERE)),
                               n_sample_reference=int(J["n_sample"]), n_sample_here=int(n_sample_flows))
    worst = 0.0
    for g in ("k12", "pub_k12"):
        mine, ref = groups[g]["universes"]["all"], J[g]
        pairs = {"gap_res": (mine["gap_res"], ref["gap_res"]), "gap_work": (mine["gap_work"], ref["gap_work"]),
                 "net_commute_gap": (mine["net_commute_gap"], ref["net_commute_gap"]),
                 "share_cross_class": (mine["share_cross_class"], ref["share_cross_class"]),
                 "share_metro_cross": (mine["share_cross_city_suburb"], ref["share_metro_cross"])}
        for b in ("city", "suburb"):
            for a in CLS:
                pairs[f"live_{a}_given_work_{b}"] = (mine["live_given_work"][b][a], ref["live_given_work"][f"work_{b}"][f"live_{a}"])
        rows = {}
        for name, (m, r) in pairs.items():
            dv, ds_ = abs(m["value"] - r["value"]), abs(m["se"] - r["se"])
            worst = max(worst, dv, ds_)
            rows[name] = dict(here=m["value"], here_se=m["se"], reference=r["value"], reference_se=r["se"],
                              abs_diff_value=dv, abs_diff_se=ds_)
        out[g] = rows
    out["max_abs_diff"] = worst
    if worst > 1e-8 or out["n_sample_reference"] != out["n_sample_here"]:
        raise SystemExit(f"reproduction of teacher_residence.json flows failed: max abs diff {worst:.3g}, "
                         f"n {out['n_sample_here']} vs {out['n_sample_reference']}")
    return out


def build(outdir: Path) -> dict[str, Any]:
    t0 = time.time()
    d, rep, geo, sel, diag = load_sample()
    LOG.info("sample: %d persons (%.2f m weighted) of %d employed persons in households; %d unclassified dropped",
             diag["n_sample"], diag["weighted_sample_millions"], diag["n_employed_in_households"],
             diag["n_unclassified_dropped"])
    sub = d[sel].reset_index(drop=True)
    rs = rep.sub(sel)
    rc, wc = sub.res_pow_cls.to_numpy(), sub.wk_cls.to_numpy()
    gm, um = group_masks(sub), universe_masks(sub)
    metro = np.isin(rc, METRO_ROWS)                                    # residence POWPUMA class city or suburb
    cmask = {f"{g}:{u}": mg & mu for g, mg in gm.items() for u, mu in um.items()}
    specs = cell_specs(rc, wc, cmask)
    one = np.ones(len(sub))
    for key, m in cmask.items():
        specs[f"{key}|pop"] = (m.astype(float), one)                   # the population's share of all workers
    R = tr.Ratios(rs, specs)
    cnt = dict(cmask)                                                  # head counts; metro_pow cells come from the `all` cells
    cnt.update({f"{g}:metro_pow": mg & metro for g, mg in gm.items()})
    M = np.column_stack([m.astype(float) for m in cnt.values()])
    est_n, rep_n = rs.run(lambda w: w @ M)
    se_n = tr.sdr_se(est_n, rep_n)
    tot_n = {k: dict(value=float(est_n[i]) / 1e6, se=float(se_n[i]) / 1e6) for i, k in enumerate(cnt)}
    LOG.info("ratios done in %.0fs", time.time() - t0)

    code = 3 * pd.Categorical(rc, categories=list(CLS)).codes + pd.Categorical(wc, categories=list(CLS)).codes
    groups: dict[str, Any] = {}
    for g in gm:
        gd: dict[str, Any] = dict(label=GROUP_LABEL[g], definition=GROUP_NOTE[g], universes={})
        for u in list(um) + ["metro_pow"]:
            key = f"{g}:{u}"
            blk = flow_block(R, f"{g}:all" if u == "metro_pow" else key, METRO_ROWS if u == "metro_pow" else CLS)
            blk["n_sample"] = int(cnt[key].sum())
            blk["workers_millions"] = tot_n[key]
            nc = np.bincount(code[cnt[key]], minlength=9).reshape(3, 3)
            blk["n_cells"] = {f"{a}->{b}": int(nc[i, j]) for i, a in enumerate(CLS) for j, b in enumerate(CLS)}
            if u == "metro_pow":                                       # share of the group's workers who live in a city / suburb POWPUMA
                blk["share_of_group"] = R.derive(lambda gt, g=g: sum(gt(cell_key(f"{g}:all", a, b)) for a in METRO_ROWS for b in CLS))
            else:
                blk["share_of_all_workers"] = R.get(f"{key}|pop")
            gd["universes"][u] = blk
        groups[g] = gd

    differences: dict[str, Any] = {}
    for u in ("all", "away"):
        differences[u] = {f"{g}_minus_nonteacher": diff_block(R, f"{g}:{u}", f"nonteacher:{u}") for g in ("k12", "pub_k12")}
    differences["metro_pow"] = {f"{g}_minus_nonteacher": diff_block(R, f"{g}:all", "nonteacher:all", METRO_ROWS)
                                for g in ("k12", "pub_k12")}
    LOG.info("blocks done in %.0fs", time.time() - t0)

    repro = reproduction_check(groups, diag["n_us_powpuma"])
    LOG.info("reproduction of teacher_residence.json flows: max abs diff %.2g", repro.get("max_abs_diff", float("nan")))
    sens = class_cuts(d, rep, geo, sel)
    LOG.info("class cuts done in %.0fs", time.time() - t0)

    keep = ("res", "work", "x", "p_work_suburb_given_live_city", "p_work_city_given_live_suburb", "share_cross_city_suburb",
            "share_cross_class", "block", "marginal")
    headline = {g: dict(all_workers={k: groups[g]["universes"]["all"][k] for k in keep},
                        metropolitan_residents={k: groups[g]["universes"]["metro_pow"][k] for k in keep})
                for g in groups}

    R_out = dict(
        meta=dict(generated=dt.datetime.now().isoformat(timespec="seconds"), script="data/spatial/commuting_all.py",
                  imports="data/spatial/teacher_residence.py (unchanged)", year=tr.YEAR,
                  pums_vintage="ACS 2018 1-year PUMS", replicate_weights=tr.N_REP,
                  se_formula="sqrt(4/80 * sum_r (x_r - x)^2) over the 80 successive-difference replicate weights; derived "
                             "statistics (ratios, differences, logs) are evaluated on every replicate",
                  pc_cut=geo["info"]["pc_cut"], metro_cut=geo["info"]["metro_cut"],
                  powpuma_classes=dict(n_powpuma_by_class=geo["info"]["n_powpuma_by_class"],
                                       n_puma_by_class=geo["info"]["n_puma_by_class"],
                                       pop_by_class_2010=geo["info"]["pop_by_class"]),
                  classes=list(CLS), universes=UNIV_LABEL, groups={g: GROUP_NOTE[g] for g in gm},
                  layout="groups.<group>.universes.<universe>: res[a], work[b], x['a->b'] (joint shares), work_given_live[a][b] = "
                         "P(work b | live a), live_given_work[b][a] = P(live a | work b), p_work_suburb_given_live_city, "
                         "p_work_city_given_live_suburb, share_cross_city_suburb, share_cross_class, gap_res, gap_work, "
                         "net_commute_gap, block{share_of_group, city_share_jobs, city_share_residents, net_inflow_to_city}, "
                         "marginal{city_share_residents, city_share_jobs}, n_sample, n_cells, workers_millions; each statistic "
                         "is {value, se}. differences.<universe>.<group>_minus_nonteacher.<statistic>; class_cuts[i].groups.<group>; "
                         "headline = the main statistics of groups.<group>.universes.all and .metro_pow",
                  methods=METHODS),
        sample=diag, headline=headline, groups=groups, differences=differences, class_cuts=sens,
        reproduction=repro)
    clean = ds.clean(R_out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "commuting_all.json").write_text(json.dumps(clean, indent=1))
    write_md(outdir / "commuting_all.md", clean)
    LOG.info("wrote commuting_all.{md,json} to %s (%.0fs)", outdir, time.time() - t0)
    return clean


METHODS = [
    "Data and sample: the ACS 2018 1-year PUMS cached by teacher_residence.py (employed persons 16+, ESR 1-2, in households; 80 "
    "replicate weights), restricted as in its `flows_block` to persons with a US place-of-work POWPUMA (POWSP 1-56, POWPUMA > 1). "
    "Workers at home carry their own POWPUMA as the place of work (checked: all of them do), so they count as working where they "
    "live; persons employed but not at work last week and workplaces abroad or in Puerto Rico are out.",
    "Classes: `teacher_residence.build_puma_geo()` at the main cut (0.5 / 0.5): a POWPUMA is `city` if at least half of its 2010 "
    "population lives in principal cities of OMB (Sept. 2018) metropolitan areas, `suburb` if at least half lives in metropolitan "
    "counties and it is not `city`, `other` otherwise. Residence class = class of the POWPUMA that contains the residence PUMA; "
    "workplace class = class of the place-of-work POWPUMA (one geography at both ends).",
    "Groups: all workers; non-teachers = everyone except OCCP 2310-2330 (the `k12` flag of teacher_residence.py). K-12 and public "
    "K-12 (COW 3-4) teachers are recomputed by the same code as a check and reproduce `flows.k12` / `flows.pub_k12` of "
    "teacher_residence.json.",
    "Estimator: person-weighted cell shares; every probability is a ratio of cells. Standard errors from the 80 ACS replicate "
    "weights, SE = sqrt(4/80 sum_r (x_r - x)^2), with derived statistics and between-group differences evaluated on every "
    "replicate (covariances respected). Sampling error only.",
    "Definitions: `block` = workers who live AND work in a city or suburb POWPUMA (the closed two-location system); crossing share "
    "= (city->suburb + suburb->city) / block (`share_metro_cross` of teacher_residence.py); `marginal` = workers who live (work) in "
    "a city or suburb POWPUMA, whatever the other end. Caveat: a POWPUMA is a whole county or group of counties, so city-suburb "
    "commuting inside one county is invisible and the crossing probabilities are lower bounds; workers at home count as "
    "non-crossers (see `away`); the classes depend on the cuts (`class_cuts`).",
]


# --------------------------------------------------------------------------
# Report
# --------------------------------------------------------------------------

def p1(s: dict[str, float]) -> str:
    return f"{100 * s['value']:.1f}% ({100 * s['se']:.2f})"


def p2(s: dict[str, float]) -> str:
    return f"{100 * s['value']:.2f}% ({100 * s['se']:.2f})"


def f3(s: dict[str, float]) -> str:
    return f"{s['value']:.3f} ({s['se']:.3f})"


def pp(s: dict[str, float]) -> str:
    """a difference of shares in percentage points"""
    return f"{100 * s['value']:+.2f} pp ({100 * s['se']:.2f})"


def write_md(path: Path, R: dict[str, Any]) -> None:
    L: list[str] = []
    a = L.append
    G, Dif, S, M = R["groups"], R["differences"], R["sample"], R["meta"]
    order = ["all", "nonteacher", "k12", "pub_k12"]
    U = {g: G[g]["universes"] for g in order}
    gl = {g: GROUP_LABEL[g] for g in order}

    a("# Residence against workplace of all workers and of non-teachers (D1b)\n")
    a(f"Generated {M['generated']} by `data/spatial/commuting_all.py`, which imports `teacher_residence.py` (unchanged) and builds, "
      "for all workers and for non-teachers, the table of `teacher_residence.md` section \"Residence against workplace of K-12 "
      "teachers on one geography\" (JSON key `flows`): ACS 2018 1-year PUMS, employed persons 16+ with a US POWPUMA of work, residence "
      "and workplace classified on the same POWPUMA classes (city / suburb / other), replicate-weight standard errors in parentheses "
      "(percentage points for shares). Machine-readable: `commuting_all.json` (`groups.<group>.universes.<universe>`).\n")

    # ---- headline
    a("## Headline\n")
    a(f"Sample: {S['n_sample']:,} persons, {S['weighted_sample_millions']:.1f} million weighted; K-12 teachers (OCCP 2310-2330) are "
      f"{100 * S['k12_share_of_sample']:.2f}% of them.\n")
    a("| group | live: city | suburb | other | work: city | suburb | other |")
    a("|---|---|---|---|---|---|---|")
    for g in ("all", "nonteacher"):
        u = U[g]["all"]
        a(f"| **{gl[g]}** | " + " | ".join(p1(u["res"][c]) for c in CLS) + " | " + " | ".join(p1(u["work"][c]) for c in CLS) + " |")
    a("")
    a("| group | P(work in suburb \\| live in city) | P(work in city \\| live in suburb) | city<->suburb crossing share |")
    a("|---|---|---|---|")
    for g in order:
        u = U[g]["all"]
        a(f"| {'**' + gl[g] + '**' if g in ('all', 'nonteacher') else gl[g]} | {p2(u['p_work_suburb_given_live_city'])} | "
          f"{p2(u['p_work_city_given_live_suburb'])} | {p2(u['share_cross_city_suburb'])} |")
    a("")
    a("| group | block: city share of jobs | city share of residents | net inflow into the city | marginal: city share of jobs | "
      "city share of residents |")
    a("|---|---|---|---|---|---|")
    for g in order:
        u = U[g]["all"]
        a(f"| {'**' + gl[g] + '**' if g in ('all', 'nonteacher') else gl[g]} | {p2(u['block']['city_share_jobs'])} | "
          f"{p2(u['block']['city_share_residents'])} | {pp(u['block']['net_inflow_to_city'])} | {p2(u['marginal']['city_share_jobs'])} | "
          f"{p2(u['marginal']['city_share_residents'])} |")
    a("")
    a("The crossing share is the share of workers who live and work in a city or suburb POWPUMA and live in one and work in the other "
      "(`share_metro_cross` of teacher_residence.py). `block` shares are among those same workers (a closed city-suburb economy); "
      "`marginal` shares are among workers who live (work) in a city or suburb POWPUMA whatever the other end. The K-12 rows are "
      "recomputed by the same code and reproduce `teacher_residence.json` (section 6).\n")
    ua, un = U["all"]["all"], U["nonteacher"]["all"]
    dk = Dif["all"]["k12_minus_nonteacher"]

    def rel(s_: dict[str, float]) -> str:
        return "about as often" if abs(s_["value"]) < 2 * s_["se"] else ("less often" if s_["value"] < 0 else "more often")

    a(f"- **All workers.** {p1(ua['p_work_suburb_given_live_city'])} of city residents work in a suburban POWPUMA and "
      f"{p1(ua['p_work_city_given_live_suburb'])} of suburban residents work in a city POWPUMA; {p1(ua['share_cross_city_suburb'])} of the "
      f"workers who live and work in the city or suburb cross between them. The city holds {p1(ua['block']['city_share_jobs'])} of those "
      f"workers' jobs and {p1(ua['block']['city_share_residents'])} of their residents (net inflow {pp(ua['block']['net_inflow_to_city'])}).")
    a(f"- **Non-teachers** are almost the same as all workers (teachers are {100 * S['k12_share_of_sample']:.1f}% of them): "
      f"{p1(un['p_work_suburb_given_live_city'])}, {p1(un['p_work_city_given_live_suburb'])}, crossing share {p1(un['share_cross_city_suburb'])}.")
    a(f"- **K-12 teachers against non-teachers** (difference, SE): P(work in suburb | live in city) {pp(dk['p_work_suburb_given_live_city'])}, "
      f"{rel(dk['p_work_suburb_given_live_city'])}; P(work in city | live in suburb) {pp(dk['p_work_city_given_live_suburb'])}, "
      f"{rel(dk['p_work_city_given_live_suburb'])}; crossing share {pp(dk['share_cross_city_suburb'])}; net commuting gap "
      f"{f3(dk['net_commute_gap'])} log points (section 4).\n")

    # ---- methods
    a("## Methods note\n")
    for t in M["methods"]:
        a(f"- {t}")
    a("")

    # ---- 1 shares
    a("## 1. Where workers live and where they work\n")
    a("Shares of each group's workers (the three residence shares and the three workplace shares each sum to 100%). The last four "
      "columns are the log gaps of teacher_residence.py, log(suburb / city) by residence, by workplace and their difference (the net "
      "commuting gap), and the share working in a different class from the one they live in.\n")
    a("| group | n | workers (millions) | live: city | suburb | other | work: city | suburb | other | gap by residence | gap by workplace | "
      "net commuting gap | works in the other class |")
    a("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for g in order:
        u = U[g]["all"]
        a(f"| {gl[g]} | {u['n_sample']:,} | {u['workers_millions']['value']:.2f} | " + " | ".join(p1(u["res"][c]) for c in CLS) + " | "
          + " | ".join(p1(u["work"][c]) for c in CLS) + f" | {f3(u['gap_res'])} | {f3(u['gap_work'])} | {f3(u['net_commute_gap'])} | "
          f"{p1(u['share_cross_class'])} |")
    a("")
    a("City and suburb shares among workers who live (work) in one of the two, whatever the other end (`marginal`), and among workers who "
      "live and work in one of the two (`block`):\n")
    a("| group | marginal: city share of residents | city share of jobs | block: share of the group | city share of residents | city share of jobs |")
    a("|---|---|---|---|---|---|")
    for g in order:
        u = U[g]["all"]
        a(f"| {gl[g]} | {p2(u['marginal']['city_share_residents'])} | {p2(u['marginal']['city_share_jobs'])} | "
          f"{p1(u['block']['share_of_group'])} | {p2(u['block']['city_share_residents'])} | {p2(u['block']['city_share_jobs'])} |")
    a("")

    # ---- 2 matrices
    a("## 2. The 3x3 residence x workplace matrices\n")
    a("Joint shares of the group's workers; the last column is the residence share and the last row the workplace share. Below each, the "
      "row-conditional probabilities P(work in b | live in a) and the column-conditional P(live in a | work in b).\n")
    for g in order:
        u = U[g]["all"]
        a(f"**{gl[g]}** ({u['n_sample']:,} persons, {u['workers_millions']['value']:.1f} million; the smallest cell holds "
          f"{min(u['n_cells'].values()):,} persons)\n")
        a("| live \\ work | city | suburb | other | all |")
        a("|---|---|---|---|---|")
        for r_ in CLS:
            a(f"| {r_} | " + " | ".join(p2(u["x"][f"{r_}->{c}"]) for c in CLS) + f" | {p2(u['res'][r_])} |")
        a("| all | " + " | ".join(p2(u["work"][c]) for c in CLS) + " | 100% |")
        a("")
        a("| P(work in ... \\| live in ...) | city | suburb | other |")
        a("|---|---|---|---|")
        for r_ in CLS:
            a(f"| live in {r_} | " + " | ".join(p2(u["work_given_live"][r_][c]) for c in CLS) + " |")
        a("")
        a("| P(live in ... \\| work in ...) | live: city | suburb | other |")
        a("|---|---|---|---|")
        for c in CLS:
            a(f"| work in {c} | " + " | ".join(p2(u["live_given_work"][c][r_]) for r_ in CLS) + " |")
        a("")

    # ---- 3 metro
    a("## 3. Residents of metropolitan areas only\n")
    a("The script already distinguishes metropolitan residents (`metro` = residence class city or suburb, used in its commute block on PUMA "
      "classes). Two versions: (a) on the geography of this table, residence POWPUMA class city or suburb (rows `city` and `suburb` of the "
      "matrix; the conditional probabilities, the crossing share and the block shares are then identical to section 2 by construction, only "
      "the margins change); (b) residence PUMA class city or suburb (the PUMA-level `metro` flag), cells on the POWPUMA classes.\n")
    a("| group | version | n | workers (millions) | live: city | suburb | other | work: city | suburb | other (outside the two classes) |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for g in ("all", "nonteacher", "k12", "pub_k12"):
        for ver, lab in (("metro_pow", "(a) POWPUMA"), ("metro_puma", "(b) PUMA flag")):
            u = U[g][ver]
            oth = p1(u["res"]["other"]) if "other" in u["res"] else "-"
            a(f"| {gl[g]} | {lab} | {u['n_sample']:,} | {u['workers_millions']['value']:.2f} | {p1(u['res'][CITY])} | {p1(u['res'][SUB])} | "
              f"{oth} | {p1(u['work'][CITY])} | {p1(u['work'][SUB])} | {p1(u['work']['other'])} |")
    a("")
    a("| group | version | P(work in suburb \\| live in city) | P(work in city \\| live in suburb) | city<->suburb crossing share | "
      "city share of jobs (block) | city share of residents (block) | works in the other class |")
    a("|---|---|---|---|---|---|---|---|")
    for g in ("all", "nonteacher", "k12", "pub_k12"):
        for ver, lab in (("metro_pow", "(a) POWPUMA"), ("metro_puma", "(b) PUMA flag")):
            u = U[g][ver]
            a(f"| {gl[g]} | {lab} | {p2(u['p_work_suburb_given_live_city'])} | {p2(u['p_work_city_given_live_suburb'])} | "
              f"{p2(u['share_cross_city_suburb'])} | {p2(u['block']['city_share_jobs'])} | {p2(u['block']['city_share_residents'])} | "
              f"{p1(u['share_cross_class'])} |")
    a("")
    a("Metropolitan-resident matrices, version (a): joint shares of the group's metropolitan residents.\n")
    for g in ("all", "nonteacher"):
        u = U[g]["metro_pow"]
        a(f"**{gl[g]}**\n")
        a("| live \\ work | city | suburb | other | all |")
        a("|---|---|---|---|---|")
        for r_ in METRO_ROWS:
            a(f"| {r_} | " + " | ".join(p2(u["x"][f"{r_}->{c}"]) for c in CLS) + f" | {p2(u['res'][r_])} |")
        a("| all | " + " | ".join(p2(u["work"][c]) for c in CLS) + " | 100% |")
        a("")

    # ---- 4 differences
    a("## 4. K-12 teachers minus non-teachers\n")
    a("Differences of the statistics between the groups (same replicate weights, so the standard error respects the covariance). "
      "Percentage points for the shares, log points for the gaps.\n")
    a("| statistic | K-12 minus non-teachers | public K-12 minus non-teachers | K-12 minus non-teachers, workers away from home |")
    a("|---|---|---|---|")
    names = {"p_work_suburb_given_live_city": "P(work in suburb \\| live in city)", "p_work_city_given_live_suburb": "P(work in city \\| live in suburb)",
             "share_cross_city_suburb": "city<->suburb crossing share", "share_cross_class": "works in the other class",
             "block.city_share_jobs": "city share of jobs (block)", "block.city_share_residents": "city share of residents (block)",
             "block.net_inflow_to_city": "net inflow into the city (block)", "marginal.city_share_residents": "city share of residents (marginal)",
             "marginal.city_share_jobs": "city share of jobs (marginal)", "gap_res": "gap by residence (log)", "gap_work": "gap by workplace (log)",
             "net_commute_gap": "net commuting gap (log)"}
    for k, lab in names.items():
        fmt = f3 if k in ("gap_res", "gap_work", "net_commute_gap") else pp
        a(f"| {lab} | {fmt(Dif['all']['k12_minus_nonteacher'][k])} | {fmt(Dif['all']['pub_k12_minus_nonteacher'][k])} | "
          f"{fmt(Dif['away']['k12_minus_nonteacher'][k])} |")
    a("")

    # ---- 5 robustness
    a("## 5. Robustness\n")
    a("**Workers away from home** (JWTR != 11: work-at-home workers carry their own area as the workplace and never cross; dropping them "
      "raises the crossing probabilities).\n")
    a("| group | n | P(work in suburb \\| live in city) | P(work in city \\| live in suburb) | city<->suburb crossing share | works in the other class |")
    a("|---|---|---|---|---|---|")
    for g in order:
        u = U[g]["away"]
        a(f"| {gl[g]} | {u['n_sample']:,} | {p2(u['p_work_suburb_given_live_city'])} | {p2(u['p_work_city_given_live_suburb'])} | "
          f"{p2(u['share_cross_city_suburb'])} | {p1(u['share_cross_class'])} |")
    a("")
    a("**Class cuts** (`teacher_residence.CUTS`: a POWPUMA is `city` if its principal-city population share is at least pc, `suburb` if its "
      "metropolitan share is at least metro and it is not `city`; the main line is 0.50 / 0.50).\n")
    a("| pc cut | metro cut | POWPUMAs city / suburb / other | group | P(work in suburb \\| live in city) | P(work in city \\| live in suburb) | "
      "crossing share | city share of jobs (block) | city share of residents (block) |")
    a("|---|---|---|---|---|---|---|---|---|")
    for r_ in R["class_cuts"]:
        n3 = " / ".join(str(r_["n_powpuma_by_class"].get(c, 0)) for c in CLS)
        for g in ("all", "nonteacher", "k12"):
            s_ = r_["groups"][g]
            a(f"| {r_['pc_cut']:.2f} | {r_['metro_cut']:.2f}{' (main)' if r_['main'] else ''} | {n3} | {gl[g]} | "
              f"{p2(s_['p_work_suburb_given_live_city'])} | {p2(s_['p_work_city_given_live_suburb'])} | {p2(s_['share_cross_city_suburb'])} | "
              f"{p2(s_['block.city_share_jobs'])} | {p2(s_['block.city_share_residents'])} |")
    a("")

    # ---- 6 checks
    rp = R["reproduction"]
    a("## 6. Checks\n")
    if rp.get("checked"):
        a(f"- **Reproduction of `teacher_residence.json`.** The K-12 and public K-12 rows computed here (gap by residence, gap by workplace, net "
          f"commuting gap, share working in the other class, city<->suburb crossing share, and the six `live_given_work` probabilities, values and "
          f"standard errors) equal `flows.k12` and `flows.pub_k12` of the existing file: largest absolute difference {rp['max_abs_diff']:.1e}; "
          f"sample {rp['n_sample_here']:,} persons here and {rp['n_sample_reference']:,} there.")
    else:
        a(f"- Reproduction check skipped: {rp.get('note')}.")
    a(f"- **Sample.** {S['n_employed_in_households']:,} employed persons 16+ in households ({S['weighted_employed_millions']:.1f} million weighted); "
      f"{S['n_us_powpuma']:,} ({S['weighted_us_powpuma_millions']:.1f} million) have a US place-of-work POWPUMA; {S['n_unclassified_dropped']} "
      f"of them could not be classified at one end and are dropped. The {S['n_not_in_sample_no_workplace']:,} persons without a US POWPUMA of work "
      f"({S['weighted_not_in_sample_millions']:.1f} million weighted) are employed-but-not-at-work persons and workers with a workplace abroad "
      f"or in Puerto Rico.")
    a(f"- **Work at home.** {100 * S['work_at_home_share_of_sample']:.1f}% of the sample works at home ({S['n_work_at_home']:,} persons); "
      f"{S['n_work_at_home_other_powpuma']} of them report a place of work in a POWPUMA other than the residence's.")
    a(f"- **Teacher definition.** OCCP codes found in 2310-2330: {S['occp_codes_in_2310_2330']} (equal to `k12`, asserted); the non-teacher group "
      f"keeps OCCP 2300 (preschool and kindergarten teachers) and the other education occupations.")
    nat = S.get("national_check") or {}
    if nat:
        a(f"- **Pipeline check** (PUMS loader of teacher_residence.py): mean travel time {nat['mean_travel_time_all_persons']:.1f} minutes against "
          f"{nat['published_mean_travel_time']} published, work-at-home share {100 * nat['work_at_home_share_all_persons']:.1f}% against "
          f"{100 * nat['published_work_at_home_share']:.1f}%.")
    a("- Standard errors are sampling errors of the PUMS only; the classification of the areas is taken as given.")
    path.write_text("\n".join(L) + "\n")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    sub.add_parser("selftest", help="unit check of the statistics registry (no data)")
    r = sub.add_parser("render", help="rewrite the .md report from the existing .json")
    r.add_argument("--outdir", default=None, help="default: <here>/estimates/settings")
    b = sub.add_parser("build", help="compute the tables and write estimates/settings/commuting_all.{md,json}")
    b.add_argument("--outdir", default=None, help="default: <here>/estimates/settings")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("data_estimate").setLevel(logging.WARNING)
    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "render":
        outdir = Path(args.outdir).resolve() if args.outdir else HERE / "estimates" / "settings"
        write_md(outdir / "commuting_all.md", json.loads((outdir / "commuting_all.json").read_text()))
        LOG.info("rewrote %s", outdir / "commuting_all.md")
    else:
        outdir = Path(args.outdir).resolve() if getattr(args, "outdir", None) else HERE / "estimates" / "settings"
        build(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
