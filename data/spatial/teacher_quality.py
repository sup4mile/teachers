#!/usr/bin/env python3
"""How much teacher quality is there in SEDA, and what must the model supply?

Companion to `data_estimate.py`.  That script builds the spatial moments; this
one asks the narrower question raised by the data-gap list at the end of its
report -- *nothing here measures teacher ability* -- and works out exactly how
far the public data can be pushed before the answer has to come from outside.

    uv run teacher_quality.py            # writes estimates/teacher_quality.txt

Every number in `notes/teacher_quality_notes.md` is printed by this script.
It reuses `data_estimate.assemble`, so the sample, the outlier screen and the
commuting-zone definition are identical to the main report's.

Six blocks:

  A  reliability   is the district growth measure precise enough to use?
  B  composition   how much of it is just who enrolls?
  C  inputs        do observable teacher inputs explain any of it?
  D  gaps          the within-CZ growth gap, raw and net of composition
  E  residual      what the model must then attribute to teacher quality
  F  aggregation   why no district-level series can identify teacher VA

Written with help from Claude Code.
"""
from __future__ import annotations

import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SPATIAL = Path(__file__).resolve().parent
sys.path.insert(0, str(SPATIAL))
os.chdir(SPATIAL)
import data_estimate as de           # noqa: E402

YEAR = 2018
GRADE_SPAN = de.SEDA_GRADE_SPAN      # 5 grades, 3-8

# Teacher value-added yardsticks the residual is measured against.  All are
# SDs of student achievement gained in ONE year from ONE teacher.
VA_SD_MATH_CFR = 0.116               # Chetty-Friedman-Rockoff I, elementary, shrunk
VA_SD_ELA_CFR = 0.080
VA_SD_MATH_TRUE = 0.163              # CFR I, underlying (unshrunk) elementary
VA_WITHIN_SCHOOL_SHARE = 0.85        # CFR I: 85% of VA variance is within school
ETG_MATH = 0.024                     # Isenberg et al. (2013), FRL gap, math
ETG_ELA = 0.034                      # ... and ELA

LOG = logging.getLogger("teacher_quality")
_out: list[str] = []


def say(*a: object) -> None:
    s = " ".join(str(x) for x in a)
    print(s)
    _out.append(s)


def wmean(x: pd.Series, w: pd.Series) -> float:
    m = x.notna() & w.notna() & (w > 0)
    return float((x[m] * w[m]).sum() / w[m].sum()) if m.any() else np.nan


def load_panel() -> pd.DataFrame:
    """The main report's 2018 panel, plus the SEDA columns it does not read.

    `data_estimate` keeps one learning rate (cs / ol).  Judging whether that
    number can carry the weight the calibration puts on it needs the standard
    errors, the shrunken estimate, and the grade-equivalent scale as well.
    """
    store = de.RawStore(SPATIAL / "raw")
    panel = de.assemble(store, [YEAR], "czone")
    panel = panel[panel.year == YEAR]

    cs = pd.read_parquet(
        SPATIAL / "raw/seda_geodist_pool_cs_6_0.parquet",
        columns=["sedalea", "subgroup", "gap", "cs_mn_avg_ol", "cs_mn_lrn_ol",
                 "cs_mn_avg_ol_se", "cs_mn_lrn_ol_se", "cs_mn_lrn_ol_se_adj",
                 "cs_mn_lrn_eb"])
    gys = pd.read_parquet(
        SPATIAL / "raw/seda_geodist_pool_gys_6_0.parquet",
        columns=["sedalea", "subgroup", "gap", "gys_mn_lrn_ol",
                 "gys_mn_lrn_ol_se"])
    keep = lambda d: d[(d.gap == 0) & (d.subgroup == "all")].drop(
        columns=["subgroup", "gap"])
    sd = keep(cs).merge(keep(gys), on="sedalea", how="outer")
    sd["leaid"] = de.id7(sd.sedalea)
    sd = sd.drop_duplicates("leaid").drop(columns=["sedalea"])

    panel = panel.merge(sd, on="leaid", how="left")
    panel["learn_gain"] = panel.cs_mn_lrn_ol * GRADE_SPAN
    panel["se_gain"] = panel.cs_mn_lrn_ol_se * GRADE_SPAN
    panel["log_Q"] = np.log(panel.gys_mn_lrn_ol.where(panel.gys_mn_lrn_ol > 0))
    return panel


def demeaned(d: pd.DataFrame, var: str, w: pd.Series, fe: str | None
             ) -> pd.Series:
    if fe is None:
        return d[var] - (d[var] * w).sum() / w.sum()
    num = (d[var] * w).groupby(d[fe]).transform("sum")
    den = w.groupby(d[fe]).transform("sum")
    return d[var] - num / den


def reliability(d: pd.DataFrame, var: str, se: str, w: str | None,
                fe: str | None) -> tuple[float, float, float, int]:
    """Share of the observed dispersion that is signal rather than sampling error.

    SEDA publishes a standard error for every estimate, so the noise term is
    known rather than assumed: var(signal) = var(estimate) - mean(se^2).
    """
    m = d[[var, se]].notna().all(axis=1)
    if w:
        m &= d[w].notna() & (d[w] > 0)
    if fe:
        m &= d[fe].notna()
    d = d[m]
    ww = d[w].astype(float) if w else pd.Series(1.0, index=d.index)
    x = demeaned(d, var, ww, fe)
    tot = float((ww * x ** 2).sum() / ww.sum())
    err = float((ww * d[se] ** 2).sum() / ww.sum())
    return tot, err, (max(tot - err, 0) / tot if tot > 0 else np.nan), len(d)


def block_a(panel: pd.DataFrame) -> float:
    say("\n" + "=" * 78)
    say("A. RELIABILITY -- is the district growth measure precise enough to use?")
    say("=" * 78)
    say("   SEDA suppresses any estimate whose OLS reliability is below 0.7, so")
    say("   this is a selected sample: the imprecise districts are already gone.")
    say(f"\n   {'measure':34s} {'weight':11s} {'':10s} {'sd':>7s} {'noise sd':>9s}"
        f" {'reliab.':>8s} {'n':>7s}")
    rows = [("cs_mn_avg_ol", "cs_mn_avg_ol_se", "level  mn_avg (cs)"),
            ("cs_mn_lrn_ol", "cs_mn_lrn_ol_se", "growth mn_lrn (cs)"),
            ("cs_mn_lrn_ol", "cs_mn_lrn_ol_se_adj", "growth mn_lrn (cs, adj. se)"),
            ("gys_mn_lrn_ol", "gys_mn_lrn_ol_se", "growth mn_lrn (gys, grade eq.)")]
    rel_within = np.nan
    for var, se, lab in rows:
        for wgt in ("pupils", None):
            for fe in (None, "czone"):
                tot, err, rel, n = reliability(panel, var, se, wgt, fe)
                say(f"   {lab:34s} {'enrol-wtd' if wgt else 'unweighted':11s} "
                    f"{'within CZ' if fe else 'national':10s} {np.sqrt(tot):7.4f} "
                    f"{np.sqrt(err):9.4f} {rel:8.3f} {n:7d}")
                if var == "cs_mn_lrn_ol" and se.endswith("_se") and wgt and fe:
                    rel_within = rel
    say("\n   Read: sampling noise is not what is keeping the quality gap small.")
    return rel_within


def block_b(panel: pd.DataFrame) -> None:
    say("\n" + "=" * 78)
    say("B. COMPOSITION -- how much of each measure is just who enrolls?")
    say("=" * 78)
    say("   within-CZ WLS, enrolment weights, CZ-clustered SEs\n")
    say(f"   {'outcome':16s} {'on':26s} {'R2 within':>9s} {'b(SES)':>9s}"
        f" {'se':>8s} {'n':>7s}")
    for y in ("cs_mn_avg_ol", "cs_mn_lrn_ol", "gys_mn_lrn_ol"):
        for xs in (["seda_ses"], ["seda_ses", "pov_rate_5_17"]):
            d = panel.dropna(subset=[y] + xs + ["pupils", "czone"])
            r = de.ols_fe(d, y, xs, fe="czone", w="pupils", cluster="czone")
            if r:
                say(f"   {y:16s} {'+'.join(xs):26s} {r['r2_within']:9.3f} "
                    f"{r['coef']['seda_ses']:+9.4f} {r['se']['seda_ses']:8.4f} "
                    f"{r['n']:7d}")
    say("\n   Read: the level is composition; the growth rate mostly is not.")
    say("   That is the case for using the learning rate as the Q_l analogue --")
    say("   and the reason it cannot also stand in for teacher quality.")


def block_c(panel: pd.DataFrame) -> None:
    say("\n" + "=" * 78)
    say("C. INPUTS -- do observable teacher inputs explain district growth?")
    say("=" * 78)
    base = ["seda_ses", "pov_rate_5_17"]
    tin = ["log_teachers_pp", "log_salary_real", "log_exp_pp"]
    for y in ("cs_mn_lrn_ol", "cs_mn_avg_ol"):
        say(f"\n   outcome: {y}")
        for xs, lab in ((tin, "teacher inputs only"), (base, "SES only"),
                        (base + tin, "both")):
            d = panel.dropna(subset=[y] + xs + ["pupils", "czone"])
            r = de.ols_fe(d, y, xs, fe="czone", w="pupils", cluster="czone")
            if not r:
                continue
            coefs = "  ".join(f"{k.replace('log_', '')}={r['coef'][k]:+.4f}"
                              f"({r['se'][k]:.4f})" for k in tin if k in r["coef"])
            say(f"     {lab:20s} R2w={r['r2_within']:.4f} n={r['n']:6d}  {coefs}")
    say("\n   Read: teacher inputs the CCD and F-33 can see move district growth")
    say("   hardly at all.  Whatever separates districts, it is not measured by")
    say("   headcount, salary or spending -- which is the Chetty-Friedman-Rockoff")
    say("   finding one aggregation level up.")


def block_d(panel: pd.DataFrame) -> dict[str, pd.DataFrame]:
    say("\n" + "=" * 78)
    say("D. GAPS -- the within-CZ growth gap, raw and net of composition")
    say("=" * 78)
    VARS = ["cs_mn_lrn_ol", "gys_mn_lrn_ol", "cs_mn_avg_ol", "log_teachers_pp"]
    res: dict[str, pd.DataFrame] = {}
    say(f"\n   {'scheme':8s} {'variable':16s} {'loc 1':>9s} {'loc 2':>9s}"
        f" {'gap':>9s} {'se':>8s} {'t':>6s} {'CZs':>5s}")
    for scheme in ("locale", "revenue", "ses"):
        part = de.partition(panel, scheme)
        part = part[part.year == YEAR]
        rows = []
        for cz, g in part.groupby("czone"):
            rec = {"czone": cz, "pupils": g.pupils.sum()}
            for v in VARS:
                for l in (1, 2):
                    gl = g[g.loc_id == l]
                    rec[f"{v}_{l}"] = wmean(gl[v], gl.pupils)
                rec[f"gap_{v}"] = rec[f"{v}_2"] - rec[f"{v}_1"]
            rows.append(rec)
        G = pd.DataFrame(rows)
        res[scheme] = G
        for v in VARS:
            m = G[f"gap_{v}"].notna()
            x, w = G.loc[m, f"gap_{v}"], G.loc[m, "pupils"]
            mu = float((x * w).sum() / w.sum())
            ws = w / w.sum()
            n = len(x)
            se = float(np.sqrt(((ws ** 2) * (x - mu) ** 2).sum() * n / max(n - 1, 1)))
            say(f"   {scheme:8s} {v:16s} {wmean(G.loc[m, f'{v}_1'], w):9.4f} "
                f"{wmean(G.loc[m, f'{v}_2'], w):9.4f} {mu:+9.4f} {se:8.4f} "
                f"{mu / se:6.2f} {n:5d}")
        say("")

    say("   The same gap as a coefficient on a location-2 dummy, so that student")
    say("   composition can be controlled for (CZ FE, enrolment-weighted):\n")
    say(f"   {'scheme':8s} {'specification':22s} {'gap':>9s} {'se':>8s}"
        f" {'t':>6s} {'n':>7s}")
    for scheme in ("locale", "revenue", "ses"):
        p = de.partition(panel, scheme)
        p = p[p.year == YEAR].assign(loc2=lambda d: (d.loc_id == 2).astype(float))
        for ctl, lab in (([], "raw gap"),
                         (["seda_ses"], "+ SES"),
                         (["seda_ses", "pov_rate_5_17"], "+ SES + poverty"),
                         (["seda_ses", "pov_rate_5_17", "log_teachers_pp",
                           "log_exp_pp"], "+ SES + inputs")):
            d = p.dropna(subset=["cs_mn_lrn_ol", "loc2", "pupils", "czone"] + ctl)
            r = de.ols_fe(d, "cs_mn_lrn_ol", ["loc2"] + ctl, fe="czone",
                          w="pupils", cluster="czone")
            if r:
                b, s = r["coef"]["loc2"], r["se"]["loc2"]
                say(f"   {scheme:8s} {lab:22s} {b:+9.4f} {s:8.4f} {b / s:6.2f}"
                    f" {r['n']:7d}")
        say("")
    say("   Read: on the paper's preferred (city/suburb) split the growth")
    say("   advantage of the good location does not survive controlling for who")
    say("   enrols -- it changes sign.  Q_2/Q_1 is 1.03 raw and about 1.00 net.")
    return res


def block_e(res: dict[str, pd.DataFrame], panel: pd.DataFrame) -> None:
    say("\n" + "=" * 78)
    say("E. RESIDUAL -- what the model must then attribute to teacher quality")
    say("=" * 78)
    say("   Q_l = (2 H_T,l / M_l)^sigma  and  H_T,l = N_T,l * E_l[h^beta], so")
    say("       dlog E[h^beta] = dlog Q / sigma - dlog(teachers per pupil).")
    say("   dlog Q comes from the grade-equivalent learning rate, which is a")
    say("   genuine ratio; the CS learning rate is a deviation around the")
    say("   national average and has no level to take a ratio of.\n")
    say(f"   {'scheme':8s} {'dlogQ':>9s} {'dlog(N_T/M)':>12s}"
        + "".join(f" {'sigma=' + str(s):>11s}" for s in (0.3, 0.5, 0.7)))
    implied = {}
    for scheme, G in res.items():
        m = G[["gap_gys_mn_lrn_ol", "gap_log_teachers_pp"]].notna().all(axis=1)
        w = G.loc[m, "pupils"]
        q1 = wmean(G.loc[m, "gys_mn_lrn_ol_1"], w)
        q2 = wmean(G.loc[m, "gys_mn_lrn_ol_2"], w)
        dlogQ = float(np.log(q2 / q1))
        dNT = wmean(G.loc[m, "gap_log_teachers_pp"], w)
        implied[scheme] = {s: dlogQ / s - dNT for s in (0.3, 0.5, 0.7)}
        say(f"   {scheme:8s} {dlogQ:+9.4f} {dNT:+12.4f}"
            + "".join(f" {dlogQ / s - dNT:+11.4f}" for s in (0.3, 0.5, 0.7)))
    say("   (columns: the implied within-CZ log gap in beta-weighted teacher quality)")

    say("\n   Is that gap big or small?  Convert the growth gap into the units the")
    say("   teacher value-added literature reports -- SDs of student achievement")
    say("   gained from one teacher in one year -- and compare:\n")
    say(f"   {'scheme':8s} {'growth gap/grade':>17s} {'as share of VA sd':>19s}"
        f" {'vs FRL access gap':>19s}")
    for scheme, G in res.items():
        w = G.pupils
        gap = wmean(G.gap_cs_mn_lrn_ol, w)
        say(f"   {scheme:8s} {gap:+17.4f} {gap / VA_SD_MATH_CFR:18.2f}x"
            f" {gap / ETG_MATH:18.2f}x")
    say(f"\n   VA sd = {VA_SD_MATH_CFR} (CFR I, elementary math, shrunk);")
    say(f"   FRL access gap = {ETG_MATH} SD math (Isenberg et al. 2013, 29 districts).")
    say("   Read: even attributing the ENTIRE within-CZ growth gap to teachers,")
    say("   the district-mean teacher-quality gap the model needs is a few percent")
    say("   of a teacher-VA standard deviation on the city/suburb split.")


def block_f(panel: pd.DataFrame, rel_within: float) -> None:
    say("\n" + "=" * 78)
    say("F. AGGREGATION -- why no district series can identify teacher VA")
    say("=" * 78)
    d = panel.dropna(subset=["teachers", "pupils", "cs_mn_lrn_ol"])
    q = d.teachers.quantile([.1, .25, .5, .75, .9])
    say("   teachers per district, unweighted: "
        + ", ".join(f"p{int(k * 100)}={v:,.0f}" for k, v in q.items()))
    say(f"   enrolment-weighted median = {de.wquantile(d.teachers, d.pupils, .5):,.0f}"
        f", mean = {wmean(d.teachers, d.pupils):,.0f}")
    tot, err, _, _ = reliability(panel, "cs_mn_lrn_ol", "cs_mn_lrn_ol_se",
                                 "pupils", "czone")
    sd_obs = np.sqrt(max(tot - err, 0))
    say(f"\n   observed within-CZ sd of the district learning rate = {sd_obs:.4f}"
        " SD/grade (noise removed)")
    say(f"   {VA_WITHIN_SCHOOL_SHARE:.0%} of teacher-VA variance is within school, "
        "not between (CFR I),")
    say("   and a district mean averages over hundreds of teachers:\n")
    say(f"   {'sd of teacher VA':>17s} {'teachers':>9s} {'sd of district mean':>20s}"
        f" {'observed / that':>16s}")
    for sdva in (VA_SD_ELA_CFR, VA_SD_MATH_CFR, VA_SD_MATH_TRUE):
        for N in (50, 200, 1000):
            samp = sdva / np.sqrt(N)
            say(f"   {sdva:17.3f} {N:9d} {samp:20.4f} {sd_obs / samp:15.1f}x")
    say("\n   Read: district growth varies several times more than random")
    say("   allocation of teachers would produce, so there IS systematic sorting")
    say("   in it -- but what survives aggregation is teacher sorting bundled")
    say("   with curriculum, leadership, peers and everything else the district")
    say("   does.  No district-level series separates the bundle.")

    say("\n   Two regressions the main report runs on the locale partition only.")
    say("   The sample, not the data, is what makes the kappa gradient look dead:\n")
    loc = de.partition(panel, "locale")
    loc = loc[loc.year == YEAR]
    say(f"   {'':5s} {'sample':22s} {'controls':16s} {'b':>9s} {'se':>8s} {'t':>6s}"
        f" {'n':>7s}")
    for nm, d_ in (("full within-CZ panel", panel), ("locale partition", loc)):
        for ctl, lab in (([], "none"),
                         (["log_str_ratio", "pov_rate_5_17"], "str + poverty")):
            dd = d_.dropna(subset=["log_salary_real", "learn_gain", "pupils",
                                   "czone"] + ctl)
            r = de.ols_fe(dd, "log_salary_real", ["learn_gain"] + ctl, fe="czone",
                          w="pupils", cluster="czone")
            if r:
                b, s = r["coef"]["learn_gain"], r["se"]["learn_gain"]
                say(f"   kappa {nm:22s} {lab:16s} {b:+9.4f} {s:8.4f}"
                    f" {b / s:6.2f} {r['n']:7d}")
    for nm, d_ in (("full within-CZ panel", panel), ("locale partition", loc)):
        for ctl, lab in (([], "none"),
                         (["seda_ses", "pov_rate_5_17"], "SES + poverty")):
            dd = d_.dropna(subset=["log_Q", "log_teachers_pp", "pupils",
                                   "czone"] + ctl)
            r = de.ols_fe(dd, "log_Q", ["log_teachers_pp"] + ctl, fe="czone",
                          w="pupils", cluster="czone")
            if r:
                b, s = r["coef"]["log_teachers_pp"], r["se"]["log_teachers_pp"]
                say(f"   sigma {nm:22s} {lab:16s} {b:+9.4f} {s:8.4f}"
                    f" {b / s:6.2f} {r['n']:7d}")
    say(f"\n   Errors-in-variables: the learning rate's within-CZ reliability is")
    say(f"   {rel_within:.3f}, so putting it on the right-hand side attenuates a")
    say(f"   coefficient by only {(1 / rel_within - 1) * 100:.1f}%.  Measurement error is not the")
    say("   story; sigma is genuinely not identified in this cross-section, and")
    say("   the kappa gradient is genuinely positive once all districts are used.")


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)-7s %(message)s")
    panel = load_panel()
    say(f"Teacher quality in SEDA -- {YEAR} CCD x SEDA 6.0 pooled 2009-19")
    say(f"panel: {len(panel):,} districts, "
        f"{int(panel.cs_mn_lrn_ol.notna().sum()):,} with a SEDA learning rate, "
        f"{panel.czone.nunique()} commuting zones")
    rel = block_a(panel)
    block_b(panel)
    block_c(panel)
    res = block_d(panel)
    block_e(res, panel)
    block_f(panel, rel)
    outdir = SPATIAL / "estimates"
    outdir.mkdir(exist_ok=True)
    (outdir / "teacher_quality.txt").write_text("\n".join(_out) + "\n")
    LOG.info("wrote %s", outdir / "teacher_quality.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
