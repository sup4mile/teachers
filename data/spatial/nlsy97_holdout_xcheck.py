#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nlsy97_holdout_xcheck.py -- independent re-computation of a subset of the SEALED NLSY97 hold-out statistics,
compared with the sealed file.  Prints only PASS/FAIL and maximum absolute differences, never a value.

`nlsy97_holdout.py` writes the hold-out statistics to estimates/settings/sealed/nlsy97_holdout.json without
printing them, so the only way to know they were computed as intended is to compute some of them a second
time, in code written separately (pandas on the raw columns, no `nlsy_composition` / `nlsy_teachers` helper
except the age-normed ability z-scores of `nlsy_ability.build_n97`, which reproduce the published
R_comp_97).  This script does that for 66 statistics -- the transition statistics by own-ability tercile
(origin city share, city to suburb, suburb to city, move rate, with top minus bottom), the teacher gaps A, B
and C of the main specification (employed, all K-12 codes, majority rule; pooled, men, women; BA+ by degree
and 16+ years), the BA+ / 16+ share by ability tercile and the BA+ vs non-BA gap in the whole-sample and in
the NLSY79-matching universe (highest grade >= 9, terciles cut over that sample; pooled, men, women), and
R_comp_97 -- and, since the first version, derives independently from the raw columns the things the first
version took from `Frame97`: the household clusters (breadth-first search on the reciprocated round-1 roster
links, households numbered by their lowest respondent index), the PSU index and the Rao-Wu draws
(`VSTRAT`, `VPSU`; `rng.integers(0, 2, strata)` per draw, one PSU per stratum with weight 2), and the
armed-forces rounds (a military roster employer held at the interview and no civilian one).  It compares:

  (0) the independent household partition, PSU multipliers of the first draw and armed-forces rounds with
      those of `Frame97` (equal);
  (a) the statistics at the sample, against the sealed point estimates;
  (b) at the first draw of each bootstrap (household: seed 20260929; PSU: seed 20260929 + 97), against
      `nlsy97_holdout.evaluate97` at that draw's multipliers;
  (c) with `--boot`, both whole bootstraps (B = 500) against the sealed household and PSU standard errors
      and 95% percentile limits.

What it does NOT validate: how the raw variables are read (the codebook meanings of `CV_MSA`, `CV_MAINJOB_FLG`,
`YEMP_CURFLAG`, `YEMP_MILFLAG`, `CV_HGC_EVER_EDT`, `CV_HIGHEST_DEGREE_EVER_EDT`: the independent code uses the
same `nlsy_ability.Cohort` accessor and the same readings, documented in the methods file) or the choice of
definitions; it validates the arithmetic, the bootstrap protocols and the derived structures.

    cd data/spatial && python3 nlsy97_holdout_xcheck.py [path/to/nlsy97_holdout.json] [--boot [B]]

The default path is estimates/settings/sealed/nlsy97_holdout.json.  The file is read by this program only;
nothing from it is displayed.  (A failure lists the keys with the largest differences, not their values.)

Written with help from Claude Code.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
logging.basicConfig(level=logging.ERROR)
import nlsy_ability as na  # noqa: E402
import nlsy97_holdout as h  # noqa: E402

TOL = 1e-9
SEED = 20260929
PSU_SEED = SEED + 97


def wmean(x: np.ndarray, w: np.ndarray) -> float:
    return float((x * w).sum() / w.sum())


def wsd(x: np.ndarray, w: np.ndarray) -> float:
    m = wmean(x, w)
    return float(np.sqrt((w * (x - m) ** 2).sum() / w.sum()))


def wq(x: np.ndarray, w: np.ndarray, q: float) -> float:
    """Weighted quantile with the project's interpolation convention (`nlsy_ability.wquantile`)."""
    o = np.argsort(x, kind="stable")
    xs, ws = x[o], w[o]
    cum = np.cumsum(ws) - 0.5 * ws
    return float(np.interp(q, cum / ws.sum(), xs))


def terc(x: np.ndarray, wcut: np.ndarray, mask: np.ndarray) -> np.ndarray:
    q1, q2 = wq(x[mask], wcut[mask], 1 / 3), wq(x[mask], wcut[mask], 2 / 3)
    return np.where(~np.isfinite(x), -1, np.where(x <= q1, 0, np.where(x <= q2, 1, 2)))


def households_bfs(pid: np.ndarray, sib: np.ndarray) -> np.ndarray:
    """Household label of every respondent: connected components of the RECIPROCATED round-1 roster links by
    breadth-first search, numbered by the lowest respondent index in the component."""
    n = len(pid)
    row_of = {int(p): i for i, p in enumerate(pid)}
    directed = set()
    for i in range(n):
        for s in sib[i]:
            if np.isfinite(s):
                j = row_of.get(int(s))
                if j is not None and j != i:
                    directed.add((i, j))
    nbr: list[list[int]] = [[] for _ in range(n)]
    for (i, j) in directed:
        if (j, i) in directed:
            nbr[i].append(j)
    lab = np.full(n, -1)
    k = 0
    for s0 in range(n):
        if lab[s0] >= 0:
            continue
        lab[s0] = k
        stack = [s0]
        while stack:
            u = stack.pop()
            for v in nbr[u]:
                if lab[v] < 0:
                    lab[v] = k
                    stack.append(v)
        k += 1
    return lab


class Independent:
    """The fixed inputs of the independent computation, `psu_mult(rng)` and `compute(mult)`."""

    def __init__(self, c: na.Cohort) -> None:
        self.c = c
        n = self.n = c.df.shape[0]
        R = na.build_n97(c)                    # age-normed CAT-ASVAB z-scores (validated by the published R_comp_97)
        self.w0 = c.get("SAMPLING_WEIGHT_CC", 1997) / 100.0
        self.zn = R.zn
        self.female = R.female.astype(bool)
        age97 = c.get("CV_AGE_INT_DATE", 1997)
        self.origin_age = (age97 >= 12) & (age97 <= 16) & (self.w0 > 0)
        # households (reciprocated links) and the design (variance stratum x PSU), from the raw columns
        pid = c.get("PUBID", 1997)
        sib = np.column_stack([c.get(f"YOUTH_SIBID0{k}.01", 1997) for k in range(1, 5)])
        self.hh = households_bfs(pid, sib)
        self.ncl = int(self.hh.max()) + 1
        vs, vp = c.get("VSTRAT"), c.get("VPSU")
        strata = np.unique(vs)
        self.H = len(strata)
        self.srank = np.searchsorted(strata, vs)
        self.vp = vp.astype(int)
        # transition rows: origin = 1997 CV_MSA (3 city, 2 suburb, 1 not in MSA), adult = rounds at 25-34
        msa97 = c.get("CV_MSA", 1997)
        orig = np.select([msa97 == 3, msa97 == 2, msa97 == 1], [0, 1, 2], default=-1)
        orig = np.where(self.origin_age, orig, -1)
        rows = []
        for Y in sorted(int(y) for y in c.by_q["CV_MSA"] if str(y).isdigit()):
            age = c.get("CV_AGE_INT_DATE", Y)
            m = c.get("CV_MSA", Y)
            w = c.get("SAMPLING_WEIGHT_CC", Y) / 100.0
            a = np.select([m == 3, m == 2, m == 1], [0, 1, 2], default=-1)
            sel = np.isfinite(age) & (age >= 25) & (age <= 34) & np.isfinite(m) & (w > 0) & (a >= 0) & (orig >= 0)
            rows += [(i, a[i], w[i]) for i in np.flatnonzero(sel)]
        tr = pd.DataFrame(rows, columns=["p", "a", "w"])
        tr["rw"] = tr.w / tr.groupby("p").p.transform("size")      # a person carries the mean of his or her round weights
        tr["o"] = orig[tr.p]
        self.tr = tr
        # person-rounds at 25-34: main-job occupation, current-job flag, schooling, armed-forces status
        pr_rows = []
        for Y in sorted(int(y) for y in c.by_q["CV_AGE_INT_DATE"] if str(y).isdigit()):
            age = c.get("CV_AGE_INT_DATE", Y)
            main = c.try_get("CV_MAINJOB_FLG", Y)
            deg = c.try_get("CV_HIGHEST_DEGREE_EVER_EDT", Y)
            hgc = c.try_get("CV_HGC_EVER_EDT", Y)
            mfl = {nn: c.try_get(f"YEMP_MILFLAG.{nn:02d}", Y) for nn in range(1, 16)}
            cfl = {nn: c.try_get(f"YEMP_CURFLAG.{nn:02d}", Y) for nn in range(1, 16)}
            occ_ = {nn: c.try_get(f"YEMP_OCCODE-2002.{nn:02d}", Y) for nn in range(1, 16)}
            for i in np.flatnonzero(np.isfinite(age) & (age >= 25) & (age <= 34)):
                occ = cur = np.nan
                if main is not None and np.isfinite(main[i]):
                    nn = int(main[i])
                    occ = occ_[nn][i] if occ_.get(nn) is not None else np.nan
                    cur = cfl[nn][i] if cfl.get(nn) is not None else np.nan
                mil_cur = civ_cur = False
                for nn in range(1, 16):
                    if cfl[nn] is None or not cfl[nn][i] == 1:
                        continue
                    if mfl[nn] is not None and mfl[nn][i] == 1:
                        mil_cur = True
                    else:
                        civ_cur = True
                g = np.nan if hgc is None or hgc[i] >= 95 else hgc[i]      # 95 = ungraded
                pr_rows.append((i, Y, occ, cur, np.nan if deg is None else deg[i], g, mil_cur and not civ_cur))
        pr = pd.DataFrame(pr_rows, columns=["p", "Y", "occ", "cur", "deg", "hgc", "armed"])
        self.pr = pr
        pr["valid"] = pr.occ.between(10, 9750)
        pr["k12"] = pr.valid & pr.occ.isin([2300, 2310, 2320, 2330, 2340])
        emp = pr[(~pr.armed) & pr.valid & (pr.cur == 1)]
        g = emp.groupby("p").agg(nr=("k12", "size"), nk=("k12", "sum"))
        self.teacher = pd.Series(2 * g.nk > g.nr, index=g.index)       # majority of the employed rounds
        self.observed = np.zeros(n, bool)
        self.observed[pr[~pr.armed].p.unique()] = True                 # at least one round that is not an armed-forces round
        self.deg_max = pr.groupby("p").deg.max().reindex(np.arange(n)).to_numpy()
        self.hgc_max = pr.groupby("p").hgc.max().reindex(np.arange(n)).to_numpy()

    def psu_mult(self, rng: np.random.Generator) -> np.ndarray:
        """Person multipliers of one Rao-Wu draw: one of the two PSUs of each stratum (equal probability) gets weight 2."""
        pick = rng.integers(0, 2, self.H)
        return np.where(self.vp - 1 == pick[self.srank], 2.0, 0.0)

    def hh_mult(self, rng: np.random.Generator) -> np.ndarray:
        """Person multipliers of one household-bootstrap draw."""
        cm = np.bincount(rng.integers(0, self.ncl, self.ncl), minlength=self.ncl).astype(float)
        return cm[self.hh]

    def compute(self, mult: np.ndarray) -> dict[str, float]:
        """The independent versions of the statistics for person multiplicities `mult`."""
        n, w0, zn = self.n, self.w0, self.zn
        w = w0 * mult
        ok = np.isfinite(zn).all(axis=1) & (w > 0)
        zs = np.column_stack([(zn[ok, k] - wmean(zn[ok, k], w[ok])) / wsd(zn[ok, k], w[ok]) for k in range(4)])
        comp = zs.mean(axis=1)
        comp = (comp - wmean(comp, w[ok])) / wsd(comp, w[ok])
        af = np.full(n, np.nan)
        af[ok] = comp
        Rm = na.wcorr_cross(zn[ok], zn[ok], w[ok])
        fit = na.fit_factor_uls(Rm, [0] * 4)
        rc = (fit.lam.sum() ** 2) / Rm.sum()
        sR = np.sqrt(rc)
        t_origin = terc(af, w, self.origin_age & ok)
        t_full = terc(af, w, ok)
        out = {"R_comp_97": rc}
        # items 1-3
        d = self.tr.copy()
        d["rwm"] = d.rw * mult[d.p]
        d["t"] = t_origin[d.p]
        d2 = d[(d.o <= 1) & (d.a <= 1) & (d.t >= 0)]
        sts = {k: {} for k in range(3)}
        for k in range(3):
            gg = d2[d2.t == k]
            co, so = gg[gg.o == 0], gg[gg.o == 1]
            sts[k]["origin_share_city"] = co.rwm.sum() / gg.rwm.sum()
            sts[k]["city_to_suburb"] = co[co.a == 1].rwm.sum() / co.rwm.sum()
            sts[k]["suburb_to_city"] = so[so.a == 0].rwm.sum() / so.rwm.sum()
            sts[k]["move_rate"] = (co[co.a == 1].rwm.sum() + so[so.a == 0].rwm.sum()) / gg.rwm.sum()
            for key, v in sts[k].items():
                out[f"trans|afqt|{key}|{k + 1}"] = v
        for key in sts[0]:
            out[f"trans|afqt|{key}|diff"] = sts[2][key] - sts[0][key]
        # item 4: teachers (employed, all K-12, majority); college = BA+ by degree or 16+ years
        pers = pd.DataFrame({"p": np.arange(n), "af": af, "w": w, "female": self.female, "deg_max": self.deg_max,
                             "hgc_max": self.hgc_max})
        pers["teacher"] = self.teacher.reindex(pers.p).to_numpy()
        members = pers[pers.teacher.notna() & pers.af.notna() & (pers.w > 0)].copy()
        members["teacher"] = members.teacher.astype(bool)
        thr = {"ba": ("deg_max", 4), "c16": ("hgc_max", 16)}
        for grp, mask in (("pooled", np.ones(len(members), bool)), ("men", ~members.female.to_numpy()),
                          ("women", members.female.to_numpy())):
            m = members[mask]
            t, nt_ = m[m.teacher], m[~m.teacher]
            mt = wmean(t.af.to_numpy(), t.w.to_numpy())
            out[f"c02-emp-all-maj|{grp}|nt|gap"] = (mt - wmean(nt_.af.to_numpy(), nt_.w.to_numpy())) / sR
            for d_, (col, lim) in thr.items():
                nb, tb = nt_[nt_[col] >= lim], t[t[col] >= lim]
                out[f"c02-emp-all-maj|{grp}|col_{d_}|gap"] = (mt - wmean(nb.af.to_numpy(), nb.w.to_numpy())) / sR
                out[f"c02-emp-all-maj|{grp}|cc_{d_}|gap"] = (wmean(tb.af.to_numpy(), tb.w.to_numpy())
                                                              - wmean(nb.af.to_numpy(), nb.w.to_numpy())) / sR
        # item 5: share with BA+ (16+ years) by ability tercile and the gap in latent ability, whole-sample universe (`edu`)
        # and the NLSY79-matching universe (`edu9`: highest grade >= 9; `edu9_cutm`: terciles cut over that sample)
        pers["t_full"] = t_full
        pers["obs"] = self.observed
        masks = {"pooled": np.ones(n, bool), "men": ~self.female, "women": self.female}
        for d_, (col, lim) in thr.items():
            known = pers[col].notna().to_numpy() & pers.obs.to_numpy() & pers.af.notna().to_numpy() & (pers.w > 0).to_numpy()
            mem = pers[known]
            flag = (mem[col] >= lim)
            for k in range(3):
                g_ = mem[mem.t_full == k]
                out[f"edu|{d_}|share|{k + 1}"] = g_.w[g_[col] >= lim].sum() / g_.w.sum()
            b1, b0 = mem[flag], mem[~flag]
            out[f"edu|{d_}|pooled|gap"] = (wmean(b1.af.to_numpy(), b1.w.to_numpy()) - wmean(b0.af.to_numpy(), b0.w.to_numpy())) / sR
            # grade >= 9 universe
            sel9 = (pers.hgc_max >= 9).to_numpy()
            mem9 = pers[known & sel9]
            t9 = terc(af, w, known & sel9)
            mem9 = mem9.assign(t9=t9[mem9.p.to_numpy()], men=~mem9.female.to_numpy())
            for grp, gm in (("pooled", np.ones(len(mem9), bool)), ("men", mem9.men.to_numpy()), ("women", ~mem9.men.to_numpy())):
                g9 = mem9[gm]
                sh = []
                for k in range(3):
                    gq = g9[g9.t9 == k]
                    sh.append(gq.w[gq[col] >= lim].sum() / gq.w.sum())
                    out[f"edu9_cutm|{d_}|share|{k + 1}" if grp == "pooled" else f"edu9_cutm|{d_}|{grp}|share|{k + 1}"] = sh[-1]
                if grp == "pooled":
                    out[f"edu9_cutm|{d_}|share|diff"] = sh[2] - sh[0]
                f1, f0 = g9[g9[col] >= lim], g9[g9[col] < lim]
                out[f"edu9|{d_}|{grp}|gap"] = (wmean(f1.af.to_numpy(), f1.w.to_numpy()) - wmean(f0.af.to_numpy(), f0.w.to_numpy())) / sR
        return out


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("--") and not a.isdigit()]
    path = Path(args[0]) if args else HERE / "estimates" / "settings" / "sealed" / "nlsy97_holdout.json"
    S = json.loads(path.read_text())["statistics"]
    c = na.Cohort("nlsy97", HERE)
    ind = Independent(c)
    fr = h.Frame97(c)                          # for the first-draw comparisons and the structures checked in (0)
    good = True

    def verdict(name: str, worst: float, d: dict[str, float] | None = None) -> None:
        nonlocal good
        flag = worst < TOL
        good &= flag
        print(f"{name}: max abs difference {worst:.1e}  [{'PASS' if flag else 'FAIL'}]")
        if not flag and d:
            for k in sorted(d, key=d.get, reverse=True)[:6]:
                print(f"  {k}: {d[k]:.1e}")

    # (0) independently derived structures equal those of Frame97
    same_hh = bool(np.array_equal(ind.hh, fr.hh_idx))
    same_armed = bool(np.array_equal(ind.pr.armed.to_numpy(), fr.pr.armed))
    same_p = bool(np.array_equal(ind.pr.p.to_numpy(), fr.pr.p))
    m_main = h.rao_wu_multipliers(fr.psu_pairs, fr.n_psu, np.random.default_rng(PSU_SEED))[fr.psu_idx]
    same_psu = bool(np.array_equal(ind.psu_mult(np.random.default_rng(PSU_SEED)), m_main))
    ok0 = same_hh and same_armed and same_p and same_psu
    good &= ok0
    print(f"independent household partition ({ind.ncl} households), person-round rows, armed-forces rounds and first Rao-Wu multipliers "
          f"equal those of Frame97: households {same_hh}, rows {same_p}, armed {same_armed}, PSU multipliers {same_psu}  [{'PASS' if ok0 else 'FAIL'}]")

    mine = ind.compute(np.ones(ind.n))
    d_s = {k: abs(v - S[k]["value"]) for k, v in mine.items()}
    verdict(f"sample: {len(mine)} statistics (transitions, teacher gaps A/B/C, BA+/16+ shares and gaps in both universes, R_comp_97) "
            "vs the sealed point estimates", max(d_s.values()), d_s)

    for lab, mult, key_se, key_lo, key_hi in (
            ("household", ind.hh_mult(np.random.default_rng(SEED)), "se", "ci95_lo", "ci95_hi"),
            ("PSU", ind.psu_mult(np.random.default_rng(PSU_SEED)), "se_psu", "ci95_lo_psu", "ci95_hi_psu")):
        mine = ind.compute(mult)
        ref = h.evaluate97(fr, mult)
        d_b = {k: abs(v - ref[k]) for k, v in mine.items()}
        verdict(f"first {lab} bootstrap draw: the same statistics vs nlsy97_holdout.evaluate97 at that draw's multipliers", max(d_b.values()), d_b)

    if "--boot" in argv:
        nums = [a for a in argv if a.isdigit()]
        B = int(nums[0]) if nums else 500
        keys = list(mine)
        rng_h, rng_p = np.random.default_rng(SEED), np.random.default_rng(PSU_SEED)
        arr_h, arr_p = np.full((B, len(keys)), np.nan), np.full((B, len(keys)), np.nan)
        for b in range(B):
            d = ind.compute(ind.hh_mult(rng_h))
            arr_h[b] = [d[k] for k in keys]
            d = ind.compute(ind.psu_mult(rng_p))
            arr_p[b] = [d[k] for k in keys]
        for lab, arr, ks, kl, kh in (("household", arr_h, "se", "ci95_lo", "ci95_hi"), ("PSU", arr_p, "se_psu", "ci95_lo_psu", "ci95_hi_psu")):
            se = np.nanstd(arr, axis=0, ddof=1)
            lo, hi = np.nanpercentile(arr, [2.5, 97.5], axis=0)
            dd = {k: max(abs(se[i] - S[k][ks]), abs(lo[i] - S[k][kl]), abs(hi[i] - S[k][kh])) for i, k in enumerate(keys)}
            verdict(f"independent {lab} bootstrap (B = {B}) of the same statistics vs the sealed standard errors and 95% percentile limits",
                    max(dd.values()), dd)
    print(f"independent re-computation: {'PASS' if good else 'FAIL'}")
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
