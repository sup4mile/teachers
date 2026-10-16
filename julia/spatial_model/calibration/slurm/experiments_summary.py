#!/usr/bin/env python3
"""
experiments_summary.py -- the fixed validation table for a set of refits (2026-10 experiments).

Reads <dir>/{theta.toml, moments.toml, report.txt} as written by `spatial_diagnose.jl --out`
for each case directory and writes a markdown report: the criterion and ϑ of each case, the
standardized residuals of its targets, and the pre-registered validation table (district
gaps, NLSY composition and flows by parental income, teacher pay and schooling, the gender
gap, the first-stage consistency moments) with the composition diagnostics and a block of model
levels (population shares, tax rates, quality, home production, WTP elasticity, the floor-mass
audit), each with its data comparator where there is one. A moment in a
case's criterion is marked "†"; residuals are (model − data)/scale. ϑ is labelled by the
`theta_names` of each case's theta.toml (the NAMES list for files without them), so cases with
different numbers of parameters share one table: an extension parameter a case does not list is
off (0), as `pad_theta` reads it.

Case directories are given as --case NAME=DIR (repeatable) or found as the subdirectories of
--runs ROOT (repeatable). The district data values default to the locale targets; --scheme
reads another partition's values from data/spatial/estimates/settings/district_settings_targets.toml.
Run from the repository root:

    python3 julia/spatial_model/calibration/slurm/experiments_summary.py \\
        --case baseline=julia/spatial_model/calibration/runs/exp-final/p0_baseline \\
        --runs julia/spatial_model/calibration/runs/panel-exact --out <file.md>

Written with help from Claude Code.
"""
import argparse
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SETTINGS = ROOT / "data/spatial/estimates/settings/district_settings_targets.toml"
# names of ϑ in order (spatial_calibrate.jl THETA_NAMES: the base eight, r_p, ζ, c_T, c_Tz, then ρw
# (the split), the M3 education block, ψ_a (peer effects, 2026-10-07), the switches of 2026-10-15,
# κ_c, θ_ν, ζ_BA and b_T, and ζ_y of round 2026-10-15c); a theta.toml without `theta_names` is read against this list
NAMES = ["logκ̃", "δκ", "ΔB", "β", "λ", "r_m", "ωf", "γ", "r_p", "ζ", "c_T", "c_Tz",
         "ρw", "ρ_c", "c_C", "c_Cz", "ρ_T0", "ϑ_s", "s_z", "σ_ϵ", "dlogφ", "c_T0", "ψ_a",
         "κ_c", "θ_ν", "ζ_BA", "b_T", "ζ_y"]

# key: (label, data value, scale); scales are the calibration's (SEs or tolerances)
DISTRICT = {"gap_income": ("earnings gap", 0.1122, 0.0168),
            "gap_teachers_pp": ("FTE/pupil gap", -0.0177, 0.0131),
            "payroll_pp_gap": ("payroll/pupil gap", -0.0091, 0.0188),
            "poverty_gap": ("poverty gap", -0.0888, 0.0090),
            "score_gap": ("SEDA score gap", 0.2887, 0.0326),
            "seda_growth_gap": ("SEDA growth gap", 0.0031, 0.0030)}
# residence/workplace split (`--external split=1`; the commuting shares, the net commuting gap and the
# supply elasticity are absent for the bundled model: NaN there, dropped from moments.toml). K-12 teacher
# values from data/spatial/estimates/settings/teacher_residence_targets.toml [k12]; the national PUMS net
# commuting gap +0.108 (SE 0.008) is only a cross-check
SPLIT = {"commute_share_all": ("crossing share, all workers", 0.107, 0.03),
         "commute_share_teach": ("crossing share, teachers", 0.084, 0.03),
         "net_commute_gap": ("net commuting gap", 0.075, 0.038),
         "gap_teach_res_emp": ("teachers per employed resident gap", 0.298, 0.041),
         "teacher_supply_elasticity": ("teacher supply elasticity (Ransom–Sims; Falch 1.4)", 3.7, 1.0)}
OTHER = {"city_share_gradient": ("city share top−bottom", -0.295, 0.037),
         "teacher_premium_m": ("teacher premium m", 0.0557, 0.05),
         "teacher_premium_f": ("teacher premium f", 0.0860, 0.05),
         "s_ratio": ("s_T/s_O", 1.1761, 0.0588),
         "gender_gap_market": ("gender gap", 0.1046, 0.05),
         "teacher_z_gap_market_sd": ("teachers' z vs market (SD)", 0.614, 0.062),
         "city_share_gradient_z": ("city share by own z, top−bottom", -0.222, 0.036),
         "wage_slope": ("wage slope", 0.2127, 0.0077),
         "rho_pc": ("mother–child ρ", 0.5766, 0.0140)}
TARGETS = ["male_teach_share", "female_teach_share", "teacher_p9010", "gap_salary", "gap_pupils",
           "cfr_effect", "wtp_share", "move_rate"]
# the joint education and occupation choice (M3, `--external college=1`): label, data value, scale
M3 = {"ba_share_nonteach": ("BA share, non-teachers", 0.3231, 0.0162),
      "ba_share_teach": ("BA share, teachers", 0.8744, 0.01),
      "ba_z_gap_sd": ("BA − non-BA mean log z (SD)", 1.2189, 0.019),
      "ba_premium": ("BA premium, within occupation", 0.3610, 0.02),
      "teacher_premium_ba_m": ("teachers − BA non-teachers, men", -0.2610, 0.05),
      "teacher_premium_ba_f": ("teachers − BA non-teachers, women", -0.2081, 0.05),
      "teacher_premium_within_ba_m": ("BA teachers − BA non-teachers, men", -0.2374, 0.05),
      "teacher_premium_within_ba_f": ("BA teachers − BA non-teachers, women", -0.1368, 0.05),
      "teacher_z_gap_ba_sd": ("teachers' z − BA non-teachers (SD)", -0.3335, 0.0625),
      "teacher_z_gap_within_ba_sd": ("BA teachers' z − BA non-teachers (SD)", -0.2384, 0.0617),
      "s_T_mean_ba": ("mean s, BA teachers", 0.6573, 0.01),
      "s_O_mean": ("mean s, non-teachers", 0.5457, 0.005),
      "wage_slope": ("wage slope", 0.2127, 0.0077),
      "wage_slope_ba": ("wage slope holding BA fixed (convention)", 0.1608, 0.0083),
      "pooled_p9010": ("pooled 90/10", 3.737, 0.187),
      "ma_premium_implied": ("master's premium implied (17/16)", 0.1564, 0.0068),
      "teacher_ba_wage_gap": ("BA − non-BA teachers' log wage", 0.4819, 0.0094),
      "s_T_mean_nonba": ("mean s, non-BA teachers", 0.5195, 0.01),
      "teacher_p9010_ba": ("teacher 90/10, BA teachers only", 2.4166, 0.1208),
      "ba_share_z_gap": ("P(BA | top z tercile) − P(bottom)", 0.653, 0.05),
      "ba_share_score_gap": ("P(BA | top score tercile) − P(bottom), observed-score terciles", 0.6243, 0.05),
      "teacher_slope_ratio": ("teachers' / BA non-teachers' wage slope", 0.5565, 0.3188)}
VECTORS = {"city_share_by_parent_income": (0.458, 0.360, 0.163),
           "c2s_by_parent_income": (0.284, 0.429, 0.407),
           "s2c_by_parent_income": (0.168, 0.112, 0.168),
           "move_rate_by_parent_income": (0.221, 0.226, 0.207)}


def payroll_gap(tab):
    """Payroll per pupil is pay per FTE times FTE per pupil: its gap is the sum of the salary and FTE
    gaps of the scheme's table, with the SE of the sum from their bootstrap correlation."""
    ck = tab["corr_keys"]
    rho = tab["corr"][ck.index("gap_salary")][ck.index("gap_teachers_pp")]
    s, f = tab["gap_salary"], tab["gap_teachers_pp"]
    return s["value"] + f["value"], (s["se"] ** 2 + f["se"] ** 2 + 2 * rho * s["se"] * f["se"]) ** 0.5


def district_values(scheme):
    if scheme in (None, "locale"):
        return DISTRICT
    tab = tomllib.loads(SETTINGS.read_text())[scheme]
    out = {k: (lab, tab[k]["value"], tab[k]["se"]) if k in tab else (lab, v, s)
           for k, (lab, v, s) in DISTRICT.items()}
    out["payroll_pp_gap"] = (DISTRICT["payroll_pp_gap"][0], *payroll_gap(tab))
    return out


def read_case(d: Path):
    th, mo, rep = d / "theta.toml", d / "moments.toml", d / "report.txt"
    if not (th.exists() and mo.exists()):
        return None
    t = tomllib.loads(th.read_text())
    m = tomllib.loads(mo.read_text())
    if "payroll_pp_gap" not in m and "gap_salary" in m and "gap_teachers_pp" in m:
        m["payroll_pp_gap"] = m["gap_salary"] + m["gap_teachers_pp"]     # files from before the key
    if "option_value_gap" not in m and "option_value_teach" in m and "option_value_nonteach" in m:
        m["option_value_gap"] = m["option_value_teach"] - m["option_value_nonteach"]   # files from before the key
    text = rep.read_text() if rep.exists() else ""
    values = list(t["theta"])
    names = list(t.get("theta_names") or NAMES[:len(values)])
    theta = dict(zip(names, values))
    g = re.search(r"pooled non-teacher 90/10: GE ([\d.]+)", text)
    it = re.search(r"polish: (\d+) iterations, ([\d.]+) minutes", text)
    opts = re.search(r"options: (.*)", text)
    return {"theta": theta, "crit": t["criterion"], "critw": t.get("criterion_weighted"),
            "m": m, "active": set(m.get("active", [])), "p9010": float(g.group(1)) if g else None,
            "iters": f"{it.group(1)} ({float(it.group(2)):.0f} min)" if it else "—",
            "options": opts.group(1) if opts else ""}


def cell(r, key, val, se, digits=3, own=True):
    """Model value and residual against (val, se); with `own`, a moment that is a target of the case
    uses the case's own value and scale instead (district columns pass own=False, so that they are
    scored against the header's data, the --scheme values)."""
    x = r["m"].get(key)
    if x is None:
        return "—"
    t = r["m"].get("targets", {}).get(key) if own else None
    if t:
        val, se = t["value"], t["se"]
    mark = "†" if key in r["active"] else ""
    return f"{x:.{digits}f}{mark} ({(x - val) / se:+.1f})"


def vec(r, key, digits=3):
    x = r["m"].get(key)
    return "—" if x is None else ", ".join(f"{v:.{digits}f}" for v in x)


def theta_columns(rows):
    """Names of ϑ for the parameter table: the names of the cases in order of first appearance,
    starting from the eight base parameters."""
    cols = NAMES[:8]
    for _, r in rows:
        cols += [n for n in r["theta"] if n not in cols]
    return cols


def is_split(r):
    """The case has a workplace: its crossing share is present and not the bundled model's zero."""
    c = r["m"].get("commute_share_all")
    return c is not None and c != 0.0


def level(r, keys, fmt):
    """`fmt` applied to the values of `keys` in the case's moments (scalars, or vectors joined by
    commas), or "—" if the case lacks one of them (moments.toml files from before the keys existed)."""
    if any(k not in r["m"] for k in keys):
        return "—"
    flat = [v for k in keys for v in (r["m"][k] if isinstance(r["m"][k], list) else [r["m"][k]])]
    return ", ".join(format(v, fmt) for v in flat)


# model levels: (header, moment keys, format); a column no case has the keys of is left out. A data
# comparator is in brackets; "model only" levels have none
LEVELS = [("population share (city, suburb) [suburb 0.579]", ["pop_share"], ".3f"),
          ("t on all income (city, suburb; model only)", ["t"], ".4f"),
          ("t on market income, salaries only (all; city, suburb) [≈ 0.016]", ["t_market", "t_market_l"], ".4f"),
          ("Q (city, suburb; model only)", ["Q"], ".4f"),
          ("ΔlogQ", ["gap_logQ"], "+.4f"),
          ("students/teacher (model only)", ["students_per_teacher"], ".1f"),
          ("teaching share (men, women)", ["male_teach_share", "female_teach_share"], ".4f"),
          ("home production (men, women)", ["hp_share_m", "hp_share_f"], ".3f"),
          ("mean s (teachers, non-teachers)", ["s_T_mean", "s_O_mean"], ".4f"),
          ("P(1→2), P(2→1)", ["p12", "p21"], ".3f"),
          ("mean C (model only)", ["mean_C"], ".5f"),
          ("WTP income elasticity, all [BFM 0.379]", ["wtp_income_elasticity"], ".3f"),
          ("WTP income elasticity within z [≈ 1 under log utility]", ["wtp_income_elasticity_z"], ".3f"),
          ("WTP share by parental-income tercile", ["wtp_share_by_parent_income"], ".4f"),
          ("mean ΔΛ (all; due to Q alone)", ["dLambda_mean", "dLambda_quality_mean"], ".2e"),
          ("teachers by residence; by workplace (model only)", ["n_teach_res", "n_teach_work"], ".4f"),
          ("teacher supply elasticity, split only [Ransom–Sims 3.7; Falch ≈ 1.4: upper bounds]", ["teacher_supply_elasticity"], ".2f"),
          ("population share at the floor (feasible location exists; none)", ["floor_mass", "infeasible_mass"], ".1e")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append", default=[], help="NAME=DIR")
    ap.add_argument("--runs", action="append", default=[], help="directory of case directories")
    ap.add_argument("--scheme", default=None, help="district scheme of the data values")
    ap.add_argument("--title", default="Calibration experiments: fixed validation table")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cases = [(n, ROOT / d) for n, d in (c.split("=", 1) for c in a.case)]
    for root in a.runs:
        cases += sorted((p.name, p) for p in (ROOT / root).iterdir() if p.is_dir())
    rows = [(n, r) for n, r in ((n, read_case(d)) for n, d in cases) if r is not None]
    missing = [n for n, d in cases if read_case(d) is None]
    dist = district_values(a.scheme)

    L = [f"# {a.title}", "",
         "Generated by `julia/spatial_model/calibration/slurm/experiments_summary.py`. Cells are model value "
         "(standardized residual against the data, scale in the header); † marks a moment in that case's "
         "criterion. District data: " + (a.scheme or "locale") + " partition.", ""]
    if missing:
        L += [f"Cases without results: {', '.join(missing)}.", ""]

    cols = theta_columns(rows)
    L += ["## Fit and parameters", "",
          "| case | criterion (diag) | weighted | LM its | " + " | ".join(cols) + " |",
          "|" + "---|" * (4 + len(cols))]
    for n, r in rows:
        cw = f"{r['critw']:.3g}" if r["critw"] is not None else ""
        L.append(f"| {n} | {r['crit']:.3g} | {cw} | {r['iters']} | " +
                 " | ".join(f"{r['theta'][c]:.4f}" if c in r["theta"] else "0.0000" if c in NAMES[8:] else "—"
                            for c in cols) + " |")
    extra = sorted({k for _, r in rows for k in r["active"]} - set(TARGETS))
    L += ["", "Targets: model value (standardized residual against the case's own target value and scale; "
          "† = in the criterion). Extra columns are validation moments some case moved into its criterion.", "",
          "| case | " + " | ".join(TARGETS + extra) + " |", "|" + "---|" * (1 + len(TARGETS) + len(extra))]
    for n, r in rows:
        tg = r["m"].get("targets", {})
        def tcell(k):
            if k not in r["m"]:
                return "—"
            x, t = r["m"][k], tg.get(k)
            res = f" ({(x - t['value']) / t['se']:+.2f})" if t else ""
            return f"{x:.4f}" + ("†" if k in r["active"] else "") + res
        L.append(f"| {n} | " + " | ".join(tcell(k) for k in TARGETS + extra) + " |")

    L += ["", "## Validation: district and SEDA gaps", "",
          "| case | " + " | ".join(f"{lab} [{v:+.4f}; {s:.4f}]" for lab, v, s in dist.values()) + " | ΔlogQ |",
          "|" + "---|" * (2 + len(dist))]
    for n, r in rows:
        L.append(f"| {n} | " + " | ".join(cell(r, k, v, s, 4, own=False) for k, (lab, v, s) in dist.items()) +
                 f" | {r['m'].get('gap_logQ', float('nan')):.4f} |")

    L += ["", "## Validation: residence and workplace (the split)", "",
          "Teachers' crossing, net commuting and staffing against the ACS commuting data (`commuting_all.md`, "
          "`teacher_residence.md`): the share of workers living in one of city and suburb and working in the "
          "other, the log gap of teachers living in the suburb over the city less the log gap of teachers working "
          "there, and public K–12 teachers per employed resident (data: teacher_residence_targets.toml [k12]; the "
          "national PUMS net gap +0.108 is a cross-check). The teacher supply elasticity is d log(n_work,2/n_work,1) / "
          "d log(κ₂/κ₁) by a one-sided household re-solve holding the aggregates fixed (the literature's 3.7 and 1.4 are "
          "single-district estimates, so upper bounds for a location). A dash is a bundled model (these moments do not "
          "exist without a workplace) or a moments file from before the keys. The last column is the inclusive-value gap "
          "I_r − W_rr relative to the own-location value, in utils (not the value of the workplace choice: it contains the "
          "position normalisation), for teachers, non-teachers and their difference.", "",
          "| case | " + " | ".join(f"{lab} [{v:+.3f}; {s:.3f}]" for lab, v, s in SPLIT.values()) +
          " | crossing, non-teachers | payroll/pupil gap | inclusive-value gap (teachers, non-teachers, difference) |",
          "|" + "---|" * (4 + len(SPLIT))]
    for n, r in rows:
        pp = dist["payroll_pp_gap"]
        L.append(f"| {n} | " + " | ".join(cell(r, k, v, s, 3, own=False) for k, (lab, v, s) in SPLIT.items()) +
                 f" | {level(r, ['commute_share_nonteach'], '.3f')} | {cell(r, 'payroll_pp_gap', pp[1], pp[2], 4, own=False)} | "
                 f"{level(r, ['option_value_teach', 'option_value_nonteach', 'option_value_gap'], '+.4f') if is_split(r) else '—'} |")

    L += ["", "## Validation: NLSY composition and flows by parental-income tercile", "",
          "Data (NLSY79): children in the city 0.458, 0.360, 0.163; city→suburb 0.284, 0.429, 0.407; "
          "suburb→city 0.168, 0.112, 0.168; move rate 0.221, 0.226, 0.207; overall 1→2 0.361, 2→1 0.153.", "",
          "| case | city share gradient [−0.295; 0.037] | children in 1 | 1→2 | 2→1 | move rate | P(1→2), P(2→1) | adults in 2 by own income | by own z: gradient [−0.222; 0.036] | children in 1 by own z |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for n, r in rows:
        lab, v, s = OTHER["city_share_gradient"]
        p12, p21 = r["m"].get("p12"), r["m"].get("p21")
        L.append(f"| {n} | {cell(r, 'city_share_gradient', v, s)} | {vec(r, 'city_share_by_parent_income')} | "
                 f"{vec(r, 'c2s_by_parent_income')} | {vec(r, 's2c_by_parent_income')} | "
                 f"{vec(r, 'move_rate_by_parent_income')} | {p12:.3f}, {p21:.3f} | {vec(r, 'suburb_share_by_own_income')} | "
                 f"{cell(r, 'city_share_gradient_z', -0.222, 0.036)} | {vec(r, 'city_share_by_own_z')} |")

    keys = ["teacher_premium_m", "teacher_premium_f", "s_ratio", "teacher_z_gap_market_sd", "gender_gap_market",
            "wage_slope", "rho_pc"]
    L += ["", "## Validation: teacher pay and schooling, wages, first-stage consistency", "",
          "| case | " + " | ".join(f"{OTHER[k][0]} [{OTHER[k][1]:+.4f}; {OTHER[k][2]:.3f}]" for k in keys) +
          " | pooled 90/10 [3.737] |", "|" + "---|" * (2 + len(keys))]
    for n, r in rows:
        p = f"{r['p9010']:.3f}" if r["p9010"] else "—"
        L.append(f"| {n} | " + " | ".join(cell(r, k, OTHER[k][1], OTHER[k][2]) for k in keys) + f" | {p} |")

    m3 = [k for k in M3 if any(k in r["m"] for _, r in rows)]
    if any("ba_share" in r["m"] for _, r in rows) and m3:
        L += ["", "## Education block (M3: joint BA and occupational choice)", "",
              "Data: ACS 2009–13 (BA shares, premia, schooling), NLSY79 (z gaps, wage slopes; slope holding BA fixed: "
              "convention 0.161, EIV 0.169), NCES (master's premium). Teacher premia are log hourly wages; z gaps are in SDs "
              "of log z. BA premium within occupation × gender (raw 0.481). P(BA) by own-z tercile, latent at the ACS "
              "share: 0.052, 0.262, 0.704.", "",
              "| case | BA share [all 0.339] | " + " | ".join(f"{M3[k][0]} [{M3[k][1]:+.4g}; {M3[k][2]:.3g}]" for k in m3) +
              " | P(BA) by z tercile |", "|" + "---|" * (3 + len(m3))]
        for n, r in rows:
            ba = r["m"].get("ba_share")
            L.append(f"| {n} | " + (f"{ba:.4f}" if ba is not None else "—") + " | " +
                     " | ".join(cell(r, k, M3[k][1], M3[k][2], 4) for k in m3) + f" | {vec(r, 'ba_share_by_z')} |")

    L += ["", "## Composition diagnostics (model only)", "",
          "| case | young Δ mean log z (SD) | adults Δ mean log z (SD) | score gap | score gap without Q | "
          "mean-log wage gap (market) | Q-adjusted | mover − stayer (born 1, 2) | teaching share |",
          "|---|---|---|---|---|---|---|---|---|"]
    for n, r in rows:
        m = r["m"]
        g = lambda k, d=3: f"{m[k]:+.{d}f}" if k in m else "—"
        L.append(f"| {n} | {g('dz_young_sd')} | {g('dz_adult_sd')} | {g('score_gap', 4)} | {g('score_gap_noQ', 4)} | "
                 f"{g('gap_meanlog_mkt', 4)} | {g('gap_meanlog_mkt_Qadj', 4)} | {vec(r, 'mover_premium')} | {g('teach_share', 4)} |")

    lev = [(h, k, f) for h, k, f in LEVELS if any(all(x in r["m"] for x in k) for _, r in rows)]
    L += ["", "## Model levels (model only)", "",
          "Population shares (data: 0.579 of the partitioned pupils in the suburb), tax rates (`t` is on all income, "
          "home-production output included; `t on market income` is teacher salaries over market labour income, "
          "like-for-like data ≈ 0.016: the §1 0.023 includes benefits), school quality, staffing, home production, "
          "schooling, flows, mean consumption and the BFM Table 8 WTP income elasticity (data 0.379, SE ≈ 0.17 including the SE of mean WTP; "
          "the within-z version holds ability fixed as BFM hold education fixed, and is about 1 under log "
          "utility). ΔΛ is the utility value of the suburb's school quality to the average parent: "
          "all of it, and the part due to Q alone (the rest is chiefly the Ξ offset in children's schooling). The "
          "last column is the consumption-floor audit. Students per teacher, Q and mean C are model levels "
          "without a data counterpart. A dash is a case whose moments file predates the key.", "",
          "| case | " + " | ".join(h for h, _, _ in lev) + " |", "|" + "---|" * (1 + len(lev))]
    for n, r in rows:
        L.append(f"| {n} | " + " | ".join(level(r, k, f) for _, k, f in lev) + " |")

    L += ["", "## Case options", ""]
    for n, r in rows:
        L.append(f"- **{n}**: {r['options']}")
    Path(a.out).write_text("\n".join(L) + "\n")
    print(f"wrote {a.out} ({len(rows)} cases)")


if __name__ == "__main__":
    main()
