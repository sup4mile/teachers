"""Bin scatters of model moments against ϑ over the Sobol screening points.

Reads one or more screen.jsonl files written by spatial_screen.jl and draws, in
<out>/:
  targeted.png     the 8 targeted moments against the 8 parameters
  validation.png   Table 4 validation moments and diagnostics against the 8 parameters
  params_best.png  pairwise scatter of ϑ among the points with the lowest criterion
  failures.png     share of failed solves by parameter decile
Each moment panel shows every solved point (gray), the lowest-criterion 10% (blue),
decile-bin means (black), the data target (dashed) and the baseline estimate (dotted).

    python3 julia/spatial_model/calibration/screen_plots.py julia/spatial_model/calibration/runs/screen*/screen.jsonl \
        --out julia/spatial_model/calibration/runs/screen/plots
"""
import argparse
import json
import pathlib
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).parent
ESTIMATE = HERE / "runs" / "exact-polish-base" / "theta.toml"

PARAM_LABELS = [r"$\log\tilde\kappa$", r"$\delta_\kappa$", r"$\Delta B$", r"$\beta$",
                r"$\lambda$", r"$r_m$", r"$\omega_f$", r"$\gamma$"]
TARGETED = [("male_teach_share", "Male teaching share"), ("female_teach_share", "Female teaching share"),
            ("teacher_p9010", "Teacher 90/10"), ("gap_salary", "Salary gap"),
            ("gap_pupils", "Enrollment gap"), ("cfr_effect", "CFR effect"),
            ("wtp_share", "WTP share"), ("move_rate", "Move rate")]
# Table 4 data values; the income and FTE gaps are the inactive targets in default_targets().
VALIDATION = [("gap_income", "Earnings gap", 0.1122), ("gap_teachers_pp", "FTE/pupil gap", -0.0177),
              ("score_gap", "SEDA score gap", 0.2887), ("seda_growth_gap", "SEDA growth gap", 0.0031),
              ("poverty_gap", "Child-poverty gap", -0.0888), ("s_ratio", "Teacher/non-teacher\nschooling", 1.176),
              ("wage_slope", "NLSY wage slope", 0.2127), ("rho_pc", "Mother–child ρ", 0.5766)]

INK, MUTED, POINT, BEST = "#1f2328", "#6e7781", "#c4c9cf", "#2f6fdb"


def load(paths):
    header, rows = None, {}
    for p in paths:
        for line in open(p):
            r = json.loads(line)
            if r.get("header"):
                if header is not None and r["problem_id"] != header["problem_id"]:
                    raise SystemExit(f"{p}: problem {r['problem_id']} differs from {header['problem_id']}")
                header = r
            else:
                rows[r["index"]] = r  # a Sobol index evaluated twice keeps the later record
    names = header["theta_names"]
    recs = []
    for i, r in sorted(rows.items()):
        d = {"index": i, "ok": r["status"] == "ok"}
        d.update(dict(zip(names, r["theta"])))
        d.update(r["moments"] or {})
        recs.append(d)
    df = pd.DataFrame(recs)
    act = [t for t in header["targets"] if t["active"]]
    ok = df["ok"]
    df["criterion"] = np.where(ok, sum(((df[t["key"]] - t["value"]) / t["se"]) ** 2 for t in act), np.nan)
    return header, df


def binned(x, y, nbins=10):
    q = np.unique(np.quantile(x, np.linspace(0, 1, nbins + 1)))
    b = np.clip(np.searchsorted(q, x, side="right") - 1, 0, len(q) - 2)
    g = pd.DataFrame({"b": b, "x": x, "y": y}).groupby("b")
    return g["x"].mean().to_numpy(), g["y"].mean().to_numpy()


def moment_grid(df, names, bounds, rows, targets, theta_hat, path, title):
    ok = df[df["ok"]]
    best = ok[ok["criterion"] <= ok["criterion"].quantile(0.10)]
    nr, nc = len(rows), len(names)
    fig, axes = plt.subplots(nr, nc, figsize=(2.1 * nc, 1.75 * nr), sharex="col", sharey="row",
                             squeeze=False, constrained_layout=True)
    for i, (key, label) in enumerate(rows):
        y = ok[key].to_numpy()
        # Clip the y range to the central 98% so a few extreme points don't flatten the rest.
        lo, hi = np.nanquantile(y, [0.01, 0.99])
        pad = 0.08 * (hi - lo)
        for j, name in enumerate(names):
            ax = axes[i, j]
            ax.scatter(ok[name], y, s=3, color=POINT, linewidths=0, rasterized=True)
            ax.scatter(best[name], best[key], s=4, color=BEST, linewidths=0, rasterized=True)
            bx, by = binned(ok[name].to_numpy(), y)
            ax.plot(bx, by, color=INK, lw=1.4, marker="o", ms=2.5)
            if targets.get(key) is not None:
                ax.axhline(targets[key], color=INK, lw=0.9, ls="--")
            ax.axvline(theta_hat[j], color=MUTED, lw=0.8, ls=":")
            ax.set_ylim(min(lo, targets.get(key, lo)) - pad, max(hi, targets.get(key, hi)) + pad)
            ax.set_xlim(*bounds[j])
            ax.tick_params(labelsize=6, length=2, color=MUTED)
            for s in ax.spines.values():
                s.set_color(MUTED)
                s.set_linewidth(0.5)
            if i == nr - 1:
                ax.set_xlabel(PARAM_LABELS[j], fontsize=9)
            if j == 0:
                ax.set_ylabel(label, fontsize=8)
    handles = [plt.Line2D([], [], color=POINT, marker="o", ls="", ms=4, label=f"solved point (n = {len(ok)})"),
               plt.Line2D([], [], color=BEST, marker="o", ls="", ms=4, label="lowest-criterion 10%"),
               plt.Line2D([], [], color=INK, marker="o", ms=3, lw=1.4, label="decile-bin mean"),
               plt.Line2D([], [], color=INK, ls="--", lw=0.9, label="data"),
               plt.Line2D([], [], color=MUTED, ls=":", lw=0.8, label="baseline estimate")]
    fig.legend(handles=handles, loc="outside lower center", ncol=5, fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11, color=INK)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def params_best(df, names, bounds, theta_hat, path, frac=0.10):
    ok = df[df["ok"]]
    best = ok[ok["criterion"] <= ok["criterion"].quantile(frac)].sort_values("criterion", ascending=False)
    c = np.log10(best["criterion"])
    n = len(names)
    fig, axes = plt.subplots(n, n, figsize=(13, 12.5), constrained_layout=True)
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                ax.hist(ok[names[j]], bins=20, range=bounds[j], color=POINT)
                ax.hist(best[names[j]], bins=20, range=bounds[j], color=BEST)
                ax.axvline(theta_hat[j], color=INK, lw=0.8, ls=":")
                ax.set_yticks([])
            elif i > j:
                sc = ax.scatter(best[names[j]], best[names[i]], c=c, cmap="Blues_r", s=7, linewidths=0,
                                vmin=c.min(), vmax=c.max() + 0.3 * (c.max() - c.min()))
                ax.plot(theta_hat[j], theta_hat[i], marker="+", color="#d1242f", ms=9, mew=1.5)
                ax.set_ylim(*bounds[i])
            else:
                ax.axis("off")
                continue
            ax.set_xlim(*bounds[j])
            ax.tick_params(labelsize=6, length=2)
            if i == n - 1:
                ax.set_xlabel(PARAM_LABELS[j], fontsize=10)
            if j == 0 and i > 0:
                ax.set_ylabel(PARAM_LABELS[i], fontsize=10)
    fig.colorbar(sc, ax=axes[:3, -3:], shrink=0.6, label=r"$\log_{10}$ criterion")
    fig.suptitle(f"ϑ at the lowest-criterion {frac:.0%} of screening points (n = {len(best)} of {len(ok)}); "
                 "+ = baseline estimate", fontsize=11)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def failures(df, names, bounds, path):
    fig, axes = plt.subplots(2, 4, figsize=(11, 4.6), sharey=True, constrained_layout=True)
    for j, (ax, name) in enumerate(zip(axes.flat, names)):
        edges = np.linspace(*bounds[j], 11)
        b = np.clip(np.digitize(df[name], edges) - 1, 0, 9)
        rate = (~df["ok"]).groupby(b).mean().reindex(range(10), fill_value=np.nan)
        ax.bar((edges[:-1] + edges[1:]) / 2, rate, width=0.9 * np.diff(edges), color=BEST)
        ax.set_xlabel(PARAM_LABELS[j])
        ax.tick_params(labelsize=7)
    axes[0, 0].set_ylabel("share failed")
    axes[1, 0].set_ylabel("share failed")
    fig.suptitle(f"Failed solves by parameter decile ({(~df['ok']).sum()} of {len(df)} points)", fontsize=11)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", default=str(HERE / "runs" / "screen" / "plots"))
    a = ap.parse_args()
    header, df = load(a.files)
    names, bounds = header["theta_names"], [tuple(b) for b in header["bounds"]]
    theta_hat = tomllib.load(open(ESTIMATE, "rb"))["theta"]
    targets = {t["key"]: t["value"] for t in header["targets"]}
    targets.update({k: v for k, _, v in VALIDATION})
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "screen_points.csv", index=False)
    moment_grid(df, names, bounds, TARGETED, targets, theta_hat, out / "targeted.png",
                "Targeted moments over the Sobol screening points")
    moment_grid(df, names, bounds, [(k, l) for k, l, _ in VALIDATION], targets, theta_hat, out / "validation.png",
                "Validation moments and diagnostics over the Sobol screening points")
    params_best(df, names, bounds, theta_hat, out / "params_best.png")
    failures(df, names, bounds, out / "failures.png")
    ok = df[df["ok"]]
    print(f"{len(df)} points, {len(ok)} solved; criterion quantiles (10/50/90%):",
          np.round(ok["criterion"].quantile([0.1, 0.5, 0.9]).to_numpy(), 2))
    print("best 5 points:")
    print(ok.nsmallest(5, "criterion")[["index", *names, "criterion"]].round(4).to_string(index=False))
    print("figures in", out)


if __name__ == "__main__":
    main()
