# =============================================================================
# Validation outcomes with teacher composition held fixed (calibration experiment 6,
# 2026-10; julia/spatial_model/calibration/calibration_log.md).
#
# At a fitted ϑ, solve the general equilibrium and then re-solve the household block in
# partial equilibrium with the teacher-quality index Q_T,l replaced, holding taxes, the
# teacher stocks and, with peer effects (ψa ≠ 0), the peer factor ā_l^ψa of the solution
# fixed (school quality is Q_l = Q_T,l ā_l^ψa):
#   q_base   Q_T,l of a reference solution (--q-from, default the exactly identified
#            baseline, runs/exp-final/p0_baseline of the repository, also when run from a code
#            snapshot; its moments.toml "QT", or "Q" for runs that predate peer effects): teacher
#            composition as in the baseline, so what moves the validation outcomes is the
#            residential side of the specification;
#   q_equal  Q_T,1 = Q_T,2 = (2 ΣH̃_T / ΣM)^σ: no teacher-quality gap, so the change from
#            the equilibrium is the teacher channel's contribution.
# Households, location choices and the stationary distribution are re-solved; the
# reported outcomes are the residential ones (income, poverty, scores, composition),
# not the teacher moments, which a fixed Q leaves inconsistent.
#
# With E = julia/spatial_model/calibration/estimation and C = julia/spatial_model/calibration/spatial_counterfactual.jl:
#   julia --project=E C --theta julia/spatial_model/calibration/runs/exp-sorting/taste_income \
#         --free zeta --out julia/spatial_model/calibration/runs/exp-counterfactual/taste_income
# Options: --theta (a run directory or theta.toml, or v1,v2,…), --q-from DIR (a run
# directory with moments.toml), --out DIR, and the problem options of spatial_diagnose.jl
# (--external, --measure, --fix, --floor-tol, …; --fix overrides the read ϑ; --floor-tol is the
# largest population share at the consumption floor, for the equilibrium and the re-solves alike).
# =============================================================================

using Printf
using TOML

const HERE = @__DIR__
include(joinpath(HERE, "spatial_cli.jl"))
include(joinpath(HERE, "spatial_calibrate.jl"))

"""
The repository root: the parent of `data/`, reached as `DISTRICT_SETTINGS` reaches the data (three levels up
from this directory). A code snapshot (runs/_code/<tag>/) links `data` to the repository's own, so the root
resolved through the link is the repository's whatever copy of the code runs.
"""
const REPO_ROOT = dirname(realpath(joinpath(HERE, "..", "..", "..", "data")))

"The default `--q-from`: the exactly identified baseline of the repository's runs, also from a snapshot."
const DEFAULT_Q_FROM = joinpath(REPO_ROOT, "julia", "spatial_model", "calibration", "runs", "exp-final", "p0_baseline")

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "theta" => "", "q-from" => DEFAULT_Q_FROM, "out" => nothing,
        "reference" => nothing, "wtp-consumption" => nothing, "district" => nothing, "weight" => "diag",
        "floor-tol" => nothing);
    floats  = ("wtp-consumption", "floor-tol"),
    strings = ("theta", "q-from", "out", "reference", "district", "weight"),
    lists   = ("fix", "external", "measure", "target", "target-se"),
    multi   = ("activate", "deactivate", "free"),
    script  = "spatial_counterfactual.jl")

"ϑ from a run directory, a theta.toml file or a comma list (eight values are padded with the extensions off)."
function read_run_theta(spec)
    isempty(spec) && error("--theta is required: a run directory, a theta.toml file or v1,v2,…")
    isdir(spec) && return read_run_theta(joinpath(spec, "theta.toml"))
    isfile(spec) && return pad_theta(Float64.(TOML.parsefile(spec)["theta"]))
    return pad_theta(parse.(Float64, split(spec, ',')))
end

"""
    fixed_q_solution(sol, QT; hh_tol, hh_maxit) -> sol′

The household block, location choices and stationary distribution re-solved at the
aggregates of `sol` (taxes, teacher stocks, peer ability) with the teacher-quality index
replaced by `QT`, so that school quality is `QT` times the solution's peer factor; under
the split (`sol.p.split`) the teachers' positions stay at the pupils' shares of `sol`.
"""
function fixed_q_solution(sol, QT; hh_tol = SOLVER.hh_tol, hh_maxit = SOLVER.hh_maxit)
    (; p, gr) = sol
    QT  = collect(Float64, QT)
    gr2 = merge(gr, (; QT, Q = iszero(p.ψa) ? copy(QT) : QT .* exp.(p.ψa .* gr.la)))
    hh  = solve_household(p, gr2; tol = hh_tol, maxit = hh_maxit, hh0 = sol.hh)
    p.split && return merge(sol, (; hh, split_probs(hh, p, gr2)..., gr = gr2, s_resid = hh.s_resid,
                                    s_passes = hh.s_passes, s_relaxed = hh.s_relaxed))
    Pi_T, Pi_O = location_probs(hh, p, gr2)
    Φ, πbar = stationary_phi(Pi_T, Pi_O, hh, p, gr2)
    return merge(sol, (; hh, Pi_T, Pi_O, Φ, πbar, gr = gr2, s_resid = hh.s_resid, s_passes = hh.s_passes,
                         s_relaxed = hh.s_relaxed))
end

const ROWS = [(:gap_income, "earnings gap (median, market)", 0.1122), (:gap_meanlog_mkt, "mean-log wage gap (market)", NaN),
              (:gap_meanlog_mkt_Qadj, "  Q-adjusted", NaN), (:poverty_gap, "child-poverty gap", -0.0888),
              (:score_gap, "SEDA score gap", 0.2887), (:score_gap_noQ, "  without the Q term", NaN),
              (:seda_growth_gap, "SEDA growth gap", 0.0031), (:city_share_gradient, "city share top − bottom", -0.295),
              (:dz_young_sd, "young Δ mean log z (SD)", NaN), (:dz_adult_sd, "adults Δ mean log z (SD)", NaN),
              (:gap_pupils, "enrollment gap", 0.3764), (:move_rate, "move rate", 0.2156), (:gap_logQ, "Δ log Q", NaN)]

function main(opts)
    ref = frozen_reference(something(opts["reference"], REFERENCE_FILE), SOLVER)
    (; base, meas, floor_tol) = problem_setup(opts, ref)
    θ = read_run_theta(opts["theta"])
    for (name, v) in opts["fix"]
        θ[theta_index(name)] = v
    end
    model = SpatialModel(base, Symbol[]; meas, floor_tol)
    sol, m = solve_theta(model, θ)
    qm  = TOML.parsefile(joinpath(opts["q-from"], "moments.toml"))
    Qb  = Float64.(get(qm, "QT", qm["Q"]))       # the teacher index; Q itself before peer effects
    Qeq = fill((2 * sum(sol.HT) / sum(sol.M))^sol.p.σ, 2)
    pe = [("Q at reference", fixed_q_solution(sol, Qb)), ("Q equal", fixed_q_solution(sol, Qeq))]
    cases = vcat([("equilibrium", m)], [(name, model_moments(s, meas)) for (name, s) in pe])
    # the partial-equilibrium re-solves get the feasibility checks `solve_theta` applies in GE
    checks = [@sprintf("%s: (e, s) residual %.1e, floor mass %.1e (no feasible location %.1e)%s", name, s.s_resid,
                       c[2].floor_mass, c[2].infeasible_mass,
                       (s.s_resid ≤ 1e-8 && c[2].floor_mass + c[2].infeasible_mass ≤ model.floor_tol) ? "" : "  ⚠ FAILS the solve_theta checks")
              for ((name, s), c) in zip(pe, cases[2:end])]
    out = opts["out"]
    io  = out === nothing ? stdout : (mkpath(out); open(joinpath(out, "counterfactual.txt"), "w"))
    @printf(io, "spatial_counterfactual.jl  %s\n  theta %s\n  ϑ = %s\n  Q: equilibrium %s, reference %s (%s), equal %s\n",
            Libc.strftime("%Y-%m-%d %H:%M:%S", time()), opts["theta"], fmtvec(θ), fmtvec(sol.gr.QT), fmtvec(Qb),
            opts["q-from"], fmtvec(Qeq))
    iszero(sol.p.ψa) || @printf(io, "  (teacher index Q_T; peers held at the equilibrium: ψa = %.3f, log ā = %s)\n",
                                sol.p.ψa, fmtvec(sol.gr.la))
    println(io)
    foreach(c -> println(io, "  ", c), checks)
    @printf(io, "\n  %-30s %9s %15s %15s %15s\n", "outcome", "data (locale)", cases[1][1], cases[2][1], cases[3][1])
    for (k, lab, d) in ROWS
        @printf(io, "  %-30s %9s %15.4f %15.4f %15.4f\n", lab, isnan(d) ? "" : @sprintf("%.4f", d),
                (c[2][k] for c in cases)...)
    end
    for (k, lab) in ((:city_share_by_parent_income, "children in 1 by parental tercile"),
                     (:c2s_by_parent_income, "1→2 by parental tercile"), (:s2c_by_parent_income, "2→1 by parental tercile"))
        @printf(io, "  %-34s %s\n", lab, join((fmtvec(c[2][k]) for c in cases), "  "))
    end
    if out !== nothing
        d = Dict{String,Any}("theta" => θ, "Q_equilibrium" => sol.gr.QT, "Q_reference" => Qb, "Q_equal" => Qeq)
        for (name, mm) in cases, (k, _, _) in ROWS
            d[replace(name, " " => "_") * "." * string(k)] = Float64(mm[k])
        end
        open(f -> TOML.print(f, d; sorted = true), joinpath(out, "counterfactual.toml"), "w")
        close(io)
        println("written to ", out)
    end
end

main(OPTS)
