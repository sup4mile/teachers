# =============================================================================
# TikTak estimation of the spatial block (julia/spatial_model/calibration/spatial_calibration.md §4).
#
# Fits ϑ = (log κ̃, δκ, ΔB, β, λ, r_m, ωf, γ) to the active targets of `default_targets()`
# by a TikTak multistart search over Distributed workers. The library code is in
# spatial_calibrate.jl; this file parses options, starts workers, loads the frozen
# reference levels C̄ and h̄, and runs `TikTak.minimize`.
#
# It runs in its own environment, estimation/, which pins the model packages to
# the versions in julia/Manifest.toml and adds TikTak by relative path. (The julia/
# environment cannot load TikTak: Plots and PyCall hold JSON at 0.21, TikTak needs 1.x.)
# The library files still load from julia/ for interactive work. Once per machine:
#   julia --project=julia/spatial_model/calibration/estimation -e 'using Pkg; Pkg.instantiate()'
#
# With E = julia/spatial_model/calibration/estimation and S = julia/spatial_model/calibration/spatial_estimate.jl:
#   julia --project=E S --calibrate-phi
#   julia --project=E S --freeze-reference
#   julia --project=E S --workers 8
#   julia --project=E S --run-dir <dir> --resume --max-evals 8000
#   julia --project=E S --fix r_m=0.10
#   julia --project=E S --external theta_nu=0.20 --measure chi=0.10
#   julia --project=E S --smoke --workers 2
#
# Options (defaults in brackets)
#   --workers N            local worker processes [SLURM_CPUS_PER_TASK − 1, else 0 = serial]
#   --run-dir DIR          TikTak run directory [runs/<problem id>]
#   --resume               continue the run in --run-dir (same problem id required)
#   --max-evals N          hard cap on model calls across resumes [4000]
#   --max-seconds S        soft deadline for this invocation; set it below the SLURM wall time
#   --n-samples N          Sobol screening points [256];  --n-local N  retained seeds [10%]
#   --local-max-evals N    model calls per local search [200];  --seed N  Sobol seed [0]
#   --fix name=value       hold a parameter of ϑ fixed (repeatable): r_m scenarios and profiles
#   --external name=value  override the external block (repeatable): the §4 sensitivity panel
#   --measure name=value   override a `Measurement` field (repeatable), e.g. chi=0.10, K=3,
#                          income=log_median_market
#   --wtp-consumption X    mean monthly consumption, 1990 dollars, for the BFM target [3,309]
#   --target key=value     override a target's value (repeatable; NaN drops it); --target-se key=value its scale
#   --activate key         put a validation target in the criterion (repeatable); --deactivate key the reverse
#   --reference FILE       frozen C̄ and h̄ [spatial_reference.toml]
#   --freeze-reference     (re)compute the reference levels at the baseline, save, and exit
#   --calibrate-phi        report the φ that matches non-teachers' schooling at the reference
#                          (copy it into EXTERNAL.φ, then --freeze-reference), and exit
#   --problem-id ID        override the fingerprint-based TikTak problem id
#   --warm-start SPEC      extra warm starts besides ϑ₀ (repeatable): a TikTak run directory
#                          (its 20 best points, re-evaluated) or a theta.toml from spatial_diagnose.jl
#   --verbose              log every model evaluation
#   --smoke                coarse grid, tiny budget, temporary directory: a pipeline check
# Names accept ASCII aliases: logkappa dkappa dB beta lambda r_m sigma theta_nu psi eta
# phi mu gamma alpha sigma_eps rho_z s_z omega_f chi. ωf and γ are in ϑ: hold them with
# --fix, not --external. The external block has 21
# occupations (spatial_calibrate.jl §1); --external sigma_eps=… re-inverts the
# occupational shares at the new dispersion.
#
# SLURM. One node: request --cpus-per-task=N and let --workers default to N − 1.
# Several nodes: replace `addprocs(n; ...)` below with ClusterManagers.SlurmManager
# (see TikTak.jl/README.md). Freeze the reference once before launching several
# runs, e.g. an r_m profile:
#   for r in 0 0.05 0.10 0.1813 0.30; do sbatch job.sh --fix r_m=$r; done
# All runs then share the units of C̄. Sensitivity runs (--external, --measure)
# keep the baseline reference, as §4 requires.
# =============================================================================

using Distributed
using Printf

const HERE = @__DIR__

Base.find_package("TikTak") === nothing &&
    error("TikTak is not in the active environment $(Base.active_project()). Run with " *
          "--project=$(joinpath(HERE, "estimation")) after Pkg.instantiate() there.")

include(joinpath(HERE, "spatial_cli.jl"))

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "workers" => default_workers(), "run-dir" => nothing, "resume" => false, "max-evals" => 4000,
        "max-seconds" => nothing, "n-samples" => 256, "n-local" => nothing, "local-max-evals" => 200,
        "seed" => 0, "wtp-consumption" => nothing, "reference" => nothing, "freeze-reference" => false,
        "calibrate-phi" => false, "problem-id" => nothing, "verbose" => false, "smoke" => false);
    flags   = ("resume", "freeze-reference", "calibrate-phi", "verbose", "smoke"),
    ints    = ("workers", "max-evals", "n-samples", "n-local", "local-max-evals", "seed"),
    floats  = ("max-seconds", "wtp-consumption"),
    strings = ("run-dir", "reference", "problem-id"),
    lists   = ("fix", "external", "measure", "target", "target-se"),
    multi   = ("warm-start", "activate", "deactivate"),
    script  = "spatial_estimate.jl")
OPTS["workers"] > 0 && !OPTS["freeze-reference"] && !OPTS["calibrate-phi"] &&
    addprocs(OPTS["workers"]; exeflags = "--project=$(Base.active_project())")

@everywhere begin
    using TikTak, LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($HERE, "spatial_calibrate.jl"))
end

const SMOKE_SOLVER = merge(SOLVER, (; Nz = 3, nϵT = 24, nXO = 24))

function progress(r)
    @printf("  [%s] local searches %d, model calls %d (%d failed points), best criterion %.6g\n",
            Libc.strftime("%H:%M:%S", time()), r.n_local_completed, r.n_evals, r.n_failed, r.fun)
    flush(stdout)
end

function report(io::IO, result, model, targets, opts, pid)
    println(io, "problem id: ", pid)
    println(io, "options: ", join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                                   if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "))
    println(io, result)
    println(io, result.message)
    if has_solution(result)
        n_conv = count(r -> r.converged, result.local_results)
        @printf(io, "local searches: %d completed, %d report convergence\n", length(result.local_results), n_conv)
        evaluate(model, result.x, targets; io)
    end
end

function main(opts)
    smoke  = opts["smoke"]
    solver = smoke ? SMOKE_SOLVER : SOLVER
    if opts["freeze-reference"]
        freeze_reference(something(opts["reference"], REFERENCE_FILE), solver)
        return
    end
    if opts["calibrate-phi"]
        path = something(opts["reference"], REFERENCE_FILE)
        ext0 = external_params()
        isfile(path) && (old = load_reference(path); ext0 = with_params(ext0; Cbar = old.Cbar, href = old.href))
        r = calibrate_phi(ext0, THETA0; solver)
        @printf("\nφ = %.6f gives mean s_O = %.6f (target %.6f) at ϑ₀; C̄ = %.8f, h̄ = %.8f.\n",
                r.φ, r.sO, schooling_target(), r.Cbar, r.href)
        println("Set EXTERNAL.φ to this value (4 decimals) in spatial_calibrate.jl, then run --freeze-reference.")
        return
    end

    run_dir = opts["run-dir"]
    smoke && run_dir === nothing && (run_dir = mktempdir(; prefix = "spatial-smoke-", cleanup = false))
    ref = smoke && opts["reference"] === nothing ?
          freeze_reference(joinpath(run_dir, "reference.toml"), solver) :
          frozen_reference(something(opts["reference"], REFERENCE_FILE), solver)

    (; base, meas, targets, active, bounds, θ0) = problem_setup(opts, ref)
    for t in setdiff(targets, active)
        println(t.active ? "  ⚠ target $(t.key) is pending" : "  target $(t.key) is validation",
                " and left out of the criterion: $(t.note)")
    end

    model     = SpatialModel(base, [t.key for t in active]; meas, solver, verbose = opts["verbose"])
    objective = MomentObjective(model, [t.value for t in active]; scales = [t.se for t in active])
    pid       = something(opts["problem-id"], "spatial-" * problem_fingerprint(model, active, bounds))
    run_dir   = something(run_dir, joinpath(HERE, "runs", pid))
    config    = smoke ?
        TikTakConfig(n_samples = 8, n_local = 2, local_max_evals = 6, max_evals = 24, x_tol = 1e-3, f_tol = 1e-4,
                     failure_exceptions = (ModelEvaluationError, CalibrationFailure, DomainError)) :
        TikTakConfig(n_samples = opts["n-samples"], n_local = opts["n-local"], max_evals = opts["max-evals"],
                     local_max_evals = opts["local-max-evals"], seed = opts["seed"], max_seconds = opts["max-seconds"],
                     x_tol = 1e-4, f_tol = 1e-6,
                     failure_exceptions = (ModelEvaluationError, CalibrationFailure, DomainError))

    println("\nSpatial calibration: problem $pid, $(nworkers() == 1 && workers() == [1] ? "serial" : "$(nworkers()) workers")")
    println("  run directory: ", run_dir)
    println("  free: ", join((n for (n, (lo, hi)) in zip(THETA_NAMES, bounds) if lo < hi), ", "),
            isempty(opts["fix"]) ? "" : ";  fixed: " * join(("$k = $v" for (k, v) in opts["fix"]), ", "))
    isempty(opts["external"]) || println("  external overrides: ", join(("$k = $v" for (k, v) in opts["external"]), ", "))
    println("  measurement: ", meas)
    println("  grid: ", solver)
    flush(stdout)

    starts = [θ0]
    for path in opts["warm-start"]
        pts  = isdir(path) && !isfile(joinpath(path, "theta.toml")) ? load_estimates(path; limit = 20) :
               [Float64.(TOML.parsefile(isdir(path) ? joinpath(path, "theta.toml") : path)["theta"])]
        for x in pts
            length(x) == length(θ0) || error("--warm-start $path: $(length(x)) values, ϑ has $(length(θ0))")
            for (name, v) in opts["fix"]
                x[theta_index(name)] = v
            end
            push!(starts, clamp.(x, first.(bounds), last.(bounds)))
        end
    end
    length(starts) > 1 && println("  warm starts: ϑ₀ and $(length(starts) - 1) more")
    t0 = time()
    result = minimize(objective, bounds; config, run_dir, problem_id = pid, resume = opts["resume"],
                      warm_start = opts["resume"] ? nothing : starts, callback = progress)
    @printf("\nTikTak finished in %.1f minutes.\n", (time() - t0) / 60)

    report(stdout, result, model, targets, opts, pid)
    open(joinpath(run_dir, "fit_summary.txt"), "w") do io
        report(io, result, model, targets, opts, pid)
    end
    println("\nSummary written to ", joinpath(run_dir, "fit_summary.txt"))
    return result
end

main(OPTS)
