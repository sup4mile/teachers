# =============================================================================
# TikTak estimation of the spatial block (data/spatial/spatial_calibration.md §4).
#
# Fits ϑ = (log κ̄, δκ, ΔB, β, λ, r_m) to the active targets of `default_targets()`
# by a TikTak multistart search over Distributed workers. The library code is in
# spatial_calibrate.jl; this file parses options, starts workers, loads the frozen
# reference levels C̄ and h̄, and runs `TikTak.minimize`.
#
# It runs in its own environment, estimation/, which pins the model packages to
# the versions in julia/Manifest.toml and adds TikTak by relative path. (The julia/
# environment cannot load TikTak: Plots and PyCall hold JSON at 0.21, TikTak needs 1.x.)
# The library files still load from julia/ for interactive work. Once per machine:
#   julia --project=julia/spatial_model/estimation -e 'using Pkg; Pkg.instantiate()'
#
# With E = julia/spatial_model/estimation and S = julia/spatial_model/spatial_estimate.jl:
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
#   --measure name=value   override a `Measurement` field (repeatable), e.g. chi=0.10, K=3
#   --wtp-consumption X    mean monthly consumption, 1990 dollars; activates the BFM target
#   --reference FILE       frozen C̄ and h̄ [spatial_reference.toml]
#   --freeze-reference     (re)compute the reference levels at the baseline, save, and exit
#   --problem-id ID        override the fingerprint-based TikTak problem id
#   --verbose              log every model evaluation
#   --smoke                coarse grid, tiny budget, temporary directory: a pipeline check
# Names accept ASCII aliases: logkappa dkappa dB beta lambda r_m sigma theta_nu psi eta
# phi mu gamma alpha sigma_eps rho_z sigma_xi chi.
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

const ALIASES = Dict("logkappa" => "logκ", "dkappa" => "δκ", "dB" => "ΔB", "beta" => "β",
                     "lambda" => "λ", "sigma" => "σ", "theta_nu" => "θν", "psi" => "ψ",
                     "eta" => "η", "phi" => "φ", "mu" => "μ", "gamma" => "γ", "alpha" => "α",
                     "sigma_eps" => "σϵ", "rho_z" => "ρz", "sigma_xi" => "σξ", "chi" => "χ")

"`name=value` → Symbol => value; the value is a Float64 when it parses as one, else a Symbol."
function parse_assignment(s)
    occursin('=', s) || error("expected name=value, got \"$s\"")
    k, v = split(s, '='; limit = 2)
    return Symbol(get(ALIASES, k, k)) => something(tryparse(Float64, v), Symbol(v))
end

function parse_options(args)
    o = Dict{String,Any}(
        "workers" => max(parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1")) - 1, 0),
        "run-dir" => nothing, "resume" => false, "max-evals" => 4000, "max-seconds" => nothing,
        "n-samples" => 256, "n-local" => nothing, "local-max-evals" => 200, "seed" => 0,
        "fix" => Pair{Symbol,Any}[], "external" => Pair{Symbol,Any}[], "measure" => Pair{Symbol,Any}[],
        "wtp-consumption" => nothing, "reference" => nothing, "freeze-reference" => false,
        "problem-id" => nothing, "verbose" => false, "smoke" => false)
    flags   = ("resume", "freeze-reference", "verbose", "smoke")
    ints    = ("workers", "max-evals", "n-samples", "n-local", "local-max-evals", "seed")
    floats  = ("max-seconds", "wtp-consumption")
    strings = ("run-dir", "reference", "problem-id")
    lists   = ("fix", "external", "measure")
    i = 1
    while i ≤ length(args)
        key = lstrip(args[i], '-')
        if key in flags
            o[key] = true
            i += 1
            continue
        end
        i < length(args) || error("option --$key needs a value")
        v = args[i+1]
        key in ints    ? (o[key] = parse(Int, v)) :
        key in floats  ? (o[key] = parse(Float64, v)) :
        key in strings ? (o[key] = String(v)) :
        key in lists   ? push!(o[key], parse_assignment(v)) :
        error("unknown option --$key (see the header of spatial_estimate.jl)")
        i += 2
    end
    for (name, v) in o["fix"]
        v isa Float64 || error("--fix $name needs a number")
    end
    return o
end

const OPTS = parse_options(ARGS)
OPTS["workers"] > 0 && !OPTS["freeze-reference"] &&
    addprocs(OPTS["workers"]; exeflags = "--project=$(Base.active_project())")

@everywhere begin
    using TikTak, LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($HERE, "spatial_calibrate.jl"))
end

const SMOKE_SOLVER = merge(SOLVER, (; Nz = 3, nϵT = 24, nXO = 24))

"Compute the reference levels at the baseline external block and ϑ₀, and save them."
function freeze_reference(path, solver)
    println("Freezing C̄ and h̄ at the baseline external block and ϑ₀ → $path")
    ext0 = external_params()
    ref  = reference_levels(ext0, THETA0; solver)
    save_reference(path, ref, ext0; θ = THETA0, solver)
    return load_reference(path)
end

"Load the frozen reference (freezing it first if the file is missing) and flag configuration drift."
function frozen_reference(path, solver)
    ref = isfile(path) ? load_reference(path) : freeze_reference(path, solver)
    drift = reference_mismatch(ref, external_params(), solver)
    @printf("Reference levels (%s): C̄ = %.8f, h̄ = %.8f, frozen %s\n", path, ref.Cbar, ref.href,
            get(ref.meta, "created", "?"))
    isempty(drift) || println("  ⚠ measured at a different baseline or grid: ", join(drift, "; "),
                              ". Rerun --freeze-reference unless this is deliberate.")
    return ref
end

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

    run_dir = opts["run-dir"]
    smoke && run_dir === nothing && (run_dir = mktempdir(; prefix = "spatial-smoke-", cleanup = false))
    ref = smoke && opts["reference"] === nothing ?
          freeze_reference(joinpath(run_dir, "reference.toml"), solver) :
          frozen_reference(something(opts["reference"], REFERENCE_FILE), solver)

    base    = calibration_base(ref; opts["external"]...)
    meas    = Measurement(; opts["measure"]...)
    targets = default_targets()
    opts["wtp-consumption"] === nothing || (targets = with_wtp_target(targets, opts["wtp-consumption"]))
    active  = active_targets(targets)
    for t in setdiff(targets, active)
        println("  ⚠ target $(t.key) is pending and left out of the criterion: $(t.note)")
    end

    bounds = copy(THETA_BOUNDS)
    θ0     = copy(THETA0)
    for (name, v) in opts["fix"]
        j = theta_index(name)
        bounds[j] = (v, v)
        θ0[j] = v
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

    t0 = time()
    result = minimize(objective, bounds; config, run_dir, problem_id = pid, resume = opts["resume"],
                      warm_start = opts["resume"] ? nothing : [θ0], callback = progress)
    @printf("\nTikTak finished in %.1f minutes.\n", (time() - t0) / 60)

    report(stdout, result, model, targets, opts, pid)
    open(joinpath(run_dir, "fit_summary.txt"), "w") do io
        report(io, result, model, targets, opts, pid)
    end
    println("\nSummary written to ", joinpath(run_dir, "fit_summary.txt"))
    return result
end

main(OPTS)
