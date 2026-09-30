# =============================================================================
# Sobol screening of the spatial block: model moments over the ϑ search box.
#
# Evaluates points of TikTak's scrambled Sobol sequence (the sequence `spatial_estimate.jl`
# screens, same seed) over THETA_BOUNDS and records every scalar field of
# `model_moments` at each point, validation moments and diagnostics included, not only
# the active targets TikTak journals. Output is one JSON line per point in
# <run-dir>/screen.jsonl after a header line with the targets, bounds and problem id;
# a rerun with the same run directory skips indices already recorded, so preempted
# jobs resume. Plots: screen_plots.py.
#
# With E = julia/spatial_model/calibration/estimation and S = julia/spatial_model/calibration/spatial_screen.jl:
#   julia --project=E S --run-dir runs/screen --first 513 --n 1024 --workers 63
#
# Options (defaults in brackets)
#   --run-dir DIR   output directory [runs/screen]
#   --first I       first index of the Sobol sequence to evaluate [513: after exact-tiktak's 512]
#   --n N           number of points [1024]
#   --seed N        Sobol scrambling seed [0, as spatial_estimate.jl]
#   --workers N     local worker processes [SLURM_CPUS_PER_TASK − 1]
#   --fix, --external, --measure, --reference: as in spatial_estimate.jl
# =============================================================================

using Distributed
using Printf

const HERE = @__DIR__
include(joinpath(HERE, "spatial_cli.jl"))

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "workers" => default_workers(), "run-dir" => joinpath(HERE, "runs", "screen"), "first" => 513,
        "n" => 1024, "seed" => 0, "reference" => nothing, "wtp-consumption" => nothing);
    ints    = ("workers", "first", "n", "seed"),
    floats  = ("wtp-consumption",),
    strings = ("run-dir", "reference"),
    lists   = ("fix", "external", "measure", "target", "target-se"),
    multi   = ("activate", "deactivate"),
    script  = "spatial_screen.jl")
OPTS["workers"] > 0 && addprocs(OPTS["workers"]; exeflags = "--project=$(Base.active_project())")

@everywhere begin
    using TikTak, LinearAlgebra
    using TikTak: JSON   # JSON is TikTak's dependency, not a direct one of estimation/
    BLAS.set_num_threads(1)
    include(joinpath($HERE, "spatial_calibrate.jl"))

    "Finite reals and short real vectors of `m`, flattened to name => value (vectors as name_1, name_2, …)."
    function flat_moments(m)
        out = Dict{String,Any}()
        for (k, v) in pairs(m)
            if v isa Real
                out[string(k)] = isfinite(v) ? Float64(v) : nothing
            elseif v isa AbstractVector{<:Real} && length(v) ≤ 25
                for (i, x) in enumerate(v)
                    out["$(k)_$i"] = isfinite(x) ? Float64(x) : nothing
                end
            end
        end
        return out
    end

    "Solve at ϑ; the moments, or the failure message for an expected failure."
    function screen_point(model, θ)
        t0 = time()
        try
            _, m = solve_theta(model, θ)
            return (; status = "ok", seconds = time() - t0, error = nothing, moments = flat_moments(m))
        catch e
            e isa Union{CalibrationFailure,DomainError} || rethrow()
            return (; status = "failed", seconds = time() - t0, error = sprint(showerror, e), moments = nothing)
        end
    end
end

function main(opts)
    ref = frozen_reference(something(opts["reference"], REFERENCE_FILE), SOLVER)
    (; base, meas, targets, active, bounds) = problem_setup(opts, ref)
    model = SpatialModel(base, [t.key for t in active]; meas)
    pid   = "spatial-" * problem_fingerprint(model, active, bounds)

    transform = BoxTransform(bounds)
    first_i, n = opts["first"], opts["n"]
    units = TikTak.sobol_points(transform.dimension, first_i + n - 1, opts["seed"])

    dir  = mkpath(opts["run-dir"])
    path = joinpath(dir, "screen.jsonl")
    done = Set{Int}()
    if isfile(path)
        for line in eachline(path)
            r = JSON.parse(line)
            haskey(r, "index") && push!(done, r["index"])
        end
    else
        header = Dict("header" => true, "problem_id" => pid, "seed" => opts["seed"],
                      "theta_names" => string.(THETA_NAMES), "bounds" => [collect(b) for b in bounds],
                      "targets" => [Dict("key" => string(t.key), "value" => isfinite(t.value) ? t.value : nothing,
                                         "se" => isfinite(t.se) ? t.se : nothing, "active" => t.active,
                                         "param" => string(t.param)) for t in targets],
                      "theta0" => THETA0)
        open(io -> println(io, JSON.json(header)), path, "w")
    end
    todo = [i for i in first_i:(first_i + n - 1) if i ∉ done]
    @printf("Screening problem %s: %d of %d points left (Sobol %d–%d, seed %d), %d workers\n  output %s\n",
            pid, length(todo), n, first_i, first_i + n - 1, opts["seed"], nworkers(), path)
    flush(stdout)

    pool = workers() == [1] ? nothing : WorkerPool(workers())
    lk, n_done, n_fail, t0 = ReentrantLock(), Ref(0), Ref(0), time()
    asyncmap(todo; ntasks = max(nworkers(), 1)) do i
        θ = to_parameters(transform, units[i])
        r = pool === nothing ? screen_point(model, θ) : remotecall_fetch(screen_point, pool, model, θ)
        rec = Dict("index" => i, "unit" => units[i], "theta" => θ, "status" => r.status,
                   "seconds" => r.seconds, "error" => r.error, "moments" => r.moments)
        lock(lk) do
            open(io -> println(io, JSON.json(rec)), path, "a")
            n_done[] += 1
            r.status == "ok" || (n_fail[] += 1)
            if n_done[] % 50 == 0 || n_done[] == length(todo)
                @printf("  [%s] %d/%d points, %d failed, %.1f min\n", Libc.strftime("%H:%M:%S", time()),
                        n_done[], length(todo), n_fail[], (time() - t0) / 60)
                flush(stdout)
            end
        end
    end
    println("Done: ", path)
end

main(OPTS)
