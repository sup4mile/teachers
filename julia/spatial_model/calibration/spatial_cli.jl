# =============================================================================
# Command-line plumbing shared by spatial_estimate.jl and spatial_diagnose.jl.
#
# Loaded on the main process before the model files. `parse_cli` needs nothing
# else; the reference and problem helpers call into spatial_calibrate.jl, which
# must be loaded by the time they run.
# =============================================================================

using Printf

const ALIASES = Dict("logkappa" => "logκ̃", "dkappa" => "δκ", "dB" => "ΔB", "beta" => "β",
                     "lambda" => "λ", "sigma" => "σ", "theta_nu" => "θν", "psi" => "ψ",
                     "eta" => "η", "phi" => "φ", "mu" => "μ", "gamma" => "γ", "alpha" => "α",
                     "sigma_eps" => "σϵ", "rho_z" => "ρz", "omega_f" => "ωf", "chi" => "χ",
                     "neT" => "nϵT")

"`name=value` → Symbol => value; the value is a Float64 when it parses as one, else a Symbol."
function parse_assignment(s)
    occursin('=', s) || error("expected name=value, got \"$s\"")
    k, v = split(s, '='; limit = 2)
    return Symbol(get(ALIASES, k, k)) => something(tryparse(Float64, v), Symbol(v))
end

"""
    parse_cli(args, defaults; flags, ints, floats, strings, lists, multi, script) -> Dict

Parse `--key value` and `--flag` arguments into a copy of `defaults`. Keys in
`lists` accumulate `name=value` pairs (`parse_assignment`) and keys in `multi` raw
strings; `--fix` values must be numbers. `script` names the file whose header
documents the options.
"""
function parse_cli(args, defaults::AbstractDict; flags = (), ints = (), floats = (), strings = (),
                   lists = (), multi = (), script = "the script")
    o = Dict{String,Any}(defaults)
    for k in lists
        o[k] = Pair{Symbol,Any}[]
    end
    for k in multi
        o[k] = String[]
    end
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
        key in multi   ? push!(o[key], String(v)) :
        error("unknown option --$key (see the header of $script)")
        i += 2
    end
    for (name, v) in get(o, "fix", ())
        v isa Float64 || error("--fix $name needs a number")
    end
    return o
end

"Worker count when --workers is not given: SLURM_CPUS_PER_TASK − 1, else 0 (serial)."
default_workers() = max(parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1")) - 1, 0)

# -----------------------------------------------------------------------------
# Reference levels and the estimation problem (need spatial_calibrate.jl)
# -----------------------------------------------------------------------------
"""
Compute the reference levels at the baseline external block and ϑ₀, and save them.
The fixed point starts from the levels already in `path`, else from the main frozen
reference (`REFERENCE_FILE`, e.g. for a --smoke run's own reference), else from the
`Params` defaults.
"""
function freeze_reference(path, solver)
    println("Freezing C̄ and h̄ at the baseline external block and ϑ₀ → $path")
    ext0  = external_params()
    start = isfile(path) ? path : isfile(REFERENCE_FILE) ? REFERENCE_FILE : nothing
    if start !== nothing
        old  = load_reference(start)
        ext0 = with_params(ext0; Cbar = old.Cbar, href = old.href)
    end
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
    names = get(ref.meta, "theta_names", nothing)
    names == collect(string.(THETA_NAMES)) ||
        println("  ⚠ the reference's ϑ is saved under names $names, not $(collect(string.(THETA_NAMES))); ",
                "refreeze so that its ϑ is in the current units")
    return ref
end

"""
    problem_setup(opts, ref) -> (; base, meas, targets, active, bounds, θ0)

The estimation problem the options describe: the fixed block with `--external`
overrides at the frozen reference `ref`, the measurement choices with `--measure`,
the targets (the BFM denominator from `--wtp-consumption`; values and scales from
`--target key=value` and `--target-se key=value`, where a NaN value drops the target;
`--activate key` and `--deactivate key` move a target into or out of the criterion),
and the search box and start with `--fix` parameters pinned. `--external` may not
name a parameter of ϑ.
"""
function problem_setup(opts, ref)
    for (name, _) in opts["external"]
        name in THETA_NAMES && error("--external $name: $name is in ϑ; hold it with --fix $name=… instead")
    end
    base    = calibration_base(ref; opts["external"]...)
    meas    = calibration_measurement(; opts["measure"]...)
    targets = default_targets()
    opts["wtp-consumption"] === nothing || (targets = with_wtp_target(targets, opts["wtp-consumption"]))
    for (field, list) in ((:value, get(opts, "target", ())), (:se, get(opts, "target-se", ())))
        for (key, v) in list
            i = findfirst(t -> t.key === key, targets)
            i === nothing && error("--target $key: no such target; targets are $([t.key for t in targets])")
            t = targets[i]
            targets[i] = with_target(t; (field => v,)..., note = t.note * "  [overridden: $field = $v]")
        end
    end
    for (flag, list) in ((true, get(opts, "activate", ())), (false, get(opts, "deactivate", ())))
        for key in list
            i = findfirst(t -> t.key === Symbol(key), targets)
            i === nothing && error("--(de)activate $key: no such target; targets are $([t.key for t in targets])")
            targets[i] = with_target(targets[i]; active = flag)
        end
    end
    active  = active_targets(targets)
    bounds  = copy(THETA_BOUNDS)
    θ0      = copy(THETA0)
    for (name, v) in opts["fix"]
        j = theta_index(name)
        bounds[j] = (v, v)
        θ0[j] = v
    end
    return (; base, meas, targets, active, bounds, θ0)
end
