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
                     "neT" => "nϵT", "zeta" => "ζ", "rho_w" => "ρw", "tau_k" => "τk", "sigma_w" => "σw",
                     # M3: the education block (Params fields `ρc` … `σd`, `ωz`; ϑ entries `ρ_c`, `ρ_T0`, `ϑ_s`, `σ_ϵ`, `dlogφ`)
                     "rho_c" => "ρ_c", "rho_T0" => "ρ_T0", "vartheta_s" => "ϑ_s", "theta_s" => "ϑ_s",
                     "dlogphi" => "dlogφ", "sigma_eps_theta" => "σ_ϵ",
                     "sigma_d" => "σd", "omega_z" => "ωz",
                     # the local peer effect (ϑ entry `ψ_a`, Params field `ψa`); `psi` is the warm-glow curvature ψ
                     "psi_a" => "ψ_a",
                     # the switches of 2026-10-15: commuting friction κ_c (Params `κc`), BA suburban amenity ζ_BA (Params `Bd`);
                     # `theta_nu` stays the external θν, and `THETA_SYNONYMS` reads it as ϑ's θ_ν under --free and --fix;
                     # closure R: teachers' suburban workplace taste b_T (Params `BT`), already ASCII ("b_T"; "bT" also reads)
                     "kappa_c" => "κ_c", "zeta_BA" => "ζ_BA", "zeta_ba" => "ζ_BA", "bT" => "b_T",
                     # round 2026-10-15c: the income-graded suburban taste ζ_y (Params `By`) and its frozen reference `lyref`
                     "zeta_y" => "ζ_y", "ly_ref" => "lyref")

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
    # a reference saved before some extensions existed (a leading part of THETA_NAMES) is
    # current if ϑ₀'s missing extensions are off
    n = names isa AbstractVector ? length(names) : 0
    (N_THETA_BASE ≤ n ≤ length(THETA_NAMES) && names == collect(string.(THETA_NAMES[1:n])) &&
     all(iszero, THETA0[n+1:end])) ||
        println("  ⚠ the reference's ϑ is saved under names $names, not $(collect(string.(THETA_NAMES))); ",
                "refreeze so that its ϑ is in the current units")
    return ref
end

"Default file of the audited district targets by partition scheme (data/spatial/district_settings.py)."
const DISTRICT_SETTINGS = joinpath(@__DIR__, "..", "..", "..", "data", "spatial", "estimates", "settings",
                                   "district_settings_targets.toml")

"The four district targets a `--district` scheme replaces (values, SEs and correlations)."
const DISTRICT_KEYS = (:gap_pupils, :gap_salary, :gap_teachers_pp, :gap_income)

"""
The `Params` fields that `theta_to_params` sets from entries of ϑ (M3, the peer effect ψa, the switches κc and Bd
of 2026-10-15, and BT of closure R): `--external` may not set them. (σν is set from ϑ only while θ_ν ≠ 0 and is not
an external key: the external block's θν is, so `--external theta_nu` stays valid.)
"""
const THETA_SET_FIELDS = (:ρc, :cC, :cCz, :ρT0, :ϑs, :cT0, :ψa, :κc, :Bd, :BT, :By)

"The entries of ϑ that act only with the joint education and occupation choice (`--external college=1`)."
const COLLEGE_PARAMETERS = (:ρ_c, :c_C, :c_Cz, :ρ_T0, :ϑ_s, :c_T0, :ζ_BA)

"The entries of ϑ that act only under the residence/workplace split (`--external split=1`)."
const SPLIT_PARAMETERS = (:κ_c, :b_T)

"""
    external_start!(θ, bounds, base) -> θ

M3: the entries s_z and σ_ϵ of ϑ mean the external block's value while they are 0, which lies outside their
box once `--free` releases them (`EXTENSION_BOUNDS`); so does θ_ν (2026-10-15), the location-taste scale σν/μ.
Such an entry starts at the external block's value of `base`, which is the same model (`theta_to_params`),
instead of at the nearest end of the box. Entries that are fixed, or whose value is inside the box, are left alone.
"""
function external_start!(θ, bounds, base)
    for (name, v) in ((:s_z, base.σξ / sqrt(1 - base.ρz^2)), (:σ_ϵ, base.σϵ), (:θ_ν, base.σν / base.μ))
        j = theta_index(name)
        lo, hi = bounds[j]
        (lo < hi && !(lo ≤ θ[j] ≤ hi)) && (θ[j] = clamp(v, lo, hi))
    end
    return θ
end

"""
    problem_setup(opts, ref) -> (; base, meas, targets, active, bounds, θ0, W, floor_tol)

The estimation problem the options describe: the fixed block with `--external`
overrides at the frozen reference `ref`, the measurement choices with `--measure`,
the targets (the BFM denominator from `--wtp-consumption`; the four district targets
of partition scheme `--district SCHEME` from `DISTRICT_SETTINGS`; values and scales from
`--target key=value` and `--target-se key=value`, where a NaN value drops the target;
`--activate key` and `--deactivate key` move a target into or out of the criterion),
the search box and start with `--free name` extension parameters released to
`EXTENSION_BOUNDS` (starting at the bound if its off value 0 lies outside the box, as for ρw)
and `--fix` parameters pinned, and the whitening matrix `W` of the
criterion (`--weight corr`: the bootstrap correlations of the district targets of the
`--district` scheme, default locale; otherwise the identity, i.e. diagonal weights).
`--external` may not name a parameter of ϑ, except s_z (an external input whose value holds
while ϑ's s_z is 0), nor σw (= ρw σν, with ρw in ϑ: `external_params` errors), nor a `Params` field that ϑ sets
(`ρc`, `cC`, `cCz`, `ρT0`, `ϑs`, `cT0`, `ψa`, `κc`, `Bd`, `BT`: `theta_to_params` would overwrite it, `BT` only
while b_T ≠ 0, so an external `BT` would hold at b_T = 0 alone); the M3 education block is external
(`--external college=1`, `sBA`, `sT0`, `Fc`, `sigma_d`, `omega_z`) except its ϑ entries (`rho_c`, `c_C`, `c_Cz`,
`rho_T0`, `vartheta_s`, `s_z`, `sigma_eps`, `dlogphi`, `c_T0`), which `--fix`/`--free` set. `--activate
pooled_p9010` turns on the pooled 90/10 in `model_moments`. A target that exists only under the residence/workplace
split (`SPLIT_ONLY_MOMENTS`) may not be active without `--external split=1`. An active target that is an education
moment (`COLLEGE_MOMENTS`) or a freed education parameter (ρ_c, c_C, c_Cz, ρ_T0, ϑ_s, c_T0, ζ_BA) needs
`--external college=1`: without it the moment is NaN and the parameter has no effect, which the search would
report only as 'nonfinite moments' or a singular Jacobian, so it is an error here. The same holds for the
switches of 2026-10-15 when freed or fixed at a nonzero value: κ_c and b_T (`SPLIT_PARAMETERS`) need `--external
split=1` and ζ_BA `--external college=1` (a fixed zero is the off value and is accepted). `--free theta_nu` and
`--fix theta_nu=…` address ϑ's θ_ν (`THETA_SYNONYMS`), the location-taste scale σν = θ_ν μ that overrides the
external `--external theta_nu` while nonzero; a freed θ_ν starts at the external block's value. Two diagnostic switches of
round 2026-10-16 are off by default: `--external lyT_fixed=x` (a `Params` field: every teacher node's income-graded taste
is graded at ln ŷ = x, inert while By = 0) and `--measure cfr_rho=ρ` (a `Measurement` field, 0 < ρ ≤ 1: the CFR
bridge's SD of teacher effectiveness is sd(log h_eff)/ρ); each is an error under the other option. `floor_tol` is the largest
population share
allowed at the consumption floor, `--floor-tol` (default `FLOOR_MASS_TOL`), which the caller
hands to `SpatialModel`. With `--district`, the payroll-per-pupil validation target is the
sum of the scheme's salary and FTE-per-pupil gaps, with their bootstrap correlation.
"""
function problem_setup(opts, ref)
    for (name, _) in opts["external"]
        # s_z is also an external input: its value holds while ϑ's s_z entry is 0 (M3)
        name in THETA_NAMES && name !== :s_z && error("--external $name: $name is in ϑ; hold it with --fix $name=… instead")
        # the Params fields that theta_to_params sets from ϑ (M3): an external value would be silently overwritten
        name in THETA_SET_FIELDS && error("--external $name: theta_to_params sets $name from ϑ; hold it with --fix instead")
        # round 2026-10-16: the CFR effectiveness concept is a measurement choice, not a Params field
        name === :cfr_rho && error("--external cfr_rho: the CFR effectiveness concept is a measurement choice; use --measure cfr_rho=…")
    end
    for (name, _) in opts["measure"]
        name === :lyT_fixed && error("--measure lyT_fixed: it is a model switch (Params); use --external lyT_fixed=…")
    end
    base    = calibration_base(ref; opts["external"]...)
    meas    = calibration_measurement(; opts["measure"]...)
    # round 2026-10-16: catch a bad diagnostic switch here, not in the first solve (a search would record every point as failed)
    (0 < meas.cfr_rho ≤ 1) ||
        error("--measure cfr_rho=$(meas.cfr_rho): the correlation of teachers' log effectiveness with log h_eff must lie in (0, 1]")
    isinf(base.lyT_fixed) && error("--external lyT_fixed=$(base.lyT_fixed): must be finite (NaN: off)")
    targets = default_targets()
    opts["wtp-consumption"] === nothing || (targets = with_wtp_target(targets, opts["wtp-consumption"]))
    scheme  = get(opts, "district", nothing)
    if scheme !== nothing
        tab = TOML.parsefile(DISTRICT_SETTINGS)[scheme]
        for key in DISTRICT_KEYS
            i = findfirst(t -> t.key === key, targets)
            v = tab[string(key)]
            targets[i] = with_target(targets[i]; value = Float64(v["value"]), se = Float64(v["se"]),
                                     note = targets[i].note * "  [district scheme $scheme]")
        end
        i  = findfirst(t -> t.key === :payroll_pp_gap, targets)
        ck = tab["corr_keys"]
        ρ  = tab["corr"][findfirst(==("gap_salary"), ck)][findfirst(==("gap_teachers_pp"), ck)]
        sal, fte = tab["gap_salary"], tab["gap_teachers_pp"]
        targets[i] = with_target(targets[i]; value = Float64(sal["value"] + fte["value"]),
                                 se = sum_se(sal["se"], fte["se"], ρ),
                                 note = targets[i].note * "  [district scheme $scheme]")
    end
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
    split_only = [t.key for t in active if t.key in SPLIT_ONLY_MOMENTS]
    isempty(split_only) || base.split ||
        error("target(s) $(join(split_only, ", ")) exist only under the residence/workplace split: add --external split=1")
    # M3: the education moments are NaN without `college`; a freed education parameter has no effect
    college_active = [t.key for t in active if t.key in COLLEGE_MOMENTS]
    (isempty(college_active) || base.college) ||
        error("--activate: $(join(college_active, ", ")) are education moments, which the model has only with the joint education and occupation choice; add --external college=1")
    base.college || begin
        freed = [string(n) for n in get(opts, "free", ()) if get(THETA_SYNONYMS, Symbol(get(ALIASES, n, n)), Symbol(get(ALIASES, n, n))) in COLLEGE_PARAMETERS]
        isempty(freed) || error("--free $(join(freed, ", ")): education parameters have no effect without --external college=1")
    end
    # the switches κ_c, b_T (split) and ζ_BA (college) are inert without their extension: freed, or fixed at a nonzero value, an error
    for (names, on, flag) in ((SPLIT_PARAMETERS, base.split, "split=1"), ((:ζ_BA,), base.college, "college=1"))
        on && continue
        theta_key(n) = (k = Symbol(get(ALIASES, n, n)); get(THETA_SYNONYMS, k, k))
        bad = [string(theta_key(n)) for n in get(opts, "free", ()) if theta_key(n) in names]
        append!(bad, string(theta_key(n)) for (n, v) in opts["fix"] if theta_key(n) in names && !iszero(v))
        isempty(bad) || error("--free/--fix $(join(unique(bad), ", ")): the parameter has no effect without --external $flag")
    end
    # the teacher-demand closure (round 2026-10-15c, `--external fte_rule=…`): b_T is solved in GE, and the FTE gap is
    # imposed, so b_T may not be freed or fixed at a nonzero value and gap_teachers_pp may not be a target
    if !isnan(base.fte_rule)
        any(t.key === :gap_teachers_pp for t in active) &&
            error("--external fte_rule: the closure imposes the FTE-per-pupil gap; it cannot also be a target (drop --activate gap_teachers_pp)")
        bk(n) = (k = Symbol(get(ALIASES, n, n)); get(THETA_SYNONYMS, k, k))
        (any(bk(n) === :b_T for n in get(opts, "free", ())) || any(bk(n) === :b_T && !iszero(v) for (n, v) in opts["fix"])) &&
            error("--external fte_rule: b_T is the closure's shadow price, solved in equilibrium; do not free or fix it")
    end
    # the pooled 90/10 needs a second pseudo-sample: computed only if it is in the criterion
    any(t.key === :pooled_p9010 for t in active) && !meas.pooled_p9010 &&
        (meas = with_measurement(meas; pooled_p9010 = true))
    bounds  = copy(THETA_BOUNDS)
    θ0      = copy(THETA0)
    for name in get(opts, "free", ())
        s = Symbol(get(ALIASES, name, name))
        s = get(THETA_SYNONYMS, s, s)
        haskey(EXTENSION_BOUNDS, s) || error("--free $name: only the extension parameters $(collect(keys(EXTENSION_BOUNDS))) are fixed by default")
        bounds[theta_index(s)] = EXTENSION_BOUNDS[s]
        # an off value 0 outside the box starts at the bound (ρw); s_z, σ_ϵ and θ_ν start at the external block's values instead
        # (`external_start!` below), which the same model gives: clamping them here would send them to the box's ends
        s in (:s_z, :σ_ϵ, :θ_ν) || (θ0[theta_index(s)] = clamp(θ0[theta_index(s)], EXTENSION_BOUNDS[s]...))
    end
    for (name, v) in opts["fix"]
        j = theta_index(name)
        bounds[j] = (v, v)
        θ0[j] = v
    end
    external_start!(θ0, bounds, base)
    weight = get(opts, "weight", "diag")
    weight in ("diag", "corr") || error("--weight must be diag or corr, got $weight")
    corr = weight == "corr" ? load_target_corr(DISTRICT_SETTINGS, something(scheme, "locale")) : nothing
    W    = target_whitener(active, corr)
    floor_tol = something(get(opts, "floor-tol", nothing), FLOOR_MASS_TOL)
    floor_tol ≥ 0 || error("--floor-tol must be nonnegative, got $floor_tol")
    return (; base, meas, targets, active, bounds, θ0, W, floor_tol)
end
