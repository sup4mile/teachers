# =============================================================================
# Delta-method uncertainty of the model's validation moments at a fitted ϑ̂
# (julia/spatial_model/calibration/calibration_log.md §5, remaining work: calibration bootstrap).
#
# A linearized first step toward the full calibration bootstrap, cheap at an exactly
# identified fit. At ϑ̂, with free parameters F, active targets A and every scalar moment
# of `model_moments` (vectors flattened to key[i]):
#   1. central differences of all moments in each free parameter, J = ∂m/∂ϑ_F; the target
#      rows J_A give the refit map G = ∂ϑ̂/∂t_A = J_A⁻¹ (square J_A: exact identification),
#      else (J_A′ΩJ_A)⁻¹J_A′Ω with the criterion weights Ω;
#   2. the targets' covariance Σ_A = D R D: D the scales (literature targets' scales times
#      --transport), R the correlations of the whole-CZ bootstrap draws of the district targets
#      (--district scheme, default locale), identity elsewhere;
#   3. Var(ϑ̂) = G Σ_A G′, the parameters' SEs;
#   4. external inputs e (--externals): ∂m/∂e by central differences at fixed ϑ̂ (σϵ
#      re-inverts the occupational block, as `--external sigma_eps=` does), the implicit
#      refit dϑ/de = −G ∂m_A/∂e and the total effect dm/de = ∂m/∂e + J dϑ/de; the
#      externals are independent of each other and of the targets;
#   5. each validation moment v: the model value, its SE from the targets (J_v Var(ϑ̂) J_v′)
#      and from the externals (Σ_e (dm_v/de SE_e)²), the data value and SE, and
#          z = (model − data) / √(SE_model² + SE_data² − 2 cov),   cov = J_v G Σ_{A,v},
#      where Σ_{A,v} is the covariance between the targets and the data moment through
#      shared district sampling: the district moments of the scheme's bootstrap draws (the four
#      district targets, the poverty, SEDA growth and score gaps) share one resampling of CZs, so
#      their draws are correlated (`load_draw_corr`); moments from other samples are independent.
#      cov = Cov(m_v, d_v) is the covariance of the model prediction and the data moment: the
#      prediction moves with the targets' sampling error ε_A as J_v G ε_A and the data moment
#      by ε_v, so Var(m_v − d_v) = Var(m_v) + Var(d_v) − 2 J_v G Σ_{A,v};
#   6. report.md (tables), uncertainty.toml (every number) and jacobian.toml (J, ∂m/∂e).
# Solves: 1 + 2|F| + 2|E| (25 at the baseline with four externals), spread over the workers
# with pmap. The validation set is the inactive targets with data values, the district
# validation moments of the --district scheme (poverty, SEDA growth and score gaps), the NLSY first-stage moments
# (wage slope, mother–child ρ) and the payroll per pupil gap (salary + FTE per pupil).
#
# With E = julia/spatial_model/calibration/estimation and U = julia/spatial_model/calibration/spatial_uncertainty.jl:
#   julia --project=E U --theta julia/spatial_model/calibration/runs/exact-polish-base --workers 16 --out DIR
#   julia --project=E U --theta DIR0 --externals sigma_eps,sigma,theta_nu=0.05 --transport 2 --out DIR
#   julia --project=E U --theta DIR0 --reuse DIR --workers 0 --target-se move_rate=0.011 --out DIR2
#
# Options (defaults in brackets)
#   --theta SPEC           theta0; a TikTak run directory (its best point); a TOML file or a
#                          directory holding theta.toml (from --out); or v1,v2,… [theta0]
#   --workers N            local worker processes [SLURM_CPUS_PER_TASK − 1, else 0]
#   --out DIR              write report.md, uncertainty.toml and jacobian.toml to DIR
#   --externals LIST       external inputs to perturb, comma-separated; `name=se` sets or
#                          overrides the SE [sigma_eps,rho_z,s_z,sigma]. Default SEs: ρz 0.014
#                          and s_z 0.0071 (NLSY, first_stage.md), σϵ 0.02 (a specification SD:
#                          a quarter of the first-stage sensitivity range 0.766–0.844, not a
#                          sampling SE), σ the FOO (2013) SE at η and K (0.296 at K = 12).
#                          With the M3 entries s_z and σ_ϵ of ϑ nonzero (free or fixed) the externals
#                          of those names are overridden by ϑ: the default list drops them (with a
#                          note), an explicit list naming them is an error
#   --ext-step X           central-difference step of an external, in its SEs [0.1]. At the
#                          baseline 0.03 and 0.1 give moment SEs within about 2% (σ's contribution
#                          within 1.5%); 0.3 is too coarse (σ's contribution moves by up to 12%, the external SE of gap_income by 9% and of the FTE gap by 7%)
#   --transport K          literature SE inflation factor for cfr_effect and wtp_share [1];
#                          the report also gives K = 1 and K = 2
#   --reuse DIR            skip the solves: take J, ∂m/∂e and the model values from DIR/jacobian.toml
#                          and redo the delta method at other target scales, weights, active sets or
#                          transport factors; use with --workers 0. Refused unless ϑ, the free
#                          parameters and their steps, the externals' values and steps and the model
#                          specification (fixed block, measurement, solver, --floor-tol, model source
#                          files) are those of the solves in DIR
#   --district SCHEME      the four district targets of a partition scheme, and the bootstrap
#                          draws whose correlations Σ_A and cov use (locale, revenue, ses, …) [locale]
#   --weight diag|corr     criterion weights Ω of the refit map when overidentified (as
#                          spatial_diagnose.jl); irrelevant at an exact fit [diag]
#   --floor-tol X          largest population share allowed at the consumption floor, with or
#                          without a feasible location (as spatial_diagnose.jl; needs model code
#                          that has it) [FLOOR_MASS_TOL]; used in every solve, perturbed ones too
#   --grid name=value      override a SOLVER setting (repeatable)
#   --fix, --free, --external, --measure, --target, --target-se, --activate, --deactivate,
#   --wtp-consumption, --reference:
#                          as in spatial_estimate.jl and spatial_diagnose.jl
# =============================================================================

using Distributed
using Printf
using SHA: sha256
using Statistics: cor, std
using TOML
using TikTak: load_estimates

const HERE = @__DIR__
include(joinpath(HERE, "spatial_cli.jl"))

"The default of `--externals`."
const DEFAULT_EXTERNALS = "sigma_eps,rho_z,s_z,sigma"

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "workers" => default_workers(), "theta" => "theta0", "out" => nothing,
        "externals" => DEFAULT_EXTERNALS, "ext-step" => 0.1, "transport" => 1.0,
        "reuse" => nothing, "wtp-consumption" => nothing, "reference" => nothing,
        "district" => nothing, "weight" => "diag", "floor-tol" => nothing);
    ints    = ("workers",),
    floats  = ("wtp-consumption", "ext-step", "transport", "floor-tol"),
    strings = ("theta", "out", "reference", "district", "weight", "externals", "reuse"),
    lists   = ("fix", "external", "measure", "target", "target-se", "grid"),
    multi   = ("activate", "deactivate", "free"),
    script  = "spatial_uncertainty.jl")
OPTS["workers"] > 0 && OPTS["reuse"] === nothing &&
    addprocs(OPTS["workers"]; exeflags = "--project=$(Base.active_project())")

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($HERE, "spatial_calibrate.jl"))
end

@everywhere begin
    """
        flatten_moments(m) -> (names, values)

    The scalar and short-vector fields of `model_moments` as flat lists, a vector `key` as
    `key[1]`, `key[2]`, …; other fields are dropped.
    """
    function flatten_moments(m)
        names, vals = String[], Float64[]
        for (k, v) in pairs(m)
            if v isa Real
                push!(names, string(k)); push!(vals, Float64(v))
            elseif (v isa AbstractVector || v isa Tuple) && all(x -> x isa Real, v)
                for (i, x) in enumerate(v)
                    push!(names, string(k, "[", i, "]")); push!(vals, Float64(x))
                end
            end
        end
        return names, vals
    end

    """
        solve_point(model, θ) -> (names, values) | String

    Solve at ϑ and flatten every moment, or return the message of a `CalibrationFailure` or
    `DomainError` (the failures the search records as failed points). Defined on every worker.
    """
    function solve_point(model::SpatialModel, θ)
        try
            _, m = solve_theta(model, θ)
            return flatten_moments(m)
        catch e
            e isa Union{CalibrationFailure,DomainError} || rethrow()
            return sprint(showerror, e)
        end
    end
end

# -----------------------------------------------------------------------------
# The point ϑ
# -----------------------------------------------------------------------------
"""
ϑ from a --theta specification (see the header). A TOML file must list `theta_names`
equal to `THETA_NAMES` or to a leading part of it of at least eight names (files saved
before some extension parameters existed, read with those extensions off; older files
still hold log κ̄ as `logκ`, not log κ̃); TikTak journals carry no names, so only their
length is checked. Short vectors are padded with the missing extensions off (`pad_theta`).
Same as `read_theta` of spatial_diagnose.jl.
"""
function read_theta(spec::AbstractString)
    spec == "theta0" && return copy(THETA0)
    checklen(θ) = (N_THETA_BASE ≤ length(θ) ≤ length(THETA_NAMES) ||
                   error("--theta $spec: $(length(θ)) values, ϑ has $(length(THETA_NAMES))"); pad_theta(θ))
    if isdir(spec)
        isfile(joinpath(spec, "theta.toml")) && return read_theta(joinpath(spec, "theta.toml"))
        return checklen(load_estimates(spec; limit = 1)[1])
    end
    if isfile(spec)
        d = TOML.parsefile(spec)
        names = get(d, "theta_names", nothing)
        (names isa AbstractVector && N_THETA_BASE ≤ length(names) ≤ length(THETA_NAMES) &&
         names == collect(string.(THETA_NAMES[1:length(names)])) && length(d["theta"]) == length(names)) ||
            error("--theta $spec: theta_names $names differ from $(collect(string.(THETA_NAMES)))")
        return checklen(Float64.(d["theta"]))
    end
    return checklen(parse.(Float64, split(spec, ',')))
end

# -----------------------------------------------------------------------------
# External inputs
# -----------------------------------------------------------------------------
"""
Default SEs of the external inputs `--externals` perturbs (σ is computed from FOO's SE at the
base's η and K, `sigma_from_wages`: 0.296 at η = 0.080, K = 12; 95% interval [−0.015, 1.144]).
The ρz and s_z SEs are the NLSY79 sampling SEs of data/spatial/estimates/first_stage.md
(s_z = slope × (1 − η)/α, slope 0.2127 with SE 0.0077). σϵ has no sampling SE (the ACS 90/10
is precise): 0.02 is a specification SD, a quarter of the first-stage sensitivity range
0.766–0.844 (age-residualized 90/10 to b s_z = χ).
"""
const EXTERNAL_SE = Dict{Symbol,NamedTuple{(:se, :basis),Tuple{Float64,String}}}(
    :σϵ  => (se = 0.02,   basis = "specification SD (0.844 − 0.766)/4, a quarter of the first-stage range; not a sampling SE"),
    :ρz  => (se = 0.014,  basis = "NLSY79/CNLSY latent mother–child correlation 0.577 (SE 0.014), first_stage.md"),
    :s_z => (se = 0.0071, basis = "NLSY79 latent wage slope 0.2127 (SE 0.0077) × (1 − η), first_stage.md"))

"""
    overridden_externals(θ) -> Vector{Symbol}

The external inputs that `theta_to_params` replaces at ϑ (M3): σϵ when ϑ's σ_ϵ is nonzero and s_z when ϑ's s_z is
nonzero, and (2026-10-15) θν when ϑ's θ_ν is nonzero (freed entries start at the external values, `external_start!`,
so they are nonzero). Perturbing such an external changes nothing (identically zero columns of ∂m/∂e), and its
first-stage link to a wage slope no longer holds.
"""
overridden_externals(θ) = Symbol[n for (n, j) in ((:σϵ, theta_index(:σ_ϵ)), (:s_z, theta_index(:s_z)), (:θν, theta_index(:θ_ν)))
                                 if !iszero(θ[j])]

"The ϑ entry that overrides an external input of `overridden_externals`."
const OVERRIDING_ENTRY = Dict(:σϵ => "σ_ϵ", :s_z => "s_z", :θν => "θ_ν")

"""
    external_inputs(opts, base, meas; θ = nothing) -> Vector{NamedTuple}

The externals of `--externals` as `(; name, value, se, step, basis)`: the value in `base` (a
`--external` override, else `EXTERNAL`, else the `Params` field), the SE (`EXTERNAL_SE`, FOO's
for σ, or `name=se` on the command line) and the central-difference step `--ext-step` × SE.
With ϑ (`θ`) the externals it overrides (`overridden_externals`) are dropped from the default list, with a
note, and are an error if `--externals` names them.
"""
function external_inputs(opts, base::Params, meas::Measurement; θ = nothing)
    over = (; opts["external"]...)
    out  = NamedTuple[]
    overridden = θ === nothing ? Symbol[] : overridden_externals(θ)
    for item in split(opts["externals"], ','; keepempty = false)
        nm, given = occursin('=', item) ? split(item, '='; limit = 2) : (item, nothing)
        nm   = String(strip(nm))
        name = Symbol(get(ALIASES, nm, nm))
        if name in overridden
            opts["externals"] == DEFAULT_EXTERNALS ||
                error("--externals $nm: ϑ's $(OVERRIDING_ENTRY[name]) is nonzero (free or fixed), so theta_to_params overrides the external $name; " *
                      "perturbing it changes nothing. Leave it out of --externals (its uncertainty is that of the estimate)")
            println("  note: external $name left out of the perturbations: ϑ's ", OVERRIDING_ENTRY[name],
                    " is nonzero, so theta_to_params overrides it (its columns would be identically zero)")
            continue
        end
        name in THETA_NAMES && name !== :s_z && error("--externals $nm: $name is in ϑ, not an external input")
        value = haskey(over, name) ? over[name] : hasproperty(EXTERNAL, name) ? EXTERNAL[name] : getfield(base, name)
        value isa Real || error("--externals $nm: the value $value is not a number")
        # round 2026-10-15c: a switch that is off (fte_rule = NaN) or inert (lyref with By = 0) has nothing to perturb
        isfinite(value) || error("--externals $nm: its value is $value (the switch is off), so perturbing it changes nothing")
        name === :lyref && all(iszero, base.By) &&
            error("--externals lyref: the income-graded taste is off (By = 0), so lyref is inert")
        if given !== nothing
            se, basis = parse(Float64, given), "--externals $item"
        elseif name === :σ
            se, basis = sigma_from_wages(; η = base.η, K = meas.K).se, "FOO (2013) Table V, SE 0.0033 per pupil mapped at η, K"
        elseif haskey(EXTERNAL_SE, name)
            se, basis = EXTERNAL_SE[name]
        else
            error("--externals $nm: no default SE; give it as $nm=se")
        end
        se > 0 || error("--externals $nm: the SE must be positive")
        push!(out, (; name, value = Float64(value), se, step = opts["ext-step"] * se, basis))
    end
    allunique(e.name for e in out) || error("--externals lists a parameter twice")
    return out
end

"The `--external` overrides with `name` set to `value`, as keywords for `calibration_base`."
external_overrides(list, name::Symbol, value) =
    Pair{Symbol,Any}[(p for p in list if first(p) !== name)..., name => value]

# -----------------------------------------------------------------------------
# Data moments and their covariance
# -----------------------------------------------------------------------------
"The literature targets: their scales are the papers' SEs and are inflated by --transport."
const LITERATURE_KEYS = (:cfr_effect, :wtp_share)

"Moments whose `se` is a modelling scale (a tolerance), not a sampling SE."
const SCALE_KEYS = (:male_teach_share, :female_teach_share, :teacher_p9010, :move_rate,
                    :teacher_premium_m, :teacher_premium_f, :s_ratio, :gender_gap_market,
                    :teacher_supply_elasticity,
                    # M3 (spatial_calibrate.jl `default_targets`): 5% of a BA share or of a 90/10, 0.01 for the BA share of
                    # teachers and the teachers' schooling, 0.02 for the BA premium, 0.005 for s_O, 0.05 for the teacher premia and the
                    # BA-share gaps; the others (ba_z_gap_sd, wage_slope_ba, teacher_ba_wage_gap, teacher z gaps, ma_premium_implied) are SEs
                    :ba_share_nonteach, :ba_share, :ba_share_teach, :ba_premium, :s_T_mean_ba, :s_T_mean_nonba, :s_O_mean,
                    :pooled_p9010, :teacher_p9010_ba, :teacher_premium_ba_m, :teacher_premium_ba_f, :teacher_premium_within_ba_m,
                    :teacher_premium_within_ba_f, :ba_share_z_gap, :ba_share_score_gap)

"""
NLSY first-stage moments the model's `wage_slope` and `rho_pc` are checked against: the latent
log-wage slope per latent AFQT SD and the latent mother–child score correlation, with SEs
(`NLSY` in spatial_first_stage.jl; first_stage.md). They are also the inputs of s_z and ρz.
"""
const NLSY_MOMENTS = (; wage_slope = (0.21266045167782227, 0.007708547769715874),
                        rho_pc     = (0.5766299079601815, 0.014023593002952637))

"""
Validation moments in the order of the pre-registered validation table; inactive targets with
data values that are not listed follow. A moment that `model_moments` does not return or leaves
undefined (NaN: the split's commuting moments and supply elasticity in a bundled run), has no
data value, or is an active target is skipped.
"""
const VALIDATION_KEYS = (:gap_income, :gap_teachers_pp, :payroll_pp_gap, :poverty_gap, :score_gap,
                         :seda_growth_gap, :city_share_gradient, :city_share_gradient_z,
                         :teacher_premium_m, :teacher_premium_f, :s_ratio, :gender_gap_market,
                         :teacher_z_gap_market_sd, :wage_slope, :rho_pc)

"The external input each NLSY first-stage validation moment is the data source of (s_z from the wage slope, ρz from the mother–child correlation)."
const FIRST_STAGE_LINK = Dict(:wage_slope => :s_z, :rho_pc => :ρz)

"Validation moments that are linear combinations of moments: log payroll per pupil = log salary + log teachers per pupil."
const DERIVED_VALIDATION = Dict{Symbol,Vector{Pair{Symbol,Float64}}}(
    :payroll_pp_gap => [:gap_salary => 1.0, :gap_teachers_pp => 1.0])

"""
Model levels shown after the validation moments (those `model_moments` returns), with labels:
population shares, tax rates, class sizes, home production, schooling, WTP. The levels of
consumption and human capital are left out: they are in the frozen units of C̄ and h̄ and move with σ
by factors, not by SEs.
"""
const LEVELS = ("pop_share[1]" => "population share, location 1", "pop_share[2]" => "population share, location 2",
                "t[1]" => "tax rate, location 1", "t[2]" => "tax rate, location 2",
                "t_market" => "teacher payroll over market income", "teach_share" => "teaching share (all adults)",
                "students_per_teacher[1]" => "students per teacher, location 1",
                "students_per_teacher[2]" => "students per teacher, location 2",
                "hp_share_m" => "home production, men", "hp_share_f" => "home production, women",
                "s_T_mean" => "mean time share s, teachers", "s_O_mean" => "mean time share s, non-teachers",
                "p12" => "P(move 1 → 2)", "p21" => "P(move 2 → 1)",
                "wtp_income_elasticity" => "WTP income elasticity",
                "teacher_supply_elasticity" => "teacher supply elasticity (split only)")

"District validation moments of a scheme's table in `DISTRICT_SETTINGS` (data values, SEs and bootstrap draws) that are not targets."
const DISTRICT_VALIDATION = (:poverty_gap, :seda_growth_gap, :score_gap)

"""
    data_moments(targets, scheme) -> Dict{Symbol,NamedTuple}

Data value and SE of every moment the tool compares with: the targets that have both (active
or not) and the validation moments that are not targets: poverty, SEDA growth and score gaps
from the `scheme` table of `DISTRICT_SETTINGS`, and the NLSY wage slope and mother–child ρ.
"""
function data_moments(targets, scheme)
    dat = Dict{Symbol,NamedTuple{(:value, :se, :scale_only),Tuple{Float64,Float64,Bool}}}()
    for t in targets
        (isfinite(t.value) && isfinite(t.se) && t.se > 0) && (dat[t.key] = (; value = t.value, se = t.se, scale_only = t.key in SCALE_KEYS))
    end
    tab = TOML.parsefile(DISTRICT_SETTINGS)[scheme]
    for key in DISTRICT_VALIDATION
        haskey(dat, key) || !haskey(tab, string(key)) ||
            (dat[key] = (; value = Float64(tab[string(key)]["value"]), se = Float64(tab[string(key)]["se"]), scale_only = false))
    end
    for (key, (v, se)) in pairs(NLSY_MOMENTS)
        haskey(dat, key) || (dat[key] = (; value = v, se, scale_only = false))
    end
    return dat
end

"""
    load_draw_corr(scheme) -> (; keys, matrix, source)

Correlations of the whole-CZ bootstrap draws of the district moments of partition scheme
`scheme`: the four district targets (`corr_keys`) and the validation moments `DISTRICT_VALIDATION`,
read from `district_settings_draws_<scheme>.csv` next to `DISTRICT_SETTINGS`
(data/spatial/district_settings.py). One multiplicity matrix of CZ resamples is shared by every
moment, so the correlations are coherent across the targets and the validation moments;
`load_target_corr` returns only the 4 × 4 block of the targets, which leaves the covariance
of the model prediction with the poverty, SEDA growth and score gaps out. Correlations are pairwise
over finite draws. Errors if the draws do not reproduce the table's `corr` block and each SE
(the files are from different runs); without a draws file it falls back to `load_target_corr`
(the validation moments are then uncorrelated with the targets) and says so in `source`.
"""
function load_draw_corr(scheme)
    path = joinpath(dirname(DISTRICT_SETTINGS), "district_settings_draws_$scheme.csv")
    tab  = TOML.parsefile(DISTRICT_SETTINGS)[scheme]
    tgt  = load_target_corr(DISTRICT_SETTINGS, scheme)
    if !isfile(path)
        @warn "$path not found: the validation moments are taken as uncorrelated with the targets"
        return (; tgt..., source = "$(basename(DISTRICT_SETTINGS)), the four district targets only: " *
                                   "no draws file $(basename(path)); the validation moments are taken as uncorrelated with the targets")
    end
    lines = filter(!isempty ∘ strip, readlines(path))
    head  = Symbol.(strip.(split(lines[1], ',')))
    X     = [something(tryparse(Float64, strip(f)), NaN) for l in lines[2:end] for f in split(l, ',')]
    X     = permutedims(reshape(X, length(head), :))                    # draws × columns
    keys  = Symbol[tgt.keys..., (k for k in DISTRICT_VALIDATION if k in head && haskey(tab, string(k)) && !(k in tgt.keys))...]
    all(k -> k in head, keys) || error("$path has no column for $(join(setdiff(keys, head), ", "))")
    col   = Dict(k => X[:, findfirst(==(k), head)] for k in keys)
    finite(x, y) = [i for i in eachindex(x) if isfinite(x[i]) && isfinite(y[i])]
    matrix = [a === b ? 1.0 : (ok = finite(col[a], col[b]); cor(col[a][ok], col[b][ok])) for a in keys, b in keys]
    n = length(tgt.keys)
    maximum(abs, matrix[1:n, 1:n] .- tgt.matrix) < 1e-6 ||
        error("$path: the correlations of the draws differ from `corr` of $(basename(DISTRICT_SETTINGS)) [$scheme]; rerun district_settings.py")
    for k in keys
        haskey(tab, string(k)) || continue
        x = filter(isfinite, col[k])
        isapprox(std(x), tab[string(k)]["se"]; rtol = 1e-6) ||
            error("$path: the SD of the draws of $k ($(std(x))) differs from its SE in $(basename(DISTRICT_SETTINGS)) ($(tab[string(k)]["se"]))")
    end
    return (; keys, matrix, source = "$(basename(path)), $(size(X, 1)) draws")
end

"""
    data_cov(dat, corr, a, b, k) -> Float64

Covariance of the data estimates of moments `a` and `b`: the variance of one, or the bootstrap
correlation times the two SEs when both are in `corr` (`load_draw_corr`), else 0. The SEs of
the literature moments are multiplied by the transport factor `k`.
"""
function data_cov(dat, corr, a::Symbol, b::Symbol, k)
    sd(x) = dat[x].se * (x in LITERATURE_KEYS ? k : 1.0)
    a === b && return sd(a)^2
    i, j = findfirst(==(a), corr.keys), findfirst(==(b), corr.keys)
    (i === nothing || j === nothing) && return 0.0
    return corr.matrix[i, j] * sd(a) * sd(b)
end

"""
    validation_rows(targets, dat, names, m0) -> Vector{NamedTuple}

The validation moments to report: `(; label, comps)` with `comps` the moments (and
coefficients) whose linear combination it is. Skips what `model_moments` does not return, what
is not defined at the point (`m0`, the moments at ϑ in the order of `names`, is NaN there: the
commuting moments and the supply elasticity exist only under the split, and are zeros by
construction and not predictions in a bundled run), what has no data value, and active targets.
"""
function validation_rows(targets, dat, names, m0)
    active = [t.key for t in active_targets(targets)]
    labels = Symbol[VALIDATION_KEYS...]
    for t in targets
        (!t.active && haskey(dat, t.key) && !(t.key in labels)) && push!(labels, t.key)
    end
    at = Dict(n => i for (i, n) in enumerate(names))
    defined(key) = haskey(at, string(key)) && isfinite(m0[at[string(key)]])
    rows = NamedTuple[]
    for label in labels
        label in active && continue
        comps = get(DERIVED_VALIDATION, label, Pair{Symbol,Float64}[label => 1.0])
        all(c -> haskey(dat, first(c)) && defined(first(c)), comps) && push!(rows, (; label, comps))
    end
    return rows
end

# -----------------------------------------------------------------------------
# Finite differences
# -----------------------------------------------------------------------------
"""
    spatial_model(base, keys; meas, solver, floor_tol = nothing) -> SpatialModel

`SpatialModel` with the consumption-floor tolerance where the model code has one (`floor_tol` is
`nothing`: its default). Model code without `--floor-tol` (a snapshot from before 2026-10-05)
is an error if a tolerance is given, not silently ignored.
"""
function spatial_model(base, keys; meas, solver, floor_tol = nothing)
    floor_tol === nothing && return SpatialModel(base, keys; meas, solver)
    hasfield(SpatialModel, :floor_tol) || error("--floor-tol: this model code has no consumption-floor tolerance")
    return SpatialModel(base, keys; meas, solver, floor_tol)
end

"`model` with another fixed block `base` (an external moved): same keys, measurement, solver and floor tolerance."
with_base(model::SpatialModel, base) =
    spatial_model(base, model.keys; meas = model.meas, solver = model.solver,
                  floor_tol = hasfield(SpatialModel, :floor_tol) ? model.floor_tol : nothing)

"Central-difference steps of the parameters of ϑ: 0.001 of each box width (as `moment_jacobian`)."
fd_steps(bounds) = 1e-3 .* [hi - lo for (lo, hi) in bounds]

"""
    spec_fingerprint(model, bounds) -> String

Short hash of what the finite differences depend on besides ϑ and the externals' values: the
fixed block (frozen C̄ and h̄, externals, `--external` overrides), measurement choices, solver
settings, consumption-floor tolerance, the search box (the steps) and the source of the three model
files (as `problem_fingerprint`, so that edited model code cannot reuse old derivatives).
`--reuse` compares it with the one saved in jacobian.toml. The active set (`--activate`, `--deactivate`) is deliberately not part of it: the
Jacobian holds every moment, and the targets enter only the algebra (`--reuse` redoes it at other active sets).
"""
function spec_fingerprint(model::SpatialModel, bounds)
    io = IOBuffer()
    print(io, model.base, "|", model.meas, "|", model.solver, "|",
          hasfield(SpatialModel, :floor_tol) ? model.floor_tol : "no floor_tol", "|", bounds)
    for f in (joinpath("..", "spatial_continuous.jl"), "spatial_moments.jl", "spatial_calibrate.jl")
        write(io, read(joinpath(HERE, f)))
    end
    return bytes2hex(sha256(take!(io)))[1:12]
end

"""
    finite_differences(opts, ref, model, θ, free, bounds, ext) -> (; names, m0, J, dmde, step, solves, minutes)

Solve at ϑ, at ϑ ± one step in each free parameter (0.001 box widths, one-sided inside the box
at a bound, as `moment_jacobian`) and at ϑ with each external of `ext` moved by ± its step,
spread over the workers with `pmap`. Returns the flattened moment names, the moments at ϑ,
`J[v, j] = ∂m_v/∂ϑ_free[j]` and `dmde[v, e] = ∂m_v/∂e` at fixed ϑ (the external block rebuilt
by `calibration_base`, so σϵ re-inverts the occupational shares). A failed solve is an error.
"""
function finite_differences(opts, ref, model::SpatialModel, θ, free, bounds, ext)
    step = fd_steps(bounds)
    offs(j) = θ[j] + step[j] > bounds[j][2] ? (0, -2) :
              θ[j] - step[j] < bounds[j][1] ? (2, 0) : (1, -1)
    shifted(j, k) = (x = collect(Float64, θ); x[j] += k * step[j]; x)
    pts, labels = Tuple{SpatialModel,Vector{Float64}}[(model, collect(Float64, θ))], ["ϑ"]
    for j in free
        push!(pts, (model, shifted(j, offs(j)[1]))); push!(labels, "$(THETA_NAMES[j]) upper")
    end
    for j in free
        push!(pts, (model, shifted(j, offs(j)[2]))); push!(labels, "$(THETA_NAMES[j]) lower")
    end
    for e in ext, s in (1, -1)
        base_e = calibration_base(ref; external_overrides(opts["external"], e.name, e.value + s * e.step)...)
        push!(pts, (with_base(model, base_e), collect(Float64, θ)))
        push!(labels, "$(e.name) $(s > 0 ? "upper" : "lower")")
    end
    minutes = @elapsed res = pmap(p -> solve_point(p[1], p[2]), pts)
    for (i, r) in enumerate(res)
        r isa String && error("the solve at \"$(labels[i])\" failed: $r")
    end
    names = res[1][1]
    all(r[1] == names for r in res) || error("the solves return different moment lists")
    vals = [r[2] for r in res]
    nf   = length(free)
    J    = reduce(hcat, [(vals[1 + i] .- vals[1 + nf + i]) ./ ((offs(j)[1] - offs(j)[2]) * step[j])
                         for (i, j) in enumerate(free)])
    dmde = isempty(ext) ? zeros(length(names), 0) :
           reduce(hcat, [(vals[1 + 2nf + 2k - 1] .- vals[1 + 2nf + 2k]) ./ (2 * e.step) for (k, e) in enumerate(ext)])
    return (; names, m0 = vals[1], J, dmde, step = step[free], solves = length(pts), minutes = minutes / 60)
end

"""
Write the finite differences (`finite_differences`) to a TOML file, read back by `--reuse`: ϑ,
the free parameters and their steps, the moment names and values at ϑ, `J` (one row per moment),
the externals with their values and steps and `dmde` (one row per moment), and the fingerprint of
the model specification (`spec_fingerprint`).
"""
function save_jacobian(path, θ, free, fd, ext, fingerprint)
    d = Dict{String,Any}("theta" => collect(θ), "theta_names" => collect(string.(THETA_NAMES)),
                         "free" => collect(string.(THETA_NAMES[free])), "step" => fd.step,
                         "moments" => fd.names, "model" => fd.m0, "J" => [fd.J[v, :] for v in axes(fd.J, 1)],
                         "external_names" => [string(e.name) for e in ext], "external_value" => [e.value for e in ext],
                         "external_step" => [e.step for e in ext],
                         "dmde" => [fd.dmde[v, :] for v in axes(fd.dmde, 1)],
                         "fingerprint" => fingerprint,
                         "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    open(io -> TOML.print(io, d; sorted = true), path, "w")
end

"""
Read a `save_jacobian` file as `finite_differences` returns it. Errors unless it was taken at
the same ϑ, free parameters and steps, externals (names, values and steps) and model
specification (`fingerprint`, from `spec_fingerprint`): the derivatives are reused as they are.
A legacy file without a fingerprint is accepted with a warning (the externals are still checked); one that holds a
leading part of ϑ (before the extensions) is accepted if the later entries of ϑ are zero.
"""
function load_jacobian(path, θ, free, step, ext, fingerprint)
    d = TOML.parsefile(path)
    same(a, b) = length(a) == length(b) && all(isapprox.(Float64.(a), Float64.(b); rtol = 1e-10, atol = 1e-12))
    # a file from before ϑ was extended holds a leading part of ϑ; the later (extension) entries must be zero (off) now
    nθ = length(d["theta"])
    (nθ ≤ length(θ) && same(d["theta"], θ[1:nθ]) && all(iszero, θ[nθ+1:end])) ||
        error("--reuse: $path was taken at another ϑ")
    get(d, "theta_names", collect(string.(THETA_NAMES[1:nθ]))) == collect(string.(THETA_NAMES[1:nθ])) ||
        error("--reuse: $path has other parameter names than ϑ")
    d["free"] == collect(string.(THETA_NAMES[free])) || error("--reuse: $path has free parameters $(d["free"])")
    same(d["step"], step[free]) || error("--reuse: $path has finite-difference steps $(d["step"]), not $(step[free])")
    d["external_names"] == [string(e.name) for e in ext] || error("--reuse: $path has externals $(d["external_names"]); use --externals")
    same(d["external_value"], [e.value for e in ext]) ||
        error("--reuse: $path has the externals at $(d["external_value"]), not $([e.value for e in ext]) (--external)")
    same(d["external_step"], [e.step for e in ext]) ||
        error("--reuse: $path has external steps $(d["external_step"]), not $([e.step for e in ext]) (--ext-step, --externals)")
    if !haskey(d, "fingerprint")
        @warn "--reuse: $path predates the specification fingerprint: ϑ, free parameters, steps and the externals' names, values and steps " *
              "agree, but --measure, --grid, --fix, --floor-tol, --reference and the model code are NOT checked"
    elseif d["fingerprint"] != fingerprint
        error("--reuse: $path was taken with another model specification (fingerprint $(d["fingerprint"]), now $fingerprint): ",
              "fixed block (--external, --fix, --activate), --measure, --grid, --floor-tol, --reference, the search box or the model code differ; recompute without --reuse")
    end
    mat(rows, n) = isempty(rows) ? zeros(0, n) : reduce(vcat, permutedims.(Vector{Float64}.(rows)))
    nf, ne = length(free), length(ext)
    J = mat(d["J"], nf)
    dmde = ne == 0 ? zeros(size(J, 1), 0) : mat(d["dmde"], ne)
    return (; names = Vector{String}(d["moments"]), m0 = Float64.(d["model"]), J, dmde, step = Float64.(d["step"]),
            solves = 0, minutes = 0.0)
end

# -----------------------------------------------------------------------------
# The delta method
# -----------------------------------------------------------------------------
"""
    refit_matrix(JA, Ω) -> G

G = ∂ϑ̂/∂t_A, the response of the fitted parameters to the targets: J_A⁻¹ when J_A is square (the
exactly identified fit), else (J_A′ΩJ_A)⁻¹J_A′Ω for the criterion weights Ω (the linearized
first-order condition J_A′Ω(m_A(ϑ̂) − t_A) = 0).
"""
function refit_matrix(JA, Ω)
    n, f = size(JA)
    n ≥ f || error("ϑ is not identified: $n active targets for $f free parameters")
    return n == f ? inv(JA) : (JA' * Ω * JA) \ (JA' * Ω)
end

"""
    analyze(; θ, free, bounds, act, W, names, m0, J, dmde, ext, corr, dat, ks, rows) -> NamedTuple

The delta method (header, steps 3–5) for each literature transport factor in `ks`. Returns the
refit maps, covariances and SEs of the parameters and of every moment, the validation rows (model
value; SE from the targets, the externals and combined; data value and SE; the covariance `cov`
through shared district sampling; z; the contribution |J_v G[:, a]| × scale of each active target `a`
and |∂m_v/∂e| × SE of each external `e`) and the identification diagnostics.
"""
function analyze(; θ, free, bounds, act, W, names, m0, J, dmde, ext, corr, dat, ks, rows)
    idx   = Dict(n => i for (i, n) in enumerate(names))
    rowsA = [idx[string(t.key)] for t in act]
    JA    = J[rowsA, :]
    all(isfinite, JA) || error("the target rows of the Jacobian are not finite")
    se0   = [t.se for t in act]
    width = [bounds[j][2] - bounds[j][1] for j in free]
    sv    = svdvals((JA ./ se0) .* width')
    sv[end] > 1e-10 * sv[1] ||
        error("the target rows of the Jacobian are singular: standardized singular values $(join(sv, ", "))")
    Ω     = Diagonal(1 ./ se0) * (W' * W) * Diagonal(1 ./ se0)
    G     = refit_matrix(JA, Ω)
    nf, ne = length(free), length(ext)
    se_e  = Float64[e.se for e in ext]
    Σ     = [[data_cov(dat, corr, a.key, b.key, k) for a in act, b in act] for k in ks]
    all(S -> isposdef(Symmetric(S)), Σ) || error("a target covariance matrix is not positive definite")
    Vθ    = [Symmetric(G * S * G') for S in Σ]
    dθde  = -G * dmde[rowsA, :]                       # refit shift of ϑ̂ per unit of each external
    dmtot = dmde .+ J * dθde                         # total effect of each external on every moment
    se_θ_T = [sqrt.(max.(diag(V), 0.0)) for V in Vθ]
    all(s -> all(x -> isfinite(x) && x > 0, s), se_θ_T) ||
        error("a parameter SE is not positive and finite: $(join(se_θ_T[1], ", "))")
    se_θ_E = vec(sqrt.(sum((dθde .* se_e') .^ 2; dims = 2)))
    se_θ   = [sqrt.(s .^ 2 .+ se_θ_E .^ 2) for s in se_θ_T]
    se_m_T = [sqrt.(max.(vec(sum((J * V) .* J; dims = 2)), 0.0)) for V in Vθ]
    se_m_E = vec(sqrt.(sum((dmtot .* se_e') .^ 2; dims = 2)))
    se_m   = [sqrt.(s .^ 2 .+ se_m_E .^ 2) for s in se_m_T]
    sdk    = [[t.se * (t.key in LITERATURE_KEYS ? k : 1.0) for t in act] for k in ks]   # the targets' scales used
    vrows  = map(rows) do r
        ii   = [idx[string(key)] for (key, _) in r.comps]
        c    = [cf for (_, cf) in r.comps]
        Jv   = vec(sum(c[i] .* J[ii[i], :] for i in eachindex(ii)))
        dv   = vec(sum(c[i] .* dmtot[ii[i], :] for i in eachindex(ii)))
        mv   = sum(c .* m0[ii])
        dval = sum(cf * dat[key].value for (key, cf) in r.comps)
        vd   = [sum(ca * cb * data_cov(dat, corr, a, b, k) for (a, ca) in r.comps, (b, cb) in r.comps) for k in ks]
        cv   = [dot(vec(G' * Jv), [sum(cf * data_cov(dat, corr, t.key, key, k) for (key, cf) in r.comps) for t in act])
                for k in ks]
        sT   = [sqrt(max(Jv' * V * Jv, 0.0)) for V in Vθ]
        sE   = sqrt(sum((dv .* se_e) .^ 2))
        z    = [(v2 = sT[i]^2 + sE^2 + vd[i] - 2cv[i]; v2 > 0 ? (mv - dval) / sqrt(v2) : NaN) for i in eachindex(ks)]
        (; r.label, r.comps, model = mv, data = dval, scale_only = any(key in SCALE_KEYS for (key, _) in r.comps),
           se_T = sT, se_E = sE, se = [sqrt(a^2 + sE^2) for a in sT], se_d = sqrt.(vd), cov = cv, z,
           z_naive = (mv - dval) / sqrt(vd[1]), contrib_E = abs.(dv) .* se_e,
           contrib_T = [abs.(G' * Jv) .* sdk[i] for i in eachindex(ks)])
    end
    levels = [(; name = n, label = l, model = m0[idx[n]], se_T = [s[idx[n]] for s in se_m_T], se_E = se_m_E[idx[n]],
               se = [s[idx[n]] for s in se_m]) for (n, l) in LEVELS if haskey(idx, n) && isfinite(m0[idx[n]])]
    fit    = maximum(abs, [(m0[rowsA[i]] - act[i].value) / act[i].se for i in eachindex(act)])
    selfse = maximum(abs, [se_m_T[1][rowsA[i]] - sqrt(Σ[1][i, i]) for i in eachindex(act)])
    return (; θ, free, act, ks, corr, names, idx, m0, J, dmde, ext, JA, sv, G, Σ, Vθ, dθde, dmtot, se_θ_T, se_θ_E, se_θ,
            se_m_T, se_m_E, se_m, rows = vrows, levels, fit, selfse, exact = length(act) == nf)
end

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
g4(x) = isfinite(x) ? @sprintf("%.4g", x) : "—"
g3(x) = isfinite(x) ? @sprintf("%.3g", x) : "—"
zf(x) = isfinite(x) ? @sprintf("%+.2f", x) : "—"
mdrow(io, cells) = println(io, "| ", join(cells, " | "), " |")
mdtable(io, head, rows) = (mdrow(io, head); mdrow(io, fill("---", length(head))); foreach(r -> mdrow(io, r), rows))
klabel(k) = k == round(k) ? string(Int(k)) : string(k)

"""
    report_text(A, opts, meas, solver, fd; floor_tol = nothing) -> String

The Markdown report: identification, parameter SEs, target and external inputs with the refit
shifts they cause, the validation table (the first transport factor in full, the others as SE
and z), the target and external contributions to the model SEs and model levels.
"""
function report_text(A, opts, meas, solver, fd; floor_tol = nothing)
    io  = IOBuffer()
    ks  = A.ks
    pn  = string.(THETA_NAMES[A.free])
    nf  = length(A.free)
    println(io, "# Delta-method uncertainty at ϑ̂\n")
    println(io, "`spatial_uncertainty.jl`, ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()), ". Options: ",
            join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                  if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "), ".  \n",
            "Measurement: ", meas, ".  \nSolver: ", solver, floor_tol === nothing ? "" : ", floor tolerance $floor_tol", ".  \n",
            fd.solves > 0 ? @sprintf("%d solves on %d worker(s): %.1f minutes.", fd.solves, nworkers(), fd.minutes) :
                            "Jacobians reused from --reuse (no solves).", "\n")

    println(io, "## Identification\n")
    println(io, "- ", length(A.act), " active targets, ", nf, " free parameters (",
            A.exact ? "exactly identified: G = J_A⁻¹" : "overidentified: G = (J_A′ΩJ_A)⁻¹J_A′Ω, Ω from --weight $(opts["weight"])", ").")
    println(io, "- Standardized J_A (rows by target scale, columns by box width) has singular values ", join((g3(s) for s in A.sv), ", "),
            @sprintf("; condition number %.3g.", A.sv[1] / A.sv[end]))
    println(io, @sprintf("- Fit at ϑ̂: max |(m − t)/scale| over the active targets = %.2e.", A.fit),
            @sprintf(" Model SE of each active target from J_A G Σ_A G′ J_A′ against its scale: max abs difference %.2e%s.", A.selfse,
                     A.exact ? " (zero when exactly identified)" : " (not zero when overidentified)"))
    A.exact && A.fit > 1e-4 &&
        println(io, "- ⚠ ϑ does not fit the active targets exactly (max residual above): the SEs and z are the linearization at this point, not at an estimate; refit it first (`spatial_diagnose.jl --polish`).")
    println(io, "- Transport factors ", join(klabel.(ks), ", "), " (literature targets ", join(string.(LITERATURE_KEYS), ", "),
            "; the first is --transport and heads the tables).")
    println(io, "- Whole-CZ bootstrap correlations of the district moments in Σ_A and in cov (scheme ", something(opts["district"], "locale"),
            "; ", A.corr.source, "), below.\n")
    used = Symbol[t.key for t in A.act]
    foreach(r -> append!(used, first.(r.comps)), A.rows)
    ci = [i for (i, a) in enumerate(A.corr.keys) if a in used]
    mdtable(io, vcat([""], string.(A.corr.keys[ci])),
            [vcat([string(A.corr.keys[i])], [@sprintf("%+.2f", A.corr.matrix[i, j]) for j in ci]) for i in ci])
    println(io)

    println(io, "## Parameters\n")
    head = vcat(["parameter", "estimate", "SE externals"], vcat([["SE targets k=$(klabel(k))", "SE total k=$(klabel(k))"] for k in ks]...))
    mdtable(io, head, [vcat([pn[j], g4(A.θ[A.free[j]]), g3(A.se_θ_E[j])],
                            vcat([[g3(A.se_θ_T[i][j]), g3(A.se_θ[i][j])] for i in eachindex(ks)]...)) for j in 1:nf])
    println(io, "\nTargets' sampling or scale error and external uncertainty are independent; SE total = √(SE targets² + SE externals²).\n")

    println(io, "## Targets, externals and the refit\n")
    kt = ks[1]
    mdtable(io, ["target", "value", "scale", "scale used (k=$(klabel(kt)))", "Δϑ̂ per +1 scale: " * join(pn, ", ")],
            [[string(t.key) * (t.key in SCALE_KEYS ? " †" : ""), g4(t.value), g3(t.se),
              g3(t.se * (t.key in LITERATURE_KEYS ? kt : 1.0)),
              join((@sprintf("%+.3g", A.G[j, a] * t.se * (t.key in LITERATURE_KEYS ? kt : 1.0)) for j in 1:nf), ", ")]
             for (a, t) in enumerate(A.act)])
    println(io, "\n† the scale is a modelling tolerance (10% of a teaching share, 5% of the teacher 90/10 or the schooling ratio, 0.05 for the move rate and wage gaps; M3: 5% of a BA share or of the pooled 90/10, 0.01 for the BA share of teachers and the teachers' schooling, 0.02 for the BA premium, 0.005 for the non-teachers' schooling, 0.05 for the teacher premia and BA-share gaps), not a sampling SE; the parameter SEs inherit it. The district correlations enter Σ_A (gap_pupils, gap_salary).\n")
    if !isempty(A.ext)
        mdtable(io, ["external", "value", "SE", "step", "Δϑ̂ per +1 SE: " * join(pn, ", "), "basis"],
                [[string(e.name), g4(e.value), g3(e.se), g3(e.step),
                  join((@sprintf("%+.3g", A.dθde[j, k] * e.se) for j in 1:nf), ", "), e.basis] for (k, e) in enumerate(A.ext)])
        println(io)
    end

    println(io, "## Validation moments\n")
    head = vcat(["moment", "model", "SE targets", "SE externals", "SE model", "data", "SE data", "ρ(m, d)", "z (naive)", "z"],
                vcat([["SE model k=$(klabel(k))", "z k=$(klabel(k))"] for k in ks[2:end]]...))
    linked(r) = get(FIRST_STAGE_LINK, r.label, nothing) in (e.name for e in A.ext)
    rho(r, i) = r.cov[i] == 0 ? "0" : @sprintf("%+.2f", r.cov[i] / (r.se[i] * r.se_d[i]))
    mdtable(io, head, [vcat([string(r.label) * (r.scale_only ? " †" : linked(r) ? " ‡" : ""), g4(r.model),
                             g3(r.se_T[1]), g3(r.se_E), g3(r.se[1]), g4(r.data), g3(r.se_d[1]), rho(r, 1), zf(r.z_naive), zf(r.z[1])],
                            vcat([[g3(r.se[i]), zf(r.z[i])] for i in 2:length(ks)]...)) for r in A.rows])
    corrd = [string(r.label) for r in A.rows if r.cov[1] != 0]
    println(io, "\nk = $(klabel(ks[1])) in the first block. z (naive) = (model − data)/SE data; z = (model − data)/√(SE model² + SE data² − 2 cov), ",
            "cov = J_v G Σ_{A,v} and ρ(m, d) = cov/(SE model × SE data) at k = $(klabel(ks[1])). cov is nonzero for ",
            isempty(corrd) ? "no moment" : join(corrd, ", "), ": the data moments that share the whole-CZ bootstrap draws with the district targets ",
            "(correlations above); the NLSY, ACS and literature moments come from other samples and have cov = 0.")
    println(io, "† the data SE is a modelling scale, not a sampling SE.", any(linked, A.rows) ?
            " ‡ the data moment is the first-stage input of the external s_z or ρz: independence overstates the SE of the model–data difference, so z is a conservative consistency check, not an independent validation." : "")
    println(io, "payroll_pp_gap = gap_salary + gap_teachers_pp (log payroll per pupil), data SE from the bootstrap correlation of the two; at an exact fit its z equals the FTE gap's (the matched salary parts cancel).\n")

    println(io, "### SE of the model value contributed by each active target (|J_v G[:, a]| × scale used, k = $(klabel(ks[1])))\n")
    mdtable(io, vcat(["moment"], [string(t.key) * (t.key in SCALE_KEYS ? " †" : "") for t in A.act]),
            [vcat([string(r.label)], [g3(c) for c in r.contrib_T[1]]) for r in A.rows])
    println(io, "\nOne target at a time: they add in quadrature only if the targets are uncorrelated (the district targets are correlated, see above), ",
            "so SE targets is not their root sum of squares. † columns are tolerances, not sampling SEs: the SE of a model value that rests on them is a conservative scale, not a sampling error.\n")

    if !isempty(A.ext)
        println(io, "### SE of the model value contributed by each external (|dm/de| × SE)\n")
        mdtable(io, vcat(["moment"], [string(e.name) for e in A.ext]),
                [vcat([string(r.label)], [g3(c) for c in r.contrib_E]) for r in A.rows])
        println(io)
    end

    if !isempty(A.levels)
        println(io, "## Model levels\n")
        mdtable(io, ["moment", "model", "SE targets", "SE externals", "SE model"],
                [[l.label * " (" * l.name * ")", g4(l.model), g3(l.se_T[1]), g3(l.se_E), g3(l.se[1])] for l in A.levels])
        println(io)
    end
    println(io, "## Notes\n")
    println(io, "- First order: central differences at ϑ̂; refits at e ± SE (`spatial_diagnose.jl --external … --polish`) check it. For σ (SE 0.30) the refit of r_m and λ is strongly nonlinear: ",
            "at exact-polish-base the refits at σ = 0.30 / 0.85 (−0.9 / +1.0 SE) move r_m by +0.504 / −0.153 and λ by −0.088 / +0.362 against the linear +0.252 / −0.281 and −0.153 / +0.171, ",
            "and σ = 0.15 (−1.4 SE) has no exact fit (r_m at its bound), so the σ-driven SEs of log κ̃, λ and r_m are local and one tail is about twice the linear SE; ",
            "σϵ and s_z agree with refits within a few percent at one SD.")
    println(io, "- The scales of the ACS targets and of the move rate (†) are tolerances, not sampling SEs, so the SEs of the parameters they pin and of the model values that rest on them are conservative scales (see the target contributions); `--target-se` overrides them and `--reuse` redoes the algebra without solves.")
    println(io, "- The externals are independent of each other and of the targets. σϵ's SD is a specification range; the first stage in fact ties σϵ̂ to s_z (a lower s_z raises σϵ̂).")
    println(io, "- The covariance of the data moments with the targets uses the correlations of the whole-CZ bootstrap draws of the scheme (the four district targets, the poverty, SEDA growth and score gaps); moments from one NLSY79 sample (move_rate, the city-share gradients, the teacher AFQT gap, ρz, s_z) and the ACS and literature moments are treated as independent of the targets and of each other.")
    return String(take!(io))
end

"""
    uncertainty_dict(A, opts, dat) -> Dict

Every number of the analysis for `uncertainty.toml`. Quantities that depend on the transport factor
are lists over `transport`, in the order of the factors. `parameters` and `refit` follow
the free parameters, `targets` the active targets, `moments` every flattened moment, and
`validation.<key>` the validation rows (`contrib_externals` follows `externals.names`, `contrib_targets`
`targets.names` and is a list over `transport`); `correlations` holds the bootstrap correlation matrix used.
"""
function uncertainty_dict(A, opts, dat; fingerprint = "", floor_tol = nothing)
    nf = length(A.free)
    byk(f) = [f(i) for i in eachindex(A.ks)]
    d = Dict{String,Any}(
        "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()), "script" => "spatial_uncertainty.jl",
        "theta_spec" => opts["theta"], "transport" => collect(A.ks), "exactly_identified" => A.exact,
        "fingerprint" => fingerprint, "floor_tol" => something(floor_tol, NaN),
        "correlations" => Dict{String,Any}("keys" => string.(A.corr.keys), "source" => A.corr.source,
                                           "matrix" => [A.corr.matrix[i, :] for i in axes(A.corr.matrix, 1)]),
        "options" => join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                           if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "),
        "fit_max_abs_residual" => A.fit, "singular_values_standardized" => A.sv, "condition_standardized" => A.sv[1] / A.sv[end],
        "parameters" => Dict{String,Any}("names" => collect(string.(THETA_NAMES[A.free])), "estimate" => A.θ[A.free],
                                         "se_externals" => A.se_θ_E, "se_targets" => byk(i -> A.se_θ_T[i]),
                                         "se" => byk(i -> A.se_θ[i])),
        "refit" => Dict{String,Any}("targets" => [string(t.key) for t in A.act], "G" => [A.G[j, :] for j in 1:nf],
                                    "dtheta_de" => [A.dθde[j, :] for j in 1:nf],
                                    "cov_theta_targets" => byk(i -> [Matrix(A.Vθ[i])[j, :] for j in 1:nf])),
        "targets" => Dict{String,Any}("names" => [string(t.key) for t in A.act], "value" => [t.value for t in A.act],
                                      "scale" => [t.se for t in A.act],
                                      "cov" => byk(i -> [A.Σ[i][a, :] for a in eachindex(A.act)])),
        "externals" => Dict{String,Any}("names" => [string(e.name) for e in A.ext], "value" => [e.value for e in A.ext],
                                        "se" => [e.se for e in A.ext], "step" => [e.step for e in A.ext],
                                        "basis" => [e.basis for e in A.ext]),
        "moments" => Dict{String,Any}("names" => A.names, "model" => A.m0, "se_externals" => A.se_m_E,
                                      "se_targets" => byk(i -> A.se_m_T[i]), "se" => byk(i -> A.se_m[i])),
        "validation" => Dict{String,Any}(string(r.label) => Dict{String,Any}(
                "model" => r.model, "data" => r.data, "data_se_is_scale" => r.scale_only,
                "se_targets" => r.se_T, "se_externals" => r.se_E, "se_model" => r.se, "se_data" => r.se_d,
                "cov_model_data" => r.cov, "z" => r.z, "z_naive" => r.z_naive, "contrib_externals" => r.contrib_E,
                "contrib_targets" => r.contrib_T)
            for r in A.rows))
    return d
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
function run_uncertainty(opts)
    ref = frozen_reference(something(opts["reference"], REFERENCE_FILE), SOLVER)
    setup = problem_setup(opts, ref)
    (; base, meas, targets, active, bounds, W) = setup
    opts["floor-tol"] === nothing || hasfield(SpatialModel, :floor_tol) ||
        error("--floor-tol: this model code has no consumption-floor tolerance")
    floor_tol = get(setup, :floor_tol, nothing)             # older model code: no tolerance in the setup
    grid   = NamedTuple(k => (k in (:Nz, :nϵT, :nXO, :maxit, :hh_maxit) ? Int(v) : v) for (k, v) in opts["grid"])
    solver = merge(SOLVER, grid)
    model  = spatial_model(base, [t.key for t in active]; meas, solver, floor_tol)
    θ = read_theta(opts["theta"])
    for (name, v) in opts["fix"]
        θ[theta_index(name)] = v
    end
    external_start!(θ, bounds, base)       # a freed s_z, σ_ϵ or θ_ν at its off value 0 starts at the external block's value (as spatial_diagnose.jl)
    free   = findall(first.(bounds) .< last.(bounds))
    scheme = something(opts["district"], "locale")
    corr   = load_draw_corr(scheme)
    dat    = data_moments(targets, scheme)
    ext    = external_inputs(opts, base, meas; θ)
    ks     = vcat(opts["transport"], [k for k in (1.0, 2.0) if k != opts["transport"]])
    out    = opts["out"]
    out === nothing || mkpath(out)
    println("spatial_uncertainty.jl  ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()), "\n  free parameters: ",
            join(THETA_NAMES[free], ", "), "\n  active targets: ", join((t.key for t in active), ", "),
            "\n  externals: ", join(("$(e.name) = $(e.value) ± $(e.se)" for e in ext), ", "), "\n  workers: ", nworkers())

    fp = spec_fingerprint(model, bounds)
    fd = opts["reuse"] === nothing ? finite_differences(opts, ref, model, θ, free, bounds, ext) :
         load_jacobian(joinpath(opts["reuse"], "jacobian.toml"), θ, free, fd_steps(bounds), ext, fp)
    out !== nothing && opts["reuse"] === nothing && save_jacobian(joinpath(out, "jacobian.toml"), θ, free, fd, ext, fp)
    rows = validation_rows(targets, dat, fd.names, fd.m0)
    A = analyze(; θ, free, bounds, act = active, W, fd.names, fd.m0, fd.J, fd.dmde, ext, corr, dat, ks, rows)
    ft = hasfield(SpatialModel, :floor_tol) ? model.floor_tol : nothing
    text = report_text(A, opts, meas, solver, fd; floor_tol = ft)
    println(text)
    if out !== nothing
        write(joinpath(out, "report.md"), text)
        open(io -> TOML.print(io, uncertainty_dict(A, opts, dat; fingerprint = fp, floor_tol = ft); sorted = true),
             joinpath(out, "uncertainty.toml"), "w")
        println("Report, uncertainty.toml", opts["reuse"] === nothing ? " and jacobian.toml" : "", " written to ", out)
    end
end

run_uncertainty(OPTS)
