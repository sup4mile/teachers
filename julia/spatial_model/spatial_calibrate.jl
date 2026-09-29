# =============================================================================
# Calibration harness for the spatial block of data/spatial/spatial_calibration.md.
#
# The internal parameters ϑ = (log κ̄, δκ, ΔB, β, λ, r_m) are fitted to the seven
# moments of §2, conditional on the external block of §1:
#     κ = (κ̄, κ̄ e^{δκ}),  B = (0, ΔB),  m₁₂ = m₂₁ = r_m C̄,
# with C̄ and h̄ frozen at a reference solution (§1, "Units"). The external block
# has 21 occupations: 19 market occupations, home production and K–12 teaching.
#
# Layout
#   1  External inputs (§1), the occupational block (Table 2) and the calculations behind them
#   2  Frozen reference levels C̄, h̄
#   3  Targets (§2)
#   4  The internal parameter vector ϑ
#   5  Model evaluation ϑ → moments (the callable handed to TikTak)
#   6  Serial diagnostics: moment table, criterion, moment Jacobian (T6)
#
# Nothing here depends on TikTak: `spatial_estimate.jl` wraps `SpatialModel` in a
# TikTak `MomentObjective` and runs the global search.
# =============================================================================

isdefined(@__MODULE__, :model_moments) || include(joinpath(@__DIR__, "spatial_moments.jl"))

using TOML
using SHA

# -----------------------------------------------------------------------------
# 1. External inputs
#
# Table 1 values; * marks those still provisional in §1. The ability rows (σϵ, ρz,
# s_z) and the occupational block are placeholders until Table 2 (T2) delivers
# estimates. ωf is internal in Table 3 but not yet part of ϑ.
# -----------------------------------------------------------------------------

"""
Occupations of the calibration (§1, Table 2): the 19 market occupations of the ACS
file in its order, home production (`HP_OCC`) and K–12 teaching (`T_OCC`).
`men`/`women` are ACS 2009–13 weighted counts, ages 25–34, from
data/LaborMarketData/wages_occ_shares_v2.xlsx, sheet `moments_shares`. Military
and unemployed are dropped; 15–29-hour workers are split equally between home
production and their occupation.
"""
const ACS_OCCUPATIONS = (;
    names = ["Executives, Administrative, and Managerial", "Management Related",
             "Architects, Engineers, Math, and Computer Science",
             "Natural and Social Scientists, Recreation, Religious, Arts, Athletes",
             "Doctors and Lawyers", "Nurses, Therapists, and Other Health Service",
             "Teachers, Postsecondary", "Teachers, Non-Postsecondary and Librarians",
             "Health and Science Technicians", "Sales, All", "Administrative Support, Clerks, Record",
             "Fire, Police, and Guards", "Food, Cleaning, and Personal Services and Private Household",
             "Farm, Related Agriculture, Logging, and Extraction", "Mechanics and Construction",
             "Precision Manufacturing", "Manufacturing Operators",
             "Fabricators, Inspectors, and Material Handlers", "Vehicle Operators",
             "Home Production", "Kindergarten - Secondary Teachers"],
    men   = [54751, 24884, 36021.5, 24503, 10932, 12358, 6059, 3596, 30792, 61715.5, 48822,
             25086, 45172, 23031, 93691, 12867, 20492.5, 33825.5, 32518, 178848.5, 15468.5],
    women = [49517.5, 31581, 10938, 29176, 10617, 73932, 6163.5, 12179, 26232.5, 53598, 110723,
             6224.5, 58229.5, 3666, 2554, 4041.5, 6063.5, 8800, 2926, 279887, 50120.5])

const N_OCC  = length(ACS_OCCUPATIONS.names)
const HP_OCC = N_OCC - 1
const T_OCC  = N_OCC

"ACS counts of gender g (1 = men, 2 = women), all 21 occupations."
acs_counts(g) = g == 1 ? ACS_OCCUPATIONS.men : ACS_OCCUPATIONS.women

"Shares over the 20 non-teaching occupations: the Table 2 share targets."
function occ_shares(g)
    c = acs_counts(g)[1:N_OCC .!= T_OCC]
    return c ./ sum(c)
end

"K–12 teaching share over all 25–34-year-olds of gender g, home production included (Table 3)."
teach_share_data(g) = acs_counts(g)[T_OCC] / sum(acs_counts(g))

"""
    roy_shares(logΘ, s)

Occupation shares when X_i = Θ_i ϵ_i^α and log X_i ~ N(log Θ_i + c, s²) iid, s = ασϵ:
P_i = ∫ φ(u) ∏_{j≠i} Φ(u + (log Θ_i − log Θ_j)/s) du.
"""
function roy_shares(logΘ, s)
    n = length(logΘ)
    share(i) = quadgk(u -> pdf(Normal(), u) *
                           prod((cdf(Normal(), u + (logΘ[i] - logΘ[j]) / s) for j in 1:n if j != i); init = 1.0),
                      -10.0, 10.0; rtol = 1e-10)[1]
    return [share(i) for i in 1:n]
end

"""
    invert_roy_shares(shares, s; ref = length(shares), tol = 1e-10, maxit = 10_000)

log Θ reproducing `shares` in `roy_shares`, normalized to log Θ_ref = 0, by the damped
fixed point log Θ ← log Θ + 0.7 s (log shares − log model shares).
"""
function invert_roy_shares(shares, s; ref = length(shares), tol = 1e-10, maxit = 10_000)
    logΘ = zeros(length(shares))
    for _ in 1:maxit
        d = log.(shares) .- log.(roy_shares(logΘ, s))
        maximum(abs, d) < tol && return logΘ
        logΘ .+= 0.7 * s .* d
        logΘ .-= logΘ[ref]
    end
    error("share inversion did not converge in $maxit iterations at s = $s")
end

"""
    occupation_block(; σϵ, α) -> (; A, r)

Provisional Table 2 values at ability dispersion ασϵ, by inverting the static Roy
model over the 20 non-teaching occupations with teaching ignored:
- `A`: A_i/A_HP from men's shares (length `N_OCC`, A_HP = 1, A[T] = NaN);
- `r`: relative female wedges r_i = Θ_{i,f}/(ωf A_i) from women's shares, r_HP = 1,
  so that 1 − τω_{i,f} = ωf r_i (`female_wedges`).
Placeholders until T2 fits σϵ jointly with the ability block and corrects women's
shares for selection into teaching (§1, Table 2).
"""
function occupation_block(; σϵ, α)
    s  = α * σϵ
    lm = invert_roy_shares(occ_shares(1), s; ref = HP_OCC)
    lf = invert_roy_shares(occ_shares(2), s; ref = HP_OCC)
    return (; A = vcat(exp.(lm), NaN), r = vcat(exp.(lf .- lm), NaN))
end

"τω with no male or teaching wedges and 1 − τω_{i,f} = ωf r_i for i ≠ T."
function female_wedges(r, ωf)
    τω = zeros(N_OCC, 2)
    for i in 1:N_OCC
        i == T_OCC || (τω[i, 2] = 1 - ωf * r[i])
    end
    return τω
end

const EXTERNAL = (;
    σ   = 0.56,     # * class-size curvature: FOO (2013) wage mapping at η = 0.080, K = 12
    θν  = 0.25,     # * relative location-taste scale σν/μ: Eckert–Kleineberg (2024)
    ψ   = 0.0,      # warm-glow curvature: linear warm glow
    η   = 0.080,    # goods elasticity: 2018 spending net of teacher pay (`eta_2018`)
    φ   = 1.12,     # * time elasticity: s_O = 0.55 at μ = 1 (placeholder years, T2)
    μ   = 1.0,      # weight on log consumption: normalization
    γ   = 0.83,     # teacher wage curvature: starting value; internal (Table 3)
    α   = 1.0,      # ability elasticity: normalization
    σϵ  = 0.46,     # placeholder (T2): pooled 90/10 net of s_z, no selection (`sigma_eps_placeholder`)
    ρz  = 0.9,      # placeholder (T2): old code value
    s_z = 0.104,    # placeholder (T2): stationary log-z SD at b s_z = χ, i.e. χ(1 − η)
    ωf  = 0.934)    # placeholder (Table 3, internal): female level, 6.0% female teaching share at ϑ₀

"""
    external_params(; kwargs...) -> Params

The fixed block of §1 as a 21-occupation, two-location `Params`: 19 market
occupations, home production (A_HP = 1) and K–12 teaching (`T_OCC`). Keywords
override entries of `EXTERNAL` or any other `Params` field. Unless `A` or `τω` is
passed, the occupational block comes from `occupation_block` at the given σϵ and α,
with female wedges scaled by the level ωf. The internal fields (κ, B, β, λ, mcost)
are placeholders that `theta_to_params` overwrites, and `Cbar`/`href` are replaced
by the frozen reference levels (`calibration_base`). The location-taste scale
enters as θν, so σν = θν μ, and ability as its stationary SD s_z, so σξ = s_z √(1 − ρz²).
"""
function external_params(; kwargs...)
    x = merge(EXTERNAL, NamedTuple(kwargs))
    occ = haskey(x, :A) && haskey(x, :τω) ? nothing : occupation_block(; x.σϵ, x.α)
    A  = haskey(x, :A)  ? x.A  : occ.A
    τω = haskey(x, :τω) ? x.τω : female_wedges(occ.r, x.ωf)
    rest = Base.structdiff(x, NamedTuple{(:θν, :s_z, :ωf, :A, :τω)})
    return Params(N_OCC, 2; T = T_OCC, A, τω, τmove_off = 0.0, mcost_off = 0.0, rest...,
                  σν = taste_scale(x.μ, x.θν), σξ = innovation_sd(x.s_z, x.ρz))
end

"Location-taste scale σν = θν μ for the relative scale θν = σν/μ."
taste_scale(μ, θν = 0.25) = θν * μ

"AR(1) innovation SD for a stationary log-ability SD s_z: σξ = s_z √(1 − ρz²)."
innovation_sd(s_z, ρz) = s_z * sqrt(1 - ρz^2)

"""
    eta_2018(; spend_gdp, labshare, teachers, salary, benefits, gdp) -> (; S, t, η)

§1 goods elasticity from 2018 aggregates. Spending over labour income is
S = η(1 − t) + t, where t is the teacher-payroll tax, so η = (S − t)/(1 − t).
- `spend_gdp`: education institutions, all levels, 6.0% of GDP (Digest 2021, Table 605.20);
- `labshare`: PWT 11 labour share, 0.5906 (FRED LABSHPUSA156NRUG);
- `teachers`: 3.2 million public FTE teachers, fall 2018 (Digest 2018);
- `salary`: average public teacher salary 2017–18, \$60,483 (Digest 2018, Table 211.60);
- `benefits`: instruction benefits over salaries 2017–18, 103.8/241.5 (Digest 2020, Table 236.20);
- `gdp`: nominal GDP, \$20,656.5 billion (BEA, FRED GDPA).
Pass `payroll` to replace teachers × salary × (1 + benefits), e.g. all instruction
salaries and benefits.
"""
function eta_2018(; spend_gdp = 0.060, labshare = 0.5906, teachers = 3.2e6, salary = 60_483.0,
                  benefits = 103.821 / 241.516, gdp = 20_656.516e9,
                  payroll = teachers * salary * (1 + benefits))
    S = spend_gdp / labshare
    t = payroll / (labshare * gdp)
    return (; S, t, η = (S - t) / (1 - t))
end

"""
    sigma_eps_placeholder(; p9010, η, s_z) -> σϵ

Placeholder ability dispersion: the pooled within-occupation 90/10 of hourly wages
(3.737, non-teachers, ACS 2009–13) read as lognormal with log-wage SD
√(σϵ² + s_z²)/(1 − η) at α = 1, ignoring selection within occupations. Selection
lowers within-occupation dispersion, so this is a lower bound; T2 replaces it.
"""
sigma_eps_placeholder(; p9010 = 3.7374, η = 0.080, s_z = 0.104) =
    sqrt((log(p9010) * (1 - η) / (2 * quantile(Normal(), 0.9)))^2 - s_z^2)

"""
    sigma_from_wages(; η, K, d, N̄, b, se_b) -> (σ, se)

FOO (2013) direct-wage mapping σ = (1−η)(K/d) N̄ b for a log-wage effect −b per
pupil over d school years at mean class size N̄ (§3).
"""
function sigma_from_wages(; η = 0.080, K = 12.0, d = 3.0, N̄ = 24.357, b = 0.0063, se_b = 0.0033)
    factor = (1 - η) * K / d * N̄
    return (σ = factor * b, se = factor * se_b)
end

"STAR first-year cross-check σ = (1−η) χ d_S K / (d log(N₀/N₁)) (§3)."
sigma_from_star(; η = 0.080, χ = log(1.12), K = 12.0, dS = 0.20, N0 = 22.4, N1 = 15.1, d = 1.0) =
    (1 - η) * χ * dS * K / (d * log(N0 / N1))

"Zero-goods-cost time shares (s_O, s_T) and s_T/s_O (§1)."
function schooling_shares(; μ, φ, η, γ)
    sO = μ * φ / (μ * φ + 1 - η)
    sT = μ * φ * γ / (μ * φ * γ + 1 - γ * η)
    return (; sO, sT, ratio = sT / sO)
end

"Time elasticity φ matching the non-teacher time share s_O at zero goods cost (§1)."
phi_from_schooling(sO; μ = 1.0, η = 0.080) = sO * (1 - η) / (μ * (1 - sO))

"Analytic CFR mapping (§3): Δlog y ≈ β sd / (K(1−η)) for non-teachers without goods costs."
cfr_analytic(β, sd; K = 12.0, η = 0.080) = β * sd / (K * (1 - η))

"Invert `cfr_analytic` for β: an initialization, not an estimate."
beta_from_cfr(effect, sd; K = 12.0, η = 0.080) = effect * K * (1 - η) / sd

"Move probability at equal attributes and gross consumption, goods ratio r_m = m/C (audit §3)."
move_probability(r_m, θν = 0.25) = 1 / (1 + exp(-log1p(-r_m) / θν))

"BFM WTP as a consumption share: \$/month over mean monthly consumption, same price year (T3)."
wtp_share_target(wtp_month, consumption_month) = wtp_month / consumption_month

"""
    external_report(io = stdout; x = EXTERNAL, K = 12.0, χ = log(1.12))

Print the intermediate calculations behind §1 at the external block `x`. The
note's values are in brackets.
"""
function external_report(io::IO = stdout; x = EXTERNAL, K = 12.0, χ = log(1.12))
    e  = eta_2018()
    ea = eta_2018(; payroll = (241.516 + 103.821) * 1e9)
    w  = sigma_from_wages(; η = x.η, K)
    sh = schooling_shares(; x.μ, x.φ, x.η, x.γ)
    occ = occupation_block(; x.σϵ, x.α)
    mkt = [i for i in 1:N_OCC if i ∉ (HP_OCC, T_OCC)]
    println(io, "\n===== External inputs: intermediate calculations =====")
    @printf(io, "  η  2018: S = 0.060/0.5906 = %.4f, t = %.4f, (S − t)/(1 − t) = %.4f   [0.102, 0.023, 0.080]\n",
            e.S, e.t, e.η)
    @printf(io, "     all instruction salaries and benefits: t = %.4f, η = %.4f   [0.030, 0.073]\n", ea.t, ea.η)
    @printf(io, "  σ  FOO wage mapping (η=%.3f, K=%g)       = %.4f (SE %.4f), 95%% [%.3f, %.3f]   [0.565; −0.015, 1.144]\n",
            x.η, K, w.σ, w.se, w.σ - 1.96 * w.se, w.σ + 1.96 * w.se)
    @printf(io, "     σ at K = 5: %.4f   [0.235];  STAR first-year check at χ = %.4f: %.4f   [0.687]\n",
            sigma_from_wages(; η = x.η, K = 5.0).σ, χ, sigma_from_star(; η = x.η, χ, K))
    @printf(io, "  σν = θν μ = %.2f × %.3f                  = %.4f\n", x.θν, x.μ, taste_scale(x.μ, x.θν))
    @printf(io, "  φ  matching s_O = 13.75/25 at μ = %.2f   = %.4f   [1.12]\n", x.μ,
            phi_from_schooling(13.75 / 25; x.μ, x.η))
    @printf(io, "  s_O, s_T at zero goods cost              = %.4f, %.4f;  s_T/s_O = %.4f   [0.55, 0.50; 0.91 vs. 1.23]\n",
            sh.sO, sh.sT, sh.ratio)
    @printf(io, "  χ = log(1.12) = %.4f;  s_z = χ(1 − η) = %.4f;  σξ at ρz = %.2f: %.4f\n",
            χ, χ * (1 - x.η), x.ρz, innovation_sd(x.s_z, x.ρz))
    @printf(io, "  σϵ placeholder from the pooled 90/10 (3.737) at s_z = %.3f: %.4f\n",
            x.s_z, sigma_eps_placeholder(; η = x.η, x.s_z))
    @printf(io, "  CFR target log(1.013) = %.5f;  β ≈ %.4f / sd(log h_T) by the analytic map\n",
            log(1.013), beta_from_cfr(log(1.013), 1.0; K, η = x.η))
    @printf(io, "  occupations: %d (%d market, home production = %d, K–12 teaching = %d)\n",
            N_OCC, length(mkt), HP_OCC, T_OCC)
    @printf(io, "     teaching shares: men %.4f, women %.4f   [1.9%%, 6.0%%];  home production: men %.4f, women %.4f\n",
            teach_share_data(1), teach_share_data(2), occ_shares(1)[HP_OCC], occ_shares(2)[HP_OCC])
    @printf(io, "     at ασϵ = %.3f: log A_i/A_HP ∈ [%.3f, %.3f];  log r_i ∈ [%.3f, %.3f];  female level ωf = %.4f\n",
            x.α * x.σϵ, extrema(log.(occ.A[mkt]))..., extrema(log.(occ.r[mkt]))..., x.ωf)
    print(io, "  move probability at θν = $(x.θν), r_m ∈ {0, .05, .10, .1813, .30}: ")
    println(io, join((@sprintf("%.3f", move_probability(r, x.θν)) for r in (0.0, 0.05, 0.10, 0.1813, 0.30)), ", "))
    println(io, "======================================================")
end

# -----------------------------------------------------------------------------
# 2. Frozen reference levels
#
# The goods cost is m = r_m C̄ and the warm glow is normalized by h̄. Both are
# measured once at a reference configuration and then frozen, so r_m keeps its
# units across trial parameters, profile runs, and counterfactuals.
# -----------------------------------------------------------------------------
const REFERENCE_FILE = joinpath(@__DIR__, "spatial_reference.toml")

"""
    reference_levels(ext, θ = THETA0; solver = SOLVER, maxit = 5, rtol = 1e-4, verbose = true)

Mean consumption C̄ and mean child human capital h̄ at the reference configuration
`theta_to_params(θ, ext)`. The goods cost r_m C̄ and the warm-glow reference h̄
feed back into the solution, so iterate (C̄, h̄) to a fixed point, starting from
`ext.Cbar` and `ext.href`. Returns `(; Cbar, href, sol, err)`, where `err` is the
last relative change.
"""
function reference_levels(ext::Params, θ = THETA0; solver = SOLVER, maxit = 5, rtol = 1e-4,
                          verbose = true)
    Cbar, href = ext.Cbar, ext.href
    sol, err = nothing, Inf
    for it in 1:maxit
        p   = theta_to_params(θ, with_params(ext; Cbar, href))
        sol = solve_ge(p; solver..., verbose = false)
        sol.converged || error("reference solve did not converge at iteration $it")
        C1, h1 = mean_consumption(sol), mean_child_h(sol)
        err = max(abs(C1 / Cbar - 1), abs(h1 / href - 1))
        Cbar, href = C1, h1
        verbose && @printf("  reference %d: C̄ = %.8f  h̄ = %.8f  rel. change %.2e\n", it, Cbar, href, err)
        err < rtol && break
    end
    err < rtol || @printf("  ⚠ reference levels moved %.2e at the last iteration (rtol %.1e)\n", err, rtol)
    return (; Cbar, href, sol, err)
end

"The external block and solver settings a reference is measured at, as stored in its TOML file."
function reference_config(ext::Params, solver)
    scalars  = (:α, :φ, :η, :σ, :γ, :μ, :ψ, :σν, :ρz, :σξ, :σϵ, :Mtot)
    nonteach = [i for i in eachindex(ext.A) if i != ext.T]
    external = merge(Dict{String,Any}(string(k) => getfield(ext, k) for k in scalars),
                     Dict{String,Any}("T" => ext.T, "A" => ext.A[nonteach],
                                      "tau_omega_female" => ext.τω[nonteach, 2]))
    return (; external, solver = Dict{String,Any}(string(k) => v for (k, v) in pairs(solver)))
end

"Write frozen reference levels and the configuration they were measured at to `path`."
function save_reference(path, ref, ext::Params; θ = THETA0, solver = SOLVER)
    cfg = reference_config(ext, solver)
    d = Dict{String,Any}(
        "Cbar" => ref.Cbar, "href" => ref.href, "fixed_point_change" => ref.err,
        "theta_names" => collect(string.(THETA_NAMES)), "theta" => collect(Float64, θ),
        "solver" => cfg.solver, "external" => cfg.external,
        "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
        "note" => "Frozen C̄ and h̄ (spatial_calibration.md §1). Recompute only with a deliberate " *
                  "change to the reference configuration; r_m is in units of this C̄.")
    mkpath(dirname(path))
    open(io -> TOML.print(io, d; sorted = true), path, "w")
    return path
end

"Read frozen reference levels; returns `(; Cbar, href, meta)` with the saved TOML as `meta`."
function load_reference(path = REFERENCE_FILE)
    d = TOML.parsefile(path)
    return (; Cbar = Float64(d["Cbar"]), href = Float64(d["href"]), meta = d)
end

"""
Entries where a loaded reference was measured at a configuration other than
(`ext`, `solver`). Iteration caps and damping are skipped: they do not move a
converged solution. Numbers and vectors match to a relative 1e-10, since the
occupational block is itself computed by quadrature.
"""
function reference_mismatch(ref, ext::Params, solver)
    cfg = reference_config(ext, solver)
    numeric(x) = x isa Real || (x isa AbstractVector && all(y -> y isa Real, x))
    same(a, b) = a == b || (numeric(a) && numeric(b) && length(a) == length(b) &&
                            isapprox(Float64.(a), Float64.(b); rtol = 1e-10))
    out = String[]
    for (section, now) in (("external", cfg.external), ("solver", cfg.solver)), (k, v) in now
        k in ("maxit", "hh_maxit", "damping") && continue
        saved = get(get(ref.meta, section, Dict()), k, nothing)
        same(saved, v) || push!(out, v isa AbstractVector ? "$section.$k differs" : "$section.$k = $saved (now $v)")
    end
    return sort!(out)
end

"The calibration's fixed block: `external_params(; kwargs...)` with frozen C̄ and h̄."
calibration_base(ref; kwargs...) = with_params(external_params(; kwargs...); Cbar = ref.Cbar, href = ref.href)

# -----------------------------------------------------------------------------
# 3. Targets
#
# Scales enter the criterion as ((model − target)/se)², i.e. diagonal
# inverse-variance weights. The district SEs are cross-CZ SD/√n placeholders until
# the whole-CZ bootstrap (T5); literature scales include no transport uncertainty
# yet. Treat every scale as provisional.
# -----------------------------------------------------------------------------
"""
    Target(key, value, se, param, source, note)

A calibration target. `key` names a field of `model_moments`; `value` is NaN
while the target is pending, which drops it from the criterion; `se` is its scale
in the criterion; `param` is the principal parameter it informs (§2).
"""
struct Target
    key::Symbol
    value::Float64
    se::Float64
    param::Symbol
    source::String
    note::String
end

"The seven moments of §2, in the note's order."
function default_targets()
    cfr = log(1.013)
    return [
        Target(:male_teach_share, teach_share_data(1), 0.1 * teach_share_data(1), :logκ,
               "ACS 2009–13, men 25–34 (wages_occ_shares_v2.xlsx: moments_shares)",
               "K–12 teachers over all men, home production included; scale = 10% of value"),
        Target(:gap_salary, 0.01623, 0.1020 / sqrt(176), :δκ, "spatial_moments.json: gap_salary_real_locale",
               "CWIFT-deflated salary per FTE, 176 CZs; scale = cross-CZ SD/√n pending bootstrap (T5)"),
        Target(:gap_pupils, 0.3764, 0.9691 / sqrt(188), :ΔB, "spatial_moments.json: gap_pupils_locale",
               "enrollment, 188 CZs; scale = cross-CZ SD/√n pending bootstrap (T5)"),
        Target(:cfr_effect, cfr, 0.25 * cfr, :β, "CFR (2014b): ≈1.3% earnings at 28 per teacher-VA SD-year",
               "log(1.013); scale = 25% of value, provisional transport uncertainty"),
        Target(:wtp_share, NaN, NaN, :λ, "BFM (2007) Table 7 col. 4: \$19.70/month (SE \$7.40)",
               "PENDING: divide by mean monthly consumption in 1990 dollars (T3); see `with_wtp_target`"),
        Target(:gap_income, 0.2708, 0.2057 / sqrt(187), :r_m, "spatial_moments.json: gap_seda_lninc_locale",
               "log median household income proxy, 187 CZs; measurement map pending (T3)"),
        Target(:gap_teachers_pp, -0.01216, 0.1306 / sqrt(188), :none, "spatial_moments.json: gap_teachers_pp_locale",
               "FTE per pupil vs. model headcount per student, 188 CZs; scale = cross-CZ SD/√n (T5)"),
    ]
end

"""
    with_wtp_target(targets, consumption_month; wtp = 19.70, se = 7.40)

Fill the BFM target as a consumption share, given mean monthly consumption in the
same price year as the \$19.70 estimate (1990 dollars).
"""
function with_wtp_target(targets, consumption_month; wtp = 19.70, se = 7.40)
    return [t.key === :wtp_share ?
            Target(t.key, wtp_share_target(wtp, consumption_month), se / consumption_month, t.param, t.source,
                   @sprintf("\$%.2f / \$%.2f monthly consumption; scale from BFM's SE", wtp, consumption_month)) :
            t for t in targets]
end

"Targets with a finite value; pending ones are dropped from the criterion."
active_targets(targets) = filter(t -> isfinite(t.value) && isfinite(t.se) && t.se > 0, targets)

# -----------------------------------------------------------------------------
# 4. The internal parameter vector
# -----------------------------------------------------------------------------
const THETA_NAMES = (:logκ, :δκ, :ΔB, :β, :λ, :r_m)

"""
Starting values. κ̄ = 0.2365 gives the 1.9% male teaching share and β = 0.387 the
analytic CFR map at the reference sd(log h_T) = 0.368 (K = 12), jointly with ωf in
`EXTERNAL` for the 6.0% female share. The rest are the code's κ₂/κ₁ = 1.2,
B = (0, 0.1), λ = 0.70 and r_m = 0.1813 (§2).
"""
const THETA0 = [log(0.2365), log(0.9 / 0.75), 0.1, 0.387, 0.70, 0.1813]

"""
Search box for ϑ. §4 imposes β ∈ (0, 1), λ ≥ 0 and r_m ≥ 0; the other limits
bound the Sobol screening region. Zero is included for r_m as the mechanism-off
boundary. Near ϑ₀ the male teaching share is 0.3% at κ̄ = 0.15 and 9.6% at 0.40;
both limits solve.
Widen a limit if estimates pile up against it.
"""
const THETA_BOUNDS = [(log(0.15), log(0.40)), # log κ̄
                      (-0.5, 0.5),            # δκ
                      (-0.5, 0.5),            # ΔB
                      (0.0, 0.6),             # β
                      (0.0, 1.5),             # λ
                      (0.0, 0.5)]             # r_m

"Index of a parameter in ϑ, by name."
function theta_index(name::Symbol)
    i = findfirst(==(name), THETA_NAMES)
    i === nothing && throw(ArgumentError("unknown parameter $name; use one of $THETA_NAMES"))
    return i
end

"""
    theta_to_params(θ, base) -> Params

κ = (κ̄, κ̄ e^{δκ}), B = (0, ΔB), β, λ, and the symmetric goods cost r_m C̄, with C̄ =
`base.Cbar` frozen. There is no utility moving cost.
"""
function theta_to_params(θ, base::Params)
    logκ, δκ, ΔB, β, λ, r_m = θ
    κ̄ = exp(logκ)
    return with_params(base; κ = [κ̄, κ̄ * exp(δκ)], B = [0.0, ΔB], β, λ,
                       τmove = zeros(2, 2), mcost = offdiag(r_m * base.Cbar, 2))
end

"ϑ implied by a two-location `Params` (inverse of `theta_to_params`)."
params_to_theta(p::Params) =
    [log(p.κ[1]), log(p.κ[2] / p.κ[1]), p.B[2] - p.B[1], p.β, p.λ, p.mcost[1, 2] / p.Cbar]

# -----------------------------------------------------------------------------
# 5. Model evaluation
# -----------------------------------------------------------------------------
"""
Expected failure at a trial ϑ: GE non-convergence, a tax rate on its clamp, an
unconverged (e, s) fixed point, or choice mass at the consumption floor. The
driver records these as failed evaluations instead of stopping the search.
"""
struct CalibrationFailure <: Exception
    msg::String
end
Base.showerror(io::IO, e::CalibrationFailure) = print(io, "CalibrationFailure: ", e.msg)

"""
GE solver settings for estimation. With 21 occupations a converged solve near ϑ₀,
moments included, takes 15–30 s; `maxit` caps the time a failing point, typically
a location-emptying corner, can take.
"""
const SOLVER = (; Nz = 5, nϵT = 48, nXO = 48, damping = 0.75, tol = 1e-6, maxit = 120,
                  hh_tol = 1e-7, hh_maxit = 1000)

"""
    SpatialModel(base, keys; meas = Measurement(), solver = SOLVER, verbose = false)

Callable ϑ → moments: solves the GE at `theta_to_params(ϑ, base)` and returns the
`model_moments` fields named in `keys`, in order. Every call solves from the same
cold start, so evaluations are deterministic. Throws `CalibrationFailure` at an
infeasible or unconverged point. With `verbose`, logs one line per call.
"""
struct SpatialModel
    base::Params
    meas::Measurement
    solver::NamedTuple
    keys::Vector{Symbol}
    verbose::Bool
end

SpatialModel(base::Params, keys; meas = Measurement(), solver = SOLVER, verbose = false) =
    SpatialModel(base, meas, solver, collect(Symbol, keys), verbose)

"Largest population share allowed to choose a location with consumption at the smooth_log floor."
const FLOOR_MASS_TOL = 1e-8

"""
    solve_theta(model, θ) -> (sol, moments)

Solve at ϑ and check feasibility, or throw `CalibrationFailure`. Unaffordable
moves chosen with probability ≈ 0 are allowed: smooth_log then behaves like log,
and they become common when the frozen goods cost is large relative to income.
"""
function solve_theta(model::SpatialModel, θ)
    p   = theta_to_params(θ, model.base)
    sol = solve_ge(p; model.solver..., verbose = false)
    sol.converged || throw(CalibrationFailure("GE did not converge in $(model.solver.maxit) iterations"))
    all(0 .< sol.t .< TAX_CAP - 1e-9) || throw(CalibrationFailure("tax rate on its clamp, t = $(fmtvec(sol.t))"))
    sol.s_resid ≤ 1e-8 || throw(CalibrationFailure(@sprintf("(e, s) residual %.1e", sol.s_resid)))
    m = model_moments(sol, model.meas)
    m.floor_mass ≤ FLOOR_MASS_TOL ||
        throw(CalibrationFailure(@sprintf("population share %.1e chooses consumption at the smooth_log floor", m.floor_mass)))
    return sol, m
end

function (model::SpatialModel)(θ)
    t0 = time()
    try
        _, m = solve_theta(model, θ)
        v = [Float64(m[k]) for k in model.keys]
        model.verbose && (@printf("  ϑ = %s  %.1fs  moments %s\n", fmtvec(θ), time() - t0, fmtvec(v)); flush(stdout))
        return v
    catch e
        model.verbose && (@printf("  ϑ = %s  %.1fs  failed: %s\n", fmtvec(θ), time() - t0, sprint(showerror, e)); flush(stdout))
        rethrow()
    end
end

"""
    problem_fingerprint(model, targets, bounds) -> String

Short hash of everything that defines the estimation problem: bounds, the fixed
block (including frozen C̄, h̄), measurement choices, solver settings, targets,
and the source of the three model files. Used as TikTak's `problem_id`, so a
changed model cannot silently resume an old run.
"""
function problem_fingerprint(model::SpatialModel, targets, bounds)
    io = IOBuffer()
    print(io, THETA_NAMES, bounds, model.base, model.meas, model.solver, model.keys,
          [(t.key, t.value, t.se) for t in targets])
    for f in ("spatial_continuous.jl", "spatial_moments.jl", "spatial_calibrate.jl")
        write(io, read(joinpath(@__DIR__, f)))
    end
    return bytes2hex(sha256(take!(io)))[1:12]
end

# -----------------------------------------------------------------------------
# 6. Serial diagnostics
# -----------------------------------------------------------------------------
"Standardized residuals (model − target)/se over the active targets."
residuals(m, targets) = [(m[t.key] - t.value) / t.se for t in active_targets(targets)]

"Σ((model − target)/se)²: TikTak's `MomentObjective` with `scales = se`."
criterion(m, targets) = sum(abs2, residuals(m, targets))

"Print targets, model moments and standardized residuals; pending targets show the model value only."
function moment_table(io::IO, m, targets)
    @printf(io, "  %-18s %11s %11s %9s  %-6s %s\n", "moment", "target", "model", "(m−t)/se", "param", "source")
    for t in targets
        pending = !(isfinite(t.value) && isfinite(t.se))
        @printf(io, "  %-18s %11s %11.5f %9s  %-6s %s\n", t.key,
                pending ? "pending" : @sprintf("%.5f", t.value), m[t.key],
                pending ? "—" : @sprintf("%+.2f", (m[t.key] - t.value) / t.se), t.param, t.source)
    end
    @printf(io, "  criterion over %d active targets = %.4f\n", length(active_targets(targets)), criterion(m, targets))
end
moment_table(m, targets) = moment_table(stdout, m, targets)

"Validation outcomes and diagnostics from `model_moments` (everything not targeted)."
function diagnostics_table(io::IO, m)
    @printf(io, "  validation: SEDA growth gap %.4f SD/grade [data 0.0031 (0.0030)];  Δlog Q = %.4f\n",
            m.seda_growth_gap, m.gap_logQ)
    @printf(io, "  occupations: teaching share %.4f (male %.4f, female %.4f);  s_T/s_O = %.4f\n",
            m.teach_share, m.male_teach_share, m.female_teach_share, m.s_ratio)
    @printf(io, "  mobility: gross move rate %.4f;  P(1→2) = %.4f, P(2→1) = %.4f\n", m.move_rate, m.p12, m.p21)
    @printf(io, "  CFR pieces: sd(log h_T) = %.4f;  non-teachers %.5f;  analytic %.5f\n",
            m.teacher_sd, m.cfr_nonteach, m.cfr_analytic)
    @printf(io, "  WTP share by location %s;  students/teacher %s;  t = %s;  M = %s\n",
            fmtvec(m.wtp_share_l), fmtvec(m.students_per_teacher), fmtvec(m.t), fmtvec(m.M))
    @printf(io, "  levels: mean C = %.6f, mean h = %.6f;  stationarity dev %.1e\n",
            m.mean_C, m.mean_h, m.stationarity_dev)
    @printf(io, "  consumption support: min C = %.3e;  %d (node, l′) pairs at the floor, chosen by a population share %.1e\n",
            m.minC, m.n_binding, m.floor_mass)
end
diagnostics_table(m) = diagnostics_table(stdout, m)

"""
    evaluate(model, θ, targets; io = stdout) -> (sol, m)

Solve at ϑ and print the moment and diagnostics tables.
"""
function evaluate(model::SpatialModel, θ, targets; io::IO = stdout)
    t = @elapsed (sol, m) = solve_theta(model, θ)
    println(io, "\nϑ = (", join((@sprintf("%s = %.5f", n, x) for (n, x) in zip(THETA_NAMES, θ)), ", "), ")")
    @printf(io, "  solved in %.1fs\n", t)
    moment_table(io, m, targets)
    diagnostics_table(io, m)
    return sol, m
end

"""
    moment_jacobian(model, θ, targets; step, free = eachindex(θ), map = map)

Central-difference Jacobian of the standardized moments m_k/se_k with respect to
the free parameters, for the T6 identification checks. Columns are also scaled
by the width of the search box, so singular values compare like units.
- `sv`, `cond`: singular values and condition number of the scaled Jacobian;
- `R2[j]`: share of column j reproduced by a combination of the other columns
  (§4's projection check). With seven moments and five other columns it is
  mechanically close to 1, so also read
- `unique[j]`: the norm of the part of column j the others cannot reproduce, in
  standardized moments per box width: the effect no other parameter can mimic.
Pass `map = pmap` to spread the 2·length(free) solves across workers.
"""
function moment_jacobian(model::SpatialModel, θ, targets;
                         step = 1e-3 .* [hi - lo for (lo, hi) in THETA_BOUNDS],
                         free = eachindex(θ), map = map)
    act  = active_targets(targets)
    keys = [t.key for t in act]
    se   = [t.se for t in act]
    free = collect(free)
    shifted(j, s) = (x = collect(Float64, θ); x[j] += s * step[j]; x)
    pts  = vcat([shifted(j, +1) for j in free], [shifted(j, -1) for j in free])
    vals = map(x -> (m = solve_theta(model, x)[2]; [Float64(m[k]) for k in keys]), pts)
    nf   = length(free)
    J    = reduce(hcat, [(vals[i] .- vals[i + nf]) ./ (2 * step[j]) for (i, j) in enumerate(free)])
    Js   = (J ./ se) .* [hi - lo for (lo, hi) in THETA_BOUNDS[free]]'
    sv   = svdvals(Js)
    resid = [begin
                 others = Js[:, setdiff(1:nf, c)]
                 Js[:, c] .- others * (others \ Js[:, c])
             end for c in 1:nf]
    R2     = [1 - sum(abs2, resid[c]) / sum(abs2, Js[:, c]) for c in 1:nf]
    unique = norm.(resid)
    return (; J, Js, sv, cond = sv[1] / sv[end], R2, unique, names = collect(THETA_NAMES[free]), keys)
end

function print_jacobian(io::IO, jac)
    println(io, "\n  Standardized moment Jacobian ∂(m/se)/∂ϑ × (box width):")
    @printf(io, "  %-18s", "")
    foreach(n -> @printf(io, " %9s", n), jac.names); println(io)
    for (k, key) in enumerate(jac.keys)
        @printf(io, "  %-18s", key)
        foreach(x -> @printf(io, " %9.3g", x), jac.Js[k, :]); println(io)
    end
    println(io, "  singular values: ", join((@sprintf("%.3g", s) for s in jac.sv), ", "),
            @sprintf("   condition number %.3g", jac.cond))
    println(io, "  R² of each column on the others: ",
            join((@sprintf("%s %.3f", n, r) for (n, r) in zip(jac.names, jac.R2)), ", "))
    println(io, "  unique norm (effect the others cannot mimic): ",
            join((@sprintf("%s %.3g", n, u) for (n, u) in zip(jac.names, jac.unique)), ", "))
end
print_jacobian(jac) = print_jacobian(stdout, jac)
