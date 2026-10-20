# =============================================================================
# Calibration harness for the spatial block of julia/spatial_model/calibration/spatial_calibration.md.
#
# The internal parameters ϑ = (log κ̃, δκ, ΔB, β, λ, r_m, ωf, γ) are fitted to the
# moments of §2 (Table 3), conditional on the external block of §1:
#     κ = (κ̄, κ̄ e^{δκ}) with κ̄ = κ̃ h̄^{1−γ},  B = (0, ΔB),  m₁₂ = m₂₁ = r_m C̄,
#     1 − τω_{i,f} = ωf r_i,
# with C̄ and h̄ frozen at a reference solution (§1, "Units"). The external block
# has 21 occupations: 19 market occupations, home production and K–12 teaching.
# Four extension parameters follow the eight (calibration experiments, 2026-10): the
# suburb's recurring residence cost r_p (p = (0, r_p C̄)), the ability-graded taste
# for the suburb ζ (utils per SD of log z), a non-pecuniary cost of teaching c_T (utils)
# and its ability gradient c_Tz (utils per SD of log z). All are fixed at 0 unless freed
# (`--free r_p`), which reproduces the eight-parameter model exactly. A fifth extension, the
# Stone–Geary subsistence `csub` C̄ paid by every adult (`Params.csub`), is not in ϑ: it is set
# like any other external input (`--external csub=0.1`) and is off at its default 0. A sixth,
# the residence/workplace split (`Params.split`, `--external split=1`), has one parameter in ϑ,
# the workplace taste scale relative to the location taste, ρw = σw/σν (`--fix ρw=…` or
# `--free ρw`); the commuting time cost τk and the non-teacher job share of location 1,
# `jobs_city` (the baseline population shares if not given), are external inputs
# (`--external τk=0.03 --external jobs_city=0.45`). A freed ρw is identified by the all-worker crossing
# share (`--activate commute_share_all`) or, instead, by the teacher-supply elasticity
# (`--activate teacher_supply_elasticity`, Ransom–Sims 3.7); σw itself is not an external input.
#
# The joint education and occupation choice (`Params.college`, M3; model §3c; it composes with the
# split) extends ϑ to 22 entries: ρw (the split's parameter, read only with `split`), the credential
# premium ρ_c, the college cost c_C and its gradient c_Cz
# (per SD of log z, like c_Tz), the non-BA teachers' pay penalty ρ_T0, the salary-schedule
# elasticity ϑ_s, and three entries that re-estimate external inputs jointly with them: the
# stationary SD of log ability s_z, the ability-shock dispersion σ_ϵ (the occupational block is
# re-inverted at it: `occupation_block_at`) and the log deviation dlogφ of φ from `EXTERNAL.φ`;
# the last entry, c_T0, is the non-pecuniary cost of teaching without a BA (the size of the
# non-BA teaching track, so that ρ_T0 can carry its pay gap to BA teachers: ACS 0.48).
# The first-stage values hold while an entry is 0 (s_z, σ_ϵ: use the external block's; dlogφ:
# external φ). External inputs of M3: college, sBA, Fc, σd (`--external college=1`, …), ωz and
# sT0, the non-BA teachers' schooling floor (`--external sT0=0.52`; 0: none).
# The 23rd entry (2026-10-07) is the local peer effect ψ_a (`Params.ψa`): school quality
# Q_l = (2H̃_T/M_l)^σ ā_l^{ψ_a} with ā_l the students' mean ability in l. Like the other extensions
# it is 0 (off) unless fixed (`--fix psi_a=0.25`) or freed (`--free psi_a`).
# Three further entries (2026-10-15) are switches, off at 0 and bounded by `EXTENSION_BOUNDS`: the
# 24th, κ_c, is the utility friction of commuting (`Params.κc`; needs `--external split=1`), the 25th,
# θ_ν, is the location-taste scale σν = θ_ν μ (0: the external block's θν; it replaces `--external
# theta_nu` while nonzero, and `--free theta_nu` / `--fix theta_nu=…` address it), and the 26th,
# ζ_BA, is the BA holders' extra suburban amenity (`Params.Bd = (0, ζ_BA)`; needs `--external college=1`).
# The 27th (2026-10-15, closure R), b_T, is teachers' non-pecuniary value of a suburban teaching job
# (`Params.BT = (0, b_T)`, utils; needs `--external split=1`): with δκ on the salary gap it lets the FTE gap be
# matched as well, a reduced form of rationed suburban positions at near-equal pay.
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
# s_z) and the occupational block are the first-pass Table 2 estimates (T2b,
# `spatial_first_stage.jl`): Q_l and t_l fixed, selection into teaching ignored.
# φ matches non-teachers' schooling at the reference (T2, `calibrate_phi`). ωf and γ
# are internal (Table 3, in ϑ); their `EXTERNAL` entries are the ϑ₀ values.
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

"""
The same block on pooled ACS 1-year PUMS 2016–19 (T7), from
data/spatial/estimates/acs_occupations.json (`acs_occupations.py`, which reproduces
every 2009–13 workbook count exactly), in `ACS_OCCUPATIONS` order: counts `men`,
`women` (21 occupations); market-cell 90/10s and wage-sample counts (19 each); the
pooled non-teacher and K–12 teacher 90/10s; mean years of schooling of non-teachers
(home production included) and K–12 teachers; PWGTP-weighted mean log hourly wages of
the market cells and K–12 teachers (men, women). Selected by `occ_sample = :acs2016_19`.
"""
const ACS_2016_19 = (;
    men   = [56691, 25044.5, 39521.5, 24398, 9571, 14264.5, 5102.5, 3930.5, 32045, 52378.5, 45518.5, 21539.5, 41398.5, 18729, 79436, 12207.5, 18117.5, 31781.5, 29972, 145318.5, 13048.5],
    women = [54779, 28117, 12876.5, 29903, 9860.5, 73053, 5197.5, 13269.5, 24888.5, 45694.5, 87920.5, 5351.5, 51805, 4127, 2721.5, 4278.5, 5390.5, 9174, 3363, 214286, 41514.5],
    p9010_men   = [4.256833826, 3.448061015, 3.063285212, 4.180555556, 6.818181818, 5.521475363, 4.335177186, 3.261030092, 4.604036858, 4.862184512, 3.362068966, 3.547387875, 3.604952074, 4.299050709, 3.534045167, 3.170233437, 3.431340082, 3.363575466, 3.453157674],
    p9010_women = [3.684900608, 3.174603175, 3.141830987, 3.357831867, 5.794140813, 4.475510204, 4.1125, 3.6875, 3.794642857, 4.796189913, 2.865932524, 3.585400249, 3.771316735, 4.292307441, 3.917362951, 3.244928491, 3.224288032, 3.245698298, 3.699042408],
    n_men   = [53404, 23404, 36803, 21098, 8445, 12342, 3687, 3145, 29200, 46931, 39453, 19973, 33025, 15218, 69540, 11320, 16203, 26618, 25735],
    n_women = [50147, 25760, 11714, 25351, 8446, 61176, 3673, 10106, 21434, 36797, 74755, 4749, 37258, 3088, 2204, 3668, 4368, 7257, 2461],
    pooled = 3.821426487, teacher_p9010 = (; men = 2.587365743, women = 2.85029759, pooled = 2.786129844),
    schooling = (; nonteach = 13.8662, teach = 15.9337),
    mean_log_men   = [2.86148, 2.95191, 2.99456, 2.61913, 3.15297, 2.67241, 2.47832, 2.42731, 2.98295, 2.56774,
                      2.40102, 2.59559, 2.13482, 2.19399, 2.47649, 2.45677, 2.40764, 2.30031, 2.35687],
    mean_log_women = [2.73119, 2.81502, 2.93148, 2.60164, 3.10615, 2.58143, 2.48930, 2.36819, 2.64903, 2.32501,
                      2.33018, 2.46087, 2.02882, 2.04592, 2.34223, 2.22352, 2.12176, 2.12094, 2.13808],
    k12_mean_log   = (2.56254, 2.49013))

"The occupational sample `occ_sample` (:acs2009_13, the baseline, or :acs2016_19, T7)."
acs_data(occ_sample::Symbol = :acs2009_13) =
    occ_sample === :acs2009_13 ? ACS_OCCUPATIONS :
    occ_sample === :acs2016_19 ? ACS_2016_19 :
    throw(ArgumentError("occ_sample must be :acs2009_13 or :acs2016_19, got :$occ_sample"))

"ACS counts of gender g (1 = men, 2 = women), all 21 occupations."
acs_counts(g; data = ACS_OCCUPATIONS) = g == 1 ? data.men : data.women

"Shares over the 20 non-teaching occupations: the Table 2 share targets."
function occ_shares(g; data = ACS_OCCUPATIONS)
    c = acs_counts(g; data)[1:N_OCC .!= T_OCC]
    return c ./ sum(c)
end

"K–12 teaching share over all 25–34-year-olds of gender g, home production included (Table 3)."
teach_share_data(g; data = ACS_OCCUPATIONS) = acs_counts(g; data)[T_OCC] / sum(acs_counts(g; data))

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
    occupation_block(; σϵ, α, data = ACS_OCCUPATIONS) -> (; A, r)

Table 2 values at ability dispersion ασϵ, by inverting the static Roy model over
the 20 non-teaching occupations with teaching ignored:
- `A`: A_i/A_HP from men's shares (length `N_OCC`, A_HP = 1, A[T] = NaN);
- `r`: relative female wedges r_i = Θ_{i,f}/(ωf A_i) from women's shares, r_HP = 1,
  so that 1 − τω_{i,f} = ωf r_i (`female_wedges`).
At the `EXTERNAL` σϵ these are the first-pass T2b estimates. Selection into teaching
is ignored for both genders until the §4 consistency check (§1, Table 2).
"""
function occupation_block(; σϵ, α, data = ACS_OCCUPATIONS)
    s  = α * σϵ
    lm = invert_roy_shares(occ_shares(1; data), s; ref = HP_OCC)
    lf = invert_roy_shares(occ_shares(2; data), s; ref = HP_OCC)
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

"The share inversion at ασϵ = 1 by occupational sample, `occupation_block_at`'s cache (log Θ/(ασϵ), men and women)."
const OCC_UNIT_INVERSION = Dict{Symbol,NamedTuple{(:men, :women),Tuple{Vector{Float64},Vector{Float64}}}}()

"""
    occupation_block_at(σϵ, α, occ_sample) -> (; A, r)

`occupation_block` at another ability dispersion σϵ (σ_ϵ in ϑ, M3), from the inversion at ασϵ = 1.
Roy shares depend on log Θ only through log Θ/(ασϵ) (`roy_shares`; spatial_first_stage.jl), so the
inversion at any σϵ is ασϵ times the inversion at ασϵ = 1, with log Θ_HP = 0: it is done once per
occupational sample and cached, and a trial σϵ costs nothing. The result agrees with
`occupation_block(; σϵ, α, data)` to the inversion's tolerance (1e-10 in log shares).
"""
function occupation_block_at(σϵ, α, occ_sample::Symbol)
    c = get!(OCC_UNIT_INVERSION, occ_sample) do
        data = acs_data(occ_sample)
        (; men = invert_roy_shares(occ_shares(1; data), 1.0; ref = HP_OCC),
           women = invert_roy_shares(occ_shares(2; data), 1.0; ref = HP_OCC))
    end
    s  = α * σϵ
    lm, lf = s .* c.men, s .* c.women
    return (; A = vcat(exp.(lm), NaN), r = vcat(exp.(lf .- lm), NaN))
end

"""
Mean years of schooling, ages 25–34, ACS 2009–13 (data/Mincer_YrsSchool.xlsx, sheet
`years_of_school`, 2013 column; min 9, max 17 years), averaged over occupation ×
gender cells with the `ACS_OCCUPATIONS` counts: non-teachers (home production
included) 13.6419 years, K–12 teachers 16.0437. As in HHJK, s = years/25, so
s_O = 0.54567 (the φ target) and s_T/s_O = 1.1761 (untargeted). By gender: men
13.4417 and 16.0405, women 13.8402 and 16.0447.
"""
const SCHOOLING = (; nonteach = 13.6419, teach = 16.0437, horizon = 25.0)

"Non-teachers' time share s_O = years/25, the φ target (§1, *Time investment*)."
schooling_target(sch = SCHOOLING) = sch.nonteach / sch.horizon

const EXTERNAL = (;
    σ   = 0.56,     # * class-size curvature: FOO (2013) wage mapping at η = 0.080, K = 12
    θν  = 0.25,     # * relative location-taste scale σν/μ: Eckert–Kleineberg (2024)
    ψ   = 0.0,      # warm-glow curvature: linear warm glow
    η   = 0.080,    # goods elasticity: 2018 spending net of teacher pay (`eta_2018`)
    φ   = 1.0584,   # time elasticity: mean s_O = 13.64/25 at the reference Ξ (`calibrate_phi`, T2)
    μ   = 1.0,      # weight on log consumption: normalization
    γ   = 0.83,     # teacher wage curvature: internal (Table 3, in ϑ); the ϑ₀ value
    α   = 1.0,      # ability elasticity: normalization
    σϵ  = 0.7861,   # T2b: pooled within-occupation 90/10 (3.737), Roy-selected cells, Nz = 5 grid
    ρz  = 0.5766,   # T2b: latent mother–child score correlation (NLSY79/CNLSY)
    s_z = 0.1956,   # T2b: stationary log-z SD, NLSY79 latent wage slope 0.2127 × (1 − η)
    ωf  = 0.9014)   # female outside-option level: internal (Table 3, in ϑ); the ϑ₀ value

"""
Population shares of the frozen baseline (city, suburb): 0.407 and 0.593, the `pop_share` of the exactly
identified fit (exact-polish-base). The non-teachers' position shares `jobs` of the split when `jobs_city`
is not given.
"""
const BASELINE_POP_SHARES = (0.407, 0.593)

"""
The `split` option as a Bool: `true`/`false`, the symbols `Symbol("true")`/`Symbol("false")` (what the
command line passes for `split=true`; the literal `:true` is the Bool), or the number 1 or 0 (`split=1`
arrives as 1.0). Anything else, a misspelling or `split=2`, is an error and not the bundled model.
`parse_flag(v, name)` is the same for any switch of the external block (`move_util`).
"""
parse_split(v) = parse_flag(v, "split")

function parse_flag(v, name)
    v isa Bool && return v
    v === Symbol("true")  && return true
    v === Symbol("false") && return false
    v isa Real && v == 1 && return true
    v isa Real && v == 0 && return false
    throw(ArgumentError("$name must be 1 or 0 (true or false), got $(repr(v))"))
end

"""
    external_params(; kwargs...) -> Params

The fixed block of §1 as a 21-occupation, two-location `Params`: 19 market
occupations, home production (A_HP = 1) and K–12 teaching (`T_OCC`). Keywords
override entries of `EXTERNAL` or any other `Params` field (e.g. `csub`, the
Stone–Geary subsistence in units of C̄). Unless `A` or `τω` is
passed, the occupational block comes from `occupation_block` at the given σϵ and α,
on the occupational sample `occ_sample` (`acs_data`; T7 passes `:acs2016_19`), with
female wedges scaled by the level ωf. The internal fields (κ, B, β, λ, mcost)
are placeholders that `theta_to_params` overwrites, and `Cbar`/`href` are replaced
by the frozen reference levels (`calibration_base`). The location-taste scale
enters as θν, so σν = θν μ, and ability as its stationary SD s_z, so σξ = s_z √(1 − ρz²). The
residence/workplace split is `split` (1 or true, 0 or false; the command line passes a number or a
symbol, and any other value is an error), with the commuting time cost `τk` and, in place of `jobs`,
`jobs_city`, the non-teachers' job share of location 1 (jobs = (jobs_city, 1 − jobs_city)); with the
split on and no `jobs_city` the jobs are the baseline population shares (`BASELINE_POP_SHARES`). The
workplace scale σw is in ϑ (ρw = σw/σν) and may not be passed: `theta_to_params` would overwrite it, so
it is an error. `move_util` (1 or 0, as `split`; 2026-10-15, PREREG amendment 2, candidate U) makes ϑ's r_m a
utility cost of moving (`τmove`, utils) instead of a goods cost (`mcost` = r_m C̄): see `theta_to_params`.
"""
function external_params(; kwargs...)
    x = merge(EXTERNAL, (; occ_sample = :acs2009_13), NamedTuple(kwargs))
    haskey(x, :σw) && throw(ArgumentError("σw is not an external input: σw = ρw σν and ρw is in ϑ; use --fix ρw=… (or --free ρw)"))
    haskey(x, :jobs_city) && haskey(x, :jobs) && throw(ArgumentError("pass jobs_city or jobs, not both"))
    occ = haskey(x, :A) && haskey(x, :τω) ? nothing : occupation_block(; x.σϵ, x.α, data = acs_data(x.occ_sample))
    A  = haskey(x, :A)  ? x.A  : occ.A
    τω = haskey(x, :τω) ? x.τω : female_wedges(occ.r, x.ωf)
    rest = Base.structdiff(x, NamedTuple{(:θν, :s_z, :ωf, :A, :τω, :occ_sample, :jobs_city)})
    split = haskey(x, :split) && parse_split(x.split)
    haskey(x, :split) && (rest = merge(rest, (; split)))
    haskey(x, :move_util) && (rest = merge(rest, (; move_util = parse_flag(x.move_util, "move_util"))))
    haskey(x, :jobs_city) && (rest = merge(rest, (; jobs = [x.jobs_city, 1 - x.jobs_city])))
    split && !haskey(x, :jobs_city) && !haskey(x, :jobs) && (rest = merge(rest, (; jobs = collect(BASELINE_POP_SHARES))))
    return Params(N_OCC, 2; T = T_OCC, A, τω, τmove_off = 0.0, mcost_off = 0.0, rest...,
                  σν = taste_scale(x.μ, x.θν), σξ = innovation_sd(x.s_z, x.ρz), occ_sample = x.occ_sample)
end

"Location-taste scale σν = θν μ for the relative scale θν = σν/μ."
taste_scale(μ, θν = 0.25) = θν * μ

"AR(1) innovation SD for a stationary log-ability SD s_z: σξ = s_z √(1 − ρz²)."
innovation_sd(s_z, ρz) = s_z * sqrt(1 - ρz^2)

"""
Score loading of log Q_l, c = b s_z/χ with b = 1/(1 − η): a school-driven score SD
is worth χ in log earnings, while an ability SD is worth b s_z (§1, *Scale and χ*).
"""
score_loading(; s_z, η, χ = log(1.12)) = s_z / ((1 - η) * χ)

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

No-selection benchmark for the ability dispersion: the pooled within-occupation
90/10 of hourly wages (3.737, non-teachers, ACS 2009–13) read as lognormal with
log-wage SD √(σϵ² + s_z²)/(1 − η) at α = 1, ignoring selection within occupations.
Selection lowers within-occupation dispersion, so this is a lower bound on the T2b
estimate (`spatial_first_stage.jl`), which models the Roy-selected cells.
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

"""
Analytic CFR mapping (§3): Δlog y ≈ β sd / (K(1−η)) for non-teachers without goods costs. With the effectiveness
concept `ρ` of `Measurement.cfr_rho` (round 2026-10-16) the SD of log effectiveness is sd/ρ; ρ = 1 is the mapping above.
"""
cfr_analytic(β, sd; K = 12.0, η = 0.080, ρ = 1.0) = β * (sd / ρ) / (K * (1 - η))

"Invert `cfr_analytic` for β: an initialization, not an estimate (`ρ` as there)."
beta_from_cfr(effect, sd; K = 12.0, η = 0.080, ρ = 1.0) = effect * K * (1 - η) / (sd / ρ)

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
    @printf(io, "  φ  = %.4f (reference Ξ, `calibrate_phi`); at Ξ = 0, s_O = %.4f/%g needs %.4f\n", x.φ,
            SCHOOLING.nonteach, SCHOOLING.horizon, phi_from_schooling(schooling_target(); x.μ, x.η))
    @printf(io, "  s_O, s_T at zero goods cost              = %.4f, %.4f;  s_T/s_O = %.4f   [data %.4f, %.4f; %.3f]\n",
            sh.sO, sh.sT, sh.ratio, schooling_target(), SCHOOLING.teach / SCHOOLING.horizon,
            SCHOOLING.teach / SCHOOLING.nonteach)
    @printf(io, "  ability (T2b): ρz = %.4f, s_z = %.4f, σξ = %.4f;  b s_z = %.4f vs. χ = log(1.12) = %.4f   [0.577, 0.1956]\n",
            x.ρz, x.s_z, innovation_sd(x.s_z, x.ρz), x.s_z / (1 - x.η), χ)
    @printf(io, "     score loading c = b s_z/χ = %.4f   [1.88];  model latent sibling correlation ρz² = %.4f   [0.606 in data]\n",
            score_loading(; x.s_z, x.η, χ), x.ρz^2)
    @printf(io, "  σϵ = %.4f (Roy-selected cells, spatial_first_stage.jl);  no-selection benchmark at s_z: %.4f\n",
            x.σϵ, sigma_eps_placeholder(; η = x.η, x.s_z))
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

"""
    calibrate_phi(ext, θ = THETA0; target = schooling_target(), solver = SOLVER, maxit = 8,
                  tol = 1e-5, verbose = true) -> (; φ, Cbar, href, sO, sol)

φ such that non-teachers' mean time share at the reference solution equals `target`
(§1, *Time investment*: with r_m > 0 the time share depends on Ξ, so φ is set at the
reference Ξ rather than by the Ξ = 0 formula `phi_from_schooling`). Iterates φ jointly
with the reference levels C̄ and h̄, starting from `ext`; each step rescales μφ by the
ratio of the Ξ = 0 inversions b(s) = s(1 − η)/(1 − s) at the target and at the model's
mean s_O. Copy the result into `EXTERNAL.φ`, then refreeze the reference.
"""
function calibrate_phi(ext::Params, θ = THETA0; target = schooling_target(), solver = SOLVER,
                       maxit = 8, tol = 1e-5, verbose = true)
    φ, Cbar, href = ext.φ, ext.Cbar, ext.href
    binv(s) = s * (1 - ext.η) / (1 - s)
    sol, sO = nothing, NaN
    for it in 1:maxit
        p   = theta_to_params(θ, with_params(ext; φ, Cbar, href))
        sol = solve_ge(p; solver..., verbose = false)
        sol.converged || error("φ calibration: GE did not converge at iteration $it")
        sO  = s_summary(sol).meanO
        C1, h1 = mean_consumption(sol), mean_child_h(sol)
        dref = max(abs(C1 / Cbar - 1), abs(h1 / href - 1))
        verbose && @printf("  φ %d: φ = %.6f  s_O = %.6f (target %.6f)  C̄ = %.8f  h̄ = %.8f  rel. change %.2e\n",
                           it, φ, sO, target, C1, h1, dref)
        Cbar, href = C1, h1
        abs(sO - target) < tol && dref < 1e-4 && return (; φ, Cbar, href, sO, sol, converged = true)
        it < maxit && (φ *= binv(target) / binv(sO))
    end
    @printf("  ⚠ φ calibration stopped after %d iterations: s_O = %.6f vs %.6f\n", maxit, sO, target)
    return (; φ, Cbar, href, sO, sol, converged = false)
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
# inverse-variance weights. District SEs are whole-CZ bootstrap SEs (B = 1999, seed
# 20260929; data/spatial/estimates/district_bootstrap.md); the diagonal weighting
# ignores their correlation (pupils–income 0.53). Literature scales are the papers'
# SEs without transport uncertainty; ACS scales are modelling tolerances.
# -----------------------------------------------------------------------------
"""
    Target(key, value, se, param, source, note, active = true)

A calibration target. `key` names a field of `model_moments`; `value` is NaN
while the target is pending, which drops it from the criterion; `se` is its scale
in the criterion; `param` is the principal parameter it informs (§2). An inactive
target is a validation moment: reported against its data value, left out of the
criterion (`--activate key` restores it).
"""
struct Target
    key::Symbol
    value::Float64
    se::Float64
    param::Symbol
    source::String
    note::String
    active::Bool
end
Target(key, value, se, param, source, note) = Target(key, value, se, param, source, note, true)

"Copy of target `t` with the named fields replaced."
with_target(t::Target; kwargs...) =
    Target((get(kwargs, f, getfield(t, f)) for f in fieldnames(Target))...)

"""
Teachers' 90/10 of hourly wages, ACS 2009–13 (wages_occ_shares_v2.xlsx, sheet
`90_10_hr_wages_weighted`): men 2.596925, women 2.829657, and their wage-sample-
weighted average 2.772778 (weights 0.2444/0.7556, `Measurement().teacher_weights`).
"""
const TEACHER_P9010 = (; men = 2.596925, women = 2.829657, pooled = 2.7727776991105824)

"""
Literature targets (T3; sources and arithmetic in data/spatial/estimates/literature_targets.md):
- CFR (2014b), Table 3 col. 2: a 1-SD better teacher (value added) for one year raises
  earnings at 28 by \$285.55 (SE 87.64) on a mean of \$21,256, i.e. 1.34% (0.41 pp);
  col. 1 gives 1.65% (0.43 pp), a Table 5 case. In log points: 0.013344 (SE 0.004069).
- BFM (2007), Table 7 panel B col. 4: \$19.70 (SE 7.40) a month of housing user cost,
  1990, for a one-SD better school-average score (all households), over mean monthly
  consumption of \$3,309 (CE 1990–91, San Francisco area, \$39,707 a year): 0.595%
  (SE 0.224%). US (\$2,365) and West (\$2,705) denominators are Table 5 cases
  (`--wtp-consumption`). The score unit is `Measurement().wtp_score_sd`.
"""
const LITERATURE = (;
    cfr = log(1 + 285.55 / 21_256), cfr_se = (87.64 / 21_256) / (1 + 285.55 / 21_256),
    wtp = 19.70, wtp_se = 7.40, consumption_month = 39_707 / 12)

"SE of the sum of two estimates with SEs `s1`, `s2` and correlation `ρ`: the payroll-per-pupil gap is the salary gap plus the FTE-per-pupil gap."
sum_se(s1, s2, ρ) = sqrt(s1^2 + s2^2 + 2 * ρ * s1 * s2)

"The moments of Table 3 (§2), in the note's order, and the Table 4 validation moments the criterion leaves out."
function default_targets()
    cfr = LITERATURE.cfr
    return [
        Target(:male_teach_share, teach_share_data(1), 0.1 * teach_share_data(1), :logκ̃,
               "ACS 2009–13, men 25–34 (wages_occ_shares_v2.xlsx: moments_shares)",
               "K–12 teachers over all men, home production included; scale = 10% of value"),
        Target(:female_teach_share, teach_share_data(2), 0.1 * teach_share_data(2), :ωf,
               "ACS 2009–13, women 25–34 (wages_occ_shares_v2.xlsx: moments_shares)",
               "K–12 teachers over all women, home production included; scale = 10% of value"),
        Target(:teacher_p9010, TEACHER_P9010.pooled, 0.05 * TEACHER_P9010.pooled, :γ,
               "ACS 2009–13 (wages_occ_shares_v2.xlsx: 90_10_hr_wages_weighted)",
               "men 2.597, women 2.830, wage-sample weights 0.244/0.756; scale = 5% of value"),
        Target(:gap_salary, 0.0086, 0.0092, :δκ, "district_audit.md, variant 28 (spatial_moments.json: 0.01623 on 176 CZs)",
               "CWIFT-deflated teacher salary per K–12 FTE, all 188 CZs (zero-coded IL/NM/AK wage bills filled with 0.881 × instruction salaries, NYC as one district, pre-K teachers out of the FTE); bootstrap SE"),
        Target(:gap_pupils, 0.3764, 0.1154, :ΔB, "spatial_moments.json: gap_pupils_locale",
               "enrollment, 188 CZs; bootstrap SE"),
        Target(:cfr_effect, cfr, LITERATURE.cfr_se, :β, "CFR (2014b) Table 3 col. 2: 1.34% earnings at 28 per teacher-VA SD-year",
               "log(1.0134); scale = CFR's SE (0.41 pp), no transport uncertainty added"),
        Target(:wtp_share, LITERATURE.wtp / LITERATURE.consumption_month, LITERATURE.wtp_se / LITERATURE.consumption_month,
               :λ, "BFM (2007) Table 7 col. 4: \$19.70/month (SE \$7.40) per school-average SD",
               "over \$3,309/month consumption (CE 1990–91, SF area); scale from BFM's SE"),
        Target(:gap_income, 0.1122, 0.0168, :r_m, "district_bootstrap.md: pm_ftyr_all (ACS 2014–18 school districts, B20005)",
               "log median earnings of full-time year-round workers 16+, pooled over each location's districts, 188 CZs; vs. model log median market earnings; bootstrap SE. " *
               "Validation: it barely moves with r_m (Jacobian at the first fit), so r_m is identified by the move rate (§4)", false),
        Target(:gap_teachers_pp, -0.0177, 0.0131, :none, "district_audit.md, FTE variant 11 (spatial_moments.json: −0.01216)",
               "K–12 FTE teachers per K–12 pupil (pre-K teachers and pupils removed), 188 CZs, vs. model headcount per student; bootstrap SE. " *
               "Validation: the one overidentifying moment, left out so that the eight parameters are exactly identified (§4)", false),
        Target(:move_rate, 0.2156, 0.05, :r_m, "NLSY79 (nlsy_extensions.json T3/nlsy79_main): city↔suburb, ages 14–16 → 25–34",
               "SE 0.011; scale 0.05 spans the definitional range (0.168–0.314 across NLSY79/97 variants)"),
        # Validation moments added for the 2026-10 experiments (target tradeoffs, teacher pay, sorting).
        Target(:teacher_premium_m, 0.0557, 0.05, :none, "ACS PUMS 2009–13 (acs_occupations.py), men",
               "mean log hourly wage of K–12 teachers minus that of non-teacher market workers; PUMS cell means weighted by " *
               "occupation counts, as the model weights cells by its shares (+0.042 with wage-sample weights, the §4 check); scale 0.05", false),
        Target(:teacher_premium_f, 0.0860, 0.05, :none, "ACS PUMS 2009–13 (acs_occupations.py), women",
               "as for men (+0.070 with wage-sample weights); scale 0.05", false),
        Target(:s_ratio, SCHOOLING.teach / SCHOOLING.nonteach, 0.05 * SCHOOLING.teach / SCHOOLING.nonteach, :none,
               "ACS 2009–13 (Mincer_YrsSchool.xlsx): 16.04 / 13.64 years",
               "mean time share of teachers over non-teachers; scale 5% of value", false),
        Target(:gender_gap_market, 0.1046, 0.05, :none, "ACS PUMS 2009–13 (acs_occupations.py)",
               "men minus women, mean log hourly wage of non-teacher market workers (occupation-count weights; 0.103 with " *
               "wage-sample weights); scale 0.05", false),
        Target(:city_share_gradient, -0.295, 0.037, :none, "NLSY79 (nlsy_extensions.md, T3): MSA children 14–16 in the central city, by parental-income tercile",
               "top minus bottom tercile, 16.3 − 45.8 pp (NLSY97 −26.0); model: share of the young schooled in location 1 by tercile " *
               "of the parent's pre-tax income; household-bootstrap SE (settings/nlsy_composition.md)", false),
        Target(:teacher_z_gap_market_sd, 0.614, 0.062, :none, "NLSY79 (settings/nlsy_teachers.md): K–12 teachers at 25–34 vs employed non-teachers",
               "gap in mean latent AFQT, latent SDs (household-bootstrap SE); model: mean log z of teachers minus non-teacher market " *
               "workers, SDs of log z (`teacher_selection`); −0.334 against non-teacher college graduates", false),
        Target(:city_share_gradient_z, -0.222, 0.036, :none, "NLSY79 (settings/nlsy_composition.md): MSA children 14–16 in the central city, by own-AFQT tercile",
               "top minus bottom tercile, 21.6 − 43.8 pp; model: the young in location 1 by tercile of own z (`own_z_composition`); " *
               "household-bootstrap SE", false),
        Target(:payroll_pp_gap, -0.0091, sum_se(0.0092, 0.0131, 0.4058), :none, "district_audit.md: the salary (variant 28) and FTE (variant 11) gap targets",
               "Δ log teacher payroll per K–12 pupil, the sum of the two gaps (0.0086 − 0.0177); scale from the whole-CZ bootstrap, " *
               "sqrt(0.0092² + 0.0131² + 2ρ·0.0092·0.0131) with their locale correlation ρ = 0.406 (district_settings_targets.toml); " *
               "model: Δ log of teacher pay per student", false),
        # Validation moments of the residence/workplace split (2026-10; data/spatial/estimates/settings/commuting_all.md,
        # teacher_residence.md, teacher_residence_targets.toml). They need `--external split=1`; in the bundled model
        # the commuting shares, the net commuting gap and the supply elasticity are NaN (zero by construction, not
        # predictions) and are left out of the moment tables.
        Target(:commute_share_all, 0.107, 0.03, :ρw, "ACS 2018 PUMS (commuting_all.md): workers living and working in a city or suburb POWPUMA",
               "share crossing between the city and suburb classes, all workers, 10.70% (sampling SE 0.03 pp; a lower bound, since a POWPUMA is a " *
               "county or more, and 5–22% across class cuts), so the scale is 0.03; vs. model share of market workers with r ≠ w. " *
               "Meant to identify ρw: `--free ρw --activate commute_share_all`", false),
        Target(:commute_share_teach, 0.084, 0.03, :none, "ACS 2018 PUMS (commuting_all.md): K–12 teachers (OCCP 2310-2330)",
               "share crossing between the city and suburb classes, 8.41% (SE 0.15 pp), same lower-bound scale 0.03; vs. model share of teachers with r ≠ w", false),
        Target(:net_commute_gap, 0.075, 0.038, :none, "ACS 2018 (teacher_residence_targets.toml [k12]): public K–12 teacher residents per K–12 FTE",
               "Δ log teachers living in the suburb over the city less Δ log teachers working there, +0.075 (SE 0.038; across specifications " *
               "−0.05 to 0.14): more teachers live in the suburb than work there. The national POWPUMA-class PUMS figure +0.108 (sampling SE " *
               "0.008) is only a cross-check; model: Δ log n_res − Δ log n_work", false),
        Target(:gap_teach_res_emp, 0.298, 0.041, :none, "ACS 2014–18 (teacher_residence_targets.toml [k12]): public K–12 teachers per employed resident",
               "Δ log (public K–12 teachers per employed resident), suburb minus city, +0.298 (SE 0.041; across specifications 0.18 to 0.37); " *
               "model: Δ log (teachers residing / market workers residing)", false),
        Target(:teacher_supply_elasticity, 3.7, 1.0, :ρw, "Ransom–Sims (2010) 3.7; Falch (2010) ≈ 1.4",
               "own-wage elasticity of teacher supply to a district: single-district estimates, so upper bounds for a location whose neighbours' pay " *
               "moves with it; scale 1.0 is a tolerance, not an SE. Model: d log(n_work,2/n_work,1)/d log(κ₂/κ₁), a one-sided household re-solve " *
               "at κ₂ e^{0.01} holding (H̃_T, M, t) fixed (`teacher_supply_elasticity`). Meant to identify ρw in place of the crossing " *
               "share: `--free ρw --activate teacher_supply_elasticity`", false),
        # The joint education and occupation choice (M3; `Params.college`, `--external college=1`), 2026-10-05: the targets of the
        # eight M3 parameters and the validation moments, all inactive until `--activate`d (values: settings/education_targets.toml,
        # D3; settings/nlsy_ba_slopes_targets.toml, D3b; the model counterparts are NaN without college).
        Target(:ba_share_nonteach, 0.3231, 0.0162, :c_C, "ACS 2009–13 (education_targets.toml ba_share_nonteacher)",
               "BA or more among non-teachers 25–34, home production included (T2 convention); scale 5% of the value", false),
        Target(:ba_z_gap_sd, 1.2189, 0.019, :c_Cz, "NLSY79 (education_targets.toml afqt_gap_ba, at the ACS BA share)",
               "mean latent AFQT of 16+ years minus the rest, latent SDs (1.227 at the NLSY79 share 0.258, SE 0.020); model: mean log z of BA " *
               "minus non-BA adults in SDs of log z; AFQT at 15–23 partly reflects schooling (D3 suggests scale 0.05)", false),
        Target(:ba_premium, 0.3610, 0.02, :ρ_c, "ACS 2009–13 (education_targets.toml ba_premium_within_occupation)",
               "BA minus non-BA log wage within occupation × gender and age (raw 0.481: the model has no BA-specific occupational sorting, " *
               "so ρ_c would absorb it); scale 0.02", false),
        Target(:ba_share_teach, 0.8744, 0.01, :c_T0, "ACS 2009–13 (education_targets.toml teacher_ba_share)",
               "BA share of K–12 teachers, codes 2300–2340 (core codes 0.961); the rest is the non-BA teaching track; scale 0.01. " *
               "With c_T0 (the cost of the non-BA track) it is c_T0's moment and ρ_T0 is pinned by `teacher_ba_wage_gap`; without it " *
               "(the first release, 16 parameters) it is ρ_T0's, which then overpays the penalty (the BA–non-BA teacher gap 1.05 vs 0.48)", false),
        Target(:teacher_ba_wage_gap, 0.4819, 0.0094, :ρ_T0, "ACS PUMS 2009–13 (education_targets.toml teacher_ba_wage_gap)",
               "mean log hourly wage of BA+ minus non-BA K–12 teachers (codes 2300–2340), men and women pooled, convention P (men 0.202, " *
               "women 0.566; with age and census-code controls 0.354); the non-BA track is the pre-K and 'other instructors' codes, so the " *
               "raw gap is its counterpart; model: BA minus non-BA teachers' mean log wage; sampling SE 0.0094", false),
        Target(:s_T_mean_ba, 0.6573, 0.01, :ϑ_s, "ACS 2009–13 (education_targets.toml s_T_mean_ba)",
               "mean years of BA+ teachers / 25 = 16.43/25 (43% hold a graduate degree); model: mean time share of BA teachers; scale 0.01", false),
        Target(:s_O_mean, schooling_target(), 0.005, :dlogφ, "ACS 2009–13 (Mincer_YrsSchool.xlsx; education_targets.toml s_O_mean)",
               "mean time share of non-teachers, 13.64/25: the φ target (re-fitted by dlogφ when college is on); scale 0.005 (0.125 years)", false),
        Target(:wage_slope_ba, 0.1608, 0.0083, :s_z, "NLSY79 (nlsy_ba_slopes_targets.toml wage_slope_ba, D3b)",
               "log wage per SD of latent AFQT holding a BA indicator and gender fixed, observed slope/√R (the convention of the 0.2127 wage_slope it " *
               "pairs with; EIV 0.1691, split-half IV 0.1635); bootstrap SE 0.0083 (the file's criterion scale, which adds half the spread of the three " *
               "versions, is 0.0093: `--target-se wage_slope_ba=0.0093`); with the BA control s_z ≈ (1 − η) 0.161 = 0.148", false),
        Target(:pooled_p9010, 3.7374, 0.05 * 3.7374, :σ_ϵ, "ACS 2009–13 (wages_occ_shares_v2.xlsx: 90_10_hr_wages_weighted)",
               "pooled within-occupation 90/10 of non-teachers' hourly wages, wage-sample weights; scale 5% of the value; needs " *
               "`--measure pooled_p9010=1` (set by --activate)", false),
        # validation of the education block
        Target(:teacher_z_gap_ba_sd, -0.3335, 0.0625, :none, "NLSY79 (settings/nlsy_teachers.md): teachers vs non-teacher college graduates",
               "gap in mean latent AFQT, latent SDs (household-bootstrap SE); model: all teachers' mean log z minus that of BA non-teacher " *
               "market workers, SDs of log z", false),
        Target(:teacher_z_gap_within_ba_sd, -0.2384, 0.0617, :none, "NLSY79 (education_targets.toml teacher_z_gap_within_ba)",
               "teachers with 16+ years minus non-teachers with 16+ years; model: BA teachers vs BA non-teacher market workers", false),
        Target(:teacher_premium_ba_m, -0.2610, 0.05, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_premium_ba_m)",
               "mean log hourly wage of all K–12 teachers (codes 2300–2340) minus BA non-teacher workers, men; scale 0.05", false),
        Target(:teacher_premium_ba_f, -0.2081, 0.05, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_premium_ba_f)",
               "as for men, women; scale 0.05", false),
        Target(:teacher_premium_within_ba_m, -0.2374, 0.05, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_premium_within_ba_m)",
               "BA+ teachers minus BA+ non-teacher workers, men; scale 0.05", false),
        Target(:teacher_premium_within_ba_f, -0.1368, 0.05, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_premium_within_ba_f)",
               "BA+ teachers minus BA+ non-teacher workers, women; scale 0.05", false),
        Target(:teacher_ba_wage_gap_m, 0.2017, 0.0219, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_ba_wage_gap_m)",
               "BA+ minus non-BA K–12 teachers, mean log hourly wage, men (age and census-code controlled 0.151); sampling SE", false),
        Target(:teacher_ba_wage_gap_f, 0.5664, 0.0095, :none, "ACS PUMS 2009–13 (education_targets.toml teacher_ba_wage_gap_f)",
               "as for men, women (age and census-code controlled 0.418); sampling SE", false),
        Target(:teacher_p9010_ba, 2.4166, 0.05 * 2.4166, :γ, "ACS PUMS 2009–13, computed in the 2026-10-05 review of M3 (tests/data_teacher_p9010_by_ba.py: the sample of education_occupation.py, the workbook's record-level rule)",
               "90/10 of hourly wages of BA+ K–12 teachers, men 2.4486 and women 2.4062 (10,963 and 33,407 records) pooled with the teacher_p9010 weights 0.2444/0.7556; " *
               "the same computation reproduces the workbook's all-teacher 2.597/2.830 to 0.1% (2.596/2.833). Core codes 2310–2330: 2.347 (men 2.413, women 2.326); " *
               "core and BA+: 2.291; non-BA teachers: 3.461 (men 4.004, women 3.285). The all-teacher target `teacher_p9010` mixes the two tracks; this is " *
               "γ's moment when the non-BA track is not to enter its identification (`--deactivate teacher_p9010 --activate teacher_p9010_ba`); scale 5% of the value (a tolerance; no bootstrap SE)", false),
        Target(:s_T_mean_nonba, 12.9878 / 25, 0.01, :none, "ACS 2009–13 (education_occupation.md §1.1; education_occupation.json years_nba|pooled|teacher|P)",
               "mean years of schooling (capped 9–17) of K–12 teachers without a BA over 25: 12.99 (0.010) years; the schooling floor `sT0` " *
               "of the non-BA teaching track (`--external sT0=0.52`), which would otherwise choose s_T* ≈ 0.30 (7.5 years); scale 0.01", false),
        Target(:ma_premium_implied, 0.1564, 0.0068, :none, "NCES Digest 2022 Table 211.20, NTPS 2017–18 (education_targets.toml ma_ba_salary_premium)",
               "log salary of teachers with a master's over a bachelor's at 1–5 years of experience; model: ϑ_s log(17/16), in the units of the " *
               "s_T_mean_ba target (ACS years capped at 17: a graduate degree is one year beyond the BA; with a master's at 17.5 years, " *
               "`ma_premium_implied_175`, the premium is ϑ_s log(17.5/16))", false),
        Target(:wage_slope, 0.21266, 0.0077, :none, "NLSY79 (nlsy_ability.json): latent wage slope, no BA control",
               "unconditional slope per SD of latent AFQT (the first-stage s_z (1 − η)⁻¹ target); validation once s_z is in ϑ", false),
        Target(:ba_share, 0.3391, 0.017, :none, "ACS 2009–13 (education_targets.toml ba_share_all)",
               "BA or more among all 25–34 (home production included); scale 5%", false),
        Target(:ba_share_z_gap, 0.653, 0.05, :none, "NLSY79 (nlsy_ba_slopes.md §4; education_occupation.md §3, latent terciles at the ACS share)",
               "P(BA | top tercile of latent AFQT) − P(BA | bottom): 0.704 − 0.052; the latent scale, as the model has no score noise; scale 0.05", false),
        Target(:ba_share_score_gap, 0.6243, 0.05, :none, "NLSY79 (nlsy_ba_slopes_targets.toml p_ba_afqt_tercile_*_at_acs_share; nlsy_ba_slopes.md §4)",
               "P(BA | top tercile of the OBSERVED AFQT) − P(BA | bottom) = 0.6925 − 0.0682 (middle 0.2566): all persons, cell-shifted like for like to the ACS " *
               "all-persons BA share 0.3391 (the first release's mixed version, shifted to the non-teacher share 0.323, was 0.670 − 0.062 = 0.608, " *
               "0.237 in the middle: not like for like); model: the same by tercile of a noisy score of reliability 0.916 (`ba_share_by_score`); scale 0.05", false),
        Target(:teacher_slope_ratio, 0.5565, 0.3188, :none, "NLSY79 (nlsy_ba_slopes_targets.toml teacher_slope_ratio, D3b)",
               "wage–AFQT slope of all K–12 teachers (BA control) over that of non-teacher BA workers, observed/√R; the closed form γ(1−η)/(1−γη) " *
               "is 0.600 at the baseline γ = 0.619; a weak check (147 teachers, scale = the bootstrap SE)", false),
    ]
end

"""
    with_wtp_target(targets, consumption_month; wtp = 19.70, se = 7.40)

The BFM target as a consumption share with another denominator: mean monthly
consumption in the same price year as the \$19.70 estimate (1990 dollars). The
default is `LITERATURE.consumption_month`.
"""
function with_wtp_target(targets, consumption_month; wtp = LITERATURE.wtp, se = LITERATURE.wtp_se)
    return [t.key === :wtp_share ?
            with_target(t; value = wtp_share_target(wtp, consumption_month), se = se / consumption_month,
                        note = @sprintf("\$%.2f / \$%.2f monthly consumption; scale from BFM's SE", wtp, consumption_month)) :
            t for t in targets]
end

"Active targets with a finite value and scale: the criterion's moments."
active_targets(targets) = filter(t -> t.active && isfinite(t.value) && isfinite(t.se) && t.se > 0, targets)

"Moments that exist only under the residence/workplace split (NaN in the bundled model): they cannot be targets without `split`."
const SPLIT_ONLY_MOMENTS = (:commute_share_all, :commute_share_teach, :commute_share_nonteach, :net_commute_gap,
                            :teacher_supply_elasticity, :commute_dir_teach, :commute_dir_nonteach,
                            :teacher_own_elasticity, :teacher_cross_elasticity)

# -----------------------------------------------------------------------------
# 4. The internal parameter vector
# -----------------------------------------------------------------------------
"""
Names of ϑ. `logκ̃` is the log teacher pay scale normalized at the reference human
capital h̄, κ̄ = κ̃ h̄^{1−γ} (`theta_to_params`); ϑ vectors saved before 2026-09-29
16:00 hold log κ̄ under the name `logκ` instead. Everything after γ is an extension
parameter (`EXTENSION_BOUNDS`), appended in the order it was added; shorter vectors
saved earlier (eight values before 2026-10-02, eleven before c_Tz, twelve before ρw, twenty-one
before c_T0) are extended by `pad_theta` with the missing extensions off. ρw = σw/σν is the workplace
taste scale of the residence/workplace split (`Params.split`, `--external split=1`) relative to the
location taste: it is read only with `split`, which needs it fixed or freed, and is inert otherwise.
The ten entries added with the joint education and occupation choice (M3, 2026-10-05):
`ρw` (the split's parameter, 13th so that the layout does not shift), `ρ_c`,
`c_C`, `c_Cz`, `ρ_T0`, `ϑ_s` (the education block of `Params.college`), and `s_z`, `σ_ϵ`, `dlogφ`,
which re-estimate the stationary SD of log ability, the ability-shock dispersion and φ
(`theta_to_params`); their off value 0 means the external block's `s_z`, `σϵ`, `φ`; and `c_T0`,
the non-pecuniary cost of teaching without a BA (appended after the first review: it frees
`ρ_T0` to carry the BA–non-BA teacher pay gap while `c_T0` carries the track's size).
`ψ_a` (2026-10-07; ϑ vectors saved before it have twenty-two values) is the elasticity of human capital
in local peer ability, Q_l = (2H̃_T/M_l)^σ ā_l^{ψ_a} (`Params.ψa`; `build_grids`, `peer_log_ability`).
Three switches were added on 2026-10-15 (ϑ vectors saved before them have twenty-three values; all are
0, off, unless fixed or freed): `κ_c`, the utility friction of commuting (`Params.κc`, utils per commuting
pair; it needs `split`), `θ_ν`, the location-taste scale relative to μ, σν = θ_ν μ (0: the external block's
`θν`, so that `--external theta_nu=…` still sets it; the workplace scale is then σw = ρw σν with the
new σν), and `ζ_BA`, the extra suburban amenity of BA holders (`Params.Bd = (0, ζ_BA)`, utils; it needs
`college`). `θ_ν` is also the ϑ-name of the external `θν` (`THETA_SYNONYMS`). `b_T` (closure R, added later on
2026-10-15; ϑ vectors saved before it have twenty-six values) is teachers' non-pecuniary value of teaching in the
suburb, `Params.BT = (0, b_T)` (utils, inside the workplace nest; it needs `split`). `ζ_y` (round 2026-10-15c;
ϑ vectors saved before it have twenty-seven values) is the income-graded suburban taste, `Params.By = (0, ζ_y)`, utils
per log point of the node's predicted pre-tax income relative to the frozen reference `lyref` (`income_taste`).
"""
const THETA_NAMES = (:logκ̃, :δκ, :ΔB, :β, :λ, :r_m, :ωf, :γ, :r_p, :ζ, :c_T, :c_Tz,
                     :ρw, :ρ_c, :c_C, :c_Cz, :ρ_T0, :ϑ_s, :s_z, :σ_ϵ, :dlogφ, :c_T0, :ψ_a,
                     :κ_c, :θ_ν, :ζ_BA, :b_T, :ζ_y)

"""
The names under which an external input also reads as an entry of ϑ (`--fix sigma_eps=…`,
`--free theta_nu`). `--external` keeps addressing the external block (σϵ, θν), and a nonzero ϑ entry overrides it.
"""
const THETA_SYNONYMS = Dict(:σϵ => :σ_ϵ, :θν => :θ_ν)

"Number of parameters of the eight-parameter model, before the extensions."
const N_THETA_BASE = 8

"""
A ϑ saved before some extensions existed, with the missing ones appended at their off values (0):
eight, eleven, twelve, twenty-one, twenty-two, twenty-three and twenty-six values (before r_p, c_Tz, ρw, c_T0,
ψ_a, the switches κ_c, θ_ν, ζ_BA, b_T, and ζ_y) all load, and a full vector is returned as it is.
"""
pad_theta(θ) = N_THETA_BASE ≤ length(θ) ≤ length(THETA_NAMES) ?
               vcat(collect(Float64, θ), zeros(length(THETA_NAMES) - length(θ))) :
               throw(ArgumentError("ϑ has $(length(θ)) values; expected $N_THETA_BASE to $(length(THETA_NAMES))"))

"""
Starting values, re-tuned after T2b at the frozen reference. The first entry is
log κ̃, the pay scale at the reference human capital: κ̄ = κ̃ h̄^{1−γ}
(`theta_to_params`), so κ̃ = 0.2971 is κ̄ = 0.1980 at h̄ = 0.09187 and γ = 0.83. It gives the
1.9% male teaching share and β = 0.2314 the analytic CFR map at the reference
sd(log h_T) = 0.616 (K = 12), jointly with ωf = 0.9014 for the 6.0% female share
(model: 1.94% and 5.98%). The rest are the code's κ₂/κ₁ = 1.2, B = (0, 0.1),
λ = 0.70, r_m = 0.1813 (§2) and γ = 0.83. ωf and γ equal their `EXTERNAL` entries,
so the reference configuration is the same whether read from ϑ₀ or `EXTERNAL`.
"""
const THETA0 = [log(0.2971), log(0.9 / 0.75), 0.1, 0.2314, 0.70, 0.1813, EXTERNAL.ωf, EXTERNAL.γ,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

"""
Search box for ϑ. §4 imposes β ∈ (0, 1), λ ≥ 0 and r_m ≥ 0; the other limits
bound the Sobol screening region. γ ≤ 1 keeps the teacher wage concave in h_T.
Zero is included for r_m as the mechanism-off boundary. The pay scale κ̃ is
normalized at h̄, so its box does not shift with γ: at the T2b block the male
teaching share near ϑ₀ is 1.8% at κ̄ = 0.190 and 3.9% at 0.2365 (κ̃ = 0.285, 0.355),
so the share is steep in κ̃. Widen a limit if estimates pile up against it: the
first Table 5 panel (2026-09-29) hit κ̃ = 0.18 and 0.50 and r_m = 0.5 in its σ, K and
r_m-profile cases, so these limits were widened to κ̃ ∈ [0.10, 1.0], r_m ∈ [0, 1].
"""
const THETA_BOUNDS = [(log(0.10), log(1.00)), # log κ̃ = log κ̄ − (1 − γ) log h̄
                      (-0.5, 0.5),            # δκ
                      (-0.5, 0.5),            # ΔB
                      (0.0, 0.6),             # β
                      (0.0, 1.5),             # λ
                      (0.0, 1.0),             # r_m
                      (0.6, 1.3),             # ωf
                      (0.4, 1.0),             # γ
                      (0.0, 0.0),             # r_p: off; `--free r_p` uses EXTENSION_BOUNDS
                      (0.0, 0.0),             # ζ: off
                      (0.0, 0.0),             # c_T: off
                      (0.0, 0.0),             # c_Tz: off
                      (0.0, 0.0),             # ρw: off (read only with `split`; `--free ρw` uses EXTENSION_BOUNDS)
                      (0.0, 0.0),             # ρ_c: off
                      (0.0, 0.0),             # c_C: off
                      (0.0, 0.0),             # c_Cz: off
                      (0.0, 0.0),             # ρ_T0: off
                      (0.0, 0.0),             # ϑ_s: off
                      (0.0, 0.0),             # s_z: off (the external block's)
                      (0.0, 0.0),             # σ_ϵ: off (the external block's)
                      (0.0, 0.0),             # dlogφ: off (the external φ)
                      (0.0, 0.0),             # c_T0: off
                      (0.0, 0.0),             # ψ_a: off (no peer effect)
                      (0.0, 0.0),             # κ_c: off (no commuting friction; needs `split`)
                      (0.0, 0.0),             # θ_ν: off (the external block's θν)
                      (0.0, 0.0),             # ζ_BA: off (needs `college`)
                      (0.0, 0.0),             # b_T: off (needs `split`)
                      (0.0, 0.0)]             # ζ_y: off (income-graded suburban taste)

"""
Search box of an extension parameter once freed (`--free name`). r_p is the suburb's
residence cost in units of C̄ (1 = mean consumption): at 0.3 the poorest children born in
the suburb can afford neither staying nor moving (0.8% of the population at the
consumption floor, 2026-10-02 smoke test) and trial solves near it take tens of minutes,
while the housing data put r_p near 0.03, so the box stops at 0.15. ζ (the suburb taste),
c_T (the non-pecuniary cost of teaching) and c_Tz (its gradient per SD of log z; negative
selects teachers positively on ability) are in utils (the location-taste scale is σν = 0.25).
ρw = σw/σν lies in (0, 1], RUM-consistent, and is bounded below at 0.01 (σw = 0.0025: the logit arguments
stay at or below about 40); its off value 0 is outside the box, so `problem_setup` starts it at the bound.
With a 3% commuting cost and a ≈ 10% crossing share the data point to a small σw, a few hundredths of σν
(ρw ≈ 0.066 at τk = 0.03); the τk = 0.015 sensitivity ended at the former bound 0.05 with the crossing share
unmatched.
The M3 entries (`Params.college` on): ρ_c, the log credential premium of a BA, may be negative (a BA
also raises h through the time share, so the premium within occupation is ρ_c plus that return plus
selection); c_C and c_Cz (per SD of log z) are the college cost in utils; ρ_T0 the log pay
penalty of a non-BA teacher (ACS: 0.35–0.48 raw BA gap among teachers); ϑ_s the salary-schedule
elasticity (NCES master's premium 0.156 = ϑ_s log(17.5/16) at ϑ_s = 1.74); s_z and σ_ϵ the
ability scale and shock dispersion (first-stage values 0.1956, 0.7861); dlogφ the log deviation
of φ from its external value; c_T0 the extra utility cost of teaching without a BA (it may be
negative: a net amenity of the track). ψ_a, the elasticity of h in the students' mean ability, lies in
[0, 1]: 1 makes the peers' ability count as much as one's own (α = 1).
The switches of 2026-10-15: κ_c, the utility friction of commuting, is in utils and lies in [0, 3]. It acts
inside the workplace nest, so its scale is σw = ρw σν: 3 utils is 12 σw at ρw = 1 and σν = 0.25, which
leaves a crossing pair all but empty, but at the C4′ value ρw = 0.068 (σw ≈ 0.017) κ_c = 0.5 is already
about 29 σw and shuts commuting. The intended use (candidate S1) therefore fixes ρw high (ρw = 1: the
crossing share 0.107 is reached at κ_c ≈ 0.48 at the C4′ ϑ) and frees κ_c in its place. θ_ν = σν/μ lies in
[0.03, 0.60] (0.25 is the Eckert–Kleineberg value of the external block; a small θ_ν makes location choice
nearly deterministic, so trial solves near the lower end are stiff). σw = ρw θ_ν μ falls with both, so ρw's
lower bound no longer bounds σw once θ_ν is free: the corner ρw = 0.01, θ_ν = 0.03 gives σw ≈ 3e−4; do not
free both near it. ζ_BA, the BA holders' extra suburban amenity, is in utils, in [−1.5, 1.5] (the
ability-graded ζ is bounded by ±1 per SD of log z); it also moves the college choice, so it is read jointly
with c_C (`education_amenity`). The off value of θ_ν, 0, lies outside its box, so `problem_setup` starts a
freed θ_ν at the external block's θν (`external_start!`), as for s_z and σ_ϵ; the other two start at 0.
b_T (closure R), teachers' non-pecuniary value of a suburban teaching job, is in utils, in [−1.5, 1.5]. Like κ_c it
acts inside the workplace nest (scale σw) next to μ log C, so a shift of b_T utils moves teachers between
workplaces about as a suburban pay change of b_T C/(μ y) log points does: roughly 0.85 b_T/μ at S1 and C4′ (the
teacher mean of μy/C is 1.17; the gap is goods spending: education, the moving cost r_m C̄ and the residence cost
r_p C̄); it starts at 0. ζ_y (round 2026-10-15c), the income-graded suburban taste, is in utils per log point of
predicted income and lies in [−1, 1]: matching the parental-income gradient needs a tilt of about 0.2 utils between the
bottom and top income terciles at σν = 0.25, about 1.3 log points apart, so ζ_y ≈ 0.15–0.3; it starts at 0.
"""
const EXTENSION_BOUNDS = Dict(:r_p => (0.0, 0.15), :ζ => (-1.0, 1.0), :c_T => (-1.0, 2.0), :c_Tz => (-1.0, 1.0),
                              :ρw => (0.01, 1.0),
                              :ρ_c => (-0.6, 0.8), :c_C => (-2.0, 3.0), :c_Cz => (-2.0, 1.0), :ρ_T0 => (0.0, 1.0),
                              :ϑ_s => (0.0, 3.0), :s_z => (0.05, 0.40), :σ_ϵ => (0.40, 1.40),
                              :dlogφ => (-0.5, 0.5), :c_T0 => (-1.0, 3.0), :ψ_a => (0.0, 1.0),
                              :κ_c => (0.0, 3.0), :θ_ν => (0.03, 0.60), :ζ_BA => (-1.5, 1.5), :b_T => (-1.5, 1.5),
                              :ζ_y => (-1.0, 1.0))

"Index of a parameter in ϑ, by name (`THETA_SYNONYMS`: `σϵ` reads as `σ_ϵ`)."
function theta_index(name::Symbol)
    name = get(THETA_SYNONYMS, name, name)
    i = findfirst(==(name), THETA_NAMES)
    i === nothing && throw(ArgumentError("unknown parameter $name; use one of $THETA_NAMES"))
    return i
end

"""
    theta_to_params(θ, base) -> Params

κ = (κ̄, κ̄ e^{δκ}) with κ̄ = κ̃ h̄^{1−γ} and h̄ = `base.href` frozen, so that the city
wage of a teacher with human capital h is κ̃ h̄ (h/h̄)^γ and κ̃ keeps its meaning as γ
moves; B = (0, ΔB), β, λ, the symmetric goods cost r_m C̄ with C̄ =
`base.Cbar` frozen, γ, and the female wedges 1 − τω_{i,f} = ωf r_i for i ≠ T. The
relative wedges r_i = (1 − τω_{i,f})/(1 − τω_{HP,f}) are read from `base` (r_HP = 1),
so only their level ωf moves. There is no utility moving cost, unless `base.move_util` (2026-10-15, candidate U of
PREREG amendment 2): then r_m is a utility cost of moving, τmove = r_m utils off the diagonal, and the goods cost
is zero (mcost = 0); the move rate identifies r_m in either reading. The extensions give the
suburb a residence cost p = (0, r_p C̄) and a taste ζ per SD of log z, Bz = (0, ζ/s_z)
with s_z = σξ/√(1 − ρz²), and teaching a non-pecuniary cost c_T + c_Tz log z/s_z; an
eight-value ϑ leaves them off. The workplace taste scale is σw = ρw σν (read only if `base.split`; with the
split off ρw is inert and `σw` stays at `base`'s). The M3 entries (`THETA_NAMES`) set the education block of `Params`
(ρc = ρ_c, cC = c_C, cCz = c_Cz/s_z per unit of log z, ρT0 = ρ_T0, ϑs = ϑ_s, cT0 = c_T0; they act only if
`base.college`) and re-estimate three external inputs when nonzero: s_z (σξ = s_z√(1 − ρz²), and
the per-SD conversions above use it), σ_ϵ (the occupational block, A_i and the relative female
wedges, is re-inverted at it: `occupation_block_at`, cached) and φ = `base.φ` e^{dlogφ}. At 0 they
leave the external block's values. ψa = ψ_a is the peer effect (0: off). The last three entries are the switches of
2026-10-15: κc = κ_c, the utility friction of commuting (read only with `split`), σν = θ_ν μ when θ_ν is nonzero (0:
`base.σν`, the external block's θν; σw = ρw σν then uses this σν), and the BA holders' extra suburban amenity,
Bd = (0, ζ_BA) (read only with `college`). At zero all three leave `base` as it is. The 27th, b_T (closure R), is
teachers' non-pecuniary value of teaching in the suburb, BT = (0, b_T) (read only with `split`; 0 leaves `base`).
The 28th, ζ_y (round 2026-10-15c), is the income-graded suburban taste, By = (0, ζ_y) per log point of predicted income
relative to `base.lyref` (0 leaves `base`).
"""
function theta_to_params(θ, base::Params)
    logκ̃, δκ, ΔB, β, λ, r_m, ωf, γ, r_p, ζ, c_T, c_Tz = pad_theta(θ)
    ρw, ρ_c, c_C, c_Cz, ρ_T0, ϑ_s, s_zθ, σϵθ, dlogφ, c_T0, ψ_a, κ_c, θ_ν, ζ_BA, b_T, ζ_y = pad_theta(θ)[13:28]
    if !iszero(σϵθ)       # M3: the occupational block at σ_ϵ (shares unchanged; cached inversion)
        occ  = occupation_block_at(σϵθ, base.α, base.occ_sample)
        base = with_params(base; σϵ = σϵθ, A = occ.A, τω = female_wedges(occ.r, 1 - base.τω[HP_OCC, 2]))
    end
    iszero(s_zθ) || (base = with_params(base; σξ = innovation_sd(s_zθ, base.ρz)))   # M3: internal s_z
    κ̄  = exp(logκ̃) * base.href^(1 - γ)
    τω = copy(base.τω)
    for i in eachindex(base.A)
        i == base.T && continue
        τω[i, 2] = 1 - ωf * (1 - base.τω[i, 2]) / (1 - base.τω[HP_OCC, 2])
    end
    s_z = base.σξ / sqrt(1 - base.ρz^2)
    φ   = iszero(dlogφ) ? base.φ : base.φ * exp(dlogφ)            # M3: φ re-fitted to non-teachers' schooling
    σν  = iszero(θ_ν) ? base.σν : taste_scale(base.μ, θ_ν)        # Θ: internal location-taste scale
    τmove, mcost = base.move_util ? (offdiag(r_m, 2), zeros(2, 2)) :      # U: r_m in utils
                                    (zeros(2, 2), offdiag(r_m * base.Cbar, 2))
    return with_params(base; κ = [κ̄, κ̄ * exp(δκ)], B = [0.0, ΔB], β, λ, γ, τω, φ, σν, τmove, mcost,
                       pres = [0.0, r_p * base.Cbar], Bz = [0.0, ζ / s_z], cT = c_T, cTz = c_Tz / s_z,
                       σw = base.split ? ρw * σν : base.σw,
                       ρc = ρ_c, cC = c_C, cCz = c_Cz / s_z, ρT0 = ρ_T0, ϑs = ϑ_s, cT0 = c_T0, ψa = ψ_a,
                       κc = κ_c, Bd = [0.0, ζ_BA], BT = iszero(b_T) ? base.BT : [0.0, b_T],
                       By = iszero(ζ_y) ? base.By : [0.0, ζ_y])
end

"ϑ implied by a two-location `Params` (inverse of `theta_to_params`; the entries that re-estimate external inputs, s_z, σ_ϵ, dlogφ and θ_ν, read as off)."
params_to_theta(p::Params) =
    [log(p.κ[1]) - (1 - p.γ) * log(p.href), log(p.κ[2] / p.κ[1]), p.B[2] - p.B[1], p.β, p.λ,
     p.move_util ? p.τmove[1, 2] : p.mcost[1, 2] / p.Cbar,
     1 - p.τω[HP_OCC, 2], p.γ, (p.pres[2] - p.pres[1]) / p.Cbar, (p.Bz[2] - p.Bz[1]) * p.σξ / sqrt(1 - p.ρz^2), p.cT,
     p.cTz * p.σξ / sqrt(1 - p.ρz^2),
     p.split ? p.σw / p.σν : 0.0, p.ρc, p.cC, p.cCz * p.σξ / sqrt(1 - p.ρz^2), p.ρT0, p.ϑs, 0.0, 0.0, 0.0, p.cT0, p.ψa,
     p.κc, 0.0, p.Bd[2] - p.Bd[1], p.BT[2] - p.BT[1], p.By[2] - p.By[1]]

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

"Largest population share allowed at the smooth_log floor (the default `floor_tol`; `--floor-tol`)."
const FLOOR_MASS_TOL = 1e-8

"""
    SpatialModel(base, keys; meas = calibration_measurement(), solver = SOLVER, verbose = false,
                 floor_tol = FLOOR_MASS_TOL)

Callable ϑ → moments: solves the GE at `theta_to_params(ϑ, base)` and returns the
`model_moments` fields named in `keys`, in order. Every call solves from the same
cold start, so evaluations are deterministic. Throws `CalibrationFailure` at an
infeasible or unconverged point, where `floor_tol` is the largest population share
allowed at the consumption floor (`solve_theta`). With `verbose`, logs one line per call.
"""
struct SpatialModel
    base::Params
    meas::Measurement
    solver::NamedTuple
    keys::Vector{Symbol}
    verbose::Bool
    floor_tol::Float64
end

SpatialModel(base::Params, keys; meas = calibration_measurement(), solver = SOLVER, verbose = false,
             floor_tol = FLOOR_MASS_TOL) =
    SpatialModel(base, meas, solver, collect(Symbol, keys), verbose, floor_tol)

"""
    calibration_measurement(; kwargs...) -> Measurement

The calibration's measurement choices: `Measurement` defaults with home production
(`HP_OCC`) split out of the wage and earnings statistics, and the income gap measured
as the gap in log median earnings of market workers, the model object closest to a
district median of earnings (T3). Keywords override fields.
"""
calibration_measurement(; kwargs...) = Measurement(; hp = HP_OCC, income = :log_median_market, kwargs...)

"""
    solve_theta(model, θ) -> (sol, moments)

Solve at ϑ and check feasibility, or throw `CalibrationFailure`. Unaffordable
moves chosen with probability ≈ 0 are allowed: smooth_log then behaves like log,
and they become common when the frozen goods cost is large relative to income. The
population share at the consumption floor, `floor_mass + infeasible_mass` (it chooses a
floored location although another is feasible, or has no feasible one), may not exceed
`model.floor_tol`.
"""
function solve_theta(model::SpatialModel, θ)
    p   = theta_to_params(θ, model.base)
    sol = solve_ge(p; model.solver..., verbose = false)
    sol.converged || throw(CalibrationFailure("GE did not converge in $(model.solver.maxit) iterations"))
    all(0 .< sol.t .< TAX_CAP - 1e-9) || throw(CalibrationFailure("tax rate on its clamp, t = $(fmtvec(sol.t))"))
    sol.s_resid ≤ 1e-8 || throw(CalibrationFailure(@sprintf("(e, s) residual %.1e", sol.s_resid)))
    m = model_moments(sol, model.meas; hh_tol = get(model.solver, :hh_tol, 1e-7), hh_maxit = get(model.solver, :hh_maxit, 1000))
    m.floor_mass + m.infeasible_mass ≤ model.floor_tol ||
        throw(CalibrationFailure(@sprintf("population share %.1e chooses consumption at the smooth_log floor (%.1e with a feasible location, %.1e with none; --floor-tol %.1e)",
                                          m.floor_mass + m.infeasible_mass, m.floor_mass, m.infeasible_mass, model.floor_tol)))
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
          [(t.key, t.value, t.se) for t in targets], model.floor_tol)
    for f in (joinpath("..", "spatial_continuous.jl"), "spatial_moments.jl", "spatial_calibrate.jl")
        write(io, read(joinpath(@__DIR__, f)))
    end
    return bytes2hex(sha256(take!(io)))[1:12]
end

# -----------------------------------------------------------------------------
# 6. Serial diagnostics
# -----------------------------------------------------------------------------
"Standardized residuals (model − target)/se over the active targets."
residuals(m, targets) = [(m[t.key] - t.value) / t.se for t in active_targets(targets)]

"""
    std_residuals(model, θ, act) -> Vector

Standardized residuals (m − target)/se at ϑ over the targets `act`, or `Inf`s where
the solve fails with a `CalibrationFailure` or a `DomainError` (the failures TikTak
records as failed points). A failed point logs the message and ϑ (in `--theta` form, to
re-evaluate it): the search otherwise sees only the `Inf`s. Defined here so that it exists on every worker.
"""
function std_residuals(model::SpatialModel, θ, act)
    try
        _, m = solve_theta(model, θ)
        return [(Float64(m[t.key]) - t.value) / t.se for t in act]
    catch e
        e isa Union{CalibrationFailure,DomainError} || rethrow()
        @printf("  ⚠ point fails: %s;  ϑ = %s\n", sprint(showerror, e), join((@sprintf("%.10g", x) for x in θ), ","))
        flush(stdout)
        return fill(Inf, length(act))
    end
end

"Σ((model − target)/se)²: TikTak's `MomentObjective` with `scales = se`."
criterion(m, targets) = sum(abs2, residuals(m, targets))

"""
    target_whitener(act, corr) -> Matrix

L⁻¹ for the correlation matrix R = LL′ of the standardized residuals of the active
targets `act`, so that Σ(L⁻¹r)² = r′R⁻¹r is the criterion with the full inverse
covariance of the targets in place of the diagonal weights (the scales are SEs, so
Σ = D R D with D = diag(se)). `corr` is `nothing` (identity: the diagonal weights) or
`(; keys, matrix)` from `load_target_corr`; targets it does not cover are uncorrelated
with the rest, as are targets from different data sources.
"""
function target_whitener(act, corr)
    n = length(act)
    R = Matrix{Float64}(I, n, n)
    if corr !== nothing
        idx = Dict(k => i for (i, k) in enumerate(corr.keys))
        for a in 1:n, b in 1:n
            i, j = get(idx, act[a].key, 0), get(idx, act[b].key, 0)
            (a != b && i > 0 && j > 0) && (R[a, b] = corr.matrix[i, j])
        end
    end
    return Matrix(inv(cholesky(Symmetric(R)).L))
end

"""
    load_target_corr(path, scheme) -> (; keys, matrix)

Bootstrap correlations of a partition scheme's district targets from a
`district_settings_targets.toml` file (`data/spatial/district_settings.py`): its table
`scheme` holds `corr_keys` (model keys) and `corr` (rows).
"""
function load_target_corr(path, scheme)
    d = TOML.parsefile(path)[scheme]
    return (; keys = Symbol.(d["corr_keys"]), matrix = reduce(vcat, permutedims.(Vector{Float64}.(d["corr"]))))
end

"""
Print targets, model moments and standardized residuals; pending targets show the model value only. A
validation moment the specification does not define (NaN: the split's moments in the bundled model, which are
zero by construction and not predictions) is left out; an active target that is NaN is shown.
"""
function moment_table(io::IO, m, targets)
    @printf(io, "  %-18s %11s %11s %9s  %-6s %s\n", "moment", "target", "model", "(m−t)/se", "param", "source")
    for t in targets
        !t.active && !isfinite(m[t.key]) && continue      # a validation moment the specification does not define (the split's without split, M3's without college)
        pending = !(isfinite(t.value) && isfinite(t.se))
        @printf(io, "  %-18s %11s %11.5f %9s  %-6s %s\n", t.key,
                pending ? "pending" : @sprintf("%.5f", t.value), m[t.key],
                pending ? "—" : @sprintf("%+.2f%s", (m[t.key] - t.value) / t.se, t.active ? "" : "*"), t.param, t.source)
    end
    any(t -> !t.active, targets) && println(io, "  * validation: not in the criterion")
    @printf(io, "  criterion over %d active targets = %.4f\n", length(active_targets(targets)), criterion(m, targets))
end
moment_table(m, targets) = moment_table(stdout, m, targets)

"Validation outcomes and diagnostics from `model_moments` (everything not targeted)."
function diagnostics_table(io::IO, m)
    @printf(io, "  validation: SEDA growth gap %.4f SD/grade [data 0.0031 (0.0030)];  Δlog Q = %.4f\n",
            m.seda_growth_gap, m.gap_logQ)
    iszero(m.psi_a) || @printf(io, "  peer effect ψa = %.3f: Δlog Q = teachers %.4f + peers %.4f (Δlog ā = %+.4f, log ā %s);  fixed-point residual %.1e\n",
                               m.psi_a, m.gap_logQ_teach, m.gap_logQ_peer, m.log_abar[2] - m.log_abar[1], fmtvec(m.log_abar), m.peer_resid)
    @printf(io, "              score-level gap %.4f SD [0.2887];  child-poverty gap %+.4f [−0.0888] (poor shares %s;\n",
            m.score_gap, m.poverty_gap, fmtvec(m.poverty_share_l))
    println(io, "              provisional: relative line, home-production output counted as income)")
    @printf(io, "  income gaps: mean log %.4f, log mean %.4f, log median %.4f;  market workers: mean log %.4f, log median %.4f\n",
            m.gap_income_mean_log, m.gap_income_log_mean, m.gap_income_log_median,
            m.gap_income_mean_log_market, m.gap_income_log_median_market)
    @printf(io, "  occupations: teaching share %.4f (male %.4f, female %.4f);  teacher 90/10 men %.3f, women %.3f [2.597, 2.830]\n",
            m.teach_share, m.male_teach_share, m.female_teach_share, m.teacher_p9010_g...)
    @printf(io, "  ability (§4 check): wage slope %.4f [NLSY 0.2127];  mother–child ρ %.4f [0.5766];  c = %.4f;  s_T/s_O = %.4f [1.176]\n",
            m.wage_slope, m.rho_pc, m.score_loading, m.s_ratio)
    @printf(io, "  mobility: gross move rate %.4f;  P(1→2) = %.4f, P(2→1) = %.4f [NLSY79 0.361, 0.153];  by parental-income tercile %s [NLSY79 0.221, 0.226, 0.207]\n",
            m.move_rate, m.p12, m.p21, fmtvec(m.move_rate_by_parent_income))
    hasproperty(m, :move_rate_by_own_z) &&
        @printf(io, "            by own-z tercile: move rate %s [NLSY79 own AFQT 0.195, 0.184, 0.273], 1→2 %s [0.284, 0.378, 0.508], 2→1 %s [0.126, 0.118, 0.208]\n",
                fmtvec(m.move_rate_by_own_z), fmtvec(m.c2s_by_own_z), fmtvec(m.s2c_by_own_z))
    @printf(io, "  by parental-income tercile: children in location 1 %s [NLSY79 0.458, 0.360, 0.163; top − bottom %+.3f vs −0.295];\n",
            fmtvec(m.city_share_by_parent_income), m.city_share_gradient)
    @printf(io, "              1→2 %s [0.284, 0.429, 0.407];  2→1 %s [0.168, 0.112, 0.168];  adults in 2 by own-income tercile %s\n",
            fmtvec(m.c2s_by_parent_income), fmtvec(m.s2c_by_parent_income), fmtvec(m.suburb_share_by_own_income))
    @printf(io, "  composition: mean log z gap young %+.3f SD, adults %+.3f SD;  score gap without Q %.4f;  market mean-log wage gap %+.4f, Q-adjusted %+.4f;  mover − stayer by birth location %s\n",
            m.dz_young_sd, m.dz_adult_sd, m.score_gap_noQ, m.gap_meanlog_mkt, m.gap_meanlog_mkt_Qadj, fmtvec(m.mover_premium))
    @printf(io, "  by own-z tercile: children in location 1 %s [NLSY79 own AFQT 0.438, 0.251, 0.216; top − bottom %+.3f vs −0.222]\n",
            fmtvec(m.city_share_by_own_z), m.city_share_gradient_z)
    @printf(io, "  teacher log-wage premium men %+.3f, women %+.3f [+0.056, +0.086];  market gender gap %.3f [0.105];  teachers' mean log z vs non-teachers %+.3f SD (vs market workers %+.3f);  Δ mean log h_T by work location %+.4f\n",
            m.teacher_premium_m, m.teacher_premium_f, m.gender_gap_market, m.teacher_z_gap_sd, m.teacher_z_gap_market_sd,
            m.gap_logh_teacher)
    @printf(io, "  CFR pieces: sd(log h_T) = %.4f;  non-teachers %.5f;  analytic %.5f",
            m.teacher_sd, m.cfr_nonteach, m.cfr_analytic)
    (hasproperty(m, :cfr_rho) && m.cfr_rho != 1.0) &&      # round 2026-10-16: effectiveness concept ρ (Measurement.cfr_rho)
        @printf(io, ";  ρ = %.3f, so the intervention's SD of log effectiveness is sd/ρ = %.4f", m.cfr_rho, m.cfr_sd_eff)
    println(io)
    @printf(io, "  WTP share by location %s;  students/teacher %s (model only);  t = %s;  M = %s\n",
            fmtvec(m.wtp_share_l), fmtvec(m.students_per_teacher), fmtvec(m.t), fmtvec(m.M))
    @printf(io, "  levels (model only): mean C = %.6f, mean h = %.6f;  stationarity dev %.1e\n",
            m.mean_C, m.mean_h, m.stationarity_dev)
    # 2026-10-18: the fertility shock (Params.pkid < 1): pupils (`M`, `pop_share`) and adults by location come apart
    (hasproperty(m, :pkid) && m.pkid < 1) &&
        @printf(io, "  fertility shock: P(child) = %.4f;  adults by location %s (log gap %.4f) vs pupils %s (gap_pupils %.4f);  parents' share of adults %s (log gap %.4f);  parents' WTP share %.5f [BFM households with children ≈ 0.0074]\n",
                m.pkid, fmtvec(m.adult_share), m.gap_adults, fmtvec(m.pop_share), m.gap_pupils,
                fmtvec(m.parent_share_l), m.gap_parent_share, m.wtp_share_parents)
    (hasproperty(m, :pkid_z) && m.pkid < 1 && m.pkid_z != 0) &&
        @printf(io, "  ability-graded fertility: logit slope %.4f per ability SD;  P(child | z) %s;  young's mean log z %+.4f ergodic SDs\n",
                m.pkid_z, fmtvec(m.pk_z), m.fert_tilt_sd)
    @printf(io, "  model levels: population shares %s [suburb 0.579 of partitioned pupils];  t %s, salaries over market income %.4f %s [≈ 0.016; the §1 0.023 includes benefits];  Q %s\n",
            fmtvec(m.pop_share), fmtvec(m.t), m.t_market, fmtvec(m.t_market_l), fmtvec(m.Q))
    @printf(io, "              home production men %.3f, women %.3f;  mean s: teachers %.4f, non-teachers %.4f;  WTP income elasticity %.3f [BFM 0.379], within z %.3f [≈ 1 under log utility] (slope %.4f), share by parental-income tercile %s\n",
            m.hp_share_m, m.hp_share_f, m.s_T_mean, m.s_O_mean, m.wtp_income_elasticity,
            m.wtp_income_elasticity_z, m.wtp_income_slope_z, fmtvec(m.wtp_share_by_parent_income))
    @printf(io, "              value of the suburb's school quality ΔΛ: mean %.2e, by z %s; due to Q alone: mean %.2e, by z %s\n",
            m.dLambda_mean, fmtvec(m.dLambda_z), m.dLambda_quality_mean, fmtvec(m.dLambda_quality))
    if isfinite(m.commute_share_all)         # the residence/workplace split
        @printf(io, "  residence/workplace: crossing share all %.4f [0.107], teachers %.4f [0.084], non-teachers %.4f;  net commuting gap %+.4f [+0.075 (0.038); national PUMS cross-check +0.108];  teachers per employed resident gap %+.4f [+0.298 (0.041)];  payroll per pupil gap %+.4f [−0.0091]\n",
                m.commute_share_all, m.commute_share_teach, m.commute_share_nonteach, m.net_commute_gap, m.gap_teach_res_emp,
                m.payroll_pp_gap)
        @printf(io, "              teachers by residence %s, by workplace %s;  teacher supply elasticity %.2f [Ransom–Sims 3.7, Falch ≈ 1.4: single-district, so upper bounds for a location];  inclusive-value gap I_r − W_rr (not the value of the workplace choice): teachers %+.4f, non-teachers %+.4f, difference %+.4f utils\n",
                fmtvec(m.n_teach_res), fmtvec(m.n_teach_work), m.teacher_supply_elasticity, m.option_value_teach,
                m.option_value_nonteach, m.option_value_gap)
    else
        @printf(io, "  residence/workplace: bundled model (no commuting, supply elasticity or inclusive-value gap);  teachers per employed resident gap %+.4f [+0.298 (0.041)];  payroll per pupil gap %+.4f [−0.0091];  teachers %s\n",
                m.gap_teach_res_emp, m.payroll_pp_gap, fmtvec(m.n_teach_work))
    end
    if hasproperty(m, :teacher_district_elasticity)      # the switches of 2026-10-15 (κ_c, θ_ν, ζ_BA) and their diagnostics
        isfinite(m.commute_dir_teach[1]) &&
            @printf(io, "  commuting direction (live in 1 → work in 2, live in 2 → work in 1): teachers %s [ACS, city–suburb block: 0.0862, 0.0832], non-teachers %s [0.0919, 0.1164]\n",
                    fmtvec(m.commute_dir_teach), fmtvec(m.commute_dir_nonteach))
        @printf(io, "  teacher supply (PREREG A6): small-district own-pay elasticity by workplace %s [A6: ≤ 10];  suburb's own-location elasticity ε22 %.3f [A6: ≤ 3.7, Ransom–Sims], cross ε12 %.3f (household partial equilibrium; ratio ε22 − ε12 above)\n",
                fmtvec(m.teacher_district_elasticity), m.teacher_own_elasticity, m.teacher_cross_elasticity)
    end
    if isfinite(m.ba_share)
        hasproperty(m, :ba_share_by_residence) &&
            @printf(io, "  BA share of adults by residence %s, suburb − city %+.4f (the counterpart of the sealed hold-out B1)\n",
                    fmtvec(m.ba_share_by_residence), m.ba_share_gap)
        @printf(io, "  education (M3): BA share %.4f [0.339] (non-teachers %.4f [0.323], teachers %.4f [0.874]), by own-z tercile %s [NLSY79 latent 0.052, 0.262, 0.704], by tercile of a noisy score %s [observed 0.068, 0.257, 0.692, like for like at 0.339];  BA − non-BA mean log z %+.3f SD [1.219]\n",
                m.ba_share, m.ba_share_nonteach, m.ba_share_teach, fmtvec(m.ba_share_by_z), fmtvec(m.ba_share_by_score), m.ba_z_gap_sd)
        @printf(io, "              BA premium %+.3f (unadjusted %+.3f) [0.361 within occupation, 0.481 raw];  teachers vs BA non-teachers: men %+.3f, women %+.3f [−0.261, −0.208], BA teachers %+.3f, %+.3f [−0.237, −0.137];  teachers' log z vs BA non-teachers %+.3f SD, BA teachers %+.3f [−0.334, −0.238]\n",
                m.ba_premium, m.ba_premium_raw, m.teacher_premium_ba_m, m.teacher_premium_ba_f,
                m.teacher_premium_within_ba_m, m.teacher_premium_within_ba_f, m.teacher_z_gap_ba_sd, m.teacher_z_gap_within_ba_sd)
        @printf(io, "              mean s: BA teachers %.4f [0.657], non-BA teachers %.4f [0.520], BA non-teachers %.4f [0.652], non-BA %.4f [0.497];  master's premium implied %.3f (%.3f at 17.5 years) [0.156];  wage slope %.4f, holding BA fixed %.4f [0.213, 0.161 convention, 0.169 EIV];  teacher/BA non-teacher slope ratio %.3f [0.556 (0.319), closed form %.3f]%s\n",
                m.s_T_mean_ba, m.s_T_mean_nonba, m.s_O_mean_ba, m.s_O_mean_nonba, m.ma_premium_implied, m.ma_premium_implied_175, m.wage_slope, m.wage_slope_ba,
                m.teacher_slope_ratio, m.slope_ratio_closed_form,
                isfinite(m.pooled_p9010) ? @sprintf(";  pooled 90/10 %.3f [3.737]", m.pooled_p9010) : "")
        @printf(io, "              non-BA teaching track: BA − non-BA teachers' mean log wage %+.3f (men %+.3f, women %+.3f) [0.482 (0.0094), men 0.202, women 0.566; controlled 0.354];  teacher 90/10 all %.3f (men %.3f, women %.3f) [2.773 (2.597, 2.830)], BA teachers only %.3f (men %.3f, women %.3f) [2.417 (2.449, 2.406)]\n",
                m.teacher_ba_wage_gap, m.teacher_ba_wage_gap_m, m.teacher_ba_wage_gap_f, m.teacher_p9010, m.teacher_p9010_g[1], m.teacher_p9010_g[2],
                m.teacher_p9010_ba, m.teacher_p9010_ba_g[1], m.teacher_p9010_ba_g[2])
    end
    # round 2026-10-16: who teaches where (diagnostics from `teacher_sorting_stats`, `sorting_moments`; none is a target)
    if hasproperty(m, :gap_logh_teacher_work) && isfinite(m.gap_logh_teacher_work) && isfinite(m.gap_logh_teacher_res)
        sg(x) = @sprintf("%+.4f", x)
        rw(a, b) = sg(a) * " | " * sg(b)
        @printf(io, "  teacher sorting, suburb − city, residence | workplace: mean log h_T %s;  mean log h_eff %s;  log E[h_eff^(β/σ)] %s;  mean log pay %s;  mean ln ŷ %s;  log mean pay %s | %s\n",
                rw(m.gap_logh_teacher_res, m.gap_logh_teacher_work), rw(m.gap_logheff_teacher_res, m.gap_logheff_teacher_work),
                rw(m.gap_Heff_teacher_res, m.gap_Heff_teacher_work), rw(m.gap_logpay_teacher_res, m.gap_logpay_teacher_work),
                rw(m.gap_lyhat_teacher_res, m.gap_lyhat_teacher_work), sg(m.gap_salary_res), sg(m.gap_salary))
        @printf(io, "              Δ mean log h_T = log Q_l + α log z + α log ϵ + φ log s + η log e: residence %s;  workplace %s  [quadrature residual %.1e];  born in the suburb: residence %s, workplace %s\n",
                join((sg(m[Symbol("gap_logh_teacher_res_", t)]) for t in (:Q, :z, :eps, :s, :e)), " + "),
                join((sg(m[Symbol("gap_logh_teacher_work_", t)]) for t in (:Q, :z, :eps, :s, :e)), " + "),
                m.decomp_resid_max, fmtvec(m.share_born_suburb_res), fmtvec(m.share_born_suburb_work))
        isfinite(m.gap_logQ_teach_quality) && isfinite(m.gap_logQ_teach_staffing) &&
            @printf(io, "              Δ log Q_T = %+.4f = quality σ Δlog E[h_eff^(β/σ)] %+.4f + staffing σ Δlog(teachers per pupil) %+.4f  (check %.1e)\n",
                    m.gap_logQ_teach, m.gap_logQ_teach_quality, m.gap_logQ_teach_staffing, m.gap_logQ_teach_check)
        isfinite(m.gap_logpay_teacher_res_fixedw) && isfinite(m.gap_logh_teacher_res_fixedw) &&
            @printf(io, "              at a given workplace, suburban − city residents: mean log pay %+.4f, mean log h_T %+.4f;  (r, w) cells 11, 21, 12, 22: teachers %s, mean log h_T %s, mean log pay %s\n",
                    m.gap_logpay_teacher_res_fixedw, m.gap_logh_teacher_res_fixedw, fmtvec(m.teacher_rw_mass),
                    fmtvec(m.teacher_rw_logh), fmtvec(m.teacher_rw_logpay))
        isfinite(m.gap_logy_mkt_nonteach_res) && isfinite(m.gap_lyhat_mkt_nonteach_res) &&
            @printf(io, "              non-teaching market workers, suburb − city residents: mean log income %+.4f, mean ln ŷ %+.4f;  mean ln ŷ of teachers %.4f, of market workers %.4f, of non-teaching market workers %.4f\n",
                    m.gap_logy_mkt_nonteach_res, m.gap_lyhat_mkt_nonteach_res, m.mean_log_yhat_teach, m.mean_log_yhat_mkt,
                    m.mean_log_yhat_nonteach_mkt)
    end
    if hasproperty(m, :lyT_fixed) && hasproperty(m, :cfr_rho) && (isfinite(m.lyT_fixed) || m.cfr_rho != 1.0)
        sw = String[]
        isfinite(m.lyT_fixed) &&
            push!(sw, @sprintf("lyT_fixed = %.4f (every teacher node's income taste is graded at this ln ŷ; inert while By = 0)", m.lyT_fixed))
        m.cfr_rho != 1.0 &&
            push!(sw, @sprintf("cfr_rho = %.3f (CFR bridge: sd of log effectiveness = sd(log h_eff)/ρ = %.4f)", m.cfr_rho, m.cfr_sd_eff))
        println(io, "  diagnostic switches of round 2026-10-16: ", join(sw, ";  "))
    end
    @printf(io, "  consumption support: min C = %.3e;  %d (node, l′) pairs at the floor, chosen by a population share %.1e (no feasible location: %.1e)\n",
            m.minC, m.n_binding, m.floor_mass, m.infeasible_mass)
end
diagnostics_table(m) = diagnostics_table(stdout, m)

"""
    evaluate(model, θ, targets; io = stdout) -> (sol, m)

Solve at ϑ and print the moment and diagnostics tables.
"""
function evaluate(model::SpatialModel, θ, targets; io::IO = stdout)
    t = @elapsed (sol, m) = solve_theta(model, θ)
    println(io, "\nϑ = (", join((@sprintf("%s = %.5f", n, x) for (n, x) in zip(THETA_NAMES, θ)), ", "), ")")
    @printf(io, "  solved in %.1fs;  (e, s) residual %.1e, at most %d passes in the final sweep; relaxed step at up to %d nodes in a sweep\n",
            t, sol.s_resid, sol.s_passes, sol.s_relaxed)
    moment_table(io, m, targets)
    diagnostics_table(io, m)
    return sol, m
end

"""
    moment_jacobian(model, θ, targets; bounds, step, free, map = map)

Central-difference Jacobian of the standardized moments m_k/se_k with respect to
the free parameters (those with lo < hi in `bounds`; steps are 0.001 box widths), for
the T6 identification checks. Columns are also scaled
by the width of the search box, so singular values compare like units.
- `sv`, `cond`: singular values and condition number of the scaled Jacobian;
- `R2[j]`: share of column j reproduced by a combination of the other columns
  (§4's projection check). With seven moments and five other columns it is
  mechanically close to 1, so also read
- `unique[j]`: the norm of the part of column j the others cannot reproduce, in
  standardized moments per box width: the effect no other parameter can mimic.
Pass `map = pmap` to spread the 2·length(free) solves across workers. At a bound of
`bounds` the difference is one-sided (the point and a double step inward), so no
solve leaves the box.
"""
function moment_jacobian(model::SpatialModel, θ, targets;
                         bounds = THETA_BOUNDS, step = 1e-3 .* [hi - lo for (lo, hi) in bounds],
                         free = findall(first.(bounds) .< last.(bounds)), map = map)
    act  = active_targets(targets)
    keys = [t.key for t in act]
    se   = [t.se for t in act]
    free = collect(free)
    # (upper, lower) offsets in steps: central, or one-sided at a bound
    offs(j) = θ[j] + step[j] > bounds[j][2] ? (0, -2) :
              θ[j] - step[j] < bounds[j][1] ? (2, 0) : (1, -1)
    shifted(j, k) = (x = collect(Float64, θ); x[j] += k * step[j]; x)
    pts  = vcat([shifted(j, offs(j)[1]) for j in free], [shifted(j, offs(j)[2]) for j in free])
    vals = map(x -> (m = solve_theta(model, x)[2]; [Float64(m[k]) for k in keys]), pts)
    nf   = length(free)
    J    = reduce(hcat, [(vals[i] .- vals[i + nf]) ./ ((offs(j)[1] - offs(j)[2]) * step[j])
                         for (i, j) in enumerate(free)])
    Js   = (J ./ se) .* [hi - lo for (lo, hi) in bounds[free]]'
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
