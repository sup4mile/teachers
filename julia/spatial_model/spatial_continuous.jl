# =============================================================================
# Spatial occupational-choice model with teacher spillovers + altruism.
# General-equilibrium stationary solution: L locations, I occupations.
# Specialized to ε = 0 (parental cost share) and solved with the continuous
# idiosyncratic shocks integrated analytically rather than on a discrete grid.
#
# Two structural facts make the high-dimensional choice problem tractable:
#
#   1. NON-TEACHING COLLAPSE.  For a non-teaching occupation i≠T the whole
#      sub-problem (optimal goods investment e, location probabilities, value)
#      depends on (i, z, ε_i) only through the income capacity z^α·X_{O,i}, where
#          X_{O,i} ≡ Θ_{i,g} ε_i^α ,   Θ_{i,g} ≡ (1−τ^ω_{i,g}) A_i .
#      (Income is (1−τ^ω)A_i·h = Q_l z^α X_{O,i} s_O^φ e^η, so consumption, the
#      e-FOC, and the value all see i only through X_{O,i}.) Hence the best
#      non-teaching option is argmax_i X_{O,i}, and the scalar
#          X_O* ≡ max_{i≠T} Θ_{i,g} ε_i^α
#      is an EXACT sufficient statistic.  With ε_i iid lognormal each Θ_{i,g}ε_i^α
#      is lognormal, so X_O* has the closed-form max-CDF F_{X_O*}(x)=∏_{i≠T}Φ_i(x).
#      This collapses the I-dimensional choice integral to ONE dimension in X_O*
#      plus the scalar teaching shock ε_T, killing the N_ε^I curse of dimensionality.
#
#   2. CONTINUOUS OCCUPATION MARGIN.  A young agent teaches iff W_T(ε_T) >
#      W_O(X_O*).  Both occupational values are monotone increasing in their
#      scalar argument, so the indifference locus is found by inverting the value
#      functions (spline inversion).  Integrating the hard teach/don't-teach
#      indicator against the atomless density makes every aggregate (altruism term
#      Λ, teaching share, H̃_T, tax base/bill, migration kernel) a smooth function
#      of the threshold, so the household map contracts without occupation-choice tremble.
#
#   3. z STAYS DISCRETE.  The persistent ability state is Rouwenhorst-discretised;
#      Φ_l(z), the spatial distribution, is the dominant left eigenvector of the
#      migration × AR(1) kernel.
#
# The remaining pieces are the closed-form time share s, the 1-D goods margin e,
# the endogenous Φ_l(z), the teacher aggregator H̃_T, the local balanced budgets
# t_l, and the outer (H̃_T, M, t) general-equilibrium fixed point.
#
# Symbol note:  ε  (\varepsilon) = parental cost share — fixed at 0 here.
#               ϵ  (\epsilon)    = idiosyncratic ability shock (continuous).
# =============================================================================

using LinearAlgebra
using Printf
using Optim
using Dierckx
using QuadGK
using Distributions

# Smooth log barrier used in place of log(consumption): equals log(Y) for Y > δ,
# then switches to a linear extension below δ so the optimiser always sees a
# finite, nonzero gradient instead of log's −∞ / blow-up near zero.
@inline smooth_log(Y, δ = 1e-6) = Y > δ ? log(Y) : log(δ) + (Y - δ) / δ

# Warm-glow kernel f(h') in the altruism term.  ψ = 1 is the log kernel (the
# ψ→1 limit of the power form, and the pre-sorting nesting case), ψ < 1 makes f
# CONVEX in log h' so the marginal value of a better school rises with the child's
# ability — the §3.2 sorting mechanism, and ψ = 0 at the BASELINE (see below).
# href is a fixed normalisation (see `mean_child_h`):
# h'^{1-ψ} is not scale-free the way log h' is, so href must be set to the
# benchmark mean h̄ or ψ silently rescales the whole altruism term.
# The isapprox guard handles ψ EXACTLY 1 (the 0/0); `expm1` handles the band
# just below it, where (x^(1-ψ) − 1) is a cancelling difference of two numbers
# both ≈ 1 — expm1 of the combined exponent keeps full precision there.  Same
# treatment in `child_Ef`'s power branch, which evaluates the same kernel.
@inline warm(h, ψ, href) = isapprox(ψ, 1.0; atol = 1e-10) ? log(h / href) :
                           expm1((1 - ψ) * log(h / href)) / (1 - ψ)

# -----------------------------------------------------------------------------
# 1. Parameters.  All structural constants of the model live in this immutable
#    struct; a single instance `p::Params` is threaded through every routine.
#    occ_θ is retained as a field for interface compatibility but is UNUSED here
#    (the continuous integration removes the need for the occupation tremble).
#
# THE BASELINE SORTS ON ABILITY.  `Params()` ships BOTH sorting mechanisms on:
#   §3.1  the moving cost is denominated in GOODS (`mcost`), not utility
#         (`τmove = 0`).  A fixed goods loss is regressive in ability, so the
#         mobility–ability slope is positive in BOTH occupational branches, not
#         just among teachers.
#   §3.2  the warm-glow kernel is the power form at ψ = 0 (linear in h'), not
#         the log form at ψ = 1, so the marginal value of a better school rises
#         with the child's ability.
# Together these make the post-migration mean-ability gap POSITIVE
# (`idx_work_signed ≈ +0.08`, ~150 % of it within-branch) where the pre-sorting
# parameterization had it negative and entirely compositional.  Turning both off
# — `nesting_params()` — recovers the original model exactly.
#
# The three constants below are FROZEN normalizations, not equilibrium objects.
# They were computed once from `nesting_params()` solved at Nz = 5,
# nϵT = nXO = 48, and `sorting_mechanism_tests` re-derives and re-checks them:
#   HBAR_BENCH  = mean_child_h    of that solution  → the kernel's `href`
#   CBAR_BENCH  = mean_consumption of that solution → `Cbar`
#   MCOST_BENCH = (1 − e^{−τmove/μ})·C̄ at τmove = 0.20, μ = 1 — the goods loss
#                 whose utility value at C̄ is exactly the old utility cost.  So
#                 the switch to goods denomination adds NO new free parameter.
# -----------------------------------------------------------------------------
const TAUMOVE_BENCH = 0.20                                   # the pre-sorting utility moving cost
const HBAR_BENCH    = 0.1685483614                           # h̄ at the nesting benchmark
const CBAR_BENCH    = 0.1619724217                           # C̄ at the nesting benchmark
const MCOST_BENCH   = (1 - exp(-TAUMOVE_BENCH / 1.0)) * CBAR_BENCH   # μ = 1 ⇒ 0.0293606

Base.@kwdef struct Params
    T::Int = 1                                   # index of the teaching occupation within 1:I
    A::Vector{Float64} = [NaN, 1.5]              # length I; non-teaching productivity A_i (A[T] unused)
    α::Float64 = 0.30                            # ability (z, ε) elasticity of human capital
    φ::Float64 = 0.40                            # time-investment (s) elasticity of human capital
    η::Float64 = 0.20                            # goods-investment (e) elasticity of human capital
    σ::Float64 = 0.25                            # teacher-spillover curvature (enters the quality index Q)
    γ::Float64 = 0.80                            # human-capital elasticity of the teaching wage
    κ::Vector{Float64} = [0.75, 0.9]             # length L; local teaching-wage shifter per location
    μ::Float64 = 1.00                            # weight on log consumption in flow utility
    λ::Float64 = 0.70                            # altruism strength on the warm-glow kernel f(h'); see `warm` and ψ below
    occ_θ::Float64 = 0.0                         # occupation-choice tremble scale — ignored in continuous solver
    B::Vector{Float64}      = [0.0, 0.1]         # length L; location amenity value
    σν::Float64             = 0.20               # scale of the Gumbel location-taste shock
    τmove::Matrix{Float64}  = zeros(2, 2)        # L×L; UTILITY cost of moving born-location → work-location. 0 at the baseline: superseded by `mcost`
    τω::Matrix{Float64} = [0.0 0.00; 0.0 0.1]    # I×2; labor-income wedge by (occupation, gender), income keeps (1−τω)
    τe::Matrix{Float64} = [0.0 0.00; 0.0 0.00]   # I×2; education barrier by (occupation, gender), cost scales as (1+τe)
    ρz::Float64 = 0.9                            # AR(1) persistence of log ability z
    σξ::Float64 = 0.20                           # std of the AR(1) innovation to log z
    σϵ::Float64 = 0.30                           # std of the iid log idiosyncratic shock ε
    β::Float64 = 0.15                            # exponent linking human capital to the teacher aggregator H̃_T
    Mtot::Float64 = 2.0                          # total population mass (both genders, all locations)
    # --- sorting mechanisms; both are ON at the baseline, and `nesting_params()`
    #     turns both off to recover the original model exactly ---
    ψ::Float64    = 0.0                          # warm-glow curvature; ψ = 1 ⇒ log kernel (the nesting case)
    href::Float64 = HBAR_BENCH                   # fixed reference h̄ in the kernel; a NORMALIZATION, not re-solved in GE
    mcost::Matrix{Float64} = [0.0 MCOST_BENCH; MCOST_BENCH 0.0]  # L×L; GOODS cost of moving l → l' (mcost[l,l] = 0)
    Cbar::Float64 = CBAR_BENCH                   # fixed reference consumption `mcost` is denominated in (see redenominate_move_cost)
end


"""
    Params(I, L; ...)

Convenience constructor for `I` occupations and `L` locations (location 1 = base;
locations 2..L share `B_other`, `κ_other`, and pairwise moving costs `τmove_off`
(utility-denominated) and `mcost_off` (goods-denominated); non-teaching
occupations share productivity `A_other`).  Any field can still be overridden by
keyword via `kwargs...`.

Defaults match the 2×2 baseline: `τmove_off = 0`, `mcost_off = MCOST_BENCH`.
Note `mcost` is an ABSOLUTE goods amount, so it is NOT scale-free across
economies — a parameterization whose consumption level differs materially from
`CBAR_BENCH` should re-derive its own cost with `redenominate_move_cost` rather
than inherit this default (see `generalized_dims_test`).
"""
function Params(I::Int, L::Int;
                 T::Int = 1,
                 A_other::Float64 = 1.5,
                 κ_base::Float64 = 0.75, κ_other::Float64 = 0.9,
                 B_base::Float64 = 0.0,  B_other::Float64 = 0.1,
                 τmove_off::Float64 = 0.0,
                 mcost_off::Float64 = MCOST_BENCH,
                 τω_other::Float64 = 0.1,
                 kwargs...)
    A = fill(A_other, I)
    A[T] = NaN
    κ = vcat(κ_base, fill(κ_other, L - 1))
    B = vcat(B_base, fill(B_other, L - 1))
    τmove = fill(τmove_off, L, L)
    mcost = fill(mcost_off, L, L)
    for l in 1:L
        τmove[l, l] = 0.0
        mcost[l, l] = 0.0
    end
    τω = zeros(I, 2)
    I ≥ 2 && (τω[2, 2] = τω_other)
    τe = zeros(I, 2)
    defaults = (; T, A, κ, B, τmove, mcost, τω, τe)
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

"""
    nesting_params(; kwargs...)

The PRE-SORTING parameterization both mechanisms nest at: the moving cost back in
UTILITY units (`τmove = TAUMOVE_BENCH` off-diagonal, `mcost = 0`) and the log warm
glow (`ψ = 1`, `href = 1`).  Everything else is a baseline default.

This is the regression baseline — it reproduces the model as it stood before the
sorting mechanisms were switched on — and the point the frozen normalizations
`HBAR_BENCH` / `CBAR_BENCH` are measured at.  Sorting here is NEGATIVE and
entirely compositional (`idx_work_signed < 0`, within-branch gradient ≈ 0), which
is exactly what the baseline moves away from.
"""
function nesting_params(; L::Int = 2, kwargs...)
    τmove = fill(TAUMOVE_BENCH, L, L)
    for l in 1:L
        τmove[l, l] = 0.0
    end
    defaults = (; τmove, mcost = zeros(L, L), ψ = 1.0, href = 1.0, Cbar = 1.0)
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

# -----------------------------------------------------------------------------
# 1b. Optimal time share s (Proposition 2′).  Combining the s-FOC with the e-FOC
#     (both evaluated under the logit location average) cancels every state
#     variable except ONE scalar,
#         Ξ ≡ Σ_{l'} π_{l'|l} · m_{l,l'} / C_{l'},
#     the choice-weighted moving cost per unit of consumption.  So s stays an
#     EXPLICIT policy — no extra numerical optimisation — but is now
#     state-dependent whenever the goods-denominated moving cost m is nonzero.
#     At m ≡ 0 (Ξ = 0) this collapses to the old state-independent closed form.
#     Teaching and non-teaching differ because teaching income runs through the
#     extra human-capital exponent γ.
# -----------------------------------------------------------------------------
@inline function s_of(i::Int, p::Params, Ξ::Float64 = 0.0)
    b = p.μ * p.φ * (1 + Ξ)
    return i == p.T ? b * p.γ / (b * p.γ + 1 - p.γ * p.η) : b / (b + 1 - p.η)
end

"Ξ = 0 reference time share — the state-independent value the solver warm-starts
 from, and the ONLY value s takes when there is no goods moving cost (i.e. at
 `nesting_params()`, not at the baseline). Retained under its old name for callers
 written against the pre-Proposition-2′ interface."
@inline s_closed(i::Int, p::Params) = s_of(i, p, 0.0)

# -----------------------------------------------------------------------------
# 2. Discretization helpers
# -----------------------------------------------------------------------------

"Rouwenhorst discretization of log z' = ρ log z + σ ξ, ξ~N(0,1). Returns z, Π."
function rouwenhorst(N::Int, ρ::Float64, σ::Float64)
    p = (1 + ρ) / 2
    Π = [p 1-p; 1-p p]
    for n in 3:N
        Πp = Π
        Π = zeros(n, n)
        @views Π[1:n-1, 1:n-1] .+= p     .* Πp
        @views Π[1:n-1, 2:n  ] .+= (1-p) .* Πp
        @views Π[2:n  , 1:n-1] .+= (1-p) .* Πp
        @views Π[2:n  , 2:n  ] .+= p     .* Πp
        @views Π[2:n-1, :]     ./= 2
    end
    σz   = σ / sqrt(1 - ρ^2)
    ψ    = sqrt(N - 1) * σz
    logz = range(-ψ, ψ; length = N)
    return exp.(logz), Π
end

"Stationary distribution π (πΠ = π) of a row-stochastic matrix Π via power iteration."
function stationary(Π; tol = 1e-12, max_iter = 10000)
    π = fill(1 / size(Π, 1), size(Π, 1))       # start from the uniform distribution
    for _ in 1:max_iter
        π_new = vec(π' * Π)
        if maximum(abs.(π_new .- π)) < tol
            return π_new ./ sum(π_new)
        end
        π = π_new
    end
    return π ./ sum(π)
end

"Force a vector to be strictly increasing (for safe spline inversion)."
function make_increasing(v::AbstractVector)
    w = collect(float.(v))
    @inbounds for i in 2:length(w)
        w[i] = max(w[i], w[i-1] + 1e-12)
    end
    return w
end

# -----------------------------------------------------------------------------
# 3. Grids + continuous-shock machinery.
#
#   * z, Πz             : persistent state (discrete, Rouwenhorst).
#   * dT                : law of the teaching shock ε_T ~ LogNormal(-σϵ²/2, σϵ).
#   * ϵTgrid            : quantile grid for ε_T on which teaching policies live.
#   * dO[i,g]           : law of X_{O,i} = Θ_{i,g} ε_i^α  (lognormal).
#   * XOgrid[:,g]       : log-spaced grid for X_O* = max_i X_{O,i}, per gender.
#   * logΘ[i,g]         : log Θ_{i,g}, used to reconstruct log h_O from X_O*.
# -----------------------------------------------------------------------------
function build_grids(p::Params; H̃T, M, t, Nz = 5, nϵT = 64, nXO = 64,
                     q_lo = 1e-5, q_hi = 1 - 1e-6)
    I = length(p.A); L = length(p.B)
    z, Πz = rouwenhorst(Nz, p.ρz, p.σξ)        # ability nodes z and their transition matrix Πz
    Mfloor = 1e-8
    # Floor population mass in Q to avoid 0/0 -> NaN when a location collapses.
    Q = (2 .* H̃T ./ max.(M, Mfloor)) .^ p.σ    # per-location teacher quality index Q_l (from H̃_T and mass M)

    # Teaching shock ε_T: unit-mean lognormal (mean of log is −σϵ²/2). ϵTgrid is a
    # quantile grid spanning (q_lo, q_hi) on which all teaching policies are stored.
    dT     = LogNormal(-p.σϵ^2 / 2, p.σϵ)
    ϵTgrid = quantile.(dT, range(q_lo, q_hi, nϵT))

    nonteach = [i for i in 1:I if i != p.T]    # indices of the non-teaching occupations

    # The X_O* sufficient-statistic collapse requires the education barrier τe to
    # be COMMON across non-teaching occupations: heterogeneous (1−τω)A_i are fine
    # (they enter Θ_{i,g}, hence X), but a τe that varies with i≠T would make the
    # value depend on i beyond X_{O,i}, and W_nonteach (which uses τe of the
    # first non-teaching occupation) would silently be wrong.  Fail loudly.
    i0_nonteach = nonteach[1]
    common_τe = all(p.τe[i, g] == p.τe[i0_nonteach, g] for i in nonteach, g in 1:2)
    @assert common_τe "this solver requires a common τe across non-teaching occupations " *
                      "(a τe that varies with i≠T breaks the X_O* sufficient-statistic collapse)"
    # dO[i,g] is the law of the per-occupation income capacity X_{O,i}=Θ_{i,g}ε_i^α,
    # and logΘ[i,g] stores log Θ_{i,g} (needed later to back out log h_O from X_O*).
    dO   = Array{LogNormal{Float64},2}(undef, I, 2)
    logΘ = fill(NaN, I, 2)
    for g in 1:2, i in nonteach
        Θ        = (1 - p.τω[i, g]) * p.A[i]    # Θ_{i,g} = (1−τ^ω) A_i, the non-stochastic capacity shifter
        logΘ[i, g] = log(Θ)
        # log X_{O,i} = log Θ + α log ε ~ N(log Θ − ασϵ²/2, (ασϵ)²)
        dO[i, g] = LogNormal(log(Θ) - p.α * p.σϵ^2 / 2, p.α * p.σϵ)
    end

    # XOgrid[:,g]: log-spaced grid for the sufficient statistic X_O* = max_i X_{O,i},
    # spanning the (q_lo, q_hi) quantiles across the non-teaching occupations.
    XOgrid = Array{Float64,2}(undef, nXO, 2)
    for g in 1:2
        lo = minimum(quantile(dO[i, g], q_lo) for i in nonteach)
        hi = maximum(quantile(dO[i, g], q_hi) for i in nonteach)
        XOgrid[:, g] = exp.(range(log(lo), log(hi), nXO))
    end

    # The returned named tuple `gr` is the read-only grid/shock bundle passed everywhere:
    #   z, Πz, Nz   ability nodes, transition, count      Q, t   quality index and tax rates
    #   I, L        #occupations, #locations              dT, ϵTgrid  teaching-shock law and grid
    #   dO, XOgrid, logΘ  non-teaching-shock laws, X_O* grid, log Θ
    #   nonteach    non-teaching occupation indices        nϵT, nXO   shock grid sizes
    return (; z, Πz, Nz, Q, t, I, L, dT, ϵTgrid, dO, XOgrid, logΘ, nonteach,
              nϵT, nXO)
end

# CDF of X_O* = max_{i≠T} X_{O,i}: since the X_{O,i} are independent, the max-CDF
# is just the product of the individual lognormal CDFs.
@inline Fxo(gr, g, x) = prod(cdf(gr.dO[i, g], x) for i in gr.nonteach)
# pdf of X_O*, by differentiating the product CDF:  f = F · Σ_i f_i/F_i.
function fxo(gr, g, x)
    F = Fxo(gr, g, x)
    F <= 0.0 && return 0.0
    ratio_sum = 0.0
    for i in gr.nonteach
        Fi = cdf(gr.dO[i, g], x)
        Fi > 0.0 && (ratio_sum += pdf(gr.dO[i, g], x) / Fi)     # accumulate f_i / F_i
    end
    return F * ratio_sum
end

# E[log Θ_{i*} | X_O* = x]: which non-teaching occupation attained the max?
# (For a single non-teaching occupation this is just logΘ of that occupation.)
function ElogΘ(gr, g, x)
    isone(length(gr.nonteach)) && return gr.logΘ[gr.nonteach[1], g]
    weight_sum    = 0.0
    weighted_logΘ = 0.0
    for i in gr.nonteach
        Fi     = max(cdf(gr.dO[i, g], x), 1e-300)
        weight = pdf(gr.dO[i, g], x) / Fi     # ∝ P(i = argmax | max = x)
        weight_sum    += weight
        weighted_logΘ += weight * gr.logΘ[i, g]
    end
    return weight_sum > 0.0 ? weighted_logΘ / weight_sum : gr.logΘ[gr.nonteach[1], g]
end

# E[Θ_{i*}^{−(1−ψ)} | X_O* = x]: the power-kernel analogue of ElogΘ, using the
# identical argmax weights.  Needed because the warm-glow kernel of §6 is applied
# to h_O itself, not to log h_O, so with I > 2 non-teaching occupations averaging
# a point estimate of log h_O would carry a Jensen bias.
function EΘpow(gr, g, x, ψ)
    isone(length(gr.nonteach)) && return exp(-(1 - ψ) * gr.logΘ[gr.nonteach[1], g])
    weight_sum = 0.0
    weighted   = 0.0
    for i in gr.nonteach
        Fi     = max(cdf(gr.dO[i, g], x), 1e-300)
        weight = pdf(gr.dO[i, g], x) / Fi     # ∝ P(i = argmax | max = x)
        weight_sum += weight
        weighted   += weight * exp(-(1 - ψ) * gr.logΘ[i, g])
    end
    return weight_sum > 0.0 ? weighted / weight_sum :
           exp(-(1 - ψ) * gr.logΘ[gr.nonteach[1], g])
end

# E[ε_{i*} | X_O* = x]: the idiosyncratic draw in the occupation the agent ACTUALLY
# works in, i.e. the ε behind the ability a_{i*} = z·ε_{i*} of §2.  The sufficient
# statistic is X_{O,i*} = Θ_{i*} ε_{i*}^α, so ε_{i*} = (x/Θ_{i*})^{1/α} and the
# conditional mean is x^{1/α}·E[Θ_{i*}^{−1/α} | x] under the SAME argmax weights as
# ElogΘ.  With a single non-teaching occupation Θ_{i*} is degenerate and this is
# exact; with several it is the conditional mean, which is what the ability
# aggregates (`student_ability_moments`) integrate.
function Eϵ_nonteach(gr, g, x, α)
    xpow = x^(1 / α)
    isone(length(gr.nonteach)) &&
        return xpow * exp(-gr.logΘ[gr.nonteach[1], g] / α)
    weight_sum = 0.0
    weighted   = 0.0
    for i in gr.nonteach
        Fi     = max(cdf(gr.dO[i, g], x), 1e-300)
        weight = pdf(gr.dO[i, g], x) / Fi     # ∝ P(i = argmax | max = x)
        weight_sum += weight
        weighted   += weight * exp(-gr.logΘ[i, g] / α)
    end
    return xpow * (weight_sum > 0.0 ? weighted / weight_sum :
                   exp(-gr.logΘ[gr.nonteach[1], g] / α))
end

# -----------------------------------------------------------------------------
# 4. Stage-2 location values (ε = 0) and the Gumbel log-sum, for each branch.
#
# Teaching: wage ω_{T,l'} = κ_l' h^γ, h = Q_l (z ε_T)^α s_T^φ e^η.
# Non-teaching: income capacity (1−τ^ω)A_i h = Q_l z^α X_O s_O^φ e^η is a
# function of X_O alone (the collapse), so we never name a specific occupation.
# Both add the precomputed altruism term Λ[zi, l'].
# -----------------------------------------------------------------------------
# Teaching branch: given the teaching shock ϵT, ability node zi, birth location l,
# gender g, goods investment e and time share sT, return the per-work-location
# value vector W (utility of living in each l', including the altruism term Λ)
# plus the agent's human capital h.  Wage in l' is κ_l' h^γ; Y is post-tax income
# net of schooling AND of the goods moving cost m_{l,l'} (zero on the diagonal).
function W_teach(ϵT, zi, l, g, e, sT, Λ, p::Params, gr)
    h  = gr.Q[l] * (gr.z[zi] * ϵT)^p.α * sT^p.φ * e^p.η   # human capital produced when teaching
    W  = Vector{Float64}(undef, gr.L)
    for lp in 1:gr.L
        wage = p.κ[lp] * h^p.γ                            # teaching wage in candidate work location lp
        Y    = (1 - gr.t[lp]) * (1 - p.τω[p.T, g]) * wage -
               (1 + p.τe[p.T, g]) * e - p.mcost[l, lp]    # net consumption
        W[lp] = p.μ * smooth_log(Y) + Λ[zi, lp]           # flow utility + altruism value of location lp
    end
    return W, h
end

# Non-teaching branch: same shape as W_teach but driven by the sufficient
# statistic XO = X_O* rather than a named occupation.  `cap` is the income
# capacity (1−τ^ω)A_i·h; net income needs no h^γ term (γ applies only to teaching).
function W_nonteach(XO, zi, l, g, e, sO, Λ, p::Params, gr)
    i0  = gr.nonteach[1]
    τeO = p.τe[i0, g]                          # common τ^e across non-teaching occupations
    cap = gr.Q[l] * gr.z[zi]^p.α * XO * sO^p.φ * e^p.η   # income capacity (1−τ^ω)A_i·h
    W   = Vector{Float64}(undef, gr.L)
    for lp in 1:gr.L
        Y = (1 - gr.t[lp]) * cap - (1 + τeO) * e - p.mcost[l, lp]   # post-tax net consumption
        W[lp] = p.μ * smooth_log(Y) + Λ[zi, lp]
    end
    return W, cap
end

# Gumbel (logit) smoothing over the location choice: given the per-location value
# vector W and birth location l, return (V̄, π) where V̄ is the expected value
# (the σν-scaled log-sum-exp inclusive value) and π is the choice-probability
# vector over work locations.  Amenity B and moving cost τmove enter here; the
# max-subtraction m is the standard log-sum-exp overflow guard.
function logsum_probs(W, l, p::Params)
    L = length(W)
    scores = [(W[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in 1:L]
    smax   = maximum(scores)
    denom  = sum(exp.(scores .- smax))
    return p.σν * (smax + log(denom)), exp.(scores .- smax) ./ denom
end

# Inner 1-D goods margin: choose goods investment e to maximise the location
# log-sum (inclusive value).  The objective V̄(u), u = log e, is smooth and
# single-peaked, and its first two derivatives in u are available in closed form,
# so we solve the FOC with a WARM-STARTED, SAFEGUARDED NEWTON iteration instead of
# a cold golden-section search — the same optimum, but 1–2 evaluations per node
# once the previous sweep's e is used as the seed.  With a ≡ the effective goods
# exponent in income (η for non-teaching, ηγ for teaching because the teaching
# wage runs through h^γ), mc ≡ (1+τe)·e, m_{l'} ≡ mcost[l,l'] the goods moving
# cost, y^g_{l'} the gross income term and Y_{l'} = y^g_{l'} − mc − m_{l'}:
#     w_{l'}  ≡ dW_{l'}/du = μ [ a − ((1−a) mc − a m_{l'}) / Y_{l'} ]
#     w'_{l'} ≡ dw_{l'}/du = μ (a² y^g_{l'} − mc) / Y_{l'} − w_{l'}²/μ
#     g  ≡ dV̄/du   = Σ_{l'} π_{l'} w_{l'}
#     g' ≡ d²V̄/du² = Var_π(w)/σν + Σ_{l'} π_{l'} w'_{l'}
# (at m ≡ 0 the second line reduces to the earlier −(μa − w)(μ − w)/μ.)
# where π is the logit location weight (logsum_probs).  `income(e)` returns
# (W, Y, mc, ygross): the per-location value vector, the net-consumption vector,
# the common marginal schooling cost, and the per-location gross income.  The
# analytic w uses the smooth_log slope 1/Y, valid only where the barrier is slack
# (Y > δ); if any Y ≤ δ, or a Newton step leaves the box, or the objective is
# locally non-concave, we fall back to the robust bracketed search.  At a sensible
# solution the barrier is slack, so the fallback essentially never triggers once
# warm-started.  Both solvers return (e, V̄, Y, π): the location weights π and the
# net consumption Y are needed by the outer s iteration below.
function solve_e_bracket(income, l, p::Params; lo = log(1e-8), hi = log(1e4))
    f(le) = -logsum_probs(income(exp(le))[1], l, p)[1]
    res = optimize(f, lo, hi)
    e = exp(Optim.minimizer(res))
    W, Y, _, _ = income(e)
    _, π = logsum_probs(W, l, p)
    return e, -Optim.minimum(res), Y, π                  # e, V̄, Y, π
end

function solve_e_newton(income, a, l, p::Params, u0;
                        lo = log(1e-8), hi = log(1e4), δ = 1e-6,
                        maxit = 30, gtol = 1e-10)
    u = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    mvec = @view p.mcost[l, :]
    for _ in 1:maxit
        e = exp(u)
        W, Y, mc, ygross = income(e)
        Vbar, π = logsum_probs(W, l, p)
        any(≤(δ), Y) && return solve_e_bracket(income, l, p; lo, hi)   # barrier active → 1/Y invalid
        w    = @. p.μ * (a - ((1 - a) * mc - a * mvec) / Y)   # dW_{l'}/du
        g    = dot(π, w)
        varw = dot(π, w .^ 2) - g^2
        wp   = @. p.μ * (a^2 * ygross - mc) / Y - w^2 / p.μ   # dw_{l'}/du
        gp   = varw / p.σν + dot(π, wp)
        abs(g) < gtol && return e, Vbar, Y, π
        gp ≥ -1e-14 && return solve_e_bracket(income, l, p; lo, hi)    # not concave here
        unew = u - g / gp
        (isfinite(unew) && lo ≤ unew ≤ hi) || return solve_e_bracket(income, l, p; lo, hi)
        u = unew
    end
    return solve_e_bracket(income, l, p; lo, hi)         # no convergence → robust fallback
end

# Outer fixed point in the time share s (Proposition 2′).  Given s the goods
# margin is the Newton solve above; given the resulting (π, Y) the s-FOC is the
# explicit `s_of(i, p, Ξ)`.  The feedback is stabilising (higher s ⇒ higher h ⇒
# higher C ⇒ lower Ξ ⇒ lower s) and contracts by a factor of ~0.03 per pass, so no
# damping is needed: a WARM-STARTED sweep converges in 2 passes and only the very
# first (cold) sweep of a solve needs ~6–8.  With mcost ≡ 0, Ξ ≡ 0 and the loop
# exits after ONE pass at the old constant s, which is what makes the whole
# mechanism exactly nesting.
#
# S_PASSES_MAX reports the COST (pass count, dominated by the cold first sweep);
# S_RESID_MAX reports the ACCURACY (|s_{n+1} − s_n| at exit) and is what
# verify_solution actually audits.  They are running maxima over every cell, so
# solve_ge zeroes them a SECOND time just before the final clean solve and
# stores the result in the returned sol as (s_resid, s_passes): the audit then
# describes the policies it is handed, not the worst of the GE path, and stays
# correct when several solutions are verified after all of them have solved.
const S_PASSES_MAX = Ref(0)
const S_RESID_MAX  = Ref(0.0)
const S_MAXIT = 12
const S_TOL   = 1e-10

# Ξ = Σ_{l'} π_{l'|l} m_{l,l'} / C_{l'}.  The floor on Y only guards against a
# division blow-up if the smooth_log barrier is binding; `verify_solution` is what
# reports that situation.
@inline function xi_from(π, Y, mvec)
    Ξ = 0.0
    @inbounds for lp in eachindex(π)
        Ξ += π[lp] * mvec[lp] / max(Y[lp], 1e-8)
    end
    return Ξ
end

function solve_e_s(make_income, i::Int, a, l, p::Params, u0, s0)
    s = (isfinite(s0) && 0.0 < s0 < 1.0) ? s0 : s_of(i, p, 0.0)
    mvec = @view p.mcost[l, :]
    local e, Vbar
    passes = 0
    resid  = 0.0
    for _ in 1:S_MAXIT
        passes += 1
        e, Vbar, Y, π = solve_e_newton(make_income(s), a, l, p, u0)
        snew  = s_of(i, p, xi_from(π, Y, mvec))
        resid = abs(snew - s)
        s = snew
        resid < S_TOL && break
    end
    passes > S_PASSES_MAX[] && (S_PASSES_MAX[] = passes)
    resid  > S_RESID_MAX[]  && (S_RESID_MAX[]  = resid)
    return e, Vbar, s
end

# Teaching branch: h = base·e^η, wage_{l'} = κ_{l'} h^γ, so net income runs
# through e^{ηγ}; effective goods exponent a = ηγ.  Seeded at u0 = log(e_prev)
# and s0 = s_prev.  Returns (e, V̄, s).
function solve_e_teach(ϵT, zi, l, g, Λ, p::Params, gr, u0, s0)
    τeT  = p.τe[p.T, g]
    coef = 1 - p.τω[p.T, g]
    cap0 = gr.Q[l] * (gr.z[zi] * ϵT)^p.α
    make_income(s) = let base = cap0 * s^p.φ
        function income(e)
            h  = base * e^p.η
            mc = (1 + τeT) * e
            W  = Vector{Float64}(undef, gr.L)
            Y  = Vector{Float64}(undef, gr.L)
            yg = Vector{Float64}(undef, gr.L)
            for lp in 1:gr.L
                ygl = (1 - gr.t[lp]) * coef * p.κ[lp] * h^p.γ
                y   = ygl - mc - p.mcost[l, lp]
                yg[lp] = ygl
                Y[lp]  = y
                W[lp]  = p.μ * smooth_log(y) + Λ[zi, lp]
            end
            return W, Y, mc, yg
        end
    end
    return solve_e_s(make_income, p.T, p.η * p.γ, l, p, u0, s0)
end

# Non-teaching branch: income capacity cap = K·e^η; effective goods exponent a = η.
function solve_e_nonteach(XO, zi, l, g, Λ, p::Params, gr, u0, s0)
    i0  = gr.nonteach[1]
    τeO = p.τe[i0, g]
    K0  = gr.Q[l] * gr.z[zi]^p.α * XO
    make_income(s) = let K = K0 * s^p.φ
        function income(e)
            cap = K * e^p.η
            mc  = (1 + τeO) * e
            W   = Vector{Float64}(undef, gr.L)
            Y   = Vector{Float64}(undef, gr.L)
            yg  = Vector{Float64}(undef, gr.L)
            for lp in 1:gr.L
                ygl = (1 - gr.t[lp]) * cap
                y   = ygl - mc - p.mcost[l, lp]
                yg[lp] = ygl
                Y[lp]  = y
                W[lp]  = p.μ * smooth_log(y) + Λ[zi, lp]
            end
            return W, Y, mc, yg
        end
    end
    return solve_e_s(make_income, i0, p.η, l, p, u0, s0)
end

# -----------------------------------------------------------------------------
# 5. Occupation margin: teach iff W_T(ε_T) > W_O(X_O*).
#
# Both W_T(·) and W_O(·) are monotone increasing in their scalar argument, so the
# indifference locus is found by inverting the value functions.
#   teach_wt(ε_T)     = P(X_O* < threshold) = F_{X_O*}( W_O^{-1}(W_T(ε_T)) )
#   nonteach_wt(X_O*) = P(ε_T  < threshold) = F_{ε_T}( W_T^{-1}(W_O(X_O*)) )
# i.e. the probability that a young agent with the given shock chooses to teach
# (resp. not teach), obtained by mapping the value through the OTHER branch's
# inverse and evaluating that branch's shock CDF.  Clamped to {0,1} outside the
# value ranges.
# -----------------------------------------------------------------------------
# Bundles everything needed to evaluate the occupation margin for one (gender g)
# at one household state: the value functions and their inverses as splines, the
# raw value vectors (for range clamping), and the grids object gr.
struct ChoiceMaps
    splWT::Spline1D; splWO::Spline1D          # W_T(ε_T) and W_O(X_O*) as splines over their grids
    splWTinv::Spline1D; splWOinv::Spline1D    # the inverse maps value → ε_T and value → X_O*
    WTv::Vector{Float64}; WOv::Vector{Float64} # raw value vectors on the grids (define the invertible range)
    g::Int                                    # gender index (1 or 2) these maps were built for
    gr                                        # grids/shock-distribution bundle (see build_grids)
end

# Theory says W_T(ε_T) and W_O(X_O*) are strictly increasing; make_increasing
# only patches TINY numerical wiggles (1e-12 bumps).  If a value vector is
# materially non-monotone the inverse spline can misplace the occupation
# threshold, so detect it and warn instead of silently patching.
const NONMONO_WARNINGS = Ref(0)
function warn_nonmonotone(v::AbstractVector, name::String)
    drop = 0.0
    @inbounds for i in 2:length(v)
        d = v[i-1] - v[i]
        d > drop && (drop = d)
    end
    scale = max(abs(v[end] - v[1]), 1e-12)
    if drop > 1e-8 * scale && NONMONO_WARNINGS[] < 5
        NONMONO_WARNINGS[] += 1
        suffix = NONMONO_WARNINGS[] == 5 ? " (suppressing further warnings)" : ""
        @printf("  ⚠ %s materially non-monotone (max drop %.3e vs range %.3e): occupation threshold from spline inversion may be off%s\n",
                name, drop, scale, suffix)
        flush(stdout)
    end
end

function choice_maps(WTv, WOv, g, gr)
    warn_nonmonotone(WTv, "W_T(ε_T)")
    warn_nonmonotone(WOv, "W_O(X_O*)")
    splWT = Spline1D(gr.ϵTgrid, WTv)
    splWO = Spline1D(gr.XOgrid[:, g], WOv)
    splWTinv = Spline1D(make_increasing(WTv), gr.ϵTgrid; k = 1)
    splWOinv = Spline1D(make_increasing(WOv), gr.XOgrid[:, g]; k = 1)
    return ChoiceMaps(splWT, splWO, splWTinv, splWOinv, WTv, WOv, g, gr)
end

@inline function teach_wt(cm::ChoiceMaps, ϵT)
    wT = cm.splWT(ϵT)
    wT <= cm.WOv[1]   && return 0.0
    wT >= cm.WOv[end] && return 1.0
    return Fxo(cm.gr, cm.g, cm.splWOinv(wT))
end
@inline function nonteach_wt(cm::ChoiceMaps, XO)
    wO = cm.splWO(XO)
    wO <= cm.WTv[1]   && return 0.0
    wO >= cm.WTv[end] && return 1.0
    return cdf(cm.gr.dT, cm.splWTinv(wO))
end

# Integrate a teaching-side grid quantity qT(ε_T) and a non-teaching-side grid
# quantity qO(X_O*) against their densities, weighted by the choice indicator.
#
# The quadrature runs over the quantile grids; the shock mass OUTSIDE the grids
# (≈ q_lo + tail beyond q_hi per branch) is added back via a constant-
# extrapolation tail correction  (tail mass) · q(endpoint) · w(endpoint),  so no
# probability mass is silently dropped and the migration-kernel rows sum to ≈ 1.
# (Constant extrapolation slightly underestimates upper-tail quantities that
# grow in the shock, e.g. h^{β/σ}; with q_hi = 1−1e-6 the residual is ~1e-6.)
# Returns (∫ teaching, ∫ non-teaching).
function integrate_choice(qTvec, qOvec, cm::ChoiceMaps, gr)
    g = cm.g

    # Teaching side: ∫ qT(ε_T)·f_T(ε_T) over the region where the agent teaches.
    ϵlo, ϵhi = gr.ϵTgrid[1], gr.ϵTgrid[end]
    splqT = Spline1D(gr.ϵTgrid, qTvec)
    teach_integrand(ϵ) = begin
        w = teach_wt(cm, ϵ)
        w <= 0.0 ? 0.0 : pdf(gr.dT, ϵ) * splqT(ϵ) * w
    end
    IT  = quadgk(teach_integrand, ϵlo, ϵhi)[1]
    IT += cdf(gr.dT, ϵlo)  * qTvec[1]   * teach_wt(cm, ϵlo)   # lower tail
    IT += ccdf(gr.dT, ϵhi) * qTvec[end] * teach_wt(cm, ϵhi)   # upper tail

    # Non-teaching side: ∫ qO(X_O*)·f_{X_O*}(x) over the region where the agent does not teach.
    xlo, xhi = gr.XOgrid[1, g], gr.XOgrid[end, g]
    splqO = Spline1D(gr.XOgrid[:, g], qOvec)
    nonteach_integrand(x) = begin
        w = nonteach_wt(cm, x)
        w <= 0.0 ? 0.0 : fxo(gr, g, x) * splqO(x) * w
    end
    IO  = quadgk(nonteach_integrand, xlo, xhi)[1]
    IO += Fxo(gr, g, xlo)       * qOvec[1]   * nonteach_wt(cm, xlo)  # lower tail
    IO += (1 - Fxo(gr, g, xhi)) * qOvec[end] * nonteach_wt(cm, xhi)  # upper tail
    return IT, IO
end

# -----------------------------------------------------------------------------
# 6. The altruism term Λ.  For a child born in (l', g', z') the expected
# warm glow integrates the child's optimal occupation choice:
#   E[f(h')] = ∫ f_T(ε_T) f(h_T(ε_T)) teach_wt(ε_T) dε_T
#            + ∫ f_{X_O*}(x) E[f(h_O) | x] nonteach_wt(x) dx,
# with the kernel f = `warm` (log at ψ = 1; power otherwise) and
#   log h_O = log[Q_l z^α s_O^φ e^η] + log X_O* − log Θ_{i*}.
# On the teaching side h_T is on the grid directly.  On the non-teaching side
# only the DISTRIBUTION of Θ_{i*} given X_O* is known, so the expectation is
# taken in the right space for each kernel: E[log Θ_{i*}] for ψ = 1 (ElogΘ) and
# E[Θ_{i*}^{−(1−ψ)}] otherwise (EΘpow).  With a single non-teaching occupation
# both are degenerate and the two branches agree exactly.
# Then Λ[zi, l'] = ½ Σ_{g'} Σ_{z'} Π_z[zi,z'] λ E[f(h') | l', g', z'].
#
# Λ is investment- and occupation-independent for ANY kernel, so ψ leaves the
# household FOCs (Propositions 1–3) untouched: only this section changes.
# -----------------------------------------------------------------------------
function child_Ef(l, g, zi, hh, p, gr)
    cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
    # teaching: h_T directly available on the ε_T grid.
    fTvec = warm.(@view(hh.hT[g, l, zi, :]), p.ψ, p.href)
    # non-teaching: reconstruct h_O on the X_O grid.  s_O is now a per-node
    # policy, so the s term stays INSIDE the k-loop.
    base0 = log(gr.Q[l]) + p.α * log(gr.z[zi])
    if isapprox(p.ψ, 1.0; atol = 1e-10)
        fOvec = [base0 + p.φ * log(hh.sO[g, l, zi, k]) + p.η * log(hh.eO[g, l, zi, k]) +
                 log(gr.XOgrid[k, g]) - ElogΘ(gr, g, gr.XOgrid[k, g]) - log(p.href)
                 for k in 1:gr.nXO]
    else
        # E[(h_O/href)^{1−ψ}] = (Q z^α s^φ e^η X_O / href)^{1−ψ} · E[Θ_{i*}^{−(1−ψ)}]
        fOvec = [begin
                     loga = base0 + p.φ * log(hh.sO[g, l, zi, k]) +
                            p.η * log(hh.eO[g, l, zi, k]) +
                            log(gr.XOgrid[k, g]) - log(p.href)
                     # exp(u) = (h_O/href)^{1−ψ}·E[Θ^{−(1−ψ)}]; expm1 on the
                     # COMBINED exponent avoids the cancellation as ψ → 1, and
                     # makes the limit u/(1−ψ) → loga − E[logΘ] visibly the
                     # ψ = 1 branch above.  EΘpow is a positive weighted mean,
                     # so its log is always well defined.
                     u = (1 - p.ψ) * loga + log(EΘpow(gr, g, gr.XOgrid[k, g], p.ψ))
                     expm1(u) / (1 - p.ψ)
                 end
                 for k in 1:gr.nXO]
    end
    IT, IO = integrate_choice(fTvec, fOvec, cm, gr)
    return IT + IO
end

function compute_lambda(hh, p, gr)
    (; Nz, L, Πz) = gr
    Ef = zeros(L, 2, Nz)                       # [l', g', z']
    for lp in 1:L, g in 1:2, zpi in 1:Nz
        Ef[lp, g, zpi] = child_Ef(lp, g, zpi, hh, p, gr)
    end
    Λ = zeros(Nz, L)
    for zi in 1:Nz, lp in 1:L
        acc = 0.0
        for g in 1:2, zpi in 1:Nz
            acc += 0.5 * Πz[zi, zpi] * p.λ * Ef[lp, g, zpi]
        end
        Λ[zi, lp] = acc
    end
    return Λ
end

# -----------------------------------------------------------------------------
# 7. Household block.  Iterate the intergenerational fixed point
#    V → h' → Λ → V.  s and the goods margin e are both refined each sweep (s is
#    explicit given Ξ, so this costs one extra pass of the Newton solve per node
#    only when a goods moving cost is active).  No occupation tremble: the
#    atomless integration keeps Λ (hence V) continuous in V, so the map contracts.
#
#    Λ0 warm-starts the altruism term (e.g. from the previous GE iteration's
#    household solution): the sweep count is governed by how far Λ starts from
#    its fixed point, so this saves most of the household time in the GE loop.
#
#    eT0/eO0 (and sT0/sO0) likewise warm-start the goods and time policies.  Each
#    per-node solve is seeded at the previous sweep's (e, s) (and, on the first
#    sweep, at the previous GE iteration's when supplied), so the Newton e-solver
#    converges in 1–2 steps instead of re-searching from scratch every sweep.
# -----------------------------------------------------------------------------
function solve_household(p::Params, gr; tol = 1e-6, maxit = 500,
                          verbose = false, print_every = 25,
                          Λ0 = nothing, eT0 = nothing, eO0 = nothing,
                          sT0 = nothing, sO0 = nothing)
    (; Nz, L, nϵT, nXO) = gr
    i0 = gr.nonteach[1]
    # Policy/value arrays indexed [gender, birth-location, ability node, shock node].
    # Teaching side lives on the ϵT grid, non-teaching side on the X_O* grid.
    eT = zeros(2, L, Nz, nϵT); WT = zeros(2, L, Nz, nϵT); hT = zeros(2, L, Nz, nϵT)  # goods e, value W_T, human capital h_T
    eO = zeros(2, L, Nz, nXO); WO = zeros(2, L, Nz, nXO)                              # goods e, value W_O (non-teaching)
    sT = fill(s_of(p.T, p, 0.0), 2, L, Nz, nϵT)   # time share, now a per-node policy
    sO = fill(s_of(i0,  p, 0.0), 2, L, Nz, nXO)
    # Seed the goods/time policies from a previous GE iteration when their shape
    # matches; otherwise the first e-solve falls back to the bracketed search
    # (u0 non-finite) and s starts at its Ξ = 0 value.
    (eT0 !== nothing && size(eT0) == size(eT)) && (eT .= eT0)
    (eO0 !== nothing && size(eO0) == size(eO)) && (eO .= eO0)
    (sT0 !== nothing && size(sT0) == size(sT)) && (sT .= sT0)
    (sO0 !== nothing && size(sO0) == size(sO)) && (sO .= sO0)
    Λ  = (Λ0 !== nothing && size(Λ0) == (Nz, L)) ? copy(Λ0) : zeros(Nz, L)            # altruism term Λ[zi, l'] (warm-started if given)
    hh = (; eT, WT, hT, eO, WO, sT, sO, Λ)     # the household solution bundle, updated each sweep
    t0 = time()
    local err
    for it in 1:maxit
        WT0 = copy(WT); WO0 = copy(WO)     # previous-sweep values, for the convergence check
        # Refresh every (gender, location, ability, shock) policy given the current Λ.
        # u0 = log(previous e) warm-starts each Newton solve; a zero (unset) entry
        # gives a non-finite seed, which the solver treats as a cold bracketed solve.
        for g in 1:2, l in 1:L, zi in 1:Nz
            for k in 1:nϵT
                e, V̄, s = solve_e_teach(gr.ϵTgrid[k], zi, l, g, Λ, p, gr,
                                        log(eT[g, l, zi, k]), sT[g, l, zi, k])
                eT[g, l, zi, k] = e
                sT[g, l, zi, k] = s
                WT[g, l, zi, k] = log(1 - s) + V̄   # add the log time-endowment term log(1−s) to the location value
                hT[g, l, zi, k] = gr.Q[l] * (gr.z[zi] * gr.ϵTgrid[k])^p.α * s^p.φ * e^p.η
            end
            for k in 1:nXO
                e, V̄, s = solve_e_nonteach(gr.XOgrid[k, g], zi, l, g, Λ, p, gr,
                                           log(eO[g, l, zi, k]), sO[g, l, zi, k])
                eO[g, l, zi, k] = e
                sO[g, l, zi, k] = s
                WO[g, l, zi, k] = log(1 - s) + V̄
            end
        end
        Λ = compute_lambda(hh, p, gr)      # re-derive the altruism term from the updated policies
        hh = (; eT, WT, hT, eO, WO, sT, sO, Λ)
        err = max(maximum(abs, WT .- WT0), maximum(abs, WO .- WO0))
        if verbose && (it == 1 || it % print_every == 0 || err < tol)
            @printf("    HH %4d/%d  err=%.3e  elapsed=%.1fs\n", it, maxit, err, time() - t0)
            flush(stdout)
        end
        err < tol && break
    end
    return hh
end

# -----------------------------------------------------------------------------
# 8. Location-choice probabilities on the grids (used for migration/aggregates).
# -----------------------------------------------------------------------------
function location_probs(hh, p::Params, gr)
    (; Nz, L, nϵT, nXO) = gr
    πT = zeros(2, L, Nz, nϵT, L)
    πO = zeros(2, L, Nz, nXO, L)
    for g in 1:2, l in 1:L, zi in 1:Nz
        for k in 1:nϵT
            W, _ = W_teach(gr.ϵTgrid[k], zi, l, g, hh.eT[g, l, zi, k],
                           hh.sT[g, l, zi, k], hh.Λ, p, gr)
            _, πi = logsum_probs(W, l, p)
            @views πT[g, l, zi, k, :] .= πi
        end
        for k in 1:nXO
            W, _ = W_nonteach(gr.XOgrid[k, g], zi, l, g, hh.eO[g, l, zi, k],
                              hh.sO[g, l, zi, k], hh.Λ, p, gr)
            _, πi = logsum_probs(W, l, p)
            @views πO[g, l, zi, k, :] .= πi
        end
    end
    return πT, πO
end

# -----------------------------------------------------------------------------
# 9. Endogenous spatial distribution Φ_l(z): the stationary mass of agents in
#    each (location, ability) cell.  It is the dominant left eigenvector of the
#    combined kernel P = (migration probabilities) × (AR(1) ability transition),
#    scaled to total population.  The migration probabilities are the location
#    choices integrated over the occupation margin (continuous integrals here).
# -----------------------------------------------------------------------------
function stationary_phi(πT, πO, hh, p::Params, gr)
    (; Nz, L, Πz) = gr
    # πbar[l0, zi, lp] = probability an agent born in l0 with ability zi works in lp,
    # averaged over gender and integrated over the occupation choice.
    πbar = zeros(L, Nz, L)
    for l0 in 1:L, zi in 1:Nz
        for g in 1:2
            cm = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
            for lp in 1:L
                IT, IO = integrate_choice(πT[g, l0, zi, :, lp],
                                          πO[g, l0, zi, :, lp], cm, gr)
                πbar[l0, zi, lp] += 0.5 * (IT + IO)   # ½ per gender (equal masses)
            end
        end
    end
    # Assemble the joint (location, ability) transition kernel P and take its
    # stationary distribution; idx flattens (l, zi) into a single row/column index.
    n = L * Nz
    P = zeros(n, n)
    idx(l, zi) = (l - 1) * Nz + zi
    for l0 in 1:L, zi in 1:Nz, l in 1:L, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]   # migrate then age (AR(1) draws z')
    end
    φ = stationary(P); φ .*= p.Mtot / 2        # scale eigenvector to the per-gender population mass
    return reshape(φ, Nz, L), πbar
end

# -----------------------------------------------------------------------------
# 10. Aggregates by work location: the teacher aggregator H̃_T, the population
#     mass M, and the balanced-budget tax rate t.  Each is a policy quantity
#     integrated over the occupation margin and summed against the mass Φ.
# -----------------------------------------------------------------------------
function aggregates(Φ, πT, πO, hh, p::Params, gr)
    (; Nz, L) = gr
    expo = p.β / p.σ                           # exponent on h in the teacher aggregator H̃_T
    HT = zeros(L); Inc = zeros(L); Wbill = zeros(L)   # teacher aggregator, taxable income, teacher wage bill
    for l0 in 1:L, zi in 1:Nz
        cell_mass = 0.5 * Φ[zi, l0]                         # mass in this (birth-loc, ability) cell, ½ per gender
        cell_mass == 0.0 && continue
        for g in 1:2
            cm   = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
            hTv  = hh.hT[g, l0, zi, :]          # teacher human capital across the ε_T grid
            capO = [gr.Q[l0] * gr.z[zi]^p.α * gr.XOgrid[k, g] * hh.sO[g, l0, zi, k]^p.φ *
                    hh.eO[g, l0, zi, k]^p.η for k in 1:gr.nXO]   # non-teaching income capacity across X_O* grid
            for lp in 1:L
                # Build the per-shock integrand vectors weighted by the migration
                # probability into lp, then integrate each over the teach/not-teach margin.
                wageT = p.κ[lp] .* hTv .^ p.γ                        # teaching wage in lp
                qHT   = πT[g, l0, zi, :, lp] .* hTv .^ expo          # contribution to H̃_T (teachers only)
                qIncT = πT[g, l0, zi, :, lp] .* (1 - p.τω[p.T, g]) .* wageT   # taxable teaching income
                qWb   = πT[g, l0, zi, :, lp] .* wageT               # gross teacher wage bill
                qIncO = πO[g, l0, zi, :, lp] .* capO                # taxable non-teaching income
                zeroO = zeros(gr.nXO)                               # non-teaching side is zero for teacher-only terms

                ITh, _   = integrate_choice(qHT,   zeroO, cm, gr)   # ∫ over teachers only
                ITi, IOi = integrate_choice(qIncT, qIncO, cm, gr)   # income, both branches
                ITw, _   = integrate_choice(qWb,   zeroO, cm, gr)   # wage bill, teachers only

                HT[lp]    += cell_mass * ITh
                Inc[lp]   += cell_mass * (ITi + IOi)
                Wbill[lp] += cell_mass * ITw
            end
        end
    end
    M = [2 * sum(@view Φ[:, l]) for l in 1:L]   # total mass per location (×2 to undo the per-gender ½)
    # Balanced local budget: tax rate = teacher wage bill / taxable income, clamped to [0, 0.9].
    t = [Inc[lp] > 0 ? clamp(Wbill[lp] / Inc[lp], 0.0, 0.9) : 0.0 for lp in 1:L]
    return HT, M, t
end

# -----------------------------------------------------------------------------
# 11. Outer general-equilibrium fixed point over (H̃_T, M, t).
# -----------------------------------------------------------------------------
fmtvec(v) = "(" * join((@sprintf("%.4f", x) for x in v), ", ") * ")"
reldiff(nw, od; fl = 1e-3) = maximum(abs.(nw .- od) ./ max.(abs.(od), fl))

function solve_ge(p::Params = Params();
                  Nz = 5, nϵT = 64, nXO = 64, damping = 0.3, tol = 1e-5,
                  maxit = 1000, hh_tol = 1e-6, hh_maxit = 1000, verbose = true,
                  init = nothing, q_lo = 1e-5, q_hi = 1 - 1e-6)
    I = length(p.A); L = length(p.B)
    @assert length(p.κ) == L && size(p.τmove) == (L, L) && size(p.mcost) == (L, L) &&
            size(p.τω) == (I, 2) && size(p.τe) == (I, 2) && 1 <= p.T <= I "Params dims inconsistent"
    @assert all(p.mcost[l, l] == 0.0 for l in 1:L) "mcost must have a zero diagonal (staying put is free)"
    S_PASSES_MAX[] = 0            # zeroed again before the final clean solve below
    S_RESID_MAX[]  = 0.0

    if verbose
        @printf("  [ε=0, continuous]  I=%d, L=%d, Nz=%d  shock grids: nϵT=%d, nXO=%d\n",
                I, L, Nz, nϵT, nXO)
        @printf("  Ξ=0 reference s: ")
        for i in 1:I
            @printf("%s s=%.4f   ", i == p.T ? "T" : "occ$i", s_of(i, p, 0.0))
        end
        @printf("| ψ=%.3f href=%.4f  max mcost=%.4f", p.ψ, p.href, maximum(p.mcost))
        println(); flush(stdout)
    end

    HT = fill(1.0, L); M = fill(p.Mtot / L, L); t = fill(0.10, L)
    # Household warm start: Λ from the previous GE iteration (and, optionally,
    # a previously solved equilibrium passed via `init` for continuation runs).
    Λ0 = nothing; eT0 = nothing; eO0 = nothing; sT0 = nothing; sO0 = nothing
    if init !== nothing
        if length(init.HT) == L
            HT .= init.HT; M .= init.M; t .= init.t
            size(init.hh.Λ) == (Nz, L) && (Λ0 = copy(init.hh.Λ))
            eT0 = copy(init.hh.eT); eO0 = copy(init.hh.eO)   # seed goods policy for continuation
            sT0 = copy(init.hh.sT); sO0 = copy(init.hh.sO)   # and the time policy
        elseif verbose
            @printf("  ⚠ warm-start ignored: init has L=%d ≠ %d\n", length(init.HT), L)
        end
    end
    local hh, Pi_T, Pi_O, Φ, πbar
    Mfloor = 1e-8
    for it in 1:maxit
        # One GE iteration: rebuild grids at the current (H̃_T, M, t), solve the
        # household problem, derive location choices, the stationary distribution,
        # and the implied aggregates.  @elapsed times each stage for the log line.
        gr = build_grids(p; H̃T = HT, M = M, t = t, Nz, nϵT, nXO, q_lo, q_hi)
        t_hh  = @elapsed (hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit,
                                               verbose = verbose && it == 1,
                                               Λ0, eT0, eO0, sT0, sO0))
        Λ0 = hh.Λ                              # carry Λ forward to warm-start the next iteration's household solve
        eT0 = hh.eT; eO0 = hh.eO               # carry the goods policy forward too
        sT0 = hh.sT; sO0 = hh.sO               # and the time policy
        t_pi  = @elapsed ((Pi_T, Pi_O) = location_probs(hh, p, gr))
        t_phi = @elapsed ((Φ, πbar) = stationary_phi(Pi_T, Pi_O, hh, p, gr))
        t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, Pi_T, Pi_O, hh, p, gr))

        collapsed = Mn .< Mfloor
        # A collapsed location has no meaningful balanced-budget tax; hold its tax
        # fixed so a 0/undefined tax-base ratio does not create spurious oscillation.
        tn[collapsed] .= t[collapsed]
        t_err = any(.!collapsed) ? maximum(abs, tn[.!collapsed] .- t[.!collapsed]) : 0.0
        # Convergence = max relative change in H̃_T and M, and absolute change in t.
        err = max(reldiff(HTn, HT), reldiff(Mn, M), t_err)
        # Damped update: move a fraction `damping` toward the new guess for stability.
        @. HT = (1 - damping) * HT + damping * HTn
        @. M  = (1 - damping) * M  + damping * Mn
        @. t  = (1 - damping) * t  + damping * tn

        if verbose
            @printf("GE %3d  err=%.3e  H̃_T=%s  M=%s  t=%s  [hh=%.1fs pi=%.1fs phi=%.1fs agg=%.1fs]\n",
                    it, err, fmtvec(HT), fmtvec(M), fmtvec(t), t_hh, t_pi, t_phi, t_agg)
            flush(stdout)
        end
        if err < tol
            verbose && println("GE aggregates converged.")
            break
        end
    end

    gr = build_grids(p; H̃T = HT, M = M, t = t, Nz, nϵT, nXO, q_lo, q_hi)
    # Final clean solve at the converged (H̃_T, M, t) so the returned policies,
    # location probabilities, and distribution are all mutually consistent.
    # Zero the s diagnostics so they measure THIS solve only: warm-started off the
    # converged Λ, it is the tightest sweep of the run and the only one whose
    # residual describes the returned policies.
    S_PASSES_MAX[] = 0
    S_RESID_MAX[]  = 0.0
    hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit, Λ0, eT0, eO0, sT0, sO0)
    Pi_T, Pi_O = location_probs(hh, p, gr)
    Φ, πbar = stationary_phi(Pi_T, Pi_O, hh, p, gr)
    # The solution bundle `sol`: household policies hh, location probs Pi_T/Pi_O,
    # spatial distribution Φ, migration matrix πbar, aggregates HT/M/t, grids gr, params p,
    # and the s-fixed-point diagnostics of the clean solve (audited by verify_solution).
    return (; hh, Pi_T, Pi_O, Φ, πbar, HT, M, t, gr, p,
              s_resid = S_RESID_MAX[], s_passes = S_PASSES_MAX[])
end

# -----------------------------------------------------------------------------
# 12. Diagnostics
# -----------------------------------------------------------------------------

"""
    verify_solution(sol; δ = 1e-6, tclamp = 0.9)

Post-solve audit that the numerical safeguards are SLACK at the solution:
(i) no (state, work-location) cell has consumption Y ≤ δ, i.e. the smooth_log
barrier never binds and every stored value is the true log; (ii) no local tax
rate sits on the clamp [0, tclamp]; (iii) the s fixed point of Proposition 2′
converged (residual below `s_atol`).  Returns true if clean.  The s diagnostics
come from `sol` itself (`sol.s_resid`, `sol.s_passes`, recorded by `solve_ge` on
the final clean solve), so verifying several solutions after all of them have
solved reports each one's own residual; a `sol` predating those fields falls
back to the globals.  The pass count is reported for cost, not audited: it is a
max over cells and is set by the first (cold) sweep of the clean solve; a
warm-started sweep converges in 2.

The goods moving cost `mcost` subtracts a FIXED amount from a mover's
consumption, so it is the low-ability movers that push the barrier first — check
(i) is the binding constraint on how large `mcost` can be calibrated.  At the
baseline (m = MCOST_BENCH ≈ 0.0294) min C ≈ 0.025, so roughly 1.5× MCOST_BENCH is
the largest clean bump; 2× already drives cells negative.  That is why
`comparative_statics`'s high-move-cost case scales m by 1.5 and not more.
"""
function verify_solution(sol; δ = 1e-6, tclamp = 0.9, s_atol = 1e-8)
    (; hh, t, gr, p) = sol
    # Prefer the per-solve diagnostics stored in `sol`; fall back to the globals
    # for solutions saved before solve_ge started attaching them.
    s_resid  = hasproperty(sol, :s_resid)  ? sol.s_resid  : S_RESID_MAX[]
    s_passes = hasproperty(sol, :s_passes) ? sol.s_passes : S_PASSES_MAX[]
    i0 = gr.nonteach[1]
    minY = Inf
    n_binding = 0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        for k in 1:gr.nϵT                          # teaching branch
            e = hh.eT[g, l, zi, k]
            h = gr.Q[l] * (gr.z[zi] * gr.ϵTgrid[k])^p.α * hh.sT[g, l, zi, k]^p.φ * e^p.η
            for lp in 1:gr.L
                Y = (1 - gr.t[lp]) * (1 - p.τω[p.T, g]) * p.κ[lp] * h^p.γ -
                    (1 + p.τe[p.T, g]) * e - p.mcost[l, lp]
                minY = min(minY, Y); Y <= δ && (n_binding += 1)
            end
        end
        for k in 1:gr.nXO                          # non-teaching branch
            e   = hh.eO[g, l, zi, k]
            cap = gr.Q[l] * gr.z[zi]^p.α * gr.XOgrid[k, g] * hh.sO[g, l, zi, k]^p.φ * e^p.η
            for lp in 1:gr.L
                Y = (1 - gr.t[lp]) * cap - (1 + p.τe[i0, g]) * e - p.mcost[l, lp]
                minY = min(minY, Y); Y <= δ && (n_binding += 1)
            end
        end
    end
    ok = true
    if s_resid > s_atol
        @printf("  ⚠ the s fixed point left a residual of %.2e (> %.1e) after %d passes: the Ξ iteration did not converge\n",
                s_resid, s_atol, s_passes)
        ok = false
    end
    if n_binding > 0
        @printf("  ⚠ smooth_log barrier BINDING at %d (state, l') cells (min C = %.3e ≤ δ = %.1e): values there are barrier-distorted\n",
                n_binding, minY, δ)
        ok = false
    end
    for lp in 1:gr.L
        if t[lp] >= tclamp - 1e-9 || t[lp] <= 0.0
            @printf("  ⚠ tax rate t[%d] = %.4f is ON the clamp [0, %.1f]: budget not truly balanced\n",
                    lp, t[lp], tclamp)
            ok = false
        end
    end
    ok && @printf("  ✓ feasibility audit passed: min C = %.4e > δ = %.1e; taxes interior t = %s; s residual %.1e (≤ %d passes)\n",
                  minY, δ, fmtvec(t), s_resid, s_passes)
    return ok
end

"Copy of `p` with the named fields replaced — `Params` is immutable, so this is
 how a single dial is turned without respecifying the whole struct."
with_params(p::Params; kwargs...) =
    Params(; merge(NamedTuple(k => getfield(p, k) for k in fieldnames(Params)),
                   NamedTuple(kwargs))...)

"""
    mean_child_h(sol)

Φ-weighted, choice-weighted mean child human capital h̄ at `sol`.  This is the
value `href` must be FROZEN at before any ψ sweep: `warm` measures h' in units of
href, and h' ≈ 0.17 ≪ 1 at the benchmark, so leaving href = 1 makes
(h'/href)^{1−ψ} collapse toward 1 and the ψ mechanism look inert.  With
h'/h̄ ≈ 1 the kernel is log(h'/h̄) to first order for EVERY ψ, so Λ₂ − Λ₁ keeps
its ψ = 1 level and ψ does only curvature.

href is a normalization, not an equilibrium object: re-solving it inside the GE
loop would add a fixed point and destroy that interpretation.  Evaluated at
`nesting_params()` this returns `HBAR_BENCH`, the baseline's shipped `href`;
`sorting_mechanism_tests` re-derives it and checks the constant still holds.
"""
function mean_child_h(sol)
    (; hh, Φ, gr, p) = sol
    num = 0.0; den = 0.0
    onesT = ones(gr.nϵT); onesO = ones(gr.nXO)
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2
        cell = 0.5 * Φ[zi, l]
        cell == 0.0 && continue
        cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
        base0 = log(gr.Q[l]) + p.α * log(gr.z[zi])
        hOvec = [exp(base0 + p.φ * log(hh.sO[g, l, zi, k]) + p.η * log(hh.eO[g, l, zi, k]) +
                     log(gr.XOgrid[k, g]) - ElogΘ(gr, g, gr.XOgrid[k, g])) for k in 1:gr.nXO]
        IhT, IhO = integrate_choice(hh.hT[g, l, zi, :], hOvec, cm, gr)
        ImT, ImO = integrate_choice(onesT, onesO, cm, gr)
        num += cell * (IhT + IhO)
        den += cell * (ImT + ImO)
    end
    return num / den
end

"""
    mean_consumption(sol)

Φ-weighted, choice-weighted and location-choice-weighted mean consumption C̄ at
`sol` — the denominator `redenominate_move_cost` uses to convert the utility
moving cost τmove into a goods cost.  Evaluated at `nesting_params()` this returns
`CBAR_BENCH`, the baseline's shipped `Cbar`, from which `MCOST_BENCH` (the
baseline's off-diagonal `mcost`) is derived.
"""
function mean_consumption(sol)
    (; hh, Φ, gr, p) = sol
    i0 = gr.nonteach[1]
    num = 0.0; den = 0.0
    onesT = ones(gr.nϵT); onesO = ones(gr.nXO)
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2
        cell = 0.5 * Φ[zi, l]
        cell == 0.0 && continue
        cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
        # Per node, average net consumption over the location choice π_{l'|l}.
        CT = [sum(sol.Pi_T[g, l, zi, k, lp] *
                  ((1 - gr.t[lp]) * (1 - p.τω[p.T, g]) * p.κ[lp] * hh.hT[g, l, zi, k]^p.γ -
                   (1 + p.τe[p.T, g]) * hh.eT[g, l, zi, k] - p.mcost[l, lp])
                  for lp in 1:gr.L) for k in 1:gr.nϵT]
        CO = [begin
                  cap = gr.Q[l] * gr.z[zi]^p.α * gr.XOgrid[k, g] *
                        hh.sO[g, l, zi, k]^p.φ * hh.eO[g, l, zi, k]^p.η
                  sum(sol.Pi_O[g, l, zi, k, lp] *
                      ((1 - gr.t[lp]) * cap - (1 + p.τe[i0, g]) * hh.eO[g, l, zi, k] -
                       p.mcost[l, lp]) for lp in 1:gr.L)
              end for k in 1:gr.nXO]
        IcT, IcO = integrate_choice(CT, CO, cm, gr)
        ImT, ImO = integrate_choice(onesT, onesO, cm, gr)
        num += cell * (IcT + IcO)
        den += cell * (ImT + ImO)
    end
    return num / den
end

"""
    redenominate_move_cost(p, Cbar) -> Params

Re-express the UTILITY moving cost `p.τmove` as an equivalent GOODS cost, at no
new free parameter:  m = (1 − exp(−τmove/μ))·C̄, which is the goods loss whose
utility value is exactly τmove at consumption C̄.  Returns a copy of `p` with
`τmove = 0`, `mcost = m` and `Cbar` recorded.

This is the §3.1 sorting mechanism: a goods cost is regressive in ability (it
costs a low-consumption agent proportionally more), so the mobility–ability slope
turns positive for EVERYONE, not just for teachers.  Setting `mcost = 0` with
`τmove` restored — i.e. `nesting_params()` — recovers the original model exactly.

The BASELINE already ships the result of this call at `Cbar = CBAR_BENCH`:
`Params().mcost` equals `redenominate_move_cost(nesting_params(), CBAR_BENCH).mcost`
(asserted in `mcost_test`).  Call it directly only to re-derive the cost for an
economy sitting at a different consumption scale, where the absolute
`MCOST_BENCH` would be the wrong size.
"""
function redenominate_move_cost(p::Params, Cbar::Float64)
    L = length(p.B)
    m = [(1 - exp(-p.τmove[l, lp] / p.μ)) * Cbar for l in 1:L, lp in 1:L]
    for l in 1:L
        m[l, l] = 0.0
    end
    return with_params(p; τmove = zeros(L, L), mcost = m, Cbar = Cbar)
end

"Φ-weighted mean and (min, max) range of the time share s across all household
 states, split by branch.  With mcost ≡ 0 both ranges are degenerate at the
 Proposition-2 constants; a spread of a few percentage points is the signature of
 the Proposition-2′ Ξ correction."
function s_summary(sol)
    (; hh, Φ, gr) = sol
    accT = accO = wT = wO = 0.0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        w = 0.5 * Φ[zi, l]
        accT += w * sum(@view hh.sT[g, l, zi, :]) / gr.nϵT; wT += w
        accO += w * sum(@view hh.sO[g, l, zi, :]) / gr.nXO; wO += w
    end
    return (; meanT = accT / wT, meanO = accO / wO,
              rangeT = (minimum(hh.sT), maximum(hh.sT)),
              rangeO = (minimum(hh.sO), maximum(hh.sO)))
end

"Gender-pooled teaching share among young born in l, weighted by endogenous G_l(z)."
function teaching_share_endo(hh, l, Φ, gr, p::Params)
    (; Nz) = gr
    mass = sum(@view Φ[:, l])
    share = 0.0
    for zi in 1:Nz
        Gz = Φ[zi, l] / mass
        for g in 1:2
            cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
            IT, _ = integrate_choice(ones(gr.nϵT), zeros(gr.nXO), cm, gr)
            share += 0.5 * Gz * IT
        end
    end
    return share
end

"""
    student_ability_moments(sol)

Choice-weighted ability moments of the YOUNG (the students) by BIRTH location —
the location an agent is educated in, which is the location `Φ[z, l]` indexes.

Two ability concepts, both conditioned on the occupation the agent actually picks:

  * `z`  — the persistent, inherited component.  Its cross-location distribution
    is `G_l(z) = Φ[:,l]/ΣΦ[:,l]` (printed by `report_ge`) and does NOT depend on
    the occupational choice, so `mean_z` here reproduces the `mean z` row of
    `report_ge`; what the occupational split adds is `mean_z_teach` vs
    `mean_z_non`, i.e. which types teaching draws (γ < 1 ⇒ negatively selected).
  * `a`  — the ability in the CHOSEN occupation, a_i = z·ε_i (§2 of the notes).
    Teachers get a = z·ε_T on the ε_T grid; everyone else gets a = z·ε_{i*} with
    ε_{i*} recovered from the sufficient statistic by `Eϵ_nonteach`.  `a` is the
    object the occupational choice actually operates on: both branches select on
    their own ε, so E[a] exceeds E[z]·E[ε] = E[z] by the selection premium even
    though ε is unit-mean unconditionally.

Every moment is `Φ`-weighted across (ability node, gender) cells and integrated
over the teach/don't-teach margin with `integrate_choice`, so the teaching and
non-teaching pieces carry the equilibrium branch masses: `teach_share` reproduces
`teaching_share_endo` and `mass` reproduces `ΣΦ[:, l]`, both to the choice
quadrature's own accuracy (~1e-3 at nϵT = nXO = 32, the same slack
`report_checks` reports as the π̄ row-sum deviation).

Returns per-location vectors (`mean_a`, `mean_a_teach`, `mean_a_non`, `mean_z`,
`mean_z_teach`, `mean_z_non`, `teach_share`, `mass`), the branch-split ability
mass `ΦT`/`ΦO` (Nz×L, summing to `Φ`), and the population-pooled scalars of the
same names suffixed `_all`.
"""
function student_ability_moments(sol)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, L, nϵT, nXO) = gr
    onesT = ones(nϵT); onesO = ones(nXO)
    ΦT = zeros(Nz, L); ΦO = zeros(Nz, L)        # Φ split by chosen branch
    aT = zeros(L); aO = zeros(L)                # Σ mass · E[a] within branch
    for l in 1:L, zi in 1:Nz, g in 1:2
        cell = 0.5 * Φ[zi, l]                   # ½ per gender (equal masses)
        cell == 0.0 && continue
        cm  = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
        aTv = z[zi] .* gr.ϵTgrid                                        # a = z·ε_T (teaching)
        aOv = [z[zi] * Eϵ_nonteach(gr, g, gr.XOgrid[k, g], p.α) for k in 1:nXO]  # a = z·ε_{i*}
        ImT, ImO = integrate_choice(onesT, onesO, cm, gr)               # branch masses
        IaT, IaO = integrate_choice(aTv,   aOv,   cm, gr)               # branch ability
        ΦT[zi, l] += cell * ImT;  ΦO[zi, l] += cell * ImO
        aT[l]     += cell * IaT;  aO[l]     += cell * IaO
    end
    safe(num, den) = [den[l] > 0 ? num[l] / den[l] : NaN for l in 1:L]
    massT = vec(sum(ΦT; dims = 1)); massO = vec(sum(ΦO; dims = 1))
    zT = [sum(z[zi] * ΦT[zi, l] for zi in 1:Nz) for l in 1:L]           # Σ mass · z, teachers
    zO = [sum(z[zi] * ΦO[zi, l] for zi in 1:Nz) for l in 1:L]           # … and non-teachers
    mass = massT .+ massO
    return (; mass, massT, massO, ΦT, ΦO,
              teach_share  = safe(massT, mass),
              mean_a       = safe(aT .+ aO, mass),
              mean_a_teach = safe(aT, massT),
              mean_a_non   = safe(aO, massO),
              mean_z       = safe(zT .+ zO, mass),
              mean_z_teach = safe(zT, massT),
              mean_z_non   = safe(zO, massO),
              teach_share_all  = sum(massT) / sum(mass),
              mean_a_all       = sum(aT .+ aO) / sum(mass),
              mean_a_teach_all = sum(aT) / sum(massT),
              mean_a_non_all   = sum(aO) / sum(massO),
              mean_z_all       = sum(zT .+ zO) / sum(mass),
              mean_z_teach_all = sum(zT) / sum(massT),
              mean_z_non_all   = sum(zO) / sum(massO))
end

"""
    student_ability_density(sol, l, agrid)

The DENSITY of student ability a = z·ε over `agrid`, for the students born in
location `l`, split into the teaching and non-teaching branches.  Returns
`(fT, fO)`, normalised so that ∫(fT+fO) = 1 and ∫fT = `teach_share[l]`, i.e. the
branch areas are the equilibrium occupation shares (this is what makes the split
"choice-weighted" rather than two conditional densities).

Exact, not a histogram of the quadrature nodes: with z discrete and the shock
grids quantile-spaced, binning nodes produces a comb artifact, so each branch
density is written down in closed form and evaluated wherever it is wanted.

  * teaching:      a = z·ε_T, so the contribution of an ability node is
                   f_ε(a/z)·(1/z)·P(teach | ε_T = a/z).
  * non-teaching:  the agent works in i* = argmax_{i≠T} Θ_i ε_i^α, and the joint
                   density of (i* = i, ε_{i*} = ε) is f_ε(ε)·∏_{j≠i,T} F_{X_j}(Θ_i ε^α)
                   — occupation i wins iff every rival's capacity falls below its
                   own.  Summing that over i and weighting by
                   P(don't teach | X_O* = Θ_i ε^α) gives the branch density; the
                   ∏ collapses to 1 with a single non-teaching occupation.

Both branches are then mass-weighted by the cell mass ½Φ[z, l] and by gender.
"""
function student_ability_density(sol, l, agrid)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, dT, dO, nonteach, logΘ) = gr
    fT = zeros(length(agrid)); fO = zeros(length(agrid))
    for zi in 1:Nz, g in 1:2
        cell = 0.5 * Φ[zi, l]
        cell == 0.0 && continue
        cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
        for (ia, a) in enumerate(agrid)
            ϵ  = a / z[zi]
            fϵ = cell * pdf(dT, ϵ) / z[zi]          # ε ~ F_ε (same law in every occupation)
            fT[ia] += fϵ * teach_wt(cm, ϵ)
            for i in nonteach
                x  = exp(logΘ[i, g]) * ϵ^p.α        # this occupation's capacity X_{O,i}
                pr = 1.0
                for j in nonteach
                    j != i && (pr *= cdf(dO[j, g], x))   # …and it beats every rival
                end
                fO[ia] += fϵ * pr * nonteach_wt(cm, x)
            end
        end
    end
    mass = sum(@view Φ[:, l])
    return fT ./ mass, fO ./ mass
end

function report_ge(sol)
    (; hh, Φ, πbar, HT, M, t, gr, p) = sol
    (; z, Nz, Q, Πz, L) = gr

    println("\n===== GE stationary solution (ε = 0, continuous distributions) =====")
    ss = s_summary(sol)
    @printf("  s (Φ-weighted mean [min, max]):  T=%.4f [%.4f, %.4f]   nonteach=%.4f [%.4f, %.4f]\n",
            ss.meanT, ss.rangeT..., ss.meanO, ss.rangeO...)
    @printf("  H̃_T    = %s\n", fmtvec(HT))
    @printf("  M      = %s      [total %.4f]\n", fmtvec(M), sum(M))
    @printf("  t      = %s\n", fmtvec(t))
    @printf("  Q      = %s\n", fmtvec(Q))

    ergodic = stationary(Πz)
    println("\n  Endogenous ability distribution  G_l(z)   (vs ergodic G*):")
    @printf("    %10s", "z")
    for l in 1:L; @printf(" %10s", "G_$(l)(z)"); end
    @printf(" %10s\n", "G*(z)")
    for zi in 1:Nz
        @printf("    %10.4f", z[zi])
        for l in 1:L; @printf(" %10.4f", Φ[zi, l] / sum(@view Φ[:, l])); end
        @printf(" %10.4f\n", ergodic[zi])
    end
    mean_z(w) = sum(z[zi] * w[zi] for zi in 1:Nz) / sum(w)
    print("    mean z: ")
    for l in 1:L; @printf(" loc %d = %.4f  ", l, mean_z(@view Φ[:, l])); end
    @printf("ergodic = %.4f\n", mean_z(ergodic))

    println("\n  Migration matrix  Π̄[born→work]  (rows sum to 1):")
    for l0 in 1:L
        mass = sum(@view Φ[:, l0])
        row = [sum(Φ[zi, l0] * πbar[l0, zi, l] for zi in 1:Nz) / mass for l in 1:L]
        print("    born $l0 :")
        for l in 1:L; @printf("  ->%d %.4f", l, row[l]); end
        println()
    end

    println("\n  Teaching share among young born in l (endogenous G_l):")
    for l in 1:L
        @printf("    born %d :  %.4f\n", l, teaching_share_endo(hh, l, Φ, gr, p))
    end

    # Student (young) ability by birth location, split by the occupation each
    # agent chooses.  a = z·ε in the CHOSEN occupation; mean z is the same object
    # as the row above, here re-cut by branch to expose teaching's selection on z.
    sa = student_ability_moments(sol)
    println("\n  Student ability by birth location  (choice-weighted; a = z·ε in the chosen occupation):")
    @printf("    %-6s %-8s %-9s %-9s %-11s %-13s %-11s %-13s\n",
            "born", "teach", "mean z", "mean a", "a | teach", "a | nonteach",
            "z | teach", "z | nonteach")
    for l in 1:L
        @printf("    %-6d %-8.4f %-9.4f %-9.4f %-11.4f %-13.4f %-11.4f %-13.4f\n",
                l, sa.teach_share[l], sa.mean_z[l], sa.mean_a[l],
                sa.mean_a_teach[l], sa.mean_a_non[l],
                sa.mean_z_teach[l], sa.mean_z_non[l])
    end
    @printf("    %-6s %-8.4f %-9.4f %-9.4f %-11.4f %-13.4f %-11.4f %-13.4f\n",
            "all", sa.teach_share_all, sa.mean_z_all, sa.mean_a_all,
            sa.mean_a_teach_all, sa.mean_a_non_all,
            sa.mean_z_teach_all, sa.mean_z_non_all)
    println("====================================================================\n")
end

# -----------------------------------------------------------------------------
# 13. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model (ε = 0, continuous) in general equilibrium ...")
    # Solve at the baseline (both sorting mechanisms on) with a moderately fine
    # shock grid.  damping = 0.75 blends each GE update 75 % toward the new
    # fixed-point guess; the baseline tolerates large steps (damping_diagnostic
    # finds every value in 0.2–0.8 usable, with the residual falling as damping
    # rises), so this is the fast end rather than a stability compromise.
    # tol = 1e-6 is the convergence threshold on the (H̃_T, M, t) residual.
    sol = solve_ge(Params(); Nz = 3, nϵT = 48, nXO = 48, damping = 0.75,
                   tol = 1e-6, hh_tol = 1e-6, maxit = 500, hh_maxit = 1000,
                   verbose = true)
    report_ge(sol)
    verify_solution(sol)
    return sol
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
