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

# -----------------------------------------------------------------------------
# 1. Parameters.  All structural constants of the model live in this immutable
#    struct; a single instance `p::Params` is threaded through every routine.
#    occ_θ is retained as a field for interface compatibility but is UNUSED here
#    (the continuous integration removes the need for the occupation tremble).
# -----------------------------------------------------------------------------
Base.@kwdef struct Params
    T::Int = 1                                   # index of the teaching occupation within 1:I
    A::Vector{Float64} = [NaN, 1.5]              # length I; non-teaching productivity A_i (A[T] unused)
    α::Float64 = 0.30                            # ability (z, ε) elasticity of human capital
    φ::Float64 = 0.40                            # time-investment (s) elasticity of human capital
    η::Float64 = 0.20                            # goods-investment (e) elasticity of human capital
    σ::Float64 = 0.25                            # teacher-spillover curvature (enters the quality index Q)
    γ::Float64 = 0.80                            # human-capital elasticity of the teaching wage
    κ::Vector{Float64} = [0.85, 1.0]             # length L; local teaching-wage shifter per location
    μ::Float64 = 1.00                            # weight on log consumption in flow utility
    λ::Float64 = 0.50                            # altruism strength on child log-human-capital, f(h') = log h'
    occ_θ::Float64 = 0.0                         # occupation-choice tremble scale — ignored in continuous solver
    B::Vector{Float64}      = [0.0, 0.2]         # length L; location amenity value
    σν::Float64             = 0.25               # scale of the Gumbel location-taste shock
    τmove::Matrix{Float64}  = [0.0 0.30; 0.30 0.0]  # L×L; utility cost of moving born-location → work-location
    τω::Matrix{Float64} = [0.0 0.00; 0.0 0.1]    # I×2; labor-income wedge by (occupation, gender), income keeps (1−τω)
    τe::Matrix{Float64} = [0.0 0.00; 0.0 0.00]   # I×2; education barrier by (occupation, gender), cost scales as (1+τe)
    ρz::Float64 = 0.70                           # AR(1) persistence of log ability z
    σξ::Float64 = 0.20                           # std of the AR(1) innovation to log z
    σϵ::Float64 = 0.30                           # std of the iid log idiosyncratic shock ε
    β::Float64 = 0.10                            # exponent linking human capital to the teacher aggregator H̃_T
    Mtot::Float64 = 2.0                          # total population mass (both genders, all locations)
end


"""
    Params(I, L; ...)

Convenience constructor for `I` occupations and `L` locations (location 1 = base;
locations 2..L share `B_other`, `κ_other`, and pairwise moving cost `τmove_off`;
non-teaching occupations share productivity `A_other`).  Any field can still be
overridden by keyword via `kwargs...`.
"""
function Params(I::Int, L::Int;
                 T::Int = 1,
                 A_other::Float64 = 1.5,
                 κ_base::Float64 = 0.85, κ_other::Float64 = 1.0,
                 B_base::Float64 = 0.0,  B_other::Float64 = 0.15,
                 τmove_off::Float64 = 0.30,
                 τω_other::Float64 = 0.1,
                 kwargs...)
    A = fill(A_other, I)
    A[T] = NaN
    κ = vcat(κ_base, fill(κ_other, L - 1))
    B = vcat(B_base, fill(B_other, L - 1))
    τmove = fill(τmove_off, L, L)
    for l in 1:L
        τmove[l, l] = 0.0
    end
    τω = zeros(I, 2)
    I ≥ 2 && (τω[2, 2] = τω_other)
    τe = zeros(I, 2)
    defaults = (; T, A, κ, B, τmove, τω, τe)
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

# -----------------------------------------------------------------------------
# 1b. Closed-form optimal time share s.  With ε = 0 the time-investment FOC has a
#     closed-form solution that is constant per occupation (independent of the
#     state), so it is computed once here rather than in the inner optimisation.
#     Teaching and non-teaching have different formulas because teaching income
#     runs through the extra human-capital exponent γ.
# -----------------------------------------------------------------------------
@inline function s_closed(i::Int, p::Params)
    if i == p.T
        return p.μ * p.φ * p.γ / ((p.μ * p.φ - p.η) * p.γ + 1)
    else
        return p.μ * p.φ / (p.μ * p.φ + 1 - p.η)
    end
end

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
    Q = (2 .* H̃T ./ M) .^ p.σ                  # per-location teacher quality index Q_l (from H̃_T and mass M)

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

# -----------------------------------------------------------------------------
# 4. Stage-2 location values (ε = 0) and the Gumbel log-sum, for each branch.
#
# Teaching: wage ω_{T,l'} = κ_l' h^γ, h = Q_l (z ε_T)^α s_T^φ e^η.
# Non-teaching: income capacity (1−τ^ω)A_i h = Q_l z^α X_O s_O^φ e^η is a
# function of X_O alone (the collapse), so we never name a specific occupation.
# Both add the precomputed altruism term Λ[zi, l'].
# -----------------------------------------------------------------------------
# Teaching branch: given the teaching shock ϵT, ability node zi, birth location l,
# gender g, and goods investment e, return the per-work-location value vector W
# (utility of living in each l', including the altruism term Λ) plus the agent's
# human capital h.  Wage in l' is κ_l' h^γ; Y is post-tax income net of schooling.
function W_teach(ϵT, zi, l, g, e, Λ, p::Params, gr)
    sT = s_closed(p.T, p)
    h  = gr.Q[l] * (gr.z[zi] * ϵT)^p.α * sT^p.φ * e^p.η   # human capital produced when teaching
    W  = Vector{Float64}(undef, gr.L)
    for lp in 1:gr.L
        wage = p.κ[lp] * h^p.γ                            # teaching wage in candidate work location lp
        Y    = (1 - gr.t[lp]) * (1 - p.τω[p.T, g]) * wage - (1 + p.τe[p.T, g]) * e  # net consumption
        W[lp] = p.μ * smooth_log(Y) + Λ[zi, lp]           # flow utility + altruism value of location lp
    end
    return W, h
end

# Non-teaching branch: same shape as W_teach but driven by the sufficient
# statistic XO = X_O* rather than a named occupation.  `cap` is the income
# capacity (1−τ^ω)A_i·h; net income needs no h^γ term (γ applies only to teaching).
function W_nonteach(XO, zi, l, g, e, Λ, p::Params, gr)
    i0  = gr.nonteach[1]
    sO  = s_closed(i0, p)
    τeO = p.τe[i0, g]                          # common τ^e across non-teaching occupations
    cap = gr.Q[l] * gr.z[zi]^p.α * XO * sO^p.φ * e^p.η   # income capacity (1−τ^ω)A_i·h
    W   = Vector{Float64}(undef, gr.L)
    for lp in 1:gr.L
        Y = (1 - gr.t[lp]) * cap - (1 + τeO) * e          # post-tax net consumption
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
# wage runs through h^γ) and mc ≡ (1+τe)·e, per work-location l':
#     w_{l'}  ≡ dW_{l'}/du = μ [ a − (1−a) mc / Y_{l'} ]
#     w'_{l'} ≡ dw_{l'}/du = −(μa − w_{l'})(μ − w_{l'}) / μ
#     g  ≡ dV̄/du   = Σ_{l'} π_{l'} w_{l'}
#     g' ≡ d²V̄/du² = Var_π(w)/σν + Σ_{l'} π_{l'} w'_{l'}
# where π is the logit location weight (logsum_probs).  `income(e)` returns
# (W, Y, mc): the per-location value vector, the net-consumption vector, and the
# common marginal schooling cost.  The analytic w uses the smooth_log slope 1/Y,
# valid only where the barrier is slack (Y > δ); if any Y ≤ δ, or a Newton step
# leaves the box, or the objective is locally non-concave, we fall back to the
# robust bracketed search.  At a sensible solution the barrier is slack, so the
# fallback essentially never triggers once warm-started.
function solve_e_bracket(income, l, p::Params; lo = log(1e-8), hi = log(1e4))
    f(le) = -logsum_probs(income(exp(le))[1], l, p)[1]
    res = optimize(f, lo, hi)
    return exp(Optim.minimizer(res)), -Optim.minimum(res)        # e, V̄
end

function solve_e_newton(income, a, l, p::Params, u0;
                        lo = log(1e-8), hi = log(1e4), δ = 1e-6,
                        maxit = 30, gtol = 1e-10)
    u = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    for _ in 1:maxit
        e = exp(u)
        W, Y, mc = income(e)
        Vbar, π = logsum_probs(W, l, p)
        any(≤(δ), Y) && return solve_e_bracket(income, l, p; lo, hi)   # barrier active → 1/Y invalid
        w    = @. p.μ * (a - (1 - a) * mc / Y)          # dW_{l'}/du
        g    = dot(π, w)
        varw = dot(π, w .^ 2) - g^2
        wp   = @. -(p.μ * a - w) * (p.μ - w) / p.μ      # dw_{l'}/du
        gp   = varw / p.σν + dot(π, wp)
        abs(g) < gtol && return e, Vbar
        gp ≥ -1e-14 && return solve_e_bracket(income, l, p; lo, hi)    # not concave here
        unew = u - g / gp
        (isfinite(unew) && lo ≤ unew ≤ hi) || return solve_e_bracket(income, l, p; lo, hi)
        u = unew
    end
    return solve_e_bracket(income, l, p; lo, hi)         # no convergence → robust fallback
end

# Teaching branch: h = base·e^η, wage_{l'} = κ_{l'} h^γ, so net income runs
# through e^{ηγ}; effective goods exponent a = ηγ.  Seeded at u0 = log(e_prev).
function solve_e_teach(ϵT, zi, l, g, Λ, p::Params, gr, u0)
    sT   = s_closed(p.T, p)
    base = gr.Q[l] * (gr.z[zi] * ϵT)^p.α * sT^p.φ
    τeT  = p.τe[p.T, g]
    coef = 1 - p.τω[p.T, g]
    function income(e)
        h  = base * e^p.η
        mc = (1 + τeT) * e
        W  = Vector{Float64}(undef, gr.L)
        Y  = Vector{Float64}(undef, gr.L)
        for lp in 1:gr.L
            y = (1 - gr.t[lp]) * coef * p.κ[lp] * h^p.γ - mc
            Y[lp] = y
            W[lp] = p.μ * smooth_log(y) + Λ[zi, lp]
        end
        return W, Y, mc
    end
    return solve_e_newton(income, p.η * p.γ, l, p, u0)
end

# Non-teaching branch: income capacity cap = K·e^η; effective goods exponent a = η.
function solve_e_nonteach(XO, zi, l, g, Λ, p::Params, gr, u0)
    i0  = gr.nonteach[1]
    sO  = s_closed(i0, p)
    τeO = p.τe[i0, g]
    K   = gr.Q[l] * gr.z[zi]^p.α * XO * sO^p.φ
    function income(e)
        cap = K * e^p.η
        mc  = (1 + τeO) * e
        W   = Vector{Float64}(undef, gr.L)
        Y   = Vector{Float64}(undef, gr.L)
        for lp in 1:gr.L
            y = (1 - gr.t[lp]) * cap - mc
            Y[lp] = y
            W[lp] = p.μ * smooth_log(y) + Λ[zi, lp]
        end
        return W, Y, mc
    end
    return solve_e_newton(income, p.η, l, p, u0)
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
# log-human-capital integrates the child's optimal occupation choice:
#   E[log h'] = ∫ f_T(ε_T) log h_T(ε_T) teach_wt(ε_T) dε_T
#             + ∫ f_{X_O*}(x) E[log h_O | x] nonteach_wt(x) dx,
# where log h_O = log[Q_l z^α s_O^φ e^η] + log X_O* − E[log Θ_{i*} | X_O*].
# Then Λ[zi, l'] = ½ Σ_{g'} Σ_{z'} Π_z[zi,z'] λ E[log h' | l', g', z'].
# -----------------------------------------------------------------------------
function child_Elogh(l, g, zi, hh, p, gr)
    cm  = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
    sO  = s_closed(gr.nonteach[1], p)
    # teaching: log h_T directly available on the ε_T grid.
    loghTvec = log.(hh.hT[g, l, zi, :])
    # non-teaching: reconstruct log h_O on the X_O grid.
    base = log(gr.Q[l]) + p.α * log(gr.z[zi]) + p.φ * log(sO)
    loghOvec = [base + p.η * log(hh.eO[g, l, zi, k]) +
                log(gr.XOgrid[k, g]) - ElogΘ(gr, g, gr.XOgrid[k, g])
                for k in 1:gr.nXO]
    IT, IO = integrate_choice(loghTvec, loghOvec, cm, gr)
    return IT + IO
end

function compute_lambda(hh, p, gr)
    (; Nz, L, Πz) = gr
    Elogh = zeros(L, 2, Nz)                    # [l', g', z']
    for lp in 1:L, g in 1:2, zpi in 1:Nz
        Elogh[lp, g, zpi] = child_Elogh(lp, g, zpi, hh, p, gr)
    end
    Λ = zeros(Nz, L)
    for zi in 1:Nz, lp in 1:L
        acc = 0.0
        for g in 1:2, zpi in 1:Nz
            acc += 0.5 * Πz[zi, zpi] * p.λ * Elogh[lp, g, zpi]
        end
        Λ[zi, lp] = acc
    end
    return Λ
end

# -----------------------------------------------------------------------------
# 7. Household block.  Iterate the intergenerational fixed point
#    V → h' → Λ → V.  s is closed form; only the goods margin e and the value
#    are refined each sweep.  No occupation tremble: the atomless integration
#    keeps Λ (hence V) continuous in V, so the map contracts.
#
#    Λ0 warm-starts the altruism term (e.g. from the previous GE iteration's
#    household solution): the sweep count is governed by how far Λ starts from
#    its fixed point, so this saves most of the household time in the GE loop.
#
#    eT0/eO0 likewise warm-start the goods policy.  Each per-node e-solve is
#    seeded at the previous sweep's e (and, on the first sweep, at eT0/eO0 from
#    the previous GE iteration when supplied), so the Newton e-solver converges in
#    1–2 steps instead of re-searching from scratch every sweep.
# -----------------------------------------------------------------------------
function solve_household(p::Params, gr; tol = 1e-6, maxit = 500,
                          verbose = false, print_every = 25,
                          Λ0 = nothing, eT0 = nothing, eO0 = nothing)
    (; Nz, L, nϵT, nXO) = gr
    sT = s_closed(p.T, p)
    sO = s_closed(gr.nonteach[1], p)
    # Policy/value arrays indexed [gender, birth-location, ability node, shock node].
    # Teaching side lives on the ϵT grid, non-teaching side on the X_O* grid.
    eT = zeros(2, L, Nz, nϵT); WT = zeros(2, L, Nz, nϵT); hT = zeros(2, L, Nz, nϵT)  # goods e, value W_T, human capital h_T
    eO = zeros(2, L, Nz, nXO); WO = zeros(2, L, Nz, nXO)                              # goods e, value W_O (non-teaching)
    # Seed the goods policy from a previous GE iteration when its shape matches;
    # otherwise the first solve falls back to the bracketed search (u0 non-finite).
    (eT0 !== nothing && size(eT0) == size(eT)) && (eT .= eT0)
    (eO0 !== nothing && size(eO0) == size(eO)) && (eO .= eO0)
    Λ  = (Λ0 !== nothing && size(Λ0) == (Nz, L)) ? copy(Λ0) : zeros(Nz, L)            # altruism term Λ[zi, l'] (warm-started if given)
    hh = (; eT, WT, hT, eO, WO, Λ)             # the household solution bundle, updated each sweep
    t0 = time()
    local err
    for it in 1:maxit
        WT0 = copy(WT); WO0 = copy(WO)     # previous-sweep values, for the convergence check
        # Refresh every (gender, location, ability, shock) policy given the current Λ.
        # u0 = log(previous e) warm-starts each Newton solve; a zero (unset) entry
        # gives a non-finite seed, which the solver treats as a cold bracketed solve.
        for g in 1:2, l in 1:L, zi in 1:Nz
            for k in 1:nϵT
                e, V̄ = solve_e_teach(gr.ϵTgrid[k], zi, l, g, Λ, p, gr, log(eT[g, l, zi, k]))
                eT[g, l, zi, k] = e
                WT[g, l, zi, k] = log(1 - sT) + V̄   # add the log time-endowment term log(1−s) to the location value
                hT[g, l, zi, k] = gr.Q[l] * (gr.z[zi] * gr.ϵTgrid[k])^p.α * sT^p.φ * e^p.η
            end
            for k in 1:nXO
                e, V̄ = solve_e_nonteach(gr.XOgrid[k, g], zi, l, g, Λ, p, gr, log(eO[g, l, zi, k]))
                eO[g, l, zi, k] = e
                WO[g, l, zi, k] = log(1 - sO) + V̄
            end
        end
        Λ = compute_lambda(hh, p, gr)      # re-derive the altruism term from the updated policies
        hh = (; eT, WT, hT, eO, WO, Λ)
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
            W, _ = W_teach(gr.ϵTgrid[k], zi, l, g, hh.eT[g, l, zi, k], hh.Λ, p, gr)
            _, πi = logsum_probs(W, l, p)
            @views πT[g, l, zi, k, :] .= πi
        end
        for k in 1:nXO
            W, _ = W_nonteach(gr.XOgrid[k, g], zi, l, g, hh.eO[g, l, zi, k], hh.Λ, p, gr)
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
    sO   = s_closed(gr.nonteach[1], p)
    HT = zeros(L); Inc = zeros(L); Wbill = zeros(L)   # teacher aggregator, taxable income, teacher wage bill
    for l0 in 1:L, zi in 1:Nz
        cell_mass = 0.5 * Φ[zi, l0]                         # mass in this (birth-loc, ability) cell, ½ per gender
        cell_mass == 0.0 && continue
        for g in 1:2
            cm   = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
            hTv  = hh.hT[g, l0, zi, :]          # teacher human capital across the ε_T grid
            capO = [gr.Q[l0] * gr.z[zi]^p.α * gr.XOgrid[k, g] * sO^p.φ *
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

function solve_ge(p::Params = Params();
                  Nz = 5, nϵT = 64, nXO = 64, damping = 0.3, tol = 1e-5,
                  maxit = 1000, hh_tol = 1e-6, hh_maxit = 1000, verbose = true,
                  init = nothing)
    I = length(p.A); L = length(p.B)
    @assert length(p.κ) == L && size(p.τmove) == (L, L) &&
            size(p.τω) == (I, 2) && size(p.τe) == (I, 2) && 1 <= p.T <= I "Params dims inconsistent"

    if verbose
        @printf("  [ε=0, continuous]  I=%d, L=%d, Nz=%d  shock grids: nϵT=%d, nXO=%d\n",
                I, L, Nz, nϵT, nXO)
        @printf("  Closed-form s: ")
        for i in 1:I
            @printf("%s s=%.4f   ", i == p.T ? "T" : "occ$i", s_closed(i, p))
        end
        println(); flush(stdout)
    end

    HT = fill(1.0, L); M = fill(p.Mtot / L, L); t = fill(0.10, L)
    # Household warm start: Λ from the previous GE iteration (and, optionally,
    # a previously solved equilibrium passed via `init` for continuation runs).
    Λ0 = nothing; eT0 = nothing; eO0 = nothing
    if init !== nothing
        if length(init.HT) == L
            HT .= init.HT; M .= init.M; t .= init.t
            size(init.hh.Λ) == (Nz, L) && (Λ0 = copy(init.hh.Λ))
            eT0 = copy(init.hh.eT); eO0 = copy(init.hh.eO)   # seed goods policy for continuation
        elseif verbose
            @printf("  ⚠ warm-start ignored: init has L=%d ≠ %d\n", length(init.HT), L)
        end
    end
    local hh, Pi_T, Pi_O, Φ, πbar
    for it in 1:maxit
        # One GE iteration: rebuild grids at the current (H̃_T, M, t), solve the
        # household problem, derive location choices, the stationary distribution,
        # and the implied aggregates.  @elapsed times each stage for the log line.
        gr = build_grids(p; H̃T = HT, M = M, t = t, Nz, nϵT, nXO)
        t_hh  = @elapsed (hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit,
                                               verbose = verbose && it == 1, Λ0, eT0, eO0))
        Λ0 = hh.Λ                              # carry Λ forward to warm-start the next iteration's household solve
        eT0 = hh.eT; eO0 = hh.eO               # carry the goods policy forward too
        t_pi  = @elapsed ((Pi_T, Pi_O) = location_probs(hh, p, gr))
        t_phi = @elapsed ((Φ, πbar) = stationary_phi(Pi_T, Pi_O, hh, p, gr))
        t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, Pi_T, Pi_O, hh, p, gr))

        # Convergence = max relative change in H̃_T and M, and absolute change in t.
        err = max(maximum(abs, (HTn .- HT) ./ HT),
                  maximum(abs, (Mn  .- M ) ./ M ),
                  maximum(abs,  tn  .- t))
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

    gr = build_grids(p; H̃T = HT, M = M, t = t, Nz, nϵT, nXO)
    # Final clean solve at the converged (H̃_T, M, t) so the returned policies,
    # location probabilities, and distribution are all mutually consistent.
    hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit, Λ0, eT0, eO0)
    Pi_T, Pi_O = location_probs(hh, p, gr)
    Φ, πbar = stationary_phi(Pi_T, Pi_O, hh, p, gr)
    # The solution bundle `sol`: household policies hh, location probs Pi_T/Pi_O,
    # spatial distribution Φ, migration matrix πbar, aggregates HT/M/t, grids gr, params p.
    return (; hh, Pi_T, Pi_O, Φ, πbar, HT, M, t, gr, p)
end

# -----------------------------------------------------------------------------
# 12. Diagnostics
# -----------------------------------------------------------------------------

"""
    verify_solution(sol; δ = 1e-6, tclamp = 0.9)

Post-solve audit that the numerical safeguards are SLACK at the solution:
(i) no (state, work-location) cell has consumption Y ≤ δ, i.e. the smooth_log
barrier never binds and every stored value is the true log; (ii) no local tax
rate sits on the clamp [0, tclamp].  Returns true if clean.
"""
function verify_solution(sol; δ = 1e-6, tclamp = 0.9)
    (; hh, t, gr, p) = sol
    sT = s_closed(p.T, p)
    i0 = gr.nonteach[1]
    sO = s_closed(i0, p)
    minY = Inf
    n_binding = 0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        for k in 1:gr.nϵT                          # teaching branch
            e = hh.eT[g, l, zi, k]
            h = gr.Q[l] * (gr.z[zi] * gr.ϵTgrid[k])^p.α * sT^p.φ * e^p.η
            for lp in 1:gr.L
                Y = (1 - gr.t[lp]) * (1 - p.τω[p.T, g]) * p.κ[lp] * h^p.γ -
                    (1 + p.τe[p.T, g]) * e
                minY = min(minY, Y); Y <= δ && (n_binding += 1)
            end
        end
        for k in 1:gr.nXO                          # non-teaching branch
            e   = hh.eO[g, l, zi, k]
            cap = gr.Q[l] * gr.z[zi]^p.α * gr.XOgrid[k, g] * sO^p.φ * e^p.η
            for lp in 1:gr.L
                Y = (1 - gr.t[lp]) * cap - (1 + p.τe[i0, g]) * e
                minY = min(minY, Y); Y <= δ && (n_binding += 1)
            end
        end
    end
    ok = true
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
    ok && @printf("  ✓ feasibility audit passed: min C = %.4e > δ = %.1e; taxes interior t = %s\n",
                  minY, δ, fmtvec(t))
    return ok
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

function report_ge(sol)
    (; hh, Φ, πbar, HT, M, t, gr, p) = sol
    (; z, Nz, Q, Πz, L) = gr

    println("\n===== GE stationary solution (ε = 0, continuous distributions) =====")
    @printf("  s (closed): ")
    for i in 1:gr.I
        @printf("%s=%.4f  ", i == p.T ? "T" : "occ$i", s_closed(i, p))
    end
    println()
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
    println("====================================================================\n")
end

# -----------------------------------------------------------------------------
# 13. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model (ε = 0, continuous) in general equilibrium ...")
    # Solve at the default parameters with a moderately fine shock grid.
    # damping = 0.5 blends each GE update halfway toward the new fixed-point guess;
    # tol = 1e-6 is the convergence threshold on the (H̃_T, M, t) residual.
    sol = solve_ge(Params(); Nz = 3, nϵT = 48, nXO = 48, damping = 0.9,
                   tol = 1e-6, hh_tol = 1e-6, maxit = 500, hh_maxit = 1000,
                   verbose = true)
    report_ge(sol)
    verify_solution(sol)
    return sol
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
