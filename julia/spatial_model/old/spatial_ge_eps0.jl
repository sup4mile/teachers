# =============================================================================
# Spatial occupational-choice model with teacher spillovers + altruism.
# GENERAL-EQUILIBRIUM stationary solution:  L locations, I occupations.
#
#                       *** ε = 0 SPECIALIZATION ***
#
# This is a STANDALONE rewrite of spatial_ge.jl that hard-codes the parental
# cost share ε = 0 (parents finance none of their child's education; each agent
# pays the full cost of her own schooling out of future labor income).  It does
# NOT include or import spatial_ge.jl — everything needed is copied here so the
# file runs on its own (`julia spatial_ge_eps0.jl`).
#
# Why ε = 0 simplifies things (see notes A.6 in spatial_teachers_model.md):
#
#   1. CONSUMPTION DECOUPLES FROM THE CHILD.  With ε = 0 the parental-finance
#      term vanishes, so old-age consumption conditional on the work location l'
#          C_{l'} = (1−t_{l'})(1−τω) ω_i(h) − (1+τe) e
#      is DETERMINISTIC given the parent's own (state, s, e): it no longer
#      depends on the child's realized (z', ϵ', g').  The entire expectation
#      over child shocks in the consumption term disappears.  The child enters
#      ONLY through the additive altruism term
#          Λ_{l'}(z) = ½ Σ_{g'} E_{z'|z} E_{ϵ'} [ λ f(h'(z',ϵ',l',g')) ],
#      which is independent of the parent's (s, e).
#
#   2. CLOSED FORM FOR THE TIME SHARE s.  Because Λ differentiates to zero in
#      the (s,e) FOCs, the telescoping argument in A.6 pins s down in closed
#      form for EVERY state (i,g,l,z,ϵ), regardless of L, the location wedges,
#      or the altruism strength λ:
#          non-teaching:  s_O* = μφ / (μφ + 1 − η)
#          teaching:      s_T* = μφγ / ((μφ − η)γ + 1).
#      So the whole policy array S collapses to the I numbers {s_i*}.
#
#   3. THE INNER PROBLEM IS 1-D.  With s known analytically, only the goods
#      margin e remains.  e is a scalar fixed point coupled to the location
#      probabilities π_{l'|l} (it does NOT generally have a closed form — A.6),
#      so we solve a single univariate optimization for e at each state instead
#      of the 2-D (s,e) optimization in the general-ε code.
#
# Everything spatial — the endogenous ability distribution Φ_l(z), the migration
# matrix, the teacher-HC aggregator H̃_T, the local balanced budgets t_l — is
# UNCHANGED by ε = 0 and is still solved numerically exactly as in §10.
#
# Outer fixed point over the aggregate vector (H̃_T, M, t):
#   1. Household block: given (Q,t), iterate the intergenerational policy
#      fixed point  -> S (closed form), E, H, V.  The only intergenerational
#      coupling left is via Λ -> π -> e (and the child's own e feeding h').
#   2. Location-choice probabilities  π_{l'|l}(h; i,g,z).
#   3. Stationary joint mass Φ_l(z): dominant left eigenfunction of
#      K_{(l₀,z)→(l,z')} = π_z(z'|z) π̄_{l₀→l}(z).
#   4. New aggregates from policies integrated against the endogenous Φ:
#         M_l   = 2∫Φ_l
#         H̃_T,l = ½Σ_g Σ_{l₀} ∫∫ 1_T π_{l|l₀} h^{β/σ} dF_ε Φ_{l₀}
#         t_l   = 𝒲_l / 𝓘_l                      (local balanced budget)
#   5. Damp and repeat until the aggregates stop moving.
#
# Symbol note:  ε  (\varepsilon) = parental cost share — FIXED AT 0 HERE.
#               ϵ  (\epsilon)    = idiosyncratic ability shock (grid ϵgrid)
#               β  (\beta)       = teacher-HC elasticity; enters only H̃_T.
# =============================================================================

using LinearAlgebra
using Printf
using Optim

# When an (s,e) trial yields negative consumption we report it.  The GE solver
# silences this during intermediate iterations; only the final solution prints it.
const INFEASIBLE_VERBOSE = Ref(true)

# Smooth log barrier: log(Y) for Y > δ, then a linear extension below δ.
# The linear part has slope 1/δ, so the optimizer always sees a nonzero gradient
# pointing away from infeasibility (unlike log(max(Y, 1e-30)) which is flat
# for any Y ≤ 1e-30 and leaves the line search with no information).
@inline smooth_log(Y, δ = 1e-6) = Y > δ ? log(Y) : log(δ) + (Y - δ) / δ

# -----------------------------------------------------------------------------
# 1. Parameters.  Occupation T (default 1) = Teacher; all other occupations
#    1..I share the same linear wage ω_i(h) = A[i] h (A[T] unused).
#    Gender index: 1 = m, 2 = f (fixed at 2). Distortion matrices τω, τe are
#    [occupation, gender]. Location vectors B, κ have length L; τmove is L×L.
#
#    NOTE: there is no ε field — this file is the ε = 0 specialization, and the
#    parental cost share is hard-coded to zero throughout.
# -----------------------------------------------------------------------------
Base.@kwdef struct Params
    T::Int = 1                                   # index of the teaching occupation
    A::Vector{Float64} = [NaN, 1.5]              # length I; non-teaching productivity (A[T] unused)
    # Reduced-form schooling tech: h = Q_l (zϵ)^α s^φ e^η ,  Q_l = (2 H̃_T/M)^σ
    α::Float64 = 0.30                            # ability elasticity
    φ::Float64 = 0.40                            # time-investment elasticity
    η::Float64 = 0.20                            # goods-investment elasticity
    σ::Float64 = 0.25                            # teacher-spillover curvature (in Q)
    # Teaching wage:  ω_T,l'(h) = κ_l' h^γ ;  non-teaching wage = A_i h
    γ::Float64 = 0.80
    κ::Vector{Float64} = [0.85, 1.0]             # length L
    # Preferences
    μ::Float64 = 1.00                            # weight on log consumption
    λ::Float64 = 0.50                            # altruism strength, f(h') = log h'
    # Logit scale for the occupation (Roy) choice.  0.0 = exact hard argmax (the
    # original model); >0 makes the child's occupation a softmax over occupation
    # values, removing the argmax discontinuity behind the household limit cycle.
    occ_θ::Float64 = 0.05
    # Locations: amenity, Gumbel scale, utility moving-cost matrix τ[l,l']
    B::Vector{Float64}      = [0.0, 0.2]        # length L
    σν::Float64             = 0.25
    τmove::Matrix{Float64}  = [0.0 0.30; 0.30 0.0]  # L×L, symmetric, zero diagonal
    # Gender-/occupation-specific distortions ([occupation, gender], I×2)
    τω::Matrix{Float64} = [0.0 0.00; 0.0 0.1]   # labor-income wedge  (1-τω)
    τe::Matrix{Float64} = [0.0 0.00; 0.0 0.00]   # education barrier   (1+τe)
    # Ability processes:  log z' = ρz log z + σξ ξ ;  log ϵ ~ N(-σϵ²/2, σϵ²)
    ρz::Float64 = 0.70
    σξ::Float64 = 0.20
    σϵ::Float64 = 0.30
    # General-equilibrium aggregates.  After the efficient class-size reduction
    # (notes A.1) only the aggregator survives,
    #     H̃_T,l = Σ_g ∫ h^{β/σ} dF_{T,l,g},     Q_l = (2 H̃_T,l / M_l)^σ.
    # The H̃_T-elasticity of the GE map is ≈ β, so any β < 1 contracts.
    β::Float64 = 0.10
    # Total measure of agents; the split Σ_l M_l = Mtot is endogenous.
    Mtot::Float64 = 2.0
end

"""
    Params(I, L; T=1, A_other=1.5, κ_base=0.85, κ_other=1.0,
           B_base=0.0, B_other=0.15, τmove_off=0.30, τω_other=0.1, kwargs...)

Build a Params for `I` occupations and `L` locations, with every location
parameterized the same way (location 1 is the "base" location; locations
2..L all share `B_other`, `κ_other`, and pairwise moving cost `τmove_off`),
and every non-teaching occupation sharing productivity `A_other`. With
I = L = 2 (the default) this reproduces `Params()` exactly. Any other field
of `Params` (e.g. `β`, `σν`, ...) can be overridden via `kwargs`.
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
# 1b. Closed-form time share s (the ε = 0 payoff).
#
# From notes A.6: with ε = 0 the (s,e) FOCs telescope and pin s in closed form
# for EVERY (i,g,l,z,ϵ).  These are the SAME constants as the single-location,
# no-altruism limit (A.4), but here they hold for the full spatial model with
# any L, asymmetric (t_{l'}, κ_{l'}, B_{l'}, τ), and any λ.
#     non-teaching:  s_O* = μφ / (μφ + 1 − η)
#     teaching:      s_T* = μφγ / ((μφ − η)γ + 1)
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

"Rouwenhorst discretization of  log z' = ρ log z + σ ξ, ξ~N(0,1). Returns z (levels) and Π."
function rouwenhorst(N::Int, ρ::Float64, σ::Float64)
    p = (1 + ρ) / 2
    Π = [p 1-p; 1-p p]
    for n in 3:N
        Πp = Π
        Π  = zeros(n, n)
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

"Gauss–Hermite nodes/weights for ∫ f(x) e^{-x²} dx."
function gauss_hermite(n::Int)
    β = sqrt.((1:n-1) ./ 2)
    F = eigen(SymTridiagonal(zeros(n), β))
    return F.values, sqrt(π) .* (F.vectors[1, :] .^ 2)
end

"Discretize ϵ with log ϵ ~ N(-σ²/2, σ²) so that E[ϵ] = 1. Returns (nodes, probs)."
function lognormal_grid(n::Int, σ::Float64)
    x, w = gauss_hermite(n)
    ϵ = exp.(-σ^2 / 2 .+ sqrt(2) * σ .* x)
    p = w ./ sum(w)
    return ϵ, p
end

"Stationary distribution of a row-stochastic matrix Π via power iteration."
function stationary(Π; tol = 1e-10, max_iter = 10000)
    π = fill(1 / size(Π, 1), size(Π, 1))
    for iter in 1:max_iter
        π_new = vec(π' * Π)
        if maximum(abs.(π_new .- π)) < tol
            return π_new ./ sum(π_new)
        end
        π = π_new
    end
    return π ./ sum(π)
end

"Bundle grids + given aggregates (Q, t) into a NamedTuple. I = #occupations,
L = #locations are read off the lengths of p.A and p.B."
function build_grids(p::Params; H̃T, M, t, Nz = 3, Nϵ = 5)
    I = length(p.A)
    L = length(p.B)
    z, Πz     = rouwenhorst(Nz, p.ρz, p.σξ)
    ϵgrid, ϵw = lognormal_grid(Nϵ, p.σϵ)
    Q = (2 .* H̃T ./ M) .^ p.σ
    # Flatten the I-dimensional shock grid to a single concrete index sflat ∈ 1:Sϵ.
    # Edecode[sflat, i] = shock-grid index for occupation i at flat index sflat.
    # Wflat[sflat]      = ∏_i ϵw[Edecode[sflat,i]], the joint quadrature weight.
    Sϵ      = Nϵ^I
    Edecode = Matrix{Int}(undef, Sϵ, I)
    Wflat   = Vector{Float64}(undef, Sϵ)
    for (sflat, ci) in enumerate(CartesianIndices(ntuple(_ -> Nϵ, I)))
        w = 1.0
        for m in 1:I
            Edecode[sflat, m] = ci[m]
            w *= ϵw[ci[m]]
        end
        Wflat[sflat] = w
    end
    return (; z, Πz, ϵgrid, ϵw, Nz, Nϵ, Q, t, I, L, Sϵ, Edecode, Wflat)
end

# -----------------------------------------------------------------------------
# 3. Occupation choice and the child-side altruism object
# -----------------------------------------------------------------------------

# Occupation-choice WEIGHTS over the I occupations at flat shock sflat, written
# into the preallocated buffer `w` (length I). θ == 0 ⇒ one-hot at the Roy argmax
# (ties → lower index); θ > 0 ⇒ softmax with temperature θ over the occupation
# values V[i,g,l,zi,Edecode[sflat,i]], removing the argmax discontinuity.
@inline function occupation_weights!(w, V, g, l, zi, sflat, Edecode, θ)
    Iocc = size(V, 1)
    if θ == 0.0
        i★, best = 1, V[1, g, l, zi, Edecode[sflat, 1]]
        for i in 2:Iocc
            v = V[i, g, l, zi, Edecode[sflat, i]]
            v > best && ((i★, best) = (i, v))
        end
        @inbounds for i in 1:Iocc
            w[i] = (i == i★) ? 1.0 : 0.0
        end
    else
        m = -Inf
        for i in 1:Iocc
            v = V[i, g, l, zi, Edecode[sflat, i]]
            v > m && (m = v)
        end
        den = 0.0
        @inbounds for i in 1:Iocc
            w[i] = exp((V[i, g, l, zi, Edecode[sflat, i]] - m) / θ)
            den += w[i]
        end
        @inbounds for i in 1:Iocc
            w[i] /= den
        end
    end
    return w
end

# Child human capital entering a parent's ALTRUISM term.  For a child born in
# `l` with gender g, persistent state z'(zi), and joint shock vector ϵidx, the
# child optimally chooses an occupation; we store its human capital h'.
#
# NOTE vs the general-ε code: with ε = 0 the parent's consumption no longer
# contains the child's education outlay, so we ONLY need Hc (for Λ) — the Dc
# array of child education outlays is gone entirely.
function child_hc(E, H, V, p::Params, grids)
    (; Nz, L, Sϵ, Edecode, I) = grids
    Hc = zeros(L, 2, Nz, Sϵ)
    wocc = Vector{Float64}(undef, I)
    for l in 1:L, g in 1:2, zi in 1:Nz, sflat in 1:Sϵ
        occupation_weights!(wocc, V, g, l, zi, sflat, Edecode, p.occ_θ)
        for i in 1:I
            wocc[i] == 0.0 && continue
            ei = Edecode[sflat, i]
            Hc[l, g, zi, sflat] += wocc[i] * H[i, g, l, zi, ei]
        end
    end
    return Hc
end

# Precompute the (s,e)-independent altruism term Λ[zi, lp] once per HH sweep.
# Λ[zi, lp] = Σ_{gp, zpi, sflat} ½·Πz[zi,zpi]·Wflat[sflat]·λ·log(Hc[lp,gp,zpi,sflat]).
# With ε = 0 the whole stage-2 value is  V_{lp} = μ·log(C_{lp}) + Λ[zi,lp].
function precompute_lambda(Hc, p::Params, grids)
    (; Πz, Nz, Sϵ, Wflat, L) = grids
    Λ = zeros(Nz, L)
    for zi in 1:Nz, lp in 1:L
        acc = 0.0
        for gp in 1:2, zpi in 1:Nz
            pz = Πz[zi, zpi]
            pz == 0.0 && continue
            for sflat in 1:Sϵ
                acc += pz * Wflat[sflat] * 0.5 * p.λ * log(Hc[lp, gp, zpi, sflat])
            end
        end
        Λ[zi, lp] = acc
    end
    return Λ
end

# -----------------------------------------------------------------------------
# 4. Stage-2 value of each work location (ε = 0), and the Gumbel log-sum
# -----------------------------------------------------------------------------

# Expected value W_l' of working in each location, given the parent's own
# (occupation i, gender g, birth loc l, ability zi, shock ei) and choice (s,e).
#
# ε = 0 ⇒ consumption is deterministic in the child:
#     C_{l'} = (1−t_{l'})(1−τω) ω_i(h) − (1+τe) e ,
# and the stage-2 value separates additively into the consumption term and the
# precomputed altruism term Λ[zi, lp] (no expectation over child shocks):
#     W_{l'} = μ·log(C_{l'}) + Λ[zi, lp].
# This is the whole reason the inner problem collapses to 1-D in e.
function location_values(i, g, l, zi, ei, s, e, Λ, p::Params, grids)
    (; z, ϵgrid, Q, t, L) = grids
    h = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
    W = Vector{Float64}(undef, L)
    for lp in 1:L
        wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
        # Full own-education cost (1+τe)e — no (1−ε) factor because ε = 0.
        Y = (1 - t[lp]) * (1 - p.τω[i, g]) * wage - (1 + p.τe[i, g]) * e
        W[lp] = p.μ * smooth_log(Y) + Λ[zi, lp]
    end
    return W, h
end

# Gumbel log-sum over locations (folds in amenity B and moving cost τ) and the
# associated choice probabilities π_{l'|l}.
#     u(lp) = W[lp] + B[lp] − τmove[l, lp] + σν · ν_{lp},   ν_{lp} ~ iid Gumbel
#     EV     = σν · log Σ_{lp} exp( u(lp)/σν )
#     π_{lp|l} = exp( u(lp)/σν ) / Σ_{lp'} exp( u(lp')/σν )
function logsum_probs(W, l, p::Params)
    L = length(W)
    x   = ((W[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in 1:L)
    x   = collect(x)
    m   = maximum(x)                                  # max-subtract for stability
    den = sum(exp.(x .- m))
    return p.σν * (m + log(den)), exp.(x .- m) ./ den
end

# -----------------------------------------------------------------------------
# 5. Inner problem: with s = s_closed(i,p) FIXED, maximize log(1−s) + V̄(h(s,e))
#    over the single scalar e > 0.  log(1−s) is constant in e, so we maximize
#    V̄ over x = log(e) (Brent), then add log(1−s) when storing V.
# -----------------------------------------------------------------------------
function neg_Vbar_e(loge, i, g, l, zi, ei, s, Λ, p::Params, grids)
    e = exp(loge)
    W, _ = location_values(i, g, l, zi, ei, s, e, Λ, p, grids)
    V̄, _ = logsum_probs(W, l, p)
    return -V̄
end

# One sweep: re-solve the optimal e at every state given the current child HC
# (through Λ).  s is the closed-form constant; only e is optimized.
function update_policy!(S, E, H, V, Hc, p::Params, grids;
                        loge_lo = log(1e-8), loge_hi = log(1e4))
    (; z, ϵgrid, Nz, Nϵ, Q, I, L) = grids
    Λ = precompute_lambda(Hc, p, grids)
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        s = s_closed(i, p)
        # 1-D Brent maximization of V̄ over loge = log(e) (derivative-free, robust
        # to the smooth_log kink).  The objective is concave in e for the
        # interior region, so a bracketed Brent search is reliable.
        res = optimize(loge -> neg_Vbar_e(loge, i, g, l, zi, ei, s, Λ, p, grids),
                       loge_lo, loge_hi)
        e   = exp(Optim.minimizer(res))
        V̄   = -Optim.minimum(res)
        S[i, g, l, zi, ei] = s
        E[i, g, l, zi, ei] = e
        H[i, g, l, zi, ei] = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
        # Store the FULL stage-1 value so occupation choice (which compares V
        # across occupations with different s_i) is correct: V = log(1−s) + V̄.
        V[i, g, l, zi, ei] = log(1 - s) + V̄
    end
end

# -----------------------------------------------------------------------------
# 6. Location-choice probabilities  π_{l'|l}(h; i,g,z)  at every state.
# -----------------------------------------------------------------------------
function location_probs(S, E, H, V, p::Params, grids)
    (; Nz, Nϵ, I, L) = grids
    Hc = child_hc(E, H, V, p, grids)
    Λ  = precompute_lambda(Hc, p, grids)
    Pi = zeros(I, 2, L, Nz, Nϵ, L)
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        s, e = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        W, _ = location_values(i, g, l, zi, ei, s, e, Λ, p, grids)
        _, πi = logsum_probs(W, l, p)
        @views Pi[i, g, l, zi, ei, :] .= πi
    end
    return Pi
end

# -----------------------------------------------------------------------------
# 7. Endogenous spatial distribution of ability Φ_l(z).  (Unchanged by ε.)
#
#    π̄_{l₀→l'}(z) = ½ Σ_g ∫ Σ_i 1_{i,g,l₀}(z,ϵ) π_{l'|l₀} dF_ε   is the gender-
#    pooled probability that a parent born in (l₀,z) works/raises a child in l'.
#    Φ_l is the dominant left eigenvector of the row-stochastic kernel
#        P[(l₀,z),(l,z')] = π̄_{l₀→l}(z) · π_z(z'|z), rescaled so Σ Φ = Mtot/2.
#    Returns Φ as (Nz × L) array [z, location] and the migration tensor πbar.
# -----------------------------------------------------------------------------
function stationary_phi(Pi, V, p::Params, grids)
    (; Nz, Πz, L, Sϵ, Wflat, Edecode, I) = grids

    πbar = zeros(L, Nz, L)
    wocc = Vector{Float64}(undef, I)
    for l0 in 1:L, zi in 1:Nz, g in 1:2, sflat in 1:Sϵ
        occupation_weights!(wocc, V, g, l0, zi, sflat, Edecode, p.occ_θ)
        w = 0.5 * Wflat[sflat]
        for i in 1:I
            wocc[i] == 0.0 && continue
            ei = Edecode[sflat, i]
            for lp in 1:L
                πbar[l0, zi, lp] += w * wocc[i] * Pi[i, g, l0, zi, ei, lp]
            end
        end
    end

    n   = L * Nz
    P   = zeros(n, n)
    idx(l, zi) = (l - 1) * Nz + zi
    for l0 in 1:L, zi in 1:Nz, l in 1:L, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]
    end
    φ  = stationary(P)
    φ .*= p.Mtot / 2
    Φ  = reshape(φ, Nz, L)
    return Φ, πbar
end

# -----------------------------------------------------------------------------
# 8. Aggregates implied by the policies and the endogenous distribution Φ.
#    (Unchanged by ε — these are the §8–§9 objects.)
# -----------------------------------------------------------------------------
function aggregates(Φ, Pi, H, V, p::Params, grids)
    (; Nz, L, Sϵ, Wflat, Edecode, I) = grids
    expo  = p.β / p.σ                                 # exponent in H̃_T = Σ∫ h^{β/σ}
    HT    = zeros(L)
    Inc   = zeros(L)
    Wbill = zeros(L)

    wocc = Vector{Float64}(undef, I)
    for l0 in 1:L, zi in 1:Nz
        wbase = 0.5 * Φ[zi, l0]
        wbase == 0.0 && continue
        for g in 1:2, sflat in 1:Sϵ
            occupation_weights!(wocc, V, g, l0, zi, sflat, Edecode, p.occ_θ)
            wjk = wbase * Wflat[sflat]
            for i in 1:I
                wocc[i] == 0.0 && continue
                ei = Edecode[sflat, i]
                h  = H[i, g, l0, zi, ei]
                for lp in 1:L
                    pij = Pi[i, g, l0, zi, ei, lp]
                    pij == 0.0 && continue
                    wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
                    Inc[lp] += wjk * wocc[i] * pij * (1 - p.τω[i, g]) * wage
                    if i == p.T
                        Wbill[lp] += wjk * wocc[i] * pij * wage
                        HT[lp]    += wjk * wocc[i] * pij * h^expo
                    end
                end
            end
        end
    end

    M = [2 * sum(@view Φ[:, l]) for l in 1:L]
    t = [Inc[lp] > 0 ? clamp(Wbill[lp] / Inc[lp], 0.0, 0.9) : 0.0 for lp in 1:L]
    return HT, M, t
end

# -----------------------------------------------------------------------------
# 9. Household block: warm-startable solver used inside the GE loop.
# -----------------------------------------------------------------------------
function recompute_H!(H, S, E, p::Params, grids)
    (; z, ϵgrid, Nz, Nϵ, Q, I, L) = grids
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        H[i, g, l, zi, ei] = Q[l] * (z[zi] * ϵgrid[ei])^p.α *
                             S[i, g, l, zi, ei]^p.φ * E[i, g, l, zi, ei]^p.η
    end
    return H
end

function solve_household!(S, E, H, V, p::Params, grids; tol = 1e-5, maxit = 200,
                           verbose = false, print_every = 10)
    # s is the closed-form constant per occupation — set it ONCE up front; the
    # household iteration then only refines e (and V) through the altruism term.
    for i in 1:grids.I
        @views S[i, :, :, :, :] .= s_closed(i, p)
    end
    recompute_H!(H, S, E, p, grids)
    t0 = time()
    for iter in 1:maxit
        Hc = child_hc(E, H, V, p, grids)
        E0, V0 = copy(E), copy(V)
        update_policy!(S, E, H, V, Hc, p, grids)
        err = max(maximum(abs, E .- E0), maximum(abs, V .- V0))
        if verbose && (iter == 1 || iter % print_every == 0 || err < tol)
            @printf("    HH %4d/%d  err=%.3e  elapsed=%.1fs\n", iter, maxit, err, time() - t0)
            flush(stdout)
        end
        err < tol && break
    end
    return S, E, H, V
end

# -----------------------------------------------------------------------------
# 10. Outer general-equilibrium fixed point over (H̃_T, M, t).
# -----------------------------------------------------------------------------
"Format a vector as \"(x1, x2, ...)\" with 4 decimals, for variable-length printing."
fmtvec(v) = "(" * join((@sprintf("%.4f", x) for x in v), ", ") * ")"

function solve_ge(p::Params = Params();
                  Nz = 5, Nϵ = 7, damping = 0.5, tol = 1e-5, maxit = 1000,
                  hh_tol = 1e-6, hh_maxit = 1000, verbose = true, init = nothing)
    prev_verbose = INFEASIBLE_VERBOSE[]
    INFEASIBLE_VERBOSE[] = false

    I = length(p.A)
    L = length(p.B)
    @assert length(p.κ) == L && size(p.τmove) == (L, L) &&
            size(p.τω) == (I, 2) && size(p.τe) == (I, 2) &&
            1 <= p.T <= I "Params dimensions inconsistent with I=$I occupations, L=$L locations"

    if verbose
        joint_shock = Nϵ^I
        policy_states = I * 2 * L * Nz * Nϵ
        @printf("  [ε=0 model]  Grid: I=%d, L=%d, Nz=%d, Nϵ=%d  →  joint-shock grid Nϵ^I=%d,  policy states=%d\n",
                I, L, Nz, Nϵ, joint_shock, policy_states)
        @printf("  Closed-form s: ")
        for i in 1:I
            @printf("%s s=%.4f   ", i == p.T ? "T" : "occ$i", s_closed(i, p))
        end
        println()
        joint_shock > 100 &&
            @printf("  ⚠  Nϵ^I=%d is large. Consider reducing Nϵ to 5–7.\n", joint_shock)
        flush(stdout)
    end

    HT = fill(1.0, L)
    M  = fill(p.Mtot / L, L)
    t  = fill(0.10, L)
    dims = (I, 2, L, Nz, Nϵ)
    # S is initialized to the closed-form values (solve_household! also enforces
    # this every call, but seeding here keeps the very first H consistent).
    S = zeros(dims...)
    for i in 1:I
        @views S[i, :, :, :, :] .= s_closed(i, p)
    end
    E = fill(0.10, dims...)
    H = zeros(dims...)
    V = zeros(dims...)

    # Optional warm start (continuation): seed the GE aggregates and household
    # policies from a previously solved equilibrium so the fixed-point iteration
    # starts near the answer. Only used when the grid dimensions match; otherwise
    # we fall back to the cold defaults above. S is reseeded from s_closed by
    # solve_household! every call, so warm-starting it is harmless.
    if init !== nothing
        if size(init.E) == dims
            HT .= init.HT
            M  .= init.M
            t  .= init.t
            S  .= init.S
            E  .= init.E
            H  .= init.H
            V  .= init.V
        elseif verbose
            @printf("  ⚠ warm-start ignored: init dims %s ≠ %s\n",
                    size(init.E), dims)
        end
    end

    local Pi, Φ, πbar
    try
        for it in 1:maxit
            grids = build_grids(p; H̃T = HT, M = M, t = t, Nz, Nϵ)

            t_hh = @elapsed solve_household!(S, E, H, V, p, grids;
                                             tol = hh_tol, maxit = hh_maxit,
                                             verbose = verbose && it == 1,
                                             print_every = max(1, hh_maxit ÷ 20))
            t_pi  = @elapsed (Pi      = location_probs(S, E, H, V, p, grids))
            t_phi = @elapsed ((Φ, πbar) = stationary_phi(Pi, V, p, grids))
            t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, Pi, H, V, p, grids))

            err = max(maximum(abs, (HTn .- HT) ./ HT),
                      maximum(abs, (Mn  .- M ) ./ M ),
                      maximum(abs,  tn  .- t))
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

        grids   = build_grids(p; H̃T = HT, M = M, t = t, Nz, Nϵ)
        solve_household!(S, E, H, V, p, grids; tol = hh_tol, maxit = hh_maxit)
        Pi      = location_probs(S, E, H, V, p, grids)
        Φ, πbar = stationary_phi(Pi, V, p, grids)
        return (; S, E, H, V, Pi, Φ, πbar, HT, M, t, grids, p)
    finally
        INFEASIBLE_VERBOSE[] = prev_verbose
    end
end

# -----------------------------------------------------------------------------
# 10b. θ-continuation (annealing) over the occupation-tremble scale.
#
# The tremble occ_θ > 0 solves a PERTURBED model with bias O(occ_θ); occ_θ = 0
# is the true model but a cold start there hits the argmax limit cycle (A.7).
# Resolution: solve at a comfortably smooth θ first, then re-solve at smaller
# and smaller θ, warm-starting each stage from the previous solution via `init`.
# Near the fixed point the near-hard-argmax stages converge in a few sweeps,
# so the final solution carries (almost) no tremble bias without the cycling.
# -----------------------------------------------------------------------------
"Copy of `p` with the occupation-tremble scale replaced by `θ` (all else equal)."
function with_occθ(p::Params, θ::Float64)
    nt = NamedTuple{fieldnames(Params)}(ntuple(i -> getfield(p, i), fieldcount(Params)))
    return Params(; nt..., occ_θ = θ)   # rightmost keyword wins
end

"""
    solve_ge_annealed(p; θs = [0.05, 0.01, 0.001], kwargs...)

Solve the GE fixed point at each occupation-tremble scale in `θs` (decreasing),
warm-starting every stage from the previous stage's solution. `kwargs` are
passed through to `solve_ge`. Returns the final (smallest-θ) solution.
"""
function solve_ge_annealed(p::Params = Params(); θs = [0.05, 0.01, 0.001],
                           verbose = true, kwargs...)
    @assert issorted(θs; rev = true) "θs should decrease (smooth → sharp)"
    sol = nothing
    for θ in θs
        verbose && @printf("―― θ-continuation stage: occ_θ = %g ――\n", θ)
        sol = solve_ge(with_occθ(p, θ); init = sol, verbose, kwargs...)
    end
    return sol
end

# -----------------------------------------------------------------------------
# 11. Diagnostics
# -----------------------------------------------------------------------------

"""
    verify_solution(sol; δ = 1e-6, tclamp = 0.9)

Post-solve audit that the numerical safeguards are SLACK at the solution:
(i) no (state, work-location) cell has consumption Y ≤ δ, i.e. the smooth_log
barrier never binds and every stored value is the true log; (ii) no local tax
rate sits on the clamp [0, tclamp].  Returns true if clean.
"""
function verify_solution(sol; δ = 1e-6, tclamp = 0.9)
    (; E, H, t, grids, p) = sol
    (; Nz, Nϵ, I, L) = grids
    minY = Inf; nbad = 0
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        e, h = E[i, g, l, zi, ei], H[i, g, l, zi, ei]
        for lp in 1:L
            wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
            Y = (1 - grids.t[lp]) * (1 - p.τω[i, g]) * wage - (1 + p.τe[i, g]) * e
            minY = min(minY, Y); Y <= δ && (nbad += 1)
        end
    end
    ok = true
    if nbad > 0
        @printf("  ⚠ smooth_log barrier BINDING at %d (state, l') cells (min C = %.3e ≤ δ = %.1e): values there are barrier-distorted\n",
                nbad, minY, δ)
        ok = false
    end
    for lp in 1:L
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

"Gender-pooled teaching share among young born in l, weighted by the *endogenous* G_l(z)."
function teaching_share_endo(V, l, Φ, grids, p::Params)
    (; Nz, Sϵ, Wflat, Edecode, I) = grids
    mass = sum(@view Φ[:, l])
    s = 0.0
    wocc = Vector{Float64}(undef, I)
    for zi in 1:Nz
        gz = Φ[zi, l] / mass
        for g in 1:2, sflat in 1:Sϵ
            occupation_weights!(wocc, V, g, l, zi, sflat, Edecode, p.occ_θ)
            s += 0.5 * gz * Wflat[sflat] * wocc[p.T]
        end
    end
    return s
end

function report_ge(sol)
    (; V, Φ, πbar, HT, M, t, grids, p) = sol
    (; z, Nz, Q, Πz, L) = grids

    println("\n========== General-equilibrium stationary solution (ε = 0) ==========")
    @printf("  s (closed): ")
    for i in 1:grids.I
        @printf("%s=%.4f  ", i == p.T ? "T" : "occ$i", s_closed(i, p))
    end
    println()
    @printf("  H̃_T    = %s\n", fmtvec(HT))
    @printf("  M      = %s      [total %.4f]\n", fmtvec(M), sum(M))
    @printf("  t      = %s\n", fmtvec(t))
    @printf("  Q      = %s\n", fmtvec(Q))

    erg = stationary(Πz)
    println("\n  Endogenous ability distribution  G_l(z)   (vs ergodic G*):")
    @printf("    %10s", "z")
    for l in 1:L
        @printf(" %10s", "G_$(l)(z)")
    end
    @printf(" %10s\n", "G*(z)")
    for zi in 1:Nz
        @printf("    %10.4f", z[zi])
        for l in 1:L
            @printf(" %10.4f", Φ[zi, l] / sum(@view Φ[:, l]))
        end
        @printf(" %10.4f\n", erg[zi])
    end
    mean_z(w) = sum(z[zi] * w[zi] for zi in 1:Nz) / sum(w)
    print("    mean z: ")
    for l in 1:L
        @printf(" loc %d = %.4f  ", l, mean_z(@view Φ[:, l]))
    end
    @printf("ergodic = %.4f\n", mean_z(erg))

    println("\n  Migration matrix  Π̄[born→work]  (rows sum to 1):")
    for l0 in 1:L
        mass = sum(@view Φ[:, l0])
        row = [sum(Φ[zi, l0] * πbar[l0, zi, l] for zi in 1:Nz) / mass for l in 1:L]
        print("    born $l0 :")
        for l in 1:L
            @printf("  ->%d %.4f", l, row[l])
        end
        println()
    end

    println("\n  Teaching share among young born in l (endogenous G_l):")
    for l in 1:L
        @printf("    born %d :  %.4f\n", l, teaching_share_endo(V, l, Φ, grids, p))
    end
    println("=====================================================================\n")
end

# -----------------------------------------------------------------------------
# 12. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model (ε = 0) in general equilibrium with θ-continuation ...")
    # A single cold solve at a tiny occ_θ (e.g. 1e-4) is the worst of both
    # worlds: near-hard argmax (limit-cycle risk) AND leftover tremble bias.
    # Anneal instead: smooth solve first, then warm-started sharpening.
    sol = solve_ge_annealed(Params(); θs = [0.05, 0.01, 0.001],
                            Nz = 5, Nϵ = 21, damping = 0.5, tol = 1e-6,
                            hh_tol = 1e-6, maxit = 500, hh_maxit = 1000,
                            verbose = true)
    report_ge(sol)
    verify_solution(sol)
    return sol
end

# Only auto-run the headline solve when this file is executed as a script
# (`julia spatial_ge_eps0.jl`), NOT when it is `include`d by a test/diagnostic file.
if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
