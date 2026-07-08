# =============================================================================
# Spatial occupational-choice model with teacher spillovers + altruism.
# GENERAL-EQUILIBRIUM stationary solution:  L locations, I occupations.
#
# The default parameterization is L = 2 locations, I = 2 occupations
# (Teacher vs. non-teaching), exactly as in the original 2x2 model. The code
# below is written for general I and L: every location is parameterized the
# same way (amenity B, symmetric moving costs τmove, teaching-wage shifter κ),
# and every non-teaching occupation i ≠ T shares the same linear-in-h wage
# function ω_i(h) = A[i] h, but may have its own productivity A[i].
#
# The household problem is itself a fixed point across generations: a parent's
# old-age consumption includes the share ε of her child's education outlay, and
# in a stationary equilibrium the child's policy equals the parent's policy.
#
# Outer fixed point over the aggregate vector (H̃_T, M, t):
#   1. Household block: given (Q,t), iterate the intergenerational policy
#      fixed point  -> S, E, H, V.
#   2. Location-choice probabilities  π_{l'|l}(h; i,g,z).
#   3. Stationary joint mass Φ_l(z): fixed point (dominant left eigenfunction) of the
#      migration ⊗ ability-transition kernel
#      K_{(l₀,z)→(l,z')} = π_z(z'|z) π̄_{l₀→l}(z).
#   4. New aggregates from policies integrated against the endogenous Φ:
#         M_l   = 2∫Φ_l
#         H̃_T,l = ½Σ_g Σ_{l₀} ∫∫ 1_T π_{l|l₀} h^{β/σ} dF_ε Φ_{l₀}
#         t_l   = 𝒲_l / 𝓘_l                      (local balanced budget)
#   5. Damp and repeat until the aggregates stop moving.
#
# key spatial subtlety: because z is the only inherited state
# and parents choose where their child is born, the cross-sectional distribution
# of z is *endogenous and location-specific*: Φ_l(z) ≠ the ergodic G*_z, and
# generally G_l(z) ≠ G_l'(z) for l ≠ l'.  Every aggregate integrates policies
# against Φ_l(z) — NOT the ergodic z law.
#
# Occupation choice with I occupations: each occupation i draws its own
# idiosyncratic shock ϵ_i, so the joint shock state is an I-tuple
# ϵidx = (ϵ_1,...,ϵ_I) (one Nϵ-grid index per occupation), living on
# CartesianIndices((Nϵ,...,Nϵ)) (I copies). The agent picks
# i★ = argmax_i V[i,g,l,z,ϵ_i]; for I = 2 this is exactly the original
# "teacher draws ϵ1, non-teacher draws ϵ2" Roy choice. NOTE: the size of this
# joint-shock grid is Nϵ^I, so larger I needs a coarser Nϵ to stay tractable.
#
# Symbol note:  ε  (\varepsilon) = parental cost share        (p.ε, scalar)
#               ϵ  (\epsilon)    = idiosyncratic ability shock (grid ϵgrid)
#               β  (\beta)       = teacher-HC elasticity in the raw schooling tech;
#                                  enters only the aggregator H̃_T = Σ∫ h^{β/σ}.
# =============================================================================

using LinearAlgebra
using Printf
using Optim

# When an (s,e) trial yields negative consumption we report it.  The GE solver
# silences this during intermediate iterations; only the final solution prints it.
const INFEASIBLE_VERBOSE = Ref(true)

# Smooth log barrier: log(Y) for Y > δ, then a linear extension below δ.
# The linear part has slope 1/δ, so LBFGS always sees a nonzero gradient
# pointing away from infeasibility (unlike log(max(Y, 1e-30)) which is flat
# for any Y ≤ 1e-30 and leaves LBFGS with no information).
@inline smooth_log(Y, δ = 1e-6) = Y > δ ? log(Y) : log(δ) + (Y - δ) / δ

# -----------------------------------------------------------------------------
# 1. Parameters.  Occupation T (default 1) = Teacher; all other occupations
#    1..I share the same linear wage ω_i(h) = A[i] h (A[T] unused).
#    Gender index: 1 = m, 2 = f (fixed at 2). Distortion matrices τω, τe are
#    [occupation, gender]. Location vectors B, κ have length L; τmove is L×L.
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
    ε::Float64 = 0.0                            # parental share of child's education
    # Logit scale for the occupation (Roy) choice.  0.0 = exact hard argmax (the
    # original model); >0 makes the child's occupation a softmax over occupation
    # values, mirroring the Gumbel location choice. A small scale removes the
    # argmax discontinuity that otherwise creates period-2 limit cycles in the
    # household fixed point (see HOUSEHOLD_CONVERGENCE_FINDINGS.md).
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
    # General-equilibrium aggregates. Teacher human capital enters the raw
    # schooling tech as h_T^β; after the efficient class-size reduction (notes
    # A.1) only the aggregator survives,
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
    # This removes the runtime-rank CartesianIndex / abstract-typed prod() in the hot loop.
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
# 3. Stage-2 value of each work location, and the Gumbel log-sum
# -----------------------------------------------------------------------------

# Occupation choice (shared helper).
#
# sflat is the flat shock index ∈ 1:Sϵ; Edecode[sflat, i] is the shock-grid
# index for occupation i at that flat index. Returns (i*, ei*) where ei* is
# the own-shock grid index of the chosen occupation.
# Ties favour the lower-index occupation (occupation 1 = Teacher by default).
@inline function choose_occupation(V, g, l, zi, sflat, Edecode)
    Iocc = size(V, 1)
    i★, best = 1, V[1, g, l, zi, Edecode[sflat, 1]]
    for i in 2:Iocc
        v = V[i, g, l, zi, Edecode[sflat, i]]
        if v > best
            i★, best = i, v
        end
    end
    return i★, Edecode[sflat, i★]
end

# Occupation-choice WEIGHTS over the I occupations at flat shock sflat, written
# into the preallocated buffer `w` (length I). With θ == 0 this is a one-hot at
# the Roy argmax (ties → lower index), i.e. EXACTLY `choose_occupation`, so the
# hard-choice model is the nested θ = 0 case. With θ > 0 it is a softmax with
# temperature θ over the occupation values V[i,g,l,zi,Edecode[sflat,i]] — a logit
# occupation choice that is continuous in V and removes the discontinuity behind
# the household limit cycle. The own-shock index for occupation i is
# Edecode[sflat, i] at every i, so callers index H/E/Pi with that.
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

# Child-side objects entering a parent's consumption: for a child born in `l`
# with gender g, persistent state z'(zi), and joint shock vector ϵidx
# (ϵ_1,...,ϵ_I, one per occupation), the child optimally chooses an occupation;
# we store its education outlay D=(1+τe)e' and its human capital h'. Built from
# the current policy guess (E, H, V).
function child_objects(E, H, V, p::Params, grids)
    (; Nz, L, Sϵ, Edecode, I) = grids
    # Dc = child's GROSS education outlay (1+τe)e';  Hc = child's human capital h'.
    # Flat shock index sflat ∈ 1:Sϵ replaces the I trailing CartesianIndex dims.
    # The occupation is chosen with weights `wocc` (θ = 0 ⇒ hard argmax).
    Dc = zeros(L, 2, Nz, Sϵ)
    Hc = zeros(L, 2, Nz, Sϵ)
    wocc = Vector{Float64}(undef, I)
    for l in 1:L, g in 1:2, zi in 1:Nz, sflat in 1:Sϵ
        occupation_weights!(wocc, V, g, l, zi, sflat, Edecode, p.occ_θ)
        for i in 1:I
            wocc[i] == 0.0 && continue
            ei = Edecode[sflat, i]
            Dc[l, g, zi, sflat] += wocc[i] * (1 + p.τe[i, g]) * E[i, g, l, zi, ei]
            Hc[l, g, zi, sflat] += wocc[i] * H[i, g, l, zi, ei]
        end
    end
    return Dc, Hc
end

# Precompute the (s,e)-independent altruism term Λ[zi, lp] once per HH sweep.
# Λ[zi, lp] = Σ_{gp, zpi, sflat} 0.5·Πz[zi,zpi]·Wflat[sflat]·λ·log(Hc[lp,gp,zpi,sflat])
# When ε=0 the full EV collapses to μ·log(Y) + Λ[zi,lp], eliminating the Sϵ inner loop.
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

# Expected value W_l' of working in each location, given the parent's own
# (occupation i, gender g, birth loc l, ability zi, shock ei) and choice (s,e).
# Λ[zi, lp] = precomputed altruism integral (child-HC term, independent of (s,e)).
# Returns the length-L vector W and the parent's human capital h.
function location_values(i, g, l, zi, ei, s, e, Dc, Λ, p::Params, grids)
    (; z, ϵgrid, Πz, Nz, Q, t, L, Sϵ, Wflat) = grids
    h = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
    W = Vector{Float64}(undef, L)
    if p.ε == 0.0
        # Fast path: C = Y independent of child shocks; weights sum to 1 so
        # EV = μ·log(Y) + Λ[zi, lp].  Eliminates the entire Sϵ inner loop.
        # log(max(Y, 1e-30)) keeps the objective smooth across the Y=0 boundary so
        # LBFGS's line search never encounters a hard discontinuity.
        for lp in 1:L
            wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
            Y = (1 - t[lp]) * (1 - p.τω[i, g]) * wage - (1 + p.τe[i, g]) * e
            W[lp] = p.μ * smooth_log(Y) + Λ[zi, lp]
        end
    else
        for lp in 1:L
            wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
            Y = (1 - t[lp]) * (1 - p.τω[i, g]) * wage - (1 - p.ε) * (1 + p.τe[i, g]) * e
            EV = Λ[zi, lp]
            for gp in 1:2, zpi in 1:Nz
                pz = Πz[zi, zpi]
                pz == 0.0 && continue
                for sflat in 1:Sϵ
                    C = Y - p.ε * Dc[lp, gp, zpi, sflat]
                    EV += pz * Wflat[sflat] * 0.5 * p.μ * smooth_log(C)
                end
            end
            W[lp] = EV
        end
    end
    return W, h
end

# Gumbel log-sum over locations (folds in amenity B and moving cost τ) and the
# associated choice probabilities π_{l'|l}.
#
# Given the stage-2 values W[lp] (the expected utility of working in lp, already
# integrated over child shocks), the location-specific amenities B[lp], and the
# utility moving costs τmove[l, lp], the net payoff of choosing destination lp is
#     u(lp) = W[lp] + B[lp] − τmove[l, lp] + σν · ν_{lp},
# where ν_{lp} ~ i.i.d. Gumbel(0,1). The standard Gumbel extreme-value result gives
# the expected maximum (= expected utility from optimal location choice):
#     EV = σν · log Σ_{lp} exp( u(lp) / σν )
# and the choice probabilities (softmax):
#     π_{lp|l} = exp( u(lp)/σν ) / Σ_{lp'} exp( u(lp')/σν ).
# A larger σν increases both location-choice dispersion and the option value of moving.
function logsum_probs(W, l, p::Params)
    L = length(W)
    # Scaled net payoff of each work location lp: value + amenity B − moving cost τ,
    # divided by the Gumbel scale σν. (l = the location the agent is choosing FROM.)
    x   = ((W[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in 1:L)
    x   = collect(x)
    m   = maximum(x)                                  # max-subtract for numerical stability
    den = sum(exp.(x .- m))                           # Σ exp(x - m); equals Σ exp(x) / exp(m)
    # σν·(m + log(den)) = σν·log(exp(m)·den) = σν·log(Σ exp(x)), the log-sum-exp EV.
    # The softmax probabilities are exp(x-m)/den, cancelling exp(m) from numerator and denominator.
    return p.σν * (m + log(den)), exp.(x .- m) ./ den
end

# -----------------------------------------------------------------------------
# 4. Inner problem: maximize  log(1-s) + V̄(h(s,e))  over (s,e) for one state.
#    Unconstrained via transforms  s = logistic(x1) ∈ (0,1),  e = exp(x2) > 0.
# -----------------------------------------------------------------------------
function neg_objective(x, i, g, l, zi, ei, Dc, Λ, p::Params, grids)
    s = 1 / (1 + exp(-x[1]))
    e = exp(x[2])
    W, _ = location_values(i, g, l, zi, ei, s, e, Dc, Λ, p, grids)
    all(==(-Inf), W) && return 1e10
    V̄, _ = logsum_probs(W, l, p)
    return -(log(1 - s) + V̄)
end

# One sweep: re-optimize (s,e) at every state given the current child objects.
function update_policy!(S, E, H, V, Dc, Hc, p::Params, grids)
    (; z, ϵgrid, Nz, Nϵ, Q, I, L) = grids
    Λ = precompute_lambda(Hc, p, grids)
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        # Warm start from the current policy, mapped to the unconstrained space
        # (x1 = logit(s), x2 = log(e)) the optimizer works in.
        s0, e0 = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        x0  = [log(s0 / (1 - s0)), log(e0)]
        # Fallback: if the warm-start objective is already terrible (e.g. the GE
        # aggregates shifted so that (s0,e0) now yields Y < 0), reset to a safe
        # default rather than handing LBFGS a zero-gradient starting point.
        if neg_objective(x0, i, g, l, zi, ei, Dc, Λ, p, grids) > 50.0
            x0 = [log(0.4 / 0.6), log(0.05)]
        end
        res = optimize(x -> neg_objective(x, i, g, l, zi, ei, Dc, Λ, p, grids),
                       x0, LBFGS(), Optim.Options(iterations = 200, g_tol = 1e-6))
        # Map the optimizer's solution back to (s,e) and store the updated policy,
        # the implied human capital H, and the value V (= −minimised objective).
        x = Optim.minimizer(res)
        s, e = 1 / (1 + exp(-x[1])), exp(x[2])
        S[i, g, l, zi, ei] = s
        E[i, g, l, zi, ei] = e
        H[i, g, l, zi, ei] = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
        V[i, g, l, zi, ei] = -Optim.minimum(res)
    end
end

# -----------------------------------------------------------------------------
# 5. Location-choice probabilities  π_{l'|l}(h; i,g,z)  at every state.
#
#    Pi[i,g,l,zi,ei,l'] is the probability that an agent in state (i,g,l,z,ϵ)
#    chooses to work in l'.
# -----------------------------------------------------------------------------
function location_probs(S, E, H, V, p::Params, grids)
    (; Nz, Nϵ, I, L) = grids
    Dc, Hc = child_objects(E, H, V, p, grids)
    Λ = precompute_lambda(Hc, p, grids)
    Pi = zeros(I, 2, L, Nz, Nϵ, L)
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        s, e = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        W, _ = location_values(i, g, l, zi, ei, s, e, Dc, Λ, p, grids)
        _, πi = logsum_probs(W, l, p)
        @views Pi[i, g, l, zi, ei, :] .= πi
    end
    return Pi
end

# -----------------------------------------------------------------------------
# 6. Endogenous spatial distribution of ability Φ_l(z).
#
#    π̄_{l₀→l'}(z) = ½ Σ_g ∫ Σ_i 1_{i,g,l₀}(z,ϵ) π_{l'|l₀} dF_ε   is the gender-
#    pooled probability that a parent born in (l₀,z) works/raises a child in l'.
#    The stationary joint mass solves the eigen-problem
#        Φ_l(z') = Σ_{l₀} Σ_z π_z(z'|z) π̄_{l₀→l}(z) Φ_{l₀}(z),
#    i.e. the dominant left eigenvector of the (L·Nz)×(L·Nz) row-stochastic kernel
#        P[(l₀,z),(l,z')] = π̄_{l₀→l}(z) · π_z(z'|z).
#    Mass is conserved; we rescale Φ so Σ_l Σ_z Φ_l = Mtot/2.
#
#    Returns Φ as an (Nz × L) array [z, location] and the migration tensor
#    πbar[l₀, z, l'].
# -----------------------------------------------------------------------------
function stationary_phi(Pi, V, p::Params, grids)
    (; Nz, Πz, L, Sϵ, Wflat, Edecode, I) = grids

    # πbar[l0, z, l'] = gender- and shock-averaged probability that a parent born in
    # (l0, z) ends up raising a child in l'. Integrate the location probabilities Pi
    # over gender (the 0.5), the flat joint shock index (Wflat quadrature weights),
    # and the occupation choice (weights wocc; θ = 0 ⇒ hard argmax).
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

    # Assemble the (L·Nz)×(L·Nz) row-stochastic kernel on the joint state (location, z):
    #   P[(l0,z) -> (l,z')] = πbar[l0,z,l] · Πz[z,z']  (migrate, then ability evolves).
    # idx flattens (location, z-index) into a single row/column index.
    n   = L * Nz
    P   = zeros(n, n)
    idx(l, zi) = (l - 1) * Nz + zi
    for l0 in 1:L, zi in 1:Nz, l in 1:L, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]
    end
    # Stationary distribution = dominant left eigenvector of P, rescaled so the total
    # mass is Mtot/2 (one young agent per parent). Reshaped to Φ[z, location].
    φ  = stationary(P)
    φ .*= p.Mtot / 2
    Φ  = reshape(φ, Nz, L)
    return Φ, πbar
end

# -----------------------------------------------------------------------------
# 7. Aggregates implied by the policies and the endogenous distribution Φ.
# -----------------------------------------------------------------------------
function aggregates(Φ, Pi, H, V, p::Params, grids)
    (; Nz, L, Sϵ, Wflat, Edecode, I) = grids
    # expo = β/σ is the exponent in the teacher-HC aggregator
    #   H̃_{T,l} = Σ_g ∫ h^{β/σ} dF_{T,l,g}
    # (see notes A.1 and header comment). β < σ means the aggregator is concave
    # in teacher HC, so doubling teacher quality raises Q less than proportionally.
    expo  = p.β / p.σ
    HT    = zeros(L)                                  # accumulates H̃_{T,l} = Σ∫ h^{β/σ} for teachers
    Inc   = zeros(L)                                  # taxable labor income Σ∫ (1−τω)·wage dΦ, per work location
    Wbill = zeros(L)                                  # teacher wage bill Σ∫ wage·1_T dΦ, per work location

    # Integrate policies against the ENDOGENOUS distribution Φ (not the ergodic z law).
    # Φ[zi, l0] is the mass of agents born in location l0 with ability z[zi]; it is
    # determined jointly with migration (see stationary_phi). Using Φ rather than the
    # ergodic G* matters because high-amenity locations attract high-ability parents,
    # so the local ability distribution is endogenously skewed (header, §key subtlety).
    #
    # Outer loops sum over birth state (l0, z); inner loops over gender and the flat
    # shock index sflat (concrete Int, avoids abstract CartesianIndex overhead).
    wocc = Vector{Float64}(undef, I)
    for l0 in 1:L, zi in 1:Nz
        wbase = 0.5 * Φ[zi, l0]                       # mass at (l0,z), halved because each gender gets equal weight
        wbase == 0.0 && continue
        for g in 1:2, sflat in 1:Sϵ
            # Occupation choice weights (θ = 0 ⇒ one-hot at the Roy argmax).
            occupation_weights!(wocc, V, g, l0, zi, sflat, Edecode, p.occ_θ)
            # wjk = mass element for this (l0, zi, g, shock) cell; Wflat[sflat] is
            # the joint quadrature weight ∏_i ϵw[Edecode[sflat,i]] for the I-shock draw.
            wjk = wbase * Wflat[sflat]
            for i in 1:I
                wocc[i] == 0.0 && continue
                ei = Edecode[sflat, i]
                h  = H[i, g, l0, zi, ei]
                # Agents born in (l0, z) choose a work location lp with probability Pi[...].
                # Pi already integrates the Gumbel shocks, so this spreads each birth-state
                # mass cell across all potential work destinations.
                for lp in 1:L
                    pij = Pi[i, g, l0, zi, ei, lp]
                    pij == 0.0 && continue
                    # After-distortion wage: teacher earns κ_lp h^γ, non-teacher earns A_i h.
                    # The labor-income wedge τω reduces the income subject to taxation.
                    wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
                    Inc[lp] += wjk * wocc[i] * pij * (1 - p.τω[i, g]) * wage
                    if i == p.T                          # only teachers feed H̃_T and the wage bill
                        Wbill[lp] += wjk * wocc[i] * pij * wage     # gross (pre-tax) teacher wage for budget
                        HT[lp]    += wjk * wocc[i] * pij * h^expo   # teacher-HC aggregator contribution
                    end
                end
            end
        end
    end

    # M_l = 2 × (total mass of parents born in l), because each parent-child pair
    # contributes one young agent and one old agent to the headcount.
    M = [2 * sum(@view Φ[:, l]) for l in 1:L]
    # Local balanced-budget tax rate: t_l = Wbill_l / Inc_l (teacher wages financed
    # entirely by a proportional tax on local labor income). Clamped to [0, 0.9]
    # to keep consumption positive even at extreme calibrations.
    t = [Inc[lp] > 0 ? clamp(Wbill[lp] / Inc[lp], 0.0, 0.9) : 0.0 for lp in 1:L]
    return HT, M, t
end

# -----------------------------------------------------------------------------
# 8. Household block: warm-startable solver used inside the GE loop.
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
    recompute_H!(H, S, E, p, grids)
    t0 = time()
    for iter in 1:maxit
        Dc, Hc = child_objects(E, H, V, p, grids)
        E0, V0 = copy(E), copy(V)
        update_policy!(S, E, H, V, Dc, Hc, p, grids)
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
# 9. Outer general-equilibrium fixed point over (H̃_T, M, t).
# -----------------------------------------------------------------------------
"Format a vector as \"(x1, x2, ...)\" with 4 decimals, for variable-length printing."
fmtvec(v) = "(" * join((@sprintf("%.4f", x) for x in v), ", ") * ")"

function solve_ge(p::Params = Params();
                  Nz = 5, Nϵ = 7, damping = 0.5, tol = 1e-5, maxit = 1000,
                  hh_tol = 1e-6, hh_maxit = 1000, verbose = true)
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
        @printf("  Grid: I=%d, L=%d, Nz=%d, Nϵ=%d  →  joint-shock grid Nϵ^I=%d,  policy states=%d\n",
                I, L, Nz, Nϵ, joint_shock, policy_states)
        joint_shock > 100 &&
            @printf("  ⚠  Nϵ^I=%d is large (was Nϵ=%d pre-generalisation). Consider reducing Nϵ to 5–7.\n",
                    joint_shock, Nϵ)
        flush(stdout)
    end

    # Initial guess for the aggregate vector the outer loop solves for.
    HT = fill(1.0, L)                                # teacher-HC aggregator H̃_T
    M  = fill(p.Mtot / L, L)                         # location headcounts (start even)
    t  = fill(0.10, L)                               # local tax rates
    # Policy/value arrays, shared (and warm-started) across GE iterations.
    dims = (I, 2, L, Nz, Nϵ)                         # [occ, gender, loc, z, ϵ]
    S = fill(0.40, dims...)                          # schooling time s
    E = fill(0.10, dims...)                          # education goods e
    H = zeros(dims...)                               # human capital h
    V = zeros(dims...)                               # value function

    local Pi, Φ, πbar
    try
        for it in 1:maxit
            # One GE iteration: rebuild grids at current aggregates, solve the
            # household block, then recompute choice probs, the distribution, and
            # the implied new aggregates (HTn, Mn, tn).
            grids = build_grids(p; H̃T = HT, M = M, t = t, Nz, Nϵ)

            t_hh = @elapsed solve_household!(S, E, H, V, p, grids;
                                             tol = hh_tol, maxit = hh_maxit,
                                             verbose = verbose && it == 1,
                                             print_every = max(1, hh_maxit ÷ 20))
            t_pi  = @elapsed (Pi      = location_probs(S, E, H, V, p, grids))
            t_phi = @elapsed ((Φ, πbar) = stationary_phi(Pi, V, p, grids))
            t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, Pi, H, V, p, grids))

            # Convergence measured on the relative move of (H̃_T, M) and absolute move of t.
            err = max(maximum(abs, (HTn .- HT) ./ HT),
                      maximum(abs, (Mn  .- M ) ./ M ),
                      maximum(abs,  tn  .- t))
            # Dampened update toward the new aggregates (stabilises the fixed point).
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
# 10. Diagnostics
# -----------------------------------------------------------------------------

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
            s += 0.5 * gz * Wflat[sflat] * wocc[p.T]   # probability mass on teaching
        end
    end
    return s
end

function report_ge(sol)
    (; V, Φ, πbar, HT, M, t, grids, p) = sol
    (; z, Nz, Q, Πz, L) = grids

    println("\n========== General-equilibrium stationary solution ==========")
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
    println("=============================================================\n")
end

# -----------------------------------------------------------------------------
# 11. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model in general equilibrium ...")
    sol = solve_ge(Params(); Nz = 3, Nϵ = 11, damping = 0.5, tol = 1e-5)
    report_ge(sol)
    return sol
end


# Only auto-run the headline solve when this file is executed as a script
# (`julia spatial_ge.jl`), NOT when it is `include`d by a test/diagnostic file.
if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
