# =============================================================================
# Spatial occupational-choice model with teacher spillovers + altruism.
# GENERAL-EQUILIBRIUM stationary solution:  L = 2 locations, I = 2 occupations.
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
# generally G_1(z) ≠ G_2(z).  Every aggregate integrates policies against Φ_l(z)
# — NOT the ergodic z law.
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

# -----------------------------------------------------------------------------
# 1. Parameters.  Occupation 1 = Teacher (T); occupation 2 = non-teaching.
#    Gender index: 1 = m, 2 = f.  Distortion matrices are [occupation, gender].
# -----------------------------------------------------------------------------
Base.@kwdef struct Params
    T::Int = 1                                   # index of the teaching occupation
    A::Vector{Float64} = [NaN, 1.5]              # non-teaching productivity (A[T] unused)
    # Reduced-form schooling tech: h = Q_l (zϵ)^α s^φ e^η ,  Q_l = (2 H̃_T/M)^σ
    α::Float64 = 0.30                            # ability elasticity
    φ::Float64 = 0.40                            # time-investment elasticity
    η::Float64 = 0.20                            # goods-investment elasticity
    σ::Float64 = 0.25                            # teacher-spillover curvature (in Q)
    # Teaching wage:  ω_T,l'(h) = κ_l' h^γ ;  non-teaching wage = A_i h
    γ::Float64 = 0.80
    κ::Vector{Float64} = [0.85, 1.0]
    # Preferences
    μ::Float64 = 1.00                            # weight on log consumption
    λ::Float64 = 0.50                            # altruism strength, f(h') = log h'
    ε::Float64 = 0.0                            # parental share of child's education
    # Locations: amenity, Gumbel scale, utility moving-cost matrix τ[l,l']
    B::Vector{Float64}      = [0.0, 0.10]
    σν::Float64             = 0.50
    τmove::Matrix{Float64}  = [0.0 0.30; 0.30 0.0]
    # Gender-/occupation-specific distortions
    τω::Matrix{Float64} = [0.0 0.00; 0.0 0.1]   # labor-income wedge  (1-τω)
    τe::Matrix{Float64} = [0.0 0.00; 0.0 0.00]   # education barrier   (1+τe)
    # Ability processes:  log z' = ρz log z + σξ ξ ;  log ϵ ~ N(-σϵ²/2, σϵ²)
    ρz::Float64 = 0.70
    σξ::Float64 = 0.20
    σϵ::Float64 = 0.30
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

"Bundle grids + given aggregates (Q, t) into a NamedTuple."
function build_grids(p::Params; H̃T, M, t, Nz = 3, Nϵ = 5)
    z, Πz     = rouwenhorst(Nz, p.ρz, p.σξ)
    ϵgrid, ϵw = lognormal_grid(Nϵ, p.σϵ)
    Q = (2 .* H̃T ./ M) .^ p.σ
    return (; z, Πz, ϵgrid, ϵw, Nz, Nϵ, Q, t)
end

# -----------------------------------------------------------------------------
# 3. Stage-2 value of each work location, and the Gumbel log-sum
# -----------------------------------------------------------------------------

# Child-side objects entering a parent's consumption: for a child born in `l`
# with gender g, persistent state z'(zi), and shock pair (ϵ1=j, ϵ2=k), the child
# optimally chooses an occupation; we store its education outlay D=(1+τe)e' and
# its human capital h'. Built from the current policy guess (E, H, V).
function child_objects(E, H, V, p::Params, grids)
    (; Nz, Nϵ) = grids
    # Dc = child's GROSS education outlay (1+τe)e';  Hc = child's human capital h'.
    # Both are indexed [birth-loc, gender, z', ϵ1(teacher shock), ϵ2(non-teach shock)].
    Dc = zeros(2, 2, Nz, Nϵ, Nϵ)                 # [birth-loc, gender, z', ϵ1, ϵ2]
    Hc = zeros(2, 2, Nz, Nϵ, Nϵ)
    for l in 1:2, g in 1:2, zi in 1:Nz, j in 1:Nϵ, k in 1:Nϵ
        # Child picks the occupation with the higher value; each occupation draws
        # its OWN shock (teacher uses ϵ1=j, non-teaching uses ϵ2=k). ei = chosen index.
        i′, ei = V[1, g, l, zi, j] ≥ V[2, g, l, zi, k] ? (1, j) : (2, k)
        Dc[l, g, zi, j, k] = (1 + p.τe[i′, g]) * E[i′, g, l, zi, ei]   # what parent co-funds
        Hc[l, g, zi, j, k] = H[i′, g, l, zi, ei]                      # what parent's altruism values
    end
    return Dc, Hc
end

# Expected value W_l' of working in each location, given the parent's own
# (occupation i, gender g, birth loc l, ability zi, shock ei) and choice (s,e).
# Returns the length-2 vector W and the parent's human capital h.
function location_values(i, g, l, zi, ei, s, e, Dc, Hc, p::Params, grids)
    (; z, ϵgrid, ϵw, Πz, Nz, Nϵ, Q, t) = grids
    # Parent's own human capital from the reduced-form schooling tech.
    h = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
    W = fill(-Inf, 2)
    for lp in 1:2                                    # lp = candidate WORK location
        # Wage: teachers earn κ·h^γ; others earn linear A·h.
        wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
        # Y = disposable resources after tax, income wedge, and the parent's OWN
        #     share (1-ε) of this period's education spending e.
        Y = (1 - t[lp]) * (1 - p.τω[i, g]) * wage - (1 - p.ε) * (1 + p.τe[i, g]) * e
        # EV = expected continuation value, integrating over the child's gender gp,
        #      next-period ability z', and the two shock draws (j,k). The 0.5 is the
        #      gender probability; pz·ϵw[j]·ϵw[k] are the remaining quadrature weights.
        EV = 0.0; feasible = true
        for gp in 1:2, zpi in 1:Nz
            pz = Πz[zi, zpi]                          # ability transition prob z -> z'
            pz == 0.0 && continue
            for j in 1:Nϵ, k in 1:Nϵ
                # Old-age consumption: resources Y less the parent's share ε of the
                # child's education outlay. C ≤ 0 means this (s,e) is infeasible.
                C = Y - p.ε * Dc[lp, gp, zpi, j, k]
                if C ≤ 0
                    INFEASIBLE_VERBOSE[] && println("  infeasible (s,e) = ($s, $e) for parent state (i=$i, g=$g, l=$l, z=$zi, ϵ=$ei) and child state (gp=$gp, zpi=$zpi, j=$j, k=$k)")
                    feasible = false; break
                end
                # Flow utility: μ·log(consumption) + λ·log(child's human capital).
                EV += pz * ϵw[j] * ϵw[k] * 0.5 *
                      (p.μ * log(C) + p.λ * log(Hc[lp, gp, zpi, j, k]))
            end
            feasible || break
        end
        W[lp] = feasible ? EV : -Inf                 # value of working in lp
    end
    return W, h
end

# Gumbel log-sum over locations (folds in amenity B and moving cost τ) and the
# associated choice probabilities π_{l'|l}.
function logsum_probs(V, l, p::Params)
    # Scaled net payoff of each work location lp: value + amenity B − moving cost τ,
    # divided by the Gumbel scale σν. (l = the location the agent is choosing FROM.)
    x   = ((V[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in 1:2)
    x   = collect(x)
    m   = maximum(x)                                  # max-subtract for numerical stability
    den = sum(exp.(x .- m))
    # Returns (expected value = σν·log-sum-exp,  softmax choice probabilities π_{l'|l}).
    return p.σν * (m + log(den)), exp.(x .- m) ./ den
end

# -----------------------------------------------------------------------------
# 4. Inner problem: maximize  log(1-s) + V̄(h(s,e))  over (s,e) for one state.
#    Unconstrained via transforms  s = logistic(x1) ∈ (0,1),  e = exp(x2) > 0.
# -----------------------------------------------------------------------------
function neg_objective(x, i, g, l, zi, ei, Dc, Hc, p::Params, grids)
    s = 1 / (1 + exp(-x[1]))
    e = exp(x[2])
    W, _ = location_values(i, g, l, zi, ei, s, e, Dc, Hc, p, grids)
    all(==(-Inf), W) && return 1e10
    V̄, _ = logsum_probs(W, l, p)
    return -(log(1 - s) + V̄)
end

# One sweep: re-optimize (s,e) at every state given the current child objects.
function update_policy!(S, E, H, V, Dc, Hc, p::Params, grids)
    (; z, ϵgrid, Nz, Nϵ, Q) = grids
    for i in 1:2, g in 1:2, l in 1:2, zi in 1:Nz, ei in 1:Nϵ
        # Warm start from the current policy, mapped to the unconstrained space
        # (x1 = logit(s), x2 = log(e)) the optimizer works in.
        s0, e0 = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        x0  = [log(s0 / (1 - s0)), log(e0)]
        res = optimize(x -> neg_objective(x, i, g, l, zi, ei, Dc, Hc, p, grids),
                       x0, LBFGS())
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
# 5. General-equilibrium primitives layered on top of the household Params.
#
#    β and the total measure M do not appear in the PE household problem (which
#    takes H̃_T as given), so they live here rather than in Params.
# -----------------------------------------------------------------------------
Base.@kwdef struct GEParams
    p::Params = Params()
    # Teacher human capital enters the raw schooling tech as h_T^β; after the
    # efficient class-size reduction (notes A.1) only the aggregator survives,
    #     H̃_T,l = Σ_g ∫ h^{β/σ} dF_{T,l,g},     Q_l = (2 H̃_T,l / M_l)^σ.
    # The H̃_T-elasticity of the GE map is ≈ β, so any β < 1 contracts.
    β::Float64 = 0.10
    # Total measure of agents; the split M₁ + M₂ = Mtot is endogenous.
    Mtot::Float64 = 2.0
end

# -----------------------------------------------------------------------------
# 6. Occupation choice (shared helper).
#
#    With I = 2, occupation 1 = Teacher draws on its own shock ϵ₁ (=ϵgrid[j]) and
#    occupation 2 = non-teaching draws on ϵ₂ (=ϵgrid[k]).  Returns (i*, eᵢ*),
#    where eᵢ* is the grid index of the chosen occupation's own shock.
# -----------------------------------------------------------------------------
@inline function choose_occupation(V, g, l, zi, j, k)
    return V[1, g, l, zi, j] ≥ V[2, g, l, zi, k] ? (1, j) : (2, k)
end

# -----------------------------------------------------------------------------
# 7. Location-choice probabilities  π_{l'|l}(h; i,g,z)  at every state.
#
#    Pi[i,g,l,zi,ei,l'] is the probability that an agent in state (i,g,l,z,ϵ)
#    chooses to work in l'.
# -----------------------------------------------------------------------------
function location_probs(S, E, H, V, p::Params, grids)
    (; Nz, Nϵ) = grids
    Dc, Hc = child_objects(E, H, V, p, grids)
    Pi = zeros(2, 2, 2, Nz, Nϵ, 2)
    for i in 1:2, g in 1:2, l in 1:2, zi in 1:Nz, ei in 1:Nϵ
        s, e = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        W, _ = location_values(i, g, l, zi, ei, s, e, Dc, Hc, p, grids)
        _, πi = logsum_probs(W, l, p)
        @views Pi[i, g, l, zi, ei, :] .= πi
    end
    return Pi
end

# -----------------------------------------------------------------------------
# 8. Endogenous spatial distribution of ability Φ_l(z).
#
#    π̄_{l₀→l'}(z) = ½ Σ_g ∫ Σ_i 1_{i,g,l₀}(z,ϵ) π_{l'|l₀} dF_ε   is the gender-
#    pooled probability that a parent born in (l₀,z) works/raises a child in l'.
#    The stationary joint mass solves the eigen-problem
#        Φ_l(z') = Σ_{l₀} Σ_z π_z(z'|z) π̄_{l₀→l}(z) Φ_{l₀}(z),
#    i.e. the dominant left eigenvector of the (2Nz)×(2Nz) row-stochastic kernel
#        P[(l₀,z),(l,z')] = π̄_{l₀→l}(z) · π_z(z'|z).
#    Mass is conserved; we rescale Φ so Σ_l Σ_z Φ_l = M/2.
#
#    Returns Φ as an (Nz × 2) array [z, location] and the migration tensor
#    πbar[l₀, z, l'].
# -----------------------------------------------------------------------------
function stationary_phi(Pi, V, gp::GEParams, grids)
    (; Nz, Nϵ, ϵw, Πz) = grids

    # πbar[l0, z, l'] = gender- and shock-averaged probability that a parent born in
    # (l0, z) ends up raising a child in l'. Integrate the location probabilities Pi
    # over gender (the 0.5) and the two ability shocks (the ϵw quadrature weights),
    # using each state's optimally chosen occupation i★.
    πbar = zeros(2, Nz, 2)
    for l0 in 1:2, zi in 1:Nz, g in 1:2, j in 1:Nϵ, k in 1:Nϵ
        i★, ei = choose_occupation(V, g, l0, zi, j, k)
        w = 0.5 * ϵw[j] * ϵw[k]
        for lp in 1:2
            πbar[l0, zi, lp] += w * Pi[i★, g, l0, zi, ei, lp]
        end
    end

    # Assemble the (2Nz)×(2Nz) row-stochastic kernel on the joint state (location, z):
    #   P[(l0,z) -> (l,z')] = πbar[l0,z,l] · Πz[z,z']  (migrate, then ability evolves).
    # idx flattens (location, z-index) into a single row/column index.
    n   = 2 * Nz
    P   = zeros(n, n)
    idx(l, zi) = (l - 1) * Nz + zi
    for l0 in 1:2, zi in 1:Nz, l in 1:2, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]
    end
    # Stationary distribution = dominant left eigenvector of P, rescaled so the total
    # mass is Mtot/2 (one young agent per parent). Reshaped to Φ[z, location].
    φ  = stationary(P)
    φ .*= gp.Mtot / 2
    Φ  = reshape(φ, Nz, 2)
    return Φ, πbar
end

# -----------------------------------------------------------------------------
# 9. Aggregates implied by the policies and the endogenous distribution Φ.
# -----------------------------------------------------------------------------
function aggregates(Φ, Pi, H, V, gp::GEParams, grids)
    p = gp.p
    (; Nz, Nϵ, ϵw) = grids
    expo  = gp.β / p.σ                                # exponent in the H̃_T aggregator
    HT    = zeros(2)                                  # teacher-HC aggregator, per location
    Inc   = zeros(2)                                  # taxable labor income, per location
    Wbill = zeros(2)                                  # teacher wage bill, per location

    # Integrate policies against the ENDOGENOUS distribution Φ (not the ergodic z law).
    # Outer loops sum over birth state (l0, z); inner loops over gender and shocks.
    for l0 in 1:2, zi in 1:Nz
        wbase = 0.5 * Φ[zi, l0]                       # mass at (l0,z), split by gender
        wbase == 0.0 && continue
        for g in 1:2, j in 1:Nϵ, k in 1:Nϵ
            i★, ei = choose_occupation(V, g, l0, zi, j, k)
            wjk = wbase * ϵw[j] * ϵw[k]               # full quadrature weight for this draw
            h   = H[i★, g, l0, zi, ei]
            for lp in 1:2                             # spread mass over WORK locations via Pi
                pij = Pi[i★, g, l0, zi, ei, lp]
                pij == 0.0 && continue
                wage = i★ == p.T ? p.κ[lp] * h^p.γ : p.A[i★] * h
                Inc[lp] += wjk * pij * (1 - p.τω[i★, g]) * wage
                if i★ == p.T                          # only teachers feed H̃_T and the wage bill
                    Wbill[lp] += wjk * pij * wage
                    HT[lp]    += wjk * pij * h^expo
                end
            end
        end
    end

    # Headcount M_l = 2·(mass in l); tax rate t_l balances the local budget
    # (teacher wage bill / taxable income), clamped to a sane [0, 0.9] range.
    M = [2 * sum(@view Φ[:, l]) for l in 1:2]
    t = [Inc[lp] > 0 ? clamp(Wbill[lp] / Inc[lp], 0.0, 0.9) : 0.0 for lp in 1:2]
    return HT, M, t
end

# -----------------------------------------------------------------------------
# 10. Household block: warm-startable solver used inside the GE loop.
# -----------------------------------------------------------------------------
function recompute_H!(H, S, E, p::Params, grids)
    (; z, ϵgrid, Nz, Nϵ, Q) = grids
    for i in 1:2, g in 1:2, l in 1:2, zi in 1:Nz, ei in 1:Nϵ
        H[i, g, l, zi, ei] = Q[l] * (z[zi] * ϵgrid[ei])^p.α *
                             S[i, g, l, zi, ei]^p.φ * E[i, g, l, zi, ei]^p.η
    end
    return H
end

function solve_household!(S, E, H, V, p::Params, grids; tol = 1e-5, maxit = 200)
    recompute_H!(H, S, E, p, grids)
    for _ in 1:maxit
        Dc, Hc = child_objects(E, H, V, p, grids)
        E0, V0 = copy(E), copy(V)
        update_policy!(S, E, H, V, Dc, Hc, p, grids)
        err = max(maximum(abs, E .- E0), maximum(abs, V .- V0))
        err < tol && break
    end
    return S, E, H, V
end

# -----------------------------------------------------------------------------
# 11. Outer general-equilibrium fixed point over (H̃_T, M, t).
# -----------------------------------------------------------------------------
function solve_ge(gp::GEParams = GEParams();
                  Nz = 5, Nϵ = 7, damping = 0.5, tol = 1e-4, maxit = 100,
                  hh_tol = 1e-6, hh_maxit = 200, verbose = true)
    p = gp.p
    prev_verbose = INFEASIBLE_VERBOSE[]
    INFEASIBLE_VERBOSE[] = false

    # Initial guess for the aggregate vector the outer loop solves for.
    HT = [1.0, 1.0]                                  # teacher-HC aggregator H̃_T
    M  = [gp.Mtot / 2, gp.Mtot / 2]                  # location headcounts (start even)
    t  = [0.10, 0.10]                                # local tax rates
    # Policy/value arrays, shared (and warm-started) across GE iterations.
    dims = (2, 2, 2, Nz, Nϵ)                         # [occ, gender, loc, z, ϵ]
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
            solve_household!(S, E, H, V, p, grids; tol = hh_tol, maxit = hh_maxit)
            Pi          = location_probs(S, E, H, V, p, grids)
            Φ, πbar     = stationary_phi(Pi, V, gp, grids)
            HTn, Mn, tn = aggregates(Φ, Pi, H, V, gp, grids)

            # Convergence measured on the relative move of (H̃_T, M) and absolute move of t.
            err = max(maximum(abs, (HTn .- HT) ./ HT),
                      maximum(abs, (Mn  .- M ) ./ M ),
                      maximum(abs,  tn  .- t))
            # Dampened update toward the new aggregates (stabilises the fixed point).
            @. HT = (1 - damping) * HT + damping * HTn
            @. M  = (1 - damping) * M  + damping * Mn
            @. t  = (1 - damping) * t  + damping * tn

            verbose && @printf("GE %3d  err=%.3e  H̃_T=(%.4f,%.4f)  M=(%.4f,%.4f)  t=(%.4f,%.4f)\n",
                               it, err, HT[1], HT[2], M[1], M[2], t[1], t[2])
            if err < tol
                verbose && println("GE aggregates converged.")
                break
            end
        end

        grids   = build_grids(p; H̃T = HT, M = M, t = t, Nz, Nϵ)
        solve_household!(S, E, H, V, p, grids; tol = hh_tol, maxit = hh_maxit)
        Pi      = location_probs(S, E, H, V, p, grids)
        Φ, πbar = stationary_phi(Pi, V, gp, grids)
        return (; S, E, H, V, Pi, Φ, πbar, HT, M, t, grids, gp)
    finally
        INFEASIBLE_VERBOSE[] = prev_verbose
    end
end

# -----------------------------------------------------------------------------
# 12. Diagnostics
# -----------------------------------------------------------------------------

"Gender-pooled teaching share among young born in l, weighted by the *endogenous* G_l(z)."
function teaching_share_endo(V, l, Φ, grids)
    (; Nz, Nϵ, ϵw) = grids
    mass = sum(@view Φ[:, l])
    s = 0.0
    for zi in 1:Nz
        gz = Φ[zi, l] / mass
        for g in 1:2, j in 1:Nϵ, k in 1:Nϵ
            i★, _ = choose_occupation(V, g, l, zi, j, k)
            i★ == 1 && (s += 0.5 * gz * ϵw[j] * ϵw[k])
        end
    end
    return s
end

function report_ge(sol)
    (; V, Φ, πbar, HT, M, t, grids) = sol
    (; z, Nz, Q, Πz) = grids

    println("\n========== General-equilibrium stationary solution ==========")
    @printf("  H̃_T    = (%.4f, %.4f)\n", HT[1], HT[2])
    @printf("  M      = (%.4f, %.4f)      [total %.4f]\n", M[1], M[2], sum(M))
    @printf("  t      = (%.4f, %.4f)\n", t[1], t[2])
    @printf("  Q      = (%.4f, %.4f)\n", Q[1], Q[2])

    erg = stationary(Πz)
    println("\n  Endogenous ability distribution  G_l(z)   (vs ergodic G*):")
    @printf("    %10s %10s %10s %10s\n", "z", "G_1(z)", "G_2(z)", "G*(z)")
    for zi in 1:Nz
        g1 = Φ[zi, 1] / sum(@view Φ[:, 1])
        g2 = Φ[zi, 2] / sum(@view Φ[:, 2])
        @printf("    %10.4f %10.4f %10.4f %10.4f\n", z[zi], g1, g2, erg[zi])
    end
    mean_z(w) = sum(z[zi] * w[zi] for zi in 1:Nz) / sum(w)
    @printf("    mean z:  loc 1 = %.4f   loc 2 = %.4f   ergodic = %.4f\n",
            mean_z(@view Φ[:, 1]), mean_z(@view Φ[:, 2]), mean_z(erg))

    println("\n  Migration matrix  Π̄[born→work]  (rows sum to 1):")
    for l0 in 1:2
        mass = sum(@view Φ[:, l0])
        row = [sum(Φ[zi, l0] * πbar[l0, zi, l] for zi in 1:Nz) / mass for l in 1:2]
        @printf("    born %d :  ->1 %.4f   ->2 %.4f\n", l0, row[1], row[2])
    end

    println("\n  Teaching share among young born in l (endogenous G_l):")
    for l in 1:2
        @printf("    born %d :  %.4f\n", l, teaching_share_endo(V, l, Φ, grids))
    end
    println("=============================================================\n")
end

# -----------------------------------------------------------------------------
# 13. Run
# -----------------------------------------------------------------------------
function main()
    gp = GEParams()
    println("Solving spatial model in general equilibrium ...")
    sol = solve_ge(gp; Nz = 3, Nϵ = 15, damping = 0.5, tol = 1e-5)
    report_ge(sol)
    return sol
end


# @time main()