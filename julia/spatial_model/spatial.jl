# =============================================================================
# Spatial occupational-choice model with teacher spillovers + altruism.
# PARTIAL-EQUILIBRIUM household problem:  L = 2 locations, I = 2 occupations.
#
# We take the aggregates as given:
#       H̃_T[l]  teacher human-capital stock,
#       M[l]     population,
#       t[l]     local labor-income tax
#
# The household problem is itself a fixed point across generations: a parent's
# old-age consumption includes the share ε of her child's education outlay, and
# in a stationary equilibrium the child's policy equals the parent's policy.
# We therefore iterate the policy functions (s, e, occupation, location) until
# they converge.
#
# Symbol note:  ε  (\varepsilon) = parental cost share  (scalar, p.ε)
#               ϵ  (\epsilon)    = idiosyncratic ability shock (grid ϵgrid).
# =============================================================================

using LinearAlgebra
using Printf
using Optim

# -----------------------------------------------------------------------------
# 1. Parameters.  Occupation 1 = Teacher (T); occupation 2 = non-teaching.
#    Gender index: 1 = m, 2 = f.  Distortion matrices are [occupation, gender].
# -----------------------------------------------------------------------------
Base.@kwdef struct Params
    T::Int = 1                                   # index of the teaching occupation
    A::Vector{Float64} = [NaN, 1.0]              # non-teaching productivity (A[T] unused)
    # Reduced-form schooling tech: h = Q_l (zϵ)^α s^φ e^η ,  Q_l = (2 H̃_T/M)^σ
    α::Float64 = 0.30                            # ability elasticity
    φ::Float64 = 0.40                            # time-investment elasticity
    η::Float64 = 0.20                            # goods-investment elasticity
    σ::Float64 = 0.10                            # teacher-spillover curvature (in Q)
    # Teaching wage:  ω_T,l'(h) = κ_l' h^γ ;  non-teaching wage = A_i h
    γ::Float64 = 0.80
    κ::Vector{Float64} = [0.90, 1.10]
    # Preferences
    μ::Float64 = 1.00                            # weight on log consumption
    λ::Float64 = 0.50                            # altruism strength, f(h') = log h'
    ε::Float64 = 0.50                            # parental share of child's education
    # Locations: amenity, Gumbel scale, utility moving-cost matrix τ[l,l']
    B::Vector{Float64}      = [0.0, 0.10]
    σν::Float64             = 0.50
    τmove::Matrix{Float64}  = [0.0 0.30; 0.30 0.0]
    # Gender-/occupation-specific distortions
    τω::Matrix{Float64} = [0.0 0.00; 0.0 0.15]   # labor-income wedge  (1-τω)
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
    for n in 3:N                                 # build the (n×n) matrix recursively
        Πp = Π
        Π  = zeros(n, n)
        @views Π[1:n-1, 1:n-1] .+= p     .* Πp
        @views Π[1:n-1, 2:n  ] .+= (1-p) .* Πp
        @views Π[2:n  , 1:n-1] .+= (1-p) .* Πp
        @views Π[2:n  , 2:n  ] .+= p     .* Πp
        @views Π[2:n-1, :]     ./= 2             # halve interior rows (avoid double count)
    end
    σz   = σ / sqrt(1 - ρ^2)                     # unconditional sd of log z
    ψ    = sqrt(N - 1) * σz                      # grid half-width
    logz = range(-ψ, ψ; length = N)
    return exp.(logz), Π                         # z in levels, transition matrix
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
function stationary(Π)
    π = fill(1 / size(Π, 1), size(Π, 1))
    for _ in 1:2000
        π = vec(π' * Π)
    end
    return π ./ sum(π)
end

"Bundle grids + given aggregates (Q, t) into a NamedTuple."
function build_grids(p::Params; H̃T, M, t, Nz = 3, Nϵ = 5)
    z, Πz     = rouwenhorst(Nz, p.ρz, p.σξ)
    ϵgrid, ϵw = lognormal_grid(Nϵ, p.σϵ)
    Q = (2 .* H̃T ./ M) .^ p.σ                    # teacher-quality shifter per location
    return (; z, Πz, ϵgrid, ϵw, Nz, Nϵ, Q, t)
end

# -----------------------------------------------------------------------------
# 3. Stage-2 value of each work location, and the Gumbel log-sum
# -----------------------------------------------------------------------------

# Child-side objects entering a parent's consumption: for a child born in `l`
# with gender g, persistent state z'(zi), and shock pair (ϵ1=j, ϵ2=k), the child
# optimally chooses an occupation; we store its education outlay D=(1+τe)e' and
# its human capital h'. Built from the current policy guess (E, H, W).
function child_objects(E, H, W, p::Params, grids)
    (; Nz, Nϵ) = grids
    Dc = zeros(2, 2, Nz, Nϵ, Nϵ)                 # [birth-loc, gender, z', ϵ1, ϵ2]
    Hc = zeros(2, 2, Nz, Nϵ, Nϵ)
    for l in 1:2, g in 1:2, zi in 1:Nz, j in 1:Nϵ, k in 1:Nϵ
        # Compare teaching (own shock j) vs non-teaching (own shock k).
        i′, ei = W[1, g, l, zi, j] ≥ W[2, g, l, zi, k] ? (1, j) : (2, k)
        Dc[l, g, zi, j, k] = (1 + p.τe[i′, g]) * E[i′, g, l, zi, ei]
        Hc[l, g, zi, j, k] = H[i′, g, l, zi, ei]
    end
    return Dc, Hc
end

# Expected value V_l' of working in each location, given the parent's own
# (occupation i, gender g, birth loc l, ability zi, shock ei) and choice (s,e).
# Returns the length-2 vector V and the parent's human capital h.
function location_values(i, g, l, zi, ei, s, e, Dc, Hc, p::Params, grids)
    (; z, ϵgrid, ϵw, Πz, Nz, Nϵ, Q, t) = grids
    h = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
    V = fill(-Inf, 2)
    for lp in 1:2
        wage = i == p.T ? p.κ[lp] * h^p.γ : p.A[i] * h
        # Part of consumption independent of the child's realized state:
        Y = (1 - t[lp]) * (1 - p.τω[i, g]) * wage - (1 - p.ε) * (1 + p.τe[i, g]) * e
        EV = 0.0; feasible = true
        for gp in 1:2, zpi in 1:Nz                # child gender + persistent state
            pz = Πz[zi, zpi]
            pz == 0.0 && continue
            for j in 1:Nϵ, k in 1:Nϵ              # child's fresh idiosyncratic vector
                C = Y - p.ε * Dc[lp, gp, zpi, j, k]
                if C ≤ 0
                    feasible = false; break       # infeasible (s,e): bail out
                end
                EV += pz * ϵw[j] * ϵw[k] * 0.5 *
                      (p.μ * log(C) + p.λ * log(Hc[lp, gp, zpi, j, k]))
            end
            feasible || break
        end
        V[lp] = feasible ? EV : -Inf
    end
    return V, h
end

# Gumbel log-sum over locations (folds in amenity B and moving cost τ) and the
# associated choice probabilities π_{l'|l}.
function logsum_probs(V, l, p::Params)
    x   = ((V[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in 1:2)
    x   = collect(x)
    m   = maximum(x)
    den = sum(exp.(x .- m))
    return p.σν * (m + log(den)), exp.(x .- m) ./ den
end

# -----------------------------------------------------------------------------
# 4. Inner problem: maximize  log(1-s) + V̄(h(s,e))  over (s,e) for one state.
#    Unconstrained via transforms  s = logistic(x1) ∈ (0,1),  e = exp(x2) > 0.
# -----------------------------------------------------------------------------
function neg_objective(x, i, g, l, zi, ei, Dc, Hc, p::Params, grids)
    s = 1 / (1 + exp(-x[1]))
    e = exp(x[2])
    V, _ = location_values(i, g, l, zi, ei, s, e, Dc, Hc, p, grids)
    all(==(-Inf), V) && return 1e10              # penalize fully infeasible choices
    V̄, _ = logsum_probs(V, l, p)
    return -(log(1 - s) + V̄)                     # Optim minimizes
end

# One sweep: re-optimize (s,e) at every state given the current child objects.
# Warm-starts from the incumbent policy for speed.
function update_policy!(S, E, H, W, Dc, Hc, p::Params, grids)
    (; z, ϵgrid, Nz, Nϵ, Q) = grids
    for i in 1:2, g in 1:2, l in 1:2, zi in 1:Nz, ei in 1:Nϵ
        s0, e0 = S[i, g, l, zi, ei], E[i, g, l, zi, ei]
        x0  = [log(s0 / (1 - s0)), log(e0)]
        res = optimize(x -> neg_objective(x, i, g, l, zi, ei, Dc, Hc, p, grids),
                       x0, NelderMead())
        x = Optim.minimizer(res)
        s, e = 1 / (1 + exp(-x[1])), exp(x[2])
        S[i, g, l, zi, ei] = s
        E[i, g, l, zi, ei] = e
        H[i, g, l, zi, ei] = Q[l] * (z[zi] * ϵgrid[ei])^p.α * s^p.φ * e^p.η
        W[i, g, l, zi, ei] = -Optim.minimum(res)
    end
end

# -----------------------------------------------------------------------------
# 5. Outer intergenerational fixed point of the policy functions
# -----------------------------------------------------------------------------
function solve_household(p::Params, grids; tol = 1e-5, maxit = 300)
    (; z, ϵgrid, Nz, Nϵ, Q) = grids
    dims = (2, 2, 2, Nz, Nϵ)                      # [occ, gender, birth-loc, z, ϵ]
    S = fill(0.40, dims...)                       # time investment   s ∈ (0,1)
    E = fill(0.10, dims...)                       # goods investment  e > 0
    H = zeros(dims...)                            # human capital     (derived)
    W = zeros(dims...)                            # value of occupation i at state
    for I in CartesianIndices(H)                  # initialize H consistent with (S,E)
        i, g, l, zi, ei = Tuple(I)
        H[I] = Q[l] * (z[zi] * ϵgrid[ei])^p.α * S[I]^p.φ * E[I]^p.η
    end

    for it in 1:maxit
        Dc, Hc = child_objects(E, H, W, p, grids) # child policy enters parent's C
        E0, W0 = copy(E), copy(W)
        update_policy!(S, E, H, W, Dc, Hc, p, grids)
        err = max(maximum(abs, E .- E0), maximum(abs, W .- W0))
        @printf("  iter %3d   ‖ΔE,ΔW‖∞ = %.3e\n", it, err)
        if err < tol
            println("  household policy converged.")
            break
        end
    end
    return S, E, H, W
end

# -----------------------------------------------------------------------------
# 6. Diagnostics
# -----------------------------------------------------------------------------

"Population-weighted teaching share for (gender g, birth loc l) using the ergodic z."
function teaching_share(W, g, l, grids, πz)
    (; Nz, Nϵ, ϵw) = grids
    s = 0.0
    for zi in 1:Nz, j in 1:Nϵ, k in 1:Nϵ
        W[1, g, l, zi, j] ≥ W[2, g, l, zi, k] && (s += πz[zi] * ϵw[j] * ϵw[k])
    end
    return s
end

function report(S, E, H, W, p::Params, grids)
    (; Nz, Nϵ, Q, t, Πz) = grids
    πz   = stationary(Πz)
    zmed, ϵmed = (Nz + 1) ÷ 2, (Nϵ + 1) ÷ 2
    gnames, lnames = ("m", "f"), (1, 2)

    println("\nTeacher-quality shifter  Q = ", round.(Q; digits = 4),
            "   |  tax t = ", t)

    println("\nTeaching share (ergodic-z weighted):")
    for l in 1:2, g in 1:2
        @printf("   born in %d, %s : %.3f\n", l, gnames[g], teaching_share(W, g, l, grids, πz))
    end

    println("\nNon-teacher (g=m) at median (z,ϵ): policy + location choice probs")
    for l in 1:2
        # rebuild π_{l'|l} from the converged policy for a representative worker
        s, e = S[2, 1, l, zmed, ϵmed], E[2, 1, l, zmed, ϵmed]
        V, h = location_values(2, 1, l, zmed, ϵmed, s, e, child_objects(E, H, W, p, grids)..., p, grids)
        _, π = logsum_probs(V, l, p)
        @printf("   born %d :  s=%.3f  e=%.3f  h=%.3f  ->  π(stay,move)=(%.3f, %.3f)\n",
                l, s, e, h, π[l], π[3 - l])
    end
end

# -----------------------------------------------------------------------------
# 7. Run the example
# -----------------------------------------------------------------------------
function main()
    p     = Params()
    grids = build_grids(p; H̃T = [1.0, 1.5], M = [1.0, 1.0], t = [0.10, 0.12],
                        Nz = 3, Nϵ = 9)
    println("Solving partial-equilibrium household fixed point ...")
    S, E, H, W = solve_household(p, grids)
    report(S, E, H, W, p, grids)
end

@time main()