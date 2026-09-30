# =============================================================================
# Spatial occupational-choice model with teacher spillovers and altruism.
# Stationary general equilibrium with L locations and I occupations (T = teaching),
# at zero parental cost share (ε = 0), with continuous idiosyncratic shocks.
#
# Solution method
#   1. Non-teaching collapse.  For i ≠ T the household problem depends on (i, ϵ_i)
#      only through income capacity X_{O,i} = Θ_{i,g} ϵ_i^α, Θ_{i,g} = (1−τω_{i,g})A_i.
#      X_O* = max_{i≠T} X_{O,i} is therefore a sufficient statistic, and since each
#      X_{O,i} is lognormal its CDF is ∏_i F_i(x).  The I-dimensional shock integral
#      reduces to two scalars, (ϵ_T, X_O*).
#   2. Continuous occupation margin.  Agents teach iff W_T(ϵ_T) > W_O(X_O*).  Both
#      values are increasing, so the threshold follows from inverting splines, and
#      every aggregate is a quadrature integral against the atomless shock densities.
#   3. Ability z is Rouwenhorst-discretised; the spatial distribution Φ_l(z) is the
#      stationary distribution of the migration × AR(1) kernel.
#   The outer loop is a damped fixed point in (H̃_T, M, t).
#
# Notation: ε (\varepsilon) is the parental cost share; ϵ (\epsilon) is the
# idiosyncratic ability shock.
# =============================================================================

using LinearAlgebra
using Printf
using Optim
using Dierckx
using QuadGK
using Distributions

const C_FLOOR = 1e-6    # consumption below which smooth_log is linearised
const M_FLOOR = 1e-8    # population floor in Q_l (guards 0/0 if a location empties)
const TAX_CAP = 0.9     # local tax rates are clamped to [0, TAX_CAP]

"log(C), continued linearly below δ so the optimiser always sees a finite gradient."
@inline smooth_log(C, δ = C_FLOOR) = C > δ ? log(C) : log(δ) + (C - δ) / δ

is_log_kernel(ψ) = isapprox(ψ, 1.0; atol = 1e-10)

"""
    warm(h, ψ, href)

Warm-glow kernel f(h) = ((h/href)^{1−ψ} − 1)/(1−ψ), with the log limit at ψ = 1.
ψ < 1 makes f convex in log h, so the value of a better school rises with the
child's ability. `expm1` avoids cancellation as ψ → 1.
"""
@inline warm(h, ψ, href) = is_log_kernel(ψ) ? log(h / href) :
                           expm1((1 - ψ) * log(h / href)) / (1 - ψ)

# -----------------------------------------------------------------------------
# 1. Parameters
#
# `Params()` is the baseline, with both sorting mechanisms on: a goods-denominated
# moving cost (`mcost`; `τmove = 0`) and a linear warm glow (ψ = 0).
# `nesting_params()` switches both off and recovers the pre-sorting model.
#
# Frozen normalisations, measured once at `nesting_params()` (Nz = 5,
# nϵT = nXO = 48) and re-checked by `sorting_mechanism_tests` in the script:
#   HBAR_BENCH  mean child human capital h̄ → warm-glow reference `href`
#   CBAR_BENCH  mean consumption C̄         → `Cbar`
#   MCOST_BENCH goods cost whose utility value at C̄ equals TAUMOVE_BENCH (μ = 1)
# -----------------------------------------------------------------------------
const TAUMOVE_BENCH = 0.20
const HBAR_BENCH    = 0.1685483614
const CBAR_BENCH    = 0.1619724217
const MCOST_BENCH   = (1 - exp(-TAUMOVE_BENCH / 1.0)) * CBAR_BENCH

"L×L matrix with `x` off the diagonal and zeros on it."
offdiag(x, L) = [l == lp ? 0.0 : Float64(x) for l in 1:L, lp in 1:L]

Base.@kwdef struct Params
    # occupations and technology
    T::Int             = 1             # index of teaching within 1:I
    A::Vector{Float64} = [NaN, 1.5]    # productivity of non-teaching occupations (A[T] unused)
    α::Float64 = 0.30                  # ability elasticity of human capital
    φ::Float64 = 0.40                  # time (s) elasticity of human capital
    η::Float64 = 0.20                  # goods (e) elasticity of human capital
    σ::Float64 = 0.25                  # curvature of the teacher-quality index Q
    γ::Float64 = 0.80                  # human-capital elasticity of the teaching wage
    β::Float64 = 0.15                  # H̃_T aggregates teachers' h^{β/σ}
    κ::Vector{Float64} = [0.75, 0.9]   # teaching-wage shifter by location
    # preferences
    μ::Float64    = 1.00               # weight on log consumption
    λ::Float64    = 0.70               # altruism weight on the warm glow f(h')
    ψ::Float64    = 0.0                # warm-glow curvature (ψ = 1: log kernel)
    href::Float64 = HBAR_BENCH         # warm-glow normalisation (fixed, not solved in GE)
    # locations
    B::Vector{Float64}     = [0.0, 0.1]              # amenity by location
    σν::Float64            = 0.20                    # scale of the Gumbel location taste shock
    τmove::Matrix{Float64} = zeros(2, 2)             # utility cost of moving l → l'
    mcost::Matrix{Float64} = offdiag(MCOST_BENCH, 2) # goods cost of moving l → l'
    Cbar::Float64          = CBAR_BENCH              # consumption level `mcost` is denominated at
    # wedges by (occupation, gender)
    τω::Matrix{Float64} = [0.0 0.0; 0.0 0.1]   # labour-income wedge; income keeps (1 − τω)
    τe::Matrix{Float64} = [0.0 0.0; 0.0 0.0]   # education barrier; goods cost is (1 + τe)e
    # ability process and population
    ρz::Float64   = 0.9                # AR(1) persistence of log z
    σξ::Float64   = 0.20               # AR(1) innovation std of log z
    σϵ::Float64   = 0.30               # std of log ϵ (iid across occupations)
    Mtot::Float64 = 2.0                # total population, both genders
end

"""
    Params(I, L; ...)

`I` occupations and `L` locations. Location 1 is the base; locations 2..L share
`B_other`, `κ_other` and the off-diagonal moving costs `τmove_off` (utility) and
`mcost_off` (goods). Non-teaching occupations share `A_other`; `τω_other` is the
wedge on occupation 2 for gender 2. Other fields can be passed as keywords.

`mcost_off` is an absolute goods amount. An economy whose consumption level is far
from `CBAR_BENCH` should re-derive it with `redenominate_move_cost`.
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
    τω = zeros(I, 2)
    I ≥ 2 && (τω[2, 2] = τω_other)
    defaults = (; T, A, τω, τe = zeros(I, 2),
                κ = vcat(κ_base, fill(κ_other, L - 1)),
                B = vcat(B_base, fill(B_other, L - 1)),
                τmove = offdiag(τmove_off, L), mcost = offdiag(mcost_off, L))
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

"""
    nesting_params(; kwargs...)

The pre-sorting model: utility moving cost (`τmove = TAUMOVE_BENCH`, `mcost = 0`)
and log warm glow (`ψ = 1`, `href = 1`). This is the regression benchmark and the
point where `HBAR_BENCH` and `CBAR_BENCH` are measured.
"""
function nesting_params(; kwargs...)
    defaults = (; τmove = offdiag(TAUMOVE_BENCH, 2), mcost = zeros(2, 2),
                ψ = 1.0, href = 1.0, Cbar = 1.0)
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

"Copy of `p` with the named fields replaced."
with_params(p::Params; kwargs...) =
    Params(; merge(NamedTuple(k => getfield(p, k) for k in fieldnames(Params)),
                   NamedTuple(kwargs))...)

"""
    redenominate_move_cost(p, Cbar) -> Params

Convert the utility moving cost `p.τmove` into the goods cost m = (1 − e^{−τmove/μ})·C̄,
the loss whose utility value at consumption C̄ equals τmove. Returns `p` with
`τmove = 0`, `mcost = m` and `Cbar` recorded. `Params().mcost` is this conversion
applied to `nesting_params()` at `CBAR_BENCH`.
"""
function redenominate_move_cost(p::Params, Cbar::Float64)
    L = length(p.B)
    m = [(1 - exp(-p.τmove[l, lp] / p.μ)) * Cbar for l in 1:L, lp in 1:L]
    m[diagind(m)] .= 0.0
    return with_params(p; τmove = zeros(L, L), mcost = m, Cbar)
end

function check_params(p::Params)
    I, L = length(p.A), length(p.B)
    @assert I ≥ 2 && 1 ≤ p.T ≤ I "need a teaching and at least one non-teaching occupation"
    @assert length(p.κ) == L && size(p.τmove) == (L, L) && size(p.mcost) == (L, L) &&
            size(p.τω) == (I, 2) && size(p.τe) == (I, 2) "Params dimensions are inconsistent"
    @assert all(iszero, diag(p.mcost)) "mcost must have a zero diagonal"
    # The X_O* collapse needs a common τe across non-teaching occupations.
    nonteach = [i for i in 1:I if i != p.T]
    @assert all(p.τe[i, g] == p.τe[nonteach[1], g] for i in nonteach, g in 1:2) "τe must be common across non-teaching occupations"
    return nothing
end

# -----------------------------------------------------------------------------
# 2. Discretisation
# -----------------------------------------------------------------------------
"Rouwenhorst discretisation of log z' = ρ log z + σξ, ξ ~ N(0,1). Returns (z, Π)."
function rouwenhorst(N::Int, ρ::Float64, σ::Float64)
    q = (1 + ρ) / 2
    Π = [q 1-q; 1-q q]
    for n in 3:N
        Πp = Π
        Π = zeros(n, n)
        @views Π[1:n-1, 1:n-1] .+= q     .* Πp
        @views Π[1:n-1, 2:n  ] .+= (1-q) .* Πp
        @views Π[2:n  , 1:n-1] .+= (1-q) .* Πp
        @views Π[2:n  , 2:n  ] .+= q     .* Πp
        @views Π[2:n-1, :]     ./= 2
    end
    zmax = sqrt(N - 1) * σ / sqrt(1 - ρ^2)
    return exp.(range(-zmax, zmax; length = N)), Π
end

"Stationary distribution π = πΠ of a row-stochastic matrix, by power iteration."
function stationary(Π; tol = 1e-12, max_iter = 10_000)
    π = fill(1 / size(Π, 1), size(Π, 1))
    for _ in 1:max_iter
        π_new = vec(π' * Π)
        converged = maximum(abs, π_new .- π) < tol
        π = π_new
        converged && break
    end
    return π ./ sum(π)
end

"Force a vector to be strictly increasing (for spline inversion)."
function make_increasing(v::AbstractVector)
    w = collect(float.(v))
    @inbounds for i in 2:length(w)
        w[i] = max(w[i], w[i-1] + 1e-12)
    end
    return w
end

"""
    build_grids(p; H̃T, M, t, Nz, nϵT, nXO, q_lo, q_hi)

Grids and shock laws at aggregates (H̃_T, M, t):
- `z`, `Πz`: Rouwenhorst ability nodes and transitions;
- `Q`: teacher-quality index Q_l = (2H̃_T/M)^σ;
- `dT`, `ϵTgrid`: law of ϵ_T ~ LogNormal(−σϵ²/2, σϵ) and a log-spaced grid between its
  `q_lo` and `q_hi` quantiles. Teachers come from the upper tail of ϵ_T; a grid spaced
  in quantiles would leave that whole tail (above the 0.98 quantile at nϵT = 48) in one
  interval;
- `dO[i,g]`: law of X_{O,i} = Θ_{i,g} ϵ_i^α; `XOgrid[:,g]`: log-spaced grid between the
  `q_lo` and `q_hi` quantiles of X_O*;
- `lϵT`, `lXO`: the grids in logs, the coordinate of every spline and integral;
- `logΘ[i,g]`: log Θ_{i,g}, to recover h_O from X_O*.
"""
function build_grids(p::Params; H̃T, M, t, Nz = 5, nϵT = 64, nXO = 64,
                     q_lo = 1e-5, q_hi = 1 - 1e-6)
    I, L = length(p.A), length(p.B)
    z, Πz = rouwenhorst(Nz, p.ρz, p.σξ)
    # H̃_T is an integral of a spline and can undershoot zero when a location has
    # almost no teachers; the floor keeps Q real there.
    Q = (2 .* max.(H̃T, M_FLOOR) ./ max.(M, M_FLOOR)) .^ p.σ

    dT     = LogNormal(-p.σϵ^2 / 2, p.σϵ)
    ϵTgrid = exp.(range(log(quantile(dT, q_lo)), log(quantile(dT, q_hi)), nϵT))

    nonteach = [i for i in 1:I if i != p.T]
    dO   = Matrix{LogNormal{Float64}}(undef, I, 2)
    logΘ = fill(NaN, I, 2)
    for g in 1:2, i in nonteach
        logΘ[i, g] = log((1 - p.τω[i, g]) * p.A[i])
        dO[i, g]   = LogNormal(logΘ[i, g] - p.α * p.σϵ^2 / 2, p.α * p.σϵ)
    end

    XOgrid = Matrix{Float64}(undef, nXO, 2)
    for g in 1:2
        lo = quantile_xo(dO, nonteach, g, q_lo)
        hi = quantile_xo(dO, nonteach, g, q_hi)
        XOgrid[:, g] = exp.(range(log(lo), log(hi), nXO))
    end
    return (; z, Πz, Nz, Q, t, I, L, dT, ϵTgrid, dO, XOgrid, lϵT = log.(ϵTgrid), lXO = log.(XOgrid),
              logΘ, nonteach, nϵT, nXO)
end

"""
    quantile_xo(dO, nonteach, g, q)

q-quantile of X_O* = max_{i≠T} X_{O,i}, by bisection in log x. F* ≤ min_i F_i puts
it above every occupation's own q-quantile, and F* ≥ 1 − Σ_i(1 − F_i) puts it below
their (1 − (1−q)/n)-quantiles. With several occupations the lower tail of X_O* is far
above the lowest occupation's, so bounding the grid by individual quantiles would
leave many nodes where X_O* has no mass. Exact with one non-teaching occupation.
"""
function quantile_xo(dO, nonteach, g, q)
    isone(length(nonteach)) && return quantile(dO[nonteach[1], g], q)
    n = length(nonteach)
    a = log(maximum(quantile(dO[i, g], q) for i in nonteach))
    b = log(maximum(quantile(dO[i, g], 1 - (1 - q) / n) for i in nonteach))
    while b - a > 1e-12 * max(1.0, abs(a))
        c = (a + b) / 2
        prod(cdf(dO[i, g], exp(c)) for i in nonteach) < q ? (a = c) : (b = c)
    end
    return exp((a + b) / 2)
end

"CDF of X_O* = max_{i≠T} X_{O,i}: the product of independent lognormal CDFs."
@inline Fxo(gr, g, x) = prod(cdf(gr.dO[i, g], x) for i in gr.nonteach)

"Density of X_O*: f = F · Σ_i f_i/F_i, in one pass over the occupations."
function fxo(gr, g, x)
    F = 1.0
    ratio_sum = 0.0
    for i in gr.nonteach
        Fi = cdf(gr.dO[i, g], x)
        F *= Fi
        Fi > 0.0 && (ratio_sum += pdf(gr.dO[i, g], x) / Fi)
    end
    return F <= 0.0 ? 0.0 : F * ratio_sum
end

"""
    argmax_mean(f, gr, g, x)

E[f(log Θ_{i*}) | X_O* = x], where i* is the non-teaching occupation attaining the
max. P(i* = i | X_O* = x) ∝ f_i(x)/F_i(x). Exact (degenerate) with one
non-teaching occupation.
"""
function argmax_mean(f, gr, g, x)
    (; nonteach, dO, logΘ) = gr
    i0 = nonteach[1]
    isone(length(nonteach)) && return f(logΘ[i0, g])
    num = den = 0.0
    for i in nonteach
        w = pdf(dO[i, g], x) / max(cdf(dO[i, g], x), 1e-300)
        num += w * f(logΘ[i, g])
        den += w
    end
    return den > 0.0 ? num / den : f(logΘ[i0, g])
end

"E[Θ_{i*}^{−(1−ψ)} | X_O* = x]."
EΘpow(gr, g, x, ψ) = argmax_mean(lθ -> exp(-(1 - ψ) * lθ), gr, g, x)

"E[ϵ_{i*} | X_O* = x] = x^{1/α} E[Θ_{i*}^{−1/α} | x]: the shock in the occupation actually chosen."
Eϵ_nonteach(gr, g, x, α) = x^(1 / α) * argmax_mean(lθ -> exp(-lθ / α), gr, g, x)

# -----------------------------------------------------------------------------
# 3. The household's location problem at one node
#
# In both branches post-tax income in work location l' is
#     y_{l'} = c_{l'} · (s^φ e^η)^{γ_b},
# with γ_b = γ for teaching (wage κ h^γ) and γ_b = 1 otherwise:
#     teaching:      c_{l'} = (1−t_{l'})(1−τω_{T,g}) κ_{l'} [Q_l (z ϵ_T)^α]^γ
#     non-teaching:  c_{l'} = (1−t_{l'}) Q_l z^α X_O*
# Consumption is C_{l'} = y_{l'} − (1+τe)e − m_{l,l'}, and the value of l' is
# W_{l'} = μ log C_{l'} + Λ(z, l'). A Gumbel shock over l' gives logit choice.
# -----------------------------------------------------------------------------
struct Node
    c::Vector{Float64}    # income coefficient c_{l'}
    γb::Float64           # income elasticity of s^φ e^η
    τe::Float64           # education barrier
    m::Vector{Float64}    # goods moving cost m_{l,l'}
    Λ::Vector{Float64}    # altruism term Λ(z, l')
    l::Int                # birth location
end

function teach_node(ϵT, zi, l, g, Λ, p::Params, gr)
    k = (1 - p.τω[p.T, g]) * (gr.Q[l] * (gr.z[zi] * ϵT)^p.α)^p.γ
    c = [(1 - gr.t[lp]) * p.κ[lp] * k for lp in 1:gr.L]
    return Node(c, p.γ, p.τe[p.T, g], p.mcost[l, :], Λ[zi, :], l)
end

function nonteach_node(XO, zi, l, g, Λ, p::Params, gr)
    k = gr.Q[l] * gr.z[zi]^p.α * XO
    c = [(1 - gr.t[lp]) * k for lp in 1:gr.L]
    return Node(c, 1.0, p.τe[gr.nonteach[1], g], p.mcost[l, :], Λ[zi, :], l)
end

consumption(nd::Node, s, e, p::Params) =
    nd.c .* (s^p.φ * e^p.η)^nd.γb .- (1 + nd.τe) * e .- nd.m

location_values(nd::Node, C, p::Params) = p.μ .* smooth_log.(C) .+ nd.Λ

"Teacher human capital h_T = Q_l (z ϵ_T)^α s^φ e^η."
teacher_h(ϵT, zi, l, s, e, p::Params, gr) = gr.Q[l] * (gr.z[zi] * ϵT)^p.α * s^p.φ * e^p.η

"Non-teaching income capacity (1−τω)A_{i*} h_O = Q_l z^α X_O* s^φ e^η."
nonteach_cap(XO, zi, l, s, e, p::Params, gr) = gr.Q[l] * gr.z[zi]^p.α * XO * s^p.φ * e^p.η

"Location values W and teacher human capital at a teaching node."
function W_teach(ϵT, zi, l, g, e, s, Λ, p::Params, gr)
    nd = teach_node(ϵT, zi, l, g, Λ, p, gr)
    return location_values(nd, consumption(nd, s, e, p), p), teacher_h(ϵT, zi, l, s, e, p, gr)
end

"Location values W and income capacity at a non-teaching node."
function W_nonteach(XO, zi, l, g, e, s, Λ, p::Params, gr)
    nd = nonteach_node(XO, zi, l, g, Λ, p, gr)
    return location_values(nd, consumption(nd, s, e, p), p), nonteach_cap(XO, zi, l, s, e, p, gr)
end

"Inclusive value σν·log Σ exp((W + B − τmove)/σν) and logit probabilities over work locations."
function logsum_probs(W, l, p::Params)
    v = [(W[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in eachindex(W)]
    vmax = maximum(v)
    ev = exp.(v .- vmax)
    S = sum(ev)
    return p.σν * (vmax + log(S)), ev ./ S
end

"""
    time_share(γb, p, Ξ = 0)

Optimal time share (Proposition 2′). Combining the s- and e-FOCs leaves one
state variable, Ξ = Σ_{l'} π_{l'} m_{l,l'}/C_{l'}. With b = μφ(1+Ξ),
s = bγ_b / (bγ_b + 1 − γ_b η). Without a goods moving cost Ξ = 0 and s is constant.
"""
@inline function time_share(γb, p::Params, Ξ = 0.0)
    b = p.μ * p.φ * (1 + Ξ) * γb
    return b / (b + 1 - γb * p.η)
end

"Ξ = 0 time share of occupation `i`."
s_of(i::Int, p::Params, Ξ = 0.0) = time_share(i == p.T ? p.γ : 1.0, p, Ξ)

"Robust fallback for the goods margin: bounded Brent search over u = log e."
function solve_e_bracket(nd::Node, s, p::Params; lo, hi)
    V(u) = logsum_probs(location_values(nd, consumption(nd, s, exp(u), p), p), nd.l, p)[1]
    e = exp(Optim.minimizer(optimize(u -> -V(u), lo, hi)))
    C = consumption(nd, s, e, p)
    V̄, π = logsum_probs(location_values(nd, C, p), nd.l, p)
    return e, V̄, C, π
end

"""
    solve_e(nd, s, u0, p)

Maximise the inclusive value V̄ over u = log e at fixed s. Newton's method,
warm-started at `u0`, with analytic derivatives. With a = ηγ_b and mc = (1+τe)e,
    w_{l'} = ∂W_{l'}/∂u = μ(a y_{l'} − mc)/C_{l'}
    ∂V̄/∂u  = E_π[w]
    ∂²V̄/∂u² = Var_π(w)/σν + E_π[μ(a² y − mc)/C − w²/μ].
Falls back to `solve_e_bracket` if the log barrier binds, the objective is not
locally concave, or a step leaves [lo, hi]. Returns (e, V̄, C, π).
"""
function solve_e(nd::Node, s, u0, p::Params;
                 lo = log(1e-8), hi = log(1e4), maxit = 30, gtol = 1e-10)
    a  = p.η * nd.γb
    y1 = nd.c .* s^(p.φ * nd.γb)        # income at e = 1
    u  = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    for _ in 1:maxit
        e  = exp(u)
        y  = y1 .* e^a
        mc = (1 + nd.τe) * e
        C  = y .- mc .- nd.m
        V̄, π = logsum_probs(location_values(nd, C, p), nd.l, p)
        any(≤(C_FLOOR), C) && break
        w = @. p.μ * (a * y - mc) / C
        g = dot(π, w)
        abs(g) < gtol && return e, V̄, C, π
        H = (dot(π, w .^ 2) - g^2) / p.σν + dot(π, @. p.μ * (a^2 * y - mc) / C - w^2 / p.μ)
        H ≥ -1e-14 && break
        u -= g / H
        (isfinite(u) && lo ≤ u ≤ hi) || break
    end
    return solve_e_bracket(nd, s, p; lo, hi)
end

"""
    solve_node(nd, p, u0, s0; maxit = 12, tol = 1e-10)

Joint (e, s) solution at one node: alternate the Newton e-solve with the explicit
s(Ξ) update. The map contracts by about 0.03 per pass and exits after one pass when
there is no goods moving cost. Returns (e, V̄, s, residual, passes).
"""
function solve_node(nd::Node, p::Params, u0, s0; maxit = 12, tol = 1e-10)
    s = (isfinite(s0) && 0.0 < s0 < 1.0) ? s0 : time_share(nd.γb, p)
    e = V̄ = NaN
    resid = Inf
    passes = 0
    while passes < maxit && resid ≥ tol
        passes += 1
        e, V̄, C, π = solve_e(nd, s, u0, p)
        Ξ = sum(π[lp] * nd.m[lp] / max(C[lp], 1e-8) for lp in eachindex(π))
        s_new = time_share(nd.γb, p, Ξ)
        resid = abs(s_new - s)
        s = s_new
    end
    return e, V̄, s, resid, passes
end

# -----------------------------------------------------------------------------
# 4. Occupation margin: teach iff W_T(ϵ_T) > W_O(X_O*)
#
#   teach_wt(ϵ_T)    = P(teach | ϵ_T)  = F_{X_O*}( W_O^{-1}(W_T(ϵ_T)) )
#   nonteach_wt(X_O*) = P(don't | X_O*) = F_{ϵ_T}( W_T^{-1}(W_O(X_O*)) )
#
# Values are splined in the log shocks, where they are close to linear; in levels a
# cubic spline between widely spaced tail nodes misses the values by several units.
# -----------------------------------------------------------------------------
struct ChoiceMaps{G}
    splWT::Spline1D;    splWO::Spline1D      # W_T(log ϵ_T), W_O(log X_O*)
    splWTinv::Spline1D; splWOinv::Spline1D   # their inverses, returning logs
    WTv::Vector{Float64}; WOv::Vector{Float64}
    g::Int
    gr::G
end

# The values should be strictly increasing; `make_increasing` only patches
# round-off. Warn (at most 5 times) if a value vector is materially non-monotone.
const NONMONO_WARNINGS = Ref(0)
function warn_nonmonotone(v::AbstractVector, name::String)
    drop  = maximum(v[i-1] - v[i] for i in 2:length(v); init = 0.0)
    scale = max(abs(v[end] - v[1]), 1e-12)
    if drop > 1e-8 * scale && NONMONO_WARNINGS[] < 5
        NONMONO_WARNINGS[] += 1
        suffix = NONMONO_WARNINGS[] == 5 ? " (suppressing further warnings)" : ""
        @printf("  ⚠ %s non-monotone (max drop %.3e vs range %.3e): occupation threshold may be off%s\n",
                name, drop, scale, suffix)
        flush(stdout)
    end
end

function choice_maps(WTv, WOv, g, gr)
    warn_nonmonotone(WTv, "W_T(ϵ_T)")
    warn_nonmonotone(WOv, "W_O(X_O*)")
    return ChoiceMaps(Spline1D(gr.lϵT, WTv), Spline1D(gr.lXO[:, g], WOv),
                      Spline1D(make_increasing(WTv), gr.lϵT; k = 1),
                      Spline1D(make_increasing(WOv), gr.lXO[:, g]; k = 1),
                      WTv, WOv, g, gr)
end

choice_maps(hh, g, l, zi, gr) = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)

"P(teach | log ϵ_T = u)."
@inline function teach_wt_log(cm::ChoiceMaps, u)
    wT = cm.splWT(u)
    wT <= cm.WOv[1]   && return 0.0
    wT >= cm.WOv[end] && return 1.0
    return Fxo(cm.gr, cm.g, exp(cm.splWOinv(wT)))
end

"P(don't teach | log X_O* = v)."
@inline function nonteach_wt_log(cm::ChoiceMaps, v)
    wO = cm.splWO(v)
    wO <= cm.WTv[1]   && return 0.0
    wO >= cm.WTv[end] && return 1.0
    return cdf(cm.gr.dT, exp(cm.splWTinv(wO)))
end

@inline teach_wt(cm::ChoiceMaps, ϵT) = teach_wt_log(cm, log(ϵT))
@inline nonteach_wt(cm::ChoiceMaps, XO) = nonteach_wt_log(cm, log(XO))

# Integrals of a grid quantity over the teachers / non-teachers of one cell:
# quadrature in the log shock (density f(eᵘ)eᵘ) over the grid, plus the shock mass
# outside the grid at the endpoint values (constant extrapolation; the residual is
# O(1 − q_hi)).
function integrate_teach(q, cm::ChoiceMaps, gr)
    (; dT, ϵTgrid, lϵT) = gr
    lo, hi = ϵTgrid[1], ϵTgrid[end]
    spl = Spline1D(lϵT, q)
    f(u) = (w = teach_wt_log(cm, u); w <= 0.0 ? 0.0 : (ϵ = exp(u); pdf(dT, ϵ) * ϵ * spl(u) * w))
    return quadgk(f, lϵT[1], lϵT[end])[1] +
           cdf(dT, lo)  * q[1]   * teach_wt(cm, lo) +
           ccdf(dT, hi) * q[end] * teach_wt(cm, hi)
end

function integrate_nonteach(q, cm::ChoiceMaps, gr)
    g = cm.g
    lx = view(gr.lXO, :, g)
    lo, hi = gr.XOgrid[1, g], gr.XOgrid[end, g]
    spl = Spline1D(lx, q)
    f(v) = (w = nonteach_wt_log(cm, v); w <= 0.0 ? 0.0 : (x = exp(v); fxo(gr, g, x) * x * spl(v) * w))
    return quadgk(f, lx[1], lx[end])[1] +
           Fxo(gr, g, lo)       * q[1]   * nonteach_wt(cm, lo) +
           (1 - Fxo(gr, g, hi)) * q[end] * nonteach_wt(cm, hi)
end

"(∫ over teachers of qT(ϵ_T), ∫ over non-teachers of qO(X_O*)) for one cell."
integrate_choice(qT, qO, cm::ChoiceMaps, gr) =
    (integrate_teach(qT, cm, gr), integrate_nonteach(qO, cm, gr))

# -----------------------------------------------------------------------------
# 5. Altruism term
#
#   Λ(z, l') = λ Σ_{z'} Πz(z, z') · ½ Σ_{g'} E[f(h') | l', g', z'],
# where the expectation integrates the child's occupation choice. On the
# non-teaching side h_O = cap/Θ_{i*} and only the law of Θ_{i*} given X_O* is known,
# so E[f(h_O) | X_O*] is taken over Θ_{i*} in the space the kernel needs.
# -----------------------------------------------------------------------------
function Ewarm_nonteach(cap, x, g, p::Params, gr)
    is_log_kernel(p.ψ) && return log(cap / p.href) - argmax_mean(identity, gr, g, x)
    u = (1 - p.ψ) * log(cap / p.href) + log(EΘpow(gr, g, x, p.ψ))
    return expm1(u) / (1 - p.ψ)
end

function child_Ef(l, g, zi, hh, p::Params, gr)
    cm = choice_maps(hh, g, l, zi, gr)
    fT = warm.(hh.hT[g, l, zi, :], p.ψ, p.href)
    fO = [Ewarm_nonteach(hh.capO[g, l, zi, k], gr.XOgrid[k, g], g, p, gr) for k in 1:gr.nXO]
    return sum(integrate_choice(fT, fO, cm, gr))
end

function compute_lambda(hh, p::Params, gr)
    (; Nz, L, Πz) = gr
    iszero(p.λ) && return zeros(Nz, L)
    Ef = [(child_Ef(l, 1, zi, hh, p, gr) + child_Ef(l, 2, zi, hh, p, gr)) / 2
          for zi in 1:Nz, l in 1:L]
    return p.λ .* (Πz * Ef)
end

# -----------------------------------------------------------------------------
# 6. Household block: iterate V → policies → Λ → V to a fixed point.
# -----------------------------------------------------------------------------
"""
    solve_household(p, gr; tol, maxit, hh0, Λ0)

Solve every (gender, birth location, ability, shock) node given Λ, update Λ,
and repeat until the values converge. `hh0` warm-starts the policies (e, s) and,
through `Λ0`, the altruism term. Arrays are indexed [g, l, zi, k].

Returns the policies and values plus `s_resid`/`s_passes`, the largest (e, s)
fixed-point residual and pass count over the nodes in the final sweep.
"""
function solve_household(p::Params, gr; tol = 1e-6, maxit = 500,
                         verbose = false, print_every = 25,
                         hh0 = nothing, Λ0 = hh0 === nothing ? nothing : hh0.Λ)
    (; Nz, L, nϵT, nXO) = gr
    szT, szO = (2, L, Nz, nϵT), (2, L, Nz, nXO)
    eT = zeros(szT); sT = fill(s_of(p.T, p), szT);           WT = zeros(szT); hT   = zeros(szT)
    eO = zeros(szO); sO = fill(s_of(gr.nonteach[1], p), szO); WO = zeros(szO); capO = zeros(szO)
    if hh0 !== nothing && size(hh0.eT) == szT && size(hh0.eO) == szO
        eT .= hh0.eT; sT .= hh0.sT; eO .= hh0.eO; sO .= hh0.sO
    end
    Λ  = (Λ0 !== nothing && size(Λ0) == (Nz, L)) ? copy(Λ0) : zeros(Nz, L)
    hh = (; eT, sT, WT, hT, eO, sO, WO, capO, Λ)

    s_resid, s_passes = 0.0, 0
    t0 = time()
    for it in 1:maxit
        WT_old, WO_old = copy(WT), copy(WO)
        s_resid, s_passes = 0.0, 0
        # A zero e (unset) gives u0 = −Inf, which the e-solver treats as a cold start.
        for g in 1:2, l in 1:L, zi in 1:Nz
            for k in 1:nϵT
                ϵ  = gr.ϵTgrid[k]
                nd = teach_node(ϵ, zi, l, g, Λ, p, gr)
                e, V̄, s, r, n = solve_node(nd, p, log(eT[g, l, zi, k]), sT[g, l, zi, k])
                eT[g, l, zi, k] = e
                sT[g, l, zi, k] = s
                WT[g, l, zi, k] = log(1 - s) + V̄
                hT[g, l, zi, k] = teacher_h(ϵ, zi, l, s, e, p, gr)
                s_resid, s_passes = max(s_resid, r), max(s_passes, n)
            end
            for k in 1:nXO
                x  = gr.XOgrid[k, g]
                nd = nonteach_node(x, zi, l, g, Λ, p, gr)
                e, V̄, s, r, n = solve_node(nd, p, log(eO[g, l, zi, k]), sO[g, l, zi, k])
                eO[g, l, zi, k]   = e
                sO[g, l, zi, k]   = s
                WO[g, l, zi, k]   = log(1 - s) + V̄
                capO[g, l, zi, k] = nonteach_cap(x, zi, l, s, e, p, gr)
                s_resid, s_passes = max(s_resid, r), max(s_passes, n)
            end
        end
        Λ .= compute_lambda(hh, p, gr)
        err = max(maximum(abs, WT .- WT_old), maximum(abs, WO .- WO_old))
        if verbose && (it == 1 || it % print_every == 0 || err < tol)
            @printf("    HH %4d/%d  err=%.3e  elapsed=%.1fs\n", it, maxit, err, time() - t0)
            flush(stdout)
        end
        err < tol && break
    end
    return (; hh..., s_resid, s_passes)
end

# -----------------------------------------------------------------------------
# 7. Location choice, stationary distribution and aggregates
# -----------------------------------------------------------------------------
"Logit work-location probabilities πT[g, l, zi, k, l'] and πO[g, l, zi, k, l'] at the policies."
function location_probs(hh, p::Params, gr)
    (; Nz, L, nϵT, nXO) = gr
    πT = zeros(2, L, Nz, nϵT, L)
    πO = zeros(2, L, Nz, nXO, L)
    for g in 1:2, l in 1:L, zi in 1:Nz
        for k in 1:nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p)
            πT[g, l, zi, k, :] .= logsum_probs(location_values(nd, C, p), l, p)[2]
        end
        for k in 1:nXO
            nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p)
            πO[g, l, zi, k, :] .= logsum_probs(location_values(nd, C, p), l, p)[2]
        end
    end
    return πT, πO
end

"""
    stationary_phi(πT, πO, hh, p, gr) -> (Φ, πbar)

`πbar[l, zi, l']` is the migration probability born-in-l → work-in-l', integrated
over the occupation choice and averaged over gender. `Φ[zi, l]` is the per-gender
stationary mass of the kernel P((l,z) → (l',z')) = πbar[l, z, l'] Πz[z, z'].
"""
function stationary_phi(πT, πO, hh, p::Params, gr)
    (; Nz, L, Πz) = gr
    πbar = zeros(L, Nz, L)
    for l in 1:L, zi in 1:Nz, g in 1:2
        cm = choice_maps(hh, g, l, zi, gr)
        for lp in 1:L
            IT, IO = integrate_choice(πT[g, l, zi, :, lp], πO[g, l, zi, :, lp], cm, gr)
            πbar[l, zi, lp] += (IT + IO) / 2
        end
    end
    idx(l, zi) = (l - 1) * Nz + zi
    P = zeros(L * Nz, L * Nz)
    for l in 1:L, zi in 1:Nz, lp in 1:L, zpi in 1:Nz
        P[idx(l, zi), idx(lp, zpi)] = πbar[l, zi, lp] * Πz[zi, zpi]
    end
    Φ = stationary(P) .* (p.Mtot / 2)
    return reshape(Φ, Nz, L), πbar
end

"""
    aggregates(Φ, πT, πO, hh, p, gr) -> (H̃_T, M, t)

By work location: the teacher aggregator H̃_T = ∫ h_T^{β/σ} over teachers,
population M, and the balanced-budget tax rate t = teacher wage bill / taxable income.
"""
function aggregates(Φ, πT, πO, hh, p::Params, gr)
    (; Nz, L) = gr
    expo = p.β / p.σ
    HT = zeros(L); income = zeros(L); wagebill = zeros(L)
    for l in 1:L, zi in 1:Nz
        mass = Φ[zi, l] / 2                      # per gender
        iszero(mass) && continue
        for g in 1:2
            cm   = choice_maps(hh, g, l, zi, gr)
            hTv  = hh.hT[g, l, zi, :]
            capO = hh.capO[g, l, zi, :]
            for lp in 1:L
                πTv = πT[g, l, zi, :, lp]
                wT  = integrate_teach(πTv .* p.κ[lp] .* hTv .^ p.γ, cm, gr)
                HT[lp]       += mass * integrate_teach(πTv .* hTv .^ expo, cm, gr)
                wagebill[lp] += mass * wT
                income[lp]   += mass * ((1 - p.τω[p.T, g]) * wT +
                                        integrate_nonteach(πO[g, l, zi, :, lp] .* capO, cm, gr))
            end
        end
    end
    M = 2 .* vec(sum(Φ; dims = 1))
    t = [income[lp] > 0 ? clamp(wagebill[lp] / income[lp], 0.0, TAX_CAP) : 0.0 for lp in 1:L]
    return HT, M, t
end

# -----------------------------------------------------------------------------
# 8. General equilibrium: damped fixed point in (H̃_T, M, t)
# -----------------------------------------------------------------------------
fmtvec(v) = "(" * join((@sprintf("%.4f", x) for x in v), ", ") * ")"
reldiff(new, old; floor = 1e-3) = maximum(abs.(new .- old) ./ max.(abs.(old), floor))

"""
    solve_ge(p = Params(); Nz, nϵT, nXO, damping, tol, maxit, hh_tol, hh_maxit,
             verbose, init, q_lo, q_hi)

Each iteration rebuilds the grids at the current (H̃_T, M, t), solves the household
problem (warm-started from the previous iteration), computes location choices, Φ
and the implied aggregates, and takes a damped step. Convergence is on the max
relative change in H̃_T and M and absolute change in t. `init` is a previous
solution to continue from. A final solve at the converged aggregates makes the
returned objects mutually consistent.
"""
function solve_ge(p::Params = Params();
                  Nz = 5, nϵT = 64, nXO = 64, damping = 0.3, tol = 1e-5,
                  maxit = 1000, hh_tol = 1e-6, hh_maxit = 1000, verbose = true,
                  init = nothing, q_lo = 1e-5, q_hi = 1 - 1e-6)
    check_params(p)
    I, L = length(p.A), length(p.B)
    grid_kw = (; Nz, nϵT, nXO, q_lo, q_hi)

    if verbose
        @printf("  [ε=0, continuous]  I=%d, L=%d, Nz=%d  shock grids: nϵT=%d, nXO=%d\n",
                I, L, Nz, nϵT, nXO)
        print("  Ξ=0 reference s: ")
        for i in 1:I
            @printf("%s s=%.4f   ", i == p.T ? "T" : "occ$i", s_of(i, p))
        end
        @printf("| ψ=%.3f href=%.4f  max mcost=%.4f\n", p.ψ, p.href, maximum(p.mcost))
        flush(stdout)
    end

    HT = fill(1.0, L); M = fill(p.Mtot / L, L); t = fill(0.10, L)
    hh0 = nothing
    if init !== nothing
        if length(init.HT) == L
            HT .= init.HT; M .= init.M; t .= init.t
            hh0 = init.hh
        elseif verbose
            @printf("  ⚠ warm start ignored: init has L=%d ≠ %d\n", length(init.HT), L)
        end
    end

    converged = false
    for it in 1:maxit
        gr = build_grids(p; H̃T = HT, M, t, grid_kw...)
        t_hh = @elapsed (hh0 = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit,
                                                verbose = verbose && it == 1, hh0))
        t_pi  = @elapsed ((πT, πO) = location_probs(hh0, p, gr))
        t_phi = @elapsed ((Φ, _) = stationary_phi(πT, πO, hh0, p, gr))
        t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, πT, πO, hh0, p, gr))

        # An empty location has no meaningful budget: hold its tax fixed.
        live = Mn .≥ M_FLOOR
        tn[.!live] .= t[.!live]
        t_err = any(live) ? maximum(abs, tn[live] .- t[live]) : 0.0
        err = max(reldiff(HTn, HT), reldiff(Mn, M), t_err)

        @. HT = (1 - damping) * HT + damping * HTn
        @. M  = (1 - damping) * M  + damping * Mn
        @. t  = (1 - damping) * t  + damping * tn

        if verbose
            @printf("GE %3d  err=%.3e  H̃_T=%s  M=%s  t=%s  [hh=%.1fs pi=%.1fs phi=%.1fs agg=%.1fs]\n",
                    it, err, fmtvec(HT), fmtvec(M), fmtvec(t), t_hh, t_pi, t_phi, t_agg)
            flush(stdout)
        end
        if err < tol
            converged = true
            verbose && println("GE aggregates converged.")
            break
        end
    end
    verbose && !converged && println("  ⚠ GE did not converge in $maxit iterations.")

    gr = build_grids(p; H̃T = HT, M, t, grid_kw...)
    hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit, hh0)
    Pi_T, Pi_O = location_probs(hh, p, gr)
    Φ, πbar = stationary_phi(Pi_T, Pi_O, hh, p, gr)
    return (; hh, Pi_T, Pi_O, Φ, πbar, HT, M, t, gr, p, converged,
              s_resid = hh.s_resid, s_passes = hh.s_passes)
end

# -----------------------------------------------------------------------------
# 9. Diagnostics and summary statistics
# -----------------------------------------------------------------------------
"""
    verify_solution(sol; δ = C_FLOOR, tclamp = TAX_CAP, s_atol = 1e-8)

Check that the numerical safeguards are slack at the solution: (i) consumption
exceeds δ in every (node, work location), so smooth_log equals log; (ii) no tax
rate is on the clamp [0, tclamp]; (iii) the (e, s) fixed point converged.
Returns true if all hold. Condition (i) is what limits how large `mcost` can be:
the goods cost bites first on low-ability movers.
"""
function verify_solution(sol; δ = C_FLOOR, tclamp = TAX_CAP, s_atol = 1e-8)
    (; hh, t, gr, p, s_resid, s_passes) = sol
    minC = Inf
    n_binding = 0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        for k in 1:gr.nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p)
            minC = min(minC, minimum(C)); n_binding += count(≤(δ), C)
        end
        for k in 1:gr.nXO
            nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p)
            minC = min(minC, minimum(C)); n_binding += count(≤(δ), C)
        end
    end
    ok = true
    if s_resid > s_atol
        @printf("  ⚠ s fixed point residual %.2e (> %.1e) after %d passes\n", s_resid, s_atol, s_passes)
        ok = false
    end
    if n_binding > 0
        @printf("  ⚠ smooth_log barrier BINDING at %d (node, l') cells (min C = %.3e ≤ δ = %.1e)\n",
                n_binding, minC, δ)
        ok = false
    end
    for lp in 1:gr.L
        if t[lp] >= tclamp - 1e-9 || t[lp] <= 0.0
            @printf("  ⚠ tax rate t[%d] = %.4f is on the clamp [0, %.1f]: budget not balanced\n",
                    lp, t[lp], tclamp)
            ok = false
        end
    end
    ok && @printf("  ✓ feasibility audit passed: min C = %.4e > δ = %.1e; taxes interior t = %s; s residual %.1e (≤ %d passes)\n",
                  minC, δ, fmtvec(t), s_resid, s_passes)
    return ok
end

"""
    cell_mean(sol, qT, qO)

Population mean of a quantity given on the grids, weighted by Φ and integrated
over the occupation choice. `qT(g, l, zi)` and `qO(g, l, zi)` return the vectors
on the ϵ_T and X_O* grids for that cell.
"""
function cell_mean(sol, qT, qO)
    (; hh, Φ, gr) = sol
    onesT, onesO = ones(gr.nϵT), ones(gr.nXO)
    num = den = 0.0
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        num += mass * sum(integrate_choice(qT(g, l, zi), qO(g, l, zi), cm, gr))
        den += mass * sum(integrate_choice(onesT, onesO, cm, gr))
    end
    return num / den
end

"""
    mean_child_h(sol)

Mean human capital h̄ at `sol`. The warm-glow reference `href` is frozen at this
value (`HBAR_BENCH` at `nesting_params()`) so that ψ changes only the curvature of
the kernel, not its level.
"""
function mean_child_h(sol)
    (; hh, gr) = sol
    hO(g, l, zi) = [hh.capO[g, l, zi, k] * EΘpow(gr, g, gr.XOgrid[k, g], 0.0) for k in 1:gr.nXO]
    return cell_mean(sol, (g, l, zi) -> hh.hT[g, l, zi, :], hO)
end

"""
    mean_consumption(sol)

Mean consumption C̄, averaged over location choices, at `sol` (`CBAR_BENCH` at
`nesting_params()`). Used by `redenominate_move_cost`.
"""
function mean_consumption(sol)
    (; hh, gr, p, Pi_T, Pi_O) = sol
    CT(g, l, zi) = map(1:gr.nϵT) do k
        nd = teach_node(gr.ϵTgrid[k], zi, l, g, hh.Λ, p, gr)
        dot(Pi_T[g, l, zi, k, :], consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p))
    end
    CO(g, l, zi) = map(1:gr.nXO) do k
        nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hh.Λ, p, gr)
        dot(Pi_O[g, l, zi, k, :], consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p))
    end
    return cell_mean(sol, CT, CO)
end

"Mean time share among teachers and non-teachers, and its (min, max) over all nodes."
function s_summary(sol)
    (; hh, Φ, gr) = sol
    numT = numO = denT = denO = 0.0
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        sT, sO = integrate_choice(hh.sT[g, l, zi, :], hh.sO[g, l, zi, :], cm, gr)
        mT, mO = integrate_choice(ones(gr.nϵT), ones(gr.nXO), cm, gr)
        numT += mass * sT; denT += mass * mT
        numO += mass * sO; denO += mass * mO
    end
    return (; meanT = numT / denT, meanO = numO / denO,
              rangeT = extrema(hh.sT), rangeO = extrema(hh.sO))
end

"Teaching share among the young born in `l` (both genders)."
function teaching_share_endo(hh, l, Φ, gr, p::Params)
    mass = sum(@view Φ[:, l])
    return sum(Φ[zi, l] / mass / 2 * integrate_teach(ones(gr.nϵT), choice_maps(hh, g, l, zi, gr), gr)
               for zi in 1:gr.Nz, g in 1:2)
end

"""
    student_ability_moments(sol)

Ability moments of the young by birth location, split by the occupation chosen:
- `z`: persistent ability. Its mean equals the `mean z` row of `report_ge`; the split
  shows which types select into teaching.
- `a = z·ϵ`: ability in the chosen occupation (teachers: ϵ_T; others: ϵ_{i*} via
  `Eϵ_nonteach`). E[a] exceeds E[z] because both branches select on their own ϵ.

Returns per-location vectors (`mean_a`, `mean_a_teach`, `mean_a_non`, `mean_z`,
`mean_z_teach`, `mean_z_non`, `teach_share`, `mass`, `massT`, `massO`), the
branch split of Φ (`ΦT`, `ΦO`), and population-pooled scalars suffixed `_all`.
"""
function student_ability_moments(sol)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, L, nϵT, nXO) = gr
    onesT, onesO = ones(nϵT), ones(nXO)
    ΦT = zeros(Nz, L); ΦO = zeros(Nz, L)
    aT = zeros(L);     aO = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm  = choice_maps(hh, g, l, zi, gr)
        aTv = z[zi] .* gr.ϵTgrid
        aOv = [z[zi] * Eϵ_nonteach(gr, g, gr.XOgrid[k, g], p.α) for k in 1:nXO]
        mT, mO = integrate_choice(onesT, onesO, cm, gr)
        IaT, IaO = integrate_choice(aTv, aOv, cm, gr)
        ΦT[zi, l] += mass * mT; ΦO[zi, l] += mass * mO
        aT[l]     += mass * IaT; aO[l]    += mass * IaO
    end
    ratio(num, den) = [den[l] > 0 ? num[l] / den[l] : NaN for l in 1:L]
    massT, massO = vec(sum(ΦT; dims = 1)), vec(sum(ΦO; dims = 1))
    mass = massT .+ massO
    zT, zO = vec(z' * ΦT), vec(z' * ΦO)
    return (; mass, massT, massO, ΦT, ΦO,
              teach_share  = ratio(massT, mass),
              mean_a       = ratio(aT .+ aO, mass),
              mean_a_teach = ratio(aT, massT),
              mean_a_non   = ratio(aO, massO),
              mean_z       = ratio(zT .+ zO, mass),
              mean_z_teach = ratio(zT, massT),
              mean_z_non   = ratio(zO, massO),
              teach_share_all  = sum(massT) / sum(mass),
              mean_a_all       = sum(aT .+ aO) / sum(mass),
              mean_a_teach_all = sum(aT) / sum(massT),
              mean_a_non_all   = sum(aO) / sum(massO),
              mean_z_all       = sum(zT .+ zO) / sum(mass),
              mean_z_teach_all = sum(zT) / sum(massT),
              mean_z_non_all   = sum(zO) / sum(massO))
end

"""
    student_ability_density(sol, l, agrid) -> (fT, fO)

Density of student ability a = z·ϵ over `agrid` for those born in `l`, split by
branch and normalised so ∫(fT + fO) = 1 and ∫fT is the teaching share. Evaluated
in closed form, which avoids the comb artefacts of binning quadrature nodes:
- teaching: f_ϵ(a/z)/z · P(teach | ϵ_T = a/z);
- non-teaching: Σ_i f_ϵ(ϵ) ∏_{j≠i} F_j(Θ_i ϵ^α) · P(don't teach | X_O* = Θ_i ϵ^α),
  where the product is the probability that occupation i attains the max.
"""
function student_ability_density(sol, l, agrid)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, dT, dO, nonteach, logΘ) = gr
    fT = zeros(length(agrid)); fO = zeros(length(agrid))
    for zi in 1:Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        for (ia, a) in enumerate(agrid)
            ϵ  = a / z[zi]
            fϵ = mass * pdf(dT, ϵ) / z[zi]
            fT[ia] += fϵ * teach_wt(cm, ϵ)
            for i in nonteach
                x = exp(logΘ[i, g]) * ϵ^p.α
                wins = prod((cdf(dO[j, g], x) for j in nonteach if j != i); init = 1.0)
                fO[ia] += fϵ * wins * nonteach_wt(cm, x)
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
    @printf("  s (choice-weighted mean [min, max]):  T=%.4f [%.4f, %.4f]   nonteach=%.4f [%.4f, %.4f]\n",
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
    mean_z(w) = dot(z, w) / sum(w)
    print("    mean z: ")
    for l in 1:L; @printf(" loc %d = %.4f  ", l, mean_z(@view Φ[:, l])); end
    @printf("ergodic = %.4f\n", mean_z(ergodic))

    println("\n  Migration matrix  Π̄[born→work]  (rows sum to 1):")
    for l in 1:L
        row = [dot(Φ[:, l], πbar[l, :, lp]) / sum(@view Φ[:, l]) for lp in 1:L]
        print("    born $l :")
        for lp in 1:L; @printf("  ->%d %.4f", lp, row[lp]); end
        println()
    end

    println("\n  Teaching share among young born in l (endogenous G_l):")
    for l in 1:L
        @printf("    born %d :  %.4f\n", l, teaching_share_endo(hh, l, Φ, gr, p))
    end

    sa = student_ability_moments(sol)
    println("\n  Student ability by birth location  (choice-weighted; a = z·ϵ in the chosen occupation):")
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
# 10. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model (ε = 0, continuous) in general equilibrium ...")
    sol = solve_ge(Params(); Nz = 3, nϵT = 48, nXO = 48, damping = 0.75,
                   tol = 1e-6, hh_tol = 1e-6, maxit = 500, hh_maxit = 1000)
    report_ge(sol)
    verify_solution(sol)
    return sol
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
