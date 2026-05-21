#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "spatial.jl"))

# ── helpers ──────────────────────────────────────────────────────────────────

function make_small_params(; τmove = [0.0 0.2; 0.2 0.0], B = [0.0, 0.0], κ = [1.0, 1.05], t = [0.2, 0.25])
    A_O = [1.1, 0.95, 1.05]
    p = Params(
        2,
        [:f, :m],
        0.6,   # α
        0.8,   # β
        0.3,   # η
        1.0,   # σ
        1.0,   # μ
        0.4,   # ϕ
        0.9,   # γ
        κ,
        A_O,
        Dict(:f => [0.1, 0.05, 0.0], :m => [0.0, 0.0, 0.0]),  # τw_O
        Dict(:f => [0.05, 0.02, 0.0], :m => [0.0, 0.0, 0.0]), # τe_O
        Dict(:f => 0.0, :m => 0.0),  # τw_T
        Dict(:f => 0.0, :m => 0.0),  # τe_T
        0.2,   # ε_parent
        0.1,   # λ
        B,
        τmove,
        0.2,   # σν
        0.6,   # ρz
        0.25,  # σξ
        0.4,   # σeps
        7,     # Neps
        3,     # Nz
        0.2, 0.2, 0.2, 0.3,  # damp_m, damp_H, damp_t, damp_E
        1e-5,
        300,
    )
    xz, Pz = tauchen(p.Nz, p.ρz, p.σξ)
    z = exp.(xz)
    eps_dist = LogNormal(-0.5 * p.σeps^2, p.σeps)
    eps_nodes, eps_w = quantile_histogram(eps_dist, p.Neps)
    grids = Grids(xz, z, Pz, eps_nodes, eps_w)
    m = fill(1.0 / (p.L * p.Nz), p.L, p.Nz)
    eqm = Eqm(
        m,
        [1.0, 1.0],
        [1.0, 1.0],
        t,
        zeros(p.Nz, p.L, 2),
        zeros(p.Nz, p.L, 2),
        zeros(p.Nz, p.L, 2, p.Neps),
        zeros(p.Nz, p.L, 2, p.Neps),
        zeros(p.Nz, p.L, 2, p.Neps),
        zeros(p.Nz, p.L, 2, p.Neps),
    )
    return p, grids, eqm
end

function compute_child_objects(p::Params, grids::Grids, eqm::Eqm)
    child_logh = Matrix{Float64}(undef, p.Nz, p.L)
    child_cost = Matrix{Float64}(undef, p.Nz, p.L)
    for k in 1:p.Nz
        for ldest in 1:p.L
            elog = 0.0
            eco = 0.0
            for kp in 1:p.Nz
                wkp = grids.Pz[k, kp]
                elog += wkp * 0.5 * (eqm.E_logh[kp, ldest, 1] + eqm.E_logh[kp, ldest, 2])
                eco += wkp * 0.5 * (eqm.E_cost[kp, ldest, 1] + eqm.E_cost[kp, ldest, 2])
            end
            child_logh[k, ldest] = elog
            child_cost[k, ldest] = eco
        end
    end
    return child_logh, child_cost
end

# ── test sets ─────────────────────────────────────────────────────────────────

@testset "Core numerics" begin
    @test logsumexp2(1.0, 1.0) ≈ 1.0 + log(2.0) atol = 1e-12
    p1, p2 = softmax2(1.5, 0.5, 0.3)
    @test p1 + p2 ≈ 1.0 atol = 1e-12
    @test p1 > p2

    xz, Pz = tauchen(5, 0.7, 0.2)
    @test length(xz) == 5
    @test size(Pz) == (5, 5)
    @test all(abs.(sum(Pz, dims = 2) .- 1.0) .< 1e-12)
    π = stationary_dist(Pz)
    @test abs(sum(π) - 1.0) < 1e-12
    @test maximum(abs.((π' * Pz)[:] .- π)) < 1e-10
end

@testset "Closed-form s* limits (Appendix A.4)" begin
    # Single location, no altruism, no parental finance, no ability persistence:
    # the FOCs have closed-form solutions independent of state.
    # s_O* = μϕ / (μϕ + 1 - η)   (eq. A.4 non-teaching)
    # s_T* = μϕγ / ((μϕ - η)γ + 1) (eq. A.4 teaching)
    p, grids, eqm = make_small_params(
        τmove = [0.0 0.0; 0.0 0.0],   # no moving cost → symmetric locations
        B = [0.0, 0.0],
        κ = [1.0, 1.0],               # identical locations
        t = [0.2, 0.2],
    )
    # Override to single-location-like conditions: ε_parent=0, λ=0
    p2 = Params(
        p.L, p.genders,
        p.α, p.β, p.η, p.σ, p.μ, p.ϕ, p.γ,
        p.κ, p.A_O,
        p.τw_O, p.τe_O,
        Dict(:f => 0.0, :m => 0.0),  # τw_T
        Dict(:f => 0.0, :m => 0.0),  # τe_T
        0.0,   # ε_parent = 0
        0.0,   # λ = 0
        p.B, p.τmove, p.σν,
        0.0,   # ρz = 0 → iid z
        p.σξ, p.σeps,
        p.Neps, p.Nz,
        p.damp_m, p.damp_H, p.damp_t, p.damp_E,
        p.tol, p.maxit,
    )

    s_O_expected = s_O_const(p2)
    s_T_expected = s_T_const(p2)

    @test 0.0 < s_O_expected < 1.0
    @test 0.0 < s_T_expected < 1.0

    # Verify the closed-form formula values
    μ, ϕ, η, γ = p2.μ, p2.ϕ, p2.η, p2.γ
    @test s_O_expected ≈ μ * ϕ / (μ * ϕ + 1 - η) atol = 1e-14
    @test s_T_expected ≈ μ * ϕ * γ / ((μ * ϕ - η) * γ + 1) atol = 1e-14

    # With symmetric locations, zero moving costs, ε=0, λ=0, the solver
    # should converge to s values near the closed-form constants.
    eqm2 = Eqm(
        eqm.m_young, eqm.M, eqm.Htilde, [0.2, 0.2],
        eqm.E_logh, eqm.E_cost,
        eqm.UT, eqm.UO_best, eqm.costT, eqm.costO_best,
    )
    child_logh = zeros(p2.Nz, p2.L)
    k = ceil(Int, p2.Nz / 2)  # middle z node
    _, _, _, _, _, _, sT_solved = solve_T(p2, grids, eqm2, child_logh, 1, k, :m, grids.eps_nodes[4])
    _, _, _, _, _, _, sO_solved = solve_O(p2, grids, eqm2, child_logh, 1, k, :m, grids.eps_nodes[4], 1)

    @test sT_solved ≈ s_T_expected atol = 0.02
    @test sO_solved ≈ s_O_expected atol = 0.02
end

@testset "Location-choice logic" begin
    p, grids, eqm = make_small_params()
    child_logh = zeros(p.Nz, p.L)
    k = 2

    _, πO, hO, eO, _, _, _ = solve_O(p, grids, eqm, child_logh, 1, k, :f, grids.eps_nodes[3], 1)
    @test πO[1] + πO[2] ≈ 1.0 atol = 1e-10
    @test 0.0 <= πO[1] <= 1.0
    @test hO > 0.0
    @test eO > 0.0

    _, _, hO_low, _, _, _, _ = solve_O(p, grids, eqm, child_logh, 1, k, :f, grids.eps_nodes[3], argmin(p.A_O))
    @test hO > hO_low

    _, πT, hT, eT, _, _, _ = solve_T(p, grids, eqm, child_logh, 1, k, :m, grids.eps_nodes[3])
    @test πT[1] + πT[2] ≈ 1.0 atol = 1e-10
    @test 0.0 <= πT[1] <= 1.0
    @test hT > 0.0
    @test eT > 0.0

    p_sym, grids_sym, eqm_sym = make_small_params(τmove = [0.0 0.0; 0.0 0.0], B = [0.0, 0.0], κ = [1.0, 1.0], t = [0.2, 0.2])
    _, π_sym, _, _, _, _, _ = solve_O(p_sym, grids_sym, eqm_sym, child_logh, 1, k, :m, grids_sym.eps_nodes[3], 1)
    @test π_sym[1] ≈ 0.5 atol = 1e-6

    p_sticky, grids_sticky, eqm_sticky = make_small_params(τmove = [0.0 8.0; 8.0 0.0], B = [0.0, 0.0], κ = [1.0, 1.0], t = [0.2, 0.2])
    _, π_sticky, _, _, _, _, _ = solve_O(p_sticky, grids_sticky, eqm_sticky, child_logh, 1, k, :m, grids_sticky.eps_nodes[3], 1)
    @test π_sticky[1] > 0.999
end

@testset "state_moments migration sums to 1" begin
    p, grids, eqm = make_small_params()
    child_logh = zeros(p.Nz, p.L)
    for l in 1:p.L, k in 1:p.Nz, gsym in p.genders
        _, _, _, mig, _, _, _, _, _, _, _ = state_moments(p, grids, eqm, child_logh, l, k, gsym)
        @test mig[1] + mig[2] ≈ 1.0 atol = 1e-8
        @test 0.0 <= mig[1] <= 1.0
        @test 0.0 <= mig[2] <= 1.0
    end
end

@testset "Htilde aggregation consistency" begin
    # After calling state_moments, accumulating Hcontrib with proper weights
    # should reproduce eqm.Htilde (up to the fact eqm is not at fixed point).
    # We just check the formula is internally consistent: Htilde > 0 and finite.
    p, grids, eqm = make_small_params()
    child_logh = zeros(p.Nz, p.L)
    H = zeros(p.L)
    for l0 in 1:p.L, k in 1:p.Nz, (ig, gsym) in enumerate(p.genders)
        _, _, _, _, Hc, _, _, _, _, _, _ = state_moments(p, grids, eqm, child_logh, l0, k, gsym)
        wgt = 0.5 * eqm.m_young[l0, k]
        for ldest in 1:p.L
            H[ldest] += wgt * Hc[ldest]
        end
    end
    @test all(isfinite.(H))
    @test all(H .>= 0.0)
    @test sum(H) > 0.0
end

@testset "Stationary equilibrium integration (small grid)" begin
    d = load_baseline("benchmark", 1970)
    τw = Array{Float64}(d["τ_w"])
    τe = Array{Float64}(d["τ_e"])
    κ0 = Float64(d["κ"])
    p = Params(
        2,
        [:f, :m],
        Float64(d["α"]),
        Float64(d["β"]),
        Float64(d["η"]),
        Float64(d["σ"]),
        Float64(d["μ"]),
        Float64(d["ϕ"]),
        Float64(d["γ"]),
        [κ0, 1.03 * κ0],
        Vector{Float64}(d["a_by_occ"]),
        Dict(:f => Vector{Float64}(τw[:, 1]), :m => Vector{Float64}(τw[:, 2])),  # τw_O
        Dict(:f => Vector{Float64}(τe[:, 1]), :m => Vector{Float64}(τe[:, 2])),  # τe_O
        Dict(:f => 0.0, :m => 0.0),  # τw_T
        Dict(:f => 0.0, :m => 0.0),  # τe_T
        0.2,   # ε_parent
        0.1,   # λ
        [0.0, 0.0],
        [0.0 0.2; 0.2 0.0],
        0.2,   # σν
        0.7,   # ρz
        0.4,   # σξ
        0.6,   # σeps
        9,     # Neps
        5,     # Nz
        0.2, 0.2, 0.2, 0.3,  # damp_m, damp_H, damp_t, damp_E
        1e-4,
        500,
    )
    @test n_nonteach_occupations(p) == size(τw, 1)
    @test length(p.τw_O[:f]) == n_nonteach_occupations(p)
    @test length(p.τe_O[:m]) == n_nonteach_occupations(p)

    eqm, grids = solve_stationary(p)

    @test size(eqm.m_young) == (p.L, p.Nz)
    @test all(isfinite.(eqm.m_young))
    @test sum(eqm.m_young) ≈ 1.0 atol = 1e-10
    @test all(eqm.m_young .>= 0.0)
    @test eqm.M[1] ≈ 2.0 * sum(eqm.m_young[1, :]) atol = 1e-8
    @test eqm.M[2] ≈ 2.0 * sum(eqm.m_young[2, :]) atol = 1e-8
    @test sum(eqm.M) ≈ 2.0 atol = 1e-8
    @test all(eqm.Htilde .> 0.0)
    @test all((0.0 .<= eqm.t) .& (eqm.t .<= 0.95))

    child_logh, child_cost = compute_child_objects(p, grids, eqm)

    I = zeros(p.L)
    W = zeros(p.L)
    H = zeros(p.L)
    for l0 in 1:p.L, k in 1:p.Nz, (ig, gsym) in enumerate(p.genders)
        _, _, _, mig, Hc, tb, wb, _, _, _, _ = state_moments(p, grids, eqm, child_logh, l0, k, gsym)
        @test mig[1] + mig[2] ≈ 1.0 atol = 1e-8
        @test (0.0 <= mig[1] <= 1.0) && (0.0 <= mig[2] <= 1.0)
        wgt = 0.5 * eqm.m_young[l0, k]
        for ldest in 1:p.L
            H[ldest] += wgt * Hc[ldest]
            I[ldest] += wgt * tb[ldest]
            W[ldest] += wgt * wb[ldest]
        end
    end

    # H̃_T aggregation consistency (§8 aggregator should match stored Htilde)
    @test maximum(abs.(H .- eqm.Htilde)) < 5e-4

    # Local budget balance: t_l * I_l ≈ W_l (§9)
    @test maximum(abs.(eqm.t .* I .- W) ./ max.(abs.(W), 1e-12)) < 5e-3
end

@testset "Symmetric equilibrium" begin
    # With identical primitives in both locations (same κ, same t, zero moving costs,
    # zero amenity gap) the equilibrium must be symmetric: M[1] ≈ M[2], Htilde[1] ≈ Htilde[2].
    d = load_baseline("benchmark", 1970)
    τw = Array{Float64}(d["τ_w"])
    τe = Array{Float64}(d["τ_e"])
    κ0 = Float64(d["κ"])
    p_sym = Params(
        2,
        [:f, :m],
        Float64(d["α"]),
        Float64(d["β"]),
        Float64(d["η"]),
        Float64(d["σ"]),
        Float64(d["μ"]),
        Float64(d["ϕ"]),
        Float64(d["γ"]),
        [κ0, κ0],                 # identical κ
        Vector{Float64}(d["a_by_occ"]),
        Dict(:f => Vector{Float64}(τw[:, 1]), :m => Vector{Float64}(τw[:, 2])),
        Dict(:f => Vector{Float64}(τe[:, 1]), :m => Vector{Float64}(τe[:, 2])),
        Dict(:f => 0.0, :m => 0.0),
        Dict(:f => 0.0, :m => 0.0),
        0.2,
        0.1,
        [0.0, 0.0],               # identical amenities
        [0.0 0.0; 0.0 0.0],       # zero moving costs
        0.2,
        0.7,
        0.4,
        0.6,
        9,
        5,
        0.2, 0.2, 0.2, 0.3,
        1e-4,
        500,
    )
    eqm_sym, _ = solve_stationary(p_sym)

    @test eqm_sym.M[1] ≈ eqm_sym.M[2] atol = 0.05
    @test eqm_sym.Htilde[1] ≈ eqm_sym.Htilde[2] atol = 0.05
    @test eqm_sym.t[1] ≈ eqm_sym.t[2] atol = 0.05
    @test sum(eqm_sym.m_young[1, :]) ≈ sum(eqm_sym.m_young[2, :]) atol = 0.05
end
