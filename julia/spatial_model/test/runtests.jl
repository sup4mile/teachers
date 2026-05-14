#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "..", "spatial.jl"))

function make_small_params(; τmove = [0.0 0.2; 0.2 0.0], B = [0.0, 0.0], κ = [1.0, 1.05], t = [0.2, 0.25])
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
        1.1,   # A_O
        Dict(:f => 0.1, :m => 0.0),
        Dict(:f => 0.05, :m => 0.0),
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
        0.2, 0.2, 0.2, 0.3,
        1e-5,
        300
    )
    xz, Pz = tauchen(p.Nz, p.ρz, p.σξ)
    z = exp.(xz)
    eps_dist = LogNormal(-0.5 * p.σeps^2, p.σeps)
    eps_nodes, eps_w = quantile_histogram(eps_dist, p.Neps)
    grids = Grids(xz, z, Pz, eps_nodes, eps_w)
    m = fill(1.0 / (p.L * p.Nz), p.L, p.Nz)
    eqm = Eqm(m, [1.0, 1.0], [1.0, 1.0], t, zeros(p.Nz, p.L, 2), zeros(p.Nz, p.L, 2))
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

@testset "Location-choice logic" begin
    p, grids, eqm = make_small_params()
    child_logh = zeros(p.Nz, p.L)
    child_cost = zeros(p.Nz, p.L)
    k = 2

    _, πO, hO, eO, _, _ = solve_O(p, grids, eqm, child_logh, child_cost, 1, k, :f, grids.eps_nodes[3])
    @test πO[1] + πO[2] ≈ 1.0 atol = 1e-10
    @test 0.0 <= πO[1] <= 1.0
    @test hO > 0.0
    @test eO > 0.0

    _, πT, hT, eT, _, _ = solve_T(p, grids, eqm, child_logh, child_cost, 1, k, :m, grids.eps_nodes[3])
    @test πT[1] + πT[2] ≈ 1.0 atol = 1e-10
    @test 0.0 <= πT[1] <= 1.0
    @test hT > 0.0
    @test eT > 0.0

    p_sym, grids_sym, eqm_sym = make_small_params(τmove = [0.0 0.0; 0.0 0.0], B = [0.0, 0.0], κ = [1.0, 1.0], t = [0.2, 0.2])
    _, π_sym, _, _, _, _ = solve_O(p_sym, grids_sym, eqm_sym, child_logh, child_cost, 1, k, :m, grids_sym.eps_nodes[3])
    @test π_sym[1] ≈ 0.5 atol = 1e-6

    p_sticky, grids_sticky, eqm_sticky = make_small_params(τmove = [0.0 8.0; 8.0 0.0], B = [0.0, 0.0], κ = [1.0, 1.0], t = [0.2, 0.2])
    _, π_sticky, _, _, _, _ = solve_O(p_sticky, grids_sticky, eqm_sticky, child_logh, child_cost, 1, k, :m, grids_sticky.eps_nodes[3])
    @test π_sticky[1] > 0.999
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
        mean(Vector{Float64}(d["a_by_occ"])),
        Dict(:f => mean(τw[:, 1]), :m => mean(τw[:, 2])),
        Dict(:f => mean(τe[:, 1]), :m => mean(τe[:, 2])),
        0.2,
        0.1,
        [0.0, 0.0],
        [0.0 0.2; 0.2 0.0],
        0.2,
        0.7,
        0.4,
        0.6,
        9,   # Neps
        5,   # Nz
        0.2, 0.2, 0.2, 0.3,
        1e-4,
        500
    )

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
        _, _, _, mig, Hc, tb, wb = state_moments(p, grids, eqm, child_logh, child_cost, l0, k, gsym)
        @test mig[1] + mig[2] ≈ 1.0 atol = 1e-8
        @test (0.0 <= mig[1] <= 1.0) && (0.0 <= mig[2] <= 1.0)
        wgt = 0.5 * eqm.m_young[l0, k]
        for ldest in 1:p.L
            H[ldest] += wgt * Hc[ldest]
            I[ldest] += wgt * tb[ldest]
            W[ldest] += wgt * wb[ldest]
        end
    end

    @test maximum(abs.(H .- eqm.Htilde)) < 5e-4
    @test maximum(abs.(eqm.t .* I .- W) ./ max.(abs.(W), 1e-12)) < 5e-3
end
