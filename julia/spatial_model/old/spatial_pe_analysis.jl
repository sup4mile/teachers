using Printf
using Plots

include(joinpath(@__DIR__, "spatial.jl"))

gr()

function fixed_point_error(S, E, H, V, p::Params, grids)
    S1, E1, H1, V1 = copy(S), copy(E), copy(H), copy(V)
    Dc, Hc = child_objects(E, H, V, p, grids)
    update_policy!(S1, E1, H1, V1, Dc, Hc, p, grids)
    return max(maximum(abs, E1 .- E), maximum(abs, V1 .- V))
end

function monotone_rate(A, grids; tol = 1e-8, either = false)
    Nz, Nϵ = grids.Nz, grids.Nϵ
    total_z = 2 * 2 * 2 * Nϵ
    total_ϵ = 2 * 2 * 2 * Nz
    pass_z = 0
    pass_ϵ = 0
    for i in 1:2, g in 1:2, l in 1:2
        for ei in 1:Nϵ
            v = vec(A[i, g, l, :, ei])
            d = diff(v)
            (all(d .>= -tol) || (either && all(d .<= tol))) && (pass_z += 1)
        end
        for zi in 1:Nz
            v = vec(A[i, g, l, zi, :])
            d = diff(v)
            (all(d .>= -tol) || (either && all(d .<= tol))) && (pass_ϵ += 1)
        end
    end
    return pass_z, total_z, pass_ϵ, total_ϵ
end

function report_policy_checks(S, E, V, grids)
    s_ok = all(isfinite.(S)) && all((S .> 0) .& (S .< 1))
    e_ok = all(isfinite.(E)) && all(E .> 0)
    v_ok = all(isfinite.(V))

    println("\nPolicy bounds / finiteness:")
    @printf("  s in (0,1): %s\n", s_ok ? "ok" : "FAIL")
    @printf("  e > 0     : %s\n", e_ok ? "ok" : "FAIL")
    @printf("  V finite  : %s\n", v_ok ? "ok" : "FAIL")

    println("\nMonotonicity checks (share monotone):")
    ps_z, ts_z, ps_ϵ, ts_ϵ = monotone_rate(S, grids; either = true)
    @printf("  s monotone in z: %3d/%3d\n", ps_z, ts_z)
    @printf("  s monotone in ϵ: %3d/%3d\n", ps_ϵ, ts_ϵ)

    pe_z, te_z, pe_ϵ, te_ϵ = monotone_rate(E, grids)
    @printf("  e increasing in z: %3d/%3d\n", pe_z, te_z)
    @printf("  e increasing in ϵ: %3d/%3d\n", pe_ϵ, te_ϵ)

    pv_z, tv_z, pv_ϵ, tv_ϵ = monotone_rate(V, grids)
    @printf("  V increasing in z: %3d/%3d\n", pv_z, tv_z)
    @printf("  V increasing in ϵ: %3d/%3d\n", pv_ϵ, tv_ϵ)
end

function plot_key_outcomes(S, E, H, V, grids; outdir)
    mkpath(outdir)
    z, ϵ = grids.z, grids.ϵgrid
    g, l = 1, 1
    zmed, ϵmed = (grids.Nz + 1) ÷ 2, (grids.Nϵ + 1) ÷ 2

    p1 = plot(layout = (2, 2), size = (900, 650))
    plot!(p1[1], z, S[1, g, l, :, ϵmed], label = "Teacher", xlabel = "z", ylabel = "s")
    plot!(p1[1], z, S[2, g, l, :, ϵmed], label = "Non-teacher")
    plot!(p1[2], z, E[1, g, l, :, ϵmed], label = "Teacher", xlabel = "z", ylabel = "e")
    plot!(p1[2], z, E[2, g, l, :, ϵmed], label = "Non-teacher")
    plot!(p1[3], z, H[1, g, l, :, ϵmed], label = "Teacher", xlabel = "z", ylabel = "h")
    plot!(p1[3], z, H[2, g, l, :, ϵmed], label = "Non-teacher")
    plot!(p1[4], z, V[1, g, l, :, ϵmed], label = "Teacher", xlabel = "z", ylabel = "V")
    plot!(p1[4], z, V[2, g, l, :, ϵmed], label = "Non-teacher")
    savefig(p1, joinpath(outdir, "policies_values_vs_z.png"))

    p2 = plot(layout = (2, 2), size = (900, 650))
    plot!(p2[1], ϵ, S[1, g, l, zmed, :], label = "Teacher", xlabel = "ϵ", ylabel = "s")
    plot!(p2[1], ϵ, S[2, g, l, zmed, :], label = "Non-teacher")
    plot!(p2[2], ϵ, E[1, g, l, zmed, :], label = "Teacher", xlabel = "ϵ", ylabel = "e")
    plot!(p2[2], ϵ, E[2, g, l, zmed, :], label = "Non-teacher")
    plot!(p2[3], ϵ, H[1, g, l, zmed, :], label = "Teacher", xlabel = "ϵ", ylabel = "h")
    plot!(p2[3], ϵ, H[2, g, l, zmed, :], label = "Non-teacher")
    plot!(p2[4], ϵ, V[1, g, l, zmed, :], label = "Teacher", xlabel = "ϵ", ylabel = "V")
    plot!(p2[4], ϵ, V[2, g, l, zmed, :], label = "Non-teacher")
    savefig(p2, joinpath(outdir, "policies_values_vs_eps.png"))
end

function sweep_parameterizations()
    println("\nParameter sweep (ε = 0 throughout):")
    BASE_GRIDS = (H̃T = [1.0, 1.5], M = [1.0, 1.0], t = [0.10, 0.12])
    cases = [
        # --- GE aggregates / shock dispersion ---
        (label = "baseline-ε0",      p = Params(ε = 0.0),
         grids = BASE_GRIDS),
        (label = "high-tax",         p = Params(ε = 0.0),
         grids = (H̃T = [1.0, 1.5], M = [1.0, 1.0], t = [0.25, 0.30])),
        (label = "low-Q",            p = Params(ε = 0.0),
         grids = (H̃T = [0.6, 0.6], M = [1.0, 1.0], t = [0.10, 0.12])),
        (label = "high-σϵ",          p = Params(ε = 0.0, σϵ = 0.60),
         grids = BASE_GRIDS),
        # --- Altruism ---
        (label = "high-altruism",    p = Params(ε = 0.0, λ = 0.80),
         grids = BASE_GRIDS),
        (label = "low-altruism",     p = Params(ε = 0.0, λ = 0.20),
         grids = BASE_GRIDS),
        # --- Teaching wage schedule ---
        (label = "high-teach-wage",  p = Params(ε = 0.0, κ = [1.20, 1.40]),
         grids = BASE_GRIDS),
        (label = "high-γ",           p = Params(ε = 0.0, γ = 0.95),
         grids = BASE_GRIDS),
        # --- Gender / occupation distortions ---
        (label = "gender-wedge",     p = Params(ε = 0.0, τω = [0.0 0.0; 0.0 0.20]),
         grids = BASE_GRIDS),
        # --- Spatial frictions ---
        (label = "high-move-cost",   p = Params(ε = 0.0, τmove = [0.0 0.60; 0.60 0.0]),
         grids = BASE_GRIDS),
        (label = "low-mobility",     p = Params(ε = 0.0, σν = 0.25),
         grids = BASE_GRIDS),
        # --- Ability process ---
        (label = "high-persistence", p = Params(ε = 0.0, ρz = 0.90),
         grids = BASE_GRIDS),
        (label = "low-persistence",  p = Params(ε = 0.0, ρz = 0.40),
         grids = BASE_GRIDS),
    ]

    for c in cases
        grids = build_grids(c.p; c.grids..., Nz = 3, Nϵ = 5)
        S, E, H, V = solve_household(c.p, grids; tol = 1e-4, maxit = 150)
        err = fixed_point_error(S, E, H, V, c.p, grids)
        s_ok = all(isfinite.(S)) && all((S .> 0) .& (S .< 1))
        e_ok = all(isfinite.(E)) && all(E .> 0)
        v_ok = all(isfinite.(V))
        @printf("  %-16s  err=%.2e  bounds=%s\n", c.label, err,
                (s_ok && e_ok && v_ok) ? "ok" : "FAIL")
    end
end

function main()
    p = Params()
    grids = build_grids(p; H̃T = [1.0, 1.5], M = [1.0, 1.0], t = [0.10, 0.12],
                        Nz = 7, Nϵ = 21)
    println("Solving baseline partial-equilibrium household fixed point ...")
    S, E, H, V = solve_household(p, grids)
    report(S, E, H, V, p, grids)

    report_policy_checks(S, E, V, grids)
    plot_key_outcomes(S, E, H, V, grids; outdir = joinpath(@__DIR__, "figures"))
    sweep_parameterizations()
end

@time main()
