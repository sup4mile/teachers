# =============================================================================
# Diagnostics & tests for the GENERAL-EQUILIBRIUM spatial solution (spatial_ge.jl).
#
# Mirrors spatial_pe_analysis.jl, but where the PE script checks the *household*
# fixed point taking aggregates as given, this script checks the *closed* economy:
# the stationary equilibrium in (H̃_T, M, t, Φ).  It verifies
#
#   1. Internal consistency of a converged equilibrium:
#        • GE-map residual          — feeding the aggregates back reproduces them
#        • Φ stationarity           — Φ is the fixed point of migration + ability evolution
#        • mass conservation        — Σ_l M_l = Mtot, M_l = 2∫Φ_l
#        • row-stochasticity        — π_{l'|l} and π̄_{l₀→l'} rows sum to 1
#        • Q consistency            — Q_l = (2H̃_T,l/M_l)^σ
#        • household fixed point     — one extra policy sweep barely moves (E,V)
#        • policy bounds / shape    — s∈(0,1), e>0, V finite, monotone in (z,ϵ)
#   2. Economic-sense properties via controlled experiments:
#        • symmetry   — identical locations ⇒ symmetric equilibrium
#        • scale      — doubling Mtot doubles (M,H̃_T), leaves (Q,t,Gₗ) fixed
#        • comparative statics across a battery of parameter changes (signs sane,
#          policies in-bounds; a `conv?` flag marks the few hard cases — large β,
#          σν→0 — that the simple damped iteration cannot settle in the budget)
#        • damping diagnostic — the outer iteration is only weakly contracting;
#          this reports the residual reached at several damping values.
#
# Run:  julia --project=.. spatial_ge_analysis.jl
# =============================================================================

using Printf
using Plots

include(joinpath(@__DIR__, "spatial_ge.jl"))

gr()

# -----------------------------------------------------------------------------
# 0. Generic helpers (re-implemented here so this file does not pull in
#    spatial_pe_analysis.jl, whose top level runs the PE example on include).
# -----------------------------------------------------------------------------

"Fraction of (i,g,l) slices along which A is monotone in z and in ϵ.  `either`
 accepts non-decreasing OR non-increasing (used for s, whose sign can flip)."
function monotone_rate(A, grids; tol = 1e-8, either = false)
    Nz, Nϵ = grids.Nz, grids.Nϵ
    total_z, total_ϵ = 2 * 2 * 2 * Nϵ, 2 * 2 * 2 * Nz
    pass_z = pass_ϵ = 0
    for i in 1:2, g in 1:2, l in 1:2
        for ei in 1:Nϵ
            d = diff(vec(A[i, g, l, :, ei]))
            (all(d .>= -tol) || (either && all(d .<= tol))) && (pass_z += 1)
        end
        for zi in 1:Nz
            d = diff(vec(A[i, g, l, zi, :]))
            (all(d .>= -tol) || (either && all(d .<= tol))) && (pass_ϵ += 1)
        end
    end
    return pass_z, total_z, pass_ϵ, total_ϵ
end

policies_in_bounds(S, E, V) =
    all(isfinite, S) && all((S .> 0) .& (S .< 1)) &&
    all(isfinite, E) && all(E .> 0) && all(isfinite, V)

mark(ok) = ok ? "ok" : "FAIL"

# -----------------------------------------------------------------------------
# 1. Equilibrium-consistency residuals.
#
#    A converged `sol` should satisfy each equilibrium condition exactly (up to
#    the tolerance of the outer iteration).  Each function returns a residual we
#    can threshold; structural identities (row sums, mass, Q) hold to ~1e-12,
#    while the GE-map / household residuals inherit the outer/inner tolerances.
# -----------------------------------------------------------------------------

"Residual of the aggregate fixed point: re-evaluate the GE map at `sol` and
 compare to the stored aggregates.  This is the *honest* convergence measure —
 the `err` printed inside solve_ge can stall in a small limit cycle at high
 damping, but this residual is zero only at a true fixed point."
function ge_map_residual(sol)
    HTn, Mn, tn = aggregates(sol.Φ, sol.Pi, sol.H, sol.V, sol.gp, sol.grids)
    rH = maximum(abs, (HTn .- sol.HT) ./ sol.HT)
    rM = maximum(abs, (Mn  .- sol.M ) ./ sol.M )
    rt = maximum(abs,  tn  .- sol.t)
    return (; res = max(rH, rM, rt), rH, rM, rt, HTn, Mn, tn)
end

"Residual of the Φ eigen-problem (A.5): rebuild the joint kernel K from π̄ and
 Πz and check that vec(Φ) is its dominant left eigenvector with eigenvalue 1.
 Also returns the kernel row-sum deviation (mass conservation of K) and Σ Φ."
function phi_stationarity_residual(sol)
    (; Nz, Πz) = sol.grids
    πbar = sol.πbar
    n = 2 * Nz
    idx(l, zi) = (l - 1) * Nz + zi
    P = zeros(n, n)
    for l0 in 1:2, zi in 1:Nz, l in 1:2, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]
    end
    φ = zeros(n)
    for l in 1:2, zi in 1:Nz
        φ[idx(l, zi)] = sol.Φ[zi, l]
    end
    eig_res     = maximum(abs, vec(φ' * P) .- φ) / sum(φ)
    row_dev     = maximum(abs, vec(sum(P, dims = 2)) .- 1)
    return (; eig_res, row_dev, mass = sum(φ))
end

"Worst deviation from 1 of the location-choice rows Σ_{l'} π_{l'|l}."
function pi_rowsum_dev(sol)
    (; Nz, Nϵ) = sol.grids
    dev = 0.0
    for i in 1:2, g in 1:2, l in 1:2, zi in 1:Nz, ei in 1:Nϵ
        dev = max(dev, abs(sum(@view sol.Pi[i, g, l, zi, ei, :]) - 1))
    end
    return dev
end

"Worst deviation from 1 of the migration rows Σ_{l'} π̄_{l₀→l'}(z)."
pibar_rowsum_dev(sol) =
    maximum(abs, [sum(@view sol.πbar[l0, zi, :]) - 1
                  for l0 in 1:2, zi in 1:sol.grids.Nz])

"One extra household policy sweep at the converged grids; the residual measures
 how close the stored (E,V) are to the household fixed point."
function household_residual(sol)
    (; S, E, H, V, grids) = sol
    p = sol.gp.p
    E1, V1 = copy(E), copy(V)
    S1, H1 = copy(S), copy(H)
    Dc, Hc = child_objects(E, H, V, p, grids)
    update_policy!(S1, E1, H1, V1, Dc, Hc, p, grids)
    return max(maximum(abs, E1 .- E), maximum(abs, V1 .- V))
end

"Consistency of the teacher-quality shifter:  Q_l = (2 H̃_T,l / M_l)^σ."
q_consistency_dev(sol) =
    maximum(abs, sol.grids.Q .- (2 .* sol.HT ./ sol.M) .^ sol.gp.p.σ)

function report_ge_checks(sol)
    println("\n========== Equilibrium-consistency checks ==========")

    # (a) Structural identities — hold to machine precision at *any* iterate,
    #     converged or not, because they are built into the GE map itself.
    println("  Structural identities (should be ~0 always):")
    ph = phi_stationarity_residual(sol)
    @printf("    Φ eigen-residual             : %.2e   [%s]\n", ph.eig_res, mark(ph.eig_res < 1e-8))
    @printf("    K row-sum deviation          : %.2e   [%s]\n", ph.row_dev, mark(ph.row_dev < 1e-10))
    mass_dev = abs(ph.mass - sol.gp.Mtot / 2)
    Mtot_dev = abs(sum(sol.M) - sol.gp.Mtot)
    @printf("    Σ Φ = Mtot/2                 : %.2e   [%s]\n", mass_dev, mark(mass_dev < 1e-8))
    @printf("    Σ_l M_l = Mtot               : %.2e   [%s]\n", Mtot_dev, mark(Mtot_dev < 1e-8))
    pir = pi_rowsum_dev(sol)
    pbr = pibar_rowsum_dev(sol)
    @printf("    Σ π_{l'|l} = 1               : %.2e   [%s]\n", pir, mark(pir < 1e-10))
    @printf("    Σ π̄_{l₀→l'} = 1              : %.2e   [%s]\n", pbr, mark(pbr < 1e-10))
    qd = q_consistency_dev(sol)
    @printf("    Q_l = (2H̃_T,l/M_l)^σ         : %.2e   [%s]\n", qd, mark(qd < 1e-10))

    # (b) Equilibrium conditions that close only at a true fixed point; their
    #     residuals inherit the outer/inner iteration tolerances.
    gm = ge_map_residual(sol)
    Ml_dev = maximum(abs, sol.M .- [2 * sum(@view sol.Φ[:, l]) for l in 1:2])
    hh = household_residual(sol)
    println("\n  Convergence-sensitive equilibrium conditions:")
    @printf("    GE-map residual              : %.2e   (H̃_T %.1e, M %.1e, t %.1e)  [%s]\n",
            gm.res, gm.rH, gm.rM, gm.rt, mark(gm.res < 5e-3))
    @printf("    M_l = 2 ∫ Φ_l                : %.2e   [%s]\n", Ml_dev, mark(Ml_dev < 5e-3))
    @printf("    household policy residual    : %.2e   [%s]\n", hh, mark(hh < 1e-3))
    gm.res < 5e-3 || println("    ⚠ outer iteration not tightly converged — " *
                             "lower `damping` and/or raise `maxit`.")

    # Budget: solve_ge clamps t to (0,0.9); if a rate hits the clamp the local
    # balanced budget t·I = W does NOT hold exactly.  Flag it.
    clamp_hit = any(sol.t .>= 0.9 - 1e-9) || any(sol.t .<= 1e-9)
    @printf("  budget clamp inactive          : %s   (t = %.3f, %.3f)\n",
            mark(!clamp_hit), sol.t[1], sol.t[2])

    println("\n  Policy bounds / finiteness:")
    S, E, V, grids = sol.S, sol.E, sol.V, sol.grids
    @printf("    s∈(0,1) : %s    e>0 : %s    V finite : %s\n",
            mark(all((S .> 0) .& (S .< 1))), mark(all(E .> 0)), mark(all(isfinite, V)))

    println("\n  Monotonicity (household policies):")
    ps_z, ts_z, ps_ϵ, ts_ϵ = monotone_rate(S, grids; either = true)
    pe_z, te_z, pe_ϵ, te_ϵ = monotone_rate(E, grids)
    pv_z, tv_z, pv_ϵ, tv_ϵ = monotone_rate(V, grids)
    @printf("    s monotone in z: %3d/%-3d  in ϵ: %3d/%-3d\n", ps_z, ts_z, ps_ϵ, ts_ϵ)
    @printf("    e increasing in z: %3d/%-3d in ϵ: %3d/%-3d\n", pe_z, te_z, pe_ϵ, te_ϵ)
    @printf("    V increasing in z: %3d/%-3d in ϵ: %3d/%-3d\n", pv_z, tv_z, pv_ϵ, tv_ϵ)
    println("====================================================")
end

# -----------------------------------------------------------------------------
# 2. Symmetry test.
#
#    Make the two locations literally identical (B and κ symmetric; τmove, τω,
#    τe, A are already location-symmetric).  Because the solver starts from a
#    symmetric guess, the unique stable stationary equilibrium must be symmetric:
#    H̃_T, M, t equal across locations and G_1(z) = G_2(z).
# -----------------------------------------------------------------------------
function symmetry_test(; Nz = 3, Nϵ = 5, damping = 0.2, maxit = 80)
    println("\n========== Symmetry test (identical locations) ==========")
    p  = Params(B = [0.0, 0.0], κ = [1.0, 1.0])
    sol = solve_ge(GEParams(p = p); Nz, Nϵ, damping, tol = 1e-3, maxit,
                   hh_tol = 1e-5, hh_maxit = 40, verbose = false)

    g1 = sol.Φ[:, 1] ./ sum(@view sol.Φ[:, 1])
    g2 = sol.Φ[:, 2] ./ sum(@view sol.Φ[:, 2])
    aH, aM, at, aG = abs(sol.HT[1] - sol.HT[2]), abs(sol.M[1] - sol.M[2]),
                     abs(sol.t[1] - sol.t[2]), maximum(abs, g1 .- g2)

    @printf("  H̃_T = (%.4f, %.4f)   |Δ| = %.2e\n", sol.HT[1], sol.HT[2], aH)
    @printf("  M   = (%.4f, %.4f)   |Δ| = %.2e\n", sol.M[1],  sol.M[2],  aM)
    @printf("  t   = (%.4f, %.4f)   |Δ| = %.2e\n", sol.t[1],  sol.t[2],  at)
    @printf("  max |G_1(z) - G_2(z)|            = %.2e\n", aG)
    ok = max(aH, aM, at, aG) < 1e-3
    @printf("  symmetric equilibrium recovered : %s\n", mark(ok))
    println("=========================================================")
    return ok
end

# -----------------------------------------------------------------------------
# 3. Scale-invariance test.
#
#    Q_l depends on H̃_T,l/M_l only, and Φ rescales linearly with Mtot, so
#    doubling the total measure must double (M, H̃_T) while leaving the intensive
#    objects (Q, t, and the shape Gₗ) untouched.
# -----------------------------------------------------------------------------
function scale_invariance_test(; Nz = 3, Nϵ = 5, damping = 0.2, maxit = 80)
    println("\n========== Scale-invariance test (Mtot: 2 → 4) ==========")
    kw = (; Nz, Nϵ, damping, tol = 1e-3, maxit, hh_tol = 1e-5, hh_maxit = 40, verbose = false)
    s1 = solve_ge(GEParams(Mtot = 2.0); kw...)
    s2 = solve_ge(GEParams(Mtot = 4.0); kw...)

    dM  = maximum(abs, s2.M  .- 2 .* s1.M)               # extensive: should double
    dHT = maximum(abs, s2.HT .- 2 .* s1.HT)
    dQ  = maximum(abs, s2.grids.Q .- s1.grids.Q)         # intensive: should match
    dt  = maximum(abs, s2.t .- s1.t)
    g(s, l) = s.Φ[:, l] ./ sum(@view s.Φ[:, l])
    dG  = maximum(maximum(abs, g(s2, l) .- g(s1, l)) for l in 1:2)

    @printf("  M(4) vs 2·M(2)   : Δ = %.2e   [%s]\n", dM,  mark(dM  < 5e-3))
    @printf("  H̃_T(4) vs 2·H̃_T(2): Δ = %.2e   [%s]\n", dHT, mark(dHT < 5e-3))
    @printf("  Q  unchanged     : Δ = %.2e   [%s]\n", dQ,  mark(dQ  < 5e-3))
    @printf("  t  unchanged     : Δ = %.2e   [%s]\n", dt,  mark(dt  < 5e-3))
    @printf("  Gₗ unchanged     : Δ = %.2e   [%s]\n", dG,  mark(dG  < 5e-3))
    println("=========================================================")
    return max(dM, dHT, dQ, dt, dG) < 5e-3
end

# -----------------------------------------------------------------------------
# 4. Comparative statics.
#
#    Each case is a one-parameter deviation from the baseline.  We solve the
#    full GE and report the converged aggregates, the GE-map residual, whether
#    policies stay in bounds, and two summary statistics:
#       outflow  = Φ-weighted average prob. of leaving the birth location,
#       teach    = population-weighted teaching share (endogenous Gₗ).
# -----------------------------------------------------------------------------
"Φ-weighted gross outflow rate Σ_{l₀} Σ_z Φ[z,l₀] π̄_{l₀→l'≠l₀}(z) / Σ Φ."
function outflow_rate(sol)
    (; Nz) = sol.grids
    num = den = 0.0
    for l0 in 1:2, zi in 1:Nz
        w = sol.Φ[zi, l0]
        num += w * sol.πbar[l0, zi, 3 - l0]
        den += w
    end
    return num / den
end

"Population-weighted teaching share across both birth locations (endogenous Gₗ)."
function teach_share_total(sol)
    num = sum(teaching_share_endo(sol.V, l, sol.Φ, sol.grids) * sum(@view sol.Φ[:, l]) for l in 1:2)
    return num / sum(sol.Φ)
end

function comparative_statics(; Nz = 3, Nϵ = 5, damping = 0.2, maxit = 120, hh_maxit = 120)
    println("\n========== Comparative statics (one-parameter deviations) ==========")
    cases = [
        (label = "baseline",        gp = GEParams()),
        # Teacher-spillover strength: larger β amplifies the H̃_T feedback.
        (label = "high-β",          gp = GEParams(β = 0.50)),
        (label = "zero-β",          gp = GEParams(β = 0.0)),
        # Altruism.
        (label = "high-altruism",   gp = GEParams(p = Params(λ = 0.90))),
        (label = "low-altruism",    gp = GEParams(p = Params(λ = 0.10))),
        # Spatial frictions: higher move cost / lower taste noise ⇒ less mobility.
        (label = "high-move-cost",  gp = GEParams(p = Params(τmove = [0.0 1.50; 1.50 0.0]))),
        (label = "low-σν",          gp = GEParams(p = Params(σν = 0.20))),
        # Amenity pull toward location 2 ⇒ M₂ should rise vs baseline.
        (label = "amenity-loc2",    gp = GEParams(p = Params(B = [0.0, 0.50]))),
        # Teaching wage schedule.
        (label = "high-teach-wage", gp = GEParams(p = Params(κ = [1.20, 1.40]))),
        # Gender wage wedge on women in non-teaching.
        (label = "gender-wedge",    gp = GEParams(p = Params(τω = [0.0 0.0; 0.0 0.30]))),
        # Ability persistence.
        (label = "high-persistence",gp = GEParams(p = Params(ρz = 0.90))),
        (label = "low-persistence", gp = GEParams(p = Params(ρz = 0.40))),
    ]

    @printf("  %-16s %-9s %-6s %-19s %-19s %-7s %-7s %s\n",
            "case", "res", "conv?", "H̃_T", "M", "out", "teach", "bounds")
    sols  = Dict{String,Any}()
    stuck = String[]
    for c in cases
        sol = solve_ge(c.gp; Nz, Nϵ, damping, tol = 1e-3, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        sols[c.label] = sol
        res  = ge_map_residual(sol).res
        conv = res < 1e-2                         # treat >1e-2 as "did not settle"
        conv || push!(stuck, c.label)
        ok   = policies_in_bounds(sol.S, sol.E, sol.V)
        @printf("  %-16s %.1e  %-6s (%.3f,%.3f)       (%.3f,%.3f)      %.3f  %.3f  %s\n",
                c.label, res, conv ? "yes" : "NO", sol.HT[1], sol.HT[2],
                sol.M[1], sol.M[2], outflow_rate(sol), teach_share_total(sol), mark(ok))
    end
    if !isempty(stuck)
        println("\n  ⚠ did not converge in the test budget (read their aggregates with care): ",
                join(stuck, ", "))
        println("    expected: high β weakens the contraction, and σν→0 / large move costs")
        println("    make the Φ-eigenproblem near-degenerate.  The household block is")
        println("    already solved to the headline inner budget here, so a remaining")
        println("    stall is in the *outer* map: re-solve with a smaller `damping`")
        println("    (see the damping diagnostic) and a larger `maxit`.")
    end

    # Directional sanity checks against the baseline.
    println("\n  Sign checks:")
    base = sols["baseline"]
    chk(name, cond) = @printf("    %-40s %s\n", name, cond ? "ok" : "FAIL")
    chk("amenity in loc 2 raises M₂", sols["amenity-loc2"].M[2] > base.M[2])
    chk("high move cost lowers outflow", outflow_rate(sols["high-move-cost"]) < outflow_rate(base))
    chk("low σν lowers outflow",        outflow_rate(sols["low-σν"]) < outflow_rate(base))
    chk("every case in bounds",
        all(policies_in_bounds(s.S, s.E, s.V) for s in values(sols)))
    println("====================================================================")
    return sols
end

# -----------------------------------------------------------------------------
# 5. Damping diagnostic.
#
#    The outer map is only weakly contracting (its H̃_T-elasticity is ≈ β), so the
#    successive-relaxation step size matters.  Report the residual reached at a
#    few damping values to guide the choice for production runs.
# -----------------------------------------------------------------------------
function damping_diagnostic(; Nz = 3, Nϵ = 15, maxit = 150, hh_maxit = 120)
    println("\n========== Damping diagnostic (residual vs step size) ==========")
    println("  The outer GE map is only weakly contracting (its H̃_T-elasticity")
    println("  is ≈ β), so the successive-relaxation step size matters.  To keep")
    println("  the comparison fair, this re-solves the baseline at the headline")
    @printf("  grid and inner budget (Nϵ=%d, hh_maxit=%d), holds the outer\n", Nϵ, hh_maxit)
    @printf("  budget fixed at %d iterations, and varies only `damping`.\n", maxit)
    @printf("  %-9s %-12s %s\n", "damping", "GE residual", "usable(<5e-3)")
    results = Tuple{Float64,Float64}[]
    for d in (0.50, 0.30, 0.20)
        sol = solve_ge(GEParams(); Nz, Nϵ, damping = d, tol = 1e-4, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        res = ge_map_residual(sol).res
        push!(results, (d, res))
        @printf("  %-9.2f %.3e    %s\n", d, res, res < 5e-3 ? "yes" : "no")
    end
    _, i = findmin(last.(results))
    best_d, best_res = results[i]
    @printf("  lowest residual in this budget: damping = %.2f (%.2e).\n", best_d, best_res)
    println("  The headline solve uses damping = 0.2; under-relaxation (smaller")
    println("  steps) costs iterations but is more robust for a weakly contracting")
    println("  map, which is why the production solve favours a smaller step.")
    println("================================================================")
end

# -----------------------------------------------------------------------------
# 6. Plots.
# -----------------------------------------------------------------------------
function plot_ge_outcomes(sol; outdir)
    mkpath(outdir)
    (; z, ϵgrid, Nz, Nϵ, Πz) = sol.grids
    S, E, H, V = sol.S, sol.E, sol.H, sol.V
    g = 1
    zmed, ϵmed = (Nz + 1) ÷ 2, (Nϵ + 1) ÷ 2

    # (a) Endogenous ability distribution by location vs the ergodic law.
    erg = stationary(Πz)
    g1  = sol.Φ[:, 1] ./ sum(@view sol.Φ[:, 1])
    g2  = sol.Φ[:, 2] ./ sum(@view sol.Φ[:, 2])
    pd  = plot(z, g1, marker = :circle, label = "G₁(z)", xlabel = "z",
               ylabel = "mass", title = "Endogenous ability distribution")
    plot!(pd, z, g2, marker = :square,  label = "G₂(z)")
    plot!(pd, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pd, joinpath(outdir, "ge_ability_distribution.png"))

    # (b) Policies / value vs z (location 1, men, median ϵ).
    p1 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        plot!(p1[idx], z, A[1, g, 1, :, ϵmed], label = "Teacher",     xlabel = "z", ylabel = lab)
        plot!(p1[idx], z, A[2, g, 1, :, ϵmed], label = "Non-teacher")
    end
    savefig(p1, joinpath(outdir, "ge_policies_vs_z.png"))

    # (c) Policies / value vs ϵ (location 1, men, median z).
    p2 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        plot!(p2[idx], ϵgrid, A[1, g, 1, zmed, :], label = "Teacher",     xlabel = "ϵ", ylabel = lab)
        plot!(p2[idx], ϵgrid, A[2, g, 1, zmed, :], label = "Non-teacher")
    end
    savefig(p2, joinpath(outdir, "ge_policies_vs_eps.png"))

    # (d) Migration matrix Π̄[born → work] as a heatmap (base Plots only).
    Π̄ = [sum(sol.Φ[zi, l0] * sol.πbar[l0, zi, l] for zi in 1:Nz) / sum(@view sol.Φ[:, l0])
         for l0 in 1:2, l in 1:2]
    pm = heatmap(["work 1", "work 2"], ["born 1", "born 2"], Π̄,
                 clims = (0, 1), c = :viridis, yflip = true,
                 title = "Migration matrix Π̄[born→work]")
    annotate!(pm, vec([(j, i, text(@sprintf("%.3f", Π̄[i, j]), 9, :white))
                       for i in 1:2, j in 1:2]))
    savefig(pm, joinpath(outdir, "ge_migration_matrix.png"))

    println("\nSaved GE figures to ", outdir)
end

# -----------------------------------------------------------------------------
# 7. Run.   The headline solve (Nz=3, Nϵ=15) plus the symmetry, scale,
#    comparative-statics, and damping batteries.  Each GE solve nests a household
#    fixed point inside every outer iteration, so cost scales with grid size ×
#    outer iters × inner sweeps.  To keep the convergence flags *fair*, the
#    comparative-statics and damping batteries solve the household block to the
#    same inner budget as the headline (hh_maxit=120) instead of a truncated one:
#    this is slower, but a remaining "did not converge" flag then reflects the
#    outer GE map, not an under-solved household problem.  (Expect the run to take
#    noticeably longer than a truncated-budget pass for this reason.)
# -----------------------------------------------------------------------------
function main_test()
    println("Solving spatial model in general equilibrium (baseline) ...")
    sol = solve_ge(GEParams(); Nz = 3, Nϵ = 15, damping = 0.2, tol = 1e-4,
                   maxit = 150, hh_tol = 1e-5, hh_maxit = 120, verbose = true)
    report_ge(sol)
    report_ge_checks(sol)
    plot_ge_outcomes(sol; outdir = joinpath(@__DIR__, "figures"))

    symmetry_test()
    scale_invariance_test()
    comparative_statics()
    damping_diagnostic()
    return sol
end

@time main_test()
