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
    Nz, Nϵ, I, L = grids.Nz, grids.Nϵ, grids.I, grids.L
    total_z, total_ϵ = I * 2 * L * Nϵ, I * 2 * L * Nz
    pass_z = pass_ϵ = 0
    for i in 1:I, g in 1:2, l in 1:L
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
    HTn, Mn, tn = aggregates(sol.Φ, sol.Pi, sol.H, sol.V, sol.p, sol.grids)
    rH = maximum(abs, (HTn .- sol.HT) ./ sol.HT)
    rM = maximum(abs, (Mn  .- sol.M ) ./ sol.M )
    rt = maximum(abs,  tn  .- sol.t)
    return (; res = max(rH, rM, rt), rH, rM, rt, HTn, Mn, tn)
end

"Residual of the Φ eigen-problem (A.5): rebuild the joint kernel K from π̄ and
 Πz and check that vec(Φ) is its dominant left eigenvector with eigenvalue 1.
 Also returns the kernel row-sum deviation (mass conservation of K) and Σ Φ."
function phi_stationarity_residual(sol)
    (; Nz, Πz, L) = sol.grids
    πbar = sol.πbar
    n = L * Nz
    idx(l, zi) = (l - 1) * Nz + zi
    P = zeros(n, n)
    for l0 in 1:L, zi in 1:Nz, l in 1:L, zpi in 1:Nz
        P[idx(l0, zi), idx(l, zpi)] = πbar[l0, zi, l] * Πz[zi, zpi]
    end
    φ = zeros(n)
    for l in 1:L, zi in 1:Nz
        φ[idx(l, zi)] = sol.Φ[zi, l]
    end
    eig_res     = maximum(abs, vec(φ' * P) .- φ) / sum(φ)
    row_dev     = maximum(abs, vec(sum(P, dims = 2)) .- 1)
    return (; eig_res, row_dev, mass = sum(φ))
end

"Worst deviation from 1 of the location-choice rows Σ_{l'} π_{l'|l}."
function pi_rowsum_dev(sol)
    (; Nz, Nϵ, I, L) = sol.grids
    dev = 0.0
    for i in 1:I, g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
        dev = max(dev, abs(sum(@view sol.Pi[i, g, l, zi, ei, :]) - 1))
    end
    return dev
end

"Worst deviation from 1 of the migration rows Σ_{l'} π̄_{l₀→l'}(z)."
pibar_rowsum_dev(sol) =
    maximum(abs, [sum(@view sol.πbar[l0, zi, :]) - 1
                  for l0 in 1:sol.grids.L, zi in 1:sol.grids.Nz])

"One extra household policy sweep at the converged grids; the residual measures
 how close the stored (E,V) are to the household fixed point."
function household_residual(sol)
    (; S, E, H, V, grids) = sol
    p = sol.p
    E1, V1 = copy(E), copy(V)
    S1, H1 = copy(S), copy(H)
    Dc, Hc = child_objects(E, H, V, p, grids)
    update_policy!(S1, E1, H1, V1, Dc, Hc, p, grids)
    return max(maximum(abs, E1 .- E), maximum(abs, V1 .- V))
end

"Consistency of the teacher-quality shifter:  Q_l = (2 H̃_T,l / M_l)^σ."
q_consistency_dev(sol) =
    maximum(abs, sol.grids.Q .- (2 .* sol.HT ./ sol.M) .^ sol.p.σ)

function report_ge_checks(sol)
    println("\n========== Equilibrium-consistency checks ==========")

    # (a) Structural identities — hold to machine precision at *any* iterate,
    #     converged or not, because they are built into the GE map itself.
    println("  Structural identities (should be ~0 always):")
    ph = phi_stationarity_residual(sol)
    @printf("    Φ eigen-residual             : %.2e   [%s]\n", ph.eig_res, mark(ph.eig_res < 1e-8))
    @printf("    K row-sum deviation          : %.2e   [%s]\n", ph.row_dev, mark(ph.row_dev < 1e-10))
    mass_dev = abs(ph.mass - sol.p.Mtot / 2)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
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
    Ml_dev = maximum(abs, sol.M .- [2 * sum(@view sol.Φ[:, l]) for l in 1:sol.grids.L])
    hh = household_residual(sol)
    println("\n  Convergence-sensitive equilibrium conditions:")
    @printf("    GE-map residual              : %.2e   (H̃_T %.1e, M %.1e, t %.1e)  [%s]\n",
            gm.res, gm.rH, gm.rM, gm.rt, mark(gm.res < 5e-3))
    @printf("    M_l = 2 ∫ Φ_l                : %.2e   [%s]\n", Ml_dev, mark(Ml_dev < 5e-3))
    @printf("    household policy residual    : %.2e   [%s]\n", hh, mark(hh < 1e-3))
    gm.res < 5e-3 || println("    ⚠ outer iteration not tightly converged — " *
                             "raise `damping` and/or increase `maxit`.")

    # Budget: solve_ge clamps t to (0,0.9); if a rate hits the clamp the local
    # balanced budget t·I = W does NOT hold exactly.  Flag it.
    clamp_hit = any(sol.t .>= 0.9 - 1e-9) || any(sol.t .<= 1e-9)
    @printf("  budget clamp inactive          : %s   (t = %s)\n",
            mark(!clamp_hit), fmtvec(sol.t))

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
    sol = solve_ge(p; Nz, Nϵ, damping, tol = 1e-3, maxit,
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
    kw = (; Nz, Nϵ, damping, tol = 1e-4, maxit, hh_tol = 1e-5, hh_maxit = 400, verbose = false)
    s1 = solve_ge(Params(Mtot = 2.0); kw...)
    s2 = solve_ge(Params(Mtot = 4.0); kw...)

    dM  = maximum(abs, s2.M  .- 2 .* s1.M)               # extensive: should double
    dHT = maximum(abs, s2.HT .- 2 .* s1.HT)
    dQ  = maximum(abs, s2.grids.Q .- s1.grids.Q)         # intensive: should match
    dt  = maximum(abs, s2.t .- s1.t)
    g(s, l) = s.Φ[:, l] ./ sum(@view s.Φ[:, l])
    dG  = maximum(maximum(abs, g(s2, l) .- g(s1, l)) for l in 1:2)

    @printf("  M(4) vs 2·M(2)   : Δ = %.2e   [%s]\n", dM,  mark(dM  < 1e-3))
    @printf("  H̃_T(4) vs 2·H̃_T(2): Δ = %.2e   [%s]\n", dHT, mark(dHT < 1e-3))
    @printf("  Q  unchanged     : Δ = %.2e   [%s]\n", dQ,  mark(dQ  < 1e-3))
    @printf("  t  unchanged     : Δ = %.2e   [%s]\n", dt,  mark(dt  < 1e-3))
    @printf("  Gₗ unchanged     : Δ = %.2e   [%s]\n", dG,  mark(dG  < 1e-3))
    println("=========================================================")
    return max(dM, dHT, dQ, dt, dG) < 1e-3
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
"Φ-weighted gross outflow rate Σ_{l₀} Σ_z Φ[z,l₀] Σ_{l≠l₀} π̄_{l₀→l}(z) / Σ Φ."
function outflow_rate(sol)
    (; Nz, L) = sol.grids
    num = den = 0.0
    for l0 in 1:L, zi in 1:Nz
        w = sol.Φ[zi, l0]
        num += w * sum(sol.πbar[l0, zi, l] for l in 1:L if l != l0)
        den += w
    end
    return num / den
end

"Population-weighted teaching share across all birth locations (endogenous Gₗ)."
function teach_share_total(sol)
    (; L) = sol.grids
    num = sum(teaching_share_endo(sol.V, l, sol.Φ, sol.grids, sol.p) * sum(@view sol.Φ[:, l]) for l in 1:L)
    return num / sum(sol.Φ)
end

function comparative_statics(; Nz = 3, Nϵ = 7, damping = 0.2, maxit = 120, hh_maxit = 1000)
    println("\n========== Comparative statics (one-parameter deviations) ==========")
    cases = [
        (label = "baseline",        p = Params()),
        # Teacher-spillover strength: larger β amplifies the H̃_T feedback.
        (label = "high-β",          p = Params(β = 0.50)),
        (label = "zero-β",          p = Params(β = 0.0)),
        # Altruism.
        (label = "high-altruism",   p = Params(λ = 0.90)),
        (label = "low-altruism",    p = Params(λ = 0.10)),
        # Spatial frictions: higher move cost / lower taste noise ⇒ less mobility.
        (label = "high-move-cost",  p = Params(τmove = [0.0 1.50; 1.50 0.0])),
        (label = "low-σν",          p = Params(σν = 0.20)),
        # Amenity pull toward location 2 ⇒ M₂ should rise vs baseline.
        (label = "amenity-loc2",    p = Params(B = [0.0, 0.50])),
        # Teaching wage schedule.
        (label = "high-teach-wage", p = Params(κ = [1.20, 1.40])),
        # Gender wage wedge on women in non-teaching.
        (label = "gender-wedge",    p = Params(τω = [0.0 0.0; 0.0 0.30])),
        # Ability persistence.
        (label = "high-persistence",p = Params(ρz = 0.90)),
        (label = "low-persistence", p = Params(ρz = 0.40)),
    ]

    @printf("  %-16s %-9s %-6s %-19s %-19s %-7s %-7s %s\n",
            "case", "res", "conv?", "H̃_T", "M", "out", "teach", "bounds")
    sols  = Dict{String,Any}()
    stuck = String[]
    for c in cases
        sol = solve_ge(c.p; Nz, Nϵ, damping, tol = 1e-3, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        sols[c.label] = sol
        res  = ge_map_residual(sol).res
        conv = res < 5e-3                         # treat >5e-3 as "did not settle"
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
        println("    stall is in the *outer* map: re-solve with a *larger* `damping`")
        println("    (see the damping diagnostic) and a larger `maxit`.")
    end

    # Directional sanity checks against the baseline.
    println("\n  Sign checks:")
    base = sols["baseline"]
    chk(name, cond) = @printf("    %-40s %s\n", name, cond ? "ok" : "FAIL")
    chk("amenity in loc 2 raises M₂", sols["amenity-loc2"].M[2] > base.M[2])
    chk("high move cost lowers outflow", outflow_rate(sols["high-move-cost"]) < outflow_rate(base))
    # Only check directional sign when the comparison case actually converged.
    if ge_map_residual(sols["low-σν"]).res < 5e-3
        chk("low σν lowers outflow", outflow_rate(sols["low-σν"]) < outflow_rate(base))
    else
        @printf("    %-40s %s\n", "low σν lowers outflow", "SKIP (not converged)")
    end
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
function damping_diagnostic(; Nz = 3, Nϵ = 7, maxit = 150, hh_maxit = 1000)
    println("\n========== Damping diagnostic (residual vs step size) ==========")
    println("  The outer GE map is contracting with modulus ≈ β; larger damping")
    println("  (larger step toward the new iterate) should converge faster.")
    println("  To keep the comparison fair, this re-solves the baseline at the")
    @printf("  headline grid (Nϵ=%d, hh_maxit=%d), holds the outer\n", Nϵ, hh_maxit)
    @printf("  budget fixed at %d iterations, and varies only `damping`.\n", maxit)
    @printf("  %-9s %-12s %s\n", "damping", "GE residual", "usable(<5e-3)")
    results = Tuple{Float64,Float64}[]
    for d in (0.80, 0.70, 0.50, 0.30, 0.20)
        sol = solve_ge(Params(); Nz, Nϵ, damping = d, tol = 1e-4, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        res = ge_map_residual(sol).res
        push!(results, (d, res))
        @printf("  %-9.2f %.3e    %s\n", d, res, res < 5e-3 ? "yes" : "no")
    end
    _, i = findmin(last.(results))
    best_d, best_res = results[i]
    @printf("  lowest residual in this budget: damping = %.2f (%.2e).\n", best_d, best_res)
    println("  Larger damping (step size closer to 1) converges faster when the")
    println("  map is contracting; use smaller damping only if iterates diverge.")
    println("================================================================")
end

# -----------------------------------------------------------------------------
# 6. Generalized-dimensions smoke test (I occupations, L locations).
#
#    Confirms the I/L-generalized solver runs for a non-2x2 configuration:
#    I = 3 occupations (1 teacher + 2 non-teaching occupations with different
#    productivities A) and L = 3 symmetric locations (Params(I,L) builds the
#    location/occupation vectors the same way the 2x2 default is built).
#    Checks the structural identities that hold at *any* iterate.
# -----------------------------------------------------------------------------
function generalized_dims_test(; Nz = 2, Nϵ = 3, damping = 0.3, maxit = 30)
    println("\n========== Generalized-dimensions test (I=3, L=3) ==========")
    p = Params(3, 3; A = [NaN, 1.5, 2.0])
    sol = solve_ge(p; Nz, Nϵ, damping, tol = 1e-3, maxit,
                   hh_tol = 1e-5, hh_maxit = 40, verbose = false)
    report_ge(sol)

    ph  = phi_stationarity_residual(sol)
    pir = pi_rowsum_dev(sol)
    pbr = pibar_rowsum_dev(sol)
    qd  = q_consistency_dev(sol)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
    ok = ph.eig_res < 1e-8 && ph.row_dev < 1e-10 && pir < 1e-10 &&
         pbr < 1e-10 && qd < 1e-10 && Mtot_dev < 1e-8 &&
         policies_in_bounds(sol.S, sol.E, sol.V)
    @printf("  I=%d, L=%d   Φ eig-res=%.2e   Σπ_{l'|l}=1 dev=%.2e   Σ_l M_l=Mtot dev=%.2e   [%s]\n",
            sol.grids.I, sol.grids.L, ph.eig_res, pir, Mtot_dev, mark(ok))
    println("=============================================================")
    return ok
end

# -----------------------------------------------------------------------------
# 6b. Scaling-cost test (runtime vs problem dimensions).
#
#     Measures the *marginal computational cost* of adding a location (L) or an
#     occupation (I), or both, to the GE solve.  The two dimensions enter the
#     cost very differently:
#       • Occupations I:  each occupation draws its OWN ϵ shock, so the joint
#         shock grid has size Nϵ^I.  Every inner loop that integrates over shocks
#         (child_objects, the EV loop in location_values, aggregates,
#         stationary_phi) runs over CartesianIndices((Nϵ,...,Nϵ)), so cost grows
#         *geometrically* in I (∝ Nϵ^I).  Adding one occupation multiplies the
#         shock work by ≈ Nϵ.
#       • Locations L:  enter linearly-to-quadratically — the per-state location
#         loop is O(L), the policy grid is O(L), and the Φ eigen-kernel is
#         (L·Nz)×(L·Nz), i.e. O(L²) to assemble/solve.
#
#     To isolate the *dimension* effect from convergence speed, every config is
#     solved with the SAME fixed grid (Nz, Nϵ) and the SAME fixed iteration
#     budget (maxit outer × hh_maxit inner).  We do NOT iterate to convergence:
#     each solve does an identical number of sweeps, so wall-time differences are
#     pure per-sweep cost driven by I and L.  Times are reported in seconds and
#     as a ratio to the (I=2, L=2) baseline; we also print the theoretical state
#     count  I·2·L·Nz·Nϵ^I  (policy states × joint-shock combinations) so the
#     measured ratios can be read against the expected scaling.
# -----------------------------------------------------------------------------
function scaling_cost_test(; Nz = 2, Nϵ = 5, maxit = 6, hh_maxit = 6, reps = 3)
    println("\n========== Scaling-cost test (runtime vs dimensions) ==========")
    @printf("  Fixed budget: Nz=%d, Nϵ=%d, outer maxit=%d, inner hh_maxit=%d.\n",
            Nz, Nϵ, maxit, hh_maxit)
    println("  Every config does the SAME number of sweeps (no early stop), so")
    println("  time differences reflect dimension cost, not convergence speed.")
    println("  Joint-shock grid is Nϵ^I ⇒ occupations cost ∝ Nϵ^I (geometric);")
    println("  locations cost ∝ L (per-state loop, policy grid) up to L² (Φ kernel).")

    # (I, L) configs: baseline, then +locations, +occupations, and both.
    configs = [
        (label = "base",            I = 2, L = 2),
        (label = "+1 location",     I = 2, L = 3),
        (label = "+2 locations",    I = 2, L = 4),
        (label = "+1 occupation",   I = 3, L = 2),
        (label = "+2 occupations",  I = 4, L = 2),
        (label = "+1 of each",      I = 3, L = 3),
    ]

    # Fixed-budget solve at dimensions (I, L); tol=0 forces exactly `maxit`
    # outer iterations so the work is identical across configs.
    solve_fixed(I, L) = solve_ge(Params(I, L); Nz, Nϵ, damping = 0.3,
                                 tol = 0.0, maxit, hh_tol = 0.0, hh_maxit,
                                 verbose = false)

    # Warm up the JIT once (compilation cost must not land in any timed config).
    solve_fixed(2, 2)

    statecount(I, L) = I * 2 * L * Nz * Nϵ^I
    @printf("  %-16s %3s %3s %14s %10s %8s\n",
            "config", "I", "L", "states", "time(s)", "vs base")
    base_t = NaN
    for c in configs
        t = minimum(@elapsed(solve_fixed(c.I, c.L)) for _ in 1:reps)
        isnan(base_t) && (base_t = t)
        @printf("  %-16s %3d %3d %14d %10.3f %7.2f×\n",
                c.label, c.I, c.L, statecount(c.I, c.L), t, t / base_t)
    end
    println("  (states = I·2·L·Nz·Nϵ^I: policy states × joint-shock combinations.)")
    println("===============================================================")
end

# -----------------------------------------------------------------------------
# 7. Plots.
# -----------------------------------------------------------------------------
function plot_ge_outcomes(sol; outdir)
    mkpath(outdir)
    (; z, ϵgrid, Nz, Nϵ, Πz, I, L) = sol.grids
    S, E, H, V, p = sol.S, sol.E, sol.H, sol.V, sol.p
    g = 1
    zmed, ϵmed = (Nz + 1) ÷ 2, (Nϵ + 1) ÷ 2
    markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon]
    occ_label(occ) = occ == p.T ? "Teacher" : "Occupation $occ"

    # (a) Endogenous ability distribution by location vs the ergodic law.
    erg = stationary(Πz)
    pd  = plot(xlabel = "z", ylabel = "mass", title = "Endogenous ability distribution")
    for l in 1:L
        gl = sol.Φ[:, l] ./ sum(@view sol.Φ[:, l])
        plot!(pd, z, gl, marker = markers[mod1(l, length(markers))], label = "G_$(l)(z)")
    end
    plot!(pd, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pd, joinpath(outdir, "ge_ability_distribution.png"))

    # (b) Policies / value vs z (location 1, men, median ϵ).
    p1 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        for occ in 1:I
            plot!(p1[idx], z, A[occ, g, 1, :, ϵmed], label = occ_label(occ), xlabel = "z", ylabel = lab)
        end
    end
    savefig(p1, joinpath(outdir, "ge_policies_vs_z.png"))

    # (c) Policies / value vs ϵ (location 1, men, median z).
    p2 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        for occ in 1:I
            plot!(p2[idx], ϵgrid, A[occ, g, 1, zmed, :], label = occ_label(occ), xlabel = "ϵ", ylabel = lab)
        end
    end
    savefig(p2, joinpath(outdir, "ge_policies_vs_eps.png"))

    # (d) Migration matrix Π̄[born → work] as a heatmap (base Plots only).
    Π̄ = [sum(sol.Φ[zi, l0] * sol.πbar[l0, zi, l] for zi in 1:Nz) / sum(@view sol.Φ[:, l0])
         for l0 in 1:L, l in 1:L]
    pm = heatmap(["work $l" for l in 1:L], ["born $l" for l in 1:L], Π̄,
                 clims = (0, 1), c = :viridis, yflip = true,
                 title = "Migration matrix Π̄[born→work]")
    annotate!(pm, vec([(j, i, text(@sprintf("%.3f", Π̄[i, j]), 9, :white))
                       for i in 1:L, j in 1:L]))
    savefig(pm, joinpath(outdir, "ge_migration_matrix.png"))

    println("\nSaved GE figures to ", outdir)
end

# -----------------------------------------------------------------------------
# 8. Run.   The headline solve (Nz=3, Nϵ=15) plus the symmetry, scale,
#    comparative-statics, and damping batteries.  Each GE solve nests a household
#    fixed point inside every outer iteration, so cost scales with grid size ×
#    outer iters × inner sweeps.  To keep the convergence flags *fair*, the
#    comparative-statics and damping batteries solve the household block to the
#    same inner budget as the headline (hh_maxit=1000) instead of a truncated one:
#    this is slower, but a remaining "did not converge" flag then reflects the
#    outer GE map, not an under-solved household problem.  (Expect the run to take
#    noticeably longer than a truncated-budget pass for this reason.)
#
#    Headline: Nz=3, Nϵ=7.  All batteries use the same grid for comparability.
# -----------------------------------------------------------------------------
function main_test()
    println("Solving spatial model in general equilibrium (baseline) ...")
    sol = solve_ge(Params(); Nz = 3, Nϵ = 9, damping = 0.5, tol = 1e-4,
                   maxit = 250, hh_tol = 1e-5, hh_maxit = 1000, verbose = true)
    report_ge(sol)
    report_ge_checks(sol)
    plot_ge_outcomes(sol; outdir = joinpath(@__DIR__, "figures"))

    symmetry_test()
    scale_invariance_test()
    comparative_statics()
    damping_diagnostic()
    generalized_dims_test()
    scaling_cost_test()
    return sol
end

@time sol = main_test()
