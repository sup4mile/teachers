# =============================================================================
# Diagnostics & tests for the ε = 0 GENERAL-EQUILIBRIUM spatial solution
# (spatial_ge_eps0.jl).
#
# Standalone counterpart of spatial_ge_analysis.jl: it `include`s the ε = 0
# solver (NOT the general-ε spatial_ge.jl) and runs the same equilibrium-
# consistency and economic-sense batteries.  The only structural differences vs
# the general-ε analysis are:
#
#   • household_residual uses child_hc (not child_objects) and the 1-arg-lighter
#     update_policy! (no Dc), matching the ε = 0 solver's signatures;
#   • an extra check verifies the closed-form time share: s equals s_closed(i,p)
#     and is CONSTANT across (z, ϵ, l) — the headline simplification of ε = 0.
#
# Run:  julia --project=.. spatial_ge_eps0_analysis.jl
# =============================================================================

using Printf
using Plots

include(joinpath(@__DIR__, "spatial_ge_eps0.jl"))

gr()

# -----------------------------------------------------------------------------
# 0. Generic helpers.
# -----------------------------------------------------------------------------

"Fraction of (i,g,l) slices along which A is monotone in z and in ϵ.  `either`
 accepts non-decreasing OR non-increasing (used for s, whose sign can flip; with
 ε = 0, s is constant, so it passes trivially)."
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
# -----------------------------------------------------------------------------

"Residual of the aggregate fixed point: re-evaluate the GE map at `sol` and
 compare to the stored aggregates.  Zero only at a true fixed point."
function ge_map_residual(sol)
    HTn, Mn, tn = aggregates(sol.Φ, sol.Pi, sol.H, sol.V, sol.p, sol.grids)
    rH = maximum(abs, (HTn .- sol.HT) ./ sol.HT)
    rM = maximum(abs, (Mn  .- sol.M ) ./ sol.M )
    rt = maximum(abs,  tn  .- sol.t)
    return (; res = max(rH, rM, rt), rH, rM, rt, HTn, Mn, tn)
end

"Residual of the Φ eigen-problem: rebuild the joint kernel K from π̄ and Πz and
 check that vec(Φ) is its dominant left eigenvector with eigenvalue 1."
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
 how close the stored (E,V) are to the household fixed point.  Uses the ε = 0
 signatures: child_hc (no Dc) and the e-only update_policy!."
function household_residual(sol)
    (; S, E, H, V, grids) = sol
    p = sol.p
    E1, V1 = copy(E), copy(V)
    S1, H1 = copy(S), copy(H)
    Hc = child_hc(E, H, V, p, grids)
    update_policy!(S1, E1, H1, V1, Hc, p, grids)
    return max(maximum(abs, E1 .- E), maximum(abs, V1 .- V))
end

"Consistency of the teacher-quality shifter:  Q_l = (2 H̃_T,l / M_l)^σ."
q_consistency_dev(sol) =
    maximum(abs, sol.grids.Q .- (2 .* sol.HT ./ sol.M) .^ sol.p.σ)

"Verify the ε = 0 closed form for s: S[i,...] == s_closed(i,p) at every state,
 i.e. the time share is constant across (g, l, z, ϵ).  Returns the worst |Δ|."
function s_closedform_dev(sol)
    (; S, grids, p) = sol
    (; Nz, Nϵ, I, L) = grids
    dev = 0.0
    for i in 1:I
        target = s_closed(i, p)
        for g in 1:2, l in 1:L, zi in 1:Nz, ei in 1:Nϵ
            dev = max(dev, abs(S[i, g, l, zi, ei] - target))
        end
    end
    return dev
end

function report_ge_checks(sol)
    println("\n========== Equilibrium-consistency checks (ε = 0) ==========")

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
    # The headline ε = 0 identity: s collapses to the closed-form constants.
    sd = s_closedform_dev(sol)
    @printf("    s = s_closed(i,p) (const)    : %.2e   [%s]\n", sd, mark(sd < 1e-12))

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
    @printf("    s monotone in z: %3d/%-3d  in ϵ: %3d/%-3d   (ε=0 ⇒ s constant ⇒ trivially passes)\n",
            ps_z, ts_z, ps_ϵ, ts_ϵ)
    @printf("    e increasing in z: %3d/%-3d in ϵ: %3d/%-3d\n", pe_z, te_z, pe_ϵ, te_ϵ)
    @printf("    V increasing in z: %3d/%-3d in ϵ: %3d/%-3d\n", pv_z, tv_z, pv_ϵ, tv_ϵ)
    println("============================================================")
end

# -----------------------------------------------------------------------------
# 2. Symmetry test.
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
# -----------------------------------------------------------------------------
function scale_invariance_test(; Nz = 3, Nϵ = 5, damping = 0.2, maxit = 80)
    println("\n========== Scale-invariance test (Mtot: 2 → 4) ==========")
    kw = (; Nz, Nϵ, damping, tol = 1e-4, maxit, hh_tol = 1e-5, hh_maxit = 400, verbose = false)
    s1 = solve_ge(Params(Mtot = 2.0); kw...)
    s2 = solve_ge(Params(Mtot = 4.0); kw...)

    dM  = maximum(abs, s2.M  .- 2 .* s1.M)
    dHT = maximum(abs, s2.HT .- 2 .* s1.HT)
    dQ  = maximum(abs, s2.grids.Q .- s1.grids.Q)
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
# 3b. Tremble (occupation-softmax) invariance test.
#
# The child's occupation is smoothed by a softmax with temperature occ_θ — the
# "trembles over occupations".  occ_θ > 0 is a NUMERICAL device: it removes the
# Roy-argmax discontinuity that otherwise drives a limit cycle in the household
# fixed point, so it buys convergence.  Economically it should be innocuous: the
# true equilibrium is the θ → 0 (hard-Roy) limit, and the trembles must not move
# it.  This test solves the GE at a ladder of occ_θ values and checks that the
# aggregates (H̃_T, M, t) and the occupational shares
#   (i)  sit close to the θ → 0 reference at the operating temperature, and
#   (ii) shrink toward that reference monotonically as occ_θ → 0,
# which together say the trembles are not changing the equilibrium.
#
# Numerics.  The tremble effect is O(θ), so at θ = 0.001 the true deviation is
# ~1e-3.  An earlier version solved each θ COLD at tol = 1e-4; the relative
# deviations were then dominated by the GE-map convergence floor — a few×1e-4
# absolute error on H̃_T,1 ≈ 0.04 is already ~1e-2 relative — so the O(θ) signal
# sat *below* the noise and monotonicity failed spuriously (every metric dipped
# together at one θ, the signature of independent solves landing at slightly
# different residuals, not of a real tremble effect).  Two changes fix it:
#   • CONTINUATION: solve θ from largest (smoothest, easiest) down to 0, warm-
#     starting each solve from its larger-θ neighbour.  All solves stay on one
#     branch and reach a tight fixed point in a few iterations — which is also
#     what makes this (the script's heaviest test) cheap.
#   • TIGHT TOLERANCE: tol = hh_tol = 1e-8 with a generous iteration budget, so
#     the convergence floor (~1e-8) is far below the smallest real signal.
# -----------------------------------------------------------------------------

"Economy-wide population share of each occupation i (length I, sums to 1),
 integrated against the endogenous Φ and the SAME softmax weights occ_θ the
 solver uses internally.  occ_θ = 0 ⇒ hard Roy counts."
function occupation_shares(sol)
    (; V, Φ, grids, p) = sol
    (; Nz, Sϵ, Wflat, Edecode, I, L) = grids
    num  = zeros(I)
    wocc = Vector{Float64}(undef, I)
    for l in 1:L, zi in 1:Nz, g in 1:2, sflat in 1:Sϵ
        occupation_weights!(wocc, V, g, l, zi, sflat, Edecode, p.occ_θ)
        w = 0.5 * Φ[zi, l] * Wflat[sflat]      # gender-pool ½, joint-shock weight
        for i in 1:I
            num[i] += w * wocc[i]
        end
    end
    return num ./ sum(num)                      # sum(num) = Σ Φ = Mtot/2
end

"Per-location, gender-pooled teaching shares as a length-L vector."
teach_shares_by_loc(sol) =
    [teaching_share_endo(sol.V, l, sol.Φ, sol.grids, sol.p) for l in 1:sol.grids.L]

function tremble_invariance_test(; Nz = 2, Nϵ = 9,
                                  θs = (0.0, 0.0001, 0.001, 0.005, 0.01),
                                  tol = 1e-8, hh_tol = 1e-8,
                                  maxit = 4000, hh_maxit = 8000,
                                  atol_agg = 1e-2, atol_share = 1e-2,
                                  mono_slack = 1e-6, htfloor = 1e-6, damp=.01)
    println("\n========== Tremble (occupation-softmax) invariance test ==========")
    θs = sort(collect(Float64, θs))

    # Damping schedule: smaller θ is more argmax-like and prone to the household
    # limit cycle, so it gets heavier damping; larger θ is smooth and converges
    # fast at light damping. The warm start (below) compensates for the slower
    # per-step progress at small θ, so the iteration count stays low throughout.
    #gedamp(θ) = θ <= 0.0   ? 0.01 :
                #θ <  0.01  ? 0.10  :
                #             0.20

    # CONTINUATION: solve in DESCENDING θ (smoothest → hardest), warm-starting
    # each solve from its larger-θ neighbour. This keeps every solve on the same
    # equilibrium branch and lets the tight (tol = 1e-8) fixed point be reached in
    # a handful of iterations — the key optimization for this heavy test.
    sols, conv, resd = Dict{Float64,Any}(), Dict{Float64,Bool}(), Dict{Float64,Float64}()
    prev = nothing
    t_sweep = @elapsed for θ in reverse(θs)
        sol      = solve_ge(Params(occ_θ = θ); Nz, Nϵ, damping = damp, # gedamp(θ),
                            tol, maxit, hh_tol, hh_maxit, verbose = false, init = prev)
        sols[θ]  = sol
        resd[θ]  = ge_map_residual(sol).res
        conv[θ]  = resd[θ] < max(1e-3, 100 * tol)
        prev     = sol                     # warm-start the next (smaller) θ
    end

    # Reference = SMALLEST CONVERGED θ — the best available proxy for the θ → 0
    # hard-Roy limit.  (θ = 0 itself may not converge; that's exactly why the
    # trembles exist.)
    converged = sort([θ for θ in θs if conv[θ]])
    if isempty(converged)
        println("  ⚠ no occ_θ converged in the test budget (res ≥ tol for all);")
        println("    cannot assess invariance — raise maxit / hh_maxit or damping.")
        println("==================================================================")
        return false
    end
    θref = first(converged)
    ref  = sols[θref]
    rHT, rM, rt           = ref.HT, ref.M, ref.t
    rocc, rtch            = occupation_shares(ref), teach_shares_by_loc(ref)

    @printf("  reference: occ_θ = %.3f (smallest converged)   I=%d, L=%d   [sweep %.1fs, tol=%.0e]\n",
            θref, ref.grids.I, ref.grids.L, t_sweep, tol)
    @printf("  %-7s %-6s %-9s %-9s %-9s %-9s %-10s %-9s\n",
            "occ_θ", "conv?", "res", "ΔH̃_T", "ΔM", "Δt", "Δocc-shr", "Δteach")
    adev = Dict{Float64,Float64}()      # aggregate block (H̃_T, M, t)
    sdev = Dict{Float64,Float64}()      # occupational-share block
    devs = Dict{Float64,Float64}()      # overall
    for θ in θs
        s   = sols[θ]
        # Relative deviations, but with a small denominator floor so a near-zero
        # aggregate component can't blow the ratio up.
        dHT = maximum(abs, (s.HT .- rHT) ./ max.(abs.(rHT), htfloor))
        dM  = maximum(abs, (s.M  .- rM ) ./ max.(abs.(rM ), htfloor))
        dt  = maximum(abs,  s.t  .- rt)
        doc = maximum(abs, occupation_shares(s) .- rocc)
        dtc = maximum(abs, teach_shares_by_loc(s) .- rtch)
        adev[θ] = max(dHT, dM, dt)
        sdev[θ] = max(doc, dtc)
        devs[θ] = max(adev[θ], sdev[θ])
        @printf("  %-7.3f %-6s %.2e  %.2e  %.2e  %.2e  %.2e   %.2e%s\n",
                θ, conv[θ] ? "yes" : "NO", resd[θ], dHT, dM, dt, doc, dtc,
                conv[θ] ? "" : "   (residual high — read with care)")
    end

    # The property the trembles MUST have is not zero effect at every θ — they do
    # perturb the equilibrium — but that the effect VANISHES in the θ → 0 limit.
    # Verify the deviation from the reference shrinks monotonically as occ_θ → θref
    # (the tremble effect is O(θ) → 0), the defensible "they don't change the
    # equilibrium" claim. mono_slack absorbs residual-level (~tol) jitter between
    # two separately-converged solves so monotonicity is judged on signal, not
    # noise; it stays far below the smallest real O(θ) deviation in the ladder.
    cdev     = [(θ, devs[θ]) for θ in converged]          # ascending θ
    monotone = all(cdev[k][2] <= cdev[k+1][2] + mono_slack for k in 1:length(cdev)-1)

    # Practical ceiling: the largest converged tremble that still keeps a block
    # within tolerance of the θ → 0 limit (θref is 0-deviation by construction).
    ceiling(dev, atol) = maximum(θ for θ in converged if dev[θ] < atol)
    agg_ceil   = ceiling(adev, atol_agg)
    share_ceil = ceiling(sdev, atol_share)
    fmtceil(c) = c <= θref ? @sprintf("only at/below reference occ_θ = %.3f", θref) :
                             @sprintf("occ_θ ≤ %.3f", c)

    println()
    @printf("  converges to θ → 0 limit (monotone in occ_θ)   [%s]\n", mark(monotone))
    @printf("  aggregates within %.0e of limit : %s\n", atol_agg,  fmtceil(agg_ceil))
    @printf("  occ shares within %.0e of limit : %s\n", atol_share, fmtceil(share_ceil))
    θref == 0.0 ||
        @printf("  ⚠ reference is occ_θ = %.3f, not 0 (hard-Roy θ=0 did not converge) — deviations are vs that proxy, not the true θ → 0 limit.\n", θref)
    ok = monotone && θref == 0.0
    @printf("  trembles vanish in the θ → 0 limit             : %s\n", mark(ok))
    println("==================================================================")
    return ok
end

# -----------------------------------------------------------------------------
# 4. Comparative statics.
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
        (label = "high-β",          p = Params(β = 0.50)),
        (label = "zero-β",          p = Params(β = 0.0)),
        (label = "high-altruism",   p = Params(λ = 0.90)),
        (label = "low-altruism",    p = Params(λ = 0.10)),
        (label = "zero-altruism",   p = Params(λ = 0.0)),   # ε=0 AND λ=0 ⇒ no intergenerational coupling
        (label = "high-move-cost",  p = Params(τmove = [0.0 1.50; 1.50 0.0])),
        (label = "low-σν",          p = Params(σν = 0.20)),
        (label = "amenity-loc2",    p = Params(B = [0.0, 0.50])),
        (label = "high-teach-wage", p = Params(κ = [1.20, 1.40])),
        (label = "gender-wedge",    p = Params(τω = [0.0 0.0; 0.0 0.30])),
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
        conv = res < 5e-3
        conv || push!(stuck, c.label)
        ok   = policies_in_bounds(sol.S, sol.E, sol.V)
        @printf("  %-16s %.1e  %-6s (%.3f,%.3f)       (%.3f,%.3f)      %.3f  %.3f  %s\n",
                c.label, res, conv ? "yes" : "NO", sol.HT[1], sol.HT[2],
                sol.M[1], sol.M[2], outflow_rate(sol), teach_share_total(sol), mark(ok))
    end
    if !isempty(stuck)
        println("\n  ⚠ did not converge in the test budget (read their aggregates with care): ",
                join(stuck, ", "))
    end

    println("\n  Sign checks:")
    base = sols["baseline"]
    chk(name, cond) = @printf("    %-40s %s\n", name, cond ? "ok" : "FAIL")
    chk("amenity in loc 2 raises M₂", sols["amenity-loc2"].M[2] > base.M[2])
    chk("high move cost lowers outflow", outflow_rate(sols["high-move-cost"]) < outflow_rate(base))
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
# -----------------------------------------------------------------------------
function damping_diagnostic(; Nz = 3, Nϵ = 7, maxit = 150, hh_maxit = 1000)
    println("\n========== Damping diagnostic (residual vs step size) ==========")
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
    println("================================================================")
end

# -----------------------------------------------------------------------------
# 6. Generalized-dimensions smoke test (I occupations, L locations).
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
    sd  = s_closedform_dev(sol)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
    ok = ph.eig_res < 1e-8 && ph.row_dev < 1e-10 && pir < 1e-10 &&
         pbr < 1e-10 && qd < 1e-10 && sd < 1e-12 && Mtot_dev < 1e-8 &&
         policies_in_bounds(sol.S, sol.E, sol.V)
    @printf("  I=%d, L=%d   Φ eig-res=%.2e   Σπ_{l'|l}=1 dev=%.2e   s-closed dev=%.2e   Σ_l M_l=Mtot dev=%.2e   [%s]\n",
            sol.grids.I, sol.grids.L, ph.eig_res, pir, sd, Mtot_dev, mark(ok))
    println("=============================================================")
    return ok
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

    erg = stationary(Πz)
    pd  = plot(xlabel = "z", ylabel = "mass", title = "Endogenous ability distribution (ε=0)")
    for l in 1:L
        gl = sol.Φ[:, l] ./ sum(@view sol.Φ[:, l])
        plot!(pd, z, gl, marker = markers[mod1(l, length(markers))], label = "G_$(l)(z)")
    end
    plot!(pd, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pd, joinpath(outdir, "ge_eps0_ability_distribution.png"))

    p1 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        for occ in 1:I
            plot!(p1[idx], z, A[occ, g, 1, :, ϵmed], label = occ_label(occ), xlabel = "z", ylabel = lab)
        end
    end
    savefig(p1, joinpath(outdir, "ge_eps0_policies_vs_z.png"))

    p2 = plot(layout = (2, 2), size = (900, 650))
    for (idx, (A, lab)) in enumerate(((S, "s"), (E, "e"), (H, "h"), (V, "V")))
        for occ in 1:I
            plot!(p2[idx], ϵgrid, A[occ, g, 1, zmed, :], label = occ_label(occ), xlabel = "ϵ", ylabel = lab)
        end
    end
    savefig(p2, joinpath(outdir, "ge_eps0_policies_vs_eps.png"))

    Π̄ = [sum(sol.Φ[zi, l0] * sol.πbar[l0, zi, l] for zi in 1:Nz) / sum(@view sol.Φ[:, l0])
         for l0 in 1:L, l in 1:L]
    pm = heatmap(["work $l" for l in 1:L], ["born $l" for l in 1:L], Π̄,
                 clims = (0, 1), c = :viridis, yflip = true,
                 title = "Migration matrix Π̄[born→work] (ε=0)")
    annotate!(pm, vec([(j, i, text(@sprintf("%.3f", Π̄[i, j]), 9, :white))
                       for i in 1:L, j in 1:L]))
    savefig(pm, joinpath(outdir, "ge_eps0_migration_matrix.png"))

    println("\nSaved ε=0 GE figures to ", outdir)
end

# -----------------------------------------------------------------------------
# 8. Run.
# -----------------------------------------------------------------------------
function main_test()
    println("Solving spatial model (ε = 0) in general equilibrium (baseline) ...")
    sol = solve_ge(Params(); Nz = 3, Nϵ = 9, damping = 0.5, tol = 1e-4,
                   maxit = 250, hh_tol = 1e-5, hh_maxit = 1000, verbose = true)
    report_ge(sol)
    report_ge_checks(sol)
    plot_ge_outcomes(sol; outdir = joinpath(@__DIR__, "figures"))

    symmetry_test()
    scale_invariance_test()
    tremble_invariance_test()
    comparative_statics()
    damping_diagnostic()
    generalized_dims_test()
    return sol
end

@time sol = main_test()
