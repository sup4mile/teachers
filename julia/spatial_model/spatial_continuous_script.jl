# =============================================================================
# Diagnostics & tests for the ε = 0, CONTINUOUS-DISTRIBUTION GE spatial solution
# (spatial_continuous.jl).
#
# Standalone counterpart of spatial_ge_eps0_analysis.jl: it `include`s the
# continuous solver and runs the same equilibrium-consistency and economic-
# sense batteries, adapted to the continuous solver's data structures
# (sol.hh.{eT,WT,hT,eO,WO,Λ}, sol.Pi_T / sol.Pi_O, sol.gr in place of S/E/H/V,
# Pi, grids). The structural differences vs the discrete ε = 0 analysis are:
#
#   • occ_θ does not exist here (the atomless integration removes the need for
#     an occupation tremble), so there is no tremble-invariance test. Its role
#     — checking that a numerical stabilization device doesn't move the
#     equilibrium — is played by a GRID-REFINEMENT test: the solution should be
#     stable as the (nϵT, nXO) quadrature grids are refined.
#   • household_residual re-solves ONE extra sweep from the converged Λ (via
#     solve_household(...; maxit=1, Λ0=hh.Λ)) instead of calling update_policy!
#     directly.
#   • an extra monotonicity check on W_T(ε_T) and W_O(X_O*): the choice-map
#     spline inversion (teach_wt/nonteach_wt) silently mislocates the
#     occupation threshold if these are not increasing, so this check guards
#     the assumption spatial_continuous.jl's own warn_nonmonotone also probes.
#
# Run:  julia --project=.. spatial_continuous_script.jl
# =============================================================================

using Printf
using Plots

include(joinpath(@__DIR__, "spatial_continuous.jl"))

gr()  # Plots.jl GR backend (unrelated to the solver's `gr` grids NamedTuple)

# -----------------------------------------------------------------------------
# 0. Generic helpers.
# -----------------------------------------------------------------------------

"Fraction of leading-index slices along which A is non-decreasing in its LAST
 dimension (the shock-grid axis). A has shape (2, L, Nz, nshock)."
function monotone_rate(A; tol = 1e-8)
    front = size(A)[1:end-1]
    pass = 0
    total = 0
    for idx in CartesianIndices(front)
        total += 1
        d = diff(vec(A[Tuple(idx)..., :]))
        all(d .>= -tol) && (pass += 1)
    end
    return pass, total
end

policies_in_bounds(hh) =
    all(isfinite, hh.eT) && all(hh.eT .> 0) &&
    all(isfinite, hh.eO) && all(hh.eO .> 0) &&
    all(isfinite, hh.hT) && all(hh.hT .> 0) &&
    all(isfinite, hh.WT) && all(isfinite, hh.WO)

mark(ok) = ok ? "ok" : "FAIL"

# -----------------------------------------------------------------------------
# 1. Equilibrium-consistency residuals.
# -----------------------------------------------------------------------------

"Residual of the aggregate fixed point: re-evaluate the GE map at `sol` and
 compare to the stored aggregates. Zero only at a true fixed point."
function ge_map_residual(sol)
    HTn, Mn, tn = aggregates(sol.Φ, sol.Pi_T, sol.Pi_O, sol.hh, sol.p, sol.gr)
    rH = maximum(abs, (HTn .- sol.HT) ./ sol.HT)
    rM = maximum(abs, (Mn  .- sol.M ) ./ sol.M )
    rt = maximum(abs,  tn  .- sol.t)
    return (; res = max(rH, rM, rt), rH, rM, rt, HTn, Mn, tn)
end

"Residual of the Φ eigen-problem: rebuild the joint kernel K from π̄ and Πz and
 check that vec(Φ) is its dominant left eigenvector with eigenvalue 1."
function phi_stationarity_residual(sol)
    (; Nz, Πz, L) = sol.gr
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
    eig_res = maximum(abs, vec(φ' * P) .- φ) / sum(φ)
    row_dev = maximum(abs, vec(sum(P, dims = 2)) .- 1)
    return (; eig_res, row_dev, mass = sum(φ))
end

"Worst deviation from 1 of the location-choice rows Σ_{l'} π_{l'|l}, across both
 the teaching (Pi_T) and non-teaching (Pi_O) branches."
function pi_rowsum_dev(sol)
    (; Nz, L, nϵT, nXO) = sol.gr
    dev = 0.0
    for g in 1:2, l in 1:L, zi in 1:Nz
        for k in 1:nϵT
            dev = max(dev, abs(sum(@view sol.Pi_T[g, l, zi, k, :]) - 1))
        end
        for k in 1:nXO
            dev = max(dev, abs(sum(@view sol.Pi_O[g, l, zi, k, :]) - 1))
        end
    end
    return dev
end

"Worst deviation from 1 of the migration rows Σ_{l'} π̄_{l₀→l'}(z)."
pibar_rowsum_dev(sol) =
    maximum(abs, [sum(@view sol.πbar[l0, zi, :]) - 1
                  for l0 in 1:sol.gr.L, zi in 1:sol.gr.Nz])

"One extra household sweep FROM THE CONVERGED Λ: solve_household with maxit=1
 recomputes (WT,WO) at the stored Λ and returns a fresh Λ. The residual is how
 far the recomputed policies are from the stored ones — small only if the
 stored household solution is close to its fixed point."
function household_residual(sol)
    (; hh, gr, p) = sol
    hh1 = solve_household(p, gr; tol = 0.0, maxit = 1, Λ0 = hh.Λ)
    return max(maximum(abs, hh1.WT .- hh.WT), maximum(abs, hh1.WO .- hh.WO))
end

"Consistency of the teacher-quality shifter: Q_l = (2 H̃_T,l / M_l)^σ."
q_consistency_dev(sol) =
    maximum(abs, sol.gr.Q .- (2 .* sol.HT ./ sol.M) .^ sol.p.σ)

"Worst non-monotonicity (largest backward step) in W_T(ε_T) and W_O(X_O*)
 across all (g,l,z) states. The choice-map spline inversion silently mislocates
 the occupation threshold if these are not increasing, so this should be ~0."
function choice_value_monotonicity(sol)
    (; hh, gr) = sol
    (; Nz, L) = gr
    worst = 0.0
    for g in 1:2, l in 1:L, zi in 1:Nz
        wT = @view hh.WT[g, l, zi, :]
        wO = @view hh.WO[g, l, zi, :]
        worst = max(worst, maximum(max.(0.0, .-diff(wT))))
        worst = max(worst, maximum(max.(0.0, .-diff(wO))))
    end
    return worst
end

function report_checks(sol)
    println("\n========== Equilibrium-consistency checks (ε = 0, continuous) ==========")

    println("  Structural identities (exact by construction; should be ~0 always):")
    ph = phi_stationarity_residual(sol)
    mass_dev = abs(ph.mass - sol.p.Mtot / 2)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
    @printf("    Σ Φ = Mtot/2                 : %.2e   [%s]\n", mass_dev, mark(mass_dev < 1e-8))
    @printf("    Σ_l M_l = Mtot               : %.2e   [%s]\n", Mtot_dev, mark(Mtot_dev < 1e-8))
    pir = pi_rowsum_dev(sol)
    @printf("    Σ_{l'} π_{l'|l} = 1  (softmax)  : %.2e   [%s]\n", pir, mark(pir < 1e-6))
    qd = q_consistency_dev(sol)
    @printf("    Q_l = (2H̃_T,l/M_l)^σ         : %.2e   [%s]\n", qd, mark(qd < 1e-10))
    cvm = choice_value_monotonicity(sol)
    @printf("    W_T, W_O increasing in shock  : %.2e   [%s]\n", cvm, mark(cvm < 1e-6))

    # Unlike the discrete solver, π̄ and the Φ eigen-residual come from a
    # QUADRATURE integral (integrate_choice's spline + tail correction), not an
    # exact sum, so they are only accurate to the current (nϵT, nXO) resolution
    # — they shrink under grid refinement (see grid_refinement_test) rather than
    # holding at machine precision. Tolerance here is set to the ballpark the
    # solver's own docstring quotes for its tail correction (~1e-3 at nϵT,nXO~64).
    println("\n  Quadrature-accuracy identities (resolution-dependent, not exact):")
    pbr = pibar_rowsum_dev(sol)
    @printf("    Σ_{l'} π̄_{l₀→l'} = 1  (quadrature) : %.2e   [%s]\n", pbr, mark(pbr < 1e-3))
    @printf("    Φ eigen-residual                   : %.2e   [%s]\n", ph.eig_res, mark(ph.eig_res < 1e-3))
    @printf("    K row-sum deviation                : %.2e   [%s]\n", ph.row_dev, mark(ph.row_dev < 1e-3))

    gm = ge_map_residual(sol)
    Ml_dev = maximum(abs, sol.M .- [2 * sum(@view sol.Φ[:, l]) for l in 1:sol.gr.L])
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
    ok_bounds = policies_in_bounds(sol.hh)
    @printf("    e>0, h>0, W finite (teach & non-teach) : %s\n", mark(ok_bounds))

    println("\n  Monotonicity (household policies, in the shock-grid direction):")
    pe_z, te_z = monotone_rate(sol.hh.eT)
    ph_z, th_z = monotone_rate(sol.hh.hT)
    pw_z, tw_z = monotone_rate(sol.hh.WT)
    @printf("    teach:     e_T in ε_T: %3d/%-3d   h_T in ε_T: %3d/%-3d   W_T in ε_T: %3d/%-3d\n",
            pe_z, te_z, ph_z, th_z, pw_z, tw_z)
    peo_z, teo_z = monotone_rate(sol.hh.eO)
    pwo_z, two_z = monotone_rate(sol.hh.WO)
    @printf("    non-teach: e_O in X_O*: %3d/%-3d                     W_O in X_O*: %3d/%-3d\n",
            peo_z, teo_z, pwo_z, two_z)
    println("==========================================================================")
end

# -----------------------------------------------------------------------------
# 2. Symmetry test.
# -----------------------------------------------------------------------------
function symmetry_test(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.2, maxit = 80)
    println("\n========== Symmetry test (identical locations) ==========")
    p = Params(B = [0.0, 0.0], κ = [1.0, 1.0])
    sol = solve_ge(p; Nz, nϵT, nXO, damping, tol = 1e-3, maxit,
                   hh_tol = 1e-5, hh_maxit = 200, verbose = false)

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
function scale_invariance_test(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.2, maxit = 80)
    println("\n========== Scale-invariance test (Mtot: 2 → 4) ==========")
    kw = (; Nz, nϵT, nXO, damping, tol = 1e-4, maxit, hh_tol = 1e-5, hh_maxit = 400, verbose = false)
    s1 = solve_ge(Params(Mtot = 2.0); kw...)
    s2 = solve_ge(Params(Mtot = 4.0); kw...)

    dM  = maximum(abs, s2.M  .- 2 .* s1.M)
    dHT = maximum(abs, s2.HT .- 2 .* s1.HT)
    dQ  = maximum(abs, s2.gr.Q .- s1.gr.Q)
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
# 3b. Grid-refinement (quadrature resolution) invariance test.
#
# spatial_continuous.jl has no occupation tremble: the atomless integration
# over ε_T and X_O* is what keeps the household map continuous, so there is no
# analogue of the discrete file's occ_θ → 0 continuation test. What DOES need
# checking is that the household/aggregate quadrature grids (nϵT, nXO) are fine
# enough that the reported equilibrium is not a grid artifact — i.e. refining
# the grids should leave the solution (numerically) unchanged. This plays the
# same role the tremble-invariance test plays in the discrete file: confirming
# a purely numerical device (there, occ_θ; here, grid resolution) is not moving
# the economics.
#
# Solves are warm-started from the coarser grid's solution (continuation), both
# for speed and so all resolutions sit on the same equilibrium branch.
# -----------------------------------------------------------------------------
function grid_refinement_test(; Nz = 3, damping = 0.3, tol = 1e-6, maxit = 200,
                               hh_tol = 1e-7, hh_maxit = 800,
                               grids = [(24, 24), (48, 48), (96, 96)],
                               atol = 5e-3)
    println("\n========== Grid-refinement test (nϵT, nXO resolution) ==========")
    sols = Dict{Tuple{Int,Int},Any}()
    prev = nothing
    for (nϵT, nXO) in grids
        sol = solve_ge(Params(); Nz, nϵT, nXO, damping, tol, maxit,
                       hh_tol, hh_maxit, verbose = false, init = prev)
        sols[(nϵT, nXO)] = sol
        prev = sol
    end
    ref = sols[grids[end]]

    @printf("  %-12s %-9s %-9s %-9s %-9s\n", "grid", "res", "ΔH̃_T", "ΔM", "Δt")
    devs = Float64[]
    for (nϵT, nXO) in grids
        s = sols[(nϵT, nXO)]
        dHT = maximum(abs, (s.HT .- ref.HT) ./ max.(abs.(ref.HT), 1e-8))
        dM  = maximum(abs, (s.M  .- ref.M ) ./ max.(abs.(ref.M ), 1e-8))
        dt  = maximum(abs,  s.t  .- ref.t)
        push!(devs, max(dHT, dM, dt))
        res = ge_map_residual(s).res
        @printf("  (%3d,%3d)    %.1e   %.2e   %.2e   %.2e\n", nϵT, nXO, res, dHT, dM, dt)
    end

    # Coarsest-vs-finest gap should already be small, and shrink monotonically
    # (excluding the trivially-zero comparison of the reference to itself).
    shrinking = length(devs) < 3 || issorted(reverse(devs[1:end-1]))
    stable    = devs[end-1] < atol
    ok = stable && shrinking
    @printf("  deviation shrinks with resolution : %s\n", mark(shrinking))
    @printf("  second-finest grid within %.0e of finest : %s\n", atol, mark(stable))
    println("==================================================================")
    return ok
end

# -----------------------------------------------------------------------------
# 4. Comparative statics.
# -----------------------------------------------------------------------------
"Φ-weighted gross outflow rate Σ_{l₀} Σ_z Φ[z,l₀] Σ_{l≠l₀} π̄_{l₀→l}(z) / Σ Φ."
function outflow_rate(sol)
    (; Nz, L) = sol.gr
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
    (; L) = sol.gr
    num = sum(teaching_share_endo(sol.hh, l, sol.Φ, sol.gr, sol.p) * sum(@view sol.Φ[:, l]) for l in 1:L)
    return num / sum(sol.Φ)
end

function comparative_statics(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.2, maxit = 120, hh_maxit = 600)
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
        sol = solve_ge(c.p; Nz, nϵT, nXO, damping, tol = 1e-3, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        sols[c.label] = sol
        res  = ge_map_residual(sol).res
        conv = res < 5e-3
        conv || push!(stuck, c.label)
        ok   = policies_in_bounds(sol.hh)
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
        all(policies_in_bounds(s.hh) for s in values(sols)))
    println("====================================================================")
    return sols
end

# -----------------------------------------------------------------------------
# 5. Damping diagnostic.
# -----------------------------------------------------------------------------
function damping_diagnostic(; Nz = 3, nϵT = 32, nXO = 32, maxit = 150, hh_maxit = 600)
    println("\n========== Damping diagnostic (residual vs step size) ==========")
    @printf("  %-9s %-12s %s\n", "damping", "GE residual", "usable(<5e-3)")
    results = Tuple{Float64,Float64}[]
    for d in (0.80, 0.70, 0.50, 0.30, 0.20)
        sol = solve_ge(Params(); Nz, nϵT, nXO, damping = d, tol = 1e-4, maxit,
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
# 6. Generalized-dimensions smoke test (I=3 occupations, L=3 locations).
#
# The X_O* collapse requires a common τe across non-teaching occupations
# (spatial_continuous.jl asserts this); Params(3,3;...)'s default τe = 0
# everywhere satisfies it.
# -----------------------------------------------------------------------------
function generalized_dims_test(; Nz = 2, nϵT = 48, nXO = 48, damping = 0.3, maxit = 40)
    println("\n========== Generalized-dimensions test (I=3, L=3) ==========")
    p = Params(3, 3; A = [NaN, 1.5, 2.0])
    sol = solve_ge(p; Nz, nϵT, nXO, damping, tol = 1e-3, maxit,
                   hh_tol = 1e-5, hh_maxit = 200, verbose = false)
    report_ge(sol)

    ph  = phi_stationarity_residual(sol)
    pir = pi_rowsum_dev(sol)
    pbr = pibar_rowsum_dev(sol)
    qd  = q_consistency_dev(sol)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
    # ph.eig_res, ph.row_dev, pbr are quadrature-accuracy (not exact) — see
    # report_checks; pir, qd, Mtot_dev are exact structural identities.
    ok = ph.eig_res < 1e-3 && ph.row_dev < 1e-3 && pir < 1e-6 &&
         pbr < 1e-3 && qd < 1e-10 && Mtot_dev < 1e-8 &&
         policies_in_bounds(sol.hh)
    @printf("  I=%d, L=%d   Φ eig-res=%.2e   Σπ_{l'|l}=1 dev=%.2e   Σ_l M_l=Mtot dev=%.2e   [%s]\n",
            sol.gr.I, sol.gr.L, ph.eig_res, pir, Mtot_dev, mark(ok))
    println("=============================================================")
    return ok
end

# -----------------------------------------------------------------------------
# 7. Plots.
# -----------------------------------------------------------------------------
function plot_ge_outcomes(sol; outdir)
    mkpath(outdir)
    (; z, ϵTgrid, XOgrid, Nz, L, Πz) = sol.gr
    hh = sol.hh
    g = 1
    zmed = (Nz + 1) ÷ 2

    erg = stationary(Πz)
    markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon]
    pd = plot(xlabel = "z", ylabel = "mass", title = "Endogenous ability distribution (ε=0, continuous)")
    for l in 1:L
        gl = sol.Φ[:, l] ./ sum(@view sol.Φ[:, l])
        plot!(pd, z, gl, marker = markers[mod1(l, length(markers))], label = "G_$(l)(z)")
    end
    plot!(pd, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pd, joinpath(outdir, "continuous_ability_distribution.png"))

    pT = plot(layout = (1, 3), size = (1200, 350))
    for l in 1:L
        plot!(pT[1], ϵTgrid, hh.eT[g, l, zmed, :], xlabel = "ε_T", ylabel = "e", title = "Teaching e(ε_T)", label = "loc $l")
        plot!(pT[2], ϵTgrid, hh.hT[g, l, zmed, :], xlabel = "ε_T", ylabel = "h", title = "Teaching h(ε_T)", label = "loc $l")
        plot!(pT[3], ϵTgrid, hh.WT[g, l, zmed, :], xlabel = "ε_T", ylabel = "W", title = "Teaching value", label = "loc $l")
    end
    savefig(pT, joinpath(outdir, "continuous_teach_policies.png"))

    pO = plot(layout = (1, 2), size = (900, 350))
    for l in 1:L
        plot!(pO[1], XOgrid[:, g], hh.eO[g, l, zmed, :], xlabel = "X_O*", ylabel = "e", title = "Non-teach e(X_O*)", label = "loc $l")
        plot!(pO[2], XOgrid[:, g], hh.WO[g, l, zmed, :], xlabel = "X_O*", ylabel = "W", title = "Non-teach value", label = "loc $l")
    end
    savefig(pO, joinpath(outdir, "continuous_nonteach_policies.png"))

    Π̄ = [sum(sol.Φ[zi, l0] * sol.πbar[l0, zi, l] for zi in 1:Nz) / sum(@view sol.Φ[:, l0])
         for l0 in 1:L, l in 1:L]
    pm = heatmap(["work $l" for l in 1:L], ["born $l" for l in 1:L], Π̄,
                 clims = (0, 1), c = :viridis, yflip = true,
                 title = "Migration matrix Π̄[born→work] (continuous)")
    annotate!(pm, vec([(j, i, text(@sprintf("%.3f", Π̄[i, j]), 9, :white))
                       for i in 1:L, j in 1:L]))
    savefig(pm, joinpath(outdir, "continuous_migration_matrix.png"))

    println("\nSaved continuous GE figures to ", outdir)
end

# -----------------------------------------------------------------------------
# 7b. Household-choice & distribution INTUITION plots.
#
# plot_ge_outcomes above shows the raw policy/value objects. The figures here are
# aimed at the ECONOMICS a reader wants to take away:
#   (1) the occupation margin — teach iff W_T(ε_T) > W_O(X_O*) — as a smooth
#       teaching PROBABILITY in each shock, and the value crossing that generates it;
#   (2) how teaching SELECTS on ability z (and how location shifts that selection);
#   (3) the shock DISTRIBUTIONS agents draw from, overlaid with who ends up teaching,
#       plus the stationary joint mass over (ability, location).
# All are read-only functions of the converged `sol`; nothing here feeds back into
# the solve, so they are safe to add/remove freely.
# -----------------------------------------------------------------------------

"Gender-pooled teaching probability at a single (location l, ability node zi):
 integrate the teach indicator over BOTH shock margins (ε_T and X_O*), pooling
 the two genders. This is teaching_share_endo's per-cell kernel, exposed by z."
function teach_prob_cell(hh, l, zi, gr, p)
    s = 0.0
    for g in 1:2
        cm = choice_maps(hh.WT[g, l, zi, :], hh.WO[g, l, zi, :], g, gr)
        IT, _ = integrate_choice(ones(gr.nϵT), zeros(gr.nXO), cm, gr)
        s += 0.5 * IT
    end
    return s
end

function plot_household(sol; outdir)
    mkpath(outdir)
    grd = sol.gr
    (; z, ϵTgrid, XOgrid, Nz, L, dT) = grd
    hh = sol.hh
    p  = sol.p
    g  = 1                                     # gender slice for the single-state panels
    zmed = (Nz + 1) ÷ 2
    lmed = (L + 1) ÷ 2
    markers = [:circle, :square, :diamond, :utriangle, :star5, :hexagon]
    zcolors = palette(:viridis, Nz)

    # --- (1) Occupation margin: teaching PROBABILITY in each shock. -----------
    # Left : P(teach | ε_T) = teach_wt — rises in the teaching shock.
    # Right: P(teach | X_O*) = 1 − nonteach_wt — falls in the outside option.
    # One curve per ability node z (at the median location), so the reader sees
    # both how the shock and how ability move the teaching decision.
    p1 = plot(layout = (1, 2), size = (1000, 380))
    for zi in 1:Nz
        cm  = choice_maps(hh.WT[g, lmed, zi, :], hh.WO[g, lmed, zi, :], g, grd)
        ptε = [teach_wt(cm, ε) for ε in ϵTgrid]
        ptx = [1 - nonteach_wt(cm, x) for x in XOgrid[:, g]]
        plot!(p1[1], ϵTgrid, ptε, color = zcolors[zi], lw = 2,
              label = @sprintf("z=%.2f", z[zi]))
        plot!(p1[2], XOgrid[:, g], ptx, color = zcolors[zi], lw = 2, label = "")
    end
    plot!(p1[1], xlabel = "teaching shock ε_T", ylabel = "P(teach)",
          title = "Teaching prob. rises in ε_T", ylims = (-0.02, 1.02),
          legend = :bottomright)
    plot!(p1[2], xlabel = "outside option X_O*", ylabel = "P(teach)",
          title = "Teaching prob. falls in X_O*", ylims = (-0.02, 1.02))
    savefig(p1, joinpath(outdir, "continuous_occupation_margin.png"))

    # --- (2) The value crossing that GENERATES that margin. ------------------
    # W_T(ε_T) and W_O(X_O*) live on a common utility axis but different shock
    # supports; plotting each against its shock QUANTILE (probability rank) puts
    # them on one [0,1] axis, so the marginal teacher (where the curves meet the
    # same choice) is directly visible. Median (location, ability) cell.
    qεgrid = [cdf(dT, ε)     for ε in ϵTgrid]
    qxgrid = [Fxo(grd, g, x) for x in XOgrid[:, g]]
    p2 = plot(size = (760, 420), xlabel = "shock quantile (probability rank)",
              ylabel = "occupation value", legend = :bottomright,
              title = @sprintf("Occupation values vs shock rank  (loc %d, z=%.2f)", lmed, z[zmed]))
    plot!(p2, qεgrid, hh.WT[g, lmed, zmed, :], lw = 2, marker = :circle, ms = 2,
          label = "W_T(ε_T): teach")
    plot!(p2, qxgrid, hh.WO[g, lmed, zmed, :], lw = 2, marker = :square, ms = 2,
          label = "W_O(X_O*): don't teach")
    savefig(p2, joinpath(outdir, "continuous_value_crossing.png"))

    # --- (3) Teaching selection on ability, per location. --------------------
    # Conditional teaching share at each (location, ability) cell: does teaching
    # draw high- or low-z types, and how does location shift the whole schedule?
    p3 = plot(size = (760, 420), xlabel = "ability z", ylabel = "teaching share",
              title = "Teaching share by ability & location", legend = :topright,
              ylims = (-0.02, 1.02))
    for l in 1:L
        ts = [teach_prob_cell(hh, l, zi, grd, p) for zi in 1:Nz]
        plot!(p3, z, ts, lw = 2, marker = markers[mod1(l, length(markers))],
              label = "loc $l")
    end
    savefig(p3, joinpath(outdir, "continuous_teach_share_by_ability.png"))

    # --- (4) Shock distributions, overlaid with who teaches. -----------------
    # Left : ε_T density (teaching shock, lognormal). Right: X_O* density (max
    # over non-teaching occupations). Each density is scaled to unit peak so it
    # shares the [0,1] axis with the dashed P(teach | shock) curve at the median
    # (location, ability): the reader sees which PART of each distribution teaches.
    cm  = choice_maps(hh.WT[g, lmed, zmed, :], hh.WO[g, lmed, zmed, :], g, grd)
    fε  = [pdf(dT, ε)    for ε in ϵTgrid]
    fx  = [fxo(grd, g, x) for x in XOgrid[:, g]]
    ptε = [teach_wt(cm, ε) for ε in ϵTgrid]
    ptx = [1 - nonteach_wt(cm, x) for x in XOgrid[:, g]]
    fεn = fε ./ maximum(fε)
    fxn = fx ./ maximum(fx)
    p4 = plot(layout = (1, 2), size = (1000, 380))
    plot!(p4[1], ϵTgrid, fεn, lw = 2, color = :steelblue, fillrange = 0,
          fillalpha = 0.15, label = "f(ε_T) (scaled)", xlabel = "ε_T",
          ylabel = "scaled density / P(teach)", title = "Teaching-shock distribution",
          ylims = (-0.02, 1.05), legend = :right)
    plot!(p4[1], ϵTgrid, ptε, lw = 2, color = :firebrick, ls = :dash, label = "P(teach | ε_T)")
    plot!(p4[2], XOgrid[:, g], fxn, lw = 2, color = :steelblue, fillrange = 0,
          fillalpha = 0.15, label = "f(X_O*) (scaled)", xlabel = "X_O*",
          ylabel = "scaled density / P(teach)", title = "Outside-option distribution",
          ylims = (-0.02, 1.05), legend = :right)
    plot!(p4[2], XOgrid[:, g], ptx, lw = 2, color = :firebrick, ls = :dash, label = "P(teach | X_O*)")
    savefig(p4, joinpath(outdir, "continuous_shock_distributions.png"))

    # --- (5) Stationary joint mass over (ability, location). -----------------
    # Where the population actually sits in the discrete state space Φ[z, l].
    p5 = heatmap(["loc $l" for l in 1:L],
                 [@sprintf("z=%.2f", z[zi]) for zi in 1:Nz],
                 sol.Φ, c = :viridis, title = "Stationary mass Φ[z, location]")
    savefig(p5, joinpath(outdir, "continuous_joint_distribution.png"))

    println("Saved household-choice intuition figures to ", outdir)
end

# -----------------------------------------------------------------------------
# 8. Run.
# -----------------------------------------------------------------------------
function main_test()
    println("Solving spatial model (ε = 0, continuous) in general equilibrium (baseline) ...")
    sol = solve_ge(Params(); Nz = 5, nϵT = 48, nXO = 48, damping = 0.5, tol = 1e-6,
                   maxit = 250, hh_tol = 1e-6, hh_maxit = 1000, verbose = true)
    report_ge(sol)
    verify_solution(sol)
    report_checks(sol)
    plot_ge_outcomes(sol; outdir = joinpath(@__DIR__, "figures"))
    plot_household(sol; outdir = joinpath(@__DIR__, "figures"))

    symmetry_test()
    scale_invariance_test()
    grid_refinement_test()
    comparative_statics()
    damping_diagnostic()
    generalized_dims_test()
    return sol
end

@time sol = main_test()
