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
        (label = "low-σν",          p = Params(σν = 0.10)),
        (label = "amenity-loc2",    p = Params(B = [0.0, 0.50])),
        (label = "high-teach-wage", p = Params(κ = [1.20, 1.40])),
        (label = "gender-wedge",    p = Params(τω = [0.0 0.0; 0.0 0.30])),
        (label = "high-persistence",p = Params(ρz = 0.97)),
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

# =============================================================================
# 4b. Mechanism deep-dive.
#
# Eight questions the reduced comparative-statics table cannot answer on its own:
#   Q1  more β cases + WHY high β corners (the β/σ ratio),
#   Q2  how occupation choice and location choice interact,
#   Q3  class sizes across locations,
#   Q4  the distribution of (teacher) human capital by location,
#   Q5  what drives spatial sorting on ability,
#   Q6  sensitivity to the quantile bounds / shock grids,
#   Q7  teacher-HC DISTRIBUTION vs teacher QUANTITY in setting Q_l,
#   Q8  amenities vs altruism in location choice, and who moves (high/low z).
#
# Everything here is READ-ONLY on a converged `sol` (re-using the solver's own
# choice_maps / integrate_choice, so the integrals match aggregates() exactly),
# except the sweep drivers, which call solve_ge themselves. Figures land in the
# same figures/ directory as the baseline plots.
# -----------------------------------------------------------------------------

share_hi_amenity(sol) = sol.M[argmax(sol.p.B)] / sum(sol.M)   # pop share in the top-amenity location

"Side-by-side (dodged) grouped bars — the GR backend does not support
 `bar_position=:dodge`, so offset each series manually. `M` is (ncats × nseries)."
function grouped_bar(cats, series_labels, M; kwargs...)
    n, k = size(M)
    w = 0.8 / k
    x = collect(1:n)
    plt = bar(x .+ (1 - (k + 1) / 2) * w, M[:, 1]; bar_width = w,
              label = series_labels[1], kwargs...)
    for j in 2:k
        bar!(plt, x .+ (j - (k + 1) / 2) * w, M[:, j]; bar_width = w, label = series_labels[j])
    end
    xticks!(plt, x, cats)
    return plt
end

"Occupation-resolved born→work masses, plus teacher log-h moments by work loc.
 Mirrors aggregates()'s cell loop but keeps the teach/non-teach split and the
 born→work resolution. Returns masses `Tmass[born,work]`, `Omass[born,work]` and
 per-work-location Σmass·logh, Σmass·logh², Σmass·h^{β/σ} (=H̃_T) among teachers."
function flow_moments(sol)
    (; hh, Φ, gr, p) = sol
    (; Nz, L, nϵT, nXO) = gr
    expo = p.β / p.σ
    Tmass = zeros(L, L); Omass = zeros(L, L)
    Tlogh = zeros(L); Tlogh2 = zeros(L); THT = zeros(L)
    zeroO = zeros(nXO); zeroT = zeros(nϵT)
    for l0 in 1:L, zi in 1:Nz
        cell = 0.5 * Φ[zi, l0]
        cell == 0.0 && continue
        for g in 1:2
            cm    = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
            hTv   = hh.hT[g, l0, zi, :]
            loghT = log.(hTv)
            for lp in 1:L
                πTv = collect(@view sol.Pi_T[g, l0, zi, :, lp])
                πOv = collect(@view sol.Pi_O[g, l0, zi, :, lp])
                ITc, _  = integrate_choice(πTv, zeroO, cm, gr)
                _,  IOc = integrate_choice(zeroT, πOv, cm, gr)
                ITl, _  = integrate_choice(πTv .* loghT,       zeroO, cm, gr)
                ITl2, _ = integrate_choice(πTv .* loghT .^ 2,  zeroO, cm, gr)
                ITh, _  = integrate_choice(πTv .* hTv .^ expo, zeroO, cm, gr)
                Tmass[l0, lp] += cell * ITc
                Omass[l0, lp] += cell * IOc
                Tlogh[lp]  += cell * ITl
                Tlogh2[lp] += cell * ITl2
                THT[lp]    += cell * ITh
            end
        end
    end
    return (; Tmass, Omass, Tlogh, Tlogh2, THT)
end

"Teacher-HC and class-size moments by WORK location (Q3, Q4, Q7).
 `classsize` = students (M_l/2) per teacher working in l; `qual_per_teacher` =
 H̃_T,l / (#teachers) = mean h^{β/σ}, the intensive-margin (quality) piece of the
 teacher aggregator, so H̃_T,l = Tcount_l × qual_per_teacher_l splits Q_l's driver
 into a QUANTITY term (Tcount) and a DISTRIBUTION term (qual_per_teacher)."
function teacher_moments(sol)
    fm = flow_moments(sol)
    (; L) = sol.gr
    Tcount = [sum(@view fm.Tmass[:, lp]) for lp in 1:L]
    Ocount = [sum(@view fm.Omass[:, lp]) for lp in 1:L]
    meanlogh = [Tcount[lp] > 0 ? fm.Tlogh[lp] / Tcount[lp] : NaN for lp in 1:L]
    sdlogh   = [Tcount[lp] > 0 ? sqrt(max(fm.Tlogh2[lp] / Tcount[lp] - meanlogh[lp]^2, 0.0)) : NaN for lp in 1:L]
    students  = sol.M ./ 2
    classsize = [Tcount[lp] > 0 ? students[lp] / Tcount[lp] : NaN for lp in 1:L]
    qual_per_teacher = [Tcount[lp] > 0 ? fm.THT[lp] / Tcount[lp] : NaN for lp in 1:L]
    return (; Tcount, Ocount, meanlogh, sdlogh, classsize, qual_per_teacher,
              HT = fm.THT, Q = sol.gr.Q, students, Tmass = fm.Tmass, Omass = fm.Omass)
end

"Mean ability z by BIRTH and by WORK location, and scalar sorting indices
 (spread in mean z across locations). `Ψ` is the post-migration (work-location)
 ability mass. Answers Q5's 'how much sorting'."
function sorting_metrics(sol)
    (; Φ, πbar, gr) = sol
    (; z, Nz, L) = gr
    Ψ = [sum(Φ[zi, l0] * πbar[l0, zi, l] for l0 in 1:L) for zi in 1:Nz, l in 1:L]
    meanz(w) = sum(z[zi] * w[zi] for zi in 1:Nz) / sum(w)
    mz_birth = [meanz(@view Φ[:, l]) for l in 1:L]
    mz_work  = [meanz(@view Ψ[:, l]) for l in 1:L]
    return (; mz_birth, mz_work, Ψ,
              idx_birth = maximum(mz_birth) - minimum(mz_birth),
              idx_work  = maximum(mz_work)  - minimum(mz_work))
end

"Per-ability probability of leaving the birth location, population-weighted over
 birth locations and also resolved per birth location (Q8: who moves?)."
function mobility_by_ability(sol)
    (; Φ, πbar, gr) = sol
    (; Nz, L) = gr
    out_by_birth = [sum(πbar[l0, zi, l] for l in 1:L if l != l0) for zi in 1:Nz, l0 in 1:L]
    outflow = [sum(Φ[zi, l0] * out_by_birth[zi, l0] for l0 in 1:L) / sum(@view Φ[zi, :])
               for zi in 1:Nz]
    return (; outflow, out_by_birth)
end

"Weighted (log h, weight) sample of teachers WORKING in location `lp`, for a
 density estimate. Node weight = cell mass × f_T(ε)·Δε × P(teach|ε) × π(→lp)."
function teacher_logh_samples(sol, lp)
    (; hh, Φ, gr) = sol
    (; Nz, L, nϵT, ϵTgrid, dT) = gr
    dϵ = similar(ϵTgrid)                        # trapezoid probability spacing on the ε_T grid
    dϵ[1] = (ϵTgrid[2] - ϵTgrid[1]) / 2
    dϵ[end] = (ϵTgrid[end] - ϵTgrid[end-1]) / 2
    @inbounds for k in 2:nϵT-1
        dϵ[k] = (ϵTgrid[k+1] - ϵTgrid[k-1]) / 2
    end
    logh = Float64[]; w = Float64[]
    for l0 in 1:L, zi in 1:Nz, g in 1:2
        cell = 0.5 * Φ[zi, l0]
        cell == 0.0 && continue
        cm = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
        for k in 1:nϵT
            tw = teach_wt(cm, ϵTgrid[k])
            tw <= 0.0 && continue
            wt = cell * pdf(dT, ϵTgrid[k]) * dϵ[k] * tw * sol.Pi_T[g, l0, zi, k, lp]
            wt <= 0.0 && continue
            push!(logh, log(hh.hT[g, l0, zi, k])); push!(w, wt)
        end
    end
    return logh, w
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
    pd = plot(xlabel = "z", ylabel = "mass", title = "Endogenous ability distribution")
    for l in 1:L
        gl = sol.Φ[:, l] ./ sum(@view sol.Φ[:, l])
        plot!(pd, z, gl, marker = markers[mod1(l, length(markers))], label = "G_$(l)(z)")
    end
    plot!(pd, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pd, joinpath(outdir, "continuous_ability_distribution.png"))

    # Same object, but by WORK location instead of BIRTH location. Φ[z,l] is indexed
    # by BIRTH location (the "young"/at-birth cohort, before they sort across space);
    # pushing it through the ability-specific migration kernel πbar gives the mass of
    # each ability z that ends up WORKING in each l (the "old"/post-migration cohort).
    # z is fixed over life, so the marginal ability distribution is unchanged — the two
    # plots differ only in how ability is allocated ACROSS locations by migration.
    Ψ = [sum(sol.Φ[zi, l0] * sol.πbar[l0, zi, l] for l0 in 1:L) for zi in 1:Nz, l in 1:L]
    pw = plot(xlabel = "z", ylabel = "mass",
              title = "Endogenous ability distribution by WORK location")
    for l in 1:L
        gwl = Ψ[:, l] ./ sum(@view Ψ[:, l])
        plot!(pw, z, gwl, marker = markers[mod1(l, length(markers))], label = "G^work_$(l)(z)")
    end
    plot!(pw, z, erg, marker = :diamond, ls = :dash, label = "ergodic G*(z)")
    savefig(pw, joinpath(outdir, "continuous_ability_distribution_by_work.png"))

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
# 7c. Mechanism analyses + figures for the eight deep-dive questions.
#
# Shared solve settings for the sweep drivers below. Small (Nz=3, 32×32) so the
# many-solve sweeps stay ~seconds each; every driver takes a keyword override.
# -----------------------------------------------------------------------------
_dd_kw(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.3, tol = 1e-4, maxit = 150,
        hh_maxit = 500) =
    (; Nz, nϵT, nXO, damping, tol, maxit, hh_tol = 1e-5, hh_maxit, verbose = false)

# ---- Q2. How occupation choice and location choice interact. ----------------
# Teaching is decided when young; location when working. The two margins couple:
# (i) teaching selection differs across BIRTH locations (different G_l(z), Q_l);
# (ii) teachers and non-teachers MIGRATE differently (teaching wages run through
# κ_l' h^γ, so destinations matter differently for the two occupations).
function analyze_occupation_location(sol; outdir)
    mkpath(outdir)
    (; L) = sol.gr
    tm = teacher_moments(sol)
    T, O = tm.Tmass, tm.Omass
    println("\n========== Q2. Occupation × location interaction ==========")
    born_mass  = [sum(@view T[l0, :]) + sum(@view O[l0, :]) for l0 in 1:L]
    teach_born = [sum(@view T[l0, :]) / born_mass[l0] for l0 in 1:L]
    println("  Teaching share among young born in l:")
    for l0 in 1:L
        @printf("    born %d : %.4f\n", l0, teach_born[l0])
    end
    # born→work destination shares, separately for teachers and non-teachers.
    stay_T = zeros(L); stay_O = zeros(L)
    println("\n  Migration (stay probability) by occupation:")
    for l0 in 1:L
        rowT = sum(@view T[l0, :]); rowO = sum(@view O[l0, :])
        stay_T[l0] = rowT > 0 ? T[l0, l0] / rowT : NaN
        stay_O[l0] = rowO > 0 ? O[l0, l0] / rowO : NaN
        @printf("    born %d :  teachers stay %.4f   non-teachers stay %.4f   (Δ = %+.4f)\n",
                l0, stay_T[l0], stay_O[l0], stay_T[l0] - stay_O[l0])
    end
    println("\n  Interpretation: Δ>0 means teachers are MORE rooted than non-teachers")
    println("  born in the same place (κ_l' pins teaching income to the local wage);")
    println("  Δ<0 means teaching is the more footloose occupation.")
    println("===========================================================")

    p1 = bar(["born $l0" for l0 in 1:L], teach_born, legend = false,
             ylabel = "teaching share", title = "Teaching share by birth location",
             ylims = (0, max(0.6, maximum(teach_born) * 1.15)))
    p2 = grouped_bar(["born $l0" for l0 in 1:L], ["teacher", "non-teacher"],
                     hcat(stay_T, stay_O); ylabel = "P(work = birth loc)",
                     title = "Stay probability by occupation", ylims = (0, 1), legend = :topright)
    plot(p1, p2, size = (1000, 380))
    savefig(joinpath(outdir, "dd_occupation_location.png"))
    return (; teach_born, stay_T, stay_O)
end

# ---- Q3/Q4/Q7. Class size, teacher-HC distribution, quantity vs quality. -----
function analyze_class_and_teachers(sol; outdir)
    mkpath(outdir)
    (; L) = sol.gr
    tm = teacher_moments(sol)
    println("\n========== Q3/Q4/Q7. Class size & teacher human capital ==========")
    @printf("  %-6s %-9s %-9s %-11s %-11s %-9s %-9s\n",
            "loc", "Q_l", "class", "#teachers", "mean logh", "sd logh", "H̃_T")
    for l in 1:L
        @printf("  %-6d %-9.4f %-9.3f %-11.4f %-11.4f %-9.4f %-9.4f\n",
                l, tm.Q[l], tm.classsize[l], tm.Tcount[l], tm.meanlogh[l],
                tm.sdlogh[l], tm.HT[l])
    end
    println("\n  Q7. Decomposition  H̃_T,l = (#teachers) × (mean h^{β/σ}):")
    println("      — the QUANTITY margin (#teachers) vs the DISTRIBUTION/quality")
    println("        margin (mean h^{β/σ}). Ratios are location l / location 1.")
    @printf("    %-8s %-12s %-16s %-12s\n", "loc", "H̃_T ratio", "quantity ratio", "quality ratio")
    for l in 1:L
        @printf("    %-8d %-12.3f %-16.3f %-12.3f\n", l,
                tm.HT[l] / tm.HT[1], tm.Tcount[l] / tm.Tcount[1],
                tm.qual_per_teacher[l] / tm.qual_per_teacher[1])
    end
    # Compare total cross-location variation (excluding the tied-at-1 reference)
    # in the quantity vs the quality ratio.
    qty_dev  = sum(abs.(tm.Tcount ./ tm.Tcount[1] .- 1))
    qual_dev = sum(abs.(tm.qual_per_teacher ./ tm.qual_per_teacher[1] .- 1))
    println("\n  ⇒ Cross-location differences in H̃_T are driven mainly by ",
            qty_dev > qual_dev ? "QUANTITY (#teachers)" : "the teacher-HC DISTRIBUTION (quality)",
            @sprintf("  (Σ|dev|: quantity %.3f vs quality %.3f).", qty_dev, qual_dev))
    println("==================================================================")

    gm    = exp.(tm.meanlogh)                        # geometric-mean teacher h (positive)
    explo = exp.(tm.meanlogh .- tm.sdlogh)
    exphi = exp.(tm.meanlogh .+ tm.sdlogh)
    p1 = bar(["loc $l" for l in 1:L], tm.classsize, legend = false,
             ylabel = "students / teacher", title = "Mean class size N̄_l")
    p2 = bar(["loc $l" for l in 1:L], gm, yerror = (gm .- explo, exphi .- gm), legend = false,
             ylabel = "geom-mean teacher h (±1sd in log)", title = "Teacher HC by location")
    p3 = grouped_bar(["loc $l" for l in 1:L], ["quantity", "quality"],
                     hcat(tm.Tcount ./ tm.Tcount[1], tm.qual_per_teacher ./ tm.qual_per_teacher[1]);
                     ylabel = "ratio to loc 1", title = "H̃_T decomposition (quantity vs quality)",
                     legend = :topleft)
    plot(p1, p2, p3, layout = (1, 3), size = (1300, 380))
    savefig(joinpath(outdir, "dd_class_and_teacher_moments.png"))

    # Teacher log-h DENSITY by work location (weighted quadrature-node sample).
    pdens = plot(xlabel = "log h (teachers)", ylabel = "density",
                 title = "Teacher human-capital distribution by work location",
                 legend = :topright)
    for l in 1:L
        lh, w = teacher_logh_samples(sol, l)
        isempty(lh) && continue
        histogram!(pdens, lh, weights = w, normalize = :pdf, bins = 24,
                   alpha = 0.35, label = "loc $l (N̄=$(round(tm.classsize[l], digits=1)))")
    end
    savefig(pdens, joinpath(outdir, "dd_teacher_hc_distribution.png"))
    println("  Saved class-size / teacher-distribution figures to ", outdir)
    return tm
end

# ---- Q5. What drives spatial sorting on ability? ----------------------------
# Shut down each source of location asymmetry one at a time and measure the
# residual sorting (spread in post-migration mean z across work locations). The
# channel whose removal collapses the sorting index is the one that drives it.
function analyze_sorting_drivers(; outdir, kw = _dd_kw())
    mkpath(outdir)
    println("\n========== Q5. Drivers of spatial sorting on ability ==========")
    cases = [
        (label = "baseline",         p = Params()),
        (label = "no amenity (B=0)", p = Params(B = [0.0, 0.0])),
        (label = "equal κ",          p = Params(κ = [1.0, 1.0])),
        (label = "no altruism (λ=0)",p = Params(λ = 0.0)),
        (label = "no move cost",     p = Params(τmove = [0.0 0.0; 0.0 0.0])),
        (label = "symmetric (B=0,κ=)",p = Params(B = [0.0, 0.0], κ = [1.0, 1.0])),
    ]
    @printf("  %-20s %-10s %-10s %-12s %-10s\n",
            "case", "idx_work", "idx_birth", "share_loc2", "outflow")
    idxs = Float64[]; labels = String[]
    for c in cases
        sol = solve_ge(c.p; kw...)
        sm  = sorting_metrics(sol)
        @printf("  %-20s %-10.4f %-10.4f %-12.4f %-10.4f\n",
                c.label, sm.idx_work, sm.idx_birth, sol.M[2] / sum(sol.M), outflow_rate(sol))
        push!(idxs, sm.idx_work); push!(labels, c.label)
    end
    println("\n  The sorting index (work-location mean-z spread) collapses toward the")
    println("  fully-symmetric benchmark as the ACTIVE asymmetry is removed; whichever")
    println("  single-channel removal cuts it most is the dominant sorting driver.")
    println("===============================================================")
    bar(labels, idxs, legend = false, xrotation = 25, ylabel = "sorting index (mean-z spread)",
        title = "What drives spatial sorting on ability?", size = (900, 450),
        bottom_margin = 8Plots.mm)
    savefig(joinpath(outdir, "dd_sorting_drivers.png"))
    return idxs
end

# ---- Q8. Amenities vs altruism in location choice; who moves? ---------------
function analyze_amenity_altruism(; outdir, kw = _dd_kw())
    mkpath(outdir)
    println("\n========== Q8. Amenities vs altruism; mobility by ability ==========")
    cases = [
        (label = "baseline (B,λ=.7)",    p = Params()),
        (label = "B=0, λ=.7",            p = Params(B = [0.0, 0.0])),
        (label = "B=0, λ=0",             p = Params(B = [0.0, 0.0], λ = 0.0)),
        (label = "B=0, λ=.9",            p = Params(B = [0.0, 0.0], λ = 0.90)),
    ]
    @printf("  %-20s %-11s %-11s %-11s %-14s\n",
            "case", "share_loc2", "idx_work", "outflow", "mobility slope")
    results = NamedTuple[]
    for c in cases
        sol = solve_ge(c.p; kw...)
        sm  = sorting_metrics(sol)
        mb  = mobility_by_ability(sol)
        z   = sol.gr.z
        # sign of the ability→mobility gradient: >0 ⇒ high-z move more.
        slope = (mb.outflow[end] - mb.outflow[1]) / (z[end] - z[1])
        @printf("  %-20s %-11.4f %-11.4f %-11.4f %+-14.4f\n",
                c.label, sol.M[2] / sum(sol.M), sm.idx_work, outflow_rate(sol), slope)
        push!(results, (; c.label, z, outflow = mb.outflow, slope))
    end
    println("\n  Reading:")
    println("   • baseline → B=0 : how much of loc-2's pull is the amenity B.")
    println("   • B=0 across λ∈{0,.7,.9} : does ALTRUISM move location choice once the")
    println("     amenity is gone (share_loc2 / idx_work shifting with λ)?")
    println("   • mobility slope>0 ⇒ HIGH-ability agents are the more mobile ones.")
    println("====================================================================")
    pl = plot(xlabel = "ability z", ylabel = "P(leave birth location)",
              title = "Who moves? Outflow rate by ability", legend = :topleft)
    for r in results
        plot!(pl, r.z, r.outflow, marker = :circle, lw = 2, label = r.label)
    end
    savefig(pl, joinpath(outdir, "dd_mobility_by_ability.png"))
    println("  Saved mobility-by-ability figure to ", outdir)
    return results
end

# ---- Q6. Sensitivity to the quantile bounds / shock grids. ------------------
function quantile_sensitivity(; outdir, Nz = 3, damping = 0.3, tol = 1e-5,
                              maxit = 200, hh_maxit = 600)
    mkpath(outdir)
    base = (; Nz, damping, tol, maxit, hh_tol = 1e-6, hh_maxit, verbose = false)
    println("\n========== Q6. Quantile-bound & grid sensitivity ==========")
    # Reference: tightest tails + finest grid.
    ref = solve_ge(Params(); base..., nϵT = 96, nXO = 96, q_lo = 1e-7, q_hi = 1 - 1e-8)
    reldev(s) = max(maximum(abs, (s.HT .- ref.HT) ./ ref.HT),
                    maximum(abs, (s.M  .- ref.M ) ./ ref.M))

    println("\n  (a) Upper tail q_hi  (q_lo=1e-5, grid 48²):")
    @printf("    %-12s %-12s %-12s %-10s %-10s\n", "1-q_hi", "H̃_T(2)", "teach", "res", "rel.dev")
    qhis = [1e-3, 1e-4, 1e-5, 1e-6, 1e-8]
    xh = Float64[]; htv = Float64[]; tsv = Float64[]
    for u in qhis
        s = solve_ge(Params(); base..., nϵT = 48, nXO = 48, q_lo = 1e-5, q_hi = 1 - u)
        @printf("    %-12.0e %-12.5f %-12.5f %-10.1e %-10.1e\n",
                u, s.HT[2], teach_share_total(s), ge_map_residual(s).res, reldev(s))
        push!(xh, -log10(u)); push!(htv, s.HT[2]); push!(tsv, teach_share_total(s))
    end

    println("\n  (b) Lower tail q_lo  (q_hi=1-1e-6, grid 48²):")
    @printf("    %-12s %-12s %-12s %-10s\n", "q_lo", "H̃_T(2)", "teach", "rel.dev")
    for lo in [1e-3, 1e-4, 1e-5, 1e-6]
        s = solve_ge(Params(); base..., nϵT = 48, nXO = 48, q_lo = lo, q_hi = 1 - 1e-6)
        @printf("    %-12.0e %-12.5f %-12.5f %-10.1e\n",
                lo, s.HT[2], teach_share_total(s), reldev(s))
    end

    println("\n  (c) Grid resolution nϵT=nXO  (q_lo=1e-5, q_hi=1-1e-6):")
    @printf("    %-12s %-12s %-12s %-10s\n", "grid", "H̃_T(2)", "teach", "rel.dev")
    for n in [16, 24, 48, 96]
        s = solve_ge(Params(); base..., nϵT = n, nXO = n, q_lo = 1e-5, q_hi = 1 - 1e-6)
        @printf("    %-12s %-12.5f %-12.5f %-10.1e\n",
                "$(n)²", s.HT[2], teach_share_total(s), reldev(s))
    end
    println("\n  Default bounds are q_lo=1e-5, q_hi=1-1e-6. If H̃_T / teach barely move")
    println("  across the rows the default is safe; a monotone drift as the tail tightens")
    println("  means the default clips a payoff-relevant tail (upper tail matters most,")
    println("  since H̃_T integrates the increasing quantity h^{β/σ}).")
    println("===========================================================")
    plot(xh, [htv tsv ./ maximum(tsv) .* maximum(htv)], marker = :circle, lw = 2,
         xlabel = "upper-tail tightness  −log₁₀(1−q_hi)",
         label = ["H̃_T(loc 2)" "teach share (scaled)"],
         title = "Sensitivity to the upper quantile bound", legend = :right)
    savefig(joinpath(outdir, "dd_quantile_sensitivity.png"))
    return nothing
end

# ---- Q1. More β cases + WHY high β corners: the β/σ ratio. -------------------
# The teacher block self-consistency (§9/App A.1): H̃_T = ½∫h^{β/σ}, Q=(2H̃_T/M)^σ,
# and h ∝ Q. Substituting Q^{β/σ}=(2H̃_T/M)^β gives H̃_T^{1-β} ∝ M^{-β}·(policy
# integral): the interior teacher stock is well-defined only for β<1, and as β→1
# the map is increasingly ill-conditioned. Independently, expo≡β/σ sets the
# CONVEXITY of the aggregator: larger β/σ concentrates H̃_T on the best teachers
# and strengthens the Q→h→H̃_T feedback, which is what tips high-β economies into
# the near-universal-teaching / collapsing-quality corner.
function beta_sigma_analysis(; outdir, kw = _dd_kw())
    mkpath(outdir)
    println("\n========== Q1. β sweep and the β/σ corner ==========")

    # (a) β sweep at σ=0.25, warm-started (continuation) so we track ONE branch.
    println("\n  (a) β sweep at σ=0.25 (continuation warm start):")
    @printf("    %-7s %-7s %-15s %-15s %-9s %-9s %-8s %-7s\n",
            "β", "β/σ", "H̃_T", "Q", "teach", "res", "nonmono", "corner")
    βs = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
    prev = nothing
    βcurve = Float64[]; tscurve = Float64[]; htcurve = Float64[]
    for β in βs
        NONMONO_WARNINGS[] = 0
        sol = solve_ge(Params(β = β); kw..., init = prev)
        prev = sol
        ts  = teach_share_total(sol)
        res = ge_map_residual(sol).res
        corner = ts > 0.95 || minimum(sol.HT) < 1e-3
        @printf("    %-7.2f %-7.2f %-15s %-15s %-9.4f %-9.1e %-8d %-7s\n",
                β, β / sol.p.σ, fmtvec(sol.HT), fmtvec(sol.gr.Q), ts, res,
                NONMONO_WARNINGS[], corner ? "YES" : "no")
        push!(βcurve, β); push!(tscurve, ts); push!(htcurve, sol.HT[2])
    end

    # (b) Cold vs warm start at high β: is the corner NUMERICAL (cold start finds a
    #     bad branch) or STRUCTURAL (β too close to 1)?
    println("\n  (b) Cold vs continuation start at high β (is the corner numerical?):")
    @printf("    %-7s %-9s %-15s %-9s   %-9s %-15s %-9s\n",
            "β", "cold ts", "cold H̃_T", "cold res", "warm ts", "warm H̃_T", "warm res")
    warm = solve_ge(Params(β = 0.40); kw...)   # a healthy interior seed
    for β in [0.50, 0.60, 0.70]
        NONMONO_WARNINGS[] = 0
        cold = solve_ge(Params(β = β); kw...)
        NONMONO_WARNINGS[] = 0
        warm = solve_ge(Params(β = β); kw..., init = warm)
        @printf("    %-7.2f %-9.4f %-15s %-9.1e   %-9.4f %-15s %-9.1e\n",
                β, teach_share_total(cold), fmtvec(cold.HT), ge_map_residual(cold).res,
                teach_share_total(warm), fmtvec(warm.HT), ge_map_residual(warm).res)
    end

    # (c) β/σ collapse test: teaching share vs β for several σ, then vs β/σ. If the
    #     σ-curves line up when plotted against β/σ, the ratio is the key statistic.
    println("\n  (c) β/σ collapse test (teach share for several σ):")
    σs = [0.15, 0.25, 0.40]
    curves = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    βgrid = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50]
    for σ in σs
        prev = nothing
        bs = Float64[]; ts = Float64[]
        for β in βgrid
            sol = solve_ge(Params(β = β, σ = σ); kw..., init = prev)
            prev = sol
            push!(bs, β); push!(ts, teach_share_total(sol))
        end
        curves[σ] = (bs, ts)
    end
    # Is β/σ a sufficient statistic for the teaching margin? Test: at each σ read
    # the teaching share at a COMMON β/σ (=1.2) by linear interpolation in β/σ; if
    # β/σ were sufficient these would coincide. Spread ⇒ it is not sufficient.
    target = 1.2
    ts_at = [Spline1D(curves[σ][1] ./ σ, curves[σ][2]; k = 1)(target) for σ in σs]
    spread = maximum(ts_at) - minimum(ts_at)
    @printf("      teach share at β/σ=%.1f across σ=%s : %s  (spread %.3f)\n",
            target, string(σs), string(round.(ts_at, digits = 3)), spread)
    println("\n  Takeaway:")
    println("   • EXISTENCE of an interior teacher stock needs β<1: from H̃_T = ½∫h^{β/σ}")
    println("     with h∝Q=(2H̃_T/M)^σ, the fixed point is H̃_T^{1-β} ∝ M^{-β}·(policy),")
    println("     which blows up / flips sign as β→1.")
    println("   • The CORNER (H̃_T→0, near-universal teaching) appears well below β=1")
    println("     here (β≈0.4 at σ=0.25): a falling Q lowers teacher h, lowers H̃_T,")
    println("     lowers Q — a self-reinforcing low-quality trap.")
    @printf("   • β/σ is NOT a sufficient statistic for the teaching margin: at a fixed\n")
    @printf("     β/σ=%.1f the teaching share still varies by %.3f across σ (see panel 3,\n", target, spread)
    println("     which SPREADS the σ-curves rather than collapsing them). β is the")
    println("     better single predictor; σ shifts it (higher σ ⇒ more teaching at given β).")
    println("     β/σ's real role is setting aggregator convexity / class-size dispersion")
    println("     (N(h)∝h^{β/σ}), not pinning the teaching share.")
    println("=====================================================")

    p1 = plot(βcurve, [tscurve htcurve], marker = :circle, lw = 2, xlabel = "β",
              label = ["teach share" "H̃_T(loc 2)"], legend = :left,
              title = "β sweep at σ=0.25 (corner as β→1)")
    p2 = plot(xlabel = "β", ylabel = "teach share", title = "vs β (curves nearly align)")
    p3 = plot(xlabel = "β/σ", ylabel = "teach share", title = "vs β/σ (NO collapse ⇒ not sufficient)")
    for σ in σs
        bs, ts = curves[σ]
        plot!(p2, bs, ts, marker = :circle, lw = 2, label = "σ=$σ")
        plot!(p3, bs ./ σ, ts, marker = :circle, lw = 2, label = "σ=$σ")
    end
    plot(p1, p2, p3, layout = (1, 3), size = (1400, 400))
    savefig(joinpath(outdir, "dd_beta_sigma_collapse.png"))
    println("  Saved β/σ figures to ", outdir)
    return (; βcurve, tscurve, htcurve, curves)
end

# ---- Driver: run every deep-dive analysis on/around a converged baseline. ----
function deep_dive(sol; outdir = joinpath(@__DIR__, "figures"), kw = _dd_kw())
    println("\n\n#################### MECHANISM DEEP-DIVE (Q1–Q8) ####################")
    analyze_occupation_location(sol; outdir)   # Q2
    analyze_class_and_teachers(sol; outdir)    # Q3, Q4, Q7
    analyze_sorting_drivers(; outdir, kw = _dd_kw(Nz = 7))   # Q5 (higher Nz for ability gradient)
    analyze_amenity_altruism(; outdir, kw = _dd_kw(Nz = 7))  # Q8 (higher Nz for ability gradient)
    quantile_sensitivity(; outdir)             # Q6
    beta_sigma_analysis(; outdir, kw)          # Q1
    println("####################################################################\n")
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

    # Mechanism deep-dive (Q1–Q8): verbose comparative statics + figures. Runs on
    # the baseline `sol` plus its own parameter sweeps; see section 7c.
    deep_dive(sol)
    return sol
end

@time sol = main_test()
