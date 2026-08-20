# =============================================================================
# Diagnostics, tests and figures for the GE spatial model in spatial_continuous.jl.
#
# Nothing here feeds back into the solver: every routine either READS a converged
# solution `sol` or re-solves the model at perturbed parameters. Figures are
# written to figures/.
#
# TWO PARAMETERIZATIONS run through this file, and every test states which one it
# uses:
#   `Params()`         — the BASELINE. Both sorting mechanisms on: the moving cost
#                        is denominated in goods (`mcost = MCOST_BENCH`,
#                        `τmove = 0`) and the warm glow is the power kernel at
#                        ψ = 0. Post-migration ability sorting is POSITIVE and
#                        mostly within-branch.
#   `nesting_params()` — the pre-sorting benchmark the two mechanisms nest at
#                        (utility `τmove = 0.20`, log warm glow ψ = 1). Sorting
#                        here is negative and purely compositional, and this is
#                        the point at which the frozen normalizations
#                        HBAR_BENCH / CBAR_BENCH are measured.
# Diagnostics whose EXPECTED reading differs between the two say so in place.
#
# Layout:
#   0  generic helpers
#   1  equilibrium-consistency residuals
#   2  invariance / robustness tests (symmetry, scale, grids, quantile bounds, damping)
#   3  generalized dimensions: the I = 3, L = 3 economy
#   4  read-outs on a converged solution (moments, sorting, migration)
#   5  comparative statics
#   6  mechanism analyses and parameter sweeps
#   7  figures
#   8  run
#
# Run:  julia --project=.. spatial_continuous_script.jl
# =============================================================================

using Printf
using Plots

include(joinpath(@__DIR__, "spatial_continuous.jl"))

gr()  # Plots.jl GR backend (unrelated to the solver's `gr` grids NamedTuple)

const FIGDIR = joinpath(@__DIR__, "figures")

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

"Marker/colour cycles reused across the figures."
const MARKERS = [:circle, :square, :diamond, :utriangle, :star5, :hexagon]
mk(i) = MARKERS[mod1(i, length(MARKERS))]

"""
    node_mass(x, F)

Probability mass carried by each node of a quadrature grid under the CDF `F`:
node k gets the mass of the midpoint cell around it, with the two end cells
extended to the tails, so the masses sum to EXACTLY 1. This is what turns a shock
grid into a weighted node sample (see `teacher_logh_samples`).

Use this rather than a trapezoid rule on the density: the shock grids are
QUANTILE grids, so their spacing is wildly non-uniform and the trapezoid rule
overstates the sparse upper tail badly (at nϵT = 32 it inflates the total mass by
~7 % and the teaching share — which lives in that tail — by ~35 %).
"""
function node_mass(x, F)
    n = length(x)
    m = zeros(n)
    lo = 0.0
    for k in 1:n
        hi = k == n ? 1.0 : F((x[k] + x[k+1]) / 2)
        m[k] = max(hi - lo, 0.0)
        lo = hi
    end
    return m
end

"Cumulative trapezoid integral of `f` sampled on the grid `x` (∫ from x[1]),
 same length as `x`. Turns a density into a CDF for range-trimming."
function cumtrapz(x, f)
    F = zeros(length(x))
    @inbounds for k in 2:length(x)
        F[k] = F[k-1] + (f[k] + f[k-1]) * (x[k] - x[k-1]) / 2
    end
    return F
end

"Unweighted mean of `f` over an index range."
mean_over(f, rng) = sum(f, rng) / length(rng)

"Weighted covariance and OLS slope of y on x under weights w."
function _wcov(x, y, w)
    sw = sum(w)
    mx = sum(w .* x) / sw
    my = sum(w .* y) / sw
    return sum(w .* (x .- mx) .* (y .- my)) / sw
end
_wslope(x, y, w) = _wcov(x, y, w) / max(_wcov(x, x, w), 1e-300)

"Index range of a shock grid trimmed to the central mass, `F` its CDF. The
 quadrature grids run out to the 1e-5 / 1−1e-6 quantiles, where there is no
 visible mass, so plots against a shock axis need trimming to stay readable."
function bulk_range(x, F; lo = 0.02, hi = 0.98)
    q = F.(x)
    i1 = something(findfirst(≥(lo), q), 1)
    i2 = something(findlast(≤(hi), q), length(x))
    return i1 < i2 ? (i1:i2) : eachindex(x)
end

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

"""
    ability_neutrality_dev(sol)

Proposition A: with NO goods-denominated moving cost (`p.mcost ≡ 0`) the location
probabilities π_{l'|l} are invariant to the idiosyncratic shock index k — the
shock scales income multiplicatively, so it shifts every location's log
consumption by the same amount and cancels out of the logit.  This is exactly the
homogeneity that makes sorting a pure OCCUPATIONAL-COMPOSITION effect with no
within-branch ability gradient.

Reported, not asserted, and the sign of the reading flips with the
parameterization:
  • `nesting_params()` (and any `mcost = 0` variant, including a pure ψ sweep —
    that mechanism works through `zi`, not `k`): ~1e-12.
  • the BASELINE, where `mcost = MCOST_BENCH ≠ 0`: ≈ 0.4. Large is the PASS
    condition here — it is the headline evidence the goods moving cost is live.
"""
function ability_neutrality_dev(sol)
    (; Pi_T, Pi_O) = sol
    dT = maximum(abs, Pi_T .- Pi_T[:, :, :, 1:1, :])
    dO = maximum(abs, Pi_O .- Pi_O[:, :, :, 1:1, :])
    return max(dT, dO)
end

function report_checks(sol)
    println("\n========== Equilibrium-consistency checks ==========")

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

    # π̄ and the Φ eigen-residual come from a QUADRATURE integral (integrate_choice's
    # spline + tail correction), not an exact sum, so they are only accurate to the
    # current (nϵT, nXO) resolution — they shrink under grid refinement (see
    # grid_refinement_test) rather than holding at machine precision. Tolerance is
    # the ballpark the solver's own docstring quotes for its tail correction
    # (~1e-3 at nϵT, nXO ~ 64).
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

    # Reported, NOT a pass/fail: ≈0 means location choice is ability-neutral
    # within a branch (nesting_params(), or any pure ψ variant); a large value
    # means the goods moving cost has broken that homogeneity, which is the
    # mechanism working as designed and the BASELINE reading.
    and = ability_neutrality_dev(sol)
    @printf("\n  Prop. A: max_k |π_{l'|l}(k) − π_{l'|l}(k=1)| : %.2e   (%s)\n", and,
            maximum(sol.p.mcost) > 0 ? "mcost ≠ 0 (baseline) ⇒ expected to be LARGE" :
                                       "mcost = 0 ⇒ expected ≈ 0")

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
    println("====================================================")
end

# -----------------------------------------------------------------------------
# 2. Invariance and robustness tests.
# -----------------------------------------------------------------------------

"Identical locations ⇒ identical aggregates and identical ability distributions."
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

"Doubling total population doubles the extensive aggregates and leaves the
 intensive ones (Q, t, G_l) alone."
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

"""
    grid_refinement_test(; grids, ...)

Refining the (nϵT, nXO) quadrature grids must leave the equilibrium numerically
unchanged — i.e. the reported solution is not a grid artifact. This is the
continuous solver's analogue of a tremble-invariance test: it confirms that a
purely numerical device (grid resolution) is not moving the economics.

Solves are warm-started from the coarser grid's solution (continuation), both for
speed and so all resolutions sit on the same equilibrium branch.
"""
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

"""
    quantile_sensitivity(; outdir, ...)

Sensitivity to where the shock grids are TRUNCATED (`q_lo`, `q_hi`) and to their
resolution. If H̃_T and the teaching share barely move across the rows the
defaults are safe; a monotone drift as a tail tightens means the default clips a
payoff-relevant tail. The upper tail matters most, since H̃_T integrates the
increasing quantity h^{β/σ}.
"""
function quantile_sensitivity(; outdir = FIGDIR, Nz = 3, damping = 0.3, tol = 1e-5,
                              maxit = 200, hh_maxit = 600)
    mkpath(outdir)
    base = (; Nz, damping, tol, maxit, hh_tol = 1e-6, hh_maxit, verbose = false)
    println("\n========== Quantile-bound & grid sensitivity ==========")
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
    println("\n  Defaults are q_lo=1e-5, q_hi=1-1e-6.")
    println("=======================================================")
    plot(xh, [htv tsv ./ maximum(tsv) .* maximum(htv)], marker = :circle, lw = 2,
         xlabel = "upper-tail tightness  −log₁₀(1−q_hi)",
         label = ["H̃_T(loc 2)" "teach share (scaled)"],
         title = "Sensitivity to the upper quantile bound", legend = :right)
    savefig(joinpath(outdir, "continuous_quantile_sensitivity.png"))
    return nothing
end

"How far the outer fixed point gets, as a function of the update step size."
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
# 3. Generalized dimensions: an I = 3, L = 3 economy.
#
# Two jobs. First a SMOKE TEST that nothing in the solver is hard-wired to 2×2:
# with I = 3 there are two non-teaching occupations, so the baseline's ψ ≠ 1 also
# exercises the `EΘpow` branch of `child_Ef` (the Jensen-correct
# E[Θ_{i*}^{−(1−ψ)} | X_O*]) that the 2-occupation baseline never reaches.
#
# Second, and the reason the three locations are NOT symmetric, an economics
# experiment. Both non-base locations are attractive, but for DIFFERENT reasons:
#
#     loc 1  the plain location: base teaching wage, no amenity
#     loc 2  AMENITY pull      (B₂ > 0,  κ₂ = κ₁)  — an EXOGENOUS pull
#     loc 3  TEACHING-WAGE pull (κ₃ > κ₁, B₃ = 0)  — an ENDOGENOUS pull
#
# The amenity enters every agent's location value directly and changes nothing
# about the local technology, so its only effect on schools runs through the
# DENOMINATOR of Q_l = (2H̃_T,l/M_l)^σ: every extra resident is another student.
# The teaching-wage shifter has no direct effect on the location value at all —
# it raises teacher income, draws more and better teachers, and so raises H̃_T and
# Q_l; the improved school is then what pulls everyone else in through the
# altruism term, paid for by a much higher local tax rate.
#
# So the pair separates a pull that DILUTES schools from a pull that IS better
# schools. The read-out below is the comparison, and the headline is that the
# amenity location can be the most populous and still have the worst schools.
#
# The goods moving cost is RE-DERIVED for this economy rather than inherited.
# `mcost` is an ABSOLUTE goods amount, so it is not scale-free: this economy's own
# C̄ differs from the 2×2 benchmark's, and inheriting MCOST_BENCH would misprice
# the move in whichever direction the two scales differ — too cheap to bind, or
# cheap enough at the benchmark and yet negative-consumption here (which is what
# happens at a lower κ, where the smooth_log barrier starts binding). Solving once
# with the utility cost, reading that economy's own C̄ and h̄, then re-denominating
# is exactly the procedure that produced the baseline's constants, and is what any
# new parameterization at a different scale should do.
# -----------------------------------------------------------------------------
function generalized_dims_test(; Nz = 3, nϵT = 48, nXO = 48, damping = 0.3, maxit = 80,
                                κ = [0.90, 0.90, 1.20], B = [0.0, 0.25, 0.0])
    println("\n========== Generalized dimensions (I=3, L=3) ==========")
    kw = (; Nz, nϵT, nXO, damping, tol = 1e-3, maxit, hh_tol = 1e-5,
            hh_maxit = 300, verbose = false)
    @printf("  loc 1 = plain (κ=%.2f, B=%.2f) | loc 2 = AMENITY (κ=%.2f, B=%.2f) | loc 3 = TEACHING WAGE (κ=%.2f, B=%.2f)\n",
            κ[1], B[1], κ[2], B[2], κ[3], B[3])

    # Pass 1: utility moving cost, log warm glow — the nesting parameterization at
    # these dimensions, used only to read off this economy's own scale.
    p_nest = Params(3, 3; A = [NaN, 1.5, 2.0], κ, B, τmove_off = TAUMOVE_BENCH,
                    mcost_off = 0.0, ψ = 1.0, href = 1.0)
    s_nest = solve_ge(p_nest; kw...)
    href3, Cbar3 = mean_child_h(s_nest), mean_consumption(s_nest)
    @printf("  own-scale normalizations: h̄₃ = %.6f   C̄₃ = %.6f\n", href3, Cbar3)

    # Pass 2: the baseline structure (goods cost + power kernel) at this scale.
    p = redenominate_move_cost(with_params(p_nest; ψ = Params().ψ, href = href3), Cbar3)
    @printf("  re-denominated m₃ = %.6f  (vs the inherited MCOST_BENCH = %.6f — the wrong absolute size at this scale)\n",
            p.mcost[1, 2], MCOST_BENCH)
    sol = solve_ge(p; kw...)
    report_ge(sol)

    # ---- The amenity-vs-teaching-wage comparison. ---------------------------
    tm = teacher_moments(sol)
    sm = sorting_metrics(sol)
    println("  Amenity (loc 2) vs teaching-wage (loc 3) pull:")
    @printf("    %-6s %-7s %-7s %-9s %-9s %-9s %-9s %-9s %-9s %-9s\n",
            "loc", "κ", "B", "M", "Q", "t", "#teach", "class N̄", "teach sh", "mean z")
    tshare = tm.Tcount ./ (tm.Tcount .+ tm.Ocount)
    for l in 1:3
        @printf("    %-6d %-7.2f %-7.2f %-9.4f %-9.4f %-9.4f %-9.4f %-9.3f %-9.4f %-9.4f\n",
                l, p.κ[l], p.B[l], sol.M[l], sol.gr.Q[l], sol.t[l],
                tm.Tcount[l], tm.classsize[l], tshare[l], sm.mz_work[l])
    end
    M, Q, t = sol.M, sol.gr.Q, sol.t
    @printf("\n    headcount   M₂/M₁ = %.2f (amenity)      M₃/M₁ = %.2f (teaching wage)\n",
            M[2] / M[1], M[3] / M[1])
    @printf("    schools     Q₂/Q₁ = %.2f (amenity)      Q₃/Q₁ = %.2f (teaching wage)\n",
            Q[2] / Q[1], Q[3] / Q[1])
    @printf("    tax rate    t₂ = %.4f                 t₃ = %.4f   (teacher wage bill / taxable income)\n",
            t[2], t[3])
    @printf("    class size  N̄₂ = %.1f                   N̄₃ = %.1f\n",
            tm.classsize[2], tm.classsize[3])
    chk(name, cond) = @printf("      %-52s %s\n", name, mark(cond))
    println("\n    Sign checks (the point of making locations 2 and 3 different):")
    chk("the amenity draws population (M₂ > M₁)", M[2] > M[1])
    chk("… but does NOT buy schools (Q₂ ≲ Q₁)", Q[2] < 1.02 * Q[1])
    chk("the teaching wage buys the best schools (Q₃ = max Q)", Q[3] ≈ maximum(Q))
    chk("… and the highest local tax rate (t₃ = max t)", t[3] ≈ maximum(t))
    chk("high-z agents sort to the school-quality location", sm.mz_work[3] ≈ maximum(sm.mz_work))
    println("\n    Reading: an amenity and a teaching-wage premium both attract people, but")
    println("    only one of them attracts them BY improving schools. The amenity enters")
    println("    the location value directly and its residents only enlarge the Q")
    println("    denominator, so loc 2 can be the most populous location and still have")
    println("    the worst schools; κ₃ works entirely through H̃_T, so loc 3's pull comes")
    println("    with smaller classes, better teachers, higher-ability residents — and a")
    println("    tax rate several times loc 2's to pay the teacher wage bill.")

    # ---- Structural identities at these dimensions. -------------------------
    ph  = phi_stationarity_residual(sol)
    pir = pi_rowsum_dev(sol)
    pbr = pibar_rowsum_dev(sol)
    qd  = q_consistency_dev(sol)
    Mtot_dev = abs(sum(sol.M) - sol.p.Mtot)
    feas = verify_solution(sol)
    # ph.eig_res, ph.row_dev, pbr are quadrature-accuracy (not exact) — see
    # report_checks; pir, qd, Mtot_dev are exact structural identities.
    ok = ph.eig_res < 1e-3 && ph.row_dev < 1e-3 && pir < 1e-6 &&
         pbr < 1e-3 && qd < 1e-10 && Mtot_dev < 1e-8 &&
         policies_in_bounds(sol.hh) && feas
    @printf("\n  I=%d, L=%d   Φ eig-res=%.2e   Σπ_{l'|l}=1 dev=%.2e   Σ_l M_l=Mtot dev=%.2e   [%s]\n",
            sol.gr.I, sol.gr.L, ph.eig_res, pir, Mtot_dev, mark(ok))
    println("=========================================================")
    return (; sol, tm, sm, ok)
end

# -----------------------------------------------------------------------------
# 4. Read-outs on a converged solution.
#
# All of these are pure functions of a converged `sol`, re-using the solver's own
# choice_maps / integrate_choice so their integrals match aggregates() exactly.
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

"Teacher-HC and class-size moments by WORK location.
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

"""
    sorting_metrics(sol)

Mean ability z by BIRTH and by WORK location, and scalar sorting indices (spread
in mean z across locations). `Ψ` is the post-migration (work-location) ability
mass.

`idx_work` is the |max − min| SPREAD, which hides the sign, and the sign is the
whole point. At `nesting_params()` sorting is NEGATIVE (≈ −0.041: the
high-amenity, high-κ location attracts LOW-z agents, because location choice is
ability-neutral within a branch and teaching — the branch that moves — selects on
low z). At the BASELINE the two mechanisms flip it POSITIVE (≈ +0.077); neither
alone suffices to flip it (warm glow alone ≈ +0.050, goods moving cost alone
≈ −0.011). `idx_work_signed` = mz_work[L] − mz_work[1] keeps the sign and is the
headline number.
"""
function sorting_metrics(sol)
    (; Φ, πbar, gr) = sol
    (; z, Nz, L) = gr
    Ψ = [sum(Φ[zi, l0] * πbar[l0, zi, l] for l0 in 1:L) for zi in 1:Nz, l in 1:L]
    meanz(w) = sum(z[zi] * w[zi] for zi in 1:Nz) / sum(w)
    mz_birth = [meanz(@view Φ[:, l]) for l in 1:L]
    mz_work  = [meanz(@view Ψ[:, l]) for l in 1:L]
    return (; mz_birth, mz_work, Ψ,
              idx_birth = maximum(mz_birth) - minimum(mz_birth),
              idx_work  = maximum(mz_work)  - minimum(mz_work),
              idx_birth_signed = mz_birth[L] - mz_birth[1],
              idx_work_signed  = mz_work[L]  - mz_work[1])
end

"""
    within_branch_gradient(sol; lo = 1, hi = L)

The within-branch ability gradient

    g_b = σν · ∂ log( π^b_{hi|l}(z) / π^b_{lo|l}(z) ) / ∂ log z ,   b ∈ {T, O},

i.e. how much the UTILITY gap between working in `hi` and in `lo` rises with the
agent's own ability, holding the occupational branch fixed. Estimated as the
Φ-weighted OLS slope of log(π_hi/π_lo) on log z over the `gr.z` grid, averaged
over shock nodes and genders (by Prop. A every shock node gives the same answer
whenever `mcost = 0`; with the baseline's `mcost ≠ 0` they genuinely differ and
the average is the summary).

This is the statistic that separates the two sorting channels. At
`nesting_params()` g ≈ 0 in BOTH branches and all of `idx_work` comes from
occupational composition. Turning the mechanisms on, one at a time and then both,
gives

    parameterization                g_T      g_O
    nesting (both off)              0.0006   0.0006
    warm glow only (ψ = 0, τmove)   0.059    0.060
    moving cost only (mcost, ψ = 1) 0.115    0.144
    BASELINE (both)                 0.171    0.200

so the two are additive to first order (0.059 + 0.115 ≈ 0.171).

Returns per-birth-location vectors `(; gT, gO)`.
"""
function within_branch_gradient(sol; lo = 1, hi = sol.gr.L)
    (; Φ, gr, p) = sol
    (; z, Nz, L, nϵT, nXO) = gr
    L < 2 && return (; gT = fill(NaN, L), gO = fill(NaN, L))
    logz = log.(z)
    gT = zeros(L); gO = zeros(L)
    for l in 1:L
        w = collect(@view Φ[:, l])
        yT = [mean_over(g -> mean_over(k -> log(sol.Pi_T[g, l, zi, k, hi] /
                                                sol.Pi_T[g, l, zi, k, lo]), 1:nϵT), 1:2)
              for zi in 1:Nz]
        yO = [mean_over(g -> mean_over(k -> log(sol.Pi_O[g, l, zi, k, hi] /
                                                sol.Pi_O[g, l, zi, k, lo]), 1:nXO), 1:2)
              for zi in 1:Nz]
        gT[l] = p.σν * _wslope(logz, yT, w)
        gO[l] = p.σν * _wslope(logz, yO, w)
    end
    return (; gT, gO)
end

"""
    sorting_decomposition(sol; lo = 1, hi = L)

Split the ability sorting into a WITHIN-BRANCH and a COMPOSITION term:

    Cov_z(z, π_hi(z)) = Σ_b θ̄_b Cov_z(z, π^b_hi(z))   (within-branch gradient)
                      + Σ_b π̄^b_hi Cov_z(z, θ_b(z))   (occupational composition)
                      + interaction

with θ_b(z) the branch share at ability z and π^b_hi(z) the branch-conditional
probability of working in `hi`. All covariances are Φ-weighted over the z grid at
a fixed birth location; the aggregate row weights birth locations by their mass.

This decomposition, not `idx_work` itself, is the result. At `nesting_params()`
term 1 ≈ 0 and term 2 (negative) carries essentially all of `idx_work`. At the
BASELINE term 1 is ≈ +0.014 and term 2 is still ≈ −0.005, so the within-branch
share of the total reads ≈ 150 %: the mechanisms did not merely add sorting, they
overturned a composition effect that still points the other way.
Returns per-birth-location vectors plus mass-weighted totals.
"""
function sorting_decomposition(sol; lo = 1, hi = sol.gr.L)
    (; hh, Φ, gr) = sol
    (; z, Nz, L, nϵT, nXO) = gr
    within = zeros(L); comp = zeros(L); inter = zeros(L); tot = zeros(L)
    zeroT = zeros(nϵT); zeroO = zeros(nXO)
    for l0 in 1:L
        θT = zeros(Nz); θO = zeros(Nz); pT = zeros(Nz); pO = zeros(Nz)
        for zi in 1:Nz
            mT = mO = qT = qO = 0.0
            for g in 1:2
                cm = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
                iT, iO = integrate_choice(ones(nϵT), ones(nXO), cm, gr)   # branch masses
                jT, _  = integrate_choice(collect(@view sol.Pi_T[g, l0, zi, :, hi]), zeroO, cm, gr)
                _,  jO = integrate_choice(zeroT, collect(@view sol.Pi_O[g, l0, zi, :, hi]), cm, gr)
                mT += 0.5 * iT; mO += 0.5 * iO; qT += 0.5 * jT; qO += 0.5 * jO
            end
            θT[zi] = mT; θO[zi] = mO
            pT[zi] = mT > 0 ? qT / mT : 0.0        # π_hi conditional on the branch
            pO[zi] = mO > 0 ? qO / mO : 0.0
        end
        w = collect(@view Φ[:, l0]); sw = sum(w)
        wmean(v) = sum(w .* v) / sw
        θTb, θOb = wmean(θT), wmean(θO)
        pTb, pOb = wmean(pT), wmean(pO)
        π_hi = θT .* pT .+ θO .* pO
        tot[l0]    = _wcov(z, π_hi, w)
        within[l0] = θTb * _wcov(z, pT, w) + θOb * _wcov(z, pO, w)
        comp[l0]   = pTb * _wcov(z, θT, w) + pOb * _wcov(z, θO, w)
        inter[l0]  = tot[l0] - within[l0] - comp[l0]
    end
    mass = [sum(@view Φ[:, l0]) for l0 in 1:L]; mass ./= sum(mass)
    return (; within, comp, inter, tot,
              within_tot = sum(mass .* within), comp_tot = sum(mass .* comp),
              inter_tot  = sum(mass .* inter),  tot_tot  = sum(mass .* tot))
end

"""
    location_choice_by_ability(sol, l0)

Branch-resolved location-choice probabilities out of birth location `l0`:
`πT[zi, lp]` and `πO[zi, lp]` are the probability of WORKING in `lp` conditional
on ability node `zi` AND on the occupational branch. The per-shock `Pi_T` / `Pi_O`
are integrated over the occupation margin (so the branch weights are the
endogenous teaching probabilities) and pooled over genders; each row is then
normalised by its branch mass, so both matrices have rows summing to 1.

The gap between the two is what `within_branch_gradient` summarises as a slope:
πT and πO answer "where does a teacher / a non-teacher of ability z go?", while
`sol.πbar` is their branch-share-weighted average.
"""
function location_choice_by_ability(sol, l0)
    (; hh, gr) = sol
    (; Nz, L, nϵT, nXO) = gr
    zeroT = zeros(nϵT); zeroO = zeros(nXO)
    πT = zeros(Nz, L); πO = zeros(Nz, L)
    for zi in 1:Nz
        mT = mO = 0.0
        for g in 1:2
            cm = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
            iT, iO = integrate_choice(ones(nϵT), ones(nXO), cm, gr)
            mT += 0.5 * iT; mO += 0.5 * iO
            for lp in 1:L
                jT, _ = integrate_choice(collect(@view sol.Pi_T[g, l0, zi, :, lp]), zeroO, cm, gr)
                _, jO = integrate_choice(zeroT, collect(@view sol.Pi_O[g, l0, zi, :, lp]), cm, gr)
                πT[zi, lp] += 0.5 * jT; πO[zi, lp] += 0.5 * jO
            end
        end
        mT > 0 && (πT[zi, :] ./= mT)
        mO > 0 && (πO[zi, :] ./= mO)
    end
    return (; πT, πO)
end

"Per-ability probability of leaving the birth location, population-weighted over
 birth locations and also resolved per birth location."
function mobility_by_ability(sol)
    (; Φ, πbar, gr) = sol
    (; Nz, L) = gr
    out_by_birth = [sum(πbar[l0, zi, l] for l in 1:L if l != l0) for zi in 1:Nz, l0 in 1:L]
    outflow = [sum(Φ[zi, l0] * out_by_birth[zi, l0] for l0 in 1:L) / sum(@view Φ[zi, :])
               for zi in 1:Nz]
    return (; outflow, out_by_birth)
end

"Weighted (log h, weight) sample of teachers WORKING in location `lp`, for a
 density estimate. Node weight = cell mass × P(ε in the node's cell) × P(teach|ε)
 × π(→lp)."
function teacher_logh_samples(sol, lp)
    (; hh, Φ, gr) = sol
    (; Nz, L, nϵT, ϵTgrid, dT) = gr
    dϵ = node_mass(ϵTgrid, ε -> cdf(dT, ε))     # exact per-node probability mass on the ε_T grid
    logh = Float64[]; w = Float64[]
    for l0 in 1:L, zi in 1:Nz, g in 1:2
        cell = 0.5 * Φ[zi, l0]
        cell == 0.0 && continue
        cm = choice_maps(hh.WT[g, l0, zi, :], hh.WO[g, l0, zi, :], g, gr)
        for k in 1:nϵT
            tw = teach_wt(cm, ϵTgrid[k])
            tw <= 0.0 && continue
            wt = cell * dϵ[k] * tw * sol.Pi_T[g, l0, zi, k, lp]
            wt <= 0.0 && continue
            push!(logh, log(hh.hT[g, l0, zi, k])); push!(w, wt)
        end
    end
    return logh, w
end

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

"""
    report_sorting(sol)

The sorting battery in one block: Prop. A deviation, the signed sorting index,
the within-branch gradient per branch and birth location, and the
within-vs-composition decomposition. This is the read-out both sorting mechanisms
are judged against.
"""
function report_sorting(sol)
    (; L) = sol.gr
    sm = sorting_metrics(sol)
    wg = within_branch_gradient(sol)
    dc = sorting_decomposition(sol)
    println("\n========== Sorting diagnostics ==========")
    @printf("  Prop. A deviation (π invariant to shock node) : %.3e\n", ability_neutrality_dev(sol))
    @printf("  mean z by work location : %s\n", fmtvec(sm.mz_work))
    @printf("  idx_work (spread) = %+.5f    idx_work_signed (loc %d − loc 1) = %+.5f\n",
            sm.idx_work, L, sm.idx_work_signed)
    println("\n  Within-branch ability gradient  g_b = σν·dlog(π_L/π_1)/dlog z")
    println("  (≈0 at nesting_params(); ≈0.17 / 0.20 at the baseline, which runs both mechanisms):")
    @printf("    %-14s %-12s %-12s\n", "birth loc", "g_T", "g_O")
    for l in 1:L
        @printf("    %-14d %+-12.5f %+-12.5f\n", l, wg.gT[l], wg.gO[l])
    end
    println("\n  Cov_z(z, π_L(z)) decomposition  (within-branch | composition | interaction):")
    @printf("    %-14s %-14s %-14s %-14s %-14s\n", "birth loc", "within", "composition", "interaction", "total")
    for l in 1:L
        @printf("    %-14d %+-14.6f %+-14.6f %+-14.6f %+-14.6f\n",
                l, dc.within[l], dc.comp[l], dc.inter[l], dc.tot[l])
    end
    @printf("    %-14s %+-14.6f %+-14.6f %+-14.6f %+-14.6f\n",
            "mass-wtd", dc.within_tot, dc.comp_tot, dc.inter_tot, dc.tot_tot)
    share_within = abs(dc.tot_tot) > 1e-12 ? dc.within_tot / dc.tot_tot : NaN
    @printf("  ⇒ within-branch share of total sorting : %.1f%%\n", 100 * share_within)
    println("=========================================")
    return (; sm, wg, dc)
end

# -----------------------------------------------------------------------------
# 5. Comparative statics.
# -----------------------------------------------------------------------------
"""
    comparative_statics(; ...)

One-parameter deviations from the BASELINE, reported in two tables: the
per-location equilibrium objects (H̃_T, M, Q, t) and the behavioural summaries
(teaching share, outflow, sorting index, class size, teacher share of employment).

The last three rows turn the sorting mechanisms off — one at a time, then both —
so the table doubles as a mechanism ablation.

`high-move-cost` scales the GOODS cost rather than reinstating a utility one:
1.5× MCOST_BENCH is the ceiling the smooth_log barrier allows, since 2× drives
mover consumption negative (see verify_solution's docstring).
"""
function comparative_statics(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.2, maxit = 120,
                              hh_maxit = 600)
    println("\n========== Comparative statics (one-parameter deviations) ==========")
    τb = TAUMOVE_BENCH
    cases = [
        (label = "baseline",        p = Params()),
        (label = "high-β",          p = Params(β = 0.50)),
        (label = "zero-β",          p = Params(β = 0.0)),
        (label = "high-altruism",   p = Params(λ = 0.90)),
        (label = "low-altruism",    p = Params(λ = 0.10)),
        (label = "zero-altruism",   p = Params(λ = 0.0)),   # ε=0 AND λ=0 ⇒ no intergenerational coupling
        (label = "high-move-cost",  p = Params(mcost = 1.5 .* Params().mcost)),
        (label = "no-move-cost",    p = Params(mcost = zeros(2, 2))),
        (label = "low-σν",          p = Params(σν = 0.10)),
        (label = "amenity-loc2",    p = Params(B = [0.0, 0.50])),
        (label = "high-teach-wage", p = Params(κ = [1.20, 1.40])),
        (label = "gender-wedge",    p = Params(τω = [0.0 0.0; 0.0 0.30])),
        (label = "high-persistence",p = Params(ρz = 0.97)),
        (label = "low-persistence", p = Params(ρz = 0.40)),
        # --- mechanism ablations (see nesting_params) ---
        (label = "log-warm-glow",   p = Params(ψ = 1.0, href = 1.0)),                    # warm glow off
        (label = "utility-move-cost", p = Params(τmove = [0.0 τb; τb 0.0],
                                                 mcost = zeros(2, 2))),                  # goods cost off
        (label = "nesting (both off)", p = nesting_params()),
    ]

    println("\n  Equilibrium objects by location:")
    @printf("  %-19s %-8s %-6s %-18s %-18s %-18s %-18s\n",
            "case", "res", "conv?", "H̃_T", "M", "Q", "t")
    sols  = Dict{String,Any}()
    stuck = String[]
    for c in cases
        sol = solve_ge(c.p; Nz, nϵT, nXO, damping, tol = 1e-3, maxit,
                       hh_tol = 1e-5, hh_maxit, verbose = false)
        sols[c.label] = sol
        res  = ge_map_residual(sol).res
        conv = res < 5e-3
        conv || push!(stuck, c.label)
        @printf("  %-19s %.1e  %-6s %-18s %-18s %-18s %-18s\n",
                c.label, res, conv ? "yes" : "NO", fmtvec(sol.HT), fmtvec(sol.M),
                fmtvec(sol.gr.Q), fmtvec(sol.t))
    end

    println("\n  Behaviour and composition:")
    @printf("  %-19s %-8s %-8s %-12s %-22s %-18s %-7s\n",
            "case", "teach", "out", "idx_work±", "class size N̄", "teacher share", "bounds")
    for c in cases
        sol = sols[c.label]
        tm  = teacher_moments(sol)
        sm  = sorting_metrics(sol)
        tshare = tm.Tcount ./ (tm.Tcount .+ tm.Ocount)     # teachers as a share of local employment
        @printf("  %-19s %-8.4f %-8.4f %+-12.5f %-22s %-18s %-7s\n",
                c.label, teach_share_total(sol), outflow_rate(sol), sm.idx_work_signed,
                fmtvec(tm.classsize), fmtvec(tshare), mark(policies_in_bounds(sol.hh)))
    end
    if !isempty(stuck)
        println("\n  ⚠ did not converge in the test budget (read their aggregates with care): ",
                join(stuck, ", "))
    end

    println("\n  Sign checks:")
    base = sols["baseline"]
    chk(name, cond) = @printf("    %-44s %s\n", name, cond ? "ok" : "FAIL")
    chk("amenity in loc 2 raises M₂", sols["amenity-loc2"].M[2] > base.M[2])
    chk("high move cost lowers outflow", outflow_rate(sols["high-move-cost"]) < outflow_rate(base))
    chk("no move cost raises outflow", outflow_rate(sols["no-move-cost"]) > outflow_rate(base))
    chk("higher teaching wage raises Q", minimum(sols["high-teach-wage"].gr.Q .- base.gr.Q) > 0)
    if ge_map_residual(sols["low-σν"]).res < 5e-3
        chk("low σν lowers outflow", outflow_rate(sols["low-σν"]) < outflow_rate(base))
    else
        @printf("    %-44s %s\n", "low σν lowers outflow", "SKIP (not converged)")
    end

    # The mechanism ablation. The baseline must sort POSITIVELY on ability and the
    # nesting parameterization negatively, with each single-mechanism case between
    # them — the property the baseline was re-parameterized for, so it is a hard
    # sign check, not just a reported number. The within-branch gradient is printed
    # alongside because it is the cleaner statistic (no composition noise) and it
    # shows the two mechanisms are additive to first order:
    # g(moving cost) + g(warm glow) ≈ g(baseline).
    println("\n  Mechanism ablation (sorting on ability):")
    @printf("    %-26s %-12s %-11s %-11s\n", "case", "idx_work±", "g_T", "g_O")
    abl = ["nesting (both off)", "log-warm-glow", "utility-move-cost", "baseline"]
    idxs = Dict{String,Float64}()
    for lab in abl
        sm = sorting_metrics(sols[lab]); wg = within_branch_gradient(sols[lab])
        idxs[lab] = sm.idx_work_signed
        @printf("    %-26s %+-12.5f %+-11.5f %+-11.5f\n",
                lab, sm.idx_work_signed, wg.gT[1], wg.gO[1])
    end
    ib, inest = idxs["baseline"], idxs["nesting (both off)"]
    chk("baseline sorts POSITIVELY on ability", ib > 0)
    chk("nesting sorts negatively on ability", inest < 0)
    chk("each single mechanism sits between the two",
        all(inest < idxs[lab] < ib for lab in ("log-warm-glow", "utility-move-cost")))

    chk("every case in bounds",
        all(policies_in_bounds(s.hh) for s in values(sols)))
    println("====================================================================")
    return sols
end

# -----------------------------------------------------------------------------
# 6. Mechanism analyses and parameter sweeps.
#
# Two sorting mechanisms, BOTH ON at the baseline, each nesting the original model
# exactly:
#   warm glow    power kernel  f(h') = ((h'/h̄)^{1−ψ} − 1)/(1−ψ)   [ψ = 1 nests]
#   moving cost  goods-denominated  m_{l,l'}                      [m = 0 nests]
# `nesting_params()` turns both off at once and is the regression point.
# -----------------------------------------------------------------------------

"Solve settings for the mechanism tests. Nz = 5 so the within-branch gradient has
 enough z nodes to fit a slope; grids stay small so a sweep is ~seconds a point."
_mech_kw(; Nz = 5, nϵT = 32, nXO = 32, damping = 0.3, tol = 1e-5, maxit = 200,
          hh_maxit = 800) =
    (; Nz, nϵT, nXO, damping, tol, maxit, hh_tol = 1e-6, hh_maxit, verbose = false)

"Solve settings for the many-solve sweeps below: small enough (Nz = 3, 32²) that
 a sweep point costs seconds. Every driver takes a keyword override."
_sweep_kw(; Nz = 3, nϵT = 32, nXO = 32, damping = 0.3, tol = 1e-4, maxit = 150,
           hh_maxit = 500) =
    (; Nz, nϵT, nXO, damping, tol, maxit, hh_tol = 1e-5, hh_maxit, verbose = false)

"""
    href_invariance_test(sol; href = 0.25, nsweep = 60)

At ψ = 1 the kernel is log(h'/href), so href only shifts f by the CONSTANT
−log href. That constant propagates into Λ as exactly −λ·log href, uniform across
(z, l'), hence cancels out of every choice probability and every FOC. Changing
href must therefore leave the policies e, h and the whole altruism GRADIENT
Λ − mean(Λ) alone, and shift Λ itself by −λ log href.

Run at the HOUSEHOLD level on `sol`'s grids, with a fixed sweep count: the
household and GE stopping rules compare successive W's, which are not invariant
to a level shift, so two GE solves would stop at slightly different points and
mask the identity behind solver tolerance.

The identity is exact in exact arithmetic but only holds up to the QUADRATURE
MASS DEFECT here: Λ integrates f against the choice density, and that density
integrates to 1 − O(1e-3) at these grid sizes (see `report_checks`'s
quadrature-accuracy block), so a constant c in f comes back as c·(1 − defect).
The tolerance is therefore scaled by the measured defect; the test still catches
any *first-order* leak of href into behaviour, which is what it exists to rule out.
"""
function href_invariance_test(sol; href = 0.25, nsweep = 60)
    println("\n========== href-invariance test (ψ = 1) ==========")
    (; gr) = sol
    p1 = with_params(sol.p; ψ = 1.0, href = 1.0)
    p2 = with_params(sol.p; ψ = 1.0, href = href)
    h1 = solve_household(p1, gr; tol = 0.0, maxit = nsweep)
    h2 = solve_household(p2, gr; tol = 0.0, maxit = nsweep)
    # Measured mass defect of the choice quadrature, at h1's policies.
    defect = 0.0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        cm = choice_maps(h1.WT[g, l, zi, :], h1.WO[g, l, zi, :], g, gr)
        iT, iO = integrate_choice(ones(gr.nϵT), ones(gr.nXO), cm, gr)
        defect = max(defect, abs(iT + iO - 1))
    end
    shift = sol.p.λ * abs(log(href))
    tolΛ = max(5 * shift * defect, 1e-12)
    tolp = 1e-5                       # policy response is second order in the defect
    dev(a, b) = maximum(abs, a .- b)
    dΛ = maximum(abs, (h2.Λ .- h1.Λ) .+ sol.p.λ * log(href))
    dp = max(dev(h1.eT, h2.eT), dev(h1.eO, h2.eO), dev(h1.hT, h2.hT),
             dev(h1.sT, h2.sT), dev(h1.sO, h2.sO))
    @printf("  quadrature mass defect               : %.2e  (⇒ Λ tolerance %.1e)\n", defect, tolΛ)
    @printf("  policies (e, h, s) unchanged by href : %.2e   [%s]\n", dp, mark(dp < tolp))
    @printf("  Λ shifted by −λ log href             : %.2e   [%s]\n", dΛ, mark(dΛ < tolΛ))
    println("==================================================")
    return dp < tolp && dΛ < tolΛ
end

"""
    prop2prime_test(sol; nsample = 6, atol = 1e-5)

Independent check of Proposition 2′ (the explicit s policy). At a sample of
household states, maximise the true objective

    F(s, e) = log(1 − s) + V̄( W(e, s) )

JOINTLY in (s, e) by 2-D numerical search, and compare the optimum to the stored
(`hh.sT`/`hh.eT`, `hh.sO`/`hh.eO`) policy. This is the one place the s analytics
could be wrong in a way no other test would catch: every other diagnostic takes
the FOC as given.

The search is deliberately started AWAY from the stored policy (`pert`): started
at it, a derivative-free method would report the start point back and the test
would be vacuous.
"""
function prop2prime_test(sol; nsample = 6, atol = 1e-5, pert = 0.25)
    (; hh, gr, p) = sol
    println("\n========== Proposition 2′ test (joint (s,e) optimum) ==========")
    sig(x) = 1 / (1 + exp(-x))
    logit(s) = log(s / (1 - s))
    opts = Optim.Options(g_tol = 1e-14, x_abstol = 1e-14, f_abstol = 1e-16,
                         iterations = 100_000)
    worst_s = 0.0; worst_e = 0.0; ncheck = 0
    # Objective: the true (s, e) problem, log(1−s) plus the location log-sum.
    obj(Wof, l) = v -> -(log(1 - sig(v[1])) +
                         logsum_probs(Wof(exp(v[2]), sig(v[1])), l, p)[1])
    states = [(g, l, zi) for g in 1:2 for l in 1:gr.L for zi in 1:gr.Nz]
    for (g, l, zi) in states[1:min(nsample, length(states))]
        for k in (1, max(1, gr.nϵT ÷ 2), gr.nϵT)
            s★, e★ = hh.sT[g, l, zi, k], hh.eT[g, l, zi, k]
            WofT = (e, s) -> W_teach(gr.ϵTgrid[k], zi, l, g, e, s, hh.Λ, p, gr)[1]
            v = Optim.minimizer(optimize(obj(WofT, l),
                                         [logit(s★ * (1 - pert)), log(e★ * (1 + pert))],
                                         NelderMead(), opts))
            worst_s = max(worst_s, abs(sig(v[1]) - s★))
            worst_e = max(worst_e, abs(exp(v[2]) - e★) / e★)
            ncheck += 1
        end
        for k in (1, max(1, gr.nXO ÷ 2), gr.nXO)
            s★, e★ = hh.sO[g, l, zi, k], hh.eO[g, l, zi, k]
            WofO = (e, s) -> W_nonteach(gr.XOgrid[k, g], zi, l, g, e, s, hh.Λ, p, gr)[1]
            v = Optim.minimizer(optimize(obj(WofO, l),
                                         [logit(s★ * (1 - pert)), log(e★ * (1 + pert))],
                                         NelderMead(), opts))
            worst_s = max(worst_s, abs(sig(v[1]) - s★))
            worst_e = max(worst_e, abs(exp(v[2]) - e★) / e★)
            ncheck += 1
        end
    end
    ok = worst_s < atol && worst_e < 1e-3
    @printf("  states checked (started %.0f%% off the policy) : %d\n", 100 * pert, ncheck)
    @printf("  worst |Δs| vs 2-D optimum      : %.2e   [%s]\n", worst_s, mark(worst_s < atol))
    @printf("  worst rel |Δe| vs 2-D optimum  : %.2e   [%s]\n", worst_e, mark(worst_e < 1e-3))
    println("==============================================================")
    return ok
end

"""
    psi_sweep(; href, ψs, kw)

Sweep the warm-glow curvature ψ as a CONTINUATION, reporting the sorting
diagnostics at each point. Everything else is held at the BASELINE, so the goods
moving cost is on throughout and this sweep traces the warm-glow mechanism's
marginal contribution on top of it. ψ = 1 is the warm-glow-off ablation;
`Params().ψ = 0` is the baseline and the far end of the sweep.

`href` must be frozen at the benchmark mean child human capital (see
`mean_child_h`) or ψ silently rescales the whole altruism term instead of only
bending it — hence the default, the baseline's own frozen `HBAR_BENCH`.

Expected: g rises monotonically as ψ falls; `idx_work_signed` crosses zero between
ψ = 1 (≈ −0.011) and ψ = 0.5 (≈ +0.035), reaching ≈ +0.077 at the baseline ψ = 0;
and the ψ > 1 point pushes it further negative — the cleanest falsification of the
mechanism.
"""
function psi_sweep(; href = Params().href, ψs = [1.5, 1.0, 0.75, 0.5, 0.25, 0.0],
                     kw = _mech_kw(), outdir = nothing)
    println("\n========== warm-glow curvature ψ sweep (baseline ψ = $(Params().ψ)) ==========")
    @printf("  href (frozen h̄) = %.6f\n", href)
    @printf("  %-6s %-11s %-11s %-13s %-11s %-11s %-9s\n",
            "ψ", "g_T(loc1)", "g_O(loc1)", "idx_work±", "within", "composition", "PropA")
    rows = NamedTuple[]
    prev = nothing
    for ψ in ψs
        NONMONO_WARNINGS[] = 0
        sol = solve_ge(Params(ψ = ψ, href = href); kw..., init = prev)
        prev = sol
        wg = within_branch_gradient(sol); sm = sorting_metrics(sol)
        dc = sorting_decomposition(sol)
        @printf("  %-6.2f %+-11.5f %+-11.5f %+-13.5f %+-11.6f %+-11.6f %-9.1e\n",
                ψ, wg.gT[1], wg.gO[1], sm.idx_work_signed, dc.within_tot, dc.comp_tot,
                ability_neutrality_dev(sol))
        push!(rows, (; ψ, gT = wg.gT[1], gO = wg.gO[1], idx = sm.idx_work_signed,
                       within = dc.within_tot, comp = dc.comp_tot,
                       nonmono = NONMONO_WARNINGS[]))
    end
    nm = sum(r.nonmono for r in rows)
    @printf("  NONMONO_WARNINGS across the sweep : %d   [%s]\n", nm, mark(nm == 0))
    println("  ψ works through `zi`, not the shock index k, so it does not move Prop. A:")
    println("  the level here (≈0.4) is set by the baseline's goods moving cost and is")
    println("  flat across the sweep. It is ≈0 only when mcost is also switched off.")
    println("==========================================================")
    if outdir !== nothing
        mkpath(outdir)
        ψv = [r.ψ for r in rows]
        p1 = plot(ψv, [[r.gT for r in rows] [r.gO for r in rows]], marker = :circle, lw = 2,
                  xlabel = "ψ", ylabel = "within-branch gradient g",
                  label = ["g_T" "g_O"], title = "Ability gradient vs warm-glow curvature")
        hline!(p1, [0.0], ls = :dash, color = :black, label = "")
        p2 = plot(ψv, [[r.idx for r in rows] [r.within for r in rows] [r.comp for r in rows]],
                  marker = :circle, lw = 2, xlabel = "ψ",
                  label = ["idx_work (signed)" "within-branch" "composition"],
                  title = "Sorting and its decomposition")
        hline!(p2, [0.0], ls = :dash, color = :black, label = "")
        plot(p1, p2, layout = (1, 2), size = (1100, 400))
        savefig(joinpath(outdir, "continuous_psi_sweep.png"))
    end
    return rows
end

"""
    mcost_test(; Cbar = Params().Cbar, kw)

The moving-cost mechanism in isolation: hold the warm glow at the baseline ψ and
compare the cost denominated in UTILITY (`τmove = TAUMOVE_BENCH`, `mcost = 0`)
against the same cost denominated in GOODS — which is what the baseline ships.

Two things are checked. First a CALIBRATION assertion: the baseline's `mcost` must
be exactly the re-denomination m = (1 − e^{−τmove/μ})·C̄ of the old utility cost at
the frozen `Cbar`, so the switch really did cost no new free parameter and
`MCOST_BENCH` has not drifted from the formula that defines it. Then the
behavioural read-out: the s policy spreads by a few points, the s inner loop stays
cheap, Prop. A breaks (that is the mechanism), the gradient appears in BOTH
branches, and the feasibility audit still passes.
"""
function mcost_test(; Cbar = Params().Cbar, kw = _mech_kw())
    println("\n========== goods- vs utility-denominated moving cost ==========")
    τb   = TAUMOVE_BENCH
    base = Params()                                             # goods cost — the baseline
    putil = Params(τmove = [0.0 τb; τb 0.0], mcost = zeros(2, 2))  # same cost, utility units

    # Calibration assertion: the shipped mcost IS the re-denominated τmove.
    m_expect = redenominate_move_cost(putil, Cbar).mcost
    ok_cal = isapprox(m_expect, base.mcost; atol = 1e-12)
    @printf("  C̄ (frozen) = %.6f   τmove_off = %.3f  ⇒  m = %.6f   (m/C̄ = %.3f)\n",
            Cbar, τb, base.mcost[1, 2], base.mcost[1, 2] / Cbar)
    @printf("  baseline mcost == redenominate(τmove, C̄)  : %s   (max dev %.2e)\n",
            mark(ok_cal), maximum(abs, m_expect .- base.mcost))

    sols = (("utility cost τmove (mechanism off)", solve_ge(putil; kw...)),
            ("goods cost mcost  (BASELINE)",       solve_ge(base;  kw...)))
    for (lab, s) in sols
        wg = within_branch_gradient(s); sm = sorting_metrics(s); dc = sorting_decomposition(s)
        ss = s_summary(s)
        println("\n  ", lab, ":")
        @printf("    g_T = %+.5f   g_O = %+.5f   idx_work± = %+.5f\n",
                wg.gT[1], wg.gO[1], sm.idx_work_signed)
        @printf("    decomposition: within %+.6f | composition %+.6f | interaction %+.6f\n",
                dc.within_tot, dc.comp_tot, dc.inter_tot)
        @printf("    Prop. A deviation = %.2e\n", ability_neutrality_dev(s))
        @printf("    s_T ∈ [%.4f, %.4f]   s_O ∈ [%.4f, %.4f]\n", ss.rangeT..., ss.rangeO...)
        verify_solution(s)
    end
    println("  (s is constant across states in the first block and spreads in the")
    println("   second: that spread IS the Proposition-2′ Ξ correction switching on.)")
    println("====================================================================")
    return (; util = sols[1][2], base = sols[2][2], ok_cal)
end

"""
    sorting_mechanism_tests(; outdir)

The sorting battery end to end:

  0. `nesting_params()` — the pre-sorting regression point. Read out its sorting
     diagnostics (they should be the "no sorting" case: g ≈ 0, `idx_work_signed`
     negative, ~100 % compositional) and RE-DERIVE the two frozen normalizations
     h̄ and C̄ from it, asserting they still equal the `HBAR_BENCH` / `CBAR_BENCH`
     the baseline ships. That assertion is what keeps `Params()`'s hard-coded
     constants honest: if the model changes underneath them, this fails rather
     than silently shifting the baseline's calibration.
  1. ψ sweep — the warm-glow mechanism, on top of the goods moving cost.
  2. `mcost_test` — the moving-cost mechanism, on top of the warm glow, plus the
     calibration assertion that `Params().mcost` is the re-denominated `τmove`.
  3. the BASELINE itself, with both mechanisms on — the headline read-out.

Proposition 2′ is checked twice: at `nesting_params()`, where `s` is constant and
the test is a regression, and at the baseline, where `s` is genuinely
state-dependent and the test is the real one.
"""
function sorting_mechanism_tests(; outdir = FIGDIR, kw = _mech_kw())
    println("\n\n############ SORTING MECHANISMS ############")

    # ---- The pre-sorting regression point + the frozen normalizations.
    println("\n---- nesting_params() read-out (mechanisms OFF) ----")
    nest = solve_ge(nesting_params(); kw...)
    ok_href = href_invariance_test(nest)
    report_sorting(nest)
    ok_p2 = prop2prime_test(nest)

    href = mean_child_h(nest)
    Cbar = mean_consumption(nest)
    dh, dc = abs(href - HBAR_BENCH), abs(Cbar - CBAR_BENCH)
    # Loose tolerance: these are re-derived at the caller's grids/tolerances, not
    # at the Nz=5, 48² solve the constants were minted from, so a few 1e-4 of
    # quadrature drift is expected. A shift big enough to matter for the
    # calibration would be orders of magnitude larger.
    atol_norm = 1e-3
    ok_norm = dh < atol_norm && dc < atol_norm
    println("\n  Frozen normalizations re-derived at nesting_params():")
    @printf("    h̄ = %.6f  vs HBAR_BENCH = %.6f   (dev %.2e)  [%s]\n",
            href, HBAR_BENCH, dh, mark(dh < atol_norm))
    @printf("    C̄ = %.6f  vs CBAR_BENCH = %.6f   (dev %.2e)  [%s]\n",
            Cbar, CBAR_BENCH, dc, mark(dc < atol_norm))
    ok_norm || println("    ⚠ the baseline's frozen constants no longer match the model — " *
                       "re-mint HBAR_BENCH/CBAR_BENCH in spatial_continuous.jl.")

    # ---- Each mechanism's marginal contribution.
    rows = psi_sweep(; href = Params().href, kw, outdir)
    mt   = mcost_test(; Cbar = Params().Cbar, kw)

    # ---- The baseline, both mechanisms on. Gradients additive to first order.
    # `mcost_test` already solved `Params()` at these settings — reuse it rather
    # than paying for an identical solve.
    println("\n---- the BASELINE (ψ = $(Params().ψ), goods mcost) ----")
    base = mt.base
    report_sorting(base)
    ok_p2b = prop2prime_test(base)            # s is state-dependent here — the real test
    verify_solution(base)

    sm_n, sm_b = sorting_metrics(nest), sorting_metrics(base)
    @printf("\n  idx_work_signed:  nesting %+.5f  →  baseline %+.5f   [%s]\n",
            sm_n.idx_work_signed, sm_b.idx_work_signed,
            mark(sm_n.idx_work_signed < 0 < sm_b.idx_work_signed))
    println("############################################\n")
    return (; nest, base, rows, mt,
              ok = ok_href && ok_p2 && ok_p2b && ok_norm && mt.ok_cal &&
                   sm_n.idx_work_signed < 0 < sm_b.idx_work_signed)
end

# ---- How occupation choice and location choice interact. --------------------
# Teaching is decided when young; location when working. The two margins couple:
# (i) teaching selection differs across BIRTH locations (different G_l(z), Q_l);
# (ii) teachers and non-teachers MIGRATE differently (teaching wages run through
# κ_l' h^γ, so destinations matter differently for the two occupations).
function analyze_occupation_location(sol; outdir = FIGDIR)
    mkpath(outdir)
    (; L) = sol.gr
    tm = teacher_moments(sol)
    T, O = tm.Tmass, tm.Omass
    println("\n========== Occupation × location interaction ==========")
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
    println("=======================================================")

    p1 = bar(["born $l0" for l0 in 1:L], teach_born, legend = false,
             ylabel = "teaching share", title = "Teaching share by birth location",
             ylims = (0, max(0.6, maximum(teach_born) * 1.15)))
    p2 = grouped_bar(["born $l0" for l0 in 1:L], ["teacher", "non-teacher"],
                     hcat(stay_T, stay_O); ylabel = "P(work = birth loc)",
                     title = "Stay probability by occupation", ylims = (0, 1), legend = :topright)
    plot(p1, p2, size = (1000, 380))
    savefig(joinpath(outdir, "continuous_occupation_location.png"))
    return (; teach_born, stay_T, stay_O)
end

# ---- Class size, teacher-HC distribution, quantity vs quality. --------------
function analyze_class_and_teachers(sol; outdir = FIGDIR)
    mkpath(outdir)
    (; L) = sol.gr
    tm = teacher_moments(sol)
    println("\n========== Class size & teacher human capital ==========")
    @printf("  %-6s %-9s %-9s %-11s %-11s %-9s %-9s\n",
            "loc", "Q_l", "class", "#teachers", "mean logh", "sd logh", "H̃_T")
    for l in 1:L
        @printf("  %-6d %-9.4f %-9.3f %-11.4f %-11.4f %-9.4f %-9.4f\n",
                l, tm.Q[l], tm.classsize[l], tm.Tcount[l], tm.meanlogh[l],
                tm.sdlogh[l], tm.HT[l])
    end
    println("\n  Decomposition  H̃_T,l = (#teachers) × (mean h^{β/σ}):")
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
    println("========================================================")

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
    savefig(joinpath(outdir, "continuous_class_and_teacher_moments.png"))

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
    savefig(pdens, joinpath(outdir, "continuous_teacher_hc_distribution.png"))
    println("  Saved class-size / teacher-distribution figures to ", outdir)
    return tm
end

# ---- What drives spatial sorting on ability? --------------------------------
# Shut down each source of location asymmetry one at a time and measure the
# residual sorting (spread in post-migration mean z across work locations). The
# channel whose removal collapses the sorting index is the one that drives it.
function analyze_sorting_drivers(; outdir = FIGDIR, kw = _sweep_kw())
    mkpath(outdir)
    println("\n========== Drivers of spatial sorting on ability ==========")
    cases = [
        (label = "baseline",         p = Params()),
        (label = "no amenity (B=0)", p = Params(B = [0.0, 0.0])),
        (label = "equal κ",          p = Params(κ = [1.0, 1.0])),
        (label = "no altruism (λ=0)",p = Params(λ = 0.0)),
        (label = "no move cost",     p = Params(mcost = zeros(2, 2))),
        (label = "log warm glow ψ=1",p = Params(ψ = 1.0, href = 1.0)),
        (label = "nesting (both off)",p = nesting_params()),
        (label = "symmetric (B=0,κ=)",p = Params(B = [0.0, 0.0], κ = [1.0, 1.0])),
    ]
    # The SIGNED index is what is plotted and returned: the sign is the result,
    # and the |max−min| spread would hide it.
    @printf("  %-20s %-12s %-10s %-10s %-12s %-10s\n",
            "case", "idx_work±", "idx_work", "idx_birth", "share_loc2", "outflow")
    idxs = Float64[]; labels = String[]
    for c in cases
        sol = solve_ge(c.p; kw...)
        sm  = sorting_metrics(sol)
        @printf("  %-20s %+-12.5f %-10.4f %-10.4f %-12.4f %-10.4f\n",
                c.label, sm.idx_work_signed, sm.idx_work, sm.idx_birth,
                sol.M[2] / sum(sol.M), outflow_rate(sol))
        push!(idxs, sm.idx_work_signed); push!(labels, c.label)
    end
    println("\n  The sorting index (work-location mean-z gap, loc L − loc 1) collapses")
    println("  toward the fully-symmetric benchmark as the ACTIVE asymmetry is removed;")
    println("  whichever single-channel removal cuts it most is the dominant driver.")
    println("  The baseline value is POSITIVE — the attractive location draws HIGH-z")
    println("  agents. The 'no move cost' / 'log warm glow' / 'nesting' rows are the")
    println("  mechanism ablation: dropping either sorting mechanism shrinks the index,")
    println("  and dropping both ('nesting') drives it NEGATIVE, back to the pre-sorting")
    println("  regime where the attractive location drew LOW-z agents through")
    println("  occupational composition alone. Note κ remains the dominant LOCATION")
    println("  asymmetry: 'equal κ' collapses the index further than either ablation.")
    println("===========================================================")
    bar(labels, idxs, legend = false, xrotation = 25,
        ylabel = "signed sorting index (mean-z gap, loc L − loc 1)",
        title = "What drives spatial sorting on ability?", size = (900, 450),
        bottom_margin = 8Plots.mm)
    savefig(joinpath(outdir, "continuous_sorting_drivers.png"))
    return idxs
end

# ---- Amenities vs altruism in location choice; who moves? -------------------
function analyze_amenity_altruism(; outdir = FIGDIR, kw = _sweep_kw())
    mkpath(outdir)
    println("\n========== Amenities vs altruism; mobility by ability ==========")
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
    println("================================================================")
    pl = plot(xlabel = "ability z", ylabel = "P(leave birth location)",
              title = "Who moves? Outflow rate by ability", legend = :topleft)
    for r in results
        plot!(pl, r.z, r.outflow, marker = :circle, lw = 2, label = r.label)
    end
    savefig(pl, joinpath(outdir, "continuous_mobility_by_ability.png"))
    println("  Saved mobility-by-ability figure to ", outdir)
    return results
end

# ---- The β sweep and the β/σ corner. ----------------------------------------
# The teacher block self-consistency: H̃_T = ½∫h^{β/σ}, Q = (2H̃_T/M)^σ, and h ∝ Q.
# Substituting Q^{β/σ} = (2H̃_T/M)^β gives H̃_T^{1-β} ∝ M^{-β}·(policy integral): the
# interior teacher stock is well-defined only for β<1, and as β→1 the map is
# increasingly ill-conditioned. Independently, expo ≡ β/σ sets the CONVEXITY of the
# aggregator: larger β/σ concentrates H̃_T on the best teachers and strengthens the
# Q→h→H̃_T feedback, which is what tips high-β economies into the
# near-universal-teaching / collapsing-quality corner.
function beta_sigma_analysis(; outdir = FIGDIR, kw = _sweep_kw())
    mkpath(outdir)
    println("\n========== β sweep and the β/σ corner ==========")

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
    println("================================================")

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
    savefig(joinpath(outdir, "continuous_beta_sigma_collapse.png"))
    println("  Saved β/σ figures to ", outdir)
    return (; βcurve, tscurve, htcurve, curves)
end

"Run every mechanism analysis on/around a converged baseline."
function mechanism_analyses(sol; outdir = FIGDIR, kw = _sweep_kw())
    println("\n\n#################### MECHANISM ANALYSES ####################")
    analyze_occupation_location(sol; outdir)
    analyze_class_and_teachers(sol; outdir)
    analyze_sorting_drivers(; outdir, kw = _sweep_kw(Nz = 7))   # higher Nz for the ability gradient
    analyze_amenity_altruism(; outdir, kw = _sweep_kw(Nz = 7))
    quantile_sensitivity(; outdir)
    beta_sigma_analysis(; outdir, kw)
    println("############################################################\n")
end

# -----------------------------------------------------------------------------
# 7. Figures.
# -----------------------------------------------------------------------------

"Aggregate GE objects: the ability distributions by birth and by work location,
 the raw policy/value functions, and the migration matrix."
function plot_ge_outcomes(sol; outdir = FIGDIR)
    mkpath(outdir)
    (; z, ϵTgrid, XOgrid, Nz, L, Πz) = sol.gr
    hh = sol.hh
    g = 1
    zmed = (Nz + 1) ÷ 2

    erg = stationary(Πz)
    pd = plot(xlabel = "z", ylabel = "mass", title = "Endogenous ability distribution")
    for l in 1:L
        gl = sol.Φ[:, l] ./ sum(@view sol.Φ[:, l])
        plot!(pd, z, gl, marker = mk(l), label = "G_$(l)(z)")
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
        plot!(pw, z, gwl, marker = mk(l), label = "G^work_$(l)(z)")
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
                 title = "Migration matrix Π̄[born→work]")
    annotate!(pm, vec([(j, i, text(@sprintf("%.3f", Π̄[i, j]), 9, :white))
                       for i in 1:L, j in 1:L]))
    savefig(pm, joinpath(outdir, "continuous_migration_matrix.png"))

    println("\nSaved GE figures to ", outdir)
end

"""
    plot_student_ability(sol; outdir, na = 400, trim = 0.005)

The ability distribution of STUDENTS (the young) in each location, in both
ability concepts, with the occupational choice folded in.

Students are educated where they are born, so the relevant location index is the
BIRTH location — `Φ[z, l]` — not the work location the `_by_work` figure uses.
Four panels:

  (1) density of a = z·ε in the CHOSEN occupation, one curve per birth location
      (`student_ability_density`). This is the distribution the occupation margin
      actually operates on, and it sits to the right of z: both branches select
      on their own ε, so E[a] > E[z] even though ε is unit-mean unconditionally.
  (2) the same density split into teachers and non-teachers, each carrying its
      branch mass so the two curves ADD UP to panel (1)'s — the visual form of
      "taking the occupational choice into account".
  (3) the persistent component: G_l(z) overall (solid), among teachers (dashed)
      and among non-teachers (dotted), against the ergodic G*(z). z is discrete
      (Rouwenhorst), so these are mass points, not a density.
  (4) mean ability a by location, overall and by branch — the numbers
      `student_ability_moments` returns and `report_ge` prints.

`a` is evaluated on `na` log-spaced points and the axis is trimmed to the central
1 − 2·`trim` of the pooled distribution (the ability support runs out to z_max·ε
at the far tail of ε, where there is no visible mass).
"""
function plot_student_ability(sol; outdir = FIGDIR, na = 400, trim = 0.005)
    mkpath(outdir)
    (; z, L, Πz, dT) = sol.gr
    sa  = student_ability_moments(sol)
    erg = stationary(Πz)
    lcolors = palette(:tab10, max(L, 3))

    # One common a-grid for every location, spanning the full (z, ε) support.
    agrid = exp.(range(log(z[1] * quantile(dT, 1e-4)),
                       log(z[end] * quantile(dT, 1 - 1e-4)), na))
    dens  = [student_ability_density(sol, l, agrid) for l in 1:L]
    # Trim the axis to the bulk, pooling locations by their population mass.
    wl    = sa.mass ./ sum(sa.mass)
    fpool = sum(wl[l] .* (dens[l][1] .+ dens[l][2]) for l in 1:L)
    Fpool = cumtrapz(agrid, fpool)
    Fpool ./= Fpool[end]
    alo = agrid[max(searchsortedfirst(Fpool, trim), 1)]
    ahi = agrid[min(searchsortedfirst(Fpool, 1 - trim), na)]

    pa = plot(xlabel = "ability a = z·ε (chosen occupation)", ylabel = "density",
              title = "Student ability by birth location", legend = :topright,
              xlims = (alo, ahi))
    pb = plot(xlabel = "ability a = z·ε (chosen occupation)", ylabel = "density",
              title = "… split by occupation (areas = branch shares)",
              legend = :topright, xlims = (alo, ahi))
    for l in 1:L
        fT, fO = dens[l]
        c = lcolors[mod1(l, L)]
        plot!(pa, agrid, fT .+ fO, lw = 2, color = c,
              label = @sprintf("loc %d  (mean a = %.3f)", l, sa.mean_a[l]))
        plot!(pb, agrid, fT, lw = 2, color = c, ls = :dash,
              label = @sprintf("loc %d teach (%.0f%%)", l, 100 * sa.teach_share[l]))
        plot!(pb, agrid, fO, lw = 2, color = c, ls = :dot,
              label = @sprintf("loc %d nonteach (%.0f%%)", l, 100 - 100 * sa.teach_share[l]))
    end

    pz = plot(xlabel = "persistent ability z", ylabel = "mass",
              title = "G_l(z) overall and by occupation", legend = :topright)
    for l in 1:L
        c = lcolors[mod1(l, L)]
        plot!(pz, z, sol.Φ[:, l] ./ sum(@view sol.Φ[:, l]), lw = 2, color = c,
              marker = mk(l), label = "loc $l")
        sum(@view sa.ΦT[:, l]) > 0 && plot!(pz, z, sa.ΦT[:, l] ./ sum(@view sa.ΦT[:, l]),
              lw = 2, color = c, ls = :dash, label = "loc $l | teach")
        sum(@view sa.ΦO[:, l]) > 0 && plot!(pz, z, sa.ΦO[:, l] ./ sum(@view sa.ΦO[:, l]),
              lw = 2, color = c, ls = :dot, label = "loc $l | nonteach")
    end
    plot!(pz, z, erg, lw = 2, color = :gray, ls = :dashdot, label = "ergodic G*(z)")

    # Means as a DOT plot, not bars: the cross-location differences are a few
    # percent, so a zero-baseline bar chart hides exactly what is being compared.
    # mean z is on the same axis on purpose — the gap a − z is the SELECTION
    # premium the occupational choice creates out of a unit-mean ε.
    pm = plot(xlabel = "birth location", ylabel = "mean ability",
              title = @sprintf("Mean ability by location (pooled a = %.3f)", sa.mean_a_all),
              legend = :outertopright, xticks = (1:L, ["loc $l" for l in 1:L]),
              xlims = (0.5, L + 0.5))
    for (lab, v, mkr, c) in (("mean z",          sa.mean_z,       :diamond,   :gray),
                             ("mean a (all)",    sa.mean_a,       :circle,    1),
                             ("a | teachers",    sa.mean_a_teach, :utriangle, 2),
                             ("a | non-teach",   sa.mean_a_non,   :square,    3))
        plot!(pm, 1:L, v, lw = 2, marker = mkr, ms = 7, color = c, label = lab)
    end

    plot(pa, pb, pz, pm, layout = (2, 2), size = (1200, 800),
         left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
    savefig(joinpath(outdir, "continuous_student_ability.png"))
    println("Saved student-ability figures to ", outdir)
    return sa
end

"""
    plot_household(sol; outdir)

Household-choice intuition, as opposed to the raw policy objects in
`plot_ge_outcomes`:
  (1) the occupation margin — teach iff W_T(ε_T) > W_O(X_O*) — as a smooth
      teaching PROBABILITY in each shock, and the value crossing that generates it;
  (2) how teaching SELECTS on ability z, and how location shifts that selection;
  (3) the shock DISTRIBUTIONS agents draw from, overlaid with who ends up teaching,
      plus the stationary joint mass over (ability, location).
"""
function plot_household(sol; outdir = FIGDIR)
    mkpath(outdir)
    grd = sol.gr
    (; z, ϵTgrid, XOgrid, Nz, L, dT) = grd
    hh = sol.hh
    p  = sol.p
    g  = 1                                     # gender slice for the single-state panels
    zmed = (Nz + 1) ÷ 2
    lmed = (L + 1) ÷ 2
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
        plot!(p3, z, ts, lw = 2, marker = mk(l), label = "loc $l")
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
    p5 = heatmap(["loc $l" for l in 1:L],
                 [@sprintf("z=%.2f", z[zi]) for zi in 1:Nz],
                 sol.Φ, c = :viridis, title = "Stationary mass Φ[z, location]")
    savefig(p5, joinpath(outdir, "continuous_joint_distribution.png"))

    println("Saved household-choice figures to ", outdir)
end

"""
    plot_location_choice(sol; outdir)

WHERE agents go, by ability and by occupational branch. One column per birth
location:

  top row     π_{l'|l₀}(z) against the persistent ability z, teachers (solid)
              versus non-teachers (dashed). The slope of each curve is the
              within-branch ability gradient `within_branch_gradient` summarises;
              the gap between the branches is the occupational-composition channel
              that `sorting_decomposition` separates out.
  bottom row  the same probabilities against realised ability a = z·ε at the
              median z node, so the shock margin is visible too: teachers are
              indexed by a = z·ε_T, non-teachers by a = z·E[ε_{i*} | X_O*], which
              puts both branches on one ability axis. Axes are trimmed to the
              central 96 % of each shock distribution.

Location choice is ability-neutral within a branch whenever `mcost = 0` (Prop. A),
so at `nesting_params()` every curve here is FLAT and only the teach/non-teach gap
survives; the baseline's goods moving cost is what tilts them.
"""
function plot_location_choice(sol; outdir = FIGDIR)
    mkpath(outdir)
    grd = sol.gr
    (; z, ϵTgrid, XOgrid, Nz, L, dT) = grd
    p = sol.p
    g = 1                                       # gender slice for the shock-axis row
    zmed = (Nz + 1) ÷ 2
    dcolors = palette(:tab10, max(L, 3))

    # Headroom above 1 so the legend sits clear of the curves in the first panel.
    ylim = (-0.05, 1.32)
    leg(l0) = l0 == 1 ? :top : false

    panels = Any[]
    # --- top row: probability of working in l', by persistent ability z.
    for l0 in 1:L
        lc = location_choice_by_ability(sol, l0)
        pl = plot(xlabel = "ability z", ylabel = l0 == 1 ? "P(work in l′)" : "",
                  title = "Born in loc $l0", ylims = ylim,
                  legend = leg(l0), legend_columns = 2, legendfontsize = 7)
        for lp in 1:L
            c = dcolors[mod1(lp, L)]
            plot!(pl, z, lc.πT[:, lp], lw = 2, color = c, marker = mk(lp), ms = 4,
                  label = "→ $lp | teach")
            plot!(pl, z, lc.πO[:, lp], lw = 2, color = c, ls = :dash, marker = mk(lp),
                  ms = 4, markerstrokecolor = c, label = "→ $lp | nonteach")
        end
        push!(panels, pl)
    end

    # --- bottom row: the same, against realised ability a = z·ε at the median z.
    iT = bulk_range(ϵTgrid, ε -> cdf(dT, ε))
    iO = bulk_range(view(XOgrid, :, g), x -> Fxo(grd, g, x))
    aT = z[zmed] .* ϵTgrid[iT]
    aO = [z[zmed] * Eϵ_nonteach(grd, g, XOgrid[k, g], p.α) for k in iO]
    for l0 in 1:L
        pl = plot(xlabel = @sprintf("ability a = z·ε  (z = %.2f)", z[zmed]),
                  ylabel = l0 == 1 ? "P(work in l′)" : "",
                  title = "Born in loc $l0 — by realised ability",
                  ylims = ylim, xlims = (min(aT[1], aO[1]), max(aT[end], aO[end])),
                  legend = leg(l0), legend_columns = 2, legendfontsize = 7)
        for lp in 1:L
            c = dcolors[mod1(lp, L)]
            plot!(pl, aT, sol.Pi_T[g, l0, zmed, iT, lp], lw = 2, color = c,
                  label = "→ $lp | teach")
            plot!(pl, aO, sol.Pi_O[g, l0, zmed, iO, lp], lw = 2, color = c, ls = :dash,
                  label = "→ $lp | nonteach")
        end
        push!(panels, pl)
    end

    plot(panels..., layout = (2, L), size = (520 * L, 760),
         left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
    savefig(joinpath(outdir, "continuous_location_choice.png"))
    println("Saved location-choice figures to ", outdir)
    return nothing
end

"""
    plot_e_ratio(sol; outdir, lnum = 2, lden = 1)

The goods-investment policy in location `lnum` RELATIVE to location `lden`, for
teachers and for non-teachers separately: e(a; lnum) / e(a; lden), one curve per
ability node z, plotted against realised ability a = z·ε.

The ratio, rather than the two levels, is the object of interest because the
locations differ only through what a parent faces there — the teacher-quality
index Q_l, the teaching-wage shifter κ_l, the tax rate t_l, the amenity B_l and
the moving costs out of l. A ratio above 1 says parents in `lnum` buy MORE goods
investment for the same child; a ratio that varies with a says the location gap is
not a pure level shift but changes the return to investing in a better child.

Teachers and non-teachers are shown separately because κ_l enters only the
teaching branch: any gap between the two panels is the part of the location effect
that runs through teaching income rather than through Q_l and t_l, which both
branches share.
"""
function plot_e_ratio(sol; outdir = FIGDIR, lnum = 2, lden = 1)
    mkpath(outdir)
    grd = sol.gr
    (; z, ϵTgrid, XOgrid, Nz, L, dT) = grd
    p  = sol.p
    hh = sol.hh
    g  = 1
    L < 2 && return nothing
    zcolors = palette(:viridis, Nz)

    iT = bulk_range(ϵTgrid, ε -> cdf(dT, ε))
    iO = bulk_range(view(XOgrid, :, g), x -> Fxo(grd, g, x))
    εO = [Eϵ_nonteach(grd, g, XOgrid[k, g], p.α) for k in iO]

    pT = plot(xlabel = "ability a = z·ε_T", ylabel = @sprintf("e(loc %d) / e(loc %d)", lnum, lden),
              title = "Teachers", legend = :topright)
    pO = plot(xlabel = "ability a = z·ε_{i*}", ylabel = "",
              title = "Non-teachers", legend = false)
    rmin = Inf; rmax = -Inf
    for zi in 1:Nz
        rT = hh.eT[g, lnum, zi, iT] ./ hh.eT[g, lden, zi, iT]
        rO = hh.eO[g, lnum, zi, iO] ./ hh.eO[g, lden, zi, iO]
        rmin = min(rmin, minimum(rT), minimum(rO))
        rmax = max(rmax, maximum(rT), maximum(rO))
        plot!(pT, z[zi] .* ϵTgrid[iT], rT, lw = 2, color = zcolors[zi],
              label = @sprintf("z=%.2f", z[zi]))
        plot!(pO, z[zi] .* εO, rO, lw = 2, color = zcolors[zi], label = "")
    end
    hline!(pT, [1.0], ls = :dash, color = :black, label = "")
    hline!(pO, [1.0], ls = :dash, color = :black, label = "")
    # Shared y-axis (including the 1.0 reference) so the two branches are
    # directly comparable — the teacher/non-teacher gap is the κ_l channel.
    ylim = (min(0.99, 0.99 * rmin), max(1.01, 1.01 * rmax))
    plot!(pT; ylims = ylim); plot!(pO; ylims = ylim)

    plot(pT, pO, layout = (1, 2), size = (1100, 430),
         plot_title = @sprintf("Goods investment e: location %d relative to location %d",
                               lnum, lden),
         plot_titlevspan = 0.12,
         left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 3Plots.mm)
    savefig(joinpath(outdir, "continuous_e_ratio_by_location.png"))

    @printf("\n  e(loc %d)/e(loc %d) over the plotted range: [%.4f, %.4f]   (Q = %s, t = %s, κ = %s)\n",
            lnum, lden, rmin, rmax, fmtvec(grd.Q), fmtvec(sol.t), fmtvec(p.κ))
    println("Saved goods-investment ratio figure to ", outdir)
    return (; rmin, rmax)
end

"Every figure that reads a single converged solution."
function plot_all(sol; outdir = FIGDIR)
    plot_ge_outcomes(sol; outdir)
    plot_student_ability(sol; outdir)
    plot_household(sol; outdir)
    plot_location_choice(sol; outdir)
    plot_e_ratio(sol; outdir)
    return nothing
end

# -----------------------------------------------------------------------------
# 8. Run.
# -----------------------------------------------------------------------------
function main_test()
    println("Solving the spatial model in general equilibrium (baseline) ...")
    sol = solve_ge(Params(); Nz = 5, nϵT = 48, nXO = 48, damping = 0.5, tol = 1e-6,
                   maxit = 250, hh_tol = 1e-6, hh_maxit = 1000, verbose = true)
    report_ge(sol)
    verify_solution(sol)
    report_checks(sol)
    report_sorting(sol)
    plot_all(sol)

    symmetry_test()
    scale_invariance_test()
    grid_refinement_test()
    comparative_statics()
    damping_diagnostic()
    generalized_dims_test()

    # Nesting assertions, the ψ sweep and the goods moving cost (section 6).
    sorting_mechanism_tests()

    # Verbose mechanism analyses + figures: runs on the baseline `sol` plus its
    # own parameter sweeps (section 6).
    mechanism_analyses(sol)
    return sol
end

# Running the file as a script executes the whole battery; `include`ing it from a
# REPL or another script just brings the diagnostics into scope.
if abspath(PROGRAM_FILE) == @__FILE__
    @time sol = main_test()
end
