# Notes on `spatial_ge_analysis.jl`

Review of the current GE test output and the test code. The suite now runs clean:
every structural identity holds to machine precision, all twelve comparative-statics
cases converge, the damping diagnostic behaves as the contraction theory predicts,
and symmetry / scale-invariance / generalized-dims / monotonicity all pass. The
earlier convergence problems (grid-dependent stalls, backwards damping advice,
comparative-statics "NO" flags) are gone. What remains are a few modeling
observations and minor test-hygiene items, none of them failures.

## Current headline (Nz=3, Nϵ=9, damping=0.5)

```
H̃_T = (0.0378, 0.1121)     M = (0.7980, 1.2020)   [total 2.0000]
t   = (0.1424, 0.3130)     Q = (0.5548, 0.6571)
mean z: loc 1 = 1.0445   loc 2 = 1.0365   ergodic = 1.0397
migration: born 1 → (0.7219, 0.2781)   born 2 → (0.1846, 0.8154)
teaching share: born 1 = 0.2297   born 2 = 0.3332
```

Equilibrium-consistency checks: all `[ok]`. Structural identities ~1e-11 or
tighter (Φ eigen-residual 3.15e-11, row sums / mass / Q at machine ε). The
convergence-sensitive conditions are well inside tolerance — GE-map residual
4.16e-05, M_l = 2∫Φ_l 3.32e-05, household policy residual 2.93e-09 — and the
budget clamp is inactive. Monotonicity is full (s, e, V all 72/72 in z and 24/24
in ϵ).

## Resolved since the previous notes

1. **Convergence is now uniform.** All twelve comparative-statics cases report
   `conv? = yes` (residuals ~7e-4, except high-β at 1.4e-3). The old failures
   (`zero-β`, `high-persistence`, `low-persistence`) no longer appear; the
   "read these aggregates with care" caveat no longer fires.
2. **Damping diagnostic now matches the theory.** Residual falls monotonically as
   damping rises — 0.20→7.77e-5, 0.30→7.24e-5, 0.50→4.15e-5, 0.70→1.54e-5,
   0.80→1.26e-5 — confirming a contracting outer map where a larger step converges
   faster. Every damping value is usable (<5e-3), and the script's advice ("use
   smaller damping only if iterates diverge") is now the correct direction.
3. **Grid / threshold reconciliation.** The damping diagnostic and comparative
   statics both run at Nϵ=7, hh_maxit=1000, and the "usable" / "conv?" cutoffs are
   both 5e-3. The earlier Nϵ=15-vs-7 mismatch and the 1e-2-vs-5e-3 threshold split
   are gone.
4. **Scale-invariance is tight.** Solves run to tol=1e-4 and the deltas are now
   ~1e-6 to 1e-8 (M 8.51e-06, H̃_T 7.66e-07, Q 2.69e-08, t 3.54e-07, Gₗ 7.06e-09),
   comfortably below the 1e-3 band rather than just scraping under it.
5. **Sign checks guard against non-converged inputs.** The low-σν directional
   check is skipped unless that case converged (it does here), so a stalled solve
   can no longer silently feed a sign test.

## Still worth noting

1. **Degenerate corners pass silently.** `high-β` drives H̃_T → (0.000, 0.000)
   with teaching share 0.997, and `high-teach-wage` reaches 0.979. These are
   genuine equilibria but economically extreme (near-universal teaching with
   teacher human capital collapsing), and nothing in the sign checks flags the
   corner. Worth a guard if these parameterizations matter.
2. **Endogenous sorting is weak.** At the headline, G_1 / G_2 / ergodic differ by
   only ~0.002–0.006 (e.g. at z=1.486: 0.2558 / 0.2461 / 0.2500). In the I=3/L=3
   test G_l is essentially flat against the ergodic law (0.4978–0.5022 vs 0.5000).
   The location-specific ability distribution is barely active — confirm this is
   intended rather than a muted channel.
3. **Headline-grid comment is stale.** `main_test` solves at Nϵ=9, but the comment
   block above it still says "Nz=3, Nϵ=7" (and the older line says Nϵ=15). Cosmetic,
   but the docstring no longer matches the code.
4. **Grids still vary across batteries.** Headline Nϵ=9, comp-stats / damping Nϵ=7,
   symmetry / scale Nϵ=5. Defensible for runtime, but residuals across sections
   aren't strictly comparable.
5. **Scaling-cost test remains small.** With reps=3 the ratios are steadier than
   before (base 1.00×, +1 loc 1.70×, +2 loc 2.46×, +1 occ 1.55×, +2 occ 2.51×,
   +1 each 2.45×), but absolute times are still 0.05–0.13s and the geometric-in-I
   signal is only loosely visible (+2 occupations at 20000 states ≈ +2 locations at
   800 states). Bigger dimensions would make the Nϵ^I scaling clearer.

## Generalized dimensions (I=3, L=3)

Runs and passes the structural identities (Φ eig-res 9.77e-12, Σπ and Σ_l M_l
deviations at machine ε). H̃_T = (0.0008, 0.0062, 0.0062), M = (0.2348, 0.8826,
0.8826) — the two symmetric locations 2 and 3 coincide exactly, as expected.
