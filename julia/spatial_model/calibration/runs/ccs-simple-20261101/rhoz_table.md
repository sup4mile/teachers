# ρ_pc ceiling of the simple model: re-estimating ρz (2026-11-01)

Fixed-ϑ evaluations (PANEL_MAX_ITER=0) at the simple model's exact refit (runs/exp-ccs/refit/ccs1_psi0: K_d 1, ρ 1, ψ 0) with the AR(1) persistence ρz of log ability raised; `--external rho_z` holds the stationary SD s_z (σξ = s_z √(1 − ρz²)). s_z is read back from the score loading, s_z = c(1 − η)χ. Prediction: with the draw independent across generations and of z, ρ_pc ≈ ρz √((1 − a)(1 − b)), a the adult draw share and b the child's at grade 4 (b = (4/12)² c²v̄/(s_z² + (4/12)² c²v̄)). The criterion is that of the 8 targets at the unchanged ϑ (the refit is exact at ρz = 0.5766). Data: ρ_pc 0.5766 (se 0.014), wage slope 0.2127 (se 0.0077); the latent sibling correlation in the NLSY is 0.606 vs the model's ρz².

| case | ρz | s_z | σξ | ρz² | ρ_pc (z) | predicted ρ_pc | ρ_pc all parents | wage slope (z) | draw share | SD(ln q) | criterion (8 targets) | max \|z\| targeted | first-stage Σz² | S13 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| refit (ρz 0.5766) | 0.5766 | 0.1956 | 0.1598 | 0.332 | 0.2976 (-19.9) | 0.2978 | 0.2975 | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | 1.51e-17 | 3.04e-09 | 444.6 | 1472.8 |
| rhoz_05766 (control) | 0.5766 | 0.1956 | 0.1598 | 0.332 | 0.2976 (-19.9) | 0.2978 | 0.2975 | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | 1.51e-17 | 3.04e-09 | 444.6 | 1472.8 |
| rhoz_075 | 0.75 | 0.1956 | 0.1294 | 0.562 | 0.3870 (-13.5) | 0.3874 | 0.3869 | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | 0.001 | 0.02 | 230.9 | 1261.6 |
| rhoz_090 | 0.9 | 0.1956 | 0.0853 | 0.810 | 0.4643 (-8.0) | 0.4649 | 0.4642 | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | 0.005 | 0.06 | 111.8 | 1145.4 |
| rhoz_097 | 0.97 | 0.1956 | 0.0476 | 0.941 | 0.5004 (-5.4) | 0.5010 | 0.5003 | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | 0.008 | 0.08 | 77.1 | 1112.6 |

Ceiling: at ρz → 1, ρ_pc → √((1 − a)(1 − b)) ≈ 0.517, z = -4.3; matching 0.5766 would need ρz ≈ 1.116 > 1.
