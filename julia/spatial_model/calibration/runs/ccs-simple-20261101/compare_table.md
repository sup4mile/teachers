# CCS simple-baseline comparison (2026-11-01)

All rows are exact refits (8 parameters ϑ = (logκ̃, δκ, ΔB, β, λ, r_m, ωf, γ), 8 targets) from the September fit, ψ = 0, default solver (Nz 5, 48/48). ρ = corr(ln θ, ln h_T) of the CFR effectiveness concept θ = h_T e^ν (--measure cfr_rho; under CCS also the draw, --external ccs_rho); K_d = independent teacher draws over K = 12 years. Built by compare_table.py from the run directories named in LOG.md.

## 1. Parameters and fit

| model | K_d | ρ | logκ̃ | δκ | ΔB | β | λ | r_m | ωf | γ | criterion (8 targets) | max \|z\| targeted | condition no. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| efficient (September) | eff | 1 | -1.2038 | 0.0045 | 0.0662 | 0.2278 | 0.1245 | 0.2059 | 0.9664 | 0.6194 | 1.61e-19 | 3.99e-10 | 30.3 |
| efficient, ρ 0.3 | eff | 0.3 | -0.9920 | 0.0048 | 0.0662 | 0.0683 | 0.0713 | 0.3594 | 0.9664 | 0.6194 | 5.40e-18 | 1.62e-09 | 18.7 |
| **simple CCS** (K_d 1, ρ 1) | 1 | 1 | -1.2553 | 0.0032 | 0.0662 | 0.2108 | 0.1241 | 0.2059 | 0.9812 | 0.5694 | 1.51e-17 | 3.04e-09 | 28.1 |
| CCS + θ (K_d 1, ρ 0.3) | 1 | 0.3 | -1.0129 | 0.0037 | 0.0662 | 0.0631 | 0.0706 | 0.3621 | 0.9815 | 0.5684 | 3.58e-18 | 1.85e-09 | 19.2 |
| CCS + K (K_d 12, ρ 1) | 12 | 1 | -1.2182 | 0.0042 | 0.0662 | 0.2261 | 0.1278 | 0.2007 | 0.9680 | 0.6145 | 4.42e-19 | 5.53e-10 | 30.2 |
| CCS + θ + K (K_d 12, ρ 0.3) | 12 | 0.3 | -0.9942 | 0.0047 | 0.0662 | 0.0678 | 0.0714 | 0.3590 | 0.9680 | 0.6145 | 5.03e-18 | 1.62e-09 | 18.7 |

## 2. Targeted moments (the 8; matched exactly in every row)

| model | male_teach_share | female_teach_share | teacher_p9010 | gap_salary | gap_pupils | cfr_effect | wtp_share | move_rate |
|---|---|---|---|---|---|---|---|---|
| Data (se) | 0.0194 (0.0019) | 0.0599 (0.0060) | 2.7728 (0.1386) | 0.0086 (0.0092) | 0.3764 (0.1154) | 0.0133 (0.0041) | 0.0060 (0.0022) | 0.2156 (0.0500) |
| efficient (September) | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |
| efficient, ρ 0.3 | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |
| **simple CCS** (K_d 1, ρ 1) | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |
| CCS + θ (K_d 1, ρ 0.3) | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |
| CCS + K (K_d 12, ρ 1) | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |
| CCS + θ + K (K_d 12, ρ 0.3) | 0.0194 | 0.0599 | 2.7728 | 0.0086 | 0.3764 | 0.0133 | 0.0060 | 0.2156 |

## 3. Untargeted validation (S13): model value (z = (model − data)/se)

None of these is a target in any row. Block subtotals are Σz² over the block; S13 is their sum.

### 3a. Sorting and teacher allocation

| model | V1 earnings gap | V3 child-poverty gap | V4 SEDA score gap | V6a city share by parental income (ACS within-CZ) | V6z city share by own AFQT (NLSY79) | V2 FTE/pupil gap | V5 SEDA growth gap | payroll per pupil gap |
|---|---|---|---|---|---|---|---|---|
| Data (se) | 0.1122 (0.0168) | -0.0888 (0.0090) | 0.2887 (0.0326) | -0.187 (0.0170) | -0.222 (0.0360) | -0.0177 (0.0131) | 0.0031 (0.0030) | -0.0091 (0.0188) |
| efficient (September) | -0.0102 (-7.3) | 0.0024 (10.1) | 0.0605 (-7.0) | 0.0108 (11.6) | 0.0027 (6.2) | 0.0234 (3.1) | 0.0119 (2.9) | 0.0320 (2.2) |
| efficient, ρ 0.3 | -0.0104 (-7.3) | 0.0026 (10.2) | 0.0592 (-7.0) | 0.0110 (11.6) | 0.0027 (6.2) | 0.0250 (3.3) | 0.0117 (2.9) | 0.0336 (2.3) |
| **simple CCS** (K_d 1, ρ 1) | -0.0123 (-7.4) | 0.0031 (10.2) | 0.0416 (-7.6) | 0.0123 (11.7) | 0.0028 (6.2) | 0.0171 (2.7) | 0.0100 (2.3) | 0.0257 (1.9) |
| CCS + θ (K_d 1, ρ 0.3) | -0.0127 (-7.4) | 0.0033 (10.2) | 0.0395 (-7.6) | 0.0127 (11.7) | 0.0028 (6.2) | 0.0199 (2.9) | 0.0095 (2.1) | 0.0285 (2.0) |
| CCS + K (K_d 12, ρ 1) | -0.0102 (-7.3) | 0.0024 (10.1) | 0.0580 (-7.1) | 0.0109 (11.6) | 0.0027 (6.2) | 0.0220 (3.0) | 0.0118 (2.9) | 0.0306 (2.1) |
| CCS + θ + K (K_d 12, ρ 0.3) | -0.0106 (-7.3) | 0.0026 (10.2) | 0.0562 (-7.1) | 0.0112 (11.7) | 0.0027 (6.2) | 0.0245 (3.2) | 0.0115 (2.8) | 0.0331 (2.2) |

### 3b. Teacher pay and first stage

| model | V9 teacher premium, men | V9 teacher premium, women | V11 market gender gap | V12 latent wage slope | V12 mother–child ρ |
|---|---|---|---|---|---|
| Data (se) | 0.0557 (0.0500) | 0.086 (0.0500) | 0.1046 (0.0500) | 0.2127 (0.0077) | 0.5766 (0.0140) |
| efficient (September) | -0.7194 (-15.5) | -0.6311 (-14.3) | 0.3023 (4.0) | 0.2138 (0.1) | 0.5760 (-0.0) |
| efficient, ρ 0.3 | -0.7194 (-15.5) | -0.6311 (-14.3) | 0.3023 (4.0) | 0.2138 (0.1) | 0.5760 (-0.0) |
| **simple CCS** (K_d 1, ρ 1) | -0.8587 (-18.3) | -0.7581 (-16.9) | 0.2824 (3.6) | 0.2581 (5.9) | 0.2976 (-19.9) |
| CCS + θ (K_d 1, ρ 0.3) | -0.8612 (-18.3) | -0.7602 (-16.9) | 0.2820 (3.5) | 0.2582 (5.9) | 0.2976 (-19.9) |
| CCS + K (K_d 12, ρ 1) | -0.7319 (-15.8) | -0.6428 (-14.6) | 0.3003 (3.9) | 0.2157 (0.4) | 0.5178 (-4.2) |
| CCS + θ + K (K_d 12, ρ 0.3) | -0.7319 (-15.8) | -0.6427 (-14.6) | 0.3003 (3.9) | 0.2158 (0.4) | 0.5178 (-4.2) |

### 3c. Block subtotals (Σz², untargeted)

| model | sorting | teacher allocation | teacher pay | first stage | **S13** |
|---|---|---|---|---|---|
| efficient (September) | 379.2 | 23.2 | 446.0 | 15.7 | **864.1** |
| efficient, ρ 0.3 | 380.5 | 24.0 | 446.0 | 15.7 | **866.2** |
| **simple CCS** (K_d 1, ρ 1) | 393.0 | 15.8 | 619.5 | 444.6 | **1472.8** |
| CCS + θ (K_d 1, ρ 0.3) | 395.5 | 16.8 | 622.7 | 444.7 | **1479.6** |
| CCS + K (K_d 12, ρ 1) | 380.4 | 22.2 | 460.6 | 33.1 | **896.2** |
| CCS + θ + K (K_d 12, ρ 0.3) | 382.4 | 23.2 | 460.6 | 33.1 | **899.3** |

## 4. CCS mechanism and the teacher-quality gap (untargeted)

ρ_pc, wage slope: the first-stage rows of S13 (z in brackets). Draw share: c² E_l Var(ln r | l)/sd_S², the teacher draw's share of the adult latent score variance (`ccs_score_draw_share`). SD(ln q) = SD(ln r | l) by location (city, suburb). Δ log Q̄_T: suburb − city gap in the log teacher index, split into teacher quality (β Δμ) and staffing (−σ Δ ln N). ΔΛ: the suburb's quality value (Φ-weighted mean) and, under CCS, its dispersion part.

| model | ρ_pc | wage slope | draw share | SD(ln q) | class size N_l | ln Q^eff − ln Q̄ | Δ log Q̄_T (quality / staffing) | ΔΛ quality | ΔΛ disp | children's Δ log z (SD) | teachers' z vs market (SD) | WTP income elast. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Data | 0.5766 (se 0.014) | 0.2127 (se 0.0077) | — | — | — | — | — | — | — | — | 0.614 | 0.379 |
| efficient (September) | 0.5760 (-0.0) | 0.2138 (0.1) | — | — | — | — | 0.0149 (0.0018 / 0.0131) | 0.0017 | — | -0.0050 | -0.2864 | 1.094 |
| efficient, ρ 0.3 | 0.5760 (-0.0) | 0.2138 (0.1) | — | — | — | — | 0.0146 (6.24e-04 / 0.0140) | 0.0016 | — | -0.0050 | -0.2863 | 1.094 |
| **simple CCS** (K_d 1, ρ 1) | 0.2976 (-19.9) | 0.2581 (5.9) | 0.672 | (0.1508, 0.1487) | (25.47, 25.04) | (0.0205, 0.0199) | 0.0125 (0.0029 / 0.0096) | 0.0013 | -3.99e-05 | -0.0051 | -0.3363 | 1.091 |
| CCS + θ (K_d 1, ρ 0.3) | 0.2976 (-19.9) | 0.2582 (5.9) | 0.672 | (0.1508, 0.1487) | (25.52, 25.01) | (0.0018, 0.0018) | 0.0119 (8.07e-04 / 0.0111) | 0.0013 | -3.91e-05 | -0.0052 | -0.3373 | 1.091 |
| CCS + K (K_d 12, ρ 1) | 0.5178 (-4.2) | 0.2157 (0.4) | 0.146 | (0.0435, 0.0429) | (25.55, 24.99) | (0.0205, 0.0199) | 0.0148 (0.0025 / 0.0123) | 0.0016 | -3.37e-06 | -0.0050 | -0.2912 | 1.094 |
| CCS + θ + K (K_d 12, ρ 0.3) | 0.5178 (-4.2) | 0.2158 (0.4) | 0.146 | (0.0435, 0.0429) | (25.59, 24.97) | (0.0018, 0.0018) | 0.0144 (6.89e-04 / 0.0137) | 0.0016 | -3.32e-06 | -0.0050 | -0.2912 | 1.094 |

Other: CFR SD of log effectiveness `cfr_sd_eff` efficient (September) 0.6565, efficient, ρ 0.3 2.1883, simple CCS (K_d 1, ρ 1) 0.7093, CCS + θ (K_d 1, ρ 0.3) 2.3681, CCS + K (K_d 12, ρ 1) 0.6613, CCS + θ + K (K_d 12, ρ 0.3) 2.2043.

## 5. Draw-share formula check

Analytic: CFR pins SD(ln q) ≈ K(1 − η)e/√K_d (e = 0.013344 the CFR target, K = 12, η = 0.08) and the loading c = s_z/((1 − η)χ), so the draw's share of the adult score variance is (Ke)²/(K_d χ² + (Ke)²), independent of s_z (χ = log 1.12). Semi-analytic: c² v̄/(s_z² + c² v̄) with the model's v̄ = Σ_l pop_l SD(ln q | l)² and s_z = c(1 − η)χ (so the non-draw score variance is set to s_z²). Model: `ccs_score_draw_share`.

| model | K(1−η)e/√K_d | model SD(ln q) | analytic share | semi-analytic share | model share |
|---|---|---|---|---|---|
| **simple CCS** (K_d 1, ρ 1) | 0.1473 | (0.1508, 0.1487) | 0.666 | 0.673 | 0.672 |
| CCS + θ (K_d 1, ρ 0.3) | 0.1473 | (0.1508, 0.1487) | 0.666 | 0.673 | 0.672 |
| CCS + K (K_d 12, ρ 1) | 0.0425 | (0.0435, 0.0429) | 0.143 | 0.146 | 0.146 |
| CCS + θ + K (K_d 12, ρ 0.3) | 0.0425 | (0.0435, 0.0429) | 0.143 | 0.146 | 0.146 |

## 6. Grid-pad check (padcheck.jl at the converged point)

| case | max \|x\| | pad | covered | nodes x (city; suburb) |
|---|---|---|---|---|
| ccs12_psi0 | 0.0768 | 0.30 | True | -0.0737, 8.90e-04, 0.0768; -0.0727, 8.84e-04, 0.0758 |
| ccs12_rho03 | 0.0754 | 0.30 | True | -0.0754, 2.37e-05, 0.0754; -0.0743, 2.37e-05, 0.0744 |
| ccs1_psi0 | 0.2722 | 0.60 | True | -0.2393, 0.0096, 0.2722; -0.2384, 0.0089, 0.2691 |
| ccs1_rho03 | 0.2616 | 0.60 | True | -0.2607, 2.54e-04, 0.2616; -0.2571, 2.36e-04, 0.2579 |
| rhoz_075 | 0.2721 | 0.60 | True | -0.2393, 0.0096, 0.2721; -0.2384, 0.0089, 0.2692 |
| rhoz_090 | 0.2721 | 0.60 | True | -0.2393, 0.0096, 0.2721; -0.2384, 0.0089, 0.2692 |
| rhoz_097 | 0.2721 | 0.60 | True | -0.2393, 0.0096, 0.2721; -0.2384, 0.0089, 0.2692 |
