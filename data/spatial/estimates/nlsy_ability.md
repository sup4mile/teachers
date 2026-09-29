# NLSY ability moments for the first-stage indirect inference (T2a)

Generated 2026-09-28T18:15:41 by `nlsy_ability.py`; cluster bootstrap B = 200 (seed 20260928), SE = bootstrap SD, CI = 95% percentile. Data side of T2b: every statistic is to be recomputed on simulated model panels as defined in the last section.

## What feeds T2b

| key | quantity | estimate (SE) | 95% CI | n |
|---|---|---|---|---|
| `rho_pc_latent` | Latent mother-child ability correlation (two-factor ULS) | 0.577 (0.014) | [0.548, 0.601] | 8137 |
| `wage_slope_latent_pooled` | Log-wage slope per SD of latent ability (pooled) | 0.196 (0.007) | [0.185, 0.213] | 34118 |
| `R_comp_mother` | Reliability of the equal-weight mother composite (4 IRT z) | 0.915 (0.002) | [0.910, 0.919] | 3480 |
| `R_comp_child` | Reliability of the child composite actually used (mean of available tests) | 0.892 (0.003) | [0.886, 0.898] | 8137 |
| `R_comp_pop` | Reliability of the 4-subtest AFQT composite, NLSY79 population | 0.916 (0.002) | [0.912, 0.919] | 9033 |
| `implied_s_z` | Implied stationary SD of log ability, s_z = slope_latent (1−η) | 0.181 (0.007) | [0.170, 0.195] | 34118 |
| `gap_to_chi` | Pooled latent slope minus χ = log 1.12 | 0.083 (0.007) | [0.072, 0.099] | 34118 |

η = 0.0807, χ = log 1.12 = 0.11333. `rho_pc_latent` is the target for ρz; `wage_slope_latent_pooled` is the target for b·s_z; the R_comp's set var(u) in the score measurement equation (mother composite, child composite, population composite).

## A. NLSY79 population measurement model

Sample: SAMPLE_ID 1-8, 10, 11, 13, 14 (drops poor-white supplement 9/12 and military 15-20), all four IRT z valid, N = 9033. Weight = 1979 SAMPWEIGHT. IRT z-scores are already age-normed by NLS (four-month groups within birth year).

| subtest | loading | SE | communality | KR-20 α (number-right) |
|---|---|---|---|---|
| AR | 0.882 | 0.004 | 0.778 | 0.908 |
| WK | 0.857 | 0.004 | 0.734 | 0.924 |
| PC | 0.805 | 0.005 | 0.647 | 0.809 |
| MK | 0.876 | 0.004 | 0.767 | 0.895 |

R_comp_pop = 0.916 (0.002); omega = 0.916 (0.002); RMSR = 0.0373; Heywood: False. KR-20 is for the number-right score, not the IRT score: a benchmark for the communalities. Own-norming sensitivity (section standard scores ranked within birth-year x 4-month cells): corr(afqt_c, afqt_c_own) = 0.993 (0.001), R_comp_own = 0.909.

## B. Mother-child measurement system

Mothers: SAMPLE_ID 5-8, 13, 14 with four valid IRT z; children: CNLSY, PIAT math / recognition / comprehension and PPVT-R at ages 5-14 (60-179 months), rank-normed within 3-month age cells over all CNLSY assessments (reference weight = round child weight), averaged within child. Dyad weight = mother's 1979 weight / number of her children in the sample. N mothers = 3480, N children = 8137; mean assessment rounds per child = 3.86; share with PPVT = 0.93. Child age is `CSAGE` (child supplement), with `MSAGE` only where CSAGE is missing (0 assessments). **Imputed comprehension:** NLS gives children with RECOG < 19 no comprehension test and copies the recognition score into COMP (5131 in-window assessments, 15.5% of in-window COMP assessments). Those rows stay in the COMP norming reference (so real scores are ranked against the whole age cell) but are excluded from child-level COMP means and from the COMP two-year stabilities. `R_comp_child` is the reliability of the mean of the available tests.

| key | quantity | estimate (SE) | 95% CI | n |
|---|---|---|---|---|
| `rho_pc_latent` | Latent mother-child ability correlation (two-factor ULS) | 0.577 (0.014) | [0.548, 0.601] | 8137 |
| `R_comp_mother` | Reliability of the equal-weight mother composite (4 IRT z) | 0.915 (0.002) | [0.910, 0.919] | 3480 |
| `R_comp_child` | Reliability of the child composite actually used (mean of available tests) | 0.892 (0.003) | [0.886, 0.898] | 8137 |
| `r_pc_composite` | Weighted correlation of the two equal-weight composites | 0.516 (0.013) | [0.485, 0.538] | 8137 |

| indicator | loading | SE | communality |
|---|---|---|---|
| AR | 0.876 | 0.007 | 0.767 |
| WK | 0.877 | 0.007 | 0.770 |
| PC | 0.802 | 0.009 | 0.643 |
| MK | 0.860 | 0.007 | 0.740 |
| MATH | 0.832 | 0.008 | 0.692 |
| RECOG | 0.826 | 0.007 | 0.682 |
| COMP | 0.860 | 0.007 | 0.740 |
| PPVT | 0.789 | 0.010 | 0.623 |

RMSR = 0.0418; Heywood: False; omega (mother, child) = 0.915, 0.896; check r/sqrt(R_M R_C) = 0.571 (SE 0.014) vs the ULS rho above.

Two-year child stabilities (not used):

| test | corr |
|---|---|
| MATH | 0.661 (0.007) |
| RECOG | 0.740 (0.006) |
| COMP | 0.635 (0.008) |
| PPVT | 0.750 (0.015) |

### Sensitivities

| variant | rho_pc_latent | R_comp_mother | R_comp_child | N mothers | N children |
|---|---|---|---|---|---|
| Main | 0.577 (0.014) | 0.915 (0.002) | 0.892 (0.003) | 3480 | 8137 |
| Cross-sectional mothers only (SAMPLE_ID 5-8) | 0.529 (0.018) | 0.905 (0.003) | 0.886 (0.004) | 2140 | 4796 |
| Unweighted | 0.630 (0.011) | 0.922 (0.002) | 0.905 (0.002) | 3480 | 8137 |
| Child weights (mean CSAMWT_REV), no 1/n_c | 0.611 (0.013) | 0.920 (0.003) | 0.901 (0.003) | 3480 | 8137 |
| Children born 1981-2000 | 0.584 (0.015) | 0.914 (0.003) | 0.893 (0.003) | 3180 | 6436 |
| One child per mother (first-born) | 0.591 (0.017) | 0.915 (0.002) | 0.897 (0.004) | 3480 | 3480 |
| Child scores from ages 10-14 only | 0.593 (0.013) | 0.915 (0.002) | 0.875 (0.004) | 3258 | 7393 |
| Linear age-norming | 0.569 (0.014) | 0.915 (0.002) | 0.888 (0.003) | 3480 | 8137 |
| Drop within-domain pairs AR-MK, WK-PC, RECOG-COMP | 0.605 (0.015) | 0.887 (0.004) | 0.841 (0.005) | 3480 | 8137 |
| Mothers' own-normed section scores | 0.574 (0.014) | 0.905 (0.003) | 0.893 (0.003) | 3480 | 8137 |
| Drop only the RECOG-COMP pair | 0.595 (0.014) | 0.915 (0.002) | 0.841 (0.005) | 3480 | 8137 |

## S. Siblings

Families with ≥ 2 children in the dyad sample: 2688; sibling pairs: 7720. Ordered pairs, each family's weight = mother's weight.

| key | quantity | estimate (SE) | 95% CI | n |
|---|---|---|---|---|
| `rho_sib_latent` | Latent sibling correlation (off-diagonal cross-test) | 0.606 (0.017) | [0.573, 0.635] | 7720 |
| `S.sib_comp_corr` | Raw sibling correlation of the child composite | 0.547 (0.016) | [0.515, 0.576] | 7720 |
| `implied_rho_from_sib` | Implied rho = rho_sib_latent / rho_pc_latent | 1.051 (0.031) | [0.994, 1.112] | 7720 |

Including k = l terms: rho_sib_latent = 0.631 (0.016).

Cross-sibling cross-test correlations C[k, l] (child k, sibling l):

| | MATH | RECOG | COMP | PPVT |
|---|---|---|---|---|
| MATH | 0.508 | 0.415 | 0.416 | 0.433 |
| RECOG | 0.415 | 0.462 | 0.434 | 0.378 |
| COMP | 0.416 | 0.434 | 0.447 | 0.410 |
| PPVT | 0.433 | 0.378 | 0.410 | 0.524 |

## C. Log-wage slope on latent AFQT (NLSY79)

Person-years of the block A sample, survey rounds with ages 25-34; income > 0, weeks ≥ 50, hours/week ≥ 35, not armed forces (ESR ≠ 4), highest grade ≥ 9; log(income / annual hours); weighted 1/99 trim within year; weight = 1979 SAMPWEIGHT; controls: age and year dummies (+ female).

| slope | pooled | male | female |
|---|---|---|---|
| observed, per SD of afqt_c | 0.188 (0.007) | 0.171 (0.009) | 0.215 (0.009) |
| latent, per SD of f | 0.196 (0.007) | 0.179 (0.009) | 0.224 (0.009) |
| N person-years | 34118 | 19245 | 14873 |
| N persons | 6904 | 3667 | 3237 |
| 90/10 (year-demeaned) | 3.40 (0.06) | 3.43 (0.08) | 3.20 (0.05) |

Implied s_z = 0.181 (0.007); gap to χ = 0.083 (0.007).

Participation (LPM slope per SD of afqt_c; all block A person-years, hgc ≥ 9, not armed forces, no labor-supply filter):

| outcome | pooled | male | female |
|---|---|---|---|
| 1{FYFT wage worker} | 0.079 (0.005) | 0.077 (0.005) | 0.080 (0.006) |
| 1{hours < 780} | -0.063 (0.004) | -0.043 (0.004) | -0.084 (0.006) |

Sensitivities (pooled latent slope):

| variant | slope (SE) |
|---|---|
| Main | 0.196 (0.007) |
| Cross-sectional sample only (SAMPLE_ID 1-8) | 0.188 (0.008) |
| Fully unweighted (incl. standardization and R_comp_pop) | 0.209 (0.006) |
| Own-normed afqt (A4), Block A sample | 0.196 (0.007) |
| Hourly wage from HRP1 | 0.184 (0.006) |

## V. NLSY97 vintage check

CAT-ASVAB thetas age-normed within birth-year x quarter (1997 weight, fixed across draws), one-factor ULS: R_comp_97 = 0.912 (0.002) (N = 7093). Wage sample as block C with `YINC-1700` and CVC hours/weeks of calendar year Y−1; both NLSY97 samples, weighted by 1997 SAMPLING_WEIGHT_CC; bootstrap by person.

| | pooled | male | female |
|---|---|---|---|
| latent slope | 0.175 (0.009) | 0.159 (0.011) | 0.203 (0.014) |
| observed slope | 0.168 (0.008) | 0.152 (0.011) | 0.194 (0.014) |
| N person-years | 17543 | 9660 | 7883 |

| subtest | loading | posterior-variance reliability |
|---|---|---|
| AR | 0.878 | 0.768 |
| WK | 0.812 | 0.751 |
| PC | 0.853 | 0.713 |
| MK | 0.856 | 0.777 |

## E. Earnings checks (not targets)

Permanent earnings = mean over ages 25-34 of log wage-and-salary income residualized on age and year dummies within generation (and, for young adults, sex), after 1/99 trimming within year.

| statistic | estimate (SE) | n |
|---|---|---|
| E1 mother-child correlation | 0.135 (0.018) | 5529 |
| E1 covariance | 0.098 (0.014) | |
| E1 slope (child on mother) | 0.107 (0.014) | |
| E1 correlation, male / female children | 0.113 (0.025) / 0.155 (0.024) | 2769 / 2760 |
| E1 ≥ 2 years each: correlation, slope | 0.129 (0.022), 0.103 (0.018) | |
| E2 sibling correlation | 0.281 (0.023) | 4319 |
| E3 implied rho (E2 / E1) | 2.090 (0.393) | |

## Definitions for T2b

In every case the model panel must use the same sample, ages, weights and composites. 'Score signal' means `θ0 + θ1 (log z + c log Q_l)` without the noise u; var(u) comes from the R_comp's.

- **`rho_pc_latent`** (ρz). Data: mothers' scores are the 1980 ASVAB taken at ages 15-23 (four IRT z-scores, age-normed by NLS within birth cohort); children's scores are PIAT math, recognition, comprehension and PPVT-R at ages 5-14, each age-normed within 3-month cells and averaged over the child's assessment rounds (about 3.9 rounds per child, so the child factor is a multi-round average; imputed comprehension scores excluded). Estimate: two-factor ULS on the 8x8 weighted pairwise-complete correlation matrix, weights = mother's 1979 weight / n_c. Model counterpart: the correlation between the latent score signals of a parent and a child, one child per parent, parents weighted by population mass. Includes assortative mating and the Q_l channel because the data factor does too.
- **`wage_slope_latent_pooled`** (b·s_z). `afqt_c` = equal-weight mean of the four standardized IRT z-scores, re-standardized (weighted) in the Block A POPULATION sample (all persons with four valid scores, any labor-market state; NOT the wage sample). Pooled person-years at ages 25-34 at interview, with income, weeks and hours for the previous calendar year; filters income > 0, weeks ≥ 50, hours/weeks ≥ 35, ESR ≠ 4, highest grade ≥ 9; weighted 1st/99th percentile trim of log hourly wage within survey year; weighted OLS on afqt_c with female, age and year dummies; slope divided by sqrt(R_comp_pop). Model counterpart: the same regression on the noisy simulated composite standardized in the simulated population, divided by the square root of its reliability, or equivalently the slope on the noise-free signal standardized in the population. Standardizing in the wage sample instead would move the slope by about 7% (SD of afqt_c in the full-year full-time sample is 0.94). `implied_s_z` = slope × (1 − η).
- **`R_comp_mother`, `R_comp_pop`**: (Σλ)² / (1'S1) for the equal-weight composite of the four standardized IRT z-scores, S the sample correlation matrix, in the dyad sample and the population sample respectively. **`R_comp_child`**: reliability of the composite actually used, the mean of each child's available standardized tests, R = Σ_i w_i (mean_{k∈K_i} λ_k)² / Σ_i w_i (1'S_{K_i}1 / |K_i|²) over dyad children (`R_comp_child_all4` in the JSON assumes all four tests present). The model composite's noise variance is set so that the simulated composite's reliability matches these.
- **`rho_sib_latent`**: Σ_{k≠l} C_kl λ_k λ_l / Σ_{k≠l} (λ_k λ_l)² over ordered sibling pairs (mother-level weights); a check, expected ≥ ρ² but sibling-shared inputs are outside the model.
- **Earnings checks (E1-E3)**: permanent earnings = mean of age-and-year-residualized log wage-and-salary income over ages 25-34 per person; correlations weighted with mother's weight / (number of children with earnings).
- **Age norming**: data scores are age-normed to N(0,1) before averaging. In the model, scores are already age-free.
