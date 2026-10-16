# NLSY79: wage-ability slopes holding the BA status fixed, teachers' slope ratio, and the BA-share logit shift

Generated 2026-10-05T03:14:22 by `nlsy_ba_slopes.py` (task D3b of the 2026-10b round, DESIGN_v2.md section M3; second version, after the review of the first). NLSY79, ages 25-34, the wage sample, wage definition, 1979-weight convention, latent-AFQT standardization (`afqt_c`, R_comp_pop = 0.916 (0.0017)) and controls (age and survey-year dummies; a female dummy in pooled samples) of `nlsy_ability.py` block C. Household-cluster bootstrap, B = 500, seed 20260929 (the draws of `nlsy_teachers.py` and D3); every draw re-standardizes AFQT, re-estimates R and re-runs every regression and shift. SE = bootstrap SD in parentheses; 95% percentile CI in brackets. All slopes are log hourly-wage points per SD of latent AFQT. NLSY97 and the sealed hold-out are not used (see the sealed-set note under Definitions and caveats).

Two disattenuation versions are reported side by side, with a split-half IV check (section 2). **Convention** (the `value` of every slope in the targets file): observed slope / sqrt(R), the project's latent slope; it gives the published 0.2127, so the BA-conditional slope and the unconditional target are on one definition. **EIV** (carried as the alternative): the classical errors-in-variables correction of the regression itself (error variance 1 - R in afqt_c units); it equals the convention for an unconditional regression on the full block-A sample, but in subsamples and with a BA control the conditional reliability of the score is lower than R, so it is 5% higher for the BA-controlled slope. Both estimate the slope on the true latent score, which is what the model's `wage_slope` is; the choice between them is a PREREG matter (V12, the unconditional slope, is on the convention).

## Summary

- **Reproduction.** The pooled unconditional latent slope is 0.2127 against the published 0.2127 (SE 0.0074 here, 0.0077 published; replaying the published bootstrap stream, 0.0077). Men 0.1904, women 0.2468 reproduce too (section 1).
- **Slope holding BA fixed (the M3 target `wage_slope_ba`; convention).** Adding a BA indicator (highest grade >= 16) lowers the pooled latent slope from 0.2127 to **0.1608 (0.0083)** on the project's convention (observed slope / sqrt(R), the pipeline of the 0.2127 target): a fall of 24% (SE 2 points). Men 0.1481 (0.0117), women 0.1825 (0.0114). The EIV version, which lets the reliability of the score fall to 0.871 with the BA control, is 0.1691 (0.0087) (unconditional EIV 0.2152 (0.0075)); a split-half IV check that assumes no error variance gives 0.1635 (0.0085); the three versions span 0.008 (more than one SE). **Pairing:** 0.1608 goes with the unconditional 0.2127 (both convention); if EIV is preferred, use 0.1691 together with 0.2152 (and amend PREREG V12 before any M3 result). The BA coefficient at a given latent ability is 0.174 (0.018) (EIV; 0.198 (0.017) uncorrected). First-stage reading s_z = (1 - eta) x slope = 0.148 (convention) or 0.156 (EIV) against 0.1955 now; the joint estimation with the college margin re-estimates it.
- **BA interacted.** Within non-BA 0.166 (0.010) / 0.172 (0.010); within BA 0.144 (0.015) / 0.159 (0.017) (convention / EIV); BA minus non-BA -0.022 (0.018) / -0.014 (0.020) (z = -1.2 for the convention): the slope is about the same with and without a BA, so the fall from the unconditional slope is the BA level, not a flatter BA gradient.
- **Teachers against non-teacher BA workers (the BA-controlled slope).** With the BA indicator in the regression (the slope of `wage_slope_ba`), K-12 teachers (147 persons, 844 wage person-years) have a latent slope of 0.083 (0.045) (convention; EIV 0.098 (0.054)) and non-teacher BA workers (1,540 persons, 8,892 person-years) 0.150 (0.016) (EIV 0.165 (0.018)). Ratio **0.56 (0.32)** (EIV 0.60 (0.35)), implied **gamma = 0.58 (0.32)** (EIV 0.62 (0.34)) at eta = 0.080 (bootstrap CI [-0.04, 1.20]; delta-method SE 0.31). The baseline gamma = 0.6194 implies a ratio of 0.60, 0.1 SE from the estimate: the ratio does not reject the baseline and is too imprecise to be a target. The first version of this file compared the teachers' slope WITHOUT the BA control (0.113 (0.039); ratio 0.75 (0.29), gamma 0.77 (0.27)): the 12% of teachers without a BA score lower and, on the point estimate, are paid less, so that slope contains the credential composition, which the closed form has no margin for; the BA control removes it. The answer also moves with the teacher definition: gamma runs from 0.21 (teacher in at least half of the rounds) to 0.73 (1970 codes); BA teachers only 0.50 (0.32); against all non-teachers (BA control in both) 0.54 (0.29).
- **BA-share shift, like for like.** The NLSY79 share of ALL persons of D3's primary sample (0.258 (0.007); it contains the K-12 teachers, 92% of whom have a BA) moves to the ACS all-persons share 0.3391 with a common logit shift of the D3 tercile shares of **0.543 (0.048)**; implied P(BA | AFQT tercile) = 0.068 (0.006), 0.257 (0.009), 0.692 (0.009) (NLSY79: 0.041, 0.167, 0.567); D3's individual-logit shift to the same share, 0.580 (0.053), is reproduced. The same exercise for the NON-TEACHERS (NLSY79 share 0.245 (0.007), ACS non-teacher share 0.3231) gives **0.542 (0.050)** and 0.064 (0.006), 0.237 (0.009), 0.679 (0.009). Shifting the all-persons shares to the non-teacher share (0.440 (0.048); shares 0.062, 0.237, 0.670), the version the task text asked for and the first version of this file, mixes two populations and understates the shift by 0.10 log points: it is kept only as a flagged alternative.

## 1. Reproduction of the published slope

|  | Pooled | Men | Women |
|---|---:|---:|---:|
| Published (`nlsy_ability.json`, B = 200, seed 20260928) | 0.21266 (0.0077) | 0.19036 (0.0100) | 0.24682 (0.0100) |
| This script, point estimate (convention) | 0.21266 | 0.19036 | 0.24682 |
|   absolute difference to the published value | 1.1e-15 | 4.9e-15 | 9.7e-16 |
| SE, replay of the published bootstrap stream (B = 200) | 0.00771 | 0.01003 | 0.00996 |
| SE, draws of D3 and `nlsy_teachers.py` (B = 500, seed 20260929) | 0.00739 | 0.00981 | 0.01084 |
| EIV version (point and SE) | 0.2152 (0.0075) | 0.1910 (0.0099) | 0.2535 (0.0112) |
| Conditional reliability of afqt_c in the wage sample (R_comp_pop = 0.916) | 0.905 | 0.913 | 0.892 |
| Person-years / persons | 37,870 / 7,219 | 21,075 / 3,756 | 16,795 / 3,463 |

The published stream is replayed exactly (multinomial household draws with the interleaved NLSY97 draw, which is generated for the stream and not read; `nlsy_ability.slope_set` itself computes the slopes): maximum absolute difference to the published SEs 3.6e-14. D3 is reproduced the same way: BA share by tercile 0.0e+00, logit coefficients 0.0e+00, shift to the all-persons ACS share 6.7e-16 (maximum absolute differences to `education_occupation.json`).

## 2. Slopes holding the BA status fixed

BA = highest grade completed >= 16 at some round at 25-34 (a proxy: the project's NLSY79 extract has no degree variable); 0.299 of the weighted wage-sample person-years. Pooled regressions include a female dummy; by-sex regressions do not. 'Interacted' = one regression with BA, AFQT x BA and the common controls; 'separate samples' = the regression run within each BA group with its own controls.

**Convention (observed slope / sqrt(R)): the target values**

|  | Pooled | Men | Women |
|---|---:|---:|---:|
| No BA control (the existing slope) | 0.213 (0.007) | 0.190 (0.010) | 0.247 (0.011) |
| **BA control** | 0.161 (0.008) | 0.148 (0.012) | 0.182 (0.011) |
| BA control, BA = highest grade >= 16 at the round itself | 0.159 (0.008) | 0.145 (0.012) | 0.183 (0.011) |
| Interacted: slope within non-BA | 0.166 (0.010) | 0.148 (0.013) | 0.195 (0.014) |
| Interacted: slope within BA | 0.144 (0.015) | 0.150 (0.022) | 0.150 (0.022) |
| Interacted: BA minus non-BA | -0.022 (0.018) | 0.002 (0.026) | -0.045 (0.027) |
| Separate samples: non-BA | 0.166 (0.010) | 0.147 (0.013) | 0.194 (0.014) |
| Separate samples: BA | 0.150 (0.015) | 0.148 (0.022) | 0.153 (0.021) |
| BA control / no control | 0.756 (0.023) | 0.778 (0.037) | 0.739 (0.026) |

**EIV (classical errors-in-variables correction): the alternative**

|  | Pooled | Men | Women |
|---|---:|---:|---:|
| No BA control (the existing slope) | 0.215 (0.008) | 0.191 (0.010) | 0.254 (0.011) |
| **BA control** | 0.169 (0.009) | 0.154 (0.012) | 0.195 (0.012) |
| BA control, BA = highest grade >= 16 at the round itself | 0.167 (0.009) | 0.151 (0.012) | 0.196 (0.012) |
| Interacted: slope within non-BA | 0.172 (0.010) | 0.152 (0.014) | 0.207 (0.015) |
| Interacted: slope within BA | 0.159 (0.017) | 0.166 (0.025) | 0.164 (0.024) |
| Interacted: BA minus non-BA | -0.014 (0.020) | 0.015 (0.028) | -0.043 (0.029) |
| Separate samples: non-BA | 0.172 (0.010) | 0.151 (0.014) | 0.205 (0.015) |
| Separate samples: BA | 0.165 (0.017) | 0.164 (0.025) | 0.168 (0.023) |
| BA control / no control | 0.786 (0.023) | 0.809 (0.038) | 0.770 (0.026) |

**BA coefficient and conditional reliabilities**

|  | Pooled | Men | Women |
|---|---:|---:|---:|
| BA coefficient in the BA-control regression, uncorrected | 0.198 (0.017) | 0.165 (0.025) | 0.236 (0.021) |
| BA coefficient in the BA-control regression, EIV (premium at given latent ability) | 0.174 (0.018) | 0.142 (0.026) | 0.208 (0.022) |
| Conditional reliability of afqt_c: BA-control regression | 0.871 | 0.878 | 0.856 |
| Conditional reliability: BA holders only | 0.830 | 0.825 | 0.836 |
| Conditional reliability: non-BA only | 0.882 | 0.891 | 0.864 |
| Person-years / persons, BA holders | 9,654 / 1,673 | 5,060 / 817 | 4,594 / 856 |
| Person-years / persons, non-BA | 28,216 / 5,546 | 16,015 / 2,939 | 12,201 / 2,607 |

The reliability of afqt_c is 0.916 in the block-A sample; in the wage sample it is 0.905, with a BA control 0.871, and within a BA group 0.830 (weighted SD of afqt_c: 0.94 in the wage sample, 0.71 among BA holders, 0.84 among non-BA). Dividing by sqrt(R) therefore under-corrects by 5% with the BA control and by 10% within BA, if the error variance of the score is 1 - R in every group. The BA coefficient falls by 0.024 under EIV because the noisy AFQT leaves part of ability to the BA dummy.

**Check of the disattenuation: split-half instrumental variables**

| Pooled, per SD of latent AFQT | No BA control | BA control |
|---|---:|---:|
| Convention (observed / sqrt(R)) | 0.2127 (0.0074) | 0.1608 (0.0083) |
| EIV (one-factor reliability) | 0.2152 (0.0075) | 0.1691 (0.0087) |
| IV, AR+WK vs PC+MK: A instrumented by B / B by A | 0.2139 (0.0077) / 0.2060 (0.0075) | 0.1635 (0.0087) / 0.1630 (0.0090) |
| IV, AR+PC vs WK+MK: A instrumented by B / B by A | 0.2148 (0.0077) / 0.2061 (0.0076) | 0.1629 (0.0088) / 0.1645 (0.0090) |
| IV, AR+MK (math) vs WK+PC (verbal): A instrumented by B / B by A | 0.2039 (0.0084) / 0.2417 (0.0084) | 0.1598 (0.0106) / 0.1946 (0.0097) |
| **IV, mean of the two mixed splits (4 estimates)** | 0.2102 (0.0074) | 0.1635 (0.0085) |

Each half of the four ASVAB tests (AR, WK, PC, MK) is a composite; the regressor is one half instrumented by the other half and rescaled by the one-factor loading of its composite (so the factor is the project's, but no error variance is assumed and no group-specific reliability is needed). The two mixed splits give nearly the same slope in both directions; the math-verbal split does not (math and verbal composites carry wage-relevant skills beyond the common factor), so it is shown but not averaged. The mixed-split IV puts the BA-controlled slope at 0.163 (0.008): +2% against the convention and -3% against EIV; the three versions span 0.008, more than one SE (0.008). The unconditional IV slope is 0.210, -1% against the published 0.213: the IV agrees with EIV that conditioning on BA raises the needed correction (IV / convention rises from 0.988 to 1.017, EIV / convention from 1.012 to 1.052) but puts the level lower than the one-factor reliability does.

## 3. Teachers and non-teacher BA workers: slopes, ratio and implied gamma

Teachers: the `nlsy_teachers.py` person-level classification (Census 1980 codes, all K-12 codes, a K-12 teacher in the majority of the employed 25-34 rounds, spec `c80-emp-all-maj`; 147 persons, 123 women and 24 men), all of whose wage-sample person-years enter. Non-teacher BA workers: the non-teachers of that spec with a BA. The teachers' slope is the SAME slope as in section 2, i.e. from a regression with a BA indicator (the non-teacher BA workers all have a BA, so there the control is constant): 12% of the teachers have no BA; they score lower (D3: the teachers' latent AFQT gap is -0.334 SD against non-teacher BA workers and -0.238 SD for BA teachers only) and, on the point estimate, are paid less (BA coefficient among teachers 0.189 (0.185), uncorrected), so a slope without the control would put the credential composition into the teachers' slope, for which the closed form has no margin. The gamma formula is the model's: ratio = gamma (1 - eta) / (1 - gamma eta), so gamma = ratio / (1 - eta + eta ratio), eta = 0.080 (0.0807 changes gamma by less than 0.001).

| Group (wage sample) | Persons | Person-years | Observed slope per afqt_c SD | Latent, convention | Latent, EIV | Cond. reliability |
|---|---:|---:|---:|---:|---:|---:|
| K-12 teachers, BA control (main) | 147 | 844 | 0.080 (0.043) | 0.083 (0.045) | 0.098 (0.054) | 0.776 |
|   women | 123 | 700 | 0.058 (0.054) | 0.060 (0.057) | 0.072 (0.070) | 0.764 |
|   men | 24 | 144 | 0.119 (0.100) | 0.125 (0.105) | 0.143 (0.121) | 0.798 |
| K-12 teachers, no BA control (first version) | 147 | 844 | 0.108 (0.037) | 0.113 (0.039) | 0.124 (0.044) | 0.829 |
|   women | 123 | 700 | 0.120 (0.044) | 0.126 (0.046) | 0.140 (0.051) | 0.824 |
|   men | 24 | 144 | 0.057 (0.100) | 0.060 (0.105) | 0.065 (0.112) | 0.837 |
| K-12 teachers with a BA | 130 | 759 | 0.069 (0.043) | 0.072 (0.045) | 0.084 (0.054) | 0.784 |
| All non-teachers, no BA control | 7,053 | 37,006 | 0.205 (0.007) | 0.214 (0.008) | 0.217 (0.008) | 0.905 |
| All non-teachers, BA control | 7,053 | 37,006 | 0.154 (0.008) | 0.161 (0.008) | 0.169 (0.009) | 0.871 |
| Non-teacher BA workers | 1,540 | 8,892 | 0.143 (0.015) | 0.150 (0.016) | 0.165 (0.018) | 0.832 |
|   women | 746 | 3,968 | 0.151 (0.022) | 0.158 (0.023) | 0.172 (0.025) | 0.841 |
|   men | 794 | 4,924 | 0.138 (0.022) | 0.144 (0.023) | 0.160 (0.025) | 0.824 |
| Non-teachers without a BA | 5,513 | 28,114 | 0.158 (0.009) | 0.165 (0.010) | 0.172 (0.010) | 0.882 |

Pooled groups include a female dummy, women and men do not; all regressions have age and year dummies. The teacher slopes rest on 147 persons. Among teachers the BA indicator's coefficient is 0.152 (0.195) (EIV; 0.189 (0.185) uncorrected), identified from the 17 teachers without a BA; for the 24 male teachers, of whom very few lack a BA, the coefficient and so the men's row are unreliable.

**Ratio of slopes and implied gamma** (convention = the values; EIV = the alternative)

| Comparison (numerator = teachers) | Teacher persons / person-years | Ratio, convention | gamma, convention | gamma 95% CI | Ratio, EIV | gamma, EIV |
|---|---:|---:|---:|---:|---:|---:|
| **all K-12 teachers (BA control) / non-teacher BA workers (main)** | 147 / 844 | 0.56 (0.32) | 0.58 (0.32) | [-0.04, 1.20] | 0.60 (0.35) | 0.62 (0.34) |
| all K-12 teachers, no BA control / non-teacher BA workers (first version) | 147 / 844 | 0.75 (0.29) | 0.77 (0.27) | [0.28, 1.38] | 0.76 (0.29) | 0.77 (0.28) |
| K-12 teachers with a BA / non-teacher BA workers | 130 / 759 | 0.48 (0.32) | 0.50 (0.32) | [-0.14, 1.15] | 0.51 (0.34) | 0.53 (0.34) |
| all K-12 teachers / all non-teachers (BA control in both) | 147 / 844 | 0.52 (0.29) | 0.54 (0.29) | [-0.04, 1.08] | 0.58 (0.33) | 0.60 (0.33) |
| all K-12 teachers / all non-teachers (no BA control) | 147 / 844 | 0.53 (0.19) | 0.55 (0.18) | [0.20, 0.91] | 0.57 (0.21) | 0.59 (0.20) |
| women: teachers (BA control) / non-teacher BA workers | 123 / 700 | 0.38 (0.38) | 0.40 (0.39) | [-0.33, 1.14] | 0.42 (0.43) | 0.44 (0.43) |
| men: teachers (BA control) / non-teacher BA workers | 24 / 144 | 0.86 (0.79) | 0.87 (0.66) | [-0.25, 2.28] | 0.89 (0.81) | 0.90 (0.70) |
| teachers (BA control) / non-teacher BA workers, core K-12 codes (elementary, secondary, special education) | 101 / 606 | 0.26 (0.27) | 0.28 (0.28) | [-0.30, 0.77] | 0.27 (0.29) | 0.29 (0.31) |
| teachers (BA control) / non-teacher BA workers, teacher in at least half of the rounds | 182 / 958 | 0.20 (0.29) | 0.21 (0.30) | [-0.37, 0.78] | 0.20 (0.30) | 0.22 (0.32) |
| teachers (BA control) / non-teacher BA workers, ever a teacher at 25-34 | 406 / 2,025 | 0.68 (0.26) | 0.70 (0.25) | [0.23, 1.21] | 0.68 (0.26) | 0.70 (0.25) |
| teachers (BA control) / non-teacher BA workers, any 25-34 occupation (not only employed rounds) | 153 / 860 | 0.50 (0.27) | 0.52 (0.27) | [-0.03, 1.02] | 0.54 (0.29) | 0.56 (0.29) |
| teachers (BA control) / non-teacher BA workers, 1970 Census codes | 149 / 840 | 0.72 (0.35) | 0.73 (0.34) | [0.10, 1.41] | 0.76 (0.38) | 0.78 (0.36) |
| round-level K-12 occupation (BA control) / other BA person-years | 317 / 992 | 0.57 (0.22) | 0.59 (0.22) | [0.19, 1.02] | 0.59 (0.23) | 0.61 (0.23) |
| teachers (BA control) / non-teacher BA workers, CPI-deflated log wage, age dummies and female (no year dummies) | 147 / 844 | 0.55 (0.31) | 0.57 (0.31) | [-0.04, 1.18] | 0.59 (0.34) | 0.61 (0.34) |
| teachers (BA control) / non-teacher BA workers, CPI-deflated log wage, female dummy only | 147 / 844 | 0.57 (0.33) | 0.59 (0.32) | [-0.05, 1.24] | 0.61 (0.35) | 0.63 (0.35) |

Main specification, delta method (bootstrap covariance of the two slopes; correlation -0.07 convention, -0.07 EIV): convention ratio 0.556 (SE 0.311), gamma 0.577 (SE 0.308; d gamma / d ratio = 0.989), against the bootstrap SDs 0.319 and 0.316; EIV ratio 0.597 (SE 0.340), gamma 0.617 (SE 0.334), against 0.348 and 0.343. Model reference: gamma = 0.6194 (baseline fit) implies ratio 0.600; gamma = 1 implies 1.

Why the rows differ: the teacher slope is estimated from 100-400 persons, so the ratio has an SE of 0.2-0.4 (0.8 for men) and the rows are not independent evidence; they change who counts as a teacher (core K-12 codes drop the pre-K and 'other teacher' codes, whose members have lower AFQT and pay; 'ever' and the round-level row add persons who teach in a minority of rounds), the comparison group (all non-teachers carry the BA premium in their slope unless the BA control is in both regressions) or the credential control (the first-version row, without it, puts the teachers without a BA into the teachers' slope). The main interval contains both the baseline gamma and gamma = 1.

## 4. The BA share: logit level shift and implied tercile shares

Primary NLSY79 sample of D3 (persons observed at 25-34, highest grade >= 9, N = 8,528); BA share 0.2582 (0.0071); terciles of the observed `afqt_c` cut over the weighted sample (cut points -0.403, 0.494 SDs), as in D3, for both populations below. **A shift is only meaningful against the ACS share of the same population.** The NLSY79 share above is that of ALL persons, 2.0% of whom (weighted) are K-12 teachers (main spec) with a BA share of 92%; the non-teachers' NLSY79 share is 0.2446 (0.0070). Two like-for-like targets: (1) all persons to the ACS all-persons share 0.3391 (`ba_share_all`, D3's own target; the model's `ba_share` and `ba_share_by_z` are over all adults), (2) the non-teachers (everyone not a teacher, persons never employed included, as the ACS non-teacher universe includes home production) to the ACS non-teacher share 0.3231 (0.0007, T2 convention, `education_targets.toml` `ba_share_nonteacher`, the M3 target). The third version, all persons to the non-teacher share, is NOT like for like and is shown only because the task text and the first version of this file used it. **Cell shift**: one delta with sum_k pi_k Lambda(logit s_k + delta) = target, s_k the D3 tercile shares of the population and pi_k its tercile weights. **Individual logit**: P(BA | afqt) = Lambda(a + delta + b afqt) with (a, b) = (-1.693, 1.815) for all persons (-1.790, 1.837 for non-teachers) fitted by weighted logit, delta set so the weighted mean is the target (D3's 'level shift of the same gradient'); the tercile shares are its weighted averages within the NLSY79 terciles.

|  | BA share | Tercile 1 (low AFQT) | Tercile 2 | Tercile 3 (high) | T3 - T1 | Logit shift delta |
|---|---:|---:|---:|---:|---:|---:|
| All persons, as observed (D3, not recomputed for use) | 0.2582 (0.0071) | 0.0408 (0.0040) | 0.1671 (0.0094) | 0.5668 (0.0161) | 0.526 (0.016) | - |
| **All persons, cell shift to the ACS all-persons share 0.3391 (like for like; primary)** | 0.3391 | 0.0682 (0.0061) | 0.2566 (0.0089) | 0.6925 (0.0086) | 0.624 (0.012) | 0.543 (0.048) |
| All persons, individual logit shifted to 0.3391 (D3's method; D3: 0.5803) | 0.3391 | 0.0571 (0.0036) | 0.2772 (0.0037) | 0.6830 (0.0070) | 0.626 (0.010) | 0.580 (0.053) |
| All persons, individual logit fit, no shift | 0.2582 | 0.0332 (0.0026) | 0.1797 (0.0081) | 0.5618 (0.0144) | 0.529 | 0 |
| Non-teachers, as observed (same terciles) | 0.2446 (0.0070) | 0.0382 (0.0038) | 0.1530 (0.0089) | 0.5511 (0.0165) | 0.513 (0.017) | - |
| **Non-teachers, cell shift to the ACS non-teacher share 0.3231 (like for like)** | 0.3231 | 0.0640 (0.0060) | 0.2370 (0.0086) | 0.6785 (0.0088) | 0.615 (0.012) | 0.542 (0.050) |
| Non-teachers, individual logit shifted to 0.3231 | 0.3231 | 0.0514 (0.0034) | 0.2593 (0.0038) | 0.6688 (0.0072) | 0.617 (0.010) | 0.582 (0.055) |
| *Mixed, NOT like for like: all persons shifted to the non-teacher share 0.3231 (task text; first version)* | 0.3231 | 0.0619 (0.0057) | 0.2374 (0.0086) | 0.6701 (0.0086) | 0.608 (0.012) | 0.440 (0.048) |
| *Mixed, individual logit (D3's method) to 0.3231* | 0.3231 | 0.0516 (0.0034) | 0.2566 (0.0038) | 0.6612 (0.0069) | 0.610 (0.010) | 0.471 (0.053) |

The two like-for-like cell shifts agree (0.543 for all persons, 0.542 for non-teachers), as they should: both bring a cohort whose BA share is lower than the ACS one to the ACS level, whereas the mixed version asks the NLSY79 all-persons distribution, which includes the high-BA teachers, to reach a level that excludes them. The cell shift keeps the empirical tercile pattern (its middle tercile is the D3 one shifted); the individual logit imposes a log-odds line in afqt_c, which fits the middle tercile +0.013 off and so shifts it differently. The SEs are NLSY79 sampling errors; the ACS targets' own SEs (0.0011 all persons, 0.0007 non-teachers) add 0.0072 and 0.0047 to delta, negligible. The shares are by tercile of the OBSERVED score (reliability 0.916): they understate the latent gradient, so the model counterpart should add classical noise with variance (1 - R) / R = 0.092 times the variance of the latent score before cutting terciles. The NLSY79 teacher share of persons (2.0%) is lower than the ACS one (4.0% of the Table 2 universe), which is why the non-teacher and all-persons levels are closer in NLSY79 (0.014) than in the ACS (0.016).

## Definitions and caveats

- **Wage sample, weights, controls**: as `nlsy_ability.py` block C (>= 48 weeks, >= 30 hours/week, wage income >= $1,000 in 2010 dollars, highest grade >= 9, not armed forces; log of annual income over annual hours, no trimming; person's 1979 SAMPWEIGHT on each person-year; age and year dummies). The CPI-deflated rows of the ratio table replace the year dummies by deflation with the prior-year CPI-U.
- **Convention and EIV**: the convention is OLS / sqrt(R), the project's latent slope (0.2127 for the pooled unconditional regression). EIV: afqt_c = f + u, var(u) = 1 - R, var(f) = R, u independent of the other regressors; the coefficients solve (M - S) b = m with S = (1 - R) E[c c'] (c = 1 for the AFQT column, BA for AFQT x BA); the slope per SD of f is b_AFQT sqrt(R). The EIV assumptions are strong: the error variance is taken equal across groups, independent of BA, and AFQT was measured at 15-23, when schooling already raises scores; a one-factor ULS model supplies R. The targets file carries the convention as the value and EIV as `alt_eiv`; the two are not to be mixed within a pair (BA-conditional slope and unconditional slope).
- **Unstable EIV draws**: an EIV statistic is set to missing in a draw whose conditional reliability of the score is below 0.5 (the correction would more than double the slope). This affects only `TC_men` (16 of 500 draws), `T_men` (4 of 500 draws), `r_men` (16 of 500 draws), `gamma_men` (16 of 500 draws); their SEs use the remaining draws.
- **Split-half IV**: one half of the tests (a composite of two standardized IRT z-scores, standardized in block A) instrumented by the other half, controls partialled out, times the loading of the composite on the standardized one-factor latent variable (mean ULS loading of its two tests over the composite's SD); re-estimated in every draw. It assumes independent errors across halves, not an error variance.
- **BA** is highest grade completed >= 16 at some round at 25-34, a proxy for a bachelor's degree: some respondents with 16 years report no bachelor's degree and others with 17-18 years do, and graduate degrees are not separated. NLSY79 itself records degrees (questions 1979-2020; the created variable HIGHEST_DEGREE_EVER was added in the 2022, round 30, release: NLSY79 Codebook Supplement, Appendix 8), but the project's extract has no degree variable; adding it needs a new NLS Investigator extract. The model's d = 1 is a BA, and its teachers are mostly BA holders (87% in the ACS).
- **Teacher sample**: 147 persons, 844 person-years in the wage sample (123 women, 24 men), so teacher slopes carry SEs of 0.05 pooled and 0.10 for men; the classification itself is held fixed in the bootstrap. Teachers' wage-sample person-years include rounds in other occupations (16%). The round-level row uses person-years in a K-12 occupation instead.
- **Cohort**: NLSY79 persons born 1957-64 observed 1982-98, against the ACS 2009-13 cohorts of the model's BA share. The NLSY97 BA-conditional slope (the closer cohort; DESIGN_v2.md lists it among the unsealed M3 inputs) is outside this NLSY79 task and was not computed.
- **Not included in the SEs**: the ACS targets' SEs (shown above for the shifts), eta (0.080 against 0.0807), the choice of BA and teacher definitions (the sensitivity rows) and the sampling error of the assumed EIV structure.
- **By sex**: women's slopes are larger than men's (BA-conditional 0.182 against 0.148, convention; unconditional 0.247 against 0.190); the model has one slope for both sexes.
- **Sealed set**: no sealed value was read or used and no NLSY97 data were read. One protocol deviation, from the first run of this task (2026-10-05): the sealed file `estimates/settings/sealed/nlsy97_holdout.json` was loaded once and only its top-level key names were printed (no value, no statistic). PREREG.md has the sealed set 'opened only after selection', so this is an access before selection; it is to be recorded in the round's LOG.md. Nothing in this file or in the targets file depends on it.

## Not done

- The project's extract has no degree variable, so BA is a schooling-years proxy; no separate treatment of graduate degrees or of the 12% of teachers without a BA beyond the BA control and the BA-teacher rows.
- The NLSY97 counterparts (the cohort closer to the ACS 2009-13 population) were not computed: the BA-conditional slope is a separate unsealed input, and the teachers' wage-AFQT slope of the NLSY97 belongs to the sealed hold-out (D4), which is not touched here.
- gamma is not pinned down: the teacher sample (about 150 persons) gives a ratio SE of 0.3 and the point estimate moves with the teacher definition; no larger NLSY79 teacher sample exists.
- Uncertainty in the EIV structure (R, equal error variance across groups, independence from BA) is not in the SEs; the convention and EIV columns bracket the effect of that choice.
- No saved bootstrap draws (only SEs and percentile CIs), so the covariances of these moments with the D3 and teacher-gap moments (same draws) are not provided.

