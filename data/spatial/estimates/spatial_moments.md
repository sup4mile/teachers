# Spatial moments for `spatial_continuous.jl`

Generated 2026-08-31T11:17:39 by `data_estimate.py`. Base year **2018**, partition scheme **all**, labor market **czone (Dorn 1990 commuting zones)**; a location is a group of school districts inside one labor market, and location 2 is the advantaged group.


## Inputs

| table | status | columns |
|---|---|---|
| `acs_school_districts` | ok | 15 |
| `czone_crosswalk` | ok | 2 |
| `edge_cwift` | ok | 15 |
| `edge_geocode_lea` | ok | 65 |
| `seda_*` | ok | 98 tables |
| `urban_ccd_directory` | ok | 69 |
| `urban_ccd_enrollment` | ok | 7 |
| `urban_ccd_finance` | ok | 163 |
| `urban_edfacts_assessments` | ok | 23 |
| `urban_edfacts_grad_rates` | ok | 15 |
| `urban_saipe` | ok | 10 |

`within-CZ var share` is the fraction of the 2018 cross-district variance in that variable that survives commuting-zone fixed effects. It is the descriptive case for the whole design: where it is high, a within-CZ two-location model is where the variation is, and a cross-metro model (Option 3 in the calibration notes) would be explaining the smaller half.


Panel: 85497 district-years, 16477 districts, 738 CZs, 621 CZ-years surviving the split.


## Estimates from the data

| key | model object | value | ratio 2/1 | within-CZ var share | n |
|---|---|---|---|---|---|
| `gap_pupils_locale` | M_2/M_1 | 0.3764 | 1.457 | 0.58 | 188 |
| `gap_teachers_locale` | H̃_T,2/H̃_T,1 (before the β weighting) | 0.3632 | 1.438 | 0.59 | 188 |
| `gap_teachers_pp_locale` | H̃_T,l / M_l | -0.01216 | 0.988 | 0.18 | 188 |
| `gap_str_ratio_locale` | class size N(h) | 0.01216 | 1.012 | 0.18 | 188 |
| `gap_salary_per_teacher_locale` | κ_2/κ_1 (undeflated) | -0.01705 | 0.983 | 0.13 | 177 |
| `gap_salary_real_locale` | κ_2/κ_1 | 0.01623 | 1.016 | 0.21 | 176 |
| `gap_exp_pp_locale` | W_l / M_l | -0.09564 | 0.909 | 0.20 | 188 |
| `gap_rev_pp_local_locale` | t_l I_l / M_l | -0.07082 | 0.932 | 0.44 | 188 |
| `gap_rev_pp_state_locale` | G_l (notes §5.2) | -0.0003101 | 1.000 | 0.44 | 188 |
| `gap_rev_pp_total_locale` | W_l / M_l | -0.1211 | 0.886 | 0.25 | 188 |
| `gap_proptax_pp_locale` | t_l | -0.1041 | 0.901 | 0.38 | 166 |
| `gap_t_eff_locale` | t_l | 0.0006032 | — | 0.38 | 188 |
| `gap_share_local_locale` | the pure local budget of §10 | 0.01055 | — | 0.52 | 188 |
| `gap_share_state_locale` | G_l | 0.01963 | — | 0.49 | 188 |
| `gap_pov_rate_5_17_locale` | ability/income sorting | -0.08878 | — | 0.70 | 188 |
| `gap_learn_rate_locale` | Q_2/Q_1 | 0.003072 | — | 0.67 | 187 |
| `gap_score_mean_locale` | Q level (composition-laden) | 0.2887 | — | 0.82 | 187 |
| `gap_seda_ses_locale` | ability/income sorting | 0.8299 | — | 0.78 | 187 |
| `gap_seda_lninc_locale` | I_2/I_1 (the income base) | 0.2708 | — | 0.67 | 187 |
| `trend_gap_salary_real_locale` | comparative statics across steady states | 0.02697 | — | — | — |
| `trend_gap_rev_pp_local_locale` | comparative statics across steady states | -0.04694 | — | — | — |
| `trend_gap_exp_pp_locale` | comparative statics across steady states | 0.03924 | — | — | — |
| `trend_gap_pupils_locale` | comparative statics across steady states | 0.08806 | — | — | — |
| `trend_gap_share_local_locale` | comparative statics across steady states | -0.03157 | — | — | — |
| `kappa_quality_gradient` | κ_l alongside Q_l | 0.04033 (0.0241) | — | — | 2721 |
| `revenue_to_salary_passthrough` | the balanced budget t_l I_l = W_l | 0.05598 (0.0133) | — | — | 2885 |
| `local_revenue_share` | the pure local budget of §10 vs. the G_l extension | 0.4726 | — | — | 3249 |
| `gap_pupils_revenue` | M_2/M_1 | -0.02474 | 0.976 | 0.48 | 621 |
| `gap_teachers_revenue` | H̃_T,2/H̃_T,1 (before the β weighting) | 0.01313 | 1.013 | 0.49 | 621 |
| `gap_teachers_pp_revenue` | H̃_T,l / M_l | 0.03661 | 1.037 | 0.20 | 621 |
| `gap_str_ratio_revenue` | class size N(h) | -0.03661 | 0.964 | 0.20 | 621 |
| `gap_salary_per_teacher_revenue` | κ_2/κ_1 (undeflated) | 0.04258 | 1.043 | 0.16 | 594 |
| `gap_salary_real_revenue` | κ_2/κ_1 | 0.01543 | 1.016 | 0.26 | 593 |
| `gap_exp_pp_revenue` | W_l / M_l | 0.07589 | 1.079 | 0.20 | 621 |
| `gap_rev_pp_local_revenue` | t_l I_l / M_l | 0.5636 | 1.757 | 0.41 | 621 |
| `gap_rev_pp_state_revenue` | G_l (notes §5.2) | -0.3193 | 0.727 | 0.44 | 621 |
| `gap_rev_pp_total_revenue` | W_l / M_l | 0.1143 | 1.121 | 0.24 | 621 |
| `gap_proptax_pp_revenue` | t_l | 0.6012 | 1.824 | 0.37 | 556 |
| `gap_t_eff_revenue` | t_l | 0.005521 | — | 0.45 | 620 |
| `gap_share_local_revenue` | the pure local budget of §10 | 0.1926 | — | 0.50 | 621 |
| `gap_share_state_revenue` | G_l | -0.1711 | — | 0.48 | 621 |
| `gap_pov_rate_5_17_revenue` | ability/income sorting | -0.03599 | — | 0.61 | 620 |
| `gap_learn_rate_revenue` | Q_2/Q_1 | 0.007128 | — | 0.65 | 613 |
| `gap_score_mean_revenue` | Q level (composition-laden) | 0.1994 | — | 0.77 | 615 |
| `gap_seda_ses_revenue` | ability/income sorting | 0.5735 | — | 0.69 | 615 |
| `gap_seda_lninc_revenue` | I_2/I_1 (the income base) | 0.1573 | — | 0.55 | 615 |
| `trend_gap_salary_real_revenue` | comparative statics across steady states | -0.007311 | — | — | — |
| `trend_gap_rev_pp_local_revenue` | comparative statics across steady states | 0.01498 | — | — | — |
| `trend_gap_exp_pp_revenue` | comparative statics across steady states | -0.03062 | — | — | — |
| `trend_gap_pupils_revenue` | comparative statics across steady states | -0.05651 | — | — | — |
| `trend_gap_share_local_revenue` | comparative statics across steady states | 0.009807 | — | — | — |
| `gap_pupils_ses` | M_2/M_1 | -0.001287 | 0.999 | 0.47 | 619 |
| `gap_teachers_ses` | H̃_T,2/H̃_T,1 (before the β weighting) | -0.00613 | 0.994 | 0.49 | 618 |
| `gap_teachers_pp_ses` | H̃_T,l / M_l | -0.003611 | 0.996 | 0.20 | 618 |
| `gap_str_ratio_ses` | class size N(h) | 0.003611 | 1.004 | 0.20 | 618 |
| `gap_salary_per_teacher_ses` | κ_2/κ_1 (undeflated) | 0.02959 | 1.030 | 0.16 | 592 |
| `gap_salary_real_ses` | κ_2/κ_1 | 0.01875 | 1.019 | 0.26 | 591 |
| `gap_exp_pp_ses` | W_l / M_l | -0.05091 | 0.950 | 0.19 | 619 |
| `gap_rev_pp_local_ses` | t_l I_l / M_l | 0.2321 | 1.261 | 0.41 | 619 |
| `gap_rev_pp_state_ses` | G_l (notes §5.2) | -0.2573 | 0.773 | 0.43 | 619 |
| `gap_rev_pp_total_ses` | W_l / M_l | -0.0455 | 0.956 | 0.23 | 619 |
| `gap_proptax_pp_ses` | t_l | 0.2326 | 1.262 | 0.37 | 553 |
| `gap_t_eff_ses` | t_l | 0.001293 | — | 0.44 | 618 |
| `gap_share_local_ses` | the pure local budget of §10 | 0.1212 | — | 0.50 | 619 |
| `gap_share_state_ses` | G_l | -0.08085 | — | 0.48 | 619 |
| `gap_pov_rate_5_17_ses` | ability/income sorting | -0.1018 | — | 0.61 | 619 |
| `gap_learn_rate_ses` | Q_2/Q_1 | 0.01455 | — | 0.65 | 617 |
| `gap_score_mean_ses` | Q level (composition-laden) | 0.4047 | — | 0.77 | 619 |
| `gap_seda_ses_ses` | ability/income sorting | 1.135 | — | 0.69 | 619 |
| `gap_seda_lninc_ses` | I_2/I_1 (the income base) | 0.3542 | — | 0.55 | 619 |
| `trend_gap_salary_real_ses` | comparative statics across steady states | 0.01081 | — | — | — |
| `trend_gap_rev_pp_local_ses` | comparative statics across steady states | 0.005593 | — | — | — |
| `trend_gap_exp_pp_ses` | comparative statics across steady states | -0.02365 | — | — | — |
| `trend_gap_pupils_ses` | comparative statics across steady states | 0.03478 | — | — | — |
| `trend_gap_share_local_ses` | comparative statics across steady states | 0.01482 | — | — | — |

- **`gap_pupils_locale`** — within-CZ gap, enrollment, group total (location 2 - location 1). M_2/M_1, the enrollment split. B_l is the residual backed out to match it, so this is a target, never a check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_teachers_locale`** — within-CZ gap, FTE teachers, group total (location 2 - location 1). The teacher stock by location. Together with the enrollment gap it says whether the good location's advantage is more teachers or better ones -- the model's answer is 'better', through β.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_teachers_pp_locale`** — within-CZ gap, FTE teachers per pupil (location 2 - location 1). This is H̃_T,l/M_l up to the β-weighting of teacher human capital, so it maps into Q_l = (2H̃_T,l/M_l)^σ directly. If it is near zero while the score gap is not, the quality gap is composition, not class size.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_str_ratio_locale`** — within-CZ gap, student-teacher ratio (location 2 - location 1). Class size. Prop. 1 has better teachers take LARGER classes, which is not what districts do -- check the mapping (notes §5.3) before targeting this.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_salary_per_teacher_locale`** — within-CZ gap, teacher salary per FTE, nominal (location 2 - location 1). Compare with the deflated version: the wedge between the two is how much of the raw salary gap is cost of living rather than pay.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_salary_real_locale`** — within-CZ gap, teacher salary per FTE, CWIFT-deflated (location 2 - location 1). The direct target for κ_2/κ_1; CWIFT-deflated, so it is a teaching-wage gap rather than a local price level.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_exp_pp_locale`** — within-CZ gap, current expenditure per pupil (location 2 - location 1). The resource gap the budget constraint has to deliver.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_rev_pp_local_locale`** — within-CZ gap, local revenue per pupil (location 2 - location 1). The local half of the budget; with the state half it pins how far pure local finance is from the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_rev_pp_state_locale`** — within-CZ gap, state revenue per pupil (location 2 - location 1). Motivates the exogenous transfer G_l of notes §5.2 -- state aid is compensatory, so expect this gap to be NEGATIVE even where the local gap is positive.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_rev_pp_total_locale`** — within-CZ gap, total revenue per pupil (location 2 - location 1). Local plus state plus federal. If this gap has the opposite sign to the local one, redistribution is doing more than local finance, and a pure local budget gets the sign of the resource gap wrong.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_proptax_pp_locale`** — within-CZ gap, property-tax revenue per pupil (location 2 - location 1). The property-tax base per pupil, the closest thing to t_l I_l in the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_t_eff_locale`** — within-CZ gap, local revenue / local income base (location 2 - location 1). The effective local tax rate. Keep it endogenous from the budget and use this as an untargeted check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_share_local_locale`** — within-CZ gap, local share of revenue (location 2 - location 1). If this is far below 1 the pure local-finance budget in §10 is counterfactual -- see G_l.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_share_state_locale`** — within-CZ gap, state share of revenue (location 2 - location 1). The size of the transfer the model is currently missing.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_pov_rate_5_17_locale`** — within-CZ gap, child poverty rate 5-17 (location 2 - location 1). The observable counterpart of sorting on z: how much richer the good location's families are.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_learn_rate_locale`** — within-CZ gap, SEDA learning rate (location 2 - location 1). The closest public analogue to Q_2/Q_1; growth rather than levels, so it nets out the SES composition the model is supposed to generate endogenously.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_score_mean_locale`** — within-CZ gap, SEDA mean score (location 2 - location 1). NOT a Q target -- it mixes school quality with who enrolls. Useful as the untargeted check that the model's sorting reproduces a level gap it did not fit.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_seda_ses_locale`** — within-CZ gap, SEDA composite SES (location 2 - location 1). Same as the poverty gap, on SEDA's composite. Reported in SEDA index units, whose cross-district SD is about 0.9 -- divide by that before setting it beside the model's mean-ability gap (`idx_work_signed`, ≈ +0.08 at the baseline), which is in units of z. The sign is the part that compares directly.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`gap_seda_lninc_locale`** — within-CZ gap, log median household income (location 2 - location 1). The income base I_l, in logs. This is what t_l is levied on, so together with the local-revenue gap it says whether the advantaged location raises more because it is richer or because it taxes itself harder.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (locale split), enrollment-weighted across 188 CZs
- **`trend_gap_salary_real_locale`** — change in the within-CZ salary_real gap, 2014->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, locale split, CZ-years
- **`trend_gap_rev_pp_local_locale`** — change in the within-CZ rev_pp_local gap, 2006->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, locale split, CZ-years
- **`trend_gap_exp_pp_locale`** — change in the within-CZ exp_pp gap, 2006->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, locale split, CZ-years
- **`trend_gap_pupils_locale`** — change in the within-CZ pupils gap, 2006->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, locale split, CZ-years
- **`trend_gap_share_local_locale`** — change in the within-CZ share_local gap, 2006->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, locale split, CZ-years
- **`kappa_quality_gradient`** — log real teacher salary on district learning rate, within CZ. The gradient is NOT distinguishable from zero (t = +1.67), so this is an upper bound on how much the cross-section knows. Sign tells you whether κ_l and Q_l co-move (agglomeration of advantage) or offset (compensating differentials). The model sets both exogenously; if they co-move strongly in the data, κ_2 > κ_1 is doing the sorting work the paper attributes to the spillover, and the two need separating.  
  *Source:* F-33 teacher salaries / CCD FTE, CWIFT-deflated, x SEDA, 2018, CZ fixed effects
- **`revenue_to_salary_passthrough`** — elasticity of teacher salary to local revenue per pupil. The elasticity is significant (t = +4.20). The model spends the whole local budget on teachers, i.e. an elasticity of 1. Anything well below 1 says the budget also buys class-size reduction, facilities and non-teaching staff, and that the mapping from t_l to κ_l needs a wedge.  
  *Source:* F-33 2018, district within CZ, CZ fixed effects
- **`local_revenue_share`** — local share of district revenue (enrollment-weighted mean). Benchmark: NCES reports 0.44 nationally. A pure local budget describes no year in the sample, which is the case for adding the exogenous transfer G_l (notes §5.2) and makes finance equalization the natural counterfactual.  
  *Source:* F-33 2018, all sample districts, enrollment-weighted
- **`gap_pupils_revenue`** — within-CZ gap, enrollment, group total (location 2 - location 1). M_2/M_1, the enrollment split. B_l is the residual backed out to match it, so this is a target, never a check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_teachers_revenue`** — within-CZ gap, FTE teachers, group total (location 2 - location 1). The teacher stock by location. Together with the enrollment gap it says whether the good location's advantage is more teachers or better ones -- the model's answer is 'better', through β.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_teachers_pp_revenue`** — within-CZ gap, FTE teachers per pupil (location 2 - location 1). This is H̃_T,l/M_l up to the β-weighting of teacher human capital, so it maps into Q_l = (2H̃_T,l/M_l)^σ directly. If it is near zero while the score gap is not, the quality gap is composition, not class size.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_str_ratio_revenue`** — within-CZ gap, student-teacher ratio (location 2 - location 1). Class size. Prop. 1 has better teachers take LARGER classes, which is not what districts do -- check the mapping (notes §5.3) before targeting this.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_salary_per_teacher_revenue`** — within-CZ gap, teacher salary per FTE, nominal (location 2 - location 1). Compare with the deflated version: the wedge between the two is how much of the raw salary gap is cost of living rather than pay.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_salary_real_revenue`** — within-CZ gap, teacher salary per FTE, CWIFT-deflated (location 2 - location 1). The direct target for κ_2/κ_1; CWIFT-deflated, so it is a teaching-wage gap rather than a local price level.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_exp_pp_revenue`** — within-CZ gap, current expenditure per pupil (location 2 - location 1). The resource gap the budget constraint has to deliver.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_rev_pp_local_revenue`** — within-CZ gap, local revenue per pupil (location 2 - location 1). The local half of the budget; with the state half it pins how far pure local finance is from the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_rev_pp_state_revenue`** — within-CZ gap, state revenue per pupil (location 2 - location 1). Motivates the exogenous transfer G_l of notes §5.2 -- state aid is compensatory, so expect this gap to be NEGATIVE even where the local gap is positive.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_rev_pp_total_revenue`** — within-CZ gap, total revenue per pupil (location 2 - location 1). Local plus state plus federal. If this gap has the opposite sign to the local one, redistribution is doing more than local finance, and a pure local budget gets the sign of the resource gap wrong.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_proptax_pp_revenue`** — within-CZ gap, property-tax revenue per pupil (location 2 - location 1). The property-tax base per pupil, the closest thing to t_l I_l in the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_t_eff_revenue`** — within-CZ gap, local revenue / local income base (location 2 - location 1). The effective local tax rate. Keep it endogenous from the budget and use this as an untargeted check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_share_local_revenue`** — within-CZ gap, local share of revenue (location 2 - location 1). If this is far below 1 the pure local-finance budget in §10 is counterfactual -- see G_l.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_share_state_revenue`** — within-CZ gap, state share of revenue (location 2 - location 1). The size of the transfer the model is currently missing.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_pov_rate_5_17_revenue`** — within-CZ gap, child poverty rate 5-17 (location 2 - location 1). The observable counterpart of sorting on z: how much richer the good location's families are.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_learn_rate_revenue`** — within-CZ gap, SEDA learning rate (location 2 - location 1). The closest public analogue to Q_2/Q_1; growth rather than levels, so it nets out the SES composition the model is supposed to generate endogenously.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_score_mean_revenue`** — within-CZ gap, SEDA mean score (location 2 - location 1). NOT a Q target -- it mixes school quality with who enrolls. Useful as the untargeted check that the model's sorting reproduces a level gap it did not fit.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_seda_ses_revenue`** — within-CZ gap, SEDA composite SES (location 2 - location 1). Same as the poverty gap, on SEDA's composite. Reported in SEDA index units, whose cross-district SD is about 0.9 -- divide by that before setting it beside the model's mean-ability gap (`idx_work_signed`, ≈ +0.08 at the baseline), which is in units of z. The sign is the part that compares directly.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`gap_seda_lninc_revenue`** — within-CZ gap, log median household income (location 2 - location 1). The income base I_l, in logs. This is what t_l is levied on, so together with the local-revenue gap it says whether the advantaged location raises more because it is richer or because it taxes itself harder.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (revenue split), enrollment-weighted across 621 CZs
- **`trend_gap_salary_real_revenue`** — change in the within-CZ salary_real gap, 2014->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, revenue split, CZ-years
- **`trend_gap_rev_pp_local_revenue`** — change in the within-CZ rev_pp_local gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, revenue split, CZ-years
- **`trend_gap_exp_pp_revenue`** — change in the within-CZ exp_pp gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, revenue split, CZ-years
- **`trend_gap_pupils_revenue`** — change in the within-CZ pupils gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, revenue split, CZ-years
- **`trend_gap_share_local_revenue`** — change in the within-CZ share_local gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, revenue split, CZ-years
- **`gap_pupils_ses`** — within-CZ gap, enrollment, group total (location 2 - location 1). M_2/M_1, the enrollment split. B_l is the residual backed out to match it, so this is a target, never a check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_teachers_ses`** — within-CZ gap, FTE teachers, group total (location 2 - location 1). The teacher stock by location. Together with the enrollment gap it says whether the good location's advantage is more teachers or better ones -- the model's answer is 'better', through β.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_teachers_pp_ses`** — within-CZ gap, FTE teachers per pupil (location 2 - location 1). This is H̃_T,l/M_l up to the β-weighting of teacher human capital, so it maps into Q_l = (2H̃_T,l/M_l)^σ directly. If it is near zero while the score gap is not, the quality gap is composition, not class size.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_str_ratio_ses`** — within-CZ gap, student-teacher ratio (location 2 - location 1). Class size. Prop. 1 has better teachers take LARGER classes, which is not what districts do -- check the mapping (notes §5.3) before targeting this.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_salary_per_teacher_ses`** — within-CZ gap, teacher salary per FTE, nominal (location 2 - location 1). Compare with the deflated version: the wedge between the two is how much of the raw salary gap is cost of living rather than pay.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_salary_real_ses`** — within-CZ gap, teacher salary per FTE, CWIFT-deflated (location 2 - location 1). The direct target for κ_2/κ_1; CWIFT-deflated, so it is a teaching-wage gap rather than a local price level.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_exp_pp_ses`** — within-CZ gap, current expenditure per pupil (location 2 - location 1). The resource gap the budget constraint has to deliver.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_rev_pp_local_ses`** — within-CZ gap, local revenue per pupil (location 2 - location 1). The local half of the budget; with the state half it pins how far pure local finance is from the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_rev_pp_state_ses`** — within-CZ gap, state revenue per pupil (location 2 - location 1). Motivates the exogenous transfer G_l of notes §5.2 -- state aid is compensatory, so expect this gap to be NEGATIVE even where the local gap is positive.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_rev_pp_total_ses`** — within-CZ gap, total revenue per pupil (location 2 - location 1). Local plus state plus federal. If this gap has the opposite sign to the local one, redistribution is doing more than local finance, and a pure local budget gets the sign of the resource gap wrong.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_proptax_pp_ses`** — within-CZ gap, property-tax revenue per pupil (location 2 - location 1). The property-tax base per pupil, the closest thing to t_l I_l in the data.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_t_eff_ses`** — within-CZ gap, local revenue / local income base (location 2 - location 1). The effective local tax rate. Keep it endogenous from the budget and use this as an untargeted check.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_share_local_ses`** — within-CZ gap, local share of revenue (location 2 - location 1). If this is far below 1 the pure local-finance budget in §10 is counterfactual -- see G_l.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_share_state_ses`** — within-CZ gap, state share of revenue (location 2 - location 1). The size of the transfer the model is currently missing.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_pov_rate_5_17_ses`** — within-CZ gap, child poverty rate 5-17 (location 2 - location 1). The observable counterpart of sorting on z: how much richer the good location's families are.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_learn_rate_ses`** — within-CZ gap, SEDA learning rate (location 2 - location 1). The closest public analogue to Q_2/Q_1; growth rather than levels, so it nets out the SES composition the model is supposed to generate endogenously.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_score_mean_ses`** — within-CZ gap, SEDA mean score (location 2 - location 1). NOT a Q target -- it mixes school quality with who enrolls. Useful as the untargeted check that the model's sorting reproduces a level gap it did not fit.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_seda_ses_ses`** — within-CZ gap, SEDA composite SES (location 2 - location 1). Same as the poverty gap, on SEDA's composite. Reported in SEDA index units, whose cross-district SD is about 0.9 -- divide by that before setting it beside the model's mean-ability gap (`idx_work_signed`, ≈ +0.08 at the baseline), which is in units of z. The sign is the part that compares directly.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`gap_seda_lninc_ses`** — within-CZ gap, log median household income (location 2 - location 1). The income base I_l, in logs. This is what t_l is levied on, so together with the local-revenue gap it says whether the advantaged location raises more because it is richer or because it taxes itself harder.  
  *Source:* CCD directory + F-33 + SAIPE + CWIFT + SEDA, 2018; district x CZ (ses split), enrollment-weighted across 619 CZs
- **`trend_gap_salary_real_ses`** — change in the within-CZ salary_real gap, 2014->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, ses split, CZ-years
- **`trend_gap_rev_pp_local_ses`** — change in the within-CZ rev_pp_local gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, ses split, CZ-years
- **`trend_gap_exp_pp_ses`** — change in the within-CZ exp_pp gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, ses split, CZ-years
- **`trend_gap_pupils_ses`** — change in the within-CZ pupils gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, ses split, CZ-years
- **`trend_gap_share_local_ses`** — change in the within-CZ share_local gap, 1994->2018. Whether the spatial gaps widened tells you whether a single calibrated steady state is enough or the spatial block needs its own time series. See notes §7: F-33 is thin before 1990.  
  *Source:* F-33 / CCD panel, ses split, CZ-years

## Estimates that rest on a stated assumption

| key | model object | value | ratio 2/1 | within-CZ var share | n |
|---|---|---|---|---|---|
| `sigma_congestion` | σ in Q_l = (2 H̃_T,l / M_l)^σ | -0.002236 (0.00782) | — | — | 3006 |
| `sigma_nu` | σν (Gumbel scale on the location choice) | — | — | — | 3037 |

- **`sigma_congestion`** — class-size / congestion elasticity σ. The slope is NOT distinguishable from zero (t = -0.29), so this is an upper bound on how much the cross-section knows. σ is a free 0.25 in the code. This converts the raw slope (-0.017 SD of achievement gained over grades 3-8 per log point of teachers per pupil) into log-h units with 1 SD ≈ 0.13 log earnings -- the conversion is an assumption, and OLS here is a correlation, not a class-size experiment. Use as a prior, then discipline with STAR / Angrist-Lavy.  
  *Source:* SEDA learning rates x CCD teachers per pupil, 2018, district level, CZ fixed effects, 187 CZ clusters
- **`sigma_nu`** — location-taste dispersion σ_ν. Slope -2.203 = 1/σ_ν under the logit, and it is NOT distinguishable from zero (t = -0.62), so this is an upper bound on how much the cross-section knows; a non-positive slope has no σ_ν behind it, hence the blank value. OLS is biased here -- quality is jointly determined with the share through congestion -- so treat it as an upper bound on 1/σ_ν and keep Eckert-Kleineberg's 0.2 as the prior. Needs an instrument (boundary discontinuities, state-aid formula kinks).  
  *Source:* CCD enrollment shares x SEDA learning rates, 2018, district within CZ, CZ fixed effects

## Literature anchors (not estimated here)

| key | model object | value | ratio 2/1 | within-CZ var share | n |
|---|---|---|---|---|---|
| `lambda_wtp_level` | λ (altruism strength) | — | — | — | — |
| `psi_wtp_gradient` | ψ (warm-glow curvature) | — | — | — | — |
| `moving_cost_m` | mcost / τmove (the goods cost of moving) | — | — | — | — |
| `rho_z_sigma_xi` | ρz, σξ | — | — | — | — |
| `sigma_nu_prior` | σν | 0.2 | — | — | — |
| `exposure_effect` | untargeted validation of the Q_l mechanism | 0.04 | — | — | — |

- **`lambda_wtp_level`** — average WTP for school quality. ~$33/month per SD of test score with boundary fixed effects, ~$17 once neighborhood sociodemographics are controlled -- 1-2% of monthly housing user cost. This pins the LEVEL of λ; ψ is pinned by how WTP varies with parental income (BFM Table 7). Target them separately or they fight for the same moment (notes §4).  
  *Source:* Bayer, Ferreira & McMillan (2007), NBER 13236, Tables 5-7
- **`psi_wtp_gradient`** — income gradient in WTP for school quality. The baseline sets ψ = 0 (linear in h'), which is the strongest sorting the power kernel allows. The gradient in BFM says how much of that is real. This is the single most valuable number still missing from the empirical section.  
  *Source:* Bayer, Ferreira & McMillan (2007) Table 7; ACS district composition by income
- **`moving_cost_m`** — school-quality capitalization at district borders. The boundary-discontinuity house-price premium is the observed price of moving between locations, so `mcost` stops being a free normalization (currently MCOST_BENCH ≈ 0.0294, reverse-engineered from a 0.20 utility cost) and becomes a measured object.  
  *Source:* Black (1999, QJE); Bayer, Ferreira & McMillan (2007)
- **`rho_z_sigma_xi`** — parent-child ability persistence. Must be the parent-child COGNITIVE correlation, not the income IGE: z is ability. Currently ρz = 0.9, σξ = 0.20 by assumption. Cross-check against the stationary AFQT dispersion already matched in the aspatial calibration.  
  *Source:* CNLSY (NLSY79 mothers x Child/Young Adult PIAT and AFQT)
- **`sigma_nu_prior`** — spatial taste dispersion, literature prior. Equals the code's current σν = 0.20 exactly. Reassuring, but it means σν is presently an import rather than an estimate -- the within-CZ share regression above is the way to earn it.  
  *Source:* Eckert & Kleineberg, 'Saving the American Dream?', spatial GE with education policy
- **`exposure_effect`** — neighborhood exposure effect on adult income rank. The single best untargeted test available. Simulate a child moved from location 1 to location 2 at each age and read off the implied convergence rate; ~4%/year without having targeted it would validate the whole Q_l channel (notes §6).  
  *Source:* Chetty & Hendren (2018), NBER 23001; Opportunity Atlas

## What the empirical section still needs

- **Residence-workplace split for teachers** (notes §5.1). Stage 2 gives one Gumbel draw for where you live, pay tax and teach. Within a metro that is the binding assumption: city teachers commonly live in the suburbs. Two independent logits plus a commuting cost buy two independent moment sets (family sorting AND teacher flows).  
  *Source:* SASS/NTPS restricted-use; state administrative files (WI via Biasi-Fu-Stromme, NC, TX); Bates-Dinerstein-Johnston-Sorkin on commute as the dominant teacher preference
- **Teacher ability composition by location** (β, γ, the H̃_T aggregator). The paper's distinctive claim is about WHO teaches, but every moment above is a quantity or a dollar. Nothing here measures teacher ability, so β and γ are untouched by the data.  
  *Source:* NLSY79/97 restricted geocode files (county/MSA of residence) crossed with the AFQT-by-occupation moments already in hand
- **An instrument for the quality-share relationship** (σ, σν). σ_ν and σ are both estimated off an equilibrium relationship between school quality and enrollment that the model says is simultaneous.  
  *Source:* district-boundary discontinuities; state-aid formula kinks (Lafortune-Rothstein-Schanzenbach); court-ordered finance reforms
- **A pre-1990 cross-section** (notes §7). The aspatial calibration has 1970/1990/2010. F-33 is thin before 1990 and SEDA starts in 2009, so the spatial block realistically calibrates to one date.  
  *Source:* Census of Governments school district finances (published back to 1957); Digest of Education Statistics
- **Tract -> school district crosswalk for the Opportunity Atlas** (the untargeted validation in notes §6). OI is tract/county/CZ, never district. The crosswalk is what turns the exposure-effect validation from a hand-wave into a number.  
  *Source:* NCES EDGE district boundaries, population-weighted onto tracts
- **Class size vs. teacher quality within districts** (notes §5.3, and whether gap_str_ratio is a legitimate target). Prop. 1 implies N(h) = h^(β/σ) Q^(-1/σ): better teachers get BIGGER classes. If districts do the opposite, the model's class sizes cannot be read as observed student-teacher ratios.  
  *Source:* CRDC school-level teacher experience/certification (--urban-topics crdc_teachers); state admin data
