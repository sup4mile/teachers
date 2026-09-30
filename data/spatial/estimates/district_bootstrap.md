# District moments: whole-CZ bootstrap

Generated 2026-09-29T16:17:44 by `district_bootstrap.py`. Base year 2018, `locale` partition (suburb NCES locale 21-23 minus city 11-13), labor market: czone (Dorn 1990 commuting zones).

**Method.** The point estimates are computed with `data_estimate.py`'s own panel, screen, partition and `czone_gaps` (enrollment-weighted mean across CZs of the within-CZ gap, weights = partitioned CZ pupils). The bootstrap resamples whole CZs with replacement: **B = 1999, seed = 20260929** (`numpy.random.default_rng`; one multiplicity matrix shared by all statistics), 188 CZs / 3291 districts in the partition. Each draw recomputes the pupil-weighted mean over the resampled CZs that have the statistic (so n varies by draw). Held fixed: the district panel after the within-year [0.5%, 99.5%] screen, the partition and the CZ filters (500-pupil and 10%-group-share cuts). `SE` = SD of the draws; the interval is the 2.5-97.5 percentile. `placeholder` = unweighted cross-CZ SD / sqrt(n) as currently used in `spatial_calibrate.jl`; `lin. SE` is the analytic cluster-robust (linearised ratio) SE, a cross-check; `n_eff` = 1 / sum of squared normalised CZ weights; `top5` = share of the weight on the five largest CZs.

**Reproduction.** Rebuilt from the raw 2018 tables, the 20 statistics that exist in the tracked `spatial_moments.json` match it (max absolute difference 1.1e-16; identical n for 20 of 20).


## Targets

`salary_real` is the CWIFT-deflated teacher salary per FTE (teacher-weighted within location); `seda_lninc` is log median household income from SEDA covariates (the current income proxy).

| statistic | value | n CZs | boot SE | 95% pct interval | placeholder SE | boot / placeholder | lin. SE | n_eff | top5 |
|---|---|---|---|---|---|---|---|---|---|
| `pupils` | 0.3764 | 188 | 0.1154 | [0.1577, 0.6041] | 0.0707 | 1.63 | 0.1177 | 42 | 26% |
| `salary_real` | 0.0162 | 176 | 0.0090 | [-0.0005, 0.0339] | 0.0077 | 1.17 | 0.0089 | 41 | 24% |
| `teachers_pp` | -0.0122 | 188 | 0.0139 | [-0.0398, 0.0138] | 0.0095 | 1.46 | 0.0138 | 42 | 26% |
| `seda_lninc` | 0.2708 | 187 | 0.0268 | [0.2192, 0.3247] | 0.0150 | 1.78 | 0.0279 | 42 | 26% |

- `pupils`: enrollment, group total. Tracked value 0.37635 (n = 188), difference +0.00e+00.
- `salary_real`: teacher salary per FTE, CWIFT-deflated. Tracked value 0.01623 (n = 176), difference +0.00e+00.
- `teachers_pp`: FTE teachers per pupil. Tracked value -0.01216 (n = 188), difference +0.00e+00.
- `seda_lninc`: log median household income. Tracked value 0.27080 (n = 187), difference +5.55e-17.

## Validation moments

`local_revenue_share` is a level (pupil-weighted mean of the local revenue share over 3249 sample districts); its `n` column counts CZ clusters and its placeholder is the SD of CZ-level shares / sqrt(CZs).

| statistic | value | n CZs | boot SE | 95% pct interval | placeholder SE | boot / placeholder | lin. SE | n_eff | top5 |
|---|---|---|---|---|---|---|---|---|---|
| `pov_rate_5_17` | -0.0888 | 188 | 0.0090 | [-0.1074, -0.0721] | 0.0061 | 1.47 | 0.0093 | 42 | 26% |
| `learn_rate` | 0.0031 | 187 | 0.0030 | [-0.0027, 0.0092] | 0.0021 | 1.41 | 0.0030 | 42 | 26% |
| `score_mean` | 0.2887 | 187 | 0.0326 | [0.2269, 0.3549] | 0.0225 | 1.45 | 0.0335 | 42 | 26% |
| `local_revenue_share` | 0.4726 | 188 | 0.0206 | [0.4323, 0.5113] | 0.0096 | 2.14 | 0.0210 | 43 | 24% |

- `pov_rate_5_17`: child poverty rate 5-17. Tracked value -0.08878 (n = 188), difference +0.00e+00.
- `learn_rate`: SEDA learning rate. Tracked value 0.00307 (n = 187), difference +0.00e+00.
- `score_mean`: SEDA mean score. Tracked value 0.28873 (n = 187), difference +0.00e+00.
- `local_revenue_share`: local share of district revenue, pupil-weighted over sample districts (level; CZ-level shares for the placeholder SD). Tracked value 0.47259 (n = 3249), difference +0.00e+00.

## Income alternatives (ACS 5-year 2014-18, school districts)

Same partition, same CZ weights. Within a location: pupil-weighted mean of the district log median (`ln_*`), or of the level then logged (`lvl_*`), or weighted by earners (`ew_*`). Suburb minus city. Districts with a bottom-coded (2,499) or missing median are dropped; a top code (250,001) is kept.

| statistic | value | n CZs | boot SE | 95% pct interval | placeholder SE | boot / placeholder | lin. SE | n_eff | top5 |
|---|---|---|---|---|---|---|---|---|---|
| `ln_earn_all` | 0.1578 | 188 | 0.0214 | [0.1174, 0.2007] | 0.0138 | 1.55 | 0.0221 | 42 | 26% |
| `ln_earn_male` | 0.2180 | 188 | 0.0229 | [0.1739, 0.2625] | 0.0152 | 1.51 | 0.0235 | 42 | 26% |
| `ln_earn_female` | 0.0996 | 188 | 0.0195 | [0.0621, 0.1386] | 0.0129 | 1.50 | 0.0200 | 42 | 26% |
| `ln_earn_ftyr` | 0.1288 | 188 | 0.0183 | [0.0949, 0.1653] | 0.0099 | 1.85 | 0.0188 | 42 | 26% |
| `ln_earn_ftyr_male` | 0.1624 | 188 | 0.0203 | [0.1239, 0.2027] | 0.0111 | 1.84 | 0.0208 | 42 | 26% |
| `ln_earn_ftyr_female` | 0.0752 | 188 | 0.0174 | [0.0436, 0.1107] | 0.0087 | 2.02 | 0.0179 | 42 | 26% |
| `ln_earn_25p` | 0.1439 | 188 | 0.0194 | [0.1072, 0.1833] | 0.0111 | 1.75 | 0.0199 | 42 | 26% |
| `ln_earn_25p_male` | 0.1982 | 188 | 0.0214 | [0.1569, 0.2402] | 0.0131 | 1.64 | 0.0221 | 42 | 26% |
| `ln_earn_25p_female` | 0.0767 | 188 | 0.0182 | [0.0428, 0.1140] | 0.0099 | 1.84 | 0.0187 | 42 | 26% |
| `ln_earn_emp` | 0.1453 | 188 | 0.0203 | [0.1069, 0.1856] | 0.0119 | 1.71 | 0.0209 | 42 | 26% |
| `ln_hhinc_2544` | 0.2493 | 188 | 0.0276 | [0.1979, 0.3070] | 0.0162 | 1.70 | 0.0288 | 42 | 26% |
| `ln_hhinc_acs` | 0.2627 | 188 | 0.0280 | [0.2102, 0.3185] | 0.0155 | 1.81 | 0.0292 | 42 | 26% |
| `ln_earn_hh_mean` | 0.1418 | 188 | 0.0246 | [0.0964, 0.1927] | 0.0139 | 1.77 | 0.0257 | 42 | 26% |
| `lvl_earn_all` | 0.1671 | 188 | 0.0215 | [0.1253, 0.2101] | 0.0135 | 1.59 | 0.0222 | 42 | 26% |
| `lvl_earn_ftyr` | 0.1417 | 188 | 0.0185 | [0.1070, 0.1785] | 0.0102 | 1.82 | 0.0190 | 42 | 26% |
| `lvl_hhinc_2544` | 0.2662 | 188 | 0.0272 | [0.2159, 0.3218] | 0.0165 | 1.65 | 0.0284 | 42 | 26% |
| `ew_ln_earn_all` | 0.1528 | 188 | 0.0196 | [0.1159, 0.1921] | 0.0147 | 1.33 | 0.0201 | 42 | 26% |
| `ew_ln_earn_ftyr` | 0.1228 | 188 | 0.0166 | [0.0921, 0.1567] | 0.0098 | 1.70 | 0.0170 | 42 | 26% |
| `pm_earn_all` | 0.1435 | 188 | 0.0196 | [0.1068, 0.1829] | 0.0140 | 1.41 | 0.0201 | 42 | 26% |
| `pm_earn_male` | 0.1944 | 188 | 0.0213 | [0.1534, 0.2370] | 0.0153 | 1.39 | 0.0219 | 42 | 26% |
| `pm_earn_female` | 0.0907 | 188 | 0.0184 | [0.0569, 0.1283] | 0.0134 | 1.37 | 0.0189 | 42 | 26% |
| `pm_ftyr_all` | 0.1122 | 188 | 0.0168 | [0.0808, 0.1466] | 0.0093 | 1.80 | 0.0172 | 42 | 26% |
| `pm_ftyr_male` | 0.1437 | 188 | 0.0180 | [0.1089, 0.1797] | 0.0104 | 1.73 | 0.0185 | 42 | 26% |
| `pm_ftyr_female` | 0.0613 | 188 | 0.0163 | [0.0325, 0.0951] | 0.0084 | 1.95 | 0.0167 | 42 | 26% |

- `ln_earn_all`: log median earnings, 16+ with earnings (B20002_001).
- `ln_earn_male`: log median earnings, men 16+ with earnings (B20002_002).
- `ln_earn_female`: log median earnings, women 16+ with earnings (B20002_003).
- `ln_earn_ftyr`: log median earnings, full-time year-round 16+ (B20018_001).
- `ln_earn_ftyr_male`: log median earnings, men full-time year-round (B20017_003).
- `ln_earn_ftyr_female`: log median earnings, women full-time year-round (B20017_006).
- `ln_earn_25p`: log median earnings, 25+ with earnings (B20004_001).
- `ln_earn_25p_male`: log median earnings, men 25+ (B20004_007).
- `ln_earn_25p_female`: log median earnings, women 25+ (B20004_013).
- `ln_earn_emp`: log median earnings, civilian employed 16+ (B24011_001).
- `ln_hhinc_2544`: log median household income, householder 25-44 (B19049_003).
- `ln_hhinc_acs`: log median household income, ACS 2014-18 (B19013_001).
- `ln_earn_hh_mean`: log mean earnings per earning household (B19061_001/B19051_002).
- `lvl_earn_all`: log ratio of location-average median earnings, 16+ with earnings (B20002_001).
- `lvl_earn_ftyr`: log ratio of location-average median earnings, full-time year-round 16+ (B20018_001).
- `lvl_hhinc_2544`: log ratio of location-average median household income, householder 25-44 (B19049_003).
- `ew_ln_earn_all`: log median earnings, 16+ with earnings (B20002_001), earner-weighted within location.
- `ew_ln_earn_ftyr`: log median earnings, full-time year-round 16+ (B20018_001), earner-weighted within location.
- `pm_earn_all`: log pooled location median earnings, all workers 16+ with earnings (B20001 bins summed over the districts of each location).
- `pm_earn_male`: log pooled location median earnings, men 16+ with earnings (B20001 bins summed over the districts of each location).
- `pm_earn_female`: log pooled location median earnings, women 16+ with earnings (B20001 bins summed over the districts of each location).
- `pm_ftyr_all`: log pooled location median earnings, full-time year-round workers 16+ (B20005 bins summed over the districts of each location).
- `pm_ftyr_male`: log pooled location median earnings, men, full-time year-round (B20005 bins summed over the districts of each location).
- `pm_ftyr_female`: log pooled location median earnings, women, full-time year-round (B20005 bins summed over the districts of each location).

### Reading the income alternatives

- The model's statistic (`income = :log_median_market` in `spatial_calibrate.jl`) is the gap in log median earnings of market workers. The closest ACS school-district objects are the median of earnings for everyone 16+ with earnings and, for full-time year-round workers, B20018/B20017. No school-district table cuts earnings by age, so 25-34 cannot be isolated; the only age-specific object is household income of householders 25-44 (B19049_003), which is a household, not a labor-income, concept.
- The current proxy (`seda_lninc`, log median household income) is 0.271 (SE 0.027). Earnings-based gaps are about half as large: 0.158 (SE 0.021) for the average of district median earnings (B20002_001), 0.144 (SE 0.020) for the pooled location median of everyone with earnings, and 0.112 (SE 0.017) for the pooled median of full-time year-round workers. Household income of 25-44 householders is 0.249 (SE 0.028) against 0.263 for all householders (ACS B19013), close to the SEDA proxy, while mean earnings per earning household (B19061/B19051) is 0.142. So the wedge is between household income (a median over all households, with non-labor income) and earnings, not age.
- `pm_*` (pooled location median) is the closest analogue of the model object: the median is taken over the pooled earnings distribution of a location (summed B20001/B20005 bin counts over the location's districts, linear interpolation inside bins) rather than averaging district medians. Interpolation reproduces the published district medians with a mean log error of +0.008 (SD 0.012).
- The men/women rows show the source of the level difference: women's earnings gap (city vs suburb) is smaller than men's, and FTYR restrictions shrink both.

## All other gaps in `data_estimate.py`

| statistic | value | n CZs | boot SE | 95% pct interval | placeholder SE | boot / placeholder | lin. SE | n_eff | top5 |
|---|---|---|---|---|---|---|---|---|---|
| `teachers` | 0.3632 | 188 | 0.1126 | [0.1544, 0.5898] | 0.0708 | 1.59 | 0.1148 | 42 | 26% |
| `str_ratio` | 0.0122 | 188 | 0.0139 | [-0.0138, 0.0398] | 0.0095 | 1.46 | 0.0138 | 42 | 26% |
| `salary_per_teacher` | -0.0171 | 177 | 0.0180 | [-0.0536, 0.0140] | 0.0080 | 2.26 | 0.0181 | 41 | 26% |
| `exp_pp` | -0.0956 | 188 | 0.0157 | [-0.1291, -0.0659] | 0.0103 | 1.52 | 0.0154 | 42 | 26% |
| `rev_pp_local` | -0.0708 | 188 | 0.0558 | [-0.1718, 0.0503] | 0.0330 | 1.69 | 0.0558 | 42 | 26% |
| `rev_pp_state` | -0.0003 | 188 | 0.0678 | [-0.1356, 0.1417] | 0.0247 | 2.74 | 0.0670 | 42 | 26% |
| `rev_pp_total` | -0.1211 | 188 | 0.0163 | [-0.1549, -0.0906] | 0.0109 | 1.49 | 0.0158 | 42 | 26% |
| `proptax_pp` | -0.1041 | 166 | 0.0473 | [-0.1884, -0.0066] | 0.0342 | 1.38 | 0.0479 | 33 | 31% |
| `t_eff` | 0.0006 | 188 | 0.0008 | [-0.0008, 0.0024] | 0.0005 | 1.50 | 0.0008 | 42 | 26% |
| `share_local` | 0.0105 | 188 | 0.0214 | [-0.0286, 0.0557] | 0.0105 | 2.04 | 0.0215 | 42 | 26% |
| `share_state` | 0.0196 | 188 | 0.0203 | [-0.0227, 0.0571] | 0.0096 | 2.11 | 0.0203 | 42 | 26% |
| `seda_ses` | 0.8299 | 187 | 0.1013 | [0.6377, 1.0350] | 0.0576 | 1.76 | 0.1052 | 42 | 26% |

- `teachers`: FTE teachers, group total. Tracked value 0.36315 (n = 188), difference +0.00e+00.
- `str_ratio`: student-teacher ratio. Tracked value 0.01216 (n = 188), difference +0.00e+00.
- `salary_per_teacher`: teacher salary per FTE, nominal. Tracked value -0.01705 (n = 177), difference +0.00e+00.
- `exp_pp`: current expenditure per pupil. Tracked value -0.09564 (n = 188), difference +0.00e+00.
- `rev_pp_local`: local revenue per pupil. Tracked value -0.07082 (n = 188), difference +0.00e+00.
- `rev_pp_state`: state revenue per pupil. Tracked value -0.00031 (n = 188), difference +0.00e+00.
- `rev_pp_total`: total revenue per pupil. Tracked value -0.12114 (n = 188), difference +0.00e+00.
- `proptax_pp`: property-tax revenue per pupil. Tracked value -0.10412 (n = 166), difference +1.39e-17.
- `t_eff`: local revenue / local income base. Tracked value 0.00060 (n = 188), difference +0.00e+00.
- `share_local`: local share of revenue. Tracked value 0.01055 (n = 188), difference +0.00e+00.
- `share_state`: state share of revenue. Tracked value 0.01963 (n = 188), difference +0.00e+00.
- `seda_ses`: SEDA composite SES. Tracked value 0.82991 (n = 187), difference -1.11e-16.

## Bootstrap correlation of the targets and validation moments

| | `pupils` | `salary_real` | `teachers_pp` | `seda_lninc` | `pov_rate_5_17` | `learn_rate` | `score_mean` | `local_revenue_share` |
|---|---|---|---|---|---|---|---|---|
| `pupils` | 1.00 | -0.04 | -0.26 | 0.53 | -0.59 | 0.29 | 0.36 | 0.32 |
| `salary_real` | -0.04 | 1.00 | 0.26 | 0.32 | -0.29 | -0.12 | 0.34 | 0.16 |
| `teachers_pp` | -0.26 | 0.26 | 1.00 | 0.07 | 0.08 | -0.65 | 0.10 | 0.25 |
| `seda_lninc` | 0.53 | 0.32 | 0.07 | 1.00 | -0.92 | 0.05 | 0.91 | 0.52 |
| `pov_rate_5_17` | -0.59 | -0.29 | 0.08 | -0.92 | 1.00 | -0.23 | -0.86 | -0.49 |
| `learn_rate` | 0.29 | -0.12 | -0.65 | 0.05 | -0.23 | 1.00 | -0.03 | -0.19 |
| `score_mean` | 0.36 | 0.34 | 0.10 | 0.91 | -0.86 | -0.03 | 1.00 | 0.37 |
| `local_revenue_share` | 0.32 | 0.16 | 0.25 | 0.52 | -0.49 | -0.19 | 0.37 | 1.00 |

Full draws: `district_bootstrap_draws.csv` (one row per draw).
