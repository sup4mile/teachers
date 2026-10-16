# District moments under three settings: locale, revenue and SES splits

Generated 2026-10-02T23:48:11 by `district_settings.py` (year 2018, whole-CZ bootstrap B = 1999, seed 20260929). For each way of cutting a commuting zone in two it computes the audited moment set used for the locale calibration targets (`district_audit.md` variant 28 and FTE variant 11, `district_bootstrap.md` pooled earnings medians) with the same corrections, the validation moments, and a whole-CZ bootstrap of every one of them. The scripts it builds on are unchanged.

## Summary

1. *Reproduction.* The locale full-sample targets match the adopted values: gap_salary 0.0086 (0.0092), gap_teachers_pp -0.0177 (0.0131), gap_pupils 0.3764 (0.1154), gap_income 0.1122 (0.0168). 19 moments are also compared with the full-precision values stored in `district_audit.json` / `district_bootstrap.json`: largest absolute difference 0.0e+00 in the value and 1.4e-17 in the bootstrap SE. The tracked-definition moments of all three partitions equal `spatial_moments.json` (locale: 14 moments, max |diff| 1.1e-16, same n for 14; revenue: 14 moments, max |diff| 5.6e-17, same n for 14; ses: 14 moments, max |diff| 6.9e-18, same n for 14).
2. *Coverage.* The locale split needs a city and a suburb in the same CZ (188 CZs, 3,260 districts, 28.6 m pupils); the revenue and ses splits need only two groups of districts (621 CZs, 9,245 districts, 44.7 m pupils; 619 CZs, 8,869 districts, 43.4 m). 186 CZs are in all three. They hold 99.9% of the locale pupils but 75.2% and 74.8% of the revenue and ses pupils, so the full-sample revenue and ses averages are weighted toward many smaller CZs (effective number of CZs 42 locale, 89 revenue, 94 ses). New York City (32 geographic districts, 960,484 pupils) is in the locale partition and in neither of the others: its split variables (F-33 revenue, SEDA SES) are filed under the Chancellor's Office id and are missing for the geographic districts.
3. *The four targets, same 186 CZs* (value, bootstrap SE):

| target | locale | revenue | ses |
|---|---|---|---|
| `gap_pupils` | 0.3772 (0.1194) | -0.0029 (0.0215) | 0.0104 (0.0203) |
| `gap_salary` | 0.0087 (0.0092) | 0.0211 (0.0126) | 0.0136 (0.0072) |
| `gap_teachers_pp` | -0.0169 (0.0127) | 0.0356 (0.0068) | -0.0088 (0.0074) |
| `gap_income` | 0.1123 (0.0169) | 0.1516 (0.0140) | 0.2725 (0.0135) |

   The pupil-weighted median split equalises enrollment within a CZ, so the enrollment gap is about zero under revenue (-0.003) and ses (0.010) against 0.377 for locale: it no longer measures a size difference. The earnings gap is 0.112 (locale), 0.152 (revenue), 0.272 (ses); the poverty gap -0.089, -0.042, -0.113; the SEDA growth gap 0.0030, 0.0080, 0.0158.
4. *State-aid offsets (revenue split, common CZs).* The local revenue gap is $4,362 per pupil (log 0.587); state aid offsets 0.52 of each local dollar (SE 0.04) and federal aid 0.06, leaving a total-revenue gap of $1,864 (log 0.118), a current-expenditure gap of $1,149 and a teacher-payroll gap of $431 per K-12 pupil (log 0.084). The within-CZ district regression gives a state-on-local slope of -0.47 (SE 0.05, 9,206 districts, 621 CZs).
5. *NYC and pre-K.* With NYC filled in as one fully observed district, the revenue-split enrollment gap moves from -0.025 to 0.020 (SE 0.047) and the locale earnings gap from 0.112 to 0.123 (section 9). Splitting on local revenue per K-12 rather than per CCD pupil moves 1.5% of pupils to the other location and changes the revenue-split targets by at most 0.0059 (section 9).
6. *Hold-out.* The pre-registered split puts 94 of the 188 locale CZs in half A (values in section 8) and 94 in the sealed half B; B's results are in `sealed/holdout_B.json` and are not reported here.

## 1. Method

**Settings.** `data_estimate.partition` at its tracked sample rules (districts with >= 500 pupils in a CZ, each location >= 10% of the CZ's partitioned pupils). Location 2 is the advantaged group: suburb (locale), above the CZ's pupil-weighted median of local revenue per pupil (revenue), above the median SEDA SES (ses). The revenue and ses cut points are computed on the districts that have the split variable, after the panel's [0.5%, 99.5%] screen of `rev_pp_local`.

**Moments.** Each is the pupil-weighted mean across CZs (weights = the CZ's partitioned pupils) of a within-CZ gap, location 2 minus location 1; a CZ missing the statistic drops out of the mean (and from that statistic's `n CZs`). Log gaps are `log(x_2/x_1)` of location aggregates, level gaps are `x_2 - x_1`.

| moment | definition |
|---|---|
| `gap_pupils` | log ratio of location totals of CCD enrollment (pre-K included), the tracked definition |
| `gap_salary` | audit variant 28: teacher wage bill per K-12 FTE, divided by CWIFT, teacher-FTE-weighted mean over a location's districts, [0.5%, 99.5%] screen panel-wide. Zero-coded F-33 wage bills (IL, AK, NM and others) are filled with 0.881 x instruction salaries (the median teacher / instruction ratio of the locale sample, the same constant in every setting), the 32 NYC geographic districts are one district with the Chancellor's Office F-33 record, and the FTE denominator excludes pre-K teacher FTE |
| `gap_teachers_pp` | audit FTE variant 11: log ratio of location totals of (teachers - pre-K teacher FTE) / (pupils - pre-K enrollment) |
| `gap_income` | `pm_ftyr_all`: log of the location median of full-time year-round earnings of workers 16+, read off ACS 5-year 2014-18 B20005 bin counts summed over the location's districts, linear interpolation inside bins; `pm_earn_all` is the same from B20001 (all workers with earnings) |
| `poverty_gap` | difference in the SAIPE child poverty rate 5-17 (ratio of location totals), levels |
| `seda_growth_gap`, `score_gap` | SEDA 6.0 learning rate and mean score, enrollment-weighted district means, difference in levels (SEDA's published units) |
| `seda_ses`, `seda_lninc` | SEDA composite SES and log median household income, difference |
| `rev_pp_*`, `exp_pp` | log ratio of location totals of local, state, federal, total revenue and of current expenditure per CCD pupil; `rev_pp_fed` is F-33 federal revenue (total = local + state + federal holds exactly) |
| `payroll_pp_k12` | log ratio of location totals of the audited (hybrid) teacher wage bill per K-12 pupil, unscreened like the other per-pupil ratios, districts with a positive wage bill |
| `share_local`, `share_state` | difference in the location ratio of local (state) to total revenue |
| `local_revenue_share` | LEVEL: pupil-weighted mean of the districts' local revenue share over the partition's districts, as in `district_bootstrap.py` |
| `salary_tracked`, `teachers_pp_tracked` | memo: the tracked definitions the corrections replace |
| `usd_*` | memo: location 2 minus location 1 in dollars per pupil, which add up across sources |

**Bootstrap.** Whole CZs resampled with replacement, B = 1999, seed 20260929, `numpy.random.default_rng`; one (B x n_CZ) multiplicity matrix per CZ set shared by all of its statistics, so the correlations are coherent. SE = SD of the draws; the interval is the 2.5-97.5 percentile. As in `district_bootstrap.py` the panel, screen, partition and CZ filters are held fixed and only CZ membership of the average is resampled; a restricted CZ set (common CZs, hold-out half) is resampled over its own CZs.

## 2. Reproduction check (locale, full sample)

The four adopted targets must reproduce the digits printed in `district_audit.md` / `district_bootstrap.md`; every other moment that exists in `district_audit.json` / `district_bootstrap.json` is compared at full precision.

| moment | reported (value, SE) | this script | n CZs | value - stored | SE - stored |
|---|---|---|---|---|---|
| `gap_salary` | 0.0086 (0.0092) | 0.0086 (0.0092) | 188 | +0.0e+00 | +0.0e+00 |
| `gap_teachers_pp` | -0.0177 (0.0131) | -0.0177 (0.0131) | 188 | +0.0e+00 | +0.0e+00 |
| `gap_pupils` | 0.3764 (0.1154) | 0.3764 (0.1154) | 188 | +0.0e+00 | -1.4e-17 |
| `gap_income` | 0.1122 (0.0168) | 0.1122 (0.0168) | 188 | +0.0e+00 | +0.0e+00 |
| `pm_earn_all` |  | 0.1435 (0.0196) | 188 | +0.0e+00 | +0.0e+00 |
| `poverty_gap` |  | -0.0888 (0.0090) | 188 | +0.0e+00 | -1.7e-18 |
| `seda_growth_gap` |  | 0.0031 (0.0030) | 187 | +0.0e+00 | +0.0e+00 |
| `score_gap` |  | 0.2887 (0.0326) | 187 | +0.0e+00 | +0.0e+00 |
| `seda_ses` |  | 0.8299 (0.1013) | 187 | +0.0e+00 | -1.4e-17 |
| `seda_lninc` |  | 0.2708 (0.0268) | 187 | +0.0e+00 | +3.5e-18 |
| `rev_pp_local` |  | -0.0708 (0.0558) | 188 | +0.0e+00 | +0.0e+00 |
| `rev_pp_state` |  | -0.0003 (0.0678) | 188 | +0.0e+00 | +0.0e+00 |
| `rev_pp_total` |  | -0.1211 (0.0163) | 188 | +0.0e+00 | +0.0e+00 |
| `exp_pp` |  | -0.0956 (0.0157) | 188 | +0.0e+00 | +0.0e+00 |
| `share_local` |  | 0.0105 (0.0214) | 188 | +0.0e+00 | +0.0e+00 |
| `share_state` |  | 0.0196 (0.0203) | 188 | +0.0e+00 | +0.0e+00 |
| `local_revenue_share` |  | 0.4726 (0.0206) | 188 | +0.0e+00 | +0.0e+00 |
| `salary_tracked` |  | 0.0162 (0.0090) | 176 | +0.0e+00 | +0.0e+00 |
| `teachers_pp_tracked` |  | -0.0122 (0.0139) | 188 | +0.0e+00 | +0.0e+00 |

All four targets match to the reported digits: **yes** (standard errors checked: yes). The tracked-definition moments of all three partitions equal `spatial_moments.json` (locale: 14 moments, max |diff| 1.1e-16, same n for 14; revenue: 14 moments, max |diff| 5.6e-17, same n for 14; ses: 14 moments, max |diff| 6.9e-18, same n for 14).

## 3. Coverage

| | locale | revenue | ses | locale, common CZs | revenue, common CZs | ses, common CZs |
|---|---|---|---|---|---|---|
| CZs | 188 | 621 | 619 | 186 | 186 | 186 |
| districts | 3,260 | 9,245 | 8,869 | 3,252 | 6,186 | 5,904 |
| pupils (millions) | 28.55 | 44.66 | 43.45 | 28.52 | 33.57 | 32.49 |
| share of all regular-district pupils | 59.1% | 92.4% | 89.9% | 59.0% | 69.5% | 67.2% |
| districts per CZ, median | 6 | 7 | 7 | 6 | 22 | 21 |
| CZs with <= 3 districts | 27.1% | 19.5% | 19.5% | 26.9% | 0.5% | 0.5% |
| effective number of CZs (1 / sum of squared weights) | 42 | 89 | 94 | 42 | 55 | 57 |
| weight on the 5 largest CZs | 25.8% | 15.7% | 15.1% | 25.8% | 20.8% | 20.1% |
| pupils in location 2 | 57.9% | 49.5% | 50.0% | 57.9% | 49.9% | 50.2% |
| NYC pupils in the partition | 960,484 | 0 | 0 | 960,484 | 0 | 0 |
| pupils with a valid audited salary (variant 28) | 98.9% | 99.2% | 99.4% | 98.9% | 99.2% | 99.4% |
| pupils with a valid tracked salary (positive, screened) | 89.9% | 94.2% | 94.8% | 90.0% | 93.8% | 94.5% |
| CZs with an audited salary gap | 188 | 619 | 615 | 186 | 186 | 186 |
| ... their share of the partition's pupils | 100.0% | 100.0% | 99.9% | 100.0% | 100.0% | 100.0% |
| pupils with valid K-12 FTE and K-12 pupils | 99.9% | 99.9% | 99.9% | 99.9% | 99.9% | 99.9% |
| pupils with ACS full-time year-round earnings bins | 96.1% | 99.5% | 99.8% | 96.1% | 99.6% | 99.8% |
| pupils with a SAIPE poverty rate | 96.2% | 99.6% | 99.9% | 96.1% | 99.6% | 99.8% |
| pupils with a SEDA learning rate | 92.9% | 97.1% | 99.8% | 92.9% | 96.6% | 99.8% |
| pupils with a SEDA mean score | 93.0% | 97.2% | 99.9% | 93.0% | 96.7% | 99.8% |
| pupils with F-33 local and total revenue | 96.4% | 100.0% | 99.9% | 96.4% | 100.0% | 99.9% |

The 186 common CZs are those in all three partitions. The locale partition has CZ 34115 (Fairbanks North Star Borough School District, 0.07% of locale weight) which is not in the revenue partition: 59% of its pupils in districts of >= 500 pupils have the split variable. The locale partition has CZ 20200 (Champlain Valley Unified School District, 0.05% of locale weight) which is not in the ses partition: 0% of its pupils in districts of >= 500 pupils have the split variable. 434 revenue-only and 432 ses-only CZs are outside the locale partition (they lack a city group or a suburb group with >= 10% of the CZ's partitioned pupils). Hybrid share (median teacher / instruction salaries, non-zero districts, NYC excluded) by partition: locale 0.881, revenue 0.872, ses 0.872; the 0.881 of the audit is used everywhere.

## 4. Audited moments by setting (full samples)

### locale

locale split (city vs suburb), the 188 CZs of the tracked moments: 188 CZs, 3,260 districts, 28.55 m pupils.

| moment | what | role | value | boot SE | 95% interval | n CZs |
|---|---|---|---|---|---|---|
| `gap_pupils` | enrollment, log(M2/M1) (CCD enrollment incl. pre-K) | target | 0.3764 | 0.1154 | [0.1577, 0.6041] | 188 |
| `gap_salary` | CWIFT-deflated teacher salary per K-12 FTE, audit variant 28 | target | 0.0086 | 0.0092 | [-0.0096, 0.0263] | 188 |
| `gap_teachers_pp` | K-12 FTE teachers per K-12 pupil, audit FTE variant 11 | target | -0.0177 | 0.0131 | [-0.0436, 0.0077] | 188 |
| `gap_income` | log pooled median earnings, full-time year-round workers (`pm_ftyr_all`) | target | 0.1122 | 0.0168 | [0.0808, 0.1466] | 188 |
| `pm_earn_all` | log pooled median earnings, all workers with earnings (`pm_earn_all`) | target | 0.1435 | 0.0196 | [0.1068, 0.1829] | 188 |
| `poverty_gap` | child poverty rate 5-17 (SAIPE), difference in levels | validation | -0.0888 | 0.0090 | [-0.1074, -0.0721] | 188 |
| `seda_growth_gap` | SEDA learning rate (`learn_rate`), difference | validation | 0.0031 | 0.0030 | [-0.0027, 0.0092] | 187 |
| `score_gap` | SEDA mean score (`score_mean`), difference | validation | 0.2887 | 0.0326 | [0.2269, 0.3549] | 187 |
| `seda_ses` | SEDA composite SES, difference | validation | 0.8299 | 0.1013 | [0.6377, 1.0350] | 187 |
| `seda_lninc` | log median household income (SEDA), difference | validation | 0.2708 | 0.0268 | [0.2192, 0.3247] | 187 |
| `rev_pp_local` | local revenue per pupil | diagnostic | -0.0708 | 0.0558 | [-0.1718, 0.0503] | 188 |
| `rev_pp_state` | state revenue per pupil | diagnostic | -0.0003 | 0.0678 | [-0.1356, 0.1417] | 188 |
| `rev_pp_fed` | federal revenue per pupil | diagnostic | -0.5677 | 0.0546 | [-0.6764, -0.4614] | 188 |
| `rev_pp_total` | total revenue per pupil | diagnostic | -0.1211 | 0.0163 | [-0.1549, -0.0906] | 188 |
| `exp_pp` | current expenditure per pupil | diagnostic | -0.0956 | 0.0157 | [-0.1291, -0.0659] | 188 |
| `payroll_pp_k12` | audited teacher wage bill per K-12 pupil | diagnostic | -0.0267 | 0.0201 | [-0.0681, 0.0106] | 188 |
| `share_local` | local share of revenue, difference in levels | diagnostic | 0.0105 | 0.0214 | [-0.0286, 0.0557] | 188 |
| `share_state` | state share of revenue, difference in levels | diagnostic | 0.0196 | 0.0203 | [-0.0227, 0.0571] | 188 |
| `local_revenue_share` | local share of revenue, LEVEL (pupil-weighted over the partition's districts) | diagnostic | 0.4726 | 0.0206 | [0.4323, 0.5113] | 188 |
| `salary_tracked` | tracked `salary_real` gap (zero-coded wage bills kept, NYC absent, pre-K in FTE) | memo (tracked definition) | 0.0162 | 0.0090 | [-0.0005, 0.0339] | 176 |
| `teachers_pp_tracked` | tracked FTE-per-pupil gap (pre-K teachers and pupils included) | memo (tracked definition) | -0.0122 | 0.0139 | [-0.0398, 0.0138] | 188 |
| `usd_rev_pp_local` | local revenue, $ per pupil, 2 minus 1 | $ per pupil | -768 | 622 | [-2023, 463] | 188 |
| `usd_rev_pp_state` | state revenue, $ per pupil, 2 minus 1 | $ per pupil | -772 | 517 | [-1867, 127] | 188 |
| `usd_rev_pp_fed` | federal revenue, $ per pupil, 2 minus 1 | $ per pupil | -643 | 67 | [-781, -518] | 188 |
| `usd_rev_pp_total` | total revenue, $ per pupil, 2 minus 1 | $ per pupil | -2184 | 346 | [-2916, -1537] | 188 |
| `usd_exp_pp` | current expenditure, $ per pupil, 2 minus 1 | $ per pupil | -1301 | 258 | [-1826, -817] | 188 |
| `usd_payroll_k12` | teacher payroll (K-12 pupils), $ per pupil, 2 minus 1 | $ per pupil | -80 | 131 | [-340, 172] | 188 |

### revenue

revenue split (below / above the CZ's pupil-weighted median of local revenue per pupil): 621 CZs, 9,245 districts, 44.66 m pupils.

| moment | what | role | value | boot SE | 95% interval | n CZs |
|---|---|---|---|---|---|---|
| `gap_pupils` | enrollment, log(M2/M1) (CCD enrollment incl. pre-K) | target | -0.0247 | 0.0218 | [-0.0673, 0.0156] | 621 |
| `gap_salary` | CWIFT-deflated teacher salary per K-12 FTE, audit variant 28 | target | 0.0178 | 0.0098 | [0.0001, 0.0374] | 619 |
| `gap_teachers_pp` | K-12 FTE teachers per K-12 pupil, audit FTE variant 11 | target | 0.0348 | 0.0057 | [0.0233, 0.0458] | 621 |
| `gap_income` | log pooled median earnings, full-time year-round workers (`pm_ftyr_all`) | target | 0.1286 | 0.0114 | [0.1053, 0.1508] | 620 |
| `pm_earn_all` | log pooled median earnings, all workers with earnings (`pm_earn_all`) | target | 0.1000 | 0.0117 | [0.0758, 0.1230] | 620 |
| `poverty_gap` | child poverty rate 5-17 (SAIPE), difference in levels | validation | -0.0360 | 0.0055 | [-0.0472, -0.0260] | 620 |
| `seda_growth_gap` | SEDA learning rate (`learn_rate`), difference | validation | 0.0071 | 0.0014 | [0.0043, 0.0099] | 613 |
| `score_gap` | SEDA mean score (`score_mean`), difference | validation | 0.1994 | 0.0225 | [0.1570, 0.2442] | 615 |
| `seda_ses` | SEDA composite SES, difference | validation | 0.5735 | 0.0577 | [0.4648, 0.6887] | 615 |
| `seda_lninc` | log median household income (SEDA), difference | validation | 0.1573 | 0.0184 | [0.1224, 0.1952] | 615 |
| `rev_pp_local` | local revenue per pupil | diagnostic | 0.5636 | 0.0170 | [0.5305, 0.5969] | 621 |
| `rev_pp_state` | state revenue per pupil | diagnostic | -0.3193 | 0.0304 | [-0.3805, -0.2616] | 621 |
| `rev_pp_fed` | federal revenue per pupil | diagnostic | -0.2137 | 0.0431 | [-0.3028, -0.1351] | 621 |
| `rev_pp_total` | total revenue per pupil | diagnostic | 0.1143 | 0.0084 | [0.0968, 0.1303] | 621 |
| `exp_pp` | current expenditure per pupil | diagnostic | 0.0759 | 0.0083 | [0.0594, 0.0924] | 621 |
| `payroll_pp_k12` | audited teacher wage bill per K-12 pupil | diagnostic | 0.0785 | 0.0118 | [0.0548, 0.1016] | 621 |
| `share_local` | local share of revenue, difference in levels | diagnostic | 0.1926 | 0.0077 | [0.1774, 0.2081] | 621 |
| `share_state` | state share of revenue, difference in levels | diagnostic | -0.1711 | 0.0074 | [-0.1859, -0.1570] | 621 |
| `local_revenue_share` | local share of revenue, LEVEL (pupil-weighted over the partition's districts) | diagnostic | 0.4491 | 0.0144 | [0.4212, 0.4765] | 621 |
| `salary_tracked` | tracked `salary_real` gap (zero-coded wage bills kept, NYC absent, pre-K in FTE) | memo (tracked definition) | 0.0154 | 0.0099 | [-0.0039, 0.0350] | 593 |
| `teachers_pp_tracked` | tracked FTE-per-pupil gap (pre-K teachers and pupils included) | memo (tracked definition) | 0.0366 | 0.0060 | [0.0245, 0.0478] | 621 |
| `usd_rev_pp_local` | local revenue, $ per pupil, 2 minus 1 | $ per pupil | 3960 | 265 | [3464, 4472] | 621 |
| `usd_rev_pp_state` | state revenue, $ per pupil, 2 minus 1 | $ per pupil | -1989 | 205 | [-2414, -1630] | 621 |
| `usd_rev_pp_fed` | federal revenue, $ per pupil, 2 minus 1 | $ per pupil | -224 | 43 | [-311, -142] | 621 |
| `usd_rev_pp_total` | total revenue, $ per pupil, 2 minus 1 | $ per pupil | 1747 | 164 | [1416, 2061] | 621 |
| `usd_exp_pp` | current expenditure, $ per pupil, 2 minus 1 | $ per pupil | 1054 | 142 | [773, 1335] | 621 |
| `usd_payroll_k12` | teacher payroll (K-12 pupils), $ per pupil, 2 minus 1 | $ per pupil | 384 | 66 | [258, 514] | 621 |

### ses

ses split (below / above the CZ's pupil-weighted median SEDA SES): 619 CZs, 8,869 districts, 43.45 m pupils.

| moment | what | role | value | boot SE | 95% interval | n CZs |
|---|---|---|---|---|---|---|
| `gap_pupils` | enrollment, log(M2/M1) (CCD enrollment incl. pre-K) | target | -0.0013 | 0.0235 | [-0.0499, 0.0428] | 619 |
| `gap_salary` | CWIFT-deflated teacher salary per K-12 FTE, audit variant 28 | target | 0.0089 | 0.0061 | [-0.0028, 0.0211] | 615 |
| `gap_teachers_pp` | K-12 FTE teachers per K-12 pupil, audit FTE variant 11 | target | -0.0070 | 0.0063 | [-0.0186, 0.0056] | 618 |
| `gap_income` | log pooled median earnings, full-time year-round workers (`pm_ftyr_all`) | target | 0.2317 | 0.0125 | [0.2072, 0.2556] | 618 |
| `pm_earn_all` | log pooled median earnings, all workers with earnings (`pm_earn_all`) | target | 0.2392 | 0.0110 | [0.2169, 0.2603] | 618 |
| `poverty_gap` | child poverty rate 5-17 (SAIPE), difference in levels | validation | -0.1018 | 0.0043 | [-0.1102, -0.0929] | 619 |
| `seda_growth_gap` | SEDA learning rate (`learn_rate`), difference | validation | 0.0145 | 0.0017 | [0.0111, 0.0177] | 617 |
| `score_gap` | SEDA mean score (`score_mean`), difference | validation | 0.4047 | 0.0194 | [0.3654, 0.4406] | 619 |
| `seda_ses` | SEDA composite SES, difference | validation | 1.1353 | 0.0450 | [1.0467, 1.2251] | 619 |
| `seda_lninc` | log median household income (SEDA), difference | validation | 0.3542 | 0.0152 | [0.3245, 0.3833] | 619 |
| `rev_pp_local` | local revenue per pupil | diagnostic | 0.2321 | 0.0235 | [0.1878, 0.2789] | 619 |
| `rev_pp_state` | state revenue per pupil | diagnostic | -0.2573 | 0.0274 | [-0.3110, -0.2036] | 619 |
| `rev_pp_fed` | federal revenue per pupil | diagnostic | -0.6572 | 0.0341 | [-0.7255, -0.5934] | 619 |
| `rev_pp_total` | total revenue per pupil | diagnostic | -0.0455 | 0.0113 | [-0.0682, -0.0248] | 619 |
| `exp_pp` | current expenditure per pupil | diagnostic | -0.0509 | 0.0106 | [-0.0713, -0.0299] | 619 |
| `payroll_pp_k12` | audited teacher wage bill per K-12 pupil | diagnostic | 0.0136 | 0.0086 | [-0.0035, 0.0310] | 619 |
| `share_local` | local share of revenue, difference in levels | diagnostic | 0.1212 | 0.0105 | [0.1012, 0.1428] | 619 |
| `share_state` | state share of revenue, difference in levels | diagnostic | -0.0808 | 0.0101 | [-0.1017, -0.0613] | 619 |
| `local_revenue_share` | local share of revenue, LEVEL (pupil-weighted over the partition's districts) | diagnostic | 0.4474 | 0.0139 | [0.4198, 0.4742] | 619 |
| `salary_tracked` | tracked `salary_real` gap (zero-coded wage bills kept, NYC absent, pre-K in FTE) | memo (tracked definition) | 0.0187 | 0.0062 | [0.0058, 0.0305] | 591 |
| `teachers_pp_tracked` | tracked FTE-per-pupil gap (pre-K teachers and pupils included) | memo (tracked definition) | -0.0036 | 0.0068 | [-0.0162, 0.0102] | 618 |
| `usd_rev_pp_local` | local revenue, $ per pupil, 2 minus 1 | $ per pupil | 1824 | 274 | [1323, 2387] | 619 |
| `usd_rev_pp_state` | state revenue, $ per pupil, 2 minus 1 | $ per pupil | -1812 | 248 | [-2305, -1340] | 619 |
| `usd_rev_pp_fed` | federal revenue, $ per pupil, 2 minus 1 | $ per pupil | -655 | 35 | [-724, -587] | 619 |
| `usd_rev_pp_total` | total revenue, $ per pupil, 2 minus 1 | $ per pupil | -642 | 193 | [-1031, -289] | 619 |
| `usd_exp_pp` | current expenditure, $ per pupil, 2 minus 1 | $ per pupil | -521 | 157 | [-819, -217] | 619 |
| `usd_payroll_k12` | teacher payroll (K-12 pupils), $ per pupil, 2 minus 1 | $ per pupil | 100 | 49 | [7, 201] | 619 |

## 5. The same CZs

All three settings restricted to the 186 CZs in every partition (bootstrap over those CZs, CZ weights unchanged). `full` columns repeat section 4.

| moment | locale, full | locale, common | revenue, full | revenue, common | ses, full | ses, common |
|---|---|---|---|---|---|---|
| `gap_pupils` | 0.3764 (0.1154) | 0.3772 (0.1194) | -0.0247 (0.0218) | -0.0029 (0.0215) | -0.0013 (0.0235) | 0.0104 (0.0203) |
| `gap_salary` | 0.0086 (0.0092) | 0.0087 (0.0092) | 0.0178 (0.0098) | 0.0211 (0.0126) | 0.0089 (0.0061) | 0.0136 (0.0072) |
| `gap_teachers_pp` | -0.0177 (0.0131) | -0.0169 (0.0127) | 0.0348 (0.0057) | 0.0356 (0.0068) | -0.0070 (0.0063) | -0.0088 (0.0074) |
| `gap_income` | 0.1122 (0.0168) | 0.1123 (0.0169) | 0.1286 (0.0114) | 0.1516 (0.0140) | 0.2317 (0.0125) | 0.2725 (0.0135) |
| `pm_earn_all` | 0.1435 (0.0196) | 0.1440 (0.0197) | 0.1000 (0.0117) | 0.1243 (0.0142) | 0.2392 (0.0110) | 0.2769 (0.0121) |
| `poverty_gap` | -0.0888 (0.0090) | -0.0890 (0.0094) | -0.0360 (0.0055) | -0.0420 (0.0072) | -0.1018 (0.0043) | -0.1126 (0.0050) |
| `seda_growth_gap` | 0.0031 (0.0030) | 0.0030 (0.0030) | 0.0071 (0.0014) | 0.0080 (0.0018) | 0.0145 (0.0017) | 0.0158 (0.0022) |
| `score_gap` | 0.2887 (0.0326) | 0.2892 (0.0337) | 0.1994 (0.0225) | 0.2391 (0.0276) | 0.4047 (0.0194) | 0.4726 (0.0201) |
| `seda_ses` | 0.8299 (0.1013) | 0.8314 (0.1052) | 0.5735 (0.0577) | 0.6670 (0.0728) | 1.1353 (0.0450) | 1.2746 (0.0506) |
| `seda_lninc` | 0.2708 (0.0268) | 0.2712 (0.0278) | 0.1573 (0.0184) | 0.1850 (0.0237) | 0.3542 (0.0152) | 0.4042 (0.0163) |
| `rev_pp_local` | -0.0708 (0.0558) | -0.0689 (0.0560) | 0.5636 (0.0170) | 0.5874 (0.0210) | 0.2321 (0.0235) | 0.2574 (0.0284) |
| `rev_pp_state` | -0.0003 (0.0678) | -0.0001 (0.0666) | -0.3193 (0.0304) | -0.3501 (0.0392) | -0.2573 (0.0274) | -0.2943 (0.0359) |
| `rev_pp_fed` | -0.5677 (0.0546) | -0.5677 (0.0560) | -0.2137 (0.0431) | -0.2483 (0.0572) | -0.6572 (0.0341) | -0.7603 (0.0386) |
| `rev_pp_total` | -0.1211 (0.0163) | -0.1209 (0.0163) | 0.1143 (0.0084) | 0.1184 (0.0105) | -0.0455 (0.0113) | -0.0506 (0.0138) |
| `exp_pp` | -0.0956 (0.0157) | -0.0952 (0.0155) | 0.0759 (0.0083) | 0.0793 (0.0105) | -0.0509 (0.0106) | -0.0569 (0.0128) |
| `payroll_pp_k12` | -0.0267 (0.0201) | -0.0258 (0.0197) | 0.0785 (0.0118) | 0.0841 (0.0152) | 0.0136 (0.0086) | 0.0131 (0.0107) |
| `share_local` | 0.0105 (0.0214) | 0.0107 (0.0216) | 0.1926 (0.0077) | 0.2056 (0.0099) | 0.1212 (0.0105) | 0.1379 (0.0132) |
| `share_state` | 0.0196 (0.0203) | 0.0195 (0.0204) | -0.1711 (0.0074) | -0.1835 (0.0095) | -0.0808 (0.0101) | -0.0933 (0.0131) |
| `local_revenue_share` | 0.4726 (0.0206) | 0.4730 (0.0201) | 0.4491 (0.0144) | 0.4614 (0.0179) | 0.4474 (0.0139) | 0.4592 (0.0176) |
| `salary_tracked` | 0.0162 (0.0090) | 0.0162 (0.0090) | 0.0154 (0.0099) | 0.0128 (0.0123) | 0.0187 (0.0062) | 0.0234 (0.0071) |
| `teachers_pp_tracked` | -0.0122 (0.0139) | -0.0113 (0.0136) | 0.0366 (0.0060) | 0.0376 (0.0074) | -0.0036 (0.0068) | -0.0047 (0.0083) |

## 6. Bootstrap correlation of the four targets

**locale**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.12 | -0.37 | 0.28 |
| `gap_salary` | -0.12 | 1.00 | 0.41 | 0.15 |
| `gap_teachers_pp` | -0.37 | 0.41 | 1.00 | -0.07 |
| `gap_income` | 0.28 | 0.15 | -0.07 | 1.00 |

**revenue**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.04 | -0.08 | 0.17 |
| `gap_salary` | -0.04 | 1.00 | 0.24 | 0.06 |
| `gap_teachers_pp` | -0.08 | 0.24 | 1.00 | 0.23 |
| `gap_income` | 0.17 | 0.06 | 0.23 | 1.00 |

**ses**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.07 | -0.14 | 0.07 |
| `gap_salary` | -0.07 | 1.00 | 0.06 | 0.31 |
| `gap_teachers_pp` | -0.14 | 0.06 | 1.00 | 0.05 |
| `gap_income` | 0.07 | 0.31 | 0.05 | 1.00 |

**locale_common**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.09 | -0.37 | 0.31 |
| `gap_salary` | -0.09 | 1.00 | 0.42 | 0.17 |
| `gap_teachers_pp` | -0.37 | 0.42 | 1.00 | -0.04 |
| `gap_income` | 0.31 | 0.17 | -0.04 | 1.00 |

**revenue_common**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.08 | 0.07 | 0.05 |
| `gap_salary` | -0.08 | 1.00 | 0.32 | 0.07 |
| `gap_teachers_pp` | 0.07 | 0.32 | 1.00 | 0.29 |
| `gap_income` | 0.05 | 0.07 | 0.29 | 1.00 |

**ses_common**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.05 | 0.05 | -0.09 |
| `gap_salary` | -0.05 | 1.00 | 0.15 | 0.25 |
| `gap_teachers_pp` | 0.05 | 0.15 | 1.00 | 0.24 |
| `gap_income` | -0.09 | 0.25 | 0.24 | 1.00 |

The correlation among the targets and the validation moments (every statistic, draw by draw) is in `district_settings.json` (`corr`) and the draws in `district_settings_draws_<set>.csv`.

## 7. State-aid offsets under the revenue split

### 7a. Decomposition of the gaps

Within-CZ gap between location 2 and location 1 (above- and below-median local revenue per pupil for the revenue split; suburb and city, high and low SES for the others) of each component, in logs (ratio of location totals) and in dollars per pupil (location 2 minus 1, CZ-weighted); the dollar gaps add up (total = local + state + federal). `per $ of local gap` is the component's dollar gap over the local dollar gap, so a negative number is the offset (n.m. where the local gap is within two SEs of zero).

| set | component | log gap | boot SE | $ per pupil | boot SE | per $ of local gap | boot SE |
|---|---|---|---|---|---|---|---|
| revenue | local revenue | 0.564 | 0.017 | 3,960 | 265 | 1 |  |
| revenue | state revenue | -0.319 | 0.030 | -1,989 | 205 | -0.502 | 0.033 |
| revenue | federal revenue | -0.214 | 0.043 | -224 | 43 | -0.056 | 0.010 |
| revenue | total revenue | 0.114 | 0.008 | 1,747 | 164 | 0.441 | 0.036 |
| revenue | current expenditure | 0.076 | 0.008 | 1,054 | 142 | 0.266 | 0.026 |
| revenue | teacher payroll, K-12 pupils | 0.078 | 0.012 | 384 | 66 | 0.097 | 0.013 |
| revenue_common | local revenue | 0.587 | 0.021 | 4,362 | 341 | 1 |  |
| revenue_common | state revenue | -0.350 | 0.039 | -2,252 | 272 | -0.516 | 0.038 |
| revenue_common | federal revenue | -0.248 | 0.057 | -245 | 57 | -0.056 | 0.012 |
| revenue_common | total revenue | 0.118 | 0.010 | 1,864 | 208 | 0.427 | 0.042 |
| revenue_common | current expenditure | 0.079 | 0.010 | 1,149 | 182 | 0.263 | 0.030 |
| revenue_common | teacher payroll, K-12 pupils | 0.084 | 0.015 | 431 | 85 | 0.099 | 0.016 |
| locale | local revenue | -0.071 | 0.056 | -768 | 622 | 1 |  |
| locale | state revenue | -0.000 | 0.068 | -772 | 517 | n.m. |  |
| locale | federal revenue | -0.568 | 0.055 | -643 | 67 | n.m. |  |
| locale | total revenue | -0.121 | 0.016 | -2,184 | 346 | n.m. |  |
| locale | current expenditure | -0.096 | 0.016 | -1,301 | 258 | n.m. |  |
| locale | teacher payroll, K-12 pupils | -0.027 | 0.020 | -80 | 131 | n.m. |  |
| ses | local revenue | 0.232 | 0.024 | 1,824 | 274 | 1 |  |
| ses | state revenue | -0.257 | 0.027 | -1,812 | 248 | -0.993 | 0.102 |
| ses | federal revenue | -0.657 | 0.034 | -655 | 35 | -0.359 | 0.053 |
| ses | total revenue | -0.045 | 0.011 | -642 | 193 | -0.352 | 0.139 |
| ses | current expenditure | -0.051 | 0.011 | -521 | 157 | -0.285 | 0.121 |
| ses | teacher payroll, K-12 pupils | 0.014 | 0.009 | 100 | 49 | 0.055 | 0.023 |

### 7b. Within-CZ district regressions

Each aggregate per pupil on local revenue per pupil, district level, CZ fixed effects, pupil weights, CZ-clustered SEs (`data_estimate.ols_fe`); per-pupil ratios screened at the panel's [0.5%, 99.5%] within-year tails (state revenue and local revenue as `data_estimate` screens them; federal revenue and payroll screened here the same way). The slope in dollars is dollars of the left-hand side per extra dollar of local revenue per pupil: the state-on-local rows are the offset rate.

| sample | regression | slope | SE | districts | CZs | within R2 |
|---|---|---|---|---|---|---|
| revenue | state on local, $ per pupil | -0.466 | 0.053 | 9,206 | 621 | 0.38 |
| revenue | state on local, logs | -0.537 | 0.035 | 9,206 | 621 | 0.43 |
| revenue | federal on local, $ per pupil | -0.044 | 0.008 | 9,233 | 621 | 0.05 |
| revenue | total revenue on local, $ per pupil | 0.510 | 0.058 | 9,197 | 621 | 0.34 |
| revenue | total revenue on local, logs | 0.166 | 0.021 | 9,197 | 621 | 0.22 |
| revenue | current expenditure on local, $ per pupil | 0.310 | 0.032 | 9,201 | 621 | 0.26 |
| revenue | current expenditure on local, logs | 0.107 | 0.015 | 9,201 | 621 | 0.12 |
| revenue | teacher payroll per K-12 pupil on local, $ | 0.109 | 0.013 | 9,201 | 621 | 0.30 |
| revenue | teacher payroll per K-12 pupil on local, logs | 0.116 | 0.017 | 9,201 | 621 | 0.14 |
| revenue_common | state on local, $ per pupil | -0.481 | 0.059 | 6,164 | 186 | 0.39 |
| revenue_common | state on local, logs | -0.561 | 0.041 | 6,164 | 186 | 0.44 |
| revenue_common | federal on local, $ per pupil | -0.046 | 0.009 | 6,177 | 186 | 0.06 |
| revenue_common | total revenue on local, $ per pupil | 0.493 | 0.065 | 6,152 | 186 | 0.32 |
| revenue_common | total revenue on local, logs | 0.157 | 0.024 | 6,152 | 186 | 0.20 |
| revenue_common | current expenditure on local, $ per pupil | 0.316 | 0.036 | 6,152 | 186 | 0.27 |
| revenue_common | current expenditure on local, logs | 0.102 | 0.018 | 6,152 | 186 | 0.11 |
| revenue_common | teacher payroll per K-12 pupil on local, $ | 0.114 | 0.015 | 6,158 | 186 | 0.32 |
| revenue_common | teacher payroll per K-12 pupil on local, logs | 0.117 | 0.020 | 6,158 | 186 | 0.14 |
| locale | state on local, $ per pupil | -0.494 | 0.067 | 3,222 | 188 | 0.39 |
| locale | state on local, logs | -0.594 | 0.046 | 3,222 | 188 | 0.47 |
| locale | federal on local, $ per pupil | -0.046 | 0.011 | 3,229 | 188 | 0.07 |
| locale | total revenue on local, $ per pupil | 0.484 | 0.074 | 3,220 | 187 | 0.31 |
| locale | total revenue on local, logs | 0.150 | 0.028 | 3,220 | 187 | 0.17 |
| locale | current expenditure on local, $ per pupil | 0.321 | 0.040 | 3,220 | 187 | 0.28 |
| locale | current expenditure on local, logs | 0.103 | 0.021 | 3,220 | 187 | 0.10 |
| locale | teacher payroll per K-12 pupil on local, $ | 0.116 | 0.017 | 3,218 | 188 | 0.34 |
| locale | teacher payroll per K-12 pupil on local, logs | 0.118 | 0.024 | 3,218 | 188 | 0.15 |
| ses | state on local, $ per pupil | -0.527 | 0.049 | 8,828 | 619 | 0.45 |
| ses | state on local, logs | -0.559 | 0.035 | 8,828 | 619 | 0.46 |
| ses | federal on local, $ per pupil | -0.046 | 0.009 | 8,852 | 619 | 0.05 |
| ses | total revenue on local, $ per pupil | 0.440 | 0.053 | 8,814 | 619 | 0.29 |
| ses | total revenue on local, logs | 0.155 | 0.022 | 8,814 | 619 | 0.20 |
| ses | current expenditure on local, $ per pupil | 0.263 | 0.030 | 8,819 | 619 | 0.21 |
| ses | current expenditure on local, logs | 0.095 | 0.015 | 8,819 | 619 | 0.10 |
| ses | teacher payroll per K-12 pupil on local, $ | 0.103 | 0.013 | 8,817 | 619 | 0.26 |
| ses | teacher payroll per K-12 pupil on local, logs | 0.108 | 0.016 | 8,817 | 619 | 0.13 |

## 8. Pre-registered hold-out (locale)

Rule (fixed before any hold-out number was computed): sort the 188 locale CZs by partitioned pupils (descending, ties by ascending CZ id; 0 tied pairs), pair consecutive CZs (1-2, 3-4, ...; an unpaired last CZ goes to half A), and for each pair take one `rng.integers(2)` from `numpy.random.default_rng(20261002)`, in pair order: 0 sends the first CZ of the pair to half B, 1 sends the second; the other CZ goes to half A. Half A has 94 CZs and half B 94. **Half B's results are written only to `sealed/holdout_B.json` and are not printed or summarised in any log or report.** The weights of the two halves are public, so B's point estimates are implied by the full-sample and A values; the seal is procedural.

**Half A**

| moment | what | role | value | boot SE | 95% interval | n CZs |
|---|---|---|---|---|---|---|
| `gap_pupils` | enrollment, log(M2/M1) (CCD enrollment incl. pre-K) | target | 0.5220 | 0.1515 | [0.2491, 0.8492] | 94 |
| `gap_salary` | CWIFT-deflated teacher salary per K-12 FTE, audit variant 28 | target | 0.0090 | 0.0138 | [-0.0178, 0.0365] | 94 |
| `gap_teachers_pp` | K-12 FTE teachers per K-12 pupil, audit FTE variant 11 | target | -0.0312 | 0.0149 | [-0.0626, -0.0039] | 94 |
| `gap_income` | log pooled median earnings, full-time year-round workers (`pm_ftyr_all`) | target | 0.1183 | 0.0234 | [0.0793, 0.1700] | 94 |
| `pm_earn_all` | log pooled median earnings, all workers with earnings (`pm_earn_all`) | target | 0.1513 | 0.0312 | [0.0982, 0.2198] | 94 |
| `poverty_gap` | child poverty rate 5-17 (SAIPE), difference in levels | validation | -0.0992 | 0.0155 | [-0.1319, -0.0731] | 94 |
| `seda_growth_gap` | SEDA learning rate (`learn_rate`), difference | validation | 0.0069 | 0.0030 | [0.0015, 0.0135] | 94 |
| `score_gap` | SEDA mean score (`score_mean`), difference | validation | 0.2947 | 0.0563 | [0.2019, 0.4199] | 94 |
| `seda_ses` | SEDA composite SES, difference | validation | 0.8869 | 0.1649 | [0.6109, 1.2447] | 94 |
| `seda_lninc` | log median household income (SEDA), difference | validation | 0.2872 | 0.0440 | [0.2134, 0.3825] | 94 |
| `rev_pp_local` | local revenue per pupil | diagnostic | -0.1426 | 0.0561 | [-0.2492, -0.0309] | 94 |
| `rev_pp_state` | state revenue per pupil | diagnostic | 0.0111 | 0.1070 | [-0.1760, 0.2460] | 94 |
| `rev_pp_fed` | federal revenue per pupil | diagnostic | -0.6111 | 0.0998 | [-0.8208, -0.4360] | 94 |
| `rev_pp_total` | total revenue per pupil | diagnostic | -0.1479 | 0.0224 | [-0.1976, -0.1129] | 94 |
| `exp_pp` | current expenditure per pupil | diagnostic | -0.1144 | 0.0238 | [-0.1632, -0.0687] | 94 |
| `payroll_pp_k12` | audited teacher wage bill per K-12 pupil | diagnostic | -0.0425 | 0.0281 | [-0.1015, 0.0084] | 94 |
| `share_local` | local share of revenue, difference in levels | diagnostic | -0.0024 | 0.0270 | [-0.0537, 0.0523] | 94 |
| `share_state` | state share of revenue, difference in levels | diagnostic | 0.0349 | 0.0243 | [-0.0129, 0.0837] | 94 |
| `local_revenue_share` | local share of revenue, LEVEL (pupil-weighted over the partition's districts) | diagnostic | 0.4337 | 0.0292 | [0.3816, 0.4953] | 94 |
| `salary_tracked` | tracked `salary_real` gap (zero-coded wage bills kept, NYC absent, pre-K in FTE) | memo (tracked definition) | 0.0229 | 0.0136 | [-0.0019, 0.0514] | 90 |
| `teachers_pp_tracked` | tracked FTE-per-pupil gap (pre-K teachers and pupils included) | memo (tracked definition) | -0.0316 | 0.0152 | [-0.0637, -0.0037] | 94 |
| `usd_rev_pp_local` | local revenue, $ per pupil, 2 minus 1 | $ per pupil | -1162 | 779 | [-2863, 137] | 94 |
| `usd_rev_pp_state` | state revenue, $ per pupil, 2 minus 1 | $ per pupil | -653 | 525 | [-1704, 346] | 94 |
| `usd_rev_pp_fed` | federal revenue, $ per pupil, 2 minus 1 | $ per pupil | -709 | 117 | [-962, -511] | 94 |
| `usd_rev_pp_total` | total revenue, $ per pupil, 2 minus 1 | $ per pupil | -2524 | 488 | [-3610, -1778] | 94 |
| `usd_exp_pp` | current expenditure, $ per pupil, 2 minus 1 | $ per pupil | -1487 | 385 | [-2300, -760] | 94 |
| `usd_payroll_k12` | teacher payroll (K-12 pupils), $ per pupil, 2 minus 1 | $ per pupil | -196 | 159 | [-538, 89] | 94 |

**Bootstrap correlation of the four targets, half A**

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.17 | -0.35 | -0.03 |
| `gap_salary` | -0.17 | 1.00 | 0.64 | 0.56 |
| `gap_teachers_pp` | -0.35 | 0.64 | 1.00 | 0.35 |
| `gap_income` | -0.03 | 0.56 | 0.35 | 1.00 |

**CZ assignment** (CZ id: half), in pair order of partitioned pupils:

```
38300:A 19400:B  24300:B 32000:A  19600:B 11304:A  33100:B 19700:A  20500:B 39400:A  37800:B 11600:A
35001:B 28900:A  37400:A 20901:B  38000:B 21501:A  11302:A 33000:B  37500:B 31301:A  36100:A 38801:B
31600:A 29502:B  15200:A 31201:B  900:B 14200:A  12701:A 7100:B  37200:A 16300:B  15900:A 18000:B
1701:B 24100:A  5600:B 20401:A  5202:B 33803:A  7600:A 2000:B  2400:A 13101:B  28202:A 19100:B
37100:B 10700:A  36000:A 3500:B  12200:A 37000:B  34901:A 28401:B  19200:A 30402:B  8202:A 6900:B
35100:B 8100:A  12501:A 18600:B  20600:A 38901:B  38200:A 19000:B  35801:A 27501:B  8401:B 17700:A
29301:B 38601:A  7500:B 28800:A  2500:B 32900:A  4200:A 20800:B  1302:B 13501:A  23100:B 6000:A
18800:A 13600:B  24000:A 24400:B  39100:B 1400:A  6401:B 22500:A  17400:B 37700:A  38100:B 32100:A
100:B 16400:A  23801:A 14100:B  15000:B 29700:A  12100:A 3800:B  16600:B 3003:A  700:B 11101:A
11700:B 9701:A  22601:A 30802:B  20100:B 11900:A  19901:A 36200:B  15100:A 36800:B  36301:A 37300:B
3400:B 2700:A  23900:A 17100:B  14700:B 21400:A  22200:B 26801:A  32801:B 9600:A  39000:B 10801:A
29204:A 28502:B  8000:B 3901:A  23500:B 13700:A  6600:A 30100:B  33300:B 9800:A  24802:B 16500:A
28001:A 17900:B  14000:A 34115:B  16200:A 35401:B  33400:B 23400:A  14500:B 22800:A  13000:B 26201:A
13200:A 39302:B  22900:A 18100:B  33500:B 22700:A  29901:A 36600:B  6100:A 26002:B  15600:B 38502:A
20200:B 14600:A  33902:A 22400:B  25900:A 24200:B  38700:A 11500:B  25601:A 16000:B  20001:B 25701:A
18900:B 6200:A  15500:A 26704:B  38401:B 2900:A  13400:B 18300:A
```

## 9. Sensitivities (not part of the targets)

### 9a. NYC as one fully observed district

The headline panel collapses NYC's 32 geographic districts into one row but keeps, for the revenue, SAIPE, SEDA and ACS moments, the geographic districts' missing values, because those records are filed under the Chancellor's Office id (`3620580`, agency type 3): F-33 revenue and expenditure, SAIPE, SEDA (`New York City Public Schools`) and the ACS school-district tables all have NYC in one record. This panel attaches all of them to the collapsed NYC row (per-pupil ratios screened with the panel's own cuts). The enrollment, FTE and salary moments are unchanged under locale; the others move as shown.

NYC's attached values: 960,484 CCD pupils; local revenue $20,798, state $12,341 and total revenue $35,179 per pupil; SAIPE 1,204,282 children 5-17, 23.1% poor; SEDA SES -0.376, growth 0.0148, score -0.0016; median earnings (B20002) $38,845. It sits in location 1 (locale), 2 (revenue), 1 (ses).

| setting | moment | headline | NYC filled | difference |
|---|---|---|---|---|
| locale | `gap_pupils` | 0.3764 (0.1154) | 0.3764 (0.1154) | +0.0000 |
| locale | `gap_salary` | 0.0086 (0.0092) | 0.0086 (0.0092) | +0.0000 |
| locale | `gap_teachers_pp` | -0.0177 (0.0131) | -0.0177 (0.0131) | +0.0000 |
| locale | `gap_income` | 0.1122 (0.0168) | 0.1234 (0.0164) | +0.0111 |
| locale | `poverty_gap` | -0.0888 (0.0090) | -0.0951 (0.0091) | -0.0063 |
| locale | `seda_growth_gap` | 0.0031 (0.0030) | 0.0037 (0.0026) | +0.0006 |
| locale | `score_gap` | 0.2887 (0.0326) | 0.2889 (0.0326) | +0.0001 |
| locale | `pm_earn_all` | 0.1435 (0.0196) | 0.1507 (0.0190) | +0.0072 |
| locale | `rev_pp_local` | -0.0708 (0.0558) | -0.0594 (0.0548) | +0.0114 |
| locale | `rev_pp_state` | -0.0003 (0.0678) | -0.0363 (0.0666) | -0.0360 |
| locale | `rev_pp_total` | -0.1211 (0.0163) | -0.1244 (0.0165) | -0.0033 |
| revenue | `gap_pupils` | -0.0247 (0.0218) | 0.0198 (0.0467) | +0.0445 |
| revenue | `gap_salary` | 0.0178 (0.0098) | 0.0172 (0.0096) | -0.0007 |
| revenue | `gap_teachers_pp` | 0.0348 (0.0057) | 0.0308 (0.0058) | -0.0041 |
| revenue | `gap_income` | 0.1286 (0.0114) | 0.1205 (0.0122) | -0.0082 |
| revenue | `poverty_gap` | -0.0360 (0.0055) | -0.0313 (0.0068) | +0.0047 |
| revenue | `seda_growth_gap` | 0.0071 (0.0014) | 0.0077 (0.0016) | +0.0006 |
| revenue | `score_gap` | 0.1994 (0.0225) | 0.1925 (0.0219) | -0.0069 |
| revenue | `pm_earn_all` | 0.1000 (0.0117) | 0.0936 (0.0122) | -0.0064 |
| revenue | `rev_pp_local` | 0.5636 (0.0170) | 0.5574 (0.0178) | -0.0062 |
| revenue | `rev_pp_state` | -0.3193 (0.0304) | -0.3061 (0.0307) | +0.0131 |
| revenue | `rev_pp_total` | 0.1143 (0.0084) | 0.1176 (0.0090) | +0.0032 |
| ses | `gap_pupils` | -0.0013 (0.0235) | -0.0204 (0.0294) | -0.0191 |
| ses | `gap_salary` | 0.0089 (0.0061) | 0.0091 (0.0060) | +0.0001 |
| ses | `gap_teachers_pp` | -0.0070 (0.0063) | -0.0032 (0.0084) | +0.0038 |
| ses | `gap_income` | 0.2317 (0.0125) | 0.2294 (0.0121) | -0.0024 |
| ses | `poverty_gap` | -0.1018 (0.0043) | -0.1039 (0.0045) | -0.0022 |
| ses | `seda_growth_gap` | 0.0145 (0.0017) | 0.0133 (0.0020) | -0.0013 |
| ses | `score_gap` | 0.4047 (0.0194) | 0.3966 (0.0200) | -0.0081 |
| ses | `pm_earn_all` | 0.2392 (0.0110) | 0.2366 (0.0110) | -0.0026 |
| ses | `rev_pp_local` | 0.2321 (0.0235) | 0.2199 (0.0244) | -0.0122 |
| ses | `rev_pp_state` | -0.2573 (0.0274) | -0.2557 (0.0266) | +0.0016 |
| ses | `rev_pp_total` | -0.0455 (0.0113) | -0.0519 (0.0113) | -0.0064 |

NYC pupils in the partition with NYC filled: locale 960,484, revenue 960,484, ses 960,484.

### 9b. Revenue split on local revenue per K-12 pupil

The audited FTE and payroll moments use K-12 pupils, but the revenue split divides local revenue by CCD enrollment, which includes pre-K pupils. This variant splits on local revenue per K-12 pupil (`pupils - pre-K enrollment`), screened at the same tails; moments unchanged. 1.5% of the pupils in districts that are in both partitions (170 of 9,243 districts) switch location; the new partition has 621 CZs and 9,244 districts (baseline 621 and 9,245).

| moment | per CCD pupil (headline) | per K-12 pupil | difference |
|---|---|---|---|
| `gap_pupils` | -0.0247 (0.0218) | -0.0306 (0.0225) | -0.0059 |
| `gap_salary` | 0.0178 (0.0098) | 0.0184 (0.0099) | +0.0005 |
| `gap_teachers_pp` | 0.0348 (0.0057) | 0.0346 (0.0057) | -0.0002 |
| `gap_income` | 0.1286 (0.0114) | 0.1252 (0.0115) | -0.0034 |
| `poverty_gap` | -0.0360 (0.0055) | -0.0344 (0.0055) | +0.0016 |
| `seda_growth_gap` | 0.0071 (0.0014) | 0.0069 (0.0014) | -0.0002 |
| `score_gap` | 0.1994 (0.0225) | 0.1953 (0.0226) | -0.0041 |
| `pm_earn_all` | 0.1000 (0.0117) | 0.0968 (0.0118) | -0.0032 |
| `rev_pp_local` | 0.5636 (0.0170) | 0.5631 (0.0171) | -0.0006 |
| `rev_pp_state` | -0.3193 (0.0304) | -0.3197 (0.0306) | -0.0004 |
| `rev_pp_total` | 0.1143 (0.0084) | 0.1148 (0.0084) | +0.0005 |

## 10. Judgment calls and caveats

- **NYC.** F-33 revenue and expenditure, SAIPE, SEDA and the ACS school-district tables file New York City under the Chancellor's Office id (`3620580`), which the frame filter drops; the CCD frame has its 32 geographic districts (960,484 pupils), whose values for those sources are missing. Locale: NYC enters `gap_pupils`, `gap_teachers_pp` and, through the audit's collapse to one row with the Chancellor's wage bill and CWIFT, `gap_salary`, but not `gap_income`, the poverty and SEDA gaps or the revenue moments, so the New York CZ's city location in those is the CZ's other cities. Revenue and ses: the split variable (local revenue per pupil, SEDA SES) is missing for NYC with or without the collapse, so NYC is outside both partitions, the collapse changes nothing there, and the New York CZ is split over its other districts. Section 9a attaches NYC's own records: it moves the locale `gap_income` by +0.0111 and `poverty_gap` by -0.0063, and puts NYC into the revenue and ses partitions, moving the revenue-split `gap_pupils` by +0.0445 (to 0.0198, SE 0.0467) and `gap_teachers_pp` by -0.0041. The headline numbers follow the instruction (tracked sample rules, the audit's corrections) and are not NYC-complete.
- **Pre-K corrections.** They are district-level data corrections (pre-K FTE out of the salary denominator, K-12 FTE over K-12 pupils) and apply unchanged inside each location under every split. The partition itself, the enrollment gap and the per-pupil revenue and expenditure moments still divide by CCD enrollment, which includes pre-K pupils, so the revenue split ranks districts by local revenue per CCD pupil while the FTE and payroll moments use K-12 pupils. Section 9b splits on revenue per K-12 pupil instead: 1.5% of pupils change location and the four revenue targets move by less than their bootstrap SEs (largest change 0.0059).
- **Hybrid fill and NYC collapse.** One constant (0.881, the locale sample's median ratio of teacher to instruction salaries) fills the zero-coded wage bills in every setting, as the target is defined; the sample medians are 0.881 (locale), 0.872 (revenue) and 0.872 (ses), NYC excluded. Without the fill and the collapse a tracked salary exists for 89.9% (locale), 94.2% (revenue) and 94.8% (ses) of the partitioned pupils; with them 98.9%, 99.2% and 99.4%. The locale-sample constant is not re-estimated per setting.
- **Common CZs.** The comparison set is the 186 CZs in all three partitions, so that one set serves all three settings. Pairwise with locale the common sets would have 187 (revenue) and 187 (ses) CZs, because two locale CZs are each missing from one other partition: CZ 34115 (Fairbanks North Star Borough School District) from the revenue partition, 0.07% of the locale weight; CZ 20200 (Champlain Valley Unified School District) from the ses partition, 0.05% of the locale weight. Their weights are negligible: on the pairwise sets (`common.pairwise` in the JSON) the seven TOML moments and `pm_earn_all` differ from the three-way ones by at most 0.0009. The CZ weights are each setting's own partitioned pupils, so the common-CZ revenue and ses averages weight a CZ by how many pupils its revenue (ses) partition covers, which is more than its city + suburb pupils.
- **Weights and effective sample.** Weights are partitioned pupils, so a few large metros dominate: the effective number of CZs is 42 (locale), 89 (revenue) and 94 (ses), and 19% of the revenue CZs have three or fewer partitioned districts, where the median split is a comparison of one district with the rest. The revenue and ses locations are not geographic and include town and rural districts: location 1/2 means below/above the CZ median, wherever the districts lie.
- **Mechanical gaps.** Under the revenue split the local-revenue gap, and under the ses split the SES gap, are produced by the split and are not findings; the pupil-weighted median cut also equalises enrollment, so `gap_pupils` is about zero by construction there. They are reported as diagnostics, not targets.
- **Bootstrap scope.** As in `district_bootstrap.py`, only CZ membership of the average is resampled; the district panel, the [0.5%, 99.5%] screens, the revenue / SES cut points and the CZ filters are held fixed. A CZ missing a statistic drops out of that statistic's draws, so n varies by draw. The sets use separate multiplicity matrices (same seed, different n), so correlations are comparable within a set, not across sets.
- **Offsets are associations.** The state-on-local slopes are within-CZ district correlations; state aid formulas respond to local wealth by design, so they are the offset rate the model needs, not the causal effect of exogenous local revenue.
- **Seal.** B's weights and CZ ids are public, so its point estimates are implied by the full-sample and A numbers; the seal only keeps B out of the logs and reports.

## Files

- `district_settings.json`: everything above, machine readable (moments with SE, interval and n; the full bootstrap correlation matrix of every set; coverage; regressions; hold-out assignment and half A).
- `district_settings_targets.toml`: value / se / n_cz of the four targets and three validation moments, and the 4 x 4 target correlation, for each setting, the common-CZ versions and half A.
- `district_settings_draws_<set>.csv`: the bootstrap draws (one row per draw, one column per moment) for `locale`, `revenue`, `ses`, `locale_common`, `revenue_common`, `ses_common`, `locale_holdout_A`.
- `district_settings-<jobid>.out`: the SLURM log of the run; it carries no hold-out numbers.
- `sealed/holdout_B.json`: hold-out half B, written sealed. Nothing in the logs or in this report reads it back; open it only when the researchers decide to.
