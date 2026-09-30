# District moments: audit of pay and FTE counts (T3)

Generated 2026-09-29T16:25:07 by `district_audit.py` (year 2018, `locale` partition, whole-CZ bootstrap B = 1999, seed 20260929). Every number here is recomputed with `data_estimate.py`'s own panel, partition and `czone_gaps`; the baseline rows reproduce the tracked moments exactly.

## Summary

**Salary gap** (`gap_salary_real_locale` = 0.0162, whole-CZ bootstrap SE 0.0090, i.e. 1.8 SEs from zero; 176 CZs).

1. *Coverage.* It averages over 89.0% of the partitioned pupils: 12 CZs drop out, among them New York (F-33 files all of NYC under the Chancellor's Office, which the frame filter removes; the 32 geographic districts are `-2` in F-33), Chicago and the other Illinois CZs, Albuquerque and Alaska (F-33 reports zero in all four teacher-salary fields for those states; 8.3% of all regular districts). The other three targets use all 188 CZs.
2. *Sign and size.* The positive gap comes entirely from the CWIFT deflator: undeflated it is 0.0003 (teacher-weighted mean) and -0.0010 (ratio of totals, same 176 CZs). Across the trimming, deflator-vintage, denominator, sample-rule and CZ-dropping variants (same 176 CZs) it ranges from 0.0039 (hybrid + NYC collapsed + FTE denominator excludes pre-K teacher FTE) to 0.0261 (min group share 20% of CZ pupils (baseline 10%)); it is -0.0039 with unweighted district means and 0.0019 with unweighted CZ averaging. Treat it as a small gap, imprecisely measured (bootstrap SE 0.0090).
3. *Trimming.* Because > 5% of the panel is exactly zero, the lower [0.5%, 99.5%] cut is 0 and only the upper tail is trimmed. The upper cut ($123k deflated) removes genuine high-pay districts (Long Island, Westchester, Greenwich, Bucks County). On the common 176 CZs this changes the gap by +0.0000; but without the trim one more CZ (New York, weight 5.3%) enters with a spurious city group (White Plains and other Westchester cities, no NYC) and a gap of -0.23, moving the aggregate to 0.0024. The screen is thus what keeps that CZ out.
4. *Universe.* F-33 '2018' and CCD 2018 are the same year and the same enrollment (96% exactly equal, 98.6% within 10%), and the four teacher fields are a strict subset of instruction salaries. One mismatch matters: CCD FTE includes pre-K teachers, and the within-CZ regression of log pay per FTE on the pre-K share gives -0.73 (SE 0.28), i.e. pre-K pay is mostly missing from the numerator. Dividing by non-pre-K FTE moves the gap to 0.0043.
5. *Other numerators.* Adding instruction benefits: 0.0090; instruction salaries per teacher: 0.0116 (187 CZs; 0.0077 on the common CZs); regular-program teacher salaries only (which drops special-education, vocational and other teachers from the numerator but not the denominator): 0.0572.
6. *Full-coverage version.* Filling zero-coded wage bills with 0.88 x instruction salaries and collapsing NYC into one district gives 0.0195 (SE 0.0086) on all 188 CZs; also dividing by non-pre-K FTE gives 0.0086 (SE 0.0092). Collapsing NYC alone (zeros missing): 0.0169 on 177 CZs.

**FTE per pupil** (`gap_teachers_pp_locale` = -0.0122, bootstrap SE 0.0139, 188 CZs).

1. The screen has no effect on it: the gap is a ratio of location totals of the unscreened `teachers` and `pupils`. Applying the screen changes it by +0.0001.
2. CCD FTE includes pre-K teachers and CCD enrollment includes pre-K pupils; removing both gives -0.0177 (SE 0.0131). Removing only the teachers (mismatched) would give -0.0014.
3. The sign is carried by the largest CZs. The unweighted CZ average is -0.0481 (SE 0.0095), dropping the five largest CZs gives -0.0389, and restricting to the CZs with a salary gap gives -0.0253. Adding aides to the numerator gives -0.0288.
4. Numerator and denominator are the same CCD district universe (FTE teachers over students, both with pre-K), so the gap is not affected by the F-33 issues above.

**Sample rules matter for the enrollment gap** (section 5): the pupils gap moves from 0.376 to 0.155 when each location must hold >= 20% of CZ pupils, and to -0.015 with unweighted CZ averaging.

## 1. How the salary gap is built

`gap_salary_real_locale` (2018) = pupil-weighted mean across commuting zones of

    log( salary_real[suburb] / salary_real[city] ),

where, within a location (the CZ's suburb or city districts), the statistic is the **teacher-FTE-weighted mean of district `salary_real`** (`stat = 'mean'`, `weight = 'teachers'`), and the CZ weight is the CZ's partitioned pupils. The steps:

1. **Sample frame:** CCD directory 2018 (Urban Institute API `school-districts/ccd/directory/2018`), regular local districts only (agency type 1, or 2 = supervisory-union component). `pupils` = CCD `enrollment` (> 0; F-33 responsible enrollment if missing); `teachers` = CCD `teachers_total_fte` (> 0), i.e. FTE teachers of all levels (pre-K, K, elementary, secondary, ungraded), not aides, coordinators or counselors.
2. **Numerator (F-33 2018, Urban API `school-districts/ccd/finance/2018`):** `wagebill_teachers` = row sum (NaN-skipping, sentinel codes -1/-2/-3/-9/-99 set to missing) of the four teacher-salary fields `salaries_teachers_regular_prog`, `salaries_teachers_sped`, `salaries_teachers_vocational`, `salaries_teachers_other_ed`. These are **teacher salaries only**: no benefits, no aides, substitutes, coordinators or other instruction salaries (those sit in `salaries_instruction`, not used; the fall-back to it exists only when all four fields are missing). `salary_per_teacher = wagebill_teachers / teachers`.
3. **Deflator (NCES EDGE CWIFT):** the district-level `LEA_CWIFTEST` of the release labelled 2018 (school year 2018-19), merged on the nearest release within 2 years (exact for 2018: coverage 96% of districts; on this machine the other LEA vintages are 2013, 2015, 2019, 2021, 2022). `salary_real = salary_per_teacher / cwift`. CWIFT is a labor-cost index built from wages of comparable college-educated non-teachers, so the deflated gap is teacher pay relative to local comparable-worker pay levels, not cost-of-living-adjusted pay.
4. **Screen:** `screen_outliers` blanks `salary_real` (and `salary_per_teacher`, `teachers_pp`, `str_ratio`, ...) outside the panel-wide (all regular districts, one year) [0.5%, 99.5%] quantiles. **It acts on the derived column only**: `gap_teachers_pp` and `gap_salary_per_teacher` are ratios of location totals of the unscreened `teachers`, `pupils`, `wagebill_teachers`, so the screen never touches them (only the `mean`-type `salary_real` gap is screened).
5. **Partition:** districts with a CZ (Dorn 1990 crosswalk via CCD county), >= 500 pupils, NCES locale 11-13 (city) or 21-23 (suburb); a CZ is kept if both groups hold >= 10% of its partitioned pupils (188 CZs). A CZ has a salary gap only if both locations contain at least one district with a non-missing `salary_real`; that gives 176.
6. **Gap and average:** log of the ratio of the suburb to the city teacher-weighted mean; averaged over CZs with weights = CZ partitioned pupils (all partitioned districts, including those without a salary).

### Sample attrition

| step | districts | pupils (millions) |
|---|---|---|
| CCD directory 2018, all agencies | 19,840 | 50.96 |
| regular districts (agency type 1-2) | 13,613 | 48.33 |
| ... with a commuting zone, >= 500 pupils, city or suburb locale, in a CZ that has both | 3,510 | 31.42 |
| ... after dropping CZs where a group has < 10% of pupils (the 188 CZs of the tracked moment) | 3,291 | 28.55 |
| ... with valid CCD FTE teachers (> 0) | 3,285 | 28.52 |
| ... with an F-33 teacher wage bill (non-missing) | 3,249 | 27.53 |
| ... and wage bill > 0 (F-33 zero-coded otherwise) | 2,931 | 25.89 |
| ... and a CWIFT index | 2,885 | 25.81 |
| ... surviving the [0.5%, 99.5%] screen = districts with `salary_real` > 0 | 2,852 | 25.67 |
| `salary_real` non-missing including zeros (what enters the tracked gap) | 3,164 | 27.30 |

12 of 188 CZs have no salary gap, holding **11.0% of the partitioned pupils** (the other targets use all of them).

| group | districts | pupils (m) | FTE valid | F-33 wage bill missing | wage bill == 0 | ... share of pupils | CWIFT missing | `salary_real` > 0 |
|---|---|---|---|---|---|---|---|---|
| city | 656 | 12.03 | 99.7% | 5.8% | 5.0% | 5.5% | 7.9% | 85.7% |
| suburb | 2,635 | 16.52 | 99.8% | 0.2% | 10.8% | 6.0% | 1.3% | 86.9% |

CZs without a salary gap (largest first):

| CZ | largest district | pupils | share of partition | share of CZ pupils in districts with `salary_real` > 0 |
|---|---|---|---|---|
| 19400 | NEW YORK CITY GEOGRAPHIC DISTRICT # 2 | 1,525,410 | 5.3% | 32% |
| 24300 | City of Chicago SD 299 | 1,279,698 | 4.5% | 0% |
| 34901 | ALBUQUERQUE PUBLIC SCHOOLS | 118,768 | 0.4% | 0% |
| 24400 | Rockford SD 205 | 63,395 | 0.2% | 0% |
| 23900 | Peoria SD 150 | 32,981 | 0.1% | 0% |
| 23500 | Champaign CUSD 4 | 26,743 | 0.1% | 0% |
| 24802 | Springfield SD 186 | 22,543 | 0.1% | 0% |
| 34115 | Fairbanks North Star Borough School District | 20,264 | 0.1% | 0% |

## 2. Do numerator and denominator cover the same universe?

**Findings that change what the salary gap measures:**

1. **Zero-coded F-33 teacher salaries.** 1,124 of the regular districts (8.3%) report 0 in all four teacher-salary fields while instruction and total salaries are positive: it is missing data coded as zero, not a sentinel. It is concentrated in whole states: state 02: 53/53 districts (100%); state 35: 89/89 districts (100%); state 17: 853/894 districts (95%); state 23: 50/256 districts (20%); state 50: 23/258 districts (9%); state 33: 15/180 districts (8%) (02 = Alaska, 35 = New Mexico, 17 = Illinois, 23 = Maine, 50 = Vermont, 33 = New Hampshire). Because > 5% of the panel is exactly zero, the lower screen quantile is 0: **the screen has no lower bite**; the upper cuts are $137,686 nominal and $123,301 deflated. Zeros are kept as salaries of 0 in the location averages. They drop a CZ only when a whole group is zero (log of 0). In CZs where some districts are zero and others positive (mostly OH, CT, MA, ME/VT/NH and a few IL/AK ones) the zeros lower that location's mean (see the 'zeros missing' rows in section 3).
2. **New York City has no salary data in the panel.** F-33 reports the whole city under the NYC Chancellor's Office (`3620580`, agency type 3, dropped by the frame filter; F-33 enrollment 960,484, teacher wage bill $7.70 billion), while the CCD frame splits it into 32 geographic districts (960,484 pupils, 64,972 FTE) whose F-33 records are all `-2` (not reported). 32 of 32 have a missing wage bill, so the largest city group in the country has no salary and CZ 19400 (New York) drops from the salary gap while staying in the pupils, FTE-per-pupil and income targets.
3. **Combined effect:** the salary gap averages over 176 CZs holding 89.0% of the partitioned pupils. The excluded CZs include the second and third largest (New York, Chicago) and Albuquerque. The salary target and the other three targets are therefore averages over different populations (see 'restrict to the CZs that have a salary gap' rows for the other targets).

**Other checks.**

- *Teacher salaries are a strict subset of instruction salaries.* Among non-zero districts the ratio (four teacher fields / `salaries_instruction`) has median 0.88, 5th-95th percentile [0.76, 1.00], and exceeds 1.001 in 0 districts. The remainder (about 12%) is instruction pay for aides, substitutes and others, which the CCD teacher FTE does not count: using `salaries_instruction` over teachers alone (section 3) overstates pay per teacher by a factor that varies with aide intensity.
- *Same districts, same year.* CCD 2018 enrollment equals F-33 `enrollment_fall_responsible` exactly in 96.0% of sample districts and within 10% in 98.6% (6 differ by more). F-33 '2018' matches CCD 2018 (exact match with the 2017 and 2019 CCD in only 1.7% (2017), 1.6% (2019)), so there is no one-year offset between numerator and denominator.
- *Pre-K teachers and aides.* CCD FTE includes pre-K teachers (mean 1.8% of teacher FTE), whose pay may sit outside the four F-33 fields. Within-CZ regression of log salary per FTE on the pre-K FTE share and the aide/teacher ratio (non-zero districts, pupil-weighted, CZ-clustered, n = 2520): pre-K share -0.728 (SE 0.283); aide ratio +0.080 (SE 0.028). A coefficient near -1 on the pre-K share would indicate that pre-K pay is missing from the numerator.
- *Timing:* F-33 wage bills are fiscal-year flows; CCD FTE are fall counts; both refer to the school year labelled 2018 in the Urban API.
- *Teacher quality composition:* pay per FTE mixes the wage schedule and teacher experience/education; the target is `Delta log[kappa_l E(h_T^gamma | T, l)]` per the calibration note, so this is intended.

## 3. Sensitivity of the salary gap

Baseline 0.01623 on 176 CZs. Whole-CZ bootstrap SEs, same seed and B as `district_bootstrap.json`; a variant that blanks a district keeps the CZ weights and partition of the baseline unless the variant changes the sample rule.

`gap on baseline CZs` restricts the same variant to the CZs that have a baseline gap (176 for salary), so that a change in coverage is separated from a change in the data; the bootstrap resamples the CZs of each column.

The hybrid rows use 0.881 = median of (teacher salaries / instruction salaries) among non-zero sample districts to fill zero-coded wage bills.

| # | variant | gap | boot SE | n CZs | diff vs baseline | gap on baseline CZs | boot SE | n | diff |
|---|---|---|---|---|---|---|---|---|---|
| 0 | baseline: teacher-weighted mean of screened district `salary_real` (zeros kept; [0.5,99.5] screen) | 0.0162 | 0.0090 | 176 |  | 0.0162 | 0.0090 | 176 |  |
| 1 | pupil-weighted (not teacher-weighted) mean within location | 0.0168 | 0.0095 | 176 | +0.0006 | 0.0168 | 0.0095 | 176 | +0.0006 |
| 2 | unweighted district mean within location | -0.0039 | 0.0074 | 176 | -0.0201 | -0.0039 | 0.0074 | 176 | -0.0201 |
| 3 | ratio of totals, deflated: sum(wage bill/CWIFT) / sum(FTE) | 0.0027 | 0.0158 | 177 | -0.0135 | 0.0166 | 0.0092 | 176 | +0.0004 |
| 4 | undeflated, same aggregation (`salary_per_teacher`, screened, teacher-weighted mean) | 0.0003 | 0.0101 | 176 | -0.0159 | 0.0003 | 0.0100 | 176 | -0.0159 |
| 5 | undeflated ratio of totals (tracked `gap_salary_per_teacher_locale`) | -0.0171 | 0.0180 | 177 | -0.0333 | -0.0010 | 0.0101 | 176 | -0.0172 |
| 6 | baseline with zero-coded wage bills set to missing | 0.0159 | 0.0090 | 176 | -0.0004 | 0.0159 | 0.0090 | 176 | -0.0004 |
| 7 | zeros missing, no trimming | 0.0024 | 0.0158 | 177 | -0.0138 | 0.0162 | 0.0092 | 176 | +0.0000 |
| 8 | zeros missing, trim [0.5%, 99.5%] recomputed on non-zero districts | 0.0159 | 0.0089 | 176 | -0.0003 | 0.0159 | 0.0090 | 176 | -0.0003 |
| 9 | zeros missing, trim [1%, 99%] | 0.0146 | 0.0089 | 176 | -0.0017 | 0.0146 | 0.0089 | 176 | -0.0017 |
| 10 | zeros missing, trim [2.5%, 97.5%] | 0.0142 | 0.0087 | 176 | -0.0021 | 0.0142 | 0.0087 | 176 | -0.0021 |
| 11 | zeros missing, trim [5%, 95%] | 0.0068 | 0.0076 | 171 | -0.0094 | 0.0068 | 0.0076 | 171 | -0.0094 |
| 12 | zeros missing, drop nominal salary/FTE < $25,000 or > $125,000 | 0.0156 | 0.0091 | 176 | -0.0006 | 0.0156 | 0.0091 | 176 | -0.0006 |
| 13 | zeros missing, drop nominal salary/FTE < $25,000 only | 0.0025 | 0.0158 | 177 | -0.0137 | 0.0164 | 0.0092 | 176 | +0.0001 |
| 14 | zeros missing, drop districts whose CCD and F-33 enrollment differ by > 10% | 0.0020 | 0.0158 | 177 | -0.0142 | 0.0159 | 0.0093 | 176 | -0.0004 |
| 15 | zeros missing, drop teacher pay outside 50-100% of instruction salaries | 0.0026 | 0.0158 | 177 | -0.0136 | 0.0165 | 0.0092 | 176 | +0.0002 |
| 16 | zeros missing, drop students-per-teacher outside [6, 30] | 0.0031 | 0.0158 | 177 | -0.0131 | 0.0170 | 0.0092 | 176 | +0.0008 |
| 17 | zeros missing, FTE denominator excludes pre-K teacher FTE | 0.0043 | 0.0091 | 176 | -0.0119 | 0.0043 | 0.0090 | 176 | -0.0119 |
| 18 | zeros missing, CWIFT 2015 release | 0.0155 | 0.0087 | 176 | -0.0007 | 0.0155 | 0.0088 | 176 | -0.0007 |
| 19 | zeros missing, CWIFT 2019 release | 0.0167 | 0.0088 | 176 | +0.0004 | 0.0167 | 0.0089 | 176 | +0.0004 |
| 20 | instruction salaries (F-33 `salaries_instruction`) per FTE teacher, deflated | 0.0116 | 0.0103 | 187 | -0.0046 | 0.0077 | 0.0102 | 176 | -0.0085 |
| 21 | instruction salaries per (FTE teachers + instructional aides), deflated | 0.0181 | 0.0119 | 188 | +0.0019 | 0.0191 | 0.0116 | 176 | +0.0029 |
| 22 | total salaries per total staff FTE, deflated | -0.0323 | 0.0138 | 188 | -0.0486 | -0.0214 | 0.0127 | 176 | -0.0377 |
| 23 | regular-program teacher salaries only per FTE teacher, deflated | 0.0572 | 0.0102 | 176 | +0.0409 | 0.0572 | 0.0106 | 176 | +0.0409 |
| 24 | teacher salaries + instruction benefits per FTE teacher, deflated | 0.0090 | 0.0093 | 176 | -0.0073 | 0.0090 | 0.0093 | 176 | -0.0073 |
| 25 | hybrid: teacher salaries, else 0.88 x instruction salaries where zero-coded (recovers IL, AK, NM) | 0.0188 | 0.0090 | 187 | +0.0026 | 0.0157 | 0.0090 | 176 | -0.0006 |
| 26 | zeros missing + NYC collapsed to one district ([0.5,99.5] on non-zero) | 0.0169 | 0.0085 | 177 | +0.0006 | 0.0159 | 0.0090 | 176 | -0.0003 |
| 27 | hybrid + NYC collapsed | 0.0195 | 0.0086 | 188 | +0.0033 | 0.0157 | 0.0090 | 176 | -0.0006 |
| 28 | hybrid + NYC collapsed + FTE denominator excludes pre-K teacher FTE | 0.0086 | 0.0092 | 188 | -0.0076 | 0.0039 | 0.0090 | 176 | -0.0123 |
| 29 | instruction salaries per FTE teacher + NYC collapsed | 0.0116 | 0.0103 | 187 | -0.0046 | 0.0077 | 0.0102 | 176 | -0.0085 |
| 30 | min district size 0 pupils (baseline 500) | 0.0161 | 0.0089 | 176 | -0.0002 | 0.0161 | 0.0090 | 176 | -0.0002 |
| 31 | min district size 250 pupils (baseline 500) | 0.0161 | 0.0089 | 176 | -0.0001 | 0.0161 | 0.0090 | 176 | -0.0001 |
| 32 | min district size 1000 pupils (baseline 500) | 0.0168 | 0.0089 | 176 | +0.0006 | 0.0168 | 0.0090 | 176 | +0.0006 |
| 33 | min district size 2000 pupils (baseline 500) | 0.0198 | 0.0092 | 169 | +0.0036 | 0.0189 | 0.0088 | 168 | +0.0026 |
| 34 | min group share 5% of CZ pupils (baseline 10%) | 0.0174 | 0.0088 | 189 | +0.0012 | 0.0162 | 0.0090 | 176 | +0.0000 |
| 35 | min group share 20% of CZ pupils (baseline 10%) | 0.0261 | 0.0099 | 148 | +0.0098 | 0.0261 | 0.0097 | 148 | +0.0098 |
| 36 | min group share 30% of CZ pupils (baseline 10%) | 0.0171 | 0.0082 | 100 | +0.0009 | 0.0171 | 0.0080 | 100 | +0.0009 |
| 37 | unweighted average across CZs | 0.0019 | 0.0078 | 176 | -0.0143 | 0.0019 | 0.0074 | 176 | -0.0143 |
| 38 | drop the 1 largest CZs by pupils | 0.0172 | 0.0099 | 175 | +0.0009 | 0.0172 | 0.0097 | 175 | +0.0009 |
| 39 | drop the 5 largest CZs by pupils | 0.0142 | 0.0100 | 173 | -0.0020 | 0.0142 | 0.0103 | 173 | -0.0020 |
| 40 | drop the 10 largest CZs by pupils | 0.0108 | 0.0068 | 168 | -0.0054 | 0.0108 | 0.0066 | 168 | -0.0054 |
| 41 | partition: revenue split (tracked: 0.0154) | 0.0154 | 0.0099 | 593 | -0.0008 | 0.0128 | 0.0119 | 176 | -0.0035 |
| 42 | partition: ses split (tracked: 0.0188) | 0.0187 | 0.0062 | 591 | +0.0025 | 0.0226 | 0.0071 | 175 | +0.0064 |

## 4. Sensitivity of the FTE-per-pupil gap

Baseline -0.01216 on 188 CZs (equals minus the student-teacher-ratio gap).

| # | variant | gap | boot SE | n CZs | diff vs baseline | gap on baseline CZs | boot SE | n | diff |
|---|---|---|---|---|---|---|---|---|---|
| 0 | baseline: ratio of location totals, sum(teachers)/sum(pupils) (screen has no effect) | -0.0122 | 0.0139 | 188 |  | -0.0122 | 0.0139 | 188 |  |
| 1 | apply the [0.5%, 99.5%] screen of `teachers_pp` (drop flagged districts) | -0.0121 | 0.0139 | 187 | +0.0001 | -0.0121 | 0.0139 | 187 | +0.0001 |
| 2 | screen [1%, 99%] | -0.0127 | 0.0139 | 187 | -0.0005 | -0.0127 | 0.0139 | 187 | -0.0005 |
| 3 | screen [2.5%, 97.5%] | -0.0107 | 0.0141 | 186 | +0.0014 | -0.0107 | 0.0141 | 186 | +0.0014 |
| 4 | screen [5%, 95%] | -0.0216 | 0.0165 | 183 | -0.0094 | -0.0216 | 0.0165 | 183 | -0.0094 |
| 5 | drop students-per-teacher outside [6, 30] | -0.0132 | 0.0139 | 187 | -0.0011 | -0.0132 | 0.0139 | 187 | -0.0011 |
| 6 | mean of district ratios, teacher-weighted | -0.0087 | 0.0147 | 188 | +0.0034 | -0.0087 | 0.0147 | 188 | +0.0034 |
| 7 | mean of district ratios, unweighted | -0.0032 | 0.0146 | 188 | +0.0089 | -0.0032 | 0.0146 | 188 | +0.0089 |
| 8 | denominator: F-33 fall enrollment (responsible) instead of CCD enrollment | -0.0108 | 0.0155 | 188 | +0.0013 | -0.0108 | 0.0155 | 188 | +0.0013 |
| 9 | drop districts whose CCD and F-33 enrollment differ by > 10% | -0.0099 | 0.0154 | 188 | +0.0022 | -0.0099 | 0.0154 | 188 | +0.0022 |
| 10 | numerator excludes pre-K teacher FTE, denominator unchanged (pre-K pupils stay in: mismatched) | -0.0014 | 0.0141 | 188 | +0.0108 | -0.0014 | 0.0141 | 188 | +0.0108 |
| 11 | K-12 consistent: (teachers - pre-K teacher FTE) / (pupils - pre-K enrollment) | -0.0177 | 0.0131 | 188 | -0.0056 | -0.0177 | 0.0131 | 188 | -0.0056 |
| 12 | numerator: teachers + instructional aides | -0.0288 | 0.0152 | 184 | -0.0166 | -0.0288 | 0.0152 | 184 | -0.0166 |
| 13 | numerator: total staff FTE | 0.0142 | 0.0428 | 188 | +0.0264 | 0.0142 | 0.0428 | 188 | +0.0264 |
| 14 | restrict to districts in the salary sample (valid, positive `salary_real`) | -0.0261 | 0.0107 | 176 | -0.0139 | -0.0261 | 0.0107 | 176 | -0.0139 |
| 15 | drop IL, AK and NM (zero-coded F-33 states) | -0.0145 | 0.0143 | 177 | -0.0023 | -0.0145 | 0.0143 | 177 | -0.0023 |
| 16 | min district size 0 pupils (baseline 500) | -0.0113 | 0.0142 | 188 | +0.0009 | -0.0113 | 0.0142 | 188 | +0.0009 |
| 17 | min district size 250 pupils (baseline 500) | -0.0117 | 0.0140 | 188 | +0.0005 | -0.0117 | 0.0140 | 188 | +0.0005 |
| 18 | min district size 1000 pupils (baseline 500) | -0.0137 | 0.0139 | 186 | -0.0016 | -0.0137 | 0.0139 | 186 | -0.0016 |
| 19 | min district size 2000 pupils (baseline 500) | -0.0178 | 0.0131 | 180 | -0.0056 | -0.0171 | 0.0130 | 179 | -0.0049 |
| 20 | min group share 5% of CZ pupils (baseline 10%) | -0.0163 | 0.0138 | 202 | -0.0041 | -0.0122 | 0.0139 | 188 | +0.0000 |
| 21 | min group share 20% of CZ pupils (baseline 10%) | -0.0095 | 0.0156 | 158 | +0.0027 | -0.0095 | 0.0156 | 158 | +0.0027 |
| 22 | min group share 30% of CZ pupils (baseline 10%) | -0.0029 | 0.0171 | 108 | +0.0092 | -0.0029 | 0.0171 | 108 | +0.0092 |
| 23 | unweighted average across CZs | -0.0481 | 0.0095 | 188 | -0.0359 | -0.0481 | 0.0095 | 188 | -0.0359 |
| 24 | drop the 1 largest CZs by pupils | -0.0125 | 0.0150 | 187 | -0.0003 | -0.0125 | 0.0150 | 187 | -0.0003 |
| 25 | drop the 5 largest CZs by pupils | -0.0389 | 0.0098 | 183 | -0.0267 | -0.0389 | 0.0098 | 183 | -0.0267 |
| 26 | drop the 10 largest CZs by pupils | -0.0421 | 0.0086 | 178 | -0.0299 | -0.0421 | 0.0086 | 178 | -0.0299 |
| 27 | restrict to the CZs that have a salary gap | -0.0253 | 0.0104 | 176 | -0.0132 | -0.0253 | 0.0104 | 176 | -0.0132 |
| 28 | partition: revenue split (tracked: 0.0366) | 0.0366 | 0.0060 | 621 | +0.0488 | 0.0378 | 0.0075 | 187 | +0.0499 |
| 29 | partition: ses split (tracked: -0.0036) | -0.0036 | 0.0068 | 618 | +0.0086 | -0.0040 | 0.0085 | 187 | +0.0082 |

## 5. Other targets under the same sample rules


**`pupils`**

| # | variant | gap | boot SE | n CZs | diff vs baseline | gap on baseline CZs | boot SE | n | diff |
|---|---|---|---|---|---|---|---|---|---|
| 0 | baseline | 0.3764 | 0.1154 | 188 |  | 0.3764 | 0.1154 | 188 |  |
| 1 | min district size 0 pupils (baseline 500) | 0.3795 | 0.1156 | 188 | +0.0032 | 0.3795 | 0.1156 | 188 | +0.0032 |
| 2 | min district size 250 pupils (baseline 500) | 0.3790 | 0.1156 | 188 | +0.0026 | 0.3790 | 0.1156 | 188 | +0.0026 |
| 3 | min district size 1000 pupils (baseline 500) | 0.3673 | 0.1156 | 186 | -0.0090 | 0.3673 | 0.1156 | 186 | -0.0090 |
| 4 | min district size 2000 pupils (baseline 500) | 0.3355 | 0.1167 | 180 | -0.0408 | 0.3180 | 0.1143 | 179 | -0.0584 |
| 5 | min group share 5% (baseline 10%) | 0.4278 | 0.1300 | 202 | +0.0515 | 0.3764 | 0.1154 | 188 | +0.0000 |
| 6 | min group share 20% (baseline 10%) | 0.1554 | 0.0900 | 158 | -0.2210 | 0.1554 | 0.0900 | 158 | -0.2210 |
| 7 | min group share 30% (baseline 10%) | 0.1139 | 0.0875 | 108 | -0.2624 | 0.1139 | 0.0875 | 108 | -0.2624 |
| 8 | unweighted average across CZs | -0.0152 | 0.0714 | 188 | -0.3915 | -0.0152 | 0.0714 | 188 | -0.3915 |
| 9 | drop the 1 largest CZs by pupils | 0.4045 | 0.1266 | 187 | +0.0282 | 0.4045 | 0.1266 | 187 | +0.0282 |
| 10 | drop the 5 largest CZs by pupils | 0.3755 | 0.1206 | 183 | -0.0009 | 0.3755 | 0.1206 | 183 | -0.0009 |
| 11 | drop the 10 largest CZs by pupils | 0.2209 | 0.1029 | 178 | -0.1554 | 0.2209 | 0.1029 | 178 | -0.1554 |
| 12 | restrict to the CZs that have a salary gap | 0.4305 | 0.1166 | 176 | +0.0541 | 0.4305 | 0.1166 | 176 | +0.0541 |
| 13 | partition: revenue split | -0.0247 | 0.0218 | 621 | -0.4011 | -0.0030 | 0.0214 | 187 | -0.3794 |
| 14 | partition: ses split | -0.0013 | 0.0235 | 619 | -0.3776 | 0.0106 | 0.0200 | 187 | -0.3657 |

**`seda_lninc`**

| # | variant | gap | boot SE | n CZs | diff vs baseline | gap on baseline CZs | boot SE | n | diff |
|---|---|---|---|---|---|---|---|---|---|
| 0 | baseline | 0.2708 | 0.0268 | 187 |  | 0.2708 | 0.0271 | 187 |  |
| 1 | min district size 0 pupils (baseline 500) | 0.2710 | 0.0268 | 187 | +0.0002 | 0.2710 | 0.0272 | 187 | +0.0002 |
| 2 | min district size 250 pupils (baseline 500) | 0.2709 | 0.0268 | 187 | +0.0001 | 0.2709 | 0.0272 | 187 | +0.0001 |
| 3 | min district size 1000 pupils (baseline 500) | 0.2698 | 0.0270 | 185 | -0.0010 | 0.2698 | 0.0267 | 185 | -0.0010 |
| 4 | min district size 2000 pupils (baseline 500) | 0.2701 | 0.0271 | 178 | -0.0007 | 0.2680 | 0.0277 | 177 | -0.0028 |
| 5 | min group share 5% (baseline 10%) | 0.2672 | 0.0270 | 200 | -0.0036 | 0.2708 | 0.0271 | 187 | +0.0000 |
| 6 | min group share 20% (baseline 10%) | 0.2494 | 0.0286 | 157 | -0.0214 | 0.2494 | 0.0279 | 157 | -0.0214 |
| 7 | min group share 30% (baseline 10%) | 0.2250 | 0.0291 | 107 | -0.0458 | 0.2250 | 0.0292 | 107 | -0.0458 |
| 8 | unweighted average across CZs | 0.2658 | 0.0149 | 187 | -0.0050 | 0.2658 | 0.0152 | 187 | -0.0050 |
| 9 | drop the 1 largest CZs by pupils | 0.2896 | 0.0238 | 186 | +0.0188 | 0.2896 | 0.0239 | 186 | +0.0188 |
| 10 | drop the 5 largest CZs by pupils | 0.2830 | 0.0231 | 182 | +0.0122 | 0.2830 | 0.0243 | 182 | +0.0122 |
| 11 | drop the 10 largest CZs by pupils | 0.2754 | 0.0225 | 177 | +0.0046 | 0.2754 | 0.0226 | 177 | +0.0046 |
| 12 | restrict to the CZs that have a salary gap | 0.2751 | 0.0296 | 175 | +0.0043 | 0.2751 | 0.0297 | 175 | +0.0043 |
| 13 | partition: revenue split | 0.1573 | 0.0184 | 615 | -0.1135 | 0.1850 | 0.0237 | 186 | -0.0858 |
| 14 | partition: ses split | 0.3542 | 0.0152 | 619 | +0.0834 | 0.4041 | 0.0157 | 187 | +0.1333 |

## 6. Influence of the largest CZs

| CZ | largest district in the CZ | weight | gap pupils | gap salary | gap FTE/pupil | gap income |
|---|---|---|---|---|---|---|
| 38300 | Los Angeles Unified | 9.1% | 0.094 | 0.008 | -0.009 | 0.083 |
| 19400 | NEW YORK CITY GEOGRAPHIC DISTRICT # 2 | 5.3% | -0.551 |  | 0.169 | 0.160 |
| 24300 | City of Chicago SD 299 | 4.5% | 0.665 |  | 0.067 | 0.327 |
| 32000 | HOUSTON ISD | 3.9% | 0.662 | 0.084 | 0.043 | 0.281 |
| 19600 | Newark Public School District | 3.0% | 2.093 | 0.003 | 0.125 | 0.634 |
