# Education and occupation: evidence for a joint BA and occupational choice

Generated 2026-10-05T02:24:43 by `education_occupation.py` (DESIGN.md section M3, data task D3). Sources are listed below; SEs are Census successive-difference replicate weights (80 replicates) for the ACS, the household bootstrap (B = 500, seed 20260929, the draws of `nlsy_teachers.py`) for the NLSY79, and the delta method for NCES. Machine-readable results: `education_occupation.json`; moments for the Julia code: `education_targets.toml`.

**Sources.**

- ACS 2009-13 5-year PUMS (data years 2009-13), Census API `https://api.census.gov/data/2013/acs/acs5/pums`, cached in `raw/acs_pums/`; sample, HHJK rules and K-12 codes of `acs_occupations.py`; the 80 replicate weights PWGTP1-PWGTP80 (`--fetch-repw`).
- NCES, *Digest of Education Statistics* 2022, Table 211.20 (average base salary of full-time public school teachers by highest degree and years of experience), <https://nces.ed.gov/programs/digest/d22/tables/dt22_211.20.asp> (2020 edition, `d20/tables/dt20_211.20.asp`, for the check of the 2017-18 rows).
- NCES, *Digest* 2018, Table 330.10 (average undergraduate tuition and fees and room and board), <https://nces.ed.gov/programs/digest/d18/tables/dt18_330.10.asp>; Table 303.70 (undergraduate enrollment), <https://nces.ed.gov/programs/digest/d18/tables/dt18_303.70.asp>; *Digest* 2019, Table 331.30 (net price), <https://nces.ed.gov/programs/digest/d19/tables/dt19_331.30.asp>.
- College Board, *Trends in College Pricing 2018* (Ma, Pender and Welch), Figure 9 (public four-year, in-state; printed p. 18) and Figure 10 (private nonprofit four-year; p. 19), <https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf>.
- BLS, Consumer Expenditure Survey 2018, Table R-1 (mean annual expenditures per consumer unit, people, children and earners per unit, college tuition), <https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf>.
- NLSY79 through the `nlsy_teachers.py` / `nlsy_ability.py` machinery (latent AFQT, 1979 weights); the published teacher gaps are in `nlsy_teachers.json`.

## Summary

- **BA shares (ACS, ages 25-34, Table 2 universe, home production included).** BA or more (SCHL >= 21): 0.339 (0.0011) of everyone (men 0.301 (0.0011), women 0.376 (0.0013)); 0.323 (0.0007) of non-teachers (T2 convention; men 0.294 (0.0007), women 0.352 (0.0010)); 0.874 (0.0017) of the K-12 teachers of the project's definition but 0.961 (0.0011) of core K-12 teachers (elementary, secondary, special education; pre-K/kindergarten 0.472 (0.0084), other instructors 0.594 (0.0085)). Graduate degree: 0.378 (0.0021) of teachers.
- **Schooling by BA status.** Mean years (capped 9-17): BA+ non-teachers 16.29 (0.001), non-BA non-teachers 12.43 (0.002); all non-teachers 13.68 (0.004) (workbook 13.64), teachers 16.00 (0.006) (BA+ teachers 16.43 (0.002)).
- **BA wage premium** among non-teacher market workers (log hourly wage, BA+ minus non-BA): 0.481 (0.0012) raw (men 0.479 (0.0019), women 0.508 (0.0016)), 0.490 (0.0013) controlling for age and sex, 0.361 (0.0013) also within occupation. The model has no BA-specific occupational sorting, so the raw (total return) and the within-occupation premium are alternative targets with different meanings; the choice is the M3 designer's (section 1.2).
- **Teachers' pay (section 1.3).** All K-12 teachers (codes 2300-2340) earn +0.0557 (0.0049) (men) and +0.0860 (0.0031) (women) against all non-teacher market workers, the existing targets +0.0557 / +0.0860 reproduced, and -0.261 (0.0049) and -0.208 (0.0033) against non-teacher BA+ workers. The 12.6% of them without a BA (pre-K/kindergarten, other instructors) earn 0.48 (0.009) less than those with one, so for a model whose teachers all have a BA the like-for-like rows are BA+ teachers: +0.079 (0.0052) and +0.157 (0.0029) against all non-teachers, -0.237 (0.0053) and -0.137 (0.0032) against BA+ non-teachers (core codes -0.237 (0.0054) and -0.125 (0.0030)). Teachers with a graduate degree earn +0.195 (0.0043) more than BA-only teachers (age-controlled); the NTPS 2017-18 salary schedule gives +0.156 (0.0068) (master's over bachelor's, 1-5 years of experience).
- **Weeks rule (section 1.4).** The Table 2 wage sample (48+ weeks worked) leaves out 20.6% of full-time K-12 teachers against 12.1% of other full-time workers, mostly teachers reporting 40-47 weeks (12.7%), whose hourly wage is higher for the same salary. Admitting those records raises the teacher premium against all non-teachers from +0.056 / +0.086 (men / women) to +0.074 / +0.112, and the BA+-teacher premium against BA+ non-teachers from -0.237 / -0.137 to -0.217 / -0.113; the existing teacher-premium targets, which use the rule, are probably 0.02-0.04 lower than a sample that keeps those teachers would give.
- **Ability (NLSY79).** Latent AFQT of persons with 16+ years of schooling minus the rest: +1.227 (0.0197) SD (employed +1.208 (0.0201); men +1.30 (0.031), women +1.16 (0.029)). P(BA | observed AFQT tercile) = 0.041 (0.0040), 0.167 (0.0094), 0.567 (0.0161) at the NLSY79 BA share (0.258); 0.057, 0.277, 0.683 at the ACS share (0.339); the latent logit slope is 2.00 (0.068) per latent SD (observed scale 1.81 (0.057)). The existing teacher gaps are confirmed to six digits: +0.614 (0.0616) against employed non-teachers and -0.334 (0.0625) against non-teacher college graduates.
- **Cost of a BA (2017-18).** Published tuition and fees $9,037 a year at public 4-year institutions (in-state) and $15,924 at all 4-year institutions; net of grants and tax benefits $3,780 (public) to $14,500 (private nonprofit). A one-off cost paid by one person is measured against consumption *per adult* ($32,223: the CE's $61,224 per consumer unit over 1.9 adults), over the working life. Central F_c = 0.0208 of C-bar (net tuition and fees, undiscounted); 0.0360 if annuitized at 3%; an aggregate-flow check gives 0.0384. The first release's per-consumer-unit figures (0.011 undiscounted, 0.019 at 3%) understate F_c by the 1.9 adults per unit (section 4).

*Changes from the first release (2026-10-05 01:14).* (1) tuition is measured against expenditures per adult, not per consumer unit (a factor 1.9), and the undiscounted figure is central (F_c 0.0189 -> 0.0208); (2) the teacher comparisons are extended to BA+ teachers and the core codes, which are the like-for-like counterparts of a model whose teachers all have a BA; (3) the wage-sample weeks rule and its effect on the teacher premia are quantified; (4) the BA-premium choice (raw against within occupation) is left to the M3 designer; (5) the ability gradient is reported on the observed and the latent scale and at the ACS BA share; (6) the TOML carries `model_moment` names, `choice_group`s, scales for every target, correct universe sizes and URLs for every source. No ACS, NLSY or NCES statistic of the first release changed; only the tuition ratios did.

## 1. ACS PUMS 2009-13: education, occupation and wages at 25-34

**Sample.** `acs_occupations.build_sample(rules="workbook")` on the 2013 5-year PUMS (data years 2009-13): ages 25-34, SCHL >= 12 (at least grade 9), armed forces and unemployed persons with 15+ usual hours dropped, 50 states and DC, group quarters in: 1,632,604 persons (795,434 men, 837,170 women), 980,062 in the wage sample. Home production is out of the labor force or under 15 hours; 15-29 hours are split 0.5/0.5 between home production and the occupation, so every person has total weight 1 (Table 2 counts are sums of these unweighted weights). **K-12 teachers** are Census codes 2300-2340. **BA+** is SCHL >= 21 (bachelor's, master's, professional or doctoral degree); **graduate degree** SCHL >= 22. **Years of schooling** are the workbook's mapping (grade 9 = 9, high school or GED 12, some college 13, associate 14, bachelor 16, graduate 19) capped at 9-17, so a graduate degree counts 17. **Wages:** the wage sample of Table 2 (employed, 30+ hours, 48+ weeks, income >= $1,000 in 2010 dollars), log hourly wage in 1999 dollars (annual wage and business income over usual hours times weeks), PWGTP-weighted.

**Two aggregation conventions.** *P*: every row weighted by PWGTP (the plain population mean). *T2*: the convention of the existing premium targets (`teacher_premium_m/f`): occupation x gender cell means (PWGTP-weighted) combined with the *unweighted* Table 2 counts of the population in question (for a BA subset, the counts of BA holders in each cell), because the model's occupation shares match those counts. They differ by up to 0.01 (K-12 teachers are 1.94% of men by counts and 1.76% PWGTP-weighted). The pooled teacher premia average the two gender premia with the K-12 wage-sample weights 0.2444 / 0.7556 (the `teacher_p9010` convention; used for every teacher group: the men's share of the BA+ and core-code teacher records, 0.247 and 0.252, would change the pooled figures by at most 0.0009), which removes the gender-mix effect (the unadjusted pooled gap is +0.041 (0.0027) against all non-teachers and -0.262 (0.0028) against BA+). SEs hold the Table 2 counts fixed.

### 1.1 BA shares and years of schooling

|  | Men | Women | Pooled |
|---|---:|---:|---:|
| BA+ share, all (P) | 0.301 (0.0011) | 0.376 (0.0013) | 0.339 (0.0011) |
| BA+ share, non-teachers (T2) | 0.294 (0.0007) | 0.352 (0.0010) | 0.323 (0.0007) |
| BA+ share, non-teachers (P) | 0.291 (0.0011) | 0.347 (0.0012) | 0.319 (0.0011) |
| BA+ share, non-teacher market workers, home production out (T2) | 0.333 (0.0007) | 0.412 (0.0011) | 0.369 (0.0008) |
| BA+ share, K-12 teachers (P) | 0.881 (0.0031) | 0.872 (0.0020) | 0.874 (0.0017) |
| Graduate degree, share of K-12 teachers | 0.354 (0.0037) | 0.385 (0.0026) | 0.378 (0.0021) |
| Graduate degree, share of BA+ K-12 teachers | 0.402 (0.0042) | 0.442 (0.0028) | 0.432 (0.0023) |
| Teachers' share of all 25-34 (T2 counts) | 0.0194 | 0.0599 | 0.0402 |
| Mean years, non-teachers (T2) | 13.50 (0.003) | 13.85 (0.005) | 13.68 (0.004) |
| Mean years, non-teachers (P) | 13.49 (0.005) | 13.83 (0.006) | 13.66 (0.005) |
| Mean years, K-12 teachers (P) | 15.99 (0.011) | 16.00 (0.008) | 16.00 (0.006) |
| Mean years, BA+ non-teachers (T2) | 16.27 (0.001) | 16.30 (0.001) | 16.29 (0.001) |
| Mean years, non-BA non-teachers (T2) | 12.34 (0.002) | 12.52 (0.002) | 12.43 (0.002) |
| Mean years, BA+ teachers (P) | 16.40 (0.004) | 16.44 (0.003) | 16.43 (0.002) |
| Mean years, non-BA teachers (P) | 12.92 (0.024) | 13.01 (0.011) | 12.99 (0.010) |
| Teachers' years over non-teachers' (T2) | 1.184 (0.0008) | 1.155 (0.0005) | 1.170 (0.0004) |

Dividing years by 25 gives the model's time shares: non-BA non-teachers 0.497, BA+ non-teachers 0.652, BA+ teachers 0.657 (a BA is 16/25 = 0.64; graduate degrees, 17 years after the cap, are 29% of BA+ non-teachers and 43% of BA+ teachers).

### 1.2 The BA wage premium, and which premium to match

|  | Men | Women | Pooled |
|---|---:|---:|---:|
| BA+ minus non-BA, non-teacher market workers: raw (T2) | 0.479 (0.0019) | 0.508 (0.0016) | 0.481 (0.0012) |
|   raw (P) | 0.486 (0.0020) | 0.507 (0.0018) | 0.483 (0.0013) |
|   controlling for single-year age (pooled: and sex) (P) | 0.475 (0.0020) | 0.505 (0.0017) | 0.490 (0.0013) |
|   also occupation x gender fixed effects (P) | 0.326 (0.0021) | 0.396 (0.0020) | 0.361 (0.0013) |

Holding occupation x gender fixed lowers the age-controlled premium from 0.490 to 0.361: about 26% of it is BA holders' sorting into better-paid occupations. **Which premium to match.** The model has one credential premium rho_c in every non-teaching occupation, and its non-teaching occupation choice does not depend on d, so the model's raw d = 1 minus d = 0 premium among non-teaching market workers equals its within-occupation premium (plus selection on z) by construction. Matching the data's *raw* premium (0.481) assigns the between-occupation sorting (about 0.12 log points) to rho_c: it is the total return to the college decision, which is what college-going responds to, but it overstates the wage structure within occupations. Matching the *within-occupation* premium (0.361) fits that wage structure but leaves the model's total return to college below the data's, which could understate the incentive to attend. Both are in the TOML with the same `choice_group`; a criterion should contain only one, and which one is primary is the M3 designer's choice.

### 1.3 Teachers' pay: which teachers, which comparison

The K-12 group of the existing targets (Census codes 2300-2340) includes pre-K/kindergarten teachers and other teachers and instructors, 12.6% of the group in all lack a BA (section 1.5). A model in which teaching requires a BA has only d = 1 teachers, so the premium of the model's teachers has several data counterparts (T2 convention, log hourly wage; the first two rows are the existing figures):

| Teachers | Compared with | Men (T2) | Women (T2) | Pooled (T2) | Men (P) | Women (P) | Pooled (P) | Pooled, age-controlled (P) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| all K-12 teachers (2300-2340) | all non-teacher workers | +0.0557 (0.0049) | +0.0860 (0.0031) | +0.0786 (0.0026) | +0.051 (0.0050) | +0.081 (0.0032) | +0.074 (0.0027) | +0.072 (0.0027) |
| all K-12 teachers (2300-2340) | BA+ non-teacher workers | -0.2610 (0.0049) | -0.2081 (0.0033) | -0.2211 (0.0028) | -0.270 (0.0050) | -0.214 (0.0034) | -0.227 (0.0028) | -0.226 (0.0029) |
| BA+ K-12 teachers | all non-teacher workers | +0.0793 (0.0052) | +0.1574 (0.0029) | +0.1383 (0.0024) | +0.075 (0.0053) | +0.152 (0.0031) | +0.133 (0.0026) | +0.131 (0.0026) |
| BA+ K-12 teachers | BA+ non-teacher workers | -0.2374 (0.0053) | -0.1368 (0.0032) | -0.1614 (0.0026) | -0.246 (0.0054) | -0.142 (0.0032) | -0.168 (0.0026) | -0.168 (0.0027) |
| core-code teachers (2310-2330) | all non-teacher workers | +0.0794 (0.0054) | +0.1691 (0.0028) | +0.1472 (0.0023) | +0.075 (0.0055) | +0.164 (0.0030) | +0.142 (0.0024) | +0.138 (0.0024) |
| core-code teachers (2310-2330) | BA+ non-teacher workers | -0.2373 (0.0054) | -0.1250 (0.0030) | -0.1524 (0.0024) | -0.246 (0.0055) | -0.131 (0.0031) | -0.159 (0.0024) | -0.161 (0.0025) |

Wage-sample record-count weights for the first row give +0.0420 (0.0049) / +0.0704 (0.0032) (the calibration log's check, +0.042 / +0.070). Which row is the counterpart of 'teachers versus d = 1 non-teachers' depends on the model: (a) if every teacher is forced to d = 1, the model's teachers are the BA+ teachers, and the counterpart is **BA+ teachers against BA+ non-teachers** (fourth row), the wage analogue of NLSY comparison C (teachers with 16+ years against non-teachers with 16+ years, section 3), with BA+ teachers against all non-teachers (third row) as the counterpart of the existing teacher_premium targets; core-code teachers (last two rows) are a robustness variant that drops the non-BA codes; (b) if the model also has a non-BA teaching track, all K-12 teachers against d = 1 non-teachers (second row, the figure of the first release) is the counterpart of the all-teacher comparison and the fourth row that of the d = 1 teachers. Moving from the second to the fourth row raises teachers' relative pay by 0.060 (pooled), 27% of the shortfall against BA+ workers (women 0.071, men 0.024): that part is the low pay of the non-BA teachers, and it matters most for women, who are most of the non-BA teachers.

|  | Men | Women | Pooled |
|---|---:|---:|---:|
| BA+ minus non-BA K-12 teachers: raw (P) | +0.202 (0.0219) | +0.566 (0.0095) | +0.482 (0.0094) |
|   controlling for single-year age (pooled: and sex) | +0.189 (0.0213) | +0.561 (0.0097) | +0.474 (0.0092) |
|   also Census code (2300-2340; pooled: code x gender) | +0.151 (0.0240) | +0.418 (0.0125) | +0.354 (0.0114) |
| Teachers with a graduate degree minus BA-only teachers: raw | +0.196 (0.0103) | +0.235 (0.0042) | +0.225 (0.0042) |
|   controlling for single-year age (pooled: and sex) | +0.162 (0.0105) | +0.206 (0.0044) | +0.195 (0.0043) |

Records in the wage sample: non-teachers 930,479 (369,451 BA+), K-12 teachers 49,583 (44,370 BA+, 5,213 without a BA, 40,761 in the core codes; 20,040 with a graduate degree, 24,330 with exactly a BA). Mean log wages (T2): teachers 2.599 (0.0048) (men) and 2.525 (0.0031) (women); BA+ teachers 2.623 (0.0052) and 2.596 (0.0028); all non-teachers 2.543 (0.0011) and 2.439 (0.0010); BA+ non-teachers 2.860 (0.0014) and 2.733 (0.0012). The BA-minus-non-BA teacher gap is the data counterpart of the pay penalty of a non-BA teaching track (rho_T0 in the model's terms); the graduate-degree rows are the ACS counterpart of the salary-schedule premium of section 2 (different objects: hourly wages at ages 25-34, including professional and doctoral degrees and district composition, against base salary at equal experience).

### 1.4 Weeks worked and the wage sample

The Table 2 wage sample keeps full-time workers (30+ usual hours) who report 48 or more weeks worked in the past year (bracket 48-49 or 50-52 weeks) and at least $1,000 of income. The weeks bracket is selective for teachers (plausibly because many are on school-year contracts): of the 63,220 full-time K-12 teachers in the sample, 20.6% are outside the wage sample, against 12.1% of the 1,057,768 other full-time workers (PWGTP shares of full-time workers):

| Full-time workers | K-12 teachers | BA+ K-12 teachers | Non-teachers | BA+ non-teachers |
|---|---:|---:|---:|---:|
| in the wage sample | 0.794 (0.0018) | 0.789 (0.0020) | 0.879 (0.0004) | 0.902 (0.0006) |
| 48+ weeks but outside (income floor) | 0.000 (0.0001) | 0.000 (0.0001) | 0.001 (<0.0001) | 0.001 (<0.0001) |
| 40-47 weeks (WKW 3) | 0.127 (0.0014) | 0.134 (0.0016) | 0.042 (0.0003) | 0.039 (0.0004) |
| 27-39 weeks (WKW 4) | 0.055 (0.0011) | 0.055 (0.0011) | 0.037 (0.0002) | 0.029 (0.0003) |
| 14-26 weeks (WKW 5) | 0.014 (0.0006) | 0.013 (0.0006) | 0.022 (0.0002) | 0.016 (0.0002) |
| 13 or fewer weeks (WKW 6) | 0.010 (0.0005) | 0.009 (0.0005) | 0.019 (0.0002) | 0.013 (0.0002) |
| weeks missing | 0.000 | 0.000 | 0.000 | 0.000 |

Teachers who report fewer weeks (40-47, a school-year contract, often with the summer unpaid or not counted) have a higher implied hourly wage for the same annual salary, so dropping them removes relatively well-paid teachers from the hourly sample and the teacher premia of the Table 2 rule are probably biased down. Sensitivity: the wage sample also admits the 40-47 week bracket (weeks set to the midpoint 43.5), for teachers and non-teachers alike, or also the 27-39 bracket (weeks 33); SEs are recomputed with the replicate weights.

| Premium (log points, T2) | 48+ weeks (Table 2 rule) | also 40-47 weeks | also 27-47 weeks |
|---|---:|---:|---:|
| All K-12 teachers minus all non-teachers (existing targets): men | +0.0557 (0.0049) | +0.0743 (0.0043) | +0.0865 (0.0041) |
| All K-12 teachers minus all non-teachers (existing targets): women | +0.0860 (0.0031) | +0.1122 (0.0030) | +0.1283 (0.0029) |
| All K-12 teachers minus all non-teachers (existing targets): pooled | +0.0786 (0.0026) | +0.1031 (0.0024) | +0.1184 (0.0023) |
| All K-12 teachers minus BA+ non-teachers: men | -0.2610 (0.0049) | -0.2422 (0.0045) | -0.2313 (0.0042) |
| All K-12 teachers minus BA+ non-teachers: women | -0.2081 (0.0033) | -0.1816 (0.0032) | -0.1677 (0.0031) |
| All K-12 teachers minus BA+ non-teachers: pooled | -0.2211 (0.0028) | -0.1962 (0.0026) | -0.1829 (0.0025) |
| BA+ teachers minus all non-teachers: men | +0.0793 (0.0052) | +0.0994 (0.0048) | +0.1123 (0.0046) |
| BA+ teachers minus all non-teachers: women | +0.1574 (0.0029) | +0.1810 (0.0028) | +0.1977 (0.0026) |
| BA+ teachers minus all non-teachers: pooled | +0.1383 (0.0024) | +0.1614 (0.0023) | +0.1774 (0.0021) |
| BA+ teachers minus BA+ non-teachers: men | -0.2374 (0.0053) | -0.2172 (0.0050) | -0.2055 (0.0047) |
| BA+ teachers minus BA+ non-teachers: women | -0.1368 (0.0032) | -0.1129 (0.0030) | -0.0984 (0.0029) |
| BA+ teachers minus BA+ non-teachers: pooled | -0.1614 (0.0026) | -0.1379 (0.0025) | -0.1238 (0.0024) |
| BA+ minus non-BA non-teachers, raw: men | 0.4795 (0.0019) | 0.4791 (0.0019) | 0.4807 (0.0018) |
| BA+ minus non-BA non-teachers, raw: women | 0.5081 (0.0016) | 0.5081 (0.0017) | 0.5111 (0.0017) |
| BA+ minus non-BA non-teachers, raw: pooled | 0.4806 (0.0012) | 0.4806 (0.0012) | 0.4830 (0.0012) |

The wage sample grows from 980,062 to 1,035,863 and 1,077,693 records. Admitting the 40-47 week records raises every teacher premium, by 0.019 to 0.027 (0.030 to 0.042 with the 27-39 bracket as well): well above the sampling SEs (0.003-0.005) and within the 0.05 modelling scale; the BA premium is unchanged (raw) or moves by less than 0.003. The sign is the one argued above (the 48+ week rule drops teachers whose implied hourly wage is high), and the teacher premium targets, which use the rule, are probably 0.02-0.04 lower than a sample that keeps those teachers would give (the hours of a 40-47 week worker and the weeks within the bracket are measured with error, so the size is indicative). The rule is inherited from the workbook (HHJK) and the existing targets use it, so the main tables keep it; both variants are in the JSON and as `alt_weeks_40plus` / `alt_weeks_27plus` in the TOML.

### 1.5 K-12 teachers by Census occupation code

| Code | Share of K-12 teachers (pooled) | BA+ men | BA+ women | BA+ pooled |
|---|---:|---:|---:|---:|
| 2300 pre-K and kindergarten teachers | 0.106 (0.0015) | 0.572 (0.0443) | 0.469 (0.0083) | 0.472 (0.0084) |
| 2310 elementary and middle school teachers | 0.612 (0.0023) | 0.965 (0.0025) | 0.963 (0.0014) | 0.963 (0.0013) |
| 2320 secondary school teachers | 0.145 (0.0017) | 0.955 (0.0044) | 0.976 (0.0027) | 0.967 (0.0024) |
| 2330 special education teachers | 0.043 (0.0011) | 0.885 (0.0194) | 0.912 (0.0083) | 0.907 (0.0068) |
| 2340 other teachers and instructors | 0.095 (0.0015) | 0.517 (0.0134) | 0.648 (0.0110) | 0.594 (0.0085) |
| core K-12 (2310 + 2320 + 2330) | 0.800 (0.0022) | 0.958 (0.0022) | 0.961 (0.0013) | 0.961 (0.0011) |

Licensure binds for the core K-12 codes (96% BA+). The project's K-12 group also contains pre-K/kindergarten teachers and 'other teachers and instructors' (20% of teachers by PWGTP, 47% and 59% BA+), so only 87% of the group has a BA. A model that imposes s >= 16/25 on every teacher is matched to teaching shares that include these codes.

### 1.6 Reproduction of the published block

The sample reproduces `acs_occupations.json` exactly: maximum absolute differences are 0.0 persons in the 42 Table 2 counts, 0.0 wage-sample records, 6.0e-14 in cell mean log wages and 0.0e+00 years in cell mean schooling. The T2 teacher premia are 0.0557 / 0.0860 (coded targets 0.0557 / 0.0860); with wage-sample weights 0.0420 / 0.0704. Years of schooling: K-12 teachers 15.988 (men) and 16.003 (women), non-teachers 13.662 pooled, against 16.04 and 13.64 in the Mincer workbook (the calibration's constants 16.0437 and 13.6419 come from the workbook, which cannot be matched exactly).

## 2. Teacher salary schedules (NCES Digest Table 211.20)

Source: Digest of Education Statistics 2022, Table 211.20 (<https://nces.ed.gov/programs/digest/d22/tables/dt22_211.20.asp>), average base salary (current dollars) of full-time public school teachers by highest degree and years of teaching experience; the 2017-18 rows are the National Teacher and Principal Survey (NTPS) and are identical in the 2020 and 2022 editions (max absolute difference 0 over 114 cells). Base salary excludes supplemental contracts and nonschool income; experience counts all school years, full or part time. Standard errors are the table's; the SE of the log ratio treats the two means as independent (the two groups are disjoint teachers).

| 2017-18, experience | Bachelor's | Master's | Master's - bachelor's ($) | log(M/B) (SE) |
|---|---:|---:|---:|---:|
| 1 year or less | 41,520 (363) | 47,150 (785) | 5,630 | 0.127 (0.019) |
| 2 years | 43,000 (372) | 51,720 (702) | 8,720 | 0.185 (0.016) |
| 3 years | 43,510 (415) | 50,840 (628) | 7,330 | 0.156 (0.016) |
| 4 years | 43,620 (309) | 50,900 (530) | 7,280 | 0.154 (0.013) |
| 5 years | 44,350 (356) | 51,330 (526) | 6,980 | 0.146 (0.013) |
| 6 to 9 years | 45,990 (272) | 54,200 (339) | 8,210 | 0.164 (0.009) |
| **1-3 years** (weighted by teachers) | 42,815 | 50,308 | 7,493 | **0.161** (0.010) |
| **1-5 years** (weighted by teachers) | 43,298 | 50,641 | 7,343 | **0.156** (0.007) |
| **1-9 years** (weighted by teachers) | 44,319 | 51,991 | 7,672 | **0.159** (0.005) |
| All experience (composition differs) | 49,890 | 63,120 | 13,230 | 0.235 |

| Survey year (survey) | log(M/B), 1-5 years (SE) | dollars, 1-5 years | log(M/B), first year (SE) |
|---|---:|---:|---:|
| 2011-12 (SASS) | 0.139 (0.012) | 5,838 | 0.115 (0.030) |
| 2015-16 (NTPS) | 0.160 (0.008) | 7,068 | 0.175 (0.020) |
| 2017-18 (NTPS) | 0.156 (0.007) | 7,343 | 0.127 (0.019) |
| 2020-21 (NTPS) | 0.152 (0.007) | 7,647 | 0.115 (0.023) |

The master's premium at the same experience is 0.13 in the first year and 0.16 over years 1-5, and the 1-5 year figure is stable across surveys (0.14-0.16); the 2011-12 row is from the Schools and Staffing Survey (SASS), the others from the NTPS. It is below the ACS hourly-wage premium of teachers aged 25-34 with a graduate degree over BA-only teachers (0.195 (0.0043) age-controlled, 0.225 (0.0042) raw), which also contains professional and doctoral degrees, district composition and hours; at a given age the graduate teachers also have fewer years of experience, which works the other way.

## 3. NLSY79: ability, BA status and teachers

Block A (population sample, N = 9,033), persons observed at 25-34 (survey years 1982-1998; 8,776 persons). Ability is `afqt_c`; a gap in means divided by sqrt(R_comp_pop) = 0.957 is a gap in latent SDs. BA+ = highest grade completed >= 16 at some round at 25-34 (the extract has no degree variable; the convention of `nlsy_teachers.py`). 1979 weights, one observation per person, household-cluster bootstrap with the draws of `nlsy_teachers.py`. Failed draws: 0.

| Sample | Persons | BA+ share | Gap BA+ minus non-BA (SD), pooled | men | women |
|---|---:|---:|---:|---:|---:|
| all persons observed at 25-34, highest grade >= 9 (the ACS Table 2 universe): primary | 8,528 | 0.258 (0.0071) | +1.227 (0.0197) | +1.298 (0.0313) | +1.157 (0.0286) |
| all persons observed at 25-34, no grade rule | 8,776 | 0.253 (0.0070) | +1.259 (0.0195) | +1.332 (0.0314) | +1.188 (0.0286) |
| persons with an employed round at 25-34, highest grade >= 9 | 8,102 | 0.263 (0.0073) | +1.208 (0.0201) | +1.286 (0.0312) | +1.130 (0.0294) |
| persons with an employed round at 25-34 (the 'employed' sample of nlsy_teachers.py) | 8,319 | 0.259 (0.0072) | +1.238 (0.0199) | +1.319 (0.0313) | +1.156 (0.0294) |
| employed non-K-12-teachers, highest grade >= 9 | 7,942 | 0.249 (0.0073) | +1.225 (0.0207) | +1.293 (0.0311) | +1.150 (0.0313) |
| employed non-K-12-teachers (comparison group A of nlsy_teachers.py) | 8,159 | 0.245 (0.0072) | +1.255 (0.0205) | +1.327 (0.0312) | +1.177 (0.0311) |

The primary sample (first row) is the NLSY counterpart of the ACS Table 2 universe. The gap is large (about 1.2 latent SDs) and insensitive to the sample; at the ACS BA share (0.339 against 0.258 (0.0071) here) a logit level shift of the same gradient gives +1.219 (0.0188), so the cohort difference in BA attainment does not move it.

**BA status by AFQT tercile and the logit slope.** The observed score has reliability R = 0.916, so statistics of its terciles and of the logit of BA on the observed score are attenuated relative to the latent ability that the model's z represents. Three forms are reported. *Observed* (the data as they are): terciles of the weighted sample, pooled over gender (cut points -0.403 and 0.494 observed SDs), at the NLSY79 BA share. *At the ACS share*: the fitted observed-score logit with its level shifted so that the BA share is the ACS's (0.339; the NLSY79 cohort has fewer BAs). *Latent*: the logit of BA on the latent factor, estimated by maximum likelihood with a normal factor and reliability R (Gauss-Hermite integration of the logistic over the factor given the score; assumes a logistic P(BA | factor) and a normal factor), with P(BA | tercile of the factor) at both shares.

| P(BA+ | tercile) | Observed, pooled | men | women | Logit fit to the observed, pooled | Observed logit, at ACS share | Latent terciles | Latent, at ACS share |
|---|---:|---:|---:|---:|---:|---:|---:|
| T1 (lowest) | 0.041 (0.0040) | 0.037 (0.0054) | 0.045 (0.0055) | 0.033 (0.0026) | 0.057 (0.0036) | 0.027 (0.0025) | 0.052 (0.0036) |
| T2 | 0.167 (0.0094) | 0.153 (0.0132) | 0.179 (0.0126) | 0.180 (0.0081) | 0.277 (0.0037) | 0.154 (0.0073) | 0.262 (0.0034) |
| T3 (highest) | 0.567 (0.0161) | 0.569 (0.0207) | 0.564 (0.0191) | 0.562 (0.0144) | 0.683 (0.0070) | 0.568 (0.0140) | 0.704 (0.0070) |
| T3 - T1 | +0.526 (0.0163) | +0.532 (0.0216) | +0.520 (0.0192) |  | +0.626 (0.0105) | +0.541 (0.0143) | +0.653 (0.0106) |

| Logit of BA+ on ability | Pooled | Men | Women |
|---|---:|---:|---:|
| Slope per observed AFQT SD (logit, observed score) | 1.81 (0.057) | 1.84 (0.086) | 1.80 (0.070) |
| Slope per latent SD, first-order correction (observed / sqrt(R)) | 1.90 (0.059) | 1.93 (0.090) | 1.88 (0.073) |
| Slope per latent SD, errors-in-variables logit (MLE) | 2.00 (0.068) | 2.04 (0.105) | 1.97 (0.084) |

The first-order correction understates the latent slope (1.90 against 2.00, 5% more) because measurement error also flattens the non-linearity of the logistic (the attenuation of a logit slope exceeds the linear factor sqrt(R)). The latent logit is a parametric deconvolution: it implies a latent gap of 1.274 against 1.227 measured directly (men 1.289 against 1.298, women 1.264 against 1.157), and the observed-score logit fits the observed terciles to within 0.013. **Model counterparts.** A model without score noise should be compared with the *latent* terciles and slope. To compare with the observed-scale statistics instead, add noise to the model's standardized log z (theta): x = sqrt(R) theta + sqrt(1 - R) u with u standard normal and R = 0.916, and form terciles of x (over adults) and the logit slope of P(d = 1) on x. The first columns are at the NLSY79 BA share; the model's BA share will be near the ACS's, so the 'at ACS share' columns (a level shift of the same logit) are the relevant levels.

**Teacher gaps recomputed (confirmation).** Same Frame, spec `c80-emp-all-maj` (employed, 1980 Census codes, all K-12 codes, majority of rounds), same bootstrap draws:

| Comparison (latent AFQT SDs) | Here, pooled | Published `nlsy_teachers.json`, pooled | men | women |
|---|---:|---:|---:|---:|
| A. teachers minus employed non-teachers | +0.614 (0.0616) | +0.614 (0.0616) | +0.641 (0.1668) | +0.619 (0.0670) |
| B. teachers minus non-teacher persons with 16+ years | -0.334 (0.0625) | -0.334 (0.0625) | -0.351 (0.1667) | -0.278 (0.0687) |
| C. teachers with 16+ years minus non-teachers with 16+ years | -0.238 (0.0617) | -0.238 (0.0617) | -0.281 (0.1637) | -0.178 (0.0670) |

Teachers have 16+ years of schooling in 0.918 (0.0225) of cases, employed non-teachers in 0.245 (0.0072). Point estimates and SEs equal the published ones (max abs. difference 3.3e-16 in values, 1.1e-13 in SEs).

## 4. The cost of a BA, in model units

Sources: published prices from NCES Digest 2018 Table 330.10 (<https://nces.ed.gov/programs/digest/d18/tables/dt18_330.10.asp>; average undergraduate tuition and required fees, dormitory and board charges of full-time students, 2017-18, FTE-weighted, in-state rates at public institutions; current = constant dollars in the base year); net prices from the College Board (*Trends in College Pricing 2018*, Figures 9 and 10, printed pp. 18-19: published price minus grant aid from all sources and federal education tax benefits, 2017-18, 2018 dollars; <https://research.collegeboard.org/media/pdf/trends-college-pricing-2018-full-report.pdf>) and NCES Digest 2019 Table 331.30 (<https://nces.ed.gov/programs/digest/d19/tables/dt19_331.30.asp>; cost of attendance minus grants for first-time full-time students awarded Title IV aid, which includes room, board, books and other expenses and excludes the 17% with no Title IV aid). Enrollment weights for the blend: fall 2017 undergraduates at 4-year institutions, Table 303.70 (<https://nces.ed.gov/programs/digest/d18/tables/dt18_303.70.asp>; 72.8% public among public and private nonprofit; for-profits, 6%, excluded). The College Board's published public in-state price ($10,270, 2018 dollars) is above NCES's ($9,037): different coverage and enrollment weights; the net prices are the College Board's own (re-read from the report at run time: 10 figures checked, all equal). A BA costs 4 years of the annual amount (5 years in the JSON).

**Units: two choices.** The model is a two-period OLG whose adult period stands for the whole working life; its flows were converted as annual dollars over CE mean annual expenditures ($61,224.13 per consumer unit, CE 2018 Table R-1, <https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf>, as in `housing_costs.md`).

1. *Per consumer unit or per person.* The CE unit has 2.5 people, 1.9 adults (2.5 minus 0.6 children under 18; the one-decimal inputs leave a range of about 1.8-2.0) and 1.3 earners. A housing cost is per household and scales with the household, so its ratio to expenditures per unit is unit-free (`housing_costs.md`, *Units*). Tuition is paid by one student and does not scale with household size. The model's agent is one adult (every person of the Table 2 universe, home production included, is in the model), so the consistent denominator is expenditures **per adult**, $61,224 / 1.9 = $32,223. Per earner ($47,095) is the alternative if only workers count; it gives 0.68 times the per-adult F_c. The first release divided by the per-unit figure, which understates F_c by the factor 1.9 for a one-adult model household.
2. *Annual or one-off.* A *recurring* annual cost such as the residence cost r_p enters as annual dollars over annual expenditures (column i). A *one-off* goods cost paid from adult consumption (the education cost e and the moving cost m are paid this way, `spatial_teachers_model_eps0.md` section 2, *Timing*) is spread over the 40 working years: the consistent counterpart is **(ii) total cost over working-life consumption**. Undiscounted, total / (40 x annual consumption) (ii-a), is the counterpart in a model without an interest rate and with the other education parameters calibrated on annual aggregate flows (eta from education spending over labor income; the schooling time share s from years over 25): in a stationary economy an annual aggregate flow ratio equals the undiscounted lifetime ratio per person. Discounting the 40 years at 3% (ii-b, total / (23.115 x annual consumption), the annuitized debt service as a share of annual consumption) or also compounding the four college payments to the start of working life (ii-c, factor 4.309 instead of 4) adds the cost of financing and gives 1.73 and 1.86 times (ii-a); it is the right figure only if the model's period is meant to carry the interest on tuition. Column (i) charges the cost every year of adult life and overstates a four-year cost by 40/4 = 10 undiscounted (23.1/4 = 5.8 at 3%). A single-year reading of the adult period would reverse this; it is not the model's.

Per adult (denominator $32,223):

| Measure (2017-18) | Annual $ | 4-year total $ | (i) annual / CE | (ii-a) total / 40 CE | (ii-b) at 3% | (ii-c) compounded, at 3% |
|---|---:|---:|---:|---:|---:|---:|
| published tuition and required fees, public 4-year (in-state) | 9,037 | 36,148 | 0.2804 | 0.0280 | 0.0485 | 0.0523 |
| published tuition and required fees, all 4-year institutions | 15,924 | 63,696 | 0.4942 | 0.0494 | 0.0855 | 0.0921 |
| published tuition and required fees, private nonprofit 4-year | 33,755 | 135,020 | 1.0475 | 0.1048 | 0.1813 | 0.1953 |
| published tuition, fees, room and board, public 4-year (in-state) | 20,050 | 80,200 | 0.6222 | 0.0622 | 0.1077 | 0.1160 |
| published tuition, fees, room and board, all 4-year institutions | 27,357 | 109,428 | 0.8490 | 0.0849 | 0.1469 | 0.1583 |
| net tuition and fees, public 4-year in-state (after grants and tax benefits) | 3,780 | 15,120 | 0.1173 | 0.0117 | 0.0203 | 0.0219 |
| net tuition and fees, private nonprofit 4-year | 14,500 | 58,000 | 0.4500 | 0.0450 | 0.0779 | 0.0839 |
| net tuition and fees, public and private nonprofit 4-year averaged with fall-2017 undergraduate enrollment weights (0.728 public) | 6,700 | 26,802 | 0.2079 | 0.0208 | 0.0360 | 0.0388 |
| net price (cost of attendance minus grants), public 4-year, Title IV aid recipients | 13,460 | 53,840 | 0.4177 | 0.0418 | 0.0723 | 0.0779 |
| net price (cost of attendance minus grants), all 4-year, Title IV aid recipients | 18,120 | 72,480 | 0.5623 | 0.0562 | 0.0973 | 0.1048 |

Denominator sensitivity for the central measure (net tuition and fees, enrollment-weighted):

| Denominator | (i) annual / CE | (ii-a) total / 40 CE | (ii-b) at 3% | (ii-c) compounded, at 3% |
|---|---:|---:|---:|---:|
| per consumer unit (the first release) | 0.1094 | 0.0109 | 0.0189 | 0.0204 |
| per adult (central) | 0.2079 | 0.0208 | 0.0360 | 0.0388 |
| per earner | 0.1423 | 0.0142 | 0.0246 | 0.0265 |

The same for the central measure and the two ends of the range, (ii-a) / (ii-b):

| Measure | CU (ii-a) | CU (ii-b) | Adult (ii-a) | Adult (ii-b) | Earner (ii-a) | Earner (ii-b) |
|---|---:|---:|---:|---:|---:|---:|
| net tuition and fees, public and private nonprofit 4-year averaged with fall-2017 undergraduate enrollment weights (0.728 public) [central] | 0.0109 | 0.0189 | 0.0208 | 0.0360 | 0.0142 | 0.0246 |
| net tuition and fees, public 4-year in-state (after grants and tax benefits) | 0.0062 | 0.0107 | 0.0117 | 0.0203 | 0.0080 | 0.0139 |
| published tuition and required fees, all 4-year institutions | 0.0260 | 0.0450 | 0.0494 | 0.0855 | 0.0338 | 0.0585 |

**Aggregate-flow cross-check (unit-free).** The CE records $798.09 of college tuition per consumer unit in 2018 (SE 55.80), 1.30% of expenditures (SE 0.09 points, ignoring the covariance). In a stationary economy this share equals P(d = 1) x F_c, so with the ACS BA+ share 0.339 it implies F_c = 0.0384 (SE about 0.0027) of consumption, with no per-unit or per-adult conversion and no discounting. It is larger than the central BA-only figure because it counts all postsecondary tuition households pay (graduate and professional schools, two-year colleges, non-completers) over the BA+ share, and close to the annuitized per-adult figure (0.0360); it supports a range of about 0.02-0.04 per adult.

**Recommended.** The education goods cost is tuition and fees net of grants, not room and board (living costs are paid whether or not one attends; the model has no young-age consumption). Central F_c = **0.0208** of C-bar per adult (net tuition and fees, enrollment-weighted public/private nonprofit, four years, (ii-a), undiscounted: $26,802 over 40 x $32,223); 0.0360 annuitized at 3% (ii-b), 0.0142 / 0.0246 per earner; range 0.0117 (public 4-year net, undiscounted) to 0.0855 (published tuition, all 4-year, 3%). The choice between undiscounted and annuitized is a judgment about whether the model's period carries the interest on tuition, and the level is small next to the 0.48 log-point BA wage premium: tuition is a minor part of the cost of a BA next to the time cost (the s margin) and any utility cost c_C. Note that eta is calibrated on total education spending, 6.0% of GDP in 2018 (NCES Digest 2021 Table 605.20, government and private expenditures on educational institutions at all levels, postsecondary included, <https://nces.ed.gov/programs/digest/d21/tables/dt21_605.20.asp>): if eta is left unchanged, adding F_c counts the college part of that spending twice (F_c = 0 with eta unchanged is the alternative).

## 5. Moments for the model

`education_targets.toml` carries each moment with its definition, SE, a suggested criterion scale for targets, the model counterpart and the name of the model-side moment (`model_moment`). Roles follow DESIGN.md M3: *target* (BA share of non-teachers, BA wage premium in its two forms, NLSY79 AFQT gap, s_O mean), alternative forms by gender or with controls, and *validation* (schooling ratio, teacher premia for each teacher group and comparison, teacher z gaps, BA status by AFQT tercile, the latent and observed logit slopes). `ba_premium` and `ba_premium_within_occupation` (and the age-controlled form) are alternatives of one moment (`choice_group`).

| Key | Role | Value (SE) | Units | Aggregation | Model moment |
|---|---|---:|---|---|---|
| `ba_share_nonteacher` | target | 0.3231 (0.0007) | share | T2 | ba_share_nonteach |
| `ba_premium` | target | 0.4806 (0.0012) | log points | T2 | ba_premium |
| `ba_premium_within_occupation` | target | 0.3610 (0.0013) | log points | P | ba_premium |
| `afqt_gap_ba` | target | 1.2272 (0.0197) | latent AFQT SDs | weighted means | ba_z_gap_sd |
| `ba_share_nonteacher_m` | target_by_gender | 0.2938 (0.0007) | share | T2 | ba_share_nonteach |
| `ba_premium_m` | target_by_gender | 0.4795 (0.0019) | log points | T2 | ba_premium |
| `ba_premium_within_occupation_m` | target_by_gender | 0.3257 (0.0021) | log points | P | ba_premium |
| `afqt_gap_ba_m` | target_by_gender | 1.2982 (0.0313) | latent AFQT SDs | weighted means | ba_z_gap_sd |
| `ba_share_nonteacher_f` | target_by_gender | 0.3522 (0.0010) | share | T2 | ba_share_nonteach |
| `ba_premium_f` | target_by_gender | 0.5081 (0.0016) | log points | T2 | ba_premium |
| `ba_premium_within_occupation_f` | target_by_gender | 0.3957 (0.0020) | log points | P | ba_premium |
| `afqt_gap_ba_f` | target_by_gender | 1.1574 (0.0286) | latent AFQT SDs | weighted means | ba_z_gap_sd |
| `ba_share_all` | validation | 0.3391 (0.0011) | share | P | ba_share |
| `ba_premium_age_adj` | target_alt | 0.4896 (0.0013) | log points | P | ba_premium |
| `s_O_mean` | target | 0.5457 | share of 25 years | workbook | s_O_mean |
| `s_ratio` | validation | 1.1761 | ratio | workbook | s_ratio |
| `s_ba_nonteacher` | validation | 0.6515 (<0.0001) | share of 25 years | T2 | none (mean s of d = 1 non-teachers) |
| `s_nonba_nonteacher` | validation | 0.4971 (0.0001) | share of 25 years | T2 | none (mean s of d = 0 non-teachers) |
| `s_T_mean_ba` | validation | 0.6573 (0.0001) | share of 25 years | P | s_T_mean_ba |
| `teacher_ba_share` | validation | 0.8744 (0.0017) | share | P | ba_share_teach |
| `teacher_graduate_share` | validation | 0.3779 (0.0021) | share | P | none |
| `teacher_premium_all` | validation | 0.0786 (0.0026) | log points | T2 | teacher_premium_m/f |
| `teacher_premium_all_m` | validation | 0.0557 (0.0049) | log points | T2 | teacher_premium_m/f |
| `teacher_premium_all_f` | validation | 0.0860 (0.0031) | log points | T2 | teacher_premium_m/f |
| `teacher_premium_ba` | validation | -0.2211 (0.0028) | log points | T2 | teacher_premium_ba_m/f |
| `teacher_premium_ba_m` | validation | -0.2610 (0.0049) | log points | T2 | teacher_premium_ba_m/f |
| `teacher_premium_ba_f` | validation | -0.2081 (0.0033) | log points | T2 | teacher_premium_ba_m/f |
| `teacher_premium_ba_teachers_all` | validation | 0.1383 (0.0024) | log points | T2 | teacher_premium_m/f when every teacher has d = 1 |
| `teacher_premium_ba_teachers_all_m` | validation | 0.0793 (0.0052) | log points | T2 | teacher_premium_m/f when every teacher has d = 1 |
| `teacher_premium_ba_teachers_all_f` | validation | 0.1574 (0.0029) | log points | T2 | teacher_premium_m/f when every teacher has d = 1 |
| `teacher_premium_within_ba` | validation | -0.1614 (0.0026) | log points | T2 | teacher_premium_ba_m/f when every teacher has d = 1 |
| `teacher_premium_within_ba_m` | validation | -0.2374 (0.0053) | log points | T2 | teacher_premium_ba_m/f when every teacher has d = 1 |
| `teacher_premium_within_ba_f` | validation | -0.1368 (0.0032) | log points | T2 | teacher_premium_ba_m/f when every teacher has d = 1 |
| `teacher_premium_core_all` | validation | 0.1472 (0.0023) | log points | T2 | none (robustness) |
| `teacher_premium_core_all_m` | validation | 0.0794 (0.0054) | log points | T2 | none (robustness) |
| `teacher_premium_core_all_f` | validation | 0.1691 (0.0028) | log points | T2 | none (robustness) |
| `teacher_premium_core_ba` | validation | -0.1524 (0.0024) | log points | T2 | none (robustness) |
| `teacher_premium_core_ba_m` | validation | -0.2373 (0.0054) | log points | T2 | none (robustness) |
| `teacher_premium_core_ba_f` | validation | -0.1250 (0.0030) | log points | T2 | none (robustness) |
| `teacher_ba_wage_gap` | validation | 0.4819 (0.0094) | log points | P | none (pay penalty of a non-BA teaching track) |
| `teacher_ba_wage_gap_m` | validation | 0.2017 (0.0219) | log points | P | none (pay penalty of a non-BA teaching track) |
| `teacher_ba_wage_gap_f` | validation | 0.5664 (0.0095) | log points | P | none (pay penalty of a non-BA teaching track) |
| `grad_teacher_premium` | validation | 0.1953 (0.0043) | log points | P | none (compare with ma_premium_implied) |
| `ma_ba_salary_premium` | validation | 0.1564 (0.0068) | log points | delta method | ma_premium_implied |
| `teacher_z_gap_employed_nonteachers` | validation | 0.6141 (0.0616) | latent AFQT SDs | weighted means | teacher_z_gap_market_sd |
| `teacher_z_gap_nonteacher_ba` | validation | -0.3335 (0.0625) | latent AFQT SDs | weighted means | teacher_z_gap_ba_sd |
| `teacher_z_gap_within_ba` | validation | -0.2384 (0.0617) | latent AFQT SDs | weighted means | teacher_z_gap_within_ba_sd |
| `p_ba_afqt_tercile_1` | validation | 0.0408 (0.0040) | share | weighted means | ba_share_by_z |
| `p_ba_afqt_tercile_2` | validation | 0.1671 (0.0094) | share | weighted means | ba_share_by_z |
| `p_ba_afqt_tercile_3` | validation | 0.5668 (0.0161) | share | weighted means | ba_share_by_z |
| `p_ba_afqt_top_minus_bottom` | validation | 0.5260 (0.0163) | share difference | weighted means | ba_share_by_z |
| `afqt_logit_slope` | validation | 1.8148 (0.0572) | log-odds per observed AFQT SD | weighted means | ba_share_by_z (slope form) |
| `afqt_logit_slope_latent` | validation | 1.9986 (0.0685) | log-odds per latent AFQT SD | EIV logit | ba_share_by_z (slope form) |
| `tuition_Fc` | cost | 0.0208 | units of Cbar | (ii-a) | Fc |

## Caveats

- **Sampling error is not the error that matters.** ACS replicate-weight SEs are 0.001-0.005 on shares and premia; the conventions (T2 against P, comparison group, wage-sample rules) move the numbers by 0.005-0.01 (the weeks rule by more, section 1.4), and the existing targets use a 0.05 modelling scale for premia. The suggested scales in the TOML follow that practice. SEs treat the Table 2 counts as fixed.
- **Degrees against years.** ACS BA+ is a degree (SCHL >= 21); NLSY79 BA+ is 16+ completed years (no degree variable in the extract). The NLSY79 cohort (born 1957-64, observed 1982-98) had a BA share of 0.258 (0.0071) against 0.339 (0.0011) in the ACS 2009-13; the gap is not sensitive to this (section 3), the tercile shares are (hence the 'at ACS share' columns).
- **AFQT timing and measurement.** The ASVAB was taken at ages 15-23 and partly reflects schooling already attained, so the BA gap overstates the ability difference at the schooling decision (the caveat of `nlsy_teachers.md`); latent SDs assume classical error and the same error variance in every group; terciles of the observed score flatten the latent gradient. The latent logit slope and terciles are model-based (a logistic P(BA | factor), a normal factor and the composite reliability R) and not direct data moments.
- **Teachers.** The K-12 definition (codes 2300-2340) includes pre-K/kindergarten and other instructors; 'teaching requires a BA' holds for the core codes only. The existing teacher premia and the teaching-share targets use the whole group; the BA+ and core-code rows of section 1.3 are the like-for-like counterparts of a model whose teachers all have d = 1. The +0.614 / -0.334 NLSY79 gaps rest on 160 teachers (24 men, 136 women).
- **Wages.** Hourly wages are annual wage and business income over usual hours x weeks (weeks bracketed to 49 or 51). The 48+ week rule is selective for teachers (section 1.4): it leaves out 20.6% of full-time K-12 teachers against 12.1% of other full-time workers, mostly teachers who report 40-47 weeks and whose implied hourly wage is higher for the same salary, so the teacher premia are probably biased down; the sensitivity is reported, and the rule is inherited from the workbook and the existing targets. The salary-schedule premium (base salary, same experience) and the ACS hourly premium (same age) are different objects.
- **NCES salary table.** Experience-specific SEs are the table's (survey design included); the log-ratio SEs assume independence; the aggregates weight rows by all-teacher counts because counts by degree are not published.
- **Costs.** Published prices are averages of sticker prices (in-state at public institutions); net prices average over students including those paying nothing; the IPEDS net price covers Title IV recipients only and includes living costs. The CE denominator is per consumer unit (2.5 people, 1.9 adults, 1.3 earners; one-decimal inputs, so adults per unit are known to about +/-5%); the central F_c divides by adults, the per-earner and per-unit figures are in section 4. 2017-18 academic-year and 2018 calendar-year dollars differ by about 1-2% and are not adjusted. The 3% rate and 40 working years are the task's assumptions. BA duration is the nominal four years; completers take longer and non-completers also pay. As for the other goods costs, mapping the CE's expenditures per adult to the model's C-bar is a unit convention, not a calibration of C-bar.
- **Untouched sample.** The NLSY97 counterparts are reserved for the sealed hold-out (DESIGN.md D4) and were not computed here.

Written with help from Claude Code.
