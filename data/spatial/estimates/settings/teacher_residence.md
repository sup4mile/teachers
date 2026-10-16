# Teacher residence versus workplace: evidence for the M2 split

Generated 2026-10-05T02:50:20 by `teacher_residence.py` (DESIGN.md section M2, M2_SPEC_v2.md). Residence-based teachers: ACS 5-year 2014-18 school-district table C24010; workplace-based teachers: the audited K-12 FTE of `district_settings.py`; locale partition, 188 CZs, whole-CZ bootstrap B = 1999, seed 20260929 (the draws of `district_settings.py`). Composition and commuting: ACS 2018 1-year PUMS (replicate-weight standard errors) with the Census 2010 PUMA equivalency files and the OMB September-2018 principal-city list. EDU is the occupation group 'Educational instruction, and library occupations' (SOC 25-0000, education, training and library). Sign convention: suburb minus city, in logs. **Names.** The model's `gap_teach_res_emp` and `net_commute_gap` are about K-12 teachers; the district table can only give the EDU group, so the district gaps are EDU proxies (`gap_edu_res_emp`, `net_commute_gap_edu`) and the K-12 counterparts are built in section 3 (`[k12]` of the TOML). Residents are counted once.

## Summary

- **K-12 counterparts of the model moments (188 CZs; approximate).** The ACS district table holds the whole EDU group, of which public K-12 teachers are the share s (PUMS). Adding the log ratio of s in suburban to city PUMAs of the same CZ (Delta log s = 0.271, survey SE 0.017, 118 CZs with both classes, 93% of the weight) to the EDU district gaps gives **`gap_teach_res_emp` = 0.298 (0.041)** (public K-12 teachers per employed resident; range over specifications [0.18, 0.37]) and **`net_commute_gap` = 0.075 (0.038)** (K-12 residents per K-12 FTE; range [-0.05, 0.14]); the K-12 residents' gap is 0.456 (0.120) against a K-12 FTE gap of 0.381 (0.114). The SE adds the CZ-bootstrap variance of the district gap and the composition term (together, correlation -0.24), the PUMS survey variance and the ACS sampling variance; it is conservative. Cross-checks without the district table: public K-12 teachers per employed resident, suburban against city PUMAs of the same CZ, 0.317 (0.041); national residence against workplace of public K-12 teachers on POWPUMA classes, 0.108 (0.008) (by residence 0.787 (0.016), by workplace 0.679 (0.016)). Reading: more public K-12 teachers per resident worker live in the suburb (positive in every specification, [0.18, 0.37]); the net commuting gap is small and positive in the headline (more teachers live in the suburb than work there: a net flow of teachers from suburban homes to city schools) but its sign is not secure across specifications [-0.05, 0.14].
- **EDU proxies (not the model counterparts).** `gap_edu_res_emp` (EDU residents / employed residents) 0.028 (0.025); `gap_edu_res` 0.185 (0.112); `gap_fte` (K-12 FTE) 0.381 (0.114); **`net_commute_gap_edu` = `gap_edu_res` - `gap_fte` = -0.196 (0.030)** (95% interval [-0.257, -0.138]), the log gap of EDU residents per K-12 FTE: 3.99 in the city, 3.18 in the suburb. Its negative sign is composition: public K-12 teachers are 31.0% (0.4) of EDU workers in principal-city PUMAs and 38.7% (0.2) in suburban PUMAs (postsecondary teachers 19.6% (0.2) and 13.0% (0.2)). The four adopted targets reproduce (gap_pupils 0.3764, gap_salary 0.0086, gap_teachers_pp -0.0177, gap_income 0.1122); the bootstrap correlations of the new moments with them are in section 2 (`net_commute_gap_edu` with `gap_teachers_pp` 0.28, `gap_salary` 0.09).
- **Residents are counted once.** The earnings rule (each LEA takes its own layer's record) counts residents of an elementary district inside a secondary district twice when both sit in one location (employed residents overcounted by 8.6% in the matched sample, more in the suburbs). The headline uses the unified and elementary records, which tile the states, and keeps the secondary LEAs' FTE and pupils; the earnings rule would give `net_commute_gap_edu` -0.153 (0.035) and `gap_edu_res` 0.228 (0.114) (bias +0.043), residents from unified + secondary records -0.199 (0.030); the ratio `gap_edu_res_emp` is insensitive (0.029 (0.025) under the earnings rule).
- **Commuting (PUMS, commuters 25-64).** K-12 teachers commute 22.9 (0.1) minutes one way against 28.2 (0.0) for other workers; 18.6% (0.2) commute more than 30 minutes against 27.1% (0.1) (30.2% (0.3) against 41.3% (0.1) for 30 minutes or more: reported times heap at exactly 30, 11.6% (0.2) and 14.2% (0.0)); 20.7% (0.2) work outside the POWPUMA they live in against 27.1% (0.1), and 8.0% (0.2) of public K-12 teachers who live and work in a city or suburban POWPUMA cross between the two classes. Working in another POWPUMA adds **17.1 minutes** (SE 0.3) one way for teachers in metropolitan areas (residence-PUMA fixed effects) and working in a POWPUMA of the other class, city or suburb, **16.0** (0.5; residence and workplace classified on the same POWPUMA geography). Both measure crossings of a county line, so they may overstate the increment within a county; Nashville teachers who live out of the county commute 13.8 minutes longer (Santelli and Grissom 2024).
- **Model units.** Extra commute 16.6 minutes one way (range 13.8-20.7), 185 days, value of an hour $14.41 (USDOT: 50% of hourly household income; range 35%-60%): an annual time cost of $1,473 (102 hours). **r_k per adult = 0.0457 (range 0.0259 to 0.0704)** in units of mean consumption per adult (CE 2018 $61,224.13 per consumer unit / 1.9 adults = $32,223); per CE consumer unit (the unit of r_p) it is 0.0241 (0.0136 to 0.0370). **tau_k** (share of income lost, M2_SPEC_v2) is **0.023** for teachers (range 0.013-0.035; 185 days) and 0.034 (0.019-0.045) for other workers (250 days, 18.1 extra minutes). Time cost only (section 5).

## 1. Data, definitions and units

- **Residents (a).** ACS 5-year 2014-18 table C24010, 'Sex by Occupation for the Civilian Employed Population 16 Years and Over' (`teacher_residence.py fetch`, school-district layers unified, elementary and secondary, one state per request as in `acs_earnings.py`): C24010_001 total, C24010_014 (men) + C24010_050 (women) 'Educational instruction, and library occupations' (the 2018 vintage's label for SOC 25-0000, Education, Training, and Library; vintages through 2017 call it 'Education, training, and library occupations'). Annotation codes (negative) are masked (none occur); 'Remainder of <state>' rows (code 99999) are dropped; one record per `leaid`, as `district_bootstrap.load_income_columns` does for B20005. 2,534 unified, 518 elementary and 152 secondary LEAs of the partition have their own layer's record; 56 have none. **Residents are counted once** (headline): the unified and elementary records tile each state (section 2) and the secondary layer overlays them, so the secondary LEAs contribute their FTE and pupils but no residents; the earnings rule of the B20005 moments (every LEA its own record) is the variant `earnings_rule`.
- **Workplace (b).** `teachers_k12` (CCD teacher FTE minus pre-K teacher FTE) and `pupils_k12` (pupils minus pre-K enrollment) of the audited panel of `district_settings.py`: the FTE variant 11 behind `gap_teachers_pp`; New York City is one district (the collapsed row carries the Chancellor's id 3620580, under which the ACS files New York City). Regular districts only (CCD agency types 1-2), so charter LEAs are not in the FTE.
- **Matched sample.** Every district moment uses the districts with a C24010 EDU count, a C24010 employed total, an audited K-12 FTE and K-12 pupils (and, for the secondary LEAs of the headline, an FTE and pupils): 3,199 of 3,260 partition districts, 99.4% of partition pupils, so that `net_commute_gap_edu` is exactly the log gap of EDU residents per FTE on one set of districts. Each statistic on its own maximal sample is the `unmatched` variant in section 2.
- **Moments.** As the audited targets: within each CZ, location totals (or ratios of totals) summed over the districts of the location, then log(suburb / city); pupil-weighted mean over CZs (weights = the CZ's partitioned pupils); whole-CZ bootstrap with `district_bootstrap.draw_counts` (B = 1999, seed 20260929), SE = SD of the draws, interval = 2.5-97.5 percentiles. The draws are the ones `district_settings.py` uses, so the correlations in section 2 can weight the new moments with the four existing targets.
- **PUMS (d, e).** ACS 2018 1-year PUMS through the Census API, one state at a time (`acs_occupations.py`'s `Api`): employed civilians (ESR 1-2) 16+ in households, with OCCP, SOCP, INDP, COW, PUMA, POWPUMA, POWSP, JWMNP, JWTR, AGEP, SEX, SCHL, wage variables and the 80 replicate weights; 1,474,819 records, 155.7 million weighted. Standard errors are the Census successive-difference replication, SE = sqrt(4/80 sum_r (x_r - x)^2), applied to every statistic including differences and regressions. K-12 teachers are OCCP 2310, 2320 and 2330 (elementary and middle, secondary and special-education teachers; SOC 25-2020, 25-2030, 25-2050); adding preschool and kindergarten teachers (2300) is a sensitivity (the PUMS cannot separate kindergarten from preschool, and the FTE count kindergarten teachers). 'Public' is a local or state government employer (COW 3-4). Pipeline check: the national mean travel time to work of workers 16+ who do not work at home is 27.15 minutes (published: 27.1) and 5.2% work at home (published: 5.3%), all persons, ACS 2018 1-year (Census Bureau, commuting guidance 'Mean travel time to work', www2.census.gov/programs-surveys/commuting/guidance/acs-1yr/Mean-travel-time.pdf).
- **PUMA classes.** 2010 PUMAs from the Census PUMA equivalency files (`PUMSEQ10_<ST>.txt`: 2010 population of each place in each PUMA) and the OMB September-2018 list of principal cities of *metropolitan* statistical areas (list 2; the NCES city locale is territory in a principal city of an urbanized area). pc_share = population in principal-city places / PUMA population; **city** if pc_share >= 0.50; **suburb** if at least 50% of the population lives in counties of metropolitan areas (list 1) and pc_share is below the cut; **other** otherwise. 689 of the 690 principal cities of metropolitan areas outside Puerto Rico are found in the files (1 unmatched: Macon-Bibb County (state 13, place 49008; Macon-Bibb County, GA)). 733 city, 1273 suburban and 345 other PUMAs hold 96, 168 and 45 million people (2010). POWPUMAs (county based; the PUMA to POWPUMA conversion of the 2018 PUMS code lists) take the population-weighted shares of the PUMAs they aggregate: 146 city, 523 suburban and 313 other. Spot checks: the PUMAs of New York City, Chicago, Milwaukee, Hartford and Boston are all city; those of Nassau, Suffolk and Westchester (NY), Lake (IL) and Oakland (MI) counties are suburb (White Plains, 35% principal city, is suburb); DuPage County's six PUMAs are suburb and a seventh, which holds Naperville, a principal city, is city.
- **Units.** r_k is a goods cost per period in units of mean consumption: the annual time cost of the extra commute divided by CE 2018 mean annual expenditures, $61,224.13 per consumer unit (BLS Table R-1, the denominator of `housing_costs.py`; 2.5 persons, 0.6 children under 18, so 1.9 adults, BLS Table 1400) or, because the model's agents are adults and the cost is per worker, per adult, $32,223. tau_k is the same cost over the worker's annual earnings. The value of an hour is the USDOT local-personal rate (50% of hourly median household income, $13.60 in 2015 dollars) uprated by the CPI-U ratio 251.107 / 237.017 = 1.0594 to 2018 dollars.

## 2. Districts: residents against workplace (a)-(c)

Locale partition (city = NCES locale 11-13, suburb = 21-23; districts of >= 500 pupils; each location >= 10% of the CZ's partitioned pupils): 188 CZs, 3,260 districts, 28.55 million pupils (effective number of CZs 42, the five largest carry 26% of the weight). Headline matched sample (residents counted once): 3,199 districts (602 city, 2597 suburban), 99.4% of the pupils; the 61 others have no ACS record (56 under the earnings rule) or no audited FTE or pupil count (6). In the headline sample: 5.76 million EDU residents, 94.7 million employed residents, 1.69 million K-12 FTE teachers (57% of the 2.97 million audited K-12 FTE of all regular districts; 152 secondary LEAs contribute FTE and pupils only).

### EDU moments (suburb minus city, log gaps)

| moment | what | value | boot SE | ACS sampling SE | 95% interval | n CZs |
|---|---|---|---|---|---|---|
| `gap_edu_res_emp` | EDU proxy for gap_teach_res_emp (NOT the K-12 model counterpart): Delta log (EDU residents / employed residents) | 0.0277 | 0.0247 | 0.0040 | [-0.0215, 0.0736] | 188 |
| `gap_edu_res` | Delta log EDU residents (workers by residence) | 0.1851 | 0.1120 | 0.0041 | [-0.0186, 0.4190] | 188 |
| `gap_fte` | Delta log K-12 FTE teachers (workers by workplace), matched sample | 0.3807 | 0.1138 |  | [0.1662, 0.6052] | 188 |
| `net_commute_gap_edu` | EDU proxy for net_commute_gap (NOT the K-12 model counterpart): gap_edu_res - gap_fte = Delta log (EDU residents per K-12 FTE) | -0.1955 | 0.0301 | 0.0041 | [-0.2566, -0.1378] | 188 |
| `gap_edu_per_pupil` | Delta log (EDU residents per K-12 pupil) | -0.2140 | 0.0360 |  | [-0.2888, -0.1463] | 188 |
| `gap_emp_res` | Delta log employed residents | 0.1575 | 0.1139 |  | [-0.0588, 0.3867] | 188 |
| `gap_pupils_k12` | Delta log K-12 pupils, matched sample | 0.3992 | 0.1180 |  | [0.1769, 0.6355] | 188 |
| `gap_fte_per_pupil` | Delta log (K-12 FTE per K-12 pupil), matched sample (the gap_teachers_pp definition) | -0.0185 | 0.0131 |  | [-0.0444, 0.0072] | 188 |

`net_commute_gap_edu` equals `gap_edu_res` - `gap_fte` CZ by CZ (largest absolute difference 3.6e-16), so it is also the log gap of EDU residents per K-12 FTE. `gap_fte` = `gap_pupils_k12` + `gap_fte_per_pupil` (the matched-sample pupils gap is 0.3992 against the tracked `gap_pupils` 0.3764), and the matched `gap_fte_per_pupil` -0.0185 is the adopted `gap_teachers_pp` (-0.0177) on the districts that also have an ACS record. The ACS sampling SE is computed from the published margins of error of C24010 (90% MOE / 1.645, men's and women's EDU counts added, location counts by the root-sum-of-squares rule, independence across districts and CZs, delta method for the log gaps; FTE are exact counts): small next to the bootstrap SE, which already contains most of it.

### Levels (CZ-weighted means of the location values, on the CZs that have the gap)

| level | city | suburb |
|---|---|---|
| EDU residents per K-12 FTE teacher | 3.99 (0.16) | 3.18 (0.09) |
| EDU residents / employed residents | 0.0603 (0.0015) | 0.0616 (0.0015) |
| EDU residents per K-12 pupil | 0.2443 (0.0119) | 0.1913 (0.0068) |
| K-12 FTE per K-12 pupil | 0.0620 (0.0021) | 0.0610 (0.0023) |

EDU residents per FTE is about 3-4: the group includes postsecondary teachers, preschool teachers, teaching assistants, librarians and private-school staff that are not in the public K-12 FTE; only the gaps are comparable with the model.

### Bootstrap correlation with the four adopted targets

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` | `k12_gap_teach_res_emp` | `k12_net_commute_gap` | `k12_gap_k12_res` | `gap_edu_res_emp` | `gap_edu_res` | `gap_fte` | `net_commute_gap_edu` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `gap_pupils` | 1.00 | -0.12 | -0.37 | 0.28 | -0.01 | -0.06 | 0.94 | -0.20 | 0.95 | 0.99 | -0.23 |
| `gap_salary` | -0.12 | 1.00 | 0.41 | 0.15 | -0.08 | -0.06 | -0.09 | 0.08 | -0.05 | -0.07 | 0.09 |
| `gap_teachers_pp` | -0.37 | 0.41 | 1.00 | -0.07 | 0.52 | 0.35 | -0.16 | 0.63 | -0.20 | -0.27 | 0.28 |
| `gap_income` | 0.28 | 0.15 | -0.07 | 1.00 | 0.10 | 0.47 | 0.41 | 0.18 | 0.44 | 0.29 | 0.54 |
| `k12_gap_teach_res_emp` | -0.01 | -0.08 | 0.52 | 0.10 | 1.00 | 0.76 | 0.27 | 0.77 | 0.13 | 0.06 | 0.27 |
| `k12_net_commute_gap` | -0.06 | -0.06 | 0.35 | 0.47 | 0.76 | 1.00 | 0.27 | 0.67 | 0.19 | -0.01 | 0.73 |
| `k12_gap_k12_res` | 0.94 | -0.09 | -0.16 | 0.41 | 0.27 | 0.27 | 1.00 | 0.06 | 0.98 | 0.96 | 0.02 |
| `gap_edu_res_emp` | -0.20 | 0.08 | 0.63 | 0.18 | 0.77 | 0.67 | 0.06 | 1.00 | 0.03 | -0.13 | 0.62 |
| `gap_edu_res` | 0.95 | -0.05 | -0.20 | 0.44 | 0.13 | 0.19 | 0.98 | 0.03 | 1.00 | 0.96 | 0.07 |
| `gap_fte` | 0.99 | -0.07 | -0.27 | 0.29 | 0.06 | -0.01 | 0.96 | -0.13 | 0.96 | 1.00 | -0.19 |
| `net_commute_gap_edu` | -0.23 | 0.09 | 0.28 | 0.54 | 0.27 | 0.73 | 0.02 | 0.62 | 0.07 | -0.19 | 1.00 |

`k12_*` are the K-12 counterparts of section 3 (district gap and composition term resampled together). `gap_edu_res` and `gap_fte` correlate 0.9 or more with `gap_pupils`: both are size gaps. The net gaps are nearly uncorrelated with the size gap and correlate with `gap_teachers_pp` and `gap_salary`, the information that diagonal weights discard. The full draws are in `teacher_residence_draws.csv`.

### Overlap of the ACS layers and other variants

Layers in the earnings-rule matched sample (share of employed residents): unified 84.2%, elementary 7.9%, secondary 7.9%; pupils: unified 90.1%, elementary 6.3%, secondary 3.1%. 28 CZs, carrying 42% of the weight, have a secondary LEA. Summing every layer's record counts 1.086 times the residents of the unified and elementary records. In the CZ locations that hold a secondary LEA the elementary and secondary layers hold nearly the same number of employed residents (city 2.79 against 2.89 million, suburb 5.10 against 5.24; the elementary layer is 3.5% and 2.7% smaller), so counting either layer once leaves at most a small coverage gap (elementary districts below the 500-pupil screen or outside the partition).

The unified and elementary records tile the states: their employed residents are 99.98% of the state totals of the same table, and the secondary layer adds 6.7% nationally (Montana 89%, Illinois 38%, Arizona 33%, California 27%). The state sums of EDU residents (9.277 million) match the PUMS count (9.279 million).

| variant | `gap_edu_res_emp` | `gap_edu_res` | `gap_fte` | `net_commute_gap_edu` | n CZs |
|---|---|---|---|---|---|
| **headline: residents counted once (unified + elementary records); FTE and pupils from every LEA** | 0.028 (0.025) | 0.185 (0.112) | 0.381 (0.114) | -0.196 (0.030) | 188 |
| residents once from unified + secondary records (elementary residents dropped where a secondary LEA is present) | 0.026 (0.025) | 0.182 (0.113) | 0.381 (0.114) | -0.199 (0.030) | 188 |
| earnings rule: every LEA takes its own layer's record (overlapping residents counted twice) | 0.029 (0.025) | 0.228 (0.114) | 0.381 (0.114) | -0.153 (0.035) | 188 |
| unified districts only (elementary and secondary LEAs leave every sum) | 0.021 (0.024) | 0.097 (0.123) | 0.313 (0.121) | -0.216 (0.030) | 183 |
| New York City dropped from both the resident and the FTE / pupil counts | 0.014 (0.022) | 0.432 (0.218) | 0.644 (0.232) | -0.212 (0.033) | 188 |
| New York CZ left out | 0.015 (0.023) | 0.227 (0.111) | 0.423 (0.112) | -0.196 (0.032) | 187 |
| each statistic on its own maximal sample |  | 0.186 (0.112) | 0.374 (0.113) | -0.188 (0.030) | 188 |

The earnings rule suits ratios and medians, not counts: it adds 8.6% to the employed residents, more in the suburbs (secondary records add 6.6% to the city's and 10.3% to the suburb's employed residents in the matched sample), so it inflates the suburban count: it moves `net_commute_gap_edu` by +0.043 (1.4 SE) and `gap_edu_res` by +0.043; the two single-count constructions agree within 0.004. The ratio `gap_edu_res_emp` moves by +0.0016 (numerator and denominator are overcounted alike). The largest change in `net_commute_gap_edu` across the variants is 0.043 (1.4 SE; `earnings_rule`). The unified-only sample drops 5 CZs; leaving New York out or dropping it from both counts changes the net gap by 0.000 and 0.016.

### Reproduction check

| target | adopted (value, SE) | this script (value, SE) | n CZs |
|---|---|---|---|
| `gap_pupils` | 0.3764 (0.1154) | 0.3764 (0.1154) | 188 |
| `gap_salary` | 0.0086 (0.0092) | 0.0086 (0.0092) | 188 |
| `gap_teachers_pp` | -0.0177 (0.0131) | -0.0177 (0.0131) | 188 |
| `gap_income` | 0.1122 (0.0168) | 0.1122 (0.0168) | 188 |

All four match to the printed digits: **yes** (standard errors checked: yes).

## 3. Composition of SOC 25-0000 (d) and the K-12 counterparts

ACS 2018 1-year PUMS, employed 16+ in households: 9.28 million workers (95,593 records) in education, training and library occupations, 4.83 million of them K-12 teachers (OCCP 2310, 2320, 2330), 3.42 million in the public sector.

| occupation (OCCP; SOCP in the records) | workers (millions) | share of EDU (SE) |
|---|---|---|
| 2205 postsecondary teachers (SOCP 251000) | 1.381 | 14.88% (0.14) |
| 2300 preschool and kindergarten teachers (SOCP 252010) | 0.586 | 6.31% (0.11) |
| 2310 elementary and middle school teachers (SOCP 252020) | 3.826 | 41.23% (0.19) |
| 2320 secondary school teachers (SOCP 252030) | 0.751 | 8.10% (0.11) |
| 2330 special education teachers (SOCP 252050) | 0.254 | 2.73% (0.08) |
| 2350 tutors (SOCP 253041) | 0.185 | 2.00% (0.06) |
| 2360 other teachers and instructors (SOCP 2530XX) | 0.698 | 7.53% (0.12) |
| 2400 archivists, curators and museum technicians (SOCP 254010) | 0.053 | 0.57% (0.03) |
| 2435 librarians and media collections specialists (SOCP 254022) | 0.175 | 1.88% (0.06) |
| 2440 library technicians (SOCP 254031) | 0.038 | 0.41% (0.02) |
| 2545 teaching assistants (SOCP 259040) | 1.215 | 13.10% (0.17) |
| 2555 other educational instruction and library workers (SOCP 2590XX) | 0.117 | 1.26% (0.04) |
| **K-12 teachers (2310 + 2320 + 2330)** | 4.831 | **52.06% (0.20)** |
| K-12 plus preschool and kindergarten (+ 2300) | 5.416 | 58.37% (0.20) |
| public-sector K-12 (COW 3-4) | 3.421 | **36.87% (0.21)** |
| public K-12 in elementary and secondary schools (INDP 7860) | | 36.74% (0.21) |

As a check on the count, the PUMS public K-12 teachers (headcount) are 3.42 million against 2.97 million audited FTE of regular districts.

### By class of residence PUMA

Share of EDU workers (by residence PUMA) in each occupation group; SE in parentheses.

| | city | suburb | other | metropolitan | suburb - city | log ratio suburb / city |
|---|---|---|---|---|---|---|
| EDU workers' distribution over classes | 32.0% (0.2) | 55.8% (0.2) | 12.2% (0.1) | | | |
| K-12 teachers | 46.4% (0.4) | 54.3% (0.3) | 56.9% (0.5) | 51.4% (0.2) | 7.9 pp (0.5) | 0.158 (0.009) |
| K-12 + preschool / kindergarten | 52.4% (0.3) | 60.9% (0.3) | 62.5% (0.5) | 57.8% (0.2) | 8.6 pp (0.5) | 0.152 (0.008) |
| public K-12 teachers | 31.0% (0.4) | 38.7% (0.2) | 43.9% (0.6) | 35.9% (0.2) | 7.7 pp (0.4) | 0.221 (0.013) |
| public K-12 + preschool / kindergarten | 31.8% (0.4) | 39.5% (0.2) | 45.1% (0.6) | 36.7% (0.2) | 7.7 pp (0.4) | 0.218 (0.013) |
| 2205 postsecondary teachers | 19.6% (0.2) | 13.0% (0.2) | 11.1% (0.3) | | -6.6 pp (0.3) | -0.410 (0.019) |
| 2300 preschool and kindergarten teachers | 6.0% (0.2) | 6.6% (0.1) | 5.6% (0.3) | | 0.6 pp (0.2) | 0.102 (0.034) |
| 2310 elementary and middle school teachers | 36.6% (0.4) | 43.0% (0.3) | 45.1% (0.5) | | 6.4 pp (0.5) | 0.161 (0.012) |
| 2320 secondary school teachers | 7.5% (0.2) | 8.3% (0.1) | 8.8% (0.3) | | 0.8 pp (0.2) | 0.096 (0.029) |
| 2330 special education teachers | 2.2% (0.1) | 3.0% (0.1) | 3.0% (0.2) | | 0.8 pp (0.1) | 0.296 (0.059) |
| 2350 tutors | 2.4% (0.1) | 1.9% (0.1) | 1.2% (0.1) | | -0.5 pp (0.1) | -0.220 (0.059) |
| 2360 other teachers and instructors | 8.0% (0.2) | 7.7% (0.2) | 5.7% (0.2) | | -0.3 pp (0.2) | -0.036 (0.032) |
| 2400 archivists, curators and museum technicians | 0.9% (0.1) | 0.4% (0.0) | 0.4% (0.1) | | -0.5 pp (0.1) | -0.763 (0.105) |
| 2435 librarians and media collections specialists | 1.9% (0.1) | 1.8% (0.1) | 2.1% (0.2) | | -0.1 pp (0.1) | -0.052 (0.053) |
| 2440 library technicians | 0.3% (0.0) | 0.4% (0.0) | 0.5% (0.1) | | 0.1 pp (0.0) | 0.264 (0.138) |
| 2545 teaching assistants | 13.1% (0.3) | 12.7% (0.2) | 15.1% (0.4) | | -0.4 pp (0.4) | -0.028 (0.028) |
| 2555 other educational instruction and library workers | 1.5% (0.1) | 1.1% (0.1) | 1.4% (0.1) | | -0.4 pp (0.1) | -0.294 (0.077) |

The city's EDU workforce is tilted toward postsecondary teachers (colleges and universities are in principal cities) and away from elementary and secondary teachers, so residents' EDU gaps understate the suburb's teacher residents.

**Is the PUMA classification defensible?** It uses only public Census files (the 2010 population of every place within every PUMA and the OMB list of principal cities), uses the principal-city criterion that underlies the NCES city locale and passes the spot checks of section 1. It is not the NCES locale partition of the districts: a PUMA holds 100,000 people or more, so PUMAs near the cut mix city and suburb, and PUMAs in outlying counties of a metropolitan area count as suburb whether or not they are in the urbanized area. The composition result does not depend on it in sign: the log ratio of the public K-12 share is 0.20-0.28 on national classes and 0.25-0.34 inside CZs across the five cuts of the table below.

### Dispersion across PUMAs

Share of a PUMA's EDU workers who are K-12 teachers, 2,351 PUMAs with on average 41 sampled EDU workers (minimum 4), so a PUMA's share has a sampling SD of about 0.09. 'Net of sampling' subtracts the average replicate-weight sampling variance of the PUMA shares from the variance across PUMAs (weights: the PUMA's EDU workers).

| | K-12 teachers | public K-12 teachers |
|---|---|---|
| mean (EDU-weighted) | 0.521 (0.002) | 0.369 (0.002) |
| SD across PUMAs (EDU-weighted) | 0.138 (0.014) | 0.133 (0.014) |
| SD net of sampling noise (EDU-weighted) | 0.102 (0.019) | 0.098 (0.019) |
| SD across PUMAs (unweighted) | 0.136 (0.016) | 0.135 (0.015) |
| SD net of sampling noise (unweighted) | 0.092 (0.023) | 0.094 (0.021) |
| 10th percentile | 0.335 (0.020) | 0.194 (0.019) |
| 25th percentile | 0.435 (0.015) | 0.275 (0.010) |
| median | 0.528 (0.005) | 0.370 (0.005) |
| 75th percentile | 0.618 (0.011) | 0.461 (0.009) |
| 90th percentile | 0.690 (0.022) | 0.539 (0.019) |

The composition of EDU workers differs across areas by about 0.1 (one SD, net of sampling error) around a mean share of 0.52, of the same order as the 8 point city / suburb difference in the K-12 share, so one area's EDU count is a noisy proxy for its teachers; averages over the CZs of the partition are less noisy.

### The K-12 counterparts of the model moments (approximate)

Inside each CZ of the partition the PUMAs are assigned to the CZ that holds most of their 2010 population (Dorn's crosswalk, as in the district panel); the share s of EDU workers who are public K-12 teachers is computed for the city and the suburban PUMAs of the CZ (class cut as above), log(suburb / city) per CZ, and averaged with the CZ pupil weights over the 118 CZs with PUMAs of both classes and a positive share in both (93.5% of the weight; 64,665 sampled EDU workers). Since EDU residents = public K-12 residents / s, Delta log (public K-12 residents) = Delta log (EDU residents) + Delta log s: the K-12 moments are the headline EDU moments of section 2 plus this term (the district gaps of all CZs plus the mean of the term over the CZs that have it).

Standard errors. The district gap and the composition term are resampled TOGETHER over whole CZs with the draws of `district_settings.py` (the two are correlated: -0.24 for the net gap; the CZ-bootstrap SE of the composition term is 0.023 against its replicate-weight survey SE 0.017); the survey variance of the composition term and the ACS sampling variance of C24010 are added. The CZ bootstrap already contains part of that sampling noise (each CZ's PUMS and ACS estimates are noisy and the bootstrap resamples them), so the total SE is conservative; `se_boot` alone is the lower figure. Adding the district SE and the survey SE in quadrature instead would ignore the CZ variance of the composition term and its covariance with the district gap.

| moment | what | value | SE | bootstrap | survey | ACS | 95% interval | range over specifications | n CZs |
|---|---|---|---|---|---|---|---|---|---|
| `gap_teach_res_emp` | public K-12 teacher residents / employed residents (gap_teach_res_emp) | **0.2984** | 0.0407 | 0.0365 | 0.0174 | 0.0040 | [0.22, 0.38] | [0.18, 0.37] (15 specifications) | 188 |
| `net_commute_gap` | public K-12 teacher residents per K-12 FTE (net_commute_gap) | **0.0752** | 0.0379 | 0.0334 | 0.0174 | 0.0041 | [0.00, 0.15] | [-0.05, 0.14] (15 specifications) | 188 |
| `gap_k12_res` | public K-12 teacher residents (gap_edu_res + composition) | **0.4559** | 0.1197 | 0.1183 | 0.0174 | 0.0041 | [0.22, 0.69] | [0.33, 0.52] (15 specifications) | 188 |

Value of each K-12 moment under the specifications (district gap of the variant plus the composition term of the specification):

| specification | `net_commute_gap` | `gap_teach_res_emp` | `gap_k12_res` |
|---|---|---|---|
| **headline: residents counted once, public K-12 / EDU within CZs, PUMA classes 0.5 / 0.5** | 0.075 | 0.298 | 0.456 |
| earnings rule: every LEA takes its own layer's record (residents counted twice) | 0.118 | 0.300 | 0.499 |
| residents once from unified + secondary records | 0.072 | 0.297 | 0.452 |
| unified districts only | 0.055 | 0.292 | 0.367 |
| New York City dropped from both counts (not in the range: the New York CZ keeps its weight with almost no city location left, so its gap is dominated by what remains; `drop_ny_cz` is the clean version) | 0.059 | 0.285 | 0.703 |
| New York CZ left out | 0.075 | 0.285 | 0.498 |
| composition: all-sector K-12 teachers / EDU | -0.011 | 0.212 | 0.370 |
| composition: public K-12 plus preschool / kindergarten (OCCP 2300) / EDU | 0.070 | 0.293 | 0.450 |
| composition: national city / suburb PUMA split (not CZ by CZ) | 0.026 | 0.249 | 0.407 |
| composition: POWPUMA classes (county based, national) | -0.046 | 0.177 | 0.335 |
| composition: group quarters included (as in C24010) | 0.094 | 0.317 | 0.474 |
| composition: PUMA class cut pc >= 0.25, metro >= 0.5 | 0.058 | 0.281 | 0.438 |
| composition: PUMA class cut pc >= 0.75, metro >= 0.5 | 0.142 | 0.365 | 0.523 |
| composition: PUMA class cut pc >= 0.9, metro >= 0.5 | 0.142 | 0.366 | 0.523 |
| composition: PUMA class cut pc >= 0.5, metro >= 0.9 | 0.074 | 0.297 | 0.455 |
| district gap plus composition added CZ by CZ on the CZs that have both | 0.070 | 0.301 | 0.477 |

Cross-checks that do not use the district table: **`gap_teach_res_emp` measured directly in the PUMS** (public K-12 teachers / all employed residents, suburban against city PUMAs of the same CZ, CZ-weighted, 118 CZs) is 0.317 (0.041) (survey and CZ-bootstrap variance added); on the national POWPUMA classes it is 0.184 (0.016); the national net commuting gap of public K-12 teachers on POWPUMA classes (section 4) is 0.108 (0.008). The composition term with group quarters included (C24010 counts them, the household PUMS sample does not): 0.289 against 0.271; 1.2% of EDU workers live in group quarters (1.7% in city PUMAs, 0.8% in suburban ones).

The composition term and the extra commute of working in a POWPUMA of the other class under other cuts of the PUMA classes (pc_share is the share of a PUMA's 2010 population living in principal-city places; the metropolitan cut is the share in metropolitan counties):

| city if pc_share >= | suburb if metro share >= | PUMAs city / suburb | Delta log s, national classes | Delta log s, inside the CZs (n CZs) | K-12 `net_commute_gap` | extra minutes, other class of POWPUMA | share crossing |
|---|---|---|---|---|---|---|---|
| 0.25 | 0.50 | 990 / 1019 | 0.199 (0.014) | 0.253 (0.017) (131) | 0.058 | 15.33 (0.42) | 10.2% |
| 0.50 (main) | 0.50 | 733 / 1273 | 0.221 (0.013) | 0.271 (0.017) (118) | 0.075 | 16.04 (0.47) | 8.5% |
| 0.75 | 0.50 | 542 / 1464 | 0.268 (0.014) | 0.337 (0.020) (98) | 0.142 | 19.00 (0.83) | 4.4% |
| 0.90 | 0.50 | 432 / 1574 | 0.284 (0.013) | 0.338 (0.020) (82) | 0.142 | 20.74 (0.96) | 3.2% |
| 0.50 | 0.90 | 733 / 1158 | 0.216 (0.013) | 0.269 (0.018) (114) | 0.074 | 15.79 (0.50) | 8.6% |

## 4. Commuting (e)

ACS 2018 1-year PUMS, employed 25-64 in households. Commuters are workers with a positive one-way travel time (JWMNP > 0) who do not work at home; 'outside' is a place of work in another POWPUMA (or state) than the POWPUMA that contains the worker's residence PUMA (Census PUMA to POWPUMA conversion of the 2018 PUMS code lists), among commuters with a US place of work. Teachers: OCCP 2310, 2320, 2330. Reported times heap at multiples of five: 11.6% (0.2) of teachers and 14.2% (0.0) of other workers report exactly 30 minutes, so 'more than 30' and '30 or more' differ by that much; the literature figure (Santelli and Grissom, 19% of Nashville teachers) is a '30 or more' share.

| group | n | mean one-way minutes | share > 30 min | share >= 30 min | share working at home | share working outside the residence POWPUMA |
|---|---|---|---|---|---|---|
| K-12 teachers | 46,103 | 22.94 (0.12) | 18.6% (0.2) | 30.2% (0.3) | 1.5% (0.1) | 20.7% (0.2) |
|   women | 35,270 | 22.45 (0.13) | 17.6% (0.3) | 29.1% (0.3) | 1.6% (0.1) | 19.8% (0.2) |
|   men | 10,833 | 24.45 (0.25) | 21.5% (0.5) | 33.6% (0.6) | 1.4% (0.1) | 23.3% (0.4) |
|   public-sector | 32,979 | 22.58 (0.15) | 18.0% (0.3) | 29.3% (0.4) | 1.1% (0.1) | 19.9% (0.3) |
| other workers | 1,159,207 | 28.20 (0.03) | 27.1% (0.1) | 41.3% (0.1) | 5.4% (0.0) | 27.1% (0.1) |
|   women | 544,648 | 26.44 (0.04) | 24.3% (0.1) | 38.0% (0.1) | 6.1% (0.0) | 24.4% (0.1) |
|   men | 614,559 | 29.70 (0.04) | 29.4% (0.1) | 44.1% (0.1) | 4.9% (0.0) | 29.3% (0.1) |
|   with a bachelor's degree or more | 441,221 | 29.62 (0.04) | 30.6% (0.1) | 44.6% (0.1) | 7.7% (0.1) | 30.4% (0.1) |
|   public-sector (non-teachers) | 120,189 | 26.25 (0.08) | 24.1% (0.2) | 37.2% (0.2) | 2.0% (0.1) | 23.0% (0.2) |
| K-12 teachers, metropolitan PUMAs | 38,096 | 23.61 (0.13) | 19.4% (0.3) | 31.6% (0.3) | 1.5% (0.1) | 21.4% (0.2) |
| other workers, metropolitan PUMAs | 972,023 | 28.72 (0.03) | 27.8% (0.1) | 42.5% (0.1) | 5.5% (0.0) | 27.3% (0.1) |
| K-12 teachers living in city PUMAs | 11,922 | 23.35 (0.22) | 17.9% (0.4) | 30.3% (0.6) | 1.5% (0.1) | 15.6% (0.3) |
| K-12 teachers living in suburban PUMAs | 26,174 | 23.73 (0.15) | 20.1% (0.3) | 32.2% (0.4) | 1.5% (0.1) | 24.2% (0.3) |
| other workers, city PUMAs | 345,611 | 27.53 (0.05) | 24.7% (0.1) | 39.7% (0.1) | 5.4% (0.0) | 19.4% (0.1) |
| other workers, suburban PUMAs | 626,412 | 29.42 (0.04) | 29.6% (0.1) | 44.2% (0.1) | 5.6% (0.0) | 32.0% (0.1) |

Teachers minus other workers (SE):

| sample | mean minutes | share > 30 min | share >= 30 min | share at home | share outside POWPUMA |
|---|---|---|---|---|---|
| all | -5.27 (0.12) | -8.5 pp (0.2) | -11.1 pp (0.3) | -3.9 pp (0.1) | -6.4 pp (0.2) |
| metropolitan PUMAs | -5.11 (0.13) | -8.4 pp (0.3) | -11.0 pp (0.3) | -4.0 pp (0.1) | -5.9 pp (0.2) |
| living in city PUMAs | -4.18 (0.22) | -6.8 pp (0.4) | -9.4 pp (0.6) | -3.9 pp (0.1) | -3.8 pp (0.3) |
| living in suburban PUMAs | -5.69 (0.15) | -9.5 pp (0.3) | -12.0 pp (0.4) | -4.1 pp (0.1) | -7.8 pp (0.3) |
| vs workers with a bachelor's degree or more | -6.69 (0.11) | -12.0 pp (0.2) | -14.4 pp (0.3) | -6.2 pp (0.1) | -9.8 pp (0.2) |

### Extra commute from working elsewhere

Weighted regression of one-way minutes on the crossing indicator, sex, age and age squared with residence-PUMA fixed effects (workers who live in the same PUMA are compared); replicate-weight SEs. Two definitions of crossing: another POWPUMA (or state) than the one that contains the residence PUMA, and a POWPUMA of the other class (city / suburb) than the one that contains the residence PUMA (both ends on the POWPUMA geography; both ends city or suburb).

| sample | crossing | extra minutes (SE) | n | n crossing | share crossing |
|---|---|---|---|---|---|
| K-12 teachers | other POWPUMA or state | 17.51 (0.28) | 43,338 | 8,875 | 20.7% |
| K-12 teachers, metropolitan PUMAs | other POWPUMA or state | 17.11 (0.31) | 35,738 | 7,597 | 21.4% |
| public K-12 teachers, metropolitan | other POWPUMA or state | 17.10 (0.40) | 25,319 | 5,251 | 20.7% |
| women teachers, metropolitan | other POWPUMA or state | 16.30 (0.33) | 27,235 | 5,577 | 20.5% |
| men teachers, metropolitan | other POWPUMA or state | 19.43 (0.75) | 8,503 | 2,020 | 24.2% |
| teachers living in suburban PUMAs | other POWPUMA or state | 17.31 (0.37) | 24,603 | 5,856 | 24.2% |
| teachers living in city PUMAs | other POWPUMA or state | 16.53 (0.71) | 11,135 | 1,741 | 15.6% |
| other workers | other POWPUMA or state | 19.43 (0.08) | 1,073,580 | 296,969 | 27.1% |
| other workers, metropolitan | other POWPUMA or state | 18.83 (0.09) | 898,494 | 250,580 | 27.3% |
| other workers with a bachelor's degree, metropolitan | other POWPUMA or state | 20.12 (0.11) | 362,485 | 112,747 | 30.8% |
| K-12 teachers living and working in a city or suburban POWPUMA | POWPUMA of the other class (city / suburb) | 16.04 (0.47) | 35,671 | 2,968 | 8.5% |
| other workers living and working in a city or suburban POWPUMA | POWPUMA of the other class (city / suburb) | 17.42 (0.12) | 894,476 | 107,424 | 12.0% |
| K-12 teachers, as above, residence PUMAs whose own class equals their POWPUMA's | POWPUMA of the other class (city / suburb) | 16.90 (0.52) | 28,639 | 2,486 | 8.9% |
| K-12 teachers, SUPERSEDED: class of the residence PUMA against the class of the work POWPUMA | mixed geographies (see below) | 12.17 (0.53) | 35,066 | 8,419 | 24.3% |

Without controls, K-12 teachers in metropolitan PUMAs commute 36.4 minutes when they work outside the residence POWPUMA and 20.1 when they work inside: 16.3 minutes (SE 0.32); other workers 18.0. **Bracketing.** Classifying the residence by its PUMA and the workplace by its POWPUMA (two geographies; the SUPERSEDED row of the table) codes 24% of metropolitan teachers as crossing against 9% on one geography, most of them working in their own residence POWPUMA, and biases the estimate down (12.2 against 16.0 minutes); it is not used. Both definitions used now are crossings of a county line (POWPUMAs are whole counties or groups of counties): the PUMS cannot see a move between a central city and its suburbs inside one county, and a county-line crossing is on average longer than a district-to-district move inside a county, so both can OVERSTATE the increment of working in the other location. Santelli and Grissom's out-of-county figure for Nashville (13.8 minutes) anchors the low end of the range.

### Residence against workplace of K-12 teachers on one geography (national)

Both residence and workplace are classified by the POWPUMA (the finest area on which the PUMS reports the workplace), all employed persons 16+ with a US POWPUMA; shares of all workers, so the log gaps are those of counts.

| | residence: city | suburb | other | workplace: city | suburb | other | gap by residence | gap by workplace | net commuting gap | gap, residents per employed resident | works in the other class | city <-> suburb |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| K-12 teachers | 0.86% | 1.79% | 0.38% | 0.93% | 1.72% | 0.38% | 0.731 (0.013) | 0.611 (0.014) | **0.120 (0.008)** | 0.128 (0.013) | 9.9% (0.1) | 8.4% (0.2) |
| public K-12 teachers | 0.58% | 1.27% | 0.30% | 0.62% | 1.23% | 0.29% | 0.787 (0.016) | 0.679 (0.016) | **0.108 (0.008)** | 0.184 (0.016) | 9.4% (0.2) | 8.0% (0.2) |

EDU workers by residence against public K-12 teachers by workplace (the PUMS version of the district measure): gap by residence 0.637 (0.009), by workplace 0.679 (0.016), net -0.042 (0.013); of which composition (EDU against public K-12 residents) -0.150 (0.012). The sign agrees with the district measure and the size is smaller; the classes are national and county-based.

Where the teachers who work in a class of area live (public K-12 teachers; share of those working in the class):

| work in | live in city | live in suburb | live in other |
|---|---|---|---|
| city POWPUMAs | 83.6% (0.4) | 15.3% (0.4) | 1.1% (0.1) |
| suburb POWPUMAs | 4.1% (0.2) | 94.0% (0.2) | 1.9% (0.1) |

**Check against the IBO count for New York City.** Public K-12 teachers (COW 3-4) with a workplace in the five NYC POWPUMAs (93 thousand weighted, 873 records): **72.7% (1.6) live in the NYC PUMAs** (74.0% (1.3) for all K-12 teachers, 76.5% (0.2) for all commuters) against 67% of the 70,328 teachers in the IBO count for 2011-12; 51.2% (1.9) live in the borough group (POWPUMA) where they work. The PUMS counts persons who report a teaching occupation (aides and substitutes included) in 2018, the IBO the Department of Education's teachers in schools, so the two differ by a few points.

### Literature

| source | what it shows |
|---|---|
| [Boyd, Lankford, Loeb and Wyckoff (2005), 'The Draw of Home: How Teachers' Preferences for Proximity Disadvantage Urban Schools', Journal of Policy Analysis and Management 24(1): 113-132 (NBER WP 9953, Sept. 2003)](https://www.nber.org/papers/w9953) | New York State, teachers entering public school teaching 1999-2002: 61% (60.8%) took their first job within 15 miles of their hometown and 85% within 40 miles (Table 2; 90% within 40 miles in New York City); urban districts import teachers from the suburbs: 36% of urban teachers have a suburban hometown, 6% of suburban teachers an urban hometown; the ratio of suburban-hometown teachers taking a city job to city-hometown teachers taking a suburban job runs from 2.6 (Buffalo) to 7.5 (New York City) (p. 7) (read: NBER working paper PDF, text and Tables 2-5) |
| [Reininger (2012), 'Hometown Disadvantage? It Depends on Where You're From', Educational Evaluation and Policy Analysis 34(2): 127-145](https://eric.ed.gov/?id=EJ963419) | NELS:88: the median distance between the high school a teacher attended in 10th grade and her residence eight years later is 13 miles, about a quarter of the distance moved by other professionals; 4 miles for teachers who grew up in dense cities (read: not read directly: the figure is as summarised at https://www.deansforimpact.org/about/news-and-blog/2017/03/06/13-miles-the-inherent-localism-of-teaching) |
| [New York City Independent Budget Office, 'Where Do NYC's Teachers and Principals Live Compared With Where They Work?'](https://www.ibo.nyc.ny.us/iboreports/printnycbtn19.pdf) | school year 2011-12: 67% of 70,328 New York City teachers live within the five boroughs; 81% of Staten Island teachers live in the borough where they work, 57% of Brooklyn teachers; teachers in Queens and the Bronx are the most likely to live outside the five boroughs (read: the PDF could not be opened from this environment (the old IBO address now redirects to the IBO home page, ibo.nyc.gov); figures and the date (Aug. 15, 2013) as given in the search summary of the report) |
| [Santelli and Grissom (2024), 'A Bad Commute: Travel Time to Work Predicts Teacher Turnover and Other Workplace Outcomes', AERA Open 10; working paper: EdWorkingPaper No. 22-691 (Dec. 2022), 'A Bad Commute: Does Travel Time to Work Predict Teacher and Leader Turnover and Other Workplace Outcomes?'](https://edworkingpapers.com/ai22-691) | Metropolitan Nashville Public Schools (countywide), 2008-09 to 2019-20, home and school addresses with HERE travel times: mean one-way commute 21.8 minutes for teachers (p. 10); 19% of teachers travel 30 minutes or more (p. 15; compare the PUMS share of 30 minutes or more, section 4, not the share above 30); teachers who live out of the county commute 13.8 minutes longer than those who live in it (p. 15); working paper: teachers who commute 45 minutes or more each way are 10 percentage points more likely to transfer schools (overall rate 12%) and 3 percentage points more likely to exit the district (overall rate 7%) than teachers in the same school who commute 0-5 minutes (pp. 16-17); predicted transfer probability 8% against 18%, exit 8% against 11% (p. 22). Published version: each 5 minutes of one-way commute raises the transfer probability by 0.8-1.0 percentage points over most of the distribution (read: EdWorkingPaper No. 22-691 (Dec. 2022), PDF text, pp. 10, 15, 16, 17, 22 (printed page numbers); the 0.8-1.0 percentage points per 5 minutes is as stated in the published AERA Open abstract (via search; the published article was not opened)) |
| [U.S. Department of Transportation (2016), 'The Value of Travel Time Savings: Departmental Guidance for Conducting Economic Evaluations, Revision 2 (2016 Update)', memorandum of Sept. 27, 2016](https://www.transportation.gov/sites/dot.gov/files/docs/2016%20Revised%20Value%20of%20Travel%20Time%20Guidance.pdf) | local personal travel (commuting included) is valued at 50% of hourly median household income ($56,516 / 2,080 = $27.20 in 2015), $13.60 an hour (Tables 1, 3, 4; text p. 11); plausible range 35%-60%, $9.50-$16.30 (Tables 2 and 5, pp. 14 and 19); business travel 100% of gross compensation, $25.40 (read: read from a web.archive.org copy (transportation.gov refuses scripted requests)) |
| [NCES, Digest of Education Statistics 2017, Table 203.90, 'Average daily attendance (ADA) as a percentage of total enrollment, school day length, and school year length in public schools, by school level and state: 2007-08 and 2011-12' (source: SASS Public School Data File, 2007-08 and 2011-12)](https://nces.ed.gov/programs/digest/d17/tables/dt17_203.90.asp) | United States average of 179 days in the school year (SE 0.1) and 1,203 hours, 2011-12 (the table reports school year length for 2011-12 only; 2007-08 has attendance and hours per day, 6.6, only) (read: web page (United States row)) |
| [NCES SASS / NTPS tabulations](https://nces.ed.gov/surveys/ntps/) | no NCES SASS / NTPS table of the share of public school teachers living in the district where they teach was found: the surveys report teacher mobility between schools and districts, not residence (read: searched; nothing found) |

## 5. The commuting cost: r_k and tau_k (f)

Annual time cost = 2 x (extra one-way minutes / 60) x days a year x (value of an hour). r_k = annual time cost / Cbar (M2_SPEC.md, a goods cost in units of mean consumption); tau_k = annual time cost / annual earnings (M2_SPEC_v2.md: income x (1 - tau_k) when r != w; the value of an hour is a share of the worker's own wage, so tau_k = share x 2 x extra hours x days / annual hours).

| input | low | central | high | source |
|---|---|---|---|---|
| extra one-way minutes of working in the other location, teachers | 13.8 | 16.6 | 20.7 | central = average of the two main teacher estimates in section 4 (crossing a POWPUMA 17.1, crossing to a POWPUMA of the other class 16.0; SE 0.35); low / high = the smallest and largest teacher estimate over the specifications (POWPUMA crossing, the POWPUMA-class crossing under five PUMA class cuts, Nashville); the 95% interval of the central is 15.1-17.7 |
| extra one-way minutes, other workers | 15.6 | 18.1 | 19.3 | the same two estimators on other workers 25-64 (18.8 and 17.4) and the class cuts |
| days a year, teachers | 180 | 185 | 190 | NCES: 179 school days (Digest 2017, Table 203.90, 2011-12), plus about a week of contract days (assumption) |
| days a year, other workers | 240 | 250 | 260 | 50 weeks of 5 days less holidays and days at home (assumption) |
| value of an hour (2018 dollars) | 10.06 | 14.41 | 17.27 | USDOT 2016: 35% / 50% / 60% of hourly median household income ($27.20 in 2015), x CPI-U 1.0594 |
| annual hours worked (full-time year-round workers 25-64) | | teachers 2,248, others 2,248 | | PUMS: usual weekly hours x weeks worked (51 for 50-52 weeks, 48.5 for 48-49), mean |
| Cbar per consumer unit | | 61,224.13 | | CE 2018 (BLS Table R-1; as r_p) |
| Cbar per adult | 34,013 | 32,223 | 30,612 | CE 2018 Table 1400: 2.5 persons less 0.6 children under 18 = 1.9 adults per consumer unit (one decimal: 1.8-2.0) |
| **r_k per adult** | **0.0259** | **0.0457** | **0.0704** | all-low and all-high combinations; central SE from the extra minutes alone 0.0010 |
| r_k per consumer unit | 0.0136 | 0.0241 | 0.0370 | the unit of r_p; central SE from the extra minutes alone 0.0005 |
| **tau_k, teachers** | **0.0129** | **0.0227** | **0.0351** | value of time 35% / 50% / 60% of the own wage; central SE from the extra minutes alone 0.0005; with 2,080 hours 0.0246 |
| **tau_k, other workers** | **0.0194** | **0.0336** | **0.0447** | same shares; central SE 0.0002; with 2,080 hours 0.0363 |

Central: 102 hours a year (33 extra minutes a day) x $14.41 = $1,473 a year: 2.4% of Cbar per consumer unit, 4.6% of Cbar per adult. For scale, full-time year-round K-12 teachers 25-64 earn a median $23.84 an hour on reported hours (a heap: 2% of them are exactly at it, so its SE is not informative; geometric mean $23.23, SE 0.08; other full-time year-round workers' median $22.35), so half of the teacher median is $11.92. The model's agents are adults, so the per-adult r_k is the one to use with the per-adult consumption scale (DESIGN_v2 'Units'); r_p (housing, `housing_costs.md`, 0.028) is a household cost over household spending and unit-free, a per-worker cost is not. M2_SPEC_v2 applies one tau_k to every worker: the other-worker value fits the 90%+ of the population that are not teachers; the teacher value fits teachers (fewer commuting days a year).

One-at-a-time sensitivity (other inputs at their central values; tau_k for teachers):

| change | extra minutes | days | $ an hour | r_k per adult | r_k per consumer unit | tau_k |
|---|---|---|---|---|---|---|
| extra minutes at its low end (smallest estimate) | 13.8 | 185 | 14.41 | 0.0381 | 0.0200 | 0.0189 |
| extra minutes at its high end (largest estimate) | 20.7 | 185 | 14.41 | 0.0572 | 0.0301 | 0.0284 |
| value of time at 35% of hourly household income | 16.6 | 185 | 10.06 | 0.0319 | 0.0168 | 0.0159 |
| value of time at 60% of hourly household income | 16.6 | 185 | 17.27 | 0.0548 | 0.0288 | 0.0273 |
| value of time at 50% of the teacher median hourly wage | 16.6 | 185 | 11.92 | 0.0378 | 0.0199 | 0.0227 |
| days 180 | 16.6 | 180 | 14.41 | 0.0445 | 0.0234 | 0.0221 |
| days 190 | 16.6 | 190 | 14.41 | 0.0469 | 0.0247 | 0.0234 |
| days 250 (a full-year worker) | 16.6 | 250 | 14.41 | 0.0618 | 0.0325 | 0.0307 |
| extra minutes: crossing a POWPUMA (teachers, metro) | 17.1 | 185 | 14.41 | 0.0472 | 0.0248 | 0.0235 |
| extra minutes: crossing a POWPUMA of the other class, city / suburb (teachers) | 16.0 | 185 | 14.41 | 0.0442 | 0.0233 | 0.0220 |
| extra minutes: Santelli-Grissom out-of-county, 13.8 | 13.8 | 185 | 14.41 | 0.0381 | 0.0200 | 0.0189 |
| central with the class-of-area minutes at PUMA class cut pc >= 0.25, metro >= 0.5 | 16.2 | 185 | 14.41 | 0.0447 | 0.0235 | 0.0223 |
| central with the class-of-area minutes at PUMA class cut pc >= 0.75, metro >= 0.5 | 18.1 | 185 | 14.41 | 0.0498 | 0.0262 | 0.0248 |
| central with the class-of-area minutes at PUMA class cut pc >= 0.9, metro >= 0.5 | 18.9 | 185 | 14.41 | 0.0522 | 0.0275 | 0.0260 |
| central with the class-of-area minutes at PUMA class cut pc >= 0.5, metro >= 0.9 | 16.5 | 185 | 14.41 | 0.0454 | 0.0239 | 0.0226 |
| adults per consumer unit 1.8 (CE counts are rounded to one decimal) | 16.6 | 185 | 14.41 | 0.0433 | 0.0241 | 0.0227 |
| adults per consumer unit 2.0 (CE counts are rounded to one decimal) | 16.6 | 185 | 14.41 | 0.0481 | 0.0241 | 0.0227 |

## 6. Caveats

- **The K-12 moments are approximate.** The ACS table has no K-12 detail at the school-district level; the group (EDU) also holds preschool and postsecondary teachers, aides, librarians and private-school staff, and the K-12 share differs by area (section 3). The K-12 versions add a PUMS composition term (a different survey, PUMA classes instead of the NCES locale partition, 2018 only, CZs assigned by population share) and are reported with a specification range ([0.18, 0.37] for `gap_teach_res_emp`, [-0.05, 0.14] for `net_commute_gap`): the sign of the net gap is not secure, the level of `gap_teach_res_emp` is. The EDU proxies must not be targeted against a K-12 teacher object. Charter-school teachers are in the residents but not in the FTE (charter LEAs are not regular districts). The PUMS K-12 count includes aides who report a teaching occupation, excludes kindergarten teachers (OCCP 2300 mixes them with preschool teachers; the FTE count kindergarten and drop pre-K; the `incl_prek` specification brackets it), and the household PUMS sample excludes group quarters while C24010 includes them (the `group_quarters` specification adds them).
- **Timing.** C24010 is ACS 2014-18 (occupation coded on the 2010 list for 2014-17, the 2018 list for 2018); FTE and pupils are fall 2018; PUMS is 2018 only. School-district boundaries are those of the ACS 2018 vintage.
- **Who is in the locations.** Only districts of the partition (>= 500 pupils, city or suburb locale) enter either count, so teachers who live in town or rural districts, in small districts or outside the CZ and work in the partition count in the FTE and not in the residents (and the reverse). The net gap is a net of all such flows, not just city-suburb commuting.
- **Overlap and NYC.** Residents are counted once (unified + elementary records); the earnings rule that counts overlapping residents twice is a variant (it biases the net gap upward by the amount in section 2). New York City enters with the Chancellor's record, matching the FTE, pupil and salary targets; the `nyc_dropped` variant drops it from both the resident and the FTE / pupil counts (the earnings moments drop NYC from the ACS-based statistics only).
- **Residents are persons, FTE are full-time equivalents.** Levels of EDU residents per FTE are not comparable with the model's teacher counts; only the gaps are.
- **POWPUMA geography.** POWPUMAs are whole counties or groups of counties, so the PUMS cannot see commuting between a central city and its suburbs inside one county, and both extra-commute estimates are crossings of a county line, which may overstate the increment of working in the other location (the Nashville out-of-county figure anchors the low end). PUMAs near the class cut are classified by a threshold (section 3 shows the cuts); the classes are not the NCES locale partition.
- **r_k and tau_k are time costs only.** Vehicle operating costs of the extra miles are not included. For scale (assumptions: 35 mph, 10-30 cents a mile of variable cost; the IRS 2018 standard business rate is 54.5 cents a mile, 25 of them depreciation, Notice 2018-03): 3,578 extra miles a year would add $358-$1,073, i.e. 0.006-0.018 per consumer unit and 0.011-0.033 per adult to r_k. The value of an hour is the USDOT guidance for personal travel (50% of household income), not an estimate for teachers; days a year are assumptions; the extra minutes are estimated on teachers and on other workers separately.
- **Units.** r_p is a cost per household over spending per household and is unit-free; a commuting cost is per worker, so it is divided by CE spending per ADULT (1.9 adults per consumer unit from the rounded published counts: 2.5 persons less 0.6 children; the per-adult scale is uncertain by about 5%). The per-consumer-unit r_k (the unit the task fixed) understates the per-adult value by that factor (about 1.9). tau_k is unit-free.
- **Standard errors.** District: whole-CZ bootstrap with the audited panel, screens and partition held fixed; the ACS sampling error of C24010 (margins of error) is reported separately (small). PUMS: Census replicate weights; the PUMA classes, OMB list and POWPUMA conversion are treated as known. The K-12 moments resample the district gap and the composition term together over CZs and add the survey and ACS variances: conservative.
- **Literature.** Boyd et al. (2005) and Reininger (2012) are about the distance between a teacher's hometown and her school, not residence and workplace; the IBO figure was not verified against the PDF; no national NCES table of teachers' residence relative to the district of employment was found in SASS / NTPS.

## Files

- `teacher_residence.json`: everything above, machine readable.
- `teacher_residence_targets.toml`: the K-12 counterparts of the model moments with their SEs and specification ranges (`[k12.*]`), tau_k and r_k (`[tau_k.*]`, `[r_k.*]`), the EDU district gaps (`[locale.*]`: proxies), their bootstrap correlation with the four adopted targets, the overlap variants, the composition terms and the commuting statistics.
- `teacher_residence_draws.csv`: the bootstrap draws (one row per draw) of the four adopted targets, the EDU moments and the K-12 moments (`k12_*`).
- Raw inputs under `raw/` (gitignored): `acs_c24010_school_districts_2018.parquet`, `acs_pums/acs1_2018_workers/<ST>.parquet`, `puma_geo/`.
