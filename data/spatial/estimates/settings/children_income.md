# Within-CZ parental-income composition of school-age children (task D5, DESIGN_v2.md)

Generated 2026-10-05T02:42:32 by `children_income.py`. NCES EDGE ACS-ED Children's Tabulation, ACS 5-year 2014-18 (Profile Tables CDP03.2, CDP05.1, CDP02.1, CDP03.4; group 102 'Relevant Children' is the primary, 104 and 105 are variants), merged to the audited district panel of `district_settings.py` and cut by its locale partition (city = NCES locale 11-13, suburb = 21-23; districts of at least 500 pupils; each location at least 10% of the CZ's partitioned pupils; 188 CZs, 3,260 districts, 28.55 million pupils). Statistics are within-CZ and are averaged over CZs with the partitioned pupils as weights; whole-CZ bootstrap, B = 1999, seed 20260929, the multiplicity matrix of `district_settings.py` (so the bootstrap correlations with the four district targets are valid). Hold-out half B (94 CZs) is in `sealed/children_income_holdout_B.json` and is not read, printed or summarised anywhere here; half A (94 CZs) is reported in section 7. Dollars are 2018 dollars.

## Summary

- **Within-CZ analogue of V6.** The share of a CZ's relevant children who live in the city is **55.5 (2.5)% in the bottom, 43.9 (2.6)% in the middle and 36.8 (2.5)% in the top national family-income tercile** (cut points $44,743 and $100,063). **Top minus bottom = -0.187 (SE 0.017)**, 95% interval [-0.220, -0.155], on 188 CZs. The NLSY79 V6 is 45.8 / 36.0 / 16.3% and -0.295 (0.037): the within-CZ gradient of 2014-18 is 63% as steep, 2.7 combined SEs flatter. Pooling the CZs (ratio of sums, the NLSY79 statistic's structure) gives 57.4 (2.4) / 44.5 (2.4) / 34.5 (2.8)% and -0.229 (0.020).
- **Levels differ, so the gradient should be read against the level.** The city holds 44.9 (2.5)% of the partition's children (NLSY79: 30.2% of SMSA children with known status), so the gradient relative to the overall city share is **-0.416 (0.045)** here against -0.977 for the NLSY79. The model's location 1 holds 0.407 of the population (experiments-2026-10b LOG.md, I0), close to the 42% of pupils and 45% of children here, not to the NLSY79's 30%.
- **The two locations.** Composition by national tercile (bottom / middle / top): city 40.9 (1.1) / 30.8 (0.4) / 28.3 (1.2)%, suburb 24.8 (0.9) / 32.3 (0.7) / 42.8 (1.5)%. Median family income of children (pooled brackets) $59,716 in the city, $89,025 in the suburb: **log gap 0.410 (0.042)** (mean family income 0.261 (0.037); average of district medians 0.398 (0.039)). Child poverty (related children 5-17) 25.4 (0.8)% against 13.0 (0.7)%, gap -0.124 (0.010).
- **Robustness.** Across the 11 variants of section 4 the gradient runs from -0.207 (`public`) to -0.164 (`city_vs_rest`); the largest moves from -0.187 are `city_vs_rest` (-0.164) and `public` (-0.207). Terciles at the bracket edges ($50,000 and $100,000) give -0.181; each CZ's own terciles -0.180; families as published -0.183. The Census B17024 cross-check (children 6-17 by income-to-poverty ratio, residence-based, a different income concept) gives -0.198 (0.018).
- **Correlation with the four district targets** (bootstrap): the gradient with `gap_income` -0.38, `gap_teachers_pp` -0.36, `gap_salary` -0.26, `gap_pupils` -0.06. The tercile shares themselves move with the size of the city (`city_share_t1` with `gap_pupils` -0.91); the difference does not.
- **Half A** (94 CZs): 52.7 (2.6) / 40.3 (3.0) / 34.2 (3.6)%, top minus bottom -0.185 (0.027) (full sample -0.187 (0.017)); log median gap 0.399 (0.062) (full 0.410 (0.042)). Half B is sealed.
- **For the calibration.** `city_share_gradient` -0.187 (0.017) is the within-CZ counterpart of the scored V6 (-0.295); it is a different object, not a re-estimate of the same number (section 1, caveats). For orientation the baseline model's V6 is +0.011 (calibration log, section 4), 12 SEs from it. The model's own counterpart is `parent_income_composition` (terciles of parental pre-tax income over all adults of one CZ); the `within_cz_cuts` variant is the closest data analogue of that single-CZ definition.

## 1. Data, definitions and units

- **Source.** NCES EDGE, ACS-ED Children's Tabulation, 'ACS 2014-2018 Profile' (school year 2017-18), Table Viewer <https://nces.ed.gov/programs/edge/TableViewer/acsProfile/2018>: 'All Districts' export, one table and one population group per file ('Download Text File'; pipe-delimited, one row per district record). 11,910 district records: 10,037 unified (GeoId prefix 97000US), 1,437 elementary (95000US), 436 secondary (96000US), keyed by the 7-digit NCES LEAID. Cached under `raw/children_income/` with a manifest (`manifest_acsed.json`: URLs, bytes, SHA-256) together with the ReadMe and the record-layout workbook (`ACS-ED_2014-2018_RecordLayouts.xlsx`) that ship in every export (the workbook's variable descriptions say 'in 2017 inflation-adjusted dollars', the Table Viewer and the ACS 2014-18 documentation say 2018 dollars; terciles and log gaps do not depend on it).
- **Population groups** ('iterations'): 102 *Relevant Children* (every child 0-17, or 18-19 without a high school diploma, who lives in the district and whose assigned grade is within the district's grade range, each child relevant to one district only, so overlapping elementary and secondary districts do not double count children); 104 *Relevant Children - Enrolled (public and private)*; 105 *Relevant Children - Enrolled Public*. **Primary: 102**, the population the NLSY79 statistic and the model's children correspond to (every child, whatever the school type); 105 is the 'children enrolled in public school' population of the task statement and is carried as a variant. The ACS-ED has no Children's detailed tables after the 2008-12 vintage, only the Profile Tables used here.
- **Variables** (record-layout names). CDP03.2 *Income and benefits*: `CDP03_34est` families; `CDP03_35est`-`CDP03_44est` families by **family income** in ten brackets (< $10k, $10-15k, $15-25k, $25-35k, $35-50k, $50-75k, $75-100k, $100-150k, $150-200k, $200k+); `CDP03_45est` / `CDP03_46est` median / mean family income. CDP05.1 *Sex and age*: `CDP05_1est` total population of the group = the number of children. CDP02.1: `CDP02_9est` average number of children in families with children (a check on the scaling). CDP03.4: `CDP03_67pct` related children 5-17 below the poverty level, percent. The household-income brackets (`CDP03_6`-`15`) are not used.
- **Unit and scaling.** The income table counts the *families in which the group's children live*; it is not a child-by-income tabulation. Each district's ten counts are scaled by its children per family (`CDP05_1est` / `CDP03_34est`; national mean 1.81; correlation with the published `CDP02_9est` 0.981 across the 8,020 districts of at least 500 families, 0.32 over all 11,909 where the rounding of tiny counts dominates), so that the district's counts add up to its children; this assumes constant family size across brackets inside a district. A family with children in two overlapping districts is in both records. The unscaled version is the `family` variant.
- **National terciles.** The national distribution of relevant children by family income is the sum of the 11,910 district records (all valid): 60.90 million children. The cut points at one third and two thirds are interpolated linearly inside the bracket that holds them (uniform density inside a bracket, the same split in every district): **$44,743** (inside $35-50k) and **$100,063** (inside $75-100k). Families as published: $45,290 and $101,627. The terciles are of the nation, as for V6 (NLSY79 cut points are over all 14-16-year-olds with income).
- **Statistics.** For CZ c, location l (1 = city, 2 = suburb) and tercile t, `T[c,l,t]` is the sum over the districts of the location of each bracket's count times the share of the bracket in t. `city_share_t` = T[c,1,t] / (T[c,1,t] + T[c,2,t]); `city_share_gradient` = city_share_3 - city_share_1 (the within-CZ analogue of V6); both averaged over CZs with the CZ weights (partitioned pupils). `pooled_*` is the ratio of sums over CZs. `comp_*` is the tercile composition of a location. Medians come from the pooled bracket counts of a location by linear interpolation (NaN if the median is in the open $200k+ bracket); the mean gap is the exact aggregate of the published district means (weights = children); child poverty is weighted by children.
- **Panel and sample.** `district_settings.build_audited_panel` (New York City one district; its collapsed row carries the Chancellor's id 3620580, under which the ACS-ED files it), locale partition, base year 2018: 3,260 districts. 3,192 have an ACS-ED record (99.64% of pupils; city 99.66%, suburb 99.63%): layers {'unified': 2527, 'elementary': 513, 'secondary': 152}. The four district targets reproduce (gap_pupils 0.3764 (0.1154), gap_salary 0.0086 (0.0092), gap_teachers_pp -0.0177 (0.0131), gap_income 0.1122 (0.0168)); the bootstrap draws of `gap_pupils` and `gap_income` rebuilt here equal those of `district_settings.py` (max abs difference 0e+00); the hold-out split is identical to the one stored in `district_settings.json` (seed 20261002).

## 2. City share by national family-income tercile

City share of the CZ's relevant children (group 102, children by scaling), CZ-weighted mean of the within-CZ ratios, whole-CZ bootstrap SE in parentheses; 'pooled' is the ratio of sums over CZs. NLSY79: central-city share of SMSA children 14-16 by tercile of 1978 family income (cut points $14,000 and $23,000, nominal).

| tercile | family income cut | share of the nation's children | city share (SE) | 95% interval | pooled over CZs (SE) | NLSY79 V6 |
|---|---|---|---|---|---|---|
| T1 (bottom) | below $44,743 | 33.3% | 55.5 (2.5)% | [51.0, 60.4] | 57.4 (2.4)% | 45.8% |
| T2 | $44,743 to $100,063 | 33.3% | 43.9 (2.6)% | [39.1, 49.0] | 44.5 (2.4)% | 36.0% |
| T3 (top) | $100,063 and above | 33.3% | 36.8 (2.5)% | [32.0, 41.6] | 34.5 (2.8)% | 16.3% |
| **top minus bottom** | | | **-0.187 (0.017)** | [-0.220, -0.155] | -0.229 (0.020) | -0.295 (0.037) |
| middle minus bottom | | | -0.116 (0.008) | [-0.132, -0.100] | | |
| top minus middle | | | -0.071 (0.011) | [-0.094, -0.050] | | |
| all children | | | 44.9 (2.5)% | [40.1, 49.6] | | 30.2% (1.5) |
| gradient / overall city share | | | -0.416 (0.045) | [-0.505, -0.334] | | -0.977 |

The within-CZ gradient is -0.187; the NLSY79 V6 is -0.295 (0.037); the difference is +0.108 (+2.7 combined SEs; the two samples are independent). n = 188 of 188 CZs (a CZ drops out only if a location has no district with an ACS-ED record).

### Composition of the two locations by national tercile

| tercile | city: share of its children (SE) | suburb (SE) | suburb minus city, pp |
|---|---|---|---|
| T1 | 40.9 (1.1)% | 24.8 (0.9)% | -16.1 |
| T2 | 30.8 (0.4)% | 32.3 (0.7)% | +1.6 |
| T3 | 28.3 (1.2)% | 42.8 (1.5)% | +14.5 |

## 3. Children by family income in the two locations

Shares of each location's children in the ten brackets (CZ-weighted means of the within-CZ shares), suburb minus city, and the city share of the children in each bracket (CZ-weighted mean of the within-CZ ratios).

| family income | city, % (SE) | suburb, % (SE) | suburb minus city, pp (SE) | city share of the bracket's children, % (SE) | nation, % |
|---|---|---|---|---|---|
| < $10k | 8.0 (0.4) | 3.9 (0.2) | -4.1 (0.4) | 59.6 (2.3) | 6.0 |
| $10-15k | 4.5 (0.2) | 2.3 (0.1) | -2.2 (0.2) | 59.3 (2.4) | 3.4 |
| $15-25k | 10.2 (0.3) | 5.6 (0.2) | -4.6 (0.3) | 57.7 (2.4) | 7.9 |
| $25-35k | 10.1 (0.2) | 6.6 (0.2) | -3.4 (0.3) | 53.9 (2.5) | 8.5 |
| $35-50k | 12.4 (0.2) | 9.8 (0.3) | -2.6 (0.3) | 50.0 (2.6) | 11.5 |
| $50-75k | 15.6 (0.2) | 15.5 (0.4) | -0.1 (0.3) | 45.2 (2.6) | 16.4 |
| $75-100k | 10.8 (0.2) | 13.4 (0.3) | 2.6 (0.2) | 40.5 (2.6) | 12.9 |
| $100-150k | 13.3 (0.3) | 20.0 (0.4) | 6.7 (0.5) | 36.9 (2.5) | 16.9 |
| $150-200k | 6.3 (0.3) | 10.3 (0.4) | 4.0 (0.4) | 35.5 (2.5) | 7.7 |
| $200k+ | 8.7 (0.7) | 12.6 (1.0) | 3.8 (0.9) | 38.0 (2.6) | 8.8 |

| statistic | city | suburb | suburb minus city (SE) | 95% interval |
|---|---|---|---|---|
| pooled median family income of children, $ | 59,716 | 89,025 | log gap **0.4100 (0.0424)** | [0.329, 0.488] |
| mean family income of children, $ | 88,758 | 114,118 | log gap 0.2612 (0.0373) | [0.190, 0.333] |
| average district log median family income (published medians) | | | 0.3979 (0.0390) | [0.322, 0.471] |
| child poverty rate, related children 5-17 | 25.4 (0.8)% | 13.0 (0.7)% | -0.1244 (0.0103) (levels) | [-0.144, -0.104] |

For scale: the audited earnings gap `gap_income` (log pooled median earnings of full-time year-round workers, suburb minus city) is 0.1122 (0.0168); the family income of *children* differs by 0.41 in the median and 0.26 in the mean because family income aggregates two earners, transfers and capital income, and children are concentrated in high-income suburban families.

## 4. Sensitivity

Full sample and half A, same national terciles unless stated. T1-T3: city share by tercile; the last columns are the log suburb-minus-city gap in median family income and the child poverty gap.

| variant | what | T1, % | T2, % | T3, % | top minus bottom (SE) | half A | median gap (SE) | poverty gap | n CZs |
|---|---|---|---|---|---|---|---|---|---|
| `primary` | primary: relevant children (102), children by scaling, national terciles, all layers | 55.5 | 43.9 | 36.8 | -0.187 (0.017) | -0.185 | 0.410 (0.042) | -0.124 | 188 |
| `public` | relevant children enrolled in public school (105) | 55.5 | 43.8 | 34.9 | -0.207 (0.018) | -0.201 | 0.456 (0.045) | -0.136 | 188 |
| `enrolled` | relevant children enrolled in public or private school (104) | 55.2 | 43.5 | 36.6 | -0.186 (0.017) | -0.184 | 0.408 (0.043) | -0.123 | 188 |
| `family` | families as published (no scaling to children), group 102 | 54.7 | 43.1 | 36.4 | -0.183 (0.017) | -0.181 | 0.406 (0.042) | -0.124 | 188 |
| `unified` | unified districts only | 57.1 | 45.3 | 38.3 | -0.189 (0.018) | -0.186 | 0.410 (0.043) | -0.126 | 183 |
| `nyc_no_record` | New York City without its ACS-ED record | 51.1 | 40.1 | 34.1 | -0.170 (0.017) | -0.185 | 0.382 (0.038) | -0.118 | 188 |
| `detroit_no_record` | Detroit without its ACS-ED record (no id crosswalk) | 55.1 | 43.7 | 36.8 | -0.183 (0.017) | -0.177 | 0.401 (0.042) | -0.122 | 188 |
| `drop_ny_cz` | New York CZ left out | 53.9 | 42.3 | 36.0 | -0.179 (0.016) | -0.185 | 0.390 (0.040) | -0.124 | 187 |
| `bracket_cuts` | terciles at the bracket edges closest to 1/3 and 2/3 (no interpolation) | 54.9 | 43.1 | 36.8 | -0.181 (0.017) | -0.179 | 0.410 (0.042) | -0.124 | 188 |
| `within_cz_cuts` | terciles of each CZ's own distribution (cuts in every CZ) | 54.9 | 42.9 | 36.9 | -0.180 (0.017) | -0.177 | 0.410 (0.042) | -0.124 | 188 |
| `equal_cz_weights` | every CZ weighted equally | 62.3 | 51.5 | 45.0 | -0.173 (0.012) | -0.189 | 0.350 (0.023) | -0.113 | 188 |
| `city_vs_rest` | city against every other district of the CZ (suburb, town and rural), CZs where the city holds at least 10% of the pupils | 48.3 | 37.3 | 31.8 | -0.164 (0.017) | -0.161 | 0.376 (0.044) | -0.116 | 181 |
| `b17024` | Census B17024, children 6-17 by income-to-poverty ratio, unified + elementary records, national IPR terciles (cut points 1.65 and 3.74) | 55.5 | 42.7 | 35.7 | -0.198 (0.018) | -0.195 | 0.446 (0.043) (median IPR) | -0.124 | 188 |

The `equal_cz_weights` variant weights every CZ equally (levels rise because the city share of children is larger in small CZs); `within_cz_cuts` uses each CZ's own terciles, so T1-T3 are within-CZ terciles and the gradient is free of between-CZ income differences. The gradient steepens with CZ size: -0.173 with equal weights, -0.187 with pupil weights, -0.229 pooled over CZs (the pooled ratio is dominated by the largest CZs).

## 5. Child poverty: three sources on the same districts

Suburb minus city gap in the poverty rate of children, districts with an ACS-ED record, a B17024 record and a SAIPE record (3,038 districts, 93.0% of pupils; New York City, whose SAIPE record is filed under another id, is not in it), CZ-weighted, same draws.

| source | city | suburb | gap (SE) |
|---|---|---|---|
| ACS-ED, related children 5-17 (CDP03_67pct) x children | 24.6% | 13.1% | -0.1159 (0.0100) |
| Census B17024, children 6-17 with IPR below 1.00 | 24.8% | 13.2% | -0.1161 (0.0098) |
| SAIPE, children 5-17 in poverty / children 5-17 | 21.0% | 12.1% | -0.0888 (0.0091) |

The two ACS sources agree; SAIPE, the estimator behind the scored `poverty_gap` (V3), shows a gap 23% smaller on the same districts (the ACS gap is 1.31 times the SAIPE gap). SAIPE district estimates are partly model-based; whether that explains the difference was not tested here. On the full ACS-ED sample the gap is -0.1244 (0.0103) (with New York City's record) and the audited SAIPE-based `poverty_gap` is -0.0888 (0.0090).

## 6. Bootstrap correlation with the four district targets

Correlation of the whole-CZ bootstrap draws of the new statistics with the draws of the four audited locale targets (the same multiplicity matrix, B = 1999, seed 20260929).

| | `gap_pupils` | `gap_salary` | `gap_teachers_pp` | `gap_income` |
|---|---|---|---|---|
| `city_share_gradient` | -0.06 | -0.26 | -0.36 | -0.38 |
| `city_share_t1` | -0.91 | +0.16 | +0.50 | -0.15 |
| `city_share_t2` | -0.96 | +0.08 | +0.48 | -0.27 |
| `city_share_t3` | -0.94 | -0.02 | +0.25 | -0.41 |
| `gap_median_family_income` | +0.23 | +0.13 | +0.40 | +0.47 |
| `gap_mean_family_income` | +0.19 | +0.15 | +0.36 | +0.59 |
| `gap_child_poverty` | -0.32 | -0.12 | -0.25 | -0.45 |
| `comp_city_t1` | +0.27 | +0.13 | +0.15 | +0.64 |
| `comp_suburb_t3` | +0.15 | +0.02 | +0.39 | -0.06 |

Targets (value, SE): `gap_pupils` 0.3764 (0.1154); `gap_salary` 0.0086 (0.0092); `gap_teachers_pp` -0.0177 (0.0131); `gap_income` 0.1122 (0.0168).

## 7. Hold-out half A

Pre-registered split of the locale CZs (`district_settings.holdout_split`: Sort by partitioned pupils (descending, ties by ascending CZ id), pair consecutive CZs (1-2, 3-4, ...; an unpaired last CZ goes to A); for each pair one `rng.integers(2)` from `default_rng(seed)`, in pair order: 0 sends the first CZ of the pair to hold-out B, 1 sends the second; the other goes to A.; seed 20261002); half A has 94 CZs, half B 94. Bootstrap on half A's own multiplicity matrix (n = 94, B = 1999, seed 20260929). Half B was computed by the same code and written only to `sealed/children_income_holdout_B.json`.

| statistic | full sample (SE) | half A (SE) |
|---|---|---|
| city share, tercile 1, % | 55.5 (2.5) | 52.7 (2.6) |
| tercile 2, % | 43.9 (2.6) | 40.3 (3.0) |
| tercile 3, % | 36.8 (2.5) | 34.2 (3.6) |
| all children, % | 44.9 (2.5) | 42.0 (3.0) |
| **top minus bottom** | -0.187 (0.017) | -0.185 (0.027) |
| top minus bottom, pooled over CZs | -0.229 (0.020) | -0.239 (0.032) |
| gradient / overall city share | -0.416 (0.045) | -0.442 (0.083) |
| log gap, median family income of children | 0.410 (0.042) | 0.399 (0.062) |
| log gap, mean family income | 0.261 (0.037) | 0.242 (0.056) |
| child poverty gap | -0.124 (0.010) | -0.130 (0.017) |

Half-A district targets (value, SE): `gap_pupils` 0.5220 (0.1515); `gap_salary` 0.0090 (0.0138); `gap_teachers_pp` -0.0312 (0.0149); `gap_income` 0.1183 (0.0234).

## 8. Checks

- **National row, group 104.** The EDGE preview's national table has 31,882,600 families; the sum of the district records has 31,955,899 (+0.23%: families with children in two overlapping districts are in both records); bracket shares differ by at most 0.23 pp; national median / mean family income $71,399 / $98,062.
- **National row, group 105.** The EDGE preview's national table has 28,351,170 families; the sum of the district records has 28,598,479 (+0.87%: families with children in two overlapping districts are in both records); bracket shares differ by at most 0.21 pp; national median / mean family income $67,594 / $91,281.
- **National row, group 102.** Not defined by the EDGE viewer: for the nation 'relevant' children are all children.
- **Interpolated against published medians** (11,717 districts): mean log error +0.0067 (children-weighted +0.0098), SD 0.049, 1st / 99th percentile -0.129 / +0.143; 114 districts have no published median (open top interval or suppressed).
- **Rounding.** 93.0% of bracket counts are multiples of 5, none equals 1-3, 7.0% equal 4 (a count of 1-4 is evidently published as 4: 0.10% of national family mass; 0.007% in the partition).
- **Unmatched partition districts.** 68 (102,148 pupils, 0.36%); largest: Achievement School District (11,031); NORTH COLONIE CSD (5,844); Central Community School District (4,714); Essex-Westford Educational Community Unified School District (3,980).
- **Children against pupils.** Matched city districts: 16.32 million relevant children, 12.03 million CCD pupils; suburb 20.11 / 16.52 million.

## 9. Caveats and what could not be done

- **The income table counts families, not children.** The ACS-ED has no child-by-family-income cross-tabulation: CDP03.2 gives the families (households) in which the group's children live, by family income. Children are obtained by scaling every district's ten counts by its children per family (CDP05_1est / CDP03_34est), which reproduces the number of children in each district exactly and the distribution under the assumption that the number of children per family does not vary with family income inside a district (it does vary across districts, which the scaling captures). A family with children in two overlapping districts (an elementary and a secondary district) is counted in both records, a child in one; the unscaled variant (families as published) gives -0.183 against -0.187.
- **Relevant children, not pupils.** Group 102 is every child who lives in the district and whose grade it serves, whatever the school type; it is the population the NLSY79 and the model's children correspond to. Group 105 (enrolled in public school) drops private-school and unenrolled children, who are richer and more urban on average: the gradient is -0.207 (SE 0.018). Group 102 is the primary for that reason; the public-school and enrolled groups are reported beside it. The CZ weights are the partitioned public pupils of the CCD, as for every district target, while the children counted here are residents: 16.3 million relevant children in the matched city districts against 12.0 million pupils (charter and other independent LEAs are not in the pupil count; their children are counted in their district of residence).
- **Not like-for-like with the NLSY79 statistic.** V6 is the 1979 central-city share of SMSA children aged 14-16 with known central-city status (20-40% of SMSA residents have none), by tercile of 1978 family income of all 14-16-year-olds, pooled across SMSAs. Here: 2014-18, city = the NCES locale of the child's district, terciles of the national distribution of children by family income in 2018 dollars, averaged over CZs. The levels differ for that reason (the city holds 44.9% of the partition's children, 30.2% of the NLSY79's MSA children); the comparison of gradients is suggestive, not an estimate of a change over time. The children also differ in age: relevant children are of school age (about 5-17, and 18-19 without a diploma; pre-school children only where a district serves them), the NLSY79's are 14-16.
- **City and suburb are district-level.** A district is in one location, whatever its internal geography: large central-city districts include suburban neighborhoods and the reverse, so the city share by income is a district-composition statistic. Only districts of at least 500 pupils in CZs with both locations are used (the effective number of CZs is 42, the five largest carry 26% of the weight, New York 5.3% alone). Town and rural districts (NCES locale 31-43) are outside the partition, whereas the NLSY79's 'suburb' is everything in the SMSA outside the central city; the `city_vs_rest` variant puts them in the second location (against the first, the city) and gives -0.164 on 181 CZs (the city then holds 38.5% of the children).
- **Matching.** 3,192 of the 3,260 partition districts (99.64% of pupils) have an ACS-ED record. The 68 that do not (102 thousand pupils) are mostly state-wide, county-office, regional vocational, virtual and charter LEAs with no residential boundary, plus a few the id match misses (e.g. North Colonie CSD, NY, 5,844 pupils: its ACS-ED record appears to carry the id 3600971 against the CCD's 3621000; not crosswalked). Detroit's record is filed under the pre-2016 id 2612000; it is matched to 2601103 (Detroit Public Schools Community District, 49,931 pupils) by an explicit crosswalk: without it the gradient is -0.1829 (see section 4). New York City's record (filed under the Chancellor's id 3620580, which the audited panel carries for its collapsed row) is used, unlike the earnings, poverty and SEDA targets, whose New York CZ has no city record; without it the gradient is -0.170, the largest of the variants.
- **Interpolation.** The national cut points fall inside the $35-50k bracket (first) and at the edge of the $75-100k bracket (second, $100,063 for group 102), so the first tercile boundary assumes a uniform density inside $35-50k; the variant with terciles at the bracket edges closest to 1/3 and 2/3 ($50k and $100k, national group shares 37.4, 29.3, 33.4%) gives -0.181. Medians interpolated from a district's own ten brackets exceed its published median by 0.7% on average (SD 4.9%), which shifts the level of the pooled-median gap little because both locations are interpolated.
- **Rounding.** ACS-ED counts are rounded to the nearest 5, and the values 1-3 never occur while 4 occurs in 7.0% of cells (so a count of 1-4 is evidently published as 4); these cells hold 0.10% of all families nationally and 0.007% in the partition: negligible for pooled counts.
- **Standard errors.** The whole-CZ bootstrap resamples CZs, so it carries the sampling error of the ACS estimates (it enters the CZ-level statistics) and the dispersion of the statistic across CZs, but not the uncertainty of the national cut points (a fraction of a percent of income) or of the scaling assumption.
- **Child poverty: two ACS sources agree, SAIPE shows a smaller gap.** On the districts that have all three, the suburb-minus-city gap in the child poverty rate is -0.116 (ACS-ED), -0.116 (B17024) and -0.089 (SAIPE, the estimator behind the scored `poverty_gap` = V3): the ACS gap is 1.31 times the SAIPE gap. SAIPE district estimates are partly model-based (they may smooth differences between districts of one CZ; not tested here). If the direct ACS estimates are the better measure of the within-CZ gap, V3 understates the sorting by that factor; V3 is not changed here.
- **Not done.** (1) No child-weighted income distribution: the ACS-ED cross-tabulation does not exist (children are scaled, see above), and a PUMS version (children by exact family income, PUMAs classified as principal-city or suburban as in `teacher_residence.py`) was not attempted: PUMAs do not nest in school districts, so it would be a different city/suburb definition. (2) Income detail is limited to ten brackets (the Census detailed table B19101 has 16, but for all families, not families with children). (3) The national row of the EDGE preview is not defined for group 102 (for the nation 'relevant' children are all children), so the national cross-check of the district sum covers groups 104 and 105 only. (4) Census B17024 is a cross-check in another income concept (income relative to the poverty line, family size adjusted), for children 6-17 by residence (unified and elementary records, each child once), not a substitute. (5) Nothing here is a model-implied moment: the model's moment is the share of the young schooled in location 1 by tercile of the parent's pre-tax income across all adults (`parent_income_composition`), the data counterpart of which is `city_share_t1..3` with a single CZ's terciles replaced by national ones (the `within_cz_cuts` variant uses each CZ's own).
- **Seal.** The CZ ids and CZ weights of half B are public (`district_settings.json`), so its point estimates are in principle implied by the full-sample and half-A numbers; the seal only keeps them out of the logs and reports. Nothing in this script's log, this report, the JSON or the TOML reads or prints half B.

## Files

- `children_income.py` (fetch / build / selftest); `children_income.json` (all statistics, full sample and half A); `children_income_targets.toml` (full-sample values and SEs); `children_income_draws.csv` and `children_income_draws_holdout_A.csv` (bootstrap draws with the four district targets); `sealed/children_income_holdout_B.json` (half B, sealed).
- Raw: `raw/children_income/` (ACS-ED exports for groups 102, 104, 105, the ReadMe and record layouts, `manifest_acsed.json`, B17024 parquet).
