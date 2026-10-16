# B1 hold-out (BA-or-more share gap, suburb minus city): methods (no values)

Generated 2026-10-15T08:49:22 by `ba_share_holdout.py`. **This file contains definitions, samples, coverage and PASS/FAIL checks only; it reports no value of B1 or of any variant.** The statistics are in `estimates/settings/sealed/ba_share_holdout_B1.json` (git-ignored), written with `json.dump` and neither printed nor logged. It is to be opened by a person only after the model specification of the 2026-10-15 round is chosen (`julia/spatial_model/calibration/runs/experiments-2026-10d/PREREG.md`, A7); no specification is re-selected on it.

## 1. What is sealed

**B1**, pre-registered under "Sealed holdout" in the 2026-10-15 PREREG: the within-CZ gap in the share of adults 25+ with a bachelor's degree or more, suburb minus city, ACS 5-year 2014-18 by school district, on the audited locale partition and CZ weights of `district_settings.py` (the CZs of the four district targets). The model counterpart is the BA share of adults by residence (`ba_share_by_residence` / `ba_share_gap`, M3 only). The sealed file holds, for the full sample and for the two pre-registered hold-out halves, the level gap (fraction and percentage points), the log ratio, the CZ-weighted city and suburb shares, the version pooled over CZs, whole-CZ bootstrap SEs and 95% intervals, the bootstrap correlations with the four district targets (`gap_income`, `gap_pupils`, `gap_salary`, `gap_teachers_pp`), the bootstrap draws of the headline, a per-CZ table, and the same statistics for 9 variants (240 sample x statistic cells in all).

## 2. What `ba_holders` measures (and why the data were pulled again)

The PREREG names `raw/acs_school_districts.parquet` (`ba_holders`, `pop_25plus`) as the input. `fetch_data.ACS_SD_VARS` maps `B15003_022E` to `ba_holders` and `B15003_001E` to `pop_25plus`. B15003 is *Educational attainment for the population 25 years and over*: one line per **highest** level completed. The Census metadata labels (read from the API at fetch time and stored in `raw/acs_b15003_labels_2018.json`) are: line 1 Total; line 22 **Bachelor's degree**; line 23 Master's degree; line 24 Professional school degree; line 25 Doctorate degree. `ba_holders` is therefore the number of adults whose highest degree is *exactly* a bachelor's: it **excludes** holders of graduate and professional degrees, so `ba_holders / pop_25plus` is *not* the "BA-or-more" share the PREREG defines. The headline therefore uses lines 22-25 of the same table, pulled for the same vintage (ACS 2018 5-year) and the same three school-district layers (`python3 ba_share_holdout.py fetch`, Census API, state by state). The exactly-bachelor's share is the variant `ba_only`.

Check of the new pull against the old file (13,389 records each): 13,389 of 13,389 records have `B15003_001E` = `pop_25plus` and 13,389 have `B15003_022E` = `ba_holders` (records only in one file: 0 / 0); lines 2-25 add up to line 1 in 13,389 records; 0 records carry an annotation code. **Departure from the PREREG wording:** the data file named there does not measure the concept it names; the concept (BA or more) is kept and the data are replaced by the full table.

## 3. Definitions

**Records.** One row per `leaid` (state FIPS + 5-digit district code) and layer (unified, elementary, secondary). The "Remainder of <state>" rows (code 99999; 75 rows) are not districts and are dropped (as `district_bootstrap.load_income_columns`); 13,314 records remain, none with a repeated id. `adults` = `B15003_001E`; `BA+` = `B15003_022E + 023E + 024E + 025E`; a record is valid if the adult count and the four degree lines are published (no annotation code), there is at least one adult, and BA+ does not exceed the adults. No record carries an annotation code and lines 2-25 add up to line 1 in every record; the 11 records that are not used are districts with no adults 25+ ("School District Not Defined", unorganized territories, a few elementary and secondary districts without residents).

**Panel and partition.** `district_settings.build_audited_panel` (New York City collapsed to one district under the Chancellor's id 3620580, under which the ACS files it) and `scheme_tables(A, "locale")`: city = NCES locale 11-13, suburb = 21-23, districts with at least 500 pupils, each location holding at least 10% of the CZ's partitioned pupils: 188 CZs, 3,260 districts (625 city, 2,635 suburb), 28.55 million pupils. The B15003 records are merged on `leaid` after the NYC collapse. CZ weights are the partitioned pupils, as for the four targets.

**Location shares and the CZ gap.** For CZ c and location l (1 = city, 2 = suburb), s_{c,l} = (sum of BA+ over the location's districts with a valid record) / (sum of adults 25+ over the same districts): a ratio of location totals, as every per-pupil object of `data_estimate.Gap` (not a mean of district shares). The CZ gap is s_{c,2} - s_{c,1}; a CZ with no valid record in one location has no gap and drops out of every statistic. **B1 = the CZ-weighted mean of the CZ gaps** (weights: partitioned pupils), a fraction; the sealed file also has it in percentage points, the CZ-weighted mean of log(s_{c,2} / s_{c,1}), the CZ-weighted means of the city and suburb shares over the CZs that have the gap (so that their difference is the gap), and the version pooled over CZs (sum of BA+ over sum of adults, per location, adult-weighted, then suburb minus city).

**Overlapping layers.** An adult who lives in an elementary district inside a secondary district appears in both records, so summing every layer double-counts adults (the unified and elementary records tile each state; the secondary layer overlays the elementary one). The headline counts adults **once, from the unified and elementary records** (the rule of `teacher_residence.py` for employed residents and of `children_income.py` for children); the 152 secondary LEAs of the partition (3.15% of its pupils) keep their pupils in the CZ weights and contribute no adults. This rule was fixed on those precedents before any share was computed. The rule of the four earnings targets (`load_income_columns`: every LEA takes the record of its own layer) is the variant `all_layers`; the two other layer rules (`unified_secondary`, `unified_only`) are in section 4.

**Matching.** ACS record to LEA on the 7-digit `leaid`. Detroit's record is filed under the pre-2016 id 2612000 (Detroit City School District) while the CCD frame carries its successor Detroit Public Schools Community District (2601103; the same boundary): it is refiled to 2601103, as in `children_income.py` (variant `no_detroit_crosswalk` leaves it out). No other crosswalk is applied: the other unmatched LEAs are listed in section 6 (non-geographic LEAs such as county offices, regional vocational districts and virtual schools have no ACS record).

## 4. Specifications (all sealed)

| key | description |
|---|---|
| `headline` | HEADLINE: BA or more (B15003 lines 22-25), adults counted once (unified + elementary records), Detroit refiled to 2601103, CZ weights = partitioned pupils |
| `all_layers` | every LEA takes the record of its own layer (the rule of the earnings targets): adults in an elementary district inside a secondary district are counted twice when both sit in one location |
| `unified_secondary` | unified + secondary records: the elementary LEAs' adults are dropped in the CZ locations that hold a valid secondary record |
| `unified_only` | unified districts only: elementary and secondary LEAs leave the adult sums |
| `ba_only` | bachelor's degree as the HIGHEST degree (line 22 = the old `ba_holders`, over `pop_25plus`), adults counted once |
| `nyc_dropped` | New York City's adults leave the sums (its pupils stay in the CZ weights) |
| `drop_ny_cz` | the New York CZ leaves the CZ mean |
| `no_detroit_crosswalk` | Detroit's record is not refiled (Detroit has no adults) |
| `adult_weights` | CZ weights = adults 25+ with a valid record in the CZ's partitioned districts |
| `equal_weights` | every CZ weighted equally |

Every specification is evaluated on the full sample and on each hold-out half. The headline is `headline`; the variants are reported for robustness and are not used to choose among specifications.

## 5. Inference and hold-out

Whole-CZ bootstrap, B = 1999, seed 20260929, the multiplicity matrix of `district_settings.py` (`district_bootstrap.draw_counts(188, B, seed)`): SE = standard deviation (ddof 1) of the draws, 95% interval = 2.5 and 97.5 percentiles, CZs with an undefined statistic drop out of the draw. The draws are therefore the same draws as those of the four district targets, and the bootstrap correlations with them are valid. The pre-registered split of the locale CZs (`district_settings.holdout_split`, seed 20261002; pairs of CZs sorted by pupils, one to each half) is reproduced here (94 CZs in A, 94 in B) and compared with the assignment stored in `district_settings.json`; each half is bootstrapped on its own number of CZs with the same seed, as `district_settings.py` does. Only the ACS sampling variance of the district counts is not added (section 9).

## 6. Sample and coverage

Partition: 188 CZs, 3,260 districts, 28.55 million partitioned pupils (CZ weights: effective number of CZs 41.7, five largest CZs 25.8% of the weight). ACS record layers of the partitioned districts: elementary 518 (6.3% of pupils), none 55 (0.3% of pupils), secondary 152 (3.1% of pupils), unified 2,535 (90.2% of pupils).

| specification | districts with a valid record | share of districts | share of partition pupils | city / suburb share of pupils | CZs with a gap (of 188) | share of the CZ weight | in half A / B |
|---|---|---|---|---|---|---|---|
| `headline` | 3,053 | 93.7% | 96.5% | 96.9% / 96.2% | 188 | 100.0% | 94 / 94 |
| `all_layers` | 3,205 | 98.3% | 99.7% | 99.7% / 99.7% | 188 | 100.0% | 94 / 94 |
| `unified_secondary` | 2,710 | 83.1% | 93.6% | 94.3% / 93.1% | 188 | 100.0% | 94 / 94 |
| `unified_only` | 2,535 | 77.8% | 90.2% | 91.5% / 89.3% | 183 | 98.7% | 93 / 90 |
| `ba_only` | 3,053 | 93.7% | 96.5% | 96.9% / 96.2% | 188 | 100.0% | 94 / 94 |
| `nyc_dropped` | 3,052 | 93.6% | 93.2% | 88.9% / 96.2% | 188 | 100.0% | 94 / 94 |
| `drop_ny_cz` | 3,053 | 93.7% | 96.5% | 96.9% / 96.2% | 187 | 94.7% | 94 / 93 |
| `no_detroit_crosswalk` | 3,052 | 93.6% | 96.4% | 96.5% / 96.2% | 188 | 100.0% | 94 / 94 |
| `adult_weights` | 3,053 | 93.7% | 96.5% | 96.9% / 96.2% | 188 | 100.0% | 94 / 94 |
| `equal_weights` | 3,053 | 93.7% | 96.5% | 96.9% / 96.2% | 188 | 100.0% | 94 / 94 |

Districts of the partition outside the headline sample: 207 (3.47% of partition pupils). By reason: no ACS record 55 (91,550 pupils); secondary record (adults counted from the elementary records) 152 (899,308 pupils). Largest by pupils (public CCD names):

| LEA | leaid | pupils | locale | reason |
|---|---|---|---|---|
| Kern High | 0619540 | 40,340 | city | secondary record (adults counted from the elementary records) |
| Sweetwater Union High | 0638640 | 39,165 | suburb | secondary record (adults counted from the elementary records) |
| Anaheim Union High | 0602630 | 30,292 | city | secondary record (adults counted from the elementary records) |
| Phoenix Union High School District (4286) | 0406330 | 27,573 | city | secondary record (adults counted from the elementary records) |
| Chaffey Joint Union High | 0608160 | 23,883 | suburb | secondary record (adults counted from the elementary records) |
| East Side Union High | 0611820 | 22,606 | city | secondary record (adults counted from the elementary records) |
| William S. Hart Union High | 0642510 | 22,416 | suburb | secondary record (adults counted from the elementary records) |
| Antelope Valley Union High | 0602820 | 21,642 | suburb | secondary record (adults counted from the elementary records) |
| Grossmont Union High | 0616230 | 16,760 | suburb | secondary record (adults counted from the elementary records) |
| Oxnard Union High | 0629270 | 16,701 | city | secondary record (adults counted from the elementary records) |

Largest districts with no ACS record (non-geographic LEAs -- county offices, regional vocational and virtual districts, statewide turnaround districts -- and districts whose CCD and Census ids differ; none is refiled except Detroit):

| LEA | leaid | pupils | locale |
|---|---|---|---|
| Achievement School District | 4700147 | 11,031 | city |
| NORTH COLONIE CSD | 3621000 | 5,844 | suburb |
| Central Community School District | 2200125 | 4,714 | suburb |
| Essex-Westford Educational Community Unified School District | 5000395 | 3,980 | suburb |
| Pinkerton Academy School District | 3307360 | 3,230 | suburb |
| Orange County Department of Education | 0691024 | 3,093 | city |
| Kern County Office of Education | 0691012 | 2,508 | city |
| RSU 22 | 2314815 | 2,429 | suburb |
| Mastery Schools of Camden Inc. | 3400772 | 2,317 | city |
| Greater Lowell Regional Vocational Technical | 2505480 | 2,256 | suburb |

Coverage by CZ location (share of the location's partitioned pupils whose district carries a valid headline record; 376 CZ locations): minimum 45.0%, fifth percentile 89.6%, median 100.0%; 1 below 50% and 13 below 80% (their CZs carry 8.7% of the CZ weight). In those locations 97.0% of the unmatched pupils are in secondary LEAs (whose pupils count as unmatched although the same households' adults are counted through the elementary districts) and the rest in LEAs without a record; the per-CZ coverage is in the sealed `cz_table` (`pupil_cov_city`, `pupil_cov_suburb`).

Overlap: the adults of all valid records of the partition exceed those counted once by the factor 1.0869 (city 1.0665, suburb 1.1046); the headline counts 130,960,933 adults, the `all_layers` variant 142,336,363.

Tiling of the states by the unified and elementary records (adults 25+, `B15003_001E`): national sum of unified + elementary records over the state totals of the same line = 0.9997; with the secondary layer added 1.0661; 51 of 52 states within 1% and 52 within 5% (smallest 0.99, largest 1.00). States with the largest double-count when all layers are summed: Montana 1.87; Illinois 1.38; Arizona 1.29; California 1.27; New Jersey 1.19.

New York City (`3620580`): in the partition (location 1, 960,484 pupils, 63.0% of its CZ's partitioned pupils, its CZ carries 5.34% of the CZ weight); ACS record layer unified; in the headline sample: True.
Detroit (`2601103`): in the partition (location 1, 49,931 pupils, 8.3% of its CZ's partitioned pupils, its CZ carries 2.12% of the CZ weight); ACS record layer unified; in the headline sample: True.

## 7. Checks run before sealing (no sealed value involved)

`seal` runs every check first and refuses to seal if one fails; `check` prints the same lines.

```
  Census labels of lines 1 and 22-25: PASS  [total; bachelor's; master's; professional school; doctorate]
  new pull = old file (same records, pop_25plus, ba_holders): PASS  [13,389 of 13,389 records]
  B15003 lines 2-25 add up to line 1; no annotation codes: PASS
  one record per leaid after dropping the 'Remainder' rows: PASS  [75 remainder rows dropped, 13,314 records]
  every record that is not used is an empty district (0 adults 25+): PASS  [11 records with no adults, 11 not used]
  the four locale targets reproduce the reported values: PASS  [value and SE]
  hold-out split = district_settings.json: PASS  [94 / 94 CZs, seed 20261002]
  bootstrap draws = district_settings.py (recomputed; stored csv): PASS  [max abs diff 0.0e+00; csv 5.0e-11]
  statistic machinery reproduces the income and enrollment targets: PASS  [max abs diff 0.0e+00]
  B1 tables: the 188 CZs and weights of the audited table, every specification: PASS  [188 CZs x 10 specifications]
  headline: valid gap in >= 90% of the CZs and of the CZ weight: PASS  [188 of 188 CZs, 100.0% of the weight]
  New York City and Detroit have a record in the headline sample: PASS
```

Reproduction of the four district targets on the same draws (value (bootstrap SE); the digits of `district_audit.md` / `district_bootstrap.md`):

| target | computed here | reported | n CZs |
|---|---|---|---|
| `gap_pupils` | 0.3764 (0.1154) | 0.3764 (0.1154) | 188 |
| `gap_salary` | 0.0086 (0.0092) | 0.0086 (0.0092) | 188 |
| `gap_teachers_pp` | -0.0177 (0.0131) | -0.0177 (0.0131) | 188 |
| `gap_income` | 0.1122 (0.0168) | 0.1122 (0.0168) | 188 |

The 1999 x 188 multiplicity matrix reproduces the draws of the four targets: maximum absolute difference between the draws recomputed from the CZ table and those of `moment_set` 0.0e+00, against the stored `district_settings_draws_locale.csv` (1999 draws) 5.0e-11 (the file is rounded to 10 decimals). The two code paths that B1 uses are validated on public quantities: the weighted-mean/bootstrap routine reproduces the income target (difference 0.0e+00) and the CZ-table routine the enrollment gap column (0.0e+00). `ba_share_holdout_xcheck.py` recomputes every sealed statistic in separately written code (including the `ba_only` variant from the old parquet) and prints only PASS/FAIL and the maximum absolute differences.

## 8. The sealed file

`estimates/settings/sealed/ba_share_holdout_B1.json`: `meta` (B, seeds, code and input SHA-256), `definition`, `specs`, `headline` (`ba_share_gap` and `pp` for `full`, `A`, `B`: value, se, ci95_lo, ci95_hi, n_cz), `samples` (`full`, `A`, `B`: `n_cz`, `cz_ids`, `stats` with keys `ba_share_gap`, `ba_share_gap_pp`, `ba_log_ratio`, `ba_share_city`, `ba_share_suburb`, `pooled_gap`, `pooled_city`, `pooled_suburb`, and `corr`: the bootstrap correlation matrix of the draws of `gap_pupils`, `gap_salary`, `gap_teachers_pp`, `gap_income`, `ba_share_gap`, `ba_log_ratio`, `pooled_gap`), `variants` (the same `stats` per sample and specification, plus the correlation of the gap with the four targets), `draws_full` (the B draws of the headline gap, log ratio and pooled gap) and `cz_table` (per CZ: weight, location shares, BA+ and adult totals, gap, pupil coverage of the two locations). `python3 ba_share_holdout.py verify-seal` prints the structure and consistency checks (finite, shares in [0, 1], pp = 100 x fraction, city + gap = suburb, halves partition the CZs), never a value.

## 9. Caveats

- The ACS 5-year school-district tabulation is a sample (2014-18, about 5% of households over five years): the district counts have sampling error (margins of error are not pulled). The SE reported is the whole-CZ bootstrap SE, as for the four district targets; it does not add the ACS sampling variance of the district counts, which is small relative to the variation across CZs but not zero for small districts.
- The location of an adult is the locale of the school district he lives in, not of his household's neighbourhood; adults without children are included (B15003 is the whole population 25+, including group quarters). B1 is the education sorting of *residents*, not of parents.
- The headline counts adults from the unified and elementary records. A CZ location made up only of secondary-layer LEAs has no adults under the headline (the counts of such CZs are in section 6; `unified_secondary` and `all_layers` bracket the rule). Dual-system states (elementary + secondary districts) are where the overlap matters.
- The partition is on 2017-18 CCD districts of at least 500 pupils; the ACS boundaries are those of 2017. Non-geographic LEAs and LEAs whose id differs between the CCD and the ACS (apart from Detroit) have no record and leave the sample (section 6).
- CZ weights are the partitioned pupils (the four targets' weights), not adults. The `adult_weights` and `equal_weights` variants show the sensitivity.
- `ba_holders` of the existing parquet is the exactly-bachelor's count (section 2). `fetch_data.py` writes `ba_holders` and `pop_25plus`; no other script of the repository reads them (checked by `grep`), so no earlier estimate used them as a college share.
- This file reports no value; the checks above use only public quantities (the four reported district targets, counts, and Census totals of adults).

## 10. Disclosure log

What was looked at before sealing, and what the runs printed. (a) Structure inspection of `raw/acs_school_districts.parquet`, run twice: column names, dtypes, layer counts and its first three rows (Choctaw County, Claiborne County and Clarksdale Municipal, Mississippi: median income, home value, population, `ba_holders` and `pop_25plus`, i.e. the exactly-bachelor's counts of three small districts); no aggregate of any BA share was computed from it by location, CZ or sample. (b) The Census metadata labels of B15003 lines 1 and 22-25 (section 2). (c) `fetch` printed record counts only. (d) The comparison of the new pull with the old file printed counts of equal records. (e) The partition was matched to the ACS records and the unmatched LEAs listed by pupils (names and pupils: public inputs, section 6); the adult count (`B15003_001E`, a total, not a degree count) of the 2612000 record was printed to confirm that it is Detroit's before the refile was fixed; no degree count entered that decision. (f) `check` and `seal` print sample sizes, coverage shares of pupils, the reproduction of the four public targets and PASS/FAIL lines; nothing computed from a degree count is printed, logged or written outside the sealed file. (g) The sealed file was verified by `verify-seal` and the independent `ba_share_holdout_xcheck.py` (structure, PASS/FAIL, maximum absolute differences). (h) `seal` scans this file, and `leak-scan` the run logs, for renderings of the sealed numbers: the high-precision ones (fractions at 5-6 decimals, percentage points at 4) must not appear, and the coarse ones (fractions at 4 decimals, percentage points at 2-3, which collide with public numbers by chance) must not appear more often than for random decoys of the sealed numbers; the counts are in the run record. This file is generated from `prepare` output, never from a statistic.

SHA-256 of the sealed file at the time of writing: `de2fcdba2a0dcd57910b7417c774bda4415ca099b85d18760f5074f6a1eea94f`.
