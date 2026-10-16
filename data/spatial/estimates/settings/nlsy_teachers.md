# NLSY79: latent ability (AFQT) of K-12 teachers and non-teachers at ages 25-34

Generated 2026-10-03T00:41:14 by `nlsy_teachers.py`. NLSY79 population sample (block A of `nlsy_ability.py`, N = 9033), person-rounds at ages 25-34 (survey years 1982-1998). Ability = `afqt_c` of `nlsy_ability.py` (standardized, age-normed composite of the four ASVAB IRT z-scores, 1979 weights); a gap is a difference of group means of `afqt_c` divided by sqrt(R_comp_pop), i.e. in SDs of the latent factor (R_comp_pop = 0.916). Weights: 1979 `SAMPWEIGHT`, one observation per person. Household-cluster bootstrap, B = 500, seed 20260929; SE = bootstrap SD, in parentheses; 95% percentile CI in brackets; AFQT standardization and R_comp_pop are re-estimated in every draw. Comparison A: all non-teachers; B: non-teachers with 16+ years of schooling (against all teachers); C: teachers with 16+ years against non-teachers with 16+ years.

## Summary

- **Relative to all non-teachers, K-12 teachers are positively selected on latent AFQT: 0.614 (0.062) SD** (employed; teacher = K-12 teacher in the majority of 25-34 rounds; 160 teachers, 8159 non-teachers). Men 0.641 (0.167) (24 teachers), women 0.619 (0.067) (136 teachers). All respondents with any 25-34 occupation: 0.583 (0.065) (165 teachers). The gap runs through schooling: 92% of these teachers and 24% of employed non-teachers have 16+ years.
- **Relative to non-teacher college graduates, teachers are negatively selected: -0.334 (0.063) SD** (B; employed; 1608 comparison persons); men -0.351 (0.167), women -0.278 (0.069). Any 25-34 occupation: -0.376 (0.065). Within college graduates only (C): -0.238 (0.062) (136 teachers with 16+ years against 1608).
- In the model's log-z units a gap of g latent SDs is g x s_z log points, s_z = 0.1955 (`implied_s_z` of `nlsy_ability.py`; its sampling error is not propagated): comparison A +0.614 SD = +0.120; B -0.334 SD = -0.065; C -0.238 SD = -0.047.
- **Ever a K-12 teacher at 25-34** (employed sample, 449 persons, 103 men and 346 women): A 0.453 (0.046), B -0.531 (0.048), C -0.255 (0.048).
- **Teaching shares at 25-34** (employed K-12 teachers over all person-rounds with highest grade >= 9, not armed forces, home production included): men 0.98 (0.15)%, women 3.45 (0.29)%, against ACS 2009-13 1.94% and 5.99% (project headline, unweighted persons; PWGTP-weighted 1.76% and 5.43%): 56% and 64% of the weighted ACS shares. Giving the NLSY rates the ACS share of bachelor's holders closes 72% of the gap for women and 18% for men (section 5).

## 1. Occupation variables, code vintages and K-12 teacher codes

The extract `raw/nlsy79_ability.parquet` already contains the CPS-job occupation variables (block `work` of `fetch_data.NLSY_BLOCKS`), so nothing was cut from the full-cohort CSV. Ages 25-34 occur only in survey years 1982-1998 (born 1957-64). Per round, with reference numbers and, for block-A person-rounds at ages 25-34 not in the armed forces, the number with a valid occupation code and with a K-12 teacher code (`80` = `CPSOCC80`, 1980 Census codes, the main classification; `70` = `CPSOCC70` through 1993 and `OCCALL-EMP.01`, job #1, 1970 codes, from 1994, the sensitivity):

| survey year | person-rounds | CPSOCC80 | 1970-code variable | valid (80) | K-12 (80) | valid (70) | K-12 (70) |
|---|---|---|---|---|---|---|---|
| 1982 | 125 | R08280.00 | R07021.00 | 110 | 1 | 110 | 1 |
| 1983 | 963 | R09453.00 | R09450.00 | 822 | 31 | 822 | 30 |
| 1984 | 1874 | R12557.00 | R12554.00 | 1598 | 42 | 1598 | 40 |
| 1985 | 2823 | R16505.00 | R16502.00 | 2430 | 77 | 2430 | 76 |
| 1986 | 3948 | R19231.00 | R19228.00 | 3447 | 95 | 3447 | 94 |
| 1987 | 5241 | R23179.00 | R23176.00 | 4622 | 117 | 4622 | 117 |
| 1988 | 6702 | R25257.00 | R25254.00 | 5932 | 142 | 5932 | 136 |
| 1989 | 7841 | R29247.00 | R29244.00 | 6874 | 173 | 6873 | 158 |
| 1990 | 8118 | R31274.00 | R31271.00 | 7099 | 162 | 7124 | 167 |
| 1991 | 8175 | R35231.00 | R35228.00 | 7087 | 153 | 7071 | 188 |
| 1992 | 7736 | R37281.00 | R37278.00 | 6661 | 178 | 6705 | 170 |
| 1993 | 6881 | R41937.00 | R41821.00 | 5865 | 177 | 5864 | 173 |
| 1994 | 5976 | R45879.01 | R45879.04 | 5020 | 152 | 5007 | 148 |
| 1996 | 3868 | R52709.00 | R52706.00 | 3449 | 113 | 3448 | 111 |
| 1998 | 1573 | R64737.00 | R64726.00 | 1401 | 41 | 1410 | 39 |

Other rounds (nobody aged 25-34; vintage from the codebook title): `CPSOCC70` 1979-1981 (1970 Census, 3-digit); `OCCALL-EMP.01` 1980-1981, 2000 (1970 Census, 3-digit); `OCCALL-EMP.01` 2002 (2000 Census, 3-digit); `OCCALL-EMP.01` 2004-2022 (2000/2002 Census, 4-digit). `OCCALL-EMP.01` before 1994 covers only a first employer other than the current job (14-144 respondents a round in 1982-92, against 7338 in 1994) and is not used. `CPSOCC80` for 1993 is the occupation of job #1 (the codebook notes that 14 first jobs are not the CPS job and carry 1970 codes; `CPSOCC70` for 1993 is the CPS item). From 1994 the CPS occupation exists for respondents with an employer for >= 10 hours a week and >= 9 weeks since the last interview.

Codebook check (58 variable-rounds): the vintage in the variable title agrees with the vintage of the code bands in the `.cdb` entry (`PROFESSIONAL, TECHNICAL AND KINDRED` = 1970, `MANAGERIAL AND PROFESSIONAL SPECIALITIES` = 1980, `Business and Financial Operations` second group = 2000 3-digit, `Executive, Administrative and Managerial` = 2002 4-digit) for 58; the extract's reference number equals the codebook's for 58 and its number of valid values equals the codebook's for 58; `ESR_COL` reads 1 = employed, 4 = in active forces in 19 of 19 rounds. The `.cdb` lists only code bands, not individual codes: the individual teacher codes below are checked against Attachment 3 and against the two codings of the same jobs.

K-12 teacher codes with the official labels of NLSY79 Codebook Supplement Attachment 3 (person-rounds = block A, ages 25-34, valid occupation; the 2000 and 2002 vintages do not occur at these ages):

| vintage | code | type | Attachment 3 label | person-rounds 25-34 |
|---|---|---|---|---|
| 1970 Census, 3-digit | 143 | pre-K and kindergarten | Prekindergarten and kindergarten teachers | 308 |
| 1970 Census, 3-digit | 142 | elementary | Elementary school teachers | 849 |
| 1970 Census, 3-digit | 144 | secondary | Secondary school teachers | 210 |
| 1970 Census, 3-digit | 145 | teachers n.e.c. / other | Teachers, except college and university, n.e.c. | 281 |
| 1980 Census, 3-digit | 155 | pre-K and kindergarten | Teachers, prekindergarten and kindergarten (231) | 303 |
| 1980 Census, 3-digit | 156 | elementary | Teachers, elementary school (232) | 776 |
| 1980 Census, 3-digit | 157 | secondary | Teachers, secondary school (233) | 229 |
| 1980 Census, 3-digit | 158 | special education | Teachers, special education (235) | 62 |
| 1980 Census, 3-digit | 159 | teachers n.e.c. / other | Teachers, n.e.c. (234, 239) | 284 |
| 2000 Census, 3-digit | 230 | pre-K and kindergarten | Preschool and Kindergarten Teachers | - |
| 2000 Census, 3-digit | 231 | elementary | Elementary and Middle School Teachers | - |
| 2000 Census, 3-digit | 232 | secondary | Secondary School Teachers | - |
| 2000 Census, 3-digit | 233 | special education | Special Education Teachers | - |
| 2000 Census, 3-digit | 234 | teachers n.e.c. / other | Other Teachers and Instructors | - |
| 2000/2002 Census, 4-digit | 2300 | pre-K and kindergarten | Preschool and Kindergarten Teachers | - |
| 2000/2002 Census, 4-digit | 2310 | elementary | Elementary and Middle School Teachers | - |
| 2000/2002 Census, 4-digit | 2320 | secondary | Secondary School Teachers | - |
| 2000/2002 Census, 4-digit | 2330 | special education | Special Education Teachers | - |
| 2000/2002 Census, 4-digit | 2340 | teachers n.e.c. / other | Other Teachers and Instructors | - |

Not counted (neighbouring codes): 1970 Census, 3-digit: 102 Agriculture teachers, 140 Teachers, college and university, subject not specified, 141 Adult education teachers; 1980 Census, 3-digit: 113 Earth, environmental, and marine science teachers (2212), 154 Postsecondary teachers, subject not specified; 2000 Census, 3-digit: 220 Postsecondary Teachers, 254 Teacher Assistants.

Concordance of the two codings on the same person-rounds (62330 person-rounds with a valid code in both, years 1982-1998; rows = 1970 code, columns = 1980 code). K-12 under both 1581, only 1970 63, only 1980 71: the codings agree on 92% of the person-rounds that either calls K-12 teaching.

| 1970 \ 1980 | 155 | 156 | 157 | 158 | 159 | other |
|---|---|---|---|---|---|---|
| 141 | 0 | 0 | 0 | 0 | 22 | 4 |
| 142 | 5 | 754 | 21 | 20 | 0 | 46 |
| 143 | 291 | 7 | 1 | 1 | 0 | 7 |
| 144 | 0 | 7 | 201 | 2 | 0 | 0 |
| 145 | 3 | 4 | 2 | 36 | 226 | 10 |
| other | 4 | 4 | 4 | 2 | 35 | 60611 |

1970 elementary (142) maps to 1980 elementary (156), pre-K/kindergarten (143) to 155, secondary (144) to 157, and teachers n.e.c. (145) and adult education (141) to 159; special education (158) has no 1970 code and comes mostly from 142 and 145.

## 2. Gaps in latent AFQT: teachers minus non-teachers

Main specification: 1980 Census codes, all K-12 teacher codes, a teacher is a K-12 teacher in the majority of his or her 25-34 rounds with a valid occupation. Latent SDs; mean latent AFQT is measured from the block-A population mean (0). Counts are unweighted persons.

**A. Teachers against all non-teachers**

| sample | group | teachers | comparison persons | mean latent AFQT, teachers | mean latent AFQT, comparison | gap (teachers - comparison) |
|---|---|---|---|---|---|---|
| employed | pooled | 160 | 8159 | 0.621 (0.060) | 0.007 (0.004) | 0.614 (0.062) [0.492, 0.731] |
|  | men | 24 | 4136 | 0.661 (0.167) | 0.019 (0.013) | 0.641 (0.167) [0.327, 0.977] |
|  | women | 136 | 4023 | 0.613 (0.065) | -0.007 (0.014) | 0.619 (0.067) [0.487, 0.756] |
| any 25-34 occupation | pooled | 165 | 8378 | 0.579 (0.064) | -0.004 (0.004) | 0.583 (0.065) [0.454, 0.715] |
|  | men | 24 | 4202 | 0.661 (0.167) | 0.011 (0.013) | 0.650 (0.167) [0.335, 0.987] |
|  | women | 141 | 4176 | 0.563 (0.070) | -0.018 (0.013) | 0.581 (0.071) [0.441, 0.723] |
| all respondents observed at 25-34 | pooled | 165 | 8611 | 0.579 (0.064) | -0.019 (0.003) | 0.598 (0.065) [0.470, 0.730] |
|  | men | 24 | 4248 | 0.661 (0.167) | 0.004 (0.013) | 0.656 (0.167) [0.340, 0.993] |
|  | women | 141 | 4363 | 0.563 (0.070) | -0.043 (0.013) | 0.605 (0.071) [0.464, 0.750] |

**B. All teachers against non-teacher college graduates (16+ years)**

| sample | group | teachers | comparison persons | mean latent AFQT, teachers | mean latent AFQT, comparison | gap (teachers - comparison) |
|---|---|---|---|---|---|---|
| employed | pooled | 160 | 1608 | 0.621 (0.060) | 0.954 (0.017) | -0.334 (0.063) [-0.455, -0.220] |
|  | men | 24 | 817 | 0.661 (0.167) | 1.011 (0.026) | -0.351 (0.167) [-0.667, -0.013] |
|  | women | 136 | 791 | 0.613 (0.065) | 0.891 (0.027) | -0.278 (0.069) [-0.415, -0.142] |
| any 25-34 occupation | pooled | 165 | 1625 | 0.579 (0.064) | 0.955 (0.017) | -0.376 (0.065) [-0.508, -0.244] |
|  | men | 24 | 821 | 0.661 (0.167) | 1.012 (0.026) | -0.351 (0.167) [-0.668, -0.014] |
|  | women | 141 | 804 | 0.563 (0.070) | 0.893 (0.027) | -0.330 (0.073) [-0.467, -0.190] |
| all respondents observed at 25-34 | pooled | 165 | 1638 | 0.579 (0.064) | 0.951 (0.017) | -0.372 (0.065) [-0.503, -0.241] |
|  | men | 24 | 826 | 0.661 (0.167) | 1.008 (0.026) | -0.347 (0.166) [-0.663, -0.013] |
|  | women | 141 | 812 | 0.563 (0.070) | 0.889 (0.027) | -0.327 (0.073) [-0.462, -0.188] |

**C. College-graduate teachers against non-teacher college graduates**

| sample | group | teachers | comparison persons | mean latent AFQT, teachers | mean latent AFQT, comparison | gap (teachers - comparison) |
|---|---|---|---|---|---|---|
| employed | pooled | 136 | 1608 | 0.716 (0.059) | 0.954 (0.017) | -0.238 (0.062) [-0.360, -0.119] |
|  | men | 23 | 817 | 0.730 (0.164) | 1.011 (0.026) | -0.281 (0.164) [-0.608, 0.034] |
|  | women | 113 | 791 | 0.713 (0.062) | 0.891 (0.027) | -0.178 (0.067) [-0.314, -0.052] |
| any 25-34 occupation | pooled | 140 | 1625 | 0.703 (0.059) | 0.955 (0.017) | -0.252 (0.062) [-0.372, -0.127] |
|  | men | 23 | 821 | 0.730 (0.164) | 1.012 (0.026) | -0.282 (0.164) [-0.608, 0.033] |
|  | women | 117 | 804 | 0.698 (0.064) | 0.893 (0.027) | -0.196 (0.068) [-0.329, -0.062] |
| all respondents observed at 25-34 | pooled | 140 | 1638 | 0.703 (0.059) | 0.951 (0.017) | -0.248 (0.062) [-0.369, -0.122] |
|  | men | 23 | 826 | 0.730 (0.164) | 1.008 (0.026) | -0.278 (0.163) [-0.605, 0.037] |
|  | women | 117 | 812 | 0.698 (0.064) | 0.889 (0.027) | -0.192 (0.068) [-0.326, -0.057] |

## 3. Sensitivities

**Employed sample**: gaps in latent SDs; A = against all non-teachers, B = against non-teacher college graduates, C = college-graduate teachers against non-teacher college graduates.

| variant | teachers (pooled) | A pooled | A men | A women | B pooled | B men | B women | C pooled |
|---|---|---|---|---|---|---|---|---|
| all K-12 codes; majority of rounds (main) | 160 | 0.614 (0.062) | 0.641 (0.167) | 0.619 (0.067) | -0.334 (0.063) | -0.351 (0.167) | -0.278 (0.069) | -0.238 (0.062) |
| all K-12 codes; at least half of rounds | 205 | 0.584 (0.063) | 0.637 (0.159) | 0.587 (0.069) | -0.368 (0.064) | -0.356 (0.159) | -0.319 (0.071) | -0.220 (0.061) |
| all K-12 codes; ever at 25-34 | 449 | 0.453 (0.046) | 0.490 (0.098) | 0.454 (0.053) | -0.531 (0.048) | -0.530 (0.099) | -0.484 (0.056) | -0.255 (0.048) |
| core (elementary, secondary, special education); majority of rounds | 106 | 0.711 (0.073) | 0.695 (0.175) | 0.725 (0.079) | -0.236 (0.074) | -0.297 (0.174) | -0.172 (0.081) | -0.224 (0.073) |
| core (elementary, secondary, special education); ever at 25-34 | 283 | 0.577 (0.051) | 0.595 (0.117) | 0.586 (0.057) | -0.406 (0.053) | -0.425 (0.117) | -0.352 (0.060) | -0.277 (0.052) |
| all K-12 codes; majority of rounds; 1970 codes | 162 | 0.602 (0.065) | 0.704 (0.162) | 0.592 (0.072) | -0.348 (0.066) | -0.289 (0.162) | -0.309 (0.073) | -0.262 (0.065) |
| all K-12 codes; majority of rounds; highest grade >= 9 | 160 | 0.586 (0.061) | 0.609 (0.167) | 0.595 (0.067) | -0.334 (0.063) | -0.351 (0.167) | -0.278 (0.069) | -0.238 (0.062) |
| all K-12 codes; majority of rounds; cross-sectional sample (SAMPLE_ID 1-8) | 110 | 0.633 (0.069) | 0.667 (0.194) | 0.642 (0.072) | -0.295 (0.070) | -0.293 (0.194) | -0.250 (0.075) | -0.229 (0.070) |
| all K-12 codes; majority of rounds; unweighted | 160 | 0.616 (0.072) | 0.752 (0.163) | 0.593 (0.081) | -0.430 (0.073) | -0.381 (0.161) | -0.365 (0.083) | -0.252 (0.064) |
| all K-12 codes; majority of rounds; mean round weights | 160 | 0.625 (0.061) | 0.649 (0.166) | 0.632 (0.067) | -0.330 (0.062) | -0.349 (0.166) | -0.273 (0.068) | -0.236 (0.061) |

**Any 25-34 occupation**: gaps in latent SDs; A = against all non-teachers, B = against non-teacher college graduates, C = college-graduate teachers against non-teacher college graduates.

| variant | teachers (pooled) | A pooled | A men | A women | B pooled | B men | B women | C pooled |
|---|---|---|---|---|---|---|---|---|
| all K-12 codes; majority of rounds (main) | 165 | 0.583 (0.065) | 0.650 (0.167) | 0.581 (0.071) | -0.376 (0.065) | -0.351 (0.167) | -0.330 (0.073) | -0.252 (0.062) |
| all K-12 codes; at least half of rounds | 195 | 0.581 (0.064) | 0.585 (0.159) | 0.594 (0.070) | -0.383 (0.065) | -0.419 (0.159) | -0.324 (0.072) | -0.255 (0.063) |
| all K-12 codes; ever at 25-34 | 505 | 0.428 (0.046) | 0.526 (0.098) | 0.415 (0.053) | -0.566 (0.048) | -0.503 (0.100) | -0.535 (0.057) | -0.237 (0.046) |
| core (elementary, secondary, special education); majority of rounds | 108 | 0.698 (0.075) | 0.704 (0.175) | 0.708 (0.083) | -0.260 (0.075) | -0.297 (0.174) | -0.204 (0.085) | -0.242 (0.075) |
| core (elementary, secondary, special education); ever at 25-34 | 318 | 0.569 (0.052) | 0.642 (0.118) | 0.565 (0.059) | -0.425 (0.055) | -0.386 (0.119) | -0.385 (0.063) | -0.252 (0.052) |
| all K-12 codes; majority of rounds; 1970 codes | 164 | 0.561 (0.066) | 0.705 (0.164) | 0.544 (0.073) | -0.400 (0.067) | -0.296 (0.164) | -0.370 (0.074) | -0.284 (0.063) |
| all K-12 codes; majority of rounds; highest grade >= 9 | 165 | 0.554 (0.065) | 0.617 (0.167) | 0.556 (0.071) | -0.376 (0.065) | -0.351 (0.167) | -0.330 (0.073) | -0.252 (0.062) |
| all K-12 codes; majority of rounds; cross-sectional sample (SAMPLE_ID 1-8) | 114 | 0.599 (0.072) | 0.672 (0.194) | 0.601 (0.078) | -0.338 (0.072) | -0.294 (0.194) | -0.302 (0.080) | -0.243 (0.070) |
| all K-12 codes; majority of rounds; unweighted | 165 | 0.577 (0.074) | 0.763 (0.163) | 0.548 (0.082) | -0.483 (0.074) | -0.381 (0.161) | -0.426 (0.084) | -0.261 (0.066) |
| all K-12 codes; majority of rounds; mean round weights | 165 | 0.586 (0.065) | 0.658 (0.166) | 0.585 (0.072) | -0.381 (0.066) | -0.349 (0.166) | -0.335 (0.073) | -0.250 (0.061) |

**No occupation required** (every respondent observed at 25-34; a person with no reported occupation is a non-teacher):

| variant | teachers (pooled) | A pooled | A men | A women | B pooled | B men | B women | C pooled |
|---|---|---|---|---|---|---|---|---|
| majority of rounds | 165 | 0.598 (0.065) | 0.656 (0.167) | 0.605 (0.071) | -0.372 (0.065) | -0.347 (0.166) | -0.327 (0.073) | -0.248 (0.062) |
| at least half of rounds | 195 | 0.597 (0.064) | 0.591 (0.159) | 0.618 (0.070) | -0.379 (0.065) | -0.415 (0.158) | -0.321 (0.072) | -0.251 (0.063) |
| ever at 25-34 | 505 | 0.444 (0.046) | 0.532 (0.098) | 0.440 (0.053) | -0.561 (0.048) | -0.498 (0.099) | -0.530 (0.057) | -0.232 (0.046) |

Ties: 45 employed-sample persons have exactly half of their rounds in K-12 teaching (non-teachers under the majority rule, teachers under 'at least half').

## 4. By type of teacher

Employed sample, majority rule within the type (a person is a type-T teacher if more than half of his or her counted rounds are in type-T codes); the comparison groups are those of the all-K-12 definition, so persons who teach mostly but not mainly in one type are in neither group. 1980 codes (main):

| type | teachers | of which women | A pooled | A women | B pooled | B women | C pooled |
|---|---|---|---|---|---|---|---|
| pre-K and kindergarten (155) | 21 | 21 | 0.050 (0.211) | 0.064 (0.212) | -0.897 (0.210) | -0.834 (0.210) | -0.634 (0.205) |
| elementary (156) | 71 | 62 | 0.715 (0.086) | 0.686 (0.095) | -0.233 (0.088) | -0.212 (0.097) | -0.215 (0.086) |
| secondary (157) | 13 | 6 | 0.941 (0.215) | 1.087 (0.273) | -0.007 (0.215) | 0.189 (0.272) | -0.007 (0.215) |
| special education (158) | 2 | 2 | 0.358 (0.260) | 0.371 (0.260) | -0.589 (0.261) | -0.526 (0.260) | -0.589 (0.261) |
| teachers n.e.c. / other (159) | 7 | 6 | 0.424 (0.399) | 0.654 (0.416) | -0.523 (0.400) | -0.243 (0.419) | -0.135 (0.418) |
| core (elementary + secondary + special education) (156+157+158) | 106 | 84 | 0.711 (0.073) | 0.725 (0.079) | -0.236 (0.074) | -0.172 (0.081) | -0.224 (0.073) |

1970 codes (`CPSOCC70` through 1993, `OCCALL-EMP.01` from 1994; no separate special education code):

| type | teachers | of which women | A pooled | A women | B pooled | B women | C pooled |
|---|---|---|---|---|---|---|---|
| pre-K and kindergarten (143) | 21 | 21 | -0.091 (0.213) | -0.078 (0.214) | -1.040 (0.213) | -0.979 (0.213) | -0.636 (0.205) |
| elementary (142) | 75 | 65 | 0.697 (0.083) | 0.691 (0.091) | -0.253 (0.084) | -0.210 (0.093) | -0.237 (0.084) |
| secondary (144) | 11 | 5 | 0.920 (0.257) | 1.079 (0.313) | -0.029 (0.256) | 0.178 (0.312) | -0.029 (0.256) |
| teachers n.e.c. / other (145) | 3 | 3 | -0.376 (0.256) | -0.363 (0.255) | -1.325 (0.255) | -1.264 (0.256) | -1.174 (0.018) |
| core (elementary + secondary) (142+144) | 108 | 86 | 0.746 (0.071) | 0.743 (0.078) | -0.203 (0.072) | -0.158 (0.079) | -0.192 (0.072) |

## 5. Teaching shares at 25-34 (sanity check against the ACS)

Round shares: weighted person-rounds of the block-A sample in an employed K-12 teacher occupation (1980 codes) over all person-rounds, employed or not (home production included); the denominator is restricted to highest grade >= 9 and not armed forces, as in the ACS sample (`acs`), or only to not armed forces (`no grade rule`). Person shares: persons classified as teachers (employed sample) over persons observed at 25-34 with highest grade >= 9. Percent.

|  | NLSY: employed K-12 / all (acs) | NLSY: no grade rule | NLSY: K-12 occupation, any status | NLSY: core only | NLSY persons, majority | NLSY persons, ever | NLSY: K-12 / employed | ACS 2009-13: project unweighted / project PWGTP / PWGTP simple rule (here) |
|---|---|---|---|---|---|---|---|---|
| men | 0.98 (0.15) | 0.95 (0.15) | 1.02 (0.16) | 0.70 (0.13) | 0.69 (0.16) | 2.75 (0.33) | 1.10 (0.17) | 1.94 / 1.76 / 1.78 |
| women | 3.45 (0.29) | 3.39 (0.29) | 4.02 (0.32) | 2.26 (0.23) | 3.30 (0.34) | 8.34 (0.53) | 4.74 (0.40) | 5.99 / 5.43 / 5.61 |

By code (employed K-12 / all, acs denominator, percent): men: pre-K and kindergarten 0.01 (0.01), elementary 0.37 (0.08), secondary 0.28 (0.07), special education 0.04 (0.02), teachers n.e.c. / other 0.27 (0.05); women: pre-K and kindergarten 0.63 (0.09), elementary 1.72 (0.19), secondary 0.40 (0.07), special education 0.14 (0.04), teachers n.e.c. / other 0.56 (0.07).

Split by schooling (NLSY: 16+ years at the round; ACS 2009-13 pooled 1-year PUMS in `raw/acs_pums`: bachelor's or more, simple rule, PWGTP weights; percent). 'NLSY rates at the ACS BA+ share' applies the NLSY K-12 rates by schooling to the ACS share of bachelor's holders; 'gap closed' is the fraction of the distance from the NLSY to the ACS simple-rule share that this closes:

|  | NLSY: 16+ share | NLSY: K-12 given 16+ | NLSY: K-12 given < 16 | ACS: BA+ share | ACS: K-12 given BA+ | ACS: K-12 given < BA | NLSY rates at the ACS BA+ share | ACS simple-rule share | gap closed |
|---|---|---|---|---|---|---|---|---|---|
| men | 24.6 (0.9) | 3.3 (0.6) | 0.22 (0.05) | 29.3 | 5.3 | 0.34 | 1.12 | 1.78 | 18% |
| women | 23.6 (0.9) | 12.2 (1.1) | 0.75 (0.09) | 37.1 | 13.0 | 1.24 | 5.00 | 5.61 | 72% |

For women the education mix closes about 72% of the gap (the K-12 rate among college graduates, 12.2 (1.1)% against 13.0%, is close to the ACS); for men it closes about 18%, and the rate among college graduates itself is lower (3.3 (0.6)% against 5.3%). The NLSY79 shares are not evidence of undercounted teachers: the occupation codes recover the ACS rate for college-educated women. What remains (the male rate among college graduates and the non-college rates) is not explained here; it may be a genuine difference between the 1957-64 cohort in 1982-98 and the 2009-13 cohorts, or a difference between the 1980 and 2002/2010 Census codings of 'other' teachers.

## Notes and caveats

- **What the three comparisons mean.** Comparison A is positive because teaching requires a college degree (92% of majority teachers against 24% of employed non-teachers have 16+ years) and college graduates score far higher on the AFQT than everyone else. The model has no schooling margin, so its non-teachers include people who would be non-college workers in the data and `z` carries whatever the AFQT measures, schooling-related or not; A is the literal counterpart of 'teachers minus non-teachers'; whether B or C is the better counterpart depends on how much of the schooling difference the model's z is meant to carry. B and C isolate selection among college graduates, where teachers are below other graduates.
- **Timing of the test.** The ASVAB was taken in 1980 at ages 15-23 and is age-normed, not grade-normed, so it partly reflects schooling attained by the test date (respondents already in college score higher than their latent ability alone implies). A compares groups with very different final schooling and may overstate the ability difference; B and C compare groups with (almost) the same final schooling and are less exposed.
- **Latent scale.** Dividing a difference of means of `afqt_c` by sqrt(R_comp_pop) assumes classical measurement error with the same error variance in every group (the convention of the wage slope in `nlsy_ability.py`).
- **Persistent teachers.** The majority rule keeps persistent teachers: 449 persons teach in at least one employed round at 25-34 but only 160 in a majority of them (teaching spells are short for many). 'Ever a teacher' gives a smaller A (0.453 (0.046) against 0.614 (0.062)) and a more negative B (-0.531 (0.048) against -0.334 (0.063)), with C about unchanged (-0.255 (0.048) against -0.238 (0.062)); the person-round teaching share in section 5 is the cross-sectional analogue of the model's equilibrium teachers.
- **Men.** Only 24 men are majority teachers in the employed sample (103 are ever teachers), so the male gaps are imprecise; by type the cells are smaller still (section 4).
- **Occupation measure.** The CPS job is the current job or, for the non-employed, the most recent job, so the 'any 25-34 occupation' sample mixes current and past jobs. Occupations are 3-digit Census codes assigned by the NLS from verbatim answers; the 1970 and 1980 codings agree on 92% of K-12 teacher person-rounds (section 1). Teacher aides, postsecondary and adult education teachers are not counted; private-school teachers are (occupation-based, like the ACS). 'Teachers n.e.c.' and pre-K/kindergarten teachers are in the main definition because they are in the ACS definition (2300-2340) that the model's teaching shares use; sections 3 and 4 show the core and the single types.
- **Weights and attrition.** 1979 weights, no correction for later attrition; a person enters the person-level classification with whatever rounds he or she was observed at 25-34 (1 to 10). SEs cluster on the NLSY79 household.
- **Teaching shares.** NLSY79 person-rounds cover 1982-98 and the 1957-64 birth cohort; the ACS benchmark is 2009-13 and the project's own shares use unweighted persons with hours-based home-production weights, while the NLSY shares are weighted person-rounds employed at the interview.
