# NLSY extensions: CFR replication in the CNLSY (T6a) and residential transitions (T3)

Generated 2026-09-29T16:13:28 by `nlsy_extensions.py`. T3: cluster bootstrap, B = 500. T6a: analytic CR1 SEs clustered by the NLSY79 mother's household, plus a household bootstrap (B = 250) for the disattenuated slopes. Seed 20260929. Data: `raw/nlsy79_ability`, `raw/nlscya_ability`, `raw/nlsy97_ability` from `fetch_data.py --sources nlsy`.

## A. T6a: CFR replication in the CNLSY

Children of NLSY79 women, Young Adult survey, earnings (wages and salary, prior calendar year) at the interview at 27-29 closest to 28 (ties to the older age), regressed on **observed** PIAT/PPVT composites at ages 8-14. Scores are rank-normed within 3-month age cells (as in `nlsy_ability.py`); the composite is the mean of the available tests (at least 2 of MATH, RECOG, COMP, PPVT), then standardized to weighted SD 1 in each regression sample, so a coefficient is per observed score SD. *Single-year*: one observation per child-round score (all rounds at ages 8-14, so the child's earnings row repeats), the CFR design. *Round average*: the child's average over rounds (no lag defined). Standard errors cluster on the mother's household. Weights: mother's 1979 weight split so each child counts equally (the `nlsy_ability.py` dyad convention); unweighted below.

### Which CFR / NLSY79 number each column is comparable to

CFR (2014b, Appendix Table 3): earnings at 28 are individual W-2 wages, zeros included (33.1% have none at 28), capped at $100,000 and expressed relative to the mean ($21,622). Unconditionally, one student SD of the score goes with $7,709, **36% of the mean** (0.305 in logs, = log 1.357). Conditioning on cubics in prior-year math and English scores and teacher fixed effects gives $2,585, **12% of the mean** (13.9% math, 10.1% English); the model's chi = log 1.12 = 0.113 is that conditional number.

- **Unconditional columns** (`base`, `controls`, round average; Table 1) have no lagged-score control. Compare them with CFR's unconditional **36%** (steps 1-2, whose definitions are CFR's) and, for the log-wage steps 5-6 disattenuated, with the NLSY79 latent slope **0.213** (controls: age, year, sex; ACS sample; no lagged score). The conditional CFR 12% is not the comparison for any of these. The `controls` columns add parental characteristics, which CFR's unconditional number does not have; they fall between CFR's unconditional and conditional numbers.
- **Lag-conditional columns** (Table 2) add a cubic in the child's PIAT math and reading scores from the previous assessment round (two years earlier when the child was tested in consecutive biennial rounds), the analogue of CFR's cubics in prior-year scores (there is no teacher fixed effect in the NLSY). Compare these with CFR's conditional **12%** (steps 1-2 only, which share CFR's definitions), and with chi = 0.113 in logs (steps 3-6 are logs, not earnings/mean, so the comparison to 12% is loose). The reference for the lag sample is `base (lag sample)`, the unconditional coefficient on the same rows, so the effect of the control is separated from the effect of the sample change. Disattenuation is not applied to lag-conditional coefficients: the lag absorbs part of the measurement error and the raw reliability no longer applies.

**Controls.** *Base*: earnings-year and age-at-earnings dummies, female (and age-at-test dummies for single-year), the specification of the NLSY79 slope. *Controls*: adds race (Hispanic, Black, other), mother's highest grade, mother's age at birth and log mother's family income (mean over the child's ages 8-14, 2010$, with a missing indicator), CFR-style parental characteristics; mother's AFQT is deliberately excluded.

**Sample.** 7016 children have a YA interview at 25-34 (19301 person-rounds); 5800 have one at 27-29. 8534 children have a composite at ages 8-14 (mean 2.73 rounds); 23221 child-round scores. Interview ages at the chosen round: {'27': 917, '28': 2579, '29': 2304}.

**Reliability** (one-factor ULS on the 4 tests in the CFR-definition sample; reliability of the available-test mean as a measure of the common factor, `comp_reliability_avail`): single-year composite 0.835, round-average composite 0.895 (weighted); 0.848 and 0.904 (unweighted). `nlsy_ability.py` reports 0.892 for the child composite over ages 5-14. The disattenuated slope is the observed slope divided by the square root of the reliability.

### Table 1 (weighted): no lagged-score control

| step | definition | single-year, base | single-year, controls | round avg, base | round avg, controls | N obs / children (single; avg) | compare with |
|---|---|---|---|---|---|---|---|
| 1 | CFR definition: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28 | 0.230 (0.014) | 0.158 (0.015) | 0.239 (0.015) | 0.166 (0.017) | 14380 / 5106; 5106 / 5106 | CFR unconditional 36% |
| 2 | Drop zeros (earnings > 0), capped earnings / mean | 0.171 (0.012) | 0.106 (0.013) | 0.180 (0.013) | 0.113 (0.015) | 12240 / 4349; 4349 / 4349 | (no CFR analogue) |
| 3 | Log earnings, earnings > 0 | 0.229 (0.018) | 0.144 (0.019) | 0.243 (0.020) | 0.154 (0.022) | 12240 / 4349; 4349 / 4349 | (no CFR analogue) |
| 4 | Log earnings, ACS-style sample | 0.171 (0.013) | 0.109 (0.015) | 0.180 (0.014) | 0.117 (0.016) | 7970 / 2831; 2831 / 2831 | (no analogue) |
| 5 | Log hourly wage, ACS-style sample | 0.160 (0.013) | 0.102 (0.014) | 0.168 (0.014) | 0.108 (0.015) | 7970 / 2831; 2831 / 2831 | NLSY79 0.213 after disattenuation |
| 6 | Log hourly wage, ACS-style, ages 25-34 pooled | 0.166 (0.011) | 0.110 (0.011) | 0.177 (0.011) | 0.120 (0.013) | 26286 / 4366; 9278 / 4366 | NLSY79 0.213 after disattenuation |
| 5d | 5, per latent SD (observed / sqrt R) | 0.175 (0.014) | 0.111 (0.015) | 0.177 (0.015) | 0.114 (0.016) |  | NLSY79 0.213 (SE 0.008) |
| 6d | 6, per latent SD (observed / sqrt R) | 0.181 (0.012) | 0.120 (0.012) | 0.187 (0.012) | 0.127 (0.013) |  | NLSY79 0.213 (SE 0.008) |

### Table 2 (weighted): with cubics in the previous round's math and reading scores (single-year design)

| step | definition | base, lag sample | lag control | controls + lag control | N obs / children (lag sample) | compare with |
|---|---|---|---|---|---|---|
| 1 | CFR definition: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28 | 0.229 (0.014) | 0.174 (0.016) | 0.118 (0.016) | 13537 / 4952 | CFR conditional 12% |
| 2 | Drop zeros (earnings > 0), capped earnings / mean | 0.173 (0.013) | 0.133 (0.014) | 0.084 (0.014) | 11552 / 4225 | (no CFR analogue) |
| 3 | Log earnings, earnings > 0 | 0.228 (0.018) | 0.175 (0.018) | 0.109 (0.018) | 11552 / 4225 | chi = 0.113 (loose) |
| 4 | Log earnings, ACS-style sample | 0.172 (0.014) | 0.126 (0.015) | 0.079 (0.015) | 7556 / 2745 | chi = 0.113 (loose) |
| 5 | Log hourly wage, ACS-style sample | 0.162 (0.013) | 0.119 (0.015) | 0.074 (0.016) | 7556 / 2745 | chi = 0.113 (loose) |
| 6 | Log hourly wage, ACS-style, ages 25-34 pooled | 0.167 (0.011) | 0.130 (0.014) | 0.088 (0.015) | 24698 / 4242 | chi = 0.113 (loose) |

### Table 1 (unweighted): no lagged-score control

| step | definition | single-year, base | single-year, controls | round avg, base | round avg, controls | N obs / children (single; avg) | compare with |
|---|---|---|---|---|---|---|---|
| 1 | CFR definition: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28 | 0.277 (0.012) | 0.185 (0.012) | 0.289 (0.012) | 0.199 (0.013) | 14380 / 5106; 5106 / 5106 | CFR unconditional 36% |
| 2 | Drop zeros (earnings > 0), capped earnings / mean | 0.199 (0.010) | 0.121 (0.010) | 0.208 (0.011) | 0.129 (0.011) | 12240 / 4349; 4349 / 4349 | (no CFR analogue) |
| 3 | Log earnings, earnings > 0 | 0.263 (0.016) | 0.166 (0.016) | 0.280 (0.017) | 0.185 (0.018) | 12240 / 4349; 4349 / 4349 | (no CFR analogue) |
| 4 | Log earnings, ACS-style sample | 0.192 (0.012) | 0.124 (0.012) | 0.205 (0.013) | 0.137 (0.015) | 7970 / 2831; 2831 / 2831 | (no analogue) |
| 5 | Log hourly wage, ACS-style sample | 0.181 (0.012) | 0.117 (0.012) | 0.193 (0.013) | 0.128 (0.014) | 7970 / 2831; 2831 / 2831 | NLSY79 0.213 after disattenuation |
| 6 | Log hourly wage, ACS-style, ages 25-34 pooled | 0.180 (0.009) | 0.118 (0.009) | 0.192 (0.010) | 0.132 (0.010) | 26286 / 4366; 9278 / 4366 | NLSY79 0.213 after disattenuation |
| 5d | 5, per latent SD (observed / sqrt R) | 0.197 (0.013) | 0.127 (0.013) | 0.203 (0.013) | 0.135 (0.015) |  | NLSY79 0.213 (SE 0.008) |
| 6d | 6, per latent SD (observed / sqrt R) | 0.196 (0.010) | 0.128 (0.010) | 0.202 (0.010) | 0.139 (0.011) |  | NLSY79 0.213 (SE 0.008) |

### Table 2 (unweighted): with cubics in the previous round's math and reading scores (single-year design)

| step | definition | base, lag sample | lag control | controls + lag control | N obs / children (lag sample) | compare with |
|---|---|---|---|---|---|---|
| 1 | CFR definition: earnings capped at $100,000 (2010$) / mean, zeros included, age ~28 | 0.274 (0.012) | 0.203 (0.012) | 0.130 (0.012) | 13537 / 4952 | CFR conditional 12% |
| 2 | Drop zeros (earnings > 0), capped earnings / mean | 0.197 (0.011) | 0.145 (0.010) | 0.083 (0.010) | 11552 / 4225 | (no CFR analogue) |
| 3 | Log earnings, earnings > 0 | 0.257 (0.016) | 0.195 (0.015) | 0.116 (0.014) | 11552 / 4225 | chi = 0.113 (loose) |
| 4 | Log earnings, ACS-style sample | 0.190 (0.012) | 0.139 (0.011) | 0.085 (0.012) | 7556 / 2745 | chi = 0.113 (loose) |
| 5 | Log hourly wage, ACS-style sample | 0.181 (0.012) | 0.130 (0.011) | 0.078 (0.011) | 7556 / 2745 | chi = 0.113 (loose) |
| 6 | Log hourly wage, ACS-style, ages 25-34 pooled | 0.181 (0.009) | 0.123 (0.009) | 0.074 (0.009) | 24698 / 4242 | chi = 0.113 (loose) |

### Disattenuated log-wage slopes with the reliability re-estimated (weighted, base spec)

| step | design | per latent SD | SE (analytic, R fixed) | SE (household bootstrap, R re-estimated) | 95% percentile CI |
|---|---|---|---|---|---|
| 5d | single | 0.175 | 0.014 | 0.014 | [0.146, 0.200] |
| 5d | avg | 0.177 | 0.015 | 0.015 | [0.147, 0.204] |
| 6d | single | 0.181 | 0.012 | 0.011 | [0.158, 0.201] |
| 6d | avg | 0.187 | 0.012 | 0.012 | [0.164, 0.207] |
B = 250. NLSY79: 0.213 (SE 0.008).

Steps 1-3 and 6 use the whole age-28 (or 25-34) earnings sample the definition allows; step 4 restricts to the ACS-style sample without changing the outcome; step 5 divides by constructed hours.

### Subject-specific single-year scores under CFR's definition (step 1)

| weights | score | base | controls | base, lag sample | lag control | controls + lag |
|---|---|---|---|---|---|---|
| weighted | composite | 0.230 (0.014) | 0.158 (0.015) | 0.229 (0.014) | 0.174 (0.016) | 0.118 (0.016) |
| weighted | PIAT math | 0.240 (0.013) | 0.174 (0.013) | 0.243 (0.013) | 0.166 (0.012) | 0.129 (0.012) |
| weighted | reading (RECOG, COMP) | 0.192 (0.014) | 0.121 (0.014) | 0.189 (0.014) | 0.103 (0.016) | 0.061 (0.015) |
| unweighted | composite | 0.277 (0.012) | 0.185 (0.012) | 0.274 (0.012) | 0.203 (0.012) | 0.130 (0.012) |
| unweighted | PIAT math | 0.270 (0.011) | 0.185 (0.011) | 0.268 (0.011) | 0.164 (0.010) | 0.121 (0.009) |
| unweighted | reading (RECOG, COMP) | 0.239 (0.011) | 0.148 (0.011) | 0.236 (0.012) | 0.134 (0.011) | 0.080 (0.011) |

### Sensitivities (weighted; specification as labelled: base, controls, lag control, controls + lag control)

| variant | design | spec | observed | per latent SD | N obs / children |
|---|---|---|---|---|---|
| rung1_uncapped | single | base | 0.241 (0.016) | 0.264 (0.017) | 14380 / 5106 |
| rung1_uncapped | single | ctrl | 0.166 (0.017) | 0.182 (0.018) | 14380 / 5106 |
| rung1_uncapped | single | lag | 0.179 (0.017) | n/a | 13537 / 4952 |
| rung1_uncapped | single | ctrl_lag | 0.120 (0.017) | n/a | 13537 / 4952 |
| rung1_uncapped | avg | base | 0.251 (0.017) | 0.265 (0.018) | 5106 / 5106 |
| rung1_uncapped | avg | ctrl | 0.175 (0.019) | 0.185 (0.020) | 5106 / 5106 |
| rung1_age26_30 | single | base | 0.228 (0.013) | 0.250 (0.015) | 16176 / 5710 |
| rung1_age26_30 | single | ctrl | 0.159 (0.014) | 0.174 (0.016) | 16176 / 5710 |
| rung1_age26_30 | single | lag | 0.170 (0.015) | n/a | 15203 / 5544 |
| rung1_age26_30 | single | ctrl_lag | 0.115 (0.015) | n/a | 15203 / 5544 |
| rung1_age26_30 | avg | base | 0.238 (0.015) | 0.252 (0.015) | 5710 / 5710 |
| rung1_age26_30 | avg | ctrl | 0.169 (0.016) | 0.178 (0.017) | 5710 / 5710 |
| rung1_interview_age28_only | single | base | 0.233 (0.018) | 0.255 (0.020) | 6637 / 2276 |
| rung1_interview_age28_only | single | ctrl | 0.165 (0.020) | 0.181 (0.022) | 6637 / 2276 |
| rung1_interview_age28_only | single | lag | 0.155 (0.022) | n/a | 6179 / 2206 |
| rung1_interview_age28_only | single | ctrl_lag | 0.107 (0.022) | n/a | 6179 / 2206 |
| rung1_interview_age28_only | avg | base | 0.243 (0.020) | 0.257 (0.021) | 2276 / 2276 |
| rung1_interview_age28_only | avg | ctrl | 0.175 (0.023) | 0.185 (0.024) | 2276 / 2276 |
| rung5_age26_30 | single | base | 0.160 (0.012) | 0.175 (0.014) | 8871 / 3131 |
| rung5_age26_30 | single | ctrl | 0.103 (0.013) | 0.112 (0.015) | 8871 / 3131 |
| rung5_age26_30 | single | lag | 0.122 (0.014) | n/a | 8399 / 3038 |
| rung5_age26_30 | single | ctrl_lag | 0.077 (0.015) | n/a | 8399 / 3038 |
| rung5_age26_30 | avg | base | 0.169 (0.013) | 0.179 (0.014) | 3131 / 3131 |
| rung5_age26_30 | avg | ctrl | 0.110 (0.015) | 0.117 (0.016) | 3131 / 3131 |
| rung5_12_months | single | base | 0.159 (0.013) | 0.174 (0.014) | 7784 / 2768 |
| rung5_12_months | single | ctrl | 0.102 (0.014) | 0.111 (0.016) | 7784 / 2768 |
| rung5_12_months | single | lag | 0.120 (0.015) | n/a | 7383 / 2684 |
| rung5_12_months | single | ctrl_lag | 0.076 (0.016) | n/a | 7383 / 2684 |
| rung5_12_months | avg | base | 0.166 (0.014) | 0.176 (0.015) | 2768 / 2768 |
| rung5_12_months | avg | ctrl | 0.108 (0.016) | 0.114 (0.017) | 2768 / 2768 |
| rung5_interview_age28_only | single | base | 0.152 (0.017) | 0.166 (0.019) | 3663 / 1251 |
| rung5_interview_age28_only | single | ctrl | 0.092 (0.019) | 0.101 (0.020) | 3663 / 1251 |
| rung5_interview_age28_only | single | lag | 0.110 (0.019) | n/a | 3430 / 1212 |
| rung5_interview_age28_only | single | ctrl_lag | 0.067 (0.019) | n/a | 3430 / 1212 |
| rung5_interview_age28_only | avg | base | 0.161 (0.019) | 0.170 (0.020) | 1251 / 1251 |
| rung5_interview_age28_only | avg | ctrl | 0.099 (0.021) | 0.104 (0.022) | 1251 / 1251 |
| rung1_lag_composite_cubic | single | lag | 0.167 (0.015) | n/a | 13633 / 4969 |
| rung1_lag_composite_cubic | single | ctrl_lag | 0.117 (0.015) | n/a | 13633 / 4969 |

### Constructed wage sample at ~28

N = 5800; zero earnings 0.147 (CFR: 0.331); job in >= 11 of 12 reference-year months 0.668; usual hours known at the interview 0.744 (of those, >= 30: 0.909); ACS-style sample 3241 (0.559); mean earnings $28107 ($27533 capped at $100,000; 0.015 above the cap), 2010$; zero earnings with a job in the year 0.021; positive earnings with no job in the history 0.047 (unweighted shares).

### Caveats for T6a

- **Annual weeks and hours are not released for the CNLSY Young Adults**, unlike NLSY79 (`WKSWK-PCY`, `HRSWK-PCY`) and NLSY97 (`CVC_*_YR_ALL`). ACS-style full-year is proxied by a job in at least 11 of the 12 reference-year months (job-history start/stop dates; about 48 weeks), full-time by usual weekly hours at all jobs current on the interview date (`TOTHOURS`, >= 30), and hourly wage by income / (52 x months/12 x usual hours). Hours at the interview can differ from hours in the reference year; the sample is therefore restricted to workers still employed at the interview.
- Earnings are wages and salary (as in NLSY79), top-coded from 2006; CFR use W-2 wages from tax records. Zeros are reported zeros (`Q15-5` universe is all YAs).
- Test ages 8-14 are PIAT/PPVT raw scores normed within age cells, not grade-level state tests. CNLSY assessments are biennial, so the 'previous round' is two years earlier (see `share_lag_gap_2` in the JSON), not the prior grade; there is no teacher fixed effect.
- The CNLSY is children of women aged 14-22 in 1979 (born mostly 1975-1992), a different population from CFR's NYC 1989-2009 cohorts. Mother-family weights assume the NLSY79 weights carry over to children.

## B. T3: residential transitions between the city and the suburb of an MSA

**Summary.** Among children in an MSA at both ages, the share living in the other type of location as adults (the model's move rate) is **31.4% (0.5) in the NLSY97** (city to suburb 32.4%, suburb to city 30.9%) and **21.6% (1.1) in the NLSY79** (city to suburb 36.1%, suburb to city 15.3%). The two are not comparable measurements of one number: the NLSY97 rate contains a mechanical reclassification from the change of MSA standard between origin and adult rounds (the switch rate roughly doubles at the 2003-04 break; about 6-8 pp), while the NLSY79 rate drops the 30-40% of SMSA residents with unknown central-city status (33-41% at ages 25-34 through 1996). Read the NLSY79 figure as the cleaner point and the NLSY97 figure as an upper bound. By parental-income tercile the move rate is nearly flat in both cohorts (below); the composition of moves by origin is not.

Origin: residence at 12-16, the round-1 (1997) NLSY97 interview (`CV_MSA`, respondents aged 12-16 at interview) or the 1979 NLSY79 interview (`SMSARES`, aged 14-16). Adult: interview rounds at 25-34; each person carries the mean of his or her round weights over the rounds observed at 25-34 in an MSA (the transition share is the weighted share of person-rounds). *City* = MSA central city, *suburb* = in MSA, not in central city (labelled 'not central city' in the codebook). The 2x2 conditions on being in an MSA at both ages and on a known central-city status. Weights: NLSY97 round `SAMPLING_WEIGHT_CC`; NLSY79 round `SAMPWEIGHT`. SEs are cluster bootstrap (persons in NLSY97, households in NLSY79); income tercile cuts are recomputed in each draw. Parental income terciles: NLSY97 `CV_INCOME_GROSS_YR` (1996 gross household income, topcoded at 2%), NLSY79 `TNFI_TRUNC` (1978 family income), over all respondents of the origin age with valid income.

### NLSY97 (main)

| row % | adult city | adult suburb |
|---|---|---|
| origin city | 67.6 | 32.4 (0.9) |
| origin suburb | 30.9 (0.7) | 69.1 |

Person-years in the 2x2: 37502 (6221 persons). Move rate (share of origin-city plus origin-suburb children living in the other type as adults): **31.4% (0.5)**. Origin share city 33.3% (0.6); adult share city 43.2% (0.6).

Origin coverage (all eligible respondents, weighted share; N): not in MSA 19.6% (1494); MSA, suburb 53.3% (4136); MSA, central city 26.0% (2697); MSA, unknown 1.1% (90); not in country 0.0% (0).

By parental-income tercile (cut points 30,500 and 57,538 in nominal dollars):

| tercile | city to suburb % | suburb to city % | move rate % | origin city share % | person-years / persons |
|---|---|---|---|---|---|
| T1 | 33.4 (1.5) | 28.0 (1.5) | 30.6 (1.0) | 47.7 | 11353 / 1864 |
| T2 | 38.6 (2.1) | 28.7 (1.3) | 31.8 (1.1) | 31.0 | 8103 / 1328 |
| T3 | 29.5 (2.2) | 34.4 (1.3) | 33.3 (1.1) | 21.7 | 7979 / 1335 |

Top minus bottom tercile: move rate 2.8 (1.6) pp, city to suburb -3.9 (2.6) pp, suburb to city 6.4 (2.0) pp.

Full 3x3 with non-MSA (row %, rows origin city / suburb / non-MSA; columns adult city / suburb / non-MSA): [66.7, 31.9, 1.4]; [30.3, 67.7, 2.0]; [30.0, 52.6, 17.4].

### NLSY79 (check)

| row % | adult city | adult suburb |
|---|---|---|
| origin city | 63.9 | 36.1 (2.2) |
| origin suburb | 15.3 (1.2) | 84.7 |

Person-years in the 2x2: 8161 (1498 persons). Move rate (share of origin-city plus origin-suburb children living in the other type as adults): **21.6% (1.1)**. Origin share city 30.2% (1.5); adult share city 30.0% (1.2).

Origin coverage (all eligible respondents, weighted share; N): not in SMSA 29.1% (1070); SMSA, suburb 34.2% (985); SMSA, central city not known 20.6% (823); SMSA, central city 15.3% (738).

By parental-income tercile (cut points 14,000 and 23,000 in nominal dollars):

| tercile | city to suburb % | suburb to city % | move rate % | origin city share % | person-years / persons |
|---|---|---|---|---|---|
| T1 | 28.4 (3.1) | 16.8 (2.5) | 22.1 (2.0) | 45.8 | 3236 / 593 |
| T2 | 42.9 (4.8) | 11.2 (2.3) | 22.6 (2.6) | 36.0 | 1827 / 329 |
| T3 | 40.7 (5.4) | 16.8 (2.2) | 20.7 (2.1) | 16.3 | 1840 / 345 |

Top minus bottom tercile: move rate -1.4 (3.0) pp, city to suburb 12.4 (6.3) pp, suburb to city 0.0 (3.4) pp.

Full 3x3 with non-MSA (row %, rows origin city / suburb / non-MSA; columns adult city / suburb / non-MSA): [58.9, 33.3, 7.9]; [13.6, 75.7, 10.6]; [11.4, 32.9, 55.8].

### Variants (move rate; city to suburb; suburb to city; in %)

| variant | move rate | city to suburb | suburb to city | person-years / persons |
|---|---|---|---|---|
| nlsy97_origin_age12_parent | 32.7 (0.6) | 34.0 (1.0) | 32.1 (0.8) | 31805 / 5202 |
| nlsy97_origin_age12_youth | 28.9 (1.8) | 28.0 (3.3) | 29.3 (2.1) | 3255 / 546 |
| nlsy97_same_region | 29.2 (0.6) | 30.5 (1.0) | 28.5 (0.7) | 33167 / 5768 |
| nlsy97_adult_first25_27 | 30.4 (0.6) | 29.7 (1.0) | 30.8 (0.8) | 16174 / 5923 |
| nlsy97_adult_near30 | 32.1 (0.6) | 33.9 (1.0) | 31.1 (0.8) | 16179 / 5865 |
| nlsy97_adult_years_2004_2011 | 30.8 (0.6) | 30.8 (1.0) | 30.7 (0.7) | 25403 / 6086 |
| nlsy97_adult_years_2013_2019 | 32.3 (0.7) | 35.5 (1.1) | 30.7 (0.8) | 12099 / 5602 |
| nlsy79_adult_years_le_1996 | 20.3 (1.1) | 36.9 (2.2) | 13.2 (1.1) | 7454 / 1412 |
| nlsy79_adult_years_1987_1990 | 16.8 (1.3) | 25.8 (2.4) | 13.0 (1.4) | 2434 / 1119 |

### Definitional breaks: round-to-round switching between city and suburb

Persons classified as city or suburb at two consecutive rounds (>= 18 at the first), weighted by the second round's weight. **NLSY97 changes MSA standard between the 2003 and 2004 rounds (to the 2000 standards; the earlier rounds carry no label in the codebook, presumably the 1990 standards) and between 2011 and 2013 (2010 standards)**; the switch rate jumps at both, and the city share jumps by 6 pp at the first. The 1997 origin is on the earlier standard and every adult round at 25-34 (2005 onward) is on a later one, so the NLSY97 transition rate includes this reclassification (the excess switching at the 2003-04 break is about 6-8 pp against the 9-12% annual switching either side).

NLSY97:

| rounds | gap (yrs) | N | switch rate % (approx SE) | share city at second round % |
|---|---|---|---|---|
| 1998-1999 | 1 | 1062 | 6.8 (0.9) | 32.1 |
| 1999-2000 | 1 | 2146 | 8.6 (0.7) | 35.2 |
| 2000-2001 | 1 | 3379 | 8.5 (0.5) | 35.2 |
| 2001-2002 | 1 | 4571 | 9.8 (0.5) | 36.3 |
| 2002-2003 | 1 | 5854 | 8.9 (0.4) | 36.6 |
| 2003-2004 | 1 | 5746 | 16.7 (0.5) | 42.9 |
| 2004-2005 | 1 | 6381 | 12.7 (0.5) | 42.6 |
| 2005-2006 | 1 | 6372 | 11.5 (0.4) | 42.4 |
| 2006-2007 | 1 | 6378 | 11.6 (0.4) | 43.4 |
| 2007-2008 | 1 | 6374 | 12.2 (0.5) | 43.5 |
| 2008-2009 | 1 | 6368 | 9.8 (0.4) | 42.7 |
| 2009-2010 | 1 | 6458 | 10.2 (0.4) | 41.6 |
| 2010-2011 | 1 | 6665 | 8.5 (0.4) | 41.9 |
| 2011-2013 | 2 | 6340 | 16.0 (0.5) | 41.6 |
| 2013-2015 | 2 | 6251 | 11.7 (0.5) | 40.7 |
| 2015-2017 | 2 | 5955 | 11.2 (0.5) | 40.1 |
| 2017-2019 | 2 | 5834 | 8.4 (0.4) | 38.1 |
| 2019-2021 | 2 | 5853 | 8.7 (0.4) | 36.9 |
| 2021-2023 | 2 | 5703 | 6.7 (0.4) | 35.9 |

NLSY79 (SMSA definitions: no spike in switching, but the 1998 round has no 'central city not known' category at all, so the classified sample changes there; the falling city share is the cohort moving to the suburbs with age):

| rounds | gap (yrs) | N | switch rate % (approx SE) | share city at second round % |
|---|---|---|---|---|
| 1979-1980 | 1 | 1931 | 5.5 (0.6) | 36.1 |
| 1980-1981 | 1 | 2469 | 8.2 (0.7) | 34.6 |
| 1981-1982 | 1 | 3017 | 6.9 (0.6) | 35.5 |
| 1982-1983 | 1 | 3031 | 5.0 (0.5) | 30.0 |
| 1983-1984 | 1 | 3261 | 5.2 (0.5) | 33.1 |
| 1984-1985 | 1 | 3325 | 3.2 (0.4) | 31.5 |
| 1985-1986 | 1 | 3138 | 4.8 (0.5) | 31.3 |
| 1986-1987 | 1 | 3105 | 6.4 (0.5) | 30.0 |
| 1987-1988 | 1 | 3246 | 3.2 (0.4) | 29.2 |
| 1988-1989 | 1 | 3171 | 5.0 (0.5) | 28.4 |
| 1989-1990 | 1 | 3260 | 3.8 (0.4) | 27.5 |
| 1990-1991 | 1 | 3379 | 2.4 (0.3) | 26.3 |
| 1991-1992 | 1 | 3470 | 1.5 (0.2) | 27.1 |
| 1992-1993 | 1 | 3241 | 5.0 (0.5) | 25.0 |
| 1993-1994 | 1 | 3764 | 3.2 (0.3) | 22.6 |
| 1994-1996 | 2 | 3533 | 5.0 (0.4) | 21.8 |
| 1996-1998 | 2 | 3833 | 8.3 (0.5) | 21.0 |

NLSY79: share of SMSA residents aged 25-34 whose central-city status is 'not known' (excluded from the 2x2):

| year | N | unknown share % |
|---|---|---|
| 1982 | 108 | 41.0 |
| 1983 | 762 | 37.5 |
| 1984 | 1539 | 40.6 |
| 1985 | 2259 | 39.8 |
| 1986 | 3124 | 38.8 |
| 1987 | 4342 | 38.0 |
| 1988 | 5373 | 40.2 |
| 1989 | 6256 | 40.6 |
| 1990 | 6514 | 39.8 |
| 1991 | 6561 | 40.4 |
| 1992 | 6162 | 37.7 |
| 1993 | 5759 | 35.9 |
| 1994 | 4913 | 34.9 |
| 1996 | 3268 | 32.7 |
| 1998 | 1325 | 0.0 |

### Caveats for T3

- **Cross-MSA movers cannot be separated** in the public files (no county/CBSA code). A child who leaves a metro area and lands in another metro's central city or suburb counts as a city/suburb transition (or a stay). The model's move rate is within one commuting zone, so a cross-MSA move that changes type inflates the rate and one that keeps type deflates it (a cross-MSA 'stayer' counts as no move). The same-region variant drops cross-region movers only (a region change requires leaving the MSA); same-region cross-MSA moves remain.
- **Central-city definitions change over time.** The NLSY97 origin (1997) is on the pre-2004 (unlabelled, presumably 1990) standard; adult rounds use 2000 standards through 2011 and 2010 standards for 2013-2019 (splits by adult year above). New CBSA definitions reclassify places between 'central city', 'suburb' and 'not in MSA', mechanically moving people across categories without moves.
- NLSY79: 'SMSA, central city not known' is 29% of 1979 SMSA residents aged 14-16 and 33-41% of SMSA residents at 25-34 through 1996 (none in 1998, when the classification became complete); those are excluded from the 2x2 (selection into unclassifiable SMSAs), and 1979 age-14-16 residence is the interview residence (at home), not residence at 14.
- Transition shares depend on the ages: residence at 25-34 is the average over person-rounds; parental income is a single noisy year, missing for a share of respondents, and terciles use the whole origin-age population, so the MSA sample within a tercile is not balanced.
