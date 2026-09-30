# Literature targets for the spatial calibration (BFM, CFR, FOO)

Checks the literature inputs behind Table 3 (BFM WTP, CFR earnings effect), Table 1 (FOO class size, CFR $\chi$) and sections 2-3 of [spatial_calibration.md](../../../julia/spatial_model/calibration/spatial_calibration.md). All numbers below were read from the papers' PDFs (text extracted with `pdftotext`) or from the BLS/Census tables. Page numbers are journal pages for BFM (JPE 115(4)) and FOO (QJE 128(1)), and page numbers of the Opportunity Insights PDF for CFR (its running heads match the PDF page). No existing file was modified.

## Flags against the note (details in the sections)

1. **Test is CLAS, not CAP.** BFM's school-quality measure is the 1992-93 California Learning Assessment System (CLAS) fourth-grade math and reading score, averaged over two years. Its scale is never stated, and I could not find CLAS scale or student-level SD documentation.
2. **74 is a household-weighted SD of school averages in the full Bay Area sample (242,100 households)**, not an SD within the boundary sample where the $19.70 is estimated. BFM report no student-level SD and no between-school variance share. The conversion needs an outside ICC (section 1c).
3. **CFR's 1.3% is the smallest of three specifications.** Table 3 col. 2 gives 1.34% (SE 0.41 pp); the baseline col. 1 gives 1.65% (SE 0.43 pp). CFR's SD of teacher VA is 0.163/0.124 (math/English, elementary) and 0.134/0.098 (middle), about 0.13 on average, not 0.14/0.10.
4. **CFR's 12% per SD conditions on lagged scores.** The controls include cubics in prior-year math and English scores and teacher fixed effects. The unconditional association at age 28 is $7,709, about 36% of mean earnings. T6a (CNLSY replication) should say whether the observed-score regression includes lagged-score controls. Without them the CFR benchmark is 36% (log 0.305), not 12%.
5. **$19.70 is a mean over all households.** Households with children have a higher WTP (Table 8: +$7.41, SE 3.58, relative to no children), so "reweight to parents" raises it to about $24.6 (my computation; the SE is not computable from published tables).
6. **BFM's table numbering.** In the JPE the structural mean WTP is Table 7 col. 4 (the note is right). In the NBER working paper it is Table 6 col. 4. The JPE text says the $17.30 hedonic is in "column 2 of table 3"; it is actually Table 3, Panel B, col. 4.
7. **FOO confirmed**, but Table V reports the effect of a one-pupil reduction (positive sign) and the wage effect is significant only at 10% (t = 1.91).

## 1. Bayer, Ferreira and McMillan (2007, JPE 115(4): 588-638)

### 1a. Table and column; what the dollar amount is

- **Structural mean MWTP.** JPE **Table 7**, "Delta Regressions: Implied Mean Willingness to Pay", sample "Within 0.20 Mile of Boundary (N = 27,458)", **Panel B (including neighborhood sociodemographics), column (4) (boundary fixed effects included): "Average test score (in standard deviations) 19.7 (7.4)"** (p. 624). The same number is the first row of **Table 8**, "Mean MWTP 19.69 (7.41)" (p. 625). Text (p. 622): "the estimated mean MWTP for school quality is $19.70 per month compared with the estimated effect of $17.30 on housing prices in the analogous hedonic price regression".
- **Hedonic.** JPE Table 3 (p. 605), Panel B, col. (4): 17.3 (5.9), 0.20-mile sample, boundary FE, neighborhood sociodemographics. The text calls it "$17 per month" (p. 606).
- **Unit of the dollar amount.** Monthly user cost of housing. Table 7 note: "The dependent variable is the monthly user cost of housing, which equals monthly rent for renter-occupied units and a monthly user cost for owner-occupied housing, calculated as described in the text. Standard errors corrected for clustering at the school level are reported in parentheses." Owner values are converted to a rent equivalent with a PUMA-specific hedonic value-to-rent ratio; "the average estimate of the ratio of house values to rents is 264.1" (p. 604 fn. 18 and Data Appendix). So $19.70 per month corresponds to a house-value equivalent of $19.70 x 264.1 = $5,203 (the JPE says $17 is "approximately $4,500 in house value in 1990", and 17.3 x 264.1 = $4,569).
- **Price year.** The data are the restricted-access 1990 Census: April 1990 house values and rents, and 1989 calendar-year income. BFM say "in 1990" (p. 606) and never state a deflator. I read the $ as 1990 dollars, with income slightly earlier (1989).
- **Sample and estimator.** Six Bay Area counties (Alameda, Contra Costa, Marin, San Mateo, San Francisco, Santa Clara), boundary discontinuity design embedded in the heterogeneous sorting model. Elementary attendance zones cover 195 schools (about a third), and San Francisco is absent from the boundary samples.

### 1b. Unit of school quality and the 74-point SD

- Measure (p. 595): "As our primary measure of school quality, we use the average fourth-grade mathematics and reading test score for each school, averaged over two years".
- SD (p. 595): "The average test score, our measure of school quality, has a mean of 527 and a standard deviation of 74." Table 1 (pp. 596-597), full sample (N = 242,100 households), col. 1-2: mean 527, SD 74. In the 0.20-mile boundary sample (N = 27,548) the mean is 507, high side 544, low side 471, gap 74 (Table 1 col. 6); the 0.10-mile sample has a gap of 75 (Table 2). Text (p. 599): "the magnitude of the discontinuity is around 75 points (which is approximately a standard deviation)".
- So one SD is the household-weighted cross-school SD of the school-average score assigned to households in the full sample, not an SD within the boundary sample where the coefficient is estimated.
- Test (Data Appendix, p. 632): "The 1992-93 California Learning Assessment System (CLAS) data set provides detailed data about school performance and peer group measures. The CLAS was a test administered in the early 1990s that gives us information on student performance in mathematics, reading, and writing for grades 4, 8, and 10." It is not the California Assessment Program.
- The scale is not stated. A mean of 527 and SD of 74 with no stated scale suggest a composite across subjects (my inference only, unverified). The ICC conversion below does not need the scale.
- The Table 7 coefficient is labelled "(in standard deviations)", i.e. per 74 points.

### 1c. Converting a school-average SD to student-level SDs

- BFM report **no student-level SD** of CLAS and **no between-school variance share**. A search for CLAS scale-score documentation (1993-94) turned up nothing usable (the one ERIC report found is an image scan).
- The conversion is scale-free. If SD(school mean) / SD(student) = $r=\sqrt{\rho}$, where $\rho$ is the between-school share of variance (ICC), then one school-average SD = $r$ student-level SDs, and WTP per student-level SD = 19.70 / $r$.
- ICC evidence for grade 4 (none is California-specific; I did not find a California or Bay Area ICC):
  - Hedges and Hedberg (2007, EEPA 29(1): 60-87), national school-level ICCs for grade 4: **0.232 (math), 0.242 (reading)**, as quoted in Hedges and Hedberg (2013), ERIC ED557573, pp. 22-23. That paper's own state averages (AR, AZ, FL, KY, MA, NC, WI, not CA): **0.175 (math), 0.170 (reading)**. I read the 2013 paper; I did not read the 2007 paper's tables (image-only PDF).
  - A 2025 J. Res. Educ. Effectiveness abstract, seen only through a search snippet (page returned 403, unverified): average ICC 0.194 (math) and 0.168 (reading), grades 3-8, US states, 2009-19.
- Implied $r$: $\rho$ = 0.17 gives 0.41; 0.20 gives 0.45; 0.232 gives 0.48; 0.242 gives 0.49; 0.30 gives 0.55. A working central value is $r\approx0.45$ (range 0.41-0.49; up to 0.55 if the Bay Area is more segregated than the typical state).
- The SD of observed two-year school means also includes sampling noise, of order $(1-\rho)/n$ with $n\approx 60$-100 tested students. That raises $r$ by about 2%.
- Caveat: applying $r$ treats a school-mean difference as a difference in the average child's causal score outcome. BFM hold neighborhood sociodemographics fixed, so their coefficient is a school-quality valuation net of composition, not a peer-composition effect.

### 1d. Sample means reported (1990 Census, restricted-access)

| Variable | Value | Source |
|---|---|---|
| Mean house value (owners), full sample | $297,700 (SD $178,479) | JPE Table 1, p. 596 |
| Mean monthly rent (renters), full sample | $744 (SD $316) | JPE Table 1; text p. 595 "rents approximately $750 per month" (the NBER WP typo says "per week") |
| Owner share; rooms | 0.60; 5.11 | JPE Table 1 |
| Mean block-group income | $54,742/yr (SD $26,075), "just under $55,000 per year" | JPE Table 1, p. 595 |
| Households; persons | 242,100 households, "around 650,000 people" | JPE p. 594; 650,000 / 242,100 = 2.68 persons per household (my computation) |
| Boundary sample (0.20 mi): house value; rent; block-group income | $250,005; $678; $46,271 | JPE Table 1 col. 3 |
| Implied mean monthly user cost | about $974 | My computation: 0.60 x 297,700/264.1 + 0.40 x 744 = 676 + 298. The JPE says $17.30 is "roughly 1.8 percent of the average monthly user cost" (p. 606), and 17.3/974 = 1.78%. The JPE does not report the mean user cost. |
| Household income (mean; SD) | **$54,103/yr (SD $50,719)** | Not in the JPE. From the earlier working paper, BFM (2004), NBER w10871 ("Tiebout Sorting, Social Multipliers, and the Demand for School Quality"), Table 1 col. 1-2, same 242,100 households |
| Share of households with children under 18 | **0.333** | Same source (w10871 Table 1); not in the JPE |
| Monthly house price (rent-equivalent user cost) | $1,087 (SD $755) | Same source; probably built slightly differently. It implies 19.7/1,087 = 1.81%, so the "1.8%" also matches $19.7 |

The JPE reports no household size, share with children or household income. The values from w10871 use the same sample (mean test score 527, SD 74; block-group income $54,744 vs $54,742 in the JPE), but they are from a different version of the paper.

### 1e. Per household overall, or households with children?

- The $19.70 is a **mean over all households**. Model (p. 615): "When the household characteristics included in the model are constructed to have mean zero, $\delta_h$ is the mean indirect utility provided by housing choice h." And (p. 620) averaging over the sample "corresponds roughly to the mean MWTP of all households". The paper says the interaction characteristics are demeaned only conditionally ("when ... constructed to have mean zero"); it does not confirm the demeaning was done, so "mean over the sample" is the natural reading.
- Table 8 (p. 625) heterogeneity in MWTP per +1 SD: household income (+$10,000) **1.38 (0.33)**; children under 18 vs. no children **7.41 (3.58)**; black vs. white -14.31 (7.36); college degree or more vs. some college or less 13.03 (3.57). Text (p. 626): "The presence of children increases demand for school quality."
- Implied group means with $p=0.333$ (w10871 Table 1), my computation: with children about 19.69 + 7.41 x 0.667 = **$24.6**; without about 19.69 - 7.41 x 0.333 = $17.2. The SE for the with-children level is not computable from published tables (needs the covariance).
- The mean and the with-children WTP also embed the sample's income distribution (mean $54,103), which matters for the consumption denominator.

### 1f. Robustness of the hedonic coefficient (upside uncertainty)

JPE Table 4 (p. 609), monthly user cost per SD, 0.20-mile sample with boundary FE and sociodemographics: baseline 17.3 (5.9); with school peer and teacher measures 22.6 (8.5); block and block-group measures 19.8 (5.7); alternative block-group measure 23.8 (5.6); dropping top-coded houses 16.1 (5.7). Owner-occupied only, dependent variable house value: census value $9,376 (2,460) and transaction price $9,176 (2,738) per SD; text (p. 612): "equivalent to approximately $35 in monthly user costs, which is roughly twice the baseline estimate". The 0.10-mile sample gives 14.6 (6.3) (Table 3, col. 8). Footnote 20 (p. 605): "The low estimated value may partly reflect the informational problem households face in attempting to distinguish the quality of a school." Footnote 38 (p. 622): the direct effect "appears low", and an increase in school quality "may have an additional indirect effect on prices as households re-sort".

## 2. Household consumption and income, 1990

All BLS figures are annual, "consumer units" (households plus financially independent members), Interview + Diary combined. BLS files were read with `pdftotext` (BLS blocks direct downloads).

| Universe | Avg annual expenditures | Personal insurance and pensions | Income before taxes | Persons / children under 18 | Source |
|---|---|---|---|---|---|
| (i) All US consumer units, 1990 | **$28,381** | $2,592 | $31,889 | 2.6 / 0.7 | BLS CE 1990, Table 8 (region of residence), https://www.bls.gov/cex/1990/Standard/region.pdf |
| (ii) West region, 1990 | **$32,461** | $3,042 | $35,385 | 2.6 / 0.8 | Same Table 8 |
| (iii) "San Francisco" MSA, 1990-91 (2,520 thousand consumer units) | **$39,707** | $3,624 | $42,215 (after taxes $36,748) | 2.6 / 0.6 | BLS CE 1990-91, Table 24 (selected western MSAs), https://www.bls.gov/cex/1991/msas/west.pdf |
| (iii') "San Francisco" MSA, 1989-90 (2,550 thousand) | $38,927 | $3,533 | $40,795 (after taxes $36,543) | 2.6 / 0.6 | BLS CE 1989-90, Table 24, https://www.bls.gov/cex/1990/MSAs/west.PDF |
| West, all MSAs, 1990-91 | $32,797 | $3,129 | $36,278 | 2.7 / 0.8 | Same 1990-91 Table 24 |

- The CE MSA tables are two-year pooled samples ("1990-91", "1989-90"), not single years. There is no separate 1990-only Bay Area table.
- The "San Francisco" area is the San Francisco-Oakland-San Jose CMSA. The BLS's current published-area definition (A422) is the ten counties Alameda, Contra Costa, Marin, Napa, San Francisco, San Mateo, Santa Clara, Santa Cruz, Solano and Sonoma (https://www.bls.gov/regions/ce_areadef.pdf). I did not confirm the 1990-91 definition. It is wider than BFM's six counties, and 2.52 million consumer units is consistent with the wider area.
- Shelter and housing (SF 1990-91): housing $13,883, shelter $9,417, owned dwellings $6,113, rented $2,852; estimated monthly rental value of owned home $709; market value of owned home $151,199.
- Bay Area income, 1989. The 1990 Census (CPH-L-124, medians, 1989 income): median household income San Francisco PMSA (SF, Marin, San Mateo) **$40,494**, Oakland PMSA (Alameda, Contra Costa) $40,621, San Jose PMSA (Santa Clara) $48,115; per capita income $22,049, $18,782 and $20,423. Source: https://www2.census.gov/programs-surveys/decennial/tables/cph/cph-l/cph-l-124-h.csv (and `-p.csv` for per capita), landing page https://www.census.gov/data/tables/time-series/dec/cph-series/cph-l/cph-l-124.html. I did not find a published Census mean household income for the six counties. The best mean is BFM's own sample: **$54,103** (section 1d). Cross-check: about $20,000 per capita x 2.68 persons per household is about $54,000.
- BFM's sample mean income ($54,103) is 28% above the CE "San Francisco" pre-tax income ($42,215). CE incomes are known to be lower than Census money income, and CE consumer units include single-person units.

## 3. Chetty, Friedman and Rockoff (2014b, AER 104(9): 2633-2679)

### 3a. Earnings effect of a 1 SD teacher

- Abstract-level statement (p. 2): "At age 28, the oldest age at which we currently have a sufficiently large sample size to estimate earnings impacts, a 1 SD increase in teacher quality in a single grade raises annual earnings by 1.3%."
- **Table 3** (p. 49), OLS, one observation per student-subject-year, SEs clustered by school-cohort. The dependent variable is W-2 wage earnings at age 28 in dollars; the effect is stated as a percent by dividing by the mean of the dependent variable in the regression sample:

| Col. | Spec | Teacher VA ($) | SE | Mean of dep. var. | % of mean | SE (pp) | 95% CI |
|---|---|---|---|---|---|---|---|
| 1 | Baseline controls | 349.84 | 91.92 | 21,256 | **1.65%** | 0.43 | [0.80, 2.49] |
| 2 | + parent characteristics | 285.55 | 87.64 | 21,256 | **1.34%** | 0.41 | [0.54, 2.15] |
| 3 | + twice-lagged score controls | 308.98 | 110.17 | 21,468 | 1.44% | 0.51 | [0.43, 2.45] |
| 5 | Total income (wages + self-employment) | 353.83 | 88.62 | 22,108 | 1.60% | 0.40 | [0.81, 2.39] |

- Quotes (p. 19): "A 1 SD increase in teacher VA in a single grade increases earnings at age 28 by $350, 1.65% of mean earnings in the regression sample." "The smallest of the three estimates implies that a 1 SD increase in teacher VA raises earnings by 1.34%." CFR's Section VI uses 1.34% (Table 3, col. 2, p. 29).
- The 1.3% in the note is therefore col. 2 (or its rounding), the lowest of the three. Regression VA is the leave-year-out normalized estimate $\hat m_{jt}$, scaled so a 1 unit change is 1 SD of true teacher quality.
- The effect is on mean earnings including zeros (below).

### 3b. 1 SD of teacher VA in student test-score SDs

- p. 14: "We define $\sigma$ as the standard deviation of teacher effects for the corresponding subject and school-level using the estimates in Table 2 of our companion paper: **0.163 for math and 0.124 for English in elementary school and 0.134 for math and 0.098 for English in middle school**."
- p. 19: "A 1 SD increase in teacher quality raises end-of-year scores by **0.13 SD** of the student test score distribution on average across grades and subjects." The test scores are normalized to mean 0, SD 1 by year and grade (p. 7).
- So the prompt's 0.14/0.10 is not CFR's number. Elementary-school math is the closest to 0.14 (0.163), and middle-school math 0.134.
- The effect on scores fades: coefficient on teacher VA on the score in year t, t+1, t+2, t+3, t+4 is 0.993, 0.533, 0.362, 0.255, 0.221 (Appendix Table 10). The persistent 0.22 at t+4 is relevant to the choice of exposure $K$.

### 3c. Cross-sectional association of age-28 earnings with test scores

- **Appendix Table 3** (p. 55), column 3 (earnings at age 28), each cell a separate OLS, test scores in SD units, pooled subjects and grades, one observation per student-subject-school year:

| Row | Coefficient ($) | SE | Mean of dep. var. | % of mean |
|---|---|---|---|---|
| No controls | 7,709 | 23 | 21,622 | 35.7% (CFR text: "$7,700 (36%)", p. 12) |
| **With controls (row 2)** | **2,585** | **59** | 21,622 | **11.96%** (SE 0.27 pp) |
| Math, full controls | 2,998 | 83 | 21,622 | 13.9% |
| English, full controls | 2,192 | 88 | 21,622 | 10.1% |

- Quote (p. 19): "A 1 SD increase in student test scores, controlling for the student- and class-level characteristics $X_{it}$, is associated with a 12% increase in earnings at age 28 (Appendix Table 3, Column 3, Row 2)."
- **Controls** (p. 15): student-level "cubic polynomials in prior-year math and English scores, interacted with the student's grade level", plus ethnicity, gender, age, lagged suspensions and absences, grade repetition, free or reduced-price lunch, special education and limited English; class-level: class size and type, class and school-grade means of prior scores, class and school-year means of the covariates, grade and year dummies; **teacher fixed effects**. Row 2 of Appendix Table 3 uses this full vector. So the 12% is the association with the current score conditional on prior scores, i.e. closer to a score-gain association.
- **Zeros are included**: "Individuals with no W-2 are coded as having 0 earnings. 33.1% of individuals have 0 wage earnings at age 28 in our sample" (p. 9). Earnings are capped at $100,000 (1.3% of individuals are above it, p. 9), in 2010 dollars, and are individual W-2 wage earnings. Mean age-28 earnings in the analysis sample are $20,885 (p. 12), and $21,256 in the Table 3 regression sample.
- **By age** (Appendix Table 4, p. 56, cohorts 1979-80, with controls): the percentage effect rises 6.1% (age 20), 7.6% (22), 11.6% (23), 12.6% (24), 12.2% (25), 12.0% (26), 12.6% (27), 13.1% (28); at age 28, $2,784 (SE 171) on mean $21,320. So the age-28 association has not stopped rising, and lifetime associations may be larger.
- The note's description of the 12% (single-year observed scores, zeros in the mean, controls) matches CFR. CFR also report that the unconditional association is three times as large (36%).

### 3d. Robustness range for the earnings effect

- Across the three main specifications: 1.34-1.65% (Table 3 cols. 1-3), total income 1.60% (col. 5). All lie within one SE of each other.
- Age profile (Appendix Table 9 Panel B, p. 61), effect on earnings in dollars (SE), divided by mean earnings: age 25 $141 (44) = 0.80%; 26 $230 (47) = 1.22%; 27 $254 (63) = 1.26%; 28 $350 (92) = 1.65%. At ages 20-22 the effect is negative ($-32, $-35, $-18) because higher-VA teachers raise college attendance.
- Extensive margin (col. 4): +0.38 pp in the probability of working (SE 0.16), which CFR say is at most 81/350 = 23% of the total earnings effect (p. 20).
- Earnings growth ages 22-28 (col. 6): $286 (SE 82), "2.5%", and (p. 20) "teachers' impacts on lifetime earnings could be larger than the 1.34% impact observed at age 28".
- Consistency check by CFR (p. 20): 0.13 x 12% = 1.55%, "similar to the observed impact of 1.34%".
- Quasi-experimental estimates for earnings exist only as "very imprecise and fragile" (p. 25); CFR rely on the cross-class design.
- No separate subject or school-level earnings effects ("inadequate precision", p. 29 fn. 35).

## 4. Fredriksson, Ockert and Oosterbeek (2013, QJE 128(1): 249-285)

**Table V** (p. 273; PDF p. 26 of the UvA copy), "IV Estimates of Class Size in Fourth-Sixth Grade", one-school districts, column (1) (baseline covariates): **ln(wage) 0.0063\* (0.0033), N = 3,185**; column (2) (no baseline covariates): 0.0043 (0.0037). Confirmed.

- **Sign convention.** The table reports the effect of a one-pupil reduction in class size: "We find a 0.6% increase in wages for a pupil reduction in class size" (p. 275). So the log wage falls by 0.0063 per extra pupil, as the note says. The asterisk means significance at the 10% level (t = 1.91); the 95% CI is [-0.0002, 0.0128].
- **What is measured.** ln(wage) in full-time equivalents, wage earners only, ages 27-42, labor market outcomes averaged over 2007-09; the instrument is "Above threshold" (Maimonides-type rule), with average class size in grades 4-6. Text (p. 272): "We interpret these IV estimates as the effects of one pupil change throughout upper primary school (grades 4-6)". Same table: earnings relative to the average 0.0117\* (0.0061), P(earnings > 0) 0.0016 (0.0024).
- **Mean class size, grades 4-6: 24.357** (SD 3.489) in the one-school-district sample, N = 5,920 individuals, 191 districts; 24.066 in the full sample (Appendix Table A.1, p. 283). Confirmed.
- The note's $\sigma$ arithmetic reproduces (section 5).

## 5. Conversions

### 5a. BFM WTP as a share of consumption and income

WTP = $19.70 per month (SE $7.40), so $236.4 per year (SE $88.8). With log utility in the model, compensating consumption over consumption is the WTP over consumption, so these are the $a/C$ targets for one school-average SD. Monthly consumption is annual expenditures / 12.

| Denominator | Monthly $ | Share = 19.70 / denominator | SE of share (7.40 / denominator) |
|---|---|---|---|
| (i) US 1990 expenditures: 28,381 / 12 | 2,365.1 | **0.833%** | 0.313% |
| (ii) West 1990 expenditures: 32,461 / 12 | 2,705.1 | **0.728%** | 0.274% |
| (iii) SF 1990-91 expenditures: 39,707 / 12 | 3,308.9 | **0.595%** | 0.224% |
| (iii') SF 1989-90 expenditures: 38,927 / 12 | 3,243.9 | 0.607% | 0.228% |
| US, net of personal insurance and pensions: (28,381 - 2,592) / 12 | 2,149.1 | 0.917% | 0.344% |
| West, net: (32,461 - 3,042) / 12 | 2,451.6 | 0.804% | 0.302% |
| SF 1990-91, net: (39,707 - 3,624) / 12 | 3,006.9 | 0.655% | 0.246% |
| Income: BFM mean household income 54,103 / 12 | 4,508.6 | **0.437%** | 0.164% |
| Income: Census median household income, SF PMSA 40,494 / 12 | 3,374.5 | 0.584% | 0.219% |
| Income: CE SF pre-tax income 1990-91, 42,215 / 12 | 3,517.9 | 0.560% | 0.210% |
| Cross-check: BFM mean income times the CE SF expenditure-to-income ratio 39,707 / 42,215 = 0.941, i.e. 54,103 x 0.941 / 12 | 4,240.7 | 0.465% | 0.175% |

- The Bay Area (iii) row is the natural match for BFM's sample, but its mean expenditure is from a wider CMSA and two pooled years, and CE consumer units differ from Census households. The share ranges 0.46-0.92% across all rows; 0.6% (SF, total expenditures) to 0.73% (West) is the plausible band. The sampling SE of the WTP itself is about 38% of the estimate and dominates.
- Shares of housing cost: 19.70 / 974 = 2.0% (my implied mean user cost) and 19.70 / 1,087 = 1.8% (w10871). The hedonic 17.30 is 0.52% of SF 1990-91 expenditures (17.3 / 3,308.9).
- With children (about $24.6, section 1e), share of SF expenditures is 24.6 / 3,308.9 = 0.74%. That is only indicative, since CE consumption for households with children is higher than the all-household mean and is not tabulated here.

### 5b. Per student-level SD (uses the ICC-based $r$, section 1c)

WTP per student-level SD = 19.70 / $r$; as a share of SF 1990-91 consumption = 0.595% / $r$; implied earnings effect of a school-average SD through the bridge $\chi=0.113$ is $\chi r$.

| ICC $\rho$ | $r=\sqrt\rho$ | WTP per student SD ($/mo) | Share of SF consumption | $\chi r$ (earnings, log) |
|---|---|---|---|---|
| 0.17 | 0.412 | 47.8 | 1.44% | 4.7% |
| 0.20 | 0.447 | 44.1 | 1.33% | 5.1% |
| 0.232 | 0.482 | 40.9 | 1.24% | 5.4% |
| 0.242 | 0.492 | 40.0 | 1.21% | 5.6% |
| 0.30 | 0.548 | 36.0 | 1.09% | 6.2% |

### 5c. CFR effect and SE in log points

Level percentage $x=\hat\beta/\bar y$ (percent of mean earnings including zeros). Log points: $\ln(1+x)$, SE by the delta method $\mathrm{se}(x)/(1+x)$. Arithmetic: $x=285.55/21{,}256=0.013434$; $\ln(1.013434)=0.013344$; $\mathrm{se}(x)=87.64/21{,}256=0.004123$; $\mathrm{se}/(1+x)=0.004069$.

| Estimate | $x$ | SE of $x$ | log points | SE (log) |
|---|---|---|---|---|
| Table 3 col. 2 (the note's 1.3% basis, 1.34%) | 0.013434 | 0.004123 | **0.01334** | **0.00407** |
| Table 3 col. 1 (baseline) | 0.016458 | 0.004325 | 0.01632 | 0.00425 |
| Table 3 col. 3 | 0.014392 | 0.005131 | 0.01429 | 0.00506 |
| Note's rounded 1.3% | 0.013 | n/a | 0.01292 (the 0.01291 in the note) | n/a |

- Using 1.34% instead of the rounded 1.3% raises the log effect by about 3%. Using the baseline 1.65% raises it by about 26%.
- Per student SD of teacher-induced score gain: 1.34 / 0.13 = 10.3% (baseline 1.65 / 0.13 = 12.7%), against the 12% conditional association (0.13 x 11.96% = 1.55%).
- Cross-sectional bridge: conditional $x=2{,}585/21{,}622=0.11955$, $\ln(1.11955)=0.11293$ (SE $0.00273/1.11955=0.00244$), consistent with $\chi=\ln1.12=0.1133$. Unconditional: $x=7{,}709/21{,}622=0.3565$, $\ln(1.3565)=0.305$. Math with controls: 13.9% (0.130); English: 10.1% (0.097).

### 5d. FOO check of $\sigma$ (Table 1)

$\sigma=(1-\eta)(K/3)(24.357)(0.0063)$ at $\eta=0.080$, $K=12$: $0.92\times4\times24.357\times0.0063=0.5647$ (note: 0.565). 95% CI: $0.92\times4\times24.357\times(0.0063\mp1.96\times0.0033)$ = $[-0.0151,\,1.1444]$ (note: $[-0.015,1.144]$). At $K=5$: $0.92\times(5/3)\times24.357\times0.0063=0.2353$ (note: 0.235). All reproduce.

## 6. Open items

- No CLAS scale or student-level SD source found; the school-average to student-level conversion relies on non-California ICCs (0.17-0.24 for grade 4). A California grade-4 ICC (for example from SEDA or the CDE) would tighten $r$.
- The JPE does not report household income, household size or share with children; the household values above come from the earlier NBER working paper w10871 (same sample) and from BFM's household count.
- The CE Bay Area denominators are two-year CMSA averages, not a single 1990 or six-county figure.
- The 1.8% and mean user-cost figures are my reconstruction from JPE Table 1; the JPE reports only the 1.8% ratio.

## Sources

- BFM (2007, JPE), JSTOR copy: https://matthewturner.org/ec2410/readings/Bayer_Ferreira_Mcmillan_JPE_2007.pdf ; DOI https://doi.org/10.1086/522381
- BFM working paper (NBER w13236, July 2007): https://www.nber.org/papers/w13236 (PDF https://www.nber.org/system/files/working_papers/w13236/w13236.pdf)
- BFM (2004, NBER w10871), household characteristics in Table 1: https://www.nber.org/papers/w10871 (PDF https://www.nber.org/system/files/working_papers/w10871/w10871.pdf)
- BLS CE 1990 Table 8 (region): https://www.bls.gov/cex/1990/Standard/region.pdf
- BLS CE 1990-91 Table 24 (western MSAs): https://www.bls.gov/cex/1991/msas/west.pdf ; 1989-90: https://www.bls.gov/cex/1990/MSAs/west.PDF
- BLS CE published-area definitions: https://www.bls.gov/regions/ce_areadef.pdf
- Census CPH-L-124 (1989 income, metro areas): https://www.census.gov/data/tables/time-series/dec/cph-series/cph-l/cph-l-124.html ; median household income CSV https://www2.census.gov/programs-surveys/decennial/tables/cph/cph-l/cph-l-124-h.csv ; per capita https://www2.census.gov/programs-surveys/decennial/tables/cph/cph-l/cph-l-124-p.csv
- CFR (2014b, AER): https://opportunityinsights.org/wp-content/uploads/2018/03/teachers2.pdf ; DOI https://doi.org/10.1257/aer.104.9.2633
- FOO (2013, QJE): https://pure.uva.nl/ws/files/1912854/124460_394782.pdf ; DOI https://doi.org/10.1093/qje/qjs048
- Hedges and Hedberg (2013), ICC values for state samples and quotation of the 2007 national values: https://files.eric.ed.gov/fulltext/ED557573.pdf ; Hedges and Hedberg (2007): https://journals.sagepub.com/doi/10.3102/0162373707299706 (tables not read)
