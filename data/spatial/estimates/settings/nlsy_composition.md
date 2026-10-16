# NLSY79: residential composition by own ability and own income

Generated 2026-10-02T23:43:01 by `nlsy_composition.py`. NLSY79 (`raw/nlsy79_ability`), the transition machinery of `nlsy_extensions.py` (T3): same samples, weights and definitions. Household-cluster bootstrap, B = 500, seed 20260931 (= 20260929 + 2, the seed `nlsy_extensions.t3` uses for the NLSY79 main table); every table uses the same household multiplicities in every draw, re-standardizes AFQT with the draw weights and recomputes all tercile cut points. SEs are bootstrap SDs. *City* = central city of an SMSA, *suburb* = in an SMSA, not in the central city; person-rounds with 'central city not known' are excluded.

## Summary

- **Adult suburb share (25-34, SMSA with known status), overall 70.5 (0.7)%.** By own AFQT tercile: 64.2 (1.2), 76.8 (1.1), 71.1 (1.2); top minus bottom **6.9 (1.7) pp** (middle minus bottom 12.6 (1.6), top minus middle -5.7 (1.6)). By own-earnings tercile: 67.6 (1.2), 68.0 (1.2), 74.7 (1.0); top minus bottom **7.1 (1.6) pp**.
- **Transitions by own AFQT tercile** (origin 14-16, adult 25-34): top minus bottom, move rate 7.8 (2.9) pp, city to suburb 22.4 (6.0) pp, suburb to city 8.2 (2.9) pp; origin city share -22.2 (3.6) pp. The same four differences by parental-income tercile (the published table; origin city share from the same draws): -1.4 (3.0), 12.4 (6.3), 0.0 (3.4) and -29.5 (3.7) pp.
- Reproduction of the published parental-income table (statistics and SEs): maximum absolute difference 0.0e+00 (rerun of `t3_estimate` with the published seed vs the saved table), 1.0e-17 (this script's generic code vs the rerun).

## 1. Reproduction: transitions by parental-income tercile (NLSY79 main table)

Terciles of 1978 family income (`TNFI_TRUNC`) over respondents aged 14-16 in 1979 with valid income; cut points 14,000 and 23,000 (nominal dollars). Rebuilt here / published (`estimates/nlsy_extensions.md`, section B):

| tercile | city to suburb % | suburb to city % | move rate % | origin city share % | person-years / persons |
|---|---|---|---|---|---|
| T1 | 28.4 (3.1) / 28.4 (3.1) | 16.8 (2.5) / 16.8 (2.5) | 22.1 (2.0) / 22.1 (2.0) | 45.8 / 45.8 | 3236 / 593 |
| T2 | 42.9 (4.8) / 42.9 (4.8) | 11.2 (2.3) / 11.2 (2.3) | 22.6 (2.6) / 22.6 (2.6) | 36.0 / 36.0 | 1827 / 329 |
| T3 | 40.7 (5.4) / 40.7 (5.4) | 16.8 (2.2) / 16.8 (2.2) | 20.7 (2.1) / 20.7 (2.1) | 16.3 / 16.3 | 1840 / 345 |

Top minus bottom, published: move rate -1.4 (3.0) pp, city to suburb 12.4 (6.3) pp, suburb to city 0.0 (3.4) pp; rebuilt: move rate -1.4 (3.0), city to suburb 12.4 (6.3), suburb to city 0.0 (3.4).

## 2. Transitions by tercile of own AFQT

Own AFQT = `afqt_c` of `nlsy_ability.py` (equal-weight mean of the four standardized ASVAB IRT z-scores, re-standardized; 1979 weights; block-A sample, N = 9033); terciles cut over respondents aged 14-16 in 1979 with a valid score (N = 3459), cut points -0.435 and 0.501 SD. Same origin (1979 residence at 14-16), adult rounds (25-34) and weights as the table above.


| tercile | city to suburb % | suburb to city % | move rate % | origin city share % | person-years / persons |
|---|---|---|---|---|---|
| T1 | 28.4 (2.8) | 12.6 (1.9) | 19.5 (1.7) | 43.8 (2.6) | 3677 / 656 |
| T2 | 37.8 (4.6) | 11.8 (1.8) | 18.4 (1.9) | 25.1 (2.5) | 2345 / 435 |
| T3 | 50.8 (5.3) | 20.8 (2.4) | 27.3 (2.4) | 21.6 (2.6) | 1813 / 341 |

Top minus bottom tercile (pp): city to suburb 22.4 (6.0); suburb to city 8.2 (2.9); move rate 7.8 (2.9); origin city share -22.2 (3.6).

## 3. Adult (25-34) suburb share by own AFQT and by own adult earnings

Share of person-rounds in the suburb among person-rounds in an SMSA with known central-city status; every NLSY79 population-sample respondent observed at 25-34 (no origin requirement); weights = round `SAMPWEIGHT` / the person's number of classified rounds at 25-34 (city, suburb or non-SMSA; as in T3). AFQT terciles: cut over the whole block-A sample (cut points -0.433 and 0.484 SD; the origin-age cut classifies 99.4% of 14-16-year-olds identically). Earnings: mean over the person's 25-34 rounds of wage-and-salary income (Q13-5; CPI-U of the reference year, 2010 dollars; zeros included), terciles over respondents with at least one earnings observation (N = 9467), weighted by their mean round weight; cut points $17,256 and $35,707. The weighted correlation of own AFQT and own mean earnings is 0.38 (N = 8830).

Main tables:

| grouping | T1 % | T2 % | T3 % | T3 - T1 (pp) | T2 - T1 (pp) | T3 - T2 (pp) | overall % | person-years-persons by tercile |
|---|---|---|---|---|---|---|---|---|
| own AFQT | 64.2 (1.2) | 76.8 (1.1) | 71.1 (1.2) | 6.9 (1.7) | 12.6 (1.6) | -5.7 (1.6) | 70.5 (0.7) | 15146-2848 / 9674-1838 / 8072-1555 |
| own adult earnings | 67.6 (1.2) | 68.0 (1.2) | 74.7 (1.0) | 7.1 (1.6) | 0.4 (1.6) | 6.7 (1.5) | 70.5 (0.7) | 11932-2349 / 11732-2222 / 11263-2076 |

Sensitivities:

| grouping | T1 % | T2 % | T3 % | T3 - T1 (pp) | T2 - T1 (pp) | T3 - T2 (pp) | overall % | person-years-persons by tercile |
|---|---|---|---|---|---|---|---|---|
| earnings, positive earners (positive rounds only) | 67.0 (1.2) | 69.3 (1.2) | 75.0 (1.0) | 8.0 (1.6) | 2.3 (1.5) | 5.7 (1.5) | 70.5 (0.7) | 11381-2255 / 11546-2165 / 10701-1978 |
| earnings, full-year full-time wage workers (ACS-style rounds) | 68.4 (1.3) | 71.1 (1.2) | 73.6 (1.2) | 5.3 (1.8) | 2.8 (1.8) | 2.5 (1.6) | 70.5 (0.7) | 9577-1882 / 9971-1869 / 9369-1719 |
| earnings, terciles within sex | 63.0 (1.2) | 73.4 (1.1) | 73.7 (1.1) | 10.7 (1.6) | 10.4 (1.5) | 0.3 (1.5) | 70.5 (0.7) | 12074-2407 / 11143-2093 / 11710-2147 |
| own AFQT, 2x2 sample (origin at 14-16 known) | 63.4 (1.8) | 76.9 (1.9) | 72.9 (2.0) | 9.6 (2.7) | 13.5 (2.7) | -3.9 (2.7) | 71.0 (1.1) | 4620-913 / 2884-583 / 2499-501 |
| own earnings, 2x2 sample | 65.0 (2.0) | 70.1 (1.8) | 76.4 (1.7) | 11.4 (2.6) | 5.1 (2.7) | 6.3 (2.5) | 71.0 (1.1) | 3528-733 / 3583-716 / 3262-621 |

Earnings tercile cut points (2010 dollars): all rounds 17,256 / 35,707; positive 21,281 / 38,365; FYFT 28,285 / 43,818.

## Notes and caveats

- **Same machinery, different conditioning.** Sections 1-2 condition on a known SMSA status at 14-16 (in 1979 29% of 14-16-year-olds are outside any SMSA and 21% are in an SMSA whose central-city status is not known) and on being classified at 25-34; section 3 does not require an origin.
- **Origin at 14-16 is the 1979 interview residence.** Cross-SMSA movers cannot be separated (see the T3 caveats in `nlsy_extensions.md`); AFQT was measured in 1980 when respondents were 15-23, after the 1979 residence, and its terciles are age-normed scores.
- **'Central city not known'** is 33-41% of SMSA residents at 25-34 through 1996 and is dropped; if unknown status is related to ability or income the shares are selected.
- **Own earnings** are wage-and-salary income reports (top-coded in the revised variable) at 25-34 for men and women: women's zeros enter the bottom tercile unless the within-sex or positive-earner variants are used. Terciles use the person's mean over rounds, which averages out transitory income but mixes survey years.
- **SEs** cluster on the NLSY79 household (siblings share a draw); AFQT standardization and all cut points are recomputed in every draw. Differences are computed draw by draw, so their SEs include the covariance between terciles.
