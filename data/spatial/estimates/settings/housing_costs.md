# Housing costs by school district: suburb-minus-city gaps and the residence cost r_p

Generated 2026-10-02T23:45:11 by `housing_costs.py`. ACS 5-year 2014-18 school-district tables (`acs_housing.py`, dollars of 2018), merged to the district panel and the CZ partitions of `district_bootstrap.py`; whole-CZ bootstrap, B = 1999, seed 20260929.

## Summary

- **Locale partition (suburb minus city), pooled-location median monthly housing cost (B25104 bins).** $1,224 in the city and $1,365 in the suburb (CZ-weighted means of the CZ-level medians): a gap of **$141 a month (SE 22)**, **0.114 log points (SE 0.017)**, over 188 CZs.
- **In model units.** r_p = 12 x gap / CE-2018 mean annual expenditures per consumer unit ($61,224.13, BLS) = **0.0276 (SE 0.0044)**, i.e. 2.8% of Cbar. **Range across the all-tenure dollar measures, locale partition: 0.013 (average gross rent, B25064, renter weights (b)) to 0.033 (average of district medians, B25105, occupied-unit weights (b))**; the pooled median sits inside. Owner subsets (with / without a mortgage) are lower, see the scale table. Pooled median under the revenue and ses partitions: 0.0288 and 0.0523; across the same measures and the three partitions 0.013 to 0.058.
- **Relative to earnings.** 12 x the pooled-median dollar gap is **3.3% (SE 0.5) of the CZ-pooled median full-time year-round (FTYR) earnings** ($51,323; ratio of CZ-weighted averages), 3.2% (SE 0.5) as the average of the within-CZ ratios.
- **District medians, locale (log gaps, SE).** monthly housing costs (B25105) 0.123 (0.017); gross rent (B25064) 0.073 (0.012); owner costs (B25088) 0.087 (0.025); home value (B25077) 0.096 (0.039). **Owner share** 0.509 in the city, 0.689 in the suburb: gap 0.180 (0.010).
- **Upper bound.** The whole level gap in housing spending is not a fixed cost: it mixes (i) a price per unit of housing, (ii) more or better housing (rooms, lot size, quality, tenure mix) in the suburb, which scales with income, and (iii) for owners, mortgage principal, which is saving. Treating the entire gap as the fixed, non-homothetic cost r_p x Cbar therefore gives an **upper bound** on that component; the homothetic part (richer households buy more housing wherever they live) belongs in the budget share, not in r_p.

## Data and method

- **Tables** (ACS 5-year 2014-18, school districts, three layers): B25105 median monthly housing costs; B25064 median gross rent; B25088 median selected monthly owner costs (all, with, without a mortgage); B25077 median home value; B25003 tenure (occupied, owner, renter units); B25104 distribution of monthly housing costs (15 cost bins from 'less than $100' to '$3,000 or more', plus 'No cash rent', which is excluded from the cost distribution). Census annotation codes (-666666666 etc.) are masked as missing; bottom-coded medians (<= $99, <= $9,999 of value) are dropped and top codes (4,001; 3,501 rent; 1,501 owner costs without a mortgage; 2,000,001 value) kept as lower bounds, as the earnings code does.
- **Merge and partition** exactly as the earnings tables in `district_bootstrap.py`: `leaid`, 'Remainder of <state>' rows dropped, `data_estimate.partition` (min 500 pupils, min 10% group share), `czone_gaps`, CZ weight = partitioned pupils. Within a location the median is pooled from summed bin counts (a), or the district medians are averaged with occupied-unit weights (b), or tenure counts are summed (c). Districts without a table are skipped.
- **Bin interpolation check.** The median interpolated from a district's B25104 bins reproduces its published B25105 median with a mean log error of 0.0019 (SD 0.0152, 12956 districts).
- **Bootstrap:** whole CZs resampled with replacement (`district_bootstrap.draw_counts`), B = 1999, seed 20260929; one multiplicity matrix per partition shared by every statistic; SE = SD of the draws, interval = 2.5-97.5 percentiles. The locale draws are the ones `district_bootstrap.py` uses (same n and seed). Reproduction: rebuilt `pm_ftyr_all` = 0.11223 (SE 0.01678) vs the published 0.11223 (SE 0.01678).
- **Model units.** r_p = 12 x (monthly dollar gap) / (mean annual expenditures per consumer unit, CE 2018) with the denominator **$61,224.13**: BLS Consumer Expenditure Survey, Table R-1 'All consumer units: annual detailed expenditure means', calendar year 2018 (131,439 thousand consumer units, 2.5 persons and 1.3 earners per unit, income before taxes $78,635; housing $20,090.59, of which shelter $11,747.43). Source: <https://www.bls.gov/cex/tables/calendar-year/mean/cu-all-detail-2018.pdf>; narrative report <https://www.bls.gov/opub/reports/consumer-expenditures/2018/>.

- **Denominator sensitivity.** CE total expenditures include personal insurance and pensions ($7,295.55, 11.9% of the total), saving-like items; dropping them from the denominator ($53,928.58) raises every r_p by 13.5% (pooled-median, locale: 0.0313). For scale, 12 x the pooled-median gap is 14.4% of CE shelter spending ($11,747.43).
- `locale` (locale (suburb 21-23 minus city 11-13)): 188 CZs, 3291 districts; ACS housing tables matched for 97.3% of districts (96.1% of pupils).
- `revenue` (revenue (above minus below the CZ median local revenue per pupil)): 621 CZs, 9245 districts; ACS housing tables matched for 98.9% of districts (99.6% of pupils).
- `ses` (ses (above minus below the CZ median SEDA SES)): 619 CZs, 8869 districts; ACS housing tables matched for 99.9% of districts (99.9% of pupils).

## (a) Pooled-location median monthly housing cost (B25104)

Location levels are CZ-weighted means of the CZ-level pooled medians (dollars of 2018 per month); the gap is the CZ-weighted mean of the within-CZ gap (location 2 minus location 1; city/low = 1).

| partition | loc 1 $ | loc 2 $ | gap $ (SE) | gap log (SE) | 95% interval (log) | r_p = 12 x gap / CE (SE) | 12 x gap / FTYR earnings % (SE) | n CZs |
|---|---|---|---|---|---|---|---|---|
| locale | 1,224 | 1,365 | 141 (22) | 0.114 (0.017) | [0.083, 0.147] | 0.0276 (0.0044) | 3.3 (0.5) | 188 |
| revenue | 1,060 | 1,207 | 147 (17) | 0.122 (0.011) | [0.102, 0.143] | 0.0288 (0.0033) | 3.7 (0.4) | 620 |
| ses | 1,010 | 1,277 | 267 (17) | 0.222 (0.010) | [0.203, 0.241] | 0.0523 (0.0033) | 6.7 (0.4) | 618 |

Sensitivity, pooled-location MEAN monthly housing cost from bin midpoints (top bin '$3,000 or more' valued at $3,500 / $4,000). The CE denominator is a mean, but the median gap is the better proxy for a cost every household pays; the mean gap also contains the upper tail of larger suburban homes.

| partition | top bin | gap $ (SE) | gap log (SE) | r_p (SE) |
|---|---|---|---|---|
| locale | $3,500 | 136 (24) | 0.102 (0.016) | 0.0267 (0.0047) |
| locale | $4,000 | 145 (26) | 0.105 (0.017) | 0.0283 (0.0051) |
| revenue | $3,500 | 170 (15) | 0.130 (0.009) | 0.0334 (0.0030) |
| revenue | $4,000 | 189 (17) | 0.140 (0.010) | 0.0370 (0.0034) |
| ses | $3,500 | 273 (14) | 0.211 (0.009) | 0.0536 (0.0027) |
| ses | $4,000 | 296 (16) | 0.221 (0.009) | 0.0579 (0.0031) |

## (b) Occupied-unit-weighted means of district medians

Weights B25003_001 (occupied units) within a location; CZ weights = pupils. `ln` = gap in the average log median (the earnings convention), `lvl` = log ratio of the average level, `$` = gap in the average level; `r_p` as above (not defined for the home-value stock). Weights 'tenure' (`_tw` in the JSON) weight rent by renter units and owner costs and value by owner units.

| series | partition | weights | ln gap (SE) | lvl gap (SE) | $ gap (SE) | loc 1 $ | loc 2 $ | r_p (SE) | n CZs |
|---|---|---|---|---|---|---|---|---|---|
| `B25105_001` | locale | occupied | 0.123 (0.017) | 0.132 (0.017) | 169 (24) | 1,230 | 1,399 | 0.0332 (0.0046) | 188 |
| `B25105_001` | revenue | occupied | 0.121 (0.011) | 0.122 (0.011) | 149 (17) | 1,074 | 1,223 | 0.0292 (0.0034) | 620 |
| `B25105_001` | ses | occupied | 0.226 (0.010) | 0.229 (0.010) | 278 (18) | 1,012 | 1,290 | 0.0544 (0.0035) | 618 |
| `B25064_001` | locale | occupied | 0.073 (0.012) | 0.078 (0.012) | 87 (14) | 1,149 | 1,235 | 0.0170 (0.0028) | 188 |
| `B25064_001` | revenue | occupied | 0.089 (0.009) | 0.090 (0.009) | 101 (13) | 1,005 | 1,107 | 0.0199 (0.0025) | 620 |
| `B25064_001` | ses | occupied | 0.160 (0.008) | 0.163 (0.009) | 183 (15) | 965 | 1,149 | 0.0360 (0.0030) | 618 |
| `B25064_001` | locale | tenure | 0.060 (0.012) | 0.065 (0.012) | 69 (15) | 1,139 | 1,207 | 0.0134 (0.0029) | 188 |
| `B25064_001` | revenue | tenure | 0.088 (0.009) | 0.088 (0.009) | 99 (12) | 993 | 1,092 | 0.0194 (0.0024) | 620 |
| `B25064_001` | ses | tenure | 0.155 (0.008) | 0.158 (0.008) | 178 (14) | 967 | 1,145 | 0.0349 (0.0028) | 618 |
| `B25088_001` | locale | occupied | 0.087 (0.025) | 0.095 (0.024) | 119 (36) | 1,412 | 1,531 | 0.0233 (0.0071) | 188 |
| `B25088_001` | revenue | occupied | 0.147 (0.012) | 0.147 (0.012) | 191 (22) | 1,151 | 1,341 | 0.0374 (0.0043) | 620 |
| `B25088_001` | ses | occupied | 0.229 (0.013) | 0.229 (0.013) | 287 (16) | 1,107 | 1,394 | 0.0562 (0.0032) | 618 |
| `B25088_001` | locale | tenure | 0.092 (0.024) | 0.099 (0.023) | 129 (36) | 1,421 | 1,550 | 0.0253 (0.0070) | 188 |
| `B25088_001` | revenue | tenure | 0.147 (0.012) | 0.147 (0.011) | 193 (22) | 1,160 | 1,353 | 0.0378 (0.0042) | 620 |
| `B25088_001` | ses | tenure | 0.235 (0.012) | 0.235 (0.012) | 296 (17) | 1,099 | 1,395 | 0.0580 (0.0033) | 618 |
| `B25088_002` | locale | occupied | 0.059 (0.021) | 0.066 (0.020) | 92 (38) | 1,815 | 1,907 | 0.0181 (0.0075) | 188 |
| `B25088_002` | revenue | occupied | 0.141 (0.011) | 0.143 (0.011) | 253 (28) | 1,492 | 1,746 | 0.0496 (0.0055) | 620 |
| `B25088_002` | ses | occupied | 0.189 (0.008) | 0.191 (0.008) | 322 (18) | 1,466 | 1,787 | 0.0631 (0.0036) | 618 |
| `B25088_003` | locale | occupied | 0.045 (0.022) | 0.051 (0.021) | 32 (15) | 617 | 650 | 0.0064 (0.0030) | 188 |
| `B25088_003` | revenue | occupied | 0.135 (0.011) | 0.138 (0.011) | 82 (8) | 508 | 590 | 0.0160 (0.0015) | 620 |
| `B25088_003` | ses | occupied | 0.153 (0.011) | 0.156 (0.011) | 94 (8) | 505 | 599 | 0.0184 (0.0016) | 618 |
| `B25077_001` | locale | occupied | 0.096 (0.039) | 0.114 (0.037) | 13,784 (11,941) | 300,484 | 314,268 |  | 188 |
| `B25077_001` | revenue | occupied | 0.277 (0.019) | 0.282 (0.019) | 80,046 (11,613) | 204,908 | 284,953 |  | 620 |
| `B25077_001` | ses | occupied | 0.363 (0.016) | 0.365 (0.016) | 91,694 (7,983) | 201,357 | 293,051 |  | 618 |
| `B25077_001` | locale | tenure | 0.104 (0.039) | 0.120 (0.037) | 17,195 (12,106) | 302,253 | 319,448 |  | 188 |
| `B25077_001` | revenue | tenure | 0.277 (0.019) | 0.283 (0.019) | 81,163 (11,966) | 207,233 | 288,396 |  | 620 |
| `B25077_001` | ses | tenure | 0.367 (0.015) | 0.370 (0.015) | 93,144 (8,313) | 198,735 | 291,879 |  | 618 |

## (c) Owner share

Owner-occupied over occupied units (B25003), pooled over the districts of a location.

| partition | loc 1 | loc 2 | gap (SE) | 95% interval | n CZs |
|---|---|---|---|---|---|
| locale | 0.509 | 0.689 | 0.180 (0.010) | [0.162, 0.199] | 188 |
| revenue | 0.652 | 0.650 | -0.002 (0.012) | [-0.024, 0.025] | 620 |
| ses | 0.590 | 0.714 | 0.123 (0.008) | [0.107, 0.140] | 618 |

## Scale: r_p from each dollar measure (locale partition)

| measure | $ gap per month (SE) | 12 x gap per year | r_p (SE) | 95% interval for r_p |
|---|---|---|---|---|
| pooled-location median, B25104 bins (a) | 141 (22) | 1,691 | 0.0276 (0.0044) | [0.0194, 0.0360] |
| average of district medians, B25105, occupied-unit weights (b) | 169 (24) | 2,034 | 0.0332 (0.0046) | [0.0243, 0.0424] |
| average gross rent, B25064, occupied-unit weights (b) | 87 (14) | 1,043 | 0.0170 (0.0028) | [0.0115, 0.0225] |
| average gross rent, B25064, renter weights (b) | 69 (15) | 822 | 0.0134 (0.0029) | [0.0080, 0.0194] |
| average owner costs, B25088_001, occupied-unit weights (b) | 119 (36) | 1,429 | 0.0233 (0.0071) | [0.0095, 0.0368] |
| average owner costs, B25088_001, owner weights (b) | 129 (36) | 1,551 | 0.0253 (0.0070) | [0.0118, 0.0386] |
| pooled-location MEAN from bin midpoints, top bin $3,500 | 136 (24) | 1,636 | 0.0267 (0.0047) | [0.0181, 0.0361] |
| pooled-location MEAN from bin midpoints, top bin $4,000 | 145 (26) | 1,736 | 0.0283 (0.0051) | [0.0189, 0.0386] |
| average owner costs with a mortgage, B25088_002 (b) | 92 (38) | 1,106 | 0.0181 (0.0075) | [0.0042, 0.0328] |
| average owner costs without a mortgage, B25088_003 (b) | 32 (15) | 389 | 0.0064 (0.0030) | [0.0009, 0.0127] |

Pooled-median dollar gap against the earnings scale (locale; B20005 FTYR bins), CZ-weighted means with whole-CZ bootstrap SEs:

| denominator | 12 x gap / earnings % (SE) | 95% interval % | n CZs |
|---|---|---|---|
| 12 x pooled-median dollar gap (a) / CZ-pooled median FTYR earnings, computed within each CZ and averaged with the CZ weights | 3.2 (0.5) | [2.3, 4.1] | 188 |
| 12 x (average pooled-median dollar gap) / (average CZ-pooled median FTYR earnings): ratio of the two CZ-weighted averages | 3.3 (0.5) | [2.4, 4.2] | 188 |
| 12 x (average pooled-median dollar gap) / (average city median FTYR earnings) | 3.5 (0.6) | [2.5, 4.6] | 188 |
| 12 x (average pooled-median dollar gap) / (average suburb median FTYR earnings) | 3.2 (0.5) | [2.3, 4.1] | 188 |

Earnings levels (CZ-weighted means, dollars of 2018): median FTYR earnings city $48,074, suburb $53,484, the two pooled $51,323.

## Caveats

- **Upper bound, not a point estimate of the fixed cost.** The gap in total housing *spending* bundles the price of a unit of housing, the quantity and quality of housing (suburban homes have more rooms and land), tenure (owner costs include mortgage principal, which is saving, and property taxes, part of which pay for the local schools) and the income gradient of housing demand. Only the first is the model's p_2 - p_1; the homothetic part of the second is in the preferences. Hedonic or repeat-sales price indices would give a smaller number; district medians cannot.
- **Units.** B25105/B25064/B25088 are per housing unit (household); the CE denominator is per consumer unit (2.5 persons, 1.3 earners), the model's Cbar is per adult. The ratio is a unit-free share if housing cost and consumption scale together with household size, which is the assumption used here.
- **Medians against a mean.** The numerator is a median (or an average of medians) and the denominator a mean; the share of expenditures would be larger against a median expenditure. The mean-based sensitivity in (a) brackets the numerator side.
- **Timing and geography.** ACS 2014-18 (5-year pooled, dollars of 2018) merged to the 2018 CCD panel by district id. Locations are groups of whole districts inside a 1990 commuting zone, so a 'city' district mixes in suburban neighborhoods and vice versa; the gaps are within-CZ averages weighted by pupils, so they are dominated by large CZs (see `n_eff_cz` in the JSON).
- **Revenue and ses partitions** compare the advantaged with the disadvantaged half of a CZ, not the suburb with the city; their housing gaps are price differences between rich and poor school districts.
- **B25077 (home value)** is a stock; it is shown for reference and not converted to r_p.
