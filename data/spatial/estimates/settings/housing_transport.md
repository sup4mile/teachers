# Housing market and transport evidence for income-driven demand (task D2, DESIGN.md M1 / M1b)

Generated 2026-10-05T02:56:00 by `housing_transport.py`. ACS 5-year 2014-18 school-district tables (dollars of 2018), merged to the district panel and the locale partition of `district_bootstrap.py` (city = NCES locale 11-13, suburb = 21-23; base year 2018, at least 500 pupils per district, each location at least 10% of the CZ's partitioned pupils; 188 CZs, 3291 districts). Location statistics pool the districts of a location within a CZ; gaps are suburb minus city per CZ and are averaged over CZs with the partitioned pupils as weights; whole-CZ bootstrap, B = 1999, seed 20260929: the draws of `housing_costs.py` and `district_settings.py`, shared by every statistic here. CE = BLS Consumer Expenditure Survey, calendar year 2018 (mean annual expenditures per consumer unit $61,224.13, the unit of the model's Cbar).

**Revision after review (2026-10-05).** Compared with the first version of this report: (b) the central `alphaH` / `hmin` come from the fit on expenditure per equivalent adult (square-root scale), 0.178 / 0.0971 against 0.201 / 0.0725 at the consumer-unit level, with the ranges widened to 0.12-0.27 / 0.04-0.14; `housing_share` is defined as a ratio of aggregates and the recommended terciles come from the deciles; (c) eps_H applies one horizon convention to every anchor: city 1.93 / suburb 2.69 (ranges 0.71-3.37 / 1.10-4.70) against 1.26 / 1.76; (d) r_car is standardized for household size, 0.0132 (0.0065-0.0406) against 0.0189; (a) the validation moments for the lower tail are the renters' (the bottom of the all-units distribution is mostly outright owners); (e) associate-or-more and children-term variants are added. The partition, the bootstrap draws and every other statistic are unchanged.

## Summary

- **(a) Housing-cost distribution.** The median gap reproduces `housing_costs.md`: **0.1137 (0.0169)** log points. All units, log gap (SE) at the 10th / 25th / 50th / 75th / 90th percentile: 0.111 (0.033) / 0.082 (0.018) / 0.114 (0.017) / 0.119 (0.018) / 0.111 (0.020); p10 minus median -0.003 (0.021): **pooled over tenure the suburb does not lack cheap units**. The pooled distribution mixes three groups with very different cost levels, and its bottom is mostly outright owners, not renters: of the units below the CZ-pooled 10th percentile 60% (city) and 77% (suburb) are owned without a mortgage, 37% and 19% rented; their cash costs (taxes, insurance, utilities) are not the price of the housing consumed, so the all-units and owner profiles do not measure R H at the bottom. **Gross rent is the validation moment for the lower tail of R_r H_i**: renters' log gap 0.198 (0.027) / 0.089 (0.012) / 0.057 (0.013) / 0.038 (0.013) / 0.023 (0.013) (p10 minus median 0.140 (0.023); dollar gap $115 at p10, $22 at p90): the suburb lacks cheap rentals, not a uniform fixed charge; rental units below the CZ-pooled 25th percentile 28.4% city vs 21.8% suburb (gap -6.5 (0.9) pp). The all-units 90th percentile is partly extrapolated (19 of 376 CZ-location cells fall in the open $3,000+ bin): 0.096-0.132 across treatments.
- **(b) Housing share and Stone-Geary fit.** CE 2018 shares of average annual expenditures by quintile of income before taxes (lowest to highest): shelter 24.1 / 20.5 / 19.9 / 18.9 / 17.5%; utilities, fuels and public services 9.1 / 8.5 / 7.6 / 6.7 / 4.9%; **shelter + utilities 33.2 / 29.1 / 27.5 / 25.6 / 22.4%** (all consumer units 25.8%); housing 40.0 / 35.8 / 34.5 / 32.1 / 29.9% (32.8%). Household size rises with income across the CE groups (people per consumer unit 1.7 / 2.1 / 2.5 / 2.9 / 3.2), so the fall of the share with expenditure partly reflects economies of scale in shared housing, not only nonhomotheticity; the model's household is one adult. **Central fit, expenditure per equivalent adult (square-root scale): aH = 0.178 (SE 0.007), hmin = 0.0971 (SE 0.0066) of mean expenditure per equivalent adult ($39,118; h = $3,799 per equivalent adult, $5,945 for a household of mean size)**, RMSE 0.09 pp (deciles 0.178 and 0.0968). Other scales (aH; hmin): none 0.201, 0.0725; square root of persons 0.178, 0.0971; OECD-modified 0.170, 0.1062; per adult 0.154, 0.1226; per capita 0.141, 0.1356 (none = the consumer-unit fit of the first version). The curve fits best at an equivalence elasticity of 0.58 (quintiles) / 0.42 (deciles) (two degrees of freedom: suggestive), in line with the square root. Envelope over scales, definitions, quintiles / deciles and dropping the lowest group: aH 0.123-0.274, hmin 0.036-0.140. **The model's `housing_share` must be a ratio of aggregates** (sum R H / sum Y'): **0.258 (0.004)** in the data, terciles of income 0.316 / 0.272 / 0.233 (deciles); a mean of household-level ratios is larger (the mean of the five quintile ratios, 0.276, is a lower bound). ACS cross-check (national, 2018 1-year, household income, no size adjustment): 31.9% of households pay 30% or more of income; the same form fitted to mean burdens gives aH = 0.196 (0.183-0.220 over the open-class assumptions).
- **(c) Supply elasticities.** One horizon convention for every anchor: the long-run stock elasticity eps_T / (1 - (1 - delta)^T) at delta = 3.5% (Baum-Snow and Duranton, slide 76), a factor 1.52 for Saiz's 30-year and 3.34 for the 10-year estimates. Saiz (2010): population-weighted mean 1.55 over the 95 MSAs of Table VI (2.54 unweighted over 269), long run **2.37** (1.00 from Baum-Snow and Han's 10-year tract elasticities, 3.86 from the unweighted mean). Period check (Baum-Snow and Han, Online Appendix D): Saiz's method on 2000-10 data gives a mean of 1.25 against 2.54 for 1970-2000, a ratio 0.49; the horizon alone predicts 0.46 at delta = 3.5% (4.5% would match), so the lower 2000-10 estimates are mostly the shorter horizon and no further period adjustment is made. Suburban tracts are 1.40 times as elastic as central-city tracts (all variants 1.25-1.73: units / floorspace, four estimators, city shares 0.25-0.5, within-region medians). **Recommended long-run eps_H: city 1.93 (0.71-3.37), suburb 2.69 (1.10-4.70)**, metro average 2.37 (1.00-3.86); constructed from the literature, not estimated. Without the horizon conversion (Saiz read as long run, the first version) 1.26 / 1.76.
- **(d) Transport.** Zero-vehicle households: 13.3% city vs 5.3% suburb (gap -8.0 (1.0) pp); vehicles per household 1.56 vs 1.94 (gap 0.377 (0.030)). Within CZ, controlling for log median household income, the suburb coefficient is 0.245 (0.025) vehicles, **but suburban households are larger** (one-person households 32% of city and 24% of suburban households; 4+ persons 22% vs 27%), and the model's household is one adult. Within household-size classes the income-adjusted gaps are 0.130 / 0.167 / 0.193 / 0.213 vehicles (1 / 2 / 3 / 4+ persons); **size-standardized and income-adjusted: 0.172 (0.022)**; one-person households 0.130 (0.017); with household-size shares as controls 0.085 (0.014). **Car requirement as a fixed goods cost: r_car = 0.0132 (SE 0.0017) of mean consumption** ($811 a year; 0.0065 to 0.0406); the first version's 0.0189 did not adjust for household size (-30%); a single adult reading (one-person households, one-person consumption) gives 0.0169. Vehicles are a normal good (vehicle count elasticity 0.65 with respect to expenditures across CE income groups).
- **(e) BFM income gradient.** Elasticity of WTP for school quality with respect to household income at the means: **0.379 (SE 0.091; 0.169 if the SE of the mean MWTP is also counted)**. Folding in the college-degree term through the 1990 Bay Area association (PUMS: 35.5% of householders with a bachelor's degree or more, 43.1% with an associate degree or more (BFM's block groups: 44%), correlation with household income 0.32, slope 0.0358 per $10,000): **0.508 (SE 0.097)**; associate-or-more 0.507; adding the children term 0.530; 0.49-0.54 across all variants.

## Recommended inputs for the model

| Input | Value (SE) | Range | Basis |
|---|---|---|---|
| `gap_housing_cost` (target for a_H) | 0.1137 (0.0169) |  | log gap of the pooled median monthly housing cost, suburb - city (B25104 bins): reproduces housing_costs.md (0.114, SE 0.017) |
| lower-tail validation moments for R_r H_i: renters' log gap at p10 / p25 / p50 / p75 / p90 | 0.198 (0.027) / 0.089 (0.012) / 0.057 (0.013) / 0.038 (0.013) / 0.023 (0.013) |  | gross rent, B25063; dollar gap 115 / 77 / 59 / 43 / 22; gross rent is the price of the unit; the bottom of the all-units and owner distributions is mostly outright owners whose cash costs (taxes, insurance, utilities) do not measure the housing consumed |
| `alphaH` (marginal budget share of housing) | 0.178 (0.007) | 0.123 - 0.274 | Stone-Geary NLS on CE 2018 quintile means of shelter + utilities / total expenditures, expenditure per equivalent adult (square-root scale); no size adjustment (consumer-unit level) 0.201 |
| `hmin` (units of Cbar at R = 1) | 0.0971 (0.0066) | 0.036 - 0.140 | subsistence housing expenditure per equivalent adult over mean expenditure per equivalent adult (= that of a household of any size if housing and consumption scale together with size): the model's hmin in units of Cbar at R = 1 (h = $3,799 per equivalent adult of $39,118; $5,945 for a household of mean size); no size adjustment 0.0725 |
| `housing_share` (sum R H / sum Y'), ratio of aggregates | 0.258 (0.004) | household-level mean of ratios larger; lower bound 0.276 | ratio of aggregates: aggregate shelter + utilities / aggregate total expenditures, CE 2018, consumer-unit weights (sum R H / sum Y' over adults in the model); one-person consumer units 0.317, two or more persons 0.246 (CE Table 1400) |
| `housing_share_by_income`, terciles of income before taxes (ratio of aggregates within tercile) | 0.316 (0.006) / 0.272 (0.004) / 0.233 (0.004) |  | CE deciles with fractional weights (quintile-based: 0.311 / 0.272 / 0.233); SE: CE sampling error, upper bound; consumer-unit income terciles, so household-size differences are included |
| `epsH` city | 1.93 | 0.71 - 3.37 | long run (delta = 3.5%): Saiz's population-weighted mean over 95 MSAs x 1.52 = 2.37, split by the Baum-Snow-Han suburb/city ratio; without the horizon conversion 1.26 |
| `epsH` suburb | 2.69 | 1.10 - 4.70 | ratio suburb/city 1.40 (all variants 1.25-1.73); sensitivity grid 0.5, 1, 2, 4, elastic; without the horizon conversion 1.76 |
| `r_car` (fixed goods cost of the suburban car requirement, units of Cbar) | 0.0132 (0.0017) | 0.0065 - 0.0406 | size-standardized, income-adjusted extra vehicles x CE private-vehicle spending per vehicle / mean expenditures per consumer unit; single-adult reading 0.0169; extensive margin only 0.0031; first version (no household-size adjustment) 0.0189 |
| `wtp_income_elasticity` | 0.379 (0.091) |  | BFM Table 8: 1.38 per $10,000 at mean income $54,103 and mean WTP $19.69 |
| `wtp_income_elasticity`, college term folded in | 0.508 (0.097) | 0.490 - 0.538 | plus 13.03 x slope of the householder college share on income (1990 Bay Area PUMS); associate-or-more 0.507, with the children term 0.530 |

How the pieces map into M1 / M1b (DESIGN_v2.md): the supply shifter a_H is fitted to `gap_housing_cost` (the all-units median, as in `housing_costs.md`); the **renter** percentile profile and renters' units-below-the-cutoff shares of (a) are the validation moments for the lower tail of R_r H_i (at the bottom the all-units and owner profiles are dominated by outright owners' cash costs, which do not measure housing consumed); `alphaH`, `hmin` and the tercile `housing_share` come from (b), all for shelter + utilities and, for the first two, in units of expenditure per equivalent adult (the tercile shares are validation of the fitted curve; they are consumer-unit income terciles and so include household-size differences); `epsH` is a long-run elasticity for the stationary equilibrium, with the grid of (c) for sensitivity; `r_car` is a separate suburban fixed cost (ACS housing costs do not contain vehicles, so it does not double count the housing gap), standardized for household size because suburban households are larger; the WTP elasticity targets the model's `wtp_income_elasticity`.

**Model definitions (what the model must compute to be comparable with these inputs).** Preferences (1 - aH) log c + aH log(H - hmin) over non-housing consumption c and housing H; Y' is the household's outlay on c and housing (income after taxes, education goods, moving, residence-charge and other fixed goods costs), C = Y' - R hmin the supernumerary expenditure. Then H = hmin + aH C / R, R H = R hmin + aH C, c = (1 - aH) C and Y' = C + R hmin (DESIGN_v2: the line of DESIGN.md M1b 'Y' = (1 + aH) C + R hmin' double counts aH C; it holds only for preferences log c + aH log(H - hmin), whose marginal housing share is aH / (1 + aH)). aH is the MARGINAL budget share of housing, the quantity estimated by the Stone-Geary fit s(y) = aH + (1 - aH) h / y; hmin is h over mean expenditure per equivalent adult (units of Cbar at R = 1). `housing_share` must be computed as a RATIO OF AGGREGATES, sum_i R H_i / sum_i Y'_i over the adults (and `housing_share_by_income` as the same ratio within terciles of pre-tax income): the CE targets are aggregate shelter + utilities over aggregate expenditures. A mean over adults of R H_i / Y'_i (M1B_SPEC) is a different moment: because the share falls with expenditure it is larger (the consumer-unit-weighted mean of the five CE quintile shares, 0.276, is a lower bound for a household-level mean of ratios, against the aggregate 0.258). Units: dollar amounts per household (r_car, the housing-cost gap) are divided by the mean expenditure per consumer unit, $61,224, as in housing_costs.md, which is unit-free if per-household costs and consumption scale together with household size; hmin is already a ratio of two amounts in the same unit (h per equivalent adult over mean expenditure per equivalent adult), so the dollar h per equivalent adult must not be divided by $61,224 (a household of mean size has h = hmin x $61,224).

## (a) Housing-cost distribution by location

Pooled monthly housing cost (dollars of 2018) by location: bin counts summed over the districts of a location within each CZ, percentiles by linear interpolation inside the bin (the rule of `housing_costs.grouped_median`), gap = log(suburb) - log(city) per CZ, then the CZ-weighted mean. All units: B25104 (renters' gross rent and owners' selected monthly owner costs in 15 bins, top bin $3,000 or more; units with no cash rent excluded). Renters: B25063 gross rent (24 bins, top $3,500 or more). Owners: B25087 selected monthly owner costs, with a mortgage (17 bins, top $4,000 or more) and without (19 bins, top $1,500 or more), pooled on the union of the two grids; the $1,500-or-more no-mortgage count is spread over the union bins at or above $1,500 like the with-mortgage owners of the partitioned districts (all CZs pooled) (shares 0.351, 0.238, 0.154, 0.094, 0.057, 0.106 over 1,500, 2,000, 2,500, 3,000, 3,500, 4,000). When a percentile falls in an open top bin it is extrapolated with a Pareto tail through the survival function at the lower edge of the top bin and of the nearest non-empty bin below it (shape at least 1); the number of CZ-location cells that needed this is in the last column and the alternatives are in the sensitivity table. Bin check: B25104 bins sum to 72.91 million units = renters (B25063) 27.39 + owners (B25087) 45.51 million.

**Reproduction.** The all-units median gap is 0.11374 (0.01686) here and 0.11374 (0.01686) in `housing_costs.json` (max CZ-level difference 0.0e+00); the earnings target `gap_income` of `district_settings.py` is 0.11223 (0.01678) on the same draws (reported 0.1122, 0.0168).

### all occupied units (B25104)

| percentile | city $ | suburb $ | gap $ (SE) | log gap (SE) | 95% interval | log gap minus median gap (SE) | top-bin cells (city / suburb) |
|---|---|---|---|---|---|---|---|
| 10 | 480 | 543 | 63 (20) | 0.111 (0.033) | [0.047, 0.174] | -0.003 (0.021) | 0 / 0 |
| 25 | 806 | 877 | 71 (17) | 0.082 (0.018) | [0.048, 0.117] | -0.032 (0.005) | 0 / 0 |
| 50 | 1,224 | 1,365 | 141 (22) | 0.114 (0.017) | [0.083, 0.147] |  | 0 / 0 |
| 75 | 1,782 | 1,993 | 211 (36) | 0.119 (0.018) | [0.085, 0.156] | 0.006 (0.005) | 1 / 3 |
| 90 | 2,541 | 2,816 | 275 (58) | 0.111 (0.020) | [0.074, 0.152] | -0.003 (0.010) | 8 / 11 |

### renters, gross rent (B25063)

| percentile | city $ | suburb $ | gap $ (SE) | log gap (SE) | 95% interval | log gap minus median gap (SE) | top-bin cells (city / suburb) |
|---|---|---|---|---|---|---|---|
| 10 | 559 | 674 | 115 (16) | 0.198 (0.027) | [0.147, 0.254] | 0.140 (0.023) | 0 / 0 |
| 25 | 837 | 915 | 77 (11) | 0.089 (0.012) | [0.067, 0.114] | 0.032 (0.007) | 0 / 0 |
| 50 | 1,133 | 1,192 | 59 (15) | 0.057 (0.013) | [0.033, 0.083] |  | 0 / 0 |
| 75 | 1,504 | 1,547 | 43 (21) | 0.038 (0.013) | [0.014, 0.066] | -0.019 (0.004) | 0 / 0 |
| 90 | 1,951 | 1,973 | 22 (27) | 0.023 (0.013) | [-0.001, 0.052] | -0.034 (0.007) | 0 / 0 |

### owners, selected monthly owner costs (B25087)

| percentile | city $ | suburb $ | gap $ (SE) | log gap (SE) | 95% interval | log gap minus median gap (SE) | top-bin cells (city / suburb) |
|---|---|---|---|---|---|---|---|
| 10 | 455 | 512 | 57 (17) | 0.106 (0.029) | [0.051, 0.162] | 0.018 (0.014) | 0 / 0 |
| 25 | 778 | 845 | 67 (23) | 0.089 (0.028) | [0.035, 0.146] | 0.001 (0.008) | 0 / 0 |
| 50 | 1,402 | 1,513 | 110 (38) | 0.088 (0.025) | [0.040, 0.136] |  | 0 / 0 |
| 75 | 2,130 | 2,213 | 84 (51) | 0.059 (0.022) | [0.017, 0.105] | -0.029 (0.012) | 0 / 1 |
| 90 | 2,999 | 3,041 | 42 (87) | 0.034 (0.025) | [-0.012, 0.085] | -0.054 (0.019) | 6 / 4 |

### Owners by mortgage status (log gap, SE)

| percentile | owners with a mortgage (B25087) | owners without a mortgage (B25087) |
|---|---|---|
| 10 | 0.093 (0.029) | 0.113 (0.039) |
| 25 | 0.079 (0.021) | 0.076 (0.023) |
| 50 | 0.058 (0.020) | 0.043 (0.021) |
| 75 | 0.037 (0.023) | 0.010 (0.022) |
| 90 | 0.010 (0.027) | -0.024 (0.030) |

### Units below the CZ-pooled percentile

For each CZ the cut is the 10th / 25th / 50th percentile of the two locations pooled; the entry is the share of each location's units below it (equal shares of 10 / 25 / 50% would mean the same distribution in both locations). Gap = suburb minus city, in percentage points.

| distribution | cut (percentile) | city % | suburb % | gap pp (SE) |
|---|---|---|---|---|
| all occupied units (B25104) | 10 | 11.9 | 9.5 | -2.4 (0.7) |
| all occupied units (B25104) | 25 | 28.0 | 23.8 | -4.2 (0.9) |
| all occupied units (B25104) | 50 | 55.1 | 47.1 | -8.0 (1.2) |
| renters, gross rent (B25063) | 10 | 12.1 | 8.3 | -3.8 (0.4) |
| renters, gross rent (B25063) | 25 | 28.4 | 21.8 | -6.5 (0.9) |
| renters, gross rent (B25063) | 50 | 52.7 | 46.6 | -6.1 (1.3) |
| owners, selected monthly owner costs (B25087) | 10 | 12.2 | 9.4 | -2.8 (0.7) |
| owners, selected monthly owner costs (B25087) | 25 | 27.7 | 24.1 | -3.6 (1.0) |
| owners, selected monthly owner costs (B25087) | 50 | 53.9 | 48.4 | -5.4 (1.4) |

### Sensitivity to the top-bin treatment (log gap)

| distribution | treatment | p10 | p25 | p50 | p75 | p90 | CZs |
|---|---|---|---|---|---|---|---|
| all occupied units (B25104) | Pareto tail (main) | 0.111 | 0.082 | 0.114 | 0.119 | 0.111 | 188 |
|  | CZs with an undetermined percentile dropped | 0.111 | 0.082 | 0.114 | 0.120 | 0.132 | 177 (p90) |
|  | percentile set at the bin's lower edge | 0.111 | 0.082 | 0.114 | 0.114 | 0.096 | 188 |
|  | common sample, no extrapolated cell (177 CZs) | 0.104 | 0.089 | 0.131 | 0.139 | 0.132 | 177 |
|  | finer grid: renters (B25063) + owners (B25087) pooled, Pareto tail above $4,000 | 0.112 | 0.080 | 0.116 | 0.126 | 0.108 | 188 |
| renters, gross rent (B25063) | Pareto tail (main) | 0.198 | 0.089 | 0.057 | 0.038 | 0.023 | 188 |
|  | CZs with an undetermined percentile dropped | 0.198 | 0.089 | 0.057 | 0.038 | 0.023 | 188 (p90) |
|  | percentile set at the bin's lower edge | 0.198 | 0.089 | 0.057 | 0.038 | 0.023 | 188 |
|  | common sample, no extrapolated cell (188 CZs) | 0.198 | 0.089 | 0.057 | 0.038 | 0.023 | 188 |
| owners, selected monthly owner costs (B25087) | Pareto tail (main) | 0.106 | 0.089 | 0.088 | 0.059 | 0.034 | 188 |
|  | CZs with an undetermined percentile dropped | 0.106 | 0.089 | 0.088 | 0.057 | 0.068 | 181 (p90) |
|  | percentile set at the bin's lower edge | 0.106 | 0.089 | 0.088 | 0.058 | 0.044 | 188 |
|  | common sample, no extrapolated cell (181 CZs) | 0.114 | 0.117 | 0.115 | 0.094 | 0.068 | 181 |

### Who occupies the cheapest units

The pooled distribution mixes three groups with very different cost levels, and their mix differs by location. For each CZ the cut is the 10th / 25th / 50th percentile of the all-units distribution of the two locations pooled; the entries are the shares of the units below the cut that are rented (cash rent, B25063), owned with a mortgage or owned without a mortgage (B25087), each group's count below the cut taken by linear interpolation inside its own bins; the occupancy rows are shares of occupied units (cash renters plus owners). CZ-weighted means, whole-CZ bootstrap, gap = suburb minus city in percentage points.

| group | population | city % | suburb % | gap pp (SE) |
|---|---|---|---|---|
| renters (B25063 gross rent, cash renters) | all occupied units | 48.4 | 30.2 | -18.2 (1.0) |
| owners with a mortgage (B25087) | all occupied units | 34.0 | 47.7 | 13.8 (0.7) |
| owners without a mortgage (B25087) | all occupied units | 17.7 | 22.1 | 4.4 (0.4) |
| renters (B25063 gross rent, cash renters) | units below the pooled 10th percentile (about $528) | 36.6 | 19.3 | -17.3 (1.3) |
| owners with a mortgage (B25087) | units below the pooled 10th percentile (about $528) | 3.1 | 3.6 | 0.5 (0.3) |
| owners without a mortgage (B25087) | units below the pooled 10th percentile (about $528) | 60.3 | 77.2 | 16.8 (1.3) |
| renters (B25063 gross rent, cash renters) | units below the pooled 25th percentile (about $853) | 44.8 | 25.4 | -19.5 (1.1) |
| owners with a mortgage (B25087) | units below the pooled 25th percentile (about $853) | 7.3 | 8.5 | 1.2 (0.5) |
| owners without a mortgage (B25087) | units below the pooled 25th percentile (about $853) | 47.9 | 66.1 | 18.3 (1.0) |
| renters (B25063 gross rent, cash renters) | units below the pooled 50th percentile (about $1,312) | 55.3 | 37.4 | -17.9 (0.8) |
| owners with a mortgage (B25087) | units below the pooled 50th percentile (about $1,312) | 15.5 | 20.2 | 4.8 (0.6) |
| owners without a mortgage (B25087) | units below the pooled 50th percentile (about $1,312) | 29.2 | 42.4 | 13.2 (0.7) |

Monthly cost by group, city / suburb in dollars (CZ-weighted means of the percentiles):

| group | p10 | p25 | p50 | p75 | p90 |
|---|---|---|---|---|---|
| renters (B25063) | 559 / 674 | 837 / 915 | 1,133 / 1,192 | 1,504 / 1,547 | 1,951 / 1,973 |
| owners with a mortgage (B25087) | 975 / 1,071 | 1,313 / 1,413 | 1,807 / 1,884 | 2,471 / 2,515 | 3,392 / 3,359 |
| owners without a mortgage (B25087) | 291 / 337 | 429 / 470 | 612 / 639 | 860 / 862 | 1,231 / 1,174 |

**Reading.** For the pooled distribution (all units) the log gap is about the same at the 10th, 50th and 90th percentile and lower at the 25th: no bottom-heavy gap. But the pooled gap is not a clean measure of the price of housing at the bottom: the cheapest units are mostly owned outright, not rented. Of the units below the CZ-pooled 10th percentile, 60% (city) and 77% (suburb) are owned without a mortgage and 37% and 19% are rented (below the 25th percentile: 48% / 66% and 45% / 25%); mortgage owners are a few percent of the bottom. The suburb has more outright owners (22% of occupied units against 18% in the city) and more mortgage owners (48% against 34%) but far fewer renters (30% against 48%), so a gap between pooled percentiles combines within-group price differences with composition. The cash cost of an outright owner is property tax, insurance and utilities, not the price of the housing consumed, so the bottom of the all-units and owner distributions does not measure R_r H_i. **Gross rent does** (a rent is the price of the unit): the renter log gap falls from 0.198 at the 10th to 0.023 at the 90th percentile and in dollars from $115 to $22, so the suburb lacks cheap rentals rather than facing a uniform fixed charge. **The renter (gross rent) profile and renters' units-below-the-cutoff shares are the validation moments for the lower tail of R_r H_i**. The all-units median gap remains the target for the supply shifter (it is the `housing_costs.md` target), but it mixes the tenure groups too (units owned without a mortgage are 29% of the city's and 42% of the suburb's units below the pooled median): the like-for-like renter median gap is 0.057 (0.013), about half of 0.114. The all-units and owner profiles are descriptive only.

## (b) Housing expenditure shares by income and the Stone-Geary fit

Source: BLS Consumer Expenditure Survey, calendar year 2018, **Table 1101** (quintiles of income before taxes) and **Table 1110** (deciles), 'Annual expenditure means, shares, standard errors, and coefficients of variation' (https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error/cu-income-quintiles-before-taxes-2018.xlsx and ...-deciles-before-taxes-2018.xlsx; PDFs cached beside them). Shares are item means over mean total expenditures of the group (they reproduce the published one-decimal shares). Consumer units (CUs), not adults: see the household-size subsection; dollars of 2018.

### Quintiles of income before taxes

| group | CUs (000) | income before taxes $ | expenditures $ | shelter $ | shelter % | utilities $ | utilities % | shelter + utilities $ | shelter + utilities % (SE) | shelter + utilities ex. telephone % | housing $ | housing % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Lowest 20 percent | 26,333 | 11,285 | 26,399 | 6,355 | 24.1 | 2,415 | 9.1 | 8,770 | 33.2 (0.9) | 30.6 | 10,553 | 40.0 |
| Second 20 percent | 26,315 | 31,237 | 39,968 | 8,212 | 20.5 | 3,404 | 8.5 | 11,616 | 29.1 (0.7) | 26.3 | 14,293 | 35.8 |
| Third 20 percent | 26,196 | 54,900 | 51,729 | 10,270 | 19.9 | 3,955 | 7.6 | 14,225 | 27.5 (0.5) | 24.8 | 17,860 | 34.5 |
| Fourth 20 percent | 26,221 | 90,478 | 69,131 | 13,069 | 18.9 | 4,633 | 6.7 | 17,702 | 25.6 (0.5) | 23.1 | 22,217 | 32.1 |
| Highest 20 percent | 26,374 | 204,975 | 118,781 | 20,814 | 17.5 | 5,836 | 4.9 | 26,650 | 22.4 (0.5) | 20.7 | 35,501 | 29.9 |
| All consumer units | 131,439 | 78,635 | 61,224 | 11,747 | 19.2 | 4,049 | 6.6 | 15,796 | 25.8 (0.4) | 23.5 | 20,091 | 32.8 |

The SE of a share (percentage points, in parentheses) is the delta-method SE of item mean over total with the published SEs of the CE means and independent sampling errors (an upper bound; with a correlation of 0.5 between item and total the SE of the shelter + utilities share is 0.7 pp for the lowest quintile). The computed shares reproduce the published one-decimal shares to within 0.05 pp. The shelter + utilities SE treats the two items as independent. Mean of the five group shares (a mean-of-ratios counterpart weighted by consumer units): shelter + utilities 27.6%.

CE 'Utilities, fuels, and public services' includes telephone services (about a third of it); 'shelter + utilities ex. telephone' (shelter, natural gas, electricity, fuel oil and other fuels, water and other public services) covers what an ACS gross rent or owner cost covers and is the closest match to the ACS-based `gap_housing_cost` target; the requested shelter + utilities remains the central definition.

### Deciles of income before taxes (shares of expenditures, %)

| group | income before taxes $ | expenditures $ | shelter | utilities | shelter + utilities | shelter + utilities ex. telephone | housing |
|---|---|---|---|---|---|---|---|
| Lowest 10 percent | 5,724 | 25,309 | 23.6 | 8.4 | 32.0 | 29.4 | 39.0 |
| Second 10 percent | 16,848 | 27,488 | 24.5 | 9.8 | 34.4 | 31.6 | 40.9 |
| Third 10 percent | 26,306 | 37,164 | 22.1 | 8.6 | 30.7 | 28.0 | 37.2 |
| Fourth 10 percent | 36,167 | 42,771 | 19.2 | 8.4 | 27.6 | 24.9 | 34.6 |
| Fifth 10 percent | 47,738 | 49,241 | 19.7 | 7.9 | 27.6 | 24.9 | 34.8 |
| Sixth 10 percent | 62,083 | 54,223 | 19.9 | 7.4 | 27.4 | 24.7 | 34.3 |
| Seventh 10 percent | 79,250 | 64,029 | 19.2 | 7.0 | 26.1 | 23.6 | 32.7 |
| Eighth 10 percent | 101,729 | 74,236 | 18.7 | 6.5 | 25.2 | 22.7 | 31.6 |
| Ninth 10 percent | 138,383 | 95,056 | 16.9 | 5.7 | 22.6 | 20.4 | 29.9 |
| Highest 10 percent | 271,773 | 142,554 | 17.9 | 4.4 | 22.3 | 20.8 | 29.9 |

### Stone-Geary fit on consumer-unit expenditure (no household-size adjustment): s(Y) = h/Y + aH (1 - h/Y), Y = total expenditures of the group

Nonlinear least squares in (h, aH) on the group shares, equal weight per group. The form is linear in (aH, (1 - aH) h) after writing s = aH + (1 - aH) h / Y, so the NLS solution equals that regression; a Levenberg-Marquardt fit from a different start reproduces it (largest relative difference 2e-15). h/mean = h over CE mean expenditures per consumer unit ($61,224.13). SE (CE) = parametric bootstrap of the CE sampling error (item mean and total drawn normal with the published SEs, independent across groups and of each other: conservative); SE (resid) = delta method from the residual variance (n - 2 degrees of freedom). Elasticity = aH / s(mean): the elasticity of housing expenditure with respect to total expenditures at the mean.

| groups | definition | aH (SE, CE) | SE resid | h $ | h / mean exp. (SE, CE) | SE resid | s at mean % | elasticity at mean | RMSE pp | max abs resid pp | R2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| quintiles | shelter | 0.156 (0.005) | 0.003 | 2,570 | 0.0420 (0.0048) | 0.0027 | 19.2 | 0.82 | 0.26 | 0.50 | 0.986 |
| quintiles | utilities | 0.045 (0.001) | 0.007 | 1,440 | 0.0235 (0.0014) | 0.0052 | 6.7 | 0.67 | 0.54 | 0.74 | 0.866 |
| quintiles | shelter + utilities | 0.201 (0.006) | 0.006 | 4,436 | 0.0725 (0.0054) | 0.0047 | 25.9 | 0.78 | 0.44 | 0.66 | 0.985 |
| quintiles | shelter + utilities ex. telephone | 0.181 (0.005) | 0.003 | 4,036 | 0.0659 (0.0050) | 0.0020 | 23.5 | 0.77 | 0.19 | 0.27 | 0.997 |
| quintiles | housing | 0.273 (0.008) | 0.004 | 4,659 | 0.0761 (0.0080) | 0.0039 | 32.9 | 0.83 | 0.33 | 0.64 | 0.991 |
| deciles | shelter | 0.156 (0.005) | 0.006 | 2,583 | 0.0422 (0.0044) | 0.0049 | 19.1 | 0.81 | 0.77 | 1.45 | 0.892 |
| deciles | utilities | 0.046 (0.001) | 0.006 | 1,385 | 0.0226 (0.0013) | 0.0041 | 6.8 | 0.68 | 0.70 | 1.42 | 0.784 |
| deciles | shelter + utilities | 0.202 (0.005) | 0.008 | 4,387 | 0.0717 (0.0050) | 0.0067 | 25.9 | 0.78 | 1.03 | 2.05 | 0.922 |
| deciles | shelter + utilities ex. telephone | 0.182 (0.005) | 0.008 | 4,005 | 0.0654 (0.0046) | 0.0061 | 23.5 | 0.77 | 0.95 | 1.69 | 0.924 |
| deciles | housing | 0.274 (0.007) | 0.007 | 4,623 | 0.0755 (0.0077) | 0.0061 | 32.9 | 0.83 | 0.86 | 1.62 | 0.940 |

Fit of the central definition (shelter + utilities, quintiles):

| group | expenditures $ | data share % | fitted share % | residual pp | subsistence h/Y % |
|---|---|---|---|---|---|
| Lowest 20 percent | 26,399 | 33.22 | 33.54 | -0.32 | 16.8 |
| Second 20 percent | 39,968 | 29.06 | 28.98 | 0.08 | 11.1 |
| Third 20 percent | 51,729 | 27.50 | 26.97 | 0.53 | 8.6 |
| Fourth 20 percent | 69,131 | 25.61 | 25.24 | 0.37 | 6.4 |
| Highest 20 percent | 118,781 | 22.44 | 23.10 | -0.66 | 3.7 |

Robustness (aH; h / mean expenditures):

| definition | quintiles | quintiles, lowest dropped | quintiles, net of pensions | deciles | deciles, lowest dropped | deciles, net of pensions |
|---|---|---|---|---|---|---|
| shelter | 0.156; 0.0420 | 0.161; 0.0357 | 0.195; 0.0278 | 0.156; 0.0422 | 0.153; 0.0459 | 0.194; 0.0280 |
| shelter + utilities | 0.201; 0.0725 | 0.194; 0.0809 | 0.250; 0.0580 | 0.202; 0.0717 | 0.192; 0.0839 | 0.250; 0.0570 |
| shelter + utilities ex. telephone | 0.181; 0.0659 | 0.179; 0.0684 | 0.226; 0.0528 | 0.182; 0.0654 | 0.173; 0.0754 | 0.226; 0.0522 |
| housing | 0.273; 0.0761 | 0.269; 0.0816 | 0.340; 0.0486 | 0.274; 0.0755 | 0.266; 0.0861 | 0.340; 0.0481 |

'Net of pensions' deflates by expenditures excluding personal insurance and pensions (mean $53,928), a closer match to consumption (these entries are in units of that smaller mean). Dropping the lowest group removes consumer units whose expenditures exceed income (students, retirees, business losses).

### Household size and the income gradient: fits on equivalized expenditure

The CE counts consumer units (CUs), and household size rises with income: people per CU 1.7 / 2.1 / 2.5 / 2.9 / 3.2 and earners 0.5 / 0.8 / 1.3 / 1.8 / 2.1 from the lowest to the highest quintile. The fall of the housing share with expenditure therefore mixes nonhomotheticity with economies of scale in shared housing, and the model's household is one adult. With an equivalence scale e = n^theta the Stone-Geary share of a household is s = aH + (1 - aH) h1 / y with y = Y / e (equivalized expenditure) and h1 the subsistence housing expenditure per equivalent adult. h1 over the mean expenditure per equivalent adult (CE mean expenditures per CU $61,224.13 over the CU-weighted mean scale) is the subsistence share of a household of ANY size if housing and consumption scale together with household size, the unit-free assumption of `housing_costs.md`; it is the model's hmin in units of Cbar. Scales (theta = 0 no adjustment, 1 per capita; Foerster and Mira d'Ercole (OECD), Table 2): square root of persons (the OECD income-distribution convention, elasticity 0.5; central); OECD-modified (1 + 0.5 (adults - 1) + 0.3 children, elasticity about 0.53; adults = persons minus children under 18, the CE's child age limit instead of the OECD's 14); per adult (children count zero: the model's literal 'one household per adult', DESIGN_v2 mean expenditure about $32,000); per capita. People per CU are published to one decimal and treated as exact. Shelter + utilities, equal weight on the groups, SE (CE) as above.

| group | people | children under 18 | adults | earners | expenditure per sqrt-equivalent $ | per OECD-modified $ | per adult $ | per capita $ |
|---|---|---|---|---|---|---|---|---|
| Lowest 20 percent | 1.7 | 0.3 | 1.4 | 0.5 | 20,247 | 20,464 | 18,856 | 15,529 |
| Second 20 percent | 2.1 | 0.5 | 1.6 | 0.8 | 27,581 | 27,564 | 24,980 | 19,032 |
| Third 20 percent | 2.5 | 0.6 | 1.9 | 1.3 | 32,716 | 31,736 | 27,226 | 20,692 |
| Fourth 20 percent | 2.9 | 0.7 | 2.2 | 1.8 | 40,595 | 38,194 | 31,423 | 23,838 |
| Highest 20 percent | 3.2 | 0.8 | 2.4 | 2.1 | 66,401 | 61,227 | 49,492 | 37,119 |
| mean scale (CU-weighted) / mean expenditure per equivalent adult $ |  |  |  |  | 1.565 / 39,118 | 1.624 / 37,701 | 1.900 / 32,224 | 2.480 / 24,688 |

Central definition (shelter + utilities), Stone-Geary fit by scale:

| groups | scale | aH (SE, CE) | h $ per equivalent adult | mean expenditure per equivalent adult $ | h / mean (SE, CE) | s at mean % | elasticity at mean | RMSE pp | R2 |
|---|---|---|---|---|---|---|---|---|---|
| quintiles | none (consumer unit) | 0.201 (0.006) | 4,436 | 61,224 | 0.0725 (0.0054) | 25.9 | 0.78 | 0.44 | 0.985 |
| quintiles | square root of persons (central) | 0.178 (0.007) | 3,799 | 39,118 | 0.0971 (0.0066) | 25.8 | 0.69 | 0.09 | 0.999 |
| quintiles | OECD-modified | 0.170 (0.008) | 4,005 | 37,701 | 0.1062 (0.0070) | 25.8 | 0.66 | 0.04 | 1.000 |
| quintiles | per adult (children count zero) | 0.154 (0.009) | 3,950 | 32,224 | 0.1226 (0.0079) | 25.8 | 0.60 | 0.28 | 0.994 |
| quintiles | per capita | 0.141 (0.010) | 3,348 | 24,688 | 0.1356 (0.0085) | 25.7 | 0.55 | 0.51 | 0.980 |
| deciles | none (consumer unit) | 0.202 (0.005) | 4,387 | 61,224 | 0.0717 (0.0050) | 25.9 | 0.78 | 1.03 | 0.922 |
| deciles | square root of persons (central) | 0.178 (0.007) | 3,809 | 39,359 | 0.0968 (0.0062) | 25.8 | 0.69 | 0.96 | 0.932 |
| deciles | OECD-modified | 0.170 (0.007) | 4,006 | 38,005 | 0.1054 (0.0066) | 25.8 | 0.66 | 0.98 | 0.929 |
| deciles | per adult (children count zero) | 0.155 (0.008) | 3,936 | 32,567 | 0.1209 (0.0072) | 25.7 | 0.60 | 1.06 | 0.916 |
| deciles | per capita | 0.143 (0.009) | 3,341 | 24,991 | 0.1337 (0.0080) | 25.7 | 0.55 | 1.25 | 0.885 |

All housing definitions by scale (aH; h / mean expenditure per equivalent adult):

| groups | definition | none | square root of persons | OECD-modified | per adult | per capita |
|---|---|---|---|---|---|---|
| quintiles | shelter | 0.156; 0.0420 | 0.143; 0.0561 | 0.138; 0.0616 | 0.129; 0.0714 | 0.123; 0.0781 |
| quintiles | shelter + utilities | 0.201; 0.0725 | 0.178; 0.0971 | 0.170; 0.1062 | 0.154; 0.1226 | 0.141; 0.1356 |
| quintiles | shelter + utilities ex. telephone | 0.181; 0.0659 | 0.161; 0.0880 | 0.153; 0.0963 | 0.139; 0.1110 | 0.128; 0.1222 |
| quintiles | housing | 0.273; 0.0761 | 0.252; 0.1014 | 0.244; 0.1107 | 0.230; 0.1273 | 0.218; 0.1403 |
| deciles | shelter | 0.156; 0.0422 | 0.143; 0.0562 | 0.138; 0.0611 | 0.131; 0.0692 | 0.124; 0.0760 |
| deciles | shelter + utilities | 0.202; 0.0717 | 0.178; 0.0968 | 0.170; 0.1054 | 0.155; 0.1209 | 0.143; 0.1337 |
| deciles | shelter + utilities ex. telephone | 0.182; 0.0654 | 0.161; 0.0878 | 0.153; 0.0956 | 0.140; 0.1091 | 0.129; 0.1202 |
| deciles | housing | 0.274; 0.0755 | 0.252; 0.1013 | 0.244; 0.1101 | 0.230; 0.1258 | 0.219; 0.1385 |

Fit by equivalence elasticity theta (persons as the size measure; shelter + utilities; aH; h / mean expenditure per equivalent adult; RMSE in pp):

| groups | theta = 0 | 0.25 | 0.5 | 0.75 | 1 | best fit on a 0.01 grid |
|---|---|---|---|---|---|---|
| quintiles | 0.201; 0.072; 0.44 | 0.191; 0.083; 0.28 | 0.178; 0.097; 0.09 | 0.162; 0.114; 0.19 | 0.141; 0.136; 0.51 | theta = 0.58: 0.174; 0.102; 0.06 |
| deciles | 0.202; 0.072; 1.03 | 0.191; 0.083; 0.97 | 0.178; 0.097; 0.96 | 0.162; 0.114; 1.04 | 0.143; 0.134; 1.25 | theta = 0.42: 0.183; 0.092; 0.96 |

Fit of the central equivalized specification (shelter + utilities, quintiles, square-root scale):

| group | equivalized expenditure $ | data share % | fitted share % | residual pp | subsistence h/y % |
|---|---|---|---|---|---|
| Lowest 20 percent | 20,247 | 33.22 | 33.25 | -0.03 | 18.8 |
| Second 20 percent | 27,581 | 29.06 | 29.15 | -0.09 | 13.8 |
| Third 20 percent | 32,716 | 27.50 | 27.37 | 0.13 | 11.6 |
| Fourth 20 percent | 40,595 | 25.61 | 25.52 | 0.08 | 9.4 |
| Highest 20 percent | 66,401 | 22.44 | 22.53 | -0.10 | 5.7 |

**Central choice: the square-root scale, aH = 0.178 and hmin = 0.0971.** (i) It is the OECD's income-distribution convention (OECD-modified, elasticity about 0.53, gives 0.170 and 0.106); (ii) the Stone-Geary curve fits the group means best at an equivalence elasticity of 0.58 for the quintiles (RMSE 0.06 pp, against 0.44 pp with no adjustment and 0.51 pp per capita) and 0.42 for the deciles (two degrees of freedom: suggestive, not an estimate of the scale); (iii) no adjustment and per capita bracket it, and they are the ends of the range together with the other definitions. In dollars a household of mean size has h = $5,945 a year ($495 a month; $3,799 per equivalent adult on mean expenditure per equivalent adult of $39,118), against $4,436 with no adjustment. The marginal share aH does not depend on the unit of the model, only on the adjustment for size: it falls from 0.201 (no adjustment) to 0.141 (per capita). Per adult, the model's literal unit, gives 0.154 and 0.123. The 'mean expenditure per equivalent adult' denominator is aggregate expenditure over aggregate equivalent adults; the CU-weighted mean of the groups' equivalized expenditures would give a slightly larger hmin (0.101 at the square-root scale).

**Cross-check with CE Table 1400 (size of consumer unit).** One-person CUs (expenditures $36,087) spend 31.7% (SE 0.7 pp) of expenditures on shelter + utilities, CUs of two or more people ($71,623, 3.1 persons) 24.6%. The CU-level income curve predicts 29.9% and 25.1%, the square-root curve 26.5% and 25.5%: one-person CUs spend 5.2 pp more than the equivalized curve implies, i.e. economies of scale in housing are stronger than the square root (a second person adds little to subsistence housing). With the central aH the one-person share implies a subsistence share of 0.169 of their expenditure ($6,094 a year) against 0.082 ($5,892) for CUs of two or more. Read literally for a model of single adults the subsistence share would be nearer 0.17 than 0.097, and the aggregate share nearer 0.317 than 0.258. This is a one-point cross-check (one-person CUs are older and more often retired), not a fit; it is reported, not included in the range.

### Terciles of income before taxes (for `housing_share_by_income`)

Terciles are built from the equal-sized groups with fractional weights (from the deciles: tercile 1 = the lowest three deciles plus one third of the fourth, and so on; from the quintiles: tercile 1 = the lowest quintile plus two thirds of the second). The decile-based terciles are closer to true tercile cuts and are the recommended `housing_share_by_income` values. Ratio of means (aggregate dollars over aggregate expenditures, the moment the model must compute) in %, with the parametric-bootstrap SE of the CE sampling error (independent errors, an upper bound) in parentheses; the last three columns are the consumer-unit-weighted means of the groups' ratios (a mean-of-ratios counterpart). These are terciles of consumer-unit income before taxes, so they include household-size differences (people per CU rise from 1.7 to 3.2 across the quintiles), unlike the equivalized fit above:

| groups | definition | tercile 1 (ratio of means) | tercile 2 | tercile 3 | T1 (mean of group shares) | T2 | T3 |
|---|---|---|---|---|---|---|---|
| quintiles | shelter | 22.3 (0.5) | 19.7 (0.3) | 17.9 (0.4) | 22.7 | 19.8 | 18.1 |
| quintiles | utilities | 8.8 (0.2) | 7.5 (0.1) | 5.4 (0.1) | 8.9 | 7.6 | 5.6 |
| quintiles | shelter + utilities | 31.1 (0.6) | 27.2 (0.3) | 23.3 (0.4) | 31.6 | 27.4 | 23.7 |
| quintiles | shelter + utilities ex. telephone | 28.4 (0.5) | 24.6 (0.3) | 21.3 (0.4) | 28.9 | 24.8 | 21.6 |
| quintiles | housing | 37.9 (0.7) | 34.1 (0.5) | 30.5 (0.6) | 38.3 | 34.3 | 30.8 |
| quintiles | mean expenditures $ | 31,824 | 52,849 | 98,990 |  |  |  |
| deciles | shelter | 22.7 (0.5) | 19.6 (0.3) | 17.9 (0.4) | 23.0 | 19.6 | 18.0 |
| deciles | utilities | 8.9 (0.2) | 7.6 (0.1) | 5.4 (0.1) | 8.9 | 7.7 | 5.7 |
| deciles | shelter + utilities | 31.6 (0.6) | 27.2 (0.4) | 23.3 (0.4) | 31.9 | 27.3 | 23.6 |
| deciles | shelter + utilities ex. telephone | 28.9 (0.5) | 24.5 (0.3) | 21.3 (0.4) | 29.2 | 24.6 | 21.5 |
| deciles | housing | 38.2 (0.7) | 34.1 (0.6) | 30.5 (0.5) | 38.6 | 34.2 | 30.7 |
| deciles | mean expenditures $ | 31,263 | 52,393 | 100,005 |  |  |  |

### Cross-check: ACS national housing cost burdens by household income

ACS 2018 1-year, United States (tables B25074 renters' gross rent and B25095 owners' selected monthly owner costs as a percentage of household income, households with a computed burden, B25106 for the coarse 30%-or-more shares; mean income by bin from B19001 midpoints and, for $100,000 and over, from aggregate income B19025). The 5-year 2014-18 file is in the JSON and agrees to about 0.01. Mean burden uses class midpoints, with 12% for 'less than 20%' and 70% for '50% or more'; the median is interpolated inside classes (shown when it lies in a bounded class). ACS costs include utilities and, for owners, mortgage principal, taxes and insurance (not maintenance); CE shelter excludes principal and CE utilities include telephone service. Burdens are household-level ratios to household income with no adjustment for household size, so the comparison is with the consumer-unit-level (no adjustment) CE fit.

| household income | households (M) | mean income $ | burden >= 30% (%) | burden >= 50% (%) | mean burden (%) | median burden (%) |
|---|---|---|---|---|---|---|
| $0 - $10,000 | 5.60 | 5,000 | 92.2 | 81.4 | 62.9 |  |
| $10,000 - $20,000 | 10.14 | 14,994 | 78.3 | 54.4 | 51.9 |  |
| $20,000 - $35,000 | 16.06 | 27,460 | 62.8 | 27.2 | 39.9 | 36.3 |
| $35,000 - $50,000 | 14.74 | 42,350 | 40.8 | 10.2 | 29.6 | 26.8 |
| $50,000 - $75,000 | 20.79 | 62,132 | 23.5 | 4.0 | 23.2 | 21.3 |
| $75,000 - $100,000 | 15.10 | 87,500 | 12.3 | 1.4 | 19.0 |  |
| $100,000 and over | 35.36 | 189,795 | 4.6 | 0.4 | 15.4 |  |

B25106 (coarse bins): share paying 30% or more of income: Less than $20,000 83.2% (90% margin of error +-0.13 pp); $20,000-34,999 62.8% (90% margin of error +-0.21 pp); $35,000-49,999 40.8% (90% margin of error +-0.22 pp); $50,000-74,999 23.5% (90% margin of error +-0.14 pp); $75,000 or more 6.9% (90% margin of error +-0.05 pp). All households with a computed burden: 31.9% (B25106: 31.9%, margin of error +-0.06 pp): national ACS sampling error is negligible here, the uncertainty is the class values for the open classes.

Stone-Geary form on income, s(y) = a + b / y, fitted to the mean burdens (weights = households): a = 0.196, h = $3,488 = 0.039 of mean household income ($89,731), RMSE 7.1 pp (the lowest bins are fitted poorly); a ranges 0.183-0.220 and h/mean income 0.033-0.044 when the class values are varied (less than 20%: 10-16%; 50% or more: 60-80%). In dollars h is 0.057 of CE mean expenditures, against 0.072 from the CE fit. CE shelter + utilities over income before taxes, by quintile: 77.7 / 37.2 / 25.9 / 19.6 / 13.0% (mean incomes 11,285 / 31,237 / 54,900 / 90,478 / 204,975), against ACS mean burdens of 39.9 / 29.6 / 23.2 / 19.0 / 15.4% at mean incomes 27,460 / 42,350 / 62,132 / 87,500 / 189,795 for the bins from $20,000 up.

## (c) Housing supply elasticities

### Saiz (2010)

A. Saiz, 'The Geographic Determinants of Housing Supply', *Quarterly Journal of Economics* 125(3): 1253-1296. Price elasticity of the supply of housing units by metropolitan area, estimated from 1970-2000 changes (a 30-year horizon) with instrumented demand shocks interacted with the share of land unavailable for development (geography) and the Wharton regulation index (Table V column 2; estimates in Table VI, pp. 1282-1284, for the 95 MSAs with population over 500,000 in 2000). Text (p. 1281): population-weighted average 1.75, unweighted 2.5. The author's file (269 MSAs; downloaded from the link on urbaneconomics.mit.edu) reproduces both and Table VI:

| sample | MSAs | mean | population-weighted mean | median | p25 - p75 | p10 - p90 | min - max |
|---|---|---|---|---|---|---|---|
| 95 MSAs over 500,000 (Table VI) | 95 | 1.87 | 1.55 | 1.65 | 1.09 - 2.37 | 0.78 - 3.09 | 0.60 - 5.45 |
| all 269 MSAs | 269 | 2.54 | 1.76 | 2.26 | 1.56 - 3.19 | 1.07 - 4.37 | 0.60 - 12.15 |

Selected MSAs (Table VI): San Francisco 0.66, Oakland 0.70, San Jose 0.76, New York 0.76, Los Angeles 0.63, Chicago 0.81, Boston 0.86, Seattle 0.88, Washington 1.61, Atlanta 2.55, Dallas 2.18, Houston 2.30, Phoenix 1.61, Minneapolis 1.45, Indianapolis 4.00. Land-constrained coastal metros are below 1; large interior metros are above 2. These are metro-wide averages of housing units, with no separation of the city from the suburbs.

### Baum-Snow and Han (2024)

N. Baum-Snow and L. Han, 'The Microgeography of Housing Supply', *Journal of Political Economy* 132(6): 1897-1946. Tract-level supply elasticities for 2000-2010 changes (a ten-year horizon; the authors attribute most of the gap to Saiz's much larger 1970-2000 elasticities to the later, shorter period, introduction p. 4), identified from Bartik-type labor demand shocks to commuting destinations; the tract file below is the authors' predicted-elasticity release of September 2023 (50,410 tracts in 306 metro regions; 63,897 tract-region rows, `wgt` = 1 / number of regions of a tract). Text (accepted manuscript, September 5, 2023): in the estimation sample the average elasticities are 0.35 (housing units), 0.42 (floorspace) and 0.09 (land development); predicted for all 50,410 urban tracts, unit and floorspace elasticities average 0.29 and 0.51 (quartiles 0.14-0.44 and 0.33-0.70), with more variation within than between metro regions (introduction, p. 2); at the CBD the implied unit elasticity is 0.11, rising to 0.59 near halfway to the region edge (Table 5 column 1; section 5.2, p. 25); supply elasticities increase with CBD distance, flatten in the suburbs (about 30 to 90 percent of the way out) and rise again at the fringe because the developed fraction falls with distance while regulation is tighter in the suburbs than in central cities (Figure 3 discussion, p. 30). The file reproduces the text (FMM-IV linear specification, tract-weighted): units 0.289 (quartiles 0.14-0.44), floorspace 0.511 (0.33-0.70). Aggregated to the metro-region level the FMM-IV unit elasticities average 0.38-0.41 (s.d. 0.11; main text p. 33 and Online Appendix D), against 2.54 for Saiz's 1970-2000 estimates for the same metros; the authors attribute the gap to the later and shorter period, to aggregation (their aggregation weights tracts with large initial stocks, which respond less, more heavily) and to differences in the price and quantity measures.

| estimator (2001 developed fraction) | units | floorspace | land development | units, housing-unit weighted | floorspace, housing-unit weighted |
|---|---|---|---|---|---|
| FMM lin | 0.289 | 0.511 | 0.010 | 0.301 | 0.525 |
| FMM quad | 0.339 | 0.601 | -0.017 | 0.346 | 0.609 |
| IV lin | 0.367 | 0.447 | 0.083 | 0.381 | 0.459 |
| IV quad | 0.400 | 0.481 | 0.053 | 0.405 | 0.484 |

**By distance to the CBD** (fraction of the way from the CBD to the region edge; housing-unit weighted):

| distance | housing units % | developed fraction | FMM lin units | FMM lin space | FMM quad units | FMM quad space | IV lin units | IV lin space |
|---|---|---|---|---|---|---|---|---|
| 0.00-0.05 | 8.7 | 0.57 | 0.13 | 0.31 | 0.22 | 0.46 | 0.16 | 0.25 |
| 0.05-0.10 | 14.6 | 0.48 | 0.21 | 0.41 | 0.26 | 0.49 | 0.29 | 0.36 |
| 0.10-0.20 | 23.9 | 0.38 | 0.29 | 0.51 | 0.31 | 0.55 | 0.40 | 0.47 |
| 0.20-0.30 | 15.9 | 0.29 | 0.35 | 0.59 | 0.37 | 0.64 | 0.47 | 0.54 |
| 0.30-0.50 | 19.8 | 0.24 | 0.36 | 0.61 | 0.41 | 0.69 | 0.46 | 0.53 |
| 0.50-0.70 | 11.2 | 0.22 | 0.34 | 0.58 | 0.42 | 0.72 | 0.38 | 0.47 |
| 0.70-1.00 | 5.9 | 0.12 | 0.38 | 0.64 | 0.48 | 0.80 | 0.38 | 0.49 |

**By developed fraction of the tract in 2001:**

| developed fraction | housing units % | mean CBD distance | FMM lin units | FMM lin space | FMM quad units | FMM quad space |
|---|---|---|---|---|---|---|
| 0.00-0.10 | 19.4 | 0.49 | 0.54 | 0.80 | 0.76 | 1.10 |
| 0.10-0.30 | 26.5 | 0.30 | 0.41 | 0.66 | 0.42 | 0.69 |
| 0.30-0.60 | 40.1 | 0.21 | 0.22 | 0.44 | 0.16 | 0.39 |
| 0.60-1.00 | 14.0 | 0.14 | -0.00 | 0.13 | 0.16 | 0.40 |

**Central city versus suburb.** Within each region the tracts closest to the CBD holding 42.1% of occupied units (the city share of partitioned pupils in the locale partition) are the 'city' and the rest the 'suburb'; means are housing-unit weighted pooled over regions. The classification is by distance, not by the principal-city boundary, so it is a proxy for the locale split.

| estimator | quantity | city (10-year) | suburb (10-year) | suburb / city | within-region ratio, median / housing-weighted median |
|---|---|---|---|---|---|
| FMM lin | units | 0.236 | 0.347 | 1.47 | 1.35 / 1.53 |
| FMM lin | space | 0.439 | 0.588 | 1.34 | 1.27 / 1.36 |
| FMM quad | units | 0.278 | 0.396 | 1.43 | 1.60 / 1.37 |
| FMM quad | space | 0.515 | 0.677 | 1.31 | 1.48 / 1.29 |
| IV lin | units | 0.314 | 0.429 | 1.37 |  |
| IV lin | space | 0.390 | 0.510 | 1.31 |  |
| IV quad | units | 0.341 | 0.451 | 1.32 |  |
| IV quad | space | 0.417 | 0.533 | 1.28 |  |

Mean developed fraction in 2001: city 0.45, suburb 0.26. Ratio of suburban to central-city elasticity (FMM linear, units) at other city shares: 0.25: 1.62, 0.33: 1.54, 0.421: 1.47, 0.5: 1.43.

### Recommended long-run eps_H

**Convention.** The model's housing market is a steady state, so the relevant elasticity is that of the housing stock after full adjustment. For a stock with annual depreciation delta, a response measured over T years is the fraction 1 - (1 - delta)^T of the long-run response: eps_LR = eps_T / (1 - (1 - delta)^T) (N. Baum-Snow and G. Duranton, 'Housing Supply and Housing Affordability', chapter for the Handbook of Regional and Urban Economics vol. 6, slides of the Urban Economics Association summer school 2025, slide 76 of 90). At delta = 3.5% the long run is 3.34 times a 10-year and 1.52 times a 30-year elasticity (5.47 and 2.20 at 2%; 7.13 and 2.74 at 1.5%). **One convention is applied to every anchor**: Saiz (2010) is a 30-year estimate (1970-2000 changes), Baum-Snow and Han (2024) and the Saiz-method replication of their Online Appendix D are 10-year estimates (2000-10 changes). The slide adds that supply elasticities have fallen 'in US supply-constrained markets only' when Baum-Snow and Han are compared with Saiz. The slide and the Appendix D statements quoted below were checked against the text of the cached files (slides: True; manuscript: True).

**Period check (Online Appendix D, pp. 17-18, and Table A8, p. 26, of the accepted manuscript).** Baum-Snow and Han replicate Saiz's estimation for 2000-2010 'using the same set of metros, specifications and variables': 'The average implied metro supply elasticity declines to 1.25 for 2000-2010, indicating that many metropolitan areas have become increasingly supply constrained over time', against a mean of 2.54 (s.d. 1.43) for Saiz's 1970-2000 estimates (Table A8 reports means of 0.80-1.30 across its 2000-10 specifications and 1.54-2.56 for the 1970-2000 ones). The ratio is 0.492. The horizon alone predicts (1 - (1 - delta)^10) / (1 - (1 - delta)^30) = 0.456 at delta = 3.5%, i.e. a 2000-10 mean of 1.16, and a depreciation rate of 4.5% would match the ratio exactly. So the lower 2000-10 estimates are mostly the shorter horizon, and no additional period adjustment is made to the central anchor; the 2000-10 vintage converted with the same factor gives 2.55 for the population-weighted 95-MSA mean (Saiz's mean x 0.49 x 3.34), close to the central anchor. Baum-Snow and Duranton's reading, that the fall is concentrated in supply-constrained markets, and the post-2010 zip-code evidence of Gorback and Keys (2024, 2010-20; slide 74: 'Find even smaller supply elasticities'), are reasons for a lower value in constrained coastal CZs, which the model's CZ-average locations do not distinguish; they are covered by the low anchor, and the unconverted 10-year values (0.30-0.76) are not long-run elasticities but are covered by the grid's 0.5.

| metro anchor (housing units) | horizon (years) | as published | long run (delta = 3.5%) | role |
|---|---|---|---|---|
| Saiz (2010), population-weighted mean over the 95 MSAs of Table VI | 30 | 1.55 | 2.37 | central |
| Saiz (2010), population-weighted mean over all 269 MSAs | 30 | 1.76 | 2.68 |  |
| Saiz (2010), unweighted mean over all 269 MSAs | 30 | 2.54 | 3.86 | high |
| Saiz (2010), 25th percentile of the 95 MSAs (dispersion across metros) | 30 | 1.09 | 1.67 | dispersion |
| Saiz (2010), 75th percentile of the 95 MSAs (dispersion across metros) | 30 | 2.37 | 3.61 | dispersion |
| Saiz's method on 2000-10 data (Baum-Snow and Han, Online Appendix D), mean over metros | 10 | 1.25 (0.80-1.30 across Table A8 specifications) | 4.17 | period check |
| Saiz's population-weighted mean x (1.25 / 2.54), 10-year equivalent | 10 | 0.76 | 2.55 | period check |
| Baum-Snow and Han, metro-region aggregates of the tract elasticities (units, FMM-IV, mean) | 10 | 0.38-0.41 | 1.27-1.37 | reference |
| Baum-Snow and Han, tract elasticities, housing-unit-weighted mean (units, FMM linear) | 10 | 0.30 | 1.00 | low |

Construction: eps_H(city) and eps_H(suburb) follow from the metro anchor and the within-metro ratio suburb / city of Baum-Snow and Han's predicted tract elasticities, pooled housing-unit-weighted means on the city / suburb split matched to the pupil share of the locale partition (city share s = 0.421); with equal price changes in the two places the metro elasticity is s eps_c + (1 - s) eps_s (perfect substitutes, slide 77), so eps_c = anchor / (s + (1 - s) ratio) and eps_s = ratio x eps_c. The central ratio 1.40 is the mean of the four estimators for units; the range 1.25-1.73 covers all 64 variants: units and floorspace, the four estimators, city shares 0.25-0.5, and the within-region medians (matched split only: 1.28-1.47; the extremes are the within-region median of FMM lin space at s = 0.5 and the within-region median of FMM quad units at s = 0.25). Mean ratio for units by city share: 0.25: 1.49, 0.33: 1.44, 0.421: 1.40, 0.5: 1.37.

|  | central | low - high | ten-year Baum-Snow-Han value (units; floorspace) |
|---|---|---|---|
| eps_H city | 1.93 | 0.71 - 3.37 | 0.24; 0.44 |
| eps_H suburb | 2.69 | 1.10 - 4.70 | 0.35; 0.59 |
| metro average | 2.37 | 1.00 - 3.86 |  |
| memo: no horizon conversion (Saiz's 30-year estimate read as long run: the first version of this report) | 1.26 city / 1.76 suburb |  |  |

Low - high is the extreme over the three metro anchors (low: Baum-Snow and Han tract elasticities, housing-unit weighted, x 3.34; high: Saiz's unweighted mean over 269 MSAs x 1.52) and the lowest and highest suburb / city ratio. **Why the central anchor is Saiz's population-weighted mean:** the model's CZs are pupil-weighted and mostly large metros, which the population weights approximate; it is the long-run object under the one convention; and the period check shows that the 2000-10 literature does not call for a lower long-run value once horizons are aligned. Baum-Snow and Han's tract elasticities, a within-metro identification closer to the model's city-versus-suburb margin, give a much lower value (low anchor; 1.00 for the metro, 0.71 for the city at the extremes), and they are a floor for sensitivity runs rather than the central case because the tract estimates are within-region elasticities identified from relative changes and aggregate to metro values (0.38-0.41 x 3.34 = 1.27-1.37) well below Saiz's. Sensitivity to the depreciation rate (central anchor with each delta, the same ratio):

| depreciation delta | factor, 30 years | factor, 10 years | central metro anchor | eps_H city | eps_H suburb | Baum-Snow-Han tract anchor |
|---|---|---|---|---|---|---|
| 1.5% | 2.74 | 7.13 | 4.26 | 3.47 | 4.84 | 2.14 |
| 2.0% | 2.20 | 5.47 | 3.42 | 2.78 | 3.88 | 1.64 |
| 3.5% | 1.52 | 3.34 | 2.37 | 1.93 | 2.69 | 1.00 |

Suggested sensitivity grid for the model: eps_H in {0.5, 1, 2, 4} in both places and the elastic limit (which nests the fixed-cost model). Floorspace elasticities are about 1.8 times the unit ones in the Baum-Snow-Han data; if the model's H is read as floorspace or quality-adjusted housing services the values scale up by that factor. The recommended values are constructed from the literature (a metro anchor split by a within-metro ratio), not estimated for the model's two locations.

## (d) Transport: the suburban car requirement

ACS 5-year 2014-18 school-district tables, same panel, locale partition and draws as (a): B25044 tenure by vehicles available (occupied units), B25046 aggregate vehicles available (vehicles per household = B25046_001 / B25044_001), B08201 household size by vehicles available (a second zero-vehicle measure, one-person households and the household-size classes used below). Pooled over the districts of a location within each CZ, suburb minus city, CZ-weighted mean, whole-CZ bootstrap.

| statistic | city | suburb | gap (SE) | 95% interval | CZs |
|---|---|---|---|---|---|
| households with no vehicle available, share of occupied units (B25044) | 0.133 | 0.053 | -0.0795 (0.0103) | [-0.101, -0.061] | 188 |
| vehicles per occupied unit: aggregate vehicles (B25046) / occupied units (B25044) | 1.562 | 1.940 | 0.3774 (0.0297) | [0.322, 0.439] | 188 |
| log ratio, vehicles per occupied unit (suburb / city) | 1.562 | 1.940 | 0.2304 (0.0231) | [0.188, 0.279] | 188 |
| owner-occupied units with no vehicle (B25044) | 0.045 | 0.022 | -0.0232 (0.0050) | [-0.034, -0.014] | 188 |
| renter-occupied units with no vehicle (B25044) | 0.217 | 0.124 | -0.0934 (0.0115) | [-0.117, -0.073] | 188 |
| vehicles per owner-occupied unit (B25046_002 / B25044_002) | 1.917 | 2.164 | 0.2466 (0.0292) | [0.192, 0.306] | 188 |
| vehicles per renter-occupied unit (B25046_003 / B25044_009) | 1.199 | 1.439 | 0.2405 (0.0208) | [0.200, 0.283] | 188 |
| renter share of occupied units (B25044_009 / B25044_001) | 0.491 | 0.311 | -0.1804 (0.0095) | [-0.199, -0.162] | 188 |
| units with two or more vehicles (B25044), share | 0.482 | 0.649 | 0.1668 (0.0132) | [0.142, 0.194] | 188 |
| households with no vehicle (B08201_002 / B08201_001) | 0.133 | 0.053 | -0.0795 (0.0103) | [-0.101, -0.061] | 188 |
| one-person households with no vehicle (B08201_008 / B08201_007) | 0.238 | 0.137 | -0.1015 (0.0129) | [-0.128, -0.078] | 188 |

### Within-CZ district regressions

District outcome on a suburb indicator, CZ fixed effects (weighted within-transformation), weights = occupied units (households) or pupils, log median household income of the district (`median_hh_income`, ACS B19013) as the income control, and (a further check) the district's renter share, which is a choice correlated with location, or the district's shares of 1-, 3- and 4-or-more-person households (the 2-person share is the omitted category). SE (cluster) = CZ-clustered; SE (boot) = whole-CZ bootstrap with the same draws as everywhere else. 'log income quadratic' adds the squared deviation of log income from its mean.

| outcome | weights | controls | suburb (SE cluster; SE boot) | other coefficients (SE cluster; SE boot) | districts | CZs | within R2 |
|---|---|---|---|---|---|---|---|
| zero-vehicle share | households | no controls | -0.0786 (0.0111; 0.0104) |  | 3200 | 188 | 0.337 |
| zero-vehicle share | households | log income | -0.0492 (0.0087; 0.0083) | log income -0.1176 (0.0114; 0.0108) | 3200 | 188 | 0.554 |
| zero-vehicle share | households | log income quadratic | -0.0472 (0.0091; 0.0086) | log income -0.1160 (0.0116; 0.0111); log income squared (centered) 0.0388 (0.0149; 0.0146) | 3200 | 188 | 0.561 |
| zero-vehicle share | households | log income renter share | -0.0197 (0.0059; 0.0056) | log income -0.0330 (0.0180; 0.0170); renter share 0.2915 (0.0484; 0.0460) | 3200 | 188 | 0.672 |
| zero-vehicle share | households | log income size shares | -0.0268 (0.0066; 0.0063) | log income -0.0639 (0.0175; 0.0164); share of 1-person households 0.6529 (0.0939; 0.0868); share of 3-person households 0.2511 (0.0938; 0.0920); share of 4+-person households 0.2575 (0.0459; 0.0423) | 3200 | 188 | 0.627 |
| zero-vehicle share | pupils | no controls | -0.0682 (0.0100; 0.0095) |  | 3200 | 188 | 0.293 |
| zero-vehicle share | pupils | log income | -0.0393 (0.0071; 0.0068) | log income -0.1152 (0.0121; 0.0112) | 3200 | 188 | 0.562 |
| zero-vehicle share | pupils | log income quadratic | -0.0368 (0.0074; 0.0070) | log income -0.1151 (0.0127; 0.0120); log income squared (centered) 0.0444 (0.0117; 0.0116) | 3200 | 188 | 0.575 |
| zero-vehicle share | pupils | log income renter share | -0.0143 (0.0051; 0.0049) | log income -0.0381 (0.0131; 0.0124); renter share 0.2621 (0.0405; 0.0385) | 3200 | 188 | 0.670 |
| zero-vehicle share | pupils | log income size shares | -0.0203 (0.0052; 0.0051) | log income -0.0677 (0.0180; 0.0167); share of 1-person households 0.6021 (0.0891; 0.0813); share of 3-person households 0.2744 (0.0717; 0.0710); share of 4+-person households 0.2429 (0.0426; 0.0394) | 3200 | 188 | 0.628 |
| vehicles per household | households | no controls | 0.3686 (0.0311; 0.0292) |  | 3200 | 188 | 0.392 |
| vehicles per household | households | log income | 0.2454 (0.0262; 0.0250) | log income 0.4915 (0.0386; 0.0370) | 3200 | 188 | 0.594 |
| vehicles per household | households | log income quadratic | 0.2401 (0.0272; 0.0259) | log income 0.4872 (0.0392; 0.0376); log income squared (centered) -0.1065 (0.0587; 0.0560) | 3200 | 188 | 0.597 |
| vehicles per household | households | log income renter share | 0.0772 (0.0160; 0.0153) | log income 0.0087 (0.0599; 0.0568); renter share -1.6641 (0.1354; 0.1301) | 3200 | 188 | 0.798 |
| vehicles per household | households | log income size shares | 0.0847 (0.0149; 0.0143) | log income 0.1780 (0.0589; 0.0563); share of 1-person households -3.6423 (0.3341; 0.3154); share of 3-person households -1.5788 (0.2970; 0.2947); share of 4+-person households -0.5504 (0.1833; 0.1774) | 3200 | 188 | 0.810 |
| vehicles per household | pupils | no controls | 0.3441 (0.0288; 0.0274) |  | 3200 | 188 | 0.360 |
| vehicles per household | pupils | log income | 0.2215 (0.0223; 0.0213) | log income 0.4879 (0.0367; 0.0344) | 3200 | 188 | 0.594 |
| vehicles per household | pupils | log income quadratic | 0.2142 (0.0232; 0.0221) | log income 0.4876 (0.0377; 0.0355); log income squared (centered) -0.1318 (0.0449; 0.0435) | 3200 | 188 | 0.599 |
| vehicles per household | pupils | log income renter share | 0.0672 (0.0140; 0.0136) | log income 0.0119 (0.0469; 0.0449); renter share -1.6183 (0.1066; 0.1024) | 3200 | 188 | 0.792 |
| vehicles per household | pupils | log income size shares | 0.0700 (0.0123; 0.0121) | log income 0.1873 (0.0584; 0.0555); share of 1-person households -3.5912 (0.3249; 0.3055); share of 3-person households -1.6799 (0.2467; 0.2455); share of 4+-person households -0.5758 (0.2052; 0.1955) | 3200 | 188 | 0.806 |

Income explains 33% of the within-CZ suburb-city difference in vehicles per household and 37% of the difference in the zero-vehicle share; adding the renter share removes most of the rest, but tenure is itself a location-related choice (within tenure, the pooled gaps are +0.25 vehicles for owners and +0.24 for renters). Adding the district household-size shares cuts the suburb coefficient on vehicles per household from 0.245 to 0.085: most of what remains after income is household size.

### Household size: the vehicles gap holding household size fixed

Suburban households are larger and larger households own more vehicles, while the model's household is one adult, so the vehicles gap per household is separated from household composition with B08201 (household size by vehicles available), the school-district table of 1-, 2-, 3- and 4-or-more-person households by vehicles (0, 1, 2, 3, 4+). Vehicles per household within a size class use the midpoints 0, 1, 2, 3 and, for '4 or more', a mean calibrated by location to the aggregate B25046 / B25044 (4.37 in the city, 4.36 in the suburb); the class-based total reproduces the B25046-based one (largest difference over CZ x location cells 0.034; pooled raw gap 0.377 against 0.377). The 4-or-more-person class is open-ended, so larger suburban families inside it are not removed. Pooled over the districts of a location within each CZ, CZ-weighted, whole-CZ bootstrap:

| household size | share of households, city % | suburb % | gap pp (SE) | vehicles per household, city | suburb | gap (SE) | gap, income-adjusted (SE) | zero-vehicle share gap pp (SE) | zero-vehicle gap, income-adjusted pp (SE) |
|---|---|---|---|---|---|---|---|---|---|
| 1 person | 31.8 | 23.6 | -8.2 (0.6) | 0.90 | 1.07 | 0.169 (0.015) | 0.130 (0.017) | -10.2 (1.3) | -7.0 (1.3) |
| 2 persons | 31.4 | 32.7 | 1.3 (0.3) | 1.65 | 1.93 | 0.285 (0.031) | 0.167 (0.025) | -6.0 (1.0) | -3.4 (0.7) |
| 3 persons | 15.2 | 16.8 | 1.7 (0.2) | 1.90 | 2.26 | 0.363 (0.037) | 0.193 (0.026) | -5.6 (0.8) | -3.0 (0.5) |
| 4 or more persons | 21.6 | 26.9 | 5.3 (0.5) | 2.16 | 2.51 | 0.350 (0.034) | 0.213 (0.023) | -4.6 (0.7) | -2.3 (0.4) |

| standardized suburb-minus-city gap | vehicles per household (SE) | zero-vehicle share, pp (SE) |
|---|---|---|
| raw gap, all households (B25046 / B25044) | 0.377 (0.030) |  |
| size-standardized gap, CZ-pooled size distribution (within-class gaps averaged with each CZ's pooled household-size shares) | 0.281 (0.027) | -6.7 (0.9) |
|   with the city's own household-size shares | 0.269 (0.025) |  |
|   with the suburb's own household-size shares | 0.286 (0.028) |  |
| size-standardized regression gap, no income control (district regressions within class, national class shares 27 / 32 / 16 / 24%) | 0.272 (0.027) | -6.7 (0.9) |
| **size-standardized and income-adjusted regression gap** (log median household income control) | 0.172 (0.022) | -4.1 (0.7) |

### Cost of a vehicle and the fixed goods cost

AAA, *Your Driving Costs*, 2018 edition (composite average vehicle, nine categories weighted by sales; new vehicle; brochure p. 9; every number below was found in the cached brochure text: True): ownership ('fixed') costs $5,960 a year = insurance $1,189 + license, registration and taxes $738 + depreciation $3,289 + finance $744; operating costs 19.26 cents a mile (fuel 11.05, maintenance, repair and tires 8.21); total per year $7,531 at 10,000 miles, **$8,849 at 15,000** (58.99 cents a mile), $10,213 at 20,000. These are costs of a new vehicle held five years, so they overstate what the typical household pays.

CE 2018 (quintiles of income before taxes) transportation lines and vehicles per consumer unit; private-vehicle spending = vehicle purchases (net outlay) + gasoline, other fuels and motor oil + other vehicle expenses (finance charges, maintenance and repairs, insurance, rentals, leases, licenses):

| group | vehicles per CU | with a vehicle % | transportation $ | transportation % of expenditures | private-vehicle spending $ | per vehicle $ |
|---|---|---|---|---|---|---|
| Lowest 20 percent | 1.0 | 66 | 3,718 | 14.1 | 3,464 | 3,464 |
| Second 20 percent | 1.6 | 89 | 6,761 | 16.9 | 6,322 | 3,951 |
| Third 20 percent | 1.9 | 94 | 8,636 | 16.7 | 8,087 | 4,256 |
| Fourth 20 percent | 2.3 | 95 | 11,288 | 16.3 | 10,418 | 4,530 |
| Highest 20 percent | 2.7 | 97 | 18,387 | 15.5 | 16,411 | 6,078 |

Across the five groups the elasticity of the vehicle count with respect to total expenditures is 0.65 (deciles 0.66) and that of private-vehicle spending 1.01: spending on vehicles rises one for one with total expenditures, so a car is only partly a fixed cost. Cost per vehicle in the CE (all consumer units): private-vehicle spending $4,707; ownership-like items (purchases, finance, insurance, leases and licenses) $3,129. **The CE publishes vehicles per consumer unit to one decimal** (1.9 overall; 1.0-2.7 by quintile), so every CE cost per vehicle carries a rounding error of up to +-2.6% (all-CU private-vehicle spending $4,586-$4,834), and the vehicle-count elasticity above is computed from the rounded counts.

**Conversion.** r_car = (extra vehicles per household) x (annual cost of a vehicle) / (mean expenditures per consumer unit, $61,224.13), a fixed goods cost of living in the suburb in units of mean consumption. The ratio of a per-household cost to per-household consumption is unit-free if both scale with household size (`housing_costs.md`), so the extra vehicles are standardized for household size; for one-person households the denominator is the mean expenditure of one-person consumer units ($36,087, CE Table 1400). Extra vehicles per household, from the unadjusted to the composition-adjusted: raw pooled 0.377 (0.030); income-adjusted regression, not adjusted for household size (the first version) 0.245 (0.025); size-standardized 0.281 (0.027); **size-standardized and income-adjusted 0.172 (0.022) (central)**; one-person households 0.169 (0.015) raw and 0.130 (0.017) income-adjusted; income and household-size shares as controls 0.085 (0.014) (an over-control, since the size mix is itself a location outcome: the low end). On the extensive margin the fall in the zero-vehicle share is 0.080 (0.010) raw, 0.049 (0.008) income-adjusted and 0.041 (0.007) size-standardized and income-adjusted, times the cost of one vehicle.

| extra | value | cost per vehicle-year | $ | annual $ per household | denominator $ | r_car (SE from the extra-vehicle SE only) |
|---|---|---|---|---|---|---|
| extra vehicles, raw (not adjusted) | 0.377 | CE private-vehicle spending per vehicle | 4,707 | 1,776 | 61,224 | 0.0290 (0.0023) |
| extra vehicles, raw (not adjusted) | 0.377 | CE ownership-like spending per vehicle | 3,129 | 1,181 | 61,224 | 0.0193 (0.0015) |
| extra vehicles, raw (not adjusted) | 0.377 | AAA total, 15,000 miles | 8,849 | 3,340 | 61,224 | 0.0545 (0.0043) |
| extra vehicles, income-adjusted (not size-adjusted: first version) | 0.245 | CE private-vehicle spending per vehicle | 4,707 | 1,155 | 61,224 | 0.0189 (0.0019) |
| extra vehicles, income-adjusted (not size-adjusted: first version) | 0.245 | CE ownership-like spending per vehicle | 3,129 | 768 | 61,224 | 0.0125 (0.0013) |
| extra vehicles, income-adjusted (not size-adjusted: first version) | 0.245 | AAA total, 15,000 miles | 8,849 | 2,172 | 61,224 | 0.0355 (0.0036) |
| extra vehicles, size-standardized | 0.281 | CE private-vehicle spending per vehicle | 4,707 | 1,323 | 61,224 | 0.0216 (0.0021) |
| extra vehicles, size-standardized | 0.281 | CE ownership-like spending per vehicle | 3,129 | 879 | 61,224 | 0.0144 (0.0014) |
| extra vehicles, size-standardized | 0.281 | AAA total, 15,000 miles | 8,849 | 2,487 | 61,224 | 0.0406 (0.0039) |
| extra vehicles, size-standardized and income-adjusted | 0.172 | CE private-vehicle spending per vehicle | 4,707 | 811 | 61,224 | 0.0132 (0.0017) |
| extra vehicles, size-standardized and income-adjusted | 0.172 | CE ownership-like spending per vehicle | 3,129 | 539 | 61,224 | 0.0088 (0.0011) |
| extra vehicles, size-standardized and income-adjusted | 0.172 | AAA total, 15,000 miles | 8,849 | 1,524 | 61,224 | 0.0249 (0.0031) |
| extra vehicles, income and household-size shares as controls | 0.085 | CE private-vehicle spending per vehicle | 4,707 | 399 | 61,224 | 0.0065 (0.0011) |
| extra vehicles, income and household-size shares as controls | 0.085 | CE ownership-like spending per vehicle | 3,129 | 265 | 61,224 | 0.0043 (0.0007) |
| extra vehicles, income and household-size shares as controls | 0.085 | AAA total, 15,000 miles | 8,849 | 749 | 61,224 | 0.0122 (0.0021) |
| fewer zero-vehicle households, income-adjusted (not size-adjusted) | 0.049 | CE private-vehicle spending per vehicle | 4,707 | 231 | 61,224 | 0.0038 (0.0006) |
| fewer zero-vehicle households, income-adjusted (not size-adjusted) | 0.049 | CE ownership-like spending per vehicle | 3,129 | 154 | 61,224 | 0.0025 (0.0004) |
| fewer zero-vehicle households, income-adjusted (not size-adjusted) | 0.049 | AAA total, 15,000 miles | 8,849 | 435 | 61,224 | 0.0071 (0.0012) |
| fewer zero-vehicle households, size-standardized and income-adjusted | 0.041 | CE private-vehicle spending per vehicle | 4,707 | 191 | 61,224 | 0.0031 (0.0006) |
| fewer zero-vehicle households, size-standardized and income-adjusted | 0.041 | CE ownership-like spending per vehicle | 3,129 | 127 | 61,224 | 0.0021 (0.0004) |
| fewer zero-vehicle households, size-standardized and income-adjusted | 0.041 | AAA total, 15,000 miles | 8,849 | 359 | 61,224 | 0.0059 (0.0011) |
| extra vehicles, one-person households, raw | 0.169 | CE private-vehicle spending per vehicle | 4,707 | 797 | 36,087 | 0.0221 (0.0020) |
| extra vehicles, one-person households, raw | 0.169 | CE ownership-like spending per vehicle | 3,129 | 530 | 36,087 | 0.0147 (0.0013) |
| extra vehicles, one-person households, raw | 0.169 | AAA total, 15,000 miles | 8,849 | 1,499 | 36,087 | 0.0415 (0.0037) |
| extra vehicles, one-person households, income-adjusted | 0.130 | CE private-vehicle spending per vehicle | 4,707 | 611 | 36,087 | 0.0169 (0.0022) |
| extra vehicles, one-person households, income-adjusted | 0.130 | CE ownership-like spending per vehicle | 3,129 | 406 | 36,087 | 0.0113 (0.0015) |
| extra vehicles, one-person households, income-adjusted | 0.130 | AAA total, 15,000 miles | 8,849 | 1,149 | 36,087 | 0.0318 (0.0041) |

(The JSON has all five cost measures, including AAA ownership costs and the AAA total at 10,000 miles.)

**Central: r_car = 0.0132 (SE 0.0017)** = size-standardized, income-adjusted extra vehicles x CE private-vehicle spending per vehicle ($811 a year), against 0.0189 in the first version, which adjusted for income but not for household size (a 30% reduction). Range 0.0065 (income and household-size shares as controls, CE private-vehicle cost) to 0.0406 (size-standardized, no income adjustment, AAA total cost at 15,000 miles; the unadjusted raw gap with the AAA cost gives 0.0545, the first version's upper end). **Single-adult reading** (one-person households, income-adjusted, divided by one-person consumer units' expenditures): 0.0169 (0.0022); it is higher than the central value because vehicles scale less than proportionally with household expenditure (economies of scale in car use). Extensive margin only (size-standardized, income-adjusted): 0.0031. For scale, the housing-cost gap of `housing_costs.md` is r_p = 0.0276 (0.0044), so the car requirement adds about 48% to the residence charge if the whole housing gap were also a fixed cost (the same household-size caveat applies to the housing-cost gap: suburban households are larger and buy more housing).

## (e) Income gradient of willingness to pay (BFM 2007)

Bayer, Ferreira and McMillan (2007), JPE 115(4): 588-638, Table 8 (p. 625), preferred specification (boundary fixed effects, neighborhood sociodemographics), 'Average Test Score +1 SD' column; household income and college row are interaction terms of the heterogeneous-choice model. Mean MWTP $19.69 (SE 7.41) per month per school-average SD of test scores (monthly user cost of housing, 1990 dollars); household income (+$10,000) 1.38 (SE 0.33); college degree or more vs some college or less 13.03 (SE 3.57); children under 18 vs none 7.41 (SE 3.58). Mean household income $54,103 (SD $50,719; BFM (2004), NBER w10871, Table 1 (same 242,100 households), via literature_targets.md section 1d). The education variable is the household head's.

**Elasticity at the means.** (1.38 / 10,000) x 54,103 / 19.69 = **0.3792**; SE 0.0907 holding the means fixed (0.33 / 10,000 x 54,103 / 19.69), 0.169 if the SE of the mean MWTP (7.41) is added independently (Table 8 gives no covariance; with demeaned interaction characteristics the mean and the slopes are approximately uncorrelated). The elasticity of the WTP *share* of income is -0.621: school quality is a necessity in BFM's data. In the linear Table 8 gradient, WTP as a share of monthly income is 0.71% at 0.5 times the mean income, 0.44% at 1.0 times the mean income, 0.30% at 2.0 times the mean income.

### Folding in the college-degree term

Table 8 holds the other characteristics fixed. If households with more income are also more likely to be college educated, the total gradient of WTP with respect to income is 1.38 per $10,000 plus 13.03 times the slope of the college probability on income. That slope is measured in the Census 1990 PUMS 5% sample for the six Bay Area counties (California file (https://www2.census.gov/census_1990/1990_PUMS_A/PUMSAXCA.zip); households in the Oakland, San Francisco and San Jose PMSAs; householder with a bachelor's degree or more; total household income in 1989).

| statistic | PUMS 1990 Bay Area | BFM sample |
|---|---|---|
| households (records / weighted) | 93,954 / 1,939,709 | 242,100 |
| mean household income $ | 50,775 | 54,103 |
| SD of household income $ | 42,390 | 50,719 |
| householders with a bachelor's degree or more | 35.5% | 44% (block-group average) |
| mean income, college / non-college $ | 68,901 / 40,786 |  |
| correlation of college with income | 0.317 |  |
| slope of college on income, per $10,000 (linear projection) | 0.0358 (bootstrap SE 0.00041) |  |
| local slope at $54,103: logit / local linear, per $10,000 | 0.0402 / 0.0445 |  |

| slope of college on income | total WTP gradient per $10,000 | of which college term | elasticity | SE | college part of the elasticity |
|---|---|---|---|---|---|
| linear projection (Cov / Var), PUMS | 1.847 | 0.467 | 0.508 | 0.097 | 0.128 |
| local logit at the BFM mean, PUMS | 1.904 | 0.524 | 0.523 | 0.099 | 0.144 |
| local linear, +-$20,000 window, PUMS | 1.960 | 0.580 | 0.538 | 0.101 | 0.159 |
| linear projection with BFM's moments (44% college, SD $50,719) and the PUMS correlation | 1.785 | 0.405 | 0.490 | 0.096 | 0.111 |

Assumptions: independent sampling errors of the two Table 8 coefficients (income and college are positively correlated, so their covariance is probably negative and the SE conservative); the children, race, work-status and age terms are held fixed; the PUMS Bay Area sample differs from BFM's (mean income $50,775 vs $54,103, 35.5% vs 44% college heads), so the last row rescales the slope to BFM's reported moments. **Approximate elasticity with the college term: 0.508 (SE 0.097), 0.49-0.54 across the four slope measures.**

### Education definition and the children term

Table 8 does not define 'college degree or more' and BFM's Table 1 has 44% of the heads of household of a census block group with a college degree or more (block-group average); the PUMS share of Bay Area householders with a bachelor's degree or more is 35.5% and with an associate degree or more (YEARSCH >= 12) 43.1%, so both are shown. Table 8 also has a children-under-18 coefficient (7.41, SE 3.58); the 1990 PUMS share of households with a person under 18 is 32.7% and its linear slope on household income 0.0113 per $10,000 (bootstrap SE 0.00043). Slopes are linear projections Cov / Var; SEs combine the independent sampling errors of the Table 8 coefficients and the bootstrap SEs of the slopes (300 resamples of the PUMS households).

| college definition (and children term) | PUMS share | slope of college on income per $10,000 (SE) | total WTP gradient per $10,000 | college (+ children) part of the elasticity | elasticity | SE |
|---|---|---|---|---|---|---|
| bachelor's degree or more | 35.5% | 0.0358 (0.00041) | 1.847 | 0.128 | 0.508 | 0.097 |
| associate degree or more | 43.1% | 0.0357 (0.00041) | 1.845 | 0.128 | 0.507 | 0.097 |
| bachelor's degree or more + children | 35.5% | 0.0358 (0.00041) | 1.930 | 0.128 + 0.023 | 0.530 | 0.098 |
| associate degree or more + children | 43.1% | 0.0357 (0.00041) | 1.928 | 0.128 + 0.023 | 0.530 | 0.098 |

The associate-or-more definition reproduces BFM's 44% (and the slope, 0.0357 against 0.0358 per $10,000) and leaves the elasticity unchanged at 0.507; adding the children term raises it by 0.023 to 0.530. The range across the slope measures and these variants is 0.490-0.538.

## Caveats

- **(a) Cost, not price.** ACS monthly housing costs are gross rent (rent plus utilities) and owners' selected monthly costs (mortgage payments including principal, taxes, insurance, utilities; no maintenance); a percentile gap mixes the price of housing services, quantity and quality, and, for owners, saving. Outright owners' cash costs (taxes, insurance, utilities) are not the price of the housing consumed, and they fill the bottom of the pooled distribution, so only the renter (gross rent) profile is recommended for the lower tail. Units with no cash rent are excluded. Percentiles come from binned counts (linear interpolation inside bins; the median reproduces the published district medians to a mean log error of 0.002, `housing_costs.md`). Elementary and secondary district layers overlap geographically with unified ones in a few places, as in `housing_costs.py`. The tenure split is within the same CZ partition, so the renter / owner gaps include the tenure-specific mix of structures, and none of the percentile gaps is adjusted for household size (suburban households are larger and buy more housing).
- **(b) Units, household size and income.** The CE counts consumer units (2.5 persons, 1.3 earners), not adults, and expenditures, not income; its income quintiles mix permanent and transitory income (expenditures of the lowest quintile are 2.3 times its income). Household size rises with income (1.7 to 3.2 persons per CU), which the equivalence scales address by assumption: the scale is a convention (the OECD table of elasticities runs from 0 to 1), the central square root is supported by the fit of the curve (two degrees of freedom) and not estimated, people and children per CU are published to one decimal and applied to the group means (the within-group dispersion of household size is ignored), 'children' are under 18, and the CE Table 1400 cross-check suggests that housing economies of scale exceed the square root (one-person households spend more than the equivalized curve implies), which would raise the subsistence share of single adults above the central value. The ratio h / (mean expenditure per equivalent adult) is the model's hmin only under the assumption that housing and consumption scale together with household size (`housing_costs.md`). Total expenditures include personal insurance and pensions (11.9%), which the model treats as saving; the net-of-pensions fit is reported at the consumer-unit level only. CE shelter excludes mortgage principal; utilities include telephone service. The Stone-Geary shares ignore the price of housing (relative prices differ across areas), and two parameters (three with the scale) are fitted to five or ten group means, so the residual SEs have 3 (8) degrees of freedom. The SE (CE) treats item and total as independent draws (conservative). The tercile shares are consumer-unit income terciles, so they include household-size differences. ACS burdens are relative to household income with no size adjustment, include utilities and principal, and the mean burden needs class values for the open classes.
- **(c) Elasticities.** The recommended values are constructed (a metro anchor split by a within-metro ratio), not estimated for the model's two locations. The horizon conversion assumes full stock adjustment with constant depreciation (3.5% a year, from the Baum-Snow and Duranton slide; 2% would give factors of 2.20 and 5.47) and equal price changes in city and suburb; it is applied to Saiz's 30-year estimates (1970-2000, metro-wide, housing units, no separation of city and suburb) and to Baum-Snow and Han's 2000-2010 tract predictions (identified from within-region variation in access to jobs) alike. The Online Appendix D numbers are read from the accepted manuscript of September 2023 (the text refers to 237 metros, the Table A8 notes to predicted values for the full Saiz sample of 269 metros; the 1.25 is the text's figure, Table A8 reports 0.80-1.30 across specifications) and may differ in the published JPE version. The city / suburb split of the tracts is by distance to the CBD at the pupil share of the city, not by principal-city boundaries. The model's H is a composite of quantity and quality: floorspace elasticities are higher than unit elasticities.
- **(d) Vehicles.** B25044, B25046 and B08201 count vehicles available to occupied housing units (households, not adults). The household-size classes use midpoints for 0-3 vehicles and a calibrated mean for '4 or more'; the 4-or-more-person class is open-ended; the within-class regressions use district medians of household income, so they control for income imperfectly (within-district dispersion), and the standardized gap uses national class shares; the renter share absorbs most of the remaining difference but tenure is itself a location choice, and the regression with household-size shares as controls over-controls (the size mix is a location outcome). A car is partly a normal good (vehicle count elasticity about 0.65 in the CE, spending elasticity about 1), so the fixed-cost reading applies to the first vehicle; the CE vehicles per consumer unit are published to one decimal (+-2.6% on the cost per vehicle, and the elasticity is computed from rounded counts); the AAA figures are for new vehicles driven 10,000-20,000 miles. The conversion divides household costs by expenditures per consumer unit as `housing_costs.py` does; the one-person variant divides by one-person consumer units' expenditures.
- **(e) BFM.** The WTP gradient is a coefficient of a heterogeneous-choice model with 135 interactions, conditional on the other household characteristics; Table 8 gives no covariances; mean income comes from the working paper version. The 1990 association uses the public PUMS 5% sample for the six counties, not BFM's restricted long-form data, and the householder's education; BFM do not define 'college degree or more' (bachelor's and associate-or-more are both shown) and the children term is the household-level presence of a person under 18.

## Sources

- ACS 5-year 2014-18 (vintage 2018) school-district tables B25104 (via `acs_housing.py`), B25063, B25087, B25044, B25046, B08201 (Census API, `housing_transport.py fetch`); national B25074, B25095, B25106, B19001, B19025, B25003 (acs5 and acs1 2018).
- BLS Consumer Expenditure Survey 2018: Table 1101 (quintiles of income before taxes), Table 1110 (deciles) and Table 1400 (size of consumer unit) at https://www.bls.gov/cex/tables/calendar-year/mean-item-share-average-standard-error.htm; Table R-1 mean expenditures ($61,224.13) via `housing_costs.py`.
- M. F. Foerster and M. Mira d'Ercole (OECD), 'The OECD approach to measuring income distribution and poverty: strengths, limits and statistical issues', Table 2 (equivalence scales and elasticities: per capita 1, 'old' OECD 0.73, OECD-modified 0.53, square root 0.50, household 0), https://cipe.umd.edu/conferences/oecdumd/conf_papers/Papers/The_OECD_Approach_to_Measuring_Income_Distribution_and_Poverty.pdf; B. Buhmann, L. Rainwater, G. Schmaus and T. Smeeding (1988), Review of Income and Wealth 34(2): 115-142 (sensitivity of inequality and poverty to the equivalence elasticity, LIS).
- AAA, Your Driving Costs, 2018 edition: https://exchange.aaa.com/wp-content/uploads/2018/09/18-0090_2018-Your-Driving-Costs-Brochure_FNL-Lo-5-2.pdf
- A. Saiz (2010), QJE 125(3): 1253-1296, https://doi.org/10.1162/qjec.2010.125.3.1253; MSA elasticities file via https://urbaneconomics.mit.edu/research/data.
- N. Baum-Snow and L. Han (2024), JPE 132(6): 1897-1946, https://doi.org/10.1086/728110; accepted version (September 5, 2023; introduction p. 4, Figure 4 discussion p. 33, Online Appendix D pp. 17-18 and Table A8 p. 26) and tract elasticity estimates (gammas_hat_all.dta, 'Explanation of Housing Supply Elasticity Estimates', September 2023) from the authors' page https://sites.google.com/site/baumsnow/research; replication package https://doi.org/10.7910/DVN/LFGAAW (not downloaded).
- N. Baum-Snow and G. Duranton, 'Housing Supply and Housing Affordability', chapter for the Handbook of Regional and Urban Economics vol. 6; slides of the Urban Economics Association summer school, July 2025, https://urbaneconomics.org/workshops/summerschool2025/slides/baumsnow.pdf (slide 73 Saiz, slide 74 Baum-Snow and Han, slide 76 long-run conversion, slide 77 aggregation).
- P. Bayer, F. Ferreira and R. McMillan (2007), JPE 115(4): 588-638, Table 8 (p. 625), https://doi.org/10.1086/522381; Table 1 (44% of block-group heads with a college degree or more); mean household income from NBER w10871 Table 1 (`literature_targets.md`).
- Census 1990 PUMS 5% sample (PUMS A), California: https://www2.census.gov/census_1990/1990_PUMS_A/PUMSAXCA.zip, data dictionary PUMSUSDD.doc in the same directory's DOCUMENT folder (YEARSCH 12-13 associate, 14 bachelor's, R18UNDR).
- Internal: `housing_costs.md` (r_p, median gap), `district_settings.py` (draws), `literature_targets.md` (BFM), `calibration_log.md`, `runs/experiments-2026-10b/DESIGN.md` and `DESIGN_v2.md` (M1, M1b), `M1B_SPEC.md`.

