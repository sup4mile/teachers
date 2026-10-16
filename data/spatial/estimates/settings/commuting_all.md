# Residence against workplace of all workers and of non-teachers (D1b)

Generated 2026-10-05T02:08:08 by `data/spatial/commuting_all.py`, which imports `teacher_residence.py` (unchanged) and builds, for all workers and for non-teachers, the table of `teacher_residence.md` section "Residence against workplace of K-12 teachers on one geography" (JSON key `flows`): ACS 2018 1-year PUMS, employed persons 16+ with a US POWPUMA of work, residence and workplace classified on the same POWPUMA classes (city / suburb / other), replicate-weight standard errors in parentheses (percentage points for shares). Machine-readable: `commuting_all.json` (`groups.<group>.universes.<universe>`).

## Headline

Sample: 1,445,532 persons, 152.6 million weighted; K-12 teachers (OCCP 2310-2330) are 3.03% of them.

| group | live: city | suburb | other | work: city | suburb | other |
|---|---|---|---|---|---|---|
| **all workers** | 31.2% (0.03) | 57.0% (0.03) | 11.8% (0.02) | 35.2% (0.04) | 54.0% (0.04) | 10.8% (0.02) |
| **non-teacher workers** | 31.3% (0.03) | 57.0% (0.03) | 11.8% (0.02) | 35.4% (0.04) | 53.9% (0.04) | 10.7% (0.02) |

| group | P(work in suburb \| live in city) | P(work in city \| live in suburb) | city<->suburb crossing share |
|---|---|---|---|
| **all workers** | 9.09% (0.06) | 11.38% (0.04) | 10.70% (0.03) |
| **non-teacher workers** | 9.11% (0.06) | 11.48% (0.04) | 10.77% (0.03) |
| K-12 teachers (reference) | 8.52% (0.31) | 8.19% (0.19) | 8.41% (0.15) |
| public K-12 teachers (reference) | 8.65% (0.37) | 7.49% (0.20) | 7.97% (0.19) |

| group | block: city share of jobs | city share of residents | net inflow into the city | marginal: city share of jobs | city share of residents |
|---|---|---|---|---|---|
| **all workers** | 39.67% (0.04) | 35.48% (0.03) | +4.19 pp (0.03) | 39.50% (0.04) | 35.36% (0.03) |
| **non-teacher workers** | 39.80% (0.04) | 35.57% (0.03) | +4.23 pp (0.03) | 39.63% (0.04) | 35.45% (0.03) |
| K-12 teachers (reference) | 35.38% (0.31) | 32.59% (0.30) | +2.79 pp (0.17) | 35.19% (0.31) | 32.50% (0.30) |
| public K-12 teachers (reference) | 33.83% (0.35) | 31.35% (0.36) | +2.48 pp (0.17) | 33.65% (0.35) | 31.27% (0.35) |

The crossing share is the share of workers who live and work in a city or suburb POWPUMA and live in one and work in the other (`share_metro_cross` of teacher_residence.py). `block` shares are among those same workers (a closed city-suburb economy); `marginal` shares are among workers who live (work) in a city or suburb POWPUMA whatever the other end. The K-12 rows are recomputed by the same code and reproduce `teacher_residence.json` (section 6).

- **All workers.** 9.1% (0.06) of city residents work in a suburban POWPUMA and 11.4% (0.04) of suburban residents work in a city POWPUMA; 10.7% (0.03) of the workers who live and work in the city or suburb cross between them. The city holds 39.7% (0.04) of those workers' jobs and 35.5% (0.03) of their residents (net inflow +4.19 pp (0.03)).
- **Non-teachers** are almost the same as all workers (teachers are 3.0% of them): 9.1% (0.06), 11.5% (0.04), crossing share 10.8% (0.03).
- **K-12 teachers against non-teachers** (difference, SE): P(work in suburb | live in city) -0.58 pp (0.32), about as often; P(work in city | live in suburb) -3.30 pp (0.19), less often; crossing share -2.36 pp (0.16); net commuting gap -0.058 (0.008) log points (section 4).

## Methods note

- Data and sample: the ACS 2018 1-year PUMS cached by teacher_residence.py (employed persons 16+, ESR 1-2, in households; 80 replicate weights), restricted as in its `flows_block` to persons with a US place-of-work POWPUMA (POWSP 1-56, POWPUMA > 1). Workers at home carry their own POWPUMA as the place of work (checked: all of them do), so they count as working where they live; persons employed but not at work last week and workplaces abroad or in Puerto Rico are out.
- Classes: `teacher_residence.build_puma_geo()` at the main cut (0.5 / 0.5): a POWPUMA is `city` if at least half of its 2010 population lives in principal cities of OMB (Sept. 2018) metropolitan areas, `suburb` if at least half lives in metropolitan counties and it is not `city`, `other` otherwise. Residence class = class of the POWPUMA that contains the residence PUMA; workplace class = class of the place-of-work POWPUMA (one geography at both ends).
- Groups: all workers; non-teachers = everyone except OCCP 2310-2330 (the `k12` flag of teacher_residence.py). K-12 and public K-12 (COW 3-4) teachers are recomputed by the same code as a check and reproduce `flows.k12` / `flows.pub_k12` of teacher_residence.json.
- Estimator: person-weighted cell shares; every probability is a ratio of cells. Standard errors from the 80 ACS replicate weights, SE = sqrt(4/80 sum_r (x_r - x)^2), with derived statistics and between-group differences evaluated on every replicate (covariances respected). Sampling error only.
- Definitions: `block` = workers who live AND work in a city or suburb POWPUMA (the closed two-location system); crossing share = (city->suburb + suburb->city) / block (`share_metro_cross` of teacher_residence.py); `marginal` = workers who live (work) in a city or suburb POWPUMA, whatever the other end. Caveat: a POWPUMA is a whole county or group of counties, so city-suburb commuting inside one county is invisible and the crossing probabilities are lower bounds; workers at home count as non-crossers (see `away`); the classes depend on the cuts (`class_cuts`).

## 1. Where workers live and where they work

Shares of each group's workers (the three residence shares and the three workplace shares each sum to 100%). The last four columns are the log gaps of teacher_residence.py, log(suburb / city) by residence, by workplace and their difference (the net commuting gap), and the share working in a different class from the one they live in.

| group | n | workers (millions) | live: city | suburb | other | work: city | suburb | other | gap by residence | gap by workplace | net commuting gap | works in the other class |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all workers | 1,445,532 | 152.56 | 31.2% (0.03) | 57.0% (0.03) | 11.8% (0.02) | 35.2% (0.04) | 54.0% (0.04) | 10.8% (0.02) | 0.603 (0.001) | 0.426 (0.002) | 0.177 (0.001) | 12.5% (0.04) |
| non-teacher workers | 1,396,928 | 147.94 | 31.3% (0.03) | 57.0% (0.03) | 11.8% (0.02) | 35.4% (0.04) | 53.9% (0.04) | 10.7% (0.02) | 0.599 (0.001) | 0.421 (0.002) | 0.179 (0.001) | 12.5% (0.04) |
| K-12 teachers (reference) | 48,604 | 4.62 | 28.4% (0.26) | 59.0% (0.29) | 12.7% (0.17) | 30.8% (0.28) | 56.7% (0.30) | 12.4% (0.18) | 0.731 (0.013) | 0.611 (0.014) | 0.120 (0.008) | 9.9% (0.14) |
| public K-12 teachers (reference) | 34,338 | 3.26 | 26.9% (0.31) | 59.2% (0.33) | 13.9% (0.21) | 29.1% (0.31) | 57.3% (0.34) | 13.7% (0.23) | 0.787 (0.016) | 0.679 (0.016) | 0.108 (0.008) | 9.4% (0.16) |

City and suburb shares among workers who live (work) in one of the two, whatever the other end (`marginal`), and among workers who live and work in one of the two (`block`):

| group | marginal: city share of residents | city share of jobs | block: share of the group | city share of residents | city share of jobs |
|---|---|---|---|---|---|
| all workers | 35.36% (0.03) | 39.50% (0.04) | 87.1% (0.02) | 35.48% (0.03) | 39.67% (0.04) |
| non-teacher workers | 35.45% (0.03) | 39.63% (0.04) | 87.2% (0.03) | 35.57% (0.03) | 39.80% (0.04) |
| K-12 teachers (reference) | 32.50% (0.30) | 35.19% (0.31) | 86.1% (0.18) | 32.59% (0.30) | 35.38% (0.31) |
| public K-12 teachers (reference) | 31.27% (0.35) | 33.65% (0.35) | 84.9% (0.22) | 31.35% (0.36) | 33.83% (0.35) |

## 2. The 3x3 residence x workplace matrices

Joint shares of the group's workers; the last column is the residence share and the last row the workplace share. Below each, the row-conditional probabilities P(work in b | live in a) and the column-conditional P(live in a | work in b).

**all workers** (1,445,532 persons, 152.6 million; the smallest cell holds 3,712 persons)

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 28.08% (0.03) | 2.84% (0.02) | 0.27% (0.01) | 31.19% (0.03) |
| suburb | 6.49% (0.02) | 49.74% (0.04) | 0.79% (0.01) | 57.01% (0.03) |
| other | 0.67% (0.01) | 1.41% (0.01) | 9.71% (0.02) | 11.80% (0.02) |
| all | 35.25% (0.04) | 53.98% (0.04) | 10.77% (0.02) | 100% |

| P(work in ... \| live in ...) | city | suburb | other |
|---|---|---|---|
| live in city | 90.04% (0.06) | 9.09% (0.06) | 0.87% (0.02) |
| live in suburb | 11.38% (0.04) | 87.24% (0.04) | 1.38% (0.01) |
| live in other | 5.72% (0.06) | 11.95% (0.09) | 82.33% (0.10) |

| P(live in ... \| work in ...) | live: city | suburb | other |
|---|---|---|---|
| work in city | 79.68% (0.06) | 18.41% (0.05) | 1.91% (0.02) |
| work in suburb | 5.25% (0.03) | 92.14% (0.04) | 2.61% (0.02) |
| work in other | 2.52% (0.05) | 7.32% (0.07) | 90.15% (0.09) |

**non-teacher workers** (1,396,928 persons, 147.9 million; the smallest cell holds 3,566 persons)

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 28.16% (0.04) | 2.85% (0.02) | 0.27% (0.01) | 31.28% (0.03) |
| suburb | 6.54% (0.02) | 49.63% (0.04) | 0.79% (0.01) | 56.95% (0.03) |
| other | 0.69% (0.01) | 1.42% (0.01) | 9.66% (0.02) | 11.77% (0.02) |
| all | 35.38% (0.04) | 53.89% (0.04) | 10.72% (0.02) | 100% |

| P(work in ... \| live in ...) | city | suburb | other |
|---|---|---|---|
| live in city | 90.03% (0.06) | 9.11% (0.06) | 0.87% (0.02) |
| live in suburb | 11.48% (0.04) | 87.14% (0.04) | 1.38% (0.02) |
| live in other | 5.82% (0.06) | 12.06% (0.09) | 82.11% (0.10) |

| P(live in ... \| work in ...) | live: city | suburb | other |
|---|---|---|---|
| work in city | 79.58% (0.06) | 18.48% (0.05) | 1.94% (0.02) |
| work in suburb | 5.29% (0.03) | 92.08% (0.04) | 2.63% (0.02) |
| work in other | 2.52% (0.05) | 7.33% (0.08) | 90.15% (0.09) |

**K-12 teachers (reference)** (48,604 persons, 4.6 million; the smallest cell holds 146 persons)

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 25.65% (0.27) | 2.42% (0.09) | 0.31% (0.03) | 28.39% (0.26) |
| suburb | 4.83% (0.11) | 53.25% (0.28) | 0.89% (0.05) | 58.96% (0.29) |
| other | 0.33% (0.03) | 1.08% (0.04) | 11.24% (0.17) | 12.65% (0.17) |
| all | 30.81% (0.28) | 56.75% (0.30) | 12.44% (0.18) | 100% |

| P(work in ... \| live in ...) | city | suburb | other |
|---|---|---|---|
| live in city | 90.37% (0.32) | 8.52% (0.31) | 1.11% (0.11) |
| live in suburb | 8.19% (0.19) | 90.31% (0.18) | 1.51% (0.09) |
| live in other | 2.62% (0.21) | 8.55% (0.38) | 88.83% (0.39) |

| P(live in ... \| work in ...) | live: city | suburb | other |
|---|---|---|---|
| work in city | 83.26% (0.35) | 15.66% (0.35) | 1.08% (0.09) |
| work in suburb | 4.26% (0.15) | 93.83% (0.16) | 1.91% (0.08) |
| work in other | 2.53% (0.23) | 7.13% (0.42) | 90.34% (0.46) |

**public K-12 teachers (reference)** (34,338 persons, 3.3 million; the smallest cell holds 105 persons)

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 24.28% (0.30) | 2.33% (0.10) | 0.33% (0.03) | 26.94% (0.31) |
| suburb | 4.44% (0.12) | 53.85% (0.33) | 0.91% (0.07) | 59.20% (0.33) |
| other | 0.33% (0.04) | 1.11% (0.05) | 12.42% (0.22) | 13.87% (0.21) |
| all | 29.05% (0.31) | 57.29% (0.34) | 13.66% (0.23) | 100% |

| P(work in ... \| live in ...) | city | suburb | other |
|---|---|---|---|
| live in city | 90.14% (0.38) | 8.65% (0.37) | 1.21% (0.12) |
| live in suburb | 7.49% (0.20) | 90.97% (0.21) | 1.54% (0.12) |
| live in other | 2.40% (0.25) | 8.02% (0.41) | 89.59% (0.42) |

| P(live in ... \| work in ...) | live: city | suburb | other |
|---|---|---|---|
| work in city | 83.59% (0.42) | 15.27% (0.41) | 1.14% (0.12) |
| work in suburb | 4.07% (0.17) | 93.99% (0.17) | 1.94% (0.09) |
| work in other | 2.38% (0.24) | 6.67% (0.49) | 90.95% (0.51) |

## 3. Residents of metropolitan areas only

The script already distinguishes metropolitan residents (`metro` = residence class city or suburb, used in its commute block on PUMA classes). Two versions: (a) on the geography of this table, residence POWPUMA class city or suburb (rows `city` and `suburb` of the matrix; the conditional probabilities, the crossing share and the block shares are then identical to section 2 by construction, only the margins change); (b) residence PUMA class city or suburb (the PUMA-level `metro` flag), cells on the POWPUMA classes.

| group | version | n | workers (millions) | live: city | suburb | other | work: city | suburb | other (outside the two classes) |
|---|---|---|---|---|---|---|---|---|---|
| all workers | (a) POWPUMA | 1,221,375 | 134.56 | 35.4% (0.03) | 64.6% (0.03) | - | 39.2% (0.04) | 59.6% (0.04) | 1.2% (0.01) |
| all workers | (b) PUMA flag | 1,207,838 | 133.63 | 35.6% (0.03) | 64.0% (0.03) | 0.3% (0.00) | 39.5% (0.04) | 59.1% (0.04) | 1.4% (0.01) |
| non-teacher workers | (a) POWPUMA | 1,180,714 | 130.53 | 35.5% (0.03) | 64.5% (0.03) | - | 39.3% (0.04) | 59.5% (0.04) | 1.2% (0.01) |
| non-teacher workers | (b) PUMA flag | 1,167,750 | 129.63 | 35.7% (0.03) | 64.0% (0.03) | 0.3% (0.00) | 39.6% (0.04) | 59.0% (0.04) | 1.4% (0.01) |
| K-12 teachers (reference) | (a) POWPUMA | 40,661 | 4.04 | 32.5% (0.30) | 67.5% (0.30) | - | 34.9% (0.31) | 63.7% (0.31) | 1.4% (0.07) |
| K-12 teachers (reference) | (b) PUMA flag | 40,088 | 4.00 | 32.8% (0.29) | 66.8% (0.30) | 0.4% (0.04) | 35.2% (0.31) | 63.2% (0.31) | 1.7% (0.08) |
| public K-12 teachers (reference) | (a) POWPUMA | 28,331 | 2.81 | 31.3% (0.35) | 68.7% (0.35) | - | 33.3% (0.35) | 65.2% (0.35) | 1.4% (0.09) |
| public K-12 teachers (reference) | (b) PUMA flag | 27,914 | 2.78 | 31.6% (0.34) | 68.0% (0.35) | 0.4% (0.05) | 33.6% (0.35) | 64.6% (0.34) | 1.8% (0.10) |

| group | version | P(work in suburb \| live in city) | P(work in city \| live in suburb) | city<->suburb crossing share | city share of jobs (block) | city share of residents (block) | works in the other class |
|---|---|---|---|---|---|---|---|
| all workers | (a) POWPUMA | 9.09% (0.06) | 11.38% (0.04) | 10.70% (0.03) | 39.67% (0.04) | 35.48% (0.03) | 11.8% (0.04) |
| all workers | (b) PUMA flag | 9.09% (0.06) | 11.50% (0.04) | 10.76% (0.03) | 40.03% (0.04) | 35.84% (0.03) | 11.9% (0.04) |
| non-teacher workers | (a) POWPUMA | 9.11% (0.06) | 11.48% (0.04) | 10.77% (0.03) | 39.80% (0.04) | 35.57% (0.03) | 11.8% (0.04) |
| non-teacher workers | (b) PUMA flag | 9.11% (0.06) | 11.60% (0.04) | 10.83% (0.03) | 40.16% (0.04) | 35.93% (0.03) | 11.9% (0.04) |
| K-12 teachers (reference) | (a) POWPUMA | 8.52% (0.31) | 8.19% (0.19) | 8.41% (0.15) | 35.38% (0.31) | 32.59% (0.30) | 9.7% (0.15) |
| K-12 teachers (reference) | (b) PUMA flag | 8.52% (0.31) | 8.31% (0.19) | 8.50% (0.16) | 35.80% (0.31) | 32.99% (0.30) | 9.8% (0.15) |
| public K-12 teachers (reference) | (a) POWPUMA | 8.65% (0.37) | 7.49% (0.20) | 7.97% (0.19) | 33.83% (0.35) | 31.35% (0.36) | 9.3% (0.19) |
| public K-12 teachers (reference) | (b) PUMA flag | 8.65% (0.37) | 7.61% (0.20) | 8.05% (0.19) | 34.26% (0.35) | 31.77% (0.35) | 9.4% (0.19) |

Metropolitan-resident matrices, version (a): joint shares of the group's metropolitan residents.

**all workers**

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 31.84% (0.04) | 3.21% (0.02) | 0.31% (0.01) | 35.36% (0.03) |
| suburb | 7.36% (0.02) | 56.39% (0.04) | 0.89% (0.01) | 64.64% (0.03) |
| all | 39.19% (0.04) | 59.60% (0.04) | 1.20% (0.01) | 100% |

**non-teacher workers**

| live \ work | city | suburb | other | all |
|---|---|---|---|---|
| city | 31.92% (0.04) | 3.23% (0.02) | 0.31% (0.01) | 35.45% (0.03) |
| suburb | 7.41% (0.02) | 56.25% (0.04) | 0.89% (0.01) | 64.55% (0.03) |
| all | 39.33% (0.04) | 59.47% (0.04) | 1.20% (0.01) | 100% |

## 4. K-12 teachers minus non-teachers

Differences of the statistics between the groups (same replicate weights, so the standard error respects the covariance). Percentage points for the shares, log points for the gaps.

| statistic | K-12 minus non-teachers | public K-12 minus non-teachers | K-12 minus non-teachers, workers away from home |
|---|---|---|---|
| P(work in suburb \| live in city) | -0.58 pp (0.32) | -0.46 pp (0.38) | -0.95 pp (0.33) |
| P(work in city \| live in suburb) | -3.30 pp (0.19) | -3.99 pp (0.20) | -3.84 pp (0.19) |
| city<->suburb crossing share | -2.36 pp (0.16) | -2.80 pp (0.20) | -2.85 pp (0.16) |
| works in the other class | -2.69 pp (0.14) | -3.10 pp (0.17) | -3.24 pp (0.14) |
| city share of jobs (block) | -4.42 pp (0.32) | -5.98 pp (0.36) | -4.70 pp (0.32) |
| city share of residents (block) | -2.98 pp (0.31) | -4.22 pp (0.36) | -3.06 pp (0.31) |
| net inflow into the city (block) | -1.44 pp (0.17) | -1.75 pp (0.18) | -1.64 pp (0.18) |
| city share of residents (marginal) | -2.95 pp (0.30) | -4.18 pp (0.36) | -3.02 pp (0.31) |
| city share of jobs (marginal) | -4.45 pp (0.32) | -5.99 pp (0.36) | -4.71 pp (0.32) |
| gap by residence (log) | 0.132 (0.014) | 0.188 (0.017) | 0.135 (0.014) |
| gap by workplace (log) | 0.190 (0.014) | 0.258 (0.016) | 0.201 (0.014) |
| net commuting gap (log) | -0.058 (0.008) | -0.070 (0.008) | -0.066 (0.008) |

## 5. Robustness

**Workers away from home** (JWTR != 11: work-at-home workers carry their own area as the workplace and never cross; dropping them raises the crossing probabilities).

| group | n | P(work in suburb \| live in city) | P(work in city \| live in suburb) | city<->suburb crossing share | works in the other class |
|---|---|---|---|---|---|
| all workers | 1,365,545 | 9.59% (0.06) | 12.04% (0.04) | 11.31% (0.03) | 13.2% (0.04) |
| non-teacher workers | 1,317,763 | 9.62% (0.06) | 12.16% (0.04) | 11.40% (0.03) | 13.3% (0.04) |
| K-12 teachers (reference) | 47,782 | 8.66% (0.32) | 8.32% (0.19) | 8.55% (0.16) | 10.0% (0.14) |
| public K-12 teachers (reference) | 33,928 | 8.74% (0.37) | 7.58% (0.20) | 8.06% (0.19) | 9.6% (0.16) |

**Class cuts** (`teacher_residence.CUTS`: a POWPUMA is `city` if its principal-city population share is at least pc, `suburb` if its metropolitan share is at least metro and it is not `city`; the main line is 0.50 / 0.50).

| pc cut | metro cut | POWPUMAs city / suburb / other | group | P(work in suburb \| live in city) | P(work in city \| live in suburb) | crossing share | city share of jobs (block) | city share of residents (block) |
|---|---|---|---|---|---|---|---|---|
| 0.25 | 0.50 | 347 / 325 / 310 | all workers | 5.51% (0.03) | 22.22% (0.07) | 11.60% (0.04) | 68.83% (0.04) | 64.41% (0.03) |
| 0.25 | 0.50 | 347 / 325 / 310 | non-teacher workers | 5.50% (0.04) | 22.42% (0.07) | 11.65% (0.04) | 68.99% (0.04) | 64.52% (0.03) |
| 0.25 | 0.50 | 347 / 325 / 310 | K-12 teachers (reference) | 5.68% (0.16) | 16.36% (0.37) | 10.02% (0.17) | 63.69% (0.27) | 60.67% (0.31) |
| 0.50 | 0.50 (main) | 146 / 523 / 313 | all workers | 9.09% (0.06) | 11.38% (0.04) | 10.70% (0.03) | 39.67% (0.04) | 35.48% (0.03) |
| 0.50 | 0.50 (main) | 146 / 523 / 313 | non-teacher workers | 9.11% (0.06) | 11.48% (0.04) | 10.77% (0.03) | 39.80% (0.04) | 35.57% (0.03) |
| 0.50 | 0.50 (main) | 146 / 523 / 313 | K-12 teachers (reference) | 8.52% (0.31) | 8.19% (0.19) | 8.41% (0.15) | 35.38% (0.31) | 32.59% (0.30) |
| 0.75 | 0.50 | 50 / 619 / 313 | all workers | 13.27% (0.11) | 4.48% (0.02) | 5.51% (0.02) | 13.63% (0.03) | 11.08% (0.02) |
| 0.75 | 0.50 | 50 / 619 / 313 | non-teacher workers | 13.29% (0.12) | 4.51% (0.02) | 5.55% (0.02) | 13.68% (0.03) | 11.10% (0.02) |
| 0.75 | 0.50 | 50 / 619 / 313 | K-12 teachers (reference) | 12.55% (0.57) | 3.32% (0.11) | 4.35% (0.11) | 12.24% (0.23) | 10.56% (0.20) |
| 0.90 | 0.50 | 29 / 640 / 313 | all workers | 14.33% (0.15) | 2.96% (0.02) | 3.78% (0.02) | 8.65% (0.02) | 6.84% (0.02) |
| 0.90 | 0.50 | 29 / 640 / 313 | non-teacher workers | 14.36% (0.16) | 2.97% (0.02) | 3.80% (0.02) | 8.67% (0.03) | 6.86% (0.02) |
| 0.90 | 0.50 | 29 / 640 / 313 | K-12 teachers (reference) | 13.17% (0.78) | 2.50% (0.09) | 3.22% (0.10) | 7.88% (0.19) | 6.34% (0.17) |
| 0.50 | 0.90 | 146 / 402 / 434 | all workers | 8.49% (0.05) | 11.84% (0.04) | 10.74% (0.03) | 42.66% (0.04) | 38.57% (0.03) |
| 0.50 | 0.90 | 146 / 402 / 434 | non-teacher workers | 8.52% (0.06) | 11.94% (0.04) | 10.81% (0.04) | 42.79% (0.04) | 38.66% (0.03) |
| 0.50 | 0.90 | 146 / 402 / 434 | K-12 teachers (reference) | 7.63% (0.30) | 8.68% (0.20) | 8.49% (0.16) | 38.55% (0.32) | 35.61% (0.30) |

## 6. Checks

- **Reproduction of `teacher_residence.json`.** The K-12 and public K-12 rows computed here (gap by residence, gap by workplace, net commuting gap, share working in the other class, city<->suburb crossing share, and the six `live_given_work` probabilities, values and standard errors) equal `flows.k12` and `flows.pub_k12` of the existing file: largest absolute difference 7.8e-16; sample 1,445,532 persons here and 1,445,532 there.
- **Sample.** 1,474,819 employed persons 16+ in households (155.7 million weighted); 1,445,532 (152.6 million) have a US place-of-work POWPUMA; 0 of them could not be classified at one end and are dropped. The 29,287 persons without a US POWPUMA of work (3.2 million weighted) are employed-but-not-at-work persons and workers with a workplace abroad or in Puerto Rico.
- **Work at home.** 5.3% of the sample works at home (79,987 persons); 0 of them report a place of work in a POWPUMA other than the residence's.
- **Teacher definition.** OCCP codes found in 2310-2330: [2310, 2320, 2330] (equal to `k12`, asserted); the non-teacher group keeps OCCP 2300 (preschool and kindergarten teachers) and the other education occupations.
- **Pipeline check** (PUMS loader of teacher_residence.py): mean travel time 27.1 minutes against 27.1 published, work-at-home share 5.2% against 5.3%.
- Standard errors are sampling errors of the PUMS only; the classification of the areas is taken as given.
