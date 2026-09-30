# Spatial data pipeline — `fetch_data.py`

Downloads and packages the public data used to motivate and calibrate the
two-location teachers model. **Acquisition only**: every table is fetched,
parsed, and written to `raw/` as it came. Merging districts to commuting zones,
splitting each CZ into two locations, and computing the within-CZ gaps happen
downstream — see `notes/spatial_calibration_notes.md` for the moment map these
data feed.

Model calibration code, notes, the Julia environment, SLURM helpers and runs live
in [`julia/spatial_model/calibration/`](../../julia/spatial_model/calibration/).
Start with the [calibration note](../../julia/spatial_model/calibration/spatial_calibration.md)
and [calibration log](../../julia/spatial_model/calibration/calibration_log.md).
Data acquisition scripts, raw inputs and generated data reports remain here.

## Quick start

```bash
# uv resolves the environment from pyproject.toml/uv.lock on first run

# see every URL that would be requested, download nothing
uv run fetch_data.py --dry-run

# full run: everything lands in raw/ next to this script
uv run fetch_data.py

# just refresh the finance data
uv run fetch_data.py --sources urban --urban-topics ccd_finance

# the NLSY test batteries for the ability block (opt-in, ~1.1 GB)
uv run fetch_data.py --sources nlsy
```

Expect the full download to take 30–90 minutes and roughly 5–15 GB on disk,
most of it the CCD directory and F-33 panels. `--years 1990 2000 2010 2018`
cuts that to a few minutes if you only want the benchmark cross-sections.

## What it pulls

| Source | What | Coverage | Key required |
|---|---|---|---|
| Urban Institute Education Data API | CCD directory (enrollment, FTE teachers, locale, county), CCD enrollment, **F-33 finance** (revenue by source, expenditure, salaries), SAIPE child poverty, EDFacts assessments and grad rates | 1986–2022 (F-33: 1991, 1994–2018) | no |
| SEDA (Stanford Digital Repository) | district / county / metro / **commuting-zone** mean scores and **learning rates**, plus covariates and the school–district crosswalk | 2009–2019 | no |
| NCES EDGE | **CWIFT** (comparable wage index for teachers, at district / county / state level), district geocode+locale files | 2013–2023 | no |
| David Dorn | 1990 county → 1990 commuting-zone crosswalk (3,141 counties) | — | no |
| *optional* Census ACS | school-district income, home value, property tax | 2009–2022 | yes (`--census-key` or `CENSUS_API_KEY`) |
| *optional* NLS (BLS) | **NLSY79** ASVAB/AFQT (section scores, SEs, item responses), **NLSY79 Child/YA** PIAT/PPVT with mother links, **NLSY97** CAT-ASVAB with posterior variances; plus ages, weights, parental and own schooling, occupation, wages, annual wage income, weeks and hours (young adults: income and weights) | 1979–2023 | no |

ACS is off by default: `--sources urban seda edge czones acs`. Nothing else
needs a key. A free Census key takes about a minute at
<https://api.census.gov/data/key_signup.html>; put it in `.env` as
`CENSUS_API_KEY` and `source` that before the run.

The district layers refuse `in=state:*` before the 2022 vintage, so ACS is
fetched one state at a time (`ACS_WORKERS` in parallel; about 20 minutes for
2009-2022). Variables a vintage does not publish are dropped after a probe
request rather than failing the year: `B25103` (median property tax) starts in
2010 and `B15003` (detailed attainment) in 2012.

NLSY is also off by default: `--sources nlsy`, with `--nlsy-cohorts` to pick
among `nlsy79 nlscya nlsy97`. It feeds the occupations-and-ability block (Table 2 and item
T2 in `spatial_calibration.md`), not the spatial moments. See the NLSY
section below.

`--urban-topics crdc_teachers` additionally pulls school-level CRDC teacher
experience/certification (2011, 2013, 2015, 2017); it is not in the default set.

## What it produces

```
data/spatial/
  raw/             one table per source-endpoint (parquet by default; --format csv|both)
    urban_ccd_directory, urban_ccd_finance, urban_saipe, urban_edfacts_*
    seda_*                (one per SEDA file matched at the requested levels)
    edge_cwift, edge_geocode_lea
    czone_crosswalk
    acs_school_districts  (only with --sources acs)
    nlsy79_ability, nlscya_ability, nlsy97_ability   (only with --sources nlsy)
    _downloads/           the original csv / zip / dta files, untouched
  metadata/
    urban_api_*.csv       endpoint, variable and bulk-download lists
    seda_file_list.txt    what the SDR reported for the release
    SEDA_documentation_*.pdf
    <cohort>_ability_vars.csv   NLSY variable dictionary (reference number,
                                year, question name, title, block)
    <cohort>_ability.NLSY79 / .CHILDYA / .NLSY97   the same selection as an
                                NLS Investigator tagset
    manifest.json         every URL attempted, its status and row count
```

Parquet is the default because `_downloads/` already holds a verbatim text copy;
pass `--format both` if you want CSV alongside it.

## Notes on the sources

- **Urban bulk CSVs.** For a full annual panel the script prefers the portal's
  one-file-per-endpoint bulk CSV over paging the API. `--no-bulk` forces paging;
  `--years` implies paging, since a bulk file covers all years at once.
- **SEDA file names.** The file list comes from the Stanford Digital Repository
  at run time, so a new release only needs its druid (`--seda-druid`, or a new
  entry in `SEDA_DRUIDS`). Name patterns are only a fallback.
- **Archive quirks.** Several source zips were built on a Mac and carry
  `__MACOSX` resource forks that parse as one-row tables; NCES ships
  tab- and pipe-delimited files under a `.txt` extension; and the EDGE geocode
  archive gives the same table as a named-column `.xlsx` and a header-less
  `.TXT`. The zip reader skips the forks, sniffs the delimiter, and prefers the
  spreadsheet — which is why `openpyxl` is a hard requirement, not an optional.
- **NLSY without Investigator.** NLS Investigator needs a login, but the NLS
  also publishes each cohort whole at `nlsinfo.org/cohort-data/` and links the
  files from its [access page](https://www.nlsinfo.org/accessing-data-cohorts).
  The script reads the current file names from that page, so a new release
  needs no code change. Each zip holds a CSV of every variable, 2.9–8.1 GB
  uncompressed, with reference numbers (`R0618301`) as column names, plus a
  codebook and a variable index (`.sdf`). The script reads the index, selects
  variables by **question name** with the regexes in `NLSY_BLOCKS`, and cuts
  only those columns from the CSV. The selection is 493 NLSY79, 354 Child/YA
  and 340 NLSY97 variables. Selecting by question name (`MATH1996`, `ASVAB-3`,
  `KEY!SEX`) means new survey rounds come in without editing a list of
  R-numbers. To add a variable, add its question name to a block. The zips are
  kept in `_downloads/` and reused when a rerun only changes the selection.
- **NLSY values are left as released.** Negative values are NLS missing-data
  codes (-1 refused, -2 don't know, -3 invalid skip, -4 valid skip,
  -5 non-interview), not scores. Many fields carry implied decimals. NLSY79
  `AFQT-3` and NLSY97 `ASVAB_MATH_VERBAL_SCORE_PCT` run 0–100000, so they are
  percentiles to three decimals. Sampling weights carry two. The NLSY97
  CAT-ASVAB ability estimate is split into a `_POS` and a `_NEG` field. The
  codebook (`.cdb`) in each zip gives every variable's scaling. To link a child to its mother, match Child/YA `MPUBID`
  (`C0000200`) to NLSY79 `CASEID` (`R0000100`).
- **Reruns are cheap to scope.** Every fetch is independent, so re-pulling one
  source (`--sources edge`) overwrites only that source's tables.

---

# Moment construction — `data_estimate.py`

Reads whatever `fetch_data.py` has left in `raw/` and turns it into the spatial
moments that discipline `julia/spatial_model/spatial_continuous.jl`. Every
headline number is a **within-commuting-zone gap**, enrollment-weighted across
CZs, following the Option 1 interpretation argued in
`notes/spatial_calibration_notes.md`: a *location* is a group of school
districts inside one labor market, and location 2 is always the advantaged
group, matching the model's `κ = [0.75, 0.9]`, `B = [0.0, 0.1]`.

```bash
uv run data_estimate.py                       # base year 2018, locale split
uv run data_estimate.py --scheme all          # all three CZ partitions
uv run data_estimate.py --base-year 2010 --years 1994 2000 2010 2018
uv run data_estimate.py --geo cbsa            # force the metro fallback
```

One partition takes about 10 seconds; `--scheme all` takes about 70, almost all
of it in the CZ-by-year loop over the two median splits.

## Partitioning a CZ in two

| `--scheme` | location 1 / location 2 | needs |
|---|---|---|
| `locale` (default) | city (NCES locale 11–13) / suburb (21–23) | CCD, 2006 on |
| `revenue` | below / above the CZ enrollment-weighted median of local revenue per pupil | F-33 |
| `ses` | below / above the CZ median SEDA SES (`sesavgall`, from the district covariate file; falls back to SAIPE child poverty) | SEDA or SAIPE |

A CZ is kept only when both groups clear 10% of its enrollment
(`--min-group-share`), and districts below 500 pupils are dropped
(`--min-pupils`). `--geo` picks the labor market: Dorn's commuting zones by
default, falling back to the CCD's CBSA code with a warning when
`--sources czones` has not been run.

## What it estimates

- **Gaps** in the objects the model has: `M_2/M_1`, `H̃_T,l/M_l`, `κ_2/κ_1`
  (CWIFT-deflated), `Q_2/Q_1` (SEDA learning rates), per-pupil revenue by
  source, the effective local tax rate, and the poverty / SES / median-income
  gaps that stand in for sorting on ability. Location totals are summed and per-pupil objects are
  ratios of totals — *not* averages of district ratios, which would answer a
  question the model does not ask.
- **Regressions** with CZ fixed effects and clustered SEs: σ (congestion),
  σ_ν (mobility), the κ–quality gradient, and the local-revenue → teacher-pay
  pass-through.
- **A within-CZ variance share** for every gap, on the base-year cross-section —
  the descriptive case for the whole design. Pooling the years instead would
  count two decades of nominal drift in every dollar variable as within-CZ
  variance, since each CZ appears in every year.
- **Literature anchors** for what no public district table can identify (λ, ψ,
  the moving cost, ρ_z), carried with their source rather than silently omitted.

Derived ratios are blanked outside their within-year [0.5%, 99.5%] range:
F-33 wage bills matched to a CCD teacher count from a different reporting
universe produce $2M per FTE and class sizes of 1000.

Three things the raw tables force, which are worth knowing before reading a
number:

- **SEDA files are long in subgroup.** Race, gender and economic-disadvantage
  rows sit under the same district, and rows with `gap == 1` are *differences*
  between two subgroups rather than levels. Only the all-students level rows are
  read — averaging over the lot pulls the district mean toward zero (it moved
  the city/suburb score gap from 0.09 to 0.29 SD).
- **CWIFT starts with the 2013-14 district universe** and stacks district,
  county and state rows in one file. Only the district rows are used, merged on
  the nearest release within two years, so `κ_2/κ_1` exists for 2014 and 2018
  and is genuinely absent before that rather than back-cast.
- **SEDA learning rates are SDs gained per grade**, fitted over grades 3-8. The
  reported gap keeps that published unit; the σ and σ_ν regressions work with
  the gain cumulated over the five grades, since `Q_l` is a level.

## Output

```
data/spatial/estimates/
  spatial_moments.json     every estimate: value, source, and a note on how it
                           disciplines the model; plus the input inventory and
                           the data-gap list
  spatial_moments.md       the same, as a readable report
  czone_gaps_<scheme>.csv  the CZ x year gap panel behind the headline numbers
  gap_trends_<scheme>.csv  each gap by year
```

Nothing is silently skipped. A moment that cannot be built becomes a
`missing-input` record naming the `fetch_data.py` flag that would supply it, so
the report doubles as a to-do list for the empirical section — as does its
closing "what the empirical section still needs".

**Status: run against the full download (2026-08-31).** Every input in the
report's inventory now resolves against a real table — CCD, F-33, SAIPE,
SEDA 6.0, CWIFT, Dorn's crosswalk, EDFacts and ACS — and the commuting zone,
not the CBSA fallback, is the labor market. 83 estimates, none blocked on
missing data.

ACS closed the last gap: `t_eff` (local revenue over the ACS income base) is
3.0% in aggregate for 2018, and the within-CZ gap in it is +0.06pp on the
locale split, +0.55pp on the revenue split. Note the denominator is ACS money
income, which runs well below BEA personal income, so the level is roughly
1.6x what a national-accounts base would give; the within-CZ gap is the object
to read, not the level.

Two results worth reading before trusting them: `σ` and `σ_ν` both come back
wrong-signed and insignificant (t = -0.29 and -0.62). That is the within-CZ
cross-section saying it cannot identify either without an instrument, which is
already the first item on the data-gap list — not a bug in the pipeline.

---

# Ability moments — `nlsy_ability.py`

Reads the three NLSY tables that `fetch_data.py --sources nlsy` leaves in
`raw/` and builds the data side of the ability block (Table 2, item T2a in
`spatial_calibration.md`). T2b uses the latent correlation and wage slope in
closed form, then fits wage dispersion outside equilibrium; it does not need
simulated panels or score noise. Full-equilibrium consistency checks follow
the first internal fit (T6).

```bash
uv run nlsy_ability.py              # B = 200 household bootstrap, about a minute
uv run nlsy_ability.py --quick      # B = 20
uv run nlsy_ability.py --selftest   # factor fit, norming and disattenuation on synthetic data
```

- **Mother–child measurement system.** NLSY79 women's four AFQT subtests (IRT
  z-scores, age-normed by NLS) against their CNLSY children's PIAT math,
  reading recognition, reading comprehension and PPVT at ages 5–14. The
  children's scores are rank-normed within 3-month age cells and averaged over
  rounds. A two-factor model fit by least squares to the off-diagonal
  correlations gives the latent mother–child correlation, loadings and
  composite reliabilities. Sibling pairs give the latent sibling correlation.
- **Wage slope.** Log hourly wage (annual wage income over annual hours) on the
  standardized AFQT composite, on an ACS-style sample: ages 25–34, full-year
  full-time, at least some high school. The slope is divided by the square
  root of the composite's reliability to put it per latent SD. NLSY97 repeats
  it for the cohorts the ACS 2009–13 sample covers.
- **Checks.** Mother–child and sibling correlations of permanent log earnings;
  participation gradients in ability.

Two data traps are handled explicitly. PIAT comprehension is not
administered below a recognition raw score of 19; the file copies the
recognition score instead. Those values stay in the norming reference but are
kept out of the child means. PIAT/PPVT ages come from the child supplement
(`CSAGE`), not the mother supplement.

Writes `estimates/nlsy_ability.json` and `estimates/nlsy_ability.md`. The
report ends with the exact definition of each moment's model counterpart.

## First-stage estimation — `spatial_first_stage.jl` (T2b)

**Complete as a first pass; verified 2026-09-29.** From the repository root:

```bash
julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_first_stage.jl
```

Add `--write` to regenerate `estimates/first_stage.md` and
`estimates/first_stage.toml`. The script combines the ACS occupational shares
and wage 90/10s with the T2a latent moments, holding school quality and taxes
fixed and ignoring teaching selection for both genders. On the model's Nz = 5
grid, it estimates ρz = 0.5766, s_z = 0.1956 and σϵ = 0.7861, matching the pooled
within-occupation wage 90/10 of 3.7374. It also reports grid convergence,
untargeted cell 90/10s and first-stage sensitivity refits.

The rounded estimates are installed in `spatial_calibrate.jl`; the reference
levels in `spatial_reference.toml` have been refrozen after retuning the starting
teaching margin and CFR map. This is not a joint internal fit. Schooling inputs
remain under T2, and full-equilibrium counterparts and consistency checks remain
under T4/T6. See the [first-stage report](estimates/first_stage.md) and
[calibration status and verification](../../julia/spatial_model/calibration/spatial_calibration.md).

## District bootstrap and audit — `district_bootstrap.py`, `district_audit.py`, `acs_earnings.py` (T3, T5)

The 2018 district moments were rebuilt on the SSCC machines on 2026-09-29 and
match the tracked `spatial_moments.json` exactly (68 non-trend estimates):

```bash
python3 fetch_data.py --sources urban seda edge czones acs \
    --urban-topics ccd_directory ccd_finance saipe --years 2018 --seda-levels geodist_pool
python3 data_estimate.py --scheme all --years 2018 --base-year 2018 --outdir <dir> --dump-panel
set -a; source .env; set +a; python3 acs_earnings.py   # ACS 2014-18 earnings tables by school district
python3 district_bootstrap.py                          # about 1 minute
python3 district_audit.py                              # about 4 minutes
```

Keep `--outdir` away from `estimates/` for a 2018-only run: the trend estimates
need the earlier years and would be overwritten. `district_bootstrap.py`
resamples whole commuting zones (B = 1999, seed 20260929) for every targeted and
validation gap and for the earnings alternatives; one resampling matrix is shared,
so `estimates/district_bootstrap_draws.csv` holds joint draws. `district_audit.py`
documents how the salary and FTE gaps are built and recomputes them under
alternatives; its variant 28 (zero-coded wage bills filled, New York City as one
district, K–12 FTE) is the salary target and FTE variant 11 (K–12 teachers per
K–12 pupil) the FTE-per-pupil gap, validation since the exactly identified fit. `acs_earnings.py` adds the earnings tables and the binned
earnings distributions (B20001, B20005) behind the pooled-median earnings gap (validation).
Reports: `estimates/district_bootstrap.md`, `estimates/district_audit.md`.

## NLSY extensions — `nlsy_extensions.py` (T6a, T3)

```bash
python3 fetch_data.py --sources nlsy   # NLSY_BLOCKS includes geography, family income and YA job history
python3 nlsy_extensions.py             # B = 500; --quick, --only t3|t6a
```

T6a regresses CNLSY children's earnings at about 28 on observed scores at ages
8–14, from CFR's definitions to the ACS-style latent slope, with and without
lagged-score controls. T3 measures city–suburb transitions between ages 12–16 and
25–34 in the NLSY97 and NLSY79, by parental-income tercile; the NLSY79 move rate
(0.216) is the calibration's residential-transition target. Report:
`estimates/nlsy_extensions.md`.

## Occupational block from ACS microdata — `acs_occupations.py` (T4, T7)

```bash
set -a; source .env; set +a
python3 acs_occupations.py fetch     # Census API PUMS, cached in raw/acs_pums/ (~100 MB)
python3 acs_occupations.py build     # about 6 minutes; --selftest, --check-bls-table
```

Rebuilds `data/LaborMarketData/wages_occ_shares_v2.xlsx` from Census PUMS microdata
(ACS 2013 5-year file): every group count and wage-sample count matches exactly and
the cell 90/10s to 0.2%, with occupation codes mapped Census → occ1990 → HHJK groups
through the BLS/IPUMS table, the Census 2002→2010→2018 crosswalks and eleven inferred
code departures, and with rules the workbook does not state (listed in the report).
It adds mean hourly and log wages and schooling by occupation × gender (T4) and the
same block on pooled ACS 1-year 2016–19 (T7), which `spatial_calibrate.jl` carries as
`ACS_2016_19` and `spatial_first_stage.jl --t7 --write` refits
(`estimates/first_stage_t7.md`). Report: `estimates/acs_occupations.md`.

## Literature targets

`estimates/literature_targets.md` documents the CFR, BFM and FOO numbers behind
Table 3 (tables, pages, quotes) and the conversion of BFM's willingness to pay to
a consumption share.

## Internal calibration — `spatial_estimate.jl`, `spatial_diagnose.jl` (T6)

From the repository root, on SLURM (`econ-grad` preempts the other partitions):

```bash
E=julia/spatial_model/calibration/estimation
julia --project=$E julia/spatial_model/calibration/spatial_estimate.jl --calibrate-phi       # φ from schooling, after a solver or Table 2 change
julia --project=$E julia/spatial_model/calibration/spatial_estimate.jl --freeze-reference    # C̄, h̄ at ϑ₀
sbatch -c 48 julia/spatial_model/calibration/slurm/estimate.sh --run-dir julia/spatial_model/calibration/runs/<name> \
       --n-samples 470 --n-local 47 --local-max-evals 150 --max-evals 8000     # TikTak
sbatch julia/spatial_model/calibration/slurm/diagnose.sh --theta julia/spatial_model/calibration/runs/<name> \
       --polish --jacobian --consistency --out julia/spatial_model/calibration/runs/<name>/polish   # LM polish and checks
n=$(grep -vc '^#' julia/spatial_model/calibration/slurm/panel_cases.tsv)
sbatch --array=1-$n%10 julia/spatial_model/calibration/slurm/panel.sh julia/spatial_model/calibration/runs/polish-noincome   # Table 5
python3 julia/spatial_model/calibration/slurm/panel_summary.py                               # estimates/sensitivity_panel.md
```

`spatial_diagnose.jl` prints the moment table (validation moments starred), the
validation and §4 consistency diagnostics and the standardized Jacobian at any ϑ
(`--theta theta0`, a TikTak run directory, a `theta.toml`, or comma-separated
values); `--grid Nz=9` re-solves on another grid. Both scripts take `--fix`,
`--external`, `--measure`, `--target key=value`, `--target-se key=value`,
`--activate key` and `--deactivate key`. TikTak's first local search runs alone
before the parallel phase opens, so a run spends its first ~40 minutes on one
worker; the Levenberg–Marquardt polish is the fast finisher. A requeued TikTak job
resumes its journal only if the model files are unchanged (the problem id hashes
them). The baseline is the exactly identified fit in
`julia/spatial_model/calibration/runs/exact-polish-base/` (eight moments, eight parameters; the
FTE-per-pupil gap is validation); the first, overidentified fit is in
`julia/spatial_model/calibration/runs/polish-noincome/`. Both are described in the calibration
log's §4.
