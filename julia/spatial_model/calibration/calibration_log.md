# Spatial calibration log

Companion to [spatial_calibration.md](spatial_calibration.md): target definitions, evidence mappings, validation and next experiments. Tables 1–5 retain their original numbering.

**Current baseline:** the exactly identified 2026-09-29 fit in [`runs/exact-polish-base/`](runs/exact-polish-base/). The October experiments leave it unchanged; their extensions are diagnostic prototypes, not an adopted specification. The log-shock solver correction supersedes all earlier equilibrium results; the outside-equilibrium first stage is unaffected. Verification below refers to saved runs, not new estimation during this documentation update.

**Accepted direction (2026-10-15, after both rounds):** the user provisionally accepted the recommendations, conditional on preserving previous work. Carry M3, estimated wage compression and the S1 supply repair into a small comparison of income-dependent residential sorting and explicit teacher demand. Strengthen residence/commuting acceptance criteria and use finer ability grids. The [accepted development plan](#accepted-development-plan-2026-10-15) supersedes the historical task lists below. Source and selected fit evidence have been checkpointed; no baseline has been replaced and B1 remains sealed.

**2026-10-05 round:** a pre-registered comparison of the 10-03 mechanisms selects **C4′** (residence/workplace split + college choice + priced suburban access; [§4](#model-experiments-2026-10-05)). Its validation score is about a third of the baseline's, but the household-sorting misses remain, so the 09-29 fit stays the reference until a specification is adopted.

**2026-10-07 round:** two further mechanisms ([§4](#model-experiments-2026-10-07)), neither adopted. The regular Roy model ($\gamma=1$) fixes the baseline's teacher premia but removes teachers' selection within occupation and raises their wage dispersion to the market's; on C4′ it worsens the validation score. Local peer effects in school quality create no sorting. Once sorting is imposed, SEDA growth bounds their elasticity at $\psi_a\lesssim0.06$.

**2026-10-15 recommendations (historical, before the rounds):** use C4′ as the development candidate, beginning with an audit of commuting and teacher-supply estimands, then a bounded experiment with heterogeneous housing access. The September fit remains the reference. [Recommended sequence](#recommended-sequence-2026-10-15); the completed rounds and accepted follow-up are recorded separately.

**2026-10-15 rounds (done; [§4](#model-experiments-2026-10-15)).** The sequence above was carried out as a pre-registered round and an informed follow-up. Files: `runs/experiments-2026-10d/` and `-10e/` (untracked); code before the round is frozen in `runs/_code/pre-20261015/`. New mechanisms are switches, off by default. Scores are $S_{1015}$, which replaces the NLSY79 parental-income gradient with its ACS 2014–18 within-CZ counterpart, −0.187. C4′ scores 315.1.

- *Supply.* C4′'s teacher supply is hypersensitive: the small-district elasticity is ≈68, because σw ≈ 0.017. The repair S1 (ρw = 1 plus a utility commuting friction κc fitted to the crossing share) gives a full-equilibrium own-location elasticity of 1.6. Non-teachers' commuting direction stays reversed in every fit.
- *Staffing closures.* Under S1, a salary-matched δκ over-staffs the suburb (FTE +0.064). An FTE-matched δκ cuts suburban pay (−0.040) and removes the score gap. A non-pecuniary closure, R, matches salary and FTE together: teachers' value of a suburban job, $b_T$ (the existing `Params.BT`), fitted to FTE.
- *Sorting.* Four forces create almost no household sorting: θν, bounded housing access, persistence (ρz up to 0.90) and peer effects. School quality is worth only ≈0.002 utils to parents. A utility moving cost (U) flattens the own-AFQT mobility gradient. Only the BA-graded suburban amenity ζ_BA (E) reaches the parental-income gradient V6a.
- *Best candidate.* S1·R·E (S1 + R + E) scores 105 against C4′'s 315.
  - It matches V6a, V6z (−0.225 against −0.222), poverty (−1.1 SE), staffing, salary and supply, and passes A1–A4 and A6.
  - It fails A5's grid check: re-solved at Nz = 9 with the same parameters, the score gap moves 0.32 SE. Re-fitted at
    Nz = 9, its moments are unchanged. **C4′ therefore stays the development candidate.**
  - It is a reduced form. ζ_BA ≈ 1.39 utils is weakly identified, because V6a saturates. Earnings and scores overshoot (+5.3 and +3.4 SE). Teachers' residence and net-commuting gaps are 5 and 20 times the data.
- *Hold-out.* B1, the BA share by residence, stays sealed until the user instructs otherwise (amendment 5). S1·R·E predicts +0.57.

- [External inputs and first stage](#1-external-inputs)
- [Internal targets and validation](#2-internal-calibration)
- [Evidence-to-model mappings](#3-mapping-evidence-to-model-objects)
- [Estimation and experiments](#4-estimation-and-sensitivity)
- [Next experiments and remaining work](#5-to-do)

## Known issues at the current fit

| Issue | Evidence and implication |
|---|---|
| Household sorting | Earnings gap −0.010 versus 0.112; score gap 0.061 versus 0.289 SD; poverty +0.2 versus −8.9 pp. Tested parameter and weighting changes do not resolve these misses. C4′ reaches only 0.034, 0.082 and −2.3 pp. The own-AFQT sorting miss replicates in NLSY97 (2026-10-05). Local peer effects in school quality add no sorting (2026-10-07). Persistence does not help either: raising ρz to 0.90 moves V6a by −0.004, because school quality is worth only ≈0.002 utils to parents against σν = 0.25. A BA-graded suburban amenity reaches V6a only by also moving teachers (2026-10-15). With the non-pecuniary staffing closure R it reaches V6a, V6z and poverty with staffing matched (S1·R·E, $S_{1015}$ 105). But it does so as a weakly identified reduced form, with ζ_BA ≈ 1.4 utils. It overshoots earnings and scores, puts teachers' residences in the suburb at five times the data gap, and fails the Nz = 9 grid check (2026-10-15b). |
| Teacher allocation | Salary fit implies FTE/pupil +0.023 versus −0.018 and SEDA growth 0.012 versus 0.003 SD/grade. With a residence/workplace split and positively selected teachers (C4′): +0.020 and 0.013. The split's staffing-ratio response to suburban pay is ≈31. Audited on 2026-10-15: C4′'s workplace nest is hypersensitive, with a small-district elasticity ≈68 against 1.4–3.7 in the literature. The repair S1 (ρw = 1 plus a commuting friction) gives a full-GE own-location elasticity of 1.6. With realistic supply, however, a salary-matched δκ over-staffs the suburb (FTE +0.064), and an FTE-matched δκ cuts suburban pay (−0.040) and removes the suburb's school-quality edge. The model lacks a demand side or position rationing ([§4](#model-experiments-2026-10-15)). A non-pecuniary workplace taste $b_T$ (closure R) matches FTE and salary together. On its own, like SC, it removes the suburb's school-quality edge (2026-10-15b). |
| Teacher education and selection | Baseline: schooling ratio 0.76 versus 1.18; teacher pay well below market pay; teachers' ability −0.29 SD versus +0.61. The college choice with licensure (M3) gives 1.14, premia −0.10/+0.01 versus +0.06/+0.09, and ability +0.55 SD (−0.29 versus BA non-teachers, data −0.33), without a schooling floor or teacher utility costs. Still missed: MA premium 0.07 versus 0.16; non-BA teachers' schooling 0.37 versus 0.52. The regular Roy model ($\gamma=1$) brings baseline teacher pay to market pay (premia −0.09/+0.01) but removes selection (0.00 SD) and raises teachers' 90/10 to 3.75 versus 2.77. Under M3 it overpays teachers (+0.25/+0.37) (2026-10-07). |
| Occupational wages | Cell wage levels and dispersion are missed. The market gender log-wage gap is 0.30 versus 0.10, and 0.36 under M3. |
| Family inputs and participation | Latent sibling correlation 0.33 versus 0.606; non-teaching participation is independent of inherited ability, unlike the data. Parents finance none of their children's investment. A family component shared by siblings would fix the sibling check, but it would not create sorting (2026-10-15 persistence test). |
| Measurement and inference | ACS/NLSY samples differ; public mobility data cannot isolate within-MSA moves; BFM score units use an assumed ICC; poverty uses a provisional threshold. Teacher payroll is 2.5% of model market income (4.5% under M3), versus ≈1.6% in the data; the comparator is not like-for-like (age and income definitions). Full calibration bootstrap and literature transport uncertainty remain open. |

The [model](../spatial_continuous.jl) ([writeup](../../../notes/spatial_teachers_model_eps0.md)) represents city and suburb within a CZ. Spatial gaps are enrollment-weighted within-CZ contrasts, suburb minus city, in logs unless stated otherwise. District inputs are from 2018, SEDA pools 2009–19, and occupations use ACS 2009–13. Revenue and SES partitions are robustness checks. Residence, workplace and school are currently bundled; housing and teacher effort are omitted, and home production is taxed.

Code: `spatial_calibrate.jl` (inputs and targets), `spatial_first_stage.jl` (Table 2), `spatial_moments.jl` (counterparts), `spatial_estimate.jl` (search), `spatial_diagnose.jl` (polish, Jacobian and checks), `spatial_cli.jl` (options); batch scripts in `slurm/`.

## 1. External inputs

### Table 1: Fixed inputs

Asterisks mark provisional values not yet rebuilt from the 2018 sample.

| Parameter | Value | Basis | Source |
|---|---|---|---|
| Ability elasticity $\alpha$ | 1 | Normalization; only $\alpha\sigma_\epsilon$, $\alpha s_z$ identified | — |
| Population $M$, amenity $B_1$, home-production productivity $A_{HP}$ | 2, 0, 1 | Normalizations | — |
| Consumption weight $\mu$ | 1 | Normalization (below) | — |
| Goods-investment elasticity $\eta$ | 0.080 | $(S-t)/(1-t)$, $S=0.102$, $t\approx0.023$ | [NCES 605.20](https://nces.ed.gov/programs/digest/d21/tables/dt21_605.20.asp); [PWT 11](https://fred.stlouisfed.org/series/LABSHPUSA156NRUG); [NCES 211.60](https://nces.ed.gov/programs/digest/d18/tables/dt18_211.60.asp), [236.20](https://nces.ed.gov/programs/digest/d20/tables/dt20_236.20.asp) |
| Time-investment elasticity $\phi$ | 1.0584 | Non-teacher time share $s_O=13.64/25=0.5457$ at the reference $\Xi$ | ACS 2009–13 (`Mincer_YrsSchool.xlsx`); HHJK (2019) mapping |
| Class-size curvature $\sigma$ | 0.56* | Class-size wage effect | [FOO (2013)](https://pure.uva.nl/ws/files/1912854/124460_394782.pdf#page=26) |
| Location taste scale $\sigma_\nu/\mu$ | 0.25* | Literature | [Eckert–Kleineberg (2024)](https://www.fpeckert.me/papers/ek_2024.pdf#page=17) |
| Warm-glow curvature $\psi$ | 0 | Linear warm glow | — |
| Earnings–achievement bridge $\chi$ | 0.113 | $\log 1.12$; converts school-driven score changes only, not $s_z$ (section 3) | [CFR (2014b)](https://opportunityinsights.org/wp-content/uploads/2018/03/teachers2.pdf) |
| School years of exposure $K$ | 12 | Grades 1–12, equal per-year exposure | Model timing |

**Goods investment.** Spending over labor income is $S=\eta(1-t)+t$, where $t$ funds teacher payroll, so $\eta=(S-t)/(1-t)$. For 2018, $S=0.060/0.591=0.102$ and $t\approx0.023$ (3.2 million FTE teachers at \$60,483, benefits at 0.43 of salary, against GDP of \$20.66 trillion), giving $\eta=0.0807$ (`eta_2018`). Using gross $S$ would count teacher pay twice. Counting all instructional staff instead gives $\eta=0.073$. Model tax rates also include home-production output in the base, so they are not directly comparable to the market-income payroll share.

**Time investment.** With zero goods moving costs, $s_O=\mu\phi/(\mu\phi+1-\eta)$ and $s_T=\mu\phi\gamma/(\mu\phi\gamma+1-\gamma\eta)$, which depend only on $\mu\phi$. Location choice and WTP depend only on ratios to $\mu$, and $\bar\kappa$ absorbs its remaining role, so $\mu=1$ and $\phi$ matches non-teachers' schooling ($s=$ years/25). With $r_m>0$, evaluate at the reference $\Xi$.

**Schooling (T2).** From `data/Mincer_YrsSchool.xlsx` (sheet `years_of_school`, 2013 column: ACS 2009–13, ages 25–34, years capped at 9–17), averaged over occupation-by-gender cells with the Table 2 counts: non-teachers, home production included, 13.64 years (men 13.44, women 13.84) and K–12 teachers 16.04 (both genders), so $s_O=0.5457$ and $s_T/s_O=1.176$. The $\Xi=0$ formula gives $\phi=1.105$; the goods moving cost raises time shares through $\Xi$, and `calibrate_phi` sets $\phi=1.0584$ so that non-teachers' mean $s_O$ at the reference solution is 0.5457 (C̄, h̄ iterated jointly). The ACS microdata pull (T4) cross-checks the workbook.

**Known miss.** A common $\phi$ with $\gamma<1$ gives teachers less schooling: $s_T/s_O=0.76$ at the baseline versus 1.18 in the data. The implied leisure advantage lets the model match teaching shares at low pay. Non-teacher schooling is close to target ($s_O=0.550$ versus 0.546).

**Exposure.** Per-year evidence pins $\sigma/K$ and $\beta/K$, while $\beta/\sigma$ is free of $K$, so $K$ sets how much schools matter. $K=5$ (SEDA's tested grades) is a sensitivity case.

**Location tastes.** Without housing, effective dispersion is roughly $\sigma_\nu/\mu$ plus housing share over supply elasticity ($\approx0.45$), so 0.25 likely overstates population responses; sensitivity extends to 0.50.

**Units.** Moving costs and warm glow use the fixed reference in `spatial_reference.toml`: $\bar C=0.044114$, $\bar h=0.095780$, set after the 2026-09-29 solver correction and schooling calibration. Keep these fixed across comparative statics; rebuilding the reference requires a consistent external block and $\phi$.

### Table 2: Occupations and ability, estimated outside the equilibrium

Targets come from the ACS 2009–13 sample in `data/LaborMarketData/wages_occ_shares_v2.xlsx`: ages 25–34 with at least one year of high school, military and unemployed dropped. Home production is out of the labor force or under 15 hours a week; 15–29 hours are split equally between home production and the reported occupation. Wages are hourly, for workers employed at least 48 weeks/year and 30 hours/week, with annual wage-and-salary plus business/farm income of at least \$1,000 in 2010 dollars (workbook `readme`, cells A3:B3). No percentile trimming is specified in the workbook. The model has 21 occupations: the file's 19 market occupations, home production and K–12 teaching (postsecondary and other teachers are market occupations).

| Parameter | Target | Value | Estimate (T2b) | Source |
|---|---|---|---|---|
| Relative productivities $A_i/A_{HP}$, 19 market occupations | Men's shares over the 20 non-teaching occupations | `moments_shares` | $\log A_i/A_{HP}\in[-1.53,-0.33]$ | ACS 2009–13 |
| Relative female wedges $\tau^\omega_{i,f}$, 19 free | Women's shares over the 20 non-teaching occupations | `moments_shares` | $\log r_i\in[-1.58,0.44]$ | ACS 2009–13 |
| Ability dispersion $\sigma_\epsilon$ | Weighted within-occupation 90/10 of hourly wages, non-teachers | 3.74 | 0.786 (Nz = 5) | ACS 2009–13 |
| Ability persistence $\rho_z$ | Mother–child correlation of latent test scores (net of the $Q_l$ channel) | 0.577 (SE 0.014); 0.53–0.63 across samples | 0.577 | NLSY79 AFQT, CNLSY PIAT/PPVT |
| Stationary dispersion $s_z$ | Slope of log wage on standardized latent AFQT, $\approx s_z/(1-\eta)$ | 0.213 (SE 0.008) | 0.196 (SE 0.007) | NLSY79; NLSY97 0.199 |
| *Check only* | Occupation-specific 90/10s | `90_10_hr_wages_weighted` | Missed (below) | ACS 2009–13 |
| *Check only* | Mean log wages by occupation, and teachers' | ACS 2009–13 PUMS (`acs_occupations.py`) | Missed (§4, *Mean wages*) | ACS 2009–13 |
| *Check only* | Mother–child and sibling correlations: earnings 0.135 and 0.281; latent scores 0.606 (siblings) | `nlsy_ability.md` | $\rho_z^2=0.33$ (siblings) | NLSY79, CNLSY |

Forty-one parameters, forty-one moments. No education wedges and no wage wedges for men or teaching. The levels of $A_i$ and of the female outside option relative to teaching are internal (Table 3).

**Ability block.** The first stage fixes $Q_l,t_l$ and ignores teaching selection consistently in shares and wage distributions. Non-teacher log wages are $b[\alpha\log z+\log(\Theta_{i,g}\epsilon_i^\alpha)]$ plus a constant, with $b=1/(1-\eta)$. Latent mother–child correlation pins $\rho_z=0.577$; the latent wage slope pins $s_z=0.21266(1-\eta)/\alpha=0.1956$; the pooled wage 90/10 then pins $\sigma_\epsilon$. Scores measure $z$, not occupational shocks, separating the two dispersions. Composite reliabilities (0.92 mothers, 0.89 children) matter for observed-score checks. Read $\rho_z$ as reduced form, including assortative mating; the model has no spouses. The sibling correlation $\rho_z^2=0.33$ versus 0.606 remains a miss. Section 4 checks teaching selection and spatial variation in the full equilibrium.

**Ability grid.** Use the same Rouwenhorst grid in the first stage and equilibrium. It matches the variance and autocorrelation of $\log z$ exactly; quantiles still require convergence checks. Estimates of $\sigma_\epsilon$ are 0.7861, 0.7866 and 0.7867 at Nz = 5, 9 and 15, versus 0.7870 under continuous $z$. Nz = 5 stays within 0.1% of the continuous pooled 90/10; equilibrium grid checks are in section 4.

**Scale and $\chi$ (decided).** $s_z$ targets the NLSY latent-score/log-wage slope; $\chi$ converts school-driven achievement changes only. The score loading is $c=b s_z/\chi\approx1.88$. The [CNLSY replication](../../../data/spatial/estimates/nlsy_extensions.md) supports distinguishing the slopes: earnings at 28 relative to the mean rise 0.230 per observed score SD without lagged-score controls and 0.118 with lagged-score cubics and parental controls, close to CFR's conditional 12%. The alternative $b s_z=\chi$ implies $s_z=0.104$, about 13 SEs below the NLSY estimate, and is a sensitivity case. Neither association identifies the causal return to inherited ability.

**NLSY moments (T2a).** From `nlsy_ability.py`; rerun with 200 bootstrap draws (seed 20260928), clustered over NLSY79 households and NLSY97 persons, with no failed draws. Full results: [NLSY report](../../../data/spatial/estimates/nlsy_ability.md).

**Wage sample.** Both NLSY cohorts require 48 weeks/year, 30 hours/week and \$1,000 annual wage-and-salary income in 2010 dollars, with no percentile trimming. Deflate using the prior calendar year (the income reference year); year dummies absorb deflation in log-wage regressions. The NLSY79 slope is 0.213 on 37,870 person-years (7,219 people; 95% CI [0.200, 0.229]); 1/99 trimming gives 0.196 (SE 0.007). NLSY97 gives 0.199 (SE 0.009).

These are **ACS-matched cutoffs, not identical samples**. The ACS workbook includes business/farm income; the NLSY extracts contain wages and salary only. NLSY uses age at interview and prior-calendar-year income/hours/weeks, with weekly hours constructed as annual hours divided by weeks worked. Existing status and education filters also differ: NLSY79 excludes armed forces (ESR = 4) but does not drop current unemployment or missing ESR; NLSY97 has no military-status filter in this extract and retains missing highest grade. These differences remain explicit limitations of the comparison.

- *Persistence:* the latent mother–child correlation is 0.577. Sample choices move it more than sampling error (0.53 with cross-sectional mothers, 0.63 unweighted), so both variants go to Table 5. Sorting on $Q_l$ raises the score correlation, so expect $\rho_z$ at or somewhat below 0.58. The sibling ratio $\mathrm{cov}_{sib}/\mathrm{cov}_{pc}=1.05$ is a check only (siblings share inputs outside the model).
- *Scale:* the log-wage slope is 0.213 per latent SD (men 0.190, women 0.247), so $s_z\approx0.1955$ (SE 0.007); NLSY97, the ACS cohorts, gives 0.199, so vintage does not close the gap to $\chi$. The NLSY79 slope exceeds $\chi$ by 0.099 (SE 0.008).
- *Selection into work:* each observed composite SD raises full-year full-time wage work (including the real-income floor) by 0.08 and lowers working under 15 hours by 0.06 (0.08 for women). The model makes participation independent of $z$, so mothers' earnings are a selected measure.

**Pooled 90/10.** The target averages cell ratios with occupation-by-gender wage-sample weights (`occ_gender_weights`; home production has zero weight); build the model counterpart the same way. Age-residualized, it is 3.64.

**First-stage results (T2b).** `spatial_first_stage.jl` matches pooled 90/10 = 3.737 at $\hat\sigma_\epsilon=0.786$. Shares depend on $\log\Theta/(\alpha\sigma_\epsilon)$, so the inversion is solved once and rescaled. Roy selection is essential: an unselected-lognormal calculation gives only 0.431. Sensitivity estimates range from 0.766 (age-residualized 90/10) to 0.844 ($b s_z=\chi$). Full estimates and refits: [first-stage report](../../../data/spatial/estimates/first_stage.md).

**First-stage verification (2026-09-29).** The script reproduces the saved tables, workbook counts, 38 wage ratios and weights, and NLSY constants. Independent Python integration gives $\hat\sigma_\epsilon=0.786108$ versus 0.786107 in [first_stage.toml](../../../data/spatial/estimates/first_stage.toml), with shares agreeing to $4.4\times10^{-13}$. `EXTERNAL` uses four-decimal estimates. These checks validate the first stage; full-equilibrium consistency is reported separately.

*Cell 90/10s (check only) are missed.* Common $\sigma_\epsilon$ gives cell ratios [3.21, 4.15] versus [2.80, 6.59], with weighted correlation −0.19 across 38 cells. The model gives men/women 3.65/3.85 versus 3.82/3.63. Occupation-specific dispersion would need additional identification; matching the pooled target does not validate the cells.

**Why no equilibrium.** Among non-teaching alternatives, agents maximize $(1-\tau^\omega_{i,g})A_i\epsilon_i^\alpha$; $z$, $Q_l$ and spatial terms cancel. Conditional non-teaching shares therefore identify this Roy block when teaching selection is ignored. Small aggregate teaching shares do not guarantee small errors within occupations or wage tails; recheck after each substantive model change.

**Caveats.**

1. If $z$ alone gives a 90/10 above 3.74, $\sigma_\epsilon$ is at a corner. Fix $\eta$ first, since wages load on ability through $\alpha/(1-\eta)$.
2. Measurement error and transitory shocks in ACS wages bias $\sigma_\epsilon$ up given $s_z$.
3. Table 2 and the teaching margin use ACS 2009–13, against late-2010s spatial targets. Wage dispersion was rising (90/10 of 3.36 in 2000), and men's home-production share (22.9%, against 18.6% in 2000) reflects the recession, lowering market $A_i$ and the levels set against them, including $\bar\kappa$. T7 measures it: on 2016–19 the pooled 90/10 is 3.82 and $\hat\sigma_\epsilon$ 0.806, market $A_i$ rise 5% against home production, and the internal fit barely moves.
4. Mean wages and cell-level 90/10s are untargeted checks on a common $\sigma_\epsilon$.

## 2. Internal calibration

Write $\kappa_1=\bar\kappa$, $\kappa_2=\bar\kappa e^{\delta_\kappa}$, $B=(0,\Delta B)$, $m_{12}=m_{21}=r_m\bar C$ and $1-\tau^\omega_{i,f}=\omega_f r_i$. The pay level enters as $\bar\kappa=\tilde\kappa\,\bar h^{1-\gamma}$, so a teacher with the reference human capital $\bar h$ earns $\tilde\kappa\bar h$ in the city at any $\gamma$; without this normalization the admissible $\bar\kappa$ shifts with $\gamma$ (since $h<1$), a diagonal valley for the search. Fit all parameters jointly; each row names the main source of identification. District moments are from [spatial_moments.md](../../../data/spatial/estimates/spatial_moments.md).

### Table 3: Internal parameters and targets

| Parameter | Target | Value | Source |
|---|---|---|---|
| ***Teaching margin*** | | | |
| Teacher pay scale $\tilde\kappa$ | Male teaching share | 1.9% | ACS 2009–13 |
| Female outside-option level $\omega_f$ | Female teaching share | 6.0% | ACS 2009–13 |
| Teacher wage curvature $\gamma$ | Teachers' 90/10 of hourly wages | 2.77 | ACS 2009–13 |
| ***Spatial block*** | | | |
| Spatial pay gap $\delta_\kappa$ | Adjusted salary log gap | 0.0086 (SE 0.0092) | District, 188 CZs (audited) |
| Amenity gap $\Delta B$ | Enrollment log gap | 0.3764 (SE 0.1154) | District, 188 CZs |
| Teacher spillover $\beta$ | Earnings at 28 per teacher-quality SD-year | 1.34% (SE 0.41 pp) | CFR (2014b), Table 3 col. 2 |
| Altruism $\lambda$ | WTP per SD of school-average scores, share of consumption | 0.595% (SE 0.224%): \$19.70 (SE \$7.40) over \$3,309/month | BFM (2007), Table 7; CE 1990–91 |
| Goods moving cost $r_m=m/\bar C$ | City↔suburb move rate, ages 14–16 to 25–34 | 0.216 (SE 0.011; scale 0.05) | NLSY79 |

Eight parameters, eight moments: scales affect the search, not the exact fit. FTE/pupil was the first fit's ninth target and is now validation; `--activate gap_teachers_pp` restores it. Earnings barely respond to $r_m$, so mobility identifies moving costs and earnings remain validation; `--activate gap_income` adds them to the criterion.

**Target definitions.**

- **Teaching shares:** K–12 teachers over all 25–34-year-olds of each gender, home production included in the denominator (Table 2 sample). Scale: 10% of the value.
- **Teacher 90/10:** wage-sample-weighted average of men's (2.597) and women's (2.830) ratios (weights 0.2444 and 0.7556, the K–12 wage-sample counts 12,118 and 37,465), 2.773. The model counterpart pools teachers across locations by gender (pre-tax $(1-\tau^\omega_{T,g})\kappa_{l'}h_T^\gamma$) and averages the two ratios with the same weights (`teacher_wage_stats`). Scale: 5% of the value.
- **District scales (T5):** whole-CZ bootstrap SEs, B = 1999, seed 20260929, holding the screen, partition and filters fixed ([bootstrap report](../../../data/spatial/estimates/district_bootstrap.md), `district_bootstrap.py`). They are 1.2–1.8 times the old SD/√n placeholders: the pupil weights leave about 42 effective CZs of 188 (the top 5 carry 26%). Target draws are correlated (pupils–income 0.53, salary–income 0.32, FTE–pupils −0.26), which the diagonal weighting ignores. The 2018 district moments were rebuilt on this machine from fresh downloads and reproduce the tracked `spatial_moments.json` exactly (68 non-trend estimates; the trends need earlier years).
- **Salary gap:** CWIFT-adjusted; matches $\Delta\log[\kappa_l\,\mathbb E(h_T^\gamma\mid T,l)]$, which includes teacher composition. Prefer standardized schedules where available. Audit ([district audit](../../../data/spatial/estimates/district_audit.md)): the published 0.0162 covers 176 CZs, losing New York (F-33 files the city under one agency) and every CZ in Illinois, New Mexico and Alaska (zero-coded teacher-salary fields), and counts pre-K teachers in the FTE denominator while their pay is mostly outside the four F-33 teacher-salary fields (within-CZ elasticity of pay per FTE to the pre-K share −0.73, SE 0.28). The target fills the zero-coded wage bills with 0.881 × instruction salaries, treats New York City as one district and divides by K–12 FTE: 0.0086 (SE 0.0092) on 188 CZs; the published value is a Table 5 case. The sign comes from the CWIFT deflator alone (undeflated 0.0003); across trimming, deflator and sample variants the gap runs from 0.004 to 0.026, so read it as small and imprecise.
- **Enrollment:** equals the population-share gap with equal cohorts.
- **WTP:** BFM's \$19.70 is monthly housing user cost (rent; house value/264.1 for owners), 1990 Census, mean over all households, per SD of school-average 1992–93 CLAS grade-4 scores (Table 7 panel B col. 4; also Table 8 row 1). Over mean monthly consumption in the San Francisco area (CE 1990–91, \$39,707 a year) it is 0.595% (SE 0.224%); US (\$2,365) and West (\$2,705) denominators give 0.833% and 0.728% (Table 5). One school-average SD is $\sqrt{\mathrm{ICC}}\approx0.45$ student SDs (grade-4 ICCs near 0.2; range 0.41–0.49; BFM report no student-level SD), `Measurement().wtp_score_sd`. Parents: BFM's Table 8 puts households with children \$7.41 (SE 3.58) above the mean, about \$24.6 (a Table 5 case). Sources and arithmetic: [literature targets](../../../data/spatial/estimates/literature_targets.md).
- **CFR:** Table 3 col. 2, \$285.55 (SE 87.64) on mean earnings \$21,256 at 28, 0.01334 (SE 0.00407) in logs, with CFR's SE as the scale. Col. 1 (1.65%, SE 0.43 pp) is a Table 5 case.
- **Income:** the model counterpart is the gap in log median pre-tax earnings of market workers, teachers included and home producers excluded (`income = :log_median_market`); the model has no hours margin, so the data counterpart is full-time year-round earnings. It is the log of each location's median over its pooled ACS 2014–18 earnings distribution (B20005 bins summed over the location's districts, linear interpolation within bins, which reproduces published district medians to a mean log error of 0.008): 0.1122 (SE 0.0168). No school-district table cuts earnings by age, so this is all workers 16+. The household-income proxy it replaces (0.2708) is twice every earnings-based gap: all earners 0.1435, men 0.194, women 0.091, the mean of district log medians 0.158; household income is 0.249 even for householders 25–44, so the wedge is households versus individuals, not age. Both are Table 5 cases. The gap does not identify $r_m$ (§4, *First internal fit*), so it is validation.
- **FTE per pupil:** matches teacher headcount per student $n_l/S_l$. K–12 teachers per K–12 pupil (pre-K teachers and pupils removed, consistent with the salary gap): −0.0177 (SE 0.0131), against −0.0122 with pre-K. The sign is carried by the largest CZs (unweighted CZ mean −0.048). Validation (Table 4): it was the one overidentifying moment, and no parameter vector found matches it together with the salary gap (§4, *Exactly identified fit*).
- **Move rate:** share of youths living in an MSA central city (suburb: in an MSA, outside the central city) at the 1979 interview, ages 14–16, who live in the other type in the rounds at 25–34, among those in an MSA at both ages; NLSY79, weighted, 1,498 persons, 0.216 (SE 0.011) (`nlsy_extensions.py`, [report](../../../data/spatial/estimates/nlsy_extensions.md)). Variants span 0.168 (adult years 1987–90) to 0.314 (NLSY97, whose origin and adult rounds straddle a change in MSA standards, an upper bound), hence the scale of 0.05. Cross-MSA moves cannot be separated in the public files. The model counterpart is the share of the young schooled in one location who live in the other, $\sum_l S_l(1-P_{ll})/\sum_l S_l$. By parental-income tercile the rate is flat (NLSY79 0.221, 0.226, 0.207; NLSY97 0.306, 0.318, 0.333), a check on $r_m$, which makes moving costlier for poorer children (`move_by_parent_income`).

### Table 4: Untargeted validation

| Outcome | Gap (data) | Model counterpart |
|---|---|---|
| Earnings gap (full-time year-round median) | 0.1122 (SE 0.0168) | Δ log median pre-tax earnings of market workers (`gap_income`, `income = :log_median_market`) |
| K–12 teacher FTE per K–12 pupil | −0.0177 (SE 0.0131) | Δ log teacher headcount per student, $n_l/S_l$ (`gap_teachers_pp`); a Table 3 target in the first fit |
| Child poverty | −8.88 pp (SE 0.90) | Share of children whose parent's pre-tax income (home production included) is below the economy-wide 18th percentile (`poverty_gap`; the line is provisional) |
| SEDA learning rate | 0.0031 SD/grade (SE 0.0030) | $\Delta\log Q/(\chi K(1-\eta))$ |
| SEDA score level | 0.2887 SD (SE 0.0327) | $\Delta\mathbb E[\log z+(g/K)c\log Q_l]$ at $g=5.5$ over the young, in cohort SDs of that signal (`score_gap`); composition plus $Q$, descriptive |
| Local revenue share | 0.4726 (SE 0.0206) | Diagnostic only |
| Move rate by parental-income tercile | NLSY79 0.221, 0.226, 0.207 (flat) | `move_by_parent_income` |

## 3. Mapping evidence to model objects

**Teacher quality (CFR).** With teacher headcount $n_l$ and student mass $S_l=M_l/2$,

$$
Q_l=\left[\frac{n_l}{S_l}\,\mathbb E(h_T^{\beta/\sigma}\mid T,l)\right]^\sigma .
$$

The efficient class-size allocation leaves no equilibrium variation in teacher value added, so match CFR with an auxiliary intervention: raise one teacher's log human capital by one SD for one year, holding class size and location fixed. For a non-teaching child,

$$
\Delta\log y'\simeq\frac{\beta\,s_{\log h_T}}{K(1-\eta)}\;\longrightarrow\;\log(1.013).
$$

This assumes one SD of teacher human capital equals one SD of causal effectiveness.

**Class size.** FOO (2013, Table V, col. 1) find log wages fall 0.0063 (SE 0.0033) per extra pupil over grades 4–6. With three years of exposure and mean class size 24.357,

$$
\sigma\simeq(1-\eta)\frac{K}{3}(24.357)(0.0063)=0.565
\quad\text{at }\eta=0.080,\ K=12,
$$

with a 95% interval of $[-0.015,\,1.144]$; it is 0.235 at $K=5$. Proposition 1 requires $\sigma<1$. Cross-checks: Krueger's STAR kindergarten comparison gives 0.687, and [JJP (2016)](https://www.nber.org/papers/w20847)'s spending elasticity of about 0.7 is close to the implied $\sigma/(1-\eta)\simeq0.61$. Recompute $\sigma$ whenever $K$ or $\eta$ changes.

**Earnings–achievement bridge.** $\Delta\log y_{\rm adult}=\chi\,\Delta a$, with $a$ in student test-score SDs. CFR (2014b) report about 12% higher earnings at 28 per score SD, a conditional association; $\chi=\log1.12$. It is used only for school-driven score changes (SEDA, WTP), read as $b\,s_z/c$; the ability scale $s_z$ comes from the NLSY (section 1, *Scale and $\chi$*).

**SEDA.** The weighted learning-rate gap is 0.0031 SD/grade (SE 0.0030) and the score-level gap 0.2887 SD; with CZ effects, 0.0046 (0.0022) raw and −0.0042 (0.0022) with SES controls. Compare via $\Delta g\simeq\Delta\log Q/(\chi K(1-\eta))$. Growth also reflects families and peers, so the gap is neither a causal teacher effect nor an upper bound on one.

**Willingness to pay.** Compute compensating consumption $a$ from

$$
\mu\log C+\lambda\,\mathbb E[f(h'_0)]=\mu\log(C-a)+\lambda\,\mathbb E[f(h'_1)],
$$

converting BFM's school-average SD (74 points) to student-level SDs before applying $\chi$. Use the \$19.70 structural estimate, not the \$17.30 hedonic one (Table 3 panel B col. 4, SE 5.9). The test is the 1992–93 CLAS, not the CAP, and the 74 points is the SD over BFM's full sample of 242,100 households, not the boundary sample behind the \$19.70. With no student-level SD reported, one school-average SD is $\sqrt{\mathrm{ICC}}\approx0.45$ student SDs, so $x=0.45\chi$ (`wtp_score_sd`). The model computes $a$ for every parent and divides mean $a$ by mean consumption (`school_quality_wtp`).

## 4. Estimation and sensitivity

Minimize $[\mathcal M(\vartheta)-\widehat{\mathcal M}]'W[\mathcal M(\vartheta)-\widehat{\mathcal M}]$. The baseline uses diagonal inverse fit-scale weights; correlated district weights are available for overidentified cases. With eight parameters and eight matched moments, weights affect the search but not the root. Full bootstrap inference, including Table 2 re-estimation and literature transport uncertainty, remains to be done.

- **Solver correction (2026-09-29).** Quantile-spaced shock grids with level splines under-resolved the teacher tail after $\sigma_\epsilon$ increased. Grids, splines and quadrature now use log shocks. All earlier equilibrium results, including `runs/fit1-quantile-grid-invalid`, are superseded. The first stage is unaffected; $\phi$ and the reference units were recalibrated.
- **Solution checks.** Verify equilibrium and household residuals, consumption support, nonnegative probability masses, stationarity, multiple starts and finer ability/shock grids. `spatial_diagnose.jl --grid Nz=9` (or `nXO=96`) re-solves with the reference fixed. Nonmonotonicity warnings occurred at distant trial points, not near the baseline.
- **External-block consistency.** Recompute occupation shares, cell and pooled 90/10s, wage–score slope and mother–child score correlation with teaching selection, $Q_l$ and taxes operating. Re-estimate Table 2 if distortions are material. A selection correction depending on $z$ must use the selected joint distribution in both shares and wages.
- **Measurement implementation.** `micro_sample` integrates adult choices on 600 log-shock midpoints; 300–1,200 midpoint comparisons move quantile gaps by less than $10^{-5}$. Earlier Monte Carlo checks agreed with the counterparts but preceded the solver correction. Score exposure is $(g/K)c\log Q_l$: adults $g=K$, CNLSY children $g=4$, SEDA $g=5.5$. `--consistency` reports the first-stage comparisons.

### Exactly identified fit (2026-09-29): the baseline

Eight Table 3 targets match at criterion $1.6\times10^{-19}$ in [`runs/exact-polish-base/`](runs/exact-polish-base/). After a 512-point Sobol screen, Levenberg–Marquardt from 15 dispersed starts converges to the same root (`runs/exact-ms/`); no alternative root was found. FTE and earnings remain validation. Distant search points produced nonmonotonicity warnings; the final polish and grid checks did not.

| Parameter | $\log\tilde\kappa$ | $\delta_\kappa$ | $\Delta B$ | $\beta$ | $\lambda$ | $r_m$ | $\omega_f$ | $\gamma$ |
|---|---|---|---|---|---|---|---|---|
| Estimate | −1.2038 ($\bar\kappa=0.123$) | 0.0045 | 0.0662 | 0.2278 | 0.1245 | 0.206 | 0.9664 | 0.6194 |
| First fit (nine moments) | −1.2080 ($\bar\kappa=0.120$) | −0.0018 | 0.0557 | 0.2240 | 0.1236 | 0.248 | 0.9682 | 0.6107 |

All eight targets are matched to at least six significant digits (teachers' 90/10 2.773: men 2.65, women 2.81). The equilibrium local tax rates are $t_1=0.016346$ (city, 1.635%) and $t_2=0.017115$ (suburb, 1.712%); their definition and experimental comparisons are below.

| Validation moment | Data | Model | $(m-t)/se$ | First fit |
|---|---|---|---|---|
| FTE per pupil gap | −0.0177 (SE 0.0131) | +0.0234 | +3.14 | +0.0048 (+1.72, targeted) |
| Earnings gap | 0.1122 (SE 0.0168) | −0.0102 | −7.28 | −0.0152 (−7.6) |
| SEDA learning-rate gap | 0.0031 SD/grade (SE 0.0030) | 0.0119 | +2.93 | 0.0021 (−0.3) |
| SEDA score-level gap | 0.2887 SD (SE 0.0327) | 0.061 | −7.0 | 0.007 |
| Child-poverty gap | −8.88 pp (SE 0.90) | +0.24 pp | | +0.52 pp |
| Move rate by parental-income tercile | 0.221, 0.226, 0.207 | 0.205, 0.216, 0.226 | | 0.170, 0.182, 0.191 |

- **Pay–staffing tradeoff.** Dropping the FTE target lets the salary gap match exactly, but raises FTE/pupil to +0.023 versus −0.018. The overidentified first fit split the discrepancy between salary (−1.23 SE) and FTE (+1.72 SE).
- **Sorting still fails.** The earnings gap is −0.010, the poverty gap +0.2 pp and the score-level gap 0.061 SD against 0.289. The larger suburban pay premium raises $\Delta\log Q$ from 0.003 to 0.015, so the SEDA growth gap, which the first fit matched, is now 0.0119 SD/grade, 2.9 SEs above the data. The move rate by tercile now matches the data in level but still rises with parental income, while the data are flat.
- **Identification.** The standardized 8 × 8 Jacobian has full rank, singular values 6.94–210 and condition number 30.3. Pay scale and $\beta$ overlap through teaching shares, but the full set of moments separates them. The Jacobian and normalization are reported in [spatial_calibration.md](spatial_calibration.md).
- **External-block consistency.** Occupation shares move by at most 0.15 pp (men) and 0.75 pp (women); the wage slope is 0.214 versus 0.213 and mother–child correlation 0.5760 versus 0.5766. Pooled non-teacher 90/10 falls to 3.622 versus 3.737. Absorbing this would raise $\sigma_\epsilon$ roughly 3%, within the tested sensitivity range; the first stage is retained.
- **Mean wages (untargeted).** Essentially unchanged: cell deviations from the market mean correlate 0.09 with the data (SD 0.17 vs. 0.26), the market gender gap is 0.30 against 0.10, teachers earn 0.72 (men) and 0.63 (women) log points below the market mean against 0.04 and 0.07 above it, and $s_T/s_O=0.76$.
- **Solution checks.** Stationarity deviation $1.5\times10^{-5}$; the population share choosing a location at the consumption floor is $1\times10^{-23}$. Re-solving at the estimate on finer grids keeps every target matched: the criterion is $4.6\times10^{-6}$ at Nz = 9, $9.3\times10^{-6}$ at Nz = 15 and $2.4\times10^{-6}$ with nϵT = nXO = 96 (largest residual about 0.003 SE), and the validation moments move by at most 0.0001 (`runs/exact-grid-*`).

### Calibration experiments (2026-10)

Six experiments compared with the unchanged baseline on a pre-registered validation table: district earnings, staffing, poverty and SEDA gaps; NLSY parental-income composition and flows; teacher pay and schooling; gender wages and first-stage consistency. Own-AFQT location gradients and teacher ability were added during the experiments. Detailed cases and tables: `slurm/cases_*.tsv`, `runs/experiments-2026-10/` (untracked); final-code evaluations: `runs/exp-final/`, generated by `slurm/experiments_summary.py`. Late housing runs are separate (below; fixed-charge cases revalidated 2026-10-05).

**Prototype switches.** All are off in the baseline: recurring suburban residence cost $r_p\bar C$; suburban taste $\zeta\log z/s_z$; teacher schooling floor $s_T\ge16.04/25$; teaching utility cost $c_T+c_{Tz}\log z/s_z$. The floor uses observed **mean** schooling, not an independently measured minimum. The teaching cost changes occupational utility, not a training technology. These switches diagnose missing mechanisms; they are not proposed baseline assumptions ([model note §12](../../../notes/spatial_teachers_model_eps0.md)).

Experimental teacher-premium targets use occupation-count weights (+0.056/+0.086 men/women; gender gap 0.105), unlike the wage-sample-weighted checks above (+0.042/+0.070; 0.103). Other added targets are schooling ratio 1.176, city-share gradients by parental income −0.295 (SE 0.037) and own AFQT −0.222 (0.036), and teachers' latent AFQT relative to employed non-teachers +0.614 (0.062). All are inactive unless explicitly targeted. `--district SCHEME` changes district targets; `--weight corr` uses their bootstrap correlations. Data: `data/spatial/estimates/settings/`.

**Headline comparison (October 2–3 experiments).** This is the comparison discussed in the calibration review, with equilibrium local tax rates added. Earnings and staffing gaps are suburb minus city in logs; scores are in SDs and learning rates in SD/grade. The two taste specifications target the own-AFQT location gradient, not these four validation outcomes.

| Specification | Earnings gap | Teachers/pupil gap | Score gap, SD | Learning-rate gap, SD/grade | City tax $t_1$, % | Suburb tax $t_2$, % |
|---|---|---|---|---|---|---|
| Data | 0.112 | −0.018 | 0.289 | 0.0031 | — | — |
| Original baseline | −0.010 | 0.023 | 0.060 | 0.0119 | 1.635 | 1.712 |
| Schooling floor + teaching costs, including ability gradient | −0.011 | 0.022 | 0.051 | 0.0102 | 3.267 | 3.419 |
| Ability-dependent suburb taste, fitted to AFQT sorting | 0.072 | −0.147 | 0.115 | −0.0533 | 1.942 | 1.562 |
| Same taste + teacher extensions | 0.173 | 0.112 | 0.765 | 0.0732 | 3.461 | 3.295 |

The teacher extensions alone leave household sorting essentially absent. Adding the suburb taste improves income sorting but displaces teachers from suburbs under the baseline's negative ability selection; with the teacher extensions, it instead concentrates teachers there and overshoots both educational gaps. These are diagnostic results, not an adopted specification.

**Local taxes.** The reported rates are endogenous model outcomes, $t_l=\mathcal W_l/\mathcal I_l$: local teacher payroll divided by the model's taxable income, including home-production output. They fund teacher wages only and are not observed total local tax rates. The baseline is **1.635% in the city and 1.712% in the suburb**; the teacher extensions roughly double both rates. Values are $100\times t$ from each saved `moments.toml`, ordered city then suburb. No empirical tax-rate target is supplied here.

**Additional composition and housing diagnostics.** Poverty and city-share gradients are fractions; the score composition term is the gap with the $Q$ term removed, in the same cohort SD units as the total score.

| Specification | Poverty gap | Score gap without $Q$ | City-share gradient by parental income |
|---|---|---|---|
| Data | −0.089 | Not separately observed | −0.295 |
| Original baseline | +0.002 | −0.005 | +0.011 |
| Teacher floor and costs | +0.003 | −0.005 | +0.012 |
| … + AFQT-targeted taste | −0.097 | 0.384 | −0.165 |
| … at baseline $Q$ (partial equilibrium) | −0.058 | 0.388 | −0.110 |
| AFQT-targeted taste, baseline teacher block | −0.025 | 0.413 | −0.067 |
| Schooling floor + $c_T$ + ACS residence cost | −0.019 | 0.009 | −0.025 |

At baseline $Q$, the bundled taste case gives earnings 0.115, scores 0.453 and growth 0.0119; taxes remain fixed at that case's equilibrium rates (3.461%, 3.295%), not the baseline's tax rates. With the schooling floor, $c_T$ and ACS residence cost, earnings are 0.024, FTE/pupil 0.029, scores 0.086 and growth 0.0140; city/suburb taxes are 3.422%/3.470%.

Sources: `runs/exp-final/{p0_baseline,pr3_teacher_sel,pc_taste_afqt,pr3_teacher_sel_taste_afqt,pc_teacher_housing_acs}/moments.toml` and `runs/exp-counterfactual/teacher_sel_taste_afqt`. Teacher costs target premia and ability selection. Fixed-$Q$ re-solves hold taxes and teacher stocks fixed; they do not satisfy the full teacher-allocation equilibrium.

1. **Revenue partition.** On 621 CZs, salary and FTE gaps are both positive (0.018 and 0.035), unlike the locale split. Holding $r_m$ at its locale value (NLSY mobility does not measure revenue-group transitions), the model fits seven targets exactly. Adding FTE gives criterion 2.0, or 2.3 with correlation weights. Sorting still fails: earnings 0.011 versus 0.129, poverty −0.007 versus −0.036, scores 0.059 versus 0.199. Signs persist on common CZs. State aid offsets 0.47 of each local dollar, so local revenue is not directly a model object.
2. **Target tradeoffs.** Restoring FTE gives criterion 5.1. Adding earnings leaves the gap at or below 0.013 under either weighting. Adding teacher pay and schooling pushes $\gamma$ to 1 and the teacher 90/10 about 6 SE above target. Combining these targets gives weighted criterion 107.
3. **Focused sensitivity.** Refits of class-size curvature, exposure, location-taste dispersion and WTP at the exact baseline change no validation conclusion (`runs/panel-exact/`). Most fit exactly; $\sigma=0.15$ and $K=5$ hit the moving-cost bound. WTP mainly changes $\lambda$; taste dispersion changes moving-cost units and other parameters without generating sorting.
4. **Teacher education/pay.** The schooling floor raises $s_T/s_O$ from 0.76 to 1.17. Adding $c_T$ brings pay near the market; adding $c_{Tz}<0$ changes teacher ability from −0.29 to +0.61 SD. The gender wage gap worsens to 0.39 versus 0.10. NLSY teachers are +0.61 SD against employed non-teachers but −0.24 within college graduates (−0.33 for all teachers against non-teacher graduates). ASVAB at ages 15–23 may already reflect schooling. See [teacher-selection evidence](../../../data/spatial/estimates/settings/nlsy_teachers.md).
5. **Residential sorting.** ACS housing spending implies $r_p\approx0.028$ even if the entire suburb–city difference is a fixed charge. It produces little sorting and displaces the baseline's underpaid teachers. Ability-targeted taste produces strong $z$ sorting, but sends teachers away from suburbs under negative selection and into suburbs under positive selection. Targeting FTE in the latter case still leaves growth +9.6 SE above the data. With $Q$ equalized, the bundled AFQT case gives earnings 0.110, scores 0.400 and parental-income city-share gradient −0.106: commuting alone would not solve the remaining composition mismatch.
6. **Validation.** Half-A refits repeat the failures on half B: the bundled taste case overshoots scores by 14 SE and growth by 16 SE. Although half-B results were sealed, model selection used full-sample moments; this is a stability check, not independent out-of-sample validation. CRDC shows fewer novice teachers in suburbs (10.2% versus 13.1%) and more certified teachers (+2.1 pp). Converting the novice gap into about +0.01 log teacher human capital requires assumed experience effects and the CFR mapping; it is suggestive, not direct validation of teacher productivity.

#### Explicit housing / cost-of-living experiment (2026-10-02–03)

**Question and specification.** Could a recurring suburban cost generate persistent income sorting and improve the earnings, poverty and SEDA validation gaps? The prototype adds $p=(0,r_p\bar C)$ to the household budget: $C_{l'}=y_{l'}-(1+\tau^e)e-m_{l,l'}-p_{l'}$. Every suburban resident pays it each model period, including stayers. A fixed charge takes a larger share of poorer households' resources, so it can sort the resident stock as well as movers. This is a reduced-form housing/access or COL differential; it does not introduce housing quantities, endogenous rents, supply or a general COL price index. The baseline keeps $r_p=0$. See [model note §12](../../../notes/spatial_teachers_model_eps0.md).

**Which location costs more in the data? Suburbs.** The [housing-cost estimates](../../../data/spatial/estimates/settings/housing_costs.md) use ACS 2014–18 housing tables merged to the 2018 school-district panel, with city locales 11–13 and suburb locales 21–23. Across 188 commuting zones, pupil-weighted means of within-CZ pooled median monthly housing costs are **$1,224 in cities and $1,365 in suburbs**, in 2018 dollars: **+$141/month (SE $22)** and a mean within-CZ log gap of **0.114 (SE 0.017)**. Renter-weighted district median gross rent is also higher in suburbs, by $69/month (SE $15). These are housing-spending comparisons within this sample, not quality-adjusted prices or a universal urban/suburban ranking.

Converting the annual housing-spending gap using CE-2018 mean expenditures per consumer unit ($61,224.13) gives $r_p=0.0276$ (SE 0.0044), rounded to **0.028** in the experiment. Treating the entire spending gap as a fixed charge is an upper-bound interpretation: it includes housing quantity/quality, tenure and mortgage principal, and is not adjusted for household size. It also assumes housing costs and consumption scale together when mapping household data to the per-adult model. Locale measures imply $r_p$ around 0.013–0.033; **0.058** is a broad sensitivity value from the range across measures and partitions, including SES, rather than the upper confidence limit for the city–suburb estimate.

**Saved results.** The first two housing cases fix the charge and refit the original eight targets. The third frees $r_p$ in [0, 0.15] and adds earnings to the criterion, so earnings is targeted there. Gaps are suburb minus city; earnings and staffing are logs, scores are SDs and growth is SD/grade.

| Case | $r_p$ | Earnings gap | FTE/pupil gap | Score gap | Learning-rate gap | Move rate |
|---|---|---|---|---|---|---|
| Data | — | 0.112 | −0.018 | 0.289 | 0.0031 | 0.216 |
| Original baseline | 0 | −0.010 | 0.023 | 0.060 | 0.0119 | 0.216 |
| ACS charge, baseline teacher block (`housing_acs`) | 0.028 | −0.004 | −0.063 | −0.120 | −0.0230 | 0.216 |
| High charge, baseline teacher block (`housing_acs_hi`, revalidated) | 0.058 | 0.005 | −0.152 | −0.306 | −0.0597 | 0.216 |
| Charge fitted with earnings (`housing_income`) | 0.0798 | 0.077 | −0.169 | −0.370 | −0.0772 | 0.314 |
| ACS charge + teacher schooling floor + $c_T$ | 0.028 | 0.024 | 0.029 | 0.086 | 0.0140 | 0.219 |

Sources: [case definitions](slurm/cases_housing.tsv), saved `theta.toml`, `moments.toml` and `report.txt` in [`runs/exp-housing/`](runs/exp-housing/), and [`runs/exp-final/pc_teacher_housing_acs/`](runs/exp-final/pc_teacher_housing_acs/). The last row also targets teacher wage premia and is the combined case reported above. `housing_nlsy` (charge fitted to parental-income location sorting) and `housing_taste` (charge plus ability-dependent taste fitted to earnings and parental-income sorting) have report headers but no saved final parameters or moments; no completed estimates are reported for them.

**Interpretation and numerical status.** At the ACS charge, income sorting remains nearly absent and suburban staffing and school-quality gaps turn negative under the baseline's teacher selection. Teacher extensions restore positive school gaps but leave earnings sorting weak. The fixed-charge cases match their eight calibration targets; the earnings-targeted fit has criterion 10.82 and sacrifices mobility while still missing earnings. **Revalidation (2026-10-05).** The runs were refit on code with the floored-node spline cap, the relaxed solver and `--floor-tol 1e-5` ([`runs/exp-housing-reval/`](runs/exp-housing-reval/)).
- **ACS charge:** results unchanged; floor mass $3\times10^{-36}$.
- **High charge:** still an exact fit, now with nonnegative floor mass ($1\times10^{-28}$; previously −0.000015). Its row above is the revalidated one; the original was 0.007, −0.134, −0.273, −0.053.
- **Earnings-targeted row:** superseded. Spurious teachers at floored nodes had moved its salary gap from 0.014 to 0.076. Its refit stopped at criterion 7.0 after 26 iterations, still improving slowly (originally 10.8), when its node lost the ResearchDrive mount.
- **Free-charge fits:** none is exact.
  - Fitted to parental-income sorting: criterion 34, $r_p$ 0.119, gradient −0.091.
  - Charge plus taste: criterion 21, with ζ at its upper bound of 1.0; earnings 0.138.
  - Both empty suburban schools: FTE −0.72/−0.78, score −1.20/−0.76, growth −0.29.
- **Reading:** under the bundled baseline teacher block, a uniform charge cannot sort households without draining suburban schools. This matches the 2026-10-05 results below. These saved prototypes do not establish a validation improvement or justify replacing the baseline. Failure of a uniform residence charge does not reject an endogenous housing market.

**Assessment.** The experiments identify structural restrictions rather than a replacement calibration. Household sorting must coexist with modest school-effectiveness differences, and teacher schooling/selection need an economic explanation. Within the current score mapping, observed growth contributes approximately $5.5\times0.0031=0.017$ SD to the 0.289 score gap; the remaining roughly 0.272 SD is composition or other accumulated differences. This is model accounting, not a causal decomposition of SEDA. The next experiments are in §5.

### Model experiments (2026-10-05)

The 10-03 mechanisms were tested under a pre-registration written before any result (`runs/experiments-2026-10b/PREREG.md`, amendments 1–3; specs, reviews, working `LOG.md` and `scores_final.md` in the same untracked directory; runs in `runs/exp-{m2,m3,c3}/`).
- **Score:** $S=\sum z^2+2k$ over a fixed validation set: V1 earnings, V2 FTE, V3 poverty, V4 score, V5 growth, V6/V6z city-share gradients, V9 teacher premia, V11 gender gap, V12 wage slope and mother–child ρ, and payroll per pupil. $k$ counts parameters beyond the base eight. Lowest $S$ is selected.
- **Sealed set:** NLSY97 and half-B district moments, opened only after selection.
- **Numerics:** every fit is exact (criterion ≤ $2\times10^{-19}$; floor mass ≤ $10^{-28}$) on main-tree code checked bitwise against the prototypes and the baseline.

| Candidate | Added mechanism (free parameters ↔ added targets) | $k$ | $S$, data SE | $S$, + model SE |
|---|---|---|---|---|
| C0 | Baseline | 0 | 797 | 352 |
| C1 | M2: all workers choose a workplace (position-weighted logit; commuting cost $\tau_k=0.03$ of the wage); $\rho_w$ ↔ crossing share 0.107 | 1 | 801 | 300 |
| C2 | M3: BA choice; licensure with a non-BA teaching track and MA margin; $s_z,\sigma_\epsilon$ internal. 9 parameters ↔ BA shares, premium and AFQT gap, BA wage slope, schooling, pooled 90/10, BA − non-BA teacher pay | 9 | 384 | 274 |
| C3 | M2 + M3 | 10 | 437 | 215 |
| **C4′** | C3 + suburban access charge $r_p=0.04\bar C$ (fixed; replaces M1, below) | 10 | **284** | **110** |

Unscored variants:
- $r_p=0.075$: $S$ = 266.
- Staffing closure ($\delta_\kappa$ ↔ FTE, salary as validation): C3 535, C4′ 376.

Model SEs come from the delta method (`spatial_uncertainty.jl`). They leave the order unchanged. For C4′ they come mostly from the taste dispersion $\sigma$, whose linearization is unreliable, so the data-SE score is primary.

C4′ parameters:
- Base eight: $\log\tilde\kappa$ −0.778, $\delta_\kappa$ −0.002, $\Delta B$ 0.109, $\beta$ 0.239, $\lambda$ 0.095, $r_m$ 0.245, $\omega_f$ 0.932, $\gamma$ 0.610.
- Workplace and college blocks: $\rho_w$ 0.068, $\rho_c$ −0.110, $c_C$ 0.084, $c_{Cz}$ −0.468, $\rho_{T0}$ 0.410, $\vartheta_s$ 1.143, $s_z$ 0.151, $\sigma_\epsilon$ 0.803, $d\log\phi$ −0.230, $c_{T0}$ 0.510.

| Validation († targeted) | Earnings | FTE/pupil | Poverty | Score | Growth | V6 | V6z | Premia m/f | Gender gap | Teachers' z vs all |
|---|---|---|---|---|---|---|---|---|---|---|
| Data | 0.112 | −0.018 | −0.089 | 0.289 | 0.0031 | −0.295 | −0.222 | 0.06/0.09 | 0.105 | +0.61 |
| C0 | −0.010 | 0.023 | 0.002 | 0.061 | 0.0119 | 0.011 | 0.003 | −0.72/−0.63 | 0.302 | −0.29 |
| C1 | −0.001 | 0.043 | −0.004 | 0.117 | 0.0220 | 0.002 | 0.002 | −0.72/−0.63 | 0.302 | −0.29 |
| C2 | −0.001 | 0.036 | −0.002 | 0.092 | 0.0173 | 0.003 | 0.002 | −0.10/0.01 | 0.357 | +0.55 |
| C3 | 0.013 | 0.065 | −0.011 | 0.179 | 0.0330 | −0.010 | 0.001 | −0.10/0.01 | 0.357 | +0.55 |
| **C4′** | 0.034 | 0.020 | −0.023 | 0.082 | 0.0127 | −0.034 | −0.007 | −0.10/0.00 | 0.356 | +0.55 |
| C4′ $z$: data SE / + model SE | −4.6/−2.1 | +2.9/+0.5 | +7.3/+3.2 | −6.3/−1.4 | +3.2/+0.3 | +7.0/+4.9 | +6.0/+5.7 | m −3.2/−1.5, f −1.6/−0.8 | +5.0/+3.5 | −1.0/−0.9 |
| ζ diagnostic (below) | 0.119 | −0.018† | −0.058 | 0.436 | 0.0057 | −0.112 | −0.222† | −0.10/0.01 | 0.358 | +0.55 |

**Model levels.** Some levels are the same in every case: population 0.407/0.593 city/suburb, pinned by the pupil-gap target; teaching shares 1.9%/6.0% (men/women); home production 22.6%/34.1%; and $P(1\to2)$, $P(2\to1)$ = 0.265, 0.182. Taxes are $t_l$ as defined above.

| | C0 | C1 | C2 | C3 | C4′ | Data |
|---|---|---|---|---|---|---|
| Local tax $t$, city/suburb (%) | 1.63/1.71 | 1.62/1.72 | 2.94/3.10 | 2.91/3.12 | 3.04/3.02 | — |
| Teacher payroll / market income (%) | 2.53 | 2.53 | 4.54 | 4.54 | 4.54 | ≈1.6 (not like-for-like) |
| $\Delta\log Q$ | 0.015 | 0.028 | 0.022 | 0.041 | 0.016 | — |
| Pupils per teacher, city/suburb | 25.6/25.0 | 25.9/24.8 | 25.8/24.9 | 26.2/24.6 | 25.5/25.0 | — |
| Mean schooling $s$, teachers/others | 0.42/0.55 | 0.42/0.55 | 0.62/0.55 | 0.62/0.55 | 0.62/0.55 | ratio 1.18 |
| BA share | — | — | 0.345 | 0.345 | 0.345 | 0.339 |
| Children's mean log $z$, suburb − city (SD) | −0.005 | −0.004 | −0.003 | −0.002 | +0.013 | — |
| WTP income elasticity | 1.09 | 1.09 | 1.10 | 1.10 | 1.12 | 0.38 (SE 0.17) |
| Teachers' crossing share | — | 0.118 | — | 0.114 | 0.113 | 0.084 |
| Teachers per employed resident, gap | — | −0.007 | — | −0.018 | 0.001 | +0.298 |
| Staffing-ratio elasticity to suburban pay (PE) | — | 31 | — | 32 | 31 | No matched counterpart yet |

**Elasticity qualification (2026-10-15).** The last row reports $\partial\log(n_2/n_1)/\partial\log\kappa_2$ at fixed aggregates, with $\kappa_1$ unchanged. The previously displayed literature range 1.4–3.7 concerns own-employer teacher supply. It is a comparator requiring a mapping, not a direct data counterpart; see the [proposed audit](#recommended-sequence-2026-10-15).

- **M1, income-driven demand: not built.** The baseline's WTP income elasticity is already 1.09 against BFM's 0.38. The suburb's school quality is worth 4.6×10⁻⁴ utils to the average parent, against ≈0.4 needed. The endogenous-rent block was gated on a data price wedge ≥ 0.08 $\bar C$. The ACS/CE wedge is ≈0.036 $\bar C$ (missing cheap suburban rentals 0.023 + car 0.013), with an upper value of 0.075, so the gate failed. C4′ instead prices access at 0.04.
- **M2, workplace split.** The crossing share pins a small workplace dispersion: most teachers work where they live, yet staffing becomes very pay-elastic. Alone (C1), it worsens FTE and growth. The model's teachers are negatively selected and are not suburbanites (net commuting −0.05 against +0.075). The split decouples staffing from residential sorting only with positive teacher selection.
- **M3, education and occupation.** M3 produces positive teacher selection overall and negative selection within college graduates (−0.29 against −0.33). It raises teacher pay close to market without teacher-specific utility costs. The licensure barrier on the non-BA track ($c_{T0}$) carries the teachers' BA share: without it the share is 0.53, against 0.87 in the data. Misses: BA teachers within BA −0.12 (data −0.24), the MA premium and the gender gap. Sorting is unchanged. First-stage recheck (`--consistency`): $s_z$ falls to 0.151, as the BA-conditional wage slope predicts (0.148). The unconditional wage slope is 0.199 against 0.213 (−1.7 SE), and mother–child ρ is 0.576 against 0.577.
- **Combined (C4′).** The charge adds a little sorting: earnings 0.034, V6 −0.034, children's ability gap +0.013 SD. Unlike in the bundled model, it does not drain suburban schools (FTE +0.020, growth 0.013). Sorting and the gender gap remain the main misses.
- **Diagnostic: M2 + M3 absorb a sorting force.** The ζ row replaces the charge with the AFQT-graded suburb taste (ζ = 0.21, fitted to V6z) under the staffing closure ($S$ = 125). Earnings, score and growth come close to the data, and teachers per employed resident is +0.36 (data +0.30). Net commuting overshoots: +0.39. With salary targeted instead, FTE is −0.16. ζ is a stand-in, not a mechanism: it shows what a sorting mechanism must deliver, and that the extended teacher block can absorb it.
- **Sealed set (opened after selection).** The sorting miss replicates in a new cohort, and M3's selection pattern holds out of sample.
  - NLSY97 city share by own AFQT, top − bottom tercile: −0.157 (design SE 0.026). NLSY79 −0.222; C4′ −0.007.
  - NLSY97 teachers' latent AFQT: vs employed non-teachers +0.44 (0.06; model +0.55); vs BA non-teachers −0.28 (0.06; model −0.29); BA teachers within BA −0.20 (0.06; model −0.12).
  - Half-B city-share gradient of children by family income: −0.189 (0.022), against −0.185 in half A.
  - The teachers' wage–AFQT slope is not in the sealed file.
  - Disclosure: the read also displayed NLSY97 flows, BA shares and adult suburb shares, which the pre-registration had set aside. They are not used.

### Model experiments (2026-10-07)

Items 1–2 of the 2026-10-07 list (§5), on main-tree code: working log `runs/experiments-2026-10c/LOG.md` (untracked), fits in `runs/exp-roy/` and `runs/exp-peer/`. These are diagnostics, not a pre-registered comparison, because the NLSY97 sealed set is spent. Every case is refitted exactly (criterion ≤ $10^{-15}$) and scored with `score_candidates.py` on the 10-05 set, with data SEs.

**Regular Roy model, $\gamma=1$.** Teaching pays $\kappa_l h$, linear like every other occupation, so $z$ and $Q$ cancel from the teaching threshold and $s_T=s_O$ ([note §12](../../../notes/spatial_teachers_model_eps0.md)). γ is fixed at 1 (`--fix gamma=1`). Its target, teachers' 90/10, becomes validation (`--deactivate teacher_p9010`), which leaves 7 × 7 in the baseline and 17 × 17 in C4′; $k$ falls by one.

| | $S$ | Teachers' 90/10 | Premia m/f | Teachers' $z$ vs all | $s_T/s_O$ | Gender gap | FTE/pupil | Growth | Score | Earnings |
|---|---|---|---|---|---|---|---|---|---|---|
| Data | | 2.77 | 0.06/0.09 | +0.61 | 1.18 | 0.105 | −0.018 | 0.0031 | 0.289 | 0.112 |
| C0 | 797 | 2.77† | −0.72/−0.63 | −0.29 | 0.76 | 0.302 | 0.023 | 0.0119 | 0.061 | −0.010 |
| C0, γ = 1 | 408 | 3.75 | −0.09/0.01 | 0.00 | 1.00 | 0.482 | 0.047 | 0.0213 | 0.113 | −0.001 |
| C4′ | 284 | 2.77† | −0.10/0.00 | +0.55 | 1.14 | 0.356 | 0.020 | 0.0127 | 0.082 | 0.034 |
| C4′, γ = 1 | 372 | 3.79 | +0.25/+0.37 | +0.66 | 1.17 | 0.533 | 0.002 | 0.0038 | 0.033 | 0.033 |

- **Baseline.** ϑ: $\log\tilde\kappa$ −1.061, $\delta_\kappa$ 0.010, $\Delta B$ 0.068, β 0.273, λ 0.118, $r_m$ 0.198, $\omega_f$ 0.828.
  - Teacher pay moves to market pay (premia −0.09/0.01). A lower $\omega_f$ holds the female teaching share, which widens the market gender gap to 0.48.
  - Selection on ability disappears: teachers' $z$ is the population's (data +0.61 SD), and their 90/10 rises to the market's, 3.75 (+7 at the 5% scale).
  - The salary gap buys more staffing: FTE +0.047 and growth 0.021 ($z$ = +5.0 and +6.1). Sorting is unchanged.
  - Levels: $t$ 3.0/3.2%, teacher payroll 4.6% of market income.
- **Identification.** Condition number 47.8 (C0 30.3). With $Q$ cancelling, β no longer moves the teaching shares: the $R^2$ of the $\tilde\kappa$ and β columns on the others falls from 0.90 to 0.43 and 0.70.
- **Score.** 408 against 797. All of the gain is the teacher premia (V9 446 → 12); FTE, growth and the gender gap worsen. Counting the released 90/10 gives ≈ 458.
- **C4′.** ϑ stays near C4′ except $\omega_f$ 0.79 (0.93), $\rho_{T0}$ 0.72 (0.41), $\vartheta_s$ 0.78 (1.14) and $c_{T0}$ 0.28 (0.51).
  - The licensure track still selects teachers positively (+0.66 SD; data +0.61), and the schooling ratio fits (1.17). Within the BA, selection disappears: BA teachers −0.01 SD against −0.24.
  - Teachers are therefore overpaid (premia +0.25/+0.37), and their 90/10 is the market's (3.79; BA teachers 3.65 against 2.42).
  - The staffing gradient fits: FTE +0.002, growth 0.0038 ($z$ = +1.5, +0.2), payroll per pupil +0.011, with $\Delta\log Q$ 0.005.
  - The gender gap widens to 0.53 and the score gap falls to 0.033. Taxes rise to 4.4%, and salaries reach 6.4% of market income.
  - Score 372 against 284, or ≈ 426 counting the 90/10.
- **Reading.** In the bundled baseline, γ = 1 is a cheap fix for the teacher premia, which dominate its score, but it gets selection and dispersion wrong. M3 already fixes the premia at γ ≈ 0.61 while keeping the negative selection within the BA and targeting the 90/10. Setting γ = 1 there improves the staffing gradient and SEDA growth, but worsens pay, dispersion and the gender gap. C4′ stays preferred: γ < 1 is what draws the weaker BA graduates into teaching.

**Local peer effects, $\psi_a$.** School quality carries the students' mean ability: $Q_l=(2\tilde H_{T,l}/M_l)^\sigma\,\bar a_l^{\psi_a}$, where $\ln\bar a_l$ is the mean log $z$ of the young schooled in $l$ minus its ergodic mean (`Params.ψa`; ϑ entry `ψ_a`, alias `psi_a`, box [0, 1]; note §12).
- ā is a fourth aggregate of the outer fixed point. $\psi_a=0$ reproduces the earlier fits bitwise (C0 and C4′ checked), and new moments split $\Delta\log Q$ into teacher and peer parts.
- $\psi_a$ has no target of its own, so it is fixed at 0.25 or 1 and each specification is refitted to its own targets.

| | $\psi_a$ | Children's $z$ gap (SD) | $\Delta\log Q$: teachers, peers | FTE/pupil | Growth | Score | Earnings | V6 | $S$ |
|---|---|---|---|---|---|---|---|---|---|
| Data | | | | −0.018 | 0.0031 | 0.289 | 0.112 | −0.295 | |
| C0 | 0 | −0.005 | 0.015, 0 | 0.023 | 0.0119 | 0.061 | −0.010 | 0.011 | 797 |
| | 0.25 | −0.005 | 0.015, −0.0002 | 0.024 | 0.0119 | 0.060 | −0.010 | 0.011 | 798 |
| | 1 | −0.005 | 0.016, −0.0010 | 0.025 | 0.0117 | 0.059 | −0.010 | 0.011 | 799 |
| C4′ | 0 | +0.013 | 0.016, 0 | 0.020 | 0.0127 | 0.082 | 0.034 | −0.034 | 284 |
| | 0.25 | +0.013 | 0.016, 0.0005 | 0.019 | 0.0128 | 0.083 | 0.034 | −0.034 | 283 |
| ζ benchmark | 0 | +0.41 | 0.007, 0 | −0.018† | 0.0057 | 0.436 | 0.119 | −0.112 | 125 |
| | 0.25 | +0.41 | 0.010, 0.015 | −0.018† | 0.0199 | 0.509 | 0.129 | −0.122 | 176 |
| | 1, not refitted | +0.41 | 0.001, 0.062 | −0.045 | 0.0503 | 0.662 | 0.158 | −0.149 | 471 |

- **Peers do not create sorting.** At the baseline and C4′ the children's ability gap is near zero. The peer part of $\Delta\log Q$ is therefore at most 0.002, and the refits move only $\delta_\kappa$. Starting the fixed point from sorted states (log ā = ±0.1, half an SD of log $z$) returns to the same equilibrium. With λ pinned by the BFM WTP, parents value school quality far too little (M1 above) for peer ability to start a sorting spiral.
- **With sorting imposed, peers add quality, not sorting.** In the ζ benchmark (children's ability gap +0.41 SD, held by ζ ↔ V6z), $\psi_a=0.25$ adds 0.015 to $\Delta\log Q$, and refitting $\delta_\kappa$ to FTE raises the teacher part as well (0.007 → 0.010). SEDA growth rises from 0.006 to 0.020 (data 0.003, SE 0.003). The score and earnings gaps, already above the data here, rise further (0.509, 0.129). ζ barely moves (0.215 → 0.214): peers do not substitute for the sorting force.
- **Bound.** Growth is linear in $\psi_a$. If sorting is as large as the data require, the SEDA growth gap bounds $\psi_a\lesssim0.06$ at +2 SE (0.08 without the refit), with no room at the point estimate. In the score mapping, that is at most 0.04 score SDs per SD of the peers' mean ability. The bound inherits the mapping, which counts school quality, peers included, as growth over 5.5 grades; peer effects that act before grade 3 would escape it.

### Model experiments (2026-10-15)

This is the [recommended sequence](#recommended-sequence-2026-10-15), run as a pre-registered round. The pre-registration (amendments 1–5), case lists, scores and working log are in `runs/experiments-2026-10d/` (untracked); fits are in `runs/exp-1015/`. Code before the round is frozen in `runs/_code/pre-20261015/`. Every new mechanism is a switch that is off by default, checked bitwise against C0 and C4′, and now in the main tree (note §12).

**Scoring.** $S_{1015}$ is the 10-05 set with V6 replaced by its ACS within-CZ counterpart V6a, −0.187 (SE 0.017), scored with data SEs plus $2k$. The acceptance criteria are measured against C4′ ($S_{1015}$ = 315.1):
- A1: lower $S$.
- A2: no sorting moment worse by more than 1 SE, and at least three better by 1 SE.
- A3: FTE, growth and payroll $|z|$ at most C4′'s + 1.
- A4: teachers' $z$ at least 0.40 SD and a within-BA gap of at most 0.
- A5: an exact fit and grid checks.
- A6: small-district supply elasticity at most 10 and own-location elasticity at most 3.7.
- A8: for utility moving costs, the own-AFQT move gradient within 2 SE of NLSY79's +7.8 pp.

**Supply audit** (`spatial_supply_audit.jl`; [table](runs/exp-1015/audit/README.md)). The ratio elasticity 31 at C4′ is the household-PE response of $n_2/n_1$, and 91% of it is workplace switching inside $\sigma_w\approx0.017$. In full equilibrium, $\varepsilon_{22}$ is 5.9 and the ratio 12.3. The small-district elasticity $\mu y/(\sigma_wC)$ is ≈68; reaching 3.7 would need $\sigma_w>\sigma_\nu$, which a RUM does not allow.

The repair **S1** fixes $\rho_w=1$ and fits a utility commuting friction $\kappa_c$ (0.49) to the crossing share. Under S1 the full-GE $\varepsilon_{22}$ is 1.6 and the small-district elasticity 4.7, within the Falch–Ransom–Sims range. The literature note (`audit/literature.md`) supports a flat residence–workplace nest. The direction of non-teachers' commuting stays reversed in every fit (city→suburb 0.13 and suburb→city 0.09, against ACS 0.09 and 0.12).

**Candidates** (all exact unless stated; $z$ in parentheses):

| | $S_{1015}$ | sorting / allocation | Earnings, score, V6a, V6z | FTE, growth, salary | own-AFQT move gap | fails |
|---|---|---|---|---|---|---|
| Data | | | 0.112, 0.289, −0.187, −0.222 | −0.018, 0.0031, 0.009 | +7.8 pp | |
| C4′ | 315.1 | 232 / 23 | 0.034, 0.082, −0.034, −0.007 | 0.020 (2.9), 0.0127 (3.2), † | +9.2 | — |
| S1 | 393.7 | 143 / 190 | 0.056, 0.222, −0.055, −0.009 | 0.064 (6.2), 0.038 (11.5), † | +9.1 | A1, A3 |
| U (C4′, utility moving cost) | 404.0 | 299 / 45 | 0.037, −0.064, −0.036, −0.013 | −0.041, −0.016, † | −0.2 | A1–A3, A6, A8 |
| S1·U | 279.3 | 134 / 85 | 0.067, 0.170, −0.065, −0.014 | 0.037 (4.2), 0.026 (7.7), † | −0.2 | A3, A8 |
| S1·E (ζ_BA ↔ V6a) | 1708.5 | 234 / 1412 | 0.194, 0.741, −0.187†, −0.123 | 0.169, 0.103, † | +5.3 | A1–A3 |
| S1·U·E | 1523.1 | 193 / 1267 | 0.194, 0.694, −0.187†, −0.111 | 0.166, 0.097, † | −1.5 | A1–A3, A8 |
| S1·SC (δκ ↔ FTE) | 448.8* | 342 / 46* | 0.014, −0.028, −0.016, −0.005 | †, −0.0068, −0.040 (−5.3) | +9.2 | A1–A3 |
| S1·U·SC | 310.7* | 235 / 16* | 0.045, 0.013, −0.045, −0.013 | †, −0.0020, −0.021 (−3.2) | −0.2 | A2, A8 |
| S1·R (b_T ↔ FTE; round b) | 407.1 | 338 / 12 | 0.017, −0.030, −0.017, −0.005 | †, −0.0072 (−3.4), † | +9.3 | A1, A2 |
| **S1·R·E** (round b) | **105.1** | 40 / 3 | 0.201, 0.398, −0.187†, −0.225 | †, −0.0021 (−1.7), † | +3.0 | A5 (Nz = 9) |

\* With FTE targeted, the salary gap takes V2's place in the score (amendment 4).

Five fits did not reach an exact fit and were stopped. Their best iterates were re-evaluated (`exp-1015/reeval/*_last`):
- **Θ (θ_ν ↔ V6a), E and U·E on C4′.** Their criteria stalled at 9.1, 5.0 and 4.1, with θ_ν 0.044 and ζ_BA ≈ 1. With ρw free near C4′'s, they also fail A6 by construction.
- **S1·SC·E.** It stalled at criterion 0.25, holding FTE against the amenity with a suburban pay cut: salary −0.198 (z −22.5), S* 752.
- **S1·Θ.** It stalled at 12.1. It reaches V6a only through λ 0.45 (S1: 0.095) and fails A1, A3 and A6 (own elasticity 4.15).

**Reading.**
- *Teacher supply.* S1 makes it realistic, and then staffing exposes the missing demand side. With δκ on salary, FTE follows teachers' residence (+0.064). With δκ on FTE (SC), the suburban pay cut (salary −0.040) also removes the suburb's school-quality advantage, and the score gap (−0.03) with it. The data's score gap, 0.289 SD at FTE −0.018 and growth 0.003, must therefore come from household composition.
- *U.* Moving costs in utils release poor households and give the best score (S1·U, 279). But they remove the rise of mobility with own ability, which the criterion written for this risk rejects.
- *E.* The BA amenity matches V6a but carries BA teachers (87% BA) to the suburb, so it fails allocation by a wide margin.
- *Outcome of round a.* No candidate meets the criteria, and C4′ stays the development candidate. In the follow-up round, S1·R·E meets every criterion except A5's grid check ([below](#follow-up-round-2026-10-15b-a-non-pecuniary-staffing-closure)).
- *Hold-out.* The sealed B1 (BA share by residence) has not been opened: the repository's standing instruction requires the user's explicit go-ahead (amendment 5). Model counterparts range from +0.007 (C4′) to +0.30 (S1·E) and +0.57 (S1·R·E).

**Persistence and family inputs** (exploratory, `exp-1015/explore/`). Raising $\rho_z$ from 0.58 to 0.90 at fixed ϑ, with $s_z$ held, moves V6a by −0.004 at C4′ and at S1, V6z by −0.008 and −0.012, and the score gap by 0.02–0.03. A family component known to parents works through the same predictable-ability term, so this bounds it from above. The suburb's school quality is worth about 0.002 utils to parents, against a taste scale $\sigma_\nu$ = 0.25. The reason is that λ is pinned by the BFM WTP and $\Delta\log Q$ by SEDA growth. No transmission mechanism that runs through school quality can amplify a force that small, so shared family inputs are not the next candidate.

#### Follow-up round 2026-10-15b: a non-pecuniary staffing closure

The round is informed by round a: its pre-registration cites those results. Files: `runs/experiments-2026-10e/` (PREREG, LOG, `scores_all.md`, `scores_profile.md`, `scores_nz9.md`; untracked); fits in `runs/exp-1015b/`; code `runs/_code/r1015r`.

**Design.** The data fill suburban positions at equal pay (salary +0.009, FTE −0.018), even though teachers live disproportionately in the suburb (+0.298 per employed resident). The closure **R** therefore puts the staffing wedge in non-pecuniary terms. It is teachers' value of a suburban teaching job, $b_T$: the existing `Params.BT`, now the 27th entry of ϑ (note §12), fitted to FTE while $\delta_\kappa$ stays on salary. It is a reduced form of rationed positions or of workplace preferences. Candidates: S1·R and S1·R·E (with ζ_BA ↔ V6a), under round a's criteria.

**S1·R** (exact, b_T −0.059) fails A1 and A2, as pre-registered. Like SC it removes the excess suburban staffing, and with it the school-quality part of the score gap (−0.030). It creates no sorting. Unlike SC, it needs no suburban pay cut.

**S1·R·E** is exact at ζ_BA = 1.3945, b_T −0.345, δκ −0.108 and c_C 1.23.
- *How the fit was found.* The LM crept along a ridge (ζ_BA with the college cost c_C that holds the BA share). So the fit was found through a profile: ζ_BA fixed on a grid, the other 19 targets fitted exactly, and V6a(ζ_BA) read off.
- *V6a saturates.* It is −0.175, −0.180, −0.183, −0.186, −0.187 and −0.188 at ζ_BA 0.90, 1.00, 1.10, 1.25, 1.40 and 1.50. Its slope at the root is −0.0075 per util, so one SE of V6a spans about 2 utils of ζ_BA.
- *Confirming the fit.* A secant step and a re-evaluation of the pre-registered case at the root give criterion 2.3e−15 over its 20 targets.

The other targeted and scored results:
- **Sorting** (z in parentheses): earnings 0.201 (+5.3), poverty −0.099 (−1.1), scores 0.398 (+3.4), V6z −0.225 (−0.1).
- **Allocation:** growth −0.0021 (−1.7). FTE and salary are targeted, and payroll follows from them.
- **Supply:** own-location elasticity 3.47, small-district elasticity 4.67.
- **Teachers:** latent z 0.546; within-BA gap −0.13.
- **Score:** $S_{1015}$ = 105.1, against C4′'s 315.1.

S1·R·E meets A1–A4, A6 and A5's convergence conditions (criterion 2.3e−15, floor mass 1e−94, stationarity 2.1e−5). **It fails A5's grid check.** Re-solved at Nz = 9 with the same ϑ, the score gap moves −0.32 SE and SEDA growth −0.26 SE, against a limit of 0.25; V6z moves +0.24 and V6a +0.13. Its strong ability sorting needs a finer ability grid than five nodes.

Round b's PREREG writes the grid check into A5 "for a selected candidate". Round a ran it only "before adoption". Read strictly, S1·R·E fails A5, so no candidate meets A1–A6 and **C4′ stays the development candidate**. S1·R·E is the best-scoring candidate, and the only criterion it fails is A5's grid check.

**Checks of S1·R·E** (§5 item 5, run for the best candidate):
- *Grid.* It fails at Nz = 9, as above. C4′ passes the same re-solve, with no scored moment moving more than 0.069 SE, so the sensitivity comes from S1·R·E's strong ability sorting: the children's log-ability gap is 0.41 SD, against 0.01 at C4′.
- *Identification.* The standardized Jacobian's condition number is 2,320, against C4′'s 56. The smallest singular value is 0.50, against 7.0 for the next. ζ_BA and c_C are almost collinear with the other columns (R² = 1.000), and the effect the others cannot mimic is only 0.58 and 0.99 box widths; every other parameter is at 7.5 or more.
- *Multistart.* Eight LM polishes at the root, with every other parameter perturbed by ±15%, all return to the same solution: ϑ within 2.6e−10 box widths, scored moments within 1e−8 SE. With the strictly decreasing profile, the fit is unique in these neighbourhoods. Its weakness is a flat direction, not a second solution.
- *Finer grid.* The ζ_BA profile was re-fitted at Nz = 9, with the other 19 targets exact. V6a is −0.18657, −0.18690 and −0.18716 at ζ_BA 1.3945, 1.45 and 1.50, so the Nz = 9 root lies at ζ_BA ≈ 1.47, inside the box but close to its bound.
  - Along that profile the scored moments match the Nz = 5 fit: V4 0.397–0.398, growth −0.0021 to −0.0027, V1 0.200–0.201, V6z −0.220 to −0.222. $S_{1015}$ is 104.7–106.7, and the points pass A1–A4 and A6.
  - So the grid failure reflects how the coarse grid maps ϑ to moments, and it moves the weakly identified ζ_BA by about 0.08 utils. The fitted moments are robust to the grid. That covers both S1·R·E's successes and its validation failures: teachers +1.57 and BA share +0.57 at Nz = 9.
  - Re-fitting at Nz = 9 and checking at Nz = 13 would repair A5 numerically. That step is informed by the failure, so it was not run as part of the pre-registered selection.

**What its validation says.** E and R are reduced forms, and the untargeted moments show them doing implausible work.
- *Teachers.* They are 87% BA, so they live with other BA holders in the suburb and commute back to city schools.
  - Teachers per employed resident +1.54 (data +0.298).
  - Teachers' net commuting +1.57 (+0.075).
  - Teacher crossing share 0.29 (0.084).
  - To hold staffing, b_T must be worth about 30 log points of suburban pay (b_T C/(μy)). The PREREG expected +0.4 and +0.2 to +0.35 for the two residence moments.
- *Sorting overshoots in income and scores.* BA-graded sorting carries more income and ability than the data's parental-income gradient implies.
- *Teacher premia rise to +0.13/+0.24* (data +0.06/+0.09). The pay level offsets the suburban workplace penalty.
- *Directional mobility by own AFQT becomes too sorted.* City→suburb is 0.12 and 0.57 for the bottom and top terciles (data 0.28 and 0.51). The own-AFQT move gap is +3.0 pp against +7.8.
- *ζ_BA ≈ 1.39 utils* is 5.6 σν and weakly identified.
- *The BA-share counterpart is +0.57,* against C4′'s +0.007. The sealed B1 tests exactly this; it has not been opened (amendment 5; A7 waits for the user's instruction).

**Reading.** S1·R·E shows what the model needs: a sorting force that reaches V6a and V6z with poverty, together with staffing held by something other than pay. It also shows what the next mechanism must avoid. It should sort households by income without drawing teachers' residences five times past the data, and without overshooting earnings and scores.

#### Round 2026-10-15c: income-graded sorting (Y) and teacher demand (D)

The round follows the [accepted development plan](#accepted-development-plan-2026-10-15).
- **Files.** Pre-registration (amendments 1 and 2), case lists, working log and the bounded-test script are in `runs/experiments-2026-10f/` (untracked). Outputs are in `runs/exp-1015c/`; final scores in `runs/experiments-2026-10f/scores_final.md`. Code snapshots, all in `runs/_code/`:
  - `r1015c1`: the fits.
  - `r1015c2`: the robust closure update; S1·Y·D00's fit, the Nz = 13 re-solves and the re-evaluations.
  - `r1015c3`: r1015c2 plus the review fixes; its bitwise and closure tests, the C4′ closure check and the smoke test.

  The main tree equals r1015c3; the pre-round tree is frozen in `runs/_code/pre-20261015c`, and r1015c3 is bitwise identical to it with the switches off.
- **Switches.** Both are off by default and bitwise-inert at C0, C4′ and S1·R (every moment key identical).
  - **Y**, `Params.By` with `lyref`, ϑ entry 28 `ζ_y`: the suburb's value rises by ζ_y per log point of the node's predicted pre-tax income ŷ. ŷ is constant within a node, so Proposition 2′ and the X_O* collapse hold.
  - **D**, `Params.fte_rule`: positions per pupil follow a rule, and b_T becomes the rationing price that `solve_ge` solves so the FTE gap equals the rule. With the rule at −0.0177 it reproduces S1·R (b_T −0.0593 both ways). D0 sets the rule to 0, so FTE becomes validation.
- **Designated criteria, added to A1–A8:** teachers' net commuting (+0.075, |z| ≤ 3), teacher crossing share (×½ to ×2 of 0.084), commuting directions (×½ to ×2) and the school-quality WTP elasticity (≤ 1.22). The per-employed-resident gap (+0.298) is reported, not designated. The identity net commuting = residents gap − FTE gap + Δlog(employed residents per pupil) puts the last term at −0.24 in the data, and the model, with one child per adult, cannot reproduce it.
- **Bounded test** (S1·R's aggregates fixed; ΔB, r_m, b_T, κ_c and c_C re-matched; `runs/exp-1015c/bounded/`). At V6a = −0.187:
  - **Y** needs ζ_y = 0.27. It leaves teachers where S1·R has them (net commuting −0.01, crossing 0.107), but V6z is only −0.05 and V4 0.05, and V1 overshoots (0.166).
  - **E** saturates at −0.167 and moves teachers' residences 4.6 times past the data.
  - **Z**, ability-graded, overshoots V4 (0.52).
  - Income grading sorts on income, not ability. Part of the V1 overshoot is a concept gap: V6a is measured on family income, so household structure and race are in the data's gradient but not in the model's.
- **GE fits.** All are exact at Nz = 9 (criteria below 1e−16) and pass A5's re-solve at Nz = 13 (largest move 0.08 SE). D cases run with the closure; † marks a targeted moment.

| | C4′ | S1·D0 | S1·Y | S1·Y·D0 | S1·Y·R | S1·Y·D00 |
|---|---|---|---|---|---|---|
| S_1015 (k) | 315.1 (10) | 341.8 (10) | 810.4 (11) | **151.6** (11) | 133.9 (12) | 153.9 (10) |
| ζ_y; b_T | —; 0 | —; −0.047 | 0.353; 0 | 0.231; +0.058 | 0.244; +0.051 | 0.229; −0.029 |
| V1 earnings 0.112 | 0.034 | 0.025 | 0.159 | 0.169 | 0.168 | 0.170 |
| V3 poverty −0.089 | −0.023 | −0.017 | −0.050 | −0.074 | −0.071 | −0.074 |
| V4 scores 0.289 | 0.082 | 0.025 | −0.182 | 0.196 | 0.152 | 0.187 |
| V6a −0.187 | −0.034 | −0.025 | † | † | † | † |
| V6z −0.222 | −0.007 | −0.006 | −0.054 | −0.043 | −0.045 | −0.043 |
| FTE −0.018 | +0.020 | 0 (rule) | −0.153 | 0 (rule) | † | 0 (rule) |
| growth 0.003 | 0.013 | 0.0025 | −0.052 | 0.021 | 0.012 | 0.019 |
| payroll −0.009 | 0.029 | 0.009 | −0.145 | 0.009 | −0.009 | 0.069 |
| T1 net commuting 0.075 | −0.023 | 0.006 | 0.019 | −0.044 | −0.037 | −0.024 |
| ε22 (≤ 3.7) | 13.2 | 3.62 | 3.87 | 3.70 | 3.71 | 3.70 |
| fails | baseline | A1, A2 | A1–A3, A6‼ | A3 (growth z +6.0), T1 (z −3.1) | A6‼ (3.715) | A3 (growth, payroll) |

Under D the FTE gap equals the rule by construction. Its z (+1.4) tests the common positions-per-pupil norm, and S1·D0's lower growth comes from setting FTE to 0 rather than to −0.018. S1·Y·D00's salary gap is not scored: +0.069 (z +6.5).

- **Selection (pre-registered).** No candidate meets every criterion, so C4′ stays the development candidate. The best-scoring candidate without a severe failure is S1·Y·D0: S_1015 151.6 against C4′'s 315.1. B1 stays sealed.
  - The outcome is close to a knife-edge. S1·Y·R scores lowest (133.9) and fails only A6, by 0.015 (ε22 3.715 > 3.7), which counts as a severe failure under the pre-registration.
  - C4′ stays by default, although it fails A6 itself (ε22 13.2).
- **Reading.**
  - *Y needs a demand closure.* Alone, it prices teachers out of the suburb (FTE −0.153, growth −0.052). PREREG expected the opposite sign, as under E.
  - *With D0 or R, Y matches the targeted V6a and, untargeted, most of V3.* V1 overshoots by 3.4 SE, and V6z stays at z ≈ +5 against the all-race target. The race concept is discussed below.
  - *Growth overshoots through teachers.* Of S1·Y·D0's 0.196 score gap, 0.081 is composition and 0.115 is school quality.
    - With ψ_a = 0 and FTE fixed by the rule, Q can differ only through teachers' human capital. Suburban teachers' is 0.10 log points higher; without Y the gap is at most 0.02, and SEDA growth allows about 0.015.
    - The likely route is that higher-h teachers have higher predicted income, live in the suburb and mostly work there. It is not isolated yet.
    - The same gap shows in pay: holding the salary gap needs δκ ≈ −0.066.
    - The closure's shadow price is positive (+0.058), not the expected −0.1 to −0.4.
  - *Pay.* S1·Y·D00 shows that, under Y, the suburb's salary gap would be 0.069 without δκ ≈ −0.066.
- **Numerics.**
  - The first closure (r1015c1) could lock into a 2-cycle (S1·Y·D00's start). r1015c2 replaces it with a robust update: the secant is used only when the aggregates have settled, the slope doubles on growing sign flips, and the slope is bounded at 48.
  - The r1015c1 fits re-evaluated on r1015c2 agree to 4e−6.
  - From a cold start at C4′ (∂FTE/∂b ≈ 27), r1015c3 converges in 57 and 50 iterations (rules −0.0177 and 0).
  - With the switches off, the default model at L = 1 and 2 is bitwise identical under r1015c3 and pre-20261015c (smoke test, `tests/smoke_l12_*.txt`).
- **Review.** An Opus review with adversarial verification confirmed four findings and refuted three. The confirmed ones are the closure cycle, an L = 1 crash, inert `--externals` and a missing ζ_y display. All are fixed in r1015c3.
- **Diagnostic after the round: race composition.** The sorting targets include between-race composition, and the model has no race. The round's working log has the tables, scripts and caveats (`runs/experiments-2026-10f/race_concept/`; they read no sealed file, compute no BA statistic, and were checked by a second review).
  - *V6z.* NLSY79's own-AFQT gradient is ≈ 0 among non-Black non-Hispanic youth (−0.020, SE 0.039) and −0.14 to −0.19 among Black and Hispanic youth. Pooled within race it is −0.05 with 1979 shares and −0.08 to −0.10 with 2018 shares: 22–45% of −0.222. With 2018 shares it still separates C4′ (≈ +2.3 SE) from the Y fits (≈ +1 SE).
  - *V4.* SEDA's score gap is ≈ 0 for white students and 0.16–0.21 SD for Black and Hispanic students. Pooled within race it is ≈ 0.10 (SE 0.023): about a third of 0.289. S1·Y·D0's 0.196 overshoots that by ≈ 4 SE.
  - *V5.* Growth is ≈ 0 pooled within race.
  - *Covariates.* Household income and child-poverty gaps (total-population covariates) keep about three-quarters of their size within race.

  Part of the V4 and V6z shortfall of a raceless income-sorting model is therefore one of concept. Race composition does not explain the growth overshoot, which misses by about 6 SE against either concept. S1·Y·D0's composition part (0.081) is close to the within-race composition the data imply (≈ 0.09–0.10). Whether race explains the V1 overshoot depends on a within-race V6a, which has not been measured; household structure is the more direct candidate. The pre-registered selection is unaffected; [§5](#next-round-proposed-after-2026-10-15c) has the proposal.

### First internal fit, overidentified (2026-09-29)

Historical reference for Table 5: [`runs/polish-noincome/`](runs/polish-noincome/), nine moments including FTE, criterion 5.15. Salary misses by −1.23 SE and FTE by +1.72 SE; earnings are −0.015 versus 0.112 and scores 0.007 versus 0.289. A second minimum near $r_m=0$ has criterion 28.6. Removing FTE from the criterion yields the single root found above, with a larger FTE validation miss. Detailed reports and search journals remain in `runs/`.

### Table 5: Sensitivity panel

The original 40-case panel uses the **overidentified first fit**, including FTE (`runs/panel/`, `slurm/panel_cases.tsv`; [results](../../../data/spatial/estimates/sensitivity_panel.md)). Only $\sigma$, $K$, location-taste and WTP variants have been rerun at the exact baseline (`runs/panel-exact/`). Do not compare criteria across different target sets.

| Experiment | Main finding at the first fit |
|---|---|
| Moving-cost profile | Criterion 28.7 at $r_m=0$, 5.15 near 0.248, 5.44 at 0.30. Mobility identifies moving costs; income does not. |
| Location-taste scale, 1/7–0.50 | $r_m$ ranges roughly 0.14–0.44, with similar fit. The cost level depends on assumed taste dispersion. |
| Warm glow, score bridge, ability, persistence, investment, WTP and CFR variants | Criteria roughly 5.1–5.4; validation misses persist. WTP changes $\lambda$, CFR changes $\beta$. First-stage refits: [report](../../../data/spatial/estimates/first_stage.md). |
| Earnings targets (FTYR, all-earner, household) | Criteria 54, 63 and 112; model earnings gap remains near zero. |
| Class-size curvature and exposure | At $\sigma=0.30,0.56,0.85$, $\lambda\approx0.036,0.124,0.49$ and $r_m\approx0.89,0.25,0.06$. $\sigma=0.15$ and $K=5$ hit $r_m=1$. Some parameter movement reflects frozen reference units rather than different outcomes. |
| ACS 2016–19 occupational block (T7) | Pooled 90/10 rises to 3.82 and $\sigma_\epsilon$ to 0.806; refit criterion 5.03, validation misses unchanged. This refit also predates the exact baseline. [Report](../../../data/spatial/estimates/first_stage_t7.md). |

For any remaining exact-baseline reruns: explicitly activate FTE when targeted; dropping mobility leaves seven moments for eight parameters unless another restriction is supplied; use a separate `PANEL_OUT`. The baseline grid checks do not establish numerical accuracy of every extended specification.

## 5. To-do

### Next round (proposed after 2026-10-15c)

This is a proposal. It needs the user's decisions on how the model treats race, and on whether to model the group dimension. The accepted plan's preservation rules apply.

1. **Decide how the model treats race before adding mechanisms.** Every sorting target mixes in between-race composition, and the model has no race ([round 2026-10-15c](#round-2026-10-15c-income-graded-sorting-y-and-teacher-demand-d)). A raceless model with the full income distribution should produce the part of the between-race location differences that runs through income, but not the rest: segregation, wealth, discrimination and preferences.
   - **(a) Race-consistent re-scoring: a diagnostic and data work only; first.**
     - Re-measure V1, V3, V4, V5, V6a and V6z under a race-consistent concept, on the same partition, sample and bootstrap. Prefer a concept that removes only the between-race differences not explained by income, for example by reweighting groups to a common income distribution. Pooled within-race gaps are a bound. Sources: SEDA subgroups, ACS district tables by race, NLSY79 by race.
     - Register every construction choice before computing V1, V3 and V6a by race: population shares (1979 or 2018, national or CZ), groups, tercile convention, `wa` or `wp`, CZ sample and weights, and the income concept (family income or individual earnings, i.e. household structure). The registration must state that within-race V4, V5 and V6z have already been seen.
     - Re-score every fit (C4′, S1·R·E and all round-c candidates) under both concepts in one permanent table. This needs no new runs.
     - This is not a replacement target set unless the first stage (ability, wage dispersion, teacher shares) is re-estimated on the same concept. Any selection under new targets is a new round.
   - **(b) Model the group dimension: likely needed for this model's purpose.**
     - Teacher mobility responds to student composition (Hanushek, Kain and Rivkin 2004; Jackson 2009). Without the group dimension the model absorbs that in b_T, and pay or composition counterfactuals inherit the error.
     - A minimal version adds an exogenous household type whose location composition comes from the data: white students are 29% of city and 54% of suburban assessments. Teachers' workplace value may then depend on composition.
     - Frame group differences as resources and circumstances, not as the model's inherited ability.
   - *Recommendation:* run (a) first, because it is cheap and shows how much of each sorting miss is concept. Then decide on (b).
2. **Then teacher-quality sorting.** Growth misses by about 6 SE against either race concept.
   - First isolate the route: add a teacher-human-capital-by-residence moment and check ζ_y switched off for teachers.
   - Then test a hiring margin that limits how far teacher quality follows teachers' residence, for example filling positions from a common applicant pool by credentials or seniority.
   - Discipline it with teacher-quality gaps by locale, such as the CRDC shares of novice and uncertified teachers.
3. **Keep the demand closure.** Income-graded sorting without it prices teachers out of the suburb (S1·Y). D0 and R give the same sorting; D0 leaves FTE as validation.
4. **B1.** It stays sealed.
   - When the user instructs it to be opened, it is compared as sealed with the model's `ba_share_gap`, as an all-race moment.
   - Under race-consistent targets, its between-race component is a stated caveat.
   - A within-race B1 needs the user's instruction and a definition fixed before opening.
5. **Unchanged.** M3 with γ < 1 and S1; fits at Nz ≥ 9 with Nz = 13 checks; non-homothetic housing (CGHH-style quality tiers) as the structural successor to the reduced-form ζ_y once the concepts are settled.

### Accepted development plan (2026-10-15)

**Decision and scope.** After reviewing both October 15 rounds, the user provisionally accepted the recommendations on condition that previous work remains preserved. Development can proceed without another general approval request. This accepts a research direction, not a fitted specification: retain the September exact baseline, C4′ as the designated development candidate and S1·R·E as a diagnostic benchmark. The historical experiment scores and selection decisions stand. This documentation update launches no model runs.

**Preservation is a requirement of the accepted plan.** Before this update, current source, configuration, spatial notes and selected fit evidence were copied into the verified local archive [`runs/_code/accepted-direction-20261015/`](runs/_code/accepted-direction-20261015/README.md). It contains 505 files with a SHA-256 manifest, the Git HEAD and working-tree status, including uncommitted source changes. The archive is Git-ignored and local to ResearchDrive; it is not a separate off-site backup or a complete optimizer/raw-data checkpoint. All existing run directories, empirical inputs, preregistrations, working logs and earlier code snapshots remain in place. Sealed hold-out contents were not read or archived.

- Never overwrite historical fitted parameters, moments, diagnostics, case lists, protocols or snapshots. Use a new named round and separate output directories; keep reference units and input vintages explicit.
- Add new mechanisms through configurations/switches that retain existing behavior when disabled. Verify the September baseline and C4′ configurations before substantive new runs, and preserve S1·R·E for comparison.
- Snapshot the new round's code, configuration and input provenance before estimation. Keep raw/generated/hold-out ignore rules and the untracked status of local `AGENTS.md` intact. Do not reset, clean or discard existing uncommitted work.

**Next sequence.**

1. **Write a small comparison before running it.** Specify the economic mechanisms, parameter-to-moment map, candidate configurations, scoring components and numerical acceptance thresholds. The main comparison should distinguish the contribution of income-dependent sorting from explicit teacher demand, then evaluate their combination. Keep new parameters tied to independent evidence where available; any outcome used to estimate a parameter must be labelled targeted.
2. **Retain the successful education block and repair supply.** Carry forward M3's college/occupation choice and licensure, with estimated $\gamma<1$, and use S1's commuting-friction repair in new candidates. Preserve the original C4′ configuration for comparison. Report own-location, cross-location, ratio and small-district elasticities with their adjustment margins; matching all-worker crossing alone is insufficient.
3. **Test income-dependent residential sorting.** Prioritize non-homothetic housing or amenity demand, meaning expenditure shares vary with income, rather than an amenity exclusively for BA holders. Discipline it using housing budget shares and rental distributions, recognizing that spending is not a quality-adjusted price. A bounded calculation with explicitly fixed aggregates should precede full equilibrium. The failed uniform charge and bounded-access test do not reject housing quantities/endogenous rents. Do not try to obtain sorting simply by raising school-quality WTP, already too income-elastic. Additional broad Roy, contemporaneous peer-effect or persistence searches are lower priority.
4. **Replace the staffing reduced form with teacher demand.** Test district budgets, authorized positions, hiring/rationing rules or class-size constraints that jointly explain salary and FTE per pupil. Define how available positions are filled and how selection affects school quality. Keep $b_T$ as a diagnostic comparator, not an unexplained permanent substitute for demand.
5. **Strengthen joint acceptance criteria.** Require explicit checks of teacher residence (+0.298 log gap per employed resident), teacher crossing (0.084), net commuting (+0.075) and both commuting directions, with sample/geography limitations recorded. Evaluate these alongside earnings, parental-income and own-AFQT location gradients, poverty, staffing, teacher selection, score levels and growth. Retain own-AFQT mobility, the gender wage gap, MA premium and non-BA schooling checks. Fix thresholds before new results; a lower aggregate score cannot offset a severe failure on the designated criteria. Report target status, component discrepancies and added parameters so fitted moments are not presented as independent validation.
6. **Resolve numerics and identification before inference.** Estimate strongly sorted candidates at $N_z\ge9$ and re-solve at a finer grid, initially $N_z=13$, using a pre-specified fixed-parameter tolerance. Check equilibrium residuals, multistarts, Jacobians and profiles for weak parameters. After choosing a satisfactory specification, align samples/payroll definitions, recheck the first stage in equilibrium, and complete the full calibration bootstrap with Table 2 re-estimated and literature-input sensitivities. Parallel or potentially long computations go through SLURM, with refreshed ResearchDrive credentials, existing log directories and recorded jobs/configurations/output paths.

**Hold-out.** B1 stays sealed for evaluation of the structural candidate. Acceptance of this plan is not an instruction to open it; explicit user instruction is still required. Previously used NLSY97 evidence is not a fresh hold-out, and overlapping ACS vintages are not automatically independent validation.

**Readiness.** Clear to proceed with design and implementation under this accepted plan. The next deliverable is the mechanism specification and pre-specified comparison, not another commuting audit or immediate adoption of S1·R·E. Exact parameterizations and numerical thresholds remain implementation work; their absence does not require repeating the user's general authorization. New baseline adoption follows the joint empirical and numerical assessment.

### Recommended sequence (2026-10-15)

**Historical sequence, before the October 15 rounds.** Retained to explain the completed experiments; the accepted development plan above now sets priorities. Keep `runs/exact-polish-base/` as the reference and develop from C4′ (`runs/exp-c3/c4p/`). C4′ substantially improves teacher pay and selection, but still gives earnings 0.034 against 0.112, scores 0.082 against 0.289 SD, and growth 0.0127 against 0.0031 SD/grade. Further optimization of the already exact baseline is not the next priority. The milestone is household sorting that coexists with credible teacher allocation and small learning-rate gaps.

1. **Audit commuting and teacher supply before changing the block.** The current `teacher_supply_elasticity` in [spatial_moments.jl](spatial_moments.jl) computes $\varepsilon_{\rm ratio}=\partial\log(n_2/n_1)/\partial\log\kappa_2$, holding $\kappa_1$, teacher human-capital aggregates, population, taxes, school quality, position weights and the stationary distribution fixed. Here $n_w$ is teacher mass by workplace. Household occupation, residence and workplace choices adjust; under M3 education also adjusts. Thus $\varepsilon_{\rm ratio}=\varepsilon_{22}-\varepsilon_{12}$, where $\varepsilon_{w2}=\partial\log n_w/\partial\log\kappa_2$. The saved value near 31 is not the own-location elasticity $\varepsilon_{22}$.

   [Falch (2010)](https://www.journals.uchicago.edu/doi/abs/10.1086/649905) estimates supply to individual Norwegian schools at about 1.4; [Ransom–Sims (2010)](https://www.journals.uchicago.edu/doi/abs/10.1086/649904) infer supply to Missouri school districts from wage responses of quits, obtaining about 3.7. Before using either as a restriction, report the model's own-location, cross-location and ratio responses separately, distinguish fixed-education/residence responses from those allowing these choices to adjust, and report full-equilibrium responses separately. Align the pay change, geographic aggregation and adjustment horizon. The code comment interpreting the employer estimates as upper bounds also needs this justification.

   Audit both commuting directions and teacher-specific crossing shares alongside the all-worker share currently targeted. The [ACS flow report](../../../data/spatial/estimates/settings/commuting_all.md) uses POWPUMA classes: within-county city/suburb crossings can be invisible, so measured crossing shares are lower bounds. Reconcile public versus all K–12 teachers and sample ages across flow and staffing moments. If excessive responsiveness remains after this audit, test workplace attachment or heterogeneous commuting costs with independent identifying evidence. Keep the staffing closure ($\delta_\kappa$ ↔ FTE, salary as validation) as a sensitivity case; changing the closure does not repair the supply mechanism.

2. **Test heterogeneous housing access next.** The [housing report](../../../data/spatial/estimates/settings/housing_transport.md) gives suburban gross-rent gaps of 0.198 log points at the 10th percentile and 0.023 at the 90th. This motivates a test of limited inexpensive suburban housing through quantities, minimum bundles or availability, with endogenous rents. These spending distributions do not identify quality-adjusted prices; the mechanism remains a hypothesis. Failure of the uniform residence charge does not reject it.

   Discipline housing parameters with rental distributions and housing budget shares. Start with a bounded calculation holding school quality fixed and state which other aggregates remain fixed; ask whether empirically supported parameters produce meaningful income and ability sorting. If they do, proceed to full equilibrium and joint recalibration. Compare with the existing $\zeta$ benchmark, retaining its diagnostic status. If a sorting moment identifies a new parameter, label it targeted and assess the remaining sorting outcomes as validation. If housing access cannot deliver the required sorting, parental investment/shared family inputs remain a subsequent candidate, motivated by the sibling-correlation miss.

3. **Carry forward M3 and estimated wage compression.** Retain college/occupation choice, licensure and the residence/workplace split in the candidate, with $\gamma$ estimated rather than fixed at 1. M3 improves teacher schooling and pay and produces positive overall selection with negative selection within college graduates. The completed regular-Roy and contemporaneous peer-effect experiments do not justify another broad round of either: the former worsens dispersion and selection under C4′, while the latter adds little sorting and raises growth once sorting is imposed. Keep the market gender gap, MA premium and non-BA teacher schooling as unresolved checks.

4. **Pre-specify acceptance criteria.** Fix the candidate set, identifying moments, validation set and numerical checks before inspecting new results. Require joint improvement in earnings, parental-income and AFQT location gradients, staffing, and both score levels and growth while preserving teacher-selection gains. Retain poverty as a check with its provisional measurement caveat. Report component discrepancies alongside the aggregate score so that improvements in teacher pay cannot conceal failures elsewhere. Every new free parameter needs an explicit source of identification, and a moment promoted to a target must lose its validation label. Reserve genuinely unused evidence: the earlier NLSY97 holdout is spent, and overlapping ACS vintages provide stability checks rather than automatically independent validation.

5. **Complete calibration and inference after specification choice.** Align occupational and spatial sample definitions and build the like-for-like teacher-payroll comparator. Recheck first-stage moments inside the chosen equilibrium and revisit the occupational wage block where discrepancies persist. Run the extended specification's own grid, equilibrium-residual, multistart and identification checks; the baseline's numerical checks do not establish its accuracy. Then complete the full calibration bootstrap, re-estimating Table 2, and literature-input sensitivities. Use SLURM for parallel or potentially long workloads, refresh ResearchDrive credentials before submission, and record job IDs, commands, configurations and outputs.

**Immediate deliverable:** the commuting-elasticity audit, deciding whether teacher supply needs repair before the housing-access experiment. Adoption of a new baseline follows the joint validation assessment, not this documentation update.

**Status after the 2026-10-15 rounds** ([§4](#model-experiments-2026-10-15)):
1. *Supply: done.* S1 repairs teacher supply (full-equilibrium own-location elasticity 1.6). No fit reproduces the direction of non-teachers' commuting. Under S1, staffing needs a closure that pay does not provide.
2. *Housing access: done.* It fails its bounded gate. Persistence bounds shared family inputs from above (ρz 0.90 moves V6a by −0.004), so they are not the next candidate either.
3. *Carry forward M3 with γ < 1: kept.* γ ≈ 0.61 in every fit.
4. *Pre-specified criteria: done.* A1–A8, $S_{1015}$ and the sealed B1 are in place. B1 is still sealed.
5. *Inference after specification choice: not reached,* because no candidate was selected. For the best candidate, S1·R·E, the grid, identification and multistart checks were run:
   - The grid check fails at fixed ϑ but holds after a re-fit.
   - The condition number is 2,320, with ζ_BA weakly identified.
   - The multistart finds a unique solution.
   - The bootstrap and the in-equilibrium first-stage recheck remain open.

**Recommendations at round close (subsequently accepted and expanded above).**
- **A sorting mechanism graded by income, not education.** S1·R·E shows the sorting and staffing targets can be met together once staffing is held by something other than pay. But a BA-graded amenity overshoots earnings and scores, moves teachers' residences fivefold, and is weakly identified.
  - Proposition 4 implies that log utility with proportional location premia does not sort. A non-homothetic location demand is the structural candidate, for example income-elastic demand for housing quantity or for a suburban amenity.
  - Discipline it jointly by V1, V6a and V6z, with the teacher-residence moments as validation (+0.298 per employed resident; net commuting +0.075).
  - It must not raise the school-quality WTP income elasticity, already 1.14 against BFM's 0.38.
- **A teacher-demand side** to replace b_T. In S1·R·E, b_T ≈ −0.35 is worth about 30 log points of pay. District hiring rules, budgets or class-size floors would make the staffing closure testable against FTE and salary together.
- **Finer ability grids for sorted candidates.** Score any candidate with a large composition gap at Nz ≥ 9 (S1·R·E: 0.41 SD), and check it at a finer grid still.
- **B1.** Keep it for the structural mechanism unless the user instructs otherwise. S1·R·E predicts +0.57.

### Additional model experiments (2026-10-07)

Historical task list; the sequence above now sets the order of remaining work.

1. **Regular Roy model.** Solve and calibrate the regular Roy model with $\gamma = 1$. *Done* ([§4](#model-experiments-2026-10-07)): exact fits of C0 and C4′, neither adopted.
2. **Local peer effects.** Add average ability in location $l$, $\bar a_l$, raised to a power $\psi$ as a multiplicative component of $Q_l$ in the human capital production function: for example, $h = \bar a_l^{\psi} z \epsilon \cdots$. The current specification sets $\psi = 0$. *Done* ([§4](#model-experiments-2026-10-07)): implemented as $\psi_a$ (`Params.ψa`; $\psi$ is the warm glow). It creates no sorting at the fitted λ, and once sorting is imposed SEDA growth bounds $\psi_a\lesssim0.06$.
3. **Teacher commuting.** Allow teachers to commute across locations so their residence can differ from their workplace. Build on the residence/workplace prototypes in §4; the reference baseline still requires teachers to live where they work.

### Next model experiments (2026-10-05)

Historical task list; the 2026-10-15 sequence above refines the commuting comparison and prioritizes housing access.

The 2026-10-03 list (income-driven demand, residence vs. workplace, joint education and occupation) is done: see §4, 2026-10-05. Household sorting by income and ability is now the binding problem. It replicates in NLSY97, and the ζ diagnostic shows that the extended teacher block can absorb a sorting force.

1. **A sorting mechanism.** It must deliver:
   - earnings ≈ 0.11, poverty ≈ −0.09 and city-share gradients of −0.2 to −0.3;
   - small quality and growth gaps;
   - neither a higher WTP income elasticity (already 3× BFM) nor a uniform price wedge above the measured ≈0.04 $\bar C$.

   Candidates:
   - Income-graded access with endogenous rents. Suburban rents exceed city rents by 0.20 log points at the 10th percentile but 0.02 at the 90th: the suburb lacks cheap rentals.
   - Peer composition in school quality or children's outcomes. This gives self-reinforcing sorting at a small $\Delta\log Q$; check it against SEDA growth. Tested in school quality on 2026-10-07 ([§4](#model-experiments-2026-10-07)): no self-reinforcing sorting at the fitted λ, and growth bounds the elasticity.
   - Parental investment and shared family inputs (sibling correlation 0.33 against 0.61).

   Compare each with the ζ benchmark row.
2. **Teacher supply elasticity.** The saved staffing-ratio diagnostic gives 31 at the fitted $\rho_w$, and ≈4 only with no crossing or 46% crossing. The original comparison to 1.4–3.7 treated this as an own-employer elasticity; audit the mapping before concluding that the commuting and supply evidence cannot be matched. The [2026-10-15 recommendation](#recommended-sequence-2026-10-15) makes workplace attachment conditional on that audit and retains the staffing closure as sensitivity rather than the default.
3. **Wages under M3.** Market gender gap 0.36 against 0.10; MA premium 0.07 against 0.16; non-BA teachers' schooling 0.37 against 0.52; BA teachers within BA −0.12 against −0.24.
4. **Level comparators.** Build a like-for-like teacher-payroll share (model ages and income definition) and check the tax levels against it.

**Evaluation.** Keep the protocol: pre-registration, fixed scored set, `score_candidates.py`, and sealed data opened only after selection. The pre-registered NLSY97 sealed moments are now used, and the other NLSY97 statistics were partly displayed. The next round needs a fresh holdout, for example district moments from another ACS vintage, or the teachers' wage–AFQT slope once computed.

**Other remaining work:**
- The full calibration bootstrap with Table 2 re-estimated, and literature transport uncertainty. The delta-method tool covers targets and externals at $k=1$.
- Exact-baseline sensitivity reruns where informative.

Configuration, first-stage estimation, target construction, moment implementation and the historical T7 vintage check are complete; their sources and results are linked above.
