# Calibration of the spatial teachers model

The baseline is the exactly identified September 29, 2026 calibration in [`runs/exact-polish-base/`](runs/exact-polish-base/). It matches eight targeted moments and passes numerical and local-identification checks, but falls short on household income sorting, teacher allocation and several teacher labor-market outcomes. Subsequent experiments improve teacher schooling and pay without resolving the spatial sorting failures; none has replaced the baseline.

## Baseline strategy

The [model](../spatial_continuous.jl) ([writeup](../../../notes/spatial_teachers_model_eps0.md)) represents a city and suburb within a commuting zone (CZ). Spatial gaps are enrollment-weighted within-CZ contrasts, suburb minus city, in logs unless stated otherwise. District inputs are from 2018, SEDA outcomes pool 2009–19, and occupations and teaching shares use ACS 2009–13.

Calibration proceeds in three stages: fix aggregate and literature inputs (Table 1); estimate occupations and ability outside equilibrium (Table 2); then jointly estimate eight internal parameters against eight equilibrium moments (Table 3). Recheck the first stage in the fitted equilibrium and evaluate untargeted outcomes separately (Table 4). This separation preserves income sorting, staffing and achievement as tests of the model rather than additional fitting targets. The [calibration log](calibration_log.md) provides definitions, derivations, implementation details and the experimental history.

## Table 1: Fixed inputs

| Parameter | Value | Basis |
|---|---|---|
| Ability elasticity $\alpha$; consumption weight $\mu$ | 1; 1 | Normalizations |
| Population $M$; amenity $B_1$; home productivity $A_{HP}$ | 2; 0; 1 | Normalizations |
| Goods-investment elasticity $\eta$ | 0.080 | 2018 education spending net of teacher payroll, relative to labor income (NCES/PWT) |
| Time-investment elasticity $\phi$ | 1.0584 | Non-teacher schooling: 13.64 years/25 at the reference equilibrium (ACS) |
| Class-size curvature $\sigma$ | 0.56* | Class-size wage effect, FOO (2013) |
| Location taste scale $\sigma_\nu/\mu$ | 0.25* | Eckert–Kleineberg (2024) |
| Warm-glow curvature $\psi$ | 0 | Linear warm glow |
| Earnings–achievement bridge $\chi$ | 0.113 | $\log 1.12$, CFR (2014b); school-driven score changes only |
| School exposure $K$ | 12 years | Grades 1–12, equal yearly exposure |

\*Provisional literature inputs. Reference consumption $\bar C$ and human capital $\bar h$ fix moving-cost and warm-glow units. Sources and calculations: [log, section 1](calibration_log.md#1-external-inputs).

## Table 2: Occupations and ability

The model has 19 market occupations, home production and K–12 teaching. The first stage holds school quality and taxes fixed and ignores teaching selection; it matches 41 moments with 41 parameters. There are no education wedges or wage wedges for men or teachers.

| Parameter | Target | Estimate |
|---|---|---|
| Relative productivities $A_i/A_{HP}$, 19 | Men's non-teaching occupation shares, ACS | $\log(A_i/A_{HP})\in[-1.53,-0.33]$ |
| Relative female wedges $r_i$, 19 free | Women's non-teaching occupation shares, ACS | $\log r_i\in[-1.58,0.44]$ |
| Occupational ability dispersion $\sigma_\epsilon$ | Weighted non-teacher within-cell hourly-wage 90/10: 3.737, ACS | 0.786 |
| Ability persistence $\rho_z$ | Latent mother–child score correlation: 0.577 (SE 0.014), NLSY79/CNLSY | 0.577 |
| Stationary ability dispersion $s_z$ | Log-wage slope per latent score SD: 0.213 (SE 0.008), NLSY79 | 0.196 |

Shares condition on not teaching, including home production. Wage ratios use occupation-by-gender wage-sample weights. The NLSY slope pins $s_z$ independently of $\chi$; $\chi$ converts school-driven achievement changes. Full occupational estimates: [first-stage report](../../../data/spatial/estimates/first_stage.md).

## Table 3: Internal calibration

Write $\kappa_1=\bar\kappa=\tilde\kappa\bar h^{1-\gamma}$, $\kappa_2=\bar\kappa e^{\delta_\kappa}$, $B=(0,\Delta B)$, $m_{12}=m_{21}=r_m\bar C$ and $1-\tau^\omega_{i,f}=\omega_f r_i$. Rows indicate the main identifying moments; estimation is joint.

| Parameter | Value | Target | Source | Data | Model |
|---|---|---|---|---|---|
| Teacher pay scale $\log\tilde\kappa$ | −1.2038 ($\bar\kappa=0.123$) | Male teaching share | ACS | 1.95% | 1.95% |
| Female outside option $\omega_f$ | 0.9664 | Female teaching share | ACS | 5.99% | 5.99% |
| Teacher wage curvature $\gamma$ | 0.6194 | Teachers' hourly-wage 90/10 | ACS | 2.773 | 2.773 |
| Spatial pay gap $\delta_\kappa$ | 0.0045 | Adjusted salary gap | Audited districts | 0.0086 (SE 0.0092) | 0.0086 |
| Amenity gap $\Delta B$ | 0.0662 | Enrollment gap | Districts | 0.3764 (SE 0.1154) | 0.3764 |
| Teacher spillover $\beta$ | 0.2278 | Earnings effect per teacher-quality SD-year | CFR (2014b) | 1.34% (SE 0.41 pp) | 1.34% |
| Altruism $\lambda$ | 0.1245 | WTP per school-average score SD, share of consumption | BFM (2007) | 0.595% (SE 0.224 pp) | 0.595% |
| Goods moving cost $r_m$ | 0.206 | City↔suburb move rate | NLSY79 | 0.216 (SE 0.011) | 0.216 |

Teaching-share denominators include home production; the teacher 90/10 averages gender-specific ratios with wage-sample weights. District targets cover 188 CZs. CFR maps to $\beta\,\mathrm{sd}(\log h_T)/[K(1-\eta)]$; BFM maps to compensating consumption, with one school-average SD taken as 0.45 student SDs. [Definitions and sources](calibration_log.md#2-internal-calibration).

The fit matches all eight targets (criterion $1.6\times10^{-19}$). Fifteen multistart searches find no other root, finer grids preserve the fit, and the scaled moment Jacobian has full rank (condition number 30.3). These checks support numerical stability and local identification; they do not establish the model's fit to untargeted outcomes. [Estimation and identification diagnostics](calibration_log.md#exactly-identified-fit-2026-09-29-the-baseline).

## Table 4: Untargeted validation

| Outcome | Data | Model at the fit |
|---|---|---|
| Median full-time year-round earnings gap | 0.1122 (SE 0.0168) | −0.0102 |
| K–12 teacher FTE per pupil gap | −0.0177 (SE 0.0131) | +0.0234 |
| Child-poverty gap | −8.88 pp (SE 0.90) | ≈+0.2 pp; provisional poverty line |
| SEDA score-level gap | 0.2887 SD (SE 0.0327) | 0.061 SD |
| SEDA learning-rate gap | 0.0031 SD/grade (SE 0.0030) | 0.0119 SD/grade |
| Move rate by parental-income tercile | 0.221, 0.226, 0.207 | ≈0.20, 0.22, 0.23 |
| Local revenue-share gap | 0.4726 (SE 0.0206) | Diagnostic only |

The first-stage consistency checks pass, but the baseline does not reproduce household income sorting: it predicts a slightly negative suburban earnings gap and almost no poverty gap. Matching the salary gap puts more teachers per pupil in the suburb, opposite to the data (3.1 SEs above the observed FTE gap). The score-level gap is too small, while learning-rate growth is too high (2.9 SEs). Earnings and staffing remain validation outcomes, not targets used to identify moving costs or spatial pay.

Teacher selection and labor-market outcomes also remain weak. Teacher/non-teacher schooling is 0.76 versus 1.18 in the data; relative teacher pay is too low; and the market gender log-wage gap is 0.30 versus 0.10. Occupational wage levels and dispersion remain additional diagnostics. [Known issues](calibration_log.md#known-issues-at-the-current-fit).

**Sensitivity.** The class-size, exposure, location-taste and WTP cases have been rerun at the exactly identified baseline and change no validation conclusion; the rest of the panel and the ACS 2016–19 refit have not. [Table 5 details](calibration_log.md#table-5-sensitivity-panel).

## Additional mechanisms

The accepted scope for extensions comprises commuting, local housing costs, fertility/parenthood shocks, education-based teacher salary schedules, teacher education and licensure requirements, and peer effects. These mechanisms remain experimental additions to the baseline. Their implementation and magnitude require independent empirical or institutional support; improved fit alone is insufficient.

The October 19 comparison refitted all 64 combinations to the baseline's eight targets, adding the all-worker crossing share as a target when commuting is enabled. Extension magnitudes were fixed externally, apart from the commuting friction fitted to that share. All combinations retained exact target fit and passed finer-grid checks. Their effects on calibration are:

| Mechanism | Calibration implication |
|---|---|
| Separate residence and workplace choices | Adds a commuting friction identified by the crossing share (10.7%). Raises spatial income and score gaps mainly through school quality, but over-staffs suburban schools and overstates learning-rate growth. |
| Local housing costs | A suburban access charge of 2% of reference consumption generates modest income sorting and offsets some over-staffing when combined with salary schedules. Results depend on the charge; housing expenditure gaps are not quality-adjusted price gaps. |
| Fertility/parenthood shocks | At a parenthood probability of 0.802, assuming no school valuation among childless adults raises altruism from 0.125 to 0.155 to preserve all-household WTP. Other baseline parameters and outcomes barely change; tested ability gradients also have small effects. |
| Education-based salary schedules | An externally set schooling-pay elasticity of 1.74, based on the observed MA salary premium, raises relative teacher schooling to 1.20 (data: 1.18) and substantially reduces the teacher-pay shortfall. |
| Education and licensure requirements | A BA entry floor improves schooling and pay as an alternative to the salary schedule. It does not bind under the tested schedule, so combining the two adds no benefit. Observed mean schooling is not an entry requirement. |
| Peer effects | At the tested magnitude, effects are negligible because the model generates little ability sorting. Peer effects do not initiate residential sorting; much stronger effects and alternative peer inputs amplify or create school-quality gaps but do not sort households by ability ([below](#peer-effect-strength-and-specification)). |

**Strongest result.** Housing costs plus the salary schedule are the simplest combination near the best untargeted fit: the sum of squared standardized errors across 13 validation moments falls from 864 to 364. Most of the improvement comes from teacher pay. Income and ability sorting, teacher selection and the gender wage gap remain unresolved. At the main externally fixed magnitudes, no combination achieves even one-third of the observed income-sorting index; larger housing charges strengthen sorting but make the interpretation and empirical support for the charge consequential. [Comparison, sensitivity checks and evidence](calibration_log.md#factorial-round-the-six-eligible-mechanisms-in-every-combination-2026-10-19).

### Untargeted fit of the leading combinations

**HS** combines housing costs ($r_p=0.020$) and the salary schedule (schooling-pay elasticity 1.74). **HSP** adds peer effects ($\psi_a=0.10$) and has the lowest S13; HS is the simplest combination within 10 points of it. Both match the eight calibration targets exactly, so every outcome below is untargeted; the wage slope and mother–child correlation also check inputs estimated in the first stage. The baseline column is the no-extension refit on the factorial's grid ($N_z=9$). Spatial gaps are suburb minus city, in logs unless noted. City-share gradients are top-minus-bottom tercile differences, by parental income or own AFQT.

| Outcome | Data | Baseline | HS | HSP |
|---|---|---|---|---|
| Median full-time year-round earnings gap | 0.1122 | −0.0102 | 0.0122 | 0.0122 |
| Child-poverty gap (pp; provisional model poverty line) | −8.88 | +0.24 | −1.18 | −1.18 |
| SEDA score-level gap (SD) | 0.2887 | 0.0605 | 0.0685 | 0.0686 |
| City-share gradient by parental income (ACS within-CZ) | −0.1870 | +0.0108 | −0.0123 | −0.0123 |
| City-share gradient by own AFQT (NLSY79) | −0.2220 | +0.0026 | −0.0024 | −0.0024 |
| K–12 teacher FTE per pupil gap | −0.0177 | +0.0234 | +0.0226 | +0.0224 |
| SEDA learning-rate gap (SD/grade) | 0.0031 | 0.0119 | 0.0117 | 0.0117 |
| Teacher payroll per pupil gap | −0.0091 | +0.0320 | +0.0312 | +0.0310 |
| Male teacher log-wage premium over market mean | +0.0557 | −0.7195 | −0.1490 | −0.1490 |
| Female teacher log-wage premium over market mean | +0.0860 | −0.6311 | −0.0497 | −0.0497 |
| Market gender log-wage gap | 0.1046 | 0.3024 | 0.3281 | 0.3281 |
| Log-wage slope per latent score SD | 0.2127 | 0.2138 | 0.2058 | 0.2058 |
| Mother–child latent score correlation | 0.5766 | 0.5760 | 0.5761 | 0.5761 |
| **S13: sum of squared standardized errors, above 13 rows** | — | **864.13** | **364.11** | **363.96** |
| Teachers' latent ability relative to market workers (SD) | +0.614 | −0.287 | −0.225 | −0.225 |
| Teacher/non-teacher schooling ratio | 1.176 | 0.761 | 1.203 | 1.203 |
| Move rates by parental-income tercile (low, middle, high) | 0.221, 0.226, 0.207 | 0.205, 0.216, 0.226 | 0.206, 0.216, 0.225 | 0.206, 0.216, 0.225 |

S13 standardizes by empirical SEs where available (score-level SE 0.0326, against 0.0327 in Table 4) and by a 0.05 fit scale for each teacher premium and the gender gap. The last three rows are outside S13. Sources: saved moments for the [baseline](runs/exp-factorial/none/moments.toml), [HS](runs/exp-factorial/HS/moments.toml) and [HSP](runs/exp-factorial/HSP/moments.toml), and the [factorial scores](runs/factorial-20261019/scores.md).

The schedule brings relative schooling close to the data and removes most of the teacher-pay shortfall, but teachers remain less able than other market workers, not more. Income sorting improves only modestly, staffing and payroll gaps keep the wrong sign, and the gender wage gap worsens. Peers lower S13 by just 0.15.

### Peer-effect strength and specification

The October 25 round asks whether peer effects deliver the missing sorting when they are much stronger or use a different peer input. School quality is $Q_l=Q_{T,l}P_l^{\psi_a}$, where $Q_{T,l}$ is the teacher index and $P_l$ is one of four peer inputs, each relative to an average location:

- **A** (the existing form): students' mean ability;
- **B**: students' total ability, their number times their mean;
- **C**: adult residents' mean human capital;
- **D**: adult residents' total human capital.

At each fixed $\psi_a$, the eight baseline parameters were refitted exactly to the eight targets at $N_z=9$, from the baseline fit and from HS, so every outcome below is untargeted. SI, the income-sorting index, averages the earnings gap, the child-poverty gap and the parental-income city-share gradient, each as a fraction of its data value. The children's z gap is the gap in students' mean log ability.

| Fit | S13 | SI | Earnings gap | Score gap (SD) | FTE/pupil gap | Growth gap (SD/grade) | Children's z gap (SD) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Data | — | 1 | 0.112 | 0.289 | −0.018 | 0.0031 |  |
| Baseline | 864.1 | −0.06 | −0.010 | 0.061 | 0.023 | 0.0119 | −0.0050 |
| A: mean ability, ψa 32 | 971.4 | −0.11 | −0.018 | 0.017 | 0.080 | 0.0041 | −0.0059 |
| C: mean adult h, ψa 1 | 885.8 | −0.07 | −0.012 | 0.049 | 0.038 | 0.0098 | −0.0052 |
| B: total ability, ψa 0.1 | 830.5 | −0.00 | −0.002 | 0.105 | −0.033 | 0.0199 | −0.0040 |
| B: total ability, ψa 1 | 3949.2 | 0.50 | 0.071 | 0.493 | −0.547 | 0.0915 | 0.0047 |
| D: total adult h, ψa 1 | 5671.0 | 0.62 | 0.089 | 0.580 | −0.675 | 0.1089 | 0.0068 |
| HS | 364.1 | 0.10 | 0.012 | 0.069 | 0.023 | 0.0117 | 0.0045 |
| HS + A, ψa 32 | 335.6 | 0.15 | 0.019 | 0.115 | −0.027 | 0.0199 | 0.0056 |
| HS + B, ψa 0.1 | 336.2 | 0.16 | 0.020 | 0.119 | −0.031 | 0.0206 | 0.0056 |

- **Mean forms (A, C) amplify existing sorting; they do not create it.** At the baseline, students are very slightly abler in the city, so stronger peers move every sorting outcome away from the data. Under HS the sign flips and S13 falls from 364 to 336, but only at $\psi_a=32$, hundreds of times the literature's 0.03–0.16, and children's ability sorting barely moves. C behaves like A and feeds on itself, because adults carry the school quality of the place where they grew up; no equilibrium converges beyond $\psi_a\approx2$.
- **Total forms (B, D) create a quality gap from enrollment.** They add $\psi_a$ times the targeted enrollment gap (0.376) to the suburb's quality advantage. Income sorting rises with $\psi_a$, to half the data's index or more at $\psi_a=1$, but the score, growth and staffing gaps are then 6–50 SEs off.
- **No specification sorts households by ability.** Children's ability gap is at most 0.008 SD in every refit, and the own-AFQT city-share gradient stays 6 SEs from the data. Given the BFM WTP, parents value quality too little for quality gaps of the size SEDA growth allows to sort households.

The baseline is unchanged. Details: [calibration log](calibration_log.md#peer-effect-strength-and-specification-2026-10-25) and [scores of every refit](runs/peer-20261025/refit_scores.md).

### Constant class size with random teacher assignment

The baseline assumes efficient class sizes within a location (Proposition 1), so every student there receives the same school quality. The constant-class-size (CCS) model ([theory note](../../../notes/spatial_teachers_model_ccs.md)) instead gives every class in a location the same size and assigns each student at random to one teacher for all of school. Her school input $q=h_T^\beta N_l^{-\sigma}$ becomes part of her state, and parents value a location's whole teacher distribution. The draw rescales the young's occupational shocks, so the model needs no new state dimension, and it nests the baseline when a location's teachers are identical. An extension adds noisy teacher effectiveness θ, whose log has correlation ρ with $\log h_T$ as in the CFR measurement, and $K_d$ independent teacher draws; $K_d=12$ is a new teacher each year.

The switch is `--external ccs=1`, off by default and then bit-for-bit inert. `ccs_draws` sets $K_d$ and `ccs_rho` sets ρ, with the matching CFR concept `--measure cfr_rho=ρ`. Every row below refits the eight baseline parameters exactly to the eight targets, so every outcome shown is untargeted; the mother–child correlation and the wage slope also check the first-stage estimates of $\rho_z$ and $s_z$.

| Fit | β | γ | S13 | Sorting / allocation / pay / first stage (Σz²) | Mother–child corr. (z) | Wage slope (z) | Share of score variance from the draw | Earnings gap |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Data |  |  | — |  | 0.577 | 0.213 |  | 0.112 |
| Baseline (efficient class sizes) | 0.228 | 0.619 | 864.1 | 379.2 / 23.2 / 446.0 / 15.7 | 0.576 (−0.0) | 0.214 (+0.1) | — | −0.010 |
| Efficient, ρ 0.3 | 0.068 | 0.619 | 866.2 | 380.5 / 24.0 / 446.0 / 15.7 | 0.576 (−0.0) | 0.214 (+0.1) | — | −0.010 |
| **Simple CCS** (one teacher, ρ 1) | 0.211 | 0.569 | 1472.8 | 393.0 / 15.8 / 619.5 / 444.6 | 0.298 (−19.9) | 0.258 (+5.9) | 0.672 | −0.012 |
| CCS + θ (ρ 0.3) | 0.063 | 0.568 | 1479.6 | 395.5 / 16.8 / 622.7 / 444.7 | 0.298 (−19.9) | 0.258 (+5.9) | 0.672 | −0.013 |
| CCS + K (a new teacher each year) | 0.226 | 0.614 | 896.2 | 380.4 / 22.2 / 460.6 / 33.1 | 0.518 (−4.2) | 0.216 (+0.4) | 0.146 | −0.010 |
| CCS + θ + K | 0.068 | 0.614 | 899.3 | 382.4 / 23.2 / 460.6 / 33.1 | 0.518 (−4.2) | 0.216 (+0.4) | 0.146 | −0.011 |

Rows use the linear warm glow ($\psi=0$). Block scores and S13 are sums of squared standardized errors; z is the standardized error.

- **Random assignment adds luck, not sorting.** The two locations' teacher draws are almost equally dispersed (log-quality SD 0.151 vs 0.149 in the simple model), so dispersion adds about 3% as much to the suburb's value as mean quality does. Sorting outcomes move by about 0.002 at most. Teachers select on worse draws, so the teacher premia worsen.
- **With one teacher, the draw dominates the first stage.** The CFR target fixes the draw's dispersion, so the draw accounts for two thirds of the adult score variance with one teacher and 15% with yearly draws. The draw is not inherited: the mother–child correlation falls to 0.298 and the wage slope misses by 5.9 SEs. Re-estimating $\rho_z$ cannot close the gap, since even $\rho_z\to1$ leaves the correlation at 0.517 ([ρz diagnostic](runs/ccs-simple-20261101/rhoz_table.md)).
- **Noisy effectiveness is almost a relabelling.** At ρ = 0.3, β falls to about ρβ, while the draw's dispersion, the first stage and sorting barely move. Only $K_d$ matters: yearly draws raise the correlation to 0.518 (−4.2 SEs) and bring S13 back to 896.
- **Warm-glow curvature changes who values quality, not where they live.** Log and CRRA-2 kernels ($\psi=1,2$) change only λ and move S13 by less than one point. As the theory predicts, the value of teacher dispersion changes sign at log utility, but it is negligible here.

**CCS combines with every eligible mechanism.** The November 2 round crossed CCS with the six factorial mechanisms (C commuting, H housing cost, F fertility, S salary schedule, L licensure, P peers) in 161 exact refits: all 64 combinations under the simple model and under $K_d=12$, and 11 combinations under each ρ = 0.3 variant. The CCS cost is S13 minus that of the efficient fit with the same mechanisms.

| Model | None: S13 (CCS cost) | Best fit: S13 (CCS cost) | Highest income sorting (SI) | Mother–child corr. |
|---|---|---|---|---|
| Efficient (factorial) | 864.1 | HS/HSP 364.0 | CHSP 0.320 | 0.567–0.576 |
| Simple CCS | 1473.9 (+609.8) | HS/HSP 827.1 (+463.1) | CHSP 0.303 | 0.295–0.300 |
| CCS + K | 896.3 (+32.2) | HS/HSP 383.2 (+19.2) | CHSP 0.322 | 0.512–0.518 |
| Efficient, ρ 0.3 | 866.2 | HS/HSP 392.5 | CHS 0.261 | 0.569–0.577 |
| CCS + θ | 1480.6 (+614.4) | HS/HSP 854.9 (+462.4) | CHS 0.243 | 0.296–0.299 |
| CCS + θ + K | 899.4 (+33.2) | HS/HSP 413.5 (+21.0) | CHS 0.260 | 0.513–0.518 |

At ρ 0.3 the CCS cost is against the efficient ρ 0.3 fit. The correlation is the range across each model's fits (data 0.577).

With yearly draws, CCS is nearly neutral and additive: it costs 18–34 S13 points in every combination, and each mechanism's effect is the factorial's to within 7 points. With one teacher, the first-stage cost is 406–429 points in every combination, and no mechanism repairs it. CCS amplifies no mechanism's income sorting (interaction effect on SI at most 0.006).

**Conclusion.** The baseline remains the September efficient fit. One-teacher CCS cannot be reconciled with the NLSY first stage under the CFR/χ bridge, with or without the eligible mechanisms. $K_d=12$ is the variant consistent with it, at 32 S13 points over the baseline, and it preserves the factorial's ranking (best HS/HSP, 383 against 364). No CCS version helps with sorting; which, if any, to carry forward is open.

Sources: comparison tables for [November 1](runs/ccs-simple-20261101/compare_table.md) and [October 26](runs/ccs-20261026/refit_table.md), the November 2 [results](runs/ccs-mech-20261102/RESULTS.md) and [scores](runs/ccs-mech-20261102/scores.md), and the calibration log ([October 26](calibration_log.md#constant-class-size-and-random-teacher-assignment-2026-10-26), [November 1](calibration_log.md#simple-constant-class-size-baseline-and-the-θk_d-comparison-2026-11-01), [November 2](calibration_log.md#ccs-with-the-six-eligible-mechanisms-2026-11-02)).
