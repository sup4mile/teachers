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
| Peer effects | At the tested magnitude, effects are negligible because the model generates little ability sorting. Peer effects do not initiate residential sorting. |

**Strongest result.** Housing costs plus the salary schedule are the simplest combination near the best untargeted fit: the sum of squared standardized errors across 13 validation moments falls from 864 to 364. Most of the improvement comes from teacher pay. Income and ability sorting, teacher selection and the gender wage gap remain unresolved. At the main externally fixed magnitudes, no combination achieves even one-third of the observed income-sorting index; larger housing charges strengthen sorting but make the interpretation and empirical support for the charge consequential. [Comparison, sensitivity checks and evidence](calibration_log.md#factorial-round-the-six-eligible-mechanisms-in-every-combination-2026-10-19).
