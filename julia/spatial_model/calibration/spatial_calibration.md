# Calibration of the spatial teachers model

Calibrate the [spatial model](../spatial_continuous.jl) ([writeup](../../../notes/spatial_teachers_model_eps0.md)) to a city and suburb within a commuting zone (CZ). Spatial gaps are enrollment-weighted within-CZ contrasts, suburb minus city, in logs unless stated otherwise. District inputs are from 2018 and SEDA outcomes pool 2009–19; occupations and teaching shares use ACS 2009–13.

**Strategy.** Fix aggregate and literature inputs (Table 1); estimate occupations and ability outside equilibrium (Table 2); then solve jointly for eight internal parameters that match eight moments exactly in equilibrium (Table 3). Check the first stage in the fitted equilibrium and assess untargeted outcomes (Table 4). The [calibration log](calibration_log.md) contains definitions, derivations, code references, run history and known issues.

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

Teaching-share denominators include home production; the teacher 90/10 averages gender-specific ratios with wage-sample weights. District targets cover 188 CZs. CFR maps to $\beta\,\mathrm{sd}(\log h_T)/[K(1-\eta)]$; BFM maps to compensating consumption, with one school-average SD taken as 0.45 student SDs. Teaching shares use 10% relative fit scales and the teacher 90/10 uses 5%; other scales are the reported SEs except the move rate (0.05). With exact identification the scales guide the search but do not affect the estimate. [Definitions and sources](calibration_log.md#2-internal-calibration).

The exactly identified fit (2026-09-29) matches all eight targets (criterion $1.6\times10^{-19}$). Levenberg–Marquardt from 15 starts across the search box finds no other root, and finer grids keep every target matched ([details](calibration_log.md#exactly-identified-fit-2026-09-29-the-baseline)). The K–12 FTE-per-pupil gap is validation (Table 4). With the salary gap matched, the model puts more teachers per pupil in the suburb, the opposite sign from the data. The earnings gap does not identify moving costs and is also validation.

**Moment Jacobian at the chosen parameters.** The table reports $J^s_{kj}=(w_j/s_k)\,\partial m_k/\partial\vartheta_j$ at the unrounded Table 3 estimate, where $s_k$ is the fit scale defined above and $w_j$ is the parameter's search-box width. In column order, $w=(\log 10,\,0.7,\,0.6,\,1,\,1,\,0.6,\,1.5,\,1)$. Derivatives use central differences with steps $0.001w_j$, re-solving equilibrium with Tables 1–2 and the reference units fixed ([saved diagnostics](runs/exact-polish-base/report.txt)). Thus a local increase of 1% of a parameter's box width changes moment $k$ by approximately $0.01J^s_{kj}$ fit-scale units; the unscaled derivative is $s_kJ^s_{kj}/w_j$.

| Moment | $\log\tilde\kappa$ | $\omega_f$ | $\gamma$ | $\delta_\kappa$ | $\Delta B$ | $\beta$ | $\lambda$ | $r_m$ |
|---|---|---|---|---|---|---|---|---|
| Male teaching share | 46.0 | 9.65 | −9.55 | 12.1 | 0.375 | 20.2 | 0.00209 | 1.91 |
| Female teaching share | 35.7 | −15.6 | −14.4 | 9.33 | 0.210 | 15.8 | 0.00119 | 1.70 |
| Teachers' hourly-wage 90/10 | 2.94 | −0.690 | 11.1 | 1.11 | −0.501 | 0.987 | −0.00264 | −0.453 |
| Adjusted salary gap | −2.13 | −0.112 | −0.982 | 209 | 0.118 | 0.263 | 0.0560 | 1.87 |
| Enrollment gap | −1.80 | 0.212 | 0.284 | −6.49 | 49.6 | 0.672 | 0.278 | 2.85 |
| CFR earnings effect | 0.00570 | −0.382 | −1.24 | −0.0355 | −0.0601 | 8.86 | −0.000359 | 0.413 |
| BFM WTP share | 7.52 | −1.56 | −1.23 | 1.95 | −0.0276 | −3.45 | 32.0 | 0.215 |
| City↔suburb move rate | 11.3 | −0.884 | −1.91 | 2.63 | −3.14 | −5.36 | −0.0180 | −19.4 |

The standardized Jacobian has rank eight, singular values from 6.94 to 210 and condition number 30.3, supporting local identification. Salary, enrollment, WTP and mobility respond strongly to $\delta_\kappa$, $\Delta B$, $\lambda$ and $r_m$, respectively. Pay scale and $\beta$ overlap through teaching shares (each column's projection on the other seven has $R^2\approx0.90$), but the full set of moments separates them. These magnitudes and the condition number depend on the stated scaling.

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

The first-stage consistency and numerical grid checks pass. The main failures are income sorting, the staffing gradient, teacher schooling and relative pay, occupational wage levels/dispersion, and the gender wage gap. Matching the salary gap puts the FTE gap 3.1 SEs above the data and the SEDA growth gap 2.9 SEs above it. Teacher/non-teacher schooling is 0.76 versus 1.18 in the data; the market gender log-wage gap is 0.30 versus 0.10. See [known issues](calibration_log.md#known-issues-at-the-current-fit).

**Sensitivity.** The deterministic panel and ACS 2016–19 refit remain to be run at the exactly identified baseline. [Table 5 details](calibration_log.md#table-5-sensitivity-panel).

**Remaining work.** Rerun the Table 5 panel at the exactly identified baseline; bootstrap the full calibration with Table 2 re-estimated, adding literature transport uncertainty; address the sorting, staffing and wage/schooling misses. [Detailed task log](calibration_log.md#5-to-do).
