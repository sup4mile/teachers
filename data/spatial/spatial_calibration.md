# Calibration of the spatial teachers model

This note calibrates [spatial_continuous.jl](../../julia/spatial_model/spatial_continuous.jl) ([writeup](../../notes/spatial_teachers_model_eps0.md)) to a within-commuting-zone economy: location 1 is the city, location 2 the suburb. Inputs are 2018 values plus SEDA's pooled 2009–19 outcomes, read as a late-2010s cross-section. The exception is the occupational block (Table 2) and the teaching margin, which use ACS 2009–13 (caveat 3 under Table 2). Spatial gaps are enrollment-weighted within-CZ contrasts, suburb minus city, in logs unless noted; revenue and SES partitions are robustness checks.

**No parameter has been estimated yet.** Values below are targets, literature inputs or starting values.

| Block | Parameters | Disciplined by | Solves equilibrium? |
|---|---|---|---|
| Fixed inputs (Table 1) | Normalizations; $\eta$; $\phi$; $\sigma$; $\sigma_\nu/\mu$; $\psi$; $\chi$; $K$ | Aggregates, schooling, literature | No |
| Occupations and ability (Table 2) | Relative $A_i$; $\sigma_\epsilon$; $\rho_z$, $s_z$; relative female wedges | ACS 2009–13; NLSY79/CNLSY (first stage, $Q_l$ fixed) | No; checked after the internal fit |
| Internal (Table 3) | $\bar\kappa$; female outside option; $\gamma$; $\delta_\kappa$, $\Delta B$, $\beta$, $\lambda$, $r_m$ | Teaching shares, district moments, CFR, BFM | Yes |

Code: `spatial_calibrate.jl` (inputs, targets, diagnostics), `spatial_moments.jl` (model moments), `spatial_estimate.jl` (TikTak search).

## 1. External inputs

### Table 1: Fixed inputs

Asterisks mark provisional values not yet rebuilt from the 2018 sample.

| Parameter | Value | Basis | Source |
|---|---|---|---|
| Ability elasticity $\alpha$ | 1 | Normalization; only $\alpha\sigma_\epsilon$, $\alpha s_z$ identified | — |
| Population $M$, amenity $B_1$, home-production productivity $A_{HP}$ | 2, 0, 1 | Normalizations | — |
| Consumption weight $\mu$ | 1 | Normalization (below) | — |
| Goods-investment elasticity $\eta$ | 0.080 | $(S-t)/(1-t)$, $S=0.102$, $t\approx0.023$ | [NCES 605.20](https://nces.ed.gov/programs/digest/d21/tables/dt21_605.20.asp); [PWT 11](https://fred.stlouisfed.org/series/LABSHPUSA156NRUG); [NCES 211.60](https://nces.ed.gov/programs/digest/d18/tables/dt18_211.60.asp), [236.20](https://nces.ed.gov/programs/digest/d20/tables/dt20_236.20.asp) |
| Time-investment elasticity $\phi$ | 1.12* | Non-teacher time share $s_O=0.55$ | ACS; HHJK (2019) mapping |
| Class-size curvature $\sigma$ | 0.56* | Class-size wage effect | [FOO (2013)](https://pure.uva.nl/ws/files/1912854/124460_394782.pdf#page=26) |
| Location taste scale $\sigma_\nu/\mu$ | 0.25* | Literature | [Eckert–Kleineberg (2024)](https://www.fpeckert.me/papers/ek_2024.pdf#page=17) |
| Warm-glow curvature $\psi$ | 0 | Linear warm glow | — |
| Earnings–achievement bridge $\chi$ | 0.113 | $\log 1.12$; converts school-driven score changes only, not $s_z$ (section 3) | [CFR (2014b)](https://opportunityinsights.org/wp-content/uploads/2018/03/teachers2.pdf) |
| School years of exposure $K$ | 12 | Grades 1–12, equal per-year exposure | Model timing |

**Goods investment.** Spending over labor income is $S=\eta(1-t)+t$, where $t$ is the tax funding teacher payroll, so $\eta=(S-t)/(1-t)$. For 2018, $S=0.060/0.591=0.102$ (education spending over the PWT labor share), and $t\approx0.023$ (3.2 million FTE teachers at \$60,483, benefits at 0.43 of salary, against 2018 GDP of \$20.66 trillion); unrounded, $\eta=0.0807$ (`eta_2018`). The old 0.103 is the gross $S$, which counts teacher pay twice. Counting all instructional staff gives $t=0.030$, $\eta=0.073$. At the reference solution, teacher payroll is 4.2% of market income, against $t=0.023$: the model's 25–34-year-olds teach at about twice the economy-wide rate. The solver's own $t$ (2.9%) also taxes home-production output, which is 31% of model income.

**Time investment.** With zero goods moving costs, $s_O=\mu\phi/(\mu\phi+1-\eta)$ and $s_T=\mu\phi\gamma/(\mu\phi\gamma+1-\gamma\eta)$, which depend only on $\mu\phi$. Location choice and WTP depend only on ratios to $\mu$, and $\bar\kappa$ absorbs its remaining role, so $\mu=1$ and $\phi$ matches non-teachers' schooling ($s=$ years/25). With $r_m>0$, evaluate at the reference $\Xi$.

**Known miss.** With a common $\phi$ and $\gamma<1$, teachers get 0.91 of non-teachers' schooling ($s_T=0.50$), against 1.23 in the data. Report this as untargeted. The implied teacher leisure advantage (≈0.11 log consumption) is absorbed by $\bar\kappa$.

**Exposure.** Per-year evidence pins $\sigma/K$ and $\beta/K$, while $\beta/\sigma$ is free of $K$, so $K$ sets how much schools matter. $K=5$ (SEDA's tested grades) is a sensitivity case.

**Location tastes.** Without housing, effective dispersion is roughly $\sigma_\nu/\mu$ plus housing share over supply elasticity ($\approx0.45$), so 0.25 likely overstates population responses; sensitivity extends to 0.50.

**Units.** Freeze benchmark $\bar C$ and $\bar h$, which scale the moving cost and warm glow. Keep consumption, moving costs and WTP in the same units. They are frozen at $\vartheta_0$ in `spatial_reference.toml` ($\bar C=0.012554$, $\bar h=0.021000$) and depend on Table 2, so refreeze after T2.

### Table 2: Occupations and ability, estimated outside the equilibrium

Targets come from the ACS 2009–13 sample in `data/LaborMarketData/wages_occ_shares_v2.xlsx`: ages 25–34 with at least one year of high school, military and unemployed dropped. Home production is out of the labor force or under 15 hours a week; 15–29 hours are split equally between home production and the reported occupation. Wages are hourly, for full-year, full-time workers. The model has 21 occupations: the file's 19 market occupations, home production and K–12 teaching (postsecondary and other teachers are market occupations).

| Parameter | Target | Value | Source |
|---|---|---|---|
| Relative productivities $A_i/A_{HP}$, 19 market occupations | Men's shares over the 20 non-teaching occupations | `moments_shares` | ACS 2009–13 |
| Relative female wedges $\tau^\omega_{i,f}$, 19 free | Women's shares over the 20 non-teaching occupations | `moments_shares` | ACS 2009–13 |
| Ability dispersion $\sigma_\epsilon$ | Weighted within-occupation 90/10 of hourly wages, non-teachers | 3.74 | ACS 2009–13 |
| Ability persistence $\rho_z$ | Mother–child correlation of latent test scores (net of the $Q_l$ channel) | 0.577 (SE 0.014); 0.53–0.63 across samples | NLSY79 AFQT, CNLSY PIAT/PPVT |
| Stationary dispersion $s_z$ | Slope of log wage on standardized latent AFQT, $\approx s_z/(1-\eta)$ | 0.196 (SE 0.007) | NLSY79; NLSY97 0.175 |
| *Check only* | Occupation-specific 90/10s; mean wages by occupation | `90_10_hr_wages_weighted` | ACS 2009–13 |
| *Check only* | Mother–child and sibling correlations: earnings 0.135 and 0.281; latent scores 0.606 (siblings) | `nlsy_ability.md` | NLSY79, CNLSY |

Forty-one parameters, forty-one moments. No education wedges and no wage wedges for men or teaching. The levels of $A_i$ and of the female outside option relative to teaching are internal (Table 3).

**Ability block.** Estimate $(\sigma_\epsilon,\rho_z,s_z)$ in a first stage outside equilibrium, holding $Q_l$ and $t_l$ fixed; within a cell they are constants. Non-teachers' log wages are then $b[\alpha\log z+\log(\Theta_{i,g}\epsilon_i^\alpha)]$ plus a constant, and $z$ is independent of occupation, gender and $\epsilon$. Each moment then pins one parameter. The model's latent mother–child score correlation is $\rho_z$, so $\rho_z=0.577$. The slope per latent SD is $b\,\alpha s_z$, so $s_z=0.196(1-\eta)/\alpha=0.181$ (`implied_s_z`). $\sigma_\epsilon$ solves the pooled 90/10 given $s_z$, with the shares block inverted at each value; the 90/10 depends on both $z$ and $\epsilon$, while scores measure $z$ but not $\epsilon_i$, which separates them. Both data targets are latent (the factor-model correlation and the slope over $\sqrt{R}$), so no simulated panels or score noise are needed; $\mathrm{var}(u)$, from the composite reliabilities (0.92 for mothers, 0.89 for children), matters only for composite-level checks. The model's latent sibling correlation is $\rho_z^2\approx0.33$, against 0.606 in the data (check only). Use individual earnings (no spouses in the model) and read $\rho_z$ as reduced form, including assortative mating. What the first stage leaves out goes to the section 4 consistency check: sorting on $Q_l$ in scores (signal $\log z+c\log Q_l$) raises the model correlation above $\rho_z$; cross-location variation in $Q_l$ and $t_l$ adds wage dispersion, lowering $\sigma_\epsilon$; and teachers (about 4% of the slope sample) load on $z$ with a different elasticity and select on $z$ in GE.

**Ability grid.** Rouwenhorst matches the stationary variance and autocorrelation of $\log z$ exactly at any Nz, so $\rho_z$ and $s_z$ do not depend on it; the 90/10, a quantile, does. Compute the first stage on the full model's Rouwenhorst grid, not a continuous normal, so it hands the GE the distribution the GE uses. Compare Nz = 5, 9 and 15 with the continuous-$z$ limit and set the full model's Nz (`SOLVER` in `spatial_calibrate.jl`, now 5) where the 90/10 and $\hat\sigma_\epsilon$ stop moving; confirm it in the GE (section 4), since solve time grows with Nz.

**Scale and $\chi$ (decided).** $s_z$ is free and targets the NLSY slope, which matches the model object: latent scores, log hourly wages at 25–34, the 90/10's full-year full-time sample. CFR's 12% is a descriptive association, and each of its definitions (single-year observed scores in grades 3–8, earnings at 28 with zeros as a share of the mean, controls) pushes it below $b\,s_z$. $\chi$ only converts school-driven score changes (SEDA, WTP): a score SD from $Q_l$ is worth $b\,s_z/c$, so $c=b\,s_z/\chi\approx1.7$, also used in the score counterparts. Whether CFR's cross-sectional 12% supports this reading is open (T6a); otherwise report the $\chi$ inconsistency. The alternative $b\,s_z=\chi$ ($s_z=0.104$, about 12 SEs below the slope) is in Table 5. It cuts $z$'s share of within-occupation log-wage variance from about 15% to 5%, and sorting and teacher selection scale with it.

**NLSY moments (T2a).** From `nlsy_ability.py`; SEs from a 200-draw bootstrap over NLSY79 households.

- *Persistence:* the latent mother–child correlation is 0.577. Sample choices move it more than sampling error (0.53 with cross-sectional mothers, 0.63 unweighted), so both variants go to Table 5. Sorting on $Q_l$ raises the score correlation, so expect $\rho_z$ at or somewhat below 0.58. The sibling ratio $\mathrm{cov}_{sib}/\mathrm{cov}_{pc}=1.05$ is a check only (siblings share inputs outside the model).
- *Scale:* the log-wage slope is 0.196 per latent SD (men 0.179, women 0.224), so $s_z\approx0.18$; NLSY97, the ACS cohorts, gives 0.175, so vintage does not close the gap to $\chi$.
- *Selection into work:* each latent SD raises full-year full-time work by 0.08 and lowers working under 15 hours by 0.06 (0.08 for women). The model makes participation independent of $z$, so mothers' earnings are a selected measure.

**Pooled 90/10.** The target averages cell ratios with occupation-by-gender wage-sample weights (`occ_gender_weights`; home production has zero weight); build the model counterpart the same way. Age-residualized, it is 3.64.

**Why no equilibrium.** Non-teachers choose the occupation maximizing $(1-\tau^\omega_{i,g})A_i\epsilon_i^\alpha$, so $z$, $Q_l$, taxes, teacher pay, amenities, altruism and moving costs drop out of the shares. Selection into teaching biases $\log A_i$ by at most 0.013 at a 2% teaching share (0.039 at 10%): ignore it for men (1.9%), correct it for women (6.0%). Scores still carry $Q_l$; the external-block consistency check (section 4) covers it.

**Caveats.**

1. If $z$ alone gives a 90/10 above 3.74, $\sigma_\epsilon$ is at a corner. Fix $\eta$ first, since wages load on ability through $\alpha/(1-\eta)$.
2. Measurement error and transitory shocks in ACS wages bias $\sigma_\epsilon$ up given $s_z$.
3. Table 2 and the teaching margin use ACS 2009–13, against late-2010s spatial targets. Wage dispersion was rising (90/10 of 3.36 in 2000), and men's home-production share (22.9%, against 18.6% in 2000) reflects the recession, lowering market $A_i$ and the levels set against them, including $\bar\kappa$. Within-CZ gaps should be less sensitive (T7).
4. Mean wages and cell-level 90/10s are untargeted checks on a common $\sigma_\epsilon$.

## 2. Internal calibration

Write $\kappa_1=\bar\kappa$, $\kappa_2=\bar\kappa e^{\delta_\kappa}$, $B=(0,\Delta B)$ and $m_{12}=m_{21}=r_m\bar C$. Fit all parameters jointly; each row names the main source of identification. District moments are from [spatial_moments.md](estimates/spatial_moments.md).

### Table 3: Internal parameters and targets

| Parameter | Target | Value | Source |
|---|---|---|---|
| ***Teaching margin*** | | | |
| Teacher pay level $\bar\kappa$ | Male teaching share | 1.9% | ACS 2009–13 |
| Female outside-option level | Female teaching share | 6.0% | ACS 2009–13 |
| Teacher wage curvature $\gamma$ | Teachers' 90/10 of hourly wages | 2.77; start $\gamma$ at 0.83 | ACS 2009–13 |
| ***Spatial block*** | | | |
| Spatial pay gap $\delta_\kappa$ | Adjusted salary log gap | 0.01623 | District, 176 CZs |
| Amenity gap $\Delta B$ | Enrollment log gap | 0.3764 | District, 188 CZs |
| Teacher spillover $\beta$ | Earnings at 28 per teacher-quality SD-year | ≈1.3% | CFR (2014b) |
| Altruism $\lambda$ | WTP per SD of school-average scores | \$19.70/month (SE \$7.40) | BFM (2007), Table 7 |
| Goods moving cost $r_m=m/\bar C$ | Household-income log gap | 0.2708 | District, 187 CZs |
| — (overidentifying) | Teacher FTE per pupil log gap | −0.01216 | District, 188 CZs |

Eight parameters, nine moments.

**Target definitions.**

- **Teaching shares:** K–12 teachers over all 25–34-year-olds of each gender, home production included in the denominator (Table 2 sample).
- **Teacher 90/10:** wage-sample-weighted average of men's (2.60) and women's (2.83) ratios (24% men, 76% women); build the model counterpart the same way.
- **Salary gap:** CWIFT-adjusted; matches $\Delta\log[\kappa_l\,\mathbb E(h_T^\gamma\mid T,l)]$, which includes teacher composition. Prefer standardized schedules where available.
- **Enrollment:** equals the population-share gap with equal cohorts.
- **WTP:** BFM report 1990 housing user-cost dollars; convert to a consumption share (see section 3) and reweight to parents where feasible.
- **Income:** district median household income is only a proxy for model labor income. Build a comparable statistic or move it to validation, and check whether it identifies $r_m$ with the Jacobian and profiles.
- **FTE per pupil:** matches teacher headcount per student $n_l/S_l$.

### Table 4: Untargeted validation

| Outcome | Gap (data) | Model counterpart |
|---|---|---|
| Child poverty | −8.88 pp | Low-income child share (mapping to be defined) |
| SEDA learning rate | 0.0031 SD/grade (SE 0.0030) | $\Delta\log Q/(\chi K(1-\eta))$ |
| SEDA score level | 0.2887 SD | Composition plus $Q$; descriptive |
| Local revenue share | 0.4726 | Diagnostic only |

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

converting BFM's school-average SD (74 points) to student-level SDs before applying $\chi$. Use the \$19.70 structural estimate, not the \$17.30 hedonic one.

## 4. Estimation and sensitivity

Minimize $[\mathcal M(\vartheta)-\widehat{\mathcal M}]'W[\mathcal M(\vartheta)-\widehat{\mathcal M}]$ with diagonal inverse-variance weights, adding transport uncertainty to literature moments. Bootstrap whole CZs, re-estimating Table 2 in the same loop. Impose $0<\beta<1$, $\lambda\ge0$, $r_m\ge0$.

- **Identification:** Jacobian rank and singular values; whether $\lambda$, $\Delta B$ and $\beta$ substitute for each other; a profile of $r_m$ including $r_m=0$.
- **Solution:** equilibrium residuals, positive consumption, quadrature accuracy, Nz convergence (5, 9, 15; Table 2, *Ability grid*), multiple starts.
- **External-block consistency:** recompute Table 2 moments in the full model at the optimum: the mother–child score correlation with $c\log Q_l$ (parents' birth locations from the migration kernel, children's from mothers' location choices), wages including $Q_l$ and $t_l$, and teachers in the slope sample. If they differ beyond sampling error, re-estimate Table 2 net of the model-implied distortion and refit.
- **Residential transitions:** a childhood-to-adulthood transition rate, conditioned on pre-move resources, would help identify $r_m$.

### Table 5: Sensitivity panel

| Input | Values | Notes |
|---|---|---|
| $r_m$ | 0, 0.05, 0.10, 0.1813, 0.30 | Refit the rest at each value |
| $\psi$ | 0, 0.5, 1 | |
| $\sigma_\nu/\mu$ | 1/7, 0.20, 0.25, 0.30, 0.40, 0.50 | Upper end allows for housing |
| $\sigma$ ($K=12$) | 0.15, 0.30, 0.60, 0.85 | Refit $\beta$ |
| $\chi$ | 0.10–0.20 | Rescales $c$ at the baseline $s_z$ |
| $s_z$ | Free (baseline); $b\,s_z=\chi$ ($s_z=0.104$) | Refit $\sigma_\epsilon$ and the shares block |
| $\rho_z$ target | 0.577; 0.53 (cross-sectional mothers); 0.63 (unweighted) | Sampling choice moves it more than the SE |
| $K$ | 5, 9, 12 | Rescale $\sigma$, $\beta$ |
| $\eta$ | 0.073, 0.080, 0.103 | Recompute $\sigma$ |
| Table 2 parameters | Bootstrap draws | |

**Interpretation limits.** The model bundles residence, work and school, and local taxes fund teacher payroll only. It omits effort responses, teacher–student match effects and childhood exposure timing. The test is whether it jointly explains pay, sorting and a small school-effectiveness gap; a large teacher contribution to the score-level gap should be tested, not imposed.

## 5. To-do

- [x] **T1 — Model configuration:** 21 occupations (19 market, home production, K–12 teaching), passing occupation-specific `A` and `τω` to `Params`; rebuild $\eta$ for 2018; fix $\bar C$, $\bar h$. `external_params()` builds the 21-occupation model from provisional Table 2 values; $\vartheta_0$ matches both teaching shares ($\bar\kappa=0.2365$, female level 0.934) and the CFR map ($\beta=0.387$). The female level and $\gamma$ are not yet in $\vartheta$.
- [ ] **T2 — Table 2:** ACS 2009–13 shares and 90/10s are in `wages_occ_shares_v2.xlsx`. Remaining: schooling for teachers and non-teachers (sets $\phi$; independent of T2b); T2b.
  - [x] **T2a — NLSY data:** `nlsy_ability.py` links 8,137 CNLSY children to 3,480 NLSY79 mothers, age-norms the tests, fits the two-factor measurement system, estimates the log-wage slope on an ACS-matched sample (with an NLSY97 vintage check), computes the sibling and earnings checks and bootstraps by household. Main weights are the mother's 1979 weight split across her children; sensitivities cover the CNLSY weights and a 1981–2000 birth-cohort restriction. Results are under Table 2.
  - [ ] **T2b — First-stage estimation:** outside equilibrium, with $Q_l$ and $t_l$ fixed (Table 2, *Ability block*). (i) Set $\rho_z=0.577$ and $s_z=0.181$ in closed form. (ii) Fit $\sigma_\epsilon$ to the pooled 90/10 by a one-dimensional search: at each value, invert the shares block ($A_i$, $\tau^\omega_{i,f}$, with the women's teaching-selection correction); in each occupation–gender cell, convolve $b\,\alpha\log z$ on the Rouwenhorst grid with the Roy-selected $b\log(\Theta_{i,g}\epsilon_i^\alpha)$; average cell 90/10s with `occ_gender_weights`. (iii) Nz: run 5, 9 and 15 against the continuous-$z$ limit and pick the full model's Nz (*Ability grid*). (iv) Set $c=b\,s_z/\chi\approx1.7$ and report the $\chi$ gap and the sibling check. (v) Replace `σϵ`, `ρz`, `s_z` and the `occupation_block` placeholders in `spatial_calibrate.jl`; refreeze $\bar C$, $\bar h$ (`--freeze-reference`). (vi) Table 5 variants: $b\,s_z=\chi$ and the two $\rho_z$ targets change only the closed-form rows, then refit $\sigma_\epsilon$. No simulated panels or score noise. The $Q_l$, $t_l$ and teacher pieces go to the section 4 consistency check.
- [ ] **T3 — Spatial targets:** income statistic; audit pay and FTE counts; convert CFR and BFM targets; residential-transition target.
- [ ] **T4 — Model counterparts** for every target in Tables 2–4.
- [ ] **T5 — Freeze targets:** definition, sample, weights and uncertainty for each.
- [ ] **T6 — Verify and run:** Section 4 checks, joint fit, consistency check, sensitivity panel.
  - [ ] **T6a — CFR replication in CNLSY:** on CNLSY children past 28, regress earnings at 28 (zeros included, as a share of mean) on observed, not latent, scores at about ages 8–14. Near 0.11 means the gap is definitional and supports the baseline. If the CNLSY gives about 0.11 even with ACS-style definitions, the NLSY slope is the outlier, and $b\,s_z=\chi$ becomes the stronger case.
- [ ] **T7 — 2018 occupational block (if time):** rerun the ACS extract on a sample centered on 2018 (e.g., pooled 2016–19) with the same definitions; refit Table 2 and the teaching margin; compare $A_i$, wedges and $\sigma_\epsilon$ with 2009–13.

**Later extensions:** teacher geocodes; residence/workplace split; payroll-specific transfers; pre-1990 cross-sections; Opportunity Atlas crosswalk.
