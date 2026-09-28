# Calibration of the spatial teachers model

This note describes how to calibrate the spatial teachers model for a within-commuting-zone city/suburb economy: location 1 is the city and location 2 the suburb. The model is [spatial_continuous.jl](../../julia/spatial_model/spatial_continuous.jl), documented in [its writeup](../../notes/spatial_teachers_model_eps0.md). Inputs are 2018 values combined with SEDA's pooled 2009–19 outcomes, read as a late-2010s cross-section. Spatial gaps are enrollment-weighted averages of within-CZ contrasts, suburb minus city, in logs unless noted. Revenue and SES partitions are robustness checks.

**No parameter has been estimated yet.** Values below are targets, literature inputs, or starting values.

The parameters fall into three blocks, and only the last requires solving the equilibrium:

| Block | Parameters | Disciplined by | Solves equilibrium? |
|---|---|---|---|
| Fixed inputs (Table 1) | Normalizations; $\eta,\phi,\mu$; $\sigma$; $\sigma_\nu/\mu$; $\psi$; $\chi$; $K$ | Aggregates, literature, assumptions | No |
| Occupations and ability (Table 2) | Relative $A_i$; $\sigma_\epsilon$; $\rho_z$, $s_z$; relative female wedges | ACS shares and wages; NLSY parents, children, siblings | No; checked after the internal fit |
| Internal (Table 3) | $\bar\kappa$; female outside-option level; $\gamma$; $\delta_\kappa$, $\Delta B$, $\beta$, $\lambda$, $r_m$ | Teaching shares, district moments, CFR, BFM | Yes |

Code: [spatial_calibrate.jl](../../julia/spatial_model/spatial_calibrate.jl) holds inputs, targets and diagnostics; [spatial_moments.jl](../../julia/spatial_model/spatial_moments.jl) the model counterparts; [spatial_estimate.jl](../../julia/spatial_model/spatial_estimate.jl) the TikTak search. The [external-parameter audit](external_parameter_audit.md) derives $\sigma$, $\sigma_\nu/\mu$ and the moving-cost treatment.

## 1. External inputs

### Table 1: Fixed inputs

Asterisks mark provisional values carried over from the old draft.

| Parameter | Value | Basis | Source |
|---|---|---|---|
| Ability elasticity $\alpha$ | 1 | Normalization: only $\alpha\sigma_\epsilon$ and $\alpha s_z$ are identified. The code still uses 0.30 | — |
| Population $M$, amenity $B_1$, one $A_i$ | 2, 0, 1 | Normalizations | — |
| Goods-investment elasticity $\eta$ | 0.103* | Education spending/GDP (0.066) ÷ labor share (0.641) | U.S. aggregates, 1995/2000/2010 |
| Time-investment elasticity $\phi$ | 2.745* | Mincer schooling returns | Old calibration |
| Consumption weight $\mu$ | 0.714* | Schooling moments; see the feasibility note below | Old calibration |
| Class-size curvature $\sigma$ | 0.25* | Direct-wage mapping gives 0.2294 at $\eta=0.103$, $K=5$ (§3) | [FOO (2013), Table V](https://pure.uva.nl/ws/files/1912854/124460_394782.pdf#page=26); [audit §1](external_parameter_audit.md#1-class-size-025-has-a-defensible-magnitude-but-the-original-citation-was-incomplete) |
| Location taste scale $\sigma_\nu/\mu$ | 0.25* | Literature relative scale; set $\sigma_\nu=0.25\mu$ | [Eckert–Kleineberg (2024), §2.2, Table 1](https://www.fpeckert.me/papers/ek_2024.pdf#page=17); [audit §2](external_parameter_audit.md#2-location-tastes-correct-the-source-normalization-and-interpretation) |
| Warm-glow curvature $\psi$ | 0 | Specification: linear in child human capital | — |
| Earnings–achievement bridge $\chi$ | 0.13 | Measurement assumption (§3) | [CFR (2014b), §III.B](https://opportunityinsights.org/wp-content/uploads/2018/03/teachers2.pdf) |
| Effective school years $K$ | 5 | Grade intervals, grades 3–8 | SEDA grade span |

**Investment block.** $\eta$ is measured from aggregate education spending, so publicly financed teacher payroll must not also enter private goods investment. There is a schooling feasibility problem to resolve before fitting $\mu$ and $\phi$. With zero goods moving costs,

$$
s_O=\frac{\mu\phi}{\mu\phi+1-\eta},\qquad
s_T=\frac{\mu\phi\gamma}{\mu\phi\gamma+1-\gamma\eta},
$$

so $\gamma<1$ implies $s_T<s_O$. The old values imply $s_T/s_O\simeq0.933$, not the historical schooling ratio of 1.23. Goods moving costs make $s$ state-dependent, so this is a diagnostic rather than an impossibility result. Check the map from $s$ to years of schooling before keeping that target.

**Units and reference levels.** Compute benchmark $\bar C$ and $\bar h$ once and freeze them during fitting and counterfactuals: the goods moving cost and the warm glow are defined relative to them. Consumption, moving costs and WTP must share lifetime or annual units. $\chi$ and $K$ are measurement assumptions, not estimates.

### Table 2: Occupations and ability, estimated outside the equilibrium

| Parameter | Target | Data |
|---|---|---|
| Relative productivities $A_i/A_1$, $i\ne T$ | Men's occupation shares among non-teachers | ACS |
| Ability dispersion $\sigma_\epsilon$ (common to all occupations, including teaching) | Wage differences across occupations, given the shares (caveat 2) | ACS |
| Ability persistence $\rho_z$ | Sibling covariance ÷ parent–child covariance | NLSY79, Child/Young Adult |
| Stationary dispersion $s_z$, with $\sigma_\xi=s_z\sqrt{1-\rho_z^2}$ | Parent–child covariance of permanent log earnings, $b^2\rho_z s_z^2$ with $b=1/(1-\eta)$ | NLSY79, Child/Young Adult |
| Relative female wedges $\tau^\omega_{i,f}$, $i\ne T$ | Women's occupation shares among non-teachers | ACS |
| *Check, not a target* | Within-occupation wage variance | ACS, NLSY |

Restrictions kept from the old calibration: no education wedges, no wage wedges for men or for teaching, and one $A_i$ normalized. Two levels stay internal (Table 3): the level of $A_i$ relative to teaching pay, and the female outside option relative to teaching. The NLSY test and link files are pulled by `fetch_data.py --sources nlsy`.

**Why this block does not need the equilibrium.**

- **The non-teaching occupation choice is a static Roy model.** Non-teachers pick the occupation with the highest $(1-\tau^\omega_{i,g})A_i\epsilon_i^\alpha$ ([Proposition 3](../../notes/spatial_teachers_model_eps0.md#8-occupational-choice)). Ability $z$, school quality $Q_l$, taxes, teacher pay, amenities, altruism and moving costs all cancel. Non-teaching shares depend only on relative $A_i$, $\sigma_\epsilon$ and the female wedges.
- **Non-teacher log wages are linear in log ability.** $\log y=[\log X_O^*+\alpha\log z+\log Q_l]/(1-\eta)+\text{const}$. This is essentially exact without goods moving costs and approximate with them. Spatial and teacher parameters shift every non-teaching occupation equally, so they drop out of comparisons across occupations.
- **Teaching only removes agents from the non-teaching pool.** Teachers come from agents with low $X_O^*$ and, because $\gamma<1$, somewhat lower $z$. Given $X_O^*$, which non-teaching occupation wins does not depend on teaching.
- **The resulting bias is small.** In a stylized simulation (five occupations, with the teaching threshold of [App. A.3](../../notes/spatial_teachers_model_eps0.md#a3-occupational-thresholds-under-destination-homogeneity)), ignoring teachers biases $\log A_i$ by at most 0.013 at a 2% teaching share (about the male share). At a 10% share the bias is 0.039. $\mathrm{Var}(\log z)$ among non-teachers is unchanged to three decimals. Ignore the bias for men and correct it for women.

**Caveats.**

1. **$z$ does not affect which non-teaching occupation people choose,** because it is common to all of them. Wage variances therefore separate $s_z$ from $\sigma_\epsilon$ only through the lognormal functional form, and $s_z$ becomes a residual. That residual also absorbs transitory earnings, measurement error, age and hours.
   - **Use covariances instead.** The parent–child covariance is $b^2\rho_z s_z^2$ and the sibling covariance is $b^2\rho_z^2 s_z^2$. So $\rho_z=\mathrm{cov}_{sib}/\mathrm{cov}_{pc}$ and $b^2s_z^2=\mathrm{cov}_{pc}^2/\mathrm{cov}_{sib}$. Covariances are robust to independent transitory noise.
   - **Siblings share more than $z$,** such as family and neighborhood, which biases $\rho_z$ upward.
   - **Parent–child earnings also carry the school-quality channel** that the spatial block estimates.
   - **Test scores avoid that channel,** but need a measurement model to put them in model units.
   - **Occupational mobility tables are not targets.** The model predicts no inheritance of *which* non-teaching occupation a person holds, a prediction that is counterfactual.
2. **Shares and mean wages may conflict.** With one $\sigma_\epsilon$ and a common education barrier, a higher $A_i$ raises both an occupation's share and its mean wage. In the data, large occupations such as office support, sales and production typically pay less than small professional ones. Time investment $s$ is also the same across non-teaching occupations ([Proposition 2′](../../notes/spatial_teachers_model_eps0.md#7-human-capital-investment)). Raw mean-wage gaps, which largely reflect schooling, would then load entirely onto $A_i$ and selection, and inflate $\sigma_\epsilon$. Check the ACS pattern before choosing among three options:
   - **Occupation-specific $\sigma_{\epsilon,i}$.** Proposition 3 still holds, but the code needs a vector `σϵ`.
   - **Shares plus within-occupation dispersion,** as the old Fréchet calibration used.
   - **Schooling-adjusted wage gaps.**
3. **The model must use the same occupations.** The benchmark code has two occupations, teaching and other. That leaves no relative $A_i$ to estimate, and shares cannot pin down $\sigma_\epsilon$. Estimates from $J$ non-teaching occupations also do not carry over to a two-occupation model, whose outside option is one lognormal rather than a maximum of $J$ lognormals. The code supports $I>2$.
4. **Order.** Wages load on ability through $\alpha/(1-\eta)$, so fix $\eta$ first.

## 2. Internal calibration

Write $\kappa_1=\bar\kappa$, $\kappa_2=\bar\kappa e^{\delta_\kappa}$, $B=(0,\Delta B)$ and $m_{12}=m_{21}=r_m\bar C$, with no separate utility moving cost. Fit all parameters jointly, because their effects overlap. Each row names the principal source of identification, not a one-to-one mapping.

### Table 3: Internal parameters and targets

| Parameter | Target | Value | Source |
|---|---|---|---|
| ***Teaching margin*** | | | |
| Teacher pay level $\bar\kappa$ | Male teaching share | 1.9% in 2010; update to the 2018 sample | ACS |
| Female outside-option level (old $\lambda_f$; not altruism) | Female teaching share | Pending | ACS |
| Teacher wage curvature $\gamma$ | Model-generated teacher wage dispersion | Pending; start at the old Fréchet value 0.83 | ACS |
| ***Spatial block*** | | | |
| Spatial pay gap $\delta_\kappa$ | Adjusted salary log gap | 0.01623 | [District moments](estimates/spatial_moments.md), 176 CZs |
| Amenity gap $\Delta B$ | Enrollment log gap | 0.3764 | District moments, 188 CZs |
| Teacher spillover $\beta$ | Earnings at 28 per teacher-quality SD-year | ≈1.3% | [CFR (2014b)](https://www.nber.org/system/files/working_papers/w19424/revisions/w19424.rev2.pdf) |
| Altruism $\lambda$ | WTP per SD of school-average scores | \$19.70/month (SE \$7.40) | [BFM (2007)](../../literature/Bayer-UnifiedFrameworkMeasuring-2007.pdf), Table 7, col. 4 |
| Goods moving cost $r_m=m/\bar C$ | Household-income proxy log gap | 0.2708 | District moments, 187 CZs |
| — (additional restriction) | Teacher FTE per pupil log gap | −0.01216 | District moments, 188 CZs |

**Status.** When complete, the internal block has 8 parameters and 9 moments. The harness currently fits $\vartheta=(\log\bar\kappa,\delta_\kappa,\Delta B,\beta,\lambda,r_m)$ to the male teaching share and the six spatial moments. $\gamma$ and the female level sit at code placeholders until Table 2 is estimated. The starting values are $\kappa=(0.75,0.9)$, $B=(0,0.1)$, $\beta=0.15$, $\lambda=0.70$ and $r_m=0.1813$; none has been validated.

**Target definitions.**

- **Teaching shares.** Use the benchmark age, employment and occupation rules. The old draft used ages 25–34 and split part-time workers' weights between market work and home production. Include home production in the denominator if it remains a choice.
- **Salary gap.** Adjusted for outside wages (CWIFT). It matches $\Delta\log[\kappa_l\,\mathbb E(h_T^\gamma\mid T,l)]$, which includes teacher composition, so it is not the pay-intercept gap. Prefer experience- or qualification-standardized schedules where available.
- **Enrollment.** With equal cohorts, enrollment shares equal population shares. Match the weighted log contrast, not pooled national totals.
- **CFR.** Requires an auxiliary teacher intervention that holds class size and location fixed (§3).
- **WTP.** BFM report 1990 housing user-cost dollars. Convert to a consumption share and the matching model quality change (§3). Households with children pay \$7.41 more (Table 8), so reweight to parents where feasible.
- **Income.** District median household income is only a proxy for model labor income: it differs in household vs. individual, median vs. mean, and total vs. labor income. Construct a comparable statistic or move this moment to validation. Population shares and capitalization do not measure the moving cost. Whether income sorting identifies $r_m$ must be checked with the Jacobian and profiles, not inferred from the parameter count.
- **FTE per pupil.** Matches model teacher headcount per student $n_l/S_l$, not class size or teacher human capital.

### Table 4: Untargeted validation

| Outcome | Gap (data) | Model counterpart |
|---|---|---|
| Child poverty | −8.88 pp | Low-income child share by location (mapping to be defined) |
| SEDA learning rate | 0.0031 SD/grade (SE 0.0030) | $\Delta\log Q/(\chi K(1-\eta))$ (§3) |
| SEDA score level | 0.2887 SD | Composition plus $Q$; descriptive |
| Local revenue share | 0.4726 | Diagnostic only: model taxes fund teacher payroll, not total spending |

## 3. Mapping evidence to model objects

**Teacher quality (CFR).** With teacher headcount $n_l$ and student mass $S_l=M_l/2$, the implemented technology is

$$
Q_l=\left[\frac{n_l}{S_l}\,\mathbb E(h_T^{\beta/\sigma}\mid T,l)\right]^\sigma .
$$

The relevant exponent is $\beta/\sigma$; teacher counts capture only the extensive margin. Efficient class-size allocation equalizes $h_T^\beta N(h_T)^{-\sigma}$ within a location, which leaves no equilibrium variation in teacher value added. So match CFR with an auxiliary intervention: raise one teacher's simulated log human capital by one SD for one year, holding class size and location fixed. For a non-teaching child who reoptimizes investment,

$$
\Delta\log y'\simeq\frac{\beta\,s_{\log h_T}}{K(1-\eta)}\;\longrightarrow\;\log(1.013).
$$

Implement the finite intervention and average over children. The mapping assumes one SD of teacher human capital equals one SD of causal effectiveness. Test scores, salary residuals and VA are not interchangeable measures of it ([Wiswall 2013](https://doi.org/10.1016/j.jpubeco.2013.01.006); [CFR 2014a](https://www.aeaweb.org/articles?id=10.1257/aer.104.9.2593)).

**Class size.** The anchor is [FOO (2013), Table V, col. 1](https://pure.uva.nl/ws/files/1912854/124460_394782.pdf#page=26): log wages fall 0.0063 (SE 0.0033) per extra pupil in average class size over grades 4–6. With three years of exposure and a mean class size of 24.357,

$$
\sigma_{\rm wage}\simeq(1-\eta)\frac{K}{3}(24.357)(0.0063)=0.2294
\quad\text{at }\eta=0.103,\ K=5.
$$

The mapped SE is 0.1202, so the 95% interval includes zero. Use wages rather than annual earnings, which include an hours response absent from the model. As a cross-check with STAR: after $d$ years of reducing class size from $N_0$ to $N_1$, a gain of $d_S$ student SDs implies $\sigma_{\rm ext}\simeq(1-\eta)\chi d_S K/[d\log(N_0/N_1)]$. Krueger's kindergarten comparison ($N_0=22.4$, $N_1=15.1$, $d=1$, $d_S\simeq0.20$) gives 0.2957. Recompute $\sigma$ whenever $K$ or $\eta$ changes. The local class-size association (−0.00224, SE 0.00782) does not identify $\sigma$.

**Earnings–achievement bridge.** $\Delta\log y_{\rm adult}=\chi\,\Delta a$, with $a$ in student test-score SDs. CFR (2014b, §III.B; online App. Table 3, col. 3) report about 12% higher earnings at 28 per score SD. That is a conditional association, not a causal return, and it suggests $\chi\approx0.12$ (or $\log1.12\simeq0.113$). CFR's nearby 0.13 figure is score SDs per SD of teacher VA, not $\chi$. The current 0.13 is therefore a provisional choice.

**SEDA.** The [district estimates](estimates/teacher_quality.txt) give a weighted learning-rate gap of 0.0031 SD/grade (SE 0.0030) and a score-level gap of 0.2887 SD. Within-CZ learning-rate reliability is 0.974, or 0.915 with adjusted SEs. Regressions with CZ effects give 0.0046 (0.0022) raw and −0.0042 (0.0022) with SES controls; report these separately. Learning rates measure relative growth, so compare the model through

$$
\Delta g\simeq\frac{\Delta\log Q}{\chi K(1-\eta)}.
$$

Growth also reflects families, peers and other inputs, so the observed gap is neither a causal teacher effect nor an upper bound on one. See the [SEDA documentation](metadata/SEDA_documentation_6.0.pdf) and the [validation study](https://edopportunity.org/papers/learning_rate_validation_nontechnical_summary.pdf).

**Willingness to pay.** Compute compensating consumption $a$, holding other location attributes fixed:

$$
\mu\log C+\lambda\,\mathbb E[f(h'_0)]=\mu\log(C-a)+\lambda\,\mathbb E[f(h'_1)].
$$

Translate BFM's score unit with the same $\chi$ bridge. Use the \$19.70 structural estimate. The \$17.30 hedonic estimate measures capitalization, and neither figure measures moving costs. BFM's income gradient (\$1.38/month per \$10,000) does not identify $\psi<1$, because log altruism already makes dollar WTP rise with income.

## 4. Estimation and sensitivity

Minimize $[\mathcal M(\vartheta)-\widehat{\mathcal M}]'W[\mathcal M(\vartheta)-\widehat{\mathcal M}]$ with diagonal inverse-variance weights, adding transport uncertainty to literature moments. Bootstrap whole CZs, rebuilding the moments. Bootstrap the Table 2 estimation in the same loop so its uncertainty carries through. Impose $0<\beta<1$, $\lambda\ge0$ and $r_m\ge0$.

- **Identification checks.**
  - Rank and scaled singular values of the moment Jacobian.
  - Whether $\lambda$, $\Delta B$ or $\beta$ can reproduce each other's moment effects.
  - A profile of $r_m$ that refits everything else, including $r_m=0$. A flat profile means weak identification, even with a positive point estimate.
- **Solution checks.** Equilibrium residuals, positive consumption, quadrature accuracy and multiple starting values.
- **External-block consistency.** At the internal optimum, recompute the Table 2 moments inside the full model, including teacher selection, the moving cost, the school-quality gap and sorting. If they differ from the Table 2 fit by more than sampling error, re-estimate Table 2 against the data net of the model-implied distortion, then refit.
- **Residential transitions.** Add a matched childhood-to-adulthood gross transition rate as a further moment. Income-conditioned transitions would also test the goods-cost mechanism ([audit, Route B](external_parameter_audit.md#route-b-add-residential-transitions-by-pre-move-resources)). Directional flows and stationary shares are not independent restrictions. Condition on pre-move resources, not realized adult earnings.

### Table 5: Sensitivity panel

| Input | Values | Notes |
|---|---|---|
| $r_m$ | 0, 0.05, 0.10, 0.1813, 0.30 | Until $r_m$ is estimated, refit the other parameters at each value; 0 turns the mechanism off. These are scenarios, not bounds |
| $\psi$ | 0, 0.5, 1 | |
| $\sigma_\nu/\mu$ | 1/7, 0.20, 0.25, 0.30, 0.40 | |
| $\sigma$ (at $K=5$) | 0.05, 0.10, 0.25, 0.50 | Refit $\beta$ too, since the aggregator depends on $\beta/\sigma$. $\sigma=0$ needs a different allocation rule |
| $\chi$ | 0.10–0.20 | |
| $K$, $\eta$ | Alternatives | Recompute the class-size mapping |
| Table 2 parameters | Bootstrap draws | |

**Interpretation limits.**

- **Geography.** The model bundles residence, work and school, so residential moves and teacher job switches are different margins. [CWIFT](https://nces.ed.gov/programs/edge/economic/teacherwage) adjusts outside wages, not living costs.
- **Finance.** Local taxes fund teacher payroll only. The salary/local-revenue association (0.056) is not a pass-through elasticity. A finance extension would need payroll-specific transfers and a budget-to-pay rule.
- **Missing mechanisms.** The model omits effort responses ([Biasi 2021](https://www.aeaweb.org/articles?id=10.1257/pol.20200295)), teacher–student comparative advantage ([Biasi–Fu–Stromme](https://www.barbarabiasi.com/uploads/1/0/1/2/101280322/biasi_fu_stromme.pdf)) and childhood exposure timing ([Chetty–Hendren 2018](https://opportunityinsights.org/paper/neighborhoodsi/)).

The empirical test is whether the model jointly explains pay, sorting and a small school-effectiveness gap. A large teacher contribution to the score-level gap should be tested, not imposed.

## 5. To-do

- [ ] **T1 — Fix the model configuration.**
  - Choose the occupation set, which must match Table 2, and decide whether home production is a choice.
  - Settle investment accounting and the map from $s$ to schooling (the §1 feasibility check).
  - Fix reference units: $\bar C$, $\bar h$, $\chi$, $K$.
- [ ] **T2 — Estimate the occupations-and-ability block (Table 2).**
  - Build ACS shares and wage moments by gender and occupation, and state how these national targets relate to the within-CZ benchmark.
  - Build NLSY parent–child and sibling covariances.
  - Choose the wage moments (caveat 2) and fit the Roy simulator.
- [ ] **T3 — Finish the spatial targets.**
  - Resolve the income statistic.
  - Audit district pay and FTE counts.
  - Convert the literature targets (CFR intervention, BFM units).
  - Build a residential-transition target.
- [ ] **T4 — Implement matching model moments** for every target in Tables 2–4, including the CFR intervention and the WTP calculation.
- [ ] **T5 — Freeze the target file.**
  - For each target, record its definition, sample, weights and uncertainty. Use the whole-CZ bootstrap, and keep literature transport uncertainty separate from sampling error.
  - Remove stale structural labels from `data_estimate.py`, the generated reports and the README.
- [ ] **T6 — Verify and run.** Run the §4 checks on a trial parameterization, then the joint fit, the external-block consistency check and the sensitivity panel.

**Extensions that need not block the benchmark:** teacher geocodes and ability by location; a residence/workplace split; payroll-specific transfers; pre-1990 cross-sections; an Opportunity Atlas crosswalk with an exposure extension.
