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

**Sensitivity.** The class-size, exposure, location-taste and WTP cases have been rerun at the exactly identified baseline and change no validation conclusion; the rest of the panel and the ACS 2016–19 refit have not. [Table 5 details](calibration_log.md#table-5-sensitivity-panel).

**Experiments (2026-10).** Three rounds leave the baseline unchanged.
- *First round.* Parameter and weighting changes do not resolve sorting; imposed residential sorting moves teacher allocation too strongly because teachers live where they teach. Schooling floors and teaching utility costs improve selected teacher moments but remain diagnostic prototypes, not an adopted specification. [Results and limitations](calibration_log.md#calibration-experiments-2026-10).
- *10-05.* A pre-registered comparison selects C4′: a residence/workplace split, a joint college and occupation choice, and priced suburban access. Household sorting remains the main miss. [Results](calibration_log.md#model-experiments-2026-10-05).
- *10-07.* The regular Roy model ($\gamma=1$) removes teachers' selection within occupation and raises their wage dispersion to the market's. It closes the baseline's teacher pay gap but overpays teachers under C4′. Local peer effects in school quality create no sorting; once sorting is imposed, SEDA growth bounds them. [Results](calibration_log.md#model-experiments-2026-10-07).

## Current direction and user preferences (2026-10-16)

The user provisionally accepts the latest recommendations, conditional on preserving all previous work, and explicitly prefers to err on the side of simplicity. College choice and licensure are substantial departures from the earlier model that the user is skeptical of; their fit improvements justify keeping them as experimental benchmarks, not treating them as required ingredients of the eventual model. The user remains open to new experiments to make progress.

Proceed first with a pre-specified measurement comparison (race composition and family income versus individual earnings), re-scoring existing fits under both the original and alternative concepts. Then isolate teacher-quality sorting and test the smallest useful allocation/hiring change. Include a simpler specification close to the September structure and comparisons with extensions removed where feasible. Retaining a mechanism should depend on what it explains, its independent evidence and its complexity; a lower aggregate score alone is insufficient. Explicit household groups, structural housing demand, college choice and licensure remain options to assess, not commitments. The September baseline and historical candidate designations remain unchanged.

**Mechanisms in plain terms.** Income-dependent sorting (Y) currently adds a suburban preference that rises with predicted log income: richer households are more inclined to live in the suburb. This is a reduced-form preference interaction, not a solved housing market. Teacher demand (D) currently fixes the relative number of teaching positions per pupil and solves an access/rationing term so workplace choices fill those positions. D0 imposes equal positions per pupil across locations; it does not derive that rule from an optimizing district budget or specify a separate hiring process for teacher quality. These are experimental devices whose necessity and interpretation remain under review.

**Preservation and readiness.** Clear to proceed with this diagnostic and experimental sequence without another general approval request. No new specification, target replacement or explicit group structure is adopted. Existing source and results remain intact; new experiments use separate directories and preserve the old configurations. The three documents changed here were archived and hash-verified in [`runs/_code/accepted-preferences-20261016/`](runs/_code/accepted-preferences-20261016/README.md), in addition to the earlier source/evidence checkpoint. This update launches no model runs. B1 remains sealed. [Detailed preferences and sequence](calibration_log.md#accepted-preferences-and-next-sequence-2026-10-16).

**Status after round 2026-10-16** ([log](calibration_log.md#round-2026-10-16-measurement-concepts-the-teacher-quality-route-and-a-simpler-comparison)). Steps 1 and 2 and step 3's simple comparison were carried out in a pre-registered round, and step 4 is under way. Step 3's hiring margin was not built because its pre-registered trigger failed: at S1·Y·D0, where the pre-registration evaluates it, teachers' residence gap in human capital (0.128) exceeds their workplace gap (0.101); at the selection they are 0.096 and 0.076. A margin that changes who fills positions would therefore have to act on where teachers live. As the pre-registration notes, that is household structure, which is the user's decision.
- *Selection.* Under the original targets, the pre-registered selection is **P·TB·Y·D0·ρ**, the first candidate of these rounds to meet every criterion (under the round's pre-registered A4 and A5 rules for fits without college and for overidentified fits). It is now the development candidate; it is not adopted, and the September fit stays the reference.
  - It is the core model with the S1 split, income-graded sorting (Y), the D0 staffing closure and the simple teacher block TB (a schooling requirement and an ability-graded teaching cost). College and licensure are dropped. CFR is read with ρ = 0.3, the correlation of teacher effectiveness with pay-relevant skill.
  - It has 12 free parameters and S_1015 114.6, against C4′'s 18 and 315.1. It passes the Nz = 13 check, the identification check (condition number 22.4) and an 8-start multistart.
- *Simplicity.* Replacing college and licensure with TB improves the score by 19.1 with 7 fewer parameters. Of that, 14.1 is in the teacher-pay block, which TB targets; without that block the score still improves by 5.0. Removing TB costs 404.5, all of it in teacher pay: teachers become negatively selected and underpaid. Removing Y costs 348 without college (190 with it), and removing D0 costs 659 (from S1·Y·D0). Removing ρ costs 17 at P·TB·Y·D0 and 14 in the college model with D0 (M·Y·D0·ρ → S1·Y·D0), but with R (M·Y·R·ρ → S1·Y·R) it lowers the score by 22.
- *Measurement.* Per-adult family income (I1) barely moves the targets. Ranking each child by a co-resident parent's own earnings (I2, a pre-registered sensitivity) flattens V6a to 0.59 of its target, as much as race does. Within race, V1, V3 and V6a are 0.44–0.59 of their all-race values, V4 is 0.35 and V6z 0.45.
  - Under the original targets the selection misses V4 and V6z by about 5 SE.
  - Within race it fits them, as do the other two ρ = 0.3 fits (at ρ = 1 the fits with a staffing closure overshoot V4 and growth instead). It overshoots V1 (+8.8 SE) and V6a (−5.5 SE, by construction, since V6a is targeted at its original value).
  - The race concept therefore decides how much income sorting the targets require.
- *Teacher quality.* With Y, teachers sort by their own predicted income, mostly through the occupational shock. At a fixed workplace, the selection's teacher pay contrast by residence is 3.4× the PUMS. Exempting teachers from Y reverses the sorting rather than removing it. Linear interpolation between the base and that exemption at the ρ = 1 fits S1·Y·D0 and P·TB·Y·D0 (arithmetic, not a model run) puts the workplace gap near zero, and growth and the fixed-workplace contrast near the data, at a teacher grading of roughly 60–80% of non-teachers'.
- *Decisions for the user.*
  - How to treat race: a within-race first stage, or a group dimension in the model.
  - Whether to build household structure, with spouses' income in location choice, as the teacher-sorting margin.
  - When to open B1.
- *Remaining step-4 checks* before adoption: sample and payroll alignment, the in-equilibrium first stage (the selection holds the first stage's 90/10, occupation shares and NLSY moments, but misses the gender gap and cell wages), and the full bootstrap.

## Accepted development direction (2026-10-15, after both rounds)

Historical plan and results, retained for provenance. The October 16 preferences above control new work, including reconsideration of the earlier instruction to retain college choice and licensure.

The user provisionally accepted the next-step recommendations, conditional on preserving previous work. The September exact fit remains the reference baseline, C4′ remains the designated development candidate, and S1·R·E remains a diagnostic benchmark. Acceptance authorizes development along the following sequence; it does not adopt a new specification or authorize opening sealed B1.

1. Retain college/occupation choice, licensure and estimated wage compression ($\gamma\approx0.61$), and carry the S1 teacher-supply repair into new candidates.
2. Test a small set of income-dependent residential sorting mechanisms: housing or amenity demand whose expenditure share varies with income, with housing quantities and endogenous rents where appropriate. The failed fixed-charge and bounded-access tests do not reject this broader mechanism.
3. Develop explicit teacher demand or hiring constraints (budgets, authorized positions or class-size rules) to replace the reduced-form workplace penalty $b_T$.
4. Before new runs, fix identifying moments and acceptance thresholds. Include teacher residence, teacher crossing and commuting direction as explicit criteria alongside income/ability sorting, staffing, teacher selection, achievement levels and growth. Keep the gender wage gap visible; distinguish targets from validation and report component errors alongside the aggregate score.
5. Estimate strongly sorted candidates at $N_z\ge9$ and check a finer grid (initially $N_z=13$). Check identification before the full bootstrap; complete sample/payroll alignment and in-equilibrium first-stage checks before final inference.

**Preservation and readiness.** A verified local checkpoint of current source, pre-update documentation and selected fit evidence is in [`runs/_code/accepted-direction-20261015/`](runs/_code/accepted-direction-20261015/README.md). Existing outputs, historical protocols and older snapshots remain intact. New mechanisms must preserve the existing configurations and use new output directories. We are clear to proceed with specification design and implementation under this plan; the immediate deliverable is the next round's mechanism definitions, target map and pre-specified comparison, followed by checks and SLURM runs. [Detailed plan and preservation rules](calibration_log.md#accepted-development-plan-2026-10-15).

**Status after round 2026-10-15c** ([log](calibration_log.md#round-2026-10-15c-income-graded-sorting-y-and-teacher-demand-d)). Items 1–5 above were followed in a pre-registered comparison. Item 5's identification checks wait for a selected specification.
- *What was tested.* Income-graded residential demand (Y, a reduced form of non-homothetic demand) and explicit teacher demand (D: positions per pupil follow a rule, and $b_T$ becomes the rationing price). Both are switches that are bitwise-inert when off. Five candidates were fitted exactly at $N_z=9$ and checked at $N_z=13$.
- *Selection.* No candidate meets every criterion, so **C4′ stays the development candidate**. The best without a severe failure is S1·Y·D0, at $S_{1015}$ 151.6 against C4′'s 315.1.
  - It matches the targeted parental-income gradient ($\zeta_y$ ↔ V6a) and, untargeted, most of the poverty gap.
  - 59% of its score gap (0.196) is school quality, which is the part SEDA growth rejects: it fails growth by z +6.0, through teacher quality following income sorting. It also fails teachers' net commuting marginally (z −3.1).
  - S1·Y·R scores lower (133.9) but fails A6 by 0.015 (ε22 3.715 > 3.7), a severe failure under the pre-registration. C4′ stays by default, although it fails A6 itself (ε22 13.2).
- *Teacher demand.* Y needs the demand closure: without it, suburban staffing collapses (FTE −0.153).
- *Race composition.* A diagnostic after the round found that much of the sorting block is between-race composition, which the model lacks.
  - The own-AFQT gradient pooled within race is −0.05 to −0.10, against −0.222 overall; it is ≈ 0 among non-Black non-Hispanic youth.
  - The score gap pooled within race is ≈ 0.10, against 0.289 overall. It is ≈ 0 for white students and 0.16–0.21 for Black and Hispanic students.
  - Household income and child-poverty gaps keep about three-quarters of their size.
  - This explains part of the V4 and V6z shortfall, but not the growth overshoot.
- *Next step at the close of the round.* Decide how the model treats race before adding mechanisms. The [proposal](calibration_log.md#next-round-proposed-after-2026-10-15c) is retained below; the October 16 acceptance and simplicity preference now govern its implementation.
  1. Re-score all fits under a race-consistent concept registered in advance. This is a diagnostic and data work only.
  2. Decide whether to model the group dimension. Teacher allocation likely needs it, because teacher mobility responds to student composition.

  Teacher-quality sorting is the next model margin after that. B1 stays sealed.

## Recommended next round (2026-10-15)

Historical plan written before the two October 15 rounds; the accepted direction above now controls new work. The audit described as immediate below is complete.

Keep the September exact fit as the reference and use **C4′ as the development candidate**, retaining its college/occupation choice and residence/workplace separation. The next milestone is credible household sorting alongside plausible teacher allocation. These are recommendations for new work, not completed experiments or adoption of a new baseline.

1. **Audit commuting and teacher supply first.** The saved elasticity near 31 is the response of the suburb/city teacher-count ratio to suburban pay; the literature's 1.4–3.7 estimates concern supply to individual employers. Separate own-location and cross-location responses, align the adjustment horizon and geography, and account for unobserved within-POWPUMA commuting before imposing a target. If excessive responsiveness remains, test workplace attachment or heterogeneous commuting costs using independent evidence.
2. **Test heterogeneous housing access.** Use the rental distribution to discipline access to inexpensive suburban housing, housing quantities and endogenous rents. Start with school quality fixed; proceed to full equilibrium if the mechanism produces meaningful sorting at empirically supported housing parameters. Keep the ability-dependent taste case as a diagnostic benchmark.
3. **Retain the college-choice mechanism and estimate wage compression.** Carry forward M3 with estimated $\gamma<1$. The completed $\gamma=1$ and contemporaneous peer-effect experiments do not justify making either the next priority.
4. **Specify acceptance criteria before running.** Evaluate earnings, parental-income and AFQT sorting, staffing, teacher selection, and both score levels and growth jointly. Keep the gender wage gap visible. Give new parameters explicit identifying moments, distinguish targets from validation, and reserve genuinely unused evidence.
5. **Complete inference after specification choice.** Align samples and payroll definitions, recheck the first stage in equilibrium, and validate the extended model's grids and identification before the full calibration bootstrap and literature-input sensitivities.

The [detailed sequence](calibration_log.md#recommended-sequence-2026-10-15) gives the elasticity definitions, evidence, decision points and remaining work. The immediate deliverable is the commuting-elasticity audit.

**Status (2026-10-15, both rounds done; [log](calibration_log.md#model-experiments-2026-10-15)).** Items 1, 2 and 4 are done; item 3 is kept; item 5 waits for a specification choice.

- *Audit.* C4′'s supply is hypersensitive: own-location elasticity 13 and small-district elasticity ≈68, because σw ≈ 0.017. The repair S1 (ρw = 1 plus a utility commuting friction fitted to the crossing share) gives a full-equilibrium own-location elasticity of 1.6.
- *Housing.* The bounded housing-access test fails its gate.
- *Refits.* No pre-registered candidate meets all acceptance criteria, and C4′ remains the development candidate.
  - Under realistic supply, matching the salary gap over-staffs the suburb.
  - Matching FTE instead removes the school-quality contribution to the score gap.
  - A utility moving cost flattens the own-AFQT mobility gradient.
  - A BA-graded suburban amenity matches the parental-income gradient but moves teachers with other BA holders.
  - A non-pecuniary staffing closure (teachers' value of a suburban job, $b_T$) matches salary and FTE together.
  - With that closure, the BA amenity gives the best score: S1·R·E, 105 against C4′'s 315.
    - It matches the parental-income and AFQT gradients, poverty, staffing and pay.
    - It fails the Nz = 9 grid check at its fitted parameters. Re-fitted at Nz = 9 its moments are unchanged, but the
      amenity moves from 1.39 to ≈1.47.
    - It is a weakly identified reduced form (amenity ≈ 1.4 utils). It overshoots earnings and scores, and it puts
      teachers' residences in the suburb at five times the data gap.
- *Persistence.* Raising intergenerational persistence barely moves sorting. Shared family inputs are therefore not the next lever; household sorting must come from location forces graded by income or education outside school quality.
- *Next.* Replace the reduced forms with a sorting force graded by income, not education, such as non-homothetic housing or
  amenity demand, judged against the teacher-residence moments. Add a teacher-demand side in place of $b_T$. Score
  sorted candidates at Nz ≥ 9. The sealed hold-out B1 is opened only on the user's instruction
  ([log §5](calibration_log.md#recommended-sequence-2026-10-15)).
