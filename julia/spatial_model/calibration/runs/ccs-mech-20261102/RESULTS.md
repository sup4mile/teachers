# CCS × mechanisms round 2026-11-02: results

## Sources
- Design: PREREG.md, including the §9 amendments, all made before any fit.
- Jobs and checks: LOG.md.
- Tables:
  - `scores.md`: every fit by variant, with blocks, Δ against the efficient comparator, the CCS and draw-luck diagnostics, factorial effects, parameters and sensitivities;
  - `scores.csv`;
  - `compare.md`: the compact comparison;
  - `selection.json`.
- Outputs: `runs/exp-ccs-mech/<variant>/<case>`.
- Code: `runs/_code/ccs-mech-r1`.

Conventions:
- Every scored moment is untargeted.
- S13 = Σ z² over the 13 untargeted moments.
- Δ = CCS fit − its efficient comparator. The comparator is the factorial fit for ρ = 1 and `effr03` for ρ = 0.3. `effr03` itself is compared with the factorial fit.
- All differences are computed from unrounded S13.
- An independent results review recomputed every score, selection, OLS effect and check from the raw outputs. This version applies its corrections.

**HSP and HSLP are one fit.** Wherever the licensure floor is slack at schedT = 1.74, they coincide: |Δϑ| ≤ 4.5e−11 and |Δz| ≤ 3.5e−9 in ccs1, ccs12 and the factorial. The same holds for every pair that differs only by L when S is present.
- In ccs12 the registered rule resolves a 9.7e−9 difference in S13 literally, so it selects HSLP.
- In ccs1 the same rule selects HSP.
- The checks and sensitivities therefore ran on **HSP in ccs1** and **HSLP in ccs12**. That changes only the low-schedT sensitivity, where the floor binds in HSLP alone (§ Checks).

## Headline

1. **Does CCS work technically with every mechanism? Yes.**
   - All 161 refits are exact: criterion ≤ 9.7e−17, except ccs1/HFP at 2.6e−13, which is within the 1e−10 rule.
   - Every fit passes the numerics and the pad. The draw rules reach at most 0.47 of the pad.
   - At the selections, Nz = 13 moves no scored moment by more than 0.0055 SE; the rule is 0.25.
   - All 15 pad/smoke checks pass.
   - The multistarts return to the same ϑ: max relative |Δϑ| ≤ 1.1e−9 over 8 starts.
   - The Jacobians are as well conditioned as the efficient model's: condition numbers 30.6 (ccs1/HSP) and 35.1 (ccs12/HSLP), against 35.5 for the factorial HSP.
   - The code tests pass. The per-draw Proposition 5 test (ii) failed 110/120 on its first run (6811184_1); it passed 120/120 after a warm second node sweep and a comparator fix (h̄ and mean consumption compared by numerator and denominator; jobs 6811390 and 6813362). See Departures.

2. **CCS + K (K_d 12) is close to neutral and combines additively with every mechanism.**
   - The CCS cost is +17.6 to +34.5 in S13 across all 64 combinations.
   - The mechanisms' effects are the efficient model's. The largest gap in a main effect is 6.8 (S: −176.6 against −169.8). The largest gap in an interaction is 10.8 (S×L +390.5 against +379.7).
   - Best fit HSLP (= HSP) 383.2; the simplest within 10 is HS, at 383.3. The factorial best is HSP 364.0.
   - The cost is mostly the first stage: ρ_pc 0.517 against 0.577.
     - At none, the first stage is 17.5 of the 32.2 cost and teacher pay 14.6.
     - At the best fit, the first stage is 17.3 of 19.2 (90%) and teacher pay 2.5.

3. **Simple CCS (K_d 1) costs 423–629 S13 points in every combination.**
   - The cost is mostly the first-stage block: Δfirst runs from +406.5 to +429.0. Of that block, ρ_pc alone (0.295–0.300 against 0.577, SE 0.014) contributes 390.6–403.5 in z², and the wage slope (0.25–0.26 against 0.213) 19.4–35.3.
   - The first-stage cost varies only a little across mechanisms. Main effects on Δfirst: S −9.3 and L −6.2 (through the wage slope), C −6.3 (through ρ_pc). No mechanism repairs it.
   - Best fit HSP 827.1 (simplest within 10: HS 827.2).
   - **S, L and C meet the pre-registered "mitigates the CCS cost" rule**, with main effects on ΔS13 of −86.9, −71.9 and −24.6. **The mechanism is z² convexity, not a repair.**
   - *Teacher pay.*
     - CCS lowers the teacher premia (men/women) by 0.140/0.128 without S or L, by 0.110/0.103 with S, and by 0.097–0.100/0.090–0.092 with L alone: a shift of 0.09–0.14 in every combination.
     - The same shift costs far more in z² when the premia are far from the data (efficient none: −0.72/−0.63 against 0.056/0.086) than once S or L has brought them close (efficient S: −0.15/−0.05).
     - By block, S's −86.9 is pay −66.7, first stage −9.3 and allocation −8.8. L's −71.9 is pay −62.3, first stage −6.2 and sorting −4.6. The pay block carries 77% (S) and 87% (L).
     - Holding the none shift fixed (−0.140/−0.128 in every combination), convexity alone gives pay main effects of S −66.1 and L −52.1, against −66.7 and −62.3 actual. So convexity explains essentially all of S's pay effect and about 85% of L's; the rest of L's is its smaller shift.
   - *Allocation.*
     - C's −24.6 and H's +11.4 are the same convexity in the allocation block, with main effects on Δalloc of C −23.5 and H +12.7.
     - CCS lowers the FTE, growth and payroll gaps by similar amounts: about 0.006, 0.002 and 0.006 at none and with C.
     - C's efficient gaps are far above the data (FTE +0.078 against −0.018), so the shift helps. H's are below it (FTE −0.038), so the shift hurts.
     - Across the 32 C combinations, Δalloc runs from −34.3 (CS = CSLP) to −1.5 (CHFL).

4. **CCS + θ (ρ = 0.3) adds nothing new to CCS.**
   - CCS + θ behaves like simple CCS: ΔS13 +436 to +622 against effr03.
   - CCS + θ + K behaves like CCS + K: +21 to +33.
   - The ρ = 0.3 concept in the efficient model (effr03 against the factorial) is a separate effect. There is no θ draw in effr03: the change comes from the CFR measurement acting through the recalibrated β (0.228 → 0.068–0.074).
     - Effect on S13: +2.0 at none, F and P; +1.2 to +1.6 at S and L; +5.2 at C; **−30.0 at H (an improvement)**; +28.5 and +28.6 at HS and HSP; +46.2 and +46.3 at CHFSLP and CHS.
     - At HSP, the +28.6 is sorting +19.7 (V3 +8.3, V6a +8.7, V1 +4.7) and allocation +8.8.
     - It is a real effect, not another local solution. HSP's β (0.073) lies within the range of the other effr03 fits, the fit took 4 LM iterations, and Nz 13 moves it by ≤ 0.0005 SE. No multistart was run at effr03.
   - Mechanism, at H: the concept shrinks the suburb's teacher-quality deficit (gap_logQ −0.016 → −0.004).
   - Mechanism, at HS: β falls 0.246 → 0.073 and composition sorting roughly halves (SI 0.103 → 0.055; Q-adjusted wage gap +0.0064 → −0.0025).

5. **CCS does not amplify any mechanism's income sorting.**
   - The largest CCS × m effect on SI is +0.0055 (L, ccs1) in ccs1 and +0.0031 (H) in ccs12, against the +0.05 threshold.
   - Simple CCS slightly *reduces* sorting: ΔSI averages −0.012 and ranges −0.023 to +0.003. Combinations with SI ≥ 0.10 fall from 40 (factorial) to 32 (ccs1; 40 in ccs12); HS/HSP drop from 0.103 to 0.099.
   - The highest SI is CHSP (= CHSLP) in both ρ = 1 variants: 0.303 (ccs1) and 0.322 (ccs12), against 0.320 in the factorial.
   - Simple CCS also lowers the SEDA score gap V4: none 0.060 → 0.042, CHS 0.260 → 0.209.

6. **Draw luck sorts only through H, and only a little.**
   - The adults' own-draw gap (suburb − city, in SD of ln r) is positive in every H combination and negative in every combination without H:
     - ccs1: +0.002 to +0.010 with H, −0.004 to −0.001 without;
     - ccs12: +0.001 to +0.003 with H;
     - ρ 0.3 variants: ≤ +0.004 with H.
   - With a fixed suburban goods cost, adults with luckier draws (higher incomes) are slightly more likely to live in the suburb.
   - The effect on wages is negligible. The Q-adjusted wage gap, adjusted for k(ln Q̄ + ln r), and the Q̄-adjusted gap, adjusted for k ln Q̄ only, differ by ≤ 0.0016 in ccs1 and ≤ 0.00014 in ccs12.

## Per variant (selections by the registered rules)

Best S15 equals best S13 in every variant.

| variant | exact | best S13 = best S15 (simplest within 10) | ΔS13 at best | highest SI | none: S13 (ΔS13; Δ blocks sort/alloc/pay/first) | mitigates (ΔS13 main < −10) | amplifies (CCS×m on SI ≥ +0.05) |
|---|---|---|---|---|---|---|---|
| ccs1 (simple) | 64/64 | HSP 827.1 (HS) | +463.1 | CHSP 0.303 | 1473.9 (+609.8; +13.8/−7.5/+174.5/+428.9) | C, S, L (convexity; headline 3) | none |
| ccs12 (+K) | 64/64 | HSLP = HSP 383.2 (HS) | +19.2 | CHSP 0.322 | 896.3 (+32.2; +1.2/−1.1/+14.6/+17.5) | none | none |
| ccs1r03 (+θ) | 11/11 | HSP 854.9 (HS) | +462.4 | CHS 0.243 | 1480.6 (+614.4; +15.0/−7.2/+177.6/+429.0) | (11 combinations: no OLS) | — |
| ccs12r03 (+θ+K) | 11/11 | HSP 413.5 (HS) | +21.0 | CHS 0.260 | 899.4 (+33.2; +1.9/−0.7/+14.6/+17.5) | (no OLS) | — |
| effr03 (efficient, ρ 0.3) | 11/11 | HSP 392.5 (HS) | +28.6 vs factorial | CHS 0.261 | 866.2 (+2.0 vs factorial) | — | — |

## Per mechanism: the four questions

- Columns (2) and (3): **OLS** = effect-coded main effect over the 64 fits of ccs1 or ccs12; **1−0** = single mechanism minus none.
- Column (1) is "exact in every variant" for all six mechanisms.

| m | (2) S13 effect: ccs1 OLS / 1−0; ccs12 OLS; ccs1r03 1−0; ccs12r03 1−0 (factorial OLS) | (3) CCS-cost effect on ΔS13: ccs1 OLS / 1−0; ccs12 OLS / 1−0; ccs1r03 1−0; ccs12r03 1−0 | (4) SI main: ccs1 / ccs12 / factorial; CCS×m on SI: ccs1 / ccs12 |
|---|---|---|---|
| C commuting | +65.7 / +61.2; +89.2; +64.1; +88.7 (+90.4) | −24.6 / −24.5; −1.1 / +0.3; −24.8; −0.3 | +0.193 / +0.203 / +0.204; −0.010 / −0.000 |
| H housing | −87.6 / +31.3; −100.4; −12.2; −20.3 (−99.0) | +11.4 / +19.3; −1.4 / +2.3; +7.8; −0.4 | +0.115 / +0.116 / +0.113; +0.002 / +0.003 |
| F fertility | −0.5 / −0.0; −0.8; −0.0; −0.1 (−0.8) | +0.3 / +0.0; +0.0 / +0.0; +0.0; +0.0 | −0.001 / −0.001 / −0.001; 0 / 0 |
| S schedule | −256.7 / −572.5; −176.6; −576.2; −434.3 (−169.8) | −86.9 / −150.8; −6.8 / −11.7; −154.0; −12.2 | +0.071 / +0.072 / +0.071; −0.000 / +0.000 |
| L licensure | −261.8 / −535.3; −195.2; −539.4; −401.7 (−189.8) | −71.9 / −146.8; −5.4 / −11.9; −150.0; −12.4 | +0.031 / +0.026 / +0.026; +0.006 / +0.001 |
| P peers | −0.0 / +0.2; −0.0; +0.2; +0.2 (−0.0) | −0.0 / −0.0; −0.0 / −0.0; −0.0; −0.0 | 0 / 0 / 0; 0 / 0 |

Reading:
- **F and P** are as inert under CCS as without it. ψa = 0.10 changes nothing, and neither does 0.25 (sensitivity: −0.2 S13).
- **C** keeps its efficient-model profile: an over-staffed suburb, a worse allocation block, and the largest single source of income sorting. Under simple CCS it trims the cost slightly, through allocation-block convexity.
- **H** is as useful under CCS + K as without CCS. Under simple CCS it raises the CCS cost a little (+11.4), through the same convexity with the opposite sign. With ρ 0.3 its single effect turns favourable (−12.2 / −20.3), as in effr03 (−20.0 single − none), so this comes from the ρ concept, not from CCS.
- **S and L** are by far the most useful mechanisms in every variant. Under simple CCS they are worth more, mostly through pay-block convexity.
- **Interactions.** In ccs12 they are the factorial's: S×L +390.5 against +379.7; C×S +161 against +161; C×H −89 against −87. Not so in ccs1: S×L is +523.5, 38% larger. The extra +143.8 is the S×L effect on ΔS13; L adds no pay-block convexity once S has closed the premia gap. C×S (+150) and C×H (−80) are close to the factorial's.

## Checks (PREREG §6) and sensitivity (§7)

- **Nz = 13** (15 re-solves, job 6822945): the largest scored-moment shift is 0.0055 SE.

  | variant | largest shift (SE) |
  |---|---|
  | ccs1 | 0.0046–0.0055 |
  | ccs1r03 | 0.0042–0.0051 |
  | ccs12 | 0.0013–0.0023 |
  | ccs12r03 | 0.0012–0.0016 |
  | effr03 | 0.0005–0.0013 |

  All pass, so no follow-up was needed.
- **Pad checks** (15, job 6822946): all PASS.
  - The stored ϑ re-solves to criterion ≤ 7.1e−17.
  - max |x| is 0.262–0.276 at pad 0.6 and 0.075–0.077 at pad 0.3.
  - No key is lost against the comparator.
  - The only ⚠ is the known smooth-floor barrier, with chosen mass ≤ 5.8e−23.
  - The effr03 checks have no CCS and are verification solves only.
- **Jacobian** (job 6822948):

  | fit | singular values | condition number |
  |---|---|---|
  | ccs1/HSP | 219 … 7.16 | 30.6 |
  | ccs12/HSLP | 233 … 6.64 | 35.1 |
  | factorial HSP (reference) | — | 35.5 |

  No direction is lost.
- **Multistart** (8 refits, job 6822947): all 8 converge (criterion ≤ 6.4e−18) to the central ϑ, with max relative |Δϑ| 1.1e−9 for ccs1/HSP and 7.5e−10 for ccs12/HSLP.
- **Sensitivity** (10 refits, job 6822949). All are exact and pass the numerics. ΔS13 against the central fit, from unrounded S13, with SI:

  | case | r_p 0.003 | r_p 0.035 | schedT 0.96 | schedT 2.57 | ψa 0.25 |
  |---|---|---|---|---|---|
  | ccs1/HSP | +62.4 (SI −0.01) | −45.5 (0.20) | +339.8 (0.01) | −57.0 (0.15) | −0.2 (0.10) |
  | ccs12/HSLP | +66.8 (−0.00) | −46.8 (0.20) | +34.4 (0.08) | −12.1 (0.15) | −0.2 (0.11) |
  | factorial HSP (reference) | +65.7 (−0.01) | −46.4 (0.20) | +270.6 (0.02) | −9.1 (0.15) | −0.2 (0.10) |

  - The r_p and ψa responses and all the SI responses match the efficient model's.
  - The schedT-low contrast between ccs1 (+340) and ccs12 (+34) is partly the selection artefact above: the floor binds only in HSLP. ccs12's +34.4 matches the factorial's HSLP (+33.5).
  - The like-for-like comparison is ccs1 HSP +339.8 against factorial HSP +270.6. Simple CCS makes schedT matter more in both directions (+340 / −57 against +271 / −9), again through the shifted premia.

## Surprises and caveats

- **The mitigation flags (C, S, L under simple CCS) satisfy the pre-registered rule, but they are z² convexity**, in the pay block for S and L and in the allocation block for C. No mechanism addresses CCS's actual failure: the mother–child ρ, and the wage slope, under K_d = 1. Those vary by at most about 22 in z² across all 64 combinations.
- **Simple CCS slightly lowers income sorting and the SEDA score gap.** It does not raise them. With K_d = 12 these effects are ≈ 0.
- **The ρ 0.3 concept alone changes the efficient model.** It is not uniformly worse: −30.0 at H, about +2 elsewhere, but +28.6 at HS/HSP and +46 at CHS/CHFSLP. At the H combinations it halves composition sorting (headline 4).
- **Numerical flags.** In 3 of 179 LM runs, a Jacobian point failed to solve (GE did not converge in 120 iterations) and its column was zeroed.
  - In ccs1/HFP the ΔB + h point failed at every Jacobian evaluation from iteration 2 on (5 times), so ΔB was effectively frozen. LM stalled at 2.6e−13 and stopped at iteration 10. The fit is still exact by the rule, and its ϑ equals ccs1/HF to 8.8e−6, since P is inert, so nothing changes.
  - ccs12/HSLP_ms4 and ccs1/HSP__sched_low each lost the γ column once; both ended exact.

## Departures from PREREG

- None in the design or the scoring.
- The §9 amendments were made after the phase-1 review and before any fit:
  - a C-aware retry rule;
  - Δ computed only against comparators that pass the numerics;
  - a machine-readable pad check;
  - guards against empty case lists.
- Code test (ii) failed 110/120 on its first run (6811184_1). Its two causes, both in the test, were a single cold node sweep and a ratio-of-mixtures comparator for h̄ and mean consumption. It passed 120/120 with a warm second sweep (6811390), and 120/120 again with numerator and denominator compared at 1e−8 on the run tree (6813362).
- No retry was needed: all 161 fits were exact at the first attempt.
- No Nz 13 or pad follow-up wave was needed.
- The checks and sensitivities ran at the literal rule's selections, HSP in ccs1 and HSLP in ccs12 (one fit with HSP; see the top).
