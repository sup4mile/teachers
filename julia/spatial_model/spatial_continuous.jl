# =============================================================================
# Spatial occupational-choice model with teacher spillovers and altruism.
# Stationary general equilibrium with L locations and I occupations (T = teaching),
# at zero parental cost share (ε = 0), with continuous idiosyncratic shocks.
#
# Solution method
#   1. Non-teaching collapse.  For i ≠ T the household problem depends on (i, ϵ_i)
#      only through income capacity X_{O,i} = Θ_{i,g} ϵ_i^α, Θ_{i,g} = (1−τω_{i,g})A_i.
#      X_O* = max_{i≠T} X_{O,i} is therefore a sufficient statistic, and since each
#      X_{O,i} is lognormal its CDF is ∏_i F_i(x).  The I-dimensional shock integral
#      reduces to two scalars, (ϵ_T, X_O*).
#   2. Continuous occupation margin.  Agents teach iff W_T(ϵ_T) > W_O(X_O*).  Both
#      values are increasing, so the threshold follows from inverting splines, and
#      every aggregate is a quadrature integral against the atomless shock densities.
#   3. Ability z is Rouwenhorst-discretised; the spatial distribution Φ_l(z) is the
#      stationary distribution of the migration × AR(1) kernel.
#   4. Education before occupation (`Params.college`, off by default; §3c).  Given (g, l, z)
#      and before the shocks ϵ are drawn, the young choose d ∈ {0, 1} (no BA / BA) by a logit
#      with scale σd on EV₁ − EV₀ − cC − cCz log z, EV_d = E_ϵ[max_i W_i | d].  Each branch is
#      solved as in 1–3 with its own returns and schooling constraint, so the household block
#      carries two sets of policies and P₁(g, l, z); everything that integrates over a cell
#      (location probabilities, Φ, aggregates, Λ, the moments) mixes them with the weights
#      P(d | g, l, z).  With `college` off nothing of this is reached.
#   5. Local peer effects (`Params.ψa`, off by default).  School quality carries a peer input P_l,
#      Q_l = (2H̃_T/M_l)^σ P_l^{ψa}.  By default P_l is the students' mean ability ā_l, log ā_l the mean
#      log z of the young schooled in l relative to its ergodic mean (`peer_log_ability`); `peer_input`
#      = :hc makes it the adults' human capital by residence, and `peer_total` the group's total
#      (number × geometric mean) instead of its mean (`peer_log_input`, specifications A–D).  log P
#      is a fourth outer aggregate.
#   6. Constant class size (`Params.ccs`, off by default; §3e, 2026-10-26).  Classes are equal-sized within a
#      location, N_l = S_l/n_T,l, and students are assigned to teachers at random, so school quality
#      q = Q̄_l r is a state of the young: ln Q̄_l = β E[ln h_eff | teacher in l] − σ ln N_l and ln r the
#      mean-zero teacher draw.  By Proposition 5 the node at (Q̄_l r, ϵ) is the node at Q̄_l with the effective
#      shocks ϵ̃_T = r^{1/α} ϵ_T, X̃_O = r X_O, so the nodes are solved once, on effective-shock grids widened by
#      a pad, and each node x = ln r of the n_q-point Gauss rule R_l (matched to the moments of ln r) only
#      shifts the shock laws (`ccs_law`); every integral over a cell is the R_l-mixture of the draw integrals
#      (`MixedMaps`, `over_draws`).  The outer state replaces H̃_T by (ln Q̄_l, the teachers' central moments of
#      ln h_eff) by location (`ccs_teacher_moments`).
#   The outer loop is a damped fixed point in (H̃_T, M, t), and log P when ψa ≠ 0; under `ccs`, in
#   (ln Q̄, teacher moments, M, t), with H̃_T still computed as the efficient-equivalent index.
#
# Notation: ε (\varepsilon) is the parental cost share; ϵ (\epsilon) is the
# idiosyncratic ability shock.
# =============================================================================

using LinearAlgebra
using Printf
using Optim
using Dierckx
using QuadGK
using Distributions

const C_FLOOR = 1e-6    # consumption below which smooth_log is linearised
const M_FLOOR = 1e-8    # population floor in Q_l (guards 0/0 if a location empties)
const TAX_CAP = 0.9     # local tax rates are clamped to [0, TAX_CAP]

"log(C), continued linearly below δ so the optimiser always sees a finite gradient."
@inline smooth_log(C, δ = C_FLOOR) = C > δ ? log(C) : log(δ) + (C - δ) / δ

is_log_kernel(ψ) = isapprox(ψ, 1.0; atol = 1e-10)

"""
    warm(h, ψ, href)

Warm-glow kernel f(h) = ((h/href)^{1−ψ} − 1)/(1−ψ), with the log limit at ψ = 1.
ψ < 1 makes f convex in log h, so the value of a better school rises with the
child's ability. `expm1` avoids cancellation as ψ → 1.
"""
@inline warm(h, ψ, href) = is_log_kernel(ψ) ? log(h / href) :
                           expm1((1 - ψ) * log(h / href)) / (1 - ψ)

# -----------------------------------------------------------------------------
# 1. Parameters
#
# `Params()` is the baseline, with both sorting mechanisms on: a goods-denominated
# moving cost (`mcost`; `τmove = 0`) and a linear warm glow (ψ = 0).
# `nesting_params()` switches both off and recovers the pre-sorting model.
#
# Seven extensions (calibration experiments, 2026-10), each off at its default, so the
# model above is reproduced exactly (an eighth, the joint education and occupation choice
# `college`, follows the list; it composes with `split`):
#   pres    recurring goods cost of living in l′, paid by every resident each period
#           (housing/access); it enters consumption and Ξ exactly like `mcost`
#   Bz      taste for location l′ per unit of log z (ability-graded location taste)
#   sT_min  minimum time share in teaching (a schooling requirement); binds when s_T* is lower
#   cT      non-pecuniary cost of teaching, in utils, subtracted from W_T
#   cTz     its ability gradient: W_T falls by cT + cTz log z (cTz < 0: teacher training
#           is less costly for higher z, which selects teachers positively on ability)
#   csub    Stone–Geary subsistence c̲ = csub·Cbar, a goods cost paid by every adult in every
#           location (utility μ log(C − c̲) over total consumption C); like `pres` it does
#           not scale with (s, e), so it enters consumption and Ξ, and Proposition 2′ holds.
#           `consumption` is the supernumerary C − c̲, the argument of the log
#   split   residence and workplace apart (§3b): children go to school where the household
#           lives, and every worker also chooses a workplace w, in a nested logit with the
#           residence outer (σν) and the workplace inner (σw = ρw σν ≤ σν). A worker's positions
#           are spread over the workplaces with shares P_w: pupils for teachers (`gr.PT`), `jobs`
#           for non-teachers. Teachers' pay and their taste `BT` depend on w, taxes on r, and
#           commuting (r ≠ w) costs the share `τk` of consumption-relevant income, a time cost.
#           The solution then carries the joint probabilities instead of `Pi_T`, `Pi_O`.
#           Teacher demand has no general-equilibrium closure of its own: the positions are fixed
#           shares of the pupils, and the calibration handles the demand side by fitting the pay
#           premium δκ to the FTE-per-pupil gap (`--activate gap_teachers_pp --deactivate
#           gap_salary`), which is equivalent at a calibrated point to a demand-side closure
#           (a workplace vacancy value, or the finance rule κ_w). With τk > 0 the position
#           weights are normalised within each residence, so even at equal pay (δκ = 0) teachers
#           do not locate in proportion to pupils (FTE gap +0.055 at the baseline ϑ, τk = 0.03,
#           ρw = 0.066, against 0 at τk = 0): the equal-pay gap is part of what δκ offsets.
#           Home producers (non-teachers whose best occupation is home production) enter the
#           workplace nest like every other non-teacher, an approximation: they pick a workplace
#           over `jobs` and lose τk of their output when they cross, about 0.3% of their income
#           at τk = 0.03 and a crossing share of 0.107
#
# The eighth extension, `college` (§3c), makes the BA a choice taken before the occupation:
#   sBA     BA schooling, as a time share (16/25): a floor for every BA holder and a cap for the rest
#           (0: no constraint); ρc  log credential premium of a BA in every non-teaching occupation;
#   ρT0     log pay penalty of a teacher without a BA;  ϑs  elasticity of a BA teacher's pay in the
#           time share, pay κ h^γ (s/sBA)^ϑs (the salary schedule);  Fc  tuition, in units of Cbar
#   cC, cCz utility cost of a BA and its gradient per unit of log z, common to every BA holder;
#   sT0     time-share floor of a teacher without a BA (0: none, the default; the non-BA teaching track's
#           mean schooling is 12.99 years = 0.52 in the ACS);  cT0  its non-pecuniary cost, in utils,
#           on top of cT + cTz log z (0: none): the margin of the track's size, so that ρT0 can carry
#           its pay gap to BA teachers (ACS 0.48), not its size
#   σd      scale of the college taste shock;  ωz  loading of teacher effectiveness on z:
#           h_eff = h_T z^{−(1−ωz)α} enters Q and the CFR experiment (1: as h_T; any `college`)
# A ninth extension, local peer effects (2026-10-07), composes with all of the above:
#   ψa      elasticity of human capital in the local peer ability ā_l: every student schooled in l
#           acquires h = Q_T,l ā_l^{ψa} (z ϵ)^α s^φ e^η, with Q_T,l = (2H̃_T/M_l)^σ the teacher index (`QT`),
#           so the peer factor is part of the school-quality index Q = Q_T ā^{ψa} (`build_grids`) and
#           reaches everything that Q reaches (the technology, the scores and SEDA growth of the
#           moments, the value of a location to parents). ā_l is the geometric
#           mean of the students' z relative to the economy's (`peer_log_ability`), so it moves quality
#           across locations, not its overall level. 0: off, and Q is the teacher index alone
#   peer_input, peer_total (2026-10-25)  the peer input P_l of Q = Q_T P^{ψa} (`peer_log_input`; the solver's
#           fourth aggregate `la` is log P_l). `peer_input = :ability` (the default) is the students' ā_l above;
#           `:hc` is the adults' human capital by residence, log P_l = E[log h | adult resident of l] −
#           E[log h | all adults]: every old agent who lives in l (teachers, market non-teachers and home producers,
#           both genders, both education branches under `college`, by residence r under the split, the childless
#           under the fertility shock), with her own h, acquired when young where she was born (h_T, or h_O =
#           Q (z ϵ_{i*})^α s^φ e^η; `adult_log_h`). The reference is the endogenous economy-wide mean, so the
#           factor moves quality across locations as ā does (a level reference would feed Q back into its own
#           level through h, explosively for ψa ≥ 1). `peer_total = true` makes P the group's total, number ×
#           geometric mean, relative to an average location's: it adds the scale term log(n_l/n̄) of the students
#           n_l = Σ_z Φ[z, l] (or log(m_l/m̄) of the adult residents). Specifications A (:ability, mean; the
#           default), B (:ability, total), C (:hc, mean), D (:hc, total). At ψa = 0 both switches are inert, bit
#           for bit, and so is the default at any ψa. Unlike A and C, the scale term also moves the level: the
#           students' mean log Q rises by ψa Σ_l s_l log(L s_l) ≥ 0 (s_l the shares; +0.0174 ψa at the baseline's
#           0.407/0.593). Adults equal students by location in a stationary equilibrium without the fertility shock,
#           so there D − C = B − A exactly. Under :hc the adults carry their origin's Q, so Δlog P ≈ κ Δlog Q +
#           composition with κ the net stayer share (≈ 0.5 at the baseline): Q gaps multiply by 1/(1 − ψa κ), and the
#           damped fixed point cannot converge once ψa κ ≥ 1 (ψa ≈ 2), whatever the damping (review 2026-10-25)
# Two further switches (2026-10-15), each off at its default, compose with all of the above:
#   κc      utility friction of commuting (`split` only): a worker who lives in r and works at w ≠ r loses
#           κc utils on top of the time cost τk. It is an additive term of the location values of the
#           nested node (`location_values`), constant in (s, e), so it does not touch consumption, taxes,
#           Ξ or the per-destination value derivatives (the node's g and H move only through the choice
#           probabilities π), and Proposition 2′ and the X_O* collapse hold (0: off). It acts inside the
#           workplace nest, so its scale is σw = ρw σν, not σν
#   Bd      extra amenity of residence r for BA holders, in utils, by location (`college` only; the
#           education branch d = 1 of `with_education`): added to the residence-level additive term of
#           the BA branch's node, so within a branch it enters location choice only, like `Bz`. It also
#           raises the inclusive value of every BA node, hence W_T, W_O and EV₁, so it moves the college
#           choice: a uniform shift of Bd is a shift of cC, which is why it is normalised to (0, ζ_BA) and
#           read jointly with cC. `hh.Λ`, `dLambda_z` and `school_quality_wtp` keep their meaning; the BA's
#           valuation of the suburb is ΔB + Bd[2] − Bd[1] (zeros: off)
# Two switches of round 2026-10-15c, each off at its default, compose with all of the above:
#   By      income-graded location taste, an income interaction of the location preference (as in Bayer,
#           Ferreira and McMillan 2007): a worker values residence r by an extra By_r (ln ŷ − lyref) utils,
#           where ŷ is the node's predicted pre-tax income (`predicted_log_income`): its income at its own
#           location (r = w = l), before the local tax, at the time share of Proposition 2′ with Ξ = 0 (within
#           the node's bounds) and the goods investment of the destination-homogeneous first-order condition.
#           ŷ is a constant of the node, so the term enters location choice only, like Bz and Bd: consumption,
#           taxes, Ξ, Proposition 2′ and the X_O* collapse are untouched (`income_taste`). `lyref` is a frozen
#           reference (the calibration's population mean of ln ŷ), so that By tilts the location values across
#           incomes and B keeps its meaning; it has no default (NaN) and must be given whenever By ≠ 0. A uniform
#           By would be a shift of B; it is normalised to (0, ζ_y) (zeros: off)
#   fte_rule  teacher-demand closure (`split` only, two locations; NaN: off): districts' positions per pupil
#           differ by the log gap `fte_rule` (0: a common pupil–teacher norm), and teachers' value of the
#           suburban workplace, BT = (0, b), is the shadow price that rations the positions: `solve_ge` solves b
#           as an aggregate of the outer fixed point, so that teachers working in each location per pupil differ
#           by exactly `fte_rule`. The solution carries the equilibrium b in `p.BT`. Off, BT is a parameter
# One diagnostic switch of round 2026-10-16, off at its default (NaN), composes with all of the above:
#   lyT_fixed  teachers' income-graded taste at one common income (`teach_node`, `income_taste`): when finite, the taste
#           term of every TEACHER node is By_r (lyT_fixed − lyref) instead of By_r (ln ŷ − lyref) with the node's own
#           ln ŷ, so By cannot sort teachers by their own pay (it is a constant per residence and branch for them);
#           non-teacher nodes are unchanged. It touches the location values of teacher nodes only: consumption, taxes,
#           Proposition 2′, `predicted_log_income` and `mean_predicted_log_income` do not move. Inert when By = 0.
#           A diagnostic (`--external lyT_fixed=-3.04`), not a calibrated parameter; NaN: off, and the model is
#           bit for bit the one without it
# One switch of 2026-10-18 (the fertility shock, §3d), off at its default (1), composes with all of the above:
#   pkid    probability that an adult has a child. The shock d ∈ {0, 1} is drawn at the transition, with the
#           location taste shocks ν and independently of everything else, so occupation, (s, e) and income are
#           chosen before it is known; it multiplies the altruism term: a parent values location l′ by
#           μ log C + Λ(z, l′) + (taste terms), a childless adult by μ log C + (taste terms). The young therefore
#           maximise the expected inclusive value pkid·V̄₁ + (1 − pkid)·V̄₀ (`fert_value_derivs`), whose derivatives
#           weight destinations by the mixed probabilities π = pkid·π₁ + (1 − pkid)·π₀, so Proposition 2′ holds
#           with Ξ at π and the X_O* collapse is untouched. Adults (workers, taxpayers, teachers) locate with π;
#           children are born to parents only, 1/pkid each, so the population is stationary and the young of the
#           next generation locate with the parents' π₁ (`parent_probs`, `stationary_phi`): Φ, M and the pupils
#           follow π₁, every worker-side aggregate follows π. The childless keep the taste terms (Bz, Bd, By) and
#           the amenity B. 1: everyone has a child, the model without the shock bit for bit
#   pkid_z  ability gradient of the fertility shock (2026-10-18b; needs pkid < 1): P(child | z) =
#           logistic(a + pkid_z x), x = (log z − m)/s_z the adult's own log ability in units of its ergodic SD (m, s_z
#           the ergodic mean and SD of the Rouwenhorst chain), so pkid_z is the logit slope per ability SD. z is set
#           at birth, before every choice, so the node's weight on the parent's inclusive value is P(child | z) and the
#           warm glow's level, common to every choice of the cell, still drops out of all choices. Every parent has
#           1/pkid children whatever her z, so children per adult rise with ability (children as a normal good,
#           through ability's loading on income). The intercept a solves ρ(diag(P(child | z)) Πz) = pkid, so the
#           stationary population is constant and the share of parents among adults is exactly pkid. The young's
#           ability distribution is then the Perron vector of diag(P(child | z)) Πz, tilted toward the more fertile
#           (`fertility_probs`, `gr.pk`), and the parents' arrays carry the weight P(child | z)/pkid (`parent_probs`).
#           0: P(child) = pkid for everyone, bit for bit the shock above
# One switch of 2026-10-19 (the education-based teacher salary schedule), off at its default (NaN), composes with
# everything above except `college` (whose BA teachers have their own schedule, `ϑs`):
#   schedT  total elasticity of a teacher's pay in her time share s under a salary schedule that pays for credentials:
#           pay κ_w h^γ (s/sBA)^{schedT − γφ}, so that pay ∝ s^schedT at given (Q, z, ϵ_T, e): schooling is paid by the
#           schedule's lane elasticity, the rest of h by the compression γ. The pay is normalised at the BA, s = sBA
#           (16/25), where it is κ_w h^γ as without the schedule. It acts through the node's extra elasticity
#           ϑs = schedT − γφ (`teach_schedule`, either sign), so the time share solves Proposition 2′ with b = μ schedT
#           (1 + Ξ) (`time_share`) and `pay_factor` / `teacher_logpay_shift` carry it into the wage moments. NaN: off,
#           and the model is bit for bit the one without it
# One switch of 2026-10-26 (constant class size, §3e), off at its default, composes with everything above except `college`:
#   ccs     constant class size with random teacher assignment, in place of the efficient class sizes of Proposition 1:
#           classes are equal-sized within a location, N_l = S_l/n_T,l (S_l = Σ_z Φ[z, l] the pupils, by residence =
#           school location under the split; n_T,l the teachers who work in l), and every student of l is assigned at
#           random, independently of (z, ϵ, g), to K_d = `ccs_draws` teachers, one per school block, so that
#           ln q = ln Q̄_l + ln r with ln Q̄_l = β μ_l − σ ln N_l, μ_l = E[ln h_eff | teacher works in l] (`effective_h`), and
#           ln r = (β/K_d) Σ_t (ln h_eff,t − μ_l), mean zero with cumulants κ_m(ln r) = β^m κ_m(ln h_eff | l)/K_d^{m−1}
#           (`draw_cumulants`). q is a state of the young, seen with (z, ϵ) before (s, e, i); Q = Q̄ × the peer factor.
#           Proposition 5: the problem at (g, l, z, ϵ, q) is the one at Q̄_l with the effective shocks ϵ̃_T = r^{1/α} ϵ_T and
#           X̃_O = r X_O, because q enters only through h = q (z ϵ)^α s^φ e^η. The nodes are therefore solved once, on the
#           effective-shock grids widened by a pad that covers |x| (|x|/α for ϵ_T), with `gr.Q` = Q̄ (no node function
#           changes), and the draw x = ln r only shifts the laws of the effective shocks (`ccs_law`): every integral over a
#           cell is the mixture over the n_q nodes of R_l, the Gauss rule matched to the moments of ln r (`gauss_from_moments`,
#           Golub–Welsch), and an integrand that depends on the law (E[Θ_{i*} | X̃_O*]) or on the true shock (ϵ = ϵ̃ e^{−x/α})
#           goes through `over_draws`. Parents value a location by the whole teacher distribution (`child_Ef` mixes R_l).
#           The outer fixed point carries (ln Q̄_l, the teachers' central moments of ln h_eff to order 2n_q) by location
#           (`ccs_teacher_moments`); H̃_T is still computed, the efficient-equivalent index (2H̃_T/M)^σ of the same teachers.
#           Under CCS `gr.dT` and `gr.dO` are `nothing`, so a consumer that has not been taught the draws throws instead of
#           integrating at r = 1. false: Proposition 1, the model without it bit for bit (`college` is rejected: its
#           choice would need P₁ by draw)
#   ccs_nq  number n_q of Gauss nodes of R_l (3; 1 is the mean-field model q = Q̄_l, which still differs from Proposition 1
#           in the level of Q̄ and in N_l)
#   ccs_draws  K_d, the independent teacher draws that compose q (1: one teacher for the whole of school; 12: one a year)
#   ccs_rho the CFR reliability ρ of measured teacher effectiveness (1: none): with ρ < 1 effectiveness is ln a = ln h_eff + ν,
#           which the randomly assigned students receive too, so the per-teacher variance gains κ₂(ν) = κ₂(ln h_eff | l)
#           (1/ρ² − 1), Gaussian, before the β/K_d scaling (the CCS cases set it to the calibration's `cfr_rho`)
# `occ_sample` records the occupational sample of the calibration's occupational block and is
# not used here. Nor is `move_util` (2026-10-15): it tells the calibration (`theta_to_params`) to read its
# moving-cost parameter r_m as a utility cost `τmove` (utils) instead of the goods cost `mcost` (r_m C̄); the
# model has both fields and treats them as always.
#
# Frozen normalisations, measured once at `nesting_params()` (Nz = 5,
# nϵT = nXO = 48) and re-checked by `sorting_mechanism_tests` in the script:
#   HBAR_BENCH  mean child human capital h̄ → warm-glow reference `href`
#   CBAR_BENCH  mean consumption C̄         → `Cbar`
#   MCOST_BENCH goods cost whose utility value at C̄ equals TAUMOVE_BENCH (μ = 1)
# -----------------------------------------------------------------------------
const TAUMOVE_BENCH = 0.20
const HBAR_BENCH    = 0.1685483614
const CBAR_BENCH    = 0.1619724217
const MCOST_BENCH   = (1 - exp(-TAUMOVE_BENCH / 1.0)) * CBAR_BENCH

"L×L matrix with `x` off the diagonal and zeros on it."
offdiag(x, L) = [l == lp ? 0.0 : Float64(x) for l in 1:L, lp in 1:L]

Base.@kwdef struct Params
    # occupations and technology
    T::Int             = 1             # index of teaching within 1:I
    A::Vector{Float64} = [NaN, 1.5]    # productivity of non-teaching occupations (A[T] unused)
    α::Float64 = 0.30                  # ability elasticity of human capital
    φ::Float64 = 0.40                  # time (s) elasticity of human capital
    η::Float64 = 0.20                  # goods (e) elasticity of human capital
    σ::Float64 = 0.25                  # curvature of the teacher-quality index Q
    γ::Float64 = 0.80                  # human-capital elasticity of the teaching wage
    β::Float64 = 0.15                  # H̃_T aggregates teachers' h^{β/σ}
    ψa::Float64 = 0.0                  # peer effect: elasticity of h in the local peer input P_l (0: off)
    peer_input::Symbol = :ability      # P_l: the students' ability (:ability, ā_l) or the adults' human capital by residence (:hc)
    peer_total::Bool   = false         # P_l is the peer group's total, number × geometric mean (false: its mean)
    κ::Vector{Float64} = [0.75, 0.9]   # teaching-wage shifter by location
    # preferences
    μ::Float64    = 1.00               # weight on log consumption
    λ::Float64    = 0.70               # altruism weight on the warm glow f(h')
    ψ::Float64    = 0.0                # warm-glow curvature (ψ = 1: log kernel)
    href::Float64 = HBAR_BENCH         # warm-glow normalisation (fixed, not solved in GE)
    # locations
    B::Vector{Float64}     = [0.0, 0.1]              # amenity by location
    σν::Float64            = 0.20                    # scale of the Gumbel location taste shock
    τmove::Matrix{Float64} = zeros(2, 2)             # utility cost of moving l → l'
    mcost::Matrix{Float64} = offdiag(MCOST_BENCH, 2) # goods cost of moving l → l'
    Cbar::Float64          = CBAR_BENCH              # consumption level `mcost` is denominated at
    pres::Vector{Float64}  = [0.0, 0.0]              # recurring goods cost of residence in l'
    Bz::Vector{Float64}    = [0.0, 0.0]              # taste for l' per unit of log z
    sT_min::Float64        = 0.0                     # minimum time share in teaching
    cT::Float64            = 0.0                     # non-pecuniary cost of teaching (utils)
    cTz::Float64           = 0.0                     # its gradient per unit of log z
    csub::Float64          = 0.0                     # subsistence c̲ paid by every adult, in units of Cbar
    # wedges by (occupation, gender)
    τω::Matrix{Float64} = [0.0 0.0; 0.0 0.1]   # labour-income wedge; income keeps (1 − τω)
    τe::Matrix{Float64} = [0.0 0.0; 0.0 0.0]   # education barrier; goods cost is (1 + τe)e
    # ability process and population
    ρz::Float64   = 0.9                # AR(1) persistence of log z
    σξ::Float64   = 0.20               # AR(1) innovation std of log z
    σϵ::Float64   = 0.30               # std of log ϵ (iid across occupations)
    Mtot::Float64 = 2.0                # total population, both genders
    # residence/workplace split (§3b); the fields below matter only if `split`
    split::Bool           = false      # workers choose a residence r and a workplace w
    σw::Float64           = 0.25       # scale of the workplace taste shock, 0 < σw ≤ σν (σw = ρw σν)
    τk::Float64           = 0.0        # share of consumption-relevant income lost to commuting (r ≠ w)
    jobs::Vector{Float64} = [0.5, 0.5] # non-teaching positions by workplace (shares of the jobs)
    BT::Vector{Float64}   = zeros(2)   # non-pecuniary value of teaching at workplace w (utils)
    κc::Float64           = 0.0        # utility friction of commuting (r ≠ w), on top of τk (0: none)
    fte_rule::Float64     = NaN        # teacher-demand closure: log gap in positions per pupil; BT[2] solved in GE (NaN: off)
    # joint education and occupation choice (§3c); `sBA` … `σd` matter only if `college`
    college::Bool = false              # the young choose a BA (d = 1) or not (d = 0) before ϵ
    sBA::Float64  = 16 / 25            # BA schooling: time-share floor of BA holders, cap of the rest (0: none)
    ρc::Float64   = 0.0                # log credential premium of a BA, every non-teaching occupation
    ρT0::Float64  = 0.0                # log pay penalty of a teacher without a BA
    ϑs::Float64   = 0.0                # elasticity of BA teachers' pay in s: κ h^γ (s/sBA)^ϑs
    Fc::Float64   = 0.0                # tuition of a BA, goods, in units of Cbar
    cC::Float64   = 0.0                # utility cost of a BA, common to every BA holder
    cCz::Float64  = 0.0                # its gradient per unit of log z
    σd::Float64   = 0.25               # scale of the college taste shock
    sT0::Float64  = 0.0                # time-share floor of a teacher without a BA (0: none)
    cT0::Float64  = 0.0                # extra non-pecuniary cost of teaching without a BA (utils)
    Bd::Vector{Float64} = [0.0, 0.0]   # extra amenity of residence r for BA holders (d = 1), utils; location choice only
    By::Vector{Float64} = [0.0, 0.0]   # taste for residence r per unit of predicted log income ln ŷ − lyref (utils); location choice only
    lyref::Float64      = NaN          # frozen reference of ln ŷ for `By`; must be set when By ≠ 0 (NaN: unset)
    lyT_fixed::Float64  = NaN          # diagnostic: ln ŷ of every teacher node in the income taste `By` (NaN: off, the node's own)
    pkid::Float64       = 1.0          # fertility shock (§3d): probability that an adult has a child, who alone values Λ (1: off)
    pkid_z::Float64     = 0.0          # its logit slope per ergodic SD of log z, mean held at pkid (§3d; 0: off, P(child) = pkid)
    schedT::Float64     = NaN          # salary schedule without college: teachers' total pay elasticity in s, pay κ h^γ (s/sBA)^{schedT − γφ} (NaN: off)
    # constant class size with random teacher assignment (§3e, 2026-10-26); `ccs_nq` … `ccs_rho` matter only if `ccs`
    ccs::Bool           = false        # equal classes within a location, N_l = S_l/n_T,l; q = Q̄_l r a state of the young (false: Proposition 1)
    ccs_nq::Int         = 3            # Gauss nodes of the draw rule R_l (1: mean field)
    ccs_draws::Int      = 1            # K_d, independent teacher draws composing q (1: one teacher; 12: one a year)
    ccs_rho::Float64    = 1.0          # CFR reliability ρ: the draw's teacher variance is κ₂(ln h_eff)/ρ² (1: no ν)
    # teacher effectiveness h_T z^{−(1−ωz)α} for Q and the CFR experiment (any `college`; 1 = h_T)
    ωz::Float64   = 1.0
    occ_sample::Symbol = :acs2009_13   # occupational sample of the calibration's block (metadata)
    move_util::Bool = false            # calibration metadata: ϑ's r_m sets τmove (utils), not mcost (goods)
end

"""
    Params(I, L; ...)

`I` occupations and `L` locations. Location 1 is the base; locations 2..L share
`B_other`, `κ_other` and the off-diagonal moving costs `τmove_off` (utility) and
`mcost_off` (goods). Non-teaching occupations share `A_other`; `τω_other` is the
wedge on occupation 2 for gender 2. Other fields can be passed as keywords; `jobs`
defaults to equal shares and `BT` to zero.

`mcost_off` is an absolute goods amount. An economy whose consumption level is far
from `CBAR_BENCH` should re-derive it with `redenominate_move_cost`.
"""
function Params(I::Int, L::Int;
                T::Int = 1,
                A_other::Float64 = 1.5,
                κ_base::Float64 = 0.75, κ_other::Float64 = 0.9,
                B_base::Float64 = 0.0,  B_other::Float64 = 0.1,
                τmove_off::Float64 = 0.0,
                mcost_off::Float64 = MCOST_BENCH,
                τω_other::Float64 = 0.1,
                kwargs...)
    A = fill(A_other, I)
    A[T] = NaN
    τω = zeros(I, 2)
    I ≥ 2 && (τω[2, 2] = τω_other)
    defaults = (; T, A, τω, τe = zeros(I, 2),
                κ = vcat(κ_base, fill(κ_other, L - 1)),
                B = vcat(B_base, fill(B_other, L - 1)),
                τmove = offdiag(τmove_off, L), mcost = offdiag(mcost_off, L),
                pres = zeros(L), Bz = zeros(L), Bd = zeros(L), By = zeros(L), jobs = fill(1.0 / L, L), BT = zeros(L))
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

"""
    nesting_params(; kwargs...)

The pre-sorting model: utility moving cost (`τmove = TAUMOVE_BENCH`, `mcost = 0`)
and log warm glow (`ψ = 1`, `href = 1`). This is the regression benchmark and the
point where `HBAR_BENCH` and `CBAR_BENCH` are measured.
"""
function nesting_params(; kwargs...)
    defaults = (; τmove = offdiag(TAUMOVE_BENCH, 2), mcost = zeros(2, 2),
                ψ = 1.0, href = 1.0, Cbar = 1.0)
    return Params(; merge(defaults, NamedTuple(kwargs))...)
end

"Copy of `p` with the named fields replaced."
with_params(p::Params; kwargs...) =
    Params(; merge(NamedTuple(k => getfield(p, k) for k in fieldnames(Params)),
                   NamedTuple(kwargs))...)

"""
    redenominate_move_cost(p, Cbar) -> Params

Convert the utility moving cost `p.τmove` into the goods cost m = (1 − e^{−τmove/μ})·C̄,
the loss whose utility value at consumption C̄ equals τmove. Returns `p` with
`τmove = 0`, `mcost = m` and `Cbar` recorded. `Params().mcost` is this conversion
applied to `nesting_params()` at `CBAR_BENCH`.
"""
function redenominate_move_cost(p::Params, Cbar::Float64)
    L = length(p.B)
    m = [(1 - exp(-p.τmove[l, lp] / p.μ)) * Cbar for l in 1:L, lp in 1:L]
    m[diagind(m)] .= 0.0
    return with_params(p; τmove = zeros(L, L), mcost = m, Cbar)
end

function check_params(p::Params)
    I, L = length(p.A), length(p.B)
    @assert I ≥ 2 && 1 ≤ p.T ≤ I "need a teaching and at least one non-teaching occupation"
    @assert length(p.κ) == L && size(p.τmove) == (L, L) && size(p.mcost) == (L, L) &&
            length(p.pres) == L && length(p.Bz) == L && length(p.Bd) == L && length(p.By) == L &&
            size(p.τω) == (I, 2) && size(p.τe) == (I, 2) "Params dimensions are inconsistent"
    @assert all(iszero, diag(p.mcost)) "mcost must have a zero diagonal"
    @assert 0.0 ≤ p.sT_min < 1.0 "sT_min must lie in [0, 1)"
    @assert p.csub ≥ 0.0 "csub must be nonnegative"
    @assert isfinite(p.ψa) "ψa must be finite"
    @assert p.peer_input in (:ability, :hc) "peer_input must be :ability (the students' ability) or :hc (the adults' human capital), got :$(p.peer_input)"
    @assert 0.0 ≤ p.τk < 1.0 "τk must lie in [0, 1)"
    @assert p.κc ≥ 0.0 && isfinite(p.κc) "κc must be nonnegative and finite"
    @assert all(isfinite, p.Bd) "Bd must be finite"
    @assert all(isfinite, p.By) "By must be finite"
    @assert all(iszero, p.By) || isfinite(p.lyref) "the income-graded taste By ≠ 0 needs its frozen reference lyref (--external lyref=…)"
    @assert !isinf(p.lyT_fixed) "lyT_fixed must be finite (NaN: off)"
    @assert 0.0 < p.pkid ≤ 1.0 "pkid, the probability of having a child, must lie in (0, 1] (1: off)"
    @assert isfinite(p.pkid_z) "pkid_z must be finite"
    @assert iszero(p.pkid_z) || p.pkid < 1.0 "pkid_z, the ability gradient of the fertility shock, needs pkid < 1"
    # the peer reference of `peer_log_ability` is the ergodic distribution of Πz, which pkid_z ≠ 0 tilts; the
    # adults' specifications (peer_input = :hc) are not implemented with it either, so the switch is excluded for all
    @assert iszero(p.pkid_z) || iszero(p.ψa) "pkid_z ≠ 0 is not implemented with peer effects (ψa ≠ 0)"
    if !isnan(p.fte_rule)
        @assert isfinite(p.fte_rule) "fte_rule must be finite (NaN: off)"
        @assert p.split && L == 2 "the teacher-demand closure (fte_rule) needs the split and two locations"
    end
    @assert all(≥(0.0), p.jobs) && isapprox(sum(p.jobs), 1.0; atol = 1e-12) "jobs must be nonnegative shares summing to one"
    if p.split
        @assert length(p.jobs) == L && length(p.BT) == L "the split needs one `jobs` and one `BT` entry per location"
        @assert 0.0 < p.σw ≤ p.σν "σw must lie in (0, σν]: σw = ρw σν with 0 < ρw ≤ 1"
    end
    @assert 0.0 ≤ p.sBA < 1.0 "sBA must lie in [0, 1)"
    @assert p.σd > 0.0 && p.Fc ≥ 0.0 && p.ϑs ≥ 0.0 && p.ωz ≥ 0.0 "σd must be positive; Fc, ϑs and ωz nonnegative"
    p.college && !iszero(p.ϑs) && @assert p.sBA > 0.0 "the salary schedule (s/sBA)^ϑs needs sBA > 0"
    if !isnan(p.schedT)
        @assert isfinite(p.schedT) && p.schedT > 0.0 "schedT, the schedule's pay elasticity in s, must be positive and finite (NaN: off)"
        @assert !p.college "schedT is the schedule without college; under `college` the BA teachers' schedule is ϑs"
        @assert p.sBA > 0.0 "the salary schedule schedT is normalised at s = sBA, which must be positive"
    end
    # §3e (2026-10-26): constant class size
    @assert !(p.ccs && p.college) "constant class size (ccs) is not implemented with `college`: the college choice would need P₁ by teacher draw"
    @assert p.ccs_nq ≥ 1 "ccs_nq, the number of Gauss nodes of the teacher-draw rule, must be at least 1"
    @assert p.ccs_draws ≥ 1 "ccs_draws, the number K_d of independent teacher draws, must be at least 1"
    @assert 0.0 < p.ccs_rho ≤ 1.0 "ccs_rho, the CFR reliability of teacher effectiveness, must lie in (0, 1] (1: no ν)"
    @assert 0.0 ≤ p.sT0 < 1.0 "sT0 must lie in [0, 1)"
    p.college && p.sBA > 0.0 && @assert p.sT0 ≤ p.sBA "the floor of non-BA teachers, sT0, may not exceed the BA schooling sBA (their cap)"
    # The X_O* collapse needs a common τe across non-teaching occupations.
    nonteach = [i for i in 1:I if i != p.T]
    @assert all(p.τe[i, g] == p.τe[nonteach[1], g] for i in nonteach, g in 1:2) "τe must be common across non-teaching occupations"
    return nothing
end

# -----------------------------------------------------------------------------
# 2. Discretisation
# -----------------------------------------------------------------------------
"Rouwenhorst discretisation of log z' = ρ log z + σξ, ξ ~ N(0,1). Returns (z, Π)."
function rouwenhorst(N::Int, ρ::Float64, σ::Float64)
    q = (1 + ρ) / 2
    Π = [q 1-q; 1-q q]
    for n in 3:N
        Πp = Π
        Π = zeros(n, n)
        @views Π[1:n-1, 1:n-1] .+= q     .* Πp
        @views Π[1:n-1, 2:n  ] .+= (1-q) .* Πp
        @views Π[2:n  , 1:n-1] .+= (1-q) .* Πp
        @views Π[2:n  , 2:n  ] .+= q     .* Πp
        @views Π[2:n-1, :]     ./= 2
    end
    zmax = sqrt(N - 1) * σ / sqrt(1 - ρ^2)
    return exp.(range(-zmax, zmax; length = N)), Π
end

"Stationary distribution π = πΠ of a row-stochastic matrix, by power iteration."
function stationary(Π; tol = 1e-12, max_iter = 10_000)
    π = fill(1 / size(Π, 1), size(Π, 1))
    for _ in 1:max_iter
        π_new = vec(π' * Π)
        converged = maximum(abs, π_new .- π) < tol
        π = π_new
        converged && break
    end
    return π ./ sum(π)
end

"""
    fertility_probs(p, z, Πz) -> Vector

P(child | z) at the ability nodes (§3d). Without the ability gradient (`p.pkid_z` = 0) every node has `p.pkid`. With it,
P(child | z) = logistic(a + pkid_z x), x = (log z − m)/s with m and s the mean and SD of log z under the ergodic
distribution of `Πz`, and the intercept a solves ρ(diag(P(child | z)) Πz) = pkid by bisection. Since every parent has
1/pkid children, ρ is the growth factor of the population (the z-marginal of the next generation's kernel is
diag(P(child | z)) Πz/pkid whatever the location choices), so that root keeps the population stationary and makes the
share of parents among adults, the mean of P(child | z) under the Perron vector, exactly pkid.
"""
function fertility_probs(p::Params, z, Πz)
    iszero(p.pkid_z) && return fill(p.pkid, length(z))
    erg = stationary(Πz)
    lz  = log.(z)
    m   = dot(erg, lz)
    x   = (lz .- m) ./ sqrt(dot(erg, (lz .- m) .^ 2))
    pk(a) = @. 1 / (1 + exp(-(a + p.pkid_z * x)))
    root(a) = maximum(real, eigvals(Diagonal(pk(a)) * Πz))
    lo, hi = -50.0, 50.0                          # ρ rises in a from 0 to 1
    while hi - lo > 1e-13
        mid = (lo + hi) / 2
        root(mid) < p.pkid ? (lo = mid) : (hi = mid)
    end
    return pk((lo + hi) / 2)
end

"Force a vector to be strictly increasing (for spline inversion)."
function make_increasing(v::AbstractVector)
    w = collect(float.(v))
    @inbounds for i in 2:length(w)
        w[i] = max(w[i], w[i-1] + 1e-12)
    end
    return w
end

"""
    build_grids(p; H̃T, M, t, la, Nz, nϵT, nXO, q_lo, q_hi)

Grids and shock laws at aggregates (H̃_T, M, t) and the local peer input `la` = log P_l (zeros by
default):
- `z`, `Πz`: Rouwenhorst ability nodes and transitions;
- `Q`: school quality, the index of the technology h = Q_l (z ϵ)^α s^φ e^η: the teacher-quality
  index `QT` = (2H̃_T/M)^σ times the peer factor P_l^{ψa}, whatever the specification of P
  (`peer_log_input`: by default the students' mean ability ā_l, `peer_log_ability`). With ψa = 0 (the
  default) `Q` equals `QT` and `la` has no effect; `la` is returned as given;
- `dT`, `ϵTgrid`: law of ϵ_T ~ LogNormal(−σϵ²/2, σϵ) and a log-spaced grid between its
  `q_lo` and `q_hi` quantiles. Teachers come from the upper tail of ϵ_T; a grid spaced
  in quantiles would leave that whole tail (above the 0.98 quantile at nϵT = 48) in one
  interval;
- `dO[i,g]`: law of X_{O,i} = Θ_{i,g} ϵ_i^α; `XOgrid[:,g]`: log-spaced grid between the
  `q_lo` and `q_hi` quantiles of X_O*;
- `lϵT`, `lXO`: the grids in logs, the coordinate of every spline and integral;
- `logΘ[i,g]`: log Θ_{i,g}, to recover h_O from X_O*;
- `PT`: the pupils' shares M_l/ΣM, the teachers' positions of the residence/workplace split
  (§3b; unused otherwise);
- `pk`: P(child | z) at the ability nodes (§3d, `fertility_probs`; `p.pkid` at every node without the
  ability gradient, and 1 without the shock).

Constant class size (`p.ccs`, §3e, 2026-10-26; the keywords `ccs` and `pad` are read only then): `QT` = Q̄ = exp(`ccs.lQT`)
(the teacher index of equal classes, ln Q̄_l = β μ_l − σ ln N_l, a state of the solver; without `ccs` it starts at the efficient
σ ln(2H̃_T/M)), and `ccs.c[l]` are the central moments of ln h_eff among the teachers who work in l, orders 0..2n_q in
`c[l][1:2n_q+1]` with `c[l][1]` = 1 (`ccs_teacher_moments`; `nothing`, or a zero variance, gives the degenerate draw).
The shock grids are the effective-shock grids, those above widened at both ends by `pad` (log units of r) with the base
spacing: ⌈pad/α/Δ⌉ nodes for ϵ_T and a common ⌈pad/Δ_g⌉ (max over g) for X_O*, each gender keeping its own Δ_g, so that
`pad` = 0 is the base grid itself and `nϵT`, `nXO` count the widened grids. The pad must cover max |x| of the rules
(`solve_ge`'s `ccs_pad`). The top-level `dT` and `dO` are `nothing` (§3e's safety convention), and `gr.ccs` =
(; x, w, laws, base, pad, lQT, c, nq, draws, rho) carries, by location l, the nodes `x[l]` (of ln r) and weights `w[l]` of
R_l (`ccs_rule`), the laws `laws[l][k]` = `ccs_law(…, x[l][k], p)` of the effective shocks under draw k, and `base`, the
law at r = 1 (that of the true shocks). With `ccs` off the grids are built by exactly the operations they always were.
"""
function build_grids(p::Params; H̃T, M, t, la = zeros(length(p.B)), Nz = 5, nϵT = 64, nXO = 64,
                     q_lo = 1e-5, q_hi = 1 - 1e-6, ccs = nothing, pad = 0.0)
    I, L = length(p.A), length(p.B)
    z, Πz = rouwenhorst(Nz, p.ρz, p.σξ)
    # H̃_T is an integral of a spline and can undershoot zero when a location has
    # almost no teachers; the floor keeps Q real there. §3e: under constant class size Q_T is Q̄ = exp(lQT).
    QT = p.ccs ? exp.(ccs_lQT(p, ccs, H̃T, M)) : (2 .* max.(H̃T, M_FLOOR) ./ max.(M, M_FLOOR)) .^ p.σ
    Q  = iszero(p.ψa) ? copy(QT) : QT .* exp.(p.ψa .* la)

    dT     = LogNormal(-p.σϵ^2 / 2, p.σϵ)
    ϵTgrid = exp.(range(log(quantile(dT, q_lo)), log(quantile(dT, q_hi)), nϵT))

    nonteach = [i for i in 1:I if i != p.T]
    dO   = Matrix{LogNormal{Float64}}(undef, I, 2)
    logΘ = fill(NaN, I, 2)
    for g in 1:2, i in nonteach
        logΘ[i, g] = log((1 - p.τω[i, g]) * p.A[i])
        dO[i, g]   = LogNormal(logΘ[i, g] - p.α * p.σϵ^2 / 2, p.α * p.σϵ)
    end

    XOgrid = Matrix{Float64}(undef, nXO, 2)
    for g in 1:2
        lo = quantile_xo(dO, nonteach, g, q_lo)
        hi = quantile_xo(dO, nonteach, g, q_hi)
        XOgrid[:, g] = exp.(range(log(lo), log(hi), nXO))
    end
    p.ccs && return ccs_grids(p, (; z, Πz, Nz, Q, QT, la = collect(Float64, la), t, I, L, ϵTgrid, XOgrid, logΘ, nonteach,
                                    PT = max.(M, M_FLOOR) ./ sum(max.(M, M_FLOOR)), pk = fertility_probs(p, z, Πz)), ccs, pad)
    return (; z, Πz, Nz, Q, QT, la = collect(Float64, la), t, I, L, dT, ϵTgrid, dO, XOgrid,
              lϵT = log.(ϵTgrid), lXO = log.(XOgrid),
              logΘ, nonteach, nϵT, nXO, PT = max.(M, M_FLOOR) ./ sum(max.(M, M_FLOOR)),
              pk = fertility_probs(p, z, Πz))
end

"""
    quantile_xo(dO, nonteach, g, q)

q-quantile of X_O* = max_{i≠T} X_{O,i}, by bisection in log x. F* ≤ min_i F_i puts
it above every occupation's own q-quantile, and F* ≥ 1 − Σ_i(1 − F_i) puts it below
their (1 − (1−q)/n)-quantiles. With several occupations the lower tail of X_O* is far
above the lowest occupation's, so bounding the grid by individual quantiles would
leave many nodes where X_O* has no mass. Exact with one non-teaching occupation.
"""
function quantile_xo(dO, nonteach, g, q)
    isone(length(nonteach)) && return quantile(dO[nonteach[1], g], q)
    n = length(nonteach)
    a = log(maximum(quantile(dO[i, g], q) for i in nonteach))
    b = log(maximum(quantile(dO[i, g], 1 - (1 - q) / n) for i in nonteach))
    while b - a > 1e-12 * max(1.0, abs(a))
        c = (a + b) / 2
        prod(cdf(dO[i, g], exp(c)) for i in nonteach) < q ? (a = c) : (b = c)
    end
    return exp((a + b) / 2)
end

"CDF of X_O* = max_{i≠T} X_{O,i}: the product of independent lognormal CDFs."
@inline Fxo(gr, g, x) = prod(cdf(gr.dO[i, g], x) for i in gr.nonteach)

"Density of X_O*: f = F · Σ_i f_i/F_i, in one pass over the occupations."
function fxo(gr, g, x)
    F = 1.0
    ratio_sum = 0.0
    for i in gr.nonteach
        Fi = cdf(gr.dO[i, g], x)
        F *= Fi
        Fi > 0.0 && (ratio_sum += pdf(gr.dO[i, g], x) / Fi)
    end
    return F <= 0.0 ? 0.0 : F * ratio_sum
end

"""
    argmax_mean(f, gr, g, x)

E[f(log Θ_{i*}) | X_O* = x], where i* is the non-teaching occupation attaining the
max. P(i* = i | X_O* = x) ∝ f_i(x)/F_i(x). Exact (degenerate) with one
non-teaching occupation.
"""
function argmax_mean(f, gr, g, x)
    (; nonteach, dO, logΘ) = gr
    i0 = nonteach[1]
    isone(length(nonteach)) && return f(logΘ[i0, g])
    num = den = 0.0
    for i in nonteach
        w = pdf(dO[i, g], x) / max(cdf(dO[i, g], x), 1e-300)
        num += w * f(logΘ[i, g])
        den += w
    end
    return den > 0.0 ? num / den : f(logΘ[i0, g])
end

"E[Θ_{i*}^{−(1−ψ)} | X_O* = x]."
EΘpow(gr, g, x, ψ) = argmax_mean(lθ -> exp(-(1 - ψ) * lθ), gr, g, x)

"E[ϵ_{i*} | X_O* = x] = x^{1/α} E[Θ_{i*}^{−1/α} | x]: the shock in the occupation actually chosen."
Eϵ_nonteach(gr, g, x, α) = x^(1 / α) * argmax_mean(lθ -> exp(-lθ / α), gr, g, x)

# -----------------------------------------------------------------------------
# 3. The household's location problem at one node
#
# In both branches post-tax income in work location l' is
#     y_{l'} = c_{l'} · (s^φ e^η)^{γ_b},
# with γ_b = γ for teaching (wage κ h^γ) and γ_b = 1 otherwise:
#     teaching:      c_{l'} = (1−t_{l'})(1−τω_{T,g}) κ_{l'} [Q_l (z ϵ_T)^α]^γ
#     non-teaching:  c_{l'} = (1−t_{l'}) Q_l z^α X_O*
# Consumption is C_{l'} = y_{l'} − (1+τe)e − m_{l,l'} − p_{l'} − c̲ (p = `pres`,
# c̲ = `csub` C̄), and the value of l' is W_{l'} = μ log C_{l'} + Λ(z, l') + Bz_{l'} log z.
# A Gumbel shock over l' gives logit choice. Under `Params.split` the node is a `NestedNode`
# over residence and workplace (§3b) instead, and everything that follows is generic in the node.
# -----------------------------------------------------------------------------
"The household's location problem at one node: a `Node`, or a `NestedNode` under the split."
abstract type AbstractNode end

struct Node <: AbstractNode
    c::Vector{Float64}    # income coefficient c_{l'}
    γb::Float64           # income elasticity of s^φ e^η
    τe::Float64           # education barrier
    m::Vector{Float64}    # goods costs m_{l,l'} + p_{l'} + c̲ (moving, residence, subsistence)
    Λ::Vector{Float64}    # altruism term Λ(z, l') plus the location taste Bz_{l'} log z
    l::Int                # birth location
    smin::Float64         # lower bound on the time share (teaching: `sT_min`)
    smax::Float64         # upper bound on the time share (§3c, no BA: `sBA`; Inf = none)
    ϑs::Float64           # extra elasticity of income in s (§3c, BA teachers' salary schedule)
    Λ0::Vector{Float64}   # §3d: the additive term of the childless, `Λ` without the altruism term (unused at pk = 1)
    pk::Float64           # §3d: probability of having a child, P(child | z) = `gr.pk[zi]` (`Params.pkid`, `pkid_z`; 1: everyone, the term `Λ`)
end

"A node without the §3c fields: no upper bound on s and no extra elasticity (the baseline)."
Node(c, γb, τe, m, Λ, l, smin) = Node(c, γb, τe, m, Λ, l, smin, Inf, 0.0)

"A node without the fertility shock (§3d): everyone has a child, so the childless term is never used."
Node(c, γb, τe, m, Λ, l, smin, smax, ϑs) = Node(c, γb, τe, m, Λ, l, smin, smax, ϑs, Λ, 1.0)

"""
    fertile(p) -> Bool

Whether the fertility shock is on (§3d, `Params.pkid` < 1): some adults have no child and do not value Λ.
"""
fertile(p::Params) = p.pkid < 1.0

"""
    childless_term(Λz, zi, p, gr) -> Vector

The additive location term of a childless adult at ability node `zi` (§3d): the node's taste term `Λz` without the
altruism term Λ(z, ·), i.e. the ability-graded taste Bz log z alone (the BA amenity `Bd` and the income taste `By` are
added later to both terms, by `with_education` and `income_taste`). Without the shock it is `Λz` itself, never read.
"""
childless_term(Λz, zi, p::Params, gr) = fertile(p) ? p.Bz .* log(gr.z[zi]) : Λz

"""
    goods_cost(l, p)

Goods cost of each work location l' that does not scale with (s, e): the move from l,
the residence cost p_{l'} and the subsistence c̲ = `csub` C̄, which every adult pays in
every location. Every `Node` takes its `m` from here, so each of them enters consumption
and Ξ.
"""
goods_cost(l, p::Params) = p.mcost[l, :] .+ p.pres .+ p.csub * p.Cbar

"""
    teach_node(ϵT, zi, l, g, Λ, p, gr; d)
    nonteach_node(XO, zi, l, g, Λ, p, gr; d)

The node of the teaching problem at shock ϵ_T and of the non-teaching problem at X_O*: the
baseline node (`plain_teach_node`, `plain_nonteach_node`; under `p.split` the `NestedNode` of
§3b, `nested_teach_node`, `nested_nonteach_node`), and under `p.college` the node of education
branch `d` ∈ {0, 1} derived from it (`with_education`, §3c); `d` is ignored otherwise. With
`p.lyT_fixed` finite (round 2026-10-16, a diagnostic) the income-graded taste of a teacher node is graded at
that common ln ŷ instead of the node's own (`income_taste`); every teacher node of the model is built here.
"""
function teach_node(ϵT, zi, l, g, Λ, p::Params, gr; d = nothing)
    nd = p.split ? nested_teach_node(ϵT, zi, l, g, Λ, p, gr) : plain_teach_node(ϵT, zi, l, g, Λ, p, gr)
    nd = p.college ? with_education(nd, p, d, :teach) : nd
    return isnan(p.lyT_fixed) ? income_taste(nd, p, gr) : income_taste(nd, p, gr; ly = p.lyT_fixed)
end

function nonteach_node(XO, zi, l, g, Λ, p::Params, gr; d = nothing)
    nd = p.split ? nested_nonteach_node(XO, zi, l, g, Λ, p, gr) : plain_nonteach_node(XO, zi, l, g, Λ, p, gr)
    nd = p.college ? with_education(nd, p, d, :nonteach) : nd
    return income_taste(nd, p, gr)
end

"""
    teach_schedule(p) -> (scale, ϑ)

The education-based salary schedule without college (`p.schedT`, 2026-10-19): the extra elasticity ϑ = schedT − γφ of a
teacher's pay in s and the income scale sBA^{−ϑ}, so that pay κ h^γ (s/sBA)^ϑ ∝ s^schedT. (1.0, 0.0) when off, and
under `college`, whose nodes take their schedule from `with_education`.
"""
teach_schedule(p::Params) = isnan(p.schedT) || p.college ? (1.0, 0.0) :
    (ϑ = p.schedT - p.γ * p.φ; (p.sBA^(-ϑ), ϑ))

"The baseline teaching node (§3)."
function plain_teach_node(ϵT, zi, l, g, Λ, p::Params, gr)
    k = (1 - p.τω[p.T, g]) * (gr.Q[l] * (gr.z[zi] * ϵT)^p.α)^p.γ
    c = [(1 - gr.t[lp]) * p.κ[lp] * k for lp in 1:gr.L]
    Λz = Λ[zi, :] .+ p.Bz .* log(gr.z[zi])
    isnan(p.schedT) && return Node(c, p.γ, p.τe[p.T, g], goods_cost(l, p), Λz, l,
                                   p.sT_min, Inf, 0.0, childless_term(Λz, zi, p, gr), gr.pk[zi])
    scale, ϑ = teach_schedule(p)
    return Node(c .* scale, p.γ, p.τe[p.T, g], goods_cost(l, p), Λz, l,
                p.sT_min, Inf, ϑ, childless_term(Λz, zi, p, gr), gr.pk[zi])
end

"The baseline non-teaching node (§3)."
function plain_nonteach_node(XO, zi, l, g, Λ, p::Params, gr)
    k = gr.Q[l] * gr.z[zi]^p.α * XO
    c = [(1 - gr.t[lp]) * k for lp in 1:gr.L]
    Λz = Λ[zi, :] .+ p.Bz .* log(gr.z[zi])
    return Node(c, 1.0, p.τe[gr.nonteach[1], g], goods_cost(l, p),
                Λz, l, 0.0, Inf, 0.0, childless_term(Λz, zi, p, gr), gr.pk[zi])
end

function consumption(nd::Node, s, e, p::Params)
    y = nd.c .* (s^p.φ * e^p.η)^nd.γb
    iszero(nd.ϑs) || (y = y .* s^nd.ϑs)       # §3c: BA teachers' salary schedule, y ∝ s^{γφ + ϑs}
    return y .- (1 + nd.τe) * e .- nd.m
end

location_values(nd::Node, C, p::Params) = p.μ .* smooth_log.(C) .+ nd.Λ

"Teacher human capital h_T = Q_l (z ϵ_T)^α s^φ e^η."
teacher_h(ϵT, zi, l, s, e, p::Params, gr) = gr.Q[l] * (gr.z[zi] * ϵT)^p.α * s^p.φ * e^p.η

"Non-teaching income capacity (1−τω)A_{i*} h_O = Q_l z^α X_O* s^φ e^η."
nonteach_cap(XO, zi, l, s, e, p::Params, gr) = gr.Q[l] * gr.z[zi]^p.α * XO * s^p.φ * e^p.η

"Location values W and teacher human capital at a teaching node (`d`: the education branch under `college`)."
function W_teach(ϵT, zi, l, g, e, s, Λ, p::Params, gr; d = nothing)
    nd = teach_node(ϵT, zi, l, g, Λ, p, gr; d)
    return location_values(nd, consumption(nd, s, e, p), p), teacher_h(ϵT, zi, l, s, e, p, gr)
end

"""
Location values W and income capacity at a non-teaching node (`d`: the education branch under
`college`, where the capacity includes the credential premium of a BA, e^{ρc}).
"""
function W_nonteach(XO, zi, l, g, e, s, Λ, p::Params, gr; d = nothing)
    nd = nonteach_node(XO, zi, l, g, Λ, p, gr; d)
    return location_values(nd, consumption(nd, s, e, p), p), nonteach_cap(XO, zi, l, s, e, p, gr) * credential(p, d)
end

"Inclusive value σν·log Σ exp((W + B − τmove)/σν) and logit probabilities over work locations."
function logsum_probs(W, l, p::Params)
    v = [(W[lp] + p.B[lp] - p.τmove[l, lp]) / p.σν for lp in eachindex(W)]
    vmax = maximum(v)
    ev = exp.(v .- vmax)
    S = sum(ev)
    return p.σν * (vmax + log(S)), ev ./ S
end

"""
    time_share(γb, p, Ξ = 0)

Optimal time share (Proposition 2′). Combining the s- and e-FOCs leaves one
state variable, Ξ = Σ_{l'} π_{l'} m_{l,l'}/C_{l'}, where m_{l,l'} is every goods cost
that does not scale with (s, e): the moving cost plus the residence cost p_{l'} and
the subsistence c̲ (`goods_cost`).
With b = μφ(1+Ξ), s = bγ_b / (bγ_b + 1 − γ_b η). Without goods costs Ξ = 0 and s is
constant.
"""
@inline function time_share(γb, p::Params, Ξ = 0.0)
    b = p.μ * p.φ * (1 + Ξ) * γb
    return b / (b + 1 - γb * p.η)
end

"""
    time_share(γb, p, Ξ, ϑs)

Time share with a salary schedule (§3c): pay κ h^γ (s/sBA)^ϑs makes income ∝ s^{γ_bφ + ϑs}
e^{γ_bη}, and the same two first-order conditions as in Proposition 2′ give
s = b/(b + 1 − γ_bη) with b = μ(γ_bφ + ϑs)(1 + Ξ).
"""
@inline function time_share(γb, p::Params, Ξ, ϑs)
    b = p.μ * (p.φ * γb + ϑs) * (1 + Ξ)
    return b / (b + 1 - γb * p.η)
end

"""
    node_share(nd, p, Ξ) -> s

The time share of node `nd` at Ξ, with its constraints: the optimum of Proposition 2′ (with the
node's salary-schedule elasticity `ϑs`, if any) clamped to [`smin`, `smax`]. With `smin` and
`smax` as the baseline has them (0 or `sT_min`; Inf) this is `max(time_share, smin)`.
"""
@inline node_share(nd::AbstractNode, p::Params, Ξ = 0.0) =
    min(max(iszero(nd.ϑs) ? time_share(nd.γb, p, Ξ) : time_share(nd.γb, p, Ξ, nd.ϑs), nd.smin), nd.smax)

"Ξ = 0 time share of occupation `i`."
s_of(i::Int, p::Params, Ξ = 0.0) = time_share(i == p.T ? p.γ : 1.0, p, Ξ)

"""
    predicted_log_income(nd, p, gr) -> Float64

ln ŷ, the predicted pre-tax income of a node, the grading variable of the income-graded location taste
`Params.By` (round 2026-10-15c). It is the node's income at its own location (r = w = l: no commuting) before
the local tax, at the time share of Proposition 2′ with Ξ = 0 within the node's bounds (`node_share`) and the
goods investment of the destination-homogeneous first-order condition, (1 + τe) e = a y with a = ηγ_b. Income
is y = c s^{γ_bφ + ϑs} e^a (`consumption`), so
    (1 − a) ln ŷ = ln c_own + (γ_bφ + ϑs) ln s + a ln(a/(1 + τe)),
with c_own the node's income coefficient at its own location over 1 − t_l. c_own carries the wedge 1 − τω,
Q_l, z and the shock, κ_l and the pay factor for teachers, and the credential premium e^{ρc} for BA
non-teachers; ŷ depends on none of (s, e), so it is a constant of the node.
"""
function predicted_log_income(nd::AbstractNode, p::Params, gr)
    l   = nd.l
    own = nd isa NestedNode ? (l - 1) * gr.L + l : l
    a   = p.η * nd.γb
    s   = node_share(nd, p, 0.0)
    return (log(nd.c[own] / (1 - gr.t[l])) + (p.φ * nd.γb + nd.ϑs) * log(s) + a * log(a / (1 + nd.τe))) / (1 - a)
end

"""
    income_taste(nd, p, gr; ly = nothing) -> node

The node with the income-graded location taste `Params.By` (round 2026-10-15c): its residence-level additive
term `Λ` (by residence r) gains By_r (ln ŷ − lyref), with ln ŷ the node's predicted log income
(`predicted_log_income`), a constant of the node. Like `Bz` and `Bd` it enters location choice only. Returns
`nd` itself when `By` is zero (off), so the model without it is reproduced bit for bit.
The keyword `ly` (round 2026-10-16, `teach_node` with `Params.lyT_fixed`) replaces the node's own ln ŷ by a
given number; its default `nothing` evaluates `predicted_log_income` and performs the same floating-point
operations as before the keyword existed (a default of `predicted_log_income(...)` itself would be evaluated even
when `By` is zero).
"""
function income_taste(nd::Node, p::Params, gr; ly = nothing)
    all(iszero, p.By) && return nd
    τy = p.By .* ((ly === nothing ? predicted_log_income(nd, p, gr) : ly) - p.lyref)
    Λ = nd.Λ .+ τy
    # §3d: the childless have the income taste too (their term without the shock is never read)
    return Node(nd.c, nd.γb, nd.τe, nd.m, Λ, nd.l, nd.smin, nd.smax, nd.ϑs,
                nd.pk < 1.0 ? nd.Λ0 .+ τy : Λ, nd.pk)
end

# The `NestedNode` of §3b, defined below (the only other subtype of `AbstractNode`).
function income_taste(nd::AbstractNode, p::Params, gr; ly = nothing)
    all(iszero, p.By) && return nd
    τy = p.By .* ((ly === nothing ? predicted_log_income(nd, p, gr) : ly) - p.lyref)
    Λ = nd.Λ .+ τy
    return NestedNode(nd.c, nd.γb, nd.τe, nd.m, Λ, nd.bw, nd.lP, nd.l, nd.smin, nd.smax, nd.ϑs,
                      nd.pk < 1.0 ? nd.Λ0 .+ τy : Λ, nd.pk)
end

"Robust fallback for the goods margin: bounded Brent search over u = log e."
function solve_e_bracket(nd::Node, s, p::Params; lo, hi)
    V(u) = logsum_probs(location_values(nd, consumption(nd, s, exp(u), p), p), nd.l, p)[1]
    e = exp(Optim.minimizer(optimize(u -> -V(u), lo, hi)))
    C = consumption(nd, s, e, p)
    V̄, π = logsum_probs(location_values(nd, C, p), nd.l, p)
    return e, V̄, C, π
end

"""
    solve_e(nd, s, u0, p)

Maximise the inclusive value V̄ over u = log e at fixed s. Newton's method,
warm-started at `u0`, with analytic derivatives. With a = ηγ_b and mc = (1+τe)e,
    w_{l'} = ∂W_{l'}/∂u = μ(a y_{l'} − mc)/C_{l'}
    ∂V̄/∂u  = E_π[w]
    ∂²V̄/∂u² = Var_π(w)/σν + E_π[μ(a² y − mc)/C − w²/μ].
Falls back to `solve_e_bracket` if the log barrier binds, the objective is not
locally concave, or a step leaves [lo, hi]. Returns (e, V̄, C, π).
"""
function solve_e(nd::Node, s, u0, p::Params;
                 lo = log(1e-8), hi = log(1e4), maxit = 30, gtol = 1e-10)
    nd.pk < 1.0 && return solve_e_fert(nd, s, u0, p; lo, hi, maxit, gtol)     # §3d
    a  = p.η * nd.γb
    y1 = nd.c .* s^(p.φ * nd.γb)        # income at e = 1
    iszero(nd.ϑs) || (y1 = y1 .* s^nd.ϑs)   # §3c: BA teachers' salary schedule
    u  = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    for _ in 1:maxit
        e  = exp(u)
        y  = y1 .* e^a
        mc = (1 + nd.τe) * e
        C  = y .- mc .- nd.m
        V̄, π = logsum_probs(location_values(nd, C, p), nd.l, p)
        any(≤(C_FLOOR), C) && break
        w = @. p.μ * (a * y - mc) / C
        g = dot(π, w)
        abs(g) < gtol && return e, V̄, C, π
        H = (dot(π, w .^ 2) - g^2) / p.σν + dot(π, @. p.μ * (a^2 * y - mc) / C - w^2 / p.μ)
        H ≥ -1e-14 && break
        u -= g / H
        (isfinite(u) && lo ≤ u ≤ hi) || break
    end
    return solve_e_bracket(nd, s, p; lo, hi)
end

"""
    solve_node(nd, p, u0, s0; maxit = 12, relax_maxit = 200, tol = 1e-10)

Joint (e, s) solution at one node: alternate the Newton e-solve with the explicit
s(Ξ) update. The map contracts by about 0.03 per pass and exits after one pass when
there is no goods moving cost. A goods cost paid in every location (`csub`, `pres`) makes
its slope about −sΞ/((1+Ξ)μ), which tends to −1 as consumption nears zero, and the
undamped passes then oscillate. If `maxit` of them leave the residual at or above `tol`,
up to `relax_maxit` further passes take the relaxed step s ← s + ω (s(Ξ) − s) with
ω = 1/2, halved whenever the residual grows; e is re-solved at each pass as before. A
binding schooling requirement (`nd.smin` above the optimal s) holds s at `nd.smin`, and e
is optimal given it; likewise a binding cap (`nd.smax` below it, §3c) holds s at `nd.smax`.
Returns (e, V̄, s, residual, passes, relaxed), where `relaxed` says
whether the relaxed phase was needed (it never is at the baseline). Generic in the node:
under the split the destinations of π, C and m are the (residence, workplace) pairs.
"""
function solve_node(nd::AbstractNode, p::Params, u0, s0; maxit = 12, relax_maxit = 200, tol = 1e-10)
    s = min(max((isfinite(s0) && 0.0 < s0 < 1.0) ? s0 : node_share(nd, p), nd.smin), nd.smax)
    e = V̄ = NaN
    resid = Inf
    passes = 0
    while passes < maxit && resid ≥ tol
        passes += 1
        e, V̄, C, π = solve_e(nd, s, u0, p)
        Ξ = sum(π[lp] * nd.m[lp] / max(C[lp], 1e-8) for lp in eachindex(π))
        s_new = node_share(nd, p, Ξ)
        resid = abs(s_new - s)
        s = s_new
    end
    relaxed = resid ≥ tol
    ω = 0.5
    while passes < maxit + relax_maxit && resid ≥ tol
        passes += 1
        e, V̄, C, π = solve_e(nd, s, u0, p)
        Ξ = sum(π[lp] * nd.m[lp] / max(C[lp], 1e-8) for lp in eachindex(π))
        s_new = node_share(nd, p, Ξ)
        r = abs(s_new - s)
        r > resid && (ω /= 2)
        resid = r
        s = r < tol ? s_new : s + ω * (s_new - s)
    end
    return e, V̄, s, resid, passes, relaxed
end

# -----------------------------------------------------------------------------
# 3b. Residence and workplace: the nested location problem (`Params.split`)
#
# Children attend school where the household lives, but a worker chooses a residence r and a
# workplace w inside the CZ. The choice is a nested logit: the residence outer (scale σν) and
# the workplace inner (scale σw = ρw σν with 0 < ρw ≤ 1, which keeps it consistent with random
# utility; ρw = 1 is a flat logit over the L² pairs). The inner nest is weighted by positions: a
# worker of group o has positions spread over the workplaces with shares P^o_w, each with an iid
# taste shock, so that with W_{rw} = μ log C_{rw} + b_w + Λ_r (the last constant within a nest)
#     π_{w|r} = P^o_w exp(W_{rw}/σw) / Σ_{w'} P^o_{w'} exp(W_{rw'}/σw),
#     I_r     = σw log Σ_w P^o_w exp(W_{rw}/σw),
#     V̄       = σν log Σ_r exp((I_r + B_r − τmove_{l,r})/σν),        π_{rw} = π_r π_{w|r}.
# Teachers' positions are proportional to the pupils, P^T_w = S_w/ΣS (`gr.PT`), so that with equal
# pay, no commuting cost and b = 0 they locate in proportion to pupils; non-teachers' follow
# the data's job shares (`Params.jobs`). b_w = BT_w for teachers and 0 for non-teachers.
# With a commuting cost the weights, normalised within each residence, draw commuters toward the
# larger location, and the staffing at equal pay departs from the pupils': there is no
# general-equilibrium teacher-demand closure, and the calibration fits δκ to the FTE-per-pupil gap
# instead, which is equivalent at a calibrated point (see `split` in §1).
# The home producers (the non-teachers whose best occupation is home production) are in the
# workplace nest like every other non-teacher, an approximation: the X_O* collapse leaves no
# separate node for them, and they lose the share τk of their output when they cross (about 0.3% of
# home-production income at τk = 0.03 and a crossing share of 0.107).
# Teaching pay κ_w is paid at the workplace and the tax t_r levied at the residence.
# Commuting (r ≠ w) costs the share τk of consumption-relevant income, a time cost that leaves
# taxable income and measured wages alone:
#     teaching:      c_{rw} = (1−t_r)(1−τω_{T,g}) κ_w [Q_l (z ϵ_T)^α]^γ (1 − τk 1{r≠w})
#     non-teaching:  c_{rw} = (1−t_r) Q_l z^α X_O* (1 − τk 1{r≠w})
# so income stays homogeneous in (s, e), the goods costs m_{rw} = m_{l,r} + p_r + c̲ do not
# involve commuting, and Proposition 2′ holds with Ξ = Σ_{rw} π_{rw} m_{rw}/C_{rw}. Destinations
# are d = (r − 1)L + w: `reshape(x, L, L)` of a vector over d is indexed [w, r].
# -----------------------------------------------------------------------------
"Choice probability below which a destination at the consumption floor drops out of the Newton step of a `NestedNode`."
const P_FLOOR = 1e-14

"""
    NestedNode

The `Node` of the residence/workplace split. Its vectors over destinations d = (r − 1)L + w
(`c`, `m`, and the C, π of `solve_node`) have L² entries; `Λ` is by residence r, and `bw`, `lP`
by workplace w.
"""
struct NestedNode <: AbstractNode
    c::Vector{Float64}    # income coefficient c_d: tax of r, pay of w, commuting time
    γb::Float64           # income elasticity of s^φ e^η
    τe::Float64           # education barrier
    m::Vector{Float64}    # goods costs m_{l,r} + p_r + c̲ of destination d (by residence)
    Λ::Vector{Float64}    # altruism term Λ(z, r) plus the location taste Bz_r log z
    bw::Vector{Float64}   # non-pecuniary value of working at w (teachers: `BT`; non-teachers: 0)
    lP::Vector{Float64}   # log P_w: the position weights of the inner nest
    l::Int                # birth location
    smin::Float64         # lower bound on the time share (teaching: `sT_min`)
    smax::Float64         # upper bound on the time share (§3c, no BA: `sBA`; Inf = none)
    ϑs::Float64           # extra elasticity of income in s (§3c, BA teachers' salary schedule)
    Λ0::Vector{Float64}   # §3d: the childless's term by residence, `Λ` without the altruism term (unused at pk = 1)
    pk::Float64           # §3d: probability of having a child, P(child | z) = `gr.pk[zi]` (1: everyone)
end

"A `NestedNode` without the §3c fields: no upper bound on s and no extra elasticity (no college)."
NestedNode(c, γb, τe, m, Λ, bw, lP, l, smin) = NestedNode(c, γb, τe, m, Λ, bw, lP, l, smin, Inf, 0.0)

"A `NestedNode` without the fertility shock (§3d)."
NestedNode(c, γb, τe, m, Λ, bw, lP, l, smin, smax, ϑs) = NestedNode(c, γb, τe, m, Λ, bw, lP, l, smin, smax, ϑs, Λ, 1.0)

"A vector over residences r, stacked over the destinations d = (r − 1)L + w."
by_residence(x) = repeat(x; inner = length(x))

"A vector over workplaces w, stacked over the destinations d = (r − 1)L + w."
by_workplace(x) = repeat(x, length(x))

"Income coefficients c_d = (1 − t_r) f_w k (1 − τk 1{r ≠ w}) over the destinations d = (r − 1)L + w."
commuting_income(t, f, k, τk) =
    vec([(1 - t[r]) * f[w] * k * (r == w ? 1.0 : 1 - τk) for w in eachindex(f), r in eachindex(t)])

"""
    nested_teach_node(ϵT, zi, l, g, Λ, p, gr) -> NestedNode

The node of a teacher under the split: the workplace nest with the pupils' shares `gr.PT` as position
weights, pay κ_w at the workplace, the taste `p.BT`, and the commuting time cost `p.τk`. Teacher
demand has no equilibrium closure of its own (the calibration fits δκ to the FTE-per-pupil gap,
which is equivalent at a calibrated point); see §3b.
"""
function nested_teach_node(ϵT, zi, l, g, Λ, p::Params, gr)
    k = (1 - p.τω[p.T, g]) * (gr.Q[l] * (gr.z[zi] * ϵT)^p.α)^p.γ
    Λz = Λ[zi, :] .+ p.Bz .* log(gr.z[zi])
    isnan(p.schedT) && return NestedNode(commuting_income(gr.t, p.κ, k, p.τk), p.γ, p.τe[p.T, g],
                                         by_residence(goods_cost(l, p)), Λz,
                                         p.BT, log.(gr.PT), l, p.sT_min, Inf, 0.0, childless_term(Λz, zi, p, gr), gr.pk[zi])
    scale, ϑ = teach_schedule(p)       # 2026-10-19: the salary schedule without college
    return NestedNode(commuting_income(gr.t, p.κ, k * scale, p.τk), p.γ, p.τe[p.T, g],
                      by_residence(goods_cost(l, p)), Λz,
                      p.BT, log.(gr.PT), l, p.sT_min, Inf, ϑ, childless_term(Λz, zi, p, gr), gr.pk[zi])
end

"""
    nested_nonteach_node(XO, zi, l, g, Λ, p, gr) -> NestedNode

The node of a non-teaching worker under the split: the workplace nest with the position weights
`p.jobs`. Every non-teaching node is in it, home producers included: the X_O* collapse leaves no
separate node for the occupation that attains the max, so a home producer chooses a workplace like
any other non-teacher and loses the share τk of her output when she crosses. This is an
approximation (about 0.3% of home-production income at τk = 0.03 and a crossing share of 0.107). The
market-worker moments (`commute_shares`) leave home producers out; their consumption, their
job-access term in the residential choice and the inclusive-value gap (`option_values`) do not.
"""
function nested_nonteach_node(XO, zi, l, g, Λ, p::Params, gr)
    k = gr.Q[l] * gr.z[zi]^p.α * XO
    Λz = Λ[zi, :] .+ p.Bz .* log(gr.z[zi])
    return NestedNode(commuting_income(gr.t, ones(gr.L), k, p.τk), 1.0, p.τe[gr.nonteach[1], g],
                      by_residence(goods_cost(l, p)), Λz,
                      zeros(gr.L), log.(p.jobs), l, 0.0, Inf, 0.0, childless_term(Λz, zi, p, gr), gr.pk[zi])
end

function consumption(nd::NestedNode, s, e, p::Params)
    y = nd.c .* (s^p.φ * e^p.η)^nd.γb
    iszero(nd.ϑs) || (y = y .* s^nd.ϑs)       # §3c: BA teachers' salary schedule, y ∝ s^{γφ + ϑs}
    return y .- (1 + nd.τe) * e .- nd.m
end

"""
Indicator of commuting, 1{r ≠ w}, over the destinations d = (r − 1)L + w of a `NestedNode` with `L`
locations: the pairs that bear the utility friction `Params.κc`.
"""
commuting_indicator(L) = vec([r == w ? 0.0 : 1.0 for w in 1:L, r in 1:L])

"""
W_d = μ log C_d + b_w + Λ_r − κc 1{r ≠ w} of living in r and working at w, before the amenity B_r and
the moving cost. The friction `κc` (utils, `Params.κc`; 0 by default, when W is unchanged) is the same
on every commuting pair and constant in (s, e), so the per-destination derivatives w_d and dw_d of
`nested_value_derivs` are those of W without it (the node's g and H change only through the choice
probabilities π, which it moves), and consumption, taxes and Ξ do not know about it. It acts inside the
workplace nest, whose scale is σw = ρw σν: κc / σw, not κc / σν, measures how much commuting it removes.
"""
function location_values(nd::NestedNode, C, p::Params)
    W = p.μ .* smooth_log.(C) .+ by_workplace(nd.bw) .+ by_residence(nd.Λ)
    return iszero(p.κc) ? W : W .- p.κc .* commuting_indicator(length(nd.Λ))
end

"""
    nested_logsum_probs(W, l, nd, p) -> (V̄, π, πr, πw, I)

Inclusive value V̄ and choice probabilities of a `NestedNode` born in `l`, at the values `W` of
`location_values`: the position-weighted logsum over workplaces within each residence (inner
values `I[r]`, probabilities `πw[w, r]` = π_{w|r}), then over residences (`πr[r]`); `π =
vec(πw .* πr')` is the joint over the pairs d. Both logsums are shifted by their maximum, and a
workplace without positions has probability zero.
"""
function nested_logsum_probs(W, l, nd::NestedNode, p::Params)
    L    = length(nd.Λ)
    x    = reshape(W ./ p.σw, L, L) .+ nd.lP         # [w, r]
    xmax = maximum(x; dims = 1)
    ex   = exp.(x .- xmax)
    Sw   = sum(ex; dims = 1)
    I    = vec(p.σw .* (xmax .+ log.(Sw)))
    πw   = ex ./ Sw
    u    = (I .+ p.B .- p.τmove[l, :]) ./ p.σν
    umax = maximum(u)
    eu   = exp.(u .- umax)
    Su   = sum(eu)
    πr   = eu ./ Su
    return p.σν * (umax + log(Su)), vec(πw .* πr'), πr, πw, I
end

"""
    nested_value_derivs(nd, s, u, p) -> (; V̄, g, H, C, π, ok)

V̄ and its first two derivatives in u = log e at fixed s. With w_d = ∂W_d/∂u as in `solve_e`
and w̄_r = E_{π_{w|r}}[w], the derivative of the inner value I_r,
    g = ∂V̄/∂u   = E_π[w],
    H = ∂²V̄/∂u² = E_π[μ(a² y − mc)/C − w²/μ] + Σ_r π_r Var_{π_{w|r}}(w_{rw})/σw
                   + Σ_r π_r (w̄_r − g)²/σν.
A pair at the consumption floor with probability below `P_FLOOR` drops out of the sums (its
terms are of that order); with a larger probability the log barrier binds, and `ok` is false.
"""
function nested_value_derivs(nd::NestedNode, s, u, p::Params)
    L  = length(nd.Λ)
    a  = p.η * nd.γb
    e  = exp(u)
    y1 = nd.c .* s^(p.φ * nd.γb)                          # income at e = 1
    iszero(nd.ϑs) || (y1 = y1 .* s^nd.ϑs)                 # §3c: BA teachers' salary schedule
    y  = y1 .* e^a
    mc = (1 + nd.τe) * e
    C  = y .- mc .- nd.m
    V̄, π, πr, πw = nested_logsum_probs(location_values(nd, C, p), nd.l, nd, p)
    floored = C .≤ C_FLOOR
    any(floored .& (π .> P_FLOOR)) && return (; V̄, g = NaN, H = NaN, C, π, ok = false)
    w  = @. p.μ * (a * y - mc) / C
    dw = @. p.μ * (a^2 * y - mc) / C - w^2 / p.μ
    w[floored] .= 0.0
    dw[floored] .= 0.0
    wr = vec(sum(πw .* reshape(w, L, L); dims = 1))
    g  = dot(π, w)
    H  = dot(π, dw) + dot(πr, vec(sum(πw .* (reshape(w, L, L) .- wr') .^ 2; dims = 1))) / p.σw +
         dot(πr, (wr .- g) .^ 2) / p.σν
    return (; V̄, g, H, C, π, ok = true)
end

"`solve_e_bracket` for a `NestedNode`."
function solve_e_bracket(nd::NestedNode, s, p::Params; lo, hi)
    V(u) = nested_logsum_probs(location_values(nd, consumption(nd, s, exp(u), p), p), nd.l, nd, p)[1]
    e = exp(Optim.minimizer(optimize(u -> -V(u), lo, hi)))
    C = consumption(nd, s, e, p)
    V̄, π = nested_logsum_probs(location_values(nd, C, p), nd.l, nd, p)
    return e, V̄, C, π
end

"`solve_e` for a `NestedNode`: the same Newton step and fallback with the derivatives of `nested_value_derivs`."
function solve_e(nd::NestedNode, s, u0, p::Params;
                 lo = log(1e-8), hi = log(1e4), maxit = 30, gtol = 1e-10)
    nd.pk < 1.0 && return solve_e_fert(nd, s, u0, p; lo, hi, maxit, gtol)     # §3d
    u = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    for _ in 1:maxit
        d = nested_value_derivs(nd, s, u, p)
        d.ok || break
        abs(d.g) < gtol && return exp(u), d.V̄, d.C, d.π
        d.H ≥ -1e-14 && break
        u -= d.g / d.H
        (isfinite(u) && lo ≤ u ≤ hi) || break
    end
    return solve_e_bracket(nd, s, p; lo, hi)
end

# -----------------------------------------------------------------------------
# 3d. The fertility shock (`Params.pkid`, 2026-10-18)
#
# At the transition, with the location taste shocks ν and independently of them and of the agent's state, an adult
# learns whether she has a child (probability pkid) or not. A parent values destination l′ by W₁ = μ log C + Λ (the
# node's `Λ`: the altruism term plus the taste terms), a childless adult by W₀ = μ log C + Λ0 (the node's `Λ0`: the
# taste terms alone). Consumption is the same for both, since the occupation and (s, e) are chosen before the shock.
# The young therefore maximise the expected inclusive value
#     V̄ = pkid V̄₁ + (1 − pkid) V̄₀,      V̄_d = σν log Σ_l′ exp((W_d + B − τmove)/σν)   (under the split: §3b's logsums),
# whose derivatives in u = log e are the same mixtures of the two types' derivatives: g = E_π[w] at the mixed
# probabilities π = pkid π₁ + (1 − pkid) π₀, and H = pkid H₁ + (1 − pkid) H₀. The per-destination derivatives w do not
# depend on the type, so Proposition 2′ holds with Ξ = Σ π m/C at the mixed π, and nothing in the X_O* collapse changes.
# Adults (workers, taxpayers, teachers) locate with π (`node_probs`, `location_probs`); children are born to parents
# only, 1/pkid each, and go to school where their parent lives, so the next generation locates with π₁ (`parent_probs`).
# Off (pkid = 1) no node reaches this section.
# With the ability gradient (`Params.pkid_z` ≠ 0, 2026-10-18b) the node's weight is P(child | z) = `gr.pk[zi]`
# (`fertility_probs`) instead of pkid: z is set at birth, so the weight is common to the occupation and location choices
# of a cell and everything above holds node by node. Every parent still has 1/pkid children, so a node's children per
# adult are P(child | z)/pkid, and the parents' arrays (`parent_probs`) carry that weight on π₁: the next generation's
# kernel is then diag(P(child | z)/pkid) times the old one, whose Perron root is 1 by the choice of the intercept, and
# Φ is its Perron vector (`stationary_phi`).
# -----------------------------------------------------------------------------
"The node of a childless adult (§3d): the node with its childless term `Λ0` as its additive term, and pk = 1."
childless(nd::Node) = Node(nd.c, nd.γb, nd.τe, nd.m, nd.Λ0, nd.l, nd.smin, nd.smax, nd.ϑs, nd.Λ0, 1.0)
childless(nd::NestedNode) =
    NestedNode(nd.c, nd.γb, nd.τe, nd.m, nd.Λ0, nd.bw, nd.lP, nd.l, nd.smin, nd.smax, nd.ϑs, nd.Λ0, 1.0)

"""
    plain_value_derivs(nd, s, u, p) -> (; V̄, g, H, C, π, ok)

V̄ and its first two derivatives in u = log e at fixed s of a `Node` for one fertility type (its additive term `Λ`):
the terms of the Newton step of `solve_e`. `ok` is false when a destination is at the consumption floor, where
`solve_e` falls back to the bracket. `nested_value_derivs` is the `NestedNode` version.
"""
function plain_value_derivs(nd::Node, s, u, p::Params)
    a  = p.η * nd.γb
    y1 = nd.c .* s^(p.φ * nd.γb)            # income at e = 1
    iszero(nd.ϑs) || (y1 = y1 .* s^nd.ϑs)   # §3c: BA teachers' salary schedule
    e  = exp(u)
    y  = y1 .* e^a
    mc = (1 + nd.τe) * e
    C  = y .- mc .- nd.m
    V̄, π = logsum_probs(location_values(nd, C, p), nd.l, p)
    any(≤(C_FLOOR), C) && return (; V̄, g = NaN, H = NaN, C, π, ok = false)
    w = @. p.μ * (a * y - mc) / C
    g = dot(π, w)
    H = (dot(π, w .^ 2) - g^2) / p.σν + dot(π, @. p.μ * (a^2 * y - mc) / C - w^2 / p.μ)
    return (; V̄, g, H, C, π, ok = true)
end

"The value derivatives of one fertility type: `plain_value_derivs` or `nested_value_derivs`."
value_derivs(nd::Node, s, u, p::Params) = plain_value_derivs(nd, s, u, p)
value_derivs(nd::NestedNode, s, u, p::Params) = nested_value_derivs(nd, s, u, p)

"""
    fert_value_derivs(nd, s, u, p) -> (; V̄, g, H, C, π, π1, ok)

The expected inclusive value pk V̄₁ + (1 − pk) V̄₀ of a node with the fertility shock and its first two derivatives in
u = log e (§3d): the mixtures, with the node's weight pk, of the parent's (`nd` itself) and the childless's
(`childless(nd)`) `value_derivs`. `π` is the mixed probability vector, `π1` the parent's.
"""
function fert_value_derivs(nd::AbstractNode, s, u, p::Params)
    d1 = value_derivs(nd, s, u, p)
    d0 = value_derivs(childless(nd), s, u, p)
    q  = nd.pk
    return (; V̄ = q * d1.V̄ + (1 - q) * d0.V̄, g = q * d1.g + (1 - q) * d0.g, H = q * d1.H + (1 - q) * d0.H,
              C = d1.C, π = q .* d1.π .+ (1 - q) .* d0.π, π1 = d1.π, ok = d1.ok && d0.ok)
end

"Inclusive value and probabilities of one fertility type of a node at consumption `C` (`logsum_probs`, `nested_logsum_probs`)."
type_logsum(nd::Node, C, p::Params) = logsum_probs(location_values(nd, C, p), nd.l, p)
function type_logsum(nd::NestedNode, C, p::Params)
    V̄, π = nested_logsum_probs(location_values(nd, C, p), nd.l, nd, p)
    return V̄, π
end

"""
    fert_logsum(nd, C, p) -> (V̄, π, π1, π0)

The expected inclusive value of a node with the fertility shock at consumption `C` and its mixed, parents' and
childless's location probabilities (§3d).
"""
function fert_logsum(nd::AbstractNode, C, p::Params)
    V1, π1 = type_logsum(nd, C, p)
    V0, π0 = type_logsum(childless(nd), C, p)
    q = nd.pk
    return q * V1 + (1 - q) * V0, q .* π1 .+ (1 - q) .* π0, π1, π0
end

"`solve_e_bracket` with the fertility shock (§3d): a bounded Brent search over u = log e of the expected inclusive value."
function solve_e_bracket_fert(nd::AbstractNode, s, p::Params; lo, hi)
    V(u) = fert_logsum(nd, consumption(nd, s, exp(u), p), p)[1]
    e = exp(Optim.minimizer(optimize(u -> -V(u), lo, hi)))
    C = consumption(nd, s, e, p)
    V̄, π = fert_logsum(nd, C, p)
    return e, V̄, C, π
end

"""
`solve_e` with the fertility shock (§3d): the Newton step of `solve_e` on the expected inclusive value
(`fert_value_derivs`), with the same exits and fallback (`solve_e_bracket_fert`). Returns (e, V̄, C, π) with the
mixed π, which is what Ξ needs.
"""
function solve_e_fert(nd::AbstractNode, s, u0, p::Params; lo, hi, maxit, gtol)
    u = clamp(isfinite(u0) ? u0 : log(1e-3), lo, hi)
    for _ in 1:maxit
        d = fert_value_derivs(nd, s, u, p)
        d.ok || break
        abs(d.g) < gtol && return exp(u), d.V̄, d.C, d.π
        d.H ≥ -1e-14 && break
        u -= d.g / d.H
        (isfinite(u) && lo ≤ u ≤ hi) || break
    end
    return solve_e_bracket_fert(nd, s, p; lo, hi)
end

"""
    node_probs(nd, C, p; who = :all) -> Vector

Location probabilities of node `nd` at consumption `C`, over the destinations l′ (under the split, the pairs
d = (r − 1)L + w): `who = :all` those of all adults (with the fertility shock the mixed π of §3d), `:parents` the
parents' π₁ (where the children are born and go to school), `:childless` π₀. Without the shock the three are the
node's logit probabilities, computed as `location_probs` always has.
"""
function node_probs(nd::AbstractNode, C, p::Params; who::Symbol = :all)
    who in (:all, :parents, :childless) || throw(ArgumentError("who must be :all, :parents or :childless, got :$who"))
    nd.pk < 1.0 || return type_logsum(nd, C, p)[2]
    who === :parents   && return type_logsum(nd, C, p)[2]
    who === :childless && return type_logsum(childless(nd), C, p)[2]
    return fert_logsum(nd, C, p)[2]
end

# -----------------------------------------------------------------------------
# 3c. Education before occupation (`Params.college`)
#
# Given (g, l, z) and before the occupational shocks are drawn, a young agent chooses d ∈ {0, 1}
# (no BA, BA) by a logit with scale σd,
#     P₁(g, l, z) = Λ([EV₁ − EV₀ − cC − cCz log z]/σd),
#     EV_d = E_ϵ[max(W_T^d, W_O^d)] = ∫_teachers W_T^d + ∫_non-teachers W_O^d  (branch d's choice maps),
# so the cost cC + cCz log z falls on every BA holder, teachers and non-teachers alike. Both
# branches keep the full teach / don't-teach structure of §4:
#   d = 1  non-teaching income × e^{ρc}, a credential premium common to every non-teaching
#          occupation (home production included, so X_O* stays a sufficient statistic and the
#          first-stage shares are unchanged); teaching pay κ h^γ (s/sBA)^ϑs, the salary schedule
#          as a lane in s, so income ∝ s^{γφ+ϑs} and Proposition 2′ holds with the elasticity
#          γφ + ϑs (`time_share`); s ≥ sBA for teachers (max(sT_min, sBA)) and non-teachers
#          (Kuhn–Tucker, as `smin`); tuition Fc C̄ is a goods cost like m (in `Node.m`, so in Ξ);
#   d = 0  non-teaching as in §3; a teaching track without a BA at pay κ h^γ e^{−ρT0} and no
#          licensing floor by default (`sT0`: a floor of its own, e.g. the track's mean schooling 0.52);
#          s ≤ sBA, a cap (Kuhn–Tucker, `smax`; none if sBA = 0); the cost cT0 (utils) of teaching
#          without a BA, on top of cT + cTz log z, makes the track's size a margin of its own.
# The warm glow mixes the branches, E f(h′) = P₁ E[f | d = 1] + (1 − P₁) E[f | d = 0], and uses
# human capital net of the credential premium (a BA non-teacher's h is her income capacity over
# e^{ρc}; a teacher's h is her human capital, the pay factor being a wage).
#
# Household solution under `college`: `hh = (; college = true, br1, br0, P1, Λ, …)`. Branch d
# (`br1`, `br0`) carries the arrays of §6 (eT, sT, WT, hT, eO, sO, WO, capO, Λ, infeasT, infeasO),
# where `capO` is the INCOME capacity of non-teachers (with e^{ρc} for d = 1), and `P1[g, l, zi]` is
# P₁. The work-location probabilities are pairs (d = 1, d = 0). Everything that integrates over a cell
# loops over `cell_branches` with the weights P(d | g, l, z); the arrays have no top-level `eT`, …, so
# a consumer that has not been taught the branches throws instead of reading half of the economy.
# With `college` off none of this is reached, and cC = −∞ with ρc = ρT0 = ϑs = Fc = sBA = 0
# reproduces the baseline.
#
# Under `Params.split` (§3b) the two extensions compose: each branch's nodes are `NestedNode`s with
# the branch's income scale, tuition, time-share bounds and salary-schedule elasticity
# (`with_education`), the branch's location probabilities are the joint residence × workplace
# arrays [g, l, zi, k, r, w], the EV of the college choice is the nested inclusive value inside W
# (nothing changes in `college_values`), and the aggregates and the tax base mix the branches with
# P(d | g, l, z) by residence and workplace (`aggregates`). The floored-node cap of `choice_maps`
# applies in each branch (the branch arrays carry `infeasT`, `infeasO`).
# -----------------------------------------------------------------------------
"Education branches of the model: (1,) without college (the baseline), (1, 0) with it."
branch_ids(p::Params) = p.college ? (1, 0) : (1,)

"Household arrays of education branch `d`: `hh` itself without college."
branch_of(hh, p::Params, d) = !p.college ? hh : d == 1 ? hh.br1 : hh.br0

"The entry of education branch `d` in a pair `π` of per-branch arrays (`location_probs` under college): `π` itself without college."
branch_probs(π, p::Params, d) = !p.college ? π : d == 1 ? π[1] : π[2]

"""
    cell_branches(hh, p, g, l, zi)
    cell_branches(hh, πT, πO, p, g, l, zi)

The education branches of the cell (gender g, birth location l, ability node zi) as named
tuples `(; h, w, d)`: household arrays, weight P(d | g, l, z) and label. The second form adds
the branch's work-location probabilities `πT`, `πO`. Without college there is one branch of
weight 1 holding `hh` (and `πT`, `πO` as passed), so a loop over it is the baseline loop.
"""
@inline cell_branches(hh, p::Params, g, l, zi) = p.college ?
    ((; h = hh.br1, w = hh.P1[g, l, zi], d = 1), (; h = hh.br0, w = 1 - hh.P1[g, l, zi], d = 0)) :
    ((; h = hh, w = 1.0, d = 1),)
@inline cell_branches(hh, πT, πO, p::Params, g, l, zi) = p.college ?
    ((; h = hh.br1, πT = πT[1], πO = πO[1], w = hh.P1[g, l, zi], d = 1),
     (; h = hh.br0, πT = πT[2], πO = πO[2], w = 1 - hh.P1[g, l, zi], d = 0)) :
    ((; h = hh, πT, πO, w = 1.0, d = 1),)

"The credential premium e^{ρc} of a BA (d = 1) in every non-teaching occupation; 1 otherwise and without college."
@inline credential(p::Params, d) = p.college && d == 1 ? exp(p.ρc) : 1.0

"""
Multiplier of a teacher's pay κ h^γ in education branch `d` at time share `s`: the salary
schedule (s/sBA)^ϑs for BA teachers, e^{−ρT0} without a BA; without college 1, or the schedule
(s/sBA)^{schedT − γφ} when `p.schedT` is set (2026-10-19, `teach_schedule`).
"""
@inline pay_factor(s, p::Params, d) =
    !p.college ? (isnan(p.schedT) ? 1.0 : (s / p.sBA)^teach_schedule(p)[2]) :
    d == 1 ? (iszero(p.ϑs) ? 1.0 : (s / p.sBA)^p.ϑs) : exp(-p.ρT0)

"`pay_factor` at each time share of the vector `sT`."
teacher_pay_factor(sT, p::Params, d) = [pay_factor(s, p, d) for s in sT]

"""
    teacher_logpay_shift(hb, g, l, zi, luT, u, p, d)

Log of the pay factor at the log shocks `u` for the teachers of cell (g, l, zi) of branch `hb`
(the time-share spline over the log-shock grid `luT`, as `micro_sample` splines log h): 0 without
college (or (schedT − γφ)(log s − log sBA) with the schedule `p.schedT`, 2026-10-19), −ρT0 without a BA,
ϑs (log s − log sBA) for BA teachers.
"""
function teacher_logpay_shift(hb, g, l, zi, luT, u, p::Params, d)
    if !p.college
        isnan(p.schedT) && return 0.0
        ls = Spline1D(luT, log.(hb.sT[g, l, zi, :]))
        return teach_schedule(p)[2] .* (ls(u) .- log(p.sBA))
    end
    d == 0 && return -p.ρT0
    iszero(p.ϑs) && return 0.0
    ls = Spline1D(luT, log.(hb.sT[g, l, zi, :]))
    return p.ϑs .* (ls(u) .- log(p.sBA))
end

"""
Teacher effectiveness h_T z^{−(1−ωz)α} of the teachers' human capital `h` at ability node `zi`:
what a teacher contributes to H̃_T (hence Q) and to the CFR experiment. ωz = 1 (the default)
leaves h_T; ωz = 0 makes effectiveness independent of academic ability z. Pay stays on h_T.
"""
effective_h(h, zi, p::Params, gr) = p.ωz == 1.0 ? h : h .* gr.z[zi]^(-(1 - p.ωz) * p.α)

"""
    with_education(nd, p, d, kind) -> Node or NestedNode

The node of education branch `d` (§3c) of the baseline node `nd` of the same problem, `kind` being
`:teach` or `:nonteach`. Income, goods costs and the time-share bounds change as follows (the
location problem, the vectors over destinations and Λ do not, except the BA amenity `Bd` below):
- d = 1: non-teaching income × e^{ρc}; teaching income × sBA^{−ϑs}, which makes the pay
  κ h^γ (s/sBA)^ϑs = (κ h^γ sBA^{−ϑs}) s^ϑs with the extra income elasticity ϑs in s (`Node.ϑs`);
  tuition Fc C̄ added to every goods cost; s ≥ sBA (teaching: max(sT_min, sBA));
- d = 0: teaching income × e^{−ρT0}, a floor sT0 (none by default: no licensing) and s ≤ sBA, a cap
  (none if sBA = 0); the extra cost cT0 of teaching without a BA is in the value, not the node
  (`solve_household_college`).
Only the node's income vector `c`, `m`, `smin`, `smax`, `ϑs` are touched, so the same function serves the
`NestedNode` of the residence/workplace split (income scales on every destination, tuition on every
pair): the education block does not know where the household lives or works. The one exception is
the residence-level additive term of the node, which the BA branch's extra amenity `Bd` raises
(`education_amenity`): within the branch it enters location choice only, but it also raises the BA
node's inclusive value and so the college choice (see `education_amenity`).
"""
function with_education(nd::Node, p::Params, d, kind::Symbol)
    scale, m, smin, smax, ϑs = education_terms(nd.m, p, d, kind)
    Λ = education_amenity(nd.Λ, p, d)
    return Node(nd.c .* scale, nd.γb, nd.τe, m, Λ, nd.l, smin, smax, ϑs,
                nd.pk < 1.0 ? education_amenity(nd.Λ0, p, d) : Λ, nd.pk)       # §3d: the childless have the BA amenity too
end

function with_education(nd::NestedNode, p::Params, d, kind::Symbol)
    scale, m, smin, smax, ϑs = education_terms(nd.m, p, d, kind)
    Λ = education_amenity(nd.Λ, p, d)
    return NestedNode(nd.c .* scale, nd.γb, nd.τe, m, Λ, nd.bw, nd.lP, nd.l,
                      smin, smax, ϑs, nd.pk < 1.0 ? education_amenity(nd.Λ0, p, d) : Λ, nd.pk)
end

"""
    education_amenity(Λ, p, d)

The residence-level additive term of a node in education branch `d`: `Λ` (altruism plus the ability
taste `Bz log z`) plus the extra amenity `Params.Bd` of each residence for BA holders (d = 1), in
utils. It is a taste for where to live, so it is in the location values of the node and nowhere else:
not in consumption, taxes or Ξ, so Proposition 2′ and the X_O* collapse hold, and not in `hh.Λ` (the
altruism term, `dLambda_z`, `school_quality_wtp`). Through the location values it raises the inclusive
value of every BA node, hence W_T, W_O and EV₁, so P(BA) moves with it: a uniform shift of `Bd` is
exactly a shift of the BA cost cC, so `Bd` is normalised to (0, ζ_BA) and ζ_BA is read jointly with cC.
With `Bd` zero, or for d = 0, it returns `Λ` itself.
"""
education_amenity(Λ, p::Params, d) = d == 1 && !all(iszero, p.Bd) ? Λ .+ p.Bd : Λ

"""
    education_terms(m, p, d, kind) -> (scale, m′, smin, smax, ϑs)

What branch `d` changes at a node whose goods costs are `m`: the factor on its income vector, its goods costs,
its time-share bounds and its extra income elasticity in s (`with_education`).
"""
function education_terms(m, p::Params, d, kind::Symbol)
    d in (0, 1) || throw(ArgumentError("college is on: pass the education branch d ∈ {0, 1}, got d = $d"))
    teach = kind === :teach
    scale = d == 1 ? (teach ? (iszero(p.ϑs) ? 1.0 : p.sBA^(-p.ϑs)) : credential(p, 1)) :
                     (teach ? exp(-p.ρT0) : 1.0)
    return (scale, d == 1 ? m .+ p.Fc * p.Cbar : m,
            d == 1 ? (teach ? max(p.sT_min, p.sBA) : p.sBA) : (teach ? p.sT0 : 0.0),
            d == 1 || iszero(p.sBA) ? Inf : p.sBA, d == 1 && teach ? p.ϑs : 0.0)
end

"""
    college_values(hh, gr) -> (EV1, EV0)

EV_d[g, l, zi] = ∫_teachers W_T + ∫_non-teachers W_O of branch d at its own choice maps, the
expected value E_ϵ[max(W_T, W_O)] of the occupational problem given education d. The values are
those of the choice maps (`cm.WTv`, `cm.WOv`): `choice_maps` caps the values at the nodes with no
feasible destination (`cap_infeasible`), and the expectation of the maximum is taken at the values
the choice itself uses; a cell without such a node has the raw values, as the baseline does. Under
the split W is the nested inclusive value (§3b), so nothing else changes.
"""
function college_values(hh, gr)
    EV1 = zeros(2, gr.L, gr.Nz); EV0 = zeros(2, gr.L, gr.Nz)
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        for (EV, hb) in ((EV1, hh.br1), (EV0, hh.br0))
            cm = choice_maps(hb, g, l, zi, gr)
            EV[g, l, zi] = sum(integrate_choice(cm.WTv, cm.WOv, cm, gr))
        end
    end
    return EV1, EV0
end

"""
    college_choice!(hh, p, gr) -> hh

The college choice P₁[g, l, zi] = Λ([EV₁ − EV₀ − cC − cCz log z]/σd), in place, from the branch values.
"""
function college_choice!(hh, p::Params, gr)
    EV1, EV0 = college_values(hh, gr)
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz
        Δ = EV1[g, l, zi] - EV0[g, l, zi] - p.cC - p.cCz * log(gr.z[zi])
        hh.P1[g, l, zi] = 1 / (1 + exp(-Δ / p.σd))
    end
    return hh
end

"""
    solve_household_college(p, gr; tol, maxit, verbose, print_every, hh0, Λ0)

`solve_household` with the college choice: each sweep solves every node of both education
branches given Λ, computes P₁ from the new values (`college_choice!`) and then Λ with the
mixed warm glow (`child_Ef`), until the values converge. Returns `(; college = true, br1, br0,
P1, Λ, s_resid, s_passes, s_relaxed)` (§3c); `hh0`, a solution of the same layout, warm-starts
the policies and P₁. `s_resid`, `s_passes` and `s_relaxed` are over both branches.
"""
function solve_household_college(p::Params, gr; tol = 1e-6, maxit = 500,
                                 verbose = false, print_every = 25,
                                 hh0 = nothing, Λ0 = hh0 === nothing ? nothing : hh0.Λ)
    (; Nz, L, nϵT, nXO) = gr
    szT, szO = (2, L, Nz, nϵT), (2, L, Nz, nXO)
    Λ  = (Λ0 !== nothing && size(Λ0) == (Nz, L)) ? copy(Λ0) : zeros(Nz, L)
    branch() = (; eT = zeros(szT), sT = fill(s_of(p.T, p), szT), WT = zeros(szT), hT = zeros(szT),
                  eO = zeros(szO), sO = fill(s_of(gr.nonteach[1], p), szO), WO = zeros(szO),
                  capO = zeros(szO), Λ, infeasT = falses(szT), infeasO = falses(szO))
    br1, br0 = branch(), branch()
    P1 = fill(0.5, 2, L, Nz)
    if hh0 !== nothing && haskey(hh0, :br1) && size(hh0.br1.eT) == szT && size(hh0.br1.eO) == szO
        for (b, b0) in ((br1, hh0.br1), (br0, hh0.br0))
            b.eT .= b0.eT; b.sT .= b0.sT; b.eO .= b0.eO; b.sO .= b0.sO
        end
        P1 .= hh0.P1
    end
    hh = (; college = true, br1, br0, P1, Λ)

    # one sweep over the nodes of branch d; the worst (e, s) residual, passes, relaxed-step nodes
    function solve_branch!(b, d)
        s_resid, s_passes, n_relaxed = 0.0, 0, 0
        cred = credential(p, d)
        cT0  = d == 0 ? p.cT0 : 0.0                  # the extra cost of teaching without a BA
        for g in 1:2, l in 1:L, zi in 1:Nz
            for k in 1:nϵT
                ϵ  = gr.ϵTgrid[k]
                nd = teach_node(ϵ, zi, l, g, Λ, p, gr; d)
                e, V̄, s, r, n, rx = solve_node(nd, p, log(b.eT[g, l, zi, k]), b.sT[g, l, zi, k])
                b.eT[g, l, zi, k] = e
                b.sT[g, l, zi, k] = s
                b.WT[g, l, zi, k] = log(1 - s) + V̄ - p.cT - p.cTz * log(gr.z[zi]) - cT0
                b.hT[g, l, zi, k] = teacher_h(ϵ, zi, l, s, e, p, gr)
                b.infeasT[g, l, zi, k] = all(≤(C_FLOOR), consumption(nd, s, e, p))
                s_resid, s_passes, n_relaxed = max(s_resid, r), max(s_passes, n), n_relaxed + rx
            end
            for k in 1:nXO
                x  = gr.XOgrid[k, g]
                nd = nonteach_node(x, zi, l, g, Λ, p, gr; d)
                e, V̄, s, r, n, rx = solve_node(nd, p, log(b.eO[g, l, zi, k]), b.sO[g, l, zi, k])
                b.eO[g, l, zi, k]   = e
                b.sO[g, l, zi, k]   = s
                b.WO[g, l, zi, k]   = log(1 - s) + V̄
                b.capO[g, l, zi, k] = nonteach_cap(x, zi, l, s, e, p, gr) * cred
                b.infeasO[g, l, zi, k] = all(≤(C_FLOOR), consumption(nd, s, e, p))
                s_resid, s_passes, n_relaxed = max(s_resid, r), max(s_passes, n), n_relaxed + rx
            end
        end
        return s_resid, s_passes, n_relaxed
    end

    s_resid, s_passes, s_relaxed = 0.0, 0, 0
    t0 = time()
    for it in 1:maxit
        old = (copy(br1.WT), copy(br1.WO), copy(br0.WT), copy(br0.WO))
        s_resid, s_passes, n_relaxed = 0.0, 0, 0
        for (b, d) in ((br1, 1), (br0, 0))
            r, n, rx = solve_branch!(b, d)
            s_resid, s_passes, n_relaxed = max(s_resid, r), max(s_passes, n), n_relaxed + rx
        end
        s_relaxed = max(s_relaxed, n_relaxed)
        college_choice!(hh, p, gr)
        Λ .= compute_lambda(hh, p, gr)
        err = max(maximum(abs, br1.WT .- old[1]), maximum(abs, br1.WO .- old[2]),
                  maximum(abs, br0.WT .- old[3]), maximum(abs, br0.WO .- old[4]))
        if verbose && (it == 1 || it % print_every == 0 || err < tol)
            @printf("    HH %4d/%d  err=%.3e  BA share %.4f  elapsed=%.1fs\n", it, maxit, err,
                    sum(P1) / length(P1), time() - t0)
            flush(stdout)
        end
        err < tol && break
    end
    return (; hh..., s_resid, s_passes, s_relaxed)
end

"""
    college_child_Ef(l, g, zi, hh, p, gr)

E[f(h′)] of a child born in l with gender g and ability node zi, mixing the education branches:
P₁ E[f | d = 1] + (1 − P₁) E[f | d = 0], each integrated over its own occupation choice, with
the non-teachers' human capital net of the credential premium (`child_Ef`).
"""
function college_child_Ef(l, g, zi, hh, p::Params, gr)
    P1 = hh.P1[g, l, zi]
    function Ef(hb, d)
        cred = credential(p, d)
        cm = choice_maps(hb, g, l, zi, gr)
        fT = warm.(hb.hT[g, l, zi, :], p.ψ, p.href)
        fO = [Ewarm_nonteach(hb.capO[g, l, zi, k] / cred, gr.XOgrid[k, g], g, p, gr) for k in 1:gr.nXO]
        return sum(integrate_choice(fT, fO, cm, gr))
    end
    return P1 * Ef(hh.br1, 1) + (1 - P1) * Ef(hh.br0, 0)
end

"The pairs ((πT₁, πT₀), (πO₁, πO₀)) of work-location probabilities of the two education branches (`location_probs`; `who` as there)."
function college_location_probs(hh, p::Params, gr; who::Symbol = :all)
    pr1 = location_probs(hh.br1, p, gr; d = 1, who)
    pr0 = location_probs(hh.br0, p, gr; d = 0, who)
    return (pr1[1], pr0[1]), (pr1[2], pr0[2])
end

"Teaching share among the young born in `l` (both genders), mixing the education branches (`teaching_share_endo`)."
function college_teaching_share(hh, l, Φ, gr, p::Params)
    mass = sum(@view Φ[:, l])
    return sum(Φ[zi, l] / mass / 2 * br.w * integrate_teach(ones(gr.nϵT), choice_maps(br.h, g, l, zi, gr), gr)
               for zi in 1:gr.Nz, g in 1:2 for br in cell_branches(hh, p, g, l, zi))
end

# -----------------------------------------------------------------------------
# 4. Occupation margin: teach iff W_T(ϵ_T) > W_O(X_O*)
#
#   teach_wt(ϵ_T)    = P(teach | ϵ_T)  = F_{X_O*}( W_O^{-1}(W_T(ϵ_T)) )
#   nonteach_wt(X_O*) = P(don't | X_O*) = F_{ϵ_T}( W_T^{-1}(W_O(X_O*)) )
#
# Values are splined in the log shocks, where they are close to linear; in levels a
# cubic spline between widely spaced tail nodes misses the values by several units.
# -----------------------------------------------------------------------------
struct ChoiceMaps{G}
    splWT::Spline1D;    splWO::Spline1D      # W_T(log ϵ_T), W_O(log X_O*)
    splWTinv::Spline1D; splWOinv::Spline1D   # their inverses, returning logs
    WTv::Vector{Float64}; WOv::Vector{Float64}
    g::Int
    gr::G
end

# The values should be strictly increasing; `make_increasing` only patches
# round-off. Warn (at most 5 times) if a value vector is materially non-monotone.
const NONMONO_WARNINGS = Ref(0)
function warn_nonmonotone(v::AbstractVector, name::String)
    drop  = maximum(v[i-1] - v[i] for i in 2:length(v); init = 0.0)
    scale = max(abs(v[end] - v[1]), 1e-12)
    if drop > 1e-8 * scale && NONMONO_WARNINGS[] < 5
        NONMONO_WARNINGS[] += 1
        suffix = NONMONO_WARNINGS[] == 5 ? " (suppressing further warnings)" : ""
        @printf("  ⚠ %s non-monotone (max drop %.3e vs range %.3e): occupation threshold may be off%s\n",
                name, drop, scale, suffix)
        flush(stdout)
    end
end

function choice_maps(WTv, WOv, g, gr)
    warn_nonmonotone(WTv, "W_T(ϵ_T)")
    warn_nonmonotone(WOv, "W_O(X_O*)")
    return ChoiceMaps(Spline1D(gr.lϵT, WTv), Spline1D(gr.lXO[:, g], WOv),
                      Spline1D(make_increasing(WTv), gr.lϵT; k = 1),
                      Spline1D(make_increasing(WOv), gr.lXO[:, g]; k = 1),
                      WTv, WOv, g, gr)
end

"""
    cap_infeasible(WTv, WOv, infT, infO) -> (WTv′, WOv′)

Values at the nodes of a cell with no feasible destination (every destination at the consumption
floor: `hh.infeasT`, `hh.infeasO`) sit on `smooth_log`'s linear branch. They are not utilities but
a penalty that grows with the shortfall in consumption (hundreds of utils within a few nodes of the
feasible ones), and a cubic spline through them rings across the jump to the first feasible node.
At the m2_rp075 point the spline of W_T peaked at −5.7 (min W_O −6.3), and after a step of 0.001 in
δκ at +3.7, above every W_O: teach weights of one on a fifth of the interval where nobody teaches,
which moved the male teaching share by 1.6%, the teacher 90/10 from 2.84 to 2.90 and the
population share at the floor from 5e−7 to 2e−4. The values at those nodes are replaced by one
level for all the infeasible nodes of a branch,
    c_T = min(a_O − 1, a_T) for the teacher nodes,    c_O = min(a_T − 1, a_O) for the non-teacher nodes,
with a_T, a_O the lowest values among the feasible nodes of each branch (∞ for a branch without
one): the lower of one util below the other branch's lowest feasible value, so that nobody takes
the infeasible branch there against a feasible alternative (the teach weight is exactly zero), and
of the branch's own lowest feasible value, so that the values stay nondecreasing in the shock and the
jump at the feasibility boundary is max(0, a_T − a_O + 1), not the hundreds of the extrapolation. A
cell with no feasible node on either branch is left alone. Only nodes with no feasible destination
change, so a cell without any, as at the bundled baseline, is untouched bit for bit.
"""
function cap_infeasible(WTv, WOv, infT, infO)
    aT = minimum(view(WTv, .!infT); init = Inf)
    aO = minimum(view(WOv, .!infO); init = Inf)
    cT, cO = min(aO - 1, aT), min(aT - 1, aO)
    WT, WO = copy(WTv), copy(WOv)
    isfinite(cT) && (WT[infT] .= cT)
    isfinite(cO) && (WO[infO] .= cO)
    return WT, WO
end

"""
    choice_maps(hh, g, l, zi, gr) -> ChoiceMaps (or MixedMaps under CCS)

The choice maps of cell (g, l, zi) of a household solution. At the nodes with no feasible destination
(`hh.infeasT`, `hh.infeasO`) the values are capped first (`cap_infeasible`); a cell without such a
node, and an `hh` without the flags, use the raw values. `choice_maps(WTv, WOv, g, gr)` builds the
maps from value vectors as they are. When `gr` carries the teacher-draw laws of constant class size (`ccs_on`, §3e)
the cell form returns a `MixedMaps` over the draws of R_l, l being the cell's school location (its birth location).
"""
function choice_maps(hh, g, l, zi, gr)
    WTv, WOv = hh.WT[g, l, zi, :], hh.WO[g, l, zi, :]
    if haskey(hh, :infeasT)
        infT, infO = hh.infeasT[g, l, zi, :], hh.infeasO[g, l, zi, :]
        (any(infT) || any(infO)) && ((WTv, WOv) = cap_infeasible(WTv, WOv, infT, infO))
    end
    ccs_on(gr) && return mixed_maps(WTv, WOv, g, l, gr)       # §3e
    return choice_maps(WTv, WOv, g, gr)
end

"P(teach | log ϵ_T = u)."
@inline function teach_wt_log(cm::ChoiceMaps, u)
    wT = cm.splWT(u)
    wT <= cm.WOv[1]   && return 0.0
    wT >= cm.WOv[end] && return 1.0
    return Fxo(cm.gr, cm.g, exp(cm.splWOinv(wT)))
end

"P(don't teach | log X_O* = v)."
@inline function nonteach_wt_log(cm::ChoiceMaps, v)
    wO = cm.splWO(v)
    wO <= cm.WTv[1]   && return 0.0
    wO >= cm.WTv[end] && return 1.0
    return cdf(cm.gr.dT, exp(cm.splWTinv(wO)))
end

@inline teach_wt(cm::ChoiceMaps, ϵT) = teach_wt_log(cm, log(ϵT))
@inline nonteach_wt(cm::ChoiceMaps, XO) = nonteach_wt_log(cm, log(XO))

# Integrals of a grid quantity over the teachers / non-teachers of one cell:
# quadrature in the log shock (density f(eᵘ)eᵘ) over the grid, plus the shock mass
# outside the grid at the endpoint values (constant extrapolation; the residual is
# O(1 − q_hi)). The grid quantity is interpolated by a spline of degree `k` in the log
# shock: cubic by default; `k = 1` for a quantity that must stay nonnegative (a cubic
# spline through nonnegative values can undershoot zero between nodes). A linear spline of
# a 0/1 indicator is a ramp one grid interval wide, which an adaptive rule can step over
# (it missed 10% of a tail mass of 6e−6 and half of a mass of 0.04), so `k = 1` splits the
# quadrature at the grid nodes.
function integrate_teach(q, cm::ChoiceMaps, gr; k = 3)
    (; dT, ϵTgrid, lϵT) = gr
    lo, hi = ϵTgrid[1], ϵTgrid[end]
    spl = Spline1D(lϵT, q; k)
    f(u) = (w = teach_wt_log(cm, u); w <= 0.0 ? 0.0 : (ϵ = exp(u); pdf(dT, ϵ) * ϵ * spl(u) * w))
    I = k == 1 ? quadgk(f, lϵT...)[1] : quadgk(f, lϵT[1], lϵT[end])[1]
    return I +
           cdf(dT, lo)  * q[1]   * teach_wt(cm, lo) +
           ccdf(dT, hi) * q[end] * teach_wt(cm, hi)
end

function integrate_nonteach(q, cm::ChoiceMaps, gr; k = 3)
    g = cm.g
    lx = view(gr.lXO, :, g)
    lo, hi = gr.XOgrid[1, g], gr.XOgrid[end, g]
    spl = Spline1D(lx, q; k)
    f(v) = (w = nonteach_wt_log(cm, v); w <= 0.0 ? 0.0 : (x = exp(v); fxo(gr, g, x) * x * spl(v) * w))
    I = k == 1 ? quadgk(f, lx...)[1] : quadgk(f, lx[1], lx[end])[1]
    return I +
           Fxo(gr, g, lo)       * q[1]   * nonteach_wt(cm, lo) +
           (1 - Fxo(gr, g, hi)) * q[end] * nonteach_wt(cm, hi)
end

"(∫ over teachers of qT(ϵ_T), ∫ over non-teachers of qO(X_O*)) for one cell; `k` as above."
integrate_choice(qT, qO, cm::ChoiceMaps, gr; k = 3) =
    (integrate_teach(qT, cm, gr; k), integrate_nonteach(qO, cm, gr; k))

# -----------------------------------------------------------------------------
# 3e. Constant class size (`Params.ccs`, 2026-10-26)
#
# (Placed after §4 because it builds on the choice maps.) Classes are equal-sized within a location, N_l = S_l/n_T,l,
# and each student of l is assigned at random to K_d teachers (`ccs_draws`), so school quality is
#     ln q = ln Q̄_l + ln r,   ln Q̄_l = β μ_l − σ ln N_l,   μ_l = E[ln h_eff | teacher works in l],
#     ln r = (β/K_d) Σ_t (ln h_eff,t − μ_l),   κ_m(ln r) = β^m κ_m(ln h_eff | l)/K_d^{m−1}   (`draw_cumulants`),
# and q is a state of the young, seen with (z, ϵ) before (s, e, i). With `ccs_rho` = ρ < 1 the students also receive
# the CFR noise ν of teacher effectiveness, which adds κ₂(ν) = κ₂(ln h_eff | l)(1/ρ² − 1) to κ₂ before the scaling.
# Proposition 5 (effective shocks): q enters the young's problem only through h = q (z ϵ_i)^α s^φ e^η, so the problem at
# (g, l, z, ϵ, q = Q̄ r) is the problem at Q̄ with ϵ̃_T = r^{1/α} ϵ_T and X̃_O = r X_O: W_T(ϵ_T; q) = W̄_T(r^{1/α} ϵ_T),
# W_O(X_O*; q) = W̄_O(r X_O*), and likewise every policy, h_T and the income capacity. The nodes are therefore solved once,
# at gr.Q = Q̄ on the effective-shock grids (widened by a pad, `build_grids`), with no change to any node function, and a
# draw x = ln r changes only the laws of the effective shocks (`ccs_law`):
#     ϵ̃_T ~ LogNormal(−σϵ²/2 + x/α, σϵ),    X̃_{O,i} ~ LogNormal(log Θ_{i,g} − ασϵ²/2 + x, ασϵ),
# under which every integral of §4 is exact for that draw: the occupation weights, the densities, the end masses and
# `argmax_mean`, since P(i* = i | X̃* = x̃, r) ∝ f_i(x̃/r)/F_i(x̃/r) is the ratio of the shifted laws. The law of x in
# location l is represented by R_l, the n_q-point Gauss rule matched to the moments of ln r (`gauss_from_moments`;
# n_q = 1 is the mean field), and the integral over a cell is the R_l-mixture of the draw integrals (`MixedMaps`,
# the cell form of `choice_maps` under CCS). Integrands that depend on the draw itself, through the law (E over Θ_{i*}
# given X̃_O*: `argmax_mean`, `EΘpow`, `Eϵ_nonteach`) or through the true shock ϵ = ϵ̃ e^{−x/α}, go through `over_draws`.
# Safety convention: under CCS the top-level `gr.dT` and `gr.dO` are `nothing`, so a consumer that has not been taught
# the draws (a ChoiceMaps built from the top-level grid, a direct `pdf(gr.dT, …)`) throws instead of integrating at r = 1.
# -----------------------------------------------------------------------------
"""
    ccs_on(gr) -> Bool

Whether the grid `gr` carries the teacher-draw laws of constant class size (`build_grids` under `Params.ccs`, §3e). The
switch is read from the grid, not from `Params`, because `choice_maps` and the integrals have no `p`; a draw's law
(`ccs_law`) has no `ccs` field and reads false.
"""
ccs_on(gr) = hasproperty(gr, :ccs) && gr.ccs !== nothing

"ln Q̄ by location for `build_grids` under CCS: `ccs.lQT`, or without a state the efficient σ ln(2H̃_T/M) (as the baseline's Q_T)."
ccs_lQT(p::Params, ccs, H̃T, M) = ccs === nothing ? p.σ .* log.(2 .* max.(H̃T, M_FLOOR) ./ max.(M, M_FLOOR)) :
                                                    collect(Float64, ccs.lQT)

"""
    moments_to_cumulants(m) -> κ

Cumulants κ_1..κ_n of the distribution with moments `m` of orders 0..n (`m[j + 1]` the order j, `m[1]` the mass; raw or
central), by κ_n = μ_n − Σ_{k=1}^{n−1} C(n−1, k−1) κ_k μ_{n−k} with μ = m/m[1]. From order 2 on they do not depend on
the centre.
"""
function moments_to_cumulants(m::AbstractVector)
    μ = m ./ m[1]
    n = length(m) - 1
    κ = zeros(n)
    for j in 1:n
        s = μ[j + 1]
        for k in 1:j-1
            s -= binomial(j - 1, k - 1) * κ[k] * μ[j - k + 1]
        end
        κ[j] = s
    end
    return κ
end

"Moments of orders 0..n (mass 1) of the distribution with cumulants κ_1..κ_n: μ_n = Σ_{k=1}^{n} C(n−1, k−1) κ_k μ_{n−k}."
function cumulants_to_moments(κ::AbstractVector)
    n = length(κ)
    μ = zeros(n + 1)
    μ[1] = 1.0
    for j in 1:n
        s = 0.0
        for k in 1:j
            s += binomial(j - 1, k - 1) * κ[k] * μ[j - k + 1]
        end
        μ[j + 1] = s
    end
    return μ
end

"""
    draw_cumulants(c, β, K_d, ρ = 1.0) -> κ
    draw_moments(c, β, K_d, ρ = 1.0) -> m

The law of the teacher draw ln r of one location (§3e) from the central moments `c` of ln h_eff among its teachers
(orders 0..n in `c[1:n+1]`; `ccs_teacher_moments`): their cumulants (`moments_to_cumulants`), κ₂ raised by the CFR
noise κ₂(ν) = κ₂(1/ρ² − 1) (that is, κ₂/ρ²; ν Gaussian, so the higher cumulants are unchanged), then scaled to the mean
of K_d independent draws times β, κ_m(ln r) = β^m κ_m/K_d^{m−1}, with κ_1 = 0 (ln r is centred at μ_l, which is in Q̄).
`draw_moments` returns the moments of ln r, orders 0..n (`cumulants_to_moments`): central, with `m[1]` = 1 and `m[2]` = 0.
"""
function draw_cumulants(c::AbstractVector, β, K_d, ρ = 1.0)
    κ = moments_to_cumulants(c)
    κ[1] = 0.0
    length(κ) ≥ 2 && (κ[2] /= ρ^2)
    for m in eachindex(κ)
        κ[m] *= β^m / Float64(K_d)^(m - 1)
    end
    return κ
end

draw_moments(c::AbstractVector, β, K_d, ρ = 1.0) = cumulants_to_moments(draw_cumulants(c, β, K_d, ρ))

"""
    gauss_from_moments(m, n; sd_tol = 1e-12, diag_tol = 1e-7) -> (x, w)

The n-point Gauss rule (nodes `x`, weights `w` summing to 1) of the distribution with moments `m` of orders 0..K
(`m[j + 1]` the order j, K ≥ 2n − 1; raw or central, any mass), by Golub–Welsch: the moments are standardized
(s_j = E[(X − μ)^j]/sd^j), H = [s_{i+j−2}] is the (n+1)×(n+1) Hankel matrix with upper Cholesky factor R, and the
Jacobi matrix has the diagonal α_j = r_{j,j+1}/r_{j,j} − r_{j−1,j}/r_{j−1,j−1} and the off-diagonal β_j = r_{j+1,j+1}/r_{j,j};
its eigenvalues are the standardized nodes and the squared first components of its eigenvectors the weights. The rule
matches the moments of orders 0..2n − 1. Only the first n rows of R enter the recurrence, and they need neither
s_{2n} nor the last pivot r_{n+1,n+1} (the Cholesky of the leading n×n block plus one forward substitution for the
last column), so a measure with exactly n support points (r_{n+1,n+1} = 0) still gets its exact n-point rule. If the
leading block is not positive definite, or a pivot r_{j,j} falls below `diag_tol` (the measure has fewer than n points
up to rounding), a node is dropped and the rule rebuilt. n = 1 and a standard deviation below `sd_tol` (the test is on
the variable whose moments are given: for the teacher draw, ln r after the β scaling, so a β sweep does not switch
the order of the rule) return the one-point rule at the mean, ([μ], [1.0]).
"""
function gauss_from_moments(m::AbstractVector, n::Int; sd_tol = 1e-12, diag_tol = 1e-7)
    K = length(m) - 1
    μ = m[2] / m[1]
    # central moments about μ, normalised to mass 1
    cm = [sum(binomial(j, i) * (m[i + 1] / m[1]) * (-μ)^(j - i) for i in 0:j) for j in 0:K]
    sd = K ≥ 2 && cm[3] > 0.0 ? sqrt(cm[3]) : 0.0
    (n ≤ 1 || !(sd ≥ sd_tol)) && return ([μ], [1.0])
    s = [cm[j + 1] / sd^j for j in 0:K]
    for nn in min(n, (K + 1) ÷ 2):-1:2
        H = [s[i + j - 1] for i in 1:nn, j in 1:nn]          # orders 0..2nn−2
        F = cholesky(Symmetric(H); check = false)
        issuccess(F) || continue
        R = F.U
        all(R[j, j] > diag_tol for j in 1:nn) || continue
        r = zeros(nn)                                         # the last column r_{j,nn+1}: orders nn−1..2nn−1
        for j in 1:nn
            r[j] = (s[j + nn] - sum(R[k, j] * r[k] for k in 1:j-1; init = 0.0)) / R[j, j]
        end
        a = [(j < nn ? R[j, j + 1] : r[nn]) / R[j, j] - (j > 1 ? R[j - 1, j] / R[j - 1, j - 1] : 0.0) for j in 1:nn]
        b = [R[j + 1, j + 1] / R[j, j] for j in 1:nn-1]
        E = eigen(SymTridiagonal(a, b))
        w = E.vectors[1, :] .^ 2
        (all(isfinite, E.values) && all(isfinite, w) && all(>(0.0), w)) || continue
        return (μ .+ sd .* E.values, w ./ sum(w))
    end
    return ([μ], [1.0])
end

"""
    ccs_rule(c, p) -> (; x, w)

R_l, the rule of the teacher draw of a location with teacher central moments `c` (orders 0..2n_q of ln h_eff; `nothing`:
none): the `p.ccs_nq`-point Gauss rule (`gauss_from_moments`) of the moments of ln r (`draw_moments` with `p.β`,
`p.ccs_draws`, `p.ccs_rho`). Nodes `x` are values of x = ln r; a degenerate draw (β = 0, homogeneous teachers, n_q = 1) is
([0.0], [1.0]).
"""
function ccs_rule(c, p::Params)
    c === nothing && return (; x = [0.0], w = [1.0])
    x, w = gauss_from_moments(draw_moments(c, p.β, p.ccs_draws, p.ccs_rho), p.ccs_nq)
    return (; x, w)
end

"""
    ccs_law(g0, x, p) -> NamedTuple

The law-only grid object of the teacher draw x = ln r (§3e): the laws of the effective shocks, `dT` = LogNormal(−σϵ²/2 +
x/α, σϵ) and `dO[i, g]` = LogNormal(log Θ_{i,g} − ασϵ²/2 + x, ασϵ), with the (shared, not copied) effective-shock grids of
`g0` (`ϵTgrid`, `lϵT`, `XOgrid`, `lXO`, `nϵT`, `nXO`), `nonteach`, `logΘ` and `x`. It is what the integrals of §4 read as
`gr` (`integrate_teach`, `integrate_nonteach`, `Fxo`, `fxo`, `argmax_mean`) and what the per-draw `ChoiceMaps` of a
`MixedMaps` hold, and it has no Q, z, t or `ccs` field, so passing it where a full grid is needed throws. At x = 0 the
laws are those of the true shocks, built by the same operations as `build_grids` without CCS. (The plan's `ccs_law(gr0,
x)` takes `p` as well: the ϵ_T shift is x/α.)
"""
function ccs_law(g0, x::Real, p::Params)
    I = size(g0.logΘ, 1)
    dT = LogNormal(-p.σϵ^2 / 2 + x / p.α, p.σϵ)
    dO = Matrix{LogNormal{Float64}}(undef, I, 2)
    for g in 1:2, i in g0.nonteach
        dO[i, g] = LogNormal(g0.logΘ[i, g] - p.α * p.σϵ^2 / 2 + x, p.α * p.σϵ)
    end
    return (; dT, dO, ϵTgrid = g0.ϵTgrid, lϵT = g0.lϵT, XOgrid = g0.XOgrid, lXO = g0.lXO, nonteach = g0.nonteach,
              logΘ = g0.logΘ, nϵT = g0.nϵT, nXO = g0.nXO, x = Float64(x))
end

"`n` nodes below and above the increasing vector `v` of log nodes, at the spacing `Δ`, kept exact in the middle: exp of the result."
pad_grid(v, Δ, n) = n == 0 ? copy(v) : vcat(exp.(log(v[1]) .- Δ .* (n:-1:1)), v, exp.(log(v[end]) .+ Δ .* (1:n)))

"""
    ccs_grids(p, g0, ccs, pad) -> gr

The grids of `build_grids` under CCS (§3e) from the base grids `g0`: the effective-shock grids widened at both ends by
`pad` with the base spacing (`pad_grid`; ⌈pad/α/Δ⌉ nodes for ϵ_T, a common ⌈pad/Δ_g⌉ = max over g for X_O*), the rule
R_l of every location from `ccs.c` (`ccs_rule`), the laws of its draws (`ccs_law`), and `dT = dO = nothing`.
"""
function ccs_grids(p::Params, g0, ccs, pad)
    pad ≥ 0.0 || throw(ArgumentError("the CCS grid pad must be nonnegative, got $pad"))
    L = g0.L
    lT0 = log.(g0.ϵTgrid)
    ΔT  = (lT0[end] - lT0[1]) / (length(lT0) - 1)
    nT  = pad > 0.0 ? ceil(Int, pad / p.α / ΔT) : 0
    ϵTgrid = pad_grid(g0.ϵTgrid, ΔT, nT)
    ΔO  = [(log(g0.XOgrid[end, g]) - log(g0.XOgrid[1, g])) / (size(g0.XOgrid, 1) - 1) for g in 1:2]
    nO  = pad > 0.0 ? maximum(ceil(Int, pad / ΔO[g]) for g in 1:2) : 0
    XOgrid = nO == 0 ? copy(g0.XOgrid) : hcat((pad_grid(g0.XOgrid[:, g], ΔO[g], nO) for g in 1:2)...)
    gl = (; ϵTgrid, lϵT = log.(ϵTgrid), XOgrid, lXO = log.(XOgrid), nonteach = g0.nonteach, logΘ = g0.logΘ,
            nϵT = length(ϵTgrid), nXO = size(XOgrid, 1))
    c     = ccs === nothing ? nothing : ccs.c
    rules = [ccs_rule(c === nothing ? nothing : c[l], p) for l in 1:L]
    laws  = [[ccs_law(gl, x, p) for x in r.x] for r in rules]
    lQT   = log.(g0.QT)
    return (; g0.z, g0.Πz, g0.Nz, g0.Q, g0.QT, g0.la, g0.t, g0.I, g0.L, dT = nothing, gl.ϵTgrid, dO = nothing,
              gl.XOgrid, gl.lϵT, gl.lXO, g0.logΘ, g0.nonteach, gl.nϵT, gl.nXO, g0.PT, g0.pk,
              ccs = (; x = [r.x for r in rules], w = [r.w for r in rules], laws, base = ccs_law(gl, 0.0, p),
                       pad = Float64(pad), lQT = ccs === nothing ? lQT : collect(Float64, ccs.lQT),
                       c = c === nothing ? nothing : [copy(v) for v in c], nq = p.ccs_nq, draws = p.ccs_draws, rho = p.ccs_rho))
end

"""
    MixedMaps

The choice maps of one cell under CCS (§3e): the per-draw `ChoiceMaps` `maps[k]`, which share the splines of the cell's
values (the nodes are the same for every draw, Proposition 5) and each hold the law of draw k (`maps[k].gr`, a
`ccs_law`), the weights `w` and nodes `x` (of ln r) of the cell's rule R_l, the gender `g`, and the values `WTv`, `WOv`
(capped as by `choice_maps`). `integrate_teach`, `integrate_nonteach` and `integrate_choice` take it in place of a
`ChoiceMaps` and return Σ_k w_k × the draw-k integral, for a quantity given at the effective-shock nodes that does not
depend on the draw (masses, location probabilities, policies, values, h_T, the income capacity); a quantity that does
goes through `over_draws`. `teach_wt` and `nonteach_wt` have no method for it: P(teach) at a shock depends on whether
the shock is effective or true, so it is taken per draw.
"""
struct MixedMaps{C<:ChoiceMaps}
    maps::Vector{C}
    w::Vector{Float64}
    x::Vector{Float64}
    g::Int
    WTv::Vector{Float64}; WOv::Vector{Float64}
end

"The `MixedMaps` of a cell born in `l` from its (capped) values: the splines once, then one `ChoiceMaps` per law of R_l."
function mixed_maps(WTv, WOv, g, l, gr)
    laws = gr.ccs.laws[l]
    cm1  = choice_maps(WTv, WOv, g, laws[1])
    maps = [k == 1 ? cm1 : ChoiceMaps(cm1.splWT, cm1.splWO, cm1.splWTinv, cm1.splWOinv, cm1.WTv, cm1.WOv, g, laws[k])
            for k in eachindex(laws)]
    return MixedMaps(maps, gr.ccs.w[l], gr.ccs.x[l], g, cm1.WTv, cm1.WOv)
end

"`integrate_teach` over the draws of a `MixedMaps`: Σ_k w_k ∫ (draw k's law); the top-level `gr` is not read."
function integrate_teach(q, cm::MixedMaps, gr; k = 3)
    I = 0.0
    for j in eachindex(cm.w)
        I += cm.w[j] * integrate_teach(q, cm.maps[j], cm.maps[j].gr; k)
    end
    return I
end

"`integrate_nonteach` over the draws of a `MixedMaps` (as `integrate_teach`)."
function integrate_nonteach(q, cm::MixedMaps, gr; k = 3)
    I = 0.0
    for j in eachindex(cm.w)
        I += cm.w[j] * integrate_nonteach(q, cm.maps[j], cm.maps[j].gr; k)
    end
    return I
end

"`integrate_choice` over the draws of a `MixedMaps`."
integrate_choice(qT, qO, cm::MixedMaps, gr; k = 3) =
    (integrate_teach(qT, cm, gr; k), integrate_nonteach(qO, cm, gr; k))

"""
    over_draws(f, cm, gr)

The draw mixture of `f(cm_k, law_k, x_k)` (§3e): Σ_k w_k f(cm.maps[k], law k, x_k) for a `MixedMaps`, where law k
(`cm.maps[k].gr`, a `ccs_law`) is what the integrals and `argmax_mean` / `EΘpow` / `Eϵ_nonteach` must be given as `gr`,
and x_k = ln r_k recovers the true shocks, ϵ_T = ϵ̃_T e^{−x/α} and X_O = X̃_O e^{−x}; `f(cm, gr, 0.0)` for a plain
`ChoiceMaps`, so code written with it is the baseline's at r = 1. `f` may return a number, a tuple or an array (the sum
broadcasts). Use it for any integrand that depends on the law of the draw or on the true shock, and with the do-block
syntax: `over_draws(cm, gr) do cmk, law, x … end`.
"""
over_draws(f, cm::ChoiceMaps, gr) = f(cm, gr, 0.0)
function over_draws(f, cm::MixedMaps, gr)
    acc = cm.w[1] .* f(cm.maps[1], cm.maps[1].gr, cm.x[1])
    for k in 2:length(cm.w)
        acc = acc .+ cm.w[k] .* f(cm.maps[k], cm.maps[k].gr, cm.x[k])
    end
    return acc
end

"SD of ln r under the rule (x, w) of a location (the matched SD for n_q ≥ 2; 0 for the one-point rule)."
rule_sd(x, w) = sqrt(max(sum(w .* x .^ 2) - sum(w .* x)^2, 0.0))

# -----------------------------------------------------------------------------
# 5. Altruism term
#
#   Λ(z, l') = λ Σ_{z'} Πz(z, z') · ½ Σ_{g'} E[f(h') | l', g', z'],
# where the expectation integrates the child's occupation choice. On the
# non-teaching side h_O = cap/Θ_{i*} and only the law of Θ_{i*} given X_O* is known,
# so E[f(h_O) | X_O*] is taken over Θ_{i*} in the space the kernel needs.
# -----------------------------------------------------------------------------
function Ewarm_nonteach(cap, x, g, p::Params, gr)
    is_log_kernel(p.ψ) && return log(cap / p.href) - argmax_mean(identity, gr, g, x)
    u = (1 - p.ψ) * log(cap / p.href) + log(EΘpow(gr, g, x, p.ψ))
    return expm1(u) / (1 - p.ψ)
end

"""
    child_Ef(l, g, zi, hh, p, gr)

E[f(h′)] of a child born in l with gender g and ability node zi, over the child's occupation choice. Under constant
class size (§3e) it is also the expectation over the teacher draw, Σ_k w_k E[f(h′) | r_k] with R_l (`over_draws`): the
non-teachers' E[f(h_O) | X̃_O*] is taken under each draw's law of Θ_{i*} (h_O = capO/Θ_{i*} at the effective node holds
for every draw), so parents value a location by its whole teacher distribution. Without CCS it is one pass at `gr`.
"""
function child_Ef(l, g, zi, hh, p::Params, gr)
    p.college && return college_child_Ef(l, g, zi, hh, p, gr)
    cm = choice_maps(hh, g, l, zi, gr)
    fT = warm.(hh.hT[g, l, zi, :], p.ψ, p.href)
    return over_draws(cm, gr) do cmk, law, _
        fO = [Ewarm_nonteach(hh.capO[g, l, zi, k], law.XOgrid[k, g], g, p, law) for k in 1:law.nXO]
        sum(integrate_choice(fT, fO, cmk, law))
    end
end

function compute_lambda(hh, p::Params, gr)
    (; Nz, L, Πz) = gr
    iszero(p.λ) && return zeros(Nz, L)
    Ef = [(child_Ef(l, 1, zi, hh, p, gr) + child_Ef(l, 2, zi, hh, p, gr)) / 2
          for zi in 1:Nz, l in 1:L]
    return p.λ .* (Πz * Ef)
end

# -----------------------------------------------------------------------------
# 6. Household block: iterate V → policies → Λ → V to a fixed point.
# -----------------------------------------------------------------------------
"""
    solve_household(p, gr; tol, maxit, hh0, Λ0)

Solve every (gender, birth location, ability, shock) node given Λ, update Λ,
and repeat until the values converge. `hh0` warm-starts the policies (e, s) and,
through `Λ0`, the altruism term. Arrays are indexed [g, l, zi, k].

Returns the policies and values plus `infeasT`, `infeasO` (flags, true at the nodes with no
destination above the consumption floor at the policies, which `choice_maps` caps), and
`s_resid`/`s_passes`, the largest (e, s)
fixed-point residual and pass count over the nodes in the final sweep, and `s_relaxed`,
the largest number of nodes that needed the relaxed update (`solve_node`) in any one sweep:
a cold sweep may need it at many nodes, while the warm-started final sweep rarely does.
Under `p.college` the household block has two education branches and the college choice
(`solve_household_college`, §3c); its arrays are then `hh.br1[…]`, `hh.br0[…]` and `hh.P1`.
"""
function solve_household(p::Params, gr; tol = 1e-6, maxit = 500,
                         verbose = false, print_every = 25,
                         hh0 = nothing, Λ0 = hh0 === nothing ? nothing : hh0.Λ)
    p.college && return solve_household_college(p, gr; tol, maxit, verbose, print_every, hh0, Λ0)
    (; Nz, L, nϵT, nXO) = gr
    szT, szO = (2, L, Nz, nϵT), (2, L, Nz, nXO)
    eT = zeros(szT); sT = fill(s_of(p.T, p), szT);           WT = zeros(szT); hT   = zeros(szT)
    eO = zeros(szO); sO = fill(s_of(gr.nonteach[1], p), szO); WO = zeros(szO); capO = zeros(szO)
    infeasT, infeasO = falses(szT), falses(szO)       # no destination above the consumption floor
    if hh0 !== nothing && size(hh0.eT) == szT && size(hh0.eO) == szO
        eT .= hh0.eT; sT .= hh0.sT; eO .= hh0.eO; sO .= hh0.sO
    end
    Λ  = (Λ0 !== nothing && size(Λ0) == (Nz, L)) ? copy(Λ0) : zeros(Nz, L)
    hh = (; eT, sT, WT, hT, eO, sO, WO, capO, Λ, infeasT, infeasO)

    s_resid, s_passes, s_relaxed = 0.0, 0, 0
    t0 = time()
    for it in 1:maxit
        WT_old, WO_old = copy(WT), copy(WO)
        s_resid, s_passes, n_relaxed = 0.0, 0, 0
        # A zero e (unset) gives u0 = −Inf, which the e-solver treats as a cold start.
        for g in 1:2, l in 1:L, zi in 1:Nz
            for k in 1:nϵT
                ϵ  = gr.ϵTgrid[k]
                nd = teach_node(ϵ, zi, l, g, Λ, p, gr)
                e, V̄, s, r, n, rx = solve_node(nd, p, log(eT[g, l, zi, k]), sT[g, l, zi, k])
                eT[g, l, zi, k] = e
                sT[g, l, zi, k] = s
                WT[g, l, zi, k] = log(1 - s) + V̄ - p.cT - p.cTz * log(gr.z[zi])
                hT[g, l, zi, k] = teacher_h(ϵ, zi, l, s, e, p, gr)
                infeasT[g, l, zi, k] = all(≤(C_FLOOR), consumption(nd, s, e, p))
                s_resid, s_passes, n_relaxed = max(s_resid, r), max(s_passes, n), n_relaxed + rx
            end
            for k in 1:nXO
                x  = gr.XOgrid[k, g]
                nd = nonteach_node(x, zi, l, g, Λ, p, gr)
                e, V̄, s, r, n, rx = solve_node(nd, p, log(eO[g, l, zi, k]), sO[g, l, zi, k])
                eO[g, l, zi, k]   = e
                sO[g, l, zi, k]   = s
                WO[g, l, zi, k]   = log(1 - s) + V̄
                capO[g, l, zi, k] = nonteach_cap(x, zi, l, s, e, p, gr)
                infeasO[g, l, zi, k] = all(≤(C_FLOOR), consumption(nd, s, e, p))
                s_resid, s_passes, n_relaxed = max(s_resid, r), max(s_passes, n), n_relaxed + rx
            end
        end
        s_relaxed = max(s_relaxed, n_relaxed)
        Λ .= compute_lambda(hh, p, gr)
        err = max(maximum(abs, WT .- WT_old), maximum(abs, WO .- WO_old))
        if verbose && (it == 1 || it % print_every == 0 || err < tol)
            @printf("    HH %4d/%d  err=%.3e  elapsed=%.1fs\n", it, maxit, err, time() - t0)
            flush(stdout)
        end
        err < tol && break
    end
    return (; hh..., s_resid, s_passes, s_relaxed)
end

# -----------------------------------------------------------------------------
# 7. Location choice, stationary distribution and aggregates
# -----------------------------------------------------------------------------
"""
Logit work-location probabilities πT[g, l, zi, k, l'] and πO[g, l, zi, k, l'] at the policies.
Under `Params.split` the joint residence × workplace probabilities πT[g, l, zi, k, r, w] and
πO[g, l, zi, k, r, w] (`nested_location_probs`). Under `p.college` a pair ((πT₁, πT₀), (πO₁, πO₀)) for
the two education branches (`college_location_probs`), each of them the arrays of the baseline or,
under the split, the joint ones; `d` selects a branch's arrays `hh` (§3c). `who` (§3d, `node_probs`) picks
the adults' probabilities (`:all`, the default: with the fertility shock, the mixture of parents and the
childless), the parents' (`:parents`) or the childless's (`:childless`); without the shock all three are the same.
"""
function location_probs(hh, p::Params, gr; d = nothing, who::Symbol = :all)
    p.college && d === nothing && return college_location_probs(hh, p, gr; who)
    p.split && return nested_location_probs(hh, p, gr; d, who)
    (; Nz, L, nϵT, nXO) = gr
    πT = zeros(2, L, Nz, nϵT, L)
    πO = zeros(2, L, Nz, nXO, L)
    for g in 1:2, l in 1:L, zi in 1:Nz
        for k in 1:nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hh.Λ, p, gr; d)
            C  = consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p)
            πT[g, l, zi, k, :] .= node_probs(nd, C, p; who)
        end
        for k in 1:nXO
            nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hh.Λ, p, gr; d)
            C  = consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p)
            πO[g, l, zi, k, :] .= node_probs(nd, C, p; who)
        end
    end
    return πT, πO
end

"""
    parent_probs(hh, p, gr, πT, πO) -> (πTk, πOk)

The parents' location probabilities (§3d): where the next generation is born and goes to school, in the layout of
`location_probs`. Without the fertility shock every adult is a parent, and they are the adults' `πT`, `πO` as given.
With the ability gradient (`p.pkid_z` ≠ 0) each node's π₁ carries the weight P(child | z)/pkid (`gr.pk`), so that an
entry is P(parent in l′ | node)/pkid, the node's children in l′ per adult: everything that reads the parents' arrays
(Φ, the parents' and children's samples, the parents' WTP) then counts children and parents. The weight is 1 without
the gradient, which is not applied, so the arrays are bit for bit π₁.
"""
function parent_probs(hh, p::Params, gr, πT, πO)
    fertile(p) || return (πT, πO)
    πTk, πOk = location_probs(hh, p, gr; who = :parents)
    iszero(p.pkid_z) && return (πTk, πOk)
    c = gr.pk ./ p.pkid
    return fertility_weighted(πTk, c), fertility_weighted(πOk, c)
end

"Location probabilities [g, l, zi, …] (or the pair of the education branches) times `c[zi]` (§3d, `parent_probs`)."
function fertility_weighted(π::AbstractArray, c)
    out = copy(π)
    for zi in axes(out, 3)
        selectdim(out, 3, zi) .*= c[zi]
    end
    return out
end
fertility_weighted(π::Tuple, c) = map(x -> fertility_weighted(x, c), π)

"""
    migration_probs(πT, πO, hh, p, gr) -> πbar

`πbar[l, zi, l']`, the probability born-in-l → work-in-l' (under the split: lives-in-l′) at the location
probabilities `πT`, `πO`, integrated over the occupation choice and averaged over gender; under `p.college` it
mixes the education branches with the weights P(d | g, l, z) (§3c). With the adults' probabilities it is the
migration of the young; with the parents' (§3d, `parent_probs`) the kernel of the next generation's birthplace.
"""
function migration_probs(πT, πO, hh, p::Params, gr)
    (; Nz, L) = gr
    πbar = zeros(L, Nz, L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
        cm = choice_maps(br.h, g, l, zi, gr)
        for lp in 1:L
            IT, IO = integrate_choice(br.πT[g, l, zi, :, lp], br.πO[g, l, zi, :, lp], cm, gr)
            πbar[l, zi, lp] += br.w * (IT + IO) / 2
        end
    end
    return πbar
end

"""
    stationary_phi(πT, πO, hh, p, gr) -> (Φ, πbar)

`πbar[l, zi, l']` is the migration probability born-in-l → work-in-l' at `πT`, `πO` (`migration_probs`).
`Φ[zi, l]` is the per-gender stationary mass of the kernel P((l,z) → (l',z')) = πbar[l, z, l'] Πz[z, z'],
the young by birthplace. With the fertility shock (§3d) pass the parents' probabilities (`parent_probs`): children are
born where their parent lives, so the kernel is the parents' and the returned πbar is the parents' migration. With
its ability gradient (`p.pkid_z` ≠ 0) those carry the weight P(child | z)/pkid, so the rows of the kernel sum to it
(not to 1) and πbar[l, zi, l′] is the children in l′ per young of (l, zi); Φ is the kernel's Perron vector, by a power
iteration normalised at every step (the root is 1 by `fertility_probs`, up to the quadrature error of the rows).
"""
function stationary_phi(πT, πO, hh, p::Params, gr)
    (; Nz, L, Πz) = gr
    πbar = migration_probs(πT, πO, hh, p, gr)
    idx(l, zi) = (l - 1) * Nz + zi
    P = zeros(L * Nz, L * Nz)
    for l in 1:L, zi in 1:Nz, lp in 1:L, zpi in 1:Nz
        P[idx(l, zi), idx(lp, zpi)] = πbar[l, zi, lp] * Πz[zi, zpi]
    end
    Φ = (fertile(p) && !iszero(p.pkid_z) ? perron_left(P) : stationary(P)) .* (p.Mtot / 2)
    return reshape(Φ, Nz, L), πbar
end

"Left Perron vector of a nonnegative matrix K, summing to 1: power iteration v ← v′K/Σ(v′K) (§3d's tilted kernel)."
function perron_left(K; tol = 1e-13, max_iter = 100_000)
    v = fill(1 / size(K, 1), size(K, 1))
    for _ in 1:max_iter
        v_new = vec(v' * K)
        v_new ./= sum(v_new)
        converged = maximum(abs, v_new .- v) < tol
        v = v_new
        converged && break
    end
    return v
end

"""
    peer_log_ability(Φ, gr) -> Vector

Local peer ability by location, log ā_l = E[log z | young in l] − E*[log z]: the mean log persistent
ability of the students of l (the young schooled there, Φ[:, l]; under the split, where the household
lives) relative to its ergodic mean under Πz. ā_l is the geometric mean of the students' z relative to
the economy's. Since z′ ~ Πz(z, ·) whatever the location, the young as a whole have the ergodic
distribution of z in a stationary equilibrium (without the fertility shock's ability gradient), so the log ā_l average to zero over the young (up to
the quadrature error of Φ, of the order of its stationarity deviation) and the peer factor ā_l^{ψa}
of `build_grids` shifts school quality between locations without changing its overall level. An
empty location gets 0. With the fertility shock's ability gradient (§3d, `pkid_z` ≠ 0) the young are tilted away
from the ergodic distribution, which is why that switch is not allowed with peer effects (`check_params`).
It is the peer input of the default specification A (`peer_log_input`), and the children's ability by location in
every specification (`log_abar` of the moments).
"""
function peer_log_ability(Φ, gr)
    lz  = log.(gr.z)
    ref = dot(stationary(gr.Πz), lz)
    return [(m = sum(view(Φ, :, l)); m > 0 ? dot(view(Φ, :, l), lz) / m - ref : 0.0) for l in 1:gr.L]
end

"""
    adult_log_h(Φ, πT, πO, hh, p, gr) -> (; m, H, Hbar)

The adults by residence and their own human capital (2026-10-25): `m[l]` the mass of the adults who live in l,
`H[l]` = E[log h | adult resident of l] (NaN for an empty location) and `Hbar` = E[log h | all adults], the
m-weighted mean of `H`, so that Σ_l m_l (H_l − Hbar) = 0 up to rounding. The adults are the young of Φ at the
location choices `πT`, `πO`, which must be the adults' (`location_probs` with `who = :all`: under the fertility shock,
§3d, the childless are adults too), over every birth location, ability node, gender, occupation (teachers, market
non-teachers and home producers) and, under `college`, education branch with the weights P(d | g, l, z)
(`cell_branches`); under the split (joint probabilities) by residence r, from the residence marginals (method below).
h is the adult's own, acquired when young in her birth location: h_T for teachers (`hh.hT`); for non-teachers
h_O = Q_l (z ϵ_{i*})^α s^φ e^η = capO/(e^{ρc d} Θ_{i*}) (`capO` carries the credential premium of a BA, `credential`),
and given X_O* = x every choice depends on x alone, so the occupation that attains the max has P(i* = i | x) ∝
f_i(x)/F_i(x) among the choosers too and E[log h_O | x, choices] = log(capO/e^{ρc d}) − E[log Θ_{i*} | x]
(`argmax_mean`, exact for any number of occupations); this is `mean_child_h`'s h_O in logs. The masses and integrals
are those of the solver's choice maps (`integrate_choice`, cubic splines in the log shock, as `migration_probs` and
`location_stats` integrate), so `m` is the location (residence) mass of the adults that `migration_probs` implies.
Under constant class size (§3e) h is at the effective node, the true h for every draw, and E[log Θ_{i*} | X̃_O*] is taken
under each draw's law (`over_draws`).
"""
function adult_log_h(Φ, πT, πO, hh, p::Params, gr)
    (; Nz, L, nXO) = gr
    ccs = ccs_on(gr)
    EΘ = ccs ? nothing :
         [[argmax_mean(identity, gr, g, gr.XOgrid[k, g]) for k in 1:nXO] for g in 1:2]   # E[log Θ_{i*} | X_O*], by gender
    m = zeros(L); S = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm  = choice_maps(br.h, g, l, zi, gr)
        lhT = log.(br.h.hT[g, l, zi, :])
        if ccs                                   # §3e: log(capO/e^{ρc d}); E[log Θ_{i*}] by draw below
            lcO = log.(br.h.capO[g, l, zi, :] ./ credential(p, br.d))
            for lp in 1:L
                πTv, πOv = br.πT[g, l, zi, :, lp], br.πO[g, l, zi, :, lp]
                mT, mO = integrate_choice(πTv, πOv, cm, gr)
                sT = integrate_teach(πTv .* lhT, cm, gr)
                sO = over_draws(cm, gr) do cmk, law, _
                    EΘk = [argmax_mean(identity, law, g, law.XOgrid[k, g]) for k in 1:nXO]
                    integrate_nonteach(πOv .* (lcO .- EΘk), cmk, law)
                end
                m[lp] += mass * (mT + mO)
                S[lp] += mass * (sT + sO)
            end
            continue
        end
        lhO = log.(br.h.capO[g, l, zi, :] ./ credential(p, br.d)) .- EΘ[g]
        for lp in 1:L
            πTv, πOv = br.πT[g, l, zi, :, lp], br.πO[g, l, zi, :, lp]
            mT, mO = integrate_choice(πTv, πOv, cm, gr)
            sT, sO = integrate_choice(πTv .* lhT, πOv .* lhO, cm, gr)
            m[lp] += mass * (mT + mO)
            S[lp] += mass * (sT + sO)
        end
    end
    return (; m, H = [m[l] > 0 ? S[l] / m[l] : NaN for l in 1:L], Hbar = sum(S) / sum(m))
end

"""
    peer_spec(p) -> String

The peer specification of `p` (2026-10-25), "A" … "D" (`peer_log_input`): A the students' mean ability (the default),
B their total, C the adults' mean human capital, D their total.
"""
peer_spec(p::Params) = p.peer_input === :ability ? (p.peer_total ? "B" : "A") : (p.peer_total ? "D" : "C")

"Description of the peer specification of `p` (`peer_spec`)."
peer_spec_name(p::Params) = string(p.peer_total ? "total " : "mean ",
                                   p.peer_input === :ability ? "ability of the students" : "human capital of the adult residents")

"""
    peer_terms(Φ, πT, πO, hh, p, gr; adults) -> (; mean, scale, mass)

The two terms of the peer input log P_l (`peer_log_input`) by `p.peer_input`, whatever `p.peer_total`: `mass` the
peer group's size by location (the students n_l = Σ_z Φ[z, l] for :ability, the adult residents m_l for :hc),
`scale` = log(mass_l/(Σ_l mass_l/L)) its log relative to an average location's, and `mean` its mean term, log ā_l
(`peer_log_ability`) for :ability and H_l − Hbar (`adult_log_h`) for :hc. An empty location gets 0 in both terms.
`adults` passes `adult_log_h(Φ, πT, πO, hh, p, gr)` if it is already at hand (read only for :hc).
"""
function peer_terms(Φ, πT, πO, hh, p::Params, gr;
                    adults = p.peer_input === :hc ? adult_log_h(Φ, πT, πO, hh, p, gr) : nothing)
    if p.peer_input === :ability
        mterm = peer_log_ability(Φ, gr)
        mass  = vec(sum(Φ; dims = 1))
    else
        mass  = adults.m
        mterm = [mass[l] > 0 ? adults.H[l] - adults.Hbar : 0.0 for l in 1:gr.L]
    end
    avg = sum(mass) / gr.L
    return (; mean = mterm, scale = [mass[l] > 0 ? log(mass[l] / avg) : 0.0 for l in 1:gr.L], mass)
end

"""
    peer_log_input(Φ, πT, πO, hh, p, gr; adults) -> Vector

The peer input log P_l of school quality Q_l = Q_T,l P_l^{ψa} (`build_grids`; the solver's aggregate `la`) at the
stationary distribution Φ and the adults' location probabilities `πT`, `πO` (`location_probs`; the joint ones under
the split, pairs under `college`) of the household block `hh`, by `p.peer_input` and `p.peer_total` (2026-10-25):
- A (`:ability`, mean; the default): log ā_l = E[log z | student of l] − E*[log z], `peer_log_ability(Φ, gr)` itself,
  so that the default is the model without the switches bit for bit;
- B (`:ability`, total): log(n_l/n̄) + log ā_l, n_l = Σ_z Φ[z, l] the students of l and n̄ = Σ_l n_l/L;
- C (`:hc`, mean): E[log h | adult resident of l] − E[log h | all adults] (`adult_log_h`), against the endogenous
  adult-mass-weighted mean, so that Σ_l m_l log P_l = 0 and the factor moves quality across locations, not its
  level (as ā does in A; a fixed reference would feed Q back into its own level through h, explosively for ψa ≥ 1);
- D (`:hc`, total): log(m_l/m̄) + C, m_l the adult residents of l and m̄ = Σ_l m_l/L.
"Total" is number × geometric mean relative to an average location's (n̄ times the reference mean), so that B − A and
D − C are exactly the scale term (`peer_terms`). `πT`, `πO` and `hh` are read by C and D only. An empty location gets
0, as in `peer_log_ability`. `adults` as in `peer_terms`.
"""
function peer_log_input(Φ, πT, πO, hh, p::Params, gr;
                        adults = p.peer_input === :hc ? adult_log_h(Φ, πT, πO, hh, p, gr) : nothing)
    p.peer_input === :ability && !p.peer_total && return peer_log_ability(Φ, gr)    # A: exactly the pre-switch input
    t = peer_terms(Φ, πT, πO, hh, p, gr; adults)
    return p.peer_total ? t.mean .+ t.scale : t.mean
end

"""
    aggregates(Φ, πT, πO, hh, p, gr) -> (H̃_T, M, t)

By work location: the teacher aggregator H̃_T = ∫ h_T^{β/σ} over teachers,
population M, and the balanced-budget tax rate t = teacher wage bill / taxable income.
Under `p.college` the cells mix the education branches with the weights P(d | g, l, z),
teachers' pay carries the pay factor of its branch (`pay_factor`), non-teachers' income is
their capacity `capO` (with the credential premium), and H̃_T aggregates the teachers'
effectiveness (`effective_h`, h_T unless ωz ≠ 1) (§3c).
"""
function aggregates(Φ, πT, πO, hh, p::Params, gr)
    (; Nz, L) = gr
    expo = p.β / p.σ
    HT = zeros(L); income = zeros(L); wagebill = zeros(L)
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2                     # per gender
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
            mass = mass0 * br.w
            cm   = choice_maps(br.h, g, l, zi, gr)
            hTv  = br.h.hT[g, l, zi, :]
            hEv  = effective_h(hTv, zi, p, gr)
            pf   = teacher_pay_factor(br.h.sT[g, l, zi, :], p, br.d)
            capO = br.h.capO[g, l, zi, :]
            for lp in 1:L
                πTv = br.πT[g, l, zi, :, lp]
                wT  = integrate_teach(πTv .* p.κ[lp] .* hTv .^ p.γ .* pf, cm, gr)
                HT[lp]       += mass * integrate_teach(πTv .* hEv .^ expo, cm, gr)
                wagebill[lp] += mass * wT
                income[lp]   += mass * ((1 - p.τω[p.T, g]) * wT +
                                        integrate_nonteach(br.πO[g, l, zi, :, lp] .* capO, cm, gr))
            end
        end
    end
    M = 2 .* vec(sum(Φ; dims = 1))
    t = [income[lp] > 0 ? clamp(wagebill[lp] / income[lp], 0.0, TAX_CAP) : 0.0 for lp in 1:L]
    return HT, M, t
end

# Residence and workplace apart (§3b). `location_probs` returns the joint probabilities
# π[g, l, zi, k, r, w] (under `college`: a pair of such arrays, one per education branch), and
# `stationary_phi` and `aggregates` have methods for such arrays.
"The joint probabilities of the split: an array [g, l, zi, k, r, w], or the pair of them of the two education branches."
const JointProbs = Union{AbstractArray{<:Real,6}, Tuple{AbstractArray{<:Real,6}, AbstractArray{<:Real,6}}}

"Residence marginal Σ_w π[g, l, zi, k, r, w] of joint probabilities (of each branch under `college`)."
residence_probs(π::AbstractArray{<:Real,6}) = dropdims(sum(π; dims = 6); dims = 6)
residence_probs(π::Tuple) = map(residence_probs, π)

"Workplace marginal Σ_r π[g, l, zi, k, r, w] of joint probabilities (of each branch under `college`)."
workplace_probs(π::AbstractArray{<:Real,6}) = dropdims(sum(π; dims = 5); dims = 5)
workplace_probs(π::Tuple) = map(workplace_probs, π)

"The joint probabilities of node k in cell (g, l, zi) as a vector over the pairs d = (r − 1)L + w."
pair_probs(π, g, l, zi, k) = vec(permutedims(π[g, l, zi, k, :, :]))

"""
    nested_location_probs(hh, p, gr; d = nothing) -> (πTW, πOW)

Joint residence × workplace probabilities πTW[g, l, zi, k, r, w] and πOW[g, l, zi, k, r, w] at
the policies: `location_probs` under the split (`d`: the arrays of education branch `d` under `college`; `who`: the
adults', the parents' or the childless's, §3d).
"""
function nested_location_probs(hh, p::Params, gr; d = nothing, who::Symbol = :all)
    (; Nz, L, nϵT, nXO) = gr
    πT = zeros(2, L, Nz, nϵT, L, L)
    πO = zeros(2, L, Nz, nXO, L, L)
    for g in 1:2, l in 1:L, zi in 1:Nz
        for k in 1:nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hh.Λ, p, gr; d)
            C  = consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p)
            πT[g, l, zi, k, :, :] .= reshape(node_probs(nd, C, p; who), L, L)'
        end
        for k in 1:nXO
            nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hh.Λ, p, gr; d)
            C  = consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p)
            πO[g, l, zi, k, :, :] .= reshape(node_probs(nd, C, p; who), L, L)'
        end
    end
    return πT, πO
end

"""
`stationary_phi` for joint probabilities: children go to school where the household lives, so
only the residence marginals matter (`πbar[l, zi, r]` is the probability born-in-l → lives-in-r).
Under `college` the marginals are taken branch by branch and the baseline method mixes them.
"""
stationary_phi(πT::JointProbs, πO::JointProbs, hh, p::Params, gr) =
    stationary_phi(residence_probs(πT), residence_probs(πO), hh, p, gr)

"`migration_probs` for joint probabilities: born-in-l → lives-in-r, from the residence marginals."
migration_probs(πT::JointProbs, πO::JointProbs, hh, p::Params, gr) =
    migration_probs(residence_probs(πT), residence_probs(πO), hh, p, gr)

"`adult_log_h` for joint probabilities: the adults by residence r, from the residence marginals (wherever they work)."
adult_log_h(Φ, πT::JointProbs, πO::JointProbs, hh, p::Params, gr) =
    adult_log_h(Φ, residence_probs(πT), residence_probs(πO), hh, p, gr)

"""
`aggregates` for joint probabilities. The tax base is by residence, the teachers' side by workplace:
H̃_T[w] and the wage bill[w] integrate the teachers who work at w (wherever they live), the taxable
income[r] is the pre-commuting income of the residents of r, and t_r = wage bill[r]/income[r] is the
tax rate that funds the teachers who work in r with the income of those who live there. Under
`p.college` the cells mix the education branches with the weights P(d | g, l, z), with the pay
factor of the branch on the teachers' pay, the credential premium in the non-teachers' capacity
`capO` and the teachers' effectiveness in H̃_T, as in the baseline method (§3c).
"""
function aggregates(Φ, πT::JointProbs, πO::JointProbs, hh, p::Params, gr)
    (; Nz, L) = gr
    expo = p.β / p.σ
    HT = zeros(L); income = zeros(L); wagebill = zeros(L)
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2                     # per gender
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
            mass = mass0 * br.w
            cm   = choice_maps(br.h, g, l, zi, gr)
            hTv  = br.h.hT[g, l, zi, :]
            hEv  = effective_h(hTv, zi, p, gr)
            pf   = teacher_pay_factor(br.h.sT[g, l, zi, :], p, br.d)
            capO = br.h.capO[g, l, zi, :]
            for r in 1:L
                πOr = vec(sum(br.πO[g, l, zi, :, r, :]; dims = 2))        # non-teachers: residence only
                income[r] += mass * integrate_nonteach(πOr .* capO, cm, gr)
                for w in 1:L
                    πTv = br.πT[g, l, zi, :, r, w]
                    wT  = integrate_teach(πTv .* p.κ[w] .* hTv .^ p.γ .* pf, cm, gr)
                    HT[w]       += mass * integrate_teach(πTv .* hEv .^ expo, cm, gr)
                    wagebill[w] += mass * wT
                    income[r]   += mass * (1 - p.τω[p.T, g]) * wT
                end
            end
        end
    end
    M = 2 .* vec(sum(Φ; dims = 1))
    t = [income[r] > 0 ? clamp(wagebill[r] / income[r], 0.0, TAX_CAP) : 0.0 for r in 1:L]
    return HT, M, t
end

"""
    split_probs(hh, p, gr) -> (; Pi_T, Pi_O, Pi_Tr, Pi_Or, Pi_TW, Pi_OW, Φ, πbar)

The location-choice objects of a solution under the split: the joint probabilities `Pi_TW`,
`Pi_OW`, their residence marginals `Pi_Tr`, `Pi_Or` (where the children go to school) and the
stationary distribution. `Pi_T` and `Pi_O` are `nothing`, so that a consumer that has not been
migrated to the split throws instead of reading the wrong object. Under `p.college` the four
probability objects are pairs, one entry per education branch (`cell_branches`).
With the fertility shock (§3d) the probabilities are the adults' and `πbar` is their migration born-in-l →
lives-in-r, while Φ is the stationary distribution of the parents' kernel; the parents' objects are `Pi_TWk`,
`Pi_OWk`, `Pi_Trk`, `Pi_Ork` and `πbar_kid` (without the shock, the adults' own; `Pi_Tk = Pi_Ok = nothing`). With the
shock's ability gradient (`pkid_z` ≠ 0) they carry the weight P(child | z)/pkid (`parent_probs`): children per adult,
not probabilities.
"""
function split_probs(hh, p::Params, gr)
    Pi_TW, Pi_OW = location_probs(hh, p, gr)
    Pi_TWk, Pi_OWk = parent_probs(hh, p, gr, Pi_TW, Pi_OW)
    Φ, πbar_kid = stationary_phi(Pi_TWk, Pi_OWk, hh, p, gr)
    πbar = fertile(p) ? migration_probs(Pi_TW, Pi_OW, hh, p, gr) : πbar_kid
    Pi_Tr, Pi_Or = residence_probs(Pi_TW), residence_probs(Pi_OW)
    Pi_Trk, Pi_Ork = fertile(p) ? (residence_probs(Pi_TWk), residence_probs(Pi_OWk)) : (Pi_Tr, Pi_Or)
    return (; Pi_T = nothing, Pi_O = nothing, Pi_Tr, Pi_Or, Pi_TW, Pi_OW, Φ, πbar,
              Pi_Tk = nothing, Pi_Ok = nothing, Pi_TWk, Pi_OWk, Pi_Trk, Pi_Ork, πbar_kid)
end

"""
    teachers_by_workplace(Φ, πT, πO, hh, p, gr) -> Vector

Teacher mass by workplace w under the split (joint probabilities; both genders, every birth location, ability
node and education branch), integrated exactly as the calibration's `location_stats` integrates it: the
numerator of the FTE-per-pupil gap that the teacher-demand closure (`Params.fte_rule`) holds.
"""
function teachers_by_workplace(Φ, πT::JointProbs, πO::JointProbs, hh, p::Params, gr)
    (; Nz, L) = gr
    n = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        for r in 1:L, w in 1:L
            n[w] += mass * integrate_teach(br.πT[g, l, zi, :, r, w], cm, gr)
        end
    end
    return n
end

"""
    fte_gap(Φ, πT, πO, hh, p, gr) -> Float64

Δlog teachers per pupil, log(n_2/S_2) − log(n_1/S_1), with teachers by workplace (`teachers_by_workplace`) and
pupils by residence, S_l = Σ_z Φ[z, l]: the calibration's `gap_teachers_pp` at these choices.
"""
function fte_gap(Φ, πT::JointProbs, πO::JointProbs, hh, p::Params, gr)
    n = teachers_by_workplace(Φ, πT, πO, hh, p, gr)
    S = vec(sum(Φ; dims = 1))
    return log(n[2] / S[2]) - log(n[1] / S[1])
end

"""
    ccs_teacher_moments(Φ, πT, πO, hh, p, gr; c0 = nothing, order = 2 p.ccs_nq) -> (; n, μ, c)

The teachers by workplace w and their effectiveness (§3e, 2026-10-26): `n[w]` the teacher mass, `μ[w]` = E[ln h_eff |
teacher works in w] (`effective_h`) and `c[w]` the central moments of ln h_eff about μ[w], orders 0..`order` in
`c[w][1:order+1]`, normalised to mass 1 (`c[w][1]` = 1, `c[w][2]` = 0); NaN for a workplace without teachers. These are the
state of constant class size: ln Q̄_w = β μ_w − σ ln(S_w/n_w) and the rule R_w (`ccs_rule`). Cells are mass-weighted by
Φ[zi, l]/2 per gender (both genders, every birth location and, with `cell_branches`, branch), as `aggregates` does, and
`πT` are the work-location probabilities (the adults'; the joint ones under the split, method below, by workplace).
For each (cell, draw, workplace) ONE vector-valued `quadgk` pass integrates, in the log effective shock u and split at
the grid nodes,
    pdf(dT, eᵘ) eᵘ · teach_wt_log(cm_k, u) · π_w(u) · (lh(u) − c0_w)^j,    j = 0..order,
with π_w the linear (k = 1, so nonnegative) spline of the location probabilities and lh the cubic spline of ln h_eff, plus
the end-point masses as `integrate_teach` adds them; the draws are mixed with `over_draws` (one pass at `gr` without CCS,
so a baseline solution works too: the moments of its teachers). `n` and `μ` themselves come from `integrate_teach` (cubic
splines of π and π ln h_eff, as `aggregates` integrates H̃_T), so at β = 0 the CCS n_T equals the baseline's H̃_T; the
pass's own mass and mean differ from them at the spline error (≈ 1e-6) and are used only to normalise and centre
`c`. Every power is integrated against the same nonnegative
discrete measure (the Gauss–Kronrod nodes and weights the pass chose, the end masses), so the raw moments about c0 are a
moment sequence and their Hankel matrix is positive semidefinite by construction; the shift to central moments about μ
(binomial) is exact algebra on that measure. (Separate `integrate_teach` calls per power would build a spline per power and
need not be.) `c0` is the centre, ideally the previous μ; `nothing` (or any non-finite entry) first runs an order-1 pass
about 0 for it.
"""
function ccs_teacher_moments(Φ, πT, πO, hh, p::Params, gr; c0 = nothing, order::Int = 2 * p.ccs_nq)
    (; Nz, L, lϵT) = gr
    if c0 === nothing || !all(isfinite, c0)
        c0 = order ≤ 1 ? zeros(L) : ccs_teacher_moments(Φ, πT, πO, hh, p, gr; c0 = zeros(L), order = 1).μ
        c0 = [isfinite(v) ? v : 0.0 for v in c0]
    end
    J = order
    raw = [zeros(J + 1) for _ in 1:L]
    n, s1 = zeros(L), zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, πT, πO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm  = choice_maps(br.h, g, l, zi, gr)
        lh  = log.(effective_h(br.h.hT[g, l, zi, :], zi, p, gr))
        slh = Spline1D(lϵT, lh)
        for w in 1:L
            πv = br.πT[g, l, zi, :, w]
            all(iszero, πv) && continue
            # n and μ as `aggregates` and `location_stats` integrate them (cubic): at β = 0 the CCS n_T is the baseline's
            # H̃_T, so the two models nest to rounding
            n[w]  += mass * integrate_teach(πv, cm, gr)
            s1[w] += mass * integrate_teach(πv .* lh, cm, gr)
            sπ = Spline1D(lϵT, πv; k = 1)
            v  = over_draws(cm, gr) do cmk, law, _
                teacher_moment_pass(cmk, law, slh, sπ, lh, πv, c0[w], J)
            end
            raw[w] .+= mass .* v
        end
    end
    μ = [n[w] > 0.0 ? s1[w] / n[w] : NaN for w in 1:L]
    c = map(1:L) do w
        raw[w][1] > 0.0 || return fill(NaN, J + 1)
        m = raw[w] ./ raw[w][1]                  # raw moments about c0, mass 1 (the pass's own measure)
        d = m[2]                                  # μ − c0
        [sum(binomial(j, i) * m[i + 1] * (-d)^(j - i) for i in 0:j) for j in 0:J]
    end
    return (; n, μ, c)
end

"""
    teacher_moment_pass(cm, law, slh, sπ, lh, πv, c0, J) -> Vector

One vector-valued quadrature pass of `ccs_teacher_moments` for a cell, a draw (`cm`, `law`: a `ChoiceMaps` and the
law it holds) and a workplace: ∫ pdf(dT, eᵘ) eᵘ P(teach | u) π_w(u) (lh(u) − c0)^j du, j = 0..J, over the grid (split at
its nodes) plus the end-point masses at the grid's end values, as `integrate_teach` does.
"""
function teacher_moment_pass(cm::ChoiceMaps, law, slh, sπ, lh, πv, c0, J)
    (; dT, ϵTgrid, lϵT) = law
    function f(u)
        wt = teach_wt_log(cm, u)
        wt <= 0.0 && return zeros(J + 1)
        ϵ = exp(u)
        a = pdf(dT, ϵ) * ϵ * wt * max(sπ(u), 0.0)
        d = slh(u) - c0
        return [a * d^j for j in 0:J]
    end
    I  = quadgk(f, lϵT...)[1]
    lo, hi = ϵTgrid[1], ϵTgrid[end]
    a1 = cdf(dT, lo) * πv[1] * teach_wt(cm, lo)
    a2 = ccdf(dT, hi) * πv[end] * teach_wt(cm, hi)
    return I .+ a1 .* (lh[1] - c0) .^ (0:J) .+ a2 .* (lh[end] - c0) .^ (0:J)
end

"`ccs_teacher_moments` for joint probabilities (the split): the teachers by workplace w, from the workplace marginals."
ccs_teacher_moments(Φ, πT::JointProbs, πO::JointProbs, hh, p::Params, gr; kw...) =
    ccs_teacher_moments(Φ, workplace_probs(πT), workplace_probs(πO), hh, p, gr; kw...)

"Initial slope ∂(FTE gap)/∂b of the teacher-demand closure, per util (refined by secant steps in `solve_ge`)."
const FTE_SLOPE0 = 3.0
"Bounds on the closure's slope. Too small a slope overshoots, so the lower bound is the binding safeguard."
const FTE_SLOPE_MIN, FTE_SLOPE_MAX = 1.5, 48.0
"Largest step of the closure's b in one outer iteration (utils): far from the fixed point the FTE gap is no guide."
const FTE_STEP_MAX = 0.25
"b moves only once the other aggregates change by less than this (far from the fixed point the FTE gap is no guide)."
const FTE_SETTLE = 1e-2
"A secant pair is used only if both iterates had the other aggregates within this and b moved by more than `FTE_SECANT_DB`."
const FTE_SECANT_ERR, FTE_SECANT_DB = 1e-4, 1e-4

# -----------------------------------------------------------------------------
# 8. General equilibrium: damped fixed point in (H̃_T, M, t), and ā with peer effects
# -----------------------------------------------------------------------------
fmtvec(v) = "(" * join((@sprintf("%.4f", x) for x in v), ", ") * ")"
reldiff(new, old; floor = 1e-3) = maximum(abs.(new .- old) ./ max.(abs.(old), floor))

# Constant class size in the outer loop (§3e, 2026-10-26): the state (ln Q̄, teacher moments c, centre μ) by location,
# its start, its implied update, the change of the rules, the grid pad and the warm start across grid sizes.
"The degenerate teacher moments of order 0..2n_q (mass 1, no dispersion): the cold start of `c`, a one-point rule."
ccs_degenerate(p::Params) = vcat(1.0, zeros(2 * p.ccs_nq))

"""
    ccs_start(p, HT, M, init) -> (; lQT, c, μ)

The CCS state at the start of `solve_ge`. With a warm start `init` of the same L (a CCS or a baseline solution): ln Q̄
from `init.ccs.lQT` if it carries it, otherwise the efficient σ ln(2H̃_T/M) of `init`; the moments `c` and centres `μ`
recomputed from `init` with `ccs_teacher_moments` at `p`'s order (at `init`'s own grid and laws, so a different n_q or
pad, or a baseline `init`, is never misread), degenerate where `init` has no teachers. Without it: the baseline's initial
quality σ ln(2H̃_T/M) at the default H̃_T, M (`HT`, `M` as passed), degenerate moments and no centre (NaN).
"""
function ccs_start(p::Params, HT, M, init)
    L = length(HT)
    if init !== nothing && length(init.HT) == L
        lQT = hasproperty(init, :ccs) && init.ccs !== nothing ? collect(Float64, init.ccs.lQT) :
              p.σ .* log.(2 .* max.(init.HT, M_FLOOR) ./ max.(init.M, M_FLOOR))
        tm = ccs_teacher_moments(init.Φ, sol_probs(init)..., init.hh, p, init.gr)
        ok = [tm.n[l] > M_FLOOR && isfinite(tm.μ[l]) && all(isfinite, tm.c[l]) for l in 1:L]
        return (; lQT, c = [ok[l] ? tm.c[l] : ccs_degenerate(p) for l in 1:L], μ = [ok[l] ? tm.μ[l] : NaN for l in 1:L])
    end
    return (; lQT = p.σ .* log.(2 .* max.(HT, M_FLOOR) ./ max.(M, M_FLOOR)),
              c = [ccs_degenerate(p) for _ in 1:L], μ = fill(NaN, L))
end

"""
    ccs_implied(tm, Φ, live, lQT, c, μ, p) -> (; lQT, c, μ, ok)

The undamped CCS state implied by the teacher moments `tm` (`ccs_teacher_moments`) at the stationary distribution Φ:
ln Q̄_l = β μ_l − σ ln(S_l/n_T,l) with S_l = Σ_z Φ[z, l] the pupils, the moments `tm.c` and the centres `tm.μ`. A location
that is not `live` or has no teachers (n_T ≤ M_FLOOR, or non-finite moments) holds its current state, as `solve_ge` holds
t and log P there (`ok` false).
"""
function ccs_implied(tm, Φ, live, lQT, c, μ, p::Params)
    S  = vec(sum(Φ; dims = 1))
    L  = length(S)
    ok = [live[l] && tm.n[l] > M_FLOOR && isfinite(tm.μ[l]) && all(isfinite, tm.c[l]) for l in 1:L]
    return (; lQT = [ok[l] ? p.β * tm.μ[l] - p.σ * log(S[l] / tm.n[l]) : lQT[l] for l in 1:L],
              c = [ok[l] ? copy(tm.c[l]) : copy(c[l]) for l in 1:L],
              μ = [ok[l] ? tm.μ[l] : μ[l] for l in 1:L], ok)
end

"""
    ccs_rule_change(rules, c, p) -> Float64

The largest change, over the locations, between the rules `rules` (a grid's `gr.ccs`, with fields `x`, `w` by location)
and those of the moments `c` (`ccs_rule`): max |Δx|, |Δw| when a location's rule keeps its order, 1.0 when it changes
order (a node dropped or restored), which holds the outer loop until the order settles.
"""
function ccs_rule_change(rules, c, p::Params)
    Δ = 0.0
    for l in eachindex(c)
        r = ccs_rule(c[l], p)
        if length(r.x) == length(rules.x[l])
            Δ = max(Δ, maximum(abs, r.x .- rules.x[l]), maximum(abs, r.w .- rules.w[l]))
        else
            Δ = max(Δ, 1.0)
        end
    end
    return Δ
end

"max_l,k |x_k| over the rules of the moments `c` by location (`ccs_rule`)."
ccs_max_x(c, p::Params) = maximum(maximum(abs, ccs_rule(v, p).x) for v in c)

"The automatic CCS grid pad of the moments `c`: 1.25 × max_l,k |x_k| of their rules, rounded up to a multiple of 0.05."
ccs_auto_pad(c, p::Params) = max(0.0, 0.05 * ceil(1.25 * ccs_max_x(c, p) / 0.05 - 1e-9))

"The damped moments (1 − d) c + d cn by location: the moments of the mixture of the two centred laws, so PSD stays PSD."
ccs_damp(c, cn, damping) = [(1 - damping) .* c[l] .+ damping .* cn[l] for l in eachindex(c)]

"SD(ln r) of the rules of a CCS grid, by location (`rule_sd`)."
ccs_sds(gr) = [rule_sd(gr.ccs.x[l], gr.ccs.w[l]) for l in 1:gr.L]

"Linear interpolation at `xn` of the values `v` on the increasing nodes `xo`, flat beyond both ends."
function interp_flat(xo, v, xn)
    n = length(xo)
    out = zeros(length(xn))
    for (i, x) in enumerate(xn)
        if x <= xo[1]
            out[i] = v[1]
        elseif x >= xo[n]
            out[i] = v[n]
        else
            j = searchsortedlast(xo, x)
            θ = (x - xo[j]) / (xo[j + 1] - xo[j])
            out[i] = (1 - θ) * v[j] + θ * v[j + 1]
        end
    end
    return out
end

"""
    ccs_regrid_warm(hh0, gr_old, gr_new) -> hh0′

The warm start of `solve_household` carried across a change of the CCS grid (a wider pad, or a warm start from a
solution on other grids): the policies (eT, sT) and (eO, sO) of `hh0`, solved on `gr_old`, interpolated onto `gr_new`
in the log shock, linearly (log e for e; flat beyond the ends), with Λ as it is. Returns `hh0` itself when the grid
sizes already match, when either is missing, or when the (g, l, z) layout differs (`solve_household` then cold-starts).
"""
function ccs_regrid_warm(hh0, gr_old, gr_new)
    (hh0 === nothing || gr_old === nothing) && return hh0
    (; Nz, L, nϵT, nXO) = gr_new
    (size(hh0.eT) == (2, L, Nz, nϵT) && size(hh0.eO) == (2, L, Nz, nXO)) && return hh0
    (size(hh0.eT, 4) == length(gr_old.lϵT) && size(hh0.eO, 4) == size(gr_old.lXO, 1)) || return hh0
    (size(hh0.eT)[1:3] == (2, L, Nz) && size(hh0.eO)[1:3] == (2, L, Nz)) || return hh0
    eT = zeros(2, L, Nz, nϵT); sT = zeros(2, L, Nz, nϵT)
    eO = zeros(2, L, Nz, nXO); sO = zeros(2, L, Nz, nXO)
    for g in 1:2, l in 1:L, zi in 1:Nz
        e = hh0.eT[g, l, zi, :]
        all(>(0.0), e) && (eT[g, l, zi, :] .= exp.(interp_flat(gr_old.lϵT, log.(e), gr_new.lϵT)))
        sT[g, l, zi, :] .= interp_flat(gr_old.lϵT, hh0.sT[g, l, zi, :], gr_new.lϵT)
        e = hh0.eO[g, l, zi, :]
        all(>(0.0), e) && (eO[g, l, zi, :] .= exp.(interp_flat(gr_old.lXO[:, g], log.(e), gr_new.lXO[:, g])))
        sO[g, l, zi, :] .= interp_flat(gr_old.lXO[:, g], hh0.sO[g, l, zi, :], gr_new.lXO[:, g])
    end
    return (; eT, sT, eO, sO, Λ = hh0.Λ)
end

"""
    ccs_solution(hh, p, gr, HT, M, t, la, converged, s_relaxed, lQT, c, μ, pad) -> sol

The solution of `solve_ge` under CCS: the baseline's (or, under the split, `split_probs`') objects at the final household
solve, plus `ccs` = (; lQT, c, nT, μ, N, x, w, pad, nq, draws, rho, sd, lQT_now, c_now): `lQT` and `c`, the state the
returned grid was built at (ln Q̄ by location and the teachers' central moments of ln h_eff, orders 0..2n_q); `nT`, `μ`
and `c_now` the teacher mass, mean and central moments by workplace measured at the returned solution
(`ccs_teacher_moments`), `N` = S/n_T the class size, `lQT_now` = β μ − σ ln N the ln Q̄ they imply (= `lQT` at a fixed
point); `x`, `w` the nodes and weights of the rules R_l (as `gr.ccs`), `sd` their SD(ln r) (`rule_sd`), and `pad`, `nq`,
`draws`, `rho` the grid pad and the switches; `lQeff` = σ ln(H̃_T/S) = ln (2H̃_T/M)^σ with H̃_T measured at the returned
solution (`aggregates`), the efficient-equivalent quality the same teachers would give under Proposition 1. The
solution's `HT` is that measured H̃_T, not the loop's damped `HT` (under CCS H̃_T is not a state and not in the
convergence error, so the damped value can trail it).
"""
function ccs_solution(hh, p::Params, gr, HT, M, t, la, converged, s_relaxed, lQT, c, μ, pad)
    if p.split
        sol = (; hh, split_probs(hh, p, gr)..., HT, M, t, la, gr, p, converged,
                 s_resid = hh.s_resid, s_passes = hh.s_passes, s_relaxed)
        πT, πO = sol.Pi_TW, sol.Pi_OW
    else
        Pi_T, Pi_O = location_probs(hh, p, gr)
        Pi_Tk, Pi_Ok = parent_probs(hh, p, gr, Pi_T, Pi_O)
        Φ, πbar_kid = stationary_phi(Pi_Tk, Pi_Ok, hh, p, gr)
        πbar = fertile(p) ? migration_probs(Pi_T, Pi_O, hh, p, gr) : πbar_kid
        sol = (; hh, Pi_T, Pi_O, Φ, πbar, HT, M, t, la, gr, p, converged,
                 s_resid = hh.s_resid, s_passes = hh.s_passes, s_relaxed, Pi_Tk, Pi_Ok, πbar_kid)
        πT, πO = Pi_T, Pi_O
    end
    tm = ccs_teacher_moments(sol.Φ, πT, πO, hh, p, gr; c0 = μ)
    S  = vec(sum(sol.Φ; dims = 1))
    N  = [tm.n[l] > M_FLOOR ? S[l] / tm.n[l] : NaN for l in 1:gr.L]
    HTnow = aggregates(sol.Φ, πT, πO, hh, p, gr)[1]
    return merge(sol, (; HT = HTnow, ccs = (; lQT = copy(lQT), c = [copy(v) for v in c], nT = tm.n, μ = tm.μ, N,
                                  x = gr.ccs.x, w = gr.ccs.w, pad, nq = p.ccs_nq, draws = p.ccs_draws, rho = p.ccs_rho,
                                  sd = ccs_sds(gr),
                                  lQT_now = p.β .* tm.μ .- p.σ .* log.(N), c_now = tm.c,
                                  lQeff = p.σ .* log.(max.(HTnow, M_FLOOR) ./ max.(S, M_FLOOR)))))
end

"""
    solve_ge(p = Params(); Nz, nϵT, nXO, damping, tol, maxit, hh_tol, hh_maxit,
             verbose, init, q_lo, q_hi, la_start, ccs_pad)

Each iteration rebuilds the grids at the current (H̃_T, M, t), solves the household
problem (warm-started from the previous iteration), computes location choices, Φ
and the implied aggregates, and takes a damped step. Convergence is on the max
relative change in H̃_T and M and absolute change in t. `init` is a previous
solution to continue from. A final solve at the converged aggregates makes the
returned objects mutually consistent. `s_resid` and `s_passes` are those of the final
solve, `s_relaxed` the largest number of nodes that needed the relaxed (e, s) update in
any household sweep of the whole solve (`solve_household`). Under `p.split` the solution
carries `Pi_T = Pi_O = nothing` and, in their place, the residence marginals `Pi_Tr`, `Pi_Or`
and the joint residence × workplace probabilities `Pi_TW`, `Pi_OW` (`split_probs`); the
aggregates H̃_T and the wage bill are by workplace, M and the tax rate t by residence.
With peer effects (`p.ψa ≠ 0`) the local peer input `la` = log P_l (`peer_log_input`: by default the students'
mean ability log ā_l, `peer_log_ability`; specifications B–D by `p.peer_input`, `p.peer_total`) is a
fourth aggregate: the grids are built at it, it is recomputed from Φ and the iteration's adults' location
probabilities and household block, damped like the others and converges in absolute change; it starts at
`init.la` or, without a warm start that carries it, at `la_start`·(−1, +1) (city, suburb; two locations): under A
and C 0, the default, is no sorting, and a positive value starts from a suburb-sorted state, which can select that
branch when a strong peer effect makes the stationary equilibrium multiple (2026-10-25). Under B and D the scale
term alone puts the gap in log P near the pupils' log gap (0.376 at the baseline), so a suburb-sorted start needs
la_start above half of it, and small values only change the path. With ψa = 0 it stays at 0 whatever
`la_start`, and the iterations are those of the model without peers. The solution carries `la`, the value its
grid (and so its Q) was built at.
With the teacher-demand closure (`p.fte_rule` finite, round 2026-10-15c) the teachers' value of the suburban
workplace b = BT[2] is a further aggregate: each iteration measures the FTE-per-pupil gap at the choices
(`fte_gap`) and, once the other aggregates move by less than `FTE_SETTLE`, moves b by the damped step
damping·(fte_rule − gap)/slope. The slope starts at `FTE_SLOPE0` and is refined by secant steps between iterates whose
other aggregates had settled to `FTE_SECANT_ERR` and whose b differ by more than `FTE_SECANT_DB`. It is held within
[`FTE_SLOPE_MIN`, `FTE_SLOPE_MAX`] and doubled when the residual changes sign and grows, an overshoot. (The first
version refined the slope from any consecutive iterates with a relative clamp: near the fixed point, aggregate drift
contaminated the tiny b steps and the slope fell until b cycled at ±`FTE_STEP_MAX`, 2026-10-15.) Convergence also
requires |fte_rule − gap| < tol. b starts at `p.BT[2]` (or at `init.p.BT[2]`), and the returned `p` carries the
equilibrium b in `BT`, so that everything downstream (moments, partial-equilibrium re-solves) holds it fixed.
Without the closure the iterations are those of the model without it.
With constant class size (`p.ccs`, §3e, 2026-10-26) H̃_T is no longer a state (it is still computed and damped, the
efficient-equivalent index (2H̃_T/M)^σ of the same teachers): the outer state adds ln Q̄_l and the teachers' central
moments c_l of ln h_eff (orders 0..2n_q) by location. Each iteration builds the grids at them (`build_grids` with
`ccs = (; lQT, c)` and the pad), measures the teachers by workplace (`ccs_teacher_moments`, centred at the previous μ) and
takes the damped step to ln Q̄ = β μ − σ ln(S/n_T) and to their moments (`ccs_implied`; a convex combination of centred
moment vectors is the mixture's, so the Hankel matrices stay PSD), holding an empty location's state as t and log P are.
The error drops the H̃_T term and adds max |Δ ln Q̄| and the change of the rules (`ccs_rule_change`), BEFORE the closure's
`FTE_SETTLE` gate and its secant record (so b does not move while ln Q̄ and the moments do). Start: a warm start (CCS or
baseline) gives ln Q̄ (`init.ccs.lQT`, or its efficient σ ln(2H̃_T/M)) and moments recomputed from it (`ccs_start`); a cold
start has the baseline's initial quality and degenerate moments. `ccs_pad` (log units of r) is the grid pad: `:auto` =
1.25 max_l,k |x| rounded up to 0.05 and only ever widened, or a fixed number; when the grid grows the household warm
start is interpolated onto it (`ccs_regrid_warm`). The solution carries `ccs` (`ccs_solution`); `verify_solution` checks
max |x| ≤ pad. Without `ccs` none of this is reached and the iterations are those of the model without it.
"""
function solve_ge(p::Params = Params();
                  Nz = 5, nϵT = 64, nXO = 64, damping = 0.3, tol = 1e-5,
                  maxit = 1000, hh_tol = 1e-6, hh_maxit = 1000, verbose = true,
                  init = nothing, q_lo = 1e-5, q_hi = 1 - 1e-6, la_start = 0.0, ccs_pad = :auto)
    check_params(p)
    I, L = length(p.A), length(p.B)
    (ccs_pad === :auto || (ccs_pad isa Real && ccs_pad ≥ 0)) ||
        throw(ArgumentError("ccs_pad must be :auto or a nonnegative number (log units of r), got $ccs_pad"))
    (iszero(p.ψa) || iszero(la_start)) || @assert L == 2 "la_start starts the peer input at la_start·(−1, +1): it needs two locations (city, suburb)"
    grid_kw = (; Nz, nϵT, nXO, q_lo, q_hi)

    if verbose
        @printf("  [ε=0, continuous]  I=%d, L=%d, Nz=%d  shock grids: nϵT=%d, nXO=%d\n",
                I, L, Nz, nϵT, nXO)
        print("  Ξ=0 reference s: ")
        for i in 1:I
            @printf("%s s=%.4f   ", i == p.T ? "T" : "occ$i", s_of(i, p))
        end
        @printf("| ψ=%.3f href=%.4f  max mcost=%.4f\n", p.ψ, p.href, maximum(p.mcost))
        p.split && @printf("  split: σw=%.4f (σw/σν=%.3f), τk=%.4f, jobs=%s, BT=%s\n", p.σw, p.σw / p.σν, p.τk,
                           fmtvec(p.jobs), fmtvec(p.BT))
        p.college && @printf("  college: σd=%.3f sBA=%.3f ρc=%.3f ρT0=%.3f ϑs=%.3f Fc=%.4f cC=%.3f cCz=%.3f ωz=%.2f sT0=%.3f cT0=%.3f\n",
                             p.σd, p.sBA, p.ρc, p.ρT0, p.ϑs, p.Fc, p.cC, p.cCz, p.ωz, p.sT0, p.cT0)
        isnan(p.schedT) || @printf("  salary schedule: schedT=%.4f (pay κ h^γ (s/sBA)^%.4f, sBA=%.4f; teachers' s at Ξ = 0: %.4f)\n",
                                   p.schedT, teach_schedule(p)[2], p.sBA, time_share(p.γ, p, 0.0, teach_schedule(p)[2]))
        iszero(p.ψa) || @printf("  peers: ψa=%.3f, specification %s, the %s (Q_l = (2H̃_T/M_l)^σ P_l^ψa)%s\n", p.ψa,
                                peer_spec(p), peer_spec_name(p), iszero(la_start) ? "" : @sprintf("; cold start log P = %.4f·(−1, +1)", la_start))
        iszero(p.κc) || @printf("  commuting friction: κc=%.4f utils per commuting pair (r ≠ w)\n", p.κc)
        all(iszero, p.Bd) || @printf("  BA amenity: Bd=%s utils (residence r, BA holders)\n", fmtvec(p.Bd))
        all(iszero, p.By) || @printf("  income-graded taste: By=%s utils per log point of predicted income (lyref=%.4f)\n",
                                     fmtvec(p.By), p.lyref)
        isnan(p.fte_rule) || @printf("  teacher-demand closure: positions per pupil gap %.4f; BT[2] solved in GE\n", p.fte_rule)
        isnan(p.lyT_fixed) || @printf("  diagnostic: teachers' income taste graded at ln ŷ = lyT_fixed = %.4f (non-teachers: their own)%s\n",
                                      p.lyT_fixed, all(iszero, p.By) ? "  [By = 0: inert]" : "")
        fertile(p) && @printf("  fertility shock: P(child) = pkid = %.4f; only parents value Λ, children follow the parents (§3d)\n", p.pkid)
        fertile(p) && !iszero(p.pkid_z) && @printf("    ability-graded: logit slope pkid_z = %.4f per ergodic SD of log z; 1/pkid children per parent\n", p.pkid_z)
        p.ccs && @printf("  constant class size (§3e): N_l = S_l/n_T,l, random assignment; n_q = %d Gauss nodes, K_d = %d teacher draws, ρ = %.3f; pad %s\n",
                         p.ccs_nq, p.ccs_draws, p.ccs_rho, ccs_pad === :auto ? ":auto" : @sprintf("%.3f", ccs_pad))
        flush(stdout)
    end

    HT = fill(1.0, L); M = fill(p.Mtot / L, L); t = fill(0.10, L)
    la = zeros(L)                                 # log P_l, the peer input: an aggregate only with peer effects
    peer = !iszero(p.ψa)
    la_warm = false                               # whether a warm start supplies la
    closure = p.split && !isnan(p.fte_rule)       # teacher-demand closure: b = BT[2] is an aggregate
    b = closure ? p.BT[2] : NaN                  # (BT has one entry per location: read it only under the closure)
    hh0 = nothing
    if init !== nothing
        if length(init.HT) == L
            HT .= init.HT; M .= init.M; t .= init.t
            peer && hasproperty(init, :la) && (la .= init.la; la_warm = true)
            closure && hasproperty(init, :p) && (b = init.p.BT[2])
            hh0 = init.hh
        elseif verbose
            @printf("  ⚠ warm start ignored: init has L=%d ≠ %d\n", length(init.HT), L)
        end
    end
    peer && !la_warm && !iszero(la_start) && (la .= la_start .* [-1.0, 1.0])   # a sorted cold start (2026-10-25)
    closure && (p = with_params(p; BT = [p.BT[1], b]))
    slope, b_prev, f_prev, e_prev, r_prev = FTE_SLOPE0, NaN, NaN, NaN, NaN
    ccs = p.ccs                                   # §3e: constant class size, state (ln Q̄, moments, centres) by location
    if ccs
        cs = ccs_start(p, HT, M, init)            # HT, M are init's when it is a warm start of the same L
        lQT, cmom, cμ = cs.lQT, cs.c, cs.μ
        pad = ccs_pad === :auto ? 0.0 : Float64(ccs_pad)
        gr_prev = init !== nothing && length(init.HT) == L ? init.gr : nothing   # the grid hh0 was solved on
        pad_warned = false
    end

    converged = false
    s_relaxed = 0
    for it in 1:maxit
        if ccs                                    # §3e: the pad covers the rules of the current moments (:auto only widens it)
            ccs_pad === :auto && (pad = max(pad, ccs_auto_pad(cmom, p)))
            if verbose && !pad_warned && ccs_pad !== :auto
                xmax = ccs_max_x(cmom, p)
                if xmax > pad
                    @printf("  ⚠ CCS: max |x| = %.4f exceeds the fixed pad %.4f (ccs_pad): the draws' laws leave the grid\n", xmax, pad)
                    pad_warned = true
                end
            end
        end
        gr = ccs ? build_grids(p; H̃T = HT, M, t, la, grid_kw..., ccs = (; lQT, c = cmom), pad) :
                   build_grids(p; H̃T = HT, M, t, la, grid_kw...)
        ccs && (hh0 = ccs_regrid_warm(hh0, gr_prev, gr); gr_prev = gr)
        t_hh = @elapsed (hh0 = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit,
                                                verbose = verbose && it == 1, hh0))
        s_relaxed = max(s_relaxed, hh0.s_relaxed)
        t_pi  = @elapsed ((πT, πO) = location_probs(hh0, p, gr))
        # §3d: the young are the parents' children (the adults' πT, πO without the fertility shock)
        t_phi = @elapsed ((Φ, _) = stationary_phi(parent_probs(hh0, p, gr, πT, πO)..., hh0, p, gr))
        t_agg = @elapsed ((HTn, Mn, tn) = aggregates(Φ, πT, πO, hh0, p, gr))

        # An empty location has no meaningful budget: hold its tax fixed.
        live = Mn .≥ M_FLOOR
        tn[.!live] .= t[.!live]
        t_err = any(live) ? maximum(abs, tn[live] .- t[live]) : 0.0
        err = ccs ? max(reldiff(Mn, M), t_err) : max(reldiff(HTn, HT), reldiff(Mn, M), t_err)   # §3e: H̃_T is not a state
        if peer                                   # the peer input log P (A: the students' mean ability), in absolute change
            lan = peer_log_input(Φ, πT, πO, hh0, p, gr)
            lan[.!live] .= la[.!live]             # nor peers: hold an empty location's log P too
            err = max(err, maximum(abs, lan .- la))
            @. la = (1 - damping) * la + damping * lan
        end
        if ccs          # §3e: ln Q̄ and the teacher moments, in absolute change of ln Q̄ and of the rules (before the closure's gate)
            tm   = ccs_teacher_moments(Φ, πT, πO, hh0, p, gr; c0 = cμ)
            imp  = ccs_implied(tm, Φ, live, lQT, cmom, cμ, p)
            err  = max(err, maximum(abs, imp.lQT .- lQT), ccs_rule_change(gr.ccs, imp.c, p))
            @. lQT = (1 - damping) * lQT + damping * imp.lQT
            cmom = ccs_damp(cmom, imp.c, damping)   # the mixture's moments: stays PSD
            cμ   = imp.μ
        end

        fte = NaN
        if closure                                # positions per pupil: b rations the suburban positions
            fte = fte_gap(Φ, πT, πO, hh0, p, gr)
            res = p.fte_rule - fte
            if err < FTE_SETTLE                   # the other aggregates have settled: move b
                if max(err, e_prev) < FTE_SECANT_ERR && abs(b - b_prev) > FTE_SECANT_DB
                    sec = (fte - f_prev) / (b - b_prev)
                    isfinite(sec) && sec > 0 && (slope = clamp(sec, FTE_SLOPE_MIN, FTE_SLOPE_MAX))
                end
                isfinite(r_prev) && res * r_prev < 0 && abs(res) > abs(r_prev) &&
                    (slope = min(2 * slope, FTE_SLOPE_MAX))   # a growing overshoot: halve the step
                b_prev, f_prev, e_prev, r_prev = b, fte, err, res
                b  += clamp(damping * res / slope, -FTE_STEP_MAX, FTE_STEP_MAX)
                p   = with_params(p; BT = [p.BT[1], b])
            end
            err = max(err, abs(res))
        end

        @. HT = (1 - damping) * HT + damping * HTn
        @. M  = (1 - damping) * M  + damping * Mn
        @. t  = (1 - damping) * t  + damping * tn

        if verbose
            @printf("GE %3d  err=%.3e  H̃_T=%s  M=%s  t=%s%s%s  [hh=%.1fs pi=%.1fs phi=%.1fs agg=%.1fs]\n",
                    it, err, fmtvec(HT), fmtvec(M), fmtvec(t), peer ? "  log P=" * fmtvec(la) : "",
                    closure ? @sprintf("  FTE gap=%.5f b=%.5f (slope %.2f)", fte, b, slope) : "",
                    t_hh, t_pi, t_phi, t_agg)
            ccs && @printf("        CCS: ln Q̄=%s  SD(ln r)=%s  nodes=%s  pad=%.2f (nϵT=%d, nXO=%d)\n", fmtvec(lQT),
                           fmtvec(ccs_sds(gr)), "(" * join(length.(gr.ccs.x), ", ") * ")", pad, gr.nϵT, gr.nXO)
            flush(stdout)
        end
        if err < tol
            converged = true
            verbose && println("GE aggregates converged.")
            break
        end
    end
    verbose && !converged && println("  ⚠ GE did not converge in $maxit iterations.")

    ccs && ccs_pad === :auto && (pad = max(pad, ccs_auto_pad(cmom, p)))
    gr = ccs ? build_grids(p; H̃T = HT, M, t, la, grid_kw..., ccs = (; lQT, c = cmom), pad) :
               build_grids(p; H̃T = HT, M, t, la, grid_kw...)
    ccs && (hh0 = ccs_regrid_warm(hh0, gr_prev, gr))
    hh = solve_household(p, gr; tol = hh_tol, maxit = hh_maxit, hh0)
    ccs && return ccs_solution(hh, p, gr, HT, M, t, la, converged, max(s_relaxed, hh.s_relaxed), lQT, cmom, cμ, pad)
    p.split && return (; hh, split_probs(hh, p, gr)..., HT, M, t, la, gr, p, converged,
                         s_resid = hh.s_resid, s_passes = hh.s_passes, s_relaxed = max(s_relaxed, hh.s_relaxed))
    Pi_T, Pi_O = location_probs(hh, p, gr)
    Pi_Tk, Pi_Ok = parent_probs(hh, p, gr, Pi_T, Pi_O)          # §3d: the adults' own without the fertility shock
    Φ, πbar_kid = stationary_phi(Pi_Tk, Pi_Ok, hh, p, gr)
    πbar = fertile(p) ? migration_probs(Pi_T, Pi_O, hh, p, gr) : πbar_kid
    return (; hh, Pi_T, Pi_O, Φ, πbar, HT, M, t, la, gr, p, converged,
              s_resid = hh.s_resid, s_passes = hh.s_passes, s_relaxed = max(s_relaxed, hh.s_relaxed),
              Pi_Tk, Pi_Ok, πbar_kid)     # §3d: with pkid_z ≠ 0 these three carry P(child | z)/pkid (`parent_probs`)
end

# -----------------------------------------------------------------------------
# 9. Diagnostics and summary statistics
# -----------------------------------------------------------------------------
"""
    verify_solution(sol; δ = C_FLOOR, tclamp = TAX_CAP, s_atol = 1e-8)

Check that the numerical safeguards are slack at the solution: (i) consumption
exceeds δ in every (node, work location), so smooth_log equals log (under the split, in every
(node, residence, workplace)); (ii) no tax
rate is on the clamp [0, tclamp]; (iii) the (e, s) fixed point converged. It also
reports how many nodes needed the relaxed (e, s) update (`solve_node`) in the largest
sweep, which is not a failure. Returns true if all hold. Condition (i) is what limits how large `mcost` can be:
the goods cost bites first on low-ability movers. Under constant class size (§3e) it also checks that the grid pad
covers every node of the draw rules, max |x| ≤ pad (otherwise a draw's law leaves the effective-shock grid), and
reports it; the pad nodes are in (i) and (iii) like any other.
"""
function verify_solution(sol; δ = C_FLOOR, tclamp = TAX_CAP, s_atol = 1e-8)
    (; hh, t, gr, p, s_resid, s_passes, s_relaxed) = sol
    minC = Inf
    n_binding = 0
    for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz, d in branch_ids(p)
        hb = branch_of(hh, p, d)
        for k in 1:gr.nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p)
            minC = min(minC, minimum(C)); n_binding += count(≤(δ), C)
        end
        for k in 1:gr.nXO
            nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sO[g, l, zi, k], hb.eO[g, l, zi, k], p)
            minC = min(minC, minimum(C)); n_binding += count(≤(δ), C)
        end
    end
    ok = true
    if s_resid > s_atol
        @printf("  ⚠ s fixed point residual %.2e (> %.1e) after %d passes\n", s_resid, s_atol, s_passes)
        ok = false
    end
    @printf("  (e, s) update: relaxed step needed at up to %d of %d nodes in a sweep\n", s_relaxed,
            length(branch_ids(p)) * 2 * gr.L * gr.Nz * (gr.nϵT + gr.nXO))
    if n_binding > 0
        @printf("  ⚠ smooth_log barrier BINDING at %d (node, l') cells (min C = %.3e ≤ δ = %.1e)\n",
                n_binding, minC, δ)
        ok = false
    end
    for lp in 1:gr.L
        if t[lp] >= tclamp - 1e-9 || t[lp] <= 0.0
            @printf("  ⚠ tax rate t[%d] = %.4f is on the clamp [0, %.1f]: budget not balanced\n",
                    lp, t[lp], tclamp)
            ok = false
        end
    end
    if ccs_on(gr)                                 # §3e: the pad must cover the draws
        xmax = maximum(maximum(abs, x) for x in gr.ccs.x)
        if xmax > gr.ccs.pad + 1e-12
            @printf("  ⚠ CCS: max |x| = %.4f of the draw rules exceeds the grid pad %.4f: the shifted laws leave the grid\n",
                    xmax, gr.ccs.pad)
            ok = false
        else
            @printf("  CCS: draw rules within the grid pad, max |x| = %.4f ≤ pad = %.4f (nϵT = %d, nXO = %d)\n",
                    xmax, gr.ccs.pad, gr.nϵT, gr.nXO)
        end
    end
    ok && @printf("  ✓ feasibility audit passed: min C = %.4e > δ = %.1e; taxes interior t = %s; s residual %.1e (≤ %d passes)\n",
                  minC, δ, fmtvec(t), s_resid, s_passes)
    return ok
end

"""
    sol_probs(sol) -> (πT, πO)

The work-location probabilities of a solution as `cell_branches` takes them: `Pi_T`, `Pi_O`, or the
joint residence × workplace probabilities `Pi_TW`, `Pi_OW` under the split (a pair of arrays, one per
education branch, under `college`).
"""
sol_probs(sol) = sol.p.split ? (sol.Pi_TW, sol.Pi_OW) : (sol.Pi_T, sol.Pi_O)

"""
    cell_mean(sol, qT, qO)

Population mean of a quantity given on the grids, weighted by Φ and integrated
over the occupation choice. `qT(br, g, l, zi)` and `qO(br, g, l, zi)` return the vectors
on the ϵ_T and X_O* grids for that cell and education branch `br` (a `cell_branches` entry:
`br.h` the household arrays, `br.πT`, `br.πO` the location probabilities (the joint ones under the
split), `br.d` the label); under `p.college` the branches are mixed with the weights
P(d | g, l, z) (§3c).
With `law = true` (§3e, 2026-10-26) the quantity may depend on the teacher draw: the closures are called as
`qT(br, g, l, zi, law, x)`, `qO(…, law, x)` for every draw (x = ln r, `law` its `ccs_law`, to be passed wherever an
integral or `argmax_mean` takes `gr`) and the cell's integral is their R_l-mixture (`over_draws`); without CCS they are
called once with (`gr`, 0.0). With `law = false` (the default) the closures take (br, g, l, zi) as always, and the
quantity must not depend on the draw.
"""
function cell_mean(sol, qT, qO; law::Bool = false)
    (; hh, Φ, gr, p) = sol
    PT, PO = sol_probs(sol)
    onesT, onesO = ones(gr.nϵT), ones(gr.nXO)
    num = den = 0.0
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2, br in cell_branches(hh, PT, PO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        if law
            v = over_draws(cm, gr) do cmk, lw, x
                sum(integrate_choice(qT(br, g, l, zi, lw, x), qO(br, g, l, zi, lw, x), cmk, lw))
            end
            num += mass * v
        else
            num += mass * sum(integrate_choice(qT(br, g, l, zi), qO(br, g, l, zi), cm, gr))
        end
        den += mass * sum(integrate_choice(onesT, onesO, cm, gr))
    end
    return num / den
end

"""
    mean_predicted_log_income(sol)

Population mean of the predicted log income ln ŷ (`predicted_log_income`), the grading variable of the
income-graded location taste `Params.By` (round 2026-10-15c), over every node weighted by Φ and the occupation
choice (`cell_mean`). The calibration freezes the reference `lyref` at this mean at a reference solution; ln ŷ does
not depend on `By`, so the mean is the same with the taste on or off at the same aggregates.
"""
function mean_predicted_log_income(sol)
    (; gr, p) = sol
    qT(br, g, l, zi) = [predicted_log_income(teach_node(gr.ϵTgrid[k], zi, l, g, br.h.Λ, p, gr; d = br.d), p, gr)
                        for k in 1:gr.nϵT]
    qO(br, g, l, zi) = [predicted_log_income(nonteach_node(gr.XOgrid[k, g], zi, l, g, br.h.Λ, p, gr; d = br.d), p, gr)
                        for k in 1:gr.nXO]
    return cell_mean(sol, qT, qO)
end

"""
    mean_child_h(sol)

Mean human capital h̄ at `sol`. The warm-glow reference `href` is frozen at this
value (`HBAR_BENCH` at `nesting_params()`) so that ψ changes only the curvature of
the kernel, not its level. Non-teachers' capacity is net of the credential premium (§3c).
Under constant class size (§3e) E[Θ_{i*}^{−1} | X̃_O*] is taken under each draw's law (`cell_mean` with `law = true`).
"""
function mean_child_h(sol)
    (; hh, gr, p) = sol
    hO(br, g, l, zi, lw = gr, x = 0.0) = [br.h.capO[g, l, zi, k] / credential(p, br.d) * EΘpow(lw, g, lw.XOgrid[k, g], 0.0)
                                         for k in 1:lw.nXO]
    ccs_on(gr) && return cell_mean(sol, (br, g, l, zi, _...) -> br.h.hT[g, l, zi, :], hO; law = true)
    return cell_mean(sol, (br, g, l, zi) -> br.h.hT[g, l, zi, :], hO)
end

"""
    mean_consumption(sol)

Mean consumption C̄, averaged over location choices, at `sol` (`CBAR_BENCH` at
`nesting_params()`). Used by `redenominate_move_cost`. It is net of the subsistence
c̲ = `csub` C̄, like `consumption`; under the split it averages over the (residence, workplace)
pairs.
"""
function mean_consumption(sol)
    sol.p.split && return mean_consumption_split(sol)
    (; gr, p) = sol
    CT(br, g, l, zi) = map(1:gr.nϵT) do k
        nd = teach_node(gr.ϵTgrid[k], zi, l, g, br.h.Λ, p, gr; d = br.d)
        dot(br.πT[g, l, zi, k, :], consumption(nd, br.h.sT[g, l, zi, k], br.h.eT[g, l, zi, k], p))
    end
    CO(br, g, l, zi) = map(1:gr.nXO) do k
        nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, br.h.Λ, p, gr; d = br.d)
        dot(br.πO[g, l, zi, k, :], consumption(nd, br.h.sO[g, l, zi, k], br.h.eO[g, l, zi, k], p))
    end
    return cell_mean(sol, CT, CO)
end

function mean_consumption_split(sol)
    (; gr, p) = sol
    CT(br, g, l, zi) = map(1:gr.nϵT) do k
        nd = teach_node(gr.ϵTgrid[k], zi, l, g, br.h.Λ, p, gr; d = br.d)
        dot(pair_probs(br.πT, g, l, zi, k), consumption(nd, br.h.sT[g, l, zi, k], br.h.eT[g, l, zi, k], p))
    end
    CO(br, g, l, zi) = map(1:gr.nXO) do k
        nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, br.h.Λ, p, gr; d = br.d)
        dot(pair_probs(br.πO, g, l, zi, k), consumption(nd, br.h.sO[g, l, zi, k], br.h.eO[g, l, zi, k], p))
    end
    return cell_mean(sol, CT, CO)
end

"""
Mean time share among teachers and non-teachers, and its (min, max) over all nodes. Under
`p.college` also the means by education: BA and non-BA teachers (`meanT_ba`, `meanT_nonba`), BA and
non-BA non-teachers (`meanO_ba`, `meanO_nonba`); NaN without college.
"""
function s_summary(sol)
    (; hh, Φ, gr, p) = sol
    numT = numO = denT = denO = 0.0
    num = zeros(2, 2); den = zeros(2, 2)           # [teacher / non-teacher, d = 1 / d = 0]
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        sT, sO = integrate_choice(br.h.sT[g, l, zi, :], br.h.sO[g, l, zi, :], cm, gr)
        mT, mO = integrate_choice(ones(gr.nϵT), ones(gr.nXO), cm, gr)
        numT += mass * sT; denT += mass * mT
        numO += mass * sO; denO += mass * mO
        j = br.d == 1 ? 1 : 2
        num[1, j] += mass * sT; den[1, j] += mass * mT
        num[2, j] += mass * sO; den[2, j] += mass * mO
    end
    ext(name) = extrema(vcat((vec(getproperty(branch_of(hh, p, d), name)) for d in branch_ids(p))...))
    return (; meanT = numT / denT, meanO = numO / denO, rangeT = ext(:sT), rangeO = ext(:sO),
              meanT_ba = p.college ? num[1, 1] / den[1, 1] : NaN,
              meanT_nonba = p.college ? num[1, 2] / den[1, 2] : NaN,
              meanO_ba = p.college ? num[2, 1] / den[2, 1] : NaN,
              meanO_nonba = p.college ? num[2, 2] / den[2, 2] : NaN)
end

"The BA share P₁ of the young (Φ-weighted over cells, both genders), overall and by ability node: `p.college` only."
function ba_share_young(sol)
    (; hh, Φ, gr) = sol
    w = [hh.P1[g, l, zi] * Φ[zi, l] / 2 for g in 1:2, l in 1:gr.L, zi in 1:gr.Nz]
    return sum(w) / sum(Φ), [sum(w[:, :, zi]) / sum(Φ[zi, :]) for zi in 1:gr.Nz]
end

"Teaching share among the young born in `l` (both genders)."
function teaching_share_endo(hh, l, Φ, gr, p::Params)
    p.college && return college_teaching_share(hh, l, Φ, gr, p)
    mass = sum(@view Φ[:, l])
    return sum(Φ[zi, l] / mass / 2 * integrate_teach(ones(gr.nϵT), choice_maps(hh, g, l, zi, gr), gr)
               for zi in 1:gr.Nz, g in 1:2)
end

"""
    student_ability_moments(sol)

Ability moments of the young by birth location, split by the occupation chosen:
- `z`: persistent ability. Its mean equals the `mean z` row of `report_ge`; the split
  shows which types select into teaching.
- `a = z·ϵ`: ability in the chosen occupation (teachers: ϵ_T; others: ϵ_{i*} via
  `Eϵ_nonteach`). E[a] exceeds E[z] because both branches select on their own ϵ.

Returns per-location vectors (`mean_a`, `mean_a_teach`, `mean_a_non`, `mean_z`,
`mean_z_teach`, `mean_z_non`, `teach_share`, `mass`, `massT`, `massO`), the
branch split of Φ (`ΦT`, `ΦO`), and population-pooled scalars suffixed `_all`.
"""
function student_ability_moments(sol)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, L, nϵT, nXO) = gr
    onesT, onesO = ones(nϵT), ones(nXO)
    ΦT = zeros(Nz, L); ΦO = zeros(Nz, L)
    aT = zeros(L);     aO = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm  = choice_maps(br.h, g, l, zi, gr)
        if cm isa MixedMaps                      # §3e: the true shocks, ϵ_T = ϵ̃_T e^{−x/α} and ϵ_{i*} = ϵ̃_{i*} e^{−x/α}
            mT, mO = integrate_choice(onesT, onesO, cm, gr)
            IaT, IaO = over_draws(cm, gr) do cmk, law, x
                sc = exp(-x / p.α)
                integrate_choice(z[zi] .* law.ϵTgrid .* sc,
                                 [z[zi] * Eϵ_nonteach(law, g, law.XOgrid[k, g], p.α) * sc for k in 1:nXO], cmk, law)
            end
            ΦT[zi, l] += mass * mT; ΦO[zi, l] += mass * mO
            aT[l]     += mass * IaT; aO[l]    += mass * IaO
            continue
        end
        aTv = z[zi] .* gr.ϵTgrid
        aOv = [z[zi] * Eϵ_nonteach(gr, g, gr.XOgrid[k, g], p.α) for k in 1:nXO]
        mT, mO = integrate_choice(onesT, onesO, cm, gr)
        IaT, IaO = integrate_choice(aTv, aOv, cm, gr)
        ΦT[zi, l] += mass * mT; ΦO[zi, l] += mass * mO
        aT[l]     += mass * IaT; aO[l]    += mass * IaO
    end
    ratio(num, den) = [den[l] > 0 ? num[l] / den[l] : NaN for l in 1:L]
    massT, massO = vec(sum(ΦT; dims = 1)), vec(sum(ΦO; dims = 1))
    mass = massT .+ massO
    zT, zO = vec(z' * ΦT), vec(z' * ΦO)
    return (; mass, massT, massO, ΦT, ΦO,
              teach_share  = ratio(massT, mass),
              mean_a       = ratio(aT .+ aO, mass),
              mean_a_teach = ratio(aT, massT),
              mean_a_non   = ratio(aO, massO),
              mean_z       = ratio(zT .+ zO, mass),
              mean_z_teach = ratio(zT, massT),
              mean_z_non   = ratio(zO, massO),
              teach_share_all  = sum(massT) / sum(mass),
              mean_a_all       = sum(aT .+ aO) / sum(mass),
              mean_a_teach_all = sum(aT) / sum(massT),
              mean_a_non_all   = sum(aO) / sum(massO),
              mean_z_all       = sum(zT .+ zO) / sum(mass),
              mean_z_teach_all = sum(zT) / sum(massT),
              mean_z_non_all   = sum(zO) / sum(massO))
end

"""
    student_ability_density(sol, l, agrid) -> (fT, fO)

Density of student ability a = z·ϵ over `agrid` for those born in `l`, split by
branch and normalised so ∫(fT + fO) = 1 and ∫fT is the teaching share. Evaluated
in closed form, which avoids the comb artefacts of binning quadrature nodes:
- teaching: f_ϵ(a/z)/z · P(teach | ϵ_T = a/z);
- non-teaching: Σ_i f_ϵ(ϵ) ∏_{j≠i} F_j(Θ_i ϵ^α) · P(don't teach | X_O* = Θ_i ϵ^α),
  where the product is the probability that occupation i attains the max.
Under constant class size (§3e) a and the densities are those of the true shocks (the laws at r = 1, `gr.ccs.base`), and
the occupation weights the R_l-mixture of each draw's, at the effective shocks ϵ e^{x/α} and X_O e^{x} (`over_draws`).
"""
function student_ability_density(sol, l, agrid)
    (; hh, Φ, gr, p) = sol
    ccs_on(gr) && return ccs_student_ability_density(sol, l, agrid)
    (; z, Nz, dT, dO, nonteach, logΘ) = gr
    fT = zeros(length(agrid)); fO = zeros(length(agrid))
    for zi in 1:Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        for (ia, a) in enumerate(agrid)
            ϵ  = a / z[zi]
            fϵ = mass * pdf(dT, ϵ) / z[zi]
            fT[ia] += fϵ * teach_wt(cm, ϵ)
            for i in nonteach
                x = exp(logΘ[i, g]) * ϵ^p.α
                wins = prod((cdf(dO[j, g], x) for j in nonteach if j != i); init = 1.0)
                fO[ia] += fϵ * wins * nonteach_wt(cm, x)
            end
        end
    end
    mass = sum(@view Φ[:, l])
    return fT ./ mass, fO ./ mass
end

"`student_ability_density` under constant class size (§3e): true-shock densities, draw-mixed occupation weights."
function ccs_student_ability_density(sol, l, agrid)
    (; hh, Φ, gr, p) = sol
    (; z, Nz, nonteach, logΘ) = gr
    (; dT, dO) = gr.ccs.base                    # the laws of the true shocks (r = 1)
    fT = zeros(length(agrid)); fO = zeros(length(agrid))
    for zi in 1:Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        for (ia, a) in enumerate(agrid)
            ϵ  = a / z[zi]
            fϵ = mass * pdf(dT, ϵ) / z[zi]
            xs = [exp(logΘ[i, g]) * ϵ^p.α for i in nonteach]
            ws = [prod((cdf(dO[j, g], xs[ii]) for j in nonteach if j != i); init = 1.0) for (ii, i) in enumerate(nonteach)]
            tw, ow = over_draws(cm, gr) do cmk, _, x
                (teach_wt(cmk, ϵ * exp(x / p.α)), sum(ws[ii] * nonteach_wt(cmk, xs[ii] * exp(x)) for ii in eachindex(xs)))
            end
            fT[ia] += fϵ * tw
            fO[ia] += fϵ * ow
        end
    end
    mass = sum(@view Φ[:, l])
    return fT ./ mass, fO ./ mass
end

function report_ge(sol)
    (; hh, Φ, πbar, HT, M, t, gr, p) = sol
    (; z, Nz, Q, Πz, L) = gr

    println("\n===== GE stationary solution (ε = 0, continuous distributions) =====")
    ss = s_summary(sol)
    @printf("  s (choice-weighted mean [min, max]):  T=%.4f [%.4f, %.4f]   nonteach=%.4f [%.4f, %.4f]\n",
            ss.meanT, ss.rangeT..., ss.meanO, ss.rangeO...)
    if p.college
        ba, ba_z = ba_share_young(sol)
        @printf("  s by education: BA teachers %.4f, non-BA teachers %.4f, BA non-teachers %.4f, non-BA non-teachers %.4f;  BA share of the young %.4f, by z %s\n",
                ss.meanT_ba, ss.meanT_nonba, ss.meanO_ba, ss.meanO_nonba, ba, fmtvec(ba_z))
    end
    @printf("  H̃_T    = %s\n", fmtvec(HT))
    @printf("  M      = %s      [total %.4f]\n", fmtvec(M), sum(M))
    @printf("  t      = %s\n", fmtvec(t))
    @printf("  Q      = %s\n", fmtvec(Q))
    if ccs_on(gr) && hasproperty(sol, :ccs)     # §3e: constant class size
        cc   = sol.ccs
        Qeff = exp.(cc.lQeff)                       # (2H̃_T/M)^σ at the solution's teachers
        @printf("  CCS    : n_q = %d nodes, K_d = %d teacher draws, ρ = %.3f, pad = %.2f;  teachers n_T = %s, μ = E[ln h_eff | l] = %s\n",
                cc.nq, cc.draws, cc.rho, cc.pad, fmtvec(cc.nT), fmtvec(cc.μ))
        @printf("           class size N_l = S_l/n_T,l = %s;  Q̄_l = %s (ln Q̄ = %s; implied now %s);  SD(ln q | l) = %s\n",
                fmtvec(cc.N), fmtvec(gr.QT), fmtvec(cc.lQT), fmtvec(cc.lQT_now), fmtvec(cc.sd))
        @printf("           efficient-equivalent Q^eff = (2H̃_T/M)^σ = %s;  ln Q^eff − ln Q̄ = %s\n",
                fmtvec(Qeff), fmtvec(cc.lQeff .- cc.lQT))
        for l in 1:L
            @printf("           R_%d: x = %s  w = %s\n", l, fmtvec(cc.x[l]), fmtvec(cc.w[l]))
        end
    end
    if !iszero(p.ψa)         # 2026-10-25: the specification, and log P at Φ split into its mean and scale terms
        PT, PO = sol_probs(sol)
        ad = adult_log_h(Φ, PT, PO, hh, p, gr)
        pt = peer_terms(Φ, PT, PO, hh, p, gr; adults = ad)
        @printf("  peers  : ψa = %.3f, specification %s (the %s);  log P = %s (at Φ: %s = mean %s + scale %s);  teacher part Q_T = %s\n",
                p.ψa, peer_spec(p), peer_spec_name(p), fmtvec(gr.la), fmtvec(peer_log_input(Φ, PT, PO, hh, p, gr; adults = ad)),
                fmtvec(pt.mean), fmtvec(p.peer_total ? pt.scale : zeros(L)), fmtvec(gr.QT))
        @printf("           students' log ā = %s;  adult residents' mean log h − economy mean = %s (adults %s)\n",
                fmtvec(peer_log_ability(Φ, gr)), fmtvec(ad.H .- ad.Hbar), fmtvec(ad.m))
    end

    ergodic = stationary(Πz)
    println("\n  Endogenous ability distribution  G_l(z)   (vs ergodic G*):")
    @printf("    %10s", "z")
    for l in 1:L; @printf(" %10s", "G_$(l)(z)"); end
    @printf(" %10s\n", "G*(z)")
    for zi in 1:Nz
        @printf("    %10.4f", z[zi])
        for l in 1:L; @printf(" %10.4f", Φ[zi, l] / sum(@view Φ[:, l])); end
        @printf(" %10.4f\n", ergodic[zi])
    end
    mean_z(w) = dot(z, w) / sum(w)
    print("    mean z: ")
    for l in 1:L; @printf(" loc %d = %.4f  ", l, mean_z(@view Φ[:, l])); end
    @printf("ergodic = %.4f\n", mean_z(ergodic))

    println("\n  Migration matrix  Π̄[born→work]  (rows sum to 1):")
    for l in 1:L
        row = [dot(Φ[:, l], πbar[l, :, lp]) / sum(@view Φ[:, l]) for lp in 1:L]
        print("    born $l :")
        for lp in 1:L; @printf("  ->%d %.4f", lp, row[lp]); end
        println()
    end
    if fertile(p)            # §3d: the parents' rows (where the next generation is born) and the adults by location
        (; πbar_kid) = sol
        @printf("  parents (P(child) = %.4f)%s:\n", p.pkid,
                iszero(p.pkid_z) ? "" : @sprintf(", children per young (P(child | z) = %s, slope %.3f)", fmtvec(gr.pk), p.pkid_z))
        for l in 1:L
            row = [dot(Φ[:, l], πbar_kid[l, :, lp]) / sum(@view Φ[:, l]) for lp in 1:L]
            print("    born $l :")
            for lp in 1:L; @printf("  ->%d %.4f", lp, row[lp]); end
            println()
        end
        adults = [sum(Φ[zi, l] * πbar[l, zi, lp] for l in 1:L, zi in 1:Nz) for lp in 1:L]
        @printf("  adults by location %s (shares %s) vs the young %s\n", fmtvec(2 .* adults),
                fmtvec(adults ./ sum(adults)), fmtvec(vec(sum(Φ; dims = 1)) ./ sum(Φ)))
    end

    println("\n  Teaching share among young born in l (endogenous G_l):")
    for l in 1:L
        @printf("    born %d :  %.4f\n", l, teaching_share_endo(hh, l, Φ, gr, p))
    end

    sa = student_ability_moments(sol)
    println("\n  Student ability by birth location  (choice-weighted; a = z·ϵ in the chosen occupation):")
    @printf("    %-6s %-8s %-9s %-9s %-11s %-13s %-11s %-13s\n",
            "born", "teach", "mean z", "mean a", "a | teach", "a | nonteach",
            "z | teach", "z | nonteach")
    for l in 1:L
        @printf("    %-6d %-8.4f %-9.4f %-9.4f %-11.4f %-13.4f %-11.4f %-13.4f\n",
                l, sa.teach_share[l], sa.mean_z[l], sa.mean_a[l],
                sa.mean_a_teach[l], sa.mean_a_non[l],
                sa.mean_z_teach[l], sa.mean_z_non[l])
    end
    @printf("    %-6s %-8.4f %-9.4f %-9.4f %-11.4f %-13.4f %-11.4f %-13.4f\n",
            "all", sa.teach_share_all, sa.mean_z_all, sa.mean_a_all,
            sa.mean_a_teach_all, sa.mean_a_non_all,
            sa.mean_z_teach_all, sa.mean_z_non_all)
    println("====================================================================\n")
end

# -----------------------------------------------------------------------------
# 10. Run
# -----------------------------------------------------------------------------
function main()
    println("Solving spatial model (ε = 0, continuous) in general equilibrium ...")
    sol = solve_ge(Params(); Nz = 3, nϵT = 48, nXO = 48, damping = 0.75,
                   tol = 1e-6, hh_tol = 1e-6, maxit = 500, hh_maxit = 1000)
    report_ge(sol)
    verify_solution(sol)
    return sol
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end
