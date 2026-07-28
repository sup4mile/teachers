# Spatial Model with Altruism

This document presents the spatial occupational-choice / teacher-spillover economy under the assumption that the parental cost share is zero. Parents still choose the location in which their children are born and educated, and still care about their children's human capital through a warm-glow term, but they finance none of the children's educational expenditures.

Relative to the previous iteration, imposing no parental cost sharing ($\varepsilon = 0$) from the beginning buys three structural simplifications, each developed below:

1. **Consumption decouples from the child**. Old-age consumption is deterministic given the agent's own state, so the location value separates additively into a consumption term and a warm-glow term. Altruism shapes behavior only through location choice.
2. **The time investment is closed form**. The optimal time share $s$ is a constant per occupation — independent of ability, birth location, destination wedges, and the altruism strength. Under the sorting extension of §5.1 this weakens only to a *conditional* closed form (Proposition 2′).
3. **Non-teaching occupations collapse to a scalar**. The best non-teaching option is summarized by a one-dimensional sufficient statistic with known distribution, and the teach/not-teach margin is a single-crossing threshold. The $I$-dimensional choice problem becomes two-dimensional.

### Two sorting extensions

In its baseline form the model generates essentially **no sorting of ability across space**, and for a structural rather than a calibration reason: location-choice probabilities are exactly invariant to the agent's ability draws (Proposition 4 below). Two minimal departures break that invariance, each chosen to introduce no new free parameter and to leave the tractability above intact:

- **§5.1 — Goods-denominated moving cost.** The utility moving cost $\tau_{l,l'}$ is re-denominated into a goods cost $m_{l,l'}$ paid out of old-age consumption. A fixed goods cost is a proportionally larger burden on the low-ability, so richer agents are less deterred by distance and concentrate wherever the pull is. Proposition 3 survives untouched; Proposition 2 weakens to Proposition 2′.
- **§5.2 — Power warm-glow kernel.** The kernel $f(h') = \log h'$ is replaced by $f(h') = \big[(h'/\bar h)^{1-\psi} - 1\big]/(1-\psi)$ with a *fixed* reference level $\bar h$. Better schools are then worth more in utils to parents with better children, so the move-to-opportunity motive becomes ability-graded. Nothing structural breaks: Propositions 1–3 hold verbatim.

Both extensions nest the baseline exactly: $m_{l,l'} \equiv 0$ (with $\tau_{l,l'}$ restored) and $\psi = 1$ return every result in its original form. Throughout, statements that hold only in the nested baseline are flagged.

## 1 Environment

There are $L=2$ locations indexed by $l \in \{1,2\}$. In each location the same set of $I$ occupations is available, with $i=T$ denoting teaching. Locations share common occupational productivities $\{A_i\}_{i\neq T}$, labor-market distortions $\{\tau^{\omega}_{i,g}\}$, and educational barriers $\{\tau^{e}_{i,g}\}$. Locations differ along:

- the (endogenous) stock of teacher human capital $\widetilde{H}_{T,l}$,
- the (endogenous) population or mass of agents in each location $M_l$
- the local labor-income tax $t_l$, which funds teachers under a local balanced-budget rule,
- an exogenous amenity value $B_l$,
- the teaching wage intercept $\kappa_l$.

A young agent born in $l$ who later works in $l'$ pays a moving cost. In the baseline this is a utility cost $\tau_{l,l'} \geq 0$ with $\tau_{l,l}=0$; under the §5.1 extension it is instead a **goods** cost $m_{l,l'} \geq 0$ with $m_{l,l}=0$, paid out of old-age consumption. The two are set at parity through the consumption-equivalent
$$
m_{l,l'} \;=\; \big(1 - e^{-\tau_{l,l'}/\mu}\big)\,\bar C ,
$$
where $\bar C$ is a **fixed** reference consumption level (benchmark mean $C$), so that at the benchmark the goods cost is worth exactly the utils the original cost was worth. This is a re-denomination, not a new parameter.

The economy is populated by overlapping generations of finite-lived agents (young and old) with constant total measure $M$. Population in each location $l$ has equal cohort sizes $M_l/2$ within each location in stationary equilibrium; the split $\{M_l\}$ is pinned down by location flows. Groups $g \in \{m,f\}$ denote gender. A parent's child is male or female with equal probability, independently of the parent's gender, so young agents are evenly split by gender in every location. Because distortions $\tau^\omega_{i,g}$, barriers $\tau^e_{i,g}$, and hence policies are gender-specific, the coin flip over the child's gender is payoff-relevant to the parent, and the gender mix among *old* (working-age) agents in each location is endogenous — men and women sort differently across locations.

## 2 Timing

Each cohort's life cycle proceeds in three stages.

**Stage 1 — Young.** Born in location $l$ (chosen by the parent), with persistent ability $z$ inherited from the parent. The vector of occupation-specific idiosyncratic shocks $\vec\epsilon = (\epsilon_1,\ldots,\epsilon_I)$ is drawn iid from $F_\epsilon$, giving ability $a_i = z\,\epsilon_i$ in each occupation. The agent chooses time investment $s$, goods investment $e$, and an occupation $i$. The full education cost $(1+\tau^{e}_{i,g})\,e$ is borne by the agent herself, paid out of future labor income.

**Stage 2 — Transition.** The agent reaches working age, draws iid Gumbel taste shocks $\{\nu_{l'}\}_{l'=1}^{L}$, and chooses a work location $l'$ in which to earn income, pay local taxes, and raise her own child.

**Stage 3 — Old.** She earns occupation-specific labor income in $l'$, pays the local tax, her education debt, and (under §5.1) the goods moving cost $m_{l,l'}$, and consumes the residual. Her child's persistent ability $z' \sim \pi_z(\cdot\mid z)$, idiosyncratic vector $\vec\epsilon'$, and gender $g'$ realize, and she receives a warm-glow payoff $\lambda f(h')$ from the child's accumulated human capital. There is no financial transfer between generations.

## 3 Ability Transmission

Ability in occupation $i$ is the product $a_i = z\,\epsilon_i$ of a persistent scalar and an occupation-specific idiosyncratic shock.

**Persistent component.** The scalar $z$ is the only ability state inherited across generations, following a log-AR(1):
$$
\log z' \;=\; \rho_z \log z \;+\; \sigma_{\xi}\,\xi, \qquad \xi \sim N(0,1)\ \text{iid},\quad |\rho_z|<1,
$$
with transition density $\pi_z(z'\mid z)$.

**Idiosyncratic component.** The shocks $\{\epsilon_i\}_{i=1}^I$ are iid across occupations and generations, lognormal with unit mean:
$$
\log \epsilon_i \;\sim\; N\!\left(-\tfrac{1}{2}\sigma_\epsilon^2,\; \sigma_\epsilon^2\right).
$$
The vector $\vec\epsilon$ is not inherited; the child draws a fresh $\vec\epsilon'$ independent of the parent's $\vec\epsilon$ and of $z'$. Both components are continuously distributed, and this continuity does real work below: it makes every aggregate a smooth function of the occupational threshold.

**Information.** The young agent observes $(z,\vec\epsilon)$ before choosing. The parent's expectation over her child factors as $\mathbb E_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}$, all mutually independent.

**Endogenous spatial distribution of $z$.** Because $z$ is inherited and parents choose where children are born, the cross-sectional distribution of $z$ among young agents differs across locations and is endogenous to migration. We track it by the joint mass $\Phi_l(z)$.

Setting $\rho_z = 0$, $\lambda = 0$, and $L=1$ collapses the model to the standard occupational-choice previous baseline, with $a_i = z\epsilon_i$ iid each generation.

## 4 Technologies

**Human capital.** A young agent of group $g$ in location $l$ who chooses occupation $i$ and is matched with a teacher of human capital $h_{T,l}$ acquires
$$
h \;=\; (h_{T,l})^{\beta}\,(z\,\epsilon_i)^{\alpha}\,s^{\phi}\,e^{\eta}\,N(h_{T,l})^{-\sigma},
$$
where $N(\cdot)$ is class size. Student-teacher matching is random within $l$; class sizes are allocated efficiently across the active teachers subject to the local slot constraint $\tfrac{M_l}{2} = \sum_g \int N(h)\,dF_{T,l,g}(h)$, where $F_{T,l,g}$ is the endogenous measure of human capital among teachers working in $l$. This yields:

**Proposition 1 (Class size and teacher quality).** *Efficient class-size allocation equalizes $h^{\beta}N(h)^{-\sigma}$ across all active teachers in $l$, with common value*
$$
h^{\beta}\,N(h)^{-\sigma} \;=\; \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^{\sigma} \;\equiv\; Q_l,
\qquad
\widetilde{H}_{T,l} \;\equiv\; \sum_g \int_0^\infty h^{\beta/\sigma}\,dF_{T,l,g}(h),
$$
*so better teachers take larger classes, $N(h) = h^{\beta/\sigma}\,Q_l^{-1/\sigma}$, and the student technology reduces to*
$$
h_{i,g,l}(z,\epsilon_i) \;=\; Q_l\,(z\,\epsilon_i)^{\alpha}\,s^{\phi}\,e^{\eta}.
$$

(Derivation in Appendix A.1.) Only the single moment $\widetilde{H}_{T,l}$ of the teacher distribution enters the rest of the model.

**Production and wages.** Non-teaching output is linear in human capital, $y_i = A_i h$, so a non-teacher earns $\omega_{i,l'}(h) = A_i h$ regardless of work location — training location matters only through $h$. Teaching wages take the power form
$$
\omega_{T,l'}(h) \;=\; \kappa_{l'}\,h^{\gamma},
$$
with location-specific intercepts $\kappa_{l'}$ and common elasticity $\gamma$. The nonlinearity ($\gamma \neq 1$) is what lets absolute advantage, and not just comparative advantage, shape selection into teaching.

## 5 Preferences and the Separation of Altruism

A young agent born in $l$, group $g$, with state $(z,\vec\epsilon)$, who chooses occupation $i$, values
$$
U \;=\; \ln(1-s) \;+\; \mathbb{E}\!\left[\, \mu \ln C \;+\; B_{l'} \;+\; \sigma_{\nu}\nu_{l'} \;+\; \lambda\,f(h')\,\right],
$$
with altruism strength $\lambda \geq 0$ and warm-glow kernel $f$ specified in §5.2. With $\varepsilon = 0$ there is no parental-finance term, and old-age consumption conditional on working in $l'$ is
$$
C_{l'} \;=\; (1-t_{l'})\,(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l}(z,\epsilon_i)\big) \;-\; (1+\tau^{e}_{i,g})\,e \;-\; m_{l,l'},
$$
which is **deterministic** given the agent's own state and choices: it does not depend on the child's realized $(z',\vec\epsilon',g')$. The expectation over the child therefore applies only to the warm-glow term, and the Stage-2 value of working in $l'$ separates additively:
$$
V_{l'} \;=\; \mu \ln C_{l'} \;+\; \Lambda_{l'}(z),
\qquad
\Lambda_{l'}(z) \;\equiv\; \lambda\,\tfrac{1}{2}\sum_{g'}\mathbb{E}_{z'\mid z}\,\mathbb{E}_{\vec\epsilon'}\!\left[f\big(h'(z',\vec\epsilon',l',g')\big)\right].
$$
The altruism term $\Lambda_{l'}(z)$ depends on the child's birth location $l'$ (through $Q_{l'}$ and the child's policies) and on the parent's $z$ (through the AR(1)), but **not** on the parent's investments $(s,e)$ or occupation. Two consequences follow, and both are preserved by everything in §5.1–§5.2. First, $\Lambda$ differentiates to zero in the investment first-order conditions: altruism affects investment only indirectly, through the location-choice probabilities. Second, $\Lambda$ is occupation-independent, so it will cancel from every pairwise occupational comparison at fixed location weights — but not from location choice itself, where it is a "move-to-opportunity" motive alongside the amenity $B_{l'}$.

Note that the moving cost has been moved *inside* $C_{l'}$ and correspondingly deleted from the utils bracket. Setting $m \equiv 0$ and restoring $-\tau_{l,l'}$ recovers the prior baseline.

### 5.1 The moving cost in goods

Two things change relative to a utils-denominated cost, one substantive and one technical.

*Substantive.* Under $\mu\ln C$ a proportional wage gap is valued identically at every income level, and every location premium in the baseline model is proportional — non-teaching wages $A_i h$ do not depend on $l'$ at all, teaching wages differ by the multiplicative $\kappa_{l'}$, and taxes are proportional. A *fixed* goods charge is the exception: it is a larger share of a poor agent's consumption than of a rich one's, so
$$
\frac{\partial \ln C_{l'}}{\partial \ln z}
\;=\; \frac{\alpha}{1-\eta}\left(1 + \frac{m_{l,l'}}{C_{l'}}\right)
\qquad (i \neq T),
$$
which is destination-*specific*. High-$z$ agents are less deterred by distance and therefore concentrate wherever the pull is (high $B$, high $\kappa$, high $Q$). This is the mechanism of Proposition 4(b) below and the canonical urban-sorting force (Diamond 2016).

*Technical.* The moving cost is not homogeneous in the agent's ability index, so the cancellation behind Proposition 2 no longer runs to a constant. What survives is Proposition 2′ in §7: $s$ remains **explicit given a single scalar**, the choice-weighted cost-to-consumption ratio
$$
\Xi \;\equiv\; \sum_{l'}\pi_{l'\mid l}\,\frac{m_{l,l'}}{C_{l'}} \;\geq\; 0 .
$$
Because $m_{l,l}=0$, $\Xi$ is bounded by the outflow rate times $m/\bar C$ and is small in any sensible calibration; $s$ deviates from its baseline constant only at $O(\Xi)$.

### 5.2 The warm-glow kernel

Replace $f(h') = \log h'$ with the CRRA-type kernel
$$
f(h') \;=\; \frac{(h'/\bar h)^{1-\psi} - 1}{1-\psi},
\qquad \psi \geq 0,
$$
where $\bar h > 0$ is a **fixed constant** — the benchmark mean child human capital, set once and held fixed across any $\psi$ experiment, so that $\psi$ changes only curvature and not the strength of the altruism motive. Since $f \to \log(h'/\bar h)$ as $\psi \to 1$ and $\lambda\log\bar h$ is a constant common to all $(l',z)$, the case $\psi = 1$ reproduces the baseline exactly: a constant shift in $\Lambda$ drops out of both location choice and occupational choice.

Why it sorts. Write the child's human capital as $h' = H_{l'}\cdot u$, where $H_{l'} \propto Q_{l'}^{1/(1-\eta)}$ collects the location component and $u$ the child's own ability and shocks. Under $\log$, $\Lambda_{2}-\Lambda_{1} = \tfrac{\lambda}{1-\eta}\log(Q_2/Q_1)$ is a **constant** — a better school is worth the same number of utils to every parent, which is precisely why the baseline generates no altruism-driven sorting. Under the power kernel the gap is **multiplicative** in the child's ability,
$$
\Lambda_{2}(z)-\Lambda_{1}(z)
\;=\; \frac{\lambda}{(1-\psi)\,\bar h^{\,1-\psi}}\Big(H_2^{1-\psi} - H_1^{1-\psi}\Big)\,
\mathbb{E}\big[u^{1-\psi}\,\big|\,z\big],
\qquad
\mathbb{E}\big[u^{1-\psi}\mid z\big] \;\propto\; z^{\,m_\psi\rho_z}\,e^{\,m_\psi^{2}\sigma_\xi^{2}/2},
$$
with $m_\psi \equiv \dfrac{(1-\psi)\alpha}{1-\eta}$ (derivation in Appendix A.4). A better school is worth *more* to a parent whose child is better; $\log$ is exactly the knife-edge where it is worth the same to everyone. The sign of the induced sorting is the sign of $m_\psi$: $\psi < 1$ gives positive assortative matching, $\psi > 1$ negative. Note that $\Lambda_2 > \Lambda_1$ whenever $H_2 > H_1$ for every $\psi$ — the better location is always preferred; only the *gradient* changes sign.

Because $\Lambda$ is still investment-independent and occupation-independent, **nothing in Propositions 1, 2, or 3 is affected by $\psi$.** The kernel touches only the definition of $\Lambda$ and the aggregation that produces it (§8).

## 6 Location Choice

Having chosen occupation $i$ when young, the agent observes the Gumbel shocks at Stage 2 and picks
$$
l^{*} \in \arg\max_{l'}\; V_{l'} + B_{l'} + \sigma_\nu \nu_{l'},
$$
yielding logit choice probabilities and a log-sum value:
$$
\pi_{l'\mid l} \;=\; \frac{\exp\!\big[(V_{l'} + B_{l'})/\sigma_\nu\big]}{\sum_{l''}\exp\!\big[(V_{l''} + B_{l''})/\sigma_\nu\big]},
\qquad
\bar V \;=\; \sigma_\nu \ln\!\sum_{l'}\exp\!\big[(V_{l'} + B_{l'})/\sigma_\nu\big],
$$
where the moving cost now enters through $V_{l'} = \mu\ln C_{l'} + \Lambda_{l'}(z)$ rather than the utils bracket. The log-sum $\bar V$ folds the Stage-2 problem into Stage-1 choices. As $\sigma_\nu \to 0$ location choice becomes deterministic with the moving cost acting as a wedge in net values; as $m_{l,l'} \to C_{l'}$ agents stay put.

**Whether migration probabilities depend on ability is a central question**, and in the baseline the answer is no.

**Proposition 4 (Ability-neutrality and how it breaks).** *Fix $(g,l,z)$ and the altruism term $\Lambda$. Take the non-teaching branch and write the ability index $K \equiv Q_l\,z^{\alpha}X_O^{*}(s_O^{*})^{\phi}$, so that income capacity is $Ke^{\eta}$.*

*(a) Baseline neutrality. If $m_{l,l'}\equiv 0$, then for any $\theta>0$*
$$
e^{*}(\theta K) = \theta^{\frac{1}{1-\eta}}e^{*}(K),
\qquad
C_{l'}(\theta K) = \theta^{\frac{1}{1-\eta}}C_{l'}(K)\ \ \forall l',
\qquad
\bar V(\theta K) = \bar V(K) + \tfrac{\mu}{1-\eta}\ln\theta,
$$
*and hence $\pi^{*}(\theta K) = \pi^{*}(K)$ exactly: location probabilities are **invariant** to the agent's ability draws. (Proof: substituting $e \mapsto \theta^{1/(1-\eta)}e$ scales both $Ke^{\eta}$ and $(1+\tau^e)e$ by $\theta^{1/(1-\eta)}$, so every $C_{l'}$ scales by the same factor and every $V_{l'}$ shifts by the same constant; a common shift leaves logit weights unchanged. The teaching branch is identical with $\eta \mapsto \gamma\eta$.) With $f = \log$, $\Lambda_{2}-\Lambda_{1}$ is also constant in $z$, so the model delivers zero within-branch sorting: all cross-location variation in mean ability comes from occupational composition, and since $\gamma<1$ makes teaching negatively selected on $z$, that residual channel has the empirically wrong sign.*

*(b) Both extensions break it, with gradients*
$$
g \;\equiv\; \frac{\partial\,[V_2 - V_1]}{\partial \ln z}
\;=\;
\underbrace{\mu\,\frac{\alpha}{1-\eta}\left(\frac{m_{l,2}}{C_2} - \frac{m_{l,1}}{C_1}\right)}_{\text{§5.1: goods cost}}
\;+\;
\underbrace{m_\psi\,\rho_z\,\big[\Lambda_2(z)-\Lambda_1(z)\big]}_{\text{§5.2: power warm glow}} ,
$$
*(non-teaching branch; replace $\alpha/(1-\eta)$ by $\gamma\alpha/(1-\gamma\eta)$ for teaching). Both terms vanish exactly at $m\equiv0$, $\psi=1$.*

The two channels are complementary in what they reach: the cost term applies to **everyone** and is largest for the agents most exposed to distance; the warm-glow term applies through the **teacher-quality** primitive $Q_{l'}$ that the paper is about, and is itself increasing in $z$, so sorting concentrates at the top of the ability distribution. Neither requires a new free parameter.

*Magnitude.* In the logit, the swing in the probability of working in location 2 across the ability support is approximately $\Delta\pi_2 \approx \bar\pi_1\bar\pi_2\,g\cdot\mathrm{range}(\ln z)/\sigma_\nu$. At the benchmark calibration ($\alpha=0.3$, $\eta=0.2$, $\rho_z=0.9$, $\sigma_\xi=0.2$, $\sigma_\nu=0.2$, $\mu=1$) the ability support spans $\mathrm{range}(\ln z)\approx1.84$ and $\bar\pi_1\bar\pi_2\approx0.25$, so a gradient of $g\approx0.06$–$0.10$ moves $\pi_2$ by 14–23 pp across the support — enough to dominate the composition channel, since teachers are only about a fifth of agents. Both extensions land there at parameter values that are *already in the model*: re-denominating $\tau_{l,l'} = 0.20$ utils gives $m/\bar C = 1-e^{-0.20/\mu} = 0.18$ and hence $g \approx 0.068$; and at $\psi = 0$ with the benchmark $\Lambda_2-\Lambda_1 \approx 0.19$, $m_\psi\rho_z = 0.34$ gives $g \approx 0.06$.

## 7 Human Capital Investment

When young, the agent solves, for each candidate occupation,
$$
W_{i} \;=\; \max_{s,e}\;\Big[\ln(1 - s) + \bar V\big(h_{i,g,l}(z,\epsilon_i);\, i,g,l,z\big)\Big],
$$
then picks the occupation with the highest $W_i$. Because $\bar V$ is a log-sum, its derivative with respect to any investment is the $\pi$-weighted average of the per-location derivatives; and because $\Lambda_{l'}$ is investment-independent, only the consumption terms matter. The two margins behave very differently.

**Proposition 2′ (Conditional closed-form time investment).** *With $\varepsilon = 0$, and with $\Xi \equiv \sum_{l'}\pi_{l'\mid l}\,m_{l,l'}/C_{l'}$ the choice-weighted moving-cost-to-consumption ratio, the optimal time share is*
$$
s_O^{*} \;=\; \frac{\mu\phi\,(1+\Xi)}{\mu\phi\,(1+\Xi) + 1 - \eta}
\qquad (i \neq T),
\qquad\qquad
s_T^{*} \;=\; \frac{\mu\phi\gamma\,(1+\Xi)}{\mu\phi\gamma\,(1+\Xi) + 1 - \gamma\eta},
$$
*independent of $(z,\vec\epsilon)$, the birth location, the destination wedges $(t_{l'},\kappa_{l'},B_{l'})$, and the altruism strength $\lambda$ **except** through the scalar $\Xi$. At $m\equiv0$ we have $\Xi=0$ and Proposition 2 is recovered:*
$$
s_O^{*} = \frac{\mu\phi}{\mu\phi + 1-\eta},
\qquad
s_T^{*} = \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1}.
$$

The proof (Appendix A.2) rests on the same cancellation as before: income is a common power of $(s,e)$ across destinations, so the $e$-FOC pins the $\pi$-weighted average of $1/C_{l'}$ down to $e$ and $\Xi$, and substituting into the $s$-FOC eliminates every location-specific object *except* $\Xi$. Teaching differs from non-teaching only because income runs through the extra exponent $\gamma$; $\psi$ never appears, because $\Lambda$ is investment-independent for any kernel.

The practical content of Proposition 2′ is that the time margin does **not** become a second numerical optimization. $\Xi$ is a scalar function of $(e,\pi)$, so $s \leftarrow s(\Xi(e,\pi))$ is one extra line inside the existing $(e,\pi)$ fixed point rather than a joint $(s,e)$ solve, and the feedback is stabilizing (a higher $s$ raises $h$, raises $C$, lowers $\Xi$, lowers $s$). At the benchmark, $\Xi \approx 0.07$ and $s_O$ moves from $0.333$ to $\approx 0.349$: $s$ becomes genuinely state-dependent, but only mildly, so the baseline constants are excellent warm starts.

**The goods margin.** No such collapse applies to $e$. The FOC for a non-teaching occupation is
$$
\sum_{l'}\pi_{l'\mid l}\;\frac{\eta\,y_{l'}/e \;-\; (1+\tau^{e}_{i,g})}{C_{l'}} \;=\; 0,
\qquad
y_{l'} \equiv (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}(h),
$$
(for teaching, $\eta \mapsto \gamma\eta$), a one-dimensional fixed point in $e$ coupled to the choice probabilities $\pi_{l'\mid l}$, which themselves depend on $e$ through $C_{l'}$. Equivalently — and this is how the solver treats it — $e$ directly maximizes the inclusive value $\bar V$, which is a smooth, single-peaked one-dimensional problem. The FOC collapses to a closed form only when consumption is common across active destinations ($L=1$, or $\sigma_\nu\to 0$, or symmetric net returns *and* $m\equiv0$), in which case $C = \tfrac{1-\eta}{\eta}(1+\tau^e_{i,g})\,e$ for non-teaching and the familiar constant investment rates of the baseline draft reappear.

## 8 Occupational Choice

With $s$ pinned down by Proposition 2′, occupational choice compares the values $W_i$ across the $I$ occupations, each solved at its own optimal $(s,e)$. Two observations reduce this to a one-dimensional threshold problem.

**Proposition 3 (Sufficient statistic for non-teaching choice).** *Suppose the educational barrier $\tau^e_{i,g}$ is common across non-teaching occupations. For $i \neq T$, define the income capacity for a given state, occupation pair*
$$
X_{O,i} \;\equiv\; \Theta_{i,g}\,\epsilon_i^{\alpha},
\qquad
\Theta_{i,g} \;\equiv\; (1-\tau^{\omega}_{i,g})\,A_i .
$$
*Then the non-teaching value $W_i$ depends on $(i,\epsilon_i)$ only through $X_{O,i}$ and is strictly increasing in it, so the best non-teaching option is $\arg\max_{i\neq T} X_{O,i}$, and*
$$
X_O^{*} \;\equiv\; \max_{i \neq T}\,\Theta_{i,g}\,\epsilon_i^{\alpha}
$$
*is a sufficient statistic for the entire non-teaching branch. Since $\log X_{O,i} \sim N\big(\log\Theta_{i,g} - \tfrac{\alpha\sigma_\epsilon^2}{2},\, \alpha^2\sigma_\epsilon^2\big)$ independently across $i$, the distribution of $X_O^*$ is the product of lognormal cdfs:*
$$
F_{X^*_O}(x) \;=\; \prod_{i \neq T} F_{X_{O,i}}(x).
$$

The logic: non-teaching after-tax income in $l'$ is $(1-t_{l'})\,\Theta_{i,g}\,h = (1-t_{l'})\,Q_l\,z^\alpha X_{O,i}\,(s_O^*)^\phi e^\eta$, so consumption, the $e$-FOC, and the value all see $(i,\epsilon_i)$ only through $X_{O,i}$; with a common $\tau^e$ nothing else varies with $i$. The agent's occupational problem therefore reduces to the pair of scalars $(\epsilon_T, X_O^*)$ — the ability draws in the $I-2$ dominated non-teaching occupations are irrelevant.

**Both extensions preserve Proposition 3.** The moving cost $m_{l,l'}$ depends on $(l,l')$ but not on $i$ or $\epsilon_i$, so it enters $C_{l'}$ without disturbing the collapse; and although $s_O^*$ is now state-dependent through $\Xi$, $\Xi$ itself is a function of $(C_{l'})_{l'}$, hence of $(i,\epsilon_i)$ only through $X_{O,i}$. So $s_O^*, e^*, C_{l'}, W_i$ all remain functions of $X_{O,i}$ alone, and $W_i$ remains strictly increasing in it by the envelope theorem. The warm-glow kernel $\psi$ does not appear in the occupational comparison at all, since $\Lambda$ is occupation-independent. Monotonicity of $W_T(\epsilon_T)$ and $W_O(X_O^*)$ — which the threshold inversion below relies on — is likewise unaffected.

**The teaching margin.** Write $W_T(\epsilon_T)$ and $W_O(X_O^*)$ for the two branch values at a fixed $(g,l,z)$ (each inclusive of its $\ln(1-s^*)$ term and its own optimal $e$). Both are strictly increasing in their scalar argument. The agent teaches iff $W_T(\epsilon_T) > W_O(X_O^*)$, so the indifference locus is a threshold: for each $\epsilon_T$ there is a cutoff $x^*(\epsilon_T) = W_O^{-1}\big(W_T(\epsilon_T)\big)$ such that the agent teaches iff $X_O^* < x^*(\epsilon_T)$. Conditional on $(g,l,z)$, the probability of teaching given the teaching draw is
$$
P\big(\text{teach} \mid \epsilon_T\big) \;=\; F_{X^*_O}\!\big(x^*(\epsilon_T)\big),
$$
and the unconditional teaching share integrates this against the density of $\epsilon_T$. This locus is the spatial, $\varepsilon=0$ analogue of the baseline draft's teacher/other cutoff $a^*_{T,g}(a_O)$: because $\gamma \neq 1$, the threshold depends on the *level* $z$ and not only on relative draws — absolute advantage matters for selection into teaching. Explicit closed forms for the threshold obtain in the destination-homogeneous case (Appendix A.3); in general it is computed by inverting the two monotone value functions.

**Smoothness.** Because $\epsilon_T$ and $X_O^*$ are continuously distributed, the indifference locus carries zero probability mass, and every aggregate below — the teaching share, $\widetilde H_T$, the tax base, the migration kernel, the altruism term $\Lambda$ — is a smooth functional of the value functions. No smoothing device (occupational taste shocks) is needed for the household fixed point to be well-behaved, in contrast to a discretized-$\vec\epsilon$ implementation where the hard argmax creates discontinuities.

**The household fixed point.** The warm-glow term closes a loop: values determine the child's occupational threshold and policies, which determine $\mathbb E[f(h')]$ and hence $\Lambda$, which enters the values. For the expectation, the non-teaching branch uses
$$
\log h_O \;=\; \log\!\big(Q_l\,z^{\alpha}(s_O^*)^{\phi} e^{\eta}\big) \;+\; \log X_O^* \;-\; \log \Theta_{i^*,g},
$$
so that under $\psi=1$ the relevant object is $\mathbb E\big[\log\Theta_{i^*,g}\mid X_O^*\big]$ — a conditional expectation over which occupation attained the max — while under a general $\psi$ it is
$$
h_O^{1-\psi} \;=\; \big(Q_l z^{\alpha}(s_O^*)^{\phi}e^{\eta}\big)^{1-\psi}\,(X_O^*)^{1-\psi}\;\mathbb E\big[\Theta_{i^*,g}^{-(1-\psi)}\,\big|\,X_O^*\big].
$$
The two coincide at $\psi=1$ and both are exact when there is a single non-teaching occupation ($I=2$), where $\Theta_{i^*,g}$ is degenerate; for $I>2$ the same $\mathrm{pdf}/\mathrm{CDF}$ argmax weights apply, and averaging $\Theta_i^{-(1-\psi)}$ rather than $\log\Theta_i$ is what removes a Jensen bias. The household block iterates $V \mapsto h' \mapsto \Lambda \mapsto V$ to convergence; by the smoothness property this map is a well-behaved contraction in practice, though with $\psi<1$ the term $\Lambda_{l'}(z) \propto z^{m_\psi\rho_z}$ grows as a power rather than a log of $z$, so location choice is more elastic at the top of the grid and the loop may need more sweeps.

## 9 Distributions and Aggregation

Let $\Phi_l(z)$ be the joint mass of young agents in location $l$ with persistent ability $z$, gender-pooled (per-gender masses are $\tfrac12\Phi_l$ by the gender coin flip), and let $\mathbb 1_{i,g,l}(z,\vec\epsilon)$ indicate the optimal occupation. All aggregates are integrals of policy quantities against $dF_\epsilon^{\otimes I}\,\Phi_l(dz)$; by §8 each such integral reduces to one-dimensional integrals over $\epsilon_T$ (teaching side) and $X_O^*$ (non-teaching side), split at the threshold.

**Migration kernel.** The probability that an agent born in $l_0$ with ability $z$ works in $l'$, averaging over gender and integrating over occupation choice and shocks:
$$
\bar\pi_{l_0 \to l'}(z) \;=\; \tfrac12 \sum_g \int \sum_i \mathbb 1_{i,g,l_0}(z,\vec\epsilon)\;\pi_{l'\mid l_0}\big(h_{i,g,l_0};\,i,g,z\big)\; dF_\epsilon^{\otimes I}(\vec\epsilon).
$$

**Law of motion for $\Phi_l$.** A parent in $(l_0,z)$ locates her child in $l$ with probability $\bar\pi_{l_0\to l}(z)$, and the child draws $z' \sim \pi_z(\cdot \mid z)$:
$$
\Phi_{l,t+1}(z') \;=\; \sum_{l_0}\int \pi_z(z'\mid z)\;\bar\pi_{l_0\to l}(z)\;\Phi_{l_0,t}(dz).
$$
In stationary equilibrium $\Phi$ is the dominant left eigenfunction (eigenvalue 1) of the kernel $\pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)$ on $\{1,\ldots,L\}\times\mathbb R_+$, scaled so that $\sum_l 2\int\Phi_l(dz) = M$. Total population per location is $M_l = 2\int\Phi_l(dz)$, and the gender composition of workers in each location follows from the gender-specific version of $\bar\pi$.

**Teacher aggregator.** Teachers trained anywhere who work in $l'$ contribute to $\widetilde H_{T,l'}$:
$$
\widetilde{H}_{T,l'} \;=\; \tfrac{1}{2}\sum_g\sum_{l_0}\!\int\!\!\int \mathbb 1_{T,g,l_0}(z,\vec\epsilon)\;\pi_{l'\mid l_0}\;h_{T,g,l_0}(z,\epsilon_T)^{\beta/\sigma}\; dF_\epsilon^{\otimes I}(\vec\epsilon)\,\Phi_{l_0}(dz).
$$
The full teacher distribution $F_{T,l',g}$ is the corresponding push-forward, but only this moment enters the model.

## 10 Government Budget

Each location runs a separate balanced budget: payroll taxes on labor income earned in $l'$ fund teacher wages paid in $l'$,
$$
t_{l'}\,\mathcal{I}_{l'} \;=\; \mathcal{W}_{l'},
$$
where $\mathcal I_{l'}$ is aggregate taxable income earned in $l'$ (net of the wedges $\tau^\omega_{i,g}$, summed over both teachers and non-teachers, all genders, and all training locations) and $\mathcal W_{l'}$ is the gross teaching wage bill $\int \omega_{T,l'}(h)$ over teachers working in $l'$. Both integrate against the same objects as §9: the occupation indicator, the migration probabilities, and $\Phi_{l_0}$.

## 11 Stationary Equilibrium

A stationary equilibrium is a tuple
$$
\Big\{\,t_l,\; \widetilde{H}_{T,l},\; M_l,\; \Phi_l,\; s_{i,g,l}(z,\epsilon_i),\; e_{i,g,l}(z,\epsilon_i),\; i^*(\,\cdot\,;g,l)\,\Big\}_{l,\,g}
$$
such that:

1. **(Optimization)** Given $(\widetilde H_{T,l}, M_l, t_l)_l$, the goods policies solve the $e$-FOC of §7 jointly with the location probabilities and with $s$ from Proposition 2′, occupational choice follows the threshold rule of §8, and the altruism term $\Lambda$ is consistent with the children's policies (the household fixed point).
2. **(Spatial distribution)** $\Phi_l$ is the stationary fixed point of the law of motion in §9, with $M_l = 2\int\Phi_l(dz)$.
3. **(Teacher aggregate)** $\widetilde H_{T,l}$ equals the aggregator in §9 evaluated at the equilibrium policies and $\Phi$.
4. **(Class size)** $h^\beta N(h)^{-\sigma} = Q_l$ for all teachers in $l$, with $\sum_g\int N\,dF_{T,l,g} = M_l/2$.
5. **(Budget)** $t_l\,\mathcal I_l = \mathcal W_l$ in each location.

Computationally, the outer fixed point runs over the small vector of aggregates $(\widetilde H_{T,l}, M_l, t_l)_l$: guess, solve the household block, integrate to update, and iterate with damping. The inner household block is the $V \mapsto \Lambda \mapsto V$ loop of §8, with the per-node $(s,e,\pi)$ fixed point of Proposition 2′ nested inside it.

**A caution on multiplicity.** Positive assortative matching introduces a complementarity between $\Phi_l$ and $Q_l$ — better students attract better teachers, which raises $Q_l$, which attracts better students — so the outer loop may require heavier damping and multiple equilibria become a live possibility. Sweeps over $m$ or $\psi$ should be run as continuations from the neighbouring solution so that a single branch is tracked.

**A caution on units.** Under $\mu\ln C$ with a utils-denominated moving cost the levels of $A_i$, $\kappa_l$, and $\bar h$ are innocuous: the model is invariant to a common rescaling of consumption. Both extensions break that invariance deliberately — $m_{l,l'}/C_{l'}$ and $h'/\bar h$ are the ratios that do the work — so a consumption normalization ($\bar C$ at the benchmark, with $\bar h$ the benchmark mean child human capital) must be fixed before any comparative static is run. Invariance to scaling the *population* $M$, under which $Q$, $h$, and $C$ are all unchanged, is unaffected and should continue to hold.

---

## Appendix A — Derivations

### A.1 Class size and the teacher-quality shifter $Q_l$

A teacher with human capital $h$ working in $l$ is paid $\omega_{T,l}(h)$ regardless of class size. Class sizes are allocated across active teachers to maximize aggregate student human capital subject to the slot constraint $\sum_j N_j = M_l/2$. Random matching makes each class's contribution proportional to $N_j^{1-\sigma}h_j^{\beta}$ (times a common average of student inputs), so the planner's FOC equalizes $h^{\beta}N(h)^{-\sigma}$ across teachers at some constant $\mathcal Q_l$. Then $N(h) = h^{\beta/\sigma}\mathcal Q_l^{-1/\sigma}$, and substituting into the slot constraint,
$$
\frac{M_l}{2} \;=\; \mathcal Q_l^{-1/\sigma}\sum_g\int h^{\beta/\sigma}dF_{T,l,g}(h) \;=\; \mathcal Q_l^{-1/\sigma}\,\widetilde H_{T,l}
\quad\Longrightarrow\quad
\mathcal Q_l = \left(\frac{2\widetilde H_{T,l}}{M_l}\right)^{\sigma} = Q_l .
$$
Substituting $h_{T}^{\beta}N^{-\sigma} = Q_l$ into the raw technology gives the reduced form of Proposition 1; dependence on the specific teacher drops out.

### A.2 Closed form for the time share $s$ (Proposition 2′)

Fix a state $(i,g,l,z,\vec\epsilon)$ with $i \neq T$ and write
$$
C_{l'} \;=\; y_{l'} - mc - m_{l,l'},
\qquad
mc \equiv (1+\tau^e)e,
\qquad
y_{l'} = (1-t_{l'})(1-\tau^\omega_{i,g})A_i h,
\qquad
h = Q_l(z\epsilon_i)^\alpha s^\phi e^\eta .
$$
Since $\Lambda$ is investment-independent — for **any** warm-glow kernel, so $\psi$ plays no role here — the FOCs weight per-destination consumption terms by $\pi_{l'\mid l}$ (the log-sum envelope). Define
$$
R \;\equiv\; \sum_{l'}\pi_{l'\mid l}\frac{1}{C_{l'}},
\qquad
\Xi \;\equiv\; \sum_{l'}\pi_{l'\mid l}\frac{m_{l,l'}}{C_{l'}} .
$$

*Step 1 (the $e$-FOC).* Since $h \propto e^{\eta}$, $\partial y_{l'}/\partial e = \eta y_{l'}/e$, and the $e$-FOC reads $\sum_{l'}\pi_{l'\mid l}\,[\eta y_{l'}/e - (1+\tau^e)]/C_{l'} = 0$. Substituting $y_{l'} = C_{l'} + mc + m_{l,l'}$ and using $\sum_{l'}\pi_{l'\mid l}=1$, then multiplying by $e$:
$$
\eta \;+\; \eta\big(mc\,R + \Xi\big) \;-\; mc\,R \;=\; 0
\qquad\Longleftrightarrow\qquad
mc\,R \;=\; \frac{\eta\,(1+\Xi)}{1-\eta}. \tag{$\star'$}
$$
At $m\equiv0$ this is the original $(\star)$: $R = \eta/[(1-\eta)(1+\tau^e)e]$.

*Step 2 (the $s$-FOC).* Since $h \propto s^{\phi}$, the $s$-FOC reads $\tfrac{1}{1-s} = \tfrac{\mu\phi}{s}\sum_{l'}\pi_{l'\mid l}\,y_{l'}/C_{l'}$. But $y_{l'}/C_{l'} = 1 + mc/C_{l'} + m_{l,l'}/C_{l'}$, so the weighted sum is $1 + mc\,R + \Xi$, and by $(\star')$
$$
\sum_{l'}\pi_{l'\mid l}\frac{y_{l'}}{C_{l'}}
\;=\; (1+\Xi)\left[1 + \frac{\eta}{1-\eta}\right]
\;=\; \frac{1+\Xi}{1-\eta}.
$$
Every location-specific object has cancelled *except* through the single scalar $\Xi$. Hence $s(1-\eta) = \mu\phi(1+\Xi)(1-s)$:
$$
s_O^{*} \;=\; \frac{\mu\phi(1+\Xi)}{\mu\phi(1+\Xi) + 1 - \eta}.
$$
For teaching, $\omega_{T,l'} = \kappa_{l'}h^{\gamma}$ implies $\partial y_{l'}/\partial e = \gamma\eta\,y_{l'}/e$ and $\partial y_{l'}/\partial s = \gamma\phi\,y_{l'}/s$; the same two steps give $mc\,R = \gamma\eta(1+\Xi)/(1-\gamma\eta)$, then $\sum\pi\,y/C = (1+\Xi)/(1-\gamma\eta)$ and
$$
s_T^{*} \;=\; \frac{\mu\phi\gamma(1+\Xi)}{\mu\phi\gamma(1+\Xi) + 1 - \gamma\eta},
$$
which at $\Xi=0$ equals $\mu\phi\gamma/[(\mu\phi-\eta)\gamma+1]$ as before.

The cancellation used only $\sum_{l'}\pi_{l'\mid l}=1$ and the homogeneity of *income* in $(s,e)$ — never the symmetry of destinations — so these expressions hold for any $L$, any wedges, any $\lambda$, and any $\psi$. What the goods cost destroys is the homogeneity of *consumption*: $m_{l,l'}$ does not scale with $(s,e)$, which is exactly why it survives into the answer as $\Xi$ and, by the same token, why it breaks the ability-neutrality of Proposition 4(a).

Two things $(\star')$ still does not do. It does not solve for $e$: with heterogeneous $C_{l'}$ the $\pi$-weighted sum of reciprocals does not invert to a single power of $e$, so the goods margin remains a one-dimensional numerical problem. And it does not make $s$ explicit in primitives, since $\Xi$ depends on $(e,\pi)$ — but the dependence is through one scalar, so it is resolved by a short fixed point rather than a second optimization.

### A.3 Occupational thresholds under destination homogeneity

This subsection characterizes the baseline case $m\equiv0$; the two extensions leave it intact in the following sense. The kernel $\psi$ never enters, because $\Lambda$ is occupation-independent and so is absorbed into $D(z,l)$ below. The goods cost enters only through $\Xi$, which under destination homogeneity is the *common* constant $m/C$ shared by both branches; it shifts $s_O^*$ and $s_T^*$ but leaves the structure of both cutoffs — in particular which variables cancel — exactly as stated.

When active destinations deliver a common consumption level $C_{l'} \equiv C$ ($L=1$, or $\sigma_\nu \to 0$, or symmetric net returns), the log-sum factors as $\bar V = \mu\ln C(h_i) + D(z,l)$ with $D$ occupation-independent, and $(\star')$ solves in closed form: $C^{(i)} = (1-\eta)(1-t)(1-\tau^\omega_{i,g})A_i h_i$ for non-teaching and $C^{(T)} = (1-\gamma\eta)(1-t)(1-\tau^\omega_{T,g})\kappa_l h_T^{\gamma}$ for teaching, with
$$
h_i \propto (z\epsilon_i)^{\frac{\alpha}{1-\eta}} \quad (i\neq T),
\qquad
h_T \propto (z\epsilon_T)^{\frac{\alpha}{1-\gamma\eta}} .
$$
Two thresholds follow from $W_i = W_j$ (values affine in log ability):

*Non-teaching vs. non-teaching — comparative advantage only.* The common $s_O^*$ and exponents cancel $z$, $Q_l$, and every spatial wedge, leaving a cutoff in relative draws:
$$
\ln\frac{\epsilon_i}{\epsilon_j}
\;=\;
\frac{1}{\alpha}\ln\frac{(1-\tau^\omega_{j,g})A_j}{(1-\tau^\omega_{i,g})A_i}
\;-\;
\frac{\eta}{\alpha}\ln\frac{1+\tau^e_{j,g}}{1+\tau^e_{i,g}} .
$$
With a common $\tau^e$ this is exactly the statement that the agent ranks non-teaching occupations by $X_{O,i} = (1-\tau^\omega_{i,g})A_i\,\epsilon_i^{\alpha}$ — Proposition 3.

*Teaching vs. non-teaching — absolute advantage.* The ability exponents differ ($\tfrac{\gamma\alpha}{1-\gamma\eta} \neq \tfrac{\alpha}{1-\eta}$ unless $\gamma=1$), so $z$ does not cancel and the cutoff solves
$$
\frac{\gamma\alpha}{1-\gamma\eta}\,\ln(z\,\epsilon_T^{*})
\;=\;
\frac{\alpha}{1-\eta}\,\ln(z\,\epsilon_O)
\;+\;
\frac{1}{\mu}\ln\frac{1-s^*_O}{1-s^*_T}
\;+\;
\vartheta_O-\vartheta_T,
$$
where $\vartheta_O,\vartheta_T$ collect the (log) productivity, wedge, $Q_l$, and parameter constants of the two branches. The cutoff rises with the persistent level $z$ — the absolute-advantage channel — and depends on $\widetilde H_{T,l}$ through $Q_l$ except in the linear-wage knife-edge $\gamma = 1$, where the $Q_l$ terms cancel (the spatial analogue of the baseline's $\beta=\sigma$ case). Away from destination homogeneity these are leading-order characterizations; the exact threshold adds an occupation-specific destination-mixing correction and is found by inverting the monotone branch values as in §8.


### A.4 The altruism term under the power kernel

Take the non-teaching branch in the destination-homogeneous benchmark, where the goods policy is $e \propto \big[Q_{l'}(z'\epsilon)^{\alpha}s^{\phi}\big]^{1/(1-\eta)}$ and hence
$$
h' \;=\; H_{l'}\cdot u,
\qquad
H_{l'} \;\propto\; Q_{l'}^{\frac{1}{1-\eta}},
\qquad
u \;\propto\; (z')^{\frac{\alpha}{1-\eta}}\cdot(\text{shock terms}),
$$
with the teaching analogue replacing $1-\eta$ by $1-\gamma\eta$ and adding the $\gamma$ exponent. Applying the kernel,
$$
\Lambda_{l'}(z)
\;=\; \frac{\lambda}{1-\psi}\Big[\bar h^{-(1-\psi)}H_{l'}^{1-\psi}\,\mathbb E\big[u^{1-\psi}\mid z\big] - 1\Big].
$$
Under the log-AR(1) $\log z' = \rho_z\log z + \sigma_\xi\xi$, for any exponent $m$,
$$
\mathbb E\big[(z')^{m}\,\big|\,z\big] \;=\; z^{m\rho_z}\,e^{m^{2}\sigma_\xi^{2}/2}
\qquad\Longrightarrow\qquad
\frac{\partial \ln \mathbb E\big[(z')^{m}\mid z\big]}{\partial \ln z} \;=\; m\rho_z ,
$$
so with $m_\psi = (1-\psi)\alpha/(1-\eta)$ the *location gap* inherits a constant log-elasticity in $z$:
$$
\Lambda_2(z)-\Lambda_1(z) \;\propto\; z^{\,m_\psi\rho_z}
\qquad\Longrightarrow\qquad
\frac{\partial\big[\Lambda_2 - \Lambda_1\big]}{\partial \ln z} \;=\; m_\psi\,\rho_z\,\big[\Lambda_2(z)-\Lambda_1(z)\big],
$$
which is the second term of Proposition 4(b). Three observations. First, $m_\psi \to 0$ as $\psi\to1$, recovering the constant gap of the baseline: the log kernel is the unique knife-edge at which school quality is worth the same in utils to every parent. Second, the elasticity is constant but the *level* $\Lambda_2-\Lambda_1$ is increasing in $z$ when $\psi<1$, so the gradient itself steepens with ability and sorting concentrates at the top. Third, the multiplicative structure is what makes the reference $\bar h$ matter: $\log$ is scale-free in differences but $h'^{1-\psi}$ is not, so with $h' \ll 1$ in the benchmark calibration a literal $\bar h = 1$ would shrink $\Lambda_2-\Lambda_1$ by nearly an order of magnitude and the mechanism would appear to do nothing. Setting $\bar h$ to benchmark mean child human capital makes $h'/\bar h \approx 1$, so $f(h') \approx \log(h'/\bar h)$ to first order for every $\psi$ and the $\psi$-sweep isolates curvature.
