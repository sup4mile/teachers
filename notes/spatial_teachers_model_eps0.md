# Spatial Model with Altruism ($\varepsilon = 0$)

This document presents the spatial occupational-choice / teacher-spillover economy under the assumption that the parental cost share is zero. Parents still choose the location in which their children are born and educated, and still care about their children's human capital through a warm-glow term, but they finance none of the children's educational expenditures.

Imposing $\varepsilon = 0$ from the beginning buys three structural simplifications, each developed below:

1. **Consumption decouples from the child**. Old-age consumption is deterministic given the agent's own state, so the location value separates additively into a consumption term and a warm-glow term. Altruism shapes behavior only through location choice.
2. **The time investment is closed form**. The optimal time share $s$ is a constant per occupation — independent of ability, birth location, destination wedges, and the altruism strength.
3. **Non-teaching occupations collapse to a scalar**. The best non-teaching option is summarized by a one-dimensional sufficient statistic with known distribution, and the teach/not-teach margin is a single-crossing threshold. The $I$-dimensional choice problem becomes two-dimensional.

## 1 Environment

There are $L=2$ locations indexed by $l \in \{1,2\}$. In each location the same set of $I$ occupations is available, with $i=T$ denoting teaching. Locations share common occupational productivities $\{A_i\}_{i\neq T}$, labor-market distortions $\{\tau^{\omega}_{i,g}\}$, and educational barriers $\{\tau^{e}_{i,g}\}$. Locations differ along:

- the (endogenous) stock of teacher human capital $\widetilde{H}_{T,l}$,
- the (endogenous) population or mass of agents in each location $M_l$
- the local labor-income tax $t_l$, which funds teachers under a local balanced-budget rule,
- an exogenous amenity value $B_l$,
- the teaching wage intercept $\kappa_l$.

A young agent born in $l$ who later works in $l'$ pays a utility moving cost $\tau_{l,l'} \geq 0$, with $\tau_{l,l}=0$.

The economy is populated by overlapping generations of finite-lived agents (young and old) with constant total measure $M$. Population in each location $l$ has equal cohort sizes $M_l/2$ within each location in stationary equilibrium; the split $\{M_l\}$ is pinned down by location flows. Groups $g \in \{m,f\}$ denote gender. A parent's child is male or female with equal probability, independently of the parent's gender, so young agents are evenly split by gender in every location. Because distortions $\tau^\omega_{i,g}$, barriers $\tau^e_{i,g}$, and hence policies are gender-specific, the coin flip over the child's gender is payoff-relevant to the parent, and the gender mix among *old* (working-age) agents in each location is endogenous — men and women sort differently across locations.

## 2 Timing

Each cohort's life cycle proceeds in three stages.

**Stage 1 — Young.** Born in location $l$ (chosen by the parent), with persistent ability $z$ inherited from the parent. The vector of occupation-specific idiosyncratic shocks $\vec\epsilon = (\epsilon_1,\ldots,\epsilon_I)$ is drawn iid from $F_\epsilon$, giving ability $a_i = z\,\epsilon_i$ in each occupation. The agent chooses time investment $s$, goods investment $e$, and an occupation $i$. The full education cost $(1+\tau^{e}_{i,g})\,e$ is borne by the agent herself, paid out of future labor income.

**Stage 2 — Transition.** The agent reaches working age, draws iid Gumbel taste shocks $\{\nu_{l'}\}_{l'=1}^{L}$, and chooses a work location $l'$ in which to earn income, pay local taxes, and raise her own child.

**Stage 3 — Old.** She earns occupation-specific labor income in $l'$, pays the local tax and her education debt, and consumes the residual. Her child's persistent ability $z' \sim \pi_z(\cdot\mid z)$, idiosyncratic vector $\vec\epsilon'$, and gender $g'$ realize, and she receives a warm-glow payoff $\lambda f(h')$ from the child's accumulated human capital. There is no financial transfer between generations.

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
U \;=\; \ln(1-s) \;+\; \mathbb{E}\!\left[\, \mu \ln C \;+\; B_{l'} - \tau_{l,l'} \;+\; \sigma_{\nu}\nu_{l'} \;+\; \lambda\,f(h')\,\right],
$$
with warm-glow kernel $f(h') = \log h'$ and altruism strength $\lambda \geq 0$. With $\varepsilon = 0$ there is no parental-finance term, and old-age consumption conditional on working in $l'$ is
$$
C_{l'} \;=\; (1-t_{l'})\,(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l}(z,\epsilon_i)\big) \;-\; (1+\tau^{e}_{i,g})\,e,
$$
which is **deterministic** given the agent's own state and choices: it does not depend on the child's realized $(z',\vec\epsilon',g')$. The expectation over the child therefore applies only to the warm-glow term, and the Stage-2 value of working in $l'$ separates additively:
$$
V_{l'} \;=\; \mu \ln C_{l'} \;+\; \Lambda_{l'}(z),
\qquad
\Lambda_{l'}(z) \;\equiv\; \lambda\,\tfrac{1}{2}\sum_{g'}\mathbb{E}_{z'\mid z}\,\mathbb{E}_{\vec\epsilon'}\!\left[\log h'(z',\vec\epsilon',l',g')\right].
$$
The altruism term $\Lambda_{l'}(z)$ depends on the child's birth location $l'$ (through $Q_{l'}$ and the child's policies) and on the parent's $z$ (through the AR(1)), but **not** on the parent's investments $(s,e)$ or occupation. Two consequences follow. First, $\Lambda$ differentiates to zero in the investment first-order conditions: altruism affects investment only indirectly, through the location-choice probabilities. Second, $\Lambda$ is occupation-independent, so it will cancel from every pairwise occupational comparison at fixed location weights — but not from location choice itself, where it is a "move-to-opportunity" motive alongside the amenity $B_{l'}$.

## 6 Location Choice

Having chosen occupation $i$ when young, the agent observes the Gumbel shocks at Stage 2 and picks
$$
l^{*} \in \arg\max_{l'}\; V_{l'} + B_{l'} - \tau_{l,l'} + \sigma_\nu \nu_{l'},
$$
yielding logit choice probabilities and a log-sum value:
$$
\pi_{l'\mid l} \;=\; \frac{\exp\!\big[(V_{l'} + B_{l'} - \tau_{l,l'})/\sigma_\nu\big]}{\sum_{l''}\exp\!\big[(V_{l''} + B_{l''} - \tau_{l,l''})/\sigma_\nu\big]},
\qquad
\bar V \;=\; \sigma_\nu \ln\!\sum_{l'}\exp\!\big[(V_{l'} + B_{l'} - \tau_{l,l'})/\sigma_\nu\big].
$$
The log-sum $\bar V$ folds the Stage-2 problem into Stage-1 choices. As $\sigma_\nu \to 0$ location choice becomes deterministic with $\tau_{l,l'}$ acting as a wedge in net values; as $\tau_{l,l'} \to \infty$ agents stay put. Note that migration probabilities depend on the agent's ability: through $V_{l'}$, higher-$h$ agents weigh the net-of-tax wage differences differently, and through $\Lambda_{l'}(z)$, higher-$z$ parents weigh teacher-quality differences for their children.

## 7 Human Capital Investment

When young, the agent solves, for each candidate occupation,
$$
W_{i} \;=\; \max_{s,e}\;\Big[\ln(1 - s) + \bar V\big(h_{i,g,l}(z,\epsilon_i);\, i,g,l,z\big)\Big],
$$
then picks the occupation with the highest $W_i$. Because $\bar V$ is a log-sum, its derivative with respect to any investment is the $\pi$-weighted average of the per-location derivatives; and because $\Lambda_{l'}$ is investment-independent, only the consumption terms matter. The two margins behave very differently.

**Proposition 2 (Closed-form time investment).** *With $\varepsilon = 0$ the optimal time share is a constant per occupation,*
$$
s_O^{*} \;=\; \frac{\mu\phi}{\mu\phi + 1 - \eta}
\qquad (i \neq T),
\qquad\qquad
s_T^{*} \;=\; \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1},
$$
*independent of $(z,\vec\epsilon)$, the birth location, the destination wedges $(t_{l'},\kappa_{l'},B_{l'},\tau_{l,l'})$, and the altruism strength $\lambda$.*

The proof (Appendix A.2) rests on a cancellation: income is a common power of $(s,e)$ across destinations, so the $e$-FOC pins down the $\pi$-weighted average of $1/C_{l'}$ in terms of $e$ alone, and substituting into the $s$-FOC eliminates every location-specific object. Teaching differs from non-teaching only because income runs through the extra exponent $\gamma$.

**The goods margin.** No such collapse applies to $e$. The FOC for a non-teaching occupation is
$$
\sum_{l'}\pi_{l'\mid l}\;\frac{\eta\,g_{l'}/e \;-\; (1+\tau^{e}_{i,g})}{C_{l'}} \;=\; 0,
\qquad
g_{l'} \equiv (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}(h),
$$
(for teaching, $\eta \mapsto \gamma\eta$), a one-dimensional fixed point in $e$ coupled to the choice probabilities $\pi_{l'\mid l}$, which themselves depend on $e$ through $C_{l'}$. Equivalently — and this is how the solver treats it — $e$ directly maximizes the inclusive value $\bar V$, which is a smooth, single-peaked one-dimensional problem. The FOC collapses to a closed form only when consumption is common across active destinations ($L=1$, or $\sigma_\nu\to 0$, or symmetric net returns), in which case $C = \tfrac{1-\eta}{\eta}(1+\tau^e_{i,g})\,e$ for non-teaching and the familiar constant investment rates of the baseline draft reappear.

## 8 Occupational Choice

With $s$ pinned down by Proposition 2, occupational choice compares the values $W_i$ across the $I$ occupations, each solved at its own optimal $e$. Two observations reduce this to a one-dimensional threshold problem.

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

**The teaching margin.** Write $W_T(\epsilon_T)$ and $W_O(X_O^*)$ for the two branch values at a fixed $(g,l,z)$ (each inclusive of its $\ln(1-s^*)$ term and its own optimal $e$). Both are strictly increasing in their scalar argument. The agent teaches iff $W_T(\epsilon_T) > W_O(X_O^*)$, so the indifference locus is a threshold: for each $\epsilon_T$ there is a cutoff $x^*(\epsilon_T) = W_O^{-1}\big(W_T(\epsilon_T)\big)$ such that the agent teaches iff $X_O^* < x^*(\epsilon_T)$. Conditional on $(g,l,z)$, the probability of teaching given the teaching draw is
$$
P\big(\text{teach} \mid \epsilon_T\big) \;=\; F_{X^*_O}\!\big(x^*(\epsilon_T)\big),
$$
and the unconditional teaching share integrates this against the density of $\epsilon_T$. This locus is the spatial, $\varepsilon=0$ analogue of the baseline draft's teacher/other cutoff $a^*_{T,g}(a_O)$: because $\gamma \neq 1$, the threshold depends on the *level* $z$ and not only on relative draws — absolute advantage matters for selection into teaching. Explicit closed forms for the threshold obtain in the destination-homogeneous case (Appendix A.3); in general it is computed by inverting the two monotone value functions.

**Smoothness.** Because $\epsilon_T$ and $X_O^*$ are continuously distributed, the indifference locus carries zero probability mass, and every aggregate below — the teaching share, $\widetilde H_T$, the tax base, the migration kernel, the altruism term $\Lambda$ — is a smooth functional of the value functions. No smoothing device (occupational taste shocks) is needed for the household fixed point to be well-behaved, in contrast to a discretized-$\vec\epsilon$ implementation where the hard argmax creates discontinuities.

**The household fixed point.** The warm-glow term closes a loop: values determine the child's occupational threshold and policies, which determine $\mathbb E[\log h']$ and hence $\Lambda$, which enters the values. For the expectation, the non-teaching branch uses
$$
\log h_O \;=\; \log\!\big(Q_l\,z^{\alpha}(s_O^*)^{\phi} e^{\eta}\big) \;+\; \log X_O^* \;-\; \mathbb E\big[\log \Theta_{i^*,g} \,\big|\, X_O^*\big],
$$
where the conditional expectation accounts for which occupation attained the max. The household block iterates $V \mapsto h' \mapsto \Lambda \mapsto V$ to convergence; by the smoothness property this map is a well-behaved contraction in practice.

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
\Big\{\,t_l,\; \widetilde{H}_{T,l},\; M_l,\; \Phi_l,\; e_{i,g,l}(z,\epsilon_i),\; i^*(\,\cdot\,;g,l)\,\Big\}_{l,\,g}
$$
(with $s$ given by Proposition 2) such that:

1. **(Optimization)** Given $(\widetilde H_{T,l}, M_l, t_l)_l$, the goods policies solve the $e$-FOC of §7 jointly with the location probabilities, occupational choice follows the threshold rule of §8, and the altruism term $\Lambda$ is consistent with the children's policies (the household fixed point).
2. **(Spatial distribution)** $\Phi_l$ is the stationary fixed point of the law of motion in §9, with $M_l = 2\int\Phi_l(dz)$.
3. **(Teacher aggregate)** $\widetilde H_{T,l}$ equals the aggregator in §9 evaluated at the equilibrium policies and $\Phi$.
4. **(Class size)** $h^\beta N(h)^{-\sigma} = Q_l$ for all teachers in $l$, with $\sum_g\int N\,dF_{T,l,g} = M_l/2$.
5. **(Budget)** $t_l\,\mathcal I_l = \mathcal W_l$ in each location.

Computationally, the outer fixed point runs over the small vector of aggregates $(\widetilde H_{T,l}, M_l, t_l)_l$: guess, solve the household block, integrate to update, and iterate with damping. The inner household block is the $V \mapsto \Lambda \mapsto V$ loop of §8.

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

### A.2 Closed form for the time share $s$

Fix a state $(i,g,l,z,\vec\epsilon)$ with $i \neq T$ and write $C_{l'} = g_{l'} - (1+\tau^e)e$ with after-tax income $g_{l'} = (1-t_{l'})(1-\tau^\omega_{i,g})A_i h$ and $h = Q_l(z\epsilon_i)^\alpha s^\phi e^\eta$. Since $\Lambda$ is investment-independent, the FOCs weight per-destination consumption terms by $\pi_{l'\mid l}$ (the log-sum envelope).

*Step 1.* Since $h \propto e^{\eta}$, $\partial g_{l'}/\partial e = \eta g_{l'}/e$, and the $e$-FOC reads $\sum_{l'}\pi_{l'\mid l}\,[\eta g_{l'}/e - (1+\tau^e)]/C_{l'} = 0$. Substituting $g_{l'} = C_{l'} + (1+\tau^e)e$ and using $\sum_{l'}\pi_{l'\mid l}=1$:
$$
\sum_{l'}\pi_{l'\mid l}\,\frac{1}{C_{l'}} \;=\; \frac{\eta}{(1-\eta)(1+\tau^e)\,e}. \tag{$\star$}
$$

*Step 2.* Since $h \propto s^{\phi}$, the $s$-FOC reads $\tfrac{1}{1-s} = \tfrac{\mu\phi}{s}\sum_{l'}\pi_{l'\mid l}\,g_{l'}/C_{l'}$. But $g_{l'}/C_{l'} = 1 + (1+\tau^e)e/C_{l'}$, so by $(\star)$ the weighted sum equals $1 + \tfrac{\eta}{1-\eta} = \tfrac{1}{1-\eta}$ — every location-specific object has cancelled. Hence $s(1-\eta) = \mu\phi(1-s)$:
$$
s_O^{*} \;=\; \frac{\mu\phi}{\mu\phi + 1 - \eta}.
$$
For teaching, $\omega_{T,l'} = \kappa_{l'}h^{\gamma}$ implies $\partial g_{l'}/\partial e = \gamma\eta\,g_{l'}/e$ and $\partial g_{l'}/\partial s = \gamma\phi\,g_{l'}/s$; the same two steps give $\sum\pi\,g/C = 1/(1-\gamma\eta)$ and
$$
s_T^{*} \;=\; \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1}.
$$
The cancellation used only $\sum_{l'}\pi_{l'\mid l}=1$ and the homogeneity of income in $(s,e)$ — never the symmetry of destinations — so the constants hold for any $L$, any wedges, and any $\lambda$. It does not, however, solve $(\star)$ for $e$: with heterogeneous $C_{l'}$ the $\pi$-weighted sum of reciprocals does not invert to a single power of $e$, which is why the goods margin remains a one-dimensional numerical problem.

### A.3 Occupational thresholds under destination homogeneity

When active destinations deliver a common consumption level $C_{l'} \equiv C$ ($L=1$, or $\sigma_\nu \to 0$, or symmetric net returns), the log-sum factors as $\bar V = \mu\ln C(h_i) + D(z,l)$ with $D$ occupation-independent, and $(\star)$ solves in closed form: $C^{(i)} = (1-\eta)(1-t)(1-\tau^\omega_{i,g})A_i h_i$ for non-teaching and $C^{(T)} = (1-\gamma\eta)(1-t)(1-\tau^\omega_{T,g})\kappa_l h_T^{\gamma}$ for teaching, with
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
