# Spatial Model with Altruism

This document extends the baseline occupational-choice / teacher-spillover economy to a setting with $L=2$ locations and a "move-to-opportunity" altruistic link between generations. Parents choose the location in which their children are born and educated, and finance an exogenous share $\varepsilon$ of their children's educational expenditures. Locations differ in teacher quality, fiscal policy, and amenities; agents face Gumbel taste shocks over locations and an additive utility moving cost.

## 1 Environment

There are $L=2$ locations indexed by $l \in \{1,2\}$. In each location the same set of $I$ occupations is available, with $i=1$ (or $i=T$) denoting teaching. Locations share common occupational productivities $\{A_i\}_{i=2}^{I}$, common labor-market distortions $\{\tau^{\omega}_{i,g}\}$, and common educational barriers $\{\tau^{e}_{i,g}\}$. Locations differ along:

- the (endogenous) stock of teacher human capital $\widetilde{H}_{T,l}$,
- the local labor-income tax $t_l$, which funds teachers under a local balanced-budget rule (one per location),
- an exogenous amenity value $B_l$,
- (possibly) a location-specific teaching wage profile $\omega_{T,l}(\cdot)$.

In addition, a young agent born in $l$ who later works in $l'$ pays a utility moving cost $\tau_{l,l'} \geq 0$, with $\tau_{l,l}=0$. Moving costs attenuate the dispersion induced by the Gumbel taste shocks.

The economy is populated by overlapping generations of finite-lived agents (young and old) with constant total measure $M$. Mass in location $l$ is $M_l$, with $M_1+M_2=M$ and equal cohort sizes $M_l/2$ (young) and $M_l/2$ (old) within each location in stationary equilibrium. The split $\{M_l\}$ is pinned down by location flows. Groups $g \in \{m, f\}$ denote gender. A parent's child is male or female with equal probability, independently of the parent's gender. Labor-market distortions $\tau^\omega_{i,g}$, educational barriers $\tau^e_{i,g}$, and investment policies $(s_{i,g,l}, e_{i,g,l})$ are gender-specific, so the coin-flip over the child's gender $g'$ is payoff-relevant from the parent's perspective. Because gender is an independent coin flip, the gender share of young agents is $\tfrac{1}{2}$ in every location. The gender composition of old (working-age) agents in each location is endogenous, however: men and women face different distortions and sort differently across locations, so the gender mix among workers in $l'$ is determined in equilibrium by the flow equations.

## 2 Timing

Each cohort's life cycle proceeds in three stages.

**Stage 1 — Young.** Born in location $l$ (chosen by the parent in Stage 2 of the previous cohort), with persistent ability $z$ inherited from the parent. Conditional on $l$ and $z$:

1. The vector of occupation-specific idiosyncratic shocks $\vec\epsilon = (\epsilon_1,\ldots,\epsilon_I)$ is drawn iid from $F_\epsilon$, giving ability $a_i = z\,\epsilon_i$ in each occupation.
2. The agent chooses time investment $s_{i,g,l}$ and goods investment $e_{i,g,l}$ in occupation-specific human capital, and selects an occupation $i$.
3. The parent (currently old) finances a share $\varepsilon$ of $(1+\tau^{e}_{i,g})\,e_{i,g,l}$. The remaining share $1-\varepsilon$ is paid by the agent out of her own future labor income.

**Stage 2 — Transition.** The agent reaches working age. She draws i.i.d. Gumbel taste shocks $\{\nu_{l'}\}_{l'=1}^{L}$ and chooses a location $l' \in \{1,2\}$ in which to (i) earn labor income, (ii) pay local taxes, and (iii) raise and educate her own child. She forms expectations over the child's gender $g' \in \{m,f\}$ (each equally likely), the child's persistent ability $z' \sim \pi_z(\cdot\mid z)$ (the AR(1) transition), and the child's iid idiosyncratic vector $\vec\epsilon' \sim F_\epsilon^{\otimes I}$, and the resulting human capital $h'$.

**Stage 3 — Old.** Conditional on $l'$:

1. She earns occupation-specific labor income.
2. Her child's persistent state $z'$ and idiosyncratic vector $\vec\epsilon'$ are realized.
3. She pays the parental share $\varepsilon$ of the child's (now-realized) educational expenses.
4. She consumes the residual and receives a warm-glow payoff $\lambda\,f(h')$ from the child's accumulated human capital.

## 3 Ability Transmission

Ability is summarized by a persistent scalar $z \geq 0$ and a vector of occupation-specific idiosyncratic shocks $\vec\epsilon = (\epsilon_1,\ldots,\epsilon_I)$. The ability the agent draws on in occupation $i$ is the product
$$
a_i \;=\; z \cdot \epsilon_i.
$$

**Persistent component.** The scalar $z$ is the only ability state inherited across generations. It follows a log-AR(1):
$$
\log z' \;=\; \rho_z \log z \;+\; \sigma_{\xi}\,\xi, \qquad \xi \sim N(0,1)\ \text{iid},
$$
with $|\rho_z|<1$. Let $\pi_z(z'\mid z)$ denote the implied transition density and $G_z^{*}$ the (ergodic) unconditional distribution.

**Idiosyncratic component.** The shocks $\{\epsilon_i\}_{i=1}^I$ are iid across occupations and generations:
$$
\log \epsilon_i \;\sim\; N\!\left(-\tfrac{1}{2}\sigma_\epsilon^2,\; \sigma_\epsilon^2\right) \implies \mathbb E[\epsilon_i] = 1.
$$
Denote the common scalar law by $F_\epsilon$ and the joint distribution over $\vec\epsilon$ by $F_\epsilon^{\otimes I}$. The vector $\vec\epsilon$ is not inherited; the child draws a fresh $\vec\epsilon'$ independent of the parent's $\vec\epsilon$ and of $z'$.

**Information.** When the young agent solves her Stage-1 problem, both $z$ and $\vec\epsilon$ are observed. When the parent (in Stage 2 of the previous cohort) forms expectations over her child's future outcomes, the relevant expectation factors as $\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\,\mathbb E_{g'}$, with $z'\sim\pi_z(\cdot\mid z)$, $\vec\epsilon'\sim F_\epsilon^{\otimes I}$, and $g'\sim\mathrm{Unif}\{m,f\}$ all mutually independent.

**Endogenous spatial distribution of $z$.** Because $z$ is the only inherited component and parents choose where their children are born, the unconditional (cross-sectional) distribution of $z$ among young agents differs across locations and is endogenous to migration flows. We track this distribution by the joint mass $\Phi_l(z)$.

**Recovery of baseline.** Setting $\rho_z = 0$ makes $z$ iid across generations and the model with $L=1$, $\lambda = 0$, $\varepsilon = 0$ collapses to the standard occupational-choice baseline (with $a_i = z\,\epsilon_i$ iid each generation).

## 4 Technologies

**Human capital.** A young agent of group $g$ in location $l$ with persistent state $z$ and idiosyncratic vector $\vec\epsilon$ acquires occupation-specific human capital
$$
h_{i,g,l}(z,\epsilon_i) \;=\; (h_{T,l})^{\beta}\,(z\,\epsilon_i)^{\alpha}\,(s_{i,g,l})^{\phi}\,(e_{i,g,l})^{\eta}\,N(h_{T,l})^{-\sigma},
$$
with $h_{T,l}$ the human capital of the (randomly matched) teacher working in $l$. Student-teacher matching is random within $l$. The teacher's optimal class-size choice yields the indifference condition $h^{\beta} N(h)^{-\sigma}=\text{const}$ across all active teachers in $l$ (derived in Appendix A.1). Combined with the local class-size resource constraint
$$
\frac{M_l}{2} \;=\; \sum_g \int_0^\infty N(h)\,dF_{T,l,g}(h),
$$
this yields
$$
h^{\beta}\,N(h)^{-\sigma} \;=\; \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^{\sigma},
\qquad
\widetilde{H}_{T,l} \;\equiv\; \sum_g \int_0^\infty h^{\beta/\sigma}\,dF_{T,l,g}(h),
$$
where $F_{T,l,g}$ is the (endogenous) cdf of human capital among teachers working in $l$ from group $g$. Substituting back gives the reduced-form student technology
$$
\quad h_{i,g,l}(z,\epsilon_i) \;=\; \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^{\sigma}\,(z\,\epsilon_i)^{\alpha}\,s_{i,g,l}^{\phi}\,e_{i,g,l}^{\eta}.\quad
$$
For notational convenience define the *teacher-quality shifter* $Q_l \equiv (2\widetilde{H}_{T,l}/M_l)^{\sigma}$.

**Production.** Non-teaching output is $y_{i,g}=A_i\,h_{i,g,l}$, with $\omega_{i,l'}(h)=A_i h$ for $i\neq T$ — i.e., non-teaching wages depend on training location $l$ only through $h_{i,g,l}$, not on the work location $l'$. Teaching wages $\omega_{T,l'}(h)$ are taken as exogenous primitives in shape (strictly increasing, continuously differentiable). We adopt the parametric form
$$
\omega_{T,l'}(h) \;=\; \kappa_{l'}\,h^{\gamma},
$$
with location-specific intercepts $\kappa_{l'}$ and a common elasticity $\gamma$; location-specific tax rates $t_{l'}$ adjust endogenously to satisfy the local balanced-budget condition. A teacher trained in $l$ with human capital $h$ who works in $l'$ earns $\omega_{T,l'}(h)$.

## 5 Preferences

Let $\lambda \geq 0$ index the strength of altruism, and let $f(\cdot)$ be the warm-glow kernel (for example, $f(h')=\log h'$, so $\lambda\,f$ enters as $\lambda\,\mathbb E[\log h']$). A young agent born in $l$, group $g$, with persistent state $z$ and idiosyncratic vector $\vec\epsilon$, who chooses occupation $i$ when young, values
$$
U \;=\; \ln(1-s_{i,g,l}) \;+\; \mathbb{E}\!\left[\, \mu \ln C \;+\; B_{l'} - \tau_{l,l'} \;+\; \sigma_{\nu}\nu_{l'} \;+\; \lambda\,f(h')\,\right],
$$
where the expectation is over: the Gumbel shocks $\{\nu_{l'}\}$ (which determine $l'$), the child's gender $g' \in \{m,f\}$ (each with probability $\tfrac{1}{2}$), the child's persistent state $z' \sim \pi_z(\cdot\mid z)$, and the child's idiosyncratic vector $\vec\epsilon' \sim F_\epsilon^{\otimes I}$ — which determine $h'$, the child's optimal occupation $i'$, and education spending $e'_{i',g',l'}$. The moving cost $\tau_{l,l'}\geq 0$ is paid in utils; $\tau_{l,l}=0$.

Conditional on realized $l'$, $(z',\vec\epsilon')$, and child's group $g'$, old-age consumption is:
$$
\begin{aligned}
C(l',z',\vec\epsilon', g' \mid i,g,l,z,\vec\epsilon) \;=\;& (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l}(z,\epsilon_i)\big) \\
&-\;(1-\varepsilon)(1+\tau^{e}_{i,g})\,e_{i,g,l} \\
&-\;\varepsilon\,(1+\tau^{e}_{i'(z',\vec\epsilon',l',g'),g'})\,e'_{i'(z',\vec\epsilon',l',g'),g',l'}(z',\epsilon'_{i'}).
\end{aligned}
$$
The first net term is the parent's after-tax labor income. The second is her $(1-\varepsilon)$ share of her own childhood education. The third is the parental $\varepsilon$ share of her child's (now-realized) educational outlay; the child's optimal $i'$ and $e'$ are functions of $(z',\vec\epsilon',g')$ and the child's birth location $l'$, with corresponding aggregates $(\widetilde{H}_{T,l'},M_{l'},t_{l'})$. The parent's own labor-market distortion $(1-\tau^\omega_{i,g})$ and own educational barrier $(1+\tau^e_{i,g})$ retain $g$, the parent's gender; only the child's educational outlay carries the uncertain $g'$.

Setting $\varepsilon=0$, $\lambda=0$, $L=1$, and $\rho_z = 0$ recovers the baseline model.

## 6 Location Choice

Having chosen occupation $i$ when young, the agent observes $\{\nu_{l'}\}$ at Stage 2 and picks
$$
l^{*} \in \arg\max_{l' \in \{1,2\}}\; V_{l'}(h_{i,g,l}; i,g,l,z) + B_{l'} - \tau_{l,l'} + \sigma_\nu \nu_{l'},
$$
where $V_{l'}$ is the expected value of working in $l'$, integrating out the child's gender, persistent state, and idiosyncratic vector:
$$
V_{l'}(h; i,g,l,z) \;=\; \tfrac{1}{2}\sum_{g'\in\{m,f\}}\mathbb{E}_{z'\mid z}\,\mathbb{E}_{\vec\epsilon'}\!\left[\, \mu \ln C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon) \;+\; \lambda\,f\!\big(h'(z',\vec\epsilon',l',g')\big)\,\right].
$$
The moving cost $\tau_{l,l'}$ enters additively in the location-choice utility, parallel to the amenity $B_{l'}$. The Gumbel structure delivers closed-form choice probabilities,
$$
\pi_{l'\mid l}(h; i,g,z) \;=\; \frac{\exp\!\big[(V_{l'}(h; i,g,l,z) + B_{l'} - \tau_{l,l'})/\sigma_\nu\big]}{\sum_{l''}\exp\!\big[(V_{l''}(h; i,g,l,z) + B_{l''} - \tau_{l,l''})/\sigma_\nu\big]},
$$
and the log-sum value
$$
\bar V(h; i,g,l,z) \;=\; \sigma_\nu \ln\!\sum_{l'}\exp\!\big[(V_{l'}(h; i,g,l,z) + B_{l'} - \tau_{l,l'})/\sigma_\nu\big]
$$
folds the Stage-2 problem cleanly into Stage-1 occupational choice. As $\sigma_\nu \to 0$, $\pi_{l'\mid l}$ collapses to a deterministic indicator with $\tau_{l,l'}$ acting as a wedge in net values; as $\tau_{l,l'} \to \infty$ for $l'\neq l$, agents stay in their birth location regardless of $\nu$.

## 7 Occupation and Human Capital Choice

When young, the agent solves
$$
i^*(z,\vec\epsilon; g, l) \;=\; \arg\max_{i \in \{1,\dots,I\}}\; \Big[\ln(1 - s^*_{i,g,l}) + \bar V\!\big(h^*_{i,g,l}(z,\epsilon_i); i, g, l, z\big)\Big],
$$
with $\big(s^*_{i,g,l}, e^*_{i,g,l}\big)$ given by the FOCs below, taking the location-choice probabilities $\pi_{l'\mid l}$ and the child's policy functions as given. (Verbose derivations are in Appendix A.2–A.3.)

For non-teaching ($i \neq T$, with $\omega_{i,l'}(h)=A_i h$), the FOCs reduce to:

**FOC for $e_{i,g,l}$:**
$$
\sum_{l'} \pi_{l'\mid l}\;\tfrac{1}{2}\!\sum_{g'}\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{(1-t_{l'})(1-\tau^\omega_{i,g})A_i\eta\,h_{i,g,l}/e_{i,g,l} \;-\; (1-\varepsilon)(1+\tau^e_{i,g})}{C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon)}\right]\;=\;0.
$$

**FOC for $s_{i,g,l}$:**
$$
\frac{1}{1-s_{i,g,l}} \;=\; \mu\,\phi\,(1-\tau^\omega_{i,g})A_i\,\frac{h_{i,g,l}}{s_{i,g,l}}\,\sum_{l'} \pi_{l'\mid l}\;\tfrac{1}{2}\!\sum_{g'}\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{1-t_{l'}}{C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon)}\right].
$$

For teaching ($i=T$), the FOCs are analogous with $\partial\omega_{T,l'}(h)/\partial h$ replacing $A_i$ and $\omega_{T,l'}(h)$ replacing $A_i h$ in the consumption term — i.e., $A_i\eta h/e$ becomes $\eta\,(\partial \omega_{T,l'}/\partial h)\,h/e$ and the wage entering $C$ is $\omega_{T,l'}(h_{T,g,l})$. Because $\omega_{T,l'}$ depends on the work location $l'$, optimal $(s_{T,g,l}, e_{T,g,l})$ depend on the full vector $\{\pi_{l'\mid l}\}_{l'}$ rather than only own-location returns. For the power-form $\omega_{T,l'}(h)=\kappa_{l'} h^{\gamma}$ used in the code, $\partial\omega_{T,l'}/\partial h = \gamma\kappa_{l'}h^{\gamma-1}$.

Note two things. First, the moving cost $\tau_{l,l'}$ enters the FOCs only through the choice probabilities $\pi_{l'\mid l}$: $\tau_{l,l'}$ shifts the weights $\pi_{l'\mid l}$ on each location's return. Second, since $\pi_{l'\mid l}$ depends on $V_{l'}$, which depends on $C$, which depends on $(s,e)$, the FOCs constitute a fixed-point problem in $(s, e, \pi)$ for each $(i,g,l,z,\vec\epsilon)$.

## 8 Distributions and Aggregation

**Endogenous distribution of $z$ by location.** Let $\Phi_l(z)$ denote the joint mass of young agents in location $l$ with persistent ability $z$, gender-pooled. By the iid gender coin flip, the per-gender mass is symmetric: $\Phi_l^{g}(z) = \tfrac{1}{2}\Phi_l(z)$ for $g\in\{m,f\}$. The marginal of $z$ for young in $l$ is therefore $G_l(z) \equiv \Phi_l(\{z'\leq z\})/\big(M_l/2\big)$. The iid idiosyncratic shocks $\vec\epsilon$ are distributed exogenously according to $F_\epsilon^{\otimes I}$ within every cohort and location. So every population-weighted integral over $(z,\vec\epsilon)$ in $(g,l)$ factors as
$$
\int\!\!\int (\cdot)\, dF_\epsilon^{\otimes I}(\vec\epsilon)\, \Phi_l^{g}(dz) \;=\; \tfrac{1}{2}\int\!\!\int (\cdot)\, dF_\epsilon^{\otimes I}(\vec\epsilon)\, \Phi_l(dz),
$$
with total young-in-$(g,l)$ mass $\int\Phi_l^{g}(dz) = M_l/4$.

**Optimal-occupation indicator.** Let $\mathbb 1_{i,g,l}(z,\vec\epsilon) \equiv \mathbb 1[i^*(z,\vec\epsilon; g, l) = i]$.

**Mass flowing through each $(i,g,l)\to l'$.**
$$
m^t_{i,g,l, l'} \;=\; \tfrac{1}{2}\!\int\!\!\int \mathbb 1_{i,g,l}(z,\vec\epsilon)\,\pi_{l'\mid l}(h_{i,g,l}(z,\epsilon_i); i, g, z)\,dF_\epsilon^{\otimes I}(\vec\epsilon)\,\Phi_{l,t}(dz).
$$

**Aggregate teacher human capital** working in $l'$ at $t+1$ (used to teach the cohort born in $l'$ at $t+1$):
$$
\widetilde{H}_{T,l',t+1} \;=\; \tfrac{1}{2}\sum_g\sum_{l_0}\!\int\!\!\int \mathbb 1_{T,g,l_0}(z,\vec\epsilon)\,\pi_{l'\mid l_0}\!\big(h_{T,g,l_0}; T, g, z\big)\,h_{T,g,l_0}(z,\epsilon_T)^{\beta/\sigma}\,dF_\epsilon^{\otimes I}(\vec\epsilon)\,\Phi_{l_0,t}(dz).
$$
This expression uses the individual $(z,\vec\epsilon,g)$ policies and the endogenous joint mass $\Phi_{l_0}$ on the birth-location side, together with the exogenous iid $\vec\epsilon$ distribution.

**Migration probability of parents born in $l_0$ with persistent state $z$.** Integrating over occupation choice and idiosyncratic shocks at fixed gender,
$$
\bar\pi_{l_0 \to l'}^{g}(z) \;=\; \int \sum_i \mathbb 1_{i,g,l_0}(z,\vec\epsilon)\,\pi_{l'\mid l_0}(h_{i,g,l_0};i,g,z)\,dF_\epsilon^{\otimes I}(\vec\epsilon),
$$
and gender-pooling,
$$
\bar\pi_{l_0 \to l'}(z) \;=\; \tfrac{1}{2}\sum_g \bar\pi_{l_0 \to l'}^{g}(z).
$$

**Law of motion for the endogenous joint mass $\Phi_l$.** A parent in $(l_0,z)$ migrates to $l'=l$ with prob $\bar\pi_{l_0\to l}(z)$ and her child's persistent state is drawn from $\pi_z(\cdot\mid z)$:
$$
\Phi_{l,t+1}(z') \;=\; \sum_{l_0}\int \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)\,\Phi_{l_0,t}(dz).
$$

**Population flows.** Total population in $l'$ in $t+1$ is
$$
M_{l',t+1} \;=\; 2\!\int\!\Phi_{l',t+1}(dz) \;=\; \sum_{l_0}\,M_{l_0,t}\,\bar\Pi_{l_0,l',t},
\qquad
\bar\Pi_{l_0,l',t} \;=\; \frac{\int \bar\pi_{l_0\to l'}(z)\,\Phi_{l_0,t}(dz)}{\int \Phi_{l_0,t}(dz)},
$$
i.e., $\bar\Pi_{l_0,l',t}$ is the $\Phi_{l_0,t}$-weighted average over $z$ of the conditional migration probability. Rows of $[\bar\Pi_{l_0,l',t}]$ sum to 1. Young agents are evenly split by gender in every location by construction; the gender composition of old agents working in $l'$ is endogenous, determined by the gender-specific migration matrix $\{\bar\pi_{l_0\to l'}^{g}(z)\}$.

**Teacher distribution.** The (unnormalized) mass-weighted measure of human capital among teachers working in $l'$ from group $g$, $F_{T,l',g}(\cdot)$, is the push-forward of the joint distribution of $(z,\vec\epsilon,l_0)$ via the map $(z,\vec\epsilon,l_0) \mapsto h_{T,g,l_0}(z,\epsilon_T)$, weighted by $\mathbb 1_{T,g,l_0}\,\pi_{l'\mid l_0}$ and by the joint mass $\Phi_{l_0}(dz)\,dF_\epsilon$. The total mass of $F_{T,l',g}$ equals the number of teachers from group $g$ working in $l'$, so that the class-size constraint $M_l/2 = \sum_g\int N\,dF_{T,l,g}$ and the aggregator $\widetilde{H}_{T,l'} = \sum_g\int h^{\beta/\sigma}\,dF_{T,l',g}$ both balance. Only the moment $\widetilde{H}_{T,l'}$ enters the rest of the model.

## 9 Government Budget

Each location runs a separate balanced budget: payroll taxes on labor income earned in $l'$ fund teacher wages paid in $l'$. With taxable income
$$
\mathcal{I}_{l'} \;=\; \tfrac{1}{2}\sum_g\sum_{l_0}\,\sum_i\!\int\!\!\int \mathbb 1_{i,g,l_0}(z,\vec\epsilon)\,\pi_{l'\mid l_0}\,(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l_0}(z,\epsilon_i)\big)\,dF_\epsilon^{\otimes I}(\vec\epsilon)\,\Phi_{l_0}(dz)
$$
and teacher wage bill
$$
\mathcal{W}_{l'} \;=\; \tfrac{1}{2}\sum_g\sum_{l_0}\!\int\!\!\int \mathbb 1_{T,g,l_0}(z,\vec\epsilon)\,\pi_{l'\mid l_0}\,\omega_{T,l'}\!\big(h_{T,g,l_0}(z,\epsilon_T)\big)\,dF_\epsilon^{\otimes I}(\vec\epsilon)\,\Phi_{l_0}(dz),
$$
the local rate $t_{l'}$ solves
$$
t_{l'}\,\mathcal{I}_{l'} \;=\; \mathcal{W}_{l'},\qquad l' \in \{1,2\}.
$$
Both the wage bill and tax base sum across training locations $l_0$, since teachers and other workers can migrate, and they integrate over the endogenous distributions $\Phi_{l_0}$ on the birth-location side.

## 10 Stationary Equilibrium

A stationary equilibrium is a tuple
$$
\Big\{\,t_l,\; \widetilde{H}_{T,l},\; M_l,\; \Phi_l,\; F_{T,l,g},\; \big(s_{i,g,l},\, e_{i,g,l}\big)_{i,g},\; i^*(\,\cdot\,;g,l)\,\Big\}_{l\in\{1,2\},\,g}
$$
such that:

1. **(Optimization)** For each $(g,l,z,\vec\epsilon)$, $\big(s_{i,g,l}, e_{i,g,l}\big)$ satisfy the FOCs given $\big(\widetilde{H}_{T,1}, \widetilde{H}_{T,2}, M_1, M_2, t_1, t_2\big)$ and the location-choice probabilities $\pi_{l'\mid l}$. Occupational choice $i^*(z,\vec\epsilon;g,l) = \arg\max_i\big[\ln(1-s_{i,g,l}) + \bar V(h_{i,g,l};i,g,l,z)\big]$.

2. **(Endogenous spatial distribution)** $\Phi_l$ is a stationary fixed point of the law of motion in §8:
$$
\Phi_l(z') \;=\; \sum_{l_0}\int \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)\,\Phi_{l_0}(dz),
$$
and total population $M_l = 2\int\Phi_l(dz)$ is consistent with the row-stochastic migration matrix $[\bar\Pi_{l_0,l'}]$. Young agents are evenly split by gender ($\tfrac{1}{2}$ each) in every location by construction.

3. **(Teacher distribution and aggregates)** $F_{T,l,g}$ is the push-forward described in §8, and
$$
\widetilde{H}_{T,l} \;=\; \sum_g \int h^{\beta/\sigma}\,dF_{T,l,g}(h)
$$
coincides with the explicit aggregator in §8 (using policies and $\Phi$).

4. **(Class size)** $h^{\beta} N(h)^{-\sigma} = (2\widetilde{H}_{T,l}/M_l)^{\sigma}$ for all teachers in $l$, with $\sum_g\int N\,dF_{T,l,g} = M_l/2$.

5. **(Budget)** Each $t_l$ satisfies the local balanced-budget condition $t_{l}\mathcal{I}_{l}=\mathcal{W}_{l}$.

---

## Appendix A — Derivations

This appendix collects verbose derivations of the results stated in the main text.

### A.1 Reduced-form human capital and the teacher-quality shifter $Q_l$

**Setup.** The raw technology states that a young agent of group $g$ in $l$ with abilities $(z,\vec\epsilon)$, who chooses occupation $i$, $(s,e)$, and is matched with a teacher of human capital $h_{T,l}$, accumulates
$$
h \;=\; h_{T,l}^{\beta}\,(z\epsilon_i)^{\alpha}\,s^{\phi}\,e^{\eta}\,N(h_{T,l})^{-\sigma}, \tag{A.1}
$$
where $N(h_{T,l})$ is the class size the teacher is willing to take. The dependence on $N$ captures congestion in instruction.

**Teacher's indifference condition.** A teacher of human capital $h_{T,l}$ working in $l$ is paid $\omega_{T,l}(h_{T,l})$ regardless of class size, but exerts effort that disutility scales with $N$. In equilibrium, all active teachers in $l$ must be indifferent over choosing different $N$. Standard arguments (see baseline derivation) yield that the locally relevant student-side product $h^{\beta}N(h)^{-\sigma}$ is constant across all $h$ among active teachers in $l$:
$$
h_{T,l}^{\beta}\,N(h_{T,l})^{-\sigma} \;=\; \mathcal{Q}_l \quad \text{(constant in $l$).} \tag{A.2}
$$
This is the indifference condition referenced in §4.

**Pinning down $\mathcal Q_l$ via the class-size resource constraint.** Total student–teacher slots in $l$ must absorb the local young cohort of size $M_l/2$:
$$
\frac{M_l}{2} \;=\; \sum_g\int_0^\infty N(h)\,dF_{T,l,g}(h). \tag{A.3}
$$
From (A.2), $N(h) = h^{\beta/\sigma}\,\mathcal Q_l^{-1/\sigma}$. Substituting into (A.3):
$$
\frac{M_l}{2} \;=\; \mathcal Q_l^{-1/\sigma}\!\sum_g\!\int h^{\beta/\sigma}\,dF_{T,l,g}(h) \;=\; \mathcal Q_l^{-1/\sigma}\,\widetilde{H}_{T,l},
$$
so
$$
\mathcal Q_l \;=\; \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^{\sigma} \;\equiv\; Q_l. \tag{A.4}
$$

**Reduced form.** Substituting (A.2)–(A.4) into (A.1):
$$
h_{i,g,l}(z,\epsilon_i) \;=\; \big[h_{T,l}^{\beta}N(h_{T,l})^{-\sigma}\big]\,(z\epsilon_i)^{\alpha}\,s_{i,g,l}^{\phi}\,e_{i,g,l}^{\eta} \;=\; Q_l\,(z\epsilon_i)^{\alpha}\,s_{i,g,l}^{\phi}\,e_{i,g,l}^{\eta}, \tag{A.5}
$$
matching the boxed expression in §4. The dependence on the specific teacher drops out, because the indifference condition forces $h^{\beta}N^{-\sigma}=Q_l$ for every $h$ in the support of $F_{T,l,g}$.

### A.2 FOCs for non-teaching occupations

Fix an occupation $i\neq T$ with $\omega_{i,l'}(h)=A_i h$. The young agent's expected utility, given choice of occupation $i$, is
$$
\mathcal U \;=\; \ln(1-s) \;+\; \mathbb E\!\left[\mu\ln C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon) + B_{l'} - \tau_{l,l'} + \sigma_\nu\nu_{l'} + \lambda f(h')\right].
$$
Taking expectations explicitly over $\nu$ (Gumbel) gives the log-sum $\bar V$; what remains, conditional on $l'$, is
$$
\mathcal{U}(s,e) \;=\; \ln(1-s) \;+\; \sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\Big[\mu\ln C(l',z',\vec\epsilon',g') + \lambda f(h'(z',\vec\epsilon',l',g'))\Big] + \text{const}_{l',B,\tau,\nu}.
$$
The continuation $\lambda f(h')$ does not depend on the parent's $(s,e)$ — it depends on the child's policies — so the FOCs for $(s,e)$ involve only the $\mu\ln C$ term.

**FOC for $e$.** Differentiating $\mu\ln C$ with respect to $e\equiv e_{i,g,l}$:
$$
\frac{\partial}{\partial e}\,\mu\ln C \;=\; \frac{\mu}{C}\cdot\frac{\partial C}{\partial e}, \qquad
\frac{\partial C}{\partial e} \;=\; (1-t_{l'})(1-\tau^{\omega}_{i,g})A_i\,\frac{\partial h}{\partial e} \;-\; (1-\varepsilon)(1+\tau^{e}_{i,g}). \tag{A.6}
$$
Using $h = Q_l (z\epsilon_i)^\alpha s^\phi e^\eta$, $\partial h/\partial e = \eta h/e$. Setting the sum across $(l',g',z',\vec\epsilon')$ to zero:
$$
\boxed{\quad
\sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{(1-t_{l'})(1-\tau^{\omega}_{i,g})A_i\eta\,h/e \;-\; (1-\varepsilon)(1+\tau^{e}_{i,g})}{C(l',z',\vec\epsilon',g')}\right] \;=\; 0.
\quad}
\tag{A.7}
$$

**FOC for $s$.** Differentiate $\ln(1-s) + \mu\ln C$:
$$
-\,\frac{1}{1-s} \;+\; \sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{\mu}{C}\cdot\frac{\partial C}{\partial s}\right] \;=\; 0.
$$
Using $\partial C/\partial s = (1-t_{l'})(1-\tau^{\omega}_{i,g})A_i\,\partial h/\partial s$ and $\partial h/\partial s = \phi\, h/s$, we get
$$
\boxed{\quad
\frac{1}{1-s} \;=\; \mu\,\phi\,(1-\tau^{\omega}_{i,g})A_i\,\frac{h}{s}\,\sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{1-t_{l'}}{C(l',z',\vec\epsilon',g')}\right].
\quad}
\tag{A.8}
$$

**Joint system.** Equations (A.7) and (A.8), together with $h=Q_l(z\epsilon_i)^\alpha s^\phi e^\eta$ and the location-choice fixed point in $\{\pi_{l'\mid l}\}$, form a coupled system for $(s,e,\pi)$ at each $(i,g,l,z,\vec\epsilon)$.

### A.3 FOCs for teaching

For $i=T$ the wage $\omega_{T,l'}(h)$ depends on the *work* location $l'$ (and is non-linear in $h$). The derivations parallel A.2 but with two changes:

(i) Wherever $A_i h$ appeared in $C$, replace with $\omega_{T,l'}(h_{T,g,l}(z,\epsilon_T))$. In particular,
$$
C \;=\; (1-t_{l'})(1-\tau^{\omega}_{T,g})\,\omega_{T,l'}(h_{T,g,l}) \;-\; (1-\varepsilon)(1+\tau^{e}_{T,g})e_{T,g,l} \;-\; \varepsilon(1+\tau^{e}_{i',g'})e'_{i',g',l'}(z',\epsilon_{i'}').
$$

(ii) Wherever $A_i \partial h/\partial x$ appeared in $\partial C/\partial x$, replace with $\omega_{T,l'}'(h)\,\partial h/\partial x$.

Concretely, with $h = Q_l(z\epsilon_T)^\alpha s^\phi e^\eta$:
$$
\frac{\partial C}{\partial e}\bigg|_{l'} \;=\; (1-t_{l'})(1-\tau^{\omega}_{T,g})\,\omega_{T,l'}'(h)\,\eta\frac{h}{e} \;-\; (1-\varepsilon)(1+\tau^{e}_{T,g}), \tag{A.9}
$$
$$
\frac{\partial C}{\partial s}\bigg|_{l'} \;=\; (1-t_{l'})(1-\tau^{\omega}_{T,g})\,\omega_{T,l'}'(h)\,\phi\frac{h}{s}. \tag{A.10}
$$
Substituting into the analogues of (A.7) and (A.8) yields the teaching FOCs. For the power form $\omega_{T,l'}(h)=\kappa_{l'}h^{\gamma}$, $\omega_{T,l'}'(h)=\gamma\kappa_{l'}h^{\gamma-1}$, so $\omega_{T,l'}'(h)\,h = \gamma\omega_{T,l'}(h)$, and the FOCs become
$$
\sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{(1-t_{l'})(1-\tau^{\omega}_{T,g})\gamma\eta\,\kappa_{l'}h^{\gamma}/e \;-\; (1-\varepsilon)(1+\tau^{e}_{T,g})}{C}\right] \;=\; 0, \tag{A.11}
$$
$$
\frac{1}{1-s} \;=\; \mu\,\phi\,\gamma\,(1-\tau^{\omega}_{T,g})\,\frac{1}{s}\,\sum_{l'}\pi_{l'\mid l}\,\tfrac{1}{2}\!\sum_{g'}\,\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{(1-t_{l'})\kappa_{l'}h^{\gamma}}{C}\right]. \tag{A.12}
$$

### A.4 Single-location, no-altruism limit: closed forms for $s$

Set $L=1$ (drop the $\sum_{l'}\pi_{l'\mid l}$), $\varepsilon=0$, $\lambda=0$, $\rho_z=0$. Drop conditioning on $(z',\vec\epsilon',g')$ since $C$ no longer depends on the child. We get
$$
C \;=\; (1-t)(1-\tau^{\omega}_{i,g})\omega_{i}(h) \;-\; (1+\tau^{e}_{i,g})e,
$$
and the FOCs in A.2 reduce to
$$
(1-t)(1-\tau^\omega)A_i\eta h/e \;=\; (1+\tau^e), \qquad \frac{1}{1-s} \;=\; \mu\phi(1-\tau^\omega)A_i\,\frac{h}{s\,C}.
$$
From the $e$-FOC: $(1-t)(1-\tau^\omega)A_i h = (1+\tau^e)e/\eta$. Substituting into $C$:
$$
C \;=\; \frac{(1+\tau^e)e}{\eta} \;-\; (1+\tau^e)e \;=\; (1+\tau^e)e\,\frac{1-\eta}{\eta}.
$$
Hence $(1-t)(1-\tau^\omega)A_i h / C = 1/(1-\eta)$ and the $s$-FOC becomes
$$
\frac{1}{1-s} \;=\; \frac{\mu\phi}{s(1-\eta)} \quad\Longrightarrow\quad s\,(1-\eta) \;=\; \mu\phi\,(1-s) \quad\Longrightarrow\quad
\boxed{\;s_O^{*} \;=\; \frac{\mu\phi}{\mu\phi + 1 - \eta}.\;}
$$
This matches `s_O_const(p) = μϕ/(μϕ + 1 - η)` in `spatial.jl`.

For teaching with $\omega_T(h)=\kappa h^{\gamma}$, the analogous algebra gives $(1-t)(1-\tau^\omega)\kappa h^\gamma / C = 1/(1-\gamma\eta)$ and
$$
\frac{1}{1-s} \;=\; \frac{\mu\phi\gamma}{s(1-\gamma\eta)} \quad\Longrightarrow\quad
\boxed{\;s_T^{*} \;=\; \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1}.\;}
$$
This matches `s_T_const(p) = μϕγ/((μϕ-η)γ + 1)` in `spatial.jl`. Both constants are independent of $(z,\vec\epsilon,l)$ in this limit; with migration, altruism, and parental finance, the full FOCs are needed.

### A.5 Stationarity of $\Phi_l$

From §8, $\Phi_{l,t+1}(z') = \sum_{l_0}\int \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)\,\Phi_{l_0,t}(dz)$. Integrating in $z'$ and using $\int\pi_z(z'\mid z)\,dz' = 1$:
$$
\int\Phi_{l,t+1}(dz') \;=\; \sum_{l_0}\int \bar\pi_{l_0\to l}(z)\,\Phi_{l_0,t}(dz),
$$
which gives the population law $M_{l,t+1} = \sum_{l_0}M_{l_0,t}\bar\Pi_{l_0,l,t}$ via $M_l = 2\int\Phi_l$. In a stationary equilibrium the joint mass $\Phi_l(\cdot)$ is the dominant left eigenfunction (eigenvalue 1) of the joint kernel $K_{(l_0,z)\to(l,z')} = \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)$ on $\{1,\ldots,L\}\times\mathbb{R}_+$.

---

## Appendix B — Computational Notes

I don't believe we recover closed-form analogues of baseline Proposition 2 (that occupational choice does not depend on aggregate teacher human capital $\widetilde H$) because teaching wages depend on the work location $l'$ rather than the training location $l$, so the prospective teacher's expected return is a $\pi_{l'\mid l}$-weighted average across locations rather than a single-location expression; and altruism plus the $\varepsilon$-share of child education couples the parent's consumption to the child's optimal $e'_{i'(z',\vec\epsilon',l',g'),g',l'}(z',\epsilon'_{i'})$, which is itself the solution of a problem in the child's location $l'$ with aggregates $(\widetilde{H}_{T,l'}, M_{l'}, t_{l'})$.

Coding the model amounts to four nested fixed points:

- **Inner:** for each $(z,\vec\epsilon, g, l, i)$, jointly solve the FOCs and the choice probabilities for $(s_{i,g,l}, e_{i,g,l}, \pi_{l'\mid l})$ given aggregates and the child's policy.
- **Child consistency:** the child's policy $\big(s'_{i',g',l'}, e'_{i',g',l'}\big)$ entering $C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon)$ must coincide with the equilibrium policy of a gender-$g'$ agent in stationary equilibrium, for each $g' \in \{m,f\}$.
- **Aggregation:** integrate individual policies against the *endogenous* joint mass $\Phi_{l}$ on the birth-location side and the *exogenous* iid $\epsilon$ distribution to compute $\widetilde{H}_{T,l}$, $M_l$, and gender shares.
- **Budget:** find $t_l$ satisfying local budget for each $l$.

A natural outer loop iterates on $\big(\Phi_1,\Phi_2,\widetilde{H}_{T,1},\widetilde{H}_{T,2},M_1,M_2,t_1,t_2\big)$ until all flow, aggregation, and budget equations balance. In `julia/spatial_model/spatial.jl` the joint mass $\Phi_l$ is represented as the matrix `m_young[l, k]`, a discrete approximation indexed by location and Tauchen node for $z$.

### B.1 Simplifications used in `spatial.jl`

The Julia implementation introduces one structural simplification relative to the exact model in the main text. All other aspects of the model are fully and exactly implemented.

**Two-draw idiosyncratic shocks.** The model in §3 specifies $I$ i.i.d. draws $(\epsilon_1,\ldots,\epsilon_I)$, one per occupation. Integrating over the full $I$-dimensional joint distribution is computationally expensive. The code instead uses two independent draws: $\epsilon_T \sim F_\epsilon$ for teaching and a single $\epsilon_O \sim F_\epsilon$ applied uniformly to all non-teaching occupations. Within the teaching-vs-non-teaching comparison the threshold method integrates over $(\epsilon_T, \epsilon_O)$ jointly — the outer loop in `state_moments` runs over $\epsilon_T$ nodes, and for each $\epsilon_T$ the function `invert_threshold` finds the $\epsilon_O^*$ such that $U_{O^*}(\epsilon_O^*) = U_T(\epsilon_T)$, so the integral over $(\epsilon_T, \epsilon_O)$ is exact conditional on this two-draw structure. The optimal non-teaching occupation at a given $\epsilon_O$ is $\arg\max_i A_i \cdot h(\epsilon_O)$, which in practice is determined by the occupation-specific productivities $\{A_i\}$ and distortions $\{\tau^\omega_{i,g}\}$ rather than by $\epsilon_O$ itself.

**The following aspects of the model are fully implemented without simplification:**

- **State-dependent time investment $s$.** Both `solve_T` and `solve_O` solve the exact FOCs (A.8) and (A.12) for $s$ via bisection at each $(z,\epsilon,g,l)$ point. The closed-form constants from Appendix A.4 serve only as initial guesses for the bisection; the converged $s$ is fully state-dependent.

- **Goods-investment FOC.** The FOCs (A.7) and (A.11) for $e$ are implemented exactly, including the $(1-\varepsilon)$ factor on the RHS and the child's education costs inside $C$ (through `expected_child_terms`, which integrates $1/C$ over $(z',\epsilon_T',\epsilon_O',g')$ using the threshold method).

- **Expected log utility.** Old-age utility $\mu\,\mathbb E[\ln C]$ is computed by explicit integration over the child's $(\epsilon_T',\epsilon_O')$ draws in `expected_logC_invC`, using the same threshold-histogram method as the occupation-choice aggregation. The Jensen approximation $\mu\ln(\mathbb E[C])$ is not used.

- **Non-teacher taxable income.** The contribution of non-teachers to the tax base in location $l'$ is computed as the joint expectation $\mathbb E_{\epsilon_T,\epsilon_O}[\mathbf 1[\text{choose O}]\cdot\pi_{l'\mid l}(\epsilon_O)\cdot(1-\tau^\omega_O)\,w(\epsilon_O)]$ via `hist_weighted_sum_gt_joint`. This avoids the otherwise tempting product of separately computed conditional means $\mathbb E[\pi_{l'}|\epsilon_O>\epsilon^*]\cdot\mathbb E[\text{taxO}_{l'}|\epsilon_O>\epsilon^*]$, which ignores the positive correlation between $\pi_{l'}$ and taxable wages through $\epsilon_O$.

- **Teaching distortions.** The `Params` struct carries gender-specific labor-market distortion $\tau^\omega_{T,g}$ (`τw_T`) and educational barrier $\tau^e_{T,g}$ (`τe_T`) for teaching. Both enter `solve_T` (consumption and FOCs), `state_moments` (the child's reported cost `costT` and taxable income `taxableT`), and the budget aggregation. In the current baseline calibration they are set to zero for both genders, but the infrastructure is in place to set them non-trivially.

### B.2 Outer iteration in the code

The outer loop in `solve_stationary` damps updates of $(m_{\text{young}}, M, \widetilde H, t, E[\log h], E[\text{cost}])$ separately, with damping weights `damp_m, damp_H, damp_t, damp_E`. A period-2 cycle detector is included: if the iterates revisit the same point every two steps (a common pathology of threshold/histogram maps), the average of the two phases is returned. The convergence criterion is a sup-norm error across $m_{\text{young}}$, $\widetilde H$, and $t$.
