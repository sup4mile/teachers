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

**Stage 1 — Young.** Born in location $l$ (chosen by the parent in Stage 2 of the previous cohort), with persistent ability $z$ inherited from the parent. Given $l$ and $z$:

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
with $h_{T,l}$ the human capital of the (randomly matched) teacher working in $l$. Student-teacher matching is random within $l$. Class sizes are allocated efficiently — to maximize aggregate student human capital subject to the local slot constraint — which yields the equalization condition $h^{\beta} N(h)^{-\sigma}=\text{const}$ across all active teachers in $l$ (derived in Appendix A.1). Combined with the local class-size resource constraint
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
For notational convenience define the teacher-quality shifter $Q_l \equiv (2\widetilde{H}_{T,l}/M_l)^{\sigma}$.

**Production.** Non-teaching output is $y_{i,g}=A_i\,h_{i,g,l}$, with $\omega_{i,l'}(h)=A_i h$ for $i\neq T$ — i.e., non-teaching wages depend on training location $l$ only through $h_{i,g,l}$, not on the work location $l'$. Teaching wages $\omega_{T,l'}(h)$ are taken as exogenous primitives in shape (strictly increasing, continuously differentiable). We adopt the parametric form
$$
\omega_{T,l'}(h) \;=\; \kappa_{l'}\,h^{\gamma},
$$
with location-specific intercepts $\kappa_{l'}$ and a common elasticity $\gamma$; location-specific tax rates $t_{l'}$ adjust endogenously to satisfy the local balanced-budget condition. A teacher trained in $l$ with human capital $h$ who works in $l'$ earns $\omega_{T,l'}(h)$. These specifics are more a question of calibration and are subject to change.

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
with $\big(s^*_{i,g,l}, e^*_{i,g,l}\big)$ given by the FOCs below, taking the location-choice probabilities $\pi_{l'\mid l}$ and the child's policy functions as given. (Derivations are in the appendix.) The *occupational threshold* separating occupations $i$ and $j$ is the level set $\{(z,\vec\epsilon): W_{i,g,l}=W_{j,g,l}\}$ of the bracketed value $W_{i,g,l}$ — the spatial counterpart of the baseline teacher/other cutoff $a^*_{T,g}(a_O)$. In the full model this threshold is only an implicit locus; it admits a closed form solely in the $\varepsilon=0$, destination-homogeneous special case derived in Appendix A.8.

The investment policies are functions of the agent's state: $s^*_{i,g,l}=s^*_{i,g,l}(z,\epsilon_i)$ and $e^*_{i,g,l}=e^*_{i,g,l}(z,\epsilon_i)$. We suppress these arguments for readability, but they are not constant across $(z,\epsilon_i)$ in the spatial model. The closed forms in Appendix A.4 are constant only in the single-location limit, where the $\sum_{l'}\pi_{l'\mid l}$ collapses; with $L=2$ the FOCs aggregate returns across locations using $\pi_{l'\mid l}$ weights that themselves depend on $h_{i,g,l}(z,\epsilon_i)$.

For non-teaching ($i \neq T$, with $\omega_{i,l'}(h)=A_i h$), the FOCs reduce to:

**FOC for $e_{i,g,l}$:**
$$
\sum_{l'} \pi_{l'\mid l}\;\tfrac{1}{2}\!\sum_{g'}\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{(1-t_{l'})(1-\tau^\omega_{i,g})A_i\eta\,h_{i,g,l}/e_{i,g,l} \;-\; (1-\varepsilon)(1+\tau^e_{i,g})}{C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon)}\right]\;=\;0.
$$

**FOC for $s_{i,g,l}$:**
$$
\frac{1}{1-s_{i,g,l}} \;=\; \mu\,\phi\,(1-\tau^\omega_{i,g})A_i\,\frac{h_{i,g,l}}{s_{i,g,l}}\,\sum_{l'} \pi_{l'\mid l}\;\tfrac{1}{2}\!\sum_{g'}\mathbb E_{z'\mid z}\,\mathbb E_{\vec\epsilon'}\!\left[\frac{1-t_{l'}}{C(l',z',\vec\epsilon',g' \mid i,g,l,z,\vec\epsilon)}\right].
$$

For teaching ($i=T$), the FOCs are analogous with $\partial\omega_{T,l'}(h)/\partial h$ replacing $A_i$ and $\omega_{T,l'}(h)$ replacing $A_i h$ in the consumption term — i.e., $A_i\eta h/e$ becomes $\eta\,(\partial \omega_{T,l'}/\partial h)\,h/e$ and the wage entering $C$ is $\omega_{T,l'}(h_{T,g,l})$. Because $\omega_{T,l'}$ depends on the work location $l'$, optimal $(s_{T,g,l}, e_{T,g,l})$ depend on the full vector $\{\pi_{l'\mid l}\}_{l'}$ through both the wage intercepts $\kappa_{l'}$ and the tax rates $t_{l'}$. (Non-teaching policies likewise depend on $\{\pi_{l'\mid l}\}_{l'}$, but only through the $t_{l'}$ — and, when $\varepsilon>0$, the child-education — terms in $C$, since $\omega_{i,l'}=A_i h$ is independent of $l'$.) For the power-form $\omega_{T,l'}(h)=\kappa_{l'} h^{\gamma}$ used in the code, $\partial\omega_{T,l'}/\partial h = \gamma\kappa_{l'}h^{\gamma-1}$.

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
\Big\{\,t_l,\; \widetilde{H}_{T,l},\; M_l,\; \Phi_l,\; F_{T,l,g},\; \big(s_{i,g,l}(z,\epsilon_i),\, e_{i,g,l}(z,\epsilon_i)\big)_{i,g},\; i^*(\,\cdot\,;g,l)\,\Big\}_{l\in\{1,2\},\,g}
$$
such that:

1. **(Optimization)** For each $(g,l,z,\vec\epsilon)$, $\big(s_{i,g,l}, e_{i,g,l}\big)$ satisfy the FOCs given $\big(\widetilde{H}_{T,1}, \widetilde{H}_{T,2}, M_1, M_2, t_1, t_2\big)$ and the location-choice probabilities $\pi_{l'\mid l}$. Occupational choice $i^*(z,\vec\epsilon;g,l) = \arg\max_i\big[\ln(1-s_{i,g,l}) + \bar V(h_{i,g,l};i,g,l,z)\big]$.

2. **(Endogenous spatial distribution)** $\Phi_l$ is a stationary fixed point of the law of motion:
$$
\Phi_l(z') \;=\; \sum_{l_0}\int \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)\,\Phi_{l_0}(dz),
$$
and total population $M_l = 2\int\Phi_l(dz)$ is consistent with the row-stochastic migration matrix $[\bar\Pi_{l_0,l'}]$. Young agents are evenly split by gender ($\tfrac{1}{2}$ each) in every location by construction.

3. **(Teacher distribution and aggregates)** $F_{T,l,g}$ is the push-forward, and
$$
\widetilde{H}_{T,l} \;=\; \sum_g \int h^{\beta/\sigma}\,dF_{T,l,g}(h)
$$
coincides with the explicit aggregator (using policies and $\Phi$).

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

**Efficient class-size allocation.** A teacher of human capital $h_{T,l}$ working in $l$ is paid $\omega_{T,l}(h_{T,l})$ regardless of class size. Class sizes $\{N_j\}$ are allocated across the active teachers $j$ in $l$ to maximize aggregate student human capital subject to the local slot constraint $\sum_j N_j = M_l/2$. Because students are randomly matched, each teacher's class contributes $N_j$ students whose human capital is proportional to $h_j^{\beta}N_j^{-\sigma}$ times the common (matching-independent) average of $(z\epsilon)^\alpha s^\phi e^\eta$; the planner therefore maximizes $\sum_j N_j^{1-\sigma}h_j^{\beta}$ subject to $\sum_j N_j = M_l/2$. The FOC $(1-\sigma)N_j^{-\sigma}h_j^{\beta}=\text{const}$ equalizes the student-side product $h^{\beta}N(h)^{-\sigma}$ across all active teachers in $l$:
$$
h_{T,l}^{\beta}\,N(h_{T,l})^{-\sigma} \;=\; \mathcal{Q}_l \quad \text{(constant in $l$).} \tag{A.2}
$$
(Equivalently, this is the allocation a competitive market for school slots would implement; it requires no class-size disutility in preferences, consistent with §5.)

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
matching the expression in section 4. The dependence on the specific teacher drops out, because the indifference condition forces $h^{\beta}N^{-\sigma}=Q_l$ for every $h$ in the support of $F_{T,l,g}$.

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

For teaching with $\omega_T(h)=\kappa h^{\gamma}$, the analogous algebra gives $(1-t)(1-\tau^\omega)\kappa h^\gamma / C = 1/(1-\gamma\eta)$ and
$$
\frac{1}{1-s} \;=\; \frac{\mu\phi\gamma}{s(1-\gamma\eta)} \quad\Longrightarrow\quad
\boxed{\;s_T^{*} \;=\; \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1}.\;}
$$
Both constants are independent of $(z,\vec\epsilon,l)$ in this limit; with migration, altruism, and parental finance, the full FOCs are needed.

### A.5 Stationarity of $\Phi_l$

From section 8, $\Phi_{l,t+1}(z') = \sum_{l_0}\int \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)\,\Phi_{l_0,t}(dz)$. Integrating in $z'$ and using $\int\pi_z(z'\mid z)\,dz' = 1$:
$$
\int\Phi_{l,t+1}(dz') \;=\; \sum_{l_0}\int \bar\pi_{l_0\to l}(z)\,\Phi_{l_0,t}(dz),
$$
which gives the population law $M_{l,t+1} = \sum_{l_0}M_{l_0,t}\bar\Pi_{l_0,l,t}$ via $M_l = 2\int\Phi_l$. In a stationary equilibrium the joint mass $\Phi_l(\cdot)$ is the dominant left eigenfunction (eigenvalue 1) of the joint kernel $K_{(l_0,z)\to(l,z')} = \pi_z(z'\mid z)\,\bar\pi_{l_0\to l}(z)$ on $\{1,\ldots,L\}\times\mathbb{R}_+$.

### A.6 The special case $\varepsilon = 0$: simplifications and a closed form for $s$

The main text carries the parental cost share $\varepsilon \in (0,1)$ throughout. This appendix records what the model collapses to when $\varepsilon = 0$ — parents finance none of their child's education and each agent pays the full cost of her own schooling out of future labor income. The upshot: the optimal time share $s$ admits a closed form for *every* $(i,g,l,z,\vec\epsilon)$, regardless of $L$, the location-specific wedges, or the altruism strength $\lambda$; the goods margin $e$ and the spatial objects do not collapse.

**Consumption decouples from the child.** With $\varepsilon = 0$ the parental-finance term in the consumption identity of §5 vanishes, and old-age consumption conditional on the work location $l'$ is
$$
C_{l'} \;=\; (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}(h_{i,g,l}) \;-\; (1+\tau^{e}_{i,g})\,e_{i,g,l},
$$
which is **deterministic** given the parent's own state and $(s,e)$: it no longer depends on the child's realized $(z',\vec\epsilon',g')$. The Stage-2 location value therefore separates additively,
$$
V_{l'}(h;i,g,l,z) \;=\; \mu \ln C_{l'}(h) \;+\; \Lambda_{l'}(z),
\qquad
\Lambda_{l'}(z) \equiv \tfrac{1}{2}\sum_{g'}\mathbb{E}_{z'\mid z}\,\mathbb{E}_{\vec\epsilon'}\!\left[\lambda\, f\!\big(h'(z',\vec\epsilon',l',g')\big)\right],
$$
where the warm-glow term $\Lambda_{l'}(z)$ depends on the child's birth location $l'$ and (through the AR(1)) on the parent's $z$, but **not** on the parent's $(s,e)$. Consequently $\Lambda_{l'}$ differentiates to zero in the $(s,e)$ first-order conditions: altruism affects investment only indirectly, through the choice probabilities $\pi_{l'\mid l}$ that embed $\Lambda_{l'}$ via $V_{l'}$. Setting $\varepsilon=0$ in (A.7)–(A.8) and dropping the $\tfrac12\sum_{g'}\mathbb{E}_{z'}\mathbb{E}_{\vec\epsilon'}$ averaging (the summand no longer depends on the child), the non-teaching FOCs reduce to
$$
\sum_{l'}\pi_{l'\mid l}\;\frac{(1-t_{l'})(1-\tau^{\omega}_{i,g})A_i\,\eta\, h/e \;-\; (1+\tau^{e}_{i,g})}{C_{l'}} \;=\; 0,
\qquad
\frac{1}{1-s} \;=\; \mu\,\phi\,(1-\tau^{\omega}_{i,g})A_i\,\frac{h}{s}\sum_{l'}\pi_{l'\mid l}\,\frac{1-t_{l'}}{C_{l'}}.
$$

**Setup for the closed form.** Fix a state $(i,g,l,z,\vec\epsilon)$ and write $C_{l'} = g_{l'} - (1+\tau^e_{i,g})e$, where the after-tax gross income is
$$
g_{l'} \equiv (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}(h),
\qquad h = Q_l\,(z\epsilon_i)^\alpha s^\phi e^\eta .
$$
For non-teaching $\omega_{i,l'}(h)=A_i h$, so $g_{l'}$ varies across $l'$ only through $1-t_{l'}$; for teaching $\omega_{T,l'}=\kappa_{l'}h^\gamma$, it varies through both $1-t_{l'}$ and $\kappa_{l'}$.

**Step 1 — collapse the $e$-FOC.** Take non-teaching first. The marginal product of $e$ in income is $\partial g_{l'}/\partial e = \eta\,g_{l'}/e$ (since $h\propto e^\eta$). The $e$-FOC is
$$
\sum_{l'}\pi_{l'\mid l}\,\frac{\eta\,g_{l'}/e - (1+\tau^e_{i,g})}{C_{l'}} \;=\; 0 .
$$
Write $g_{l'} = C_{l'} + (1+\tau^e_{i,g})e$, so $\eta g_{l'}/e = \eta C_{l'}/e + \eta(1+\tau^e_{i,g})$ and the numerator becomes $\eta C_{l'}/e - (1-\eta)(1+\tau^e_{i,g})$. Dividing by $C_{l'}$ and summing with $\sum_{l'}\pi_{l'\mid l}=1$,
$$
\frac{\eta}{e} \;=\; (1-\eta)(1+\tau^e_{i,g})\sum_{l'}\pi_{l'\mid l}\,\frac{1}{C_{l'}}
\quad\Longrightarrow\quad
\sum_{l'}\pi_{l'\mid l}\,\frac{1}{C_{l'}} \;=\; \frac{\eta}{(1-\eta)(1+\tau^e_{i,g})\,e}. \tag{$\star$}
$$

**Step 2 — substitute into the $s$-FOC.** The marginal product of $s$ gives $\partial g_{l'}/\partial s = \phi\,g_{l'}/s$, so the $s$-FOC is
$$
\frac{1}{1-s} \;=\; \mu\phi\,\frac{1}{s}\sum_{l'}\pi_{l'\mid l}\,\frac{g_{l'}}{C_{l'}} .
$$
Now $g_{l'}/C_{l'} = 1 + (1+\tau^e_{i,g})e/C_{l'}$, so by $(\star)$,
$$
\sum_{l'}\pi_{l'\mid l}\,\frac{g_{l'}}{C_{l'}}
= 1 + (1+\tau^e_{i,g})e\sum_{l'}\pi_{l'\mid l}\frac{1}{C_{l'}}
= 1 + \frac{\eta}{1-\eta}
= \frac{1}{1-\eta}.
$$
Every location-specific quantity ($t_{l'},\,\pi_{l'\mid l},\,\tau$, the level $h$) has cancelled. The $s$-FOC is $\tfrac{1}{1-s} = \tfrac{\mu\phi}{s(1-\eta)}$, i.e. $s(1-\eta)=\mu\phi(1-s)$, giving
$$
\boxed{\;s_O^{*} \;=\; \frac{\mu\phi}{\mu\phi + 1 - \eta}\;}.
$$

**Teaching.** Identical steps with $\omega_{T,l'}=\kappa_{l'}h^\gamma$, for which $\partial\omega/\partial e = \gamma\eta\,\omega/e$ and $\partial\omega/\partial s=\gamma\phi\,\omega/s$ (so $\eta\to\gamma\eta$ in $(\star)$ and an extra $\gamma$ multiplies the $s$-FOC). $(\star)$ becomes $\sum_{l'}\pi_{l'\mid l}/C_{l'} = \gamma\eta/[(1-\gamma\eta)(1+\tau^e_{T,g})e]$, $\sum_{l'}\pi_{l'\mid l}\,g_{l'}/C_{l'} = 1/(1-\gamma\eta)$, and
$$
\frac{1}{1-s} = \frac{\mu\phi\gamma}{s(1-\gamma\eta)}
\quad\Longrightarrow\quad
\boxed{\;s_T^{*} \;=\; \frac{\mu\phi\gamma}{(\mu\phi-\eta)\gamma + 1}\;}.
$$
These coincide with the single-location, no-altruism constants of A.4, but here they hold for the full spatial model with any $L$, asymmetric $(t_{l'},\kappa_{l'},B_{l'},\tau_{l,l'})$, and any $\lambda$. The time share is thus constant across $(z,\vec\epsilon)$, across birth locations $l$, and across destinations, so the entire policy array $S$ collapses to the $I$ numbers $\{s_i^*\}$.

**Why $e$ stays implicit.** The cancellation that pins down $s$ used only $\sum_{l'}\pi_{l'\mid l}=1$ and the homogeneity of $g_{l'}$ in $(s,e)$; it never required solving $(\star)$ for $e$. But $(\star)$ itself is $\sum_{l'}\pi_{l'\mid l}/(G_{l'}e^\eta - (1+\tau^e)e) = \eta/[(1-\eta)(1+\tau^e)e]$ with $G_{l'}=(1-t_{l'})(1-\tau^\omega)A_i Q_l(z\epsilon_i)^\alpha (s^*)^\phi$ — a $\pi$-weighted sum of reciprocals with heterogeneous $G_{l'}$, which does not invert to a single power of $e$. It is closed-form only when the $C_{l'}$ are common across $l'$ (symmetric destinations, or degenerate $\pi$ as $\sigma_\nu\to0$, where $(\star)$ reduces to a single location's condition $C = (1-\eta)(1+\tau^e)e/\eta$). Hence $s$ is analytic but $e$ remains a scalar fixed point coupled to $\pi$. The altruistic sorting term $\Lambda_{l'}(z)$, the endogenous $\Phi_l$, the migration matrix, and the aggregates of §8–§9 are likewise unaffected by $\varepsilon$ and must still be solved numerically as in §10.

### A.7 Gumbel taste shocks on occupational choice

In the main text the occupation is chosen by a hard argmax (§7); the numerical implementation instead adds iid Gumbel taste shocks $\{\eta_i\}_{i=1}^{I}$ (scale $\theta_o>0$) to each occupation's value, paralleling the location-choice shocks of §6. With the occupational value
$$
W_{i,g,l}(z,\vec\epsilon) \;\equiv\; \ln\!\big(1-s^*_{i,g,l}\big) \;+\; \bar V\!\big(h^*_{i,g,l}(z,\epsilon_i);\,i,g,l,z\big),
$$
choosing $\arg\max_i\{W_{i,g,l}+\theta_o\eta_i\}$ yields the closed-form choice probabilities
$$
\boxed{\;\;
\rho_{i\mid g,l}(z,\vec\epsilon) \;=\; \frac{\exp\!\big[W_{i,g,l}(z,\vec\epsilon)/\theta_o\big]}{\sum_{i'=1}^{I}\exp\!\big[W_{i',g,l}(z,\vec\epsilon)/\theta_o\big]}.
\;\;}
$$
As $\theta_o\to 0$ these collapse to the hard Roy indicator $\mathbb 1_{i,g,l}(z,\vec\epsilon)$ of §8 (the code uses $\theta_o=0.05$); everywhere §8–§9 integrates against $\mathbb 1_{i,g,l}$, it is replaced by the smooth weight $\rho_{i\mid g,l}$.

**Why this matters for computation.** The household block solves the fixed point $V \mapsto h' \mapsto \Lambda \mapsto V$ by successive approximation. Under a hard argmax, $h'$ — and hence the altruism term $\Lambda$ and the value $V$ — is *piecewise-constant and discontinuous* in $V$ on the discretized $\vec\epsilon$ grid: any child cell near teacher/non-teacher indifference flips occupation every sweep, so $V$ hops between two branches in a period-2 limit cycle whose residual floor exceeds both the household and GE tolerances, exhausting the iteration budget and capping GE convergence. The Gumbel smoothing makes $h'$ continuous in $V$, restores the contraction, and reaches machine tolerance in a few sweeps — at a distortion shrinkable to zero with $\theta_o$. The occupational argmax was the lone remaining discontinuity in the map.

### A.8 Occupational thresholds: when are they closed form?

Write the occupational value with investments folded in,
$$
W_{i,g,l}(z,\vec\epsilon) \;\equiv\; \ln(1-s^*_{i,g,l}) \;+\; \bar V\!\big(h^*_{i,g,l}(z,\epsilon_i);\,i,g,l,z\big),
$$
so $i^*=\arg\max_i W_{i,g,l}$ and the threshold between occupations $i,j$ is the locus $W_{i,g,l}=W_{j,g,l}$ (the baseline teacher/other cutoff $a^*_{T,g}(a_O)$ of the draft is the $I=2$ instance).

**General spatial model ($\varepsilon>0$): no closed form.** $W_{i,g,l}$ inherits the coupled $(s,e,\pi)$ fixed point of §7 — both investment margins are state-dependent (A.6), $C$ averages over the child's $(z',\vec\epsilon',g')$, and $\bar V$ is a destination log-sum whose weights $\pi_{l'\mid l}$ depend on $h$. The threshold is only an implicit locus.

**$\varepsilon=0$ alone: still implicit.** Setting $\varepsilon=0$ collapses the time share to the constants $s^*_O,s^*_T$ of A.6 and decouples $C_{l'}$ from the child, so $W_{i,g,l}$ becomes a function of the single scalar $h_i$ (hence of $z\epsilon_i$). But the goods margin $e$ remains a $\pi$-weighted scalar fixed point and the log-sum $\bar V$ does not collapse, so $W_{i,g,l}$ is not an explicit function of $z\epsilon_i$ and the threshold stays implicit.

**$\varepsilon=0$ with destination-homogeneous net returns: closed form.** Add the condition that the active destinations deliver a common consumption level $C_{l'}\equiv C$ — single location $L=1$; or degenerate migration $\sigma_\nu\to0$; or symmetric net returns. Then the log-sum factors,
$$
\bar V(h_i;i,g,l,z) \;=\; \mu\ln C(h_i) \;+\; D(z,l),
\qquad
D(z,l)\equiv \sigma_\nu\ln\!\sum_{l'}\exp\!\big[(\Lambda_{l'}(z)+B_{l'}-\tau_{l,l'})/\sigma_\nu\big],
$$
and the additive term $D(z,l)$ is **occupation-independent** (the child's human capital does not depend on the parent's occupation), so it cancels in every pairwise comparison $W_i-W_j$. The two margins then solve in closed form (A.6 algebra with $C_{l'}=C$): for non-teaching $C^{(i)}=(1-\eta)(1-t)(1-\tau^\omega_{i,g})A_i\,h_i$ and for teaching $C^{(T)}=(1-\gamma\eta)(1-t)(1-\tau^\omega_{T,g})\kappa_l\,h_T^{\gamma}$, with
$$
h_i \;\propto\; (z\epsilon_i)^{\frac{\alpha}{1-\eta}}\quad(i\neq T),
\qquad
h_T \;\propto\; (z\epsilon_T)^{\frac{\alpha}{1-\gamma\eta}} .
$$
Hence $\ln C^{(i)}$ is affine in $\ln(z\epsilon_i)$:
$$
\ln C^{(i)} = \tfrac{\alpha}{1-\eta}\ln(z\epsilon_i)+\Theta_i\ (i\neq T),
\qquad
\ln C^{(T)} = \tfrac{\gamma\alpha}{1-\gamma\eta}\ln(z\epsilon_T)+\Theta_T,
$$
with intercepts
$$
\Theta_i=\frac{\ln[(1-\tau^\omega_{i,g})A_i]+\ln Q_l}{1-\eta}-\frac{\eta\ln(1+\tau^e_{i,g})}{1-\eta}+c_O,
\qquad
\Theta_T=\frac{\ln[(1-\tau^\omega_{T,g})\kappa_l]+\gamma\ln Q_l}{1-\gamma\eta}-\frac{\gamma\eta\ln(1+\tau^e_{T,g})}{1-\gamma\eta}+c_T,
$$
where $c_O,c_T$ collect pure parameter constants ($\eta,\gamma,t,s^*$). Two cases follow.

*Non-teaching vs. non-teaching — comparative advantage.* $s^*_O$ is common, so $W_i=W_j\iff \ln C^{(i)}=\ln C^{(j)}$, a cutoff in **relative** ability only:
$$
\boxed{\;\;
\ln\frac{\epsilon_i}{\epsilon_j}
\;=\;
\frac{1}{\alpha}\ln\frac{(1-\tau^\omega_{j,g})A_j}{(1-\tau^\omega_{i,g})A_i}
\;-\;
\frac{\eta}{\alpha}\ln\frac{1+\tau^e_{j,g}}{1+\tau^e_{i,g}}.
\;\;}
$$
It is independent of $z$, of the teacher aggregate $\widetilde H_{T,l}$ (the shifter $Q_l$ cancels), and of every spatial wedge $(t_{l'},B_{l'},\tau_{l,l'})$.

*Teaching vs. non-teaching — absolute advantage.* The ability exponents differ ($\tfrac{\gamma\alpha}{1-\gamma\eta}\neq\tfrac{\alpha}{1-\eta}$), so $z$ no longer cancels and $W_T=W_O$ solves explicitly for the cutoff $\epsilon_T^*$:
$$
\boxed{\;\;
\frac{\gamma\alpha}{1-\gamma\eta}\,\ln(z\epsilon_T^{*})
\;=\;
\frac{\alpha}{1-\eta}\,\ln(z\epsilon_O)
\;+\;
\frac{1}{\mu}\ln\frac{1-s^*_O}{1-s^*_T}
\;+\;
\Theta_O-\Theta_T.
\;\;}
$$
The cutoff $\epsilon_T^*(\epsilon_O,z)$ rises with the persistent level $z$ — the absolute-advantage channel — and, through $Q_l$ in $\Theta_O-\Theta_T$, depends on $\widetilde H_{T,l}$ except in the linear-wage knife-edge $\gamma=1$, where the $Q_l$ terms cancel and the teaching threshold too becomes independent of $\widetilde H_{T,l}$ (the spatial analogue of the baseline's $\beta=\sigma$ / Proposition 3).

Away from destination homogeneity these are the leading-order cutoffs; the exact spatial threshold adds the destination-mixing correction $\bar V(h_i)-\mu\ln C(h_i)-D(z,l)$, which is occupation-specific and must be computed numerically.



