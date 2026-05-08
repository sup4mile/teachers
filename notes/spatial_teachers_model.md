# Spatial Model with Altruism

This document extends the baseline occupational-choice / teacher-spillover economy to a setting with $L=2$ locations and a "move-to-opportunity" altruistic link between generations. Parents choose the location in which their children are born and educated, and finance an exogenous share $\varepsilon$ of their children's educational expenditures. Locations differ in teacher quality, fiscal policy, and amenities; agents face Gumbel taste shocks over locations and an additive utility moving cost.

## 1 Environment

There are $L=2$ locations indexed by $l \in \{1,2\}$. In each location the same set of $I$ occupations is available, with $i=1$ (or $i=T$) denoting teaching. Locations share common occupational productivities $\{A_i\}_{i=2}^{I}$, common labor-market distortions $\{\tau^{\omega}_{i,g}\}$, and common educational barriers $\{\tau^{e}_{i,g}\}$. Locations differ along:

- the (endogenous) stock of teacher human capital $\widetilde{H}_{T,l}$,
- the local labor-income tax $t_l$, which funds teachers under a *local* balanced-budget rule (one per location),
- an exogenous amenity value $B_l$,
- (possibly) a location-specific teaching wage profile $\omega_{T,l}(\cdot)$.

In addition, a young agent born in $l$ who later works in $l'$ pays a utility moving cost $\tau_{l,l'} \geq 0$, with $\tau_{l,l}=0$. Moving costs attenuate the dispersion induced by the Gumbel taste shocks.

The economy is populated by overlapping generations of finite-lived agents (young and old) with constant total measure $M$. Mass in location $l$ is $M_l$, with $M_1+M_2=M$ and equal cohort sizes $M_l/2$ (young) and $M_l/2$ (old) within each location in stationary equilibrium. The split $\{M_l\}$ is pinned down by location flows. Groups $g \in \{m, f\}$ denote gender. A parent's child is male or female with equal probability, independently of the parent's gender. Labor-market distortions $\tau^\omega_{i,g}$, educational barriers $\tau^e_{i,g}$, and investment policies $(s_{i,g,l}, e_{i,g,l})$ are gender-specific, so the coin-flip over the child's gender $g'$ is payoff-relevant from the parent's perspective. Because gender is an independent coin flip, the gender share of *young* agents is $\tfrac{1}{2}$ in every location. The gender composition of *old* (working-age) agents in each location is endogenous, however: men and women face different distortions and sort differently across locations, so the gender mix among workers in $l'$ is determined in equilibrium by the flow equations.

## 2 Timing

Each cohort's life cycle proceeds in three stages.

**Stage 1 — Young.** Born in location $l$ (chosen by the parent in Stage 2 of the previous cohort). Conditional on $l$:

1. Ability vector $\vec a$ is drawn from $F(\vec a \mid \vec a^{p})$, where $\vec a^{p}$ is the parent's ability.
2. The agent chooses time investment $s_{i,g,l}$ and goods investment $e_{i,g,l}$ in occupation-specific human capital, and selects an occupation $i$.
3. The parent (currently old) finances a share $\varepsilon$ of $(1+\tau^{e}_{i,g})\,e_{i,g,l}$. The remaining share $1-\varepsilon$ is paid by the agent out of her own future labor income.

**Stage 2 — Transition.** The agent reaches working age. She draws i.i.d. Gumbel taste shocks $\{\nu_{l'}\}_{l'=1}^{L}$ and chooses a location $l' \in \{1,2\}$ in which to (i) earn labor income, (ii) pay local taxes, and (iii) raise and educate her own child. She forms expectations over the child's gender $g' \in \{m,f\}$ (each equally likely) and the child's ability $\vec a'$ (drawn from $F(\vec a' \mid \vec a)$), and the resulting human capital $h'$.

**Stage 3 — Old.** Conditional on $l'$:

1. She earns occupation-specific labor income.
2. Her child's ability $\vec a'$ is realized from a log-AR(1) process around her own ability.
3. She pays the parental share $\varepsilon$ of the child's (now-realized) educational expenses.
4. She consumes the residual and receives a warm-glow payoff $\lambda\,f(h')$ from the child's accumulated human capital.

## 3 Ability Transmission

For each occupation $i$, abilities follow a log-AR(1) across generations:
$$
\ln a_i' = \rho \ln a_i + \xi_i', \qquad \xi_i' \sim \mathcal{N}(0,\sigma_{\xi}^{2}),\qquad \rho \in [0,1).
$$
Stacking across $i$ yields $F(\vec a' \mid \vec a)$ with stationary marginal $F(\vec a)$, identical across groups. As $\rho \to 0$ we recover the i.i.d. baseline.

## 4 Technologies

**Human capital.** A young agent of group $g$ in location $l$ with abilities $\vec a$ acquires occupation-specific human capital
$$
h_{i,g,l}(\vec a) = (h_{T,l})^\beta (a_i)^\alpha (s_{i,g,l})^{\phi}\,(e_{i,g,l})^{\eta}\,N(h_{T,l})^{-\sigma},
$$
with $h_{T,l}$ the human capital of a generic teacher working in $l$. Student-teacher matching is random within $l$, and the indifference condition $h^\beta N(h)^{-\sigma}=$ const. across teachers in $l$ combined with the local class-size resource constraint
$$
\frac{M_l}{2} = \sum_g \int_0^\infty N(h_{T,l,g})\,dF_{T,l,g}(h_{T,l,g})
$$
yields
$$
h_{T, l, g}^\beta\, N(h_{T,l,g})^{-\sigma} = \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^\sigma,
\qquad
\widetilde{H}_{T,l} \equiv \sum_g \int_0^\infty h^{\beta/\sigma}\,dF_{T,l,g}(h),
$$
where $F_{T,l,g}(\cdot)$ is the (endogenous) cdf of human capital among teachers *working in* $l$ from group $g$. Substituting back gives the reduced-form
$$
h_{i,g,l}(\vec a) = \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right)^\sigma a_i^\alpha\, s_{i,g,l}^{\phi}\, e_{i,g,l}^{\eta}.
$$

**Production.** Non-teaching output is $y_{i,g} = A_i\,h_{i,g,l}$, with $\omega_{i,l'}(h)=A_i h$ for $i \neq T$ — i.e., non-teaching wages depend on training location $l$ only through $h_{i,g,l}$, not on the work location $l'$. Teaching wages $\omega_{T,l'}(h)$ are taken as exogenous primitives in shape (strictly increasing, continuously differentiable), with location-specific tax rates $t_{l'}$ adjusting endogenously to satisfy the local balanced-budget condition. A teacher trained in $l$ with human capital $h$ who works in $l'$ earns $\omega_{T,l'}(h)$.

## 5 Preferences

Let $\lambda \geq 0$ index the strength of altruism. A young agent born in $l$, group $g$, with ability vector $\vec a$, who chooses occupation $i$ when young, values
$$
U \;=\; \ln(1-s_{i,g,l}) \;+\; \mathbb{E}\!\left[\, \mu \ln C \;+\; B_{l'} - \tau_{l,l'} \;+\; \sigma_{\nu}\nu_{l'} \;+\; \lambda\,f(h')\,\right],
$$
where the expectation is over: the Gumbel shocks $\{\nu_{l'}\}$ (which determine $l'$), the child's gender $g' \in \{m,f\}$ (each with probability $\tfrac{1}{2}$), and the child's ability draw $\vec a' \sim F(\vec a' \mid \vec a)$ (which determines $h'$, the child's optimal occupation $i'$, and investment $e'_{i',g',l'}$). The moving cost $\tau_{l,l'}\geq 0$ is paid in utils and attenuates Gumbel-driven moves; $\tau_{l,l}=0$.

Conditional on realized $l'$, $\vec a'$, and child's group $g'$, old-age consumption is:
$$
C(l',\vec a', g' \mid i,g,l) \;=\; (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l}(\vec a)\big)
\;-\;(1-\varepsilon)(1+\tau^{e}_{i,g})\,e_{i,g,l}
\;-\;\varepsilon\,(1+\tau^{e}_{i'(\vec a',l'),g'})\,e'_{i'(\vec a',l'),g',l'}(\vec a').
$$
The first netting term is the agent's $(1-\varepsilon)$ share of her own childhood education. The second is the parental $\varepsilon$ share of her child's (now-realized) educational outlay; the child's optimal $i'$ and $e'$ are functions of the realized $(\vec a', g')$ and the child's birth location $l'$ (with corresponding aggregates $\widetilde{H}_{T,l'}$, $M_{l'}$, $t_{l'}$). The parent's own labor-market distortion $(1-\tau^\omega_{i,g})$ and own educational barrier $(1+\tau^e_{i,g})$ retain $g$, the parent's gender; only the child's educational outlay carries the uncertain $g'$.

Setting $\varepsilon=0$, $\lambda=0$, $L=1$, and $\rho = 0$ recovers the baseline model.

## 6 Location Choice

Having chosen occupation $i$ when young, the agent observes $\{\nu_{l'}\}$ at Stage 2 and picks
$$
l^{*} \in \arg\max_{l' \in \{1,2\}}\; V_{l'}(h_{i,g,l}; i,g,l) + B_{l'} - \tau_{l,l'} + \sigma_\nu \nu_{l'},
$$
where $V_{l'}$ is the expected value of working in $l'$, integrating out both the child's gender and ability:
$$
V_{l'}(h; i,g,l) \;=\; \tfrac{1}{2}\sum_{g'\in\{m,f\}}\mathbb{E}_{\vec a' \mid \vec a}\!\left[\, \mu \ln C(l',\vec a',g' \mid i,g,l) \;+\; \lambda\,f\!\big(h'(\vec a',l',g')\big)\,\right].
$$
The moving cost $\tau_{l,l'}$ enters additively in the location-choice utility, parallel to the amenity $B_{l'}$. The Gumbel structure delivers closed-form choice probabilities,
$$
\pi_{l'\mid l}(h; i,g) \;=\; \frac{\exp\!\big[(V_{l'}(h; i,g,l) + B_{l'} - \tau_{l,l'})/\sigma_\nu\big]}{\sum_{l''}\exp\!\big[(V_{l''}(h; i,g,l) + B_{l''} - \tau_{l,l''})/\sigma_\nu\big]},
$$
and the log-sum value
$$
\bar V(h; i,g,l) \;=\; \sigma_\nu \ln\!\sum_{l'}\exp\!\big[(V_{l'}(h; i,g,l) + B_{l'} - \tau_{l,l'})/\sigma_\nu\big]
$$
folds the Stage-2 problem cleanly into Stage-1 occupational choice. As $\sigma_\nu \to 0$, $\pi_{l'\mid l}$ collapses to a deterministic indicator with $\tau_{l,l'}$ acting as a wedge in net values; as $\tau_{l,l'} \to \infty$ for $l'\neq l$, agents stay in their birth location regardless of $\nu$.

## 7 Occupation and Human Capital Choice

When young, the agent solves
$$
i^*(\vec a; g, l) \;=\; \arg\max_{i \in \{1,\dots,I\}}\; \Big[\ln(1 - s^*_{i,g,l}) + \bar V\!\big(h^*_{i,g,l}(\vec a); i, g, l\big)\Big],
$$
with $\big(s^*_{i,g,l}, e^*_{i,g,l}\big)$ given by the FOCs below, taking the location-choice probabilities $\pi_{l'\mid l}$ and the child's policy functions as given.

For non-teaching ($i \neq T$, with $\omega_{i,l'}(h)=A_i h$), the FOCs reduce to:

**FOC for $e_{i,g,l}$:**
$$
\sum_{l'} \pi_{l'\mid l}(h_{i,g,l};i,g)\;\tfrac{1}{2}\!\sum_{g'\in\{m,f\}}\mathbb E_{\vec a'\mid \vec a}\!\left[\frac{(1-t_{l'})(1-\tau^\omega_{i,g})A_i\eta\,h_{i,g,l}/e_{i,g,l} \;-\; (1-\varepsilon)(1+\tau^e_{i,g})}{C(l',\vec a',g' \mid i,g,l)}\right]\;=\;0.
$$

**FOC for $s_{i,g,l}$:**
$$
\frac{1}{1-s_{i,g,l}} \;=\; \mu\,\phi\,(1-\tau^\omega_{i,g})A_i\,\frac{h_{i,g,l}}{s_{i,g,l}}\,\sum_{l'} \pi_{l'\mid l}(h_{i,g,l};i,g)\;\tfrac{1}{2}\!\sum_{g'\in\{m,f\}}\mathbb E_{\vec a'\mid \vec a}\!\left[\frac{1-t_{l'}}{C(l',\vec a',g' \mid i,g,l)}\right].
$$

For teaching ($i=T$), the FOCs are analogous with $\partial \omega_{T,l'}(h)/\partial h$ replacing $A_i$ and $\omega_{T,l'}(h)$ replacing $A_i h$ in the consumption term — i.e., $A_i \eta h/e$ becomes $\eta\,(\partial \omega_{T,l'}/\partial h)\,h/e$ and the wage entering $C$ is $\omega_{T,l'}(h_{T,g,l})$. Because $\omega_{T,l'}$ depends on the *work* location $l'$, optimal $(s_{T,g,l}, e_{T,g,l})$ depend on the full vector $\{\pi_{l'\mid l}\}_{l'}$ rather than only own-location returns.

Note two things. First, the moving cost $\tau_{l,l'}$ enters the FOCs only through the choice probabilities $\pi_{l'\mid l}$: $\tau_{l,l'}$ shifts the weights $\pi_{l'\mid l}$ on each location's return. Second, since $\pi_{l'\mid l}$ depends on $V_{l'}$, which depends on $C$, which depends on $(s,e)$, the FOCs constitute a fixed-point problem in $(s, e, \pi)$ for each $(i,g,l,\vec a)$.

## 8 Distributions and Aggregation

Let $\mathbb 1_{i,g,l}(\vec a) \equiv \mathbb 1[i^*(\vec a; g, l) = i]$ index optimal occupational choice when young. The mass of cohort-$t$ young in $(g,l)$ choosing occupation $i$ and working in $l'$ is
$$
m^t_{i,g,l, l'} \;=\; \tfrac{M_{l,t}}{4}\, \int \mathbb 1_{i,g,l}(\vec a)\,\pi_{l'\mid l}(h_{i,g,l}(\vec a); i, g)\,dF(\vec a).
$$

**Aggregate teacher human capital** working in $l'$ at $t+1$ (used to teach the cohort born in $l'$ at $t+1$) is
$$
\widetilde{H}_{T,l',t+1} \;=\; \sum_g \sum_{l_0}\, \tfrac{M_{l_0,t}}{4}\!\int \mathbb 1_{T,g,l_0}(\vec a)\,\pi_{l'\mid l_0}\!\big(h_{T,g,l_0}(\vec a); T, g\big)\,h_{T,g,l_0}(\vec a)^{\beta/\sigma}\,dF(\vec a).
$$

This is equivalent to the above definition of $\widetilde H$, but using policy functions and the exogenous distribution of ability rather the endogenous distribution of human capital.

**Population flows.** The mass in $l'$ next period is
$$
M_{l',t+1} \;=\; \sum_{l_0}\, M_{l_0,t}\,\bar\pi_{l_0, l',t},
\qquad
\bar\pi_{l_0, l',t} \;=\; \tfrac{1}{2}\sum_g \sum_i \int \mathbb 1_{i,g,l_0}(\vec a)\,\pi_{l'\mid l_0}(h_{i,g,l_0};i,g)\,dF(\vec a).
$$
The matrix $[\bar\pi_{l_0, l',t}]$ has rows that sum to 1 for each $t$. The $\tfrac{1}{2}$ comes from that young agents in every location are always evenly split by gender. The gender composition of old agents working in $l'$ is endogenous.

**Teacher distribution.** The cdf of human capital among teachers working in $l'$ from group $g$, $F_{T,l',g}(\cdot)$, is the (population-normalized) push-forward of the joint distribution of $\big(\vec a, l_0\big)$ via $h_{T,g,l_0}(\vec a)$, weighted by $\mathbb 1_{T,g,l_0}\,\pi_{l'\mid l_0}$. Only the moment $\widetilde{H}_{T,l'}$ enters the rest of the model.

## 9 Government Budget

Each location runs a separate balanced budget: payroll taxes on labor income earned in $l'$ fund teacher wages paid in $l'$. With taxable income
$$
\mathcal{I}_{l'} \;=\; \sum_g \sum_{l_0}\, \tfrac{M_{l_0}}{4}\, \sum_{i} \int \mathbb 1_{i,g,l_0}(\vec a)\,\pi_{l'\mid l_0}(h_{i,g,l_0}; i, g)\,(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}\!\big(h_{i,g,l_0}(\vec a)\big)\,dF(\vec a)
$$
and teacher wage bill
$$
\mathcal{W}_{l'} \;=\; \sum_g \sum_{l_0}\, \tfrac{M_{l_0}}{4}\!\int \mathbb 1_{T,g,l_0}(\vec a)\,\pi_{l'\mid l_0}(h_{T,g,l_0}; T, g)\,\omega_{T,l'}\!\big(h_{T,g,l_0}(\vec a)\big)\,dF(\vec a),
$$
the local rate $t_{l'}$ solves
$$
t_{l'}\,\mathcal{I}_{l'} \;=\; \mathcal{W}_{l'},\qquad l' \in \{1,2\}.
$$
Both the wage bill and tax base sum across training locations $l_0$, since teachers and other workers can migrate.

## 10 Stationary Equilibrium

A stationary equilibrium is a tuple
$$
\Big\{\,t_l,\; \widetilde{H}_{T,l},\; M_l,\; F_{T,l,g},\; \big(s_{i,g,l},\, e_{i,g,l}\big)_{i,g},\; i^*(\,\cdot\,;g,l)\,\Big\}_{l\in\{1,2\},\,g}
$$
such that:

1. **(Optimization)** For each $(g,l,\vec a)$, $\big(s_{i,g,l}, e_{i,g,l}\big)$ satisfy the FOCs given $\big(\widetilde{H}_{T,1}, \widetilde{H}_{T,2}, M_1, M_2, t_1, t_2\big)$ and the location-choice probabilities $\pi_{l'\mid l}$. Occupational choice $i^*(\vec a;g,l) = \arg\max_i\big[\ln(1-s_{i,g,l}) + \bar V(h_{i,g,l};i,g,l)\big]$.

2. **(Distributions)** $F_{T,l,g}$ and the population masses $M_l$ are stationary fixed points of the flow equations. Young agents are split evenly by gender ($\tfrac{1}{2}$ each) in every location by construction.

3. **(Aggregates)** $\widetilde{H}_{T,l} = \sum_g \int h^{\beta/\sigma}\,dF_{T,l,g}(h)$ is consistent with the law of motion.

4. **(Class size)** $h_{T,l,g}^{\beta} N(h_{T,l,g})^{-\sigma} = (2\widetilde{H}_{T,l}/M_l)^{\sigma}$ for all teachers in $l$, with $\sum_g\int N\,dF_{T,l,g} = M_l/2$.

5. **(Budget)** Each $t_l$ satisfies the local balanced-budget condition.


## 11 Computational Notes

The model does not admit closed-form analogues of baseline Proposition 2 (that occupational choice does not depend on aggregate teacher human capital $\widetilde H$) because teaching wages depend on the work location $l'$ rather than the training location $l$, so the prospective teacher's expected return is a $\pi_{l'\mid l}$-weighted average across locations rather than a single-location expression; and altruism plus the $\varepsilon$-share of child education couples the parent's consumption to the child's optimal $e'_{i'(\vec a',l'),g,l'}(\vec a')$, which is itself the solution of a problem in the child's location $l'$ with aggregates $(\widetilde{H}_{T,l'}, M_{l'}, t_{l'})$.

Coding the model amounts to four nested fixed points:

- **Inner:** for each $(\vec a, g, l, i)$, jointly solve the FOCs and the choice probabilities for $(s_{i,g,l}, e_{i,g,l}, \pi_{l'\mid l})$ given aggregates and the child's policy.
- **Child consistency:** the child's policy $\big(s'_{i',g',l'}, e'_{i',g',l'}\big)$ entering $C(l',\vec a',g' \mid i,g,l)$ must coincide with the equilibrium policy of a gender-$g'$ agent in stationary equilibrium, for each $g' \in \{m,f\}$.
- **Aggregation:** integrate individual policies over $(\vec a, g, l)$ to compute $\widetilde{H}_{T,l}$, $M_l$, $\theta^g_l$.
- **Budget:** find $t_l$ satisfying local budget for each $l$.

A natural outer loop iterates on $\big(\widetilde{H}_{T,1},\widetilde{H}_{T,2},M_1,M_2,t_1,t_2\big)$ until all flow, aggregation, and budget equations balance.
