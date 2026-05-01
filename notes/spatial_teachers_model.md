# Spatial Model with Altruism

This document extends the baseline occupational-choice / teacher-spillover economy to a setting with $L=2$ locations and a "move-to-opportunity" altruistic link between generations. Parents choose the location in which their children are born and educated, and finance an exogenous share $\varepsilon$ of their children's educational expenditures.

## 1 Environment

There are $L=2$ locations indexed by $l \in \{1,2\}$. In each location the same set of $I$ occupations is available, with $i=1$ denoting teaching. Locations share common occupational productivities $\{A_i\}_{i=2}^{I}$, common labor-market distortions $\{\tau^{\omega}_{i,g}\}$, and common educational barriers $\{\tau^{e}_{i,g}\}$. Locations differ only along the following dimensions:

- the (endogenous) stock of teacher human capital $\widetilde{H}_{T,l}$,
- the local labor-income tax $t_l$, which funds teachers under a *local* balanced-budget rule (one per location),
- and an exogenous amenity value $B_l$.

The economy is populated by overlapping generations of finite-lived agents (young and old). A measure $M_l$ of young agents resides in location $l$, with $M_1+M_2=M$ pinned down by location flows below.

## 2 Timing

Each cohort's life cycle proceeds in three stages.

**Stage 1 — Young.** Born in location $l$ (chosen by the parent in Stage 2 of the previous cohort). Conditional on $l$:

1. Ability vector $\vec a$ is drawn from $F(\vec a \mid \vec a^{p})$, where $\vec a^{p}$ is the parent's ability.
2. The agent chooses time investment $s_{i,g,l}$ and goods investment $e_{i,g,l}$ in occupation-specific human capital, and selects an occupation $i$.
3. The parent (currently old) finances a share $\varepsilon$ of $(1+\tau^{e}_{i,g})\,e_{i,g,l}$. The remaining share $1-\varepsilon$ is paid by the agent out of her own future labor income.

**Stage 2 — Transition.** The agent reaches working age. She draws i.i.d. Gumbel taste shocks $\{\nu_{l'}\}_{l'=1}^{L}$ and chooses a location $l' \in \{1,2\}$ in which to (i) earn labor income, (ii) pay local taxes, and (iii) raise and educate her own child. She forms expectations over the child's ability $\vec a'$ and resulting human capital $h'$.

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
Stacking across $i$ yields $F(\vec a' \mid \vec a)$ with stationary marginal $F(\vec a)$, identical across groups. As $\rho \to 0$ we recover the i.i.d. baseline, as $\rho \to 1$ ability is dynastic.

## 4 Technologies

**Human capital.** As in (3)–(4) of the baseline,
$$
h'_{i,g} = (h_{T,l})^\beta (a_i)^\alpha (s_{i,g,l})^{\phi}\,(e_{i,g,l})^{\eta}{N(h_{T,l})^{-\sigma}}
$$
with $h_{T,l}$ the human capital of the teacher to which the student is matched in location $l$. Class size in $l$ satisfies the *local* resource constraint
$$
\frac{M_l}{2} = \int_0^\infty N(h_{T,l})\,dF_{T,l}(h_{T,l}),
$$

Where $F_{T, l}(h)$ is the location-specific (endogenous) cdf of teacher human capital. Combined with the indifference condition governing the allocation of teachers across classrooms,

$$ 
N(h') = \left(\frac{h'}{h}\right)^\frac{\beta}{\sigma} N(h),
$$

yields a location-specific class size distribution characterized by

$$
h_{T, l, g} ^ \beta N(h_{T,l, g}) ^ {-\sigma} = \left(\frac{2\widetilde{H}_{T,l}}{M_l}\right) ^ \sigma
$$

$$
\widetilde{H}_{T,l} \equiv \sum_g \int_0^\infty (h_{T,l, g})^{\beta/\sigma}\,dF_{T,l,g}(h_{T,l, g}).
$$

**Production.** Non-teaching output is $y_{i,g} = A_i\,h_{i,g}$, identical across locations. Teaching wages $\omega_{T,l}(h_{T,l})$ are taken as exogenous primitives (strictly increasing, continuously differentiable in $h_{T,l}$); they are pinned down in levels in equilibrium by the local balanced-budget condition.

## 5 Preferences

Let $\lambda \geq 0$ index the strength of altruism. A young agent born in $l$ with ability vector $\vec a$ values
$$
U \;=\; \ln(1-s_{i,g,l}) \;+\; \mathbb{E}\!\left[\, \mu \ln C \;+\; B_{l'} \;+\; \sigma_{\nu}\nu_{l'} \;+\; \lambda\,f(h')\,\right],
$$
where the expectation is over the Gumbel shocks $\{\nu_{l'}\}$ (which determine $l'$) and over the child's ability draw $\vec a'$ (which determines $h'$ and the child's optimal occupation $i'$ and investment $e'_{i',g,l'}$).

Conditional on the realized $l'$ and $\vec a'$, old-age consumption is:
$$
C(l',\vec a') \;=\; (1-t_{l'})(1-\tau^{\omega}_{i,g})\,\omega_{i,l'}(h_i)
\;-\;(1-\varepsilon)(1+\tau^{e}_{i,g})\,e_{i,g,l}
\;-\;\varepsilon\,(1+\tau^{e}_{i'(\vec a'),g})\,e'_{i'(\vec a'),g,l'}(\vec a').
$$
The first netting term is the agent's $(1-\varepsilon)$ share of her own childhood education. The second is the parental $\varepsilon$ share of her child's (now-realized) educational outlay; the child's optimal $i'$ and $e'$ are functions of the realized $\vec a'$.

Setting $\varepsilon=0$, $\lambda=0$, and $L=1$ recovers the previous baseline budget constraint.

## 6 Location Choice

Having chosen occupation $i$ when young, the agent observes $\{\nu_{l'}\}$ at Stage 2 and picks
$$
l^{*} \in \arg\max_{l' \in \{1,2\}}\; V_{l'}(h_i; i,g,l) + B_{l'} + \sigma_\nu \nu_{l'},
$$
where $V_{1}$ and $V_{2}$ are the deterministic values of working in locations 1 and 2 respectively, integrating out the child's ability:
$$
V_{l'}(h_i; i,g,l) \;=\; \mathbb{E}_{\vec a' \mid \vec a}\!\left[\, \mu \ln C(l',\vec a') \;+\; \lambda\,f\!\left(h'(\vec a',l')\right) \,\right],\qquad l' \in \{1,2\}.
$$
The Gumbel structure delivers closed-form choice probabilities,
$$
\pi_{l'}(h_i,i,g,l) \;=\; \frac{\exp\!\big[(V_{l'}(h_i; i,g,l) + B_{l'})/\sigma_\nu\big]}{\sum_{l''}\exp\!\big[(V_{l''}(h_i; i,g,l) + B_{l''})/\sigma_\nu\big]},
$$
and the value
$$
\bar V(h_i; i,g,l) \;=\; \sigma_\nu \ln\!\sum_{l'}\exp\!\big[(V_{l'}(h_i; i,g,l) + B_{l'})/\sigma_\nu\big]
$$
folds the Stage-2 problem cleanly into Stage-1 occupational choice, replacing the deterministic old-age.

## 7 Occupation and Human Capital Choice

For non-teaching ($i \neq T$, with $\omega_i(h_i)=A_i h_i$), the FOCs reduce to:

**FOC for $e_{i,g,l}$:**
$$
\sum_{l'} \pi_{l'}\, \mathbb E_{\vec a'\mid \vec a}\!\left[\frac{(1-t_{l'})(1-\tau^\omega_{i,g})A_i\eta\,h_i/e \;-\; (1-\varepsilon)(1+\tau^e_{i,g})}{C(l',\vec a')}\right]\;=\;0.
$$

**FOC for $s_{i,g,l}$:**
$$
\frac{1}{1-s} \;=\; \mu\,\phi\,(1-\tau^\omega_{i,g})A_i\,\frac{h_i}{s}\,\sum_{l'}\pi_{l'}\,\mathbb E_{\vec a'\mid \vec a}\!\left[\frac{1-t_{l'}}{C(l',\vec a')}\right].
$$
