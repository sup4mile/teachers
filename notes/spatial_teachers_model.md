# A Spatial Model of Occupational Choice with Teacher Spillovers and Intergenerational Links

## 1. Overview

This note extends the baseline Roy / Hsieh–Hurst–Jones–Klenow–style model of occupational choice with teacher human-capital spillovers (hereafter "the baseline") to a two-location, overlapping-generations economy with intergenerational altruism. Relative to the baseline, three new elements are introduced:

1. **Spatial choice.** Working-age agents choose where to live from $\ell \in \{1, 2\}$ subject to idiosyncratic Gumbel preference shocks and exogenous amenities $B_\ell$. Teacher human capital, wages, labor-market wedges, and educational barriers are all indexed by $\ell$.
2. **Intergenerational altruism.** Parents receive a "warm glow" $\lambda f(h')$ from their child's human capital and fund an exogenous share $\varepsilon \in [0,1]$ of the child's goods investment in education.
3. **Skill inheritance.** A child's ability vector is drawn from a log-AR(1) process centered at the parent's ability, with innovations independent across occupations.

The human-capital production technology, class-size determination, and broad Roy-style occupational choice structure are inherited unchanged from the baseline.

---

## 2. Environment

### 2.1 Demographics and geography

Time is discrete; the population is stationary. Each period a mass $M/2$ of new students is born. An individual lives for two periods — *young* (student) and *old* (worker/parent). The economy has two locations indexed by $\ell \in \{1, 2\}$. Occupations are indexed $i \in \{1, \dots, I\}$, with $i = T$ reserved for teaching. Final-goods occupations have linear technologies $y_{i,g,\ell} = A_{i,\ell}\, h_{i,g,\ell}$.

Agents are partitioned into groups $g \in \{1, \dots, G\}$. For concreteness in the spatial model we think of $g$ as gender ($G = 2$), consistent with the "coin flip" over the child's gender described below.

### 2.2 Skills

Each agent is endowed with a vector $\vec a = (a_1, \dots, a_I) \in \mathbb{R}_+^I$. **Skills are inherited across generations** according to independent log-AR(1) processes:

$$
\log a'_i \;=\; \rho_i \log a_i + \eta'_i, \qquad \eta'_i \sim \mathcal{N}(0,\, \sigma^2_{\eta,i}),
$$

independent across $i$. Let $F(\vec a' \mid \vec a)$ denote the implied conditional distribution. The unconditional distribution is stationary, and its marginal equals the baseline $F(\vec a)$.

### 2.3 Gender

The child's gender $g' \in \{m, f\}$ is drawn i.i.d. with $P(g' = m) = 1/2$, independent of ability. Parents integrate over $g'$ when forming period-2 expectations.

### 2.4 Timing

For a single individual, the lifetime timing is:

**Period 1 (young, in birth location $\ell_1$):**

1. $(g,\, \vec a,\, \ell_1)$ is observed. ($\ell_1$ is set by the parent's period-2 location choice.)
2. Chooses occupation $o$ and investments $(s_o, e_o)$.
3. Human capital $h$ realizes deterministically from $(a_o, s_o, e_o, h_{T,\ell_1})$.
4. Experiences disutility $\ln(1 - s_o)$.

**Period 2 (old, parent):**

1. Location shocks $\{\nu_\ell\}_{\ell \in \{1,2\}}$ drawn i.i.d. Type-I EV.
2. Chooses location $\ell_2$, trading off wages, amenities, wedge regime, and shocks.
3. Child's gender $g'$ and ability $\vec a'$ realize; $\vec a' \sim F(\cdot \mid \vec a)$.
4. Child solves their own period-1 problem, producing $e'_{o'}(\vec a', g')$ and $h'$.
5. Parent consumes the residual of after-tax labor income net of own past educational outlays and the share-$\varepsilon$ contribution to the child.

The critical timing assumption — that **$\ell_2$ is chosen before $\vec a'$ is known** — is what makes the parent's consumption stochastic ex ante and drives most of the new computational content of the model.

---

## 3. Preferences

Lifetime utility from the perspective of a young agent with state $(g, \vec a, \ell_1)$ is

$$
\mathcal{U} \;=\; \ln(1 - s_o) \;+\; \mu\,\mathbb{E}\!\left[\,B_{\ell_2} + \xi\,\nu_{\ell_2} \;+\; \mathbb{E}_{g', \vec a' \mid \vec a}\!\left[\,\ln c_{\ell_2} \;+\; \lambda f(h')\,\right]\,\right],
$$

where the outer $\mathbb{E}$ is over the Gumbel shocks (and therefore over the optimal $\ell_2$), and the inner $\mathbb{E}_{g', \vec a' \mid \vec a}$ is over the child's draw. The warm-glow function $f(\cdot)$ is non-decreasing (e.g. $f(h') = \ln h'$ preserves the global log structure and is the natural default); $\lambda \ge 0$ is the altruism weight; $\xi > 0$ scales the Gumbel noise.

---

## 4. Technologies

### 4.1 Human-capital production

A young agent of type $(g, \vec a)$ in location $\ell_1$ who chooses occupation $i$ produces

$$
h_{i,g,\ell_1} \;=\; \phi(h_{T,\ell_1}, a_i)\, s_{i,g}^{\phi}\, e_{i,g}^{\eta}\, N(h_{T,\ell_1})^{\sigma}, \qquad \phi(h_{T,\ell_1}, a_i) \;=\; h_{T,\ell_1}^{\beta}\, a_i^{\alpha}.
$$

This is the baseline production function with one change: the teacher input $h_{T,\ell_1}$ is **location-specific**. Children are taught by teachers who live in their birth location. Class size $N(\cdot)$ continues to be chosen so that students are indifferent across teachers *within* a location (the baseline equation (7) holds location-by-location):

$$
N(h_{T,\ell_1,g'}) = \left(\frac{h_{T,\ell_1,g'}}{h_{T,\ell_1,g}}\right)^{-\beta/\sigma} N(h_{T,\ell_1,g}).
$$

### 4.2 Wages

Final-goods occupations: $\omega_{i,g,\ell}(h) = A_{i,\ell}\, h$ for $i \ne T$.
Teaching: $\omega_{T,g,\ell}(h_T)$ is strictly increasing, continuously differentiable, and potentially nonlinear, exactly as in the baseline.

---

## 5. Budget Constraint

Let $\tau^{\omega}_{i,g,\ell}$ denote a location-specific labor-market wedge, $\tau^{e}_{i,g,\ell}$ a location-specific educational barrier, and $t$ the flat labor-income tax. The old-age budget constraint is

$$
\boxed{\;c \;=\; (1-t)\,(1 - \tau^{\omega}_{o,g,\ell_2})\, \omega_{o,g,\ell_2}(h)
\;-\; (1-\varepsilon)(1 + \tau^{e}_{o,g,\ell_1})\, e_{o}
\;-\; \varepsilon\,(1 + \tau^{e}_{o',g',\ell_2})\, e'_{o'}(\vec a', g')\;}
$$

Three remarks on this constraint:

- Education is financed out of period-2 income, as in the baseline; the only change is the split into an own share $(1 - \varepsilon)$ and a child share $\varepsilon$.
- Because $e'_{o'}$ is realized only after $\vec a'$ and $g'$ are drawn — which is *after* the location choice — **consumption is stochastic ex ante** from the parent's perspective.
- The own educational wedge is $\tau^{e}_{o,g,\ell_1}$ (incurred when young in the birth location), whereas the child's wedge is $\tau^{e}_{o',g',\ell_2}$ (incurred when the child is young in $\ell_2$, which is the child's birth location).

---

## 6. The Household's Recursive Problem

### 6.1 The old agent's problem

Fix the state $(g, \vec a, h, \ell_1, o)$ carried from period 1. Define the location-specific expected utility **net of the Gumbel shock** as

$$
\bar V_{\ell}(g, \vec a, h, \ell_1, o) \;\equiv\; B_\ell \;+\; \mathbb{E}_{g', \vec a' \mid \vec a}\!\left[\, \ln c_{\ell} \;+\; \lambda\, f\!\left(h'(g', \vec a', \ell)\right) \,\right],
$$

where $c_\ell$ is the budget constraint evaluated at $\ell_2 = \ell$. The location choice is

$$
\ell_2^{\star} \;=\; \arg\max_{\ell \in \{1,2\}}\, \left\{\,\bar V_\ell(g, \vec a, h, \ell_1, o) + \xi\, \nu_\ell\,\right\}.
$$

With Type-I EV shocks, the ex-ante value (before the Gumbel draw) and choice probabilities take the standard log-sum-exp / logit forms:

$$
V_2(g, \vec a, h, \ell_1, o) \;=\; \xi\, \ln\!\left[\sum_{\ell=1}^{2} \exp\!\left(\frac{\bar V_\ell}{\xi}\right)\right] + \xi\,\bar{\gamma},
$$

$$
\pi_\ell(g, \vec a, h, \ell_1, o) \;=\; \frac{\exp(\bar V_\ell / \xi)}{\sum_{\ell'} \exp(\bar V_{\ell'} / \xi)},
$$

where $\bar\gamma$ is the Euler–Mascheroni constant (irrelevant for choices).

### 6.2 The young agent's problem

Given state $(g, \vec a, \ell_1)$:

$$
V_1(g, \vec a, \ell_1) \;=\; \max_{o \in \{1, \dots, I\}}\; \max_{s_o,\, e_o}\; \left\{\, \ln(1 - s_o) \;+\; \mu\, V_2\!\left(g, \vec a, h(a_o, s_o, e_o, h_{T,\ell_1}),\, \ell_1,\, o\right)\,\right\}.
$$

The young agent does **not** choose $\ell_1$ — they inherit it — but they *anticipate* the period-2 location option through the log-sum-exp $V_2$.

---

## 7. Characterization of Decisions

### 7.1 Location choice

The location choice probability $\pi_\ell$ has the standard logit form. A few observations worth flagging:

- **$h$ matters for sorting.** Since $\omega_{o,g,\ell}(h)$ typically differs across $\ell$ (through $A_{o,\ell}$ or the teacher wage profile), agents with higher $h$ sort differentially across locations. A high-wedge location for one's group induces a departure.
- **Warm glow induces teacher-driven sorting.** Because $h_{T,\ell}$ differs across locations, expected $h'$ and hence $\mathbb{E}[f(h')]$ differs; altruistic parents tilt toward the location with better expected child outcomes, independent of own-wage considerations.
- **Risk over $\vec a'$ matters.** Because $\bar V_\ell$ contains $\mathbb{E}_{\vec a'}[\ln c_\ell]$ rather than $\ln \mathbb{E}_{\vec a'}[c_\ell]$, log-risk-aversion over the child's ability realization enters the location ranking. A location with higher variance of child-investment cost $\varepsilon(1+\tau^e_{o',g',\ell}) e'_{o'}$ is penalized relative to the risk-neutral ranking.
- **Moving costs.** A term $-\kappa\, \mathbf{1}\{\ell \ne \ell_1\}$ can be added to $\bar V_\ell$ additively; this preserves the logit structure and only alters the $\bar V_\ell$ differential.

### 7.2 A commitment problem between parent and child

The parent does not choose the child's $e'_{o'}$: the child optimizes on their own behalf. Because the parent absorbs share $\varepsilon$ of the cost while the child enjoys the full human-capital return, the child faces a distorted price $(1-\varepsilon)(1+\tau^e_{o',g',\ell_2})$ rather than the full social price $(1+\tau^e_{o',g',\ell_2})$. The child therefore over-invests in goods (from the parent's perspective) relative to what an altruistic parent would privately prefer, and there is a wedge between the parent-preferred and the realized $e'_{o'}$.

Two natural alternative modeling closures of this issue:

1. **Non-cooperative (current specification).** Child chooses $e'_{o'}$; parent absorbs the residual.
2. **Parental transfer as a bequest.** Parent chooses $\varepsilon e'_{o'}$ as a direct bequest after observing $\vec a'$, which the child then uses. Requires a richer parent problem but closes the commitment gap.

We proceed with specification (1). This is the source of one of the model's most interesting implications and is worth preserving; the alternative is noted for robustness.

### 7.3 Young agent's FOCs — time and goods investment

Fix occupation $o$ in birth location $\ell_1$. Using the log-sum-exp structure of $V_2$, envelope derivatives are weighted by the logit probabilities:

$$
\frac{\partial V_2}{\partial h} \;=\; \sum_\ell \pi_\ell\, \frac{\partial \bar V_\ell}{\partial h},
\qquad
\frac{\partial V_2}{\partial e_o} \;=\; \sum_\ell \pi_\ell\, \frac{\partial \bar V_\ell}{\partial e_o}.
$$

Since $h'$ is a function of $(\vec a', g', \ell)$ alone — it does **not** depend on the parent's $h$ or $e_o$ — the warm-glow term drops out of these derivatives. Computing the remaining pieces:

$$
\frac{\partial \bar V_\ell}{\partial h} \;=\; \mathbb{E}_{g', \vec a'}\!\left[\frac{(1-t)(1-\tau^{\omega}_{o,g,\ell})\, \omega'_{o,g,\ell}(h)}{c_\ell}\right],
$$

$$
\frac{\partial \bar V_\ell}{\partial e_o} \;=\; -\,(1-\varepsilon)(1+\tau^{e}_{o,g,\ell_1})\, \mathbb{E}_{g', \vec a'}\!\left[\frac{1}{c_\ell}\right].
$$

Given $\partial h / \partial s = \phi h / s$ and $\partial h / \partial e = \eta h / e$:

**FOC for time investment $s_o$:**

$$
\frac{1}{1 - s_o} \;=\; \mu\, \frac{\phi\, h}{s_o}\, \underbrace{\sum_\ell \pi_\ell\, \mathbb{E}_{g', \vec a'}\!\left[\frac{(1-t)(1-\tau^{\omega}_{o,g,\ell})\,\omega'_{o,g,\ell}(h)}{c_\ell}\right]}_{\displaystyle \equiv\, W(h; o, g, \ell_1)}.
$$

**FOC for goods investment $e_o$:**

$$
(1-\varepsilon)(1 + \tau^{e}_{o,g,\ell_1})\, \underbrace{\sum_\ell \pi_\ell\, \mathbb{E}_{g', \vec a'}\!\left[\frac{1}{c_\ell}\right]}_{\displaystyle \equiv\, C(h; o, g, \ell_1)}
\;=\;
\frac{\eta\, h}{e_o}\, W(h; o, g, \ell_1).
$$

**Time–goods ratio:**

$$
\frac{s_o}{1 - s_o} \;=\; \frac{\mu\,\phi}{\eta}\, \cdot \, \frac{(1-\varepsilon)(1+\tau^{e}_{o,g,\ell_1})\, e_o \cdot C(h)}{W(h)} \cdot \frac{1}{1}  \;=\; \mu\,\phi\, h \cdot W(h)
$$

so that

$$
s_o \;=\; \frac{\mu\,\phi\, h\, W}{1 + \mu\,\phi\, h\, W}, \qquad
e_o \;=\; \frac{\eta\, h\, W}{(1-\varepsilon)(1+\tau^{e}_{o,g,\ell_1})\, C}.
$$

These look like closed-form rules, but they are **implicit**: $h$ depends on $(s_o, e_o)$, and $W$ and $C$ depend on $c_\ell$, which itself depends on $(e_o, h)$. The pair $(s_o, e_o)$ must be solved as a fixed point given $\pi_\ell$ (which itself depends on $h$).

**Contrast with the baseline.** In the baseline, log utility plus linear-in-$h$ wages meant that the $1/c$ terms in the FOCs cancel exactly, yielding the closed forms in baseline equations (23)–(26). Two distinct features of the spatial extension break this:

1. The probability weighting by $\pi_\ell$ combines several different location-specific wedges $(\tau^{\omega}_{o,g,\ell})$ and productivities $A_{o,\ell}$. Even in the absence of uncertainty over $\vec a'$, the ratio $W / C$ is no longer a scalar linear function of $h$.
2. Uncertainty over $\vec a'$ generates dispersion in $c_\ell$ across $\vec a'$-draws. $\mathbb{E}[1/c_\ell]$ is not the reciprocal of $\mathbb{E}[c_\ell]$, so Jensen-type terms cannot be eliminated.

Consequently, the policy functions must be computed numerically at each state.

### 7.4 Occupational choice

The young agent chooses

$$
o^{\star}(g, \vec a, \ell_1) \;=\; \arg\max_{o}\, V_1^{\,o}(g, \vec a, \ell_1),
$$

where $V_1^o$ is the inner maximum conditional on occupation $o$. The occupational-choice boundary in $\vec a$-space is implicit and, unlike the baseline's equation (30), cannot be written without computing the location-integrated values for each candidate occupation. The three forces that shape it are:

- **Location-specific returns.** Different wedges and productivity scales enter through $\bar V_\ell$, generating heterogeneous relative payoffs by $\ell$.
- **Inheritance.** Higher own $a_i$ raises expected $a'_i$ via the log-AR(1), hence expected warm glow; the effect is occupation-specific and asymmetric if $\rho_i$ varies across $i$.
- **Option value of location.** Occupations that deliver higher $h$ are more "portable" — their holders receive higher inclusive value $V_2$ because the log-sum-exp over locations benefits from higher wages in *both* locations.

---

## 8. Aggregation

### 8.1 Teacher human capital by location

The aggregate state now contains two teacher stocks, one per location:

$$
\widetilde H_{T,\ell} \;=\; \sum_{g}\, \int h_{T,g}^{-\beta / \sigma}\, d F_{T, g, \ell}(h_{T,g}),
$$

where $F_{T, g, \ell}$ is the stationary distribution of teacher human capital for old agents of group $g$ who (i) chose occupation $T$ in period 1 and (ii) chose location $\ell$ in period 2. Both margins are endogenous.

### 8.2 Dynamics

The law of motion for $\widetilde H_{T,\ell}$ depends on three joint objects:

1. The distribution of $(\vec a, g)$ among period-1 students *born* in $\ell_1 = \ell$ — determined by last period's parents' location choices.
2. The occupational selection of these students (which $o^{\star}(\vec a, g, \ell)$ picks $T$).
3. The period-2 location choices of prospective teachers — i.e., the logit probabilities $\pi_\ell$ applied to the conditional distribution of teacher human capital.

In stationary equilibrium, the migration of new teachers from $\ell_1$ to $\ell_2$ is cross-coupled between locations: the teacher stock in $\ell$ today depends on choices in $\ell'$ yesterday.

### 8.3 Government / public finance — a modeling choice

The baseline features a single local government with balanced budget,

$$
t \cdot \text{(total taxable income)} \;=\; \text{(total teacher pay)}.
$$

Three natural closures for the two-location model:

- **(a) Two local budgets.** Each location's tax base funds its own teachers; $t$ becomes location-specific. This introduces a Tiebout-type force: high-income workers cluster, fund better schools, and perpetuate across-location inequality in $h_{T,\ell}$.
- **(b) One national budget.** A single federal $t$ transfers to teachers nationally. Shuts down school-finance Tiebout sorting, but other location differences remain.
- **(c) Block grant plus local supplement.** A convex combination of (a) and (b). Useful for matching U.S. school-finance institutions.

The choice affects quantitative predictions substantially but does not alter the form of the household problem.

---

## 9. Equilibrium

**Definition.** Given exogenous objects $\{B_\ell,\ A_{i,\ell},\ \tau^{\omega}_{i,g,\ell},\ \tau^{e}_{i,g,\ell}\}$ and the skill transition kernel $F(\cdot \mid \vec a)$, a stationary equilibrium consists of:

1. Investment and occupational policies $\{s_{o,g,\ell_1}(\vec a),\, e_{o,g,\ell_1}(\vec a),\, o^{\star}_{g,\ell_1}(\vec a)\}$ solving the young agent's problem.
2. Location-choice probabilities $\{\pi_\ell(g, \vec a, h, \ell_1, o)\}$ solving the old agent's problem.
3. Stationary distributions of teachers and workers by location, $F_{T, g, \ell}$ and $F_{O, g, \ell}$.
4. Teacher wage functions $\omega_{T, g, \ell}(\cdot)$ and aggregates $\widetilde H_{T, \ell}$ satisfying the laws of motion.
5. Government budget(s) balanced under the chosen public-finance closure.

Consistency requires:

- The period-1 distribution of students by $\ell_1$ corresponds to the child distribution induced by last period's parents' location choices.
- Teacher wages clear the teaching market location-by-location.
- The joint distribution of $(\vec a, h, \ell_1, o)$ among old agents, weighted by $\pi_\ell$, reproduces the current cross-sectional population by location.

---

## 10. Open Modeling Choices

The following are flagged for future decisions, as discussed at various points above:

| Choice | Options | Tradeoff |
|---|---|---|
| Moving costs | (i) None; (ii) $-\kappa \mathbf{1}\{\ell\ne\ell_1\}$; (iii) origin-destination costs $\kappa_{\ell_1 \to \ell}$; (iv) time-cost scaling with $h$ | Higher $\kappa$ dampens Tiebout sorting and generates path dependence. |
| Public finance | (a) local, (b) national, (c) block-grant | Determines strength of school-funding sorting force. |
| Parent vs. child choice of $e'_{o'}$ | Non-cooperative (child), or parent-bequest | The bequest version removes the commitment wedge but complicates the period-2 problem. |
| Warm-glow functional form | $f(h') = \ln h'$, $f(h') = h'$, CRRA | Affects sensitivity of altruism-driven sorting to dispersion in $h'$. |
| Gender asymmetries | Group-$g$-specific wedges $\tau^{\omega}_{i,g,\ell}$ for both parent and child | The coin flip on $g'$ already generates within-family heterogeneity. |

---

## 11. Computational Considerations

Relative to the baseline, several features substantially raise the computational cost. The baseline admits a closed-form solution for $(s, e)$ given the aggregates, with the occupational boundary (30) available in explicit form. The spatial extension loses both.

**Sources of added computational cost.**

1. *No closed-form investment rules.* Consumption is stochastic ex ante (via $e'(\vec a', g')$) and location-choice probabilities weight the wage expectations. Log-utility no longer delivers cancellation; investment policies are fixed points of the FOCs.
2. *Nested expectations.* Each policy evaluation requires integrating over $(g', \vec a')$ at the parent stage, over Gumbel shocks implicitly via log-sum-exp, and over $\vec a$ at the child stage for the aggregate distributions.
3. *Larger state space.* The old agent's state is $(g, \vec a, h, \ell_1, o)$ — roughly doubled vs. the baseline by $\ell_1$. The aggregate state includes $(\widetilde H_{T,1},\, \widetilde H_{T,2})$ rather than a scalar.
4. *Cross-location fixed point.* Teacher stocks in each location depend on migration decisions of new teachers, which depend on wages and amenities in each location. Stationarity requires a two-dimensional outer fixed point.
5. *Joint inner fixed point.* The child's $(s, e)$ depend on $\pi_\ell$, which depend on the implied $h$, which depends on $(s, e)$. An inner fixed point is required at each gridpoint.

**A practical algorithm.**

1. Discretize $\vec a$ (e.g., Gauss–Hermite on each $\log a_i$, leveraging independence across $i$).
2. Guess $(\widetilde H_{T,1}, \widetilde H_{T,2})$ and teacher wage profiles.
3. For each state $(g, \vec a, \ell_1)$ and each candidate occupation $o$:
    a. Inner fixed point over $(s_o, e_o, \pi_\ell)$ using the FOCs and the logit formula.
    b. Compute $V_1^o$.
4. Set $o^\star = \arg\max_o V_1^o$.
5. Aggregate implied teacher distributions and choice probabilities.
6. Update $(\widetilde H_{T,1}, \widetilde H_{T,2})$ and wage profiles; iterate to convergence.
7. Check the government budget(s); adjust $t$ (or location-specific $t_\ell$) accordingly.

Independence of skill innovations across occupations reduces the multidimensional integration at step 3 to products of one-dimensional quadrature rules, which mitigates the curse of dimensionality substantially.
