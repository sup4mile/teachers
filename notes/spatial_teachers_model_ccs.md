# Constant Class Size and Random Teacher Assignment

The [baseline model](spatial_teachers_model_eps0.md) assumes efficient class sizes: teachers with more human capital
teach larger classes, by just enough to equalize school quality within a location (Proposition 1). This implies far
more variation in class size within locations than we observe. Here classes are equal-sized and each student is assigned at random to one teacher, whose
human capital becomes part of her state. The draw acts as a shift of the student's occupational shocks, so the model
remains tractable (Proposition 5). In a benchmark, parents value a location by a power mean of its teachers' human
capital, so dispersion matters, with a sign set by the curvature of the warm glow (Proposition 6). The model nests the
baseline with equal teachers within locations; Section 6 adds noisy effectiveness and several teachers per student.

## 1 Environment

We keep the baseline environment and change only how classes are formed. A young agent of gender $g$ in location $l$
has ability $z$, which follows a log-AR(1) across generations with persistence $\rho_z$, and iid lognormal occupational
shocks $\vec\epsilon$. She chooses a schooling time share $s$, a goods investment $e$ and an occupation $i$, where
$i = T$ is teaching. Her human capital is
$$
h = h_T^{\beta}N_l^{-\sigma}(z\epsilon_i)^{\alpha}s^{\phi}e^{\eta} \equiv q\,(z\epsilon_i)^{\alpha}s^{\phi}e^{\eta}, \tag{1}
$$
where $h_T$ is her teacher's human capital, $N_l$ her class size and $q$ her school input, with $\beta > 0$ and
$0 < \sigma < 1$. A non-teacher earns $(1-\tau^\omega_{i,g})A_ih$, and a teacher in $l$ earns $\kappa_lh^{\gamma}$ with
$\gamma < 1$. When old, she chooses where to live, consumes $C$, and has a child who is schooled there. Parents do not
pay for their children's education, so the value of location $l'$ is $\mu\ln C_{l'} + \Lambda_{l'}(z)$, plus amenities
and taste shocks, where
$$
\Lambda_{l'}(z) = \lambda\,\mathbb E\big[f(h')\,\big|\,l', z\big], \qquad f(h') = \frac{(h'/\bar h)^{1-\psi} - 1}{1-\psi}, \tag{2}
$$
is the expected warm glow from the child, with $f = \log(h'/\bar h)$ at $\psi = 1$.

When all students in a location receive the same $q$, the baseline reduces the young agent's problem to two scalars:
the teaching shock $\epsilon_T$ and the best non-teaching income capacity
$X_O^* = \max_{i \ne T}(1-\tau^\omega_{i,g})A_i\epsilon_i^{\alpha}$ (Propositions 2′ and 3). Location choice depends on
ability only through fixed goods costs and, when $\psi \ne 1$, through $\Lambda$ (Proposition 4). *Destination
homogeneity* refers to the case in which every location an agent might choose gives her the same consumption.

**Schools.** Every class in $l$ has $N_l = S_l/n_{T,l}$ pupils, where $S_l = M_l/2$ is the mass of young schooled in
$l$ and $n_{T,l} = \sum_g\int dF_{T,l,g}$ the mass of teachers who work there, $F_{T,l,g}$ being the baseline's measure
of their human capital. Let $G_l = \sum_gF_{T,l,g}/n_{T,l}$ be that distribution scaled to unit mass. Each student is
assigned at random to one teacher for all of school, so $h_T \sim G_l$, independent of $(g, z, \vec\epsilon)$, and
$q = h_T^{\beta}N_l^{-\sigma}$. Write $q = \bar Q_l\,r$, where
$$
\log\bar Q_l = \beta\bar\ell_l - \sigma\log N_l, \qquad \bar\ell_l = \mathbb E[\log h_T\mid l], \tag{3}
$$
is the geometric-mean school quality of $l$. The draw $\log r = \beta(\log h_T - \bar\ell_l)$ has mean zero, and its
law $R_l$ is $G_l$ in logs, centred and scaled by $\beta$:
$$
\kappa_m(\log r\mid l) = \beta^{m}\,\kappa_m(\log h_T\mid l), \qquad m \ge 2, \tag{4}
$$
so $\mathrm{Var}(\log r\mid l) = \beta^2\,\mathrm{Var}(\log h_T\mid l)$. The young agent observes $q$, together with
$(z, \vec\epsilon)$, before choosing $(s, e, i)$.

## 2 The Young Agent's Problem

**Proposition 5 (Effective shocks).** *The problem of a young agent with draw $r$ is the baseline problem at school
quality $\bar Q_l$ and effective shocks $\tilde\epsilon_i = r^{1/\alpha}\epsilon_i$. Hence (a) Propositions 2′ and 3
hold at every draw, with non-teaching statistic $rX_O^*$, and (b) under destination homogeneity and with no fixed goods
costs, the teaching decision depends on $(q, z)$ only through $qz^{\alpha}$, and a better draw lowers the probability
of teaching.*

The draw enters the problem only through (1), where $q(z\epsilon_i)^{\alpha} = \bar Q_l(z\tilde\epsilon_i)^{\alpha}$
(Appendix A.1). Two implications follow. First, the draw adds no dimension to the policy functions: every aggregate is
a mixture, over draws, of baseline objects under shifted shock laws. Second, the draw acts like ability. It shifts
absolute advantage as $z^{\alpha}$ does, so students who draw better teachers select out of teaching, and the agent's
own draw adds no sorting force beyond that of ability: it enters location choice only through fixed goods costs
(Proposition 4). The child's draw does add one, through $\Lambda$ (Section 3).

## 3 Parents' Valuation of Teacher Quality

A child draws her teacher where her parent lives, so the expectation in (2) now also integrates over the child's draw
$r' \sim R_{l'}$. Parents value the whole law $R_{l'}$, not only $\bar Q_{l'}$. A benchmark shows how.

**Benchmark.** Take the non-teaching branch, setting aside the child's occupational choice (which depends on $r'$ by
Proposition 5b), with no fixed goods costs and destination homogeneity. Then $h' \propto q'^{1/(1-\eta)}u$, where
$u$, the child's ability and shocks, is independent of $q'$. Let $\mathcal M_p(x\mid l) = (\mathbb E[x^p\mid l])^{1/p}$
be the power mean of order $p$ over the teachers of $l$, with $\mathcal M_0$ the geometric mean.

**Proposition 6 (Valuation of teacher quality).** *In the benchmark, let $\chi = (1-\psi)/(1-\eta)$. Then, for a
constant $c > 0$,*
$$
\mathbb E\big[f(h')\,\big|\,l', z\big] = \frac{1}{1-\psi}\Big[c\,N_{l'}^{-\sigma\chi}\,\mathcal M_{\chi\beta}(h_T\mid l')^{\chi\beta}\,\mathbb E\big[u^{1-\psi}\mid z\big] - 1\Big]. \tag{5}
$$
*Parents value a location by the power mean of order $\chi\beta$ of its teachers' human capital. Holding
$\bar Q_{l'}$ fixed, dispersion in teacher human capital raises the location's value if $\psi < 1$, leaves it unchanged
if $\psi = 1$, and lowers it if $\psi > 1$.*

With $\bar Q_{l'}$ fixed, a spread in $\log q'$ is a spread in $\log h'$, amplified by $1/(1-\eta)$ because goods
investment rises with the school input. The kernel values it if it is convex in $\log h'$ ($\psi < 1$) and penalizes it
if it is concave ($\psi > 1$); the log kernel is linear in $\log h'$. The effect is second order. If $h_T$ is lognormal
within $l'$, the certainty equivalent of the draw, in log human capital, is
$$
\frac{(1-\psi)\,\beta^2}{2(1-\eta)^2}\,\mathrm{Var}(\log h_T\mid l'). \tag{6}
$$

**Corollary (Dispersion and sorting).** *In the benchmark, let $m_\psi = (1-\psi)\alpha/(1-\eta)$, $\rho_z > 0$, and
$\Lambda_2(z) - \Lambda_1(z) = \Gamma(z) + D(z)$, where $\Gamma$ is the gap with $r \equiv 1$ in both locations and $D$
the contribution of dispersion. Then $\partial X/\partial\log z = m_\psi\rho_zX$ for $X \in \{\Gamma, D\}$, and:
(a) if $\bar Q_2 > \bar Q_1$, $\Gamma$ rises with $z$ if $\psi < 1$, is constant if $\psi = 1$, and falls with $z$ if
$\psi > 1$, as in the baseline; (b) if $\bar Q_2 = \bar Q_1$ and $\log r_2$ is a mean-preserving spread of $\log r_1$,
$D$ rises with $z$ for every $\psi \ne 1$.*

Dispersion thus favors high-ability parents whatever the curvature. For $\psi < 1$ it raises a location's value by
more for them; for $\psi > 1$ it lowers the value by less.

## 4 Relation to Efficient Class Sizes

With equal teachers in each location, $r \equiv 1$ and the two models coincide. Otherwise, fix a location's teachers
and its pupil–teacher ratio. Efficient class sizes give every student the baseline's $Q_l$, written $Q_l^{\text{eff}}$;
random assignment gives a mean $\mathbb E[q\mid l]$ and a geometric mean $\bar Q_l$, all power means of $h_T$:
$$
Q_l^{\text{eff}} = N_l^{-\sigma}\mathcal M_{\beta/\sigma}(h_T\mid l)^{\beta}, \qquad \mathbb E[q\mid l] = N_l^{-\sigma}\mathcal M_{\beta}(h_T\mid l)^{\beta}, \qquad \bar Q_l = N_l^{-\sigma}\mathcal M_0(h_T\mid l)^{\beta}. \tag{7}
$$
The first is the baseline's $Q_l = (2\widetilde H_{T,l}/M_l)^{\sigma}$, since
$\widetilde H_{T,l} = n_{T,l}\,\mathbb E[h_T^{\beta/\sigma}\mid l]$ and $M_l/2 = S_l$. Power means increase in their
order and $\sigma < 1$, so $Q_l^{\text{eff}} \ge \mathbb E[q\mid l] \ge \bar Q_l$. The gain from matching class size to
teacher human capital, in expected log school input, is $\log Q_l^{\text{eff}} - \log\bar Q_l$; to second order in the
dispersion of $\log h_T$ it is $\beta^2\,\mathrm{Var}(\log h_T\mid l)/(2\sigma)$, exactly so if $h_T$ is lognormal.

## 5 Equilibrium

**Definition.** A stationary equilibrium consists of local tax rates $t_l$, teacher masses and human-capital
distributions $(n_{T,l}, G_l)$ (in place of the baseline's $\widetilde H_{T,l}$), populations $M_l$, distributions
$\Phi_l(z)$ of the young by ability, and policies such that:

1. given $(\bar Q_l, R_l, t_l)$, policies solve the problems of Proposition 5, and $\Lambda$ satisfies (2) given the
   children's policies and draws;
2. $\Phi$ is stationary under the migration kernel, integrated over occupations, shocks and draws;
3. $(n_{T,l}, G_l)$ are generated by teachers' choices;
4. $N_l = S_l/n_{T,l}$, and $\bar Q_l$ and $R_l$ follow from $G_l$ by (3) and (4);
5. local taxes fund local teacher pay, $t_l\mathcal I_l = \mathcal W_l$, where $\mathcal I_l$ is taxable income earned
   in $l$ and $\mathcal W_l$ is the teacher wage bill.

## 6 Extension: Noisy Effectiveness and Several Teachers

The model generalizes in two directions. Teacher $j$'s effectiveness is $\log\theta_j = \log h_{T,j} + \nu_j$, with
$\nu_j$ independent mean-zero noise scaled so that $\mathrm{corr}(\log\theta, \log h_T\mid l) = \rho$, and a student
draws $K_d$ teachers independently of $(g, z, \vec\epsilon)$ and of each other (a new one each year if $K_d$ is the
number of grades), all known before she chooses $(s, e, i)$, with school input
$$
q = \prod_{t=1}^{K_d}\big(\theta_{j_t}^{\beta}N_l^{-\sigma}\big)^{1/K_d}. \tag{8}
$$
Sections 1–5 are the case $\theta = h_T$ ($\rho = 1$) and $K_d = 1$. Since $\mathbb E[\log\theta\mid l] = \bar\ell_l$,
(3) is unchanged, while (4) becomes
$$
\kappa_m(\log r\mid l) = \frac{\beta^{m}}{K_d^{\,m-1}}\,\kappa_m(\log\theta\mid l), \quad m \ge 2, \qquad \mathrm{Var}(\log\theta\mid l) = \frac{\mathrm{Var}(\log h_T\mid l)}{\rho^2}: \tag{9}
$$
noise adds luck, and more teachers average it out. Proposition 5 and the Corollary use only $q = \bar Q_lr$ with $r$
independent of $(g, z, \vec\epsilon)$ and $\mathbb E\log r = 0$, so they carry over unchanged. In (5),
$\mathcal M_{\chi\beta}(h_T\mid l')$ becomes $\mathcal M_{\chi\beta/K_d}(\theta\mid l')$, so the weight on dispersion
falls with $K_d$, and (6) is divided by $K_d$, with $\theta$ for $h_T$. In (7), with efficient class sizes matched to
effectiveness, $\theta$ replaces $h_T$ and $\mathbb E[q\mid l] = N_l^{-\sigma}\mathcal M_{\beta/K_d}(\theta\mid l)^\beta$;
since $\sigma < 1 \le K_d$, the ordering is unchanged.

**Table 1.** Efficient class sizes, constant class size (Sections 1–5) and its extension (Section 6)

| | Efficient class sizes | Constant class size | Extension |
|---|---|---|---|
| Class size in $l$ | $N(h_T) \propto h_T^{\beta/\sigma}$ | $N_l$ | $N_l$ |
| Teachers of a student | irrelevant | one, $h_T \sim G_l$ | $K_d$, effectiveness $\theta$ |
| School input in $l$ | $Q_l^{\text{eff}}$, common to all | $\bar Q_lr$, a random draw | $\bar Q_lr$, a random draw |
| $\mathrm{Var}(\log r\mid l)$ | $0$ | $\beta^2\mathrm{Var}(\log h_T\mid l)$ | $\beta^2\mathrm{Var}(\log h_T\mid l)/(\rho^2K_d)$ |
| State of the young | $(g, l, z, \vec\epsilon)$ | $(g, l, z, \vec\epsilon, q)$ | $(g, l, z, \vec\epsilon, q)$ |
| Teacher moment parents value (benchmark) | $\mathcal M_{\beta/\sigma}(h_T\mid l)$ | $\mathcal M_{\chi\beta}(h_T\mid l)$ | $\mathcal M_{\chi\beta/K_d}(\theta\mid l)$ |

## Appendix A Proofs

### A.1 Proposition 5

By (1), $h = \bar Q_l(z\tilde\epsilon_i)^{\alpha}s^{\phi}e^{\eta}$, and nothing else in the problem depends on $q$
($\Lambda$ depends on the child's draw, not the parent's), which is the main claim. (a) Non-teaching income is
$(1-\tau^\omega_{i,g})A_ih = \bar Q_lz^{\alpha}(rX_{O,i})s^{\phi}e^{\eta}$, so $rX_{O,i}$ plays the role of $X_{O,i}$.
(b) With destination homogeneity and no fixed goods costs, the branch values are affine in $\log(qz^{\alpha})$, up to
common terms, with coefficient $\mu\gamma/(1-\gamma\eta)$ in teaching and $\mu/(1-\eta)$ outside it (baseline Appendix
A.3); the first is smaller iff $\gamma < 1$. Fixed goods costs give the slopes branch-specific factors and can make the
time share depend on $z$ through the location probabilities (Proposition 2′). $\square$

### A.2 Proposition 6

With one draw from $G_l$, $\mathbb E[q^{\chi}] = N^{-\sigma\chi}\,\mathbb E[h_T^{\chi\beta}]
= N^{-\sigma\chi}\,\mathcal M_{\chi\beta}(h_T)^{\chi\beta}$. Since $(h'/\bar h)^{1-\psi} \propto q'^{\chi}u^{1-\psi}$
with $q'$ independent of $u$, this gives (5). Write $\mathbb E[q^{\chi}] = \bar Q^{\chi}\,\mathbb E[r^{\chi}]$. Because
$\mathbb E\log r = 0$ and $e^{\chi x}$ is convex, $\mathbb E[r^{\chi}] \ge 1$, and it rises under a mean-preserving
spread of $\log r$; the factor $1/(1-\psi)$ in (5) gives the sign.

*Equation (6) and the efficiency gain.* In the cumulants of $\log x$,
$\log\mathcal M_p(x) = \kappa_1 + \tfrac12p\kappa_2 + O(\kappa_3)$, with no remainder if $x$ is lognormal. The
certainty equivalent is the factor $c_r$ on $h'$ with $c_r^{1-\psi} = \mathbb E[r^{\chi}]$, so
$\log c_r = \log\mathcal M_{\chi}(r)/(1-\eta)$; with (4), this gives (6). Applied to $h_T$ with $p = \beta/\sigma$ and
$p = 0$ in (7), the same expansion gives the efficiency gain.

### A.3 Corollary

Write $\mathbb E[q_l^{\chi}] = \bar Q_l^{\chi}(1 + d_l)$ with $d_l = \mathbb E[r_l^{\chi}] - 1 \ge 0$, and let
$\upsilon(z) = \mathbb E[u^{1-\psi}\mid z]$. By (5), with $k = \lambda c > 0$,
$$
\Gamma = \frac{k\,(\bar Q_2^{\chi} - \bar Q_1^{\chi})}{1-\psi}\,\upsilon(z), \qquad D = \frac{k\,(\bar Q_2^{\chi}d_2 - \bar Q_1^{\chi}d_1)}{1-\psi}\,\upsilon(z).
$$
The log-AR(1) gives $\upsilon(z) \propto z^{m_\psi\rho_z}$ (baseline Appendix A.4), so
$\partial X/\partial\log z = m_\psi\rho_zX$. Since $m_\psi/(1-\psi) = \alpha/(1-\eta) > 0$ and $\rho_z > 0$,
$\partial\Gamma/\partial\log z$ has the sign of $\bar Q_2^{\chi} - \bar Q_1^{\chi}$, which for $\bar Q_2 > \bar Q_1$ is
the sign of $1-\psi$. At $\bar Q_1 = \bar Q_2$, $\partial D/\partial\log z$ has the sign of $d_2 - d_1$, positive for
every $\chi \ne 0$ because $e^{\chi x}$ is convex and $\log r_2$ is a mean-preserving spread of $\log r_1$. $\square$

### A.4 Section 6

The proofs go through with $h_T$ replaced by $\theta$ and one draw by $K_d$ independent ones. In (9),
$\log r = (\beta/K_d)\sum_t(\log\theta_t - \bar\ell_l)$ is a sum of independent terms, and cumulants of order
$m \ge 2$ are additive over independent terms and homogeneous of degree $m$, so
$\kappa_m(\log r) = K_d(\beta/K_d)^m\kappa_m(\log\theta)$; the variance follows from
$\mathrm{Cov}(\log\theta, \log h_T\mid l) = \mathrm{Var}(\log h_T\mid l)$. In (5) and (7),
$\mathbb E[q^{\chi}] = N^{-\sigma\chi}(\mathbb E[\theta^{\chi\beta/K_d}])^{K_d}
= N^{-\sigma\chi}\mathcal M_{\chi\beta/K_d}(\theta)^{\chi\beta}$, with $\chi = 1$ for $\mathbb E[q\mid l]$. $\square$
