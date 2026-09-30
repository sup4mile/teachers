# =============================================================================
# Model counterparts of the calibration targets in
# julia/spatial_model/calibration/spatial_calibration.md (§2–3; to-do item T4).
#
# Every function reads a converged `sol = solve_ge(p)` and integrates with the
# solver's own choice maps, so masses and budgets agree with `aggregates`.
#
# Conventions
#   - Gender g = 1 is male and g = 2 female (the female wage wedge is τω[:, 2]).
#   - An adult's location l′ is where she works, lives, pays tax and schools her
#     child: the model bundles residence, work and school.
#   - Students in l′ are the young born there, S_{l′} = Σ_z Φ[z, l′] = M_{l′}/2.
#   - Spatial gaps are location 2 (suburb) minus location 1 (city), in logs.
#
# Layout
#   1  Measurement assumptions
#   2  Work-location aggregates
#   3  CFR fixed-class-size teacher intervention
#   4  BFM willingness to pay for school quality
#   5  Wage distributions: a weighted pseudo-sample of adults
#   6  All moments
#   7  External-block consistency (Table 2 in the full equilibrium, §4)
# =============================================================================

isdefined(@__MODULE__, :solve_ge) || include(joinpath(@__DIR__, "..", "spatial_continuous.jl"))

# -----------------------------------------------------------------------------
# 1. Measurement assumptions
# -----------------------------------------------------------------------------
"""
    Measurement(; χ, K, teacher_sd, wtp_score_sd, income)

Assumptions that map model objects into the units of the data. None of these is
estimated; each is a documented choice (§1 and §3 of the calibration note).
- `χ = log(1.12)`: log adult earnings per student SD of achievement (the CFR bridge).
- `K = 12`: school years of exposure (grades 1–12), used for the teacher
  intervention and SEDA.
- `teacher_sd = :within`: SD of teacher log h within work locations (`:pooled`
  uses the economy-wide SD).
- `wtp_score_sd = 0.45`: BFM's unit, one SD of school-average scores, in student
  SDs: √ICC with a grade-4 between-school variance share near 0.2 (Hedges–Hedberg:
  0.23–0.24 nationally, 0.17–0.18 on state averages; range 0.41–0.49). BFM report
  no student-level SD for their CLAS scores (data/spatial/estimates/literature_targets.md).
- `income = :mean_log`: residence income statistic. y is pre-tax labour income,
  home-production output included, over all adults: `:mean_log` is E[log y], `:log_mean`
  log E[y], `:log_median` the log median. The `_market` variants (`:mean_log_market`,
  `:log_median_market`) restrict to market workers, teachers included and home
  producers (occupation `hp`) excluded, the population of an earnings median.
- `hp = 0`: index of the home-production occupation (0: none). Needed by the `_market`
  income statistics and the wage-sample moments, which exclude home producers.
- `teacher_weights = (0.2444, 0.7556)`: weights of men's and women's teacher 90/10s
  in the pooled teacher 90/10, as in the data (ACS 2009–13 K–12 wage-sample counts,
  12,118 men and 37,465 women).
- `pov_rate = 0.18`: economy-wide share of children counted as poor; sets the income
  line of the child-poverty validation moment (provisional, T4).
- `nquad = 600`: midpoints per shock dimension in the pseudo-sample (`micro_sample`).
- `seda_grade = 5.5`, `child_grade = 4.0`: mean grade at which SEDA scores (grades
  3–8) and CNLSY child scores (ages 5–14, averaged over rounds) are taken. A score
  at grade g carries g/K of the full-exposure loading c on log Q_l; parents' scores
  (NLSY79 AFQT at ages 15–23) carry all of it (`ability_moments`).
"""
Base.@kwdef struct Measurement
    χ::Float64            = log(1.12)
    K::Float64            = 12.0
    teacher_sd::Symbol    = :within
    wtp_score_sd::Float64 = 0.45
    income::Symbol        = :mean_log
    hp::Int               = 0
    teacher_weights::NTuple{2,Float64} = (12_118.0, 37_465.0) ./ 49_583.0
    pov_rate::Float64     = 0.18
    nquad::Int            = 600
    seda_grade::Float64   = 5.5
    child_grade::Float64  = 4.0
end

"SEDA growth gap implied by a quality gap (§3): Δg ≈ Δlog Q / (χ K (1 − η))."
seda_growth(ΔlogQ, χ, K, η) = ΔlogQ / (χ * K * (1 - η))

# -----------------------------------------------------------------------------
# 2. Work-location aggregates
# -----------------------------------------------------------------------------
"""
    location_stats(sol)

One pass over (birth location, ability, gender) cells, integrating over the
occupation choice and the work-location logit. Indexed by work location l′:
- `nT[g, l′]`, `nO[g, l′]`: masses of teachers and non-teachers;
- `pay[l′]`: teacher pay Σ (1−τω_{T,g}) κ_{l′} h_T^γ;
- `logh[l′]`, `logh2[l′]`: Σ log h_T and Σ (log h_T)² over teachers;
- `logy[l′]`, `y[l′]`: Σ log y and Σ y over all workers, y = pre-tax labour income;
- `S[l′]`: students, the young born in l′.
"""
function location_stats(sol)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L) = gr
    nT = zeros(2, L); nO = zeros(2, L)
    pay = zeros(L); logh = zeros(L); logh2 = zeros(L); logy = zeros(L); y = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm   = choice_maps(hh, g, l, zi, gr)
        lhT  = log.(hh.hT[g, l, zi, :])
        yO   = hh.capO[g, l, zi, :]
        lyO  = log.(yO)
        for lp in 1:L
            πT  = Pi_T[g, l, zi, :, lp]
            πO  = Pi_O[g, l, zi, :, lp]
            lyT = log((1 - p.τω[p.T, g]) * p.κ[lp]) .+ p.γ .* lhT
            payT = mass * integrate_teach(πT .* exp.(lyT), cm, gr)
            nT[g, lp] += mass * integrate_teach(πT, cm, gr)
            nO[g, lp] += mass * integrate_nonteach(πO, cm, gr)
            pay[lp]   += payT
            logh[lp]  += mass * integrate_teach(πT .* lhT, cm, gr)
            logh2[lp] += mass * integrate_teach(πT .* lhT .^ 2, cm, gr)
            logy[lp]  += mass * (integrate_teach(πT .* lyT, cm, gr) +
                                 integrate_nonteach(πO .* lyO, cm, gr))
            y[lp]     += payT + mass * integrate_nonteach(πO .* yO, cm, gr)
        end
    end
    return (; nT, nO, pay, logh, logh2, logy, y, S = vec(sum(Φ; dims = 1)))
end

"P(i* = i | X_O* = x): the chance that non-teaching occupation `i` attains the max."
function argmax_prob(gr, g, x, i)
    w(j) = pdf(gr.dO[j, g], x) / max(cdf(gr.dO[j, g], x), 1e-300)
    den = sum(w, gr.nonteach)
    return den > 0.0 ? w(i) / den : 0.0
end

"""
    home_production_stats(sol, hp)

Occupation `hp` (home production in the calibration) integrated over non-teachers
with P(i* = hp | X_O*). The solver taxes home-production output like market
income, so the solved `t` is not comparable to the §1 t ≈ 0.023, which is teacher
compensation over market labour income. Returns
- `share[g]`: home-production share of gender g;
- `income_share`: its share of pre-tax income;
- `t_all`: teacher wage bill over all income (the solved tax base);
- `t_market`: teacher wage bill over market income (the §1 counterpart).
"""
function home_production_stats(sol, hp::Int)
    (; hh, Φ, gr) = sol
    (; Nz, L, nXO) = gr
    hp in gr.nonteach || throw(ArgumentError("occupation $hp is not a non-teaching occupation"))
    st  = location_stats(sol)
    pHP = [[argmax_prob(gr, g, gr.XOgrid[k, g], hp) for k in 1:nXO] for g in 1:2]
    nHP = zeros(2); yHP = 0.0
    for l in 1:L, zi in 1:Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        nHP[g] += mass * integrate_nonteach(pHP[g], cm, gr)
        yHP    += mass * integrate_nonteach(pHP[g] .* hh.capO[g, l, zi, :], cm, gr)
    end
    pay, y = sum(st.pay), sum(st.y)
    return (; share = nHP ./ vec(sum(st.nT .+ st.nO; dims = 2)), income_share = yHP / y,
              t_all = pay / y, t_market = pay / (y - yHP))
end

"SD of teacher log h: `:within` work locations (mass-weighted) or `:pooled`."
function teacher_logh_sd(st, scope::Symbol)
    n = vec(sum(st.nT; dims = 1))
    if scope === :pooled
        m = sum(st.logh) / sum(n)
        return sqrt(max(sum(st.logh2) / sum(n) - m^2, 0.0))
    elseif scope === :within
        v = [n[l] > 0 ? st.logh2[l] / n[l] - (st.logh[l] / n[l])^2 : 0.0 for l in eachindex(n)]
        return sqrt(max(dot(n, v) / sum(n), 0.0))
    end
    throw(ArgumentError("teacher_sd must be :within or :pooled, got :$scope"))
end

# -----------------------------------------------------------------------------
# 3. CFR fixed-class-size teacher intervention
# -----------------------------------------------------------------------------
"""
    teacher_intervention(sol, meas; stats = location_stats(sol))

The CFR (2014b) experiment: for one of K school years a child's teacher is one SD
better in log h_T, at a fixed class size. Since h = h_T^β N^{−σ}(zϵ)^α s^φ e^η, the
child's school input rises by δ = β·sd(log h_T)/K in logs. Each child re-solves
(e, s) at Q_l e^δ with Λ, taxes and aggregates held fixed. Occupation and
work-location weights stay at the baseline. The earnings ratio at a node does not
depend on the destination, so fixing location only fixes the weights.

Returns
- `effect`: log(E y₁ / E y₀) over all children, pre-tax earnings (the target);
- `effect_nonteach`: the same over non-teachers only;
- `mean_dlog`: E[log(y₁/y₀)] over all children;
- `analytic`: β·sd/(K(1−η)), exact for non-teachers when there is no goods cost;
- `sd_logh`, `δ`.
"""
function teacher_intervention(sol, meas::Measurement; stats = location_stats(sol))
    (; hh, Φ, Pi_T, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    sd  = teacher_logh_sd(stats, meas.teacher_sd)
    δ   = p.β * sd / meas.K
    gr1 = merge(gr, (; Q = gr.Q .* exp(δ)))
    Y0 = Y1 = Y0O = Y1O = D = N = 0.0
    yT = zeros(nϵT); rT = zeros(nϵT)
    yO = zeros(nXO); rO = zeros(nXO)
    for l in 1:L, zi in 1:Nz, g in 1:2
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        for k in 1:nϵT
            e0, s0, h0 = hh.eT[g, l, zi, k], hh.sT[g, l, zi, k], hh.hT[g, l, zi, k]
            nd = teach_node(ϵTgrid[k], zi, l, g, hh.Λ, p, gr1)
            e1, _, s1, _, _ = solve_node(nd, p, log(e0), s0)
            rT[k] = (teacher_h(ϵTgrid[k], zi, l, s1, e1, p, gr1) / h0)^p.γ
            yT[k] = (1 - p.τω[p.T, g]) * h0^p.γ * dot(view(Pi_T, g, l, zi, k, :), p.κ)
        end
        for k in 1:nXO
            e0, s0 = hh.eO[g, l, zi, k], hh.sO[g, l, zi, k]
            nd = nonteach_node(XOgrid[k, g], zi, l, g, hh.Λ, p, gr1)
            e1, _, s1, _, _ = solve_node(nd, p, log(e0), s0)
            yO[k] = hh.capO[g, l, zi, k]
            rO[k] = nonteach_cap(XOgrid[k, g], zi, l, s1, e1, p, gr1) / yO[k]
        end
        y0T, y0O = integrate_choice(yT, yO, cm, gr)
        y1T, y1O = integrate_choice(yT .* rT, yO .* rO, cm, gr)
        dT, dO   = integrate_choice(log.(rT), log.(rO), cm, gr)
        nTc, nOc = integrate_choice(ones(nϵT), ones(nXO), cm, gr)
        Y0  += mass * (y0T + y0O); Y1  += mass * (y1T + y1O)
        Y0O += mass * y0O;         Y1O += mass * y1O
        D   += mass * (dT + dO);   N   += mass * (nTc + nOc)
    end
    return (; effect = log(Y1 / Y0), effect_nonteach = log(Y1O / Y0O), mean_dlog = D / N,
              analytic = δ / (1 - p.η), sd_logh = sd, δ)
end

# -----------------------------------------------------------------------------
# 4. BFM willingness to pay for school quality
# -----------------------------------------------------------------------------
"""
    school_quality_wtp(sol, meas)

BFM (2007) willingness to pay for a school whose average score is Δ =
`meas.wtp_score_sd` student SDs higher. Through the earnings bridge every child
schooled in l′ gets log h′ + x, x = χΔ, holding the child's policies fixed. For
f(h) = ((h/h̄)^{1−ψ} − 1)/(1−ψ), f(h eˣ) − f(h) = [1 + (1−ψ) f(h)]·c(x) with
c(x) = expm1((1−ψ)x)/(1−ψ), so the altruism term rises by

    ΔΛ_{l′}(z) = c(x)·[λ + (1−ψ) Λ_{l′}(z)]        (= λx at ψ = 1).

A parent in l′ with consumption C pays a = C(1 − e^{−ΔΛ/μ}), which solves
μ log(C − a) + Λ + ΔΛ = μ log C + Λ. Location choices are held fixed.

Returns `share = E[a]/E[C]` over parents (the target: \$19.70/month over mean
monthly consumption) and `share_l` by location. It also audits the consumption
support: `minC` over all (node, l′) pairs, the count of pairs at or below
`C_FLOOR` (`n_binding`), and `floor_mass`, the population share whose chosen
location has consumption at the floor. Binding pairs are usually unaffordable
moves chosen with probability ≈ 0, which smooth_log handles like log would;
`floor_mass` is what measures a real feasibility problem.
"""
function school_quality_wtp(sol, meas::Measurement)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    x  = meas.χ * meas.wtp_score_sd
    cx = is_log_kernel(p.ψ) ? x : expm1((1 - p.ψ) * x) / (1 - p.ψ)
    ΔΛ = is_log_kernel(p.ψ) ? fill(p.λ * x, Nz, L) : cx .* (p.λ .+ (1 - p.ψ) .* hh.Λ)
    frac = -expm1.(-ΔΛ ./ p.μ)                        # a/C by (parent z, l′)
    Csum = zeros(L); Asum = zeros(L)
    minC = Inf; n_binding = 0; floor_mass = 0.0
    CT = zeros(nϵT, L); CO = zeros(nXO, L)
    FT = zeros(nϵT);    FO = zeros(nXO)
    for l in 1:L, zi in 1:Nz, g in 1:2
        for k in 1:nϵT
            nd = teach_node(ϵTgrid[k], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sT[g, l, zi, k], hh.eT[g, l, zi, k], p)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            FT[k] = sum(Pi_T[g, l, zi, k, b]; init = 0.0)
            CT[k, :] .= Pi_T[g, l, zi, k, :] .* C
        end
        for k in 1:nXO
            nd = nonteach_node(XOgrid[k, g], zi, l, g, hh.Λ, p, gr)
            C  = consumption(nd, hh.sO[g, l, zi, k], hh.eO[g, l, zi, k], p)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            FO[k] = sum(Pi_O[g, l, zi, k, b]; init = 0.0)
            CO[k, :] .= Pi_O[g, l, zi, k, :] .* C
        end
        mass = Φ[zi, l] / 2
        iszero(mass) && continue
        cm = choice_maps(hh, g, l, zi, gr)
        floor_mass += mass * sum(integrate_choice(FT, FO, cm, gr))
        for lp in 1:L
            c = mass * sum(integrate_choice(CT[:, lp], CO[:, lp], cm, gr))
            Csum[lp] += c
            Asum[lp] += c * frac[zi, lp]
        end
    end
    return (; share = sum(Asum) / sum(Csum), share_l = Asum ./ Csum, minC, n_binding,
              floor_mass = floor_mass / sum(Φ), x)
end

# -----------------------------------------------------------------------------
# 5. Wage distributions: a weighted pseudo-sample of adults
#
# Quantile moments (90/10s, medians, income lines) need the distribution of wages,
# not integrals of grid quantities. Each (g, l, zi) cell is integrated on midpoints
# of a uniform grid in log ϵ_T and log X_O*, where every policy is a smooth function
# of the shock; a point carries the mass of its slice times the choice
# probabilities. Pooled over cells, the points are a weighted sample of adults.
# -----------------------------------------------------------------------------
"""
    micro_sample(sol; n = 600, groups = nothing) -> NamedTuple of vectors

The adult cross-section as a weighted pseudo-sample. Each (gender g, birth location
l, ability node zi) cell is integrated on `n` midpoints of a uniform grid in log ϵ_T
(teachers) and log X_O* (non-teachers) between the solver's grid ends; a point's
weight is its cell mass Φ[zi, l]/2 times the shock mass of its slice times:
- teachers: P(teach | ϵ_T) and the work-location logit π_{l′};
- non-teachers: P(don't teach | X_O*), P(i* ∈ group k | X_O*) and π_{l′}.
`logw` is the log pre-tax wage at the slice midpoint: log[(1−τω_{T,g}) κ_{l′} h_T^γ]
for teachers and log[(1−τω_{i*,g}) A_{i*} h_O] = log `capO` for non-teachers, which
depends neither on l′ nor on which occupation attains the max; `lo` and `hi` are its
values at the slice edges, so that a point is its mass spread uniformly over
[lo, hi] (`mquantile`). `groups` partitions the non-teaching occupations (default: one
group); `grp` is 0 for teachers and the group index otherwise. Weights sum to sum(Φ)
up to the shock mass outside the grids and the solver's quadrature error (≈1e-3).
"""
function micro_sample(sol; n = 600, groups = nothing)
    (; hh, Φ, Pi_T, Pi_O, gr) = sol
    (; Nz, L, ϵTgrid, XOgrid, dT) = gr
    groups = something(groups, [gr.nonteach])
    G   = length(groups)
    out = (; g = Int[], l = Int[], zi = Int[], lp = Int[], grp = Int[],
             logw = Float64[], lo = Float64[], hi = Float64[], w = Float64[])
    # a slice's log wage runs from e[j] to e[j+1] (edges), with midpoint value mid[j]
    function push_pts!(g, l, zi, lp, k, shift, mid, e, w)
        m = length(mid)
        append!(out.g, fill(g, m)); append!(out.l, fill(l, m)); append!(out.zi, fill(zi, m))
        append!(out.lp, fill(lp, m)); append!(out.grp, fill(k, m)); append!(out.w, w)
        append!(out.logw, shift .+ mid)
        append!(out.lo, shift .+ min.(view(e, 1:m), view(e, 2:m+1)))
        append!(out.hi, shift .+ max.(view(e, 1:m), view(e, 2:m+1)))
    end

    luT = log.(ϵTgrid)
    uT  = collect(range(luT[1], luT[end]; length = n + 1))
    umT = (uT[1:end-1] .+ uT[2:end]) ./ 2
    ϵm  = exp.(umT)
    fT  = pdf.(dT, ϵm) .* ϵm .* (uT[2] - uT[1])          # density of log ϵ_T is f(ϵ)·ϵ
    for g in 1:2
        lx = log.(XOgrid[:, g])
        v  = collect(range(lx[1], lx[end]; length = n + 1))
        vm = (v[1:end-1] .+ v[2:end]) ./ 2
        xm = exp.(vm)
        fO = [fxo(gr, g, x) * x * (v[2] - v[1]) for x in xm]
        pg = G == 1 ? ones(n, 1) :
             [sum(argmax_prob(gr, g, xm[j], i) for i in groups[k]) for j in 1:n, k in 1:G]
        for l in 1:L, zi in 1:Nz
            mass = Φ[zi, l] / 2
            iszero(mass) && continue
            cm   = choice_maps(hh, g, l, zi, gr)
            wT   = mass .* fT .* [teach_wt(cm, ϵ) for ϵ in ϵm]
            splT = Spline1D(luT, log.(hh.hT[g, l, zi, :]))
            lwT, lwTe = sol.p.γ .* splT(umT), sol.p.γ .* splT(uT)
            wO   = mass .* fO .* [nonteach_wt(cm, x) for x in xm]
            splO = Spline1D(lx, log.(hh.capO[g, l, zi, :]))
            lcO, lcOe = splO(vm), splO(v)
            for lp in 1:L
                πT = clamp.(Spline1D(luT, Pi_T[g, l, zi, :, lp]; k = 1)(umT), 0.0, 1.0)
                πO = clamp.(Spline1D(lx, Pi_O[g, l, zi, :, lp]; k = 1)(vm), 0.0, 1.0)
                push_pts!(g, l, zi, lp, 0, log((1 - sol.p.τω[sol.p.T, g]) * sol.p.κ[lp]), lwT, lwTe, wT .* πT)
                for k in 1:G
                    push_pts!(g, l, zi, lp, k, 0.0, lcO, lcOe, wO .* πO .* view(pg, :, k))
                end
            end
        end
    end
    return out
end

"""
    mquantile(ms, sel, ps)

Quantiles at probabilities `ps` of log wages over the points of `ms` selected by the
Bool vector `sel`. Each point spreads its mass uniformly over [lo, hi], its slice of
the shock grid, so the CDF F(w) = Σ w_k clamp((w − lo_k)/(hi_k − lo_k), 0, 1)/Σ w_k
is continuous and exact up to the within-slice curvature of the wage; F(w) = p is
solved by bisection.
"""
function mquantile(ms, sel, ps)
    lo, hi, w = ms.lo[sel], ms.hi[sel], ms.w[sel]
    W = sum(w)
    W > 0 || return map(_ -> NaN, ps)
    width = max.(hi .- lo, 1e-12)
    F(x) = sum(w[k] * clamp((x - lo[k]) / width[k], 0.0, 1.0) for k in eachindex(w)) / W
    a0, b0 = minimum(lo), maximum(hi)
    return map(ps) do p
        a, b = a0, b0
        for _ in 1:60
            m = (a + b) / 2
            F(m) < p ? (a = m) : (b = m)
            b - a < 1e-10 && break
        end
        (a + b) / 2
    end
end

"90/10 ratio of exp(log wage) over the points of `ms` selected by the Bool vector `sel`."
function p9010(ms, sel)
    q = mquantile(ms, sel, (0.1, 0.9))
    return exp(q[2] - q[1])
end

"Weighted mean of `x` with weights `w`."
wmean(x, w) = sum(x .* w) / sum(w)

"""
    ms_groups(gr, hp) -> groups

Non-teaching groups for `micro_sample` that split home production (`hp`) from the
market occupations: group 1 = [hp], group 2 = the rest. With `hp = 0`, one group.
"""
ms_groups(gr, hp::Int) = hp == 0 ? [gr.nonteach] : [[hp], setdiff(gr.nonteach, [hp])]

"""
Points of `ms` that are market workers: everyone but the home producers, who are the
points in non-teaching group `hpgrp` (0: home production not split out, all points).
"""
market_workers(ms, hpgrp::Int) = hpgrp == 0 ? trues(length(ms.w)) : ms.grp .!= hpgrp

"Group of home production in `micro_sample(sol; groups)`, or 0 if `hp` = 0."
hp_group(groups, hp::Int) = hp == 0 ? 0 : findfirst(grp -> grp == [hp], groups)

"""
    teacher_wage_stats(ms, weights) -> (; p9010, p9010_g)

Teachers' 90/10 of hourly (pre-tax) wages by gender (`p9010_g`, men and women, pooled
over locations as in the national ACS) and their weighted average with `weights`,
the Table 3 target.
"""
function teacher_wage_stats(ms, weights)
    r = [p9010(ms, (ms.grp .== 0) .& (ms.g .== g)) for g in 1:2]
    return (; p9010 = dot(collect(weights), r) / sum(weights), p9010_g = r)
end

"""
    income_stats(ms, hpgrp, L)

Residence (work-location) income statistics of pre-tax labour income, by l′:
- over all adults (home production included): `mean_log`, `log_mean`, `log_median`;
- over market workers (teachers included, home producers — non-teaching group
  `hpgrp` — excluded):
  `mean_log_market`, `log_median_market`.
"""
function income_stats(ms, hpgrp::Int, L)
    mkt = market_workers(ms, hpgrp)
    at(sel) = [findall(sel .& (ms.lp .== l)) for l in 1:L]
    at_mask(sel) = [sel .& (ms.lp .== l) for l in 1:L]
    all_l, mkt_l = at(trues(length(ms.w))), at(mkt)
    return (;
        mean_log          = [wmean(ms.logw[i], ms.w[i]) for i in all_l],
        log_mean          = [log(wmean(exp.(ms.logw[i]), ms.w[i])) for i in all_l],
        log_median        = [mquantile(ms, sel_l, (0.5,))[1] for sel_l in at_mask(trues(length(ms.w)))],
        mean_log_market   = [wmean(ms.logw[i], ms.w[i]) for i in mkt_l],
        log_median_market = [mquantile(ms, sel_l, (0.5,))[1] for sel_l in at_mask(mkt)])
end

"""
    score_loading(p, χ) -> c

Loading of log Q_l in the latent score signal log z + c log Q_l, c = b s_z/χ with
b = 1/(1 − η) and s_z = σξ/√(1 − ρz²) (calibration note §1, *Scale and χ*).
"""
score_loading(p::Params, χ) = p.σξ / sqrt(1 - p.ρz^2) / ((1 - p.η) * χ)

"""
    ability_moments(sol, ms, meas; hpgrp)

Moments of the latent score signal S_g = log z + (g/K) c log Q_l for a score taken at
grade g by someone schooled in l (full exposure g = K for adults), each standardized
by its SD over the cohort:
- `score_gap`: SEDA score-level gap at `meas.seda_grade`, E[S | schooled in 2] −
  E[S | schooled in 1] over the young (Φ), in cohort SDs of that latent signal (an
  observed-score gap is about √R ≈ 0.95 of it); composition plus Q (Table 4,
  descriptive);
- `wage_slope`: weighted OLS slope of log wage on the standardized adult score with
  a female dummy, over market workers (teachers included, home producers excluded),
  the model counterpart of the NLSY79 latent slope b·s_z (§4 consistency check);
- `rho_pc`: correlation between mothers' adult scores and their children's scores at
  `meas.child_grade`, one child per mother, mothers weighted by population mass; the
  child's z′ follows Πz and the child is schooled where the mother lives (§4
  consistency check);
- `rho_pc_all`: the same over parents of both genders; `c`, `sd_S` (adult scores).
"""
function ability_moments(sol, ms, meas::Measurement; hpgrp::Int)
    (; Φ, gr, p) = sol
    (; z, Q, Πz, Nz, L) = gr
    c  = score_loading(p, meas.χ)
    score(g) = [log(z[zi]) + min(g / meas.K, 1.0) * c * log(Q[l]) for zi in 1:Nz, l in 1:L]
    S, Sc, Ss = score(meas.K), score(meas.child_grade), score(meas.seda_grade)
    w  = Φ ./ sum(Φ)
    wsd(X) = (m = sum(w .* X); sqrt(sum(w .* (X .- m) .^ 2)))
    mS = sum(w .* S)
    sd = wsd(S)
    ES(l) = sum(Φ[:, l] .* Ss[:, l]) / sum(Φ[:, l])

    mkt = market_workers(ms, hpgrp)
    x   = [(S[ms.zi[k], ms.l[k]] - mS) / sd for k in eachindex(ms.w)][mkt]
    X   = hcat(ones(length(x)), x, Float64.(ms.g[mkt] .== 2))
    wk  = ms.w[mkt]
    β   = (X' * (wk .* X)) \ (X' * (wk .* ms.logw[mkt]))

    function parent_child(sel)
        # parent mass by (zi, l, l′), then the child's z′ ~ Πz and schooling location l′
        P = zeros(Nz, L, L)
        for k in findall(sel)
            P[ms.zi[k], ms.l[k], ms.lp[k]] += ms.w[k]
        end
        P ./= sum(P)
        m1 = m2 = s11 = s22 = s12 = 0.0
        for zi in 1:Nz, l in 1:L, lp in 1:L, zj in 1:Nz
            q = P[zi, l, lp] * Πz[zi, zj]
            a, b = S[zi, l], Sc[zj, lp]
            m1 += q * a; m2 += q * b; s11 += q * a^2; s22 += q * b^2; s12 += q * a * b
        end
        return (s12 - m1 * m2) / sqrt((s11 - m1^2) * (s22 - m2^2))
    end
    return (; score_gap = (ES(2) - ES(1)) / wsd(Ss), wage_slope = β[2],
              rho_pc = parent_child(ms.g .== 2), rho_pc_all = parent_child(trues(length(ms.w))), c, sd_S = sd)
end

"""
    poverty_gap(ms, pov_rate, L) -> (; gap, share, line)

Child-poverty validation (Table 4): children live with their parent in l′, and a
child is poor when the parent's pre-tax income (home production included) is below
the `pov_rate` quantile of all adults' income. Returns the suburb-minus-city gap in
the poor share (in shares, not pp), the shares by l′ and the log income line.
"""
function poverty_gap(ms, pov_rate, L)
    line  = mquantile(ms, trues(length(ms.w)), (pov_rate,))[1]
    below = clamp.((line .- ms.lo) ./ max.(ms.hi .- ms.lo, 1e-12), 0.0, 1.0)   # share of each slice below the line
    share = [sum((ms.w .* below)[ms.lp .== l]) / sum(ms.w[ms.lp .== l]) for l in 1:L]
    return (; gap = share[2] - share[1], share, line)
end

"""
    move_by_parent_income(sol, ms; nq = 3) -> Vector

Children's move rate, schooled in one location and living in the other as adults, by
`nq`-tile of the parent's pre-tax income (home production included; quantiles over
all adults): the NLSY transitions by parental income (T3). A parent at (zi, l′) has a
child of ability z′ ~ Πz[zi, :] schooled in l′, who moves with probability
1 − π̄[l′, z′, l′]. Slices straddling a cut are split uniformly (`micro_sample`).
"""
function move_by_parent_income(sol, ms; nq = 3)
    (; πbar, gr) = sol
    (; Πz, Nz) = gr
    cuts  = vcat(-Inf, collect(mquantile(ms, trues(length(ms.w)), [k / nq for k in 1:nq-1])), Inf)
    below(c) = clamp.((c .- ms.lo) ./ max.(ms.hi .- ms.lo, 1e-12), 0.0, 1.0)
    pmove = [sum(Πz[zi, zj] * (1 - πbar[lp, zj, lp]) for zj in 1:Nz) for zi in 1:Nz, lp in 1:gr.L]
    pm    = [pmove[ms.zi[k], ms.lp[k]] for k in eachindex(ms.w)]
    return map(1:nq) do q
        wq = ms.w .* (below(cuts[q+1]) .- below(cuts[q]))
        sum(wq .* pm) / sum(wq)
    end
end

# -----------------------------------------------------------------------------
# 6. All moments
# -----------------------------------------------------------------------------
"""
    model_moments(sol, meas = Measurement()) -> NamedTuple

Targeted moments (§2), in the order of the note's table:
- `male_teach_share`, `female_teach_share`: teaching share among men and women (κ̄
  and the female outside-option level ωf);
- `teacher_p9010`: teachers' 90/10 of hourly wages, the `meas.teacher_weights` average
  of men's and women's ratios (γ); `teacher_p9010_g` by gender;
- `gap_salary`: Δlog of mean teacher pay by work location, Δlog[κ_l E(h_T^γ | T, l)];
- `gap_pupils`: Δlog enrollment, log(M₂/M₁);
- `cfr_effect`: log earnings gain from the CFR intervention, vs. log(1.013);
- `wtp_share`: BFM WTP as a share of mean consumption;
- `gap_income`: Δ residence income statistic (`meas.income`); every variant is also
  returned as `gap_income_<statistic>`;
- `gap_teachers_pp`: Δlog teacher headcount per student, n_l/S_l.

Validation and diagnostics:
- `seda_growth_gap`, `gap_logQ`: the SEDA growth mapping of §3 and the quality gap;
- `score_gap`: SEDA score-level gap in cohort SDs of log z + c log Q (`ability_moments`);
- `poverty_gap`, `poverty_share_l`: Δ share of children below the `meas.pov_rate`
  income line, and the shares by location (`poverty_gap`);
- `wage_slope`, `rho_pc`, `score_loading`: the NLSY moments in the full equilibrium
  (§4 consistency check; `ability_moments`);
- `teach_share`: for the occupational block;
- `move_rate`, `p12`, `p21`: gross and directional birth → adult location flows (the
  move rate is the Table 3 residential-transition target);
- `move_rate_by_parent_income`: the move rate by parental-income tercile (validation);
- `teacher_sd`, `cfr_nonteach`, `cfr_analytic`: pieces of the CFR moment;
- `wtp_share_l`: WTP share by location;
- `s_ratio`: mean s_T / mean s_O (the §1 schooling diagnostic, T1);
- `students_per_teacher`, `t`, `Q`, `M`;
- `mean_C`, `mean_h`: drift from the frozen C̄ and h̄;
- `stationarity_dev`: max |adults_l − students_l|/students_l (quadrature accuracy);
- `minC`, `n_binding`, `floor_mass`: consumption support (see `school_quality_wtp`).
"""
function model_moments(sol, meas::Measurement = Measurement())
    (; Φ, πbar, gr, p) = sol
    (; L, Nz) = gr
    @assert L == 2 "the calibration benchmark is the city/suburb economy (L = 2)"
    gap(v) = log(v[2]) - log(v[1])

    st  = location_stats(sol)
    cfr = teacher_intervention(sol, meas; stats = st)
    wtp = school_quality_wtp(sol, meas)
    grp = ms_groups(gr, meas.hp)
    hpg = hp_group(grp, meas.hp)
    ms  = micro_sample(sol; n = meas.nquad, groups = grp)
    tw  = teacher_wage_stats(ms, meas.teacher_weights)
    inc = income_stats(ms, hpg, L)
    ab  = ability_moments(sol, ms, meas; hpgrp = hpg)
    pov = poverty_gap(ms, meas.pov_rate, L)

    nTl    = vec(sum(st.nT; dims = 1))
    adults = nTl .+ vec(sum(st.nO; dims = 1))
    meas.income in keys(inc) ||
        throw(ArgumentError("income must be one of $(keys(inc)), got :$(meas.income)"))
    (meas.hp == 0 && endswith(string(meas.income), "_market")) &&
        throw(ArgumentError("income = :$(meas.income) needs the home-production index `hp`"))
    income = inc[meas.income]

    S  = st.S
    P  = [sum(Φ[zi, l] * πbar[l, zi, lp] for zi in 1:Nz) for l in 1:L, lp in 1:L]
    P ./= sum(P; dims = 2)                      # rows integrate to 1 + O(quadrature error)
    gap_logQ = gap(gr.Q)
    ss = s_summary(sol)

    return (;
        male_teach_share   = sum(st.nT[1, :]) / sum(st.nT[1, :] .+ st.nO[1, :]),
        female_teach_share = sum(st.nT[2, :]) / sum(st.nT[2, :] .+ st.nO[2, :]),
        teacher_p9010      = tw.p9010,
        gap_salary         = gap(st.pay ./ nTl),
        gap_pupils         = gap(S),
        cfr_effect         = cfr.effect,
        wtp_share          = wtp.share,
        gap_income         = income[2] - income[1],
        gap_teachers_pp    = gap(nTl ./ S),
        seda_growth_gap    = seda_growth(gap_logQ, meas.χ, meas.K, p.η),
        gap_logQ,
        score_gap          = ab.score_gap,
        poverty_gap        = pov.gap,
        teacher_p9010_g    = tw.p9010_g,
        gap_income_mean_log          = inc.mean_log[2] - inc.mean_log[1],
        gap_income_log_mean          = inc.log_mean[2] - inc.log_mean[1],
        gap_income_log_median        = inc.log_median[2] - inc.log_median[1],
        gap_income_mean_log_market   = meas.hp == 0 ? NaN : inc.mean_log_market[2] - inc.mean_log_market[1],
        gap_income_log_median_market = meas.hp == 0 ? NaN : inc.log_median_market[2] - inc.log_median_market[1],
        wage_slope = ab.wage_slope, rho_pc = ab.rho_pc, score_loading = ab.c,
        poverty_share_l = pov.share,
        move_rate_by_parent_income = move_by_parent_income(sol, ms),
        teach_share        = sum(nTl) / sum(adults),
        move_rate          = sum(S[l] * (1 - P[l, l]) for l in 1:L) / sum(S),
        p12 = P[1, 2], p21 = P[2, 1],
        teacher_sd = cfr.sd_logh, cfr_nonteach = cfr.effect_nonteach,
        cfr_analytic = cfr.analytic,
        wtp_share_l = wtp.share_l,
        s_ratio = ss.meanT / ss.meanO,
        students_per_teacher = S ./ nTl,
        t = copy(sol.t), Q = copy(gr.Q), M = 2 .* S,
        mean_C = mean_consumption(sol), mean_h = mean_child_h(sol),
        stationarity_dev = maximum(abs.(adults .- S) ./ S),
        minC = wtp.minC, n_binding = wtp.n_binding, floor_mass = wtp.floor_mass)
end

# -----------------------------------------------------------------------------
# 7. External-block consistency (Table 2 in the full equilibrium, §4)
# -----------------------------------------------------------------------------
"""
    occupation_consistency(sol, meas; n = meas.nquad)

The Table 2 moments recomputed in the full equilibrium, where the first stage
(spatial_first_stage.jl) held them fixed: selection into teaching, Q_l, t_l and the
goods moving cost all operate (§4, *External-block consistency*).
- `shares[k, g]`: shares of the non-teaching occupations `gr.nonteach[k]` among
  non-teachers of gender g, the first stage's share targets;
- `p9010[k, g]`: within-cell 90/10 of pre-tax hourly wages (NaN for home
  production, `meas.hp`); the caller pools them with the data's wage-sample weights;
- `mean_logw[k, g]`, `teacher_mean_logw[g]`: mean log pre-tax hourly wage by cell and of
  teachers, for the untargeted mean-wage check (levels are in model units);
- `teacher_p9010_g`, `teacher_p9010`: teachers' 90/10s (`teacher_wage_stats`);
- `wage_slope`, `rho_pc`, `rho_pc_all`, `c`: ability moments (`ability_moments`),
  against the NLSY targets b·s_z and ρz.
"""
function occupation_consistency(sol, meas::Measurement; n = meas.nquad)
    (; gr) = sol
    groups = [[i] for i in gr.nonteach]
    ms  = micro_sample(sol; n, groups)
    hpg = hp_group(groups, meas.hp)
    nt  = length(groups)
    within(k, g) = (ms.grp .== k) .& (ms.g .== g)
    shares = [sum(ms.w[within(k, g)]) / sum(ms.w[(ms.grp .> 0) .& (ms.g .== g)]) for k in 1:nt, g in 1:2]
    R  = [k == hpg ? NaN : p9010(ms, within(k, g)) for k in 1:nt, g in 1:2]
    tw = teacher_wage_stats(ms, meas.teacher_weights)
    ab = ability_moments(sol, ms, meas; hpgrp = hpg)
    cellmean(sel) = (i = findall(sel); wmean(ms.logw[i], ms.w[i]))
    ml = [k == hpg ? NaN : cellmean(within(k, g)) for k in 1:nt, g in 1:2]
    tl = [cellmean((ms.grp .== 0) .& (ms.g .== g)) for g in 1:2]
    return (; occ = copy(gr.nonteach), shares, p9010 = R, mean_logw = ml, teacher_mean_logw = tl, teacher_p9010 = tw.p9010,
              teacher_p9010_g = tw.p9010_g, ab.wage_slope, ab.rho_pc, ab.rho_pc_all, ab.c)
end
