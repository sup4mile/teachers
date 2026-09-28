# =============================================================================
# Model counterparts of the calibration targets in
# data/spatial/spatial_calibration.md (§2–3; to-do item T4).
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
#   5  All moments
# =============================================================================

isdefined(@__MODULE__, :solve_ge) || include(joinpath(@__DIR__, "spatial_continuous.jl"))

# -----------------------------------------------------------------------------
# 1. Measurement assumptions
# -----------------------------------------------------------------------------
"""
    Measurement(; χ, K, teacher_sd, wtp_score_sd, income)

Assumptions that map model objects into the units of the data. None of these is
estimated; each is a documented choice (§1 and §3 of the calibration note).
- `χ = 0.13`: log adult earnings per student SD of achievement (the CFR bridge).
- `K = 5`: effective school years, used for the teacher intervention and SEDA.
- `teacher_sd = :within`: SD of teacher log h within work locations (`:pooled`
  uses the economy-wide SD).
- `wtp_score_sd = 1.0`: BFM's score unit in student SDs. Placeholder until the
  school-average-SD conversion is settled (T3).
- `income = :mean_log`: residence income statistic, E[log y] (`:log_mean` uses
  log E[y]). y is pre-tax labour income; the district proxy is a log median (T3).
"""
Base.@kwdef struct Measurement
    χ::Float64            = 0.13
    K::Float64            = 5.0
    teacher_sd::Symbol    = :within
    wtp_score_sd::Float64 = 1.0
    income::Symbol        = :mean_log
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
# 5. All moments
# -----------------------------------------------------------------------------
"""
    model_moments(sol, meas = Measurement()) -> NamedTuple

Targeted moments (§2), in the order of the note's table:
- `male_teach_share`: teaching share among men (the target for κ̄);
- `gap_salary`: Δlog of mean teacher pay by work location, Δlog[κ_l E(h_T^γ | T, l)];
- `gap_pupils`: Δlog enrollment, log(M₂/M₁);
- `cfr_effect`: log earnings gain from the CFR intervention, vs. log(1.013);
- `wtp_share`: BFM WTP as a share of mean consumption;
- `gap_income`: Δ residence income statistic (`meas.income`), vs. the 0.2708 proxy;
- `gap_teachers_pp`: Δlog teacher headcount per student, n_l/S_l.

Validation and diagnostics:
- `seda_growth_gap`, `gap_logQ`: the SEDA growth mapping of §3 and the quality gap;
- `female_teach_share`, `teach_share`: for the occupational block;
- `move_rate`, `p12`, `p21`: gross and directional birth → adult location flows;
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

    nTl    = vec(sum(st.nT; dims = 1))
    adults = nTl .+ vec(sum(st.nO; dims = 1))
    income = meas.income === :mean_log ? st.logy ./ adults :
             meas.income === :log_mean ? log.(st.y ./ adults) :
             throw(ArgumentError("income must be :mean_log or :log_mean, got :$(meas.income)"))

    S  = st.S
    P  = [sum(Φ[zi, l] * πbar[l, zi, lp] for zi in 1:Nz) for l in 1:L, lp in 1:L]
    P ./= sum(P; dims = 2)                      # rows integrate to 1 + O(quadrature error)
    gap_logQ = gap(gr.Q)
    ss = s_summary(sol)

    return (;
        male_teach_share   = sum(st.nT[1, :]) / sum(st.nT[1, :] .+ st.nO[1, :]),
        gap_salary         = gap(st.pay ./ nTl),
        gap_pupils         = gap(S),
        cfr_effect         = cfr.effect,
        wtp_share          = wtp.share,
        gap_income         = income[2] - income[1],
        gap_teachers_pp    = gap(nTl ./ S),
        seda_growth_gap    = seda_growth(gap_logQ, meas.χ, meas.K, p.η),
        gap_logQ,
        female_teach_share = sum(st.nT[2, :]) / sum(st.nT[2, :] .+ st.nO[2, :]),
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
