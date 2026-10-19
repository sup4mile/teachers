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
#     Under the residence/workplace split (`Params.split`) residence r and workplace w
#     come apart. Children, taxes, income, poverty and composition stay by residence
#     (`ms.lp` = r in the pseudo-sample; `ms.lw` = w), while teachers' counts, pay and quality
#     are by workplace (`location_stats`); the solution then has `Pi_T = Pi_O = nothing`.
#   - Students in l′ are the young born there, S_{l′} = Σ_z Φ[z, l′] = M_{l′}/2.
#   - Spatial gaps are location 2 (suburb) minus location 1 (city), in logs.
#   - With `Params.college` (spatial_continuous.jl §3c) every cell is a mixture of two education
#     branches, d = 1 (BA) and d = 0, with weights P(d | g, l, z): the integrals loop over
#     `cell_branches`, `micro_sample` carries d, and §6b adds the education moments.
#
# Layout
#   1  Measurement assumptions
#   2  Work-location aggregates (residence and workplace under the split)
#   3  CFR fixed-class-size teacher intervention
#   4  BFM willingness to pay for school quality, overall and by parental income, and the
#      value of the quality gap ΔΛ
#   5  Wage distributions: a weighted pseudo-sample of adults (one point set per residence ×
#      workplace under the split), commuting shares and the option value of the workplace
#   6  All moments
#   6b Education moments (`college_stats`) and the pooled within-occupation 90/10
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
- `pooled_p9010 = false`: compute the pooled within-occupation 90/10 of non-teachers' wages
  (`pooled_p9010`, a slow second pseudo-sample; NaN otherwise), the σϵ target of the joint
  estimation; `wage_sample = :acs2009_13` picks the wage-sample weights (`WAGE_SAMPLE_COUNTS`).
- `cfr_rho = 1.0` (round 2026-10-16), 0 < ρ ≤ 1: the correlation of teachers' log effectiveness with
  log h_eff (the model's teacher quality, `effective_h`) in the CFR bridge. CFR's one SD is one SD of value
  added, log a = log h_eff + ν with ν independent of everything in the model (absorbed in the level of Q), so
  sd(log a) = sd(log h_eff | workplace)/ρ and the intervention raises the child's school input by
  δ = β sd(log a)/K (`teacher_intervention`; only β times the loading of h_eff on a is identified, so the loading
  is one). It is a measurement concept like χ and K, not a parameter. 1.0: a is h_eff, and the model without
  it is reproduced bit for bit (`sd/1.0 == sd`). `--measure cfr_rho=0.3`.
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
    pooled_p9010::Bool    = false
    wage_sample::Symbol   = :acs2009_13
    cfr_rho::Float64      = 1.0
end

"Copy of `m` with the named fields replaced."
with_measurement(m::Measurement; kwargs...) =
    Measurement(; merge(NamedTuple(k => getfield(m, k) for k in fieldnames(Measurement)), NamedTuple(kwargs))...)

"""
Wage-sample counts of the 19 market occupations by gender, the weights of the pooled 90/10
(ACS wage sample: ages 25–34, 48+ weeks, 30+ hours; home production and teachers have none).
`acs2009_13` is `ACS_WAGES.n_men`, `.n_women` of spatial_first_stage.jl (workbook sheet
`occ_gender_weights`), `acs2016_19` is `ACS_2016_19.n_men`, `.n_women` of spatial_calibrate.jl.
"""
const WAGE_SAMPLE_COUNTS = (;
    acs2009_13 = (; men = [51028, 23030, 33560, 20877, 9642, 10527, 4292, 2850, 27754, 54653, 41631, 23058,
                           35132, 18053, 76036, 11654, 17580, 27169, 27228],
                    women = [44932, 28486, 9787, 24365, 8845, 60037, 4239, 8960, 22118, 42370, 92826, 5479,
                             40793, 2646, 2016, 3371, 4749, 6584, 2122]),
    acs2016_19 = (; men = [53404, 23404, 36803, 21098, 8445, 12342, 3687, 3145, 29200, 46931, 39453, 19973,
                           33025, 15218, 69540, 11320, 16203, 26618, 25735],
                    women = [50147, 25760, 11714, 25351, 8446, 61176, 3673, 10106, 21434, 36797, 74755, 4749,
                             37258, 3088, 2204, 3668, 4368, 7257, 2461]))

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
- `logh[l′]`, `logh2[l′]`: Σ log h_T and Σ (log h_T)² over teachers (`gap_logh_teacher`);
- `logh_eff[l′]`, `logh2_eff[l′]`: the same of their effectiveness h_T z^{−(1−ωz)α} (`effective_h`,
  which is h_T unless ωz ≠ 1): the SD behind the CFR experiment (`teacher_logh_sd`);
- `logy[l′]`, `y[l′]`: Σ log y and Σ y over all workers, y = pre-tax labour income;
- `S[l′]`: students, the young born in l′;
- `nT_res[g, l′]`, `nO_work[g, l′]`: teachers by residence and non-teachers by workplace, equal to
  `nT` and `nO` in the bundled model.
Under the split the teachers' `nT`, `pay`, `logh`, `logh2`, `logh_eff` and `logh2_eff` are by workplace w,
and `nO`, `logy`, `y` and `S` by residence r (income is the pre-commuting income of the residents, which the
tax base is), with `nT_res` and `nO_work` the other way round.
Under `p.college` the cells mix the education branches with the weights P(d | g, l, z), and a
teacher's pay carries the pay factor of her branch (`pay_factor`).
"""
function location_stats(sol)
    sol.p.split && return location_stats_split(sol)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L) = gr
    nT = zeros(2, L); nO = zeros(2, L)
    pay = zeros(L); logh = zeros(L); logh2 = zeros(L); logh_eff = zeros(L); logh2_eff = zeros(L)
    logy = zeros(L); y = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_T, Pi_O, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm   = choice_maps(br.h, g, l, zi, gr)
        lhT  = log.(br.h.hT[g, l, zi, :])
        lhE  = p.ωz == 1.0 ? lhT : log.(effective_h(br.h.hT[g, l, zi, :], zi, p, gr))
        lpf  = log.(teacher_pay_factor(br.h.sT[g, l, zi, :], p, br.d))
        yO   = br.h.capO[g, l, zi, :]
        lyO  = log.(yO)
        for lp in 1:L
            πT  = br.πT[g, l, zi, :, lp]
            πO  = br.πO[g, l, zi, :, lp]
            lyT = log((1 - p.τω[p.T, g]) * p.κ[lp]) .+ p.γ .* lhT .+ lpf
            payT = mass * integrate_teach(πT .* exp.(lyT), cm, gr)
            nT[g, lp] += mass * integrate_teach(πT, cm, gr)
            nO[g, lp] += mass * integrate_nonteach(πO, cm, gr)
            pay[lp]   += payT
            logh[lp]  += mass * integrate_teach(πT .* lhT, cm, gr)
            logh2[lp] += mass * integrate_teach(πT .* lhT .^ 2, cm, gr)
            logh_eff[lp]  += mass * integrate_teach(πT .* lhE, cm, gr)
            logh2_eff[lp] += mass * integrate_teach(πT .* lhE .^ 2, cm, gr)
            logy[lp]  += mass * (integrate_teach(πT .* lyT, cm, gr) +
                                 integrate_nonteach(πO .* lyO, cm, gr))
            y[lp]     += payT + mass * integrate_nonteach(πO .* yO, cm, gr)
        end
    end
    return (; nT, nO, pay, logh, logh2, logh_eff, logh2_eff, logy, y, S = vec(sum(Φ; dims = 1)),
              nT_res = nT, nO_work = nO)
end

"`location_stats` under the split: the joint probabilities `Pi_TW`, `Pi_OW` give workplace and residence separately."
function location_stats_split(sol)
    (; hh, Φ, Pi_TW, Pi_OW, gr, p) = sol
    (; Nz, L) = gr
    nT = zeros(2, L); nO = zeros(2, L); nT_res = zeros(2, L); nO_work = zeros(2, L)
    pay = zeros(L); logh = zeros(L); logh2 = zeros(L); logh_eff = zeros(L); logh2_eff = zeros(L)
    logy = zeros(L); y = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_TW, Pi_OW, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm   = choice_maps(br.h, g, l, zi, gr)
        lhT  = log.(br.h.hT[g, l, zi, :])
        lhE  = p.ωz == 1.0 ? lhT : log.(effective_h(br.h.hT[g, l, zi, :], zi, p, gr))
        lpf  = log.(teacher_pay_factor(br.h.sT[g, l, zi, :], p, br.d))
        yO   = br.h.capO[g, l, zi, :]
        lyO  = log.(yO)
        for r in 1:L, w in 1:L
            πT  = br.πT[g, l, zi, :, r, w]
            πO  = br.πO[g, l, zi, :, r, w]
            lyT = log((1 - p.τω[p.T, g]) * p.κ[w]) .+ p.γ .* lhT .+ lpf
            payT = mass * integrate_teach(πT .* exp.(lyT), cm, gr)
            mT   = mass * integrate_teach(πT, cm, gr)
            mO   = mass * integrate_nonteach(πO, cm, gr)
            nT[g, w] += mT;  nT_res[g, r] += mT
            nO[g, r] += mO;  nO_work[g, w] += mO
            pay[w]   += payT
            logh[w]  += mass * integrate_teach(πT .* lhT, cm, gr)
            logh2[w] += mass * integrate_teach(πT .* lhT .^ 2, cm, gr)
            logh_eff[w]  += mass * integrate_teach(πT .* lhE, cm, gr)
            logh2_eff[w] += mass * integrate_teach(πT .* lhE .^ 2, cm, gr)
            logy[r]  += mass * (integrate_teach(πT .* lyT, cm, gr) +
                                integrate_nonteach(πO .* lyO, cm, gr))
            y[r]     += payT + mass * integrate_nonteach(πO .* yO, cm, gr)
        end
    end
    return (; nT, nO, pay, logh, logh2, logh_eff, logh2_eff, logy, y, S = vec(sum(Φ; dims = 1)), nT_res, nO_work)
end

"P(i* = i | X_O* = x): the chance that non-teaching occupation `i` attains the max."
function argmax_prob(gr, g, x, i)
    w(j) = pdf(gr.dO[j, g], x) / max(cdf(gr.dO[j, g], x), 1e-300)
    den = sum(w, gr.nonteach)
    return den > 0.0 ? w(i) / den : 0.0
end

"""
    home_production_stats(sol, hp; stats = location_stats(sol))

Occupation `hp` (home production in the calibration) integrated over non-teachers
with P(i* = hp | X_O*). The solver taxes home-production output like market
income, so the solved `t` is not comparable to the §1 t ≈ 0.023, which is teacher
compensation over market labour income. Returns
- `share[g]`: home-production share of gender g;
- `income_share`: its share of pre-tax income;
- `income_l[l′]`: home-production income by work location;
- `t_all`: teacher wage bill over all income (the solved tax base);
- `t_market`: teacher wage bill over market income. The model pays salaries only, so its
  like-for-like data comparator is ≈ 0.016 (teacher salaries over market labour income);
  the §1 0.023 includes benefits;
- `t_market_l[l′]`: the same by work location (the solved `t` is the all-income version);
- `n_l[l′]`: mass of home producers by location.
Under the split the locations are residences: `Pi_Or` stands in for `Pi_O`, so that `income_l`
and `n_l` are by residence, and `t_market_l` is the pay of the teachers who work in l′ over the
market income of the residents of l′, the budget the tax rate balances.
"""
function home_production_stats(sol, hp::Int; stats = location_stats(sol))
    sol.Pi_O === nothing &&
        return home_production_stats(merge(sol, (; Pi_T = sol.Pi_Tr, Pi_O = sol.Pi_Or)), hp; stats)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L, nXO) = gr
    hp in gr.nonteach || throw(ArgumentError("occupation $hp is not a non-teaching occupation"))
    pHP = [[argmax_prob(gr, g, gr.XOgrid[k, g], hp) for k in 1:nXO] for g in 1:2]
    nHP = zeros(2); yHP = zeros(L); nHPl = zeros(L)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_T, Pi_O, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        nHP[g] += mass * integrate_nonteach(pHP[g], cm, gr)
        for lp in 1:L
            yHP[lp] += mass * integrate_nonteach(pHP[g] .* br.h.capO[g, l, zi, :] .* br.πO[g, l, zi, :, lp], cm, gr)
            nHPl[lp] += mass * integrate_nonteach(pHP[g] .* br.πO[g, l, zi, :, lp], cm, gr)
        end
    end
    pay, y = sum(stats.pay), sum(stats.y)
    return (; share = nHP ./ vec(sum(stats.nT .+ stats.nO; dims = 2)), income_share = sum(yHP) / y,
              income_l = yHP, t_all = pay / y, t_market = pay / (y - sum(yHP)),
              t_market_l = stats.pay ./ (stats.y .- yHP), n_l = nHPl)
end

"SD of teacher log effectiveness (log h_T unless ωz ≠ 1): `:within` work locations (mass-weighted) or `:pooled`."
function teacher_logh_sd(st, scope::Symbol)
    n = vec(sum(st.nT; dims = 1))
    if scope === :pooled
        m = sum(st.logh_eff) / sum(n)
        return sqrt(max(sum(st.logh2_eff) / sum(n) - m^2, 0.0))
    elseif scope === :within
        v = [n[l] > 0 ? st.logh2_eff[l] / n[l] - (st.logh_eff[l] / n[l])^2 : 0.0 for l in eachindex(n)]
        return sqrt(max(dot(n, v) / sum(n), 0.0))
    end
    throw(ArgumentError("teacher_sd must be :within or :pooled, got :$scope"))
end

"""
    teacher_sorting_stats(sol, meas) -> NamedTuple

Who teaches where (round 2026-10-16; a diagnostic of the teacher-quality route, no target uses it). One pass over the
cells (birth location l, ability z, gender g, education branch) in the loop of `location_stats_split` (and of
`location_stats` in the bundled model, where r = w = the work location, so only the diagonal is filled), accumulating
mass-weighted sums over the teachers by residence r and workplace w as L × L matrices indexed `[r, w]`:
- `m`: mass; `lh`, `le`: Σ log h_T and Σ log h_eff (`effective_h`, the effectiveness behind Q and the CFR
  experiment); `hb`: Σ h_eff^{β/σ}, what H̃_T aggregates;
- `lpay`, `pay`: Σ log and Σ level of pre-tax pay (1 − τω) κ_w h_T^γ × pay factor (the `lyT` of `location_stats`);
- `yh`: Σ ln ŷ, the predicted log income of the teacher's node (`predicted_log_income`, at her birth location,
  whatever she chooses and whatever `Params.lyT_fixed` is);
- `dQ`, `dz`, `dϵ`, `ds`, `de`: Σ of the five terms of log h_T = log Q_l + α log z + α log ϵ_T + φ log s + η log e
  (`teacher_h`): the school quality of the birth location, ability, the occupational shock, schooling and goods
  investment; `born`: the mass born in location 2 (l = 2);
and, for the non-teaching market workers by residence r (home producers, weight P(i* = hp | X_O*) when `meas.hp ≠ 0`,
are left out), the vectors `mO` (mass), `lyO` (Σ log income, the `capO` of `location_stats`) and `yhO` (Σ ln ŷ of
`nonteach_node`). `resid` is the largest node-level distance between log h_T and the sum of its five terms (round-off).
The integrals are those of `location_stats` (`integrate_teach` of the same integrands for the masses and log h_T),
so that the sums over r and over w reproduce its `nT`, `logh` and `pay` up to the order of the additions. The pieces
of the decomposition are integrated separately and sum to the total up to the quadrature tolerance, not exactly
(`sorting_moments`, `decomp_resid_max`).
"""
function teacher_sorting_stats(sol, meas::Measurement)
    (; hh, Φ, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    PT, PO = p.split ? (sol.Pi_TW, sol.Pi_OW) : (sol.Pi_T, sol.Pi_O)
    pairs_rw = p.split ? [(r, w) for r in 1:L for w in 1:L] : [(w, w) for w in 1:L]
    expo = p.β / p.σ
    hp_ok = meas.hp != 0 && meas.hp in gr.nonteach
    wmkt = [hp_ok ? [1 - argmax_prob(gr, g, XOgrid[k, g], meas.hp) for k in 1:nXO] : ones(nXO) for g in 1:2]
    m = zeros(L, L); lh = zeros(L, L); le = zeros(L, L); hb = zeros(L, L)
    lpay = zeros(L, L); pay = zeros(L, L); yh = zeros(L, L)
    dQ = zeros(L, L); dz = zeros(L, L); dϵ = zeros(L, L); ds = zeros(L, L); de = zeros(L, L); born = zeros(L, L)
    mO = zeros(L); lyO = zeros(L); yhO = zeros(L)
    resid = 0.0
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, PT, PO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        hbr, d = br.h, br.d
        cm   = choice_maps(hbr, g, l, zi, gr)
        hT   = hbr.hT[g, l, zi, :]
        lhT  = log.(hT)
        hE   = effective_h(hT, zi, p, gr)
        lhE  = p.ωz == 1.0 ? lhT : log.(hE)
        hEx  = hE .^ expo
        lpf  = log.(teacher_pay_factor(hbr.sT[g, l, zi, :], p, d))
        cQ, cz = log(gr.Q[l]), p.α * log(gr.z[zi])
        t_ϵ, t_s, t_e = p.α .* gr.lϵT, p.φ .* log.(hbr.sT[g, l, zi, :]), p.η .* log.(hbr.eT[g, l, zi, :])
        resid = max(resid, maximum(abs, lhT .- (cQ .+ cz .+ t_ϵ .+ t_s .+ t_e)))
        qyh  = [predicted_log_income(teach_node(ϵTgrid[k], zi, l, g, hbr.Λ, p, gr; d), p, gr) for k in 1:nϵT]
        for (r, w) in pairs_rw
            πT  = p.split ? br.πT[g, l, zi, :, r, w] : br.πT[g, l, zi, :, w]
            lyT = log((1 - p.τω[p.T, g]) * p.κ[w]) .+ p.γ .* lhT .+ lpf
            mT  = mass * integrate_teach(πT, cm, gr)
            m[r, w]    += mT
            lh[r, w]   += mass * integrate_teach(πT .* lhT, cm, gr)
            le[r, w]   += mass * integrate_teach(πT .* lhE, cm, gr)
            hb[r, w]   += mass * integrate_teach(πT .* hEx, cm, gr)
            lpay[r, w] += mass * integrate_teach(πT .* lyT, cm, gr)
            pay[r, w]  += mass * integrate_teach(πT .* exp.(lyT), cm, gr)
            yh[r, w]   += mass * integrate_teach(πT .* qyh, cm, gr)
            dQ[r, w]   += cQ * mT
            dz[r, w]   += cz * mT
            dϵ[r, w]   += mass * integrate_teach(πT .* t_ϵ, cm, gr)
            ds[r, w]   += mass * integrate_teach(πT .* t_s, cm, gr)
            de[r, w]   += mass * integrate_teach(πT .* t_e, cm, gr)
            l == 2 && (born[r, w] += mT)
        end
        lyOk = log.(hbr.capO[g, l, zi, :])
        qyhO = [predicted_log_income(nonteach_node(XOgrid[k, g], zi, l, g, hbr.Λ, p, gr; d), p, gr) for k in 1:nXO]
        for r in 1:L
            πO = p.split ? vec(sum(br.πO[g, l, zi, :, r, :]; dims = 2)) : br.πO[g, l, zi, :, r]
            q  = πO .* wmkt[g]
            mO[r]  += mass * integrate_nonteach(q, cm, gr)
            lyO[r] += mass * integrate_nonteach(q .* lyOk, cm, gr)
            yhO[r] += mass * integrate_nonteach(q .* qyhO, cm, gr)
        end
    end
    return (; m, lh, le, hb, lpay, pay, yh, dQ, dz, dϵ, ds, de, born, mO, lyO, yhO, resid)
end

"""
    sorting_moments(ts, sol) -> NamedTuple

The `model_moments` keys of round 2026-10-16 that come from the sums of `teacher_sorting_stats` (`ts`), each a suburb
minus city contrast (index 2 minus 1) or a level; every entry is finite or NaN (an empty cell, a location without
teachers, or the contrasts that need a workplace for both residences in the bundled model). By residence (`_res`) and
by workplace (`_work`):
- `gap_logh_teacher_res`, `gap_logh_teacher_work`: Δ mean log h_T (the work version is `gap_logh_teacher` of
  `model_moments` up to the order of the additions); `gap_logheff_teacher_res/_work`: Δ mean log h_eff;
  `gap_Heff_teacher_res/_work`: Δ log E[h_eff^{β/σ}] (log ratio of means);
- `gap_logpay_teacher_res/_work`: Δ mean log pre-tax pay; `gap_salary_res`: Δ log mean pay by residence (the workplace
  version is `gap_salary`);
- `gap_lyhat_teacher_res/_work`: Δ mean ln ŷ of the teachers' nodes; `mean_log_yhat_teach`: the mean ln ŷ of all
  teachers (the level `Params.lyT_fixed` is set to in the diagnostic T-b); `mean_log_yhat_mkt`: the mean ln ŷ of
  all market workers (teachers and non-teaching market workers: home producers excluded when `meas.hp ≠ 0`), and
  `mean_log_yhat_nonteach_mkt` that of the non-teaching ones;
- `gap_logh_teacher_{res,work}_{Q,z,eps,s,e}`: the contribution of each term of log h_T = log Q_l + α log z + α log ϵ +
  φ log s + η log e to the Δ mean log h_T (Δ of the mean of the term); `decomp_resid_max`: the largest distance, in log
  points, between a cell's mean log h_T and the sum of its five terms (integrated separately, so it is the quadrature
  tolerance, not zero); there is no 'other' term, `teacher_h` having no further factor;
- `share_born_suburb_res/_work`: the share of the teachers (by residence, by workplace) born in location 2;
- levels `teacher_logh_res/_work`, `teacher_logpay_res/_work`, `teacher_lyhat_res/_work`, `teacher_Heff_res/_work`
  (mean log h_T, mean log pay, mean ln ŷ and E[h_eff^{β/σ}] by residence and workplace) and
  `nonteach_mkt_logy_res`, `nonteach_mkt_lyhat_res`, `nonteach_mkt_mass_res` (non-teaching market workers by residence);
- `gap_lyhat_mkt_nonteach_res`, `gap_logy_mkt_nonteach_res`: Δ mean ln ŷ and Δ mean log pre-tax income of the
  non-teaching market workers by residence;
- `gap_logQ_teach_quality` = σ Δ log E[h_eff^{β/σ} | teacher, w] and `gap_logQ_teach_staffing` = σ Δ log (2 n_w/M_w), with
  n_w the teachers working at w and M_w = 2 S_w the pupils there (the normalisation of `QT` in `build_grids`; the
  second is σ times `gap_teachers_pp`): the two parts of Δ log Q_T = `gap_logQ_teach`, and
  `gap_logQ_teach_check` = `gap_logQ_teach` − quality − staffing, which is at the GE tolerance (the grid is built at the
  aggregates of the last iteration, the sums here are at the choices);
- the (r, w) cells flattened column-major, element (w − 1) L + r: `teacher_rw_mass`, `teacher_rw_logh`, `teacher_rw_logpay`;
- `gap_logpay_teacher_res_fixedw`, `gap_logh_teacher_res_fixedw`: the residence contrast at a given workplace,
  Σ_w m_w (x[2, w] − x[1, w]) / Σ_w m_w over the workplaces with teachers of both residences (m_w the teacher mass at w);
  NaN in the bundled model, where nobody crosses.
"""
function sorting_moments(ts, sol)
    (; Φ, gr, p) = sol
    L = gr.L
    gap(v) = v[2] - v[1]
    slog(x) = x > 0 && isfinite(x) ? log(x) : NaN
    rows(x) = vec(sum(x; dims = 2))                 # by residence r
    cols(x) = vec(sum(x; dims = 1))                 # by workplace w
    avg(num, den) = [d > 0 ? n / d : NaN for (n, d) in zip(num, den)]
    mr, mw = rows(ts.m), cols(ts.m)
    res_work(x) = (avg(rows(x), mr), avg(cols(x), mw))
    lh_r, lh_w = res_work(ts.lh)
    le_r, le_w = res_work(ts.le)
    lp_r, lp_w = res_work(ts.lpay)
    yh_r, yh_w = res_work(ts.yh)
    He_r, He_w = res_work(ts.hb)
    pay_r = avg(rows(ts.pay), mr)
    Q_r, Q_w = res_work(ts.dQ); z_r, z_w = res_work(ts.dz); ϵ_r, ϵ_w = res_work(ts.dϵ)
    s_r, s_w = res_work(ts.ds); e_r, e_w = res_work(ts.de)
    bo_r, bo_w = res_work(ts.born)
    parts = ts.dQ .+ ts.dz .+ ts.dϵ .+ ts.ds .+ ts.de
    cell_resid = [ts.m[r, w] > 0 ? abs(ts.lh[r, w] - parts[r, w]) / ts.m[r, w] : 0.0 for r in 1:L, w in 1:L]
    cellmean(x) = [ts.m[r, w] > 0 ? x[r, w] / ts.m[r, w] : NaN for r in 1:L, w in 1:L]
    # the residence contrast at a given workplace, averaged over the workplaces with teachers of both residences
    tol = 1e-12 * sum(ts.m)
    ok  = [w for w in 1:L if ts.m[1, w] > tol && ts.m[2, w] > tol]
    function fixedw(x)
        isempty(ok) && return NaN
        mwk = [ts.m[1, w] + ts.m[2, w] for w in ok]
        return sum(mwk[i] * (x[2, w] - x[1, w]) for (i, w) in enumerate(ok)) / sum(mwk)
    end
    # the teacher-quality index: quality of the teachers and staffing (QT = (2 H̃_T/M)^σ with H̃_T = n E[h_eff^{β/σ}])
    S    = vec(sum(Φ; dims = 1))
    stf  = [slog(2 * mw[w] / (2 * S[w])) for w in 1:L]
    g_q  = p.σ * gap(slog.(He_w))
    g_s  = p.σ * gap(stf)
    g_QT = log(gr.QT[2]) - log(gr.QT[1])
    nm_ly, nm_yh = avg(ts.lyO, ts.mO), avg(ts.yhO, ts.mO)
    tot_t, tot_o = sum(ts.m), sum(ts.mO)
    return (;
        gap_logh_teacher_res = gap(lh_r), gap_logh_teacher_work = gap(lh_w),
        gap_logheff_teacher_res = gap(le_r), gap_logheff_teacher_work = gap(le_w),
        gap_Heff_teacher_res = gap(slog.(He_r)), gap_Heff_teacher_work = gap(slog.(He_w)),
        gap_logpay_teacher_res = gap(lp_r), gap_logpay_teacher_work = gap(lp_w),
        gap_salary_res = gap(slog.(pay_r)),
        gap_lyhat_teacher_res = gap(yh_r), gap_lyhat_teacher_work = gap(yh_w),
        mean_log_yhat_teach = tot_t > 0 ? sum(ts.yh) / tot_t : NaN,
        mean_log_yhat_mkt = tot_t + tot_o > 0 ? (sum(ts.yh) + sum(ts.yhO)) / (tot_t + tot_o) : NaN,
        mean_log_yhat_nonteach_mkt = tot_o > 0 ? sum(ts.yhO) / tot_o : NaN,
        gap_lyhat_mkt_nonteach_res = gap(nm_yh), gap_logy_mkt_nonteach_res = gap(nm_ly),
        gap_logh_teacher_res_Q = gap(Q_r), gap_logh_teacher_work_Q = gap(Q_w),
        gap_logh_teacher_res_z = gap(z_r), gap_logh_teacher_work_z = gap(z_w),
        gap_logh_teacher_res_eps = gap(ϵ_r), gap_logh_teacher_work_eps = gap(ϵ_w),
        gap_logh_teacher_res_s = gap(s_r), gap_logh_teacher_work_s = gap(s_w),
        gap_logh_teacher_res_e = gap(e_r), gap_logh_teacher_work_e = gap(e_w),
        decomp_resid_max = max(maximum(cell_resid), ts.resid),
        share_born_suburb_res = bo_r, share_born_suburb_work = bo_w,
        teacher_logh_res = lh_r, teacher_logh_work = lh_w, teacher_logpay_res = lp_r, teacher_logpay_work = lp_w,
        teacher_lyhat_res = yh_r, teacher_lyhat_work = yh_w, teacher_Heff_res = He_r, teacher_Heff_work = He_w,
        nonteach_mkt_logy_res = nm_ly, nonteach_mkt_lyhat_res = nm_yh, nonteach_mkt_mass_res = copy(ts.mO),
        gap_logQ_teach_quality = g_q, gap_logQ_teach_staffing = g_s, gap_logQ_teach_check = g_QT - g_q - g_s,
        teacher_rw_mass = vec(copy(ts.m)), teacher_rw_logh = vec(cellmean(ts.lh)), teacher_rw_logpay = vec(cellmean(ts.lpay)),
        gap_logpay_teacher_res_fixedw = fixedw(cellmean(ts.lpay)), gap_logh_teacher_res_fixedw = fixedw(cellmean(ts.lh)))
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
With `meas.cfr_rho` = ρ < 1 (round 2026-10-16) the SD that CFR's "one SD of value added" is in the model's
units is sd(log a) = sd(log h_eff)/ρ, with log a the teacher's effectiveness and ρ its correlation with log h_eff
(`Measurement`): δ = β sd(log h_eff)/(ρK). ρ = 1 (the default) is the experiment above, bit for bit.

Returns
- `effect`: log(E y₁ / E y₀) over all children, pre-tax earnings (the target);
- `effect_nonteach`: the same over non-teachers only;
- `mean_dlog`: E[log(y₁/y₀)] over all children;
- `analytic`: β·sd(log a)/(K(1−η)), exact for non-teachers when there is no goods cost;
- `sd_logh`: sd(log h_eff) (`teacher_logh_sd`); `sd_eff`: the SD of log a that the intervention uses, `sd_logh/ρ`
  (equal to `sd_logh` at ρ = 1); `δ`.
Under the split a teacher's pay depends on where she works, so the work-location weights are the
workplace marginal of the joint probabilities.
Under `p.college` both education branches are re-solved with P(d | g, l, z) held fixed (§3c);
a teacher's income carries her pay factor, which moves with her time share (`pay_factor`), and
the SD of teacher quality is that of the effectiveness `effective_h` (h_T unless ωz ≠ 1).
"""
function teacher_intervention(sol, meas::Measurement; stats = location_stats(sol))
    0 < meas.cfr_rho ≤ 1 || throw(ArgumentError("cfr_rho must lie in (0, 1], got $(meas.cfr_rho)"))
    sol.Pi_T === nothing &&
        return teacher_intervention(merge(sol, (; Pi_T = workplace_probs(sol.Pi_TW), Pi_O = workplace_probs(sol.Pi_OW))),
                                    meas; stats)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    sd  = teacher_logh_sd(stats, meas.teacher_sd)
    sd_a = sd / meas.cfr_rho                    # SD of log effectiveness a, corr(log a, log h_eff) = ρ; ρ = 1: sd_a === sd
    δ   = p.β * sd_a / meas.K
    gr1 = merge(gr, (; Q = gr.Q .* exp(δ)))
    Y0 = Y1 = Y0O = Y1O = D = N = 0.0
    yT = zeros(nϵT); rT = zeros(nϵT)
    yO = zeros(nXO); rO = zeros(nXO)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_T, Pi_O, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        hb, d = br.h, br.d
        cm = choice_maps(hb, g, l, zi, gr)
        for k in 1:nϵT
            e0, s0, h0 = hb.eT[g, l, zi, k], hb.sT[g, l, zi, k], hb.hT[g, l, zi, k]
            nd = teach_node(ϵTgrid[k], zi, l, g, hb.Λ, p, gr1; d)
            e1, _, s1, _, _ = solve_node(nd, p, log(e0), s0)
            rT[k] = (teacher_h(ϵTgrid[k], zi, l, s1, e1, p, gr1) / h0)^p.γ * (pay_factor(s1, p, d) / pay_factor(s0, p, d))
            yT[k] = (1 - p.τω[p.T, g]) * h0^p.γ * dot(view(br.πT, g, l, zi, k, :), p.κ) * pay_factor(s0, p, d)
        end
        for k in 1:nXO
            e0, s0 = hb.eO[g, l, zi, k], hb.sO[g, l, zi, k]
            nd = nonteach_node(XOgrid[k, g], zi, l, g, hb.Λ, p, gr1; d)
            e1, _, s1, _, _ = solve_node(nd, p, log(e0), s0)
            yO[k] = hb.capO[g, l, zi, k]
            rO[k] = nonteach_cap(XOgrid[k, g], zi, l, s1, e1, p, gr1) * credential(p, d) / yO[k]
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
              analytic = δ / (1 - p.η), sd_logh = sd, sd_eff = sd_a, δ)
end

# -----------------------------------------------------------------------------
# 4. BFM willingness to pay for school quality, overall and by parental income, and the
#    value of the quality gap ΔΛ
# -----------------------------------------------------------------------------
"""
    wtp_fraction(sol, meas) -> (; ΔΛ, frac, x)

The BFM school-quality step of `school_quality_wtp` and `wtp_by_parent_income`: through
the earnings bridge every child schooled in l′ gets log h′ + x, x = χΔ with Δ =
`meas.wtp_score_sd` student SDs, holding the child's policies fixed. For
f(h) = ((h/h̄)^{1−ψ} − 1)/(1−ψ), f(h eˣ) − f(h) = [1 + (1−ψ) f(h)]·c(x) with
c(x) = expm1((1−ψ)x)/(1−ψ), so the altruism term rises by

    ΔΛ_{l′}(z) = c(x)·[λ + (1−ψ) Λ_{l′}(z)]        (= λx at ψ = 1).

A parent in l′ with consumption C pays a = C(1 − e^{−ΔΛ/μ}), which solves
μ log(C − a) + Λ + ΔΛ = μ log C + Λ; `frac[z, l′]` is a/C. Location choices are held
fixed. C is consumption net of subsistence, the argument of the log (`csub`).
"""
function wtp_fraction(sol, meas::Measurement)
    (; hh, gr, p) = sol
    x  = meas.χ * meas.wtp_score_sd
    cx = is_log_kernel(p.ψ) ? x : expm1((1 - p.ψ) * x) / (1 - p.ψ)
    ΔΛ = is_log_kernel(p.ψ) ? fill(p.λ * x, gr.Nz, gr.L) : cx .* (p.λ .+ (1 - p.ψ) .* hh.Λ)
    return (; ΔΛ, frac = -expm1.(-ΔΛ ./ p.μ), x)
end

"""
    school_quality_wtp(sol, meas)

BFM (2007) willingness to pay for a school whose average score is Δ =
`meas.wtp_score_sd` student SDs higher: a parent in l′ with consumption C pays
a = C(1 − e^{−ΔΛ/μ}) (`wtp_fraction`).

Returns `share = E[a]/E[C + p + c̲]` over parents (the target: \$19.70/month over mean
monthly consumption, which includes housing, so the residence cost p_{l′} of the
2026-10 extension and the subsistence c̲ = `csub` C̄ are added back to the denominator;
both are zero in the baseline) and
`share_l` by location. It also audits the consumption
support: `minC` over all (node, l′) pairs, the count of pairs at or below
`C_FLOOR` (`n_binding`), and two population shares at the floor: `floor_mass`, the
share that chooses a location with consumption at the floor although another location
is feasible, and `infeasible_mass`, the share with no location above the floor at its
policies (a goods cost paid everywhere, `csub` or `pres`, that exceeds what the poorest
earn). Binding pairs are usually unaffordable moves chosen with probability ≈ 0, which
smooth_log handles like log would; the two masses are what measure a real feasibility
problem, and `solve_theta` bounds their sum. Their integrands are choice indicators, so
they are interpolated linearly in the log shock (`k = 1`): a cubic spline through
nonnegative values can undershoot zero, and did report negative masses. The quadrature is
split at the grid nodes (`integrate_teach`), which a step-like integrand needs. Under the split
the pairs are (residence, workplace), a parent's WTP is by her residence (where the child goes
to school) and the audit counts the pairs.
"""
function school_quality_wtp(sol, meas::Measurement)
    sol.p.split && return school_quality_wtp_split(sol, meas)
    (; hh, Φ, Pi_T, Pi_O, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    (; frac, x) = wtp_fraction(sol, meas)             # a/C by (parent z, l′)
    Csum = zeros(L); Asum = zeros(L); Hsum = zeros(L)
    minC = Inf; n_binding = 0; floor_mass = 0.0; infeasible_mass = 0.0
    CT = zeros(nϵT, L); CO = zeros(nXO, L)
    FT = zeros(nϵT);    FO = zeros(nXO)
    XT = zeros(nϵT);    XO = zeros(nXO)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_T, Pi_O, p, g, l, zi)
        hb, d = br.h, br.d
        for k in 1:nϵT
            nd = teach_node(ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            XT[k] = all(b)                       # no feasible location at this node
            FT[k] = all(b) ? 0.0 : sum(br.πT[g, l, zi, k, b]; init = 0.0)
            CT[k, :] .= br.πT[g, l, zi, k, :] .* C
        end
        for k in 1:nXO
            nd = nonteach_node(XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sO[g, l, zi, k], hb.eO[g, l, zi, k], p)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            XO[k] = all(b)
            FO[k] = all(b) ? 0.0 : sum(br.πO[g, l, zi, k, b]; init = 0.0)
            CO[k, :] .= br.πO[g, l, zi, k, :] .* C
        end
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(hb, g, l, zi, gr)
        floor_mass += mass * sum(integrate_choice(FT, FO, cm, gr; k = 1))
        infeasible_mass += mass * sum(integrate_choice(XT, XO, cm, gr; k = 1))
        for lp in 1:L
            c = mass * sum(integrate_choice(CT[:, lp], CO[:, lp], cm, gr))
            Csum[lp] += c
            Asum[lp] += c * frac[zi, lp]
            h = p.pres[lp] + p.csub * p.Cbar
            iszero(h) || (Hsum[lp] += h * mass *
                          sum(integrate_choice(br.πT[g, l, zi, :, lp], br.πO[g, l, zi, :, lp], cm, gr)))
        end
    end
    return (; share = sum(Asum) / sum(Csum .+ Hsum), share_l = Asum ./ (Csum .+ Hsum), minC, n_binding,
              floor_mass = floor_mass / sum(Φ), infeasible_mass = infeasible_mass / sum(Φ), x)
end

"""
`school_quality_wtp` under the split: consumption and the floor audit run over the (residence, workplace) pairs.
Under `p.college` the cells are mixtures of the education branches, each with its own nodes and joint
probabilities (the residence marginals `Pi_Tr`, `Pi_Or` of the branch carry the housing costs).
"""
function school_quality_wtp_split(sol, meas::Measurement)
    (; hh, Φ, Pi_TW, Pi_OW, Pi_Tr, Pi_Or, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    (; frac, x) = wtp_fraction(sol, meas)             # a/C by (parent z, residence r)
    Csum = zeros(L); Asum = zeros(L); Hsum = zeros(L)
    minC = Inf; n_binding = 0; floor_mass = 0.0; infeasible_mass = 0.0
    CT = zeros(nϵT, L); CO = zeros(nXO, L)            # Σ_w π_{rw} C_{rw} by residence r
    FT = zeros(nϵT);    FO = zeros(nXO)
    XT = zeros(nϵT);    XO = zeros(nXO)
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, Pi_TW, Pi_OW, p, g, l, zi)
        hb, d = br.h, br.d
        for k in 1:nϵT
            nd = teach_node(ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p)
            π  = pair_probs(br.πT, g, l, zi, k)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            XT[k] = all(b)                       # no feasible pair at this node
            FT[k] = all(b) ? 0.0 : sum(π[b]; init = 0.0)
            CT[k, :] .= vec(sum(reshape(π .* C, L, L); dims = 1))
        end
        for k in 1:nXO
            nd = nonteach_node(XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d)
            C  = consumption(nd, hb.sO[g, l, zi, k], hb.eO[g, l, zi, k], p)
            π  = pair_probs(br.πO, g, l, zi, k)
            b  = C .≤ C_FLOOR
            minC = min(minC, minimum(C)); n_binding += count(b)
            XO[k] = all(b)
            FO[k] = all(b) ? 0.0 : sum(π[b]; init = 0.0)
            CO[k, :] .= vec(sum(reshape(π .* C, L, L); dims = 1))
        end
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(hb, g, l, zi, gr)
        floor_mass += mass * sum(integrate_choice(FT, FO, cm, gr; k = 1))
        infeasible_mass += mass * sum(integrate_choice(XT, XO, cm, gr; k = 1))
        πTr, πOr = branch_probs(Pi_Tr, p, d), branch_probs(Pi_Or, p, d)
        for r in 1:L
            c = mass * sum(integrate_choice(CT[:, r], CO[:, r], cm, gr))
            Csum[r] += c
            Asum[r] += c * frac[zi, r]
            h = p.pres[r] + p.csub * p.Cbar
            iszero(h) || (Hsum[r] += h * mass *
                          sum(integrate_choice(πTr[g, l, zi, :, r], πOr[g, l, zi, :, r], cm, gr)))
        end
    end
    return (; share = sum(Asum) / sum(Csum .+ Hsum), share_l = Asum ./ (Csum .+ Hsum), minC, n_binding,
              floor_mass = floor_mass / sum(Φ), infeasible_mass = infeasible_mass / sum(Φ), x)
end

"""
    wtp_by_parent_income(sol, ms, meas; nq = 3) -> (; elasticity, share, slope_z, elasticity_z)

Heterogeneity of the BFM willingness to pay over the adult pseudo-sample `ms`
(`micro_sample`), the Table 8 construct: adult i pays a = C (1 − e^{−ΔΛ(z, l′)/μ}) for the
school of `school_quality_wtp` (`wtp_fraction`), with C her consumption at the chosen
location (`ms.C`) and income y = exp(`logw`), pre-tax labour income with home-production
output included (the convention of `move_by_parent_income`).
- `elasticity`: the weighted OLS slope of a on y over all adults times ȳ/ā, the income
  elasticity of WTP. BFM's Table 8 has +\$1.38 (SE 0.33) per \$10,000 at mean income \$54,103
  and mean WTP \$19.69, 0.379 (SE ≈ 0.17 including the SE of mean WTP), conditional on education and children;
- `slope_z`, `elasticity_z`: the same construct conditional on ability, as BFM's is on
  education: the weighted OLS slope of a on y with a fixed effect for each z node of the
  adult (the within-node slope), and that slope times ȳ/ā, with the overall means of the
  unconditional version. Within a (z, l′) cell a is proportional to C, so under log utility
  this elasticity is about 1;
- `share[q]`: mean a over mean total consumption C + p_{l′} + c̲, the denominator of
  `school_quality_wtp`, within `nq`-tile q of y; slices straddling a cut are split
  uniformly, as in `move_by_parent_income`.
"""
function wtp_by_parent_income(sol, ms, meas::Measurement; nq = 3)
    (; p) = sol
    frac = wtp_fraction(sol, meas).frac
    y    = exp.(ms.logw)
    a    = ms.C .* [frac[ms.zi[k], ms.lp[k]] for k in eachindex(ms.w)]
    ȳ, ā = wmean(y, ms.w), wmean(a, ms.w)
    slope = sum(ms.w .* (y .- ȳ) .* (a .- ā)) / sum(ms.w .* (y .- ȳ) .^ 2)
    # z-node fixed effects: demean y and a within each ability node, then the same slope
    function demean_z(x)
        num = zeros(sol.gr.Nz); den = zeros(sol.gr.Nz)
        for k in eachindex(ms.w)
            num[ms.zi[k]] += ms.w[k] * x[k]
            den[ms.zi[k]] += ms.w[k]
        end
        return x .- [d > 0 ? n / d : 0.0 for (n, d) in zip(num, den)][ms.zi]
    end
    yz, az = demean_z(y), demean_z(a)
    slope_z = sum(ms.w .* yz .* az) / sum(ms.w .* yz .^ 2)
    Ctot  = ms.C .+ p.pres[ms.lp] .+ p.csub * p.Cbar
    cuts  = vcat(-Inf, collect(mquantile(ms, trues(length(ms.w)), [k / nq for k in 1:nq-1])), Inf)
    below(c) = clamp.((c .- ms.lo) ./ max.(ms.hi .- ms.lo, 1e-12), 0.0, 1.0)
    share = map(1:nq) do q
        wq = ms.w .* (below(cuts[q+1]) .- below(cuts[q]))
        sum(wq .* a) / sum(wq .* Ctot)
    end
    return (; elasticity = slope * ȳ / ā, share, slope_z, elasticity_z = slope_z * ȳ / ā)
end

"""
    quality_value(sol) -> (; dΛ, dΛ_equal, dΛ_quality)

The value of the suburb's school quality to a parent of ability z, Λ(z, 2) − Λ(z, 1) over
the z nodes, and the part of it that school quality alone accounts for. `dΛ` is the gap in
the solution's own Λ. `dΛ_equal` is the same gap with the altruism term re-evaluated at one
common quality, Q̄ = the population-weighted mean of Q: at fixed policies and occupation
choices h ∝ Q, so the human capital of every child born in l (`hT` and `capO`) is scaled by
Q̄/Q_l. `dΛ_quality = dΛ − dΛ_equal` is the part due to Q alone. `dΛ` also contains
what differs between the locations' children besides Q, chiefly the Ξ (moving-cost) offset
in their schooling (children born in different locations face different goods costs and
choose different (e, s)) and their occupational choices; that rest is `dΛ_equal`. With peer
effects Q includes the peer factor ā_l^{ψa} (`build_grids`), so "school quality" here is teachers
and peers together.
"""
function quality_value(sol)
    (; hh, Φ, gr, p) = sol
    w  = vec(sum(Φ; dims = 1))                         # population by location
    r  = reshape(dot(w, gr.Q) / sum(w) ./ gr.Q, 1, gr.L, 1, 1)   # Q̄/Q_l for the cells born in l
    scale(b) = merge(b, (; hT = b.hT .* r, capO = b.capO .* r))  # both education branches under college
    Λe = compute_lambda(p.college ? merge(hh, (; br1 = scale(hh.br1), br0 = scale(hh.br0))) : scale(hh), p, gr)
    dΛ, dΛe = hh.Λ[:, 2] .- hh.Λ[:, 1], Λe[:, 2] .- Λe[:, 1]
    return (; dΛ, dΛ_equal = dΛe, dΛ_quality = dΛ .- dΛe)
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
    cell_consumption(sol, g, l, zi) -> (CT, CO)

Consumption at the policies of cell (gender g, birth location l, ability node zi) on the
solver grids: `CT[k, l′]` at the teacher node k and `CO[k, l′]` at the non-teaching node k,
for work location l′. Under the split the last index runs over the (residence, workplace) pairs
d = (r − 1)L + w. Under `p.college` pass the education branch's arrays `hb` and its label `d`
(`cell_branches`).
"""
function cell_consumption(sol, g, l, zi; hb = sol.hh, d = nothing)
    sol.p.split && return cell_consumption_split(sol, g, l, zi; hb, d)
    (; gr, p) = sol
    CT = zeros(gr.nϵT, gr.L); CO = zeros(gr.nXO, gr.L)
    for k in 1:gr.nϵT
        nd = teach_node(gr.ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
        CT[k, :] .= consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p)
    end
    for k in 1:gr.nXO
        nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d)
        CO[k, :] .= consumption(nd, hb.sO[g, l, zi, k], hb.eO[g, l, zi, k], p)
    end
    return CT, CO
end

"`cell_consumption` under the split: `CT[k, d]`, `CO[k, d]` over the pairs d = (r − 1)L + w (`hb`, `d` as above: the education branch)."
function cell_consumption_split(sol, g, l, zi; hb = sol.hh, d = nothing)
    (; gr, p) = sol
    CT = zeros(gr.nϵT, gr.L^2); CO = zeros(gr.nXO, gr.L^2)
    for k in 1:gr.nϵT
        nd = teach_node(gr.ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
        CT[k, :] .= consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p)
    end
    for k in 1:gr.nXO
        nd = nonteach_node(gr.XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d)
        CO[k, :] .= consumption(nd, hb.sO[g, l, zi, k], hb.eO[g, l, zi, k], p)
    end
    return CT, CO
end

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
[lo, hi] (`mquantile`). `C` is consumption (`consumption`, after the moving, residence and
subsistence costs) at the point's work location l′, splined in the log shock from the solver
grid like `logw`. `groups` partitions the non-teaching occupations (default: one
group); `grp` is 0 for teachers and the group index otherwise. Weights sum to sum(Φ)
up to the shock mass outside the grids and the solver's quadrature error (≈1e-3).
Under the split every cell has one point set per (residence r, workplace w) pair, with the
joint probability π_{rw} in the weight: `lp` is the residence r (where the household
lives, pays tax and sends its child to school) and a new field `lw` the workplace w; the
teachers' log wage shifts by log[(1−τω_{T,g}) κ_w], the non-teachers' wages do not depend on w,
and `C` is consumption C_{rw} of the pair (after the commuting time).
Under `p.college` (§3c) the cells are mixtures of the education branches: a point's weight
carries P(d | g, l, z), its field `d` is its education (1 = BA; 1 for all without college),
a teacher's `logw` carries her pay factor (`teacher_logpay_shift`: ϑs (log s − log sBA) for
BA teachers, −ρT0 without a BA) and a BA non-teacher's `logw` is her capacity with the
credential premium e^{ρc}.
"""
function micro_sample(sol; n = 600, groups = nothing)
    sol.p.split && return micro_sample_split(sol; n, groups)
    (; hh, Φ, Pi_T, Pi_O, gr) = sol
    (; Nz, L, ϵTgrid, XOgrid, dT) = gr
    groups = something(groups, [gr.nonteach])
    G   = length(groups)
    out = (; g = Int[], l = Int[], zi = Int[], lp = Int[], grp = Int[],
             logw = Float64[], lo = Float64[], hi = Float64[], C = Float64[], w = Float64[], d = Int[])
    # a slice's log wage runs from e[j] to e[j+1] (edges), with midpoint value mid[j]
    function push_pts!(g, l, zi, lp, k, shift, mid, e, C, w, d)
        m = length(mid)
        append!(out.g, fill(g, m)); append!(out.l, fill(l, m)); append!(out.zi, fill(zi, m))
        append!(out.lp, fill(lp, m)); append!(out.grp, fill(k, m)); append!(out.w, w); append!(out.d, fill(d, m))
        append!(out.logw, shift .+ mid)
        append!(out.lo, shift .+ min.(view(e, 1:m), view(e, 2:m+1)))
        append!(out.hi, shift .+ max.(view(e, 1:m), view(e, 2:m+1)))
        append!(out.C, C)
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
        for l in 1:L, zi in 1:Nz, br in cell_branches(hh, Pi_T, Pi_O, sol.p, g, l, zi)
            mass = Φ[zi, l] / 2 * br.w
            iszero(mass) && continue
            hb, d = br.h, br.d
            cm   = choice_maps(hb, g, l, zi, gr)
            wT   = mass .* fT .* [teach_wt(cm, ϵ) for ϵ in ϵm]
            splT = Spline1D(luT, log.(hb.hT[g, l, zi, :]))
            lwT, lwTe = sol.p.γ .* splT(umT), sol.p.γ .* splT(uT)
            lwT  = lwT .+ teacher_logpay_shift(hb, g, l, zi, luT, umT, sol.p, d)
            lwTe = lwTe .+ teacher_logpay_shift(hb, g, l, zi, luT, uT, sol.p, d)
            wO   = mass .* fO .* [nonteach_wt(cm, x) for x in xm]
            splO = Spline1D(lx, log.(hb.capO[g, l, zi, :]))
            lcO, lcOe = splO(vm), splO(v)
            CT, CO = cell_consumption(sol, g, l, zi; hb, d)
            for lp in 1:L
                πT = clamp.(Spline1D(luT, br.πT[g, l, zi, :, lp]; k = 1)(umT), 0.0, 1.0)
                πO = clamp.(Spline1D(lx, br.πO[g, l, zi, :, lp]; k = 1)(vm), 0.0, 1.0)
                cT, cO = Spline1D(luT, CT[:, lp])(umT), Spline1D(lx, CO[:, lp])(vm)
                push_pts!(g, l, zi, lp, 0, log((1 - sol.p.τω[sol.p.T, g]) * sol.p.κ[lp]), lwT, lwTe, cT, wT .* πT, d)
                for k in 1:G
                    push_pts!(g, l, zi, lp, k, 0.0, lcO, lcOe, cO, wO .* πO .* view(pg, :, k), d)
                end
            end
        end
    end
    return out
end

"`micro_sample` under the split: one point set per (residence, workplace) pair, with the field `lw` (and `d`, as `micro_sample`)."
function micro_sample_split(sol; n = 600, groups = nothing)
    (; hh, Φ, Pi_TW, Pi_OW, gr) = sol
    (; Nz, L, ϵTgrid, XOgrid, dT) = gr
    groups = something(groups, [gr.nonteach])
    G   = length(groups)
    out = (; g = Int[], l = Int[], zi = Int[], lp = Int[], lw = Int[], grp = Int[],
             logw = Float64[], lo = Float64[], hi = Float64[], C = Float64[], w = Float64[], d = Int[])
    # a slice's log wage runs from e[j] to e[j+1] (edges), with midpoint value mid[j]
    function push_pts!(g, l, zi, r, w, k, shift, mid, e, C, wt, d)
        m = length(mid)
        append!(out.g, fill(g, m)); append!(out.l, fill(l, m)); append!(out.zi, fill(zi, m))
        append!(out.lp, fill(r, m)); append!(out.lw, fill(w, m)); append!(out.grp, fill(k, m)); append!(out.w, wt)
        append!(out.d, fill(d, m))
        append!(out.logw, shift .+ mid)
        append!(out.lo, shift .+ min.(view(e, 1:m), view(e, 2:m+1)))
        append!(out.hi, shift .+ max.(view(e, 1:m), view(e, 2:m+1)))
        append!(out.C, C)
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
        for l in 1:L, zi in 1:Nz, br in cell_branches(hh, Pi_TW, Pi_OW, sol.p, g, l, zi)
            mass = Φ[zi, l] / 2 * br.w
            iszero(mass) && continue
            hb, d = br.h, br.d
            cm   = choice_maps(hb, g, l, zi, gr)
            wT   = mass .* fT .* [teach_wt(cm, ϵ) for ϵ in ϵm]
            splT = Spline1D(luT, log.(hb.hT[g, l, zi, :]))
            lwT, lwTe = sol.p.γ .* splT(umT), sol.p.γ .* splT(uT)
            lwT  = lwT .+ teacher_logpay_shift(hb, g, l, zi, luT, umT, sol.p, d)
            lwTe = lwTe .+ teacher_logpay_shift(hb, g, l, zi, luT, uT, sol.p, d)
            wO   = mass .* fO .* [nonteach_wt(cm, x) for x in xm]
            splO = Spline1D(lx, log.(hb.capO[g, l, zi, :]))
            lcO, lcOe = splO(vm), splO(v)
            CT, CO = cell_consumption(sol, g, l, zi; hb, d)
            for r in 1:L, w in 1:L
                dest = (r - 1) * L + w
                πT = clamp.(Spline1D(luT, br.πT[g, l, zi, :, r, w]; k = 1)(umT), 0.0, 1.0)
                πO = clamp.(Spline1D(lx, br.πO[g, l, zi, :, r, w]; k = 1)(vm), 0.0, 1.0)
                cT, cO = Spline1D(luT, CT[:, dest])(umT), Spline1D(lx, CO[:, dest])(vm)
                push_pts!(g, l, zi, r, w, 0, log((1 - sol.p.τω[sol.p.T, g]) * sol.p.κ[w]), lwT, lwTe, cT, wT .* πT, d)
                for k in 1:G
                    push_pts!(g, l, zi, r, w, k, 0.0, lcO, lcOe, cO, wO .* πO .* view(pg, :, k), d)
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
    teacher_wage_stats(ms, weights; d = nothing) -> (; p9010, p9010_g)

Teachers' 90/10 of hourly (pre-tax) wages by gender (`p9010_g`, men and women, pooled
over locations as in the national ACS) and their weighted average with `weights`,
the Table 3 target. With `d` (`p.college`: 1 = BA, 0 = no BA) the teachers of that education only
(`teacher_p9010_ba`, which shows how much of the all-teacher 90/10 is the gap between the two tracks).
"""
function teacher_wage_stats(ms, weights; d = nothing)
    seld = d === nothing ? trues(length(ms.w)) : (ms.d .== d)
    r = [p9010(ms, (ms.grp .== 0) .& (ms.g .== g) .& seld) for g in 1:2]
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
    wls_slope(X, y, w; ba = nothing) -> Float64

Coefficient of the second column of `X` (the standardized score; the first is the constant) in the weighted
least-squares regression of `y` on the columns of `X` and, if given, the BA indicator `ba`. A BA indicator without
weighted variation (everyone or no one in the sample has a BA, as at a corner of the college choice, where the
weight of the other branch is zero) is collinear with the constant: it is dropped, and the slope is that of the
regression without it, the limit of the BA-conditional slope as the BA share tends to 0 or 1. NaN without weight
or if X′WX is singular.
"""
function wls_slope(X, y, w; ba = nothing)
    W = sum(w)
    W > 0 || return NaN
    if ba !== nothing
        sb = sum(w .* ba) / W
        sb * (1 - sb) > 1e-10 && (X = hcat(X, ba))
    end
    F = lu(X' * (w .* X); check = false)
    return issuccess(F) ? (F \ (X' * (w .* y)))[2] : NaN
end

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
- `rho_pc_all`: the same over parents of both genders; `c`, `sd_S` (adult scores);
- `wage_slope_ba` (`p.college`; NaN otherwise): the same regression with a BA indicator (the
  point's education `ms.d`), the model counterpart of the NLSY79 BA-conditional slope (`wls_slope`:
  without any weighted variation in the BA indicator it is the plain slope);
- `teacher_slope_ratio` (`p.college`; NaN otherwise): the slope among all teachers (with the BA
  indicator) over the slope among BA non-teacher market workers, the NLSY79 check of γ
  (closed form γ(1−η)/(1−γη) without selection on the education margin).
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
    if sol.p.college      # §3c: the slopes with the education branches (`wls_slope`)
        yk, dk = ms.logw[mkt], Float64.(ms.d[mkt] .== 1)
        selT, selC = ms.grp[mkt] .== 0, (ms.grp[mkt] .> 0) .& (dk .== 1)
        slope_ba = wls_slope(X, yk, wk; ba = dk)
        slope_T  = wls_slope(X[selT, :], yk[selT], wk[selT]; ba = dk[selT])
        slope_C  = wls_slope(X[selC, :], yk[selC], wk[selC])
    else
        slope_ba = slope_T = slope_C = NaN
    end

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
    # the score gap's composition part: the same gap with the c log Q term removed, in
    # the same SD units, so score_gap − score_gap_noQ is the school-quality part (teachers and,
    # with ψa ≠ 0, peers)
    EZ(l) = sum(Φ[:, l] .* log.(z)) / sum(Φ[:, l])
    return (; score_gap = (ES(2) - ES(1)) / wsd(Ss), score_gap_noQ = (EZ(2) - EZ(1)) / wsd(Ss),
              wage_slope = β[2], rho_pc = parent_child(ms.g .== 2), rho_pc_all = parent_child(trues(length(ms.w))),
              c, sd_S = sd, wage_slope_ba = slope_ba, teacher_slope_ratio = slope_T / slope_C)
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

"""
    parent_income_composition(sol, ms; nq = 3) -> (; city_share, c2s, s2c)

Residential composition and directional flows of children by `nq`-tile of the
parent's pre-tax income (home production included; quantiles over all adults), the
NLSY counterparts of T3 (nlsy_extensions.md, section B): `city_share[q]`, the share of
children schooled in location 1; `c2s[q]` and `s2c[q]`, the shares of children
schooled in 1 (in 2) who live in 2 (in 1) as adults. Children are schooled where the
parent lives and draw z′ ~ Πz, as in `move_by_parent_income`.
"""
function parent_income_composition(sol, ms; nq = 3)
    (; πbar, gr) = sol
    (; Πz, Nz) = gr
    gr.L == 2 || throw(ArgumentError("parent_income_composition is for two locations"))
    cuts  = vcat(-Inf, collect(mquantile(ms, trues(length(ms.w)), [k / nq for k in 1:nq-1])), Inf)
    below(c) = clamp.((c .- ms.lo) ./ max.(ms.hi .- ms.lo, 1e-12), 0.0, 1.0)
    out(l)   = [sum(Πz[zi, zj] * πbar[l, zj, 3 - l] for zj in 1:Nz) for zi in 1:Nz]
    p12, p21 = out(1), out(2)
    in1 = ms.lp .== 1
    mv  = [ms.lp[k] == 1 ? p12[ms.zi[k]] : p21[ms.zi[k]] for k in eachindex(ms.w)]
    city = zeros(nq); c2s = zeros(nq); s2c = zeros(nq)
    for q in 1:nq
        wq = ms.w .* (below(cuts[q+1]) .- below(cuts[q]))
        city[q] = sum(wq[in1]) / sum(wq)
        c2s[q]  = sum(wq[in1] .* mv[in1]) / sum(wq[in1])
        s2c[q]  = sum(wq[.!in1] .* mv[.!in1]) / sum(wq[.!in1])
    end
    return (; city_share = city, c2s, s2c)
end

"""
    own_income_suburb_share(ms, hpgrp; nq = 3) -> Vector

Share of market workers (teachers included, home producers excluded) living in
location 2, by `nq`-tile of their own pre-tax income: the NLSY adult suburb share by
own-earnings tercile (closest to its positive-earner variant, since home producers,
the NLSY's zero earners, are excluded here).
"""
function own_income_suburb_share(ms, hpgrp::Int; nq = 3)
    mkt  = market_workers(ms, hpgrp)
    cuts = vcat(-Inf, collect(mquantile(ms, mkt, [k / nq for k in 1:nq-1])), Inf)
    below(c) = clamp.((c .- ms.lo) ./ max.(ms.hi .- ms.lo, 1e-12), 0.0, 1.0)
    return map(1:nq) do q
        wq = ms.w .* (below(cuts[q+1]) .- below(cuts[q])) .* mkt
        sum(wq[ms.lp .== 2]) / sum(wq)
    end
end

"""
    own_z_composition(Φ; nq = 3) -> Vector

Share of the young schooled in location 1 by `nq`-tile of their own persistent
ability z: the NLSY79 origin city share by own-AFQT tercile. z is discrete, so each
node's mass is spread uniformly over its interval of the cumulative distribution and a
node that straddles a cut is split (on the Nz = 5 grid the middle tercile is the middle
node).
"""
function own_z_composition(Φ; nq = 3)
    w   = vec(sum(Φ; dims = 2))
    cum = cumsum(w) ./ sum(w)
    lo  = vcat(0.0, cum[1:end-1])
    return map(1:nq) do q
        f = [max(0.0, min(cum[i], q / nq) - max(lo[i], (q - 1) / nq)) / max(cum[i] - lo[i], 1e-300)
             for i in eachindex(w)]
        sum(f .* Φ[:, 1]) / sum(f .* w)
    end
end

"""
    move_by_own_z(sol; nq = 3) -> (; move, c2s, s2c)

Mobility of the young by `nq`-tile of their own z (2026-10-15, PREREG amendment 3; the terciles of
`own_z_composition`, so the counterpart of V6z's): `move[q]`, the share of the young schooled in one location who
live in the other as adults; `c2s[q]` and `s2c[q]`, the shares of those schooled in location 1 (in 2) who live in
2 (in 1). A young person of own z node i schooled in l (mass Φ[i, l]) lives in l′ with probability
πbar[l, i, l′]. NLSY79 counterparts (`settings/nlsy_composition.md` §2, by own-AFQT tercile): move rate 0.195,
0.184, 0.273; city → suburb 0.284, 0.378, 0.508; suburb → city 0.126, 0.118, 0.208. AFQT is a noisy score, so the
data's gradients are attenuated relative to these latent ones. Two locations only.
"""
function move_by_own_z(sol; nq = 3)
    (; Φ, πbar, gr) = sol
    gr.L == 2 || throw(ArgumentError("move_by_own_z is for two locations"))
    w   = vec(sum(Φ; dims = 2))
    cum = cumsum(w) ./ sum(w)
    lo  = vcat(0.0, cum[1:end-1])
    m12 = [Φ[i, 1] * πbar[1, i, 2] for i in eachindex(w)]
    m21 = [Φ[i, 2] * πbar[2, i, 1] for i in eachindex(w)]
    move = zeros(nq); c2s = zeros(nq); s2c = zeros(nq)
    for q in 1:nq
        f = [max(0.0, min(cum[i], q / nq) - max(lo[i], (q - 1) / nq)) / max(cum[i] - lo[i], 1e-300)
             for i in eachindex(w)]
        move[q] = sum(f .* (m12 .+ m21)) / sum(f .* w)
        c2s[q]  = sum(f .* m12) / sum(f .* Φ[:, 1])
        s2c[q]  = sum(f .* m21) / sum(f .* Φ[:, 2])
    end
    return (; move, c2s, s2c)
end

"""
    teacher_selection(sol, hp) -> (; z_gap_sd, z_gap_market_sd, T, O, H)

Mean log z of teachers minus that of non-teachers, in SDs of log z: against all
non-teachers (home production included) and against non-teacher market workers
(home producers, occupation `hp`, excluded). The NLSY79 counterpart is K–12 teachers'
latent AFQT against other workers at 25–34. `T`, `O`, `H` are the masses of teachers,
non-teachers and home producers by ability node and education branch, `[zi, d = 1 / d = 0]`
(the second column is empty without `p.college`), which `college_stats` uses.
"""
function teacher_selection(sol, hp::Int)
    (; hh, Φ, gr, p) = sol
    s_z = p.σξ / sqrt(1 - p.ρz^2)
    lz  = log.(gr.z)
    Tb = zeros(gr.Nz, 2); Ob = zeros(gr.Nz, 2); Hb = zeros(gr.Nz, 2)
    pHP = hp == 0 ? nothing : [[argmax_prob(gr, g, gr.XOgrid[k, g], hp) for k in 1:gr.nXO] for g in 1:2]
    for l in 1:gr.L, zi in 1:gr.Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        j  = br.d == 1 ? 1 : 2
        mT, mO = integrate_choice(ones(gr.nϵT), ones(gr.nXO), cm, gr)
        Tb[zi, j] += mass * mT; Ob[zi, j] += mass * mO
        pHP === nothing || (Hb[zi, j] += mass * integrate_nonteach(pHP[g], cm, gr))
    end
    T, O, H = vec(sum(Tb; dims = 2)), vec(sum(Ob; dims = 2)), vec(sum(Hb; dims = 2))
    mean_lz(x) = dot(x, lz) / sum(x)
    return (; z_gap_sd = (mean_lz(T) - mean_lz(O)) / s_z,
              z_gap_market_sd = hp == 0 ? NaN : (mean_lz(T) - mean_lz(O .- H)) / s_z,
              T = Tb, O = Ob, H = Hb)
end

"""
    composition_stats(sol, ms, hpgrp) -> NamedTuple

Residential composition behind the spatial gaps (calibration experiments, 2026-10):
- `dz_young_sd`, `dz_adult_sd`: gap in mean log z of the young (by schooling location)
  and of adults (by work location), in SDs of log z (stocks; inherited ability);
- `gap_meanlog_mkt`: gap in mean log pre-tax wage of market workers (teachers included),
  and `gap_meanlog_mkt_Qadj` the same with each worker's own schooling-Q effect removed (Q with
  its peer factor),
  log Q_birth/(1 − η) for non-teachers and γ log Q_birth/(1 − γη) for teachers (the
  destination-homogeneous elasticities): the composition part of the earnings gap;
- `mover_premium`: mean log wage of movers minus stayers among market workers born in l
  (flows; selection of movers);
- `teacher_premium[g]`: mean log wage of teachers minus non-teacher market workers of
  gender g; `gender_gap_market`: men minus women among non-teacher market workers.
"""
function composition_stats(sol, ms, hpgrp::Int)
    (; Φ, πbar, gr, p) = sol
    (; z, Q, Nz, L) = gr
    s_z = p.σξ / sqrt(1 - p.ρz^2)
    lz  = log.(z)
    young = [dot(Φ[:, l], lz) / sum(Φ[:, l]) for l in 1:L]
    A     = [sum(Φ[zi, l] * πbar[l, zi, lp] for l in 1:L) for zi in 1:Nz, lp in 1:L]
    adult = [dot(A[:, lp], lz) / sum(A[:, lp]) for lp in 1:L]
    mkt = market_workers(ms, hpgrp)
    k   = ifelse.(ms.grp .== 0, p.γ / (1 - p.γ * p.η), 1 / (1 - p.η))
    adj = ms.logw .- k .* log.(Q[ms.l])
    mean_at(x, sel) = wmean(x[sel], ms.w[sel])
    raw  = [mean_at(ms.logw, mkt .& (ms.lp .== lp)) for lp in 1:L]
    qadj = [mean_at(adj, mkt .& (ms.lp .== lp)) for lp in 1:L]
    mover = [mean_at(ms.logw, mkt .& (ms.l .== l) .& (ms.lp .!= l)) -
             mean_at(ms.logw, mkt .& (ms.l .== l) .& (ms.lp .== l)) for l in 1:L]
    nonT(g) = mkt .& (ms.grp .> 0) .& (ms.g .== g)
    prem = [mean_at(ms.logw, (ms.grp .== 0) .& (ms.g .== g)) - mean_at(ms.logw, nonT(g)) for g in 1:2]
    return (; dz_young_sd = (young[2] - young[1]) / s_z, dz_adult_sd = (adult[2] - adult[1]) / s_z,
              gap_meanlog_mkt = raw[2] - raw[1], gap_meanlog_mkt_Qadj = qadj[2] - qadj[1],
              mover_premium = mover, teacher_premium = prem,
              gender_gap_market = mean_at(ms.logw, nonT(1)) - mean_at(ms.logw, nonT(2)))
end

"""
    commute_shares(ms, hpgrp) -> (; all, teach, nonteach)

Shares of workers who work outside the location where they live (`ms.lp ≠ ms.lw`, the
residence/workplace split): among market workers (teachers included, home producers —
non-teaching group `hpgrp` — excluded), among teachers, and among non-teaching market workers.
The data counterparts are the shares crossing between the city and suburb classes of
POWPUMAs, 0.107 for all workers and 0.084 for K–12 teachers (`commuting_all.md`), a lower
bound that misses commuting within a county. NaN for the bundled model, whose points have no
workplace: the shares are zero there by construction, not predictions, and a NaN moment is left
out of the moment tables and the validation rows.
"""
function commute_shares(ms, hpgrp::Int)
    hasproperty(ms, :lw) || return (; all = NaN, teach = NaN, nonteach = NaN)
    mkt   = market_workers(ms, hpgrp)
    cross = ms.lp .!= ms.lw
    share(sel) = sum(ms.w[sel .& cross]) / sum(ms.w[sel])
    return (; all = share(mkt), teach = share(ms.grp .== 0), nonteach = share(mkt .& (ms.grp .> 0)))
end

"""
    commute_directions(ms, hpgrp) -> (; teach, nonteach)

The direction of commuting, as two-vectors (P(work in 2 | live in 1), P(work in 1 | live in 2)): the share of
the residents of the city (location 1) who work in the suburb and of the residents of the suburb who work in
the city, among teachers (`teach`) and among non-teaching market workers (`nonteach`; home producers, group
`hpgrp`, are out as in `commute_shares`). The ACS 2018 counterparts (`commuting_all.md`, city and suburb
classes of POWPUMAs) are, as conditionals within the closed city–suburb block like the model's (and like the
crossing-share target 0.107), 0.0862, 0.0832 for K–12 teachers and 0.0919, 0.1164 for non-teachers (as row
conditionals whose denominators also hold the 1–2% who work in neither class: 0.0852, 0.0819 and 0.0911, 0.1148).
NaN for the bundled model, which has no workplace.
"""
function commute_directions(ms, hpgrp::Int)
    hasproperty(ms, :lw) || return (; teach = fill(NaN, 2), nonteach = fill(NaN, 2))
    mkt = market_workers(ms, hpgrp)
    function direction(sel)
        from1, from2 = sel .& (ms.lp .== 1), sel .& (ms.lp .== 2)
        return [sum(ms.w[from1 .& (ms.lw .== 2)]) / sum(ms.w[from1]),
                sum(ms.w[from2 .& (ms.lw .== 1)]) / sum(ms.w[from2])]
    end
    return (; teach = direction(ms.grp .== 0), nonteach = direction(mkt .& (ms.grp .> 0)))
end

"""
    option_values(sol) -> (; teach, nonteach, gap)

Mean over the teachers and over the non-teachers (home producers included) of the inclusive-value
gap relative to the own-location value, Σ_r π_r (I_r − W_rr) in utils: the inner value I_r of the
position-weighted nest less the value W_rr of working where one lives; `gap` is `teach − nonteach`.
This is not the value of the workplace choice. I_r − W_rr = σw log Σ_w P_w e^{(W_rw − W_rr)/σw} is at
most zero when commuting costs and pay is equal, and exactly zero when W_rw = W_rr for every w (τk = 0
at equal pay) although about half of the workers cross: it contains the position normalisation
(the weights P_w need not favour the own location) as well as the commuting penalty. Every worker
is in the nest, so the two means are close, and only their difference, `gap`, bears on the
occupational choice: a positive gap tilts it toward teaching. All zero for the bundled model. Under
`p.college` the means are over both education branches, each at its own nodes, with the weights P(d | g, l, z).
"""
function option_values(sol)
    sol.p.split || return (; teach = 0.0, nonteach = 0.0, gap = 0.0)
    (; hh, Φ, gr, p) = sol
    (; Nz, L, nϵT, nXO, ϵTgrid, XOgrid) = gr
    diag_d = [(r - 1) * L + r for r in 1:L]                # the pairs (r, r)
    function value(nd, s, e)
        W = location_values(nd, consumption(nd, s, e, p), p)
        _, _, πr, _, I = nested_logsum_probs(W, nd.l, nd, p)
        return dot(πr, I .- W[diag_d])
    end
    oT = zeros(nϵT); oO = zeros(nXO)
    numT = numO = denT = denO = 0.0
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        hb, d = br.h, br.d
        for k in 1:nϵT
            oT[k] = value(teach_node(ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d), hb.sT[g, l, zi, k], hb.eT[g, l, zi, k])
        end
        for k in 1:nXO
            oO[k] = value(nonteach_node(XOgrid[k, g], zi, l, g, hb.Λ, p, gr; d), hb.sO[g, l, zi, k], hb.eO[g, l, zi, k])
        end
        cm = choice_maps(hb, g, l, zi, gr)
        IT, IO = integrate_choice(oT, oO, cm, gr)
        mT, mO = integrate_choice(ones(nϵT), ones(nXO), cm, gr)
        numT += mass * IT; denT += mass * mT
        numO += mass * IO; denO += mass * mO
    end
    return (; teach = numT / denT, nonteach = numO / denO, gap = numT / denT - numO / denO)
end

"Step of log κ₂ in `teacher_supply_elasticity`."
const SUPPLY_STEP = 0.01

"""
    teacher_work_counts(hh, Φ, Pi_TW, gr, p = Params()) -> Vector

Mass of teachers by workplace w under the split, Σ_{g,l,z} Φ[z, l]/2 ∫ teach(ϵ_T) Σ_r π_{rw}: the
`nT` of `location_stats` summed over gender, for the household solution `hh` and its joint
probabilities `Pi_TW`, with the stationary distribution `Φ` as it is. Under `p.college` (pass `p`)
`hh` and `Pi_TW` are the two education branches, mixed with the weights P(d | g, l, z).
"""
function teacher_work_counts(hh, Φ, Pi_TW, gr, p::Params = Params())
    (; Nz, L) = gr
    πw = workplace_probs(Pi_TW)                       # [g, l, zi, k, w]; a pair of them under college
    n  = zeros(L)
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, πw, πw, p, g, l, zi)
            mass = mass0 * br.w
            cm = choice_maps(br.h, g, l, zi, gr)
            for w in 1:L
                n[w] += mass * integrate_teach(br.πT[g, l, zi, :, w], cm, gr)
            end
        end
    end
    return n
end

"""
    teacher_supply_elasticity(sol; step = SUPPLY_STEP, hh_tol = 1e-7, hh_maxit = 1000) -> Float64

d log(n_work,2/n_work,1) / d log(κ₂/κ₁) under the split: the response of the teachers' workplace
shares, location 2 against location 1, to the suburb's pay premium, by a one-sided household re-solve
at κ₂ e^{step} with κ₁ unchanged, holding (H̃_T, M, t), the grids (Q, the position weights) and the
stationary distribution Φ at the solution's: partial equilibrium, the response of occupation,
residence and workplace choices at given aggregates (under `p.college` also the education choice,
which the re-solve updates). n_work is the mass of teachers working at w
(`teacher_work_counts`). The literature's own-wage elasticities of teacher supply are Ransom–Sims
(2010) 3.7 and Falch (2010) ≈ 1.4; both are single-district estimates, so they are upper bounds for
a location whose neighbours' pay moves with it (the inactive target `teacher_supply_elasticity`, which
with `--free ρw` can identify ρw instead of the crossing share). NaN for the bundled model, where
workers do not choose a workplace. The ratio of `teacher_supply_response`, which also returns its two parts.
"""
teacher_supply_elasticity(sol; step = SUPPLY_STEP, hh_tol = 1e-7, hh_maxit = 1000) =
    teacher_supply_response(sol; step, hh_tol, hh_maxit).ratio

"""
    teacher_supply_response(sol; step = SUPPLY_STEP, hh_tol = 1e-7, hh_maxit = 1000) -> (; ratio, own, cross)

The re-solve of `teacher_supply_elasticity` (κ₂ e^{step}; Φ, H̃_T, M, t and the grids fixed) and the three
responses of the teachers' workplace masses n_w it gives (2026-10-15): `ratio` = d log(n₂/n₁)/d log κ₂ (the
moment `teacher_supply_elasticity`), `own` = ε₂₂ = d log n₂/d log κ₂, the suburb's own-location elasticity at
household partial equilibrium (PREREG A6: ≤ 3.7, Ransom–Sims 2010), and `cross` = ε₁₂ = d log n₁/d log κ₂, so
that ratio = own − cross. NaNs for the bundled model. `spatial_supply_audit.jl` decomposes them by margin.
"""
function teacher_supply_response(sol; step = SUPPLY_STEP, hh_tol = 1e-7, hh_maxit = 1000)
    sol.p.split || return (; ratio = NaN, own = NaN, cross = NaN)
    (; hh, Φ, gr, p) = sol
    p2  = with_params(p; κ = [p.κ[1], p.κ[2] * exp(step)])
    hh2 = solve_household(p2, gr; tol = hh_tol, maxit = hh_maxit, hh0 = hh)
    n1  = teacher_work_counts(hh, Φ, sol.Pi_TW, gr, p)
    n2  = teacher_work_counts(hh2, Φ, location_probs(hh2, p2, gr)[1], gr, p2)
    return (; ratio = (log(n2[2] / n2[1]) - log(n1[2] / n1[1])) / step,
              own = (log(n2[2]) - log(n1[2])) / step, cross = (log(n2[1]) - log(n1[1])) / step)
end

"""
    teacher_district_elasticity(sol) -> Vector

The own-pay supply elasticity of teachers to a school district of negligible weight inside workplace w, for
each w. It is the teacher-mass-weighted mean, over teacher nodes and residences r, of
    μ y_{rw} / (σ_in C_{rw}),
with y_{rw} the node's income at the pair (r, w) before goods costs (what pay scales: C_{rw} = y_{rw} −
(1 + τe) e − m_r) and C_{rw} its consumption, both at the solution's policies; σ_in is the workplace taste
scale σw under the split and the location taste scale σν in the bundled model, where the destination is
the location. NaN for a workplace without teachers.

Derivation. A district j inside w raises its pay by a proportion x. Its pay scales the income of a worker
who lives in r, ∂ log y_{rj}/∂ log κ_j = 1, so ∂W_{rj}/∂ log κ_j = μ y_{rj}/C_{rj}, with W = μ log C + b_w + Λ_r
(− κc 1{r ≠ w}) as in `location_values`: the friction, the taste and the amenities are constants of the node.
The policies (e, s) maximise the node's inclusive value, which one district of negligible weight moves by a
negligible amount, so they do not respond (the envelope theorem); nor do occupation, the aggregates (H̃_T, M,
t) or the stationary distribution. Each position has its own iid taste shock of scale σ_in, so given r the
probability of choosing a position of j is π_{j|r} ∝ P_j exp(W_{rj}/σ_in) and
    ∂ log π_{rj}/∂ log κ_j = μ y_{rj} / (σ_in C_{rj}) (1 − π_{j|r}),
the response of the residence choice π_r through the inclusive value I_r being of the order of π_{j|r}. For a
district of negligible weight π_{j|r} → 0, so the factor (1 − π_{j|r}) drops. The district's teachers come from
the pairs (r, w) of the nodes in proportion to the mass π^T_{rw} the nodes put on w (at equal pay the district
is a typical position of w), so the elasticity of its supply is the average of the pairs' elasticities with
those masses as weights. It is the elasticity of the supply to one district, not the response of a workplace's
teacher mass to its own pay (`teacher_supply_elasticity`, a ratio of two workplaces' masses, in which a
rise in κ_w also draws from the other workplace through the workplace nest and the district's neighbours rise
with it): PREREG A6's bar for a school district is 10 (for the own-location elasticity ε₂₂ at household
partial equilibrium, `teacher_supply_response`, it is 3.7: Ransom–Sims 2010). Under `p.college` the cells mix
the two education branches with P(d | g, l, z), each at its own nodes. Pairs at the consumption floor
(C ≤ C_FLOOR) are dropped from the weights and from the integrand alike: they carry no choice mass unless the
whole node is infeasible, and their integrand, of the order of 1/C_FLOOR, would leak into neighbouring nodes
through the cubic spline of `integrate_teach`.
"""
function teacher_district_elasticity(sol)
    (; hh, Φ, gr, p) = sol
    (; Nz, L, nϵT, ϵTgrid) = gr
    σin = p.split ? p.σw : p.σν
    PT, PO = p.split ? (sol.Pi_TW, sol.Pi_OW) : (sol.Pi_T, sol.Pi_O)
    pairs_rw = p.split ? [(r, w) for r in 1:L for w in 1:L] : [(w, w) for w in 1:L]
    num = zeros(L); den = zeros(L)
    ρ  = zeros(nϵT, L, L)                                # the elasticity of node k at the pair (r, w)
    ok = zeros(nϵT, L, L)                                # 1 where that pair is above the consumption floor, else 0
    for l in 1:L, zi in 1:Nz, g in 1:2, br in cell_branches(hh, PT, PO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        hb, d = br.h, br.d
        cm = choice_maps(hb, g, l, zi, gr)
        for k in 1:nϵT
            nd = teach_node(ϵTgrid[k], zi, l, g, hb.Λ, p, gr; d)
            e  = hb.eT[g, l, zi, k]
            C  = consumption(nd, hb.sT[g, l, zi, k], e, p)
            x  = p.μ .* (C .+ (1 + nd.τe) * e .+ nd.m) ./ (σin .* max.(C, C_FLOOR))   # over the destinations
            f  = Float64.(C .> C_FLOOR)                  # floored pairs drop out
            if p.split
                ρ[k, :, :]  .= reshape(x .* f, L, L)'    # [w, r] → [r, w]
                ok[k, :, :] .= reshape(f, L, L)'
            else
                for w in 1:L
                    ρ[k, w, w], ok[k, w, w] = x[w] * f[w], f[w]
                end
            end
        end
        for (r, w) in pairs_rw
            πT = p.split ? br.πT[g, l, zi, :, r, w] : br.πT[g, l, zi, :, w]
            den[w] += mass * integrate_teach(πT .* ok[:, r, w], cm, gr)
            num[w] += mass * integrate_teach(πT .* ρ[:, r, w], cm, gr)
        end
    end
    return num ./ den
end

# -----------------------------------------------------------------------------
# 6. All moments
# -----------------------------------------------------------------------------
"""
    model_moments(sol, meas = Measurement(); hh_tol = 1e-7, hh_maxit = 1000) -> NamedTuple

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
- `payroll_pp_gap`: Δlog teacher payroll per student, log(pay₂/S₂) − log(pay₁/S₁), which equals
  `gap_salary + gap_teachers_pp` (district data −0.0091, an inactive target);
- `seda_growth_gap`, `gap_logQ`: the SEDA growth mapping of §3 and the quality gap;
- `gap_logQ_teach`, `gap_logQ_peer`, `psi_a`: the two parts of `gap_logQ`, the teacher-quality index
  Δlog Q_T and the peer factor ψa Δlog ā (`build_grids`; zero without peer effects, ψa = 0);
  `log_abar`: the students' log ā_l at the solution (`peer_log_ability`, also without peers; its
  gap is s_z times `dz_young_sd`) and `peer_resid`, its largest distance from the log ā_l the grid
  was built at (the convergence of the peer fixed point; 0 without peers);
- `score_gap`: SEDA score-level gap in cohort SDs of log z + c log Q (`ability_moments`);
- `poverty_gap`, `poverty_share_l`: Δ share of children below the `meas.pov_rate`
  income line, and the shares by location (`poverty_gap`);
- `wage_slope`, `rho_pc`, `score_loading`: the NLSY moments in the full equilibrium
  (§4 consistency check; `ability_moments`);
- `teach_share`: for the occupational block;
- `move_rate`, `p12`, `p21`: gross and directional birth → adult location flows (the
  move rate is the Table 3 residential-transition target);
- `move_rate_by_parent_income`: the move rate by parental-income tercile (validation);
- `city_share_by_parent_income`, `city_share_gradient` (top minus bottom tercile),
  `c2s_by_parent_income`, `s2c_by_parent_income`: children's location and directional
  flows by parental-income tercile (`parent_income_composition`, NLSY T3);
  `suburb_share_by_own_income`: adults (market workers) by own-income tercile;
  `city_share_by_own_z`, `city_share_gradient_z`: the young by own-z tercile
  (`own_z_composition`, NLSY own-AFQT counterpart);
- `move_rate_by_own_z`, `c2s_by_own_z`, `s2c_by_own_z` (2026-10-15, validation; PREREG amendment 3): the young's
  move rate and directional moves by own-z tercile (`move_by_own_z`; NLSY79 0.195, 0.184, 0.273; 0.284, 0.378,
  0.508; 0.126, 0.118, 0.208);
- `teacher_z_gap_sd`, `teacher_z_gap_market_sd`: teachers' mean log z against non-teachers
  (`teacher_selection`); `gap_logh_teacher`: Δ mean log h_T of teachers by work location;
- `score_gap_noQ`, `dz_young_sd`, `dz_adult_sd`, `gap_meanlog_mkt`, `gap_meanlog_mkt_Qadj`,
  `mover_premium`: composition behind the gaps (`composition_stats`, `ability_moments`);
- `teacher_premium_m`, `teacher_premium_f`, `gender_gap_market`: teachers' and women's
  mean log wages against non-teacher market workers (`composition_stats`);
- `teacher_sd`, `cfr_nonteach`, `cfr_analytic`: pieces of the CFR moment;
- `wtp_share_l`: WTP share by location;
- `s_ratio`: mean s_T / mean s_O (the §1 schooling diagnostic, T1);
- `students_per_teacher`, `t`, `Q`, `QT` (its teacher part), `M` (model levels without a data counterpart);
- `mean_C`, `mean_h`: drift from the frozen C̄ and h̄ (`mean_C` is net of the subsistence `csub` C̄;
  model units only);
- `stationarity_dev`: max |adults_l − students_l|/students_l (quadrature accuracy);
- `minC`, `n_binding`, `floor_mass`, `infeasible_mass`: consumption support (see
  `school_quality_wtp`): the population shares that choose a floored location although a feasible
  one exists, and that have no feasible location.

Model levels, reported alongside the calibration targets (2026-10):
- `pop_share`: M_l/ΣM, the population shares of city and suburb (data: 0.579 of the
  partitioned pupils are in the suburb, `district_settings.md`);
- `t_market`, `t_market_l`: teacher wage bill over market income (home-production output
  out), economy-wide and by location, while the solved `t` is on all income
  (`home_production_stats`; NaN without `meas.hp`). Its data comparator is ≈ 0.016, teacher
  salaries over market labour income; the §1 t ≈ 0.023 includes benefits;
- `hp_share_m`, `hp_share_f`: home-production shares of men and women (NaN without `meas.hp`);
- `s_O_mean`, `s_T_mean`: mean time shares of non-teachers and teachers (`s_summary`);
- `dLambda_z`, `dLambda_mean`: Λ(z, 2) − Λ(z, 1), the utility value of the suburb's
  school quality to a parent of ability z over the z nodes, and its Φ-weighted mean. It also
  contains the Ξ (moving-cost) offset in children's schooling, not only school quality;
- `dLambda_quality`, `dLambda_quality_mean`: the part of it due to Q alone, the gap minus the
  same gap with every child's h′ rescaled to the population-weighted mean Q (`quality_value`);
- `wtp_income_elasticity`, `wtp_share_by_parent_income`: the BFM Table 8 income
  elasticity of WTP (data 0.379) and the WTP share by parental-income tercile
  (`wtp_by_parent_income`);
- `wtp_income_slope_z`, `wtp_income_elasticity_z`: the within-ability slope of WTP on income, with a
  fixed effect for each z node as BFM's slope holds education fixed, and its elasticity
  (about 1 under log utility);
- `teacher_supply_elasticity`: the split's partial-equilibrium supply elasticity of teachers' workplace
  shares, reported with these levels and against the literature's 3.7 (Ransom–Sims) and 1.4 (Falch);
  defined below with the residence/workplace moments, NaN when bundled.

Residence and workplace (`Params.split`). In the bundled model a worker works where she lives, and
the moments that exist only with a workplace are NaN: the zeros they would be are not predictions, and a
NaN moment is left out of the moment tables and the validation rows. Under the split `gap_salary`,
`gap_teachers_pp`, `students_per_teacher`, `gap_logh_teacher`, `teacher_sd` and `payroll_pp_gap` are by
workplace (teachers per student residing there), `gap_income`, the poverty and composition moments and
`stationarity_dev` by residence:
- `commute_share_all`, `commute_share_teach`, `commute_share_nonteach`: shares of market workers,
  of teachers and of non-teaching market workers with r ≠ w (`commute_shares`; data 0.107 for all
  workers and 0.084 for K–12 teachers, ACS city/suburb POWPUMA crossing shares); NaN if bundled;
- `n_teach_res`, `n_teach_work`: masses of teachers by residence and by workplace;
- `gap_teach_res_emp`: Δlog (teachers residing per market worker residing), data +0.298 (SE 0.041) for
  public K–12 teachers per employed resident (teacher_residence_targets.toml, [k12]); defined in both
  models;
- `net_commute_gap`: Δlog n_res − Δlog n_work of teachers, positive when they live in the suburb more
  than they work there (data +0.075, SE 0.038, [k12] of teacher_residence_targets.toml; the national
  PUMS cross-check is +0.108 with sampling SE 0.008); NaN if bundled;
- `teacher_supply_elasticity`: d log(n_work,2/n_work,1)/d log(κ₂/κ₁), by a one-sided household re-solve at
  κ₂ e^{0.01} holding (H̃_T, M, t) fixed (`teacher_supply_elasticity`; Ransom–Sims 2010: 3.7, Falch 2010:
  ≈ 1.4, single-district estimates and so upper bounds for a location); NaN if bundled;
- `option_value_teach`, `option_value_nonteach`, `option_value_gap`: the inclusive-value gap relative to
  the own-location value, I_r − W_rr, averaged over teachers and over non-teachers, and their difference
  (`option_values`; not the value of the workplace choice); zero if bundled;
- `commute_dir_teach`, `commute_dir_nonteach` (2026-10-15, validation): (P(work in 2 | live in 1), P(work in 1 | live
  in 2)) among teachers and among non-teaching market workers (`commute_directions`; ACS 2018, within the
  city–suburb block: 0.0862, 0.0832 and 0.0919, 0.1164); NaN vectors if bundled;
- `teacher_district_elasticity` (2026-10-15, diagnostic, both models): by workplace, the own-pay supply
  elasticity of a district of negligible weight inside it, the teacher-mass-weighted mean of
  μ y_{rw}/(σ_in C_{rw}) (`teacher_district_elasticity`: σ_in = σw under the split, σν bundled; the
  derivation is in its docstring). The audit bar is 10 (PREREG A6); `teacher_supply_elasticity` is the ratio
  response of the two workplaces' masses to a rise in κ₂, which is not it;
- `teacher_own_elasticity`, `teacher_cross_elasticity` (2026-10-15, diagnostic): the two parts of
  `teacher_supply_elasticity` from the same re-solve (`teacher_supply_response`), ε₂₂ = d log n_work,2/d log κ₂,
  the own-location elasticity at household partial equilibrium (PREREG A6: ≤ 3.7, Ransom–Sims 2010), and
  ε₁₂ = d log n_work,1/d log κ₂; NaN if bundled.

Education (§3c, `Params.college`; NaN without it; `college_stats`, `ability_moments`):
- `ba_share`, `ba_share_nonteach`, `ba_share_teach`: BA shares of adults, of non-teachers (home
  production included) and of teachers;
- `ba_z_gap_sd`: mean log z of BA minus non-BA adults (SDs of log z); `ba_share_by_z`: P(BA) by
  tercile of own z, `ba_share_z_gap` its top minus bottom; `ba_share_by_score`, `ba_share_score_gap`
  the same by tercile of the noisy NLSY79-like score (`ba_share_by_score_terciles`);
- `ba_premium`, `ba_premium_raw`: mean log wage of BA minus non-BA non-teacher market workers, at
  fixed gender (the data's within-occupation premium) and unadjusted;
- `teacher_premium_ba_m`, `teacher_premium_ba_f`: teachers minus BA non-teacher market workers;
  `teacher_premium_within_ba_m`, `_f`: BA teachers minus the same;
- `teacher_z_gap_ba_sd`, `teacher_z_gap_within_ba_sd`: mean log z of all teachers, of BA teachers,
  minus BA non-teacher market workers (SDs of log z);
- `s_T_mean_ba`, `s_O_mean_ba`, `s_O_mean_nonba`, `s_T_mean_nonba`: mean time share of BA teachers, of BA and
  non-BA non-teachers and of non-BA teachers; `ma_premium_implied`: ϑs log(17/16), the salary schedule's
  master's premium in the units of the `s_T_mean_ba` target (`ma_premium_implied_175`: with 17.5 years);
  `slope_ratio_closed_form`: γ(1−η)/(1−γη), the `teacher_slope_ratio` without selection effects;
- `teacher_ba_wage_gap`, `_m`, `_f`: mean log wage of BA minus non-BA teachers (the data counterpart of
  the non-BA track's pay penalty ρT0); `teacher_p9010_ba`, `teacher_p9010_ba_g`: the teachers' 90/10 of BA
  teachers alone, beside `teacher_p9010` (all teachers);
- `wage_slope_ba`: `wage_slope` with a BA indicator; `teacher_slope_ratio`: teachers' wage–score slope
  (with the BA indicator) over that of BA non-teacher market workers (`ability_moments`);
  `pooled_p9010`: the pooled within-occupation 90/10 of non-teachers' wages (`meas.pooled_p9010` only;
  NaN otherwise);
- `ba_share_by_residence`, `ba_share_gap` (2026-10-15, validation): the BA share of adults by residence
  and its gap, suburb minus city, in levels (`ba_share_by_residence`); the counterpart of the sealed
  hold-out B1 (the data: the share of adults 25+ with a BA or more, suburb minus city). NaN without `college`.

Round 2026-10-15c:
- `teacher_shadow_b`: teachers' value of the suburban workplace, BT[2] − BT[1] (utils): the shadow price of the
  teacher-demand closure (`Params.fte_rule`) when it is on, else the parameter b_T; NaN when bundled;
- `mean_log_yhat`: population mean of the predicted log income ln ŷ, the grading variable of the income-graded
  taste `Params.By` (`mean_predicted_log_income`); the calibration freezes `lyref` at its value at a reference fit.

Round 2026-10-16 (diagnostics, unscored; appended after `mean_log_yhat`; none is a target, and the moments above are
unchanged by them). Gaps are suburb minus city (index 2 minus 1); `_res` is by residence r and `_work` by workplace w
(identical in the bundled model, r = w = the work location). They come from `teacher_sorting_stats` and
`sorting_moments`, which define them:
- `cfr_rho`, `cfr_sd_eff`: `Measurement.cfr_rho` and the SD of teachers' log effectiveness the CFR experiment uses,
  sd(log h_eff | w)/ρ (`teacher_intervention`; `teacher_sd` stays sd(log h_eff)); `lyT_fixed`: `Params.lyT_fixed` (NaN: off);
- teachers' sorting: `gap_logh_teacher_res/_work` (Δ mean log h_T; the work version is `gap_logh_teacher`),
  `gap_logheff_teacher_res/_work` (log h_eff), `gap_Heff_teacher_res/_work` (Δ log E[h_eff^{β/σ}]),
  `gap_logpay_teacher_res/_work` (mean log pre-tax pay), `gap_salary_res` (Δ log mean pay by residence),
  `gap_lyhat_teacher_res/_work`, `mean_log_yhat_teach`, `mean_log_yhat_mkt`, `mean_log_yhat_nonteach_mkt` (ln ŷ);
- the decomposition of Δ mean log h_T into log Q_l, α log z, α log ϵ, φ log s and η log e:
  `gap_logh_teacher_{res,work}_{Q,z,eps,s,e}`, with `decomp_resid_max` (quadrature residual) and
  `share_born_suburb_res/_work` (share of the teachers born in location 2);
- the teacher-quality index: `gap_logQ_teach_quality` + `gap_logQ_teach_staffing` + `gap_logQ_teach_check` = `gap_logQ_teach`;
- the (r, w) cells, column-major: `teacher_rw_mass`, `teacher_rw_logh`, `teacher_rw_logpay`, and the contrasts at a given
  workplace `gap_logpay_teacher_res_fixedw`, `gap_logh_teacher_res_fixedw` (NaN if bundled);
- non-teaching market workers by residence (home producers out): `gap_lyhat_mkt_nonteach_res`, `gap_logy_mkt_nonteach_res`;
- levels: `teacher_logh_res/_work`, `teacher_logpay_res/_work`, `teacher_lyhat_res/_work`, `teacher_Heff_res/_work`,
  `nonteach_mkt_logy_res`, `nonteach_mkt_lyhat_res`, `nonteach_mkt_mass_res`.
"""
function model_moments(sol, meas::Measurement = Measurement(); hh_tol = 1e-7, hh_maxit = 1000)
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
    pic = parent_income_composition(sol, ms)
    cmp = composition_stats(sol, ms, hpg)
    zc  = own_z_composition(Φ)
    mz  = gr.L == 2 ? move_by_own_z(sol) : (; move = fill(NaN, 3), c2s = fill(NaN, 3), s2c = fill(NaN, 3))
    tsel = teacher_selection(sol, meas.hp)
    hps = meas.hp == 0 ? nothing : home_production_stats(sol, meas.hp; stats = st)
    wti = wtp_by_parent_income(sol, ms, meas)
    cs  = commute_shares(ms, hpg)
    ov  = option_values(sol)
    tsr = teacher_supply_response(sol; hh_tol, hh_maxit)
    cdr = commute_directions(ms, hpg)
    tde = teacher_district_elasticity(sol)
    bsr = ba_share_by_residence(sol)

    nTl    = vec(sum(st.nT; dims = 1))
    adults = nTl .+ vec(sum(st.nO; dims = 1))
    nTr    = vec(sum(st.nT_res; dims = 1))                   # teachers by residence (nTl when bundled)
    p.split && (adults = nTr .+ vec(sum(st.nO; dims = 1)))   # the residents, whose number is the students'
    meas.income in keys(inc) ||
        throw(ArgumentError("income must be one of $(keys(inc)), got :$(meas.income)"))
    (meas.hp == 0 && endswith(string(meas.income), "_market")) &&
        throw(ArgumentError("income = :$(meas.income) needs the home-production index `hp`"))
    income = inc[meas.income]

    S  = st.S
    M  = 2 .* S
    P  = [sum(Φ[zi, l] * πbar[l, zi, lp] for zi in 1:Nz) for l in 1:L, lp in 1:L]
    P ./= sum(P; dims = 2)                      # rows integrate to 1 + O(quadrature error)
    gap_logQ = gap(gr.Q)
    la = peer_log_ability(Φ, gr)                # the students' log ā_l, the peer input at this Φ
    ss = s_summary(sol)
    dΛ = sol.hh.Λ[:, 2] .- sol.hh.Λ[:, 1]       # value of the suburb's school quality, by z
    qv = quality_value(sol)
    edu = college_stats(sol, ms, ss, tsel, hpg; teacher_weights = meas.teacher_weights)    # `cs` is the commuting shares
    tsm = sorting_moments(teacher_sorting_stats(sol, meas), sol)                          # round 2026-10-16 diagnostics

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
        gap_logQ_teach = gap(gr.QT), gap_logQ_peer = p.ψa * (gr.la[2] - gr.la[1]), psi_a = p.ψa,
        log_abar = la, peer_resid = iszero(p.ψa) ? 0.0 : maximum(abs, la .- gr.la),
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
        city_share_by_parent_income = pic.city_share,
        city_share_gradient = pic.city_share[end] - pic.city_share[1],
        c2s_by_parent_income = pic.c2s, s2c_by_parent_income = pic.s2c,
        suburb_share_by_own_income = own_income_suburb_share(ms, hpg),
        city_share_by_own_z = zc, city_share_gradient_z = zc[end] - zc[1],
        teacher_z_gap_sd = tsel.z_gap_sd, teacher_z_gap_market_sd = tsel.z_gap_market_sd,
        gap_logh_teacher = st.logh[2] / nTl[2] - st.logh[1] / nTl[1],
        score_gap_noQ = ab.score_gap_noQ,
        dz_young_sd = cmp.dz_young_sd, dz_adult_sd = cmp.dz_adult_sd,
        gap_meanlog_mkt = cmp.gap_meanlog_mkt, gap_meanlog_mkt_Qadj = cmp.gap_meanlog_mkt_Qadj,
        mover_premium = cmp.mover_premium,
        teacher_premium_m = cmp.teacher_premium[1], teacher_premium_f = cmp.teacher_premium[2],
        gender_gap_market = cmp.gender_gap_market,
        teach_share        = sum(nTl) / sum(adults),
        move_rate          = sum(S[l] * (1 - P[l, l]) for l in 1:L) / sum(S),
        p12 = P[1, 2], p21 = P[2, 1],
        teacher_sd = cfr.sd_logh, cfr_nonteach = cfr.effect_nonteach,
        cfr_analytic = cfr.analytic,
        wtp_share_l = wtp.share_l,
        s_ratio = ss.meanT / ss.meanO,
        students_per_teacher = S ./ nTl,
        t = copy(sol.t), Q = copy(gr.Q), QT = copy(gr.QT), M,
        mean_C = mean_consumption(sol), mean_h = mean_child_h(sol),
        stationarity_dev = maximum(abs.(adults .- S) ./ S),
        minC = wtp.minC, n_binding = wtp.n_binding, floor_mass = wtp.floor_mass,
        infeasible_mass = wtp.infeasible_mass,
        pop_share = M ./ sum(M),
        t_market = hps === nothing ? NaN : hps.t_market,
        t_market_l = hps === nothing ? fill(NaN, L) : hps.t_market_l,
        hp_share_m = hps === nothing ? NaN : hps.share[1],
        hp_share_f = hps === nothing ? NaN : hps.share[2],
        s_O_mean = ss.meanO, s_T_mean = ss.meanT,
        dLambda_z = dΛ, dLambda_mean = dot(dΛ, vec(sum(Φ; dims = 2))) / sum(Φ),
        wtp_income_elasticity = wti.elasticity, wtp_share_by_parent_income = wti.share,
        wtp_income_slope_z = wti.slope_z, wtp_income_elasticity_z = wti.elasticity_z,
        dLambda_quality = qv.dΛ_quality, dLambda_quality_mean = dot(qv.dΛ_quality, vec(sum(Φ; dims = 2))) / sum(Φ),
        payroll_pp_gap = gap(st.pay ./ S),
        commute_share_all = cs.all, commute_share_teach = cs.teach, commute_share_nonteach = cs.nonteach,
        n_teach_res = nTr, n_teach_work = nTl,
        gap_teach_res_emp = gap(nTr ./ (hps === nothing ? adults : adults .- hps.n_l)),
        net_commute_gap = p.split ? gap(nTr) - gap(nTl) : NaN,
        teacher_supply_elasticity = tsr.ratio,
        option_value_teach = ov.teach, option_value_nonteach = ov.nonteach, option_value_gap = ov.gap,
        edu..., wage_slope_ba = ab.wage_slope_ba, teacher_slope_ratio = ab.teacher_slope_ratio,
        pooled_p9010 = meas.pooled_p9010 ? pooled_p9010(sol, meas) : NaN,
        ba_share_by_residence = bsr.share, ba_share_gap = bsr.gap,
        teacher_district_elasticity = tde,
        commute_dir_teach = cdr.teach, commute_dir_nonteach = cdr.nonteach,
        teacher_own_elasticity = tsr.own, teacher_cross_elasticity = tsr.cross,
        move_rate_by_own_z = mz.move, c2s_by_own_z = mz.c2s, s2c_by_own_z = mz.s2c,
        teacher_shadow_b = p.split ? p.BT[2] - p.BT[1] : NaN, mean_log_yhat = mean_predicted_log_income(sol),
        cfr_rho = meas.cfr_rho, cfr_sd_eff = cfr.sd_eff, lyT_fixed = p.lyT_fixed,
        tsm...)
end

# -----------------------------------------------------------------------------
# 6b. Education moments (`Params.college`) and the pooled within-occupation 90/10
# -----------------------------------------------------------------------------
"Composite reliability of the NLSY79 AFQT (the project's latent-score convention: R_comp_pop = 0.9158; nlsy_ba_slopes.md)."
const AFQT_RELIABILITY = 0.916

"""
    ba_share_by_score_terciles(A1, A0, lz; R = AFQT_RELIABILITY) -> Vector

P(BA) by tercile of a noisy score x = √R θ + √(1 − R) u of standardized log ability θ, u ~ N(0, 1): the
NLSY79 BA shares by tercile of the OBSERVED AFQT (D3, D3b) are shares by terciles of such a score, whereas
the terciles of own z (`ba_share_by_z`) are the noise-free version. `A1`, `A0` are the masses of BA and
non-BA adults at the ability nodes of log z `lz`; θ is standardized over all adults and the terciles of x
are cut over all adults (bisection on the mixture CDF).
"""
function ba_share_by_score_terciles(A1, A0, lz; R = AFQT_RELIABILITY)
    w  = A1 .+ A0
    m  = dot(w, lz) / sum(w)
    sd = sqrt(dot(w, (lz .- m) .^ 2) / sum(w))
    μ  = sqrt(R) .* (lz .- m) ./ sd
    σ  = sqrt(1 - R)
    F(c) = dot(w, cdf.(Normal(), (c .- μ) ./ σ)) / sum(w)
    function cut(q)
        a, b = -12.0, 12.0
        for _ in 1:80
            c = (a + b) / 2
            F(c) < q ? (a = c) : (b = c)
        end
        return (a + b) / 2
    end
    c = [-Inf, cut(1 / 3), cut(2 / 3), Inf]
    return map(1:3) do q
        f = cdf.(Normal(), (c[q+1] .- μ) ./ σ) .- cdf.(Normal(), (c[q] .- μ) ./ σ)
        dot(A1, f) / dot(w, f)
    end
end

"""
    college_stats(sol, ms, ss, tsel, hpgrp) -> NamedTuple

The education moments of §3c, all NaN without `p.college` (listed in `model_moments`). Shares and
z gaps are integrals over the cells (`tsel`: the masses of teachers, non-teachers and home
producers by ability node and education from `teacher_selection`; `ss`: `s_summary`), the wage
moments are over the pseudo-sample `ms` (`micro_sample`, market workers: all but the home
producers, teachers included).
- `ba_premium` holds gender fixed (the within-gender differences, weighted by the genders' shares
  among non-teacher market workers): the data's within-occupation premium is a regression with
  occupation × gender cells, and a raw pooled difference would load the BA holders' different
  gender mix onto it. `ba_premium_raw` is the unadjusted difference. The model has no BA-specific
  occupational sorting, so the within-occupation premium is the data counterpart.
- `ba_share_by_z`: P(BA) by tercile of own z, cut as `own_z_composition` (a node that straddles a cut is split);
  `ba_share_by_score` the same by tercile of the NOISY score of the NLSY79 AFQT (`ba_share_by_score_terciles`),
  which is what the observed-score data terciles measure, and `ba_share_score_gap` its top minus bottom.
- `teacher_z_gap_ba_sd`, `teacher_z_gap_within_ba_sd`: all teachers' and BA teachers' mean log z minus
  that of BA non-teacher market workers, in SDs of log z (NLSY79 −0.334 and −0.238).
- `teacher_ba_wage_gap` (`_m`, `_f`): mean log wage of BA minus non-BA teachers, pooled (the data's
  plain pooled difference, convention P: 0.482) and by gender (0.202, 0.566): the data counterpart of the
  non-BA teaching track's pay penalty ρT0, with the schooling gap that the track's time share adds;
  `s_T_mean_nonba`: its mean time share (data 12.99 years/25 = 0.520); `teacher_p9010_ba`
  (`teacher_p9010_ba_g`): the teachers' 90/10 of BA teachers alone (pooled with `teacher_weights`; by
  gender), to be read beside `teacher_p9010`, whose bottom decile is mostly non-BA teachers when
  the track is large or underpaid.
- `ma_premium_implied`: ϑs log(17/16), the salary schedule's master's premium in the units of the
  `s_T_mean_ba` target (the ACS years of schooling are capped at 17, so a graduate degree is one
  year beyond the BA: 43% of BA teachers → 16.43 years); `ma_premium_implied_175` the same with a
  master's taken as 1.5 years beyond the BA (ϑs log(17.5/16), the first release's convention).
"""
function college_stats(sol, ms, ss, tsel, hpgrp::Int; teacher_weights = Measurement().teacher_weights)
    (; Φ, gr, p) = sol
    if !p.college
        return (; ba_share = NaN, ba_share_nonteach = NaN, ba_share_teach = NaN, ba_z_gap_sd = NaN,
                  ba_share_by_z = fill(NaN, 3), ba_share_z_gap = NaN,
                  ba_share_by_score = fill(NaN, 3), ba_share_score_gap = NaN, ba_premium = NaN, ba_premium_raw = NaN,
                  teacher_premium_ba_m = NaN, teacher_premium_ba_f = NaN,
                  teacher_premium_within_ba_m = NaN, teacher_premium_within_ba_f = NaN,
                  teacher_z_gap_ba_sd = NaN, teacher_z_gap_within_ba_sd = NaN,
                  s_T_mean_ba = NaN, s_O_mean_ba = NaN, s_O_mean_nonba = NaN, ma_premium_implied = NaN,
                  ma_premium_implied_175 = NaN, slope_ratio_closed_form = NaN,
                  teacher_ba_wage_gap = NaN, teacher_ba_wage_gap_m = NaN, teacher_ba_wage_gap_f = NaN,
                  s_T_mean_nonba = NaN, teacher_p9010_ba = NaN, teacher_p9010_ba_g = fill(NaN, 2))
    end
    s_z = p.σξ / sqrt(1 - p.ρz^2)
    lz  = log.(gr.z)
    (; T, O, H) = tsel                                     # [zi, d = 1 / d = 0]
    mean_lz(x) = dot(x, lz) / sum(x)
    A1, A0 = T[:, 1] .+ O[:, 1], T[:, 2] .+ O[:, 2]        # BA and non-BA adults by ability node
    M1 = O[:, 1] .- H[:, 1]                                # BA non-teacher market workers
    # P(BA) by tercile of own z
    w   = vec(sum(Φ; dims = 2))
    cum = cumsum(w) ./ sum(w)
    lo  = vcat(0.0, cum[1:end-1])
    ba_z = map(1:3) do q
        f = [max(0.0, min(cum[i], q / 3) - max(lo[i], (q - 1) / 3)) / max(cum[i] - lo[i], 1e-300) for i in eachindex(w)]
        sum(f .* A1) / sum(f .* (A1 .+ A0))
    end
    # wage moments
    mkt  = market_workers(ms, hpgrp)
    nonT = mkt .& (ms.grp .> 0)
    ba, nb = nonT .& (ms.d .== 1), nonT .& (ms.d .== 0)
    mean_at(sel) = wmean(ms.logw[sel], ms.w[sel])
    atg(sel, g) = sel .& (ms.g .== g)
    wg    = [sum(ms.w[atg(nonT, g)]) for g in 1:2]
    prem  = [mean_at(atg(ba, g)) - mean_at(atg(nb, g)) for g in 1:2]
    tprem = [mean_at(atg(ms.grp .== 0, g)) - mean_at(atg(ba, g)) for g in 1:2]
    wprem = [mean_at(atg((ms.grp .== 0) .& (ms.d .== 1), g)) - mean_at(atg(ba, g)) for g in 1:2]
    T1, T0 = (ms.grp .== 0) .& (ms.d .== 1), (ms.grp .== 0) .& (ms.d .== 0)
    tgap  = [mean_at(atg(T1, g)) - mean_at(atg(T0, g)) for g in 1:2]       # NaN if a group is empty
    tw_ba = teacher_wage_stats(ms, teacher_weights; d = 1)
    ba_sc = ba_share_by_score_terciles(A1, A0, lz)
    return (; ba_share = sum(A1) / (sum(A1) + sum(A0)), ba_share_nonteach = sum(O[:, 1]) / sum(O),
              ba_share_teach = sum(T[:, 1]) / sum(T), ba_z_gap_sd = (mean_lz(A1) - mean_lz(A0)) / s_z,
              ba_share_by_z = ba_z, ba_share_z_gap = ba_z[3] - ba_z[1],
              ba_share_by_score = ba_sc, ba_share_score_gap = ba_sc[3] - ba_sc[1],
              ba_premium = dot(wg, prem) / sum(wg), ba_premium_raw = mean_at(ba) - mean_at(nb),
              teacher_premium_ba_m = tprem[1], teacher_premium_ba_f = tprem[2],
              teacher_premium_within_ba_m = wprem[1], teacher_premium_within_ba_f = wprem[2],
              teacher_z_gap_ba_sd = (mean_lz(T[:, 1] .+ T[:, 2]) - mean_lz(M1)) / s_z,
              teacher_z_gap_within_ba_sd = (mean_lz(T[:, 1]) - mean_lz(M1)) / s_z,
              s_T_mean_ba = ss.meanT_ba, s_O_mean_ba = ss.meanO_ba, s_O_mean_nonba = ss.meanO_nonba,
              ma_premium_implied = p.ϑs * log(17 / 16), ma_premium_implied_175 = p.ϑs * log(17.5 / 16),
              slope_ratio_closed_form = p.γ * (1 - p.η) / (1 - p.γ * p.η),
              teacher_ba_wage_gap = mean_at(T1) - mean_at(T0), teacher_ba_wage_gap_m = tgap[1], teacher_ba_wage_gap_f = tgap[2],
              s_T_mean_nonba = ss.meanT_nonba, teacher_p9010_ba = tw_ba.p9010, teacher_p9010_ba_g = tw_ba.p9010_g)
end

"""
    ba_share_by_residence(sol) -> (; share, gap)

The BA share of adults by residence r, and its gap, suburb minus city (levels): the model counterpart of
the share of adults 25+ with a BA or more in a school district, suburb minus city (the sealed hold-out B1 of
`experiments-2026-10d/PREREG.md`; the data are not read here). The adults are the old generation, every
worker including the home producers, teachers and non-teachers; the stationary mass of a cell
(gender, birth location, z) is split by the education weights P(d | g, l, z) and each branch by its own
residence probabilities integrated over the occupation choice, `Pi_Tr`, `Pi_Or` under the split (the
marginals of the joint probabilities over the workplace) and `Pi_T`, `Pi_O` in the bundled model, where
a worker lives where she works. NaN without `p.college`.
"""
function ba_share_by_residence(sol)
    (; hh, Φ, gr, p) = sol
    L = gr.L
    p.college || return (; share = fill(NaN, L), gap = NaN)
    PT, PO = p.split ? (sol.Pi_Tr, sol.Pi_Or) : (sol.Pi_T, sol.Pi_O)
    ba = zeros(L); adults = zeros(L)
    for l in 1:L, zi in 1:gr.Nz, g in 1:2, br in cell_branches(hh, PT, PO, p, g, l, zi)
        mass = Φ[zi, l] / 2 * br.w
        iszero(mass) && continue
        cm = choice_maps(br.h, g, l, zi, gr)
        for r in 1:L
            m = mass * (integrate_teach(br.πT[g, l, zi, :, r], cm, gr) +
                        integrate_nonteach(br.πO[g, l, zi, :, r], cm, gr))
            adults[r] += m
            br.d == 1 && (ba[r] += m)
        end
    end
    share = ba ./ adults
    return (; share, gap = share[2] - share[1])
end

"""
The moments of the education block: `model_moments` returns them as NaN without `p.college`, so a target
among them (active, or validation in the report) needs `--external college=1` (`problem_setup` checks).
"""
const COLLEGE_MOMENTS = (:ba_share, :ba_share_nonteach, :ba_share_teach, :ba_z_gap_sd, :ba_share_by_z, :ba_share_z_gap,
                         :ba_share_by_residence, :ba_share_gap,
                         :ba_share_by_score, :ba_share_score_gap, :ba_premium, :ba_premium_raw,
                         :teacher_premium_ba_m, :teacher_premium_ba_f, :teacher_premium_within_ba_m,
                         :teacher_premium_within_ba_f, :teacher_z_gap_ba_sd, :teacher_z_gap_within_ba_sd,
                         :s_T_mean_ba, :s_O_mean_ba, :s_O_mean_nonba, :s_T_mean_nonba, :ma_premium_implied,
                         :ma_premium_implied_175, :slope_ratio_closed_form, :teacher_ba_wage_gap,
                         :teacher_ba_wage_gap_m, :teacher_ba_wage_gap_f, :teacher_p9010_ba, :teacher_p9010_ba_g,
                         :wage_slope_ba, :teacher_slope_ratio)

"""
    pooled_p9010(sol, meas; n = meas.nquad) -> Float64

The pooled within-occupation 90/10 of non-teachers' hourly wages, the σϵ target of the first stage:
the average of the 38 market cells' (occupation × gender) 90/10s of pre-tax wages with the data's
wage-sample counts as weights (`WAGE_SAMPLE_COUNTS[meas.wage_sample]`), as `occupation_consistency`
forms it (home production, `meas.hp`, is not a market cell). It needs a second pseudo-sample with
one group per non-teaching occupation, so it is slow; `model_moments` computes it only if
`meas.pooled_p9010`.
"""
function pooled_p9010(sol, meas::Measurement; n = meas.nquad)
    groups = [[i] for i in sol.gr.nonteach]
    ms  = micro_sample(sol; n, groups)
    mkt = [k for (k, grp) in enumerate(groups) if grp[1] != meas.hp]
    cnt = WAGE_SAMPLE_COUNTS[meas.wage_sample]
    length(mkt) == length(cnt.men) ||
        throw(ArgumentError("pooled_p9010: $(length(mkt)) market occupations (set `meas.hp`) but $(length(cnt.men)) wage-sample counts"))
    R = [p9010(ms, (ms.grp .== k) .& (ms.g .== g)) for k in mkt, g in 1:2]
    W = hcat(Float64.(cnt.men), Float64.(cnt.women))
    return sum(W .* R) / sum(W)
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
