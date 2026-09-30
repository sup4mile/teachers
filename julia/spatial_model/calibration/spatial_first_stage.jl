# =============================================================================
# First-stage estimation of the ability block (T2b; julia/spatial_model/calibration/spatial_calibration.md
# §1, Table 2, "Ability block").
#
# Outside the equilibrium, with Q_l and t_l fixed and selection into teaching ignored
# for both genders:
#   (i)   ρz and s_z in closed form from the NLSY moments;
#   (ii)  σϵ by a one-dimensional search on the pooled within-occupation 90/10 of
#         non-teachers' hourly wages, with the shares block (A_i, τω_{i,f}) inverted
#         at each value over the 20 non-teaching occupations;
#   (iii) the Nz convergence of (ii) against the continuous-z limit;
#   (iv)  the score loading c = b s_z/χ, the χ gap and the sibling check;
#   (vi)  the Table 5 variants of the closed-form rows, each refitting σϵ.
# Step (v), writing the estimates into `EXTERNAL` in spatial_calibrate.jl, is done by
# hand from this report.
#
# The model. A non-teacher in occupation i, gender g, earns, up to a cell constant,
#     log w = b [α log z + log(Θ_{i,g} ϵ_i^α)],   b = 1/(1 − η),
# with z independent of the chosen occupation and of ϵ (teaching ignored). Roy shares
# depend on log Θ only through log Θ/(ασϵ), so the inversion at any σϵ is ασϵ times the
# inversion at ασϵ = 1: log Θ_{i,g} = ασϵ c_{i,g}. Given that i is chosen,
#     x = log X_i/(ασϵ) has density φ(x − c_i) ∏_{j≠i} Φ(x − c_j) / P_i
# at every σϵ, and a cell's 90/10 is exp(b [q90 − q10]) of α log z + ασϵ x, with log z
# on the model's Rouwenhorst grid (or N(0, s_z²) in the continuous limit). The pooled
# 90/10 averages the 38 market cells with ACS wage-sample weights, as in the data.
#
#   julia --project=julia/spatial_model/calibration/estimation julia/spatial_model/calibration/spatial_first_stage.jl [--write] [--t7]
# `--write` saves the report to data/spatial/estimates/first_stage.{md,toml}; `--t7` adds
# the refit on the 2016–19 block (T7) and, with `--write`, saves first_stage_t7.md.
# =============================================================================

isdefined(@__MODULE__, :occupation_block) || include(joinpath(@__DIR__, "spatial_calibrate.jl"))

using LinearAlgebra
using Printf
using TOML

# -----------------------------------------------------------------------------
# 1. Data
# -----------------------------------------------------------------------------
"""
ACS 2009–13 wage moments of the 19 market occupations, in `ACS_OCCUPATIONS` order,
from data/LaborMarketData/wages_occ_shares_v2.xlsx (wage sample: ages 25–34, at least
48 weeks and 30 hours a week, at least \$1,000 in 2010 dollars):
- `p9010_men`, `p9010_women`: 90/10 of hourly wages (sheet `90_10_hr_wages_weighted`);
- `n_men`, `n_women`: wage-sample counts (sheet `occ_gender_weights`), the weights of
  the pooled 90/10; home production has none, and K–12 teachers are excluded;
- `pooled`: the weighted average of the 38 cell ratios (3.737), the σϵ target;
- `pooled_resid`: the same for age-residualized wages (`90_10_resid_weighted`, 3.635);
- `mean_log_men`, `mean_log_women`, `k12_mean_log`: PWGTP-weighted mean log hourly wages
  (1999 dollars) of the market cells and K–12 teachers (men, women), from the ACS 2009–13
  5-year PUMS with the workbook's rules (`acs_occupations.py`; not in the workbook), for
  the untargeted mean-wage check (T4).
"""
const ACS_WAGES = (;
    p9010_men   = [4.148209, 3.379131, 2.924719, 3.935401, 6.59469, 5.617173, 4.230855, 3.310792,
                   3.903154, 4.7105, 3.267888, 3.507951, 3.5, 4.510296, 3.655806, 3.191318,
                   3.522213, 3.394892, 3.466604],
    p9010_women = [3.590654, 2.974307, 2.986459, 3.238095, 5.480093, 4.527641, 3.970526, 3.581652,
                   3.436822, 4.623172, 2.802961, 3.501195, 3.859886, 4.19272, 4.000065, 3.30715,
                   3.231369, 3.241037, 3.555692],
    n_men   = [51028, 23030, 33560, 20877, 9642, 10527, 4292, 2850, 27754, 54653, 41631, 23058,
               35132, 18053, 76036, 11654, 17580, 27169, 27228],
    n_women = [44932, 28486, 9787, 24365, 8845, 60037, 4239, 8960, 22118, 42370, 92826, 5479,
               40793, 2646, 2016, 3371, 4749, 6584, 2122],
    pooled       = 3.737418491024515,
    pooled_resid = 3.635308119334235,
    mean_log_men   = [2.82123, 2.93778, 2.97909, 2.60984, 3.13275, 2.62955, 2.43289, 2.42986, 2.90428, 2.55834,
                      2.40937, 2.63397, 2.09276, 2.16848, 2.46123, 2.43587, 2.37097, 2.28475, 2.35173],
    mean_log_women = [2.68990, 2.79309, 2.91713, 2.58488, 3.05026, 2.54060, 2.42713, 2.35085, 2.65224, 2.32685,
                      2.34376, 2.52550, 1.97418, 2.01701, 2.36009, 2.19846, 2.07795, 2.11139, 2.15704],
    k12_mean_log   = (2.59897, 2.52467))

const N_MKT = N_OCC - 2    # market occupations 1:N_MKT; HP_OCC = N_MKT + 1, T_OCC = N_MKT + 2
@assert HP_OCC == N_MKT + 1 && T_OCC == N_MKT + 2 "first stage assumes market occupations, then HP, then teaching"

"Wage-sample weights of the 38 market cells, (occupation, gender); `w = ACS_2016_19` for T7."
wage_weights(w = ACS_WAGES) = hcat(Float64.(w.n_men), Float64.(w.n_women))

"Data 90/10s of the 38 market cells, (occupation, gender); `w = ACS_2016_19` for T7."
data_p9010(w = ACS_WAGES) = hcat(w.p9010_men, w.p9010_women)

"""
NLSY moments from data/spatial/estimates/nlsy_ability.json (T2a, `nlsy_ability.py`,
200 cluster-bootstrap draws):
- `rho_pc`: latent mother–child correlation (two-factor ULS), the ρz target;
  `rho_pc_xs` (cross-sectional mothers) and `rho_pc_unw` (unweighted) are Table 5 rows;
- `slope`: log-wage slope per latent SD, NLSY79, ACS cutoffs, the b α s_z target;
  `slope97` is the NLSY97 vintage check;
- `rho_sib`: latent sibling correlation (check only; the model gives ρz²).
"""
const NLSY = (;
    rho_pc  = 0.5766299079601815, rho_pc_se  = 0.014023593002952637,
    rho_pc_xs  = 0.5292689017161321,
    rho_pc_unw = 0.6299783704846894,
    slope   = 0.21266045167782227, slope_se   = 0.007708547769715874,
    slope97 = 0.19872689843795754, slope97_se = 0.009468583649377938,
    rho_sib = 0.6062739952268316,  rho_sib_se = 0.016531908660971692)

"Earnings–achievement bridge χ = log 1.12 (Table 1), as in `Measurement`."
const CHI = Measurement().χ

# -----------------------------------------------------------------------------
# 2. Closed-form rows (i) and (iv)
# -----------------------------------------------------------------------------
"""
    closed_form(; slope, ρ, η, α, χ) -> NamedTuple

(i) ρz = ρ and s_z = slope (1 − η)/α, so that the model's slope b α s_z equals the NLSY
slope; σξ = s_z √(1 − ρz²). (iv) c = b s_z/χ, the χ gap b s_z − χ, and the model's
latent sibling correlation ρz².
"""
function closed_form(; slope = NLSY.slope, ρ = NLSY.rho_pc, η = EXTERNAL.η, α = EXTERNAL.α, χ = CHI)
    b   = 1 / (1 - η)
    s_z = slope * (1 - η) / α
    return (; ρz = ρ, s_z, σξ = innovation_sd(s_z, ρ), b, c = b * s_z / χ,
              chi_gap = b * α * s_z - χ, rho_sib = ρ^2)
end

# -----------------------------------------------------------------------------
# 3. The ability distribution
# -----------------------------------------------------------------------------
"Nodes and weights of the n-point Gauss–Hermite rule for N(0, 1) (Golub–Welsch)."
function gauss_hermite_normal(n)
    F = eigen(SymTridiagonal(zeros(n), sqrt.(1.0:n-1)))
    return F.values, F.vectors[1, :] .^ 2
end

"""
    zdist(Nz, s_z; ρz) -> (; ℓ, π)

Nodes ℓ = log z and stationary probabilities π: the model's Rouwenhorst grid
(`rouwenhorst`, `stationary`) for an integer Nz, or a 64-node Gauss–Hermite rule for
N(0, s_z²) at `Nz = :continuous`. Rouwenhorst end nodes are ±√(Nz − 1) s_z and the
stationary law is binomial, so neither depends on ρz given s_z. In the GE, Φ summed
over locations is this stationary law, since each parent has one child.
"""
function zdist(Nz, s_z; ρz = NLSY.rho_pc)
    iszero(s_z) && return (; ℓ = [0.0], π = [1.0])
    if Nz === :continuous
        u, w = gauss_hermite_normal(64)
        return (; ℓ = s_z .* u, π = w)
    end
    z, Π = rouwenhorst(Nz, ρz, innovation_sd(s_z, ρz))
    return (; ℓ = log.(z), π = stationary(Π))
end

# -----------------------------------------------------------------------------
# 4. The Roy-selected occupation–gender cells
# -----------------------------------------------------------------------------
"""
    RoyCells(; h = 0.002, pad = 9.0, data = ACS_OCCUPATIONS)

The shares block at unit dispersion ασϵ = 1, on the occupational sample `data`:
- `c[i, g]`: log Θ_{i,g}/(ασϵ) over the 20 non-teaching occupations, c_HP = 0, by
  `invert_roy_shares` at s = 1 (the inversion at any s is s·c);
- `x`, `G[:, i, g]`: a uniform grid and the CDF of x = log X_i/(ασϵ) given that i is
  chosen, by cumulative trapezoid of φ(x − c_i) ∏_{j≠i} Φ(x − c_j);
- `P[i, g]`: the implied shares, a check on the inversion and the quadrature.
"""
struct RoyCells
    c::Matrix{Float64}
    x::Vector{Float64}
    G::Array{Float64,3}
    P::Matrix{Float64}
end

function RoyCells(; h = 0.002, pad = 9.0, data = ACS_OCCUPATIONS)
    n = N_OCC - 1
    c = hcat(invert_roy_shares(occ_shares(1; data), 1.0; ref = HP_OCC),
             invert_roy_shares(occ_shares(2; data), 1.0; ref = HP_OCC))
    x = collect(range(minimum(c) - pad, maximum(c) + pad; step = h))
    G = zeros(length(x), n, 2)
    P = zeros(n, 2)
    N01 = Normal()
    for g in 1:2
        lΦ = [logcdf(N01, xk - c[j, g]) for xk in x, j in 1:n]
        S  = vec(sum(lΦ; dims = 2))
        for i in 1:n
            f = exp.(logpdf.(N01, x .- c[i, g]) .+ S .- lΦ[:, i])
            cum = zeros(length(x))
            for k in 2:length(x)
                cum[k] = cum[k-1] + h * (f[k-1] + f[k]) / 2
            end
            P[i, g] = cum[end]
            G[:, i, g] = cum ./ cum[end]
        end
    end
    return RoyCells(c, x, G, P)
end

"CDF of x in cell (i, g), by linear interpolation; 0 below and 1 above the grid."
function cellcdf(rc::RoyCells, i, g, xq)
    x = rc.x
    xq <= x[1] && return 0.0
    xq >= x[end] && return 1.0
    h = x[2] - x[1]
    k = min(floor(Int, (xq - x[1]) / h) + 1, length(x) - 1)
    t = (xq - x[k]) / h
    return (1 - t) * rc.G[k, i, g] + t * rc.G[k+1, i, g]
end

"Root of f on [a, b], which must bracket a sign change, by bisection."
function bisect(f, a, b; tol = 1e-12, maxit = 200)
    fa, fb = f(a), f(b)
    fa * fb <= 0 || error("bisect: no sign change on [$a, $b] (f = $fa, $fb)")
    for _ in 1:maxit
        m  = (a + b) / 2
        fm = f(m)
        (fm == 0 || b - a < tol) && return m
        sign(fm) == sign(fa) ? ((a, fa) = (m, fm)) : ((b, fb) = (m, fm))
    end
    return (a + b) / 2
end

"""
    cell_p9010(rc, i, g, s, zd; α, b)

90/10 of wages in market cell (i, g): exp(b [q90 − q10]) of W = α log z + s x, with
x ~ G_{i,g}, log z ~ `zd` independent of x, and s = ασϵ.
"""
function cell_p9010(rc::RoyCells, i, g, s, zd; α = EXTERNAL.α, b)
    F(w) = sum(zd.π[k] * cellcdf(rc, i, g, (w - α * zd.ℓ[k]) / s) for k in eachindex(zd.π))
    lo = s * rc.x[1] + α * minimum(zd.ℓ) - 1
    hi = s * rc.x[end] + α * maximum(zd.ℓ) + 1
    q(p) = bisect(w -> F(w) - p, lo, hi; tol = 1e-10)
    return exp(b * (q(0.9) - q(0.1)))
end

"p-quantile of a discrete law with ascending nodes: the smallest node whose CDF reaches p."
dquantile(zd, p) = zd.ℓ[findfirst(≥(p - 1e-12), cumsum(zd.π))]

"All 38 market-cell 90/10s, (occupation, gender)."
cell_p9010s(rc::RoyCells, s, zd; α = EXTERNAL.α, b) =
    [cell_p9010(rc, i, g, s, zd; α, b) for i in 1:N_MKT, g in 1:2]

"Wage-sample-weighted average of cell ratios: the pooled 90/10 of the model or the data."
pooled(R; w = wage_weights()) = sum(w .* R) / sum(w)

# -----------------------------------------------------------------------------
# 5. The σϵ search (ii)
# -----------------------------------------------------------------------------
"""
    fit_sigma_eps(rc; s_z, ρz, η, α, Nz, target) -> NamedTuple

σϵ matching the pooled 90/10 `target`, given s_z, on the ability grid `Nz` (an integer
or `:continuous`). Throws if the ability component alone already exceeds the target
(§1, caveat 1). Returns σϵ, the cell 90/10s `R`, and the 90/10 of z alone: `zonly` for
continuous z (the corner check) and `zonly_grid` on the grid, where q10 and q90 fall
on lattice nodes.
"""
function fit_sigma_eps(rc::RoyCells; s_z, ρz = NLSY.rho_pc, η = EXTERNAL.η, α = EXTERNAL.α,
                       Nz = SOLVER.Nz, target = ACS_WAGES.pooled, w = wage_weights())
    b  = 1 / (1 - η)
    zd = zdist(Nz, s_z; ρz)
    gap(s) = pooled(cell_p9010s(rc, s, zd; α, b); w) - target
    zonly = exp(2 * quantile(Normal(), 0.9) * b * α * s_z)   # continuous z; the corner check
    zonly < target || error(@sprintf("σϵ at a corner: z alone gives a 90/10 of %.3f ≥ %.3f", zonly, target))
    zonly_grid = exp(b * α * (dquantile(zd, 0.9) - dquantile(zd, 0.1)))
    s = bisect(gap, 1e-3, 5.0; tol = 1e-9)
    R = cell_p9010s(rc, s, zd; α, b)
    return (; σϵ = s / α, s, R, pooled = pooled(R; w), zonly, zonly_grid, b, Nz)
end

# -----------------------------------------------------------------------------
# 6. Estimates, Nz convergence and Table 5 variants
# -----------------------------------------------------------------------------
"The Table 5 variants of the first stage: each changes one input and refits σϵ."
function table5_variants()
    base = (; slope = NLSY.slope, ρ = NLSY.rho_pc, η = EXTERNAL.η, target = ACS_WAGES.pooled)
    return [
        ("Baseline",                            base),
        ("ρz: cross-sectional mothers (0.529)", merge(base, (; ρ = NLSY.rho_pc_xs))),
        ("ρz: unweighted (0.630)",              merge(base, (; ρ = NLSY.rho_pc_unw))),
        ("b s_z = χ",                           merge(base, (; slope = CHI))),
        ("NLSY97 slope (0.199)",                merge(base, (; slope = NLSY.slope97))),
        ("Age-residualized 90/10 (3.635)",      merge(base, (; target = ACS_WAGES.pooled_resid))),
        ("η = 0.073",                           merge(base, (; η = 0.073))),
        ("η = 0.103",                           merge(base, (; η = 0.103))),
    ]
end

"""
    first_stage(; Nz = SOLVER.Nz, Nz_grid = (3, 5, 7, 9, 15, 25)) -> NamedTuple

Run (i)–(iv) and (vi): the baseline estimates on the full model's grid `Nz`, the Nz
convergence table, and the Table 5 variants.
"""
function first_stage(; Nz = SOLVER.Nz, Nz_grid = (3, 5, 7, 9, 15, 25))
    rc  = RoyCells()
    cf  = closed_form()
    fit = fit_sigma_eps(rc; s_z = cf.s_z, cf.ρz, Nz)

    # (iii) Nz convergence: σϵ̂ on each grid, and each grid's pooled 90/10 at the
    # continuous-limit σϵ̂
    lim   = fit_sigma_eps(rc; s_z = cf.s_z, cf.ρz, Nz = :continuous)
    grids = [Nz_grid..., :continuous]
    conv  = map(grids) do N
        f = N === :continuous ? lim : fit_sigma_eps(rc; s_z = cf.s_z, cf.ρz, Nz = N)
        p = pooled(cell_p9010s(rc, lim.s, zdist(N, cf.s_z; cf.ρz); α = EXTERNAL.α, cf.b))
        (; Nz = N, σϵ = f.σϵ, p9010_at_lim = p)
    end

    variants = map(table5_variants()) do (name, v)
        c = closed_form(; v.slope, v.ρ, v.η)
        f = fit_sigma_eps(rc; s_z = c.s_z, c.ρz, v.η, Nz, target = v.target)
        (; name, v..., c..., σϵ = f.σϵ, zonly = f.zonly)
    end

    occ = occupation_block(; σϵ = fit.σϵ, α = EXTERNAL.α)
    return (; rc, cf, fit, lim, conv, variants, occ, Nz)
end

# -----------------------------------------------------------------------------
# 7. Report
# -----------------------------------------------------------------------------
"Short occupation labels for tables, in `ACS_OCCUPATIONS` order."
const SHORT_NAMES = ["Executives", "Management related", "Architects, engineers, CS", "Scientists, arts",
                     "Doctors and lawyers", "Nurses, therapists", "Postsecondary teachers",
                     "Other teachers, librarians", "Technicians", "Sales", "Administrative support",
                     "Fire, police, guards", "Food, cleaning, personal", "Farm, extraction",
                     "Mechanics, construction", "Precision manufacturing", "Manufacturing operators",
                     "Fabricators, handlers", "Vehicle operators", "Home production", "K–12 teachers"]
short_name(i) = SHORT_NAMES[i]

function report(io::IO, r)
    (; rc, cf, fit, lim, conv, variants, occ, Nz) = r
    mkt  = 1:N_MKT
    D    = data_p9010()
    w    = wage_weights()
    s    = fit.s
    b    = cf.b
    zd0  = zdist(Nz, 0.0)
    Rε   = cell_p9010s(rc, s, zd0; b)                      # ϵ-only cell ratios at σϵ̂

    println(io, "# First-stage estimates of the ability block (T2b)\n")
    println(io, "Generated ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()), " by `julia/spatial_model/calibration/spatial_first_stage.jl`. ",
            "First pass outside the equilibrium: Q_l and t_l fixed, selection into teaching ignored for both genders ",
            "(spatial_calibration.md §1, Table 2, *Ability block*). η = ", EXTERNAL.η, ", α = ", EXTERNAL.α,
            @sprintf(", b = 1/(1 − η) = %.5f, χ = log 1.12 = %.5f.\n", b, CHI))

    println(io, "## Estimates\n")
    println(io, "| parameter | estimate | target | source |")
    println(io, "|---|---|---|---|")
    @printf(io, "| ρz | %.4f | latent mother–child correlation 0.577 (SE %.3f) | NLSY79/CNLSY; closed form |\n", cf.ρz, NLSY.rho_pc_se)
    @printf(io, "| s_z | %.4f (SE %.4f) | latent log-wage slope 0.2127 (SE %.4f) = b α s_z | NLSY79; closed form |\n",
            cf.s_z, NLSY.slope_se * (1 - EXTERNAL.η), NLSY.slope_se)
    @printf(io, "| σξ = s_z √(1 − ρz²) | %.4f | — | implied |\n", cf.σξ)
    @printf(io, "| σϵ (Nz = %d) | %.4f | pooled within-occupation 90/10 %.4f | ACS 2009–13; model %.4f |\n",
            Nz, fit.σϵ, ACS_WAGES.pooled, fit.pooled)
    @printf(io, "| σϵ (continuous z) | %.4f | same | — |\n", lim.σϵ)
    @printf(io, "| c = b s_z/χ | %.4f | — | score loading of log Q_l |\n\n", cf.c)

    println(io, "## Checks\n")
    @printf(io, "- **χ gap:** b α s_z − χ = %.4f − %.4f = %.4f (SE %.4f): the NLSY slope is %.1f SEs above χ.\n",
            b * EXTERNAL.α * cf.s_z, CHI, cf.chi_gap, NLSY.slope_se, cf.chi_gap / NLSY.slope_se)
    @printf(io, "- **Siblings:** model latent sibling correlation ρz² = %.4f, against %.4f (SE %.4f) in the data (check only; siblings share inputs outside the model).\n",
            cf.rho_sib, NLSY.rho_sib, NLSY.rho_sib_se)
    @printf(io, "- **Corner (caveat 1):** z alone gives a 90/10 of %.3f, well below %.3f, so σϵ is interior.\n",
            fit.zonly, ACS_WAGES.pooled)
    @printf(io, "- **Components:** at σϵ̂ the pooled 90/10 is %.3f with ϵ alone (z shut off); z alone gives %.3f for continuous z and %.3f on the Nz = %d lattice, where q10 and q90 fall on nodes. ",
            pooled(Rε), fit.zonly, fit.zonly_grid, Nz)
    @printf(io, "The no-selection lognormal benchmark (`sigma_eps_placeholder`) gives σϵ = %.4f, so Roy selection within cells raises σϵ̂ by %.0f%%.\n",
            sigma_eps_placeholder(; η = EXTERNAL.η, s_z = cf.s_z), 100 * (fit.σϵ / sigma_eps_placeholder(; η = EXTERNAL.η, s_z = cf.s_z) - 1))
    @printf(io, "- **Shares block:** the inversion reproduces the share targets to %.1e (max abs. error, both genders, after quadrature).\n",
            maximum(abs, rc.P .- hcat(occ_shares(1), occ_shares(2))))
    lA = log.(occ.A[mkt]); lr = log.(occ.r[mkt])
    @printf(io, "- **Occupational block at σϵ̂:** log A_i/A_HP ∈ [%.3f, %.3f]; log r_i ∈ [%.3f, %.3f] (market occupations; A_HP = r_HP = 1).\n",
            extrema(lA)..., extrema(lr)...)
    Rm = fit.R
    @printf(io, "- **Cell 90/10s (untargeted):** model by gender %.3f (men), %.3f (women) vs. data %.3f, %.3f; ",
            pooled(Rm[:, 1]; w = w[:, 1]), pooled(Rm[:, 2]; w = w[:, 2]), pooled(D[:, 1]; w = w[:, 1]), pooled(D[:, 2]; w = w[:, 2]))
    @printf(io, "weighted correlation of model and data cell ratios %.3f; model range [%.2f, %.2f] vs. data [%.2f, %.2f].\n\n",
            wcorr(vec(Rm), vec(D), vec(w)), extrema(Rm)..., extrema(D)...)

    println(io, "## Ability grid (Nz)\n")
    println(io, "σϵ̂ on each Rouwenhorst grid, and each grid's pooled 90/10 at the continuous-limit σϵ̂ = ",
            @sprintf("%.4f", lim.σϵ), " (target ", @sprintf("%.4f", ACS_WAGES.pooled), ").\n")
    println(io, "| Nz | σϵ̂ | Δ vs. limit | pooled 90/10 at limit σϵ̂ |")
    println(io, "|---|---|---|---|")
    for c in conv
        @printf(io, "| %s | %.4f | %+.4f | %.4f |\n", string(c.Nz), c.σϵ, c.σϵ - lim.σϵ, c.p9010_at_lim)
    end
    println(io)

    println(io, "## Table 5 variants\n")
    println(io, "Each row changes one input and refits σϵ on the Nz = $Nz grid. Rouwenhorst nodes and weights depend on s_z but not ρz, ",
            "so the ρz rows change only ρz and σξ.\n")
    println(io, "| variant | ρz | s_z | σξ | σϵ̂ | c | b s_z − χ | ρz² |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for v in variants
        @printf(io, "| %s | %.4f | %.4f | %.4f | %.4f | %.3f | %.4f | %.3f |\n",
                v.name, v.ρz, v.s_z, v.σξ, v.σϵ, v.c, v.chi_gap, v.rho_sib)
    end
    println(io)

    println(io, "## Occupational block and cell 90/10s at σϵ̂\n")
    println(io, "A_i/A_HP from men's shares, r_i = Θ_{i,f}/(ωf A_i) from women's (1 − τω_{i,f} = ωf r_i; ωf is internal). ",
            "Shares are over the 20 non-teaching occupations. Cell 90/10s are untargeted checks on a common σϵ.\n")
    println(io, "| occupation | share m | share f | A_i/A_HP | r_i | 90/10 m: model | data | 90/10 f: model | data |")
    println(io, "|---|---|---|---|---|---|---|---|---|")
    sm, sf = occ_shares(1), occ_shares(2)
    for i in 1:N_OCC-1
        if i <= N_MKT
            @printf(io, "| %s | %.4f | %.4f | %.4f | %.4f | %.2f | %.2f | %.2f | %.2f |\n", short_name(i), sm[i], sf[i],
                    occ.A[i], occ.r[i], Rm[i, 1], D[i, 1], Rm[i, 2], D[i, 2])
        else
            @printf(io, "| %s | %.4f | %.4f | 1 | 1 | — | — | — | — |\n", short_name(i), sm[i], sf[i])
        end
    end
    println(io)

    println(io, "## Limitations of the first pass\n")
    println(io, "- Teaching selection is ignored in shares and wage distributions for both genders; the §4 consistency check measures it after the first internal fit.")
    println(io, "- Cross-location variation in Q_l and t_l and the goods moving cost (Ξ) are left out of wages; both add dispersion, so σϵ̂ is an upper bound in that respect.")
    println(io, "- ACS wages include measurement error and transitory shocks, which also load on σϵ (caveat 2).")
    println(io, "- The ACS block is 2009–13, the spatial targets late 2010s (caveat 3).")
    println(io, "- Mean wages by occupation (Table 2 check) are not in the workbook and are not compared here.")
end
report(r) = report(stdout, r)

"Weighted correlation."
function wcorr(x, y, w)
    w = w ./ sum(w)
    mx, my = sum(w .* x), sum(w .* y)
    return sum(w .* (x .- mx) .* (y .- my)) / sqrt(sum(w .* (x .- mx) .^ 2) * sum(w .* (y .- my) .^ 2))
end

"Machine-readable summary of the first stage."
function first_stage_dict(r)
    (; cf, fit, lim, conv, variants, occ, Nz) = r
    return Dict{String,Any}(
        "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
        "note" => "T2b first stage (spatial_first_stage.jl): Q_l, t_l fixed, teaching selection ignored.",
        "eta" => EXTERNAL.η, "alpha" => EXTERNAL.α, "chi" => CHI, "b" => cf.b,
        "Nz" => Nz,
        "rho_z" => cf.ρz, "s_z" => cf.s_z, "sigma_xi" => cf.σξ,
        "sigma_eps" => fit.σϵ, "sigma_eps_continuous" => lim.σϵ,
        "pooled_p9010_model" => fit.pooled, "pooled_p9010_target" => ACS_WAGES.pooled,
        "c" => cf.c, "chi_gap" => cf.chi_gap, "rho_sib_model" => cf.rho_sib,
        "A" => occ.A[1:N_OCC-1], "r" => occ.r[1:N_OCC-1],
        "nz_convergence" => [Dict("Nz" => string(c.Nz), "sigma_eps" => c.σϵ, "p9010_at_limit" => c.p9010_at_lim) for c in conv],
        "variants" => [Dict("name" => v.name, "slope" => v.slope, "rho_z" => v.ρz, "eta" => v.η,
                            "target" => v.target, "s_z" => v.s_z, "sigma_xi" => v.σξ,
                            "sigma_eps" => v.σϵ, "c" => v.c, "chi_gap" => v.chi_gap) for v in variants])
end

# -----------------------------------------------------------------------------
# 8. T7: the 2016–19 block
# -----------------------------------------------------------------------------
"""
    t7_report(io, r)

Refit the first stage on the 2016–19 block (`ACS_2016_19`, T7) with the baseline ability
block (ρz, s_z) and compare it with the 2009–13 estimates in `r = first_stage()`: σϵ̂
against each sample's pooled 90/10, the occupational block (log A_i, log r_i), teaching
shares, teachers' 90/10s and schooling.
"""
function t7_report(io::IO, r)
    (; cf, fit, occ, Nz) = r
    d7  = ACS_2016_19
    rc7 = RoyCells(; data = d7)
    w7  = wage_weights(d7)
    f7  = fit_sigma_eps(rc7; s_z = cf.s_z, cf.ρz, Nz, target = d7.pooled, w = w7)
    o7  = occupation_block(; σϵ = f7.σϵ, α = EXTERNAL.α, data = d7)
    mkt = 1:N_MKT
    lA0, lA7 = log.(occ.A[mkt]), log.(o7.A[mkt])
    lr0, lr7 = log.(occ.r[mkt]), log.(o7.r[mkt])
    println(io, "# First stage on the 2016–19 occupational block (T7)\n")
    println(io, "Generated ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()), " by `spatial_first_stage.jl --t7`. ",
            "Data: pooled ACS 1-year PUMS 2016–19 with the workbook's definitions (`acs_occupations.py`; ",
            "[report](acs_occupations.md)); the ability block (ρz, s_z) is unchanged. Same first pass as T2b: ",
            "Q_l and t_l fixed, teaching selection ignored.\n")
    println(io, "| | 2009–13 | 2016–19 |")
    println(io, "|---|---|---|")
    @printf(io, "| Pooled non-teacher 90/10 (target) | %.4f | %.4f |\n", ACS_WAGES.pooled, d7.pooled)
    @printf(io, "| σϵ̂ (Nz = %d) | %.4f | %.4f |\n", Nz, fit.σϵ, f7.σϵ)
    @printf(io, "| K–12 teaching share, men / women | %.4f / %.4f | %.4f / %.4f |\n",
            teach_share_data(1), teach_share_data(2), teach_share_data(1; data = d7), teach_share_data(2; data = d7))
    @printf(io, "| Home production share, men / women | %.4f / %.4f | %.4f / %.4f |\n",
            occ_shares(1)[HP_OCC], occ_shares(2)[HP_OCC], occ_shares(1; data = d7)[HP_OCC], occ_shares(2; data = d7)[HP_OCC])
    @printf(io, "| Teachers' 90/10, men / women / pooled | %.3f / %.3f / %.3f | %.3f / %.3f / %.3f |\n",
            TEACHER_P9010.men, TEACHER_P9010.women, TEACHER_P9010.pooled,
            d7.teacher_p9010.men, d7.teacher_p9010.women, d7.teacher_p9010.pooled)
    @printf(io, "| Years of schooling, non-teachers / K–12 teachers | %.2f / %.2f | %.2f / %.2f |\n",
            SCHOOLING.nonteach, SCHOOLING.teach, d7.schooling.nonteach, d7.schooling.teach)
    @printf(io, "| log A_i/A_HP range | [%.3f, %.3f] | [%.3f, %.3f] |\n", extrema(lA0)..., extrema(lA7)...)
    @printf(io, "| log r_i range | [%.3f, %.3f] | [%.3f, %.3f] |\n\n", extrema(lr0)..., extrema(lr7)...)
    avg(x) = sum(x) / length(x)
    sd(x)  = sqrt(sum(abs2, x .- avg(x)) / (length(x) - 1))
    eq     = ones(length(mkt))
    @printf(io, "Across the 19 market occupations, log A_i/A_HP moves by %+.3f on average (SD %.3f; correlation %.3f) and log r_i by %+.3f (SD %.3f; correlation %.3f). ",
            avg(lA7 .- lA0), sd(lA7 .- lA0), wcorr(lA0, lA7, eq), avg(lr7 .- lr0), sd(lr7 .- lr0), wcorr(lr0, lr7, eq))
    @printf(io, "Market A_i rise against home production because home-production shares fell (men %.1f pp, women %.1f pp); ",
            100 * (occ_shares(1; data = d7)[HP_OCC] - occ_shares(1)[HP_OCC]), 100 * (occ_shares(2; data = d7)[HP_OCC] - occ_shares(2)[HP_OCC]))
    println(io, "the internal levels (κ̃, ωf) move with them.\n")
    println(io, "| occupation | log A_i 09–13 | log A_i 16–19 | log r_i 09–13 | log r_i 16–19 | 90/10 m: model / data 16–19 | 90/10 f: model / data 16–19 |")
    println(io, "|---|---|---|---|---|---|---|")
    D7 = data_p9010(d7)
    for i in mkt
        @printf(io, "| %s | %.3f | %.3f | %.3f | %.3f | %.2f / %.2f | %.2f / %.2f |\n", short_name(i),
                lA0[i], lA7[i], lr0[i], lr7[i], f7.R[i, 1], D7[i, 1], f7.R[i, 2], D7[i, 2])
    end
    return (; fit7 = f7, occ7 = o7)
end

const FIRST_STAGE_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "spatial", "estimates"))

function main(args = ARGS)
    t = @elapsed r = first_stage()
    report(r)
    @printf("\n(first stage in %.1fs)\n", t)
    if "--write" in args
        open(io -> report(io, r), joinpath(FIRST_STAGE_DIR, "first_stage.md"), "w")
        open(io -> TOML.print(io, first_stage_dict(r); sorted = true), joinpath(FIRST_STAGE_DIR, "first_stage.toml"), "w")
        println("Wrote first_stage.md and first_stage.toml to ", FIRST_STAGE_DIR)
    end
    if "--t7" in args
        println()
        t7_report(stdout, r)
        "--write" in args && (open(io -> t7_report(io, r), joinpath(FIRST_STAGE_DIR, "first_stage_t7.md"), "w");
                              println("Wrote first_stage_t7.md to ", FIRST_STAGE_DIR))
    end
    return r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
