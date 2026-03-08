################################################################################
#  AUTOMATED CALIBRATION OF λf AND κ
#  ------------------------------------------------------------------
#  Finds (λf, κ) to match target shares of women and men in teaching.
#
#  ALGORITHM:
#
#    Exploits known monotonicity:
#      - female teaching share is DECREASING in λf  (for fixed κ)
#      - male   teaching share is INCREASING in κ   (for fixed λf)
#      - cross-effects are weak (GE only, through the tax rate)
#
#    PHASE 1 — Fast Gauss-Seidel root-finding
#      Alternately bisect on κ (to nail the male share) and on λf (to nail
#      the female share) using Brent's method for each 1-D sub-problem.
#      Uses loose fixed-point tolerances for speed.
#      Typically converges in 3-4 outer sweeps × ~12 Brent evals each.
#      Can be skipped if a good initial guess is available.
#
#    PHASE 2 — Broyden quasi-Newton refinement
#      2D root-finding with tight fixed-point tolerances for accuracy.
#      Builds an initial Jacobian via finite differences, then uses
#      Broyden rank-1 updates. Backtracking line search for robustness.
#
#  NOTE: This script does NOT save any results or modify any files.
#  Verify the calibrated parameters by running frechet.jl and save results.
################################################################################

using Distributions, Dierckx, QuadGK, JLD, LaTeXStrings, CSV, DataFrames,
      LinearAlgebra, Optim, Roots, Random

###############################################################################
#  USER CONFIGURATION
###############################################################################

# ── Targets ──────────────────────────────────────────────────────────────────
target_mass_T_female = 0.034   # target share of women in teaching
target_mass_T_male   = 0.012   # target share of men in teaching

# ── Initial guess for (λf, κ) ───────────────────────────────────────────────
λf_init = 0.6436
κ_init  = 0.18325

# ── Search bounds (used for bracketing) ──────────────────────────────────────
λf_lo, λf_hi = 0.30, 1.25
κ_lo,  κ_hi  = 0.11, 0.45

# ── Environment settings ─────────────────────────────────────────────────────
paramname  = "benchmark_w_growth_HP"
year       = 1990
growth_pct = 0.02
include_HP = 1

# ── Phase 1 settings (loose tolerances, fast) ───────────────────────────────
skip_phase1          = true    # skip Phase 1 and go straight to Phase 2
phase1_outer_maxiter = 8
phase1_outer_tol     = 2e-3    # relative tolerance on both shares
phase1_brent_atol    = 1e-4
phase1_brent_maxiter = 25

# Phase 1 fixed-point tolerances (loose for speed)
phase1_tolHH = 1e-5
phase1_tolT  = 1e-5
phase1_tolG  = 1e-5
phase1_ν_tax = 1.0
phase1_ν_H   = 0.8

# ── Phase 2 settings (tight tolerances, accurate) ───────────────────────────
phase2_max_newton_iter = 10
phase2_tol             = 1e-6  # convergence: max relative error on both shares
phase2_fd_step_λf      = 1e-3
phase2_fd_step_κ       = 5e-4
phase2_max_step_λf     = 0.05
phase2_max_step_κ      = 0.02
phase2_armijo_c        = 1e-4
phase2_min_step_factor = 0.1

# Phase 2 fixed-point tolerances (tight for accuracy)
phase2_tolHH = 1e-7
phase2_tolT  = 1e-7
phase2_tolG  = 1e-7
phase2_ν_tax = 0.9
phase2_ν_H   = 0.7

# ── Shared fixed-point settings ─────────────────────────────────────────────
maxiterT = 200
maxiterG = 200

# ── Diagnostic output ────────────────────────────────────────────────────────
verbose_model    = false
verbose_rootfind = true

###############################################################################
#  DERIVED SETTINGS
###############################################################################
A_idx  = include_HP == 1 ? 1 : 2
growth = (1 + growth_pct)^(year - 1970)
M      = 1.0
g      = ["female", "male"]
n_G    = length(g)

###############################################################################
#  DATA LOADING
###############################################################################

base_dir  = pwd()
data_dir  = joinpath(base_dir, "data", "LaborMarketData")
code_dir  = joinpath(base_dir, "julia", "codes_w_exo_wage")
param_dir = joinpath(code_dir, "parameterization", paramname)

import XLSX

xs1  = XLSX.readxlsx(joinpath(data_dir, "wages_occ_shares_v2.xlsx"))
tab1 = xs1["moments_shares"]
xs2  = XLSX.readxlsx(joinpath(data_dir, "wages_occ_shares_v2.xlsx"))
tab2 = xs2["90_10_hr_wages_weighted"]

occ_O = tab1["A30:A49"]
occ_T = tab1["A52"]
occ   = [occ_O; occ_T]

n_O = length(occ)

share_occ_data = Array{Float64,2}(undef, n_O - 1, 2)
w_90_10_data   = Array{Float64,1}(undef, 2)

if year == 1970
    share_occ_data[:, 1] = tab1["K30:K49"]
    share_occ_data[:, 2] = tab1["C30:C49"]
    w_90_10_data[1]      = tab2["C50"]
    w_90_10_data[2]      = tab2["C51"]
elseif year == 1990
    share_occ_data[:, 1] = tab1["M30:M49"]
    share_occ_data[:, 2] = tab1["E30:E49"]
    w_90_10_data[1]      = tab2["E50"]
    w_90_10_data[2]      = tab2["E51"]
elseif year == 2010
    share_occ_data[:, 1] = tab1["O30:O49"]
    share_occ_data[:, 2] = tab1["G30:G49"]
    w_90_10_data[1]      = tab2["G50"]
    w_90_10_data[2]      = tab2["G51"]
else
    error("Unsupported year: $year")
end

###############################################################################
#  LOAD PREVIOUS PARAMETERIZATION
###############################################################################

fyear    = string(year)
fnameJLD = joinpath(param_dir, string("previousParameterization", fyear, ".jld"))
d        = load(fnameJLD)

a_by_occ_initial = d["a_by_occ"]
τ_w_initial      = d["τ_w_opt"]
τ_e_initial      = d["τ_e"]
a_T_thresh_init  = d["a_T_thresh"]
t_init           = d["t"]
H_grid_init      = d["H_grid"]
H_O_init         = d["H_O"]
HH_fp_init       = d["HH_fp"]
α     = d["α"]
β     = d["β"]
η     = d["η"]
σ     = d["σ"]
μ     = d["μ"]
ϕ     = d["ϕ"]
γ     = d["γ"]
theta = d["theta"]
λm    = d["λm"]
a_grid_init  = d["a_grid"]
a_O_10p_init = d["a_O_10p"]
a_O_90p_init = d["a_O_90p"]
a_T_10p_init = d["a_T_10p"]
a_T_90p_init = d["a_T_90p"]
h_T_initial  = d["h_T"]

if year != 1970
    d_aggA    = load(joinpath(param_dir, "aggA_1970.jld"))
    aggA_1970 = d_aggA["aggA_1970"]
end

###############################################################################
#  ABILITY DISTRIBUTION AND GRID
###############################################################################

dist = Frechet(theta, 1)

quantile_top    = 1 - 1e-4
quantile_bottom = 1e-3
n_a = 120
n_H = 11
range_H = 0.25

i_grid = log10.(range(1e3^(quantile_bottom), 1e3^(quantile_top), n_a)) ./ 3
b_grid = quantile.(dist, i_grid)

if b_grid != a_grid_init
    spl_a_grid = Spline1D(range(1, length(a_grid_init), length(a_grid_init)), a_grid_init)
    a_grid_upd = spl_a_grid.(range(1, length(a_grid_init), n_a))
    cc     = 1
    c_grid = cc .* b_grid + (1 - cc) .* a_grid_upd
    if n_a != length(a_grid_init)
        h_T_initial_upd = Array{Float64,3}(undef, n_a, n_H, n_G)
        a_T_thresh_upd  = Array{Float64,4}(undef, n_a, n_H, n_G, n_O)
        for iH = 1:n_H, iG = 1:n_G
            spl_h_T = Spline1D(a_grid_init, h_T_initial[:, iH, iG])
            h_T_initial_upd[:, iH, iG] = spl_h_T.(c_grid)
            for iO = 1:n_O-1
                spl_a_T = Spline1D(a_grid_init, a_T_thresh_init[:, iH, iG, iO])
                a_T_thresh_upd[:, iH, iG, iO] = spl_a_T.(c_grid)
            end
        end
        h_T_initial     = h_T_initial_upd
        a_T_thresh_init = a_T_thresh_upd
    end
    a_grid = c_grid
else
    a_grid = copy(a_grid_init)
end

lowbnd = a_grid[1]
upbnd  = a_grid[end]

###############################################################################
#  GROUP MEASURES
###############################################################################

gm = zeros(n_G)
for iG in 1:n_G-1
    gm[iG] = M / (2 * n_G)
end
gm[end] = M / 2 - sum(gm)

###############################################################################
#  ONE-TIME CALIBRATION: a_by_occ (men) AND τ_w_opt (women)
#
#  These depend only on the ability distribution and occupation shares, not
#  on λf or κ, so we compute them once before the root-finding loop.
###############################################################################

τ_e_base = fill(0.0, n_O - 1, n_G)

τ_w_men_raw = zeros(n_O - 1)
τ_w_men     = ones(n_O - 1) .- λm .* (ones(n_O - 1) .- τ_w_men_raw)

marg_cal = Array{Float64,3}(undef, n_a, n_O - 1, n_G)

function share_fn(occ_id, iG, a_by_occ_loc, τ_w_loc, τ_e_loc, marg_buf)
    a_by_occ_n = a_by_occ_loc ./ a_by_occ_loc[end]
    for ia in 1:n_a
        B_tmp = (((ones(n_O - 2) .- τ_w_loc[occ_id]) ./
                  (ones(n_O - 2) .- τ_w_loc[1:end.!=occ_id])) .*
                 (((ones(n_O - 2) .+ τ_e_loc[occ_id]) ./
                   (ones(n_O - 2) .+ τ_e_loc[1:end.!=occ_id])) .^ (-η)) .*
                 (a_by_occ_n[occ_id] ./ a_by_occ_n[1:end.!=occ_id]))
        A_tmp = B_tmp .^ (1 / α)
        marg_buf[ia, occ_id, iG] = prod(cdf.(dist, A_tmp .* a_grid[ia]))
    end
    spl_m = Spline1D(a_grid, marg_buf[:, occ_id, iG])
    out1 = quadgk(aa -> spl_m(aa) * pdf(dist, aa), lowbnd, upbnd)[1]
    out2 = marg_buf[:, occ_id, iG]
    return out1, out2
end

function calibrate_A_fn(x)
    share_occ_model = Vector{Float64}(undef, n_O - 1)
    for iO = 1:(n_O-1)
        share_occ_model[iO] = share_fn(iO, 2, x, τ_w_men, τ_e_base[:, 2], marg_cal)[1]
    end
    share_occ_model ./= sum(share_occ_model)
    return sum((((share_occ_model .- share_occ_data[:, 2]) ./ share_occ_data[:, 2]) .* 1e3) .^ 4)
end

println("═══════════════════════════════════════════════════════════════")
println("  Calibrating a_by_occ to match men's occupation shares...")
println("═══════════════════════════════════════════════════════════════")
res_A = optimize(calibrate_A_fn, a_by_occ_initial ./ a_by_occ_initial[end],
    LBFGS(),
    Optim.Options(show_trace=false, iterations=200000, f_abstol=1e-12, f_reltol=1e-12))
a_by_occ_cal = Optim.minimizer(res_A)
println("  Residual: ", Optim.minimum(res_A))

function calibrate_τ_fn(x)
    share_occ_model = Vector{Float64}(undef, n_O - 1)
    for iO = 1:(n_O-1)
        share_occ_model[iO] = share_fn(iO, 1, a_by_occ_cal, x, τ_e_base[:, 1], marg_cal)[1]
    end
    share_occ_model ./= sum(share_occ_model)
    return sum((((share_occ_model .- share_occ_data[:, 1]) ./ share_occ_data[:, 1]) .* 1e3) .^ 4)
end

println("  Calibrating τ_w_opt to match women's occupation shares...")
res0 = optimize(calibrate_τ_fn, zeros(n_O - 1),
                 Optim.Options(show_trace=false, iterations=200000, f_abstol=1e-12, f_reltol=1e-12))
res_τ = optimize(calibrate_τ_fn, Optim.minimizer(res0), LBFGS(),
                 Optim.Options(show_trace=false, iterations=200000, f_abstol=1e-12, f_reltol=1e-12))
τ_w_opt_cal = Optim.minimizer(res_τ)
println("  Residual: ", Optim.minimum(res_τ))

# Store marginal distributions for both groups (used inside the model solver)
marginal_cal = Array{Float64,3}(undef, n_a, n_O - 1, n_G)
for iO = 1:(n_O-1)
    _, marginal_cal[:, iO, 2] = share_fn(iO, 2, a_by_occ_cal, τ_w_men, τ_e_base[:, 2], marg_cal)
end
for iO = 1:(n_O-1)
    _, marginal_cal[:, iO, 1] = share_fn(iO, 1, a_by_occ_cal, τ_w_opt_cal, τ_e_base[:, 1], marg_cal)
end

println("  One-time calibration complete.")
println("═══════════════════════════════════════════════════════════════")

###############################################################################
#  TEACHER WAGE PROFILE
###############################################################################

ω_fn(κ_val, h)     = κ_val * h^γ
der_ω_fn(κ_val, h) = κ_val * γ * h^(γ - 1)

###############################################################################
#  MONOTONE INVERSE SPLINE
#
#  Build Spline1D(thresh_vals, a_grid) safely by sorting and deduplicating.
#  Required because at extreme parameter values the ability thresholds may
#  not be strictly increasing across the ability grid.
###############################################################################

function monotone_inverse_spline(thresh_vals::Vector{Float64},
                                  a_vals::Vector{Float64})
    perm  = sortperm(thresh_vals)
    x_srt = thresh_vals[perm]
    y_srt = a_vals[perm]
    mask  = [true; diff(x_srt) .> 1e-14]
    return Spline1D(x_srt[mask], y_srt[mask])
end

###############################################################################
#  CORE MODEL SOLVER
#
#  compute_teaching_shares(λf_val, κ_val; tol_mode, verbose)
#
#  Solves for the steady-state teaching shares of women and men given
#  parameters (λf, κ). The tol_mode argument controls the fixed-point
#  tolerances: :loose for Phase 1, :tight for Phase 2.
###############################################################################

function compute_teaching_shares(λf_val::Float64, κ_val::Float64;
                                  tol_mode::Symbol = :tight,
                                  verbose::Bool    = false)
    # Select tolerances based on mode
    if tol_mode == :loose
        tolHH_loc = phase1_tolHH
        tolT_loc  = phase1_tolT
        tolG_loc  = phase1_tolG
        ν_tax_loc = phase1_ν_tax
        ν_H_loc   = phase1_ν_H
    else  # :tight
        tolHH_loc = phase2_tolHH
        tolT_loc  = phase2_tolT
        tolG_loc  = phase2_tolG
        ν_tax_loc = phase2_ν_tax
        ν_H_loc   = phase2_ν_H
    end

    # Build τ_w from λf and calibrated τ_w_opt
    τ_w_loc = zeros(n_O - 1, n_G)
    τ_w_loc[:, 2] = copy(τ_w_men)
    τ_w_loc[:, 1] = ones(n_O - 1) .- λf_val .* (ones(n_O - 1) .- τ_w_opt_cal)
    τ_e_loc = copy(τ_e_base)

    a_by_occ_w = copy(a_by_occ_cal)
    s_O = μ * ϕ / (μ * ϕ + 1 - η)

    iHH_loc   = convert(Int, ceil(n_H / 2))
    HH_fp_loc = HH_fp_init
    H_grid_loc = copy(H_grid_init)
    t_loc      = copy(t_init)
    h_T_loc    = copy(h_T_initial)
    s_T_loc    = similar(h_T_loc)
    e_T_loc    = similar(h_T_loc)
    ω_loc      = similar(h_T_loc)
    der_ω_loc  = similar(h_T_loc)

    a_O_thresh_loc = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)
    h_O_loc        = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)

    H_O_0_loc   = Array{Float64,3}(undef, n_O - 1, n_G, n_H)
    E_O_loc     = Array{Float64,3}(undef, n_O - 1, n_H, n_G)
    Y_O_loc     = Array{Float64,3}(undef, n_O - 1, n_H, n_G)
    HH_T_0_loc  = Array{Float64,2}(undef, n_G, n_H)
    E_T_loc     = Array{Float64,2}(undef, n_H, n_G)
    sum_E_O_loc = Array{Float64,2}(undef, n_H, n_G)
    sum_Y_O_loc = Array{Float64,2}(undef, n_H, n_G)
    mass_T_loc  = Array{Float64,2}(undef, n_H, n_G)
    mass_O_loc  = Array{Float64,3}(undef, n_H, n_G, n_O - 1)

    f1_loc = Array{Float64,2}(undef, n_O - 1, n_G)
    f2_loc = Array{Float64,2}(undef, n_H, n_G)
    f3_loc = Array{Float64,2}(undef, n_O - 1, n_G)
    f4_loc = Array{Float64,2}(undef, n_H, n_G)

    f_1_T_loc     = Array{Float64,2}(undef, n_a, n_G)
    f_1_T_tmp_loc = Array{Float64,3}(undef, n_a, n_G, n_O - 1)
    f_1_O_loc     = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)

    convHH      = 1.0
    hh_iter     = 0
    max_hh_iter = 80

    while convHH > tolHH_loc && hh_iter < max_hh_iter
        hh_iter += 1
        H_grid_loc[iHH_loc] = HH_fp_loc

        if verbose
            println("    [HH iter $hh_iter] HH_fp = ",
                    round(HH_fp_loc; digits=6))
        end

        # ── Growth loop ──────────────────────────────────────────────
        convG = 1.0
        iterG = 1
        while convG > tolG_loc && iterG < maxiterG

            # ── Tax-rate loop ────────────────────────────────────────
            convT = 1.0
            iterT = 1
            while convT > tolT_loc && iterT < maxiterT
                for iG in 1:n_G
                    for ia in 1:n_a
                        a = a_grid[ia]
                        fn_s_T(k) = μ * ϕ * der_ω_fn(κ_val, k) * k /
                                    ((μ * ϕ - η) * der_ω_fn(κ_val, k) * k +
                                     ω_fn(κ_val, k))
                        fn_h_T(k) = (η^η * (1 - t_loc[iHH_loc, 1])^η *
                                     der_ω_fn(κ_val, k)^η * a^α *
                                     fn_s_T(k)^ϕ *
                                     (2 * HH_fp_loc / M)^σ)^(1 / (1 - η)) - k
                        hh_T = find_zero(fn_h_T, h_T_loc[ia, iHH_loc, iG])
                        h_T_loc[ia, iHH_loc, iG]   = hh_T
                        s_T_loc[ia, iHH_loc, iG]   = fn_s_T(hh_T)
                        e_T_loc[ia, iHH_loc, iG]   = η * (1 - t_loc[iHH_loc, 1]) *
                                                       der_ω_fn(κ_val, hh_T) * hh_T
                        ω_loc[ia, iHH_loc, iG]     = ω_fn(κ_val, hh_T)
                        der_ω_loc[ia, iHH_loc, iG] = der_ω_fn(κ_val, hh_T)

                        for iO in 1:n_O-1
                            a_O_thresh_loc[ia, iHH_loc, iG, iO] =
                                ((1 + τ_e_loc[iO, iG]) / (1 - τ_w_loc[iO, iG]) *
                                 (1 + τ_e_loc[iO, iG])^(-(1 - η)) *
                                 ((1 - s_T_loc[ia, iHH_loc, iG]) / (1 - s_O))^((1 - η) / μ) *
                                 (s_T_loc[ia, iHH_loc, iG] / s_O)^ϕ *
                                 (der_ω_loc[ia, iHH_loc, iG] / a_by_occ_w[iO]) *
                                 ((ω_loc[ia, iHH_loc, iG] /
                                   (der_ω_loc[ia, iHH_loc, iG] *
                                    h_T_loc[ia, iHH_loc, iG]) - η) /
                                  (1 - η))^(1 - η))^(1 / α) * a
                            h_O_loc[ia, iHH_loc, iG, iO] =
                                (η^η * (1 - t_loc[iHH_loc, 1])^η *
                                 (1 - τ_w_loc[iO, iG])^η /
                                 (1 + τ_e_loc[iO, iG])^η *
                                 a_by_occ_w[iO]^η * a_grid[ia]^α *
                                 s_O^ϕ * (2 * HH_fp_loc / M)^σ)^(1 / (1 - η))
                        end
                    end

                    spl_s  = Spline1D(a_grid, s_T_loc[:, iHH_loc, iG])
                    spl_dw = Spline1D(a_grid,
                                      der_ω_fn.(κ_val, h_T_loc[:, iHH_loc, iG]))

                    for iO in 1:n_O-1
                        spl_marg_l = Spline1D(a_grid, marginal_cal[:, iO, iG])
                        spl_inv    = monotone_inverse_spline(
                            a_O_thresh_loc[:, iHH_loc, iG, iO], a_grid)

                        f3_loc[iO, iG] = quadgk(
                            aa -> spl_marg_l(aa) * pdf(dist, aa) *
                                  cdf(dist, spl_inv(aa)) *
                                  aa^(α / (1 - η)),
                            lowbnd, upbnd)[1]
                        f1_loc[iO, iG] = quadgk(
                            aa -> spl_marg_l(aa) * pdf(dist, aa) *
                                  cdf(dist, spl_inv(aa)),
                            lowbnd, upbnd)[1]

                        H_O_0_loc[iO, iG, iHH_loc] =
                            (((1 - t_loc[iHH_loc, 1]) *
                              (1 - τ_w_loc[iO, iG]) /
                              (1 + τ_e_loc[iO, iG]))^η * η^η *
                             (2 * HH_fp_loc / M)^σ * s_O^ϕ *
                             a_by_occ_w[iO]^η)^(1 / (1 - η)) *
                            f3_loc[iO, iG]

                        E_O_loc[iO, iHH_loc, iG] =
                            a_by_occ_w[iO] * (1 - τ_w_loc[iO, iG]) *
                            H_O_0_loc[iO, iG, iHH_loc]
                        Y_O_loc[iO, iHH_loc, iG] =
                            a_by_occ_w[iO] * H_O_0_loc[iO, iG, iHH_loc]

                        f_1_O_loc[:, iHH_loc, iG, iO] =
                            marginal_cal[:, iO, iG] .*
                            cdf.(dist, spl_inv(a_grid))
                        spl_f1O = Spline1D(a_grid,
                                           f_1_O_loc[:, iHH_loc, iG, iO])
                        mass_O_loc[iHH_loc, iG, iO] = quadgk(
                            aa -> pdf(dist, aa) * spl_f1O(aa),
                            lowbnd, upbnd)[1]

                        for ia in 1:n_a
                            f_1_T_tmp_loc[ia, iG, iO] = max(
                                quadgk(aa -> spl_marg_l(aa) * pdf(dist, aa),
                                       lowbnd,
                                       a_O_thresh_loc[ia, iHH_loc, iG, iO])[1],
                                0.0)
                        end
                    end

                    for ia in 1:n_a
                        f_1_T_loc[ia, iG] = sum(f_1_T_tmp_loc[ia, iG, :])
                    end

                    spl_f1T = Spline1D(a_grid, f_1_T_loc[:, iG])
                    spl_wa  = Spline1D(a_grid, ω_loc[:, iHH_loc, iG])

                    f2_loc[iHH_loc, iG] = quadgk(
                        aa -> pdf(dist, aa) * spl_f1T(aa) *
                              aa^(α * β / σ / (1 - η)) *
                              spl_s(aa)^(ϕ * β / σ / (1 - η)) *
                              spl_dw(aa)^(η * β / σ / (1 - η)),
                        lowbnd, upbnd)[1]
                    f4_loc[iHH_loc, iG] = quadgk(
                        aa -> pdf(dist, aa) * spl_f1T(aa) * spl_wa(aa),
                        lowbnd, upbnd)[1]

                    mass_T_loc[iHH_loc, iG] = quadgk(
                        aa -> pdf(dist, aa) * spl_f1T(aa),
                        lowbnd, upbnd)[1]

                    HH_T_0_loc[iG, iHH_loc] =
                        ((1 - t_loc[iHH_loc, 1])^η * η^η *
                         (2 * HH_fp_loc / M)^σ)^(β / σ / (1 - η)) *
                        f2_loc[iHH_loc, iG]
                    E_T_loc[iHH_loc, iG] = f4_loc[iHH_loc, iG]
                end

                H_grid_loc[iHH_loc] = sum(HH_T_0_loc[:, iHH_loc] .* gm)
                for iG in 1:n_G
                    sum_E_O_loc[iHH_loc, iG] = sum(E_O_loc[:, iHH_loc, iG])
                    sum_Y_O_loc[iHH_loc, iG] = sum(Y_O_loc[:, iHH_loc, iG])
                end

                t_new = sum(E_T_loc[iHH_loc, :] .* gm) /
                        (sum(E_T_loc[iHH_loc, :] .* gm) +
                         sum(sum_E_O_loc[iHH_loc, :] .* gm))
                convT = abs(t_new - t_loc[iHH_loc, 1])
                t_loc[iHH_loc, 1] = (1 - ν_tax_loc) * t_loc[iHH_loc, 1] +
                                     ν_tax_loc * t_new
                t_loc[iHH_loc, 2] = t_new
                iterT += 1
            end  # tax loop

            # ── Aggregate productivity matching ──────────────────────
            if year != 1970
                aggA_loc = zeros(2)
                aggA_loc[1] = sum(sum_Y_O_loc[iHH_loc, :] .* gm) /
                              sum(gm[1] .* mass_O_loc[iHH_loc, 1, :] .+
                                  gm[2] .* mass_O_loc[iHH_loc, 2, :])
                aggA_loc[2] = sum((sum_Y_O_loc[iHH_loc, :] .-
                                   Y_O_loc[n_O-1, iHH_loc, :]) .* gm) /
                              sum(gm[1] .* mass_O_loc[iHH_loc, 1, 1:n_O-2] .+
                                  gm[2] .* mass_O_loc[iHH_loc, 2, 1:n_O-2])
                delta = aggA_1970[A_idx] * growth / aggA_loc[A_idx]
                convG = abs(aggA_1970[A_idx] * growth - aggA_loc[A_idx])
                a_by_occ_w = a_by_occ_w .* delta
            else
                convG = 0.0
            end
            iterG += 1
        end  # growth loop

        convHH = abs(log(H_grid_loc[iHH_loc] / HH_fp_loc))
        HH_fp_loc = (1 - ν_H_loc) * HH_fp_loc + ν_H_loc * H_grid_loc[iHH_loc]

        if verbose
            println("    [HH iter $hh_iter] convHH = ",
                    round(convHH; digits=8),
                    ", mass_T = (",
                    round(mass_T_loc[iHH_loc, 1]; digits=5), ", ",
                    round(mass_T_loc[iHH_loc, 2]; digits=5), ")")
        end
    end  # HH loop

    return mass_T_loc[iHH_loc, 1], mass_T_loc[iHH_loc, 2]
end

###############################################################################
#  PHASE 1 HELPERS: 1-D BISECTION SUB-PROBLEMS
###############################################################################

"""
    bracket_and_solve_κ(λf_fixed, κ_lo_b, κ_hi_b; ...)

Find κ such that mass_T_male(λf_fixed, κ) = target_mass_T_male.
Returns (κ_star, model_female_share, model_male_share).
"""
function bracket_and_solve_κ(λf_fixed::Float64, κ_lo_b::Float64, κ_hi_b::Float64;
                              atol::Float64 = 1e-4,
                              maxiter::Int  = 25,
                              verbose::Bool = false)
    eval_count = Ref(0)
    last_mf    = Ref(0.0)
    last_mm    = Ref(0.0)

    function residual(κ)
        eval_count[] += 1
        local mf, mm
        try
            mf, mm = compute_teaching_shares(λf_fixed, κ; tol_mode=:loose)
        catch e
            @warn "Model failed at κ=$κ (λf=$λf_fixed): $e"
            mid = (κ_lo_b + κ_hi_b) / 2
            return κ < mid ? -1.0 : 1.0
        end
        last_mf[] = mf;  last_mm[] = mm
        r = mm - target_mass_T_male
        if verbose
            println("      κ-solve [$(eval_count[])]  κ=", round(κ; digits=5),
                    "  mm=", round(mm; digits=5),
                    "  residual=", round(r; digits=6))
        end
        return r
    end

    κ_star = find_zero(residual, (κ_lo_b, κ_hi_b), Bisection();
                        atol=atol, maxevals=maxiter)
    return κ_star, last_mf[], last_mm[]
end

"""
    bracket_and_solve_λf(κ_fixed, λf_lo_b, λf_hi_b; ...)

Find λf such that mass_T_female(λf, κ_fixed) = target_mass_T_female.
Note: female share is DECREASING in λf, so the residual changes sign
from positive (low λf) to negative (high λf).
Returns (λf_star, model_female_share, model_male_share).
"""
function bracket_and_solve_λf(κ_fixed::Float64, λf_lo_b::Float64, λf_hi_b::Float64;
                               atol::Float64 = 1e-4,
                               maxiter::Int  = 25,
                               verbose::Bool = false)
    eval_count = Ref(0)
    last_mf    = Ref(0.0)
    last_mm    = Ref(0.0)

    function residual(λf)
        eval_count[] += 1
        mf, mm = compute_teaching_shares(λf, κ_fixed; tol_mode=:loose)
        last_mf[] = mf
        last_mm[] = mm
        r = mf - target_mass_T_female
        if verbose
            println("      λf-solve [$(eval_count[])]  λf=", round(λf; digits=5),
                    "  mf=", round(mf; digits=5),
                    "  residual=", round(r; digits=6))
        end
        return r
    end

    λf_star = find_zero(residual, (λf_lo_b, λf_hi_b), Bisection();
                         atol=atol, maxevals=maxiter)
    return λf_star, last_mf[], last_mm[]
end

###############################################################################
#  PHASE 1 — FAST GAUSS-SEIDEL (loose tolerances)
###############################################################################

function run_phase1(;verbose::Bool = true)
    println("\n═══════════════════════════════════════════════════════════════")
    println("  PHASE 1: Fast Gauss-Seidel root-finding (loose tolerances)")
    println("  Targets:  mass_T_female = $target_mass_T_female")
    println("            mass_T_male   = $target_mass_T_male")
    println("  Initial:  λf = $λf_init,  κ = $κ_init")
    println("═══════════════════════════════════════════════════════════════\n")

    λf_cur = λf_init
    κ_cur  = κ_init

    for outer = 1:phase1_outer_maxiter
        println("  ── Outer sweep $outer ──")

        # (a) Fix λf, solve for κ to match male share
        println("    Solving for κ  (λf fixed at $(round(λf_cur; digits=5)))...")
        κ_cur, mf_a, mm_a = bracket_and_solve_κ(
            λf_cur, κ_lo, κ_hi;
            atol    = phase1_brent_atol,
            maxiter = phase1_brent_maxiter,
            verbose = verbose)
        println("    → κ = ", round(κ_cur; digits=6),
                "   (mf=", round(mf_a; digits=5),
                ", mm=", round(mm_a; digits=5), ")")

        # (b) Fix κ, solve for λf to match female share
        println("    Solving for λf  (κ fixed at $(round(κ_cur; digits=5)))...")
        λf_cur, mf_b, mm_b = bracket_and_solve_λf(
            κ_cur, λf_lo, λf_hi;
            atol    = phase1_brent_atol,
            maxiter = phase1_brent_maxiter,
            verbose = verbose)
        println("    → λf = ", round(λf_cur; digits=6),
                "   (mf=", round(mf_b; digits=5),
                ", mm=", round(mm_b; digits=5), ")")

        # Check convergence on BOTH moments
        rel_err_f = abs(mf_b - target_mass_T_female) / target_mass_T_female
        rel_err_m = abs(mm_b - target_mass_T_male)   / target_mass_T_male
        println("    Relative errors:  female = ", round(rel_err_f; digits=6),
                ",  male = ", round(rel_err_m; digits=6))

        if rel_err_f < phase1_outer_tol && rel_err_m < phase1_outer_tol
            println("\n  ✓ Phase 1 converged at outer sweep $outer.")
            break
        end
        if outer == phase1_outer_maxiter
            println("\n  ⚠ Phase 1 reached max outer iterations. Proceeding to Phase 2.")
        end
    end

    return λf_cur, κ_cur
end

###############################################################################
#  PHASE 2 — BROYDEN'S METHOD (tight tolerances)
#
#  2D quasi-Newton root-finding for (λf, κ).
#
#  Cost breakdown:
#    • Initial Jacobian: 3 model solves (base + 2 finite-difference perturbations)
#    • Each Broyden step: 1 model solve + O(1) linear algebra
#    • Typical convergence: 3-5 iterations
#    • Total: ~6-8 model solves
#
#  The Jacobian J is the 2×2 matrix:
#    J = [ ∂F_female/∂λf   ∂F_female/∂κ ]
#        [ ∂F_male/∂λf     ∂F_male/∂κ   ]
#
#  where F_female = mass_T_female(λf,κ) - target_female
#        F_male   = mass_T_male(λf,κ)   - target_male
#
#  After the initial finite-difference Jacobian, subsequent iterations use
#  the Broyden rank-1 update:  J_new = J_old + (ΔF - J_old·Δx)·Δxᵀ / (Δx·Δx)
###############################################################################

function eval_residual(λf_val::Float64, κ_val::Float64;
                       verbose::Bool = false)
    mf, mm = compute_teaching_shares(λf_val, κ_val;
                                      tol_mode=:tight, verbose=verbose)
    F = [mf - target_mass_T_female,
         mm - target_mass_T_male]
    return F, mf, mm
end

function clamp_step(λf_cur, κ_cur, δλf, δκ)
    δλf = clamp(δλf, -phase2_max_step_λf, phase2_max_step_λf)
    δκ  = clamp(δκ,  -phase2_max_step_κ,  phase2_max_step_κ)
    λf_new = clamp(λf_cur + δλf, λf_lo, λf_hi)
    κ_new  = clamp(κ_cur  + δκ,  κ_lo,  κ_hi)
    return λf_new - λf_cur, κ_new - κ_cur
end

function run_phase2_broyden(λf_start::Float64, κ_start::Float64;
                            verbose::Bool = true)
    println("\n═══════════════════════════════════════════════════════════════")
    println("  PHASE 2: Broyden's method (tight tolerances)")
    println("  Starting point:  λf = ", round(λf_start; digits=6),
            ",  κ = ", round(κ_start; digits=6))
    println("  Targets:  mass_T_female = $target_mass_T_female")
    println("            mass_T_male   = $target_mass_T_male")
    println("  Fixed-point tolerances:  tolHH=$phase2_tolHH, tolT=$phase2_tolT, tolG=$phase2_tolG")
    println("  Damping:  ν_tax=$phase2_ν_tax, ν_H=$phase2_ν_H")
    println("═══════════════════════════════════════════════════════════════\n")

    n_evals = Ref(0)

    function counted_residual(λf, κ)
        n_evals[] += 1
        F, mf, mm = eval_residual(λf, κ)
        return F, mf, mm
    end

    # ── Step 1: Evaluate at starting point ───────────────────────────────
    println("  [1/3] Evaluating residual at starting point...")
    F0, mf0, mm0 = counted_residual(λf_start, κ_start)
    println("    F = [", round(F0[1]; digits=6), ", ", round(F0[2]; digits=6), "]")
    println("    (mf=", round(mf0; digits=5), ", mm=", round(mm0; digits=5), ")")

    # ── Step 2: Build initial Jacobian via finite differences ────────────
    println("\n  [2/3] Building initial Jacobian (2 extra model solves)...")

    λf_pert = λf_start + phase2_fd_step_λf
    F_λf, _, _ = counted_residual(λf_pert, κ_start)
    dF_dλf = (F_λf .- F0) ./ phase2_fd_step_λf

    κ_pert = κ_start + phase2_fd_step_κ
    F_κ, _, _ = counted_residual(λf_start, κ_pert)
    dF_dκ = (F_κ .- F0) ./ phase2_fd_step_κ

    J = hcat(dF_dλf, dF_dκ)  # 2×2 Jacobian

    println("    Initial Jacobian:")
    println("      ∂F_female/∂λf = ", round(J[1,1]; digits=4),
            "   ∂F_female/∂κ = ", round(J[1,2]; digits=4))
    println("      ∂F_male/∂λf   = ", round(J[2,1]; digits=4),
            "   ∂F_male/∂κ   = ", round(J[2,2]; digits=4))

    if J[1,1] > 0
        @warn "Unexpected: ∂F_female/∂λf > 0 (expected negative). Jacobian may be noisy."
    end
    if J[2,2] < 0
        @warn "Unexpected: ∂F_male/∂κ < 0 (expected positive). Jacobian may be noisy."
    end

    # ── Step 3: Broyden iterations ───────────────────────────────────────
    println("\n  [3/3] Broyden iterations...")

    λf_cur = λf_start
    κ_cur  = κ_start
    F_cur  = copy(F0)

    for iter = 1:phase2_max_newton_iter
        # Newton step: J · δx = -F
        δx = -(J \ F_cur)

        # Clamp step to stay within bounds and max step sizes
        δλf_clamped, δκ_clamped = clamp_step(λf_cur, κ_cur, δx[1], δx[2])
        δx_clamped = [δλf_clamped, δκ_clamped]

        # ── Backtracking line search on ‖F‖² ────────────────────────
        merit_cur = dot(F_cur, F_cur)
        step_factor = 1.0
        λf_try = λf_cur + δx_clamped[1]
        κ_try  = κ_cur  + δx_clamped[2]
        F_try, mf_try, mm_try = counted_residual(λf_try, κ_try)
        merit_try = dot(F_try, F_try)

        while merit_try > merit_cur * (1.0 - phase2_armijo_c * step_factor) &&
              step_factor > phase2_min_step_factor
            step_factor *= 0.5
            λf_try = λf_cur + step_factor * δx_clamped[1]
            κ_try  = κ_cur  + step_factor * δx_clamped[2]
            F_try, mf_try, mm_try = counted_residual(λf_try, κ_try)
            merit_try = dot(F_try, F_try)
            if verbose
                println("      Line search: step_factor=", round(step_factor; digits=3),
                        "  ‖F‖²=", round(merit_try; digits=10))
            end
        end

        # Actual step taken
        Δx = [λf_try - λf_cur, κ_try - κ_cur]
        ΔF = F_try .- F_cur

        # Broyden rank-1 update: J_new = J + (ΔF - J·Δx)·Δxᵀ / (Δx·Δx)
        denom = dot(Δx, Δx)
        if denom > 1e-30
            J .+= ((ΔF .- J * Δx) .* Δx') ./ denom
        end

        # Update state
        λf_cur = λf_try
        κ_cur  = κ_try
        F_cur  = F_try

        # Convergence check
        rel_err_f = abs(mf_try - target_mass_T_female) / target_mass_T_female
        rel_err_m = abs(mm_try - target_mass_T_male)   / target_mass_T_male

        println("    Iter $iter:  λf=", round(λf_cur; digits=6),
                "  κ=", round(κ_cur; digits=6),
                "  mf=", round(mf_try; digits=6),
                "  mm=", round(mm_try; digits=6),
                "  rel_err=(", round(rel_err_f; digits=5),
                ", ", round(rel_err_m; digits=5), ")",
                step_factor < 1.0 ? "  [step×$(round(step_factor; digits=2))]" : "")

        if rel_err_f < phase2_tol && rel_err_m < phase2_tol
            println("\n  ✓ Phase 2 converged in $iter Broyden iterations ",
                    "($(n_evals[]) total model solves).")
            return λf_cur, κ_cur
        end
    end

    println("\n  ⚠ Phase 2 reached max iterations ($(n_evals[]) model solves).")
    return λf_cur, κ_cur
end

###############################################################################
#  RUN CALIBRATION
###############################################################################

# Phase 1: fast, loose tolerances
λf_p1, κ_p1 = λf_init, κ_init
if skip_phase1
    println("\n═══════════════════════════════════════════════════════════════")
    println("  Skipping Phase 1; starting Phase 2 from initial guess.")
    println("  Initial guess:  λf = $λf_init,  κ = $κ_init")
    println("═══════════════════════════════════════════════════════════════\n")
else
    λf_p1, κ_p1 = run_phase1(verbose = verbose_rootfind)
end

# Phase 2: Broyden with tight tolerances
λf_p2, κ_p2 = run_phase2_broyden(λf_p1, κ_p1; verbose = verbose_rootfind)

###############################################################################
#  FINAL VERIFICATION
###############################################################################

println("\n═══════════════════════════════════════════════════════════════")
println("  Final verification (tight tolerances, verbose)...")
println("═══════════════════════════════════════════════════════════════")

mf_final, mm_final = compute_teaching_shares(λf_p2, κ_p2;
                                              tol_mode=:tight, verbose=true)

println("\n")
println("╔═════════════════════════════════════════════════════════════╗")
println("║                    CALIBRATION RESULT                      ║")
println("╠═════════════════════════════════════════════════════════════╣")
println("║  Environment:                                              ║")
println("║    year       = $year")
println("║    growth     = $(growth_pct * 100)%")
println("║    include_HP = $include_HP")
println("║    paramname  = $paramname")
println("╠═════════════════════════════════════════════════════════════╣")
println("║  Optimal parameters:                                       ║")
println("║    λf = ", round(λf_p2; digits=6))
println("║    κ  = ", round(κ_p2; digits=6))
println("╠═════════════════════════════════════════════════════════════╣")
println("║  Model moments vs targets:                                 ║")
println("║    mass_T_female:  model = ", round(mf_final; digits=6),
        "   target = ", target_mass_T_female)
println("║    mass_T_male:    model = ", round(mm_final; digits=6),
        "   target = ", target_mass_T_male)
println("║    rel error (f) = ",
        round(abs(mf_final - target_mass_T_female) / target_mass_T_female * 100;
              digits=4), "%")
println("║    rel error (m) = ",
        round(abs(mm_final - target_mass_T_male) / target_mass_T_male * 100;
              digits=4), "%")
println("╚═════════════════════════════════════════════════════════════╝")