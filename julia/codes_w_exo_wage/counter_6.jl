# counter_6.jl
# This Julia script computes counterfactual steady states for 1990 and 2010 where
# aggregate teacher human capital (H̃_T) is held fixed at its 1970 benchmark value.
# Everything else about the benchmark calibration (barriers, productivities, κ, etc.) is unchanged.
# This serves as a check of partial vs general equilibrium effects:
# by fixing the aggregate state variable, we isolate the direct (partial equilibrium) effects
# of parameter changes from the indirect (general equilibrium) feedback through H̃_T.

# Key changes relative to the main solution (frechet.jl):
# (1) Load HH_fp from the 1970 benchmark and hold it fixed across all years.
# (2) No fixed-point iteration over HH_fp — it is taken as exogenous.
# (3) Load calibrated a_by_occ, τ_w, κ, etc. from the benchmark for each year.
# (4) Recompute human capital investments, tax rate, wages, and distributions given fixed HH_fp.
# (5) Skip transition dynamics — only compute counterfactual steady states.
# (6) Only run for 1990 and 2010 (1970 values are identical to the benchmark by construction).

using Distributions, Dierckx, QuadGK, JLD, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim, Roots, PyCall, Random, Plots

# Change directory, if necessary:
# cd("./GitHub/teachers/julia/codes_w_exo_wage")

# Select type of parameterization (from following list):
# (1) benchmark - no growth, no HP in TFP growth
# (2) benchmark_w_growth
# (3) benchmark_w_HP
# (4) benchmark_w_growth_HP
# (5) low_beta
# (6) high_beta 
# (7) high_beta_high_sigma
# (8) low_beta_low_sigma
paramname = "benchmark"

# No transition dynamics in counterfactual (steady states only):
transition = 0

# Population size:
M = 1.0
# Groups (no discrimination / no barrier group is last element)
g = ["female", "male"]
# Import moments for calibration:
import XLSX
# Employment shares by occupation:
cd("..")
cd("..")
cd("./data/LaborMarketData")
# Employment shares:
xs1 = XLSX.readxlsx("wages_occ_shares_v2.xlsx")
tab1 = xs1["moments_shares"]
# Wage dispersion by occupation:
xs2 = XLSX.readxlsx("wages_occ_shares_v2.xlsx")
tab2 = xs2["90_10_hr_wages_weighted"]

# Reset working directory to folder with Julia script:
cd("..")
cd("..")
cd("./julia/codes_w_exo_wage")
# Labels for occupations:
occ_O = tab1["A30:A49"]
occ_T = tab1["A52"]
occ = [occ_O; occ_T]

# Innovation step size for updates:
ν = 1 # for tax rate that balances government's budget

n_G = length(g) # number of "groups"
n_O = length(occ) # number of occupations

gm = zeros(1, n_G) # measure of in individuals in groups 1,...,n_G
for iG in 1:n_G-1
    gm[iG] = M / (2 * n_G)
end
gm[end] = M / 2 - sum(gm)
gm = gm'

# Load 1970 benchmark HH_fp (the value we hold fixed across all years):
cd(string("./parameterization/", paramname))
d_1970 = load("previousParameterization1970.jld")
HH_fp_1970 = d_1970["HH_fp"]
cd("..")
cd("..")
println("Loaded 1970 benchmark HH_fp = ", HH_fp_1970)

# Loop over counterfactual years (1990 and 2010 only; 1970 is identical to benchmark by construction)
for year in [1990, 2010]

local fyear, fnameJLD, d, df, h_T_avg_ss, H_T_total

println("  COUNTERFACTUAL 6: year = ", year, " with HH_fp fixed at 1970 benchmark value")

global share_occ_data, w_90_10_data
share_occ_data = Array{Float64,2}(undef, length(occ) - 1, 2)
w_90_10_data = Array{Float64,1}(undef, 2)

# Load data for selected year:
if year == 1990
    share_occ_data[:, 1] = tab1["M30:M49"] # Census 1990 for NLSY79 (women)
    share_occ_data[:, 2] = tab1["E30:E49"] # Census 1990 for NLSY79 (men)
    w_90_10_data[1] = tab2["E50"] # Dispersion in non-teaching occupations
    w_90_10_data[2] = tab2["E51"] # Dispersion in teaching
elseif year == 2010
    share_occ_data[:, 1] = tab1["O30:O49"] # ACS 2009-2013 for NLSY97 (women)
    share_occ_data[:, 2] = tab1["G30:G49"] # ACS 2009-2013 1990 for NLSY97 (men)
    w_90_10_data[1] = tab2["G50"] # Dispersion in non-teaching occupations
    w_90_10_data[2] = tab2["G51"] # Dispersion in teaching
end


# Load calibrated parameters from the main model solution for this year
cd(string("./parameterization/", paramname))
fyear = string(year)
fnameJLD = string("previousParameterization", fyear, ".jld")
d = load(fnameJLD)
global a_by_occ, τ_w_initial, τ_e_initial, a_T_thresh, t, H_grid, H_O, HH_fp
global α, β, η, σ, μ, ϕ, γ, κ, theta, λf, λm, iHH, a_grid
global a_O_10p, a_O_90p, a_T_10p, a_T_90p, h_T_initial
# Load calibrated occupation-specific productivities from main solution:
a_by_occ = d["a_by_occ_level"] # uses the non-normalized version
τ_w_initial = d["τ_w_opt"]
τ_e_initial = d["τ_e"]
a_T_thresh = d["a_T_thresh"]
t = (d["t"])
H_grid = d["H_grid"]
H_O = d["H_O"]
HH_fp = d["HH_fp"] # benchmark value (will be overridden below)
α = d["α"]
β = d["β"]
η = d["η"]
σ = d["σ"]
μ = d["μ"]
ϕ = d["ϕ"]
γ = d["γ"]
κ = d["κ"]
theta = d["theta"]
λf = d["λf"]
λm = d["λm"]
iHH = d["iHH"]
a_grid = d["a_grid"]
a_O_10p = d["a_O_10p"]
a_O_90p = d["a_O_90p"]
a_T_10p = d["a_T_10p"]
a_T_90p = d["a_T_90p"]
h_T_initial = d["h_T"]
cd("..")
cd("..")

# Override HH_fp with the 1970 benchmark value:
println("Benchmark HH_fp for year ", year, " = ", HH_fp)
HH_fp = HH_fp_1970
println("Setting HH_fp to 1970 value = ", HH_fp)

# Distribution of abilities:
global dist
dist = Frechet(theta, 1)

# Vector of labor market discrimination in each occupation (relative to 'T')
# Change: Keep all barriers at their benchmark-calibrated values
global τ_w, τ_e, τ_w_opt
τ_w = zeros(n_O - 1, n_G)
# Vector of education barriers in each occupation (relative to 'T')
τ_e = zeros(n_O - 1, n_G)

τ_w_opt = zeros(n_O - 1, n_G)

τ_w[:, 2] = fill!(τ_w[:, 2], 0)
τ_w[:, 2] = ones(n_O - 1) .- λm * (ones(n_O - 1) .- τ_w[:, 2])
τ_e = fill!(τ_e, 0.0)

# Keep benchmark barriers for women (no change from calibration):
τ_w_opt[:, 1] = τ_w_initial[:, 1]
τ_w[:, 1] = ones(n_O - 1) .- λf * (ones(n_O - 1) .- τ_w_opt[:, 1])
println("Keeping benchmark barriers: λf = ", λf)

global quantile_top, quantile_bottom
quantile_top = 1 - 1e-4
quantile_bottom = 1e-3

# Reset grid length for individual abilities (if necessary):
global n_a
n_a = 120
# Set up the grid for aggregate human capital in teaching:
global n_H
n_H = 11
# Set the percentage range above / below fixed point:
global range_H
range_H = 0.25
# Update the grid for abilities (if necessary):
global i_grid, b_grid
i_grid = log10.(range(1e3^(quantile_bottom), 1e3^(quantile_top), n_a)) ./ 3
b_grid = quantile.(dist, i_grid)
if b_grid != a_grid
    spl_a_grid = Spline1D(range(1, length(a_grid), length(a_grid)), a_grid)
    a_grid_upd = spl_a_grid.(range(1, length(a_grid), n_a))
    cc = 1
    c_grid = cc .* b_grid + (1 - cc) .* a_grid_upd
    # Update size of loaded arrays, if necessary:
    if n_a != length(a_grid)
        h_T_initial_upd = Array{Float64,3}(undef, n_a, n_H, n_G)
        a_T_thresh_upd = Array{Float64,4}(undef, n_a, n_H, n_G, n_O)
        for iH = 1:n_H
            for iG = 1:n_G
                spl_h_T = Spline1D(a_grid, h_T_initial[:, iH, iG])
                h_T_initial_upd[:, iH, iG] = spl_h_T.(c_grid)
                for iO = 1:n_O-1
                    spl_a_T_thresh = Spline1D(a_grid, a_T_thresh[:, iH, iG, iO])
                    a_T_thresh_upd[:, iH, iG, iO] = spl_a_T_thresh.(c_grid)
                end
            end
        end
        h_T_initial = h_T_initial_upd
        a_T_thresh = a_T_thresh_upd
        a_grid = c_grid
    end
end

# Integral bounds:
global lowbnd, upbnd
lowbnd = a_grid[1] # set to smallest element in 'a_grid'
upbnd = a_grid[end] # set to largest element in 'a_grid'

# Given a productivity vector 'a_by_occ' with productivities in non-teaching occupations, the function 'share()' calculates share of people choosing occupation 'occ_id':
global marg
marg = Array{Float64,3}(undef, n_a, n_O - 1, n_G)
function share(occ_id, iG, a_by_occ, τ_w, τ_e)
    a_by_occ = a_by_occ ./ a_by_occ[end]
    for ia in 1:n_a
        B_tmp = (((ones(n_O - 2) .- τ_w[occ_id]) ./ (ones(n_O - 2) .- τ_w[1:end.!=occ_id])) .* (((ones(n_O - 2) .+ τ_e[occ_id]) ./ (ones(n_O - 2) .+ τ_e[1:end.!=occ_id])) .^ (-η)) .* (a_by_occ[occ_id] ./ a_by_occ[1:end.!=occ_id]))
        A_tmp = B_tmp .^ (1 / α)
        marg[ia, occ_id, iG] = prod(cdf.(dist, A_tmp .* a_grid[ia]))
    end
    spl_marg = Spline1D(a_grid, marg[:, occ_id, iG])
    # Function returns two output arguments:
    out1 = quadgk(aa -> spl_marg(aa) * pdf(dist, aa), lowbnd, upbnd)[1]
    out2 = marg[:, occ_id, iG]
    return out1, out2
end

# Skip calibration of a_by_occ and τ_w — use benchmark values

# Compute occupation shares and marginals using benchmark barriers and productivities
global share_occ, marginal
share_occ = Array{Float64,2}(undef, n_O - 1, n_G)
marginal = Array{Float64,3}(undef, n_a, n_O - 1, n_G)
for iG in 1:n_G
    for iO = 1:(n_O-1)
        share_occ[iO, iG] = share(iO, iG, a_by_occ, τ_w[:, iG], τ_e[:, iG])[1]
        # Fraction of workers with abilities 'a' in (non-teaching) occupation 'iO' who actually work in occupation 'iO':
        marginal[:, iO, iG] = share(iO, iG, a_by_occ, τ_w[:, iG], τ_e[:, iG])[2]
    end
    share_occ[:, iG] = share_occ[:, iG] ./ sum(share_occ[:, iG])
end

global a_T_thresh, a_O_thresh
a_T_thresh = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)
a_O_thresh = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)

global e_T, e_O
e_T = Array{Float64,3}(undef, n_a, n_H, n_G)
e_O = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)

global s_T
s_T = Array{Float64,3}(undef, n_a, n_H, n_G)

global h_T, h_T_avg, h_O, ω, w, ω_90_10, w_90_10, der_ω
h_T = Array{Float64,3}(undef, n_a, n_H, n_G)
h_T_avg = Array{Float64,2}(undef, n_H, n_G)
h_O = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)
ω = Array{Float64,3}(undef, n_a, n_H, n_G)
w = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)
ω_90_10 = Array{Float64,2}(undef, n_H, n_G + 1)
w_90_10 = Array{Float64,3}(undef, n_H, n_G + 1, n_O)
der_ω = Array{Float64,3}(undef, n_a, n_H, n_G)

global a_T_90_10, a_O_90_10
a_T_90_10 = Array{Float64,2}(undef, n_H, n_G)
a_O_90_10 = Array{Float64,3}(undef, n_H, n_G, n_O)

global spl_T, spl_O, spl_T_inv, spl_O_inv, spl_ω, spl_w, spl_h_O, spl_pdf_T
spl_T = Array{Spline1D,1}(undef, n_G)
spl_O = Array{Spline1D,1}(undef, n_G)
spl_T_inv = Array{Spline1D,1}(undef, n_G)
spl_O_inv = Array{Spline1D,1}(undef, n_G)
spl_ω = Array{Spline1D,1}(undef, n_H)
spl_w = Array{Spline1D,1}(undef, n_G)
spl_h_O = Array{Spline1D,1}(undef, n_G)
spl_pdf_T = Array{Spline1D,1}(undef, n_G)
global wb_T, wb_O, E_T, E_O, Y_T, Y_O, sum_E_O, sum_Y_O
wb_T = Array{Float64,2}(undef, n_H, n_G)
wb_O = Array{Float64,2}(undef, n_H, n_G)
E_T = Array{Float64,2}(undef, n_H, n_G)
E_O = Array{Float64,3}(undef, n_O - 1, n_H, n_G)
Y_T = Array{Float64,2}(undef, n_H, n_G)
Y_O = Array{Float64,3}(undef, n_O - 1, n_H, n_G)
sum_E_O = Array{Float64,2}(undef, n_H, n_G)
sum_Y_O = Array{Float64,2}(undef, n_H, n_G)

global HH_T_0, HH_T, H_O_0
HH_T_0 = Array{Float64,2}(undef, n_G, n_H)
HH_T = Array{Float64,1}(undef, n_H)
H_O_0 = Array{Float64,3}(undef, n_O - 1, n_G, n_H)
H_O = Array{Float64,2}(undef, n_O - 1, n_H)

global EN, EN2
EN = Array{Float64,1}(undef, n_G)
EN2 = Array{Float64,1}(undef, n_G)

global f_2, f_O
f_2 = Array{Float64,2}(undef, n_a, n_G)
f_O = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)

global spl_T_thresh
spl_T_thresh = Array{Spline1D,1}(undef, 1)

global f1, mass_O, f3, f_1_T, f_1_O, spl_f_1_T, spl_f_1_O
f1 = Array{Float64,2}(undef, n_O - 1, n_G)
mass_O = Array{Float64,3}(undef, n_H, n_G, n_O - 1)
f3 = Array{Float64,2}(undef, n_O - 1, n_G)
f_1_T = Array{Float64,2}(undef, n_a, n_G)
f_1_O = Array{Float64,4}(undef, n_a, n_H, n_G, n_O - 1)
spl_f_1_T = Array{Spline1D,2}(undef, n_H, n_G)
spl_f_1_O = Array{Spline1D,3}(undef, n_H, n_G, n_O - 1)
global f_T, spl_f_T, spl_f_O, mass_T, f2, f4, ratio, f_1_T_tmp
f_T = Array{Float64,3}(undef, n_a, n_H, n_G)
spl_f_T = Array{Spline1D,1}(undef, n_G)
spl_f_O = Array{Spline1D,3}(undef, n_H, n_G, n_O)
mass_T = Array{Float64,2}(undef, n_H, n_G)
f2 = Array{Float64,2}(undef, n_H, n_G)
f4 = Array{Float64,2}(undef, n_H, n_G)
ratio = Array{Float64,1}(undef, n_G)
f_1_T_tmp = Array{Float64,3}(undef, n_a, n_G, n_O - 1)

#######################################
# Exogenous wage profile for teachers #
#######################################
ω_fn(κ, h_T) = κ .* h_T .^ γ
# Derivative of wage profile:
der_ω_fn(κ, h_T) = κ .* γ .* h_T .^ (γ - 1)


###############################################################
# Compute equilibrium given fixed HH_fp (no HH fixed-point)  #
###############################################################

# Set the array index to the midpoint of H_grid:
global iHH
iHH = convert(Int, ceil(n_H / 2))
# Iteration and tolerance settings for the income tax rate fixed-point:
global tolT, maxiterT
tolT = 1e-7
maxiterT = 200
global aggA
aggA = zeros(2) # Aggregate productivity: aggA[1] = w/ home production, aggA[2] = w/o home production

# Time investment doesn't depend on any endogenous variables, only on parameters:
global s_O
s_O = μ * ϕ / (μ * ϕ + 1 - η)

# No HH fixed-point loop — HH_fp is held fixed at its 1970 value.
# Only the tax rate fixed-point remains.

H_grid[iHH] = HH_fp
println("HH_fp (fixed at 1970 value) = ", HH_fp)
println("t[iHH,1] = ", t[iHH, 1])
println("")

# Initiate 'while' loop ('not indexed') over the tax rate:
convT = 1
iterT = 1
while convT > tolT && iterT < maxiterT
    for iG in 1:n_G
        for ia in 1:n_a
            a = a_grid[ia]
            fn_s_T(k) = μ * ϕ * der_ω_fn(κ, k) * k / ((μ * ϕ - η) * der_ω_fn(κ, k) * k + ω_fn(κ, k))
            fn_h_T(k) = (η^η * (1 - t[iHH, 1])^η * der_ω_fn(κ, k)^η * a^α * fn_s_T(k)^ϕ * (2 * HH_fp / M)^σ)^(1 / (1 - η)) - k
            hh_T = find_zero(fn_h_T, h_T_initial[ia, iHH, iG])
            h_T[ia, iHH, iG] = hh_T
            s_T[ia, iHH, iG] = fn_s_T(hh_T)
            e_T[ia, iHH, iG] = η * (1 - t[iHH, 1]) * der_ω_fn(κ, hh_T) * hh_T
            ω[ia, iHH, iG] = ω_fn(κ, hh_T)
            der_ω[ia, iHH, iG] = der_ω_fn(κ, hh_T)
            for iO in 1:n_O-1
                # Occupational choice threshold:
                a_O_thresh[ia, iHH, iG, iO] = ((1 + τ_e[iO, iG]) / (1 - τ_w[iO, iG]) * (1 + τ_e[iO, iG])^(-(1 - η)) * ((1 - s_T[ia, iHH, iG]) / (1 - s_O))^((1 - η) / μ) * (s_T[ia, iHH, iG] / s_O)^ϕ * (der_ω[ia, iHH, iG] / a_by_occ[iO]) * ((ω[ia, iHH, iG] / (der_ω[ia, iHH, iG] * h_T[ia, iHH, iG]) - η) / (1 - η))^(1 - η))^(1 / α) * a
                h_O[ia, iHH, iG, iO] = (η^η * (1 - t[iHH, 1])^η * (1 - τ_w[iO, iG])^η / (1 + τ_e[iO, iG])^η * a_by_occ[iO]^η * a_grid[ia]^α * s_O^ϕ * (2 * HH_fp / M)^σ)^(1 / (1 - η))
                e_O[ia, iHH, iG, iO] = η * (1 - t[iHH, 1]) * (1 - τ_w[iO, iG]) / (1 + τ_e[iO, iG]) * a_by_occ[iO] * h_O[ia, iHH, iG, iO]
            end
        end
        spl_s = Spline1D(a_grid, s_T[:, iHH, iG])
        spl_dw = Spline1D(a_grid, der_ω_fn.(κ, h_T[:, iHH, iG]))

        for iO in 1:n_O-1
            spl_marg = Spline1D(a_grid, marginal[:, iO, iG])
            # Inverse occupational threshold 
            spl_inv = Spline1D(a_O_thresh[:, iHH, iG, iO], a_grid)
            f3[iO, iG] = quadgk(aa -> spl_marg(aa) * pdf(dist, aa) * cdf(dist, spl_inv(aa)) * aa^(α / (1 - η)), lowbnd, upbnd)[1]

            f1[iO, iG] = quadgk(aa -> spl_marg(aa) * pdf(dist, aa) * cdf(dist, spl_inv(aa)), lowbnd, upbnd)[1] # total mass of other using direct threshold
            # Compute H_O:
            H_O_0[iO, iG, iHH] = (((1 - t[iHH, 1]) * (1 - τ_w[iO, iG]) / (1 + τ_e[iO, iG]))^η * η^η * (2 * HH_fp / M)^σ * s_O^ϕ * a_by_occ[iO]^η)^(1 / (1 - η)) * f3[iO, iG]
            # Total earnings of others in each group:
            E_O[iO, iHH, iG] = a_by_occ[iO] * (1 - τ_w[iO, iG]) * H_O_0[iO, iG, iHH]
            # Total output of others in each group:
            Y_O[iO, iHH, iG] = a_by_occ[iO] * H_O_0[iO, iG, iHH]

            # Inverse occupational threshold between occupation 'iO' and teaching: 
            f_1_O[:, iHH, iG, iO] = marginal[:, iO, iG] .* cdf.(dist, spl_inv(a_grid))
            # Compute p.d.f.:
            spl_f_1_O[iHH, iG, iO] = Spline1D(a_grid, f_1_O[:, iHH, iG, iO])
            mass_O[iHH, iG, iO] = quadgk(aa -> pdf(dist, aa) * spl_f_1_O[iHH, iG, iO](aa), lowbnd, upbnd)[1]

            for ia in 1:n_a
                a = a_grid[ia]
                # Fraction of teachers occupation-by-occupation, given 'a' in teaching:
                f_1_T_tmp[ia, iG, iO] = maximum([quadgk(aa -> spl_marg(aa) * pdf(dist, aa), lowbnd, a_O_thresh[ia, iHH, iG, iO])[1], 0.0])
            end
        end

        for ia in 1:n_a
            # Fraction of teachers given a_T (across all occupations):
            f_1_T[ia, iG] = sum(f_1_T_tmp[ia, iG, :])
        end
        spl_f_1_T[iHH, iG] = Spline1D(a_grid, f_1_T[:, iG])
        spl_wa = Spline1D(a_grid, ω[:, iHH, iG])

        f2[iHH, iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa) * aa^(α * β / σ / (1 - η)) * spl_s(aa)^(ϕ * β / σ / (1 - η)) * spl_dw(aa)^(η * β / σ / (1 - η)), lowbnd, upbnd)[1]
        f4[iHH, iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa) * spl_wa(aa), lowbnd, upbnd)[1]

        mass_T[iHH, iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa), lowbnd, upbnd)[1]

        # Compute tomorrow's aggregate human capital in teaching ('HH_T'):
        # Note: this is the implied HH_T from the law of motion, NOT the fixed HH_fp.
        HH_T_0[iG, iHH] = ((1 - t[iHH, 1])^η * η^η * (2 * HH_fp / M)^σ)^(β / σ / (1 - η)) * f2[iHH, iG]
        # Total earnings of teachers in each group:
        E_T[iHH, iG] = f4[iHH, iG]
        # Total output of teachers in each group:
        Y_T[iHH, iG] = sum(H_O_0[:, iG, iHH] .* a_by_occ)
    end

    H_grid[iHH] = sum(HH_T_0[:, iHH] .* gm)
    for iO in 1:n_O-1
        H_O[iO, iHH] = sum(H_O_0[iO, :, iHH] .* gm)
    end
    for iG in 1:n_G
        sum_E_O[iHH, iG] = sum(E_O[:, iHH, iG])
        sum_Y_O[iHH, iG] = sum(Y_O[:, iHH, iG])
    end
    t[iHH, 2] = sum(E_T[iHH, :] .* gm) / (sum(E_T[iHH, :] .* gm) + sum(sum_E_O[iHH, :] .* gm))
    convT = abs(t[iHH, 2] - t[iHH, 1])
    println("convT=", round(convT, digits=3))
    println("t[iHH,:]=", round.(t[iHH, :], digits=3))
    t[iHH, 1] = (1 - ν) * t[iHH, 1] + ν * t[iHH, 2]

    iterT = iterT + 1
end

# No TFP growth loop — using benchmark-calibrated productivities directly.
# No HH fixed-point loop — HH_fp is exogenously held at 1970 value.

# Report the implied next-period HH_T from the law of motion (will differ from HH_fp):
println("Implied HH_T from law of motion = ", H_grid[iHH])
println("Fixed HH_fp (1970 value) = ", HH_fp)
println("Difference = ", H_grid[iHH] - HH_fp)


#############################################################
# Compute wage dispersion at steady state (iHH)             #
#############################################################
for iG in 1:n_G
    # --- Teachers ---
    # Normalised p.d.f. of abilities among teachers:
    f_T[:, iHH, iG] = f_1_T[:, iG] .* pdf.(dist, a_grid) ./ mass_T[iHH, iG]
    spl_f_T[iG] = Spline1D(a_grid, f_T[:, iHH, iG])

    # 10th percentile of teacher ability distribution:
    fn_F_T_10p(x) = 0.1 - quadgk(aa -> spl_f_T[iG](aa), lowbnd, x)[1]
    a_T_10p[iHH, iG] = find_zero(fn_F_T_10p, (lowbnd, upbnd), Bisection())

    # 90th percentile of teacher ability distribution:
    fn_F_T_90p(x) = 0.9 - quadgk(aa -> spl_f_T[iG](aa), lowbnd, x)[1]
    a_T_90p[iHH, iG] = find_zero(fn_F_T_90p, (lowbnd, upbnd), Bisection())

    # 90-10 wage ratio for teachers:
    spl_ω[iHH] = Spline1D(a_grid, ω[:, iHH, iG])
    ω_90_10[iHH, iG] = spl_ω[iHH](a_T_90p[iHH, iG]) / spl_ω[iHH](a_T_10p[iHH, iG])

    # --- Non-teaching occupations ---
    for iO in 1:n_O-1
        spl_inv = Spline1D(a_O_thresh[:, iHH, iG, iO], a_grid)
        f_1_O[:, iHH, iG, iO] = marginal[:, iO, iG] .* cdf.(dist, spl_inv(a_grid))
        spl_f_1_O[iHH, iG, iO] = Spline1D(a_grid, f_1_O[:, iHH, iG, iO])
        mass_O[iHH, iG, iO] = quadgk(aa -> pdf(dist, aa) * spl_f_1_O[iHH, iG, iO](aa), lowbnd, upbnd)[1]

        f_O[:, iHH, iG, iO] = f_1_O[:, iHH, iG, iO] .* pdf.(dist, a_grid) ./ mass_O[iHH, iG, iO]
        spl_f_O[iHH, iG, iO] = Spline1D(a_grid, f_O[:, iHH, iG, iO])

        fn_F_O_10p(x) = 0.1 - quadgk(aa -> pdf(dist, aa) * spl_f_1_O[iHH, iG, iO](aa) / mass_O[iHH, iG, iO], lowbnd, x)[1]
        a_O_10p[iHH, iG, iO] = find_zero(fn_F_O_10p, (lowbnd, upbnd), Bisection())

        fn_F_O_90p(x) = 0.9 - quadgk(aa -> pdf(dist, aa) * spl_f_1_O[iHH, iG, iO](aa) / mass_O[iHH, iG, iO], lowbnd, x)[1]
        a_O_90p[iHH, iG, iO] = find_zero(fn_F_O_90p, (lowbnd, upbnd), Bisection())

        spl_h = Spline1D(a_grid, h_O[:, iHH, iG, iO])
        w_90_10[iHH, iG, iO] = spl_h(a_O_90p[iHH, iG, iO]) / spl_h(a_O_10p[iHH, iG, iO])
    end

    # Weighted average 90-10 ratio across non-teaching occupations for group iG:
    w_90_10[iHH, iG, n_O] = w_90_10[iHH, iG, 1:n_O-1]' * share_occ_data[:, iG]
end

# Pooled 90-10 teacher wage ratio (across groups):
function fn_f_T_all_g(x)
    sum(spl_f_T[iG](x) * mass_T[iHH, iG] / sum(mass_T[iHH, :]) for iG in 1:n_G)
end
fn_F_T_10p_all_g(x) = 0.1 - quadgk(aa -> fn_f_T_all_g(aa), lowbnd, x)[1]
fn_F_T_90p_all_g(x) = 0.9 - quadgk(aa -> fn_f_T_all_g(aa), lowbnd, x)[1]
a_T_10p[iHH, n_G+1] = find_zero(fn_F_T_10p_all_g, (lowbnd, upbnd), Bisection())
a_T_90p[iHH, n_G+1] = find_zero(fn_F_T_90p_all_g, (lowbnd, upbnd), Bisection())
ω_90_10[iHH, n_G+1] = spl_ω[iHH](a_T_90p[iHH, n_G+1]) / spl_ω[iHH](a_T_10p[iHH, n_G+1])

# Pooled weighted-average 90-10 ratio across non-teaching occupations:
w_90_10[iHH, n_G+1, n_O] = mean(w_90_10[iHH, 1:n_G, n_O])

# Compute aggregate productivity at counterfactual steady state
aggA[1] = sum(sum_Y_O[iHH, :] .* gm) / sum(gm[1] .* mass_O[iHH, 1, :] + gm[2] .* mass_O[iHH, 2, :])
aggA[2] = sum((sum_Y_O[iHH, :] .- Y_O[n_O-1, iHH, :]) .* gm) / sum(gm[1] .* mass_O[iHH, 1, 1:n_O-2] + gm[2] .* mass_O[iHH, 2, 1:n_O-2])
println("aggA (counterfactual 6 - fixed HH_fp) = ", round.(aggA, digits=6))


# Transition dynamics are skipped

###########################################################
##################  PARAMETERIZATION  #####################
###########################################################
global av_e_O, av_e_T, av_s_T, av_a_rank_T
av_e_O = zeros(n_G)
av_e_T = zeros(n_G)
av_s_T = zeros(n_G)
av_a_rank_T = zeros(n_G)
for iG in 1:n_G
    spl_e_T = Spline1D(a_grid, e_T[:, iHH, iG])
    spl_s_T = Spline1D(a_grid, s_T[:, iHH, iG])
    av_e_T[iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa) * spl_e_T(aa), lowbnd, upbnd)[1] / mass_T[iHH, iG]
    av_s_T[iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa) * spl_s_T(aa), lowbnd, upbnd)[1] / mass_T[iHH, iG]
    av_a_rank_T[iG] = quadgk(aa -> pdf(dist, aa) * spl_f_1_T[iHH, iG](aa) * cdf(dist, aa), lowbnd, upbnd)[1] / mass_T[iHH, iG]
end

h_T_avg_ss = zeros(n_G)
H_T_total = zeros(n_G)
for iG in 1:n_G
    spl_h_T_ss = Spline1D(a_grid, h_T[:, iHH, iG])
    h_T_avg_ss[iG] = quadgk(aa -> spl_h_T_ss(aa) * pdf(dist, aa) * spl_f_1_T[iHH, iG](aa), lowbnd, upbnd)[1] / mass_T[iHH, iG]
    H_T_total[iG] = h_T_avg_ss[iG] * mass_T[iHH, iG] * gm[iG]
end

# Save counterfactual parameterization in JLD file:
cd(string("./parameterization/", paramname))
fnameJLD_cf = string("counter_6_", fyear, ".jld")
save(fnameJLD_cf, "a_by_occ", a_by_occ, "τ_w", τ_w, "τ_w_opt", τ_w_opt, "τ_e", τ_e, "a_T_thresh", a_T_thresh, "t", t, "H_grid", H_grid, "H_O", H_O, "HH_fp", HH_fp, "HH_fp_1970", HH_fp_1970, "α", α, "β", β, "η", η, "σ", σ, "μ", μ, "ϕ", ϕ, "γ", γ, "κ", κ, "theta", theta, "λf", λf, "λm", λm, "iHH", iHH, "a_grid", a_grid, "a_O_10p", a_O_10p, "a_O_90p", a_O_90p, "a_T_10p", a_T_10p, "a_T_90p", a_T_90p, "h_T", h_T, "f_T", f_T, "f_O", f_O, "mass_T", mass_T[iHH, :], "mass_O", mass_O[iHH, :, :], "aggA", aggA, "share_occ", share_occ, "ω_90_10", ω_90_10, "w_90_10", w_90_10)
cd("..")
cd("..")

println("____________")
println("COUNTERFACTUAL 6 RESULTS (year = ", year, ", HH_fp fixed at 1970 value)")
println("share of teachers among women= ", mass_T[iHH, 1])
println("share of teachers among men= ", mass_T[iHH, 2])
println(" ")
println("tax rate= ", t[iHH, 1])
println("output= ", sum(Y_T[iHH, :] .* gm) + sum(sum_Y_O[iHH, :] .* gm))
println("HH_fp (fixed at 1970) = ", HH_fp)
println("Implied HH_T from law of motion = ", H_grid[iHH])
println("aggA (counterfactual 6) = ", round.(aggA, digits=6))

# Write counterfactual moments to CSV file using DataFrame:
df = DataFrame(year=year, HH_fp_1970=round(HH_fp_1970; digits=6), λf=round(λf; digits=4), share_teachers_female=round(mass_T[iHH, 1]; digits=4), κ=round(κ; digits=4), share_teachers_male=round(mass_T[iHH, 2]; digits=4), γ=round(γ; digits=4), p90_p10_ω_teachers=round(ω_90_10[iHH, end]; digits=2), θ=round(theta; digits=4), p90_p10_w_other=round(w_90_10[iHH, n_G+1, n_O]; digits=2), η=round(η; digits=4), α=round(α; digits=2),
    h_T_avg_female=round(h_T_avg_ss[1]; digits=6),
    h_T_avg_male=round(h_T_avg_ss[2]; digits=6),
    H_T_total_female=round(H_T_total[1]; digits=6),
    H_T_total_male=round(H_T_total[2]; digits=6),
    H_T_total_all=round(sum(H_T_total); digits=6),
    implied_HH_T=round(H_grid[iHH]; digits=6))
# If it exists, load previous counterfactual results, append 'df':
fnameCSV_cf = string("./results/", paramname, "/counter_6_moments.csv")
if isfile(fnameCSV_cf) == true
    moments_cf = DataFrame(CSV.File(fnameCSV_cf))
    append!(moments_cf, df, cols=:union)
else
    moments_cf = df
end
show(moments_cf)

# Write DataFrame to CSV file:
CSV.write(fnameCSV_cf, moments_cf)

println("")
println("Counterfactual 6 for year ", year, " complete.")
println("================================================================")
println("")

end # end loop over counterfactual years