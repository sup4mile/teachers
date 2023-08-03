using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim, Roots, PyCall, Random
# Change directory, if necessary:
# cd("./GitHub/teachers/julia/codes_w_exo_wage")

# Set (some of the) parameters:
# Population size:
M=1.0
# Groups (no discrimination / no barrier group is last element)
g = ["female","male"]
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
occ = [occ_O;occ_T]

share_occ_data = Array{Float64,2}(undef,length(occ)-1,2)
w_90_10_data = Array{Float64,1}(undef,2) 

# Select calendar your for calibration (1970, 1990, or 2010)
year = 1970
# Load data for selected year:
if year == 1970
    share_occ_data[:,1] = tab1["K30:K49"] # Census 1970 for Project TALENT (women)
    share_occ_data[:,2] = tab1["C30:C49"] # Census 1970 for Project TALENT (men)
    w_90_10_data[1] = tab2["C50"] # Dispersion in non-teaching occupations
    w_90_10_data[2] = tab2["C51"] # Dispersion in teaching
elseif year==1990
    share_occ_data[:,1] = tab1["M30:M49"] # Census 1990 for NLSY79 (women)
    share_occ_data[:,2] = tab1["E30:E49"] # Census 1990 for NLSY79 (men)
    w_90_10_data[1] = tab2["E50"] # Dispersion in non-teaching occupations
    w_90_10_data[2] = tab2["E51"] # Dispersion in teaching
elseif year==2010
    share_occ_data[:,1] = tab1["O30:O49"] # ACS 2009-2013 for NLSY97 (women)
    share_occ_data[:,2] = tab1["G30:G49"] # ACS 2009-2013 1990 for NLSY97 (men)
    w_90_10_data[1] = tab2["G50"] # Dispersion in non-teaching occupations
    w_90_10_data[1] = tab2["G51"] # Dispersion in teaching
else
    println("No a_by_occ for the selected year")
end

# Innovation step size for updates:
ν = 0.8 # for tax rate that balances government's budget
ν2 = 0.8 # for fixed point of aggregate human capital in teaching

n_G=length(g) # number of "groups"
n_O=length(occ) # number of occupations

gm=zeros(1,n_G) # measure of in individuals in groups 1,...,n_G
for iG in 1:n_G-1
    gm[iG] = M/(2*n_G)
end
gm[end] = M/2 - sum(gm)
# Load parameters and results from previous parameterization:
    cd("./parameterization")
    fyear = string(year)
    fnameJLD = string("previousParameterization",fyear,".jld")
    d = load(fnameJLD)
    a_by_occ_initial = d["a_by_occ"]
    # a_by_occ = d["a_by_occ"]
    τ_w_initial = d["τ_w_opt"]
    τ_e_initial = d["τ_e"]
    a_T_thresh = d["a_T_thresh"]
    t = (d["t"])
    H_grid = d["H_grid"]
    H_O = d["H_O"]
    HH_fp = d["HH_fp"]
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

# Update model parameters, if required:
λf = .895 # composite barrier for women in non-teaching occupations (note: share of female teachers is decreasing in λf)
κ = 2.015 # scale parameter of teachers' wage profile

# Distribution of abilities:
dist = Frechet(theta,1)

 # Vector of labor market discrimination in each occupation (relative to 'T')
τ_w=zeros(n_O-1,n_G)
# Vector of education barriers in each occupation (relative to 'T')
τ_e=zeros(n_O-1,n_G)

# to match shares for Other for women
τ_w_opt=zeros(n_O-1,n_G) # n_G-element vector of labor market discrimination in 'O' (relative to 'T')

τ_w[:,2]=fill!(τ_w[:,2],0)
τ_w[:,2]=ones(n_O-1).-λm*(ones(n_O-1).-τ_w[:,2])
τ_e=fill!(τ_e,0.0)

quantile_top = 1 - 1e-4
quantile_bottom = 1e-3

# Reset grid length for individual abilities (if necessary):
n_a = 120
# Set up the grid for aggregate human capital in teaching:
n_H = 11
# Set the percentage range above / below fixed point:
range_H = 0.25
# Update the grid for abilities (if necessary):
    i_grid = log10.(range(1e4^(quantile_bottom),1e4^(quantile_top),n_a))./4
    b_grid = quantile.(dist,i_grid)
    spl_a_grid = Spline1D(range(1,length(a_grid),length(a_grid)),a_grid)
    a_grid_upd = spl_a_grid.(range(1,length(a_grid),n_a))
    c_grid = 1 .* b_grid + 0 .* a_grid_upd
    # Update size of loaded arrays, if necessary:
    if n_a != length(a_grid)
        h_T_initial_upd = Array{Float64,3}(undef,n_a,n_H,n_G)
        a_T_thresh_upd = Array{Float64,4}(undef,n_a,n_H,n_G,n_O)
        for iH = 1:n_H
            for iG = 1:n_G
                spl_h_T = Spline1D(a_grid,h_T_initial[:,iH,iG])
                h_T_initial_upd[:,iH,iG] = spl_h_T.(c_grid)
                for iO = 1:n_O-1
                    spl_a_T_thresh = Spline1D(a_grid,a_T_thresh[:,iH,iG,iO])
                    a_T_thresh_upd[:,iH,iG,iO] = spl_a_T_thresh.(c_grid)
                end
            end
        end
        h_T_initial = h_T_initial_upd
        a_T_thresh = a_T_thresh_upd
    end
    a_grid = c_grid

# Integral bounds:
lowbnd = a_grid[1] # set to smallest element in 'a_grid'
upbnd = a_grid[end] # set to largest element in 'a_grid'

# Given a productivity vector 'a_by_occ' with productivities in non-teaching occupations, the function 'share()' calculates share of people choosing occupation 'occ_id':
marg=Array{Float64,3}(undef,n_a,n_O-1,n_G)
function share(occ_id,iG,a_by_occ,τ_w,τ_e)
    # if iG == 1
    # println("iO = ",occ_id,", τ_w = ",τ_w)
    # end
    for ia in 1:n_a
        B_tmp = (((ones(n_O-2).-τ_w[occ_id])./(ones(n_O-2).-τ_w[1:end.!=occ_id])).*(((ones(n_O-2).+τ_e[occ_id])./(ones(n_O-2).+τ_e[1:end.!=occ_id])).^(-η)).*(a_by_occ[occ_id]./a_by_occ[1:end .!=occ_id]))
        # if minimum(B_tmp) <= 0
        #     println("Negative base!")
        # end
        A_tmp = B_tmp.^(1/α)
        # A_tmp = (( ((ones(n_O-2).-τ_w[occ_id])./(ones(n_O-2).-τ_w[1:end.!=occ_id])).*(((ones(n_O-2).+τ_e[occ_id])./(ones(n_O-2).+τ_e[1:end.!=occ_id])).^(-η)).*(a_by_occ[occ_id]./a_by_occ[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG] = prod(cdf.(dist,A_tmp.*a_grid[ia]))
    end
    # if iG == 1
    # println("iO = ",occ_id,", τ_w = ",τ_w[occ_id])
    # end
    spl_marg = Spline1D(a_grid, marg[:,occ_id,iG])
    # Function returns two output arguments:
    out1 = quadgk(aa -> spl_marg(aa)*pdf(dist,aa),lowbnd,upbnd)[1]
    # if iG == 1
    # println("out1 = ",out1)
    # end
    out2 = marg[:,occ_id,iG]
    return out1, out2 
end

# Calibrate 'a_by_occ' to match labor market shares for men in non-teaching occupations:
function calibrate_A(x)
    # Share of individuals employed in occupation 'iO' among all non-teaching individuals (the elements of share_occ_model will add up to 1):
    share_occ_model = Array{Float64,1}(undef,n_O-1)
    for iO = 1:(n_O-1)
        share_occ_model[iO]=share(iO,2,x,τ_w[:,2],τ_e[:,2])[1]
    end
    share_occ_model=share_occ_model./sum(share_occ_model)
    # Sum of squared differences (in relative terms) between data and model:
    return sum(((share_occ_model[1:end-1].-share_occ_data[1:end-1,2])./share_occ_data[1:end-1,2]).^2)
end

# Compute occupation-specific productivies to match employment shares of men (group 2):
res=optimize(calibrate_A,a_by_occ_initial, show_trace=false, iterations=10000)
a_by_occ=Optim.minimizer(res)
println("Sum of squared distances between a_by_occ_initial and a_by_occ for men is ",sum(abs.(a_by_occ- a_by_occ_initial)))

# Calibrate 'τ_w' to match labor marker shares for women in non-teaching occupations:
function calibrate_τ(x)
    # Share of individuals employed in occupation 'iO' among all non-teaching individuals (the elements of share_occ_model will add up to 1):
    share_occ_model = Array{Float64,1}(undef,n_O-1)
    for iO = 1:(n_O-1)
        share_occ_model[iO]=share(iO,1,a_by_occ,x,τ_e[:,1])[1]
    end
    share_occ_model=share_occ_model./sum(share_occ_model)
    # Sum of squared differences (in relative terms) between data and model:
    return sum(((share_occ_model[1:end-1].-share_occ_data[1:end-1,2])./share_occ_data[1:end-1,2]).^2)
end
# Compute labor market barriers for women (group 1) with a box constraint for τ_w (which has to be smaller than 1 for all occupations):

# Optimization with box constraint (τ_w < 1):
# lower = -Inf*ones(size(τ_w_initial[:,1]))
# upper = ones(size(τ_w_initial[:,1]))
# res = optimize(calibrate_τ, lower, upper, τ_w_initial[:,1], Fminbox(NelderMead()), Optim.Options(show_trace=false,iterations=10000,outer_iterations=2))

# Unconstrained optimization (preferred, if possible; Nelder Mead algorithm doesn't always work well with box constraints):
res = optimize(calibrate_τ,τ_w_initial[:,1], show_trace=false, iterations=10000)
τ_w_opt[:,1] = Optim.minimizer(res)
println("Sum of squared distances between τ_w_initial and τ_w_opt for women is ",sum(abs.(τ_w_opt-τ_w_initial)))
τ_w[:,1]=ones(n_O-1).-λf*(ones(n_O-1).-τ_w_opt[:,1])

# Compute the share of individuals employed in non-teaching occupations in each group (indexed by 'iG'):
share_occ = Array{Float64,2}(undef,n_O-1,n_G)
marginal = Array{Float64,3}(undef,n_a,n_O-1,n_G)
for iG in 1:n_G
    for iO = 1:(n_O-1)
        share_occ[iO,iG]=share(iO,iG,a_by_occ,τ_w[:,iG],τ_e[:,iG])[1]
        # Note: the sum of 'share_occ' will add up to 1*n_G.
        # Fraction of workers with abilities 'a' in (non-teaching) occupation 'iO' who actually work in occupation 'iO':
        marginal[:,iO,iG]=share(iO,iG,a_by_occ,τ_w[:,iG],τ_e[:,iG])[2]
    end
    share_occ[:,iG]=share_occ[:,iG]./sum(share_occ[:,iG])
end

a_T_thresh = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)
a_O_thresh = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)

e_T = Array{Float64,3}(undef,n_a,n_H,n_G)
e_O = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)

s_T = Array{Float64,3}(undef,n_a,n_H,n_G)

h_T = Array{Float64,3}(undef,n_a,n_H,n_G)
h_O = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)
ω = Array{Float64,3}(undef,n_a,n_H,n_G)
w = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)
ω_90_10 = Array{Float64,2}(undef,n_H,n_G+1)
w_90_10 = Array{Float64,3}(undef,n_H,n_G+1,n_O)
der_ω = Array{Float64,3}(undef,n_a,n_H,n_G)

# a_T_10p = Array{Float64,2}(undef,n_H,n_G+1)
# a_T_90p = Array{Float64,2}(undef,n_H,n_G+1)
# a_O_10p = Array{Float64,3}(undef,n_H,n_G,n_O)
# a_O_90p = Array{Float64,3}(undef,n_H,n_G,n_O)
a_T_90_10 = Array{Float64,2}(undef,n_H,n_G)
a_O_90_10 = Array{Float64,3}(undef,n_H,n_G,n_O)

spl_T = Array{Spline1D,1}(undef,n_G)
spl_O = Array{Spline1D,1}(undef,n_G)
spl_T_inv = Array{Spline1D,1}(undef,n_G)
spl_O_inv = Array{Spline1D,1}(undef,n_G)
spl_ω = Array{Spline1D,1}(undef,n_H)
spl_w = Array{Spline1D,1}(undef,n_G)
spl_h_O = Array{Spline1D,1}(undef,n_G)
spl_pdf_T = Array{Spline1D,1}(undef,n_G)
wb_T = Array{Float64,2}(undef,n_H,n_G)
wb_O = Array{Float64,2}(undef,n_H,n_G)
E_T = Array{Float64,2}(undef,n_H,n_G)
E_O = Array{Float64,3}(undef,n_O-1,n_H,n_G)
Y_T = Array{Float64,2}(undef,n_H,n_G)
Y_O = Array{Float64,3}(undef,n_O-1,n_H,n_G)
sum_E_O = Array{Float64,2}(undef,n_H,n_G)
sum_Y_O = Array{Float64,2}(undef,n_H,n_G)

HH_T_0 = Array{Float64,1}(undef,n_G)
HH_T = Array{Float64,1}(undef,n_H)
H_O_0 = Array{Float64,3}(undef,n_O-1,n_G,n_H)
H_O = Array{Float64,2}(undef,n_O-1,n_H)

EN = Array{Float64,1}(undef,n_G)
EN2 = Array{Float64,1}(undef,n_G)

f_2 = Array{Float64,2}(undef,n_a,n_G)
f_O = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)

spl_T_thresh=Array{Spline1D,1}(undef,1)
# spl_marg = Array{Spline1D,1}(undef,1)

f1 = Array{Float64,2}(undef,n_O-1,n_G)
mass_O = Array{Float64,3}(undef,n_H,n_G,n_O-1)
f3 = Array{Float64,2}(undef,n_O-1,n_G)
f_1_T = Array{Float64,2}(undef,n_a,n_G)
f_1_O = Array{Float64,4}(undef,n_a,n_H,n_G,n_O-1)
spl_f_1_T = Array{Spline1D,2}(undef,n_H,n_G)
spl_f_1_O = Array{Spline1D,3}(undef,n_H,n_G,n_O-1)
f_T = Array{Float64,3}(undef,n_a,n_H,n_G)
spl_f_T = Array{Spline1D,1}(undef,n_G)
spl_f_O = Array{Spline1D,3}(undef,n_H,n_G,n_O)
mass_T = Array{Float64,2}(undef,n_H,n_G)
f2 = Array{Float64,2}(undef,n_H,n_G)
f4 = Array{Float64,2}(undef,n_H,n_G)
ratio = Array{Float64,1}(undef,n_G)
f_1_T_tmp=Array{Float64,3}(undef,n_a,n_G,n_O-1)

#######################################
# Exogenous wage profile for teachers #
#######################################
ω_fn(h_T) = κ.*h_T.^γ
# Derivative of wage profile:
der_ω_fn(h_T) = κ.*γ.*h_T.^(γ-1)


#######################
# Steady-state of H_T #
#######################

# Set the array index to the midpoint of H_grid:
iHH = convert(Int,ceil(n_H/2))
# Iteration and tolerance settings for fixed-point problems:
# (a) Aggregate human capital:
convHH = 1
tolHH = .0025
# (b) Income tax rate:
tolT= 2.5e-4
maxiterT = 100

# Time investment doesn't depend on any endogenous variables, only on parameters:
s_O = μ*ϕ / (μ*ϕ + 1 - η)
# Initiate 'while' loop (not indexed) over aggregate human capital in teaching:
while convHH > tolHH
    global s_O
    global HH_fp
    global convHH
    global tolHH
    global H_grid
    println("HH_fp = ", HH_fp)
    println("t[iHH,1] = ",t[iHH,1])
    println("")
    H_grid[iHH] = HH_fp
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        for iG in 1:n_G
            for ia in 1:n_a
                a=a_grid[ia]
                fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k)=(η^η*(1-t[iHH,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                hh_T = find_zero(fn_h_T,h_T_initial[ia,iHH,iG])
                h_T[ia,iHH,iG]=hh_T
                s_T[ia,iHH,iG]=fn_s_T(hh_T)
                e_T[ia,iHH,iG]=η*(1-t[iHH,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iHH,iG]=ω_fn(hh_T)
                der_ω[ia,iHH,iG]=der_ω_fn(hh_T)
                for iO in 1:n_O-1
                    # Occupational choice threshold (workers with 'a' in occupation 'iO' become teachers if 'a' < 'a_O_thresh[.,.,.]' for given teaching ability a_T = a[ia]):
                    a_O_thresh[ia,iHH,iG,iO] = ((1+τ_e[iO,iG])/(1-τ_w[iO,iG])*(1+τ_e[iO,iG])^(-(1-η))*((1-s_T[ia,iHH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iHH,iG]/s_O)^ϕ*(der_ω[ia,iHH,iG]/a_by_occ[iO])*((ω[ia,iHH,iG]/(der_ω[ia,iHH,iG]*h_T[ia,iHH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iHH,iG,iO] = (η^η*(1-t[iHH,1])^η*(1-τ_w[iO,iG])^η/(1+τ_e[iO,iG])^η*a_by_occ[iO]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                    e_O[ia,iHH,iG,iO] = η * (1-t[iHH,1]) * (1-τ_w[iO,iG]) / (1+τ_e[iO,iG]) * a_by_occ[iO] * h_O[ia,iHH,iG,iO]
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iHH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iHH,iG]))

            for iO in 1:n_O-1
                spl_marg = Spline1D(a_grid, marginal[:,iO,iG])
                # Inverse occupational threshold 
                spl_inv = Spline1D(a_O_thresh[:,iHH,iG,iO],a_grid)
                f3[iO,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[iO,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[iO,iG,iHH]=(((1-t[iHH,1])*(1-τ_w[iO,iG])/(1+τ_e[iO,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*a_by_occ[iO]^η)^(1/(1-η))*f3[iO,iG]
                # Total earnings of others in each group:
                E_O[iO,iHH,iG]=a_by_occ[iO]*H_O_0[iO,iG,iHH]
                # Total output of others in each group:
                Y_O[iO,iHH,iG]=a_by_occ[iO]*H_O_0[iO,iG,iHH]

                for ia in 1:n_a
                    a = a_grid[ia]
                    # Fraction of teachers occupation-by-occupation, given 'a' in teaching (see integration in equation 17 of manuscript for details):
                    f_1_T_tmp[ia,iG,iO]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iHH,iG,iO])[1],0.0])
                end
            end

            for ia in 1:n_a
                # Fraction of teachers given a_T (across all occupations):
                f_1_T[ia,iG]= sum(f_1_T_tmp[ia,iG,:])
            end
            spl_f_1_T[iHH,iG]=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iHH,iG])

            f2[iHH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iHH,iG](aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iHH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iHH,iG](aa)*spl_wa(aa),lowbnd,upbnd)[1]

            mass_T[iHH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iHH,iG](aa),lowbnd,upbnd)[1]

            # Compute tomorrow's aggregate human capital in teaching ('HH_T'):
            HH_T_0[iG]= ( (1-t[iHH,1])^η*η^η*(2*HH_fp/M)^σ)^(β/σ/(1-η))*f2[iHH,iG]
            #*s_T^ϕ*(β/σ)^η*(sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_O-1).-τ_w[:,iG])./(ones(n_O-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iHH,iG]=f4[iHH,iG]
            # Total output of teachers in each group:
            Y_T[iHH,iG] = sum(H_O_0[:,iG,iHH].*a_by_occ)
            #η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_O-1).-τ_w[:,iG])./(ones(n_O-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])
        end
        
        H_grid[iHH]=sum(HH_T_0.*gm)
        for iO in 1:n_O-1
            H_O[iO,iHH]=sum(H_O_0[iO,:,iHH].*gm)
        end
        for iG in 1:n_G
            sum_E_O[iHH,iG]=sum(E_O[:,iHH,iG])
            sum_Y_O[iHH,iG]=sum(Y_O[:,iHH,iG])
        end
        t[iHH,2]=sum(E_T[iHH,:].*gm)/(sum(E_T[iHH,:].*gm)+sum(sum_E_O[iHH,:].*gm))
        #Y=sum(Y_T[iHH,:].*gm)+sum(sum_Y_O[iHH,:].*gm)
        convT = abs(t[iHH,2]-t[iHH,1])
         println("convT=",round(convT,digits=3))
         println("t[iHH,:]=",round.(t[iHH,:],digits=3))
        t[iHH,1] =  (1-ν)*t[iHH,1] + ν*t[iHH,2]

        iterT = iterT+1
    end
    convHH = abs(log(H_grid[iHH]/HH_fp))
    println(convHH)
    println(H_grid[iHH])
    HH_fp = (1-ν2)*HH_fp + ν2*H_grid[iHH]
end
println("Found fixed point for human capital in teaching!")

##############################
# Transition dynamics of H_T #
##############################

# Construct grid for H_T (aggregate human capital in in teaching) around fixed point:
H_grid = HH_fp .* exp.(range(-range_H,range_H,n_H))
for iH in 1:n_H
    println("Transition dynamics for H: ", iH," of ", n_H," (H = ", round(H_grid[iH];digits=6),")")
    println("")
    # Today's aggregate human capital in teaching:
    H = H_grid[iH]
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        for iG in 1:n_G
            for ia in 1:n_a
                a=a_grid[ia]
                # Human capital investments and wages in teaching:
                fn_s_T(k) = μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k) = (η^η*(1-t[iH,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*H/M)^σ)^(1/(1-η))-k
                hh_T = find_zero(fn_h_T,h_T_initial[ia,iH,iG])
                h_T[ia,iH,iG] = hh_T
                s_T[ia,iH,iG] = fn_s_T(hh_T)
                e_T[ia,iH,iG] = η*(1-t[iH,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iH,iG] = ω_fn(hh_T)
                der_ω[ia,iH,iG] = der_ω_fn(hh_T)
                # Human capital investments and wages in other occupations:
                for iO in 1:n_O-1
                    # Occupational threshold ():
                    a_O_thresh[ia,iH,iG,iO] = ((1+τ_e[iO,iG])/(1-τ_w[iO,iG])*(1+τ_e[iO,iG])^(-(1-η))*((1-s_T[ia,iH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iH,iG]/s_O)^ϕ*(der_ω[ia,iH,iG]/a_by_occ[iO])*((ω[ia,iH,iG]/(der_ω[ia,iH,iG]*h_T[ia,iH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iH,iG,iO] = (η^η*(1-t[iH,1])^η*(1-τ_w[iO,iG])^η/(1+τ_e[iO,iG])^η*a_by_occ[iO]^η*a_grid[ia]^α*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))
                    e_O[ia,iH,iG,iO] = η * (1-t[iH,1]) * (1-τ_w[iO,iG]) / (1+τ_e[iO,iG]) * a_by_occ[iO] * h_O[ia,iHH,iG,iO]
                    w[ia,iH,iG,iO] = a_by_occ[iO] * (1-τ_w[iO,iG]) * h_O[ia,iH,iG,iO]
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iH,iG]))

            for iO in 1:n_O-1
                # Spline for fraction of workers in occupation 'iO':
                spl_marg = Spline1D(a_grid, marginal[:,iO,iG])
                # Inverse occupational threshold between occupation 'iO' and teaching: 
                spl_inv = Spline1D(a_O_thresh[:,iH,iG,iO],a_grid)
                f3[iO,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[iO,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[iO,iG,iH]=(((1-t[iH,1])*(1-τ_w[iO,iG])/(1+τ_e[iO,iG]))^η*η^η*(2*H/M)^σ*s_O^ϕ*a_by_occ[iO]^η)^(1/(1-η))*f3[iO,iG]
                # Total earnings of others in each group:
                E_O[iO,iH,iG]=a_by_occ[iO]*H_O_0[iO,iG,iH]
                # Total output of others in each group:
                Y_O[iO,iH,iG]=a_by_occ[iO]*H_O_0[iO,iG,iH]

                for ia in 1:n_a
                    a = a_grid[ia]
                    # fract of teach occ by occ, given a_T
                    f_1_T_tmp[ia,iG,iO]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iH,iG,iO])[1],0.0])
                    #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iH,iG,iO])[1]
                end
            end
            for ia in 1:n_a
                # Fraction of teachers given a_T (across all occupations):
                f_1_T[ia,iG]= sum(f_1_T_tmp[ia,iG,:])
            end
            
            # Compute 10th and 90th percentile of distribution of wages:
            
            # (a) Teachers:
            # (a.1) Compute p.d.f.:
            spl_f_1_T[iH,iG]=Spline1D(a_grid, f_1_T[:,iG])
            mass_T[iH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa),lowbnd,upbnd)[1]
            f_T[:,iH,iG] = f_1_T[:,iG].*pdf.(dist,a_grid)./mass_T[iH,iG]
            # (a.2) Compute ability at 10th and 90th percentile for each group using c.d.f.:
            spl_f_T[iG] = Spline1D(a_grid,f_T[:,iH,iG])
            function fn_F_T_10p(x)
                .1 - quadgk(aa -> spl_f_T[iG](aa),lowbnd,x)[1]
            end
            # Find ability of teacher at 10th percentile of the distribution:
            a_T_10p[iH,iG] = find_zero(fn_F_T_10p,a_T_10p[iH,iG])
            # a_T_10p[iH,iG] = find_zero(fn_F_T_10p,a_grid[round(Int,n_a/2)])
            function fn_F_T_90p(x)
                .9 - quadgk(aa -> spl_f_T[iG](aa),lowbnd,x)[1]
            end
            # Find ability of teacher at 90th percentile of the distribution:
            a_T_90p[iH,iG] = find_zero(fn_F_T_90p,a_T_90p[iH,iG])
            # a_T_90p[iH,iG] = find_zero(fn_F_T_90p,a_grid[round(Int,n_a/2)])
            # (a.3) Compute wages at 10th and 90th percentile of distribution for each group:
            spl_ω[iH] = Spline1D(a_grid,ω[:,iH,iG])
            ω_90_10[iH,iG] = spl_ω[iH](a_T_90p[iH,iG])/spl_ω[iH](a_T_10p[iH,iG])
            a_T_90_10[iH,iG] = a_T_90p[iH,iG]/a_T_10p[iH,iG]
            
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iH,iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*spl_wa(aa),lowbnd,upbnd)[1]
            
            # Compute tomorrow's aggregate human capital in teaching ('HH_T'):
            HH_T_0[iG]= ( (1-t[iH,1])^η*η^η*(2*H/M)^σ)^(β/σ/(1-η))*f2[iH,iG]
            # Total earnings of teachers in each group:
            E_T[iH,iG]=f4[iH,iG]
            # Total output of teachers in each group:
            Y_T[iH,iG] = sum(H_O_0[:,iG,iH].*a_by_occ)

            # (b) All other occupations:
            for iO = 1:n_O-1
                # Inverse occupational threshold between occupation 'iO' and teaching: 
                spl_inv = Spline1D(a_O_thresh[:,iH,iG,iO],a_grid)
                f_1_O[:,iH,iG,iO] = marginal[:,iO,iG] .* cdf.(dist,spl_inv(a_grid))
                # (b.1) Compute p.d.f.:
                spl_f_1_O[iH,iG,iO]=Spline1D(a_grid, f_1_O[:,iH,iG,iO])
                mass_O[iH,iG,iO]=quadgk(aa -> pdf(dist,aa)*spl_f_1_O[iH,iG,iO](aa),lowbnd,upbnd)[1]
                f_O[:,iH,iG,iO] = f_1_O[:,iH,iG,iO].*pdf.(dist,a_grid)./mass_O[iH,iG,iO]
                # (b.2) Compute ability at 10th and 90th percentile for each group and occupation using c.d.f.:
                spl_f_O[iH,iG,iO] = Spline1D(a_grid,f_O[:,iH,iG,iO])
                function fn_F_O_10p(x)
                    # .1 - quadgk(aa -> spl_f_O[iH,iG,iO](aa),lowbnd,x)[1]
                    .1 - quadgk(aa -> pdf(dist,aa)*spl_f_1_O[iH,iG,iO](aa)/mass_O[iH,iG,iO],lowbnd,x)[1]
                end
                # println("iH = ",iH,", iG = ",iG,", iO = ",iO)
                # Find ability of worker at 10th percentile of the distribution in each group and occupation:
                a_O_10p[iH,iG,iO] = find_zero(fn_F_O_10p,a_O_10p[iH,iG,iO])
                function fn_F_O_90p(x)
                    # .9 - quadgk(aa -> spl_f_O[iH,iG,iO](aa),lowbnd,x)[1]
                    .9 - quadgk(aa -> pdf(dist,aa)*spl_f_1_O[iH,iG,iO](aa)/mass_O[iH,iG,iO],lowbnd,x)[1]
                end
                # Find ability of worker at 90th percentile of the distribution:
                a_O_90p[iH,iG,iO] = find_zero(fn_F_O_90p,a_O_90p[iH,iG,iO])
                # (a.3) Compute ratio of wages at 10th and 90th percentile of distribution for each group and occupation:
                spl_h = Spline1D(a_grid,h_O[:,iH,iG,iO])
                w_90_10[iH,iG,iO] = spl_h(a_O_90p[iH,iG,iO])/spl_h(a_O_10p[iH,iG,iO])
                a_O_90_10[iH,iG,iO] = a_O_90p[iH,iG,iO]/a_O_10p[iH,iG,iO]
            end
            # Compute weighted average of 90-10 wage ratio within group 'iG' (and for given 'iH'):
            w_90_10[iH,iG,n_O] = w_90_10[iH,iG,1:n_O-1]' * share_occ_data[:,iG]
        end
        w_90_10[iH,n_G+1,n_O] = mean(w_90_10[iH,1:n_G,n_O])
        # Compute wages at 10th and 90th percentiles in teaching (pooled across groups):
        function fn_f_T_all_g(x)
            fn = 0
            for iG = 1:n_G
                fn = fn + spl_f_T[iG](x) * mass_T[iH,iG]/sum(mass_T[iH,:])
            end
            return fn
        end
        function fn_F_T_10p_all_g(x)
            .1 - quadgk(aa -> fn_f_T_all_g(aa),lowbnd,x)[1]
        end
        function fn_F_T_90p_all_g(x)
            .9 - quadgk(aa -> fn_f_T_all_g(aa),lowbnd,x)[1]
        end
        a_T_10p[iH,n_G+1] = find_zero(fn_F_T_10p_all_g,a_T_10p[iH,n_G+1])
        a_T_90p[iH,n_G+1] = find_zero(fn_F_T_90p_all_g,a_T_90p[iH,n_G+1])
        ω_90_10[iH,n_G+1] = spl_ω[iH](a_T_90p[iH,n_G+1])/spl_ω[iH](a_T_10p[iH,n_G+1])
        # Law of motion for aggregate human in teaching: 
        HH_T[iH]=sum(HH_T_0.*gm)
        for iO in 1:n_O-1
            H_O[iO,iH]=sum(H_O_0[iO,:,iH].*gm)
        end
        for iG in 1:n_G
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
            sum_Y_O[iH,iG]=sum(Y_O[:,iH,iG])
        end
        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        convT = abs(t[iH,2]-t[iH,1])
         println("convT=",round(convT,digits=3))
         println("t[iH,:]=",round.(t[iH,:],digits=3))
        t[iH,1] =  (1-ν)*t[iH,1] + ν*t[iH,2]# min(t[1,2],1-1e-3)

        iterT = iterT+1
    end
end

# What is 'N' used for?
N=zeros(n_a,n_H,n_G)
for iH in 1:n_H
    for iG in 1:n_G
        for ia in 1:n_a
            N[ia,iH,iG]=(M/2)/H_grid[iH]*(h_T[ia,iH,iG])^(β/σ)
        end
    end
end

###########################################################
##################  PARAMETERIZATION  #####################
###########################################################
av_e_O=zeros(n_G)
av_e_T=zeros(n_G)
av_s_T=zeros(n_G)
av_a_rank_T=zeros(n_G)
av_N=zeros(n_G)
for iH in 1:n_H
    for iG in 1:n_G
        spl_e_T = Spline1D(a_grid,e_T[:,iH,iG])
        spl_s_T = Spline1D(a_grid,s_T[:,iH,iG])
        spl_N=Spline1D(a_grid, N[:,iH,iG])
        av_e_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*spl_e_T(aa),lowbnd,upbnd)[1]/mass_T[iH,iG]
        av_s_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*spl_s_T(aa),lowbnd,upbnd)[1]/mass_T[iH,iG]
        av_a_rank_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*cdf(dist,aa),lowbnd,upbnd)[1]/mass_T[iH,iG]
        av_N[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iH,iG](aa)*spl_N(aa),lowbnd,upbnd)[1]/mass_T[iH,iG]
    end
end
# Save parameterization in JLD file:
cd("./parameterization")
save(fnameJLD,"a_by_occ",a_by_occ,"τ_w_opt",τ_w_opt,"τ_e",τ_e,"a_T_thresh",a_T_thresh,"t",t,"H_grid",H_grid,"H_O",H_O,"HH_fp",HH_fp,"HH_T",HH_T,"α",α,"β",β,"η",η, "σ",σ,"μ",μ,"ϕ",ϕ,"γ",γ,"κ",κ,"theta",theta,"λf",λf,"λm",λm,"iHH",iHH,"a_grid",a_grid,"a_O_10p",a_O_10p,"a_O_90p",a_O_90p,"a_T_10p",a_T_10p,"a_T_90p",a_T_90p,"h_T",h_T,"f_T",f_T,"f_O",f_O)
cd("..")

println("____________")
println("share of teachers among women= ",mass_T[iHH,1])
println("share of teachers among men= ",mass_T[iHH,2])
println("average class size= ",sum(av_N.*gm))
# println(" ")
# println("average a_rank_T_female= ",av_a_rank_T[1])
# println("average a_rank_T_male= ",av_a_rank_T[2])
# println("average e_T_female= ",av_e_T[1])
# println("average e_T_male= ",av_e_T[2])
# println("average s_T_female= ",av_s_T[1])
# println("average s_T_male= ",av_s_T[2])
# println("average w_T_female= ",E_T[1,1]/mass_T[1])
# println("average w_T_male= ",E_T[1,2]/mass_T[2])
# println("dispersion w_T_female= ")
# println("dispersion w_T_male= ")
# println(" ")
# println("average a_rank_O_female= ",)
# println("average a_rank_O_male= ",)
# println("average e_O_female= ",)
# println("average e_O_male= ",)
# println("average s_O_female= ",s_O)
# println("average s_O_male= ",s_O)
# println("average w_O_female= ",)
# println("average w_O_male= ",)
# println("dispersion w_O_female= ",)
# println("dispersion w_O_male= ",)
println(" ")
println("tax rate= ",t[1,1])
println("output= ",sum(Y_T[iHH,:].*gm)+sum(sum_Y_O[iHH,:].*gm))
println("HH_fp= ",HH_fp)

# Write parameter values and moments to CSV file using DataFrame:

df = DataFrame(λf = round(λf;digits=4),share_teachers_female = round(mass_T[iHH,1];digits=4),κ = round(κ;digits=4),share_teachers_male = round(mass_T[iHH,2];digits=4),γ = round(γ;digits=4),p90_p10_ω_teachers = round(ω_90_10[iHH,end];digits=2), θ = round(theta;digits=4),p90_p10_w_other = round(w_90_10[iHH,n_G+1,n_O];digits=2), η = round(η;digits=4), α = round(α;digits=2) )
# If it exists, load previous parameterization (in CSV format), convert to DataFrame, and append 'df':
fnameCSV = string("./results/moments",fyear,".csv")
if isfile(fnameCSV) == true
    moments = DataFrame(CSV.File(fnameCSV))
    append!(moments,df)
else
    # If the CSV file doesn't exist, create a new DataFrame named "moments":
    moments = df
end
show(moments)

# Write DataFrame to CSV file:
CSV.write(fnameCSV,moments);
