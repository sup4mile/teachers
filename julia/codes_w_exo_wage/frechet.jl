using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim, Roots, PyCall, Random
# Change directory, if necessary:
# cd("./")
# Switch to use previous result to initiate numerical optimization:
previous = 1;
# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
# cd("C:/Users/julia/Desktop/Research/Teachers")
# cd("C:\\Users\\julia\\Desktop\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
cd("..")
cd("..")
cd("./data/LaborMarketData")
# xf = XLSX.readxlsx("C:/Users/julia/Desktop/BOX/Teachers/New Yulia's Folder/LaborMarketData/Occupation_shares_v3.xlsx")
xf = XLSX.readxlsx("Occupation_shares_v3.xlsx")
sh = xf["mom"]

# Reset working directory to folder with Julia script:
cd("..")
cd("..")
cd("./julia/codes_w_exo_wage")
# labels for occuptions
occ = sh["A3:A23"]

share_occ_data=Array{Float64,2}(undef,length(occ)-1,2)

#year=1970
#year=1990
year=2010
if year==1970
    share_occ_data[:,1] = sh["H3:H22"] # Census 1970 for Project TALENT
    share_occ_data[:,2] = sh["B3:B22"] # Census 1970 for Project TALENT
elseif year==1990
    share_occ_data[:,1] = sh["I3:I22"] # Census 1990 for NLSY79
    share_occ_data[:,2] = sh["C3:C22"] # Census 1990 for NLSY79
elseif year==2010
    share_occ_data[:,1] = sh["J3:J22"] # Census 1990 for NLSY97
    share_occ_data[:,2] = sh["D3:D22"] # Census 1990 for NLSY97
else
    println("No a_by_occ for the selected year")
end

# Innovation step size for updates:
ν = 0.9 # for tax rate that balances government's budget
ν2 = 0.8 # for fixed point of aggregate human capital in teaching

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)

if previous == 1
    cd("./parameterization")
    d = load("previousParameterization.jld")
    a_by_occ_initial = d["a_by_occ"]
    # a_by_occ = d["a_by_occ"]
    τ_w_initial = d["τ_w_opt"]
    τ_e_initial = d["τ_e"]
    a_T_thresh = d["a_T_thresh"]
    t = transpose(d["t"])
    HH_T = d["HH_T"]
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
    cd("..")
else
    # productivity in 'Other' occupations:
    a_by_occ_initial=Array{Float64,1}(undef,n_occ-1)
    #a_by_occ=fill!(a_by_occ, 1)
    a_by_occ_initial=collect(range(1,1.1,length=n_occ-1))
    α=.089
    β=0.5
    η=.103
    σ=0.5#β*η
    μ=.714
    ϕ=0.999
    #θ=3
    γ=0.99
    κ=1.1
end

# Update model parameters, if required:
γ = .875
λf = .895
κ = .55
theta = 1.4

# Distribution of abilities:
dist=Frechet(theta,1)

τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')

# to match shares for Other for women
τ_w_opt=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')

#τ_w[:,1]=fill!(τ_w[:,1],-3)
#τ_w[1,1] = -3+3*.1

τ_w[:,2]=fill!(τ_w[:,2],0)
# scale of τ_w
# guess for initial relative τ_w for women (levels will be adjusted later)
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λm=d["λm"]
#else
    # λm=1. # to match mass of men in Teaching
#end
# scale of τ_w
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λf=d["λf"]
#else
    # λf=1.1# to match mass of women in Teaching
#end

τ_w[:,2]=ones(n_occ-1).-λm*(ones(n_occ-1).-τ_w[:,2])

#τ_w[1,2] = -4
τ_e=fill!(τ_e,0.0)

# guess for initial relative τ_w for women (levels will be adjusted later)
# if previous == 1
#     d = load("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization.jld")#2010.jld")
#     τ_w_initial=d["τ_w_opt"][:,1]
#    τ_e_initial=d["τ_e"][:,1]
# else
if previous == 0
    τ_w_initial=fill!(τ_w[:,1],0.0)
    #τ_w[1,1] = -3+3*.1

    τ_e_initial=fill!(τ_e,0.0)
end

quantile_top = .999
quantile_bottom=.001

# Set up the grid for aggregate human capital in teaching:
n_H = 11
# Set the percentage range above / below fixed point:
range_H = 0.25
# Update the grid for abilities (if necessary):
if previous == 0
    a_grid = quantile.(dist,log10.(range(1000^(quantile_bottom),1000^(quantile_top),length(a_grid)))./3)
else
    # b_grid = quantile.(dist,log10.(range(1000^(quantile_bottom),1000^(quantile_top),length(a_grid)))./3)
    # a_grid = 1 .* b_grid + 0 .* a_grid
end

# Integral bounds
lowbnd=quantile.(dist,quantile_bottom) # set to smallest element in 'a_grid'
upbnd=quantile.(dist,quantile_top) # set to largest element in 'a_grid'

# Given a productivity vector 'a_by_occ' with productivities in non-teaching occupations, the function 'share()' calculates share of people choosing occupaton occ_id:
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,a_by_occ,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(a_by_occ[occ_id]./a_by_occ[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(dist,A_tmp.*a_grid[ia]))
    end
    spl_marg = Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    # Function returns two arguments:
    return quadgk(aa -> spl_marg(aa)*pdf(dist,aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

# Calibrate a_by_occ to match labor market shares for men in non-teaching occupations:
function calibrate_A(x)
    # share of individuals employed in occ i among all non-teaching individuals among this group
    # then all share_occ will add up to 1*n_g
    share_occ_model = Array{Float64,1}(undef,n_occ-1)
    for i_occ = 1:(n_occ-1)
        share_occ_model[i_occ]=share(i_occ,2,x,τ_w[:,2],τ_e[:,2])[1]
    end
    share_occ_model=share_occ_model./sum(share_occ_model)
    # obj fn is a gap in labor market shares for non-teaching male between model and data
    return sum(abs.((share_occ_model[1:end-1].-share_occ_data[1:end-1,2])./share_occ_data[1:end-1,2]))
end

# x0=1.0.*a_by_occ_initial
res=optimize(calibrate_A,a_by_occ_initial, show_trace=false)
a_by_occ=Optim.minimizer(res)
# println("sum of abs distance between model and data for occ shares for men is ", Optim.minimum(res))
println("sum of absolute distances between a_by_occ_initial and a_by_occ for men is ",sum(abs.(a_by_occ- a_by_occ_initial)))

# Calibrate τ_w/τ_w to match labor marker shares for women in non-teaching occupations:
function calibrate_τ(x)
    # share of individuals employed in occ i among all non-teaching individuals among this group
    # then all share_occ will add up to 1*n_g
    share_occ_model = Array{Float64,1}(undef,n_occ-1)
    for i_occ = 1:(n_occ-1)
        share_occ_model[i_occ]=share(i_occ,1,a_by_occ,x,τ_e[:,1])[1]
    end
    share_occ_model=share_occ_model./sum(share_occ_model)
    # obj fn is a gap in labor market shares for non-teaching male between model and data
    return sum(abs.((share_occ_model[1:end-1].-share_occ_data[1:end-1,1])./share_occ_data[1:end-1,1]))
end

x0=1.0.*τ_w_initial[:,1]
res=optimize(calibrate_τ,x0, show_trace=false)
τ_w_opt[:,1]=Optim.minimizer(res)
# println("sum of abs distance between model and data for occ shares for women is ", Optim.minimum(res))
println("sum of absolute distances between τ_w_initial and τ_w_opt for women is ",sum(abs.(τ_w_opt-τ_w_initial)))

# scale of τ_w
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λf=d["λf"]
#else
#    λf=160 # to match mass of women in Teaching
#end
τ_w[:,1]=ones(n_occ-1).-λf*(ones(n_occ-1).-τ_w_opt[:,1])
#ones(n_occ-1).-λf/λm*(ones(n_occ-1).-τ_w_opt[:,1])

# Compute the share of individuals employed in non-teaching occupations in each group (indexed by 'iG'):
# Note: the sum of 'share_occ' will add up to 1*n_g.
share_occ = Array{Float64,2}(undef,n_occ-1,n_g)
marginal=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
for iG in 1:n_g
    for i_occ = 1:(n_occ-1)
        share_occ[i_occ,iG]=share(i_occ,iG,a_by_occ,τ_w[:,iG],τ_e[:,iG])[1]
        #τ_w[:,iG]=ones(n_occ-1).-(ones(n_occ-1).-τ_w[:,iG])./factor
        marginal[:,i_occ,iG]=share(i_occ,iG,a_by_occ,τ_w[:,iG],τ_e[:,iG])[2]
    end
    share_occ[:,iG]=share_occ[:,iG]./sum(share_occ[:,iG])
end

s_O=μ*ϕ/(μ*ϕ+1-η)

maxiterAA=1e5
tol_a=1e-6

tolT= 2.5e-4
maxiterT = 100

# t=zeros(n_H,2)

a_T_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
a_O_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
# if previous == 1
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_τW=[-0.8 -0.8]_τE=[0.0 0.0]_A=10.0_α=0.1_β=0.3_η=0.1_σ=0.4.jld")
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results/previousParameterization.jld")
#      d = load("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization.jld")
#     a_T_thresh=d["a_T_thresh"]
#     t = d["t"]
#     HH_T = d["HH_T"]
#     H_O = d["H_O"]
# else
    if previous == 0
    for ia in 1:n_occ-1
        for iG in 1:n_g
            a_T_thresh[:,iG,ia]=a_grid.^0.7
        end
    end
    t[:,1] = zeros(n_H)
    # t[:,1]=t1[:,1]
    HH_T = H_grid
    H_O = H_grid
end
e_T = Array{Float64,3}(undef,length(a_grid),n_H,n_g)
e_O = Array{Float64,4}(undef,length(a_grid),n_H,n_g,n_occ-1)

s_T = Array{Float64,3}(undef,length(a_grid),n_H,n_g)

h_T = Array{Float64,3}(undef,length(a_grid),n_H,n_g)
h_O = Array{Float64,4}(undef,length(a_grid),n_H,n_g,n_occ-1)
ω = Array{Float64,3}(undef,length(a_grid),n_H,n_g)
ω_90_10 = Array{Float64,1}(undef,n_g+1)
der_ω = Array{Float64,3}(undef,length(a_grid),n_H,n_g)

a_T_10p = Array{Float64,1}(undef,n_g+1)
a_T_90p = Array{Float64,1}(undef,n_g+1)

spl_T = Array{Spline1D,1}(undef,n_g)
spl_O = Array{Spline1D,1}(undef,n_g)
spl_T_inv = Array{Spline1D,1}(undef,n_g)
spl_O_inv = Array{Spline1D,1}(undef,n_g)
spl_ω = Array{Spline1D,1}(undef,n_H)
spl_w = Array{Spline1D,1}(undef,n_g)
spl_h_O = Array{Spline1D,1}(undef,n_g)
spl_pdf_T = Array{Spline1D,1}(undef,n_g)
wb_T = Array{Float64,2}(undef,n_H,n_g)
wb_O = Array{Float64,2}(undef,n_H,n_g)
E_T = Array{Float64,2}(undef,n_H,n_g)
E_O = Array{Float64,3}(undef,n_occ-1,n_H,n_g)
Y_T = Array{Float64,2}(undef,n_H,n_g)
Y_O = Array{Float64,3}(undef,n_occ-1,n_H,n_g)
sum_E_O = Array{Float64,2}(undef,n_H,n_g)
sum_Y_O = Array{Float64,2}(undef,n_H,n_g)

HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,2}(undef,n_occ-1,n_g)
H_O = Array{Float64,2}(undef,n_occ-1,n_H)

EN = Array{Float64,1}(undef,n_g)
EN2 = Array{Float64,1}(undef,n_g)

f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_occ-1)
f_2 = Array{Float64,2}(undef,length(a_grid),n_g)
f_O = Array{Float64,2}(undef,length(a_grid),n_g)

spl_T_thresh=Array{Spline1D,1}(undef,1)
# spl_marg = Array{Spline1D,1}(undef,1)

f1 = Array{Float64,2}(undef,n_occ-1,n_g)
mass_O = Array{Float64,2}(undef,n_occ-1,n_g)
f3 = Array{Float64,2}(undef,n_occ-1,n_g)
# spl_inv=Array{Spline1D,1}(undef,1)
f_1_T = Array{Float64,2}(undef,length(a_grid),n_g)
spl_f_1_T=Array{Spline1D,1}(undef,n_g)
f_T = Array{Float64,2}(undef,length(a_grid),n_g)
spl_f_T=Array{Spline1D,1}(undef,n_g)
mass_T = Array{Float64,1}(undef,n_g)
f2 = Array{Float64,1}(undef,n_g)
f4 = Array{Float64,1}(undef,n_g)
ratio = Array{Float64,1}(undef,n_g)
f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
a_T_90_10 = Array{Float64,1}(undef,n_g)

n_h_T_tmp=10000
h_T_tmp = Array{Float64,3}(undef,n_h_T_tmp,n_H,n_g)
for iH in 1:n_H
    for iG in 1:n_g
        h_T_tmp[:,iH,iG]=collect(range(1e-5,5,length=n_h_T_tmp)) #(η^η.*a_grid.^α.*s_O^ϕ.*(2*H_grid[1]/M)^σ).^(1/(1-η))
    end
end

########### EXOGENOUS WAGE PROFILE
ω_fn(h_T) = κ.*h_T.^γ
# derivative of wage profile
der_ω_fn(h_T) = κ.*γ.*h_T.^(γ-1)

#=spl_ω = Array{Spline1D,1}(undef,n_g)
for iH in 1:n_H
    for iG in 1:n_g
        spl_ω[iG] = Spline1D(h_T_tmp[:,iH,iG], ω_tmp[:,iH,iG],k=1, bc="extrapolate")
    end
end
spl_der_ω = Array{Spline1D,1}(undef,n_g)
der_ω(k,iG)=Dierckx.derivative(spl_ω[iG],k)
for iH in 1:n_H
    for iG in 1:n_g
        spl_der_ω[iG]=Spline1D(h_T_tmp[:,iH,iG],der_ω.(h_T_tmp[:,iH,iG],iG),k=1, bc="extrapolate")
    end
end=#

# Initial guess for fixed point of HH_T:
# Set the array index to the midpoint of H_grid:

if previous == 0
    HH_fp = 1 # HH_fp1
    iH = convert(Int,ceil(n_H/2))
    # else
    # HH_fp = HH_T[iH]
end
convHH = 1
tolHH = .0025
println("t[1,1] = ",t[1,1])
while convHH > tolHH
    global HH_fp
    global convHH
    global tolHH
    global HH_T
    println("HH_fp = ", HH_fp)
    println("t[1,1] = ",t[1,1])
    println("")
    HH_T[iHH] = HH_fp
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        # println("    T: iteration ",iterT)
        # println("HH_T[iHH]=",round(HH_T[iHH],digits=3))
        # println("t[1,:]=",round.(t[1,:],digits=3))
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k)=(η^η*(1-t[1,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iHH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iHH,iG]),maximum(h_T_tmp[:,iHH,iG])))
                h_T[ia,iHH,iG]=hh_T
                s_T[ia,iHH,iG]=fn_s_T(hh_T)
                e_T[ia,iHH,iG]=η*(1-t[1,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iHH,iG]=ω_fn(hh_T)
                der_ω[ia,iHH,iG]=der_ω_fn(hh_T)
                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])/(1-τ_w[i_occ,iG])*(1+τ_e[i_occ,iG])^(-(1-η))*((1-s_T[ia,iHH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iHH,iG]/s_O)^ϕ*(der_ω[ia,iHH,iG]/a_by_occ[i_occ])*((ω[ia,iHH,iG]/(der_ω[ia,iHH,iG]*h_T[ia,iHH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iHH,iG,i_occ]=(η^η*(1-t[1,1])^η*(1-τ_w[i_occ,iG])^η/(1+τ_e[i_occ,iG])^η*a_by_occ[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iHH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iHH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg = Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[1,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*a_by_occ[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iHH,iG]=a_by_occ[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iHH,iG]=a_by_occ[i_occ]*H_O_0[i_occ,iG]

                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    # fract of teach occ by occ, given a_T
                    f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
                    #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
                end
            end

            for ia in 1:length(a_grid)
                f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
            end
            spl_f_1_T[iG]=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iHH,iG])

            f2[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*spl_wa(aa),lowbnd,upbnd)[1]

            mass_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa),lowbnd,upbnd)[1]

            # Compute HHH_T:
            HHH_T_0[iG]= ( (1-t[1,1])^η*η^η*(2*HH_fp/M)^σ)^(β/σ/(1-η))*f2[iG]
            #*s_T^ϕ*(β/σ)^η*(sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iHH,iG]=f4[iG]
            # Total output of teachers in each group:
            Y_T[iHH,iG] = sum(H_O_0[:,iG].*a_by_occ)
            #η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])
        end
        
        HH_T[iHH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iHH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iHH,iG]=sum(E_O[:,iHH,iG])
            sum_Y_O[iHH,iG]=sum(Y_O[:,iHH,iG])
        end
        t[1,2]=sum(E_T[iHH,:].*gm)/(sum(E_T[iHH,:].*gm)+sum(sum_E_O[iHH,:].*gm))
        #Y=sum(Y_T[iHH,:].*gm)+sum(sum_Y_O[iHH,:].*gm)
        convT = abs(t[1,2]-t[1,1])
         println("convT=",round(convT,digits=3))
         println("t[1,:]=",round.(t[1,:],digits=3))
        # println("HH_T[iHH]=",round(HH_T[iHH],digits=3))
        t[1,1] =  (1-ν)*t[1,1] + ν*t[1,2]# min(t[1,2],1-1e-3)

        iterT = iterT+1
    end
    convHH = abs(log(HH_T[iHH]/HH_fp))
    println(convHH)
    println(HH_T[iHH])
    HH_fp = (1-ν2)*HH_fp + ν2*HH_T[iHH]
end
println("Found fixed point for human capital in teaching!")

# Construct grid for H_T (aggregate human capital in in teaching) around fixed point:
H_grid = HH_fp .* exp.(range(-range_H,range_H,n_H))
for iH in 1:n_H
    println("H: gridpoint ", iH," of ", n_H," (H = ", round(H_grid[iH];digits=6),")")
    println("")
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Initiate 'while' loop ('not indexed') over the tax rate:
    HH_T[iH] = HH_fp
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        # println("    T: iteration ",iterT)
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        # println("t[1,:]=",round.(t[1,:],digits=3))
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k)=(η^η*(1-t[1,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
                h_T[ia,iH,iG]=hh_T
                s_T[ia,iH,iG]=fn_s_T(hh_T)
                e_T[ia,iH,iG]=η*(1-t[1,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iH,iG]=ω_fn(hh_T)
                der_ω[ia,iH,iG]=der_ω_fn(hh_T)
                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])/(1-τ_w[i_occ,iG])*(1+τ_e[i_occ,iG])^(-(1-η))*((1-s_T[ia,iH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iH,iG]/s_O)^ϕ*(der_ω[ia,iH,iG]/a_by_occ[i_occ])*((ω[ia,iH,iG]/(der_ω[ia,iH,iG]*h_T[ia,iH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iH,iG,i_occ]=(η^η*(1-t[1,1])^η*(1-τ_w[i_occ,iG])^η/(1+τ_e[i_occ,iG])^η*a_by_occ[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg = Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[1,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*a_by_occ[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG]=a_by_occ[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iH,iG]=a_by_occ[i_occ]*H_O_0[i_occ,iG]

                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    # fract of teach occ by occ, given a_T
                    f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
                    #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
                end
            end

            # Compute 10th and 90th percentile of distribution of wages:
            
            # (a) Teachers:
            # (a.1) Compute p.d.f.:
            mass_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa),lowbnd,upbnd)[1]
            f_T[:,iG] = f_1_T[:,iG].*pdf.(dist,a_grid)./mass_T[iG]
            # (a.2) Compute ability at 10th and 90th percentile for each group using c.d.f.:
            spl_f_T[iG] = Spline1D(a_grid,f_T[:,iG])
            function fn_F_T_10p(x)
                .1 - quadgk(aa -> spl_f_T[iG](aa),lowbnd,x)[1]
            end
            a_T_10p[iG] = find_zero(fn_F_T_10p,a_grid[round(Int,length(a_grid)/2)])
            function fn_F_T_90p(x)
                .9 - quadgk(aa -> spl_f_T[iG](aa),lowbnd,x)[1]
            end
            a_T_90p[iG] = find_zero(fn_F_T_90p,a_grid[round(Int,length(a_grid)/2)])
            # (a.3) Compute wages at 10th and 90th percentile of distribution for each group:
            spl_ω[iH] = Spline1D(a_grid,ω[:,iH,iG])
            ω_90_10[iG] = spl_ω[iH](a_T_90p[iG])/spl_ω[iH](a_T_10p[iG])
            a_T_90_10[iG] = a_T_90p[iG]/a_T_10p[iG]
            
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*spl_wa(aa),lowbnd,upbnd)[1]
            
            # (b) All other occupations:
            # TO BE COMPLETED!
            # Compute HHH_T:
            HHH_T_0[iG]= ( (1-t[1,1])^η*η^η*(2*HH_fp/M)^σ)^(β/σ/(1-η))*f2[iG]
            #*s_T^ϕ*(β/σ)^η*(sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iH,iG]=f4[iG]
            # Total output of teachers in each group:
            Y_T[iH,iG] = sum(H_O_0[:,iG].*a_by_occ)
            #η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( a_by_occ.^(1/(1-η)).*((1-t[1,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])
        end
        # Compute wages at 10th and 90th percentiles in teaching (pooled across groups):
        function fn_f_T_all_g(x)
            fn = 0
            for iG = 1:n_g
                fn = fn + spl_f_T[iG](x) * mass_T[iG]/sum(mass_T)
            end
            return fn
        end
        function fn_F_T_10p_all_g(x)
            .1 - quadgk(aa -> fn_f_T_all_g(aa),lowbnd,x)[1]
        end
        function fn_F_T_90p_all_g(x)
            .9 - quadgk(aa -> fn_f_T_all_g(aa),lowbnd,x)[1]
        end
        a_T_10p[n_g+1] = find_zero(fn_F_T_10p_all_g,a_grid[round(Int,length(a_grid)/2)])
        a_T_90p[n_g+1] = find_zero(fn_F_T_90p_all_g,a_grid[round(Int,length(a_grid)/2)])
        ω_90_10[n_g+1] = spl_ω[iH](a_T_90p[n_g+1])/spl_ω[iH](a_T_10p[n_g+1])
    
        HH_T[iH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
            sum_Y_O[iH,iG]=sum(Y_O[:,iH,iG])
        end
        t[1,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        #Y=sum(Y_T[iH,:].*gm)+sum(sum_Y_O[iH,:].*gm)
        convT = abs(t[1,2]-t[1,1])
         println("convT=",round(convT,digits=3))
         println("t[1,:]=",round.(t[1,:],digits=3))
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[1,1] =  (1-ν)*t[1,1] + ν*t[1,2]# min(t[1,2],1-1e-3)

        iterT = iterT+1
    end
end

#filename1 = string("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization.jld")
#save(filename1,"τ_w_opt",τ_w_opt,"a_by_occ",a_by_occ,"a_T_thresh",a_T_thresh,"t",t,"HH_T",HH_T,"H_O",H_O)

#=
## reduced form HHH_T
HHH_T_0_tmp=Array{Float64,1}(undef,n_g)
for iG in 1:n_g
    for ia in 1:length(a_grid)
        a=a_grid[ia]
        fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
        fn_h_T(k)=(η^η*(1-t[1,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
        #(2*H_grid[iHH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
        hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iHH,iG]),maximum(h_T_tmp[:,iHH,iG])))
        h_T[ia,iHH,iG]=hh_T
        s_T[ia,iHH,iG]=fn_s_T(hh_T)
    end
    spl_s=Spline1D(a_grid, s_T[:,iHH,iG])
    spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iHH,iG]), bc="extrapolate")

    for i_occ in 1:n_occ-1
        spl_marg = Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
        for ia in 1:length(a_grid)
            a = a_grid[ia]
            # fract of teach occ by occ, given a_T
            f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
        end
    end

    for ia in 1:length(a_grid)
        f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
    end
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
    
    f2[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
    HHH_T_0_tmp[iG]=(1-t[1,1])^(β*η/σ/(1-η))*f2[iG]
end
HH_T_tmp=η^(β*η/σ/(1-η-β))*(2/M)^(β/(1-η-β))*sum(HHH_T_0_tmp.*gm)^((1-η)/(1-η-β))
=#


N=zeros(length(a_grid),n_H,n_g)
n_H=1
for iH in 1:n_H
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            N[ia,iH,iG]=(M/2)/HH_fp*(h_T[ia,iH,iG])^(β/σ)
        end
    end
end


###########################################################
##################  PARAMETERIZATION  #####################
###########################################################
av_e_O=zeros(n_g)
av_e_T=zeros(n_g)
av_s_T=zeros(n_g)
av_a_rank_T=zeros(n_g)
av_N=zeros(n_g)

for iG in 1:n_g
    spl_e_T = Spline1D(a_grid,e_T[:,1,iG])
    spl_s_T = Spline1D(a_grid,s_T[:,1,iG])
    spl_N=Spline1D(a_grid, N[:,1,iG])
    av_e_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*spl_e_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_s_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*spl_s_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_a_rank_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*cdf(dist,aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_N[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T[iG](aa)*spl_N(aa),lowbnd,upbnd)[1]/mass_T[iG]
end

# Save parameterization in JLD file:
cd("./parameterization")
save("previousParameterization.jld","a_by_occ",a_by_occ,"τ_w_opt",τ_w_opt,"τ_e",τ_e,"a_T_thresh",a_T_thresh,"t",t[1,:],"HH_T",HH_T,"H_O",H_O,"HH_fp",HH_fp,"α",α,"β",β,"η",η, "σ",σ,"μ",μ,"ϕ",ϕ,"γ",γ,"κ",κ,"theta",theta,"λf",λf,"λm",λm,"iHH",iHH,"a_grid",a_grid)
cd("..")

println("____________")
println("share of teachers among women= ",mass_T[1])
println("share of teachers among men= ",mass_T[2])
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

# Write parameter values and moments to CSV file.

# Start with λf / κ (parameters) and male / female teachers' share (moments).
df = DataFrame(λf = round(λf;digits=4),share_teachers_female = round(mass_T[1];digits=4),κ = round(κ;digits=4),share_teachers_male = round(mass_T[2];digits=4),γ = round(γ;digits=4),p90_p10_ω_teachers = round(ω_90_10[end];digits=4), θ = round(theta;digits=4),p90_p10_w_other = NaN)
# If it exists, load previous parameterization (in CSV format), convert to DataFrame, and append 'df':
if isfile("./results/moments.csv") == true
    moments = DataFrame(CSV.File("./results/moments.csv"))
    append!(moments,df)
else
    # If the CSV file doesn't exist, create a new DataFrame named "moments":
    moments = df
end
println(moments)
# Write DataFrame to CSV file:
CSV.write("./results/moments.csv",moments)