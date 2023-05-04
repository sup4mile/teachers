#using Pkg;Pkg.add("SciPy")
using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim, Roots
using PyCall#, SciPy

# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
# cd("C:/Users/julia/Desktop/Research/Teachers")
 cd("C:\\Users\\julia\\Desktop\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
#cd("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/codes_w_exo_wage")
xf = XLSX.readxlsx("C:/Users/julia/Desktop/BOX/Teachers/New Yulia's Folder/LaborMarketData/Occupation_shares_v3.xlsx")
sh = xf["mom"]
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
    println("No Agrid for the selected year")
end

# step for the update
ν = 0.1
ν2 = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)

α=.089
β=0.5
η=.103
σ=0.5#β*η
μ=.714
ϕ=0.999
θ=2
γ=1.01
κ=1.1

# guess for initial relative A for men (levels will be adjusted later)
d = load("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization_1.01.jld")
τ_w_opt=d["τ_w_opt"]
Agrid=d["Agrid"]
a_T_thresh=d["a_T_thresh"]
a_O_thresh=d["a_O_thresh"]
t=d["t"]
HH_T=d["HH_T"]
H_O=d["H_O"]

τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_w_opt=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[:,2]=fill!(τ_w[:,2],0)
λm=1. # to match mass of men in Teaching
λf=1.1# to match mass of women in Teaching
τ_w[:,2]=ones(n_occ-1).-λm*(ones(n_occ-1).-τ_w[:,2])
τ_w[:,1]=ones(n_occ-1).-λf*(ones(n_occ-1).-τ_w_opt[:,1])
τ_e=fill!(τ_e,0.0)

quantile_top = .999
quantile_bottom=.001

dist=Frechet(θ,1)
a_grid=append!(quantile.(dist,quantile_bottom:.001:.004), append!(quantile.(dist,.005:.015:.979), quantile.(dist,.980:.002:quantile_top)))

# Integral bounds
lowbnd=quantile.(dist,quantile_bottom)
upbnd=quantile.(dist,quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(dist,A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(dist,aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

# share of individuals employed in occ i among all non-teaching individuals among this group
# then all share_occ will add up to 1*n_g
share_occ = Array{Float64,2}(undef,n_occ-1,n_g)
marginal=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
for iG in 1:n_g
    for i_occ = 1:(n_occ-1)
        share_occ[i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[1]
        #τ_w[:,iG]=ones(n_occ-1).-(ones(n_occ-1).-τ_w[:,iG])./factor
        marginal[:,i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[2]
    end
    share_occ[:,iG]=share_occ[:,iG]./sum(share_occ[:,iG])
end

s_O=μ*ϕ/(μ*ϕ+1-η)

########### EXOGENOUS WAGE PROFILE
ω_fn(h_T) = κ*h_T^γ
der_ω_fn(h_T) = κ*γ*h_T^(γ-1)

f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
f_1_T=Array{Float64,2}(undef,length(a_grid),n_g)
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
        for ia in 1:length(a_grid)
            a = a_grid[ia]
            # fract of teach occ by occ, given a_T
            f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
        end
    end
    for ia in 1:length(a_grid)
        f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
    end
end

p1=plot(cdf.(dist,a_grid), [pdf.(dist,a_grid).*f_1_T[:,1] pdf.(dist,a_grid).*f_1_T[:,2]], label=["female,λ=1.01" "male,λ=1.01"])

# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
# cd("C:/Users/julia/Desktop/Research/Teachers")
 cd("C:\\Users\\julia\\Desktop\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
#cd("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/codes_w_exo_wage")
xf = XLSX.readxlsx("C:/Users/julia/Desktop/BOX/Teachers/New Yulia's Folder/LaborMarketData/Occupation_shares_v3.xlsx")
sh = xf["mom"]
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
    println("No Agrid for the selected year")
end

# step for the update
ν = 0.1
ν2 = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)

α=.089
β=0.5
η=.103
σ=0.5#β*η
μ=.714
ϕ=0.999
θ=2
γ=1.0
κ=1.1

# guess for initial relative A for men (levels will be adjusted later)
d = load("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization_1.0.jld")
τ_w_opt=d["τ_w_opt"]
Agrid=d["Agrid"]
a_T_thresh=d["a_T_thresh"]
a_O_thresh=d["a_O_thresh"]
t=d["t"]
HH_T=d["HH_T"]
H_O=d["H_O"]

τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_w_opt=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[:,2]=fill!(τ_w[:,2],0)
λm=1. # to match mass of men in Teaching
λf=1.1# to match mass of women in Teaching
τ_w[:,2]=ones(n_occ-1).-λm*(ones(n_occ-1).-τ_w[:,2])
τ_w[:,1]=ones(n_occ-1).-λf*(ones(n_occ-1).-τ_w_opt[:,1])
τ_e=fill!(τ_e,0.0)

quantile_top = .999
quantile_bottom=.001

dist=Frechet(θ,1)
a_grid=append!(quantile.(dist,quantile_bottom:.001:.004), append!(quantile.(dist,.005:.015:.979), quantile.(dist,.980:.002:quantile_top)))

# Integral bounds
lowbnd=quantile.(dist,quantile_bottom)
upbnd=quantile.(dist,quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(dist,A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(dist,aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

# share of individuals employed in occ i among all non-teaching individuals among this group
# then all share_occ will add up to 1*n_g
share_occ = Array{Float64,2}(undef,n_occ-1,n_g)
marginal=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
for iG in 1:n_g
    for i_occ = 1:(n_occ-1)
        share_occ[i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[1]
        #τ_w[:,iG]=ones(n_occ-1).-(ones(n_occ-1).-τ_w[:,iG])./factor
        marginal[:,i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[2]
    end
    share_occ[:,iG]=share_occ[:,iG]./sum(share_occ[:,iG])
end

s_O=μ*ϕ/(μ*ϕ+1-η)

########### EXOGENOUS WAGE PROFILE
ω_fn(h_T) = κ*h_T^γ
der_ω_fn(h_T) = κ*γ*h_T^(γ-1)

f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
f_1_T=Array{Float64,2}(undef,length(a_grid),n_g)
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
        for ia in 1:length(a_grid)
            a = a_grid[ia]
            # fract of teach occ by occ, given a_T
            f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
        end
    end
    for ia in 1:length(a_grid)
        f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
    end
end

p2=plot(p1,cdf.(dist,a_grid), [pdf.(dist,a_grid).*f_1_T[:,1] pdf.(dist,a_grid).*f_1_T[:,2]], label=["female,λ=1.0" "male,λ=1.0"])



# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
# cd("C:/Users/julia/Desktop/Research/Teachers")
 cd("C:\\Users\\julia\\Desktop\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
#cd("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/codes_w_exo_wage")
xf = XLSX.readxlsx("C:/Users/julia/Desktop/BOX/Teachers/New Yulia's Folder/LaborMarketData/Occupation_shares_v3.xlsx")
sh = xf["mom"]
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
    println("No Agrid for the selected year")
end

# step for the update
ν = 0.1
ν2 = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)

α=.089
β=0.5
η=.103
σ=0.5#β*η
μ=.714
ϕ=0.999
θ=2
γ=0.99
κ=1.1

# guess for initial relative A for men (levels will be adjusted later)
d = load("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization_0.99.jld")
τ_w_opt=d["τ_w_opt"]
Agrid=d["Agrid"]
a_T_thresh=d["a_T_thresh"]
a_O_thresh=d["a_O_thresh"]
t=d["t"]
HH_T=d["HH_T"]
H_O=d["H_O"]

τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_w_opt=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[:,2]=fill!(τ_w[:,2],0)
λm=1. # to match mass of men in Teaching
λf=1.1# to match mass of women in Teaching
τ_w[:,2]=ones(n_occ-1).-λm*(ones(n_occ-1).-τ_w[:,2])
τ_w[:,1]=ones(n_occ-1).-λf*(ones(n_occ-1).-τ_w_opt[:,1])
τ_e=fill!(τ_e,0.0)

quantile_top = .999
quantile_bottom=.001

dist=Frechet(θ,1)
a_grid=append!(quantile.(dist,quantile_bottom:.001:.004), append!(quantile.(dist,.005:.015:.979), quantile.(dist,.980:.002:quantile_top)))

# Integral bounds
lowbnd=quantile.(dist,quantile_bottom)
upbnd=quantile.(dist,quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(dist,A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(dist,aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

# share of individuals employed in occ i among all non-teaching individuals among this group
# then all share_occ will add up to 1*n_g
share_occ = Array{Float64,2}(undef,n_occ-1,n_g)
marginal=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
for iG in 1:n_g
    for i_occ = 1:(n_occ-1)
        share_occ[i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[1]
        #τ_w[:,iG]=ones(n_occ-1).-(ones(n_occ-1).-τ_w[:,iG])./factor
        marginal[:,i_occ,iG]=share(i_occ,iG,Agrid,τ_w[:,iG],τ_e[:,iG])[2]
    end
    share_occ[:,iG]=share_occ[:,iG]./sum(share_occ[:,iG])
end

s_O=μ*ϕ/(μ*ϕ+1-η)

########### EXOGENOUS WAGE PROFILE
ω_fn(h_T) = κ*h_T^γ
der_ω_fn(h_T) = κ*γ*h_T^(γ-1)

f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
f_1_T=Array{Float64,2}(undef,length(a_grid),n_g)
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
        for ia in 1:length(a_grid)
            a = a_grid[ia]
            # fract of teach occ by occ, given a_T
            f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
        end
    end
    for ia in 1:length(a_grid)
        f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
    end
end

p3=plot(p2,cdf.(dist,a_grid), [pdf.(dist,a_grid).*f_1_T[:,1] pdf.(dist,a_grid).*f_1_T[:,2]], label=["female,λ=0.99" "male,λ=0.99"], title="Density of Teachers")

filename1 = string("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/density.png")
savefig(p3,filename1)

















# Initial guess for fixed point of HH_T:
# Set the array index to the midpoint of H_grid:
iH = convert(Int,ceil(H_grid_length/2))
if previous == 0
    HH_fp = 1 # HH_fp1
else
    HH_fp = HH_T[iH]
end
gapHH = 1
tolHH = .0025
while gapHH > tolHH
    global HH_fp
    global gapHH
    global tolHH
    global HH_T
    println("HH_fp = ", HH_fp)
    println("")
    HH_T[iH] = HH_fp
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        # println("    T: iteration ",iterT)
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        # println("t[iH,:]=",round.(t[iH,:],digits=3))
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k)=(η^η*(1-t[iH,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
                h_T[ia,iH,iG]=hh_T
                s_T[ia,iH,iG]=fn_s_T(hh_T)
                e_T[ia,iH,iG]=η*(1-t[iH,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iH,iG]=ω_fn(hh_T)
                der_ω[ia,iH,iG]=der_ω_fn(hh_T)
                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])/(1-τ_w[i_occ,iG])*(1+τ_e[i_occ,iG])^(-(1-η))*((1-s_T[ia,iH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iH,iG]/s_O)^ϕ*(der_ω[ia,iH,iG]/Agrid[i_occ])*((ω[ia,iH,iG]/(der_ω[ia,iH,iG]*h_T[ia,iH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iH,iG,i_occ]=(η^η*(1-t[iH,1])^η*(1-τ_w[i_occ,iG])^η/(1+τ_e[i_occ,iG])^η*Agrid[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]

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
            spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*spl_wa(aa),lowbnd,upbnd)[1]

            mass_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]

            # Compute HHH_T:
            HHH_T_0[iG]= ( (1-t[iH,1])^η*η^η*(2*HH_fp/M)^σ)^(β/σ/(1-η))*f2[iG]
            #*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iH,iG]=f4[iG]
            # Total output of teachers in each group:
            Y_T[iH,iG] = sum(H_O_0[:,iG].*Agrid)
            #η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])
        end
        
        HH_T[iH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
            sum_Y_O[iH,iG]=sum(Y_O[:,iH,iG])
        end
        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        #Y=sum(Y_T[iH,:].*gm)+sum(sum_Y_O[iH,:].*gm)
        convT = abs(t[iH,2]-t[iH,1])
         println("convT=",round(convT,digits=3))
         println("t[iH,:]=",round.(t[iH,:],digits=3))
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[iH,1] =  ν*t[iH,1]+(1-ν)*t[iH,2]# min(t[iH,2],1-1e-3)

        iterT = iterT+1
    end
    gapHH = abs(log(HH_T[iH]/HH_fp))
    println(gapHH)
    println(HH_T[iH])
    HH_fp = ν2*HH_T[iH] + (1-ν2)*HH_fp
end
println("Found fixed point for human capital in teaching!")

# Update H_grid:
H_grid = HH_fp#collect(range((1-H_grid_range)*HH_fp,stop=(1+H_grid_range)*HH_fp,length=H_grid_length))
H_grid_length=1
for iH in 1:H_grid_length
    println("H: gridpoint ", iH," of ", H_grid_length," (H = ", H_grid[iH],")")
    println("")
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Initiate 'while' loop ('not indexed') over the tax rate:
    HH_T[iH] = HH_fp
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        # println("    T: iteration ",iterT)
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        # println("t[iH,:]=",round.(t[iH,:],digits=3))
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
                fn_h_T(k)=(η^η*(1-t[iH,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
                h_T[ia,iH,iG]=hh_T
                s_T[ia,iH,iG]=fn_s_T(hh_T)
                e_T[ia,iH,iG]=η*(1-t[iH,1])*der_ω_fn(hh_T)*hh_T
                ω[ia,iH,iG]=ω_fn(hh_T)
                der_ω[ia,iH,iG]=der_ω_fn(hh_T)
                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])/(1-τ_w[i_occ,iG])*(1+τ_e[i_occ,iG])^(-(1-η))*((1-s_T[ia,iH,iG])/(1-s_O))^((1-η)/μ)*(s_T[ia,iH,iG]/s_O)^ϕ*(der_ω[ia,iH,iG]/Agrid[i_occ])*((ω[ia,iH,iG]/(der_ω[ia,iH,iG]*h_T[ia,iH,iG])-η)/(1-η))^(1-η))^(1/α)*a
                    h_O[ia,iH,iG,i_occ]=(η^η*(1-t[iH,1])^η*(1-τ_w[i_occ,iG])^η/(1+τ_e[i_occ,iG])^η*Agrid[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(dist,aa)*cdf(dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]

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
            spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*spl_wa(aa),lowbnd,upbnd)[1]

            mass_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]

            # Compute HHH_T:
            HHH_T_0[iG]= ( (1-t[iH,1])^η*η^η*(2*HH_fp/M)^σ)^(β/σ/(1-η))*f2[iG]
            #*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iH,iG]=f4[iG]
            # Total output of teachers in each group:
            Y_T[iH,iG] = sum(H_O_0[:,iG].*Agrid)
            #η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])
        end
        HH_T[iH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
            sum_Y_O[iH,iG]=sum(Y_O[:,iH,iG])
        end
        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        #Y=sum(Y_T[iH,:].*gm)+sum(sum_Y_O[iH,:].*gm)
        convT = abs(t[iH,2]-t[iH,1])
         println("convT=",round(convT,digits=3))
         println("t[iH,:]=",round.(t[iH,:],digits=3))
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[iH,1] =  ν*t[iH,1]+(1-ν)*t[iH,2]# min(t[iH,2],1-1e-3)

        iterT = iterT+1
    end
end

filename1 = string("C:/Users/julia/Documents/GitHub/teachers/julia/codes_w_exo_wage/parameterization/previousParameterization_1.01.jld")
save(filename1,"τ_w_opt",τ_w_opt,"Agrid",Agrid,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"t",t,"HH_T",HH_T,"H_O",H_O)

#=
## reduced form HHH_T
HHH_T_0_tmp=Array{Float64,1}(undef,n_g)
for iG in 1:n_g
    for ia in 1:length(a_grid)
        a=a_grid[ia]
        fn_s_T(k)=μ*ϕ*der_ω_fn(k)*k/((μ*ϕ-η)*der_ω_fn(k)*k+ω_fn(k))
        fn_h_T(k)=(η^η*(1-t[iH,1])^η*der_ω_fn(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
        #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
        hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
        h_T[ia,iH,iG]=hh_T
        s_T[ia,iH,iG]=fn_s_T(hh_T)
    end
    spl_s=Spline1D(a_grid, s_T[:,iH,iG])
    spl_dw=Spline1D(a_grid, der_ω_fn.(h_T[:,iH,iG]), bc="extrapolate")

    for i_occ in 1:n_occ-1
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
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
    HHH_T_0_tmp[iG]=(1-t[iH,1])^(β*η/σ/(1-η))*f2[iG]
end
HH_T_tmp=η^(β*η/σ/(1-η-β))*(2/M)^(β/(1-η-β))*sum(HHH_T_0_tmp.*gm)^((1-η)/(1-η-β))
=#


N=zeros(length(a_grid),H_grid_length,n_g)
H_grid_length=1
for iH in 1:H_grid_length
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
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
    spl_N=Spline1D(a_grid, N[:,1,iG])
    av_e_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*spl_e_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_s_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*spl_s_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_a_rank_T[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*cdf(dist,aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_N[iG]=quadgk(aa -> pdf(dist,aa)*spl_f_1_T(aa)*spl_N(aa),lowbnd,upbnd)[1]/mass_T[iG]
end

println("____________")
println("share of teachers among women= ",mass_T[1])
println("share of teachers among men= ",mass_T[2])
println("average class size= ",sum(av_N.*gm))
println(" ")
println("average a_rank_T_female= ",av_a_rank_T[1])
println("average a_rank_T_male= ",av_a_rank_T[2])
println("average e_T_female= ",av_e_T[1])
println("average e_T_male= ",av_e_T[2])
println("average s_T_female= ",av_s_T[1])
println("average s_T_male= ",av_s_T[2])
println("average w_T_female= ",E_T[1,1]/mass_T[1])
println("average w_T_male= ",E_T[1,2]/mass_T[2])
println("dispersion w_T_female= ")
println("dispersion w_T_male= ")
println(" ")
#println("average a_rank_O_female= ",)
#println("average a_rank_O_male= ",)
#println("average e_O_female= ",)
#println("average e_O_male= ",)
println("average s_O_female= ",s_O)
println("average s_O_male= ",s_O)
#println("average w_O_female= ",)
#println("average w_O_male= ",)
#println("dispersion w_O_female= ",)
#println("dispersion w_O_male= ",)
println(" ")
println("tax rate= ",t[1,1])
println("output= ",sum(Y_T[iH,:].*gm)+sum(sum_Y_O[iH,:].*gm))
println("HH_fp= ",HH_fp)
