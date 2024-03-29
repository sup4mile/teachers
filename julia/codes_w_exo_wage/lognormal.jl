#using Pkg;Pkg.add("SciPy")
using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim, Roots
using PyCall#, SciPy

# Load occupational threshold from previous parameterization:
previous = 0;
# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
cd("C:/Users/julia/Desktop/Research/Teachers")
#cd("C:\\Users\\julia\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
xf = XLSX.readxlsx("Occupation_shares_v3.xlsx")
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
#θ=3

# guess for initial relative A for men (levels will be adjusted later)
if previous == 1
    d = load("C:/Users/julia/Desktop/teacher/previousParameterization_.jld")#2010.jld")
    Agrid_initial=d["Agrid"]
else
    # productivity in 'Other' occupations
    Agrid_initial=Array{Float64,1}(undef,n_occ-1)
    #Agrid=fill!(Agrid, 1)
    Agrid_initial=collect(range(1,1.1,length=n_occ-1))
end

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
    λm=1. # to match mass of men in Teaching
#end
# scale of τ_w
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λf=d["λf"]
#else
    λf=1.1# to match mass of women in Teaching
#end

τ_w[:,2]=ones(n_occ-1).-λm*(ones(n_occ-1).-τ_w[:,2])

#τ_w[1,2] = -4
τ_e=fill!(τ_e,0.0)

# guess for initial relative τ_w for women (levels will be adjusted later)
if previous == 1
    d = load("C:/Users/julia/Desktop/teacher/previousParameterization_.jld")#2010.jld")
    τ_w_initial=d["τ_w_opt"][:,1]
#    τ_e_initial=d["τ_e"][:,1]
else
    τ_w_initial=fill!(τ_w[:,1],0.0)
    #τ_w[1,1] = -3+3*.1

    τ_e_initial=fill!(τ_e,0.0)
end

quantile_top = .999
quantile_bottom=.001

# Set up the grid for aggregate human capital in teaching:
H_grid_length = 21
# Set the percentage range above / below fixed point:
H_grid_range = 0.2
H_grid = collect(range((1-H_grid_range)*1,stop=(1+H_grid_range)*1,length=H_grid_length))


# a_grid=quantile.(Frechet(θ),.005:.015:1)
mean_a=0
std_a=1

dist=Normal(mean_a,std_a)
log_dist=LogNormal(mean_a,std_a)

a_grid=exp.(append!(quantile.(dist,quantile_bottom:.001:.004), append!(quantile.(dist,.005:.015:.979), quantile.(dist,.980:.002:quantile_top))))

# Integral bounds
lowbnd=quantile.(log_dist,quantile_bottom)
upbnd=quantile.(log_dist,quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(log_dist,A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(log_dist,aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

###### Calibrate Agrid to match labor market shares for men in Other ######
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

x0=1.0.*Agrid_initial
res=optimize(calibrate_A,x0, show_trace=false)
Agrid=Optim.minimizer(res)
println("sum of abs distance between model and data for occ shares for men is ", Optim.minimum(res))
##############

###### Calibrate τ_w/τ_w to match labor marker shares for women in Other ######
function calibrate_τ(x)
    # share of individuals employed in occ i among all non-teaching individuals among this group
    # then all share_occ will add up to 1*n_g
    share_occ_model = Array{Float64,1}(undef,n_occ-1)
    for i_occ = 1:(n_occ-1)
        share_occ_model[i_occ]=share(i_occ,1,Agrid,x,τ_e[:,1])[1]
    end
    share_occ_model=share_occ_model./sum(share_occ_model)
    # obj fn is a gap in labor market shares for non-teaching male between model and data
    return sum(abs.((share_occ_model[1:end-1].-share_occ_data[1:end-1,1])./share_occ_data[1:end-1,1]))
end

x0=1.0.*τ_w_initial
res=optimize(calibrate_τ,x0, show_trace=false)
τ_w_opt[:,1]=Optim.minimizer(res)
println("sum of abs distance between model and data for occ shares for women is ", Optim.minimum(res))
##############

# scale of τ_w
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λf=d["λf"]
#else
#    λf=160 # to match mass of women in Teaching
#end
τ_w[:,1]=ones(n_occ-1).-λf*(ones(n_occ-1).-τ_w_opt[:,1])
#ones(n_occ-1).-λf/λm*(ones(n_occ-1).-τ_w_opt[:,1])

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

maxiterAA=1e5
tol_a=1e-6

tolT= 2.5e-4
maxiterT = 100

t=zeros(H_grid_length,2)

a_T_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
a_O_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
if previous == 1
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_τW=[-0.8 -0.8]_τE=[0.0 0.0]_A=10.0_α=0.1_β=0.3_η=0.1_σ=0.4.jld")
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/previousParameterization.jld")
    d = load("C:/Users/iuliia/Desktop/teachers/julia/results_new/previousParameterization.jld")
    a_T_thresh=d["a_T_thresh"]
    t = d["t"]
    HH_T = d["HH_T"]
    H_O = d["H_O"]
else
    for ia in 1:n_occ-1
        for iG in 1:n_g
            a_T_thresh[:,iG,ia]=a_grid.^0.7
        end
    end
    t[:,1] = zeros(H_grid_length)
    t[:,1]=t1[:,1]
    HH_T = H_grid
    H_O = H_grid
end
e_T = Array{Float64,3}(undef,length(a_grid),H_grid_length,n_g)
e_O = Array{Float64,4}(undef,length(a_grid),H_grid_length,n_g,n_occ-1)

s_T = Array{Float64,3}(undef,length(a_grid),H_grid_length,n_g)

h_T = Array{Float64,3}(undef,length(a_grid),H_grid_length,n_g)
h_O = Array{Float64,4}(undef,length(a_grid),H_grid_length,n_g,n_occ-1)
ω = Array{Float64,3}(undef,length(a_grid),H_grid_length,n_g)

spl_T = Array{Spline1D,1}(undef,n_g)
spl_O = Array{Spline1D,1}(undef,n_g)
spl_T_inv = Array{Spline1D,1}(undef,n_g)
spl_O_inv = Array{Spline1D,1}(undef,n_g)
spl_ω = Array{Spline1D,1}(undef,n_g)
spl_w = Array{Spline1D,1}(undef,n_g)
spl_h_O = Array{Spline1D,1}(undef,n_g)
spl_pdf_T = Array{Spline1D,1}(undef,n_g)

wb_T = Array{Float64,2}(undef,H_grid_length,n_g)
wb_O = Array{Float64,2}(undef,H_grid_length,n_g)
E_T = Array{Float64,2}(undef,H_grid_length,n_g)
E_O = Array{Float64,3}(undef,n_occ-1,H_grid_length,n_g)
Y_T = Array{Float64,2}(undef,H_grid_length,n_g)
Y_O = Array{Float64,3}(undef,n_occ-1,H_grid_length,n_g)
sum_E_O = Array{Float64,2}(undef,H_grid_length,n_g)
sum_Y_O = Array{Float64,2}(undef,H_grid_length,n_g)

HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,2}(undef,n_occ-1,n_g)
H_O = Array{Float64,2}(undef,n_occ-1,H_grid_length)

EN = Array{Float64,1}(undef,n_g)
EN2 = Array{Float64,1}(undef,n_g)

f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_occ-1)
f_2 = Array{Float64,2}(undef,length(a_grid),n_g)
f_O = Array{Float64,2}(undef,length(a_grid),n_g)

spl_T_thresh=Array{Spline1D,1}(undef,1)
spl_marg=Array{Spline1D,1}(undef,1)

f1 = Array{Float64,2}(undef,n_occ-1,n_g)
mass_O = Array{Float64,2}(undef,n_occ-1,n_g)
f3 = Array{Float64,2}(undef,n_occ-1,n_g)
spl_inv=Array{Spline1D,1}(undef,1)
spl_f_1_T=Array{Spline1D,1}(undef,1)
f_1_T = Array{Float64,2}(undef,length(a_grid),n_g)
f_T = Array{Float64,2}(undef,length(a_grid),n_g)
mass_T = Array{Float64,1}(undef,n_g)
f2 = Array{Float64,1}(undef,n_g)
f4 = Array{Float64,1}(undef,n_g)
ratio = Array{Float64,1}(undef,n_g)
f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)

n_h_T_tmp=10000
h_T_tmp = Array{Float64,3}(undef,n_h_T_tmp,H_grid_length,n_g)
for iH in 1:H_grid_length
    for iG in 1:n_g
        h_T_tmp[:,iH,iG]=collect(range(1e-5,5,length=n_h_T_tmp)) #(η^η.*a_grid.^α.*s_O^ϕ.*(2*H_grid[1]/M)^σ).^(1/(1-η))
    end
end

# EXOGENOUS WAGE PROFILE
ω_tmp = Array{Float64,3}(undef,n_h_T_tmp,H_grid_length,n_g)
ω_tmp =0.2*h_T_tmp.+0.8.*h_T_tmp.^0.8
spl_ω = Array{Spline1D,1}(undef,n_g)
for iH in 1:H_grid_length
    for iG in 1:n_g
        spl_ω[iG] = Spline1D(h_T_tmp[:,iH,iG], ω_tmp[:,iH,iG],k=1, bc="extrapolate")
    end
end

# derivative of wage profile
spl_der_ω = Array{Spline1D,1}(undef,n_g)
der_ω(k,iG)=Dierckx.derivative(spl_ω[iG],k)
for iH in 1:H_grid_length
    for iG in 1:n_g
        spl_der_ω[iG]=Spline1D(h_T_tmp[:,iH,iG],der_ω.(h_T_tmp[:,iH,iG],iG),k=1, bc="extrapolate")
    end
end

# Initial guess for fixed point of HH_T:
# Set the array index to the midpoint of H_grid:
iH = convert(Int,ceil(H_grid_length/2))
if previous == 0
    HH_fp =HH_fp1# 1
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
                dw_dh(k)=spl_der_ω[iG](k)
                w_h(k)=spl_ω[iG](k)
                fn_s_T(k)=μ*ϕ*dw_dh(k)*k/((μ*ϕ-η)*dw_dh(k)*k+w_h(k))
                fn_h_T(k)=(η^η*(1-t[iH,1])^η*dw_dh(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
                h_T[ia,iH,iG]=hh_T
                s_T[ia,iH,iG]=fn_s_T(hh_T)
                e_T[ia,iH,iG]=η*(1-t[iH,1])*dw_dh(hh_T)*hh_T
                ω[ia,iH,iG]=w_h(hh_T)

                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])*((1-s_O)/(1-s_T[ia,iH,iG]))^μ*(((1-t[iH,1])*ω[ia,iH,iG]-e_T[ia,iH,iG])/ ((1-η)*(η^η*(1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG])*Agrid[i_occ]*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))) ))^((1-η)/α)
                    # ((1+τ_e[i_occ,iG])*((1-s_O)/(1-s_T[ia,iH,iG]))^μ*((ω[ia,iH,iG]-e_T[ia,iH,iG])/ ((1-η)*(η^η*τ_w[i_occ,iG]*Agrid[i_occ]*s_O*(2*H_grid[iH]/M)^σ)^(1/(1-η))) ))^((1-η)/α)
                    h_O[ia,iH,iG,i_occ]=(η^η*(1-t[iH,1])^η*(1-τ_w[iH,1])^η/(1+τ_e[iH,1])^η*Agrid[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, spl_der_ω[iG].(h_T[:,iH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(log_dist,aa)*cdf(log_dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(log_dist,aa)*cdf(log_dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]

                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    # fract of teach occ by occ, given a_T
                    f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(log_dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
                    #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(log_dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
                end
            end

            for ia in 1:length(a_grid)
                f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
            end
            spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*spl_wa(aa),lowbnd,upbnd)[1]

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
                dw_dh(k)=spl_der_ω[iG](k)
                w_h(k)=spl_ω[iG](k)
                fn_s_T(k)=μ*ϕ*dw_dh(k)*k/((μ*ϕ-η)*dw_dh(k)*k+w_h(k))
                fn_h_T(k)=(η^η*(1-t[iH,1])^η*dw_dh(k)^η*a^α*fn_s_T(k)^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))-k
                #(2*H_grid[iH]/M)^σ*a_grid[ia]^α*fn_s_T(k)^ϕ*fn_e_T(k)^η-k
                hh_T=find_zero(fn_h_T,(minimum(h_T_tmp[:,iH,iG]),maximum(h_T_tmp[:,iH,iG])))
                h_T[ia,iH,iG]=hh_T
                s_T[ia,iH,iG]=fn_s_T(hh_T)
                e_T[ia,iH,iG]=η*(1-t[iH,1])*dw_dh(hh_T)*hh_T
                ω[ia,iH,iG]=w_h(hh_T)

                for i_occ in 1:n_occ-1
                    a_O_thresh[ia,iG,i_occ]=((1+τ_e[i_occ,iG])*((1-s_O)/(1-s_T[ia,iH,iG]))^μ*(((1-t[iH,1])*ω[ia,iH,iG]-e_T[ia,iH,iG])/ ((1-η)*(η^η*(1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG])*Agrid[i_occ]*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))) ))^((1-η)/α)
                    # ((1+τ_e[i_occ,iG])*((1-s_O)/(1-s_T[ia,iH,iG]))^μ*((ω[ia,iH,iG]-e_T[ia,iH,iG])/ ((1-η)*(η^η*τ_w[i_occ,iG]*Agrid[i_occ]*s_O*(2*H_grid[iH]/M)^σ)^(1/(1-η))) ))^((1-η)/α)
                    h_O[ia,iH,iG,i_occ]=(η^η*(1-t[iH,1])^η*(1-τ_w[iH,1])^η/(1+τ_e[iH,1])^η*Agrid[i_occ]^η*a_grid[ia]^α*s_O^ϕ*(2*HH_fp/M)^σ)^(1/(1-η))
                end
            end
            spl_s=Spline1D(a_grid, s_T[:,iH,iG])
            spl_dw=Spline1D(a_grid, spl_der_ω[iG].(h_T[:,iH,iG]), bc="extrapolate")

            for i_occ in 1:n_occ-1
                spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG], bc="extrapolate")
                # inverse threshold
                spl_inv = Spline1D(a_O_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
                f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(log_dist,aa)*cdf(log_dist,spl_inv(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

                f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(log_dist,aa)*cdf(log_dist,spl_inv(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
                # Compute H_O:
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]
                # Total output of others in each group:
                Y_O[i_occ,iH,iG]=Agrid[i_occ]*H_O_0[i_occ,iG]

                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    # fract of teach occ by occ, given a_T
                    f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(log_dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
                    #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(log_dist,aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
                end
            end

            for ia in 1:length(a_grid)
                f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
            end
            spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])
            spl_wa=Spline1D(a_grid, ω[:,iH,iG])

            f2[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*aa^(α*β/σ/(1-η))*spl_s(aa)^(ϕ*β/σ/(1-η))*spl_dw(aa)^(η*β/σ/(1-η)),lowbnd,upbnd)[1]
            f4[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*spl_wa(aa),lowbnd,upbnd)[1]
            mass_T[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]

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
    av_e_T[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*spl_e_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_s_T[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*spl_s_T(aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_a_rank_T[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*cdf(log_dist,aa),lowbnd,upbnd)[1]/mass_T[iG]
    av_N[iG]=quadgk(aa -> pdf(log_dist,aa)*spl_f_1_T(aa)*spl_N(aa),lowbnd,upbnd)[1]/mass_T[iG]
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
