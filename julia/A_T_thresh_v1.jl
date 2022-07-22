#using Pkg;Pkg.add("SciPy")
using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings, CSV, DataFrames, LinearAlgebra, Optim
using PyCall#, SciPy

# Load occupational threshold from previous parameterization:
previous = 0;
# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
#occ = ["other1", "other2", "teaching"]
#occ = ["other1", "other2", "other3", "other4", "teaching"]
import XLSX
cd("C:\\Users\\julia\\Desktop\\teacher")
#cd("C:\\Users\\julia\\Box\\Teachers\\New Yulia's Folder\\LaborMarketData")
xf = XLSX.readxlsx("Occupation_shares.xlsx")
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
v = 0.1
ν2 = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)

α=.089
η=.103
β=.6#.4#.4#.5
σ=.4#β*η
μ=.5
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
    λm=8.5 # to match mass of men in Teaching
#end
# scale of τ_w
#if previous == 1
#    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
#    λf=d["λf"]
#else
    λf=8.# to match mass of women in Teaching
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

a_grid=exp.(append!(quantile.(Normal(mean_a,std_a),quantile_bottom:.001:.004), append!(quantile.(Normal(mean_a,std_a),.005:.015:.979), quantile.(Normal(mean_a,std_a),.980:.002:quantile_top))))

# Integral bounds
lowbnd=quantile.(LogNormal(mean_a,std_a),quantile_bottom)
upbnd=quantile.(LogNormal(mean_a,std_a),quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid,τ_w,τ_e)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id])./(ones(n_occ-2).-τ_w[1:end.!=occ_id])).*(((ones(n_occ-2).+τ_e[occ_id])./(ones(n_occ-2).+τ_e[1:end.!=occ_id])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(LogNormal(mean_a,std_a),A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG])#,bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
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

s_T=μ*ϕ/(μ*ϕ+σ/β-η)
s_O=μ*ϕ/(μ*ϕ+1-η)
cnst=(1-η)/(1-η*β/σ)*((1-s_O)/(1-s_T))^(1/μ)

maxiterAA=1e5
tol_a=1e-6

tolT= 2.5e-4
maxiterT = 100

t=zeros(H_grid_length,2)

a_T_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
a_O_thresh = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
a_T_thresh_new = Array{Float64,3}(undef,length(a_grid),n_g,n_occ-1)
if previous == 1
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_τW=[-0.8 -0.8]_τE=[0.0 0.0]_A=10.0_α=0.1_β=0.3_η=0.1_σ=0.4.jld")
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/previousParameterization.jld")
    d = load("C:/Users/julia/Desktop/teacher/previousParameterization.jld")
    a_T_thresh=d["a_T_thresh"]
    t = d["t"]
    HH_T = d["HH_T"]
    H_O = d["H_O"]
else
    for ia in 1:n_occ-1
        for iG in 1:n_g
            a_T_thresh[:,iG,ia]=1.0.*a_T_thresh_start[:,iG,ia]#a_grid.^0.7
        end
    end
    t[:,1] = zeros(H_grid_length)
    HH_T = H_grid
    H_O = H_grid
end
e_T = Array{Float64,3}(undef,length(a_grid),H_grid_length,n_g)
e_O = Array{Float64,4}(undef,length(a_grid),H_grid_length,n_g,n_occ-1)

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
sum_E_O = Array{Float64,2}(undef,H_grid_length,n_g)

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
#ratio = Array{Float64,1}(undef,n_g)
f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
f1_1=zeros(n_occ-1)
f3_1=zeros(n_occ-1)

# Find a_thresh that does not depend on H for each other occ
iAA=0
conv_a=1000
while conv_a>tol_a && iAA<1000#maxiterAA
    global a_T_thresh
    iAA+=1
    println("Iteration for a_T_thresh=",iAA)

    for iG in 1:n_g
        for i_occ in 1:n_occ-1
            f_2[:,iG] = marginal[:,i_occ,iG].*cdf.(LogNormal(mean_a,std_a),a_T_thresh[:,iG,i_occ]) # Fraction of other given a_O
            f_O[:,iG] = pdf.(LogNormal(mean_a,std_a),a_grid).*f_2[:,iG] # mass of other given a_O

            spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ])
            spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG])

            f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
            mass_O[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold

            f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1] # int(a_O^(α/(1-η))f_O da_O)
            #println(maximum([spl_marg(lowbnd),lowbnd]))
            # inverse threshold
            spl_inv = Spline1D(a_T_thresh[:,iG,i_occ],a_grid, bc="extrapolate")
            for ia in 1:length(a_grid)
                a = a_grid[ia]
                a_O_thresh[ia,iG,i_occ] = spl_inv(a_grid[ia])#maximum([spl_inv(a_grid[ia]),lowbnd])
                # fract of teach occ by occ, given a_T
                f_1_T_tmp[ia,i_occ,iG]=maximum([quadgk(aa ->spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1],0.0])
                #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
            end
        end

        for ia in 1:length(a_grid)
            f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
        end
        spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG])

        f_T[:,iG]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,iG]  # mass of teachers given a_T
        mass_T[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]

        f2[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa)*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]
    end

    #ratio[iG]=sum(f1[:,iG])*f2[iG]/sum( Agrid.^(1/(1-η)).*((ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*f3[:,iG] )

    f1_1=zeros(n_occ-1)
    f3_1=zeros(n_occ-1)
    for i_occ in 1:n_occ-1
        f1_1[i_occ]=sum(f1[i_occ,:].*gm./sum(gm))
        f3_1[i_occ]=sum( ((ones(n_g).-τ_w[i_occ,:])./(ones(n_g).+τ_e[i_occ,:])).^(η/(1-η)).*Agrid[i_occ]^(1/(1-η)).*f3[i_occ,:].*gm./sum(gm))
    end
    ratio=sum(f2.*gm./sum(gm))*sum(f1_1)/sum(f3_1)

    for iG in 1:n_g
        for i_occ in 1:n_occ-1
            a_T_thresh_new[:,iG,i_occ]=( (((1-τ_w[i_occ,iG])^(1/(1-η)))*((1+τ_e[i_occ,iG])^(-η/(1-η)))*cnst*ratio*Agrid[i_occ]^(1/(1-η))).*(a_grid.^(α/(1-η)))).^((σ/β-η)/α)
        end
    end

    abs_err = abs.(a_T_thresh_new[:,:,:] .- a_T_thresh[:,:,:])
    ind_max_err = findall(x->x==maximum(abs_err),abs_err)
    sign_max_err = sign.(a_T_thresh_new[ind_max_err[1][1],ind_max_err[1][2],ind_max_err[1][3]] .- a_T_thresh[ind_max_err[1][1],ind_max_err[1][2],ind_max_err[1][3]])
    conv_a=maximum(abs_err)
    println("error=", conv_a.*sign_max_err[1])
    a_T_thresh[:,:,:]=v.*a_T_thresh_new[:,:,:] + (1-v).*a_T_thresh[:,:,:]
end

filename1 = string("C:/Users/julia/Desktop/teacher/a_T_thresh_year=",year,"_",n_g,"_groups_α=",round(α,digits=2),"_β=",round(β,digits=2),"_η=",round(η,digits=2),"_σ=",round(σ,digits=2),"_ϕ=",round(ϕ,digits=2),"_μ=",round(μ,digits=2),"_λf=",round.(λf,digits=2),"_λm=",round.(λm,digits=2),".jld")
save(filename1,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"α",α,"β",β,"σ",σ,"η",η,"ϕ",ϕ,"μ",μ,"Agrid",Agrid,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"h_T",h_T,"h_O",h_O,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1_T",f_1_T,"N", N,"gm",gm,"av_a",av_a, "year",year, "λm",λm,"λf",λf, "τ_w_opt", τ_w_opt)


###########################
########## PLOTS ##########
###########################
#a_T_thresh_start=1.0.*a_T_thresh

#=
d1 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.5_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=8.55_λm=9.3.jld")
a_T_thresh = d1["a_T_thresh"]
β = d1["β"]
σ = d1["σ"]
λf= d1["λf"]
λm = d1["λm"]
mass_O = d1["mass_O"]
mass_T = d1["mass_T"]
av_a = d1["av_a"]
f_1_T=d1["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d1["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])


plt25 =plot(a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot(a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot(cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot(cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot(a_grid[1:50],(a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot(1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")


d6 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.5_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=8.98_λm=9.3.jld")
a_T_thresh = d6["a_T_thresh"]
β = d6["β"]
σ = d6["σ"]
λf= d6["λf"]
λm = d6["λm"]
mass_O = d6["mass_O"]
mass_T = d6["mass_T"]
av_a = d6["av_a"]
f_1_T=d6["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d6["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])

plt25 =plot!(plt25,a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot!(plt27,a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot!(plt28,cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot!(plt35,cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot!(plt5,(a_grid[1:50],a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")


d7 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.5_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=9.41_λm=9.3.jld")
a_T_thresh = d7["a_T_thresh"]
β = d7["β"]
σ = d7["σ"]
λf= d7["λf"]
λm = d7["λm"]
mass_O = d7["mass_O"]
mass_T = d7["mass_T"]
av_a = d7["av_a"]
f_1_T=d7["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d7["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])

plt25 =plot!(plt25,a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot!(plt27,a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot!(plt28,cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot!(plt35,cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot!(plt5,(a_grid[1:50],a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,2],grid=:none,label=latexstring("men"),legendfontsize=5, legend=:topleft,xticks = ([1:1:20;], occ[1:end-1]),xrotation = 45,xtickfontsize=4,ylabel="Labor Market Barriers")


#savefig(plt25,"C:/Users/julia/Desktop/teacher/a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
#savefig(plt35,"C:/Users/julia/Desktop/teacher/rank_a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
savefig(plt25,"C:/Users/julia/Desktop/teacher/a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt27,"C:/Users/julia/Desktop/teacher/a_T_frac_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt28,"C:/Users/julia/Desktop/teacher/rank_a_T_frac_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt35,"C:/Users/julia/Desktop/teacher/rank_a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")

#savefig(plt5,"C:/Users/julia/Desktop/teacher/a_T_thresh_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
#savefig(plt55,"C:/Users/julia/Desktop/teacher/τ_w_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
savefig(plt5,"C:/Users/julia/Desktop/teacher/a_T_thresh_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt55,"C:/Users/julia/Desktop/teacher/τ_w_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")








d1 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.4_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=10.95_λm=11.75.jld")
a_T_thresh = d1["a_T_thresh"]
β = d1["β"]
σ = d1["σ"]
λf= d1["λf"]
λm = d1["λm"]
mass_O = d1["mass_O"]
mass_T = d1["mass_T"]
av_a = d1["av_a"]
f_1_T=d1["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d1["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])


plt25 =plot(a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot(a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot(cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot(cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot(a_grid[1:50],(a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot(1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")


d6 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.4_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=11.5_λm=11.75.jld")
a_T_thresh = d6["a_T_thresh"]
β = d6["β"]
σ = d6["σ"]
λf= d6["λf"]
λm = d6["λm"]
mass_O = d6["mass_O"]
mass_T = d6["mass_T"]
av_a = d6["av_a"]
f_1_T=d6["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d6["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])

plt25 =plot!(plt25,a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot!(plt27,a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot!(plt28,cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot!(plt35,cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot!(plt5,(a_grid[1:50],a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")


d7 = load("C:/Users/julia/Desktop/teacher/a_T_thresh_year=2010_2_groups_α=0.09_β=0.4_η=0.1_σ=0.4_ϕ=1.0_μ=0.5_λf=12.04_λm=11.75.jld")
a_T_thresh = d7["a_T_thresh"]
β = d7["β"]
σ = d7["σ"]
λf= d7["λf"]
λm = d7["λm"]
mass_O = d7["mass_O"]
mass_T = d7["mass_T"]
av_a = d7["av_a"]
f_1_T=d7["f_1_T"]
spl_f_1_T=Spline1D(a_grid, f_1_T[:,1], bc="extrapolate")
f_T[:,1]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,1]
τ_w=d7["τ_w"]
τ_w2=ones(n_occ-1,n_g).-(ones(n_occ-1,n_g).-τ_w)./(ones(n_occ-1,n_g).-hcat(τ_w[:,2],τ_w[:,2]))

# average ability in other
av_a=Array{Float64,3}(undef,n_occ,n_g, H_grid_length)
 H_grid_length=1
for iG in 1:n_g
    for i_occ in 1:n_occ-1
        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
        spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")
        for iH in 1:H_grid_length
            av_a[i_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
        end
    end
end

# average ability in teaching
for iG in 1:n_g
    spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")
    for iH in 1:H_grid_length
        av_a[n_occ,iG,iH]=quadgk(aa -> cdf(LogNormal(mean_a,std_a),aa)*pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]
    end
end

println(mass_T)
println(av_a[end,:,1])

plt25 =plot!(plt25,a_grid, f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topright,xlabel="a_T",ylabel="f_T")
plt27 =plot!(plt27,a_grid, f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:bottomright,xlabel="a_T",ylabel="Fraction of Teachers")
plt28 =plot!(plt28,cdf.(LogNormal.(mean_a,std_a),a_grid), f_1_T[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="Fraction of Teachers")

plt35 =plot!(plt35,cdf.(LogNormal.(mean_a,std_a),a_grid), f_T[:,1]./mass_T[1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="rank_a_T",ylabel="f_T")

plt5 = plot!(plt5,(a_grid[1:50],a_T_thresh[1:50,1,1]),grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="a_occ1",ylabel="a_T_thresh")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,1],grid=:none,label=latexstring("\\frac{\\beta}{\\sigma}_f = "*string(round(β/σ,digits=2))*", \\lambda_f = "*string(round(λf,digits=2))),legend=:topleft,xlabel="Other Occupations",ylabel="Labor Market Barriers")
plt55 = plot!(plt55,1:n_occ-1,τ_w2[:,2],grid=:none,label=latexstring("men"),legendfontsize=5, legend=:topleft,xticks = ([1:1:20;], occ[1:end-1]),xrotation = 45,xtickfontsize=4,ylabel="Labor Market Barriers")


#savefig(plt25,"C:/Users/julia/Desktop/teacher/a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
#savefig(plt35,"C:/Users/julia/Desktop/teacher/rank_a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
savefig(plt25,"C:/Users/julia/Desktop/teacher/a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt27,"C:/Users/julia/Desktop/teacher/a_T_frac_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt28,"C:/Users/julia/Desktop/teacher/rank_a_T_frac_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt35,"C:/Users/julia/Desktop/teacher/rank_a_T_f_T_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")

#savefig(plt5,"C:/Users/julia/Desktop/teacher/a_T_thresh_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
#savefig(plt55,"C:/Users/julia/Desktop/teacher/τ_w_beta_sigma= "*string(round(β/σ,digits=2))*".eps")
savefig(plt5,"C:/Users/julia/Desktop/teacher/a_T_thresh_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")
savefig(plt55,"C:/Users/julia/Desktop/teacher/τ_w_beta_sigma= "*string(round(β/σ,digits=2))*"_2010_v2.pdf")



=#
