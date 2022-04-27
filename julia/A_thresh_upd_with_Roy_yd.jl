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
occ =["other1", "other2", "other3","other4", "other5", "other6","other7", "other8", "other9","other10","other11", "other12", "other13","other14", "other15", "other16","other17", "other18", "other19", "home production", "teaching"]
ν = 0.1#0.2
ν2 = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations

τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')

τ_w[:,1]=fill!(τ_w[:,1],-3)
#τ_w[1,1] = -3
τ_w[:,2]=fill!(τ_w[:,2],-4)
#τ_w[1,2] = -4
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_e=fill!(τ_e,0.0)
#τ_e[1,1] = 0.0
#τ_e[1,2] = τ_e[1,1]
gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)
α=.1
η=.1
β=.5
σ=.4
μ=1/2
ϕ=1/3
θ=3

# productivity in 'Other' occupations
Agrid=Array{Float64,1}(undef,n_occ-1)
#Agrid=fill!(Agrid, 1)
Agrid=collect(range(1,1.1,length=n_occ-1))
#Agrid=collect(range(2,2.2,length=n_occ-1))
#Agrid[1]=0.98
#Agrid[end]=0.98

quantile_top = .999

# Set up the grid for aggregate human capital in teaching:
H_grid_length = 21
# Set the percentage range above / below fixed point:
H_grid_range = 0.2
H_grid = collect(range((1-H_grid_range)*1,stop=(1+H_grid_range)*1,length=H_grid_length))
# a_grid=quantile.(Frechet(θ),.005:.015:1)
mean_a=0
std_a=1

a_grid=exp.(append!(quantile.(Normal(mean_a,std_a),.001:.001:.004), append!(quantile.(Normal(mean_a,std_a),.005:.015:.979), quantile.(Normal(mean_a,std_a),.980:.002:quantile_top))))

# Integral bounds
lowbnd=0.0
upbnd=quantile.(LogNormal(mean_a,std_a),quantile_top) # set to largest element in 'a_grid'

# given a productivity vector Agrid
# this function calculates share of people choosing occupaton occ_id
marg=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
function share(occ_id,iG,Agrid)
    for ia in 1:length(a_grid)
        A_tmp=(( ((ones(n_occ-2).-τ_w[occ_id,iG])./(ones(n_occ-2).-τ_w[1:end.!=occ_id,iG])).*(((ones(n_occ-2).+τ_e[occ_id,iG])./(ones(n_occ-2).+τ_e[1:end.!=occ_id,iG])).^(-η)).*(Agrid[occ_id]./Agrid[1:end .!=occ_id])).^(1/α))
        marg[ia,occ_id,iG]=prod(cdf.(LogNormal(mean_a,std_a),A_tmp.*a_grid[ia]))
    end
    spl_marg=Spline1D(a_grid, marg[:,occ_id,iG],bc="extrapolate")
    return quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa),lowbnd,upbnd)[1], marg[:,occ_id,iG]
end

# share of individuals employed in occ i among all non-teaching individuals among this group
# then all share_occ will add up to 1*n_g
share_occ = Array{Float64,2}(undef,n_occ-1,n_g)
marginal=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)
for iG in 1:n_g
    for i_occ = 1:(n_occ-1)
        share_occ[i_occ,iG]=share(i_occ,iG,Agrid)[1]
        marginal[:,i_occ,iG]=share(i_occ,iG,Agrid)[2]
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
ratio = Array{Float64,1}(undef,n_g)
f_1_T_tmp=Array{Float64,3}(undef,length(a_grid),n_occ-1,n_g)

# Find a_thresh that does not depend on H for each other occ
for iG in 1:n_g
    iAA=0
    conv_a=1000
    while conv_a>tol_a && iAA<500#maxiterAA
        global a_T_thresh
        iAA+=1
        println("Iteration for a_T_thresh=",iAA)

        for i_occ in 1:n_occ-1
            f_2[:,iG] = marginal[:,i_occ,iG].*cdf.(LogNormal(mean_a,std_a),a_T_thresh[:,iG,i_occ]) # Fraction of other given a_O
            f_O[:,iG] = pdf.(LogNormal(mean_a,std_a),a_grid).*f_2[:,iG] # mass of other given a_O

            spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG,i_occ],bc="extrapolate")
            spl_marg=Spline1D(a_grid, marginal[:,i_occ,iG],bc="extrapolate")

            f1[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
            mass_O[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold

            f3[i_occ,iG]=quadgk(aa -> spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1] # int(a_O^(α/(1-η))f_O da_O)
            #println(maximum([spl_marg(lowbnd),lowbnd]))
            # inverse threshold
            spl_inv = Spline1D(a_T_thresh[:,iG,i_occ],a_grid,bc="extrapolate")
            for ia in 1:length(a_grid)
                a = a_grid[ia]
                a_O_thresh[ia,iG,i_occ] = maximum([spl_inv(a_grid[ia]),lowbnd])
                # fract of teach occ by occ, given a_T
                f_1_T_tmp[ia,i_occ,iG]=quadgk(aa ->spl_marg(aa)*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
                #quadgk(aa ->maximum([spl_marg(aa),lowbnd])*pdf(LogNormal(mean_a,std_a),aa),lowbnd,a_O_thresh[ia,iG,i_occ])[1]
            end
        end

        for ia in 1:length(a_grid)
            f_1_T[ia,iG]= sum(f_1_T_tmp[ia,:,iG]) # Fraction of teachers given a_T
        end
        spl_f_1_T=Spline1D(a_grid, f_1_T[:,iG], bc="extrapolate")

        f_T[:,iG]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,iG]  # mass of teachers given a_T
        mass_T[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa),lowbnd,upbnd)[1]

        f2[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*spl_f_1_T(aa)*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]

        ratio[iG]=sum(f1[:,iG])*f2[iG]/sum( Agrid.^(1/(1-η)).*((ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*f3[:,iG] )

        for i_occ in 1:n_occ-1
            a_T_thresh_new[:,iG,i_occ]=( (((1-τ_w[i_occ,iG])^(1/(1-η)))*((1+τ_e[i_occ,iG])^(-η/(1-η)))*cnst*ratio[iG]*Agrid[i_occ]^(1/(1-η))).*(a_grid.^(α/(1-η)))).^((σ/β-η)/α)
        end

        abs_err = abs.(a_T_thresh_new[:,iG,:] .- a_T_thresh[:,iG,:])
        ind_max_err = findall(x->x==maximum(abs_err),abs_err)
        sign_max_err = sign.(a_T_thresh_new[ind_max_err[1][1],iG,ind_max_err[1][2]] .- a_T_thresh[ind_max_err[1][1],iG,ind_max_err[1][2]])
        conv_a=maximum(abs_err)
        println("error=", conv_a.*sign_max_err[1])
        a_T_thresh[:,iG,:]=ν.*a_T_thresh_new[:,iG,:] + (1-ν).*a_T_thresh[:,iG,:]
    end
end

# Initial guess for fixed point of HH_T:
# Set the array index to the midpoint of H_grid:
iH = convert(Int,ceil(H_grid_length/2))
if previous == 0
    HH_fp = 1
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
                e_T[ia,iH,iG]=(1-t[iH,1])*β/σ*η^(1/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])*a^(α/(σ/β-η))/f2[iG]
            end
            # Compute HHH_T:
            h_T[:,iH,iG] = (2*HH_fp/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
            HHH_T_0[iG]= ( (1-t[iH,1])^η*η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iH,iG] = η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])

            for i_occ in 1:n_occ-1
                for ia in 1:length(a_grid)
                    a=a_grid[ia]
                    e_O[ia,iH,iG,i_occ]=((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG])*Agrid[i_occ]*η*s_O^ϕ*(2*HH_fp/M)^σ*a^α)^(1/(1-η))
                end
                # Compute H_O:
                h_O[:,iH,iG,i_occ] = (2*HH_fp/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG,i_occ].^η
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG] = (((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ])^(1/(1-η))*f3[i_occ,iG]
            end
        end
        HH_T[iH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
        end

        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        convT = abs(t[iH,2]-t[iH,1])
         println("convT=",round(convT,digits=3))
         println("t[iH,:]=",round.(t[iH,:],digits=3))
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[iH,1] = min(t[iH,2],1-1e-3) # ν*t[iH,1]+(1-ν)*t[iH,2]

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
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        println("    T: iteration ",iterT)
        println("HH_T[iH]=",round(HH_T[iH],digits=3))
        println("t[iH,:]=",round.(t[iH,:],digits=3))
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                e_T[ia,iH,iG]=(1-t[iH,1])*β/σ*η^(1/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])*a^(α/(σ/β-η))/f2[iG]
            end
            # Compute HHH_T:
            h_T[:,iH,iG] = (2*HH_fp/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
            HHH_T_0[iG]= ( (1-t[iH,1])^η*η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
            # Total earnings of teachers in each group:
            E_T[iH,iG] = η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])

            for i_occ in 1:n_occ-1
                for ia in 1:length(a_grid)
                    a=a_grid[ia]
                    e_O[ia,iH,iG,i_occ]=((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG])*Agrid[i_occ]*η*s_O^ϕ*(2*HH_fp/M)^σ*a^α)^(1/(1-η))
                end
                # Compute H_O:
                h_O[:,iH,iG,i_occ] = (2*HH_fp/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG,i_occ].^η
                H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
                # Total earnings of others in each group:
                E_O[i_occ,iH,iG] = (((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ])^(1/(1-η))*f3[i_occ,iG]
            end
        end
        HH_T[iH]=sum(HHH_T_0.*gm)
        for i_occ in 1:n_occ-1
            H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
        end
        for iG in 1:n_g
            sum_E_O[iH,iG]=sum(E_O[:,iH,iG])
        end

        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(sum_E_O[iH,:].*gm))
        convT = abs(t[iH,2]-t[iH,1])
         println("convT=",round(convT,digits=3))
         println("t[iH,:]=",round.(t[iH,:],digits=3))
        # println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[iH,1] = min(t[iH,2],1-1e-3) # ν*t[iH,1]+(1-ν)*t[iH,2]

        iterT = iterT+1
    end
end

for iH in 1:H_grid_length
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Initiate 'while' loop ('not indexed') over the tax rate:
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            a=a_grid[ia]
            e_T[ia,iH,iG]=(1-t[iH,1])*β/σ*η^(1/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])*a^(α/(σ/β-η))/f2[iG]
        end
        # Compute HHH_T:
        h_T[:,iH,iG] = (2*HH_fp/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
        HHH_T_0[iG]= ( (1-t[iH,1])^η*η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*(f2[iG])^(σ/β-η) )^(β/σ)
        # Total earnings of teachers in each group:
        E_T[iH,iG] = η^(η/(1-η))*(2*HH_fp/M)^(σ/(1-η))*sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG])

        for i_occ in 1:n_occ-1
            for ia in 1:length(a_grid)
                a=a_grid[ia]
                e_O[ia,iH,iG,i_occ]=((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG])*Agrid[i_occ]*η*s_O^ϕ*(2*HH_fp/M)^σ*a^α)^(1/(1-η))
            end
            # Compute H_O:
            h_O[:,iH,iG,i_occ] = (2*HH_fp/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG,i_occ].^η
            H_O_0[i_occ,iG]=(((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ]^η)^(1/(1-η))*f3[i_occ,iG]
            # Total earnings of others in each group:
            E_O[i_occ,iH,iG] = (((1-t[iH,1])*(1-τ_w[i_occ,iG])/(1+τ_e[i_occ,iG]))^η*η^η*(2*HH_fp/M)^σ*s_O.^ϕ*Agrid[i_occ])^(1/(1-η))*f3[i_occ,iG]
        end
    end
    HH_T[iH]=sum(HHH_T_0.*gm)
    for i_occ in 1:n_occ-1
        H_O[i_occ,iH]=sum(H_O_0[i_occ,:].*gm)
    end
end

#Check:
# t_check is correct for ng=2
#t_check=Array{Float64,1}(undef,n_g)
HH_T_0_check=Array{Float64,2}(undef,H_grid_length,n_g)
HH_T_check=Array{Float64,1}(undef,H_grid_length)
for iG in 1:n_g
    #t_check[iG]=1/(1+sum(f1[:,iG]))
    for iH in 1:H_grid_length
        HH_T_0_check[iH,iG]=((1-t[iH,1])^η*η^(η/(1-η))*(2/M)^(σ/(1-η))*s_T^ϕ*(β/σ)^η*(sum( Agrid.^(1/(1-η)).*((1-t[iH,1]).*(ones(n_occ-1).-τ_w[:,iG])./(ones(n_occ-1).+τ_e[:,iG])).^(η/(1-η)).*s_O^(ϕ/(1-η)).*f3[:,iG] )/sum(f1[:,iG]))^η*f2[iG]^(σ/β-η))^(β/σ)
    end
end
for iH in 1:H_grid_length
    HH_T_check[iH]=sum(HH_T_0_check[iH,:].*gm)^((1-η)/(1-η-β))
end

#### Group-specific barriers
# HH_T_cf=zeros(H_grid_length)
# H_O_cf=zeros(H_grid_length)
# for iH in 1:H_grid_length
#     H=H_grid[iH]
#     for iG in 1:n_g
#         HH_T_cf[iH]+=(gm[iG]/sum(gm))*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[iH,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[iG])/(1+τ_e[iG]))^(η*η*β/σ/(1-η))*(2*H/M)^(β/(1-η))*f2[iG]^(1-β*η/σ)*(f3[iG]/f1[iG])^(η*β/σ)
#         H_O_cf[iH]+=(gm[iG]/sum(gm))*(A*η*(1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*s_O^(ϕ/(1-η))*(2*H/M)^(σ/(1-η))*f3[iG]
#     end
# end

N=zeros(length(a_grid),H_grid_length,n_g)
for iH in 1:H_grid_length
    H=HH_T[iH] # H_grid[iH]
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            N[ia,iH,iG]=(M/2)/H*(h_T[ia,iH,iG])^(β/σ)
        end
    end
end
###################################
##### CHECK THIS
###################################
#=
# Compute mean and c.v. associated with class size distribution (calculate first and second moment
for iG in 1:n_g
    global spl_pdf_T
    # Splines for 'N' and pdf of teachers:
    spl_EN = Spline1D(a_grid,N[:,1,iG])
    spl_pdf_T = Spline1D(a_grid,(cdf.(LogNormal(mean_a,std_a),a_O_thresh[:,iG])).*pdf.(LogNormal(mean_a,std_a),a_grid))
    # First moment:
    EN[iG] = quadgk(aa -> spl_pdf_T(aa) * spl_EN(aa),lowbnd,upbnd)[1]
    # EN[iG] = quadgk(aa -> spl_pdf_T(aa),lowbnd,upbnd)[1]
    # Second moment:
    spl_EN2 = Spline1D(a_grid,N[:,1,iG].^2)
    EN2[iG] = quadgk(aa -> spl_pdf_T(aa) * spl_EN2(aa),lowbnd,upbnd)[1]
end
EN_agg = (gm./(sum(gm))*EN)[1]
EN2_agg = (gm./(sum(gm))*EN2)[1]
cvN = (EN2_agg - EN_agg.^2).^(.5)./EN_agg
#/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/
filename1 = string("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_",n_g,"_groups_τW=",τ_w,"_τE=",τ_e,"_A=",round(A,digits=2),"_α=",round(α,digits=2),"_β=",round(β,digits=2),"_η=",round(η,digits=2),"_σ=",round(σ,digits=2),".jld")
save(filename1,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"α",α,"β",β,"σ",σ,"η",η,"A",A,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2,"N", N,"EN_agg",EN_agg,"cvN",cvN,"gm",gm)
filename2 = string("//Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/previousParameterization.jld")
save(filename2,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"α",α,"β",β,"σ",σ,"η",η,"A",A,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2,"N", N,"EN_agg",EN_agg,"cvN",cvN,"gm",gm)

# Find steady state
# using Roots
# f(x)=x-gm[1]/sum(gm)*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[1,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[1])/(1+τ_e[1]))^(η*η*β/σ/(1-η))*(2/M)^(β/(1-η))*f2[1]^(1-β*η/σ)*(f3[1]/f1[1])^(η*β/σ)*x^(β/(1-η))-gm[2]/sum(gm)*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[1,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[2])/(1+τ_e[2]))^(η*η*β/σ/(1-η))*(2/M)^(β/(1-η))*f2[2]^(1-β*η/σ)*(f3[2]/f1[2])^(η*β/σ)*x^(β/(1-η))
# find_zero(f,2)


###########################################################
##################  PARAMETERIZATION  #####################
###########################################################
#=
println("__________")
println("EN_agg= ",EN_agg)
println("cvN= ",cvN)
println("mass_T_female= ",mass_T[1])
println("mass_T_male= ",mass_T[2])
println("mass_T_female/mass_T_male= ",mass_T[1]*gm[1]/sum(mass_T.*vec(gm)))
av_e_O=zeros(n_g)
av_e_T=zeros(n_g)
for iG in 1:n_g
    spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")
    spl_T_thresh=Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")
    spl_e_O = Spline1D(a_grid,e_O[:,1,iG])
    spl_e_T = Spline1D(a_grid,e_T[:,1,iG])
    av_e_O[iG]=quadgk(aa -> spl_e_O(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
    av_e_T[iG]=quadgk(aa -> spl_e_T(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]
end

println("average e_T_female= ",av_e_T[1])
println("average e_T_male= ",av_e_T[2])
println("average e_O_female= ",av_e_O[1])
println("average e_O_male= ",av_e_O[2])
println("s_T/s_O= ",s_T/s_O)

println("average w_T_female= ",E_T[1,1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1])
println("average w_T_male= ",E_T[1,2]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1])
println("average w_O_female= ",E_O[1,1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1])
println("average w_O_male= ",E_O[1,2]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1])
=#

println("__________")
println(EN_agg)
println(cvN)
println(mass_T[1])
println(mass_T[2])
println(mass_T[1]*gm[1]/sum(mass_T.*vec(gm)))
av_e_O=zeros(n_g)
av_e_T=zeros(n_g)
for iG in 1:n_g
    spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")
    spl_T_thresh=Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")
    spl_e_O = Spline1D(a_grid,e_O[:,1,iG])
    spl_e_T = Spline1D(a_grid,e_T[:,1,iG])
    av_e_O[iG]=quadgk(aa -> spl_e_O(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
    av_e_T[iG]=quadgk(aa -> spl_e_T(aa)*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]
    mass_O[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
    mass_T[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]
end
println(av_e_T[1])
println(av_e_T[2])
println(av_e_O[1])
println(av_e_O[2])
println(s_T/s_O)
println(E_T[1,1]/mass_T[1])
println(E_T[1,2]/mass_T[2])
println(E_O[1,1]/mass_O[1])
println(E_O[1,2]/mass_O[2])
av_a_O=zeros(n_g)
av_a_T=zeros(n_g)
for iG in 1:n_g
    spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")
    spl_T_thresh=Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")
    av_a_O[iG]=quadgk(aa -> aa*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1]
    av_a_T[iG]=quadgk(aa -> aa*pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]/quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1]
end
println(av_a_T[1])
println(av_a_T[2])
println(av_a_O[1])
println(av_a_O[2])

E_O_1=zeros(n_g)
E_O_2=zeros(n_g)
for iG in 1:n_g
    global spl_pdf_O
    # Splines for 'N' and pdf of teachers:
    spl_E_O = Spline1D(a_grid,(1-τ_w[iG])*A*h_O[:,1,iG])
    spl_pdf_O = Spline1D(a_grid,(cdf.(LogNormal(mean_a,std_a),a_T_thresh[:,iG])).*pdf.(LogNormal(mean_a,std_a),a_grid))
    # First moment:
    E_O_1[iG] = quadgk(aa -> spl_pdf_O(aa) * spl_E_O(aa),lowbnd,upbnd)[1]
    # Second moment:
    spl_E_O_2 = Spline1D(a_grid,((1-τ_w[iG])*A*h_O[:,1,iG]).^2)
    E_O_2[iG] = quadgk(aa -> spl_pdf_O(aa) * spl_E_O_2(aa),lowbnd,upbnd)[1]
end
E_O_1_agg = (gm./(sum(gm))*E_O_1)[1]
E_O_2_agg = (gm./(sum(gm))*E_O_2)[1]
cv_E_O = (E_O_2_agg - E_O_1_agg.^2).^(.5)./E_O_1_agg

println(cv_E_O)

E_T_1=zeros(n_g)
E_T_2=zeros(n_g)
for iG in 1:n_g
    global spl_pdf_T
    # Splines for 'N' and pdf of teachers:
    spl_E_T = Spline1D(a_grid,H_O[1]*A/HH_T[1]/f1[iG]*(2*H_grid[1]/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[1,1])*A*η*s_O^ϕ*(2*H_grid[1]/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ).*a_grid.^(α/(σ/β-η)))
    spl_pdf_T = Spline1D(a_grid,(cdf.(LogNormal(mean_a,std_a),a_O_thresh[:,iG])).*pdf.(LogNormal(mean_a,std_a),a_grid))
    # First moment:
    E_T_1[iG] = quadgk(aa -> spl_pdf_T(aa) * spl_E_T(aa),lowbnd,upbnd)[1]
    # Second moment:
    spl_E_T_2 = Spline1D(a_grid,(H_O[1]*A/HH_T[1]/f1[iG]*(2*H_grid[1]/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[1,1])*A*η*s_O^ϕ*(2*H_grid[1]/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ).*a_grid.^(α/(σ/β-η))).^2)
    E_T_2[iG] = quadgk(aa -> spl_pdf_T(aa) * spl_E_T_2(aa),lowbnd,upbnd)[1]
end
E_T_1_agg = (gm./(sum(gm))*E_T_1)[1]
E_T_2_agg = (gm./(sum(gm))*E_T_2)[1]
cv_E_T = (E_T_2_agg - E_T_1_agg.^2).^(.5)./E_T_1_agg

println(cv_E_T)

# Create array with all moments of interest:
# (1) Vector with 11 parameters:
paramExp=Vector{Float64}(undef,10)
paramExp=[α,η,β,σ,μ,ϕ,A,θ,τ_w[1,1],τ_w[1,2]]
# (2) Vector with XX moments:
momExp=Vector{Float64}(undef,21)

# Create DataFrame with parameters:
calibration_params = DataFrame(Label = ["α","η","β","σ","μ","ϕ","A","θ","τ_w (female)","τ_w (male)"],Value = [α,η,β,σ,μ,ϕ,A,θ,τ_w[1,1],τ_w[1,2]])

# Create DataFrame with moments:
calibration_moms = DataFrame(Moment = ["average class size","cv of class size","share of teachers among female","share of teachers among male","share of female vs male teachers","average good investment in T female","average good investment in T male","average good investment in O female","average good investment in O male","s_T/s_O","average wage in T female","average wage in T male","average wage in O female","average wage in O male","average ability in T female","average ability in T male","average ability in O female","average ability in O male","wage dispersion in T","wage dispersion in O"],Value=[EN_agg,cvN,mass_T[1],mass_T[2],mass_T[1]*gm[1]/sum(mass_T.*vec(gm)),av_e_T[1],av_e_T[2],av_e_O[1],av_e_O[2],s_T/s_O,E_T[1,1]/mass_T[1],E_T[1,2]/mass_T[2],E_O[1,1]/mass_O[1],E_O[1,2]/mass_O[2],av_a_T[1],av_a_T[2],av_a_O[1],av_a_O[2],cv_E_T,cv_E_O])

# Create DataFrame with parameters + moments (for export / write to file):
calibration_df = DataFrame(Label = ["α","η","β","σ","μ","ϕ","A","θ","τ_w (female)","τ_w (male)","average class size","cv of class size","share of teachers among female","share of teachers among male","share of female vs male teachers","average good investment in T female","average good investment in T male","average good investment in O female","average good investment in O male","s_T/s_O","average wage in T female","average wage in T male","average wage in O female","average wage in O male","average ability in T female","average ability in T male","average ability in O female","average ability in O male","wage dispersion in T","wage dispersion in O"],Value = [α,η,β,σ,μ,ϕ,A,θ,τ_w[1,1],τ_w[1,2],EN_agg,cvN,mass_T[1],mass_T[2],mass_T[1]*gm[1]/sum(mass_T.*vec(gm)),av_e_T[1],av_e_T[2],av_e_O[1],av_e_O[2],s_T/s_O,E_T[1,1]/mass_T[1],E_T[1,2]/mass_T[2],E_O[1,1]/mass_O[1],E_O[1,2]/mass_O[2],av_a_T[1],av_a_T[2],av_a_O[1],av_a_O[2],cv_E_T,cv_E_O])
# CSV.write(fullpath_results,calibration_df)
calibration_csv = permutedims(Matrix(calibration_df))

fullpath_results = pwd()*"/parameterizations/calibration.csv"
# For new CSV file:
# writedlm(fullpath_results,  calibration_csv, ',')
# To append to existing CSV file (ignoring header row):
io = open(fullpath_results,"a")
writedlm(io,  calibration_csv[2:end,:], ','; header=true)
close(io)
=#
