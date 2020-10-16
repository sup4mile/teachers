using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings

# Load occupational threshold from previous parameterization:
previous = 1;
# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
occ = ["teaching","other"]
ν = 0.2

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations
τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[1,1] = -3
τ_w[1,2] = -4
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_e[1,1] = 0.0
τ_e[1,2] = τ_e[1,1]
gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)
α=.1
η=.1
β=.3
σ=.4
μ=1/2
ϕ=1/3
A=10 # productivity in 'Other'
θ=3

quantile_top = .999
#H_grid=collect(range(.05,stop=1.4,length=17))
H_grid=collect(range(.0005,stop=1.0,length=20))
# a_grid=quantile.(Frechet(θ),.005:.015:1)
mean_a=0
std_a=1
a_grid=exp.(append!(quantile.(Normal(mean_a,std_a),.001:.001:.004), append!(quantile.(Normal(mean_a,std_a),.005:.015:.979), quantile.(Normal(mean_a,std_a),.980:.002:quantile_top))))
s_T=μ*ϕ/(μ*ϕ+σ/β-η)
s_O=μ*ϕ/(μ*ϕ+1-η)
cnst=(1-η)/(1-η*β/σ)*((1-s_O)/(1-s_T))^(1/μ)

maxiterAA=1e5
tol_a=1e-6

tolT= 2.5e-4
maxiterT = 100

# Integral bounds
lowbnd=0.0
upbnd=quantile.(LogNormal(mean_a,std_a),quantile_top) # set to largest element in 'a_grid'

t=zeros(length(H_grid),2)

a_T_thresh = Array{Float64,2}(undef,length(a_grid),n_g)
a_O_thresh = Array{Float64,2}(undef,length(a_grid),n_g)
a_T_thresh_new = Array{Float64,2}(undef,length(a_grid),n_g)
if previous == 1
    # d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_τW=[-0.8 -0.8]_τE=[0.0 0.0]_A=10.0_α=0.1_β=0.3_η=0.1_σ=0.4.jld")
    d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/previousParameterization.jld")
    a_T_thresh=d["a_T_thresh"]
    t = d["t"]
    HH_T = d["HH_T"]
    H_O = d["H_O"]
else
    for iG in 1:n_g
        a_T_thresh[:,iG]=a_grid.^0.7
    end
    t[1,1] = 0
    HH_T = H_grid
    H_O = H_grid
end

f_O = Array{Float64,2}(undef,length(a_grid),n_g)
f_T = Array{Float64,2}(undef,length(a_grid),n_g)

f1 = Array{Float64,1}(undef,n_g)
f1_inv = Array{Float64,1}(undef,n_g)
f2 = Array{Float64,1}(undef,n_g)
f3 = Array{Float64,1}(undef,n_g)
f4 = Array{Float64,1}(undef,n_g)
f4_inv = Array{Float64,1}(undef,n_g)

e_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
e_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

h_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
h_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
ω = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

spl_T = Array{Spline1D,1}(undef,n_g)
spl_O = Array{Spline1D,1}(undef,n_g)
spl_T_inv = Array{Spline1D,1}(undef,n_g)
spl_O_inv = Array{Spline1D,1}(undef,n_g)
spl_ω = Array{Spline1D,1}(undef,n_g)
spl_w = Array{Spline1D,1}(undef,n_g)
spl_h_O = Array{Spline1D,1}(undef,n_g)
spl_pdf_T = Array{Spline1D,1}(undef,n_g)

wb_T = Array{Float64,2}(undef,length(H_grid),n_g)
wb_O = Array{Float64,2}(undef,length(H_grid),n_g)
E_T = Array{Float64,2}(undef,length(H_grid),n_g)
E_O = Array{Float64,2}(undef,length(H_grid),n_g)

HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,1}(undef,n_g)

EN = Array{Float64,1}(undef,n_g)
EN2 = Array{Float64,1}(undef,n_g)
mass_T = Array{Float64,1}(undef,n_g)
mass_O = Array{Float64,1}(undef,n_g)
f_1 = Array{Float64,2}(undef,length(a_grid),n_g)
f_2 = Array{Float64,2}(undef,length(a_grid),n_g)
f_1_T = Array{Float64,2}(undef,length(a_grid),n_g)
f_2_T = Array{Float64,2}(undef,length(a_grid),n_g)

# Find a_thresh that does not depend on H
for iG in 1:n_g
    iAA=0
    conv_a=1000
    while conv_a>tol_a && iAA<maxiterAA
        global a_T_thresh
        iAA+=1
        println("Iteration for a_T_thresh=",iAA)

        f_1[:,iG] = ones(length(a_grid)).-cdf.(LogNormal(mean_a,std_a),a_T_thresh[:,iG]) # Fraction of teachers given a_O
        f_2[:,iG] = cdf.(LogNormal(mean_a,std_a),a_T_thresh[:,iG]) # Fraction of other given a_O
        f_O[:,iG] = pdf.(LogNormal(mean_a,std_a),a_grid).*f_2[:,iG] # mass of other given a_O

        spl_T_thresh = Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")

        f1[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
        mass_O[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa)),lowbnd,upbnd)[1] # total mass of other using direct threshold
        # test1 = quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),.5),lowbnd,upbnd)[1]
        f4[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*(1-cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))),lowbnd,upbnd)[1] # total mass of teachers using direct threshold
        mass_T[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*(1-cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))),lowbnd,upbnd)[1] # total mass of teachers using direct threshold
        # test2 = quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(1-LogNormal(mean_a,std_a),.5),lowbnd,upbnd)[1]
        f3[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^(α/(1-η)),lowbnd,upbnd)[1]

        # println(test1)
        # println(test2)

        spl_inv = Spline1D(a_T_thresh[:,iG],a_grid,bc="extrapolate")
        for ia in 1:length(a_grid)
            a = a_grid[ia]
            a_O_thresh[ia,iG] = maximum([spl_inv(a_grid[ia]),0.0])
        end
        f_1_T[:,iG]= cdf.(LogNormal(mean_a,std_a),a_O_thresh[:,iG]) # Fraction of teachers given a_T
        f_2_T[:,iG]= ones(length(a_grid)).-cdf.(LogNormal(mean_a,std_a),a_O_thresh[:,iG]) # Fraction of other given a_T
        f_T[:,iG]=pdf.(LogNormal(mean_a,std_a),a_grid).*f_1_T[:,iG]  # mass of teachers given a_T

        spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")

        f1_inv[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*(1-cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))),lowbnd,upbnd)[1] # total mass of other using direct threshold
        f4_inv[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa)),lowbnd,upbnd)[1] # total mass of teachers using inverse threshold

        f2[iG]=quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]

        ratio=f1[iG]*f2[iG]/f3[iG]

        a_T_thresh_new[:,iG]=((1-τ_w[iG])*cnst*ratio.*(a_grid.^(α/(1-η)))).^((σ/β-η)/α)
        abs_err = abs.(a_T_thresh_new[:,iG] .- a_T_thresh[:,iG])
        ind_max_err = findall(x->x==maximum(abs_err),abs_err)
        sign_max_err = sign.(a_T_thresh_new[ind_max_err,iG] .- a_T_thresh[ind_max_err,iG])
        conv_a=maximum(abs_err)
        println("error=", conv_a.*sign_max_err[1])
        a_T_thresh[:,iG]=ν.*a_T_thresh_new[:,iG] + (1-ν).*a_T_thresh[:,iG]
    end
end

for iH in 1:length(H_grid)
    println("H: gridpoint ", iH," of ", length(H_grid)," (H = ", H_grid[iH],")")
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
                e_T[ia,iH,iG]=((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*a^(α/(σ/β-η))*f3[iG]/f1[iG]/f2[iG]
                e_O[ia,iH,iG]=((1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG])*A*η*s_O^ϕ*(2*H/M)^σ*a^α)^(1/(1-η))
            end
            # Compute HHH_T:
            h_T[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
            #spl_e_T=Spline1D(a_grid,e_T[:,iH,iG],bc="extrapolate")
            spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")
            #HHH_T_0[iG]= (2*H/M)^β*s_T^(ϕ*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α*β/σ)*maximum([spl_e_T(aa),0.0])^(η*β/σ),lowbnd,upbnd)[1]
            HHH_T_0[iG]= (2*H/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]
            # Compute H_O: (note that exponent sigma/beta is absent):
            h_O[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG].^η
            spl_e_O=Spline1D(a_grid,e_O[:,iH,iG],bc="extrapolate")
            spl_T_thresh=Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")
            H_O_0[iG]=(2*H/M)^σ*s_O^ϕ*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^α*maximum([spl_e_O(aa),0.0])^η,lowbnd,upbnd)[1]

            # Total earnings of teachers in each group:
            E_T[iH,iG] =  H_O[iH]*A/HH_T[iH]/f1[iG]*(2*H/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]
            #H_O[iH]*A/HH_T[iH]/f1[iG]*(2*H/M)^β*s_T^(ϕ*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α*β/σ)*maximum([spl_e_T(aa),0.0])^(η*β/σ),lowbnd,upbnd)[1]
            # Total earnings of others in each group:
            E_O[iH,iG] = (1-τ_w[iG])*A*(2*H/M)^σ*s_O^ϕ*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^α*maximum([spl_e_O(aa),0.0])^η,lowbnd,upbnd)[1]
        end
        HH_T[iH]=sum(HHH_T_0.*gm)
        H_O[iH]=sum(H_O_0.*gm)

        t[iH,2]=sum(E_T[iH,:].*gm)/(sum(E_T[iH,:].*gm)+sum(E_O[iH,:].*gm))
        convT = abs(t[iH,2]-t[iH,1])
        println("convT=",round(convT,digits=3))
        println("t[iH,:]=",round.(t[iH,:],digits=3))
        println("HH_T[iH]=",round(HH_T[iH],digits=3))
        t[iH,1] = min(t[iH,2],1-1e-3) # ν*t[iH,1]+(1-ν)*t[iH,2]

        iterT = iterT+1
    end
end

#plt1 = plot(a_grid,a_T_thresh[:,1,1])
#plt2 = plot(a_grid,a_O_thresh[:,1,1])

for iH in 1:length(H_grid)
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Initiate 'while' loop ('not indexed') over the tax rate:
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            a=a_grid[ia]
            e_T[ia,iH,iG]=((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*a^(α/(σ/β-η))*f3[iG]/f1[iG]/f2[iG]
            e_O[ia,iH,iG]=((1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG])*A*η*s_O^ϕ*(2*H/M)^σ*a^α)^(1/(1-η))
        end
        # Compute HHH_T:
        h_T[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
        #spl_e_T=Spline1D(a_grid,e_T[:,iH,iG],bc="extrapolate")
        spl_O_thresh=Spline1D(a_grid, a_O_thresh[:,iG],bc="extrapolate")
        HHH_T_0[iG]= (2*H/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]
        #HHH_T_0[iG]= (2*H/M)^β*s_T^(ϕ*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α*β/σ)*maximum([spl_e_T(aa),0.0])^(η*β/σ),lowbnd,upbnd)[1]
        # Compute H_O: (note that exponent sigma/beta is absent):
        h_O[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG].^η
        spl_e_O=Spline1D(a_grid,e_O[:,iH,iG],bc="extrapolate")
        spl_T_thresh=Spline1D(a_grid, a_T_thresh[:,iG],bc="extrapolate")

        H_O_0[iG]= (2*H/M)^σ*s_O^ϕ*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^α*maximum([spl_e_O(aa),0.0])^η,lowbnd,upbnd)[1]

        # Total earnings of teachers in each group:
        E_T[iH,iG] = H_O[iH]*A/HH_T[iH]/f1[iG]*(2*H/M)^β*s_T^(ϕ*β/σ)*(((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*f3[iG]/f1[iG]/f2[iG])^(η*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α/(σ/β-η)),lowbnd,upbnd)[1]
        #H_O[iH]*A/HH_T[iH]/f1[iG]*(2*H/M)^β*s_T^(ϕ*β/σ)*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_O_thresh(aa))*aa^(α*β/σ)*maximum([spl_e_T(aa),0.0])^(η*β/σ),lowbnd,upbnd)[1]
        # Total earnings of others in each group:
        E_O[iH,iG] = (1-τ_w[iG])*A*(2*H/M)^σ*s_O^ϕ*quadgk(aa -> pdf(LogNormal(mean_a,std_a),aa)*cdf(LogNormal(mean_a,std_a),spl_T_thresh(aa))*aa^α*maximum([spl_e_O(aa),0.0])^η,lowbnd,upbnd)[1]
    end
    HH_T[iH]=sum(HHH_T_0.*gm)
    H_O[iH]=sum(H_O_0.*gm)
end

#### Group-specific barriers
HH_T_cf=zeros(length(H_grid))
H_O_cf=zeros(length(H_grid))
for iH in 1:length(H_grid)
    H=H_grid[iH]
    for iG in 1:n_g
        HH_T_cf[iH]+=(gm[iG]/sum(gm))*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[iH,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[iG])/(1+τ_e[iG]))^(η*η*β/σ/(1-η))*(2*H/M)^(β/(1-η))*f2[iG]^(1-β*η/σ)*(f3[iG]/f1[iG])^(η*β/σ)
        H_O_cf[iH]+=(gm[iG]/sum(gm))*(A*η*(1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*s_O^(ϕ/(1-η))*(2*H/M)^(σ/(1-η))*f3[iG]
    end
end

N=zeros(length(a_grid),length(H_grid),n_g)
for iH in 1:length(H_grid)
    H=HH_T[iH] # H_grid[iH]
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            N[ia,iH,iG]=(M/2)/H*(h_T[ia,iH,iG])^(β/σ)
        end
    end
end
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
save(filename1,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"α",α,"β",β,"σ",σ,"η",η,"A",A,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2,"HH_T_cf",HH_T_cf,"H_O_cf",H_O_cf, "N", N,"EN_agg",EN_agg,"cvN",cvN,"gm",gm)
filename2 = string("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/previousParameterization.jld")
save(filename2,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"α",α,"β",β,"σ",σ,"η",η,"A",A,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2,"HH_T_cf",HH_T_cf,"H_O_cf",H_O_cf, "N", N,"EN_agg",EN_agg,"cvN",cvN,"gm",gm)

# Find steady state
# using Roots
# f(x)=x-gm[1]/sum(gm)*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[1,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[1])/(1+τ_e[1]))^(η*η*β/σ/(1-η))*(2/M)^(β/(1-η))*f2[1]^(1-β*η/σ)*(f3[1]/f1[1])^(η*β/σ)*x^(β/(1-η))-gm[2]/sum(gm)*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[1,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[2])/(1+τ_e[2]))^(η*η*β/σ/(1-η))*(2/M)^(β/(1-η))*f2[2]^(1-β*η/σ)*(f3[2]/f1[2])^(η*β/σ)*x^(β/(1-η))
# find_zero(f,2)
