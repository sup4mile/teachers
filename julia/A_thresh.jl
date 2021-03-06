using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings

# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
occ = ["teaching","other"]
ν = 0.9

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations
τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[1,1] = 0.0
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_e[1,1] = 0.0

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)
α=.75
η=.75
β=.15
σ=.16
μ=1/2
ϕ=1/3
A=2 # productivity in 'Other'
θ=3

#H_grid=collect(range(.05,stop=1.4,length=17))
H_grid=collect(range(.005,stop=1.0,length=20))
# a_grid=quantile.(Frechet(θ),.005:.015:1)
mean_a=0
std_a=1
a_grid=exp.(append!(quantile.(Normal(mean_a,std_a),.001:.001:.004), quantile.(Normal(mean_a,std_a),.005:.015:1)))
temp=exp.(append!(quantile.(Normal(mean_a,std_a),.001:.001:.004), quantile.(Normal(mean_a,std_a),.005:.015:1)))
s_T=μ*ϕ/(μ*ϕ+σ/β-η)
s_O=μ*ϕ/(μ*ϕ+1-η)
cnst=(1-η)/(1-η*β/σ)*((1-s_O)/(1-s_T))^(1/μ)

maxiterAA=200
tol_a=1e-6

tolT= 2.5e-3
maxiterT = 100

t=zeros(length(H_grid),2)

a_T_thresh = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
a_O_thresh = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
a_T_thresh_new = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
#d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_tauW=[0.1 0.0]_tauE=[0.1 0.0]_beta-to-sigma=1.15_A=2.0.jld")
#a_T_thresh=d["a_T_thresh"]
for iH in 1:length(H_grid)
    for iG in 1:n_g
        a_T_thresh[:,iH,iG]=temp.^0.7
    end
end

f_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

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

wb_T = Array{Float64,2}(undef,length(H_grid),n_g)
wb_O = Array{Float64,2}(undef,length(H_grid),n_g)
E_T = Array{Float64,2}(undef,length(H_grid),n_g)
E_O = Array{Float64,2}(undef,length(H_grid),n_g)

HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,1}(undef,n_g)

HH_T = Array{Float64,1}(undef,length(H_grid))
H_O = Array{Float64,1}(undef,length(H_grid))

mass_T = Array{Float64,2}(undef,length(H_grid),n_g)
mass_O = Array{Float64,2}(undef,length(H_grid),n_g)
f_1 = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_2 = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_1_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_2_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

# Note that the threshold doesn't deped on H
for iH in 1:length(H_grid)
    println("H: gridpoint ", iH," of ", length(H_grid)," (H = ", H_grid[iH],")")
    println("")
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Initiate 'while' loop ('not indexed') over the tax rate:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        println("    T: iteration ",iterT)
        for iG in 1:n_g
            iAA=0
            conv_a=1000
            while conv_a>tol_a && iAA<maxiterAA
                global a_T_thresh

                iAA+=1
                println("Iteration=",iAA)

                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    f_1[ia,iH,iG]= 1-cdf(LogNormal(mean_a,std_a),a_T_thresh[ia,iH,iG]) # Fraction of teachers given a_O
                    f_2[ia,iH,iG]= cdf(LogNormal(mean_a,std_a),a_T_thresh[ia,iH,iG]) # Fraction of other given a_O
                    f_O[ia,iH,iG]=pdf(LogNormal(mean_a,std_a),a)*f_2[ia,iH,iG] # mass of other given a_O
                end

                spl_O=Spline1D(a_grid[:], f_O[:,iH,iG],bc="extrapolate")
                f1[iG]=quadgk(aa -> maximum([spl_O(aa),0.0]),0.0,1e3)[1] # total mass of other using direct threshold
                mass_O[iH,iG]=quadgk(aa -> maximum([spl_O(aa),0.0]),0.0,1e3)[1]

                spl_T=Spline1D(a_grid[:], pdf.(LogNormal(mean_a,std_a),a_grid[:]).*f_1[:,iH,iG],bc="extrapolate")
                f4[iG]=quadgk(aa -> maximum([spl_T(aa),0.0]),0.0,1e3)[1] # total mass of teachers using direct threshold

                mass_T[iH,iG]=quadgk(aa -> maximum([spl_T(aa),0.0]),0.0,1e3)[1]

                f3[iG]=quadgk(aa -> maximum([spl_O(aa),0.0])*aa^(α/(1-η)),0.0,1e3)[1]

                spl_inv = Spline1D(a_T_thresh[:,iH,iG],a_grid,bc="extrapolate")
                for ia in 1:length(a_grid)
                    a = a_grid[ia]
                    a_O_thresh[ia,iH,iG] = maximum([spl_inv(a_grid[ia]),0.0])
                    f_1_T[ia,iH,iG]= cdf(LogNormal(mean_a,std_a),a_O_thresh[ia,iH,iG]) # Fraction of teachers given a_T
                    f_2_T[ia,iH,iG]= 1-cdf(LogNormal(mean_a,std_a),a_O_thresh[ia,iH,iG]) # Fraction of other given a_T
                    f_T[ia,iH,iG]=pdf(LogNormal(mean_a,std_a),a)*f_1_T[ia,iH,iG]  # mass of teachers given a_T
                end

                spl_O_inv=Spline1D(a_grid[:], pdf.(LogNormal(mean_a,std_a),a_grid[:]).*f_2_T[:,iH,iG],bc="extrapolate")
                f1_inv[iG]=quadgk(aa -> maximum([spl_O_inv(aa),0.0]),0.0,1e3)[1] # total mass of other using inverse threshold

                spl_T_inv=Spline1D(a_grid[:], f_T[:,iH,iG],bc="extrapolate")
                f4_inv[iG]=quadgk(aa -> maximum([spl_T_inv(aa),0.0]),0.0,1e3)[1] # total mass of teachers using inverse threshold

                f2[iG]=quadgk(aa -> maximum([spl_T_inv(aa),0.0])*aa^(α/(σ/β-η)),0.0,1e3)[1]

                ratio=f1[iG]*f2[iG]/f3[iG]

                a_T_thresh_new[:,iH,iG]=((1-τ_w[iG])*cnst*ratio.*(a_grid.^(α/(1-η)))).^((σ/β-η)/α)

                conv_a=maximum(abs.(a_T_thresh_new[:,iH,iG] - a_T_thresh[:,iH,iG]))
                println("error=", conv_a)
                a_T_thresh[:,iH,iG]=a_T_thresh_new[:,iH,iG]
            end

            for ia in 1:length(a_grid)
                a=a_grid[ia]
                e_T[ia,iH,iG]=((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*a^(α/(σ/β-η))*f3[iG]/f1[iG]/f2[iG]
                e_O[ia,iH,iG]=((1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG])*A*η*s_O^ϕ*(2*H/M)^σ*a^α)^(1/(1-η))
            end

            # Compute HHH_T:
            h_T[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
            spl_h_T=Spline1D(a_grid,h_T[:,iH,iG].^(β/σ),bc="extrapolate")
            spl_T_inv=Spline1D(a_grid[:], f_T[:,iH,iG],bc="extrapolate")
            HHH_T_0[iG]= quadgk(aa -> maximum([spl_h_T(aa),0.0])*maximum([spl_T_inv(aa),0.0]),0.0,1e3)[1]
            # Compute H_O: (note that exponent sigma/beta is absent):
            h_O[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG].^η
            spl_h_O=Spline1D(a_grid,h_O[:,iH,iG],bc="extrapolate")
            spl_O=Spline1D(a_grid[:], f_O[:,iH,iG],bc="extrapolate")
            H_O_0[iG]= quadgk(aa -> maximum([spl_h_O(aa),0.0])*maximum([spl_O(aa),0.0]),0.0,1e3)[1]
        end
        HH_T[iH]=sum(HHH_T_0.*gm)
        H_O[iH]=sum(H_O_0.*gm)
        for iG in 1:n_g
            for ia in 1:length(a_grid)
                ω[ia,iH,iG] = H_O[iH]*A*h_T[ia,iH,iG].^(β/σ)/(HH_T[iH]*f1[iG])
            end
            spl_ω=Spline1D(a_grid,ω[:,iH,iG],bc="extrapolate")
            spl_T_inv=Spline1D(a_grid[:], f_T[:,iH,iG],bc="extrapolate")
            # Total earnings of teachers in each group:
            E_T[iH,iG] = quadgk(aa -> maximum([spl_ω(aa),0.0])*maximum([spl_T_inv(aa),0.0]),0.0,1e3)[1]

            spl_h_O=Spline1D(a_grid,h_O[:,iH,iG],bc="extrapolate")
            spl_O=Spline1D(a_grid[:], f_O[:,iH,iG],bc="extrapolate")
            # Total earnings of others in each group:
            E_O[iH,iG] = quadgk(aa -> (1-τ_w[iG])*A*maximum([spl_h_O(aa),0.0])*maximum([spl_O(aa),0.0]),0.0,1e3)[1]
        end
        t[iH,2]=sum(E_T[iH,:])/(sum(E_T[iH,:])+sum(E_O[iH,:]))
        convT = abs(t[iH,2]-t[iH,1])
        println("convT=",round(convT,digits=3))
        println("t[iH,:]=",round.(t[iH,:],digits=3))
        t[iH,1] = min(t[iH,2],1-1e-3) # ν*t[iH,1]+(1-ν)*t[iH,2]
        iterT = iterT+1
    end
end


#plt1 = plot(a_grid,a_T_thresh[:,1,1])
#plt2 = plot(a_grid,a_O_thresh[:,1,1])

#### Group-specific barriers
HH_T_cf=zeros(length(H_grid))
H_O_cf=zeros(length(H_grid))
for iH in 1:length(H_grid)
    H=H_grid[iH]
    for iG in 1:n_g
        HH_T_cf[iH]+=gm[iG]*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/σ/(1-η))*((1-t[iH,1])*η*A)^(η*β/σ/(1-η))*(β/σ)^(η*β/σ)*((1-τ_w[iG])/(1+τ_e[iG]))^(η*η*β/σ/(1-η))*(2*H/M)^(β/(1-η))*f2[iG]^(1-β*η/σ)*(f1[iG]/f3[iG])^(η*β/σ)
        H_O_cf[iH]+=gm[iG]*(A*η*(1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*s_O^(ϕ/(1-η))*(2*H/M)^(σ/(1-η))*f3[iG]
    end
end

filename = string("C:/Users/iuliia/Desktop/teachers/julia/results_new/results_",n_g,"_groups_tauW=",τ_w,"_tauE=",τ_e,"_beta-to-sigma=",round(β/σ,digits=2),"_A=",round(A,digits=2),".jld")
save(filename,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"β",β,"σ",σ,"A",A,"HH_T",HH_T,"H_O",H_O,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2,"HH_T_cf",HH_T_cf,"H_O_cf",H_O_cf)
