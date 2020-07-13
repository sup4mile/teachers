using Distributions, Dierckx, QuadGK, JLD, Plots, LaTeXStrings

# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
occ = ["teaching","other"]
ν = 0.9

n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations
τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
τ_w[1,1] = 0.1
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_e[1,1] = 0.1

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)
α=.75
η=.75
β=.15
σ=.13
μ=1/2
ϕ=1/3
A=2 # productivity in 'Other'
θ=3

H_grid=range(.05,stop=1.4,length=17)
a_grid=quantile.(Frechet(θ),.005:.015:1)
s_T=μ*ϕ/(μ*ϕ+σ/β-η)
s_O=μ*ϕ/(μ*ϕ+1-η)
cnst=(1-η)/(1-η*β/σ)*((1-s_O)/(1-s_T))^(1/μ)
maxiterAA=200
tol_a=1e-6

t=zeros(length(H_grid),2)

a_O_thresh = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
a_T_thresh_new = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results/results_2_groups_tauW=[0.1 0.0]_tauE=[0.1 0.0]_beta-to-sigma=1.15_A=2.0.jld")
a_T_thresh=d["a_T_thresh"]

f_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

f1 = Array{Float64,1}(undef,n_g)
f2 = Array{Float64,1}(undef,n_g)
f3 = Array{Float64,1}(undef,n_g)

e_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
e_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

h_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
h_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)

spl_T = Array{Spline1D,1}(undef,n_g)
spl_O = Array{Spline1D,1}(undef,n_g)

HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,1}(undef,n_g)

# Note that the threshold doesn't deped on H
for iH in 1:1#length(H_grid)
    println("H: gridpoint ", iH," of ", length(H_grid)," (H = ", H_grid[iH],")")
    println("")
    H = H_grid[iH] # today's aggregate human capital in teaching
    for iG in 1:n_g
        iAA=0
        conv_a=1000
        while conv_a>tol_a && iAA<maxiterAA
            global a_T_thresh

            iAA+=1
            println("Iteration=",iAA)

            for ia in 1:length(a_grid)
                f_O[ia,iH,iG]=quadgk(aa -> pdf(Frechet(θ),aa),0.0,a_T_thresh[ia,iH,iG])[1]
            end

            spl_1=Spline1D(a_grid[:],a_grid[:].^(α/(1-η)).*pdf(Frechet(θ),a_grid[:]).*f_O[:,iH,iG])
            f1[iG]=quadgk(aa -> spl_1(aa),0.0,1e3)[1]
            spl_2=Spline1D(a_grid[:], pdf(Frechet(θ),a_grid[:]).*f_O[:,iH,iG])
            f2[iG]=quadgk(aa -> spl_2(aa),0.0,1e3)[1]

            spl_inv = Spline1D(a_T_thresh[:,iH,iG],a_grid)
            for ia in 1:length(a_grid)
                a_O_thresh[ia,iH,iG] = spl_inv(a_grid[ia])
                f_T[ia,iH,iG]=quadgk(aa -> pdf(Frechet(θ),aa),0.0,a_O_thresh[ia,iH,iG])[1]
            end

            spl_3=Spline1D(a_grid[:],a_grid[:].^(α*β/(σ-β*η)).*pdf(Frechet(θ),a_grid[:]).*f_T[:,iH,iG])
            f3[iG]=quadgk(aa -> spl_3(aa),0.0,1e3)[1]

            ratio=f1[iG]/f2[iG]/f3[iG]

            a_T_thresh_new[:,iH,iG]=((1-τ_w[iG])*cnst*ratio.*(a_grid.^(α/(1-η)))).^((σ/β-η)/α)

            conv_a=maximum(abs.(a_T_thresh_new[:,iH,iG] - a_T_thresh[:,iH,iG]))
            println("error=", conv_a)
            a_T_thresh[:,iH,iG]=a_T_thresh_new[:,iH,iG]
        end
    end
end
for iG in 1:n_g
    for ia in 1:length(a_grid)
        a_T_thresh[ia,:,iG].=a_T_thresh[ia,1,iG]
    end
end

plt1 = plot(a_grid,a_T_thresh[:,1,1])
plt2 = plot(a_grid,a_O_thresh[:,1,1])

for iH in 1:length(H_grid)
    H=H_grid[iH]
    for iG in 1:n_g
        for ia in 1:length(a_grid)
            a=a_grid[ia]
            e_T[ia,iH,iG]=((1-τ_w[iG])/(1+τ_e[iG]))^(η/(1-η))*((1-t[iH,1])*A*η*s_O^ϕ*(2*H/M)^σ)^(1/(1-η))*β/σ*a^(α*σ/(σ-β*η))*f1[iG]/f2[iG]/f3[iG]
            e_O[ia,iH,iG]=((1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG])*A*η*s_O^ϕ*(2*H/M)^σ*a^α)^(1/(1-η))
        end
        # Fit spline to (h^T)^(β/σ) and 'a_grid' to compute HHH_T:
        h_T[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_T.^ϕ .* e_T[:,iH,iG].^η
        spl_T[iG] = Spline1D(a_grid,h_T[:,iH,iG].^(β/σ))
        # Fit spline to h^O and 'a_grid' (note that exponent sigma/beta is absent):
        h_O[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_O.^ϕ .* e_O[:,iH,iG].^η
        spl_O[iG] = Spline1D(a_grid,h_O[:,iH,iG])
        function f_T_pwl(aa)
            # Locate closest gridpoint:
            ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
            mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
            return (mu*f_T[ind+1,iH,iG] + (1-mu)*f_T[ind,iH,iG])*gm[iG]
        end
        function f_O_pwl(aa)
            # Locate closest gridpoint:
            ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
            mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
            return (mu*f_O[ind+1,iH,iG] + (1-mu)*f_O[ind,iH,iG])*gm[iG]
        end
        HHH_T_0[iG]= quadgk(aa -> spl_T[iG](aa)*f_T_pwl(aa),a_grid[1],a_grid[end])[1]
        H_O_0[iG]= quadgk(aa -> spl_O[iG](aa)*f_O_pwl(aa),a_grid[1],a_grid[end])[1]
    end
    HH_T[iH]=sum(HHH_T_0.*gm)
    H_O[iH]=sum(H_O_0.*gm)
end


#### Group-specific barriers!!!
HH_T=zeros(length(H_grid))
H_O=zeros(length(H_grid))
for iH in 1:length(H_grid)
    H=H_grid[iH]
    for iG in 1:n_g
        HH_T[iH]+=gm[iG]*f3[iG]^((σ-β*η)/σ)*(f1[iG]/f2[iG])^(η*β/σ)*s_T^(ϕ*β/σ)*s_O^(ϕ*β*η/(σ-σ*η))*((1-t[iH,1])*η*A)^(η*β/(σ-σ*η))*(β/σ)^(η*β/σ)*((1-τ_w[iG])/(1+τ_e[iG]))^(η*η*β/(σ-σ*η))*(2*H/M)^(β/(1-η))
        H_O[iH]+=gm[iG]*f1[iG]*((1-t[iH,1])*(1-τ_w[iG])/(1+τ_e[iG])*η*A)^(η/(1-η))*s_O^(ϕ/(1-η))*(2*H/M)^(σ/(1-η))
    end
end
