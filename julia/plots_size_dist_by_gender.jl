using JLD, Plots, Distributions, LaTeXStrings
# Number of plots:
pltnm = 10

alpha = Vector{Float64}(undef,pltnm)
beta = Vector{Float64}(undef,pltnm)
eta = Vector{Float64}(undef,pltnm)
sigma = Vector{Float64}(undef,pltnm)
fracteach = Vector{Float64}(undef,pltnm)

pyplot()
# Plot #1.1:
d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results_new/results_2_groups_τW=[-3.0 -4.0]_τE=[0.0 0.0]_A=10.0_α=0.1_β=0.3_η=0.1_σ=0.4.jld")
H_grid = d["H_grid"]
a_grid = d["a_grid"]
gm = d["gm"]
a_T_thresh = d["a_T_thresh"]
a_O_thresh = d["a_O_thresh"]
β = d["β"]
σ = d["σ"]
α = d["α"]
η = d["η"]
τ_w=d["τ_w"]
H_grid = d["H_grid"]
HH_T = d["HH_T"]
# HH_T_cf = d["HH_T_cf"]
mass_O = d["mass_O"]
mass_T = d["mass_T"]
f_1 = d["f_1"]
N=d["N"]
EN_agg=d["EN_agg"]
cvN=d["cvN"]
# plt1 = plot(plt1,a_grid,N[:,1,1],grid=:none,label=latexstring("\\alpha = "*string(round(α,digits=2))*", \\beta = "*string(round(β,digits=2))*", \\eta = "*string(round(η,digits=2))*", \\sigma = "*string(round(σ,digits=2))*", E(N) = "*string(round(EN_agg,digits=2))*", c.v. = "*string(round(cvN,digits=2))*", \\%T = "*string(round((gm./(sum(gm))*mass_T)[1],digits=2))*",S/T = "*string(round(M*((gm./(sum(gm))*mass_T)[1])^(-1),digits=2))))
plt1 = plot(a_grid,N[:,1,1],grid=:none,label=latexstring("τ_w^O = "*string(τ_w[1])*", N = "*string(round(EN[1],digits=2))*", c.v. = "*string(round(cvN[1],digits=2))*", \\%T = "*string(round((mass_T)[1],digits=2))*", S/T = "*string(round(((gm'./sum(gm).*EN)./(gm'.*mass_T))[1],digits=2))))
plt1 = plot!(plt1,a_grid,N[:,1,2],grid=:none,label=latexstring("τ_w^O = "*string(τ_w[2])*", N = "*string(round(EN[2],digits=2))*", c.v. = "*string(round(cvN[2],digits=2))*", \\%T = "*string(round((mass_T)[2],digits=2))*", S/T = "*string(round(((gm'./sum(gm).*EN)./(gm'.*mass_T))[2],digits=2))))
# plt1 = plot(a_grid,N[:,1,1],grid=:none,label=latexstring("\\% in  T = "*string(round((gm./(sum(gm))*mass_T)[1],digits=2))))
# fracteach[1]=sum(mass_T[:])/(sum(mass_T[:])+sum(mass_O[:]))
