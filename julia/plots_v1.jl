###########################
########## PLOTS ##########
###########################
#a_T_thresh_start=1.0.*a_T_thresh


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
