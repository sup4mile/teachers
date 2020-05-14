using JuMP, Ipopt, Dierckx, NLsolve, Distributions, Roots, QuadGK, JLD, Plots

include("optiminv_T.jl")
include("optiminv_O.jl")
include("opt_thresh.jl")

# Switch to get values from previous parameterization:
getvalues = "on"
if getvalues == "on"
    # Change path to match architecture on local machine:
    d = load("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results/results_2_groups_tauW=[0.0 0.0]_tauE=[0.0 0.0]_beta-to-sigma=1.0_A=1.0.jld")
    # d = load("Z:/teachers_julia/results_2_groups_tauW=[0.0 0.0]_tauE=[0.1 0.0].jld")
    # Extract arrays from dictionary 'd':
    H_grid = d["H_grid"]
    HH_T = d["HH_T"]
    H_O = d["H_O"]
    λ = d["λ"]
    a_T_thresh = d["a_T_thresh"]
    a_O_thresh = d["a_O_thresh"]
    e_T = d["e_T"]
    s_T = d["s_T"]
    e_O = d["e_O"]
    s_O = d["s_O"]
    t = d["t"]
    E_O = d["E_O"]
    E_T = d["E_T"]
    mass_O = d["mass_O"]
    mass_T = d["mass_T"]
    f_1 = d["f_1"]
    f_2 = d["f_2"]
end

# Parameters:
M=1.0 # population measure
g = ["female","male"] # no discrimination / no barrier group last
occ = ["teaching","other"]
ν = 0.8 # weight on HH_T (expected H') vs HHH_T (realized H') when we update HH_T
if getvalues == "on"
    # Check if number of groups in previous results is equal to number of elements in g:
    if size(HH_T,2) < size(g,1)
        for idim in size(HH_T,2)+1:size(g,1)
            global HH_T,a_T_thresh,a_O_thresh,e_T,s_T,e_O,s_O
            # Concatenate arrays with results from legacy groups:
            a_T_thresh = cat(a_T_thresh[:,:,idim-1],a_T_thresh[:,:,idim-1],dims=3)
            a_O_thresh = cat(a_O_thresh[:,:,idim-1],a_O_thresh[:,:,idim-1],dims=3)
            e_T = cat(e_T[:,:,idim-1],e_T[:,:,idim-1],dims=3)
            s_T = cat(s_T[:,:,idim-1],s_T[:,:,idim-1],dims=3)
            e_O = cat(e_O[:,:,idim-1],e_O[:,:,idim-1],dims=3)
            s_O = cat(s_O[:,:,idim-1],s_O[:,:,idim-1],dims=3)
        end
    end
end
n_g=length(g) # number of "groups"
n_occ=length(occ) # number of occupations
τ_w=zeros(n_occ-1,n_g) # n_g-element vector of labor market discrimination in 'O' (relative to 'T')
# τ_w[1,1] = 0.1
τ_e=zeros(n_occ-1,n_g) # n_g-element vector of education barriers in 'O' (relative to 'T')
τ_e[1,1] = 0.1

gm=zeros(1,n_g) # measure of in individuals in groups 1,...,n_g
for iG in 1:n_g-1
    gm[iG] = M/(2*n_g)
end
gm[end] = M/2 - sum(gm)
α=.75
β=.5
σ=.4 # .625
η=.5
μ=1/2
ϕ=1/3
A=1 # productivity in 'Other'
if getvalues != "on"
    t=zeros(length(H_grid),2) # lump sum tax
    t=fill(0.65,(length(H_grid),2))
end
θ=3 # skill dispersion (high θ = small dispersion)

# State space (endogenous state variables): aggregate human capital in teaching today:
H_grid=collect(0.35:0.05:1.15) # Spline doesn't work for values below 0.4. Why?
# State space (exogenous state variables): idiosyncratic ability
a_grid=quantile.(Frechet(θ),.005:.015:1)
# Initialize arrays for 'VV', 's', and 'e':
VV_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
VV_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
h_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
h_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
ω = Array{Float64,4}(undef,length(a_grid),length(H_grid),n_g,2)
spl_T = Array{Spline1D,1}(undef,n_g)
spl_O = Array{Spline1D,1}(undef,n_g)
spl_TT = Array{Spline1D,1}(undef,n_g)
spl_ω = Array{Spline1D,1}(undef,n_g)
spl_h_O = Array{Spline1D,1}(undef,n_g)
T_g = Array{Float64,1}(undef,n_g)
E_T = Array{Float64,1}(undef,n_g)
E_O = Array{Float64,1}(undef,n_g)
HHH_T_0 = Array{Float64,1}(undef,n_g)
H_O_0 = Array{Float64,1}(undef,n_g)
H_T_0 = Array{Float64,1}(undef,n_g)
H_T = Array{Float64,1}(undef,length(H_grid))
mass_T = Array{Float64,2}(undef,length(H_grid),n_g)
mass_O = Array{Float64,2}(undef,length(H_grid),n_g)
f_1 = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
f_2 = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
# Only initialize the following arrays if not already loaded from previous solution:
if getvalues != "on"
    s_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    e_T = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    a_T_thresh = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    s_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    e_O = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    a_O_thresh = Array{Float64,3}(undef,length(a_grid),length(H_grid),n_g)
    # Expected aggregate human capital in 'Teaching' tomorrow:
    HH_T = Array{Float64,1}(undef,length(H_grid))
    # Aggregate human capital in 'Other' tomorrow:
    H_O = Array{Float64,1}(undef,length(H_grid))
    λ = Array{Float64,2}(undef,length(H_grid),2) # scale parameter for teachers' wage profile (used when they are not paid the marginal product)
    # Initial guess for human capital investment decisions (if no prior results have been loaded, i.e. if 'getvalues' = "off"); initialize for all groups!
    s_T[1,1,:] = .17 * ones(n_g,1)
    e_T[1,1,:] = .011 * ones(n_g,1)
    s_O[1,1,:] = .25 * ones(n_g,1)
    e_O[1,1,:] = .034 * ones(n_g,1)
end
# Switch for targeting measure of teachers:
# 1 - update λ to match targetmass_T
# 0 - fix λ, no target and teachers are paid their marginal product
targetM = 0
# Set target measure (if needed):
# if targetM == 1
targetmass_T = 0.1 # Target measure of teachers
maxupdate = 0.01 # Maxumum update for lambda
# end

tolH = 2e-3 # convergence criterion for HH_T (expected) and HHH_T (realized)
tolω = 2e-1 #tolH
tolT= 2e-2
maxiterT = 100
maxiterHH = 100

tolM = 1e-3 # convergence criterion for mass_T (actual) and targetmass_T (targeted, if at all)
maxiterM = 100

iH = 1 # index for loop over H today
maxiterH = length(H_grid)+1
iG = 1 # index for loop over groups
ia = 1 # index for loop over idiosyncratic abilities

signH = 0
flag = 0 # Flag for steady state (i.e. H = HHH_T), change the signH = sign(HH_T-H)
# Initiate 'while' loop over 'H' (indexed by 'iH') over today's aggregate human capital in teaching; the length of the loop is limited by maxiterH (not by the size of 'H'):
@time while iH < maxiterH
    println("H: gridpoint ", iH," of ", length(H_grid)," (H = ", H_grid[iH],")")
    println("")
    H = H_grid[iH] # today's aggregate human capital in teaching
    # Given 'iG'=1 and 'ia'=1, update initial guesses for "next" value of 'H' on grid:
    if iH > 1
        s_T[1,iH,1] = s_T[1,iH-1,1]
        e_T[1,iH,1] = e_T[1,iH-1,1]
        s_O[1,iH,1] = s_O[1,iH-1,1]
        e_O[1,iH,1] = e_O[1,iH-1,1]
        # println(s_T[1,iH,1])
    end
    # Initial guess for tomorrow's expected aggregate human capital and 'λ':'
    if getvalues != "on" && iH == 1
        HH_T[iH] = H_grid[iH]
        H_O[iH] = H_grid[iH]
        λ[iH,1] = H_O[iH]*A/(M/2)
    end
    # Initiate 'while' loop ('not indexed') over the lump sum tax:
    convT = 1; iterT = 1
    while convT > tolT && iterT < maxiterT
        println("    T: iteration ",iterT)
        # Initiate 'while' loop ('not indexed') over the measure of teachers:
        convM = 1; iterM = 1
        while convM > tolM && iterM < maxiterM
            if targetM == 1
                println("  Measure of teachers: iteration ", iterM)
            end
            # Initiate 'while' loop (indexed by 'iHH') over expected (HH_T) and realized (HHH_T) human capital in teaching
            convH = 1 # convergence criterion (1)
            convω = 1 # convergence criterion (2)
            iHH = 1 # iteration counter for HH_T
            while (convH > tolH || convω > tolω) && iHH < maxiterHH
                println("    H': iteration ",iHH)#," (H = ",H,", E[H'] = ", round.(HH_T[iH],digits=3),", E[H_O'] = ",round.(H_O[iH],digits=3),")")
                # Initiate 'for' loop over groups (indexed by 'iG'):
                for iG in 1:n_g
                    if iG > 1
                        # init_s_T[iG] = s_T[1,iH,iG]
                        # init_e_T[iG] = e_T[1,iH,iG]
                        # init_s_O[iG] = s_O[1,iH,iG]
                        # init_e_O[iG] = e_O[1,iH,iG]
                    end
                    # Initiate 'for' loop (indexed by 'ia') over idiosyncratic abilities:
                    for ia in 1:length(a_grid)
                        a = a_grid[ia]
                        # For given 'iH' and 'iG', update initial guess using results from grid point 'ia-1':
                        if ia > 1
                            s_T[ia,iH,iG] = s_T[ia-1,iH,iG]
                            e_T[ia,iH,iG] = e_T[ia-1,iH,iG]
                            s_O[ia,iH,iG] = s_O[ia-1,iH,iG]
                            e_O[ia,iH,iG] = e_O[ia-1,iH,iG]
                        end
                        VV_T[ia,iH,iG],s_T[ia,iH,iG],e_T[ia,iH,iG] = optiminv_T(λ[iH,1],M,HH_T[iH],H,α,β,σ,η,μ,ϕ,0,0,t[iH,1],a,s_T[ia,iH,iG],e_T[ia,iH,iG])
                        ω[ia,iH,iG,1] = λ[iH,1]*M/(2*HH_T[iH])*a^α*s_T[ia,iH,iG]^ϕ*e_T[ia,iH,iG]^η*(2*HH_T[iH]/M)^σ
                        VV_O[ia,iH,iG],s_O[ia,iH,iG],e_O[ia,iH,iG] = optiminv_O(A,M,H,α,σ,η,μ,ϕ,τ_w[iG],τ_e[iG],t[iH,1],a,s_O[ia,iH,iG],e_O[ia,iH,iG])
                    end
                    # Compute occupational threshold using VV_O and VV_T for given ('iH','iHH','iG'):
                    # for iG in 1:n_g
                    println("      iter(group) = ", iG,": ", g[iG])
                    for ia in 1:length(a_grid) # 'a_grid[ia]' indexes ability in 'other'(?)
                        a_T_thresh[ia,iH,iG] = opt_thresh(a_grid,VV_T[:,iH,iG],a_grid[ia],VV_O[ia,iH,iG])
                    end
                    # Compute the inverse occupational threshold and density of teachers:
                    spl = Spline1D(a_T_thresh[:,iH,iG],a_grid)
                    for ia in 1:length(a_grid)
                        a = a_grid[ia] # indexes ability in teaching
                        # Inverse threshold:
                        a_O_thresh[ia,iH,iG] = spl(a)
                        # Density:
                        if a_T_thresh[ia,iH,iG]>0.0
                            f_1[ia,iH,iG], err = quadgk(aa -> pdf(Frechet(θ),aa),a_T_thresh[ia,iH,iG],1e3)
                            f_2[ia,iH,iG], err = quadgk(aa -> pdf(Frechet(θ),aa),0,a_T_thresh[ia,iH,iG])
                            #f_1 = quadgk(aa -> pdf(Frechet(θ),aa),0,a_O_thresh[ia,iH,iG])
                            #f_2 = quadgk(aa -> pdf(Frechet(θ),aa),a_O_thresh[ia,iH,iG],1e3)
                        else
                            f_1[ia,iH,iG] = 1.0
                            f_2[ia,iH,iG] = 0.0
                        end
                        f_T[ia,iH,iG] = pdf(Frechet(θ),a)*f_1[ia,iH,iG]
                        f_O[ia,iH,iG] = pdf(Frechet(θ),a)*f_2[ia,iH,iG]
                    end
                    # Fit spline to (h^T)^(β/σ) and 'a_grid' to compute HHH_T:
                    spl_T[iG] = Spline1D(a_grid,((2*H/M).^σ .* a_grid.^α .* s_T[:,iH,iG].^ϕ .* e_T[:,iH,iG].^η).^(β/σ))
                    # Fit spline to h^O and 'a_grid' (note that exponent sigma/beta is absent):
                    spl_O[iG] = Spline1D(a_grid,(2*H/M).^σ .* a_grid.^α .* s_O[:,iH,iG].^ϕ .* e_O[:,iH,iG].^η)
                    h_O[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_O[:,iH,iG].^ϕ .* e_O[:,iH,iG].^η
                    # Fit another spline to h^T and 'a_grid' to compute 'H_T':
                    spl_TT[iG] = Spline1D(a_grid,(2*H/M).^σ .* a_grid.^α .* s_T[:,iH,iG].^ϕ .* e_T[:,iH,iG].^η)
                    h_T[:,iH,iG] = (2*H/M).^σ .* a_grid.^α .* s_T[:,iH,iG].^ϕ .* e_T[:,iH,iG].^η
                    function f_T_pwl(aa)
                        # Locate closest gridpoint:
                        ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
                        mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
                        return (mu*f_T[ind+1,iH,iG] + (1-mu)*f_T[ind,iH,iG])*gm[iG]
                    end
                    # function f_T_pwl(aa)
                    #     # Locate closest gridpoint:
                    #     ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
                    #     mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
                    #     return mu*f_T[ind+1,iHH,iH] + (1-mu)*f_T[ind,iHH,iH]
                    # end
                    function f_O_pwl(aa)
                        # Locate closest gridpoint:
                        ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
                        mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
                        return (mu*f_O[ind+1,iH,iG] + (1-mu)*f_O[ind,iH,iG])*gm[iG]
                    end
                    # Investments by high-ability individuals (95th percentile):
                    # println(" Group: ",g[iG])
                    # println("   High ability T: ",[s_T[end-6,1,iG] e_T[end-6,1,iG]])
                    # println("   High ability O: ",[s_O[end-6,1,iG] e_O[end-6,1,iG]])
                    # Measure of teachers:
                    global mass_T[iH,iG], err = quadgk(aa -> f_T_pwl(aa),a_grid[1],a_grid[end])
                    # println(" Measure of teachers = ",mass_T[1])
                    # Measure of workers:
                    global mass_O[iH,iG], err = quadgk(aa -> f_O_pwl(aa),a_grid[1],a_grid[end])
                    # println(" Measure of workers = ",mas_O[1])

                    # Calculate H^T[tilde] (β/σ-th moment of distribution of h_T):
                    HHH_T_0[iG], err = quadgk(aa -> spl_T[iG](aa)*f_T_pwl(aa),a_grid[1],a_grid[end])
                    # Calculate H^O:
                    H_O_0[iG], err = quadgk(aa -> spl_O[iG](aa)*f_O_pwl(aa),a_grid[1],a_grid[end])
                    # Calculate H^T:
                    H_T_0[iG], err = quadgk(aa -> spl_TT[iG](aa)*f_O_pwl(aa),a_grid[1],a_grid[end])
                end
                HHH_T = sum(HHH_T_0)
                H_O[iH] = sum(H_O_0)
                H_T[iH] = sum(H_T_0)

                # COMPUTE λ and ω (loop over 'a')!
                println("        H[T]' implied by law of motion: ", round.(HHH_T,digits=3))
                println("        H[T]' anticipated by agents: ", round.(HH_T[iH],digits=3))
                println("        H[O]' anticipated by agents: ", round.(H_O[iH],digits=3))
                println("    __________________________")
                println("")
                # Update λ':
                λ[iH,2] = H_O[iH]*A/(M/2*sum(mass_T[iH,:]))
                λ[iH,1] = ν*λ[iH,1]+(1-ν)*λ[iH,2]
                # Check for convergence of H[T]' and update:
                convH = abs(HH_T[iH]-HHH_T[1])
                HH_T[iH] = ν*HH_T[iH]+(1-ν)*HHH_T[1]
                # Update wage profile ω':
                for ia in 1:length(a_grid)
                    ω[ia,iH,iG,2] = λ[iH,1]*M/(2*HH_T[iH])*a_grid[ia]^α*s_T[ia,iH,iG]^ϕ*e_T[ia,iH,iG]^η*(2*HH_T[iH]/M)^σ
                end
                # Check for convergence of ω' and update:
                convω = maximum(abs.(ω[10:end-10,iH,iG,2] - ω[10:end-10,iH,iG,1]))
                ω[:,iH,iG,1] = ν*ω[:,iH,iG,1] + (1-ν)*ω[:,iH,iG,2]
                # Update counter for inner loop (over H'(T)):
                iHH = iHH+1
            end
            if targetM != 1
                convM = targetM*abs(targetmass_T-sum(mass_T[iH,:]))+(1-targetM)*0
            end
            iterM = iterM+1
        end
        for iG in 1:n_g
            spl_ω[iG] = Spline1D(a_grid,ω[:,iH,iG,1])
            # Density function of teachers by groups:
            function f_T_pwl(aa)
                # Locate closest gridpoint:
                ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
                mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
                return (mu*f_T[ind+1,iH,iG] + (1-mu)*f_T[ind,iH,iG])*gm[iG]
            end
            # Total wages paid to each group
            T_g[iG], er = quadgk(aa -> f_T_pwl(aa)*spl_ω[iG](aa),a_grid[1],a_grid[end])
            # Total earnings of teachers for each group
            E_T[iG], er = quadgk(aa -> (1-τ_w[iG])*f_T_pwl(aa)*spl_ω[iG](aa),a_grid[1],a_grid[end])

            spl_h_O[iG] = Spline1D(a_grid,ω[:,iH,iG,1])
            # Density function of teachers by groups:
            function f_O_pwl(aa)
                # Locate closest gridpoint:
                ind = min(maximum(findall(x -> x <= aa, a_grid)),length(a_grid)-1)
                mu = (aa - a_grid[ind])/(a_grid[ind+1]-a_grid[ind])
                return (mu*f_O[ind+1,iH,iG] + (1-mu)*f_O[ind,iH,iG])*gm[iG]
            end
            # Total earnings of other for each group
            E_O[iG], er = quadgk(aa -> (1-τ_w[iG])*f_O_pwl(aa)*A*spl_h_O[iG](aa),a_grid[1],a_grid[end])
        end
        t[iH,2]=sum(T_g)/(sum(E_T)+sum(E_O))
        convT = abs(t[iH,2]-t[iH,1])
        println("convT=",convT)
        println("t[iH,:]=",t[iH,:])
        t[iH,1] = ν*t[iH,1]+(1-ν)*t[iH,2]
        iterT = iterT+1
    end
    # println("__________________________")
    println("")
    if iH == 1
        global signH = sign(HH_T[iH]-H)
    end
    if -signH*(HH_T[iH]-H)>0 && flag == 0 # the steady state has passed
        global maxiterH = min(iH+10, length(H_grid))
        global flag = 1
    end
    if getvalues != "on"
        # Update H'(O):
        H_O[iH+1] = H_O[iH]
        # Update initial guesses for HH_T:
        HH_T[iH+1] = HH_T[iH]
        # Update 'lambda':
        λ[iH+1,1] = λ[iH,1]
    end
    global iH = iH+1
end
# Update 'convω':
# convω = (maximum(abs(ω[:,:,:,iterω]-ω[:,:,:,iterω+1]))<tolω)
# iterω = iterω+1
# println("__________________________")
# Save grids, expected aggregate human capital in teaching, initial guesses for 's' and 'e' (in each occupation):
# filename = string("results_",n_g,"_groups_tauW=",τ_w,"_tauE=",τ_e,".jld")
filename = string("/Users/simeonalder/Dropbox/Work/Research/GitHub/teachers/julia/results/results_",n_g,"_groups_tauW=",τ_w,"_tauE=",τ_e,"_beta-to-sigma=",round(β/σ,digits=1),"_A=",round(A,digits=1),".jld")

save(filename,"H_grid",H_grid,"a_grid",a_grid,"τ_e",τ_e,"τ_w",τ_w,"β",β,"σ",σ,"A",A,"HH_T",HH_T,"H_O",H_O,"λ",λ,"a_T_thresh",a_T_thresh,"a_O_thresh",a_O_thresh,"e_T",e_T,"s_T",s_T,"e_O",e_O,"s_O",s_O,"t",t,"E_O",E_O,"E_T",E_T,"mass_O",mass_O,"mass_T",mass_T,"f_1",f_1,"f_2",f_2)
