β=.15
σ=.13
α=.75
α_tilde = α*β/σ
η=.75
η_tilde = η*β/σ
μ=1/2
ϕ=1/3
ϕ_tilde = ϕ*β/σ
M = 1 # population size
A = 2 # productivity in 'Other'
θ = 3 # skill dispersion (high θ = small dispersion)
a_grid=quantile.(Frechet(θ),.005:.015:1)
# Intitial occupational cutoff
a_bar_T_0 = a_grid
# Based on Yulia's indifference condition:
sO = μ*ϕ/(μ*ϕ+1-η)
sT = μ*ϕ/(μ*ϕ+β/σ-η)
X1 = sO^ϕ*(M/2)^(-σ)
X2 = (M/2)^σ/(A*sO^ϕ*η)
X3 = (sT^ϕ_tilde)*(M/2)^(-β)
X4 = (M/2)^β/(A*sT^ϕ_tilde*η_tilde)
X5 = X3^(1-η_tilde)*X4^(-η_tilde)*X1^η_tilde*X2^(η*η_tilde/(η-1))
X6 = (X4*X5/(X1*X2^(η/(η-1))))^(1/(η_tilde-1))
X7 = X6^(η_tilde/(η_tilde-1))*sT^(ϕ_tilde)*(M/2)^(-β)
X8 = log(1-sT)
X9 = A*X1*X2^(η/(η-1))/X5
X11 = X8 + μ*log(X7*X9 - X6)
X12 = sO^(ϕ/(1-η))*(M/2)^(σ/(η-1))*A^(1/(1-η))*η^(η/(1-η))*(1-η)
X13 = X6 + μ*log(X12) - X11
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
mom1, er1 = quadgk(aa -> f_O_pwl(aa)*aa.^(α/(1-η)),a_grid[1],a_grid[end])
mom2, er2 = quadgk(aa -> f_O_pwl(aa),a_grid[1],a_grid[end])
mom3, er3 = quadgk(aa -> f_T_pwl(aa)*aa.^(α_tilde/(1-η_tilde)),a_grid[1],a_grid[end])
F = mom1/(mom2*mom3)
for i 
a_bar_T_0[i] = (((1-sO)/(1-sT))*((1-η)/(1-η_tilde))^μ*((1-τ_O_w)/(1-τ_T_e))^μ*a_grid[i]^(α*η/(1-η))*F^(-η))^((1-η_tilde)/α_tilde)
