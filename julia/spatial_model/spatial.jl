using LinearAlgebra
using Statistics
using Distributions
using JLD
using DataFrames
using CSV

"""Stable log-sum-exp for a length-2 vector."""
@inline function logsumexp2(a::Float64, b::Float64)
	m = max(a, b)
	return m + log(exp(a - m) + exp(b - m))
end

"""Softmax probabilities for length-2 utilities with scale σ (Gumbel)."""
@inline function softmax2(u1::Float64, u2::Float64, σ::Float64)
	x1 = u1 / σ
	x2 = u2 / σ
	lse = logsumexp2(x1, x2)
	return (exp(x1 - lse), exp(x2 - lse))
end

"""Quantile-bin histogram nodes (midpoint-quantiles) and equal weights."""
function quantile_histogram(dist::UnivariateDistribution, n::Int)
	nodes = Vector{Float64}(undef, n)
	w = fill(1.0 / n, n)
	for i in 1:n
		q = (i - 0.5) / n
		nodes[i] = quantile(dist, q)
	end
	return nodes, w
end

"""Piecewise-constant CDF implied by histogram nodes/weights (nodes assumed sorted)."""
function hist_cdf(nodes::AbstractVector{<:Real}, w::AbstractVector{<:Real}, x::Real)
	s = 0.0
	@inbounds for i in eachindex(nodes)
		if nodes[i] <= x
			s += w[i]
		else
			break
		end
	end
	return min(max(s, 0.0), 1.0)
end

"""Compute conditional expectation E[X | node > thresh] using histogram."""
function hist_cond_gt(nodes::Vector{Float64}, w::Vector{Float64}, x::Vector{Float64}, thresh::Float64)
	num = 0.0
	den = 0.0
	@inbounds for i in eachindex(nodes)
		if nodes[i] > thresh
			wi = w[i]
			den += wi
			num += wi * x[i]
		end
	end
	if den <= 0
		return 0.0  # empty conditioning set; caller multiplies by wo ≈ 0
	end
	return num / den
end

"""Compute conditional expectation E[X | node > thresh] for 2-column matrix X[:,j]."""
function hist_cond_gt2(nodes::Vector{Float64}, w::Vector{Float64}, x::Matrix{Float64}, thresh::Float64)
	num1 = 0.0
	num2 = 0.0
	den = 0.0
	@inbounds for i in eachindex(nodes)
		if nodes[i] > thresh
			wi = w[i]
			den += wi
			num1 += wi * x[i, 1]
			num2 += wi * x[i, 2]
		end
	end
	if den <= 0
		return (0.0, 0.0)  # empty conditioning set; caller multiplies by wo ≈ 0
	end
	return (num1 / den, num2 / den)
end

"""Compute conditional expectation E[X | node > thresh] for 2-column matrix X with weights, also return conditional mean of scalar y."""
function hist_cond_gt2_and_scalar(nodes::Vector{Float64}, w::Vector{Float64}, x2::Matrix{Float64}, y::Vector{Float64}, thresh::Float64)
	num1 = 0.0
	num2 = 0.0
	numy = 0.0
	den = 0.0
	@inbounds for i in eachindex(nodes)
		if nodes[i] > thresh
			wi = w[i]
			den += wi
			num1 += wi * x2[i, 1]
			num2 += wi * x2[i, 2]
			numy += wi * y[i]
		end
	end
	if den <= 0
		return (0.0, 0.0, 0.0)  # empty conditioning set; caller multiplies by wo ≈ 0
	end
	return (num1 / den, num2 / den, numy / den)
end

"""Compute Σ_{node > thresh} w[i]*π[i,1]*x[i,1] and Σ w[i]*π[i,2]*x[i,2] as a joint weighted sum.
Avoids the product-of-conditional-means approximation: returns E_{ε_O}[π_{l'}*taxO_{l'} | ε_O>thresh]*P(ε_O>thresh)."""
function hist_weighted_sum_gt_joint(nodes::Vector{Float64}, w::Vector{Float64}, π::Matrix{Float64}, x::Matrix{Float64}, thresh::Float64)
	s1 = 0.0
	s2 = 0.0
	@inbounds for i in eachindex(nodes)
		if nodes[i] > thresh
			wi = w[i]
			s1 += wi * π[i, 1] * x[i, 1]
			s2 += wi * π[i, 2] * x[i, 2]
		end
	end
	return (s1, s2)
end

"""Tauchen discretization of AR(1): x' = ρ x + ε, ε ~ N(0, σ_ε^2). Returns grid x and transition matrix P."""
function tauchen(n::Int, ρ::Float64, σ_ε::Float64; m::Float64 = 3.0)
	@assert n >= 3
	σ_x = σ_ε / sqrt(1 - ρ^2)
	x_max = m * σ_x
	x_min = -x_max
	x = collect(range(x_min, x_max; length = n))
	step = x[2] - x[1]
	P = Matrix{Float64}(undef, n, n)
	ϕ = Normal(0.0, σ_ε)
	@inbounds for i in 1:n
		for j in 1:n
			if j == 1
				z = (x[1] - ρ * x[i] + step / 2) / σ_ε
				P[i, j] = cdf(Normal(), z)
			elseif j == n
				z = (x[n] - ρ * x[i] - step / 2) / σ_ε
				P[i, j] = 1.0 - cdf(Normal(), z)
			else
				z_u = (x[j] - ρ * x[i] + step / 2) / σ_ε
				z_l = (x[j] - ρ * x[i] - step / 2) / σ_ε
				P[i, j] = cdf(Normal(), z_u) - cdf(Normal(), z_l)
			end
		end
		s = sum(P[i, :])
		P[i, :] ./= s
	end
	return x, P
end

"""Stationary distribution of a Markov chain with transition matrix P (rows sum to 1)."""
function stationary_dist(P::Matrix{Float64}; tol::Float64 = 1e-14, maxit::Int = 50_000)
	n = size(P, 1)
	π = fill(1.0 / n, n)
	for _ in 1:maxit
		π_new = (π' * P)[:]  # row vector times P
		if maximum(abs.(π_new .- π)) < tol
			return π_new ./ sum(π_new)
		end
		π = π_new
	end
	return π ./ sum(π)
end

struct Params
	L::Int
	genders::Vector{Symbol}
	α::Float64
	β::Float64
	η::Float64
	σ::Float64
	μ::Float64
	ϕ::Float64
	γ::Float64
	κ::Vector{Float64}              # κ_l for teaching wages
	A_O::Vector{Float64}            # non-teaching productivities by occupation
	τw_O::Dict{Symbol, Vector{Float64}} # labor wedges in O by gender and occupation
	τe_O::Dict{Symbol, Vector{Float64}} # education barriers in O by gender and occupation
	τw_T::Dict{Symbol, Float64}         # labor wedge for teaching by gender
	τe_T::Dict{Symbol, Float64}         # education barrier for teaching by gender
	ε_parent::Float64               # parent finance share
	λ::Float64                      # altruism
	B::Vector{Float64}              # amenities
	τmove::Matrix{Float64}          # moving costs in utils
	σν::Float64                     # Gumbel scale
	ρz::Float64                     # persistence in log z
	σξ::Float64                     # innovation sd in log z
	σeps::Float64                   # dispersion of lognormal ε
	Neps::Int
	Nz::Int
	damp_m::Float64
	damp_H::Float64
	damp_t::Float64
	damp_E::Float64
	tol::Float64
	maxit::Int
end

@inline n_nonteach_occupations(p::Params) = length(p.A_O)

function validate_occupation_primitives(p::Params)
	n_occ = n_nonteach_occupations(p)
	for gsym in p.genders
		@assert haskey(p.τw_O, gsym) "Missing τw_O entry for gender $(gsym)"
		@assert haskey(p.τe_O, gsym) "Missing τe_O entry for gender $(gsym)"
		@assert length(p.τw_O[gsym]) == n_occ "τw_O[$(gsym)] must have length $(n_occ)"
		@assert length(p.τe_O[gsym]) == n_occ "τe_O[$(gsym)] must have length $(n_occ)"
		@assert haskey(p.τw_T, gsym) "Missing τw_T entry for gender $(gsym)"
		@assert haskey(p.τe_T, gsym) "Missing τe_T entry for gender $(gsym)"
	end
	return nothing
end

struct Grids
	xz::Vector{Float64}        # log z grid
	z::Vector{Float64}         # z grid (positive)
	Pz::Matrix{Float64}
	eps_nodes::Vector{Float64} # common nodes for ε
	eps_w::Vector{Float64}
end

struct Eqm
	m_young::Matrix{Float64}    # [l,k] young distribution (sums to 1)
	M::Vector{Float64}          # total population in each location (young+old)
	Htilde::Vector{Float64}
	t::Vector{Float64}
	E_logh::Array{Float64,3}    # [k,l,g]
	E_cost::Array{Float64,3}    # [k,l,g]
	UT::Array{Float64,4}        # [k,l,g,eps] utility for teaching
	UO_best::Array{Float64,4}   # [k,l,g,eps] utility for best non-teaching
	costT::Array{Float64,4}     # [k,l,g,eps] goods cost if teaching
	costO_best::Array{Float64,4}# [k,l,g,eps] goods cost if best non-teaching
end

function avg_eqm(a::Eqm, b::Eqm)
	m = 0.5 .* (a.m_young .+ b.m_young)
	m ./= sum(m)
	return Eqm(
		m,
		0.5 .* (a.M .+ b.M),
		0.5 .* (a.Htilde .+ b.Htilde),
		0.5 .* (a.t .+ b.t),
		0.5 .* (a.E_logh .+ b.E_logh),
		0.5 .* (a.E_cost .+ b.E_cost),
		0.5 .* (a.UT .+ b.UT),
		0.5 .* (a.UO_best .+ b.UO_best),
		0.5 .* (a.costT .+ b.costT),
		0.5 .* (a.costO_best .+ b.costO_best),
	)
end

@inline function s_O_const(p::Params)
	return p.μ * p.ϕ / (p.μ * p.ϕ + 1 - p.η)
end

@inline function s_T_const(p::Params)
	return p.μ * p.ϕ * p.γ / ((p.μ * p.ϕ - p.η) * p.γ + 1)
end

"""Compute teacher-quality shifter Q_l = (2*Htilde_l/M_l)^σ."""
@inline function Q_l(p::Params, Htilde_l::Float64, M_l::Float64)
	return (2.0 * Htilde_l / M_l)^p.σ
end

@inline function safe_logC(C::Float64)
	return C > 1e-12 ? log(C) : -1e12
end

@inline function safe_invC(C::Float64)
	return C > 1e-12 ? 1.0 / C : 1e12
end

"""Solve for s in 1/(1-s) = K * s^(expo-1) via bisection."""
function solve_s_bisection(K::Float64, expo::Float64; tol::Float64 = 1e-10)
	if K <= 0
		return 1e-8
	end
	lo = 1e-8
	hi = 1.0 - 1e-8
	f_lo = 1.0 / (1.0 - lo) - K * lo^(expo - 1.0)
	f_hi = 1.0 / (1.0 - hi) - K * hi^(expo - 1.0)
	if f_lo >= 0
		return lo
	end
	if f_hi <= 0
		return hi
	end
	for _ in 1:80
		mid = 0.5 * (lo + hi)
		f_mid = 1.0 / (1.0 - mid) - K * mid^(expo - 1.0)
		if f_mid > 0
			hi = mid
		else
			lo = mid
		end
		if (hi - lo) < tol
			break
		end
	end
	return 0.5 * (lo + hi)
end

"""Expected log C and 1/C for a given child state (k,l,g), conditional on A (pre-child-cost resources)."""
function expected_logC_invC(
	p::Params,
	grids::Grids,
	UT::AbstractVector{Float64},
	UO_best::AbstractVector{Float64},
	costT::AbstractVector{Float64},
	costO_best::AbstractVector{Float64},
	A::Float64,
)
	Ne = length(UT)
	eps = grids.eps_nodes
	w = grids.eps_w
	log_sum = 0.0
	inv_sum = 0.0
	@inbounds for i in 1:Ne
		thr = invert_threshold(eps, UO_best, UT[i])
		pT = hist_cdf(eps, w, thr)
		Ct = A - p.ε_parent * costT[i]
		logT = safe_logC(Ct)
		invT = safe_invC(Ct)
		logO = 0.0
		invO = 0.0
		@inbounds for j in 1:Ne
			if eps[j] > thr
				Cj = A - p.ε_parent * costO_best[j]
				logO += w[j] * safe_logC(Cj)
				invO += w[j] * safe_invC(Cj)
			end
		end
		log_sum += w[i] * (pT * logT + logO)
		inv_sum += w[i] * (pT * invT + invO)
	end
	return log_sum, inv_sum
end

"""Expected log C and 1/C for parent with z index k and destination ldest."""
function expected_child_terms(
	p::Params,
	grids::Grids,
	eqm::Eqm,
	k_parent::Int,
	ldest::Int,
	A::Float64,
)
	log_sum = 0.0
	inv_sum = 0.0
	for kp in 1:p.Nz
		wkp = grids.Pz[k_parent, kp]
		for ig in 1:length(p.genders)
			UT = view(eqm.UT, kp, ldest, ig, :)
			UO = view(eqm.UO_best, kp, ldest, ig, :)
			cT = view(eqm.costT, kp, ldest, ig, :)
			cO = view(eqm.costO_best, kp, ldest, ig, :)
			logc, invc = expected_logC_invC(p, grids, UT, UO, cT, cO, A)
			log_sum += 0.5 * wkp * logc
			inv_sum += 0.5 * wkp * invc
		end
	end
	return log_sum, inv_sum
end

"""Expected child cost (goods) for parent with z index k and destination ldest."""
function expected_child_cost(p::Params, grids::Grids, eqm::Eqm, k_parent::Int, ldest::Int)
	c = 0.0
	for kp in 1:p.Nz
		wkp = grids.Pz[k_parent, kp]
		c += wkp * 0.5 * (eqm.E_cost[kp, ldest, 1] + eqm.E_cost[kp, ldest, 2])
	end
	return c
end

"""Given a birth location l, z, eps, solve location-choice and (s,e) fixed point for occupation T."""
function solve_T(p::Params, grids::Grids, eqm::Eqm, child_logh::Matrix{Float64}, l_birth::Int, k::Int, gsym::Symbol, epsT::Float64)
	π1 = 0.5
	sT = s_T_const(p)
	eT = 0.5
	Q = Q_l(p, eqm.Htilde[l_birth], eqm.M[l_birth])
	aT = grids.z[k] * epsT

	τw = p.τw_T[gsym]
	τe = p.τe_T[gsym]

	for _ in 1:200
		π = (π1, 1.0 - π1)
		h_base = Q * aT^p.α * sT^p.ϕ
		h = h_base * eT^p.η

		# destination-specific wages
		w1 = p.κ[1] * h^p.γ
		w2 = p.κ[2] * h^p.γ

		A1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w1 - (1.0 - p.ε_parent) * (1.0 + τe) * eT
		A2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w2 - (1.0 - p.ε_parent) * (1.0 + τe) * eT

		logC1, invC1 = expected_child_terms(p, grids, eqm, k, 1, A1)
		logC2, invC2 = expected_child_terms(p, grids, eqm, k, 2, A2)

		u1 = p.μ * logC1 + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
		u2 = p.μ * logC2 + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]

		π1_new, _ = softmax2(u1, u2, p.σν)

		S0 = π[1] * invC1 + π[2] * invC2
		Sκ = π[1] * (1.0 - eqm.t[1]) * p.κ[1] * invC1 + π[2] * (1.0 - eqm.t[2]) * p.κ[2] * invC2
		Sκ = max(Sκ, 1e-12)
		S0 = max(S0, 1e-12)

		# update e via teaching FOC
		denom = (1.0 - τw) * p.γ * p.η * h_base^p.γ * Sκ
		if denom > 0
			rhs = (1.0 - p.ε_parent) * (1.0 + τe) * S0 / denom
			eT_new = rhs^(1.0 / (p.η * p.γ - 1.0))
			eT_new = max(eT_new, 1e-8)
		else
			eT_new = 1e-8
		end

		# update s via teaching FOC (K excludes s to isolate s^(ϕγ-1))
		base_no_s = Q * aT^p.α
		Ks = p.μ * p.ϕ * p.γ * (1.0 - τw) * base_no_s^p.γ * eT^(p.η * p.γ) * Sκ
		sT_new = solve_s_bisection(Ks, p.ϕ * p.γ)

		s_err = abs(sT_new - sT)
		e_err = abs(eT_new - eT)
		π_err = abs(π1_new - π1)

		# damped updates
		sT = 0.6 * sT + 0.4 * sT_new
		eT = 0.6 * eT + 0.4 * eT_new
		π1 = 0.6 * π1 + 0.4 * π1_new

		if max(s_err, e_err, π_err) < 1e-8
			inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
			πmat = (π1, 1.0 - π1)
			C1 = A1 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 1)
			C2 = A2 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 2)
			return (inc, πmat, h, eT, (w1, w2), (C1, C2), sT)
		end
	end
	h_base = Q * aT^p.α * sT^p.ϕ
	h = h_base * eT^p.η
	w1 = p.κ[1] * h^p.γ
	w2 = p.κ[2] * h^p.γ
	A1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w1 - (1.0 - p.ε_parent) * (1.0 + τe) * eT
	A2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w2 - (1.0 - p.ε_parent) * (1.0 + τe) * eT
	logC1, _ = expected_child_terms(p, grids, eqm, k, 1, A1)
	logC2, _ = expected_child_terms(p, grids, eqm, k, 2, A2)
	u1 = p.μ * logC1 + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
	u2 = p.μ * logC2 + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]
	inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
	πmat = (π1, 1.0 - π1)
	C1 = A1 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 1)
	C2 = A2 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 2)
	return (inc, πmat, h, eT, (w1, w2), (C1, C2), sT)
end

"""Solve location choice fixed point for occupation O (non-teaching)."""
function solve_O(
	p::Params,
	grids::Grids,
	eqm::Eqm,
	child_logh::Matrix{Float64},
	l_birth::Int,
	k::Int,
	gsym::Symbol,
	epsO::Float64,
	occ_idx::Int,
)
	π1 = 0.5
	sO = s_O_const(p)
	eO = 0.5
	Q = Q_l(p, eqm.Htilde[l_birth], eqm.M[l_birth])
	aO = grids.z[k] * epsO
	A_occ = p.A_O[occ_idx]
	τw = p.τw_O[gsym][occ_idx]
	τe = p.τe_O[gsym][occ_idx]

	for _ in 1:200
		π = (π1, 1.0 - π1)
		h_base = Q * aO^p.α * sO^p.ϕ
		h = h_base * eO^p.η
		w = A_occ * h

		A1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * eO
		A2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * eO

		logC1, invC1 = expected_child_terms(p, grids, eqm, k, 1, A1)
		logC2, invC2 = expected_child_terms(p, grids, eqm, k, 2, A2)

		u1 = p.μ * logC1 + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
		u2 = p.μ * logC2 + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]
		π1_new, _ = softmax2(u1, u2, p.σν)

		S0 = π[1] * invC1 + π[2] * invC2
		S1 = π[1] * (1.0 - eqm.t[1]) * invC1 + π[2] * (1.0 - eqm.t[2]) * invC2
		S1 = max(S1, 1e-12)
		S0 = max(S0, 1e-12)

		denom = (1.0 - τw) * A_occ * p.η * h_base * S1
		if denom > 0
			rhs = (1.0 - p.ε_parent) * (1.0 + τe) * S0 / denom
			eO_new = rhs^(1.0 / (p.η - 1.0))
			eO_new = max(eO_new, 1e-8)
		else
			eO_new = 1e-8
		end

		base_no_s = Q * aO^p.α
		Ks = p.μ * p.ϕ * (1.0 - τw) * A_occ * base_no_s * eO^p.η * S1
		sO_new = solve_s_bisection(Ks, p.ϕ)

		s_err = abs(sO_new - sO)
		e_err = abs(eO_new - eO)
		π_err = abs(π1_new - π1)

		sO = 0.6 * sO + 0.4 * sO_new
		eO = 0.6 * eO + 0.4 * eO_new
		π1 = 0.6 * π1 + 0.4 * π1_new

		if max(s_err, e_err, π_err) < 1e-8
			inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
			πmat = (π1, 1.0 - π1)
			C1 = A1 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 1)
			C2 = A2 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 2)
			return (inc, πmat, h, eO, w, (C1, C2), sO)
		end
	end
	h_base = Q * aO^p.α * sO^p.ϕ
	h = h_base * eO^p.η
	w = A_occ * h
	A1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * eO
	A2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * eO
	logC1, _ = expected_child_terms(p, grids, eqm, k, 1, A1)
	logC2, _ = expected_child_terms(p, grids, eqm, k, 2, A2)
	u1 = p.μ * logC1 + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
	u2 = p.μ * logC2 + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]
	inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
	πmat = (π1, 1.0 - π1)
	C1 = A1 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 1)
	C2 = A2 - p.ε_parent * expected_child_cost(p, grids, eqm, k, 2)
	return (inc, πmat, h, eO, w, (C1, C2), sO)
end

"""Invert U_O(epsO) to find threshold epsO such that U_O = U_T, assuming monotone U_O."""
function invert_threshold(epsO_nodes::AbstractVector{Float64}, UO::AbstractVector{Float64}, UT::Float64)
	n = length(epsO_nodes)
	if UT <= UO[1]
		return epsO_nodes[1] - 1e-12
	end
	if UT >= UO[n]
		return epsO_nodes[n] + 1e-12
	end
	# find j such that UO[j] <= UT < UO[j+1]
	j = 1
	@inbounds for i in 1:n-1
		if (UO[i] <= UT) && (UT < UO[i + 1])
			j = i
			break
		end
	end
	# linear interpolation in utility space
	u0 = UO[j]
	u1 = UO[j + 1]
	x0 = epsO_nodes[j]
	x1 = epsO_nodes[j + 1]
	t = (UT - u0) / (u1 - u0)
	return x0 + t * (x1 - x0)
end

"""Compute deterministic moments for a given (birth location l, z-state k, gender g)."""
function state_moments(p::Params, grids::Grids, eqm::Eqm, child_logh::Matrix{Float64}, l_birth::Int, k::Int, gsym::Symbol)
	Ne = p.Neps
	eps = grids.eps_nodes
	w = grids.eps_w

	# Solve T side over eps_T nodes
	UT = Vector{Float64}(undef, Ne)
	loghT = Vector{Float64}(undef, Ne)
	costT = Vector{Float64}(undef, Ne)
	πT = Matrix{Float64}(undef, Ne, 2)
	wT = Matrix{Float64}(undef, Ne, 2)
	taxableT = Matrix{Float64}(undef, Ne, 2)
	billT = Matrix{Float64}(undef, Ne, 2)
	hpowT = Vector{Float64}(undef, Ne) # h^{β/σ}

	@inbounds for i in 1:Ne
		inc, π, h, e, (w1, w2), _, sT = solve_T(p, grids, eqm, child_logh, l_birth, k, gsym, eps[i])
		UT[i] = log(max(1.0 - sT, 1e-12)) + inc
		loghT[i] = log(max(h, 1e-300))
		costT[i] = (1.0 + p.τe_T[gsym]) * e
		πT[i, 1] = π[1]
		πT[i, 2] = π[2]
		wT[i, 1] = w1
		wT[i, 2] = w2
		taxableT[i, 1] = (1.0 - p.τw_T[gsym]) * w1
		taxableT[i, 2] = (1.0 - p.τw_T[gsym]) * w2
		billT[i, 1] = w1
		billT[i, 2] = w2
		hpowT[i] = h^(p.β / p.σ)
	end

	# Solve O side over eps_O nodes and all non-teaching occupations.
	n_occ = n_nonteach_occupations(p)
	UO = Matrix{Float64}(undef, Ne, n_occ)
	loghO = Matrix{Float64}(undef, Ne, n_occ)
	costO = Matrix{Float64}(undef, Ne, n_occ)
	πO = Array{Float64}(undef, Ne, n_occ, 2)
	taxableO = Array{Float64}(undef, Ne, n_occ, 2)

	@inbounds for iocc in 1:n_occ
		τw_occ = p.τw_O[gsym][iocc]
		τe_occ = p.τe_O[gsym][iocc]
		for i in 1:Ne
			inc, π, h, e, wgross, _, sO = solve_O(p, grids, eqm, child_logh, l_birth, k, gsym, eps[i], iocc)
			UO[i, iocc] = log(max(1.0 - sO, 1e-12)) + inc
			loghO[i, iocc] = log(max(h, 1e-300))
			costO[i, iocc] = (1.0 + τe_occ) * e
			πO[i, iocc, 1] = π[1]
			πO[i, iocc, 2] = π[2]
			taxableO[i, iocc, 1] = (1.0 - τw_occ) * wgross
			taxableO[i, iocc, 2] = (1.0 - τw_occ) * wgross
		end
	end

	# Collapse to best non-teaching option at each eps node.
	UO_best = Vector{Float64}(undef, Ne)
	loghO_best = Vector{Float64}(undef, Ne)
	costO_best = Vector{Float64}(undef, Ne)
	πO_best = Matrix{Float64}(undef, Ne, 2)
	taxableO_best = Matrix{Float64}(undef, Ne, 2)
	@inbounds for i in 1:Ne
		jstar = argmax(@view UO[i, :])
		UO_best[i] = UO[i, jstar]
		loghO_best[i] = loghO[i, jstar]
		costO_best[i] = costO[i, jstar]
		πO_best[i, 1] = πO[i, jstar, 1]
		πO_best[i, 2] = πO[i, jstar, 2]
		taxableO_best[i, 1] = taxableO[i, jstar, 1]
		taxableO_best[i, 2] = taxableO[i, jstar, 2]
	end

	# Ensure monotonicity for inversion (small numerical noise may violate).
	@inbounds for i in 2:Ne
		if UO_best[i] < UO_best[i - 1]
			UO_best[i] = UO_best[i - 1]
			loghO_best[i] = loghO_best[i - 1]
			costO_best[i] = costO_best[i - 1]
			πO_best[i, 1] = πO_best[i - 1, 1]
			πO_best[i, 2] = πO_best[i - 1, 2]
			taxableO_best[i, 1] = taxableO_best[i - 1, 1]
			taxableO_best[i, 2] = taxableO_best[i - 1, 2]
		end
	end

	# Aggregate moments with threshold method
	probT = 0.0
	E_logh = 0.0
	E_cost = 0.0
	mig = (0.0, 0.0)
	Hcontrib = (0.0, 0.0)
	taxbase = (0.0, 0.0)
	wagebill = (0.0, 0.0)

	@inbounds for i in 1:Ne
		thr = invert_threshold(eps, UO_best, UT[i])
		pT = hist_cdf(eps, w, thr)
		probT += w[i] * pT

		# teaching part
		wt = w[i] * pT
		E_logh += wt * loghT[i]
		E_cost += wt * costT[i]
		mig = (mig[1] + wt * πT[i, 1], mig[2] + wt * πT[i, 2])
		Hcontrib = (Hcontrib[1] + wt * πT[i, 1] * hpowT[i], Hcontrib[2] + wt * πT[i, 2] * hpowT[i])
		taxbase = (taxbase[1] + wt * πT[i, 1] * taxableT[i, 1], taxbase[2] + wt * πT[i, 2] * taxableT[i, 2])
		wagebill = (wagebill[1] + wt * πT[i, 1] * billT[i, 1], wagebill[2] + wt * πT[i, 2] * billT[i, 2])

		# non-teaching part (conditional epsO > thr)
		wo = w[i] * (1.0 - pT)
		if wo > 0
			loghO_gt = hist_cond_gt(eps, w, loghO_best, thr)
			costO_gt = hist_cond_gt(eps, w, costO_best, thr)
			πO1_gt, πO2_gt = hist_cond_gt2(eps, w, πO_best, thr)
			# joint weighted sum Σ_{j>thr} w[j]*π[j,l']*taxO[j,l'] — avoids product-of-means bias
			tax_joint1, tax_joint2 = hist_weighted_sum_gt_joint(eps, w, πO_best, taxableO_best, thr)
			E_logh += wo * loghO_gt
			E_cost += wo * costO_gt
			mig = (mig[1] + wo * πO1_gt, mig[2] + wo * πO2_gt)
			taxbase = (taxbase[1] + w[i] * tax_joint1, taxbase[2] + w[i] * tax_joint2)
		end
	end

	return (probT, E_logh, E_cost, mig, Hcontrib, taxbase, wagebill, UT, UO_best, costT, costO_best)
end

function solve_stationary(p::Params)
	validate_occupation_primitives(p)
	xz, Pz = tauchen(p.Nz, p.ρz, p.σξ)
	zgrid = exp.(xz)
	eps_dist = LogNormal(-0.5 * p.σeps^2, p.σeps) # mean 1
	eps_nodes, eps_w = quantile_histogram(eps_dist, p.Neps)
	grids = Grids(xz, zgrid, Pz, eps_nodes, eps_w)

	# initial young distribution: symmetric across l, stationary in z
	πz = stationary_dist(Pz)
	m0 = Matrix{Float64}(undef, p.L, p.Nz)
	for l in 1:p.L
		m0[l, :] .= (1.0 / p.L) .* πz
	end
	m0 ./= sum(m0)
	M0 = [2.0 * sum(m0[1, :]), 2.0 * sum(m0[2, :])]
	H0 = fill(1.0, p.L)
	t0 = fill(0.2, p.L)
	Elogh0 = zeros(p.Nz, p.L, length(p.genders))
	Ecost0 = zeros(p.Nz, p.L, length(p.genders))
	UT0 = zeros(p.Nz, p.L, length(p.genders), p.Neps)
	UO0 = zeros(p.Nz, p.L, length(p.genders), p.Neps)
	costT0 = zeros(p.Nz, p.L, length(p.genders), p.Neps)
	costO0 = zeros(p.Nz, p.L, length(p.genders), p.Neps)
	eqm = Eqm(m0, M0, H0, t0, Elogh0, Ecost0, UT0, UO0, costT0, costO0)

	# main fixed point loop
	eqm_prev = eqm
	eqm_prev2 = eqm
	for it in 1:p.maxit
		# continuation objects for parents: child expected logh by (parent z k, destination l')
		child_logh = Matrix{Float64}(undef, p.Nz, p.L)
		for k in 1:p.Nz
			for ldest in 1:p.L
				# average across genders for child
				elog = 0.0
				for kp in 1:p.Nz
					wkp = grids.Pz[k, kp]
					elog_kp = 0.5 * (eqm.E_logh[kp, ldest, 1] + eqm.E_logh[kp, ldest, 2])
					elog += wkp * elog_kp
				end
				child_logh[k, ldest] = elog
			end
		end

		# compute new E objects and aggregation moments
		Elogh_new = similar(eqm.E_logh)
		Ecost_new = similar(eqm.E_cost)
		UT_new = similar(eqm.UT)
		UO_new = similar(eqm.UO_best)
		costT_new = similar(eqm.costT)
		costO_new = similar(eqm.costO_best)
		barπ = Array{Float64,3}(undef, p.L, p.Nz, p.L)       # [l_birth,k,dest]
		mig_prob = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L) # [l,k,g,dest]
		Hcontrib = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L) # [l,k,g,dest]
		taxbase = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L)
		wagebill = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L)

		for l in 1:p.L
			for k in 1:p.Nz
				# gender-specific moments
				for (ig, gsym) in enumerate(p.genders)
					probT, elog, eco, mig, Hc, tb, wb, UT, UO_best, costT, costO_best = state_moments(p, grids, eqm, child_logh, l, k, gsym)
					Elogh_new[k, l, ig] = elog
					Ecost_new[k, l, ig] = eco
					UT_new[k, l, ig, :] = UT
					UO_new[k, l, ig, :] = UO_best
					costT_new[k, l, ig, :] = costT
					costO_new[k, l, ig, :] = costO_best
					mig_prob[l, k, ig, 1] = mig[1]
					mig_prob[l, k, ig, 2] = mig[2]
					Hcontrib[l, k, ig, 1] = Hc[1]
					Hcontrib[l, k, ig, 2] = Hc[2]
					taxbase[l, k, ig, 1] = tb[1]
					taxbase[l, k, ig, 2] = tb[2]
					wagebill[l, k, ig, 1] = wb[1]
					wagebill[l, k, ig, 2] = wb[2]
					# migration probabilities by gender used for population (young are 50/50)
					# mig already integrated over occ & eps for this gender
				end
				# average migration across genders for population law
				mig1 = 0.5 * (mig_prob[l, k, 1, 1] + mig_prob[l, k, 2, 1])
				mig2 = 0.5 * (mig_prob[l, k, 1, 2] + mig_prob[l, k, 2, 2])
				s = mig1 + mig2
				if s <= 0
					barπ[l, k, 1] = 0.5
					barπ[l, k, 2] = 0.5
				else
					barπ[l, k, 1] = mig1 / s
					barπ[l, k, 2] = mig2 / s
				end
			end
		end

		# update young distribution: migrate then reproduce
		old_mass = zeros(p.L, p.Nz)
		for l in 1:p.L
			for k in 1:p.Nz
				mlk = eqm.m_young[l, k]
				old_mass[1, k] += mlk * barπ[l, k, 1]
				old_mass[2, k] += mlk * barπ[l, k, 2]
			end
		end
		m_new = zeros(p.L, p.Nz)
		for ldest in 1:p.L
			for k in 1:p.Nz
				om = old_mass[ldest, k]
				if om == 0
					continue
				end
				@inbounds for kp in 1:p.Nz
					m_new[ldest, kp] += om * grids.Pz[k, kp]
				end
			end
		end
		m_new ./= sum(m_new)
		M_new = [2.0 * sum(m_new[1, :]), 2.0 * sum(m_new[2, :])]

		# update teacher stock Htilde and taxes
		H_new = zeros(p.L)
		I_new = zeros(p.L)
		W_new = zeros(p.L)
		for l0 in 1:p.L
			for k in 1:p.Nz
				mass = eqm.m_young[l0, k]
				for ig in 1:length(p.genders)
					wgt = 0.5 * mass
					for ldest in 1:p.L
						H_new[ldest] += wgt * Hcontrib[l0, k, ig, ldest]
						I_new[ldest] += wgt * taxbase[l0, k, ig, ldest]
						W_new[ldest] += wgt * wagebill[l0, k, ig, ldest]
					end
				end
			end
		end
		t_new = similar(eqm.t)
		for l in 1:p.L
			denom = max(I_new[l], 1e-12)
			t_new[l] = min(max(W_new[l] / denom, 0.0), 0.95)
		end

		# damped updates
		m_upd = (1 - p.damp_m) .* eqm.m_young .+ p.damp_m .* m_new
		m_upd ./= sum(m_upd)
		H_upd = (1 - p.damp_H) .* eqm.Htilde .+ p.damp_H .* H_new
		t_upd = (1 - p.damp_t) .* eqm.t .+ p.damp_t .* t_new
		Elogh_upd = (1 - p.damp_E) .* eqm.E_logh .+ p.damp_E .* Elogh_new
		Ecost_upd = (1 - p.damp_E) .* eqm.E_cost .+ p.damp_E .* Ecost_new
		M_upd = (1 - p.damp_m) .* eqm.M .+ p.damp_m .* M_new

		UT_upd = UT_new
		UO_upd = UO_new
		costT_upd = costT_new
		costO_upd = costO_new

		# convergence check
		err_m = maximum(abs.(m_upd .- eqm.m_young))
		err_H = maximum(abs.(H_upd .- eqm.Htilde) ./ max.(abs.(eqm.Htilde), 1e-8))
		err_t = maximum(abs.(t_upd .- eqm.t))
		err = max(err_m, err_H, err_t)

		eqm_new = Eqm(m_upd, M_upd, H_upd, t_upd, Elogh_upd, Ecost_upd, UT_upd, UO_upd, costT_upd, costO_upd)

		if it % 25 == 0
			println("iter=", it, " err=", err, " M=", eqm_new.M, " t=", eqm_new.t, " H=", eqm_new.Htilde)
		end
		if err < p.tol
			println("Converged in ", it, " iterations with err=", err)
			return eqm_new, grids
		end

		# Detect a small deterministic 2-cycle (common with threshold/histogram maps)
		if it > 2
			err2_m = maximum(abs.(eqm_new.m_young .- eqm_prev2.m_young))
			err2_H = maximum(abs.(eqm_new.Htilde .- eqm_prev2.Htilde) ./ max.(abs.(eqm_prev2.Htilde), 1e-8))
			err2_t = maximum(abs.(eqm_new.t .- eqm_prev2.t))
			err2 = max(err2_m, err2_H, err2_t)
			if err2 < p.tol
				println("Converged to period-2 cycle in ", it, " iterations (err2=", err2, "); returning averaged equilibrium.")
				return avg_eqm(eqm_new, eqm_prev), grids
			end
		end

		eqm_prev2 = eqm_prev
		eqm_prev = eqm_new
		eqm = eqm_new
	end
	println("WARNING: did not converge in maxit")
	return eqm, grids
end

function load_baseline(paramname::String, year::Int)
	path = joinpath(@__DIR__, "..", "codes_w_exo_wage", "parameterization", paramname, "previousParameterization$(year).jld")
	d = load(path)
	return d
end

function main()
	# ---- user-facing knobs (v1) ----
	paramname = "benchmark"
	year = 1970

	# Spatial primitives
	B = [0.0, 0.0]
	τmove = [0.0 0.3; 0.3 0.0]
	σν = 0.2
	ε_parent = 0.2
	λ = 0.1

	# Ability process (log z) and idiosyncratic occupation shocks
	ρz = 0.7
	σξ = 0.4
	σeps = 0.6
	Nz = 7
	Neps = 21

	# Fixed-point controls
	damp_m = 0.15
	damp_H = 0.15
	damp_t = 0.15
	damp_E = 0.30
	tol = 1e-5
	maxit = 2000

	# ---- load calibrated primitives ----
	d = load_baseline(paramname, year)
	α = Float64(d["α"])
	β = Float64(d["β"])
	η = Float64(d["η"])
	σ = Float64(d["σ"])
	μ = Float64(d["μ"])
	ϕ = Float64(d["ϕ"])
	γ = Float64(d["γ"])

	κ0 = Float64(d["κ"])
	κ = [κ0, 1.05 * κ0]  # location-specific scale

	a_by_occ = Vector{Float64}(d["a_by_occ"])
	A_O = copy(a_by_occ)

	τ_w = Array{Float64}(d["τ_w"])
	τ_e = Array{Float64}(d["τ_e"])
	# τ_* dimensions: (nOccNonTeach, nG), kept at full occupation detail.
	τw_O = Dict(:f => Vector{Float64}(τ_w[:, 1]), :m => Vector{Float64}(τ_w[:, 2]))
	τe_O = Dict(:f => Vector{Float64}(τ_e[:, 1]), :m => Vector{Float64}(τ_e[:, 2]))
	# Teaching distortions: set to zero unless the calibration provides them.
	τw_T = Dict(:f => 0.0, :m => 0.0)
	τe_T = Dict(:f => 0.0, :m => 0.0)

	p = Params(
		2,
		[:f, :m],
		α, β, η, σ, μ, ϕ, γ,
		κ,
		A_O,
		τw_O,
		τe_O,
		τw_T,
		τe_T,
		ε_parent,
		λ,
		B,
		τmove,
		σν,
		ρz,
		σξ,
		σeps,
		Neps,
		Nz,
		damp_m,
		damp_H,
		damp_t,
		damp_E,
		tol,
		maxit,
	)

	eqm, grids = solve_stationary(p)

	# ---- outputs ----
	outdir = @__DIR__
	jld_path = joinpath(outdir, "spatial_eqm_$(paramname)_$(year).jld")
	save(jld_path,
		"Params", p,
		"m_young", eqm.m_young,
		"M", eqm.M,
		"Htilde", eqm.Htilde,
		"t", eqm.t,
		"zgrid", grids.z,
		"Pz", grids.Pz,
		"eps_nodes", grids.eps_nodes,
	)

	df = DataFrame(location = 1:p.L, M = eqm.M, Htilde = eqm.Htilde, t = eqm.t, B = p.B)
	csv_path = joinpath(outdir, "spatial_eqm_summary_$(paramname)_$(year).csv")
	CSV.write(csv_path, df)
	println("Wrote ", jld_path)
	println("Wrote ", csv_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
