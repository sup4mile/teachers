
#!/usr/bin/env julia

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
		return NaN
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
		return (NaN, NaN)
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
		return (NaN, NaN, NaN)
	end
	return (num1 / den, num2 / den, numy / den)
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
	A_O::Float64                    # non-teaching composite productivity
	τw_O::Dict{Symbol, Float64}     # labor wedge in O by gender
	τe_O::Dict{Symbol, Float64}     # education barrier in O by gender
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

"""Given a birth location l, z, eps, and a scalar p1=Pr(dest=1), solve location choice fixed point for occupation T."""
function solve_T(p::Params, grids::Grids, eqm::Eqm, child_logh::Matrix{Float64}, child_cost::Matrix{Float64}, l_birth::Int, k::Int, gsym::Symbol, epsT::Float64)
	# scalar fixed point on π1
	π1 = 0.5
	sT = s_T_const(p)
	Q = Q_l(p, eqm.Htilde[l_birth], eqm.M[l_birth])
	aT = grids.z[k] * epsT

	τw = 0.0
	τe = 0.0

	for _ in 1:200
		π = (π1, 1.0 - π1)
		κeff = π[1] * (1.0 - eqm.t[1]) * p.κ[1] + π[2] * (1.0 - eqm.t[2]) * p.κ[2]
		κeff = max(κeff, 1e-12)
		# closed form for power wage profile
		h = (p.η^p.η * (κeff * p.γ)^p.η * aT^p.α * sT^p.ϕ * Q)^(1.0 / (1.0 - p.η * p.γ))
		e = p.η * (κeff * p.γ) * h^p.γ
		# destination-specific wages
		w1 = p.κ[1] * h^p.γ
		w2 = p.κ[2] * h^p.γ
		# consumption per destination
		C1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w1 - (1.0 - p.ε_parent) * (1.0 + τe) * e - p.ε_parent * child_cost[k, 1]
		C2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w2 - (1.0 - p.ε_parent) * (1.0 + τe) * e - p.ε_parent * child_cost[k, 2]
		u1 = (C1 > 1e-12 ? p.μ * log(C1) : -1e12) + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
		u2 = (C2 > 1e-12 ? p.μ * log(C2) : -1e12) + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]
		π1_new, _ = softmax2(u1, u2, p.σν)
		if abs(π1_new - π1) < 1e-10
			inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
			πmat = (π1_new, 1.0 - π1_new)
			return (inc, πmat, h, e, (w1, w2), (C1, C2))
		end
		π1 = 0.6 * π1 + 0.4 * π1_new
	end
	inc = -1e12
	πmat = (π1, 1.0 - π1)
	return (inc, πmat, NaN, NaN, (NaN, NaN), (NaN, NaN))
end

"""Solve location choice fixed point for occupation O (non-teaching)."""
function solve_O(p::Params, grids::Grids, eqm::Eqm, child_logh::Matrix{Float64}, child_cost::Matrix{Float64}, l_birth::Int, k::Int, gsym::Symbol, epsO::Float64)
	π1 = 0.5
	sO = s_O_const(p)
	Q = Q_l(p, eqm.Htilde[l_birth], eqm.M[l_birth])
	aO = grids.z[k] * epsO
	τw = p.τw_O[gsym]
	τe = p.τe_O[gsym]

	for _ in 1:200
		π = (π1, 1.0 - π1)
		taxeff = π[1] * (1.0 - eqm.t[1]) + π[2] * (1.0 - eqm.t[2])
		taxeff = max(taxeff, 1e-12)
		scale = taxeff * (1.0 - τw) / (1.0 + τe) * p.A_O
		h = (p.η^p.η * scale^p.η * aO^p.α * sO^p.ϕ * Q)^(1.0 / (1.0 - p.η))
		e = p.η * scale * h
		w = p.A_O * h
		C1 = (1.0 - eqm.t[1]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * e - p.ε_parent * child_cost[k, 1]
		C2 = (1.0 - eqm.t[2]) * (1.0 - τw) * w - (1.0 - p.ε_parent) * (1.0 + τe) * e - p.ε_parent * child_cost[k, 2]
		u1 = (C1 > 1e-12 ? p.μ * log(C1) : -1e12) + p.B[1] - p.τmove[l_birth, 1] + p.λ * child_logh[k, 1]
		u2 = (C2 > 1e-12 ? p.μ * log(C2) : -1e12) + p.B[2] - p.τmove[l_birth, 2] + p.λ * child_logh[k, 2]
		π1_new, _ = softmax2(u1, u2, p.σν)
		if abs(π1_new - π1) < 1e-10
			inc = p.σν * logsumexp2(u1 / p.σν, u2 / p.σν)
			πmat = (π1_new, 1.0 - π1_new)
			return (inc, πmat, h, e, w, (C1, C2))
		end
		π1 = 0.6 * π1 + 0.4 * π1_new
	end
	inc = -1e12
	πmat = (π1, 1.0 - π1)
	return (inc, πmat, NaN, NaN, NaN, (NaN, NaN))
end

"""Invert U_O(epsO) to find threshold epsO such that U_O = U_T, assuming monotone U_O."""
function invert_threshold(epsO_nodes::Vector{Float64}, UO::Vector{Float64}, UT::Float64)
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
function state_moments(p::Params, grids::Grids, eqm::Eqm, child_logh::Matrix{Float64}, child_cost::Matrix{Float64}, l_birth::Int, k::Int, gsym::Symbol)
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
		inc, π, h, e, (w1, w2), _ = solve_T(p, grids, eqm, child_logh, child_cost, l_birth, k, gsym, eps[i])
		UT[i] = log(1.0 - s_T_const(p)) + inc
		loghT[i] = log(max(h, 1e-300))
		costT[i] = (1.0 + 0.0) * e
		πT[i, 1] = π[1]
		πT[i, 2] = π[2]
		wT[i, 1] = w1
		wT[i, 2] = w2
		taxableT[i, 1] = (1.0 - 0.0) * w1
		taxableT[i, 2] = (1.0 - 0.0) * w2
		billT[i, 1] = w1
		billT[i, 2] = w2
		hpowT[i] = h^(p.β / p.σ)
	end

	# Solve O side over eps_O nodes
	UO = Vector{Float64}(undef, Ne)
	loghO = Vector{Float64}(undef, Ne)
	costO = Vector{Float64}(undef, Ne)
	πO = Matrix{Float64}(undef, Ne, 2)
	taxableO = Matrix{Float64}(undef, Ne, 2)

	τw = p.τw_O[gsym]
	@inbounds for i in 1:Ne
		inc, π, h, e, wgross, _ = solve_O(p, grids, eqm, child_logh, child_cost, l_birth, k, gsym, eps[i])
		UO[i] = log(1.0 - s_O_const(p)) + inc
		loghO[i] = log(max(h, 1e-300))
		costO[i] = (1.0 + p.τe_O[gsym]) * e
		πO[i, 1] = π[1]
		πO[i, 2] = π[2]
		taxableO[i, 1] = (1.0 - τw) * wgross
		taxableO[i, 2] = (1.0 - τw) * wgross
	end

	# Ensure monotonicity for inversion (small numerical noise may violate)
	# Sort by eps (already sorted). Smooth UO to be nondecreasing.
	@inbounds for i in 2:Ne
		if UO[i] < UO[i - 1]
			UO[i] = UO[i - 1]
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
		thr = invert_threshold(eps, UO, UT[i])
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
			loghO_gt = hist_cond_gt(eps, w, loghO, thr)
			costO_gt = hist_cond_gt(eps, w, costO, thr)
			πO1_gt, πO2_gt = hist_cond_gt2(eps, w, πO, thr)
			tax1_gt, tax2_gt = hist_cond_gt2(eps, w, taxableO, thr)
			E_logh += wo * loghO_gt
			E_cost += wo * costO_gt
			mig = (mig[1] + wo * πO1_gt, mig[2] + wo * πO2_gt)
			taxbase = (taxbase[1] + wo * πO1_gt * tax1_gt, taxbase[2] + wo * πO2_gt * tax2_gt)
		end
	end

	return (probT, E_logh, E_cost, mig, Hcontrib, taxbase, wagebill)
end

function solve_stationary(p::Params)
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
	eqm = Eqm(m0, M0, H0, t0, Elogh0, Ecost0)

	# main fixed point loop
	eqm_prev = eqm
	eqm_prev2 = eqm
	for it in 1:p.maxit
		# continuation objects for parents: child expected logh and cost by (parent z k, destination l')
		child_logh = Matrix{Float64}(undef, p.Nz, p.L)
		child_cost = Matrix{Float64}(undef, p.Nz, p.L)
		for k in 1:p.Nz
			for ldest in 1:p.L
				# average across genders for child
				elog = 0.0
				eco = 0.0
				for kp in 1:p.Nz
					wkp = grids.Pz[k, kp]
					elog_kp = 0.5 * (eqm.E_logh[kp, ldest, 1] + eqm.E_logh[kp, ldest, 2])
					eco_kp = 0.5 * (eqm.E_cost[kp, ldest, 1] + eqm.E_cost[kp, ldest, 2])
					elog += wkp * elog_kp
					eco += wkp * eco_kp
				end
				child_logh[k, ldest] = elog
				child_cost[k, ldest] = eco
			end
		end

		# compute new E objects and aggregation moments
		Elogh_new = similar(eqm.E_logh)
		Ecost_new = similar(eqm.E_cost)
		barπ = Array{Float64,3}(undef, p.L, p.Nz, p.L)       # [l_birth,k,dest]
		mig_prob = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L) # [l,k,g,dest]
		Hcontrib = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L) # [l,k,g,dest]
		taxbase = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L)
		wagebill = Array{Float64,4}(undef, p.L, p.Nz, length(p.genders), p.L)

		for l in 1:p.L
			for k in 1:p.Nz
				# gender-specific moments
				for (ig, gsym) in enumerate(p.genders)
					probT, elog, eco, mig, Hc, tb, wb = state_moments(p, grids, eqm, child_logh, child_cost, l, k, gsym)
					Elogh_new[k, l, ig] = elog
					Ecost_new[k, l, ig] = eco
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

		# convergence check
		err_m = maximum(abs.(m_upd .- eqm.m_young))
		err_H = maximum(abs.(H_upd .- eqm.Htilde) ./ max.(abs.(eqm.Htilde), 1e-8))
		err_t = maximum(abs.(t_upd .- eqm.t))
		err = max(err_m, err_H, err_t)

		eqm_new = Eqm(m_upd, M_upd, H_upd, t_upd, Elogh_upd, Ecost_upd)

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
	A_O = mean(a_by_occ)  # composite non-teaching productivity (simple v1)

	τ_w = Array{Float64}(d["τ_w"])
	τ_e = Array{Float64}(d["τ_e"])
	# τ_* dimensions: (nOccNonTeach, nG). We collapse to one composite by simple mean.
	τw_f = mean(τ_w[:, 1])
	τw_m = mean(τ_w[:, 2])
	τe_f = mean(τ_e[:, 1])
	τe_m = mean(τ_e[:, 2])
	τw_O = Dict(:f => τw_f, :m => τw_m)
	τe_O = Dict(:f => τe_f, :m => τe_m)

	p = Params(
		2,
		[:f, :m],
		α, β, η, σ, μ, ϕ, γ,
		κ,
		A_O,
		τw_O,
		τe_O,
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

