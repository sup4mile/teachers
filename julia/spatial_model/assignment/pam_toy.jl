# Toy: teacher -> student assignment with equal class sizes.
# One location, two occupations (T teaching, O other), no taxes, no location choice.
# Technology and choices follow notes/spatial_teachers_model_ccs.md:
#   h = q (z ϵ_i)^α s^φ e^η,  q = h_T^β N^{-σ},  N = 1/n_T (equal classes)
#   O earns A h, T earns κ h^γ; e paid out of income; U = ln(1-s) + μ ln C.
# Assignment: Gaussian copula between the student's rank in z and her teacher's rank in h_T,
#   ζ_teacher = ϱ ζ_z + sqrt(1-ϱ²) ξ,   ϱ = 0: random (CCS),  ϱ = 1: positive assortative matching.
# Every teacher has N students for any ϱ (the teacher's rank is uniform), so masses balance.
# Fixed point: the quantile function of ln h_T among teachers (in normal scores) and ln N.
#
# Run from the repository root (about 70 seconds); the output of the last run is in pam_toy_output.txt:
#   julia --project=julia julia/spatial_model/assignment/pam_toy.jl
# `include` it to use `solve` and `moments` without running `main`.

using Distributions, Printf, Statistics

# Parameters; the defaults are those of `Params()` in ../spatial_continuous.jl.
@kwdef struct P
    α::Float64 = 0.30    # ability elasticity of human capital
    φ::Float64 = 0.40    # time (s) elasticity of human capital
    η::Float64 = 0.20    # goods (e) elasticity of human capital
    σ::Float64 = 0.25    # class-size elasticity: q = h_T^β N^{-σ}
    γ::Float64 = 0.80    # human-capital elasticity of the teaching wage κ h^γ
    β::Float64 = 0.15    # elasticity of human capital in the teacher's human capital
    μ::Float64 = 1.0     # weight on log consumption
    A::Float64 = 1.5     # productivity of the non-teaching occupation (wage A h)
    κ::Float64 = 0.80    # teaching-wage shifter (`main` sets it to hit a teacher share)
    ρz::Float64 = 0.9    # AR(1) persistence of ln z across generations
    σξ::Float64 = 0.20   # AR(1) innovation sd of ln z
    σϵ::Float64 = 0.30   # sd of the occupational shocks ln ϵ_T, ln ϵ_O (iid, unit mean)
    ϱ::Float64 = 0.0     # assignment: copula correlation of student and teacher ranks (0 random, 1 PAM)
end

# Naming: a leading `l` is a log (lq = ln q, lz = ln z, lhT = ln h_T, lN = ln N, leT = ln ϵ_T), and ζ is a normal score
# (a rank u written as Φ⁻¹(u)), so ζz is a student's rank in z and ζg the grid of teacher ranks.

const NORM = Normal()

"Uniform grid of `n` standard-normal points on [lo, hi] and their probability weights (sum to 1)."
function ugrid(lo, hi, n)
    x = collect(range(lo, hi; length = n))
    w = pdf.(NORM, x)
    return x, w ./ sum(w)
end

"Linear interpolation on a uniform grid with linear extrapolation."
function lin(xg, yg, x)
    n = length(xg); dx = xg[2] - xg[1]
    i = clamp(floor(Int, (x - xg[1]) / dx) + 1, 1, n - 1)
    return yg[i] + (yg[i+1] - yg[i]) * (x - xg[i]) / dx
end

# Schooling time share by occupation: a constant, from max ln(1-s) + (return to ln s) ln s (baseline Proposition 2)
sO(p) = (a = p.μ * p.φ / (1 - p.η); a / (1 + a))
sT(p) = (a = p.μ * p.γ * p.φ / (1 - p.γ * p.η); a / (1 + a))

"Teaching branch at school input lq = ln q and ln(zϵ_T) = lx: (ln C, ln h)."
function teach(p, lq, lx)
    b  = lq + p.α * lx + p.φ * log(sT(p))                       # ln h before goods: ln(q x^α s^φ)
    le = (log(p.γ * p.η * p.κ) + p.γ * b) / (1 - p.γ * p.η)     # ln e from the FOC of κ h^γ - e
    return log(1 / (p.γ * p.η) - 1) + le, b + p.η * le          # C = (1/(γη) - 1) e,  ln h = b + η ln e
end

"Other branch at ln ϵ_O = 0: (ln C, ln h); both have slope α/(1-η) in ln ϵ_O."
function other0(p, lq, lz)
    b  = lq + p.α * lz + p.φ * log(sO(p))                       # ln h before goods, at ϵ_O = 1
    le = (log(p.η * p.A) + b) / (1 - p.η)                       # ln e from the FOC of A h - e
    return log(1 / p.η - 1) + le, b + p.η * le                  # C = (1/η - 1) e,  ln h = b + η ln e
end

"Cell (lq, lz, ln ϵ_T): P(teach), ln h_T, the cutoff c of ln ϵ_O (teach iff below), and ln h_O at ln ϵ_O = 0."
function cell(p, lq, lz, leT)
    lCT, lhT = teach(p, lq, lz + leT)
    lCO0, lhO0 = other0(p, lq, lz)
    a1 = p.α / (1 - p.η)                                                      # slope of ln C_O and ln h_O in ln ϵ_O
    c  = ((log(1 - sT(p)) + p.μ * lCT - log(1 - sO(p))) / p.μ - lCO0) / a1    # U_T = U_O at ln ϵ_O = c
    pT = cdf(NORM, (c + p.σϵ^2 / 2) / p.σϵ)                                   # ln ϵ_O ~ N(-σϵ²/2, σϵ²)
    return pT, lhT, c, lhO0
end

"""
    solve(p; nz, nξ, nϵ, damp, tol, maxit, verbose)

Stationary equilibrium of the toy: a damped fixed point in `Qf`, the quantile function of ln h_T among teachers (on the
normal-score grid `ζg`), and `lN`, the log class size. Each iteration assigns teachers to students with the current
`Qf`, solves the students' problems, and rebuilds `Qf` and `lN` from those who choose to teach.
"""
function solve(p; nz = 121, nξ = 41, nϵ = 161, damp = 0.5, tol = 1e-9, maxit = 2000, verbose = false)
    σz = p.σξ / sqrt(1 - p.ρz^2)                              # stationary sd of ln z
    ζz, wz = ugrid(-4.5, 4.5, nz)                             # students' ranks in z (normal scores); ln z = σz ζz
    ξ, wξ  = p.ϱ ≥ 1 ? ([0.0], [1.0]) : ugrid(-4.5, 4.5, nξ)  # assignment noise (none under PAM)
    xϵ, wϵ = ugrid(-5.0, 5.0, nϵ)
    leTg   = -p.σϵ^2 / 2 .+ p.σϵ .* xϵ                        # grid of ln ϵ_T (ϵ_O is integrated in closed form)
    ζg     = collect(range(-3.25, 3.25; length = 131))        # teacher ranks (normal scores) at which Qf is stored
    Qf     = fill(0.0, length(ζg))                            # ln h_T of the teacher at rank ζg; start: identical teachers
    lN     = log(15.0)                                        # ln class size, N = 1/n_T
    n      = nz * length(ξ) * nϵ                              # number of (z, ξ, ϵ_T) cells
    val    = zeros(n); wgt = zeros(n)                         # each cell's ln h_T if she teaches, and its teacher mass
    nx     = 401; hist = zeros(nx); cdfx = zeros(nx)          # histogram and CDF of ln h_T among teachers
    it = 0; err = Inf; nT = NaN
    while it < maxit
        it += 1; m = 0
        # 1. students' problems: a student at rank ζz[i] with noise ξ[k] gets the teacher at rank ϱ ζz + √(1-ϱ²) ξ
        for i in 1:nz, k in eachindex(ξ)
            lq = p.β * lin(ζg, Qf, p.ϱ * ζz[i] + sqrt(max(1 - p.ϱ^2, 0)) * ξ[k]) - p.σ * lN   # ln q = β ln h_T - σ ln N
            for j in 1:nϵ
                pT, lhT, _, _ = cell(p, lq, σz * ζz[i], leTg[j])
                m += 1; val[m] = lhT; wgt[m] = wz[i] * wξ[k] * wϵ[j] * pT
            end
        end
        nT = sum(wgt)                                         # mass of teachers per student
        # 2. the teachers' CDF of ln h_T, on a grid centred and scaled by the cells' mean and sd, with each cell's mass
        #    split linearly between its two neighbouring grid points (continuous in the cells' values, unlike a sorted
        #    empirical quantile, which stalls the iteration)
        mT = sum(wgt .* val) / nT; sT_ = sqrt(sum(wgt .* (val .- mT) .^ 2) / nT)
        x0 = mT - 9sT_; dx = 18sT_ / (nx - 1)
        fill!(hist, 0.0)
        for m2 in 1:n
            t = clamp((val[m2] - x0) / dx, 0.0, nx - 1 - 1e-9); b0 = floor(Int, t); f = t - b0
            hist[b0+1] += wgt[m2] * (1 - f); hist[b0+2] += wgt[m2] * f
        end
        cumsum!(cdfx, hist); cdfx ./= nT
        # 3. the new quantile function: invert the CDF at the ranks Φ(ζg)
        Qn = similar(Qf); b1 = 1
        for g in eachindex(ζg)
            u = cdf(NORM, ζg[g])
            while b1 < nx - 1 && cdfx[b1+1] < u
                b1 += 1
            end
            lo = cdfx[b1]; hi = cdfx[b1+1]
            Qn[g] = x0 + dx * (b1 - 1 + (hi > lo ? clamp((u - lo) / (hi - lo), 0.0, 1.0) : 0.0))
        end
        lNn = -log(nT)
        # 4. damped update
        err = max(maximum(abs.(Qn .- Qf)), abs(lNn - lN))
        verbose && it % 20 == 0 && @printf("  it %4d  err %.2e  nT %.4f\n", it, err, nT)
        Qf .= (1 - damp) .* Qf .+ damp .* Qn; lN = (1 - damp) * lN + damp * lNn
        err < tol && break
    end
    return (; p, Qf, ζg, lN, nT, it, err, σz, ζz, wz, ξ, wξ, leTg, wϵ)
end

"""
    moments(s)

Moments of the solution `s` of `solve`, as a named tuple:

    nT, N      teacher share and class size
    sdT, mhT   sd of ln h_T and mean h_T among teachers
    b, R2      slope in ln z of the assignment map E[ln h_T of one's teacher | z], and the R² of that line
    zT         teachers' mean ln z, in sd of ln z
    sdq        sd of ln q across students;  corr_qz  its correlation with ln z
    sdh, sdy   sd of ln h and of ln income, all adults
    rho_pc     parent-child correlation of ln h
    Y          output of the non-teachers, per adult
"""
function moments(s)
    p = s.p; σz = s.σz
    a1 = p.α / (1 - p.η)          # slope of ln h_O in ln ϵ_O
    m0 = -p.σϵ^2 / 2              # mean of ln ϵ_O
    nz = length(s.ζz)
    g = zeros(nz)                 # E[ln h | z], for the parent-child covariance
    # accumulators: first and second moments of ln h (all adults), ln h_T (teachers), ln q (students) and ln income
    Eh = 0.0; Eh2 = 0.0; Y = 0.0; ElzT = 0.0; ElhT = 0.0; ElhT2 = 0.0
    Elq = 0.0; Elq2 = 0.0; Elqz = 0.0; Ely = 0.0; Ely2 = 0.0; EhT = 0.0
    for i in 1:nz, k in eachindex(s.ξ)
        lz = σz * s.ζz[i]
        lq = p.β * lin(s.ζg, s.Qf, p.ϱ * s.ζz[i] + sqrt(max(1 - p.ϱ^2, 0)) * s.ξ[k]) - p.σ * s.lN
        wq = s.wz[i] * s.wξ[k]
        Elq += wq * lq; Elq2 += wq * lq^2; Elqz += wq * lq * lz
        for j in eachindex(s.leTg)
            pT, lhT, c, lhO0 = cell(p, lq, lz, s.leTg[j])
            w = wq * s.wϵ[j]
            cs  = (c - m0) / p.σϵ                              # standardised cutoff of ln ϵ_O
            λ   = pdf(NORM, cs) / max(ccdf(NORM, cs), 1e-300)  # E[x | x > cs]
            x1  = λ; x2 = 1 + cs * λ                           # first two moments of the truncated standard normal
            mO  = lhO0 + a1 * (m0 + p.σϵ * x1)                 # E[ln h_O | not teaching]
            mO2 = (lhO0 + a1 * m0)^2 + 2 * (lhO0 + a1 * m0) * a1 * p.σϵ * x1 + (a1 * p.σϵ)^2 * x2   # E[(ln h_O)² | not teaching]
            lh1 = pT * lhT + (1 - pT) * mO                     # E[ln h] in the cell
            g[i] += s.wξ[k] * s.wϵ[j] * lh1
            Eh  += w * lh1; Eh2 += w * (pT * lhT^2 + (1 - pT) * mO2)
            # non-teachers' output A E[h_O 1{not teaching}], a lognormal partial expectation
            Y   += w * p.A * exp(lhO0 + a1 * m0 + (a1 * p.σϵ)^2 / 2) * ccdf(NORM, cs - a1 * p.σϵ)
            ElzT += w * pT * lz; ElhT += w * pT * lhT; ElhT2 += w * pT * lhT^2; EhT += w * pT * exp(lhT)
            lyT = log(p.κ) + p.γ * lhT; lA = log(p.A)          # ln income: teachers κ h^γ, others A h
            Ely += w * (pT * lyT + (1 - pT) * (lA + mO)); Ely2 += w * (pT * lyT^2 + (1 - pT) * (lA^2 + 2lA * mO + mO2))
        end
    end
    # parent-child covariance of ln h: shocks and teacher draws are fresh, so it runs through z only
    xg, wx = ugrid(-4.5, 4.5, 61); cov = 0.0
    for i in 1:nz
        Eg = sum(wx[k] * lin(s.ζz, g, p.ρz * s.ζz[i] + sqrt(1 - p.ρz^2) * xg[k]) for k in eachindex(xg))
        cov += s.wz[i] * g[i] * Eg
    end
    cov -= Eh^2
    # the assignment map: E[ln h_T of one's teacher | z], its slope in ln z and how close to log-linear it is
    mp = [sum(s.wξ[k] * lin(s.ζg, s.Qf, p.ϱ * s.ζz[i] + sqrt(max(1 - p.ϱ^2, 0)) * s.ξ[k]) for k in eachindex(s.ξ)) for i in 1:nz]
    lzg = σz .* s.ζz; mb = sum(s.wz .* mp)
    b  = sum(s.wz .* (mp .- mb) .* lzg) / σz^2
    R2 = b^2 * σz^2 / max(sum(s.wz .* (mp .- mb) .^ 2), 1e-300)
    sdT = sqrt(ElhT2 / s.nT - (ElhT / s.nT)^2)
    return (; nT = s.nT, N = exp(s.lN), sdT, b, R2, zT = ElzT / s.nT / σz, sdh = sqrt(Eh2 - Eh^2), rho_pc = cov / (Eh2 - Eh^2),
            Y, sdq = sqrt(max(Elq2 - Elq^2, 0)), corr_qz = (Elqz) / (σz * sqrt(max(Elq2 - Elq^2, 1e-300))),
            sdy = sqrt(Ely2 - Ely^2), it = s.it, err = s.err, mhT = EhT / s.nT)
end

"SD of ln h_T without selection into teaching: root of σT²[(1-γη)² - β²] - 2αβϱσzσT - α²(σz² + σϵ²) = 0."
function sdT_noselect(p)
    σz = p.σξ / sqrt(1 - p.ρz^2); a = (1 - p.γ * p.η)^2 - p.β^2
    a ≤ 0 && return Inf
    bq = 2 * p.α * p.β * p.ϱ * σz; cq = p.α^2 * (σz^2 + p.σϵ^2)
    return (bq + sqrt(bq^2 + 4a * cq)) / (2a)
end

"A copy of the parameters `p` with the fields in `kw` replaced."
with(p; kw...) = P(; (f => getfield(p, f) for f in fieldnames(P))..., kw...)

"κ that gives a teacher share `share` at the parameters `p` (coarse grids)."
function calib_κ(p; share = 0.07)
    lo, hi = 0.02, 20.0                 # bisection in ln κ; the teacher share rises with κ
    for _ in 1:36
        mid = sqrt(lo * hi)
        solve(with(p; κ = mid); nz = 61, nξ = 21, nϵ = 81, tol = 1e-7, maxit = 300).nT < share ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end

"The two experiments whose output is in pam_toy_output.txt."
function main()
    base = P()
    # 1. random assignment to PAM (ϱ from 0 to 1) at two values of β. Columns: see `moments`; "(no-sel)" is
    #    `sdT_noselect`, and Y is also shown relative to ϱ = 0.
    println("κ is set at each β to give a 7% teacher share under random assignment (ϱ = 0), then held fixed across ϱ.\n")
    println("   β     ϱ  |   κ     share    N   sd lnhT (no-sel)  slope b   R²lin  z̄_T(sd)  sd ln q | sd ln h  sd ln y  ρ_pc(ln h)    Y            | iters   err")
    for β in (0.15, 0.30)
        κ = calib_κ(with(base; β)); Y0 = NaN
        for ϱ in (0.0, 0.5, 0.8, 1.0)
            p = with(base; κ, β, ϱ)
            m = moments(solve(p))
            ϱ == 0 && (Y0 = m.Y)
            @printf("  %.2f  %.2f | %.4f  %.4f %5.1f   %.4f (%.4f)   %6.3f   %.4f   %+.3f    %.4f  |  %.4f   %.4f    %.4f    %.4f (%+.2f%%) | %4d  %.1e\n",
                    β, ϱ, κ, m.nT, m.N, m.sdT, sdT_noselect(p), m.b, m.R2, m.zT, m.sdq, m.sdh, m.sdy, m.rho_pc, m.Y, 100 * (m.Y / Y0 - 1), m.it, m.err)
        end
    end
    # 2. existence: under PAM, raise β towards the bound 1 - γη, iterating without damping
    println("\nExistence bound β < 1 - γη = ", 1 - base.γ * base.η, ": PAM (ϱ = 1), κ recalibrated to a 7% share at each β, undamped iteration")
    for β in (0.15, 0.30, 0.45, 0.60, 0.75, 0.80, 0.90)
        κ = calib_κ(with(base; β, ϱ = 1.0))
        p = with(base; κ, β, ϱ = 1.0)
        s = solve(p; damp = 1.0, maxit = 600, tol = 1e-8)
        m = moments(s)
        @printf("  β %.2f  κ %8.4f  iters %4d  err %.1e  share %.4f  sd lnhT %.4f (no-sel %.4f)  slope b %.3f  z̄_T %+.3f\n", β, κ, s.it, s.err, m.nT, m.sdT, sdT_noselect(p), m.b, m.zT)
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
