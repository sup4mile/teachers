# =============================================================================
# Diagnostics at a point ϑ (T6; julia/spatial_model/calibration/spatial_calibration.md §4).
#
#   the moment and diagnostics tables at ϑ, always;
#   --consistency  Table 2 moments in the full equilibrium against the first stage
#                  and the data (§4, *External-block consistency*);
#   --jacobian     the standardized moment Jacobian: singular values, projections;
#   --polish       Levenberg–Marquardt refinement of ϑ, run before the other
#                  diagnostics; each iteration's forward-difference Jacobian and its
#                  trial steps are spread over the workers.
#
# With E = julia/spatial_model/calibration/estimation and D = julia/spatial_model/calibration/spatial_diagnose.jl:
#   julia --project=E D --theta theta0 --consistency
#   julia --project=E D --theta julia/spatial_model/calibration/runs/<id> --polish --jacobian --workers 16
#
# Options (defaults in brackets)
#   --theta SPEC           theta0; a TikTak run directory (its best point); a TOML file or
#                          a directory holding theta.toml (from --out); or v1,v2,… [theta0]
#   --workers N            local worker processes [SLURM_CPUS_PER_TASK − 1, else 0]
#   --polish               Levenberg–Marquardt from ϑ;  --max-iter N  iterations [25]
#   --jacobian             central-difference Jacobian at the (polished) ϑ
#   --consistency          external-block consistency report
#   --out DIR              write report.txt and theta.toml to DIR
#   --grid name=value      override a SOLVER setting (repeatable), e.g. Nz=9, nXO=96: the
#                          §4 grid and quadrature checks; the frozen reference is unchanged
#   --free name            release an extension parameter of ϑ (r_p, ζ/zeta, c_T, c_Tz, ρw/rho_w and the M3 entries
#                          rho_c, c_C, c_Cz, rho_T0, vartheta_s, s_z, sigma_eps, dlogphi, c_T0, the peer
#                          effect psi_a, ψ_a in Q_l = (2H̃_T/M_l)^σ ā_l^ψ_a, and the switches of 2026-10-15
#                          kappa_c (commuting friction; split), theta_nu (σν/μ) and zeta_BA (BA suburban
#                          amenity; college), and b_T (teachers' suburban workplace taste, closure R; split))
#                          to its EXTENSION_BOUNDS (repeatable); they are fixed at 0
#                          otherwise (for s_z, σ_ϵ, dlogφ, θ_ν: at the external block's s_z, σϵ, φ, θν)
#   --district SCHEME      the four district targets of a partition scheme (locale, revenue,
#                          ses, …) from data/spatial/estimates/settings/district_settings_targets.toml
#   --weight diag|corr     criterion weights: diagonal [diag], or the bootstrap correlations of
#                          the district targets of the --district scheme (default locale)
#   --floor-tol X          largest population share allowed at the consumption floor, with or
#                          without a feasible location (`floor_mass + infeasible_mass`); the
#                          evaluation fails above it [1e-8]
#   --fix, --external, --measure, --target, --target-se, --activate, --deactivate,
#   --wtp-consumption, --reference:
#                          as in spatial_estimate.jl, which also describes the residence/workplace split:
#                          --external split=1 --external τk=0.03 --external jobs_city=0.45 --fix ρw=0.1
#                          M3 (joint education and occupation choice, model §3c; it composes with the split):
#                          `--external college=1` turns it on, with `--external sBA=…`, `sT0=…` (the
#                          non-BA teachers' schooling floor), `Fc=…`, `sigma_d=…`, `omega_z=…`;
#                          `--fix rho_c=0.2 --fix c_C=0 --fix c_Cz=-0.5 --fix rho_T0=0.3
#                          --fix vartheta_s=1.1 --fix c_T0=…` set its ϑ entries at a point; its
#                          moments are in the report, and `--activate ba_share_nonteach` etc. move
#                          them into the criterion (`--activate pooled_p9010` computes the pooled 90/10)
#                          `--external move_util=1` (2026-10-15, candidate U): ϑ's r_m is a utility moving
#                          cost (τmove, utils) instead of the goods cost r_m C̄
# =============================================================================

using Distributed
using Printf
using TOML
using TikTak: load_estimates

const HERE = @__DIR__
include(joinpath(HERE, "spatial_cli.jl"))

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "workers" => default_workers(), "theta" => "theta0", "polish" => false, "max-iter" => 25,
        "jacobian" => false, "consistency" => false, "out" => nothing, "wtp-consumption" => nothing,
        "reference" => nothing, "district" => nothing, "weight" => "diag", "floor-tol" => nothing);
    flags   = ("polish", "jacobian", "consistency"),
    ints    = ("workers", "max-iter"),
    floats  = ("wtp-consumption", "floor-tol"),
    strings = ("theta", "out", "reference", "district", "weight"),
    lists   = ("fix", "external", "measure", "target", "target-se", "grid"),
    multi   = ("activate", "deactivate", "free"),
    script  = "spatial_diagnose.jl")
OPTS["workers"] > 0 && addprocs(OPTS["workers"]; exeflags = "--project=$(Base.active_project())")

@everywhere begin
    using LinearAlgebra
    BLAS.set_num_threads(1)
    include(joinpath($HERE, "spatial_calibrate.jl"))
end
include(joinpath(HERE, "spatial_first_stage.jl"))   # ACS_WAGES, NLSY, RoyCells, … for --consistency

# -----------------------------------------------------------------------------
# The point ϑ
# -----------------------------------------------------------------------------
"""
ϑ from a --theta specification (see the header). A TOML file must list `theta_names`
equal to `THETA_NAMES` or to a leading part of it of at least eight names (files saved
before some extension parameters existed, read with those extensions off; older files
still hold log κ̄ as `logκ`, not log κ̃); TikTak journals carry no names, so only their
length is checked. Short vectors are padded with the missing extensions off (`pad_theta`).
"""
function read_theta(spec::AbstractString)
    spec == "theta0" && return copy(THETA0)
    checklen(θ) = (N_THETA_BASE ≤ length(θ) ≤ length(THETA_NAMES) ||
                   error("--theta $spec: $(length(θ)) values, ϑ has $(length(THETA_NAMES))"); pad_theta(θ))
    if isdir(spec)
        isfile(joinpath(spec, "theta.toml")) && return read_theta(joinpath(spec, "theta.toml"))
        return checklen(load_estimates(spec; limit = 1)[1])
    end
    if isfile(spec)
        d = TOML.parsefile(spec)
        names = get(d, "theta_names", nothing)
        (names isa AbstractVector && N_THETA_BASE ≤ length(names) ≤ length(THETA_NAMES) &&
         names == collect(string.(THETA_NAMES[1:length(names)])) && length(d["theta"]) == length(names)) ||
            error("--theta $spec: theta_names $names differ from $(collect(string.(THETA_NAMES)))")
        return checklen(Float64.(d["theta"]))
    end
    return checklen(parse.(Float64, split(spec, ',')))
end

# -----------------------------------------------------------------------------
# Levenberg–Marquardt polish
# -----------------------------------------------------------------------------
"""
    polish(model, targets, θ0, bounds; max_iter, h, tol, io) -> (; θ, f, r, iters)

Levenberg–Marquardt on the standardized residuals r(ϑ) = W(m(ϑ) − target)/se over the
free parameters (those with lo < hi in `bounds`), in box-width units; `W` whitens
correlated targets (`target_whitener`; the identity gives diagonal weights). After each
accepted step it takes a forward-difference Jacobian (step `h` box widths, inward at
an upper bound); every iteration tries four damping factors. These solves run in
parallel with `pmap` (`std_residuals`, which is defined on every worker and logs the message and ϑ of a
point that fails to solve, so that a zeroed Jacobian column or a failed start has its cause in the log).
A step is kept only if it lowers Σr²; steps are clipped to the box. Stops when the relative
improvement falls below `tol` or damping exceeds 1e8.
"""
function polish(model::SpatialModel, targets, θ0, bounds; max_iter = 25, h = 1e-3, tol = 1e-6, io = stdout,
                W = Matrix{Float64}(I, length(active_targets(targets)), length(active_targets(targets))))
    act   = active_targets(targets)
    lo    = first.(bounds); hi = last.(bounds)
    width = hi .- lo
    free  = findall(lo .< hi)
    whiten(v) = all(isfinite, v) ? W * v : v
    θ = clamp.(collect(Float64, θ0), lo, hi)
    r = whiten(std_residuals(model, θ, act))
    f = sum(abs2, r)
    isfinite(f) || error("polish: the start point fails to solve (the failure and ϑ are logged above)")
    μ = 1e-2
    @printf(io, "\nLevenberg–Marquardt polish over %s: start criterion %.6g\n", join(THETA_NAMES[free], ", "), f)
    iters = 0
    J = nothing
    for it in 1:max_iter
        iters = it
        if J === nothing                     # new point: new Jacobian
            steps = [(θ[j] + h * width[j] > hi[j] ? -h : h) * width[j] for j in free]
            pts   = [(x = copy(θ); x[j] += steps[i]; x) for (i, j) in enumerate(free)]
            R     = whiten.(pmap(x -> std_residuals(model, x, act), pts))
            J     = reduce(hcat, [(R[i] .- r) ./ (steps[i] / width[free[i]]) for i in eachindex(free)])
            bad   = [!all(isfinite, J[:, i]) for i in eachindex(free)]
            any(bad) && (@printf(io, "  ⚠ Jacobian columns %s failed; zeroed\n", join(THETA_NAMES[free[bad]], ", ")); J[:, bad] .= 0.0)
        end
        A, g  = J' * J, J' * r
        D     = Diagonal(max.(diag(A), 1e-12))
        μs    = μ .* [0.1, 1.0, 10.0, 100.0]
        cands = map(μs) do m
            Δ = -((A + m * D) \ g)
            x = copy(θ)
            x[free] .+= Δ .* width[free]
            clamp.(x, lo, hi)
        end
        Rc = whiten.(pmap(x -> std_residuals(model, x, act), cands))
        fc = [sum(abs2, rc) for rc in Rc]
        k  = argmin(fc)
        if fc[k] < f
            rel = (f - fc[k]) / f
            θ, r, f = cands[k], Rc[k], fc[k]
            μ = μs[k] / 10
            J = nothing
            @printf(io, "  it %2d: criterion %.6g (−%.2e rel.)  damping %.1e  ϑ = %s\n", it, f, rel, μs[k], fmtvec(θ))
            flush(io)
            (rel < tol || f < 1e-16) && break   # an exact fit keeps improving ~17% a step at f ≈ 1e-20
        else
            μ = 100 * maximum(μs)
            @printf(io, "  it %2d: no improvement (best trial %.6g); damping → %.1e\n", it, fc[k], μ)
            flush(io)
            μ > 1e8 && break
        end
    end
    return (; θ, f, r, iters)
end

# -----------------------------------------------------------------------------
# External-block consistency (§4)
# -----------------------------------------------------------------------------
"""
    consistency_report(io, sol, meas)

Table 2 in the full equilibrium at `sol` against its first-stage counterparts and the
data: occupation shares among non-teachers by gender, cell 90/10s and their
wage-sample-weighted pool, teachers' 90/10s, and the NLSY moments (latent wage slope,
mother–child score correlation).
"""
function consistency_report(io::IO, sol, meas::Measurement; oc = occupation_consistency(sol, meas),
                            occ_sample::Symbol = :acs2009_13)
    @assert oc.occ == [i for i in 1:N_OCC if i != T_OCC] "non-teaching occupations out of the ACS order"
    data  = acs_data(occ_sample)
    wdata = occ_sample === :acs2016_19 ? ACS_2016_19 : ACS_WAGES
    tdata = occ_sample === :acs2016_19 ? ACS_2016_19.teacher_p9010 : TEACHER_P9010
    p   = sol.p
    mkt = 1:N_MKT
    b   = 1 / (1 - p.η)
    s_z = p.σξ / sqrt(1 - p.ρz^2)
    R1  = cell_p9010s(RoyCells(; data), p.α * p.σϵ, zdist(sol.gr.Nz, s_z; ρz = p.ρz); α = p.α, b)   # first stage at this block
    D, w = data_p9010(wdata), wage_weights(wdata)
    Rge = oc.p9010[mkt, :]
    tgt = hcat(occ_shares(1; data), occ_shares(2; data))
    println(io, "\n===== External-block consistency: Table 2 in the full equilibrium (§4), data $occ_sample =====")
    @printf(io, "  pooled non-teacher 90/10: GE %.4f  first stage %.4f  data %.4f\n", pooled(Rge; w), pooled(R1; w), wdata.pooled)
    @printf(io, "  by gender (GE / first stage / data): men %.3f / %.3f / %.3f;  women %.3f / %.3f / %.3f\n",
            pooled(Rge[:, 1]; w = w[:, 1]), pooled(R1[:, 1]; w = w[:, 1]), pooled(D[:, 1]; w = w[:, 1]),
            pooled(Rge[:, 2]; w = w[:, 2]), pooled(R1[:, 2]; w = w[:, 2]), pooled(D[:, 2]; w = w[:, 2]))
    @printf(io, "  shares among non-teachers: max |GE − target| %.4f (men), %.4f (women); home production GE %.4f, %.4f vs %.4f, %.4f\n",
            maximum(abs, oc.shares[:, 1] .- tgt[:, 1]), maximum(abs, oc.shares[:, 2] .- tgt[:, 2]),
            oc.shares[HP_OCC, :]..., tgt[HP_OCC, :]...)
    @printf(io, "  teachers' 90/10: GE %.3f (men %.3f, women %.3f) vs data %.3f (%.3f, %.3f)\n",
            oc.teacher_p9010, oc.teacher_p9010_g..., tdata.pooled, tdata.men, tdata.women)
    @printf(io, "  NLSY moments: latent wage slope %.4f vs %.4f (SE %.4f);  mother–child ρ %.4f (all parents %.4f) vs %.4f (SE %.4f);  c = %.3f\n",
            oc.wage_slope, NLSY.slope, NLSY.slope_se, oc.rho_pc, oc.rho_pc_all, NLSY.rho_pc, NLSY.rho_pc_se, oc.c)
    # mean log wages (T4 check): model units differ from dollars, so compare deviations
    # from the wage-sample-weighted market mean
    ML  = hcat(wdata.mean_log_men, wdata.mean_log_women)
    Mge = oc.mean_logw[mkt, :]
    wbar(X) = sum(w .* X) / sum(w)
    dG, dD = Mge .- wbar(Mge), ML .- wbar(ML)
    sdw(X) = sqrt(wbar((X .- wbar(X)) .^ 2))
    @printf(io, "  mean log wages by cell (deviation from the market mean): weighted correlation %.3f, RMSE %.3f, SD GE %.3f vs data %.3f\n",
            wcorr(vec(dG), vec(dD), vec(w)), sqrt(wbar((dG .- dD) .^ 2)), sdw(Mge), sdw(ML))
    mG(g) = sum(w[:, g] .* Mge[:, g]) / sum(w[:, g]); mD(g) = sum(w[:, g] .* ML[:, g]) / sum(w[:, g])
    @printf(io, "  gender gap in market mean log wages: GE %.3f vs data %.3f;  teacher premium over the market mean: GE men %+.3f, women %+.3f vs data %+.3f, %+.3f\n",
            mG(1) - mG(2), mD(1) - mD(2), oc.teacher_mean_logw[1] - mG(1), oc.teacher_mean_logw[2] - mG(2),
            wdata.k12_mean_log[1] - mD(1), wdata.k12_mean_log[2] - mD(2))
    @printf(io, "  %-27s %7s %7s | %7s %7s | %6s %6s %6s | %6s %6s %6s | %6s %6s | %6s %6s\n", "occupation", "sh m GE", "target",
            "sh f GE", "target", "m GE", "1st", "data", "f GE", "1st", "data", "lw m", "data", "lw f", "data")
    for i in 1:N_OCC-1
        if i <= N_MKT
            @printf(io, "  %-27s %7.4f %7.4f | %7.4f %7.4f | %6.2f %6.2f %6.2f | %6.2f %6.2f %6.2f | %+6.2f %+6.2f | %+6.2f %+6.2f\n", short_name(i),
                    oc.shares[i, 1], tgt[i, 1], oc.shares[i, 2], tgt[i, 2], Rge[i, 1], R1[i, 1], D[i, 1], Rge[i, 2], R1[i, 2], D[i, 2],
                    dG[i, 1], dD[i, 1], dG[i, 2], dD[i, 2])
        else
            @printf(io, "  %-27s %7.4f %7.4f | %7.4f %7.4f |\n", short_name(i), oc.shares[i, 1], tgt[i, 1], oc.shares[i, 2], tgt[i, 2])
        end
    end
    println(io, "=============================================================================")
    return oc
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
"""
Write ϑ, its criterion and the targeted moments `m` to a TOML file (read back by
--theta). `criterion` is the diagonal Σr²; with whitening `W` (--weight corr) the
weighted criterion is `criterion_weighted`.
"""
function save_theta(path, θ, m, targets; note = "", W = nothing)
    act = active_targets(targets)
    d = Dict{String,Any}("theta" => collect(θ), "theta_names" => collect(string.(THETA_NAMES)),
                         "criterion" => criterion(m, targets), "note" => note,
                         "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
                         "moments" => Dict(string(t.key) => Float64(m[t.key]) for t in act))
    W === nothing || (d["criterion_weighted"] = sum(abs2, W * residuals(m, targets)))
    open(io -> TOML.print(io, d; sorted = true), path, "w")
end

"""
Write every finite scalar and short real vector of `model_moments` to a TOML file
(`moments.toml` next to theta.toml), with the active target keys and every target's
data value, scale and status (table `targets`, as this case defined them): the input of
the experiment summaries.
"""
function save_moments(path, m, targets)
    d = Dict{String,Any}("active" => [string(t.key) for t in active_targets(targets)],
                         "targets" => Dict(string(t.key) => Dict("value" => t.value, "se" => t.se, "active" => t.active)
                                           for t in targets if isfinite(t.value) && isfinite(t.se)))
    for (k, v) in pairs(m)
        if v isa Real && isfinite(v)
            d[string(k)] = Float64(v)
        elseif v isa AbstractVector{<:Real} && length(v) ≤ 25 && all(isfinite, v)
            d[string(k)] = Float64.(collect(v))
        end
    end
    open(io -> TOML.print(io, d; sorted = true), path, "w")
end

function run_diagnose(opts)
    ref = frozen_reference(something(opts["reference"], REFERENCE_FILE), SOLVER)
    (; base, meas, targets, active, bounds, W, floor_tol) = problem_setup(opts, ref)
    weighted = opts["weight"] == "corr"
    grid   = NamedTuple(k => (k in (:Nz, :nϵT, :nXO, :maxit, :hh_maxit) ? Int(v) : v) for (k, v) in opts["grid"])
    solver = merge(SOLVER, grid)
    model  = SpatialModel(base, [t.key for t in active]; meas, solver, floor_tol)
    θ = read_theta(opts["theta"])
    for (name, v) in opts["fix"]
        θ[theta_index(name)] = v
    end
    external_start!(θ, bounds, base)       # a freed s_z, σ_ϵ or θ_ν at its off value 0 starts at the external block's value
    out = opts["out"]
    out === nothing || mkpath(out)
    io  = out === nothing ? stdout : open(joinpath(out, "report.txt"), "w")
    tee(f) = (f(stdout); io === stdout || (f(io); flush(io)))
    tee(o -> println(o, "spatial_diagnose.jl  ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
                     "\n  options: ", join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                                           if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "),
                     "\n  measurement: ", meas, "\n  solver: ", solver, "\n  workers: ", nworkers()))

    if opts["polish"]
        t = @elapsed pr = polish(model, targets, θ, bounds; max_iter = opts["max-iter"], W)
        θ = pr.θ
        tee(o -> @printf(o, "\npolish: %d iterations, %.1f minutes, criterion %.6g\n", pr.iters, t / 60, pr.f))
    end
    sol, m = solve_theta(model, θ)
    tee(o -> begin
        println(o, "\nϑ = (", join((@sprintf("%s = %.5f", n, x) for (n, x) in zip(THETA_NAMES, θ)), ", "), ")")
        moment_table(o, m, targets)
        weighted && @printf(o, "  weighted criterion (district-target correlations, --district %s) = %.4f\n",
                            something(opts["district"], "locale"), sum(abs2, W * residuals(m, targets)))
        diagnostics_table(o, m)
        any(!iszero, θ[theta_index(:ρ_c):theta_index(:c_T0)]) && @printf(o, "  M3 block in effect: college %s, s_z %.4f, σϵ %.4f, φ %.4f; ρc %.3f, cC %.3f, cCz %.3f per unit of log z, ρT0 %.3f, cT0 %.3f, ϑs %.3f, sBA %.3f, sT0 %.3f, Fc %.4f, σd %.3f\n",
                                           sol.p.college ? "on" : "off", sol.p.σξ / sqrt(1 - sol.p.ρz^2), sol.p.σϵ, sol.p.φ,
                                           sol.p.ρc, sol.p.cC, sol.p.cCz, sol.p.ρT0, sol.p.cT0, sol.p.ϑs, sol.p.sBA, sol.p.sT0, sol.p.Fc, sol.p.σd)
        @printf(o, "  (e, s) solve: residual %.1e, at most %d passes in the final sweep; relaxed step at up to %d nodes in a sweep\n",
                sol.s_resid, sol.s_passes, sol.s_relaxed)
    end)
    if opts["consistency"]
        oc = occupation_consistency(sol, meas)
        occ_sample = Symbol(get(Dict(opts["external"]), :occ_sample, :acs2009_13))
        tee(o -> consistency_report(o, sol, meas; oc, occ_sample))
    end
    if opts["jacobian"]
        free = findall(first.(bounds) .< last.(bounds))
        try
            t = @elapsed jac = moment_jacobian(model, θ, targets; free, bounds, map = pmap)
            tee(o -> (print_jacobian(o, jac); @printf(o, "  (%.1f minutes)\n", t / 60)))
        catch e
            tee(o -> println(o, "\n  ⚠ Jacobian failed: ", sprint(showerror, e isa RemoteException ? e.captured.ex : e)))
        end
    end
    if out !== nothing
        save_theta(joinpath(out, "theta.toml"), θ, m, targets; note = "spatial_diagnose.jl --theta $(opts["theta"])",
                   W = weighted ? W : nothing)
        save_moments(joinpath(out, "moments.toml"), m, targets)
        close(io)
        println("\nReport and ϑ written to ", out)
    end
end

run_diagnose(OPTS)
