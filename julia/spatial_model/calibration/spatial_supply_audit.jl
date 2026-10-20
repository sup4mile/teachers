# =============================================================================
# Audit of the teacher-supply and commuting estimands of the residence/workplace split
# (julia/spatial_model/calibration/calibration_log.md §5, "Recommended sequence (2026-10-15)", item 1).
#
# `teacher_supply_elasticity` (spatial_moments.jl) is ε₂₂ − ε₁₂ for ε_{iw′} = ∂log n_i/∂log κ_w′, n_i the mass of
# teachers working at i, after a household re-solve at κ₂ e^{0.01} with (H̃_T, M, t), the grids and Φ held. At C4′ it
# is ≈ 31 against the literature's own-employer elasticities, Falch (2010) 1.4 and Ransom–Sims (2010) 3.7. The ratio
# mixes own and cross responses and several margins; this script separates them. At the solved equilibrium it raises
# the pay of each workplace w′ in turn (κ_w′ → κ_w′ e^{step}, the other κ unchanged) and reports the teacher masses by
# workplace n_w and their total N_T under nested margins, each holding at the base solution whatever it does not list:
#   (a)  workplace only: each teacher node keeps its (e, s), Λ, occupation choice maps, education weights P(d | g, l, z),
#        residence probabilities π_r and Φ; only the conditional workplace probabilities π_{w|r} are recomputed, at the
#        new pay (`workplace_only_probs`). Not defined for the bundled model (workplace = residence): NaN;
#   (b)  + residence: the joint probabilities π_{rw} recomputed at the base (e, s) and Λ (bundled: the location choice);
#   (c₀) + investment, Λ held: (e, s) re-solved at the new pay with Λ held at its base value (a single household sweep,
#        `solve_household(…; maxit = 1)`), integrated with the base choice maps and education weights;
#   (c)  + the altruism term Λ: the household re-solved to its fixed point (policies and Λ together, so that the nodes'
#        location probabilities are all new), integrated with the base maps and weights. This is the investment margin
#        of the audit's specification, which expects that a full re-solve cannot hold Λ; (c₀) is the narrower step
#        before it (they differ by 0.002 or less in the ratio at C4′ and C0);
#   (d)  + occupation and education: the household's own choice maps and weights, Φ and the aggregates fixed: exactly
#        what `teacher_supply_elasticity` does (its ratio at w′ = 2 is reproduced and checked against the moment);
#   (e)  + the stationary distribution Φ, re-solved from the new location probabilities (`stationary_phi`), with the
#        aggregates (Q, t, H̃_T, the positions) still fixed;
#   (f)  full general equilibrium: `solve_ge` at the new pay (positions, Q, t, Φ all move).
# For each margin and shocked workplace w′ it reports ε_{1,w′} = ∂log n_1/∂log κ_w′, ε_{2,w′}, their difference
# ε_{2,w′} − ε_{1,w′} (for w′ = 2 the moment's ratio) and ε_total = ∂log N_T/∂log κ_w′, in two tables:
#   - at `--step` (0.01, the step of the moment), margins (a)–(e) by a one-sided (forward) difference at +step, which is
#     what the moment is and what (d) reproduces; --two-sided adds the backward difference (−step) and the central one;
#   - at `--ge-step` (0.001), margins (a)–(f) by central differences, (log n₊ − log n₋)/(2 step), exact to order step²
#     where a one-sided difference is exact only to order step: the local response, the derivative that the moment
#     approximates. The responses are far from linear over 0.01. At C4′ the ratio of (d) is 30.90 at +0.01 and 29.92 as a
#     derivative, and the general equilibrium much more so: in one-off solves at +0.01, ε_total of (f) is 1.66 against
#     0.94 as a derivative and ε_1 is −5.90 against −6.40; a forward difference at 0.001 is still 7% off in ε_total
#     (1.01), which is why this table is central; (f) shows its forward and backward rows as well.
# (f) is solved at the GE tolerance `--ge-tol` (1e-8), tighter than the calibration's 1e-6 as a margin (the solver's error
# divided by the step is the error of a difference quotient; at C4′ the central differences at 1e-6 and 1e-8 agree to
# 3e-4), from a zero-shock re-solve at that tolerance, against which its one-sided differences are measured. Margins
# (a)–(c) fix the occupation maps, so N_T does not move (ε_total = 0) and the response is a reallocation between
# workplaces. Beyond the margins:
#   (g)  the small-district own-pay elasticity, the model's counterpart of the single-employer estimates: the supply
#        response of a district of negligible weight inside workplace w that raises its pay (`pay_response`), with the
#        analytic derivative of margin (a), which the finite-difference margin approaches as the step falls;
#   (h)  the directional crossing shares P(work in 2 | live in 1), P(work in 1 | live in 2) of teachers (exactly) and of
#        non-teaching and all market workers (pseudo-sample) against the ACS 2018 city/suburb flows, and the net gaps;
#   (i)  the commuting block of the parameters (σν, σw, ρw, τk, κ) and the base teacher masses.
#
# With E = julia/spatial_model/calibration/estimation and A = julia/spatial_model/calibration/spatial_supply_audit.jl:
#   julia --project=E A --theta julia/spatial_model/calibration/runs/exp-c3/c4p --external split=1 --external tau_k=0.03 \
#         --external jobs_city=0.45 --external college=1 --fix r_p=0.04 --floor-tol 1e-5 \
#         --out julia/spatial_model/calibration/runs/experiments-2026-10d/audit/c4p
# (the options of the run that produced ϑ; --free, --activate and --deactivate only define the estimation problem and
# are accepted so that a run's whole option line can be passed). The bundled model needs no options:
# --theta julia/spatial_model/calibration/runs/exact-polish-base.
# Options (defaults in brackets)
#   --theta SPEC       a run directory, a theta.toml file or v1,v2,… (eight values are padded with the extensions off)
#   --out DIR          write audit.txt and audit.toml to DIR [print to the screen only]
#   --step X           the log-pay change of the first table, margins (a)–(e), one-sided [`SUPPLY_STEP` = 0.01, the step of
#                      `teacher_supply_elasticity`, which margin (d) reproduces]
#   --ge-step X        the log-pay change of the second table, margins (a)–(f), central [`GE_STEP` = 0.001]
#   --ge-tol X         the GE tolerance of margin (f) and its zero-shock reference [`GE_TOL` = 1e-8]
#   --two-sided        also −step in the first table: the backward difference, the central difference, and a check of the
#                      analytic margin (a) against the central difference (use a small --step, 1e-4)
#   --no-ge            skip margin (f), the full-general-equilibrium response (four GE solves and the reference)
#   --external, --measure, --fix, --floor-tol, --district, --activate, --deactivate, --free, --target, --target-se,
#   --wtp-consumption, --weight, --reference: the problem options of spatial_diagnose.jl; --fix overrides the read ϑ.
# The pooled 90/10 (`--activate pooled_p9010`) is not computed here. An audit of the college model C4′ takes about six
# minutes on one core, including Julia's start-up: the equilibrium and moments 105 s, the household margins 3 s per
# shock (six shocks), each of the four GE solves 37 s and the zero-shock re-solve 8 s; the bundled model takes about two
# minutes. Run it on a compute node (srun or sbatch), not on the login node.
# =============================================================================

using Printf
using TOML

const HERE = @__DIR__
include(joinpath(HERE, "spatial_cli.jl"))
include(joinpath(HERE, "spatial_calibrate.jl"))

"""
The default `--ge-step`, the log-pay change of the second table (the local response, by central differences): the
responses are far from linear over `SUPPLY_STEP`.
"""
const GE_STEP = 1e-3

"""
The default `--ge-tol`, the GE tolerance of margin (f): a margin below the calibration's (`SOLVER.tol` = 1e-6), whose
error is divided by 2 `GE_STEP` in a central difference.
"""
const GE_TOL = 1e-8

const OPTS = parse_cli(ARGS, Dict{String,Any}(
        "theta" => "", "out" => nothing, "step" => SUPPLY_STEP, "ge-step" => GE_STEP, "ge-tol" => GE_TOL,
        "no-ge" => false, "two-sided" => false,
        "reference" => nothing, "wtp-consumption" => nothing, "district" => nothing, "weight" => "diag",
        "floor-tol" => nothing);
    flags   = ("no-ge", "two-sided"),
    floats  = ("step", "ge-step", "ge-tol", "wtp-consumption", "floor-tol"),
    strings = ("theta", "out", "reference", "district", "weight"),
    lists   = ("fix", "external", "measure", "target", "target-se"),
    multi   = ("activate", "deactivate", "free"),
    script  = "spatial_supply_audit.jl")

"ϑ from a run directory, a theta.toml file or a comma list (eight values are padded with the extensions off)."
function read_run_theta(spec)
    isempty(spec) && error("--theta is required: a run directory, a theta.toml file or v1,v2,…")
    isdir(spec) && return read_run_theta(joinpath(spec, "theta.toml"))
    isfile(spec) && return pad_theta(Float64.(TOML.parsefile(spec)["theta"]))
    return pad_theta(parse.(Float64, split(spec, ',')))
end

# -----------------------------------------------------------------------------
# Data comparators
# -----------------------------------------------------------------------------
"""
ACS 2018 1-year PUMS flows between the city and suburb POWPUMA classes (data/spatial/estimates/settings/commuting_all.md,
the 3 × 3 residence × workplace tables), in percent of the residents of the origin class, by group: (live in city → works
in city, → works in suburb; live in suburb → works in city, → works in suburb). The remainder works in a third class
('other'), so the model's closed city/suburb economy is compared on the block basis, P(work in suburb | live in city) =
the second entry over the first two. Lower bounds on the true crossing: a POWPUMA is a county or group of counties.
"""
const ACS_FLOWS = (; all = (90.04, 9.09, 11.38, 87.24), teach = (90.37, 8.52, 8.19, 90.31),
                     nonteach = (90.03, 9.11, 11.48, 87.14))

"The ACS 2018 net commuting gaps, log(suburb/city) of residents less that of jobs, on the three classes (commuting_all.md §1)."
const ACS_NET_GAP = (; all = 0.177, teach = 0.120, nonteach = 0.179)

"The literature's own-employer elasticities of teacher supply: Falch (2010), Ransom–Sims (2010)."
const LITERATURE_SUPPLY = (; falch = 1.4, ransom_sims = 3.7)

"The calibration's data value of the target `key` (`default_targets`)."
target_value(key::Symbol) = only(t.value for t in default_targets() if t.key === key)

# -----------------------------------------------------------------------------
# Teacher masses and the margins
# -----------------------------------------------------------------------------
"""
    teacher_counts(hh, Φ, Pi, gr, p) -> Vector

Mass of teachers by workplace for the household solution `hh` (its choice maps and education weights), the stationary
distribution `Φ` and the teachers' location probabilities `Pi`: the joint residence × workplace arrays under the split
(`teacher_work_counts`), else the work-location probabilities `Pi_T`, where the workplace is the residence. Under
`p.college` the probabilities are the pair of the education branches.
"""
function teacher_counts(hh, Φ, Pi, gr, p::Params)
    p.split && return teacher_work_counts(hh, Φ, Pi, gr, p)
    (; Nz, L) = gr
    n = zeros(L)
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, Pi, Pi, p, g, l, zi)
            mass = mass0 * br.w
            cm = choice_maps(br.h, g, l, zi, gr)
            for lp in 1:L
                n[lp] += mass * integrate_teach(br.πT[g, l, zi, :, lp], cm, gr)
            end
        end
    end
    return n
end

"""
    teacher_flows(hh, Φ, Pi_TW, gr, p) -> Matrix

Mass of teachers by residence r (rows) and workplace w (columns) under the split: `teacher_work_counts` before the sum
over residences, so that the column sums are the masses by workplace and the row sums those by residence.
"""
function teacher_flows(hh, Φ, Pi_TW, gr, p::Params)
    (; Nz, L) = gr
    F = zeros(L, L)
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, Pi_TW, Pi_TW, p, g, l, zi)
            mass = mass0 * br.w
            cm = choice_maps(br.h, g, l, zi, gr)
            for r in 1:L, w in 1:L
                F[r, w] += mass * integrate_teach(br.πT[g, l, zi, :, r, w], cm, gr)
            end
        end
    end
    return F
end

"`p` with the pay of workplace `w` raised by the factor e^s (the other κ unchanged), as `teacher_supply_elasticity` shocks κ₂."
shock_pay(p::Params, w, s) = with_params(p; κ = [i == w ? p.κ[i] * exp(s) : p.κ[i] for i in eachindex(p.κ)])

"""
    workplace_only_probs(hh, p2, gr, Pi_Tr) -> πTW

Margin (a). The teachers' joint residence × workplace probabilities when only the conditional workplace choice responds
to the pay of `p2`: the residence probabilities `Pi_Tr` of the base solution (`sol.Pi_Tr`) times the π_{w|r} of the
nested logit recomputed at `p2` with each node's base (e, s) and Λ (so consumption is recomputed at the new pay and
nothing else moves). The layout of `sol.Pi_TW`: an array [g, l, zi, k, r, w], a pair of them (one per education
branch) under `college`. Under the split only.
"""
function workplace_only_probs(hh, p2::Params, gr, Pi_Tr)
    (; Nz, L, nϵT) = gr
    out = map(branch_ids(p2)) do d
        hb = branch_of(hh, p2, d)
        πr = branch_probs(Pi_Tr, p2, d)
        π  = zeros(2, L, Nz, nϵT, L, L)
        for g in 1:2, l in 1:L, zi in 1:Nz, k in 1:nϵT
            nd = teach_node(gr.ϵTgrid[k], zi, l, g, hb.Λ, p2, gr; d)
            C  = consumption(nd, hb.sT[g, l, zi, k], hb.eT[g, l, zi, k], p2)
            πw = nested_logsum_probs(location_values(nd, C, p2), l, nd, p2)[4]       # π_{w|r}, [w, r]
            π[g, l, zi, k, :, :] .= πr[g, l, zi, k, :] .* πw'                        # π_r π_{w|r}, [r, w]
        end
        π
    end
    return p2.college ? out : out[1]
end

"""
`hh` with the altruism term replaced by `Λ`, in the household arrays and, under `p.college`, in both education
branches. `solve_household(…; maxit = 1)` solves the nodes at the Λ it is given and then updates Λ; margin (c₀) takes
its policies and puts the base Λ back.
"""
with_lambda(hh, p::Params, Λ) = p.college ?
    merge(hh, (; br1 = merge(hh.br1, (; Λ)), br0 = merge(hh.br0, (; Λ)), Λ)) : merge(hh, (; Λ))

"""
The margins in the order of the header: (key, code, description). The key names the margin in `audit.toml`
(`margin.<key>`): `c` is the full household re-solve of the audit's specification, `c_lambda_held` the narrower step
before it. Without the split there is no workplace apart from the location; without college no education choice.
"""
margin_labels(p::Params) = [
    (:a,             "(a) ", "workplace choice only"),
    (:b,             "(b) ", p.split ? "+ residence choice" : "+ location choice"),
    (:c_lambda_held, "(c₀)", "+ investment (e, s), Λ held"),
    (:c,             "(c) ", "+ altruism term Λ (household re-solved)"),
    (:d,             "(d) ", p.college ? "+ occupation, education" : "+ occupation"),
    (:e,             "(e) ", "+ stationary distribution Φ"),
    (:f,             "(f) ", "full general equilibrium")]

"The keys of the household margins (a)–(e), which `margin_masses` computes."
const PE_MARGINS = (:a, :b, :c_lambda_held, :c, :d, :e)

"""
    margin_masses(sol, w, s; hh_tol, hh_maxit) -> Dict{Symbol,Vector{Float64}}

Teacher masses by workplace under margins (a)–(e) (see the header) after the pay of workplace `w` changes by the log
factor `s`. Everything a margin does not list is the base solution's, `sol`. The household re-solves take the settings
of `teacher_supply_elasticity` and are warm-started at the base policies.
"""
function margin_masses(sol, w, s; hh_tol, hh_maxit)
    (; hh, Φ, gr, p) = sol
    p2  = shock_pay(p, w, s)
    out = Dict{Symbol,Vector{Float64}}()
    out[:a] = p.split ? teacher_counts(hh, Φ, workplace_only_probs(hh, p2, gr, sol.Pi_Tr), gr, p) : fill(NaN, gr.L)
    out[:b] = teacher_counts(hh, Φ, location_probs(hh, p2, gr)[1], gr, p)
    hh1 = solve_household(p2, gr; tol = hh_tol, maxit = 1, hh0 = hh)           # one sweep: (e, s) at the base Λ
    out[:c_lambda_held] = teacher_counts(hh, Φ, location_probs(with_lambda(hh1, p, hh.Λ), p2, gr)[1], gr, p)
    hh2 = solve_household(p2, gr; tol = hh_tol, maxit = hh_maxit, hh0 = hh)
    πT2, πO2 = location_probs(hh2, p2, gr)
    out[:c] = teacher_counts(hh, Φ, πT2, gr, p)
    out[:d] = teacher_counts(hh2, Φ, πT2, gr, p2)
    Φ2, _ = stationary_phi(parent_probs(hh2, p2, gr, πT2, πO2)..., hh2, p2, gr)   # the parents' kernel (model §3d)
    out[:e] = teacher_counts(hh2, Φ2, πT2, gr, p2)
    return out
end

"""
    ge_masses(sol0, w, s; solver) -> (; n, converged, dlogM, dlogHT, dt, dlogQ, sol)

Margin (f): the general equilibrium at the pay of workplace `w` raised by the log factor `s`, solved from the solution
`sol0` (`init`) with `solver`'s settings; `sol` is the new solution and `n` its teacher masses by workplace, with the
movement from `sol0` of the aggregates (log M and log H̃_T by location, the tax rates and log Q). With `s = 0` it
re-solves the equilibrium of `sol0` itself: at a tolerance tighter than that of `sol0`, this is the reference that the
shocked solves are measured against and start from.
"""
function ge_masses(sol0, w, s; solver)
    sol2 = solve_ge(shock_pay(sol0.p, w, s); solver..., verbose = false, init = sol0)
    n = teacher_counts(sol2.hh, sol2.Φ, sol_probs(sol2)[1], sol2.gr, sol2.p)
    return (; n, converged = sol2.converged, dlogM = log.(sol2.M ./ sol0.M), dlogHT = log.(sol2.HT ./ sol0.HT),
              dt = sol2.t .- sol0.t, dlogQ = log.(sol2.gr.Q ./ sol0.gr.Q), sol = sol2)
end

"`ge_masses`, or NaNs with `converged = false` and a printed warning if the solve throws (the audit goes on without that margin)."
function ge_masses_safe(sol0, w, s; solver)
    try
        return ge_masses(sol0, w, s; solver)
    catch err
        println("    ⚠ the general equilibrium failed: ", sprint(showerror, err))
        nan = fill(NaN, sol0.gr.L)
        return (; n = nan, converged = false, dlogM = nan, dlogHT = nan, dt = nan, dlogQ = nan, sol = nothing)
    end
end

"Whether a `ge_masses` run converged to finite masses."
ge_ok(g) = g.converged && all(isfinite, g.n)

"""
    responses(n0, n, s) -> (; eps, ratio, total)

Responses of the teacher masses `n` (by workplace) to a log-pay change `s` against the base masses `n0`:
ε_i = (log n_i − log n0_i)/s, `ratio` = ε_2 − ε_1 and `total` = (log ΣN − log ΣN0)/s.
"""
function responses(n0, n, s)
    eps = log.(n ./ n0) ./ s
    return (; eps, ratio = eps[2] - eps[1], total = log(sum(n) / sum(n0)) / s)
end

# -----------------------------------------------------------------------------
# The small-district elasticity and the analytic derivative of margin (a)
# -----------------------------------------------------------------------------
"""
    node_income(nd, s, e, p) -> Vector

Income y_d = c_d (s^φ e^η)^{γ_b} (times s^{ϑs} for a BA teacher's salary schedule, §3c) of a node at the policies (s, e)
over its destinations d, before the goods costs: the first line of `consumption`, which is y − (1 + τe) e − m.
"""
function node_income(nd::AbstractNode, s, e, p::Params)
    y = nd.c .* (s^p.φ * e^p.η)^nd.γb
    iszero(nd.ϑs) || (y = y .* s^nd.ϑs)
    return y
end

"""
    pay_response(sol) -> (; small, n, workplace, total, y_resid, floor_mass, floor_q_share)

How the teachers' workplace choice responds to pay at the base solution, from the nested logit's derivatives. At a
node the value of living in r and working at w is W_{rw} = μ log C_{rw} + b_w + Λ_r with ∂W_{rw}/∂log κ_w = μ y_{rw}/C_{rw}
(the income y_{rw} = c_{rw}(s^φ e^η)^γ is proportional to the pay κ_w), so with q_{rw} = μ y_{rw}/(σ_in C_{rw}), where σ_in
is the scale of the inner (workplace) nest, σw under the split and σν for the bundled model, whose location choice is
the only logit,
    ∂π_{w″|r}/∂log κ_w′ = π_{w″|r} (1{w″ = w′} − π_{w′|r}) q_{rw′}.
A destination at the consumption floor (C ≤ C_FLOOR, where `smooth_log` is linear with the slope 1/C_FLOOR, a million
times that of an interior destination at C = 1) drops out of every sum, as a pair at the floor with a negligible
probability does in `nested_value_derivs`: its q is an artefact of the linearisation. `floor_mass` is the teacher mass at
such destinations as a share of all teachers (integrated with the linear rule, `k = 1`, as the moments' `floor_mass`) and
`floor_q_share` the largest share (over w) of Σ π q that they would carry with the slope of the linearisation: both are
negligible where the floor checks of `solve_theta` bound the population share at the floor (`--floor-tol`) well below that.
- `small[w]`: the supply elasticity of a district of negligible weight inside workplace w that raises its own pay by
  dlog κ. Its position weight is negligible, so the inner logit's denominator, the household's inclusive value, (e, s),
  the residence and the occupation choice do not respond to it (envelope): only its own share π_{rj} = π_r π_{j|r}
  among the workplace's districts rises, by the factor exp(ΔW_{rw}/σ_in), whence
      ε_small(w) = Σ_nodes mass ∫ teach Σ_r π_{rw} q_{rw} / n_w,
  the mass-weighted mean over teacher nodes and residences of π_{rw} μ y_{rw}/(σ_in C_{rw}) divided by the teachers
  working in w. This is the counterpart of the single-employer estimates (Falch 2010, Ransom–Sims 2010). Under the
  bundled model it treats the district as a location with the logit scale σν, an approximation;
- `workplace[w″, w′]`: the analytic ε of margin (a), ∂log n_w″/∂log κ_w′ = Σ_nodes mass ∫ teach Σ_r π_r ∂π_{w″|r}/∂log κ_w′
  / n_w″, with the residence probabilities, (e, s), Λ, the choice maps and Φ held: the limit of margin (a) as step → 0
  (NaN for the bundled model); `total[w′]` is the ε_total of that margin, ∂log N_T/∂log κ_w′ = Σ_w″ n_w″ ε_{w″,w′}/N_T
  (the probabilities π_{w″|r} sum to one, so it is zero up to rounding; NaN for the bundled model);
- `n[w]`: the teachers working in w; `y_resid`: the largest |C − (y − (1 + τe) e − m)| over the nodes, which checks `node_income`.
"""
function pay_response(sol)
    (; hh, Φ, gr, p) = sol
    (; Nz, L, nϵT) = gr
    Pi  = sol_probs(sol)[1]
    σin = p.split ? p.σw : p.σν
    n   = zeros(L); num = zeros(L); numf = zeros(L); dn = zeros(L, L)
    floor_mass = 0.0
    y_resid = 0.0
    for l in 1:L, zi in 1:Nz
        mass0 = Φ[zi, l] / 2
        iszero(mass0) && continue
        for g in 1:2, br in cell_branches(hh, Pi, Pi, p, g, l, zi)
            mass = mass0 * br.w
            cm   = choice_maps(br.h, g, l, zi, gr)
            π  = zeros(nϵT, L, L); q = zeros(nϵT, L, L)               # [k, r, w]; the bundled model on the diagonal
            qf = zeros(nϵT, L, L); fl = zeros(nϵT, L, L)              # q dropped at the floor (at its slope); the floor indicator
            for k in 1:nϵT
                nd = teach_node(gr.ϵTgrid[k], zi, l, g, br.h.Λ, p, gr; d = br.d)
                s, e = br.h.sT[g, l, zi, k], br.h.eT[g, l, zi, k]
                C, y = consumption(nd, s, e, p), node_income(nd, s, e, p)
                y_resid = max(y_resid, maximum(abs, C .- (y .- (1 + nd.τe) * e .- nd.m)))
                floored = C .≤ C_FLOOR
                qk  = ifelse.(floored, 0.0, p.μ .* y ./ (σin .* max.(C, C_FLOOR)))
                qfk = ifelse.(floored, p.μ .* y ./ (σin * C_FLOOR), 0.0)
                if p.split
                    π[k, :, :]  .= br.πT[g, l, zi, k, :, :]
                    q[k, :, :]  .= reshape(qk, L, L)'                 # [w, r] → [r, w]
                    qf[k, :, :] .= reshape(qfk, L, L)'
                    fl[k, :, :] .= reshape(floored, L, L)'
                else
                    for lp in 1:L
                        π[k, lp, lp] = br.πT[g, l, zi, k, lp]
                        q[k, lp, lp] = qk[lp]; qf[k, lp, lp] = qfk[lp]; fl[k, lp, lp] = floored[lp]
                    end
                end
            end
            πr = dropdims(sum(π; dims = 3); dims = 3)                  # π_r, [k, r]
            for w in 1:L
                n[w]    += mass * integrate_teach(vec(sum(view(π, :, :, w); dims = 2)), cm, gr)
                num[w]  += mass * integrate_teach(vec(sum(view(π, :, :, w) .* view(q, :, :, w); dims = 2)), cm, gr)
                numf[w] += mass * integrate_teach(vec(sum(view(π, :, :, w) .* view(qf, :, :, w); dims = 2)), cm, gr)
                floor_mass += mass * integrate_teach(vec(sum(view(π, :, :, w) .* view(fl, :, :, w); dims = 2)), cm, gr; k = 1)
            end
            p.split || continue
            for wo in 1:L, ws in 1:L
                f = zeros(nϵT)
                for k in 1:nϵT, r in 1:L
                    πr[k, r] > 0 || continue
                    f[k] += π[k, r, wo] * ((wo == ws) - π[k, r, ws] / πr[k, r]) * q[k, r, ws]
                end
                dn[wo, ws] += mass * integrate_teach(f, cm, gr)
            end
        end
    end
    return (; small = num ./ n, n, workplace = p.split ? dn ./ n : fill(NaN, L, L),
              total = p.split ? vec(sum(dn; dims = 1)) ./ sum(n) : fill(NaN, L), y_resid,
              floor_mass = floor_mass / sum(n), floor_q_share = maximum(ifelse.(num .+ numf .> 0, numf ./ (num .+ numf), 0.0)))
end

# -----------------------------------------------------------------------------
# Directional crossing shares
# -----------------------------------------------------------------------------
"""
    crossing_shares(sol, meas) -> NamedTuple or nothing

Residence against workplace by group under the split, each with `p12` = P(work in 2 | live in 1), `p21` = P(work in 1 |
live in 2), `share` = the share who live and work apart and `net_gap` = the log suburb/city ratio of residents less that
of jobs (positive when the group lives in the suburb more than it works there). `teach` is exact (`teacher_flows`);
`all` market workers and `nonteach` (the non-teaching market workers) come from the pseudo-sample of `micro_sample`, home
producers out as in `commute_shares`, and `teach_ms` is the teachers' own pseudo-sample row (the moment
`commute_share_teach`). `nothing` for the bundled model.
"""
function crossing_shares(sol, meas)
    sol.p.split || return nothing
    (; hh, Φ, gr, p) = sol
    rates(F) = (; p12 = F[1, 2] / sum(F[1, :]), p21 = F[2, 1] / sum(F[2, :]), share = (F[1, 2] + F[2, 1]) / sum(F),
                  net_gap = log(sum(F[2, :]) / sum(F[1, :])) - log(sum(F[:, 2]) / sum(F[:, 1])))
    grp = ms_groups(gr, meas.hp)
    hpg = hp_group(grp, meas.hp)
    ms  = micro_sample(sol; n = meas.nquad, groups = grp)
    mkt = market_workers(ms, hpg)
    pseudo(sel) = [sum(ms.w[sel .& (ms.lp .== r) .& (ms.lw .== w)]) for r in 1:gr.L, w in 1:gr.L]
    return (; teach = rates(teacher_flows(hh, Φ, sol.Pi_TW, gr, p)), nonteach = rates(pseudo(mkt .& (ms.grp .> 0))),
              all = rates(pseudo(mkt)), teach_ms = rates(pseudo(ms.grp .== 0)))
end

# -----------------------------------------------------------------------------
# The audit
# -----------------------------------------------------------------------------
"Run `f`, print `label` and the seconds it took, return its value."
function timed(f, label)
    t0 = time()
    val = f()
    @printf("  %-64s %7.1fs\n", label, time() - t0)
    flush(stdout)
    return val
end

"""
    run_audit(sol, m, meas; step, ge_step, ge_tol, ge, two_sided, solver) -> NamedTuple

Everything the audit reports at the solved equilibrium `sol` with moments `m` (`solve_theta`): the base teacher masses
`n0`; `masses[(margin, w′, ±1, tag)]` for margins :a … :e, tag :main at `step` (−1 only with `two_sided`) and tag :fine
at `ge_step` (both signs); `ge[(w′, ±1)]` for margin (f) at `ge_step` (empty if `ge` is false), solved at the tolerance
`ge_tol` from `ge_ref`, the zero-shock re-solve at that tolerance, which is the reference of (f)'s one-sided differences
(`ge_n0` its masses; the base masses `n0` if it did not converge); the small-district elasticities and the analytic
margin (a) (`pay`), the crossing shares and the consistency checks.
"""
function run_audit(sol, m, meas; step, ge_step, ge_tol, ge, two_sided, solver)
    (; hh, Φ, gr, p) = sol
    L = gr.L
    hh_tol, hh_maxit = get(solver, :hh_tol, 1e-7), get(solver, :hh_maxit, 1000)
    n0 = teacher_counts(hh, Φ, sol_probs(sol)[1], gr, p)
    sname(sgn) = sgn > 0 ? "+" : "−"
    masses = Dict{Tuple{Symbol,Int,Int,Symbol},Vector{Float64}}()
    for (tag, s0, signs) in ((:main, step, two_sided ? (1, -1) : (1,)), (:fine, ge_step, (1, -1))), w in 1:L, sgn in signs
        res = timed(@sprintf("margins (a)–(e), κ_%d %s%g", w, sname(sgn), s0)) do
            margin_masses(sol, w, sgn * s0; hh_tol, hh_maxit)
        end
        for (k, n) in res
            masses[(k, w, sgn, tag)] = n
        end
    end
    gem = Dict{Tuple{Int,Int},Any}()
    ge_ref = nothing
    if ge
        gsolver = merge(solver, (; tol = ge_tol))
        ge_ref = timed(@sprintf("margin (f), zero-shock re-solve, tol %.0e", ge_tol)) do
            ge_masses_safe(sol, 1, 0.0; solver = gsolver)
        end
        ge_ok(ge_ref) || println("    ⚠ the zero-shock re-solve did not converge: (f) is measured against the base solution")
        init = ge_ok(ge_ref) ? ge_ref.sol : sol
        for w in 1:L, sgn in (1, -1)
            gem[(w, sgn)] = timed(@sprintf("margin (f), general equilibrium, κ_%d %s%g", w, sname(sgn), ge_step)) do
                ge_masses_safe(init, w, sgn * ge_step; solver = gsolver)
            end
        end
    end
    ge_n0 = ge_ref !== nothing && ge_ok(ge_ref) ? ge_ref.n : n0
    pay = timed("small-district elasticity, analytic margin (a)") do
        pay_response(sol)
    end
    crossing = timed("crossing shares") do
        crossing_shares(sol, meas)
    end
    # the checks: zero-shock reproduction, the masses against the moments, margin (d) against the moment
    zero_a = p.split ? teacher_counts(hh, Φ, workplace_only_probs(hh, p, gr, sol.Pi_Tr), gr, p) : fill(NaN, L)
    zero_b = teacher_counts(hh, Φ, location_probs(hh, p, gr)[1], gr, p)
    reldiff_n(x) = maximum(abs.(x .- n0) ./ n0)
    tag_s(tag) = tag === :fine ? ge_step : step
    resp(k, w, sgn = 1, tag = :main) = responses(n0, masses[(k, w, sgn, tag)], sgn * tag_s(tag))
    central(k, w, tag) = (log.(masses[(k, w, 1, tag)]) .- log.(masses[(k, w, -1, tag)])) ./ (2 * tag_s(tag))
    gap_a(est) = maximum(abs, [est(w)[i] - pay.workplace[i, w] for w in 1:L, i in 1:L])    # analytic (a) against an estimate of ε_{·,w}
    checks = (; zero_a = reldiff_n(zero_a), zero_b = reldiff_n(zero_b), n_moments = reldiff_n(m.n_teach_work),
                n_pay_response = reldiff_n(pay.n), y_resid = pay.y_resid, floor_mass = pay.floor_mass,
                floor_q_share = pay.floor_q_share, tse_moment = m.teacher_supply_elasticity,
                tse_d = resp(:d, 2).ratio, share_a = resp(:a, 2).ratio / resp(:d, 2).ratio,
                total_fixed = maximum(abs, [resp(k, w).total for k in (:a, :b, :c_lambda_held, :c), w in 1:L
                                            if p.split || k !== :a]),
                analytic_fd = p.split ? gap_a(w -> resp(:a, w).eps) : NaN,
                analytic_central = p.split && two_sided ? gap_a(w -> (resp(:a, w).eps .+ resp(:a, w, -1).eps) ./ 2) : NaN,
                analytic_fine = p.split ? gap_a(w -> central(:a, w, :fine)) : NaN,
                ge_shift = ge_ref !== nothing && ge_ok(ge_ref) ? maximum(abs, log.(ge_ref.n ./ n0)) : NaN,
                ge_converged = ge_ref === nothing ? true : ge_ok(ge_ref) && all(ge_ok, values(gem)),
                crossing_all = crossing === nothing ? NaN : abs(crossing.all.share - m.commute_share_all),
                crossing_teach = crossing === nothing ? NaN : abs(crossing.teach_ms.share - m.commute_share_teach),
                crossing_nonteach = crossing === nothing ? NaN : abs(crossing.nonteach.share - m.commute_share_nonteach),
                crossing_exact = crossing === nothing ? NaN : abs(crossing.teach.share - m.commute_share_teach),
                net_gap = crossing === nothing ? NaN : abs(crossing.teach.net_gap - m.net_commute_gap))
    return (; n0, masses, ge = gem, ge_ref, ge_n0, pay, crossing, checks, step, ge_step, ge_tol, two_sided)
end

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
"A number to three decimals (no negative zero), or `n/a` for NaN (a margin that does not exist for the bundled model)."
f3(x) = isfinite(x) ? @sprintf("%.3f", round(x; digits = 3) + 0.0) : "n/a"

"A vector in scientific notation."
sci(v) = "(" * join((@sprintf("%.2e", x) for x in v), ", ") * ")"

"The step of the table `tag`: `:main` is the audit's `step` (margins (a)–(e)), `:fine` the step of (f), `ge_step`."
tag_step(a, tag) = tag === :fine ? a.ge_step : a.step

"The base masses that margin `k` is measured against: the zero-shock GE re-solve for (f), the base solution's for the rest."
margin_base(a, k) = k === :f ? a.ge_n0 : a.n0

"""
    margin_estimates(a, k, w, tag) -> Vector{NamedTuple}

The difference quotients of margin `k` after the pay shock of workplace `w`, in display order: the forward difference
(+step), the backward (−step) if computed, and the central difference if both were. Each has `kind` (:forward,
:backward, :central), the `step` (signed for the one-sided differences, the half-width for the central one), the
responses `eps`, `ratio` and `total`, and `n`, the masses at +step (at −step for the backward difference). Margin (f)
is at `a.ge_step` whatever the `tag`; the others at the step of `tag`.
"""
function margin_estimates(a, k, w, tag)
    n0 = margin_base(a, k)
    s0 = k === :f ? a.ge_step : tag_step(a, tag)
    masses(sgn) = k === :f ? (haskey(a.ge, (w, sgn)) ? a.ge[(w, sgn)].n : nothing) : get(a.masses, (k, w, sgn, tag), nothing)
    np, nm = masses(1), masses(-1)
    out = NamedTuple[]
    np === nothing || push!(out, (; kind = :forward, step = s0, responses(n0, np, s0)..., n = np))
    nm === nothing || push!(out, (; kind = :backward, step = -s0, responses(n0, nm, -s0)..., n = nm))
    if np !== nothing && nm !== nothing
        eps = (log.(np) .- log.(nm)) ./ (2 * s0)
        push!(out, (; kind = :central, step = s0, eps, ratio = eps[2] - eps[1],
                      total = (log(sum(np)) - log(sum(nm))) / (2 * s0), n = np))
    end
    return out
end

"The GE runs (`a.ge` entries) behind a difference of margin (f) of the kind `kind` at the shocked workplace `w`."
ge_runs(a, w, kind) = [a.ge[(w, sgn)] for sgn in (kind === :forward ? (1,) : kind === :backward ? (-1,) : (1, -1))]

"The step column: the signed step of a one-sided difference, ± the half-width of a central one."
fmt_step(e) = e.kind === :central ? @sprintf("±%.4f", e.step) : @sprintf("%+.4f", e.step)

"The tables of the report as (tag, margin keys): (a)–(e) at `step`, and (a)–(f) (just (a)–(e) without GE) at the step of (f)."
table_specs(a) = [(:main, PE_MARGINS), (:fine, isempty(a.ge) ? PE_MARGINS : (PE_MARGINS..., :f))]

"""
The table of margins `keys` at the steps of `tag`: one row per shocked workplace and kind of difference. In the table at
the step of (f) the household margins show only their central row (all are in the TOML); (f) shows its three rows, the
one-sided ones being the evidence of how far the response is from linear.
"""
function margin_table(io::IO, a, p, tag, keys)
    @printf(io, "  %-44s %3s %8s %9s %9s %9s %9s  %s\n", "margin", "w′", "step", "ε_1,w′", "ε_2,w′", "ε_2−ε_1", "ε_total", "difference")
    for (k, code, text) in margin_labels(p)
        k in keys || continue
        for w in 1:length(a.n0)
            es = margin_estimates(a, k, w, tag)
            tag === :fine && k !== :f && filter!(e -> e.kind === :central, es)
            for (j, e) in enumerate(es)
                bad = k === :f && !all(g -> ge_ok(g), ge_runs(a, w, e.kind))
                @printf(io, "  %-44s %3d %8s %9s %9s %9s %9s  %s%s\n", w == 1 && j == 1 ? code * " " * text : "", w, fmt_step(e),
                        f3(e.eps[1]), f3(e.eps[2]), f3(e.ratio), f3(e.total), e.kind, bad ? "  ⚠ GE not converged" : "")
            end
        end
    end
end

"The analytic derivative of margin (a) (the limit of the finite differences as the step falls), under the split."
function analytic_rows(io::IO, a, p)
    p.split || return
    for ws in 1:length(a.n0)
        ε = a.pay.workplace[:, ws]
        @printf(io, "  %-44s %3d %8s %9s %9s %9s %9s  %s\n", ws == 1 ? "(a)  analytic derivative (step → 0)" : "", ws, "→ 0",
                f3(ε[1]), f3(ε[2]), f3(ε[2] - ε[1]), f3(a.pay.total[ws]), "analytic")
    end
end

"""
The line of the decomposition of the ratio ε_2,w′ − ε_1,w′ at the steps of `tag`: the level at the first margin of `keys`,
then the increment that each later margin adds, from the central difference if there is one, else the forward one.
"""
function decomposition_line(a, p, tag, w, keys)
    parts, prev = String[], NaN
    for (k, code, _) in margin_labels(p)
        k in keys || continue
        es = margin_estimates(a, k, w, tag)
        isempty(es) && continue
        r = last(es).ratio
        isfinite(r) || continue
        push!(parts, isnan(prev) ? @sprintf("%s %.3f", strip(code), r) : @sprintf("%s %+.3f", strip(code), r - prev))
        prev = r
    end
    return join(parts, "   ")
end

"The movement of aggregate `f` (a field of `ge_masses`) of the equilibrium per unit of log κ_w′: the central difference."
ge_per_unit(a, w, f::Symbol) = (getproperty(a.ge[(w, 1)], f) .- getproperty(a.ge[(w, -1)], f)) ./ (2 * a.ge_step)

"The report of `run_audit`'s result `a` as text, with the header (ϑ, options) and the data comparators."
function audit_report(io::IO, a, sol, m, θ, opts; solver)
    (; p, gr) = sol
    L = gr.L
    c = a.checks
    println(io, "spatial_supply_audit.jl  ", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    println(io, "  options: ", join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                                     if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "))
    println(io, "  solver: ", solver)
    println(io, "\nϑ = (", join((@sprintf("%s = %.5f", n, x) for (n, x) in zip(THETA_NAMES, θ)), ", "), ")")
    @printf(io, "  %s%s; the pay of workplace w′ raised by the factor e^step: step %g, forward%s, then step %g, central%s\n",
            p.split ? "residence/workplace split" : "bundled model (workplace = residence)", p.college ? " with college (M3)" : "",
            a.step, a.two_sided ? " (and backward, central)" : "", a.ge_step,
            isempty(a.ge) ? "" : @sprintf("; margin (f) at GE tolerance %.0e", a.ge_tol))

    println(io, "\n(i) Commuting block and base teacher masses")
    @printf(io, "  σν = %.4f, σw = %s, ρw = σw/σν = %s, τk = %.4f;  κ = %s, B = %s%s\n", p.σν,
            p.split ? @sprintf("%.5f", p.σw) : "n/a", p.split ? @sprintf("%.5f", p.σw / p.σν) : "n/a", p.τk, fmtvec(p.κ),
            fmtvec(p.B), p.split ? ", positions: pupils' shares, non-teachers' jobs " * fmtvec(p.jobs) : "")
    @printf(io, "  teachers by workplace n_work = %s (by residence n_res = %s), N_T = %.6f;  M = %s, t = %s, Q = %s\n",
            fmtvec(a.n0), fmtvec(m.n_teach_res), sum(a.n0), fmtvec(sol.M), fmtvec(sol.t), fmtvec(gr.Q))

    last_key = isempty(a.ge) ? "(e)" : "(f)"
    println(io, "\n(a)–(e) Teacher mass by workplace after a pay rise at workplace w′, step ", a.step,
            a.step == SUPPLY_STEP ? " (the step of the moment)" : "", ":  ε_{i,w′} = Δ log n_i / Δ log κ_w′")
    margin_table(io, a, p, :main, PE_MARGINS)
    println(io, "\n(a)–", last_key, " The same at the step ", a.ge_step, " by central differences: the local response (the derivative at the base)")
    margin_table(io, a, p, :fine, table_specs(a)[2][2])
    analytic_rows(io, a, p)
    println(io, "  Each margin adds to the one above; what it does not list stays at the base solution. (a) recomputes π_{w|r} only (nodes")
    println(io, "  keep their (e, s), Λ, occupation maps, education weights, π_r, Φ); (b) the joint π_{rw}; (c₀) (e, s) re-solved with Λ held at")
    println(io, "  its base value (one household sweep); (c) the household fully re-solved so that Λ moves; (c₀) and (c) are integrated with the")
    println(io, "  base maps and weights; (d) the household's own maps and weights (", p.college ? "occupation, education" : "occupation",
            "): the existing `teacher_supply_elasticity`,")
    println(io, "  the ratio column at w′ = 2 of the first table; (e) Φ re-solved; (f) full GE (positions, Q, t, H̃_T, Φ). In (a)–(c) the maps are")
    println(io, "  fixed, so ε_total = 0 and the response reallocates teachers. One-sided rows are differences at ±step, exact to order step; the")
    println(io, "  central row, (log n₊ − log n₋)/(2 step), to order step². The central rows of the second table are the elasticities; the first")
    println(io, "  table is the one-sided difference at step ", a.step, a.step == SUPPLY_STEP ? ", the moment's own." : ".")
    if !isempty(a.ge)
        g0 = a.ge_ref
        @printf(io, "\n  General equilibrium (f): tolerance %.0e (the calibration's %.0e), step %g. The solves start from a zero-shock\n",
                a.ge_tol, solver.tol, a.ge_step)
        if ge_ok(g0)
            @printf(io, "  re-solve of the base solution at that tolerance, which moves the teacher masses by %.1e at most (in log),\n", c.ge_shift)
            @printf(io, "  Δ log M %s, Δ log H̃_T %s, Δt %s, Δ log Q %s:\n", sci(g0.dlogM), sci(g0.dlogHT), sci(g0.dt), sci(g0.dlogQ))
            println(io, "  the base solution is at the fixed point to that accuracy, and the one-sided rows are measured against the re-solve.")
            @printf(io, "  That checks reproduction, not the error of the shocked solves, of the order of tolerance ÷ step = %.0e in ε.\n",
                    a.ge_tol / a.ge_step)
            println(io, "  The central row does not use the reference.")
        else
            println(io, "  ⚠ The zero-shock re-solve failed or did not converge: the one-sided rows are measured against the base solution,")
            @printf(io, "  whose own solver error (tolerance %.0e) is divided by the step %g; use the central row.\n", solver.tol, a.ge_step)
        end
        println(io, "  Movement of the aggregates per unit of log κ_w′ (location 1, 2; central difference):")
        for w in 1:L
            @printf(io, "    w′ = %d: d log M %s, d log H̃_T %s, dt %s, d log Q %s\n", w, fmtvec(ge_per_unit(a, w, :dlogM)),
                    fmtvec(ge_per_unit(a, w, :dlogHT)), fmtvec(ge_per_unit(a, w, :dt)), fmtvec(ge_per_unit(a, w, :dlogQ)))
        end
    end

    println(io, "\nDecomposition of the ratio ε_2,w′ − ε_1,w′: the level of the first margin, then the increment each margin adds (w′ = 2 is the moment)")
    for (tag, keys) in table_specs(a)
        println(io, "  step ", tag_step(a, tag), ", ", tag === :main && !a.two_sided ? "forward" : "central", " differences:")
        for w in 1:L
            println(io, "    w′ = $w:  ", decomposition_line(a, p, tag, w, keys))
        end
    end

    println(io, "\n(g) Own-pay elasticity of teacher supply of a small district (negligible weight inside workplace w) raising its pay")
    σin = p.split ? p.σw : p.σν
    for w in 1:L
        ε = a.pay.small[w]
        @printf(io, "  w = %d:  ε_small = %s  (teachers working in w %.6f;  analytic (a) own ε_%d%d %s;  σ_in = %.5f)\n",
                w, f3(ε), a.pay.n[w], w, w, f3(a.pay.workplace[w, w]), σin)
        @printf(io, "          the σ_in that gives ε_small = %.1f: %.4f;  = %.1f: %.4f\n", LITERATURE_SUPPLY.ransom_sims,
                σin * ε / LITERATURE_SUPPLY.ransom_sims, LITERATURE_SUPPLY.falch, σin * ε / LITERATURE_SUPPLY.falch)
    end
    println(io, "  ε_small(w) = Σ_nodes mass ∫ teach Σ_r π_{rw} μ y_{rw}/(σ_in C_{rw}) / n_w, the mass-weighted mean of π_{rw} μ y/(σ_in C)")
    println(io, "  over teacher nodes and residences divided by the teachers in w; σ_in = ",
            p.split ? "σw, the workplace nest" : "σν (the bundled model's only logit; an approximation)", ".")
    println(io, "  The σ_in that would give a target elasticity holds y/C and the probabilities fixed: first order. Destinations at the")
    println(io, "  consumption floor (C ≤ C_FLOOR) are left out of the sums, here and in the analytic (a); the Checks give what they would carry.")
    @printf(io, "  Literature, single employers: Falch (2010) %.1f, Ransom–Sims (2010) %.1f.\n", LITERATURE_SUPPLY.falch, LITERATURE_SUPPLY.ransom_sims)

    println(io, "\n(h) Crossing between residence and workplace, by direction")
    if a.crossing === nothing
        println(io, "  not defined for the bundled model (the workplace is the residence)")
    else
        @printf(io, "  %-30s %9s %9s %9s %9s    %s\n", "", "P(2 | 1)", "P(1 | 2)", "share", "net gap", "ACS 2018 (block basis): P(2 | 1), P(1 | 2); net gap")
        for (key, label) in ((:teach, "teachers (exact)"), (:nonteach, "non-teaching market workers"), (:all, "all market workers"))
            s, d = getproperty(a.crossing, key), getproperty(ACS_FLOWS, key)
            @printf(io, "  %-30s %9.4f %9.4f %9.4f %+9.4f    %.4f, %.4f; %+.3f\n", label, s.p12, s.p21, s.share, s.net_gap,
                    d[2] / (d[1] + d[2]), d[3] / (d[3] + d[4]), getproperty(ACS_NET_GAP, key))
        end
        t = a.crossing.teach_ms
        @printf(io, "  %-30s %9.4f %9.4f %9.4f %+9.4f    (pseudo-sample row of the teachers, the moment commute_share_teach)\n",
                "teachers (pseudo-sample)", t.p12, t.p21, t.share, t.net_gap)
        @printf(io, "  Targets: crossing share of all workers %.3f, of K–12 teachers %.3f; net commuting gap of teachers %+.3f\n",
                target_value(:commute_share_all), target_value(:commute_share_teach), target_value(:net_commute_gap))
        println(io, "  (ACS 2018 K–12 +0.120, public K–12 +0.108 on three classes). The ACS flows are lower bounds (a POWPUMA is a county or more):")
        println(io, "  the crossing within a county is invisible. Net gap = log(suburb/city) of residents less that of jobs; the data's gaps are on")
        println(io, "  three classes (city, suburb, other), the model's on its two locations.")
    end

    println(io, "\nChecks")
    @printf(io, "  base masses against the moments' n_work: max relative difference %.1e;  against the pay-response integral %.1e\n",
            c.n_moments, c.n_pay_response)
    @printf(io, "  node income against consumption(): max |C − (y − (1 + τe) e − m)| = %.1e\n", c.y_resid)
    @printf(io, "  destinations at the consumption floor, left out of (g) and the analytic (a): teacher mass share %.1e,\n", c.floor_mass)
    @printf(io, "  and the share of Σ π q that they would carry with the linearisation's slope %.1e\n", c.floor_q_share)
    @printf(io, "  zero-shock reproduction of the base masses: %s(b) %.1e\n", p.split ? @sprintf("(a) %.1e, ", c.zero_a) : "", c.zero_b)
    @printf(io, "  margin (d), w′ = 2, ratio %s against the moment teacher_supply_elasticity %s%s\n", @sprintf("%.5f", c.tse_d),
            isfinite(c.tse_moment) ? @sprintf("%.5f", c.tse_moment) : "n/a (bundled)",
            isfinite(c.tse_moment) && a.step == SUPPLY_STEP ? @sprintf(" (difference %.1e)", c.tse_d - c.tse_moment) :
            (isfinite(c.tse_moment) ? " (the moment's step is $SUPPLY_STEP, this one is $(a.step))" : ""))
    p.split && @printf(io, "  margin (a) is %.0f%% of (d)'s ratio at w′ = 2;  largest |ε_total| in (a)–(c): %.1e\n", 100 * c.share_a, c.total_fixed)
    p.split && @printf(io, "  analytic against margin (a), max |Δε|: %.1e against the one-sided difference at step %g (of order step × curvature)%s;\n",
                       c.analytic_fd, a.step, isfinite(c.analytic_central) ?
                       @sprintf(", %.1e against its central difference", c.analytic_central) : "")
    p.split && @printf(io, "  %.1e against the central difference at step %g\n", c.analytic_fine, a.ge_step)
    isempty(a.ge) || @printf(io, "  margin (f): zero-shock re-solve at tolerance %.0e moves the base masses by %s (log); all GE solves converged: %s\n",
                             a.ge_tol, isfinite(c.ge_shift) ? @sprintf("%.1e", c.ge_shift) : "n/a", c.ge_converged ? "yes" : "NO")
    p.split && @printf(io, "  crossing shares against the moments: all %.1e, teachers' pseudo-sample %.1e (exact %.1e), non-teachers %.1e;\n",
                       c.crossing_all, c.crossing_teach, c.crossing_exact, c.crossing_nonteach)
    p.split && @printf(io, "  teachers' net gap %.1e\n", c.net_gap)
    return nothing
end

"""
The audit's numbers as a nested dictionary for TOML: header, base, `margin.<key>.shock_<w′>` with the forward difference at
the top (`step`, `n`, `N_T`, `eps_1`, `eps_2`, `ratio`, `eps_total`) at `--step`, `.minus` (backward) and `.central` with
--two-sided, and `.fine` for the same at `--ge-step` (forward, `.minus`, `.central`: the local response is
`.fine.central`). Margin (f) has only that step, so its top level and `.fine` coincide; it also carries the movement of
the aggregates (`per_unit`: per unit of log κ_w′, central) and `margin.f.reference` the zero-shock re-solve; then the rest.
"""
function audit_dict(a, sol, m, θ, opts; solver)
    (; p, gr) = sol
    L = gr.L
    d = Dict{String,Any}()
    d["meta"] = Dict{String,Any}("script" => "spatial_supply_audit.jl", "created" => Libc.strftime("%Y-%m-%d %H:%M:%S", time()),
                                 "theta" => opts["theta"], "step" => a.step, "ge_step" => a.ge_step, "ge_tol" => a.ge_tol,
                                 "two_sided" => a.two_sided, "split" => p.split, "college" => p.college, "solver" => string(solver),
                                 "margin_keys" => "a, b, c_lambda_held (c0: e, s re-solved with Λ held), " *
                                                  "c (the full household re-solve; the specification's (c)), d, e, f",
                                 "options" => join(("$k = $v" for (k, v) in sort!(collect(opts); by = first)
                                                    if !(v === nothing || v == false || (v isa Vector && isempty(v)))), ", "))
    d["theta"] = Dict{String,Any}(string(n) => x for (n, x) in zip(THETA_NAMES, θ))
    d["base"] = Dict{String,Any}("sigma_nu" => p.σν, "sigma_w" => p.split ? p.σw : NaN, "rho_w" => p.split ? p.σw / p.σν : NaN,
                                 "tau_k" => p.τk, "kappa" => p.κ, "n_work" => a.n0, "n_res" => m.n_teach_res, "N_T" => sum(a.n0),
                                 "M" => sol.M, "t" => sol.t, "Q" => gr.Q)
    entry(e) = Dict{String,Any}("step" => e.step, "n" => e.n, "N_T" => sum(e.n), "eps_1" => e.eps[1], "eps_2" => e.eps[2],
                                "ratio" => e.ratio, "eps_total" => e.total)
    function estimates(es)                                   # forward at the top, backward as `minus`, `central`
        out = Dict{String,Any}()
        for e in es
            e.kind === :forward  && merge!(out, entry(e))
            e.kind === :backward && (out["minus"] = entry(e))
            e.kind === :central  && (out["central"] = Dict{String,Any}("step" => e.step, "eps_1" => e.eps[1], "eps_2" => e.eps[2],
                                                                      "ratio" => e.ratio, "eps_total" => e.total))
        end
        return out
    end
    ge_fields(g) = Dict{String,Any}("converged" => g.converged, "dlogM" => g.dlogM, "dlogHT" => g.dlogHT, "dt" => g.dt,
                                    "dlogQ" => g.dlogQ)
    margins = Dict{String,Any}()
    for (k, code, text) in margin_labels(p)
        mk = Dict{String,Any}("label" => strip(code) * " " * text)
        for w in 1:L
            if k === :f
                es = margin_estimates(a, k, w, :fine)
                isempty(es) && continue
                e = estimates(es)
                merge!(e, ge_fields(a.ge[(w, 1)]))
                merge!(e["minus"], ge_fields(a.ge[(w, -1)]))
                e["per_unit"] = Dict{String,Any}("difference" => "central", "dlogM" => ge_per_unit(a, w, :dlogM),
                                                 "dlogHT" => ge_per_unit(a, w, :dlogHT), "dt" => ge_per_unit(a, w, :dt),
                                                 "dlogQ" => ge_per_unit(a, w, :dlogQ))
                e["fine"] = deepcopy(e)
            else
                e = estimates(margin_estimates(a, k, w, :main))
                e["fine"] = estimates(margin_estimates(a, k, w, :fine))
            end
            mk["shock_$w"] = e
        end
        if k === :f && a.ge_ref !== nothing
            g0 = a.ge_ref
            mk["reference"] = Dict{String,Any}("converged" => g0.converged, "n" => g0.n, "N_T" => sum(g0.n), "ge_tol" => a.ge_tol,
                                               "max_abs_dlogn_vs_base" => a.checks.ge_shift, "dlogM" => g0.dlogM, "dlogHT" => g0.dlogHT,
                                               "dt" => g0.dt, "dlogQ" => g0.dlogQ)
        end
        length(mk) > 1 && (margins[string(k)] = mk)
    end
    d["margin"] = margins
    d["analytic_a"] = Dict{String,Any}("eps_$(wo)_$(ws)" => a.pay.workplace[wo, ws] for wo in 1:L, ws in 1:L)
    for ws in 1:L
        d["analytic_a"]["eps_total_$ws"] = a.pay.total[ws]
    end
    d["small_district"] = Dict{String,Any}("elasticity" => a.pay.small, "n" => a.pay.n, "sigma_in" => p.split ? p.σw : p.σν,
                                           "falch" => LITERATURE_SUPPLY.falch, "ransom_sims" => LITERATURE_SUPPLY.ransom_sims)
    if a.crossing !== nothing
        d["crossing"] = Dict{String,Any}(string(k) => Dict{String,Any}("p12" => v.p12, "p21" => v.p21, "share" => v.share,
                                                                       "net_gap" => v.net_gap) for (k, v) in pairs(a.crossing))
    end
    d["checks"] = Dict{String,Any}(string(k) => v for (k, v) in pairs(a.checks))
    return d
end

function main(opts)
    ref = frozen_reference(something(opts["reference"], REFERENCE_FILE), SOLVER)
    (; base, meas, floor_tol) = problem_setup(opts, ref)
    meas = with_measurement(meas; pooled_p9010 = false)       # not used here (a second pseudo-sample)
    opts["step"] > 0 || error("--step must be positive, got $(opts["step"])")
    opts["ge-step"] > 0 || error("--ge-step must be positive, got $(opts["ge-step"])")
    opts["ge-tol"] > 0 || error("--ge-tol must be positive, got $(opts["ge-tol"])")
    θ = read_run_theta(opts["theta"])
    for (name, v) in opts["fix"]
        θ[theta_index(name)] = v
    end
    model = SpatialModel(base, Symbol[]; meas, floor_tol)
    sol, m = timed("equilibrium and moments") do
        solve_theta(model, θ)
    end
    sol.gr.L == 2 || error("the audit is for the city/suburb economy (L = 2)")
    a = run_audit(sol, m, meas; step = opts["step"], ge_step = opts["ge-step"], ge_tol = opts["ge-tol"], ge = !opts["no-ge"],
                  two_sided = opts["two-sided"], solver = model.solver)
    buf = IOBuffer()
    audit_report(buf, a, sol, m, θ, opts; solver = model.solver)
    txt = String(take!(buf))
    println("\n", txt)
    if opts["out"] !== nothing
        mkpath(opts["out"])
        write(joinpath(opts["out"], "audit.txt"), txt)
        open(joinpath(opts["out"], "audit.toml"), "w") do f
            TOML.print(f, audit_dict(a, sol, m, θ, opts; solver = model.solver); sorted = true)
        end
        println("written to ", opts["out"])
    end
    return a
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(OPTS)
end
