# Generate plots for presentation slides:
#=
# Select type of parameterization (from following list):
    # (1) benchmark
    # (2) counter_# (check # of counterfactual)
    # (3) high_beta
    # (4) low_beta
    # (5) high_beta_high_sigma
    # (6) low_beta_low_sigma
    paramname = "counter_4"
    pathname = string("./plots/counterfactuals/",paramname,"/")
# Select backend:
pyplot()
# (1) 'year' = 1970
# (1.a.0) Distribution of teachers' abilities (steady-state)
plt1 = plot(a_grid[2:end-2],f_T[2:end-2,iHH,1],label=string(year))
plt2 = plot(a_grid[2:end-2],f_T[2:end-2,iHH,2],label=string(year))
# (1.a.1) Teachers' human capital (steady-state)
plt1_1 = plot(a_grid[2:end-2],h_T[2:end-2,iHH,1],label=string(year))
plt2_1 = plot(a_grid[2:end-2],h_T[2:end-2,iHH,2],label=string(year))
# (1.a.2) Distribution of teachers' human capital (steady-state)
plt1_2 = plot(h_T[2:end-2,iHH,1],f_T[2:end-2,iHH,1],label=string(year))
plt2_2 = plot(h_T[2:end-2,iHH,2],f_T[2:end-2,iHH,2],label=string(year))
# (1.b) Scatter plot of labor market barriers:
plt3 = scatter(τ_w[:,1],label=string(year),markershape=:circle)
# (1.b.1) Create plot for 1970 only:
plt3_1_70 = scatter(τ_w[:,1],label=string(year),markershape=:circle,ylims=(-1.5,1),legend=:left)
plot!(plt3_1_70,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,subplot=1)
filename = string(pathname,"tau_w_women_70.eps")
savefig(plt3_1_70,filename)
# (1.b.2) Scatter plot of occupational productivities:
plt3_2_70 = scatter(a_by_occ,label=string(year),markershape=:circle,legend=:topleft,ylims=(0,1.25))
plot!(plt3_2_70,ylabel=L"A_o",title="Occupational Productivities",grid=false,subplot=1)
filename = string(pathname,"A_men_70.eps")
savefig(plt3_2_70,filename)
# (1.c) Distribution of 'e_T'
plt4 = plot(a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))

# (1.d) Law of motion for 'H_T'
plt6 = plot(H_grid./H_grid[iHH],H_grid./H_grid[iHH],linewidth=.5,linestyle=:dash)
plot!(plt6,H_grid./H_grid[iHH],HH_T./HH_T[iHH],linewidth=1,linestyle=:solid)
# (1.e) Occupational threshold for occupation 1 ("Executives,.."):
# (1.e.1) Women
plt7w = plot(a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year),linestyle=:dash)
# (1.e.2) Men
plt7m = plot(a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year),linestyle=:dash)

# (2) 'year' = 1990
# (2.a.0) 
plot!(plt1,a_grid[2:end-2],f_T[2:end-2,iHH,1],label=string(year))
plot!(plt2,a_grid[2:end-2],f_T[2:end-2,iHH,2],label=string(year))
# (2.a.1)
plot!(plt1_1,a_grid[2:end-2],h_T[2:end-2,iHH,1],label=string(year))
plot!(plt2_1,a_grid[2:end-2],h_T[2:end-2,iHH,2],label=string(year))
# (2.a.2)
plot!(plt1_2,h_T[2:end-2,iHH,1],f_T[2:end-2,iHH,1],label=string(year))
plot!(plt2_2,h_T[2:end-2,iHH,2],f_T[2:end-2,iHH,2],label=string(year))
# (2.b) 
scatter!(plt3,τ_w[:,1],label=string(year),markershape=:diamond)
# (2.b.1) Create plot for 1970-90 only:
plt3_1_70_90 = scatter!(plt3_1_70,τ_w[:,1],label=string(year),markershape=:diamond,ylims=(-1.5,1))
filename = string(pathname,"tau_w_women_70_90.eps")
savefig(plt3_1_70_90,filename)
# (2.b.2) Create plot for 1970-90 only:
plt3_2_70_90 = scatter!(plt3_2_70,a_by_occ,label=string(year),markershape=:diamond,ylims=(0,1.25))
filename = string(pathname,"A_men_70_90.eps")
savefig(plt3_2_70_90,filename)
# (2.c)
plot!(plt4,a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))

# (2.e.1)
plot!(plt7w,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year),linestyle=:dash)
# (2.e.1)
plot!(plt7m,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year),linestyle=:dash)

# (3) 'year' = 2010
# (3.a.0)
plot!(plt1,a_grid[2:end-2],f_T[2:end-2,iHH,1],label=string(year))
plot!(plt2,a_grid[2:end-2],f_T[2:end-2,iHH,2],label=string(year))
# (3.a.1)
plot!(plt1_1,a_grid[2:end-2],h_T[2:end-2,iHH,1],label=string(year))
plot!(plt2_1,a_grid[2:end-2],h_T[2:end-2,iHH,2],label=string(year))
# (3.a.2)
plot!(plt1_2,h_T[2:end-2,iHH,1],f_T[2:end-2,iHH,1],label=string(year))
plot!(plt2_2,h_T[2:end-2,iHH,2],f_T[2:end-2,iHH,2],label=string(year))
# (3.b) 
scatter!(plt3,τ_w[:,1],label=string(year),markershape=:xcross)
# (2.b.1) Create plot for 1970-2010 only:
plt3_1_70_10 = scatter!(plt3_1_70_90,τ_w[:,1],label=string(year),markershape=:xcross)
filename = string(pathname,"tau_w_women_70_10.eps")
savefig(plt3_1_70_10,filename)
# (2.b.2) Create plot for 1970-2010 only:
plt3_2_70_10 = scatter!(plt3_2_70_90,a_by_occ,label=string(year),markershape=:xcross,ylims=(0,1.25))
filename = string(pathname,"A_men_70_10.eps")
savefig(plt3_2_70_10,filename)
# (3.c)
plot!(plt4,a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))
# plot!(plt5,a_grid[2:end-2],e_T[2:end-2,iHH,2],label=string(year))
# (3.e.1)
plot!(plt7w,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year),linestyle=:dash)
# (3.e.2)
plot!(plt7m,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year),linestyle=:dash)

# (4) Add axis labels and titles:
# (4.a.0)
plot!(plt1,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Abilities",ylims=(-0.001,.105),subplot=1)
plot!(plt2,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Abilities",ylims=(-.001,.065),subplot=1)
# (4.a.1)
plot!(plt1_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
plot!(plt2_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
# (4.a.2)
plot!(plt1_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Human Capital",subplot=1)
plot!(plt2_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Human Capital",subplot=1)
# (4.b)
plot!(plt3,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,legend=:left,subplot=1)

plot!(plt3_1_70_90,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,subplot=1)
plot!(plt3_1_70_10,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,subplot=1)
# (4.c)
plot!(plt4,xlabel="Idiosyncratic Ability",ylabel=L"e_T",grid=false,title="Teachers' Human Capital Investment",legend=:left,subplot=1)
# plot!(plt5,xlabel="Idiosyncratic Ability",ylabel=L"e_T",grid=false,title="Male Teachers' Human Capital Investment",legend=:left,subplot=1)
plot!(plt6,xlabel=L"\widetilde{H}_T",ylabel=L"\widetilde{H}_T^'",grid=false,title="Law of Motion for Teachers' Human Capital",legend=false)
# (4.e.1)
plot!(plt7w,xlabel="Idiosyncratic Ability in Teaching",ylabel=L"\bar{a}_1(a)",grid=false,title="Occupational Threshold (Women)",legend=:topleft)
# (4.e.2)
plot!(plt7m,xlabel="Idiosyncratic Ability in Teaching",ylabel=L"\bar{a}_1(a)",grid=false,title="Occupational Threshold (Men)",legend=:topleft)

# (5) Save figures in EPS file format:
savefig(plt1,string(pathname,"fT_women_steadystate.eps"))
savefig(plt2,string(pathname,"fT_men_steadystate.eps"))
savefig(plt1_1,string(pathname,"hT_women_steadystate.eps"))
savefig(plt2_1,string(pathname,"hT_men_steadystate.eps"))
savefig(plt1_2,string(pathname,"hT_fT_women_steadystate.eps"))
savefig(plt2_2,string(pathname,"hT_fT_men_steadystate.eps"))
savefig(plt3,string(pathname,"tau_w_women.eps"))
savefig(plt3_1_70_90,string(pathname,"tau_w_women_70_90.eps"))
savefig(plt3_1_70_10,string(pathname,"tau_w_women_70_10.eps"))
savefig(plt4,string(pathname,"eT_steadystate.eps"))
savefig(plt6,string(pathname,"LoM.eps"))
savefig(plt7w,string(pathname,"a_O_women.eps"))
savefig(plt7m,string(pathname,"a_O_men.eps"))

=#



#=
    plots_counter.jl
    ────────────────────────
    Plotting script for counterfactual exercises of the teacher occupational
    choice model.  Loads saved JLD files from the benchmark and each
    counterfactual, then generates PNG plots:

      (A) Within-counterfactual plots  (analogous to plots_benchmark.jl)
      (B) Baseline-vs-counterfactual comparison plots for each year
      (C) Cross-counterfactual summary plots

    Expected working directory:  ./julia/codes_w_exo_wage/

    Plots are saved to:
        ./plots/<paramname>/counter_<k>/        — per-counterfactual
        ./plots/<paramname>/counter_comparison/  — cross-counterfactual
=#

using Plots, JLD, CSV, DataFrames, LaTeXStrings
import XLSX

# ═══════════════════════════════════════════════════════════════════
#  1. CONFIGURATION
# ═══════════════════════════════════════════════════════════════════

PARAMNAME = "benchmark"
BASEPATH  = string("./parameterization/", PARAMNAME, "/")

# Counterfactual definitions:
#   id        — integer label (1–6)
#   years     — which years have JLD output
#   shortname — for legends / titles
#   desc      — one-line description
CF_SPECS = [
    (id=1, years=[1990, 2010],       shortname="CF1: 1970 barriers",
     desc="Women's barriers frozen at 1970 levels"),
    (id=2, years=[1970, 1990, 2010], shortname="CF2: No barriers",
     desc="All barriers against women eliminated"),
    (id=3, years=[1970, 1990, 2010], shortname="CF3: No barriers ex HP",
     desc="Barriers eliminated except home production"),
    (id=4, years=[1970, 1990, 2010], shortname="CF4: No barriers ex HP (1970 HP)",
     desc="Barriers eliminated; HP barrier frozen at 1970"),
    (id=5, years=[1990, 2010],       shortname="CF5: 1970 κ + barriers ex HP",
     desc="κ frozen at 1970; non-HP barriers frozen at 1970"),
    (id=6, years=[1990, 2010],       shortname="CF6: Fixed H̃_T",
     desc="Aggregate teaching HC fixed at 1970 level"),
]

ALL_YEARS = [1970, 1990, 2010]

# ── Visual style ─────────────────────────────────────────────────
# Colorblind-safe palette (Okabe–Ito)
CB_BLUE   = RGB(0/255, 114/255, 178/255)
CB_ORANGE = RGB(230/255, 159/255, 0/255)
CB_GREEN  = RGB(0/255, 158/255, 115/255)
CB_RED    = RGB(213/255, 94/255, 0/255)
CB_PURPLE = RGB(204/255, 121/255, 167/255)
CB_CYAN   = RGB(86/255, 180/255, 233/255)
CB_YELLOW = RGB(240/255, 228/255, 66/255)

# Year styling (shared with benchmark)
YEAR_COLORS = Dict(1970 => CB_BLUE, 1990 => CB_ORANGE, 2010 => CB_GREEN)
YEAR_STYLES = Dict(1970 => :solid, 1990 => :dash, 2010 => :dot)
YEAR_MARKERS = Dict(1970 => :circle, 1990 => :diamond, 2010 => :xcross)

# Comparison styling: baseline vs counterfactual
COMP_BASELINE = (color=CB_BLUE, linestyle=:solid, linewidth=2.0)
COMP_COUNTER  = (color=CB_RED,  linestyle=:dash,  linewidth=2.0)

# Cross-counterfactual colours (for CF1–CF6)
CF_COLORS = [CB_BLUE, CB_ORANGE, CB_GREEN, CB_RED, CB_PURPLE, CB_CYAN]

ε = 1e-12  # avoid log(0)

gr()
default(;
    framestyle     = :box,
    grid           = false,
    legend         = :best,
    legendfontsize = 7,
    guidefontsize  = 9,
    titlefontsize  = 10,
    tickfontsize   = 7,
    dpi            = 200,
)

SZ_NORMAL = (620, 400)
SZ_WIDE   = (960, 480)
SZ_BAR    = (560, 400)

# ═══════════════════════════════════════════════════════════════════
#  2. LOAD DATA
# ═══════════════════════════════════════════════════════════════════

# ── 2a. Benchmark (main model solution) ──────────────────────────
println("Loading benchmark parameterizations …")
BENCH = Dict{Int, Dict{String, Any}}()
for yr in ALL_YEARS
    fname = string(BASEPATH, "previousParameterization", yr, ".jld")
    if isfile(fname)
        BENCH[yr] = load(fname)
        println("  Loaded benchmark for $yr")
    else
        @warn "Benchmark file not found: $fname"
    end
end

# ── 2b. Benchmark moments CSVs ──────────────────────────────────
BENCH_MOM = Dict{Int, DataFrame}()
for yr in ALL_YEARS
    csvname = string("./results/", PARAMNAME, "/moments", yr, ".csv")
    if isfile(csvname)
        BENCH_MOM[yr] = DataFrame(CSV.File(csvname))
    end
end

# ── 2c. Counterfactual JLD files ────────────────────────────────
# CF_DATA[cf_id][year] = Dict of loaded JLD contents
CF_DATA = Dict{Int, Dict{Int, Dict{String, Any}}}()
for spec in CF_SPECS
    CF_DATA[spec.id] = Dict{Int, Dict{String, Any}}()
    for yr in spec.years
        fname = string(BASEPATH, "counter_", spec.id, "_", yr, ".jld")
        if isfile(fname)
            CF_DATA[spec.id][yr] = load(fname)
            println("  Loaded counter_$(spec.id) for $yr")
        else
            @warn "Counterfactual file not found: $fname"
        end
    end
end

# ── 2d. Counterfactual moments CSVs ─────────────────────────────
CF_MOM = Dict{Int, DataFrame}()
for spec in CF_SPECS
    csvname = string("./results/", PARAMNAME, "/counter_", spec.id, "_moments.csv")
    if isfile(csvname)
        CF_MOM[spec.id] = DataFrame(CSV.File(csvname))
    end
end

# ── 2e. Occupation labels from Excel ────────────────────────────
XLSX_PATH = "../../data/LaborMarketData/wages_occ_shares_v2.xlsx"
OCC_LABELS_FULL = if isfile(XLSX_PATH)
    xs  = XLSX.readxlsx(XLSX_PATH)
    tab = xs["moments_shares"]
    string.(vec(tab["A30:A49"]))
else
    ["Occ $i" for i in 1:20]
end
OCC_LABELS = [length(s) > 18 ? s[1:17] * "…" : s for s in OCC_LABELS_FULL]
N_OCC = length(OCC_LABELS)

# ═══════════════════════════════════════════════════════════════════
#  3. HELPERS
# ═══════════════════════════════════════════════════════════════════

"Year-indexed style kwargs."
ystyle(yr) = (color=YEAR_COLORS[yr], linestyle=YEAR_STYLES[yr], linewidth=1.8)

"Reconstruct e_T from saved parameters."
function compute_e_T(d::Dict)
    η_val = d["η"]; γ_val = d["γ"]; κ_val = d["κ"]
    t_val = d["t"]; hT = d["h_T"]; iH = d["iHH"]
    na = size(hT, 1)
    # h_T may be 3D (n_a × n_H × n_G) or we use iH
    nG = size(hT, 3)
    eT = zeros(na, nG)
    for iG in 1:nG, ia in 1:na
        eT[ia, iG] = η_val * (1 - t_val[iH, 1]) * κ_val * γ_val * hT[ia, iH, iG]^γ_val
    end
    return eT
end

"Extract teaching shares.  For counterfactuals mass_T is a vector [fem, male]."
function get_mass_T_bench(yr::Int)
    haskey(BENCH_MOM, yr) || return [NaN, NaN]
    df = BENCH_MOM[yr]
    r  = nrow(df)
    return [df[r, :share_teachers_female], df[r, :share_teachers_male]]
end

function get_mass_T_cf(cf_id::Int, yr::Int)
    haskey(CF_DATA, cf_id) || return [NaN, NaN]
    # If this year has no CF data (e.g. 1970 for CF1/5/6 which equal benchmark),
    # fall back to the benchmark moments for that year.
    haskey(CF_DATA[cf_id], yr) || return get_mass_T_bench(yr)
    d = CF_DATA[cf_id][yr]
    # Counterfactuals save mass_T as a vector (length n_G)
    return d["mass_T"]
end

"Get a_grid and trimming indices from benchmark (shared across all solutions)."
function get_grid(d::Dict)
    a_grid = d["a_grid"]
    n_a    = length(a_grid)
    trim   = 3:n_a-3
    return a_grid, n_a, trim
end

"Create a scatter plot pre-configured for occupation-label x-axes."
function occ_scatter(; title="", ylabel="", ylims_val=nothing)
    kw = Dict{Symbol,Any}(
        :title          => title,
        :ylabel         => ylabel,
        :xlabel         => "",
        :legend         => :topleft,
        :size           => SZ_WIDE,
        :bottom_margin  => 22Plots.mm,
        :left_margin    => 6Plots.mm,
        :right_margin   => 4Plots.mm,
    )
    if ylims_val !== nothing
        kw[:ylims] = ylims_val
    end
    return plot(; kw...)
end

# Reference grid from benchmark 1970:
A_GRID, N_A, TRIM = get_grid(BENCH[1970])
IHH = BENCH[1970]["iHH"]

# ═══════════════════════════════════════════════════════════════════
#  PART A:  WITHIN-COUNTERFACTUAL PLOTS
#           (analogous to plots_benchmark.jl, one set per CF)
# ═══════════════════════════════════════════════════════════════════

for spec in CF_SPECS
    cf_id = spec.id
    cf_years = sort(collect(keys(CF_DATA[cf_id])))
    isempty(cf_years) && continue

    PLOTPATH = string("./plots/", PARAMNAME, "/counter_", cf_id, "/")
    mkpath(PLOTPATH)
    println("\n── Generating plots for Counter $(cf_id): $(spec.desc) ──")

    # Gather data dicts for available years:
    DD = Dict(yr => CF_DATA[cf_id][yr] for yr in cf_years)

    # For CFs where 1970 == benchmark (CF1, CF5, CF6), inject benchmark 1970 data
    # so that within-CF plots include the 1970 reference line:
    if 1970 ∉ cf_years && haskey(BENCH, 1970)
        DD[1970] = BENCH[1970]
        cf_years = vcat([1970], cf_years)
    end

    # Use a_grid from CF (should match benchmark):
    a_grid_cf = DD[cf_years[1]]["a_grid"]
    n_a_cf    = length(a_grid_cf)
    trim_cf   = 3:n_a_cf-3

    # ── A1. Distribution of teachers' abilities (f_T) ────────────
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$(spec.shortname): $glabel Teachers' Ability Dist.",
            xlabel = "Idiosyncratic Ability", ylabel = "Density",
            legend = :topright, size = SZ_NORMAL)
        for yr in cf_years
            fT = DD[yr]["f_T"]
            iH = DD[yr]["iHH"]
            plot!(p, a_grid_cf[trim_cf], fT[trim_cf, iH, ig];
                  label=string(yr), ystyle(yr)...)
        end
        savefig(p, string(PLOTPATH, "fT_$(lowercase(glabel))_steadystate.png"))
    end

    # ── A2. Teachers' human capital vs ability (h_T) ─────────────
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$(spec.shortname): $glabel Teachers' Human Capital",
            xlabel = "Idiosyncratic Ability", ylabel = "Human Capital",
            legend = :topleft, size = SZ_NORMAL)
        for yr in cf_years
            hT = DD[yr]["h_T"]
            iH = DD[yr]["iHH"]
            plot!(p, a_grid_cf[trim_cf], hT[trim_cf, iH, ig];
                  label=string(yr), ystyle(yr)...)
        end
        savefig(p, string(PLOTPATH, "hT_$(lowercase(glabel))_steadystate.png"))
    end

    # ── A3. Distribution of teachers' human capital ──────────────
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$(spec.shortname): $glabel Teachers' HC Dist.",
            xlabel = "Human Capital", ylabel = "Density",
            legend = :topright, size = SZ_NORMAL)
        for yr in cf_years
            hT = DD[yr]["h_T"]
            fT = DD[yr]["f_T"]
            iH = DD[yr]["iHH"]
            plot!(p, hT[trim_cf, iH, ig], fT[trim_cf, iH, ig];
                  label=string(yr), ystyle(yr)...)
        end
        savefig(p, string(PLOTPATH, "hT_fT_$(lowercase(glabel))_steadystate.png"))
    end

    # ── A4. Labor market barriers (τ_w) ──────────────────────────
    if length(cf_years) >= 1
        # Sort by first year's τ_w for women
        τ_ref = DD[cf_years[1]]["τ_w"][:, 1]
        sidx  = sortperm(τ_ref)
        labels = OCC_LABELS[sidx]

        p = occ_scatter(; title="$(spec.shortname): Barriers Against Women",
                          ylabel=L"\tau_w")
        for yr in cf_years
            τ_w_yr = DD[yr]["τ_w"][:, 1]
            scatter!(p, 1:N_OCC, τ_w_yr[sidx];
                     label=string(yr), color=YEAR_COLORS[yr],
                     markershape=YEAR_MARKERS[yr], markersize=4.5, msw=0.4)
        end
        xticks!(p, 1:N_OCC, labels; rotation=50, tickfontsize=6)
        savefig(p, string(PLOTPATH, "tau_w_women.png"))
    end

    # ── A5. Teachers' HC investment (e_T) ────────────────────────
    p = plot(;
        title  = "$(spec.shortname): Teachers' HC Investment",
        xlabel = "Idiosyncratic Ability", ylabel = L"e_T",
        legend = :topleft, size = SZ_NORMAL)
    for yr in cf_years
        eT = compute_e_T(DD[yr])
        plot!(p, a_grid_cf[trim_cf], eT[trim_cf, 1];
              label=string(yr), ystyle(yr)...)
    end
    savefig(p, string(PLOTPATH, "eT_steadystate.png"))

    # ── A6. Occupational thresholds ──────────────────────────────
    if haskey(DD[cf_years[1]], "a_T_thresh") || haskey(DD[cf_years[1]], "a_O_thresh")
        # Use a_O_thresh if available (saved under that key in CFs that have it)
        thresh_key = haskey(DD[cf_years[1]], "a_O_thresh") ? "a_O_thresh" : "a_T_thresh"
        trim_th = 3:n_a_cf-2
        for (ig, glabel) in enumerate(["Women", "Men"])
            p = plot(;
                title  = "$(spec.shortname): Occ. Threshold ($glabel)",
                xlabel = "Ability in Teaching", ylabel = "ā₁(a)",
                legend = :topleft, size = SZ_NORMAL)
            for yr in cf_years
                thresh = DD[yr][thresh_key]
                iH = DD[yr]["iHH"]
                plot!(p, a_grid_cf[trim_th], thresh[trim_th, iH, ig, 1];
                      label=string(yr), ystyle(yr)...)
            end
            savefig(p, string(PLOTPATH, "a_O_$(lowercase(glabel)).png"))
        end
    end

    # ── A7. Teaching shares over time (bar chart) ────────────────
    if length(cf_years) >= 2
        fem = [get_mass_T_cf(cf_id, yr)[1] for yr in cf_years]
        mal = [get_mass_T_cf(cf_id, yr)[2] for yr in cf_years]
        n   = length(cf_years)
        bw  = 0.35
        xt  = collect(1:n)

        p = plot(;
            title  = "$(spec.shortname): Teaching Shares",
            ylabel = "Fraction", legend = :topleft,
            size = SZ_BAR, xlims = (0.3, n + 0.7),
            xticks = (xt, string.(cf_years)))
        bar!(p, xt .- bw/2, fem; bar_width=bw, color=CB_BLUE,  alpha=0.85, label="Female")
        bar!(p, xt .+ bw/2, mal; bar_width=bw, color=CB_ORANGE, alpha=0.85, label="Male")
        savefig(p, string(PLOTPATH, "mass_T_over_time.png"))
    end

    # ── A8. CDF of teachers' abilities ───────────────────────────
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$(spec.shortname): CDF $glabel Teachers' Abilities",
            xlabel = "Idiosyncratic Ability", ylabel = "Cumulative Probability",
            legend = :bottomright, size = SZ_NORMAL)
        for yr in cf_years
            fT  = DD[yr]["f_T"][trim_cf, DD[yr]["iHH"], ig]
            da  = diff(a_grid_cf[trim_cf])
            cdf_vals = zeros(length(trim_cf))
            for k in 2:length(trim_cf)
                cdf_vals[k] = cdf_vals[k-1] + 0.5 * (fT[k] + fT[k-1]) * da[k-1]
            end
            cdf_vals ./= max(cdf_vals[end], ε)
            plot!(p, a_grid_cf[trim_cf], cdf_vals; label=string(yr), ystyle(yr)...)
        end
        savefig(p, string(PLOTPATH, "cdf_fT_$(lowercase(glabel)).png"))
    end

    # ── A9. Unnormalised f_T · m_T ───────────────────────────────
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$(spec.shortname): Non-normalised Distribution of $glabel Teachers",
            xlabel = "Teaching Ability", ylabel = L"f_T \cdot m_T",
            legend = :topright, size = SZ_NORMAL)
        for yr in cf_years
            fT   = DD[yr]["f_T"]
            iH   = DD[yr]["iHH"]
            mT   = get_mass_T_cf(cf_id, yr)[ig]
            plot!(p, a_grid_cf[trim_cf], fT[trim_cf, iH, ig] .* mT;
                  label=string(yr), ystyle(yr)...)
        end
        savefig(p, string(PLOTPATH, "fT_unnorm_$(lowercase(glabel)).png"))
    end

    # ── A10. Δlog(f_T) and Δlog(f_T·m_T) ────────────────────────
    # (only meaningful when we have at least two years)
    if length(cf_years) >= 2
        yr_pairs = [(cf_years[i], cf_years[i+1]) for i in 1:length(cf_years)-1]
        pair_colors  = [CB_ORANGE, CB_GREEN, CB_PURPLE]
        pair_styles  = [:solid, :dash, :dot]

        for (ig, glabel) in enumerate(["Female", "Male"])
            # Conditional density shifts
            p = plot(;
                title  = "$(spec.shortname): Δlog f_T ($glabel)",
                xlabel = "Idiosyncratic Ability", ylabel = L"\Delta \log f_T",
                legend = :topright, size = SZ_NORMAL)
            for (ip, (yr0, yr1)) in enumerate(yr_pairs)
                fT0 = DD[yr0]["f_T"][trim_cf, DD[yr0]["iHH"], ig]
                fT1 = DD[yr1]["f_T"][trim_cf, DD[yr1]["iHH"], ig]
                Δ   = log.(max.(fT1, ε)) .- log.(max.(fT0, ε))
                plot!(p, a_grid_cf[trim_cf], Δ;
                      label="$yr1 − $yr0",
                      color=pair_colors[min(ip, length(pair_colors))],
                      linestyle=pair_styles[min(ip, length(pair_styles))],
                      linewidth=1.8)
            end
            hline!(p, [0.0]; color=:gray, linestyle=:dot, linewidth=0.6, label="")
            savefig(p, string(PLOTPATH, "delta_log_fT_$(lowercase(glabel)).png"))

            # Unconditional mass shifts
            p2 = plot(;
                title  = "$(spec.shortname): Δlog(f_T·m_T) ($glabel)",
                xlabel = "Idiosyncratic Ability",
                ylabel = L"\Delta \log(f_T \cdot m_T)",
                legend = :topright, size = SZ_NORMAL)
            for (ip, (yr0, yr1)) in enumerate(yr_pairs)
                u0 = DD[yr0]["f_T"][trim_cf, DD[yr0]["iHH"], ig] .* get_mass_T_cf(cf_id, yr0)[ig]
                u1 = DD[yr1]["f_T"][trim_cf, DD[yr1]["iHH"], ig] .* get_mass_T_cf(cf_id, yr1)[ig]
                Δ  = log.(max.(u1, ε)) .- log.(max.(u0, ε))
                plot!(p2, a_grid_cf[trim_cf], Δ;
                      label="$yr1 − $yr0",
                      color=pair_colors[min(ip, length(pair_colors))],
                      linestyle=pair_styles[min(ip, length(pair_styles))],
                      linewidth=1.8)
            end
            hline!(p2, [0.0]; color=:gray, linestyle=:dot, linewidth=0.6, label="")
            savefig(p2, string(PLOTPATH, "delta_log_fT_mass_$(lowercase(glabel)).png"))

            # Net flow
            p3 = plot(;
                title  = "$(spec.shortname): Net Flow Into Teaching ($glabel)",
                xlabel = "Idiosyncratic Ability",
                ylabel = L"f_T \cdot m_T \; \mathrm{(new)} - f_T \cdot m_T \; \mathrm{(old)}",
                legend = :topright, size = SZ_NORMAL)
            for (ip, (yr0, yr1)) in enumerate(yr_pairs)
                u0 = DD[yr0]["f_T"][trim_cf, DD[yr0]["iHH"], ig] .* get_mass_T_cf(cf_id, yr0)[ig]
                u1 = DD[yr1]["f_T"][trim_cf, DD[yr1]["iHH"], ig] .* get_mass_T_cf(cf_id, yr1)[ig]
                plot!(p3, a_grid_cf[trim_cf], u1 .- u0;
                      label="$yr1 − $yr0",
                      color=pair_colors[min(ip, length(pair_colors))],
                      fillrange=0, fillalpha=0.2, linewidth=1.4,
                      linestyle=pair_styles[min(ip, length(pair_styles))])
            end
            hline!(p3, [0.0]; color=:black, linewidth=0.5, label="")
            savefig(p3, string(PLOTPATH, "net_flow_teaching_$(lowercase(glabel)).png"))
        end
    end

    # ── A11. Occupational shares ──────────────────────────────────
    # Sort by benchmark 1970 female shares so the same occupation order
    # is used consistently across all plots.
    if haskey(DD[cf_years[1]], "share_occ")
        share_ref_1970 = (haskey(BENCH, 1970) && haskey(BENCH[1970], "share_occ")) ?
                         BENCH[1970]["share_occ"][:, 1] :
                         DD[cf_years[1]]["share_occ"][:, 1]
        sidx_occ_a  = sortperm(share_ref_1970; rev=true)
        labels_occ_a = OCC_LABELS[sidx_occ_a]

        for (ig, glabel) in enumerate(["Female", "Male"])
            p = occ_scatter(;
                title  = "$(spec.shortname): $glabel Occupational Shares",
                ylabel = "Share of $glabel Non-Teaching Workers")
            for yr in cf_years
                so = DD[yr]["share_occ"][:, ig]
                scatter!(p, 1:N_OCC, so[sidx_occ_a];
                         label=string(yr), color=YEAR_COLORS[yr],
                         markershape=YEAR_MARKERS[yr], markersize=4.5, msw=0.4)
            end
            xticks!(p, 1:N_OCC, labels_occ_a; rotation=50, tickfontsize=6)
            savefig(p, string(PLOTPATH, "occ_shares_$(lowercase(glabel)).png"))
        end
    end

    println("  ✓ Within-counterfactual plots saved to $PLOTPATH")
end


# ═══════════════════════════════════════════════════════════════════
#  PART B:  BASELINE vs COUNTERFACTUAL COMPARISON PLOTS
#           (overlay benchmark and CF for each year, per CF)
# ═══════════════════════════════════════════════════════════════════

for spec in CF_SPECS
    cf_id = spec.id
    cf_years = sort(collect(keys(CF_DATA[cf_id])))
    isempty(cf_years) && continue

    PLOTPATH = string("./plots/", PARAMNAME, "/counter_", cf_id, "/")
    mkpath(PLOTPATH)
    println("\n── Baseline vs $(spec.shortname) comparison plots ──")

    for yr in cf_years
        haskey(BENCH, yr) || continue
        d_b  = BENCH[yr]
        d_cf = CF_DATA[cf_id][yr]
        iH_b  = d_b["iHH"]
        iH_cf = d_cf["iHH"]

        tag = string(yr)

        # ── B1. f_T comparison ───────────────────────────────────
        for (ig, glabel) in enumerate(["Female", "Male"])
            p = plot(;
                title  = "$glabel Teachers' Ability Distribution ($yr)",
                xlabel = "Teaching Ability", ylabel = "Density",
                legend = :topright, size = SZ_NORMAL)
            plot!(p, A_GRID[TRIM], d_b["f_T"][TRIM, iH_b, ig];
                  label="Baseline", COMP_BASELINE...)
            plot!(p, A_GRID[TRIM], d_cf["f_T"][TRIM, iH_cf, ig];
                  label=spec.shortname, COMP_COUNTER...)
            savefig(p, string(PLOTPATH, "comp_fT_$(lowercase(glabel))_$(tag).png"))
        end

        # ── B2. h_T comparison ───────────────────────────────────
        for (ig, glabel) in enumerate(["Female", "Male"])
            p = plot(;
                title  = "$glabel Teachers' Human Capital ($yr)",
                xlabel = "Idiosyncratic Ability", ylabel = "Human Capital",
                legend = :topleft, size = SZ_NORMAL)
            plot!(p, A_GRID[TRIM], d_b["h_T"][TRIM, iH_b, ig];
                  label="Baseline", COMP_BASELINE...)
            plot!(p, A_GRID[TRIM], d_cf["h_T"][TRIM, iH_cf, ig];
                  label=spec.shortname, COMP_COUNTER...)
            savefig(p, string(PLOTPATH, "comp_hT_$(lowercase(glabel))_$(tag).png"))
        end

        # ── B3. h_T distribution comparison ──────────────────────
        for (ig, glabel) in enumerate(["Female", "Male"])
            p = plot(;
                title  = "$glabel Teachers' HC Dist. ($yr)",
                xlabel = "Human Capital", ylabel = "Density",
                legend = :topright, size = SZ_NORMAL)
            plot!(p, d_b["h_T"][TRIM, iH_b, ig], d_b["f_T"][TRIM, iH_b, ig];
                  label="Baseline", COMP_BASELINE...)
            plot!(p, d_cf["h_T"][TRIM, iH_cf, ig], d_cf["f_T"][TRIM, iH_cf, ig];
                  label=spec.shortname, COMP_COUNTER...)
            savefig(p, string(PLOTPATH, "comp_hT_fT_$(lowercase(glabel))_$(tag).png"))
        end

        # ── B4. τ_w comparison (women) ───────────────────────────
        τ_b  = d_b["τ_w"][:, 1]
        τ_cf = d_cf["τ_w"][:, 1]
        sidx = sortperm(τ_b)
        labels = OCC_LABELS[sidx]

        p = occ_scatter(; title="Barriers Against Women ($yr)",
                          ylabel=L"\tau_w")
        scatter!(p, 1:N_OCC, τ_b[sidx];
                 label="Baseline", color=CB_BLUE,
                 markershape=:circle, markersize=4.5, msw=0.4)
        scatter!(p, 1:N_OCC, τ_cf[sidx];
                 label=spec.shortname, color=CB_RED,
                 markershape=:diamond, markersize=4.5, msw=0.4)
        xticks!(p, 1:N_OCC, labels; rotation=50, tickfontsize=6)
        savefig(p, string(PLOTPATH, "comp_tau_w_women_$(tag).png"))

        # ── B5. e_T comparison ───────────────────────────────────
        eT_b  = compute_e_T(d_b)
        eT_cf = compute_e_T(d_cf)
        p = plot(;
            title  = "Teachers' HC Investment ($yr)",
            xlabel = "Idiosyncratic Ability", ylabel = L"e_T",
            legend = :topleft, size = SZ_NORMAL)
        plot!(p, A_GRID[TRIM], eT_b[TRIM, 1];
              label="Baseline", COMP_BASELINE...)
        plot!(p, A_GRID[TRIM], eT_cf[TRIM, 1];
              label=spec.shortname, COMP_COUNTER...)
        savefig(p, string(PLOTPATH, "comp_eT_$(tag).png"))

        # ── B6. Unnormalised f_T · m_T comparison ────────────────
        for (ig, glabel) in enumerate(["Female", "Male"])
            mT_b  = get_mass_T_bench(yr)[ig]
            mT_cf = get_mass_T_cf(cf_id, yr)[ig]
            p = plot(;
                title  = "Unnorm. $glabel Teacher Density ($yr)",
                xlabel = "Idiosyncratic Ability", ylabel = L"f_T \cdot m_T",
                legend = :topright, size = SZ_NORMAL)
            plot!(p, A_GRID[TRIM],
                  d_b["f_T"][TRIM, iH_b, ig] .* mT_b;
                  label="Baseline", COMP_BASELINE...)
            plot!(p, A_GRID[TRIM],
                  d_cf["f_T"][TRIM, iH_cf, ig] .* mT_cf;
                  label=spec.shortname, COMP_COUNTER...)
            savefig(p, string(PLOTPATH, "comp_fT_unnorm_$(lowercase(glabel))_$(tag).png"))
        end

        # ── B7. Occupational threshold comparison ────────────────
        if haskey(d_cf, "a_T_thresh") || haskey(d_cf, "a_O_thresh")
            thresh_key_b  = haskey(d_b,  "a_O_thresh") ? "a_O_thresh" : "a_T_thresh"
            thresh_key_cf = haskey(d_cf, "a_O_thresh") ? "a_O_thresh" : "a_T_thresh"
            trim_th = 3:N_A-2
            for (ig, glabel) in enumerate(["Women", "Men"])
                p = plot(;
                    title  = "Occ. Threshold ($glabel, $yr)",
                    xlabel = "Ability in Teaching", ylabel = "ā₁(a)",
                    legend = :topleft, size = SZ_NORMAL)
                plot!(p, A_GRID[trim_th], d_b[thresh_key_b][trim_th, iH_b, ig, 1];
                      label="Baseline", COMP_BASELINE...)
                plot!(p, A_GRID[trim_th], d_cf[thresh_key_cf][trim_th, iH_cf, ig, 1];
                      label=spec.shortname, COMP_COUNTER...)
                savefig(p, string(PLOTPATH, "comp_a_O_$(lowercase(glabel))_$(tag).png"))
            end
        end

        # ── B8. Occupational share comparison + difference ────────
        if haskey(d_b, "share_occ") && haskey(d_cf, "share_occ")
            so_b  = d_b["share_occ"]
            so_cf = d_cf["share_occ"]

            # Sort by benchmark 1970 female shares for consistency
            share_ref_b = (haskey(BENCH, 1970) && haskey(BENCH[1970], "share_occ")) ?
                          BENCH[1970]["share_occ"][:, 1] : so_b[:, 1]
            sidx_occ_b  = sortperm(share_ref_b; rev=true)
            labels_occ_b = OCC_LABELS[sidx_occ_b]

            # B8a. Comparison scatter: baseline vs CF
            for (ig, glabel) in enumerate(["Female", "Male"])
                p = occ_scatter(;
                    title  = "$glabel Occ. Shares: Baseline vs $(spec.shortname) ($yr)",
                    ylabel = "Share of $glabel Non-Teaching Workers")
                scatter!(p, 1:N_OCC, so_b[sidx_occ_b, ig];
                         label="Baseline", color=CB_BLUE,
                         markershape=:circle, markersize=4.5, msw=0.4)
                scatter!(p, 1:N_OCC, so_cf[sidx_occ_b, ig];
                         label=spec.shortname, color=CB_RED,
                         markershape=:diamond, markersize=4.5, msw=0.4)
                xticks!(p, 1:N_OCC, labels_occ_b; rotation=50, tickfontsize=6)
                savefig(p, string(PLOTPATH, "comp_occ_shares_$(lowercase(glabel))_$(tag).png"))
            end

            # B8b. Difference bar chart: CF − Baseline
            for (ig, glabel) in enumerate(["Female", "Male"])
                Δso = so_cf[:, ig] .- so_b[:, ig]
                Δso_sorted = Δso[sidx_occ_b]
                bar_colors = [Δ >= 0 ? CB_BLUE : CB_RED for Δ in Δso_sorted]

                p = occ_scatter(;
                    title  = "$glabel Occ. Share Diff: $(spec.shortname) − Baseline ($yr)",
                    ylabel = "Δ Share of $glabel Non-Teaching Workers")
                bar!(p, 1:N_OCC, Δso_sorted;
                     color=bar_colors, alpha=0.8, bar_width=0.7, label="")
                hline!(p, [0.0]; color=:black, linewidth=0.6, label="")
                xticks!(p, 1:N_OCC, labels_occ_b; rotation=50, tickfontsize=6)
                savefig(p, string(PLOTPATH, "diff_occ_shares_$(lowercase(glabel))_$(tag).png"))
            end
        end
    end

    # ── B9. Teaching shares: baseline vs CF grouped bar ──────────
    common_years = sort(intersect(cf_years, collect(keys(BENCH))))
    if length(common_years) >= 1
        n = length(common_years)
        fem_b  = [get_mass_T_bench(yr)[1] for yr in common_years]
        mal_b  = [get_mass_T_bench(yr)[2] for yr in common_years]
        fem_cf = [get_mass_T_cf(cf_id, yr)[1] for yr in common_years]
        mal_cf = [get_mass_T_cf(cf_id, yr)[2] for yr in common_years]

        # Grouped bar: 4 bars per year (baseline F/M, CF F/M)
        bw = 0.18
        xt = collect(1:n)

        p = plot(;
            title  = "Teaching Shares: Baseline vs $(spec.shortname)",
            ylabel = "Fraction", legend = :topright,
            size = SZ_BAR, xlims = (0.2, n + 0.8),
            xticks = (xt, string.(common_years)))
        bar!(p, xt .- 1.5*bw, fem_b;  bar_width=bw, color=CB_BLUE,   alpha=0.85, label="Baseline F")
        bar!(p, xt .- 0.5*bw, fem_cf; bar_width=bw, color=CB_CYAN,   alpha=0.85, label="CF F")
        bar!(p, xt .+ 0.5*bw, mal_b;  bar_width=bw, color=CB_ORANGE, alpha=0.85, label="Baseline M")
        bar!(p, xt .+ 1.5*bw, mal_cf; bar_width=bw, color=CB_YELLOW, alpha=0.85, label="CF M")
        savefig(p, string(PLOTPATH, "comp_mass_T_bar.png"))
    end

    println("  ✓ Comparison plots saved to $PLOTPATH")
end


# ═══════════════════════════════════════════════════════════════════
#  PART C:  CROSS-COUNTERFACTUAL SUMMARY PLOTS
#           (compare outcomes across CF1–CF6 for each year)
# ═══════════════════════════════════════════════════════════════════

COMPPATH = string("./plots/", PARAMNAME, "/counter_comparison/")
mkpath(COMPPATH)
println("\n── Cross-counterfactual comparison plots ──")

# ── C1. Teaching shares across counterfactuals ───────────────────
for yr in ALL_YEARS
    # Collect shares: baseline + each CF
    labels_bar = String[]
    fem_vals   = Float64[]
    mal_vals   = Float64[]

    if haskey(BENCH, yr)
        push!(labels_bar, "Baseline")
        mT = get_mass_T_bench(yr)
        push!(fem_vals, mT[1])
        push!(mal_vals, mT[2])
    end
    for spec in CF_SPECS
        if haskey(CF_DATA[spec.id], yr)
            push!(labels_bar, "CF$(spec.id)")
            mT = get_mass_T_cf(spec.id, yr)
            push!(fem_vals, mT[1])
            push!(mal_vals, mT[2])
        end
    end
    length(labels_bar) < 2 && continue

    n  = length(labels_bar)
    bw = 0.35
    xt = collect(1:n)

    p = plot(;
        title  = "Teaching Shares Across Counterfactuals ($yr)",
        ylabel = "Fraction", legend = :topright,
        size = (max(560, 60*n), 400),
        xlims = (0.3, n + 0.7),
        xticks = (xt, labels_bar),
        bottom_margin = 8Plots.mm)
    bar!(p, xt .- bw/2, fem_vals; bar_width=bw, color=CB_BLUE,   alpha=0.85, label="Female")
    bar!(p, xt .+ bw/2, mal_vals; bar_width=bw, color=CB_ORANGE, alpha=0.85, label="Male")
    savefig(p, string(COMPPATH, "mass_T_all_cf_$(yr).png"))
end

# ── C2. Aggregate productivity across counterfactuals ────────────
for yr in ALL_YEARS
    labels_bar = String[]
    aggA_vals_w_HP  = Float64[]
    aggA_vals_wo_HP = Float64[]

    # Baseline aggregate productivity
    if haskey(BENCH, yr)
        # Compute from benchmark moments if possible, otherwise skip
        # aggA is not directly stored in benchmark JLD in all configurations,
        # so we try to grab it from the 1970 aggA file or moments
        # For simplicity, check if it was stored:
        push!(labels_bar, "Baseline")
        # We compute aggA from benchmark at converged steady state using same approach
        # as the model: but it's simpler to report the CF aggA values vs each other
        push!(aggA_vals_w_HP, NaN)   # placeholder — benchmark aggA not in standard JLD
        push!(aggA_vals_wo_HP, NaN)
    end
    for spec in CF_SPECS
        if haskey(CF_DATA[spec.id], yr)
            d = CF_DATA[spec.id][yr]
            if haskey(d, "aggA")
                push!(labels_bar, "CF$(spec.id)")
                aA = d["aggA"]
                push!(aggA_vals_w_HP,  aA[1])
                push!(aggA_vals_wo_HP, aA[2])
            end
        end
    end

    # Filter out baselines with NaN if no CF data either
    valid = .!isnan.(aggA_vals_wo_HP)
    any(valid) || continue

    # Only plot CFs (drop baseline placeholder if NaN)
    idx = findall(valid)
    n  = length(idx)
    n < 1 && continue
    bw = 0.35
    xt = collect(1:n)

    p = plot(;
        title  = "Aggregate Productivity Across Counterfactuals ($yr)",
        ylabel = "Aggregate Productivity", legend = :topright,
        size = (max(560, 60*n), 400),
        xlims = (0.3, n + 0.7),
        xticks = (xt, labels_bar[idx]),
        bottom_margin = 8Plots.mm)
    bar!(p, xt .- bw/2, aggA_vals_w_HP[idx];
         bar_width=bw, color=CB_BLUE, alpha=0.85, label="With HP")
    bar!(p, xt .+ bw/2, aggA_vals_wo_HP[idx];
         bar_width=bw, color=CB_ORANGE, alpha=0.85, label="Excl. HP")
    savefig(p, string(COMPPATH, "aggA_all_cf_$(yr).png"))
end

# ── C3. Female f_T overlay: all CFs for a given year ─────────────
for yr in ALL_YEARS
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$glabel Teachers' Ability Dist. — All CFs ($yr)",
            xlabel = "Idiosyncratic Ability", ylabel = "Density",
            legend = :topright, size = SZ_NORMAL)

        # Baseline
        if haskey(BENCH, yr)
            iH = BENCH[yr]["iHH"]
            plot!(p, A_GRID[TRIM], BENCH[yr]["f_T"][TRIM, iH, ig];
                  label="Baseline", color=:black, linewidth=2.0, linestyle=:solid)
        end

        # Each counterfactual
        cf_lstyles = [:dash, :dot, :dashdot, :dashdotdot, :dash, :dot]
        for (ci, spec) in enumerate(CF_SPECS)
            haskey(CF_DATA[spec.id], yr) || continue
            d  = CF_DATA[spec.id][yr]
            iH = d["iHH"]
            plot!(p, A_GRID[TRIM], d["f_T"][TRIM, iH, ig];
                  label="CF$(spec.id)",
                  color=CF_COLORS[ci],
                  linestyle=cf_lstyles[min(ci, length(cf_lstyles))],
                  linewidth=1.4)
        end
        savefig(p, string(COMPPATH, "fT_all_cf_$(lowercase(glabel))_$(yr).png"))
    end
end

# ── C4. Female h_T overlay: all CFs for a given year ─────────────
for yr in ALL_YEARS
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "$glabel Teachers' Human Capital — All CFs ($yr)",
            xlabel = "Idiosyncratic Ability", ylabel = "Human Capital",
            legend = :topleft, size = SZ_NORMAL)

        if haskey(BENCH, yr)
            iH = BENCH[yr]["iHH"]
            plot!(p, A_GRID[TRIM], BENCH[yr]["h_T"][TRIM, iH, ig];
                  label="Baseline", color=:black, linewidth=2.0, linestyle=:solid)
        end

        cf_lstyles = [:dash, :dot, :dashdot, :dashdotdot, :dash, :dot]
        for (ci, spec) in enumerate(CF_SPECS)
            haskey(CF_DATA[spec.id], yr) || continue
            d  = CF_DATA[spec.id][yr]
            iH = d["iHH"]
            plot!(p, A_GRID[TRIM], d["h_T"][TRIM, iH, ig];
                  label="CF$(spec.id)",
                  color=CF_COLORS[ci],
                  linestyle=cf_lstyles[min(ci, length(cf_lstyles))],
                  linewidth=1.4)
        end
        savefig(p, string(COMPPATH, "hT_all_cf_$(lowercase(glabel))_$(yr).png"))
    end
end

# ── C5. τ_w overlay: all CFs for a given year ────────────────────
for yr in ALL_YEARS
    haskey(BENCH, yr) || continue
    τ_ref = BENCH[yr]["τ_w"][:, 1]
    sidx  = sortperm(τ_ref)
    labels = OCC_LABELS[sidx]

    p = occ_scatter(; title="Barriers Against Women — All CFs ($yr)",
                      ylabel=L"\tau_w")
    scatter!(p, 1:N_OCC, τ_ref[sidx];
             label="Baseline", color=:black,
             markershape=:circle, markersize=5, msw=0.5)

    cf_markers = [:diamond, :xcross, :utriangle, :dtriangle, :star5, :hexagon]
    for (ci, spec) in enumerate(CF_SPECS)
        haskey(CF_DATA[spec.id], yr) || continue
        τ_cf = CF_DATA[spec.id][yr]["τ_w"][:, 1]
        scatter!(p, 1:N_OCC, τ_cf[sidx];
                 label="CF$(spec.id)", color=CF_COLORS[ci],
                 markershape=cf_markers[min(ci, length(cf_markers))],
                 markersize=4, msw=0.3)
    end
    xticks!(p, 1:N_OCC, labels; rotation=50, tickfontsize=6)
    savefig(p, string(COMPPATH, "tau_w_all_cf_$(yr).png"))
end

# ── C6. Unnormalised f_T · m_T overlay ───────────────────────────
for yr in ALL_YEARS
    for (ig, glabel) in enumerate(["Female", "Male"])
        p = plot(;
            title  = "Unnorm. $glabel Teacher Density — All CFs ($yr)",
            xlabel = "Idiosyncratic Ability", ylabel = L"f_T \cdot m_T",
            legend = :topright, size = SZ_NORMAL)

        if haskey(BENCH, yr)
            iH  = BENCH[yr]["iHH"]
            mT  = get_mass_T_bench(yr)[ig]
            plot!(p, A_GRID[TRIM],
                  BENCH[yr]["f_T"][TRIM, iH, ig] .* mT;
                  label="Baseline", color=:black, linewidth=2.0, linestyle=:solid)
        end

        cf_lstyles = [:dash, :dot, :dashdot, :dashdotdot, :dash, :dot]
        for (ci, spec) in enumerate(CF_SPECS)
            haskey(CF_DATA[spec.id], yr) || continue
            d  = CF_DATA[spec.id][yr]
            iH = d["iHH"]
            mT = get_mass_T_cf(spec.id, yr)[ig]
            plot!(p, A_GRID[TRIM],
                  d["f_T"][TRIM, iH, ig] .* mT;
                  label="CF$(spec.id)",
                  color=CF_COLORS[ci],
                  linestyle=cf_lstyles[min(ci, length(cf_lstyles))],
                  linewidth=1.4)
        end
        savefig(p, string(COMPPATH, "fT_unnorm_all_cf_$(lowercase(glabel))_$(yr).png"))
    end
end

# ── C7. e_T overlay: all CFs for a given year ───────────────────
for yr in ALL_YEARS
    p = plot(;
        title  = "Teachers' HC Investment — All CFs ($yr)",
        xlabel = "Idiosyncratic Ability", ylabel = L"e_T",
        legend = :topleft, size = SZ_NORMAL)

    if haskey(BENCH, yr)
        eT = compute_e_T(BENCH[yr])
        plot!(p, A_GRID[TRIM], eT[TRIM, 1];
              label="Baseline", color=:black, linewidth=2.0, linestyle=:solid)
    end

    cf_lstyles = [:dash, :dot, :dashdot, :dashdotdot, :dash, :dot]
    for (ci, spec) in enumerate(CF_SPECS)
        haskey(CF_DATA[spec.id], yr) || continue
        eT = compute_e_T(CF_DATA[spec.id][yr])
        plot!(p, A_GRID[TRIM], eT[TRIM, 1];
              label="CF$(spec.id)",
              color=CF_COLORS[ci],
              linestyle=cf_lstyles[min(ci, length(cf_lstyles))],
              linewidth=1.4)
    end
    savefig(p, string(COMPPATH, "eT_all_cf_$(yr).png"))
end

# ── C8. Summary table: key moments across all CFs ───────────────
# Collect into a DataFrame and write CSV
summary_rows = []
for yr in ALL_YEARS
    # Baseline row
    if haskey(BENCH_MOM, yr)
        df_b = BENCH_MOM[yr]
        r = nrow(df_b)
        push!(summary_rows, (
            year = yr,
            scenario = "Baseline",
            share_T_fem = df_b[r, :share_teachers_female],
            share_T_mal = df_b[r, :share_teachers_male],
            p90p10_T    = df_b[r, :p90_p10_ω_teachers],
            p90p10_O    = df_b[r, :p90_p10_w_other],
            HH_fp       = haskey(BENCH[yr], "HH_fp") ? BENCH[yr]["HH_fp"] : NaN,
            aggA_w_HP   = NaN,
            aggA_wo_HP  = NaN,
        ))
    end
    # CF rows
    for spec in CF_SPECS
        haskey(CF_DATA[spec.id], yr) || continue
        d = CF_DATA[spec.id][yr]
        mT = d["mass_T"]
        aA = haskey(d, "aggA") ? d["aggA"] : [NaN, NaN]
        ω9010  = haskey(d, "ω_90_10") ? d["ω_90_10"] : fill(NaN, 1, 3)
        w9010  = haskey(d, "w_90_10") ? d["w_90_10"] : fill(NaN, 1, 3, 21)
        iH = d["iHH"]
        # ω_90_10 is stored as the full array from the CF; extract pooled value
        p90p10_T = try ω9010[iH, end] catch; NaN end
        p90p10_O = try w9010[iH, end, end] catch; NaN end

        push!(summary_rows, (
            year = yr,
            scenario = "CF$(spec.id)",
            share_T_fem = mT[1],
            share_T_mal = mT[2],
            p90p10_T    = p90p10_T,
            p90p10_O    = p90p10_O,
            HH_fp       = d["HH_fp"],
            aggA_w_HP   = aA[1],
            aggA_wo_HP  = aA[2],
        ))
    end
end

if !isempty(summary_rows)
    summary_df = DataFrame(summary_rows)
    summary_csv = string(COMPPATH, "summary_moments.csv")
    CSV.write(summary_csv, summary_df)
    println("\n✓ Summary CSV saved to $summary_csv")
    println(summary_df)
end

# ── C9. Occupational shares: all CFs for each year ───────────────
# One plot per year per gender; baseline + each available CF overlaid.
for yr in ALL_YEARS
    # Use benchmark 1970 female shares for sorting
    share_ref_c = (haskey(BENCH, 1970) && haskey(BENCH[1970], "share_occ")) ?
                  BENCH[1970]["share_occ"][:, 1] : nothing
    share_ref_c === nothing && continue

    sidx_occ_c  = sortperm(share_ref_c; rev=true)
    labels_occ_c = OCC_LABELS[sidx_occ_c]

    for (ig, glabel) in enumerate(["Female", "Male"])
        p = occ_scatter(;
            title  = "$glabel Occupational Shares — All CFs ($yr)",
            ylabel = "Share of $glabel Non-Teaching Workers")

        # Baseline
        if haskey(BENCH, yr) && haskey(BENCH[yr], "share_occ")
            so_b = BENCH[yr]["share_occ"][:, ig]
            scatter!(p, 1:N_OCC, so_b[sidx_occ_c];
                     label="Baseline", color=:black,
                     markershape=:circle, markersize=5, msw=0.5)
        end

        # Each counterfactual
        cf_markers = [:diamond, :xcross, :utriangle, :dtriangle, :star5, :hexagon]
        for (ci, spec) in enumerate(CF_SPECS)
            haskey(CF_DATA[spec.id], yr) || continue
            d_c = CF_DATA[spec.id][yr]
            haskey(d_c, "share_occ") || continue
            so_c = d_c["share_occ"][:, ig]
            scatter!(p, 1:N_OCC, so_c[sidx_occ_c];
                     label="CF$(spec.id)", color=CF_COLORS[ci],
                     markershape=cf_markers[min(ci, length(cf_markers))],
                     markersize=4, msw=0.3)
        end

        xticks!(p, 1:N_OCC, labels_occ_c; rotation=50, tickfontsize=6)
        savefig(p, string(COMPPATH, "occ_shares_all_cf_$(lowercase(glabel))_$(yr).png"))
    end
end

# ── C10. Occupational share differences vs baseline: all CFs ─────
# For each year and gender, one plot showing (CF − Baseline) for
# every counterfactual overlaid, so changes in allocation can be
# compared across counterfactuals.
for yr in ALL_YEARS
    haskey(BENCH, yr) && haskey(BENCH[yr], "share_occ") || continue

    share_ref_c = (haskey(BENCH, 1970) && haskey(BENCH[1970], "share_occ")) ?
                  BENCH[1970]["share_occ"][:, 1] : BENCH[yr]["share_occ"][:, 1]
    sidx_occ_c  = sortperm(share_ref_c; rev=true)
    labels_occ_c = OCC_LABELS[sidx_occ_c]

    for (ig, glabel) in enumerate(["Female", "Male"])
        so_base = BENCH[yr]["share_occ"][:, ig]

        p = occ_scatter(;
            title  = "$glabel Occ. Share Diff vs Baseline — All CFs ($yr)",
            ylabel = "Δ Share of $glabel Non-Teaching Workers (CF − Baseline)")
        hline!(p, [0.0]; color=:black, linewidth=0.6, label="")

        cf_lstyles = [:dash, :dot, :dashdot, :dashdotdot, :dash, :dot]
        cf_markers = [:diamond, :xcross, :utriangle, :dtriangle, :star5, :hexagon]
        for (ci, spec) in enumerate(CF_SPECS)
            haskey(CF_DATA[spec.id], yr) || continue
            d_c = CF_DATA[spec.id][yr]
            haskey(d_c, "share_occ") || continue
            Δso = (d_c["share_occ"][:, ig] .- so_base)[sidx_occ_c]
            scatter!(p, 1:N_OCC, Δso;
                     label="CF$(spec.id)", color=CF_COLORS[ci],
                     markershape=cf_markers[min(ci, length(cf_markers))],
                     markersize=4, msw=0.3)
        end

        xticks!(p, 1:N_OCC, labels_occ_c; rotation=50, tickfontsize=6)
        savefig(p, string(COMPPATH, "diff_occ_shares_all_cf_$(lowercase(glabel))_$(yr).png"))
    end
end

println("\n✓ All counterfactual plots complete.")