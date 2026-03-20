# Generate plots for presentation slides:

#= commenting but keeping the old version of this script for now
# Select type of parameterization (from following list):
    # (1) benchmark
    # (2) counter_1
    # (3) counter_2
    # (4) high_beta
    # (5) low_beta
    # (6) high_beta_high_sigma
    # (7) low_beta_low_sigma
    paramname = "benchmark"
    pathname = string("./plots/",paramname,"/")
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
# plt5 = plot(a_grid[2:end-2],e_T[2:end-2,iHH,2],label=string(year))
# (1.d) Law of motion for 'H_T'
plt6 = plot(H_grid./H_grid[iHH],H_grid./H_grid[iHH],linewidth=.5,linestyle=:dash)
plot!(plt6,H_grid./H_grid[iHH],HH_T./HH_T[iHH],linewidth=1,linestyle=:solid)
# (1.e) Occupational threshold for occupation 1 ("Executives,.."):
# (1.e.1) Women
plt7w = plot(a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year))
# (1.e.2) Men
plt7m = plot(a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year))

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
# plot!(plt5,a_grid[2:end-2],e_T[2:end-2,iHH,2],label=string(year))
# (2.e.1)
plot!(plt7w,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year))
# (2.e.1)
plot!(plt7m,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year))

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
plot!(plt7w,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,1,1],label=string(year))
# (3.e.2)
plot!(plt7m,a_grid[2:end-1],a_O_thresh[2:end-1,iHH,2,1],label=string(year))

# (4) Add axis labels and titles:
# (4.a.0)
plot!(plt1,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Abilities",ylims=(-0.001,0.11),subplot=1)
plot!(plt2,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Abilities",ylims=(-0.001,0.065),subplot=1)
# (4.a.1)
plot!(plt1_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
plot!(plt2_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
# (4.a.2)
plot!(plt1_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Human Capital",subplot=1)
plot!(plt2_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Human Capital",subplot=1)
# (4.b)
plot!(plt3,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,legend=:right,subplot=1)

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
    plots_benchmark.jl
    ───────────────────
    Script for the teacher occupational choice model.
    Loads saved JLD parameterization files and moments CSVs, then generates
    PNG plots comparing steady-state results across calibration years.

    Plots generated (all PNG):
      1.  fT_female_steadystate        — Female teachers' ability distribution
      2.  fT_male_steadystate          — Male teachers' ability distribution
      3.  hT_female_steadystate        — Female teachers' human capital vs ability
      4.  hT_male_steadystate          — Male teachers' human capital vs ability
      5.  hT_fT_female_steadystate     — Distribution of female teachers' h.c.
      6.  hT_fT_male_steadystate       — Distribution of male teachers' h.c.
      7.  tau_w_women                  — Labor market barriers (all years)
      8.  tau_w_women_70_90            — Barriers: 1970 vs 1990
      9.  tau_w_women_70_10            — Barriers: 1970 vs 2010
     10.  A_occ                        — Occupational productivities (sorted)
     11.  eT_steadystate               — Teachers' human capital investment
     12.  LoM                          — Law of motion for agg. teaching h.c.
     13.  a_O_women / a_O_men          — Occupational thresholds
     14.  delta_log_fT_female/male     — Δlog(f_T): conditional distribution shifts
     15.  delta_log_fT_mass_fem/male   — Δlog(f_T·m_T): unconditional mass shifts
     16.  net_flow_teaching_fem/male   — Net flow into/out of teaching by ability
     17.  mass_T_over_time             — Teaching shares over time (bar chart)
     18.  cdf_fT_female/male           — CDF of teachers' abilities
     19.  fT_unnorm_female/male        — Unnormalised teacher density (f_T · m_T)
=#

using Plots, JLD, CSV, DataFrames, LaTeXStrings
import XLSX

# ═══════════════════════════════════════════════════════════════════
#  1. CONFIGURATION
# ═══════════════════════════════════════════════════════════════════

PARAMNAME = "benchmark"
YEARS     = [1970, 1990, 2010]
BASEPATH  = string("./parameterization/", PARAMNAME, "/")
PLOTPATH  = string("./plots/", PARAMNAME, "/")
mkpath(PLOTPATH)

# Colorblind-safe palette (Okabe–Ito)
CB_BLUE   = RGB(0/255, 114/255, 178/255)
CB_ORANGE = RGB(230/255, 159/255, 0/255)
CB_GREEN  = RGB(0/255, 158/255, 115/255)
CB_COLORS = [CB_BLUE, CB_ORANGE, CB_GREEN]
LSTYLES   = [:solid, :dash, :dot]
MARKERS   = [:circle, :diamond, :xcross]

# Tiny constant to avoid log(0)
ε = 1e-12

# Use GR backend (ships with Plots.jl, no external deps)
gr()

# Global defaults — no grid
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

# Standard figure sizes
SZ_NORMAL = (620, 400)     # line / area plots
SZ_WIDE   = (960, 480)     # scatter with occupation labels
SZ_BAR    = (480, 380)     # bar chart

# ═══════════════════════════════════════════════════════════════════
#  2. LOAD DATA
# ═══════════════════════════════════════════════════════════════════

# --- JLD parameterization files ---
DATA = Dict{Int, Dict{String, Any}}()
for yr in YEARS
    local fname = string(BASEPATH, "previousParameterization", yr, ".jld")
    DATA[yr] = load(fname)
    println("Loaded JLD for year $yr")
end

# --- Moments CSVs ---
MOMENTS = Dict{Int, DataFrame}()
for yr in YEARS
    local csvname = string("./results/", PARAMNAME, "/moments", yr, ".csv")
    if isfile(csvname)
        MOMENTS[yr] = DataFrame(CSV.File(csvname))
    end
end

# --- Occupation labels from Excel ---
XLSX_PATH = "../../data/LaborMarketData/wages_occ_shares_v2.xlsx"
OCC_LABELS_FULL = if isfile(XLSX_PATH)
    local xs  = XLSX.readxlsx(XLSX_PATH)
    local tab = xs["moments_shares"]
    string.(vec(tab["A30:A49"]))
else
    ["Occ $i" for i in 1:20]
end

# Abbreviated labels: truncate to 18 chars for readability on x-axis
OCC_LABELS = [length(s) > 18 ? s[1:17] * "…" : s for s in OCC_LABELS_FULL]
N_OCC = length(OCC_LABELS)

# ═══════════════════════════════════════════════════════════════════
#  3. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════

"Return named-tuple of style kwargs for year index i."
ystyle(i) = (color=CB_COLORS[i], linestyle=LSTYLES[i], linewidth=1.8)

"Extract [female, male] teaching shares from moments CSV (last row)."
function get_mass_T(yr::Int)
    haskey(MOMENTS, yr) || return [NaN, NaN]
    local df = MOMENTS[yr]
    local r  = nrow(df)
    return [df[r, :share_teachers_female], df[r, :share_teachers_male]]
end

"""
Reconstruct e_T (human-capital investment in teaching) from saved parameters.
    e_T(a) = η · (1-t) · κ · γ · h_T(a)^γ
"""
function compute_e_T(d::Dict)
    local η = d["η"]; local γ_val = d["γ"]; local κ_val = d["κ"]
    local t_val = d["t"]; local hT = d["h_T"]; local iH = d["iHH"]
    local na, _, nG = size(hT)
    local eT = zeros(na, nG)
    for iG in 1:nG, ia in 1:na
        eT[ia, iG] = η * (1 - t_val[iH, 1]) * κ_val * γ_val * hT[ia, iH, iG]^γ_val
    end
    return eT
end

"Create a scatter plot pre-configured for occupation-label x-axes."
function occ_scatter(; title="", ylabel="", ylims_val=nothing)
    local kw = Dict{Symbol,Any}(
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

# Common grid and index (shared across all years)
A_GRID = DATA[YEARS[1]]["a_grid"]
N_A    = length(A_GRID)
IHH    = DATA[YEARS[1]]["iHH"]
TRIM   = 3:N_A-3          # skip 3 boundary points (safer than 2)

# ═══════════════════════════════════════════════════════════════════
#  4. DISTRIBUTION OF TEACHERS' ABILITIES  (f_T)
# ═══════════════════════════════════════════════════════════════════

for (ig, glabel) in enumerate(["Female", "Male"])
    local p = plot(;
        title  = "Distribution of $glabel Teachers' Abilities",
        xlabel = "Idiosyncratic Ability", ylabel = "Density",
        legend = :topright, size = SZ_NORMAL)
    for (i, yr) in enumerate(YEARS)
        local fT = DATA[yr]["f_T"]
        plot!(p, A_GRID[TRIM], fT[TRIM, IHH, ig]; label=string(yr), ystyle(i)...)
    end
    savefig(p, string(PLOTPATH, "fT_$(lowercase(glabel))_steadystate.png"))
end

# ═══════════════════════════════════════════════════════════════════
#  5. TEACHERS' HUMAN CAPITAL vs ABILITY  (h_T)
# ═══════════════════════════════════════════════════════════════════

for (ig, glabel) in enumerate(["Female", "Male"])
    local p = plot(;
        title  = "$glabel Teachers' Human Capital",
        xlabel = "Idiosyncratic Ability", ylabel = "Human Capital",
        legend = :topleft, size = SZ_NORMAL)
    for (i, yr) in enumerate(YEARS)
        local hT = DATA[yr]["h_T"]
        plot!(p, A_GRID[TRIM], hT[TRIM, IHH, ig]; label=string(yr), ystyle(i)...)
    end
    savefig(p, string(PLOTPATH, "hT_$(lowercase(glabel))_steadystate.png"))
end

# ═══════════════════════════════════════════════════════════════════
#  6. DISTRIBUTION OF TEACHERS' HUMAN CAPITAL  (h_T on x-axis)
# ═══════════════════════════════════════════════════════════════════

for (ig, glabel) in enumerate(["Female", "Male"])
    local p = plot(;
        title  = "Distribution of $glabel Teachers' Human Capital",
        xlabel = "Human Capital", ylabel = "Density",
        legend = :topright, size = SZ_NORMAL)
    for (i, yr) in enumerate(YEARS)
        local hT = DATA[yr]["h_T"]
        local fT = DATA[yr]["f_T"]
        plot!(p, hT[TRIM, IHH, ig], fT[TRIM, IHH, ig];
              label=string(yr), ystyle(i)...)
    end
    savefig(p, string(PLOTPATH, "hT_fT_$(lowercase(glabel))_steadystate.png"))
end

# ═══════════════════════════════════════════════════════════════════
#  7. LABOR MARKET BARRIERS  (τ_w)  — sorted by 1970 magnitude
#     Y-axis uses LaTeX label
# ═══════════════════════════════════════════════════════════════════

let τ_ref  = DATA[1970]["τ_w"][:, 1],
    sidx   = sortperm(τ_ref),
    labels = OCC_LABELS[sidx]

    # Define the progression of years to plot and their file suffixes
    local plot_setups = [
        (1970, "70"),
        (1990, "70_90"),
        (2010, "70_10")
    ]
    for (yr_end, suffix) in plot_setups
        local p = occ_scatter(; title="Labor Market Barriers Against Women",
                                ylabel=L"\tau_w")
        
        # Filter and plot only the years up to yr_end
        for yr in filter(y -> y <= yr_end, YEARS)
            local τ_w = DATA[yr]["τ_w"][:, 1]
            local yi  = findfirst(==(yr), YEARS)
            
            scatter!(p, 1:N_OCC, τ_w[sidx];
                     label=string(yr), color=CB_COLORS[yi],
                     markershape=MARKERS[yi], markersize=4.5, msw=0.4)
        end
        
        xticks!(p, 1:N_OCC, labels; rotation=50, tickfontsize=6)
        
        # Save the plot with the progressive suffix
        savefig(p, string(PLOTPATH, "tau_w_women_$suffix.png"))
    end
end

# ═══════════════════════════════════════════════════════════════════
#  8. OCCUPATIONAL PRODUCTIVITIES  — sorted by 1970 values
#     Y-axis uses LaTeX label
# ═══════════════════════════════════════════════════════════════════

let a_ref   = haskey(DATA[1970], "a_by_occ_level") ? DATA[1970]["a_by_occ_level"] : DATA[1970]["a_by_occ"],
    sidx    = sortperm(a_ref),
    labels  = OCC_LABELS[sidx]

    local p = occ_scatter(; title="Occupational Productivities", ylabel=L"A_o")
    for (i, yr) in enumerate(YEARS)
        local a = haskey(DATA[yr], "a_by_occ_level") ? DATA[yr]["a_by_occ_level"] : DATA[yr]["a_by_occ"]
        scatter!(p, 1:N_OCC, a[sidx];
                 label=string(yr), color=CB_COLORS[i],
                 markershape=MARKERS[i], markersize=4.5, msw=0.4)
    end
    xticks!(p, 1:N_OCC, labels; rotation=50, tickfontsize=6)
    savefig(p, string(PLOTPATH, "A_occ.png"))
end

# ═══════════════════════════════════════════════════════════════════
#  9. TEACHERS' HUMAN CAPITAL INVESTMENT  (e_T)
#     Y-axis uses LaTeX label
# ═══════════════════════════════════════════════════════════════════

let p = plot(;
        title  = "Teachers' Human Capital Investment",
        xlabel = "Idiosyncratic Ability", ylabel = L"e_T",
        legend = :topleft, size = SZ_NORMAL)
    for (i, yr) in enumerate(YEARS)
        local eT = compute_e_T(DATA[yr])
        plot!(p, A_GRID[TRIM], eT[TRIM, 1]; label=string(yr), ystyle(i)...)
    end
    savefig(p, string(PLOTPATH, "eT_steadystate.png"))
end

# ═══════════════════════════════════════════════════════════════════
# 10. LAW OF MOTION FOR AGGREGATE TEACHING HUMAN CAPITAL
#     Axes use LaTeX labels
# ═══════════════════════════════════════════════════════════════════

let d       = DATA[YEARS[end]],
    H_grid  = d["H_grid"],
    HH_T    = d["HH_T"],
    iH      = d["iHH"],
    x_norm  = H_grid ./ H_grid[iH],
    y_norm  = HH_T   ./ HH_T[iH]

    local p = plot(;
        title  = "Law of Motion for Teachers' Human Capital",
        xlabel = L"H_T \; \mathrm{(normalised)}",
        ylabel = L"H_T' \; \mathrm{(normalised)}",
        legend = false, size = SZ_NORMAL)
    plot!(p, x_norm, x_norm; color=:gray, linewidth=0.8, linestyle=:dash)
    plot!(p, x_norm, y_norm; color=CB_BLUE, linewidth=2)
    savefig(p, string(PLOTPATH, "LoM.png"))
end

# ═══════════════════════════════════════════════════════════════════
# 11. OCCUPATIONAL THRESHOLDS
# ═══════════════════════════════════════════════════════════════════

if haskey(DATA[YEARS[1]], "a_O_thresh")
    local trim_th = 3:N_A-2
    for (ig, glabel) in enumerate(["Women", "Men"])
        local p = plot(;
            title  = latexstring("Occupational Threshold for $glabel: \$ a_{HP}^*(a_T) \$"),
            xlabel = latexstring("Ability in Teaching: \$ a_T \$"),
            ylabel = latexstring("\$ a_{HP} \$"),
            legend = :topleft, size = SZ_NORMAL)
        for (i, yr) in enumerate(YEARS)
            local a_O_thresh = DATA[yr]["a_O_thresh"]
            plot!(p, A_GRID[trim_th], a_O_thresh[trim_th, IHH, ig, end];
                  label=string(yr), ystyle(i)...)
        end
        savefig(p, string(PLOTPATH, "a_O_$(lowercase(glabel)).png"))
    end
end

# ═══════════════════════════════════════════════════════════════════
# 12. Δlog(f_T): SHIFTS IN CONDITIONAL ABILITY DISTRIBUTION
#     Y-axis uses LaTeX label; legend placed at topright
# ═══════════════════════════════════════════════════════════════════

for (ig, glabel) in enumerate(["Female", "Male"])
    local fT_70 = DATA[1970]["f_T"][TRIM, IHH, ig]
    local fT_90 = DATA[1990]["f_T"][TRIM, IHH, ig]
    local fT_10 = DATA[2010]["f_T"][TRIM, IHH, ig]

    local Δ1 = log.(max.(fT_90, ε)) .- log.(max.(fT_70, ε))
    local Δ2 = log.(max.(fT_10, ε)) .- log.(max.(fT_90, ε))

    local p = plot(;
        title  = "Δ log Conditional Density ($glabel Teachers)",
        xlabel = "Idiosyncratic Ability", ylabel = L"\Delta \log f_T",
        legend = :topright, size = SZ_NORMAL)
    plot!(p, A_GRID[TRIM], Δ1;
          label="1990 − 1970", color=CB_ORANGE, linewidth=1.8)
    plot!(p, A_GRID[TRIM], Δ2;
          label="2010 − 1990", color=CB_GREEN, linestyle=:dash, linewidth=1.8)
    hline!(p, [0.0]; color=:gray, linestyle=:dot, linewidth=0.6, label="")

    savefig(p, string(PLOTPATH, "delta_log_fT_$(lowercase(glabel)).png"))
end

# ═══════════════════════════════════════════════════════════════════
# 13. Δlog(f_T × mass_T): UNCONDITIONAL MASS SHIFTS
#     Y-axis uses LaTeX label; legend placed at topright
# ═══════════════════════════════════════════════════════════════════

let mass = Dict(yr => get_mass_T(yr) for yr in YEARS)
    for (ig, glabel) in enumerate(["Female", "Male"])
        local u70 = DATA[1970]["f_T"][TRIM, IHH, ig] .* mass[1970][ig]
        local u90 = DATA[1990]["f_T"][TRIM, IHH, ig] .* mass[1990][ig]
        local u10 = DATA[2010]["f_T"][TRIM, IHH, ig] .* mass[2010][ig]

        local Δ1 = log.(max.(u90, ε)) .- log.(max.(u70, ε))
        local Δ2 = log.(max.(u10, ε)) .- log.(max.(u90, ε))

        local p = plot(;
            title  = "Δ log Unconditional Teacher Mass ($glabel)",
            xlabel = "Idiosyncratic Ability",
            ylabel = L"\Delta \log(f_T \cdot m_T)",
            legend = :topright, size = SZ_NORMAL)
        plot!(p, A_GRID[TRIM], Δ1;
              label="1990 − 1970", color=CB_ORANGE, linewidth=1.8)
        plot!(p, A_GRID[TRIM], Δ2;
              label="2010 − 1990", color=CB_GREEN, linestyle=:dash, linewidth=1.8)
        hline!(p, [0.0]; color=:gray, linestyle=:dot, linewidth=0.6, label="")

        savefig(p, string(PLOTPATH, "delta_log_fT_mass_$(lowercase(glabel)).png"))
    end
end

# ═══════════════════════════════════════════════════════════════════
# 14. NET FLOW INTO / OUT OF TEACHING BY ABILITY
#     Legend placed at topright
# ═══════════════════════════════════════════════════════════════════

let mass = Dict(yr => get_mass_T(yr) for yr in YEARS)
    for (ig, glabel) in enumerate(["Female", "Male"])
        local u70 = DATA[1970]["f_T"][TRIM, IHH, ig] .* mass[1970][ig]
        local u90 = DATA[1990]["f_T"][TRIM, IHH, ig] .* mass[1990][ig]
        local u10 = DATA[2010]["f_T"][TRIM, IHH, ig] .* mass[2010][ig]

        local p = plot(;
            title  = "Net Flow Into Teaching ($glabel)",
            xlabel = "Teaching Ability",
            ylabel = L"f_T \cdot m_T \; \mathrm{(new)} - f_T \cdot m_T \; \mathrm{(old)}",
            legend = :topright, size = SZ_NORMAL)
        plot!(p, A_GRID[TRIM], u90 .- u70;
              label="1990 − 1970", color=CB_ORANGE,
              fillrange=0, fillalpha=0.2, linewidth=1.4)
        plot!(p, A_GRID[TRIM], u10 .- u90;
              label="2010 − 1990", color=CB_GREEN,
              fillrange=0, fillalpha=0.2, linewidth=1.4, linestyle=:dash)
        hline!(p, [0.0]; color=:black, linewidth=0.5, label="")

        savefig(p, string(PLOTPATH, "net_flow_teaching_$(lowercase(glabel)).png"))
    end
end

# ═══════════════════════════════════════════════════════════════════
# 15. TEACHING SHARES OVER TIME  (manual grouped bar — no StatsPlots)
#     Legend placed at topleft
# ═══════════════════════════════════════════════════════════════════

let mass     = Dict(yr => get_mass_T(yr) for yr in YEARS),
    fem      = [mass[yr][1] for yr in YEARS],
    mal      = [mass[yr][2] for yr in YEARS],
    n        = length(YEARS),
    bw       = 0.35,                      # half-width of each bar
    xticks_x = collect(1:n)

    local p = plot(;
        title  = "Share of Teachers Over Time",
        ylabel = "Fraction",
        legend = :topleft,
        size   = SZ_BAR,
        xlims  = (0.3, n + 0.7),
        xticks = (xticks_x, string.(YEARS)))

    bar!(p, xticks_x .- bw/2, fem;
         bar_width=bw, color=CB_BLUE,  alpha=0.85, label="Female")
    bar!(p, xticks_x .+ bw/2, mal;
         bar_width=bw, color=CB_ORANGE, alpha=0.85, label="Male")

    savefig(p, string(PLOTPATH, "mass_T_over_time.png"))
end

# ═══════════════════════════════════════════════════════════════════
# 16. CDF OF TEACHERS' ABILITIES
# ═══════════════════════════════════════════════════════════════════

for (ig, glabel) in enumerate(["Female", "Male"])
    local p = plot(;
        title  = "CDF of $glabel Teachers' Abilities",
        xlabel = "Idiosyncratic Ability", ylabel = "Cumulative Probability",
        legend = :bottomright, size = SZ_NORMAL)
    for (i, yr) in enumerate(YEARS)
        local fT  = DATA[yr]["f_T"][TRIM, IHH, ig]
        local da  = diff(A_GRID[TRIM])
        local cdf_vals = zeros(length(TRIM))
        for k in 2:length(TRIM)
            cdf_vals[k] = cdf_vals[k-1] + 0.5 * (fT[k] + fT[k-1]) * da[k-1]
        end
        cdf_vals ./= max(cdf_vals[end], ε)
        plot!(p, A_GRID[TRIM], cdf_vals; label=string(yr), ystyle(i)...)
    end
    savefig(p, string(PLOTPATH, "cdf_fT_$(lowercase(glabel)).png"))
end

# ═══════════════════════════════════════════════════════════════════
# 17. UNNORMALISED DISTRIBUTION OF TEACHERS' ABILITIES  (f_T · m_T)
#     Shows the mass of teachers at each ability level, reflecting
#     changes in both the conditional density and the overall share
#     of teachers over time.
# ═══════════════════════════════════════════════════════════════════

let mass = Dict(yr => get_mass_T(yr) for yr in YEARS)
    for (ig, glabel) in enumerate(["Female", "Male"])
        local p = plot(;
            title  = "Non-normalised Distribution of $glabel Teachers",
            xlabel = "Teaching Ability",
            ylabel = L"f_T \cdot m_T",
            legend = :topright, size = SZ_NORMAL)
        for (i, yr) in enumerate(YEARS)
            local fT = DATA[yr]["f_T"][TRIM, IHH, ig]
            plot!(p, A_GRID[TRIM], fT .* mass[yr][ig];
                  label=string(yr), ystyle(i)...)
        end
        savefig(p, string(PLOTPATH, "fT_unnorm_$(lowercase(glabel)).png"))
    end
end

# ═══════════════════════════════════════════════════════════════════
# 18. OCCUPATIONAL SHARES  — sorted by 1970 female share (descending)
#     share_occ is (n_O-1 × n_G): rows = non-teaching occupations,
#     columns = [female, male]; each column sums to 1 within each gender.
# ═══════════════════════════════════════════════════════════════════

if haskey(DATA[1970], "share_occ")
    local share_ref  = DATA[1970]["share_occ"][:, 1]   # 1970 female shares for ordering
    local sidx_occ   = sortperm(share_ref; rev=true)    # descending: largest share first
    local labels_occ = OCC_LABELS[sidx_occ]

    for (ig, glabel) in enumerate(["Female", "Male"])
        local p = occ_scatter(;
            title  = "$glabel Occupational Shares (Non-Teaching)",
            ylabel = "Share of $glabel Non-Teaching Workers")
        for (i, yr) in enumerate(YEARS)
            local so = DATA[yr]["share_occ"][:, ig]
            scatter!(p, 1:N_OCC, so[sidx_occ];
                     label=string(yr), color=CB_COLORS[i],
                     markershape=MARKERS[i], markersize=4.5, msw=0.4)
        end
        xticks!(p, 1:N_OCC, labels_occ; rotation=50, tickfontsize=6)
        savefig(p, string(PLOTPATH, "occ_shares_$(lowercase(glabel)).png"))
    end
end

println("\n✓ All plots saved to: $PLOTPATH")