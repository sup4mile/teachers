# Generate plots for presentation slides:

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
savefig(plt3_1_70,"./plots/counterfactuals/counter_1/tau_w_women_70.eps")
# (1.b.2) Scatter plot of occupational productivities:
plt3_2_70 = scatter(a_by_occ,label=string(year),markershape=:circle,legend=:topleft,ylims=(0,1.25))
plot!(plt3_2_70,ylabel=L"A_o",title="Occupational Productivities",grid=false,subplot=1)
savefig(plt3_2_70,"./plots/counterfactuals/counter_1/A_men_70.eps")
# (1.c) Distribution of 'e_T'
plt4 = plot(a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))

# (1.d) Law of motion for 'H_T'
plt6 = plot(H_grid./H_grid[iHH],H_grid./H_grid[iHH],linewidth=.5,linestyle=:dash)
plot!(plt6,H_grid./H_grid[iHH],HH_T./HH_T[iHH],linewidth=1,linestyle=:solid)

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
savefig(plt3_1_70_90,"./plots/counterfactuals/counter_1/tau_w_women_70_90.eps")
# (2.b.2) Create plot for 1970-90 only:
plt3_2_70_90 = scatter!(plt3_2_70,a_by_occ,label=string(year),markershape=:diamond,ylims=(0,1.25))
savefig(plt3_2_70_90,"./plots/counterfactuals/counter_1/A_men_70_90.eps")
# (2.c)
plot!(plt4,a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))

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
savefig(plt3_1_70_10,"./plots/counterfactuals/counter_1/tau_w_women_70_10.eps")
# (2.b.2) Create plot for 1970-2010 only:
plt3_2_70_10 = scatter!(plt3_2_70_90,a_by_occ,label=string(year),markershape=:xcross,ylims=(0,1.25))
savefig(plt3_2_70_10,"./plots/counterfactuals/counter_1/A_men_70_10.eps")
# (3.c)
plot!(plt4,a_grid[2:end-2],e_T[2:end-2,iHH,1],label=string(year))
# plot!(plt5,a_grid[2:end-2],e_T[2:end-2,iHH,2],label=string(year))

# (4) Add axis labels and titles:
# (4.a.0)
plot!(plt1,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Abilities",subplot=1)
plot!(plt2,xlabel="Idiosyncratic Ability",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Abilities",subplot=1)
# (4.a.1)
plot!(plt1_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
plot!(plt2_1,xlabel="Idiosyncratic Ability",ylabel="Human Capital",grid=false,legend=:right,title="Teachers' Human Capital",subplot=1)
# (4.a.2)
plot!(plt1_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Female Teachers' Human Capital",subplot=1)
plot!(plt2_2,xlabel="Human Capital",ylabel="Density",grid=false,legend=:right,title="Distribution of Male Teachers' Human Capital",subplot=1)
# (4.b)
plot!(plt3,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,legend=:left,subplot=1)

plot!(plt3_1_70_10,ylabel=L"τ_w",title="Labor Market Barriers Against Women",grid=false,subplot=1)
# (4.c)
plot!(plt4,xlabel="Idiosyncratic Ability",ylabel=L"e_T",grid=false,title="Teachers' Human Capital Investment",legend=:left,subplot=1)
# plot!(plt5,xlabel="Idiosyncratic Ability",ylabel=L"e_T",grid=false,title="Male Teachers' Human Capital Investment",legend=:left,subplot=1)
plot!(plt6,xlabel=L"\widetilde{H}_T",ylabel=L"\widetilde{H}_T^'",grid=false,title="Law of Motion for Teachers' Human Capital",legend=false)

# (5) Save figures in EPS file format:
savefig(plt1,"./plots/counterfactuals/counter_1/fT_women_steadystate.eps")
savefig(plt2,"./plots/counterfactuals/counter_1/fT_men_steadystate.eps")
savefig(plt1_1,"./plots/counterfactuals/counter_1/hT_women_steadystate.eps")
savefig(plt2_1,"./plots/counterfactuals/counter_1/hT_men_steadystate.eps")
savefig(plt1_2,"./plots/counterfactuals/counter_1/hT_fT_women_steadystate.eps")
savefig(plt2_2,"./plots/counterfactuals/counter_1/hT_fT_men_steadystate.eps")
savefig(plt3,"./plots/counterfactuals/counter_1/tau_w_women.eps")
savefig(plt3_1_70_90,"./plots/counterfactuals/counter_1/tau_w_women_70_90.eps")
savefig(plt3_1_70_10,"./plots/counterfactuals/counter_1/tau_w_women_70_10.eps")
savefig(plt4,"./plots/counterfactuals/counter_1/eT_steadystate.eps")
savefig(plt6,"./plots/counterfactuals/counter_1/LoM.eps")