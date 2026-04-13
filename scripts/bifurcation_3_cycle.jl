using Revise, Plots
using BifurcationKit
const BK = BifurcationKit
using DifferentialEquations
using LinearAlgebra
using Cairo


# gr()
default(titlefont = (7, "helvetica"), legendfontsize = 6, guidefont = (6, "helvetica"), xtickfont = (6, "helvetica"), ytickfont = (6, "helvetica"), grid = false, fontfamily = "helvetica")

proj_ei = [1 1 1 0; 0 0 0 1]
proj_ctln = [1 -1 0 0; 1 1 -2 0]

function normalize_rows(mat)
	return mat ./ [norm(row) for row in eachrow(mat)]
end

proj_ei = normalize_rows(proj_ei)
proj_ctln = normalize_rows(proj_ctln)

squareplus(x, a = 0.001) = (x + sqrt(x^2 + a))/2

function cycle_3!(du, X, p, t = 0)
	(;wEE_self, wEE, wEI, wIE, wII, tau_E, tau_I, b, a) = p
	x1, x2, x3, x4 = X
	du[1] =(1/tau_E)* (-x1 + squareplus(wEE_self*x1 + wEE * x3 + wEI * x4 + b, a))
	du[2] = (1/tau_E)*(-x2 + squareplus(wEE_self*x2 + wEE * x1 + wEI * x4 + b, a))
	du[3] = (1/tau_E)* (-x3 + squareplus(wEE_self*x3 + wEE * x2 + wEI * x4 + b, a))
	du[4] = (1/tau_I)*(-x4 + squareplus(wIE * (x1 + x2 + x3) + wII * x4,  a))
	du
end

function ctln!(du, X, p, t)
	(;wEE_self, wEE, wEI, wIE, wII, tau_E, tau_I, b, a) = p
	x1, x2, x3 = X
	du[1] =(1/tau_E)* (-x1 + squareplus((wEE-wEE_self)* x3 + (-wEE_self)* x2 + b, a))
	du[2] = (1/tau_E)*(-x2 + squareplus((wEE-wEE_self)* x1 +  (-wEE_self) * x3 + b, a))
	du[3] = (1/tau_E)* (-x3 + squareplus((wEE-wEE_self)* x2 + (-wEE_self) * x1 + b, a))
	du
end

function eigs_fp(p)
	(;wEE_self, wEE, wEI, wIE, wII, tau_E, tau_I, b, a) = p
	W_mat = [wEE_self 0 wEE wEI; 
	 	wEE wEE_self 0 wEI; 
	 	0 wEE wEE_self wEI; 
	 	wIE wIE wIE wII]
	
	tau = diagm([1/tau_E; 1/tau_E; 1/tau_E; 1/tau_I])

	F = tau * (W_mat - I)
	eigs = eigen(F)
	return eigs
end

function det_fp(p)
	(;wEE_self, wEE, wEI, wIE, wII, tau_E, tau_I, b, a) = p
	W_mat = [wEE_self 0 wEE wEI; 
	 	wEE wEE_self 0 wEI; 
	 	0 wEE wEE_self wEI; 
	 	wIE wIE wIE wII]
	
	tau = diagm([1/tau_E; 1/tau_E; 1/tau_E; 1/tau_I])

	F = tau * (W_mat - I)

	return det(F)
end

wEE_self =1.5; wEE = .75; wEI  = -2.25; wIE = 2.0; wII = -2.0; b = 1.0;

W = [wEE_self 0 wEE wEI; 
	 wEE wEE_self 0 wEI; 
	 0 wEE wEE_self wEI; 
	 wIE wIE wIE wII]
	
B = b*[1; 1; 1; 0]

FP = (I-W)^-1*B

inits = [[.05; .00; 0.05; 2/30], [.05; .05; .05; .07]]


par_cycle = (wEE_self = wEE_self, wEE =  wEE , wEI  = wEI, wIE = wIE, wII =  wII, tau_E = 1.0, tau_I = 0.01, b = 1.0, a = 0.0)

z0 = 10*inits[1]

prob = BifurcationProblem(cycle_3!, z0, par_cycle, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))

prob_de = ODEProblem(cycle_3!, z0, (0,11.), par_cycle)
sol = DifferentialEquations.solve(prob_de, Rodas5())

plot(sol)



argspo = (record_from_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[4,:]),
				min = minimum(xtt[4,:]),
				period = getperiod(p.prob, x, p.p))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		plot!(xtt.t, xtt[1,:]; label = "x", k...)
		plot!(xtt.t, xtt[2,:]; label = "y", k...)
        plot!(xtt.t, xtt[3,:]; label = "x", k...)
		plot!(xtt.t, xtt[4,:]; label = "y", k...)
	end)


# function to build probtrap from sol
probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 11.)

opts_br = ContinuationPar(p_min = 0.0,p_max = 8.0, ds = 0.002, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 2000, tol_stability = 1e-2)
brpo_cycle = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )





par_cycle = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = 2.5, b = 1.0, a = 0.0)

z0 = 10*inits[2]

prob = BifurcationProblem(cycle_3!, z0, par_cycle, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))


prob_de = ODEProblem(cycle_3!, z0, (0,11.), par_cycle)
sol = DifferentialEquations.solve(prob_de, Rodas5())
#plot(sol)
#hline!(p, [sol.u[end][4]], color = colorant"black", linestyle = :dash)

probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 5.)

opts_br = ContinuationPar(p_min = 0.1,p_max = 8.0, ds = -0.002, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 2000, tol_stability = 1e-2)
brpo_ei = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )




orbit = get_periodic_orbit(brpo_ei, 10)
u = orbit[:,:]
plot(u')



# par_bif = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = 6.5, b = 1.0, a = 0.0)

# prob_de = ODEProblem(cycle_3!, x0, (0,400.), par_bif)
# sol = DifferentialEquations.solve(prob_de, Rodas5())
# plot(sol)


par_btwn_loop = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = 2.5, b = 1.0, a = 0.0001)

lambda = .225
z0 = 10*(lambda * inits[1] + (1-lambda)* inits[2])
prob_de = ODEProblem(cycle_3!, z0, (0,30.), par_btwn_loop)
sol = DifferentialEquations.solve(prob_de, Rodas5())
plot(sol)

prob = BifurcationProblem(cycle_3!, z0, par_btwn_loop, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))

probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 12.)

opts_br = ContinuationPar(p_min = 2.1,p_max = 8.0, ds = 0.0005, dsmin = 0.00005, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 750, tol_stability = 1e-2)
brpo_between_soft = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )


iter = 450
bottom_loop =  get_periodic_orbit(brpo_between_soft, iter)
plot(bottom_loop)

tau_I_bottom = brpo_between_soft.param[iter]

par_btwn_bottom = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, 
tau_I = tau_I_bottom, b = 1.0, a = 0.0)

z0_bottom = bottom_loop[:,1]

prob = BifurcationProblem(cycle_3!, z0, par_btwn_bottom, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))
prob_de = ODEProblem(cycle_3!, z0_bottom, (0,30.), par_btwn_bottom)
sol = DifferentialEquations.solve(prob_de, Rodas5())
plot(sol)

probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 15.)

opts_br = ContinuationPar(p_min = 2.1,p_max = 8.0, ds = -0.0005, dsmin = 0.00005, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 1000, tol_stability = 1e-2)
brpo_bottom_back = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )

probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 15.)

opts_br = ContinuationPar(p_min = 2.1,p_max = 8.0, ds = 0.0005, dsmin = 0.00005, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 1000, tol_stability = 1e-2)

brpo_bottom_forward = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )

bottom_loop =  get_periodic_orbit(brpo_bottom_forward, 1)
plot(bottom_loop)


iter = 380
tau_I_top = brpo_between_soft.param[iter]

top_loop =  get_periodic_orbit(brpo_between_soft, iter)
plot(top_loop)
par_btwn_top = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = tau_I_top, b = 1.0, a = 0.0)

z0_top = top_loop[:,1]

prob = BifurcationProblem(cycle_3!, z0, par_btwn_top, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))
prob_de = ODEProblem(cycle_3!, z0_top, (0,30.), par_btwn_top)
sol = DifferentialEquations.solve(prob_de, Rodas5())
plot(sol)

probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 100),
	prob, sol, 15.51)

opts_br = ContinuationPar(p_min = 2.1,p_max = 8.0, ds = 0.0005, dsmin = 0.00005, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 1000, tol_stability = 1e-2)
brpo_top_forward = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )

opts_br = ContinuationPar(p_min = 2.1,p_max = 8.0, ds = -0.0005, dsmin = 0.00005, dsmax = 0.005, n_inversion = 100, nev = 3)
opts_po_cont = ContinuationPar(opts_br, max_steps = 1000, tol_stability = 1e-2)
brpo_top_back = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    argspo...
    )


pal= distinguishable_colors(
	6+1 ,
	[RGB(1, 1, 1), RGB(0, 0, 0)],
	dropseed = true,
);
pal[7] = colorant"black";

pal = pal[[1,4,3,7, 5, 2, 6]]

bif_diagram = plot(brpo_cycle.param[1:796],brpo_cycle.min[1:796], 
label = "CTLN cycle", color = pal[5],xtickfontsize=6,ytickfontsize=6, xguidefontsize=6, yguidefontsize=6, legendfontsize=6, titlefontsize=6, size = (200, 150), grid = false, fontfamily="helvetica", xlimits = (0,8))
plot!(bif_diagram, brpo_cycle.param[1:796],brpo_cycle.max[1:796], label = "", color = pal[5])
plot!(bif_diagram, brpo_cycle.param[797:end],brpo_cycle.min[797:end],  label = "",  color = pal[5],  linestyle = :dash )
plot!(bif_diagram, brpo_cycle.param[797:end],brpo_cycle.max[797:end], label = "", color = pal[5], linestyle = :dash)
plot!(bif_diagram, brpo_ei.param,brpo_ei.min, label = "E/I cycle", color = pal[7] )
plot!(bif_diagram, brpo_ei.param,brpo_ei.max, label = "", color = pal[7])
hline!(bif_diagram, [FP[4]], color = colorant"black", linestyle = :dash, label = "Fixed point")

plot!(bif_diagram, brpo_bottom_forward.param,brpo_bottom_forward.min, label = "Mixed", color = pal[6], linestyle = :dash )
plot!(bif_diagram, brpo_bottom_forward.param,brpo_bottom_forward.max, label = "",  color = pal[6], linestyle = :dash)
plot!(bif_diagram, brpo_top_forward.param,brpo_top_forward.min, label = "",  color = pal[6], linestyle = :dash )
plot!(bif_diagram, brpo_top_forward.param,brpo_top_forward.max, label = "",  color = pal[6], linestyle = :dash)

plot!(bif_diagram, brpo_bottom_back.param, brpo_bottom_back.min, label = "Mixed",  color = pal[6], linestyle = :dash )
plot!(bif_diagram, brpo_bottom_back.param, brpo_bottom_back.max, label = "",  color = pal[6], linestyle = :dash)
plot!(bif_diagram, brpo_top_back.param, brpo_top_back.min, label = "",  color = pal[6], linestyle = :dash )
plot!(bif_diagram, brpo_top_back.param, brpo_top_back.max, label = "",  color = pal[6], linestyle = :dash)
Plots.xlabel!("tau_I/tau_E")
Plots.ylabel!("x_I")
savefig("../results/plots/bif_diagram.pdf")


ctln_cycle = get_periodic_orbit(brpo_cycle, 1)
brpo_cycle.param[1]

projplot_ei = plot(size = (150, 150), legend = false)
projplot_ctln = plot(size = (150, 150), legend = false)

plot(ctln_cycle, palette = palette, legend = false, size = (200, 150),  ylimits = (0,1))
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/inhib_fast_ctln.pdf")


ctln_cycle = get_periodic_orbit(brpo_cycle, 144)
brpo_cycle.param[144]
plot(ctln_cycle, palette = pal, legend = false, size = (200, 150),  ylimits = (0,1))


init = ctln_cycle[1:3,1]
prob_de = ODEProblem(ctln!, init, (0,11.), par_cycle)
sol = DifferentialEquations.solve(prob_de, Rodas5())

plot!(sol, palette = pal[1:4], alpha = 0.5, width = 2)
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/inhib_equal_ctln.pdf")

ctln_cycle = get_periodic_orbit(brpo_cycle,439)
brpo_cycle.param[439]
plot(ctln_cycle, palette = palette, legend=false, size = (200, 150),  ylimits = (0,1))
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/e_i_3_ctln.pdf")

projpts = proj_ctln * ctln_cycle[:,:]
plot!(projplot_ctln, projpts[1,:], projpts[2,:], color = :red)
projpts = proj_ei * ctln_cycle[:,:]
plot!(projplot_ei, projpts[1,:], projpts[2,:], color = :red)

ei_cycle = get_periodic_orbit(brpo_ei, 202)
brpo_ei.param[202]
plot(ei_cycle, palette = palette, legend=false, size = (200, 150),  ylimits = (0,1))
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/e_i_3.pdf")

projpts = proj_ctln * ei_cycle[:,:]
plot!(projplot_ctln, projpts[1,:], projpts[2,:], color = :black, seriestype = :scatter)
projpts = proj_ei * ei_cycle[:,:]
plot!(projplot_ei, projpts[1,:], projpts[2,:],  color = :black)


iter = 313
mixed_cycle = get_periodic_orbit(brpo_top_back, iter)
plot(mixed_cycle, palette = palette, legend=false, size = (200, 150),  ylimits = (0,1))
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/mixed_cycle.pdf")

projpts = proj_ctln * mixed_cycle[:,:]
plot!(projplot_ctln, projpts[1,:], projpts[2,:], color = :purple)
projpts = proj_ei * mixed_cycle[:,:]
plot!(projplot_ei, projpts[1,:], projpts[2,:], color = :purple)

projfp = proj_ctln * FP
plot!(projplot_ctln, [projfp[1]], [projfp[2]], color = :gray, seriestype = :scatter, markersize = 5)

projfp = proj_ei * FP
plot!(projplot_ei, [projfp[1]], [projfp[2]], color = :gray, seriestype = :scatter, markersize = 5)

Plots.xlabel!(projplot_ctln, "x_1 - x_2")
Plots.ylabel!(projplot_ctln, "x_1 + x_2 - 2x_3")
savefig("../results/plots/proj_ctln.pdf")

Plots.xlabel!(projplot_ei, "x_1 + x_2 + x_3")
Plots.ylabel!(projplot_ei, "x_I")
savefig("../results/plots/proj_ei.pdf")
iter  = 800

par_btwn_further = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = brpo_cycle[iter].param, b = 1.0, a = 0.0)

ctln_cycle = get_periodic_orbit( brpo_cycle, iter)

plot(ctln_cycle)

init = ctln_cycle[:,1]

prob = BifurcationProblem(cycle_3!, init, par_btwn_further, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2], u = x[3]))

prob_de = ODEProblem(cycle_3!, init, (0,30.), par_btwn_further)

sol = DifferentialEquations.solve(prob_de, Rodas5())

plot(sol)

delta_x0 = [1.0; 1.0; 1.0; 0.0]
eps = .044798981786010157873
init_1 = init + eps*delta_x0
par_btwn_further = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = brpo_cycle[iter].param, b = 1.0, a = 0.0)

prob_de = ODEProblem(cycle_3!, init_1, (0,70.), par_btwn_further)
sol = DifferentialEquations.solve(prob_de, Rodas5())
plot(sol, palette = pal[1:4], legend = false,  size = (200, 150))
Plots.xlabel!("Time")
Plots.ylabel!("Firing rate")
savefig("../results/plots/quasiperiodic.pdf")

projpts = proj_ctln* hcat(sol.u[:,:]...)
projplot_weird = plot(size = (150, 150), legend = false)
plot!(projplot_weird, projpts[1,:], projpts[2,:], legend = false, color = pal[4])
savefig("../results/plots/quasiperiodic_proj.pdf")

stab_list = [] 
tau_Is = 1:.001:3
for tau_I in tau_Is
	pars = (wEE_self =1.5, wEE = .75 , wEI  = -2.25, wIE = 2.0, wII = -2.0, tau_E = 1.0, tau_I = tau_I, b = 1.0, a = 0.0)
	eigs = eigs_fp(pars)
	stab = minimum(real(eigs.values))
	push!(stab_list, stab)
end
plot(stab_list)
tau_Is[argmin(abs.(stab_list))]
