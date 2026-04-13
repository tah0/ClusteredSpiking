using Revise, Plots
using BifurcationKit
const BK = BifurcationKit
using DifferentialEquations

squareplus(x, a = 0.001) = (x + sqrt(x^2 + a))/2



function cycle_2!(du, X, p, t = 0)
	(;wEE_self, wEE, wEI, wIE, wII, tau_E, tau_I, b) = p
	x1, x2 = X
	du[1] =(1/tau_E)* (-x1 + squareplus(wEE_self*x1 + wEI * x2 + b))
	du[2] = (1/tau_I)*(-x2 + squareplus(wIE * x1 + wII * x2))
	du
end




par_cycle = (wEE_self =1.5, wEE = .75 , wEI  = -1.5, wIE = 2.0, wII = -1.0, tau_E = 1.0, tau_I = 1.0, b = 1.0)



recordFromSolution(x, p; k...) = (u1 = x[1], u2 = x[2])
z0 = [0.0, 0.0]
prob = BifurcationProblem(cycle_2!, z0, par_cycle,
# specify the continuation parameter
(@optic _.tau_I), record_from_solution = recordFromSolution)

opts_br = ContinuationPar(p_min = 0.5, p_max = 8.0, dsmax = 0.01,
# number of eigenvalues
nev = 2,
# maximum number of continuation steps
max_steps = 10000,)



br = continuation(prob, PALC(tangent=Bordered()), opts_br; normC = norminf)

par_cycle = (wEE_self =1.5, wEE = .75 , wEI  = -1.5, wIE = 2.0, wII = -1.0, tau_E = 1.0, tau_I = 3.9, b = 1.0)

z0 = [1.5, 1.5]

prob = BifurcationProblem(cycle_2!, z0, par_cycle, (@optic _.tau_I);
	record_from_solution = (x, p; k...) -> (x = x[1], y = x[2]))


prob_de = ODEProblem(cycle_2!, z0, (0,9.8), par_cycle)
sol = DifferentialEquations.solve(prob_de, Rodas5())

# arguments for periodic orbitsœ
# one function to record information and one
# function for plotting
args_po = (	record_from_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, p.p))
	end,
	plot_solution = (x, p; k...) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		arg = (marker = :d, markersize = 1)
		plot!(xtt.t, xtt[1,:]; label = "ext", arg..., k...)
		plot!(xtt.t, xtt[2,:]; label = "inh", arg..., k...)
		plot!(br; subplot = 1, putspecialptlegend = false)
		end,
	# we use the supremum norm
	normC = norminf)


probtrap, ci = BK.generate_ci_problem(PeriodicOrbitTrapProblem(M = 300),
	prob, sol, 9.8)

opts_br = ContinuationPar(p_min = 0.01, p_max = 8.0, ds =-0.002, dsmax = 0.01)
opts_po_cont = ContinuationPar(opts_br, max_steps = 1000, tol_stability = 1e-3)
brpo_fold = continuation(probtrap, ci, PALC(), opts_po_cont;
    verbosity = 3, plot = true,
    args_po...
    )

scene = plot(br, brpo_fold)
plot!(scene, brpo_fold.param,brpo_fold.min, label = "")
savefig("../results/plots/bifurcation_one_pop.pdf")