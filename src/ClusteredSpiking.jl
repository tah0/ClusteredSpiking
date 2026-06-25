module ClusteredSpiking

using DifferentialEquations
using LinearAlgebra
using Plots
using Graphs
using GraphPlot
using GraphRecipes
using Colors
using Cairo
using Fontconfig
using Plots.PlotMeasures
using StatsBase
using ProgressMeter
using Combinatorics
using DataFrames
using StatsPlots

export tln_function, solve_TLN, tln_function_inhib, solve_TLN_inhib,
       write_sa, graph_to_weights, graph_to_weights_inhib,
       plot_sol, n_digits, make_names, running_average,
       graph_plot, graph_to_plot, random_from_sA, plot_sol_inhib,
       sim_tln_spiking, sim_v_model_spiking,
       pop_rates, mean_by_pop, fixpts,
       W_sig, b_sig, switch_times, v_to_r, r_to_v, filter_input

include(joinpath(@__DIR__, "csn_ctln.jl"))

end # module ClusteredSpiking
