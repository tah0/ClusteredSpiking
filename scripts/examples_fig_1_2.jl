using Compose
using Random
mkpath("../results/plots") 
include("../src/csn_ctln.jl")
default(titlefont = (7, "helvetica"), legendfontsize = 6, guidefont = (6, "helvetica"), xtickfont = (6, "helvetica"), ytickfont = (6, "helvetica"), grid = false, fontfamily = "helvetica")
Random.seed!(1)

Plots.scalefontsizes(1)

NE = 1000;
NI = 1000;
pEE= .2;
pEI= .8;
pIE= .8;
pII= .8;
pE_self= .8;
dt = 1;
τ_I = 20;
τ_E = 40;


sA = [0 0 1 1; 1 0 0 0; 0 1 0 0; 1 0 0 0]
n = size(sA)[1]
N = n*NE + NI;

palette= distinguishable_colors(
    n +1 ,
    [RGB(1, 1, 1), RGB(0, 0, 0)],
    dropseed = true,
)
palette[n+1] = colorant"black";
T_max = 3000;

W_small, b_small = graph_to_weights_inhib(sA,  θ=.1, WII = -3)
pts = fixpts(W_small, b_small)


x0_small_cycle = [0.03, 0.005, 0.06, 0.0, 0.05]


W_ctln, b_ctln = graph_to_weights(sA,  θ=.1)

x0_ctln_cycle = x0_small_cycle[1:4]

function  stimfunction_small(t)
    stim = zeros(n+1)
    if (t > 1000) & (t < 1040)
        stim[4] += .2
        return stim 
    elseif (t > 2000) & (t < 2040)
        stim[2] += .2
        return stim
    else
        return stim
    end
end

stimfunction_ctln(t) = stimfunction_small(t)[1:4]

filt_in_ctln  = filter_input([τ_E; τ_E; τ_E; τ_E], b_ctln, stimfunction_ctln,  collect(1:T_max))
stimfunction_filt_ctln(t) = filt_in_ctln[:,max(1, Int(ceil(t)))]
sol_ctln = solve_TLN(W_ctln, zeros(4),  tau = τ_E, T_max = T_max, x0 = x0_ctln_cycle, stimfunction = stimfunction_filt_ctln);
u_ctln = hcat(sol_ctln.u...)';
v_ctln = r_to_v(u_ctln', W_ctln, hcat(stimfunction_filt_ctln.(t for t in sol_ctln.t)...))'
r_ctln = max.(v_ctln, 0)
plot( sol_ctln.t, r_ctln[:,:], lc = palette', size=(400,200),  legend = false)
savefig("../results/plots/example_graph_ctln_cycle.pdf")

filt_in_ei_ctln  = filter_input([τ_E; τ_E; τ_E; τ_E; τ_I], b_small, stimfunction_small,  collect(1:T_max))
stimfunction_filt_ei_ctln(t) = filt_in_ei_ctln[:,max(1, Int(ceil(t)))]
sol_ei_ctln = solve_TLN_inhib(W_small,zeros(5),  τ_I =  τ_I, τ_E =  τ_E, T_max = T_max, x0 =x0_small_cycle, stimfunction = stimfunction_filt_ei_ctln);
u_ei_ctln = hcat(sol_ei_ctln.u...)';
v_ei_ctln = r_to_v(u_ei_ctln', W_small, hcat(stimfunction_filt_ei_ctln.(t for t in sol_ei_ctln.t)...))'
r_ei_ctln = max.(v_ei_ctln, 0)
plot( sol_ei_ctln.t, r_ei_ctln[:,:], lc = palette', size=(400,200),  legend = false)
savefig("../results/plots/example_graph_ei_ctln_cycle.pdf")


W, b, ids  =random_from_sA(sA, NE = NE, NI = NI, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5, θ=.1, WII = -3);


stimfunction_big(t) = stimfunction_small(t)[ids]


x0_cycle = x0_small_cycle[ids];

v0_cycle = W * x0_cycle + b
t, x, spikes = sim_v_model_spiking(W, b, ids; v0 = v0_cycle, T_max = T_max, dt = dt, tau_I =τ_I, tau_E =τ_E, stimfunction = stimfunction_big, maxspikes = 10^7);
r = max.(x, 0)
p1 = plot(t, r', lc = palette',  size=(300,150),  legend = false)
plot!( sol_ctln.t, r_ctln[:,:], lc = palette', alpha = .5, linewidth = 3)
title!("Rates in clustered spiking network vs. CTLN")

xlabel!("Time (ms)")
ylabel!("Firing rate (spikes/ms)")

p2 = plot(t, r', lc = palette',  size=(300,150),  legend = false)
plot!( sol_ei_ctln.t, r_ei_ctln[:,:], lc = palette', alpha = .5, linewidth = 3)
title!("Rates in clustered spiking network vs. EI-CTLN")
xlabel!("Time (ms)")
ylabel!("Firing rate (spikes/ms)")

plot(p2, p1, layout = (2, 1), size = (300, 300), legend = false)
savefig("../results/plots/example_graph_stim_spike_vs_rate.pdf")


idxs = rand(N) .< .05
plot_neur = sort(findall(idxs))
df =spikes[∈(plot_neur).(spikes.neurons), :]
df.yval = [findfirst(x -> x == neuron, plot_neur) for neuron in df.neurons]
gr( size=(400,200),legend=false,markerstrokewidth=0,markersize=30)
scatter(df.spktimes, df.yval, color = palette[ids[Int.(df.neurons)]], markersize = .3, markerstrokewidth=0, legend = false, markershape=:circle, size=(300,150))
savefig("../results/plots/example_graph_raster_stim.pdf")



NE_plot = 50
NI_plot = 50
W,  b, labels  = random_from_sA(sA,NE = NE_plot, NI = NI_plot,  pEE=.01, pEI=.001, pIE=.001, pII=.5, pE_self=.5)
A = W .!= 0
g = DiGraph(A')
N = size(A, 1)

colors = palette[labels]


p = gplot(g, nodefillc = colors, arrowlengthfrac = 0)

draw(PDF("../results/plots/example_graph.pdf"), p)



sA = [0 1 0 1 0 0;
      1 0 1 1 0 0;
      0 1 0 1 0 0;
      1 1 1 0 1 0;
      0 0 0 1 0 0;
      0 0 0 0 0 0]

NE = 1000;
NI = 1000;
T_max = 4000;

n = size(sA)[1]
N = n*NE + NI;


pal= distinguishable_colors(
    n +1 ,
    [RGB(1, 1, 1), RGB(0, 0, 0)],
    dropseed = true,
);
pal[n+1] = colorant"black";

W_small, b_small = graph_to_weights_inhib(sA,  θ=.1, WII = -3)
pts = fixpts(W_small, b_small)

W, b, ids  =random_from_sA(sA, NE = NE, NI = NI, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5, θ=.1, WII = -3);

function  stimfunction_small(t)
    if (t > 1000) & (t < 1050)
        stim = zeros(n+1)
        stim[4:5] .+= .25
        return stim
    elseif (t > 2000) & (t < 2050)
        stim = zeros(n+1)
        stim[1:2] .+= .25
        return stim
    elseif (t > 3000) & (t < 3050)
        stim = zeros(n+1)
        stim[3] += .25
        return stim
    else
        return zeros(n+1)
    end  
end

x0 = pts.value[1];

tau = ones(n+1) * τ_E
tau[n+1] = τ_I
filt_in_ei_ctln  = filter_input(tau, b_small, stimfunction_small,  collect(1:T_max))
stimfunction_filt_ei_ctln(t) = filt_in_ei_ctln[:,max(1, Int(ceil(t)))]
sol = solve_TLN_inhib(W_small, b_small, τ_I =  τ_I, τ_E =  τ_E, T_max = T_max, x0 =x0, stimfunction = stimfunction_filt_ei_ctln);
u = hcat(sol.u...)';
plot( sol.t, u[:,:], lc = pal',  legend = false)


x0_big = x0[ids];
stimfunction_big(t) = stimfunction_small(t)[ids]

v0_big = W * x0_big + b
t, v, spikes = sim_v_model_spiking(W, b, ids; v0 = v0_big, T_max = T_max, dt = dt, tau_I =τ_I, tau_E =τ_E,  stimfunction = stimfunction_big);
plot(t, v', lc = palette',  legend = false)
savefig("../results/plots/fix_pt_fig_stim.pdf")

idxs = rand(N) .< .05
plot_neur = sort(findall(idxs))
df =spikes[∈(plot_neur).(spikes.neurons), :]
df.yval = [findfirst(x -> x == neuron, plot_neur) for neuron in df.neurons]
scatter(df.spktimes, df.yval, color = pal[ids[Int.(df.neurons)]], markersize = .3, markerstrokewidth=0,  size=(300,150), legend=false)
savefig("../results/plots/fix_pt_fig_stim_raster.pdf")



NE_plot = 50
NI_plot = 50
W,  b, labels  = random_from_sA(sA,NE = NE_plot, NI = NI_plot,  pEE=.01, pEI=.005, pIE=.005, pII=.5, pE_self=.5)
A = W .!= 0
#A += rand(size(A)...) .< .0005
g = DiGraph(A')
N = size(A, 1)

colors = pal[labels]


p = gplot(g, nodefillc = colors, arrowlengthfrac = 0)

draw(PDF("../results/plots/fixed_point_fig_graph.pdf"), p)


      
sA = [0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0]

NE = 1000;
NI = 1000;
T_max = 3000;

n = size(sA)[1]
N = n*NE + NI;


pal= distinguishable_colors(
    n +1 ,
    [RGB(1, 1, 1), RGB(0, 0, 0)],
    dropseed = true,
)
pal[n+1] = colorant"black";

W_small, b_small = graph_to_weights_inhib(sA,  θ=.1, WII = -3)
pts = fixpts(W_small, b_small)

W, b, ids  =random_from_sA(sA, NE = NE, NI = NI, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5, θ=.1, WII = -3);

function  stimfunction_small(t)
    stim = zeros(n+1)
    if (t > 500) & (t < 600)
        stim[2] += .25
    elseif (t > 1000) & (t < 1100)
        stim[3] += .25
    elseif (t > 1500) & (t < 1600)
        stim[4] = .25
    elseif (t > 2000) & (t < 2100)
        stim[5] += .25
    elseif (t > 2500) & (t < 2600)
        stim[6] += .25
    end  
    return stim
end

x0 = pts.value[1];

tau = ones(n+1) * τ_E
tau[n+1] = τ_I
filt_in_ei_ctln  = filter_input(tau, b_small, stimfunction_small,  collect(1:T_max))
stimfunction_filt_ei_ctln(t) = filt_in_ei_ctln[:,max(1, Int(ceil(t)))]
sol = solve_TLN_inhib(W_small, b_small, τ_I =  τ_I, τ_E =  τ_E, T_max = T_max, x0 =x0, stimfunction = stimfunction_filt_ei_ctln);
u = hcat(sol.u...)';
plot( sol.t, u[:,:], lc = pal',  legend = false)


x0_big = x0[ids];
stimfunction_big(t) = stimfunction_small(t)[ids]

v0_big = W * x0_big + b
t, v, spikes = sim_v_model_spiking(W, b, ids; v0 = v0_big, T_max = T_max, dt = dt, tau_I =τ_I, tau_E =τ_E,  stimfunction = stimfunction_big);
plot(t, v', lc = palette',  legend = false)
savefig("../results/plots/fix_pt_fig_stim.pdf")

idxs = rand(N) .< .05
plot_neur = sort(findall(idxs))
df =spikes[∈(plot_neur).(spikes.neurons), :]
df.yval = [findfirst(x -> x == neuron, plot_neur) for neuron in df.neurons]
scatter(df.spktimes, df.yval, color = pal[ids[Int.(df.neurons)]], markersize = .3, markerstrokewidth=0,  size=(300,150), legend=false)
savefig("../results/plots/indep_set_fig_stim_raster.pdf")


NE_plot = 50
NI_plot = 50
W,  b, labels  = random_from_sA(sA,NE = NE_plot, NI = NI_plot,  pEE=.01, pEI=.005, pIE=.005, pII=.5, pE_self=.5)
A = W .!= 0
#A += rand(size(A)...) .< .002
g = DiGraph(A')
N = size(A, 1)

colors = pal[labels]


p = gplot(g, nodefillc = colors, arrowlengthfrac = 0 )

draw(PDF("../results/plots/indep_set_graph.pdf"), p)
