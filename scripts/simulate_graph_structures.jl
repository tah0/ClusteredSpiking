# simulate_graph_structures.jl
#
# Simulates three graph structures (ctln_cycle, fix_pt, indep_set) using the
# clustered spiking network model. For each structure:
#   - Saves time, population rates, and spike data to results/data/
#   - Saves a graph connectivity plot and a spike raster to results/plots/

using Compose, Random, CSV, DelimitedFiles
include("../src/csn_ctln.jl")

mkpath("../results/plots")
mkpath("../results/data")
Random.seed!(1)

# ── Shared simulation parameters ─────────────────────────────────────────────
NE      = 1000
NI      = 1000
pEE     = 0.2
pEI     = 0.8
pIE     = 0.8
pII     = 0.8
pE_self = 0.8
dt      = 1
τ_I     = 20
τ_E     = 40

# ── Helper functions ──────────────────────────────────────────────────────────

# Save time vector, rates matrix (T×n), and spike DataFrame
function save_sim_data(name, t, rates, spikes)
    writedlm("../results/data/$(name)_t.csv", collect(t), ',')
    writedlm("../results/data/$(name)_rates.csv", rates', ',')
    CSV.write("../results/data/$(name)_spikes.csv", spikes)
    println("  Saved data: $(name)_t.csv, $(name)_rates.csv, $(name)_spikes.csv")
end

# Raster plot of a random 5% sample of neurons
function plot_and_save_raster(spikes, ids, pal, name; frac = 0.05, sz = (400, 200))
    N = length(ids)
    idxs = rand(N) .< frac
    plot_neur = sort(findall(idxs))
    isempty(plot_neur) && (println("  Warning: no neurons sampled for raster"); return)
    df = copy(spikes[∈(plot_neur).(Int.(spikes.neurons)), :])
    df.yval = [findfirst(==(Int(neuron)), plot_neur) for neuron in df.neurons]
    gr(size = sz, legend = false, markerstrokewidth = 0, markersize = 30)
    scatter(df.spktimes, df.yval,
            color = pal[ids[Int.(df.neurons)]],
            markersize = 0.3, markerstrokewidth = 0,
            legend = false, markershape = :circle, size = sz)
    xlabel!("Time (ms)")
    ylabel!("Neuron")
    title!(name)
    savefig("../results/plots/$(name)_raster.png")
    println("  Saved plot: $(name)_raster.png")
end

# Graph connectivity plot using a small instantiation (NE=50, NI=50)
function plot_and_save_connectivity(sA, pal, name; NE_p = 50, NI_p = 50)
    W_p, _, lab_p = random_from_sA(sA, NE = NE_p, NI = NI_p,
                                    pEE = 0.01, pEI = 0.005,
                                    pIE = 0.005, pII = 0.5, pE_self = 0.5)
    A = W_p .!= 0
    g = DiGraph(A')
    p = gplot(g, nodefillc = pal[lab_p], arrowlengthfrac = 0)
    draw(PNG("../results/plots/$(name)_graph.png"), p)
    println("  Saved plot: $(name)_graph.png")
end


# ── 1. ctln_cycle ─────────────────────────────────────────────────────────────
println("\n=== Simulating: ctln_cycle ===")

sA = [0 0 1 1;
      1 0 0 0;
      0 1 0 0;
      1 0 0 0]
n     = size(sA, 1)
T_max = 3000

pal = distinguishable_colors(n + 1, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed = true)
pal[n + 1] = colorant"black"

W_small, b_small = graph_to_weights_inhib(sA, θ = 0.1, WII = -3)
W, b, ids = random_from_sA(sA, NE = NE, NI = NI,
                             pEE = pEE, pEI = pEI, pIE = pIE, pII = pII,
                             pE_self = pE_self, ϵ = 0.25, δ = 0.5, θ = 0.1, WII = -3)

function stimfunction_small(t)
    stim = zeros(n + 1)
    if (t > 1000) && (t < 1040)
        stim[4] += 0.2
    elseif (t > 2000) && (t < 2040)
        stim[2] += 0.2
    end
    return stim
end
stimfunction_big(t) = stimfunction_small(t)[ids]

x0_small = [0.03, 0.005, 0.06, 0.0, 0.05]   # length n+1 for EI system
x0 = x0_small[ids]
v0 = W * x0 + b

t, v, spikes = sim_v_model_spiking(W, b, ids;
    v0 = v0, T_max = T_max, dt = dt,
    tau_I = τ_I, tau_E = τ_E,
    stimfunction = stimfunction_big, maxspikes = 10^7)

save_sim_data("ctln_cycle", t, max.(v, 0), spikes)
plot_and_save_raster(spikes, ids, pal, "ctln_cycle")
plot_and_save_connectivity(sA, pal, "ctln_cycle")


# ── 2. fix_pt ─────────────────────────────────────────────────────────────────
println("\n=== Simulating: fix_pt ===")

sA = [0 1 0 1 0 0;
      1 0 1 1 0 0;
      0 1 0 1 0 0;
      1 1 1 0 1 0;
      0 0 0 1 0 0;
      0 0 0 0 0 0]
n     = size(sA, 1)
T_max = 4000

pal = distinguishable_colors(n + 1, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed = true)
pal[n + 1] = colorant"black"

W_small, b_small = graph_to_weights_inhib(sA, θ = 0.1, WII = -3)
pts = fixpts(W_small, b_small)
W, b, ids = random_from_sA(sA, NE = NE, NI = NI,
                             pEE = pEE, pEI = pEI, pIE = pIE, pII = pII,
                             pE_self = pE_self, ϵ = 0.25, δ = 0.5, θ = 0.1, WII = -3)

function stimfunction_small(t)
    if (t > 1000) && (t < 1050)
        stim = zeros(n + 1); stim[4:5] .+= 0.25; return stim
    elseif (t > 2000) && (t < 2050)
        stim = zeros(n + 1); stim[1:2] .+= 0.25; return stim
    elseif (t > 3000) && (t < 3050)
        stim = zeros(n + 1); stim[3] += 0.25; return stim
    else
        return zeros(n + 1)
    end
end
stimfunction_big(t) = stimfunction_small(t)[ids]

x0 = pts.value[1][ids]
v0 = W * x0 + b

t, v, spikes = sim_v_model_spiking(W, b, ids;
    v0 = v0, T_max = T_max, dt = dt,
    tau_I = τ_I, tau_E = τ_E,
    stimfunction = stimfunction_big)

save_sim_data("fix_pt", t, max.(v, 0), spikes)
plot_and_save_raster(spikes, ids, pal, "fix_pt")
plot_and_save_connectivity(sA, pal, "fix_pt")


# ── 3. indep_set ──────────────────────────────────────────────────────────────
println("\n=== Simulating: indep_set ===")

sA = zeros(Int, 6, 6)
n     = size(sA, 1)
T_max = 3000

pal = distinguishable_colors(n + 1, [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed = true)
pal[n + 1] = colorant"black"

W_small, b_small = graph_to_weights_inhib(sA, θ = 0.1, WII = -3)
pts = fixpts(W_small, b_small)
W, b, ids = random_from_sA(sA, NE = NE, NI = NI,
                             pEE = pEE, pEI = pEI, pIE = pIE, pII = pII,
                             pE_self = pE_self, ϵ = 0.25, δ = 0.5, θ = 0.1, WII = -3)

function stimfunction_small(t)
    stim = zeros(n + 1)
    if     (t > 500)  && (t < 600);  stim[2] += 0.25
    elseif (t > 1000) && (t < 1100); stim[3] += 0.25
    elseif (t > 1500) && (t < 1600); stim[4] += 0.25
    elseif (t > 2000) && (t < 2100); stim[5] += 0.25
    elseif (t > 2500) && (t < 2600); stim[6] += 0.25
    end
    return stim
end
stimfunction_big(t) = stimfunction_small(t)[ids]

x0 = pts.value[1][ids]
v0 = W * x0 + b

t, v, spikes = sim_v_model_spiking(W, b, ids;
    v0 = v0, T_max = T_max, dt = dt,
    tau_I = τ_I, tau_E = τ_E,
    stimfunction = stimfunction_big)

save_sim_data("indep_set", t, max.(v, 0), spikes)
plot_and_save_raster(spikes, ids, pal, "indep_set")
plot_and_save_connectivity(sA, pal, "indep_set")

println("\nDone. All results saved to ../results/")
