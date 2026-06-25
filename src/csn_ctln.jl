
function tln_function(W, b; tau = 1,  stimfunction)
    function update!(du, u, p, t)
        new_du = (1/tau) *(-u .+ max.(W * u .+ b .+ stimfunction(t) , 0))
        for i in eachindex(du)
            du[i] = new_du[i]
        end
    end
    return update!
end

function solve_TLN(W, b; x0 = b[1]*rand(size(W, 1)), T_max = 100, dtmax = 0.1, maxiters = 1e5, tau = 1,  stimfunction = t -> 0)
    f! = tln_function(W, b, tau = tau, stimfunction = stimfunction)
    tspan = (0.0, T_max)
    prob = ODEProblem(f!, x0, tspan, dtmax = dtmax, maxiters = maxiters)
    sol = solve(prob)
    return sol
end

function tln_function_inhib(W, b, τ_E, τ_I, NI , stimfunction)
    n = size(W,1) -NI
    function update!(du, u, p, t)
        new_du = -u .+ max.(W * u .+ b .+ stimfunction(t), 0)
        new_du[1:n] *= (1/τ_E)
        new_du[n+1:end] *= (1/τ_I)
        for i in eachindex(du)
            du[i] = new_du[i]
        end
    end
    return update!
end





function solve_TLN_inhib(W, b,NI=1; x0 = mean(b)*rand(size(W, 1)), T_max = 100, dtmax = 0.5, maxiters = 1e5, τ_E = 1, τ_I = .01, stimfunction = t -> 0)
    f! = tln_function_inhib(W, b, τ_E, τ_I, NI, stimfunction)
    tspan = (0.0, T_max)
    prob = ODEProblem(f!, x0, tspan, dtmax = dtmax, maxiters = maxiters)
    sol = solve(prob)
    return sol
end


function write_sa(sA, file, name )
    file = matopen(file, "w")
    write(file, "sA"*string(name), sA)
    close(file)
end


function graph_to_weights(sA; ϵ = 0.25, δ = 0.5, θ = 1)
    n = size(sA, 1)
    W = (-1 + ϵ) * sA + (-1 - δ) * (ones(n, n) - I - sA)
    b = θ * ones(n)
    return W, b
end

function graph_to_weights_inhib(sA; ϵ = 0.25, δ = 0.5, θ = 1, WEE_self = 1 + δ, WII = -1, WIE = 2)
    W_E_self = 1+δ
    W_EE = ϵ + δ 
    θ_E = θ
    θ_I = 0
    WEI = -WEE_self * (1-WII)/WIE
    WEE = (WEE_self - 1 + ϵ)
    n = size(sA, 1)
    W = WEE * sA + WEE_self * I
    W = hcat(W, zeros(n))
    W = vcat(W, zeros(n+1)')
    W[n+1, :] .= WIE
    W[:, n+1] .= WEI
    W[n+1, n+1] = WII
    b = vcat(θ_E * ones(n), θ_I)
    return W, b
end


function plot_sol(sol)
    t = sol.t
    u = sol.u
    u = hcat(u...)'
    n = size(u, 2)
    if n <= 20
        labels = ["neuron " * string(i) for i = 1:n]
        labels = reshape(labels, 1, n)
        plot_legend = :outertopleft
    else
        labels = nothing
        plot_legend = false
    end
    colors = distinguishable_colors(
        n + 2,
        [RGB(1, 1, 1), RGB(0, 0, 0)],
        dropseed = true,)
    return Plots.plot(
        t,
        u,
        label = labels,
        legend = plot_legend,
        palette = colors[1:n],
        size = (min(4 * length(t), 1000)  , 300))
end

n_digits(n) = Int(floor(log(10,n)))+1


function make_names(n)
    l_max = n_digits(n)
    names = ["" for i = 1:n]
    for i in 1:n
        names[i] = string(i) * ' '^(l_max-n_digits(i))
    end
    return names
end




function running_average(x; bins = 10)
    T = length(x)
    x_pad = zeros(T + bins)
    x_pad[1:T] += x
    x_avg = [mean(x_pad[i:i+bins]) for i=1:T]
    return x_avg
end


function graph_plot(sA)
    g = DiGraph(sA')
    n = size(sA, 1)
    colors = distinguishable_colors(
        n + 2,
        [RGB(1, 1, 1), RGB(0, 0, 0)],
        dropseed = true,
    )
    p = graphplot(
        g,
        nodecolor = colors,
        names = make_names(n),
        nodeshape = :circle,
        method = :circular,
        node_size = 0.3,
        legend = false,
        arrow = :closed,
        curvature_scalar= 0
    )
    return p
end

function graph_to_plot(
    sA;
    ϵ = 0.25,
    δ = 0.5,
    θ = 1,
    x0 = 2 * rand(size(sA, 1)) / size(sA, 1),
    T_max = 100,
    total_pop = false,
    graph = (size(sA,1) < 10), 
)
    n = size(sA, 1)
    W, b = graph_to_weights(sA, ϵ = ϵ, δ = δ, θ = θ)
    sol = solve_TLN(W, b, x0 = x0, T_max = T_max)
    p = plot_sol(sol)
    if total_pop
        u = sol.u
        u = hcat(u...)'
        tpa = sum(u, dims = 2)
        Plots.plot!(p, sol.t, tpa, label = "total", legend = :outertopleft)
        max_val = -θ / (-1 + ϵ)
        min_val = -θ / (-1 - δ)
        hline!([min_val, max_val], label = "bounds", legend = :outertopleft)
    end
    if graph == true
        q = graph_plot(sA)

        l = @layout [a{0.7w} b]
        final = Plots.plot(p, q, layout = l, size = (3 * T_max +300 , 200), bottom_margin = 10px,
                left_margin = 10px)
        return (final, sol, x0)
    else
        return (Plots.plot(p, size = (4 * T_max  , 200), bottom_margin = 10px,
                left_margin = 10px), sol, x0)
    end
end




function random_from_sA(sA; NE=10, NI=10, pEE=.5, pEI=.5, pIE=.5, pII=.9, pEE_noedge = 0.5, pE_self=.9, ϵ=.25, δ=.5, θ=1,  WEE_self = 1 + δ, WII = -1, WIE = 2)
    n = size(sA, 1)
    if NE == 1
       A, b = graph_to_weights_inhib(sA, ϵ = ϵ, δ = δ, θ = θ,  WEE_self =  WEE_self, WII = WII, WIE = WIE)
       return A, b, collect(1:n+1)
    end
    N = n * NE + NI
    #WE_self = (1+δ)/(pE_self*(NE-1))

    WEI = -WEE_self * (1-WII)/WIE
    WEE = (WEE_self - 1 + ϵ)
    WEE_noedge = 1 + δ - WEE_self



    WEE = WEE/(pEE*NE)
    WEE_self = WEE_self/(pE_self*(NE-1))
    WEE_noedge = WEE_noedge / (pEE_noedge*NE)
    θE = θ

    WII = WII/((NI-1)*pII)
    WIE = WIE/(NE*pIE)
    #WEI = (-1-δ)*(1-WII * pII* (NI-1))/(NE*NI*pEI*pIE*WIE)
    WEI = WEI/(NI * pEI)
    labels = vcat(vec(repeat(1:n, 1, 1NE)'), (n+1)*ones(Int, NI))

    A = zeros(N, N)
    for i = 1:n, j =1:n
        if sA[i,j] != 0
            A[labels .== i, labels .== j] .= (rand(NE, NE) .< pEE)*WEE
        elseif i ==j
            A[labels .== i, labels .== j] .=  (rand(NE, NE) .< pE_self)*WEE_self
        else
            A[labels .== i, labels .== j] .= (rand(NE, NE) .< pEE_noedge)*WEE_noedge
        end
    end
    A[labels .!= n+1, labels .== n+1] .= (rand(n*NE, NI) .< pEI)*WEI
    A[labels .== n+1,labels .!= n+1] .= (rand(NI, n*NE) .< pIE)*WIE
    A[labels .== n+1, labels .== n+1] .= (rand(NI, NI) .< pII)* WII
    
    b = vcat(θE *ones(n*NE), zeros(NI))
    A -= diagm(diag(A))
    return A, b, labels
    
end

function plot_sol_inhib(sol, NE, NI)
    t = sol.t
    u = sol.u
    u = hcat(u...)'
    n = size(u, 2)
    n_grp = Int((n-NI)/NE)
    if n <= 20
        labels = ["neuron " * string(i) for i = 1:n]
        labels = reshape(labels, 1, n)
        plot_legend = :outertopleft
    else
        labels = nothing
        plot_legend = false
    end
    colors = distinguishable_colors(
         n_grp+1,
         [RGB(1, 1, 1), RGB(0, 0, 0)],
         dropseed = true,)
    
    colors_long =  distinguishable_colors(
        n_grp*NE+NI,
        [RGB(1, 1, 1), RGB(0, 0, 0)],
        dropseed = true,)
    for i = 1:n_grp
        for j = 1:NE
            colors_long[(i-1)*NE + j] = colors[i]
        end
    end
    for j = 1:NI
        colors_long[n_grp*NE + j] = colors[n_grp+1]
    end
    p = Plots.plot(
        t,
        u,
        label = labels,
        legend = plot_legend,
        palette = colors_long,
        size = (min(4 * length(t), 600)  , 300))
    u_summary  = zeros(size(u,1), n_grp +1)
    for i = 1:n_grp
        u_summary[:, i] += sum(u[:,(i-1)*NE+1:i*NE], dims =2)
    end
    u_summary[:,n_grp+1] += sum(u[:,NE*n_grp+1:n], dims =2)
    q= Plots.plot(
        t,
        u_summary,
        label = labels,
        legend = plot_legend,
        palette = colors,
        size = (min(4 * length(t), 600) , 300))
    return plot(p,q, layout = (2,1))
end

function sim_tln_spiking(W, b, labels; x0 = rand(size(W, 1)), T_max = 100, dt = 0.1, tau_I = 1, tau_E = 1, tau_S = .1, maxspikes = 10^6, stimfunction = t -> 0)
    N = size(W, 1)
    n = length(unique(labels))
    T = Int64(floor(T_max/dt )) 
    T_s = Int64(floor(tau_S/dt))
    pop_rates = zeros(n, T)
    X_t = x0
    spikes = zeros(N, T_s)
    tau = ones(N)
    tau[labels .< max(labels...)] *= tau_E
    tau[labels .== max(labels...)]*= tau_I
    tpts = (0:dt:(T_max-dt))
    neurons = zeros(maxspikes)
    spktimes = zeros(maxspikes)
    spk_ind = 1
    @showprogress for t = 1:T
        input = copy(b)
        input .+= stimfunction(t*dt)
        if spk_ind < maxspikes
            for neuron in unique(labels)
                pop_rates[neuron, t] = sum(X_t[ids .== neuron])/(sum(ids.==neuron)) #just a record, does not feed back in anywhere 
            end
            if t >= T_s + 2
                dx_dt = (1 ./tau).*( -X_t .+  max.(W * sum(spikes[:, :], dims = 2)/tau_S .+ input ,0)) #spikes feed back 
                X_t.= X_t .+ dx_dt*dt
            else
                X_t = X_t
            end
            rates = X_t
            spikes[:, 1:T_s-1] = spikes[:, 2:T_s]
            spikes[:, T_s] .= 0
            for i = 1:N
                p = rates[i] * dt
                if rand() < p
                    spikes[i, T_s] += 1
                    neurons[spk_ind] += i
                    spktimes[spk_ind] += tpts[t]
                    spk_ind += 1
                    if spk_ind == maxspikes
                        break
                    end
                end
            end
        end
    end
    df = DataFrame(neurons = neurons[1:spk_ind], spktimes = spktimes[1:spk_ind])
    return tpts, pop_rates, df
end

function sim_v_model_spiking(W, b, labels; v0 = rand(size(W, 1)), T_max = 100, dt = 0.1, tau_I = 1, tau_E = 1, maxspikes = 10^6, stimfunction = t -> 0)
    N = size(W, 1)
    n = length(unique(labels))
    T = Int64(floor(T_max/dt )) 
    pop_V = zeros(n, T)
    V_t = v0
    spikes = zeros(N)
    tau = ones(N)
    tau[labels .< max(labels...)] *= tau_E
    tau[labels .== max(labels...)]*= tau_I
    tpts = (0:dt:(T_max-dt))
    neurons = zeros(maxspikes)
    spktimes = zeros(maxspikes)
    spk_ind = 1
    @showprogress for t = 1:T
        input = copy(b)
        input .+= stimfunction(t*dt)
        if spk_ind < maxspikes
            for neuron in unique(labels)
                pop_V[neuron, t] = sum(V_t[ids .== neuron])/(sum(ids.==neuron)) #just a record, does not feed back in anywhere 
            end
            dv_dt = (1 ./tau).*( -V_t .+  W * (spikes/dt) + input) #spikes feed back 
            V_t.= V_t .+ dv_dt*dt

            rates = max.(V_t, 0)
            spikes .= 0
            for i = 1:N
                p = rates[i] * dt
                if rand() < p
                    spikes[i] += 1
                    neurons[spk_ind] += i
                    spktimes[spk_ind] += tpts[t]
                    spk_ind += 1
                    if spk_ind == maxspikes
                        break
                    end
                end
            end
        end
    end
    df = DataFrame(neurons = neurons[1:spk_ind], spktimes = spktimes[1:spk_ind])
    return tpts, pop_V, df
end


function pop_rates(X, ids, dt)
    neurons = unique(ids)
    N = length(neurons)
    T = size(X, 2)
    result = zeros(N, T)
    for neuron in neurons
        X_loc = X[ids .== neuron, :]
        result[neuron, :] = sum(X_loc,  dims = 1)/(sum(ids.==neuron)*dt)
    end
    return result
end

function mean_by_pop(X, ids)
    neurons = unique(ids)
    N = length(neurons)
    T = size(X, 2)
    result = zeros(N, T)
    for neuron in neurons
        X_loc = X[ids .== neuron, :]
        result[neuron, :] = sum(X_loc,  dims = 1)/(sum(ids.==neuron))
    end
    return result
end

function fixpts(W, b)
    n = length(b)

    points = []
    supports = []
    stabilities = []
    for σ in powerset(1:n)
        τ = setdiff(1:n, σ)
        x_σ = (I - W[σ, σ])^(-1) * b[σ]
        x = zeros(n)
        x[σ] += x_σ
        y = (W * x + b)[τ]
        if all(>(0), x_σ ) & all( <=(0), y)
            append!(points, [x])
            append!(supports, [σ])
            if max(real.(eigvals(W[σ, σ] - I))...)<0
                append!(stabilities, true)
            else
                append!(stabilities, false)
            end
        end 
    end
    return DataFrame(support = supports, value = points, stability = stabilities)
end

function W_sig(W, σ, full = true)
    if full
        result = zero(W)
        result[σ, :] += W[σ, :]
        return result
    else
        return W[σ,σ]
    end
end

function b_sig(b, σ, full = true)
    if full
        result = zero(b)
        result[σ] += b[σ]
        return result
    else
        return b[σ]
    end
end

function switch_times(W, b, u, t; start=1, stop=size(u,1))
    times = []
    hyperplanes = []
    chambers = []
    last = sign.(max.(W * u[start,:] .+b, 0))
    append!(times, t[start])
    push!(hyperplanes,[])
    push!(chambers, last)
    for i = start:stop
        current = sign.(max.(W * u[i,:] .+b, 0))
        if current != last
            append!(times, t[i])
            diff = findall(abs.(current-last).>0)
            push!(hyperplanes, diff)
            push!(chambers, current)
            last = current
        end
    end
    return DataFrame(time = times, hyperplane = hyperplanes, chamber = chambers)
end

function v_to_r(v, W, b)
    return W^(-1) * (v .- b)
end

function r_to_v(r, W, b)
    return W * r .+ b
end

function filter_input(tau, b, s, T)
    n = size(b,1)
    n_t = length(T)
    result = zeros(n, n_t)
    for (i,t) in enumerate(T)
        if i == 1
            result[:, i] = b
        else
            dx_dt = (1 ./ tau) .*( - result[:, i-1] + b + s(t))
            dt = T[i] - T[i-1]
            result[:, i] = result[:, i-1] + dx_dt * dt
        end
    end
    return result
end
