using MultivariateStats


using ClusteredSpiking
default(titlefont = (7, "helvetica"), legendfontsize = 6, guidefont = (6, "helvetica"), xtickfont = (6, "helvetica"), ytickfont = (6, "helvetica"), grid = false, fontfamily = "helvetica")




sA = [0 0 1; 1 0 0; 0 1 0]
n = size(sA, 1)
pal= Colors.distinguishable_colors(
    n +1 ,
    [RGB(1.0, 1.0, 1.0), RGB(0.0, 0.0, 0.0)],
    dropseed = true,
);
pal[n+1] = colorant"black";



NE =1000;
NI =1000;
N = n*NE + NI;
pEE= .2;
pEI= .8;
pIE= .8;
pII= .8;
pE_self= .8;
τ_I = 20;
τ_E =  40;
dt = 1;
T_max = 1000;
T_stim = 500;
θ = 0.1

W, b, ids  =random_from_sA(sA, NE=NE, NI=NI, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5,  
WEE_self = 2 , θ=θ, WII = -3, WIE = 4);

W_small, b_small, _  =random_from_sA(sA, NE=1, NI=1, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5,  
WEE_self = 2 , θ=θ, WII = -3, WIE = 4);

function  stimfunction_small(t)
    stim = zeros(n+1)
    if (t > T_stim) 
        stim[n+1] += 1.5*θ
    end
    return stim 
end

stimfunction_big(t) = stimfunction_small(t)[ids]


x0_small= θ.*[.4; .00; 0.4; .8]

x0 = x0_small[ids]
v0 = W * x0 + b
t, v, spikes = sim_v_model_spiking(W, b, ids; v0 = v0, T_max = T_max, dt = dt, tau_I =τ_I, tau_E =τ_E, stimfunction = stimfunction_big, maxspikes = 2*10^6);
sol = solve_TLN_inhib(W_small, b_small, τ_I=τ_I, τ_E =τ_E, T_max = T_max, x0 = x0_small, stimfunction = stimfunction_small)
T = (0:dt:(T_max-dt));
r = max.(v,0)



p1 = plot(legend = false, size = (255, 200))
#plot!(sol,  lc = pal', legend = false, size = (255, 200), alpha = .5)
ylabel!("Rate (spikes/ms)")
xlabel!("Time (ms)")
vspan!([T_stim, T_max], color = colorant"gray", label = "stimulus", alpha = 0.25)
plot!(T, r', lc = pal')




sA = hcat(0)
n = size(sA, 1)
pal= distinguishable_colors(
    n +1 ,
    [RGB(1, 1, 1), RGB(0, 0, 0)],
    dropseed = true,
);
pal[n+1] = colorant"black";





W, b, ids  =random_from_sA(sA, NE=NE, NI=NI, pEE=pEE, pEI=pEI, pIE=pIE, pII=pII,  pE_self=pE_self, ϵ=.25, δ=.5,  
WEE_self = 2 , θ= θ, WII = -3, WIE = 4);

function  stimfunction_small(t)
    stim = zeros(n+1)
    if (t > 500) 
        stim[n+1] += 1.5 * θ
    end
    return stim 
end

stimfunction_big(t) = stimfunction_small(t)[ids]


x0_small= [.1; .1]
x0 = x0_small[ids]
v0 = W * x0 + b
t, v, spikes = sim_v_model_spiking(W, b, ids; v0 = v0, T_max = T_max, dt = dt, tau_I =τ_I, tau_E =τ_E, stimfunction = stimfunction_big, maxspikes = 2*10^6);
T = (0:dt:(T_max-dt));
r = max.(v,0)



p2= plot(legend = false, size = (255, 200))
#plot!(sol,  lc = pal', legend = false, size = (255, 200), alpha = .5)
ylabel!("Rate (spikes/ms)")
xlabel!("Time (ms)")
vspan!([T_stim, T_max], color = colorant"gray", label = "stimulus", alpha = 0.25)
plot!(T, r', lc = pal')

plot(p2,p1, layout = (2,1), size = (255, 350))
savefig("../results/plots/paradoxical.pdf")
