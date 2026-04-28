# ClusteredSpikingNetworksToCTLNS

Code to accompany the manuscript  [Diverse mean-field dynamics of clustered, inhibition-stabilized Hawkes networks via
combinatorial threshold-linear networks](https://arxiv.org/pdf/2506.06234).

To activate the project environment, in a Julia Repl, run 

```Julia
] activate . 
```
To generate Figures 1 and 2, 

```Julia
cd("scripts/")
include("examples_fig_1_2.jl")
```

To generate Figure 4:

```Julia
include("paradoxical.jl")
```

To generate Figure 5:

```Julia  
include("bifurcation_fp.jl")
 ```

 To generate Figure 6:

```Julia
include("bifurcation_3_cycle.jl")
```



