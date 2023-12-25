# Vehicle Routing Problem (VRP)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/anmol1104/VRP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/anmol1104/VRP.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/anmol1104/VRP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/anmol1104/VRP.jl)

### capacitated vehicle routing problem with pickup and delivery, time-windows, and heterogeneous fleet of multi-route delivery vehicles

Given, a graph `G = (D, C, A)` with set of depots `D`, set of customer nodes `C`, and set of arcs `A`, the objective is to develop least cost routes from depot nodes using select vehicles such that every customer node is visited exactly once while accounting for customer service and time-window constraints; vehicle range, capacity, and working-hours constraints; and depot operations mandate and capacity constraints.

This package uses Adaptive Large Neighborhood Search (ALNS) algorithm to find an optimal solution for the location routing problem given an initial solution (here, developed using iterated clustering method) and ALNS optimization parameters,
- `j`     :   Number of segments in the ALNS
- `k`     :   Number of segments to reset ALNS
- `n`     :   Number of iterations in an ALNS segment
- `m`     :   Number of iterations in a local search
- `Ψᵣ`    :   Vector of removal operators
- `Ψᵢ`    :   Vector of insertion operators
- `Ψₗ`    :   Vector of local search operators
- `σ₁`    :   Score for a new best solution
- `σ₂`    :   Score for a new better solution
- `σ₃`    :   Score for a new worse but accepted solution
- `μ̲`     :   Minimum removal fraction
- `c̲`     :   Minimum customer nodes removed
- `μ̅`     :   Maximum removal fraction
- `c̅`     :   Maximum customer nodes removed
- `ω̅`     :   Initial temperature deviation parameter
- `τ̅`     :   Initial temperatureprobability parameter
- `ω̲`     :   Final temperature deviation parameter
- `τ̲`     :   Final temperature probability parameter
- `θ`     :   Cooling rate
- `ρ`     :   Reaction factor

The ALNS metaheuristic iteratively removes a set of nodes using,
- Random Customer Node Removal  : `:randomcustomer!`
- Random Route Removal          : `:randomroute!`
- Random Vehicle Removal        : `:randomvehicle!`
- Random Depot Removal          : `:randomdepot!` 
- Related Customer Node Removal : `:relatedcustomer!`
- Related Route removal         : `:relatedroute!`
- Related Vehicle Removal       : `:relatedvehicle!`
- Related Depot Removal         : `:relateddepot!`
- Worst Customer Node Removal   : `:worstcustomer!`
- Worst Route Removal           : `:worstroute!`
- Worst Vehicle Removal         : `:worstvehicle!`
- Worst Depot Removal           : `:worstdepot!`

and consequently inserts removed nodes using,
- Best Insertion           : `:best!`
- Precise Greedy Insertion : `:precise!`
- Perturb Greedy insertion : `:perturb!`
- Regret-two Insertion     : `:regret2!`
- Regret-three Insertion   : `:regret3!`

In every few iterations, the ALNS metaheuristic performs local search with,
- intra-move    : `:intramove!`
- inter-move    : `:intermove!`
- intra-swap    : `:intraswap!`
- inter-swap    : `:interswap!`
- intra-opt     : `:intraopt!`
- inter-opt     : `:interopt!`
- swapdepot     : `:swapdepot!`

See benchmark.jl for usage.

Additional removal, insertion, and local search methods can be defined.