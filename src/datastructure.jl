"""
    Arc(iᵗ::Int, iʰ::Int, l::Float64)

An `Arc` is a connection between a tail node with index `iᵗ` and a head node with 
index `iʰ` spanning length `l`.
"""
struct Arc
    iᵗ::Int                                                                         # Tail node index
    iʰ::Int                                                                         # Head node index
    l::Float64                                                                      # Length
end



"""
    Route(iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, iˢ::Int, iᵉ::Int, tˢ::Float64, tᵉ::Float64, n::Int, l::Float64)

A `Route` is a connection between `n` nodes spanning length `l`, with vehicle index 
`iᵛ`, depot node index `iᵈ`, centroid coordinates `(x, y)`, start node index `iˢ` 
and end node index `iᵉ`, start time `tˢ` and end time `tᵉ`. 
"""
mutable struct Route
    iᵛ::Int                                                                         # Vehicle index
    iᵈ::Int                                                                         # Depot node index
    x::Float64                                                                      # Centroid x-coordinate
    y::Float64                                                                      # Centroid y-coordinate
    iˢ::Int                                                                         # Start node index
    iᵉ::Int                                                                         # End node index
    tˢ::Float64                                                                     # Start time
    tᵉ::Float64                                                                     # End time
    n::Int                                                                          # Customers served
    l::Float64                                                                      # Length
end



"""
    Vehicle(iᵛ::Int, jᵛ::Int, iᵈ::Int, qᵛ::Float64, lᵛ::Float64, sᵛ::Float64, θˡ::Float64, θᵘ::Float64, τᶜ::Float64, τʷ::Float64, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64, r::Route)

A `Vehicle` is a mode of delivery with index `iᵛ`, type index `jᵛ`, depot node index 
`iᵈ`, capacity `qᵛ`, range `lᵛ`, speed `sᵛ`, tank status lower threshold `θˡ` and 
upper threshold `θᵘ`, parking time `τᶜ` at customer node, driver working-hours `τʷ`, 
operational cost `πᵈ` per unit distance and `πᵗ` per unit time, fixed cost `πᶠ`, and 
associated vehicle route `r`.
"""
mutable struct Vehicle
    iᵛ::Int                                                                         # Vehicle index
    jᵛ::Int                                                                         # Vehicle type index
    iᵈ::Int                                                                         # Depot node index
    qᵛ::Float64                                                                     # Capacity
    lᵛ::Float64                                                                     # Range
    sᵛ::Float64                                                                     # Speed
    θˡ::Float64                                                                     # Tank status lower threshold
    θᵘ::Float64                                                                     # Tank status upper threshold
    τᶜ::Float64                                                                     # Parking time at customer node
    τʷ::Float64                                                                     # Driver working-hours duration
    πᵈ::Float64                                                                     # Distance-based operational cost
    πᵗ::Float64                                                                     # Time-based operational cost
    πᶠ::Float64                                                                     # Fixed cost
    r::Route                                                                        # Vehicle route
end



"""
    Node

A `Node` is a point on the graph.
"""
abstract type Node end
"""
    DepotNode(iⁿ::Int, x::Float64, y::Float64, tˢ::Float64, tᵉ::Float64, V::Vector{Vehicle}, n::Int)

A `DepotNode` is an origin point on the graph at `(x,y)` with index `iⁿ`, capacity 
`qᵈ`, working-hours start time `tˢ` and end time `tᵉ`, and a fleet of vehicles `V` 
serving `n` customers.
"""
mutable struct DepotNode <: Node
    iⁿ::Int                                                                         # Depot node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    tˢ::Float64                                                                     # Working-hours start time
    tᵉ::Float64                                                                     # Working-hours end time
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
    n::Int                                                                          # Customers served
end
"""
    CustomerNode(iⁿ::Int, jⁿ::Int, x::Float64, y::Float64, qᶜ::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, F::Vector{FuelStationNode}, iᵗ::Int, iʰ::Int, tᵃ::Float64, tᵈ::Float64, θ::Float64, q::Float64, r::Route)

A `CustomerNode` is a source/sink point on the graph at `(x,y)` with index `iⁿ`, 
associated delivery or pickup node index `jⁿ`, demand `qᶜ`, service time `τᶜ`, 
earliest service time `tᵉ` and latest service time `tˡ`, set `F` with the nearest 
fuel station for every vehicle type, tail node index `iᵗ` and head node index `iʰ`, 
vehicle arrival time `tᵃ` and departure time `tᵈ`, on-arrival vehicle tank status 
`θ` and load `q` on route `r`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int                                                                         # Customer node index
    jⁿ::Int                                                                         # Associated pickup/delivery node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location on the y-axis
    qᶜ::Float64                                                                     # Demand
    τᶜ::Float64                                                                     # Service time
    tᵉ::Float64                                                                     # Earliest service time
    tˡ::Float64                                                                     # Latest service time
    F::Vector{FuelStationNode}                                                      # Nearest fuel station nodes
    iᵗ::Int                                                                         # Tail (predecessor) node index
    iʰ::Int                                                                         # Head (successor) node index
    tᵃ::Float64                                                                     # Vehicle arrival time
    tᵈ::Float64                                                                     # Vehicle departure time
    θ::Float64                                                                      # Vehicle tank status on arrival
    q::Float64                                                                      # Vehicle load on arrival
    r::Route                                                                        # Route visiting the customer node
end
"""
    FuelStationNode(iⁿ::Int, jⁿ::Int, x::Float64, y::Float64, τᵛ::Float64, πᵒ::Float64, πᶠ::Float64, q::Float64)

A `FuelStationNode` is a re-fueling point on the graph at `(x,y)` with index `iⁿ`, 
associated with vehicle type `jᵛ`, with a re-fueling rate of `τᵛ`, operational cost 
`πᵒ` per unit fuel consumed re-fueling, fixed cost `πᶠ`, and amount re-fueled `q`.

"""
mutable struct FuelStationNode <: Node
    iⁿ::Int                                                                         # Fuel Station node index
    jⁿ::Int                                                                         # Fuel Station type (consistent with vehicle type)
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    τᵛ::Float64                                                                     # Vehicle re-fueling rate
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    q::Float64                                                                      # Amount re-fueled
end



"""
    Solution(D::Vector{DepotNode}, C::OffsetVector{CustomerNode, Vector{CustomerNode}}, F::OffsetVector{FuelStationNode, Vector{FuelStationNode}}, A::Dict{Tuple{Int,Int}, Arc}, πᶠ::Float64, πᵒ::Float64, πᵖ::Float64)

A `Solution` is a graph with depot nodes `D`, customer nodes `C`, fuel station nodes
`F`, arcs `A`, fixed cost `πᶠ`, operational cost `πᵒ`, and penalty `πᵖ`.
"""
mutable struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    F::OffsetVector{FuelStationNode, Vector{FuelStationNode}}                       # Vector of fuel station nodes
    A::Dict{Tuple{Int,Int}, Arc}                                                    # Set of arcs
    πᶠ::Float64                                                                     # Fixed cost
    πᵒ::Float64                                                                     # Opertaional cost
    πᵖ::Float64                                                                     # Penalty
    Solution(D, C, F, A) = new(D, C, F, A, 0., 0., 0.)
end