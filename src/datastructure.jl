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
    Route(iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, iˢ::Int, iᵉ::Int, tˢ::Float64, tᵉ::Float64, n::Int, q::Float64, l::Float64)

A `Route` is a connection between nodes, with vehicle index `iᵛ`, depot node index 
`iᵈ`, centroid coordinates `(x, y)`, start node index `iˢ` and end node index `iᵉ`, 
route start time `tˢ` and route end time `tᵉ`, customer served `n`, demand served 
`q`, and length `l`.
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
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Length
end



"""
    Vehicle(iᵛ::Int, jᵛ::Int, iᵈ::Int, qᵛ::Float64, lᵛ::Float64, sᵛ::Float64, ρᵛ::Float64, θˡ::Float64, θᵘ::Float64, τᶜ::Float64, τʷ::Float64, r::Route, tˢ::Float64, tᵉ::Float64, n::Int, q::Float64, l::Float64, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64)

A `Vehicle` is a mode of delivery with index `iᵛ`, type index `jᵛ`, depot node index 
`iᵈ`, capacity `qᵛ`, range `lᵛ`, speed `sᵛ`, fuel consumption rate `ρᵛ`, tank status 
lower threshold `θˡ` and upper threshold `θᵘ`, parking time `τᶜ` at customer node, 
driver working-hours `τʷ`, associated vehicle route `r`, initial departure time `tˢ`, 
final arrival time `tᵉ`, customers served `n`, demand served `q`, route length `l`, 
operational cost `πᵈ` per unit distance and `πᵗ` per unit time, and fixed cost `πᶠ`.
"""
mutable struct Vehicle
    iᵛ::Int                                                                         # Vehicle index
    jᵛ::Int                                                                         # Vehicle type index
    iᵈ::Int                                                                         # Depot node index
    qᵛ::Float64                                                                     # Capacity
    lᵛ::Float64                                                                     # Range
    sᵛ::Float64                                                                     # Speed
    ρᵛ::Float64                                                                     # Fuel consumption rate
    θˡ::Float64                                                                     # Tank status lower threshold
    θᵘ::Float64                                                                     # Tank status upper threshold
    τᶜ::Float64                                                                     # Parking time at customer node
    τʷ::Float64                                                                     # Driver working-hours duration
    r::Route                                                                        # Vehicle route
    tˢ::Float64                                                                     # Start time (departure time from the depot node)
    tᵉ::Float64                                                                     # End time (arrival time at the depot node)
    n::Int                                                                          # Customers served
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Route length
    πᵈ::Float64                                                                     # Distance-based operational cost
    πᵗ::Float64                                                                     # Time-based operational cost
    πᶠ::Float64                                                                     # Fixed cost
end



"""
    Node

A `Node` is a point on the graph.
"""
abstract type Node end
"""
    DepotNode(iⁿ::Int, x::Float64, y::Float64, tˢ::Float64, tᵉ::Float64, V::Vector{Vehicle}, n::Int, q::Float64, l::Float64, πᵒ::Float64, πᶠ::Float64)

A `DepotNode` is an origin point on the graph at `(x,y)` with index `iⁿ`, capacity 
`qᵈ`, working-hours start time `tˢ` and end time `tᵉ`, fleet of vehicles `V`, 
customer nodes visited `n`, demand served `q`, route length `l`, operational cost 
`πᵒ` per unit demand served, and fixed cost `πᶠ`.
"""
mutable struct DepotNode <: Node
    iⁿ::Int                                                                         # Depot node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    tˢ::Float64                                                                     # Working-hours start time
    tᵉ::Float64                                                                     # Working-hours end time
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
    n::Int                                                                          # Customers served
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Route length
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
end
"""
    CustomerNode(iⁿ::Int, jⁿ::Int, iᵛ::Int, iᵈ::Int, Iᶠ::Vector{Int}, x::Float64, y::Float64, qᶜ::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, r::Route, iᵗ::Int, iʰ::Int, tᵃ::Float64, tᵈ::Float64, n::Int, q::Float64, l::Float64)

A `CustomerNode` is a source/sink point on the graph at `(x,y)` with index `iⁿ`, 
associated delivery/pickup node index `jⁿ`, serviced on route `r` with vehicle 
index `iᵛ` and depot node index `iᵈ`, with nearest fueling station node indices `Iᶠ` 
associated with each vehicle type, demand `qᶜ`, service time `τᶜ`, earliest service 
time `tᵉ` and latest service time `tˡ`, tail node index `iᵗ` and head node index 
`iʰ`, vehicle arrival time `tᵃ` and departure time `tᵈ`, position `n` on the route, 
on-arrival vehicle load `q` and route length `l`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int                                                                         # Customer node index
    jⁿ::Int                                                                         # Associated pickup/delivery node index
    iᵛ::Int                                                                         # Vehicle index
    iᵈ::Int                                                                         # Depot node index
    Iᶠ::Vector{Int}                                                                 # Fuel Station node indices
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location on the y-axis
    qᶜ::Float64                                                                     # Demand
    τᶜ::Float64                                                                     # Service time
    tᵉ::Float64                                                                     # Earliest service time
    tˡ::Float64                                                                     # Latest service time
    r::Route                                                                        # Route visiting the customer node
    iᵗ::Int                                                                         # Tail (predecessor) node index
    iʰ::Int                                                                         # Head (successor) node index
    tᵃ::Float64                                                                     # Vehicle arrival time
    tᵈ::Float64                                                                     # Vehicle departure time
    n::Int                                                                          # Position on the route
    q::Float64                                                                      # Vehicle load on arrival
    l::Float64                                                                      # Route length on arrival
end
"""
    FuelStationNode(iⁿ::Int, jⁿ::Int, x::Float64, y::Float64, ρᵛ::Float64, n::Int, q::Float64, l::Float64, πᵒ::Float64, πᶠ::Float64)

A `FuelStationNode` is a re-fueling point on the graph at `(x,y)` with index `iⁿ`, 
re-fueling events `n`, amount re-fueled `q`, range restored `l`, operational cost 
`πᵒ` per unit fuel consumed re-fueling, and fixed cost `πᶠ`, that can re-fuel 
vehicles of type `jⁿ` at `ρᵛ` re-fueling rate.
"""
mutable struct FuelStationNode <: Node
    iⁿ::Int                                                                         # Fuel Station node index
    jⁿ::Int                                                                         # Fuel Station type (consistent with vehicle type)
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    ρᵛ::Float64                                                                     # Refueling rate
    n::Int                                                                          # Number of re-fueling events
    q::Float64                                                                      # Amount re-fueled
    l::Float64                                                                      # Range restored
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
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