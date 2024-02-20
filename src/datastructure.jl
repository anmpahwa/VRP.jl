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
    Route(iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, iˢ::Int, iᵉ::Int, tˢ::Float64, tᵉ::Float64, θ̲::Float64, θ::Float64, ω::Float64, δ::Float64, n::Int, l::Float64)

A `Route` is a connection between nodes with vehicle index `iᵛ`, depot node index 
`iᵈ`, centroid coordinates `(x,y)`, start node index `iˢ`, end node index `iᵉ`,
start time `tˢ`, end time `tᵉ`, threshold vehicle on-arrival tank status `θ̲`, 
vehicle on-arrival tank status `θ`, post-departure fuel re-fueled `ω`, re-fuel 
detour length `δ`, total customers served `n`, and total length `l`.
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
    θ̲::Float64                                                                      # Threshold vehicle on-arrival tank status
    θ::Float64                                                                      # Vehicle on-arrival tank status at the depot node
    ω::Float64                                                                      # Fuel re-fueled post-departure from the depot node
    δ::Float64                                                                      # Re-fuel detour length
    n::Int                                                                          # Total customers served
    l::Float64                                                                      # Total length
end



"""
    Vehicle(iᵛ::Int, jᵛ::Int, iᵈ::Int, qᵛ::Float64, ωᵛ::Float64, lᵛ::Float64, sᵛ::Float64, τᶜ::Float64, τʷ::Float64, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64, r::Route)

A `Vehicle` is a mode of delivery with index `iᵛ`, type index `jᵛ`, depot node 
index `iᵈ`, package capacity `qᵛ`, tank capacity `ωᵛ`, range `lᵛ`, speed `sᵛ`, 
customer node parking time `τᶜ`, driver working-hours `τʷ`, operational cost `πᵈ` 
per unit distance and `πᵗ` per unit time, fixed cost `πᶠ`, and associated route `r`.
"""
mutable struct Vehicle
    iᵛ::Int                                                                         # Vehicle index
    jᵛ::Int                                                                         # Vehicle type index
    iᵈ::Int                                                                         # Depot node index
    qᵛ::Float64                                                                     # Package capacity
    ωᵛ::Float64                                                                     # Tank capacity
    lᵛ::Float64                                                                     # Range
    sᵛ::Float64                                                                     # Speed
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
    FuelStationNode(iⁿ::Int, jⁿ::Int, x::Float64, y::Float64, τᵛ::Float64, πᵒ::Float64, πᶠ::Float64, ω::Float64)

A `FuelStationNode` is a re-fueling point on the graph at `(x,y)` with index `iⁿ`, 
associated vehicle type `jⁿ`, re-fueling rate `τᵛ`, operational cost `πᵒ` per unit 
fuel re-fueled, fixed cost `πᶠ`, and fuel re-fueled `ω`.

"""
mutable struct FuelStationNode <: Node
    iⁿ::Int                                                                         # Fuel Station node index
    jⁿ::Int                                                                         # Fuel Station type (consistent with vehicle type)
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    τᵛ::Float64                                                                     # Vehicle re-fueling rate
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    ω::Float64                                                                      # Fuel re-fueled
end
"""
    DepotNode(iⁿ::Int, x::Float64, y::Float64, tˢ::Float64, tᵉ::Float64, F::Vector{FuelStationNode}, V::Vector{Vehicle}, n::Int)

A `DepotNode` is an origin point on the graph at `(x,y)` with index `iⁿ`, working-
hours start time `tˢ` and end time `tᵉ`, operational cost `πᵒ` per customer, fixed 
cost `πᶠ`, set of nearest fuel station nodes for every vehicle type `F`, fleet of 
vehicles `V`, and customers served `n`.
"""
mutable struct DepotNode <: Node
    iⁿ::Int                                                                         # Depot node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    tˢ::Float64                                                                     # Working-hours start time
    tᵉ::Float64                                                                     # Working-hours end time
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    F::Vector{FuelStationNode}                                                      # Nearest fuel station nodes
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
    n::Int                                                                          # Customers served
end
"""
    CustomerNode(iⁿ::Int, jⁿ::Int, x::Float64, y::Float64, qᶜ::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, F::Vector{FuelStationNode}, iᵗ::Int, iʰ::Int, tᵃ::Float64, tᵈ::Float64, q::Float64, θ̲::Float64, θ::Float64, ω::Float64, δ::Float64, r::Route)

A `CustomerNode` is a source/sink point on the graph at `(x,y)` with index `iⁿ`, 
associated delivery/pickup node index `jⁿ`, demand `qᶜ`, service time `τᶜ`, earliest 
service time `tᵉ`, latest service time `tˡ`, set of nearest fuel station nodes for 
every vehicle type `F`, tail node index `iᵗ`, head node index `iʰ`, vehicle arrival 
time `tᵃ` and departure time `tᵈ`, vehicle load `q` on-arrival if delivery node else 
on-departure if pickup node, threshold vehicle on-arrival tank status `θ̲`, vehicle 
on-arrival tank status `θ`, post-departure fuel re-fueled `ω`, re-fuel detour length 
`δ`, and associated route `r`.
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
    q::Float64                                                                      # Vehicle on-arrival load if delivery node else on-departure node if pickup node
    θ̲::Float64                                                                      # Threshold vehicle on-arrival tank status
    θ::Float64                                                                      # Vehicle on-arrival tank status at the customer node
    ω::Float64                                                                      # Fuel re-fueled post-departure from the customer node
    δ::Float64                                                                      # Re-fuel detour length
    r::Route                                                                        # Route visiting the customer node
end



"""
    Solution(F::Vector{FuelStationNode}, D::OffsetVector{DepotNode, Vector{DepotNode}}, C::OffsetVector{CustomerNode, Vector{CustomerNode}}, A::Dict{Tuple{Int,Int}, Arc}, πᶠ::Float64, πᵒ::Float64, πᵖ::Float64)

A `Solution` is a graph with fuel station nodes `F`, depot nodes `D`, customer nodes 
`C`, arcs `A`, fixed cost `πᶠ`, operational cost `πᵒ`, and penalty `πᵖ`.
"""
mutable struct Solution
    F::Vector{FuelStationNode}                                                      # Vector of fuel station nodes
    D::OffsetVector{DepotNode, Vector{DepotNode}}                                   # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int,Int}, Arc}                                                    # Set of arcs
    πᶠ::Float64                                                                     # Fixed cost
    πᵒ::Float64                                                                     # Opertaional cost
    πᵖ::Float64                                                                     # Penalty
    Solution(F, D, C, A) = new(F, D, C, A, 0., 0., 0.)
end