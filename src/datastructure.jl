"""
    Arc(iᵗ::Int, iʰ::Int, l::Float64)

An `Arc` is a connection between a tail node with index `iᵗ` and a head node with 
index `iʰ` spanning length `l`.
"""
struct Arc
    iᵗ::Int                                                                         # Tail node index
    iʰ::Int                                                                         # Head node index
    l::Float64                                                                      # Arc length
end



"""
    Route(iʳ::Int, iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, iˢ::Int, iᵉ::Int, θⁱ::Float64, θˢ::Float64, θᵉ::Float64, tⁱ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int, q::Float64, l::Float64, φ::Int64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, centroid coordinates `x, y` , start node index `iˢ`, end node index 
`iᵉ`, vehicle tank status `θⁱ`, `θˢ`, and `θᵉ` at route initiaition time `tⁱ`, start
time `tˢ`, and end time `tᵉ`, repsectively, slack time `τ`, nodes visited `n`, demand 
served from the depot `q`, route length `l`, and route initiaition status `φ`.
"""
mutable struct Route
    iʳ::Int                                                                         # Route index
    iᵛ::Int                                                                         # Vehicle index
    iᵈ::Int                                                                         # Depot node index
    x::Float64                                                                      # Centroid x-coordinate
    y::Float64                                                                      # Centroid y-coordinate
    iˢ::Int                                                                         # Route start node index
    iᵉ::Int                                                                         # Route end node index
    θⁱ::Float64                                                                     # Vehicle tank status at the initiation time
    θˢ::Float64                                                                     # Vehicle tank status at the start time
    θᵉ::Float64                                                                     # Vehicle tank status at the end time
    tⁱ::Float64                                                                     # Route initiation time
    tˢ::Float64                                                                     # Route start time
    tᵉ::Float64                                                                     # Route end time
    τ::Float64                                                                      # Route slack time
    n::Int                                                                          # Route total nodes visited
    q::Float64                                                                      # Route total demand served from the depot
    l::Float64                                                                      # Route length
    φ::Int64                                                                        # Route initialization status
end



"""
    Vehicle(iᵛ::Int, jᵛ::Int, iᵈ::Int, qᵛ::Float64, lᵛ::Float64, sᵛ::Float64, τᶠ::Float64, τᵈ::Float64, τᶜ::Float64, r̅::Int, τʷ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int, q::Float64, l::Float64, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ`, vehicle type index `jᵛ`, depot 
node index `iᵈ`, capacity `qᵛ`, range `lᵛ`, speed `sᵛ`, refueling time `τᶠ`, service 
time `τᵈ` at depot node per unit demand, parking time `τᶜ` at customer node, maximum 
routes permitted `r̅`, working-hours `τʷ`, initial departure time `tˢ`, final arrival 
time `tᵉ`, slack time `τ`, nodes visited `n`, demand served from the depot `q`, total 
route length `l`, operational cost `πᵈ` per unit distance and `πᵗ` per unit time, 
fixed cost `πᶠ`, and set of routes `R`.
"""
mutable struct Vehicle
    iᵛ::Int                                                                         # Vehicle index
    jᵛ::Int                                                                         # Vehicle type index
    iᵈ::Int                                                                         # Depot node index
    qᵛ::Float64                                                                     # Vehicle capacity
    lᵛ::Float64                                                                     # Vehicle range
    sᵛ::Float64                                                                     # Vehicle speed
    τᶠ::Float64                                                                     # Re-fueling time
    τᵈ::Float64                                                                     # Depot node service time per unit demand
    τᶜ::Float64                                                                     # Parking time at customer stop
    τʷ::Float64                                                                     # Vehicle working-hours duration
    r̅::Int                                                                          # Maximum number of vehicle routes permitted
    tˢ::Float64                                                                     # Vehicle start time (initial departure time from the depot node)
    tᵉ::Float64                                                                     # Vehicle end time (final arrival time at the depot node)
    τ::Float64                                                                      # Vehicle slack time
    n::Int                                                                          # Vehicle total nodes visited
    q::Float64                                                                      # Vehicle total demand served from the depot
    l::Float64                                                                      # Vehicle total route length
    πᵈ::Float64                                                                     # Vehicle operational cost (distance based)
    πᵗ::Float64                                                                     # Vehicle operational cost (time based)
    πᶠ::Float64                                                                     # Vehicle fixed cost
    R::Vector{Route}                                                                # Vector of vehicle routes
end



"""
    Node

A `Node` is a point on the graph.
"""
abstract type Node end
"""
    DepotNode(iⁿ::Int, x::Float64, y::Float64, qᵈ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int, q::Float64, l::Float64, πᵒ::Float64, πᶠ::Float64, V::Vector{Vehicle})

A `DepotNode` is an origin point on the graph at `(x,y)` with index `iⁿ`, capacity 
`qᵈ`, working-hours start time `tˢ` and end tme `tᵉ`, slack time `τ`, total nodes 
visited `n`, demand served `q`, total route length `l`, operational cost `πᵒ` per 
package, fixed cost `πᶠ`, operations mandate `φ`, and fleet of vehicles `V`.
"""
mutable struct DepotNode <: Node
    iⁿ::Int                                                                         # Depot node index
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    qᵈ::Float64                                                                     # Depot capacity
    tˢ::Float64                                                                     # Depot working-hours start time
    tᵉ::Float64                                                                     # Depot working-hours end time
    τ::Float64                                                                      # Depot slack time
    n::Int                                                                          # Depot total nodes visited
    q::Float64                                                                      # Depot total demand served
    l::Float64                                                                      # Depot total route length
    πᵒ::Float64                                                                     # Depot operational cost
    πᶠ::Float64                                                                     # Depot fixed cost
    φ::Int                                                                          # Depot operations mandate
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end
"""
    CustomerNode(iⁿ::Int, iʳ::Int, iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, qᶜ::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, iᵗ::Int, iʰ::Int, tᵃ::Float64, tᵈ::Float64, r::Route, τ::Float64, n::Int, q::Float64, l::Float64)

A `CustomerNode` is a source/sink point on the graph at `(x,y)` with index `iⁿ` 
and associated delivery/pickup node index `jⁿ`, on route `r` with route index `iʳ`, 
vehicle index `iᵛ`, and depot node index `iᵈ`, customer demand `qᶜ`, service time 
`τᶜ`, earliest service time `tᵉ`, latest service time `tˡ`, tail node index `iᵗ`, 
head node index `iʰ`, arrival time `tᵃ`, departure time `tᵈ`, slack time `τ`, 
position on the route `n`, vehicle load `q` and route length `l` on arrival at 
customer node.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int                                                                         # Customer node index
    jⁿ::Int                                                                         # Associated pickup/delivery node index
    iʳ::Int                                                                         # Route index
    iᵛ::Int                                                                         # Vehicle index
    iᵈ::Int                                                                         # Depot node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    qᶜ::Float64                                                                     # Customer demand
    τᶜ::Float64                                                                     # Customer service time
    tᵉ::Float64                                                                     # Customer node earliest service time
    tˡ::Float64                                                                     # Customer node latest service time
    iᵗ::Int                                                                         # Tail (predecessor) node index
    iʰ::Int                                                                         # Head (successor) node index
    tᵃ::Float64                                                                     # Customer node arrival time
    tᵈ::Float64                                                                     # Customer node departure time
    τ::Float64                                                                      # Customer slack time
    n::Int                                                                          # Customer position on the route
    q::Float64                                                                      # Vehicle load on arrival
    l::Float64                                                                      # Route length on arrival
    r::Route                                                                        # Route visiting the customer node
end



"""
    Solution(D::Vector{DepotNode}, C::OffsetVector{CustomerNode}, A::Dict{Tuple{Int,Int}, Arc})

A `Solution` is a graph with depot nodes `D`, customer nodes `C`, and arcs `A`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int,Int}, Arc}                                                    # Set of arcs
end