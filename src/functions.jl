"""
    isdepot(n::Node)

Returns `true` if node `n` is a `DepotNode`.
"""
isdepot(n::Node) = isequal(typeof(n), DepotNode)
"""
    iscustomer(n::Node)
    
Returns `true` if node `n` is a `CustomerNode`.
"""
iscustomer(n::Node) = isequal(typeof(n), CustomerNode)
"""
    isfuelstation(n::Node)
    
Returns `true` if node `n` is a `CustomerNode`.
"""
isfuelstation(n::Node) = isequal(typeof(n), FuelStationNode)


"""
    ispickup(c::CustomerNode)

Returns `true` if node `c` is a pickup node.
"""
ispickup(c::CustomerNode) = c.qᶜ < 0.
"""
    isdelivery(c::CustomerNode)

Returns `true` if node `c` is a delivery node.
"""
isdelivery(c::CustomerNode) = c.qᶜ > 0.



"""
    isequal(n₁::Node, n₂::Node)

Return `true` if node `n₁` equals node `n₂`.
Two nodes are equal if their indices (`iⁿ`) match.
"""
Base.isequal(n₁::Node, n₂::Node) = isequal(n₁.iⁿ, n₂.iⁿ)
"""
    isequal(r₁::Route, r₂::Route)

Return `true` if route `r₁` equals route `r₂`.
Two routes are the equal if their indices (`iᵈ`, `iᵛ`) match.
"""
Base.isequal(r₁::Route, r₂::Route) = isequal(r₁.iᵈ, r₂.iᵈ) && isequal(r₁.iᵛ, r₂.iᵛ)
"""
    isequal(v₁::Vehicle, v₂::Vehicle)

Return `true` if vehicle `v₁` equals vehicle `v₂`.
Two vehicles are equal if their indices (`iᵈ`, `iᵛ`) match.
"""
Base.isequal(v₁::Vehicle, v₂::Vehicle) = isequal(v₁.iᵈ, v₂.iᵈ) && isequal(v₁.iᵛ, v₂.iᵛ)



"""
    isopt(r::Route)

Returns `true` if route `r` is operational.
A `Route` is defined operational if it serves at least one customer.
"""
isopt(r::Route) = !iszero(r.n)
"""
    isopt(v::Vehicle)

Returns `true` if vehicle `v` is operational.
A `Vehicle` is defined operational if it serves at least one customer.
"""
isopt(v::Vehicle) = !iszero(v.r.n)
"""
    isopt(d::DepotNode)
    
Returns `true` if depot node `d` is operational.
A `DepotNode` is defined operational if it serves at least one customer.
"""
isopt(d::DepotNode) = !iszero(d.n)
"""
    isopt(f::FuelStationNode)
    
Returns `true` if fuel station node `f` is operational.
A `FuelStationNode` is defined operational if it refuels at least one vehicle.
"""
isopt(f::FuelStationNode) = !iszero(f.ω)



"""
    isopen(c::CustomerNode)
    
Returns `true` if customer node `c` is open.
A `CustomerNode` is defined open if it is not being serviced.
"""
isopen(c::CustomerNode) = isequal(c.r, NullRoute)
"""
    isopen(d::DepotNode)

Returns `true` if depot node `d` is open.
A `DepotNode` is defined open if it serves at least one customer.
"""  
isopen(d::DepotNode) = !iszero(d.n)
"""
    isopen(f::FuelStationNode)
    
Returns `true` if fuel station node `f` is open.
A `FuelStationNode` is defined open if it refuels at least one vehicle.
"""
isopen(f::FuelStationNode) = !iszero(f.n)



"""
    isclose(c::CustomerNode)

Returns `true` if customer node `c` is closed.
A `CustomerNode` is defined closed if it is being serviced.
"""  
isclose(c::CustomerNode) = !isequal(c.r, NullRoute)
"""
    isclose(d::DepotNode)

Returns `true` if depot node `d` is closed.
A `DepotNode` is defined closed if it serves no customer.
"""  
isclose(d::DepotNode) = iszero(d.n)
"""
    isclose(f::FuelStationNode)
    
Returns `true` if fuel station node `f` is open.
A `FuelStationNode` is defined closed if it refuels no vehicle.
"""
isclose(f::FuelStationNode) = iszero(f.n)



"""
    relatedness(m::Symbol, c₁::CustomerNode, c₂::CustomerNode, s::Solution)

Returns a measure of similarity between customer nodes `c₁` and `c₂` based on metric `m` in solution `s`.
"""
function relatedness(m::Symbol, c₁::CustomerNode, c₂::CustomerNode, s::Solution;)
    ϵ   = 1e-5
    cᵖ₁ = isdelivery(c₁) ? s.C[c₁.jⁿ] : s.C[c₁.iⁿ] 
    cᵈ₁ = isdelivery(c₁) ? s.C[c₁.iⁿ] : s.C[c₁.jⁿ]
    cᵖ₂ = isdelivery(c₂) ? s.C[c₂.jⁿ] : s.C[c₂.iⁿ] 
    cᵈ₂ = isdelivery(c₂) ? s.C[c₂.iⁿ] : s.C[c₂.jⁿ]
    φ   = 1
    q   = isequal(m, :q) * (abs(c₁.qᶜ - c₂.qᶜ))
    l   = isequal(m, :l) * (s.A[(cᵖ₁.iⁿ,cᵖ₂.iⁿ)].l + s.A[(cᵈ₁.iⁿ,cᵈ₂.iⁿ)].l)
    t   = isequal(m, :t) * (abs(cᵖ₁.tᵉ - cᵖ₂.tᵉ) + abs(cᵖ₁.tˡ - cᵖ₂.tˡ) + abs(cᵈ₁.tᵉ - cᵈ₂.tᵉ) + abs(cᵈ₁.tˡ - cᵈ₂.tˡ))
    z   = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(m::Symbol, r₁::Route, r₂::Route, s::Solution)

Returns a measure of similarity between routes `r₁` and `r₂` based on metric `m` in solution `s`.
"""
function relatedness(m::Symbol, r₁::Route, r₂::Route, s::Solution)
    ϵ = 1e-5
    φ = 1
    q = isequal(m, :q) * (0.)
    l = isequal(m, :l) * (sqrt((r₁.x - r₂.x)^2 + (r₁.y - r₂.y)^2))
    t = isequal(m, :t) * (abs(r₁.tˢ - r₂.tˢ) + abs(r₁.tᵉ - r₂.tᵉ))
    z = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(m::Symbol, v₁::Vehicle, v₂::Vehicle, s::Solution)

Returns a measure of similarity between vehicles `v₁` and `v₂` based on metric `m` in solution `s`.
"""
function relatedness(m::Symbol, v₁::Vehicle, v₂::Vehicle, s::Solution)
    ϵ  = 1e-5
    r₁ = v₁.r
    r₂ = v₂.r
    x₁ = r₁.x
    x₂ = r₂.x
    y₁ = r₁.y
    y₂ = r₂.y
    φ = 1
    q = isequal(m, :q) * (0.)
    l = isequal(m, :l) * (sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2))
    t = isequal(m, :t) * (abs(r₁.tˢ - r₂.tˢ) + abs(r₁.tᵉ - r₂.tᵉ))
    z = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(m::Symbol, d₁::DepotNode, d₂::DepotNode, s::Solution)

Returns a measure of similarity between depot nodes `d₁` and `d₂` based on metric `m` in solution `s`.
"""
function relatedness(m::Symbol, d₁::DepotNode, d₂::DepotNode, s::Solution)
    ϵ = 1e-5
    φ = 1
    q = isequal(m, :q) * (0.)
    l = isequal(m, :l) * (s.A[(d₁.iⁿ, d₂.iⁿ)].l)
    t = isequal(m, :t) * (abs(d₁.tˢ - d₂.tˢ) + abs(d₁.tᵉ - d₂.tᵉ))
    z = φ/(q + l + t + ϵ)
    return z
end



"""
    NullRoute

A `NullRoute` is a fictitious out-of-service route.
"""           
const NullRoute = Route(0, 0, 0., 0., 0, 0, Inf, Inf, 1., 0., Inf, Inf, 0, Inf)



"""
    vectorize(s::Solution)

Returns `Solution` as a sequence of nodes in the order of visits for every depot, vehicle, and route.
"""
function vectorize(s::Solution)
    Z = [[Int[] for v ∈ d.V] for d ∈ s.D]
    for d ∈ s.D
        if !isopt(d) continue end
        iⁿ = d.iⁿ
        for v ∈ d.V
            if !isopt(v) continue end
            iᵛ = v.iᵛ
            r  = v.r
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ] 
            f  = cˢ.F[v.jᵛ]
            push!(Z[iⁿ][iᵛ], d.iⁿ)
            if r.ω > 0. push!(Z[iⁿ][iᵛ], f.iⁿ) end
            c  = cˢ
            while true
                f = c.F[v.jᵛ]
                push!(Z[iⁿ][iᵛ], c.iⁿ)
                if c.ω > 0. push!(Z[iⁿ][iᵛ], f.iⁿ) end
                if isequal(c, cᵉ) break end
                c = s.C[c.iʰ]
            end
            push!(Z[iⁿ][iᵛ], d.iⁿ)
        end
    end
    return Z
end
"""
    hash(s::Solution)

Returns hash on vectorized `Solution`.
"""
Base.hash(s::Solution) = hash(vectorize(s))



"""
    isfeasible(s::Solution)

Returns `true` if customer service and time-window constraints;
vehicle capacity, range, and working-hours constraints;
are not violated.
"""
function isfeasible(s::Solution)
    if any(isopen, s.C) return false end                                    # Service constraint
    for d ∈ s.D
        if !isopt(d) continue end
        for v ∈ d.V
            if !isopt(v) continue end
            r  = v.r
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            c  = cˢ
            while true
                if c.tᵃ > c.tˡ return false end                             # Time-window constraint
                cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
                cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
                if !isequal(cᵖ.r, cᵈ.r) return false end                    # Service constraint (order of service)
                if cᵖ.tᵃ > cᵈ.tᵃ return false end                           # Service constraint (order of service)
                if c.q > v.qᵛ return false end                              # Vehicle capacity constraint
                if c.θ̲ > c.θ return false end                               # Vehicle range constraint
                if isequal(c, cᵉ) break end
                c  = s.C[c.iʰ]
            end
            if r.θ̲ > r.θ return false end                                   # Vehicle range constraint
            if d.tˢ > r.tˢ return false end                                 # Working-hours constraint (start time)
            if r.tᵉ > d.tᵉ return false end                                 # Working-hours constraint (end time)
            if r.tᵉ - r.tˢ > v.τʷ return false end                          # Working-hours constraint (duration)
        end
    end
    return true
end



"""
    f(s::Solution)

Returns objective function evaluation for solution `s`
"""
f(s::Solution) = s.πᶠ + s.πᵒ + s.πᵖ * 10 ^ ceil(log10(s.πᶠ + s.πᵒ))