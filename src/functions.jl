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
Two routes are the equal if their indices (`iᵈ`, `iᵛ`, `iʳ`) match.
"""
Base.isequal(r₁::Route, r₂::Route) = isequal(r₁.iᵈ, r₂.iᵈ) && isequal(r₁.iᵛ, r₂.iᵛ) && isequal(r₁.iʳ, r₂.iʳ)
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
isopt(v::Vehicle) = !iszero(v.n)
"""
    isopt(d::DepotNode)
    
Returns `true` if depot node `d` is operational.
A `DepotNode` is defined operational if it serves at least one customer
unless it is mandated to be operational.
"""
isopt(d::DepotNode) = !iszero(d.n) || !iszero(d.φ)



"""
    isopen(c::CustomerNode)
    
Returns `true` if customer node `c` is open.
A `CustomerNode` is defined open if it is not being serviced.
"""
isopen(c::CustomerNode) = isequal(c.r, NullRoute)
"""
    isopen(d::DepotNode)

Returns `true` if depot node `d` is operational.
A `DepotNode` is defined operational if it serves at least one customer
unless it is mandated to be operational.
"""  
isopen(d::DepotNode) = !iszero(d.n) || !iszero(d.φ)



"""
    isclose(c::CustomerNode)

Returns `true` if customer node `c` is closed.
A `CustomerNode` is defined closed if it is being serviced.
"""  
isclose(c::CustomerNode) = !isequal(c.r, NullRoute)
"""
    isclose(d::DepotNode)

Returns `true` if depot node `d` is not operational.
A `DepotNode` is defined non-operational if it serves no customer
given it is not mandated to be operational.
"""  
isclose(d::DepotNode) = iszero(d.n) && iszero(d.φ)



"""
    isactive(r::Route)

A `Route` is defined active if it has not been initialized yet.
Usage in dynamic execution and simulation.
"""
isactive(r::Route) = isone(r.φ)
"""
    isactive(c::CustomerNode)

A `CustomerNode` is defined active if its route has not been initialized yet.
Usage in dynamic execution and simulation.
"""
isactive(c::CustomerNode) = isactive(c.r)



"""
    isdormant(r::Route)

A `Route` is defined dormant if it has been initialized.
Usage in dynamic execution and simulation.
"""
isdormant(r::Route) = iszero(r.φ)
"""
    isactive(c::CustomerNode)

A `CustomerNode` is defined dormant if its route has been initialized.
Usage in dynamic execution and simulation.
"""
isdormant(c::CustomerNode) = isdormant(c.r)



"""
    hasslack(d::DepotNode)
    
Returns `true` if depot node `d` has slack.
A `DepotNode` is defined to have slack if it has the capacity to serve an additional customer.
"""
hasslack(d::DepotNode) = d.q < d.qᵈ



"""
    Route(v::Vehicle, d::DepotNode)

Returns a non-operational `Route` traversed by vehicle `v` from depot node `d`.
"""
function Route(v::Vehicle, d::DepotNode)
    iʳ = lastindex(v.R) + 1
    iᵛ = v.iᵛ
    iᵈ = d.iⁿ
    x  = 0.
    y  = 0. 
    iˢ = iᵈ
    iᵉ = iᵈ
    θⁱ = isone(iʳ) ? 1.0 : v.R[iʳ-1].θᵉ
    θˢ = θⁱ
    θᵉ = θˢ
    tⁱ = isone(iʳ) ? d.tˢ : v.R[iʳ-1].tᵉ
    tˢ = tⁱ
    tᵉ = tⁱ
    τ  = Inf
    n  = 0 
    q  = 0.
    l  = 0.
    φ  = 1
    r  = Route(iʳ, iᵛ, iᵈ, x, y, iˢ, iᵉ, θⁱ, θˢ, θᵉ, tⁱ, tˢ, tᵉ, τ, n, q, l, φ)
    return r
end
"""
    NullRoute

A `NullRoute` is a fictitious out-of-service route.
"""           
const NullRoute = Route(0, 0, 0, 0., 0., 0, 0, 0., 0., 0., Inf, Inf, Inf, 0., 0, 0, Inf, 1)



"""
    Solution(D::Vector{DepotNode}, C::OffsetVector{CustomerNode}, A::Dict{Tuple{Int,Int}, Arc})

Returns `Solution` on graph `G = (D, C, A)`.
"""
function Solution(D::Vector{DepotNode}, C::OffsetVector{CustomerNode}, A::Dict{Tuple{Int,Int}, Arc})
    πᶠ = 0.
    πᵒ = 0.
    πᵖ = 0.
    for d ∈ D πᶠ += d.φ * d.πᶠ end
    for c ∈ C πᵖ += abs(c.qᶜ) end
    return Solution(D, C, A, πᶠ, πᵒ, πᵖ)
end



"""
    Vehicle(v::Vehicle, d::DepotNode)

Returns a non-operational `Vehicle` cloning vehicle `v` at depot node `d`.
"""
function Vehicle(v::Vehicle, d::DepotNode)
    iᵛ = lastindex(d.V) + 1
    jᵛ = v.jᵛ
    iᵈ = v.iᵈ
    qᵛ = v.qᵛ
    lᵛ = v.lᵛ
    sᵛ = v.sᵛ
    τᶠ = v.τᶠ
    τᵈ = v.τᵈ
    τᶜ = v.τᶜ
    τʷ = v.τʷ
    r̅  = v.r̅
    tˢ = d.tˢ
    tᵉ = d.tˢ
    τ  = Inf
    n  = 0
    q  = 0.
    l  = 0.
    πᵈ = v.πᵈ
    πᵗ = v.πᵗ
    πᶠ = v.πᶠ
    R  = Route[]
    v  = Vehicle(iᵛ, jᵛ, iᵈ, qᵛ, lᵛ, sᵛ, τᶠ, τᵈ, τᶜ, τʷ, r̅, tˢ, tᵉ, τ, n, q, l, πᵈ, πᵗ, πᶠ, R)
    return v
end



"""
    vectorize(s::Solution)

Returns `Solution` as a sequence of nodes in the order of visits.
"""
function vectorize(s::Solution)
    Z = [Int[] for _ ∈ s.D]
    for d ∈ s.D
        iⁿ = d.iⁿ
        if !isopt(d) continue end
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ] 
                push!(Z[iⁿ], d.iⁿ)
                c  = cˢ
                while true
                    push!(Z[iⁿ], c.iⁿ)
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
            end
        end
        push!(Z[iⁿ], d.iⁿ)
    end
    return Z
end
"""
    hash(s::Solution)

Returns hash on vectorized `Solution`.
"""
Base.hash(s::Solution) = hash(vectorize(s))



"""
    f(s::Solution)

Returns objective function evaluation for solution `s`
"""
f(s::Solution) = s.πᶠ + s.πᵒ + s.πᵖ * 10 ^ ceil(log10(s.πᶠ + s.πᵒ))



"""
    isfeasible(s::Solution)

Returns `true` if customer service and time-window constraints;
vehicle range, capacity, and working-hours constraints; and
depot capacity constraints are not violated.
"""
function isfeasible(s::Solution)
    if any(isopen, s.C) return false end                                    # Service constraint
    for d ∈ s.D
        if !isopt(d) continue end
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ] 
                c  = cˢ
                while true
                    if c.tᵃ > c.tˡ return false end                         # Time-window constraint
                    if isdelivery(c)
                        nᵖ  = s.C[c.jⁿ]
                        nᵈ  = s.C[c.iⁿ]
                        qᵒ = c.q
                    else
                        nᵖ  = s.C[c.iⁿ]
                        nᵈ  = s.C[c.jⁿ]
                        qᵒ = c.q - c.qᶜ
                    end
                    if nᵖ.tᵃ > nᵈ.tᵃ return false end                       # Service constraint (order of service)
                    if !isequal(nᵖ.r, nᵈ.r) return false end                # Service constraint (order of service)
                    if c.l > v.lᵛ return false end                          # Vehicle range constraint
                    if qᵒ > v.qᵛ return false end                           # Vehicle capacity constraint
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
                if r.l > v.lᵛ return false end                              # Vehicle range constraint
                if r.q > v.qᵛ return false end                              # Vehicle capacity constraint
            end
            if d.tˢ > v.tˢ return false end                                 # Working-hours constraint (start time)
            if v.tᵉ > d.tᵉ return false end                                 # Working-hours constraint (end time)
            if v.tᵉ - v.tˢ > v.τʷ return false end                          # Working-hours constraint (duration)
            if length(v.R) > v.r̅ return false end                           # Number of routes constraint
        end
        if d.q > d.qᵈ return false end                                      # Depot capacity constraint
    end
    return true
end



"""
    relatedness(c₁::CustomerNode, c₂::CustomerNode, s::Solution)

Returns a measure of similarity between customer nodes `c₁` and `c₂` in solution `s`.
"""
function relatedness(c₁::CustomerNode, c₂::CustomerNode, s::Solution)
    ϵ  = 1e-5
    r₁ = c₁.r
    r₂ = c₂.r
    φ  = (1 + isequal(r₁, r₂)) / 2
    q  = abs(c₁.qᶜ - c₂.qᶜ)
    l  = s.A[(c₁.iⁿ,c₂.iⁿ)].l
    t  = abs(c₁.tᵉ - c₂.tᵉ) + abs(c₁.tˡ - c₂.tˡ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(r₁::Route, r₂::Route, s::Solution)

Returns a measure of similarity between routes `r₁` and `r₂` in solution `s`.
"""
function relatedness(r₁::Route, r₂::Route, s::Solution)
    ϵ  = 1e-5
    d₁ = s.D[r₁.iᵈ]
    d₂ = s.D[r₂.iᵈ]
    v₁ = d₁.V[r₁.iᵛ]
    v₂ = d₂.V[r₂.iᵛ]
    φ  = 1
    q  = abs(v₁.qᵛ - v₂.qᵛ)
    l  = sqrt((r₁.x - r₂.x)^2 + (r₁.y - r₂.y)^2)
    t  = abs(r₁.tˢ - r₂.tˢ) + abs(r₁.tᵉ - r₂.tᵉ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(v₁::Vehicle, v₂::Vehicle, s::Solution)

Returns a measure of similarity between vehicles `v₁` and `v₂` in solution `s`.
"""
function relatedness(v₁::Vehicle, v₂::Vehicle, s::Solution)
    ϵ  = 1e-5
    x₁ = 0.
    y₁ = 0.
    for r ∈ v₁.R 
        x₁ += r.n * r.x / v₁.n
        y₁ += r.n * r.y / v₁.n 
    end
    x₂ = 0.
    y₂ = 0.
    for r ∈ v₂.R 
        x₂ += r.n * r.x / v₂.n
        y₂ += r.n * r.y / v₂.n
    end
    φ  = 1
    q  = abs(v₁.qᵛ - v₂.qᵛ)
    l  = sqrt((x₁ - x₂)^2 + (y₁ - y₂)^2)
    t  = abs(v₁.tˢ - v₂.tˢ) + abs(v₁.tᵉ - v₂.tᵉ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(d₁::DepotNode, d₂::DepotNode, s::Solution)

Returns a measure of similarity between depot nodes `d₁` and `d₂` in solution `s`.
"""
function relatedness(d₁::DepotNode, d₂::DepotNode, s::Solution)
    ϵ  = 1e-5
    φ  = 1
    q  = abs(d₁.qᵈ - d₂.qᵈ)
    l  = s.A[(d₁.iⁿ, d₂.iⁿ)].l
    t  = abs(d₁.tˢ - d₂.tˢ) + abs(d₁.tᵉ - d₂.tᵉ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(c::CustomerNode, d::DepotNode, s::Solution)

Returns a measure of similarity between customer nodes `c` and depot node `d` in solution `s`.
"""
function relatedness(c::CustomerNode, d::DepotNode, s::Solution)
    ϵ  = 1e-5
    φ  = 1
    q  = 0.
    l  = s.A[(c.iⁿ,d.iⁿ)].l
    t  = 0.
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(d::DepotNode, c::CustomerNode, s::Solution)

Returns a measure of similarity between depot node `d` and customer nodes `c` in solution `s`.
"""
relatedness(d::DepotNode, c::CustomerNode, s::Solution) = relatedness(c, d, s)