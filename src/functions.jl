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
ispickup(c::CustomerNode) = c.qᶜ < 0
"""
    isdelivery(c::CustomerNode)

Returns `true` if node `c` is a delivery node.
"""
isdelivery(c::CustomerNode) = c.qᶜ > 0



"""
    isequal(p::Route, q::Route)

Return `true` if route `p` equals route `q`.
Two routes are the equal if their indices (`iᵈ`, `iᵛ`, `iʳ`) match.
"""
Base.isequal(p::Route, q::Route) = isequal(p.iᵈ, q.iᵈ) && isequal(p.iᵛ, q.iᵛ) && isequal(p.iʳ, q.iʳ)
"""
    isequal(p::Vehicle, q::Vehicle)

Return `true` if vehicle `p` equals vehicle `q`.
Two vehicles are equal if their indices (`iᵈ`, `iᵛ`) match.
"""
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.iᵈ, q.iᵈ) && isequal(p.iᵛ, q.iᵛ)
"""
    isequal(p::Node, q::Node)

Return `true` if node `p` equals node `q`.
Two nodes are equal if their indices (`iⁿ`) match.
"""
Base.isequal(p::Node, q::Node) = isequal(p.iⁿ, q.iⁿ)



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
    iʳ = length(v.R) + 1
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
    Vehicle(v::Vehicle, d::DepotNode)

Returns a non-operational `Vehicle` cloning vehicle `v` at depot node `d`.
"""
function Vehicle(v::Vehicle, d::DepotNode)
    iᵛ = length(d.V) + 1
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
    D = s.D
    C = s.C
    Z = [Int[] for _ ∈ D]
    for d ∈ D
        iⁿ = d.iⁿ
        if !isopt(d) continue end
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
                cˢ = C[r.iˢ]
                cᵉ = C[r.iᵉ] 
                push!(Z[iⁿ], d.iⁿ)
                cᵒ = cˢ
                while true
                    push!(Z[iⁿ], cᵒ.iⁿ)
                    if isequal(cᵒ, cᵉ) break end
                    cᵒ = C[cᵒ.iʰ]
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
    f(s::Solution; fixed=true, operational=true, penalty=true)

Returns objective function evaluation for solution `s`. Includes `fixed`, 
`operational`, and `penalty` cost for constraint violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, penalty=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    φᶠ, φᵒ, φᵖ = fixed, operational, penalty
    for c ∈ s.C 
        if isopen(c) πᵖ += abs(c.qᶜ)                                # Service constraint (node visit)
        else
            d = s.D[c.iᵈ]
            v = d.V[c.iᵛ]
            if isdelivery(c)
                nᵖ = c.jⁿ ≤ length(s.D) ? s.D[c.jⁿ] : s.C[c.jⁿ]
                rᵖ = isdepot(nᵖ) ? c.r : nᵖ.r
                tᵖ = isdepot(nᵖ) ? rᵖ.tˢ : nᵖ.tᵃ
                nᵈ = c
                rᵈ = nᵈ.r
                tᵈ = nᵈ.tᵃ
                qᵒ = c.q
            end
            if ispickup(c)
                nᵖ = c
                rᵖ = nᵖ.r
                tᵖ = nᵖ.tᵃ
                nᵈ = s.C[c.jⁿ]
                rᵈ = nᵈ.r
                tᵈ = nᵈ.tᵃ
                qᵒ = c.q - c.qᶜ
            end
            πᵖ += (!isequal(rᵖ, rᵈ) || (tᵖ > tᵈ)) * nᵈ.qᶜ           # Service constraint (order of service)
            πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)                     # Time-window constraint
            πᵖ += (c.l > v.lᵛ) * (c.l - v.lᵛ)                       # Vehicle range constraint
            πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)                         # Vehicle capacity constraint
        end
    end
    for d ∈ s.D
        if !isopt(d) continue end
        πᶠ += d.πᶠ
        πᵒ += d.q * d.πᵒ
        for v ∈ d.V
            if !isopt(v) continue end
            πᶠ += v.πᶠ
            πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
            for r ∈ v.R 
                if !isopt(r) continue end
                πᶠ += 0.
                πᵒ += r.l * v.πᵈ
                πᵖ += (r.l > v.lᵛ) * (r.l - v.lᵛ)                   # Vehicle range constraint
                πᵖ += (r.q > v.qᵛ) * (r.q - v.qᵛ)                   # Vehicle capacity constraint
            end
            πᵖ += (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)                     # Working-hours constraint (start time)
            πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)                     # Working-hours constraint (end time)
            πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * (v.tᵉ - v.tˢ - v.τʷ)       # Working-hours constraint (duration)
            πᵖ += (length(v.R) > v.r̅) * (length(v.R) - v.r̅)         # Number of routes constraint
        end
        πᵖ += (d.q > d.qᵈ) * (d.q - d.qᵈ)                           # Depot capacity constraint
    end
    πᵖ *= 10 ^ ceil(log10(πᶠ + πᵒ))
    z = φᶠ * πᶠ + φᵒ * πᵒ + φᵖ * πᵖ
    return z
end




"""
    isfeasible(s::Solution)

Returns `true` if customer service and time-window constraints;
vehicle range, capacity, and working-hours constraints; and
depot capacity constraints are not violated.
"""
function isfeasible(s::Solution)
    for c ∈ s.C
        if isopen(c) return false                                   # Service constraint (node visit)
        else
            d = s.D[c.iᵈ]
            v = d.V[c.iᵛ]
            if isdelivery(c)
                nᵖ = c.jⁿ ≤ length(s.D) ? s.D[c.jⁿ] : s.C[c.jⁿ]
                rᵖ = isdepot(nᵖ) ? c.r : nᵖ.r
                tᵖ = isdepot(nᵖ) ? rᵖ.tˢ : nᵖ.tᵃ
                nᵈ = c
                rᵈ = nᵈ.r
                tᵈ = nᵈ.tᵃ
                qᵒ = c.q
            end
            if ispickup(c)
                nᵖ = c
                rᵖ = nᵖ.r
                tᵖ = nᵖ.tᵃ
                nᵈ = s.C[c.jⁿ]
                rᵈ = nᵈ.r
                tᵈ = nᵈ.tᵃ
                qᵒ = c.q - c.qᶜ
            end
            if (!isequal(rᵖ, rᵈ) || (tᵖ > tᵈ)) return false end     # Service constraint (order of service)
            if c.tᵃ > c.tˡ return false end                         # Time-window constraint
            if c.l > v.lᵛ return false end                          # Vehicle range constraint
            if qᵒ > v.qᵛ return false end                           # Vehicle capacity constraint
        end
    end
    for d ∈ s.D
        for v ∈ d.V
            for r ∈ v.R
                if r.l > v.lᵛ return false end                      # Vehicle range constraint
                if r.q > v.qᵛ return false end                      # Vehicle capacity constraint
            end
            if d.tˢ > v.tˢ return false end                         # Working-hours constraint (start time)
            if v.tᵉ > d.tᵉ return false end                         # Working-hours constraint (end time)
            if v.tᵉ - v.tˢ > v.τʷ return false end                  # Working-hours constraint (duration)
            if length(v.R) > v.r̅ return false end                   # Number of routes constraint
        end
        if d.q > d.qᵈ return false end                              # Depot capacity constraint
    end
    return true
end



"""
    relatedness(c¹::CustomerNode, c²::CustomerNode, s::Solution)

Returns a measure of similarity between customer nodes `c¹` and `c²` in solution `s`.
"""
function relatedness(c¹::CustomerNode, c²::CustomerNode, s::Solution)
    ϵ  = 1e-5
    r¹ = c¹.r
    r² = c².r
    φ  = (1 + isequal(r¹, r²)) / 2
    q  = abs(c¹.qᶜ - c².qᶜ)
    l  = s.A[(c¹.iⁿ,c².iⁿ)].l
    t  = abs(c¹.tᵉ - c².tᵉ) + abs(c¹.tˡ - c².tˡ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(r¹::Route, r²::Route, s::Solution)

Returns a measure of similarity between routes `r¹` and `r²` in solution `s`.
"""
function relatedness(r¹::Route, r²::Route, s::Solution)
    ϵ  = 1e-5
    d¹ = s.D[r¹.iᵈ]
    d² = s.D[r².iᵈ]
    v¹ = d¹.V[r¹.iᵛ]
    v² = d².V[r².iᵛ]
    φ  = 1
    q  = abs(v¹.qᵛ - v².qᵛ)
    l  = sqrt((r¹.x - r².x)^2 + (r¹.y - r².y)^2)
    t  = abs(r¹.tˢ - r².tˢ) + abs(r¹.tᵉ - r².tᵉ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(v¹::Vehicle, v²::Vehicle, s::Solution)

Returns a measure of similarity between vehicles `v¹` and `v²` in solution `s`.
"""
function relatedness(v¹::Vehicle, v²::Vehicle, s::Solution)
    ϵ  = 1e-5
    x¹ = 0.
    y¹ = 0.
    for r ∈ v¹.R 
        x¹ += r.n * r.x / v¹.n
        y¹ += r.n * r.y / v¹.n 
    end
    x² = 0.
    y² = 0.
    for r ∈ v².R 
        x² += r.n * r.x / v².n
        y² += r.n * r.y / v².n
    end
    φ  = 1
    q  = abs(v¹.qᵛ - v².qᵛ)
    l  = sqrt((x¹ - x²)^2 + (y¹ - y²)^2)
    t  = abs(v¹.tˢ - v².tˢ) + abs(v¹.tᵉ - v².tᵉ)
    z  = φ/(q + l + t + ϵ)
    return z
end
"""
    relatedness(d¹::DepotNode, d²::DepotNode, s::Solution)

Returns a measure of similarity between depot nodes `d¹` and `d²` in solution `s`.
"""
function relatedness(d¹::DepotNode, d²::DepotNode, s::Solution)
    ϵ  = 1e-5
    φ  = 1
    q  = abs(d¹.qᵈ - d².qᵈ)
    l  = s.A[(d¹.iⁿ,d².iⁿ)].l
    t  = abs(d¹.tˢ - d².tˢ) + abs(d¹.tᵉ - d².tᵉ)
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