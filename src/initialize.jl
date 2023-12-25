"""
    build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
Returns a tuple of depot nodes, customer nodes, and arcs for the `instance`.

Note, `dir` locates the the folder containing instance files as sub-folders,
as follows,

    <dir>
    |-<instance>
        |-arcs.csv
        |-depot_nodes.csv
        |-customer_nodes.csv
        |-vehicles.csv
"""
function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Depot nodes
    df = DataFrame(CSV.File(joinpath(dir, "$instance/depot_nodes.csv")))
    D  = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        x  = df[k,2]
        y  = df[k,3]
        qᵈ = df[k,4]
        tˢ = df[k,5]
        tᵉ = df[k,6]
        τ  = Inf
        n  = 0
        q  = 0.
        l  = 0.
        πᵒ = df[k,7]
        πᶠ = df[k,8]
        φ  = df[k,9]
        d  = DepotNode(iⁿ, x, y, qᵈ, tˢ, tᵉ, τ, n, q, l, πᵒ, πᶠ, φ, Vehicle[])
        D[iⁿ] = d
    end
    # Customer nodes
    df = DataFrame(CSV.File(joinpath(dir, "$instance/customer_nodes.csv")))
    I  = (df[1,1]:df[nrow(df),1])::UnitRange{Int}
    C  = OffsetVector{CustomerNode}(undef, I)
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        jⁿ = df[k,2]
        iʳ = 0
        iᵛ = 0
        iᵈ = 0
        x  = df[k,3]
        y  = df[k,4]
        qᶜ = df[k,5]
        τᶜ = df[k,6]
        tᵉ = df[k,7]
        tˡ = df[k,8]
        iᵗ = 0
        iʰ = 0
        tᵃ = qᶜ > 0. ? tˡ : tᵉ
        tᵈ = tᵃ + τᶜ
        τ  = Inf
        n  = 0
        q  = 0.
        l  = 0.
        c  = CustomerNode(iⁿ, jⁿ, iʳ, iᵛ, iᵈ, x, y, qᶜ, τᶜ, tᵉ, tˡ, iᵗ, iʰ, tᵃ, tᵈ, τ, n, q, l, NullRoute)
        C[iⁿ] = c
    end
    # Arcs
    df = DataFrame(CSV.File(joinpath(dir, "$instance/arcs.csv"), header=false))
    A  = Dict{Tuple{Int,Int},Arc}()
    n  = lastindex(C)
    for iᵗ ∈ 1:n
        for iʰ ∈ 1:n
            l = df[iᵗ,iʰ] 
            a = Arc(iᵗ, iʰ, l)
            A[(iᵗ,iʰ)] = a
        end
    end
    # Vehicles
    df = DataFrame(CSV.File(joinpath(dir, "$instance/vehicles.csv")))
    for k ∈ 1:nrow(df)
        d  = D[df[k,3]]
        iᵛ = df[k,1]
        jᵛ = df[k,2]
        iᵈ = df[k,3]
        qᵛ = df[k,4]
        lᵛ = df[k,5]
        sᵛ = df[k,6]
        τᶠ = df[k,7]
        τᵈ = df[k,8]
        τᶜ = df[k,9]
        τʷ = df[k,10]
        r̅  = df[k,11]
        tˢ = d.tˢ
        tᵉ = d.tˢ
        τ  = Inf
        n  = 0
        q  = 0.
        l  = 0.
        πᵈ = df[k,12]
        πᵗ = df[k,13]
        πᶠ = df[k,14]
        v  = Vehicle(iᵛ, jᵛ, iᵈ, qᵛ, lᵛ, sᵛ, τᶠ, τᵈ, τᶜ, τʷ, r̅, tˢ, tᵉ, τ, n, q, l, πᵈ, πᵗ, πᶠ, Route[])
        push!(d.V, v)
    end
    V  = [v for d ∈ D for v ∈ d.V]
    φᵈ = !iszero(getproperty.(D, :tˢ)) || !iszero(getproperty.(D, :tᵉ))
    φᶜ = !iszero(getproperty.(C, :tᵉ)) || !iszero(getproperty.(C, :tˡ)) || !iszero(getproperty.(C, :jⁿ))
    φᵛ = !iszero(getproperty.(V, :τʷ)) || !iszero(getproperty.(V, :πᵗ))
    global φᵉ = (φᵈ || φᶜ || φᵛ)::Bool
    G  = (D, C, A)
    return G
end



"""
    initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial VRP `Solution` for the `instance`. 

Note, `dir` locates the the folder containing instance files as sub-folders, 
as follows,
    
    <dir>
    |-<instance>
        |-arcs.csv
        |-depot_nodes.csv
        |-customer_nodes.csv
        |-vehicles.csv
        
Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function initialize(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    G = build(instance; dir=dir)
    s = Solution(G...)
    preinitialize!(s)
    for c ∈ s.C
        if ispickup(c) continue end
        d = sample(rng, s.D)
        v = d.V[lastindex(d.V)]
        r = v.R[lastindex(v.R)]
        if c.jⁿ ≤ lastindex(s.D)
            insertnode!(c, d, d, r, s)
        else
            cᵖ = s.C[c.jⁿ]
            cᵈ = s.C[c.iⁿ]
            insertnode!(cᵖ, d, d, r, s)
            insertnode!(cᵈ, cᵖ, d, r, s)
        end
        v = Vehicle(v, d)
        r = Route(v, d)
        push!(v.R, r)
        push!(d.V, v)
    end
    return s
end