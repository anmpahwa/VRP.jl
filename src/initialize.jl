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
    savings([rng::AbstractRNG], instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial `Solution` created by merging routes that render the most 
savings until no merger can render further savings. 

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
function savings(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1: Initialize
    G = VRP.build(instance; dir=dir)
    s = VRP.Solution(G...)
    D = s.D
    C = s.C
    VRP.preinitialize!(s)
    # Step 2: Initialize solution with routes to every node from the depot node (first node)
    for c ∈ C
        if VRP.ispickup(c) continue end
        d = VRP.sample(rng, D)
        v = d.V[lastindex(d.V)]
        r = v.R[lastindex(v.R)]
        cᵖ = C[c.jⁿ]
        cᵈ = C[c.iⁿ]
        VRP.insertnode!(cᵖ, d, d, r, s)
        VRP.insertnode!(cᵈ, cᵖ, d, r, s)
        v = VRP.Vehicle(v, d)
        r = VRP.Route(v, d)
        push!(v.R, r)
        push!(d.V, v)
    end
    # Step 3: Merge routes iteratively until single route traversing all nodes remains
    R = [r for d ∈ s.D for v ∈ d.V for r ∈ v.R]
    K = eachindex(R)
    X = fill(Inf, (K,K))            # X[i,j]: Savings from merging route with tail node N[i] into route with tail node N[j]
    while true
        z = f(s)
        for (k₁,r₁) ∈ pairs(R)
            for (k₂,r₂) ∈ pairs(R)
                if isequal(r₁, r₂) continue end
                if !VRP.isopt(r₁) || !VRP.isopt(r₂) continue end
                cˢ  = C[r₁.iˢ]
                cᵉ  = C[r₁.iᵉ] 
                c   = cˢ
                nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                nᵗ₂ = C[r₂.iᵉ]
                nʰ₂ = D[r₂.iᵈ]
                while true
                    VRP.removenode!(c, nᵗ₁, nʰ₁, r₁, s)
                    VRP.insertnode!(c, nᵗ₂, nʰ₂, r₂, s)
                    if isequal(c, cᵉ) break end
                    c   = C[r₁.iˢ]
                    nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                    nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                    nᵗ₂ = C[r₂.iᵉ]
                    nʰ₂ = D[r₂.iᵈ]
                end
                z′ = f(s)
                Δ  = z′ - z
                X[k₁,k₂] = Δ
                c   = cᵉ
                nᵗ₁ = D[r₁.iᵈ]
                nʰ₁ = D[r₁.iᵈ]
                nᵗ₂ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                nʰ₂ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                while true
                    VRP.removenode!(c, nᵗ₂, nʰ₂, r₂, s)
                    VRP.insertnode!(c, nᵗ₁, nʰ₁, r₁, s)
                    if isequal(c, cˢ) break end
                    c   = C[r₂.iᵉ]
                    nᵗ₁ = D[r₁.iᵈ]
                    nʰ₁ = C[r₁.iˢ]
                    nᵗ₂ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                    nʰ₂ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                end
            end
        end
        k₁,k₂ = Tuple(argmin(X))
        Δ = X[k₁,k₂]
        if Δ > 0 break end 
        r₁  = R[k₁]
        r₂  = R[k₂]
        cˢ  = C[r₁.iˢ]
        cᵉ  = C[r₁.iᵉ] 
        c   = cˢ
        nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
        nᵗ₂ = C[r₂.iᵉ]
        nʰ₂ = D[r₂.iᵈ]
        while true
            VRP.removenode!(c, nᵗ₁, nʰ₁, r₁, s)
            VRP.insertnode!(c, nᵗ₂, nʰ₂, r₂, s)
            if isequal(c, cᵉ) break end
            c   = C[r₁.iˢ]
            nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
            nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
            nᵗ₂ = C[r₂.iᵉ]
            nʰ₂ = D[r₂.iᵈ]
        end
        X .= Inf
    end
    return s
end



"""
    initialize([rng::AbstractRNG], instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

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
initialize(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances")) = savings(rng, instance; dir=dir)
initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances")) = initialize(Random.GLOBAL_RNG, instance; dir=dir)