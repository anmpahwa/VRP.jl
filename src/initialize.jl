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
        n  = 0
        q  = 0.
        l  = 0.
        πᵒ = df[k,7]
        πᶠ = df[k,8]
        d  = DepotNode(iⁿ, x, y, qᵈ, tˢ, tᵉ, n, q, l, πᵒ, πᶠ, Vehicle[])
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
        n  = 0
        q  = 0.
        l  = 0.
        c  = CustomerNode(iⁿ, jⁿ, iʳ, iᵛ, iᵈ, x, y, qᶜ, τᶜ, tᵉ, tˡ, iᵗ, iʰ, tᵃ, tᵈ, n, q, l, NullRoute)
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
        n  = 0
        q  = 0.
        l  = 0.
        πᵈ = df[k,12]
        πᵗ = df[k,13]
        πᶠ = df[k,14]
        v  = Vehicle(iᵛ, jᵛ, iᵈ, qᵛ, lᵛ, sᵛ, τᶠ, τᵈ, τᶜ, τʷ, r̅, tˢ, tᵉ, n, q, l, πᵈ, πᵗ, πᶠ, Route[])
        push!(d.V, v)
    end
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
    G = build(instance; dir=dir)
    s = Solution(G...)
    D = s.D
    C = s.C
    preinitialize!(s)
    # Step 2: Initialize solution with routes to every pickup and delivery node from the depot node
    for c ∈ C
        if ispickup(c) continue end
        d = sample(rng, D)
        v = d.V[lastindex(d.V)]
        r = v.R[lastindex(v.R)]
        cᵖ = C[c.jⁿ]
        cᵈ = C[c.iⁿ]
        insertnode!(cᵖ, d, d, r, s)
        insertnode!(cᵈ, cᵖ, d, r, s)
        v = Vehicle(v, d)
        r = Route(v, d)
        push!(v.R, r)
        push!(d.V, v)
    end
    # Step 3: Merge routes iteratively until no merger can reduce objective funciton value
    R = [r for d ∈ s.D for v ∈ d.V for r ∈ v.R]
    K = eachindex(R)
    X = fill(Inf, (K,K))            # X[k₁,k₂]: Savings from merging route R[k₁] into route R[k₂] (R[k₂] --- R[k₁])
    while true
        z = f(s)
        # Step 3.1: Estimate savings from merging every route into every other route
        for (k₁,r₁) ∈ pairs(R)
            for (k₂,r₂) ∈ pairs(R)
                if isequal(r₁, r₂) continue end
                if isopt(r₁) || isopt(r₂) continue end
                # Step 3.1.1: Merge route r₁ into r₂ (r₂ --- r₁)
                cˢ  = C[r₁.iˢ]
                cᵉ  = C[r₁.iᵉ] 
                c   = cˢ
                nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                nᵗ₂ = C[r₂.iᵉ]
                nʰ₂ = D[r₂.iᵈ]
                while true
                    removenode!(c, nᵗ₁, nʰ₁, r₁, s)
                    insertnode!(c, nᵗ₂, nʰ₂, r₂, s)
                    if isequal(c, cᵉ) break end
                    c   = C[r₁.iˢ]
                    nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                    nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                    nᵗ₂ = C[r₂.iᵉ]
                    nʰ₂ = D[r₂.iᵈ]
                end
                # Step 3.1.2: Evaluate savings
                z′  = f(s)
                Δ   = z′ - z
                X[k₁,k₂] = Δ
                # Step 3.1.3. Split routes r₁ and r₂
                c   = cᵉ
                nᵗ₁ = D[r₁.iᵈ]
                nʰ₁ = D[r₁.iᵈ]
                nᵗ₂ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                nʰ₂ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                while true
                    removenode!(c, nᵗ₂, nʰ₂, r₂, s)
                    insertnode!(c, nᵗ₁, nʰ₁, r₁, s)
                    if isequal(c, cˢ) break end
                    c   = C[r₂.iᵉ]
                    nᵗ₁ = D[r₁.iᵈ]
                    nʰ₁ = C[r₁.iˢ]
                    nᵗ₂ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
                    nʰ₂ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
                end
            end
        end
        # Step 3.2: If no merger renders savings, go to step 4. else go to step 3.3.
        k₁,k₂ = Tuple(argmin(X))
        Δ = X[k₁,k₂]
        if Δ > 0 break end 
        # Step 3.3: Merge the routes that render the most savings 
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
            removenode!(c, nᵗ₁, nʰ₁, r₁, s)
            insertnode!(c, nᵗ₂, nʰ₂, r₂, s)
            if isequal(c, cᵉ) break end
            c   = C[r₁.iˢ]
            nᵗ₁ = c.iᵗ ≤ lastindex(D) ? D[c.iᵗ] : C[c.iᵗ]
            nʰ₁ = c.iʰ ≤ lastindex(D) ? D[c.iʰ] : C[c.iʰ]
            nᵗ₂ = C[r₂.iᵉ]
            nʰ₂ = D[r₂.iᵈ]
        end
        X .= Inf
    end
    postinitialize!(s)
    # Step 4: Return solution
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