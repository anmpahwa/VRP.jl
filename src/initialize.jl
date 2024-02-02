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
        v̅  = df[k,5]
        tˢ = df[k,6]
        tᵉ = df[k,7]
        n  = 0
        q  = 0.
        l  = 0.
        πᵒ = df[k,8]
        πᶠ = df[k,9]
        d  = DepotNode(iⁿ, x, y, qᵈ, v̅, tˢ, tᵉ, n, q, l, πᵒ, πᶠ, Vehicle[])
        D[iⁿ] = d
    end
    # Customer nodes
    df = DataFrame(CSV.File(joinpath(dir, "$instance/customer_nodes.csv")))
    I  = (df[1,1]:df[nrow(df),1])
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
    regret([rng::AbstractRNG], instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial `Solution` using regret-2 insertion method. 

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
function regret(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1: Initialize
    k̅ = 2
    G = build(instance; dir=dir)
    s = Solution(G...)
    preinitialize!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c) && isdelivery(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))                 # X[i,j]: insertion cost of delivery node L[i] and its associated pickup node at their best position in route R[j]
    P = ElasticMatrix(fill(((0, 0), (0, 0)), (I,J)))    # P[i,j]: best insertion postion of associated pickup node and the delivery node L[i] in route R[j]
    Y = fill(Inf, (I,k̅))                                # Y[i,k]: insertion cost of delivery node L[i] and its associated pickup node at kᵗʰ best position
    ϕ = ones(Int, J)                                    # ϕ[j]  : binary weight for route R[j]
    N = zeros(Int, (I,k̅))                               # N[i,k]: route index of delivery node L[i] and its associated pickup node at kᵗʰ best position
    Z = fill(-Inf, I)                                   # Z[i]  : regret-N cost of delivery node L[i] and its associated pickup node
    # Step 2: Iterate until all open delivery nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open delivery nodes (and the associated pickup nodes)
        z = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            for (j,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[j]) continue end
                d   = s.D[r.iᵈ]
                nᵖˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵖᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵖᵗ = d
                nᵖʰ = nᵖˢ
                while true
                    # Step 2.1.1.1: Insert associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                    insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                    nᵈˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                    nᵈᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                    nᵈᵗ = d
                    nᵈʰ = nᵈˢ
                    while true
                        # Step 2.1.1.1.1: Insert delivery node cᵈ between tail node nᵈᵗ and head node nᵈʰ in route r
                        insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                        # Step 2.1.1.1.2: Compute the insertion cost
                        z′ = f(s)
                        Δ  = z′ - z
                        # Step 2.1.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                        if Δ < X[i,j] X[i,j], P[i,j] = Δ, ((nᵖᵗ.iⁿ, nᵖʰ.iⁿ), (nᵈᵗ.iⁿ, nᵈʰ.iⁿ)) end
                        # Step 2.1.1.1.4: Revise N least insertion costs
                        k̲ = 1
                        for k ∈ 1:k̅ 
                            k̲ = k
                            if Δ < Y[i,k] break end
                        end
                        for k ∈ k̅:-1:k̲ 
                            Y[i,k] = isequal(k, k̲) ? Δ : Y[i,k-1]
                            N[i,k] = isequal(k, k̲) ? r.iʳ : N[i,k-1]
                        end
                        # Step 2.1.1.1.5: Remove delivery node cᵈ from its position between tail node nᵈᵗ and head node nᵈʰ in route r
                        removenode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                        if isequal(nᵈᵗ, nᵈᵉ) break end
                        nᵈᵗ = nᵈʰ
                        nᵈʰ = isequal(r.iᵉ, nᵈᵗ.iⁿ) ? D[nᵈᵗ.iʰ] : C[nᵈᵗ.iʰ]
                    end
                    # Step 2.1.1.2: Remove associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                    removenode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                    if isequal(nᵖᵗ, nᵖᵉ) break end
                    nᵖᵗ = nᵖʰ
                    nᵖʰ = isequal(r.iᵉ, nᵖᵗ.iⁿ) ? D[nᵖᵗ.iʰ] : C[nᵖᵗ.iʰ]
                end
            end
            # Step 2.1.2: Compute regret cost for delivery node L[i]
            Z[i] = 0.
            for k ∈ 1:k̅ Z[i] += Y[i,k] - Y[i,1] end
        end
        # Step 2.2: Insert delivery node and the associated pickup node with highest regret cost in its best position (break ties by inserting the nodes with the lowest insertion cost)
        I̲  = findall(isequal.(Z, maximum(Z)))
        i,j= Tuple(argmin(X[I̲,:]))
        i  = I̲[i]
        c  = L[i]
        cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
        cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
        r  = R[j]
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        iᵖᵗ = P[i,j][1][1]
        iᵖʰ = P[i,j][1][2]
        iᵈᵗ = P[i,j][2][1]
        iᵈʰ = P[i,j][2][2]
        nᵖᵗ = iᵖᵗ ≤ lastindex(D) ? D[iᵖᵗ] : C[iᵖᵗ]
        nᵖʰ = iᵖʰ ≤ lastindex(D) ? D[iᵖʰ] : C[iᵖʰ]
        nᵈᵗ = iᵈᵗ ≤ lastindex(D) ? D[iᵈᵗ] : C[iᵈᵗ]
        nᵈʰ = iᵈʰ ≤ lastindex(D) ? D[iᵈʰ] : C[iᵈʰ]
        insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
        insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        X[i,:] .= Inf
        P[i,:] .= (((0 ,0), (0, 0)), )
        Y[i,:] .= Inf
        N[i,:] .= 0
        Z .= -Inf 
        for (i,c) ∈ pairs(L)
            for k ∈ 1:k̅
                if iszero(N[i,k]) break end
                j = findfirst(r -> isequal(r.iʳ, N[i,k]), R)
                r = R[j]
                if isequal(r.iᵛ, v.iᵛ) Y[i,k], N[i,k] = Inf, 0 end
            end
            K = sortperm(Y[i,:])
            Y[i,:] .= Y[i,K]
            N[i,:] .= N[i,K]
        end
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            φʳ = isequal(r, c.r)
            φᵛ = isequal(r.iᵛ, v.iᵛ) && isless(c.r.tⁱ, r.tⁱ)
            φᵈ = isequal(r.iᵈ, d.iⁿ) && !hasslack(d)
            φˢ = φʳ || φᵛ || φᵈ
            if isequal(φˢ, false) continue end
            X[:,j] .= Inf
            ϕ[j] = 1  
        end
        # Step 2.4: Update solution appropriately     
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r)
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill(((0, 0), (0, 0)), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill(((0, 0), (0, 0)), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinitialize!(s)
    # Step 3: Return solution
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
initialize(rng::AbstractRNG, instance::String; dir=joinpath(dirname(@__DIR__), "instances")) = regret(rng, instance; dir=dir)
initialize(instance::String; dir=joinpath(dirname(@__DIR__), "instances")) = initialize(Random.GLOBAL_RNG, instance; dir=dir)