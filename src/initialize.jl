"""
    build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
Returns a tuple of depot nodes, customer nodes, and arcs for the `instance`.

Note, `dir` locates the the folder containing instance files as sub-folders.

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
        tᵃ = 0.
        tᵈ = 0.
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
    φᵈ = iszero(getproperty.(D, :tˢ)) && iszero(getproperty.(D, :tᵉ))
    φᶜ = iszero(getproperty.(C, :tᵉ)) && iszero(getproperty.(C, :tˡ))
    φᵛ = iszero(getproperty.(V, :τʷ)) && iszero(getproperty.(V, :πᵗ))
    global φᵉ = !(φᵈ && φᶜ && φᵛ)::Bool
    G  = (D, C, A)
    return G
end



"""
    cluster(rng::AbstractRNG, k::Int, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns `Solution` created using `k`-means clustering algorithm.
Here, each cluster is assigned to the nearest depot node and 
then each customer in this cluster is best inserted to the 
assigned depot node.

Note, `dir` locates the the folder containing instance files as sub-folders.

    <dir>
    |-<instance>
        |-arcs.csv
        |-depot_nodes.csv
        |-customer_nodes.csv
        |-vehicles.csv
"""
function cluster(rng::AbstractRNG, k::Int, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1: Initialize
    G = build(instance; dir=dir)
    s = Solution(G...)
    preinitialize!(s)
    D = s.D
    C = s.C
    # Step 2: Clustering
    N = zeros(4, eachindex(C))
    for (iⁿ,c) ∈ pairs(C) N[:,iⁿ] = [c.x, c.y, c.tᵉ, c.tˡ] end
    K = kmeans(N.parent, k; rng=rng)
    A = OffsetVector(K.assignments, eachindex(C))
    M = K.centers
    # Step 3: Add customers from each cluster to the assigned depot
    for k ∈ 1:nclusters(K)
        Y = fill(Inf, length(D))  
        for j ∈ eachindex(Y)
            d  = D[j]
            xᵒ = M[1,k]
            yᵒ = M[2,k]
            xᵈ = d.x
            yᵈ = d.y
            Y[j] = sqrt((xᵒ-xᵈ)^2 + (yᵒ-yᵈ)^2)
        end
        d = D[argmin(Y)]
        R = [r for v ∈ d.V for r ∈ v.R]
        L = filter(c -> isequal(A[c.iⁿ], k), C)
        if isempty(L) continue end
        I = eachindex(L)
        J = eachindex(R)
        W = ones(Int, I)                        # W[j]  : selection weight for customer node L[i]
        X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
        P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
        # Step 3.1: Iterate until all open customer nodes have been inserted into the route
        for _ ∈ I
            if !hasslack(d) break end
            # Step 3.1.1: Iterate through all open customer nodes and every possible insertion position in each route
            z = f(s)
            i = sample(rng, I, Weights(W))
            c = L[i]
            for (j,r) ∈ pairs(R)
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 3.1.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 3.1.1.2: Compute the insertion cost
                    z′ = f(s)
                    Δ  = z′ - z
                    # Step 3.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 3.1.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
            # Step 3.1.2: Randomly select a customer node to insert at its best position
            j  = argmin(X[i,:])
            r  = R[j]
            d  = s.D[r.iᵈ]
            v  = d.V[r.iᵛ]
            iᵗ = P[i,j][1]
            iʰ = P[i,j][2]
            nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
            nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
            insertnode!(c, nᵗ, nʰ, r, s)
            # Step 3.1.3: Revise vectors appropriately
            W[i] = 0
            # Step 3.1.4: Update solution appropriately     
            if addroute(r, s)
                r = Route(v, d)
                push!(v.R, r)
                push!(R, r)
                append!(X, fill(Inf, (I,1)))
                append!(P, fill((0, 0), (I,1)))
            end
            if addvehicle(v, s)
                v = Vehicle(v, d)
                r = Route(v, d)
                push!(d.V, v)
                push!(v.R, r) 
                push!(R, r)
                append!(X, fill(Inf, (I,1)))
                append!(P, fill((0, 0), (I,1)))
            end
        end
    end
    if any(isopen, C) best!(rng, s) end
    postinitialize!(s)
    # Step 4: Return initial solution
    return s
end



"""
    initialize([rng::AbstractRNG], instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial VRP `Solution` developed using iterated clustering method. 
If the `method` is set to `:local` search, the number of clusters are 
increased iteratively for at most as many iterations as the number of 
depot nodes. Else if the `method` is set to `:global` search, the number 
of clusters are increased iteratively for at least as many iterations as 
the number of depot nodes and at most the number of customer nodes until a 
feasible solution is found. Finally, the solution with the least objective 
function value is returned as the initial solution.

Note, `dir` locates the the folder containing instance files as sub-folders.

    <dir>
    |-<instance>
        |-arcs.csv
        |-depot_nodes.csv
        |-customer_nodes.csv
        |-vehicles.csv

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function initialize(rng::AbstractRNG, instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1. Initialize
    s = Solution(build(instance; dir=dir)...)
    z = Inf
    # Step 2. Iteratively increase the number of clusters
    k = 0
    k̲ = length(s.D)
    k̅ = length(s.C)
    ϕ = isequal(method, :local)
    while k < k̅
        k += 1
        s′ = cluster(rng, k, instance; dir=dir)
        z′ = f(s′)
        # Step 2.1. Update solution
        if z′ < z
            z = z′ 
            s = deepcopy(s′)
        end
        # Step 2.2. Check for break conditions
        ϕ = ϕ || isfeasible(s′)
        k < k̲ ? continue : (ϕ ? break : continue)
    end
    # Step 3. Return solution
    return s
end
initialize(instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances")) = initialize(Random.GLOBAL_RNG, instance; method=:method, dir=dir)