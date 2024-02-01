"""
    insert!([rng::AbstractRNG], s::Solution, method::Symbol)

Returns solution `s` after inserting open customer nodes to the solution using the given `method`.

Available methods include,
- Best Insertion           : `:best!`
- Precise Greedy Insertion : `:precise!`
- Perturb Greedy insertion : `:perturb!`
- Regret-two Insertion     : `:regret2!`
- Regret-three Insertion   : `:regret3!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Symbol)::Solution = isdefined(VRP, method) ? getfield(VRP, method)(rng, s) : getfield(Main, method)(rng, s)
insert!(s::Solution, method::Symbol) = insert!(Random.GLOBAL_RNG, s, method)



"""
    best!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting randomly selected delivery node 
and its associated pickup node at their best position until all open 
delivery nodes have been inserted to the solution.
"""
function best!(rng::AbstractRNG, s::Solution)
    # Step 1: Initialize
    preinsert!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c) && isdelivery(c)]
    I = eachindex(L)
    J = eachindex(R)
    W = ones(Int, I)                                    # W[j]  : selection weight for delivery node L[i]
    X = ElasticMatrix(fill(Inf, (I,J)))                 # X[i,j]: insertion cost of delivery node L[i] and its associated pickup node at their best position in route R[j]
    P = ElasticMatrix(fill(((0, 0), (0, 0)), (I,J)))    # P[i,j]: best insertion postion of associated pickup node and the delivery node L[i] in route R[j]
    # Step 2: Iterate until all open delivery nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Randomly select an open delivery node (and the associated pickup node) and iterate through all possible insertion positions in each route
        z  = f(s)
        i  = sample(rng, I, Weights(W))
        c  = L[i]
        cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
        cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
        for (j,r) ∈ pairs(R)
            d   = s.D[r.iᵈ]
            nᵖˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
            nᵖᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
            nᵖᵗ = d
            nᵖʰ = nᵖˢ
            while true
                # Step 2.1.1: Insert associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                nᵈˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵈᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵈᵗ = d
                nᵈʰ = nᵈˢ
                while true
                    # Step 2.1.1.1: Insert delivery node cᵈ between tail node nᵈᵗ and head node nᵈʰ in route r
                    insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z′ = f(s)
                    Δ  = z′ - z
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, ((nᵖᵗ.iⁿ, nᵖʰ.iⁿ), (nᵈᵗ.iⁿ, nᵈʰ.iⁿ)) end
                    # Step 2.1.1.4: Remove delivery node cᵈ from its position between tail node nᵈᵗ and head node nᵈʰ in route r
                    removenode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                    if isequal(nᵈᵗ, nᵈᵉ) break end
                    nᵈᵗ = nᵈʰ
                    nᵈʰ = isequal(r.iᵉ, nᵈᵗ.iⁿ) ? D[nᵈᵗ.iʰ] : C[nᵈᵗ.iʰ]
                end
                # Step 2.1.2: Remove associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                removenode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                if isequal(nᵖᵗ, nᵖᵉ) break end
                nᵖᵗ = nᵖʰ
                nᵖʰ = isequal(r.iᵉ, nᵖᵗ.iⁿ) ? D[nᵖᵗ.iʰ] : C[nᵖᵗ.iʰ]
            end
        end
        # Step 2.2: Insert the delivery node and the associated pickup node at their best position
        j = argmin(X[i,:])
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
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
        W[i] = 0
        X[i,:] .= Inf
        P[i,:] .= (((0 ,0), (0, 0)), )
        # Step 2.4: Update solution appropriately 
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r)
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill(((0, 0), (0, 0)), (I,1)))
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill(((0, 0), (0, 0)), (I,1)))
        end
    end
    postinsert!(s)
    # Step 3: Return solution
    return s
end



"""
    greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)

Returns solution `s` after iteratively inserting delivery node 
and its associated pickup node with least insertion cost until all 
open delivery nodes have been added to the solution. 
Available modes include `:pcs` (precise estimation of insertion cost) 
and `:ptb` (perturbed estimation of insertion cost).
"""
function greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)
    # Step 1: Initialize
    preinsert!(s)
    D = s.D
    C = s.C
    φ = isequal(mode,  :ptb)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c) && isdelivery(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))                 # X[i,j]: insertion cost of delivery node L[i] and its associated pickup node at their best position in route R[j]
    P = ElasticMatrix(fill(((0, 0), (0, 0)), (I,J)))    # P[i,j]: best insertion postion of associated pickup node and the delivery node L[i] in route R[j]
    ϕ = ones(Int, J)                                    # ϕ[j]  : binary weight for route R[j]
    # Step 2: Iterate until all open delivery nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open delivery nodes (and the associated pickup nodes) through every possible insertion position in each route
        z = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d   = s.D[r.iᵈ]
                nᵖˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵖᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵖᵗ = d
                nᵖʰ = nᵖˢ
                while true
                    # Step 2.1.1: Insert associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                    insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                    nᵈˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                    nᵈᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                    nᵈᵗ = d
                    nᵈʰ = nᵈˢ
                    while true
                        # Step 2.1.1.1: Insert delivery node cᵈ between tail node nᵈᵗ and head node nᵈʰ in route r
                        insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                        # Step 2.1.1.2: Compute the insertion cost
                        z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                        Δ  = z′ - z
                        # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                        if Δ < X[i,j] X[i,j], P[i,j] = Δ, ((nᵖᵗ.iⁿ, nᵖʰ.iⁿ), (nᵈᵗ.iⁿ, nᵈʰ.iⁿ)) end
                        # Step 2.1.1.4: Remove delivery node cᵈ from its position between tail node nᵈᵗ and head node nᵈʰ in route r
                        removenode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
                        if isequal(nᵈᵗ, nᵈᵉ) break end
                        nᵈᵗ = nᵈʰ
                        nᵈʰ = isequal(r.iᵉ, nᵈᵗ.iⁿ) ? D[nᵈᵗ.iʰ] : C[nᵈᵗ.iʰ]
                    end
                    # Step 2.1.2: Remove associated pickup node cᵖ between tail node nᵖᵗ and head node nᵖʰ in route r
                    removenode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
                    if isequal(nᵖᵗ, nᵖᵉ) break end
                    nᵖᵗ = nᵖʰ
                    nᵖʰ = isequal(r.iᵉ, nᵖᵗ.iⁿ) ? D[nᵖᵗ.iʰ] : C[nᵖᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Insert the delivery node and the associated pickup node with least insertion cost at their best position
        i,j= Tuple(argmin(X))
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
    postinsert!(s)
    # Step 3: Return solution
    return s
end
"""
    precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively inserting delivery node 
and its associated pickup node with least insertion cost until all 
open delivery nodes have been added to the solution. 
Estimates insertion cost precisely.
"""
precise!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:pcs)
"""
    precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively inserting delivery node 
and its associated pickup node with least insertion cost until all 
open delivery nodes have been added to the solution.  
Estimates insertion cost with a perturbration.
"""
perturb!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:ptb)



"""
    regretk!(rng::AbstractRNG, s::Solution, k̅::Int)

Returns solution `s` after iteratively adding delivery node 
and its associated pickup node with highest regret-k cost at 
its best position until all open delivery nodes have been 
added to the solution.
Note, regretk mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
function regretk!(rng::AbstractRNG, s::Solution, k̅::Int)
    # Step 1: Initialize
    preinsert!(s)
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
    postinsert!(s)
    # Step 3: Return solution
    return s
end
"""
    regret2!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively adding customer nodes with 
highest regret-2 cost at its best position until all open customer 
nodes have been added to the solution.
Note, regret2 mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
regret2!(rng::AbstractRNG, s::Solution) = regretk!(rng, s, 2)
"""
    regret3!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively adding customer nodes with 
highest regret-3 cost at its best position until all open customer 
nodes have been added to the solution.
Note, regret3 mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
regret3!(rng::AbstractRNG, s::Solution) = regretk!(rng, s, 3)