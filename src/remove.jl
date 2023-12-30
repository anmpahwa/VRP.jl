"""
    remove!([rng::AbstractRNG], q::Int, s::Solution, method::Symbol)

Returns solution removing `q` customer nodes from solution s using the given `method`.

Available methods include,
- Random Customer Node Removal  : `:randomcustomer!`
- Related Customer Node Removal : `:relatedcustomer!`
- Worst Customer Node Removal   : `:worstcustomer!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, q::Int, s::Solution, method::Symbol)::Solution = isdefined(VRP, method) ? getfield(VRP, method)(rng, q, s) : getfield(Main, method)(rng, q, s)
remove!(q::Int, s::Solution, method::Symbol) = remove!(Random.GLOBAL_RNG, q, s, method)



"""
    randomcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes
selected randomly.
"""
function randomcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    W = ones(eachindex(C))          # W[iⁿ]: selection weight of customer node C[iⁿ]
    # Step 2: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    while n < q
        iⁿ = sample(rng, eachindex(C), OffsetWeights(W))
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        W[iⁿ] = 0
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end



"""
    relatedcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes
most related to a randomly selected pivot customer node.
"""
function relatedcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    X = fill(-Inf, eachindex(C))    # X[iⁿ]: relatedness of customer node C[iⁿ] with pivot customer node C[i]
    W = ones(eachindex(C))          # W[iⁿ]: selection weight of customer node C[iⁿ]
    # Step 2: Randomly select a pivot customer node
    i = sample(rng, eachindex(C), OffsetWeights(W))
    # Step 3: For each customer node, evaluate relatedness to this pivot customer node
    m = sample(rng, s.φ ? [:q, :l, :t] : [:q, :l])
    for iⁿ ∈ eachindex(C) X[iⁿ] = isone(W[iⁿ]) ? relatedness(m, C[iⁿ], C[i], s) : -Inf end
    # Step 4: Remove q most related customer nodes
    n = 0
    while n < q
        iⁿ = argmax(X)
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        X[iⁿ] = -Inf
    end
    # Step 5: Remove redundant vehicles and routes
    postremove!(s)
    # Step 6: Return solution
    return s
end



"""
    worstcustomer!(rng::AbstractRNG, q::Int, s::Solution)

Returns solution `s` after removing exactly `q` customer nodes 
with highest removal cost (savings).
"""
function worstcustomer!(rng::AbstractRNG, q::Int, s::Solution)
    # Step 1: Initialize
    preremove!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R if isactive(r)]
    L = [c for c ∈ C if isactive(c) && isdelivery(c)]
    X = fill(-Inf, eachindex(L))   # X[i]: removal cost of customer node L[i]
    ϕ = ones(Int, eachindex(R))    # ϕʳ[j]: binary weight for route R[j]
    # Step 2: Iterate until q customer nodes have been removed
    n = 0
    while n < q
        # Step 2.1: For every closed customer node evaluate removal cost
        z = f(s)
        for i ∈ eachindex(L)
            cᵈ = L[i]
            cᵖ = C[cᵈ.jⁿ]
            if isopen(cᵈ) || isopen(cᵖ) continue end
            r = cᵈ.r
            j = findfirst(isequal(r), R)
            if iszero(ϕ[j]) continue end
            # Step 2.1.1: Remove closed customer node c between tail node nᵗ and head node nʰ in route r
            nᵖᵗ = isequal(r.iˢ, cᵖ.iⁿ) ? D[cᵖ.iᵗ] : C[cᵖ.iᵗ]
            nᵖʰ = isequal(r.iᵉ, cᵖ.iⁿ) ? D[cᵖ.iʰ] : C[cᵖ.iʰ]
            removenode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
            nᵈᵗ = isequal(r.iˢ, cᵈ.iⁿ) ? D[cᵈ.iᵗ] : C[cᵈ.iᵗ]
            nᵈʰ = isequal(r.iᵉ, cᵈ.iⁿ) ? D[cᵈ.iʰ] : C[cᵈ.iʰ]
            removenode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
            # Step 2.1.2: Evaluate the removal cost
            z′ = f(s) * (1 + rand(rng, Uniform(-0.2, 0.2)))
            Δ  = z′ - z
            X[i] = -Δ
            # Step 2.1.3: Re-insert customer node c between tail node nᵗ and head node nʰ in route r
            insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
            insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
        end
        # Step 2.2: Remove the customer node with highest removal cost (savings)
        i  = argmax(X)
        cᵈ = L[i]
        cᵖ = C[cᵈ.jⁿ]
        r  = cᵈ.r
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        nᵖᵗ = isequal(r.iˢ, cᵖ.iⁿ) ? D[cᵖ.iᵗ] : C[cᵖ.iᵗ]
        nᵖʰ = isequal(r.iᵉ, cᵖ.iⁿ) ? D[cᵖ.iʰ] : C[cᵖ.iʰ]
        removenode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
        nᵈᵗ = isequal(r.iˢ, cᵈ.iⁿ) ? D[cᵈ.iᵗ] : C[cᵈ.iᵗ]
        nᵈʰ = isequal(r.iᵉ, cᵈ.iⁿ) ? D[cᵈ.iʰ] : C[cᵈ.iʰ]
        removenode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
        tⁱ = r.tⁱ
        n += 2
        # Step 2.3: Update cost and selection weight vectors
        X[i] = -Inf
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            φʳ = isequal(r, cᵈ.r)
            φᵛ = isequal(r.iᵛ, v.iᵛ) && isless(tⁱ, r.tⁱ) && isequal(s.φ, true)
            φᵈ = false
            φˢ = φʳ || φᵛ || φᵈ
            if isequal(φˢ, false) continue end
            ϕ[j] = 1
        end
    end
    postremove!(s)
    # Step 3: Return solution
    return s
end