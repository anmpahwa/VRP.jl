"""
    localsearch!([rng::AbstractRNG], k̅::Int, s::Solution, method::Symbol)

Returns solution `s` after performing local seach on the solution using given `method` for `k̅` iterations.

Available methods include,
- intra-move    : `:intramove!`
- inter-move    : `:intermove!`
- intra-swap    : `:intraswap!`
- inter-swap    : `:interswap!`
- intra-opt     : `:intraopt!`
- inter-opt     : `:interopt!`
- swapdepot     : `:swapdepot!`

Optionally specify a random number generator `rng` as the first argument (defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k̅::Int, s::Solution, method::Symbol)::Solution = isdefined(VRP, method) ? getfield(VRP, method)(rng, k̅, s) : getfield(Main, method)(rng, k̅, s)
localsearch!(k̅::Int, s::Solution, method::Symbol) = localsearch!(Random.GLOBAL_RNG, k̅, s, method)



"""
    intramove!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after moving a randomly selected customer node 
to its best position in the same route if the move results in a reduction 
in objective function value, repeating for `k̅` iterations.
"""
function intramove!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    D = s.D
    C = s.C
    W = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        z  = f(s)
        # Step 2.1: Randomly select a customer node
        c  = sample(rng, C, OffsetWeights(W))
        if isdormant(c) continue end
        # Step 2.2: Remove this node from its position between tail node nᵗ and head node nʰ
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
        removenode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Iterate through all position in the route
        x  = 0.
        p  = (nᵗ.iⁿ, nʰ.iⁿ)
        d  = s.D[r.iᵈ]
        nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ] 
        nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
        nᵗ = d
        nʰ = nˢ
        while true
            # Step 2.3.1: Insert customer node c between tail node nᵗ and head node nʰ
            insertnode!(c, nᵗ, nʰ, r, s)
            # Step 2.3.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 2.3.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
            # Step 2.3.4: Remove node from its position between tail node nᵗ and head node nʰ
            removenode!(c, nᵗ, nʰ, r, s)
            if isequal(nᵗ, nᵉ) break end
            nᵗ = nʰ
            nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
        end
        # Step 2.4: Move the node to its best position (this could be its original position as well)
        iᵗ = p[1]
        iʰ = p[2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end
"""
    intermove!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after moving a randomly selected customer node 
to its best position in another route if the move results in a reduction 
in objective function value, repeating for `k̅` iterations.
"""
function intermove!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    D  = s.D
    C  = s.C
    R  = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    Wᶜ = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        z  = f(s)
        # Step 2.1: Select a random customer node
        c  = sample(rng, C, OffsetWeights(Wᶜ))
        r¹ = c.r
        if isdormant(r¹) continue end
        # Step 2.2: Select a random route
        Wʳ = [isdormant(r²) || isequal(r¹, r²) ? 0 : 1 for r² ∈ R]
        r² = sample(rng, R, Weights(Wʳ))
        if isdormant(r²) continue end
        # Step 2.3: Remove this node from its position between tail node nᵗ and head node nʰ
        nᵗ = isequal(r¹.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r¹.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
        removenode!(c, nᵗ, nʰ, r¹, s)
        # Step 2.4: Iterate through all position in the route
        x  = 0.
        p  = (nᵗ.iⁿ, nʰ.iⁿ)
        r  = r¹
        d  = s.D[r².iᵈ]
        nˢ = isopt(r²) ? C[r².iˢ] : D[r².iˢ] 
        nᵉ = isopt(r²) ? C[r².iᵉ] : D[r².iᵉ]
        nᵗ = d
        nʰ = nˢ
        while true
            # Step 2.4.1: Insert customer node c between tail node nᵗ and head node nʰ
            insertnode!(c, nᵗ, nʰ, r², s)
            # Step 2.4.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 2.4.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p, r = Δ, (nᵗ.iⁿ, nʰ.iⁿ), r² end
            # Step 2.4.4: Remove node from its position between tail node nᵗ and head node nʰ
            removenode!(c, nᵗ, nʰ, r², s)
            if isequal(nᵗ, nᵉ) break end
            nᵗ = nʰ
            nʰ = isequal(r².iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
        end
        # Step 2.5: Move the node to its best position (this could be its original position as well)
        iᵗ = p[1]
        iʰ = p[2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end



"""
    intraswap!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after swapping two randomly selected customers from 
the same route if the swap results in a reduction in objective function 
value, repeating for `k̅` iterations.
"""
function intraswap!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    z  = f(s)
    D  = s.D
    C  = s.C
    W² = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Swap two randomly selected customer nodes
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶
        n² = sample(rng, C, OffsetWeights(W²))
        if isdormant(n²) continue end
        W⁵ = [isdormant(n⁵) || !isequal(n².r, n⁵.r) || isequal(n², n⁵) ? 0. : relatedness(n², n⁵, s) for n⁵ ∈ C]
        n⁵ = sample(rng, C, OffsetWeights(W⁵))
        if isdormant(n⁵) continue end
        r² = n².r
        r⁵ = n⁵.r
        n¹ = isequal(r².iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r².iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        n⁴ = isequal(r⁵.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r⁵.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ]
        if isequal(n², n⁵) continue end
        # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
        if isequal(n³, n⁵)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁵, n⁶, r⁵, s)
        # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
        elseif isequal(n², n⁶)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁴, n⁵, r⁵, s)
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
        else 
            removenode!(n², n¹, n³, r², s)
            removenode!(n⁵, n⁴, n⁶, r⁵, s)
            insertnode!(n⁵, n¹, n³, r², s)
            insertnode!(n², n⁴, n⁶, r⁵, s)
        end
        # Step 2.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reswap the two customer nodes and go to step 1.1
        else
            # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
            if isequal(n³, n⁵)
                removenode!(n², n⁵, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
            elseif isequal(n², n⁶)
                removenode!(n², n⁴, n⁵, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
            else 
                removenode!(n⁵, n¹, n³, r², s)
                removenode!(n², n⁴, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
                insertnode!(n⁵, n⁴, n⁶, r⁵, s)
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end
"""
    interswap!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after swapping two randomly selected customers from 
different routes if the swap results in a reduction in objective function 
value, repeating for `k̅` iterations.
"""
function interswap!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    z  = f(s)
    D  = s.D
    C  = s.C
    W² = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Swap two randomly selected customer nodes
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶
        n² = sample(rng, C, OffsetWeights(W²))
        if isdormant(n²) continue end
        W⁵ = [isdormant(n⁵) || isequal(n².r, n⁵.r) || isequal(n², n⁵) ? 0. : relatedness(n², n⁵, s) for n⁵ ∈ C]
        n⁵ = sample(rng, C, OffsetWeights(W⁵))
        if isdormant(n⁵) continue end
        r² = n².r
        r⁵ = n⁵.r
        n¹ = isequal(r².iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r².iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        n⁴ = isequal(r⁵.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r⁵.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ]
        if isequal(n², n⁵) continue end
        # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
        if isequal(n³, n⁵)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁵, n⁶, r⁵, s)
        # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
        elseif isequal(n², n⁶)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁴, n⁵, r⁵, s)
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
        else 
            removenode!(n², n¹, n³, r², s)
            removenode!(n⁵, n⁴, n⁶, r⁵, s)
            insertnode!(n⁵, n¹, n³, r², s)
            insertnode!(n², n⁴, n⁶, r⁵, s)
        end
        # Step 2.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reswap the two customer nodes and go to step 1.1
        else
            # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
            if isequal(n³, n⁵)
                removenode!(n², n⁵, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
            elseif isequal(n², n⁶)
                removenode!(n², n⁴, n⁵, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
            else 
                removenode!(n⁵, n¹, n³, r², s)
                removenode!(n², n⁴, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
                insertnode!(n⁵, n⁴, n⁶, r⁵, s)
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end



"""
    intraopt!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from the same route 
and reconfiguring them (total possible reconfigurations 2²-1 = 3) if the 
reconfiguration results in a reduction in objective function value, repeating 
for `k̅` iterations.
"""
function intraopt!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    z = f(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    W = [!isopt(r) || isdormant(r) ? 0 : 1 for r ∈ R]
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Iteratively take 2 arcs from the same route
        # d → ... → n¹ → n² → n³ → ... → n⁴ → n⁵ → n⁶ → ... → d
        r = sample(rng, R, Weights(W))
        if !isopt(r) || isdormant(r) continue end
        (i,j) = sample(rng, 1:r.n, 2)
        (i,j) = j < i ? (j,i) : (i,j)  
        k  = 1
        c  = C[r.iˢ]
        n² = c
        n⁵ = c
        while true
            if isequal(k, i) n² = c end
            if isequal(k, j) n⁵ = c end
            if isequal(k, j) break end
            k += 1
            c  = C[c.iʰ]
        end
        n¹ = isequal(r.iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r.iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        n⁴ = isequal(r.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ] 
        if isequal(n², n⁵) || isequal(n¹, n⁵) continue end 
        # Step 2.2: Reconfigure
        # d → ... → n¹ → n⁵ → n⁴ → ... → n³ → n² → n⁶ → ... → d
        n  = n²
        tᵒ = n¹
        hᵒ = n³
        tⁿ = n⁵
        hⁿ = n⁶
        while true
            removenode!(n, tᵒ, hᵒ, r, s)
            insertnode!(n, tⁿ, hⁿ, r, s)
            hⁿ = n
            n  = hᵒ
            hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
            if isequal(n, n⁵) break end
        end
        # Step 2.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 2.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 2.5: Reconfigure back to the original state
        else
            # d → ... → n¹ → n² → n³ → ... → n⁴ → n⁵ → n⁶ → ... → d
            n  = n⁵
            tᵒ = n¹
            hᵒ = n⁴
            tⁿ = n²
            hⁿ = n⁶
            while true
                removenode!(n, tᵒ, hᵒ, r, s)
                insertnode!(n, tⁿ, hⁿ, r, s)
                hⁿ = n
                n  = hᵒ
                hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
                if isequal(n, n²) break end
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end
"""
    interopt!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from the different 
routes and reconfiguring them (total possible reconfigurations 2²-1 = 3) 
if the reconfiguration results in a reduction in objective function value, 
repeating for `k̅` iterations.
"""
function interopt!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    z  = f(s)
    D  = s.D
    C  = s.C
    R  = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    W² = [!isopt(r²) || isdormant(r²) ? 0 : 1 for r² ∈ R]
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Iteratively take 2 arcs from different routes
        # d² → ... → n¹ → n² → n³ → ... → d² and d⁵ → ... → n⁴ → n⁵ → n⁶ → ... → d⁵
        r² = sample(rng, R, Weights(W²))
        if !isopt(r²) || isdormant(r²) continue end
        W⁵ = [!isopt(r⁵) || isdormant(r⁵) || isequal(r², r⁵)  ? 0. : relatedness(r², r⁵, s) for r⁵ ∈ R]
        r⁵ = sample(rng, R, Weights(W⁵))
        if !isopt(r⁵) || isdormant(r⁵) continue end
        d² = D[r².iᵈ]
        d⁵ = D[r⁵.iᵈ]
        i  = rand(rng, 1:r².n)
        k  = 1
        c² = C[r².iˢ]
        n² = c²
        while true
            if isequal(k, i) n² = c² end
            if isequal(k, i) break end
            k += 1
            c² = C[c².iʰ]
        end
        n¹ = isequal(r².iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r².iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        j  = rand(rng, 1:r⁵.n)
        k  = 1
        c⁵ = C[r⁵.iˢ]
        n⁵ = c⁵
        while true
            if isequal(k, j) n⁵ = c⁵ end
            if isequal(k, j) break end
            k += 1
            c⁵ = C[c⁵.iʰ]
        end
        n⁴ = isequal(r⁵.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r⁵.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ]
        # Step 2.2: Reconfigure
        # d² → ... → n¹ → n⁵ → n⁶ → ...  → d² and d⁵ → ... → n⁴ → n² → n³ → ... → d⁵
        c² = n²
        tᵒ = n¹
        hᵒ = n³
        tⁿ = n⁴
        hⁿ = n⁵
        while true
            removenode!(c², tᵒ, hᵒ, r², s)
            insertnode!(c², tⁿ, hⁿ, r⁵, s)
            if isequal(hᵒ, d²) break end
            tⁿ = c² 
            c² = C[hᵒ.iⁿ]
            hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
        end
        c⁵ = n⁵
        tᵒ = c²
        hᵒ = n⁶
        tⁿ = n¹
        hⁿ = d²
        while true
            removenode!(c⁵, tᵒ, hᵒ, r⁵, s)
            insertnode!(c⁵, tⁿ, hⁿ, r², s)
            if isequal(hᵒ, d⁵) break end
            tⁿ = c⁵
            c⁵ = C[hᵒ.iⁿ]
            hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
        end
        # Step 2.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 2.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 2.5: Reconfigure back to the original state
        else
            # d² → ... → n¹ → n² → n³ → ... → d² and d⁵ → ... → n⁴ → n⁵ → n⁶ → ... → d⁵
            c² = n⁵
            tᵒ = n¹
            hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
            tⁿ = n⁴
            hⁿ = n²
            while true
                removenode!(c², tᵒ, hᵒ, r², s)
                insertnode!(c², tⁿ, hⁿ, r⁵, s)
                if isequal(hᵒ, d²) break end
                tⁿ = c² 
                c² = C[hᵒ.iⁿ]
                hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
            end
            c⁵ = n²
            tᵒ = c²
            hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
            tⁿ = n¹
            hⁿ = d²
            while true
                removenode!(c⁵, tᵒ, hᵒ, r⁵, s)
                insertnode!(c⁵, tⁿ, hⁿ, r², s)
                if isequal(hᵒ, d⁵) break end
                tⁿ = c⁵
                c⁵ = C[hᵒ.iⁿ]
                hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end



"""
    swapdepot!(rng::AbstractRNG, k̅::Int, s::Solution)

Returns solution `s` after swapping vehicles, routes, and customer nodes
between two randomly selected depot nodes if the swap results in a reduction 
in objective function value, repeating for `k̅` iterations.
"""
function swapdepot!(rng::AbstractRNG, k̅::Int, s::Solution)
    # Step 1: Initialize
    prelocalsearch!(s)
    z  = f(s)
    D  = s.D
    C  = s.C
    W¹ = isopt.(D)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Select a random depot pair
        d¹ = sample(rng, D, Weights(W¹))
        R¹ = [r¹ for v¹ ∈ d¹.V for r¹ ∈ v¹.R]
        if any(isdormant, R¹) continue end
        W² = [isequal(d¹, d²) ? 0. : relatedness(d¹, d², s) for d² ∈ D]
        d² = sample(rng, D, Weights(W²))
        R² = [r² for v² ∈ d².V for r² ∈ v².R]
        if any(isdormant, R²) continue end
        if isequal(d¹, d²) continue end
        if !isopt(d¹) && !isopt(d²) continue end
        # Step 2.2: Swap vehicles, routes, and customer nodes
        I¹ = eachindex(d¹.V)
        I² = eachindex(d².V)
        while !isempty(d¹.V)
            v = d¹.V[1]
            removevehicle!(v, d¹, s)
            insertvehicle!(v, d², s)
        end
        for iᵛ ∈ I²
            v = d².V[1]
            removevehicle!(v, d², s)
            insertvehicle!(v, d¹, s)
        end
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reconfigure back to the original state
        else
            I¹ = eachindex(d¹.V)
            I² = eachindex(d².V)
            while !isempty(d¹.V)
                v = d¹.V[1]
                removevehicle!(v, d¹, s)
                insertvehicle!(v, d², s)
            end
            for iᵛ ∈ I²
                v = d².V[1]
                removevehicle!(v, d², s)
                insertvehicle!(v, d¹, s)
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end