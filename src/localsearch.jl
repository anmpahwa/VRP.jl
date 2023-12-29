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
        n  = sample(rng, C, OffsetWeights(Wᶜ))
        cᵖ = isdelivery(n) ? C[n.jⁿ] : C[n.iⁿ]
        cᵈ = isdelivery(n) ? C[n.iⁿ] : C[n.jⁿ]
        if !isequal(cᵖ.r, cᵈ.r) continue end
        r₁ = cᵈ.r
        if isdormant(r₁) continue end
        # Step 2.2: Select a random route
        Wʳ = [isdormant(r₂) || isequal(r₁, r₂) ? 0. : 1. for r₂ ∈ R]
        r₂ = sample(rng, R, Weights(Wʳ))
        if isdormant(r₂) continue end
        # Step 2.3: Remove this node from its position between tail node nᵗ and head node nʰ
        nᵈᵗ = isequal(r₁.iˢ, cᵈ.iⁿ) ? D[cᵈ.iᵗ] : C[cᵈ.iᵗ]
        nᵈʰ = isequal(r₁.iᵉ, cᵈ.iⁿ) ? D[cᵈ.iʰ] : C[cᵈ.iʰ]
        removenode!(cᵈ, nᵈᵗ, nᵈʰ, r₁, s)
        nᵖᵗ = isequal(r₁.iˢ, cᵖ.iⁿ) ? D[cᵖ.iᵗ] : C[cᵖ.iᵗ]
        nᵖʰ = isequal(r₁.iᵉ, cᵖ.iⁿ) ? D[cᵖ.iʰ] : C[cᵖ.iʰ]
        removenode!(cᵖ, nᵖᵗ, nᵖʰ, r₁, s)
        # Step 2.4: Iterate through all position in the route
        x  = 0.
        p  = ((nᵖᵗ.iⁿ, nᵖʰ.iⁿ), (nᵈᵗ.iⁿ, nᵈʰ.iⁿ))
        r  = r₁
        d  = s.D[r₂.iᵈ]
        nˢ = isopt(r₂) ? C[r₂.iˢ] : D[r₂.iˢ] 
        nᵉ = isopt(r₂) ? C[r₂.iᵉ] : D[r₂.iᵉ]
        nᵗ = d
        nʰ = nˢ
        while true
            # Step 2.4.1: Insert customer node c between tail node nᵗ and head node nʰ
            insertnode!(cᵖ, nᵗ, nʰ, r₂, s)
            insertnode!(cᵈ, cᵖ, nʰ, r₂, s)
            # Step 2.4.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 2.4.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p, r = Δ, ((nᵗ.iⁿ, nʰ.iⁿ), (cᵖ.iⁿ, nʰ.iⁿ)), r₂ end
            # Step 2.4.4: Remove node from its position between tail node nᵗ and head node nʰ
            removenode!(cᵈ, cᵖ, nʰ, r₂, s)
            removenode!(cᵖ, nᵗ, nʰ, r₂, s)
            if isequal(nᵗ, nᵉ) break end
            nᵗ = nʰ
            nʰ = isequal(r₂.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
        end
        # Step 2.5: Move the node to its best position (this could be its original position as well)
        iᵖᵗ = p[1][1]
        iᵖʰ = p[1][2]
        iᵈᵗ = p[2][1]
        iᵈʰ = p[2][2]
        nᵖᵗ = iᵖᵗ ≤ lastindex(D) ? D[iᵖᵗ] : C[iᵖᵗ]
        nᵖʰ = iᵖʰ ≤ lastindex(D) ? D[iᵖʰ] : C[iᵖʰ]
        nᵈᵗ = iᵈᵗ ≤ lastindex(D) ? D[iᵈᵗ] : C[iᵈᵗ]
        nᵈʰ = iᵈʰ ≤ lastindex(D) ? D[iᵈʰ] : C[iᵈʰ]
        insertnode!(cᵖ, nᵖᵗ, nᵖʰ, r, s)
        insertnode!(cᵈ, nᵈᵗ, nᵈʰ, r, s)
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
    W₂ = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Swap two randomly selected customer nodes
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆
        n₂ = sample(rng, C, OffsetWeights(W₂))
        if isdormant(n₂) continue end
        m  = sample(rng, φᵉ::Bool ? [:q, :l, :t] : [:q, :l])
        W₅ = [isdormant(n₅) || !isequal(n₂.r, n₅.r) || isequal(n₂, n₅) ? 0. : relatedness(m, n₂, n₅, s) for n₅ ∈ C]
        n₅ = sample(rng, C, OffsetWeights(W₅))
        if isdormant(n₅) continue end
        r₂ = n₂.r
        r₅ = n₅.r
        n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
        n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
        n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
        n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
        if isequal(n₂, n₅) continue end
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₅, n₆, r₅, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₅, r₅, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₂, n₁, n₃, r₂, s)
            removenode!(n₅, n₄, n₆, r₅, s)
            insertnode!(n₅, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₆, r₅, s)
        end
        # Step 2.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reswap the two customer nodes and go to step 1.1
        else
            # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
            if isequal(n₃, n₅)
                removenode!(n₂, n₅, n₆, r₅, s)
                insertnode!(n₂, n₁, n₃, r₂, s)
            # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
            elseif isequal(n₂, n₆)
                removenode!(n₂, n₄, n₅, r₅, s)
                insertnode!(n₂, n₁, n₃, r₂, s)
            # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
            else 
                removenode!(n₅, n₁, n₃, r₂, s)
                removenode!(n₂, n₄, n₆, r₅, s)
                insertnode!(n₂, n₁, n₃, r₂, s)
                insertnode!(n₅, n₄, n₆, r₅, s)
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
    W₂ = isactive.(C)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Swap two randomly selected customer nodes
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆
        n₂ = sample(rng, C, OffsetWeights(W₂))
        if isdormant(n₂) continue end
        m  = sample(rng, φᵉ::Bool ? [:q, :l, :t] : [:q, :l])
        W₅ = [isdormant(n₅) || isequal(n₂.r, n₅.r) || isequal(n₂, n₅) ? 0. : relatedness(m, n₂, n₅, s) for n₅ ∈ C]
        n₅ = sample(rng, C, OffsetWeights(W₅))
        if isdormant(n₅) continue end
        if isequal(n₂, n₅) continue end
        r₂ = n₂.r
        r₅ = n₅.r
        n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
        n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
        n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
        n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₅, n₆, r₅, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₅, r₅, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₂, n₁, n₃, r₂, s)
            removenode!(n₅, n₄, n₆, r₅, s)
            insertnode!(n₅, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₆, r₅, s)
        end
        n₂ = C[n₂.jⁿ]
        n₅ = C[n₅.jⁿ]
        r₂ = n₂.r
        r₅ = n₅.r
        n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
        n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
        n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
        n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₅, n₆, r₅, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₅, r₅, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₂, n₁, n₃, r₂, s)
            removenode!(n₅, n₄, n₆, r₅, s)
            insertnode!(n₅, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₆, r₅, s)
        end
        # Step 2.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reswap the two customer nodes and go to step 1.1
        else
            n₂ = C[n₂.iⁿ]
            n₅ = C[n₅.iⁿ]
            r₂ = n₂.r
            r₅ = n₅.r
            n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
            n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
            n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
            n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
            # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
            if isequal(n₃, n₅)
                removenode!(n₂, n₁, n₃, r₂, s)
                insertnode!(n₂, n₅, n₆, r₅, s)
            # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
            elseif isequal(n₂, n₆)
                removenode!(n₂, n₁, n₃, r₂, s)
                insertnode!(n₂, n₄, n₅, r₅, s)
            # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
            else 
                removenode!(n₂, n₁, n₃, r₂, s)
                removenode!(n₅, n₄, n₆, r₅, s)
                insertnode!(n₅, n₁, n₃, r₂, s)
                insertnode!(n₂, n₄, n₆, r₅, s)
            end
            n₂ = C[n₂.jⁿ]
            n₅ = C[n₅.jⁿ]
            r₂ = n₂.r
            r₅ = n₅.r
            n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
            n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
            n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
            n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
            # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
            if isequal(n₃, n₅)
                removenode!(n₂, n₁, n₃, r₂, s)
                insertnode!(n₂, n₅, n₆, r₅, s)
            # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
            elseif isequal(n₂, n₆)
                removenode!(n₂, n₁, n₃, r₂, s)
                insertnode!(n₂, n₄, n₅, r₅, s)
            # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
            else 
                removenode!(n₂, n₁, n₃, r₂, s)
                removenode!(n₅, n₄, n₆, r₅, s)
                insertnode!(n₅, n₁, n₃, r₂, s)
                insertnode!(n₂, n₄, n₆, r₅, s)
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
and reconfiguring them (total possible reconfigurations 2₂-1 = 3) if the 
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
        # d → ... → n₁ → n₂ → n₃ → ... → n₄ → n₅ → n₆ → ... → d
        r = sample(rng, R, Weights(W))
        if !isopt(r) || isdormant(r) continue end
        (i,j) = sample(rng, 1:r.n, 2)
        (i,j) = j < i ? (j,i) : (i,j)  
        k  = 1
        c  = C[r.iˢ]
        n₂ = c
        n₅ = c
        while true
            if isequal(k, i) n₂ = c end
            if isequal(k, j) n₅ = c end
            if isequal(k, j) break end
            k += 1
            c  = C[c.iʰ]
        end
        n₁ = isequal(r.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
        n₃ = isequal(r.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
        n₄ = isequal(r.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
        n₆ = isequal(r.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ] 
        if isequal(n₂, n₅) || isequal(n₁, n₅) continue end 
        # Step 2.2: Reconfigure
        # d → ... → n₁ → n₅ → n₄ → ... → n₃ → n₂ → n₆ → ... → d
        n  = n₂
        tᵒ = n₁
        hᵒ = n₃
        tⁿ = n₅
        hⁿ = n₆
        while true
            removenode!(n, tᵒ, hᵒ, r, s)
            insertnode!(n, tⁿ, hⁿ, r, s)
            hⁿ = n
            n  = hᵒ
            hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
            if isequal(n, n₅) break end
        end
        # Step 2.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 2.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 2.5: Reconfigure back to the original state
        else
            # d → ... → n₁ → n₂ → n₃ → ... → n₄ → n₅ → n₆ → ... → d
            n  = n₅
            tᵒ = n₁
            hᵒ = n₄
            tⁿ = n₂
            hⁿ = n₆
            while true
                removenode!(n, tᵒ, hᵒ, r, s)
                insertnode!(n, tⁿ, hⁿ, r, s)
                hⁿ = n
                n  = hᵒ
                hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
                if isequal(n, n₂) break end
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
routes and reconfiguring them (total possible reconfigurations 2₂-1 = 3) 
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
    W₂ = [!isopt(r₂) || isdormant(r₂) ? 0 : 1 for r₂ ∈ R]
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Iteratively take 2 arcs from different routes
        # d₂ → ... → n₁ → n₂ → n₃ → ... → d₂ and d₅ → ... → n₄ → n₅ → n₆ → ... → d₅
        r₂ = sample(rng, R, Weights(W₂))
        if !isopt(r₂) || isdormant(r₂) continue end
        m  = sample(rng, φᵉ::Bool ? [:q, :l, :t] : [:q, :l])
        W₅ = [!isopt(r₅) || isdormant(r₅) || isequal(r₂, r₅)  ? 0. : relatedness(m, r₂, r₅, s) for r₅ ∈ R]
        r₅ = sample(rng, R, Weights(W₅))
        if !isopt(r₅) || isdormant(r₅) continue end
        d₂ = D[r₂.iᵈ]
        d₅ = D[r₅.iᵈ]
        i  = rand(rng, 1:r₂.n)
        k  = 1
        c₂ = C[r₂.iˢ]
        n₂ = c₂
        while true
            if isequal(k, i) n₂ = c₂ end
            if isequal(k, i) break end
            k += 1
            c₂ = C[c₂.iʰ]
        end
        n₁ = isequal(r₂.iˢ, n₂.iⁿ) ? D[n₂.iᵗ] : C[n₂.iᵗ]
        n₃ = isequal(r₂.iᵉ, n₂.iⁿ) ? D[n₂.iʰ] : C[n₂.iʰ]
        j  = rand(rng, 1:r₅.n)
        k  = 1
        c₅ = C[r₅.iˢ]
        n₅ = c₅
        while true
            if isequal(k, j) n₅ = c₅ end
            if isequal(k, j) break end
            k += 1
            c₅ = C[c₅.iʰ]
        end
        n₄ = isequal(r₅.iˢ, n₅.iⁿ) ? D[n₅.iᵗ] : C[n₅.iᵗ]
        n₆ = isequal(r₅.iᵉ, n₅.iⁿ) ? D[n₅.iʰ] : C[n₅.iʰ]
        # Step 2.2: Reconfigure
        # d₂ → ... → n₁ → n₅ → n₆ → ...  → d₂ and d₅ → ... → n₄ → n₂ → n₃ → ... → d₅
        c₂ = n₂
        tᵒ = n₁
        hᵒ = n₃
        tⁿ = n₄
        hⁿ = n₅
        while true
            removenode!(c₂, tᵒ, hᵒ, r₂, s)
            insertnode!(c₂, tⁿ, hⁿ, r₅, s)
            if isequal(hᵒ, d₂) break end
            tⁿ = c₂ 
            c₂ = C[hᵒ.iⁿ]
            hᵒ = isequal(r₂.iᵉ, c₂.iⁿ) ? D[c₂.iʰ] : C[c₂.iʰ]
        end
        c₅ = n₅
        tᵒ = c₂
        hᵒ = n₆
        tⁿ = n₁
        hⁿ = d₂
        while true
            removenode!(c₅, tᵒ, hᵒ, r₅, s)
            insertnode!(c₅, tⁿ, hⁿ, r₂, s)
            if isequal(hᵒ, d₅) break end
            tⁿ = c₅
            c₅ = C[hᵒ.iⁿ]
            hᵒ = isequal(r₅.iᵉ, c₅.iⁿ) ? D[c₅.iʰ] : C[c₅.iʰ]
        end
        # Step 2.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 2.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 2.5: Reconfigure back to the original state
        else
            # d₂ → ... → n₁ → n₂ → n₃ → ... → d₂ and d₅ → ... → n₄ → n₅ → n₆ → ... → d₅
            c₂ = n₅
            tᵒ = n₁
            hᵒ = isequal(r₂.iᵉ, c₂.iⁿ) ? D[c₂.iʰ] : C[c₂.iʰ]
            tⁿ = n₄
            hⁿ = n₂
            while true
                removenode!(c₂, tᵒ, hᵒ, r₂, s)
                insertnode!(c₂, tⁿ, hⁿ, r₅, s)
                if isequal(hᵒ, d₂) break end
                tⁿ = c₂ 
                c₂ = C[hᵒ.iⁿ]
                hᵒ = isequal(r₂.iᵉ, c₂.iⁿ) ? D[c₂.iʰ] : C[c₂.iʰ]
            end
            c₅ = n₂
            tᵒ = c₂
            hᵒ = isequal(r₅.iᵉ, c₅.iⁿ) ? D[c₅.iʰ] : C[c₅.iʰ]
            tⁿ = n₁
            hⁿ = d₂
            while true
                removenode!(c₅, tᵒ, hᵒ, r₅, s)
                insertnode!(c₅, tⁿ, hⁿ, r₂, s)
                if isequal(hᵒ, d₅) break end
                tⁿ = c₅
                c₅ = C[hᵒ.iⁿ]
                hᵒ = isequal(r₅.iᵉ, c₅.iⁿ) ? D[c₅.iʰ] : C[c₅.iʰ]
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
    W₁ = isopt.(D)
    # Step 2: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 2.1: Select a random depot pair
        d₁ = sample(rng, D, Weights(W₁))
        R₁ = [r₁ for v₁ ∈ d₁.V for r₁ ∈ v₁.R]
        if any(isdormant, R₁) continue end
        W₂ = [isequal(d₁, d₂) ? 0. : 1. for d₂ ∈ D]
        d₂ = sample(rng, D, Weights(W₂))
        R₂ = [r₂ for v₂ ∈ d₂.V for r₂ ∈ v₂.R]
        if any(isdormant, R₂) continue end
        if isequal(d₁, d₂) continue end
        if !isopt(d₁) && !isopt(d₂) continue end
        # Step 2.2: Swap vehicles, routes, and customer nodes
        I₁ = eachindex(d₁.V)
        I₂ = eachindex(d₂.V)
        while !isempty(d₁.V)
            v = d₁.V[1]
            movevehicle!(v, d₁, d₂, s)
        end
        for iᵛ ∈ I₂
            v = d₂.V[1]
            movevehicle!(v, d₂, d₁, s)
        end
        z′ = f(s)
        Δ  = z′ - z
        # Step 2.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 2.4: Reconfigure back to the original state
        else
            I₁ = eachindex(d₁.V)
            I₂ = eachindex(d₂.V)
            while !isempty(d₁.V)
                v = d₁.V[1]
                movevehicle!(v, d₁, d₂, s)
            end
            for iᵛ ∈ I₂
                v = d₂.V[1]
                movevehicle!(v, d₂, d₁, s)
            end
        end
    end
    postlocalsearch!(s)
    # Step 3: Return solution
    return s
end