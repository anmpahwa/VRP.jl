"""
    ALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution; mute=false))

Adaptive Large Neighborhood Search (ALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS adaptively searches large neighborhoods in the solution domain and
returns the best found solution. Additionally, displays a convergence plot.

Takes `mute` as argument. If `true` mutes progressbar and pltcnv output.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function ALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution; mute=false)
    # Step 0: Pre-initialize
    j, k = χ.j, χ.k
    n, m = χ.n, χ.m
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    μ̲, c̲ = χ.μ̲, χ.c̲
    μ̅, c̅ = χ.μ̅, χ.c̅
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲, τ̲ = χ.ω̲, χ.τ̲
    θ, ρ = χ.θ, χ.ρ   
    R = eachindex(Ψᵣ)
    I = eachindex(Ψᵢ)
    L = eachindex(Ψₗ)
    Z = OffsetVector{Float64}(undef, 0:j*(n+1))
    H = OffsetVector{UInt}(undef, 0:j*(n+1))
    # Step 1: Initialize
    s = deepcopy(sₒ)
    s⃰ = s
    z = f(sₒ)
    z⃰ = z
    h = hash(s)
    Z[0] = z
    H[0] = h
    t = ω̅ * z⃰/log(1/τ̅)
    Cᵣ, Pᵣ, Πᵣ, Wᵣ = zeros(Int, R), zeros(R), zeros(R), ones(R)
    Cᵢ, Pᵢ, Πᵢ, Wᵢ = zeros(Int, I), zeros(I), zeros(I), ones(I)
    # Step 2: Loop over segments.
    if !mute p = Progress(n * j, desc="Computing...", color=:blue, showspeed=true) end
    for u ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        for r ∈ R Cᵣ[r], Πᵣ[r] = 0, 0. end
        for i ∈ I Cᵢ[i], Πᵢ[i] = 0, 0. end
        # Step 2.2: Update selection probability for every removal and insertion operator
        for r ∈ R Pᵣ[r] = Wᵣ[r]/sum(values(Wᵣ)) end
        for i ∈ I Pᵢ[i] = Wᵢ[i]/sum(values(Wᵢ)) end
        # Step 2.3: Loop over iterations within the segment
        for v ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators.
            r = sample(rng, 1:length(Ψᵣ), Weights(Pᵣ))
            i = sample(rng, 1:length(Ψᵢ), Weights(Pᵢ))
            Cᵣ[r] += 1
            Cᵢ[i] += 1
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution.
            η = rand(rng)
            c = sum(isactive.(s.C))
            q = Int(floor(((1 - η) * min(c̲, μ̲ * c) + η * min(c̅, μ̅ * c))))
            s′= deepcopy(s)
            remove!(rng, q, s′, Ψᵣ[r])
            insert!(rng, s′, Ψᵢ[i])
            z′ = f(s′)
            h′ = hash(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁.
            if z′ < z⃰
                s = s′
                s⃰ = s′
                z = z′
                z⃰ = z′
                Πᵣ[r] += σ₁
                Πᵢ[i] += σ₂
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂.
            elseif z′ < z
                s = s′
                z = z′
                if h′ ∉ H
                    Πᵣ[r] += σ₂
                    Πᵢ[i] += σ₂
                end
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                if η < exp(-(z′ - z)/t)
                    s = s′
                    z = z′
                    if h′ ∉ H
                        Πᵣ[r] += σ₃
                        Πᵢ[i] += σ₃
                    end
                end
            end
            h = h′
            Z[(u - 1) * (n + 1) + v] = z
            H[(u - 1) * (n + 1) + v] = h
            t = max(t * θ, ω̲ * z⃰/log(1/τ̲))
            if !mute next!(p) end
        end
        # Step 2.4: Update weights for every removal and insertion operator.
        for r ∈ R if !iszero(Cᵣ[r]) Wᵣ[r] = ρ * Πᵣ[r] / Cᵣ[r] + (1 - ρ) * Wᵣ[r] end end
        for i ∈ I if !iszero(Cᵢ[i]) Wᵢ[i] = ρ * Πᵢ[i] / Cᵢ[i] + (1 - ρ) * Wᵢ[i] end end
        # Step 2.5: Reset current solution.
        if iszero(u % k) s, z = deepcopy(s⃰), z⃰ end
        # Step 2.6: Perform local search.
        s′ = deepcopy(s)
        for l ∈ L localsearch!(rng, m, s′, Ψₗ[l]) end
        z′ = f(s′)
        h′ = hash(s′)
        if z′ < z⃰
            s = s′
            s⃰ = s′
            z = z′
            z⃰ = z′
        elseif z′ < z
            s = s′
            z = z′
        end
        h = h′
        Z[u * (n + 1)] = z
        H[u * (n + 1)] = h
    end
    # Step 3: Display the convergence plot and return the best solution
    if !mute display(pltcnv(Z)) end
    return s⃰
end
ALNS(χ::ALNSparameters, s::Solution; mute=false) = ALNS(Random.GLOBAL_RNG, χ, s; mute=mute)