using VRP
using Revise
using Random
using CPUTime
using DataFrames

let
    # Set A
    A = ["lc101", "lc103", "lc109", "lr101", "lr104", "lr105", "lrc101", "lrc106", "lrc108"]
    # Set B
    B = ["lc201", "lc204", "lc207", "lr201", "lr202", "lr208", "lrc201", "lrc203", "lrc205"]
    # Define instances
    instances = [A..., B...]
    # Define random number generators
    seeds = [1010, 1106, 1509, 1604, 1905, 2104, 2412, 2703, 2710, 2807]
    # Dataframes to store solution quality and run time
    df₁ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    df₂ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    df₃ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    df₄ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    for i ∈ eachindex(instances)
        instance = instances[i]
        # Visualize instance
        display(visualize(instance))
        for j ∈ eachindex(seeds)
            seed = seeds[j]
            println("\n instance: $instance | seed: $seed")
            rng = MersenneTwister(seed);
            # Define inital solution method and build the initial solution
            s₁ = initialize(rng, instance)
            # Visualize initial solution
            display(visualize(s₁))
            # Define ALNS parameters
            x = max(100, lastindex(s₁.C))
            χ = ALNSparameters(
                j   =   50                      ,
                k   =   5                       ,
                n   =   x                       ,
                m   =   100x                    ,
                Ψᵣ  =   [
                            :randomcustomer!    ,
                            :randomroute!       ,
                            :randomvehicle!     ,
                            :randomdepot!       ,
                            :relatedcustomer!   ,
                            :relatedroute!      ,
                            :relatedvehicle!    ,
                            :relateddepot!      ,
                            :worstcustomer!     ,
                            :worstroute!        ,
                            :worstvehicle!      ,
                            :worstdepot!
                        ]                       ,
                Ψᵢ  =   [
                            :best!              ,
                            :precise!           ,
                            :perturb!           ,
                            :regret2!           ,
                            :regret3!
                        ]                       ,
                Ψₗ  =   [
                            :intramove!         ,
                            :intraswap!         ,
                            :intraopt!          ,
                            :intermove!         ,
                            :interswap!         ,
                            :interopt!
                        ]                       ,
                σ₁  =   15                      ,
                σ₂  =   10                      ,
                σ₃  =   3                       ,
                μ̲   =   0.05                    ,
                c̲   =   2                       ,
                μ̅   =   0.2                     ,
                c̅   =   30                      ,
                ω̅   =   0.05                    ,
                τ̅   =   0.5                     ,
                ω̲   =   0.01                    ,
                τ̲   =   0.01                    ,
                θ   =   0.9993                  ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t = @CPUelapsed s₂ = ALNS(rng, χ, s₁);
            # Visualize best solution
            display(visualize(s₂))
            # Optimal solution characteristics
            println("Optimal solution characteristics:")
            nᵈ₁, nᵛ₁ = 0, 0
            nᵈ₂, nᵛ₂ = 0, 0
            for d ∈ s₁.D nᵈ₁ += VRP.isopt(d) end
            for d ∈ s₂.D nᵈ₂ += VRP.isopt(d) end
            for d ∈ s₁.D for v ∈ d.V nᵛ₁ += VRP.isopt(v) end end
            for d ∈ s₂.D for v ∈ d.V nᵛ₂ += VRP.isopt(v) end end
            # Fetch objective function values
            println("Objective function value:")
            println("   Initial: $(round(s₁.πᶠ + s₁.πᵒ, digits=3))")
            println("   Optimal: $(round(s₂.πᶠ + s₂.πᵒ, digits=3))")
            # Fetch lexicograohic objective function values
            println("Lexicographic objective function value:")
            println("   Initial: $(nᵛ₁) | $(round(s₁.πᵒ, digits=3))")
            println("   Optimal: $(nᵛ₂) | $(round(s₂.πᵒ, digits=3))")
            # Check if the solutions are feasible
            println("Solution feasibility:")
            println("   Initial: $(isfeasible(s₁)) | $(round(s₁.πᵖ, digits=3))")
            println("   Optimal: $(isfeasible(s₂)) | $(round(s₂.πᵖ, digits=3))")
            # Store Results
            df₁[i,j+1] = s₂.πᶠ + s₂.πᵒ + s₂.πᵖ
            df₂[i,j+1] = s₂.πᵒ + s₂.πᵖ
            df₃[i,j+1] = nᵛ₂
            df₄[i,j+1] = t
            println(df₁)
            println(df₂)
            println(df₃)
            println(df₄)
        end
    end
    return
end 