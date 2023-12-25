using VRP
using Revise
using Random
using CPUTime
using DataFrames

let
    # Define instances
    instances = ["bar-n100-1"]
    # Define random number generators
    seeds = [1010, 1106, 1509, 1604, 1905, 2104, 2412, 2703, 2710, 2807]
    # Dataframes to store solution quality and run time
    df₁ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    df₂ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    for i ∈ eachindex(instances)
        instance = instances[i]
        # Visualize instance
        display(visualize(instance))
        for j ∈ eachindex(seeds)
            seed = seeds[j]
            println("\n instance: $instance | seed: $seed")
            rng = MersenneTwister(seed);
            # Define inital solution method and build the initial solution
            s₁ = initialize(rng, instance);
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
                            :relatedcustomer!   ,
                            :relatedroute!      ,
                            :relatedvehicle!    ,
                            :worstcustomer!     ,
                            :worstroute!        ,
                            :worstvehicle!      ,
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
                            :intraopt!
                        ]                       ,
                σ₁  =   15                      ,
                σ₂  =   10                      ,
                σ₃  =   3                       ,
                μ̲   =   0.01                    ,
                c̲   =   1                       ,
                μ̅   =   0.1                     ,
                c̅   =   15                      ,
                ω̅   =   0.05                    ,
                τ̅   =   0.5                     ,
                ω̲   =   0.01                    ,
                τ̲   =   0.01                    ,
                θ   =   0.9985                  ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t = @CPUelapsed s₂ = ALNS(rng, χ, s₁);
            # Visualize best solution
            display(visualize(s₂))
            # Fetch objective function values
            println("Objective function value:")
            println("   Initial: $(round(s₁.πᶠ + s₁.πᵒ, digits=3))")
            println("   Optimal: $(round(s₂.πᶠ + s₂.πᵒ, digits=3))")
            # Check if the solutions are feasible
            println("Solution feasibility:")
            println("   Initial: $(isfeasible(s₁)) | $(round(s₁.πᵖ, digits=3))")
            println("   Optimal: $(isfeasible(s₂)) | $(round(s₂.πᵖ, digits=3))")
            # Optimal solution characteristics
            println("Optimal solution characteristics:")
            nᵈ, nᵛ, nʳ = 0, 0, 0
            for d ∈ s₂.D nᵈ += VRP.isopt(d) end
            for d ∈ s₂.D for v ∈ d.V nᵛ += VRP.isopt(v) end end
            for d ∈ s₂.D for v ∈ d.V for r ∈ v.R nʳ += VRP.isopt(r) end end end
            println("   Number of depots: $nᵈ")
            println("   Number of vehicles: $nᵛ")
            println("   Number of routes: $nʳ")
            # Store Results
            df₁[i,j+1] = f(s₂)
            df₂[i,j+1] = t
            println(df₁)
            println(df₂)
        end
    end
    return
end