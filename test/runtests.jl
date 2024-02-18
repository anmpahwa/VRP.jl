using VRP
using Revise
using Test
using Random

let
    # ALNS parameters
    χ = ALNSparameters(
        j   =   50                      ,
        k   =   5                       ,
        n   =   5                       ,
        m   =   500                     ,
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
        θ   =   0.9965                  ,
        ρ   =   0.1
    );
    @testset "PDPTW" begin
        instances = ["lc101", "lc103", "lc104"]
        for instance ∈ instances
            visualize(instance)
            println(instance)
            s₁ = initialize(instance)
            s₂ = ALNS(χ, s₁)
            visualize(s₂)
            @test isfeasible(s₂)
            @test f(s₂) ≤ f(s₁)
        end
    end
    return
end
        
