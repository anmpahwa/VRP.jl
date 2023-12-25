using VRP
using Revise
using Test
using Random

let
    # ALNS parameters
    χ = ALNSparameters(
        j   =   50                      ,
        k   =   5                       ,
        n   =   10                      ,
        m   =   1000                    ,
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
                ]                    ,
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
                    :interopt!          ,
                    :swapdepot!
                ]                       ,
        σ₁  =   15                      ,
        σ₂  =   10                      ,
        σ₃  =   3                       ,
        μ̲   =   0.1                     ,
        c̲   =   4                       ,
        μ̅   =   0.4                     ,
        c̅   =   60                      ,
        ω̅   =   0.05                    ,
        τ̅   =   0.5                     ,
        ω̲   =   0.01                    ,
        τ̲   =   0.01                    ,
        θ   =   0.9985                  ,
        ρ   =   0.1
    );
    @testset "VRP" begin
        instances = ["bar-n100-1"]
        for instance ∈ instances
            visualize(instance)
            println(instance)
            sₒ = initialize(instance)
            s⃰  = ALNS(χ, sₒ)
            visualize(s⃰)
            @test isfeasible(s⃰)
            @test f(s⃰) ≤ f(sₒ)
        end
    end
    return
end
        
