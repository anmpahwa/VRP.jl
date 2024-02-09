"""
    visualize(instance; backend=gr)

Plots `instance`. Uses given `backend` to plot (defaults to `gr`).
"""
function visualize(instance; backend=gr)
    backend()
    D, C, _ = build(instance)
    fig= plot(legend=:none)
    K  = lastindex(C)
    X  = zeros(Float64, K)
    Y  = zeros(Float64, K)
    M₁ = fill("color", K)
    M₂ = zeros(Int, K)
    M₃ = fill(:shape, K)
    # Depot nodes
    for (k,d) ∈ pairs(D)
        X[k]  = d.x
        Y[k]  = d.y
        M₁[k] = "#b4464b"
        M₂[k] = 6
        M₃[k] = :rect
    end
    # Customer nodes
    for (k,c) ∈ pairs(C)
        X[k]  = c.x
        Y[k]  = c.y
        M₁[k] = isdelivery(c) ? "#d1e0ec" : "#ecddd1"
        M₂[k] = 5
        M₃[k] = :circle
    end
    scatter!(X, Y, color=M₁, markersize=M₂, markershape=M₃, markerstrokewidth=0)
    return fig
end
"""
    visualize(s::Solution; backend=gr)

Plots solution `s` depicting routes and unvisited nodes (if any).
Uses given `backend` to plot (defaults to `gr`).
"""
function visualize(s::Solution; backend=gr)
    backend()
    D = s.D
    C = s.C
    fig = plot(legend=:none)
    # Operational nodes: open depot nodes and closed customer nodes
    Z = vectorize(s)
    for Zᵈ ∈ Z
        for Zᵛ ∈ Zᵈ
            K  = eachindex(Zᵛ)
            X  = zeros(Float64, K)
            Y  = zeros(Float64, K)
            M₁ = fill("color", K)
            M₂ = zeros(Int, K)
            M₃ = fill(:shape, K)
            for k ∈ K
                i = Zᵛ[k]
                n = i ≤ lastindex(D) ? D[i] : C[i]
                X[k] = n.x
                Y[k] = n.y
                if isdepot(n) 
                    M₁[k] = "#82b446"
                    M₂[k] = 6
                    M₃[k] = :rect
                else 
                    M₁[k] = isdelivery(n) ? "#4682b4" : "#b47846"
                    M₂[k] = 5
                    M₃[k] = :circle
                end
            end
            scatter!(X, Y, color=M₁, markersize=M₂, markershape=M₃, markerstrokewidth=0)
            plot!(X, Y, color="#23415a")
        end
    end
    # Non-operational nodes: closed depot nodes and open customer nodes
    Z  = Int[] 
    for d ∈ D if !isopt(d) push!(Z, d.iⁿ) end end
    for c ∈ C if isopen(c) push!(Z, c.iⁿ) end end
    K  = eachindex(Z)
    X  = zeros(Float64, K)
    Y  = zeros(Float64, K)
    M₁ = fill("color", K)
    M₂ = zeros(Int, K)
    M₃ = fill(:shape, K)
    for k ∈ K
        i = Z[k]
        n = i ≤ length(D) ? D[i] : C[i]
        X[k] = n.x
        Y[k] = n.y
        if isdepot(n) 
            M₁[k] = "#b4464b"
            M₂[k] = 6
            M₃[k] = :rect
        else 
            M₁[k] = isdelivery(n) ? "#d1e0ec" : "#ecddd1"
            M₂[k] = 5
            M₃[k] = :circle
        end
    end
    scatter!(X, Y, color=M₁, markersize=M₂, markershape=M₃, markerstrokewidth=0)
     # Annotation
     x = min(minimum(getproperty.(C, :x)), minimum(getproperty.(D, :x)))
     y = max(maximum(getproperty.(C, :y)), maximum(getproperty.(D, :y)))
     annotate!(x, y, text("f(s): $(Int(round(f(s))))", :left, 10))
    return fig
end



"""
    animate(S::OffsetVector{Solution}; fps=10)

Iteratively plots solutions in `S` to develop a gif at given `fps`.
"""
function animate(S::OffsetVector{Solution}; fps=1)
    s⃰ = S[0]
    z⃰ = f(s⃰)
    figs = []
    for k ∈ eachindex(S)
        s = S[k]
        z = f(s)
        if z < 0.99z⃰ 
            s⃰ = s
            z⃰ = z
            fig = visualize(s⃰, backend=gr)
            plot!(title="Iteration #$k", titlefontsize=11)
            push!(figs, fig)
        end
    end
    anim = @animate for fig ∈ figs
        plot(fig)
    end
    gif(anim, fps=fps, show_msg=false)
end



"""
    pltcnv(Z::OffsetVector{Float64}; backend=gr)

Plots convergence using objective function evaluations vector `Z`. 
Uses given `backend` to plot (defaults to `gr`).
"""
function pltcnv(Z::OffsetVector{Float64}; backend=gr)
    backend()
    fig= plot(legend=:none)
    Y₁ = Int[]
    z⃰  = Z[0]
    for (k, z) ∈ pairs(Z)
        if z < 0.99z⃰ 
            z⃰ = z
            push!(Y₁, k)
        end
    end
    vline!(Y₁, color=:black, linewidth=0.25)
    Y₂ = zeros(eachindex(Z))
    z⃰  = minimum(Z)
    for (k, z) ∈ pairs(Z)
        Y₂[k] = (z/z⃰ - 1) * 100 
    end
    X = eachindex(Z)
    plot!(X, Y₂, xlabel="iterations", ylabel="deviation from the best (%)", color=:steelblue)
    return fig
end