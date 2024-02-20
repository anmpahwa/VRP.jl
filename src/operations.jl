"""
    insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after inserting customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d = s.D[r.iᵈ]
    v = d.V[r.iᵛ]
    # update associated depot node, route, tail-head nodes, and the customer node
    ## fetch network features
    aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗᶜ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aᶜʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    ## update cost
    s.πᶠ -= isopt(d) * d.πᶠ + isopt(v) * v.πᶠ
    s.πᵒ -= d.πᵒ * d.n + v.πᵈ * r.l
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) ? abs(c.qᶜ) : 0.
    ## update depot node
    d.n += 1
    ## update route
    if isdepot(nᵗ) r.iˢ = c.iⁿ end
    if isdepot(nʰ) r.iᵉ = c.iⁿ end
    r.n += 1
    r.l += aᵗᶜ.l + c.δ + aᶜʰ.l - aᵗʰ.l
    ## update tail-head nodes
    if iscustomer(nᵗ) nᵗ.iʰ = c.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = c.iⁿ end
    ## update the customer node
    c.iᵗ = nᵗ.iⁿ
    c.iʰ = nʰ.iⁿ
    c.r  = r
    ## update cost
    s.πᶠ += isopt(d) * d.πᶠ + isopt(v) * v.πᶠ
    s.πᵒ += d.πᵒ * d.n + v.πᵈ * r.l
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) ? abs(c.qᶜ) : 0.
    # update en-route parameters
    if isopt(r)
        ## initiate iterated parameters
        φᵗ = true
        nᵗ = d
        nᵒ = s.C[r.iˢ]
        fᵗ = nᵒ.F[v.jᵛ]
        fᵒ = nᵒ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        while true
            ## fetch network features
            aᵒᶠ = s.A[(nᵒ.iⁿ, fᵒ.iⁿ)]
            aᵒʰ = s.A[(nᵒ.iⁿ, nᵒ.iʰ)]
            aᵗᵒ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
            aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aᶠᵒ = s.A[(fᵗ.iⁿ, nᵒ.iⁿ)]
            nᵖ  = isdelivery(nᵒ) ? s.C[nᵒ.jⁿ] : s.C[nᵒ.iⁿ]   
            nᵈ  = isdelivery(nᵒ) ? s.C[nᵒ.iⁿ] : s.C[nᵒ.jⁿ]
            ## update costs
            s.πᶠ -= isopt(fᵗ) * fᵗ.πᶠ
            s.πᵒ -= fᵗ.πᵒ * (φᵗ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ -= (nᵒ.tᵃ > nᵒ.tˡ) ? (nᵒ.tᵃ - nᵒ.tˡ) : 0.
            s.πᵖ -= (nᵖ.tᵃ > nᵈ.tᵃ) ? (nᵖ.tᵃ - nᵈ.tᵃ) : 0.
            s.πᵖ -= (nᵒ.q > v.qᵛ) ? (nᵒ.q - v.qᵛ) : 0.
            s.πᵖ -= (nᵒ.θ̲ > nᵒ.θ) ? (nᵒ.θ̲ - nᵒ.θ) : 0.
            ## check for re-fueling
            θ̲ᵒ = min(aᵒᶠ.l, aᵒʰ.l)/v.lᵛ
            θᵒ = θ - aᵗᵒ.l/v.lᵛ
            θᶠ = θ - aᵗᶠ.l/v.lᵛ
            φᶠ = (θ̲ᵒ > θᵒ) && (θᶠ ≥ 0.)
            ## update network features
            ω = φᶠ ? ((1. - θᶠ) * v.ωᵛ) : (0.)
            δ = φᶠ ? (aᵗᶠ.l + aᶠᵒ.l - aᵗᵒ.l) : (0.)
            nᵒ.tᵃ = t + (φᶠ ? (aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠᵒ.l/v.sᵛ) : (aᵗᵒ.l/v.sᵛ))
            nᵒ.tᵈ = nᵒ.tᵃ + v.τᶜ + max(0., nᵒ.tᵉ - nᵒ.tᵃ - v.τᶜ) + nᵒ.τᶜ
            nᵒ.q  = q + (isdelivery(nᵒ) ? 0. : abs(nᵒ.qᶜ))
            nᵒ.θ̲  = θ̲ᵒ
            nᵒ.θ  = φᶠ ? (1. - aᶠᵒ.l/v.lᵛ) : (θᵒ)
            fᵗ.ω -= φᵗ ? r.ω : nᵗ.ω
            r.l  -= φᵗ ? r.δ : nᵗ.δ
            φᵗ ? r.ω = ω : nᵗ.ω = ω
            φᵗ ? r.δ = δ : nᵗ.δ = δ
            fᵗ.ω += φᵗ ? r.ω : nᵗ.ω
            r.l  += φᵗ ? r.δ : nᵗ.δ
            ## update costs
            s.πᶠ += isopt(fᵗ) * fᵗ.πᶠ
            s.πᵒ += fᵗ.πᵒ * (φᵗ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ += (nᵒ.tᵃ > nᵒ.tˡ) ? (nᵒ.tᵃ - nᵒ.tˡ) : 0.
            s.πᵖ += (nᵖ.tᵃ > nᵈ.tᵃ) ? (nᵖ.tᵃ - nᵈ.tᵃ) : 0.
            s.πᵖ += (nᵒ.q > v.qᵛ) ? (nᵒ.q - v.qᵛ) : 0.
            s.πᵖ += (nᵒ.θ̲ > nᵒ.θ) ? (nᵒ.θ̲ - nᵒ.θ) : 0.
            ## update iterated parameters
            φᵗ = false
            nᵗ = nᵒ
            nᵒ = nᵗ.iʰ ≤ lastindex(s.D) ? s.D[nᵗ.iʰ] : s.C[nᵗ.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fᵒ = nᵒ.F[v.jᵛ]
            t  = nᵗ.tᵈ
            θ  = nᵗ.θ
            q  = nᵗ.q - (isdelivery(nᵗ) ? abs(nᵗ.qᶜ) : 0.)
            if isequal(nᵒ, d) break end
        end
        ## fetch network features
        aᵗᵒ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
        aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aᶠᵒ = s.A[(fᵗ.iⁿ, nᵒ.iⁿ)]
        ## update costs
        s.πᶠ -= isopt(fᵗ) * fᵗ.πᶠ
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ -= (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ -= (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ -= ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
        ## check for re-fueling
        θ̲ᵒ = 0.
        θᵒ = θ - aᵗᵒ.l/v.lᵛ
        θᶠ = θ - aᵗᶠ.l/v.lᵛ
        φᶠ = (θ̲ᵒ > θᵒ) && (θᶠ ≥ 0.)
        ## update network features
        ω = φᶠ ? ((1. - θᶠ) * v.ωᵛ) : (0.)
        δ = φᶠ ? (aᵗᶠ.l + aᶠᵒ.l - aᵗᵒ.l) : (0.)
        r.tˢ  = d.tˢ
        r.tᵉ  = t + (φᶠ ? (aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠᵒ.l/v.sᵛ) : (aᵗᵒ.l/v.sᵛ))
        r.θ̲   = θ̲ᵒ
        r.θ   = φᶠ ? (1. - aᶠᵒ.l/v.lᵛ) : (θᵒ)
        fᵗ.ω -= nᵗ.ω
        r.l  -= nᵗ.δ
        nᵗ.ω  = ω
        nᵗ.δ  = δ
        fᵗ.ω += nᵗ.ω
        r.l  += nᵗ.δ
        ## update costs
        s.πᶠ += isopt(fᵗ) * fᵗ.πᶠ
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ += (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ += (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ += ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
    end
    return s
end
"""
    removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after removing customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d = s.D[r.iᵈ]
    v = d.V[r.iᵛ]
    f = c.F[v.jᵛ]
    # update associated fuel station node, depot node, route, tail-head nodes, and the customer node
    ## fetch network features
    aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗᶜ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aᶜʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cˢ  = s.C[r.iˢ]
    cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    ## update cost
    s.πᶠ -= isopt(f) * f.πᶠ + isopt(d) * d.πᶠ + isopt(v) * v.πᶠ
    s.πᵒ -= f.πᵒ * c.ω + d.πᵒ * d.n + v.πᵈ * r.l
    s.πᵖ -= (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) ? abs(c.qᶜ) : 0.
    s.πᵖ -= (cᵖ.tᵃ > cᵈ.tᵃ) ? (cᵖ.tᵃ - cᵈ.tᵃ) : 0.
    s.πᵖ -= (c.q > v.qᵛ) ? (c.q - v.qᵛ) : 0.
    s.πᵖ -= (c.θ̲ > c.θ) ? (c.θ̲ - c.θ) : 0.
    ## update fuel station node
    f.ω -= isequal(c, cˢ) * r.ω + c.ω
    ## update depot node
    d.n -= 1
    ## update route
    if isdepot(nᵗ) r.iˢ = nʰ.iⁿ end
    if isdepot(nʰ) r.iᵉ = nᵗ.iⁿ end
    r.ω -= isequal(c, cˢ) * r.ω
    r.n -= 1
    r.l -= aᵗᶜ.l + c.δ + aᶜʰ.l - aᵗʰ.l
    ## update tail-head nodes
    if iscustomer(nᵗ) nᵗ.iʰ = nʰ.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = nᵗ.iⁿ end
    ## update the customer node
    c.iᵗ = 0
    c.iʰ = 0
    c.tᵃ = isdelivery(c) ? c.tˡ : c.tᵉ
    c.tᵈ = c.tᵃ + c.τᶜ
    c.q  = 0.
    c.θ̲  = 0.
    c.θ  = 1.
    c.ω  = 0.
    c.δ  = 0.
    c.r  = NullRoute
    ## update cost
    s.πᶠ += isopt(f) * f.πᶠ + isopt(d) * d.πᶠ + isopt(v) * v.πᶠ
    s.πᵒ += f.πᵒ * c.ω + d.πᵒ * d.n + v.πᵈ * r.l
    s.πᵖ += (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) ? abs(c.qᶜ) : 0.
    s.πᵖ += (cᵖ.tᵃ > cᵈ.tᵃ) ? (cᵖ.tᵃ - cᵈ.tᵃ) : 0.
    s.πᵖ += (c.q > v.qᵛ) ? (c.q - v.qᵛ) : 0.
    s.πᵖ += (c.θ̲ > c.θ) ? (c.θ̲ - c.θ) : 0.
    # update en-route parameters
    if isopt(r)
        ## initiate iterated parameters
        φᵗ = true
        nᵗ = d
        nᵒ = s.C[r.iˢ]
        fᵗ = nᵒ.F[v.jᵛ]
        fᵒ = nᵒ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        while true
            ## fetch network features
            aᵒᶠ = s.A[(nᵒ.iⁿ, fᵒ.iⁿ)]
            aᵒʰ = s.A[(nᵒ.iⁿ, nᵒ.iʰ)]
            aᵗᵒ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
            aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aᶠᵒ = s.A[(fᵗ.iⁿ, nᵒ.iⁿ)]
            nᵖ  = isdelivery(nᵒ) ? s.C[nᵒ.jⁿ] : s.C[nᵒ.iⁿ]   
            nᵈ  = isdelivery(nᵒ) ? s.C[nᵒ.iⁿ] : s.C[nᵒ.jⁿ]
            ## update costs
            s.πᶠ -= isopt(fᵗ) * fᵗ.πᶠ
            s.πᵒ -= fᵗ.πᵒ * (φᵗ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ -= (nᵒ.tᵃ > nᵒ.tˡ) ? (nᵒ.tᵃ - nᵒ.tˡ) : 0.
            s.πᵖ -= (nᵖ.tᵃ > nᵈ.tᵃ) ? (nᵖ.tᵃ - nᵈ.tᵃ) : 0.
            s.πᵖ -= (nᵒ.q > v.qᵛ) ? (nᵒ.q - v.qᵛ) : 0.
            s.πᵖ -= (nᵒ.θ̲ > nᵒ.θ) ? (nᵒ.θ̲ - nᵒ.θ) : 0.
            ## check for re-fueling
            θ̲ᵒ = min(aᵒᶠ.l, aᵒʰ.l)/v.lᵛ
            θᵒ = θ - aᵗᵒ.l/v.lᵛ
            θᶠ = θ - aᵗᶠ.l/v.lᵛ
            φᶠ = (θ̲ᵒ > θᵒ) && (θᶠ ≥ 0.)
            ## update network features
            ω = φᶠ ? ((1. - θᶠ) * v.ωᵛ) : (0.)
            δ = φᶠ ? (aᵗᶠ.l + aᶠᵒ.l - aᵗᵒ.l) : (0.)
            nᵒ.tᵃ = t + (φᶠ ? (aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠᵒ.l/v.sᵛ) : (aᵗᵒ.l/v.sᵛ))
            nᵒ.tᵈ = nᵒ.tᵃ + v.τᶜ + max(0., nᵒ.tᵉ - nᵒ.tᵃ - v.τᶜ) + nᵒ.τᶜ
            nᵒ.q  = q + (isdelivery(nᵒ) ? 0. : abs(nᵒ.qᶜ))
            nᵒ.θ̲  = θ̲ᵒ
            nᵒ.θ  = φᶠ ? (1. - aᶠᵒ.l/v.lᵛ) : (θᵒ)
            fᵗ.ω -= φᵗ ? r.ω : nᵗ.ω
            r.l  -= φᵗ ? r.δ : nᵗ.δ
            φᵗ ? r.ω = ω : nᵗ.ω = ω
            φᵗ ? r.δ = δ : nᵗ.δ = δ
            fᵗ.ω += φᵗ ? r.ω : nᵗ.ω
            r.l  += φᵗ ? r.δ : nᵗ.δ
            ## update costs
            s.πᶠ += isopt(fᵗ) * fᵗ.πᶠ
            s.πᵒ += fᵗ.πᵒ * (φᵗ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ += (nᵒ.tᵃ > nᵒ.tˡ) ? (nᵒ.tᵃ - nᵒ.tˡ) : 0.
            s.πᵖ += (nᵖ.tᵃ > nᵈ.tᵃ) ? (nᵖ.tᵃ - nᵈ.tᵃ) : 0.
            s.πᵖ += (nᵒ.q > v.qᵛ) ? (nᵒ.q - v.qᵛ) : 0.
            s.πᵖ += (nᵒ.θ̲ > nᵒ.θ) ? (nᵒ.θ̲ - nᵒ.θ) : 0.
            ## update iterated parameters
            φᵗ = false
            nᵗ = nᵒ
            nᵒ = nᵗ.iʰ ≤ lastindex(s.D) ? s.D[nᵗ.iʰ] : s.C[nᵗ.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fᵒ = nᵒ.F[v.jᵛ]
            t  = nᵗ.tᵈ
            θ  = nᵗ.θ
            q  = nᵗ.q - (isdelivery(nᵗ) ? abs(nᵗ.qᶜ) : 0.)
            if isequal(nᵒ, d) break end
        end
        ## fetch network features
        aᵗᵒ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
        aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aᶠᵒ = s.A[(fᵗ.iⁿ, nᵒ.iⁿ)]
        ## update costs
        s.πᶠ -= isopt(fᵗ) * fᵗ.πᶠ
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ -= (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ -= (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ -= ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
        ## check for re-fueling
        θ̲ᵒ = 0.
        θᵒ = θ - aᵗᵒ.l/v.lᵛ
        θᶠ = θ - aᵗᶠ.l/v.lᵛ
        φᶠ = (θ̲ᵒ > θᵒ) && (θᶠ ≥ 0.)
        ## update network features
        ω = φᶠ ? ((1. - θᶠ) * v.ωᵛ) : (0.)
        δ = φᶠ ? (aᵗᶠ.l + aᶠᵒ.l - aᵗᵒ.l) : (0.)
        r.tˢ  = d.tˢ
        r.tᵉ  = t + (φᶠ ? (aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠᵒ.l/v.sᵛ) : (aᵗᵒ.l/v.sᵛ))
        r.θ̲   = θ̲ᵒ
        r.θ   = φᶠ ? (1. - aᶠᵒ.l/v.lᵛ) : (θᵒ)
        fᵗ.ω -= nᵗ.ω
        r.l  -= nᵗ.δ
        nᵗ.ω  = ω
        nᵗ.δ  = δ
        fᵗ.ω += nᵗ.ω
        r.l  += nᵗ.δ
        ## update costs
        s.πᶠ += isopt(fᵗ) * fᵗ.πᶠ
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ += (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ += (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ += ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
    else
        ## update costs
        s.πᶠ -= isopt(f) * f.πᶠ
        s.πᵒ -= r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ -= (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ -= (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ -= ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
        ## update network features
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
        r.θ̲  = 0.
        r.θ  = 1.
        f.ω -= r.ω
        r.l -= r.δ
        r.ω  = 0.
        r.δ  = 0.
        f.ω += r.ω
        r.l += r.δ
        ## update costs
        s.πᶠ += isopt(f) * f.πᶠ
        s.πᵒ += r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (r.θ̲ > r.θ) ? (r.θ̲ - r.θ) : 0.
        s.πᵖ += (d.tˢ > r.tˢ) ? (d.tˢ - r.tˢ) : 0.
        s.πᵖ += (r.tᵉ > d.tᵉ) ? (r.tᵉ - d.tᵉ) : 0.
        s.πᵖ += ((r.tᵉ - r.tˢ) > v.τʷ) ? ((r.tᵉ - r.tˢ) - v.τʷ) : 0.
    end
    return s
end