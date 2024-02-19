"""
    insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after inserting customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d = s.D[r.iᵈ]
    v = d.V[r.iᵛ]
    f = c.F[v.jᵛ]
    # update associated fuel station node, depot node, route, tail-head nodes, and the customer node
    ## fetch network features
    aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗᶜ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aᶜʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    φᵖᵈ = !isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)
    ## update cost
    s.πᶠ -= (isopt(f) ? f.πᶠ : 0.) + (isopt(d) ? d.πᶠ : 0.) + (isopt(v) ? v.πᶠ : 0.)
    s.πᵒ -= (f.πᵒ * c.ω) + (d.πᵒ * d.n) + (v.πᵈ * r.l)
    s.πᵖ -= (φᵖᵈ) * abs(c.qᶜ)
    ## update fuel station node
    f.ω += c.ω
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
    ## fetch network features
    φᵖᵈ = !isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)
    ## update cost
    s.πᶠ += (isopt(f) ? f.πᶠ : 0.) + (isopt(d) ? d.πᶠ : 0.) + (isopt(v) ? v.πᶠ : 0.)
    s.πᵒ += (f.πᵒ * c.ω) + (d.πᵒ * d.n) + (v.πᵈ * r.l)
    s.πᵖ += (φᵖᵈ) * abs(c.qᶜ)
    # update en-route parameters
    if isopt(r)
        ## initiate iterated parameters
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        φ  = true
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = nʰ
        while true
            ## fetch network features
            aʰᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            aʰʰ = s.A[(nʰ.iⁿ, nʰ.iʰ)]
            aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aᶠʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ]   
            cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ  = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            θᵒ  = θ - aᵗʰ.l/v.lᵛ
            θ̲   = min(aʰᶠ.l, aʰʰ.l)/v.lᵛ
            ## update costs
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= fᵗ.πᵒ * (φ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            ## update node characteristics
            if θ̲ < θᵒ
                ## directly visit the next node
                ω = 0.
                δ = 0.
                c.tᵃ = t + aᵗʰ.l/v.sᵛ
                c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q  = q
                c.θ̲  = θ̲
                c.θ  = θᵒ
                fᵗ.ω -= φ ? r.ω : nᵗ.ω
                r.l  -= φ ? r.δ : nᵗ.δ
                φ ? r.ω = ω : nᵗ.ω = ω
                φ ? r.δ = δ : nᵗ.δ = δ
                fᵗ.ω += φ ? r.ω : nᵗ.ω
                r.l  += φ ? r.δ : nᵗ.δ
            else
                ## pre-emptively re-fuel and then visit the next node
                ω = (1. - (θ - (aᵗᶠ.l/v.lᵛ))) * v.ωᵛ
                δ = aᵗᶠ.l + aᶠʰ.l - aᵗʰ.l
                c.tᵃ = t + aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠʰ.l/v.sᵛ
                c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q  = q
                c.θ̲  = θ̲
                c.θ  = 1. - aᶠʰ.l/v.lᵛ
                fᵗ.ω -= φ ? r.ω : nᵗ.ω
                r.l  -= φ ? r.δ : nᵗ.δ
                φ ? r.ω = ω : nᵗ.ω = ω
                φ ? r.δ = δ : nᵗ.δ = δ
                fᵗ.ω += φ ? r.ω : nᵗ.ω
                r.l  += φ ? r.δ : nᵗ.δ
            end
            ## fetch network features
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            ## update costs
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += fᵗ.πᵒ * (φ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            ## update iterated parameters
            φ  = false
            nᵗ = c
            nʰ = c.iʰ ≤ lastindex(s.D) ? s.D[c.iʰ] : s.C[c.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fʰ = nʰ.F[v.jᵛ]
            t  = c.tᵈ
            θ  = c.θ
            q -= c.qᶜ
            if isequal(c, cᵉ) break end
            c  = nʰ
        end
        ## fetch network features
        aʰᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aᶠʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        θᵒ  = θ - aᵗʰ.l/v.lᵛ
        θ̲   = aʰᶠ.l/v.lᵛ
        ## update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
        ## update node characteristics
        if θ̲ < θᵒ
            ## directly visit the next node
            ω = 0.
            δ = 0.
            r.tˢ = d.tˢ
            r.tᵉ = t + aᵗʰ.l/v.sᵛ
            r.θ̲  = θ̲
            r.θ  = θᵒ
            fᵗ.ω -= nᵗ.ω
            r.l  -= nᵗ.δ
            nᵗ.ω  = ω
            nᵗ.δ  = δ
            fᵗ.ω += nᵗ.ω
            r.l  += nᵗ.δ
        else
            ## pre-emptively re-fuel and then visit the next node
            ω = (1. - (θ - (aᵗᶠ.l/v.lᵛ))) * v.ωᵛ
            δ = aᵗᶠ.l + aᶠʰ.l - aᵗʰ.l
            r.tˢ = d.tˢ
            r.tᵉ = t + aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠʰ.l/v.sᵛ
            r.θ̲  = θ̲
            r.θ  = 1. - aᶠʰ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            r.l  -= nᵗ.δ
            nᵗ.ω  = ω
            nᵗ.δ  = δ
            fᵗ.ω += nᵗ.ω
            r.l  += nᵗ.δ
        end
        ## update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
    else
        ## fetch network features
        f   = d.F[v.jᵛ]
        aʰᶠ = s.A[(d.iⁿ, f.iⁿ)]
        θ̲   = aʰᶠ.l/v.lᵛ
        ## update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
        ## update route characteristics
        ω = 0.
        δ = 0.
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
        r.θ̲  = θ̲
        r.θ  = 1.
        f.ω -= r.ω
        r.l -= r.δ
        r.ω  = ω
        r.δ  = δ
        f.ω += r.ω
        r.l += r.δ
        ## update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
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
    cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    qᵒ  = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    φᵖᵈ = !isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)
    ## update cost
    s.πᶠ -= (isopt(f) ? f.πᶠ : 0.) + (isopt(d) ? d.πᶠ : 0.) + (isopt(v) ? v.πᶠ : 0.)
    s.πᵒ -= (f.πᵒ * c.ω) + (d.πᵒ * d.n) + (v.πᵈ * r.l)
    s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (φᵖᵈ) * abs(c.qᶜ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
    ## update fuel station node
    f.ω -= c.ω
    ## update depot node
    d.n -= 1
    ## update route
    if isdepot(nᵗ) r.iˢ = nʰ.iⁿ end
    if isdepot(nʰ) r.iᵉ = nᵗ.iⁿ end
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
    c.θ̲  = 1.
    c.θ  = 1.
    c.ω  = 0.
    c.δ  = 0.
    c.r  = NullRoute
    ## fetch network features
    qᵒ  = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    φᵖᵈ = !isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)
    ## update cost
    s.πᶠ += (isopt(f) ? f.πᶠ : 0.) + (isopt(d) ? d.πᶠ : 0.) + (isopt(v) ? v.πᶠ : 0.)
    s.πᵒ += (f.πᵒ * c.ω) + (d.πᵒ * d.n) + (v.πᵈ * r.l)
    s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (φᵖᵈ) * abs(c.qᶜ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
    # update en-route parameters
    if isopt(r)
        ## initiate iterated parameters
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        φ  = true
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = nʰ
        while true
            ## fetch network features
            aʰᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            aʰʰ = s.A[(nʰ.iⁿ, nʰ.iʰ)]
            aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aᶠʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            cᵖ  = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ]   
            cᵈ  = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ  = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            θᵒ  = θ - aᵗʰ.l/v.lᵛ
            θ̲   = min(aʰᶠ.l, aʰʰ.l)/v.lᵛ
            ## update costs
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= fᵗ.πᵒ * (φ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            ## update node characteristics
            if θ̲ < θᵒ
                ## directly visit the next node
                ω = 0.
                δ = 0.
                c.tᵃ = t + aᵗʰ.l/v.sᵛ
                c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q  = q
                c.θ̲  = θ̲
                c.θ  = θᵒ
                fᵗ.ω -= φ ? r.ω : nᵗ.ω
                r.l  -= φ ? r.δ : nᵗ.δ
                φ ? r.ω = ω : nᵗ.ω = ω
                φ ? r.δ = δ : nᵗ.δ = δ
                fᵗ.ω += φ ? r.ω : nᵗ.ω
                r.l  += φ ? r.δ : nᵗ.δ
            else
                ## pre-emptively re-fuel and then visit the next node
                ω = (1. - (θ - (aᵗᶠ.l/v.lᵛ))) * v.ωᵛ
                δ = aᵗᶠ.l + aᶠʰ.l - aᵗʰ.l
                c.tᵃ = t + aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠʰ.l/v.sᵛ
                c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q  = q
                c.θ̲  = θ̲
                c.θ  = 1. - aᶠʰ.l/v.lᵛ
                fᵗ.ω -= φ ? r.ω : nᵗ.ω
                r.l  -= φ ? r.δ : nᵗ.δ
                φ ? r.ω = ω : nᵗ.ω = ω
                φ ? r.δ = δ : nᵗ.δ = δ
                fᵗ.ω += φ ? r.ω : nᵗ.ω
                r.l  += φ ? r.δ : nᵗ.δ
            end
            ## fetch network features
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            ## update costs
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += fᵗ.πᵒ * (φ ? r.ω : nᵗ.ω) + v.πᵈ * r.l
            s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (c.θ̲ > c.θ) * (c.θ̲ - c.θ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            ## update iterated parameters
            φ  = false
            nᵗ = c
            nʰ = c.iʰ ≤ lastindex(s.D) ? s.D[c.iʰ] : s.C[c.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fʰ = nʰ.F[v.jᵛ]
            t  = c.tᵈ
            θ  = c.θ
            q -= c.qᶜ
            if isequal(c, cᵉ) break end
            c  = nʰ
        end
        ## fetch network features
        aʰᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        aᵗʰ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗᶠ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aᶠʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        θᵒ  = θ - aᵗʰ.l/v.lᵛ
        θ̲   = aʰᶠ.l/v.lᵛ
        ## update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
        ## update node characteristics
        if θ̲ < θᵒ
            ## directly visit the next node
            ω = 0.
            δ = 0.
            r.tˢ = d.tˢ
            r.tᵉ = t + aᵗʰ.l/v.sᵛ
            r.θ̲  = θ̲
            r.θ  = θᵒ
            fᵗ.ω -= nᵗ.ω
            r.l  -= nᵗ.δ
            nᵗ.ω  = ω
            nᵗ.δ  = δ
            fᵗ.ω += nᵗ.ω
            r.l  += nᵗ.δ
        else
            ## pre-emptively re-fuel and then visit the next node
            ω = (1. - (θ - (aᵗᶠ.l/v.lᵛ))) * v.ωᵛ
            δ = aᵗᶠ.l + aᶠʰ.l - aᵗʰ.l
            r.tˢ = d.tˢ
            r.tᵉ = t + aᵗᶠ.l/v.sᵛ + ω * fᵗ.τᵛ + aᶠʰ.l/v.sᵛ
            r.θ̲  = θ̲
            r.θ  = 1. - aᶠʰ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            r.l  -= nᵗ.δ
            nᵗ.ω  = ω
            nᵗ.δ  = δ
            fᵗ.ω += nᵗ.ω
            r.l  += nᵗ.δ
        end
        ## update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
    else
        ## fetch network features
        f   = d.F[v.jᵛ]
        aʰᶠ = s.A[(d.iⁿ, f.iⁿ)]
        θ̲   = aʰᶠ.l/v.lᵛ
        ## update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
        ## update route characteristics
        ω = 0.
        δ = 0.
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
        r.θ̲  = θ̲
        r.θ  = 1.
        f.ω -= r.ω
        r.l -= r.δ
        r.ω  = ω
        r.δ  = δ
        f.ω += r.ω
        r.l += r.δ
        ## update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ + r.l * v.πᵈ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (r.θ̲ > r.θ) * (r.θ̲ - r.θ)
    end
    return s
end