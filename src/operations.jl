"""
    insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after inserting customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    # update associated customer nodes
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ)
    if iscustomer(nᵗ) nᵗ.iʰ = c.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = c.iⁿ end
    c.iᵗ  = nᵗ.iⁿ
    c.iʰ  = nʰ.iⁿ
    c.r   = r
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ)
    # update associated vehicle-route
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= 0.
    r.x   = (r.n * r.x + c.x)/(r.n + 1)
    r.y   = (r.n * r.y + c.y)/(r.n + 1)
    if isdepot(nᵗ) r.iˢ = c.iⁿ end
    if isdepot(nʰ) r.iᵉ = c.iⁿ end
    r.n  += 1
    r.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += 0.
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.n * d.πᵒ
    s.πᵖ -= 0.
    d.n  += 1
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.n * d.πᵒ
    s.πᵖ += 0.
    # update en-route parameters
    if isopt(r)
        # initiate iterated parameters
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = nʰ
        while true
            # fetch network features
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ]   
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            aᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            aʰ = s.A[(nʰ.iⁿ, nʰ.iʰ)]
            θˡ = max(aᶠ.l, aʰ.l)/v.lᵛ
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            θᵒ = θ - aᵒ.l/v.lᵛ
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            # update costs
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= isdepot(nᵗ) ? r.ω * fᵗ.πᵒ : nᵗ.ω * fᵗ.πᵒ
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update node characteristics
            if θˡ < θᵒ
                # directly visit the next node
                ω     = 0.
                c.tᵃ  = t + aᵒ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - c.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = θ - aᵒ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            else
                # pre-emptively re-fuel and then visit the next node
                ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
                c.tᵃ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = 1. - aʰ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            end
            # fetch network features
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            # update costs
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += isdepot(nᵗ) ? r.ω * fᵗ.πᵒ : nᵗ.ω * fᵗ.πᵒ
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update iterated parameters
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
        # fetch network features
        aᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        θˡ = aᶠ.l/v.lᵛ
        aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        θᵒ = θ - aᵒ.l/v.lᵛ
        # update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
        # update node characteristics
        if θˡ < θᵒ
            # directly visit the next node
            ω     = 0.
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵒ.l/v.sᵛ
            r.θ   = θ - aᵒ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            nᵗ.ω  = ω
            fᵗ.ω += nᵗ.ω
            r.l  += 0.
        else
            # pre-emptively re-fuel and then visit the next node
            ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
            r.θ   = 1. - aʰ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            nᵗ.ω  = ω
            fᵗ.ω += nᵗ.ω
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
        end
        # update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    else
        f = d.F[v.jᵛ]
        # update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
        # update route characteristics
        r.tˢ  = d.tˢ
        r.tᵉ  = r.tˢ
        r.θ   = 1.
        f.ω  -= r.ω
        r.ω   = 0.
        f.ω  += r.ω
        # update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    end
    return s
end
"""
    removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after removing customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    f  = c.F[v.jᵛ]
    θˡ = s.A[(c.iⁿ, f.iⁿ)].l/v.lᵛ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    # update associated customer and fuel station nodes
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ -= isopt(f) ? f.πᶠ : 0.
    s.πᵒ -= c.ω * f.πᵒ
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
    if iscustomer(nᵗ) nᵗ.iʰ = nʰ.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = nᵗ.iⁿ end
    c.iᵗ  = 0
    c.iʰ  = 0
    c.tᵃ  = isdelivery(c) ? c.tˡ : c.tᵉ
    c.tᵈ  = c.tᵃ + c.τᶜ
    c.q   = 0.
    c.θ   = 1.
    f.ω  -= c.ω
    c.ω   = 0.
    f.ω  += c.ω
    c.r   = NullRoute
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ += isopt(f) ? f.πᶠ : 0.
    s.πᵒ += isopt(f) ? f.πᶠ : 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
    # update associated vehicle-route
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= 0.
    r.x   = isone(r.n) ? 0. : (r.n * r.x - c.x)/(r.n - 1)
    r.y   = isone(r.n) ? 0. : (r.n * r.y - c.y)/(r.n - 1)
    if isdepot(nᵗ) r.iˢ = nʰ.iⁿ end
    if isdepot(nʰ) r.iᵉ = nᵗ.iⁿ end
    r.n  -= 1
    r.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += 0.
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.n * d.πᵒ
    s.πᵖ -= 0.
    d.n  -= 1
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.n * d.πᵒ
    s.πᵖ += 0.
    # update en-route parameters
    if isopt(r)
        # initiate iterated parameters
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = nʰ
        while true
            # fetch network features
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ]   
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            aᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            aʰ = s.A[(nʰ.iⁿ, nʰ.iʰ)]
            θˡ = max(aᶠ.l, aʰ.l)/v.lᵛ
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            θᵒ = θ - aᵒ.l/v.lᵛ
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            # update costs
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= isdepot(nᵗ) ? r.ω * fᵗ.πᵒ : nᵗ.ω * fᵗ.πᵒ
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update node characteristics
            if θˡ < θᵒ
                # directly visit the next node
                ω     = 0.
                c.tᵃ  = t + aᵒ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - c.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = θ - aᵒ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            else
                # pre-emptively re-fuel and then visit the next node
                ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
                c.tᵃ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = 1. - aʰ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            end
            # fetch network features
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            # update costs
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += isdepot(nᵗ) ? r.ω * fᵗ.πᵒ : nᵗ.ω * fᵗ.πᵒ
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update iterated parameters
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
        # fetch network features
        aᶠ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        θˡ = aᶠ.l/v.lᵛ
        aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        θᵒ = θ - aᵒ.l/v.lᵛ
        # update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
        # update node characteristics
        if θˡ < θᵒ
            # directly visit the next node
            ω     = 0.
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵒ.l/v.sᵛ
            r.θ   = θ - aᵒ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            nᵗ.ω  = ω
            fᵗ.ω += nᵗ.ω
            r.l  += 0.
        else
            # pre-emptively re-fuel and then visit the next node
            ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
            r.θ   = 1. - aʰ.l/v.lᵛ
            fᵗ.ω -= nᵗ.ω
            nᵗ.ω  = ω
            fᵗ.ω += nᵗ.ω
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
        end
        # update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += nᵗ.ω * fᵗ.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    else
        f = d.F[v.jᵛ]
        # update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
        # update route characteristics
        r.tˢ  = d.tˢ
        r.tᵉ  = r.tˢ
        r.θ   = 1.
        f.ω  -= r.ω
        r.ω   = 0.
        f.ω  += r.ω
        # update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ + (r.tᵉ - r.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    end
    return s
end