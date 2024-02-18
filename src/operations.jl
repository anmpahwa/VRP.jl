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
    f     = d.F[v.jᵛ]
    θˡ    = s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        # initiate iterated parameters
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = cˢ
        while true
            # fetch network features
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            aˡ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            θˡ = aˡ.l/v.lᵛ
            θᵒ = θ - aᵒ.l/v.lᵛ
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update node characteristics
            if θˡ < θᵒ
                # directly visit the customer node
                ω     = 0.
                c.tᵃ  = t + aᵒ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = θ - aᵒ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            else
                # pre-emptively re-fuel and then visit the customer node
                ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
                c.tᵃ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = 1. - aʰ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += aᵗ.l + aʰ.l - aᵒ.l
            end
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
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
            c  = s.C[c.iʰ]
        end
        # fetch network features
        aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        a′ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        θˡ = a′.l/v.lᵛ
        θᵒ = θ - aᵒ.l/v.lᵛ
        # update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
        s.πᵖ -= 0.
        # update node characteristics
        if θˡ < θᵒ
            # directly visit the depot node
            ω     = 0.
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵒ.l/v.sᵛ
            r.θ   = θ - aᵒ.l/v.lᵛ
            fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
            isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
            fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
            r.l  += 0.
        else
            # pre-emptively re-fuel and then visit the depot node
            ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
            r.θ   = 1. - aʰ.l/v.lᵛ
            fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
            isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
            fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
        end
        # update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
        s.πᵖ += 0.
    else
        f = d.F[v.jᵛ]
        # update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ
        s.πᵖ -= 0.
        # update route characteristics
        r.tˢ  = d.tˢ
        r.tᵉ  = r.tˢ
        f.ω  -= r.ω
        r.ω   = 0.
        f.ω  += r.ω
        # update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ
        s.πᵖ += 0.
    end
    f     = d.F[v.jᵛ]
    θˡ    = s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
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
    # update associated customer nodes
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
    if iscustomer(nᵗ) nᵗ.iʰ = nʰ.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = nᵗ.iⁿ end
    c.iᵗ  = 0
    c.iʰ  = 0
    c.tᵃ  = isdelivery(c) ? c.tˡ : c.tᵉ
    c.tᵈ  = c.tᵃ + c.τᶜ
    c.θ   = 1.
    c.q   = 0.
    c.r   = NullRoute
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ += 0.
    s.πᵒ += 0.
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
    # update associated fuel station node
    s.πᶠ -= isopt(f) ? f.πᶠ : 0.
    s.πᵒ -= c.ω * f.πᵒ
    s.πᵖ -= 0.
    f.ω  -= c.ω
    c.ω   = 0.
    f.ω  += c.ω
    s.πᶠ += isopt(f) ? f.πᶠ : 0.
    s.πᵒ += c.ω * f.πᵒ
    s.πᵖ += 0.
    # update en-route parameters
    f     = d.F[v.jᵛ]
    θˡ    = s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        # initiate iterated parameters
        nᵗ = d
        nʰ = cˢ
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        q  = 0.
        θ  = 1.
        c  = cˢ
        while true
            # fetch network features
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            aˡ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            θˡ = aˡ.l/v.lᵛ
            θᵒ = θ - aᵒ.l/v.lᵛ
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (θˡ > c.θ) * (θˡ - c.θ)
            # update node characteristics
            if θˡ < θᵒ
                # directly visit the customer node
                ω     = 0.
                c.tᵃ  = t + aᵒ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = θ - aᵒ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += 0.
            else
                # pre-emptively re-fuel and then visit the customer node
                ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
                c.tᵃ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                c.q   = q
                c.θ   = 1. - aʰ.l/v.lᵛ
                fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
                isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
                fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
                r.l  += aᵗ.l + aʰ.l - aᵒ.l
            end
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
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
            c  = s.C[c.iʰ]
        end
        # fetch network features
        aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        a′ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        θˡ = a′.l/v.lᵛ
        θᵒ = θ - aᵒ.l/v.lᵛ
        # update costs
        s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ -= (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
        s.πᵖ -= 0.
        # update node characteristics
        if θˡ < θᵒ
            # directly visit the depot node
            ω     = 0.
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵒ.l/v.sᵛ
            r.θ   = θ - aᵒ.l/v.lᵛ
            fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
            isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
            fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
            r.l  += 0.
        else
            # pre-emptively re-fuel and then visit the depot node
            ω     = (1. - (θ - (aᵗ.l/v.lᵛ))) * v.ωᵛ
            r.tˢ  = d.tˢ
            r.tᵉ  = t + aᵗ.l/v.sᵛ + ω * fᵗ.τᵛ + aʰ.l/v.sᵛ
            r.θ   = 1. - aʰ.l/v.lᵛ
            fᵗ.ω -= isdepot(nᵗ) ? r.ω : nᵗ.ω
            isdepot(nᵗ) ? r.ω = ω : nᵗ.ω = ω
            fᵗ.ω += isdepot(nᵗ) ? r.ω : nᵗ.ω
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
        end
        # update costs
        s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
        s.πᵒ += (isdepot(nᵗ) ? r.ω : nᵗ.ω) * fᵗ.πᵒ
        s.πᵖ += 0.
    else
        f = d.F[v.jᵛ]
        # update costs
        s.πᶠ -= isopt(f) ? f.πᶠ : 0.
        s.πᵒ -= r.ω * f.πᵒ
        s.πᵖ -= 0.
        # update route characteristics
        r.tˢ  = d.tˢ
        r.tᵉ  = r.tˢ
        f.ω  -= r.ω
        r.ω   = 0.
        f.ω  += r.ω
        # update costs
        s.πᶠ += isopt(f) ? f.πᶠ : 0.
        s.πᵒ += r.ω * f.πᵒ
        s.πᵖ += 0.
    end
    f     = d.F[v.jᵛ]
    θˡ    = s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (θˡ > r.θ) * (θˡ - r.θ)
    return s
end