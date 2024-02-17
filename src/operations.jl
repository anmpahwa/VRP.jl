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
    ωˡ    = (s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ) * v.ωᵛ
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (ωˡ > r.ω) * (ωˡ - r.ω)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        # initiate iterative parameters
        nᵗ = d
        nʰ = c
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        ω  = 1.
        q  = 0.
        c  = cˢ
        while true
            # fetch network features
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            a′ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            ωˡ = (a′.l/v.lᵛ) * v.ωᵛ
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= 0.
            s.πᵒ -= 0.
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
            # update node characteristics
            φ    = refuel(nᵗ, nʰ, r, s)
            c.tᵃ = t + φ * (aᵗ.l/v.sᵛ + (v.ωᵛ + (ωˡ - c.ω)) * fᵗ.τᵛ + aʰ.l/v.sᵛ) + (1 - φ) * (aᵒ.l/v.sᵛ)
            c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.ω  = φ * (v.ωᵛ - (aʰ.l/v.lᵛ) * v.ωᵛ) + (1 - φ) * (ω - (aᵒ.l/v.lᵛ) * v.ωᵛ)
            c.q  = q
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += 0.
            s.πᵒ += 0.
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
            # update iterative parameters
            if isequal(c, cᵉ) break end
            nᵗ = c
            nʰ = c.iʰ ≤ lastindex(s.D) ? s.D[c.iʰ] : s.C[c.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fʰ = nʰ.F[v.jᵛ]
            t  = c.tᵈ
            ω  = c.ω
            q -= c.qᶜ
            c  = s.C[c.iʰ]
        end
        aᵒ   = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ   = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ   = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        a′   = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        ωˡ   = (a′.l/v.lᵛ) * v.ωᵛ
        φ    = refuel(nᵗ, nʰ, r, s)
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + φ * (aᵗ.l/v.sᵛ + (v.ωᵛ + (ωˡ - c.ω)) * fᵗ.τᵛ + aʰ.l/v.sᵛ) + (1 - φ) * (aᵒ.l/v.sᵛ)
        r.ω  = φ * (v.ωᵛ - (aʰ.l/v.lᵛ) * v.ωᵛ) + (1 - φ) * (ω - (aᵒ.l/v.lᵛ) * v.ωᵛ)
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
        r.ω  = v.ωᵛ
    end
    f     = d.F[v.jᵛ]
    ωˡ    = (s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ) * v.ωᵛ
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (ωˡ > r.ω) * (ωˡ - r.ω)
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
    ωˡ = (s.A[(c.iⁿ, f.iⁿ)].l/v.lᵛ) * v.ωᵛ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    # update associated customer nodes
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
    if iscustomer(nᵗ) nᵗ.iʰ = nʰ.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = nᵗ.iⁿ end
    c.iᵗ  = 0
    c.iʰ  = 0
    c.tᵃ  = isdelivery(c) ? c.tˡ : c.tᵉ
    c.tᵈ  = c.tᵃ + c.τᶜ
    c.ω   = v.ωᵛ
    c.q   = 0.
    c.r   = NullRoute
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
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
    f     = d.F[v.jᵛ]
    ωˡ    = (s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ) * v.ωᵛ
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (ωˡ > r.ω) * (ωˡ - r.ω)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        # initiate iterative parameters
        nᵗ = d
        nʰ = c
        fᵗ = nᵗ.F[v.jᵛ]
        fʰ = nʰ.F[v.jᵛ]
        t  = d.tˢ
        ω  = 1.
        q  = 0.
        c  = cˢ
        while true
            # fetch network features
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
            aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
            a′ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            ωˡ = (a′.l/v.lᵛ) * v.ωᵛ
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= 0.
            s.πᵒ -= 0.
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
            # update node characteristics
            φ    = refuel(nᵗ, nʰ, r, s)
            c.tᵃ = t + φ * (aᵗ.l/v.sᵛ + (v.ωᵛ + (ωˡ - c.ω)) * fᵗ.τᵛ + aʰ.l/v.sᵛ) + (1 - φ) * (aᵒ.l/v.sᵛ)
            c.tᵈ = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.ω  = φ * (v.ωᵛ - (aʰ.l/v.lᵛ) * v.ωᵛ) + (1 - φ) * (ω - (aᵒ.l/v.lᵛ) * v.ωᵛ)
            c.q  = q
            # update costs
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += 0.
            s.πᵒ += 0.
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (ωˡ > c.ω) * (ωˡ - c.ω)
            # update iterative parameters
            if isequal(c, cᵉ) break end
            nᵗ = c
            nʰ = c.iʰ ≤ lastindex(s.D) ? s.D[c.iʰ] : s.C[c.iʰ]
            fᵗ = nᵗ.F[v.jᵛ]
            fʰ = nʰ.F[v.jᵛ]
            t  = c.tᵈ
            ω  = c.ω
            q -= c.qᶜ
            c  = s.C[c.iʰ]
        end
        aᵒ   = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
        aᵗ   = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
        aʰ   = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
        a′   = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
        ωˡ   = (a′.l/v.lᵛ) * v.ωᵛ
        φ    = refuel(nᵗ, nʰ, r, s)
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + φ * (aᵗ.l/v.sᵛ + (v.ωᵛ + (ωˡ - c.ω)) * fᵗ.τᵛ + aʰ.l/v.sᵛ) + (1 - φ) * (aᵒ.l/v.sᵛ)
        r.ω  = φ * (v.ωᵛ - (aʰ.l/v.lᵛ) * v.ωᵛ) + (1 - φ) * (ω - (aᵒ.l/v.lᵛ) * v.ωᵛ)
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
        r.ω  = v.ωᵛ
    end
    f     = d.F[v.jᵛ]
    ωˡ    = (s.A[(d.iⁿ, f.iⁿ)].l/v.lᵛ) * v.ωᵛ
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ) + (ωˡ > r.ω) * (ωˡ - r.ω)
    return s
end



"""
    refuel(nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns `true` if pre-emptive re-fueling is needed in route `r` at 
node `nᵗ` to reach node `nʰ`.
"""
function refuel(nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    fᵗ = nᵗ.F[v.jᵛ]
    fʰ = nʰ.F[v.jᵛ]
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, fᵗ.iⁿ)]
    aʰ = s.A[(fᵗ.iⁿ, nʰ.iⁿ)]
    a′ = s.A[(nʰ.iⁿ, fʰ.iⁿ)]
    ωˡ = (a′.l/v.lᵛ) * v.ωᵛ
    if isdepot(nᵗ)
        if ωˡ ≥ r.ω
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= fᵗ.ω * f.πᵒ
            s.πᵖ -= 0.
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
            fᵗ.ω += (v.ωᵛ + (ωˡ - r.w))
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += fᵗ.ω * fᵗ.πᵒ
            s.πᵖ += 0.
            return true
        end
    end
    if iscustomer(nᵗ)
        c = nᵗ
        if ωˡ ≥ c.ω
            s.πᶠ -= isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ -= fᵗ.ω * f.πᵒ
            s.πᵖ -= 0.
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
            fᵗ.ω += (v.ωᵛ + (ωˡ - c.w))
            s.πᶠ += isopt(fᵗ) ? fᵗ.πᶠ : 0.
            s.πᵒ += fᵗ.ω * fᵗ.πᵒ
            s.πᵖ += 0.
            return true
        end
    end
    return false
end