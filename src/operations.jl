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
    c.jⁿ  = c.jⁿ
    c.iᵛ  = r.iᵛ
    c.iᵈ  = r.iᵈ
    c.r   = r
    isdepot(nᵗ) ? r.iˢ = c.iⁿ : nᵗ.iʰ = c.iⁿ
    isdepot(nʰ) ? r.iᵉ = c.iⁿ : nʰ.iᵗ = c.iⁿ
    c.iʰ  = nʰ.iⁿ
    c.iᵗ  = nᵗ.iⁿ
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ)
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    r.x   = (r.n * r.x + c.x)/(r.n + 1)
    r.y   = (r.n * r.y + c.y)/(r.n + 1)
    r.n  += 1
    r.q  += isdelivery(c) ? c.qᶜ : 0.
    r.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    # update associated vehicle
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    v.n  += 1
    v.q  += isdelivery(c) ? c.qᶜ : 0.
    v.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.q * d.πᵒ
    d.n  += 1
    d.q  += isdelivery(c) ? c.qᶜ : 0.
    d.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.q * d.πᵒ
    # update en-route parameters
    s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)
    s.πᵖ -= (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)
    s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = r.tˢ
        n  = 1
        q  = 0.
        l  = s.A[(d.iⁿ,cˢ.iⁿ)].l
        c  = cˢ
        while true
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ -= (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ -= (c.l > v.lᵛ) * (c.l - v.lᵛ)
            c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
            c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.n   = n
            c.q   = q
            c.l   = l
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ += (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ += (c.l > v.lᵛ) * (c.l - v.lᵛ)
            if isequal(c, cᵉ) break end
            tᵈ = c.tᵈ
            n += 1
            q += -c.qᶜ
            l += s.A[(c.iⁿ,c.iʰ)].l
            c  = s.C[c.iʰ]
        end
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
    end
    (v.tˢ, v.tᵉ) = (r.tˢ, r.tᵉ)
    s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)
    s.πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)
    s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
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
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
    cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
    # update associated customer nodes
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ)
    s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
    s.πᵖ -= (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
    s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
    s.πᵖ -= (c.l > v.lᵛ) * (c.l - v.lᵛ)
    c.jⁿ  = c.jⁿ
    c.iᵛ  = 0
    c.iᵈ  = 0
    c.r   = NullRoute
    isdepot(nᵗ) ? r.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
    isdepot(nʰ) ? r.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
    c.iʰ  = 0
    c.iᵗ  = 0
    c.tᵃ  = isdelivery(c) ? c.tˡ : c.tᵉ
    c.tᵈ  = c.tᵃ + c.τᶜ
    c.n   = 0
    c.q   = 0.
    c.l   = 0.
    qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ)
    s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
    s.πᵖ += (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
    s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
    s.πᵖ += (c.l > v.lᵛ) * (c.l - v.lᵛ)
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    r.x   = isone(r.n) ? 0. : (r.n * r.x - c.x)/(r.n - 1)
    r.y   = isone(r.n) ? 0. : (r.n * r.y - c.y)/(r.n - 1)
    r.n  -= 1
    r.q  -= isdelivery(c) ? c.qᶜ : 0.
    r.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    # update associated vehicle
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    v.n  -= 1
    v.q  -= isdelivery(c) ? c.qᶜ : 0.
    v.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.q * d.πᵒ
    d.n  -= 1
    d.q  -= isdelivery(c) ? c.qᶜ : 0.
    d.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.q * d.πᵒ
    # update en-route parameters
    s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)
    s.πᵖ -= (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)
    s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = d.tˢ
        n  = 1
        q  = 0.
        l  = s.A[(d.iⁿ,cˢ.iⁿ)].l
        c  = cˢ
        while true
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ -= (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ -= (c.l > v.lᵛ) * (c.l - v.lᵛ)
            c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
            c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.n   = n
            c.q   = q
            c.l   = l
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ += (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ += (c.l > v.lᵛ) * (c.l - v.lᵛ)
            if isequal(c, cᵉ) break end
            tᵈ = c.tᵈ
            n += 1
            q += -c.qᶜ
            l += s.A[(c.iⁿ,c.iʰ)].l
            c  = s.C[c.iʰ]
        end
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
    end
    (v.tˢ, v.tᵉ) = (r.tˢ, r.tᵉ)
    s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)
    s.πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)
    s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
    return s
end



"""
    movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)

Returns solution `s` after moving vehicle `v` from fleet of `d₁` into fleet 
of depot node `d₂`.
"""
function movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)
    deleteat!(d₁.V, findfirst(isequal(v), d₁.V))
    push!(d₂.V, v)
    v.iᵈ = d₂.iⁿ
    r = v.r
    if isopt(r)
        nᵗ  = s.C[r.iᵉ]
        nʰ  = s.C[r.iˢ]
        aᵗ₁ = s.A[(nᵗ.iⁿ, d₁.iⁿ)]
        aᵗ₂ = s.A[(nᵗ.iⁿ, d₂.iⁿ)]
        aʰ₁ = s.A[(d₁.iⁿ, nʰ.iⁿ)]
        aʰ₂ = s.A[(d₂.iⁿ, nʰ.iⁿ)]
        # update depot d₁
        s.πᶠ -= isopt(d₁) ? d₁.πᶠ : 0.
        s.πᵒ -= d₁.q * d₁.πᵒ
        d₁.n -= r.n
        d₁.q -= r.q
        d₁.l -= r.l
        s.πᶠ += isopt(d₁) ? d₁.πᶠ : 0.
        s.πᵒ += d₁.q * d₁.πᵒ
        # update associated customer nodes
        nᵗ.iʰ = d₂.iⁿ
        nʰ.iᵗ = d₂.iⁿ
        # update associated route
        s.πᶠ -= 0.
        s.πᵒ -= r.l * v.πᵈ
        r.iᵈ = d₂.iⁿ
        r.iˢ = nʰ.iⁿ
        r.iᵉ = nᵗ.iⁿ
        r.l -= aᵗ₁.l + aʰ₁.l
        r.l += aᵗ₂.l + aʰ₂.l
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = d₂.tˢ
        n  = 1
        q  = 0.
        l  = s.A[(d₂.iⁿ,cˢ.iⁿ)].l
        c  = cˢ
        while true
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ -= (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ -= (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ -= (c.l > v.lᵛ) * (c.l - v.lᵛ)
            c.iᵈ  = d₂.iⁿ
            c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
            c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.n   = n
            c.q   = q
            c.l   = l
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ)
            s.πᵖ += (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ)
            s.πᵖ += (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ)
            s.πᵖ += (c.l > v.lᵛ) * (c.l - v.lᵛ)
            if isequal(c, cᵉ) break end
            tᵈ = c.tᵈ
            n += 1
            q += -c.qᶜ
            l += s.A[(c.iⁿ,c.iʰ)].l
            c  = s.C[c.iʰ]
        end
        r.tˢ  = d₂.tˢ
        r.tᵉ  = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d₂.iⁿ)].l/v.sᵛ
        s.πᶠ += 0.
        s.πᵒ += r.l * v.πᵈ
        # update associated vehicle
        s.πᶠ -= isopt(v) ? v.πᶠ : 0.
        s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ -= (d₁.tˢ > v.tˢ) * (d₁.tˢ - v.tˢ)
        s.πᵖ -= (v.tᵉ > d₁.tᵉ) * (v.tᵉ - d₁.tᵉ)
        s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
        v.tˢ  = r.tˢ
        v.tᵉ  = r.tᵉ
        v.l  -= aᵗ₁.l + aʰ₁.l
        v.l  += aᵗ₂.l + aʰ₂.l
        s.πᶠ += isopt(v) ? v.πᶠ : 0.
        s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ += (d₂.tˢ > v.tˢ) * (d₂.tˢ - v.tˢ)
        s.πᵖ += (v.tᵉ > d₂.tᵉ) * (v.tᵉ - d₂.tᵉ)
        s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) * ((v.tᵉ - v.tˢ) - v.τʷ)
        # update depot d₂
        s.πᶠ -= isopt(d₂) ? d₂.πᶠ : 0.
        s.πᵒ -= d₂.q * d₂.πᵒ
        d₂.n += r.n
        d₂.q += r.q
        d₂.l += r.l
        s.πᶠ += isopt(d₂) ? d₂.πᶠ : 0.
        s.πᵒ += d₂.q * d₂.πᵒ
    else
        r.iᵈ = d₂.iⁿ
        r.iˢ = d₂.iⁿ
        r.iᵉ = d₂.iⁿ
        r.tˢ = d₂.tˢ
        r.tᵉ = d₂.tˢ
        v.tˢ = d₂.tˢ
        v.tᵉ = d₂.tˢ
    end
    return s
end