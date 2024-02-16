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
    # update associated vehicle
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵒ -= 0.
    s.πᵖ -= 0.
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    s.πᵒ += 0.
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
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = r.tˢ
        θ  = 1 - s.A[(d.iⁿ,cˢ.iⁿ)].l/v.lᵛ
        q  = 0.
        c  = cˢ
        while true
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= 0.
            s.πᵒ -= 0.
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
            c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
            c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.θ   = θ
            c.q   = q
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += 0.
            s.πᵒ += 0.
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
            if isequal(c, cᵉ) break end
            tᵈ = c.tᵈ
            q -= c.qᶜ
            θ -= s.A[(c.iⁿ,c.iʰ)].l/v.lᵛ
            c  = s.C[c.iʰ]
        end
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
    end
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ)
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
    s.πᵖ -= (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
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
    s.πᵖ += (!isequal(cᵖ.r, cᵈ.r) && isclose(cᵖ) && isclose(cᵈ)) * abs(c.qᶜ) + (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
    # update associated vehicle-route
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ - 0.
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
    s.πᶠ -= 0.
    s.πᵒ -= (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ)
    if isopt(r)
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = d.tˢ
        θ  = 1 - s.A[(d.iⁿ,cˢ.iⁿ)].l/v.lᵛ
        q  = 0.
        c  = cˢ
        while true
            cᵖ = isdelivery(c) ? s.C[c.jⁿ] : s.C[c.iⁿ] 
            cᵈ = isdelivery(c) ? s.C[c.iⁿ] : s.C[c.jⁿ]
            qᵒ = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ -= 0.
            s.πᵒ -= 0.
            s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
            c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
            c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
            c.θ   = θ
            c.q   = q
            qᵒ    = isdelivery(c) ? c.q : c.q + abs(c.qᶜ)
            s.πᶠ += 0.
            s.πᵒ += 0.
            s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) + (cᵖ.tᵃ > cᵈ.tᵃ) * (cᵖ.tᵃ - cᵈ.tᵃ) + (qᵒ > v.qᵛ) * (qᵒ - v.qᵛ) + (c.θ < 0.) * abs(c.θ)
            if isequal(c, cᵉ) break end
            tᵈ = c.tᵈ
            q -= c.qᶜ
            θ -= s.A[(c.iⁿ,c.iʰ)].l/v.lᵛ
            c  = s.C[c.iʰ]
        end
        r.tˢ = d.tˢ
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
    else
        r.tˢ = d.tˢ
        r.tᵉ = r.tˢ
    end
    s.πᶠ += 0.
    s.πᵒ += (r.tᵉ - r.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > r.tˢ) * (d.tˢ - r.tˢ) + (r.tᵉ > d.tᵉ) * (r.tᵉ - d.tᵉ) + ((r.tᵉ - r.tˢ) > v.τʷ) * ((r.tᵉ - r.tˢ) - v.τʷ)
    return s
end
#=
tᶠ = 0.
θʳ = 1.

tᶠ    = 0. 
θʳ   -= s.A[(c.iᵗ, c.iⁿ)].l/v.lᵛ
if θʳ ≤ v.θˡ 
    f  = c.Iᶠ[v.jᵛ]
    aᵒ = s.A[(c.iⁿ, c.iʰ)]
    aᵗ = s.A[(c.iⁿ, f.iⁿ)]
    aʰ = s.A[(f.iⁿ, c.iʰ)]
    r.l += aᵗ.l + aʰ.l - aᵒ.l
    v.l += aᵗ.l + aʰ.l - aᵒ.l
    d.l += aᵗ.l + aʰ.l - aᵒ.l
    l   += aᵗ.l
    tᶠ  += s.A[(c.iⁿ, f.iⁿ)].l/v.sᵛ
    θʳ  -= s.A[(c.iⁿ, f.iⁿ)].l/v.lᵛ
    tᶠ  += (v.θᵘ - l/v.lᵛ) * f.τᵛ
    θʳ   = v.θᵘ
    tᶠ  += s.A[(f.iⁿ, c.iʰ)].l/v.sᵛ
    θʳ  -= s.A[(f.iⁿ, c.iʰ)].l/v.lᵛ
    l   += aʰ.l
end
=#