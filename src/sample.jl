StatsBase.@weights OffsetWeights

"""
    OffsetWeights(vs, wsum=sum(vs))

Construct an `OffsetWeights` with weight values from (offset) vector `vs`.
A precomputed sum may be provided as `wsum`.
""" OffsetWeights

"""
    sample([rng::AbstractRNG], a::OffsetVector, [wv::AbstractWeights])

Select a single random element of `a`. Sampling probabilities are proportional to
the weights given in `wv`, if provided.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
StatsBase.sample(rng::AbstractRNG, a::OffsetVector) = a[rand(rng, firstindex(a):lastindex(a))]
StatsBase.sample(a::OffsetVector) = sample(Random.GLOBAL_RNG, a)

"""
    sample([rng::AbstractRNG], wv::OffsetWeights)

Select a single random integer in `eachindex(wv)` with probabilities
proportional to the weights given in `wv`.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function StatsBase.sample(rng::AbstractRNG, wv::OffsetWeights)
    t = rand(rng) * sum(wv)
    i = firstindex(wv)
    n = lastindex(wv)
    cw = wv[i]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end
StatsBase.sample(wv::OffsetWeights) = sample(Random.GLOBAL_RNG, wv)

StatsBase.sample(rng::AbstractRNG, a::OffsetVector, wv::OffsetWeights) = a[sample(rng, wv)]
StatsBase.sample(a::OffsetVector, wv::OffsetWeights) = sample(Random.GLOBAL_RNG, a, wv)