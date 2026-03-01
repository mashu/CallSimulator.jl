# MissCall.jl: Miss-call model (allele swap with probability). Functor interface.

"""
    BernoulliMissCall

Miss-call (allele swap) probabilities for V, D, and J loci. Optional per-gene override
for V only (see MissCallConfig). Callable:
(m::BernoulliMissCall)(rng, true_allele, other_allele, locus, gene_name) -> called_allele
"""
struct BernoulliMissCall
    p_v::Float64
    p_d::Float64
    p_j::Float64
    per_gene_v::Vector{Pair{String, Float64}}
end

function BernoulliMissCall(cfg::MissCallConfig)
    BernoulliMissCall(cfg.p_v, cfg.p_d, cfg.p_j, cfg.per_gene_v)
end

"""Return probability for locus (and optional gene for V)."""
function miss_prob(m::BernoulliMissCall, locus::Symbol, gene_name::String)
    if locus === :V
        for p in m.per_gene_v
            p.first == gene_name && return p.second
        end
        return m.p_v
    elseif locus === :D
        return m.p_d
    else
        return m.p_j
    end
end

"""Functors: (m::BernoulliMissCall)(rng, true_allele, other_allele, locus, gene_name) -> called allele."""
function (m::BernoulliMissCall)(
    rng::Random.AbstractRNG,
    true_allele::String,
    other_allele::String,
    locus::Symbol,
    gene_name::String,
)
    p = miss_prob(m, locus, gene_name)
    p <= 0.0 && return true_allele
    rand(rng) < p ? other_allele : true_allele
end
