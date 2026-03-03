# Noise.jl: Noise model with dispatch on Locus; callable NoiseModel.

using Random

"""Call-corruption rates (allele swap, D dropout) per locus. Pass to `NoiseModel(cfg)`."""
struct NoiseConfig
    p_allele_v::Float64
    p_allele_d::Float64
    p_allele_j::Float64
    p_d_dropout::Float64
end

p_allele_miscall(cfg::NoiseConfig, ::Type{V}) = cfg.p_allele_v
p_allele_miscall(cfg::NoiseConfig, ::Type{D}) = cfg.p_allele_d
p_allele_miscall(cfg::NoiseConfig, ::Type{J}) = cfg.p_allele_j

"""D dropout: only applies to D locus."""
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{D}) =
    rand(rng) < cfg.p_d_dropout
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{V}) = false
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{J}) = false

function NoiseConfig(;
    p_allele_v::Real = 0.04,
    p_allele_d::Real = 0.08,
    p_allele_j::Real = 0.02,
    p_d_dropout::Real = 0.15,
)
    all(p -> 0.0 <= p <= 1.0, [p_allele_v, p_allele_d, p_allele_j, p_d_dropout]) ||
        throw(ArgumentError("All probabilities must be in [0, 1]"))
    NoiseConfig(
        Float64(p_allele_v), Float64(p_allele_d), Float64(p_allele_j),
        Float64(p_d_dropout),
    )
end

"""No corruption (all rates zero). Use with `NoiseModel(NoiseConfigNone)` when calls should equal ground truth."""
const NoiseConfigNone = NoiseConfig(0.0, 0.0, 0.0, 0.0)

# Three outcomes: no change; wrong allele (same gene); D missing.
@enum NoiseType NoNoise AlleleSwap D_Dropout

"""Callable: (rng, gt, L, chr, gene_entry, true_allele) -> (called_allele, NoiseType)."""
struct NoiseModel
    cfg::NoiseConfig
end

function (noise::NoiseModel)(
    rng::Random.AbstractRNG,
    gt::DonorGenotype,
    ::Type{L},
    chr::Int,
    gene_entry::GeneEntry{L},
    true_allele::String,
) where L<:Locus
    cfg = noise.cfg

    if maybe_dropout(cfg, rng, L)
        return ("", D_Dropout)
    end

    if p_allele_miscall(cfg, L) > 0 && rand(rng) < p_allele_miscall(cfg, L)
        other = other_allele(gene_entry, true_allele)
        other !== nothing && return (other, AlleleSwap)
    end

    (true_allele, NoNoise)
end

const NOISE_TYPE_STR = Dict(NoNoise => "none", AlleleSwap => "allele_swap", D_Dropout => "dropout")
noise_type_string(n::NoiseType) = get(NOISE_TYPE_STR, n, "unknown")
