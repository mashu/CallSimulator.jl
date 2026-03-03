# Noise.jl: Noise model with dispatch on Locus; callable NoiseModel.

using Random

"""Call-corruption rates (allele swap, gene confusion, D dropout, novel) per locus. Pass to `NoiseModel(cfg)`."""
struct NoiseConfig
    p_allele_v::Float64
    p_allele_d::Float64
    p_allele_j::Float64
    p_gene_v::Float64
    p_gene_d::Float64
    p_gene_j::Float64
    p_d_dropout::Float64
    p_novel::Float64
    v_confusion::Vector{Tuple{String, String}}
    d_confusion::Vector{Tuple{String, String}}
    j_confusion::Vector{Tuple{String, String}}
end

p_allele_miscall(cfg::NoiseConfig, ::Type{V}) = cfg.p_allele_v
p_allele_miscall(cfg::NoiseConfig, ::Type{D}) = cfg.p_allele_d
p_allele_miscall(cfg::NoiseConfig, ::Type{J}) = cfg.p_allele_j

p_gene_confusion(cfg::NoiseConfig, ::Type{V}) = cfg.p_gene_v
p_gene_confusion(cfg::NoiseConfig, ::Type{D}) = cfg.p_gene_d
p_gene_confusion(cfg::NoiseConfig, ::Type{J}) = cfg.p_gene_j

confusion_tuples(cfg::NoiseConfig, ::Type{V}) = cfg.v_confusion
confusion_tuples(cfg::NoiseConfig, ::Type{D}) = cfg.d_confusion
confusion_tuples(cfg::NoiseConfig, ::Type{J}) = cfg.j_confusion

"""D dropout: only applies to D locus."""
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{D}) =
    rand(rng) < cfg.p_d_dropout
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{V}) = false
maybe_dropout(cfg::NoiseConfig, rng::Random.AbstractRNG, ::Type{J}) = false

function NoiseConfig(;
    p_allele_v::Real = 0.04,
    p_allele_d::Real = 0.08,
    p_allele_j::Real = 0.02,
    p_gene_v::Real = 0.02,
    p_gene_d::Real = 0.05,
    p_gene_j::Real = 0.0,
    p_d_dropout::Real = 0.15,
    p_novel::Real = 0.005,
    v_confusion::Vector{Tuple{String, String}} = confusion_pairs(V),
    d_confusion::Vector{Tuple{String, String}} = confusion_pairs(D),
    j_confusion::Vector{Tuple{String, String}} = confusion_pairs(J),
)
    all(p -> 0.0 <= p <= 1.0, [p_allele_v, p_allele_d, p_allele_j, p_gene_v, p_gene_d, p_gene_j, p_d_dropout, p_novel]) ||
        throw(ArgumentError("All probabilities must be in [0, 1]"))
    NoiseConfig(
        Float64(p_allele_v), Float64(p_allele_d), Float64(p_allele_j),
        Float64(p_gene_v), Float64(p_gene_d), Float64(p_gene_j),
        Float64(p_d_dropout), Float64(p_novel),
        v_confusion, d_confusion, j_confusion,
    )
end

"""No corruption (all rates zero). Use with `NoiseModel(NoiseConfigNone)` when calls should equal ground truth."""
const NoiseConfigNone = NoiseConfig(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    Tuple{String, String}[], Tuple{String, String}[], Tuple{String, String}[],
)

function confused_allele(
    rng::Random.AbstractRNG,
    cfg::NoiseConfig,
    gt::DonorGenotype,
    ::Type{L},
    chr::Int,
    source_gene::String,
) where L<:Locus
    pairs = confusion_tuples(cfg, L)
    partner = nothing
    for (a, b) in pairs
        a == source_gene && (partner = b; break)
        b == source_gene && (partner = a; break)
    end
    partner === nothing && return nothing
    gs = genes(gt, L)
    for g in gs
        g.gene == partner && is_available(g, chr) && return allele_on(g, chr)
    end
    partner * "*01"
end

# Five distinct outcomes: no change; wrong allele (same gene); wrong gene (confusion pair); D missing; non-reference allele.
@enum NoiseType NoNoise AlleleSwap GeneConfusion D_Dropout NovelAllele

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

    if p_gene_confusion(cfg, L) > 0 && rand(rng) < p_gene_confusion(cfg, L)
        alt = confused_allele(rng, cfg, gt, L, chr, gene_entry.gene)
        alt !== nothing && return (alt, GeneConfusion)
    end

    if p_allele_miscall(cfg, L) > 0 && rand(rng) < p_allele_miscall(cfg, L)
        other = other_allele(gene_entry, true_allele)
        other !== nothing && return (other, AlleleSwap)
    end

    if cfg.p_novel > 0 && rand(rng) < cfg.p_novel
        base = something(findlast('*', true_allele), 0)
        name = base === 0 ? true_allele : true_allele[1:prevind(true_allele, base)]
        return (name * "*99", NovelAllele)
    end

    (true_allele, NoNoise)
end

const NOISE_TYPE_STR = Dict(NoNoise => "none", AlleleSwap => "allele_swap", GeneConfusion => "gene_confusion", D_Dropout => "dropout", NovelAllele => "novel_allele")
noise_type_string(n::NoiseType) = get(NOISE_TYPE_STR, n, "unknown")
