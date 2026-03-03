# Expression.jl: Per-gene, per-chromosome weights; callable ExpressionProfile samples by locus/chr.

using Random
abstract type ExpressionMethod end

struct LogNormalExpr <: ExpressionMethod
    σ::Float64
end
(m::LogNormalExpr)(rng::Random.AbstractRNG) = exp(-m.σ^2 / 2 + m.σ * randn(rng))

struct ZipfExpr <: ExpressionMethod
    α::Float64
end
(m::ZipfExpr)(rng::Random.AbstractRNG) = 1.0 / (rand(rng) * 10 + 0.5)^m.α

struct DirichletExpr <: ExpressionMethod
    α::Float64
end
(m::DirichletExpr)(rng::Random.AbstractRNG) = (α = max(m.α, 0.01); randexp(rng) * rand(rng)^(1.0 / α - 1.0))

struct UniformExpr <: ExpressionMethod end
(::UniformExpr)(rng::Random.AbstractRNG) = 1.0

struct LocusWeights
    chr1::Vector{Tuple{Int, Float64}}
    chr2::Vector{Tuple{Int, Float64}}
end

"""
    ExpressionProfile

Holds (gene_index, weight) per locus per chromosome. Gene index is within that locus.
Access via `weights(ep, L, chr)`; sample via `(ep)(rng, gt, L, chr)`.
"""
struct ExpressionProfile
    v::LocusWeights
    d::LocusWeights
    j::LocusWeights
end

weights(ep::ExpressionProfile, ::Type{V}, chr::Int) = chr == 1 ? ep.v.chr1 : ep.v.chr2
weights(ep::ExpressionProfile, ::Type{D}, chr::Int) = chr == 1 ? ep.d.chr1 : ep.d.chr2
weights(ep::ExpressionProfile, ::Type{J}, chr::Int) = chr == 1 ? ep.j.chr1 : ep.j.chr2

"""Weighted sample from (index, weight) pairs; returns index."""
function weighted_sample_index(rng::Random.AbstractRNG, pairs::Vector{Tuple{Int, Float64}})
    isempty(pairs) && error("Cannot sample from empty weights")
    total = sum(last, pairs)
    total > 0.0 || return first(pairs)[1]
    u = rand(rng) * total
    cum = 0.0
    for (idx, w) in pairs
        cum += w
        cum >= u && return idx
    end
    last(pairs)[1]
end

"""Functors: (ep::ExpressionProfile)(rng, gt, L, chr) -> gene index within genes(gt, L)."""
(ep::ExpressionProfile)(rng::Random.AbstractRNG, gt::DonorGenotype, ::Type{L}, chr::Int) where L<:Locus =
    weighted_sample_index(rng, weights(ep, L, chr))

"""Log-uniform draw in [lo, hi] for allele imbalance."""
draw_imbalance(rng::Random.AbstractRNG, lo::Float64, hi::Float64) =
    exp(log(lo) + rand(rng) * (log(hi) - log(lo)))

function build_expression(
    gt::DonorGenotype;
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
    method::ExpressionMethod = LogNormalExpr(1.0),
    allele_imbalance_range::Tuple{Float64, Float64} = (0.3, 3.0),
    cis_sigma::Float64 = 0.35,
    anchor_j_fraction::Union{Nothing, Float64} = nothing,
)
    wv1, wv2 = Tuple{Int, Float64}[], Tuple{Int, Float64}[]
    wd1, wd2 = Tuple{Int, Float64}[], Tuple{Int, Float64}[]
    wj1, wj2 = Tuple{Int, Float64}[], Tuple{Int, Float64}[]

    for (gs, out1, out2) in (
            (gt.genes_v, wv1, wv2),
            (gt.genes_d, wd1, wd2),
            (gt.genes_j, wj1, wj2),
        )
        for (i, g) in enumerate(gs)
            base = method(rng)
            cis1 = exp(cis_sigma * randn(rng))
            cis2 = exp(cis_sigma * randn(rng))
            w1 = w2 = 0.0
            if is_available(g, 1)
                w1 = base * cis1
            end
            if is_available(g, 2)
                w2 = base * cis2
            end
            if is_available(g, 1) && is_available(g, 2) && zygosity(g) == Heterozygous
                r = draw_imbalance(rng, allele_imbalance_range[1], allele_imbalance_range[2])
                w1 *= r
                w2 /= r
            end
            w1 > 0 && push!(out1, (i, w1))
            w2 > 0 && push!(out2, (i, w2))
        end
    end

    if anchor_j_fraction !== nothing
        (0.0 < anchor_j_fraction < 1.0) || throw(ArgumentError("anchor_j_fraction must be in (0,1)"))
        rescale_anchor_j!(wj1, wj2, gt, anchor_j_fraction)
    end

    ExpressionProfile(
        LocusWeights(wv1, wv2),
        LocusWeights(wd1, wd2),
        LocusWeights(wj1, wj2),
    )
end

function rescale_anchor_j!(
    wj1::Vector{Tuple{Int, Float64}},
    wj2::Vector{Tuple{Int, Float64}},
    gt::DonorGenotype,
    target_frac::Float64,
)
    het_j_idx = Set{Int}()
    for (i, g) in enumerate(gt.genes_j)
        zygosity(g) == Heterozygous && push!(het_j_idx, i)
    end
    isempty(het_j_idx) && return
    for wj in (wj1, wj2)
        total = sum(last, wj)
        anchor_sum = sum(w for (i, w) in wj if i in het_j_idx)
        other = total - anchor_sum
        (anchor_sum > 0 && other > 0) || continue
        scale = (target_frac * other / (1 - target_frac)) / anchor_sum
        for k in eachindex(wj)
            if wj[k][1] in het_j_idx
                wj[k] = (wj[k][1], wj[k][2] * scale)
            end
        end
    end
    nothing
end
