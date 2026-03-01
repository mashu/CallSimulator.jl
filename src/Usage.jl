# Usage.jl: Allele usage model (weights per locus). Type-stable sampling via functor.

"""
    WeightedUsage

Weights per allele for each locus. Sampling is via (u::WeightedUsage)(rng, locus, chr_alleles).
"""
struct WeightedUsage
    v::Vector{Pair{String, Float64}}
    d::Vector{Pair{String, Float64}}
    j::Vector{Pair{String, Float64}}
end

"""Weighted sample: returns index in 1:length(weights). Weights non-negative, sum > 0."""
function _weighted_sample(rng::Random.AbstractRNG, weights::Vector{Float64})
    s = sum(weights)
    s <= 0.0 && return 1
    u = rand(rng) * s
    cum = 0.0
    for i in eachindex(weights)
        cum += weights[i]
        cum >= u && return i
    end
    length(weights)
end

"""Functors: (u::WeightedUsage)(rng, locus, chr_alleles) -> sampled allele string."""
function (u::WeightedUsage)(rng::Random.AbstractRNG, locus::Symbol, chr_alleles::Vector{String})
    pairs = locus === :V ? u.v : locus === :D ? u.d : u.j
    weights = [get_weight(pairs, a) for a in chr_alleles]
    idx = _weighted_sample(rng, weights)
    chr_alleles[idx]
end

get_weight(pairs::Vector{Pair{String, Float64}}, allele::String) =
    let w = 1.0; for p in pairs; p.first == allele && (w = p.second; break); end; w end

"""Uniform weights from genotype (all alleles equal weight)."""
function uniform_usage(gt::Genotype)
    v = String[]; d = String[]; j = String[]
    for g in gt.genes
        als = unique([g.allele_a, g.allele_b])
        for a in als
            if g.locus === :V
                push!(v, a)
            elseif g.locus === :D
                push!(d, a)
            else
                push!(j, a)
            end
        end
    end
    unique!(v); unique!(d); unique!(j)
    WeightedUsage(
        [a => 1.0 for a in v],
        [a => 1.0 for a in d],
        [a => 1.0 for a in j],
    )
end

"""Skewed weights (lognormal per gene, allele imbalance). Optionally cap anchor J fraction."""
function skewed_usage(
    gt::Genotype;
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
    lognormal_sigma::Float64 = 1.0,
    allele_imbalance::Float64 = 1.0,
    anchor_j_fraction::Union{Nothing, Float64} = nothing,
)
    v_pairs = Pair{String, Float64}[]
    d_pairs = Pair{String, Float64}[]
    j_pairs = Pair{String, Float64}[]
    for g in gt.genes
        base = exp(lognormal_sigma * randn(rng) - lognormal_sigma^2 / 2)
        if g.allele_a == g.allele_b
            w = base
            pair = g.allele_a => w
            if g.locus === :V
                push!(v_pairs, pair)
            elseif g.locus === :D
                push!(d_pairs, pair)
            else
                push!(j_pairs, pair)
            end
        else
            r = allele_imbalance <= 1.0 ? 1.0 : 1.0 / sqrt(allele_imbalance) + rand(rng) * (sqrt(allele_imbalance) - 1.0 / sqrt(allele_imbalance))
            w_a = base * r
            w_b = base / r
            if g.locus === :V
                push!(v_pairs, g.allele_a => w_a)
                push!(v_pairs, g.allele_b => w_b)
            elseif g.locus === :D
                push!(d_pairs, g.allele_a => w_a)
                push!(d_pairs, g.allele_b => w_b)
            else
                push!(j_pairs, g.allele_a => w_a)
                push!(j_pairs, g.allele_b => w_b)
            end
        end
    end
    if anchor_j_fraction !== nothing
        (anchor_j_fraction > 0.0 && anchor_j_fraction < 1.0) ||
            throw(ArgumentError("anchor_j_fraction must be in (0, 1)"))
        anchor_alleles = Set{String}()
        for g in gt.genes
            if g.locus === :J && is_heterozygous(g)
                push!(anchor_alleles, g.allele_a)
                push!(anchor_alleles, g.allele_b)
            end
        end
        if !isempty(anchor_alleles)
            total_j = sum(last, j_pairs)
            anchor_sum = sum(w for (a, w) in j_pairs if a in anchor_alleles)
            other_sum = total_j - anchor_sum
            if anchor_sum > 0.0 && other_sum > 0.0
                f = Float64(anchor_j_fraction)
                new_anchor = f * other_sum / (1.0 - f)
                scale = new_anchor / anchor_sum
                j_pairs = [a in anchor_alleles ? (a => w * scale) : (a => w) for (a, w) in j_pairs]
            end
        end
    end
    WeightedUsage(v_pairs, d_pairs, j_pairs)
end

"""Stub when KernelDensity is not loaded; replaced by EmpiricalKDE when KernelDensity is available."""
function skewed_usage_kde(gt::Genotype, log_weight_samples::Vector{Float64}; rng = Random.GLOBAL_RNG, anchor_j_fraction = nothing)
    error("skewed_usage_kde requires KernelDensity. Add KernelDensity to your environment, then: using KernelDensity; using CallSimulator")
end
