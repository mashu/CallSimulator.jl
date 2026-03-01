# EmpiricalKDE.jl: Optional KDE-based usage weights. Loaded when KernelDensity is in the environment.
# Usage: fit KDE to log(weight) samples from real data, then build WeightedUsage from that distribution.

using KernelDensity
using Random

"""
    skewed_usage_kde(gt::Genotype, log_weight_samples::Vector{Float64}; rng = Random.GLOBAL_RNG, anchor_j_fraction = nothing) -> WeightedUsage

Build allele usage weights by sampling per-gene log-weights from a KDE fitted to `log_weight_samples`
(typically from KI donor or real data). Heterozygous genes get two draws; allele imbalance is applied.
Requires `using KernelDensity` to be loaded (optional dependency).
"""
function skewed_usage_kde(
    gt::Genotype,
    log_weight_samples::Vector{Float64};
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
    anchor_j_fraction::Union{Nothing, Float64} = nothing,
)
    length(log_weight_samples) >= 2 || throw(ArgumentError("Need at least 2 log_weight_samples for KDE"))
    k = kde(log_weight_samples)
    # Sample from KDE: pick grid point with probability ∝ density, add small noise
    n_pts = length(k.x)
    w_dens = k.density .+ 1e-10
    w_dens ./= sum(w_dens)
    bw = k.bw[1]

    function sample_log_weight()
        i = _weighted_sample(rng, Vector(w_dens))
        k.x[i] + bw * randn(rng)
    end

    v_pairs = Pair{String, Float64}[]
    d_pairs = Pair{String, Float64}[]
    j_pairs = Pair{String, Float64}[]
    for g in gt.genes
        base = exp(sample_log_weight())
        if g.allele_a == g.allele_b
            pair = g.allele_a => base
            if g.locus === :V
                push!(v_pairs, pair)
            elseif g.locus === :D
                push!(d_pairs, pair)
            else
                push!(j_pairs, pair)
            end
        else
            r = 1.0 / sqrt(2.0) + rand(rng) * (sqrt(2.0) - 1.0 / sqrt(2.0))
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
