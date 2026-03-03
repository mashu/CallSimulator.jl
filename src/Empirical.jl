# Empirical.jl: Parameters estimated or fitted from real data (or defaults).
# Used to build genotype and usage so simulations follow real data distributions.

"""
    EmpiricalParams

Parameters derived from real data (or defaults): gene pools, usage skew
(lognormal sigma, allele imbalance), default anchor J fraction, and
per-locus hemizygosity probabilities (V/D/J).
Can be extended with kernel densities or other fitted distributions.
"""
struct EmpiricalParams
    gene_pool_v::Vector{String}
    gene_pool_d::Vector{String}
    gene_pool_j::Vector{String}
    lognormal_sigma::Float64
    allele_imbalance::Float64
    anchor_j_fraction_default::Float64
    p_hemi_v::Float64
    p_hemi_d::Float64
    p_hemi_j::Float64
end

function EmpiricalParams(;
    gene_pool_v::Vector{String} = default_gene_pool(V),
    gene_pool_d::Vector{String} = default_gene_pool(D),
    gene_pool_j::Vector{String} = default_gene_pool(J),
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_default::Real = 0.2,
    p_hemi_v::Real = 0.03,
    p_hemi_d::Real = 0.05,
    p_hemi_j::Real = 0.001,
)
    lognormal_sigma >= 0.0 || throw(ArgumentError("lognormal_sigma must be non-negative"))
    allele_imbalance >= 1.0 || throw(ArgumentError("allele_imbalance must be >= 1"))
    (0.0 < anchor_j_fraction_default < 1.0) ||
        throw(ArgumentError("anchor_j_fraction_default must be in (0, 1)"))
    (0.0 <= p_hemi_v <= 1.0) || throw(ArgumentError("p_hemi_v must be in [0, 1]"))
    (0.0 <= p_hemi_d <= 1.0) || throw(ArgumentError("p_hemi_d must be in [0, 1]"))
    (0.0 <= p_hemi_j <= 1.0) || throw(ArgumentError("p_hemi_j must be in [0, 1]"))
    EmpiricalParams(
        gene_pool_v,
        gene_pool_d,
        gene_pool_j,
        Float64(lognormal_sigma),
        Float64(allele_imbalance),
        Float64(anchor_j_fraction_default),
        Float64(p_hemi_v),
        Float64(p_hemi_d),
        Float64(p_hemi_j),
    )
end

"""
    ki_donor_preset() -> EmpiricalParams

Default preset with estimated values from KI donors: gene pools (V/D/J), usage skew
(lognormal_sigma, allele_imbalance), default anchor J fraction, and baseline
hemizygosity rates (p_hemi_v/d/j). Use as the default
for simulating data comparable to real donor cohorts.
"""
function ki_donor_preset(;
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_default::Real = 0.2,
    p_hemi_v::Real = 0.03,
    p_hemi_d::Real = 0.05,
    p_hemi_j::Real = 0.001,
)
    EmpiricalParams(
        gene_pool_v = default_gene_pool_v(),
        gene_pool_d = default_gene_pool_d(),
        gene_pool_j = default_gene_pool_j(),
        lognormal_sigma = lognormal_sigma,
        allele_imbalance = allele_imbalance,
        anchor_j_fraction_default = anchor_j_fraction_default,
        p_hemi_v = p_hemi_v,
        p_hemi_d = p_hemi_d,
        p_hemi_j = p_hemi_j,
    )
end
