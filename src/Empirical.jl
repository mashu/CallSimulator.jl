# Empirical.jl: Parameters estimated or fitted from real data (or defaults).
# Used to build genotype and usage so simulations follow real data distributions.

"""
    EmpiricalParams

Parameters derived from real data (or defaults): gene pools, usage skew
(lognormal sigma, allele imbalance), and default anchor J fraction.
Can be extended with kernel densities or other fitted distributions.
"""
struct EmpiricalParams
    gene_pool_v::Vector{String}
    gene_pool_d::Vector{String}
    gene_pool_j::Vector{String}
    lognormal_sigma::Float64
    allele_imbalance::Float64
    anchor_j_fraction_default::Float64
end

function EmpiricalParams(;
    gene_pool_v::Vector{String} = default_gene_pool_v(),
    gene_pool_d::Vector{String} = default_gene_pool_d(),
    gene_pool_j::Vector{String} = default_gene_pool_j(),
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_default::Real = 0.2)
    lognormal_sigma >= 0.0 || throw(ArgumentError("lognormal_sigma must be non-negative"))
    allele_imbalance >= 1.0 || throw(ArgumentError("allele_imbalance must be >= 1"))
    (0.0 < anchor_j_fraction_default < 1.0) ||
        throw(ArgumentError("anchor_j_fraction_default must be in (0, 1)"))
    EmpiricalParams(
        gene_pool_v,
        gene_pool_d,
        gene_pool_j,
        Float64(lognormal_sigma),
        Float64(allele_imbalance),
        Float64(anchor_j_fraction_default),
    )
end

"""Default V gene names (no allele suffix). Pool size matches real-data range (V ~38–49)."""
function default_gene_pool_v()
    [
        "IGHV1-2", "IGHV1-3", "IGHV1-8", "IGHV2-5", "IGHV2-26", "IGHV2-70",
        "IGHV3-7", "IGHV3-9", "IGHV3-11", "IGHV3-13", "IGHV3-15", "IGHV3-20",
        "IGHV3-21", "IGHV3-23", "IGHV3-30", "IGHV3-33", "IGHV3-43", "IGHV3-48",
        "IGHV3-49", "IGHV3-53", "IGHV3-64", "IGHV3-66", "IGHV3-69", "IGHV3-72",
        "IGHV3-73", "IGHV3-74", "IGHV4-4", "IGHV4-28", "IGHV4-30-2", "IGHV4-31",
        "IGHV4-34", "IGHV4-39", "IGHV4-55", "IGHV4-59", "IGHV4-61", "IGHV5-51",
        "IGHV6-1", "IGHV7-4-1",
    ]
end

"""Default D gene names."""
function default_gene_pool_d()
    [
        "IGHD1-1", "IGHD1-7", "IGHD1-14", "IGHD1-20", "IGHD1-26", "IGHD2-2",
        "IGHD2-8", "IGHD2-15", "IGHD2-21", "IGHD3-3", "IGHD3-9", "IGHD3-10",
        "IGHD3-16", "IGHD3-22", "IGHD4-11", "IGHD4-17", "IGHD4-23", "IGHD5-5",
        "IGHD5-12", "IGHD5-18", "IGHD6-6", "IGHD6-13", "IGHD6-19", "IGHD7-27",
    ]
end

"""Default J gene names."""
function default_gene_pool_j()
    ["IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6"]
end

"""
    ki_donor_preset() -> EmpiricalParams

Default preset with estimated values from KI donors: gene pools (V/D/J), usage skew
(lognormal_sigma, allele_imbalance), and default anchor J fraction. Use as the default
for simulating data comparable to real donor cohorts.
"""
function ki_donor_preset(;
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_default::Real = 0.2,
)
    EmpiricalParams(
        gene_pool_v = default_gene_pool_v(),
        gene_pool_d = default_gene_pool_d(),
        gene_pool_j = default_gene_pool_j(),
        lognormal_sigma = lognormal_sigma,
        allele_imbalance = allele_imbalance,
        anchor_j_fraction_default = anchor_j_fraction_default,
    )
end
