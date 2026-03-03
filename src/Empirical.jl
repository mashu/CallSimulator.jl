# Empirical.jl: Parameters estimated or fitted from real data (or defaults).
# Used with GenePools to build genotype and expression; gene pools are separate so you can tune params without re-specifying genes.

"""
    EmpiricalParams

Tunable parameters (no gene pools): usage skew (lognormal sigma, allele imbalance),
anchor J fraction range (min, max) for sampling, and per-locus hemizygosity probabilities (V/D/J).
Gene pools are passed separately via `GenePools` (e.g. `Simulator(config, GenePools(), ki_donor_preset())`).
"""
struct EmpiricalParams
    lognormal_sigma::Float64
    allele_imbalance::Float64
    anchor_j_fraction_range::Tuple{Float64, Float64}
    p_hemi_v::Float64
    p_hemi_d::Float64
    p_hemi_j::Float64
end

function Base.show(io::IO, ::MIME"text/plain", ep::EmpiricalParams)
    print(io, "EmpiricalParams(")
    print(io, "lognormal_sigma=", ep.lognormal_sigma, ", allele_imbalance=", ep.allele_imbalance,
          ", anchor_j_fraction_range=", ep.anchor_j_fraction_range, ", ")
    print(io, "p_hemi_v=", ep.p_hemi_v, ", p_hemi_d=", ep.p_hemi_d, ", p_hemi_j=", ep.p_hemi_j, ")")
end

function EmpiricalParams(;
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_range::Tuple{Real, Real} = (0.1, 0.3),
    p_hemi_v::Real = 0.03,
    p_hemi_d::Real = 0.05,
    p_hemi_j::Real = 0.001,
)
    lognormal_sigma >= 0.0 || throw(ArgumentError("lognormal_sigma must be non-negative"))
    allele_imbalance >= 1.0 || throw(ArgumentError("allele_imbalance must be >= 1"))
    ar = (Float64(anchor_j_fraction_range[1]), Float64(anchor_j_fraction_range[2]))
    (0.0 < ar[1] < ar[2] < 1.0) ||
        throw(ArgumentError("anchor_j_fraction_range must be (min, max) with 0 < min < max < 1"))
    (0.0 <= p_hemi_v <= 1.0) || throw(ArgumentError("p_hemi_v must be in [0, 1]"))
    (0.0 <= p_hemi_d <= 1.0) || throw(ArgumentError("p_hemi_d must be in [0, 1]"))
    (0.0 <= p_hemi_j <= 1.0) || throw(ArgumentError("p_hemi_j must be in [0, 1]"))
    EmpiricalParams(
        Float64(lognormal_sigma),
        Float64(allele_imbalance),
        ar,
        Float64(p_hemi_v),
        Float64(p_hemi_d),
        Float64(p_hemi_j),
    )
end

"""
    ki_donor_preset() -> EmpiricalParams

Tunable parameters with defaults from KI donors: usage skew (lognormal_sigma, allele_imbalance),
anchor J fraction range, and hemizygosity rates (p_hemi_v/d/j).
Example: `Simulator(config, GenePools(), ki_donor_preset(allele_imbalance=2.0))`.
"""
function ki_donor_preset(;
    lognormal_sigma::Real = 1.0,
    allele_imbalance::Real = 1.0,
    anchor_j_fraction_range::Tuple{Real, Real} = (0.1, 0.3),
    p_hemi_v::Real = 0.03,
    p_hemi_d::Real = 0.05,
    p_hemi_j::Real = 0.001,
)
    EmpiricalParams(
        lognormal_sigma = lognormal_sigma,
        allele_imbalance = allele_imbalance,
        anchor_j_fraction_range = anchor_j_fraction_range,
        p_hemi_v = p_hemi_v,
        p_hemi_d = p_hemi_d,
        p_hemi_j = p_hemi_j,
    )
end
