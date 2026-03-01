# Config.jl: User-facing configuration for the call simulator.
# All fields are concrete types; no Any, no try-catch.

"""
    MissCallConfig

Per-locus miss-call probabilities for V, D, and J. With probability p the true allele
is replaced by the other allele (heterozygous) or left unchanged (homozygous).
Optional per-gene overrides for the V locus only (gene name => probability); D and J
use a single probability per locus.
"""
struct MissCallConfig
    p_v::Float64
    p_d::Float64
    p_j::Float64
    per_gene_v::Vector{Pair{String, Float64}}
end

function MissCallConfig(p_v::Real, p_d::Real, p_j::Real;
                        per_gene_v::Vector{Pair{String, Float64}} = Pair{String, Float64}[])
    (0.0 <= p_v <= 1.0 && 0.0 <= p_d <= 1.0 && 0.0 <= p_j <= 1.0) ||
        throw(ArgumentError("Miss-call rates must be in [0, 1]"))
    MissCallConfig(Float64(p_v), Float64(p_d), Float64(p_j), per_gene_v)
end

"""
    SimulatorConfig

User-controlled settings: read count, RAbHIT constraint, RNG seed, output path,
subject id, miss-call config, and optional anchor J fraction override.
"""
struct SimulatorConfig
    n_reads::Int
    min_unique_vdj::Int
    seed::UInt32
    output_path::String
    subject_id::String
    miss_call::MissCallConfig
    anchor_j_fraction::Union{Nothing, Float64}
end

function SimulatorConfig(;
    n_reads::Int,
    min_unique_vdj::Int = 0,
    seed::Integer = 0,
    output_path::String = "calls.tsv",
    subject_id::String = "sim_donor",
    miss_call::MissCallConfig = MissCallConfig(0.05, 0.12, 0.03),
    anchor_j_fraction::Union{Nothing, Real} = nothing)
    n_reads > 0 || throw(ArgumentError("n_reads must be positive"))
    min_unique_vdj >= 0 || throw(ArgumentError("min_unique_vdj must be non-negative"))
    anc = anchor_j_fraction === nothing ? nothing : Float64(anchor_j_fraction)
    if anc !== nothing
        (anc > 0.0 && anc < 1.0) || throw(ArgumentError("anchor_j_fraction must be in (0, 1)"))
    end
    SimulatorConfig(
        n_reads,
        min_unique_vdj,
        UInt32(seed & 0xffffffff),
        output_path,
        subject_id,
        miss_call,
        anc,
    )
end
