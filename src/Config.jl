# Config.jl: SimulatorConfig holds run settings (read count, noise, output paths, etc.).

"""
    SimulatorConfig(; n_reads, min_unique_vdj=0, seed=0, output_path="calls.tsv", subject_id="sim_donor", noise=NoiseConfig(), anchor_j_fraction=nothing)

Run settings for the simulator. Purpose of main fields:
- `n_reads`: how many reads (rows) to generate.
- `min_unique_vdj`: if > 0, keep sampling until this many unique (V,D,J) combinations (e.g. for RAbHIT).
- `seed`: RNG seed for reproducibility.
- `output_path`, `subject_id`: where to write and how to name the donor in output.
- `noise`: call-corruption config (allele swap, D dropout); use `NoiseConfigNone` for no corruption.
- `anchor_j_fraction`: optional fraction of reads using an anchor J gene (phasing studies).
"""
struct SimulatorConfig
    n_reads::Int
    min_unique_vdj::Int
    seed::UInt32
    output_path::String
    subject_id::String
    noise::NoiseConfig
    anchor_j_fraction::Union{Nothing, Float64}
end

function SimulatorConfig(;
    n_reads::Int,
    min_unique_vdj::Int = 0,
    seed::Integer = 0,
    output_path::String = "calls.tsv",
    subject_id::String = "sim_donor",
    noise::NoiseConfig = NoiseConfig(),
    anchor_j_fraction::Union{Nothing, Real} = nothing,
)
    n_reads > 0 || throw(ArgumentError("n_reads must be positive"))
    min_unique_vdj >= 0 || throw(ArgumentError("min_unique_vdj must be non-negative"))
    anc = anchor_j_fraction === nothing ? nothing : Float64(anchor_j_fraction)
    (anc === nothing || (0.0 < anc < 1.0)) || throw(ArgumentError("anchor_j_fraction must be in (0, 1)"))
    SimulatorConfig(
        n_reads,
        min_unique_vdj,
        UInt32(seed & 0xffffffff),
        output_path,
        subject_id,
        noise,
        anc,
    )
end

const RABHIT_MIN_UNIQUE_VDJ = 2000
