# Config.jl: SimulatorConfig holds run settings (read count, noise, output paths, etc.).

"""
    SimulatorConfig(; n_reads, min_unique_vdj=0, seed=0, output_path="calls.tsv", subject_id="sim_donor", noise=NoiseConfig())

Run settings for the simulator. Purpose of main fields:
- `n_reads`: how many reads (rows) to generate.
- `min_unique_vdj`: if > 0, keep sampling until this many unique (V,D,J) combinations (e.g. for RAbHIT).
- `seed`: RNG seed for reproducibility.
- `output_path`, `subject_id`: where to write and how to name the donor in output.
- `noise`: call-corruption config (allele swap, D dropout); use `NoiseConfigNone` for no corruption.
"""
struct SimulatorConfig
    n_reads::Int
    min_unique_vdj::Int
    seed::Int
    output_path::String
    subject_id::String
    noise::NoiseConfig
end

function Base.show(io::IO, ::MIME"text/plain", cfg::SimulatorConfig)
    print(io, "SimulatorConfig(")
    print(io, "n_reads=", cfg.n_reads, ", min_unique_vdj=", cfg.min_unique_vdj,
          ", seed=", cfg.seed, ", output_path=\"", cfg.output_path, "\"",
          ", subject_id=\"", cfg.subject_id, "\", noise=")
    show(io, MIME("text/plain"), cfg.noise)
    print(io, ")")
end

function SimulatorConfig(;
    n_reads::Int,
    min_unique_vdj::Int = 0,
    seed::Integer = 0,
    output_path::String = "calls.tsv",
    subject_id::String = "sim_donor",
    noise::NoiseConfig = NoiseConfig(),
)
    n_reads > 0 || throw(ArgumentError("n_reads must be positive"))
    min_unique_vdj >= 0 || throw(ArgumentError("min_unique_vdj must be non-negative"))
    seed_int = Int(seed & 0xffffffff)
    SimulatorConfig(n_reads, min_unique_vdj, seed_int, output_path, subject_id, noise)
end

const RABHIT_MIN_UNIQUE_VDJ = 2000
