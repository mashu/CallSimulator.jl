# CallSimulator.jl: Simulate V/D/J call tables (MIAIRR-style) with ground-truth phased genotype.

module CallSimulator

using CSV
using DataFrames
using Random
using ArgParse
using Requires

include("Config.jl")
include("Empirical.jl")
include("Genotype.jl")
include("Usage.jl")
include("MissCall.jl")
include("Output.jl")
include("Simulator.jl")
include("CLI.jl")

function __init__()
    @require KernelDensity = "5ab0869b-81e5-5117-8b39-4d8f8c2e4c2a" include("EmpiricalKDE.jl")
end

export SimulatorConfig,
    MissCallConfig,
    EmpiricalParams,
    ki_donor_preset,
    Gene,
    Genotype,
    is_heterozygous,
    allele_on_chr,
    locus_of,
    locus_alleles,
    WeightedUsage,
    BernoulliMissCall,
    miss_prob,
    Simulator,
    simulate,
    build_genotype,
    build_genotype_realistic,
    uniform_usage,
    skewed_usage,
    skewed_usage_kde,
    build_calls_df,
    genotype_table,
    truth_phase_table,
    write_calls,
    write_genotype,
    write_truth_phase,
    write_simulation_output,
    config_from_args,
    run_cli,
    RABHIT_MIN_UNIQUE_VDJ

end
