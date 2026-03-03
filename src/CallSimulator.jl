# CallSimulator.jl: Simulate V/D/J call tables (MIAIRR-style) with phased genotype.
# Idiomatic Julia: Locus types (V, D, J), AlleleState hierarchy, parametric GeneEntry{L},
# locus-indexed DonorGenotype, callable ExpressionProfile and NoiseModel.
#
# API layering:
# - SimulatorConfig = run settings only (n_reads, seed, output_path, subject_id, noise).
# - GenePools + EmpiricalParams = donor and expression recipe (gene lists, p_hemi_*, lognormal_sigma, allele_imbalance, anchor_j_fraction_range).
# - Simulator(config, preset) or Simulator(config, GenePools(), params) builds genotype and expression from that recipe, then holds config + built genotype + built expression + rng. Noise is read from config when sampling.

module CallSimulator

using CSV
using DataFrames
using Random

include("Types.jl")
include("GenePools.jl")
include("GenotypeBuilder.jl")
include("Expression.jl")
include("Noise.jl")
include("Config.jl")
include("Empirical.jl")
include("Output.jl")
include("Simulator.jl")
include("Validation.jl")

export
    # Locus & types
    Locus, V, D, J,
    AlleleState, Present, Deleted,
    GeneEntry, DonorGenotype,
    Zygosity, Homozygous, Heterozygous, Hemizygous,
    zygosity, zygosity_short, zygosity_counts, is_available, is_present, allele_on, allele_display_string, phase_allele_string, genes, has_locus,
    homozygous_gene, heterozygous_gene, hemizygous_gene,
    # Gene pools & genotype build
    default_gene_pool, GenePools,
    ZygositySpec, build_donor_genotype,
    # Expression
    ExpressionMethod, LogNormalExpr, UniformExpr,
    ExpressionProfile, build_expression, weights,
    # Noise
    NoiseConfig, NoiseModel, NoiseConfigNone,
    NoiseType, NoNoise, AlleleSwap, D_Dropout,
    noise_type_string,
    # Config & run
    SimulatorConfig,
    EmpiricalParams, ki_donor_preset,
    genotype_table, truth_phase_table, build_calls_df,
    write_calls, write_genotype, write_truth_phase, write_simulation_output,
    Simulator, simulate,
    check_airr_columns, summarize_simulation, compare_with_real, gini
end
