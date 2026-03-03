# CallSimulator.jl

[![CI](https://github.com/mashu/CallSimulator.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/mashu/CallSimulator.jl/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/mashu/CallSimulator.jl/graph/badge.svg)](https://codecov.io/gh/mashu/CallSimulator.jl/branch/main)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/CallSimulator.jl/stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Simulate MIAIRR-style V/D/J call tables with ground-truth phased genotype for testing phasing and assignment. Each read has one V, one D, one J from the same chromosome; the simulator records true alleles and optional noise (allele swap, D dropout). Genotypes support **homozygous**, **heterozygous**, and **hemizygous** genes; expression can skew toward one allele per heterozygous gene (**allele imbalance**).

## Example

```julia
using CallSimulator

config = SimulatorConfig(n_reads = 10_000, output_path = "calls.tsv")
sim = Simulator(config, ki_donor_preset())
calls_df, genotype_df, truth_phase_df = simulate(sim)
write_simulation_output(calls_df, genotype_df, truth_phase_df;
    calls_path = "calls.tsv", genotype_path = "genotype.tsv", truth_phase_path = "truth_phase.tsv")
```

Gene pools (V/D/J lists) are separate from tunable parameters: `Simulator(config, GenePools(), ki_donor_preset(allele_imbalance=2.0))` uses default genes and only overrides params.
