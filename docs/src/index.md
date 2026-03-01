# CallSimulator

Simulate V/D/J gene call tables (MIAIRR-style) with **ground-truth phased genotype**. Produces both simulated calls and the reference truth; simple and well-structured.

## Design

- **Config** (`SimulatorConfig`, `MissCallConfig`): user settings (read count, seed, output paths, miss-call rates, optional anchor J fraction).
- **Empirical** (`EmpiricalParams`, `ki_donor_preset`): default preset from KI donor–estimated values (gene pools, usage skew). Optional KDE-based usage when KernelDensity is loaded.
- **Genotype** (`Gene`, `Genotype`): ground-truth phased genotype; `build_genotype` / `build_genotype_realistic` build it from gene pools.
- **Usage** (`WeightedUsage`, `skewed_usage`, `skewed_usage_kde`): allele sampling weights; optional KDE-estimated distribution for usage.
- **MissCall** (`BernoulliMissCall`): per-locus (and optional per-gene V) swap probability.
- **Simulator** (`Simulator`, `simulate`): runs the loop; returns MIAIRR calls plus genotype and truth-phase tables.
- **Output**: `build_calls_df`, `genotype_table`, `truth_phase_table`, `write_calls`, `write_genotype`, `write_truth_phase`, `write_simulation_output` for MIAIRR and ground-truth TSVs.

Output is **MIAIRR-style call table** plus **reference ground-truth phased genotype**; you can save both via the API or CLI.

## Quick start (package API)

```julia
using CallSimulator

# Default: KI donor preset
config = SimulatorConfig(n_reads = 10_000, min_unique_vdj = 2000, seed = 42, output_path = "calls.tsv")
estimated = ki_donor_preset()
sim = Simulator(config, estimated)
calls_df, genotype_df, truth_phase_df = simulate(sim)

# Save MIAIRR calls and ground truth
write_simulation_output(calls_df, genotype_df, truth_phase_df;
    calls_path = "calls.tsv",
    genotype_path = "genotype.tsv",
    truth_phase_path = "truth_phase.tsv",
)
```

See [Standalone usage](standalone.md) for the CLI, [Miss-call model](misscall.md) for swap probabilities (V, D, J and per-gene V override), and [API reference](api.md) for the full list of functions and types.
