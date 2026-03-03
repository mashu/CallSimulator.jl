# CallSimulator

Simulate V/D/J gene call tables (MIAIRR-style) with **ground-truth phased genotype** for testing phasing and assignment methods. Each read carries one V, one D, and one J call; the simulator records the true chromosome of origin and applied noise.

## Design

**In one sentence:** Build a donor genotype and expression from gene pools, add a noise model, then run the simulator to get a MIAIRR call table plus ground-truth genotype and phase.

**Pipeline (top to bottom):**

```
   Gene pools (V,D,J)  +  Config (n_reads, seed, noise, ...)
              |
              v
   DonorGenotype   (build_donor_genotype)
              |
              v
   ExpressionProfile  (build_expression)
              |
              v
   Simulator(genotype, config, expression, rng)  -->  simulate()
              |
              v
   calls_df    genotype_df    truth_phase_df
```

**What each piece does:**

1. **Genotype** — `DonorGenotype` = which V/D/J genes the donor has and zygosity: homozygous, heterozygous, or hemizygous (gene on one chromosome only, useful as phasing anchor). Built from pools with `build_donor_genotype`.
2. **Expression** — Weights per gene per chromosome; used to sample which V, which D, which J on each read. For **heterozygous** genes only, **allele imbalance** is applied: a ratio *r* (chr1 vs chr2 weight) is sampled once per gene from a log-uniform distribution on `allele_imbalance_range` (e.g. (0.3, 3.0)). Homozygous and hemizygous genes have no imbalance. Built with `build_expression(gt; ...)`.
3. **Noise** — Corrupts calls: allele swap, D dropout. Set in `SimulatorConfig` via `noise=NoiseConfig()` (default) or `noise=NoiseConfigNone`.
4. **Simulator** — Holds genotype, config (which includes noise), expression, rng; for each read picks one chromosome, samples V,D,J from it, applies noise from config, returns the three tables.

**Shortcut:** `Simulator(config, ki_donor_preset())` builds genotype and expression from preset defaults; you still choose noise in `config`.

## Quick start

```julia
using CallSimulator
using Random

# Build a donor (typical gene counts)
gt = build_donor_genotype("donor1"; n_genes_v = 42, n_genes_d = 23, n_genes_j = 6, rng = MersenneTwister(42))
estimated = ki_donor_preset()
cfg = SimulatorConfig(n_reads = 500, output_path = "calls.tsv")
expression = build_expression(gt; rng = MersenneTwister(cfg.seed + 2),
    method = LogNormalExpr(estimated.lognormal_sigma),
    allele_imbalance_range = (1.0 / estimated.allele_imbalance, estimated.allele_imbalance))
sim = Simulator(gt, cfg, expression, MersenneTwister(cfg.seed))
calls_df, genotype_df, truth_phase_df = simulate(sim)
```

Or use the preset to build genotype and expression in one go:

```julia
config = SimulatorConfig(n_reads = 10_000, output_path = "calls.tsv")
sim = Simulator(config, ki_donor_preset())
calls_df, genotype_df, truth_phase_df = simulate(sim)
```

See [Noise model](noise.md) for configuring call corruption and [API reference](api.md) for the full list.
