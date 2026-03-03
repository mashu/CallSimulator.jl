# Noise model

The simulator can **corrupt** the called allele after sampling. Configurable via `NoiseConfig(...)`. Use the predefined `NoiseConfigNone` when you want no corruption (calls equal ground truth).

## Configurable rates

| Parameter | Meaning | Typical |
|-----------|---------|---------|
| `p_allele_v`, `p_allele_d`, `p_allele_j` | Allele swap (miscall): report the other allele of the same gene | e.g. 0.04, 0.08, 0.02 |
| `p_gene_v`, `p_gene_d`, `p_gene_j` | Gene confusion: replace with allele from a similar gene (confusion pair) on same chromosome | e.g. 0.02, 0.05, 0 |
| `p_d_dropout` | D dropout: report empty D call | e.g. 0.15 |
| `p_novel` | Novel allele: report e.g. *99 | e.g. 0.005 |

All are keyword arguments to `NoiseConfig(...)`. Optional: `v_confusion`, `d_confusion`, `j_confusion` (pairs of gene names); defaults come from `confusion_pairs(L)`.

## No corruption

```julia
sim = Simulator(gt, config, expression, NoiseModel(NoiseConfigNone), rng)
```

Simulated calls will match the sampled (true) alleles; only ground-truth columns and phasing still differ from “perfect” if you use multiple chromosomes.

## Order of application

For each call: D dropout → gene confusion → allele swap → novel allele. The first that triggers wins. Ground truth stores the *true* (pre-noise) allele and `noise_v` / `noise_d` / `noise_j` record the applied noise type.
