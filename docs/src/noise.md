# Noise model

The simulator can **corrupt** the called allele after sampling. Configurable via `NoiseConfig(...)`. Use the predefined `NoiseConfigNone` when you want no corruption (calls equal ground truth).

## Configurable rates

| Parameter | Meaning | Typical |
|-----------|---------|---------|
| `p_allele_v`, `p_allele_d`, `p_allele_j` | Allele swap (miscall): report the other allele of the same gene | e.g. 0.04, 0.08, 0.02 |
| `p_d_dropout` | D dropout: report empty D call | e.g. 0.15 |

All are keyword arguments to `NoiseConfig(...)`.

## No corruption

```julia
sim = Simulator(gt, config, expression, NoiseModel(NoiseConfigNone), rng)
```

Simulated calls will match the sampled (true) alleles; only ground-truth columns and phasing still differ from “perfect” if you use multiple chromosomes.

## Order of application

For each call: D dropout → allele swap. The first that triggers wins. Ground truth stores the *true* (pre-noise) allele and `noise_v` / `noise_d` / `noise_j` record the applied noise type.
