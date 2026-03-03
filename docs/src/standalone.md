# Standalone usage (CLI)

Run the simulator from the command line without writing Julia code.

## Invocation

From the project root:

```bash
julia -e 'using CallSimulator; run_cli()' -- [options]
```

Or with the project activated:

```bash
julia --project=. -e 'using CallSimulator; run_cli()' -- [options]
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `--n-reads` | Number of reads (calls) to generate | 10000 |
| `--min-unique-vdj` | Minimum unique (V,D,J) combinations (0 = disable; RAbHIT recommends 2000) | 0 |
| `--seed` | RNG seed | 42 |
| `--output`, `-o` | Output TSV path for MIAIRR-style calls | calls.tsv |
| `--genotype` | Optional path to save ground-truth genotype table (TSV) | (none) |
| `--truth-phase` | Optional path to save ground-truth phased genotype, heterozygous/hemizygous genes (TSV) | (none) |
| `--subject` | Subject/donor identifier | sim_donor |
| `--p-v`, `--p-d`, `--p-j` | Allele swap (miscall) probability for V, D, J | 0.05, 0.12, 0.03 |
| `--p-d-dropout` | D-call dropout probability (empty D) | 0.2 |
| `--p-gene-confusion-v`, `--p-gene-confusion-d`, `--p-gene-confusion-j` | Gene confusion probability | 0.02, 0.01, 0.0 |
| `--p-novel-allele` | Novel allele call probability | 0.002 |
| `--mean-duplicate-count` | Mean duplicate_count per read | 1.3 |
| `--cis-sigma` | Cis-effect sigma for expression | 0.35 |
| `--anchor-j-fraction` | Fraction of reads carrying anchor J (empty = use empirical default) | (default) |

The CLI builds a donor with `build_donor_genotype` and `ki_donor_preset()` by default.

## Examples

Generate 5000 calls:

```bash
julia -e 'using CallSimulator; run_cli()' -- --n-reads 5000 -o my_calls.tsv
```

Save ground-truth genotype and truth-phase table:

```bash
julia -e 'using CallSimulator; run_cli()' -- --n-reads 10000 -o calls.tsv --genotype genotype.tsv --truth-phase truth_phase.tsv
```

RAbHIT-style (at least 2000 unique VDJ):

```bash
julia -e 'using CallSimulator; run_cli()' -- --n-reads 15000 --min-unique-vdj 2000 -o calls.tsv --genotype genotype.tsv --truth-phase phase.tsv
```
