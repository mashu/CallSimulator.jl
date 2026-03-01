# Standalone usage (CLI)

Run the simulator from the command line without writing Julia code.

## From project root

```bash
julia --project=. scripts/run_cli.jl [options]
```

Or from another directory, with the package on `LOAD_PATH`:

```bash
julia -e 'using CallSimulator; CallSimulator.run_cli()' -- [options]
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `--n-reads` | Number of reads (calls) to generate | 10000 |
| `--min-unique-vdj` | Minimum unique (V,D,J) combinations (0 = disable; RAbHIT recommends 2000) | 0 |
| `--seed` | RNG seed | 42 |
| `--output`, `-o` | Output TSV path for MIAIRR-style calls | calls.tsv |
| `--genotype` | Optional path to save ground-truth genotype table (TSV) | (none) |
| `--truth-phase` | Optional path to save ground-truth phased genotype, heterozygous genes only (TSV) | (none) |
| `--subject` | Subject/donor identifier | sim_donor |
| `--p-v`, `--p-d`, `--p-j` | Miss-call probability for V, D, J | 0.05, 0.12, 0.03 |
| `--anchor-j-fraction` | Fraction of reads carrying anchor J (empty = use default) | (default) |

The CLI uses the **KI donor preset** by default (estimated values from KI donors).

## Examples

Generate 5000 calls and save only the MIAIRR table:

```bash
julia --project=. scripts/run_cli.jl --n-reads 5000 -o my_calls.tsv
```

Generate calls and also save ground-truth genotype and phased truth:

```bash
julia --project=. scripts/run_cli.jl --n-reads 10000 -o calls.tsv --genotype genotype.tsv --truth-phase truth_phase.tsv
```

RAbHIT-style run (at least 2000 unique VDJ):

```bash
julia --project=. scripts/run_cli.jl --n-reads 15000 --min-unique-vdj 2000 -o calls.tsv --genotype genotype.tsv --truth-phase phase.tsv
```
