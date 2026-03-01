# Miss-call model

The simulator can **swap** the called allele (V, D, or J) with the other allele on the same chromosome with a given probability. This happens at **all three loci**: V, D, and J. We are not only swapping V.

## Per-locus probabilities

- **`p_v`**: probability that the true V allele is replaced by the other allele (when heterozygous).
- **`p_d`**: same for D.
- **`p_j`**: same for J.

So each of V, D, and J has its own miss-call (swap) rate. Defaults in the CLI are e.g. 0.05 (V), 0.12 (D), 0.03 (J).

## Per-gene V override (optional)

For the **V locus only**, you can optionally set a **per-gene** swap probability. That is: a different probability for specific V genes (e.g. `IGHV1-2` ⇒ 0.08) instead of using the single `p_v` for all V genes. This is useful when some V genes are known to be mis-called more often in practice.

- **D and J** do not have per-gene overrides; they use a single probability per locus (`p_d`, `p_j`).
- **Per-gene V** is configured via `MissCallConfig(..., per_gene_v = ["IGHV1-2" => 0.08, ...])`; it is not exposed in the CLI by default.

Summary: miss-call applies to **V, D, and J**; only V has an optional per-gene override on top of the locus-wide `p_v`.
