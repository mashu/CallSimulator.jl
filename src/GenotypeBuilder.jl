# GenotypeBuilder.jl: Build DonorGenotype from pools and zygosity specs.

using Random

"""Indices for a random subset of size n from 1:pool_len (or all if n >= pool_len)."""
random_subset_indices(rng::Random.AbstractRNG, pool_len::Int, n::Int) =
    n >= pool_len ? collect(1:pool_len) : Random.shuffle(rng, collect(1:pool_len))[1:n]

struct ZygositySpec
    n_homozygous::Int
    n_heterozygous::Int
    n_hemizygous::Int
end

ZygositySpec(; hom::Int=0, het::Int=0, hemi::Int=0) = ZygositySpec(hom, het, hemi)
total_genes(z::ZygositySpec) = z.n_homozygous + z.n_heterozygous + z.n_hemizygous

locus_output_vector(ov::Vector{GeneEntry{V}}, od::Vector{GeneEntry{D}}, oj::Vector{GeneEntry{J}}, ::Type{V}) = ov
locus_output_vector(ov::Vector{GeneEntry{V}}, od::Vector{GeneEntry{D}}, oj::Vector{GeneEntry{J}}, ::Type{D}) = od
locus_output_vector(ov::Vector{GeneEntry{V}}, od::Vector{GeneEntry{D}}, oj::Vector{GeneEntry{J}}, ::Type{J}) = oj

function build_locus_genes!(
    out_v::Vector{GeneEntry{V}},
    out_d::Vector{GeneEntry{D}},
    out_j::Vector{GeneEntry{J}},
    ::Type{L},
    spec::ZygositySpec,
    pool::Vector{String},
    rng::Random.AbstractRNG,
    allele_suffix_hom::String = "*01",
    allele_suffix_het::Tuple{String, String} = ("*01", "*02"),
) where L<:Locus
    n_needed = total_genes(spec)
    length(pool) >= n_needed ||
        throw(ArgumentError("Pool has $(length(pool)) genes but need $n_needed"))
    selected = random_subset_indices(rng, length(pool), n_needed)
    idx = 1
    out = locus_output_vector(out_v, out_d, out_j, L)

    for _ in 1:spec.n_homozygous
        g = pool[selected[idx]]
        push!(out, homozygous_gene(g, L, g * allele_suffix_hom))
        idx += 1
    end
    for _ in 1:spec.n_heterozygous
        g = pool[selected[idx]]
        push!(out, heterozygous_gene(g, L, g * allele_suffix_het[1], g * allele_suffix_het[2]))
        idx += 1
    end
    for _ in 1:spec.n_hemizygous
        g = pool[selected[idx]]
        push!(out, hemizygous_gene(g, L, g * allele_suffix_hom; on_chr=rand(rng, 1:2)))
        idx += 1
    end
    nothing
end

"""
    build_genotype(donor_id, v_spec, d_spec, j_spec; pool_v, pool_d, pool_j, rng, allele_suffix_hom, allele_suffix_het) -> DonorGenotype

Build a diploid genotype. By default alleles use *01 (homozygous/hemizygous) and *01/*02 (heterozygous).
Use `allele_suffix_het = ("*03", "*04")` etc. to get other allele names.
"""
function build_genotype(
    donor_id::String,
    v_spec::ZygositySpec,
    d_spec::ZygositySpec,
    j_spec::ZygositySpec;
    pool_v::Vector{String} = default_gene_pool(V),
    pool_d::Vector{String} = default_gene_pool(D),
    pool_j::Vector{String} = default_gene_pool(J),
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
    allele_suffix_hom::String = "*01",
    allele_suffix_het::Tuple{String, String} = ("*01", "*02"),
)
    gv = GeneEntry{V}[]
    gd = GeneEntry{D}[]
    gj = GeneEntry{J}[]
    build_locus_genes!(gv, gd, gj, V, v_spec, pool_v, rng, allele_suffix_hom, allele_suffix_het)
    build_locus_genes!(gv, gd, gj, D, d_spec, pool_d, rng, allele_suffix_hom, allele_suffix_het)
    build_locus_genes!(gv, gd, gj, J, j_spec, pool_j, rng, allele_suffix_hom, allele_suffix_het)
    DonorGenotype(donor_id, gv, gd, gj)
end

"""
    build_donor_genotype(donor_id; pool_v, pool_d, pool_j, n_genes_v, n_genes_d, n_genes_j, n_v_het_range, n_j_het_range, p_hemi_*, rng) -> DonorGenotype

Build a donor genotype with typical gene counts and zygosity.

- **Gene counts:** Human IGH donors typically have ~38–46 V, ~20–27 D, and 6 J genes.
  By default the donor gets all genes in the provided pools. Pass `n_genes_v`, `n_genes_d`,
  `n_genes_j` to use a random subset of that size (e.g. `n_genes_v = 42`).
- **Zygosity:** Heterozygous counts are drawn from `n_v_het_range` and `n_j_het_range`.
  Default `n_j_het_range = (1, 1)` gives exactly one heterozygous J (typical); use `(0, 1)` for
  sometimes 0 sometimes 1, or `(1, 2)` for 1 or 2. At least one heterozygous J is required
  when the range minimum is ≥ 1. Hemizygosity is applied with probabilities `p_hemi_v`, `p_hemi_d`, `p_hemi_j`.
- **Co-occurrence:** Each simulated read carries exactly one V, one D, and one J; they are
  always sampled from the *same* chromosome. So allele co-occurrence is by phase (no
  cross-chromosome V–D–J pairs).
"""
function build_donor_genotype(
    donor_id::String;
    pool_v::Vector{String} = default_gene_pool(V),
    pool_d::Vector{String} = default_gene_pool(D),
    pool_j::Vector{String} = default_gene_pool(J),
    n_genes_v::Union{Int, Nothing} = nothing,
    n_genes_d::Union{Int, Nothing} = nothing,
    n_genes_j::Union{Int, Nothing} = nothing,
    n_v_het_range::Tuple{Int, Int} = (5, 15),
    n_j_het_range::Tuple{Int, Int} = (1, 1),  # typically 1 het J; use (0,1) or (1,2) for other distributions
    p_hemi_v::Real = 0.03,
    p_hemi_d::Real = 0.05,
    p_hemi_j::Real = 0.001,
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
)
    n_v_pool, n_d_pool, n_j_pool = length(pool_v), length(pool_d), length(pool_j)
    n_v = n_genes_v === nothing ? n_v_pool : min(n_genes_v, n_v_pool)
    n_d = n_genes_d === nothing ? n_d_pool : min(n_genes_d, n_d_pool)
    n_j = n_genes_j === nothing ? n_j_pool : min(n_genes_j, n_j_pool)
    n_v >= 1 || throw(ArgumentError("n_genes_v / pool_v must give at least 1 V gene"))
    n_d >= 0 || throw(ArgumentError("n_genes_d / pool_d must give at least 0 D genes"))
    n_j >= 1 || throw(ArgumentError("n_genes_j / pool_j must give at least 1 J gene"))

    use_pool_v = pool_v[random_subset_indices(rng, n_v_pool, n_v)]
    use_pool_d = n_d == 0 ? String[] : pool_d[random_subset_indices(rng, n_d_pool, n_d)]
    use_pool_j = pool_j[random_subset_indices(rng, n_j_pool, n_j)]
    n_v, n_d, n_j = length(use_pool_v), length(use_pool_d), length(use_pool_j)

    n_v_het = rand(rng, first(n_v_het_range):last(n_v_het_range))
    n_j_het = rand(rng, first(n_j_het_range):last(n_j_het_range))
    n_j_het >= 1 || throw(ArgumentError("At least one heterozygous J required for phasing"))
    n_v_hom = n_v - n_v_het
    n_d_het = n_d == 0 ? 0 : rand(rng, 0:min(8, n_d ÷ 2))
    n_d_hom = n_d - n_d_het
    n_j_hom = n_j - n_j_het
    n_hemi_v = clamp(round(Int, p_hemi_v * n_v), 0, n_v - 1)
    n_hemi_d = n_d == 0 ? 0 : clamp(round(Int, p_hemi_d * n_d), 0, n_d - 1)
    n_hemi_j = clamp(round(Int, p_hemi_j * n_j), 0, n_j - 1)
    n_v_het = min(n_v_het, n_v - n_hemi_v)
    n_v_hom = n_v - n_v_het - n_hemi_v
    n_d_het = min(n_d_het, n_d - n_hemi_d)
    n_d_hom = n_d - n_d_het - n_hemi_d
    n_j_het = min(n_j_het, n_j - n_hemi_j)
    n_j_hom = n_j - n_j_het - n_hemi_j

    build_genotype(
        donor_id,
        ZygositySpec(hom=n_v_hom, het=n_v_het, hemi=n_hemi_v),
        ZygositySpec(hom=n_d_hom, het=n_d_het, hemi=n_hemi_d),
        ZygositySpec(hom=n_j_hom, het=n_j_het, hemi=n_hemi_j);
        pool_v = use_pool_v,
        pool_d = use_pool_d,
        pool_j = use_pool_j,
        rng = rng,
    )
end
