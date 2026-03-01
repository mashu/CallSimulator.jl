# Genotype.jl: Ground-truth genotype (genes with phased alleles per chromosome).

"""
    Gene

One gene in the genotype. `allele_a` on chromosome A, `allele_b` on chromosome B.
Homozygous when allele_a == allele_b.
"""
struct Gene
    gene::String
    locus::Symbol  # :V, :D, or :J
    allele_a::String
    allele_b::String
end

Gene(gene::String, locus::Symbol, allele::String) = Gene(gene, locus, allele, allele)

"""True if the gene has two distinct alleles."""
is_heterozygous(g::Gene) = g.allele_a != g.allele_b

"""
    Genotype

Ground-truth genotype for one donor: list of genes with known phase.
"""
struct Genotype
    donor::String
    genes::Vector{Gene}
end

Base.length(gt::Genotype) = length(gt.genes)

"""Allele on chromosome `chr` (1 = A, 2 = B) for gene at index `i`."""
function allele_on_chr(gt::Genotype, i::Int, chr::Int)
    g = gt.genes[i]
    chr == 1 ? g.allele_a : g.allele_b
end

"""Locus of gene at index `i`."""
locus_of(gt::Genotype, i::Int) = gt.genes[i].locus

"""Build locus index: (v_list, d_list, j_list) of (gene_index, allele_a, allele_b)."""
function locus_alleles(gt::Genotype)
    v_list = Tuple{Int, String, String}[]
    d_list = Tuple{Int, String, String}[]
    j_list = Tuple{Int, String, String}[]
    for (i, g) in enumerate(gt.genes)
        t = (i, g.allele_a, g.allele_b)
        if g.locus === :V
            push!(v_list, t)
        elseif g.locus === :D
            push!(d_list, t)
        else
            push!(j_list, t)
        end
    end
    (v_list, d_list, j_list)
end

"""Build genotype from counts and gene pools (realistic donor-like)."""
function build_genotype(
    donor::String,
    n_v_het::Int, n_v_hom::Int,
    n_d_het::Int, n_d_hom::Int,
    n_j_het::Int, n_j_hom::Int;
    gene_pool_v::Vector{String},
    gene_pool_d::Vector{String},
    gene_pool_j::Vector{String},
    allele_suffix::String = "*01",
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
)
    need_v = n_v_het + n_v_hom
    need_d = n_d_het + n_d_hom
    need_j = n_j_het + n_j_hom
    (length(gene_pool_v) >= need_v && length(gene_pool_d) >= need_d && length(gene_pool_j) >= need_j) ||
        throw(ArgumentError("Gene pools too small for requested counts"))

    perm_v = Random.shuffle(rng, collect(eachindex(gene_pool_v)))[1:need_v]
    perm_d = Random.shuffle(rng, collect(eachindex(gene_pool_d)))[1:need_d]
    perm_j = Random.shuffle(rng, collect(eachindex(gene_pool_j)))[1:need_j]

    genes = Gene[]
    for i in 1:n_v_het
        g = gene_pool_v[perm_v[i]]
        push!(genes, Gene(g, :V, g * "*01", g * "*02"))
    end
    for i in (n_v_het + 1):need_v
        g = gene_pool_v[perm_v[i]]
        push!(genes, Gene(g, :V, g * allele_suffix))
    end
    for i in 1:n_d_het
        g = gene_pool_d[perm_d[i]]
        push!(genes, Gene(g, :D, g * "*01", g * "*02"))
    end
    for i in (n_d_het + 1):need_d
        g = gene_pool_d[perm_d[i]]
        push!(genes, Gene(g, :D, g * allele_suffix))
    end
    for i in 1:n_j_het
        g = gene_pool_j[perm_j[i]]
        push!(genes, Gene(g, :J, g * "*01", g * "*02"))
    end
    for i in (n_j_het + 1):need_j
        g = gene_pool_j[perm_j[i]]
        push!(genes, Gene(g, :J, g * allele_suffix))
    end
    Genotype(donor, genes)
end

"""RAbHIT-friendly genotype: at least one heterozygous J. Gene pools must be provided (e.g. from EmpiricalParams)."""
function build_genotype_realistic(
    donor::String,
    gene_pool_v::Vector{String},
    gene_pool_d::Vector{String},
    gene_pool_j::Vector{String};
    n_v_het_range::Tuple{Int, Int} = (5, 15),
    n_j_het_range::Tuple{Int, Int} = (1, 1),
    rng::Random.AbstractRNG = Random.GLOBAL_RNG,
)
    n_v = length(gene_pool_v)
    n_d = length(gene_pool_d)
    n_j = length(gene_pool_j)
    rv = n_v_het_range[1]:n_v_het_range[2]
    rj = n_j_het_range[1]:n_j_het_range[2]
    n_v_het = min(last(rv), max(first(rv), rand(rng, rv)))
    n_j_het = min(last(rj), max(first(rj), rand(rng, rj)))
    n_j_het >= 1 || throw(ArgumentError("RAbHIT requires at least one heterozygous J"))
    n_v_hom = n_v - n_v_het
    n_d_het = rand(rng, 0:min(8, n_d ÷ 2))
    n_d_hom = n_d - n_d_het
    n_j_hom = n_j - n_j_het
    build_genotype(
        donor, n_v_het, n_v_hom, n_d_het, n_d_hom, n_j_het, n_j_hom;
        gene_pool_v, gene_pool_d, gene_pool_j, rng,
    )
end
