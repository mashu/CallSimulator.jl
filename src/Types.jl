# Types.jl: Locus and AlleleState type hierarchy; parametric GeneEntry; DonorGenotype.
# All dispatch is on types — no symbol checks.

# ------------------------------------------------------------------------------
# Locus: singleton types for V, D, J so we dispatch on locus without symbols.
# ------------------------------------------------------------------------------
abstract type Locus end
struct V <: Locus end
struct D <: Locus end
struct J <: Locus end

# ------------------------------------------------------------------------------
# AlleleState: sum type via abstract type + concrete subtypes (not bool + string).
# ------------------------------------------------------------------------------
abstract type AlleleState end

struct Present <: AlleleState
    name::String
end

struct Deleted <: AlleleState
end

is_present(::Present) = true
is_present(::Deleted) = false

allele_name(p::Present) = p.name
allele_name(::Deleted) = error("Deleted allele has no name")

"""Display string for an allele state (show/debug). Present → name, Deleted → em dash."""
allele_display_string(p::Present) = p.name
allele_display_string(::Deleted) = "—"

"""String for phase tables: Present → name, Deleted → \"\"."""
phase_allele_string(p::Present) = p.name
phase_allele_string(::Deleted) = ""

# ------------------------------------------------------------------------------
# Zygosity: for reporting; could be types, enum is compact for tables.
# ------------------------------------------------------------------------------
@enum Zygosity Homozygous Heterozygous Hemizygous

"""Short label for zygosity (hom / het / hemi)."""
zygosity_short(z::Zygosity) = (z == Homozygous && return "hom"; z == Heterozygous && return "het"; "hemi")

# ------------------------------------------------------------------------------
# GeneEntry{L}: one gene at a locus with allele state per chromosome.
# ------------------------------------------------------------------------------
struct GeneEntry{L<:Locus}
    gene::String
    chr1::AlleleState
    chr2::AlleleState
end

function zygosity(g::GeneEntry)
    p1, p2 = is_present(g.chr1), is_present(g.chr2)
    if p1 && p2
        n1, n2 = allele_name(g.chr1), allele_name(g.chr2)
        return n1 == n2 ? Homozygous : Heterozygous
    end
    (p1 ⊻ p2) && return Hemizygous
    error("Gene $(g.gene) deleted on both chromosomes")
end

is_available(g::GeneEntry, chr::Int) = chr == 1 ? is_present(g.chr1) : is_present(g.chr2)

allele_on(g::GeneEntry, chr::Int) = allele_name(chr == 1 ? g.chr1 : g.chr2)

"""Other allele of the same gene (for allele-swap noise); nothing if no other allele."""
function other_allele(g::GeneEntry, current_allele::String)
    zygosity(g) != Heterozygous && return nothing
    if is_present(g.chr1) && allele_name(g.chr1) == current_allele && is_present(g.chr2)
        return allele_name(g.chr2)
    end
    if is_present(g.chr2) && allele_name(g.chr2) == current_allele && is_present(g.chr1)
        return allele_name(g.chr1)
    end
    nothing
end

# Convenience constructors
homozygous_gene(gene::String, ::Type{L}, allele::String) where L<:Locus =
    GeneEntry{L}(gene, Present(allele), Present(allele))

function heterozygous_gene(gene::String, ::Type{L}, a1::String, a2::String) where L<:Locus
    a1 != a2 || error("Heterozygous gene requires distinct alleles")
    GeneEntry{L}(gene, Present(a1), Present(a2))
end

function hemizygous_gene(gene::String, ::Type{L}, allele::String; on_chr::Int=1) where L<:Locus
    on_chr in (1, 2) || error("on_chr must be 1 or 2")
    on_chr == 1 ? GeneEntry{L}(gene, Present(allele), Deleted()) : GeneEntry{L}(gene, Deleted(), Present(allele))
end

# ------------------------------------------------------------------------------
# DonorGenotype: locus-indexed gene vectors — dispatch on Locus type.
# ------------------------------------------------------------------------------
"""
    DonorGenotype

Diploid genotype for one donor: two chromosomes, with genes stored per locus
(`genes_v`, `genes_d`, `genes_j`). Use `genes(gt, V)` etc. and `zygosity_counts(gt, L)`.
Display in REPL or with `show(io, MIME(\"text/plain\"), gt)` to see genes and zygosities.
"""
struct DonorGenotype
    donor_id::String
    genes_v::Vector{GeneEntry{V}}
    genes_d::Vector{GeneEntry{D}}
    genes_j::Vector{GeneEntry{J}}
end

Base.length(gt::DonorGenotype) = length(gt.genes_v) + length(gt.genes_d) + length(gt.genes_j)

genes(gt::DonorGenotype, ::Type{V}) = gt.genes_v
genes(gt::DonorGenotype, ::Type{D}) = gt.genes_d
genes(gt::DonorGenotype, ::Type{J}) = gt.genes_j

function has_locus(gt::DonorGenotype, ::Type{L}) where L<:Locus
    !isempty(genes(gt, L))
end

"""Count genes by zygosity at locus L. Returns (homozygous=, heterozygous=, hemizygous=)."""
function zygosity_counts(gt::DonorGenotype, ::Type{L}) where L<:Locus
    gs = genes(gt, L)
    (homozygous = count(g -> zygosity(g) == Homozygous, gs),
     heterozygous = count(g -> zygosity(g) == Heterozygous, gs),
     hemizygous = count(g -> zygosity(g) == Hemizygous, gs))
end

# ------------------------------------------------------------------------------
# Display: show genotype and gene entries so users can inspect structure.
# ------------------------------------------------------------------------------
function Base.show(io::IO, ::MIME"text/plain", g::GeneEntry{L}) where L<:Locus
    print(io, "GeneEntry{", nameof(L), "} ", g.gene, " [", zygosity_short(zygosity(g)), "] chr1=", allele_display_string(g.chr1), " chr2=", allele_display_string(g.chr2))
end

function Base.show(io::IO, ::MIME"text/plain", gt::DonorGenotype)
    println(io, "DonorGenotype \"", gt.donor_id, "\" (", length(gt), " genes)")
    for L in (V, D, J)
        gs = genes(gt, L)
        isempty(gs) && continue
        c = zygosity_counts(gt, L)
        println(io, "  ", nameof(L), ": ", length(gs), " genes (hom=", c.homozygous, ", het=", c.heterozygous, ", hemi=", c.hemizygous, ")")
        for g in gs
            println(io, "    ", g.gene, " [", zygosity_short(zygosity(g)), "] ", allele_display_string(g.chr1), " | ", allele_display_string(g.chr2))
        end
    end
end
