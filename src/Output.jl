# Output.jl: Build MIAIRR-style call table and genotype/truth-phase tables; write TSVs.

locus_prefix(::Type{V}) = "IGHV"
locus_prefix(::Type{D}) = "IGHD"
locus_prefix(::Type{J}) = "IGHJ"

"""Unique allele names present on either chromosome (for genotype table)."""
alleles_list(g::GeneEntry) = unique(filter(!isempty, [phase_allele_string(g.chr1), phase_allele_string(g.chr2)]))

"""Genotype table: one row per (donor, locus, gene, allele). Columns: KI_case (donor id), locus, gene, sequence_id (allele name)."""
function genotype_table(gt::DonorGenotype)
    rows = Tuple{String, String, String, String}[]
    for L in (V, D, J)
        prefix = locus_prefix(L)
        for g in genes(gt, L)
            for a in alleles_list(g)
                push!(rows, (gt.donor_id, prefix, g.gene, a))
            end
        end
    end
    DataFrame(
        KI_case = first.(rows),
        locus = getindex.(rows, 2),
        gene = getindex.(rows, 3),
        sequence_id = last.(rows),
    )
end

"""Allele on chromosome for phase table (empty string if deleted)."""
chr_allele_for_phase(g::GeneEntry, chr::Int) = phase_allele_string(chr == 1 ? g.chr1 : g.chr2)

function truth_phase_table(gt::DonorGenotype)
    rows = NamedTuple{(:gene, :locus, :zygosity, :chr_a, :chr_b), Tuple{String, String, String, String, String}}[]
    for L in (V, D, J)
        prefix = locus_prefix(L)
        for g in genes(gt, L)
            z = zygosity(g)
            z == Homozygous && continue
            push!(rows, (
                gene = g.gene,
                locus = prefix,
                zygosity = zygosity_short(z),
                chr_a = chr_allele_for_phase(g, 1),
                chr_b = chr_allele_for_phase(g, 2),
            ))
        end
    end
    isempty(rows) ? DataFrame() : DataFrame(rows)
end

"""
    build_calls_df(subject_id, v_calls, d_calls, j_calls, sequence_ids, d_start, d_end, duplicate_counts; truth columns...)
"""
function build_calls_df(
    subject_id::String,
    v_calls::Vector{String},
    d_calls::Vector{String},
    j_calls::Vector{String},
    sequence_ids::Vector{String},
    d_germline_starts::Vector{Int},
    d_germline_ends::Vector{Int},
    duplicate_counts::Vector{Int};
    true_chr::Vector{Int} = Int[],
    true_v::Vector{String} = String[],
    true_d::Vector{String} = String[],
    true_j::Vector{String} = String[],
    noise_v::Vector{String} = String[],
    noise_d::Vector{String} = String[],
    noise_j::Vector{String} = String[],
)
    n = length(v_calls)
    DataFrame(
        subject = fill(subject_id, n),
        sequence_id = sequence_ids,
        v_call = v_calls,
        d_call = d_calls,
        j_call = j_calls,
        d_germline_start = d_germline_starts,
        d_germline_end = d_germline_ends,
        duplicate_count = duplicate_counts,
        locus = fill("IGH", n),
        productive = fill(true, n),
        true_chr = isempty(true_chr) ? fill(0, n) : true_chr,
        true_v = isempty(true_v) ? fill("", n) : true_v,
        true_d = isempty(true_d) ? fill("", n) : true_d,
        true_j = isempty(true_j) ? fill("", n) : true_j,
        noise_v = isempty(noise_v) ? fill("none", n) : noise_v,
        noise_d = isempty(noise_d) ? fill("none", n) : noise_d,
        noise_j = isempty(noise_j) ? fill("none", n) : noise_j,
    )
end

write_calls(df::DataFrame, path::String) = CSV.write(path, df; delim='\t')
write_genotype(df::DataFrame, path::String) = CSV.write(path, df; delim='\t')
write_truth_phase(df::DataFrame, path::String) = CSV.write(path, df; delim='\t')

function write_simulation_output(
    calls_df::DataFrame,
    genotype_df::DataFrame,
    truth_phase_df::DataFrame;
    calls_path::Union{String, Nothing} = "calls.tsv",
    genotype_path::Union{String, Nothing} = nothing,
    truth_phase_path::Union{String, Nothing} = nothing,
)
    calls_path !== nothing && write_calls(calls_df, calls_path)
    genotype_path !== nothing && write_genotype(genotype_df, genotype_path)
    truth_phase_path !== nothing && write_truth_phase(truth_phase_df, truth_phase_path)
end
