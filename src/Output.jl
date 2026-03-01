# Output.jl: Build call table and auxiliary tables (MIAIRR-compatible columns).

"""
    build_calls_df(subject_id, v_calls, d_calls, j_calls, sequence_ids, d_start, d_end) -> DataFrame

DataFrame with columns required by downstream methods: subject, v_call, d_call, j_call,
sequence_id, d_germline_start, d_germline_end, locus, productive.
"""
function build_calls_df(
    subject_id::String,
    v_calls::Vector{String},
    d_calls::Vector{String},
    j_calls::Vector{String},
    sequence_ids::Vector{String},
    d_germline_starts::Vector{Int},
    d_germline_ends::Vector{Int},
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
        locus = fill("IGH", n),
        productive = fill(true, n),
    )
end

"""Genotype table: KI_case (donor), locus prefix, gene, sequence_id (allele)."""
function genotype_table(gt::Genotype)
    rows = Tuple{String, String, String, String}[]
    for g in gt.genes
        prefix = g.locus === :V ? "IGHV" : g.locus === :D ? "IGHD" : "IGHJ"
        for a in unique([g.allele_a, g.allele_b])
            push!(rows, (gt.donor, prefix, g.gene, a))
        end
    end
    DataFrame(
        KI_case = first.(rows),
        locus = getindex.(rows, 2),
        gene = getindex.(rows, 3),
        sequence_id = last.(rows),
    )
end

"""Truth phase table: gene, locus, zygosity, chr_a, chr_b for heterozygous genes only."""
function truth_phase_table(gt::Genotype)
    rows = NamedTuple{(:gene, :locus, :zygosity, :chr_a, :chr_b), Tuple{String, String, String, String, String}}[]
    for g in gt.genes
        is_heterozygous(g) || continue
        locus_str = g.locus === :V ? "IGHV" : g.locus === :D ? "IGHD" : "IGHJ"
        push!(rows, (gene = g.gene, locus = locus_str, zygosity = "het", chr_a = g.allele_a, chr_b = g.allele_b))
    end
    isempty(rows) ? DataFrame() : DataFrame(rows)
end

"""Write calls DataFrame to TSV (MIAIRR-style)."""
function write_calls(df::DataFrame, path::String)
    CSV.write(path, df; delim = '\t')
end

"""Write genotype table (ground-truth alleles per donor) to TSV."""
function write_genotype(df::DataFrame, path::String)
    CSV.write(path, df; delim = '\t')
end

"""Write truth-phase table (phased alleles for heterozygous genes) to TSV."""
function write_truth_phase(df::DataFrame, path::String)
    CSV.write(path, df; delim = '\t')
end

"""
    write_simulation_output(calls_df, genotype_df, truth_phase_df; calls_path, genotype_path, truth_phase_path)

Write MIAIRR calls and/or ground-truth tables. Omit a path (e.g. `genotype_path = nothing`) to skip that file.
"""
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
