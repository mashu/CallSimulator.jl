# ExpressionPlots.jl: UnicodePlots-based display of expression weights per gene.
# Requires: Types (Locus, DonorGenotype, genes), Expression (ExpressionProfile, weights, ChromosomeWeights, LocusWeights).

function expression_barplot(labels::Vector{String}, w::Vector{Float64}; title::String="Expression weights", width::Int=50)
    ord = sortperm(w; rev=true)
    barplot(
        labels[ord],
        w[ord];
        title,
        xlabel = "weight",
        width = min(width, max(20, length(labels) + 5)),
    )
end

"""
    plot_expression(ep::ExpressionProfile, L::Type{<:Locus}; chr::Int=1, gt::Union{DonorGenotype, Nothing}=nothing)

Bar plot of expression weights per gene for locus `L` (V, D, or J) on chromosome `chr`.
Bars are sorted by weight (highest first). Optional `gt` uses gene names as labels.
Returns the UnicodePlots plot (display in REPL automatically or with `println(...)`).
"""
function plot_expression(ep::ExpressionProfile, L::Type{<:Locus}; chr::Int=1, gt::Union{DonorGenotype, Nothing}=nothing)
    cw = weights(ep, L, chr)
    plot_expression_weights(cw, L, chr, gt)
end

"""
    plot_expression(lw::LocusWeights, L::Type{<:Locus}; chr::Int=1, gt::Union{DonorGenotype, Nothing}=nothing)

Bar plot of expression weights for this locus (e.g. `plot_expression(sim.expression.d, D)`).
Bars sorted by weight. Use `chr=2` for the second chromosome.
"""
function plot_expression(lw::LocusWeights, L::Type{<:Locus}; chr::Int=1, gt::Union{DonorGenotype, Nothing}=nothing)
    cw = chr == 1 ? lw.chr1 : lw.chr2
    plot_expression_weights(cw, L, chr, gt)
end

"""
    plot_expression(cw::ChromosomeWeights; title::String="Expression weights")

Bar plot of expression weights for one chromosome (e.g. `plot_expression(sim.expression.d.chr1)`).
Bars sorted by weight (highest first).
"""
function plot_expression(cw::ChromosomeWeights; title::String="Expression weights")
    isempty(cw) && error("No expression weights")
    indices = [first(p) for p in cw]
    w = [last(p) for p in cw]
    labels = string.(indices)
    expression_barplot(labels, w; title)
end

function plot_expression_weights(cw::ChromosomeWeights, L::Type{<:Locus}, chr::Int, gt::Union{DonorGenotype, Nothing})
    isempty(cw) && error("No expression weights for this locus/chromosome")
    indices = [first(p) for p in cw]
    w = [last(p) for p in cw]
    labels = if gt !== nothing
        gs = genes(gt, L)
        [i <= length(gs) ? gs[i].gene : string(i) for i in indices]
    else
        string.(indices)
    end
    locus_name = string(nameof(L))
    expression_barplot(labels, w; title = "Expression ($locus_name, chr$chr)")
end

function Base.show(io::IO, ::MIME"text/plain", cw::ChromosomeWeights)
    if isempty(cw)
        print(io, "ChromosomeWeights(0 genes)")
        return
    end
    p = plot_expression(cw; title = "Expression weights")
    print(io, string(p))
end

# So println(cw) and REPL display both show the bar chart
Base.show(io::IO, cw::ChromosomeWeights) = show(io, MIME("text/plain"), cw)
