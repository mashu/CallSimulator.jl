# Validation.jl: compare simulated vs real AIRR tables (marginals, summaries).

const REQUIRED_AIRR_COLUMNS = (
    "sequence_id",
    "v_call",
    "d_call",
    "j_call",
)

"""Check whether a DataFrame has core AIRR call columns."""
function check_airr_columns(df::DataFrame)
    missing = [c for c in REQUIRED_AIRR_COLUMNS if !(c in names(df))]
    (isempty(missing), missing)
end

"""Normalized count distribution (value -> proportion)."""
function norm_counts(vals)
    c = Dict{String, Int}()
    for v in vals
        s = String(v)
        c[s] = get(c, s, 0) + 1
    end
    n = sum(values(c))
    n == 0 && return Dict{String, Float64}()
    Dict(k => v / n for (k, v) in c)
end

"""
    gini(values) -> Float64

Gini coefficient for a distribution of non-negative values (e.g. gene usage counts).
Measures inequality: 0 = all values equal (uniform usage), 1 = maximally unequal (one category
gets everything). For V/D/J call counts across genes, higher Gini means more skewed usage.
Real repertoires often have Gini in roughly 0.3–0.7 for V (many genes, skewed); J can be lower
when there are few genes. Simulated data with `LogNormalExpr` typically give similar ranges.
"""
function gini(values::Vector{Float64})
    isempty(values) && return 0.0
    s = sum(values)
    s <= 0.0 && return 0.0
    x = sort(values)
    n = length(x)
    accum = 0.0
    for i in 1:n
        accum += (2i - n - 1) * x[i]
    end
    accum / (n * s)
end

"""Total variation distance between two probability distributions."""
function total_variation_distance(a::Dict{String, Float64}, b::Dict{String, Float64})
    keys_all = union(keys(a), keys(b))
    s = 0.0
    for k in keys_all
        s += abs(get(a, k, 0.0) - get(b, k, 0.0))
    end
    0.5 * s
end

"""Summary stats for calibration (n_reads, unique_vdj, d_dropout_rate, gini per locus, etc.)."""
function summarize_simulation(df::DataFrame)
    n = nrow(df)
    dists = [norm_counts(df[!, c]) for c in [:v_call, :d_call, :j_call]]
    ginis = [gini(collect(values(d))) for d in dists]
    has_noise = "noise_d" in names(df)
    d_dropout = has_noise ? sum(df.noise_d .== "dropout") / n : sum(df.d_call .== "") / n
    has_dup = "duplicate_count" in names(df)
    mean_dup = has_dup ? sum(df.duplicate_count) / n : 1.0
    has_truth = "true_chr" in names(df)
    chr_balance = has_truth ? sum(df.true_chr .== 1) / n : NaN

    (
        n_reads = n,
        unique_vdj = length(Set(zip(df.v_call, df.d_call, df.j_call))),
        d_dropout_rate = d_dropout,
        mean_duplicate_count = mean_dup,
        chr_balance = chr_balance,
        gini_v = ginis[1], gini_d = ginis[2], gini_j = ginis[3],
    )
end

"""Total variation distance between simulated and real marginals for V, D, J."""
function compare_with_real(sim_df::DataFrame, real_df::DataFrame)
    ok1, miss1 = check_airr_columns(sim_df)
    ok2, miss2 = check_airr_columns(real_df)
    ok1 || throw(ArgumentError("sim_df missing: " * join(miss1, ", ")))
    ok2 || throw(ArgumentError("real_df missing: " * join(miss2, ", ")))

    cols = [:v_call, :d_call, :j_call]
    values = [total_variation_distance(norm_counts(sim_df[!, c]), norm_counts(real_df[!, c])) for c in cols]
    DataFrame(metric = ["tv_v", "tv_d", "tv_j"], value = values)
end
