# Simulator.jl: Main simulation loop. Type-stable; uses Config, Empirical, Genotype, Usage, MissCall.

const RABHIT_MIN_UNIQUE_VDJ = 2000

"""
    CallSimulator

Holds genotype, user config, empirical params, usage model, miss-call model, and RNG.
Construct via `CallSimulator(config, estimated)` or use the functor `(sim)(rng)`.
"""
struct Simulator
    genotype::Genotype
    config::SimulatorConfig
    usage::WeightedUsage
    miss_call::BernoulliMissCall
    rng::Random.AbstractRNG
end

"""Build Simulator from config and estimated params. Builds genotype and usage from empirical."""
function Simulator(
    config::SimulatorConfig,
    estimated::EmpiricalParams;
    rng::Random.AbstractRNG = Random.MersenneTwister(config.seed),
)
    gt = build_genotype_realistic(
        config.subject_id,
        estimated.gene_pool_v,
        estimated.gene_pool_d,
        estimated.gene_pool_j;
        rng = Random.MersenneTwister(config.seed + 1),
    )
    anchor_frac = config.anchor_j_fraction !== nothing ? config.anchor_j_fraction : estimated.anchor_j_fraction_default
    usage = skewed_usage(gt;
        rng = Random.MersenneTwister(config.seed + 2),
        lognormal_sigma = estimated.lognormal_sigma,
        allele_imbalance = estimated.allele_imbalance,
        anchor_j_fraction = anchor_frac,
    )
    miss_call = BernoulliMissCall(config.miss_call)
    Simulator(gt, config, usage, miss_call, rng)
end

"""Return (allele, gene_idx) for one locus and chromosome."""
function _sample_allele(
    rng::Random.AbstractRNG,
    usage::WeightedUsage,
    locus::Symbol,
    list::Vector{Tuple{Int, String, String}},
    chr::Int,
)
    chr_alleles = [chr == 1 ? x[2] : x[3] for x in list]
    allele = usage(rng, locus, chr_alleles)
    for (k, t) in enumerate(list)
        (chr == 1 ? t[2] : t[3]) == allele && return (allele, t[1])
    end
    (allele, list[1][1])
end

"""Apply miss-call: return called allele (other_allele if swap, else true_allele)."""
function _apply_miss(
    rng::Random.AbstractRNG,
    miss::BernoulliMissCall,
    true_allele::String,
    other_allele::String,
    locus::Symbol,
    gene_name::String,
)
    miss(rng, true_allele, other_allele, locus, gene_name)
end

"""
    simulate(sim::Simulator) -> (calls_df, genotype_df, truth_phase_df)

Generate call table from genotype with usage and miss-call. At least `config.n_reads` rows;
if `config.min_unique_vdj > 0`, continues until that many unique (v_call, d_call, j_call).
Returns: calls DataFrame (subject, v_call, d_call, j_call, sequence_id, ...), genotype table, truth phase.
"""
function simulate(sim::Simulator)
    gt = sim.genotype
    rng = sim.rng
    cfg = sim.config
    n = cfg.n_reads
    min_uv = cfg.min_unique_vdj
    v_list, d_list, j_list = locus_alleles(gt)
    isempty(v_list) && throw(ArgumentError("Genotype must have at least one V gene"))
    isempty(j_list) && throw(ArgumentError("Genotype must have at least one J gene"))

    v_calls = String[]
    d_calls = String[]
    j_calls = String[]
    sequence_ids = String[]
    d_germline_starts = Int[]
    d_germline_ends = Int[]

    read_idx = 0
    max_reads = min_uv > 0 ? max(n, 500_000) : n
    unique_vdj = Set{Tuple{String, String, String}}()

    while true
        read_idx += 1
        chr = rand(rng) < 0.5 ? 1 : 2

        allele_v, gene_idx_v = _sample_allele(rng, sim.usage, :V, v_list, chr)
        g_v = gt.genes[gene_idx_v]
        other_v = g_v.allele_a == allele_v ? g_v.allele_b : g_v.allele_a
        allele_v = _apply_miss(rng, sim.miss_call, allele_v, other_v, :V, g_v.gene)

        allele_d = ""
        if !isempty(d_list)
            allele_d, idx_d = _sample_allele(rng, sim.usage, :D, d_list, chr)
            g_d = gt.genes[idx_d]
            other_d = g_d.allele_a == allele_d ? g_d.allele_b : g_d.allele_a
            allele_d = _apply_miss(rng, sim.miss_call, allele_d, other_d, :D, g_d.gene)
        end

        allele_j, idx_j = _sample_allele(rng, sim.usage, :J, j_list, chr)
        g_j = gt.genes[idx_j]
        other_j = g_j.allele_a == allele_j ? g_j.allele_b : g_j.allele_a
        allele_j = _apply_miss(rng, sim.miss_call, allele_j, other_j, :J, g_j.gene)

        push!(v_calls, allele_v)
        push!(d_calls, allele_d)
        push!(j_calls, allele_j)
        push!(sequence_ids, "sim_$(read_idx)")
        push!(d_germline_starts, 1)
        push!(d_germline_ends, isempty(allele_d) ? 0 : 12)
        push!(unique_vdj, (allele_v, allele_d, allele_j))

        if read_idx >= n && (min_uv <= 0 || length(unique_vdj) >= min_uv)
            break
        end
        if read_idx >= max_reads
            min_uv > 0 && length(unique_vdj) < min_uv &&
                @warn "min_unique_vdj=$min_uv not reached after $max_reads reads (unique VDJ=$(length(unique_vdj)))"
            break
        end
    end

    n_out = length(v_calls)
    calls_df = build_calls_df(
        cfg.subject_id,
        v_calls,
        d_calls,
        j_calls,
        sequence_ids,
        d_germline_starts,
        d_germline_ends,
    )
    genotype_df = genotype_table(gt)
    truth_phase_df = truth_phase_table(gt)
    (calls_df, genotype_df, truth_phase_df)
end

"""Functors: (sim::Simulator)(rng) -> same as simulate(sim) (uses sim.rng if rng not provided)."""
(sim::Simulator)(rng::Random.AbstractRNG) = simulate(Simulator(sim.genotype, sim.config, sim.usage, sim.miss_call, rng))
(sim::Simulator)() = simulate(sim)
