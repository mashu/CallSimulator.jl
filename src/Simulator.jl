# Simulator.jl: Read loop — pick chromosome, sample V/D/J from it, apply noise, record.

using Random

"""Genotype, config, expression, RNG. Call `simulate(sim)` for (calls_df, genotype_df, truth_phase_df)."""
struct Simulator
    genotype::DonorGenotype
    config::SimulatorConfig
    expression::ExpressionProfile
    rng::Random.AbstractRNG
    Simulator(gt, cfg, expr, rng) = new(gt, cfg, expr, rng)
end

function Base.show(io::IO, ::MIME"text/plain", sim::Simulator)
    nv = length(sim.genotype.genes_v)
    nd = length(sim.genotype.genes_d)
    nj = length(sim.genotype.genes_j)
    println(io, "Simulator:")
    println(io, "  genotype: ", sim.genotype.donor_id, " (", nv, " V, ", nd, " D, ", nj, " J genes)")
    print(io, "  config: ")
    show(io, MIME("text/plain"), sim.config)
    println(io)
    print(io, "  expression: ")
    show(io, MIME("text/plain"), sim.expression)
    println(io)
    print(io, "  rng: ", typeof(sim.rng))
end

"""Build Simulator from config and params; uses default gene pools. Shortcut for `Simulator(config, GenePools(), params)`."""
function Simulator(
    config::SimulatorConfig,
    params::EmpiricalParams;
    rng::Random.AbstractRNG = Random.MersenneTwister(config.seed),
)
    Simulator(config, GenePools(), params; rng = rng)
end

"""Build Simulator from config, gene pools, and params (genotype + expression + noise)."""
function Simulator(
    config::SimulatorConfig,
    gene_pools::GenePools,
    params::EmpiricalParams;
    rng::Random.AbstractRNG = Random.MersenneTwister(config.seed),
)
    gt = build_donor_genotype(
        config.subject_id;
        pool_v = gene_pools.v,
        pool_d = gene_pools.d,
        pool_j = gene_pools.j,
        p_hemi_v = params.p_hemi_v,
        p_hemi_d = params.p_hemi_d,
        p_hemi_j = params.p_hemi_j,
        rng = Random.MersenneTwister(config.seed + 1),
    )
    expression = build_expression(gt;
        rng = Random.MersenneTwister(config.seed + 2),
        method = LogNormalExpr(params.lognormal_sigma),
        allele_imbalance_range = (1.0 / params.allele_imbalance, params.allele_imbalance),
        anchor_j_fraction = params.anchor_j_fraction_range,
    )
    Simulator(gt, config, expression, rng)
end

"""Sample one locus for this chromosome: expression → true allele → noise → (call, true_allele, noise_label)."""
function sample_locus(sim::Simulator, gt::DonorGenotype, ::Type{L}, chr::Int) where L<:Locus
    rng = sim.rng
    idx = sim.expression(rng, gt, L, chr)
    g = genes(gt, L)[idx]
    true_allele = allele_on(g, chr)
    call, nt = NoiseModel(sim.config.noise)(rng, gt, L, chr, g, true_allele)
    (call, true_allele, noise_type_string(nt))
end

function simulate(sim::Simulator)
    gt, rng, cfg = sim.genotype, sim.rng, sim.config
    n, min_uv = cfg.n_reads, cfg.min_unique_vdj

    has_locus(gt, V) || throw(ArgumentError("Genotype must have at least one V gene"))
    has_locus(gt, J) || throw(ArgumentError("Genotype must have at least one J gene"))

    v_calls, d_calls, j_calls = String[], String[], String[]
    sequence_ids = String[]
    d_starts, d_ends = Int[], Int[]
    duplicate_counts = Int[]
    true_chr = Int[]
    true_v, true_d, true_j = String[], String[], String[]
    noise_v, noise_d, noise_j = String[], String[], String[]

    read_idx = 0
    max_reads = min_uv > 0 ? max(n, 500_000) : n
    unique_vdj = Set{Tuple{String, String, String}}()

    while true
        read_idx += 1
        chr = rand(rng, 1:2)

        v_call, tv, nv = sample_locus(sim, gt, V, chr)
        push!(noise_v, nv)

        if has_locus(gt, D)
            d_call, td, nd = sample_locus(sim, gt, D, chr)
            push!(true_d, td)
            push!(noise_d, nd)
        else
            d_call = ""
            push!(true_d, "")
            push!(noise_d, "none")
        end

        j_call, tj, nj = sample_locus(sim, gt, J, chr)
        push!(noise_j, nj)

        push!(v_calls, v_call)
        push!(d_calls, d_call)
        push!(j_calls, j_call)
        push!(sequence_ids, "sim_$(read_idx)")
        push!(d_starts, 1)
        push!(d_ends, isempty(d_call) ? 0 : 12)
        push!(duplicate_counts, 1)
        push!(true_chr, chr)
        push!(true_v, tv)
        push!(true_j, tj)
        push!(unique_vdj, (v_calls[end], d_calls[end], j_calls[end]))

        if read_idx >= n && (min_uv <= 0 || length(unique_vdj) >= min_uv)
            break
        end
        if read_idx >= max_reads
            min_uv > 0 && length(unique_vdj) < min_uv &&
                @warn "min_unique_vdj=$min_uv not reached after $max_reads reads (unique=$(length(unique_vdj)))"
            break
        end
    end

    calls_df = build_calls_df(
        cfg.subject_id, v_calls, d_calls, j_calls,
        sequence_ids, d_starts, d_ends, duplicate_counts;
        true_chr, true_v, true_d, true_j, noise_v, noise_d, noise_j,
    )
    (calls_df, genotype_table(gt), truth_phase_table(gt))
end

(sim::Simulator)(rng::Random.AbstractRNG) = simulate(Simulator(sim.genotype, sim.config, sim.expression, rng))
(sim::Simulator)() = simulate(sim)
