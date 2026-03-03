@testset "Phasing (ground truth vs simulator output)" begin
    gt = build_genotype(
        "phasing_test",
        ZygositySpec(hom=1, het=0, hemi=0),
        ZygositySpec(hom=0, het=0, hemi=0),
        ZygositySpec(hom=1, het=1, hemi=0);
        pool_v = ["IGHV1-2"],
        pool_d = String[],
        pool_j = ["IGHJ1", "IGHJ2"],
        rng = Random.MersenneTwister(0),
    )
    tp = truth_phase_table(gt)
    @test nrow(tp) >= 1
    @test "het" in tp.zygosity
    @test "chr_a" in names(tp)
    @test "chr_b" in names(tp)

    noise = NoiseConfigNone
    config = SimulatorConfig(
        n_reads = 500,
        min_unique_vdj = 0,
        seed = 777,
        output_path = "calls.tsv",
        noise = noise,
    )
    expression = build_expression(gt; rng = Random.MersenneTwister(777), method = UniformExpr())
    sim = Simulator(gt, config, expression, NoiseModel(noise), Random.MersenneTwister(777))
    calls_df, genotype_df, truth_phase_df = simulate(sim)

    @test truth_phase_df == tp
    v_ok = all(c -> c in genotype_df[genotype_df.locus .== "IGHV", :].sequence_id, calls_df.v_call)
    j_ok = all(c -> c in genotype_df[genotype_df.locus .== "IGHJ", :].sequence_id, calls_df.j_call)
    @test v_ok
    @test j_ok
    j_calls = unique(calls_df.j_call)
    # Both alleles of the heterozygous J must appear (which J is het depends on build order)
    het_j = tp[tp.zygosity .== "het", :]
    for r in eachrow(het_j)
        @test r.chr_a in j_calls
        @test r.chr_b in j_calls
    end
    @test all(calls_df.true_chr .>= 1)
    @test all(calls_df.noise_v .== "none")
end
