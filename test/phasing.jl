@testset "Phasing (ground truth vs simulator output)" begin
    # Build a small genotype with known het genes and phase
    genes = [
        Gene("IGHV1-2", :V, "IGHV1-2*01", "IGHV1-2*02"),
        Gene("IGHJ1", :J, "IGHJ1*01", "IGHJ1*02"),
        Gene("IGHJ2", :J, "IGHJ2*01"),
    ]
    gt = Genotype("phasing_test", genes)

    tp = truth_phase_table(gt)
    @test nrow(tp) == 2
    j_phase = tp[tp.gene .== "IGHJ1", :][1, :]
    @test j_phase.chr_a == "IGHJ1*01"
    @test j_phase.chr_b == "IGHJ1*02"

    # Run simulator with zero miss-call so calls reflect true phase
    miss = MissCallConfig(0.0, 0.0, 0.0)
    config = SimulatorConfig(
        n_reads = 500,
        min_unique_vdj = 0,
        seed = 777,
        output_path = "calls.tsv",
        miss_call = miss,
    )
    usage = uniform_usage(gt)
    sim = Simulator(gt, config, usage, BernoulliMissCall(miss), Random.MersenneTwister(777))
    calls_df, genotype_df, truth_phase_df = simulate(sim)

    @test truth_phase_df == tp
    # All V and J calls must be from genotype
    v_ok = all(c -> c in genotype_df[genotype_df.locus .== "IGHV", :].sequence_id, calls_df.v_call)
    j_ok = all(c -> c in genotype_df[genotype_df.locus .== "IGHJ", :].sequence_id, calls_df.j_call)
    @test v_ok
    @test j_ok
    # With zero miss-call, each read comes from one chromosome; J calls should include both alleles for het J
    j_calls = unique(calls_df.j_call)
    @test "IGHJ1*01" in j_calls
    @test "IGHJ1*02" in j_calls
end
