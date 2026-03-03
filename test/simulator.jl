@testset "Simulator" begin
    config = SimulatorConfig(
        n_reads = 100,
        min_unique_vdj = 0,
        seed = 12345,
        output_path = "calls.tsv",
        subject_id = "sim_donor",
        noise = NoiseConfigNone,  # no noise so all calls are from genotype
    )
    estimated = ki_donor_preset()
    sim = Simulator(config, estimated; rng = Random.MersenneTwister(42))
    calls_df, genotype_df, truth_phase_df = simulate(sim)

    @test nrow(calls_df) == 100
    @test "subject" in names(calls_df)
    @test "v_call" in names(calls_df)
    @test "d_call" in names(calls_df)
    @test "j_call" in names(calls_df)
    @test "duplicate_count" in names(calls_df)
    @test "true_chr" in names(calls_df)
    @test all(calls_df.subject .== "sim_donor")
    @test all(calls_df.locus .== "IGH")

    @test nrow(genotype_df) >= 1
    @test "KI_case" in names(genotype_df)
    all_v = Set(calls_df.v_call)
    all_j = Set(calls_df.j_call)
    gt_alleles_v = Set(genotype_df[genotype_df.locus .== "IGHV", :].sequence_id)
    gt_alleles_j = Set(genotype_df[genotype_df.locus .== "IGHJ", :].sequence_id)
    @test all_v ⊆ gt_alleles_v
    @test all_j ⊆ gt_alleles_j

    out = sim()
    @test length(out) == 3
    @test out[1] isa DataFrame
end
