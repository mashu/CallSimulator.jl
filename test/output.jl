@testset "Output" begin
    sub = "subject1"
    v = ["V1*01", "V2*01"]
    d = ["D1*01", ""]
    j = ["J1*01", "J2*01"]
    ids = ["seq1", "seq2"]
    ds = [1, 1]
    de = [12, 0]
    dup = [1, 2]
    df = build_calls_df(sub, v, d, j, ids, ds, de, dup)
    @test nrow(df) == 2
    @test df.subject == fill("subject1", 2)
    @test df.sequence_id == ids
    @test df.v_call == v
    @test df.duplicate_count == dup
    @test all(df.locus .== "IGH")

    gt = build_genotype(
        "donor1",
        ZygositySpec(hom=1, het=0, hemi=0),
        ZygositySpec(hom=0, het=0, hemi=0),
        ZygositySpec(hom=1, het=1, hemi=0);
        pool_v = ["V1"],
        pool_d = ["D1"],
        pool_j = ["J1", "J2"],
        rng = Random.MersenneTwister(0),
    )
    g_df = genotype_table(gt)
    @test names(g_df) == ["KI_case", "locus", "gene", "sequence_id"]
    @test "donor1" in g_df.KI_case
    @test "IGHV" in g_df.locus
    @test "IGHJ" in g_df.locus

    tp_df = truth_phase_table(gt)
    @test names(tp_df) == ["gene", "locus", "zygosity", "chr_a", "chr_b"]
    @test nrow(tp_df) >= 1
    @test "het" in tp_df.zygosity || "hemi" in tp_df.zygosity

    gt_hom = build_genotype("d", ZygositySpec(1, 0, 0), ZygositySpec(0, 0, 0), ZygositySpec(1, 0, 0); pool_v = ["V1"], pool_d = String[], pool_j = ["J1"], rng = Random.MersenneTwister(0))
    @test nrow(truth_phase_table(gt_hom)) == 0

    td = mktempdir()
    calls_p = joinpath(td, "calls.tsv")
    gen_p = joinpath(td, "genotype.tsv")
    phase_p = joinpath(td, "phase.tsv")
    write_simulation_output(df, g_df, tp_df; calls_path = calls_p, genotype_path = gen_p, truth_phase_path = phase_p)
    @test isfile(calls_p)
    @test isfile(gen_p)
    @test isfile(phase_p)
end
