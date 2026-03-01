@testset "Output" begin
    sub = "subject1"
    v = ["V1*01", "V2*01"]
    d = ["D1*01", ""]
    j = ["J1*01", "J2*01"]
    ids = ["seq1", "seq2"]
    ds = [1, 1]
    de = [12, 0]
    df = build_calls_df(sub, v, d, j, ids, ds, de)
    @test nrow(df) == 2
    @test df.subject == fill("subject1", 2)
    @test df.sequence_id == ids
    @test df.v_call == v
    @test df.d_call == d
    @test df.j_call == j
    @test df.d_germline_start == ds
    @test df.d_germline_end == de
    @test all(df.locus .== "IGH")
    @test all(df.productive)

    gt = Genotype("donor1", [
        Gene("V1", :V, "V1*01", "V1*02"),
        Gene("J1", :J, "J1*01"),
    ])
    g_df = genotype_table(gt)
    @test names(g_df) == ["KI_case", "locus", "gene", "sequence_id"]
    @test "donor1" in g_df.KI_case
    @test "IGHV" in g_df.locus
    @test "IGHJ" in g_df.locus

    tp_df = truth_phase_table(gt)
    @test names(tp_df) == ["gene", "locus", "zygosity", "chr_a", "chr_b"]
    @test nrow(tp_df) == 1
    @test tp_df.gene[1] == "V1"
    @test tp_df.chr_a[1] == "V1*01"
    @test tp_df.chr_b[1] == "V1*02"

    gt_hom = Genotype("d", [Gene("J1", :J, "J1*01")])
    @test nrow(truth_phase_table(gt_hom)) == 0

    # write_simulation_output
    td = mktempdir()
    calls_p = joinpath(td, "calls.tsv")
    gen_p = joinpath(td, "genotype.tsv")
    phase_p = joinpath(td, "phase.tsv")
    write_simulation_output(df, g_df, tp_df; calls_path = calls_p, genotype_path = gen_p, truth_phase_path = phase_p)
    @test isfile(calls_p)
    @test isfile(gen_p)
    @test isfile(phase_p)
    write_simulation_output(df, g_df, tp_df; calls_path = nothing, genotype_path = nothing, truth_phase_path = nothing)
    # no files created for that call
end
