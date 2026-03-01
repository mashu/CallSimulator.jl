@testset "Config" begin
    mc = MissCallConfig(0.05, 0.12, 0.03)
    @test mc.p_v == 0.05
    @test mc.p_d == 0.12
    @test mc.p_j == 0.03
    @test isempty(mc.per_gene_v)

    @test_throws ArgumentError MissCallConfig(-0.1, 0.1, 0.1)
    @test_throws ArgumentError MissCallConfig(0.1, 1.5, 0.1)

    cfg = SimulatorConfig(n_reads = 1000, output_path = "out.tsv", miss_call = mc)
    @test cfg.n_reads == 1000
    @test cfg.min_unique_vdj == 0
    @test cfg.subject_id == "sim_donor"
    @test cfg.anchor_j_fraction === nothing

    cfg2 = SimulatorConfig(n_reads = 500, output_path = "c.tsv", anchor_j_fraction = 0.3)
    @test cfg2.anchor_j_fraction == 0.3

    @test_throws ArgumentError SimulatorConfig(n_reads = 0, output_path = "x.tsv")
    @test_throws ArgumentError SimulatorConfig(n_reads = 10, output_path = "x.tsv", anchor_j_fraction = 1.0)
end
