@testset "Config" begin
    noise = NoiseConfig()
    @test noise.p_allele_v >= 0.0
    @test noise.p_d_dropout >= 0.0

    cfg = SimulatorConfig(n_reads = 1000, output_path = "out.tsv", noise = noise)
    @test cfg.n_reads == 1000
    @test cfg.min_unique_vdj == 0
    @test cfg.subject_id == "sim_donor"

    @test_throws ArgumentError SimulatorConfig(n_reads = 0, output_path = "x.tsv")
    @test CallSimulator.RABHIT_MIN_UNIQUE_VDJ == 2000
end
