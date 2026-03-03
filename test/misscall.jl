# Noise model (NoiseModel functor, NoiseConfig).
@testset "Noise (NoiseModel, NoiseConfig)" begin
    cfg = NoiseConfigNone
    noise = NoiseModel(cfg)
    rng = Random.MersenneTwister(99)
    gt = CallSimulator.build_genotype(
        "t",
        ZygositySpec(1, 1, 0),
        ZygositySpec(0, 0, 0),
        ZygositySpec(1, 0, 0);
        pool_v = ["V1", "V2"],
        pool_d = String[],
        pool_j = ["J1"],
        rng = Random.MersenneTwister(0),
    )
    g_v = gt.genes_v[2]
    @test zygosity(g_v) == Heterozygous
    called, nt = noise(rng, gt, V, 1, g_v, allele_on(g_v, 1))
    @test nt == NoNoise
    @test called == allele_on(g_v, 1)

    cfg2 = NoiseConfig(p_allele_v = 1.0, p_allele_d = 0.0, p_allele_j = 0.0, p_d_dropout = 0.0)
    noise2 = NoiseModel(cfg2)
    called2, nt2 = noise2(rng, gt, V, 1, g_v, allele_on(g_v, 1))
    @test nt2 == AlleleSwap
    @test called2 == allele_on(g_v, 2)

    @test noise_type_string(NoNoise) == "none"
    @test noise_type_string(AlleleSwap) == "allele_swap"
    @test noise_type_string(D_Dropout) == "dropout"
end
