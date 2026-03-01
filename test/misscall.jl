@testset "MissCall" begin
    cfg = MissCallConfig(0.0, 0.0, 0.0)
    m = BernoulliMissCall(cfg)
    rng = Random.MersenneTwister(99)
    @test m(rng, "A*01", "A*02", :V, "V1") == "A*01"
    @test m(rng, "A*01", "A*02", :J, "J1") == "A*01"

    cfg2 = MissCallConfig(1.0, 1.0, 1.0)
    m2 = BernoulliMissCall(cfg2)
    @test m2(rng, "A*01", "A*02", :V, "V1") == "A*02"
    @test m2(rng, "A*01", "A*02", :D, "D1") == "A*02"

    @test miss_prob(m, :V, "V1") == 0.0
    @test miss_prob(m2, :J, "J1") == 1.0
end
