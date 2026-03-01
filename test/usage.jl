@testset "Usage" begin
    g1 = Gene("V1", :V, "V1*01", "V1*02")
    g2 = Gene("J1", :J, "J1*01")
    gt = Genotype("donor", [g1, g2])

    u = uniform_usage(gt)
    @test length(u.v) == 2
    @test length(u.j) == 1
    @test all(p -> p.second == 1.0, u.v)

    rng = Random.MersenneTwister(1)
    u2 = skewed_usage(gt; rng = rng, lognormal_sigma = 0.5, allele_imbalance = 1.0)
    @test length(u2.v) == 2
    @test length(u2.j) == 1
    allele = u2(rng, :V, ["V1*01", "V1*02"])
    @test allele in ["V1*01", "V1*02"]
    allele_j = u2(rng, :J, ["J1*01"])
    @test allele_j == "J1*01"

    u3 = skewed_usage(gt; rng = Random.MersenneTwister(2), anchor_j_fraction = 0.3)
    @test length(u3.j) == 1

    # skewed_usage_kde stub when KernelDensity not loaded
    @test_throws ErrorException skewed_usage_kde(gt, [0.0, 0.5, 1.0])
end
