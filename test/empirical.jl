@testset "Empirical" begin
    ep = EmpiricalParams()
    @test ep.lognormal_sigma >= 0.0
    @test ep.allele_imbalance >= 1.0
    @test 0.0 < first(ep.anchor_j_fraction_range) < last(ep.anchor_j_fraction_range) < 1.0
    @test 0.0 <= ep.p_hemi_v <= 1.0
    @test 0.0 <= ep.p_hemi_d <= 1.0
    @test 0.0 <= ep.p_hemi_j <= 1.0

    pools = GenePools()
    @test length(pools.v) >= 30
    @test length(pools.d) >= 20
    @test length(pools.j) == 6

    ki2 = ki_donor_preset(lognormal_sigma = 0.5, anchor_j_fraction_range = (0.15, 0.35))
    @test ki2.lognormal_sigma == 0.5
    @test ki2.anchor_j_fraction_range == (0.15, 0.35)

    @test_throws ArgumentError EmpiricalParams(lognormal_sigma = -0.1)
    @test_throws ArgumentError EmpiricalParams(allele_imbalance = 0.5)
    @test_throws ArgumentError EmpiricalParams(anchor_j_fraction_range = (0.0, 0.5))
end
