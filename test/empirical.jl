@testset "Empirical" begin
    ep = EmpiricalParams()
    @test length(ep.gene_pool_v) >= 30
    @test length(ep.gene_pool_d) >= 20
    @test length(ep.gene_pool_j) == 6
    @test ep.lognormal_sigma >= 0.0
    @test ep.allele_imbalance >= 1.0
    @test 0.0 < ep.anchor_j_fraction_default < 1.0

    ki = ki_donor_preset()
    @test ki.gene_pool_v == ep.gene_pool_v
    @test ki.gene_pool_d == ep.gene_pool_d
    @test ki.gene_pool_j == ep.gene_pool_j

    ki2 = ki_donor_preset(lognormal_sigma = 0.5, anchor_j_fraction_default = 0.25)
    @test ki2.lognormal_sigma == 0.5
    @test ki2.anchor_j_fraction_default == 0.25

    @test_throws ArgumentError EmpiricalParams(lognormal_sigma = -0.1)
    @test_throws ArgumentError EmpiricalParams(allele_imbalance = 0.5)
    @test_throws ArgumentError EmpiricalParams(anchor_j_fraction_default = 0.0)
end
