# Expression profile replaces old WeightedUsage; test build_expression and sampling.
@testset "Expression (build_expression, ExpressionProfile functor)" begin
    gt = CallSimulator.build_genotype(
        "donor",
        ZygositySpec(hom=1, het=1, hemi=0),
        ZygositySpec(hom=0, het=0, hemi=0),
        ZygositySpec(hom=1, het=0, hemi=0);
        pool_v = ["V1", "V2"],
        pool_d = String[],
        pool_j = ["J1"],
        rng = Random.MersenneTwister(1),
    )
    ep = build_expression(gt; rng = Random.MersenneTwister(1), method = UniformExpr())
    @test ep isa ExpressionProfile
    @test !isempty(weights(ep, V, 1))
    @test !isempty(weights(ep, J, 1))

    rng = Random.MersenneTwister(2)
    idx_v = ep(rng, gt, V, 1)
    idx_j = ep(rng, gt, J, 1)
    @test 1 <= idx_v <= length(gt.genes_v)
    @test 1 <= idx_j <= length(gt.genes_j)

    ep2 = build_expression(gt; rng = Random.MersenneTwister(3), method = LogNormalExpr(0.5), anchor_j_fraction = (0.2, 0.4))
    @test ep2 isa ExpressionProfile
end
