@testset "Genotype" begin
    g_hom = Gene("IGHV1-2", :V, "IGHV1-2*01")
    @test g_hom.allele_a == g_hom.allele_b
    @test !is_heterozygous(g_hom)

    g_het = Gene("IGHJ1", :J, "IGHJ1*01", "IGHJ1*02")
    @test is_heterozygous(g_het)

    gt = Genotype("donor1", [g_hom, g_het])
    @test length(gt) == 2
    @test allele_on_chr(gt, 1, 1) == "IGHV1-2*01"
    @test allele_on_chr(gt, 2, 2) == "IGHJ1*02"
    @test locus_of(gt, 1) === :V
    @test locus_of(gt, 2) === :J

    v_list, d_list, j_list = locus_alleles(gt)
    @test length(v_list) == 1
    @test length(d_list) == 0
    @test length(j_list) == 1
    @test v_list[1] == (1, "IGHV1-2*01", "IGHV1-2*01")
    @test j_list[1][2] != j_list[1][3]

    pool_v = ["V1", "V2", "V3"]
    pool_d = ["D1", "D2"]
    pool_j = ["J1", "J2"]
    gt2 = build_genotype("d2", 1, 1, 0, 1, 1, 0; gene_pool_v = pool_v, gene_pool_d = pool_d, gene_pool_j = pool_j, rng = Random.MersenneTwister(123))
    @test gt2.donor == "d2"
    n_genes = length(gt2.genes)
    @test n_genes >= 4
    n_j_het = count(g -> g.locus === :J && is_heterozygous(g), gt2.genes)
    @test n_j_het == 1

    gt3 = build_genotype_realistic("d3", pool_v, pool_d, pool_j; n_j_het_range = (1, 1), n_v_het_range = (1, 1), rng = Random.MersenneTwister(456))
    @test length(gt3.genes) >= 4
    j_genes = [g for g in gt3.genes if g.locus === :J]
    @test any(is_heterozygous, j_genes)

    @test_throws ArgumentError build_genotype_realistic("x", ["V1"], ["D1"], ["J1"]; n_j_het_range = (0, 0), rng = Random.MersenneTwister(0))
end
