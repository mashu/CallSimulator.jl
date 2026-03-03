@testset "Genotype (DonorGenotype, GeneEntry{L})" begin
    g_hom = homozygous_gene("IGHV1-2", V, "IGHV1-2*01")
    @test g_hom isa GeneEntry{V}
    @test zygosity(g_hom) == Homozygous
    @test is_available(g_hom, 1) && is_available(g_hom, 2)
    @test allele_on(g_hom, 1) == "IGHV1-2*01"

    g_het = heterozygous_gene("IGHJ1", J, "IGHJ1*01", "IGHJ1*02")
    @test zygosity(g_het) == Heterozygous
    @test is_available(g_het, 1) && is_available(g_het, 2)
    @test allele_on(g_het, 2) == "IGHJ1*02"

    g_hemi = hemizygous_gene("IGHD1-1", D, "IGHD1-1*01"; on_chr=1)
    @test zygosity(g_hemi) == Hemizygous
    @test is_available(g_hemi, 1) && !is_available(g_hemi, 2)
    @test allele_on(g_hemi, 1) == "IGHD1-1*01"
    @test_throws Exception allele_on(g_hemi, 2)

    gt = build_genotype(
        "donor1",
        ZygositySpec(hom=1, het=1, hemi=0),
        ZygositySpec(hom=0, het=0, hemi=0),
        ZygositySpec(hom=1, het=1, hemi=0);
        pool_v = ["V1", "V2", "V3"],
        pool_d = ["D1"],
        pool_j = ["J1", "J2", "J3"],
        rng = Random.MersenneTwister(123),
    )
    @test gt isa DonorGenotype
    @test gt.donor_id == "donor1"
    @test length(genes(gt, V)) == 2
    @test length(genes(gt, J)) == 2
    @test has_locus(gt, V) && has_locus(gt, J)

    gt2 = build_donor_genotype("d2"; pool_v = ["V1", "V2", "V3"], pool_d = ["D1", "D2"], pool_j = ["J1", "J2"], n_j_het_range = (1, 1), rng = Random.MersenneTwister(456))
    @test length(gt2.genes_j) >= 1
    @test any(zygosity(g) == Heterozygous for g in gt2.genes_j)

    @test_throws ArgumentError build_donor_genotype("x"; pool_j = ["J1"], n_j_het_range = (0, 0), rng = Random.MersenneTwister(0))

    # Display: show(io, MIME("text/plain"), gt) and zygosity_counts
    io = IOBuffer()
    show(io, MIME("text/plain"), gt)
    s = String(take!(io))
    @test occursin("DonorGenotype", s)
    @test occursin("genes", s)
    c = zygosity_counts(gt, V)
    @test c.homozygous + c.heterozygous + c.hemizygous == length(gt.genes_v)
    show(io, MIME("text/plain"), gt.genes_v[1])
    @test occursin("GeneEntry", String(take!(io)))
end
