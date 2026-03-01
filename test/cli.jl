@testset "CLI" begin
    args = Dict{String, Any}(
        "n-reads" => 50,
        "min-unique-vdj" => 0,
        "seed" => 1,
        "output" => "calls_cli.tsv",
        "genotype" => "genotype_cli.tsv",
        "truth-phase" => "phase_cli.tsv",
        "subject" => "cli_donor",
        "p-v" => 0.05,
        "p-d" => 0.12,
        "p-j" => 0.03,
        "anchor-j-fraction" => "",
    )
    config = config_from_args(args)
    @test config.n_reads == 50
    @test config.output_path == "calls_cli.tsv"
    @test config.subject_id == "cli_donor"

    td = mktempdir()
    calls_p = joinpath(td, "calls.tsv")
    gen_p = joinpath(td, "genotype.tsv")
    phase_p = joinpath(td, "phase.tsv")
    config2 = SimulatorConfig(
        n_reads = 20,
        output_path = calls_p,
        subject_id = "cli_test",
    )
    estimated = ki_donor_preset()
    sim = Simulator(config2, estimated; rng = Random.MersenneTwister(0))
    calls_df, genotype_df, truth_phase_df = simulate(sim)
    write_simulation_output(calls_df, genotype_df, truth_phase_df; calls_path = calls_p, genotype_path = gen_p, truth_phase_path = phase_p)
    @test isfile(calls_p)
    @test isfile(gen_p)
    @test isfile(phase_p)
end
