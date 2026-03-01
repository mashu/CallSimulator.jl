# CLI.jl: ArgParse-based command-line interface. No try-catch; validation via ArgumentError.

function parse_commandline()
    s = ArgParse.ArgParseSettings(
        prog = "CallSimulator",
        description = "Simulate V/D/J call tables (MIAIRR-style) with configurable miss-call and usage.",
    )
    ArgParse.@add_arg_table! s begin
        "--n-reads"
            help = "Number of reads (calls) to generate"
            arg_type = Int
            default = 10_000
        "--min-unique-vdj"
            help = "Minimum unique (v,d,j) combinations (0 = disable; RAbHIT recommends 2000)"
            arg_type = Int
            default = 0
        "--seed"
            help = "RNG seed"
            arg_type = Int
            default = 42
        "--output", "-o"
            help = "Output TSV path for MIAIRR-style calls"
            arg_type = String
            default = "calls.tsv"
        "--genotype"
            help = "Optional path to save ground-truth genotype table (TSV)"
            arg_type = String
            default = ""
        "--truth-phase"
            help = "Optional path to save ground-truth phased genotype (heterozygous genes only, TSV)"
            arg_type = String
            default = ""
        "--subject"
            help = "Subject/donor identifier"
            arg_type = String
            default = "sim_donor"
        "--p-v"
            help = "Miss-call probability for V"
            arg_type = Float64
            default = 0.05
        "--p-d"
            help = "Miss-call probability for D"
            arg_type = Float64
            default = 0.12
        "--p-j"
            help = "Miss-call probability for J"
            arg_type = Float64
            default = 0.03
        "--anchor-j-fraction"
            help = "Fraction of reads carrying anchor J (empty = use default from empirical)"
            arg_type = String
            default = ""
    end
    ArgParse.parse_args(s)
end

"""Build SimulatorConfig from parsed CLI args."""
function config_from_args(args::Dict{String, T}) where T
    anchor = get(args, "anchor-j-fraction", "")
    anchor_frac = isempty(anchor) ? nothing : parse(Float64, anchor)
    if anchor_frac !== nothing
        (anchor_frac > 0.0 && anchor_frac < 1.0) || throw(ArgumentError("--anchor-j-fraction must be in (0, 1)"))
    end
    miss = MissCallConfig(
        args["p-v"],
        args["p-d"],
        args["p-j"],
    )
    SimulatorConfig(
        n_reads = args["n-reads"],
        min_unique_vdj = args["min-unique-vdj"],
        seed = args["seed"],
        output_path = args["output"],
        subject_id = args["subject"],
        miss_call = miss,
        anchor_j_fraction = anchor_frac,
    )
end

"""Run simulator from CLI: parse args, build config and KI donor preset, simulate, write MIAIRR + optional ground truth."""
function run_cli()
    args = parse_commandline()
    config = config_from_args(args)
    estimated = ki_donor_preset()
    sim = Simulator(config, estimated)
    calls_df, genotype_df, truth_phase_df = simulate(sim)
    genotype_path = get(args, "genotype", "")
    truth_phase_path = get(args, "truth-phase", "")
    write_simulation_output(
        calls_df, genotype_df, truth_phase_df;
        calls_path = config.output_path,
        genotype_path = isempty(genotype_path) ? nothing : genotype_path,
        truth_phase_path = isempty(truth_phase_path) ? nothing : truth_phase_path,
    )
    @info "Wrote output" calls = config.output_path genotype = genotype_path truth_phase = truth_phase_path n_rows = nrow(calls_df)
end
