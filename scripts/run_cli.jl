#!/usr/bin/env julia
# Run CallSimulator from the command line. From package root:
#   julia --project=. scripts/run_cli.jl [--n-reads 10000] [-o calls.tsv] ...
# Or: julia -e 'using CallSimulator; CallSimulator.run_cli()' -- --n-reads 500 -o out.tsv

pkg_dir = joinpath(@__DIR__, "..")
pushfirst!(LOAD_PATH, joinpath(pkg_dir, "src"))
using CallSimulator
CallSimulator.run_cli()
