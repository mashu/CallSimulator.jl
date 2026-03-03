# CallSimulator test suite: full coverage for simulator logic.

using Test
using Random
using DataFrames
using CallSimulator

include("config.jl")
include("empirical.jl")
include("genotype.jl")
include("usage.jl")
include("misscall.jl")
include("output.jl")
include("simulator.jl")
include("phasing.jl")
