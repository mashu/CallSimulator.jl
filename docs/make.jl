import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter
using CallSimulator

makedocs(
    sitename = "CallSimulator",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    modules = [CallSimulator],
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Standalone usage" => "standalone.md",
        "Noise model" => "noise.md",
        "API reference" => "api.md",
    ],
)

# Prevent GitHub Pages from using Jekyll (serve Documenter's static HTML as-is)
build_dir = joinpath(@__DIR__, "build")
open(joinpath(build_dir, ".nojekyll"), "w") do io
end

deploydocs(
    repo = "github.com/mashu/CallSimulator.jl.git",
    devbranch = "main",
    push_preview = true,
)
