using Documenter
using MolecularGraphKernels

makedocs(
    root = joinpath(dirname(pathof(MolecularGraphKernels)), "..", "docs"),
    modules = [MolecularGraphKernels],
    sitename = "MolecularGraphKernels.jl",
    clean = true,
    pages = [
            "MolecularGraphKernels" => "index.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/MolecularGraphKernels.jl.git")
