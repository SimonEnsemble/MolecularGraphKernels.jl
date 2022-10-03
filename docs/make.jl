using Documenter
using MolecularGraphKernels

makedocs(
    root = joinpath(dirname(pathof(MolecularGraphKernels)), "..", "docs"),
    modules = [MolecularGraphKernels],
    sitename = "MolecularGraphKernels.jl",
    clean = true,
    pages = [
            "MolecularGraphKernels" => "index.md"
            "Molecular Graphs" => "molecular_graphs.md"
            "Product Graphs" => "product_graphs.md"
            "Graph Kernels" => "graph_kernels.md"
            "Isomorphism" => "isomorphism.md"
            "Visualization" => "visualization.md"
            "API" => "api.md"
            ],
    format = Documenter.HTML(assets = ["assets/flux.css"])
)

deploydocs(repo = "github.com/SimonEnsemble/MolecularGraphKernels.jl.git")
