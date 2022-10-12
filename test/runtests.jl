using MolecularGraphKernels

MolecularGraphKernels.banner()

include.(
    [
        "ProductGraph.jl"
        "check_isom.jl"
        "graph_conversion.jl"
        "graph_products.jl"
        "graph_kernels.jl"
        "gram_matrix.jl"
        "visualization.jl"
        "maccs.jl"
        "misc.jl"
        "docs.jl"
    ]
)
