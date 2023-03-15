import MolecularGraphKernels: banner

banner()

include.(
    [
        "init.jl"
        "ProductGraph.jl"
        "check_isom.jl"
        "graph_conversion.jl"
        "graph_products.jl"
        "con_sub_g.jl"
        "graph_kernels.jl"
        "gram_matrix.jl"
        "visualization.jl"
        "maccs.jl"
        "misc.jl"
        "docs.jl"
    ]
)
