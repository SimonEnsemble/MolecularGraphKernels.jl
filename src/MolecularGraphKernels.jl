module MolecularGraphKernels

using Cairo,
    Colors,
    Combinatorics,
    Compose,
    Distributed,
    GraphPlot,
    Graphs,
    JLD2,
    Memoization,
    MetaGraphs,
    MolecularGraph,
    PeriodicTable,
    ProgressMeter,
    RDKitMinimalLib
import PrecompileSignatures: @precompile_signatures

function __init__()
    # RDKitMinimalLib doesn't work on Windows 
    # https://github.com/eloyfelix/RDKitMinimalLib.jl/issues/13
    if !Sys.iswindows()
        # pre-compute the MACCS queries. must ignore "?". replace with nothing.
        global _maccs_queries = [
            (smarts_pattern == "?" ? nothing : get_qmol(smarts_pattern), nb_matches) for
            (smarts_pattern, nb_matches) in maccs_queries
        ]
    end
end

include.(
    [
        "ProductGraph.jl"
        "GraphMatrix.jl"
        "graph_products.jl"
        "kernels.jl"
        "visualization.jl"
        "graph_conversion.jl"
        "misc.jl"
        "maccs.jl"
        "gram_matrix.jl"
    ]
)

@precompile_signatures(MolecularGraphKernels)

end
