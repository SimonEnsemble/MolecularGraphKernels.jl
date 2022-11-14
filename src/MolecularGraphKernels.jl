module MolecularGraphKernels

using Base.Threads,
    Cairo,
    Colors,
    Compose,
    Distributed,
    FIGlet,
    GraphPlot,
    Graphs,
    MetaGraphs,
    MolecularGraph,
    PeriodicTable,
    RDKitMinimalLib,
    SharedArrays,
    SparseArrays
using Graphs.Experimental: vf2, IsomorphismProblem
using PrecompileSignatures: @precompile_signatures

import Base: display, size
import Graphs: is_directed, SimpleEdge
import MetaGraphs: weighttype, PropDict, MetaDict, set_props!, props
import SparseArrays: AbstractSparseMatrixCSC, _checkbuffers, getcolptr, rowvals, nonzeros

function __init__()
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
        "graph_products.jl"
        "kernels.jl"
        "visualization.jl"
        "graph_conversion.jl"
        "check_isom.jl"
        "misc.jl"
        "maccs.jl"
        "gram_matrix.jl"
    ]
)

export

    # ProductGraph.jl
    ProductGraph,
    is_directed,
    weighttype,
    props,
    set_props!,

    # graph_products.jl
    product_graph_adjacency_matrix,
    Modular,
    Direct,
    Weighted,

    # graph_kernels.jl
    random_walk,
    ccsi,

    # graph_conversion.jl
    MetaGraph,

    # visualization.jl
    viz_graph,

    # check_isom.jl
    is_isomorphic,

    # gram_matrix.jl
    gram_matrix,
    gm_norm,

    # misc.jl
    display,
    isomorphic_subgraphs

@precompile_signatures(MolecularGraphKernels)

end
