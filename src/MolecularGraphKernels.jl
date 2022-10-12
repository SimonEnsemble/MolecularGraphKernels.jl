module MolecularGraphKernels

using Cairo,
    Colors,
    Compose,
    FIGlet,
    GraphPlot,
    Graphs,
    MetaGraphs,
    MolecularGraph,
    RDKitMinimalLib,
    SparseArrays,
    Xtals
using Graphs.Experimental: vf2, IsomorphismProblem
using PrecompileSignatures: @precompile_signatures

import Base: display, size
import Graphs: is_directed, SimpleEdge
import MetaGraphs: weighttype, PropDict, MetaDict, set_props!, props
import SparseArrays: AbstractSparseMatrixCSC, _checkbuffers, getcolptr, rowvals, nonzeros

"""
abstract type for product graphs
"""
abstract type AbstractProductGraph end

"""
concrete product graph types
"""
struct Modular <: AbstractProductGraph end
struct Direct <: AbstractProductGraph end

function __init__()
    if !Sys.iswindows()
        # unpack the RDKit SMARTS
        global maccs_queries =
            [(x[1] == "?" ? nothing : get_qmol(x[1])) for x in rdkit_maccs_smarts_patterns]
        global maccs_counts = [x[2] for x in rdkit_maccs_smarts_patterns]
        global query_notnothing = (!(isnothing)).(maccs_queries)
        global query_notnothing_idx = findall(query_notnothing)
    end
end

include.(
    [
        "ProductGraph.jl"
        "graph_products.jl"
        "graph_kernels.jl"
        "visualization.jl"
        "graph_conversion.jl"
        "check_isom.jl"
        "misc.jl"
        "maccs_smarts.jl"
        "maccs.jl"
    ]
)

export

    # graph_products.jl
    ProductGraph,
    product_graph_adjacency_matrix,
    Modular,
    Direct,

    # random_walks.jl
    random_walk,

    # graph_conversion.jl
    MetaGraph,

    # visualization.jl
    viz_graph,

    # check_isom.jl
    is_isomorphic,

    # misc.jl
    display

@precompile_signatures(MolecularGraphKernels)

end
