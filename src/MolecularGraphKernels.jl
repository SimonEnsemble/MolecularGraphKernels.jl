module MolecularGraphKernels

using Cairo,
    Colors,
    Compose,
    FIGlet,
    GraphPlot,
    Graphs,
    MetaGraphs,
    MolecularGraph,
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

include.(
    [
        "ProductGraph.jl"
        "graph_products.jl"
        "graph_kernels.jl"
        "visualization.jl"
        "graph_conversion.jl"
        "check_isom.jl"
        "misc.jl"
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
