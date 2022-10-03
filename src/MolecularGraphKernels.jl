module MolecularGraphKernels

using Cairo, Colors, Compose, FIGlet, GraphPlot, Graphs, MetaGraphs, MolecularGraph, SparseArrays, Xtals
using Graphs.Experimental: vf2, IsomorphismProblem

import Graphs: is_directed, SimpleEdge
import MetaGraphs: weighttype, PropDict, MetaDict

"""
abstract type for product graphs
"""
abstract type AbstractProductGraph end

"""
concrete product graph types
"""
struct Modular <: AbstractProductGraph end
struct Direct <: AbstractProductGraph end

include.([
    "ProductGraph.jl"
    "ProductGraphMatrix.jl"
    "graph_products.jl"
    "graph_kernels.jl"
    "visualization.jl"
    "graph_conversion.jl"
    "check_isom.jl"
    "misc.jl"
])

export

    # graph_products.jl
    ProductGraph, ProductGraphMatrix, Modular, Direct,

    # random_walk_graph_kernels.jl
    random_walk_graph_kernel, 
    
    # graph_conversion.jl
    MetaGraph, 
    
    # visualization.jl
    viz_graph,

    # check_isom.jl
    is_isomorphic,

    # misc.jl
    display, isequal

end
