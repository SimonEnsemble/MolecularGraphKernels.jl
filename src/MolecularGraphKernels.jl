module MolecularGraphKernels

using Cairo, Colors, Compose, FIGlet, GraphPlot, Graphs, MetaGraphs, MolecularGraph, SparseArrays, Xtals
using Graphs.Experimental: vf2, IsomorphismProblem

include.([
    "graph_products.jl"
    "graph_kernels.jl"
    "visualization.jl"
    "graph_conversion.jl"
    "check_isom.jl"
])

"""
    RWGKSVM.banner()
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])
end

export

    # graph_products.jl
    ProductGraph, ProductGraphMatrix, Factor, Direct,

    # graph_kernels.jl
    graph_kernel, 
    
    # graph_conversion.jl
    MetaGraph, 
    
    # visualization.jl
    viz_graph,

    # check_isom.jl
    is_isomorphic

end
