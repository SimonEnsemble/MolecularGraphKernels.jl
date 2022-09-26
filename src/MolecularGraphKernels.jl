module MolecularGraphKernels

using Cairo, Colors, Compose, FIGlet, GraphPlot, Graphs, MetaGraphs, MolecularGraph, SparseArrays

include.([
    "direct_product_graph.jl"
    "random_walk_graph_kernel.jl"
    "graph_conversion.jl"
    "visualization.jl"
])


"""
    RWGKSVM.banner()
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])
end


export direct_product_graph, dpg_adj_mat, fixed_length_rw_kernel, MetaGraph, viz_graph

end # module MolecularGraphKernels
