module MolecularGraphKernels

using Cairo, Colors, Compose, FIGlet, GraphPlot, Graphs, MetaGraphs, MolecularGraph, SparseArrays, Xtals

include.([
    "graph_products.jl"
    "graph_kernels.jl"
    "visualization.jl"
    "graph_conversion.jl"
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
    direct_product_graph, dpg_adj_mat, csi_product_graph, csi_adj_mat,

    # graph_kernels.jl
    random_walk_kernel, 
    
    # graph_conversion.jl
    MetaGraph, 
    
    # visualization.jl
    viz_graph

end
