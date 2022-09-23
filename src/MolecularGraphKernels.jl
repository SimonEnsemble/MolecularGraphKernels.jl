module MolecularGraphKernels

using FIGlet, Graphs, MetaGraphs, MolecularGraph, SparseArrays

include.([
    "direct_product_graph.jl"
    "random_walk_graph_kernel.jl"
    "graph_conversion.jl"
])


"""
    RWGKSVM.banner()
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])
end


export direct_product_graph, fixed_length_rw_kernel, MetaGraph

end # module MolecularGraphKernels
