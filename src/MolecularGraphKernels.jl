module MolecularGraphKernels

using FIGlet, Graphs, MetaGraphs, SparseArrays

include.([
    "direct_product_graph.jl"
    "random_walk_graph_kernel.jl"
])


"""
    RWGKSVM.banner()
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])
end


export direct_product_graph, fixed_length_rw_kernel

end # module MolecularGraphKernels
