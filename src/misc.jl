import Base: display

"""
prints the node and edge lists of a graph
"""
function display(g::AbstractMetaGraph)
    function get_props(i)
        prop_vec = ["$k:$v" for (k, v) in props(g, i)]
        return reduce(*, [" "] .* prop_vec)
    end

    println("---NODES---")
    for i in 1:nv(g)
        println("[$i]", get_props(i))
    end

    println("---EDGES---")
    for ed in edges(g)
        println("($(ed.src), $(ed.dst))", get_props(ed))
    end
end

"""
Prints the stylized ASCII console banner for the package.
"""
banner() = FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])

"""
extracts the isomorphic subgraphs of the modular product graph via clique detection
"""
function isomorphic_subgraphs(
    mpg::ProductGraph{Modular};
    min_size::Int=3
)::Vector{ProductGraph{Modular}}
    max_cliques = maximal_cliques(mpg.graph)
    cliques = filter(c -> length(c) ≥ min_size, max_cliques)
    tups =
        induced_subgraph.([MetaGraph(mpg)], cliques[sortperm(length.(cliques); rev=true)])
    return [ProductGraph{Modular}(tup[1]) for tup in tups]
end

export display, isomorphic_subgraphs
