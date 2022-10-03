import Base.display
"""
prints the node and edge lists of a graph
"""
function display(g::MetaGraph)
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

display(g::ProductGraph{T}) where T <: AbstractProductGraph = display(g.graph)

"""
Prints the stylized ASCII console banner for the package.
"""
function banner()
    FIGlet.render("MolecularGraphKernels", FIGlet.availablefonts()[449])
end

import Base.isequal
"""
determines if two product graphs are not just isomorphic, but identical
"""
function isequal(a::ProductGraph, b::ProductGraph)::Bool
    if nv(a.graph) ≠ nv(b.graph)
        return false
    end
    for v in vertices(a.graph)
        if props(a.graph, v) ≠ props(b.graph, v)
            return false
        end
    end
    if ne(a.graph) ≠ ne(b.graph)
        return false
    end
    for e in edges(a.graph)
        if props(a.graph, e) ≠ props(b.graph, e)
            return false
        end
    end
    return true
end
