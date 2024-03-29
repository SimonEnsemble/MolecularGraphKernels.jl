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

const BANNER = String(read(joinpath(dirname(pathof(MolecularGraphKernels)), "banner.txt")))

"""
Prints the stylized ASCII console banner for the package.
"""
banner() = println(BANNER)

export display
