"""
    g = MetaGraph(mol)

Convert a `GraphMol` object into the corresponding `MetaGraph`

!!! note
    
    Hydrogen atoms are generally treated implicitly.
"""
function MetaGraph(mol::GraphMol)::MetaGraph
    g = MetaGraph(length(mol.nodeattrs))
    # atoms
    for v in vertices(g)
        set_prop!(g, v, :label, atomnumber(mol.nodeattrs[v].symbol))
    end
    # bonds
    for (e, att) in enumerate(mol.edgeattrs)
        add_edge!(g, mol.edges[e]..., Dict(:label => att.order))
    end
    # correct aromaticity
    for e_idx in findall(isaromaticbond(mol))
        set_prop!(g, mol.edges[e_idx]..., :label, -1)
    end
    # set atomic coordinates
    coords, _ = coords2d(mol)
    for v in vertices(g)
        set_prop!(g, v, :coords, coords[v, :])
    end
    return g
end

"""
converts a modular product graph into the corresponding direct product graph
"""
function ProductGraph{Direct}(fpg::ProductGraph{Modular})::ProductGraph{Direct}
    dpg = ProductGraph{Direct}(fpg.graph)
    for e in edges(fpg)
        if get_prop(fpg, e, :label) == D_EDGE
            rem_edge!(dpg, e)
        end
    end
    return dpg
end

"""
convert a product graph into the corresponding simple graph
"""
SimpleGraph(g::T) where {T <: ProductGraph} = g.graph

"""
convert a product grpah into the corresponding metagraph
"""
function MetaGraph(g::T) where {T <: ProductGraph}
    return MetaGraphs.MetaGraph(
        g.graph,
        g.vprops,
        g.eprops,
        g.gprops,
        g.weightfield,
        g.defaultweight,
        g.metaindex,
        g.indices
    )
end

export MetaGraph
