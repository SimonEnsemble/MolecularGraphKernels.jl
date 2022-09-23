"""
    g = MetaGraph(mol)

Convert a `GraphMol` object into the corresponding `MetaGraph`
"""
function MetaGraph(mol::G)::MetaGraph where G <: GraphMol
    g = MetaGraph(length(mol.nodeattrs))
    # atoms
    for v in vertices(g)
        set_prop!(g, v, :symbol, mol.nodeattrs[v].symbol)
    end
    # bonds
    for (e, att) in enumerate(mol.edgeattrs)
        add_edge!(g, mol.edges[e]..., Dict(:order => att.order))
    end
    return g
end
