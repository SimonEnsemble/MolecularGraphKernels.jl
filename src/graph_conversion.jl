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
    return g
end
