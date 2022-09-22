"""
    dpg_adj_mat = direct_product_graph(graph_a, graph_b)

Returns the direct product graph of `graph_a` and `graph_b`.
"""
function direct_product_graph(mol_a::MetaGraph, mol_b::MetaGraph)::SparseMatrixCSC{Int, Int}
    # input validation (no trivial graphs)
    @assert nv(mol_a) > 0 && nv(mol_b) > 0 && ne(mol_a) > 0 && ne(mol_b) > 0 "Graphs A and B must each contain at least 1 node and 1 edge."
    # input validation (check node symbol and edge order types for type-stable and efficient comparisons)
    function _type_assertion(g::MetaGraph)::Bool
        typeof(get_prop(g, 1, :symbol)) <: Union{Symbol, Int, Vector{T} where T <: Union{Symbol, Int}} &&
        typeof(get_prop(g, first(edges(g)), :order)) <: Union{Symbol, Int}
    end
    @assert all(_type_assertion.([mol_a, mol_b])) "Graphs A and B must have appropriate node and edge attribute types (see docs)."
    @assert typeof(mol_a.eprops) == typeof(mol_b.eprops) && typeof(mol_a.vprops) == typeof(mol_b.vprops) "Graphs A and B must have type-matched node and edge attribute types."

	# determine which vertices exist in the DPG
    vertex_in_dpg = spzeros(Int, nv(mol_a), nv(mol_b))
	for b in eachindex(vertices(mol_b))
		for a in eachindex(vertices(mol_a))
			vertex_in_dpg[a, b] = props(mol_a, a)[:symbol] == props(mol_b, b)[:symbol]
		end
	end
	n_verts = sum(vertex_in_dpg)

	# map (vⱼ ∈ a, vₖ ∈ b) ↦ vᵢ ∈ g
	vertex_in_dpg[vertex_in_dpg .== 1] .= 1:n_verts

    # build graph adjacency matrix
	adj_mat = spzeros(Bool, n_verts, n_verts)
    for e_a in edges(mol_a)
        for e_b in edges(mol_b)
            # only a candidate if edge labels are the same
            if props(mol_a, e_a)[:order] != props(mol_b, e_b)[:order]
                continue
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b1)
            #   v2 = (a2, b2)
            v1 = vertex_in_dpg[src(e_a), src(e_b)]
            v2 = vertex_in_dpg[dst(e_a), dst(e_b)]
            if (v1 != 0) && (v2 != 0) # if v1 and v2 are vertices in axb...
                adj_mat[v1, v2] = 1
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b2)
            #   v2 = (a2, b1)
            v1 = vertex_in_dpg[src(e_a), dst(e_b)]
            v2 = vertex_in_dpg[dst(e_a), src(e_b)]
            if (v1 != 0) && (v2 != 0) # if v1 and v2 are vertices in axb...
                adj_mat[v1, v2] = 1
            end
        end
    end

    return adj_mat + adj_mat'
end
