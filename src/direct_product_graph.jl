"""
    g = dpg(mat, map, mol)

Returns the direct product graph based on its adjacency matrix, the map of A/B nodes to DPG nodes, and graph A
"""
function dpg(adj_mat::AbstractMatrix, vertex_map::AbstractMatrix, A::MetaGraph)
    g = MetaGraph(SimpleGraph(adj_mat))
    # each vertex in the DPG has the label of the vertex in A corresponding to some row in vertex_map
    for v in vertices(g)
        a = findfirst(isequal(v), vertex_map)[1]
        set_prop!(g, v, :label, props(A, a)[:label])
    end
    # each edge n in the DPG corresponds to (nᵢ, nⱼ) in the vertex map where nₓ ↦ vᵪ ∈ V(A); ∴ L(n) = L(vᵪ₁, vᵪ₂)
    for e in edges(g)
        vᵢ = findfirst(isequal(src(e)), vertex_map)[1]
        vⱼ = findfirst(isequal(dst(e)), vertex_map)[1]
        l = get_prop(A, vᵢ, vⱼ, :label)
        set_prop!(g, e, :label, l)
    end
    return g
end


"""
    dpg_adj_mat = direct_product_graph(graph_a, graph_b)

Returns the adjacency matrix of the direct product graph of `graph_a` and `graph_b`.

    dpg_adj_mat = direct_product_graph(graph_a, graph_b; return_graph=true)

Returns the adjacency matrix and the 
"""
function direct_product_graph(mol_a::S, mol_b::T; return_graph::Bool=false)::Union{SparseMatrixCSC{Bool, Int}, MetaGraph} where {S,T <: Union{GraphMol, MetaGraph}}
    typeof(mol_a) <: GraphMol ? mol_a = MetaGraph(mol_a) : nothing
    typeof(mol_b) <: GraphMol ? mol_b = MetaGraph(mol_b) : nothing
    # input validation (no trivial graphs)
    @assert nv(mol_a) > 0 && nv(mol_b) > 0 && ne(mol_a) > 0 && ne(mol_b) > 0 "Graphs A and B must each contain at least 1 node and 1 edge."
    # input validation (check node symbol and edge order types for type-stable and efficient comparisons)
    function _type_assertion(g::MetaGraph)::Bool
        typeof(get_prop(g, 1, :label)) <: Union{Symbol, Int, Vector{Union{Symbol, Int}}} &&
        typeof(get_prop(g, first(edges(g)), :label)) <: Union{Symbol, Int}
    end
    if !all(_type_assertion.([mol_a, mol_b]))
        @warn "Graphs A and B should have Int or Symbol node and edge attribute types for efficiency."
    end

    # determine which vertices exist in the DPG
    vertex_in_dpg = spzeros(Int, nv(mol_a), nv(mol_b))
    for b in vertices(mol_b)
        for a in vertices(mol_a)
            vertex_in_dpg[a, b] = props(mol_a, a)[:label] == props(mol_b, b)[:label]
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
            if props(mol_a, e_a)[:label] != props(mol_b, e_b)[:label]
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

    adj_mat = adj_mat .|| adj_mat'

    if return_graph
        return dpg(adj_mat, vertex_in_dpg, mol_a)
    else
        return adj_mat
    end
end
