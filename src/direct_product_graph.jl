"""
    dpg_adj_mat = direct_product_graph(graph_a, graph_b)

Returns the adjacency matrix of the direct product graph of `graph_a` and `graph_b`.

!!! note
    If only the adjacency matrix of the DPG is needed (and not the DPG itself) `@ref(dpg_adj_mat)` should be used instead for performance.

!!! warning
    Edge and vertex labels should be of consistent, comparable types, preferably `Symbol`, `Int`, `Bool`, etc., for efficiency.
    Comparing, for example, `String` labels will result in significantly poorer performance!
"""
function direct_product_graph(mol_a::MetaGraph, mol_b::MetaGraph; node_label::Symbol=:label, edge_label::Symbol=:label)::MetaGraph
    # map node pairs from input graphs to vertices of direct product graph
    vertex_in_dpg = Dict{Tuple{Int, Int}, NamedTuple{(:l, :n), Tuple{Symbol, Int}}}()
    n_verts = 0
    for b in vertices(mol_b)
        l_b = props(mol_b, b)[node_label]
        for a in vertices(mol_a)
            if props(mol_a, a)[node_label] == l_b
                n_verts += 1
                vertex_in_dpg[a, b] = (l=l_b, n=n_verts) # (vᵢ, vⱼ) ↦ (lₓ, nₓ)
            end
        end
    end

    # instantiate graph to correct size
    axb = MetaGraph(n_verts)

    # label nodes from mapping
    for v in values(vertex_in_dpg)
        set_prop!(axb, v.n, node_label, v.l)
    end

    # make edges
    for e_a in edges(mol_a)
        for e_b in edges(mol_b)
            # only a candidate if edge labels are the same
            if props(mol_a, e_a)[edge_label] != props(mol_b, e_b)[edge_label]
                continue
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b1)
            #   v2 = (a2, b2)
            if haskey(vertex_in_dpg, (src(e_a), src(e_b))) && haskey(vertex_in_dpg, (dst(e_a), dst(e_b)))
                add_edge!(axb, vertex_in_dpg[src(e_a), src(e_b)].n, vertex_in_dpg[dst(e_a), dst(e_b)].n, Dict(edge_label => props(mol_a, e_a)[edge_label]))
            end

            # check candidate edge (v1, v2) in axb:
            #   v1 = (a1, b2)
            #   v2 = (a2, b1)
            if haskey(vertex_in_dpg, (src(e_a), dst(e_b))) && haskey(vertex_in_dpg, (dst(e_a), src(e_b))) # if v1 and v2 are vertices in axb...
                add_edge!(axb, vertex_in_dpg[src(e_a), dst(e_b)].n, vertex_in_dpg[dst(e_a), src(e_b)].n, Dict(edge_label => props(mol_a, e_a)[edge_label]))
            end
        end
    end

    return axb
end

direct_product_graph(mol_a::MetaGraph, mol_b::GraphMol; kwargs...) = direct_product_graph(mol_a, MetaGraph(mol_b); kwargs...)
direct_product_graph(mol_a::GraphMol, mol_b::T; kwargs...) where T <: Union{MetaGraph, GraphMol} = direct_product_graph(MetaGraph(mol_a), mol_b; kwargs...)

"""
    A = direct_product_graph(graph_a, graph_b)

Returns the adjacency matrix of the DPG.
"""
function dpg_adj_mat(mol_a::MetaGraph, mol_b::MetaGraph; node_label::Symbol=:label, edge_label::Symbol=:label)::SparseMatrixCSC{Bool, Int}
    # map (vⱼ ∈ a, vₖ ∈ b) ↦ vᵢ ∈ a*b
    vertex_in_dpg = spzeros(Int, nv(mol_a), nv(mol_b))
    for b in vertices(mol_b)
        l_b = props(mol_b, b)[node_label]
        for a in vertices(mol_a)
            vertex_in_dpg[a, b] = props(mol_a, a)[node_label] == l_b
        end
    end
    n_verts = sum(vertex_in_dpg)
    vertex_in_dpg[vertex_in_dpg .== 1] .= 1:n_verts

    # build graph adjacency matrix
    adj_mat = spzeros(Bool, n_verts, n_verts)
    for e_a in edges(mol_a)
        for e_b in edges(mol_b)
            # only a candidate if edge labels are the same
            if props(mol_a, e_a)[edge_label] != props(mol_b, e_b)[edge_label]
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

    return adj_mat .|| adj_mat'
end

dpg_adj_mat(mol_a::MetaGraph, mol_b::GraphMol; kwargs...) = dpg_adj_mat(mol_a, MetaGraph(mol_b); kwargs...)
dpg_adj_mat(mol_a::GraphMol, mol_b::T; kwargs...) where T <: Union{MetaGraph, GraphMol} = dpg_adj_mat(MetaGraph(mol_a), mol_b; kwargs...)
