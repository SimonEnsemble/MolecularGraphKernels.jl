"""
    cpg = csi_product_graph(graph_a, graph_b)

Returns the common subgraph isomorphism product graph of `graph_a` and `graph_b`
"""
function csi_product_graph(mol_a::MetaGraph, mol_b::MetaGraph; node_label::Symbol=:label, edge_label::Symbol=:label)::MetaGraph
    # map node pairs from input graphs to vertices of direct product graph
    vertex_in_dpg = Dict{Tuple{Int, Int}, NamedTuple{(:l, :n), Tuple{Int, Int}}}()
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
    for (a, b) in [[a, b] for a in vertex_in_dpg for b in vertex_in_dpg if b[2].n > a[2].n && a[1][1] ≠ b[1][1] && a[1][2] ≠ b[1][2]]
        # `a` and `b` are key => value pairs where key::Tuple{Int,Int} and value::NamedTuple
        # if the vertex pair in the product graph maps to correctly connected edges in the source graphs, make c-type edge
        e1 = has_edge(mol_a, a[1][1], b[1][1])
        e2 = has_edge(mol_b, a[1][2], b[1][2])
        if e1 && e2 && props(mol_a, a[1][1], b[1][1])[edge_label] == props(mol_b, a[1][2], b[1][2])[edge_label]
            add_edge!(axb, a[2].n, b[2].n, Dict(edge_label => props(mol_a, a[1][1], b[1][1])[edge_label], :d_type => false))
        # if not a c-type edge, and neither source graph contains the corresponding edge, make d-type edge
        elseif !e1 && !e2
            add_edge!(axb, a[2].n, b[2].n, Dict(edge_label => 1, :d_type => true))
        end
    end

    return axb
end

csi_product_graph(mol_a::MetaGraph, mol_b::GraphMol; kwargs...) = csi_product_graph(mol_a, MetaGraph(mol_b); kwargs...)
csi_product_graph(mol_a::GraphMol, mol_b::T; kwargs...) where T <: Union{MetaGraph, GraphMol} = csi_product_graph(MetaGraph(mol_a), mol_b; kwargs...)


"""
    cpg = csi_adj_mat(graph_a, graph_b)

Returns the adjacency matrix of the CSI product graph of `graph_a` and `graph_b`
"""
function csi_adj_mat(mol_a::MetaGraph, mol_b::MetaGraph; node_label::Symbol=:label, edge_label::Symbol=:label)::SparseMatrixCSC
    # map (vⱼ ∈ a, vₖ ∈ b) ↦ vᵢ ∈ a*b
    vertex_in_dpg = spzeros(Bool, nv(mol_a), nv(mol_b))
    for b in vertices(mol_b)
        l_b = props(mol_b, b)[node_label]
        for a in vertices(mol_a)
            vertex_in_dpg[a, b] = props(mol_a, a)[node_label] == l_b
        end
    end
    n_verts = sum(vertex_in_dpg)

    # product graph nodes ↦ source node pairs
    vertex_map = 1:n_verts .=> findall(vertex_in_dpg)

    # build graph adjacency matrix
    adj_mat = spzeros(Bool, n_verts, n_verts)
    for (a, b) in [[a, b] for a in vertex_map for b in vertex_map if b[1] > a[1] && a[2][1] ≠ b[2][1] && a[2][2] ≠ b[2][2]]
        # if the vertex pair in the product graph maps to correctly connected edges in the source graphs, make edge
        # else, if neither source graph contains the corresponding edge, make edge
        e1 = has_edge(mol_a, a[2][1], b[2][1])
        e2 = has_edge(mol_b, a[2][2], b[2][2])
        adj_mat[a[1], b[1]] = (e1 && e2 && props(mol_a, a[2][1], b[2][1])[edge_label] == props(mol_b, a[2][2], b[2][2])[edge_label]) || !e1 && !e2
    end

    return adj_mat .|| adj_mat'
end

csi_adj_mat(mol_a::MetaGraph, mol_b::GraphMol; kwargs...) = csi_adj_mat(mol_a, MetaGraph(mol_b); kwargs...)
csi_adj_mat(mol_a::GraphMol, mol_b::T; kwargs...) where T <: Union{MetaGraph, GraphMol} = csi_adj_mat(MetaGraph(mol_a), mol_b; kwargs...)


"""
    dpg = direct_product_graph(graph_a, graph_b)

Returns the direct product graph of `graph_a` and `graph_b`.

!!! note
    If only the adjacency matrix of the DPG is needed (and not the DPG itself) `@ref(dpg_adj_mat)` should be used instead for performance.

!!! warning
    Edge and vertex labels should be of consistent, comparable types, preferably `Symbol`, `Int`, `Bool`, etc., for efficiency.
    Comparing, for example, `String` labels will result in significantly poorer performance!
"""
function direct_product_graph(mol_a::MetaGraph, mol_b::MetaGraph; node_label::Symbol=:label, edge_label::Symbol=:label)::MetaGraph
    # map node pairs from input graphs to vertices of direct product graph
    vertex_in_dpg = Dict{Tuple{Int, Int}, NamedTuple{(:l, :n), Tuple{Int, Int}}}()
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
    A = dpg_adj_mat(graph_a, graph_b)

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
