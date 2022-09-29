# map (u ∈ g₁, v ∈ g₂) ↦ w ∈ g₁ x g₂
function _build_g₁g₂_vertex_pair_to_g₁xg₂_vertex_map(g₁::MetaGraph, g₂::MetaGraph; node_label::Symbol=:label)
    g₁g₂_vertex_pair_to_g₁xg₂_vertex_id = spzeros(Int, nv(g₁), nv(g₂))
    n_g₁xg₂_verts = 0 # nb of vertices in product graph (TBD)
    for v in vertices(g₂) # switch order b/c COLUMN sparse matrix
        l_v = props(g₂, v)[node_label] # label on g₁
        for u in vertices(g₁)
            # add vertex pair to g₁xg₂ iff they share label.
            l_u = props(g₁, u)[node_label]
            if l_v == l_u
                n_g₁xg₂_verts += 1
                g₁g₂_vertex_pair_to_g₁xg₂_vertex_id[u, v] = n_g₁xg₂_verts
            end
        end
    end
    return g₁g₂_vertex_pair_to_g₁xg₂_vertex_id
end

function _build_g₁xg₂_vertex_to_g₁g₂_vertex_pair_map(g₁g₂_vertex_pair_to_g₁xg₂_vertex_id::SparseMatrixCSC)
    g₁g₂_vertex_pairs = findall(! iszero, g₁g₂_vertex_pair_to_g₁xg₂_vertex_id) # vector of Cartesian indices
    return [(uv[1], uv[2]) for uv in g₁g₂_vertex_pairs] # vector of tuples
end

function product_graph(g₁::MetaGraph, g₂::MetaGraph, type::String; node_label::Symbol=:label, edge_label::Symbol=:label, add_g₁xg₂_labels::Bool=true)::MetaGraph
    g₁g₂_vertex_pair_to_g₁xg₂_vertex_id = _build_g₁g₂_vertex_pair_to_g₁xg₂_vertex_map(g₁, g₂, node_label=node_label)
    g₁xg₂_vertex_to_g₁g₂_vertex_pair = _build_g₁xg₂_vertex_to_g₁g₂_vertex_pair_map(g₁g₂_vertex_pair_to_g₁xg₂_vertex_id)
    
    # pre-allocate product graph
    n_g₁xg₂_vertices = sum(g₁g₂_vertex_pair_to_g₁xg₂_vertex_id .!= 0)
    g₁xg₂ = MetaGraph(n_g₁xg₂_vertices)

    @assert type in ["factor", "direct"]
    
    # loop over pairs of vertices in the product graph
    for u = 1:n_g₁xg₂_vertices
        u_g₁, u_g₂ = g₁xg₂_vertex_to_g₁g₂_vertex_pair[u]
        for v = (u + 1):n_g₁xg₂_vertices
            # this pair of vertieces in the product graph (u, v) represents
            #  (u_g₁, u_g₂) and (v_g₁, v_g₂).
            v_g₁, v_g₂ = g₁xg₂_vertex_to_g₁g₂_vertex_pair[v]
            
            # does edge (u_g₁, v_g₁) exist in g₁?
            e₁ = has_edge(g₁, u_g₁, v_g₁)
            # does edge (u_g₂, v_g₂) exist in g₂?
            e₂ = has_edge(g₂, u_g₂, v_g₂) 
            
            # is there a common adjacency? (c-edge) "c" = connected. 
            if e₁ && e₂ 
                # do they share an edge label?
                l₁ = props(g₁, u_g₁, v_g₁)[edge_label]
                l₂ = props(g₂, u_g₂, v_g₂)[edge_label]
                if l₁ == l₂
                    add_edge!(g₁xg₂, u, v, Dict(edge_label => l₁, :d_type => false))
                end
            # is there a common non-adjacency? (d-edge) "d" = disconnected. only relevant to factor graph.
            elseif (type == "factor") && !e₁ && !e₂
                add_edge!(g₁xg₂, u, v, Dict(edge_label => 0, :d_type => true))
            end
        end
    end

    if add_g₁xg₂_labels
        for u = 1:n_g₁xg₂_vertices
            set_prop!(g₁xg₂, u, :g₁g₂_vertex_pair, (g₁xg₂_vertex_to_g₁g₂_vertex_pair[1], g₁xg₂_vertex_to_g₁g₂_vertex_pair[2]))
            set_prop!(g₁xg₂, u, node_label, get_prop(g₁, g₁xg₂_vertex_to_g₁g₂_vertex_pair[1], node_label))
        end
    end
    return g₁xg₂
end
