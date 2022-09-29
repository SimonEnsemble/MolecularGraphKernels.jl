#=
#   NOTATION IN CODE
#   v₁ is a node in g₁
#   v₂ is a node in g₂
#       w = (v₁, v₂) is a vertex in the product graph
=#

"""
map (v₁ ∈ g₁, v₂ ∈ g₂) ↦ w ∈ g₁ x g₂
"""
function _build_v₁v₂_pair_to_w_map(g₁::MetaGraph, g₂::MetaGraph)
    v₁v₂_pair_to_w = spzeros(Int, nv(g₁), nv(g₂))
    n_g₁xg₂_verts = 0 # nb of vertices in product graph (TBD)
    for v₂ in vertices(g₂) # switch order b/c COLUMN sparse matrix
        l_v₂ = get_prop(g₂, v₂, :label) # label on g₂
        for v₁ in vertices(g₁)
            # add vertex pair to g₁xg₂ iff they share label.
            l_v₁ = get_prop(g₁, v₁, :label)
            if l_v₁ == l_v₂
                n_g₁xg₂_verts += 1
                v₁v₂_pair_to_w[v₁, v₂] = n_g₁xg₂_verts
            end
        end
    end
    return v₁v₂_pair_to_w
end

"""
map w ∈ g₁ x g₂ ↦ (v₁ ∈ g₁, v₂ ∈ g₂)
"""
function _build_w_to_v₁v₂_pair_map(v₁v₂_pair_to_w::SparseMatrixCSC)
    v₁v₂_pairs = findall(! iszero, v₁v₂_pair_to_w) # vector of CartesianIndices
    return [(w[1], w[2]) for w in v₁v₂_pairs] # vector of tuples
end

function product_graph(g₁::MetaGraph, g₂::MetaGraph, type::Symbol)::MetaGraph
    # get (v₁ ∈ g₁, v₂ ∈ g₂) <--> w ∈ g₁xg₂ mappings
    v₁v₂_pair_to_w = _build_v₁v₂_pair_to_w_map(g₁, g₂)
    w_to_v₁v₂_pair = _build_w_to_v₁v₂_pair_map(v₁v₂_pair_to_w)

    # pre-allocate product graph
    n_g₁xg₂_vertices = length(w_to_v₁v₂_pair)
    g₁xg₂ = MetaGraph(n_g₁xg₂_vertices)

    @assert type in [:factor, :direct]

    # loop over pairs of vertices in the product graph
    for wᵢ = 1:n_g₁xg₂_vertices
        # product graph vertex wᵢ = (u₁, u₂)
        #    with u₁ ∈ g₁, u₂ ∈ g₂
        u₁, u₂ = w_to_v₁v₂_pair[wᵢ]
        for wⱼ = (wᵢ + 1):n_g₁xg₂_vertices
            # product graph vertex wⱼ = (v₁, v₂)
            #    with v₁ ∈ g₁, v₂ ∈ g₂
            v₁, v₂ = w_to_v₁v₂_pair[wⱼ]

            # does edge (u₁, v₁) exist in g₁?
            e₁ = has_edge(g₁, u₁, v₁)
            # does edge (u₂, v₂) exist in g₂?
            e₂ = has_edge(g₂, u₂, v₂)

            # is there a common adjacency? (c-edge) "c" = connected.
            if e₁ && e₂
                # do they share an edge label?
                l₁ = get_prop(g₁, u₁, v₁, :label)
                l₂ = get_prop(g₂, u₂, v₂, :label)
                if l₁ == l₂
                    add_edge!(g₁xg₂, wᵢ, wⱼ, Dict{Symbol, Int}(:label => l₁))
                end
            # is there a common non-adjacency? (d-edge) "d" = disconnected.
            #   (only relevant to factor graph)
            #   see "Subgraph Matching Kernels for Attributed Graphs"
            elseif (type == :factor) && !e₁ && !e₂ && (u₁ != v₁) && (u₂ != v₂)
                add_edge!(g₁xg₂, wᵢ, wⱼ, Dict{Symbol, Int}(:label => 0))
            end
        end
    end

    # loop over vertices in g₁xg₂
    for w = 1:n_g₁xg₂_vertices
        v₁, v₂ = w_to_v₁v₂_pair[w]
        set_prop!(g₁xg₂, w, :v₁v₂_pair, (v₁, v₂))
        set_prop!(g₁xg₂, w, :label, get_prop(g₁, v₁, :label))
    end
    return g₁xg₂
end
