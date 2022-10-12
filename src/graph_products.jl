#=
#   NOTATION IN CODE
#   v₁ is a node in g₁
#   v₂ is a node in g₂
#       w = (v₁, v₂) is a vertex in the product graph
=#
"""
map (v₁ ∈ g₁, v₂ ∈ g₂) ↦ w ∈ g₁ x g₂
"""
function build_v₁v₂_pair_to_w_map(g₁::MetaGraph, g₂::MetaGraph)::SparseMatrixCSC{Int, Int}
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
function build_w_to_v₁v₂_pair_map(
    v₁v₂_pair_to_w::SparseMatrixCSC{Int, Int}
)::Vector{Tuple{Int, Int}}
    v₁v₂_pairs = findall(!iszero, v₁v₂_pair_to_w) # vector of CartesianIndices
    return [(w[1], w[2]) for w in v₁v₂_pairs] # vector of tuples
end

"""
applies the direct product graph adjacency matrix marking rule: common adjacency
"""
function record_adjacency!(
    A::AbstractMatrix,
    ::Type{Direct},
    g₁::MetaGraph,
    g₂::MetaGraph,
    e₁::Bool,
    e₂::Bool,
    u₁::Int,
    u₂::Int,
    v₁::Int,
    v₂::Int,
    wᵢ::Int,
    wⱼ::Int
)::Bool
    if e₁ && e₂
        if get_prop(g₁, u₁, v₁, :label) == get_prop(g₂, u₂, v₂, :label)
            A[wᵢ, wⱼ] = 1
            A[wⱼ, wᵢ] = 1
        end
        return true
    end
    return false
end

"""
applies the modular product graph adjacency matrix marking rules: common adjacency (a.k.a. the direct product graph rule) and common non-adjacency
"""
function record_adjacency!(
    A::AbstractMatrix,
    ::Type{Modular},
    g₁::MetaGraph,
    g₂::MetaGraph,
    e₁::Bool,
    e₂::Bool,
    u₁::Int,
    u₂::Int,
    v₁::Int,
    v₂::Int,
    wᵢ::Int,
    wⱼ::Int
)
    # is there a common adjacency? (c-edge) "c" = connected. 
    if record_adjacency!(A, Direct, g₁, g₂, e₁, e₂, u₁, u₂, v₁, v₂, wᵢ, wⱼ)
        # is there a common non-adjacency? (d-edge) "d" = disconnected.
        #   (only relevant to modular graph) 
        #   see "Subgraph Matching Kernels for Attributed Graphs" 
    elseif !e₁ && !e₂ && (u₁ != v₁) && (u₂ != v₂)
        A[wᵢ, wⱼ] = 1
        A[wⱼ, wᵢ] = 1
    end
end

"""
compute the product graph adjacency matrix and mappings (in each direction) for relating the source graphs and product graph
"""
function product_graph_matrix_and_maps(
    type::Type{T},
    g₁::MetaGraph,
    g₂::MetaGraph
)::Tuple{AbstractMatrix, SparseMatrixCSC, Vector} where {T <: AbstractProductGraph}
    # get (v₁ ∈ g₁, v₂ ∈ g₂) <--> w ∈ g₁xg₂ mappings
    v₁v₂_pair_to_w = build_v₁v₂_pair_to_w_map(g₁, g₂)
    w_to_v₁v₂_pair = build_w_to_v₁v₂_pair_map(v₁v₂_pair_to_w)

    # pre-allocate adj matrix
    n_g₁xg₂ = length(w_to_v₁v₂_pair)
    A = spzeros(Bool, n_g₁xg₂, n_g₁xg₂)

    # loop over pairs of vertices in the product graph
    for wᵢ in 1:n_g₁xg₂
        # product graph vertex wᵢ = (u₁, u₂)
        #    with u₁ ∈ g₁, u₂ ∈ g₂
        u₁, u₂ = w_to_v₁v₂_pair[wᵢ]
        for wⱼ in (wᵢ + 1):n_g₁xg₂
            # product graph vertex wⱼ = (v₁, v₂)
            #    with v₁ ∈ g₁, v₂ ∈ g₂
            v₁, v₂ = w_to_v₁v₂_pair[wⱼ]
            # does edge (u₁, v₁) exist in g₁?
            e₁ = has_edge(g₁, u₁, v₁)
            # does edge (u₂, v₂) exist in g₂?
            e₂ = has_edge(g₂, u₂, v₂)
            # apply the adjacency matrix marking rules for the product graph type
            record_adjacency!(A, type, g₁, g₂, e₁, e₂, u₁, u₂, v₁, v₂, wᵢ, wⱼ)
        end
    end
    return A, v₁v₂_pair_to_w, w_to_v₁v₂_pair
end

"""
return the appropriate label for an edge in a product graph given its type
"""
function product_graph_edge_label(::Type{Direct}, g₁::MetaGraph, u₁::Int, v₁::Int)
    return get_prop(g₁, u₁, v₁, :label)
end
function product_graph_edge_label(::Type{Modular}, g₁::MetaGraph, u₁::Int, v₁::Int)
    return has_edge(g₁, u₁, v₁) ? get_prop(g₁, u₁, v₁, :label) : 0
end

"""
compute the product graph of given type between g₁ and g₂
"""
function product_graph(
    type::Type{T},
    g₁::MetaGraph,
    g₂::MetaGraph
)::ProductGraph{T} where {T <: AbstractProductGraph}
    A, _, w_to_v₁v₂_pair = product_graph_matrix_and_maps(type, g₁, g₂)
    n_g₁xg₂ = length(w_to_v₁v₂_pair)

    g₁xg₂ = ProductGraph{T}(MetaGraph(SimpleGraph(A)))

    # label vertices
    for w in 1:n_g₁xg₂
        v₁, v₂ = w_to_v₁v₂_pair[w]
        set_prop!(g₁xg₂, w, :v₁v₂_pair, (v₁, v₂))
        set_prop!(g₁xg₂, w, :label, get_prop(g₁, v₁, :label))
    end
    # label edges
    for wᵢ in 1:n_g₁xg₂
        u₁, _ = w_to_v₁v₂_pair[wᵢ]
        for wⱼ in (wᵢ + 1):n_g₁xg₂
            v₁, _ = w_to_v₁v₂_pair[wⱼ]
            if A[wᵢ, wⱼ]
                set_prop!(g₁xg₂, wᵢ, wⱼ, :label, product_graph_edge_label(type, g₁, u₁, v₁))
            end
        end
    end
    return g₁xg₂
end
function product_graph(
    type::Type{T},
    g₁::MetaGraph,
    g₂::GraphMol
) where {T <: AbstractProductGraph}
    return product_graph(type, g₁, MetaGraph(g₂))
end
function product_graph(
    type::Type{T},
    g₁::GraphMol,
    g₂::Union{GraphMol, MetaGraph}
) where {T <: AbstractProductGraph}
    return product_graph(type, MetaGraph(g₁), g₂)
end

"""
compute the adjacency matrix of the product graph of given type between g₁ and g₂, but do not explicitly construct the graph
"""
function product_graph_adjacency_matrix(
    type::Type{T},
    g₁::MetaGraph,
    g₂::MetaGraph
) where {T <: AbstractProductGraph}
    return product_graph_matrix_and_maps(type, g₁, g₂)[1]
end
function product_graph_adjacency_matrix(
    type::Type{T},
    g₁::MetaGraph,
    g₂::GraphMol
) where {T <: AbstractProductGraph}
    return product_graph_adjacency_matrix(type, g₁, MetaGraph(g₂))
end
function product_graph_adjacency_matrix(
    type::Type{T},
    g₁::GraphMol,
    g₂::Union{GraphMol, MetaGraph}
) where {T <: AbstractProductGraph}
    return product_graph_adjacency_matrix(type, MetaGraph(g₁), g₂)
end
