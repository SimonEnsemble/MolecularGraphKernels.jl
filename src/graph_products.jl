#=
#   NOTATION IN CODE
#   v₁ is a node in g₁
#   v₂ is a node in g₂
#       w = (v₁, v₂) is a vertex in the product graph
=#

abstract type AbstractProductGraph end

struct ProductGraph{T <: AbstractProductGraph}
    graph::MetaGraph
end

struct Factor <: AbstractProductGraph end

struct Direct <: AbstractProductGraph end

ProductGraph{T}() where T <: AbstractProductGraph = ProductGraph{T}(MetaGraph())

ProductGraph{T}(n::Union{Int, AbstractMatrix}) where T <: AbstractProductGraph = ProductGraph{T}(MetaGraph(n))

ProductGraph{T}(g₁::G1, g₂::G2) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, MetaGraph}} = ProductGraph{T}(product_graph(g₁, g₂, T))

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
function build_w_to_v₁v₂_pair_map(v₁v₂_pair_to_w::SparseMatrixCSC{Int, Int})::Vector{Tuple{Int, Int}}
    v₁v₂_pairs = findall(!iszero, v₁v₂_pair_to_w) # vector of CartesianIndices
    return [(w[1], w[2]) for w in v₁v₂_pairs] # vector of tuples
end

"""
compute the product graph adjacency matrix and mappings (in each direction) for relating the source graphs and product graph
"""
function product_graph_matrix_and_maps(g₁::MetaGraph, g₂::MetaGraph, type::Type{T})::Tuple{SparseMatrixCSC{Bool, Int}, SparseMatrixCSC{Int, Int}, Vector{Tuple{Int, Int}}} where T <: AbstractProductGraph
    
    # get (v₁ ∈ g₁, v₂ ∈ g₂) <--> w ∈ g₁xg₂ mappings
    v₁v₂_pair_to_w = build_v₁v₂_pair_to_w_map(g₁, g₂)
    w_to_v₁v₂_pair = build_w_to_v₁v₂_pair_map(v₁v₂_pair_to_w)

    # pre-allocate adj matrix
    n_g₁xg₂ = length(w_to_v₁v₂_pair)
    A = spzeros(Bool, n_g₁xg₂, n_g₁xg₂)

    # loop over pairs of vertices in the product graph
    for wᵢ = 1:n_g₁xg₂
        # product graph vertex wᵢ = (u₁, u₂)
        #    with u₁ ∈ g₁, u₂ ∈ g₂
        u₁, u₂ = w_to_v₁v₂_pair[wᵢ]
        for wⱼ = (wᵢ + 1):n_g₁xg₂
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
                    A[wᵢ, wⱼ] = 1
                end
            # is there a common non-adjacency? (d-edge) "d" = disconnected.
            #   (only relevant to factor graph) 
            #   see "Subgraph Matching Kernels for Attributed Graphs" 
            elseif (type == Factor) && !e₁ && !e₂ && (u₁ != v₁) && (u₂ != v₂)
                A[wᵢ, wⱼ] = 1
            end
        end
    end
    return A .|| A', v₁v₂_pair_to_w, w_to_v₁v₂_pair
end

@inline label_product_graph_edge(g₁::MetaGraph, u₁::Int, v₁::Int, ::Type{Direct}) = get_prop(g₁, u₁, v₁, :label)
@inline label_product_graph_edge(g₁::MetaGraph, u₁::Int, v₁::Int, ::Type{Factor}) = has_edge(g₁, u₁, v₁) ? get_prop(g₁, u₁, v₁, :label) : 0

"""
compute the product graph (of type `type`) between g₁ and g₂
"""
function product_graph(g₁::MetaGraph, g₂::MetaGraph, type::Type{T})::MetaGraph where T <: AbstractProductGraph
    A, _, w_to_v₁v₂_pair = product_graph_matrix_and_maps(g₁, g₂, type)
    n_g₁xg₂ = length(w_to_v₁v₂_pair)
    
    g₁xg₂ = MetaGraph(SimpleGraph(A))

    # label vertices
    for w = 1:n_g₁xg₂
        v₁, v₂ = w_to_v₁v₂_pair[w]
        set_prop!(g₁xg₂, w, :v₁v₂_pair, (v₁, v₂))
        set_prop!(g₁xg₂, w, :label, get_prop(g₁, v₁, :label))
    end
    # label edges
    for wᵢ = 1:n_g₁xg₂
        u₁, _ = w_to_v₁v₂_pair[wᵢ]
        for wⱼ = (wᵢ + 1):n_g₁xg₂
            v₁, _ = w_to_v₁v₂_pair[wⱼ]
            if A[wᵢ, wⱼ]
                edge_label = label_product_graph_edge(g₁, u₁, v₁, type)
                set_prop!(g₁xg₂, wᵢ, wⱼ, :label, edge_label)
            end
        end
    end
    return g₁xg₂
end

product_graph(g₁::MetaGraph, g₂::GraphMol, type::Symbol) = product_graph(g₁, MetaGraph(g₂), type)
product_graph(g₁::GraphMol, g₂::G, type::Symbol) where G <: Union{GraphMol, MetaGraph} = product_graph(MetaGraph(g₁), g₂, type)

"""
compute the adjacency matrix of the product graph (of type `type`) between g₁ and g₂, but do not explicitly construct the graph
"""
product_graph_adj_mat(g₁::MetaGraph, g₂::MetaGraph, type::Symbol)::SparseMatrixCSC{Bool, Int} = product_graph_matrix_and_maps(g₁, g₂, type)[1]

product_graph_adj_mat(g₁::MetaGraph, g₂::GraphMol, type::Symbol) = product_graph_adj_mat(g₁, MetaGraph(g₂), type)
product_graph_adj_mat(g₁::GraphMol, g₂::G, type::Symbol) where G <: Union{GraphMol, MetaGraph} = product_graph_adj_mat(MetaGraph(g₁), g₂, type)

