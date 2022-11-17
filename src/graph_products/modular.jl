# definition for how d-type edges will be labeled
const D_EDGE = -2

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
return the appropriate label for an edge in a product graph given its type
"""
function product_graph_edge_label(::Type{Modular}, g₁::MetaGraph, u₁::Int, v₁::Int)
    return has_edge(g₁, u₁, v₁) ? get_prop(g₁, u₁, v₁, :label) : D_EDGE
end

GraphMatrix{Modular}(g₁::MetaGraph, g₂::MetaGraph) = 
    GraphMatrix{Modular}(product_graph_matrix_and_maps(Modular, g₁, g₂)[1])

export Modular
