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
return the appropriate label for an edge in a product graph given its type
"""
function product_graph_edge_label(::Type{Direct}, g₁::MetaGraph, u₁::Int, v₁::Int)
    return get_prop(g₁, u₁, v₁, :label)
end

export Direct
