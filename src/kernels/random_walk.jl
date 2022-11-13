"""
    kernel_score = rwk(adj_mat; l=n)
    kernel_score = rwk(g₁xg₂; l=n)
    kernel_score = rwk(g₁, g₂; l=n)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel (RWK) via their direct product graph.
"""
function rwk(adj_mat::AbstractMatrix; l::Int)::Int
    return sum(adj_mat^l)
end

function rwk(g₁xg₂::ProductGraph{Direct}; kwargs...)::Int
    return rwk(adjacency_matrix(g₁xg₂); kwargs...)
end

function rwk(A::AbstractMetaGraph, B::AbstractMetaGraph; kwargs...)::Int
    return rwk(product_graph_adjacency_matrix(Direct, A, B); kwargs...)
end

function rwk(A::GraphMol, B::AbstractMetaGraph; kwargs...)::Int
    return rwk(MetaGraph(A), B; kwargs...)
end

function rwk(A::Union{AbstractMetaGraph, GraphMol}, B::GraphMol; kwargs...)::Int
    return rwk(A, MetaGraph(B); kwargs...)
end
