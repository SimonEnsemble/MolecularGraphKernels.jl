"""
    kernel_score = random_walk(adj_mat; l=n)
    kernel_score = random_walk(g₁xg₂; l=n)
    kernel_score = random_walk(g₁, g₂; l=n)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel (random_walk) via their direct product graph.
"""
function random_walk(adj_mat::AbstractMatrix; l::Int)::Int
    return sum(adj_mat^l)
end

function random_walk(g₁xg₂::ProductGraph{Direct}; kwargs...)::Int
    return random_walk(adjacency_matrix(g₁xg₂); kwargs...)
end

function random_walk(A::AbstractMetaGraph, B::AbstractMetaGraph; kwargs...)::Int
    return random_walk(GraphMatrix{Direct}(A, B); kwargs...)
end

function random_walk(A::GraphMol, B::AbstractMetaGraph; kwargs...)::Int
    return random_walk(MetaGraph(A), B; kwargs...)
end

function random_walk(A::Union{AbstractMetaGraph, GraphMol}, B::GraphMol; kwargs...)::Int
    return random_walk(A, MetaGraph(B); kwargs...)
end

export random_walk
