"""
    kernel_score = random_walk(adj_mat, l)
    kernel_score = random_walk(g₁xg₂, l)
    kernel_score = random_walk(g₁, g₂, l)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their direct product graph.
"""
function random_walk(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat^l)
end

function random_walk(g₁xg₂::T, l::Int) where {T <: AbstractGraph}
    return random_walk(adjacency_matrix(g₁xg₂), l)
end
function random_walk(A::T, B::T, l::Int) where {T <: AbstractMetaGraph}
    return random_walk(product_graph_adjacency_matrix(Direct, A, B), l)
end
