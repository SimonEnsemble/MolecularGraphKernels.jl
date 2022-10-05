"""
    kernel_score = random_walk_graph_kernel(adj_mat, l)
    kernel_score = random_walk_graph_kernel(g₁xg₂, l)
    kernel_score = random_walk_graph_kernel(g₁, g₂, l, type)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their `type` product graph.
"""
function random_walk_graph_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

random_walk_graph_kernel(g₁xg₂::AbstractMetaGraph, l::Int) = random_walk_graph_kernel(adjacency_matrix(g₁xg₂), l)
random_walk_graph_kernel(A::G1, B::G2, l::Int, type::Type{T}) where {T <: AbstractProductGraph, G1,G2 <: AbstractMetaGraph} = random_walk_graph_kernel(ProductGraphMatrix{type}(A, B), l)
