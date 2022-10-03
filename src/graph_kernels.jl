"""
    kernel_score = random_walk_graph_kernel(adj_mat, l)
    kernel_score = random_walk_graph_kernel(g₁xg₂, l)
    kernel_score = random_walk_graph_kernel(g₁, g₂, l, type)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their `type` product graph.
"""
function random_walk_graph_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

random_walk_graph_kernel(A::ProductGraphMatrix{T}, l::Int) where T <: AbstractProductGraph = random_walk_graph_kernel(A.matrix, l)
random_walk_graph_kernel(A::ProductGraph{T}, l::Int) where T <: AbstractProductGraph = random_walk_graph_kernel(A.graph, l)
random_walk_graph_kernel(g₁xg₂::MetaGraph, l::Int) = random_walk_graph_kernel(adjacency_matrix(g₁xg₂), l)
random_walk_graph_kernel(A::G1, B::G2, l::Int, type::Type{T}) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, MetaGraph}} = random_walk_graph_kernel(ProductGraphMatrix{type}(A, B), l)
