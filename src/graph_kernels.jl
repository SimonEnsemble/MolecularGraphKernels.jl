"""
    kernel_score = graph_kernel(adj_mat, l)
    kernel_score = graph_kernel(g₁xg₂, l)
    kernel_score = graph_kernel(g₁, g₂, l, type)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their `type` product graph.
"""
function graph_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

graph_kernel(A::ProductGraphMatrix{T}, l::Int) where T <: AbstractProductGraph = graph_kernel(A.matrix, l)
graph_kernel(A::ProductGraph{T}, l::Int) where T <: AbstractProductGraph = graph_kernel(A.graph, l)
graph_kernel(g₁xg₂::MetaGraph, l::Int) = graph_kernel(adjacency_matrix(g₁xg₂), l)
graph_kernel(A::G1, B::G2, l::Int, type::Type{T}) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, MetaGraph}} = graph_kernel(ProductGraphMatrix{type}(A, B), l)
