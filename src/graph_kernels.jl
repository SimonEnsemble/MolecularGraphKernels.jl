"""
    kernel_score = graph_kernel(adj_mat, l)
    kernel_score = graph_kernel(g₁xg₂, l)
    kernel_score = graph_kernel(g₁, g₂, l, type)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their `type` product graph.
"""
function graph_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

graph_kernel(g₁xg₂::AbstractGraph, l::Int) = graph_kernel(adjacency_matrix(g₁xg₂), l)
graph_kernel(A::AbstractGraph, B::AbstractGraph, l::Int, type::Symbol) = graph_kernel(product_graph_adj_mat(A, B, type), l)
graph_kernel(A::AbstractGraph, B::GraphMol, l::Int, type::Symbol) = graph_kernel(A, MetaGraph(B), l, type)
graph_kernel(A::GraphMol, B::T, l::Int, type::Symbol) where T <: Union{AbstractGraph, GraphMol} = graph_kernel(MetaGraph(A), B, l, type)
