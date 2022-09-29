"""
    kernel_score = random_walk_kernel(adj_mat, l)
    kernel_score = random_walk_kernel(dpg, l)
    kernel_score = random_walk_kernel(A, B, l)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their direct product graph.
"""
function random_walk_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

random_walk_kernel(dpg::AbstractGraph, l::Int) = random_walk_kernel(adjacency_matrix(dpg), l)
random_walk_kernel(A::AbstractGraph, B::AbstractGraph, l::Int) = random_walk_kernel(dpg_adj_mat(A, B), l)
random_walk_kernel(A::AbstractGraph, B::GraphMol, l::Int) = random_walk_kernel(A, MetaGraph(B), l)
random_walk_kernel(A::GraphMol, B::T, l::Int) where T <: Union{AbstractGraph, GraphMol} = random_walk_kernel(MetaGraph(A), B, l)
