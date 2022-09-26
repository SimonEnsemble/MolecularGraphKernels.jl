"""
    kernel_score = fixed_length_rw_kernel(adj_mat, l)
    kernel_score = fixed_length_rw_kernel(dpg, l)
    kernel_score = fixed_length_rw_kernel(A, B, l)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their direct product graph.
"""
function fixed_length_rw_kernel(adj_mat::AbstractMatrix, l::Int)::Int
    return sum(adj_mat ^ l)
end

fixed_length_rw_kernel(dpg::AbstractGraph, l::Int) = fixed_length_rw_kernel(adjacency_matrix(dpg), l)

fixed_length_rw_kernel(A::AbstractGraph, B::AbstractGraph, l::Int) = fixed_length_rw_kernel(dpg_adj_mat(A, B), l)
fixed_length_rw_kernel(A::AbstractGraph, B::GraphMol, l::Int) = fixed_length_rw_kernel(A, MetaGraph(B), l)
fixed_length_rw_kernel(A::GraphMol, B::T, l::Int) where T <: Union{AbstractGraph, GraphMol} = fixed_length_rw_kernel(MetaGraph(A), B, l)
