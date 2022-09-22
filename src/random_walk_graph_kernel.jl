"""
    kernel_score = fixed_length_rw_kernel(adj_mat, l)
    kernel_score = fixed_length_rw_kernel(dpg, l)
    kernel_score = fixed_length_rw_kernel(A, B, l)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel on their direct product graph.
"""
function fixed_length_rw_kernel(adj_mat::M, l::Int)::Int where M <: AbstractMatrix
    return sum(adj_mat ^ l)
end

fixed_length_rw_kernel(dpg::G, l::Int) where G <: AbstractGraph = fixed_length_rw_kernel(adjacency_matrix(dpg), l)

fixed_length_rw_kernel(A::G, B::G, l::Int) where G <: AbstractGraph = fixed_length_rw_kernel(direct_product_graph(A, B), l)
