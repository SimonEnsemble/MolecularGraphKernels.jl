"""
type-parameterized struct for product graph adjacency matrices
"""
struct ProductGraphMatrix{T <: AbstractProductGraph, M <: AbstractMatrix}
    matrix::M
end
ProductGraphMatrix{T}(matrix::M) where {T <: AbstractProductGraph, M <: AbstractMatrix} = ProductGraphMatrix{T, M}(matrix)
ProductGraphMatrix{T}(g₁::G1, g₂::G2) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, MetaGraph}} = product_graph_adj_mat(g₁, g₂, T)
