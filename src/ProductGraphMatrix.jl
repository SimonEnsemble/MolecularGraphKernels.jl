"""
type-parameterized struct for product graph adjacency matrices
"""
mutable struct ProductGraphMatrix{T <: AbstractProductGraph, U <: Real} <: AbstractSparseMatrixCSC{U, Int}
    m::Int                
    n::Int                
    colptr::Vector{Int}   
    rowval::Vector{Int}   
    nzval::Vector{U}      
end

function ProductGraphMatrix{T}(nv::Int) where T <: AbstractProductGraph
    A = spzeros(Bool, Int, nv, nv)
    return ProductGraphMatrix{T, Bool}(A.m, A.n, A.colptr, A.rowval, A.nzval)
end

size(m::ProductGraphMatrix{T, U}) where {T, U} = (m.m, m.n)
_checkbuffers(m::ProductGraphMatrix{T, Bool}) where T = (true; m)
getcolptr(m::ProductGraphMatrix{T, Bool}) where T = m.colptr
rowvals(m::ProductGraphMatrix{T, Bool}) where T = m.rowval
nonzeros(m::ProductGraphMatrix{T, Bool}) where T = m.nzval

function ProductGraphMatrix{T}(matrix::AbstractMatrix{Bool}) where T <: AbstractProductGraph
    @assert size(matrix)[1] == size(matrix)[2] "Adjacency matrix must be square!"
    pgm = ProductGraphMatrix{T}(size(matrix)[1])
    pgm .= matrix
    return pgm
end
ProductGraphMatrix{T}(g₁::G1, g₂::G2) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, AbstractMetaGraph}} = product_graph_adj_mat(g₁, g₂, T)
