struct GraphMatrix{U <: Union{AbstractGraph, AbstractProductGraph}, T <: Real} <: AbstractMatrix{T}
    matrix::Matrix{T}
end

import Base: getindex, size
Base.getindex(g::GraphMatrix, i...) = g.matrix[i...]
Base.size(g::GraphMatrix) = size(g.matrix)

function GraphMatrix(g::AbstractGraph)
    matrix = zeros(typeof(g.defaultweight), nv(g), nv(g))
    for edge in edges(g)
        matrix[src(edge), dst(edge)] = get_prop(g, edge, :label)
    end
    return GraphMatrix{typeof(g), Float64}(matrix + matrix')
end

GraphMatrix{T}(M::AbstractMatrix{U}) where {T, U} = GraphMatrix{T, U}(M)

function GraphMatrix{T}(
    g₁::AbstractGraph,
    g₂::GraphMol
) where T
    return GraphMatrix{T}(g₁, MetaGraph(g₂))
end

function GraphMatrix{T}(
    g₁::GraphMol,
    g₂::Union{GraphMol, AbstractGraph}
) where T
    return GraphMatrix{T}(MetaGraph(g₁), g₂)
end

export GraphMatrix, getindex, size
