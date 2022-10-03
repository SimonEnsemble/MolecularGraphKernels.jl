"""
type-parameterized struct for product graphs
"""
struct ProductGraph{T <: AbstractProductGraph, U <: Real} <: AbstractMetaGraph{Int}
    graph::SimpleGraph{Int}
    vprops::Dict{Int,PropDict}
    eprops::Dict{SimpleEdge{Int},PropDict}
    gprops::PropDict
    weightfield::Symbol
    defaultweight::U
    metaindex::MetaDict
    indices::Set{Symbol}
end

function ProductGraph(x, weightfield::Symbol, defaultweight::U) where U <: Real
    T = eltype(x)
    g = SimpleGraph(x)
    vprops = Dict{T,PropDict}()
    eprops = Dict{SimpleEdge{T},PropDict}()
    gprops = PropDict()
    metaindex = MetaDict()
    idxs = Set{Symbol}()
    return ProductGraph(g, vprops, eprops, gprops, weightfield, defaultweight, metaindex, idxs)
end

function ProductGraph{S}(x, weightfield::Symbol, defaultweight::U) where {U <: Real, S <: AbstractProductGraph}
    T = eltype(x)
    g = SimpleGraph(x)
    vprops = Dict{T, PropDict}()
    eprops = Dict{SimpleEdge{T}, PropDict}()
    gprops = PropDict()
    metaindex = MetaDict()
    idxs = Set{Symbol}()
    return ProductGraph{S}(g, vprops, eprops, gprops, weightfield, defaultweight, metaindex, idxs)
end

ProductGraph{T}(x) where T = ProductGraph{T}(x, :weight, 1.0)

function ProductGraph{T}(g::MetaGraph) where T <: AbstractProductGraph
    newg = SimpleGraph{Int}(g.graph)
    return MetaGraph(newg, g.defaultweight)
end

ProductGraph{T}(g₁::G1, g₂::G2) where {T <: AbstractProductGraph, G1,G2 <: Union{GraphMol, MetaGraph}} = product_graph(g₁, g₂, T)

is_directed(::ProductGraph) = false
weighttype(::ProductGraph) = Int
