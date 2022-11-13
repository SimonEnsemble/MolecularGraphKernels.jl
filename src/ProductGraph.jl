"""
abstract type for product graphs
"""
abstract type AbstractProductGraph end

"""
concrete product graph types
"""
struct Modular <: AbstractProductGraph end
struct Direct <: AbstractProductGraph end
struct Weighted <: AbstractProductGraph end

"""
type-parameterized struct for product graphs
"""
struct ProductGraph{T <: AbstractProductGraph, U <: Real} <: AbstractMetaGraph{Int}
    graph::SimpleGraph{Int}
    vprops::Dict{Int, PropDict}
    eprops::Dict{SimpleEdge{Int}, PropDict}
    gprops::PropDict
    weightfield::Symbol
    defaultweight::U
    metaindex::MetaDict
    indices::Set{Symbol}
end

function ProductGraph{S}(
    x::Union{AbstractMatrix, AbstractGraph},
    weightfield::Symbol=:weight,
    defaultweight::U=1.0
) where {U <: Real, S <: AbstractProductGraph}
    return ProductGraph{S, U}(
        SimpleGraph(x),
        Dict{Int, PropDict}(),
        Dict{SimpleEdge{Int}, PropDict}(),
        PropDict(),
        weightfield,
        defaultweight,
        MetaDict(),
        Set{Symbol}()
    )
end

function ProductGraph{S}(
    x::MetaGraph,
    weightfield::Symbol=:weight,
    defaultweight::U=1.0
) where {U <: Real, S <: AbstractProductGraph}
    return ProductGraph{S, U}(
        x.graph,
        x.vprops,
        x.eprops,
        x.gprops,
        weightfield,
        defaultweight,
        x.metaindex,
        x.indices
    )
end

function ProductGraph{T}(
    g₁::Union{GraphMol, MetaGraph},
    g₂::Union{GraphMol, MetaGraph}
) where {T <: AbstractProductGraph}
    return product_graph(T, g₁, g₂)
end

Graphs.is_directed(::ProductGraph) = false
Graphs.is_directed(::Type{ProductGraph{T}}) where T = false
Graphs.is_directed(::Type{ProductGraph{T, U}}) where {T, U} = false

weighttype(::ProductGraph) = Int

set_props!(g::ProductGraph, e::SimpleEdge{Int}, d::Dict) = g.eprops[e] = d

props(g::ProductGraph, e::SimpleEdge{Int}) =
    if haskey(g.eprops, e)
        return g.eprops[e]
    elseif haskey(g.eprops, reverse(e))
        return g.eprops[reverse(e)]
    else
        error("No such edge: ", e)
    end
