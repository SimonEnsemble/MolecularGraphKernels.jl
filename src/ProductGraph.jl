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
    graph = SimpleGraph(x)
    vprops = Dict{Int, PropDict}()
    eprops = Dict{SimpleEdge{Int}, PropDict}()
    gprops = PropDict()
    metaindex = MetaDict()
    indices = Set{Symbol}()
    return ProductGraph{S, U}(
        graph,
        vprops,
        eprops,
        gprops,
        weightfield,
        defaultweight,
        metaindex,
        indices
    )
end

function ProductGraph{T}(
    g₁::Union{GraphMol, MetaGraph},
    g₂::Union{GraphMol, MetaGraph}
) where {T <: AbstractProductGraph}
    return product_graph(T, g₁, g₂)
end

is_directed(::ProductGraph) = false

weighttype(::ProductGraph) = Int

set_props!(g::ProductGraph, e::SimpleEdge{Int}, d::Dict) = g.eprops[e] = d

props(g::ProductGraph, e::SimpleEdge{Int}) = if haskey(g.eprops, e)
        return g.eprops[e]
    elseif haskey(g.eprops, reverse(e))
        return g.eprops[reverse(e)]
    else
        error("No such edge: ", e)
    end
