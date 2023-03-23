import GraphMakie.NetworkLayout: AbstractLayout

Base.@kwdef struct VizGraphKwargs
    savename::String = ""
    layout::Type{<: AbstractLayout} = GraphMakie.Spring
    highlight::Vector{Int}=Int[]
    highlight_color::Symbol=:orange
    highlight_width::Int=25
    node_alpha_mask::Vector{S} where {S <: Real}
    edge_alpha_mask::Vector{T} where {T <: Real}
end

function VizGraphKwargs(mol::GraphMol; kwargs...)
    return VizGraphKwargs(;
        node_alpha_mask=ones(length(mol.nodeattrs)),
        edge_alpha_mask=ones(length(mol.edgeattrs)),
        kwargs...
    )
end

function VizGraphKwargs(graph::AbstractMetaGraph; kwargs...)
    return VizGraphKwargs(;
        node_alpha_mask=ones(nv(graph)),
        edge_alpha_mask=ones(ne(graph)),
        kwargs...
    )
end

"""
returns the node labels for the indicated type of graph visualization
"""
viz_node_labels(graph::AbstractMetaGraph) = ["$v" for v in vertices(graph)]

function viz_node_labels(graph::ProductGraph)
    return ["$(get_prop(graph, v, :v₁v₂_pair))" for v in vertices(graph)]
end

"""
returns the edge labels for the indicated type of graph visualization
"""
function viz_edge_labels(graph::AbstractMetaGraph)::Vector{String}
    # general case
    edgelabel = ["$(get_prop(graph, e, :label))" for e in edges(graph)]
    replace!(edgelabel, "-1" => "a")
    return edgelabel
end

function viz_edge_labels(graph::ProductGraph{Modular})::Vector{String}
    # modular product graph specific: relabel d-type edges
    edgelabel = viz_edge_labels(MetaGraph(graph))
    replace!(edgelabel, "$D_EDGE" => "d")
    return edgelabel
end

"""
sets the colors for edges based on graph type
"""
function viz_edge_colors(::AbstractMetaGraph, edgelabel::Vector{String})::Vector{RGBA}
    # for general metagraphs
    edgestrokec = [RGBA(DARKGRAY..., 1) for _ in eachindex(edgelabel)]
    return edgestrokec
end

function viz_edge_colors(::ProductGraph{Modular}, edgelabel::Vector{String})::Vector{RGBA}
    # specific to modular product graphs
    edgestrokec = viz_edge_colors(MetaGraph(), edgelabel)
    edgestrokec[edgelabel .== "d"] .= RGBA(LIGHTGRAY..., 1)
    return edgestrokec
end

const DARKGRAY =
    parse.(Int, ["#A09C9C"[x] for x in [[2, 3], [4, 5], [6, 7]]], base=16) ./ 255

const LIGHTGRAY =
    parse.(Int, ["#E7E3E2"[x] for x in [[2, 3], [4, 5], [6, 7]]], base=16) ./ 255

function viz_node_colors(graph::AbstractMetaGraph)::Vector{RGBA}
    return [
        RGBA(
            (
                parse.(
                    Int,
                    [elements[atom].cpk_hex[x] for x in [[2, 3], [4, 5], [6, 7]]],
                    base=16
                ) ./ 255
            )...,
            1
        ) for atom in [props(graph, v)[:label] for v in vertices(graph)]
    ]
end

"""
apply the alpha mask to edge/node colors
"""
function alpha_mask!(color::Vector{RGBA}, mask::Vector{T}) where {T <: Real}
    for (i, α) in enumerate(mask)
        color[i] = RGBA(color[i].r, color[i].g, color[i].b, α)
    end
end

"""
apply the alpha mask to edge/node labels
"""
function alpha_mask!(label::Vector{String}, mask::Vector{T}) where {T <: Real}
    for (i, α) in enumerate(mask)
        if α == 0
            label[i] = ""
        end
    end
end

struct Molecular <: AbstractLayout{2, Float64} end

function (::Molecular)(graph::AbstractMetaGraph)::Vector{<: Point}
    return [Point2(get_prop(graph, v, :coords)) for v in vertices(graph)]
end

const graph_ax_kwargs = Dict(
    vcat(
        [
            :xticksvisible,
            :yticksvisible,
            :xticklabelsvisible,
            :yticklabelsvisible,
            :xgridvisible,
            :ygridvisible,
            :topspinevisible,
            :bottomspinevisible,
            :leftspinevisible,
            :rightspinevisible
        ] .=> false, 
        [:aspect => DataAspect()]
    )
)

GraphAxis(grid_pos::GridPosition; kwargs...) = 
    Axis(grid_pos; graph_ax_kwargs..., kwargs...)

function viz_graph!(
    ax::Axis, 
    graph::AbstractMetaGraph,
    opt::VizGraphKwargs=VizGraphKwargs(graph)
)
    subgraph = induced_subgraph(SimpleGraph(graph), opt.highlight)[1]
    if subgraph ≠ MetaGraph()
        GraphMakie.graphplot!(
            ax, 
            subgraph;
            node_attr=(; color=opt.highlight_color, markersize=opt.highlight_width),
            edge_attr=(; color=opt.highlight_color, linewidth=0.7 * opt.highlight_width),
            layout=_->opt.layout()(graph)[opt.highlight]
        )
    end

    node_labels = viz_node_labels(graph)
    alpha_mask!(node_labels, opt.node_alpha_mask)
    edge_labels = viz_edge_labels(graph)
    alpha_mask!(edge_labels, opt.edge_alpha_mask)
    node_colors = viz_node_colors(graph)
    alpha_mask!(node_colors, opt.node_alpha_mask)
    edge_colors = viz_edge_colors(graph, edge_labels)
    alpha_mask!(edge_colors, opt.edge_alpha_mask)
    GraphMakie.graphplot!(
        ax, 
        graph;
        layout=opt.layout(),
        node_attr=(; color=node_colors),
        edge_attr=(; color=edge_colors),
        nlabels=node_labels,
        elabels=edge_labels
    )
    return
end

"""
Visualize a molecular or product graph
"""
function viz_graph(graph::T; kwargs...) where {T <: Union{AbstractMetaGraph, GraphMol}}
    return viz_graph(graph, VizGraphKwargs(graph; kwargs...))
end

viz_graph(mol::GraphMol, opt::VizGraphKwargs) = viz_graph(MetaGraph(mol), opt)

function viz_graph(graph::AbstractMetaGraph, opt::VizGraphKwargs)
    fig = Figure()
    
    viz_graph!(Axis(fig[1, 1]; graph_ax_kwargs...), graph, opt)

    if opt.savename ≠ ""
        the_savename = opt.savename
        if split(the_savename, ".")[end] != "pdf"
            the_savename *= ".pdf"
        end
        CairoMakie.save(the_savename, fig)
    end

    return fig
end

export viz_graph
