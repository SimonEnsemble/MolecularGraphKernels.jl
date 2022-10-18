Base.@kwdef struct VizGraphKwargs
    savename::String = ""
    C::Float64 = 0.1
    layout_style::Symbol = :spring
    node_alpha_mask::Vector{S} where {S <: Real}
    edge_alpha_mask::Vector{T} where {T <: Real}
end

VizGraphKwargs(mol::GraphMol; kwargs...) = VizGraphKwargs(; node_alpha_mask=ones(length(mol.nodeattrs)), edge_alpha_mask=ones(length(mol.edgeattrs)), kwargs...)

VizGraphKwargs(graph::AbstractMetaGraph; kwargs...) = VizGraphKwargs(; node_alpha_mask=ones(nv(graph)), edge_alpha_mask=ones(ne(graph)), kwargs...)

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
    replace!(edgelabel, "0" => "d")
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

"""
sets the gplot layout style
"""
function select_graph_layout(layout_style::Symbol, C::Float64)
    return (args...) -> if layout_style == :spring
        spring_layout(args...; C=C)
    elseif layout_style == :circular
        circular_layout(args...)
    elseif layout_style == :spectral
        spectral_layout(args...)
    else
        error("Invalid layout style: ", layout_style)
    end
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

"""
Visualize a molecular or product graph
"""
viz_graph(graph::T; kwargs...) where {T <: Union{AbstractMetaGraph, GraphMol}} = viz_graph(graph, VizGraphKwargs(graph; kwargs...))

viz_graph(mol::GraphMol, opt::VizGraphKwargs) = viz_graph(MetaGraph(mol), opt)

function viz_graph(graph::AbstractMetaGraph, opt::VizGraphKwargs=VizGraphKwargs())
    layout = select_graph_layout(opt.layout_style, opt.C)

    edgelabel = viz_edge_labels(graph)
    alpha_mask!(edgelabel, opt.edge_alpha_mask)

    edgestrokec = viz_edge_colors(graph, edgelabel)
    alpha_mask!(edgestrokec, opt.edge_alpha_mask)

    nodefillc = viz_node_colors(graph)
    alpha_mask!(nodefillc, opt.node_alpha_mask)

    nodelabel = viz_node_labels(graph)
    alpha_mask!(nodelabel, opt.node_alpha_mask)

    plot = gplot(
        graph;
        layout=layout,
        nodefillc=nodefillc,
        nodelabel=nodelabel,
        edgelabel=edgelabel,
        NODELABELSIZE=5.0,
        EDGELABELSIZE=6.0,
        NODESIZE=0.3 / sqrt(nv(graph)),
        edgestrokec=edgestrokec
    )

    if opt.savename ≠ ""
        draw(PDF(opt.savename * ".pdf", 13cm, 13cm), plot)
    end

    return plot
end
