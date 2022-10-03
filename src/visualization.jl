"""
returns the node labels for the indicated type of graph visualization
"""
viz_node_labels(graph::MetaGraph, ::Type{MetaGraph}) = ["$v" for v in vertices(graph)]
viz_node_labels(graph::MetaGraph, ::Type{T}) where T <: AbstractProductGraph = ["$(get_prop(graph, v, :v₁v₂_pair))" for v in vertices(graph)]

"""
returns the edge labels for the indicated type of graph visualization
"""
function viz_edge_labels(graph::MetaGraph, ::Type{MetaGraph})::Vector{String}
    edgelabels = ["$(get_prop(graph, e, :label))" for e in edges(graph)]
    replace!(edgelabels, "-1" => "a")
    return edgelabels
end

viz_edge_labels(graph::MetaGraph, ::Type{T}) where T <: AbstractProductGraph = viz_edge_labels(graph, MetaGraph)

function viz_edge_labels(graph::MetaGraph, ::Type{Modular})::Vector{String}
    edgelabels = viz_edge_labels(graph, MetaGraph)
    replace!(edgelabels, "0" => "d")
    return edgelabels
end

"""
Visualize a molecular or product graph
"""
function viz_graph(graph::MetaGraph, type::Type{T}=MetaGraph; savename::String="", C::Float64=0.1, layout_style=nothing) where T <: Union{MetaGraph, AbstractProductGraph}
    layout = (args...) -> 
        if isnothing(layout_style)
            spring_layout(args..., C=C)
        elseif layout_style == :circular
            circular_layout(args...)
        elseif layout_style == :spectral
            spectral_layout(args...)
        else
            error("Invalid layout style: ", layout_style)
        end

    plot = gplot(
        graph, 
        layout=layout,
        nodestrokec=RGB(0, 0, 0), 
        nodefillc=[RGB((rc[:cpk_colors][atomsymbol(atom)] ./ 255)...) for atom in [props(graph, v)[:label] for v in vertices(graph)]], 
        nodelabel=viz_node_labels(graph, type),
        edgelabel=viz_edge_labels(graph, type),
        NODELABELSIZE=5.,
        EDGELABELSIZE=6.,
        NODESIZE=0.3 / sqrt(nv(graph)),
        nodestrokelw=0.5, ##! only gives perimeter width EDGELINEWIDTH (default 1) or 0
    )

    if savename ≠ ""
        draw(PDF(savename * ".pdf", 13cm, 13cm), plot)
    end

    return plot
end

viz_graph(graph::GraphMol; kwargs...) = viz_graph(MetaGraph(graph); kwargs...)

viz_graph(g::ProductGraph{T}; kwargs...) where T <: AbstractProductGraph = viz_graph(g.graph, T; kwargs...)
