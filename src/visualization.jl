"""
returns the node labels for the indicated type of graph visualization
"""
viz_node_labels(graph::AbstractGraph) = viz_node_labels(typeof(graph), graph)
viz_node_labels(::Type{T}, graph::MetaGraph) where T <: MetaGraph = ["$v" for v in vertices(graph)]
viz_node_labels(::Type{T}, graph::ProductGraph) where T <: ProductGraph = ["$(get_prop(graph, v, :v₁v₂_pair))" for v in vertices(graph)]

"""
returns the edge labels for the indicated type of graph visualization
"""
function viz_edge_labels(graph::T)::Vector{String} where T <: AbstractMetaGraph
    # for general metagraphs
    edgelabels = ["$(get_prop(graph, e, :label))" for e in edges(graph)]
    replace!(edgelabels, "-1" => "a")
    return edgelabels
end

# modular product graph specific
function viz_edge_labels(graph::ProductGraph{Modular})::Vector{String}
    edgelabels = viz_edge_labels(MetaGraph(graph))
    replace!(edgelabels, "0" => "d")
    return edgelabels
end

"""
Visualize a molecular or product graph
"""
function viz_graph(graph::AbstractMetaGraph; savename::String="", C::Float64=0.1, layout_style::Symbol=:spring, node_borders::Bool=true)
    layout = (args...) -> 
        if layout_style == :spring
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
        nodelabel=viz_node_labels(graph),
        edgelabel=viz_edge_labels(graph),
        NODELABELSIZE=5.,
        EDGELABELSIZE=6.,
        NODESIZE=0.3 / sqrt(nv(graph)),
        nodestrokelw = node_borders ? 1. : 0.,
    )

    if savename ≠ ""
        draw(PDF(savename * ".pdf", 13cm, 13cm), plot)
    end

    return plot
end

viz_graph(graph::GraphMol; kwargs...) = viz_graph(MetaGraph(graph); kwargs...)
