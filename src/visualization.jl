function viz_graph(mol::GraphMol; savename=nothing)
    graph = MetaGraph(mol)
    locs_x, locs_y = spring_layout(graph, C=0.25)

    nodelabels = ["v<sub>$i</sub>" for i = 1:nv(graph)]

    edgelabels = ["$(get_prop(graph, e, :label))" for e in edges(graph)]
    edgelabels[isaromaticbond(mol)] .= "a"
    
    gp = gplot(
        graph, 
        locs_x, 
        locs_y, 
        nodestrokec=[RGB(rgb.r/255, rgb.g/255, rgb.b/255) for rgb in atomcolor(mol)], 
        nodesize=1.0, 
        NODESIZE=0.3 / sqrt(nv(graph)), 
        nodefillc=[RGBA(rgb.r/255, rgb.g/255, rgb.b/255, 0.075) for rgb in atomcolor(mol)], 
        NODELABELSIZE=5.0,
        EDGELABELSIZE=5.0,
        EDGELINEWIDTH=15.0/nv(graph),
        nodestrokelw=1,
        nodelabel=nodelabels,
        edgelinewidth=1,
        edgelabel=edgelabels
    )

    if ! isnothing(savename)
        draw(PDF(savename * ".pdf", 13cm, 13cm), gp)
    end
    return gp
end
