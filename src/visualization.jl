function viz_graph(mol::MetaGraph; savename=nothing)
    locs_x, locs_y = spring_layout(mol, C=0.25)

    nodelabels = ["v<sub>$i</sub>" for i = 1:nv(mol)]

    edgelabels = ["$(get_prop(mol, e, :label))" for e in edges(mol)]
    replace!(edgelabels, "-1" => "a")
    
    colors = [RGB((rc[:cpk_colors][atomsymbol(atom)] ./ 255)...) for atom in [props(mol, v)[:label] for v in vertices(mol)]]
    gp = gplot(
        mol, 
        locs_x, 
        locs_y, 
        nodestrokec=colors, 
        nodesize=1.0, 
        NODESIZE=0.3 / sqrt(nv(mol)), 
        nodefillc=colors, 
        NODELABELSIZE=5.0,
        EDGELABELSIZE=5.0,
        EDGELINEWIDTH=15.0/nv(mol),
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

viz_graph(mol::GraphMol; kwargs...) = viz_graph(MetaGraph(mol); kwargs...)
