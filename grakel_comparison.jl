using Graphs, MetaGraphs
import IOCapture

function grakel_adj_mat(g::AbstractGraph)::String
    adj_mat_str = ""
    adj_mat = adjacency_matrix(g)
    for row in eachrow(adj_mat)
        row_str = ""
        for el in row
            row_str *= "$el,"
        end
        adj_mat_str *= "[" * row_str * "],"
    end
    adj_mat_str = "[" * adj_mat_str * "]"
    return adj_mat_str
end

function grakel_node_labels(g::MetaGraph)::String
    node_attr_str = ""
    for v in vertices(g)
        node_attr_str *= "$(v-1):$(get_prop(g, v, :label)),"
    end
    node_attr_str = "{" * node_attr_str * "}"
    return node_attr_str
end

function grakel_edge_labels(g::MetaGraph)::String
    edge_label_str = ""
    for e in edges(g)
        i = src(e)
        j = dst(e)
        edge_label_str *= "($(i-1),$(j-1)):$(get_prop(g, i, j, :label)),"
        edge_label_str *= "($(j-1),$(i-1)):$(get_prop(g, i, j, :label)),"
    end
    edge_label_str = "{" * edge_label_str * "}"
    return edge_label_str
end

function grakel_graph(g::MetaGraph)::String
    graph_str = "grakel.Graph("
    graph_str *= grakel_adj_mat(g)
    graph_str *= ","
    graph_str *= "node_labels="
    graph_str *= grakel_node_labels(g)
    graph_str *= ","
    graph_str *= "edge_labels="
    graph_str *= grakel_edge_labels(g)
    graph_str *= ")"
    return graph_str
end

function grakel_script(kernel::String, n::Int, graphs::Vector{<: MetaGraph})::Vector{String}
    script_str = String[]
    push!(script_str, "import time")
    push!(script_str, "import grakel")
    push!(script_str, "G = []")
    for graph in graphs
        push!(script_str, "G.append($(grakel_graph(graph)))")
    end
    push!(script_str, "kernel=grakel.$kernel")
    push!(script_str, "n=$n")
    push!(script_str, "tic = time.time()")
    push!(script_str, "for i in range(n):")
    push!(script_str, "\tkernel.fit_transform(G)")
    push!(script_str, "btime=(time.time()-tic)/n")
    push!(script_str, "print(btime, G)")
    return script_str
end

function grakel_script(kernel::String, n::Int, g₁::MetaGraph, g₂::MetaGraph)::Vector{String}
    script_str = String[]
    push!(script_str, "import time")
    push!(script_str, "import grakel")
    push!(script_str, "g1=$(grakel_graph(g₁))")
    push!(script_str, "g2=$(grakel_graph(g₂))")
    push!(script_str, "kernel=grakel.$kernel")
    push!(script_str, "n=$n")
    push!(script_str, "tic = time.time()")
    push!(script_str, "for i in range(n):")
    push!(script_str, "\tkernel.fit([g1])")
    push!(script_str, "\tkernel.transform([g2])")
    push!(script_str, "btime=(time.time()-tic)/n")
    push!(script_str, "print(btime, kernel.transform([g2])[0][0])")
    return script_str
end

function grakel_compute(kernel::String, graphs...; n::Int=1000)#::Tuple{Float64, Float64}
    script = grakel_script(kernel, n, graphs...)
    file = tempname()
    open(file, "w") do f
        return write.(f, script .* ["\n"])
    end
    err = false
    printed = IOCapture.capture() do
        try
            return run(Cmd([
                "python3"
                file
            ]))
        catch exception
            err = true
        end
    end.output
    if err
        @error file readlines(file)
        error(printed)
    end
    # t, v = parse.(Float64, split(printed))
    # return t, v
    return printed
end
