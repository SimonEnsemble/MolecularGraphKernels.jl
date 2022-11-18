### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
begin
    import IOCapture, Pkg
    IOCapture.capture() do
        return Pkg.activate(".")
    end
    using BenchmarkTools, PlutoUI
    using Graphs, MetaGraphs, MolecularGraph
    using MolecularGraphKernels
    TableOfContents(; title="Update: 2022.11.21")
end

# ╔═╡ dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
begin
    g₁ = MetaGraph(smilestomol("NC=O"))
    g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))
end;

# ╔═╡ 177de406-366a-4356-941b-b8505f208eb9
md"""
# Benchmarking
"""

# ╔═╡ dd17bcee-c71b-4141-8620-1ab5f3c91df1
md"""
## Grakel Graph-Writer

Actually, wrote a function to build and execute an entire Grakel pairwise comparison with a specified kernel and return the output and average execution time.
"""

# ╔═╡ f47322d1-9347-4b34-b60d-efda4b9cbfa8
md"""
### Adjacency Matrix
"""

# ╔═╡ 72ba5ed1-d182-4968-89be-ce5a3a3138f5
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
end;

# ╔═╡ 81a684ab-ccbc-4cf7-b88e-e2ad3c3f914d
md"""
### Node Labels
"""

# ╔═╡ a43b9345-42a8-4bf3-b1d1-cbfc84295897
function grakel_node_labels(g::MetaGraph)::String
    node_attr_str = ""
    for v in vertices(g)
        node_attr_str *= "$(v-1):$(get_prop(g, v, :label)),"
    end
    node_attr_str = "{" * node_attr_str * "}"
    return node_attr_str
end;

# ╔═╡ 16fd2f0d-6b5d-4049-9708-9523f1cca4b1
md"""
### Edge Labels
"""

# ╔═╡ 708ff0fe-0c6e-41d4-98ec-ab801a220bca
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
end;

# ╔═╡ 121ae0ca-978f-4043-92b1-a7252dc07c28
md"""
### Graph
"""

# ╔═╡ 2f36aab7-7697-44f6-b8f2-7a37d2209c37
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
end;

# ╔═╡ 0af46daf-165e-4275-a0a2-b05ae34522a0
md"""
### Script
"""

# ╔═╡ bfff619e-632f-4a42-abfb-e12cce82aab9
function grakel_script(g₁::MetaGraph, g₂::MetaGraph, kernel::String, n::Int)::Vector{String}
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
end;

# ╔═╡ 1cdf885f-3c9b-4ef3-ab27-80e50f29aef2
md"""
An example script (4-length random_walk on g₁ and g₂, average time over 1000 runs):
"""

# ╔═╡ c59742b0-5bcf-4b5b-9c94-354f2f90b284
println.(grakel_script(g₁, g₂, "RandomWalk(p=4)", 1000));

# ╔═╡ f083d2e3-5c6f-40a2-8fe0-e1a2604a24e6
md"""
### Computation
"""

# ╔═╡ 716b38b3-80f2-4ed0-9140-f446740642b0
function grakel_compute(g₁, g₂, kernel; n::Int=1000)::Tuple{Float64, Float64}
    script = grakel_script(g₁, g₂, kernel, n)
    file = tempname()
    open(file, "w") do f
        return write.(f, script .* ["\n"])
    end
    printed = IOCapture.capture() do
        return run(Cmd([
            "python3"
            file
        ]))
    end.output
    t, v = parse.(Float64, split(printed))
    return t, v
end;

# ╔═╡ 4a4f168a-c277-4292-88bc-85b5b7d4205c
md"""
Running the example shown above:
"""

# ╔═╡ af9f0e95-22d5-4af1-a086-63717f76f69b
grakel_compute(g₁, g₂, "RandomWalk(p=4)")

# ╔═╡ 62a6f895-77f9-4f88-a936-1a25563de6b4
md"""
## Comparison
"""

# ╔═╡ 203208ae-169b-4f80-8ec4-3fe742bea001
md"""
### Random Walk
"""

# ╔═╡ cc3552ac-6327-47f4-ab23-42e010e8c34c
md"""
``A\times A``
"""

# ╔═╡ 25e5176a-5716-4cd3-bba3-4f27af205421
mgk_random_walk_aa = @btime random_walk(g₁, g₁; l=4)

# ╔═╡ 62902c49-3ba2-4cff-8f57-01bd032ce7c2
gkl_random_walk_aa1 = grakel_compute(g₁, g₁, "RandomWalk(p=4)")

# ╔═╡ ce542c4f-9b98-4c4e-b113-57ee8b91cb58
gkl_random_walk_aa2 = grakel_compute(g₁, g₁, "RandomWalkLabeled(p=4)")

# ╔═╡ 5731939b-033a-41a9-aec3-2668d7a8286b
@assert mgk_random_walk_aa == gkl_random_walk_aa1 ||
        mgk_random_walk_aa == gkl_random_walk_aa2

# ╔═╡ 388656e8-7ff2-4929-bbc4-43cba431603b
md"""
!!! note
	The grakel code is much slower, returns an incorrect value, and does not respect node labels.
"""

# ╔═╡ 58958382-7d8d-48d4-84f3-ea68f44604b4
md"""
``A\times B``
"""

# ╔═╡ 66a808cf-e28e-4e47-aef0-c460988f788e
mgk_random_walk_ab = @btime random_walk(g₁, g₂; l=4)

# ╔═╡ bb18ae2f-2aaf-43f5-8632-10cf0f5f0f12
gkl_random_walk_ab1 = grakel_compute(g₁, g₂, "RandomWalk(p=4)")

# ╔═╡ 00eae4d8-88f9-4191-b156-1a2362272ea7
gkl_random_walk_ab2 = grakel_compute(g₁, g₂, "RandomWalkLabeled(p=4)")

# ╔═╡ a92a9c1a-c276-4748-aa30-94856b0dabd6
@assert mgk_random_walk_ab == gkl_random_walk_ab1 ||
        mgk_random_walk_ab == gkl_random_walk_ab2

# ╔═╡ 39ec203d-75d6-4947-8af8-2cc09fe1a92c
md"""
!!! note
	The grakel code is much slower, returns an incorrect value, and does not respect node labels.
"""

# ╔═╡ 86578480-67da-4bde-b08e-774f50704333
md"""
### Connected CSI
"""

# ╔═╡ e9315881-5c88-4651-a468-352cd5c64de8
md"""
``A\times A``
"""

# ╔═╡ 185a82b7-6b4d-43e4-b428-88c2b0876399
mgk_ccsi_aa = @btime ccsi(g₁, g₁)

# ╔═╡ cfb0a592-2064-4f10-8745-84caa6301f93
gkl_ccsi_aa = grakel_compute(g₁, g₁, "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ db095762-bc79-416f-b568-4c824f274f66
@assert mgk_ccsi_aa == gkl_ccsi_aa[2]

# ╔═╡ f7c02678-0fa7-41de-9813-4745af529213
md"""
``A\times B``
"""

# ╔═╡ 5401eeb6-8895-4d46-864c-cca78476caca
mgk_ccsi_ab = @btime ccsi(g₁, g₂)

# ╔═╡ 722bbc80-1f7a-4bbf-8a5d-817cd722044d
gkl_ccsi_ab = grakel_compute(g₁, g₂, "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ a13fbc1e-8529-467b-99bc-d05fbea03828
@assert mgk_ccsi_ab == gkl_ccsi_ab[2]

# ╔═╡ 62cae695-8d6e-4eb7-b4f1-55beef266db8
md"""
### Connected Subgraph Matching 🚧
"""

# ╔═╡ af59b147-1200-4aa4-bb9f-abe88e58187b
md"""
!!! warning "In Progress"
	Need to implement weighted product graph computation.
"""

# ╔═╡ 50cbce84-8bc6-4ab9-9380-7d77bc2221f6
md"""
!!! warning "In Progress"
	Need to figure out how Grakel takes ``k_v`` and ``k_e``.
"""

# ╔═╡ 6cb9bef0-423a-42d5-9474-1faa2eb4b6b5
println.(
    grakel_script(
        g₁,
        g₂,
        "SubgraphMatching(k=999, lw='uniform', ke=lambda e: 1, kv=lambda v: 1)",
        1000
    )
);

# ╔═╡ c06de19a-282e-4de5-b15c-2702ef1daaf7
grakel_compute(
    g₁,
    g₂,
    "SubgraphMatching(k=999, lw='uniform', ke=lambda e: 1, kv=lambda v: 1)"
)

# ╔═╡ b8a45a85-f208-49cd-80ec-9dc70e36fae9
md"""
## Scaling vs. Nodes  🚩
"""

# ╔═╡ e110a84a-ae84-4145-af12-24cfe74d41e7
graphs =
    MetaGraph.(
        smilestomol.(
            [
                "Oc1c(c(O)cc(c1)CCCCC)[C@@H]2\\C=C(/CC[C@H]2\\C(=C)C)C" # cannabidiol
                "Oc1cc(cc(O)c1C/C=C(\\C)CC\\C=C(/C)C)CCCCC" # cannabigerol
            ]
        )
    )

# ╔═╡ 5c892fd7-8a8b-4f31-8a7c-2f5a1279915e
grakel_compute(graphs[1], graphs[2], "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ bb8c45a4-7146-4c43-b060-d6377d6556e5
@btime ccsi(graphs[1], graphs[2])

# ╔═╡ Cell order:
# ╠═9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
# ╠═dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
# ╟─177de406-366a-4356-941b-b8505f208eb9
# ╟─dd17bcee-c71b-4141-8620-1ab5f3c91df1
# ╟─f47322d1-9347-4b34-b60d-efda4b9cbfa8
# ╠═72ba5ed1-d182-4968-89be-ce5a3a3138f5
# ╟─81a684ab-ccbc-4cf7-b88e-e2ad3c3f914d
# ╠═a43b9345-42a8-4bf3-b1d1-cbfc84295897
# ╟─16fd2f0d-6b5d-4049-9708-9523f1cca4b1
# ╠═708ff0fe-0c6e-41d4-98ec-ab801a220bca
# ╟─121ae0ca-978f-4043-92b1-a7252dc07c28
# ╠═2f36aab7-7697-44f6-b8f2-7a37d2209c37
# ╟─0af46daf-165e-4275-a0a2-b05ae34522a0
# ╠═bfff619e-632f-4a42-abfb-e12cce82aab9
# ╟─1cdf885f-3c9b-4ef3-ab27-80e50f29aef2
# ╠═c59742b0-5bcf-4b5b-9c94-354f2f90b284
# ╟─f083d2e3-5c6f-40a2-8fe0-e1a2604a24e6
# ╠═716b38b3-80f2-4ed0-9140-f446740642b0
# ╟─4a4f168a-c277-4292-88bc-85b5b7d4205c
# ╠═af9f0e95-22d5-4af1-a086-63717f76f69b
# ╟─62a6f895-77f9-4f88-a936-1a25563de6b4
# ╟─203208ae-169b-4f80-8ec4-3fe742bea001
# ╟─cc3552ac-6327-47f4-ab23-42e010e8c34c
# ╠═25e5176a-5716-4cd3-bba3-4f27af205421
# ╠═62902c49-3ba2-4cff-8f57-01bd032ce7c2
# ╠═ce542c4f-9b98-4c4e-b113-57ee8b91cb58
# ╠═5731939b-033a-41a9-aec3-2668d7a8286b
# ╟─388656e8-7ff2-4929-bbc4-43cba431603b
# ╟─58958382-7d8d-48d4-84f3-ea68f44604b4
# ╠═66a808cf-e28e-4e47-aef0-c460988f788e
# ╠═bb18ae2f-2aaf-43f5-8632-10cf0f5f0f12
# ╠═00eae4d8-88f9-4191-b156-1a2362272ea7
# ╠═a92a9c1a-c276-4748-aa30-94856b0dabd6
# ╟─39ec203d-75d6-4947-8af8-2cc09fe1a92c
# ╟─86578480-67da-4bde-b08e-774f50704333
# ╟─e9315881-5c88-4651-a468-352cd5c64de8
# ╠═185a82b7-6b4d-43e4-b428-88c2b0876399
# ╠═cfb0a592-2064-4f10-8745-84caa6301f93
# ╠═db095762-bc79-416f-b568-4c824f274f66
# ╟─f7c02678-0fa7-41de-9813-4745af529213
# ╠═5401eeb6-8895-4d46-864c-cca78476caca
# ╠═722bbc80-1f7a-4bbf-8a5d-817cd722044d
# ╠═a13fbc1e-8529-467b-99bc-d05fbea03828
# ╟─62cae695-8d6e-4eb7-b4f1-55beef266db8
# ╟─af59b147-1200-4aa4-bb9f-abe88e58187b
# ╟─50cbce84-8bc6-4ab9-9380-7d77bc2221f6
# ╠═6cb9bef0-423a-42d5-9474-1faa2eb4b6b5
# ╠═c06de19a-282e-4de5-b15c-2702ef1daaf7
# ╟─b8a45a85-f208-49cd-80ec-9dc70e36fae9
# ╠═e110a84a-ae84-4145-af12-24cfe74d41e7
# ╠═5c892fd7-8a8b-4f31-8a7c-2f5a1279915e
# ╠═bb8c45a4-7146-4c43-b060-d6377d6556e5
