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
    TableOfContents(; title="Update: 2022.11.08")
end

# ╔═╡ 6a040045-0e23-4a1d-bb4e-dc8ac2eeb6c6
using Combinatorics

# ╔═╡ dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
begin
    g₁ = MetaGraph(smilestomol("NC=O"))
    g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))
end;

# ╔═╡ 8c0d01a2-e7e8-4c74-b766-1d8d4b900b78
md"""
# Connected CSI/SM Kernels 🚩
"""

# ╔═╡ 1051e362-7e28-44b4-9b64-7ff3a42c570e
md"""
## Algorithm
"""

# ╔═╡ d15a77f2-6d34-498e-8e4e-8d9f92a477c5
md"""
Algorithm ``SMK`` in the paper for calculating ``k_{CSI}`` with ``\lambda(C)=1``:

1. ``\text{while } |P| > 0 \text{ do }``

2. ``\text{ }\text{ }\text{ }\text{ }v\leftarrow\text{arbitrary element of }P``

3. ``\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C\cup {v}``

4. ``\text{ }\text{ }\text{ }\text{ }value\leftarrow value+1``

5. ``\text{ }\text{ }\text{ }\text{ }P^\prime=P\cap N(v)``

6. ``\text{ }\text{ }\text{ }\text{ }SMK(C^\prime,P^\prime)``

7. ``\text{ }\text{ }\text{ }\text{ }P\leftarrow P \setminus {v}``

To constrain the kernel to connected graphs, we thought we could change line 5 to be:

``P^\prime= \{u\in P\cap N(v):\exists k\in C^\prime\rightarrow l(u,k)\ne d\}``

but that leads to under-counting by eliminating too many candidate nodes.

I also tried, among other ideas,

``P^\prime= \{u\in P\cap N(v):\exists k\in C^\prime\cup P\rightarrow l(u,k)\ne d\}`` (overcounts)

Then it occurred to me that when Kriege wrote:

	only enumerate c-cliques by making sure that only vertices are added that are adjacent to a vertex in the current clique via at least one cedge

he *didn't* mean:

	only consider as candidate nodes those which extend the current clique while maintaining *c*-edge spanning

but rather:

	only add a node from the candidate nodes to the growing clique if it maintains  *c*-edge spanning

which means that the real focus of the change is at line 3!

Algorithm ``SMK`` as I have it now for calculating ``k_{CCSI}`` with ``\lambda(C)=1``:

1. ``\text{while } |P| > 0 \text{ do }``

2. ``\text{ }\text{ }\text{ }\text{ }v\leftarrow\text{arbitrary element of }P``

3. ``\text{ }\text{ }\text{ }\text{ }\text{if }C\cup {v}\in\mathcal{C}(G_p)\text{ do}``

4. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C\cup {v}``

5. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }value\leftarrow value+1``

6. ``\text{ }\text{ }\text{ }\text{ }\text{else}``

7. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C``

6. ``\text{ }\text{ }\text{ }\text{ }SMK(C^\prime,P\cap N(v))``

7. ``\text{ }\text{ }\text{ }\text{ }P\leftarrow P \setminus {v}``
"""

# ╔═╡ 15faa4c4-4253-454d-9e6a-e319dc5ddf38
md"""
### Implementation
"""

# ╔═╡ 0f0fbb49-a724-4804-ba8a-cc31734c17dd
function extends_clique(Gₚ, C, v)
	if C == []
		return true
	end
	for u in C
		if has_edge(Gₚ, u, v) && get_prop(Gₚ, u, v, :label) ≠ 0
			return true
		end
	end
	return false
end

# ╔═╡ a76c2d1f-faf2-453a-ab90-2b2bb72ecc13
function test_algo(g₁, g₂)
	value = 0
	Gₚ = ProductGraph{Modular}(g₁, g₂)
	Vₚ = collect(vertices(Gₚ))
	cliques = []
	
	function kernel(C, P)
		while length(P) > 0
			v = first(P)
			if extends_clique(Gₚ, C, v)
				C′ = union(C, v)
				push!(cliques, C′)
				value += 1
			else
				C′ = C
			end
			# @info "" v C′ P
			kernel(C′, intersect(P, neighbors(Gₚ, v)))
			P = setdiff(P, [v])
		end
	end

	kernel([], Vₚ)
	return value, cliques
end

# ╔═╡ 941c8f5b-5ada-4aa5-97a3-8dd0cf8b02e6
md"""
### Verification
"""

# ╔═╡ 89b4dfbd-d00e-4fc6-b66b-7a7fd007f3b8
viz_graph(MetaGraph(ProductGraph{Modular}(g₁, g₂)))

# ╔═╡ 9b411ba3-9368-4ad1-96bc-0af3ad389b07
test_algo(g₁, g₂)

# ╔═╡ f15440b5-6933-4dfe-8068-7bdee34ebf2b
md"""
!!! ok "This Is Correct... Right?"
	Are we all in agreement that each of these subgraphs is a valid and correct one, and that there are no others?
	I.e., that this algorithm definitely returns the correct result in this case?
"""

# ╔═╡ c585fd94-3cac-429e-8b91-66c56ce8433c
md"""
!!! danger "Wrong!"
	Works for the simple example, but a larger one shows a duplication problem...
"""

# ╔═╡ 177de406-366a-4356-941b-b8505f208eb9
md"""
# Benchmarking 🚧
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
An example script (4-length RWK on g₁ and g₂, average time over 1000 runs):
"""

# ╔═╡ c59742b0-5bcf-4b5b-9c94-354f2f90b284
println.(grakel_script(g₁, g₂, "RandomWalk(p=4)", 1000));

# ╔═╡ f083d2e3-5c6f-40a2-8fe0-e1a2604a24e6
md"""
### Computation
"""

# ╔═╡ 716b38b3-80f2-4ed0-9140-f446740642b0
function grakel_compute(g₁, g₂, kernel; n::Int=1000)::Tuple{Float64,Float64}
	script = grakel_script(g₁, g₂, kernel, n)
	file = tempname()
	open(file, "w") do f
		write.(f, script .* ["\n"])
	end
	printed = IOCapture.capture() do
		run(Cmd([
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
mgk_rwk_aa = @btime random_walk(g₁, g₁; l=4)

# ╔═╡ 62902c49-3ba2-4cff-8f57-01bd032ce7c2
gkl_rwk_aa1 = grakel_compute(g₁, g₁, "RandomWalk(p=4)")

# ╔═╡ ce542c4f-9b98-4c4e-b113-57ee8b91cb58
gkl_rwk_aa2 = grakel_compute(g₁, g₁, "RandomWalkLabeled(p=4)")

# ╔═╡ 5731939b-033a-41a9-aec3-2668d7a8286b
@assert mgk_rwk_aa == gkl_rwk_aa1 || mgk_rwk_aa == gkl_rwk_aa2

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
mgk_rwk_ab = @btime random_walk(g₁, g₂; l=4)

# ╔═╡ bb18ae2f-2aaf-43f5-8632-10cf0f5f0f12
gkl_rwk_ab1 = grakel_compute(g₁, g₂, "RandomWalk(p=4)")

# ╔═╡ 00eae4d8-88f9-4191-b156-1a2362272ea7
gkl_rwk_ab2 = grakel_compute(g₁, g₂, "RandomWalkLabeled(p=4)")

# ╔═╡ a92a9c1a-c276-4748-aa30-94856b0dabd6
@assert mgk_rwk_ab == gkl_rwk_ab1 || mgk_rwk_ab == gkl_rwk_ab2

# ╔═╡ 39ec203d-75d6-4947-8af8-2cc09fe1a92c
md"""
!!! note
	The grakel code is much slower, returns an incorrect value, and does not respect node labels.
"""

# ╔═╡ ff85644a-b58f-40b1-8624-0888451c0bb4
md"""
### Common Subgraph Isom.
"""

# ╔═╡ 68beef00-9ff8-480c-8867-0fc1085adde9
md"""
``A\times A``
"""

# ╔═╡ f9e900aa-3143-4f34-a559-ad7b40a5909a
@btime common_subgraph_isomorphism(g₁, g₁)

# ╔═╡ a78a1080-e1e6-4912-ac3b-4911f9e15564
@btime common_subgraph_isomorphism(g₁, g₁; λ=length)

# ╔═╡ fadbec05-4964-4dc4-8837-c563c81c24bf
md"""
``A\times B``
"""

# ╔═╡ d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
@btime common_subgraph_isomorphism(g₁, g₂) # default: weight = 1 for all cliques

# ╔═╡ a5494cac-b8b8-4985-9606-8d3b6999f2f1
@btime common_subgraph_isomorphism(g₁, g₂; λ=length) # weight = clique size

# ╔═╡ 2362ad7b-e271-49d6-8535-f1b767dbef0a
md"""
!!! note "Not Implemented In Grakel"
	Grakel doesn't have the CSI kernel, only the connected variant.
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
mgk_ccsi_aa = @btime common_subgraph_isomorphism(g₁, g₁; c_cliques=true)

# ╔═╡ cfb0a592-2064-4f10-8745-84caa6301f93
gkl_ccsi_aa = grakel_compute(g₁, g₁, "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ db095762-bc79-416f-b568-4c824f274f66
@assert mgk_ccsi_aa == gkl_ccsi_aa[2]

# ╔═╡ f7c02678-0fa7-41de-9813-4745af529213
md"""
``A\times B``
"""

# ╔═╡ 5401eeb6-8895-4d46-864c-cca78476caca
mgk_ccsi_ab = @btime common_subgraph_isomorphism(g₁, g₂; c_cliques=true)

# ╔═╡ 722bbc80-1f7a-4bbf-8a5d-817cd722044d
gkl_ccsi_ab = grakel_compute(g₁, g₂, "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ a13fbc1e-8529-467b-99bc-d05fbea03828
@assert mgk_ccsi_ab == gkl_ccsi_ab[2]

# ╔═╡ 62cae695-8d6e-4eb7-b4f1-55beef266db8
md"""
### Subgraph Matching
"""

# ╔═╡ af59b147-1200-4aa4-bb9f-abe88e58187b
md"""
!!! warning "In Progress"
	Need to implement weighted product graph computation.
"""

# ╔═╡ 7c3765cb-ba0f-475a-a675-6c26eb681f85
md"""
!!! note "Not Implemented In Grakel"
	Grakel doesn't have the regular SM kernel (only the connected variation).
"""

# ╔═╡ a41af8f2-2817-4419-897b-f203a95300ed
md"""
### Connected SM
"""

# ╔═╡ fc49c8e4-2ec1-4756-a572-4bb55c7cdf3c
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
println.(grakel_script(g₁, g₂, "SubgraphMatching(k=999, lw='uniform', ke=lambda e: 1, kv=lambda v: 1)", 1000));

# ╔═╡ c06de19a-282e-4de5-b15c-2702ef1daaf7
grakel_compute(g₁, g₂, "SubgraphMatching(k=999, lw='uniform', ke=lambda e: 1, kv=lambda v: 1)")

# ╔═╡ b8a45a85-f208-49cd-80ec-9dc70e36fae9
md"""
## Scaling vs. Nodes  ✖
"""

# ╔═╡ e110a84a-ae84-4145-af12-24cfe74d41e7
graphs = MetaGraph.(smilestomol.([
	"O=C(C)Oc1ccccc1C(=O)O"
	"O=C1C[C@@H](C\\C=C1\\C)C(C)=C"
	"Oc1c(c(O)cc(c1)CCCCC)[C@@H]2\\C=C(/CC[C@H]2\\C(=C)C)C"
]))

# ╔═╡ 655a0fff-5dab-44c0-b94c-06a4a583e8d8
test_algo(graphs[1], graphs[2])

# ╔═╡ 87b2d39d-733e-4d07-b26d-f79cebb24e9d
@btime common_subgraph_isomorphism(graphs[1], graphs[2]; c_cliques=true)

# ╔═╡ 5c892fd7-8a8b-4f31-8a7c-2f5a1279915e
grakel_compute(graphs[1], graphs[2], "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ 5e66ee62-40d3-4b79-8ce4-84ba1b8af658
md"""
!!! danger "Problem"
	There is an obvious problem with the connected CSI code...
	Inspecting the cliques visited, some are duplicates!
"""

# ╔═╡ 79994dc7-6138-4e3b-b50e-1b2769f80eb2
md"""
# Cannabinoid Clustering 🚧

	1. employ CSI kernel to create Gram matrix for cannabanoids ✔

	2. use diffusion map w/ kernel matrix to embed into 2D space ✔

	3. color points according to protein 🌟

	4. how much time does it take to compute the Gram matrix? ✔

"""

# ╔═╡ ce79c1aa-1d19-4751-9199-35535e094c67
md"""
!!! note
	Used CSI kernel w/ c-clique constraint (both `λ(C) = 1` and `λ(C) = length(C)`).
	This takes ca. 10 min on a single core @ 4.3 GHz.
"""

# ╔═╡ b4ad3fa6-660b-44ec-831f-b27a794a9ed1
md"""
!!! warning "In Progress"
	See diffmap.jl in private repo
"""

# ╔═╡ d48ad6f0-e709-400f-a6d1-f7ad6868ccf1
md"""
# Literature
"""

# ╔═╡ cbf219a0-184e-4c32-b1d0-9203e2392a1b
md"""
### Original CSI Paper

Levi (1973), via combination of:

- Corneil & Gotlieb (1970)
- Sirovich (1971)
"""

# ╔═╡ 527c224d-fa2d-4b01-b17e-501ab62f6167
md"""
### Tasks (CSI Kernels)

**Nothing!**

Closest is `CSI_GED`, which gives no application, only runtime comparison on ``k(G_1,G_1)`` for 100k molecules of average atom count 24, and 10k molecules of average atom count 45.  The time limit set was ca. 2 hours; that's 1.5 seconds per computation.  We can afford ca. 300 ms per computation.
"""

# ╔═╡ e3eb5365-4c29-4eda-b0b4-7b92d17430e8
md"""
### Tasks (MCS Algorithms)

- small molecule classification
- compound activity prediction / QSAR
- reaction mapping
- database searching
- small molecule-μRNA binding
- metabolite prediction
"""

# ╔═╡ 9eb41dc9-dfb0-489a-89d1-e9dc47f0327b
md"""
### Tasks (Graph Kernels, Generally)

- atomization energy
- pure-substance phase diagrams
- reaction yield
- RNA structure
- protein-protein docking
- solvation energy
- gene function
- metal surface adsorption energy
- QSAR
- pKa
- biomolecule receptor agonism
"""

# ╔═╡ 036bed20-ad56-4913-b9d2-473d4f6773a2
md"""
### Cool Stuff

- quantum computing
- stereochemistry
"""

# ╔═╡ a49a68fa-4032-4576-af10-a863f20bac8e
md"""
# Hmm... 🚩
"""

# ╔═╡ a3583838-05c8-4f28-b674-757a6416129a
md"""
!!! danger "Solution... But..."
	Accumulating the cliques and calculating the kernel value from the unique subset works... but that has an impact on speed.
	I don't understand the Grakel implementation--it's not clear, documented, or attributed.
"""

# ╔═╡ 7e7f3c05-5837-48a9-add6-b2566cfda88e
function calculate_value(G, cliques)
	value = 0
	for C in cliques
		node_product = prod(1 for c in C)
		edge_product = prod(1 for (s, d) in combinations(C, 2); init=1)
		C_weight = 1
		value += node_product * edge_product * C_weight
	end
	return value
end

# ╔═╡ 2022a1c3-e447-4105-a65c-4f4c94156e20
function test_algo2(g₁, g₂)
	value = 0
	Gₚ = ProductGraph{Modular}(g₁, g₂)
	Vₚ = collect(vertices(Gₚ))
	cliques = []
	
	function kernel(C, P)
		while length(P) > 0
			v = first(P)
			if extends_clique(Gₚ, C, v)
				C′ = union(C, v)
				push!(cliques, C′)
			else
				C′ = C
			end
			kernel(C′, intersect(P, neighbors(Gₚ, v)))
			P = setdiff(P, [v])
		end
	end

	kernel([], Vₚ)

	cliques = unique(cliques)
	value = calculate_value(Gₚ, cliques)
	
	return value, cliques
end

# ╔═╡ 0cdf1306-3a8a-46ee-b384-0770fb7e9dba
@btime test_algo2(graphs[1], graphs[2])

# ╔═╡ bba0313d-820f-49f7-aaac-bfe3304bbb21
@btime test_algo2(g₁, g₂)

# ╔═╡ db4bae15-d7d1-42fd-848f-f0e99c4f9a96
grakel_compute(graphs[1], graphs[2], "SubgraphMatching(k=999, lw='uniform')")

# ╔═╡ Cell order:
# ╠═9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
# ╠═dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
# ╟─8c0d01a2-e7e8-4c74-b766-1d8d4b900b78
# ╟─1051e362-7e28-44b4-9b64-7ff3a42c570e
# ╟─d15a77f2-6d34-498e-8e4e-8d9f92a477c5
# ╟─15faa4c4-4253-454d-9e6a-e319dc5ddf38
# ╠═0f0fbb49-a724-4804-ba8a-cc31734c17dd
# ╠═a76c2d1f-faf2-453a-ab90-2b2bb72ecc13
# ╟─941c8f5b-5ada-4aa5-97a3-8dd0cf8b02e6
# ╠═89b4dfbd-d00e-4fc6-b66b-7a7fd007f3b8
# ╠═9b411ba3-9368-4ad1-96bc-0af3ad389b07
# ╟─f15440b5-6933-4dfe-8068-7bdee34ebf2b
# ╠═655a0fff-5dab-44c0-b94c-06a4a583e8d8
# ╟─c585fd94-3cac-429e-8b91-66c56ce8433c
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
# ╟─ff85644a-b58f-40b1-8624-0888451c0bb4
# ╟─68beef00-9ff8-480c-8867-0fc1085adde9
# ╠═f9e900aa-3143-4f34-a559-ad7b40a5909a
# ╠═a78a1080-e1e6-4912-ac3b-4911f9e15564
# ╟─fadbec05-4964-4dc4-8837-c563c81c24bf
# ╠═d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
# ╠═a5494cac-b8b8-4985-9606-8d3b6999f2f1
# ╟─2362ad7b-e271-49d6-8535-f1b767dbef0a
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
# ╟─7c3765cb-ba0f-475a-a675-6c26eb681f85
# ╟─a41af8f2-2817-4419-897b-f203a95300ed
# ╟─fc49c8e4-2ec1-4756-a572-4bb55c7cdf3c
# ╟─50cbce84-8bc6-4ab9-9380-7d77bc2221f6
# ╠═6cb9bef0-423a-42d5-9474-1faa2eb4b6b5
# ╠═c06de19a-282e-4de5-b15c-2702ef1daaf7
# ╟─b8a45a85-f208-49cd-80ec-9dc70e36fae9
# ╠═e110a84a-ae84-4145-af12-24cfe74d41e7
# ╠═87b2d39d-733e-4d07-b26d-f79cebb24e9d
# ╠═5c892fd7-8a8b-4f31-8a7c-2f5a1279915e
# ╠═5e66ee62-40d3-4b79-8ce4-84ba1b8af658
# ╟─79994dc7-6138-4e3b-b50e-1b2769f80eb2
# ╟─ce79c1aa-1d19-4751-9199-35535e094c67
# ╟─b4ad3fa6-660b-44ec-831f-b27a794a9ed1
# ╟─d48ad6f0-e709-400f-a6d1-f7ad6868ccf1
# ╟─cbf219a0-184e-4c32-b1d0-9203e2392a1b
# ╟─527c224d-fa2d-4b01-b17e-501ab62f6167
# ╟─e3eb5365-4c29-4eda-b0b4-7b92d17430e8
# ╟─9eb41dc9-dfb0-489a-89d1-e9dc47f0327b
# ╟─036bed20-ad56-4913-b9d2-473d4f6773a2
# ╠═a49a68fa-4032-4576-af10-a863f20bac8e
# ╟─a3583838-05c8-4f28-b674-757a6416129a
# ╠═6a040045-0e23-4a1d-bb4e-dc8ac2eeb6c6
# ╠═7e7f3c05-5837-48a9-add6-b2566cfda88e
# ╠═2022a1c3-e447-4105-a65c-4f4c94156e20
# ╠═0cdf1306-3a8a-46ee-b384-0770fb7e9dba
# ╠═bba0313d-820f-49f7-aaac-bfe3304bbb21
# ╠═db4bae15-d7d1-42fd-848f-f0e99c4f9a96
