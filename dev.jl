### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
begin
    import IOCapture, Pkg
    IOCapture.capture() do
        Pkg.activate(".")
        Pkg.resolve()
        Pkg.add.(
            [
                "Documenter"
                "Aqua"
                "BenchmarkTools"
                "ProfileCanvas"
                "PlutoLinks"
                "PlutoTeachingTools"
                "PlutoTest"
                "PlutoUI"
            ]
        )
        return
    end
    using Documenter
    using Aqua, BenchmarkTools, ProfileCanvas
    using PlutoLinks, PlutoTeachingTools, PlutoTest, PlutoUI
    TableOfContents()
end

# ╔═╡ 037a4025-ad55-4fb5-a14b-8b2da83db1c7
@revise using MolecularGraphKernels

# ╔═╡ e9f5391d-4832-440e-b61c-357daf275332
using Graphs, MetaGraphs, MolecularGraph, SparseArrays

# ╔═╡ cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
md"""
# Dev Notebook Header
"""

# ╔═╡ 9188ef1e-16fe-4a79-8ba6-2b0e907d743a
md"""
### Local Development Preamble

Gives access to many tools useful for developing local packages.

Should be run from the root of the local package repo under development.
"""

# ╔═╡ fb64efc5-e959-401f-96d1-464de7d47547
md"""
### Local Package Revision

Invokes file monitoring for the local package (by name).
"""

# ╔═╡ 971586d9-266b-4dfd-97d6-dc3aed449600
md"""
### Additional Dependencies

External packages needed within the notebook for the specific project(s) at hand.
"""

# ╔═╡ f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
md"""
# Development Code
"""

# ╔═╡ 5699f8a5-11d6-453d-a867-8330134d080f
begin
    import Base.Multimedia.display
    Base.Multimedia.display(mol::GraphMol) = HTML(drawsvg(mol, 250, 250))
end

# ╔═╡ 5fec82c3-99fe-4ff0-aacd-7af622f07291
md"""
## Input Structures
"""

# ╔═╡ 1f2f45f6-57ba-4c29-845f-05685ceb299a
begin
    mol₁ = smilestomol("c1(c2)cscc1cc(c3)c2ccc3")
    mol₂ = smilestomol("c1(o2)cscc1oc(c3)c2ccc3")
    g₁ = MetaGraph(mol₁)
    g₂ = MetaGraph(mol₂)
    display.([mol₁, mol₂])
end

# ╔═╡ f2c86d66-b801-4192-965c-b0b82a5c603a
viz_graph.([g₁, g₂]; layout_style=:graphmol)

# ╔═╡ 41c83665-9cff-43c1-912f-3d820d682e09
mpg = ProductGraph{Modular}(g₁, g₂)

# ╔═╡ b96dee5e-6c6b-4a1e-938a-5a9b42a96c3b
viz_graph(mpg; layout_style=:spectral)

# ╔═╡ 39326496-e4dc-4b32-b538-feaa47066982
imsgs = isomorphic_subgraphs(mpg)

# ╔═╡ a9c0e399-397e-43a7-a4ad-19af2c650d71
viz_graph(imsgs[2]; layout_style=:spring)

# ╔═╡ 87aa9631-ccef-4532-a06b-0aaee425d908
begin
    local n = 1
    dg = deepcopy(imsgs[n])
    for e in edges(imsgs[n])
        if get_prop(imsgs[n], e, :label) == 0
            rem_edge!(dg, e)
        end
    end
    viz_graph(dg; layout_style=:circular)
end

# ╔═╡ 53a7b87c-e676-4f42-81dc-dc38450078d1
md"""
#### take above graph, extract node pair labels, alpha mask original graphs to only show c-connected subgraphs
"""

# ╔═╡ cf936eb7-020f-4e84-a284-8c7098d40cde
g₁_nodes = [get_prop(dg, v, :v₁v₂_pair)[1] for v in vertices(dg)]

# ╔═╡ a60742e4-cdbc-4cb2-a1fd-9fa03ce5cf64
g₂_nodes = [get_prop(dg, v, :v₁v₂_pair)[2] for v in vertices(dg)]

# ╔═╡ 9875f635-1528-4b55-a3a7-fce2e6bbd866
g₁_edges = [Graphs.SimpleEdge(g₁_nodes[src(e)], g₁_nodes[dst(e)]) for e in edges(dg)]

# ╔═╡ ff6ddb75-f9fd-4279-8c95-d7be0ce3e812
g₂_edges = [Graphs.SimpleEdge(g₂_nodes[src(e)], g₂_nodes[dst(e)]) for e in edges(dg)]

# ╔═╡ 0a97908e-bf4b-41a1-8cf6-e7905325393a
α₀ = 0.075

# ╔═╡ f8fd1d2c-c1e5-4173-b704-8c345096b059
g₁_node_alpha_mask = [v ∈ g₁_nodes ? 1 : α₀ for v in vertices(g₁)]

# ╔═╡ f0bbbfe6-b799-4be7-8afb-61ef65f1a37e
g₂_node_alpha_mask = [v ∈ g₂_nodes ? 1 : α₀ for v in vertices(g₂)]

# ╔═╡ 983b5c21-7618-483f-abe3-32efd6452130
g₁_edge_alpha_mask = [e ∈ g₁_edges || reverse(e) ∈ g₁_edges ? 1 : α₀ for e in edges(g₁)]

# ╔═╡ d57a3747-dcc6-46b7-946f-1828cae62fa2
g₂_edge_alpha_mask = [e ∈ g₂_edges || reverse(e) ∈ g₂_edges ? 1 : α₀ for e in edges(g₂)]

# ╔═╡ 7663b5a0-d0a9-4a72-9d75-745a91160737
viz_graph(
    g₁;
    node_alpha_mask=g₁_node_alpha_mask,
    edge_alpha_mask=g₁_edge_alpha_mask,
    layout_style=:graphmol
)

# ╔═╡ 575ccd4b-50f0-405d-9a39-48bc1265512e
viz_graph(
    g₂;
    node_alpha_mask=g₂_node_alpha_mask,
    edge_alpha_mask=g₂_edge_alpha_mask,
    layout_style=:graphmol
)

# ╔═╡ 7dcc01ec-f087-467f-be59-e5404d44946f
md"""
## Subgraph Matching Kernel λ
"""

# ╔═╡ d55eaec7-0e62-44c5-a139-986b12e731fe
md"""
### *
"""

# ╔═╡ cd3a297a-d05f-48d0-81a6-c37e1dfa3777
@btime subgraph_matching(g₁, g₂; λ=length)

# ╔═╡ 024ed879-6da4-4770-81ca-9dc4cf9b4112
md"""
## Max Clique Tinkering
"""

# ╔═╡ 6a323e50-ea97-41e1-8655-026ff5d73a00
md"""
Get the list of maximal cliques from the modular product graph (450 ms)
"""

# ╔═╡ 1806ea53-6738-4082-abd7-fc7e7ca5c463
max_cliques = maximal_cliques(mpg.graph)

# ╔═╡ 5d7efccd-9575-4145-bf7b-1c614539f733
md"""
 - Translate the MPG nodes to g₁ nodes
 - Sort each list
 - Find all unique sets
(100 ms)
"""

# ╔═╡ 10184e95-d9a8-47b6-a2f3-8c25a29bca8a
g₁_csis = unique([sort([get_prop(mpg, m, :v₁v₂_pair)[1] for m in mc]) for mc in max_cliques])

# ╔═╡ 0990663a-0970-4468-9e42-15e77c209a79
md"""
Sort the g₁ vertex subsets by size (100 μs)
"""

# ╔═╡ 5c0a8f21-36b9-45f5-80b9-6726c69dd313
sorted_g₁csis = g₁_csis[sortperm(length.(g₁_csis), rev=true)]

# ╔═╡ 00b1d071-8c5c-40c9-99a5-ae5b380e3151
md"""
Allow a given node to be present only in its largest subgraph (400 μs)
"""

# ╔═╡ 2b52cdea-7af8-4b16-97bd-5ad32538e7e3
begin
	used_nodes = falses(nv(g₁))
	maxmax_csis = falses(length(sorted_g₁csis))
	for (i, csi_indices) in enumerate(sorted_g₁csis)
		if any(used_nodes[csi_indices])
			continue
		else
			used_nodes[csi_indices] .= true
			maxmax_csis[i] = true
		end
	end
	sorted_g₁csis[maxmax_csis]
end

# ╔═╡ 6d85fe2b-e9aa-4f91-af77-4f291399bc7d
md"""
The result is the maximum common subgraph isomorphism.
"""

# ╔═╡ 7e090139-72db-4658-b4f6-dac9bad7b848
viz_graph(
	induced_subgraph(g₁, sorted_g₁csis[maxmax_csis][1])[1]; 
	layout_style=:graphmol
)

# ╔═╡ 3e18dc29-6a47-4b35-8c38-fe374484f0a0
md"""
### Efficiency

Seems like there should be a faster way to do at least some of this
"""

# ╔═╡ 7ee29c76-ecd8-4eef-a69d-506f0d171d36
md"""
The SM kernel as I implemented, w/ weight function ``\lambda(C)=1``:
"""

# ╔═╡ e5dc6300-3445-4f59-8131-1694a09dbbbe
@btime common_subgraph_isomorphism(g₁, g₂)

# ╔═╡ 8193aea0-7cbc-4c1c-bd10-e6992e7eee02
md"""
`Graphs.jl` has functions for counting [induced] subgraphisomorphisms... but they don't work?  
"""

# ╔═╡ b85c081c-220c-4a20-83b5-2d1204d5a999
import Graphs.Experimental: count_induced_subgraphisomorph, count_subgraphisomorph

# ╔═╡ d124ca26-8de5-4ac7-838e-dc37524af5c3
b = count_subgraphisomorph(
	SimpleGraph(g₁),
	SimpleGraph(g₂);
	vertex_relation=(v₁, v₂) -> get_prop(g₁, v₁, :label) == get_prop(g₂, v₂, :label),
	edge_relation=(e₁, e₂) -> get_prop(g₁, e₁, :label) == get_prop(g₂, e₂, :label)
)

# ╔═╡ 765c3a87-9bce-488f-afb3-7db637c0f6c2
c = count_induced_subgraphisomorph(
	SimpleGraph(g₁), 
	SimpleGraph(g₂);
	vertex_relation=(v₁, v₂) -> get_prop(g₁, v₁, :label) == get_prop(g₂, v₂, :label),
	edge_relation=(e₁, e₂) -> get_prop(g₁, e₁, :label) == get_prop(g₂, e₂, :label)
)

# ╔═╡ d9afe6e0-86ce-49c4-a842-73862b369a99
md"""
Find maximum commmon subgraph isomorphism size (MCSS) as implemented above:
"""

# ╔═╡ bbf82170-e0ee-49d9-8aad-22329b4f6f6f
function maximum_common_subgraph_size(g₁, g₂)
	max_cliques = maximal_cliques(
		SimpleGraph(product_graph_adjacency_matrix(Modular, g₁, g₂))
	)
	g₁_csis = unique([sort([get_prop(mpg, m, :v₁v₂_pair)[1] for m in mc]) for mc in max_cliques])
	sorted_g₁csis = g₁_csis[sortperm(length.(g₁_csis), rev=true)]
	used_nodes = falses(nv(g₁))
	maxmax_csis = falses(length(sorted_g₁csis))
	for (i, csi_indices) in enumerate(sorted_g₁csis)
		if any(used_nodes[csi_indices])
			continue
		else
			used_nodes[csi_indices] .= true
			maxmax_csis[i] = true
		end
	end
	return length(sorted_g₁csis[maxmax_csis][1])
end

# ╔═╡ e30b4937-361f-47e3-95c1-3200a1466c12
@btime maximum_common_subgraph_size(g₁, g₂)

# ╔═╡ 089eddca-c7b6-44f0-9133-681db448aa73
md"""
Optimized a bit:
"""

# ╔═╡ 465cc3ca-40ee-4b90-875c-a74c5e96f2d4
function max_csi_len(g₁, g₂)
	return maximum(
		length.(
			maximal_cliques(
				SimpleGraph(product_graph_adjacency_matrix(Modular, g₁, g₂))
			)
		)
	)
end

# ╔═╡ a476ff95-5a4f-457d-b3e7-bc9b2f22cd07
@btime max_csi_len(g₁, g₂)

# ╔═╡ e9736f26-ce94-4854-aecb-b2991be0e348
md"""
If this is what we want to use, `MolecularGraph` has something for this, too:
"""

# ╔═╡ 159e1fe2-eb50-42b8-b39f-c1d6ada69048
@btime length(disconnectedmcis(mol₁, mol₂).mapping)

# ╔═╡ 00e01fc4-7c8f-44ab-acf7-0e1593edd386
md"""
And the equivalent for the largest connected 
"""

# ╔═╡ d5a8b311-d276-42b1-bdb9-5a44981bf601
@btime length(connectedmcis(mol₁, mol₂).mapping)

# ╔═╡ 57d62a74-d13a-4333-8d90-cf0e69b0973c
md"""
And a "topological constraint" version that isn't explained very well at all:
"""

# ╔═╡ 4611beba-ab89-4d68-997d-9d84c166c90a
@bind θ Slider(0:10)

# ╔═╡ 1ab6df2d-502e-4aa8-a008-a24e25aef2d1
@btime length(tcmces(mol₁, mol₂, tolerance=θ).mapping)

# ╔═╡ e6a6832f-4f76-4885-9258-b02c1d5a91e4
md"""
### Idea

If ``C`` is a clique in ``\mathcal{G}`` and ``v`` is a node in ``\mathcal{V}(\mathcal{G})``, 

$\text{deg}(v) < |C| \rightarrow v \notin C_{max}$

where ``C_{max}`` is the maximum clique.

Therefore, I propose to sort the node indices by decreasing degree of the nodes and perform a tree search by the following algorithm:
"""

# ╔═╡ 19fd937a-20b4-4e82-9c08-7edcb4161e1f
md"""
1. ``i := 1`` and ``S := 0``
2. ``C := \{\}``
3. ``C^\prime := \{ v_i \}``
4. while ``C \ne C^\prime``:
 * 5. ``C := C^\prime``
 * 6. ``C^\prime := C\cup\bigcap_{v_j\in C}\mathcal{N}(v_j)``
7. ``S:=|C|``
8. ``i:=i+1``
9. if ``v_i\in\mathcal{G}\wedge\text{deg}(v_i)>S`` goto 2
"""

# ╔═╡ 37e534f1-e03e-47c3-9a6f-ee9cc5f10de5
function len_dmcis(g₁, g₂)
	G = product_graph_adjacency_matrix(Modular, g₁, g₂)
	deg = sum.(eachcol(G))
	sorted_nodes = sortperm(deg)
	S = 0
	all_C = Set[]
	for vᵢ in sorted_nodes
		if deg[vᵢ] ≤ S
			break
		end
		C = Set([])
		C′ = Set([vᵢ])
		while C ≠ C′
			C = C′
			C′ = union(
				C, Set(findall(reduce((x, y) -> x .& y, [G[:, vⱼ] for vⱼ in C])))
			)
		end
		S = length(C)
		push!(all_C, C′)
	end
	return S, all_C
end

# ╔═╡ 9b21ce88-3ae4-45e1-b324-61ae3faf2ae4
len_dmcis(g₁, g₂)

# ╔═╡ 38e99f16-b499-4020-b51c-462b811bf135
@btime length(unique([get_prop(mpg, v, :v₁v₂_pair)[1] for v in len_dmcis(g₁, g₂)[2][1]]))

# ╔═╡ 92db43c3-ae2c-4518-9c39-1403c887a315
@test max_csi_len(g₁, g₂) == length(unique([get_prop(mpg, v, :v₁v₂_pair)[1] for v in len_dmcis(g₁, g₂)[2][1]]))

# ╔═╡ 8b0473ee-47ab-4ad1-863e-edd2883407d0
viz_graph(
	induced_subgraph(
		g₁, 
		unique(
			[get_prop(mpg, v, :v₁v₂_pair)[1] for v in len_dmcis(g₁, g₂)[2][1]]
		)
	)[1], 
	layout_style=:graphmol
)

# ╔═╡ b41f3383-2105-4a92-8b7e-edfb00addae0
md"""
# Kernel Comp Time
"""

# ╔═╡ f2f56957-1001-4bfa-962d-2210c4e8ce67
md"""
### Calculation
"""

# ╔═╡ 4a8d177e-fced-4c68-ab4d-3916f3ea3984
md"""
The maximum allowable computation time per Gram matrix element ``\tau_{ij}`` is estimated on the basis of the number ``N`` of calculations (one half of the square of the number of inputs, due to the symmetry of the matrix), the number ``P`` of parallel processes, and the overall time limit ``T``.

``\tau_{ij} = \frac{2PT}{N^2}``

Assuming the kernel computes in parallel on 48 processes with a time limit of 24 hours, for 5500 inputs, the maximum time per Gram matrix element is roughly 270 ms.
"""

# ╔═╡ 86452985-8025-414a-8663-672452fbd760
2 * 48 * 24 * 3600 / 5500^2

# ╔═╡ Cell order:
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╠═5699f8a5-11d6-453d-a867-8330134d080f
# ╟─5fec82c3-99fe-4ff0-aacd-7af622f07291
# ╠═1f2f45f6-57ba-4c29-845f-05685ceb299a
# ╠═f2c86d66-b801-4192-965c-b0b82a5c603a
# ╠═41c83665-9cff-43c1-912f-3d820d682e09
# ╠═b96dee5e-6c6b-4a1e-938a-5a9b42a96c3b
# ╠═39326496-e4dc-4b32-b538-feaa47066982
# ╠═a9c0e399-397e-43a7-a4ad-19af2c650d71
# ╠═87aa9631-ccef-4532-a06b-0aaee425d908
# ╟─53a7b87c-e676-4f42-81dc-dc38450078d1
# ╠═cf936eb7-020f-4e84-a284-8c7098d40cde
# ╠═a60742e4-cdbc-4cb2-a1fd-9fa03ce5cf64
# ╠═9875f635-1528-4b55-a3a7-fce2e6bbd866
# ╠═ff6ddb75-f9fd-4279-8c95-d7be0ce3e812
# ╠═0a97908e-bf4b-41a1-8cf6-e7905325393a
# ╠═f8fd1d2c-c1e5-4173-b704-8c345096b059
# ╠═f0bbbfe6-b799-4be7-8afb-61ef65f1a37e
# ╠═983b5c21-7618-483f-abe3-32efd6452130
# ╠═d57a3747-dcc6-46b7-946f-1828cae62fa2
# ╠═7663b5a0-d0a9-4a72-9d75-745a91160737
# ╠═575ccd4b-50f0-405d-9a39-48bc1265512e
# ╟─7dcc01ec-f087-467f-be59-e5404d44946f
# ╠═d55eaec7-0e62-44c5-a139-986b12e731fe
# ╠═cd3a297a-d05f-48d0-81a6-c37e1dfa3777
# ╟─024ed879-6da4-4770-81ca-9dc4cf9b4112
# ╟─6a323e50-ea97-41e1-8655-026ff5d73a00
# ╠═1806ea53-6738-4082-abd7-fc7e7ca5c463
# ╟─5d7efccd-9575-4145-bf7b-1c614539f733
# ╠═10184e95-d9a8-47b6-a2f3-8c25a29bca8a
# ╟─0990663a-0970-4468-9e42-15e77c209a79
# ╠═5c0a8f21-36b9-45f5-80b9-6726c69dd313
# ╟─00b1d071-8c5c-40c9-99a5-ae5b380e3151
# ╠═2b52cdea-7af8-4b16-97bd-5ad32538e7e3
# ╟─6d85fe2b-e9aa-4f91-af77-4f291399bc7d
# ╠═7e090139-72db-4658-b4f6-dac9bad7b848
# ╟─3e18dc29-6a47-4b35-8c38-fe374484f0a0
# ╟─7ee29c76-ecd8-4eef-a69d-506f0d171d36
# ╠═e5dc6300-3445-4f59-8131-1694a09dbbbe
# ╟─8193aea0-7cbc-4c1c-bd10-e6992e7eee02
# ╠═b85c081c-220c-4a20-83b5-2d1204d5a999
# ╠═d124ca26-8de5-4ac7-838e-dc37524af5c3
# ╠═765c3a87-9bce-488f-afb3-7db637c0f6c2
# ╟─d9afe6e0-86ce-49c4-a842-73862b369a99
# ╠═bbf82170-e0ee-49d9-8aad-22329b4f6f6f
# ╠═e30b4937-361f-47e3-95c1-3200a1466c12
# ╟─089eddca-c7b6-44f0-9133-681db448aa73
# ╠═465cc3ca-40ee-4b90-875c-a74c5e96f2d4
# ╠═a476ff95-5a4f-457d-b3e7-bc9b2f22cd07
# ╟─e9736f26-ce94-4854-aecb-b2991be0e348
# ╠═159e1fe2-eb50-42b8-b39f-c1d6ada69048
# ╟─00e01fc4-7c8f-44ab-acf7-0e1593edd386
# ╠═d5a8b311-d276-42b1-bdb9-5a44981bf601
# ╟─57d62a74-d13a-4333-8d90-cf0e69b0973c
# ╠═4611beba-ab89-4d68-997d-9d84c166c90a
# ╠═1ab6df2d-502e-4aa8-a008-a24e25aef2d1
# ╟─e6a6832f-4f76-4885-9258-b02c1d5a91e4
# ╟─19fd937a-20b4-4e82-9c08-7edcb4161e1f
# ╠═37e534f1-e03e-47c3-9a6f-ee9cc5f10de5
# ╠═9b21ce88-3ae4-45e1-b324-61ae3faf2ae4
# ╠═38e99f16-b499-4020-b51c-462b811bf135
# ╠═92db43c3-ae2c-4518-9c39-1403c887a315
# ╠═8b0473ee-47ab-4ad1-863e-edd2883407d0
# ╟─b41f3383-2105-4a92-8b7e-edfb00addae0
# ╟─f2f56957-1001-4bfa-962d-2210c4e8ce67
# ╟─4a8d177e-fced-4c68-ab4d-3916f3ea3984
# ╠═86452985-8025-414a-8663-672452fbd760
