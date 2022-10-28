### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

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

# ╔═╡ 41c83665-9cff-43c1-912f-3d820d682e09
mpg = ProductGraph{Modular}(g₁, g₂)

# ╔═╡ 53a7b87c-e676-4f42-81dc-dc38450078d1
md"""
## alpha masking
"""

# ╔═╡ 39326496-e4dc-4b32-b538-feaa47066982
imsgs = isomorphic_subgraphs(mpg)

# ╔═╡ 87aa9631-ccef-4532-a06b-0aaee425d908
begin
    local g = imsgs[1]
	local α₀ = 0.075
	
    dg = deepcopy(g)
    for e in edges(g)
        if get_prop(g, e, :label) == 0
            rem_edge!(dg, e)
        end
    end
	
    g₁_nodes = [get_prop(dg, v, :v₁v₂_pair)[1] for v in vertices(dg)]
	g₁_edges = [Graphs.SimpleEdge(g₁_nodes[src(e)], g₁_nodes[dst(e)]) for e in edges(dg)]
	g₁_node_alpha_mask = [v ∈ g₁_nodes ? 1 : α₀ for v in vertices(g₁)]
	g₁_edge_alpha_mask = [e ∈ g₁_edges || reverse(e) ∈ g₁_edges ? 1 : α₀ for e in edges(g₁)]

	g₂_nodes = [get_prop(dg, v, :v₁v₂_pair)[2] for v in vertices(dg)]
	g₂_edges = [Graphs.SimpleEdge(g₂_nodes[src(e)], g₂_nodes[dst(e)]) for e in edges(dg)]
	g₂_node_alpha_mask = [v ∈ g₂_nodes ? 1 : α₀ for v in vertices(g₂)]
	g₂_edge_alpha_mask = [e ∈ g₂_edges || reverse(e) ∈ g₂_edges ? 1 : α₀ for e in edges(g₂)]
end

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
## Correct CSI
"""

# ╔═╡ d55eaec7-0e62-44c5-a139-986b12e731fe
md"""
### *
"""

# ╔═╡ 34228393-b352-4fa2-9bef-5bdb7ee317bd


# ╔═╡ 8edd056e-f975-40b4-8201-277bdcf253de
@btime common_subgraph_isomorphism(g₁, g₂)

# ╔═╡ cd3a297a-d05f-48d0-81a6-c37e1dfa3777
@btime common_subgraph_isomorphism(g₁, g₂; λ=length)

# ╔═╡ 3e18dc29-6a47-4b35-8c38-fe374484f0a0
md"""
## Max-Clique Stuff
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
				C, Set(findall(reduce((x, y) -> x .&& y, [G[:, vⱼ] for vⱼ in C])))
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
@btime length(unique([get_prop(mpg, v, :v₁v₂_pair)[2] for v in len_dmcis(g₁, g₂)[2][1]]))

# ╔═╡ 92db43c3-ae2c-4518-9c39-1403c887a315
@test max_csi_len(g₁, g₂) == length(unique([get_prop(mpg, v, :v₁v₂_pair)[2] for v in len_dmcis(g₁, g₂)[2][1]]))

# ╔═╡ 8b0473ee-47ab-4ad1-863e-edd2883407d0
viz_graph(
	induced_subgraph(
		g₁, 
		unique(
			[get_prop(mpg, v, :v₁v₂_pair)[2] for v in len_dmcis(g₁, g₂)[2][1]]
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
# ╠═41c83665-9cff-43c1-912f-3d820d682e09
# ╠═53a7b87c-e676-4f42-81dc-dc38450078d1
# ╠═39326496-e4dc-4b32-b538-feaa47066982
# ╠═87aa9631-ccef-4532-a06b-0aaee425d908
# ╠═7663b5a0-d0a9-4a72-9d75-745a91160737
# ╠═575ccd4b-50f0-405d-9a39-48bc1265512e
# ╟─7dcc01ec-f087-467f-be59-e5404d44946f
# ╠═d55eaec7-0e62-44c5-a139-986b12e731fe
# ╠═34228393-b352-4fa2-9bef-5bdb7ee317bd
# ╠═8edd056e-f975-40b4-8201-277bdcf253de
# ╠═cd3a297a-d05f-48d0-81a6-c37e1dfa3777
# ╟─3e18dc29-6a47-4b35-8c38-fe374484f0a0
# ╠═465cc3ca-40ee-4b90-875c-a74c5e96f2d4
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
