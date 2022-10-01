### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
begin
	import IOCapture, Pkg
	IOCapture.capture() do
		Pkg.activate(".")
		Pkg.add.([
			"Documenter"
			
			"Aqua"
			"BenchmarkTools"
			"ProfileCanvas"
			
			"PlutoLinks"
			"PlutoTeachingTools"
			"PlutoTest"
			"PlutoUI"
		])
	end
	using Documenter
	using Aqua, BenchmarkTools, ProfileCanvas
	using PlutoLinks, PlutoTeachingTools, PlutoTest, PlutoUI
	TableOfContents()
end

# ╔═╡ e9f5391d-4832-440e-b61c-357daf275332
using Graphs, MetaGraphs

# ╔═╡ 037a4025-ad55-4fb5-a14b-8b2da83db1c7
@revise using MolecularGraphKernels

# ╔═╡ cc80b503-8edd-4f65-bbd1-46910215c399
using MolecularGraph

# ╔═╡ cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
md"""
# Dev Block
"""

# ╔═╡ cda1019e-8970-4586-9c30-d9c5be453f58
md"""
# Product Graph Types
"""

# ╔═╡ 74ddf063-b484-42c5-98bc-8be864bdc5b1
foo(::Type{Factor}) = 1

# ╔═╡ 15e8d69f-3cc0-438e-93bf-116e8b2f8193
foo(::Direct) = 2

# ╔═╡ 10810df9-be0b-443d-bd68-91a8814218dc
foo(Factor)

# ╔═╡ d5327489-16f3-4533-88f6-a84eb2950faf
begin
	mol1 = MetaGraph(smilestomol("c1ccccc1"))
	mol2 = MetaGraph(smilestomol("c1cccnc1"))
end

# ╔═╡ bcc9e06f-d272-45bd-af94-082ce3ed2487
@btime ProductGraph{Direct}(mol1, mol2)

# ╔═╡ ae411318-6eae-4767-8668-abec1e1d90b3
@btime MolecularGraphKernels.product_graph(mol1, mol2, Direct)

# ╔═╡ 1afe17c9-9d85-4213-9d55-00017fabe23e
@btime ProductGraph{Factor}(mol1, mol2)

# ╔═╡ 52e50588-0bed-4727-b9f7-ec8edfa4ac47
md"""
# Vertex Pair Label Comparison
"""

# ╔═╡ 4c8619ee-1e8a-4f56-a7a7-903a043daa86
begin
	g1 = MetaGraph(3)
	g2 = MetaGraph(3)

	set_prop!(g1, 1, :label, 1)
	set_prop!(g1, 2, :label, 8)
	set_prop!(g1, 3, :label, 1)

	set_prop!(g2, 1, :label, 1)
	set_prop!(g2, 2, :label, 8)
	set_prop!(g2, 3, :label, 1)

	add_edge!(g1, 1, 2, Dict(:label => 1))
	add_edge!(g1, 2, 3, Dict(:label => 1))

	add_edge!(g2, 1, 2, Dict(:label => 1))
	add_edge!(g2, 2, 3, Dict(:label => 1))
end;

# ╔═╡ c0bb0d4b-c914-43e1-abd2-d87d7dd23721
begin
	h1 = MetaGraph(3)
	h2 = MetaGraph(3)

	set_prop!(h1, 1, :label, 8)
	set_prop!(h1, 2, :label, 1)
	set_prop!(h1, 3, :label, 1)

	set_prop!(h2, 1, :label, 1)
	set_prop!(h2, 2, :label, 1)
	set_prop!(h2, 3, :label, 8)

	add_edge!(h1, 1, 2, Dict(:label => 1))
	add_edge!(h1, 1, 3, Dict(:label => 1))

	add_edge!(h2, 1, 3, Dict(:label => 1))
	add_edge!(h2, 2, 3, Dict(:label => 1))
end;

# ╔═╡ f28850c4-633d-460e-8708-bb6390694f53
viz_graph(g1)

# ╔═╡ 658a2625-59eb-4493-b28d-83f27455e3e7
viz_graph(g2)

# ╔═╡ 50a30de1-467a-4e2d-84e6-9a02eca806dd
viz_graph(h1)

# ╔═╡ 25912f04-5e13-4f0e-a4b0-aa5df8a5ed46
viz_graph(h2)

# ╔═╡ c3d51178-60fb-4921-b749-5bd7a1a2176f
begin
	G = product_graph(g1, g2, :direct)
	H = product_graph(h1, h2, :direct)
end;

# ╔═╡ d3a87146-84e4-4435-b3ee-47426a125180
viz_graph(G)

# ╔═╡ 992a889e-99a9-456e-a2fe-5fcd7ef6b2ad
viz_graph(H)

# ╔═╡ 0b016bc3-eca7-495a-a43c-c5779e731a81
@test is_isomorphic(G, H)

# ╔═╡ 2d9a1890-74d0-4f16-9211-cea53fb64ad6
function display(g::MetaGraph)
	function get_props(i)
		prop_vec = ["$k:$v" for (k, v) in props(g, i)]
		return reduce(*, [" "] .* prop_vec)
	end
    println("---VERTICES---")
    for i in 1:nv(g)
        println("[$i]", get_props(i))
    end
    
    println("---EDGES---")
    for ed in edges(g)
        println("($(ed.src), $(ed.dst))", get_props(ed))
    end
end

# ╔═╡ 42d6001b-d078-426a-887b-c381f72224ab
display(G)

# ╔═╡ 7643b18e-bee2-4a6e-b20a-25a8efde345a
display(H)

# ╔═╡ 236d188c-fc31-487d-9031-c7d0d5151755
@test is_isomorphic(G, H; node_labels=[:label, :v₁v₂_pair])

# ╔═╡ Cell order:
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╠═cc80b503-8edd-4f65-bbd1-46910215c399
# ╟─cda1019e-8970-4586-9c30-d9c5be453f58
# ╠═74ddf063-b484-42c5-98bc-8be864bdc5b1
# ╠═15e8d69f-3cc0-438e-93bf-116e8b2f8193
# ╠═10810df9-be0b-443d-bd68-91a8814218dc
# ╠═d5327489-16f3-4533-88f6-a84eb2950faf
# ╠═bcc9e06f-d272-45bd-af94-082ce3ed2487
# ╠═ae411318-6eae-4767-8668-abec1e1d90b3
# ╠═1afe17c9-9d85-4213-9d55-00017fabe23e
# ╟─52e50588-0bed-4727-b9f7-ec8edfa4ac47
# ╠═4c8619ee-1e8a-4f56-a7a7-903a043daa86
# ╠═c0bb0d4b-c914-43e1-abd2-d87d7dd23721
# ╠═f28850c4-633d-460e-8708-bb6390694f53
# ╠═658a2625-59eb-4493-b28d-83f27455e3e7
# ╠═50a30de1-467a-4e2d-84e6-9a02eca806dd
# ╠═25912f04-5e13-4f0e-a4b0-aa5df8a5ed46
# ╠═c3d51178-60fb-4921-b749-5bd7a1a2176f
# ╠═d3a87146-84e4-4435-b3ee-47426a125180
# ╠═992a889e-99a9-456e-a2fe-5fcd7ef6b2ad
# ╠═0b016bc3-eca7-495a-a43c-c5779e731a81
# ╠═2d9a1890-74d0-4f16-9211-cea53fb64ad6
# ╠═42d6001b-d078-426a-887b-c381f72224ab
# ╠═7643b18e-bee2-4a6e-b20a-25a8efde345a
# ╠═236d188c-fc31-487d-9031-c7d0d5151755
