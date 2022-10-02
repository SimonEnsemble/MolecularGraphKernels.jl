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

# ╔═╡ c701affd-d895-4ec4-97d3-5d1788136afd
md"""
## Template for future dev notebooks
"""

# ╔═╡ 9188ef1e-16fe-4a79-8ba6-2b0e907d743a
md"""
### Conserved Preamble

Gives access to every tool I think is useful for developing local packages.

Expects to be run from the root of the package repo.
"""

# ╔═╡ fb64efc5-e959-401f-96d1-464de7d47547
md"""
### Variable macro call

Invokes file monitoring for the local package (by name).
"""

# ╔═╡ 971586d9-266b-4dfd-97d6-dc3aed449600
md"""
### Additional deps

External deps needed within the notebook for the specific project(s) at hand.
"""

# ╔═╡ cda1019e-8970-4586-9c30-d9c5be453f58
md"""
# Product Graph Types
"""

# ╔═╡ f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
md"""
### Development Code
"""

# ╔═╡ 2f9e0038-5365-407f-8725-eae6c35c24e0
begin
	a = MetaGraph(smilestomol("NC=O"))
	b = MetaGraph(smilestomol("C(NC=O)NC=O"))
end;

# ╔═╡ 3ff0f6ee-6f55-4442-b276-f7aece131b05
begin
	dpg = @btime ProductGraph{Direct}(a, b)
	dpg_adj_mat = @btime ProductGraphMatrix{Direct}(a, b)
	fpg = @btime ProductGraph{Factor}(a, b)
	fpg_adj_mat = @btime ProductGraphMatrix{Factor}(a, b)
end;

# ╔═╡ 13f3966f-4d67-4385-a37e-efffc4165280
md"""
### Tests
"""

# ╔═╡ a001dda5-584f-4d4d-baa4-95f9f1e3ef5a
@test is_isomorphic(dpg, ProductGraph{Direct}(fpg))

# ╔═╡ bbe96dac-fbd2-4d8d-a419-869cb4d48b26
@test is_isomorphic(ProductGraphMatrix(dpg), dpg_adj_mat)

# ╔═╡ 3dcaab69-898b-4704-83ed-18a3b863e498
@test is_isomorphic(dpg, dpg)

# ╔═╡ c5c86404-7ed5-46ef-9ed9-4ca566c80021
@test is_isomorphic(fpg, fpg)

# ╔═╡ 6fa712a6-332a-4874-a15c-c5a92df4deb8
@test is_isomorphic(dpg, ProductGraph{Direct}(fpg))

# ╔═╡ 52e50588-0bed-4727-b9f7-ec8edfa4ac47
md"""
# Vertex Pair Label Comparison
"""

# ╔═╡ 63287daf-9dcc-4a18-b416-1e11e543b042
md"""
### Development Code
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

# ╔═╡ 3b4da30c-735e-4bc1-8a95-188c1f41054f
md"""
### Tests
"""

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
# ╟─c701affd-d895-4ec4-97d3-5d1788136afd
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═cc80b503-8edd-4f65-bbd1-46910215c399
# ╟─cda1019e-8970-4586-9c30-d9c5be453f58
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╠═2f9e0038-5365-407f-8725-eae6c35c24e0
# ╠═3ff0f6ee-6f55-4442-b276-f7aece131b05
# ╟─13f3966f-4d67-4385-a37e-efffc4165280
# ╠═a001dda5-584f-4d4d-baa4-95f9f1e3ef5a
# ╠═bbe96dac-fbd2-4d8d-a419-869cb4d48b26
# ╠═3dcaab69-898b-4704-83ed-18a3b863e498
# ╠═c5c86404-7ed5-46ef-9ed9-4ca566c80021
# ╠═6fa712a6-332a-4874-a15c-c5a92df4deb8
# ╟─52e50588-0bed-4727-b9f7-ec8edfa4ac47
# ╟─63287daf-9dcc-4a18-b416-1e11e543b042
# ╠═4c8619ee-1e8a-4f56-a7a7-903a043daa86
# ╠═c0bb0d4b-c914-43e1-abd2-d87d7dd23721
# ╠═f28850c4-633d-460e-8708-bb6390694f53
# ╠═658a2625-59eb-4493-b28d-83f27455e3e7
# ╠═50a30de1-467a-4e2d-84e6-9a02eca806dd
# ╠═25912f04-5e13-4f0e-a4b0-aa5df8a5ed46
# ╠═c3d51178-60fb-4921-b749-5bd7a1a2176f
# ╠═d3a87146-84e4-4435-b3ee-47426a125180
# ╠═992a889e-99a9-456e-a2fe-5fcd7ef6b2ad
# ╟─3b4da30c-735e-4bc1-8a95-188c1f41054f
# ╠═0b016bc3-eca7-495a-a43c-c5779e731a81
# ╠═2d9a1890-74d0-4f16-9211-cea53fb64ad6
# ╠═42d6001b-d078-426a-887b-c381f72224ab
# ╠═7643b18e-bee2-4a6e-b20a-25a8efde345a
# ╠═236d188c-fc31-487d-9031-c7d0d5151755
