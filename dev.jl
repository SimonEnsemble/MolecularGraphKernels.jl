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

# ╔═╡ feb99479-6e23-4838-8e22-f6307b7ba339
ProductGraph{Direct}(a, b)

# ╔═╡ eba23f46-559b-440a-8009-16eb69d1ecc5
ProductGraph{Factor}(a, b)

# ╔═╡ 42cd2333-3776-43af-a162-ea4b65660e11
ProductGraph{Direct}(ProductGraph{Factor}(a, b))

# ╔═╡ 6c69ed44-459e-4a08-aa0f-5d05e7cd5d3e
function isequal(a::ProductGraph, b::ProductGraph)::Bool
	if nv(a.graph) ≠ nv(b.graph)
		return false
	end
	for v in vertices(a.graph)
		if props(a.graph, v) ≠ props(b.graph, v)
			return false
		end
	end
	if ne(a.graph) ≠ ne(b.graph)
		return false
	end
	for e in edges(a.graph)
		if props(a.graph, e) ≠ props(b.graph, e)
			return false
		end
	end
	return true
end

# ╔═╡ 13f3966f-4d67-4385-a37e-efffc4165280
md"""
### Tests
"""

# ╔═╡ ff5f8f53-8f2f-496f-b258-79689971db27
l = 4

# ╔═╡ a001dda5-584f-4d4d-baa4-95f9f1e3ef5a
@test is_isomorphic(dpg, ProductGraph{Direct}(fpg))

# ╔═╡ bbe96dac-fbd2-4d8d-a419-869cb4d48b26
@test is_isomorphic(ProductGraphMatrix(dpg), dpg_adj_mat)

# ╔═╡ 4f73d6a9-4174-48d4-a4bb-cfe996b2fe31
@test graph_kernel(dpg.graph, l) == graph_kernel(dpg, l)

# ╔═╡ 3dcaab69-898b-4704-83ed-18a3b863e498
@test is_isomorphic(dpg, dpg)

# ╔═╡ 5853c817-1a6a-4056-afdc-c95f2013b612
@test graph_kernel(fpg, l) == graph_kernel(fpg_adj_mat, l)

# ╔═╡ 7746cfbf-ab46-42b8-9358-7720dae20da5
@test graph_kernel(dpg, l) == graph_kernel(dpg_adj_mat, l)

# ╔═╡ dcb43b58-0121-4e18-9ee8-392ac7265fdc
@test graph_kernel(fpg, l) == graph_kernel(fpg.graph, l)

# ╔═╡ c5c86404-7ed5-46ef-9ed9-4ca566c80021
@test is_isomorphic(fpg, fpg)

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

# ╔═╡ c3d51178-60fb-4921-b749-5bd7a1a2176f
begin
	G = ProductGraph{Direct}(g1, g2)
end;

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

# ╔═╡ 3990081f-259a-40d8-b35f-a110071dba90
display(g::ProductGraph) = display(g.graph)

# ╔═╡ af102c60-af61-4219-9442-95698a855a61
display(G)

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
# ╠═feb99479-6e23-4838-8e22-f6307b7ba339
# ╠═eba23f46-559b-440a-8009-16eb69d1ecc5
# ╠═42cd2333-3776-43af-a162-ea4b65660e11
# ╠═6c69ed44-459e-4a08-aa0f-5d05e7cd5d3e
# ╟─13f3966f-4d67-4385-a37e-efffc4165280
# ╠═ff5f8f53-8f2f-496f-b258-79689971db27
# ╠═a001dda5-584f-4d4d-baa4-95f9f1e3ef5a
# ╠═bbe96dac-fbd2-4d8d-a419-869cb4d48b26
# ╠═4f73d6a9-4174-48d4-a4bb-cfe996b2fe31
# ╠═3dcaab69-898b-4704-83ed-18a3b863e498
# ╠═5853c817-1a6a-4056-afdc-c95f2013b612
# ╠═7746cfbf-ab46-42b8-9358-7720dae20da5
# ╠═dcb43b58-0121-4e18-9ee8-392ac7265fdc
# ╠═c5c86404-7ed5-46ef-9ed9-4ca566c80021
# ╟─52e50588-0bed-4727-b9f7-ec8edfa4ac47
# ╟─63287daf-9dcc-4a18-b416-1e11e543b042
# ╠═4c8619ee-1e8a-4f56-a7a7-903a043daa86
# ╠═c3d51178-60fb-4921-b749-5bd7a1a2176f
# ╠═2d9a1890-74d0-4f16-9211-cea53fb64ad6
# ╠═3990081f-259a-40d8-b35f-a110071dba90
# ╠═af102c60-af61-4219-9442-95698a855a61
