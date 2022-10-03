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

# ╔═╡ 037a4025-ad55-4fb5-a14b-8b2da83db1c7
@revise using MolecularGraphKernels

# ╔═╡ e9f5391d-4832-440e-b61c-357daf275332
using Graphs, MetaGraphs, MolecularGraph

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

# ╔═╡ cda1019e-8970-4586-9c30-d9c5be453f58
md"""
# New Feature
"""

# ╔═╡ f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
md"""
### Development Code
"""

# ╔═╡ d02a6771-32aa-4f1d-97ad-7c3380f7b9da
begin
	import Graphs: is_directed, SimpleEdge
	import MetaGraphs: weighttype, PropDict, MetaDict
	
	struct Foo{T <: MolecularGraphKernels.AbstractProductGraph, U <: Real} <: AbstractMetaGraph{Int}
		graph::SimpleGraph{Int}
	    vprops::Dict{Int,PropDict}
	    eprops::Dict{SimpleEdge{Int},PropDict}
	    gprops::PropDict
	    weightfield::Symbol
	    defaultweight::U
	    metaindex::MetaDict
	    indices::Set{Symbol}
	end

	function Foo{S}(x, weightfield::Symbol, defaultweight::U) where {U <: Real, S <: MolecularGraphKernels.AbstractProductGraph}
	    T = eltype(x)
	    g = SimpleGraph(x)
	    vprops = Dict{T,PropDict}()
	    eprops = Dict{SimpleEdge{T},PropDict}()
	    gprops = PropDict()
	    metaindex = MetaDict()
	    idxs = Set{Symbol}()
	    return Foo{S}(g, vprops, eprops, gprops, weightfield, defaultweight, metaindex, idxs)
	end

	Foo{T}(x) where T = Foo{T}(x, :weight, 1.0)

	is_directed(::Foo) = false
	weighttype(foo::Foo) = Int

	Foo{T}(g1::MetaGraph, g2::MetaGraph) where T = Foo{T}(ProductGraph{T}(g1, g2).graph)

	function Foo{T}(g::MetaGraph) where T <: MolecularGraphKernels.AbstractProductGraph
	    newg = SimpleGraph{Int}(g.graph)
	    return MetaGraph(newg, g.defaultweight)
	end
end

# ╔═╡ 2b981075-65f5-45e7-b059-721aba3895eb
begin
	g1 = MetaGraph(smilestomol("c1ccccc1"))
	g2 = MetaGraph(smilestomol("c1ncccc1"))
end

# ╔═╡ 80e0a665-3ce1-4bf8-bc20-b5b68746c31e
Foo{Direct}(g1, g2)

# ╔═╡ 20c8a3e5-2d0e-4009-889d-541e52971f15
dpg = @btime ProductGraph{Direct}(g1, g2)

# ╔═╡ e52a708b-0b02-4a59-9881-2ad101ce0ad0
nv(dpg)

# ╔═╡ 2fe031a5-6295-44fe-8024-4dae6ce26438
ne(dpg)

# ╔═╡ 722ee3fa-f0f3-44c3-9011-f6b9484885d9
edges(dpg)

# ╔═╡ 37a62afe-1daa-4e61-855e-0463da06a94e
vertices(dpg)

# ╔═╡ 0263a5bf-7fbd-4fcd-b78e-c0e985051ebf


# ╔═╡ 13f3966f-4d67-4385-a37e-efffc4165280
md"""
### Tests
"""

# ╔═╡ 509a14f1-8e0b-4a88-9cbc-8010c4992aa0
# put code for test-driven development here

# ╔═╡ Cell order:
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─cda1019e-8970-4586-9c30-d9c5be453f58
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╠═d02a6771-32aa-4f1d-97ad-7c3380f7b9da
# ╠═2b981075-65f5-45e7-b059-721aba3895eb
# ╠═80e0a665-3ce1-4bf8-bc20-b5b68746c31e
# ╠═20c8a3e5-2d0e-4009-889d-541e52971f15
# ╠═e52a708b-0b02-4a59-9881-2ad101ce0ad0
# ╠═2fe031a5-6295-44fe-8024-4dae6ce26438
# ╠═722ee3fa-f0f3-44c3-9011-f6b9484885d9
# ╠═37a62afe-1daa-4e61-855e-0463da06a94e
# ╠═0263a5bf-7fbd-4fcd-b78e-c0e985051ebf
# ╟─13f3966f-4d67-4385-a37e-efffc4165280
# ╠═509a14f1-8e0b-4a88-9cbc-8010c4992aa0
