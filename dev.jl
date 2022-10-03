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
# tinker with code for the new feature (or bug fix) here

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
# ╟─13f3966f-4d67-4385-a37e-efffc4165280
# ╠═509a14f1-8e0b-4a88-9cbc-8010c4992aa0
