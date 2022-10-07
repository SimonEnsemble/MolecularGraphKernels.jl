### A Pluto.jl notebook ###
# v0.19.12

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

# ╔═╡ cda1019e-8970-4586-9c30-d9c5be453f58
md"""
# New Feature
"""

# ╔═╡ f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
md"""
### Development Code
"""

# ╔═╡ a23ae5a3-0c7f-41ad-815d-2a05e4407da4
# tinker here

# ╔═╡ 13f3966f-4d67-4385-a37e-efffc4165280
md"""
### Tests
"""

# ╔═╡ f0bb9a05-bfbe-4b65-b001-97e9563b34c5
# test here

# ╔═╡ 509c88ec-c05a-44ca-ba0c-e26f87c80f43
md"""
# Documentation and Examples
"""

# ╔═╡ 2cb9a9ce-150b-44ca-aeca-12c694d41f90
begin
	import Base.display
	display(mol::GraphMol) = HTML(drawsvg(mol, 250, 250))
end

# ╔═╡ 64b5c570-3f5a-4c58-96b3-29d54abcedaa
md"""
Molecule 1 SMILES: $(@bind mol₁_smiles confirm(TextField(default="O=C(OCCC(C)C)C")))
"""

# ╔═╡ 9ff77179-36ff-49d3-9447-99307660ca8f
mol₁ = smilestomol(mol₁_smiles)

# ╔═╡ e9f2c9b6-a17a-42e1-96cc-1935dc1cb48e
display(mol₁)

# ╔═╡ af4c6617-2ea5-47ba-b21a-12d38e0c1025
g₁ = MetaGraph(mol₁)

# ╔═╡ e19e9fc9-40db-439e-a947-26169aa222e7
viz_graph(g₁, layout_style=nothing)

# ╔═╡ a48556ba-cd3b-473e-b051-029a008556ee
md"""
Molecule 2 SMILES: $(@bind mol₂_smiles confirm(TextField(default="O=C1C[C@@H](C\\C=C1\\C)C(C)=C")))
"""

# ╔═╡ 03ef0f25-1381-4170-843f-b0b47ae3e744
mol₂ = smilestomol(mol₂_smiles)

# ╔═╡ 454a92f8-4c82-4462-b345-3104e7c2a55f
display(mol₂)

# ╔═╡ 9acb5371-3699-41cd-a70f-97d08d846462
g₂ = MetaGraph(mol₂)

# ╔═╡ 4ad2115d-4b90-433f-8943-e51386cddd00
viz_graph(g₂, layout_style=:spectral)

# ╔═╡ 6f19658a-c06e-47ee-8e42-a6358071c170
dpg = ProductGraph{Direct}(g₁, g₂)

# ╔═╡ 90f93435-0941-44bb-97c3-852c4a1736b2
viz_graph(dpg, layout_style=:circular)

# ╔═╡ b8637b31-3fe3-42ba-9bd0-8621c710f422
product_graph_adjacency_matrix{Direct}(g₁, g₂)

# ╔═╡ 4962c727-6f19-4c2d-9dc5-338ecc2914e3
random_walk_graph_kernel(dpg, 4)

# ╔═╡ 4c0b8a27-5457-4fac-bec3-43a7becd69a9
mpg = ProductGraph{Modular}(g₁, g₂)

# ╔═╡ a37e5286-52aa-477d-8ad8-b74090a58270
viz_graph(g₁)

# ╔═╡ 797a2fda-c3b0-4a4f-8c56-dcc159b263c6
viz_graph(g₂)

# ╔═╡ 0ad1de21-6abf-4eec-8db0-620647ace465
viz_graph(mpg, layout_style=nothing)

# ╔═╡ 6df530e6-a2fa-4fbb-bc0a-098a114593ec
product_graph_adjacency_matrix{Modular}(g₁, g₂)

# ╔═╡ a04ff14a-1fc9-452b-a640-a9807ecaafe6
length(maximal_cliques(SimpleGraph(mpg)))

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
# ╠═a23ae5a3-0c7f-41ad-815d-2a05e4407da4
# ╟─13f3966f-4d67-4385-a37e-efffc4165280
# ╠═f0bb9a05-bfbe-4b65-b001-97e9563b34c5
# ╟─509c88ec-c05a-44ca-ba0c-e26f87c80f43
# ╠═2cb9a9ce-150b-44ca-aeca-12c694d41f90
# ╟─64b5c570-3f5a-4c58-96b3-29d54abcedaa
# ╠═9ff77179-36ff-49d3-9447-99307660ca8f
# ╠═e9f2c9b6-a17a-42e1-96cc-1935dc1cb48e
# ╠═af4c6617-2ea5-47ba-b21a-12d38e0c1025
# ╠═e19e9fc9-40db-439e-a947-26169aa222e7
# ╟─a48556ba-cd3b-473e-b051-029a008556ee
# ╠═03ef0f25-1381-4170-843f-b0b47ae3e744
# ╠═454a92f8-4c82-4462-b345-3104e7c2a55f
# ╠═9acb5371-3699-41cd-a70f-97d08d846462
# ╠═4ad2115d-4b90-433f-8943-e51386cddd00
# ╠═6f19658a-c06e-47ee-8e42-a6358071c170
# ╠═90f93435-0941-44bb-97c3-852c4a1736b2
# ╠═b8637b31-3fe3-42ba-9bd0-8621c710f422
# ╠═4962c727-6f19-4c2d-9dc5-338ecc2914e3
# ╠═4c0b8a27-5457-4fac-bec3-43a7becd69a9
# ╠═a37e5286-52aa-477d-8ad8-b74090a58270
# ╠═797a2fda-c3b0-4a4f-8c56-dcc159b263c6
# ╠═0ad1de21-6abf-4eec-8db0-620647ace465
# ╠═6df530e6-a2fa-4fbb-bc0a-098a114593ec
# ╠═a04ff14a-1fc9-452b-a640-a9807ecaafe6
