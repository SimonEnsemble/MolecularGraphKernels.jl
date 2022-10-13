### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# ╔═╡ 6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
begin
    import IOCapture, Pkg
    IOCapture.capture() do
        Pkg.activate(".")
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

# ╔═╡ 2cb9a9ce-150b-44ca-aeca-12c694d41f90
begin
    import Base.display
    display(mol::GraphMol) = HTML(drawsvg(mol, 250, 250))
end

# ╔═╡ a23ae5a3-0c7f-41ad-815d-2a05e4407da4
g₁ = MetaGraph(smilestomol("NC=O"))

# ╔═╡ e7ef4608-4ece-488e-9240-51460de06b50
g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))

# ╔═╡ 4892ce9c-bd45-41a2-a469-c881117fedd2
function csi_kernel1(g₁::MetaGraph, g₂::MetaGraph)::Int
	mpg = ProductGraph{Modular}(g₁, g₂)
	cliques = maximal_cliques(mpg.graph)
	return length(cliques)
end

# ╔═╡ bee34c9a-739f-45f3-a521-fdca88e3351b
k1 = @btime csi_kernel1(g₁, g₂)

# ╔═╡ 4e6d6ccc-50a2-4135-abae-eeb502a72056
md"""
## Algorithm
"""

# ╔═╡ dc345b6f-6b25-4640-85c3-7ef9ad444dae
md"""
Given an adjacency matrix $A$:
"""

# ╔═╡ 6df530e6-a2fa-4fbb-bc0a-098a114593ec
A = product_graph_adjacency_matrix(Modular, g₁, g₂)

# ╔═╡ 81eeb7b5-38d5-4b4f-8f8f-05bc1fd6ef49
md"""
a clique is represented by a square matrix consisting of a subset of rows and columns where all off-diagonal elements $A_{ij}$ are $1$:
"""

# ╔═╡ 7df7dfa4-e779-4167-8756-1c07c797523a
A[[1,2], [1,2]]

# ╔═╡ 65422c37-29cd-47e1-b98c-308af00e1bf2
A[[2,3,4], [2,3,4]]

# ╔═╡ db767f8f-77d5-43a7-bbd7-8ccced63c6d7
md"""
These subset matrices are node-order invariant:
"""

# ╔═╡ 5a9aad9a-16f7-4131-9a64-0fab741912db
A[[2,1],[2,1]]

# ╔═╡ 531bfd63-d0e7-46f0-8cf0-d9a2ab3f0ad5
A[[3,4,2], [3,4,2]]

# ╔═╡ e45d1404-3a17-488d-91e2-f285be59a9a4
md"""
A *maximal clique* is represented by such a matrix where the addition of any other indices from ``A`` does not generate a clique:
"""

# ╔═╡ 8009aed4-d85b-4802-ae64-3fee0619b3db
A[[1,2,3],[1,2,3]]

# ╔═╡ c9ed2288-f389-4505-ba9c-646dbf54d8f9
A[[1,2,4],[1,2,4]]

# ╔═╡ f2009f9e-da19-4cbd-a20d-7f8981afa917
A[[1,2,5],[1,2,5]]

# ╔═╡ 843a4d5d-ab6a-4a8b-ab4a-792f45e2272e
md"""
To find the maximal cliques of a graph, the maximal cliques containing each node may be determined by first identifying the subset of nodes which may possibly form such cliques, via inspection of a single column of the adjacency matrix:
"""

# ╔═╡ 90fa2e57-5dce-4873-ac9f-66d7cd93e9c6
findall(A[:,1])

# ╔═╡ b6bfa010-4bee-4880-9be9-4ec563e43f2a
md"""
So, the set ``{1,2}`` is a possible maximal clique, and because it has size ``2`` it must in fact be so.
"""

# ╔═╡ 507b71a4-4eb6-4903-96f7-058e4570f643
findall(A[:,2])

# ╔═╡ 2d6079ff-4ef6-48a2-a43e-3fcea2295eae
md"""
The index ``2`` may be a part of additional cliques.  We have already identified all maximal cliques for index ``1``, so we can exclude it from further searches.
"""

# ╔═╡ f44608e3-44fc-42a5-ab46-4136ac4fd336
begin
	local exhausted_indices = [1]
	local subset = setdiff(vcat(2, findall(A[:,2])), exhausted_indices)
	subset, A[subset, subset]
end

# ╔═╡ f05da4cc-aa3e-4d6c-a452-0ae48b6b4424
md"""
If index `3` is part of a clique with index `2`, then index `4` may also be in the same clique; but indices `5` and `6` may not.
"""

# ╔═╡ 0baf15b8-84a8-41ee-b572-f736946c328b
begin
	local exhausted_indices = [1]
	local subset = setdiff(vcat(2, findall(A[:,2])), exhausted_indices)
	local subset1 = [2, 3, 4]
	local subset2 = [2, 5, 6]
	subset1, A[subset1, subset1], subset2, A[subset2, subset2]
end

# ╔═╡ ac3d174a-4621-493b-89f5-a14b78c9677b
md"""
Each of the subsets ``{2,3,4}`` and ``{2,5,6}`` are (maximal) cliques, exhausting the search on index ``2``.
"""

# ╔═╡ 14dc85c0-d7ee-4aba-b7c4-03b21582e614
begin
	local exhausted_indices = [1, 2]
	local subset = setdiff(vcat(3, findall(A[:,3])), exhausted_indices)
	subset, A[subset, subset]
end

# ╔═╡ d32b9c09-1973-4d31-9a3e-18ce73589749
md"""
Ah, crap.  That's a subset of an already-discovered maximal clique.
"""

# ╔═╡ 8310e7ec-2f37-4b20-a9ee-33ef74075193
function csi_kernel2(g₁, g₂)::Int
	g = SimpleGraph(product_graph_adjacency_matrix(Modular, g₁, g₂))
	cliques = maximal_cliques(g)
	return length(cliques)
end

# ╔═╡ d5dca5dd-9e59-4f63-9928-a3d953a4b5c4
k2 = @btime csi_kernel2(g₁, g₂)

# ╔═╡ 13f3966f-4d67-4385-a37e-efffc4165280
@test k1 == k2

# ╔═╡ 853ce4ae-1d90-4eba-8101-53743f809179
@btime common_subgraph_isomorphism(g₁, g₂)

# ╔═╡ f448c7d5-9a5b-4098-b611-a88d2251a6f3
gram_matrix(common_subgraph_isomorphism, [g₁, g₂])

# ╔═╡ d6b2f2de-b879-46fb-8379-adfb0fd4f73b


# ╔═╡ Cell order:
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╠═2cb9a9ce-150b-44ca-aeca-12c694d41f90
# ╠═a23ae5a3-0c7f-41ad-815d-2a05e4407da4
# ╠═e7ef4608-4ece-488e-9240-51460de06b50
# ╠═4892ce9c-bd45-41a2-a469-c881117fedd2
# ╠═bee34c9a-739f-45f3-a521-fdca88e3351b
# ╟─4e6d6ccc-50a2-4135-abae-eeb502a72056
# ╟─dc345b6f-6b25-4640-85c3-7ef9ad444dae
# ╠═6df530e6-a2fa-4fbb-bc0a-098a114593ec
# ╟─81eeb7b5-38d5-4b4f-8f8f-05bc1fd6ef49
# ╠═7df7dfa4-e779-4167-8756-1c07c797523a
# ╠═65422c37-29cd-47e1-b98c-308af00e1bf2
# ╟─db767f8f-77d5-43a7-bbd7-8ccced63c6d7
# ╠═5a9aad9a-16f7-4131-9a64-0fab741912db
# ╠═531bfd63-d0e7-46f0-8cf0-d9a2ab3f0ad5
# ╟─e45d1404-3a17-488d-91e2-f285be59a9a4
# ╠═8009aed4-d85b-4802-ae64-3fee0619b3db
# ╠═c9ed2288-f389-4505-ba9c-646dbf54d8f9
# ╠═f2009f9e-da19-4cbd-a20d-7f8981afa917
# ╟─843a4d5d-ab6a-4a8b-ab4a-792f45e2272e
# ╠═90fa2e57-5dce-4873-ac9f-66d7cd93e9c6
# ╟─b6bfa010-4bee-4880-9be9-4ec563e43f2a
# ╠═507b71a4-4eb6-4903-96f7-058e4570f643
# ╟─2d6079ff-4ef6-48a2-a43e-3fcea2295eae
# ╠═f44608e3-44fc-42a5-ab46-4136ac4fd336
# ╟─f05da4cc-aa3e-4d6c-a452-0ae48b6b4424
# ╠═0baf15b8-84a8-41ee-b572-f736946c328b
# ╟─ac3d174a-4621-493b-89f5-a14b78c9677b
# ╠═14dc85c0-d7ee-4aba-b7c4-03b21582e614
# ╟─d32b9c09-1973-4d31-9a3e-18ce73589749
# ╠═8310e7ec-2f37-4b20-a9ee-33ef74075193
# ╠═d5dca5dd-9e59-4f63-9928-a3d953a4b5c4
# ╠═13f3966f-4d67-4385-a37e-efffc4165280
# ╠═853ce4ae-1d90-4eba-8101-53743f809179
# ╠═f448c7d5-9a5b-4098-b611-a88d2251a6f3
# ╠═d6b2f2de-b879-46fb-8379-adfb0fd4f73b
