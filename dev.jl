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
viz_graph(g₁; node_alpha_mask=g₁_node_alpha_mask, edge_alpha_mask=g₁_edge_alpha_mask, layout_style=:graphmol)

# ╔═╡ 575ccd4b-50f0-405d-9a39-48bc1265512e
viz_graph(g₂; node_alpha_mask=g₂_node_alpha_mask, edge_alpha_mask=g₂_edge_alpha_mask, layout_style=:graphmol)

# ╔═╡ b5067fb9-3543-40ea-bbad-768136438c18
length(imsgs)

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
# ╠═b5067fb9-3543-40ea-bbad-768136438c18
