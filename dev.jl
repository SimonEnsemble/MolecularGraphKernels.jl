### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
begin
    import IOCapture, Pkg
    IOCapture.capture() do
        Pkg.activate(".")
        Pkg.resolve()
        Pkg.add.([
            "PlutoLinks"
            "PlutoUI"
        ])
        return
    end
    using PlutoLinks, PlutoUI
    TableOfContents()
end

# ╔═╡ 037a4025-ad55-4fb5-a14b-8b2da83db1c7
@revise using MolecularGraphKernels

# ╔═╡ e9f5391d-4832-440e-b61c-357daf275332
using Graphs, MetaGraphs, MolecularGraph, SparseArrays

# ╔═╡ cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b

# ╔═╡ 9188ef1e-16fe-4a79-8ba6-2b0e907d743a

# ╔═╡ fb64efc5-e959-401f-96d1-464de7d47547

# ╔═╡ 971586d9-266b-4dfd-97d6-dc3aed449600

# ╔═╡ f4f182e7-e8fe-4f1e-9867-0e01c8a850b1

# ╔═╡ 3b6fe89d-7727-4098-b958-e52baefe250d
md"""
## Display Function
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
    mol₁ = smilestomol("NC=O")
    mol₂ = smilestomol("CN(C=O)C=O")
    g₁ = MetaGraph(mol₁)
    g₂ = MetaGraph(mol₂)
    display.([mol₁, mol₂])
end

# ╔═╡ 53a7b87c-e676-4f42-81dc-dc38450078d1
md"""
## alpha masking
"""

# ╔═╡ 39326496-e4dc-4b32-b538-feaa47066982
imsgs = isomorphic_subgraphs(ProductGraph{Modular}(g₁, g₂))

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
    g₂;
    node_alpha_mask=g₂_node_alpha_mask,
    edge_alpha_mask=g₂_edge_alpha_mask,
    layout_style=:molecular
)

# ╔═╡ Cell order:
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╟─3b6fe89d-7727-4098-b958-e52baefe250d
# ╠═5699f8a5-11d6-453d-a867-8330134d080f
# ╟─5fec82c3-99fe-4ff0-aacd-7af622f07291
# ╠═1f2f45f6-57ba-4c29-845f-05685ceb299a
# ╟─53a7b87c-e676-4f42-81dc-dc38450078d1
# ╠═39326496-e4dc-4b32-b538-feaa47066982
# ╠═87aa9631-ccef-4532-a06b-0aaee425d908
# ╠═7663b5a0-d0a9-4a72-9d75-745a91160737
