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

# ╔═╡ d6b2f2de-b879-46fb-8379-adfb0fd4f73b
begin
	g₁ = MetaGraph(smilestomol("c1(c2)cscc1ccc2"))
	g₂ = MetaGraph(smilestomol("c1(o2)cscc1occ2"))
end;

# ╔═╡ f2c86d66-b801-4192-965c-b0b82a5c603a
viz_graph.([g₁, g₂])

# ╔═╡ 41c83665-9cff-43c1-912f-3d820d682e09
mpg = ProductGraph{Modular}(g₁, g₂)

# ╔═╡ b96dee5e-6c6b-4a1e-938a-5a9b42a96c3b
viz_graph(mpg; layout_style=:spectral)

# ╔═╡ 39326496-e4dc-4b32-b538-feaa47066982
max_cliques = maximal_cliques(mpg.graph)

# ╔═╡ 81c19258-a7de-4845-9665-3636d0256760
cliques = filter(c -> length(c) == maximum(length.(max_cliques)), max_cliques)

# ╔═╡ 7d2efc6b-1b24-49a1-8957-4c7f19036e53
begin
	c = cliques[1]
	d = Dict(j => i for (i, j) in enumerate(c))
	g = MetaGraph(length(c))
	for (i, v) in enumerate(c)
		set_prop!(g, i, :label, get_prop(mpg, v, :label))
		for n in neighbors(mpg, v)
			haskey(d, n) && add_edge!(g, v, d[n], Dict(:label => get_prop(mpg, v, n, :label)))
		end
	end
end

# ╔═╡ a9c0e399-397e-43a7-a4ad-19af2c650d71
viz_graph(g; layout_style=:circular)

# ╔═╡ Cell order:
# ╟─cd9f1c9c-ebcd-4733-a7ec-4fd743b0d81b
# ╟─9188ef1e-16fe-4a79-8ba6-2b0e907d743a
# ╠═6fe97eb4-c85d-4c2c-b892-e1c5ec91cc61
# ╟─fb64efc5-e959-401f-96d1-464de7d47547
# ╠═037a4025-ad55-4fb5-a14b-8b2da83db1c7
# ╟─971586d9-266b-4dfd-97d6-dc3aed449600
# ╠═e9f5391d-4832-440e-b61c-357daf275332
# ╟─f4f182e7-e8fe-4f1e-9867-0e01c8a850b1
# ╠═d6b2f2de-b879-46fb-8379-adfb0fd4f73b
# ╠═f2c86d66-b801-4192-965c-b0b82a5c603a
# ╠═41c83665-9cff-43c1-912f-3d820d682e09
# ╠═b96dee5e-6c6b-4a1e-938a-5a9b42a96c3b
# ╠═39326496-e4dc-4b32-b538-feaa47066982
# ╠═81c19258-a7de-4845-9665-3636d0256760
# ╠═7d2efc6b-1b24-49a1-8957-4c7f19036e53
# ╠═a9c0e399-397e-43a7-a4ad-19af2c650d71
