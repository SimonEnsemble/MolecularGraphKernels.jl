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
    g₁;
    node_alpha_mask=g₁_node_alpha_mask,
    edge_alpha_mask=g₁_edge_alpha_mask,
    layout_style=:graphmol
),

viz_graph(
    g₂;
    node_alpha_mask=g₂_node_alpha_mask,
    edge_alpha_mask=g₂_edge_alpha_mask,
    layout_style=:graphmol
)

# ╔═╡ 7dcc01ec-f087-467f-be59-e5404d44946f
md"""
## Subgraph Matching Kernel 🚩
"""

# ╔═╡ 758e4ade-e697-4ca1-badc-03f90c0b9ec8
md"""
### w/o constraints
"""

# ╔═╡ 7954e762-afff-412a-8a72-d08f7b2c01dc
@btime common_subgraph_isomorphism(g₁, g₂)

# ╔═╡ 89974cca-2057-46bd-a7d2-e6adea92579e
@btime common_subgraph_isomorphism(g₁, g₂; λ=length)

# ╔═╡ a5d88fcb-ae3a-4119-97ae-26585d34a967
md"""
### constrained
"""

# ╔═╡ 09097b03-fff4-40f6-9af9-559799b66248
import MolecularGraphKernels.smkernel_c

# ╔═╡ 34228393-b352-4fa2-9bef-5bdb7ee317bd
function constrained_subgraph_matching(Gₚ::ProductGraph{T}, λ::Function)::Int where 
	T <: Union{Modular, Weighted}
	# Algorithm: SMKernel(w, C, P)
	# Input: Product graph Gₚ, weight function λ
	# Initial: value ← 0; SMKernel(1, ∅, Vₚ)
	# Param.: Weight w of the clique C, candidate set P
	# Output: Result of the kernel function value
	
	# initialize
	value = 0
	∅ = Int[]
	Vₚ = collect(vertices(Gₚ))

	# define recursive algorithm
	function smkernel(w::Int, C::Vector{Int}, P::Vector{Int})
		while length(P) > 0 # while |P| > 0 do
			v = first(P) # v ← arbitrary element of P
			C′ = union(C, v)
			w′ = w * smkernel_c(Gₚ, v) # multiply by vertex weight
			for u in C
				w′ *= smkernel_c(Gₚ, u, v)# multiply by edge weights
			end
			value += w′ * λ(C′)
			smkernel(w′, C′, intersect(P, [u for u in neighbors(Gₚ, v) if get_prop(Gₚ, u, v, :label) ≠ 0])) # extend c-clique 🚀
			P = setdiff(P, [v]) # P ← P \ {v}
		end
		return
	end

	# run algorithm
	smkernel(1, ∅, Vₚ)
	return value
end

# ╔═╡ bb4acf8d-89e3-47bb-ab40-4afa545e4c57
@btime constrained_subgraph_matching(
	ProductGraph{Modular}(g₁, g₂),
	_->1
)

# ╔═╡ b5e9f732-b966-4030-a6a5-a953a8ceedf2
@btime constrained_subgraph_matching(
	ProductGraph{Modular}(g₁, g₂),
	length
)

# ╔═╡ b41f3383-2105-4a92-8b7e-edfb00addae0
md"""
# Kernel Comp Time
"""

# ╔═╡ f2f56957-1001-4bfa-962d-2210c4e8ce67
md"""
### Calculation
"""

# ╔═╡ 4a8d177e-fced-4c68-ab4d-3916f3ea3984
md"""
The maximum allowable computation time per Gram matrix element ``\tau_{ij}`` is estimated on the basis of the number ``N`` of calculations (one half of the square of the number of inputs, due to the symmetry of the matrix), the number ``P`` of parallel processes, and the overall time limit ``T``.

``\tau_{ij} = \frac{2PT}{N^2}``

Assuming the kernel computes in parallel on 48 processes with a time limit of 24 hours, for 5500 inputs, the maximum time per Gram matrix element is under 300 ms.
"""

# ╔═╡ 86452985-8025-414a-8663-672452fbd760
round(Int, 1000 * 2 * 48 * 24 * 3600 / 5500^2)

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
# ╟─7dcc01ec-f087-467f-be59-e5404d44946f
# ╟─758e4ade-e697-4ca1-badc-03f90c0b9ec8
# ╠═7954e762-afff-412a-8a72-d08f7b2c01dc
# ╠═89974cca-2057-46bd-a7d2-e6adea92579e
# ╟─a5d88fcb-ae3a-4119-97ae-26585d34a967
# ╠═09097b03-fff4-40f6-9af9-559799b66248
# ╠═34228393-b352-4fa2-9bef-5bdb7ee317bd
# ╠═bb4acf8d-89e3-47bb-ab40-4afa545e4c57
# ╠═b5e9f732-b966-4030-a6a5-a953a8ceedf2
# ╟─b41f3383-2105-4a92-8b7e-edfb00addae0
# ╟─f2f56957-1001-4bfa-962d-2210c4e8ce67
# ╟─4a8d177e-fced-4c68-ab4d-3916f3ea3984
# ╠═86452985-8025-414a-8663-672452fbd760
