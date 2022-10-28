### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
begin
	import IOCapture, Pkg
	IOCapture.capture() do
		Pkg.activate(".")
	end
	using BenchmarkTools, 
		Graphs, 
		MetaGraphs, 
		MolecularGraph, 
		MolecularGraphKernels, 
		PlutoUI
	TableOfContents(title="To-Do: 10/25")
end

# ╔═╡ 8e2aa3d6-9605-4525-9c9f-7cfcfa009c8f
md"""
# Bullet Point 1

	implement the common subgraph isomorphism kernel here¹ as common_subgraph_isomorphism_kernel(g1, g2). the weight function = # of nodes in the subgraph involved. can do in a few lines of code based on what you have.

¹[Grakel SM Kernel](https://ysig.github.io/GraKeL/0.1a7/kernels/subgraph_matching.html#)
"""

# ╔═╡ 62ded8f7-fcbb-4bc1-abae-c0b006ab9e8e
md"""
### Incomplete
"""

# ╔═╡ 2b15b3cb-d55f-4812-9ff2-8d74cd4bf8c0
md"""
!!! note
	What we have been doing with the CSI kernel is the "MCSI pseudo-kernel."
	`Graphs.jl` provides `maximal_cliques` but no routine for enumerating all cliques.
"""

# ╔═╡ 1d875c4b-2b05-46d7-8710-8918d95231c3
md"""
# Bullet Point 2

	to make sure your algo is correct, construct toy example as a test, compute the kernel in grakel in Python. should match. how does the speed of grakel compare with yours? compare output of random walk kernel too.
"""

# ╔═╡ 7485cf63-653e-4917-b0e4-fc4523e46239
md"""
### Incomplete
"""

# ╔═╡ f2c67185-eb8a-4a53-bc9f-ce639038cd9f
md"""
## MolecularGraphKernels.jl
"""

# ╔═╡ dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
begin
	g₁ = MetaGraph(smilestomol("NC=O"))
	g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))
end;

# ╔═╡ a4f7d0e4-eccc-481d-a914-17524d0c7da6
md"""
### Random Walk
"""

# ╔═╡ 66a808cf-e28e-4e47-aef0-c460988f788e
@btime random_walk(g₁, g₂; l=4)

# ╔═╡ 4ddd79f5-8607-47b7-9ee6-264fa69a6241
md"""
### Maximal Common Subgraph Isomorphism
"""

# ╔═╡ d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
@btime common_subgraph_isomorphism(g₁, g₂)

# ╔═╡ a5494cac-b8b8-4985-9606-8d3b6999f2f1
@btime common_subgraph_isomorphism(g₁, g₂; λ=length)

# ╔═╡ 4a37f20f-73a5-4d34-af97-d3badbb22c2e
md"""
## Grakel
"""

# ╔═╡ d3add5a8-9599-4196-984c-b597ed1d556a
md"""
!!! danger "Inconvenience"
	Grakel won't run on Python 3.10, so can't get it to run in notebook.
"""

# ╔═╡ e07f33f1-25b2-4a1f-908f-b9942b1f652c
println.(readlines("using_grakel.py"));

# ╔═╡ f3112b3e-fd57-4cf9-98c1-32ffbe084740
md"""
!!! danger "Annoyance"
	Making graphs for Grakel is excruciatingly tedious.
"""

# ╔═╡ 31d458d7-1c1a-45d3-abe7-d39da4c5419f
md"""
!!! danger "Problem"
	While the Grakel docs state that edge labels can be applied by index, `SubgraphMatching` disagrees.

	Each edge must be encoded as two pairs of vertices, and edge labels are non-optional!
"""

# ╔═╡ 00697c83-ccdf-4746-9731-665b0e0fdc06
md"""
!!! danger "Confusion"
	Grakel won't compute kernels between graphs without "fitting" to something first.
	😕?
"""

# ╔═╡ 4417ccab-5da4-4163-a36e-090e1fb2d9f9
run(`python3 using_grakel.py`);

# ╔═╡ 2ce8da75-2cdc-4d4c-80ee-8a2ac055e74c
md"""
!!! danger "Problem"
	Results from `RandomWalk` are astronomically huge. What is this calculating?
"""

# ╔═╡ 163aa4cc-0bd1-499f-8593-ea2ddd8ac02e
md"""
!!! danger "Problem"
	Results from `RandomWalkLabeled` are on the same order as, but still not equal to ours. Why?
"""

# ╔═╡ 99704c71-4fda-41c0-9225-13c687b7636e
md"""
!!! danger "Problem"
	Grakel seems to ignore node labels in `RandomWalk` and `RandomWalkLabeled`.

	Changing node labels has no effect on results!
"""

# ╔═╡ 9827d564-7965-493d-ba41-28d62ac0be90
md"""
!!! warning "Speed Comparison"
	Can't compare our kernels directly, b/c don't compute same things...

	*But* what we do calculate takes between 1 and 2 orders of magnitude less time than their code, so that's something..?
"""

# ╔═╡ 79994dc7-6138-4e3b-b50e-1b2769f80eb2
md"""
# Bullet Point 3

	employ CSI kernel to create Gram matrix for cannabanoids. use diffusion map with custom kernel matrix (see here² for package) to embed into 2D space. create function that will color the points according to a given label (what proteins they bind to) to see if the clustering can pick up on this. (even if nothing comes of this, great way to quickly get your workflow running/figure out how diffusion maps work/ what they are doing/benchmark your Gram matrix. how much time does it take to compute the Gram matrix?)

²[ManifoldLearning.jl diffmap](https://wildart.github.io/ManifoldLearning.jl/dev/diffmap/#StatsAPI.fit-Union%7BTuple%7BT%7D,%20Tuple%7BType%7BDiffMap%7D,AbstractArray%7BT,2%7D%7D%7D%20where%20T%3C:Real)
"""

# ╔═╡ c9c8f3df-205a-4b29-adfc-fc4ef1143da1
md"""
### Not Started
"""

# ╔═╡ 7adc7b44-ee3e-41c5-9977-0a64a3fb88bb
md"""
# Bullet Point 4

	search the literature to answer:
		
	1. is finding the largest common subgraph still a proper kernel? i.e positive semi-definite? not clear to me that it is. look in reviews of graph kernels, to see if they discuss anything like this. or on grakel page? [answer = no, below. so we shouldn't use that.]
		
	2. what prediction tasks have used the common subgraph isomorphism kernel? what are the size of the graphs/number of graphs? keep these papers handy so that you can cite them in your paper. good way to find out = look at the citations to the subgraph matching kernel paper. but also Google Scholar search for common subgraph isomorphism kernel (been around longer).
"""

# ╔═╡ 59742d77-b5cc-4379-bd02-87780ce61c6e
md"""
### Incomplete
"""

# ╔═╡ 1c61675b-01c1-4040-861b-e797c62c4bb7
md"""
# Bullet Point 5

	Sec. 3.3.1 of subgraph matching kernel paper here, "Restriction to Subgraph Classes", paragraph 2. the algorithm can be made faster and perhaps make for a more reasonable kernel (similarity description) by only looking at cliques spanned by c-edges. this is done on-the-fly when enumerating the cliques, not after the cliques are detected with post-filtering. "only enumerate c-cliques by making sure that only vertices are added that are adjacent to a vertex in the current clique via at least one c-edge". so when u try to extend a current clique into a maximal clique, you check if the new potential vertex is connected via a c-edge to one of the vertices in the current (possibly, non-maximum) clique. so Xiaoli was right. this is a simple modification of the clique enumeration algorithm and worth implementing. instead of inventing your own clique detection algorithm, I think you should use an algorithm that is already out there, e.g the Bron–Kerbosch algorithm.³

³[Bron-Kerbosch](https://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm)
"""

# ╔═╡ c9190f99-73f7-4bff-9138-a1fb69a9867c
md"""
### Incomplete
"""

# ╔═╡ 614434f8-8ddc-49ae-adaf-bf69aa939944
md"""
# Bullet Point 6

	what does Graphs.jl use for clique detection? is it Bron-Kerbosch? you can look at the maximal_cliques algorithm to see. the above implementation could be, copying their code, then only extending the clique if the new proposed vertex is connected to current clique by a c-edge. or at least taking inspiration from their code already to enumerate your own cliques.
"""

# ╔═╡ 52db250d-1403-4d93-a7e4-bcd9172c35d6
md"""
### Incomplete
"""

# ╔═╡ 8d805709-d844-48c3-b72c-f5b45a009a09
md"""
!!! note
	`Graphs.jl` doesn't do the kind of clique detection we want, but maybe their algorithm can be tweaked for it.  Since it is specifically tailored for maximal clique enumeration, probably better to go ground-up with it.
"""

# ╔═╡ 9ad64087-4671-4e15-a865-5088334fff6e
md"""
# Bullet Point 7

	BTW it is faster to not store the cliques, but once a max clique is found, increment the kernel. this is what Algorithm 1 is doing in the subgraph matching kernel paper, I think.
"""

# ╔═╡ 107bc909-6f05-428e-9923-992faae0b8d5
md"""
### Incomplete
"""

# ╔═╡ d70fcbca-1e41-4ecf-8a92-d9ea337bae78
md"""
!!! note
	The idea of storing the cliques was to use their info after, e.g. for eliminating automorphic duplications.  Agree that for any algorithm that can compute the contribution to the kernel without storing the clique, should do so.
"""

# ╔═╡ 1c9ff108-d3e4-4253-8a06-abde6e25207d
md"""
# Bullet Point 8

	after above explanation, use your toy example to enumerate the connected graphs that are being counted. should be much fewer now and more reasonable/explainable, I expect.
"""

# ╔═╡ 53e3ad8b-703a-44a4-8fdc-b92b26265caf
md"""
### Incomplete
"""

# ╔═╡ 3c7311e6-3cab-4aaf-8fbf-3a4a4580675a
md"""
!!! note
	Requiring that an MPG clique be spanned by c-edges vastly reduces the search space, but the overall complexity is still pretty bad.  The MPG of two graphs with size ``n`` has a maximum size of ``n!``, and clique enumeration on a graph of ``N`` nodes has worst-case complexity ``\mathcal{O}(N^k)`` given maximum clique size ``k``.  That means clique enumeration on the MPG is at worst ``\mathcal{O}(n!^n)`` (in the case of comparing a complete graph with itself).  Constraining cliques by growing only along c-edge subgraphs may make this better, but almost anything would be...
"""

# ╔═╡ 616d284c-431f-4d4e-a32c-08844ee6520e
md"""
# Bullet Point 9

	The 2D graph kernel is based on the concept of the maximum common subgraph (MCS) of two graphs. Since determining a MCS is NP-complete, exact algorithms that guarantee to return an optimal solution are too time-consuming for large data sets of graphs of even moderate size. We therefore resort to approximate solutions using the graduate assignment algorithm. (33) Since the proposed graph kernel is not positive semidefinite...⁴

⁴[Mohr, et al.](https://pubs.acs.org/doi/full/10.1021/ci900367j)
"""

# ╔═╡ 85f17dfb-40a8-46e0-8113-7f00307be7a6
md"""
!!! note
	SVM needs a PSD kernel... but P-SVM doesn't (Mohr refs 29-31!

	Something to keep in mind...
"""

# ╔═╡ ecd3bce7-7ec3-431c-ba10-bf3869ef85f9
md"""
# Bullet Point 10

	(regarding kernel score normalization⁵) does it remain positive semi-definite? I guess that should be an option in the Gram matrix function, since then you can do that after u compute the whole matrix (with its diagonal).

⁵[`normalize` kwarg](https://ysig.github.io/GraKeL/0.1a8/documentation/core_concepts.html#what-is-the-grakel-graphkernel-class)
"""

# ╔═╡ f9347e33-e568-4308-a086-560836906369
md"""
### Not Started
"""

# ╔═╡ Cell order:
# ╠═9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
# ╟─8e2aa3d6-9605-4525-9c9f-7cfcfa009c8f
# ╟─62ded8f7-fcbb-4bc1-abae-c0b006ab9e8e
# ╠═2b15b3cb-d55f-4812-9ff2-8d74cd4bf8c0
# ╟─1d875c4b-2b05-46d7-8710-8918d95231c3
# ╟─7485cf63-653e-4917-b0e4-fc4523e46239
# ╟─f2c67185-eb8a-4a53-bc9f-ce639038cd9f
# ╠═dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
# ╟─a4f7d0e4-eccc-481d-a914-17524d0c7da6
# ╠═66a808cf-e28e-4e47-aef0-c460988f788e
# ╟─4ddd79f5-8607-47b7-9ee6-264fa69a6241
# ╠═d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
# ╠═a5494cac-b8b8-4985-9606-8d3b6999f2f1
# ╟─4a37f20f-73a5-4d34-af97-d3badbb22c2e
# ╟─d3add5a8-9599-4196-984c-b597ed1d556a
# ╠═e07f33f1-25b2-4a1f-908f-b9942b1f652c
# ╟─f3112b3e-fd57-4cf9-98c1-32ffbe084740
# ╟─31d458d7-1c1a-45d3-abe7-d39da4c5419f
# ╟─00697c83-ccdf-4746-9731-665b0e0fdc06
# ╠═4417ccab-5da4-4163-a36e-090e1fb2d9f9
# ╟─2ce8da75-2cdc-4d4c-80ee-8a2ac055e74c
# ╟─163aa4cc-0bd1-499f-8593-ea2ddd8ac02e
# ╟─99704c71-4fda-41c0-9225-13c687b7636e
# ╟─9827d564-7965-493d-ba41-28d62ac0be90
# ╟─79994dc7-6138-4e3b-b50e-1b2769f80eb2
# ╠═c9c8f3df-205a-4b29-adfc-fc4ef1143da1
# ╟─7adc7b44-ee3e-41c5-9977-0a64a3fb88bb
# ╠═59742d77-b5cc-4379-bd02-87780ce61c6e
# ╟─1c61675b-01c1-4040-861b-e797c62c4bb7
# ╠═c9190f99-73f7-4bff-9138-a1fb69a9867c
# ╟─614434f8-8ddc-49ae-adaf-bf69aa939944
# ╠═52db250d-1403-4d93-a7e4-bcd9172c35d6
# ╟─8d805709-d844-48c3-b72c-f5b45a009a09
# ╠═9ad64087-4671-4e15-a865-5088334fff6e
# ╠═107bc909-6f05-428e-9923-992faae0b8d5
# ╟─d70fcbca-1e41-4ecf-8a92-d9ea337bae78
# ╟─1c9ff108-d3e4-4253-8a06-abde6e25207d
# ╠═53e3ad8b-703a-44a4-8fdc-b92b26265caf
# ╟─3c7311e6-3cab-4aaf-8fbf-3a4a4580675a
# ╟─616d284c-431f-4d4e-a32c-08844ee6520e
# ╟─85f17dfb-40a8-46e0-8113-7f00307be7a6
# ╟─ecd3bce7-7ec3-431c-ba10-bf3869ef85f9
# ╠═f9347e33-e568-4308-a086-560836906369
