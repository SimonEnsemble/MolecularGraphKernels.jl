### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
begin
    import IOCapture, Pkg
    IOCapture.capture() do
        return Pkg.activate(".")
    end
    using BenchmarkTools, PlutoUI
	using Graphs, MetaGraphs, MolecularGraph
	using MolecularGraphKernels
    TableOfContents(; title="To-Do: 10/25")
end

# ╔═╡ 1d875c4b-2b05-46d7-8710-8918d95231c3
md"""
# Algorithm Correctness 1 🚧

	to make sure your algo is correct, construct toy example as a test, compute the kernel in grakel in Python. should match. how does the speed of grakel compare with yours? compare output of random walk kernel too.
"""

# ╔═╡ dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
begin
    g₁ = MetaGraph(smilestomol("NC=O"))
    g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))
end;

# ╔═╡ e07f33f1-25b2-4a1f-908f-b9942b1f652c
println.(readlines("using_grakel.py"));

# ╔═╡ 4417ccab-5da4-4163-a36e-090e1fb2d9f9
run(`python3 using_grakel.py g1 g1`);

# ╔═╡ 25e5176a-5716-4cd3-bba3-4f27af205421
@btime random_walk(g₁, g₁; l=4)

# ╔═╡ f9e900aa-3143-4f34-a559-ad7b40a5909a
@btime common_subgraph_isomorphism(g₁, g₁)

# ╔═╡ 185a82b7-6b4d-43e4-b428-88c2b0876399
@btime common_subgraph_isomorphism(g₁, g₁; c_cliques=true)

# ╔═╡ a78a1080-e1e6-4912-ac3b-4911f9e15564
@btime common_subgraph_isomorphism(g₁, g₁; λ=length)

# ╔═╡ 4fc6a508-4702-4f8a-a623-b3f5e002c8dc
run(`python3 using_grakel.py g1 g2`);

# ╔═╡ 66a808cf-e28e-4e47-aef0-c460988f788e
@btime random_walk(g₁, g₂; l=4)

# ╔═╡ d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
@btime common_subgraph_isomorphism(g₁, g₂) # default: weight = 1 for all cliques

# ╔═╡ a5494cac-b8b8-4985-9606-8d3b6999f2f1
@btime common_subgraph_isomorphism(g₁, g₂; λ=length) # weight = clique size

# ╔═╡ 2ce8da75-2cdc-4d4c-80ee-8a2ac055e74c
md"""
!!! danger "Grakel Problems"
	- Results from `RandomWalk` are astronomically huge. 
	- Results from `RandomWalkLabeled` are on the same order as, but still not equal to ours.
	- Grakel seems to ignore node labels in `RandomWalk` and `RandomWalkLabeled`.
"""

# ╔═╡ 79994dc7-6138-4e3b-b50e-1b2769f80eb2
md"""
# Cannabinoid Clustering 🚧

	1. employ CSI kernel to create Gram matrix for cannabanoids. 

	2. use diffusion map w/ kernel matrix (see here²) to embed into 2D space. 

	3. color points according to protein, see if the clustering can pick up on this.

	4. how much time does it take to compute the Gram matrix?)

²[ManifoldLearning.jl diffmap](https://wildart.github.io/ManifoldLearning.jl/dev/diffmap/#StatsAPI.fit-Union%7BTuple%7BT%7D,%20Tuple%7BType%7BDiffMap%7D,AbstractArray%7BT,2%7D%7D%7D%20where%20T%3C:Real)
"""

# ╔═╡ ce79c1aa-1d19-4751-9199-35535e094c67
md"""
!!! note
	Used CSI kernel w/ c-clique constraint (both `λ(C) = 1` and `λ(C) = length(C)`).
	This takes ca. 10 min on a single core @ 4.3 GHz.
"""

# ╔═╡ b4ad3fa6-660b-44ec-831f-b27a794a9ed1
md"""
!!! warning "In Progress"
	See diffmap.jl in private repo
"""

# ╔═╡ 7adc7b44-ee3e-41c5-9977-0a64a3fb88bb
md"""
# Lit Search 🚧

	search the literature to answer:
		
	1. is finding the largest common subgraph still a proper kernel? i.e positive semi-definite? not clear to me that it is. look in reviews of graph kernels, to see if they discuss anything like this. or on grakel page? [answer = no, below. so we shouldn't use that.]
		
	2. what prediction tasks have used the common subgraph isomorphism kernel? what are the size of the graphs/number of graphs? keep these papers handy so that you can cite them in your paper. good way to find out = look at the citations to the subgraph matching kernel paper. but also Google Scholar search for common subgraph isomorphism kernel (been around longer).
"""

# ╔═╡ 10923005-d271-451f-8136-2a8b2cc3b9cc
md"""
!!! note
	Largest common subgraph is not a kernel
"""

# ╔═╡ 527c224d-fa2d-4b01-b17e-501ab62f6167
md"""
### Tasks (CSI Kernels)

**Nothing!**

Closest is `CSI_GED`, which gives no application, only runtime comparison on ``k(G_1,G_1)`` for 100k molecules of average atom count 24, and 10k molecules of average atom count 45.  The time limit set was ca. 2 hours; that's 1.5 seconds per computation.  We can afford ca. 300 ms per computation.
"""

# ╔═╡ e3eb5365-4c29-4eda-b0b4-7b92d17430e8
md"""
### Tasks (MCS Algorithms)

- small molecule classification
- compound activity prediction / QSAR
- reaction mapping
- database searching
- small molecule-μRNA binding
- metabolite prediction
"""

# ╔═╡ 9eb41dc9-dfb0-489a-89d1-e9dc47f0327b
md"""
### Tasks (Graph Kernels, Generally)

- atomization energy
- pure-substance phase diagrams
- reaction yield
- RNA structure
- protein-protein docking
- solvation energy
- gene function
- metal surface adsorption energy
- QSAR
- pKa
- biomolecule receptor agonism
"""

# ╔═╡ 036bed20-ad56-4913-b9d2-473d4f6773a2
md"""
### Cool Stuff

- quantum computing
- stereochemistry
"""

# ╔═╡ b346bcd0-9039-4175-a70d-59a2856ccb23
md"""
# more stuff 🚩
"""

# ╔═╡ aeae9d17-cb45-4619-870c-a2b8ce23ee3f
md"""
Function to convert Julia graphs into python graphs for grakel automatically. Will need to do scaling studies and compare performance if you wish to publish the package eventually, anyway, and will need to test on more complicated example, as the simple graphs may not capture facets/areas where bugs could be.
	
	
	
what is the performance comparison on a realistically sized graph? I think scaling vs # of nodes is important in benchmarking in general but esp here b/c the code is actually in C for grakel, so don't want to be looking at Python to C translation necessarily. 
	
	
	
 

	
It seems the labels are not being taken to consideration in grakel for the random walk kernel. Why? is something wrong with the input? is it flat-out wrong? I'd post an issue if so.
	
	
	
 

	
Your implementation of CSI is not counting only connected subgraphs. (see photo below for example we discussed)
	
	
	
 

	
What is the original paper for the CSI kernel?
	
	
	
 

	
SM kernel is for soft matches, which your code does not support. Change name of functions to CSI kernel to reflect this? in my view, SM kernel is pretty distinctive in that it implies soft matching via non-Dirac node/edge kernel. I think this is why they gave their kernel a new name, instead of something like 'extended CSI kernel'. kind of like calling a linear regression model a neural network. technically true, but misleading. 
"""

# ╔═╡ d15a77f2-6d34-498e-8e4e-8d9f92a477c5
md"""
Algorithm ``SMK`` in the paper for calculating ``k_{CSI}`` with ``\lambda(C)=1``:

1. ``\text{while } |P| > 0 \text{ do }``

2. ``\text{ }\text{ }\text{ }\text{ }v\leftarrow\text{arbitrary element of }P``

3. ``\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C\cup {v}``

4. ``\text{ }\text{ }\text{ }\text{ }value\leftarrow value+1``

5. ``\text{ }\text{ }\text{ }\text{ }P^\prime=P\cap N(v)``

6. ``\text{ }\text{ }\text{ }\text{ }SMK(C^\prime,P^\prime)``

7. ``\text{ }\text{ }\text{ }\text{ }P\leftarrow P \setminus {v}``
"""

# ╔═╡ 742ff2e7-84ba-48a0-8ffb-6c72273ee427
md"""
We thought we could change line 5 to be:

``P^\prime= \{u\in P\cap N(v):\exists k\in C^\prime\rightarrow l(u,k)\ne d\}``

but that leads to under-counting by eliminating too many candidate nodes.
"""

# ╔═╡ bcf5f6aa-96b6-4b6c-bc14-e99bd1647280
md"""
I also tried, among other ideas,

``P^\prime= \{u\in P\cap N(v):\exists k\in C^\prime\cup P\rightarrow l(u,k)\ne d\}`` (overcounts)

Then it occurred to me that when Kriege wrote:

	only enumerate c-cliques by making sure that only vertices are added that are adjacent to a vertex in the current clique via at least one cedge

he *didn't* mean:

	only consider as candidate nodes those which extend the current clique while maintaining *c*-edge spanning

but rather:

	only add a node from the candidate nodes to the growing clique if it maintains  *c*-edge spanning

which means that the real focus of the change is at line 3!
"""

# ╔═╡ 3a3ad8f4-fd72-4b0a-9ba3-400c5136eec5
md"""
Algorithm ``SMK`` as I have it now for calculating ``k_{CCSI}`` with ``\lambda(C)=1``:

1. ``\text{while } |P| > 0 \text{ do }``

2. ``\text{ }\text{ }\text{ }\text{ }v\leftarrow\text{arbitrary element of }P``

3. ``\text{ }\text{ }\text{ }\text{ }\text{if }C\cup {v}\in\mathcal{C}(G_p)\text{ do}``

4. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C\cup {v}``

5. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }value\leftarrow value+1``

6. ``\text{ }\text{ }\text{ }\text{ }\text{else}``

7. ``\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }\text{ }C^\prime\leftarrow C``

6. ``\text{ }\text{ }\text{ }\text{ }SMK(C^\prime,P\cap N(v))``

7. ``\text{ }\text{ }\text{ }\text{ }P\leftarrow P \setminus {v}``
"""

# ╔═╡ 0f0fbb49-a724-4804-ba8a-cc31734c17dd
function extends_clique(Gₚ, C, v)
	if C == []
		return true
	end
	for u in C
		if has_edge(Gₚ, u, v) && get_prop(Gₚ, u, v, :label) ≠ 0
			return true
		end
	end
	return false
end

# ╔═╡ a76c2d1f-faf2-453a-ab90-2b2bb72ecc13
function test_algo(g₁, g₂)
	value = 0
	Gₚ = ProductGraph{Modular}(g₁, g₂)
	Vₚ = collect(vertices(Gₚ))
	cliques = []
	
	function kernel(C, P)
		while length(P) > 0
			v = first(P)
			if extends_clique(Gₚ, C, v)
				C′ = union(C, v)
				push!(cliques, C′)
				value += 1
			else
				C′ = C
			end
			kernel(C′, intersect(P, neighbors(Gₚ, v)))
			P = setdiff(P, [v])
		end
	end

	kernel([], Vₚ)
	return value, cliques
end

# ╔═╡ b1373bbd-9b4b-40ee-b0d2-cec3d5c0dcdf
@btime test_algo(g₁, g₁)

# ╔═╡ 395b37e7-b0ea-4c20-8a9a-33dd2cf195de
@btime test_algo(g₁, g₂)

# ╔═╡ 0e9d12f1-92f4-4b88-b593-d6afcca8708a
begin
	mpg = ProductGraph{Modular}(g₁, g₂)
	viz_graph(MetaGraph(mpg))
end

# ╔═╡ Cell order:
# ╠═9aa20e3a-559a-11ed-1cd0-358ee55a4bb0
# ╟─1d875c4b-2b05-46d7-8710-8918d95231c3
# ╠═dbf109fc-30cf-41a7-b5f7-f98dc1434d8f
# ╟─e07f33f1-25b2-4a1f-908f-b9942b1f652c
# ╠═4417ccab-5da4-4163-a36e-090e1fb2d9f9
# ╠═25e5176a-5716-4cd3-bba3-4f27af205421
# ╠═f9e900aa-3143-4f34-a559-ad7b40a5909a
# ╠═185a82b7-6b4d-43e4-b428-88c2b0876399
# ╠═a78a1080-e1e6-4912-ac3b-4911f9e15564
# ╠═4fc6a508-4702-4f8a-a623-b3f5e002c8dc
# ╠═66a808cf-e28e-4e47-aef0-c460988f788e
# ╠═d2f7982d-1d7a-4e29-ad04-c4d6dec5906b
# ╠═a5494cac-b8b8-4985-9606-8d3b6999f2f1
# ╟─2ce8da75-2cdc-4d4c-80ee-8a2ac055e74c
# ╟─79994dc7-6138-4e3b-b50e-1b2769f80eb2
# ╟─ce79c1aa-1d19-4751-9199-35535e094c67
# ╟─b4ad3fa6-660b-44ec-831f-b27a794a9ed1
# ╟─7adc7b44-ee3e-41c5-9977-0a64a3fb88bb
# ╟─10923005-d271-451f-8136-2a8b2cc3b9cc
# ╟─527c224d-fa2d-4b01-b17e-501ab62f6167
# ╟─e3eb5365-4c29-4eda-b0b4-7b92d17430e8
# ╟─9eb41dc9-dfb0-489a-89d1-e9dc47f0327b
# ╟─036bed20-ad56-4913-b9d2-473d4f6773a2
# ╠═b346bcd0-9039-4175-a70d-59a2856ccb23
# ╠═aeae9d17-cb45-4619-870c-a2b8ce23ee3f
# ╟─d15a77f2-6d34-498e-8e4e-8d9f92a477c5
# ╟─742ff2e7-84ba-48a0-8ffb-6c72273ee427
# ╟─bcf5f6aa-96b6-4b6c-bc14-e99bd1647280
# ╟─3a3ad8f4-fd72-4b0a-9ba3-400c5136eec5
# ╠═0f0fbb49-a724-4804-ba8a-cc31734c17dd
# ╠═a76c2d1f-faf2-453a-ab90-2b2bb72ecc13
# ╠═b1373bbd-9b4b-40ee-b0d2-cec3d5c0dcdf
# ╠═395b37e7-b0ea-4c20-8a9a-33dd2cf195de
# ╠═0e9d12f1-92f4-4b88-b593-d6afcca8708a
