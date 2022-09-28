### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 7b2b99d0-3dea-11ed-18db-0765cfcce7e1
begin
	import IOCapture, Pkg
	IOCapture.capture() do
		Pkg.activate(".")
		Pkg.instantiate()
		Pkg.add("PlutoTest")
	end
	using BenchmarkTools, MolecularGraph, MolecularGraphKernels, PlutoTest
end

# ╔═╡ 91e1dc54-e924-451b-9ba6-00f3a7bb3fad


# ╔═╡ f1ba7ec0-eabb-4f2e-a213-e35913d32aaa
begin
	local mol = smilestomol("c1ccccc1-c2ccccc2")
	g = MetaGraph(mol)
	HTML(drawsvg(mol, 250, 250))
end

# ╔═╡ 7fa366f5-4b19-42e7-a6aa-0474ce2fe2be
begin
	local mol = smilestomol("c1cccnc1-c2ncccc2")
	h = MetaGraph(mol)
	HTML(drawsvg(mol, 250, 250))
end

# ╔═╡ bcc14d30-e388-408c-ac0e-408575d7a8ae
@btime dpg_adj_mat(g, h)

# ╔═╡ 015d5374-5bd5-48b6-a43c-4d8ae538ec9b
fixed_length_rw_kernel(g, h, 4)

# ╔═╡ 5c2f65cc-7e16-4f2b-a7e5-7ded1cfc5996
@btime direct_product_graph(g, h)

# ╔═╡ Cell order:
# ╠═91e1dc54-e924-451b-9ba6-00f3a7bb3fad
# ╠═7b2b99d0-3dea-11ed-18db-0765cfcce7e1
# ╠═f1ba7ec0-eabb-4f2e-a213-e35913d32aaa
# ╠═7fa366f5-4b19-42e7-a6aa-0474ce2fe2be
# ╠═bcc14d30-e388-408c-ac0e-408575d7a8ae
# ╠═015d5374-5bd5-48b6-a43c-4d8ae538ec9b
# ╠═5c2f65cc-7e16-4f2b-a7e5-7ded1cfc5996
