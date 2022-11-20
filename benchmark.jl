using BenchmarkTools, JLD2, MolecularGraph, MolecularGraphKernels

# include("grakel_comparison.jl")

@load "data.jld2"
graphs = data[:graphs]
# graphs = vcat(graphs, graphs)

@time gram_matrix(ccsi, graphs)

# @time grakel_compute("SubgraphMatching(k=999, lw='uniform')", graphs; n=1)
