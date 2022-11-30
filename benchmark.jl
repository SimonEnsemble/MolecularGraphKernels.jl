@info "Starting up"
using DataFrames, JLD2, MolecularGraphKernels
include("grakel_comparison.jl")

using Distributed
@everywhere using ProgressMeter

@info "Loading data"
@load "data.jld2"
graphs = data[:graphs]          #  49 graphs
graphs = vcat(graphs, graphs)   #  98 graphs
graphs = vcat(graphs, graphs)   # 196 graphs

@info "Starting computation"
@time gram_matrix(ccsi, graphs)
@time grakel_compute("SubgraphMatching(k=999, lw='uniform')", graphs; n=1)

@info "Done!"
