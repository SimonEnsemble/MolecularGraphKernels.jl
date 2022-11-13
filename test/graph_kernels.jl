module Test_graph_kernels

using MolecularGraph, MolecularGraphKernels, Test

@testset verbose = true "Graph Kernels" begin
    g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])

    @testset "Random Walk" begin
        l = 4
        x = random_walk(product_graph_adjacency_matrix(Direct, g₁, g₂); l=l)
        @test x == 74
        @test x == random_walk(ProductGraph{Direct}(g₁, g₂); l=l)
        @test x == random_walk(g₁, g₂; l=l)
    end

    @testset "Connected Common Subgraph Isomorphism" begin
        mpg = ProductGraph{Modular}(g₁, g₂)
        @test ccsi(mpg) == 13
        @test ccsi(mpg; λ=length) == 22
    end
end

end
