module Test_graph_kernels

using MolecularGraph, MolecularGraphKernels, Test

@testset verbose = true "Graph Kernels" begin
    @testset "Random Walk" begin
        g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
        l = 4
        x = random_walk(product_graph_adjacency_matrix(Direct, g₁, g₂); l=l)
        @test x == 74
        @test x == random_walk(ProductGraph{Direct}(g₁, g₂); l=l)
        @test x == random_walk(g₁, g₂; l=l)
    end

    @testset "ethanol-acetic acid test (EZ to count)" begin
        # the two input molecular graphs
        g₁ = MetaGraph(smilestomol("CC(O)=O")) # acetic acid
        g₂ = MetaGraph(smilestomol("CCO")) # ethanol

        ## random walks
        # l = 1. C <--> C * 4 + O <--> O * 2
        @test random_walk(g₁, g₂; l=1) == 6
        # l = 2. C - C <--> C - C * 4 + C - O <--> C - O * 2 + O - C <--> O - C * 2
        @test random_walk(g₁, g₂; l=2) == 8

        ## matching connected subgraphs
        # C ↔ C x4
        # O ↔ O x2
        # C-C ↔ C-C x2
        # C-O ↔ C-O x1
        # C-C-O ↔ C-C-O x1
        @test ccsi(g₁, g₂) == 10
        @test ccsi(g₁, g₂; λ=length) == 15
    end

    @testset "Connected Common Subgraph Isomorphism" begin
        g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
        mpg = ProductGraph{Modular}(g₁, g₂)
        @test ccsi(mpg) == 13
        @test ccsi(mpg; λ=length) == 22
    end
end

end
