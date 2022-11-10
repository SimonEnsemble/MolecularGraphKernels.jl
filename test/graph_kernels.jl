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

    @testset "ethanol-acetic acid test (EZ to count)" begin
        # the two input molecular graphs
        g₁ = MetaGraph(smilestomol("CC(O)=O")) # acetic acid
        g₂ = MetaGraph(smilestomol("CCO")) # ethanol

        ## random walks
        # l = 1. C <--> C * 4 + O <--> O * 2
        @test random_walk(g₁, g₂; l=1) == 6
        # l = 2. C - C <--> C -- C * 4 + C - O <--> C - O * 2 + O - C <--> O - C * 2
        @test random_walk(g₁, g₂; l=2) == 8

        ## matching connected subgraphs
        @test subgraph_matching(g₁, g₂) == 1
    end

    @testset "Common Subgraph Isomorphism" begin
        x = common_subgraph_isomorphism(g₁, g₂)
        @test x == 15
        @test x == common_subgraph_isomorphism(ProductGraph{Modular}(g₁, g₂))
        @test common_subgraph_isomorphism(g₁, g₂; λ=length) == 26
    end

    @testset "Subgraph Matching" begin
        mpg = ProductGraph{Modular}(g₁, g₂)
        x = subgraph_matching(mpg, _ -> 1)
        @test x == 15
        @test x == common_subgraph_isomorphism(mpg)
        @test subgraph_matching(mpg, length) == 26
        @test subgraph_matching(mpg, _ -> 1; c_cliques=true) == 11
        @test subgraph_matching(mpg, length; c_cliques=true) == 16
        @test_throws AssertionError subgraph_matching(mpg, length; c_cliques=false)
    end
end

end
