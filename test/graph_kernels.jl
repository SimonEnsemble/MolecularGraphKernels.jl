using Graphs, MetaGraphs, MolecularGraphKernels, Test

@testset verbose = true "Graph Kernels" begin
    g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
    
    @testset "Random Walk" begin
        l = 4
        x = random_walk(product_graph_adjacency_matrix(Direct, g₁, g₂); l=l)
        @test x == 74
        @test x == random_walk(ProductGraph{Direct}(g₁, g₂); l=l)
        @test x == random_walk(g₁, g₂; l=l)
    end

    @testset "Common Subgraph Isomorphism" begin
        x = common_subgraph_isomorphism(g₁, g₂)
        @test x == 3
        @test x == common_subgraph_isomorphism(ProductGraph{Modular}(g₁, g₂))
        @test x ==
              common_subgraph_isomorphism(product_graph_adjacency_matrix(Modular, g₂, g₁))
        @test common_subgraph_isomorphism(g₁, g₂; λ=length) == 8
    end
end
