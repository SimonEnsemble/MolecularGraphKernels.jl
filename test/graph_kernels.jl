using Graphs, MetaGraphs, MolecularGraphKernels, Test

@testset verbose = true "Graph Kernels" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => 1))
    add_edge!(A, 2, 3, Dict(:label => 1))
    add_edge!(A, 3, 1, Dict(:label => 1))
    set_prop!(A, 1, :label, 6)
    set_prop!(A, 2, :label, 6)
    set_prop!(A, 3, :label, 6)

    B = MetaGraph(3)
    add_edge!(B, 1, 2, Dict(:label => 1))
    add_edge!(B, 2, 3, Dict(:label => 1))
    add_edge!(B, 3, 1, Dict(:label => 2))
    set_prop!(B, 1, :label, 6)
    set_prop!(B, 2, :label, 6)
    set_prop!(B, 3, :label, 6)

    mol = smilestomol("c1nc(C)ccc1")
    g = MetaGraph(mol)

    @testset "Random Walk" begin
        l = 4
        x = random_walk(product_graph_adjacency_matrix(Direct, A, B); l=l)

        @test x == random_walk(ProductGraph{Direct}(A, B); l=l)
        @test x == random_walk(A, B; l=l)
        @test random_walk(ProductGraph{Direct}(mol, g); l=l) == random_walk(g, g; l=l)
    end

    @testset "Common Subgraph Isomorphism" begin
        g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
        x = subgraph_matching(g₁, g₂)
        @test x == 3
        @test x == subgraph_matching(ProductGraph{Modular}(g₁, g₂))
        @test x ==
              subgraph_matching(product_graph_adjacency_matrix(Modular, g₂, g₁))
    end
end
