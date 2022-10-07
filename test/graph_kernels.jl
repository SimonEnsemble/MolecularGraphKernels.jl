using Graphs, MetaGraphs, MolecularGraphKernels, Test

@testset verbose=true "graph kernels" begin
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

    l = 4

    @testset "random walk" begin
        x = random_walk(product_graph_adjacency_matrix(Direct, A, B), l)

        @test x == random_walk(ProductGraph{Direct}(A, B), l)
        @test x == random_walk(A, B, l)
    
        @test random_walk(ProductGraph{Direct}(mol, g), l) == random_walk(g, g, l)
    end
end
