module Test_graph_conversion

using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test
include("is_isomorphic.jl")
using .IsIsomorphic

@testset verbose = true "graph conversion" begin
    @testset "GraphMol to MetaGraph" begin
        # test that graph converted from GraphMol equivalent to manually generated
        A = MetaGraph(smilestomol("C1CC=1"))
        B = MetaGraph(3)
        add_edge!(B, 1, 2, Dict(:label => 1))
        add_edge!(B, 2, 3, Dict(:label => 1))
        add_edge!(B, 3, 1, Dict(:label => 2))
        set_prop!(B, 1, :label, 6)
        set_prop!(B, 2, :label, 6)
        set_prop!(B, 3, :label, 6)
        @test is_isomorphic(A, B)

        # test that multiple SMILES inputs for picoline are all isomorphs
        A = smilestomol("n1c(C)cccc1")
        B = smilestomol("C1N=C(C)C=CC=1")
        C = smilestomol("C1=C(C)N=CC=C1")
        D = smilestomol("c1-c-c-n-c(C)-c1")
        A, B, C, D = MetaGraph.([A, B, C, D])
        @test is_isomorphic(A, B)
        @test is_isomorphic(B, C)
        @test is_isomorphic(C, D)
    end

    A = smilestomol("c1cccc(CC(C)NC)c1")
    B = smilestomol("c1c(c(CCN(C)C)c[nH]2)c2ccc1")

    @testset "ProductGraph to ProductGraph" begin
        @test is_isomorphic(
            ProductGraph{Direct}(A, B),
            ProductGraph{Direct}(ProductGraph{Modular}(A, B))
        )
    end

    @testset "ProductGraph to adjacency matrix" begin
        @test is_isomorphic(
            adjacency_matrix(ProductGraph{Direct}(A, B)),
            GraphMatrix{Direct}(A, B)
        )
    end
end

end
