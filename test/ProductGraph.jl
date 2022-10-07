using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test
import Graphs: SimpleGraphs.SimpleEdgeIter

@testset "ProductGraph Interface" begin
    g1 = MetaGraph(smilestomol("c1ccccc1"))
    g2 = MetaGraph(smilestomol("c1ncccc1"))

    for type in [Direct, Modular]
        empty = ProductGraph{type}(SimpleGraph())
        @test nv(empty) == ne(empty) == 0
        @test typeof(edges(empty)) <: SimpleEdgeIter
        @test typeof(vertices(empty)) <: Base.OneTo
        @test typeof(props(empty, 1)) <: Dict
        @test_throws ErrorException props(empty, 1, 2)
        @test weighttype(empty) == Int
    end
end
