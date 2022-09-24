using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

@testset "graph_conversion" begin
    A = smilestomol("C1CC=1")
    B = MetaGraph(3)

    add_edge!(B, 1, 2, Dict(:label => :single))
    add_edge!(B, 2, 3, Dict(:label => :single))
    add_edge!(B, 3, 1, Dict(:label => :double))

    set_prop!(B, 1, :label, :C)
    set_prop!(B, 2, :label, :C)
    set_prop!(B, 3, :label, :C)

    @test MetaGraph(A) == B
end