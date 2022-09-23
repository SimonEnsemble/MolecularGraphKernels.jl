using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

@testset "graph_conversion" begin
    A = smilestomol("C1CC=1")
    B = MetaGraph(3)

    add_edge!(B, 1, 2, Dict(:order => :single))
    add_edge!(B, 2, 3, Dict(:order => :single))
    add_edge!(B, 3, 1, Dict(:order => :double))

    set_prop!(B, 1, :symbol, :C)
    set_prop!(B, 2, :symbol, :C)
    set_prop!(B, 3, :symbol, :C)

    @test MetaGraph(A) == B
end