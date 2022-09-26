using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

include("check_isom.jl")

@testset "is_isomporphic" begin
    A = MetaGraph(smilestomol("C1CC=1"))
    B = deepcopy(A)
    
    @test is_isomorphic(A, B)
    
    B = deepcopy(A)
    set_prop!(B, 2, :label, :N)

    @test !is_isomorphic(A, B)

    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, :single)

    @test !is_isomorphic(A, B)
end

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

    A = smilestomol("c1ccccc1")
    B = smilestomol("C1C=CC=CC=1")
    C = smilestomol("C1=CC=CC=C1")
    D = smilestomol("c1-c-c-c-c-c1")

    A, B, C, D = MetaGraph.([A, B, C, D])

    @test is_isomorphic(A, B)
    @test is_isomorphic(B, C)
    @test is_isomorphic(C, D)
end