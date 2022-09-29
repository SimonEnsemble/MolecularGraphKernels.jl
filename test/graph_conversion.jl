using Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

include("check_isom.jl")

@testset "is_isomporphic" begin
    A = MetaGraph(smilestomol("C1CC=1"))
    
    # trivial test (identical graphs)
    B = deepcopy(A)
    @test is_isomorphic(A, B)

    # less trivial test (isomorphic but not identical)
    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, 1)
    set_prop!(B, 1, 2, :label, 2)
    @test is_isomorphic(A, B)
    
    # isomorphism broken by changing atom identity
    B = deepcopy(A)
    set_prop!(B, 2, :label, 7)
    @test !is_isomorphic(A, B)

    # isomorphism broken by changing bond order
    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, 1)
    @test !is_isomorphic(A, B)
end

@testset "graph_conversion" begin
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