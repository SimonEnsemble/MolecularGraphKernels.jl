using MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

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