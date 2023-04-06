module Test_check_isom

using MetaGraphs, MolecularGraph, MolecularGraphKernels, Test
import MolecularGraphKernels: is_isomorphic

@testset verbose = true "isomorphism detection" begin
    @testset "MetaGraphs" begin
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

    @testset "Adjacency matrix vs ProductGraph" begin
        g₁ = MetaGraph(smilestomol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"))
        g₂ = MetaGraph(smilestomol("O=C1[C@H](C)[C@@H]2[C@](C(C)C)(C1)C2"))
        for type in [Modular, Direct]
            @test is_isomorphic(GraphMatrix{type}(g₁, g₂), ProductGraph{type}(g₁, g₂))
        end
    end
end

end
