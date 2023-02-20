module Test_gram_matrix

using MolecularGraph, MolecularGraphKernels, Test
import MolecularGraphKernels: dump_on_error

@testset verbose = true "Gram matrix" begin
    graphs = MetaGraph.(smilestomol.([
        "c1ccccc1"
        "NC=O"
        "C(NC=O)NC=O"
    ]))

    @testset "Random Walk Kernel" begin
        l = 4
        random_walk_gm = gram_matrix(random_walk, graphs; l=l)

        @test size(random_walk_gm) == (3, 3)
        @test random_walk_gm[1, 2] ==
              random_walk_gm[2, 1] ==
              random_walk(graphs[1], graphs[2]; l=l)
        @test random_walk_gm[3, 2] ==
              random_walk_gm[2, 3] ==
              random_walk(graphs[3], graphs[2]; l=l)
        @test random_walk_gm[1, 1] == random_walk(graphs[1], graphs[1]; l=l)
    end

    @testset "Connected Common Subgraph Isomorphism Kernel" begin
        csi_gm = gram_matrix(ccsi, graphs)

        @test size(csi_gm) == (3, 3)
        @test csi_gm[1, 2] == csi_gm[2, 1] == ccsi(graphs[1], graphs[2])
        @test csi_gm[3, 2] == csi_gm[3, 2] == ccsi(graphs[3], graphs[2])
        @test csi_gm[2, 2] == ccsi(graphs[2], graphs[2])
        @test csi_gm == gram_matrix(ccsi, graphs; max_depth=10)
        csi_gm2 = gram_matrix(ccsi, graphs; max_depth=3)
        @test csi_gm ≠ csi_gm2
        idx = [[2, 2], [2, 3], [3, 2]]
        @test all(csi_gm[i...] == csi_gm2[i...] for i in idx)
    end

    @testset "Normalization" begin
        K = gram_matrix(random_walk, graphs; l=3)
        @test gram_matrix(random_walk, graphs; l=3, normalize=true) == gm_norm(K)
    end

    @testset "Time Limit" begin
        @test_throws ErrorException gram_matrix(random_walk, graphs; l=4, max_runtime=0)
    end

    @testset "Error Recovery" begin
        K = ones(20, 20)
        @test isnothing(dump_on_error(K))
        K[rand(1:20), rand(1:20)] = -1
        @test_throws ErrorException dump_on_error(K)
        @test length(filter(contains("gram_matrix_dump_jl_"), readdir())) > 0
    end
end

end
