module Test_gram_matrix

using LinearAlgebra, MolecularGraph, MolecularGraphKernels, Test

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

    @testset "Kernel vector" begin
        @test kernel_vector(random_walk, graphs[1], graphs; l=4) ==
              gram_matrix(random_walk, graphs; l=4)[:, 1]
        gm = gram_matrix(random_walk, graphs; l=4)
        ngm = gm_norm(gm)
        @test kernel_vector(random_walk, graphs[1], graphs, gm; l=4) == ngm[:, 1]
        @test kernel_vector(random_walk, graphs[1], graphs, diag(gm); l=4) == ngm[:, 1]
    end
end

end
