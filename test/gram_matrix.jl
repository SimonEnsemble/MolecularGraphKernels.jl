module Test_gram_matrix

using MolecularGraph, MolecularGraphKernels, Test

@testset verbose = true "Gram matrix" begin
    graphs = MetaGraph.(smilestomol.([
        "c1ccccc1"
        "NC=O"
        "C(NC=O)NC=O"
    ]))

    @testset "Random Walk Kernel" begin
        l = 4
        rwk_gm = gram_matrix(random_walk, graphs; l=l)

        @test size(rwk_gm) == (3, 3)
        @test rwk_gm[1, 2] == rwk_gm[2, 1] == random_walk(graphs[1], graphs[2]; l=l)
        @test rwk_gm[3, 2] == rwk_gm[2, 3] == random_walk(graphs[3], graphs[2]; l=l)
        @test rwk_gm[1, 1] == random_walk(graphs[1], graphs[1]; l=l)
    end

    @testset "Common Subgraph Isomorphism Kernel" begin
        csi_gm = gram_matrix(common_subgraph_isomorphism, graphs)

        @test size(csi_gm) == (3, 3)
        @test csi_gm[1, 2] ==
              csi_gm[2, 1] ==
              common_subgraph_isomorphism(graphs[1], graphs[2])
        @test csi_gm[3, 2] ==
              csi_gm[3, 2] ==
              common_subgraph_isomorphism(graphs[3], graphs[2])
        @test csi_gm[2, 2] == common_subgraph_isomorphism(graphs[2], graphs[2])
    end

    @testset "Normalization" begin
        K = gram_matrix(random_walk, graphs; l=3)
        @test gram_matrix(random_walk, graphs; l=3, normalize=true) == gm_norm(K)
    end
end

end
