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
        rwk_gm = gram_matrix(rwk, graphs; l=l)

        @test size(rwk_gm) == (3, 3)
        @test rwk_gm[1, 2] == rwk_gm[2, 1] == rwk(graphs[1], graphs[2]; l=l)
        @test rwk_gm[3, 2] == rwk_gm[2, 3] == rwk(graphs[3], graphs[2]; l=l)
        @test rwk_gm[1, 1] == rwk(graphs[1], graphs[1]; l=l)
    end

    @testset "Connected Common Subgraph Isomorphism Kernel" begin
        csi_gm = gram_matrix(ccsi, graphs)

        @test size(csi_gm) == (3, 3)
        @test csi_gm[1, 2] ==
              csi_gm[2, 1] ==
              ccsi(graphs[1], graphs[2])
        @test csi_gm[3, 2] ==
              csi_gm[3, 2] ==
              ccsi(graphs[3], graphs[2])
        @test csi_gm[2, 2] == ccsi(graphs[2], graphs[2])
    end
end

end
