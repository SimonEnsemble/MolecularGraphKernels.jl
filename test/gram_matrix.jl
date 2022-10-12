using MolecularGraphKernels, Test

@testset verbose = true "Gram matrix" begin
    graphs = MetaGraph.(smilestomol.([
        "c1ccccc1"
        "NC=O"
        "C(NC=O)NC=O"
    ]))

    l = 4
    rwk_gm = gram_matrix(random_walk, graphs; l=l)

    @test size(rwk_gm) == (3, 3)
    @test rwk_gm[1, 2] == rwk_gm[2, 1] == random_walk(graphs[1], graphs[2]; l=l)
    @test rwk_gm[3, 2] == rwk_gm[2, 3] == random_walk(graphs[3], graphs[2]; l=l)
    @test rwk_gm[1, 1] == random_walk(graphs[1], graphs[1]; l=l)
end
