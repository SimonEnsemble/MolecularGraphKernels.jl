using Distributed

# spin up workers
addprocs(2)

# ensure environment synchronization
@everywhere begin
    import Pkg
    Pkg.activate(".")
    Pkg.instantiate()
end

using MolecularGraph, MolecularGraphKernels, Test

@testset verbose = true "Gram matrix" begin
    graphs = MetaGraph.(smilestomol.([
        "c1ccccc1"
        "NC=O"
        "C(NC=O)NC=O"
    ]))

    @testset "Built-In Kernel (CCSI)" begin
        csi_gm = gram_matrix(ccsi, graphs)

        @test size(csi_gm) == (3, 3)
        @test csi_gm[1, 2] == csi_gm[2, 1] == ccsi(graphs[1], graphs[2])
        @test csi_gm[3, 2] == csi_gm[3, 2] == ccsi(graphs[3], graphs[2])
        @test csi_gm[2, 2] == ccsi(graphs[2], graphs[2])
    end

    @testset "Custom Kernels" begin
        kernel1(g₁, g₂) = 1
        @test sum(gram_matrix(kernel1, graphs)) == 9

        kernel2(g₁, g₂; k=2) = k
        @test sum(gram_matrix(kernel2, graphs)) == 18
        @test sum(gram_matrix(kernel2, graphs; k=4)) == 36
    end

    @testset "Normalization" begin
        K = gram_matrix(random_walk, graphs; l=3)
        @test gram_matrix(random_walk, graphs; l=3, normalize=true) == gm_norm(K)
    end
end
