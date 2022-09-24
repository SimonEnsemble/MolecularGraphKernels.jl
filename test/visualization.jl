using MolecularGraphKernels, Test

@testset "visualization" begin
    @test !isnothing(viz_graph(smilestomol("c1ccccc1"); savename="benzene"))
    @test isfile("benzene.pdf")
end
