using MolecularGraphKernels, Test

@testset "visualization" begin
    @test !isnothing(viz_graph(smilestomol("c1ccccc1"); savename="benzene"))
    @test isfile("benzene.pdf")
    @info "Benzene test visualization in benzene.pdf"

    @test !isnothing(viz_graph(ProductGraph{Direct}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")); savename="dpg"))
    @test isfile("dpg.pdf")
    @info "Direct product graph visualization in dpg.pdf"
end
