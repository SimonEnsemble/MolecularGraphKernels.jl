using MolecularGraph, MolecularGraphKernels, Test

@testset verbose=true "visualization" begin
    @testset "MetaGraph" begin
        @test !isnothing(viz_graph(smilestomol("c1ccccc1"); savename="benzene"))
        vis = isfile("benzene.pdf")
        @test vis
        if vis
            @info "Benzene test visualization in benzene.pdf"
        end
    end

    @testset "Direct Product Graph" begin
        @test !isnothing(viz_graph(ProductGraph{Direct}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")); savename="dpg"))
        vis =  isfile("dpg.pdf")
        @test vis
        if vis
            @info "Direct product graph visualization in dpg.pdf"
        end
    end

    @testset "Factor Product Graph" begin
        @test !isnothing(viz_graph(ProductGraph{Factor}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")); savename="fpg"))
        vis = isfile("fpg.pdf")
        @test vis
        if vis
            @info "Factor product graph visualization in fpg.pdf"
        end
    end
end
