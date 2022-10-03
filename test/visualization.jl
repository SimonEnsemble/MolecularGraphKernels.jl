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

    @testset "Modular Product Graph" begin
        @test !isnothing(viz_graph(ProductGraph{Modular}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")); savename="mpg"))
        vis = isfile("mpg.pdf")
        @test vis
        if vis
            @info "Modular product graph visualization in mpg.pdf"
        end
    end

    @testset "Graph Plot Styles" begin
        for style in [:circular, :spectral]
            @test !isnothing(viz_graph(smilestomol("C(NC=O)NC=O"); savename="$style", layout_style=style))
            vis = isfile("$style.pdf")
            @test vis
            if vis
                @info "Modular product graph visualization in $style.pdf"
            end
        end
    end
end
