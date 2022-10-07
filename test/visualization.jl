using MolecularGraph, MolecularGraphKernels, Test

function test_vis(graph, graph_name, set_name, style=:spring)
    @testset "$set_name" begin
        rm("$graph_name.pdf", force=true)
        @test !isnothing(viz_graph(graph; savename=graph_name, layout_style=style))
        vis = isfile("$graph_name.pdf")
        @test vis
        if vis
            @info "$set_name test visualization in $graph_name.pdf"
        end
    end
end

@testset verbose=true "Graph Types" begin
    test_vis(smilestomol("c1ccccc1"), "benzene", "MetaGraph")
    test_vis(ProductGraph{Direct}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")), "dpg", "Direct Product Graph")
    test_vis(ProductGraph{Modular}(smilestomol("NC=O"), smilestomol("C(NC=O)NC=O")), "mpg", "Modular Product Graph")
end

@testset verbose=true "Graph Plot Styles" begin
    mol = smilestomol("C(NC=O)NC=O")
    for style in [:spring, :circular, :spectral]
        test_vis(mol, "$style", "$style", style)
    end
    @test_throws ErrorException viz_graph(mol, layout_style=:bogus_style)
end