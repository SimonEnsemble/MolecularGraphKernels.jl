module Test_misc

using Graphs, IOCapture, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test

@testset verbose = true "misc" begin
    @testset "display" begin
        g₁ = smilestomol("c1ccccc1")
        g₂ = smilestomol("c1ncccc1")
        dpg = ProductGraph{Direct}(g₁, g₂)

        captured_text = IOCapture.capture() do
            display(dpg)
            return
        end.output

        @test contains(captured_text, "NODES") && contains(captured_text, "EDGES")
    end

    @testset "MPG isomorphic subgraphs" begin
        g₁ = smilestomol("NC=O")
        g₂ = smilestomol("CN(C=O)C=O")
        mpg = ProductGraph{Modular}(g₁, g₂)
        imsgs = isomorphic_subgraphs(mpg; min_size=2)

        manual_isom1 = MetaGraph(3)
        set_prop!(manual_isom1, 1, :label, 7)
        set_prop!(manual_isom1, 2, :label, 6)
        set_prop!(manual_isom1, 3, :label, 8)
        add_edge!(manual_isom1, 1, 2, Dict(:label => 1))
        add_edge!(manual_isom1, 1, 3, Dict(:label => 0))
        add_edge!(manual_isom1, 3, 2, Dict(:label => 2))
        manual_isom1 = ProductGraph{Modular}(manual_isom1)

        manual_isom3 = MetaGraph(2)
        set_prop!(manual_isom3, 1, :label, 6)
        set_prop!(manual_isom3, 2, :label, 7)
        add_edge!(manual_isom3, 1, 2, Dict(:label => 1))
        manual_isom3 = ProductGraph{Modular}(manual_isom3)

        @test length(imsgs) == 3
        @test is_isomorphic(imsgs[1], manual_isom1)
        @test is_isomorphic(imsgs[2], manual_isom1)
        @test is_isomorphic(imsgs[3], manual_isom3)
    end
end

end
