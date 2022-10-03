using IOCapture, MolecularGraphKernels, Test

@testset verbose=true "misc" begin
    @testset "display" begin
        g₁ = smilestomol("c1ccccc1")
        g₂ = smilestomol("c1ncccc1")
        dpg = ProductGraph{Modular}(g₁, g₂)

        captured_text = IOCapture.capture() do
            display(dpg)
        end.output

        @test contains(captured_text, "NODES") && contains(captured_text, "EDGES")
    end

    @testset "isequal" begin
        g₁ = MetaGraph(smilestomol("CN(CC(=O)O)C(=N)N"))
        g₂ = MetaGraph(smilestomol("C(CCN)C[C@@H](C(=O)O)N"))
        dpg = ProductGraph{Direct}(g₁, g₂)
        @test isequal(dpg, ProductGraph{Direct}(ProductGraph{Modular}(g₁, g₂)))
        g₃ = MetaGraph(smilestomol("c1ccccc1"))
        @test !isequal(dpg, ProductGraph{Direct}(g₂, g₃))
        dpg2 = deepcopy(dpg)
        set_prop!(dpg2.graph, 2, :label, 15)
        @test !isequal(dpg, dpg2)
        dpg2 = deepcopy(dpg)
        @assert add_edge!(dpg2.graph, 6, 9, Dict(:label => 1))
        @test !isequal(dpg, dpg2)
        dpg2 = deepcopy(dpg)
        set_prop!(dpg2.graph, first(edges(dpg2.graph)), :label, 0)
        @test !isequal(dpg, dpg2)
    end
end
