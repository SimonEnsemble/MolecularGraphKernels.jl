using IOCapture, MolecularGraphKernels, Test

@testset verbose=true "misc" begin
    @testset "display" begin
        g₁ = smilestomol("c1ccccc1")
        g₂ = smilestomol("c1ncccc1")
        dpg = ProductGraph{Factor}(g₁, g₂)

        captured_text = IOCapture.capture() do
            display(dpg)
        end.output

        @test contains(captured_text, "NODES") && contains(captured_text, "EDGES")
    end

    @testset "isequal" begin
        g₁ = MetaGraph(smilestomol("CN(CC(=O)O)C(=N)N"))
        g₂ = MetaGraph(smilestomol("C(CCN)C[C@@H](C(=O)O)N"))
        @test isequal(ProductGraph{Direct}(g₁, g₂), ProductGraph{Direct}(ProductGraph{Factor}(g₁, g₂)))
    end
end
