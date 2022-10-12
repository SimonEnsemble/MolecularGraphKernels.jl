using IOCapture, MolecularGraphKernels, Test

@testset verbose = true "misc" begin
    @testset "display" begin
        g₁ = smilestomol("c1ccccc1")
        g₂ = smilestomol("c1ncccc1")
        dpg = ProductGraph{Modular}(g₁, g₂)

        captured_text = IOCapture.capture() do
            display(dpg)
            return
        end.output

        @test contains(captured_text, "NODES") && contains(captured_text, "EDGES")
    end
end
