using IOCapture, MolecularGraphKernels, Test

@testset verbose=true "misc" begin
    @testset "display" begin
        reference_text = """---NODES---
        [1] v₁v₂_pair:(1, 1) label:6
        [2] v₁v₂_pair:(2, 1) label:6
        [3] v₁v₂_pair:(3, 1) label:6
        [4] v₁v₂_pair:(4, 1) label:6
        [5] v₁v₂_pair:(5, 1) label:6
        [6] v₁v₂_pair:(6, 1) label:6
        [7] v₁v₂_pair:(1, 3) label:6
        [8] v₁v₂_pair:(2, 3) label:6
        [9] v₁v₂_pair:(3, 3) label:6
        [10] v₁v₂_pair:(4, 3) label:6
        [11] v₁v₂_pair:(5, 3) label:6
        [12] v₁v₂_pair:(6, 3) label:6
        [13] v₁v₂_pair:(1, 4) label:6
        [14] v₁v₂_pair:(2, 4) label:6
        [15] v₁v₂_pair:(3, 4) label:6
        [16] v₁v₂_pair:(4, 4) label:6
        [17] v₁v₂_pair:(5, 4) label:6
        [18] v₁v₂_pair:(6, 4) label:6
        [19] v₁v₂_pair:(1, 5) label:6
        [20] v₁v₂_pair:(2, 5) label:6
        [21] v₁v₂_pair:(3, 5) label:6
        [22] v₁v₂_pair:(4, 5) label:6
        [23] v₁v₂_pair:(5, 5) label:6
        [24] v₁v₂_pair:(6, 5) label:6
        [25] v₁v₂_pair:(1, 6) label:6
        [26] v₁v₂_pair:(2, 6) label:6
        [27] v₁v₂_pair:(3, 6) label:6
        [28] v₁v₂_pair:(4, 6) label:6
        [29] v₁v₂_pair:(5, 6) label:6
        [30] v₁v₂_pair:(6, 6) label:6
        ---EDGES---
        (1, 26) label:-1
        (1, 30) label:-1
        (2, 25) label:-1
        (2, 27) label:-1
        (3, 26) label:-1
        (3, 28) label:-1
        (4, 27) label:-1
        (4, 29) label:-1
        (5, 28) label:-1
        (5, 30) label:-1
        (6, 25) label:-1
        (6, 29) label:-1
        (7, 14) label:-1
        (7, 18) label:-1
        (8, 13) label:-1
        (8, 15) label:-1
        (9, 14) label:-1
        (9, 16) label:-1
        (10, 15) label:-1
        (10, 17) label:-1
        (11, 16) label:-1
        (11, 18) label:-1
        (12, 13) label:-1
        (12, 17) label:-1
        (13, 20) label:-1
        (13, 24) label:-1
        (14, 19) label:-1
        (14, 21) label:-1
        (15, 20) label:-1
        (15, 22) label:-1
        (16, 21) label:-1
        (16, 23) label:-1
        (17, 22) label:-1
        (17, 24) label:-1
        (18, 19) label:-1
        (18, 23) label:-1
        (19, 26) label:-1
        (19, 30) label:-1
        (20, 25) label:-1
        (20, 27) label:-1
        (21, 26) label:-1
        (21, 28) label:-1
        (22, 27) label:-1
        (22, 29) label:-1
        (23, 28) label:-1
        (23, 30) label:-1
        (24, 25) label:-1
        (24, 29) label:-1"""

        captured_text = IOCapture.capture() do
            display(dpg)
        end.output

        @test contains(captured_text, reference_text)
    end

    @testset "isequal" begin
        g₁ = MetaGraph(smilestomol("CN(CC(=O)O)C(=N)N"))
        g₂ = MetaGraph(smilestomol("O=C1c2ncnc2nc(N)N1"))
        @test ProductGraph{Direct}(g₁, g₂) == ProductGraph{Direct}(ProductGraph{Factor}(g₁, g₂))
    end
end
