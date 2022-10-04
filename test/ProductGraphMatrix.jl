using MolecularGraphKernels, Test

@testset "ProductGraphMatrix Interface" begin
    for type in [Direct, Modular]
        foo = ProductGraphMatrix{type}(8)
        foo[3, 6] = true
        foo[4, 6] = true
        foo[3, 4] = true
        foo = foo .|| foo'

        @test foo[3, 6] && foo[4, 6] && foo[3, 4] && foo[6, 3] && foo[6, 4] && foo[4, 3]
        @test findall(foo[:, 3]) == [4, 6]
        @test sum(foo) == 6
        @test sum(foo^3) == 24
    end
end
