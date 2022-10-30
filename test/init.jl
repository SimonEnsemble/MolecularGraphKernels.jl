module Test_init

using Test

@testset "__init__" begin
    using MolecularGraphKernels
    @test length(MolecularGraphKernels._maccs_queries) == 166
end

end
