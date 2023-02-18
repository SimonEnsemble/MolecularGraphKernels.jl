module Test_init

using Test

@testset "__init__" begin
    using MolecularGraphKernels
    if !Sys.iswindows()
        @test length(MolecularGraphKernels._maccs_queries) == 166
    else
        @test isa(MolecularGraphKernels, Module)
    end
end

end
