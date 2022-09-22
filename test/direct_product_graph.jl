using Graphs, MetaGraphs, MolecularGraphKernels, SparseArrays, Test

@testset "direct_product_graph" begin
    A = MetaGraph(3)
    B = MetaGraph(3)

    add_edge!(A, 1, 2, Dict(:order => :single))
    add_edge!(A, 2, 3, Dict(:order => :single))
    add_edge!(A, 3, 1, Dict(:order => :single))

    add_edge!(B, 1, 2, Dict(:order => :single))
    add_edge!(B, 2, 3, Dict(:order => :single))
    add_edge!(B, 3, 1, Dict(:order => :double))

    set_prop!(A, 1, :symbol, :C)
    set_prop!(A, 2, :symbol, :C)
    set_prop!(A, 3, :symbol, :C)

    set_prop!(B, 1, :symbol, :C)
    set_prop!(B, 2, :symbol, :C)
    set_prop!(B, 3, :symbol, :C)

    dpg = direct_product_graph(A, B)

    @test dpg == sparse(
        [5, 6, 4, 6, 4, 5, 2, 3, 8, 9, 1, 3, 7, 9, 1, 2, 7, 8, 5, 6, 4, 6, 4, 5], 
        [1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9], 
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        9, 9
    )

    C = MetaGraph(2)

    add_edge!(C, 1, 2, Dict(:order => 1.))
    set_prop!(C, 1, :symbol, :H)
    set_prop!(C, 2, :symbol, :H)

    @test_throws AssertionError direct_product_graph(A, C)
end
