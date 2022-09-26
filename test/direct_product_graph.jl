using Combinatorics, Graphs, MetaGraphs, MolecularGraphKernels, SparseArrays, Test

include("check_isom.jl")

@testset "dpg_adj_mat" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => :single))
    add_edge!(A, 2, 3, Dict(:label => :single))
    add_edge!(A, 3, 1, Dict(:label => :single))
    set_prop!(A, 1, :label, :C)
    set_prop!(A, 2, :label, :C)
    set_prop!(A, 3, :label, :C)

    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, :double)
    
    @test dpg_adj_mat(A, B) == sparse(
        [5, 6, 4, 6, 4, 5, 2, 3, 8, 9, 1, 3, 7, 9, 1, 2, 7, 8, 5, 6, 4, 6, 4, 5], 
        [1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9], 
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        9, 9
    )

    C = MetaGraph(2)
    add_edge!(C, 1, 2, Dict(:label => 1))
    set_prop!(C, 1, :label, :H)
    set_prop!(C, 2, :label, :H)
    h2 = smilestomol("[H]-[H]")

    x = dpg_adj_mat(C, C)

    @test x == dpg_adj_mat(h2, h2)
    @test x == dpg_adj_mat(h2, C)
    @test x == dpg_adj_mat(C, h2)
end

@testset "direct_product_graph" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => :single))
    add_edge!(A, 2, 3, Dict(:label => :single))
    add_edge!(A, 3, 1, Dict(:label => :single))
    set_prop!(A, 1, :label, :C)
    set_prop!(A, 2, :label, :C)
    set_prop!(A, 3, :label, :C)

    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, :double)

    @test adjacency_matrix(direct_product_graph(A, B)) == dpg_adj_mat(A, B)

    C = MetaGraph(2)
    add_edge!(C, 1, 2, Dict(:label => 1))
    set_prop!(C, 1, :label, :H)
    set_prop!(C, 2, :label, :H)
    h2 = smilestomol("[H]-[H]")

    @test direct_product_graph(C, C) == direct_product_graph(MetaGraph(h2), C)
    @test direct_product_graph(h2, h2) == direct_product_graph(h2, C)
    @test direct_product_graph(h2, C) == direct_product_graph(C, h2)

    g = MetaGraph(smilestomol("COP(=O)(OC)OC(Br)C(Cl)(Cl)Br"))
    h = MetaGraph(smilestomol("COP(N)(=O)SC"))

    gxh = MetaGraph(17)
    set_prop!(gxh, 1, :label, :P)
    set_prop!(gxh, 2, :label, :O)
    set_prop!(gxh, 3, :label, :O)
    set_prop!(gxh, 4, :label, :O)
    set_prop!(gxh, 5, :label, :O)
    set_prop!(gxh, 6, :label, :O)
    set_prop!(gxh, 7, :label, :O)
    set_prop!(gxh, 8, :label, :O)
    set_prop!(gxh, 9, :label, :O)
    set_prop!(gxh, 10, :label, :C)
    set_prop!(gxh, 11, :label, :C)
    set_prop!(gxh, 12, :label, :C)
    set_prop!(gxh, 13, :label, :C)
    set_prop!(gxh, 14, :label, :C)
    set_prop!(gxh, 15, :label, :C)
    set_prop!(gxh, 16, :label, :C)
    set_prop!(gxh, 17, :label, :C)
    add_edge!(gxh, 1, 2, Dict(:label => 2))
    add_edge!(gxh, 1, 3, Dict(:label => 1))
    add_edge!(gxh, 1, 4, Dict(:label => 1))
    add_edge!(gxh, 1, 5, Dict(:label => 1))
    add_edge!(gxh, 3, 10, Dict(:label => 1))
    add_edge!(gxh, 4, 11, Dict(:label => 1))
    add_edge!(gxh, 5, 12, Dict(:label => 1))

    @test is_isomorphic(direct_product_graph(g, h), gxh)
end

@testset "SMILES flexibility" begin
    A = smilestomol("c1c(C)cccc1")
    B = smilestomol("C1C=C(C)C=CC=1")
    C = smilestomol("C1=C(C)C=CC=C1")
    D = smilestomol("c1-c-c(C)-c-c-c1")
    E = smilestomol("c1ccc(C)cc1")
    A, B, C, D, E = MetaGraph.([A, B, C, D, E])

    for (x, y) in with_replacement_combinations([A, B, C, D, E], 2)
        @test is_isomorphic(direct_product_graph(x, x), direct_product_graph(x, y))
    end
end
