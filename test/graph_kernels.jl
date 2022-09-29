using Graphs, MetaGraphs, MolecularGraphKernels, Test

@testset "random_walk_kernel" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => 1))
    add_edge!(A, 2, 3, Dict(:label => 1))
    add_edge!(A, 3, 1, Dict(:label => 1))
    set_prop!(A, 1, :label, 6)
    set_prop!(A, 2, :label, 6)
    set_prop!(A, 3, :label, 6)

    B = MetaGraph(3)
    add_edge!(B, 1, 2, Dict(:label => 1))
    add_edge!(B, 2, 3, Dict(:label => 1))
    add_edge!(B, 3, 1, Dict(:label => 2))
    set_prop!(B, 1, :label, 6)
    set_prop!(B, 2, :label, 6)
    set_prop!(B, 3, :label, 6)

    l = 4

    x = random_walk_kernel(dpg_adj_mat(A, B), l)

    @test x == random_walk_kernel(direct_product_graph(A, B), l)
    @test x == random_walk_kernel(A, B, l)

    mol = smilestomol("c1nc(C)ccc1")
    g = MetaGraph(mol)
    adj_mat = dpg_adj_mat(mol, g)
    x = random_walk_kernel(adj_mat, l)

    @test x == random_walk_kernel(g, g, l)
    @test x == random_walk_kernel(mol, mol, l)
    @test x == random_walk_kernel(g, mol, l)
    @test x == random_walk_kernel(mol, g, l)
end
