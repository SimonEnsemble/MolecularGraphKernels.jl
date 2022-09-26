using Graphs, MetaGraphs, MolecularGraphKernels, Test

@testset "fixed_length_rw_kernel" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => :single))
    add_edge!(A, 2, 3, Dict(:label => :single))
    add_edge!(A, 3, 1, Dict(:label => :single))
    set_prop!(A, 1, :label, :C)
    set_prop!(A, 2, :label, :C)
    set_prop!(A, 3, :label, :C)

    B = MetaGraph(3)
    add_edge!(B, 1, 2, Dict(:label => :single))
    add_edge!(B, 2, 3, Dict(:label => :single))
    add_edge!(B, 3, 1, Dict(:label => :double))
    set_prop!(B, 1, :label, :C)
    set_prop!(B, 2, :label, :C)
    set_prop!(B, 3, :label, :C)

    l = 4

    x = fixed_length_rw_kernel(dpg_adj_mat(A, B), l)

    @test x == fixed_length_rw_kernel(direct_product_graph(A, B), l)
    @test x == fixed_length_rw_kernel(A, B, l)

    mol = smilestomol("c1nc(C)ccc1")
    g = MetaGraph(mol)
    adj_mat = dpg_adj_mat(mol, g)
    x = fixed_length_rw_kernel(adj_mat, l)

    @test x == fixed_length_rw_kernel(g, g, l)
    @test x == fixed_length_rw_kernel(mol, mol, l)
    @test x == fixed_length_rw_kernel(g, mol, l)
    @test x == fixed_length_rw_kernel(mol, g, l)
end
