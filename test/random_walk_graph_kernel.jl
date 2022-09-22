using Graphs, MolecularGraphKernels, Test

@testset "fixed_length_rw_kernel" begin
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

    l = 3

    s1 = fixed_length_rw_kernel(A, B, l)

    adj_mat = direct_product_graph(A, B)
    
    s2 = fixed_length_rw_kernel(adj_mat, l)

    dpg = MetaGraph(adj_mat)

    s3 = fixed_length_rw_kernel(dpg, l)

    @test s1 == s2 == s3
end
