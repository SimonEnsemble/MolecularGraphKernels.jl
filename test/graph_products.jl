using Combinatorics, Graphs, MetaGraphs, MolecularGraphKernels, SparseArrays, Test

include("check_isom.jl")

@testset "dpg_adj_mat" begin
    A = MetaGraph(3)
    add_edge!(A, 1, 2, Dict(:label => 1))
    add_edge!(A, 2, 3, Dict(:label => 1))
    add_edge!(A, 3, 1, Dict(:label => 1))
    set_prop!(A, 1, :label, 6)
    set_prop!(A, 2, :label, 6)
    set_prop!(A, 3, :label, 6)

    B = deepcopy(A)
    set_prop!(B, 3, 1, :label, 2)
    
    @test adj_mat(A, B, :direct) == sparse(
        [5, 6, 4, 6, 4, 5, 2, 3, 8, 9, 1, 3, 7, 9, 1, 2, 7, 8, 5, 6, 4, 6, 4, 5], 
        [1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 8, 8, 9, 9], 
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        9, 9
    )

    C = MetaGraph(2)
    add_edge!(C, 1, 2, Dict(:label => 1))
    set_prop!(C, 1, :label, 1)
    set_prop!(C, 2, :label, 1)
    h2 = MetaGraph(smilestomol("[H]-[H]"))

    x = adj_mat(C, C, :direct)

    @test x == adj_mat(h2, h2, :direct)
    @test x == adj_mat(h2, C, :direct)
    @test x == adj_mat(C, h2, :direct)
end

@testset "direct product graph" begin
#    # TODO it seems these are repeated.
#    A = MetaGraph(3)
#    add_edge!(A, 1, 2, Dict(:label => 1))
#    add_edge!(A, 2, 3, Dict(:label => 1))
#    add_edge!(A, 3, 1, Dict(:label => 1))
#    set_prop!(A, 1, :label, 6)
#    set_prop!(A, 2, :label, 6)
#    set_prop!(A, 3, :label, 6)
#
#    B = deepcopy(A)
#    set_prop!(B, 3, 1, :label, 2)
#
#   # @test adjacency_matrix(direct_product_graph(A, B)) == dpg_adj_mat(A, B)
#
#    C = MetaGraph(2)
#    add_edge!(C, 1, 2, Dict(:label => 1))
#    set_prop!(C, 1, :label, 1)
#    set_prop!(C, 2, :label, 1)
#    h2 = smilestomol("[H]-[H]")
#
#    @test direct_product_graph(C, C) == direct_product_graph(MetaGraph(h2), C)
#    @test direct_product_graph(h2, h2) == direct_product_graph(h2, C)
#    @test direct_product_graph(h2, C) == direct_product_graph(C, h2)
   
    # Fig. 4 of "Classifying the toxicity of pesticides to honey bees viasupport vector machines with random walk graph kernels"
    g₁ = MetaGraph(smilestomol("COP(=O)(OC)OC(Br)C(Cl)(Cl)Br"))
    g₂ = MetaGraph(smilestomol("COP(N)(=O)SC"))

    g₁xg₂ = MetaGraph(17)
    set_prop!(g₁xg₂, 1, :label, 15)
    set_prop!(g₁xg₂, 2, :label, 8)
    set_prop!(g₁xg₂, 3, :label, 8)
    set_prop!(g₁xg₂, 4, :label, 8)
    set_prop!(g₁xg₂, 5, :label, 8)
    set_prop!(g₁xg₂, 6, :label, 8)
    set_prop!(g₁xg₂, 7, :label, 8)
    set_prop!(g₁xg₂, 8, :label, 8)
    set_prop!(g₁xg₂, 9, :label, 8)
    set_prop!(g₁xg₂, 10, :label, 6)
    set_prop!(g₁xg₂, 11, :label, 6)
    set_prop!(g₁xg₂, 12, :label, 6)
    set_prop!(g₁xg₂, 13, :label, 6)
    set_prop!(g₁xg₂, 14, :label, 6)
    set_prop!(g₁xg₂, 15, :label, 6)
    set_prop!(g₁xg₂, 16, :label, 6)
    set_prop!(g₁xg₂, 17, :label, 6)
    add_edge!(g₁xg₂, 1, 2, Dict(:label => 2))
    add_edge!(g₁xg₂, 1, 3, Dict(:label => 1))
    add_edge!(g₁xg₂, 1, 4, Dict(:label => 1))
    add_edge!(g₁xg₂, 1, 5, Dict(:label => 1))
    add_edge!(g₁xg₂, 3, 10, Dict(:label => 1))
    add_edge!(g₁xg₂, 4, 11, Dict(:label => 1))
    add_edge!(g₁xg₂, 5, 12, Dict(:label => 1))
    
    computed_g₁xg₂ = product_graph(g₁, g₂, :direct)
    @test is_isomorphic(computed_g₁xg₂, g₁xg₂)
    @test adjacency_matrix(computed_g₁xg₂) == adj_mat(g₁, g₂, :direct)
end

 @testset "SMILES flexibility and symmetry" begin
    A = smilestomol("c1c(C)cccc1")
    B = smilestomol("C1C=C(C)C=CC=1")
    C = smilestomol("C1=C(C)C=CC=C1")
    D = smilestomol("c1-c-c(C)-c-c-c1")
    E = smilestomol("c1ccc(C)cc1")
    A, B, C, D, E = MetaGraph.([A, B, C, D, E])

    for (x, y) in with_replacement_combinations([A, B, C, D, E], 2)
        for type in [:factor, :direct]
            @test is_isomorphic(product_graph(x, x, type), product_graph(x, y, type))
        end
    end
end

@testset "factor product graph" begin
    # Fig 1 of https://arxiv.org/ftp/arxiv/papers/1206/1206.6483.pdf
    g₁ = MetaGraph(smilestomol("NC=O"))
    g₂ = MetaGraph(smilestomol("CN(C=O)C=O"))

    g₁xg₂ = MetaGraph(6)
    set_props!(g₁xg₂, 1, Dict(:label => 7, :v₁v₂_pair => (1, 3)))
    set_props!(g₁xg₂, 2, Dict(:label => 6, :v₁v₂_pair => (2, 2)))
    set_props!(g₁xg₂, 3, Dict(:label => 6, :v₁v₂_pair => (2, 4)))
    set_props!(g₁xg₂, 4, Dict(:label => 6, :v₁v₂_pair => (2, 5)))
    set_props!(g₁xg₂, 5, Dict(:label => 8, :v₁v₂_pair => (3, 1)))
    set_props!(g₁xg₂, 6, Dict(:label => 8, :v₁v₂_pair => (3, 6)))
    add_edge!(g₁xg₂, 1, 2, Dict(:label => 1))
    add_edge!(g₁xg₂, 1, 3, Dict(:label => 1))
    add_edge!(g₁xg₂, 1, 4, Dict(:label => 1))
    add_edge!(g₁xg₂, 1, 5, Dict(:label => 0)) # non-adj
    add_edge!(g₁xg₂, 1, 6, Dict(:label => 0)) # non-adj
    add_edge!(g₁xg₂, 5, 2, Dict(:label => 2))
    add_edge!(g₁xg₂, 6, 4, Dict(:label => 2))

    # TODO assert v₁v₂_pair labels the same too.
    @test is_isomorphic(product_graph(g₁, g₂, :factor), g₁xg₂; edge_labels=[:label])

    rem_edge!(g₁xg₂, 1, 5)
    rem_edge!(g₁xg₂, 1, 6)
    @test is_isomorphic(product_graph(g₁, g₂, :direct), g₁xg₂; edge_labels=[:label])

end

#@testset "csi_adj_mat" begin
#    A = smilestomol("NC=O")
#    B = smilestomol("CN(C=O)C=O")
#    C = MetaGraph(6)
#    set_props!(C, 1, Dict(:label => 7, :v₁v₂_pair => (1, 3)))
#    set_props!(C, 2, Dict(:label => 6, :v₁v₂_pair => (2, 2)))
#    set_props!(C, 3, Dict(:label => 6, :v₁v₂_pair => (2, 4)))
#    set_props!(C, 4, Dict(:label => 6, :v₁v₂_pair => (2, 5)))
#    set_props!(C, 5, Dict(:label => 8, :v₁v₂_pair => (3, 1)))
#    set_props!(C, 6, Dict(:label => 8, :v₁v₂_pair => (3, 6)))
#    add_edge!(C, 1, 2, Dict(:label => 1, :d_type => false))
#    add_edge!(C, 1, 3, Dict(:label => 1, :d_type => false))
#    add_edge!(C, 1, 4, Dict(:label => 1, :d_type => false))
#    add_edge!(C, 1, 5, Dict(:label => 1, :d_type => true))
#    add_edge!(C, 1, 6, Dict(:label => 1, :d_type => true))
#    add_edge!(C, 5, 2, Dict(:label => 2, :d_type => false))
#    add_edge!(C, 6, 4, Dict(:label => 2, :d_type => false))
#    @test is_isomorphic(SimpleGraph(csi_adj_mat(A, B)), C; edge_labels=Symbol[], node_labels=Symbol[])
#    @test is_isomorphic(SimpleGraph(csi_adj_mat(MetaGraph(A), B)), C; edge_labels=Symbol[], node_labels=Symbol[])
#    @test is_isomorphic(SimpleGraph(csi_adj_mat(A, MetaGraph(B))), C; edge_labels=Symbol[], node_labels=Symbol[])
#end
