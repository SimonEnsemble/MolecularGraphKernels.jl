using Combinatorics, Graphs, MetaGraphs, MolecularGraphKernels, SparseArrays, Test

@testset "direct product graph" begin
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
end

@testset "factor product graph" begin
    # Test case: Fig 1 of https://arxiv.org/ftp/arxiv/papers/1206/1206.6483.pdf
    mol₁ = smilestomol("NC=O")
    mol₂ = smilestomol("CN(C=O)C=O")
    g₁ = MetaGraph(mol₁)
    g₂ = MetaGraph(mol₂)
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

    #  # test: calculated result equivalent to manual construction
    @test is_isomorphic(product_graph(g₁, g₂, :factor), g₁xg₂)

    dpg = deepcopy(g₁xg₂)
    rem_edge!(dpg, 1, 5)
    rem_edge!(dpg, 1, 6)

    # test: removing d-type edges from FPG gives DPG
    @test is_isomorphic(product_graph(g₁, g₂, :direct), dpg)

    # test: type signatures (MetaGraph/GraphMol)
    @test is_isomorphic(product_graph(mol₁, mol₂, :factor), g₁xg₂)
    @test is_isomorphic(product_graph(g₁, mol₂, :factor), g₁xg₂)
    @test is_isomorphic(product_graph(mol₁, g₂, :factor), g₁xg₂)

    ## TODO assert v₁v₂_pair labels the same too.
    ##! Labels aren't the same, b/c numbering in original graphs arbitrary and non-identical
end

@testset verbose=true "product graph adjacency matrices" begin
    mol₁ = smilestomol("CNCCC(c1ccccc1)Oc2ccc(cc2)C(F)(F)F") # fluoxetine
    mol₂ = smilestomol("c1cc(ccc1C(=O)CCCN2CCC(CC2)(c3ccc(cc3)Cl)O)F") # haloperidol
    g₁ = MetaGraph(mol₁)
    g₂ = MetaGraph(mol₂)
    
    for type in [:direct, :factor]
        @testset "$type product graph" begin
            g₁xg₂ = product_graph(g₁, g₂, type)
            # test: adjacency matrix gives same topology as explicit graph
            @test is_isomorphic(g₁xg₂, SimpleGraph(product_graph_adj_mat(g₁, g₂, type)); edge_labels=Symbol[], node_labels=Symbol[])
            # test: type signatures (MetaGraph/GraphMol)
            @test is_isomorphic(g₁xg₂, SimpleGraph(product_graph_adj_mat(g₁, mol₂, type)); edge_labels=Symbol[], node_labels=Symbol[])
            @test is_isomorphic(g₁xg₂, SimpleGraph(product_graph_adj_mat(mol₁, g₂, type)); edge_labels=Symbol[], node_labels=Symbol[])
            @test is_isomorphic(g₁xg₂, SimpleGraph(product_graph_adj_mat(mol₁, mol₂, type)); edge_labels=Symbol[], node_labels=Symbol[])
        end
    end
 end

@testset verbose=true "SMILES flexibility and symmetry" begin
    A = smilestomol("c1c(C)cccc1")
    B = smilestomol("C1C=C(C)C=CC=1")
    C = smilestomol("C1=C(C)C=CC=C1")
    D = smilestomol("c1-c-c(C)-c-c-c1")
    E = smilestomol("c1ccc(C)cc1")
    A, B, C, D, E = MetaGraph.([A, B, C, D, E])

    for type in [:factor, :direct]
        @testset "$type product graph" begin
            for (x, y) in with_replacement_combinations([A, B, C, D, E], 2)
                @test is_isomorphic(product_graph(x, x, type), product_graph(x, y, type))
            end
        end
    end
end
