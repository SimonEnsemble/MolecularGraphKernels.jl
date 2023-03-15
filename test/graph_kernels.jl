module Test_graph_kernels

using MolecularGraph, MolecularGraphKernels, Test, Graphs, MetaGraphs
import MolecularGraphKernels:
    Node,
    Tree,
    ⊗ₜ,
    k_combinations,
    k_compositions,
    con_sub_g,
    combinations_from_tree,
    combinations_with_v

@testset verbose = true "Graph Kernels" begin
    @testset "Random Walk" begin
        g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
        l = 4
        x = random_walk(GraphMatrix{Direct}(g₁, g₂); l=l)
        @test x == 74
        @test x == random_walk(ProductGraph{Direct}(g₁, g₂); l=l)
        @test x == random_walk(g₁, g₂; l=l)

        g₁ = MetaGraph(smilestomol("COP(=O)(OC)OC(Br)C(Cl)(Cl)Br"))
        g₂ = MetaGraph(smilestomol("COP(N)(=O)SC"))
        dpg = ProductGraph{Direct}(g₁, g₂)
        @test random_walk(dpg; l=0) == 17
        @test random_walk(dpg; l=1) == 14
        @test random_walk(dpg; l=2) == 32
    end

    @testset "ethanol-acetic acid test (EZ to count)" begin
        # the two input molecular graphs
        g₁ = MetaGraph(smilestomol("CC(O)=O")) # acetic acid
        g₂ = MetaGraph(smilestomol("CCO")) # ethanol

        ## random walks
        # l = 1. C <--> C * 4 + O <--> O * 2
        @test random_walk(g₁, g₂; l=1) == 6
        # l = 2. C - C <--> C - C * 4 + C - O <--> C - O * 2 + O - C <--> O - C * 2
        @test random_walk(g₁, g₂; l=2) == 8

        ## matching connected subgraphs
        # C ↔ C x4
        # O ↔ O x2
        # C-C ↔ C-C x2
        # C-O ↔ C-O x1
        # C-C-O ↔ C-C-O x1
        @test ccsi(g₁, g₂) == 10
        @test ccsi(g₁, g₂; λ=sum) == 15
    end

    @testset "Connected Common Subgraph Isomorphism" begin
        g₁, g₂ = smilestomol.(["NC=O", "CN(C=O)C=O"])
        mpg = ProductGraph{Modular}(g₁, g₂)
        @test ccsi(mpg) == 13
        @test ccsi(mpg; λ=sum) == 22
    end

    @testset "All-Connected Graphlet Kernel" begin
        G = begin
            local graph = MetaGraph(5)
            for v in vertices(graph)
                set_prop!(graph, v, :visited, false)
            end
            add_edge!(graph, 1, 2)
            add_edge!(graph, 1, 4)
            add_edge!(graph, 1, 5)
            add_edge!(graph, 2, 3)
            add_edge!(graph, 2, 4)
            add_edge!(graph, 3, 4)

            for e in edges(graph)
                set_prop!(graph, e, :label, 1)
            end
            graph
        end

        @test connected_graphlet(ProductGraph{Direct}(G); n=2:4) == 46
    end
end

end
