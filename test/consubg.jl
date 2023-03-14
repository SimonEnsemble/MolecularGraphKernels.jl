
module Test_ConSubG

using MolecularGraph, MolecularGraphKernels, Test, Graphs, MetaGraphs, Combinatorics
import MolecularGraphKernels: Node, Tree, combination_tree, ⊗ₜ, kCombinations, kCompositions, ConSubG, CombinationsFromTree, CombinationsWithV, add_node!, Graphs, MetaGraphs, Combinatorics, isomorphic_trees

    @testset "ConSubG" begin
        T_manual = begin
            tree = Tree(1,1)
            add_node!(tree,1,2)
            add_node!(tree,2,3)
            add_node!(tree,3,4)
            add_node!(tree,2,4)
            add_node!(tree,1,4)
            add_node!(tree,6,3)
            add_node!(tree,1,5)
            tree[5].new = false
            tree[6].new = false
            tree[7].new = false
            tree
        end

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

        T_from_algorithm = combination_tree(1, 4, G)

        @test isomorphic_trees(T_manual, T_from_algorithm)

        @test kCombinations(2, [2, 6, 8]) == [[2, 6], [2, 8], [6, 8]]

        @test kCompositions(3, 4) == [[1, 1, 2], [1, 2, 1], [2, 1, 1]]

        @test ⊗ₜ([[3]],[[5]],T_from_algorithm) == [[]]

        @test any([i in ⊗ₜ([[2,6,7]],[[]],T_from_algorithm) for i in collect(permutations([2,6,7]))])

        @test ⊗ₜ([[]],[[2,6,7]],T_from_algorithm) == [[]]

        @test ⊗ₜ([[2]],[[6,7]],T_from_algorithm) == [[]]

        @test ⊗ₜ([[2,3],[2,5]],[[6]],T_from_algorithm) == [[]]

        @test all(i in sort.(⊗ₜ([[2,3],[2,5]],[[8]],T_from_algorithm)) for i in sort.([[2,3,8],[2,5,8]])) && all(i in sort.([[2,3,8],[2,5,8]]) for i in sort.(⊗ₜ([[2,3],[2,5]],[[8]],T_from_algorithm)))

        @test all(i in sort.(ConSubG(4,G)) for i in sort.([[4,2,1,5],[2,3,1,5],[3,4,1,5],[3,2,1,4]])) && all(i in sort.([[4,2,1,5],[2,3,1,5],[3,4,1,5],[3,2,1,4]]) for i in sort.(ConSubG(4,G)))

        @test all(i in sort.(⊗ₜ([[6,7]],[[8]],T_from_algorithm)) for i in sort.([[6,7,8]])) && all(i in sort.([[6,7,8]]) for i in sort.(⊗ₜ([[6,7]],[[8]],T_from_algorithm)))

        @test all(i in sort.(CombinationsFromTree(T_from_algorithm,4,1)) for i in sort.([[4,3,2,1],[3,2,8,1],[5,2,8,1],[7,6,8,1]])) && all(i in sort.([[4,3,2,1],[3,2,8,1],[5,2,8,1],[7,6,8,1]]) for i in sort.(CombinationsFromTree(T_from_algorithm,4,1)))

        @test all(i in sort.(CombinationsWithV(1,4,G)) for i in sort.([[4,2,1,5],[2,3,1,5],[3,4,1,5],[3,2,1,4]])) && all(i in sort.([[4,2,1,5],[2,3,1,5],[3,4,1,5],[3,2,1,4]]) for i in sort.(CombinationsWithV(1,4,G)))

    end
end