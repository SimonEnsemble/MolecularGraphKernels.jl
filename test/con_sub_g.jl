module Test_con_sub_g

using Combinatorics, Graphs, MetaGraphs, MolecularGraph, MolecularGraphKernels, Test 
import MolecularGraphKernels:
    Tree,
    combination_tree,
    ⊗ₜ,
    k_combinations,
    k_compositions,
    con_sub_g,
    combinations_from_tree,
    combinations_with_v,
    add_node!

function isomorphic_trees(t1::Tree, t2::Tree)::Bool
    if !(length(t1.nodes) == length(t2.nodes))
        return false
    end
    for i in 1:length(t1.nodes)
        n₁ = t1[i]
        n₂ = t2[i]
        if n₁.node != n₂.node ||
            n₁.vertex != n₂.vertex ||
            sort([n₁.children[v].node for v in 1:length(n₁.children)]) != sort([n₂.children[v].node for v in 1:length(n₂.children)])
            !!
            sort([n₁.children[v].vertex for v in 1:length(n₁.children)]) != sort([n₂.children[v].vertex for v in 1:length(n₂.children)])
            return false
        end
    end
    return true
end

∅ = Vector{Int}[]

@testset "con_sub_g" begin
    tree_manual = Tree(1, 1)
    add_node!(tree_manual, 1, 2)
    add_node!(tree_manual, 2, 3)
    add_node!(tree_manual, 3, 4)
    add_node!(tree_manual, 2, 4)
    add_node!(tree_manual, 1, 4)
    add_node!(tree_manual, 6, 3)
    add_node!(tree_manual, 1, 5)
    tree_manual[5].new = false
    tree_manual[6].new = false
    tree_manual[7].new = false

    G =  MetaGraph(5)
    for v in vertices(G)
        set_prop!(G, v, :visited, false)
    end
    add_edge!(G, 1, 2)
    add_edge!(G, 1, 4)
    add_edge!(G, 1, 5)
    add_edge!(G, 2, 3)
    add_edge!(G, 2, 4)
    add_edge!(G, 3, 4)

    for e in edges(G)
        set_prop!(G, e, :label, 1)
    end

    tree_algo = combination_tree(1, 4, G)

    @test isomorphic_trees(tree_manual, tree_algo)

    @test k_combinations(2, [2, 6, 8]) == [[2, 6], [2, 8], [6, 8]]

    @test k_compositions(3, 4) == [[1, 1, 2], [1, 2, 1], [2, 1, 1]]

    @test ⊗ₜ([[3]], [[5]], tree_algo) == ∅

    @test any(
        i in sort.(⊗ₜ([[2, 6, 7]], ∅, tree_algo)) for
        i in sort.(collect(permutations([2, 6, 7])))
    )

    @test ⊗ₜ(∅, [[2, 6, 7]], tree_algo) == ∅

    @test ⊗ₜ([[2]], [[6, 7]], tree_algo) == ∅

    @test ⊗ₜ([[2, 3], [2, 5]], [[6]], tree_algo) == ∅

    @test all(
        i in sort.(⊗ₜ([[2, 3], [2, 5]], [[8]], tree_algo)) for
        i in sort.([[2, 3, 8], [2, 5, 8]])
    ) && all(
        i in sort.([[2, 3, 8], [2, 5, 8]]) for
        i in sort.(⊗ₜ([[2, 3], [2, 5]], [[8]], tree_algo))
    )

    @test all(
        i in sort.(con_sub_g(4, G)) for
        i in sort.([[4, 2, 1, 5], [2, 3, 1, 5], [3, 4, 1, 5], [3, 2, 1, 4]])
    ) && all(
        i in sort.([[4, 2, 1, 5], [2, 3, 1, 5], [3, 4, 1, 5], [3, 2, 1, 4]]) for
        i in sort.(con_sub_g(4, G))
    )

    @test all(
        i in sort.(⊗ₜ([[6, 7]], [[8]], tree_algo)) for i in sort.([[6, 7, 8]])
    ) && all(i in sort.([[6, 7, 8]]) for i in sort.(⊗ₜ([[6, 7]], [[8]], tree_algo)))

    combos = sort.(combinations_from_tree(tree_algo, 4, 1))
    test_combos = sort.([[4, 3, 2, 1], [3, 2, 8, 1], [5, 2, 8, 1], [7, 6, 8, 1]])
    @test all(i in combos for i in test_combos) && all(i in test_combos for i in combos)

    combos = sort.(combinations_with_v(1, 4, G))
    test_combos = sort.([[4, 2, 1, 5], [2, 3, 1, 5], [3, 4, 1, 5], [3, 2, 1, 4]])
    @test all(i in combos for i in test_combos) && all(i in test_combos for i in combos)
end

end
