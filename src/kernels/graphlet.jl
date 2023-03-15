
mutable struct Node
    node::Int
    vertex::Int
    new::Bool
    children::Vector{Node}

    Node(node::Int, vertex::Int) = new(node, vertex, true, Node[])
end

struct Tree
    root::Int
    nodes::Dict{Int, Node}

    Tree(root::Int, vertex::Int) = new(root, Dict(root => Node(1, vertex)))
end

import Base.getindex
Base.getindex(tree::Tree, idx::Int) = tree.nodes[idx]

function add_node!(tree::Tree, node::Int, vertex::Int)
    # create node
    new_node = Node(length(tree.nodes) + 1, vertex)
    # link to parent
    tree[node].children = vcat(tree[node].children, [new_node])
    # update dictionary
    return tree.nodes[new_node.node] = new_node
end;

function combination_tree(v, k, graph)
    A = adjacency_matrix(graph)
    vertex_labels = vertices(graph)
    visited = [get_prop(graph, v, :visited) for v in vertex_labels]
    neighbor_vertices = [
        ((A[:, v]) .* vertex_labels)[(A[:, v]) .* vertex_labels .> 0] for v in vertex_labels
    ]
    function build_tree(nₜ, depth, k)
        list[:, depth] = list[:, depth - 1]
        for v′ in neighbor_vertices[tree[nₜ].vertex]
            if !list[v′, depth]
                add_node!(tree, nₜ, v′)
                nₜ′ = length(tree.nodes)
                list[v′, depth] = true
                if !visited[v′]
                    visited[v′] = true
                else
                    tree[nₜ′].new = false
                end
                if depth ≤ k
                    build_tree(nₜ′, depth + 1, k)
                end
            end
        end
    end
    tree = Tree(1, v)
    list = zeros(Bool, nv(graph), k + 1)
    list[v, 1] = true
    build_tree(1, 2, k)
    return tree
end

kCombinations(k, S) = collect(combinations(S, k))

@memoize function kCompositions(elements::Int, sum_value::Int)::Vector{Vector{Int}}
    if elements == 0 || sum_value == 0 || elements > sum_value
        return []
    end
    if elements == 1
        return [[sum_value]]
    else
        return unique(
            reduce(
                vcat,
                collect.([
                    permutations(c, elements) for c in filter(
                        c -> sum(c) == sum_value,
                        collect(
                            with_replacement_combinations(
                                1:(sum_value - elements + 1),
                                elements
                            )
                        )
                    )
                ])
            )
        )
    end
end

function ⊗ₜ(S₁, S₂, tree)::Vector{Vector{Int}}
    UnionProduct = Dict()
    if S₁ == [[]] || S₂ == [[]]
        return S₁
    else
        for s₁ in S₁
            for s₂ in S₂
                new_ct = any([tree[v].new for v in s₂])
                for i in s₁
                    vᵢ = tree[i].vertex
                    Childrenᵢ =
                        [tree[i].children[v].vertex for v in 1:length(tree[i].children)]
                    for j in s₂
                        vⱼ = tree[j].vertex
                        if vᵢ == vⱼ || (!new_ct && in(Childrenᵢ .== vⱼ)(1) == true)
                            return [[]]
                        end
                    end
                end
                UnionProduct = UnionProduct ∪ [s₁ ∪ s₂]
            end
        end
    end
    return UnionProduct
end

@memoize function CombinationsFromTree(tree, k::Int, stRoot::Int=1)::Vector{Vector{Int}}
    t = stRoot
    lnodesets = []
    k == 1 && return [[t]]
    Childrenₜ = [tree[t].children[v].node for v in 1:length(tree[t].children)]
    for i in 1:minimum([length(Childrenₜ), k - 1])
        for NodeComb in kCombinations(i, Childrenₜ)
            for string in kCompositions(i, k - 1)
                S = Dict()
                fail = false
                for pos in 1:i
                    stRoot = NodeComb[pos]
                    size = string[pos]
                    S[pos] = CombinationsFromTree(tree, size, stRoot)
                    if S[pos] == []
                        fail = true
                        break
                    end
                end
                fail && continue
                for comProduct in
                    reduce((a, b) -> ⊗ₜ(a, b, tree), [S[i] for i in 1:length(S)])
                    lnodesets = lnodesets ∪ [comProduct ∪ [t]]
                end
            end
        end
    end
    return lnodesets[length.(lnodesets) .== k]
end

function CombinationsWithV(v, k, graph)
    tree = combination_tree(v, k, graph)
    ncombs = CombinationsFromTree(tree, k)
    return [[tree[v].vertex for v in Set] for Set in ncombs]
end

function ConSubG(k, graph)
    G = deepcopy(graph)
    for v in vertices(G)
        set_prop!(G, v, :visited, false)
    end
    list = Dict()
    queue = reverse(vertices(G))
    for v in queue
        list = list ∪ CombinationsWithV(v, k, G)
        rem_vertex!(G, v)
    end
    return list
end

function connected_graphlet(G::ProductGraph; n=2:4)::Int
    if length(n) == 1
        return n * length(ConSubG(n, G))
    else
        return sum([k * length(ConSubG(k, G)) for k in n])
    end
end

function isomorphic_trees(t1, t2)
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

function connected_graphlet(
    A::Union{GraphMol, MetaGraph},
    B::Union{GraphMol, MetaGraph};
    kwargs...
)::Int
    return connected_graphlet(ProductGraph{Direct}(A, B); kwargs...)
end

export connected_graphlet
