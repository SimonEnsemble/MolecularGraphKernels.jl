"""
struct for nodes in the combination tree
"""
mutable struct Node
    node::Int
    vertex::Int
    new::Bool
    children::Vector{Node}

    Node(node::Int, vertex::Int) = new(node, vertex, true, Node[])
end

"""
struct for the combination tree
"""
struct Tree
    root::Int
    nodes::Dict{Int, Node}

    Tree(root::Int, vertex::Int) = new(root, Dict(root => Node(1, vertex)))
end

import Base.getindex
Base.getindex(tree::Tree, idx::Int) = tree.nodes[idx]

"""
add a node to `tree`, as a child of node number `node`,
associated with graph vertex `vertex`
"""
function add_node!(tree::Tree, node::Int, vertex::Int)
    # create node
    new_node = Node(length(tree.nodes) + 1, vertex)
    # link to parent
    tree[node].children = vcat(tree[node].children, [new_node])
    # update dictionary
    tree.nodes[new_node.node] = new_node
    return
end

"""
build the combination tree for size-`k` graphlets in `graph` containing vertex `v`
"""
function combination_tree(v::Int, k::Int, graph::AbstractMetaGraph)::Tree
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

"""
get the set of combinations of size `k` from collection `S`
"""
k_combinations(k::Int, S::Vector{Int}) = collect(combinations(S, k))

"""
get the set of sets of integers of size `elements` which sum to `sum_value`
"""
@memoize function k_compositions(elements::Int, sum_value::Int)::Vector{Vector{Int}}
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

const ∅ = Vector{Int}[]

"""
union product
"""
function ⊗ₜ(S₁::T, S₂::T, tree::Tree)::T where {T <: Vector{Vector{Int}}}
    if S₁ == ∅ || S₂ == ∅
        return S₁
    else
        union_product = Vector{Int}[]
        for s₁ in S₁
            for s₂ in S₂
                new_ct = any([tree[v].new for v in s₂])
                for i in s₁
                    vᵢ = tree[i].vertex
                    children = [child.vertex for child in tree[i].children]
                    for j in s₂
                        vⱼ = tree[j].vertex
                        if vᵢ == vⱼ || (!new_ct && in(children .== vⱼ)(1) == true)
                            return ∅
                        end
                    end
                end
                union_product = union_product ∪ [s₁ ∪ s₂]
            end
        end
        return union_product
    end
end

"""
generate the size-`k` node combinations from `tree` starting at `st_root`
"""
@memoize function combinations_from_tree(
    tree::Tree,
    k::Int,
    st_root::Int=1
)::Vector{Vector{Int}}
    t = st_root
    lnodesets = Vector{Int}[]
    k == 1 && return [[t]]
    children = [tree[t].children[v].node for v in 1:length(tree[t].children)]
    for i in 1:minimum([length(children), k - 1])
        for node_combo in k_combinations(i, children)
            for string in k_compositions(i, k - 1)
                S = Dict()
                fail = false
                for pos in 1:i
                    st_root = node_combo[pos]
                    size = string[pos]
                    S[pos] = combinations_from_tree(tree, size, st_root)
                    if S[pos] == []
                        fail = true
                        break
                    end
                end
                fail && continue
                for combo_product in
                    ## I will give $20 to whoever can explain why I can't reduce 
                    ## over S or [S[i] for i in eachindex(S)].  Not even kidding.
                    reduce((a, b) -> ⊗ₜ(a, b, tree), [S[i] for i in 1:length(S)])
                    lnodesets = lnodesets ∪ Vector{Int}[combo_product ∪ Int[t]]
                end
            end
        end
    end
    return lnodesets[length.(lnodesets) .== k]
end

"""
generate the set of size-`k` sets of `graph` vertices that include vertex `v`
"""
function combinations_with_v(v::Int, k::Int, graph::AbstractMetaGraph)::Vector{Vector{Int}}
    tree = combination_tree(v, k, graph)
    ncombs = combinations_from_tree(tree, k)
    return Vector{Int}[Int[tree[v].vertex for v in set] for set in ncombs]
end

"""
enumerate the set of sets of vertices comprising the size-`k` graphlets of `graph`
"""
function con_sub_g(k::Int, graph::AbstractMetaGraph)::Vector{Vector{Int}}
    G = deepcopy(graph)
    for v in vertices(G)
        set_prop!(G, v, :visited, false)
    end
    list = Vector{Int}[]
    queue = reverse(vertices(G))
    for v in queue
        list = list ∪ combinations_with_v(v, k, G)
        rem_vertex!(G, v)
    end
    return list
end

"""
calculate the connected graphlet kernel for DPG `G` using graphlet sizes from `n`
"""
function connected_graphlet(
    G::ProductGraph{Direct};
    n::Union{Vector{Int}, UnitRange}=2:4
)::Int
    if length(n) == 1
        return n * length(con_sub_g(n, G))
    else
        return sum([k * length(con_sub_g(k, G)) for k in n])
    end
end

"""
calculate the connected graphlet kernel for the DPG of `A` and `B`
"""
function connected_graphlet(
    A::Union{GraphMol, MetaGraph},
    B::Union{GraphMol, MetaGraph};
    kwargs...
)::Int
    return connected_graphlet(ProductGraph{Direct}(A, B); kwargs...)
end

export connected_graphlet
