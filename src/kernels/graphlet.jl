using Graphs.Experimental: vf2, IsomorphismProblem

"""
    isomorphism_detected = is_isomorphic(A, B)

    compare two graphs for isomorphism using node and edge labels; or, compare the graph topologies of two adjacency matrices
"""
function is_isomorphic(
    A::AbstractGraph,
    B::AbstractGraph;
    edge_labels::Vector{Symbol}=[:label],
    node_labels::Vector{Symbol}=[:label]
)::Bool
    isomorphic = false

    vf2(
        # check graph isomorphism between A and B, comparing node and edge labels
        SimpleGraph(A),
        SimpleGraph(B),
        IsomorphismProblem();
        vertex_relation=(v, w) ->
            all([get_prop(A, v, x) == get_prop(B, w, x) for x in node_labels]),
        edge_relation=(j, k) ->
            all([get_prop(A, j, x) == get_prop(B, k, x) for x in edge_labels])
    ) do x
        isomorphic = true
        return false
    end

    return isomorphic
end

function is_isomorphic(A::SimpleGraph, B::SimpleGraph)
    return is_isomorphic(A, B; edge_labels=Symbol[], node_labels=Symbol[])
end

function is_isomorphic(A::AbstractMatrix, B::AbstractMatrix)
    return is_isomorphic(SimpleGraph(A), SimpleGraph(B))
end

is_isomorphic(A::AbstractMatrix, B::ProductGraph) = is_isomorphic(A, adjacency_matrix(B))

function is_isomorphic(A::ProductGraph, B::Union{ProductGraph, AbstractMatrix})
    return is_isomorphic(adjacency_matrix(A), B)
end

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
                if depth ≤ k - 1
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
function ⊗ₜ(S₁,S₂,tree)
    unionproduct = Vector{Int}[]
    if S₁ == ∅
        return unionproduct
    end
    if S₂ == ∅
        return S₁
    else
        for s₁ ∈ S₁
            for s₂ ∈ S₂
                no_union = false
                new_ct = any([tree[v].new for v in s₂])
                for i ∈ s₁
                    vᵢ=tree[i].vertex
                    childrenᵢ = [tree[i].children[v].vertex for v in 1:length(tree[i].children)]
                    for j ∈ s₂
                        vⱼ=tree[j].vertex
                        if vᵢ == vⱼ || (!new_ct && in(childrenᵢ.==vⱼ)(1)==true)
                            no_union = true
                        end
                    end
                end
                if !no_union
                    unionproduct = unionproduct ∪ [s₁ ∪ s₂]
                end
            end
        end
    return unionproduct
    end
end

"""
generate the size-`k` node combinations from `tree` starting at `st_root`
"""
@memoize function combinations_from_tree(tree, k::Int, stRoot::Int=1)::Vector{Vector{Int}}
    t=stRoot
    lnodesets = Vector{Int}[]
    k == 1 && return [[t]]
    childrenₜ = [tree[t].children[v].node for v in 1:length(tree[t].children)]
    for i = 1:minimum([length(childrenₜ), k-1])
        for nodecombo in k_combinations(i, childrenₜ)
            for string in k_compositions(i, k - 1)
                S = Dict{Int, Any}()
                fail = false
                for pos in 1:i
                    stRoot = nodecombo[pos]
                    size = string[pos]
                    S[pos] = combinations_from_tree(tree,size,stRoot)
                    if S[pos] == []
                        fail  = true
                        break
                    end
                end
                fail && continue
                for comproduct in reduce((a,b)->⊗ₜ(a,b,tree), [S[i] for i in 1:length(S)])
                    lnodesets = lnodesets ∪ [comproduct ∪ [t]]
                end
            end
        end
    end
    return lnodesets
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
@memoize function con_sub_g(k::Int,graph::MetaGraph)::Vector{Vector{Int}}
    G = deepcopy(graph)
    for v in vertices(G)
        set_prop!(G, v, :visited, false)
    end
    list = Vector{Int}[]
    queue = reverse(vertices(G))
    for v in queue
        list = list ∪ combinations_with_v(v,k,G)
        rem_vertex!(G,v)
    end
    return list
end

"""
calculate the connected graphlet kernel for two graphs G₁ and G₂ using graphlet sizes from `n`
"""
function connected_graphlet(G₁::MetaGraph, G₂::MetaGraph; n=2:4)::Int
    c = 0
    for k in (length(n) == 1 ? [n] : n)
        c += k * count(
            is_isomorphic(
                induced_subgraph(G₁, cgi)[1],
                induced_subgraph(G₂, cgj)[1]
            )
            for cgi in con_sub_g(k, G₁) for cgj in con_sub_g(k, G₂)
        )
    end
    return c
end

export connected_graphlet
