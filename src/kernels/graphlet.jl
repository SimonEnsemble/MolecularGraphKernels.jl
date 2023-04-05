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
	UnionProduct = Dict()
	if S₁ == ∅
		return UnionProduct
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
					Childrenᵢ = [tree[i].children[v].vertex for v in 1:length(tree[i].children)]
					for j ∈ s₂
						vⱼ=tree[j].vertex
						if vᵢ == vⱼ || (!new_ct && in(Childrenᵢ.==vⱼ)(1)==true)
							no_union = true
						end
					end
				end
				if !no_union
					UnionProduct = UnionProduct ∪ [s₁ ∪ s₂]
				end
			end
		end
	return UnionProduct
	end
end

"""
generate the size-`k` node combinations from `tree` starting at `st_root`
"""
@memoize function CombinationsFromTree(tree,k::Int,stRoot::Int=1)::Vector{Vector{Int}}
	t=stRoot
	lnodesets = []
	k==1 && return [[t]]
	Childrenₜ = [tree[t].children[v].node for v in 1:length(tree[t].children)]
	for i = 1:minimum([length(Childrenₜ),k-1])
		for NodeComb in kCombinations(i,Childrenₜ)
			for string in kCompositions(i,k-1)
				S = Dict()
				fail = false
				for pos in 1:i
					stRoot = NodeComb[pos]
					size = string[pos]
					S[pos] = CombinationsFromTree(tree,size,stRoot)
					if S[pos] == []
						fail  = true
						break
					end
				end
				fail && continue
				for comProduct in reduce((a,b)->⊗ₜ(a,b,tree), [S[i] for i in 1:length(S)])
					lnodesets = lnodesets ∪ [comProduct ∪ [t]]
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
function con_sub_g(k::Int,graph::MetaGraph)::Vector{Vector{Int}}
	G = deepcopy(graph)
	for v in vertices(G)
		set_prop!(G, v, :visited, false)
	end
	list = Dict()
	queue = reverse(vertices(G))
	for v in queue
		list = list ∪ CombinationsWithV(v,k,G)
		rem_vertex!(G,v)
	end
	return list
end

"""
calculate the connected graphlet kernel for two graphs G₁ and G₂ using graphlet sizes from `n`
"""
function connected_graphlet(G₁::MetaGraph,G₂::MetaGraph; n=2:4)::Int
	count = 0
	if length(n) == 1
		cg1 = ConSubG(n, G₁)
		cg2 = ConSubG(n, G₂)
		count = count + sum([
			is_isomorphic(
				induced_subgraph(G₁, cg1[i])[1],
				induced_subgraph(G₂, cg2[j])[1]
				)
			 for i in eachindex(cg1) for j in eachindex(cg2)
				 ])*n
	else
		for k in n
		    cg1 = ConSubG(k, G₁)
		    cg2 = ConSubG(k, G₂)
		    count = count + sum([
		        is_isomorphic(
		            induced_subgraph(G₁, cg1[i])[1],
		            induced_subgraph(G₂, cg2[j])[1]
		            )
		        for i in eachindex(cg1) for j in eachindex(cg2)
					])*k
		end
	end
	return count
end

