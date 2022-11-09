"""
    kernel_score = random_walk(adj_mat; l=n)
    kernel_score = random_walk(g₁xg₂; l=n)
    kernel_score = random_walk(g₁, g₂; l=n)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel (RWK) via their direct product graph.
"""
function random_walk(adj_mat::AbstractMatrix; l::Int)::Int
    return sum(adj_mat^l)
end

function random_walk(g₁xg₂::ProductGraph{Direct}; kwargs...)::Int
    return random_walk(adjacency_matrix(g₁xg₂); kwargs...)
end

function random_walk(A::AbstractMetaGraph, B::AbstractMetaGraph; kwargs...)::Int
    return random_walk(product_graph_adjacency_matrix(Direct, A, B); kwargs...)
end

function random_walk(A::GraphMol, B::AbstractMetaGraph; kwargs...)::Int
    return random_walk(MetaGraph(A), B; kwargs...)
end

function random_walk(A::Union{AbstractMetaGraph, GraphMol}, B::GraphMol; kwargs...)::Int
    return random_walk(A, MetaGraph(B); kwargs...)
end

"""
computes the Subgraph Matching kernel on product graph Gₚ with weight function λ
"""
function subgraph_matching(
    Gₚ::ProductGraph{T},
    λ::Function;
    c_cliques::Bool=false
)::Int where {T <: Union{Modular, Weighted}}
    # Algorithm: SMKernel(w, C, P)
    # Input: Product graph Gₚ, weight function λ
    # Initial: value ← 0; SMKernel(1, ∅, Vₚ)
    # Param.: Weight w of the clique C, candidate set P
    # Output: Result of the kernel function value

    # initialize
    value = 0
    Vₚ = collect(vertices(Gₚ))

    # define recursive algorithm
    function smkernel(w::Int, C::Vector{Int}, P::Vector{Int})
        while length(P) > 0 # while |P| > 0 do
            @inbounds v = P[1] # v ← arbitrary element of P
            w′ = w
            if !c_cliques || extends_clique(Gₚ, C, v)
                C′ = union(C, v)
                w′ *= prod(smkernel_c(Gₚ, u, v) for u in C; init=smkernel_c(Gₚ, v)) # multiply by vertex, edge weights
                value += w′ * λ(C′)
            else
                C′ = C
            end
            smkernel(w′, C′, intersect(P, neighbors(Gₚ, v))) # extend clique
            @inbounds P = P[2:end] # P ← P \ {v}
        end
        return
    end

    # run algorithm
    smkernel(1, Int[], Vₚ)
    return value
end

# node weight function for SM kernel on modular product graph (i.e. CSI kernel)
smkernel_c(::ProductGraph{Modular}, ::Int)::Int = 1

# edge weight function for SM kernel on modular product graph (i.e. CSI kernel)
smkernel_c(::ProductGraph{Modular}, ::Int, ::Int)::Int = 1

# checks if a node v extends the clique C in the graph G
function extends_clique(G::AbstractGraph, C::Vector{Int}, v::Int)::Bool
	if C == []
		return true
	end
	for u in C
		if has_edge(G, u, v) && get_prop(G, u, v, :label) ≠ 0
			return true
		end
	end
	return false
end

"""
    kernel_score = common_subgraph_isomorphism(g₁xg₂)
    kernel_score = common_subgraph_isomorphism(g₁, g₂; λ=_->1)

Returns the similarity score for two graphs by applying the common subgraph isomorphism kernel (CSI) via their modular product graph.
Currently, the node-pair and edge-pair kernel functions are set as Dirac δ on the node/edge labels.
The default λ assigns a weight of 1 to every isomorphism.
"""
function common_subgraph_isomorphism(g₁xg₂::ProductGraph{Modular}; λ::Function=_ -> 1, kwargs...)::Int
    return subgraph_matching(g₁xg₂, λ; kwargs...)
end

function common_subgraph_isomorphism(
    g₁::AbstractMetaGraph,
    g₂::AbstractMetaGraph;
    kwargs...
)::Int
    return common_subgraph_isomorphism(ProductGraph{Modular}(g₁, g₂); kwargs...)
end

function common_subgraph_isomorphism(
    g₁::GraphMol,
    g₂::Union{GraphMol, AbstractMetaGraph};
    kwargs...
)::Int
    return common_subgraph_isomorphism(MetaGraph(g₁), g₂; kwargs...)
end

function common_subgraph_isomorphism(g₁::AbstractMetaGraph, g₂::GraphMol; kwargs...)::Int
    return common_subgraph_isomorphism(g₁, MetaGraph(g₂); kwargs...)
end
