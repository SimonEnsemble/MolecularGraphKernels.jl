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
    λ::Function
)::Int where {T <: Union{Modular, Weighted}}
    # Algorithm: SMKernel(w, C, P)
    # Input: Product graph Gₚ, weight function λ
    # Initial: value ← 0; SMKernel(1, ∅, Vₚ)
    # Param.: Weight w of the clique C, candidate set P
    # Output: Result of the kernel function value

    # initialize
    value = 0
    ∅ = Int[]
    Vₚ = collect(vertices(Gₚ))

    # define recursive algorithm
    function smkernel(w::Int, C::Vector{Int}, P::Vector{Int})
        while length(P) > 0 # while |P| > 0 do
            v = first(P) # v ← arbitrary element of P
            C′ = union(C, v)
            w′ = w * smkernel_c(Gₚ, v) # multiply by vertex weight
            for u in C
                w′ *= smkernel_c(Gₚ, u, v)# multiply by edge weights
            end
            value += w′ * λ(C′)
            smkernel(w′, C′, intersect(P, neighbors(Gₚ, v))) # extend clique
            P = setdiff(P, [v]) # P ← P \ {v}
        end
        return
    end

    # run algorithm
    smkernel(1, ∅, Vₚ)
    return value
end

# node weight function for SM kernel on modular product graph (i.e. CSI kernel)
smkernel_c(::ProductGraph{Modular}, ::Int) = 1

# edge weight function for SM kernel on modular product graph (i.e. CSI kernel)
smkernel_c(::ProductGraph{Modular}, ::Int, ::Int) = 1

"""
    kernel_score = common_subgraph_isomorphism(g₁xg₂)
    kernel_score = common_subgraph_isomorphism(g₁, g₂; λ=_->1)

Returns the similarity score for two graphs by applying the common subgraph isomorphism kernel (CSI) via their modular product graph.
Currently, the node-pair and edge-pair kernel functions are set as Dirac δ on the node/edge labels.
The default λ assigns a weight of 1 to every isomorphism.
"""
function common_subgraph_isomorphism(g₁xg₂::ProductGraph{Modular}; λ::Function=_ -> 1)::Int
    return subgraph_matching(g₁xg₂, λ)
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
