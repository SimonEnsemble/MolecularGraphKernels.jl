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

"""
    kernel_score = subgraph_matching(adj_mat)
    kernel_score = subgraph_matching(g₁xg₂)
    kernel_score = subgraph_matching(g₁, g₂; λ=_->1)

Returns the similarity score for two graphs by applying the subgraph matching kernel (SM) via their modular product graph.
Currently, the node-pair and edge-pair kernel functions are set as Dirac δ on the node/edge labels.
The default λ assigns a weight of 1 to every isomorphism.
"""
function subgraph_matching(adj_mat::AbstractMatrix; kwargs...)::Int
    return _subgraph_matching(SimpleGraph(adj_mat); kwargs...)
end

function subgraph_matching(g₁xg₂::ProductGraph{Modular}; kwargs...)::Int
    return _subgraph_matching(g₁xg₂.graph; kwargs...)
end

function subgraph_matching(g₁::AbstractMetaGraph, g₂::AbstractMetaGraph; kwargs...)::Int
    return subgraph_matching(product_graph_adjacency_matrix(Modular, g₁, g₂); kwargs...)
end

function subgraph_matching(
    g₁::GraphMol,
    g₂::Union{GraphMol, AbstractMetaGraph};
    kwargs...
)::Int
    return subgraph_matching(MetaGraph(g₁), g₂; kwargs...)
end

function subgraph_matching(g₁::AbstractMetaGraph, g₂::GraphMol; kwargs...)::Int
    return subgraph_matching(g₁, MetaGraph(g₂); kwargs...)
end

# internal function -- enforces requirement that graph passed to exported function be a modular product graph
function _subgraph_matching(g::SimpleGraph; λ::Function=_ -> 1)::Int
    return sum(λ.(maximal_cliques(g)))
end
