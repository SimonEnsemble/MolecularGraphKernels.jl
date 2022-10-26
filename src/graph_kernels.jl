"""
    kernel_score = random_walk(adj_mat; l=n)
    kernel_score = random_walk(g₁xg₂; l=n)
    kernel_score = random_walk(g₁, g₂; l=n)

Returns the similarity score for two graphs by applying the `l`-length random walk graph kernel (RWK) on their direct product graph.
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
    kernel_score = subgraph_matching(g₁, g₂)

Returns the similarity score for two graphs by applying the common subgraph isomorphism kernel (CSI) on their modular product graph.
"""
function subgraph_matching(adj_mat::AbstractMatrix)::Int
    return _subgraph_matching(SimpleGraph(adj_mat))
end

function subgraph_matching(g₁xg₂::ProductGraph{Modular})::Int
    return _subgraph_matching(g₁xg₂.graph)
end

function subgraph_matching(g₁::AbstractMetaGraph, g₂::AbstractMetaGraph)::Int
    return subgraph_matching(product_graph_adjacency_matrix(Modular, g₁, g₂))
end

function subgraph_matching(
    g₁::GraphMol,
    g₂::Union{GraphMol, AbstractMetaGraph}
)::Int
    return subgraph_matching(MetaGraph(g₁), g₂)
end

function subgraph_matching(g₁::AbstractMetaGraph, g₂::GraphMol)::Int
    return subgraph_matching(g₁, MetaGraph(g₂))
end

# internal function -- enforces requirement that graph passed to exported function be a modular product graph
function _subgraph_matching(g::SimpleGraph)::Int
    return length(maximal_cliques(g))
end
