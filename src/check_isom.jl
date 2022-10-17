"""
    isomorphism_detected = is_isomorphic(A, B)
    isomorphic_topology = is_isomorphic(A, B)

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
