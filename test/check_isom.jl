using Graphs.Experimental: vf2, IsomorphismProblem

function is_isomorphic(A::AbstractGraph, B::AbstractGraph; edge_labels::Vector{Symbol}=[:label], node_labels::Vector{Symbol}=[:label])::Bool
    isomorphic = false

    vf2(
        # check graph isomorphism between A and B, comparing node and edge labels
        SimpleGraph(A), SimpleGraph(B), IsomorphismProblem();
        vertex_relation = (v, w) -> all([get_prop(A, v, x) == get_prop(B, w, x) for x in node_labels]),
        edge_relation   = (j, k) -> all([get_prop(A, j, x)  == get_prop(B, k, x) for x in edge_labels])
    ) do x
        isomorphic = true
        return false
    end

    return isomorphic
end
