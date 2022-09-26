using Graphs.Experimental: vf2, IsomorphismProblem

function is_isomorphic(A::MetaGraph, B::MetaGraph)::Bool
    isomorphic = false

    vf2(
        # check graph isomorphism between A and B, comparing node and edge labels
        SimpleGraph(A), SimpleGraph(B), IsomorphismProblem();
        vertex_relation = (v, w) -> get_prop(A, v, :label) == get_prop(B, w, :label),
        edge_relation   = (j, k) -> get_prop(A, j, :label)  == get_prop(B, k, :label)
    ) do x
        isomorphic = true
        return false
    end

    return isomorphic
end
