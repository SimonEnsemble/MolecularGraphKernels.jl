function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    kwargs...
)::Matrix{Float64}
    nb_mols = length(molecules)
    matrix = SharedArray{Float64}(nb_mols, nb_mols)

    sp = sortperm(degree.(molecules), rev=true)

    @sync @distributed for i in eachindex(sp)
        for j in i:length(sp)
            @inbounds matrix[sp[i], sp[j]] = matrix[sp[j], sp[i]] = kernel(molecules[sp[i]], molecules[sp[j]]; kwargs...)
        end
    end

    matrix = Matrix{Float64}(matrix)

    if normalize
        gm_norm!(matrix)
    end

    return matrix
end

function gm_norm!(mat::Matrix{Float64})
    @inbounds sqrts = [√(mat[i, i]) for i in axes(mat, 1)]
    for j in axes(mat, 2)
        for i in axes(mat, 1)
            @inbounds mat[i, j] /= sqrts[i] * sqrts[j]
            @inbounds mat[j, i] = mat[i, j]
        end
    end
end

function gm_norm(mat::Matrix{<:Real})::Matrix{Float64}
    m = deepcopy(mat)
    gm_norm!(m)
    return m
end

export gram_matrix, gm_norm
