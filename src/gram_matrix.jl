function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    kwargs...
)::Matrix{Float64}
    sp = sortperm(degree.(molecules); rev=true)
    pairs = reduce(
        vcat,
        [[(i, j, molecules[i], molecules[j]) for j in sp if j ≥ i] for i in sp]
    )
    f_(x) = x[1], x[2], kernel(x[3], x[4]; kwargs...)
    k = @showprogress pmap(f_, pairs)
    matrix = zeros(length(molecules), length(molecules))
    for (i, j, k) in k
        @inbounds matrix[i, j] = matrix[j, i] = k
    end

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
