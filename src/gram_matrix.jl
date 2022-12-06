function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    kwargs...
)::Matrix{Float64}
    sp = sortperm(degree.(molecules); rev=true)
    tasks = reduce(
        vcat,
        [[(i, j, molecules[i], molecules[j]) for j in sp if j ≥ i] for i in sp]
    )
    f_(x) = x[1], x[2], kernel(x[3], x[4]; kwargs...)
    k = @showprogress pmap(f_, tasks)
    matrix = -1 * ones(length(molecules), length(molecules))
    for (i, j, k) in k
        matrix[i, j] = matrix[j, i] = k
    end

    dump_on_error(matrix)

    if normalize
        gm_norm!(matrix)
    end

    return matrix
end

function dump_on_error(matrix::Matrix{Float64})
    # running out of RAM can cause values to not be assigned...
    if any(matrix .== -1)
        jldsave("gram_matrix_dump_$(split(tempname(), "/")[end]).jld2"; matrix)
        error("Zero(s) on Gram matrix diagonal. Something has gone wrong!")
    end
end

function gm_norm!(mat::Matrix{Float64})
    sqrts = [√(mat[i, i]) for i in axes(mat, 1)]
    for j in axes(mat, 2)
        for i in [i for i in axes(mat, 1) if i ≥ j]
            mat[i, j] /= sqrts[i] * sqrts[j]
            mat[j, i] = mat[i, j]
        end
    end
end

function gm_norm(mat::Matrix{<:Real})::Matrix{Float64}
    m = deepcopy(mat)
    gm_norm!(m)
    return m
end

export gram_matrix, gm_norm
