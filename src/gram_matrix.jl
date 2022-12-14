function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    resume_from_error::String="",
    kwargs...
)::Matrix{Float64}
    if resume_from_error ≠ ""
        matrix = load(resume_from_error)["matrix"]
        skip_idx = findall(matrix .≠ -1)
    else
        matrix = -1 * ones(length(molecules), length(molecules))
        skip_idx = []
    end
    sp = sortperm(nv.(molecules); rev=true)
    tasks = [(i, j, molecules[i], molecules[j]) for j in sp for i in sp if j ≥ i]
    filter!(t -> CartesianIndex(t[1], t[2]) ∉ skip_idx, tasks)
    f_(x) = x[1], x[2], kernel(x[3], x[4]; kwargs...)
    k = @showprogress pmap(f_, tasks)
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
        error("Matrix element(s) not calculated. Something has gone wrong!")
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
