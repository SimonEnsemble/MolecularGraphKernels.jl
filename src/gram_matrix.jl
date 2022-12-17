function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    local_cache::String=".gram_matrix_local_cache",
    kwargs...
)::Matrix{Float64}
    # determine list of jobs
    ordered_idx = sortperm(nv.(molecules); rev=true)
    jobs = [
        (i, j, molecules[i], molecules[j]) for j in ordered_idx for
        i in ordered_idx if j ≥ i
    ]
    # check for local cache
    if isfile(local_cache)
        @warn "Using local cache"
        # check for completed jobs in local cache file
        skip_jobs = falses(length(jobs))
        for job in read_local_cache(local_cache)
            # if found, add to list of jobs to skip
            job_idx = findfirst(j -> j[1] == job[1] && j[2] == job[2], jobs)
            skip_jobs[job_idx] = true
        end
        # filter job list for only un-completed jobs
        jobs = jobs[.!skip_jobs]
    end
    # the function to run on each job
    f_(x) = x[1], x[2], kernel(x[3], x[4]; kwargs...)
    # results matrix
    matrix = fill(NaN, length(molecules), length(molecules))
    # progress meter (but w/o macro, b/c it can squash errors)
    progress = Progress(length(matrix), 1, "Calculating Gram matrix ")
    # verify that the matrix is filled completely
    nans = count(isnan.(matrix))
    while nans > 0
        # dispatch jobs to workers
        @sync @distributed for job in jobs
            # calculate kernel score
            (i, j, k) = f_(job)
            # write to local cache
            open(local_cache, "a") do f
                return write(f, "$i,$j,$k\n")
            end
        end
        # collect results from local cache
        for (i, j, k) in read_local_cache(local_cache)
            matrix[i, j] = matrix[j, i] = k
        end
        # see if any have been missed
        nans = count(isnan.(matrix))
        # update progress meter
        update!(progress, length(matrix) - nans)
    end
    # (optional) normalize results
    if normalize
        gm_norm!(matrix)
    end
    # clear local cache
    rm(local_cache)
    # return results
    return matrix
end

function read_local_cache(local_cache::String)::Vector{Vector{Int}}
    # read the cache
    return open(local_cache, "r") do f
        # extract the entires from each line in the cache file
        line_tokens = split.(readlines(f), [","])
        # if any are incomplete (corrupted cache) drop them
        filter!(lt -> length(lt) == 3, line_tokens)
        # parse the i/j/k values
        return [parse.(Float64, tokens) for tokens in line_tokens]
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
