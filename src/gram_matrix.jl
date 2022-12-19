function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    local_cache::AbstractString=split(tempname(), "/")[end],
    kwargs...
)::Matrix{Float64}
    # build the job list
    jobs, matrix = enumerate_jobs(molecules, local_cache)

    # generate task channels
    job_channel = RemoteChannel(() -> Channel{Tuple}(length(jobs)))
    result_channel = RemoteChannel(() -> Channel{Tuple}(length(workers())))

    # push jobs to channel
    progress = Progress(length(jobs); desc="Pushing jobs to queue ")
    for job in jobs
        @async errormonitor(put!(job_channel, job))
        next!(progress)
    end

    # worker task function
    function do_work(job_channel, result_channel)
        # loop until jobs channel is empty
        while isready(job_channel)
            # take a job from the job channel
            i, j, g₁, g₂ = take!(job_channel)
            # calculate the result and put it onto the result channel
            put!(result_channel, (i, j, kernel(g₁, g₂; kwargs...)))
        end
    end

    # start running the jobs on worker processes
    for p in workers()
        remote_do(do_work, p, job_channel, result_channel)
    end

    # progress meter
    progress = Progress(length(jobs); desc="Calculating Gram matrix ")
    update!(progress, 0)

    # collect results from workers
    open(local_cache, "a") do cache
        # we need to get back as many results as we sent out jobs
        for _ in jobs
            # get the return values from a job (block until one comes in)
            i, j, k = take!(result_channel)
            # fill in the matrix
            matrix[i, j] = matrix[j, i] = k
            # write to the recovery file
            write(cache, "$i,$j,$k\n")
            # increment progress meter
            next!(progress)
        end
    end

    # verify that the matrix really is complete
    @assert !any(isnan, matrix)

    # (optional) normalize results
    if normalize
        gm_norm!(matrix)
    end

    # clean up local cache (got this far -> don't need it)
    rm(local_cache)

    # return results
    return matrix
end

function enumerate_jobs(molecules, local_cache)
    # determine full list of jobs
    ordered_idx = sortperm(nv.(molecules); rev=true)
    jobs = [
        (i, j, molecules[i], molecules[j]) for j in ordered_idx for
        i in ordered_idx if j ≥ i
    ]
    # results matrix
    matrix = fill(NaN, length(molecules), length(molecules))

    # (optionally) resume from local cache
    if isfile(local_cache)
        @warn "Resuming from local cache"
        # indices of jobs that will be skipped (already completed)
        skip_jobs = trues(length(jobs))
        # read the cache
        cache_data = open(local_cache, "r") do f
            # extract the entires from each line in the cache file
            line_tokens = split.(readlines(f), [","])
            # if any are incomplete (corrupted cache) drop them
            filter!(lt -> length(lt) == 3, line_tokens)
            # parse the i/j/k values
            return [parse.(Float64, tokens) for tokens in line_tokens]
        end
        # compare cache and jobs list
        progress = Progress(length(cache_data); desc="Loading cached results ")
        for job in cache_data
            # record cached result
            i, j = Int.(job[1:2])
            k = job[3]
            matrix[i, j] = matrix[j, i] = k
            next!(progress)
        end
        # find indices of jobs not in cache
        missing_idx = findall(isnan, matrix)
        # filter job list for only un-completed jobs
        for idx in missing_idx
            # find the job corresponding to the missing index
            i = findfirst(job -> job[1] == idx[1] && job[2] == idx[2], jobs)
            # only job (i, j) will exist; just move on if looking for duplicate (j, i)
            if !isnothing(i)
                skip_jobs[i] = false
            end
        end
        jobs = jobs[.!skip_jobs]
    end

    return jobs, matrix
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
