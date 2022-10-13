function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    kwargs...
)::Matrix{Int}
    nb_mols = length(molecules)
    matrix = SharedArray{Int}(nb_mols, nb_mols)

    @sync @distributed for i in 1:nb_mols
        for j in i:nb_mols
            matrix[j, i] = kernel(molecules[i], molecules[j]; kwargs...)
            matrix[i, j] = matrix[j, i]
        end
    end
    
    return matrix
end
