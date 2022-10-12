function gram_matrix(kernel::Function, molecules::Vector{T}; kwargs...)::Matrix{Int} where {T <: AbstractMetaGraph}
    nb_mols = length(molecules)
    matrix = zeros(Int, nb_mols, nb_mols)

    @threads for i in 1:nb_mols
        mol = molecules[i]
        matrix[i:nb_mols, i] .= [kernel(mol, molecules[j]; kwargs...) for j in i:nb_mols]
    end
    for i in 1:nb_mols
        @threads for j in i:nb_mols
            matrix[i, j] = matrix[j, i]
        end
    end

    return matrix
end
