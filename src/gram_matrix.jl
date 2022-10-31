function gram_matrix(
    kernel::Function,
    molecules::Vector{MetaGraph{Int, Float64}};
    normalize::Bool=false,
    kwargs...
)::Matrix{Float64}
    nb_mols = length(molecules)
    matrix = SharedArray{Float64}(nb_mols, nb_mols)

    @sync @distributed for i in 1:nb_mols
        for j in i:nb_mols
            matrix[j, i] = kernel(molecules[i], molecules[j]; kwargs...)
            matrix[i, j] = matrix[j, i]
        end
    end

    matrix = Matrix{Float64}(matrix)

    if normalize
        gm_norm!(matrix)
    end

    return matrix
end

function gm_norm!(mat::Matrix{Float64})
	n = size(mat)[1]
	for j in 1:n
		for i in j:n
			mat[i, j] /= √(mat[i, i] * mat[j, j])
			mat[j, i] = mat[i, j]
		end
	end
end

function gm_norm(mat::Matrix{<: Real})::Matrix{Float64}
	m = deepcopy(mat)
	gm_norm!(m)
	return m
end
