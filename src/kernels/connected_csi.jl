"""
computes the Subgraph Matching kernel on product graph Gₚ with weight function λ
"""
function ccsi(Gₚ::ProductGraph{Modular, Float64}; λ::Function=_->1)::Int
    # BKC algorithm
    function bkc_core(C::Vector{Int}, P::U, D::U, S::U) where U <: BitVector
		value += λ(C)
		if !any(P) && !any(S)
			return
		else
			N = falses(n)
			P′ = falses(n)
			D′ = falses(n)
			S′ = falses(n)
			for uᵢ in findall(P)
				P[uᵢ] = false
				P′ .= P
				D′ .= D
				S′ .= S
				N .= false
				for v in neighbors(Gₚ, uᵢ)
					N[v] = true
				end
				for v in findall(D′)
					if has_edge(Gₚ, uᵢ, v) && get_prop(Gₚ, uᵢ, v, :label) ≠ 0
						if T[v]
							S′[v] = true
						else
							P′[v] = true
						end
						D′[v] = false
					end
				end
				bkc_core(C ∪ [uᵢ], P′ .& N, D′ .& N, S′ .& N)
				S[uᵢ] = true
			end
		end
	end

    # initialize and run algorithm
    value = 0
	n = nv(Gₚ)
	T = falses(n)
	P = falses(n)
	D = falses(n)
	S = falses(n)
	for u ∈ vertices(Gₚ)
		P .= false
		D .= false
		S .= false
		N = neighbors(Gₚ, u)
		for v ∈ N
			if get_prop(Gₚ, u, v, :label) ≠ 0
				if T[v]
					S[v] = true
				else
					P[v] = true
				end
			else
				if get_prop(Gₚ, u, v, :label) == 0
					D[v] = true
				end
			end
		end
		bkc_core([u], P, D, S)
		T[u] = true
	end
	return value
end

function ccsi(A::AbstractMetaGraph, B::AbstractMetaGraph; kwargs...)::Int
    return ccsi(ProductGraph{Modular}(A, B); kwargs...)
end

function ccsi(A::GraphMol, B::AbstractMetaGraph; kwargs...)::Int
    return ccsi(MetaGraph(A), B; kwargs...)
end

function ccsi(A::Union{AbstractMetaGraph, GraphMol}, B::GraphMol; kwargs...)::Int
    return ccsi(A, MetaGraph(B); kwargs...)
end
