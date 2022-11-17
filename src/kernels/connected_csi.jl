"""
computes the Subgraph Matching kernel on product graph Gₚ with weight function λ
"""
function ccsi(M::AbstractMatrix; λ::Function=_ -> 1)::Int
    function bkc_core(C::Vector{Int}, P::U, D::U, S::U) where {U <: BitVector}
        value += λ(C)
        if !any(P) && !any(S)
            return
        else
            N = falses(n)
            P′ = falses(n)
            D′ = falses(n)
            S′ = falses(n)
            for uᵢ in findall(P)
                @inbounds P[uᵢ] = false
                P′ .= P
                D′ .= D
                S′ .= S
                N .= false
                @inbounds for v in findall(M[:, uᵢ] .≠ 0)
                    N[v] = true
                end
                for v in findall(D′)
                    @inbounds if M[uᵢ, v] ≠ 0 && M[uᵢ, v] ≠ D_EDGE
                        @inbounds if T[v]
                            @inbounds S′[v] = true
                        else
                            @inbounds P′[v] = true
                        end
                        @inbounds D′[v] = false
                    end
                end
                bkc_core(C ∪ [uᵢ], P′ .& N, D′ .& N, S′ .& N)
                @inbounds S[uᵢ] = true
            end
        end
    end

    value = 0
    n = size(M, 1)
    T = falses(n)
    P = falses(n)
    D = falses(n)
    S = falses(n)
    Cₒ = [0]
    for u in 1:n
        P .= false
        D .= false
        S .= false
        @inbounds for v in findall(M[:, u] .≠ 0)
            @inbounds if M[u, v] ≠ D_EDGE
                @inbounds if T[v]
                    @inbounds S[v] = true
                else
                    @inbounds P[v] = true
                end
            else
                @inbounds D[v] = true
            end
        end
        @inbounds Cₒ[1] = u
        bkc_core(Cₒ, P, D, S)
        @inbounds T[u] = true
    end
    return value
end

ccsi(Gₚ::ProductGraph; kwargs...) = ccsi(GraphMatrix(Gₚ); kwargs...)

function ccsi(A::AbstractMetaGraph, B::AbstractMetaGraph; kwargs...)::Int
    return ccsi(ProductGraph{Modular}(A, B); kwargs...)
end

function ccsi(A::GraphMol, B::AbstractMetaGraph; kwargs...)::Int
    return ccsi(MetaGraph(A), B; kwargs...)
end

function ccsi(A::Union{AbstractMetaGraph, GraphMol}, B::GraphMol; kwargs...)::Int
    return ccsi(A, MetaGraph(B); kwargs...)
end

export ccsi
