abstract type AbstractDisplacementGroup{S,T<:Real} <: AbstractAlgebra.Group end

struct DisplacementGroup{S,T<:Real} <: AbstractDisplacementGroup{S,T}
    d::Int # dimension of euclidean space that the group acts on. E.g R^3 -> d = 3
    dof::Int # degrees of freedom of the group. E.g SE(3) -> dof = 6
    defining_objects::Vector{}
    # gens::Vector{DisplacementGroupElem{S,T}} # generators of the group as DisplacementGroupElem
end

abstract type AbstractDisplacementGroupElem{S,T<:Real} <: AbstractAlgebra.GroupElem end

struct DisplacementGroupElem{S,T<:Real} <: AbstractDisplacementGroupElem{S,T}
    parent::DisplacementGroup{S,T}
    mat::Matrix{T}

    function DisplacementGroupElem(group::DisplacementGroup{S,T}, mat::Matrix{T}) where {S,T<:Real}
        d = group.d
        @assert mat[d+1, d+1] == one(T) && det(mat[1:d, 1:d]) ≈ one(T) "mat must be in SE($d)."

        m = SMatrix{group.d + 1,group.d + 1}(mat)
        new{S,T}(group, m)
    end
end

Base.parent(g::DisplacementGroupElem) = g.parent

(G::DisplacementGroup)(mat::Matrix{<:Real}) = DisplacementGroupElem(G, mat)