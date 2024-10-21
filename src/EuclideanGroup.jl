
######################## displacement matrices in R^3
# rotation along x-axis
rotation_x(phi::Real) = [1 0 0 0; 0 cos(phi) -sin(phi) 0; 0 sin(phi) cos(phi) 0; 0 0 0 1]

# rotation along arbitrary axis. The axis is represented by a vector and a point
function rotation(w::AbstractVector, q::Point, phi::Real)
    # compute rotation matrix by base change of the rotation matrix along x-axis
    p = q.coords
    v = normalize(w) # directional vector of axis.
    n = nullspace(transpose(v)) # orthogonal directional vectors of plane with normal vector v
    n1 = n[:, 1]
    n2 = n[:, 2]

    EidB = hcat(v, n1, n2, p)
    EidB = vcat(EidB, [0 0 0 1])
    BrhoB = rotation_x(phi::Real)
    BidE = EidB^(-1)

    return EidB * BrhoB * BidE
end

# rotation along an arbitrary axis represented by a line. The vector of axis is axis.points[2] - axis.points[1]
rotation(axis::Line{3}, phi::Real) = rotation(axis.points[2] - axis.points[1], axis.points[1], phi)

#translation along arbitrary axis represented by a vector
function translation(v::AbstractVector, t::Real)
    mat = [1 0 0; 0 1 0; 0 0 1]
    mat = hcat(mat, t * v)
    mat = vcat(mat, [0 0 0 1])
    return mat
end

# translation along arbitrary axis represented as a Line. The direction of axis is axis.points[2] - axis.points[1]
translation(axis::Line, t::Real) = translation(axis.points[2] - axis.points[1], t)


abstract type AbstractDisplacementGroup{S,T<:Real} <: AbstractAlgebra.Group end
struct DisplacementGroup{N,T<:Real} <: AbstractDisplacementGroup{N,T}
    d::Int # dimension of euclidean space that the group acts on. E.g R^3 -> d = 3
    dof::Int # degrees of freedom of the group. E.g SE(3) -> dof = 6
    defining_objects::Vector{<:AbstractGeometricObject{N}}
    # gens::Vector{DisplacementGroupElem{N,T}} # generators of the group as DisplacementGroupElem
end

abstract type AbstractDisplacementGroupElem{S,T<:Real} <: AbstractAlgebra.GroupElem end

struct DisplacementGroupElem{N,T<:Real} <: AbstractDisplacementGroupElem{N,T}
    parent::DisplacementGroup{N,T}
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