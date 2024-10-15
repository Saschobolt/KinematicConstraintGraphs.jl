######################## linear algebra
# indices of linear independent columns of a matrix
function indcols_indices(A::AbstractMatrix{<:Real}; atol::Real=1e-8)
    Q, R, perm = qr(A, ColumnNorm())
    indices = perm[findall(val -> abs(val) > atol, diag(R))]
    return indices
end

# linear independent columns of a matrix
function indcols(A::AbstractMatrix{<:Real}; atol::Real=1e-8)
    indices = indcols_indices(A, atol=atol)
    return A[:, indices]
end


######################## euclidean geometry
# supertype of all geometric objects in R^N
abstract type AbstractGeometricObject{N} end

############ point
struct Point{N,T<:Real} <: AbstractGeometricObject{N}
    coords::SVector{N,T}

    Point{N,T}(coords::AbstractVector{<:Real}) where {N,T<:Real} = new(SVector{N}(coords))
end

function Point(coords::AbstractVector{<:Real})
    N = length(coords)
    return Point{N,eltype(coords)}(coords)
end

function Point(coords...)
    return Point(collect(coords))
end

# promotion and conversion
Base.convert(::Type{Point{N,T}}, p::Point{M,S}) where {M,N,T<:Real,S<:Real} = Point{N,T}(p.coords)
Base.promote_rule(::Type{Point{N,T}}, ::Type{Point{N,S}}) where {N,T<:Real,S<:Real} = Point{N,promote_type(T, S)}

Base.length(p::Point) = length(p.coords)

############ span of points
abstract type AbstractGeometricSpan{N} <: AbstractGeometricObject{N} end

mutable struct GeometricSpan{N,T<:Real} <: AbstractGeometricSpan{N}
    points::Vector{Point{N,T}}

    function GeometricSpan{N,T}(points::AbstractVector{Point{N,<:Real}}) where {N,T<:Real}
        return new{N,T}(collect(promote(points...)))
    end
end

GeometricSpan(points::AbstractVector{Point{N,T}}) where {N,T<:Real} = GeometricSpan{N,T}(points)
GeometricSpan(points::AbstractVector{Point{N,<:Real}}) where {N} = GeometricSpan{N,eltype(vcat([point.coords for point in points]...))}(points)
GeometricSpan(points::AbstractVector{Point}) = allequal(length.(points)) ? GeometricSpan(points) : throw(ArgumentError("All points must have the same dimension."))

GeometricSpan(points...) = GeometricSpan(collect(points))

function matrix(s::GeometricSpan{N,T}) where {N,T<:Real}
    mat = hcat([point.coords for point in s.points]...)
    mat = vcat(mat, ones(1, size(mat, 2)))
    return SMatrix{N + 1,length(s.points),T}(mat)
end

function dim(s::GeometricSpan)
    return rank(matrix(s)) - 1
end

# basis of the GeometricSpan, i.e. the points that span the GeometricSpan, i.e. the points that corresponds to the generators of the col space of the matrix of the GeometricSpan
function basis(s::GeometricSpan)
    basis_mat = indcols(matrix(s))
    return [Point(basis_mat[:, i]) for i in 1:size(basis_mat, 2)]
end

############ line
mutable struct Line{N,T<:Real} <: AbstractGeometricSpan{N}
    points::SVector{2,Point{N,T}}
end

############ plane
mutable struct Plane{N,T<:Real} <: AbstractGeometricSpan{N}
    point::SVector{3,Point{N,T}}
end