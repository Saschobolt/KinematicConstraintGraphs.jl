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
abstract type AbstractGeometricSpan{N,T<:Real} <: AbstractGeometricObject{N} end

mutable struct GeometricSpan{N,T<:Real} <: AbstractGeometricSpan{N,T}
    points::Vector{Point{N,T}}

    function GeometricSpan{N,T}(points::AbstractVector{<:Point{N}}) where {N,T<:Real}
        s = new{N,T}(points)
        return new{N,T}(basis(s))
    end
end

function GeometricSpan{N}(points::AbstractVector{<:Point{N}}) where {N}
    T = eltype(vcat([point.coords for point in points]...))
    return GeometricSpan{N,T}(points)
end

function GeometricSpan(points::AbstractVector{<:Point{N}}) where {N}
    return GeometricSpan{N}(points)
end

GeometricSpan(points::AbstractVector{<:Point}) = allequal(length.(points)) ? GeometricSpan(points) : throw(ArgumentError("All points must have the same dimension."))

function GeometricSpan(points...)
    if any([!(typeof(p) <: Point) for p in points])
        throw(ArgumentError("All points supplied need to be a subtype of Point."))
    end

    return GeometricSpan(collect(points))
end
function matrix(s::AbstractGeometricSpan{N,T}) where {N,T<:Real}
    mat = hcat([point.coords for point in s.points]...)
    mat = vcat(mat, ones(1, size(mat, 2)))
    return SMatrix{N + 1,length(s.points),T}(mat)
end

function dim(s::AbstractGeometricSpan)
    return rank(matrix(s)) - 1
end

# basis of the GeometricSpan, i.e. the points that span the GeometricSpan, i.e. the points that corresponds to the generators of the col space of the matrix of the GeometricSpan
function basis(s::AbstractGeometricSpan)
    basis_mat = indcols(matrix(s))[1:end-1, :]
    return [Point(basis_mat[:, i]) for i in axes(basis_mat, 2)]
end

# basis of the translation space of the AbstractGeometricSpan{N} s. Result is a N x dim(s) - matrix with last row equal to zero.
function basis_translation_space(s::AbstractGeometricSpan{N,T}) where {N,T<:Real}
    if dim(s) == 0
        return zeros(N, 0)
    end

    mat_points = matrix(s)
    mat_translations = hcat([mat_points[:, i] - mat_points[:, 1] for i in axes(mat_points, 2)[2:end]]...)
    return indcols(mat_translations[1:end-1, :])
end

# span of multiple geometric spans
span(geom_spans...) = GeometricSpan(vcat([s.points for s in geom_spans]...))

# return whether p is in s
function Base.in(p::Point, s::AbstractGeometricSpan)
    point_s = s.points[1]

    b = p.coords - point_s.coords
    A = basis_translation_space(s)

    return A * (A \ b) ≈ b
end

# intersection of two spans. If it is empty return nothing
function intersection(s1::AbstractGeometricSpan, s2::AbstractGeometricSpan)
    d1 = dim(s1)
    d2 = dim(s2)

    p1 = s1.points[1].coords
    p2 = s2.points[1].coords

    if d1 == 0 && d2 == 0
        p1 == p2 ? (return GeometricSpan(s1.points[1])) : (return nothing)
    elseif d1 == 0 && d2 > 0
        s1.points[1] in s2 ? (return GeometricSpan(s1.points[1])) : (return nothing)
    elseif d1 > 0 && d2 == 0
        s2.points[1] in s1 ? (return GeometricSpan(s2.points[1])) : (return nothing)
    end

    t1 = basis_translation_space(s1)
    t2 = basis_translation_space(s2)

    A = hcat(t1, -t2) # matrix that defines linear system of equations for two spans to intersect. 
    A_ind = indcols(A) # For the system of eqs we take the independet columns of A to avoid a LinearAlgebra.SingularException when solving the system.
    b = p2 - p1 # right hand side of system of linear equations
    sol = A_ind \ b #  As the cols of t1 are linearly independent, A_ind has more than d1 coumns and the first d1 entries of sol describe a point in s1 by p1 + t1 * sol[1:d1].

    if !(A_ind * sol ≈ b) # sanity check that the solution is actually a solution of the system of equations and thus an element in the intersection. If it is not, then the intersection is empty.
        return nothing
    end

    n = nullspace(A) # kernel of A. Dimension of kernel is dimension of the intersection

    # points that span the intersection can be computed by taking the first point in the solution and translating this point with translations corresponding to basis elements in the nullspace of A (columns of n).
    p = p1 + t1 * sol[1:d1]
    points = [Point(p), [Point(p1 + t1 * (sol[1:d1] + n[1:d1, i])) for i in axes(n, 2)]...]
    return GeometricSpan(points)
end

############ line
mutable struct Line{N,T<:Real} <: AbstractGeometricSpan{N,T}
    points::SVector{2,Point{N,T}}
end

function Line(s::AbstractGeometricSpan)
    @assert dim(s) == 1 "Points don't span a line, but have a geometric span of dimension $(dim(s))."
    b = basis(s)

    return Line(SVector{2}([b[1], b[2]]))
end

Line(points::AbstractVector{Point}) = Line(GeometricSpan(points))
Line(points...) = Line(GeometricSpan(points...))

############ plane
mutable struct Plane{N,T<:Real} <: AbstractGeometricSpan{N,T}
    point::SVector{3,Point{N,T}}
end

function Plane(s::AbstractGeometricSpan)
    @assert dim(s) == 2 "Points don't span a plane, but have a geometric span of dimension $(dim(s))."
    b = basis(s)

    return Plane(SVector{3}([b[1], b[2], b[3]]))
end

Plane(points) = Plane(GeometricSpan(points))
Plane(points...) = Plane(GeometricSpan(points...))