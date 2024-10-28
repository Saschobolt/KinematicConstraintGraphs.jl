
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

######################## displacement groups in R^3. See Graphs of kinematic constraints by Frederico Thomas, http://www.iri.upc.edu/people/thomas/papers/Computer-Aided.pdf
"""
Abstract supertype of all displacement groups. A displacement group is a group of transformations that act on a euclidean space.
All subtypes of AbstractDisplacementGroup must have the following fields:
- d::Int: dimension of euclidean space that the group acts on. E.g R^3 -> d = 3
- dof::Int: degrees of freedom of the group. E.g SE(3) -> dof = 6
- defining_object<:Union{Nothing, AbstractGeometricObject{N}}: object that defines the group. E.g. for translational group along a plane in R^3, defining_object is the defining line
"""
abstract type AbstractDisplacementGroup{N,T<:Real} <: AbstractAlgebra.Group end
Abstract3dDisplacementGroup{T<:Real} = AbstractDisplacementGroup{3,T}

# identity group
struct IdentityGroup{N,T<:Real} <: AbstractDisplacementGroup{N,T}
    d::Int
    dof::Int
    defining_object::Nothing

    function IdentityGroup{N,T}() where {N,T<:Real}
        d = N
        dof = 0
        return new{N,T}(d, dof, nothing)
    end
end

IdentityGroup(d::Integer=3) = IdentityGroup{d,Float64}()

# special euclidean Group
struct SpecialEuclideanGroup{N,T<:Real} <: AbstractDisplacementGroup{N,T}
    d::Int
    dof::Int
    defining_object::Nothing

    function SpecialEuclideanGroup{N,T}() where {N,T<:Real}
        d = N
        dof = binomial(N + 1, 2)
        new{N,Float64}(d, dof, nothing)
    end
end

SpecialEuclideanGroup(d::Integer=3) = SpecialEuclideanGroup{d,Float64}()

# rectilinear translation (P) Prismatic
struct PrismaticGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}

    function PrismaticGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0))) where {T<:Real}
        d = 3
        dof = 1
        # translate axis to a normal form: axis is a line that passes through the origin and the points of axis have distance 1.
        ax = Line(Point(0, 0, 0), Point(normalize(axis.points[2] - axis.points[1])))
        new{T}(d, dof, ax)
    end
end

PrismaticGroup(p1::Point{3}, p2::Point{3}) = PrismaticGroup(Line(p1, p2))
PrismaticGroup(p::Point{3}, v::AbstractVector) = PrismaticGroup(Line(p, p + v))

# Roation around an axis (R) Revolution
struct RevolutionGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}

    function RevolutionGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0))) where {T<:Real}
        d = 3
        dof = 1
        new{T}(d, dof, axis)
    end
end

RevolutionGroup(p1::Point{3}, p2::Point{3}) = RevolutionGroup(Line(p1, p2))
RevolutionGroup(p::Point{3}, v::AbstractVector) = RevolutionGroup(Line(p, p + v))

# Helicoidal movement (H) Screw
struct ScrewGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}
    pitch::Float64 # pitch of the screw

    function ScrewGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0)), pitch::Real=1.0) where {T<:Real}
        d = 3
        dof = 1
        new{T}(d, dof, axis, pitch)
    end
end

ScrewGroup(p1::Point{3}, p2::Point{3}) = ScrewGroup(Line(p1, p2))
ScrewGroup(p::Point{3}, v::AbstractVector, pitch::Real) = ScrewGroup(Line(p, p + v), pitch)

# Planar translation
struct PlanarTranslationGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Plane{3,T}

    function PlanarTranslationGroup(plane::Plane{3,T}=Plane(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))) where {T<:Real}
        d = 3
        dof = 2
        # translate plane to a normal form: plane is a plane that passes through the origin and the points of plane have distance 1.
        p = Plane(Point(0, 0, 0), Point(normalize(plane.points[2] - plane.points[1])), Point(normalize(plane.points[3] - plane.points[1])))
        new{T}(d, dof, p)
    end
end

PlanarTranslationGroup(p1::Point{3}, p2::Point{3}, p3::Point{3}) = PlanarTranslationGroup(Plane(p1, p2, p3))
PlanarTranslationGroup(p::Point{3}, v1::AbstractVector, v2::AbstractVector) = PlanarTranslationGroup(Plane(p, p + v1, p + v2))

# Lock movement (C) Cylindrical
struct CylindricalGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}

    function CylindricalGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0))) where {T<:Real}
        d = 3
        dof = 2
        new{T}(d, dof, axis)
    end
end

CylindricalGroup(p1::Point{3}, p2::Point{3}) = CylindricalGroup(Line(p1, p2))
CylindricalGroup(p::Point{3}, v::AbstractVector) = CylindricalGroup(Line(p, p + v))

# Spatial translation
struct SpatialTranslationGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Nothing

    function SpatialTranslationGroup{T}() where {T<:Real}
        d = 3
        dof = 3
        return new{T}(d, dof, nothing)
    end
end

SpatialTranslationGroup() = SpatialTranslationGroup{Float64}()

# Planar sliding (E) Plane
struct PlanarSlidingGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Plane{3,T}

    function PlanarSlidingGroup(plane::Plane{3,T}=Plane(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1))) where {T<:Real}
        d = 3
        dof = 3
        # translate plane to a normal form: plane is a plane that passes through the origin and the points of plane have distance 1.
        p = Plane(Point(0, 0, 0), Point(normalize(plane.points[2] - plane.points[1])), Point(normalize(plane.points[3] - plane.points[1])))
        return new{T}(d, dof, plane)
    end
end

PlanarSlidingGroup(p1::Point{3}, p2::Point{3}, p3::Point{3}) = PlanarSlidingGroup(Plane(p1, p2, p3))
PlanarSlidingGroup(p::Point{3}, v1::AbstractVector, v2::AbstractVector) = PlanarSlidingGroup(Plane(p, p + v1, p + v2))

# Spheric rotation (S) Spherical
struct SphericalRotationGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Point{3,T}

    function SphericalRotationGroup(point::Point{3,T}=Point(0, 0, 0)) where {T<:Real}
        d = 3
        dof = 3
        return new{T}(d, dof, point)
    end
end

# Translating Screw
struct TranslatingScrewGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}
    pitch::Real

    function TranslatingScrewGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0)), pitch::Real=1) where {T<:Real}
        d = 3
        dof = 3
        # translate axis to a normal form: axis is a line that passes through the origin and the points of axis have distance 1.
        ax = Line(Point(0, 0, 0), Point(normalize(axis.points[2] - axis.points[1])))
        return new{T}(d, dof, axis, pitch)
    end
end

TranslatingScrewGroup(p1::Point{3}, p2::Point{3}, pitch::Real) = TranslatingScrewGroup(Line(p1, p2), pitch)
TranslatingScrewGroup(p::Point{3}, v::AbstractVector, pitch::Real) = TranslatingScrewGroup(Line(p, p + v), pitch)

# Translating gimbal
struct TranslatingGimbalGroup{T<:Real} <: Abstract3dDisplacementGroup{T}
    d::Int
    dof::Int
    defining_object::Line{3,T}

    function TranslatingGimbalGroup(axis::Line{3,T}=Line(Point(0, 0, 0), Point(1, 0, 0))) where {T<:Real}
        d = 3
        dof = 4
        # translate axis to a normal form: axis is a line that passes through the origin and the points of axis have distance 1.
        ax = Line(Point(0, 0, 0), Point(normalize(axis.points[2] - axis.points[1])))
        return new{T}(d, dof, axis)
    end
end

TranslatingGimbalGroup(p1::Point{3}, p2::Point{3}) = TranslatingGimbalGroup(Line(p1, p2))
TranslatingGimbalGroup(p::Point{3}, v::AbstractVector) = TranslatingGimbalGroup(Line(p, p + v))

######## subgroup relations among displacement groups

function is_subgroup(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    if G1.dof > G2.dof
        return false
    end

    if typeof(G2) <: SpecialEuclideanGroup
        return true
    end

    if typeof(G1) <: PrismaticGroup
        if typeof(G2) <: PlanarTranslationGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: CylindricalGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: SpatialTranslationGroup
            return true
        elseif typeof(G2) <: PlanarSlidingGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: TranslatingScrewGroup
            return orthogonal(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: TranslatingGimbalGroup
            return true
        end
    elseif typeof(G1) <: RevolutionGroup
        if typeof(G2) <: CylindricalGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: PlanarSlidingGroup
            return orthogonal(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: SphericalRotationGroup
            return (G2.defining_object in G1.defining_object)
        elseif typeof(G2) <: TranslatingGimbalGroup
            return parallel(G1.defining_object, G2.defining_object)
        end
    elseif typeof(G1) <: ScrewGroup
        if typeof(G2) <: CylindricalGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: TranslatingScrewGroup
            return (G1.defining_object == G2.defining_object) && G1.pitch == G2.pitch
        elseif typeof(G2) <: TranslatingGimbalGroup
            return parallel(G1.defining_object, G2.defining_object)
        end
    elseif typeof(G1) <: PlanarTranslationGroup
        if typeof(G2) <: SpatialTranslationGroup
            return true
        elseif typeof(G2) <: PlanarSlidingGroup
            return parallel(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: TranslatingScrewGroup
            return orthogonal(G1.defining_object, G2.defining_object)
        elseif typeof(G2) <: TranslatingGimbalGroup
            return true
        end
    elseif typeof(G1) <: CylindricalGroup
        if typeof(G2) <: TranslatingGimbalGroup
            return parallel(G1.defining_object, G2.defining_object)
        end
    elseif typeof(G1) <: SpatialTranslationGroup
        if typeof(G2) <: TranslatingGimbalGroup
            return true
        end
    elseif typeof(G1) <: PlanarSlidingGroup
        if typeof(G2) <: TranslatingGimbalGroup
            return orthogonal(G1.defining_object, G2.defining_object)
        end
    elseif typeof(G1) <: SphericalRotationGroup
        return false
    elseif typeof(G1) <: TranslatingScrewGroup
        if typeof(G2) <: TranslatingGimbalGroup
            return parallel(G1.defining_object, G2.defining_object)
        end
    elseif typeof(G1) <: TranslatingGimbalGroup
        return false
    elseif typeof(G1) <: IdentityGroup
        return true
    end

    return false
end

function Base.:(<)(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    return is_subgroup(G1, G2)
end

function Base.:(==)(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    return typeof(G1) == typeof(G2) && G1.defining_object == G2.defining_object
end

function Base.:(<=)(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    return G1 < G2 || G1 == G2
end
