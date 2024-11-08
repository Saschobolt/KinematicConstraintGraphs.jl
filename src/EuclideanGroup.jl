
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

Base.show(io::IO, G::AbstractDisplacementGroup{N,T}) where {N,T<:Real} = print(io, "$(typeof(G))( ", G.defining_object, " )")

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

PrismaticGroup(axis::AbstractGeometricSpan) = PrismaticGroup(Line(axis))
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

RevolutionGroup(axis::AbstractGeometricSpan) = RevolutionGroup(Line(axis))
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

ScrewGroup(axis::AbstractGeometricSpan, pitch::Real) = ScrewGroup(Line(axis), pitch)
ScrewGroup(p1::Point{3}, p2::Point{3}, pitch::Real) = ScrewGroup(Line(p1, p2), pitch)
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

PlanarTranslationGroup(plane::AbstractGeometricSpan) = PlanarTranslationGroup(Plane(plane))
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

CylindricalGroup(axis::AbstractGeometricSpan) = CylindricalGroup(Line(axis))
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
        return new{T}(d, dof, p)
    end
end

PlanarSlidingGroup(plane::AbstractGeometricSpan) = PlanarSlidingGroup(Plane(plane))
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

TranslatingScrewGroup(axis::AbstractGeometricSpan, pitch::Real) = TranslatingScrewGroup(Line(axis), pitch)
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

TranslatingGimbalGroup(axis::AbstractGeometricSpan) = TranslatingGimbalGroup(Line(axis))
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
    T1 = typeof(G1)
    T2 = typeof(G2)

    if (T1 <: Union{ScrewGroup,TranslatingScrewGroup} && T2 <: Union{ScrewGroup,TranslatingScrewGroup}) && G1.pitch != G2.pitch
        return false
    end

    return typeof(G1) == typeof(G2) && G1.defining_object == G2.defining_object
end

function Base.:(<=)(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    return G1 < G2 || G1 == G2
end

function intersection(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    # if any of the two groups is the subgroup of the other return the smaller group
    if G1 <= G2
        return G1
    elseif G2 <= G1
        return G2
    end

    T1 = typeof(G1)
    T2 = typeof(G2)

    # lookup table for intersection of displacement groups according to table 4.3 in the paper by Frederico Thomas
    if T1 <: PlanarTranslationGroup && T2 <: PlanarTranslationGroup
        # prismatic group along the intersection of the planes
        return PrismaticGroup(Line(intersection(G1.defining_object, G2.defining_object)))
    elseif (T1 <: PlanarTranslationGroup && T2 <: PlanarSlidingGroup) || (T2 <: PlanarTranslationGroup && T1 <: PlanarSlidingGroup)
        # prismatic group along the intersection of the planes
        return PrismaticGroup(Line(intersection(G1.defining_object, G2.defining_object)))
    elseif T1 <: PlanarSlidingGroup && T2 <: PlanarSlidingGroup
        # prismatic group along the intersection of the planes
        return PrismaticGroup(Line(intersection(G1.defining_object, G2.defining_object)))
    elseif T1 <: TranslatingScrewGroup && T2 <: PlanarTranslationGroup
        # prismatic translation group along a direction that is parallel to the plane of the translation group and perpendicular to the axis of the screw group
        # is a prismatic group as the case where plane and axis are perpendicular is caught by G2 being a subgroup of G1.
        axis1 = G1.defining_object
        plane2 = G2.defining_object

        v_screw = axis1.points[2] - axis1.points[1]
        v_plane1 = plane2.points[2] - plane2.points[1]
        v_plane2 = plane2.points[3] - plane2.points[1]

        n = nullspace(transpose(v_screw)) # translation space of the plane perpendicular to the screw axis
        @assert size(n, 2) == 2 # sanity check: n should span a plane
        m = hcat(v_plane1, v_plane2, n)
        @assert size(nullspace(m), 2) == 1 # sanity check: intersections of two non identical planes is a line
        coeffs = nullspace(m)[1:2, :] # coeffs of the linear combination of the vectors v_plane1, v_plane2 that are orthogonal to the screw axis

        l = Line(Point(0, 0, 0), Point(coeffs[1] * v_plane1 + coeffs[2] * v_plane2))

        return PrismaticGroup(l)
    elseif T1 <: PlanarTranslationGroup && T2 <: TranslatingScrewGroup
        # reversed case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingScrewGroup && T2 <: PlanarSlidingGroup
        # case distinction if the direction of the screw axis is orthogonal to the direction of the plane of the sliding group
        if orthogonal(G1.defining_object, G2.defining_object)
            # translation group along the plane of the sliding group
            return PlanarTranslationGroup(G2.defining_object)
        elseif !orthogonal(G1.defining_object, G2.defining_object)
            # translation group along a line that is perpendicular to the axis of the screw group and parallel to the plane of the sliding group
            plane = G2.defining_object
            return intersection(G1, PlanarTranslationGroup(plane))
        end
    elseif T1 <: PlanarSlidingGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingScrewGroup && T2 <: TranslatingScrewGroup
        if !parallel(G1.defining_object, G2.defining_object)
            # prismatic translation group along a direction that is perpendicular to the rotation axes of both screw groups
            axis1 = G1.defining_object
            axis2 = G2.defining_object
            v1 = axis1.points[2] - axis1.points[1]
            v2 = axis2.points[2] - axis2.points[1]
            v = cross(v1, v2)
            return PrismaticGroup(Line(Point(0, 0, 0), Point(v)))
        elseif G1.pitch != G2.pitch
            # if the groups have different pitch, but parallel axes the intersection is the planar translation group along the plane that is perpendicular to the axes
            v = G1.defining_object.points[2] - G1.defining_object.points[1]
            b = nullspace(transpose(v)) # base of translation space of plane with normal vector v
            @assert size(b, 2) == 2 # sanity check: b should span a plane
            return PlanarTranslationGroup(Plane(Point(0, 0, 0), Point(b[:, 1]), Point(b[:, 2])))
        end
    elseif T1 <: TranslatingScrewGroup && T2 <: CylindricalGroup
        # case distinction if the direction of the screw axis is parallel to the direction of the axis of the cylindrical group or if they are perpendicular. If they are neither the intersection is trivial.
        if parallel(G1.defining_object, G2.defining_object)
            # screw group along the common axis with the pitch of the translating screw group
            return ScrewGroup(G1.defining_object, G1.pitch)
        elseif orthogonal(G1.defining_object, G2.defining_object)
            # translation group along the direction of the cylindrical group
            return PrismaticGroup(G2.defining_object)
        end
    elseif T1 <: CylindricalGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: CylindricalGroup && T2 <: CylindricalGroup
        # if the axes are parallel (but not equal), the intersection is the translation group along the common axis
        if parallel(G1.defining_object, G2.defining_object)
            return PrismaticGroup(G1.defining_object)
        end
    elseif T1 <: PlanarTranslationGroup && T2 <: CylindricalGroup
        # if the plane and the axis are parallel the intersection is the translation group along the axis
        if parallel(G1.defining_object, G2.defining_object)
            return PrismaticGroup(G2.defining_object)
        end
    elseif T1 <: CylindricalGroup && T2 <: PlanarTranslationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SpatialTranslationGroup && T2 <: CylindricalGroup
        # translation group along the axis of the cylindrical group
        return PrismaticGroup(G2.defining_object)
    elseif T1 <: CylindricalGroup && T2 <: SpatialTranslationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: PlanarSlidingGroup && T2 <: CylindricalGroup
        # if the plane and the axis are parallel the intersection is the translation group along the axis
        if parallel(G1.defining_object, G2.defining_object)
            return PrismaticGroup(G2.defining_object)
        elseif orthogonal(G1.defining_object, G2.defining_object)
            # revolution group along the axis of the cylindrical group
            return RevolutionGroup(G2.defining_object)
        end
    elseif T1 <: CylindricalGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingGimbalGroup && T2 <: CylindricalGroup
        # if the axes are not parallel (which is the case as otherwise G2 <= G1) the intersection is the translation group along the axis of the cylindrical group
        return PrismaticGroup(G2.defining_object)
    elseif T1 <: CylindricalGroup && T2 <: TranslatingGimbalGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SphericalRotationGroup && T2 <: CylindricalGroup
        # if the point is on the axis of the cylindrical group the intersection is the rotation group along the axis
        if G1.defining_object in G2.defining_object
            return RevolutionGroup(G2.defining_object)
        end
    elseif T1 <: CylindricalGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SphericalRotationGroup && T2 <: PlanarSlidingGroup
        # revolution group around an axis that contains the defining point of the spherical group and that is perependicular to the plane of the sliding group
        plane = G2.defining_object
        n = normal_vec(plane)
        p = G1.defining_object
        return RevolutionGroup(Line(p, p + n))
    elseif T1 <: PlanarSlidingGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SphericalRotationGroup && T2 <: TranslatingGimbalGroup
        # revolution group around an axis that contains the defining point of the spherical group and that is parallel to the axis of the gimbal group
        axis = G2.defining_object
        v = axis.points[2] - axis.points[1]
        p = G1.defining_object
        return RevolutionGroup(Line(p, p + v))
    elseif T1 <: TranslatingGimbalGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SphericalRotationGroup && T2 <: SphericalRotationGroup
        # revolution group around the axis spanned by the two defining points of the spherical groups
        return RevolutionGroup(Line(G1.defining_object, G2.defining_object))
    elseif T1 <: TranslatingScrewGroup && T2 <: TranslatingGimbalGroup
        # if the defining axes are not parallel the intersection is planar translation group along the plane perpendicular to the axis of the screw group
        if !parallel(G1.defining_object, G2.defining_object)
            l = G1.defining_object
            v = l.points[2] - l.points[1]
            b = nullspace(transpose(v)) # base of translation space of plane with normal vector v
            @assert size(b, 2) == 2 # sanity check: b should span a plane
            return PlanarTranslationGroup(Plane(Point(0, 0, 0), Point(b[:, 1]), Point(b[:, 2])))
        end
    elseif T1 <: TranslatingGimbalGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: PlanarSlidingGroup && T2 <: TranslatingGimbalGroup
        # if the plane and the axis are not orthogonal, the intersection is the translation group along the plane
        if !orthogonal(G1.defining_object, G2.defining_object)
            return PlanarTranslationGroup(G1.defining_object)
        end
    elseif T1 <: TranslatingGimbalGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: PlanarSlidingGroup && T2 <: SpatialTranslationGroup
        # translation group along the plane
        return PlanarTranslationGroup(G1.defining_object)
    elseif T1 <: SpatialTranslationGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingScrewGroup && T2 <: SpatialTranslationGroup
        # translation group along the plane perpendicular to the screw axis
        l = G1.defining_object
        v = l.points[2] - l.points[1]
        b = nullspace(transpose(v)) # base of translation space of plane with normal vector v
        @assert size(b, 2) == 2 # sanity check: b should span a plane
        return PlanarTranslationGroup(Plane(Point(0, 0, 0), Point(b[:, 1]), Point(b[:, 2])))
    elseif T1 <: SpatialTranslationGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingGimbalGroup && T2 <: TranslatingGimbalGroup
        # spatial translation group
        return SpatialTranslationGroup()
    end

    # if non of the conditions above apply the intersection is trivial.
    return IdentityGroup(G1.d)
end

Base.inv(G::AbstractDisplacementGroup) = G


########### composition of DisplacementGroups
mutable struct DisplacementGroupComposition{N}
    factors::Vector{<:AbstractDisplacementGroup{N}}

    function DisplacementGroupComposition(groups::AbstractVector{<:AbstractDisplacementGroup{N}}) where {N}
        return new{N}(groups)
    end
end

Base.show(io::IO, comp::DisplacementGroupComposition) = print(io, comp.factors[1], [" ⋅ $(factor)" for factor in comp.factors[2:end]]...)

# reduce a composition of displacement group to a regular representation with the minimum number of trivially intersecting factors
function reduce!(comp::DisplacementGroupComposition{N}) where {N}
    if length(comp.factors) == 1
        return comp
    end

    setdiff!(comp.factors, [IdentityGroup(N)])

    # Go through factors from left to right. If any factors have non trivial intersection, compute their composition and get regular presentation. 
    # Call the function recursively on this new presentation of the composition.
    for i in 1:length(comp.factors)-1
        G1 = comp.factors[i]
        G2 = comp.factors[i+1]
        if intersection(G1, G2) != IdentityGroup(N)
            comp.factors = vcat(comp.factors[1:i-1], (G1 * G2).factors, comp.factors[i+2:end])
            return reduce!(comp)
        end
    end

    return comp
end

function is_trivial!(comp::DisplacementGroupComposition)
    return length(reduce!(comp).factors) == 1
end

function inv!(comp::DisplacementGroupComposition)
    reverse!(map!(f -> inv(f), comp.factors, comp.factors))
    return comp
end
Base.inv(comp::DisplacementGroupComposition) = DisplacementGroupComposition(inv.(reverse(comp.factors)))

Base.:(==)(comp1::DisplacementGroupComposition{N}, comp2::DisplacementGroupComposition{N}) where {N} = (reduce!(deepcopy(comp1)).factors == reduce!(deepcopy(comp2)).factors)

# regular presentation of composition of two DisplacementGroups, i.e. factor following each other have trivial intersections.
function composition(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup)
    if G1 <= G2
        return DisplacementGroupComposition([G2])
    elseif G2 <= G1
        return DisplacementGroupComposition([G1])
    end

    T1 = typeof(G1)
    T2 = typeof(G2)

    # cases where the two groups don't intersect trivially. See table 4.3 in Thomas 1991
    if T1 <: PlanarTranslationGroup && T2 <: PlanarTranslationGroup
        # Spatial translation group
        return DisplacementGroupComposition([SpatialTranslationGroup()])
    elseif T1 <: PlanarTranslationGroup && T2 <: PlanarSlidingGroup
        # translating gimbal with axis orthogonal to plane of the sliding group
        v = normal_vec(G2.defining_object)
        return DisplacementGroupComposition([TranslatingGimbalGroup(Line(Point(0, 0, 0), Point(v)))])
    elseif T2 <: PlanarTranslationGroup && T1 <: PlanarSlidingGroup
        # reverse case of above
        return composition(G2, G1)
    elseif T1 <: PlanarSlidingGroup && T2 <: PlanarSlidingGroup
        # R_{u_0} T R_{u_1}, where u_0 is orthogonal to first plane and u_1 is orthogonal to the second plane
        u0 = normal_vec(G1.defining_object)
        u1 = normal_vec(G2.defining_object)
        return DisplacementGroupComposition([RevolutionGroup(Line(Point(0, 0, 0), Point(u0))), SpatialTranslationGroup(), RevolutionGroup(Line(Point(0, 0, 0), Point(u1)))])
    elseif T1 <: TranslatingScrewGroup && T2 <: PlanarTranslationGroup
        # translating gimbal with same axis as screw group
        return DisplacementGroupComposition([TranslatingGimbalGroup(G1.defining_object)])
    elseif T1 <: PlanarTranslationGroup && T2 <: TranslatingScrewGroup
        # reversed case of the case before
        return intersection(G2, G1)
    elseif T1 <: TranslatingScrewGroup && T2 <: PlanarSlidingGroup
        # case distinction if the direction of the screw axis is orthogonal to the direction of the plane of the sliding group
        if orthogonal(G1.defining_object, G2.defining_object)
            # translating gimbal with axis of screw group
            return DisplacementGroupComposition([TranslatingGimbalGroup(G1.defining_object)])
        elseif !orthogonal(G1.defining_object, G2.defining_object)
            # X_{v_0} R_{u_0} where v_0 is axis of the screw group and u_0 is orthogonal to plane of sliding group
            v0 = G1.defining_object
            u0 = Line(Point(0, 0, 0), Point(normal_vec(G2.defining_object)))
            return DisplacementGroupComposition([TranslatingGimbalGroup(v0), RevolutionGroup(u0)])
        end
    elseif T1 <: PlanarSlidingGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        if orthogonal(G1.defining_object, G2.defining_object)
            return composition(G2, G1)
        elseif !orthogonal(G1.defining_object, G2.defining_object)
            # R_{u_0} X_{v_0} where v_0 is axis of the screw group and u_0 is orthogonal to plane of sliding group
            v0 = G2.defining_object
            u0 = Line(Point(0, 0, 0), Point(normal_vec(G1.defining_object)))
            return DisplacementGroupComposition([RevolutionGroup(u0), TranslatingGimbalGroup(v0)])
        end
    elseif T1 <: TranslatingScrewGroup && T2 <: TranslatingScrewGroup
        if !parallel(G1.defining_object, G2.defining_object)
            # R_{u_0} T R_{u_1} where u_0 is axis of G1 and u_1 is axis of G2
            u0 = G1.defining_object
            u1 = G2.defining_object
            return DisplacementGroupComposition([RevolutionGroup(u0), SpatialTranslationGroup(), RevolutionGroup(u1)])
        elseif G1.pitch != G2.pitch
            # translating gimbal group with same axis as one of the screw groups
            return DisplacementGroupComposition([TranslatingGimbalGroup(G1.defining_object)])
        end
    elseif T1 <: TranslatingScrewGroup && T2 <: CylindricalGroup
        # case distinction if the direction of the screw axis is parallel to the direction of the axis of the cylindrical group or if they are perpendicular. If they are neither the intersection is trivial.
        if parallel(G1.defining_object, G2.defining_object)
            # translating gimbal group with the axis of any of the two groups
            return DisplacementGroupComposition([TranslatingGimbalGroup(G1.defining_object)])
        elseif orthogonal(G1.defining_object, G2.defining_object)
            # Y_{v_0, p_0} R_{u_0} where v_0 is the axis of the screw group, p_0 is the pitch of the screw group and u_0 is the axis of the cylindrical group
            v0 = G1.defining_object
            p0 = G1.pitch
            u0 = G2.defining_object
            return DisplacementGroupComposition([TranslatingScrewGroup(v0, p0), RevolutionGroup(u0)])
        end
    elseif T1 <: CylindricalGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        if parallel(G1.defining_object, G2.defining_object)
            return composition(G2, G1)
        elseif orthogonal(G1.defining_object, G2.defining_object)
            # R_{u_0} Y_{v_0, p_0} where v_0 is the axis of the screw group, p_0 is the pitch of the screw group and u_0 is the axis of the cylindrical group
            v0 = G2.defining_object
            p0 = G2.pitch
            u0 = G1.defining_object
            return DisplacementGroupComposition([RevolutionGroup(u0), TranslatingScrewGroup(v0, p0)])
        end
    elseif T1 <: CylindricalGroup && T2 <: CylindricalGroup
        # if the axes are parallel (but not equal), the composition is C_u0 R_u1 where u0 is the axis of the first and u1 the axis of the second group
        if parallel(G1.defining_object, G2.defining_object)
            u0 = G1.defining_object
            u1 = G2.defining_object
            return DisplacementGroupComposition([CylindricalGroup(u0), RevolutionGroup(u1)])
        end
    elseif T1 <: PlanarTranslationGroup && T2 <: CylindricalGroup
        # if the plane and the axis are parallel the composition is T_P0 R_u0 where P0 is the plane of the translation group and u0 is the axis of the cylindrical group
        if parallel(G1.defining_object, G2.defining_object)
            P0 = G1.defining_object
            u0 = G2.defining_object
            return DisplacementGroupComposition([PlanarTranslationGroup(P0), RevolutionGroup(u0)])
        end
    elseif T1 <: CylindricalGroup && T2 <: PlanarTranslationGroup
        # reverse case of the case before
        if parallel(G1.defining_object, G2.defining_object)
            P0 = G2.defining_object
            u0 = G1.defining_object
            return DisplacementGroupComposition([RevolutionGroup(u0), PlanarTranslationGroup(P0)])
        end
    elseif T1 <: SpatialTranslationGroup && T2 <: CylindricalGroup
        # X_v0 where v0 is the axis of the cylindrical group
        v0 = G2.defining_object
        return DisplacementGroupComposition([TranslatingGimbalGroup(v0)])
    elseif T1 <: CylindricalGroup && T2 <: SpatialTranslationGroup
        # reverse case of the case before
        return composition(G2, G1)
    elseif T1 <: PlanarSlidingGroup && T2 <: CylindricalGroup
        # if the plane and the axis are parallel the composition is G_P0 R_u0 where P0 is the plane of the sliding group and u0 is the axis of the cylindrical group
        if parallel(G1.defining_object, G2.defining_object)
            P0 = G1.defining_object
            u0 = G2.defining_object
            return DisplacementGroupComposition([PlanarSlidingGroup(P0), RevolutionGroup(u0)])
        elseif orthogonal(G1.defining_object, G2.defining_object)
            # X_v0 where v0 is orthogonal to P0
            v0 = normal_vec(G1.defining_object)
            return DisplacementGroupComposition([TranslatingGimbalGroup(v0)])
        end
    elseif T1 <: CylindricalGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        if parallel(G1.defining_object, G2.defining_object)
            P0 = G2.defining_object
            u0 = G1.defining_object
            return DisplacementGroupComposition([PlanarSlidingGroup(P0), RevolutionGroup(u0)])
        elseif orthogonal(G1.defining_object, G2.defining_object)
            return composition(G2, G1)
        end
    elseif T1 <: TranslatingGimbalGroup && T2 <: CylindricalGroup
        # X_v0 R_u0 where v0 is axis of the gimbal group and u0 is axis of the cylindrical group
        v0 = G1.defining_object
        u0 = G2.defining_object
        return DisplacementGroupComposition([TranslatingGimbalGroup(v0), RevolutionGroup(u0)])
    elseif T1 <: CylindricalGroup && T2 <: TranslatingGimbalGroup
        # reverse case of the case before
        v0 = G2.defining_object
        u0 = G1.defining_object
        return DisplacementGroupComposition([RevolutionGroup(u0), TranslatingGimbalGroup(v0)])
    elseif T1 <: SphericalRotationGroup && T2 <: CylindricalGroup
        # S_x T_u0 where x is the point of the rotation group and u0 the axis of the cylindrical group, if x is in u_0
        x = G1.defining_object
        u0 = G2.defining_object
        if x in u0
            return DisplacementGroupComposition([SphericalRotationGroup(x), PrismaticGroup(u0)])
        end
    elseif T1 <: CylindricalGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        x = G2.defining_object
        u0 = G1.defining_object
        if x in u0
            return DisplacementGroupComposition([PrismaticGroup(u0), SphericalRotationGroup(x)])
        end
    elseif T1 <: SphericalRotationGroup && T2 <: PlanarSlidingGroup
        # S_x T_P0 where x is the point of the rotation group and P0 the plane of the sliding group
        x = G1.defining_object
        P0 = G2.defining_object
        return DisplacementGroupComposition([SphericalRotationGroup(x), PlanarTranslationGroup(P0)])
    elseif T1 <: PlanarSlidingGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        x = G2.defining_object
        P0 = G1.defining_object
        return DisplacementGroupComposition([PlanarTranslationGroup(P0), SphericalRotationGroup(x)])
    elseif T1 <: SphericalRotationGroup && T2 <: TranslatingGimbalGroup
        # SE(3)
        return SpecialEuclideanGroup(3)
    elseif T1 <: TranslatingGimbalGroup && T2 <: SphericalRotationGroup
        # reverse case of the case before
        return intersection(G2, G1)
    elseif T1 <: SphericalRotationGroup && T2 <: SphericalRotationGroup
        # S_x R_u0 R_u1 where x is the defining point of G1 and the defining point of G2 is in both axis u1 and u2
        x = G1.defining_object
        y = G2.defining_object
        u0 = Line(y, y + [1, 0, 0])
        u1 = Line(y, y + [0, 1, 0])
        return DisplacementGroupComposition([SphericalRotationGroup(x), RevolutionGroup(u0), RevolutionGroup(u1)])
    elseif T1 <: TranslatingScrewGroup && T2 <: TranslatingGimbalGroup
        # if the defining axes are not parallel the composition is X_v0 R_u1 where v0 is the axis of the screw group and u1 is the axis of the Gimbal group
        if !parallel(G1.defining_object, G2.defining_object)
            v0 = G1.defining_object
            u1 = G2.defining_object
            return DisplacementGroupComposition([TranslatingGimbalGroup(v0), RevolutionGroup(u1)])
        end
    elseif T1 <: TranslatingGimbalGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        if !parallel(G1.defining_object, G2.defining_object)
            v0 = G2.defining_object
            u1 = G1.defining_object
            return DisplacementGroupComposition([RevolutionGroup(u1), TranslatingGimbalGroup(v0)])
        end
    elseif T1 <: PlanarSlidingGroup && T2 <: TranslatingGimbalGroup
        # if the plane and the axis are not orthogonal, the composition is R_u0 T R_u1 where u0 is orthogonal to the plane and u1 is parallel to the axis of the gimbal group
        if !orthogonal(G1.defining_object, G2.defining_object)
            P0 = G1.defining_object
            v1 = G2.defining_object
            u1 = v1
            u0 = Line(Point(0, 0, 0), Point(normal_vec(P0)))
            return DisplacementGroupComposition([RevolutionGroup(u0), SpatialTranslationGroup(), RevolutionGroup(u1)])
        end
    elseif T1 <: TranslatingGimbalGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        if !orthogonal(G1.defining_object, G2.defining_object)
            P0 = G2.defining_object
            v1 = G1.defining_object
            u1 = v1
            u0 = Line(Point(0, 0, 0), Point(normal_vec(P0)))
            return DisplacementGroupComposition([RevolutionGroup(u1), SpatialTranslationGroup(), RevolutionGroup(u0)])
        end
    elseif T1 <: PlanarSlidingGroup && T2 <: SpatialTranslationGroup
        # X_v0 where v0 is orthogonal to the plane of the sliding group
        v0 = Line(Point(0, 0, 0), normal_vec(G2.defining_object))
        return DisplacementGroupComposition([TranslatingGimbalGroup(v0)])
    elseif T1 <: SpatialTranslationGroup && T2 <: PlanarSlidingGroup
        # reverse case of the case before
        return composition(G2, G1)
    elseif T1 <: TranslatingScrewGroup && T2 <: SpatialTranslationGroup
        # X_v0 where v0 is the axis of the screw group
        v0 = G1.defining_object
        return DisplacementGroupComposition(TranslatingGimbalGroup(v0))
    elseif T1 <: SpatialTranslationGroup && T2 <: TranslatingScrewGroup
        # reverse case of the case before
        return composition(G2, G1)
    elseif T1 <: TranslatingGimbalGroup && T2 <: TranslatingGimbalGroup
        # R_u0 T R_u1 where u0 is parallel to the axis of G1 and u1 is parallel to the axis of G2
        u0 = G1.defining_object
        u1 = G2.defining_object
        return DisplacementGroupComposition([RevolutionGroup(u0), SpatialTranslationGroup, RevolutionGroup(u1)])
    end

    # otherwise the intersection of the groups is empty
    return DisplacementGroupComposition([G1, G2])
end

Base.:*(G1::Abstract3dDisplacementGroup, G2::Abstract3dDisplacementGroup) = composition(G1, G2)

# composition of DisplacementGroup and DisplacementGroupComposition
function composition(comp::DisplacementGroupComposition{N}, G::AbstractDisplacementGroup{N}) where {N}
    return reduce!(DisplacementGroupComposition(push!(deepcopy(comp.factors), G)))
end

Base.:*(comp::DisplacementGroupComposition{N}, G::AbstractDisplacementGroup{N}) where {N} = composition(comp, G)

function composition(G::AbstractDisplacementGroup{N}, comp::DisplacementGroupComposition{N}) where {N}
    return reduce!(DisplacementGroupComposition(pushfirst!(deepcopy(comp.factors), G)))
end

Base.:*(G::AbstractDisplacementGroup{N}, comp::DisplacementGroupComposition{N}) where {N} = composition(G, comp)

function composition(comp1::DisplacementGroupComposition, comp2::DisplacementGroupComposition)
    return reduce!(DisplacementGroupComposition(vcat(comp1.factors, comp2.factors)))
end

Base.:*(comp1::DisplacementGroupComposition, comp2::DisplacementGroupComposition) = composition(comp1, comp2)

intersection(comp::DisplacementGroupComposition{N}, G::AbstractDisplacementGroup{N}) where {N} = is_trivial!(comp) ? (return DisplacementGroupComposition([intersection(comp.factors[1], G)])) : (throw(ArgumentError("composition has to be trivial.")))

intersection(G::AbstractDisplacementGroup{N}, comp::DisplacementGroupComposition{N}) where {N} = is_trivial!(comp) ? (return DisplacementGroupComposition([intersection(comp.factors[1], G)])) : (throw(ArgumentError("composition has to be trivial.")))

intersection(comp1::DisplacementGroupComposition{N}, comp2::DisplacementGroupComposition{N}) where {N} = (is_trivial!(comp1) && is_trivial!(comp2)) ? (return DisplacementGroupComposition([intersection(comp1.factors[1], comp2.factors[1])])) : (throw(ArgumentError("compositions have to be trivial.")))