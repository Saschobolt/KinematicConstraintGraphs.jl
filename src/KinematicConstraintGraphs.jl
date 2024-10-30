module KinematicConstraintGraphs

import AbstractAlgebra
import Graphs
using StaticArrays
using LinearAlgebra
# Write your package code here.
export AbstractGeometricSpan, Point, GeometricSpan, Line, Plane, dim, basis, matrix, parallel, collinear, orthogonal, intersection
export Abstract3dDisplacementGroup
export IdentityGroup, SpecialEuclideanGroup, PrismaticGroup, RevolutionGroup, ScrewGroup, PlanarTranslationGroup, CylindricalGroup, SpatialTranslationGroup, PlanarSlidingGroup, SphericalRotationGroup, TranslatingScrewGroup, TranslatingGimbalGroup, is_subgroup, intersection

include("Geometry.jl")
include("EuclideanGroup.jl")
end