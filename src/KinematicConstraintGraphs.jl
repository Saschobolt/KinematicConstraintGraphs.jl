module KinematicConstraintGraphs

import AbstractAlgebra
import Graphs
using StaticArrays
using LinearAlgebra
# Write your package code here.
export Point, GeometricSpan, Line, Plane, dim, basis, matrix, parallel, collinear, orthogonal
export Abstract3dDisplacementGroup
export IdentityGroup, SpecialEuclideanGroup, PrismaticGroup, RevolutionGroup, ScrewGroup, PlanarTranslationGroup, CylindricalGroup, SpatialTranslationGroup, PlanarSlidingGroup, SphericalRotationGroup, TranslatingScrewGroup, TranslatingGimbalGroup

include("Geometry.jl")
include("EuclideanGroup.jl")
end
