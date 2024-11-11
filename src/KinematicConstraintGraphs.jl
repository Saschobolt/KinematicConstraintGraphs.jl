module KinematicConstraintGraphs

import AbstractAlgebra
import Graphs
using GraphMakie
using StaticArrays
using LinearAlgebra
# Write your package code here.
export AbstractGeometricSpan, Point, GeometricSpan, Line, Plane, dim, basis, matrix, parallel, collinear, orthogonal, intersection
export Abstract3dDisplacementGroup
export inv!, IdentityGroup, SpecialEuclideanGroup, PrismaticGroup, RevolutionGroup, ScrewGroup, PlanarTranslationGroup, CylindricalGroup, SpatialTranslationGroup, PlanarSlidingGroup, SphericalRotationGroup, TranslatingScrewGroup, TranslatingGimbalGroup, is_subgroup, intersection
export ConstraintGraph, add_vertex!, add_constraint!, blocks, filter_constraints!, filter_constraints, all_equiv_constraints!, all_equiv_constraints, is_interlocking

include("Geometry.jl")
include("EuclideanGroup.jl")
include("constraint_graph.jl")
end