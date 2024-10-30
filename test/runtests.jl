using KinematicConstraintGraphs
using Test
using LinearAlgebra

@testset "geometry" begin

end

###### randomized examples -> lead to numerical errors in calculations, which let tests fail.
# l1 = Line(Point(rand(3)), Point(rand(3))) # random line in 3D
# b = nullspace(transpose(l1.points[2] - l1.points[1])) # translational basis of plane perpendicular to l1
# x = Point(rand(3)) # random point in 3D
# p1 = Plane(x, x + b[:, 1], x + b[:, 2]) # random plane orthogonal to l1
# y = Point(rand(3))
# p2 = Plane(y, y + b[:, 1], y + b[:, 2]) # random plane orthogonal to l1, parallel to p1
# l2 = Line(y, y + b[:, 1] + b[:, 2]) # random line in 3D orthogonal to l1, parallel to p1 and p2
# p3 = Plane(Point(rand(3)), Point(rand(3)), Point(rand(3))) # random plane in 3D
# l3 = Line(Point(rand(3)), Point(rand(3))) # random line in 3D
# l4 = Line(x, x + (l1.points[2] - l1.points[1])) # line parallel to l1
# p4 = Plane(x, x + (l1.points[2] - l1.points[1]), Point(rand(3))) # random plane parallel to l1

l1 = Line(Point(0.0, 0, 0), Point(1, 0, 0)) # random line in 3D
# b = nullspace(transpose(l1.points[2] - l1.points[1])) # translational basis of plane perpendicular to l1
x = Point(0, 0, 1) # random point in 3D
p1 = Plane(Point(0.0, 0, 0), Point(0, 1, 0), Point(0, 0, 1)) # plane orthogonal to l1
y = Point(0, 0, 2)
p2 = Plane(Point(1.0, 0, 0), Point(1, 1, 0), Point(1, 0, 1)) # plane orthogonal to l1, parallel to p1
l2 = Line(Point(2.0, 0, 0), Point(2, 1, 0)) # line in 3D orthogonal to l1, parallel to p1 and p2
p3 = Plane(Point(rand(3)), Point(rand(3)), Point(rand(3))) # random plane in 3D
l3 = Line(Point(rand(3)), Point(rand(3))) # random line in 3D
l4 = Line(Point(0.0, 0, 1), Point(1, 0, 1)) # line parallel to l1
p4 = Plane(Point(0.0, 0, 2), Point(1, 0, 2), Point(0, 1, 2)) # random plane parallel to l1

pitch1 = rand() # random pitch for screw groups
pitch2 = rand()

I = IdentityGroup(3)
T_l1 = PrismaticGroup(l1)
T_l2 = PrismaticGroup(l2)
T_l3 = PrismaticGroup(l3)
T_l4 = PrismaticGroup(l4)
R_l1 = RevolutionGroup(l1)
R_l2 = RevolutionGroup(l2)
R_l3 = RevolutionGroup(l3)
R_l4 = RevolutionGroup(l4)
H_l1_p1 = ScrewGroup(l1, pitch1)
H_l2_p1 = ScrewGroup(l2, pitch1)
H_l3_p1 = ScrewGroup(l3, pitch1)
H_l4_p1 = ScrewGroup(l4, pitch1)
H_l1_p2 = ScrewGroup(l1, pitch2)
H_l2_p2 = ScrewGroup(l2, pitch2)
H_l3_p2 = ScrewGroup(l3, pitch2)
H_l4_p2 = ScrewGroup(l4, pitch2)
T_p1 = PlanarTranslationGroup(p1)
T_p2 = PlanarTranslationGroup(p2)
T_p3 = PlanarTranslationGroup(p3)
T_p4 = PlanarTranslationGroup(p4)
C_l1 = CylindricalGroup(l1)
C_l2 = CylindricalGroup(l2)
C_l3 = CylindricalGroup(l3)
C_l4 = CylindricalGroup(l4)
T = SpatialTranslationGroup()
G_p1 = PlanarSlidingGroup(p1)
G_p2 = PlanarSlidingGroup(p2)
G_p3 = PlanarSlidingGroup(p3)
G_p4 = PlanarSlidingGroup(p4)
S_x = SphericalRotationGroup(x)
S_y = SphericalRotationGroup(y)
Y_l1_p1 = TranslatingScrewGroup(l1, pitch1)
Y_l2_p1 = TranslatingScrewGroup(l2, pitch1)
Y_l3_p1 = TranslatingScrewGroup(l3, pitch1)
Y_l4_p1 = TranslatingScrewGroup(l4, pitch1)
Y_l1_p2 = TranslatingScrewGroup(l1, pitch2)
Y_l2_p2 = TranslatingScrewGroup(l2, pitch2)
Y_l3_p2 = TranslatingScrewGroup(l3, pitch2)
Y_l4_p2 = TranslatingScrewGroup(l4, pitch2)
X_l1 = TranslatingGimbalGroup(l1)
X_l2 = TranslatingGimbalGroup(l2)
X_l3 = TranslatingGimbalGroup(l3)
X_l4 = TranslatingGimbalGroup(l4)

@testset "subgroup relations among displacment groups" begin
    # tests for Prismatic group
    @test T_l1 == T_l4
    @test is_subgroup(T_l1, T_p4)
    @test !is_subgroup(T_l1, T_p3)
    @test is_subgroup(T_l1, C_l4)
    @test is_subgroup(T_l1, T)
    @test is_subgroup(T_l1, G_p4)
    @test !is_subgroup(T_l1, T_p3)
    @test !is_subgroup(T_l1, S_x)
    @test is_subgroup(T_l1, Y_l2_p1)
    @test !is_subgroup(T_l1, Y_l3_p2)
    @test is_subgroup(T_l1, X_l3)

    # tests for Revolution group
    @test !is_subgroup(R_l1, T_l3)
    @test !is_subgroup(R_l1, T_p3)
    @test is_subgroup(R_l1, C_l4)
    @test !is_subgroup(R_l1, T_l3)
    @test !is_subgroup(R_l1, T)
    @test is_subgroup(R_l1, G_p1)
    @test !is_subgroup(R_l1, G_p3)
    @test is_subgroup(R_l4, S_x)
    @test !is_subgroup(R_l1, Y_l3_p1)
    @test is_subgroup(R_l1, X_l4)
    @test !is_subgroup(R_l1, X_l3)

    # tests for Screw group
    @test !(H_l1_p1 == H_l1_p2)
    @test !(H_l1_p1 == H_l2_p1)
    @test !(H_l1_p1 == H_l4_p1)
    @test !is_subgroup(H_l1_p1, T_p3)
    @test !is_subgroup(H_l1_p1, T_l3)
    @test !is_subgroup(H_l1_p1, T_l3)
    @test is_subgroup(H_l1_p1, C_l4)
    @test !is_subgroup(H_l1_p1, C_l3)
    @test !is_subgroup(H_l1_p1, T)
    @test !is_subgroup(H_l1_p1, G_p3)
    @test !is_subgroup(H_l1_p1, S_x)
    @test is_subgroup(H_l1_p1, Y_l1_p1)
    @test !is_subgroup(H_l1_p1, Y_l3_p1)
    @test !is_subgroup(H_l1_p1, Y_l1_p2)
    @test is_subgroup(H_l1_p1, X_l4)

    # tests for Planar Translation group
    @test T_p1 == T_p2
    @test !is_subgroup(T_p1, C_l3)
    @test is_subgroup(T_p1, T)
    @test is_subgroup(T_p1, G_p2)
    @test !is_subgroup(T_p1, G_p3)
    @test !is_subgroup(T_p1, S_x)
    @test is_subgroup(T_p1, Y_l1_p1)
    @test !is_subgroup(T_p1, Y_l3_p1)
    @test is_subgroup(T_p1, X_l3)

    # tests for Cylindrical group
    @test !(C_l1 == C_l4)
    @test !is_subgroup(C_l1, T_p3)
    @test !is_subgroup(C_l1, C_l3)
    @test !is_subgroup(C_l1, T)
    @test !is_subgroup(C_l1, G_p3)
    @test !is_subgroup(C_l1, S_x)
    @test !is_subgroup(C_l1, Y_l3_p1)
    @test is_subgroup(C_l1, X_l4)
    @test !is_subgroup(C_l1, X_l3)

    # tests for Spatial Translation group
    @test !is_subgroup(T, T_p3)
    @test !is_subgroup(T, C_l3)
    @test !is_subgroup(T, G_p3)
    @test !is_subgroup(T, S_x)
    @test !is_subgroup(T, Y_l3_p1)
    @test is_subgroup(T, X_l3)

    # tests for Planar Sliding group
    @test G_p1 == G_p2
    @test !is_subgroup(G_p1, T_p3)
    @test !is_subgroup(G_p1, C_l3)
    @test !is_subgroup(G_p1, T)
    @test !is_subgroup(G_p1, G_p3)
    @test !is_subgroup(G_p1, S_x)
    @test !is_subgroup(G_p1, Y_l3_p1)
    @test is_subgroup(G_p1, X_l1)
    @test !is_subgroup(G_p1, X_l3)

    # tests for Translating Screw group
    @test !is_subgroup(Y_l1_p1, G_p3)
    @test !is_subgroup(Y_l1_p1, S_x)
    @test !(Y_l1_p1 == Y_l1_p2)
    @test !(Y_l1_p1 == Y_l2_p1)
    @test is_subgroup(Y_l1_p1, X_l4)
    @test !is_subgroup(Y_l1_p1, X_l3)

    # tests for Translating Gimbal group
    @test !is_subgroup(X_l1, X_l3)
end

@testset "intersection of displacement groups" begin

end
