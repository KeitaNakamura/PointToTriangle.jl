using PointToTriangle
using Test

@testset "EdgeEquation" begin
    @test all(1:1_000) do _
        X, Y, dX, dY, x, y = rand(6) .* rand([-1,1], 6)
        E = PointToTriangle.EdgeEquation(X, Y, dX, dY)
        sign(E(Vec(x,y))) == sign((x-X)*dY - (y-Y)*dX)
    end
end

@testset "PointToTriangle.vector" begin
    @test all(1:1_000_000) do _
        a = rand(Vec{3})
        b = rand(Vec{3})
        c = rand(Vec{3})
        p = rand(Vec{3})
        tri1 = PointToTriangle.Triangle(a,b,c)
        tri2 = PointToTriangle.Triangle_3DMethod(a,b,c)
        PointToTriangle.vector(p, tri1) ≈ PointToTriangle.vector(p, tri2)
    end
end
