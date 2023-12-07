using PointToTriangle
using Test

@test all(1:1_000_000) do _
    a = rand(Vec{3})
    b = rand(Vec{3})
    c = rand(Vec{3})
    p = rand(Vec{3})
    tri1 = PointToTriangle.Triangle(a,b,c)
    tri2 = PointToTriangle.Triangle_3DMethod(a,b,c)
    PointToTriangle.vector(p, tri1) ≈ PointToTriangle.vector(p, tri2)
end
