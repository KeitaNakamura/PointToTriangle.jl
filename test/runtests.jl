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
    project(p, a, b) = a + ((b-a)⋅(p-a))/((b-a)⋅(b-a))*(b-a)
    for Tri in (PointToTriangle.Triangle, PointToTriangle.Triangle_3DMethod)
        @test all(1:1_000_000) do _
            a = rand(Vec{3})
            b = rand(Vec{3})
            c = rand(Vec{3})
            p = rand(Vec{3})
            tri = Tri(a,b,c)
            # find minimum one by one
            d = Inf
            # vertices
            d = min(d, norm(a-p))
            d = min(d, norm(b-p))
            d = min(d, norm(c-p))
            # sides
            if 0 ≤ (project(p,a,b)-a)⋅normalize(b-a) ≤ norm(b-a)
                d = min(d, norm(project(p,a,b)-p))
            end
            if 0 ≤ (project(p,b,c)-b)⋅normalize(c-b) ≤ norm(c-b)
                d = min(d, norm(project(p,b,c)-p))
            end
            if 0 ≤ (project(p,c,a)-c)⋅normalize(a-c) ≤ norm(a-c)
                d = min(d, norm(project(p,c,a)-p))
            end
            # inside
            p′ = project(a, p, p - (b-a)×(c-a))
            if normalize((a-p′)×(b-p′)) ≈ normalize((b-p′)×(c-p′)) ≈ normalize((c-p′)×(a-p′))
                d = min(d, norm(p′-p))
            end
            norm(PointToTriangle.vector(p, tri)) ≈ d
        end
    end
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
