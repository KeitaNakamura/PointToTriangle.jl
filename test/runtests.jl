using PointToTriangle
using Test

using PointToTriangle: Triangle, Triangle_3DMethod
using Tensorial
using StaticArrays

@testset "EdgeEquation" begin
    @test all(1:1_000) do _
        X, Y, dX, dY, x, y = rand(6) .* rand([-1,1], 6)
        E = PointToTriangle.EdgeEquation(X, Y, dX, dY)
        sign(E(Vec(x,y))) == sign((x-X)*dY - (y-Y)*dX)
    end
end

@testset "PointToTriangle.Triangle" begin
    for T in (Float32, Float64)
        for Tri in (Triangle, Triangle_3DMethod)
            # Array
            a = rand(T, 3)
            b = rand(T, 3)
            c = rand(T, 3)
            @test typeof(@inferred Tri(a,b,c)) == Tri{T}
            # StaticArrays
            a = rand(SVector{3,T})
            b = rand(SVector{3,T})
            c = rand(SVector{3,T})
            @test typeof(@inferred Tri(a,b,c)) == Tri{T}
            # Tensorial
            a = rand(Vec{3,T})
            b = rand(Vec{3,T})
            c = rand(Vec{3,T})
            @test typeof(@inferred Tri(a,b,c)) == Tri{T}
        end
    end
end

@testset "PointToTriangle.vector" begin
    # check inference
    for T in (Float32, Float64)
        for Tri in (Triangle, Triangle_3DMethod)
            a = rand(T, 3)
            b = rand(T, 3)
            c = rand(T, 3)
            tri = Tri(a,b,c)
            # check in various arrays
            p1 = rand(T, 3)
            p2 = rand(SVector{3,T})
            p3 = rand(Vec{3,T})
            @test typeof(@inferred PointToTriangle.vector(p1, tri)) == typeof(p1)
            @test typeof(@inferred PointToTriangle.vector(p2, tri)) == typeof(p2)
            @test typeof(@inferred PointToTriangle.vector(p3, tri)) == typeof(p3)
        end
    end

    # check computation
    project(p, a, b) = a + ((b-a)⋅(p-a))/((b-a)⋅(b-a))*(b-a)
    for Tri in (Triangle, Triangle_3DMethod)
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

    # compare two methods
    @test all(1:1_000_000) do _
        a = rand(Vec{3})
        b = rand(Vec{3})
        c = rand(Vec{3})
        p = rand(Vec{3})
        tri1 = Triangle(a,b,c)
        tri2 = Triangle_3DMethod(a,b,c)
        PointToTriangle.vector(p, tri1) ≈ PointToTriangle.vector(p, tri2)
    end
end
