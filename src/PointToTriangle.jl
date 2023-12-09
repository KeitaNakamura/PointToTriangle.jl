module PointToTriangle

using Tensorial

#############
# 2D method #
#############

struct EdgeEquation{T}
    n::Vec{2,T}
    Xn::T
end

function EdgeEquation(x::T, y::T, dx::T, dy::T) where {T}
    X = Vec(x,y)
    n = Vec(dy,-dx)
    EdgeEquation(n, X⋅n)
end

# negative: left, positive: right, zero: on the line
@inline function (E::EdgeEquation{T})(x::Vec{2,T}) where {T}
    n, Xn = E.n, E.Xn
    x⋅n - Xn
end

struct Side{T}
    # edge equations
    E_ab::EdgeEquation{T}
    E_aa′::EdgeEquation{T}
    E_bb′::EdgeEquation{T}
    # precomputed for projection
    c::Vec{2,T}
    v::Vec{2,T}
end

function Side(a::Vec{2,T}, b::Vec{2,T}) where {T}
    c = normalize(b - a)
    n = Vec(-c[2],c[1])
    E_ab = EdgeEquation(a[1], a[2], c[1], c[2])
    E_aa′ = EdgeEquation(a[1], a[2], n[1], n[2])
    E_bb′ = EdgeEquation(b[1], b[2], n[1], n[2])
    v = a - (a⋅c)*c
    Side(E_ab, E_aa′, E_bb′, c, v)
end

struct Triangle{T}
    P₁::Vec{3,T}
    P₂::Vec{3,T}
    P₃::Vec{3,T}
    # precomputed
    t::Vec{3,T}
    R::Mat{3,3,T,9}
    R⁻¹::Mat{3,3,T,9}
    side1::Side{T}
    side2::Side{T}
    side3::Side{T}
end

function Triangle(P₁::AbstractVector{T}, P₂::AbstractVector{T}, P₃::AbstractVector{T}) where {T}
    @assert length(P₁) == length(P₂) == length(P₃) == 3
    Triangle(Vec{3,T}(P₁), Vec{3,T}(P₂), Vec{3,T}(P₃))
end
function Triangle(P₁::Vec{3,T}, P₂::Vec{3,T}, P₃::Vec{3,T}) where {T}
    P₁′ = P₁ - P₁
    P₂′ = P₂ - P₁
    P₃′ = P₃ - P₁
    R₁ = rotmat(normalize(P₂′) => Vec{3,T}(0,0,1))
    P₁′ = rotate(P₁′, R₁)
    P₂′ = rotate(P₂′, R₁)
    P₃′ = rotate(P₃′, R₁)
    R₂ = rotmatz(atan(P₃′[1], P₃′[2]))
    P₁′ = rotate(P₁′, R₂)
    P₂′ = rotate(P₂′, R₂)
    P₃′ = rotate(P₃′, R₂)
    R = R₂ ⋅ R₁
    side1 = Side(@Tensor(P₁′[2:3]), @Tensor(P₂′[2:3]))
    side2 = Side(@Tensor(P₂′[2:3]), @Tensor(P₃′[2:3]))
    side3 = Side(@Tensor(P₃′[2:3]), @Tensor(P₁′[2:3]))
    Triangle(P₁, P₂, P₃, -P₁, R, inv(R), side1, side2, side3)
end

@inline function vector(P₀::AbstractVector{T}, tri::Triangle{T}) where {T}
    @assert length(P₀) == 3
    convert(typeof(P₀), vector(Vec{3,T}(P₀), tri))
end
@inline function vector(P₀::Vec{3,T}, tri::Triangle{T}) where {T}
    side1, side2, side3 = tri.side1, tri.side2, tri.side3
    P₀′ = @Tensor transform(P₀, tri)[2:3]
    # side1
    E₁_ab, E₁_aa′, E₁_bb′ = side1.E_ab(P₀′), side1.E_aa′(P₀′), side1.E_bb′(P₀′)
    E₁_ab<0 && E₁_aa′>0 && E₁_bb′<0 && return transform_inv(_project(P₀′, side1), tri)-P₀
    # side2
    E₂_ab, E₂_aa′, E₂_bb′ = side2.E_ab(P₀′), side2.E_aa′(P₀′), side2.E_bb′(P₀′)
    E₂_ab<0 && E₂_aa′>0 && E₂_bb′<0 && return transform_inv(_project(P₀′, side2), tri)-P₀
    # side3
    E₃_ab, E₃_aa′, E₃_bb′ = side3.E_ab(P₀′), side3.E_aa′(P₀′), side3.E_bb′(P₀′)
    E₃_ab<0 && E₃_aa′>0 && E₃_bb′<0 && return transform_inv(_project(P₀′, side3), tri)-P₀
    # vertex1
    E₁_aa′<0 && E₃_bb′>0 && return tri.P₁-P₀
    # vertex2
    E₂_aa′<0 && E₁_bb′>0 && return tri.P₂-P₀
    # vertex3
    E₃_aa′<0 && E₂_bb′>0 && return tri.P₃-P₀
    # inside
    transform_inv(P₀′, tri)-P₀
end

@inline _project(P₀′::Vec{2,T}, side::Side{T}) where {T} = side.v + (side.c⋅P₀′)*side.c
@inline transform(a::Vec{3,T}, tri::Triangle{T}) where {T} = tri.R ⋅ (a + tri.t)
@inline transform_inv(a::Vec{3,T}, tri::Triangle{T}) where {T} = tri.R⁻¹ ⋅ a - tri.t
@inline transform_inv(a::Vec{2,T}, tri::Triangle{T}) where {T} = transform_inv(vcat(0,a), tri)

################################
# 3D method (testing purposes) #
################################

struct Triangle_3DMethod{T}
    P₁::Vec{3,T}
    P₂::Vec{3,T}
    P₃::Vec{3,T}
    # precomputed
    Nₚ::Vec{3,T} # normalized
    V₁::Vec{3,T}
    V₂::Vec{3,T}
    V₃::Vec{3,T}
end

function Triangle_3DMethod(P₁::AbstractVector{T}, P₂::AbstractVector{T}, P₃::AbstractVector{T}) where {T}
    @assert length(P₁) == length(P₂) == length(P₃) == 3
    Triangle_3DMethod(Vec{3,T}(P₁), Vec{3,T}(P₂), Vec{3,T}(P₃))
end
function Triangle_3DMethod(P₁::Vec{3,T}, P₂::Vec{3,T}, P₃::Vec{3,T}) where {T}
    Nₚ = normalize((P₂-P₁) × (P₃-P₁))
    V₁ = Nₚ × (normalize(P₁-P₂) + normalize(P₁-P₃))
    V₂ = Nₚ × (normalize(P₂-P₃) + normalize(P₂-P₁))
    V₃ = Nₚ × (normalize(P₃-P₁) + normalize(P₃-P₂))
    Triangle_3DMethod(P₁, P₂, P₃, Nₚ, V₁, V₂, V₃)
end

@inline function vector(P₀::AbstractVector{T}, tri::Triangle_3DMethod{T}) where {T}
    @assert length(P₀) == 3
    convert(typeof(P₀), vector(Vec{3,T}(P₀), tri))
end
@inline function vector(P₀::Vec{3,T}, tri::Triangle_3DMethod{T}) where {T}
    P₁, P₂, P₃, Nₚ = tri.P₁, tri.P₂, tri.P₃, tri.Nₚ
    P₀′ = P₀ - ((P₀-P₁)⋅Nₚ) * Nₚ
    pos = _position(P₀′, tri)
    pos == 0 && return P₀′-P₀
    pos == 1 && return _point_to_side(P₀′, P₀, P₁, P₂)
    pos == 2 && return _point_to_side(P₀′, P₀, P₂, P₃)
    pos == 3 && return _point_to_side(P₀′, P₀, P₃, P₁)
    error("unreachable")
end

# 0: inside, 1: P₁P₂, 2: P₂P₃, 3: P₃P₁
@inline function _position(P₀′::Vec{3,T}, tri::Triangle_3DMethod{T}) where {T}
    P₁, P₂, P₃, Nₚ, V₁, V₂, V₃ = tri.P₁, tri.P₂, tri.P₃, tri.Nₚ, tri.V₁, tri.V₂, tri.V₃
    P₁P₀′ = P₀′ - P₁
    P₂P₀′ = P₀′ - P₂
    P₃P₀′ = P₀′ - P₃
    f₁ = V₁ ⋅ P₁P₀′
    f₂ = V₂ ⋅ P₂P₀′
    f₃ = V₃ ⋅ P₃P₀′
    f₁ > 0 && f₂ ≤ 0 && return ifelse((P₁P₀′×P₂P₀′)⋅Nₚ ≥ 0, 0, 1)
    f₂ > 0 && f₃ ≤ 0 && return ifelse((P₂P₀′×P₃P₀′)⋅Nₚ ≥ 0, 0, 2)
    f₃ > 0 && f₁ ≤ 0 && return ifelse((P₃P₀′×P₁P₀′)⋅Nₚ ≥ 0, 0, 3)
    error("unreachable")
end

@inline function _point_to_side(P₀′::Vec{3,T}, P₀::Vec{3,T}, P₁::Vec{3,T}, P₂::Vec{3,T}) where {T}
    P₁P₂ = P₂ - P₁
    P₀′P₁ = P₁ - P₀′
    R = ((P₂-P₀′) × P₀′P₁) × P₁P₂
    P₀′′ = -P₀′P₁ + ((P₀′P₁⋅R)/(R⋅R)) * R
    t = (P₀′′⋅P₁P₂) / norm2(P₁P₂)
    0 ≤ t ≤ 1 && return (P₁+t*P₁P₂)-P₀
    t < 0 ? P₁-P₀ : P₂-P₀
end
@inline norm2(x) = dot(x,x)

end # module PointToTriangle
