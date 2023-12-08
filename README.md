# PointToTriangle.jl

*Fast point-to-triangle distance computation in 3D*

This program is mostly based on the study by [Jones (1995)](https://www.researchgate.net/profile/Mark-Jones-11/publication/243787422_3D_Distance_from_a_Point_to_a_Triangle/links/54cb6c3b0cf2240c27e7da11/3D-Distance-from-a-Point-to-a-Triangle.pdf).

[![CI](https://github.com/KeitaNakamura/PointToTriangle.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/KeitaNakamura/PointToTriangle.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/KeitaNakamura/PointToTriangle.jl/graph/badge.svg?token=bxCXQulaR6)](https://codecov.io/gh/KeitaNakamura/PointToTriangle.jl)

## Installation

```jl
pkg> add https://github.com/KeitaNakamura/PointToTriangle.jl.git
```

## Usage

Using the `isbitstype` type of `AbstractVector`, such as `SVector` in [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl), is recommended for sake of performance.

```jl
julia> using PointToTriangle, StaticArrays

julia> a, b, c = rand(SVector{3}), rand(SVector{3}), rand(SVector{3});

julia> tri = PointToTriangle.Triangle(a, b, c);

julia> p = rand(SVector{3});

julia> d = PointToTriangle.vector(p, tri) # compute point-to-triangle distance vector
3-element SVector{3, Float64} with indices SOneTo(3):
  0.21019419472595702
 -0.2926549510617162
 -0.1669372773882311
```

### Visualizing results

#### Plots.jl

```jl
 julia> using Plots

 julia> plotlyjs()

 julia> plot(map(Tuple, [a,b,c,a]))

 julia> plot!(map(Tuple, [p,p+d]))
```

#### Makie.jl

```jl
 julia> using GLMakie

 julia> lines([a b c a])

 julia> lines!([p p+d])
```

<img width="415" src="https://github.com/KeitaNakamura/PointToTriangle.jl/assets/16015926/64818c53-fa8e-4bcb-81e1-29b5c88076cf">
