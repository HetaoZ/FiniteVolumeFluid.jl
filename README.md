# FiniteVolumeFluid

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HetaoZ.github.io/FiniteVolumeFluid.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HetaoZ.github.io/FiniteVolumeFluid.jl/dev)
[![Build Status](https://github.com/HetaoZ/FiniteVolumeFluid.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HetaoZ/FiniteVolumeFluid.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HetaoZ/FiniteVolumeFluid.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HetaoZ/FiniteVolumeFluid.jl)

A library

## Installation

```
] add https://github.com/HetaoZ/FiniteVolumeFluid.jl.git
```

## Examples
1d SOD problem in 3d space. The results saved into VTK(or VTR/VTI) files can be viewed with Paraview.

```
using FiniteVolumeFluid

const gamma = 1.4

function initial_condition(x, t)
    if t == 0.0
        if x[1] < 0.5
            rho, u, p = 1, [0,0,0], 1e5
        else
            rho, u, p = 10, [0,0,0], 1e6
        end
        e = p/(gamma - 1)/rho
        return rho, u, e, p
    end
end

f = Fluid(dim = 3, point1 = (0,0,0), point2 = (10,1,1), ncells = (40, 4, 4), nbound = 2);

# set parameters
set_parameters(f, "gamma"=>1.4, "CFL"=>0.5)

set_initial_condition(f, initial_condition)

set_boundaries(f, (FreeBoundary, FreeBoundary), (ReflBoundary, ReflBoundary), (ReflBoundary, ReflBoundary))

set_scheme(f, Muscl())
set_scheme(f, Ausm())
set_scheme(f, RungeKutta(3))

# show(f)

tspan = (0, 1)
nframe = 0
t = 0
N = 10000

save_mesh(f, "../../out/coord")
save_to_vtk(f, ("rho","u","e","p"), (:rho,:u,:e,:p), "../../out/fluid_"*string(N))


for frame = 1:10
    if t > tspan[2]
        break
    end
    
    dt = time_step!(f)
    t += dt

    advance!(f, dt, t)

    save_to_vtk(f, ("rho","u","e","p"), (:rho,:u,:e,:p), "../../out/fluid_"*string(N + frame))

    println((frame, t, dt))
end

(1, 0.0003340765523905305, 0.0003340765523905305)
(2, 0.000557406433502483, 0.00022332988111195253)
(3, 0.000760556131444011, 0.00020314969794152797)
(4, 0.0009523991844031206, 0.00019184305295910964)
(5, 0.001137923598610722, 0.0001855244142076014)
(6, 0.001319843262026125, 0.00018191966341540302)
(7, 0.0014997482894858719, 0.00017990502745974678)
(8, 0.0016786541677924235, 0.00017890587830655174)
(9, 0.0018572371249360232, 0.00017858295714359978)
(10, 0.0020359585932986112, 0.00017872146836258806)


```
