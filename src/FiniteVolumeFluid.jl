
module FiniteVolumeFluid

using DelimitedFiles
using Printf
using Distributed
using LinearAlgebra
using Statistics
using WriteVTK
using SharedArrays


export 
    RectangularGrid,
    IdealGas,
    
    RungeKutta,
    Muscl,
    Weno,
    LaxFriedrichs,
    Ausm,
    FVSolver,
    
    FreeBoundary,
    ReflBoundary,
    Fluid
   
export 
    time_step,
    solve!,
    save_grid,
    save,
    behind_shock,
    pressure2e


include("base.jl")

include("utils/utils.jl")
include("utils/preprocess.jl")
include("utils/postprocess_vtk.jl")
include("utils/shock_dynamics.jl")

include("grid/grid.jl")
include("material/material.jl")

include("solver/solver.jl")
include("solver/flux.jl")
include("solver/reconstruction.jl")
include("solver/boundary.jl")
include("solver/local_fitting.jl")

end