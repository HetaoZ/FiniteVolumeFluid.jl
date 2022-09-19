
module FiniteVolumeFluid

using DelimitedFiles
using Printf
using Distributed
using LinearAlgebra
using Statistics
using WriteVTK
using SharedArrays



export 
    Fluid,
    FreeBoundary,
    ReflBoundary,
    Muscl,
    Weno,
    LaxFriedrichs,
    Ausm,
    RungeKutta
export 
    set_initial_condition, 
    set_parameters,
    set_boundaries,
    set_scheme,
    insert_block,
    time_step!,
    advance!,
    save_mesh,
    save_to_vtk,
    mesh_coords,
    behind_shock


include("base.jl")

include("utils/physical.jl")
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

end