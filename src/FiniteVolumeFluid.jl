
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
include("utils.jl")
include("physical.jl")
include("preprocess.jl")
include("boundary.jl")
include("solver.jl")
include("flux.jl")
include("reconstruction.jl")
include("postprocess_vtk.jl")
include("shock_dynamics.jl")

###
end