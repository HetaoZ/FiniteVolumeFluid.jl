# ---------------------------------------------
# grid

abstract type AbstractGridPrototype end

mutable struct RectangularGrid{dim} <: AbstractGridPrototype
    start
    stop 
    nel::NTuple{dim,Int}
end

struct StructuredGrid{dim}
    prototype::AbstractGridPrototype
    indices::CartesianIndices
    domain_indices
    nbound::Int
    start
    stop
    nel
    N
    d
    x::NTuple{dim,Vector{Float64}}
end


# ---------------------------------------------------------
# material

abstract type AbstractMaterial end

struct IdealGas <: AbstractMaterial 
    γ # γ(air) = 1.4
    Γ # = γ - 1
    λ # = 1/(γ-1)
end

struct ViscousGas <: AbstractMaterial
    γ # γ(air) = 1.4
    Γ # = γ - 1
    λ # = 1/(γ-1)
    U₀
    L₀
    μ
    Re
end

# ---------------------------------------------------------
# solver

abstract type AbstractRungeKutta end

struct RungeKutta{N} <: AbstractRungeKutta
    order::Int
    coeffs::NTuple{N, NTuple{3,Float64}}
end

abstract type AbstractLimiter end

struct SuperbeeLimiter <: AbstractLimiter end
struct VanLeerMeanLimiter <: AbstractLimiter end
struct VanLeerLimiter <: AbstractLimiter end
struct VanAlbabaLimiter <: AbstractLimiter end
struct MinmodLimiter <: AbstractLimiter end

abstract type AbstractReconstruction end

struct Muscl <: AbstractReconstruction
    order::Int
    stencil_width::Int
    limiter::AbstractLimiter
end

struct Weno <: AbstractReconstruction
    order::Int
    stencil_width::Int
    smoothness_function::String
    eps
    a
    C
    p::Int
end

abstract type AbstractFlux end

struct LaxFriedrichs <: AbstractFlux end

struct Ausm <: AbstractFlux
    Kp::Float64
    Ku::Float64
    sigma::Float64
end

abstract type AbstractSolver end

"有限体积法"
struct FVSolver <: AbstractSolver
    CFL::Real
    rungekutta::AbstractRungeKutta
    reconstruction::AbstractReconstruction
    flux::AbstractFlux
end

# ---------------------------------------------------------
abstract type AbstractBoundary end

struct FlowBoundary <: AbstractBoundary
    coeff::Float64
end

const FreeBoundary = FlowBoundary(1.0)
const ReflBoundary = FlowBoundary(-1.0)


# ---------------------------------------------------------

mutable struct Fluid{dim} 
    grid::StructuredGrid{dim}
    material::AbstractMaterial
    solver::AbstractSolver
    
    initial_condition::Function # 初值条件
    wall::Function  # 允许自定义任意形状的 wall
    boundaries::Vector{NTuple{2,AbstractBoundary}} # 最外侧边界

    rho::SharedArray
    u::SharedArray
    e::SharedArray
    p::SharedArray
    w::SharedArray
    wb::SharedArray
    flux::NTuple{dim,SharedArray}
    rhs::SharedArray
    source::SharedArray
    marker::SharedArray
end