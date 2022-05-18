struct Block 
    point1::Vector{Float64}
    point2::Vector{Float64}
end

function Block(point1, point2)
    return Block(collect(Float64, point1), collect(Float64, point2))
end

function in_block(point::Vector{Float64}, block::Block)
    return betweeneq(point, block.point1, block.point2)
end

struct StructuredMesh
    dim::Int # dim = 3
    point1::Vector{Float64}
    point2::Vector{Float64}
    ncells::Vector{Int}
    nbound::Int
    nx::Int
    ny::Int
    nz::Int
    NX::Int
    NY::Int
    NZ::Int
    indices::CartesianIndices
    d::Vector{Float64}
    coords::NTuple{3,Vector{Float64}}
end

function StructuredMesh(dim::Int = 3, point1::Vector{Float64} = [0,0,0], point2::Vector{Float64} = [1,1,1], ncells::Vector{Int} = [1,1,1], nbound::Int = 1)

    @assert dim ∈ (1,2,3)

    # 网格永远是三维的，但物理量是指定维度的。
    point1 = expand(point1, 3, 0)
    point2 = expand(point2, 3, 1)
    ncells = expand(ncells, 3, 1)

    d = (point2 - point1) ./ ncells
    @assert d[1] > 0 

    nx, ny, nz = Tuple(ncells)

    NX = ncells[1]+2*nbound
    x = [point1[1] + d[1] * (i-0.5-nbound) for i = 1:NX]
    
    if dim > 1
        NY = ncells[2]+2*nbound
        y = [point1[2] + d[2] * (i-0.5-nbound) for i = 1:NY]
    else
        NY = 1
        y = [point1[2]]
    end
    if dim > 2
        NZ = ncells[3]+2*nbound
        z = [point1[3] + d[3] * (i-0.5-nbound) for i = 1:NZ]
    else
        NZ = 1
        z = [point1[3]]
    end

    indices = CartesianIndices((NX, NY, NZ))

    return StructuredMesh(dim, point1, point2, ncells, nbound, nx, ny, nz, NX, NY, NZ, indices, d, (x, y, z))
end

@inline function mesh_coords(mesh::StructuredMesh, id::CartesianIndex)
    return [mesh.coords[i][id[i]] for i = 1:3]
end

@inline function in_computational_domain(mesh::StructuredMesh, id::CartesianIndex)
    return betweeneq(mesh_coords(mesh, id), mesh.point1, mesh.point2)  
end

abstract type AbstractRungeKutta end

struct RungeKutta <: AbstractRungeKutta
    order::Int
    coeffs::NTuple{order, NTuple{3,Float64}}
end

function RungeKutta(order::Int)
    if order == 3
        coeffs = ((1.0, 0.0, 1.0), 
        (0.75, 0.25, 0.25),
        (1/3, 2/3, 2/3))
    else
        error("undef order")
    end
    return RungeKutta(order, coeffs)
end

abstract type AbstractRecoScheme end

struct Muscl <: AbstractRecoScheme
    order::Int
    stencil_width::Int
end

function Muscl(order::Int)
    if order == 2
        stencil_width = 5
    else
        error("undef order")
    end
    return Muscl(order, stencil_width)
end

struct Weno <: AbstractRecoScheme
    order::Int
    stencil_width::Int
    smoothness_function::String
end

abstract type AbstractFluxScheme end

struct LaxFriedrichs <: AbstractFluxScheme end

struct Ausm <: AbstractFluxScheme
    Kp::Float64
    Ku::Float64
    sigma::Float64
end

abstract type AbstractBoundary end

struct FlowBoundary <: AbstractBoundary
    coeff::Float64
end

const FreeBoundary = FlowBoundary(1.0)
const ReflBoundary = FlowBoundary(-1.0)

mutable struct Fluid
    dim::Int
    mesh::StructuredMesh
    rho::SharedArray{Float64,3}
    u::SharedArray{Float64,4}
    e::SharedArray{Float64,3}
    p::SharedArray{Float64,3}
    w::SharedArray{Float64,4}
    wb::SharedArray{Float64,4}
    fx::SharedArray{Float64,4}
    fy::SharedArray{Float64,4}
    fz::SharedArray{Float64,4}
    rhs::SharedArray{Float64,4}
    source::SharedArray{Float64,4} 
    marker::SharedArray{Int8,3}
    boundaries::Vector{NTuple{2,AbstractBoundary}}
    blocks::Vector{Block}
    parameters::Dict
    initial_condition::Function
    rungekutta::AbstractRungeKutta
    reco_scheme::AbstractRecoScheme
    flux_scheme::AbstractFluxScheme
end

function Fluid(; dim::Int = 3, point1::Array = [0, 0, 0], point2::Array = [1, 1, 1], ncells::Vector{Int} = [1, 1, 1], nbound::Int = 1, parameters::Dict = Dict())

    default_para = Dict("rho0"=>0.,
    "u0"=>[0., 0., 0.],
    "e0"=>0.,
    "p0"=>0.,
    "gamma"=>1.4,
    "Gamma"=>0.4,
    "inv_gamma_minus_one"=>1/(1.4-1), 
    "mu"=>1.e-6, 
    "Pr"=>1.0, 
    "L0"=>1, 
    "U0"=>1, 
    "viscous"=>false,
    "CFL"=>0.1
    )

    for key in keys(parameters)
        default_para[key] = parameters[key]
    end
    default_para["Gamma"] = default_para["gamma"] - 1
    default_para["inv_gamma_minus_one"] = 1/(default_para["gamma"]-1)

    mesh = StructuredMesh(point1, point2, ncells, nbound)

    fx = SharedArray(zeros(Float64, (5, NX+1, NY, NZ)))
    fy = SharedArray(zeros(Float64, (1,1,1,1)))
    fz = SharedArray(zeros(Float64, (1,1,1,1)))
    if dim > 1
        fy = SharedArray(zeros(Float64, (5, NX, NY+1, NZ)))
    end
    if dim > 2
        fz = SharedArray(zeros(Float64, (5, NX, NY, NZ+1)))
    end

    return Fluid(dim, mesh,
    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (dim, NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 

    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # w
    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # wb

    fx,
    fy,
    fz,

    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # rhs
    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # source

    SharedArray(ones(Int8, (NX, NY, NZ))), # marker

    [(FreeBoundary, FreeBoundary), (FreeBoundary, FreeBoundary), (FreeBoundary, FreeBoundary),],
    Block[],
    default_para,
    (x,t)->(0.,zeros(Float64,dim),0.,0.),
    RungeKutta(3),
    Muscl(2),
    Ausm(0.25, 1.0, 0.75))
end

function review(f::Fluid)
    println("-- review of fluid --")
    println("# parameters")
    for k in keys(f.parameters)
        println("  ",k," : ",f.parameters[k])
    end
    println("# mesh")
    println("  dimension : ", f.dim)
    println("  domain : ", Tuple(f.mesh.point1)," -> ", Tuple(f.mesh.point2))
    println("  number of cells : ", length(f.rho))
    println("  size of cells : ", size(f.rho))
    println("# parallelism")
    println("  number of workers : ", nworkers())
    println("# boundaries")
    println("  axis 1 : ", f.boundaries[1])
    if f.dim > 1
        println("  axis 2 : ", f.boundaries[2])
    end
    if f.dim > 2
        println("  axis 3 : ", f.boundaries[3])
    end
    
    println("# physical states")
    println("  rho ∈  ", [minimum(f.rho), maximum(f.rho)])
    println("  e ∈  ", [minimum(f.e), maximum(f.e)])
    println("  p ∈  ", [minimum(f.p), maximum(f.p)])
end
