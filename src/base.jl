struct Block 
    point1::Vector{Float64}
    point2::Vector{Float64}
    function Block(point1, point2)
        return new(collect(Float64, point1), collect(Float64, point2))
    end
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
    domain_indices::SubArray
    d::Vector{Float64}
    coords::NTuple{3,Vector{Float64}}
end

function StructuredMesh(dim::Int = 3, point1::Vector{Float64} = [0,0,0], point2::Vector{Float64} = [1,1,1], ncells::Vector{Int} = [1,1,1], nbound::Int = 1)

    @assert dim ∈ (1,2,3)

    # 网格永远是三维的，但物理量是 dim 维的。
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

    indices = CartesianIndices((NX,NY,NZ))
    if dim == 3
        domain_indices = view(indices, nbound+1:nbound+nx, nbound+1:nbound+ny, nbound+1:nbound+nz)
    elseif dim == 2
        domain_indices = view(indices, nbound+1:nbound+nx, nbound+1:nbound+ny, 1)
    else
        domain_indices = view(indices, nbound+1:nbound+nx, 1, 1)
    end

    return StructuredMesh(dim, point1, point2, ncells, nbound, nx, ny, nz, NX, NY, NZ, indices, domain_indices, d, (x, y, z))
end

@inline function mesh_coords(mesh::StructuredMesh, id::CartesianIndex)
    return [mesh.coords[i][id[i]] for i = 1:mesh.dim]
end

@inline function lower_id(mesh::StructuredMesh, point::Tuple)
    return CartesianIndex(ntuple(axis -> max(1, ceil(Int, (point[axis] - mesh.coords[axis][1]) / mesh.d[axis])), 3))
end

@inline function higher_id(mesh::StructuredMesh, point::Tuple)
    id = lower_id(mesh, point)
    for axis = 1:mesh.dim
        add_cartesian(id, axis, 1)
    end
    return id
end

@inline function region_indices!(mesh::StructuredMesh, point1::Tuple, point2::Tuple; extension::Int = 0)
    l_id, h_id = lower_id(mesh, point1), higher_id(mesh, point2)
    return view(mesh.indices, ntuple(axis -> axis <= mesh.dim ? range(max(mesh.nbound+1,l_id[axis] - extension), min(mesh.ncells[axis]+mesh.nbound, h_id[axis] + extension), step = 1) : 1, length(l_id))...)
end

@inline function in_computational_domain(mesh::StructuredMesh, id::CartesianIndex)
    return betweeneq(mesh_coords(mesh, id), mesh.point1, mesh.point2)  
end

@inline function in_computational_domain(mesh::StructuredMesh, point::Vector{Float64})
    return betweeneq(point, mesh.point1, mesh.point2)  
end

abstract type AbstractRungeKutta end

struct RungeKutta{N} <: AbstractRungeKutta
    order::Int
    coeffs::NTuple{N, NTuple{3,Float64}}
end

function RungeKutta(order::Int)
    if order == 1
        coeffs = ((1.0, 0.0, 1.0),)
    elseif order == 2
        coeffs = ((1.0, 0.0, 1.0), 
        (0.5, 0.5, 0.5))
    elseif order == 3
        coeffs = ((1.0, 0.0, 1.0), 
        (0.75, 0.25, 0.25),
        (1/3, 2/3, 2/3))
    else
        error("undef order")
    end
    return RungeKutta{order}(order, coeffs)
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

function Muscl()
    return Muscl(2)
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

function Ausm()
    return Ausm(0.25, 0.75, 1.0)
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
    flux::Vector{SharedArray{Float64,4}}
    rhs::SharedArray{Float64,4}
    source::SharedArray{Float64,4} 
    last_marker::SharedArray{Int8,3}
    marker::SharedArray{Int8,3}
    boundaries::Vector{NTuple{2,AbstractBoundary}}
    blocks::Vector{Block}
    parameters::Dict
    initial_condition::Function
    rungekutta::AbstractRungeKutta
    reco_scheme::AbstractRecoScheme
    flux_scheme::AbstractFluxScheme
end

function Fluid(; dim::Int = 3, point1 = [0., 0., 0.], point2 = [1., 1., 1.], ncells = [1, 1, 1], nbound::Int = 1, parameters::Dict = Dict(), use_marker::Bool = true)

    default_para = Dict(
    "gamma"=>1.4,
    "Gamma"=>0.4,
    "inv_gamma_minus_one"=>1/(1.4-1), 
    "CFL"=>0.5
    )

    for key in keys(parameters)
        default_para[key] = parameters[key]
    end
    default_para["Gamma"] = default_para["gamma"] - 1
    default_para["inv_gamma_minus_one"] = 1/(default_para["gamma"]-1)

    point1 = collect(Float64, point1)
    point2 = collect(Float64, point2)
    ncells = collect(Int, ncells)

    mesh = StructuredMesh(dim, point1, point2, ncells, nbound)
    NX, NY, NZ = mesh.NX, mesh.NY, mesh.NZ

    fx = SharedArray(zeros(Float64, (2+dim, NX+1, NY, NZ)))
    fy = SharedArray(zeros(Float64, (1,1,1,1)))
    fz = SharedArray(zeros(Float64, (1,1,1,1)))
    if dim > 1
        fy = SharedArray(zeros(Float64, (2+dim, NX, NY+1, NZ)))
    end
    if dim > 2
        fz = SharedArray(zeros(Float64, (2+dim, NX, NY, NZ+1)))
    end

    if use_marker
        last_marker = SharedArray(ones(Int8, (NX, NY, NZ)))
        marker = SharedArray(ones(Int8, (NX, NY, NZ)))
    else
        last_marker = SharedArray(ones(Int8, (1,1,1)))
        marker = SharedArray(ones(Int8, (1,1,1)))
    end

    return Fluid(dim, mesh,
    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (dim, NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 
    SharedArray(zeros(Float64, (NX, NY, NZ))), 

    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # w
    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # wb

    [fx, fy, fz],

    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # rhs
    SharedArray(zeros(Float64, (dim+2, NX, NY, NZ))), # source

    last_marker,
    marker,

    [(FreeBoundary, FreeBoundary), (FreeBoundary, FreeBoundary), (FreeBoundary, FreeBoundary),],
    Block[],
    default_para,
    (x,t)->(0.,zeros(Float64,dim),0.,0.),
    RungeKutta(3),
    Muscl(2),
    Ausm(0.25, 1.0, 0.75))
end

function Base.show(f::Fluid)
    println("-- review of fluid --")
    println("# parameters")
    for k in keys(f.parameters)
        println("  ",k," : ",f.parameters[k])
    end
    println("# mesh")
    println("  dimension : ", f.mesh.dim)
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
    println("# schemes") 
    println("  reco scheme : ", f.reco_scheme)
    println("  flux scheme : ", f.flux_scheme)
    println("  iter scheme : Runge Kutta of O", f.rungekutta.order)
    println("# physical states")
    println("  rho ∈  ", [minimum(f.rho[f.mesh.domain_indices]), maximum(f.rho[f.mesh.domain_indices])])
    println("  e ∈  ", [minimum(f.e[f.mesh.domain_indices]), maximum(f.e[f.mesh.domain_indices])])
    println("  p ∈  ", [minimum(f.p[f.mesh.domain_indices]), maximum(f.p[f.mesh.domain_indices])])
    println()
end

function Base.copy(f::Fluid)
    
end