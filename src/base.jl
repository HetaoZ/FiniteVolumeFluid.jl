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

struct StructuredMesh{dim}
    point1::Vector{Float64}
    point2::Vector{Float64}
    ncells::Vector{Int}
    nbound::Int
    n::NTuple{dim,Int}
    N::NTuple{dim,Int}
    indices::CartesianIndices
    domain_indices::SubArray
    d::Vector{Float64}
    coords::NTuple{dim,Vector{Float64}}
end

function StructuredMesh(dim::Int = 3, point1::Vector{Float64} = [0,0,0], point2::Vector{Float64} = [1,1,1], ncells::Vector{Int} = [1,1,1], nbound::Int = 1)

    @assert dim ∈ (1,2,3)

    d = (point2 - point1) ./ ncells
    @assert d[1] > 0 

    n = Tuple(ncells)
    N = n .+ 2*nbound

    x = ntuple(axis -> [point1[axis] + d[axis] * (i-0.5-nbound) for i = 1:N[axis]], dim)

    indices = CartesianIndices(N)
    if dim == 3
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2], nbound+1:nbound+n[3])
    elseif dim == 2
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2])
    else
        domain_indices = view(indices, nbound+1:nbound+n[1])
    end

    return StructuredMesh{dim}(point1, point2, ncells, nbound, n, N, indices, domain_indices, d, x)
end

@inline function mesh_coords(mesh::StructuredMesh{dim}, id::CartesianIndex) where dim
    return [mesh.coords[i][id[i]] for i = 1:dim]
end

@inline function lower_id(mesh::StructuredMesh{dim}, point::Tuple) where dim
    return CartesianIndex(ntuple(axis -> max(1, ceil(Int, (point[axis] - mesh.coords[axis][1]) / mesh.d[axis])), dim))
end

@inline function higher_id(mesh::StructuredMesh{dim}, point::Tuple) where dim
    id = lower_id(mesh, point)
    for axis = 1:dim
        add_cartesian(id, axis, 1)
    end
    return id
end

@inline function region_indices!(mesh::StructuredMesh{dim}, point1::Tuple, point2::Tuple; extension::Int = 0) where dim
    l_id, h_id = lower_id(mesh, point1), higher_id(mesh, point2)
    return view(mesh.indices, ntuple(axis -> range(max(mesh.nbound+1,l_id[axis] - extension), min(mesh.ncells[axis]+mesh.nbound, h_id[axis] + extension), step = 1), dim)...)
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

function Muscl(order::Int = 2)
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
    eps
    a
    C
    p::Int
end

function Weno(order::Int = 5)
    if order == 5
        stencil_width = 5
    else
        error("undef order")
    end
    return Weno(order, stencil_width, "JS", 
    1.e-10,
    ([1/3 -7/6 11/6],
    [-1/6 5/6 1/3],
    [1/3 5/6 -1/6],
    [11/6 -7/6 1/3]),
    (0.1, 0.6, 0.3),
    2  # recommended p = r
    )
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

mutable struct Scheme
    rungekutta::AbstractRungeKutta
    reco_scheme::AbstractRecoScheme
    flux_scheme::AbstractFluxScheme
end

mutable struct Fluid{dim} 
    mesh::StructuredMesh{dim}
    rho::SharedArray{Float64,dim}
    u::SharedArray
    e::SharedArray{Float64,dim}
    p::SharedArray{Float64,dim}
    w::SharedArray
    wb::SharedArray
    flux::NTuple{dim,SharedArray}
    rhs::SharedArray
    source::SharedArray
    # last_marker::SharedArray{Int8,3}
    marker::SharedArray{Int8,dim}
    boundaries::Vector{NTuple{2,AbstractBoundary}}
    blocks::Vector{Block}
    parameters::Dict
    initial_condition::Function
    scheme::Scheme
end

function add_along(v::Tuple, axis::Int, n::Int)
    v1 = ntuple(i -> i==axis ? v[i]+n : v[i], length(v))
    return v1
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
    N = mesh.N

    flux = ntuple(i -> SharedArray(zeros(Float64, (2+dim, add_along(N,i,1)...))), dim)

    if use_marker
        marker = SharedArray(ones(Int8, N))
    else
        marker = SharedArray(ones(Int8, (1,)))
    end

    return Fluid{dim}(mesh,

    SharedArray(zeros(Float64, N)), 
    SharedArray(zeros(Float64, (dim, N...))), 
    SharedArray(zeros(Float64, N)), 
    SharedArray(zeros(Float64, N)), 

    SharedArray(zeros(Float64, (dim+2, N...))), # w
    SharedArray(zeros(Float64, (dim+2, N...))), # wb

    flux,

    SharedArray(zeros(Float64, (dim+2, N...))), # rhs
    SharedArray(zeros(Float64, (dim+2, N...))), # source

    marker,

    [(FreeBoundary, FreeBoundary) for i = 1:dim],
    Block[],
    default_para,

    (x,t)->(0.,zeros(Float64,dim),0.,0.),

    Scheme(RungeKutta(3),
    Muscl(2),
    Ausm(0.25, 1.0, 0.75)))
end

function Base.show(f::Fluid)
    println("-- review of fluid --")
    println("# parameters")
    for k in keys(f.parameters)
        println("  ",k," : ",f.parameters[k])
    end
    println("# mesh")
    println("  dimension : ", typeof(f.mesh).var)
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
    println("  reco scheme : ", f.scheme.reco_scheme)
    println("  flux scheme : ", f.scheme.flux_scheme)
    println("  iter scheme : Runge Kutta of O", f.scheme.rungekutta.order)
    println("# physical states")
    println("  rho ∈  ", [minimum(f.rho[f.mesh.domain_indices]), maximum(f.rho[f.mesh.domain_indices])])
    println("  e ∈  ", [minimum(f.e[f.mesh.domain_indices]), maximum(f.e[f.mesh.domain_indices])])
    println("  p ∈  ", [minimum(f.p[f.mesh.domain_indices]), maximum(f.p[f.mesh.domain_indices])])
    println()
end

function Base.copy(f::Fluid)
    
end