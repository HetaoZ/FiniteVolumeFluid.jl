"""
The function `initial_condition(x::Vector, t)` should return primitive variables `(rho, u::Vector, e, p)` or `Nothing`.
"""
function set_initial_condition(f::Fluid, initial_condition)
    f.initial_condition = initial_condition
    initialize_fluid!(f, 0.)
end

function initialize_fluid!(f::Fluid, t::Real)
    @sync @distributed for id in f.mesh.indices # 用完整的indices以确保边界外的来流条件也可以生效
        prim = f.initial_condition(mesh_coords(f.mesh, id), t)
        if typeof(prim) <: Tuple && length(prim) == 4
            f.rho[id], f.u[:,id], f.e[id], f.p[id] = prim
            f.w[:,id] = prim2cons(Float64(prim[1]), Float64.(prim[2]), Float64(prim[3]))
        # else
        #     error("Incorrect primary variables from initial_condition")
        end
    end
end

function set_parameters(f::Fluid, parameters::Pair...)
    for p in parameters
        f.parameters[p.first] = p.second
        if p.first == "gamma"
            f.parameters["Gamma"] = p.second - 1.0
            f.parameters["inv_gamma_minus_one"] = 1.0 / (p.second - 1.0)
        end
    end
end

function insert_block(f, point1::Vector, point2::Vector)
    block = Block(point1, point2)
    push!(f.blocks, block)
    f.initial_condition = (x,t) -> in_block(x, block) ? (0., zeros(Float64, size(u)), 0., 0.) : f.initial_condition(x,t)
end

function set_boundaries(f::Fluid, axis::Int, boundary_types::NTuple{2,AbstractBoundary})
    f.boundaries[axis] = boundary_types
end

function set_boundaries(f::Fluid, boundary_types::NTuple{2,AbstractBoundary}...)
    for axis in eachindex(boundary_types)
        f.boundaries[axis] = boundary_types[axis]
    end
end

function set_scheme(f::Fluid, reco_scheme::AbstractRecoScheme)
    f.scheme.reco_scheme = reco_scheme
end

function set_scheme(f::Fluid, flux_scheme::AbstractFluxScheme)
    f.scheme.flux_scheme = flux_scheme
end

function set_scheme(f::Fluid, rungekutta::AbstractRungeKutta)
    f.scheme.rungekutta = rungekutta
end

# -------------------------------------------------------------
function Fluid{dim}(grid_prototype::AbstractGridPrototype, scheme::Scheme) where dim

end


function Fluid(; parameters::Dict = Dict(), use_marker::Bool = true)

    point1, point2 = collect(prototype.start), collect(prototype.stop)
    ncells = prototype.nel

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

    mesh = StructuredGrid{dim}(prototype)
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
