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
    f.reco_scheme = reco_scheme
end

function set_scheme(f::Fluid, flux_scheme::AbstractFluxScheme)
    f.flux_scheme = flux_scheme
end

function set_scheme(f::Fluid, rungekutta::AbstractRungeKutta)
    f.rungekutta = rungekutta
end