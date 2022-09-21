
# -------------------------------------------------------------
function Fluid(grid_prototype::CartesianGrid{dim}, material::AbstractMaterial, solver::FVSolver, initial_condition::Function,  boundaries, wall::Function = (x,t)->false) where dim

    grid = StructuredGrid{dim}(grid_prototype, getnbound(solver.reconstruction))

    return Fluid(grid, material, solver, initial_condition, boundaries, wall)
end

function Fluid(grid::StructuredGrid{dim}, material::AbstractMaterial, solver::FVSolver, initial_condition::Function, boundaries, wall::Function = (x,t)->false) where dim

    N = grid.N

    rho = SharedArray(zeros(Float64, N))
    u = SharedArray(zeros(Float64, (dim, N...)))
    e = SharedArray(zeros(Float64, N))
    p = SharedArray(zeros(Float64, N))
    marker = SharedArray(ones(Int8, N))

    w = SharedArray(zeros(Float64, (dim+2, N...)))
    wb = SharedArray(zeros(Float64, (dim+2, N...)))
    rhs = SharedArray(zeros(Float64, (dim+2, N...))) # rhs

    flux_array_func = (i) -> SharedArray(zeros(Float64, (2+dim, add_along(N,i,1)...)))
    flux = ntuple(flux_array_func, dim)

    f = Fluid{dim}(grid, material, solver, initial_condition, wall, boundaries, rho, u, e, p, w, wb, flux, rhs, marker)
    initialize_fluid!(f, 0.)
    
    return f
end


function initialize_fluid!(f::Fluid{dim}, t) where dim

    zero_state = 0., zeros(Float64, dim), 0., 0., zeros(Float64, dim+2)

    @sync @distributed for id in f.grid.indices # 用完整的indices以确保边界外的来流条件也可以生效
        if f.wall(getcoordinates(f.grid, id), t)[1]
            f.rho[id], f.u[:,id], f.e[id], f.p[id], f.w[:,id] = zero_state
        else
            prim = f.initial_condition(getcoordinates(f.grid, id), t)

            # if id == CartesianIndex(23,3)
            #     println("-- initialize_fluid: 1")
            #     println("x = ",getcoordinates(f.grid, id))
            #     println("prim = ", prim, " t = ", t)
            #     println("func = ", f.initial_condition([1.0250000000000001, 0.025], 0.01))
            # end

            if typeof(prim) <: Tuple && length(prim) == 4
                f.rho[id], f.u[:,id], f.e[id], f.p[id] = prim
                f.w[:,id] = prim2cons(Float64(prim[1]), Float64.(prim[2]), Float64(prim[3]))
            end
        end
    end
end