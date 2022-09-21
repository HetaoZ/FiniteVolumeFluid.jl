using Distributed
@everywhere include("../src/FiniteVolumeFluid.jl")
@everywhere using .FiniteVolumeFluid
@everywhere using Printf, Dates, LinearAlgebra

# conditions

# 当条件不满足时，应返回 Nothing
@everywhere function initial_condition(x, t)

    idealgas = IdealGas(1.4)

    rho0, u0, p0 = 1.0, [0.,0.,0.], 1.e5

    
    if x[1] < 0.25
        # supersonic inflow
        rho, u, p = behind_shock(rho0, u0[1], p0, 1.2, idealgas, 1); u = [u, 0, 0] 
        # rho, u, p = rho0, u0, p0

        e = pressure2e(rho, p, idealgas)
        return rho, u, e, p
    else
        if t == 0.0
            e0 = pressure2e(rho0, p0, idealgas)
            return rho0, u0, e0, p0
        end
    end
end

"定义壁面和镜像反射点"

# 球面 wall
@everywhere function wall(P, t)
    center = [1, 0.5, 0.5]
    r = 0.2
    in_wall = norm(P - center) < r
    if in_wall
        refl_func = (x) -> 2*(center + r*normalize(x-center)) - x
        return in_wall, refl_func(P)
    else
        return in_wall, 0
    end
    # return false
end

# 平面 wall
# @everywhere function wall(P, t)
#     xw = 1
#     if P[1] > xw
#         refl_func = (x)-> [2 * xw - x[1], x[2], x[3]]
#         return true, refl_func(P)
#     else
#         return false, 0
#     end
# end

@everywhere function new_fluid(initial_condition, wall)

    # grid
    grid = CartesianGrid{3}((0,0,0), (2,1,1), 2 .* (20, 10, 10))

    # material
    idealgas = IdealGas(1.4)

    # solver
    solver = FVSolver(RungeKutta(3), Muscl(), LaxFriedrichs(); CFL = 0.5)

    # 外边界
    boundaries = [(FreeBoundary, FreeBoundary), 
    (ReflBoundary, ReflBoundary), 
    (ReflBoundary, ReflBoundary) ]

    return Fluid(grid, idealgas, solver, initial_condition, boundaries, wall)
end

f = new_fluid(initial_condition, wall)

# solve
frame = 0
FRAME = 0
time = 0.
N = 1000000
path = "/home/hetao/Projects/JuliaProj/out/shock_sphere_wall3d/"

save(f, path*"fluid_"*string(N+FRAME))

println("FRAME frame                      Date         Δt          t ")
@printf "%5i %5i   %s  %.3e  %.3e\n" FRAME frame Dates.now() 0 time

while frame < 10000 && time < 1
    
    global frame, time, FRAME
        
    dt = time_step(f)
    @time solve!(f, dt, time)

    frame += 1
    time += dt

    if frame%10 == 0

        FRAME += 1

        save(f, path*"fluid_"*string(N+FRAME))
        
        @printf "%5i %5i   %s  %.3e  %.3e\n" FRAME frame Dates.now() dt time
    end
end