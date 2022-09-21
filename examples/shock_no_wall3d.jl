include("../src/FiniteVolumeFluid.jl")
using .FiniteVolumeFluid
using Printf
using Dates

# grid
grid = CartesianGrid{3}((0,0,0), (2,1,1), (20, 10, 10))

# material
idealgas = IdealGas(1.4)

# solver
solver = FVSolver(RungeKutta(3), Weno(), Ausm(); CFL = 0.5)

# conditions

# 当条件不满足时，应返回 Nothing
function initial_condition(x, t)

    rho0, u0, p0 = 1.0, [0.,0.,0.], 1.e5

    # supersonic inflow
    if x[1] < 1
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

using LinearAlgebra

"定义壁面和镜像反射点"
function wall(P, t)
    center = [2, 0.5, 0.5]
    r = 0.5
    in_wall = norm(P - center) < r
    if in_wall
        refl_func = (x) -> 2*(center + r*normalize(x-center)) - x
        return in_wall, refl_func(P)
    else
        return in_wall, 0
    end
end

# 外边界
boundaries = [(FreeBoundary, FreeBoundary), 
              (ReflBoundary, ReflBoundary), 
              (ReflBoundary, ReflBoundary) ]


f = Fluid(grid, idealgas, solver, initial_condition, boundaries)

# solve
frame = 0
FRAME = 0
time = 0.
N = 1000000
path = "/home/hetao/Projects/JuliaProj/out/shock_no_wall3d/"

save(f, path*"fluid_"*string(N+FRAME))

println("FRAME frame                      Date         Δt          t ")
@printf "%5i %5i   %s  %.3e  %.3e\n" FRAME frame Dates.now() 0 time

while frame < 2 && time < 1
    
    global frame, time, FRAME
    
    dt = time_step(f)
    solve!(f, dt, time)

    frame += 1
    time += dt

    if frame%1 == 0

        FRAME += 1

        save(f, path*"fluid_"*string(N+FRAME))
        
        @printf "%5i %5i   %s  %.3e  %.3e\n" FRAME frame Dates.now() dt time
    end
end