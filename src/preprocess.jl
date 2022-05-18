"""
The function `initial_condition(x::Vector, t)` returns primitive variables (rho, u::Vector, e, p).
"""
function initialize_fluid!(f::Fluid, t::Float64)
    @sync @distributed for id in f.mesh.indices
        rho, u, e, p = f.initial_condition(mesh_coords(f.mesh, id), t)
        f.rho[id], f.u[:,id], f.e[id], f.p[id] = rho, u, e, p
        f.w[:,id] = prim2cons(rho, u, e)
    end
end

function insert_block!(f, point1::Vector, point2::Vector)
    block = Block(point1, point2)
    push!(f.blocks, block)
    f.initial_condition = (x,t) -> in_block(x, block) ? (0., zeros(Float64, size(u)), 0., 0.) : f.initial_condition(x,t)
end