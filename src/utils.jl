@inline function expand(v::Vector{T}, new_length::Int, default_value::T) where T
    return [i>length(v) ? default_value : v[i] for i in 1:new_length]
end

@inline function betweeneq(a::Vector{T}, lo::Vector{T}, hi::Vector{T}) where T <: Real
    return all(lo .≤ a .≤ hi)
end

@inline function add_cartesian(id::CartesianIndex, axis::Int, n::Int)
    return CartesianIndex(ntuple(k -> k == axis ? k+n : id[k], length(id)))
end

@inline function reset_cartesian(id::CartesianIndex, axis::Int, i::Int)
    return CartesianIndex(ntuple(k -> k == axis ? i : id[k], length(id)))
end

# function copy_fluid!(f::Fluid)
    
#     f1 = @sync Fluid(realdim = f.realdim, point1 = f.point1, point2 = f.point2, ncells = f.ncells, nbound = f.nbound, para = f.parameters)

#     @sync @distributed for id in f.mesh.indices
#         f1.rho[id] = f.rho[id]
#         f1.e[id] = f.e[id]
#         f1.p[id] = f.p[id]
#         f1.marker[id] = f.marker[id]

#         f1.u[:,id] = f.u[:,id]
#         f1.w[:,id] = f.w[:,id]
#         f1.wb[:,id] = f.wb[:,id]
#         f1.rhs[:,id] = f.rhs[:,id]
#     end

#     f1.boundx = f.boundx
#     f1.boundy = f.boundy
#     f1.boundz = f.boundz

#     return f1
# end

# function get_point_LD!(x::Vector{Float64}, f::Fluid)
#     id = ones(Int, 3)
#     id[1:f.realdim] = [ceil(Int, (x[k]-f.point1[k])/f.d[k]-0.5+f.nbound) for k = 1:f.realdim]
#     return id
# end

# function get_point_RU!(x::Vector{Float64},f::Fluid)
#     iLD = get_point_LD!(x, f)
#     return iLD .+ 1
# end

# function get_point_located_cell!(x::Vector{Float64}, f::Fluid)
#     i = get_point_LD!(x, f)
#     ip = copy(i)
#     for k in 1:f.realdim
#         if x[k] > (i[k] - f.nbound) * f.d[k] + f.point1[k]
#             ip[k] += 1
#         end
#     end  
#     return ip
# end

# function check_marked_mass!(f::Fluid, marker::Int)
    
#     V = prod(f.d[1:f.realdim])

#     mass = @sync @distributed (+) for id in f.mesh.indices
#         local_mass = 0.
#         if f.marker[id] == marker
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_mass += f.rho[id] * V
#             end
#         end
#         local_mass
#     end

#     return mass
# end

# function check_marked_mass!(f::Fluid, markers::Vector{Int})
    
#     V = prod(f.d[1:f.realdim])

#     mass = @sync @distributed (+) for id in f.mesh.indices
#         local_mass = 0.
#         if f.marker[id] in markers
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_mass += f.rho[id] * V
#             end
#         end
#         local_mass
#     end

#     return mass
# end

# function check_mass!(f::Fluid)
    
#     V = prod(f.d[1:f.realdim])

#     mass = @sync @distributed (+) for id in f.mesh.indices
#         local_mass = 0.
#         if f.marker[id] != 0
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_mass += f.rho[id] * V
#             end
#         end
#         local_mass    # else
#     end

#     return mass
# end

# function check_field!(f::Fluid, field::Symbol)
    
#     V = prod(f.d[1:f.realdim])

#     w = @sync @distributed (+) for id in f.mesh.indices
#         local_w = 0.
#         if f.marker[id] != 0
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_w += getfield(f,field)[id] * V
#             end
#         end
#         local_w 
#     end

#     return w
# end

# function check_marked_field!(f::Fluid, field::Symbol, markers::Vector{Int})
    
#     V = prod(f.d[1:f.realdim])

#     w = @sync @distributed (+) for id in f.mesh.indices
#         local_w = 0.
#         if f.marker[id] in markers
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_w += getfield(f,field)[id] * V
#             end
#         end
#         local_w 
#     end

#     return w
# end

# function check_marked_w!(f::Fluid, k::Int, markers::Vector{Int})
    
#     V = prod(f.d[1:f.realdim])

#     w = @sync @distributed (+) for id in f.mesh.indices
#         local_w = 0.
#         if f.marker[id] in markers
#             if betweeneq([f.x[id[1]], f.y[id[2]], f.z[id[3]]], f.point1, f.point2)
#                 local_w += f.w[k, id] * V
#             end
#         end
#         local_w 
#     end

#     return w
# end

# function correct_cell_w(w, gamma, rho0, u0, e0)

#     rho, u, e, p = cons2prim(w, gamma)

#     if rho > 0 && e < 0
#         return prim2cons(rho, u, e0 * 1e-20)
#     end
#     if rho > 0 && p < 0
#         error("pressure(w, gamma) < 0. ")
#         return prim2cons(rho, u, e0 * 1e-20)
#     end   
    
#     return w
# end

