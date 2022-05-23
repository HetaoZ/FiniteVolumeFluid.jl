function advance!(f::Fluid, dt::Float64, t::Real)
    backup_w!(f)
    for rk = 1:f.rungekutta.order
        update_rhs!(f)
        update_cells!(f, rk, dt)
        update_bounds!(f)
        initialize_fluid!(f, t) # Reinitialize the fluid for "everlasting sources" when t > 0
    end
end

function backup_w!(f::Fluid)
    @sync @distributed for id in eachindex(f.w)
        f.wb[id] = f.w[id]
    end     
end

function time_step!(f::Fluid)
    smax = 0.
    smax = @sync @distributed (max) for id in f.mesh.domain_indices
        f.marker[id] > 0 ? maximum([smax; sound_speed(f.rho[id], f.p[id], f.parameters["gamma"]) .+ map(abs, f.u[:,id])]) : 0.
    end
    dt = minimum(f.mesh.d[1:f.dim]) / smax * f.parameters["CFL"]
    if isnan(dt) || dt == Inf
        println("smax = ", smax)
        error("time_step NaN")
    end
    return dt
end

function update_cells!(f::Fluid, rk::Int, dt::Float64)

    coeff = f.rungekutta.coeffs[rk]

    @sync @distributed for id in f.mesh.domain_indices
        if f.marker[id] > 0
            # try
                w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * f.rhs[:, id] * dt                
                f.w[:, id] = w    
                rho, u, e, p = cons2prim(w, f.parameters["Gamma"])     
                f.rho[id] = rho
                f.u[:,id] = u
                f.e[id] = e
                f.p[id] = p
            # catch
            #     println((id,f.w[:, id],f.wb[:, id],f.rhs[:, id], dt))
            #     error()
            # end
        end
    end
end

function update_rhs!(f::Fluid)
    dim = f.dim
    @assert dim ∈ (1,2,3)
    stencil_width = f.reco_scheme.stencil_width
    half_stencil_width = Int((stencil_width+1)*0.5)

    @sync @distributed for id in f.mesh.domain_indices
        
        if f.marker[id] > 0

            point = mesh_coords(f.mesh, id)
            rhs = zeros(Float64, dim+2)
            ws = zeros(Float64, dim+2, stencil_width)

            for axis = 1:dim
                for jj in 1:stencil_width

                    pos_direction = jj > half_stencil_width
                    wall = where_is_block_wall(point, axis, pos_direction, f.blocks)

                    if wall == Missing
                        ws[:,jj] = copy(f.w[:, add_cartesian(id, axis, jj-half_stencil_width) ])
                    else
                        # when covered by blocks
                        iL, iR, λ = image_interpolant1d(f.mesh.coords[axis][id[axis] + jj-half_stencil_width], wall, f.mesh.coords[axis])
                        image_w = (1.0 - λ) * f.w[:, reset_cartesian(id, axis, iL)] + λ * f.w[:, reset_cartesian(id, axis, iR)]
                        image_w[1+axis] *= -1.0
                        ws[:,jj] = image_w
                    end
                end

                rhs += dflux!(ws, axis, f.reco_scheme, f.flux_scheme, f.parameters) / f.mesh.d[axis]

                # if isnan(rhs)
                #     println("id, ws = ", (id, ws))
                #     error("rhs")
                # end
            end
            
            f.rhs[:,id] = rhs + f.source[:,id]

        end
    end    
end 

@inline function image_interpolant1d(ghost_x::Float64, wall::Float64, mesh_x::Vector{Float64})
    return neighbors(image1d(ghost_x, wall), mesh_x)
end

@inline function image1d(ghost_x::Float64, wall::Float64)
    return 2.0 * wall - ghost_x
end

@inline function neighbors(x::Float64, mesh_x::Vector{Float64})
    iL = ceil(Int, (x - mesh_x[1])/(mesh_x[2] - mesh_x[1]))
    return iL, iL+1, (x-mesh_x[iL])/(mesh_x[iR] - mesh_x[iL])
end

"""
要求任意两个block的交集为空。
"""
@inline function where_is_block_wall(x::Vector{Float64}, axis::Int, pos_direction::Bool, blocks::Vector{Block})
    for block in blocks
        if in_block(x, block)
            return pos_direction ? block.point1[axis] : block.point2[axis]
        end
    end
    return Missing
end