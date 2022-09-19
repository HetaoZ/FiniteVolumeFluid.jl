
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

# ------------------------------------------------------------------------------
# solve

function advance!(f::Fluid, dt::Float64, t::Real)
    backup_w!(f)
    for rk = 1:f.scheme.rungekutta.order
        update_fluxes!(f)
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

function time_step!(f::Fluid{dim}) where dim
    smax = 0.
    smax = @sync @distributed (max) for id in f.mesh.domain_indices
        f.marker[id] > 0 ? maximum([smax; sound_speed(f.rho[id], f.p[id], f.parameters["gamma"]) .+ map(abs, f.u[:,id])]) : 0.
    end
    dt = minimum(f.mesh.d) / smax * f.parameters["CFL"]
    if isnan(dt) || dt == Inf
        println("smax = ", smax)
        error("time_step NaN")
    end
    return dt
end

function update_cells!(f::Fluid, rk::Int, dt::Float64; fluid_markers = (1,))

    coeff = f.scheme.rungekutta.coeffs[rk]

    @sync @distributed for id in f.mesh.domain_indices
        if f.marker[id] in fluid_markers
            # try
                w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * f.rhs[:, id] * dt                
                f.w[:, id] = w    
                rho, u, e, p = cons2prim(w, f.parameters["Gamma"])     
                f.rho[id] = rho
                f.u[:,id] = u
                f.e[id] = e
                f.p[id] = p

                # if id == CartesianIndex(3, 3, 3)
                #     println("w = " ,w)
                #     println("rho, u, e, p = ", (rho,u,e,p))
                #     println("f.w = ", f.w[:,id])
                # end
            # catch
            # if abs(sum(w) - 3.5) > 1e-8
            #     println((id,f.w[:, id],f.wb[:, id],f.rhs[:, id],dt))
            #     for axis = 1:3
            #         println("axis = ", axis)
            #         println(f.flux[axis][:,id])
            #         println(f.flux[axis][:,add_cartesian(id, axis, 1)])
            #         println()
            #     end
            #     error()
            # end
            # end
        end
    end
end

function update_fluxes!(f::Fluid{dim}; fluid_markers = (1,)) where dim

    stencil_width = f.scheme.reco_scheme.stencil_width
    total_stencil_width = stencil_width + 2
    half_stencil_width = ceil(Int, stencil_width/2)
    
    @sync @distributed for id in f.mesh.domain_indices
        if f.marker[id] in fluid_markers

            ws = zeros(Float64, dim+2, total_stencil_width)
            for axis in 1:dim
                for jj in 1:total_stencil_width
                    ws[:,jj] = f.w[:, add_cartesian(id, axis, jj-half_stencil_width-1) ]
                end
                    
                f.flux[axis][:,id] = flux!(ws[:,1:end-1], axis, f.scheme.reco_scheme, f.scheme.flux_scheme, f.parameters)

                f.flux[axis][:,add_cartesian(id, axis, 1)] = flux!(ws[:,2:end], axis, f.scheme.reco_scheme, f.scheme.flux_scheme, f.parameters)

            end
        end
    end    
end 

function update_rhs!(f::Fluid{dim}; fluid_markers = (1,)) where dim

    @sync @distributed for id in f.mesh.domain_indices
        if f.marker[id] in fluid_markers
            
            rhs = zeros(Float64, dim+2)

            for axis = 1:dim
                rhs += (f.flux[axis][:,id] - f.flux[axis][:,add_cartesian(id,axis,1)]) / f.mesh.d[axis]
            end

            
            f.rhs[:,id] = rhs # + f.source[:,id] # source item

            # if id == CartesianIndex(3, 3, 3)
            #     println("rhs = " ,rhs)
            # end
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
    return "none"
end