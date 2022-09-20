
function RungeKutta(order::Int = 3)
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

function FVSolver(rungekutta::AbstractRungeKutta, reconstruction::AbstractReconstruction, flux::AbstractFlux; CFL = 0.5)
    return FVSolver(CFL, rungekutta, reconstruction, flux)
end

# ------------------------------------------------------------------------------
# solve

function time_step(f::Fluid{dim}) where dim
    smax = 0.
    smax = @sync @distributed (max) for id in f.grid.domain_indices
        f.marker[id] > 0 ? maximum([smax; sound_speed(f.rho[id], f.p[id], f.material) .+ map(abs, f.u[:,id])])  : 0.
    end

    dt = minimum(f.grid.d) / smax * f.solver.CFL

    if isnan(dt) || dt == Inf
        println("smax = ", smax)
        error("time_step NaN")
    end
    return dt
end

function solve!(f::Fluid, dt, t)

    backup_w!(f)
    
    for rk = 1:f.solver.rungekutta.order
        @time update_fluxes!(f, t)
        update_rhs!(f)
        update_cells!(f, rk, dt)
        update_bounds!(f)
        initialize_fluid!(f, t) # Reinitialize the fluid for "everlasting sources" when t > 0
    end
end

function backup_w!(f::Fluid)
    @sync @distributed for id in f.grid.domain_indices
        f.wb[:,id] = f.w[:,id]
    end     
end

function update_cells!(f::Fluid, rk::Int, dt; fluid_markers = (1,))

    coeff = f.solver.rungekutta.coeffs[rk]

    @sync @distributed for id in f.grid.domain_indices
        if f.marker[id] in fluid_markers

                w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * f.rhs[:, id] * dt   

                f.w[:, id] = w    
                rho, u, e, p = cons2prim(w, f.material)
                f.rho[id] = rho
                f.u[:,id] = u
                f.e[id] = e
                f.p[id] = p

        end
    end
end 

function update_rhs!(f::Fluid{dim}; fluid_markers = (1,)) where dim

    @sync @distributed for id in f.grid.domain_indices
        if f.marker[id] in fluid_markers
            
            rhs = zeros(Float64, dim+2)

            for axis = 1:dim
                rhs += (f.flux[axis][:,id] - f.flux[axis][:,add_cartesian(id,axis,1)]) / f.grid.d[axis]
            end
            
            f.rhs[:,id] = rhs # + f.source[:,id] # source item

        end
    end    
end 

function update_fluxes!(f::Fluid{dim}, t::Real; fluid_markers = (1,)) where dim

    stencil_width = f.solver.reconstruction.stencil_width
    total_stencil_width = stencil_width + 2
    half_stencil_width = ceil(Int, stencil_width/2)
    
    @sync @distributed for id in f.grid.domain_indices
        if f.marker[id] in fluid_markers

            ws = zeros(Float64, dim+2, total_stencil_width)

            for axis in 1:dim

                # println("-- update_fluxes: 1")
                # @time                 
                for jj in 1:total_stencil_width

                    w_id = add_cartesian(id, axis, jj-half_stencil_width-1)
                    
                    # check wall function
                    # println("-- update_fluxes: 1.0")
                    # @time
                    point = getcoordinates(f.grid, w_id) 
                    wall_vars = f.wall(point, t)
                    if wall_vars[1]
                        image_point = wall_vars[2]
                        # println("-- update_fluxes: 1.1")
                        # @time 
                        image_w = local_fitting!(f, image_point)
                        n = normalize(image_point - point)
                        ws[:,jj] = image2ghost(image_w, n)
                    else                    
                        # println("-- update_fluxes: 1.2")
                        # @time 
                        ws[:,jj] = f.w[:, w_id]
                    end
                end
                
                # println("-- update_fluxes: 2")

                f.flux[axis][:,id] = compute_flux(ws[:,1:end-1], axis, f.solver, f.material)
                f.flux[axis][:,add_cartesian(id, axis, 1)] = compute_flux(ws[:,2:end], axis, f.solver, f.material)
                
            end
        end
    end    
end