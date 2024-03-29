
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
        update_fluxes!(f, t)
        update_cells!(f, rk, dt)
        initialize_fluid!(f, t + dt) # Reinitialize the fluid for "everlasting sources" when t > 0 # update_bounds! 放在 initialize_fluid! 中
    end
end

function backup_w!(f::Fluid)
    @sync @distributed for id in f.grid.domain_indices
        f.wb[:,id] = f.w[:,id]
    end     
end

"合并 update_rhs! 和 update_cells!"
function update_cells!(f::Fluid{dim}, rk::Int, dt; fluid_markers = (1,)) where dim

    coeff = f.solver.rungekutta.coeffs[rk]

    @sync @distributed for id in f.grid.domain_indices
        if f.marker[id] in fluid_markers

                rhs = zeros(Float64, dim+2)
                for axis = 1:dim
                    rhs += (f.flux[axis][:,id] - f.flux[axis][:,add_cartesian(id,axis,1)]) / f.grid.d[axis]
                end
                
                # rhs += f.source[:,id] # source item

                w = coeff[1] * f.w[:, id] + coeff[2] * f.wb[:, id] + coeff[3] * rhs * dt   

                f.w[:, id] = w    
                rho, u, e, p = cons2prim(w, f.material)
                f.rho[id] = rho
                f.u[:,id] = u
                f.e[id] = e
                f.p[id] = p
        end
    end
end

function update_fluxes!(f::Fluid{dim}, t::Real; fluid_markers = (1,)) where dim

    stencil_width = f.solver.reconstruction.stencil_width
    total_stencil_width = stencil_width + 1
    half_stencil_width = ceil(Int, stencil_width/2)
    
    @sync @distributed for id in f.grid.domain_indices
        if f.marker[id] in fluid_markers

            ws = zeros(Float64, dim+2, total_stencil_width)

            for axis in 1:dim

                for jj in 1:total_stencil_width
                    w_id = add_cartesian(id, axis, jj-half_stencil_width-1)
                    
                    point = getcoordinates(f.grid, w_id) 
                    wall_vars = f.wall(point, t)
                    if wall_vars[1]
                        image_point = wall_vars[2]
                        image_w = local_fitting!(f, image_point)
                        n = normalize(image_point - point)
                        ws[:,jj] = image2ghost(image_w, n)
                    else                    
                        ws[:,jj] = f.w[:, w_id]
                    end
                end

                f.flux[axis][:,id] = compute_flux(ws, axis, f.solver, f.material)


                # if axis == 2 && id == CartesianIndex(25, 14)
                #     println("-- update_fluxes: 1")
                #     println("marker = ", f.marker[25, 11:17])
                #     println("ws = "); display(ws);println()
                # end
                
                # 对于靠近外边界或浸入边界的单元，需要额外计算另一侧界面通量
                next_id = add_cartesian(id, axis, 1)
                if f.marker[next_id] ∉ fluid_markers

                    ws = ws[:, [2:end ; 1]]

                    jj = total_stencil_width
                    w_id = add_cartesian(next_id, axis, total_stencil_width-half_stencil_width-1)
                    
                    point = getcoordinates(f.grid, w_id) 
                    wall_vars = f.wall(point, t)
                    if wall_vars[1]
                        image_point = wall_vars[2]
                        image_w = local_fitting!(f, image_point)
                        n = normalize(image_point - point)
                        ws[:,jj] = image2ghost(image_w, n)
                    else                    
                        ws[:,jj] = f.w[:, w_id]
                    end

                    f.flux[axis][:,next_id] = compute_flux(ws, axis, f.solver, f.material)

                end

            end
        end
    end    
end