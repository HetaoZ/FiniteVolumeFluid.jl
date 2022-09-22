function update_bounds!(f::Fluid{dim}) where dim

    @assert typeof(f.grid) == StructuredGrid{dim}
    
    nbound = f.grid.nbound
    nel = f.grid.nel

    @sync @distributed for id in f.grid.indices
        for axis in 1:dim
            
            if id[axis] < nbound + 1
                image_id = change_cartesian(id, axis, h->2*nbound+1-h)
                rho = f.rho[image_id]
                u = f.u[:, image_id]
                u[axis] *= f.boundaries[axis][1].coeff
                e = f.e[image_id]
                p = f.p[image_id]
                w = prim2cons(rho, u, e)

                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
            if id[axis] > nel[axis] + nbound
                image_id = change_cartesian(id, axis, h->2*(nel[axis]+nbound)+1-h)
                rho = f.rho[image_id]
                u = f.u[:, image_id]
                u[axis] *= f.boundaries[axis][2].coeff
                e = f.e[image_id]
                p = f.p[image_id]
                w = prim2cons(rho, u, e)

                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
        
        end
    end 
end

