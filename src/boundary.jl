
function set_bounds!(f::Fluid{dim}, boundaries::NTuple{dim,NTuple{2,AbstractBoundary}}) where dim
    f.boundaries = boundaries
    update_bounds!(f)
end

function update_bounds!(f::Fluid{dim}) where dim
    @assert typeof(f.mesh) == StructuredMesh{dim}

    nbound = f.mesh.nbound
    ncells = f.mesh.ncells

    @sync @distributed for id in f.mesh.indices
        for axis in 1:dim
            
            if id[axis] < nbound + 1
                image_id = change_cartesian(id, axis, h->2*nbound+1-h)
                marker = f.marker[image_id]
                rho = f.rho[image_id]
                u = f.u[:, image_id]
                u[axis] *= f.boundaries[axis][1].coeff
                e = f.e[image_id]
                p = f.p[image_id]
                w = prim2cons(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
            if id[axis] > ncells[axis] + nbound
                image_id = change_cartesian(id, axis, h->2*(ncells[axis]+nbound)+1-h)
                marker = f.marker[image_id]
                rho = f.rho[image_id]
                u = f.u[:, image_id]
                u[axis] *= f.boundaries[axis][2].coeff
                e = f.e[image_id]
                p = f.p[image_id]
                w = prim2cons(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
        
        end
    end 
end