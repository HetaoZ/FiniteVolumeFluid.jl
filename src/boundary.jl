
function set_bounds!(f, boundaries::Vector{NTuple{2,AbstractBoundary}})
    f.boundaries = boundaries
    update_bounds!(f)
end

function update_bounds!(f::Fluid)
    @assert typeof(f.mesh) == StructuredMesh

    nbound = f.mesh.nbound
    ncells = f.mesh.ncells

    @sync @distributed for id in f.mesh.indices
        i, j, k = id[1], id[2], id[3]
        if i < nbound + 1
            marker = f.marker[2*nbound+1-i, j, k]
            rho = f.rho[2*nbound+1-i, j, k]
            u = f.u[:, 2*nbound+1-i, j, k]
            u[1] *= f.boundaries[1][1].coeff
            e = f.e[2*nbound+1-i, j, k]
            p = f.p[2*nbound+1-i, j, k]
            w = prim2cons(rho, u, e)

            f.marker[id] = marker
            f.rho[id] = rho
            f.u[:, id] = u
            f.e[id] = e
            f.p[id] = p
            f.w[:, id] = w
        end
        if i > ncells[1] + nbound
            marker = f.marker[2*(ncells[1]+nbound)+1-i, j, k]
            rho = f.rho[2*(ncells[1]+nbound)+1-i, j, k]
            u = f.u[:, 2*(ncells[1]+nbound)+1-i, j, k]
            u[1] *= f.boundaries[1][2].coeff
            e = f.e[2*(ncells[1]+nbound)+1-i, j, k]
            p = f.p[2*(ncells[1]+nbound)+1-i, j, k]
            w = prim2cons(rho, u, e)

            f.marker[id] = marker
            f.rho[id] = rho
            f.u[:, id] = u
            f.e[id] = e
            f.p[id] = p
            f.w[:, id] = w
        end
        
        if f.dim > 1
            if j < nbound + 1
                marker = f.marker[i, 2*nbound+1-j, k]
                rho = f.rho[i, 2*nbound+1-j, k]
                u = f.u[:, i, 2*nbound+1-j, k]
                u[2] *= f.boundaries[2][1].coeff
                e = f.e[i, 2*nbound+1-j, k]
                p = f.p[i, 2*nbound+1-j, k]
                w = prim2cons(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w
            end
            if j > ncells[2] + nbound
                marker = f.marker[i, 2*(ncells[2]+nbound)+1-j, k]
                rho = f.rho[i, 2*(ncells[2]+nbound)+1-j, k]
                u = f.u[:, i, 2*(ncells[2]+nbound)+1-j, k]
                u[2] *= f.boundaries[2][2].coeff
                e = f.e[i, 2*(ncells[2]+nbound)+1-j, k]
                p = f.p[i, 2*(ncells[2]+nbound)+1-j, k]
                w = prim2cons(rho, u, e)

                f.marker[id] = marker
                f.rho[id] = rho
                f.u[:, id] = u
                f.e[id] = e
                f.p[id] = p
                f.w[:, id] = w                         
            end

            if f.dim > 2
                if k < nbound + 1
                    marker = f.marker[i, j, 2*nbound+1-k]
                    rho = f.rho[i, j, 2*nbound+1-k]
                    u = f.u[:, i, j, 2*nbound+1-k]
                    u[3] *= f.boundaries[3][1].coeff
                    e = f.e[i, j, 2*nbound+1-k]
                    p = f.p[i, j, 2*nbound+1-k]
                    w = prim2cons(rho, u, e)

                    f.marker[id] = marker
                    f.rho[id] = rho
                    f.u[:, id] = u
                    f.e[id] = e
                    f.p[id] = p
                    f.w[:, id] = w                     
                end
                if k > ncells[3] + nbound
                    marker = f.marker[i, j, 2*(f.mesh.ncells[3]+f.mesh.nbound)+1-k]
                    rho = f.rho[i, j, 2*(f.mesh.ncells[3]+f.mesh.nbound)+1-k]
                    u = f.u[:, i, j, 2*(f.mesh.ncells[3]+f.mesh.nbound)+1-k]
                    u[3] *= f.boundaries[3][2].coeff
                    e = f.e[i, j, 2*(f.mesh.ncells[3]+f.mesh.nbound)+1-k]
                    p = f.p[i, j, 2*(f.mesh.ncells[3]+f.mesh.nbound)+1-k]
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
end