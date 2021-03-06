function save_to_vtk(f::Fluid, datanames::T where T<:Tuple, fields::T where T<:Tuple, file_name::String; remove_boundary::Bool = true)

    vtkfile = create_vtkfile(f, file_name, remove_boundary)

    if remove_boundary
        for i in eachindex(datanames)
            if fields[i] == :u
                vtkfile["ux"] = getfield(f, :u)[1,f.mesh.domain_indices]
                if f.dim > 1
                    vtkfile["uy"] = getfield(f, :u)[2,f.mesh.domain_indices]
                    if f.dim > 2
                        vtkfile["uz"] = getfield(f, :u)[3,f.mesh.domain_indices]
                    end
                end
            else
                vtkfile[datanames[i]] = getfield(f, fields[i])[f.mesh.domain_indices]
            end
        end
    else
        for i in eachindex(datanames)
            if fields[i] == :u
                vtkfile["ux"] = getfield(f, :u)[1,:,:,:]
                if f.dim > 1
                    vtkfile["uy"] = getfield(f, :u)[2,:,:,:]
                    if f.dim > 2
                        vtkfile["uz"] = getfield(f, :u)[3,:,:,:]
                    end
                end
            else
                vtkfile[datanames[i]] = getfield(f, fields[i])
            end
        end
    end

    outfiles = vtk_save(vtkfile)
end

function create_vtkfile(f::Fluid, file_name::String, remove_boundary::Bool)
    if remove_boundary
        if f.dim == 1
            x = f.mesh.coords[1][f.mesh.nbound+1:end-f.mesh.nbound]
            file = vtk_grid(file_name, x)
        elseif f.dim == 2
            x = f.mesh.coords[1][f.mesh.nbound+1:end-f.mesh.nbound]
            y = f.mesh.coords[2][f.mesh.nbound+1:end-f.mesh.nbound]
            file = vtk_grid(file_name, x, y)
        else
            x = f.mesh.coords[1][f.mesh.nbound+1:end-f.mesh.nbound]
            y = f.mesh.coords[2][f.mesh.nbound+1:end-f.mesh.nbound]
            z = f.mesh.coords[3][f.mesh.nbound+1:end-f.mesh.nbound]
            file = vtk_grid(file_name, x, y, z)
        end
    else
        file = vtk_grid(file_name, f.mesh.coords[1:f.dim]...)
    end

    return file
end

function save_mesh(f::Fluid, file_name::String; remove_boundary::Bool = true)
    if remove_boundary 
        for k in 1:f.dim
            open(file_name*"_"*string(k)*".txt","w") do file
                writedlm(file, f.mesh.coords[k][f.mesh.nbound+1:end-f.mesh.nbound])
            end
        end 
    else   
        for k in 1:f.dim
            open(file_name*"_"*string(k)*".txt","w") do file
                writedlm(file, f.mesh.coords[k])
            end
        end  
    end
end