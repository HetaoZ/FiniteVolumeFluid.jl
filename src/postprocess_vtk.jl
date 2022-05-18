function save_to_vtk(f, datanames, fields, fname)
    vtkfile = create_vtkfile(f, fname)
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
    outfiles = vtk_save(vtkfile)
end

function create_vtkfile(f, fname)
    if f.dim == 1
        x = f.mesh.point1[1]-(f.mesh.nbound-0.5)*f.mesh.d[1]:f.mesh.d[1]:f.mesh.point2[1]+(f.mesh.nbound-0.5)*f.mesh.d[1]
        file = vtk_grid(fname, x)
    elseif f.dim == 2
        x = f.mesh.point1[1]-(f.mesh.nbound-0.5)*f.mesh.d[1]:f.mesh.d[1]:f.mesh.point2[1]+(f.mesh.nbound-0.5)*f.mesh.d[1]
        y = f.mesh.point1[2]-(f.mesh.nbound-0.5)*f.mesh.d[2]:f.mesh.d[2]:f.mesh.point2[2]+(f.mesh.nbound-0.5)*f.mesh.d[2]
        file = vtk_grid(fname, x, y)
    else
        x = f.mesh.point1[1]-(f.mesh.nbound-0.5)*f.mesh.d[1]:f.mesh.d[1]:f.mesh.point2[1]+(f.mesh.nbound-0.5)*f.mesh.d[1]
        y = f.mesh.point1[2]-(f.mesh.nbound-0.5)*f.mesh.d[2]:f.mesh.d[2]:f.mesh.point2[2]+(f.mesh.nbound-0.5)*f.mesh.d[2]
        z = f.mesh.point1[3]-(f.mesh.nbound-0.5)*f.mesh.d[3]:f.mesh.d[3]:f.mesh.point2[3]+(f.mesh.nbound-0.5)*f.mesh.d[3]
        file = vtk_grid(fname, x, y, z)
    end
    return file
end

function save_fluid_mesh(f::Fluid, fname)
    for k in 1:f.dim
        open(fname*"_"*string(k)*".txt","w") do file
            writedlm(file, [f.mesh.d[k]*(i-0.5-f.mesh.nbound)+f.mesh.point1[k]  for  i=1:f.mesh.ncells[k]+f.mesh.nbound*2])
        end
    end    
end