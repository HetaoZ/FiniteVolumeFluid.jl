function save(f::Fluid{dim}, file_name::String; remove_boundary::Bool = true) where dim
    save(f, ("rho","u","e","p"), (:rho,:u,:e,:p), file_name, remove_boundary = remove_boundary)
end

function save(f::Fluid{dim}, datanames::T where T<:Tuple, fields::T where T<:Tuple, file_name::String; remove_boundary::Bool = true) where dim
    vtkfile = create_vtkfile(f, file_name, remove_boundary)

    if remove_boundary
        for i in eachindex(datanames)
            if fields[i] == :u
                vtkfile["ux"] = getfield(f, :u)[1,f.grid.domain_indices]
                if dim > 1
                    vtkfile["uy"] = getfield(f, :u)[2,f.grid.domain_indices]
                    if dim > 2
                        vtkfile["uz"] = getfield(f, :u)[3,f.grid.domain_indices]
                    end
                end
            else
                vtkfile[datanames[i]] = getfield(f, fields[i])[f.grid.domain_indices]
            end
        end
    else
        for i in eachindex(datanames)
            if fields[i] == :u
                vtkfile["ux"] = getfield(f, :u)[1,f.grid.indices]
                if dim > 1
                    vtkfile["uy"] = getfield(f, :u)[2,f.grid.indices]
                    if dim > 2
                        vtkfile["uz"] = getfield(f, :u)[3,f.grid.indices]
                    end
                end
            else
                vtkfile[datanames[i]] = getfield(f, fields[i])
            end
        end
    end

    vtk_save(vtkfile)
end

function create_vtkfile(f::Fluid{dim}, file_name::String, remove_boundary::Bool) where dim
    if remove_boundary
        if dim == 1
            x = f.grid.x[1][f.grid.nbound+1:end-f.grid.nbound]
            file = vtk_grid(file_name, x)
        elseif dim == 2
            x = f.grid.x[1][f.grid.nbound+1:end-f.grid.nbound]
            y = f.grid.x[2][f.grid.nbound+1:end-f.grid.nbound]
            file = vtk_grid(file_name, x, y)
        else
            x = f.grid.x[1][f.grid.nbound+1:end-f.grid.nbound]
            y = f.grid.x[2][f.grid.nbound+1:end-f.grid.nbound]
            z = f.grid.x[3][f.grid.nbound+1:end-f.grid.nbound]
            file = vtk_grid(file_name, x, y, z)
        end
    else
        file = vtk_grid(file_name, f.grid.x[1:dim]...)
    end

    return file
end

function save_grid(f::Fluid{dim}, file_name::String; remove_boundary::Bool = true) where dim
    if remove_boundary 
        for k in 1:dim
            open(file_name*"_"*string(k)*".txt","w") do file
                writedlm(file, f.grid.x[k][f.grid.nbound+1:end-f.grid.nbound])
            end
        end 
    else   
        for k in 1:dim
            open(file_name*"_"*string(k)*".txt","w") do file
                writedlm(file, f.grid.x[k])
            end
        end  
    end
end

# function Base.show(f::Fluid)
#     println("-- review of fluid --")
#     println("# parameters")
#     for k in keys(f.parameters)
#         println("  ",k," : ",f.parameters[k])
#     end
#     println("# grid")
#     println("  dimension : ", typeof(f.grid).var)
#     println("  domain : ", Tuple(f.grid.point1)," -> ", Tuple(f.grid.point2))
#     println("  number of cells : ", length(f.rho))
#     println("  size of cells : ", size(f.rho))
#     println("# parallelism")
#     println("  number of workers : ", nworkers())
#     println("# boundaries")
#     println("  axis 1 : ", f.boundaries[1])
#     if f.dim > 1
#         println("  axis 2 : ", f.boundaries[2])
#     end
#     if f.dim > 2
#         println("  axis 3 : ", f.boundaries[3])
#     end
#     println("# schemes") 
#     println("  reco scheme : ", f.scheme.reco_scheme)
#     println("  flux scheme : ", f.scheme.flux_scheme)
#     println("  iter scheme : Runge Kutta of O", f.scheme.rungekutta.order)
#     println("# physical states")
#     println("  rho ∈  ", [minimum(f.rho[f.grid.domain_indices]), maximum(f.rho[f.grid.domain_indices])])
#     println("  e ∈  ", [minimum(f.e[f.grid.domain_indices]), maximum(f.e[f.grid.domain_indices])])
#     println("  p ∈  ", [minimum(f.p[f.grid.domain_indices]), maximum(f.p[f.grid.domain_indices])])
#     println()
# end