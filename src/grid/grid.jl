
function StructureGrid(prototype::RectangularGrid{dim}) where dim

    point1, point2 = collect(prototype.start), collect(prototype.stop)

    @assert dim ∈ (1,2,3)

    d = (point2 - point1) ./ ncells
    @assert d[1] > 0 

    n = prototype.nel
    nbound = 3
    N = n .+ 2*nbound

    x = ntuple(axis -> [point1[axis] + d[axis] * (i-0.5-nbound) for i = 1:N[axis]], dim)

    indices = CartesianIndices(N)
    if dim == 3
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2], nbound+1:nbound+n[3])
    elseif dim == 2
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2])
    else
        domain_indices = view(indices, nbound+1:nbound+n[1])
    end

    return StructuredGrid{dim}(prototype, point1, point2, ncells, nbound, n, N, indices, domain_indices, d, x)
end

@inline function mesh_coords(mesh::StructuredGrid{dim}, id::CartesianIndex) where dim
    return [mesh.coords[i][id[i]] for i = 1:dim]
end

@inline function lower_id(mesh::StructuredGrid{dim}, point::Tuple) where dim
    return CartesianIndex(ntuple(axis -> max(1, ceil(Int, (point[axis] - mesh.coords[axis][1]) / mesh.d[axis])), dim))
end

@inline function higher_id(mesh::StructuredGrid{dim}, point::Tuple) where dim
    id = lower_id(mesh, point)
    for axis = 1:dim
        add_cartesian(id, axis, 1)
    end
    return id
end

@inline function region_indices!(mesh::StructuredGrid{dim}, point1::Tuple, point2::Tuple; extension::Int = 0) where dim
    l_id, h_id = lower_id(mesh, point1), higher_id(mesh, point2)
    return view(mesh.indices, ntuple(axis -> range(max(mesh.nbound+1,l_id[axis] - extension), min(mesh.ncells[axis]+mesh.nbound, h_id[axis] + extension), step = 1), dim)...)
end

@inline function in_computational_domain(mesh::StructuredGrid, id::CartesianIndex)
    return betweeneq(mesh_coords(mesh, id), mesh.point1, mesh.point2)  
end

@inline function in_computational_domain(mesh::StructuredGrid, point::Vector{Float64})
    return betweeneq(point, mesh.point1, mesh.point2)  
end

@inline function in_computational_domain(mesh::StructuredGrid, point::Tuple)
    return betweeneq(collect(point), mesh.point1, mesh.point2)  
end