
function StructuredGrid{dim}(prototype::RectangularGrid{dim}, nbound) where dim
    @assert dim ∈ (1,2,3)
    start, stop = collect(prototype.start), collect(prototype.stop)
    nel = collect(prototype.nel)
    n = nel
    d = (stop - start) ./ nel
    @assert all([d[i] > 0 for i in 1:dim])
    N = Tuple(nel .+ 2*nbound)
    

    x = ntuple(axis -> [start[axis] + d[axis] * (i-0.5-nbound) for i = 1:N[axis]], dim)

    indices = CartesianIndices(N)
 
    if dim == 3
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2], nbound+1:nbound+n[3])
    elseif dim == 2
        domain_indices = view(indices, nbound+1:nbound+n[1], nbound+1:nbound+n[2])
    else
        domain_indices = view(indices, nbound+1:nbound+n[1])
    end

    return StructuredGrid{dim}(prototype, indices, domain_indices, nbound, start, stop, nel, N,  d, x)
end


@inline function getcoordinates(grid::StructuredGrid{dim}, id::CartesianIndex) where dim
    return [grid.x[i][id[i]] for i = 1:dim]
end

@inline function getcoordinates_in_region(grid::StructuredGrid{dim}, rindices) where dim
    return [grid.x[i][id[i]] for i in 1:dim, id in rindices]
end

@inline function lower_id(grid::StructuredGrid{dim}, point) where dim
    return CartesianIndex(ntuple(axis -> max(1, ceil(Int, (point[axis] - grid.x[axis][1]) / grid.d[axis])), dim))
end

@inline function higher_id(grid::StructuredGrid{dim}, point) where dim
    id = lower_id(grid, point)
    for axis = 1:dim
        add_cartesian(id, axis, 1)
    end
    return id
end

@inline function region_indices!(grid::StructuredGrid{dim}, start, stop; extension::Int = 0) where dim

    l_id, h_id = lower_id(grid, start), higher_id(grid, stop)
    idrange = ntuple(axis -> range(max(grid.nbound+1,l_id[axis] - extension), min(grid.nel[axis]+grid.nbound, h_id[axis] + extension), step = 1), dim)

    return view(grid.indices, idrange...)
end

@inline function in_computational_domain(grid::StructuredGrid, id::CartesianIndex)
    return in_computational_domain(grid, getcoordinates(grid, id))  
end

@inline function in_computational_domain(grid::StructuredGrid, point::Union{Tuple,Vector})
    return betweeneq(point, grid.start, grid.stop)  
end