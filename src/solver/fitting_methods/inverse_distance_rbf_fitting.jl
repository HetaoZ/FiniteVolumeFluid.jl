const INVERSE_DISTANCE_FIT_RADIUS = 1.5

"""
inverse distance RBF fitting
"""
function local_fitting!(f::Fluid{dim}, fitting_point) where dim
    point = collect(fitting_point)
    r = INVERSE_DISTANCE_FIT_RADIUS * maximum(f.grid.d)
    start = point .- r
    stop  = point .+ r

    rindices = region_indices!(f.grid, Tuple(start), Tuple(stop))
    sampling_indices = rindices[map(rho -> rho > 1.e-14, f.rho[rindices])]
    points = getcoordinates_in_region(f.grid, sampling_indices)
    
    weights = scatter_linear_interpolant(point, points)
    samples = f.w[:, sampling_indices]
    fitting_value = samples * weights

    return fitting_value
end

"""
在一组散点上，根据反距离函数计算插值系数。
"""
function scatter_linear_interpolant(point::Vector{Float64}, x)
    inverse_distances = [1/norm(point - x[:,j]) for j in axes(x,2)]
    return inverse_distances / sum(inverse_distances)
end