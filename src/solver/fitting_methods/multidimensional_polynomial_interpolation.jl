const MP_RADIUS = 1.5

"""
multidimensional polynomial interpolation
"""
function local_fitting!(f::Fluid{dim}, fitting_point) where dim
    point = collect(fitting_point)
    r = MP_RADIUS * maximum(f.grid.d)
    start = point .- r
    stop  = point .+ r

    rindices = region_indices!(f.grid, Tuple(start), Tuple(stop))
    sampling_indices = rindices[map(rho -> rho > 1.e-14, f.rho[rindices])]
    points = getcoordinates_in_region(f.grid, sampling_indices)

    weights = scatter_polynomial_interpolant(point, points, f.grid.d * ceil(MP_RADIUS))
    samples = f.w[:, sampling_indices]
    fitting_value = samples * weights

    return fitting_value
end


function scatter_polynomial_interpolant(x, points, d)
    weights = [polynomial_interpolant(x, points[:,j], d) for j in axes(points,2)]
    return weights / sum(weights)
end

"Polynomial basis function (PBF) 多项式基函数，区别于径向基函数 (radial basis function)。其优点是在任一象限内可以再生多项式分布。"
function polynomial_interpolant(x, p, d)
    ξ = (x - p) ./ d
    return prod([abs(1 - ξi) for ξi in ξ])
end