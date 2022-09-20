const FIT_RADIUS = 0.75

function local_fitting!(f::Fluid{dim}, fitting_point) where dim
    point = collect(fitting_point)
    r = FIT_RADIUS * maximum(f.grid.d)
    start = point .- r
    stop  = point .+ r

    rindices = region_indices!(f.grid, Tuple(start), Tuple(stop))
    sampling_indices = rindices[map(m->m ∈ (1,), f.marker[rindices])]
    points = getcoordinates_in_region(f.grid, sampling_indices)
    
    # ----------------
    # 性能瓶颈
    A = mls_basis(ntuple(i->points[i,:], dim)...)
    pinv_A = pinv(A)
    # ----------------

    samples = f.w[:, sampling_indices]'
    X = pinv_A * samples
    fitting_value = vec(mls_basis(Tuple(point)...) * X)

    return fitting_value
end

@inline function image2ghost(w::Vector{Float64}, n::Vector{Float64})
    wg = copy(w)
    v = wg[2:end-1]
    # 伽利略的无粘边界反射速度，多出的 2*ub 反映参考系变换
    wg[2:end-1] = v - 2 * (v' * n * n )
    return wg
end

# -------------------------------
# mls basis

@inline function mls_basis(x::Float64)
    return [1 x x^2]
end

@inline function mls_basis(x::Float64, y::Float64)
    return [1 x y x^2 y^2 x*y]
end

@inline function mls_basis(x::Float64, y::Float64, z::Float64)
    return [1 x y z x^2 y^2 z^2 x*y y*z z*x]
end

@inline function mls_basis(x::Vector{Float64})
    return [ones(Float64, length(x)) x x.^2]
end

@inline function mls_basis(x::Vector{Float64}, y::Vector{Float64})
    return [ones(Float64, length(x)) x y x.^2 y.^2 x.*y]
end

@inline function mls_basis(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64})
    return [ones(Float64, length(x)) x y z x.^2 y.^2 z.^2 x.*y y.*z z.*x]
end