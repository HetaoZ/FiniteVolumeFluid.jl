

include("fitting_methods/multidimensional_polynomial_interpolation.jl")
# include("fitting_methods/inverse_distance_rbf_fitting.jl")
# include("fitting_methods/moving_least_square_fitting.jl")

@inline function image2ghost(w::Vector{Float64}, n::Vector{Float64})
    wg = copy(w)
    v = wg[2:end-1]
    # 伽利略的无粘边界反射速度，多出的 2*ub 反映参考系变换
    wg[2:end-1] = v - 2 * (v' * n * n )
    return wg
end