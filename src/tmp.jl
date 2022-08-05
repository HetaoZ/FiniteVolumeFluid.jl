FL[k] = weight0*( a1[1]*f1 + a1[2]*f2 + a1[3]*f3) + weight1*( a2[1]*f2 + a2[2]*f3 + a2[3]*f4) + weight2*( a3[1]*f3 + a3[2]*f4 + a3[3]*f5)

FR[k] = weight0*( a1[1]*f1 + a1[2]*f2 + a1[3]*f3) + weight1*( a2[1]*f2 + a2[2]*f3 + a2[3]*f4) + weight2*( a3[1]*f3 + a3[2]*f4 + a3[3]*f5)

function weno_getf(weight::NTuple{3,Float64}, a1::Vector{Float64}, a2::Vector{Float64}, a3::Vector{Float64}, f1::Float64, f2::Float64, f3::Float64, f4::Float64, f5::Float64)
    return dot(weight, (a1[1]*f1 + a1[2]*f2 + a1[3]*f3, a2[1]*f2 + a2[2]*f3 + a2[3]*f4, a3[1]*f3 + a3[2]*f4 + a3[3]*f5)))
end