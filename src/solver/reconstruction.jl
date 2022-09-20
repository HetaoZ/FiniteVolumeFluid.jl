
function Muscl(order::Int = 2, limiter = MinmodLimiter())
    if order == 2
        stencil_width = 3
    else
        error("undef order")
    end
    return Muscl(order, stencil_width, limiter)
end

function Weno(order::Int = 5)
    if order == 5
        stencil_width = 5
    else
        error("undef order")
    end
    return Weno(order, stencil_width, "JS", 
    1.e-10,
    ([1/3 -7/6 11/6],
    [-1/6 5/6 1/3],
    [1/3 5/6 -1/6],
    [11/6 -7/6 1/3]),
    (0.1, 0.6, 0.3),
    2  # recommended p = r
    )
end

function Ausm()
    return Ausm(0.25, 0.75, 1.0)
end

function getnbound(reco::AbstractReconstruction)
    return ceil(Int,reco.stencil_width/2)
end

# --------------------------------------------------------

reconstuct(ws::Array{Float64,2}, axis::Int, material::AbstractMaterial, muscl::Muscl) = reco_by_muscl(ws, axis, material, muscl)
reconstuct(ws::Array{Float64,2}, axis::Int, material::AbstractMaterial, weno::Weno) = reco_by_weno(ws, axis, material, weno)

# --------------------------------------------------------------------------------
# MUSCL scheme

function reco_by_muscl(ws::Array{Float64, 2}, axis::Int, material::AbstractMaterial, muscl::Muscl)

    wL, wR = muscl_interp(ws[:,1], ws[:,2], ws[:,3], ws[:,4], material, muscl.limiter)

    fL = cons2flux(axis, wL, material)
    fR = cons2flux(axis, wR, material)

    return wL, wR, fL, fR
end

function muscl_interp(w1::Vector{Float64}, w2::Vector{Float64}, w3::Vector{Float64}, w4::Vector{Float64}, idealgas::IdealGas, ::MinmodLimiter)
    n = length(w1)
    wL = zeros(Float64, n)
    wR = zeros(Float64, n)

    for i = 1:n
        rL = limited_r(w3[i] - w2[i], w2[i] - w1[i])
        rR = limited_r(w4[i] - w3[i], w3[i] - w2[i])
    
        wL[i] =  w2[i] + 0.5 * (w2[i] - w1[i]) * minmod_limiter(rL)
        wR[i] =  w3[i] - 0.5 * (w3[i] - w2[i]) * minmod_limiter(rR)
    end

    if wL[1] < 1e-14 || pressure(wL, idealgas) < 0
        wL = w2
    end
    if wR[1] < 1e-14 || pressure(wR, idealgas) < 0
        wR = w3
    end
    return wL, wR
end

function limited_r(d1::Float64, d2::Float64)
    if d2 == 0.
        return 0.
    else
        return d1/ d2
    end
end

# ---------------------------------------------------------------------
# WENO-JS scheme

function reco_by_weno(ws::Array{Float64, 2}, axis::Int, idealgas::IdealGas, weno::Weno)
    wL, wR = weno_interp(ws[:,1], ws[:,2], ws[:,3], ws[:,4], ws[:,5], ws[:,6], weno)
    fL = cons2flux(axis, wL, idealgas)
    fR = cons2flux(axis, wR, idealgas)
    return wL, wR, fL, fR
end

"Jiang-Shu's WENO scheme"
function weno_interp(FM3::Vector{Float64},FM2::Vector{Float64},FM1::Vector{Float64},FP1::Vector{Float64},FP2::Vector{Float64},FP3::Vector{Float64}, weno::Weno)
    # WENO sample points: (FM3 - FM2 - FM1 - pipe - FP1 - FP2 - FP3)
    # Jiang-Shu's WENO Scheme
    n = length(FM3)
    FL = [weno_getf(FM3[k], FM2[k], FM1[k], FP1[k], FP2[k], weno.C, weno.a[1:3], weno.eps, weno.p) for k = 1:n]
    FR = [weno_getf(FM2[k], FM1[k], FP1[k], FP2[k], FP3[k], reverse(weno.C), weno.a[2:4], weno.eps, weno.p) for k = 1:n]
    return FL,FR
end

function weno_getbeta(f1::Float64, f2::Float64, f3::Float64, f4::Float64, f5::Float64)
    return (
        13/12*( f1-2*f2+f3 )^2 + 0.25*( f1-4*f2+3*f3 )^2,
        13/12*( f2-2*f3+f4 )^2 + 0.25*( f2       -f4 )^2,
        13/12*( f3-2*f4+f5 )^2 + 0.25*( 3*f3-4*f4+f5 )^2
    )
end

function weno_weightedsum(weight, a1, a2, a3, f1::Float64, f2::Float64, f3::Float64, f4::Float64, f5::Float64)
    return dot(weight, (a1[1]*f1 + a1[2]*f2 + a1[3]*f3, a2[1]*f2 + a2[2]*f3 + a2[3]*f4, a3[1]*f3 + a3[2]*f4 + a3[3]*f5))
end

function weno_getf(f1::Float64, f2::Float64, f3::Float64, f4::Float64, f5::Float64, C::NTuple{3,Float64}, a::NTuple{3,Array{Float64}}, eps::Float64, p::Int)
    β = weno_getbeta(f1,f2,f3,f4,f5)
    α = [C[j]/(eps + β[j])^p for j in eachindex(β)]
    s = sum(α)
    weight = α ./ s
    return weno_weightedsum(weight, a[1], a[2], a[3], f1, f2, f3, f4, f5)
end

# "WENO-Z factor"
# function tau(b1::Float64,b2::Float64,b3::Float64,p::Float64)
#     #tau=(abs(b1-b2)**p +abs(b2-b3)**p +abs(b3-b1)**p )/3.d0
#     #tau=(abs(b1-b2)**p *abs(b2-b3)**p *abs(b3-b1)**p )**(1.d0/3.d0)
#     tau=3/(1/abs(b1-b2)^p +1/abs(b2-b3)^p +1/abs(b3-b1)^p )
#     #tau = abs(b3-b1)
#     return tau
# end

# -----------------------------------------------
# Limiters in high-to-low order of resolution of shock saves
# -----------------------------------------------
limiter(r::Float64, ::SuperbeeLimiter) = superbee_limiter(r)
limiter(r::Float64, ::VanLeerMeanLimiter) = van_leer_mean_limiter(r)
limiter(r::Float64, ::VanLeerLimiter) = van_leer_limiter(r)
limiter(r::Float64, ::VanAlbabaLimiter) = van_albaba_limiter(r)
limiter(r::Float64, ::MinmodLimiter) = minmod_limiter(r)

"""
superbee limiter
"""
@inline function superbee_limiter(r::Float64)
    return max(min(2*r,1), min(r,2))
end

"""
van Leer mean (double minmod) limiter
"""
@inline function van_leer_mean_limiter(r::Float64)
    if r <= 0.
        return 0.
    else
        return minimum([2*r, 2, (1+r)/2])
    end
end

"""
van Leer limiter
"""
@inline function van_leer_limiter(r::Float64)
    if 1 + r == 0.
        return 0.
    else
        return (r + abs(r)) / (1 + r)
    end
end

"""
van Albaba limiter
"""
@inline function van_albaba_limiter(r::Float64)
    return (r^2+r)/(1+r^2)
end

"""
minmod limiter
"""
@inline function minmod_limiter(r::Float64)
    if r <= 0.
        return 0.
    else
        return min(r, 1.)
    end
end

"""
minmod function
"""
@inline function minmod(a::Float64, b::Float64)
    if abs(a) < abs(b)
        c = a
    elseif abs(a) > abs(b)
        c = b
    else
        c = 0.
    end
    return c
end