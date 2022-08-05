reconstuct!(ws::Array{Float64,2}, axis::Int, Gamma::Float64, reco_scheme::Muscl) = reco_by_muscl!(ws, axis, Gamma)
reconstuct!(ws::Array{Float64,2}, axis::Int, Gamma::Float64, reco_scheme::Weno) = reco_by_weno!(ws, axis, Gamma)

# --------------------------------------------------------------------------------
# MUSCL scheme

function reco_by_muscl!(ws::Array{Float64, 2}, axis::Int, Gamma::Float64)

    wL, wR = muscl_interp(ws[:,1], ws[:,2], ws[:,3], ws[:,4])

    fL = cons2flux(axis, Gamma, wL)
    fR = cons2flux(axis, Gamma, wR)

    return wL, wR, fL, fR
end

function muscl_interp(w1::Vector{Float64}, w2::Vector{Float64}, w3::Vector{Float64}, w4::Vector{Float64})
    rL = @. map(limited_r, w3 - w2, w2 - w1)
    rR = @. map(limited_r, w4 - w3, w3 - w2)

    wL = @. w2 + 0.5 * (w2 - w1) * limiter(rL, MinmodLimiter())
    wR = @. w3 - 0.5 * (w3 - w2) * limiter(rR, MinmodLimiter())

    if wL[1] < 1e-14 || pressure(wL, Gamma) < 0
        wL = w2
    end
    if wR[1] < 1e-14 || pressure(wR, Gamma) < 0
        wR = w3
    end
    return wL, wR
end

@inline function limited_r(d1::Float64, d2::Float64)
    if d2 == 0.
        return 0.
    else
        return d1/ d2
    end
end

# ---------------------------------------------------------------------
# WENO-JS scheme

function reco_by_weno!(ws::Array{Float64, 2}, axis::Int, Gamma::Float64)
    wL, wR = weno_interp(ws[:,1], ws[:,2], ws[:,3], ws[:,4], ws[:,5], ws[:,6])

    fL = cons2flux(axis, Gamma, wL)
    fR = cons2flux(axis, Gamma, wR)

    return wL, wR, fL, fR
end

"Jiang-Shu's WENO scheme"
function weno_interp(FM3::Vector{Float64},FM2::Vector{Float64},FM1::Vector{Float64},FP1::Vector{Float64},FP2::Vector{Float64},FP3::Vector{Float64})
    # WENO sample points: (FM3 - FM2 - FM1 - pipe - FP1 - FP2 - FP3)
    # Jiang-Shu's WENO Scheme
    eps=1.e-10
    a0 = [1/3 -7/6 11/6]
    a1 = [-1/6 5/6 1/3]
    a2 = [1/3 5/6 -1/6]
    a3 = [11/6 -7/6 1/3]
    C = [0.1 0.6 0.3]
    p = 2  # recommended p = r
    FL = Vector{Float64}(undef,4)
    FR = Vector{Float64}(undef,4)
    for k=1:4
      ## for WL
      beta0=13/12*( FM3[k]-2*FM2[k]+FM1[k] )^2 + 0.25*( FM3[k]-4*FM2[k]+3*FM1[k] )^2
      beta1=13/12*( FM2[k]-2*FM1[k]+FP1[k] )^2 + 0.25*( FM2[k]           -FP1[k] )^2
      beta2=13/12*( FM1[k]-2*FP1[k]+FP2[k] )^2 + 0.25*( 3*FM1[k]-4*FP1[k]+FP2[k] )^2
      alpha0=C[1]/(eps+beta0)^p
      alpha1=C[2]/(eps+beta1)^p
      alpha2=C[3]/(eps+beta2)^p
      weight0=alpha0/(alpha0+alpha1+alpha2)
      weight1=alpha1/(alpha0+alpha1+alpha2)
      weight2=alpha2/(alpha0+alpha1+alpha2)
      FL[k] = weight0*( a0[1]*FM3[k] + a0[2]*FM2[k] + a0[3]*FM1[k]) + weight1*( a1[1]*FM2[k] + a1[2]*FM1[k] + a1[3]*FP1[k]) + weight2*( a2[1]*FM1[k] + a2[2]*FP1[k] + a2[3]*FP2[k])


      # for WR
      beta0=13/12*( FM2[k]-2*FM1[k]+FP1[k] )^2 + 0.25*( FM2[k]-4*FM1[k]+3*FP1[k] )^2
      beta1=13/12*( FM1[k]-2*FP1[k]+FP2[k] )^2 + 0.25*( FM1[k]-FP2[k] )^2
      beta2=13/12*( FP1[k]-2*FP2[k]+FP3[k] )^2 + 0.25*( 3*FP1[k]-4*FP2[k]+FP3[k] )^2
      alpha0=C[3]/(eps+beta0)^p
      alpha1=C[2]/(eps+beta1)^p
      alpha2=C[1]/(eps+beta2)^p
      weight0=alpha0/(alpha0+alpha1+alpha2)
      weight1=alpha1/(alpha0+alpha1+alpha2)
      weight2=alpha2/(alpha0+alpha1+alpha2)
      FR[k] =  weight0*( a1[1]*FM2[k] + a1[2]*FM1[k] + a1[3]*FP1[k]) + weight1*( a2[1]*FM1[k] + a2[2]*FP1[k] + a2[3]*FP2[k]) + weight2*( a3[1]*FP1[k] + a3[2]*FP2[k] + a3[3]*FP3[k])

    end

    return FL,FR
end

"WENO-Z factor"
function tau(b1::Float64,b2::Float64,b3::Float64,p::Float64)
    #tau=(abs(b1-b2)**p +abs(b2-b3)**p +abs(b3-b1)**p )/3.d0
    #tau=(abs(b1-b2)**p *abs(b2-b3)**p *abs(b3-b1)**p )**(1.d0/3.d0)
    tau=3/(1/abs(b1-b2)^p +1/abs(b2-b3)^p +1/abs(b3-b1)^p )
    #tau = abs(b3-b1)
    return tau
end

# -----------------------------------------------
# Limiters in high-to-low order of resolution of shock saves
# -----------------------------------------------
abstract type AbstractLimiter end

struct SuperbeeLimiter <: AbstractLimiter end
struct VanLeerMeanLimiter <: AbstractLimiter end
struct VanLeerLimiter <: AbstractLimiter end
struct VanAlbabaLimiter <: AbstractLimiter end
struct MinmodLimiter <: AbstractLimiter end

limiter(r::Float64, l::SuperbeeLimiter) = superbee_limiter(r)
limiter(r::Float64, l::VanLeerMeanLimiter) = van_leer_mean_limiter(r)
limiter(r::Float64, l::VanLeerLimiter) = van_leer_limiter(r)
limiter(r::Float64, l::VanAlbabaLimiter) = van_albaba_limiter(r)
limiter(r::Float64, l::MinmodLimiter) = minmod_limiter(r)

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