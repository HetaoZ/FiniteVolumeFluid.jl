reconstuct!(ws::Array{Float64,2}, axis::Int, Gamma::Float64, reco_scheme::Muscl) = reco_by_muscl!(ws, axis, Gamma)

function reco_by_muscl!(ws::Array{Float64, 2}, axis::Int, Gamma::Float64)

    rL = @. map(limited_r, ws[:,3] - ws[:,2], ws[:,2] - ws[:,1])
    rR = @. map(limited_r, ws[:,4] - ws[:,3], ws[:,3] - ws[:,2])

    wL = @. ws[:,2] + 0.5 * (ws[:,2] - ws[:,1]) * limiter(rL, MinmodLimiter())
    wR = @. ws[:,3] - 0.5 * (ws[:,3] - ws[:,2]) * limiter(rR, MinmodLimiter())

    if wL[1] < 1e-14 || pressure(wL, Gamma) < 0
        wL = ws[:,2]
    end
    if wR[1] < 1e-14 || pressure(wR, Gamma) < 0
        wR = ws[:,3]
    end

    fL = cons2flux(axis, Gamma, wL)
    fR = cons2flux(axis, Gamma, wR)

    return wL, wR, fL, fR
end

@inline function limited_r(d1, d2)
    if d2 == 0.
        return 0.
    else
        return d1/ d2
    end
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