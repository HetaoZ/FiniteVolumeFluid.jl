@inline function prim2cons(rho::Float64, u::Vector{Float64}, e::Float64)
    w = ones(Float64, length(u)+2)
    w[2:end-1] = u
    w[end] = e + 0.5*norm(u)^2
    return rho .* w
end

@inline function pressure(rho::Float64, e::Float64, Gamma::Float64)
    return Gamma * rho * e
end

@inline function pressure(w::Vector{Float64}, Gamma::Float64)
    if w[1] > 1e-14
        return Gamma * w[1] * (w[end]/w[1] - 0.5 * norm(w[2:end-1] / w[1])^2)
    elseif w[1] >= 0.
        return 0.
    else
        println("w = ", w)
        error("wrong pressure")
    end
end

@inline function pressure2e(rho::Float64, p::Float64, inv_gamma_minus_one::Float64)
    if rho < 1e-14
        return 0.
    elseif rho > 0 && p >= 0
        return p / rho * inv_gamma_minus_one
    else
        println("(rho, p) = ", (rho, p))
        error("error when pressure2e")
    end
end

@inline function temperature(rho::Float64, p::Float64, gamma::Float64)
    if rho < 1e-14
        return 0.
    elseif rho > 0 
        if p >= 0
            return  gamma * p / rho
        else
            return 0.
        end
    else

        println("(rho, p) = ", (rho, p))
        error("error when temperature")
    end
end

@inline function sound_speed(rho::Float64, p::Float64, gamma::Float64)
    return sqrt(temperature(rho, p, gamma))
end

@inline function sound_speed(w::Vector{Float64}, gamma::Float64)
    return sqrt(temperature(w[1], pressure(w, gamma - 1.0), gamma))
end

@inline function speed(w::Vector{Float64}, axis::Int)
    dim = length(w) - 2
    w[1] < 1e-14 ? zeros(Float64, dim) : w[2:end-1]/w[1]
end

@inline function reynolds_number(rho::Float64, para::Dict)
    return rho * para["U0"] * para["L0"] / para["mu"]
end

@inline function cons2prim(w::Vector{Float64}, Gamma::Float64)
    dim = length(w) - 2
    if w[1] < 1e-14
        rho = 0.
        u = zeros(Float64, dim)
        e = 0.
        p = 0. 
    elseif w[1] > 0 
        rho = w[1]
        u = w[2:end-1] / w[1]
        e = w[end] / w[1] - 0.5 * norm(u)^2
        p = pressure(rho, e, Gamma) 
    else
        println("w = ",w)
        error("wrong w")
    end
    return rho, u, e, p
end

@inline function cons2flux(axis::Int, Gamma::Float64, w::Vector{Float64})
    
    if w[1] > 0.
        rho, u, E = w[1], w[2:end-1]/w[1], w[end]/w[1]
        e = E - 0.5 * norm(u)^2
        p = pressure(rho, e, Gamma)

        flux = zeros(Float64, length(w))
        flux[1] = rho * u[axis]
        flux[2:end-1] = rho * u[axis] .* u
        flux[1+axis] += p
        flux[end] = (rho * E + p) * u[axis]
    elseif w[1] == 0.
        flux = zeros(Float64, length(w))
    else
        println("axis, w = ", (axis, w))
        error("")
    end
    return flux
end

@inline function cons2flux(axis::Int, Gamma::Float64, ws::Vector{Float64}...)
    return map(w->cons2flux(axis, Gamma, w), ws)
end

# @inline function flux2prim(axis::Int, Gamma::Float64, flux::Vector{Float64})
#     rho = flux
# end