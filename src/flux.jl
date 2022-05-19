@inline function dflux!(ws::Array{Float64, 2}, axis::Int, reco_scheme::AbstractRecoScheme, flux_scheme::AbstractFluxScheme, parameters::Dict)
    return flux!(ws[:,1:end-1], axis, reco_scheme, flux_scheme, parameters) - flux!(ws[:,2:end], axis, reco_scheme, flux_scheme, parameters)
end

@inline function flux!(ws::Array{Float64, 2}, axis::Int, reco_scheme::AbstractRecoScheme, flux_scheme::AbstractFluxScheme, parameters::Dict)
    wL, wR, fL, fR = reconstuct!(ws, axis, parameters["Gamma"], reco_scheme)
    return numerical_flux!(wL, wR, fL, fR, axis, parameters["gamma"], flux_scheme)
end

@inline function numerical_flux!(wL::Array{Float64,1}, wR::Array{Float64,1}, fL::Array{Float64,1}, fR::Array{Float64,1}, axis::Int, gamma::Float64, flux_scheme::LaxFriedrichs)
    return get_lf_flux!(wL, wR, fL, fR, gamma)
end

@inline function numerical_flux!(wL::Array{Float64,1}, wR::Array{Float64,1}, fL::Array{Float64,1}, fR::Array{Float64,1}, axis::Int, gamma::Float64, flux_scheme::Ausm)
    return get_ausm_flux!(wL, wR, fL, fR, gamma, axis, flux_scheme)
end

## _______________________________________________________

@inline function get_lf_flux!(wL::Vector{Float64}, wR::Vector{Float64}, fL::Vector{Float64}, fR::Vector{Float64}, gamma::Float64)

    return @. 0.5 * ( fL + fR - max(norm(speed(wL, axis)) + sound_speed(wL, gamma), norm(speed(wR, axis)) + sound_speed(wR, gamma)) * ( wR - wL ) )
end

function get_ausm_flux!(wL::Vector{Float64}, wR::Vector{Float64}, fL::Vector{Float64}, fR::Vector{Float64}, gamma::Float64, axis::Int, ausm::Ausm)
    len = length(fL)
    
    rhoL, uL, EL, pL, aL = get_flux_vars(wL, gamma)
    rhoR, uR, ER, pR, aR = get_flux_vars(wR, gamma)

    vel_L, vel_R = uL[axis], uR[axis]    

    a_int = max(0.5 * (aL + aR), 1.0)

    if a_int < 1e-14 || isnan(a_int)
        a_int = 1e14
    end

    Ma_L, Ma_R = vel_L / a_int, vel_R / a_int
    
    Ma_bar_sq = (vel_L^2 + vel_R^2) / (2.0*a_int^2)

    fa = 1.0

    MMa_p4_L = get_MMa_4(Ma_L, fa,  1.0)
    MMa_m4_R = get_MMa_4(Ma_R, fa, -1.0)

    rho_int = (rhoL + rhoR) * 0.5

    if rho_int < 1e-14 || fa < 1e-14
        Ma_int = 0.0
    else
        Ma_int = MMa_p4_L + MMa_m4_R - ausm.Kp/fa*max(1.0 - ausm.sigma*Ma_bar_sq, 0.0)*(pR - pL)/(rho_int*a_int^2)
    end

    PP_p5_L = get_PP_5(Ma_L, fa,  1.0)
    PP_m5_R = get_PP_5(Ma_R, fa, -1.0)

    p_int = PP_p5_L*pL + PP_m5_R*pR - ausm.Ku*PP_p5_L*PP_m5_R*(rhoL + rhoR)*fa*a_int*(vel_R - vel_L)

    vel_int = 0.5*(vel_L + vel_R)

    dm_int = a_int * Ma_int * (Ma_int > 0 ? rhoL : rhoR)
    
    Psi_int = ones(Float64, len)
    if dm_int > 0.0
        Psi_int[2:end-1] = uL
        Psi_int[end] = EL
    else
        Psi_int[2:end-1] = uR
        Psi_int[end] = ER
    end

    P_int = zeros(Float64, len)
    P_int[end] = p_int*vel_int
    P_int[1 + axis] = p_int

    f = dm_int * Psi_int + P_int

    return f
end

function get_flux_vars(w::Vector{Float64}, gamma)
    if w[1] < 1e-14
        rho, u, E, p, a = 0., zeros(Float64, length(w)-2), 0., 0., 0.
    else
        rho = w[1]
        u = w[2:end-1] ./ rho
        E = w[end] / rho
        e = E - 0.5 * norm(u)^2
        p = pressure(rho, e, gamma - 1.0)
        if p < 0
            a = sound_speed(rho, 1e-14, gamma)
        else
            a = sound_speed(rho, p, gamma)
        end
    end
    return rho, u, E, p, a
end

function get_MMa_4(Ma::Float64, fa::Float64, s::Float64)
    if abs(Ma) >= 1
        MMa_1 = 0.5 * (Ma+s*abs(Ma))
        MMa_4 = MMa_1
    else
        MMa_2_p =   s*0.25*(Ma + s*1.0)^2
        MMa_2_m = - s*0.25*(Ma - s*1.0)^2
        MMa_4 = MMa_2_p*(1.0 - s*2.0*MMa_2_m)
    end
    return MMa_4
end

function get_PP_5(Ma::Float64, fa::Float64, s::Float64)
    if abs(Ma) >= 1 
        MMa_1 = 0.5*(Ma+s*abs(Ma))
        PP_5 = MMa_1 / Ma
    else
        MMa_2_p =   s*0.25*(Ma+s*1.0)^2
        MMa_2_m = - s*0.25*(Ma-s*1.0)^2
        alpha = 3.0/16.0 * (-4.0+5.0*fa^2)
        PP_5 = MMa_2_p * ((s*2.0 - Ma) - s*16.0*alpha*Ma*MMa_2_m)
    end
    return PP_5
end

