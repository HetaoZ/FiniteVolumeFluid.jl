function Ausm()
    return Ausm(AUSM_Kp, AUSM_Ku, AUSM_sigma)
end

# --------------------------------------------------
# flux functions

function compute_flux(ws::Array{Float64, 2}, axis::Int, solver::FVSolver, material::AbstractMaterial)
    return numerical_flux(reconstuct(ws, axis, material, solver.reconstruction, solver.flux)..., axis, material)
end

"Lax-Friedrichs flux"
function numerical_flux(wL::Array{Float64,1}, wR::Array{Float64,1}, fL::Array{Float64,1}, fR::Array{Float64,1}, ::Int, material::AbstractMaterial)
    return get_lf_flux(wL, wR, fL, fR, material)
end

"AUSM^+-up flux"
function numerical_flux(wL::Array{Float64,1}, wR::Array{Float64,1}, axis::Int, material::AbstractMaterial)
    return get_ausm_flux(axis, wL, wR, material)
end

# _______________________________________________________

"Lax-Friedrichs flux scheme"
function get_lf_flux(wL::Vector{Float64}, wR::Vector{Float64}, fL::Vector{Float64}, fR::Vector{Float64}, idealgas::IdealGas)

    return 0.5 * ( fL + fR - max(norm(speed(wL)) + sound_speed(wL, idealgas), norm(speed(wR)) + sound_speed(wR, idealgas)) * ( wR - wL ) )
end


# ------------------------------
# AUSM^+-up

const AUSM_Kp = 0.25
const AUSM_Ku = 0.75
const AUSM_sigma = 1.0

"AUMS^+-up flux scheme"
function get_ausm_flux(axis::Int, wL::Vector{Float64}, wR::Vector{Float64}, idealgas::IdealGas)

    len = length(wL)
    
    rhoL, uL, EL, pL, aL = get_flux_vars(wL, idealgas)
    rhoR, uR, ER, pR, aR = get_flux_vars(wR, idealgas)

    vel_L, vel_R = uL[axis], uR[axis]    

    a_int = max(0.5 * (aL + aR), 1.0)

    if a_int < 1e-14 || isnan(a_int)
        a_int = 1e14
    end

    Ma_L, Ma_R = vel_L / a_int, vel_R / a_int
    
    Ma_bar_sq = (vel_L^2 + vel_R^2) / (2.0*a_int^2)

    # fa = 1.0 # fa = constant

    MMa_p4_L = get_MMa_4_L(Ma_L)
    MMa_m4_R = get_MMa_4_R(Ma_R)
    rho_int = (rhoL + rhoR) * 0.5

    if rho_int < 1e-14
        Ma_int = 0.0
    else
        Ma_int = MMa_p4_L + MMa_m4_R - AUSM_Kp*max(1.0 - AUSM_sigma*Ma_bar_sq, 0.0)*(pR - pL)/(rho_int*a_int^2)
    end

    PP_p5_L = get_PP_5_L(Ma_L)
    PP_m5_R = get_PP_5_R(Ma_R)

    p_int = PP_p5_L*pL + PP_m5_R*pR - AUSM_Ku*PP_p5_L*PP_m5_R*(rhoL + rhoR)*a_int*(vel_R - vel_L)

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

function get_flux_vars(w::Vector{Float64}, idealgas::IdealGas)
    if w[1] < 1e-14
        rho, u, E, p, a = 0., zeros(Float64, length(w)-2), 0., 0., 0.
    else
        rho = w[1]
        u = w[2:end-1] ./ rho
        E = w[end] / rho
        e = E - 0.5 * norm(u)^2
        p = pressure(rho, e, idealgas)
        if p < 0
            a = sound_speed(rho, 1e-14, idealgas)
        else
            a = sound_speed(rho, p, idealgas)
        end
    end
    return rho, u, E, p, a
end

function get_MMa_4_L(Ma::Float64)
    if abs(Ma) >= 1
        MMa_1 = 0.5 * (Ma+abs(Ma))
        MMa_4 = MMa_1
    else
        MMa_2_p =   0.25*(Ma + 1.0)^2
        MMa_2_m = - 0.25*(Ma - 1.0)^2
        MMa_4 = MMa_2_p*(1.0 - 2.0*MMa_2_m)
    end
    return MMa_4
end

function get_MMa_4_R(Ma::Float64)
    if abs(Ma) >= 1
        MMa_1 = 0.5 * (Ma-abs(Ma))
        MMa_4 = MMa_1
    else
        MMa_2_p =   -0.25*(Ma - 1.0)^2
        MMa_2_m =    0.25*(Ma + 1.0)^2
        MMa_4 = MMa_2_p*(1.0 + 2.0*MMa_2_m)
    end
    return MMa_4
end

function get_PP_5_L(Ma::Float64)
    if abs(Ma) >= 1 
        MMa_1 = 0.5*(Ma+abs(Ma))
        PP_5 = MMa_1 / Ma
    else
        MMa_2_p =   0.25*(Ma+1.0)^2
        MMa_2_m = - 0.25*(Ma-1.0)^2
        PP_5 = MMa_2_p * (2.0 - Ma - 3.0*Ma*MMa_2_m)
    end
    return PP_5
end

function get_PP_5_R(Ma::Float64)
    if abs(Ma) >= 1 
        MMa_1 = 0.5*(Ma-abs(Ma))
        PP_5 = MMa_1 / Ma
    else
        MMa_2_p =   -0.25*(Ma-1.0)^2
        MMa_2_m =    0.25*(Ma+1.0)^2
        PP_5 = MMa_2_p * (-2.0 - Ma + 3.0*Ma*MMa_2_m)
    end
    return PP_5
end