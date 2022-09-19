function behind_shock(ρ1, u1, p1, Mach, γ, shock_direction)
    @assert Mach >= 1
    c1 = sound_speed(ρ1, p1, γ)
    us = Mach * c1 * shock_direction
    U1 = u1 - us
    M1 = U1 / c1
    r = (γ+1)*M1^2 / (2 + (γ-1)*M1^2)
    U2 = U1 / r
    u2 = U2 + us
    ρ2 = ρ1 * r
    p2 = p1 * (1 + 2*γ/(γ+1) * (M1^2-1))
    return ρ2, u2, p2
end