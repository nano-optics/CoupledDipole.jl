
function alpha_majic(λ, ε, ε_m, a, c)

    ε_r = (ε / ε_m)
    k = 2π / λ * √ε_m

    V = 4 / 3 * π * a * a * c
    e = sqrt(c^2 - a^2) / c
    La = -1 / (2 * e^2) * ((1 - e^2) / e * atanh(e) - 1)
    Lc = (1 - e^2) / e^2 * (atanh(e) / e - 1)

    AlphaS_x = V / 4 / π * (ε_r - 1) / (1 + (ε_r - 1) * La)
    AlphaS_z = V / 4 / π * (ε_r - 1) / (1 + (ε_r - 1) * Lc)

    Omega_x = 1 / 5 * (ε_r - 2 + 3 * e^2) / (1 + (ε_r - 1) * La) - 12 / 25 * e^2
    Omega_z = 1 / 5 * (ε_r - 2 - ε_r * e^2) / (1 + (ε_r - 1) * Lc) + 9 / 25 * e^2

    Ax = AlphaS_x / (1 - Omega_x * k^2 * c^2 - 1im * 2 / 3 * k^3 * AlphaS_x)
    Az = AlphaS_z / (1 - Omega_z * k^2 * c^2 - 1im * 2 / 3 * k^3 * AlphaS_z)

    SVector(Ax, Az)
end
