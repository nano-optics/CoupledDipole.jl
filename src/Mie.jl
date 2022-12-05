

"""
ricatti_bessel(x, nmax)

RB functions for convenience
- `x`: size parameter
- `nmax`: maximum order

# Examples

```
besselj(3, 1.0)
besselj.(3, [1.0 2.0])

ricatti_bessel(collect(1.0:0.1:2.0), 3)
```

"""
function ricatti_bessel(x, nmax)

    nu = 0.5 .+ collect(0:nmax)
    bj = hcat([besselj.(_nu, x) .* sqrt.(pi / 2 * x) for _nu in nu]...)
    bh = hcat([besselh.(_nu, x) .* sqrt.(pi / 2 * x) for _nu in nu]...)
    dpsi = [bj[i, j] - j * bj[i, j+1] / x[i]
            for i in 1:length(x), j in 1:nmax]
    dxi = [bh[i, j] - j * bh[i, j+1] / x[i]
           for i in 1:length(x), j in 1:nmax]
    psi = bj[:, 2:nmax+1]
    xi = bh[:, 2:nmax+1]
    return psi, xi, dpsi, dxi
end


"""
mie_susceptibility(x, s, nmax)

Mie susceptibility
- `x`: size parameter
- `s`: relative refractive index
- `nmax`: maximum order

"""
function mie_susceptibility(x, s, nmax)
    z = s .* x
    rbx = ricatti_bessel(x, nmax)
    rbz = ricatti_bessel(z, nmax)

    # aux functions
    pp1 = rbz[1] .* rbx[3]
    pp2 = rbx[1] .* rbz[3]
    pp3 = rbz[1] .* rbx[4]
    pp4 = rbx[2] .* rbz[3]

    # numerators
    gammanum = -pp1 + s .* pp2
    deltanum = pp2 - s .* pp1

    # denominators
    gammadenom = pp3 - s .* pp4
    deltadenom = -pp4 + s .* pp3

    Gamma = gammanum ./ gammadenom
    Delta = deltanum ./ deltadenom
    A = 1im * s ./ gammadenom
    B = 1im * s ./ deltadenom
    return Gamma, Delta, A, B

end

"""
mie_ff(radius, wavelength, epsilon, medium, nmax)

Far-field Mie scattering
- `x`: size parameter
- `s`: relative refractive index
- `nmax`: maximum order
- `medium`: maximum order
- `nmax`: maximum order

"""
function mie_ff(radius, wavelength, epsilon, medium, nmax)

    k0 = 2 * pi ./ wavelength
    k = k0 .* medium
    x = k .* radius
    s = sqrt.(epsilon) ./ medium
    Gamma, Delta, A, B = mie_susceptibility(x, s, nmax)
    Gamma2 = abs.(Gamma) .^ 2
    Delta2 = abs.(Delta) .^ 2
    scatm = (Delta2 .+ Gamma2)

    pref = 2 ./ (x .^ 2)
    cc = (2 * collect(1:nmax) .+ 1)

    qsca = pref .* (scatm * cc)
    qext = -pref .* ((real(Gamma) .+ real(Delta)) * cc)
    qabs = qext - qsca
    return hcat(qext, qsca, qabs)
end
