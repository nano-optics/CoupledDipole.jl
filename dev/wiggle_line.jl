# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFramesMeta
using DataFrames
using Rotations
using LaTeXStrings
using AlgebraOfGraphics, Makie, CairoMakie, ColorSchemes

## this example looks at a diffrative line of Au nanorods in glass,
## but the line is also wiggled with a longer period
## normal incidence

## materials
wavelengths = collect(450:1:950.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole
using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFramesMeta
using DataFrames
using Rotations
using LaTeXStrings
using AlgebraOfGraphics, Makie, CairoMakie, ColorSchemes

## this example looks at a diffrative line of Au nanorods in glass,
## but the line is also wiggled with a longer period
## normal incidence

## materials
wavelengths = vcat(collect(450:10:700.0), collect(701:1:849.0), collect(850:2:950.0))
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelengths, media)

function cluster_wiggle(N, Λ, p, A, a, b, c, α=0.0, β=0.0, γ=0.0,
    material="Au", type="particle")

    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles
    rotations = [inv(QuatRotation(RotZYZ(α, β, γ))) for _ ∈ 1:N] # identical particles

    dz = A * sin.(2π / p .* (-(N - 1)*Λ/2:Λ:(N-1)*Λ/2))
    positions = SVector.(-(N - 1)*Λ/2:Λ:(N-1)*Λ/2, zero(eltype(Λ)), dz)
    positions = SVector.(-(N - 1)*Λ/2:Λ:(N-1)*Λ/2, dz, zero(eltype(Λ)))

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)
end

## array geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl0 = cluster_single(40, 80, 40)
# cluster_chain(N, Λ, a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")
cl1 = cluster_chain(200, 600, 40, 80, 40)
p2 = 8 * 600
p2 = 15 * 700
A = 200
cl2 = cluster_wiggle(200, 600, p2, A, 40, 80, 40)
# cl1 = cluster_array(200, 600, 40, 80, 40)

# visualise_threejs(cl2)

## incidence: along z
Incidence = [RotZYZ(0.0, 0.0, 0.0)]


disp0 = spectrum_dispersion(cl0, mat, Incidence)
disp1 = spectrum_dispersion(cl1, mat, Incidence)
disp2 = spectrum_dispersion(cl2, mat, Incidence)


d0 = dispersion_df(disp0, mat.wavelengths)
d1 = dispersion_df(disp1, mat.wavelengths)
d2 = dispersion_df(disp2, mat.wavelengths)

dat0 = data(@rsubset(d0, :crosstype == "extinction", :polarisation == "pol2"))
dat1 = data(@rsubset(d1, :crosstype == "extinction", :polarisation == "pol2"))
dat2 = data(@rsubset(d2, :crosstype == "extinction", :polarisation == "pol2"))

m = mapping(:wavelength, :value, col=:polarisation)

layer0 = dat0 * m * visual(Lines, linestyle=:dot)
layer1 = dat1 * m * visual(Lines, linestyle=:dash, color=:red)
layer2 = dat2 * m * visual(Lines)
fg = draw(layer0 + layer1 + layer2, facet=(; linkyaxes=:none), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 3, categorical=true)[1:2]))

fg

# save("figure.pdf", fg, px_per_unit=3)

