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
# wavelengths = vcat(collect(450:10:700.0), collect(701:1:849.0), collect(850:2:950.0))
wavelengths = collect(400:2:800.0)
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.5)])
mat = Material(wavelengths, media)


function cluster_spokes(Narm, Nspokes, Λ, a, b, c,
    material="Au", type="particle")

    N = Narm * Nspokes
    sizes = [SVector(a, b, c) for ii in 1:N] # identical particles

    angles = [2π * n / Nspokes for n = 1:Nspokes]
    Z = [exp(im * 2π * n / Nspokes) for n = 1:Nspokes]
    R = 1:Narm
    positions = [r * Λ * SVector(real(z), imag(z), 0.0) for z in Z for r in R]

    rotations = [inv(QuatRotation(RotZYZ(a + π / 4, 0, 0))) for a in angles for _ in R]

    Cluster(positions, rotations, sizes, [material for _ ∈ 1:N], type)
end

# cl = cluster_spokes(4, 6, 100, 10, 20, 10)
# visualise_threejs(cl)


## array geometry
# N, Λ, a, b, c, φ, θ, ψ, material = "Au", type="particle"
cl0 = cluster_single(40, 50, 40)
# cluster_chain(N, Λ, a, b, c, α=0.0, β=0.0, γ=0.0, material="Au", type="particle")
cl1 = cluster_chain(100, 450, 40, 50, 40)
cl2 = cluster_spokes(50, 6, 450, 40, 50, 40)
# cl1 = cluster_array(200, 600, 40, 80, 40)

# visualise_threejs(cl2)

## incidence: along z
Incidence = [RotZYZ(0.0, 0.0, 0.0)]


disp0 = spectrum_dispersion(cl0, mat, Incidence, polarisation="circular")
disp1 = spectrum_dispersion(cl1, mat, Incidence, polarisation="circular")
disp2 = spectrum_dispersion(cl2, mat, Incidence, polarisation="circular")


d0 = dispersion_df(disp0, mat.wavelengths)
d1 = dispersion_df(disp1, mat.wavelengths)
d2 = dispersion_df(disp2, mat.wavelengths)

dat0 = data(@rsubset(d0, :crosstype == "extinction", :polarisation == "pol1"))
dat1 = data(@rsubset(d1, :crosstype == "extinction", :polarisation == "pol1"))
dat2 = data(@rsubset(d2, :crosstype == "extinction", :polarisation == "pol1"))

m = mapping(:wavelength, :value, col=:polarisation)

layer0 = dat0 * m * visual(Lines, linestyle=:dot)
layer1 = dat1 * m * visual(Lines, linestyle=:dash, color=:red)
layer2 = dat2 * m * visual(Lines)
fg = draw(layer0 + layer1 + layer2, facet=(; linkyaxes=:none), axis=(; xlabel="wavelength /nm", ylabel="cross-section σ /nm²"),
    palettes=(; color=cgrad(ColorSchemes.phase.colors, 3, categorical=true)[1:2]))

fg

gd = groupby(d2, [:crosstype, :polarisation]);
@combine(gd, :peak = :wavelength[findmax(:value)[2]])



## near field map

## materials
wavelengths = [680.0]
media = Dict([("Au", epsilon_Au), ("medium", x -> 1.5)])
mat = Material(wavelengths, media)

Incidence = [RotZ(0.0)]
res = 150
x = range(-1500.0, 1500, res)
y = range(-1500.0, 1500, res)

# x = range(-10000.0, 10000, res)
# y = range(-10000.0, 10000, res)

probes = SVector.(Iterators.product(x, y, zero(eltype(x))))[:]

E², B², 𝒞, positions = map_nf(probes, cl2, mat, Incidence, polarisation="circular"; evaluate_inside=false)

df2 = (; x=positions.x, y=positions.y, z=log10.(E²[:, 2]))
df2 = (; x=positions.x, y=positions.y, z=(𝒞[:, 1]))

layer = visual(Heatmap)
m = mapping(:x, :y, :z)
f2 = data(df2) * layer * m
draw(f2, axis=(; xlabel="x /nm", ylabel="y /nm", aspect=DataAspect()); figure=(resolution=(600, 400),))
