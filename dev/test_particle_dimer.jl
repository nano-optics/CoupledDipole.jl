# include("../src/CoupledDipole.jl")
push!(LOAD_PATH, expanduser( "~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

using LinearAlgebra
using StaticArrays
using FastGaussQuadrature
using DataFrames
using VegaLite

## materials
wavelength = collect(450:2:850.0)
medium = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
mat = Material(wavelength, medium)

## dimer geometry
cl1 = cluster_dimer(80, 20, 20, 40, π/4)
cl2 = cluster_single(20, 20, 40)

## incidence: along z, along x, along y
Incidence = [SVector(0,0,0),SVector(0,π/2,0),SVector(π/2,π/2,0)]

disp1 = spectrum_dispersion(cl1, mat, Incidence)
disp2 = spectrum_dispersion(cl2, mat, Incidence)

function dispersion_df(x, wavelength)

    e = insertcols!(
        DataFrame(x.extinction, :auto),
        :wavelength => wavelength,
        :crosstype => "extinction"
    )
    a = insertcols!(
        DataFrame(x.absorption, :auto),
        :wavelength => wavelength,
        :crosstype => "absorption"
    )

    s = insertcols!(
        DataFrame(x.scattering, :auto),
        :wavelength => wavelength,
        :crosstype => "scattering"
    )

    stack([e;a;s], Not([:wavelength,:crosstype]))

end

d1 = dispersion_df(disp1, mat.wavelength)
d2 = dispersion_df(disp2, mat.wavelength)

d = [insertcols!(d1, :cluster => "dimer");
     insertcols!(d2, :cluster => "single")]

 @vlplot(data=d,
 width= 400,
 height =  300,
     mark = {:line},
     row = "crosstype",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )



oa1 = spectrum_oa(cl1, mat)
oa2 = spectrum_oa(cl2, mat)


# Juno.@enter spectrum_oa(cl1, mat)



# aa = DataFrame(:wavelength => mat.wavelength, :value => oa1.average.absorption)
# @vlplot(data=aa,
# width= 400,
# height =  300,
#     mark = {:line},
#     encoding = {x = "wavelength:q", y = "value:q"}
# )
#



function oa_df(x, wavelength)

    a = DataFrame(:extinction => x.average.extinction,
            :absorption => x.average.absorption,
            :scattering => x.average.scattering,
            :wavelength => wavelength,
            :type => "average")

        d = DataFrame(:extinction => x.dichroism.extinction,
                :absorption => x.dichroism.absorption,
                :scattering => x.dichroism.scattering,
                :wavelength => wavelength,
                :type => "dichroism")

    stack([a;d], Not([:wavelength,:type]))

end

d3 = oa_df(oa1, mat.wavelength)
d4 = oa_df(oa2, mat.wavelength)

d5 = [insertcols!(d3, :cluster => "dimer");
     insertcols!(d4, :cluster => "single")]

# filter(row -> row[:type] == "average" , d5)
d5 |> @vlplot(
 width= 400,
 height =  300,
     mark = {:line},
     row = "type",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", color = "variable:n", strokeDash="cluster:n"}
 )
