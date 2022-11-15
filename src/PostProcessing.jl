

using DataFrames

"""
    dispersion_df(x, wavelength)

Long format dataframe for dispersion spectra

- `x`: CrossSection structure with 3 arrays (extinction, absorption, scattering)
- `wavelength`: vector of wavelengths

"""
function dispersion_df(x, wavelength)

    N_l, N_2a = size(x.extinction)

    res =  DataFrame(wavelength = repeat(wavelength, outer = 3*N_2a),
    value = vcat(vec(x.extinction), vec(x.absorption), vec(x.scattering)),
    polarisation = repeat(["s", "p"], outer = Int(3*N_2a/2), inner = N_l),
    angle = repeat(repeat(1:Int(N_2a/2),outer=2), inner=N_l, outer=3),
    crosstype = repeat(["extinction", "absorption", "scattering"], inner = N_2a*N_l))
    
    return(res)

    # e = insertcols!(
    #     DataFrame(x.extinction,  [:s, :p]),
    #     :wavelength => wavelength,
    #     :crosstype => "extinction"
    # )
    # a = insertcols!(
    #     DataFrame(x.absorption,  [:s, :p]),
    #     :wavelength => wavelength,
    #     :crosstype => "absorption"
    # )

    # s = insertcols!(
    #     DataFrame(x.scattering,  [:s, :p]),
    #     :wavelength => wavelength,
    #     :crosstype => "scattering"
    # )

    # stack([e;a;s], Not([:wavelength,:crosstype]), variable_name="polarisation")

end

"""
    oa_df(x, wavelength)

Long format dataframe for orientation-averaged spectra

- `x`: Tuple with 2 CrossSection structures (average and dichroism),
   each containing 3 arrays (extinction, absorption, scattering)
- `wavelength`: vector of wavelengths

"""
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
