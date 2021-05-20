

using DataFrames

"""
    dispersion_df(x, wavelength)

Long format dataframe for dispersion spectra

- `x`: CrossSection structure with 3 arrays (extinction, absorption, scattering)
- `wavelength`: vector of wavelengths

"""
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
