

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
    value = vcat(vec(x.extinction), vec(x.absorption), vec(x.scattering)), # ext,abs,sca
    polarisation = repeat(["s", "p"], outer = 3, inner = N_l*Int(N_2a/2)), # Nl wavelengths, Na angles | ext,abs,sca
    angle = repeat(repeat(1:Int(N_2a/2),outer=2), inner=N_l, outer=3),     # 2 pol, Nl wavelengths | ext,abs,sca
    crosstype = repeat(["extinction", "absorption", "scattering"], inner = N_2a*N_l))
    
    return(res)

end

"""
    oa_df(x, wavelength)

Long format dataframe for orientation-averaged spectra

- `x`: Tuple with 2 CrossSection structures (average and dichroism),
   each containing 3 arrays (extinction, absorption, scattering)
- `wavelength`: vector of wavelengths

"""
function oa_df(x, wavelength)

    N_l = length(wavelength)

    res =  DataFrame(wavelength = repeat(wavelength, outer = 6),
    value = vcat(x.average.extinction, 
                 x.average.absorption, 
                 x.average.scattering,
                 x.dichroism.extinction,
                 x.dichroism.absorption,
                 x.dichroism.scattering), 
    type = repeat(["average", "dichroism"],inner=3*N_l), 
    crosstype = repeat(["extinction", "absorption", "scattering"], inner = N_l, outer = 2))
    
    return(res)

end
