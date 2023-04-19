"""
    dispersion_df(x, wavelength)

Long format dataframe for dispersion spectra

- `x`: CrossSection structure with 4 arrays (extinction, absorption, scattering, ext-sca)
- `wavelength`: vector of wavelengths

"""
function dispersion_df(x, wavelength; format="long", scattering_diff=false)

    N_l, N_2a = size(x.extinction)


    if scattering_diff
        scattering2 = vec(x.extinction) - vec(x.absorption)
        res = DataFrame(wavelength=repeat(wavelength, outer=4 * N_2a),
            value=vcat(vec(x.extinction), vec(x.absorption), vec(x.scattering), scattering2), # ext,abs,sca,sca2
            polarisation=repeat(["pol1", "pol2"], outer=4, inner=N_l * Int(N_2a / 2)), # Nl wavelengths, Na angles | ext,abs,sca
            angle_id=repeat(repeat(1:Int(N_2a / 2), outer=2), inner=N_l, outer=4),     # 2 pol, Nl wavelengths | ext,abs,sca
            crosstype=repeat(["extinction", "absorption", "scattering", "scattering2"], inner=N_2a * N_l))
    else
        res = DataFrame(wavelength=repeat(wavelength, outer=3 * N_2a),
            value=vcat(vec(x.extinction), vec(x.absorption), vec(x.scattering)), # ext,abs,sca
            polarisation=repeat(["pol1", "pol2"], outer=3, inner=N_l * Int(N_2a / 2)), # Nl wavelengths, Na angles | ext,abs,sca
            angle_id=repeat(repeat(1:Int(N_2a / 2), outer=2), inner=N_l, outer=3),     # 2 pol, Nl wavelengths | ext,abs,sca
            crosstype=repeat(["extinction", "absorption", "scattering"], inner=N_2a * N_l))
    end
    if format == "long"
        return (res)
    else
        return (unstack(res, :polarisation, :value))
    end
end

"""
incidence_labels(d, Incidence)
incidence_labels(d, Incidence, labels)

Merging incidence information with dataframe of dispersion spectra

- `d`: long-format output from dispersion_df
- `Incidence`: vector of Rotation objects
- `labels`: vector of labels, to give convenient names to the rotation axes (e.g. ['x', 'y', 'z']) in addition to their cartesian coordinates

"""
function incidence_labels(d, Incidence)
    tmp = map(x -> transpose(Rotations.params(Rotations.AngleAxis(x))), Incidence)
    angles_axes = DataFrame(reduce(vcat, tmp), [:angle, :x, :y, :z])
    angles_axes[!, :angle_id] = 1:nrow(angles_axes)
    # d[!, :angle_id] = 1:nrow(d) # already present from dispersion_df
    return DataFrames.leftjoin(d, angles_axes, on=:angle_id)
end

function incidence_labels(d, Incidence, labels)
    tmp = map(x -> transpose(Rotations.params(Rotations.AngleAxis(x))), Incidence)
    angles_axes = DataFrame(reduce(vcat, tmp), [:angle, :x, :y, :z])
    angles_axes[!, :angle_id] = 1:nrow(angles_axes)
    angles_axes[!, :axis_label] = labels
    # d[!, :angle_id] = 1:nrow(d) # already present from dispersion_df
    return DataFrames.leftjoin(d, angles_axes, on=:angle_id)
end


"""
    oa_df(x, wavelength)

Long format dataframe for orientation-averaged spectra

- `x`: Tuple with 2 CrossSection structures (average and dichroism),
   each containing 4 arrays (extinction, absorption, scattering, ext-sca)
- `wavelength`: vector of wavelengths

"""
function oa_df(x, wavelength, scattering_diff=false)

    N_l = length(wavelength)

    if scattering_diff
        res = DataFrame(wavelength=repeat(wavelength, outer=8),
            value=vcat(x.average.extinction,
                x.average.absorption,
                x.average.scattering,
                x.average.extinction - x.average.absorption,
                x.dichroism.extinction,
                x.dichroism.absorption,
                x.dichroism.scattering,
                x.dichroism.extinction - x.dichroism.absorption),
            type=repeat(["average", "dichroism"], inner=4 * N_l),
            crosstype=repeat(["extinction", "absorption", "scattering", "scattering2"], inner=N_l, outer=2))
    else
        res = DataFrame(wavelength=repeat(wavelength, outer=6),
            value=vcat(x.average.extinction,
                x.average.absorption,
                x.average.scattering,
                x.dichroism.extinction,
                x.dichroism.absorption,
                x.dichroism.scattering),
            type=repeat(["average", "dichroism"], inner=3 * N_l),
            crosstype=repeat(["extinction", "absorption", "scattering"], inner=N_l, outer=2))
    end
    return (res)

end
