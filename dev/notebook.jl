### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ ed601426-53a6-4bfc-9f1c-53ebe4a95716
using Markdown

# ╔═╡ 2e80aa67-7441-4aab-ae6a-266aed6c5fb8
using CoupledDipole

# ╔═╡ 867ae0d8-fb4a-4b5b-91ef-5cb4c0a382f0
using DataFrames

# ╔═╡ 545b82ad-0e48-4b42-b2c2-884c4e3e4b2b
using VegaLite

# ╔═╡ 29823b36-159a-48ea-831a-30082315cf95
using PlutoUI

# ╔═╡ 5a561f4c-080b-415d-a21e-8e89a6a6bf50
using InteractiveUtils

# ╔═╡ 02d9a90c-a22d-4531-b0df-528497b8e29e
md"""
# Coupled dipole model of a dimer of Au spheroids

We consider 2 spheroids separated by a distance d, and oriented with a dihedral angle $\varphi$ relative to each other.
"""

# ╔═╡ 310c5403-601b-47bc-9be5-ab10fe682058
@bind d Slider(50:200, default = 80)

# ╔═╡ d8cf9cdc-7ee1-42c2-b5f0-34b46d40fbf7
@bind alpha Slider(0:90, default = 45)

# ╔═╡ 780b6df4-541d-4add-8ff9-5b27ee705251
begin
	wavelength = collect(450:5:850)
    medium = Dict([("Au", epsilon_Au), ("medium", x -> 1.33)])
    mat = Material(wavelength, medium)
end

# ╔═╡ 3411ce46-78ef-4e5c-8ecb-f8a2a347b1b4
cl = cluster_dimer(d, 20, 20, 40, alpha*pi/180)

# ╔═╡ 23d7feab-ef69-40ac-aadd-7e2fb77707cb
oa = spectrum_oa(cl, mat, "gl", 36)

# ╔═╡ 560a64f1-05ca-4e5e-ad05-68dee37704a4
function oa_df(x, wavelength, cluster = "dimer")

    a = DataFrame(:extinction => x.average.extinction,
            :absorption => x.average.absorption,
            :scattering => x.average.scattering,
            :wavelength => wavelength,
            :type => "average",
		    :cluster => cluster
	)

        d = DataFrame(:extinction => x.dichroism.extinction,
                :absorption => x.dichroism.absorption,
                :scattering => x.dichroism.scattering,
                :wavelength => wavelength,
                :type => "dichroism",
		        :cluster => cluster
	 )

    stack([a;d], Not([:wavelength,:type,:cluster]))

end

# ╔═╡ 0e312bcd-6360-4713-a38c-5ebb0e152641
s = oa_df(oa, mat.wavelength, "dimer")

# ╔═╡ 2b6d88ff-b6e6-42eb-b418-0bd792c2dd27
begin
	sref = spectrum_oa(cluster_single(20, 20, 40), mat, "gl", 36)
    ref = oa_df(sref, mat.wavelength, "single")
end

# ╔═╡ ecc47ec8-8210-4237-9287-d898ae15185a
both = [s ; ref]

# ╔═╡ d61daf28-52da-4589-9009-2946365ca543
both |> @vlplot(
 width = 400,
 height = 300,
     mark = {:line},
     row = "type",
     resolve={scale={y="independent"}},
     encoding = {x = "wavelength:q", y = "value:q", 
	             color = "variable:n", strokeDash="cluster:n"}
 )


# ╔═╡ Cell order:
# ╠═ed601426-53a6-4bfc-9f1c-53ebe4a95716
# ╠═2e80aa67-7441-4aab-ae6a-266aed6c5fb8
# ╠═867ae0d8-fb4a-4b5b-91ef-5cb4c0a382f0
# ╠═545b82ad-0e48-4b42-b2c2-884c4e3e4b2b
# ╠═29823b36-159a-48ea-831a-30082315cf95
# ╠═5a561f4c-080b-415d-a21e-8e89a6a6bf50
# ╟─02d9a90c-a22d-4531-b0df-528497b8e29e
# ╟─310c5403-601b-47bc-9be5-ab10fe682058
# ╟─d8cf9cdc-7ee1-42c2-b5f0-34b46d40fbf7
# ╠═d61daf28-52da-4589-9009-2946365ca543
# ╠═780b6df4-541d-4add-8ff9-5b27ee705251
# ╠═3411ce46-78ef-4e5c-8ecb-f8a2a347b1b4
# ╠═23d7feab-ef69-40ac-aadd-7e2fb77707cb
# ╠═560a64f1-05ca-4e5e-ad05-68dee37704a4
# ╠═0e312bcd-6360-4713-a38c-5ebb0e152641
# ╠═2b6d88ff-b6e6-42eb-b418-0bd792c2dd27
# ╠═ecc47ec8-8210-4237-9287-d898ae15185a
