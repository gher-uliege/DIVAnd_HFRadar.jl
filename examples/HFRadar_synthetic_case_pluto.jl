### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 1c056adc-2e9a-11eb-3d81-0d29e98a3551
begin
    # We set up a new environment for this notebook
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add("PyCall")
    using PyCall
    try
        pyimport("matplotlib")
    catch
        run(`$(PyCall.python) -m pip install matplotlib`)
    end
    Pkg.add("PlutoUI")
    Pkg.add("DIVAnd")
    Pkg.add("PyPlot")
    Pkg.add(url="https://github.com/gher-ulg/DIVAnd_HFRadar.jl", rev="master")

    using PlutoUI
    using DIVAnd_HFRadar: DIVAndrun_HFRadar
    using DIVAnd
    using PyPlot
end


# ╔═╡ 3e8cc514-2e9a-11eb-35ad-8f110251cb48
begin
    # size of the grid
    sz = (10,11);
    # depth (meters)
    h = 50 * ones(sz);
    # land-sea mask
    # true is sea; false is land
    mask = trues(sz);
    mask[[1, end],:] .= false;
    mask[:,[1, end]] .= false;

    # 2D grid
    xi,yi = DIVAnd.ndgrid(LinRange(-1,1,sz[1]),LinRange(-1,1,sz[2]))

    # scale factor; inverse of the resolution
    pm = ones(sz) / (xi[2,1]-xi[1,1]);
    pn = ones(sz) / (yi[1,2]-yi[1,1]);

    # radial observations
    robs = 1.

    # position of the observation
    xobs = 0.
    yobs = 0.
end;

# ╔═╡ f096b782-2eff-11eb-13a0-2981bb02584f
md"### Illustration of DIVAnd with a single observation"

# ╔═╡ 6d31768a-3af2-11eb-3c08-5533dc6520b8
@bind directionobs Slider(LinRange(0.,360,361),default=45)

# ╔═╡ dcd17c5c-3af3-11eb-0dc8-df26312b5614
md"Orientation of the observation (`directionobs` = $directionobs)"
# direction of the observation (from North counted clockwise)

# ╔═╡ b39c60e8-2e9b-11eb-3c60-7fe7082d59e0
@bind len  Slider(LinRange(0.2,1,101),default=0.6)

# ╔═╡ 43b7c380-2eff-11eb-35dc-9f0a75e6341c
md"Correlation length (`len` = $len)"

# ╔═╡ 4482ce12-2e9c-11eb-028d-2100248f0f9a
@bind epsilon2  Slider(LinRange(0.001,0.1,100),default=0.01)

# ╔═╡ 2c3feeac-2f00-11eb-24bb-2b7015c0551e
md"Uncertainty of the observation (`epsilon2` = $epsilon2)"

# ╔═╡ 2575c772-2e9c-11eb-3948-e3c991d471eb
@bind eps2_boundary_constraint  Slider(LinRange(0.0001,100,100),default=0.0001)

# ╔═╡ 41fd104e-2f00-11eb-2ed7-cd6e393e07f6
md"Uncertainty of the boundary constraint (`eps2_boundary_constraint` = $eps2_boundary_constraint)"

# ╔═╡ 6e9a6224-2e9d-11eb-17b3-af648104f09d
@bind eps2_div_constraint  Slider(LinRange(0.0001,1000,100),default=0.0001)

# ╔═╡ 6668d15c-2f00-11eb-1002-a711787e0257
md"Uncertainty of the divergence constraint (`eps2_div_constraint` = $eps2_div_constraint)"

# ╔═╡ 65d63b30-2e9a-11eb-1d21-6f02c2549c10
begin
    uri,vri = DIVAndrun_HFRadar(
        mask,h,(pm,pn),(xi,yi),([xobs],[yobs]),[robs],[directionobs],
        len,epsilon2,
        eps2_boundary_constraint = eps2_boundary_constraint,
        eps2_div_constraint = eps2_div_constraint,
    )

    figure(figsize = (6,6))
    quiver(xi,yi,uri,vri, scale = 10)
    α = directionobs*pi/180
    quiver(xobs,yobs,robs .* sin.(α), robs .* cos.(α),color = "r",scale = 10)
    contourf(xi,yi,mask,levels = [0,.5],cmap = "gray")
    gcf()
end

# ╔═╡ 1aadac36-3ac9-11eb-04fe-7731a58578b8
md"In the extrem case, where only one observation (red arrow) is available, DIVAnd with the dynamical constraints produce two conter-rotating gyres for a closed domain. Have a look to the [documentation](https://gher-ulg.github.io/DIVAnd_HFRadar.jl/dev/) for more information."

# ╔═╡ Cell order:
# ╟─1c056adc-2e9a-11eb-3d81-0d29e98a3551
# ╟─3e8cc514-2e9a-11eb-35ad-8f110251cb48
# ╟─f096b782-2eff-11eb-13a0-2981bb02584f
# ╟─dcd17c5c-3af3-11eb-0dc8-df26312b5614
# ╟─6d31768a-3af2-11eb-3c08-5533dc6520b8
# ╟─43b7c380-2eff-11eb-35dc-9f0a75e6341c
# ╟─b39c60e8-2e9b-11eb-3c60-7fe7082d59e0
# ╟─2c3feeac-2f00-11eb-24bb-2b7015c0551e
# ╟─4482ce12-2e9c-11eb-028d-2100248f0f9a
# ╟─41fd104e-2f00-11eb-2ed7-cd6e393e07f6
# ╟─2575c772-2e9c-11eb-3948-e3c991d471eb
# ╟─6668d15c-2f00-11eb-1002-a711787e0257
# ╟─6e9a6224-2e9d-11eb-17b3-af648104f09d
# ╟─65d63b30-2e9a-11eb-1d21-6f02c2549c10
# ╟─1aadac36-3ac9-11eb-04fe-7731a58578b8
