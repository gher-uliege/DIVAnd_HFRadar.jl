# DIVAnd_HFRadar.jl

[![Build Status](https://github.com/gher-ulg/DIVAnd_HFRadar.jl/workflows/CI/badge.svg)](https://github.com/gher-ulg/DIVAnd_HFRadar.jl/actions)
[![codecov.io](http://codecov.io/github/gher-ulg/DIVAnd_HFRadar.jl/coverage.svg?branch=master)](http://codecov.io/github/gher-ulg/DIVAnd_HFRadar.jl?branch=master)
[![documentation latest](https://img.shields.io/badge/docs-dev-blue.svg)](https://gher-ulg.github.io/DIVAnd_HFRadar.jl/dev/)

The package `DIVAnd_HFRadar` interpolates surface current data on a regular grid possibly taking dynamical contraints into account.
The primary use-case is for radial current measurements for high-frequency radars (like WERA or CODAR SeaSonde). But it can also be applied to any other
current data (like ADCPs or drifters).

The method is described in: Barth, A., Troupin, C., Reyes, E., Alvera-Azcárate, A., Beckers J.-M. and Tintoré J. (2021): [Variational interpolation of high-frequency radar surface currents using DIVAnd](https://doi.org/10.1007/s10236-020-01432-x). Ocean Dynamics, 71, 293–308
doi: 10.1007/s10236-020-01432-x (open access)

# Installation

Install DIVAnd_HFRadar.jl in [julia](https://julialang.org/downloads/) 1.5 or later with the following command executed in Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/gher-ulg/DIVAnd_HFRadar.jl", rev="master")
```

# Documentation

[Documentation is available here](https://gher-ulg.github.io/DIVAnd_HFRadar.jl/dev/).

# Online-demo

A online-demo is available on the free service binder.org at the following link:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fonsp/pluto-on-binder/master?urlpath=pluto/open?url=https%253A%252F%252Fraw.githubusercontent.com%252Fgher-ulg%252FDIVAnd_HFRadar.jl%252Fmaster%252Fexamples%252FHFRadar_synthetic_case_pluto.jl) (this can take 2 to 10 minutes to start). If this does not work, then you can install `Pluto.jl` in Julia by running:

```julia
using Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```

Then open the file [examples/HFRadar_synthetic_case_pluto.jl](examples/HFRadar_synthetic_case_pluto.jl)

