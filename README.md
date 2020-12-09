# DIVAnd_hfradar.jl

[![Build Status](https://github.com/gher-ulg/DIVAnd_hfradar.jl/workflows/CI/badge.svg)](https://github.com/gher-ulg/DIVAnd_hfradar.jl/actions)
[![Coverage Status](https://coveralls.io/repos/gher-ulg/DIVAnd_hfradar.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/gher-ulg/DIVAnd_hfradar.jl?branch=master)
[![codecov.io](http://codecov.io/github/gher-ulg/DIVAnd_hfradar.jl/coverage.svg?branch=master)](http://codecov.io/github/gher-ulg/DIVAnd_hfradar.jl?branch=master)
[![documentation latest](https://img.shields.io/badge/docs-dev-blue.svg)](https://gher-ulg.github.io/DIVAnd_hfradar.jl/dev/)

The package `DIVAnd_hfradar` interpolates surface current data on a regular grid possibly taking dynamical contraints into account.
The primary use-case is for radial current measurements for high-frequency radars. But it can also be applied to any other
current data (like ADCPs or drifters).

The method is described in: Alexander Barth, Charles Troupin, Emma Reyes, Aida Alvera-Azcárate, Jean-Marie Beckers and Joaquı́n Tintoré (2020): Variational interpolation of high-frequency radar surface currents using DIVAnd. Ocean Dynamics (in press)

# Installation

Install DIVAnd_hfradar.jl in julia 1.5 or later with the folling command executed in Julia:

```julia
using Pkg
Pkg.add(url="https://github.com/gher-ulg/DIVAnd_hfradar.jl", rev="master")
```

# Documentation

[Documentation is available here](https://gher-ulg.github.io/DIVAnd_hfradar.jl/dev/).
