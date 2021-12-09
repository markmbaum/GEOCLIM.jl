# GEOCLIM.jl 🌎

[![Build Status](https://github.com/markmbaum/GEOCLIM.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/GEOCLIM.jl/actions)
[![Coverage](https://codecov.io/gh/markmbaum/GEOCLIM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/markmbaum/GEOCLIM.jl)

### Global Silicate Weathering Estimation

This module replicates some features of the [GEOCLIM](https://geoclimmodel.wordpress.com/) model, originally written in Fortran, but now in Julia to make them easier to use. It also extends the original weathering equations, now including
* WHAK (named after [Walker, Hayes & Kasting](http://doi.wiley.com/10.1029/JC086iC10p09776))
	1. Ignoring direct dependence on pCO2, as seen in [Goddéris et al.](https://doi.org/10.1038/ngeo2931), [Donnadieu et al.](https://doi.org/10.1029/2006GC001278), and elsewhere
	2. Including direct pCO2 dependence, as seen in [Abbot et al.](https://doi.org/10.1088/0004-637X/756/2/178) and elsewhere
* MAC (named after [Maher & Chamberlain](https://doi.org/10.1126/science.1250770)) as seen in [Graham & Pierrehumbert](https://doi.org/10.3847/1538-4357/ab9362)

------

The module implements these formulations to estimate global silicate weathering rates from gridded climatology, typically taken from the results of a Global Climate Model like [CCSM](https://www.cesm.ucar.edu/models/ccsm4.0/) or [FOAM](https://www.mcs.anl.gov/research/projects/foam/). It is intended to estimate weathering during periods of Earth history when the continental configuration was radically different, typically more than 100 million years ago. For more information about the original GEOCLIM, see the Methods/Supplement of [Goddéris et al.](https://doi.org/10.1038/ngeo2931)

------

The module won't be put in Julia's general registry. If you want to easily install, add it with the package manager using the url.
```
julia> ] add https://github.com/markmbaum/GEOCLIM.jl
```
