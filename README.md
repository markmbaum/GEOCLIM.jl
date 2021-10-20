# GEOCLIM.jl 🌎

*work in progress*

[![Build Status](https://github.com/markmbaum/GEOCLIM.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/GEOCLIM.jl/actions)
[![Coverage](https://codecov.io/gh/markmbaum/GEOCLIM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/markmbaum/GEOCLIM.jl)

This repo is intended to replicate some features of the [GEOCLIM](https://geoclimmodel.wordpress.com/) model, originally written in Fortran, but now in Julia to make them easier to use. The model estimates global silicate weathering rates from gridded climatology, typically taken from the results of a Global Climate Model like [CCSM](https://www.cesm.ucar.edu/models/ccsm4.0/) or [FOAM](https://www.mcs.anl.gov/research/projects/foam/). It is intended to estimate weathering during periods of Earth history when the continental configuration was radically different that today, typically more than 100 million years ago. For more information about the methodology of GEOCLIM, see the first link above. A good description is also found in the Methods/Supplement of

* [Goddéris, Y., Donnadieu, Y., Carretier, S. et al. *Onset and ending of the late Palaeozoic ice age triggered by tectonically paced rock weathering.* Nature Geoscience 10, 382–386 (2017)](https://doi.org/10.1038/ngeo2931)