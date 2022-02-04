# GEOCLIM.jl 🌋🌎🪨

[![Build Status](https://github.com/markmbaum/GEOCLIM.jl/workflows/CI/badge.svg)](https://github.com/markmbaum/GEOCLIM.jl/actions)
[![Coverage](https://codecov.io/gh/markmbaum/GEOCLIM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/markmbaum/GEOCLIM.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5784232.svg)](https://doi.org/10.5281/zenodo.5784232)

## Global Silicate Weathering Estimation

This module replicates some features of the [GEOCLIM](https://geoclimmodel.wordpress.com/) model, originally written in Fortran, but now in Julia to make them easier to use. It also extends the original weathering equations, now including
* WHAK (named after [Walker, Hayes & Kasting](http://doi.wiley.com/10.1029/JC086iC10p09776))
	1. Ignoring direct dependence on pCO2, as seen in [Goddéris et al.](https://doi.org/10.1038/ngeo2931), [Donnadieu et al.](https://doi.org/10.1029/2006GC001278), and elsewhere
	2. Including direct pCO2 dependence, as seen in [Abbot et al.](https://doi.org/10.1088/0004-637X/756/2/178) and elsewhere
* MAC (named after [Maher & Chamberlain](https://doi.org/10.1126/science.1250770)) as seen in [Graham & Pierrehumbert](https://doi.org/10.3847/1538-4357/ab9362)

The module implements these formulations to estimate global silicate weathering rates from gridded climatology, typically taken from the results of a global climate model like [CCSM](https://www.cesm.ucar.edu/models/ccsm4.0/) or [FOAM](https://www.mcs.anl.gov/research/projects/foam/). It is intended to estimate weathering during periods of Earth history when the continental configuration was radically different, typically more than 100 million years ago. For more information about the original GEOCLIM, see the Methods/Supplement of [Goddéris et al.](https://doi.org/10.1038/ngeo2931)

------
### Install

The module won't be put in Julia's general registry. If you want to easily install, add it with the package manager using the url.
```
julia> ] add https://github.com/markmbaum/GEOCLIM.jl
```

It's recommended that you do so in an [environment](https://pkgdocs.julialang.org/v1.2/environments/).

------
#### Weathering Functions

There are three weathering functions corresponding to the formulations listed above
1. `godderis(r, T, k, Eₐ, T₀)`
2. `whak(r, T, pCO2, k, Tₑ, T₀, pCO2₀, β=0.2)`
3. `mac(r, T, pCO2, Tₑ, T₀, pCO2₀;
             n=0.316,
             Λ=1.4e-3,
             L=1.0,
             ϕ=0.1,
             ρ=12728.0,
             k₀=8.7e-6,
             𝐀=1e2,
             X=0.36,
             tₛ=1e5,
             m=0.27,
             μ=exp(2),
             β=0.2)`
             
where `r` is runoff and `T` is temperature. You can find explanations of all the other arguments and their units in the [main source file](https://github.com/markmbaum/GEOCLIM.jl/blob/main/src/GEOCLIM.jl) or the referenced papers. These primary functions have no type restrictions.

------
#### Climatologies

The rest of the package is focused on two dimensional grids of results from GCM simulations and is structured around the `Climatology` type. Read GCM results into a `Climatology` by calling the constructor
```
Climatology(fnr,   #runoff file name
            vr,    #runoff variable name
            nullr, #runoff value interpreted as null
            convr, #conversion factor applied to runoff values
            fnT,   #temperature file name
            vT,    #temperature variable name
            fnf,   #land fraction file name
            vf;    #land fraction variable name
            fnlat="",      #file where cell latitudes can be found (empty will use the runoff file)
            latname="lat") #name of latitude vector
```
where the file names point to NetCDF files.

Then each of the weathering functions can be called on a `Climatology` by passing the struct instead of `r` and `T`, returning a global sum of weathering in each grid cell.

For example, to compute the global `mac` weathering estimate with a climatology variable called `clim`,
```
mac(clim, pCO2, Tₑ, T₀, pCO2₀)
```
The temperature and runoff are already in the struct, so you only need to provide the other weathering arguments.

------
#### Climatology Interpolation and Weathering Equilibria

Multiple climatologies can be linked into a `ClimatologyInterpolator` to easily perform cell-wise interpolation. The constructor is
```
ClimatologyInterpolator(𝒞::AbstractVector{Climatology}, x::AbstractVector{<:Real})
```
where `x` contains the independent variable you want to interpolate over. For example, if you have several simulations that are identical except for the atmospheric CO2 concentration, `x` could be `log10(pCO2)`.

Then you can find the root of some function applied to interpolated climatology by defining the function `f(C,x)`, where `C` is  a climatology, and passing it to the `findequilibrium` function. For example, to find the CO2 concentration where global weathering balances some value, you could
```julia
#make an interpolator, assuming you already have three Climatology structs
logpCO2 = [1., 2., 3.]
I = ClimatologyInterpolator(C, logpCO2) #C is a Vector{Climatology}

#define the weathering operation using mac, where x is log10(pCO2)
w(x,c) = mac(c, exp10(x)*1e-6, 11.1, 288.15, 285e-6)

#find the log10(CO2) value where weathering balances 5e4 moles/second
findequilibrium(I, w, 5e4)
```
