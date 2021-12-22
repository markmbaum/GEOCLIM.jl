module GEOCLIM

using NetCDF
using UnPack
using MultiAssign
using Roots
using BasicInterpolators
using StaticArrays: SVector
using LinearAlgebra: ⋅ #dot product

#------------------------------------------------------------------------------
# immutable physical constants

#ideal gas constant [J/K*mole]
const 𝐑 = 8.31446262

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#seconds in a year
const 𝐲𝐫 = 31536000.0

#------------------------------------------------------------------------------
# general functions

export cellarea, landfraction

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*abs(Δϕ*(cos(θ₁) - cos(θ₂)))

#computes land fraction of a topography file
#assumes latitude ∈ [-90, 90]°
#assumes land is where topo > 0
function landfraction(fn::String;
                      latname::String="lat", #variable name
                      lonname::String="lon", #variable name
                      toponame::String="topo", #variable name
                      cut::Real=Inf) #restrict to cells where -cutlat <= lat <= lat
    #read variables from file
    lat = ncread(fn, latname)
    lon = ncread(fn, lonname)
    topo = ncread(fn, toponame)
    #transpose if needed
    n, m = length(lat), length(lon)
    if size(topo) == (m,n)
        topo = collect(transpose(topo))
    else
        @assert size(topo) == (n,m)
    end
    #compute
    L = 0.0
    A = 0.0
    @inbounds for i ∈ 1:n
        if -cut <= lat[i] <= cut
            #cell weight depends on latitude
            w = cos(lat[i]*π/180)
            for j ∈ 1:m
                A += w
                if topo[i,j] > 0
                    L += w
                end
            end
        end
    end
    return L/A
end

#------------------------------------------------------------------------------
# types

include("Climatology.jl")
include("ClimatologyInterpolator.jl")

#------------------------------------------------------------------------------
# weathering functions producing weathering per unit area

export godderis, whak, mac

#==============================================================================
As implementated by
  * Goddéris, Y. et al. Onset and ending of the late Palaeozoic ice age triggered by tectonically paced rock weathering. Nature Geosci 10, 382–386 (2017).
  * Donnadieu, Y. et al. A GEOCLIM simulation of climatic and biogeochemical consequences of Pangea breakup: SIMULATION OF PANGEA BREAKUP. Geochem. Geophys. Geosyst. 7, n/a-n/a (2006).
and in other, similar papers.

Arguments
  r - runoff [m/s]
  T - temperature [K]
  k - calibration constant [mole/m^3]
  Eₐ - activation energy [J/mole]
  T₀ - reference temperature [K]
==============================================================================#
godderis(r, T, k, Eₐ, T₀) = k*r*exp((Eₐ/𝐑)*(1/T₀ - 1/T))

#==============================================================================
Originally named after
  * Walker, J. C. G., Hays, P. B. & Kasting, J. F. A negative feedback mechanism for the long-term stabilization of Earth’s surface temperature. J. Geophys. Res. 86, 9776 (1981).
As implementated by
  * Abbot, D. S., Cowan, N. B. & Ciesla, F. J. Indication of Insensitivity of Planetary Weathering Behavior and Habitable Zone to Surface Land Fraction. ApJ 756, 178 (2012).
The important difference between this function and the godderis function is the
pCO2 dependence. The temperature dependence is also written differently.

Arguments
  r - runoff [m/s]
  T - temperature [K]
  pCO2 - carbon dioxide concentration [units must be same as pCO2₀]
  k - calibration constant [mole/m^3]
  Tₑ - scaling of temperature dependence [K]
  T₀ - reference temperature [K]
  pCO2₀ - reference carbon dioxide concentration [units must be same as pCO2]
==============================================================================#
whak(r, T, pCO2, k, Tₑ, T₀, pCO2₀, β=0.2) = k*r*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β

#==============================================================================
As implementated by
  * Maher, K. & Chamberlain, C. P. Hydrologic Regulation of Chemical Weathering and the Geologic Carbon Cycle. Science 343, 1502–1504 (2014).
  * Graham, R. J. & Pierrehumbert, R. Thermodynamic and Energetic Limits on Continental Silicate Weathering Strongly Impact the Climate and Habitability of Wet, Rocky Worlds. ApJ 896, 115 (2020).

Arguments
  r - runoff [m/s]
  T - temperature [K]
  pCO2 - carbon dioxide concentration [units must be same as pCO2₀]
  Tₑ - scaling of temperature dependence [K]
  T₀ - reference temperature [K]
  pCO2₀ - reference carbon dioxide concentration [units must be same as pCO2]
Paramters
  n - thermodynamic pCO2 dependence [-]
  Λ - thermodynamic coefficient for Ceq [-]
  L - flow path length [m] 
  ϕ - porosity [-]
  ρ - mineral mass to fluid volume ratio [kg m⁻³]
  k₀ - reference rate constant [mol m⁻² yr⁻¹]
  𝐀 - specific surface area (not weathering surface area) [m²kg⁻¹]
  X - reactive mineral conc. in fresh rock [-]
  tₛ - soil age [yr]
  m - mineral molar mass [kg/mol]
  μ - scaling constant [-]
  β - pCO2 scaling [-]
==============================================================================#
function mac(r, T, pCO2, Tₑ, T₀, pCO2₀;
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
             β=0.2)
    #defined for convenience
    α = L*ϕ*ρ*𝐀*X*μ
    #equilibrium concentration
    Ceq = 1e3*Λ*(pCO2^n) #conversion from mol/liter to mol/m3
    #temperature dependence
    a = exp((T - T₀)/Tₑ)
    #pCO2 dependence
    b = (pCO2/pCO2₀)^β
    #denominator
    d = 1/(k₀*a*b) + m*𝐀*tₛ + α/(r*𝐲𝐫*Ceq)
    #weathering per unit area 
    (α/d)/𝐲𝐫
end

#--------------------------------------
#each total weathering function has a similar form, can generalize

function totalweathering(𝒻w::F,
                         𝒸::Climatology,
                         args::Vararg{Real,N};
                         kwargs...
                         ) where {F<:Function, N}
    #get climatology fields
    @unpack mask, r, T, A, f, n, m = 𝒸
    #repackage varargs
    X = ntuple(n->Float64(args[n]), N)
    #initialize the weathering sum
    ΣW = 0.0
    #sum weathering at all grid cells
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            #area of cell [m^2]
            Aᵢⱼ = A[i,j]
            #land fraction of cell [-]
            fᵢⱼ = f[i,j]
            #total weathering in cell
            ΣW += Aᵢⱼ*fᵢⱼ*𝒻w(r[i,j], T[i,j], X...; kwargs...)
        end
    end
    return ΣW
end

#--------------------------------------
#total weathering wrappers for each method to be applied to a Climatology

function godderis(𝒸::Climatology, k, Eₐ, Tₒ)
    totalweathering(godderis, 𝒸, k, Eₐ, Tₒ)
end

function whak(𝒸::Climatology, pCO2, k, Tₑ, T₀, pCO2₀, β=0.2)
    totalweathering(whak, 𝒸, pCO2, k, Tₑ, T₀, pCO2₀, β)
end

function mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀; kwargs...)
    totalweathering(mac, 𝒸, pCO2, Tₑ, T₀, pCO2₀; kwargs...)
end

#------------------------------------------------------------------------------

#====
This is a general function to perform root finding with a
ClimatologyInterpolator. The function will find the interpolation location x
where 𝒻(x, ℐ(x)), a function applied to a Climatology, equals the value y.
Typically the function 𝒻 would be a weathering function, but it could be
anything at all.
====#

export findequilibrium

function findequilibrium(ℐ::ClimatologyInterpolator{I,𝒯},
                         𝒻::F,
                         𝓎::Real;
                         tol::Real=1e-4,
                         maxevals::Int=1000
                         ) where {I,𝒯,F}
    #===
    The function to zero is the difference between
    an operation on a Climatology (like a weathering
    estimate) and the desired value of that operation.
    The returned value is converted to the same numeric
    type as the interpolator.
    ===#
    ℱ(𝓍) = convert(𝒯, 𝒻(𝓍, ℐ(𝓍)) - 𝓎)
    #the limits of the ClimatologyInterpolator's range
    𝓍₁, 𝓍₂ = ℐ.x[1], ℐ.x[end]
    #find the root with a bracketing method
    return find_zero(
        ℱ,
        (𝓍₁,𝓍₂),
        Roots.Brent(),
        atol=convert(𝒯, tol),
        rtol=convert(𝒯, tol),
        xatol=convert(𝒯, tol),
        xrtol=convert(𝒯, tol),
        maxevals=maxevals
    )
end

end
