module GEOCLIM

using NetCDF
using UnPack
using MultiAssign
using Roots
using BasicInterpolators

#------------------------------------------------------------------------------
# immutable physical constants

#ideal gas constant [J/K*mole]
const 𝐑 = 8.31446262

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#seconds in a year
const 𝐲𝐫 = 31536000.0

#------------------------------------------------------------------------------
# types

include("Climatology.jl")
include("ClimatologyInterpolator.jl")

#------------------------------------------------------------------------------
# weathering functions

#--------------------------------------
# As implementated by Goddéris et al. 2017 (and in previous papers)

export godderis

godderis(r, T, A, k, Eₐ, T₀) = k*r*A*exp((Eₐ/𝐑)*(1/T₀ - 1/T))

function godderis(𝒸::Climatology, k, Eₐ, T₀)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += godderis(r[i,j], T[i,j], A[i,j], k, Eₐ, T₀)
        end
    end
    return ΣW
end

#--------------------------------------
# As implementated by Abbot et al. 2012 and Graham and Pierrehumbert 2020 
# pCO2 dependence is added and the temperature dependence is slightly different from godderis (original geoclim)

export whak

function whak(r, T, A, pCO2, k, Tₑ, T₀, pCO2₀, β=0.2)
    k*r*A*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β
end

function whak(𝒸::Climatology, pCO2, k, Tₑ, T₀, pCO2₀, β=0.2)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += whak(r[i,j], T[i,j], A[i,j], pCO2, k, Tₑ, T₀, pCO2₀, β)
        end
    end
    return ΣW
end

#------------------------------------------------------------------------------
# MAC, as implementated by Graham and Pierrehumbert 2020 
# following Maher and Chamberlin 2014

export mac

# r input in m/s, convert to m/yr, convert result from mol/y back to mol/s
function mac(r, T, A, pCO2, Tₑ, T₀, pCO2₀;
             n=0.316, # Thermodynamic pCO2 dependence [-]
             Λ=1.4e-3, # Thermodynamic coefficient for Ceq [-]
             L=1, # Flow path length [m] 
             ϕ=0.1, # Porosity [-]
             ρ=12728, # Mineral mass to fluid volume ratio [kg m⁻³]
             k₀=8.7e-6, # Reference rate constant [mol m⁻² yr⁻¹]
             𝐀=100, # Specific surface area, not weathering surface area [m²kg⁻¹]
             X=0.36, # Reactive mineral conc. in fresh rock [-]
             tₛ=1e5, # Soil age [yr]
             m=0.27, # Mineral molar mass [kg/mol]
             μ=exp(2), # Scaling constant [-]
             β=0.2) # pCO2 scaling [-]
    #defined for convenience
    α = L*ϕ*ρ*𝐀*X*μ
    #equilibrium concentration
    Ceq = 1e3*Λ*pCO2^n #conversion from mol/liter to mol/m3, ppm to bar
    #weathering
    A*α*((k₀*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β)^-1 + m*𝐀*tₛ + α/(r*𝐲𝐫*Ceq))^-1/𝐲𝐫
end

function mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += mac(r[i,j], T[i,j], A[i,j], pCO2, Tₑ, T₀, pCO2₀)
        end
    end
    return ΣW
end

#------------------------------------------------------------------------------

#====
This is a general function to perform root finding with a
ClimatologyInterpolator. The function will find the interpolation location x
where 𝒻(ℐ(x)), a function applied to a Climatology, equals the value y.
====#

export findequilibrium

function findequilibrium(ℐ::ClimatologyInterpolator,
                         𝒻::F,
                         𝓎::Real;
                         tol::Float64=1e-4,
                         maxevals::Int=1000
                         ) where {F}
    #===
    The function to zero is the difference between
    an operation on a Climatology (like a weathering
    estimate) and the desired value of that operation.
    ===#
    ℱ(𝓍) = 𝒻(ℐ(𝓍)) - 𝓎
    #the limits of the ClimatologyInterpolator's range
    𝓍₁, 𝓍₂ = ℐ.x[1], ℐ.x[end]
    #find the root with a bracketing method
    return find_zero(
        ℱ,
        (𝓍₁,𝓍₂),
        Roots.Brent(),
        atol=tol,
        rtol=tol,
        xatol=tol,
        xrtol=tol,
        maxevals=maxevals
    )
end

end