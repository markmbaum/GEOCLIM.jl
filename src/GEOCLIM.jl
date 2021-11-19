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

#------------------------------------------------------------------------------
# types

include("Climatology.jl")
include("ClimatologyInterpolator.jl")

#------------------------------------------------------------------------------
# weathering functions

export weathering, totalweathering, weathering_mac, totalweathering_mac

#--------------------------------------
# As implementated by Goddéris et al. 2017

weathering(r, T, A, k, Eₐ, T₀) = k*r*A*exp((Eₐ/𝐑)*(1/T₀ - 1/T))

function weathering(𝒸::Climatology, k, Eₐ, T₀)
    weathering.(𝒸.r, 𝒸.T, 𝒸.A, k, Eₐ, T₀)
end

function totalweathering(𝒸::Climatology, k, Eₐ, T₀)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering(r[i,j], T[i,j], A[i,j], k, Eₐ, T₀)
        end
    end
    return ΣW
end

#--------------------------------------
# As implementated by Abbot et al. 2012 
# pCO2 dependence is added and the temperature dependence is slightly different
const β  = 0.2

function weathering(r, T, A, pCO2, k, Tₑ, T₀, pCO2₀) 
    k*r*A*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β
end

function weathering(𝒸::Climatology, pCO2, k, Tₑ, T₀, pCO2₀)
    weathering.(𝒸.r, 𝒸.T, 𝒸.A, pCO2, k, Tₑ, T₀, pCO2₀)
end

function totalweathering(𝒸::Climatology, pCO2, k, Tₑ, T₀, pCO2₀)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering(r[i,j], T[i,j], A[i,j], pCO2, k, Tₑ, T₀, pCO2₀)
        end
    end
    return ΣW
end

#------------------------------------------------------------------------------
# MAC, as implementated by Graham and Pierrehumbert 2020 
# following Maher and Chamberlin 2014

const n = 0.316 # Thermodynamic pCO2 dependence [-]
const Λ = 1.4e-3 # Thermodynamic coefficient for Ceq [-]
const L = 1 # Flow path length [m] 
const ϕ = 0.1 # Porosity [-]
const ρ = 12728 # Mineral mass to fluid volume ratio [kg m⁻³]
const k₀ = 8.7e-6 # Reference rate constant [mol m⁻² yr⁻¹]
const 𝐀 = 100 # Specific surface area [m²kg⁻¹]
const X = 0.36 # Reactive mineral conc. in fresh rock [-]
const tₛ = 1e5 # Soil age [yr]
const m = 0.27 # Mineral molar mass [kg/mol]
const μ = ℯ^2 # Scaling constant [-]
const α = L*ϕ*ρ*𝐀*X*μ # Defined for convenience [-]
const s_y = 31536000

function Ceq(pCO2)
    return Λ*(pCO2)^n*1000 #conversion from mol/liter to mol/m3, ppm to bar
end

# r input in m/s, convert to m/yr, convert result from mol/y back to mol/s
function weathering_mac(r, T, A, pCO2, Tₑ, T₀, pCO2₀) 
    A*α*((k₀*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β)^-1 + m*𝐀*tₛ + α/(r*s_y*Ceq(pCO2)))^-1/s_y
end

function weathering_mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀) 
    weathering_mac.(𝒸.r, 𝒸.T, 𝒸.A, pCO2, Tₑ, T₀, pCO2₀)
end

function totalweathering_mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering_mac(r[i,j], T[i,j], A[i,j], pCO2, Tₑ, T₀, pCO2₀)
        end
    end
    return ΣW
end


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