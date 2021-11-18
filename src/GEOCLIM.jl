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

export weathering, totalweathering

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

function weathering(r, T, A, pCO2, k, Eₐ, T₀, pCO2₀, β) 
    k*r*A*exp((Eₐ/𝐑)*(T - T₀)/T₀^2)*(pCO2/pCO2₀)^β
end

function weathering(𝒸::Climatology, pCO2, k, Eₐ, T₀, pCO2₀, β)
    weathering.(𝒸.r, 𝒸.T, 𝒸.A, pCO2, k, Eₐ, T₀, pCO2₀, β)
end

function totalweathering(𝒸::Climatology, pCO2, k, Eₐ, T₀, pCO2₀, β)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering(r[i,j], T[i,j], A[i,j], pCO2, k, Eₐ, T₀, pCO2₀, β)
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