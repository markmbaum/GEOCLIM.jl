module GEOCLIM

using NetCDF
using UnPack
using Roots
using ForwardDiff: derivative
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
# useful functions

export weathering, totalweathering, findequilibrium

weathering(r, T, A, k, Eₐ, T₀) = k*r*A*exp((Eₐ/𝐑)*(1/T₀ - 1/T))

weathering(𝒸::Climatology, k, Eₐ, T₀) = weathering.(𝒸.r, 𝒸.T, 𝒸.A, k, Eₐ, T₀)

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

function findequilibrium(ℐ::ClimatologyInterpolator,
                         𝒻::F,
                         y::Real;
                         tol::Float64=1e-3,
                         maxevals::Int=1000
                         ) where {F}
    #===
    The function to zero is the difference between
    an operation on a Climatology (like a weathering
    estimate) and the desired value of that operation.
    ===#
    ℱ(x) = 𝒻(ℐ(x)) - y
    #the limits of the ClimatologyInterpolator's range
    x₁, x₂ = ℐ.x[1], ℐ.x[end]
    #find the root with a bracketing method
    return find_zero(
        ℱ,
        (x₁,x₂),
        Roots.Brent(),
        atol=tol,
        rtol=tol,
        xatol=tol,
        xrtol=tol,
        maxevals=maxevals
    )
end

end