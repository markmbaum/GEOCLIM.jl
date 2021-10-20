module GEOCLIM

using NetCDF
using BasicInterpolators

#------------------------------------------------------------------------------
# immutable physical constants

#gas constant [J/K*mole], equivalent to kB*Av
const 𝐑 = 8.31446262

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#------------------------------------------------------------------------------
export Climatology, weathering

cellarea(Δϕ, θ₁, θ₂) = Δϕ*(cos(θ₁) - cos(θ₂))

weathering(r, T, A, k, Eₐ, T₀) = k*r*A*exp(-(Eₐ/𝐑)*(1/T - 1/T₀))

weathering(C::Climatology, k, Eₐ, T₀) = weathering.(C.r, C.T, C.A, k, Eₐ, T₀)

struct Climatology{T}
    mask::BitMatrix
    r::Matrix{T}
    T::Matrix{T}
    A::Matrix{T}
end

function Climatology(nc_runoff::String,
                     var_runoff::String,
                     fill_runoff::Real,
                     nc_temperature::String,
                     var_temperature::String)

    #first read the runoff file
    r = ncread(nc_runoff, var_runoff)
    #squash third dimension
    r = dropdims(r, dims=3)
    #replace missing data with NaNs
    f = convert(eltype(r), fill_runoff)
    r[r .== f] .= NaN
    #also replace negative values
    r[r .< 0] .= NaN
    #make a mask to keep
    mask = r .|> isnan .|> !

    #read temperature file
    T = ncread(nc_temperature, var_temperature)
    #squash third dimension
    T = dropdims(T, dims=3)
    #demand the same shape as r
    @assert size(r) == size(T) "runoff and temperature grids have different sizes"

    #calculate cell areas
    n, m = size(T)
    Δθ = π/n
    Δϕ = 2π/m
    A = similar(T)
    for i ∈ 1:n
        θ₁ = (i-1)*Δθ
        θ₂ = θ₁ + Δθ
        A[i,:] .= cellarea(Δϕ, θ₁, θ₂)
    end
    A .*= 𝐑ₑ^2

    #construct
    Climatology(mask, r, T, A)
end

end