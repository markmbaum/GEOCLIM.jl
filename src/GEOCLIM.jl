module GEOCLIM

using NetCDF
using UnPack
using BasicInterpolators

#------------------------------------------------------------------------------
# immutable physical constants

#ideal gas constant [J/K*mole]
const 𝐑 = 8.31446262

#the Earth's mean radius [m]
const 𝐑ₑ = 6.371e6

#------------------------------------------------------------------------------
export Climatology, weathering, weathering!, totalweathering

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*Δϕ*(cos(θ₁) - cos(θ₂))

struct Climatology
    mask::BitMatrix #land mask (1 for land, 0 otherwise)
    r::Matrix{Float64} #runoff [m/s]
    T::Matrix{Float64} #temperature [K]
    A::Matrix{Float64} #area [m^2]
    n::Int64 #rows
    m::Int64 #columns
end

function Base.show(io::IO, 𝒞::Climatology) where {F}
    @unpack mask, r, T, A, n, m= 𝒞
    print(io, "$n x $m Climatology\n")
    Tmax = round(maximum(T[mask]), sigdigits=4)
    Tmin = round(minimum(T[mask]), sigdigits=4)
    print(io, "  temperature ∈ [$Tmin, $Tmax] K\n")
    rmax = round(maximum(r[mask]), sigdigits=4)
    rmin = round(minimum(r[mask]), sigdigits=4)
    print(io, "  runoff ∈ [$rmin, $rmax] m/s\n")
    f = round(sum(A[mask])/sum(A), sigdigits=4)
    print(io, "  land fraction = $f")
end

function Climatology(fnᵣ::String, #runoff file name
                     vᵣ::String,  #runoff variable name
                     fᵣ::Real,    #runoff empty/fill value
                     cᵣ::Real,    #runoff conversion factor
                     fnₜ::String,  #temperature file name
                     vₜ::String)   #temperature variable name

    #first read the runoff file
    r = ncread(fnᵣ, vᵣ)
    #squash third dimension
    r = dropdims(r, dims=3)
    #replace missing data with NaNs
    fᵣ = convert(eltype(r), fᵣ)
    r[r .== fᵣ] .= NaN
    #also replace negative values
    r[r .< 0] .= NaN
    #make a mask for future calculations
    mask = r .|> isnan .|> !
    #convert units
    r .*= cᵣ

    #read temperature file
    T = ncread(fnₜ, vₜ)
    #squash third dimension
    T = dropdims(T, dims=3)
    #demand the same shape as r
    @assert size(r) == size(T) "runoff and temperature grids have different sizes"

    #calculate cell areas, assuming equal spacing in lat & lon
    n, m = size(T)
    Δθ = π/n
    Δϕ = 2π/m
    A = similar(T)
    for i ∈ 1:n
        θ₁ = (i-1)*Δθ
        θ₂ = θ₁ + Δθ
        A[i,:] .= cellarea(𝐑ₑ, Δϕ, θ₁, θ₂)
    end

    #construct
    Climatology(
        mask,
        convert(Matrix{Float64}, r),
        convert(Matrix{Float64}, T),
        convert(Matrix{Float64}, A),
        n,
        m
    )
end

weathering(r, T, A, k, Eₐ, T₀) = k*r*A*exp((Eₐ/𝐑)*(1/T₀ - 1/T))

weathering(𝒞::Climatology, k, Eₐ, T₀) = weathering.(𝒞.r, 𝒞.T, 𝒞.A, k, Eₐ, T₀)

function totalweathering(𝒞::Climatology, k, Eₐ, T₀)
    @unpack mask, r, T, A, n, m = 𝒞
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering(r[i,j], T[i,j], A[i,j], k, Eₐ, T₀)
        end
    end
    return ΣW
end

end