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
export ClimatologyInterpolator, findequilibrium

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

Base.size(𝒞::Climatology) = (𝒞.n, 𝒞.m)

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

struct ClimatologyInterpolator{I}
    x::Vector{Float64}
    mask::BitMatrix
    r::Matrix{I}
    T::Matrix{I}
    A::Matrix{Float64}
    n::Int64
    m::Int64
end

function Base.show(io::IO, ℐ::ClimatologyInterpolator{I}) where {I}
    @unpack n, m = ℐ
    print(io, "ClimatologyInterpolator{$I}, $n x $m")
end

function ClimatologyInterpolator(𝒞::AbstractVector{Climatology},
                                 x::AbstractVector{<:Real},
                                 interpolator::Type=LinearInterpolator)
    @assert issorted(x) "x vector must be sorted in ascending order"
    @assert length(𝒞) > 1 "must have at least two Climatologies"
    n, m = size(𝒞[1])
    mask = 𝒞[1].mask
    A = 𝒞[1].A
    for 𝒸 ∈ 𝒞
        @assert size(𝒸) == (n,m) "Climatologies must all be the same size"
        @assert all(𝒸.mask .== mask) "Climatologies must all have identical masks (continental configurations)"
        @assert all(𝒸.A .≈ A)
    end
    #construct interpolators
    x = collect(Float64, x)
    r = Matrix{interpolator}(undef, n, m)
    T = Matrix{interpolator}(undef, n, m)
    for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            r[i,j] = interpolator(x, map(𝒸->𝒸.r[i,j], 𝒞), NoBoundaries())
            T[i,j] = interpolator(x, map(𝒸->𝒸.T[i,j], 𝒞), NoBoundaries())
        end
    end
    #construct unified interpolator
    ClimatologyInterpolator(x, mask, r, T, A, n, m)
end

function (ℐ::ClimatologyInterpolator)(x)
    @unpack mask, r, T, A, n, m = ℐ
    #interpolate runoff and temperature at all points
    rₓ = fill(NaN, (n, m))
    Tₓ = fill(NaN, (n, m))
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            rₓ[i,j] = r[i,j](x)
            Tₓ[i,j] = T[i,j](x)
        end
    end
    #construct a new Climatology
    Climatology(mask, rₓ, Tₓ, A, n, m)
end

function findequilibrium(ℐ::ClimatologyInterpolator,
                         y::Real, #zero value for finding 0 = 𝒻(C) - y
                         𝒻::F; #function operating on a Climatology
                         tol::Float64=1e-3,
                         maxiter::Int=1000) where {F}
    #use the secant method, which will nail linear interpolators quickly
    x₁, x₂ = ℐ.x[1], ℐ.x[end]
    δ₁, δ₂ = 𝒻(ℐ(x₁)) - y, 𝒻(ℐ(x₂)) - y
    x₃ = 0.0
    δ₃ = floatmax(float(y))
    n = 0
    while abs(δ₃) > tol || abs(δ₃)/abs(y) > tol
        #approximate zero
        x₃ = x₁ - δ₁*(x₂ - x₁)/(δ₂ - δ₁)
        δ₃ = 𝒻(ℐ(x₃)) - y
        #swap values
        x₁, x₂ = x₂, x₃
        δ₁, δ₂ = δ₂, δ₃
        #count
        n += 1
        n > maxiter && error("maximum iterations reached!")
    end
    return x₃
end

end