export ClimatologyInterpolator

struct ClimatologyInterpolator{I<:OneDimensionalInterpolator}
    x::Vector{Float64}
    mask::BitMatrix
    r::Matrix{I}
    T::Matrix{I}
    A::Matrix{Float64}
    n::Int64
    m::Int64
    L::Int64
end

function Base.show(io::IO, ℐ::ClimatologyInterpolator{I}) where {I}
    @unpack n, m = ℐ
    print(io, "ClimatologyInterpolator{$I}, $n x $m")
end

Base.size(ℐ::ClimatologyInterpolator) = (ℐ.n, ℐ.m, ℐ.L)

function ClimatologyInterpolator(𝒞::AbstractVector{Climatology},
                                 x::AbstractVector{<:Real},
                                 interpolator::Type=CubicSplineInterpolator,
                                 boundaries::Type=StrictBoundaries)
    @assert issorted(x) "x vector must be sorted in ascending order"
    @assert length(𝒞) > 1 "must have at least two Climatologies"
    n, m = size(𝒞[1])
    mask = 𝒞[1].mask
    A = 𝒞[1].A
    for 𝒸 ∈ 𝒞
        #demand identical size
        @assert size(𝒸) == (n,m) "Climatologies must all be the same size"
        #demand 
        @assert all(𝒸.A .≈ A)
        mask .*= 𝒸.mask
    end
    #construct interpolators
    x = collect(Float64, x)
    @multiassign r, T = Matrix{interpolator}(undef, n, m)
    for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            r[i,j] = interpolator(x, map(𝒸->𝒸.r[i,j], 𝒞), boundaries())
            T[i,j] = interpolator(x, map(𝒸->𝒸.T[i,j], 𝒞), boundaries())
        end
    end
    #construct unified interpolator
    ClimatologyInterpolator(x, mask, r, T, A, n, m, length(𝒞))
end

function (ℐ::ClimatologyInterpolator)(x)
    #pull out fields of struct
    @unpack mask, r, T, A, n, m = ℐ
    #assume NaN until unmasked
    @multiassign rₓ, Tₓ = fill(NaN, (n, m))
    #interpolate runoff and temperature at all points
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            rₓ[i,j] = r[i,j](x)
            Tₓ[i,j] = T[i,j](x)
        end
    end
    #construct a new Climatology
    Climatology(mask, rₓ, Tₓ, A, n, m)
end