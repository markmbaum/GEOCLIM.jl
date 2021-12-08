export ClimatologyInterpolator

struct ClimatologyInterpolator{I<:OneDimensionalInterpolator,𝒯}
    x::Vector{𝒯}
    mask::BitMatrix
    𝒻r::Matrix{I}
    𝒻T::Matrix{I}
    A::Matrix{𝒯}
    f::Matrix{𝒯}
    n::Int64
    m::Int64
    L::Int64
end

function Base.show(io::IO, ℐ::ClimatologyInterpolator{I,𝒯}) where {I,𝒯}
    @unpack n, m = ℐ
    print(io, "ClimatologyInterpolator{$I,$𝒯}, $n x $m")
end

Base.size(ℐ::ClimatologyInterpolator) = (ℐ.n, ℐ.m, ℐ.L)

function ClimatologyInterpolator(𝒞::AbstractVector{Climatology{𝒯}},
                                 x::AbstractVector{<:Real},
                                 interpolator::Type=CubicSplineInterpolator,
                                 boundaries::Type=StrictBoundaries) where {𝒯}
    @assert issorted(x) "x vector must be sorted in ascending order"
    @assert length(𝒞) > 1 "must have at least two Climatologies"
    n, m = size(𝒞[1])
    mask = 𝒞[1].mask
    A = 𝒞[1].A
    f = 𝒞[1].f
    for 𝒸 ∈ 𝒞
        #demand identical size
        @assert size(𝒸) == (n,m) "Climatologies must all be the same size"
        #demand identical cell areas
        @assert all(𝒸.A .≈ A) "cell area grids must be identical"
        #demand identical cell land fractions
        @assert all(𝒸.f .≈ f) "cell land fraction grids must be identical"
        #union of all masked cells
        mask .*= 𝒸.mask
    end
    #construct interpolators
    x = collect(𝒯, x)
    @multiassign 𝒻r, 𝒻T = Matrix{interpolator}(undef, n, m)
    for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            𝒻r[i,j] = interpolator(x, map(𝒸->𝒸.r[i,j], 𝒞), boundaries())
            𝒻T[i,j] = interpolator(x, map(𝒸->𝒸.T[i,j], 𝒞), boundaries())
        end
    end
    println(x)
    #construct unified interpolator
    ClimatologyInterpolator(x, mask, 𝒻r, 𝒻T, A, f, n, m, length(𝒞))
end

function (ℐ::ClimatologyInterpolator{I,𝒯})(x) where {I,𝒯}
    #pull out fields of struct
    @unpack mask, 𝒻r, 𝒻T, A, f, n, m = ℐ
    #assume NaN until unmasked
    @multiassign r, T = fill(convert(𝒯, NaN), (n, m))
    #interpolate runoff and temperature at all points
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            #interpolate runoff
            r[i,j] = 𝒻r[i,j](x)
            #interpolate temperature
            T[i,j] = 𝒻T[i,j](x)
        end
    end
    #construct a new Climatology
    Climatology(mask, r, T, A, f, n, m)
end
