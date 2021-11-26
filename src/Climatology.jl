export Climatology

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*abs(Δϕ*(cos(θ₁) - cos(θ₂)))

function readgrid(fn, v)::Matrix{Float64}
    #read temperature file
    X = ncread(fn, v)
    #squash third dimension if necessary
    if ndims(X) > 2
        X = dropdims(X, dims=3)
    end
    size(X,1) > size(X,2) ? collect(Float64, transpose(X)) : collect(Float64, X)
end

struct Climatology
    mask::BitMatrix #land mask (1 for land, 0 otherwise)
    r::Matrix{Float64} #cell runoff [m/s]
    T::Matrix{Float64} #cell temperature [K]
    A::Matrix{Float64} #cell area [m^2]
    f::Matrix{Float64} #cell land fraction [-]
    n::Int64 #number of rows/latitudes
    m::Int64 #number of columns/longitudes
end

function Base.show(io::IO, 𝒸::Climatology)
    @unpack mask, r, T, A, f, n, m = 𝒸
    println(io, "$n x $m Climatology")
    Tmax = round(maximum(T[mask]), sigdigits=4)
    Tmin = round(minimum(T[mask]), sigdigits=4)
    println(io, "  temperature ∈ [$Tmin, $Tmax] K")
    rmax = round(maximum(r[mask]), sigdigits=4)
    rmin = round(minimum(r[mask]), sigdigits=4)
    println(io, "  runoff ∈ [$rmin, $rmax] m/s")
    println(io, "  land fraction = $(landfraction(𝒸))")
    N = sum(mask)
    print(io, "  $N/$(n*m) non-ocean cells")
end

Base.size(𝒸::Climatology) = (𝒸.n, 𝒸.m)

function Climatology(fnr::String, #runoff file name
                     vr::String,  #runoff variable name
                     nullr::Real, #runoff empty/fill value
                     convr::Real, #runoff conversion factor
                     fnT::String, #temperature file name
                     vT::String,  #temperature variable name
                     fnf::String, #land fraction file name
                     vf::String)  #land fraction variable name
    #read runoff grid
    r = readgrid(fnr, vr)
    #nullify null values
    @. r[r ≈ nullr] = NaN
    #replace negatives
    @. r[r < 0] = 0
    #apply conversion factor
    r .*= convr

    #read temperature grid
    T = readgrid(fnT, vT)

    #read cell land fractions
    f = readgrid(fnf, vf)
    #discard nonsense
    @. f[(f < 0) | (f > 1)] = 0

    #demand everything is the same size
    @assert size(r) == size(T) == size(f)
    n, m = size(r)

    #make a mask from the non-NaN runoff values
    mask = @. r |> isnan |> !

    #calculate cell areas, assuming equal spacing in lat & lon
    Δθ = π/n
    Δϕ = 2π/m
    A = zeros(Float64, n, m)
    for i ∈ 1:n
        A[i,:] .= cellarea(𝐑ₑ, Δϕ, (i-1)*Δθ, i*Δθ)
    end

    #construct
    Climatology(mask, r, T, A, f, n, m)
end

#--------------------------------------

export landfraction
export meanlandtemperature, meanlandrunoff

landfraction(𝒸::Climatology) = sum(𝒸.f .* 𝒸.A)/sum(𝒸.A)

function landmean(X::AbstractMatrix, 𝒸::Climatology)
    @unpack mask, A, f, T, n, m = 𝒸
    @assert size(X) == (n,m)
    s = 0.0
    a = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            s += A[i,j]*f[i,j]*X[i,j]
            a += A[i,j]*f[i,j]
        end
    end
    return s/a
end

meanlandtemperature(𝒸::Climatology) = landmean(𝒸.T, 𝒸)

meanlandrunoff(𝒸::Climatology) = landmean(𝒸.r, 𝒸)
