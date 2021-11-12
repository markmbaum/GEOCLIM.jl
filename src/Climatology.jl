export Climatology

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
    @. r[r == fᵣ] = NaN
    #also replace negative values
    @. r[r < 0] = NaN
    #make a mask for future calculations
    mask = @. r |> isnan |> !
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