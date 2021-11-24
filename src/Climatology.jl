export Climatology

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*abs(Δϕ)*abs(cos(θ₁) - cos(θ₂))

function checktranspose(X, n::Int, m::Int)::Matrix{Float64}
    if size(X) == (n,m)
        return collect(Float64, X)
    elseif size(X) == (m,n)
        return collect(Float64, transpose(X))
    else
        error("Incompatible matrix dimensions, looking for ($n,$m) or ($m,$n)")
    end
end

function readgrid(fn, v, n::Int, m::Int)::Matrix{Float64}
    #read temperature file
    X = ncread(fn, v)
    #squash third dimension if necessary
    if ndims(X) > 2
        X = dropdims(X, dims=3)
    end
    #check size and if transpose is needed
    X = checktranspose(X, n, m)
    return X
end

struct Climatology
    mask::BitMatrix #land mask (1 for land, 0 otherwise)
    r::Matrix{Float64} #cell runoff [m/s]
    T::Matrix{Float64} #cell temperature [K]
    A::Matrix{Float64} #cell area [m^2]
    f::Matrix{Float64} #cell land fraction
    n::Int64 #number of rows
    m::Int64 #number of columns
end

function Base.show(io::IO, 𝒞::Climatology)
    @unpack mask, r, T, A, f, n, m = 𝒞
    println(io, "$n x $m Climatology")
    Tmax = round(maximum(T[mask]), sigdigits=4)
    Tmin = round(minimum(T[mask]), sigdigits=4)
    println(io, "  temperature ∈ [$Tmin, $Tmax] K")
    rmax = round(maximum(r[mask]), sigdigits=4)
    rmin = round(minimum(r[mask]), sigdigits=4)
    println(io, "  runoff ∈ [$rmin, $rmax] m/s")
    F = round(sum(A[mask] .* f[mask])/sum(A), sigdigits=4)
    println(io, "  land fraction = $F")
    N = sum(mask)
    print(io, "  $N/$(n*m) non-ocean cells")
end

Base.size(𝒞::Climatology) = (𝒞.n, 𝒞.m)

function Climatology(fnr::String,   #runoff file name
                     vr::String,    #runoff variable name
                     nullr::Real,   #runoff empty/fill value
                     convr::Real,   #runoff conversion factor
                     fnT::String,   #temperature file name
                     vT::String,    #temperature variable name
                     fnf::String,   #land fraction file name
                     vf::String)    #land fraction variable name

    #first read the runoff file
    r = ncread(fnr, vr)
    #squash third dimension
    r = dropdims(r, dims=3)
    #replace missing data with NaNs
    nullr = convert(eltype(r), nullr)
    @. r[r == nullr] = NaN
    #also replace negative values
    @. r[r < 0] = NaN
    #make a mask for future calculations
    mask = @. r |> isnan |> !
    #convert units
    r .*= convr
    #convert types
    r = collect(Float64, r)
    #store size
    n, m = size(r)

    #read temperature grid
    T = readgrid(fnT, vT, n, m)
    #read cell land fraction grid
    f = readgrid(fnf, vf, n, m)
    #null land fractions
    @. f[(f < 0) | (f > 1)] = 0

    #calculate cell areas, assuming equal spacing in lat & lon
    Δθ = π/n
    Δϕ = 2π/m
    A = zeros(Float64, n, m)
    for i ∈ 1:n
        θ₁ = (i-1)*Δθ
        θ₂ = θ₁ + Δθ
        A[i,:] .= cellarea(𝐑ₑ, Δϕ, θ₁, θ₂)
    end

    #construct
    Climatology(mask, r, T, A, f, n, m)
end