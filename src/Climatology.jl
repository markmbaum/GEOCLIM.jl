export Climatology

function readgrid(fn, v)::Matrix
    #read temperature file
    X = ncread(fn, v)
    #squash third dimension if necessary
    if ndims(X) > 2
        X = dropdims(X, dims=3)
    end
    size(X,1) > size(X,2) ? collect(transpose(X)) : X
end

struct Climatology{𝒯}
    mask::BitMatrix #land mask (1 for land, 0 otherwise)
    r::Matrix{𝒯} #cell runoff [m/s]
    T::Matrix{𝒯} #cell temperature [K]
    A::Matrix{𝒯} #cell area [m^2]
    f::Matrix{𝒯} #cell land fraction [-]
    lat::Vector{𝒯} #cell latitudes
    n::Int64 #number of rows/latitudes
    m::Int64 #number of columns/longitudes
end

function Base.show(io::IO, 𝒸::Climatology{𝒯}) where {𝒯}
    @unpack mask, r, T, A, f, n, m = 𝒸
    println(io, "$n x $m Climatology{$𝒯}")
    Tmax = round(maximum(T[mask]), sigdigits=4)
    Tmin = round(minimum(T[mask]), sigdigits=4)
    println(io, "  temperature ∈ [$Tmin, $Tmax] K")
    rmax = round(maximum(r[mask]), sigdigits=4)
    rmin = round(minimum(r[mask]), sigdigits=4)
    println(io, "  runoff ∈ [$rmin, $rmax] m/s")
    println(io, "  land fraction = $(landfraction(𝒸))")
    N = sum(mask)
    print(io, "  $N land cells, $(n*m) total cells")
end

Base.size(𝒸::Climatology) = (𝒸.n, 𝒸.m)

function Climatology(fnr::String, #runoff file name
                     vr::String,  #runoff variable name
                     nullr::Real, #runoff empty/fill value
                     convr::Real, #runoff conversion factor
                     fnT::String, #temperature file name
                     vT::String,  #temperature variable name
                     fnf::String, #land fraction file name
                     vf::String; #land fraction variable name
                     fnlat::String="", #empty will use runoff file
                     latname::String="lat")
    #read runoff grid
    r = readgrid(fnr, vr)
    #nullify null values
    @. r[r ≈ convert(eltype(r), nullr)] = NaN
    #replace negatives
    @. r[r < 0] = 0
    #apply conversion factor
    r .*= convr

    #get cell coordinates
    fnlat = isempty(fnlat) ? fnr : fnlat
    lat = collect(eltype(r), ncread(fnlat, latname))

    #read temperature grid
    T = readgrid(fnT, vT)

    #read cell land fractions
    f = readgrid(fnf, vf)
    #discard nonsense
    @. f[(f < 0) | (f > 1)] = 0

    #demand everything is the same size
    @assert size(r) == size(T) == size(f)
    n, m = size(r)

    #insist on the same types for grids
    𝒯 = promote_type(eltype.((r, T, f))...)
    r, T, f = convert.(Matrix{𝒯}, (r, T, f))

    #make a mask from the non-NaN runoff values
    mask = @. r |> isnan |> !

    #check that the mask agrees with land fractions
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            #should have land in the cell
            @assert 0 < f[i,j] <= 1 "Climatology mask indicates a land cell, but land fraction (f) is zero"
        else
            #should have ocean in the cell
            @assert f[i,j] == 0 "Climatology mask indicates an ocean cell, but land fraction (f) is nonzero"
        end
    end

    #cell areas
    A = zeros(𝒯, n, m)
    Δϕ = 2π/m
    θₘ = (π/180)*(lat[2:end] .+ lat[1:end-1])/2 .+ π/2
    a₁ = cellarea(𝐑ₑ, Δϕ, 0, θₘ[1])
    aₙ = cellarea(𝐑ₑ, Δϕ, θₘ[end], π)
    for j ∈ 1:m
        A[1,j] = a₁
        A[n,j] = aₙ
    end
    for i ∈ 2:n-1
        aᵢ = cellarea(𝐑ₑ, Δϕ, θₘ[i-1], θₘ[i])
        for j ∈ 1:m
            A[i,j] = aᵢ
        end
    end

    #construct
    Climatology{𝒯}(mask, r, T, A, f, lat, n, m)
end

#--------------------------------------

export landfraction, meanlandtemperature, meanlandrunoff

#already exported in main file
#landfraction(𝒸::Climatology) = sum(𝒸.f .* 𝒸.A)/sum(𝒸.A)

function landfraction(𝒸::Climatology{𝒯}; cut::Real=Inf) where {𝒯}
    @unpack A, f, lat, n, m = 𝒸
    @assert cut >= 0
    L = zero(𝒯)
    S = zero(𝒯)
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if -cut <= lat[i] <= cut
            L += A[i,j]*f[i,j]
            S += A[i,j]
        end
    end
    return L/S
end

function landmean(X::AbstractMatrix{𝒯}, 𝒸::Climatology{𝒯}, cut::Real=Inf) where {𝒯}
    @unpack mask, A, f, lat, T, n, m = 𝒸
    @assert size(X) == (n,m)
    @assert cut > 0
    s = zero(𝒯)
    a = zero(𝒯)
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j] & (-cut <= lat[i] <= cut)
            #land area of cell
            LA = A[i,j]*f[i,j]
            #contributions to averaging
            s += LA*X[i,j]
            a += LA
        end
    end
    return s/a
end

meanlandtemperature(𝒸::Climatology; cut::Real=Inf) = landmean(𝒸.T, 𝒸, cut)

meanlandrunoff(𝒸::Climatology; cut::Real=Inf) = landmean(𝒸.r, 𝒸, cut)
