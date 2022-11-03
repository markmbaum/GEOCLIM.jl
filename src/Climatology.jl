export Climatology

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*abs(Δϕ*(cos(θ₁) - cos(θ₂)))

checktranspose(A::Matrix)::Matrix = (size(A,1) > size(A,2)) ? collect(transpose(A)) : A

function readgrid(fn, v)::Matrix
    #read temperature file
    X = ncread(fn, v)
    #squash third dimension if necessary
    d = ndims(X)
    if d == 3
        return checktranspose(dropdims(X, dims=3))
    elseif d == 2
        return checktranspose(X)
    end
    error("unusual number of dimensions ($d) found in file $fn, variable $v")
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

Base.size(𝒸::Climatology)::NTuple{2,Int64} = (𝒸.n, 𝒸.m)

function Base.size(𝒸::Climatology, dim::Int)::Int64
    @assert 1 <= dim <= 2 "Climatology has only two dimensions"
    size(𝒸)[dim]
end

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

    #insist on the same types for grids
    𝒯 = promote_type(eltype.((r, T, f))...)
    r, T, f = convert.(Matrix{𝒯}, (r, T, f))

    #demand everything is the same size
    @assert size(r) == size(T) == size(f) "Climatology arrays must have the same size/shape"
    n, m = size(r)
    @assert length(lat) == n "latitude vector length ($(length(lat))) must match number of Climatology rows ($n)"

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

export meanlandtemperature
export meanlandrunoff, totallandrunoff

checkcut(cut) = @assert cut >= 0 "latitude cutoff must be positive"

checksize(𝒸, X) = @assert size(X) == size(𝒸) "size mismatch between array and Climatology"

#already exported
function landfraction(𝒸::Climatology{𝒯}, cut::Real=Inf) where {𝒯}
    @unpack mask, A, f, lat, n, m = 𝒸
    checkcut(cut)
    @multiassign num, den = zero(𝒯)
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j] & (-cut <= lat[i] <= cut)
            num += A[i,j]*f[i,j]
            den += A[i,j]
        end
    end
    return iszero(den) ? zero(𝒯) : num/den
end

function landmean(X::AbstractMatrix{𝒯}, 𝒸::Climatology{𝒯}, cut::Real=Inf) where {𝒯}
    @unpack mask, A, f, lat, n, m = 𝒸
    checksize(𝒸, X)
    checkcut(cut)
    @multiassign num, den = zero(𝒯)
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j] & (-cut <= lat[i] <= cut)
            #land area of cell
            LA = A[i,j]*f[i,j]
            #contributions to averaging
            num += LA*X[i,j]
            den += LA
        end
    end
    return iszero(den) ? zero(𝒯) : num/den
end

function landsum(X::AbstractMatrix{𝒯}, 𝒸::Climatology{𝒯}, cut::Real=Inf) where {𝒯}
    @unpack mask, A, f, lat, n, m = 𝒸
    checksize(𝒸, X)
    checkcut(cut)
    s = zero(𝒯)
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j] & (-cut <= lat[i] <= cut)
            #land area of cell
            LA = A[i,j]*f[i,j]
            #contribution to sum
            s += LA*X[i,j]
        end
    end
    return s
end

meanlandtemperature(𝒸::Climatology, cut::Real=Inf) = landmean(𝒸.T, 𝒸, cut)

meanlandrunoff(𝒸::Climatology, cut::Real=Inf) = landmean(𝒸.r, 𝒸, cut)

totallandrunoff(𝒸::Climatology, cut::Real=Inf) = landsum(𝒸.r, 𝒸, cut)

#already exported
meanlandlatitude(𝒸::Climatology) = landmean(repeat(𝒸.lat, 1, 𝒸.m), 𝒸)

#already exported
meanabslandlatitude(𝒸::Climatology) = landmean(repeat(abs.(𝒸.lat), 1, 𝒸.m), 𝒸)

#--------------------------------------
export meanoceandistance

function sph2cart(θ::T, ϕ::T) where {T}
    sₜ, cₜ = sincos(θ)
    sₚ, cₚ = sincos(ϕ)
    return SVector{3,T}(sₜ*cₚ, sₜ*sₚ, cₜ)
end

function arclength(c₁::SVector{3,T}, c₂::SVector{3,T}) where {T}
    d = c₁ ⋅ c₂
    if d > one(T)
        return zero(T)
    elseif d < -one(T)
        return convert(T,π)
    end
    acos(d)
end

function meanoceandistance(𝒸::Climatology{𝒯}, cut::Real=Inf, R::Real=𝐑ₑ) where {𝒯}
    @unpack mask, lat, n, m = 𝒸
    checkcut(cut)
    @assert any(mask .== 0) "no ocean cells, can't compute distances to ocean"
    #we assume longitude values cells are evenly spaced, as they ought to be
    Δϕ = 2π/m
    ϕ = collect(𝒯, LinRange(Δϕ/2, 2π - Δϕ/2, m))
    #convert the latitude values to radians
    θ = collect(𝒯, lat*(π/180) .+ π/2)
    #create a grid of cartesian coordinates for each cell
    C = Matrix{SVector{3,𝒯}}(undef, n, m)
    for i ∈ 1:n, j ∈ 1:m
        C[i,j] = sph2cart(θ[i], ϕ[j])
    end
    #mean arclength from land cells to ocean cells
    ℒ = zero(𝒯)
    count::Int64 = 0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        #check if its a land cell withing the desired latitude band
        if mask[i,j] & (-cut <= lat[i] <= cut)
            #find the minimum distance to the ocean for cell i,j
            𝓁 = Inf
            cᵢⱼ = C[i,j]
            for k ∈ 1:n, l ∈ 1:m
                #cell k,l must be ocean
                if !mask[k,l] & ((i != k) | (j != l)) #don't check a cell against itself
                    𝓁ₖₗ = arclength(cᵢⱼ, C[k,l])
                    if 𝓁ₖₗ < 𝓁
                        𝓁 = 𝓁ₖₗ
                    end
                end
            end
            count += 1
            ℒ += 𝓁
        end
    end
    #normalize by cell count and convert from radians to meters
    R*ℒ/count
end

#--------------------------------------
export perimeter

function perimeter(mask::BitMatrix, lat::Vector{𝒯}, R::Real=𝐑ₑ) where {𝒯}
    @assert size(mask,1) == length(lat)
    n, m = size(mask)
    B = CircularArray(mask)
    Δϕ = convert(𝒯, 2π/m)
    Δθ = convert(𝒯, π/n)
    c = cos.(LinRange(-π/2, π/2, n+1))
    p = zero(𝒯)
    for i ∈ 1:n, j ∈ 1:m
        if B[i,j]
            #look left and right
            for k ∈ (j-1,j+1)
                if !B[i,k]
                    p += Δθ*R
                end
            end
            #look "up" 
            if !B[i-1,j]
                p += Δϕ*R*c[i]
            end
            #look "down" 
            if !B[i+1,j]
                p += Δϕ*R*c[i+1]
            end
        end
    end
    return p
end

perimeter(𝒸::Climatology, R::Real=𝐑ₑ) = perimeter(𝒸.mask, 𝒸.lat, R)
