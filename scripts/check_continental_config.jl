using NetCDF

##

ncdir = "C:/Users/markm/Dropbox/Mark-and-Minmin/continental-configuration/ensemble-configurations"

##

#area of a grid box rectangular in latitude and longitude
# colatitude θ ∈ [0,π]
# longitude ϕ ∈ [0,2π]
cellarea(r, Δϕ, θ₁, θ₂) = (r^2)*abs(Δϕ)*abs(cos(θ₁) - cos(θ₂))

##

fn = joinpath(ncdir, "topo_ens_n12_p-2_L30_0001.nc")
lat = ncread(fn, "lat")
lon = ncread(fn, "lon")
topo = transpose(ncread(fn, "topo"));

@assert diff(lon)[1] == diff(lat)[1]
#grid spacing in radians
Δ = diff(lon)[1]*(π/180)
#total area
A = 0.0
#land area
L = 0.0
for i ∈ 1:size(topo,1)
    #cell colatitude edges
    θ₁ = (lat[i] + 90)*(π/180) - Δ/2
    θ₂ = (lat[i] + 90)*(π/180) + Δ/2
    for j ∈ 1:size(topo,2)
        #cell area
        a = cellarea(1.0, Δ, θ₁, θ₂)
        #contribute to total area sum
        A += a
        #contribute to land area sum
        if topo[i,j] > 0
            L += a
        end
    end
end
println("land fraction = $(L/A)")