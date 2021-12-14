using GEOCLIM
using Test

##

#test climatology directories
casedirs = [
    joinpath("climatologies","e.e12.E1850C4.T31_g37.1367_p-2.0_001_10ppm_nooht"),
    joinpath("climatologies","e.e12.E1850C4.T31_g37.1367_p-2.0_001_100ppm_nooht"),
    joinpath("climatologies","e.e12.E1850C4.T31_g37.1367_p-2.0_001_1000ppm_nooht")
]

#read climatologies
C = map(casedirs) do casedir
    Climatology(
        joinpath(casedir, "ROF_T31.nc"),
        "QRUNOFF",
        1e36,
        1e-3,
        joinpath(casedir, "TS.nc"),
        "TS",
        joinpath(casedir, "ROF_T31_landfrac.nc"),
        "landfrac"
    )
end

#print one of them
println(C[1])

#co2 concentrations for each climatology
co2 = [1e1, 1e2, 1e3]

#test topography directories
topos = [
    joinpath("topographies","topo_ens_n12_p-1_L30_0001.nc"),
    joinpath("topographies","topo_ens_n12_p-1_L30_0002.nc"),
    joinpath("topographies","topo_ens_n12_p-1_L30_0003.nc")
]

##

@test all(0.29 .< landfraction.(topos) .< 0.31)

@testset "Land Quantities" begin
    #land fractions should be close to 30 %
    @test all(0.29 .< landfraction.(C) .< 0.31)
    #mean land temperature/runoff should be reasonable
    @test all(250 .< meanlandtemperature.(C) .< 300)
    @test all(1e-9 .< meanlandrunoff.(C) .< 2e-8)
end

#basic check on total weathering function ranges
@testset "Total Weathering" begin
    @test all(1e4 .< godderis.(C, 0.043, 48200, 288.15) .< 2e5)
    @test all(1e4 .< whak.(C, co2*1e-6, 0.043, 11.1, 288.15, 285e-6) .< 2e5)
    @test all(1e4 .< mac.(C, co2*1e-6, 11.1, 288.15, 285e-6) .< 2e5)
end

#test some interpolation
I = ClimatologyInterpolator(C, log10.(co2))

@testset "Equlibrium CO2" begin
    f₁(_,c) = godderis(c, 0.043, 48200, 288.15)
    @test 1.5 < findequilibrium(I, f₁, 5e4) < 2.5
    f₂(x,c) = whak(c, exp10(x)*1e-6, 0.043, 11.1, 288.15, 285e-6)
    @test 1.5 < findequilibrium(I, f₂, 5e4) < 2.5
    f₃(x,c) = mac(c, exp10(x)*1e-6, 11.1, 288.15, 285e-6)
    @test 1.5 < findequilibrium(I, f₃, 5e4) < 2.5
end
