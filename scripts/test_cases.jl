using GEOCLIM

##

#datadir = "/Users/minminfu/Dropbox/mark-minmin-paleoweathering/sample-climatologies"
datadir = "C:/Users/markm/Dropbox/mark-minmin-paleoweathering/sample-climatologies"

##

#broken up
casedir = joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-1.0_001_1000ppm_nooht")
C1 = Climatology(
    joinpath(casedir, "ROF_T31.nc"),
    "QRUNOFF",
    1f36,
    1e-3,
    joinpath(casedir, "TS.nc"),
    "TS",
    joinpath(casedir, "ROF_T31_landfrac.nc"),
    "landfrac"
)

#sorta broken up
casedir = joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-2.0_001_1000ppm_nooht")
C2 = Climatology(
    joinpath(casedir, "ROF_T31.nc"),
    "QRUNOFF",
    1f36, 
    1e-3,
    joinpath(casedir, "TS.nc"),
    "TS",
    joinpath(casedir, "ROF_T31_landfrac.nc"),
    "landfrac"
)

#not broken up
casedir = joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-3.0_001_1000ppm_nooht")
C3 = Climatology(
    joinpath(casedir, "ROF_T31.nc"),
    "QRUNOFF",
    1f36,
    1e-3,
    joinpath(casedir, "TS.nc"),
    "TS",
    joinpath(casedir, "ROF_T31_landfrac.nc"),
    "landfrac"
)