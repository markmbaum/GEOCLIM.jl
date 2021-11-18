using GEOCLIM

##

datadir = "/Users/minminfu/Dropbox/Mark-and-Minmin/GEOCLIM"
#datadir = "C:/Users/markm/Dropbox/Mark-and-Minmin/GEOCLIM"

##

C1 = Climatology(
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-1.0_001_1000ppm_nooht", "ROF_T31.nc"),
    "QRUNOFF",
    1e36,
    1e-3,
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-1.0_001_1000ppm_nooht", "TS.nc"),
    "TS",
)

C2 = Climatology(
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-2.0_001_1000ppm_nooht", "ROF_T31.nc"),
    "QRUNOFF",
    1e36,
    1e-3,
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-2.0_001_1000ppm_nooht", "TS.nc"),
    "TS",
)

C3 = Climatology(
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-3.0_001_1000ppm_nooht", "ROF_T31.nc"),
    "QRUNOFF",
    1e36,
    1e-3,
    joinpath(datadir, "e.e12.E1850C4.T31_g37.1367_p-3.0_001_1000ppm_nooht", "TS.nc"),
    "TS",
)
