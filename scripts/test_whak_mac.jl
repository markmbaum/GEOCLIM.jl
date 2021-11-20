using GEOCLIM
using UnPack

const T0 = 288.15 #[K]
const Ea = 48200.0 #[mole/m^3]
const Eab = 42300.0 #[mole/m^3]
const W = 2.5e12/(365*24*3600) #total weathering [moles/second]
const 𝐑 = 8.31446262

##

datadir = "/Users/minminfu/Dropbox/Mark-and-Minmin/GEOCLIM"
#datadir = "C:/Users/markm/Dropbox/Mark-and-Minmin/GEOCLIM"

##
# This one is the climatology. 
C0 = Climatology(
    joinpath("/Users/minminfu/GEOCLIM.jl/calibration/calibration_climatology_285ppm", "ROF_T31.nc"),
    "QRUNOFF",
    1e36,
    1e-3,
    joinpath("/Users/minminfu/GEOCLIM.jl/calibration/calibration_climatology_285ppm", "TS.nc"),
    "TS",
)

# WHAK (Goddéris formulation)
kG  = 0.043014577424652055
kGb = 0.007376358214930387
totalweathering(C0,kG,Ea,T0)+totalweathering(C0,kGb,Eab,T0)

# WHAK (Abbot formulation)
Te = 11.1
# Calibrate by simply tuning k to match modern weathering rate.
kG2 = W/totalweathering(C0,285e-6,1,Te,T0,285e-6)
totalweathering(C0,285e-6,kG2,Te,T0,285e-6)

# MAC 
Te = 11.1
#totalweathering_mac(C0,285,Te,T0,285)
#totalweathering_mac(C0,1000,Te,T0,285)



#------------------------------------------------------------------------------
# MAC, as implementated by Graham and Pierrehumbert 2020 
# following Maher and Chamberlin 2014
##
using Plots
Rvals = LinRange(0.001,10,10000)/86400/365
plot(Rvals*86400*365,
    weathering_mac.(Rvals,T0,1,0.1,Te,T0,280e-6)*31536000,
    ylims = (5e-4,5e-1),
    yticks = [1e-3,1e-2,1e-1],
    xticks = [1e-3,1e-2,1e-1,1e0],
    xaxis=:log, yaxis=:log, label="PCO2 = 0.1 bar",legend=:bottomright)
plot!(Rvals*86400*365,
    weathering_mac.(Rvals,T0,1,1,Te,T0,280e-6)*31536000,
    ylims = (5e-4,5e-1),
    yticks = [1e-2,1e-1],
    xticks = [1e-3,1e-2,1e-1,1e0],
    xaxis=:log, yaxis=:log, label="PCO2 = 1 bar",legend=:bottomright)
plot!(Rvals*86400*365,
    weathering_mac.(Rvals,T0,1,10,Te,T0,280e-6)*31536000,
    ylims = (5e-4,5e-1),
    yticks = [1e-2,1e-1],
    xticks = [1e-3,1e-2,1e-1,1e0],
    xaxis=:log, yaxis=:log, label="PCO2 = 10 bar",legend=:bottomright)
ylabel!("w (MAC) (mol/m2/yr")
xlabel!("q (m/yr)")

# MAC total weathering
totalweathering_mac(C0,285e-6,Te,T0,285e-6)