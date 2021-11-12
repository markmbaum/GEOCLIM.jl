using GEOCLIM
using NLsolve

## read GCM output into a Climatology object

dir = joinpath("calibration", "calibration_climatology_285ppm")
const C = Climatology(
    joinpath(dir, "ROF_T31.nc"),
    "QRUNOFF",
    1e36,
    1e-3,
    joinpath(dir, "TS.nc"),
    "TS"
)

## fixed parameters

const T0 = 288.15 #[K]
const Ea = 48200.0 #[mole/m^3]
const Eab = 42300.0 #[mole/m^3]
const W = 2.5e12/(365*24*3600) #total weathering [moles/second]

## the system of equations to minimize

function f!(F, x)::Nothing
    #unpack
    kG, kGb = x
    #weathering terms
    w₁ = totalweathering(C, kG, Ea, T0)
    w₂ = totalweathering(C, kGb, Eab, T0)
    #first equation is weathering balance
    F[1] = w₁ + w₂ - W
    #second argument is percentage of terms
    F[2] = w₂/(w₁ + w₂) - 0.14
    return nothing
end

## optimize

sol = nlsolve(f!, [4.35e-4, 6.46e-5], autodiff=:forward)
kG, kGb = sol.zero
println("calibration constants:")
println("  kG  = $kG")
println("  kGb = $kGb")
