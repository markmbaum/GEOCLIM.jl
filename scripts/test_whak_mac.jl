using GEOCLIM
using UnPack

const T0 = 288.15 #[K]
const Ea = 48200.0 #[mole/m^3]
const Eab = 42300.0 #[mole/m^3]
const W = 2.5e12/(365*24*3600) #total weathering [moles/second]
const 𝐑 = 8.31446262

kG  = 0.043014577424652055
kGb = 0.007376358214930387
# WHAK (Goddéris formulation)
totalweathering(C0,kG,Ea,T0)+totalweathering(C0,kGb,Eab,T0)

# WHAK (Abbot formulation)
Te = T0^2*𝐑/Ea
# Calibrate by simply tuning k to match modern weathering rate.
kG2 = W/totalweathering(C0,285,1,Ea,T0,285,0.2)
totalweathering(C0,285,kG2,Ea,T0,285,0.2)


#------------------------------------------------------------------------------
# MAC, as implementated by Graham and Pierrehumbert 2020 
# following Maher and Chamberlin 2014

#const T0 = 288.15 #[K]
##
const n = 0.316 # Thermodynamic pCO2 dependence [-]
const Λ = 1.4e-3 # Thermodynamic coefficient for Ceq [-]
const L = 1 # Flow path length [m]
const ϕ = 0.1 # Porosity [-]
const ρ = 12728 # Mineral mass to fluid volume ratio [kg m⁻³]
const k₀ = 8.7e-6 # Reference rate constant [mol m⁻² yr⁻¹]
const 𝐀 = 100 # Specific surface area [m²kg⁻¹]
const X = 0.36 # Reactive mineral conc. in fresh rock [-]
const tₛ = 1e5 # Soil age [yr]
const m = 0.27 # Mineral molar mass [kg/mol]
const μ = ℯ^2 # Scaling constant [-]
const α = L*ϕ*ρ*𝐀*X*μ # Defined for convenience [-]
##
T0 = 288.15 #[K]
Te = 11.1   #[K]
β  = 0.2
function Ceq(pCO2)
    return Λ*(pCO2*1e-6)^n*1000 #conversion from mol/liter to mol/m3
end

function weathering_mac(r, T, A, pCO2, Tₑ, T₀, pCO2₀, β) 
    #k*A*α*((k₀*exp((Eₐ/𝐑)*(T - T₀)/T₀^2)*(pCO2/pCO2₀)^β)^-1 + m*𝐀*tₛ + α/(r*Ceq(pCO2)))^-1
    A*α*((k₀*exp((T - T₀)/Tₑ)*(pCO2/pCO2₀)^β)^-1 + m*𝐀*tₛ + α/(r*Ceq(pCO2)))^-1
    #k*A*α*r*((k₀*exp((Eₐ/𝐑)*(T - T₀)/T₀^2)*(pCO2/pCO2₀)^β)^-1 + m*𝐀*tₛ + α/(365*r*Ceq(pCO2)))^-1
end

function weathering_mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀, β) 
    weathering_mac.(𝒸.r, 𝒸.T, 𝒸.A, pCO2, Tₑ, T₀, pCO2₀, β)
end

##
using Plots
Rvals = LinRange(0.001,10,1000)
y = [1e-3,1e-2,1e-1,1e0]
x = [1e-3,1e-2,1e-1,1e0]
plot(Rvals,weathering_mac.(Rvals,T0,1,0.1,Te,T0,280e-6,β), ytick = y, xtick = x,
xaxis=:log, yaxis=:log, label="PCO2 = 0.1bar")
plot!(Rvals,weathering_mac.(Rvals,T0,1,1,Te,T0,280e-6,β), ytick = y, xtick = x,
xaxis=:log, yaxis=:log, label="PCO2 = 1bar")#,minoryticks=10)
plot!(Rvals,weathering_mac.(Rvals,T0,1,10,Te,T0,280e-6,β), ytick = y, xtick = x,
xaxis=:log, yaxis=:log, label="PCO2 = 10bar",legend=:bottomright)#,minoryticks=10)
ylabel!("w (MAC) (mol/m2/yr)")
xlabel!("q (m/yr)")

function totalweathering_mac(𝒸::Climatology, pCO2, Tₑ, T₀, pCO2₀, β)
    @unpack mask, r, T, A, n, m = 𝒸
    ΣW = 0.0
    @inbounds for i ∈ 1:n, j ∈ 1:m
        if mask[i,j]
            ΣW += weathering_mac(r[i,j], T[i,j], A[i,j], pCO2, Tₑ, T₀, pCO2₀, β)
        end
    end
    return ΣW
end



