function erosion(r, S, kₑ, a, b)
    kₑ*(r^a)*(S^b)
end

function regprodopt(T, r, kᵣₚ, Eₐ)
    kᵣₚ*r*exp(-(Eₐ/𝐑)*(1/T - 1/𝐓₀))
end

function soilprod(h, h₀)
    h₀/h
end

function eqregthick(RPₒ, E, h₀)
    max(h₀*log(RPₒ/E), 0.0)
end

function dissolutionconstant(T, r, Eₐ, kd, kw)
    kd*(1 - exp(-kw*r))*exp((Eₐ/𝐑)*(1/𝐓₀ - 1/(T + 273.15)))
end 

function eqxpsurface(h, E, K, σ)
    exp(-K*((h/E)^(σ + 1))/(σ + 1))
end