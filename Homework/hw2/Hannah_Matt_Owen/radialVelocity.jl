# Add ExoJulia/ to path
push!(LOAD_PATH, "../../../ExoJulia")
using ExoJulia

#=
TrueAnomaly()

Solve Kepler's equation for the true anomaly using Newton's method.

Input
e   The eccentricity of the orbit [0<e<1]
m   The mean anomaly in radians [0<m<2pi]

Return
f   The true anomaly [0<f<2pi] on success. -1 on failure
=#
function true_anomaly(e::Float64, m::Float64)
  #tan(f/2)=sqrt((1+e)/(1-e))*tan(ea/2)
  ea = ExoJulia.Orbit.kepler_solve(m,e)
  f = atan2(sqrt(1-e^2)*sin(ea),cos(ea)-e)
end

function radial_velocity(p::Float64, mp::Float64, ms::Float64, inc::Float64,
                         f::Float64, w::Float64, gamma::Float64)
  #calculate the semi-amplitude k
  k = (2pi * 6.67259E-8 / (p * (1 - e^2)^(1.5)) * (mp^5*sin(inc))^3 / (mp + ms)^2)^(1/3)

  vRad = k*(cos(w+f)*ecc*cos(w))+gamma

end
