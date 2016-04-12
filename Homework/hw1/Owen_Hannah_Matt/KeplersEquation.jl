#=
TrueAnomaly()

Solve Kepler's equation for the true anomaly using Newton's method.

Input
e   The eccentricity of the orbit [0<e<1]
m   The mean anomaly in radians [0<m<2pi]
t   The tolerance desired for this calculation in decimal format [i.e. 0.00001]

Return
f   The true anomaly [0<f<2pi] on success. -1 on failure
=#
function TrueAnomaly(e, m, t)
  #tan(f/2)=sqrt((1+e)/(1-e))*tan(ea/2)
  ea = EccentricAnomaly(e, m, t)
  if ea == -1
    return ea
  end
  f = atan2(sqrt(1-e^2)*sin(ea),cos(ea)-e)
end

#=
EccentricAnomaly()

Solve Kepler's equation for the true anomaly using Newton's method.

Input
e   The eccentricity of the orbit [0<e<1]
m   The mean anomaly in radians [0<m<2pi]
t   The tolerance desired for this calculation in decimal format [i.e. 0.00001]

Return
f   The eccentric anomaly [0<ea<2pi] on success. -1 on failure
=#
function EccentricAnomaly(e, m, t)
  ea = m+0.85*e*sign(sin(m)) #initial guess of eccentric anomaly
  d = ea - e*sin(ea) - m

  num_iter = 0
  iter_limit = 30
  while abs(d) > t && num_iter < iter_limit
    num_iter+=1
    deltaEa = d/(1 - e*cos(ea))
    ea = ea - deltaEa
    d = ea - e*sin(ea) - m
  end

  f = -1
  if num_iter < iter_limit
    #success
    return ea
  end
  return f
end

#@stest [[EccentricAnomaly(ecc,M,1e-12) for M in linspace(0,20pi,1000)] for ecc in linspace(0,0.999,100)]
