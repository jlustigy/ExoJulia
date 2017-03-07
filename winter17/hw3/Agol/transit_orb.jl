function transit_orb(t,x,dt,nsub)
# Computes a transit lightcurve, normalized to unity, as a
# function of time, t, usually in HJD or BJD.
#
# Input parameters (x) are:
# x[1] = P  (units of day)
# x[2] = inc = inclination angle
# x[3] = p = R_p/R_* = radius of planet in units of radius of star
# x[4] = t0 = mid-point of transit
# x[5] = u1 = linear limb-darkening coefficient
# x[6] = u2 = quadratic limb-darkening coefficient
# x[7] = f0 = uneclipsed flux
# x[8] = a/R_* = semi-major axis divided by R_*
# x[9] = e = eccentricity
# x[10] = omega = longitude of pericentre
# Compute time of pericentre passage:
# The true anomaly at the time of transit:
f1=1.5*pi-x[10]*pi/180.0
ecc=x[9]
tp=(x[4]+x[1]*sqrt(1.0-ecc^2)/2.0/pi*(ecc*sin(f1)/(1.0+ecc*cos(f1))
    -2.0/sqrt(1.0-ecc^2)*atan2(sqrt(1.0-ecc^2)*tan(0.5*f1),1.0+ecc)))
flux_prior = -1.0
fluxoft = 0.0
nsub = 1
while abs(fluxoft-flux_prior) > 1e-10 
  flux_prior = fluxoft
  fluxoft = 0.0
  for j=1:nsub
    tsub = t-dt/2.0+dt*(j-0.5)/nsub
    m=2.0*pi/x[1]*(tsub-tp)
    f=kepler(m,ecc)
    radius=x[8]*(1.0-ecc^2)/(1.0+ecc*cos(f))
# Now compute sky separation:
    z0=radius*sqrt(1.0-(sin(x[2]*pi/180.0)*sin(x[10]*pi/180.0+f))^2)
    if (sin(x[10]*pi/180.0+f) < 0.0) && (z0 <= (1.0+x[3]))
      if x[3] < 1.0
        fluxoft += occultquad(z0,x[5],x[6],x[3])
      else
# We'll assume that the smaller object is not limb-darkened:
        fluxoft += occultquad(z0,  0.,  0.,x[3])
      end 
    else
      fluxoft += 1.0
    end
  end
  fluxoft /=nsub
  nsub *=2
end
# Take mean of flux over sub-time steps:
flux = fluxoft*x[7]
return flux
end
