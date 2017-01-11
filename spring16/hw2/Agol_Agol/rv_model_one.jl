function rv_model_one(times, param)
# Computes radial velocity model for one planet
ntime = size(times)
ntime = ntime[1]
f=zeros(ntime)
period = param[1]
tp     = param[2]
ecc    = param[3]
m = (times-tp)*2.0*pi/period
for j=1:ntime
  eanom=ExoJulia.Orbit.kepler_solve(m[j],ecc)
  f[j]=2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(0.5*eanom))
end
# Add in zero-point:
rvmod = param[4]*cos(f)+param[5]*sin(f)+param[6]
end
