function rv_model(times, param)
# Computes radial velocity model
# There are three parameters per planet.
sz = size(param)
nplanet = round(Int8,(sz[1]-2)/5)
ntime = size(times)
ntime = ntime[1]
rvmod = zeros(ntime)
f=zeros(ntime)
for i=1:nplanet
  period = param[(i-1)*5+1]
  tp     = param[(i-1)*5+2]
  ecc    = param[(i-1)*5+3]
  m = (times-tp)*2pi/period
  for j=1:ntime
    eanom=kepler_solve(m[j],ecc)
    f[j]=2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(0.5*eanom))
  end
  rvmod += cos(f)*param[(i-1)*5+4]
  rvmod += sin(f)*param[(i-1)*5+5]
end
rvmod += param[(nplanet-1)*5+6]
rvmod += (times-mean(times))*param[(nplanet-1)*5+7]
end
