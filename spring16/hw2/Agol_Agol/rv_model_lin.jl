function rv_model_lin(times, param, rv, sigrv)
# Computes radial velocity model
# There are three parameters per planet; the others are linearized
sz = size(param)
nplanet = round(Int8,sz[1]/3)
#println(nplanet)
ntime = size(times)
ntime = ntime[1]
functions = spzeros(2*nplanet+1,ntime)
f=zeros(ntime)
for i=1:nplanet
  period = param[(i-1)*5+1]
  tp     = param[(i-1)*5+2]
  ecc    = param[(i-1)*5+3]
  m = (times-tp)*2pi/period
  for j=1:ntime
    eanom=ExoJulia.Orbit.kepler_solve(m[j],ecc)
    f[j]=2.0*atan(sqrt((1.0+ecc)/(1.0-ecc))*tan(0.5*eanom))
  end
  functions[1+(i-1)*2,:]=cos(f)
  functions[2+(i-1)*2,:]=sin(f)
end
functions[2*nplanet+1,:]=ones(ntime)
#functions[2*nplanet+2,:]=times-mean(times)
#println(functions)
#plot(times,vec(functions[1,:]))
#plot(times,vec(functions[2,:]))
#plot(times,vec(functions[3,:]))
coeff=regress(functions,rv,sigrv)
#println(coeff)
rvmod=zeros(ntime)
for j=1:2*nplanet+1
 rvmod += coeff[j]*vec(functions[j,:])
end
rvmod,coeff
end
