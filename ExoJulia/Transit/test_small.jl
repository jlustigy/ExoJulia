using PyPlot
include("occultsmall.jl")
include("occultquad.jl")

# Computes lightcurve for small planet approximation.
nz=2001
mu=ones(nz)
p = 0.01
z=linspace(0,1+2p,nz)
c1 = c3 = 0
u1 = u2 = 0.25
c2 = u1 + 2u2
c4 = -u2
mu  = occultsmall(p,c1,c2,c3,c4,z)
plot(z,mu)
# Now, compare with occultquad.jl:
mu_exact = ones(nz)
for i=1:nz
  mu_exact[i]=occultquad(z[i],u1,u2,p)
end
plot(z,mu_exact)
