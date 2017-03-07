# find transiting planet:
using PyPlot

include("bls.jl")
#     n    = number of data points
#     t    = array {t(i)}, containing the time values of the time series
#     x    = array {x(i)}, containing the data values of the time series
#     nf   = number of frequency points in which the spectrum is computed
#     fmin = minimum frequency (MUST be > 0)
#     df   = frequency step
#     nb   = number of bins in the folded time series at any test period
#     qmi  = minimum fractional transit length to be tested
#     qma  = maximum fractional transit length to be tested

data =readdlm("mysteryplanet3.txt")
t = vec(data[:,1])
# sampling time:
dt =  t[2]-t[1]
# total duration:
ttot = maximum(t)-minimum(t)
x = vec(data[:,2])

n = length(t)

# Minimum transit duration: 1/2 hour:
tmin = 0.5/24.
# Maximum transit duration: 10 hr:
tmax = 10./24.
# Minimum frequency is 1/2 duration of data:
fmin = 2.0/(maximum(t)-minimum(t))
# Maximum frequency of 1/2 day:
fmax = 1.0/0.5
# Number of frequencies to search (this should be fine enough such that a transit is not missed):
nf = 1000
df = (fmax-fmin)/nf
nb = ceil(Int64,30*24/.5)
qmi = tmin/ttot
qma = tmax*fmax

p,bper,bpow,depth,qtran,in1,in2,f0 = bls(n,t,x,nf,fmin,df,nb,qmi,qma) #,p,bper,bpow,depth,qtran,in1,in2)
clf()
plot(log10(1./f0),p)
# Find frequency where spectrum is maximum:
np = length(f0)
imax = 0
pow_best =  0.0
for i=1:np
 if p[i] > pow_best
   imax = i
   pow_best = p[i]
 end
end
# Make the Farey sequence:

plo = 1./(f0+df/2.)
phi = 1./(f0-df/2.)
phi_best = phi[imax]
plo = plo./phi_best
phi = phi./phi_best
spec = zeros(np)
# Loop over period & compute the Farey sequence:
for i=1:np
# For unequal bounds, do a bisection search:
# First, make sure the period ratio is between 0-1,
# and if not, run bisection on inverse periods:
  if (phi[i] < 1.0)
    x1=plo[i]
    x2=phi[i]
  else
    x1=1./phi[i] 
    x2=1./plo[i]
  end
  xint = floor(x1)
  x1=x1-xint 
  x2=x2-xint
# Set up the initial Farey sequence left & right
# bounds for numerator & denominator:
  l0n=0
  l0d=1
  r0n=1
  r0d=1
# Perform the first Farey addition a/b+c/d = (a+c)/(b+d):
  mn=l0n+r0n 
  md=l0d+r0d
  m=float(mn)/float(md)
# Now, keep doing Farey addition until a rational
# number falls in the desired range:
  while (m < x1) || (m > x2)
    if (m < x1)
# Move left bound to m:
      l0n=l0n+r0n
      l0d=l0d+r0d
    else
# Move right bound to m:
      r0n=l0n+r0n
      r0d=l0d+r0d
    end
    mn=l0n+r0n
    md=l0d+r0d
    m=float(mn)/float(md)
  end
  mn=mn+md*xint
  #println(plo[i]," ",phi[i]," ",mn," ",md)
# Compute the spectrum:
  spec[i]=1./sqrt(float(mn*md))
end
spec[imax]=1.0
read(STDIN,Char)
plot(log10(1./f0),spec/maximum(spec)*(maximum(p)-median(p))+median(p))
