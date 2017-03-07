using PyPlot
include("transit_orb.jl")
include("kepler.jl")
include("occultquad.jl")
# simulate mysteryplanet02.txt

# Do 10-minute sampling for 30 days:
texp = 10./24./60. # convert to days
ntime = 30*24*6
t = linspace(0,30,ntime)  # uniformly spaced data
nsub = 10
# Set parameters of model:
#   P[d]  i     k       t_0     a1 a2  f0 a/R_* e \omega \varpi
x = [3.0,90.0,0.0075,rand()*30.,0.,0.,1.0,10.0,0.,0.]

model = zeros(ntime)
for i=1:ntime
  model[i] = transit_orb(t[i],x,texp,nsub)
end
# random noise
# randn = random normal distribution
noise = randn(ntime)*1e-4
#noise = randn(ntime).*0.

simdata = model + noise
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

plot(t,model,color="blue")
scatter(t,simdata,color="red")

data = zeros(ntime,2)
data[:,1]=t
data[:,2]=simdata
writedlm("mysteryplanet3.txt",data)
