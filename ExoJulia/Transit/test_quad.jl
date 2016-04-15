# This script runs the quadratic limb-darkening curve computation.
# Written by Eric Agol April 2016

using PyPlot
include("occultquad.jl")

function test_quad(p,b,v,u1,u2,t1,t2,nt)

# Input variables:
# p : radius ratio of planet to star
# b : impact parameter
# v : velocity of planet in units of stellar radii per unit time
# u1,2 : linear/quadratic limb-darkening parameters
# t1,2 : start & end times
# nt : number of times
# 
# Output:
# time : time array
# z    : separation of center of planet & star in units of stellar radii
# mu   : brightness of star versus time
# mu0  : brightness of star without limb-darkening
#

time = linspace(t1,t2,nt)

# Velocity is in units of solar radii per day
z = sqrt(b^2+(v.*time).^2)

mu = zeros(nt)
mu0 = zeros(nt)

for i=1:nt
  mu[i] =occultquad(z[i],u1,u2,p)
  mu0[i]=occultquad(z[i],0.,0.,p)
end


#plot(time,mu)
#plot(time,mu0)
return time,z,mu,mu0
end

p = 0.1
b = 0.0
v = 1.0
u1 = 0.5
u2 = 0.25
t1 = -1.5
t2 = 1.5
nt = 100000
@time time,z,mu,mu0=test_quad(p,b,v,u1,u2,t1,t2,nt) #(0.1,0.,1.,0.5,0.25,-1.5,1.5,100000)

plot(time,mu)
plot(time,mu0)
dmudz=zeros(nt-1)
tc = zeros(nt-1)
for i=1:nt-1
  dmudz[i] = (mu[i+1]-mu[i])/(time[i+1]-time[i])
  tc[i] = (time[i+1]+time[i])/2
end
plot(tc,dmudz+1.)
