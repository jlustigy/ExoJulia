# Compute magnification for a finite source planet + moon + star
# Positions of bodies (in units of Einstein radius of star):

import PyPlot
using CGS

function planet_moon_star(uplanet,u1moon,u2moon,mu1,mu2)
# Computes ray tracing of planet + moon + star.
# Inputs:
#   uplanet: location of planet in image plane (along x-axis)
#     in units of Einstein angle of star
#   u1moon, u2moon: location of moon relative to planet
#     in units of Einstein angle *of planet*
#   mu1: mass ratio of planet to moon
#   mu2: mass ratio of moon to planet
#
# Thanks to Brendan Brewer's raytrace.jl:
# https://github.com/eggplantbren/Microlensing/blob/master/raytrace.jl
#
# Set up coordinates of the 3 bodies (star is at origin)
# in the horizontal direction in the lens plane:
x1body = [0,uplanet,uplanet+u1moon*sqrt(mu1)]
# An example is an Earth around an M dwarf:
# mu = MEARTH/0.1/MSUN
# Set up masses of the bodies:
mbody = [1,mu1,mu2*mu1]
# Set up positions of bodies in the vertical axis in the lens plane:
x2body = [0,0,u2moon*sqrt(mu1)]
println("Mass ratio: ",mu1)

# Set up grid of rays in the lens plane:
x1min_norm = -5; x1max_norm = 5
# Choose the rays around the planet, and set the range
# in units of the Einstein radius of the planet:
# horizontal range:
x1min = x1min_norm*sqrt(mu1)+uplanet; x1max = x1max_norm*sqrt(mu1)+uplanet; nx1 = 10000
# vertical range:
x2min_norm = -15; x2max_norm = 15
x2min = x2min_norm*sqrt(mu1); x2max = x2max_norm*sqrt(mu1); nx2 = 10000

# Now 'shoot' these to the source plane:

# We're going to collect rays over a grid in the source planet (y-coordinates):
ny1 = 1000; ny2 = 1000
# Create maps to hold results of ray shooting with & without planet/moon:
map1 = zeros(Int64,(ny1,ny2))
map2 = zeros(Int64,(ny1,ny2))
# Compute the position of the planet as microlensed by star in source plane:
yp= uplanet-1/uplanet
# Set up grid in source plane:
y1min_norm= -5; y1max_norm= 5;
y2min_norm= -2; y2max_norm= 2;
# Compute range of grid in source plane in units of Einstein angle
# of the planet (sqrt(mu1)):
y1min = sqrt(mu1)*y1min_norm+yp; y1max = sqrt(mu1)*y1max_norm+yp
y2min = sqrt(mu1)*y2min_norm; y2max = sqrt(mu1)*y2max_norm
# Spacing of grid in source plane:
dy1 = (y1max-y1min)/ny1
dy2 = (y2max-y2min)/ny2
# Spacing of grid in lens plane:
dx1 = (x1max-x1min)/nx1
dx2 = (x2max-x2min)/nx2

# Keep track of the range of rays to see if we need to extend the
# grid:
x1range1 = Inf
x1range2= -Inf
x2range1 = Inf
x2range2= -Inf
# Okay, now loop over the grid in the lens plane:
for i1 in 1:nx1
  x1ray = x1min + (i1-0.5)*dx1
  for i2 in 1:nx2
    x2ray = x2min + (i2-0.5)*dx2
# x1ray, x2ray are positions of ray in the lens plane.  Compute
# their angular offset from the 3 bodies:
    r12 = (x1ray-x1body[1])^2+(x2ray-x2body[1])^2
    r22 = (x1ray-x1body[2])^2+(x2ray-x2body[2])^2
    r32 = (x1ray-x1body[3])^2+(x2ray-x2body[3])^2
# Compute the location of ray in source plane using the (vector) lens equation:
    y1ray =  x1ray -  mbody[1] * (x1ray-x1body[1])/r12 -   mbody[2] * (x1ray-x1body[2])/r22 -   mbody[3] * (x1ray-x1body[3])/r32
    y2ray =  x2ray -  mbody[1] * (x2ray-x2body[1])/r12 -   mbody[2] * (x2ray-x2body[2])/r22 -   mbody[3] * (x2ray-x2body[3])/r32
# Find the position in the source-plane grid:
    i1ray = convert(Int64,floor((y1ray - y1min)/dy1))+1
    i2ray = convert(Int64,floor((y2ray - y2min)/dy2))+1
# See if the ray lands on the grid in the source plane:
    if (1 <= i1ray <= ny1) && (1 <= i2ray <= ny2) 
# Track the range of rays which land on grid:
      if x1ray < x1range1 
        x1range1 = x1ray
      end
      if x1ray > x1range2
        x1range2 = x1ray
      end
      if x2ray < x2range1
        x2range1 = x2ray
      end
      if x2ray > x2range2
        x2range2 = x2ray
      end
# Add a ray to the microlensing map:
      map1[i1ray,i2ray] += 1
    end
# Now, do the same calculation, but ignore planet/moon:
    y1ray =  x1ray -  mbody[1] * (x1ray-x1body[1])/r12
    y2ray =  x2ray -  mbody[1] * (x2ray-x2body[1])/r12
    i1ray = convert(Int64,floor((y1ray - y1min)/dy1))+1
    i2ray = convert(Int64,floor((y2ray - y2min)/dy2))+1
    if (1 <= i1ray <= ny1) && (1 <= i2ray <= ny2) 
      map2[i1ray,i2ray] += 1
    end
  end
# Print out progress of calculation:
  if mod(i1,1000) == 0
    println(i1,"/",nx1)
  end
end
println("Range of grid, x1: ",x1min,' ',x1range1,' ',x1max,' ',x1range2)
println("Range of grid, x2: ",x2min,' ',x2range1,' ',x2max,' ',x2range2)
# Take the ratio of the microlensing maps with & without the planet+star:
map = transpose(map1./map2)
# Now make a plot:
PyPlot.imshow(map, interpolation = "nearest", extent =[y1min_norm,y1max_norm,y2min_norm,y2max_norm])
# Create a pdf file of the map:
PyPlot.savefig("map.pdf", bbox_inches="tight")
return
end
