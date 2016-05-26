# Compute magnification for a finite source planet + moon + star
# Positions of bodies (in units of Einstein radius of star):

import PyPlot
include("cgs.jl")

function planet_moon_star(yplanet,y1moon,y2moon,eps1,eps2)
# Computes ray tracing of planet + moon + star.
# Inputs:
#   yplanet: location of planet in image plane (along x-axis)
#     in units of Einstein angle of star
#   y1moon, y2moon: location of moon relative to planet
#     in units of Einstein angle *of planet*
#   eps1: mass ratio of planet to moon
#   eps2: mass ratio of moon to planet
#
# Thanks to Brendan Brewer's raytrace.jl:
# https://github.com/eggplantbren/Microlensing/blob/master/raytrace.jl
#
# Set up coordinates of the 3 bodies (star is at origin)
# in the horizontal direction in the lens plane:
y1body = [0,yplanet,yplanet+y1moon*sqrt(eps1)]
# An example is an Earth around an M dwarf:
# eps1 = MEARTH/0.1/MSUN
# Set up masses of the bodies:
mbody = [1,eps1,eps2*eps1]
# Set up positions of bodies in the vertical axis in the lens plane:
y2body = [0,0,y2moon*sqrt(eps1)]
println("Mass ratio: ",eps1)

# Set up grid of rays in the lens plane:
y1min_norm = -5; y1max_norm = 5
# Choose the rays around the planet, and set the range
# in units of the Einstein radius of the planet:
# horizontal range:
y1min = y1min_norm*sqrt(eps1)+yplanet; y1max = y1max_norm*sqrt(eps1)+yplanet; ny1 = 10000
# vertical range:
y2min_norm = -15; y2max_norm = 15
y2min = y2min_norm*sqrt(eps1); y2max = y2max_norm*sqrt(eps1); ny2 = 10000

# Now 'shoot' these to the source plane:

# We're going to collect rays over a grid in the source planet (y-coordinates):
nu1 = 1000; nu2 = 1000
# Create maps to hold results of ray shooting with & without planet/moon:
map1 = zeros(Int64,(nu1,nu2))
map2 = zeros(Int64,(nu1,nu2))
# Compute the position of the planet as microlensed by star in source plane:
yp= yplanet-1/yplanet
# Set up grid in source plane:
u1min_norm= -5; u1max_norm= 5;
u2min_norm= -2; u2max_norm= 2;
# Compute range of grid in source plane in units of Einstein angle
# of the planet (sqrt(eps1)):
u1min = sqrt(eps1)*u1min_norm+yp; u1max = sqrt(eps1)*u1max_norm+yp
u2min = sqrt(eps1)*u2min_norm; u2max = sqrt(eps1)*u2max_norm
# Spacing of grid in source plane:
du1 = (u1max-u1min)/nu1
du2 = (u2max-u2min)/nu2
# Spacing of grid in lens plane:
dy1 = (y1max-y1min)/ny1
dy2 = (y2max-y2min)/ny2

# Keep track of the range of rays to see if we need to extend the
# grid:
y1range1 = Inf
y1range2= -Inf
y2range1 = Inf
y2range2= -Inf
# Okay, now loop over the grid in the lens plane:
for i1 in 1:ny1
  y1ray = y1min + (i1-0.5)*dy1
  for i2 in 1:ny2
    y2ray = y2min + (i2-0.5)*dy2
# y1ray, y2ray are positions of ray in the lens plane.  Compute
# their angular offset from the 3 bodies:
    r12 = (y1ray-y1body[1])^2+(y2ray-y2body[1])^2
    r22 = (y1ray-y1body[2])^2+(y2ray-y2body[2])^2
    r32 = (y1ray-y1body[3])^2+(y2ray-y2body[3])^2
# Compute the location of ray in source plane using the (vector) lens equation:
    u1ray =  y1ray -  mbody[1] * (y1ray-y1body[1])/r12 -   mbody[2] * (y1ray-y1body[2])/r22 -   mbody[3] * (y1ray-y1body[3])/r32
    u2ray =  y2ray -  mbody[1] * (y2ray-y2body[1])/r12 -   mbody[2] * (y2ray-y2body[2])/r22 -   mbody[3] * (y2ray-y2body[3])/r32
# Find the position in the source-plane grid:
    i1ray = convert(Int64,floor((u1ray - u1min)/du1))+1
    i2ray = convert(Int64,floor((u2ray - u2min)/du2))+1
# See if the ray lands on the grid in the source plane:
    if (1 <= i1ray <= nu1) && (1 <= i2ray <= nu2) 
# Track the range of rays which land on grid:
      if y1ray < y1range1 
        y1range1 = y1ray
      end
      if y1ray > y1range2
        y1range2 = y1ray
      end
      if y2ray < y2range1
        y2range1 = y2ray
      end
      if y2ray > y2range2
        y2range2 = y2ray
      end
# Add a ray to the microlensing map:
      map1[i1ray,i2ray] += 1
    end
# Now, do the same calculation, but ignore planet/moon:
    u1ray =  y1ray -  mbody[1] * (y1ray-y1body[1])/r12
    u2ray =  y2ray -  mbody[1] * (y2ray-y2body[1])/r12
    i1ray = convert(Int64,floor((u1ray - u1min)/du1))+1
    i2ray = convert(Int64,floor((u2ray - u2min)/du2))+1
    if (1 <= i1ray <= nu1) && (1 <= i2ray <= nu2) 
      map2[i1ray,i2ray] += 1
    end
  end
# Print out progress of calculation:
  if mod(i1,1000) == 0
    println(i1,"/",ny1)
  end
end
println("Range of grid, y1: ",y1min,' ',y1range1,' ',y1max,' ',y1range2)
println("Range of grid, y2: ",y2min,' ',y2range1,' ',y2max,' ',y2range2)
# Take the ratio of the microlensing maps with & without the planet+star:
map = transpose(map1./map2)
# Now make a plot:
PyPlot.imshow(map, interpolation = "nearest", extent =[u1min_norm,u1max_norm,u2min_norm,u2max_norm])
# Create a pdf file of the map:
PyPlot.savefig("map.pdf", bbox_inches="tight")
return map
end
