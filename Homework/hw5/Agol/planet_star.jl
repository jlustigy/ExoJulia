# Compute magnification for a finite source planet + star
# Positions of bodies (in units of Einstein radius of star):

using PyPlot
import PyPlot
using CGS

function planet_star(uplanet,mu)
#uplanet = 0.5
x1body = [0,uplanet]
# mu = MEARTH/0.1/MSUN
mbody = [1,mu]
x2body = [0,0]
println("Mass ratio: ",mu)

# Set up grid of rays
x1min_norm = -10; x1max_norm = 10
x1min = x1min_norm*sqrt(mu)+uplanet; x1max = x1max_norm*sqrt(mu)+uplanet; nx1 = 10000
x2min_norm = -30; x2max_norm = 30
x2min = x2min_norm*sqrt(mu); x2max = x2max_norm*sqrt(mu); nx2 = 10000

# Now 'shoot' these to the source plane:

ny1 = 300; ny2 = 300
map1 = zeros(Int64,(ny1,ny2))
map2 = zeros(Int64,(ny1,ny2))
yp= uplanet-1/uplanet
y1min_norm= -10; y1max_norm= 10;
y2min_norm= -10; y2max_norm= 10;
y1min = sqrt(mu)*y1min_norm+yp; y1max = sqrt(mu)*y1max_norm+yp
y2min = sqrt(mu)*y2min_norm; y2max = sqrt(mu)*y2max_norm
dy1 = (y1max-y1min)/ny1
dy2 = (y2max-y2min)/ny2
dx1 = (x1max-x1min)/nx1
dx2 = (x2max-x2min)/nx2
for i1 in 1:nx1
  x1ray = x1min + (i1-0.5)*dx1
  for i2 in 1:nx2
    x2ray = x2min + (i2-0.5)*dx2
    r12 = (x1ray-x1body[1])^2+(x2ray-x2body[1])^2
    r22 = (x1ray-x1body[2])^2+(x2ray-x2body[2])^2
    y1ray =  x1ray -  mbody[1] * (x1ray-x1body[1])/r12 -   mbody[2] * (x1ray-x1body[2])/r22
    y2ray =  x2ray -  mbody[1] * (x2ray-x2body[1])/r12 -   mbody[2] * (x2ray-x2body[2])/r22
    i1ray = convert(Int64,floor((y1ray - y1min)/dy1))+1
    i2ray = convert(Int64,floor((y2ray - y2min)/dy2))+1
    if (1 <= i1ray <= ny1) && (1 <= i2ray <= ny2) 
      map1[i1ray,i2ray] += 1
    end
    y1ray =  x1ray -  mbody[1] * (x1ray-x1body[1])/r12
    y2ray =  x2ray -  mbody[1] * (x2ray-x2body[1])/r12
    i1ray = convert(Int64,floor((y1ray - y1min)/dy1))+1
    i2ray = convert(Int64,floor((y2ray - y2min)/dy2))+1
    if (1 <= i1ray <= ny1) && (1 <= i2ray <= ny2) 
      map2[i1ray,i2ray] += 1
    end
  end
  println(i1,"/",nx1)
end
# Now, convolve with sources:
rsource = logspace(log10(0.25),log10(4.0),81)
nr = length(rsource)
devmax = zeros(nr)
map1_conv = zeros(Float64,(ny1,ny2))
map2_conv = zeros(Float64,(ny1,ny2))
for j=1:nr
  for i1=1:ny1 
    y1 = y1min + (i1-0.5)*dy1
    for i2=1:ny2 
      y2 = y2min + (i2-0.5)*dy2
#  Loop over points within source & add to magnification:
      rn1 = convert(Int64,ceil(rsource[j] * sqrt(mu) / dy1))
      rn2 = convert(Int64,ceil(rsource[j] * sqrt(mu) / dy2))
      area=0
      for j1 = i1-rn1:i1+rn1
        if 1 <= j1 <= ny1
          y1j = y1min + (j1-0.5)*dy1
          for j2 = i2-rn2:i2+rn2
            if 1 <= j2 <= ny2
              y2j = y2min + (j2-0.5)*dy2
              if sqrt((y1-y1j)^2+(y2-y2j)^2) < (rsource[j]*sqrt(mu))
                map1_conv[i1,i2] += float(map1[j1,j2])
                map2_conv[i1,i2] += float(map2[j1,j2])
                area += 1
              end
            end
          end
        end
      end
      if area > 0
        map1_conv[i1,i2] /= area
        map2_conv[i1,i2] /= area
      end
# Compute the deviation of planet map relative to no planet:
      if map2_conv[i1,i2] > 0
        dev = map1_conv[i1,i2]/map2_conv[i1,i2]-1
        if  dev > devmax[j]
          devmax[j] = dev
        end
      end
    end
  end
  println("Convolved ",rsource[j],' ',devmax[j])
  PyPlot.imshow(transpose(map1_conv./map2_conv), interpolation = "nearest", extent =[y1min_norm,y1max_norm,y2min_norm,y2max_norm])
#  println("Hit return")
#  read(STDIN,Char)
#  clf()
end
clf()
plot(rsource,devmax)
println("Maximum deviation: ",devmax)
#map = transpose(map1./map2)
#PyPlot.imshow(map, interpolation = "nearest", extent =[y1min_norm,y1max_norm,y2min_norm,y2max_norm])
#PyPlot.savefig("map.pdf", bbox_inches="tight")
return
end
