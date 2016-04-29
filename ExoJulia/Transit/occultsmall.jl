function occultsmall(p,c1,c2,c3,c4,z)
# This routine approximates the lightcurve for a small 
# planet assuming Claret (2000)'s non-linear limb-darkening law. 
# (See section 5 of Mandel & Agol (2002) for details; note that \pi is
# missing from denominator in paper):
# Please cite Mandel & Agol (2002) if you make use of this routine.
#
# Input:
#  p      ratio of planet radius to stellar radius
#  c1-c4  non-linear limb-darkening coefficients
#  z      impact parameters (positive number normalized to stellar
#         radius)- this is an array which MUST be input to the routine
#
# Output:
#  mu     flux relative to unobscured source for each z
#
nz = length(z)
mu = ones(nz)
norm=pi*(1-c1/5-c2/3-3*c3/7-c4/2)
for i=1:nz
  if (z[i] > (1-p)) && (z[i] < 1+p)
    x=1-(z[i]-p)^2
    tmp=(1-c1*(1-0.8*x^(1//4))-c2*(1-2//3*sqrt(x))
          -c3*(1-4//7*x^(3//4))-c4*(1-0.5*x))
    mu[i]=1-tmp*(p^2*acos((z[i]-1)/p)-(z[i]-1)*sqrt(p^2-(z[i]-1)^2))/norm
  end
  if (z[i] <= (1-p)) && (z[i] != 0)
    mu[i]=1-pi*p^2*iofr(c1,c2,c3,c4,z[i],p)/norm
  end
  if z[i] == 0
    mu[i]=1-pi*p^2/norm
  end
end
return mu
end
      
function iofr(c1,c2,c3,c4,r,p)
sig1=sqrt(sqrt(1-(r-p)^2))
sig2=sqrt(sqrt(1-(r+p)^2))
intensity=(1-c1*(1+(sig2^5-sig1^5)/5/p/r)
            -c2*(1+(sig2^6-sig1^6)/6/p/r)
            -c3*(1+(sig2^7-sig1^7)/7/p/r)
            -c4*(p^2+r^2))
return intensity
end
