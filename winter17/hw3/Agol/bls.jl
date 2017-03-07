#
#
function bls(n::Integer,t::Vector,x::Vector,nf::Int64,fmin::Float64,df::Float64,nb::Int64,qmi::Float64,qma::Float64)
#
#------------------------------------------------------------------------
#     >>>>>>>>>>>> This routine computes BLS spectrum <<<<<<<<<<<<<<
#
#        [ see Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369 ]
#------------------------------------------------------------------------
#
#     Input parameters:
#     ~~~~~~~~~~~~~~~~~
#
#     n    = number of data points
#     t    = array {t(i)}, containing the time values of the time series
#     x    = array {x(i)}, containing the data values of the time series
#     u    = temporal/work/dummy array, must be dimensioned in the 
#            calling program in the same way as  {t(i)}
#     v    = the same as  {u(i)}
#     nf   = number of frequency points in which the spectrum is computed
#     fmin = minimum frequency (MUST be > 0)
#     df   = frequency step
#     nb   = number of bins in the folded time series at any test period       
#     qmi  = minimum fractional transit length to be tested
#     qma  = maximum fractional transit length to be tested
#
#     Output parameters:
#     ~~~~~~~~~~~~~~~~~~
#
#     p    = array {p(i)}, containing the values of the BLS spectrum
#            at the i-th frequency value -- the frequency values are 
#            computed as  f = fmin + (i-1)*df
#     bper = period at the highest peak in the frequency spectrum
#     bpow = value of {p(i)} at the highest peak
#     depth= depth of the transit at   *bper*
#     qtran= fractional transit length  [ T_transit/bper ]
#     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
#     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]
#
#
#     Remarks:
#     ~~~~~~~~ 
#
#     -- *fmin* MUST be greater than  *1/total time span* 
#     -- *nb*   MUST be lower than  *nbmax* 
#     -- Dimensions of arrays {y(i)} and {ibi(i)} MUST be greater than 
#        or equal to  *nbmax*. 
#     -- The lowest number of points allowed in a single bin is equal 
#        to   MAX(minbin,qmi*N),  where   *qmi*  is the minimum transit 
#        length/trial period,   *N*  is the total number of data points,  
#        *minbin*  is the preset minimum number of the data points per 
#        bin.
#     
# ========================================================================
#
#      implicit real*8 (a-h,o-z)
#
#      dimension t(*),x(*),u(*),v(*),p(*)
u = zeros(Float64,n)
v = zeros(Float64,n)
y   = zeros(Float64,nb)
ibi = zeros(Int64,nb)
#      dimension y(2000),ibi(2000)
#
minbin = 5
nbmax  = 2000
if (nb > nbmax)
  println(" NB > NBMAX !!")
  return
end
tot=t[n]-t[1]
if (fmin < 1.0/tot) 
  println(" fmin < 1/T !!")
  return
end
#------------------------------------------------------------------------
#
rn=float(n)
#      kmi=idint(qmi*float(nb))
kmi=floor(Int64,qmi*nb)
if (kmi < 1)
  kmi=1
end
kma=floor(Int64, qma*float(nb))+1
kkmi=floor(Int64,rn*qmi)
if (kkmi < minbin) 
  kkmi=minbin
end
bpow = 0.0
bper = 0.0
depth= 0.0
qtran = 0.0
in1 = 0
in2 = 0
#
# =================================
#     Set temporal time series
# =================================
#
s=0.0
t1=t[1]
for i=1:n
  u[i]=t[i]-t1
  s += x[i]
end
s=s/rn
for i=1:n
  v[i]=x[i]-s
end
#
#******************************
#     Start period search     *
#******************************
#  
f0 = zeros(nf)
p  = zeros(nf)
for jf=1:nf
  f0[jf] = fmin + df*float(jf-1)
  p0 = 1.0/f0[jf]
#
# ======================================================
#     Compute folded time series with  *p0*  period
# ======================================================
#
  for j=1:nb
    y[j]   = 0.0
    ibi[j] = 0
  end
#
  for i=1:n  
    ph     = u[i]*f0[jf]
    ph     = ph-floor(Int64,ph)
    j      = 1 + floor(Int64,nb*ph)
    ibi[j] = ibi[j] + 1
    y[j]   =   y[j] + v[i]
  end
#
# ===============================================
#     Compute BLS statistics for this period
# ===============================================
#
  power=0.0
# 
  jn1 = 0
  jn2 = 0
  rn3 = 0.0
  s3 = 0.0
  for i=1:nb
    s     = 0.0
    k     = 0
    kk    = 0
    nb2   = i+kma
    if (nb2 > nb) 
      nb2=nb
    end
    for j=i:nb2
      k     = k+1
      kk    = kk+ibi[j]
      s     = s+y[j]
      if (k >= kmi)
        if (kk >= kkmi)
          rn1   = float(kk)
          pow   = s*s/(rn1*(rn-rn1))
          if (pow >= power)
            power = pow
            jn1   = i
            jn2   = j
            rn3   = rn1
            s3    = s
          end  
        end  
      end  
    end  
  end
#
  power = sqrt(power)
  p[jf] = power
#
  if (power >= bpow)
    bpow  =  power
    in1   =  jn1
    in2   =  jn2
    qtran =  rn3/rn
    depth = -s3*rn/(rn3*(rn-rn3))
    bper  =  p0
  end
#
end
#
#     p    = array {p(i)}, containing the values of the BLS spectrum
#            at the i-th frequency value -- the frequency values are
#            computed as  f = fmin + (i-1)*df
#     bper = period at the highest peak in the frequency spectrum
#     bpow = value of {p(i)} at the highest peak
#     depth= depth of the transit at   *bper*
#     qtran= fractional transit length  [ T_transit/bper ]
#     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
#     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]

return p,bper,bpow,depth,qtran,in1,in2,f0
end

#
#
