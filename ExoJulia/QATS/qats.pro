#+
# NAME: QATS
#
#
#
# PURPOSE: Detect box-like events of fixed width q (in index units)
# and depth in a sequence, x, of zero-median values. The interval between
# successive boxes is optimally determined constrained between
# DeltaMin and DeltaMax (in index units).  The total box signal is
# returned. The indices of the instants of the optimal transits can
# also be returned.  
#
#
# CATEGORY:
#
#
#
# CALLING SEQUENCE:
# qats(x,q,DeltaMin,DeltaMax[,SBest=SBest,MBest=MBest,indices=indices])
#
#
#
# INPUTS: x - array to be searched, q - box duration, DeltaMin -
#             minimum interval, DeltaMax - maximum interval
#
#
#
# OPTIONAL INPUTS: 
#                     
#                         
#
#
#
# KEYWORD PARAMETERS:
#
#
#
# OUTPUTS:
#
#
#
# OPTIONAL OUTPUTS:SBest - Maximum signal, MBest - number of transits
#                          corresponding to optimal signal, indices  -
#                          indices of instants of boxes (corresponding
#                          to indices of x)
#
#
#
# COMMON BLOCKS:
#
#
#
# SIDE EFFECTS:
#
#
#
# RESTRICTIONS:
#
#
#
# PROCEDURE:
#
#
#
# EXAMPLE:
#
#
#
# MODIFICATION HISTORY: First written 10/4/2012 - JAC
#
#-

function shConvol(x::Vector,q::Integer)
# Convolves the light curve with a trial integer transit duration of q.
# Input:
#   x - data vector with zero mean, Gaussian noise, and negative transits
#       of unknown, but uniform, depth
#   q - transit duration (integer)
# Output:
#  cv - convolved vector which is shorter in length than x by q-1
#
  parity = -1
  nx = length(x)
  cv = zeros(nx-q+1)
  for i=1:nx-q
    for j=0:q-1
      cv[i] += parity*x[i+j]
    end
  end
  return cv
end

function omegaBounds(ml,M,DeltaMin,DeltaMax,N,q)
# Returns the bounds on the omega function:
  return,[((ml-1L)*DeltaMin > N-DeltaMax-(M-ml)*DeltaMax ? (ml-1L)*DeltaMin : N-DeltaMax-(M-ml)*DeltaMax),
          (DeltaMax-q+(ml-1L)*DeltaMax < N-q-(M-ml)*DeltaMin ? DeltaMax-q+(ml-1L)*DeltaMax : N-q-(M-ml)*DeltaMin)]

end

function gammaBounds(m,n,DeltaMin,DeltaMax,q)
# Returns the bounds on the gamma function:
  return,[((m-1L)*DeltaMin > n-DeltaMax ? (m-1L)*DeltaMin : n-DeltaMax),
          (DeltaMax-q+(m-1L)*DeltaMax < n-DeltaMin ? DeltaMax-q+(m-1L)*DeltaMax : n-DeltaMin)]

end

function computeSmnRow!(d,M,ml,DeltaMin,DeltaMax,q,Smn)

  FORWARD_FUNCTION computeSmnRow

  N = n_elements(d)

  if ml != 2 then 
    computeSmnRow!(d,M,ml-1,DeltaMin,DeltaMax,q,Smn)
  end
  oB = omegaBounds(ml,M,DeltaMin,DeltaMax,N,q)
  for nl = oB[0]:oB[2]
     if ml == 1 then 
       Smn[ml-1L,nl]=d[nl]
     else
        gB = gammaBounds(ml-1,nl,DeltaMin,DeltaMax,q)
        Smn[ml,nl]=d[nl]+maximum(Smn[ml-1,gb[1]:gb[2]])
     end
  end
return
end
  

function computeSmn(d,M,DeltaMin,DeltaMax,q,Smn=Smn)
  # Returns SBest
  # Set Smn to keep full matrix (for index retrieval)

  N = n_elements(d)
  if arg_present(Smn) then begin
     Smn = zeros(M,N-q+1L)
     computeSmnRow(d,M,M,DeltaMin,DeltaMax,q,Smn)
     oB = omegaBounds(M,M,DeltaMin,DeltaMax,N,q)
     return,max(Smn[M,oB[1]:oB[2]])/sqrt(1.0*M*q)
  endif

  Smn = dindgen(3,N-q+1)
  r = 1
  for ml = 1L, M, 1 do begin
     oB = omegaBounds(ml,M,DeltaMin,DeltaMax,N,q)
     if ml ne 1 then begin
        for nl = oB(0), oB(1), 1 do begin
           gB = gammaBounds(ml-1,nl,DeltaMin,DeltaMax,q)
           rm1 = ((r-1) mod 3 < 0 ? 3+(r-1) mod 3 : (r-1) mod 3)
           Smn[r,nl]=d[nl]+max(Smn[rm1,gb(0):gb(1)])
        endfor
     endif else for nl = oB(0), oB(1), 1 do Smn[r,nl] = d[nl]
     r = (r+1) mod 3
  endfor

  oB = omegaBounds(M,M,DeltaMin,DeltaMax,N,q)
  rm1 = ((r-1) mod 3 < 0 ? 3+(r-1) mod 3 : (r-1) mod 3)
  return,max(Smn[rm1,oB(0):oB(1)])/sqrt(1.0*M*q)
  
end

function optIndices,M,DeltaMin,DeltaMax,N,q,Smn

  indices = lindgen(M)
  oB = omegaBounds(M,M,DeltaMin,DeltaMax,N,q)
  mx = max(Smn[M-1,oB(0):oB(1)],index) & indices[M-1] = index+oB(0)
  for ml = M-1L, 1, -1 do begin
     gB = gammaBounds(ml+0L,indices[ml],DeltaMin,DeltaMax,q)
     mx = max(Smn[ml-1,gB(0):gB(1)],index) & indices[ml-1] = index+gB(0)
  endfor

  return, indices

end

PRO qats,x,q,DeltaMin,DeltaMax,SBest=SBest,MBest=MBest,indices=indices

  # if indices set, times returned (memory intensive)
  d = shConvol(x,q)
  N = n_elements(d)
  
  MMin = floor((N+q-1.0)/DeltaMax)+0L
  MMax = floor((N-q-0.0)/DeltaMin)+1L

  
  SBest = 0d0

  for M = MMin, MMax, 1L do begin
     S = computeSmn(d,M,DeltaMin+0L,DeltaMax+0L,q+0L)
     if S gt Sbest then begin
        SBest = S
        MBest = M
     endif
  endfor
  
  if arg_present(indices) then begin
     
     S = computeSmn(d,MBest,DeltaMin+0L,DeltaMax+0L,q+0L,Smn=Smn)
     indices = optIndices(MBest,DeltaMin+0L,DeltaMax+0L,N+0L,q+0L,Smn)
  endif

end
