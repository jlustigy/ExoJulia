#=
HW 3 Problem 1
Author: Hayden Smotherman
OBJECTIVE: Compute the overlap of two circles as a function of their separation.
INPUTS: Small Circle Radius (K) in terms of occulted circle radius.
        Separation (r) in terms of occulted circle radius.
OUTPUS: Uncovered area of the large circle (0<A<1)
NOTES: Inputs must be floating point numbers. r may be an array.
=#

function CircleOverlap(K::Float64,r)

  #Turn r into an array, if it isn't one already
  if ~(isa(r,Array) | isa(r,LinSpace))
    r=[r];
  end
  #Check to make sure there isn't an error in input
  if (K<0.) | any(r.<0.)
    error("InputError: Radii and separation cannot be negative")
  end
  A = -1*ones(size(r));

  for i=1:maximum(size(r))
    #Handle the special cases first.
    if (r[i].<=(1.-K))&(K<1.) #Small circle is fully inside of big circle
      A[i]=1-K^2;
    elseif (r[i]>=(1.+K)) #No overlap at all
      A[i]=1.;
    elseif (r[i]<=(K-1.))&(K>0.) #Big circle fully covers small circle
      A[i]=0;
    else #Some level of intersection
      A[i] = K^2*acos((r[i].^2+K^2-1)/(2*r[i]*K))+acos((r[i].^2+1-K^2)/(2*r[i]))-(1/2)*sqrt((-r[i]+K+1)*(r[i]+K-1)*(r[i]-K+1)*(r[i]+K+1));
      A[i] = 1-A[i]/pi;
    end

  #Make sure A is between 0 and 1
    if A[i]<0
      A[i] = 0
      println("WARNING: Manually increased area values to fit within bounds 0<A<1")
    elseif A[i]>1
      println("WARNING: Manually reduced area values to fit within bounds 0<A<1")
      A[i] = 1
    end
  end #End For
  return A
end
