# This is a "wrapper" for ttv_nplanet.jl for calling the 
# function with the LsqFit function curve_fit.jl
# The inner two planets transit; the outer does not.

#include("../../TTVFaster/Julia/compute_ttv.jl")
include("ttv_nplanet.jl")

function ttv_wrapper2(tt,param)
# These lines need modification for different choices of parameters:
nplanet = 2
ntrans = [38,24]
jmax = 5
# Call ttv_nplanet:
ttv = ttv_nplanet(nplanet,jmax,ntrans,param)
# We measure transit times, not TTVs, so add
# back in the linear ephemeris:
n1 = ntrans[1]
t01 = param[3]
per1 = param[2]
ttv1 = collect(linspace(t01,t01+per1*(n1-1),n1))
for i=1:n1
 ttv1[i]+= ttv[1,i]
end
n2 = ntrans[2]
t02 = param[8]
per2 = param[7]
ttv2 = collect(linspace(t02,t02+per2*(n2-1),n2))
for i=1:n2
  ttv2[i] += ttv[2,i]
end
# If transit times of additional planets were observable
# these would need to be added in.
#println("param2: ",param)
return [ttv1;ttv2]
end
