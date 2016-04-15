# Wrote routine to solve Kepler's equation  Eric Agol 3/30/16
#
# Follows Murray & Dermott "Solar System Dynamics", equations (2.62)-(2.64)
# A quartic-order solver based on Danby (Fundamentals of Celestial Mechanics).
# Solves Kepler's equation for an elliptical orbit.
#

#@stest [[kepler_solve(M,ecc) for M in linspace(0,2pi,100)] for ecc in linspace(0,0.999,100)]

function kepler_solve(M::Float64,ecc::Float64)
#
# Input:
#  M  mean anomaly = n(t-t0) = 2\pi/P*(t-t0), where P is period, n
#                    mean motion, and t0 time of pericenter passage
# ecc eccentricity 0 <= ecc < 1
#
# Output:
#  E  eccentric anomaly
#
# First, reduce range of M to [0,2*pi):
  Mred =  mod(M,2pi)
# Initial guess:
  E = Mred + sign(sin(Mred))*0.85*ecc
# Set the initial parameter to make it into the while loop:
  di3 = 1.0
# Set the tolerance:
  tol = 1e-12
# Begin loop for quartic solver:
  niter = 0
  while (abs(di3) > tol) & (niter < 30)
# Define e*sin(E) and e*cos(E) so that these do not need to
# be recomputed:
    SE = ecc*sin(E); CE = ecc*cos(E)
# Evaluate function and first three derivatives:
    f_of_E = E-SE-Mred; df_of_E=1-CE ; d2f_of_E=SE; d3f_of_E=CE
# Evaluate the successive estimates to the change of E:
    di1 = -f_of_E/df_of_E
    di2 = -f_of_E/(df_of_E+0.5*di1*d2f_of_E)
    di3 = -f_of_E/(df_of_E+0.5*di2*d2f_of_E+di2^2/6.*d3f_of_E)
#    di3 = -f_of_E/df_of_E
#    di3 = -(E-ecc*sin(E))/(1.0-ecc*cos(E))
# Finally, compute next estimate of the eccentric anomaly, E:
    E+=di3
    niter = niter + 1
#    println(di3)
  end
# Check that E satisfies Kepler's equation:
#println(E - ecc*sin(E))
# Add back the 2*pi portion:
E += (M-Mred)
# Output E:
if(niter == 30)
  println("Error: Reached niter = 30")
end
#if (abs(M - E + ecc*sin(E)) > tol) | (niter > 10) then
#  println(ecc, ' ',M, ' ', E-ecc*sin(E), ' ',M-E+ecc*sin(E),' ', niter)
#end
# Output the eccentric anomaly as the result:
E
end
