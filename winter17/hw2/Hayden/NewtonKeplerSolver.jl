#Solve the Kepler Equation:
# M = E-esin(E)
#using the bisection method. Can't use Newton's method because F'(E) = 0 for E=pi/2
#
#INPUTS:  Mean Anomoly (M) (Radians)
#         Eccentricity (eccen) (Unitless)
#
#OUTPUTS: Eccentric Anamoly (E) (Radians)

function NewtonKeplerSolver(M, eccen) #Newton-based solver. Solving for E
  if (M == pi) | (M == 0) | (M == 2*pi)
    E = M;
    return E
  else
    E = M; #Initial guess
    dg = E -> 1.-eccen*cos(E);
    g =  E -> E-eccen*sin(E)-M;
    i=0;
    while (maximum(abs(g(E)./dg(E))) > 1e-6) & (i < 100000)
      E = E - g(E)./dg(E);
      i+=1;
    end
    if i == 100000 #Check to see if the max number of iterations was reached
      println("Warning: Tolerance not reached")
    end
    return E
  end

end
