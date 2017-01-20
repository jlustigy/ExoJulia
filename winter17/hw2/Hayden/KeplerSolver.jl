#Solve the Kepler Equation:
# M = E-esin(E)
#using the bisection method. Can't use Newton's method because F'(E) = 0 for E=pi/2
#
#INPUTS:  Mean Anomoly (M) (Degrees)
#         Eccentricity (eccen) (Unitless)
#
#OUTPUTS: Eccentric Anamoly (E) (Degrees)

function KeplerEval(M, eccen, E) #Just Kepler's Equation. Should equal zero.
  E-eccen*sin(E)-M
end

function KeplerSolver(M, eccen) #Bisection-based solver. Solving for E
  if (M == 0.) | (M == 180.) | (M == 360.)
    E = M;
    FE=0.;
    return E
  end

  M = deg2rad(M); #Set input to radians
  Interval = [0.,2*pi] #We want to look for zeros from 0 to 2pi
  EGuess = 0.5*(Interval[1]+Interval[2]) #Find the interval midpoint
  FE = 1. #Initialize F(E), which is the evaluated Kepler's Equation
  while abs(FE) > 1e-6
    if KeplerEval(M, eccen, Interval[1]) * KeplerEval(M, eccen, EGuess) < 0 #Check for a function sign change
      Interval[2] = EGuess;
      FE = KeplerEval(M, eccen, EGuess); #Reevaluate Kepler's Equation for our new guess
    elseif KeplerEval(M, eccen, Interval[2]) * KeplerEval(M, eccen, EGuess) < 0
      Interval[1] = EGuess
      FE = KeplerEval(M, eccen, EGuess); #Same as before
    else
      println("No zero between 0 and 2pi")
      FE = 0.; #Break out of the while loop without using "break" command
    end
    EGuess = 0.5*(Interval[1]+Interval[2])
  end
  return rad2deg(EGuess)
end
