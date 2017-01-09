#Solve the Kepler Equation:
# M = E-esin(E)
#using the bisection method
#
#INPUTS:  Mean Anomoly (M) Degrees
#         Eccentricity (eccen)
#
#OUTPUTS: Eccentric Anamoly (E)

function KeplerEval(M, eccen, E)
  E-eccen*sin(E)-M
end

function KeplerSolver(M, eccen)
  M = deg2rad(M);
  Interval = [0.,2*pi]
  MidPoint = 0.5*(Interval[1]+Interval[2])
  while abs(KeplerEval(M, eccen, MidPoint)) > 1e-6
    if KeplerEval(M, eccen, Interval[1]) * KeplerEval(M, eccen, MidPoint) < 0
      Interval[2] = MidPoint
    elseif KeplerEval(M, eccen, Interval[2]) * KeplerEval(M, eccen, MidPoint) < 0
      Interval[1] = MidPoint
    else
      print("No zero between 0 and 2pi")
      break
    end
    MidPoint = 0.5*(Interval[1]+Interval[2])
  end
  return rad2deg(MidPoint)
end
