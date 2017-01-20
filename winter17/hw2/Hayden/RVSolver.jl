#HW 1b Problem 1
#===========================
FUNCTION NAME: RVSolver
INPUTS:     Time(T), Parameter vector
            Params = [K, w, P, tp, gamma, ecc]
OUTPUTS:     Radial velocity (RV)
NOTES: All angles must be input in radians.
============================#

function RVSolver(t, Params)
  K, w, P, tp, gamma, ecc = Params;
  M = 2*pi*(t-tp)/P;
  E = NewtonKeplerSolver(M, ecc);
  if ecc >= 1
    ecc = 0.9999
  end
  f = 2*atan(((1+ecc)/(1-ecc))^(1/2)*tan(E/2));
  RV = K * cos(w)*cos(f) - K *sin(w)*sin(f)+gamma+K*ecc*cos(w); #Automatically returned
end

#Note, in practice, you pass in K as an independent variable, but I include this
#function just in case you want to solve K from physical values
function KSolver(P, ecc, mpsini, Ms, mp)
  K = (2*pi*6.67e-11/(P*(1-ecc^2)^(3/2)*(Ms+mp)^2))^(1/3)*mpsini
  #Note, K is not very dependent on mp, if Ms>>mp. Rather, it is dependent on mpsini
end
