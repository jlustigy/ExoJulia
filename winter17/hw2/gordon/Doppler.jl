module Doppler
export K, vrad

include("./KeplerSolver.jl")
using KeplerSolver

# computes the radial velocity at time(s) t
function vrad(t, P, ecc, tp, gamma, omega, K)
  h = K * cos(omega)
  c = -1 * K * sin(omega)
  v0 = gamma + K * ecc * cos(omega)
  M = (2 * pi / P) * (t - tp)
  tanf2 = ((1 + ecc) / (1 - ecc))^(0.5) * tan(E(ecc, M) ./ 2.)
  f = 2 * atan(tanf2)
  return h * cos(f) + c * sin(f) + v0
end
end
