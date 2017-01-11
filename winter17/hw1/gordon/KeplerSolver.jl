module KeplerSolver

export E

using Roots

# Kepler solver
function E(ecc, M)
  tor = zeros(length(M))
  for (i, m) in enumerate(M)
    f(x) = x - ecc * sin(x) - m
    fp(x) = 1 - ecc * cos(x)
    tor[i] = nsolve(f, fp, m, 1e-6)
  end
  return tor
end

function nsolve(f, fp, x, tolerance)
  delta = tolerance + 1
  while delta > tolerance
    last = x
    x = x - f(x)/fp(x)
    delta = abs(x - last)
  end
  return x
end
end
