module KeplerSolver

export E

using Roots

# Kepler solver
function E(ecc, M)
  tor = zeros(length(M))
  for (i, m) in enumerate(M)
    f(x) = x - ecc * sin(x) - m
    tor[i] = fzero(f, m)
  end
  return tor
end
end
