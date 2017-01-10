include("./KeplerSolver.jl")
using KeplerSolver
using PyPlot

# Test the solver
max_err = 0
for x in range(0,1000)
  ecc = x/1000.
  M = linspace(0,2*pi, 1000)
  anom = E(ecc, M)
  err = maximum(abs(M - (anom - ecc * sin(anom))))
  if err > max_err
    max_err = err
  end
end
print("Maximum Error in M: ")
println(max_err)
