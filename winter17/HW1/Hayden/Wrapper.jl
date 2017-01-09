#Wrapper code to run the bisection-based Kepler Solver Function

include("KeplerSolver.jl")

#Let's solve Kepler's Equation for several values of M and e:

#M = 25 (degrees), e=0.01
println(string("E = ",KeplerSolver(25., 0.01)), " degrees")

#M = 90 (degrees), e=0.9
println(string("E = ",KeplerSolver(90., 0.9)), " degrees")

#For loop to test many, many values. Use it if you want to.
for M=0:5:360
  for e=0:0.05:1
    println(string("E = ",KeplerSolver(M, e)), " degrees, given: M=",M," and e=",e)
  end
end
