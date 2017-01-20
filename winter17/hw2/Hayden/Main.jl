#Place this script in your home directory, or alter the cd()
#cd("ExoJulia/winter17/hw2/Hayden")

include("RVFit.jl")
tic()
BestFit,Residual=RVFit("mystery_planet.txt",[300, 1.9*pi, 200, 71, -10, 0.3],300);
toc()

println("Semiamplitude (K) = $(BestFit[1])")
println("Longitude of Pericenter = $(BestFit[2])")
println("Period (P) = $(BestFit[3])")
println("Time of Pericenter Passage = $(BestFit[4])")
println("Gamma = $(BestFit[5])")
println("Eccentricity = $(BestFit[6])")
println("R^2 = $(sum(Residual.^2))")
#Runs in ~0.5 seconds w/ plotting. ~0.25 w/o plotting. (For 10 iterations)
