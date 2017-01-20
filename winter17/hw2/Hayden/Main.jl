#Place this script in your home directory, or alter the cd()
#cd("ExoJulia/winter17/hw2/Hayden")

include("RVFit.jl")
tic()
BestFit=RVFit("mystery_planet.txt",[300, 1.9*pi, 200, 71, -10, 0.3])
toc()

#Runs in ~0.5 seconds w/ plotting. ~0.25 w/o plotting.
