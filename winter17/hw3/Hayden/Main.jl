#Main Script for hw3

include("CircleOverlap.jl")
include("TransitFit.jl")
include("LayerCake.jl")

#I just discovered the xkcd() plotting option...
xkcd() #Turn off for more "professional" plots

#Problem 1: Circle Overlap
r=linspace(0,2,100);
Circles = CircleOverlap(0.1,r);
plot(r,Circles)
title("Circle Overlap")
ylabel("Uncovered Fraction")
xlabel("Separation of Center")

#Problem 2: Transit Simulation
Time = linspace(0,11,1001);
Flux=TransitSolver(Time,[9, 0.01, 0.625, 0, 0.001, 0.0])
figure()
plot(Time,Flux)
title("Simulated transit light curve")
xlabel("Time")
ylabel("Flux")

#Problem 3: Fitting the transit data
Data=TransitFit("mysteryplanet2.txt", [9, 0.01, 0.625, 0, 0.001, 0.0]);

#Problem 4: Arbitrary Limb Darkening

figure()
title("Limb Darkening")
xlabel("Time")
ylabel("Flux")

plot(linspace(0,10,1001),TransitSolver(linspace(0,10,1001), [10,3,1.3,0,.1,0.]))

plot(linspace(0,10,1001),LayerCake(linspace(0,10,1001), [10,3,1.3,0,.1,0,100],x->x)) #Linear

plot(linspace(0,10,1001),LayerCake(linspace(0,10,1001), [10,3,1.3,0,.1,0,100],x->x^2)) #Quadratic

legend(["No Limb Darkening","Linear","Quadratic"])
