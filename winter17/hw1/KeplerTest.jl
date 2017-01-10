include("./KeplerSolver.jl")
using KeplerSolver
using PyPlot

# Test the solver

ecc = 0.1
M = linspace(0,2*pi, 1000)
anom = E(ecc, M)
plot(M, abs(M - (anom - ecc * sin(anom))))
title("Error")
ylabel(L"$|M - (E - c\ \sin\ E)|$")
xlabel(L"$M$")
