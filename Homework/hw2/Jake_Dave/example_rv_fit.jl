###########################
#
# example_rv_fit.jl
#
# Note: This file is designed to be ran from the
#       Homework/ directory and is furnished with
#       @stest pointers.
#
###########################

push!(LOAD_PATH, "../ExoJulia/")

# Load premade modules
using ExoJulia
using LsqFit
using Optim

# Load new functions
include("utils.jl")
include("orbital_utils.jl")
include("rv.jl")
include("rv_fitting.jl")

# Read-in data, parse
numbers = readdlm("../Resources/mystery_planet.txt");
time = numbers[:,1];
rv = numbers[:,2];
err = numbers[:,end];

# Use the Agol periododram method
best_period = agol_periodogram(numbers, collect(linspace(1.0, 365.0, 2000)))
#@stest agol_periodogram(numbers, collect(linspace(1.0, 365.0, 2000)))

# Solve using curve_fit()
P_best1, e_best1, tp_best1 = solve_rv([time rv err], alg="cf")
#@stest solve_rv([time rv err], alg="cf")

# Solve using optimize()
P_best2, e_best2, tp_best2 = solve_rv([time rv err], alg="opt")

#=
print("~~~~~~~~ curve_fit() ~~~~~~~~\n")
print("Best fit period: $P_best1 days.\n")
print("Best fit eccentricity: $e_best1.\n")
print("Best fit time of periastron passage: $tp_best1 days.\n")
print("~~~~~~~~ optimize() ~~~~~~~~\n")
print("Best fit period: $P_best2 days.\n")
print("Best fit eccentricity: $e_best2.\n")
print("Best fit time of periastron passage: $tp_best2 days.\n")
=#
