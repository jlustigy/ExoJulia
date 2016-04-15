push!(LOAD_PATH, "../../../ExoJulia/")
push!(LOAD_PATH, ".")

# Load premade modules
using ExoJulia
using LsqFit
using Optim

include("rv_fitting.jl")

numbers = readdlm("mystery_planet.txt")

solve_rv(numbers)
