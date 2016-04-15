push!(LOAD_PATH, "../../../ExoJulia/")
push!(LOAD_PATH, ".")

# Load premade modules
using ExoJulia
using LsqFit
using Optim

include("fit_rv.jl")

numbers = readdlm("mystery_planet.txt")

fit_rv(numbers)
