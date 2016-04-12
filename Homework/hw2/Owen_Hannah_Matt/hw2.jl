# Add ExoJulia/ to path
push!(LOAD_PATH, "../../../ExoJulia")
using ExoJulia
using PyPlot
using DataFrames

data =  readtable("mystery_planet.txt", separator = ' ', header=false)

time = data[1];
rv = data[2];
err = data[3];

# pick 59 equally spaced points for f, where f is true anamoly
flist = linspace(0,2pi,59)

F = hcat([cos(f) for f in flist], [sin(f) for f in flist], ones(Float64,59), [time - time[1]])'

W = diagm(1/err.^2)

epsilon = inv(F*W*F')