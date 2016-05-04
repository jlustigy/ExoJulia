using LsqFit
include("../compute_ttv.jl")
using TTVFaster

"""
    ttvs(data)
Computer parameters for the planets passed in data.

#Arguments
* `data`: an Nx[planet transit data] array, where N is the number of planets

#Return
* `result`: an Nx5 array where each row contains [mass_ratio, period, trans0, ecosw, esinw]
"""
function ttvs(data)
  num_planets = size(data)[1]

  combined_data = [data[1];data[2]]

  function model_func(useless_input_x, params)

    planets = []

    for i=1:div(length(params),5)
      planet = Planet_plane_hk(params[5i-4],params[5i-3],params[5i-2],params[5i-1],params[5i])
      push!(planets,planet)
    end
    ttv1 = zeros(length(data[1]))
    ttv2 = zeros(length(data[2]))
    compute_ttv!(5,planets[1],planets[2],data[1],data[2],ttv1,ttv2)

    return [ttv1;ttv2]
  end

  # errors = ones(length(data)) * (30/(24*60*60)) # error of 30 sec
  x = collect(1:1:length(data))
  p0 = [0.000003, 225, 8445, 0.0001, 0.0001, 0.000003, 365, 8461, 0.0001, 0.0001]
  fit = curve_fit(model_func, x, data, p0)

  return fit.param
end


function test()
  p1 = readdlm("../ttv_planet1.txt")
  p2 = readdlm("../ttv_planet2.txt")

  data = []
  push!(data, p1[:,1])
  push!(data, p2[:,1])

  result = ttvs(data)
  return result
end

println(test())
