using LsqFit
include("../compute_ttv.jl")
using TTVFaster

"""
    ttvs(data)
Computer parameters for the planets passed in data.

#Arguments
* `data`: pass in any number of arrays containing transit data for each planet

#Return
* `result`: an Nx5 array where each row contains [mass_ratio, period, trans0, ecosw, esinw]
            where N is the number of planets
"""
function ttvs(data...)
  num_planets = length(data)

  if num_planets < 2
    println("Error: must provide data for at least 2 planets")
    return
  end

  combined_data::Array{Float64} = [] #combine all the data into one array, which we'll fit all at once
  data_start_end = [] #store the start and end of each planet's data
  for pd in data
    combined_data = [combined_data;pd[:]]
    if length(data_start_end) ==0
      push!(data_start_end, ( 1,length(pd[:])) )
    else
      start = data_start_end[end][2]+1
      push!(data_start_end, ( start, start+length(pd[:]) -1 ) )
    end
  end

  function model_func(useless_input_x, params) #curve fit makes us take the first parameter

    planets = []

    for i=1:div(length(params),5)
      planet = Planet_plane_hk(params[5*i-4],params[5*i-3],params[5*i-2],params[5*i-1],params[5*i])
      push!(planets,planet)
    end

    results = zeros(length(combined_data))

    for i=1:num_planets-1
      ttv1 = zeros(data_start_end[i][2]-data_start_end[i][1] +1)
      for j=i+1:num_planets
        ttv2 = zeros(data_start_end[j][2]-data_start_end[j][1] +1)
        compute_ttv!(5,planets[i],planets[j],
            combined_data[data_start_end[i][1]:data_start_end[i][2]],
            combined_data[data_start_end[j][1]:data_start_end[j][2]],
            ttv1,ttv2);
        results[data_start_end[i][1]:data_start_end[i][2]] = results[data_start_end[i][1]:data_start_end[i][2]] + ttv1
        results[data_start_end[j][1]:data_start_end[j][2]] = results[data_start_end[j][1]:data_start_end[j][2]] + ttv2
      end
    end

    #convert the results array from ttvs back to periods for curve fit to compare against
    for i=1:num_planets
      d_start = data_start_end[i][1]
      d_end = data_start_end[i][2]
      for j=d_start:d_end
        results[j] = results[j] + (j - d_start)*planets[i].period + planets[i].trans0
      end
    end

    return results
  end

  p0 = [0.000003, 225, 8445, 0.0001, 0.0001, 0.000003, 365, 8461, 0.0001, 0.0001]
  fit = curve_fit(model_func, [], combined_data, p0)
  return fit.param
end


function test()
  p1 = readdlm("../ttv_planet1.txt")
  p2 = readdlm("../ttv_planet2.txt")

  result = ttvs(p1,p2)
  #println("m1=$(result[1]), period1=$(result[2]), m2=$(result[6]), p2=$(result[7])")
  return "m1=$(result[1]), period1=$(result[2]), m2=$(result[6]), p2=$(result[7])"
end

println(test())
