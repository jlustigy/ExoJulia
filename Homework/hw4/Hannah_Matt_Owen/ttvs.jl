using LsqFit
include("../compute_ttv.jl")
using TTVFaster
include("MCJulia.jl")
using PyPlot


function model_func(data_start_end, params)
  planets = []
  num_planets = div(length(params),5)
  for i=1:num_planets
    #build an array of each planet object
    planet = Planet_plane_hk(params[5*i-4],params[5*i-3],params[5*i-2],params[5*i-1],params[5*i])
    push!(planets,planet)
  end

  results = zeros(data_start_end[num_planets][2])

  #loop over all the planets comparing each one to each planet it hasn't yet
  #been compared to
  for i=1:num_planets-1
    ttv1 = zeros(data_start_end[i][2]-data_start_end[i][1] +1)
    ttv_data_p1 = collect(1:1:length(ttv1))
    ttv_data_p1 = ttv_data_p1 .* planets[i].period .+ planets[i].trans0
    for j=i+1:num_planets
      ttv2 = zeros(data_start_end[j][2]-data_start_end[j][1] +1)
      ttv_data_p2 = collect(1:1:length(ttv2))
      ttv_data_p2 = ttv_data_p2 .* planets[j].period .+ planets[j].trans0
      compute_ttv!(5,planets[i],planets[j],
          ttv_data_p1,
          ttv_data_p2,
          ttv1,ttv2);
      #add the reuslts to our result array
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


function data_start_end_array(data...)
  data_start_end = [] #store the start and end of each planet's data in the array
  for pd in data
    if length(data_start_end) == 0
      push!(data_start_end, ( 1,length(pd[:])) )
    else
      start = data_start_end[end][2]+1
      push!(data_start_end, ( start, start+length(pd[:]) -1 ) )
    end
  end
  return data_start_end
end

"""
    ttvs(data)
Computer parameters for the planets passed in data.
#Arguments
* `p0`: the initial guess for the data passed in, as an array of 5N elements
        [mass_ratio, period, trans0, ecosw, esinw] for each planet
* `data`: pass in any number of arrays containing transit data for each planet
#Return
* `fit`: fit paramter for the best parameters given the provided data
"""
function ttvs(p0, data...)
  num_planets = length(data)

  if num_planets < 2
    println("Error: must provide data for at least 2 planets")
    return
  end

  combined_data::Array{Float64} = [] #combine all the data into one array, which we'll fit all at once
  data_start_end = [] #store the start and end of each planet's data in the array
  for pd in data
    combined_data = [combined_data;pd[:]]
    if length(data_start_end) == 0
      push!(data_start_end, ( 1,length(pd[:])) )
    else
      start = data_start_end[end][2]+1
      push!(data_start_end, ( start, start+length(pd[:]) -1 ) )
    end
  end

  fit = curve_fit(model_func, data_start_end, combined_data, p0)

  return fit
end

function third_planet()
  p1 = readdlm("../ttv_planet1.txt")
  p2 = readdlm("../ttv_planet2.txt")

  lh = []
  periods = collect(660:5:700)
  for period in periods #loop over the periods
    for t0=1:10 #loop over 10 different t0 values
      t0_val = t0*period/10 + p1[1]
      time = [t0_val, t0_val+period]

      #the initial guess parameter
      param0 = [0.000003, 224, 8445, 0.0001, 0.0001,
                0.000003, 365, 8461, 0.0001, 0.0001,
                0.0000015, period, t0_val, 0.0001, 0.0001 ]
      fit = ttvs(param0, p1,p2,time)

      sigma = estimate_errors(fit, 0.95)[12]
      errors = fit.resid[12]
      lh_val = exp(-errors^2/(2*sigma^2))/(2*pi*sigma^2)^0.5
      push!(lh, (period, lh_val, fit.param))
      end
  end

  #uncomment to plot the curve fit results
  #scatter([x[1] for x in lh], [x[2] for x in lh])
  #show()

  #call the MC MC solver
  perr = [0.000002, 4, 5, 0.00005, 0.00005, #assume same error for each planet
          0.000002, 5, 5, 0.00005, 0.00005,
          0.000002, 5, 5, 0.00005, 0.00005]
  planet_3 = lh[1][3] #data for all three planets
  best = 0.0
  for item in lh
    cur = item[2] #fit value
    if cur > best
      best = cur
      planet_3 = item[3]
    end
  end
  p3 = [planet_3[13], planet_3[13]+planet_3[12]] #[t0, to+period] of planet 3
  combined_data = [p1;p2;p3]
  data_start_end = data_start_end_array(p1,p2,p3)
  ye = ones((1,length(combined_data)))*30/24/3600 #set error to 30 seconds for each measurement
  mc_results = aimc(model_func, data_start_end, planet_3, perr, combined_data, ye)

  @printf("mass of p1=%0.3e +/- %0.2e, period of p1=%4.2f +/- %0.2e\nmass of p2=%0.3e +/- %0.2e, period of p2=%4.2f +/- %0.2e\nmass of p3=%0.3e +/- %0.2e, period of p3=%4.2f +/- %0.2e\n",
          mc_results[1][1],mc_results[1][2],mc_results[2][1],mc_results[2][2],mc_results[6][1],mc_results[6][2],mc_results[7][1],mc_results[7][2],mc_results[11][1],mc_results[11][2],mc_results[12][1],mc_results[12][2])

end

function test()
  p1 = readdlm("../ttv_planet1.txt")
  p2 = readdlm("../ttv_planet2.txt")

  param0 = [0.000003, 224, 8445, 0.0001, 0.0001, 0.000003, 365, 8461, 0.0001, 0.0001]
  fit = ttvs(param0, p1,p2)

  result = fit.param

  return "m1=$(result[1][1]), period1=$(result[2][1]), m2=$(result[6][1]), p2=$(result[7][1])"
end

#println(test())
third_planet()
