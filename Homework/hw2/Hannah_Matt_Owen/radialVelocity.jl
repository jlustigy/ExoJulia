# Add ExoJulia/ to path
push!(LOAD_PATH, "../../../ExoJulia")
using ExoJulia
using LsqFit


"""
estimate_period(time, rv, min_p, max_p, num_points)

Estimate the orbital period from the Agol method.
"""
function estimate_period(time, rv, min_p, max_p, num_points)
  period = 0
  prev = Inf
    for p in linspace(min_p, max_p, num_points)
    sum = 0
    data_temp = zeros(Float64, length(time), 2)
    data_temp[:,1] = mod(time,p)
    data_temp[:,2] = rv
    sorted = sortrows(data_temp, by=x->x[1])

    for i=2:length(time)
        sum += (sorted[i,2] - sorted[i-1,2])^2
    end

    if sum < prev
      prev = sum
      period = p
    end
  end
  return period
end



"""
    true_anomaly(ecc::Float64, m::Float64)

Calculate the true anomaly from the given eccentricity and mean anomaly.

#Arguments
* e   The eccentricity of the orbit [0<e<1]
* m   The mean anomaly in radians [0<m<2pi]

#Return
* f   The true anomaly [0<f<2pi] on success. -1 on failure
"""
function true_anomaly(ecc, m)
  #tan(f/2)=sqrt((1+e)/(1-e))*tan(ea/2)
  ea = ExoJulia.Orbit.kepler_solve(m,ecc)
  f = atan2(sqrt(1-ecc^2)*sin(ea),cos(ea)-ecc)
end


"""
    get_optimal_rv_parameters(data)

Compute the optimal period, time of periastron, and eccentricity for the given
data.

#Arguments
* `data::Array{Float64,3}`: an Nx3 array of data with format [time; RV; error]
* `numPlanets::Int=1`: the number of planets in the system
* `period_guess::Float64`: the best period guess

#Return
* best_params   An array containing [period, ecc, tp]
"""
function get_optimal_rv_parameters(data, period_guess; numPlanets::Int=1)
  time = data[:,1];
  rv = data[:,2];
  err = data[:,3];


  p0 = [period_guess, 0.5, rand(time)] #[period ecc tp]
  function model_rv(time, vals)
    #vals format [period ecc tp]
    F = zeros(Float64,numPlanets*2+2,length(time));
    for j=1:length(time)
      ma = 2pi/vals[1]*(mod(time[j],vals[1])-vals[3]); #calculate the mean anomaly
      f = true_anomaly(vals[2],ma); #get the true anomaly
      for i=1:numPlanets
        F[i*2-1,j] = cos(f);
        F[i*2,j] = sin(f);
        F[end-1,j] = 1;
        F[end,j] = time[j]-time[1];
      end
    end

    W = diagm(1./err.^2); #the error matrix
    epsilon = inv(F * W * (F'));

    #calculate beta which gives {h, c, v0, d}
    beta = rv' * W * (F') * epsilon
    h = beta[1]
    c = beta[2]
    v0 = beta[end-1]

    rv = zeros(Float64, length(time))
    for i=1:length(time)
      ma = 2pi/vals[1]*(mod(time[i],vals[1])-vals[3]); #calculate the mean anomaly
      f = true_anomaly(vals[2],ma); #get the true anomaly
      rv[i] = h*cos(f)+c*sin(f)+v0
    end
    return rv
  end

  fit = curve_fit(model_rv, time, rv, 1./err.^2, p0);

  return fit.param

end


function test()
  data =  readdlm("./mystery_planet.txt")
  time = data[:,1];
  rv = data[:,2];
  p_guess = estimate_period(time, rv, 100, 150, 10000)
  result = get_optimal_rv_parameters(data,p_guess);
  # println("-----------------------------------")
  # println(p_guess)
  # println(result)
  # println("optimal values are:\nperiod=$(result[1]) days\necc=$(result[2])\ntp=$(result[3])")
end

# test()
#@stest test()
