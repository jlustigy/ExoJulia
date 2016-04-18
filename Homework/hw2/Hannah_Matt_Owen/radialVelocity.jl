# Add ExoJulia/ to path
push!(LOAD_PATH, "../../../ExoJulia")
using ExoJulia
using LsqFit


"""
    true_anomaly(ecc::Float64, m::Float64)

Calculate the true anomaly from the given eccentricity and mean anomaly.

#Arguments
* e   The eccentricity of the orbit [0<e<1]
* m   The mean anomaly in radians [0<m<2pi]

#Return
* f   The true anomaly [0<f<2pi] on success. -1 on failure
"""
function true_anomaly(ecc::Float64, m::Float64)
  #tan(f/2)=sqrt((1+e)/(1-e))*tan(ea/2)
  ea = ExoJulia.Orbit.kepler_solve!(m,ecc)
  f = atan2(sqrt(1-ecc^2)*sin(ea),cos(ea)-ecc)
end

"""
    get_optimal_rv_parameters(data)

Compute the optimal period, time of periastron, and eccentricity for the given
data.

#Arguments
* `data::Array{Float64,3}`: an Nx3 array of data with format [time; RV; error]
* `numPlanets::Int`: the number of planets in the system
* `min_period::Float64`: the minimum period to guess
* `max_period::Float64`: the maximum perdiod to guess

#Return
* best_params   An array containing [period, ecc, tp, h, c, v0]
"""
function get_optimal_rv_parameters(data, numPlanets, min_period, max_period)
  time = data[:,1];
  rv = data[:,2];
  err = data[:,3];

  F = zeros(Float64,numPlanets*2+2,length(time));

  best_params = [0,0,0,0,0,0]
  cur_best = Inf
  for p in linspace(min_period,max_period,10)
    for ecc in linspace(0,0.5,10)
      for tp in linspace(0,p,round(p)) #in day increments
        #initialize F
        for j=1:length(time)
          ma = 2pi/p*(time[j]-tp); #calculate the mean anomaly
          f = true_anomaly(ecc,ma); #get the true anomaly
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
        #TODO this should depend on the numPlanets
        beta = rv' * W * (F') * epsilon
        h = beta[1]
        c = beta[2]
        v0 = beta[end-1]

        p0 = [h,c,v0]

        """
            linear_rv(time, p)

        For each time value calculate the model RV values for the given h, c, and v0
        values in p.
        """
        function linear_rv(time, vals)
          rv = zeros(Float64, length(time))
          for i=1:length(time)
            ma = 2pi/p*(time[i]-tp); #calculate the mean anomaly
            f = true_anomaly(ecc,ma); #get the true anomaly
            rv[i] = vals[1]*cos(f)+vals[2]*sin(f)+vals[3]
          end
          return rv
        end

        fit = curve_fit(linear_rv, time, rv, 1./err.^2, p0);

        fit_amt = norm(fit.param - p0);

        if fit_amt < cur_best
          cur_best = fit_amt
          best_params = [p,ecc,tp,h,c,v0]
        end
      end
    end
  end
  return best_params
end



"""
    estimate_period(time, rv, min_p, max_p)

Estimate the orbital period from the Agol method.
"""
function estimate_period(time, rv, min_p, max_p)
  period = 0
  prev = Inf
  for p in linspace(min_p, max_p, 1000)
    sum = 0
    data_temp = zeros(Float64, length(time), 2)
    data_temp[:,1] = mod(time,p)
    data_temp[:,2] = rv
    sorted = sortrows(data_temp, by=x->x[1])

    for i=2:length(time)
      sum += (data_temp[i][1] - data_temp[i-1][1])^2
    end

    if sum < prev
      prev = sum
      period = p
    end
  end
  return period
end



function test()
  data =  readdlm("./mystery_planet.txt")
  result = get_optimal_rv_parameters(data,1,110,120);
  println("-----------------------------------")
  println("optimal values are:\nperiod=$(result[1]) days\necc=$(result[2])\ntp=$(result[3])")
end

test()
#@stest get_optimal_rv_parameters(data,1,100,150);
