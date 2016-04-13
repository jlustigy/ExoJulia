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
  ea = ExoJulia.Orbit.kepler_solve(m,ecc)
  f = atan2(sqrt(1-ecc^2)*sin(ea),cos(ea)-ecc)
end

"""
    get_optimal_rv_parameters(data)

Compute the optimal period, time of periastron, and eccentricity for the given
data.

#Arguments
* `data::Array{Float64,3}`: an Nx3 array of data with format [time; RV; error]

#Return
* best_params   An array containing [period, ecc, tp, h, c, v0]
"""
function get_optimal_rv_parameters(data, numPlanets)
  #TODO use Eric's algorithm to guess P

  time = data[:,1];
  rv = data[:,2];
  err = data[:,3];

  F = zeros(Float64,numPlanets*2+2,length(time)); #TODO, should 4 be passed in?

  best_params = [0,0,0,0,0,0]
  cur_best = Inf
  for p=100:150 #should be based on guess
    for ecc=0:0.5
      for tp=1:p
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
        global g_tp = tp
        global g_p = p
        global g_ecc = ecc

        fit = curve_fit(linear_rv, time, rv, p0);

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
    linear_rv(time, p)

For each time value calculate the model RV values for the given h, c, and v0
values in p.
"""
function linear_rv(time, p)
  rv = zeros(Float64, length(time))
  for i=1:length(time)
    ma = 2pi/g_p*(time[i]-g_tp); #calculate the mean anomaly
    f = true_anomaly(g_ecc,ma); #get the true anomaly
    rv[i] = p[1]*cos(f)+p[2]*sin(f)+p[3]
  end
  return rv
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
  result = get_optimal_rv_parameters(data,1);
  println("-----------------------------------")
  println("optimal values are:\nperiod=$(result[1]) days\necc=$(result[2])\ntp=$(result[3])")
end

#@stest test()
