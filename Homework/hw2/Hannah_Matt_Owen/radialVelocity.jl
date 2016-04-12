# Add ExoJulia/ to path
push!(LOAD_PATH, "../../../ExoJulia")
using ExoJulia

#=
TrueAnomaly()

Solve Kepler's equation for the true anomaly using Newton's method.

Input
e   The eccentricity of the orbit [0<e<1]
m   The mean anomaly in radians [0<m<2pi]

Return
f   The true anomaly [0<f<2pi] on success. -1 on failure
=#
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
"""
function get_optimal_rv_parameters(data)
  #TODO use Eric's algorithm to guess P

  time = data[:,1];
  rv = data[:,2];
  err = data[:,3];

  F = zeros(Float64,4,length(time)); #TODO, should 4 be passed in?

  for p=100:150 #should be based on guess
    for ecc=0:0.5
      for tp=1:p
        #initialize F
        for j=1:length(time)
          ma = 2pi/p*(time[j]-tp); #calculate the mean anomaly
          f = true_anomaly(ecc,ma); #get the true anomaly
          F[1,j] = cos(f);
          F[2,j] = sin(f);
          F[3,j] = 1;
          F[4,j] = time[j]-time[1];
        end

        W = diagm(1./err.^2); #the error matrix
        epsilon = inv(F * W * (F'));

        #calculate beta which gives {h, c, v0, d}
        beta = rv' * W * (F') * epsilon
        h = beta[1]
        c = beta[2]
        K = sqrt(h^2+c^2)

        #TODO check if this is the closest beta
      end
    end
  end
end






data =  readdlm("./mystery_planet.txt")
get_optimal_rv_parameters(data)
