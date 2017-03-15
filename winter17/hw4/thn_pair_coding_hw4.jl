using LsqFit
using PyPlot

include("ttv_wrapper_fixp3.jl")
include("ttv_nplanet.jl")

# Read in transit times
TT_1 = readdlm("ttv_planet1.txt")[:]
TT_2 = readdlm("ttv_planet2.txt")[:]
transit_times = [TT_1; TT_2]

# Functions
function ttv_model2(x, params)
   res = ttv_nplanet(2, 5, [length(TT_1), length(TT_2)], params)
   return [add_linear(res[1,:], params[3], params[2]); add_linear(res[2, :][1:length(TT_2)], params[8], params[7])]
end

function fit_ttv2(params)
   transit_time_fit   = curve_fit(ttv_model2, 0, [TT_1; TT_2], params)
   transit_time_model = ttv_model2(0, transit_time_fit.param)
   return transit_time_model, transit_time_fit
end

function add_linear(x, initial, period)
   return x + (collect(eachindex(x))-1)*period + initial
end

function fit_ttv3_fixed(params)
   res = curve_fit(ttv_wrapper_fixp3, 0, [TT_1; TT_2], params)
   return res
end

function log_likelihood(obs, mean, err)
    return sum(-0.5*((obs - mean).^2)./(err.^2))
end

# Find the period
period_1 = mean(TT_1[2:end]-TT_1[1:end-1])
period_2 = mean(TT_2[2:end]-TT_2[1:end-1])

# Initialize Parameters
# stellar mass ratio, planet period, center of transit time, eccentricity vectors (cos/sin)
# mass ratio, period, t_0, ecosomega, esinomega
init_1 = vec([2.e-5,period_1,TT_1[1],0.001,0.001])
init_2 = vec([3.e-5,period_2,TT_2[1],0.001,0.001])
init = [init_1;init_2]

# Fit two planet model
two_planet_model, two_planet_fit = fit_ttv2(init)
println(two_planet_fit.param)

# Fit three planet model (with two planet mode results)
n_periods = 10
periods   = collect(linspace(500,10000,n_periods))
n_phases  = 10
phases    = collect(linspace(0,2*pi))
likelihood = zeros(n_periods,n_phases)
ecc = 0.01
for i = 1:n_periods 
    global p3_cur = periods[i]
    for j = 1:n_phases 
        global ecos = ecc*cos(phases[j])
        global esin = ecc*sin(phases[j])
        params = [two_planet_fit.param; [1e-3, TT_1[1]]]
        three_planet_fit   = fit_ttv3_fixed(params)
        three_planet_model = ttv_wrapper_fixp3(0, params[1:end])  # TTVs corresp. to our fit
        likelihood[i,j]      = log_likelihood(three_planet_model,transit_times,1) 
    end
end

figure()
plot(periods,likelihood)

plot(periods[2:end],likelihood[2:end,:], ".")



