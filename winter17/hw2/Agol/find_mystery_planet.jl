# This Julia script carries out a fit to an RV dataset
# 'mystery_planet.txt' Eric Agol April 2016

# This code carries out the following:
# 1). A Lomb-Scargle periodogram of the data (this amounts
#     to RV fitting with a circular orbit);
# 2). A gridded search over eccentricity & time of periastron;
# 3). A final optimization of the single-planet parameters given the best-fit
#     chi-square from the gridded search.
# 4). A search for a second planet in the residuals of the first.
# 5). An estimate of the false-alarm probability of the second planet.

include("regress.jl")
# This is an RV model for a single planet:
include("rv_model_one.jl")
# This is the linearized RV model:
include("rv_model_lin.jl")
# These routines were used in earlier versions of the code:
#using Optim
#include("rv_model.jl")
#include("chisquare.jl")
#include("kepler_solve.jl")
# Add appropriate path for ExoJulia:
push!(LOAD_PATH,"../../../ExoJulia")
using ExoJulia
using LsqFit
# PyPlot is used for plotting below (scatter & plot):
using PyPlot

# General circular periodogram:
data = readdlm("mystery_planet1.txt")
ntimes = 174
# Number of period points to search over:
nperiod = 100
period=collect(linspace(111,113,nperiod))
# Define chi-square for a circular fit and eccentric fit:
chisq_circ = zeros(nperiod)
chisq_ecc = zeros(nperiod)
# Define the times, radial-velocities, and uncertainties.
# Make these global variables so that they are accessible from
# the model routines:
global times = vec(data[:,1])
global rv = vec(data[:,2])
global sigrv = vec(data[:,3])
#sigrv = vec(data[:,3])
#times = vec(data[:,1])
# Set up a model array for the circular fits: 
xx = zeros(3,ntimes)
# Set the third basis function to unity (this computes the
# velocity zero-point in the circular fits):
xx[3,:]=1.0
#xx[4,:]=times-mean(times)
# Define the weights:
weight = 1./sigrv.^2
# Do a grid search over eccentricity and time of periastron:
necc=100
nt0=100
# Define the initial chi-square to be large:
chibest = 1e100
# Define an initial parameter vector for the circular fits:
pbest = zeros(3)
for i=1:nperiod
# First, carry out a circular orbit fit (this is a Lomb-Scargle periodogram):
  xx[1,:]=sin(2pi/period[i].*times)
  xx[2,:]=cos(2pi/period[i].*times)
# Use the 'regress.jl' routine, which carries out a linear regression
# of a dataset against basis functions that are speicified:
  coeff = regress(xx,rv,sigrv)
# From the best-fit coefficients, compute the circular model:
  model = vec(coeff[1]*xx[1,:]+coeff[2]*xx[2,:]+coeff[3]*xx[3,:])
# Compute the chi-square for the circular model:
  chisq_circ[i]=sum((rv-model).^2./sigrv.^2)
# Next, do eccentric fits.
# Loop over eccentricity & time of periastron:
  chisq_ecc[i]=chisq_circ[i]
  for jecc=1:necc
    for jt0=1:nt0
      param = [period[i],jt0*period[i]/nt0,.8+(jecc*0.199)/necc]
# The routine rv_model_lin carries out the linearized fitting of Wright & Howard (2009):
      rvmod,coefflin = rv_model_lin(times,param,rv,sigrv)
# Compute the chi-square from the best fit model:
      chisq = sum(((rv-rvmod)./sigrv).^2)
# If the chi-square improves upon the best fit so far, then save it!
      if (chisq < chisq_ecc[i])
        chisq_ecc[i] = chisq
        if(chisq < chibest)
# Save the best chi-square:
          chibest = chisq
# Save the best-fit model parameters, both the non-linear (param) and linear (coefflin):
          pbest = [param;coefflin]
        end
      end
    end
  end
#  fit = curve_fit(rv_model, times, rv, weight, param) 
#  fit = curve_fit(rv_model_one, times, rv, weight, param) 
#  fit = curve_fit(rv_model, times, rv, weight, param) 
#  chisq_ecc[i] = sum(fit.resid.^2./sigrv.^2)
#  chisq_ecc[i] = sum(((rv-rvmod)./sigrv).^2)
#  println(period[i],' ',chisq_ecc[i])
end
clf()
plot(period,log10(chisq_circ))
plot(period,log10(chisq_ecc))
xlabel("Period (days)")
ylabel("Log(chi-square)")

println("Best-fit period: ",pbest)
# Wait until user hits Enter/Return:
println("Hit return to continue")
read(STDIN,Char)
# Clear the plotting window:
clf()

# Now that we have done a gridded search over eccentricity,
# we're going to take the best fit model and further optimize it.
# We'll use the Levenberg-Marquardt algorithm, as implemented
# in curve_fit, in the package LsqFit.
# Initialize the parameters with the best parameters from the gridded search:
param = pbest
fitbest = curve_fit(rv_model_one, times, rv, weight, pbest)
# The Optim package also has a BFGS optimizer, which we tried at some point.
#fit = optimize(chisquare, pbest, method = :l_bfgs, iterations = 21)
pbest = fitbest.param
model = rv_model_one(times,pbest)
println("Chi-square: ",sum((rv-model).^2.*weight))
# Set up a finer grid of times over which to compute the best-fit model
# for plotting:
tgrid = linspace(0,pbest[1],1000)
rvmod = rv_model_one(collect(tgrid),pbest)
scatter(mod(times,pbest[1]),rv)
plot(tgrid,rvmod)
xlabel("Time (days)")
ylabel("Radial velocity (m/s)")
# Compute the 1-sigma (68.26% confidence) interval for the parameters:
errors=estimate_errors(fitbest,0.6826)
# Define parameter names (using equation 4 in Wright & Howard):
pname = ["Period [P]               ",
         "Time of pericenter [t_p] ",
         "Eccentricity [e]         ",
         "Cosine coefficient [h]   ",
         "Sine coefficient [c]     ",
         "RV zero-point [v_0]      "]
nparam = 6
# Print out the best-fit paramters (with formatted print macro):
for i=1:nparam
   @printf "%s %7.3f +/- %5.3f\n" pname[i] pbest[i] errors[i]
end
# Now plot the residuals:
println("Hit return to continue")
read(STDIN,Char)
clf()
errorbar(mod(times,pbest[1]),fitbest.resid,yerr=sigrv,fmt="o")
xlabel("Phase [days]")
ylabel("Residuals [m/s]")
println("Hit return to continue")
read(STDIN,Char)
clf()
# What is the scatter in the residuals?
println("Scatter: ",std(fitbest.resid)," Scatter in units of error bars: ",std(fitbest.resid./sigrv))
# If the scatter is larger than the uncertainties, this indicates residual stellar jitter!
# Inflate the uncertainties by the residual scatter:
sigrv = sigrv*std(fitbest.resid./sigrv)

# Now, fit the residuals with a periodogram to see if there is an
# additional planet:

rv_resid = fitbest.resid
#nperiod2 = 1000000
nperiod2 = 1000
period2=collect(linspace(10,1000,nperiod2))
chisq_circ2 = zeros(nperiod2);
chibest2 = 1e15
period_best2  = 0.
for i=1:nperiod2
  xx[1,:]=sin(2pi/period2[i].*times)
  xx[2,:]=cos(2pi/period2[i].*times)
  coeff = regress(xx,rv_resid,sigrv)
  model = vec(coeff[1]*xx[1,:]+coeff[2]*xx[2,:]+coeff[3]*xx[3,:])
  chisq_circ2[i]=sum((rv_resid-model).^2./sigrv.^2)
  if(chisq_circ2[i] < chibest2)
    chibest2 = chisq_circ2[i]
    period_best2 = period2[i]
  end
end
# Now, carry out Monte-carlo:
ntrial = 1000
chisq_trial = zeros(nperiod2,ntrial)
for j=1:ntrial
  irandrv = sortperm(rand(ntimes))
  rv_trial = rv_resid[irandrv]
  sigrv_trial = sigrv[irandrv]
  for i=1:nperiod2
    xx[1,:]=sin(2pi/period2[i].*times)
    xx[2,:]=cos(2pi/period2[i].*times)
    coeff = regress(xx,rv_trial,sigrv_trial)
    model = vec(coeff[1]*xx[1,:]+coeff[2]*xx[2,:]+coeff[3]*xx[3,:])
    chisq_trial[i,j]=sum((rv_trial-model).^2./sigrv_trial.^2)
  end
end
# Now compute 0.001% confidence limit (this corresponds to a ~4.4-sigma outlier):
# This could be improved by optimizing simultaneous fits to both planets
# at each trial period, as well as carrying out eccentric fits for the
# second planet.
chisort = sort(vec(chisq_trial))
chisq_sig=chisort[10]
plot(log10(period2),chisq_circ2)
plot(log10(period2),chisq_sig+zeros(nperiod2))
xlabel("Log(Period [days])")
ylabel("Chi-square")
#
# If the chisq_sig is smaller than the smallest chi-square
# of the periodogram, then the potential planet is insignificant at that level.
#
println("Hit return to continue")
read(STDIN,Char)
clf()
# Plot the best-fit second planet:
scatter(mod(times,period_best2),rv_resid)
xlabel("Orbital phase [days]")
ylabel("RV residual [m/s]")
