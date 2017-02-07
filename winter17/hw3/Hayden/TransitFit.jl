#=
HW 3 Problem 3
Author: Hayden Smotherman
OBJECTIVE: Simulate a planetary transit.
INPUTS: Transit Data File Name, System parameters initial guess
        (Params) [Period, Transit Duration, T0, Impact Factor, Planetary Radius (0<RP<1), Planet Flux (0<FP<1)]
OUTPUS: Flux as a fraction of stellar flux (NetFlux)
=#

using LsqFit
using PyPlot
include("TransitSolver.jl")

function TransitFit(FileName, InitGuess; Type=1)
  DATA    = readdlm(FileName);
  Time    = DATA[:,1];
  NetFlux = DATA[:,2];

  Period = InitGuess[1];
  #plot(Time, NetFlux,".")
  sTimePermutation = sortperm(Time%Period);
  sTime = (Time%Period)[sTimePermutation]; #Sorted time values
  sNetFlux = NetFlux[sTimePermutation]; #Sorted flux values

  if Type==1
    fit = curve_fit(TransitSolver, sTime, sNetFlux, InitGuess);
    figure()
    plot(sTime, sNetFlux, ".")
    plot(sTime, TransitSolver(sTime, fit.param))
    xlabel("Time (Days)")
    ylabel("Flux")
    title("Transit light curve fit")
    return fit.param
  elseif Type==2 #Experimental method to fit without simulating a curve
  #Initialize best values
    dtBest = 0.;
    qBest = 0.;
    AvgDiffBest = 0.;
    fMinBest = -1.;
    fMaxBest = -1.;

    qMin = max(2*(Time[2]-Time[1]),0.001);

    for q=qMin:0.001:0.1
      for dt=minimum(sTime):0.001:maximum(sTime)
        fMinTest = sNetFlux[dt.<sTime.<q*Period+dt];
        fMaxTest = sNetFlux[!(dt.<sTime.<q*Period+dt)];

        AvgDiff = (sum(fMaxTest)/maximum(size(fMaxTest)))^2 - (sum(fMinTest)/maximum(size(fMinTest)))^2;
        if AvgDiff > AvgDiffBest
          AvgDiffBest = AvgDiff;
          qBest = q;
          dtBest = dt;
          fMinBest = sum(fMinTest)/maximum(size(fMinTest));
          fMaxBest = sum(fMaxTest)/maximum(size(fMaxTest));
        end
      end
    end

    FluxDiff = fMaxBest - fMinBest;
    return qBest, dtBest, fMinBest, fMaxBest, sqrt(FluxDiff)
  else
    println("Please enter Type=1 or Type=2.")
  end
  #fit = curvefit(TransitSolver, Time, NetFlux, InitGuess);
end
