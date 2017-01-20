#HW 1b Problem 2
#===================
FUNCTION NAME: RV_Fit
INPUTS:   FileName, InitialGuess [K, w, P, tp, gamma, ecc]
OUTPUTS:  ?
===================#
using LsqFit
using PyPlot
include("RVSolver.jl")
include("NewtonKeplerSolver.jl")

function RVFit(FileName, InitialGuess)
  #Read the data in from the text file
  DATA = readdlm(FileName)
  t   = DATA[:,1]; #Time
  RV  = DATA[:,2]; #Radial velocity
  eRV = DATA[:,3]; #Error in radial velocity

  BestError = 1.e300; #Initialize the BestError variable
  Bestt = 0; #Initialize the Best Time variable

  #Fold the data based on the guessed period
  for j = 1:100
  #Iteratively select smaller and smaller sampling ranges
  PeriodRange = linspace((1-0.9/j)*InitialGuess[3],(1+0.9/j)*
    InitialGuess[3],10);
  for i=1:10 #Sample the given period space
    t = DATA[:,1] - minimum(DATA[:,1]); #Reset t with the first data point being zero
    while maximum(t)>PeriodRange[i] #Fold time based on the guessed period
      t[t.>PeriodRange[i]]-=PeriodRange[i]
    end
    Order = sortperm(t);
    Error = sum((RV[Order[2:end]]-RV[Order[1:end-1]]).^2);
    if Error < BestError
      BestError = Error;
      InitialGuess[3] = PeriodRange[i];
      Bestt = t;
    end
  end
  end

  fit = curve_fit(RVSolver, Bestt, RV, 1./(eRV.^2), InitialGuess);
  plot(Bestt,RV,".")
  PlotDomain = linspace(minimum(Bestt),maximum(Bestt),1000);
  plot(PlotDomain, RVSolver(PlotDomain,fit.param))
  ylabel("Radial Velocity m/s")
  xlabel("Time")
  title("Radial Velocity best-fit")
  return fit.param
end
