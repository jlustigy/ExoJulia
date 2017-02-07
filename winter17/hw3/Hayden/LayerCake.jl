#=
HW 3 Problem 4
Author: Hayden Smotherman
OBJECTIVE: Simulate a planetary transit with limb darkening
INPUTS: Time
        System parameters  (Params) [Period, Transit Duration, T0, Impact Factor, Planetary Radius (0<RP<1), Planet Flux (0<FP<1), NumLayers]
        Limb Darkening Function
OUTPUS: Flux as a fraction of stellar flux (NetFlux)
=#

using PyPlot
include("TransitSolver.jl")

function LayerCake(t, Params, LimbDark)
  Period, Transit, T0, b, Rp, Fp, NumLayer::Int = Params;

  dy = abs((LimbDark(1)-LimbDark(0))/NumLayer); #How high should each layer be
  Mu = linspace(0,1,1001);
  Theta = acos(Mu);
  dR = sin(Theta);
  LayerRadii = zeros(NumLayer);

  Checkpoint = 0.;
  j=0;

  for i=1:maximum(size(dR))
    if (LimbDark(Mu[i])-Checkpoint)>=dy
      Checkpoint = LimbDark(Mu[i]);
      j+=1;
      #if j<=NumLayer
      LayerRadii[j]=dR[i];
      #end
    end
  end
  LayerRadii = LayerRadii[LayerRadii.!=0];
  NetFlux = zeros(size(t));
  for i=1:maximum(size(LayerRadii))
    NewTransit = Transit*sqrt(((LayerRadii[i]+Rp)^2-b^2)/((1-Rp)^2-b^2))
    NewT0 = T0+(Transit-NewTransit)/2;
    NetFlux += LayerRadii[i]^2*TransitSolver(t,[Period, NewTransit, NewT0, b, Rp, Fp]);
  end
  NetFlux = NetFlux/(maximum(NetFlux)) #Normalize to 1
  return NetFlux
end
