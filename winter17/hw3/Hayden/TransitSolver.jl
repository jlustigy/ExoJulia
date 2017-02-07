#=
HW 3 Problem 2
Author: Hayden Smotherman
OBJECTIVE: Simulate a planetary transit.
INPUTS: Time (t), System parameters
        (Params) [Period, Transit Duration, T0, Impact Factor, Planetary Radius (0<RP<1), Planet Flux (0<FP<1)]
OUTPUS: Flux as a fraction of stellar flux (NetFlux)
=#

include("CircleOverlap.jl")

function TransitSolver(t, Params)
  Period       = Params[1];
  TransTime    = Params[2];
  T0           = Params[3];
  b            = Params[4];
  PlanetRadius = Params[5];
  PlanetFlux   = Params[6];

  t = t-minimum(t); #Normalize t to zero
  t -= T0; #Offset t by T0

  NetFlux = ones(size(t)); #Initialize net stellar flux and planetary flux
  for i=1:maximum(size(t))
    if (0.)<=t[i]%Period<=(TransTime)
      x = ((t[i])%Period)*(2*sqrt((1+PlanetRadius).^2-b.^2))/(TransTime)-sqrt((1+PlanetRadius).^2-b.^2);
  #x = linspace(-sqrt((1+PlanetRadius).^2-b.^2),sqrt((1+PlanetRadius).^2-b.^2),maximum(size(NetFlux[(0.).<=t%Period.<=(TransTime)])));
      r = sqrt(b.^2+x.^2);
      NetFlux[i] = CircleOverlap(PlanetRadius,r);
    elseif (Period/2)<=t[i]%Period<=(Period/2+TransTime)
      NetFlux[i] = 1-PlanetFlux;
    end
  end

  return NetFlux
end
