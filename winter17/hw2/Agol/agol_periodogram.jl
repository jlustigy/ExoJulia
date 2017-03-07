using PyPlot

pmin = 1.0
pmax = 1000.0
nperiod = 500000

#function agol_periodogram(pmin,pmax,nperiod)
# Folds RV on trial periods from pmin to pmax
# with a grid of length nperiod.  At trial each
# period, it computes the time modulo the
# period, takes the sum of squares of the
# difference between adjacent RV values, and
# mimimizes this value, returning the period
# with the minimum scatter.

data = readdlm("mystery_planet1.txt")

period=collect(linspace(pmin,pmax,nperiod))
scat_adjacent = zeros(nperiod)
rv = data[:,2]
scat_best= 1e100
pbest = 0.0
for i=1:nperiod
  phase = mod(data[:,1],period[i])
  isort = sortperm(phase)
  diff = rv[isort[2:59]]-rv[isort[1:58]]
  scat_adjacent[i] = sum(diff.*diff)
  if scat_adjacent[i] < scat_best
    scat_best = scat_adjacent[i]
    pbest = period[i]
  end
end
scatter(period,log10(scat_adjacent))
# Return the period that minimizes the scatter
# between adjacent data points:
#return pbest
#end
