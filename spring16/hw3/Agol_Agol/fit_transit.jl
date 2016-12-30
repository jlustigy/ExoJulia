# This script fits transiting planet data, mystery_planet2.txt
# Eric Agol April 2016
# The following occur:
# 1). A transit model function is defined, transit_model, using
#     the quadratic limb-darkening model, occultquad.jl
# 2). The data are read in, and an initial fit is computed.
# 3). A Markov chain is run.
# 4). Various parameters of the fit are computed, along with uncertainties.
using PyPlot
using LsqFit
push!(LOAD_PATH,"../../../ExoJulia")
using ExoJulia

#include("occultquad.jl")
data=readdlm("mystery_planet2.txt")

time = vec(data[:,1])
flux = vec(data[:,2])
eflx = vec(data[:,3])

function transit_model(time,p)
# Models transits with a circular orbit
# p[1] = period
# p[2] = t0
# p[3] = tdur
# p[4] = b
# p[5] = k
# p[6] = q1 # limb-darkening parameters due to David Kipping (2015)
# p[7] = q2 # 0 < q_1 < 1; 0 < q_2 < 1
# p[8] = f0
  period = p[1] ; t0 = p[2]; tdur = p[3]; b = p[4]
  k = p[5]; q1 = p[6] ; q2 = p[7] ; f0 = p[8]
  if (-1 < b < 1) && (0 < q1 < 1) && (0 < q2 < 1)
    v = 2*sqrt(1-b^2)/tdur
    u1 = sqrt(q1*q2)
    u2 = sqrt(q1)*(1.0-2q2)
  else
    tdur = 0.0
    v = 0.0
  end
  fmod = ones(length(time)).*f0
  for i=1:length(time)
    tmod = mod((time[i] - t0 + period/2),period)-period/2
    if abs(tmod) < (2*tdur)
      z = sqrt(b^2 + (v*tmod)^2)
      fmod[i] = f0*ExoJulia.Transit.occultquad(z,u1,u2,k)
    end
  end
  return fmod 
end

#p = [12.2,0.8,5.0,0.0,0.03,0.2,0.2,median(flux)]
p = [12.2,0.8,0.34,0.3,0.04,0.5,0.5,median(flux)]
weight = 1./eflx.^2
nparam = 8
fit = curve_fit(transit_model,time,flux,weight,p)
println("Best-fit parameters: ",fit.param)
flx_best = transit_model(time,fit.param)
chi_best= sum(((flux-flx_best)./eflx).^2)
pname = ["Period (d)","t0 (d)","Duration (d)","Impact parameter","Planet/star radius ratio","q1","q2","f0"]
#errors = estimate_errors(fit,0.6826)
errors = [0.001,0.001,0.01,0.2,0.001,0.3,0.3,1.0]
println("Chi-square: ",chi_best)
period = fit.param[1]
phase = mod((time-fit.param[2]+period/2),period)-period/2
scatter(phase,flux)
b=fit.param[4]
tdur=fit.param[3] 
t0=fit.param[2] 
#tdur = 2*sqrt(1-b^2)/v
tmod = collect(tdur*linspace(-2,2,1000)) + t0
fmod = transit_model(tmod,fit.param)
plot(tmod-t0,fmod,color="red")

# Now, run an affine-invariant markov chain:
# Foreman-Mackey et al. (2014) - 'emc' 'mcmc hammer'

nwalkers = nparam * 3
nsteps = 10000
#nsteps = 100
# Set up arrays to hold the results:
par_mcmc = zeros(nwalkers,nsteps,nparam)
chi_mcmc = zeros(nwalkers,nsteps)
# Initialize walkers:
par_trial = fit.param
for j=1:nwalkers
# Select from within uncertainties:
  chi_trial = 1e100
# Only initiate models with reasonable chi-square values:
  while chi_trial > (chi_best + 1000)
    par_trial = fit.param + errors.*randn(nparam)
    model =transit_model(time,par_trial)
    chi_trial = sum(((flux-model)./eflx).^2)
  end
  chi_mcmc[j,1]=chi_trial
  par_mcmc[j,1,:]=par_trial
  println("Success: ",par_trial,chi_trial)
end
# Initialize scale length & acceptance counter:
ascale = 2.0
accept = 0
# Next, loop over steps in markov chain:
for i=2:nsteps
  for j=1:nwalkers
    ipartner = j
# Choose another walker to 'breed' a trial step with:
    while ipartner == j
      ipartner = ceil(Int,rand()*nwalkers)
    end
# Now, choose a new set of parameters for the trial step:
    z=(rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2
    par_trial=vec(z*par_mcmc[j,i-1,:]+(1.0-z)*par_mcmc[ipartner,i-1,:])
# Compute model & chi-square:    
    model_trial =transit_model(time,par_trial)
    chi_trial=sum(((flux-model_trial)./eflx).^2)
# Next, determine whether to accept this trial step:
    alp = z^(nparam-1)*exp(-0.5*(chi_trial - chi_mcmc[j,i-1]))
    if alp >= rand()
# If step is accepted, add it to the chains!
      par_mcmc[j,i,:] = par_trial
      chi_mcmc[j,i,:] = chi_trial
      accept = accept + 1
    else
# If step is rejected, then copy last step:
      par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
      chi_mcmc[j,i,:] = chi_mcmc[j,i-1]
    end
  end
  if mod(i,1000) == 0
    frac_acc = accept/(1000*nwalkers)
    println("Number of steps: ",i," acceptance rate: ",frac_acc)
    ascale = 1.0 + (frac_acc/0.25)*(ascale-1.0)
    accept = 0
  end
end

# Now, determine time of burn-in by calculating first time median is crossed:
iburn = 0
for i=1:nparam
  med_param=median(par_mcmc[1:nwalkers,1:nsteps,i])
  for j=1:nwalkers
    istep=2
    while (par_mcmc[j,istep,i] > med_param) == (par_mcmc[j,istep-1,i] > med_param) && (istep < nsteps)
      istep=istep+1
    end
    if istep >= iburn
      iburn = istep
    end
  end
end

println("Burn-in ends: ",iburn)

# Also need trace plots & compute autocorrelation length [ ]

# Compute the various parameters:

pavg = mean(vec(par_mcmc[1:nwalkers,iburn:nsteps,1]));
psig = std(vec(par_mcmc[1:nwalkers,iburn:nsteps,1]));
@printf "Period: %8.4f +/- %8.4f\n" pavg psig

b=vec(par_mcmc[1:nwalkers,iburn:nsteps,4]);
bavg = mean(b)
bsig = std(b)
@printf "Impact par: %5.2f +/- %5.2f\n" bavg bsig

davg = mean(vec(par_mcmc[1:nwalkers,iburn:nsteps,5].^2*100));
dsig = std(vec(par_mcmc[1:nwalkers,iburn:nsteps,5].^2*100));
@printf "Depth (pct):   %5.2f +/- %5.2f\n" davg dsig

tdur = vec(par_mcmc[1:nwalkers,iburn:nsteps,3]);
tavg = mean(tdur*24)
tsig = std(tdur*24)
@printf "Duration (hr):   %5.2f +/- %5.2f\n" tavg tsig

using CGS
rho = 3/pi^2./(tdur.*24*3600).^3.*(1-b.^2).*(period.*24*3600)/GRAV;
#rho = 3/pi^2./(tdur.*24*3600).^3.*(period.*24*3600)/GRAV;
ravg = mean(rho)
rsig = std(rho)
@printf "Stellar density (g/cc):   %5.2f +/- %5.2f\n" ravg rsig

for i=1:nparam
  for j=1:nwalkers
    plot(vec(par_mcmc[j,1:nsteps,i]))
  end
  xlabel("MCMC step")
  ylabel(pname[i])
  println("Hit return to continue")
  read(STDIN,Char)
  clf()
end

for i=2:nparam
  for j=1:i-1
    scatter(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]),vec(par_mcmc[1:nwalkers,iburn:nsteps,j]))
    xlabel(pname[i])
    ylabel(pname[j])
    println("Hit return to continue")
    read(STDIN,Char)
    clf()
  end
end
