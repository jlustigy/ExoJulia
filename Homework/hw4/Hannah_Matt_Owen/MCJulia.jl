function aimc(model, time, params, perrors, yvals, yerrs)
  # Now, run an affine-invariant markov chain:
  # Foreman-Mackey et al. (2014) - 'emc' 'mcmc hammer'


  nsteps = 10000
  nparam = length(params)
  nwalkers = nparam * 3
  #nsteps = 100
  # Set up arrays to hold the results:
  par_mcmc = zeros(nwalkers,nsteps,nparam)
  chi_mcmc = zeros(nwalkers,nsteps)

  y_best = model(time,params)
  chi_best= sum(((yvals-y_best)./yerrs).^2)

  # Initialize walkers:
  par_trial = params
  for j=1:nwalkers
  # Select from within uncertainties:
    chi_trial = 1e100
  # Only initiate models with reasonable chi-square values:
    count = 0
    while chi_trial > (chi_best + 1000) && count < 100
      count += 1
      par_trial = params + perrors.*randn(nparam)
      model_output = model(time,par_trial)
      chi_trial = sum(((yvals-model_output)./yerrs).^2)
    end
    chi_mcmc[j,1]=chi_trial
    par_mcmc[j,1,:]=par_trial
    #println("Success: ",par_trial,chi_trial)
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
      model_trial =model(time,par_trial)
      chi_trial=sum(((yvals-model_trial)./yerrs).^2)
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

  # Compute the various parameters:

  results = []

  for i=1:nparam
    pavg = mean(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]));
    psig = std(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]));
    push!(results,(pavg,psig))
  end
  println(results)
  return results
end
