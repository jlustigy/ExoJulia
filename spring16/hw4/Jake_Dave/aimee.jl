# Affine Invariant McEmcE == aimee
# This module contains routines to run Affine-Invariant MCMC

function chisq(data::Array{Float64,2},model::Function,params::Vector)
    # Compute the simple chi-squared ln likelihood
    # For simple case, data = [x y err]
        
    return sum(((data[:,2] .- model(data,params)).^2.0)./(data[:,3].^2))
end

# End Function

function affine_inv_mcmc(nsteps::Int,params::Vector,model::Function,errors::Vector;
    walk_scale::Int=3,lnlike::Function=chisq,verbose::Bool=true)
    #
    # Run an affine-invariant markov chain
    # Foreman-Mackey et al. (2014) - 'emc' 'mcmc hammer'
    # Based on code written by Eric Agol 2016 (see fit_transits.jl)
    #
    # Parameters
    # ----------
    # nsteps : Int
    #    Number of MCMC steps for walkers to take
    # params : Vector
    #    Vector of inital guesses on parameters to optimize over
    # model : Function
    #    User-defined model used to compute ln likelihood
    #    Must be of the form: model(x::Array{Float64,1},y::Array{Float64,1},err::Array{Float64,1},params::Vector)
    # errors : Vector
    #   Vector of errors on initial params guess
    # walk_scale : Int
    #    Number of walkers per parameters.  3 is a good value
    # lnlike : Function
    #    Loglikelihood function.  Defaults to simple chi-square
    #
    # Returns
    # -------
    # par_mcmc : Array{Float64,3}
    #    Array containing history of all walkers for each parameter
    # ln_mcmc : Array{Float64,2}
    #    Array containing history of ln likelihood
    # nwalkers : Int
    #    Number of walkers used
    # nparam : Int
    #    Number of params fitted for
    
    nparam = length(params)
    
    # Want walkers to be a few time the number of params
    nwalkers = nparam * walk_scale
    
    # Set up arrays to hold the results:
    par_mcmc = zeros(nwalkers,nsteps,nparam)
    ln_mcmc = zeros(nwalkers,nsteps)
    
    # Initialize walkers:
    par_trial = params 
    
    # Estimate best ln like from input params
    ln_best = lnlike(data,model,par_trial)
    
    for j=1:nwalkers
        # Select from within uncertainties:
        ln_trial = 1.0e100
    
        # Only initiate models with reasonable ln like values:
        while ln_trial > (ln_best + 1000.0)
            par_trial = params + errors.*randn(nparam) 
            ln_trial = lnlike(data,model,par_trial)
        end
      
        ln_mcmc[j,1] = ln_trial
        par_mcmc[j,1,:] = par_trial
        
        if verbose
            println("Success: ",par_trial,ln_trial)
        end
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
        z = (rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2
        par_trial = vec(z*par_mcmc[j,i-1,:]+(1.0-z)*par_mcmc[ipartner,i-1,:])
    
        # Compute model & chi-square:  
            ln_trial = lnlike(data,model,par_trial)
    
        # Next, determine whether to accept this trial step:
        alp = z^(nparam-1)*exp(-0.5*(ln_trial - ln_mcmc[j,i-1]))
        if alp >= rand()
    
            # If step is accepted, add it to the chains!
            par_mcmc[j,i,:] = par_trial
            ln_mcmc[j,i,:] = ln_trial
            accept = accept + 1
       
        # If step is rejected, then copy last step:
        else
            par_mcmc[j,i,:] = par_mcmc[j,i-1,:]
            ln_mcmc[j,i,:] = ln_mcmc[j,i-1]
        end
      end
      
      if mod(i,5000) == 0
          frac_acc = accept/(5000*nwalkers)
          
          if verbose
              println("Number of steps: ",i," acceptance rate: ",frac_acc)
          end
          
          ascale = 1.0 + (frac_acc/0.25)*(ascale-1.0)
          accept = 0
      end
    end
    
    return par_mcmc, ln_mcmc, nwalkers, nparam
end

# End function

function est_burnin(par_mcmc::Array{Float64,3},nwalkers::Int,nparam::Int,nsteps::Int)
    # Determine time of burn-in by calculating first time median is crossed
    # Algorithm by Eric Agol 2016
    #
    # Parameters
    # ----------
    # par_mcmc : Array{Float64,3}
    #    Array containing history of all walkers for each parameter
    # nwalkers : Int
    #    Number of walkers used in MCMC
    # nparam : Int
    #    Number of params fit for in MCMC
    # nsteps : Int
    #    Number of MCMC steps
    #
    # Returns
    # -------
    # iburn : Int
    #    Index corresponding to the step where the burn in approximately ended
    
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
   
    return iburn
    
end

# End function