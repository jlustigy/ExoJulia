using PyPlot
using LsqFit
using PyCall

include("../compute_ttv.jl")
include("../ttv_succinct.jl")

"""
# Calculates the total TTV at each time point for a
# set of planets.
#
# Input:
#   transit_time: A vector, with each element containing
#                 an array of transit times for one planet
#           pl: A vector of Planet_plane_hk objects
#
# Returns:
#   A vector, with each element containing the TTV offset
#   at every time point for the given planets
"""

function hw4()

function calc_ttv_total(transit_time::Vector, pl::Vector)
    N = length(pl)
    jmax = 15

    ttv_total = []
    for i=1:N
        append!(ttv_total, zeros(length(transit_time[i])))
    end

    for i=1:N-1
        for j=(i+1):N
            ttv_i = zeros(length(transit_time[i]))
            ttv_j = zeros(length(transit_time[j]))
            TTVFaster.compute_ttv!(jmax, pl[i], pl[j], transit_time[i], transit_time[j], ttv_i, ttv_j)
            ttv_total[i] .+= ttv_i
            ttv_total[j] .+= ttv_j
        end
    end

    return ttv_total
end

# Load transit timing data
transit_time1 = readdlm("../ttv_planet1.txt")[:,1]
transit_time2 = readdlm("../ttv_planet2.txt")[:,1]

# Find approximate period from transit data
sum1 = 0.0
sum1a = 0.0
sum2 = 0.0
sum2a = 0.0
for i in 1:(length(transit_time1)-1)
    sum1 += transit_time1[i+1] - transit_time1[i]
end

Per1 = sum1 / (length(transit_time1)-1)

for i in 1:(length(transit_time1)-1)
    sum1a += (transit_time1[i+1] - transit_time1[i] - Per1)^2
end

for i in 1:(length(transit_time2)-1)
    sum2 += transit_time2[i+1] - transit_time2[i]
end


Per2 = sum2 / (length(transit_time2)-1)

for i in 1:(length(transit_time2)-1)
    sum2a += (transit_time2[i+1] - transit_time2[i] - Per2)^2
end

Stdev1 = sqrt(sum1a/(length(transit_time1)-1))
Stdev2 = sqrt(sum2a/(length(transit_time2)-1))

function timing_model_2p(x, p)
    m_rat1 = p[1]
    period1 = p[2]
    t01 = p[3]
    ecosw1 = p[4]
    esinw1 = p[5]
    
    m_rat2 = p[6]
    period2 = p[7]
    t02 = p[8]
    ecosw2 = p[9]
    esinw2 = p[10]
    
      # Inner planets must have smaller periods
      if ((period1 > 0.) && (period1 < period2))
        nt1 = length(transit_time1)   
        nt2 = length(transit_time2)   
        t1 = collect(linspace(t01, t01 + period1*(nt1-1), nt1))
        t2 = collect(linspace(t02, t02 + period2*(nt2-1), nt2))
        
        pl = Array(TTVFaster.Planet_plane_hk{Float64}, 2)
        pl[1] = TTVFaster.Planet_plane_hk{Float64}(m_rat1, period1, t01, ecosw1, esinw1)
        pl[2] = TTVFaster.Planet_plane_hk{Float64}(m_rat2, period2, t02, ecosw2, esinw2)

        ttv_total = calc_ttv_total(Vector{Float64}[transit_time1, transit_time2], pl)

        t1 .+= ttv_total[1]
        t2 .+= ttv_total[2]
        
        retarr = collect([t1;t2])
        return retarr
    end

    retarr[:] = typemax(Float64)
    return retarr
end

# 2 planet lsq fit
p = [1e-3, Per1, transit_time1[1], 0.05, 0.05,
     1e-3, Per2, transit_time2[1], 0.05, 0.05]

y = collect([transit_time1; transit_time2])
x = collect(linspace(1,length(y),length(y)))
weight = Array(Float64,length(y))
errorbars = Array(Float64,length(y))
weight[:] = (30.0/86400.0)^-2.
errorbars[:] = 30.0/86400.0
fit = curve_fit(timing_model_2p, x, y , errorbars, p)
println("Two-planet fit: ")
print(fit.param)
#read(STDIN,Char)

global Period3 = 4332.59
function timing_model_3p_fixed_p3(x, p)
    m_rat1 = p[1]
    period1 = p[2]
    t01 = p[3]
    ecosw1 = p[4]
    esinw1 = p[5]
    
    m_rat2 = p[6]
    period2 = p[7]
    t02 = p[8]
    ecosw2 = p[9]
    esinw2 = p[10]
    
    m_rat3 = p[11]
    t03 = p[12]
    ecosw3 = p[13]
    esinw3 = p[14]
    nt1 = length(transit_time1)
    nt2 = length(transit_time2)
    t1 = collect(linspace(t01, t01 + period1*(nt1-1), nt1))
    t2 = collect(linspace(t02, t02 + period2*(nt2-1), nt2))

      # Inner planets must have smaller periods
      if ((period1 > 0.) && (period1 < period2) && (period2 < Period3) &&
        # Restrict ecc vector to be between -1 and 1 (and not 0)
        (ecosw1 != 0.) && (esinw1 != 0.) && (ecosw2 != 0.) && (esinw2 != 0.) && (ecosw3 != 0.) && (esinw3 != 0.) &&
        (ecosw1 > -1.) && (esinw1 > -1.) && (ecosw2 > -1.) && (esinw2 > -1.) && (ecosw3 > -1.) && (esinw3 > -1.) &&
        (ecosw1 < 1.) && (esinw1 < 1.) && (ecosw2 < 1.) && (esinw2 < 1.) && (ecosw3 < 1.) && (esinw3 < 1.) &&
        (m_rat1 > 0.) && (m_rat2 > 0.) && (m_rat3 > 0.))
        #println("Good parameter!")
        t3 = collect(linspace(0., Period3*9., 10))
        
        pl = Array(TTVFaster.Planet_plane_hk{Float64}, 3)
        pl[1] = TTVFaster.Planet_plane_hk{Float64}(m_rat1, period1, t01, ecosw1, esinw1)
        pl[2] = TTVFaster.Planet_plane_hk{Float64}(m_rat2, period2, t02, ecosw2, esinw2)
        pl[3] = TTVFaster.Planet_plane_hk{Float64}(m_rat3, Period3, t03, ecosw3, esinw3)

        ttv_total = calc_ttv_total(Vector{Float64}[transit_time1, transit_time2, t3], pl)
        
        t1 .+= ttv_total[1]
        t2 .+= ttv_total[2]

        retarr = collect([t1;t2])
        return retarr
    end

    #retarr = zeros(length(transit_time1)+length(transit_time2))
    return collect([t1;t2])
    #retarr[:] = typemax(Float64)
    #return retarr
end

# 3 planet lsq fit

#p = [1e-6, Per1, transit_time1[1], 0.05, 0.05,
#    1e-6, Per2, transit_time2[2], 0.05, 0.05,
#    1e-2,      transit_time1[1], 0.05, 0.05]
param_2planet = fit.param
nphase = 10
phase = Period3 * linspace(0,1,nphase)
param3_best = zeros(15)
chibest = 1e100
y = collect([transit_time1;transit_time2])
x = collect(linspace(1,length(y),length(y)))
model3_best = zeros(length(y))
# Loop over phase of planet 3
for i=1:nphase
  p = [param_2planet; [1e-3,  phase[i], 0.01, 0.01]]

  weight = Array(Float64,length(y))
  errorbars = Array(Float64,length(y))
  weight[:] = (30.0/86400.0)^(-2.)
  errorbars[:] = 30.0/86400.0
  fit = curve_fit(timing_model_3p_fixed_p3, x, y, errorbars, p)
  model = timing_model_3p_fixed_p3(x, fit.param)
  chi = sum(((model - y)./errorbars).^2)
  if chi < chibest
     chibest = chi
     param3_best = fit.param
     model3_best = model
     println("Three-planet fit: ")
     print(chi)
     print(fit.param)
#     read(STDIN,Char)
  end
end


function mcmc_fit!(model, x::Vector, y::Vector, y_error::Vector, lsq_fit_params::Vector, 
                   param_errors::Vector, chi_best::Float64, param_mean::Vector, param_std::Vector, )  
    nparam = length(lsq_fit_params)
    nwalkers = nparam * 3
    nsteps = 50000
    # Set up arrays to hold the results:
    par_mcmc = zeros(nwalkers,nsteps,nparam)
    chi_mcmc = zeros(nwalkers,nsteps)
    # Initialize walkers:
    par_trial = lsq_fit_params
    for j=1:nwalkers
        # Select from within uncertainties:
        chi_trial = 1e100
        # Only initiate models with reasonable chi-square values:
        while chi_trial > (chi_best + 1000)
            par_trial = lsq_fit_params + param_errors.*randn(nparam)
            y_model = model(x, par_trial)
            chi_trial = sum(((y-y_model)./y_error).^2)
        end
        chi_mcmc[j,1]=chi_trial
        par_mcmc[j,1,:]=par_trial
        println("Success: ",par_trial,chi_trial)
    end
    
    # Initialize scale length & acceptance counter:
    ascale = 1.5
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
            y_model_trial = model(x,par_trial)
            chi_trial = sum(((y-y_model_trial)./y_error).^2)
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
    
    for i=1:nparam
        param_mean[i] = mean(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]))
        param_std[i] = std(vec(par_mcmc[1:nwalkers,iburn:nsteps,i]))
    end
    return par_mcmc
end

function timing_model_3p(x, p)
    m_rat1 = p[1]
    period1 = p[2]
    t01 = p[3]
    ecosw1 = p[4]
    esinw1 = p[5]
    
    m_rat2 = p[6]
    period2 = p[7]
    t02 = p[8]
    ecosw2 = p[9]
    esinw2 = p[10]
    
    m_rat3 = p[11]
    period3 = p[12]
    t03 = p[13]
    ecosw3 = p[14]
    esinw3 = p[15]
    
    # Inner planets must have smaller periods
    if ((period1 > 0.) && (period1 < period2) && (period2 < period3))
    
        t1 = collect(linspace(t01, t01 + period1*length(transit_time1), length(transit_time1)))
        t2 = collect(linspace(t02, t02 + period2*length(transit_time2), length(transit_time2)))
        t3 = collect(linspace(0., Period3*9., 10))
        
        pl = Array(TTVFaster.Planet_plane_hk{Float64}, 3)
        pl[1] = TTVFaster.Planet_plane_hk{Float64}(m_rat1, period1, t01, ecosw1, esinw1)
        pl[2] = TTVFaster.Planet_plane_hk{Float64}(m_rat2, period2, t02, ecosw2, esinw2)
        pl[3] = TTVFaster.Planet_plane_hk{Float64}(m_rat3, period3, t03, ecosw3, esinw3)

        ttv_total = calc_ttv_total(Vector{Float64}[transit_time1, transit_time2, t3], pl)
        
        t1 .+= ttv_total[1]
        t2 .+= ttv_total[2]

        retarr = collect([t1; t2])
        return retarr
    end

    retarr = zeros(length(transit_time1)+length(transit_time2))
    retarr[:] = typemax(Float64)
    return retarr
end


#p = [1e-6, Per1, transit_time1[1], 0.05, 0.05,
#    1e-6, Per2, transit_time2[2], 0.05, 0.05,
#    1e-2, Period3, transit_time1[1], 0.05, 0.05]#

p = [param3_best[1:11];Period3;param3_best[12:14]]
#                m   per  t0   ecosw esinw
param_errors = [1e-7, 0.1, 0.1, 0.05, 0.05, 
                1e-7, 0.1, 0.1, 0.05, 0.05, 
                1e-4, 10.0, 10.0, 0.05, 0.05]
chi_best = 1e20##

param_mean = zeros(15)
param_std = zeros(15)

y = collect([transit_time1; transit_time2])
x = collect(linspace(1, length(y), length(y)))

y_error = Array(Float64,length(y))
y_error[:] = 30.0/86400.0

mcmc1 = mcmc_fit!(timing_model_3p, x, y, y_error, p, param_errors, chi_best, param_mean, param_std)


# julia throws an error if I dont place some code between comment blocks
#a = 4

# output histogrammed mcmc results to a text file for later use




@pyimport numpy

num_datapts = size(mcmc1)[1]*size(mcmc1)[2]
for j=1:15
  outfile=open("$j.txt", "w")
  data = reshape(mcmc1[:, :, j], 1, num_datapts)
  n, bins = numpy.histogram(data, normed=true, bins=60)
  for i=1:length(n)
    a = n[i]
    b = bins[i]
    write(outfile, "$a $b\n")
  end
  close(outfile)
end

return
end
