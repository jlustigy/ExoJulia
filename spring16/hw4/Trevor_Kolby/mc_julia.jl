using PyPlot

type walker
    params::Array{Float64}
    chi2::Float64
end

type walker_record
    walkers::Array{walker}
    nsteps::Int64
end
walker_record(ws::Array{walker}) = walker_record(ws,length(ws)) 

type walker_array
    recs::Array{walker_record}
    nwalkers::Int8
end
walker_array(recs::Array{walker_record}) = walker_array(recs,length(recs)) 

chi2(f::Function,params::Array{Float64},data::Array{Float64,2}) = sum(((f(data[:,1],params).-data[:,2])./data[:,3]).^2) 

type ploglike 
    model::Function
    nparams::Int64
    data::Array{Float64,2}
    chi::Function
end
ploglike(func::Function,nparams::Int64,data::Array{Float64,2}) = ploglike(func,nparams,data,chi2)

import Base.push!
function push!(rec::walker_record,w::walker)
    push!(rec.walkers,w)
    rec.nsteps += 1
end
function push!(arr::walker_array,rec::walker_record)
    push!(arr.recs,rec)
    arr.nwalkers += 1
end

function init_walkers(plog::ploglike,p0::Array{Float64},p_err::Array{Float64})
    nparam = plog.nparams
    nwalker = nparam*3.0
    chi_best = plog.chi(plog.model,p0,plog.data)
    w_array = walker_array(walker_record[])
    for j=1:nwalker
        chi_trial = 1.0e100
        par_trial = p0
        while chi_trial > (chi_best + 1000)
            par_trial = p0 + p_err.*randn(nparam)
            chi_trial = plog.chi(plog.model,par_trial,plog.data)
        end
        w = walker(par_trial,chi_trial)
        push!(w_array,walker_record([w]))
    end
    return w_array
end

function walk(w::walker,partner::walker,ascale::Float64,plog::ploglike)
    #Takes a walker and a partner, picks a new trial step based on parameters
    z=(rand()*(sqrt(ascale)-1.0/sqrt(ascale))+1.0/sqrt(ascale))^2
    par_trial = z.*w.params.+(1.0-z).*partner.params
    chi_like = plog.chi
    modl = plog.model
    # Compute chi-square:    
    chi_trial=chi_like(modl,par_trial,plog.data)
    alp = z^(plog.nparams-1.0)*exp(-0.5*(chi_trial - w.chi2)) #Assumes chi2 as well, needs to be generalized
    if alp >= rand()
        # If step is accepted, add it to the chains!
        accept = 1
        return (walker(par_trial,chi_trial),accept)
    else
        # Otherwise, copy the current walker
        accept = 0
        return (w,accept)
    end
end

function mcmc(plog::ploglike,p0::Array{Float64},p_err::Array{Float64},nstep::Int64)
    w_array = init_walkers(plog,p0,p_err)
    nparam = plog.nparams
    nwalker = w_array.nwalkers
    ascale = 2.0
    accept = 0
    for i=2:nstep
        #Like j = nwalkers, but now we can access the walkers directly!
        for (j,rec) in enumerate(w_array.recs)
            current_w = rec.walkers[end]
            ipartner = j
            while ipartner == j
                ipartner = ceil(Int,rand()*nwalker)
            end
            partner_w = w_array.recs[ipartner].walkers[end]
            # Now, make a trial walker
            (w,a) = walk(current_w,partner_w,ascale,plog)
            push!(rec,w)
            accept += a
        end
        if i%1000 == 0
            frac_acc = accept/(1000*nwalker)
            println("Number of steps: $i, acceptance rate: $frac_acc")
            ascale = 1.0 + (frac_acc/0.25)*(ascale-1.0)
            accept = 0
        end
    end
    return w_array
end

get_history(w_array::walker_array,nwalk::Int64,p1::Int64,nstart::Int64) = [w.params[p1] for w in w_array.recs[nwalk].walkers[nstart:end]]

function get_median(w_array::walker_array,p1::Int64,nstart::Int64)
    plist = Float64[]
    for i=1:w_array.nwalkers
        append!(plist,get_history(w_array,i,p1,nstart))
    end
    return median(plist)
end

function get_mean(w_array::walker_array,p1::Int64,nstart::Int64)
    plist = Float64[]
    for i=1:w_array.nwalkers
        append!(plist,get_history(w_array,i,p1,nstart))
    end
    return mean(plist)
end

function get_std(w_array::walker_array,p1::Int64,nstart::Int64)
    plist = Float64[]
    for i=1:w_array.nwalkers
        append!(plist,get_history(w_array,i,p1,nstart))
    end
    return std(plist)
end

function calc_iburn(w_array::walker_array)
    # Now, determine time of burn-in by calculating first time median is crossed:
    iburn = 0
    nparam = length(w_array.recs[1].walkers[1].params)
    nwalker = w_array.nwalkers
    nstep = w_array.recs[1].nsteps
    for i=1:nparam
        med_param=get_median(w_array,i,1)
        for rec in w_array.recs
            istep=2
            while (rec.walkers[istep].params[i] > med_param) == (rec.walkers[istep-1].params[i] > med_param) && (istep < nstep)
              istep=istep+1
            end
            if istep >= iburn
              iburn = istep
            end
        end
    end
    return iburn
end

function get_params(w_array::walker_array,iburn::Int64)
    nparam = length(w_array.recs[1].walkers[1].params)
    params = Float64[get_mean(w_array,i,iburn) for i=1:nparam]
    stds = Float64[get_std(w_array,i,iburn) for i=1:nparam]
    return params,stds
end

function plot_walkers(w_array::walker_array,p1::Int64,p2::Int64,nstart::Int64)
    pones = Float64[]
    ptwos = Float64[]
    for w=1:w_array.nwalkers
        append!(pones,get_history(w_array,w,p1,nstart))
        append!(ptwos,get_history(w_array,w,p2,nstart))
    end
    scatter(pones,ptwos)
end

function plot_trace(w_array::walker_array,p1::Int64,nstart::Int64)
    for w=1:w_array.nwalkers
        eyes = collect(1:w_array.recs[w].nsteps)
        vals = get_history(w_array,w,p1,nstart)
        plot(eyes,vals)
    end
    show()
end

function drop_the_mic(plog::ploglike,p0::Array{Float64},p_err::Array{Float64},nstep::Int64)
    final_walkers = mcmc(plog,p0,p_err,nstep)
    iburn = calc_iburn(final_walkers)
    vals,errs = get_params(final_walkers,iburn)
    for (val,err) in zip(vals,errs)
        println("$val +/- $err")
    end
    for i=2:plog.nparams
        for j=1:i-1
            plot_walkers(final_walkers,i,j,iburn)
            show()
            clf()
        end
    end
    for i=1:plog.nparams
        plot_trace(final_walkers,i,1)
        show()
        clf()
    end
end