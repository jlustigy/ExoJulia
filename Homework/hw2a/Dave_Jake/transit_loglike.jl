include("transit_model.jl")
include("../../hw2/Jake_Dave/utils.jl")
include("../../hw2/Jake_Dave/orbital_utils.jl")

function transit_loglike_observables(params::Vector)
    # params [ dF, tT, tF, P, tE ]
    
    """
    Observables
    -----------
    dF : Delta Flux = Rp/Rs
    tT : Transit Duration
    tF : Flat Duration 
    P  : Period 
    tE : Time of first transit
    """
    
    # Hard bounds
    if params[1] < 0.0 || sqrt(params[1]) >= 0.1
        return -Inf
    end
    if params[2] < 0.0 || params[2] > params[4]
        return -Inf
    end
    if params[3] < 0.0 || params[3] > params[2] || params[3] > params[4]
        return -Inf
    end
    if params[4] < (best_period - best_period*0.01) || (params[4] > best_period + best_period*0.01)
        return -Inf
    end
    if params[5] < 0.0 #|| (params[5] > best_toff + best_toff*0.1))
        return -Inf
    end
    
    model = transit_model_observables(time,params);
    
    return loglike(mean_flux, model, err);
end

function transit_loglike_observables_optim(params::Vector)
    return -transit_loglike_observables(params)
end