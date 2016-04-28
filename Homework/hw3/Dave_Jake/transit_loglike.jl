include("transit_model.jl")
include("../../hw2/Jake_Dave/utils.jl")
include("../../hw2/Jake_Dave/orbital_utils.jl")

function check_impact_params(params::Vector)
    # Returns false if params will crash the impact parameter calculation
    sinF2 = sin(pi * params[3] / params[4])^2;
    sinT2 = sin(pi * params[2] / params[4])^2;
    numerator = ((1.0 - sqrt(params[1]))^2 - (sinF2/sinT2) * (1.0 + sqrt(params[1]))^2)
    if numerator < 0.0
        return false
    else
        return true
    end 
end

function transit_loglike_observables(params::Vector)
    # params [ dF, tT, tF, P, tE ]
    
    # Hard bounds
    if params[1] < 0.0 || params[1] > 0.5
        return -Inf
    end
    if (params[2] < 0.0) || (params[2] > params[4])
        return -Inf
    end
    if (params[3] < 0.0) || (params[3] > params[2]) || (params[3] > params[4])
        return -Inf
    end
    if params[4] < (best_period - best_period*0.01) || (params[4] > best_period + best_period*0.01)
        return -Inf
    end
    if(params[5] < 0.0 || (params[5] > best_tE + best_tE*0.2))
        return -Inf
    end
    if ~check_impact_params(params)
        return -Inf
    end
        
    model = transit_model_observables(time,params);
    
    return loglike(mean_flux, model, err);
end

function transit_loglike_observables_optim(params::Vector)
    return -transit_loglike_observables(params)
end