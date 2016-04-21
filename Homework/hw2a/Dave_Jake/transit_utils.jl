include("../../hw2/Jake_Dave/utils.jl")
include("../../hw2/Jake_Dave/orbital_utils.jl")
include("../../hw2/Jake_Dave/rv.jl")

function overlap_area(sep::Float64,k::Float64)
    # r: star - planet disk center separation per rstar
    # k: rp/rstar
    
    R = 1.0
    # No transit
    if sep >= k + R
        return 0.0
    end
    
    # Full transit
    if sep <= R - k
        return (k*k)/(R*R)
    end
        
    # Partial transit
    A = k*k*acos((sep*sep + k*k - R*R)/(2.0*sep*k)) + R*R*acos((sep*sep + R*R - k*k)/(2.0*sep*R))
    A -= 0.5*sqrt((-sep + k + R)*(sep + k - R)*(sep - k + R)*(sep + k + R))
    
    return A/(pi*R*R)
end

function relative_flux(sep::Float64,k::Float64)
    # Returns the relative flux observed for a planet - star system
    # in or out of transit
    # r: star - planet disk center separation
    # k: rp/rstar
    # Works as long as units are same or nothing
    
    return 1.0 - overlap_area(sep,k)
end

function center_separation(t::Float64,t0::Float64,period::Float64,d::Float64,inc::Float64,ecc::Float64,varpi::Float64)
# Computes separation of centers for star - planet system assuming star center at origin
# From Kreidberg 2015
    f = f_from_t(period, ecc, t, t0)
    
    if (pi <= (abs(f) + varpi)) && ((abs(f) + varpi <= 2.0*pi))
        return 1.0e99
    end
    
    return (d*((1.0-ecc*ecc)/(1.0 + ecc*cos(f)))*sqrt(1.0-(sin(varpi+f)^2.)*(sin(inc)^2.)))
end

function transit_loglike(params::Vector)
    # params = [rp/rs, d, ecc, varpi, inc, period] where d = a/rstar
    
    # Hard bounds
    if params[1] < 0.0 || params[1] >= 0.1
        return Inf
    end
    if params[2] <= 0.0
        return Inf
    end
    if params[3] < 0.0 || params[3] >= 1.0
        return Inf
    end
    if(params[4] < 0.0) || (params[4] > 2.0*pi)
        return Inf
    end
    if (params[5] > (100.0*pi/180.0)) || (params[5] < (80.0*pi/180.0))
        return Inf
    end
    if(params[6] <= 0.0)
        return Inf
    end
    # Probability of transit must be between (0,1]
    if (1./params[2])*(params[1] + 1.)/(1.0-params[3]^2.) <= 0.0
        return Inf
    end
    if (1./params[2])*(params[1] + 1.)/(1.0-params[3]^2.) > 1.0
        return Inf
    end
    
    # Note: - because optim looks for MINIMUM!
    model = transit_model(time,params);
    
    return -loglike(mean_flux, model, err);
end

function transit_model(time::Array{Float64,1},params::Vector)
    # params = [rp/rs, d, ecc, varpi, inc, period] where d = a/rstar
    k = params[1]
    d = params[2]
    ecc = params[3]
    varpi = params[4] + pi/2.0 # rads
    inc = params[5]
    per = params[6]
    
    fi = 3.*pi/2 - varpi
    tp = per*sqrt(1.-ecc*ecc)/(2.*pi)*(ecc*sin(fi)/(1.+ecc*cos(fi))-2./sqrt(1.-ecc*ecc)*atan2(sqrt(1.-ecc*ecc)*tan(fi/2.),1.+ecc))
    
    input_data = [time mean_flux err]
    data_fold = copy(input_data)
    data_fold[:,1] = mod(time[1] - time[1], per)
    data_fold = fastsortrows(data_fold,[1]);
    
    rel_flux = zeros(length(time))
    distance = 0.0
    for i=1:length(time)
        distance = center_separation(time[i],time[1]+tp,per,d,inc,ecc,varpi)
        rel_flux[i] = relative_flux(distance,k)
    end
    
    return rel_flux
end