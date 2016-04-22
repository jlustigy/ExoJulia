"""
Seager+ 2003 Transit Model
Observables
-----------
dF : Delta Flux = Rp/Rs
tT : Transit Duration
tF : Flat Duration 
P  : Period 
tE : Time of first transit
"""

function transit_velocity(tT::Float64, RpRs::Float64, b::Float64)
    # Calculates transit velocity (units of Rs/day) as a function of:
    # tT : Transit Duration
    # RpRs : Radius of planet / Radius of Star
    # b : Impact Parameter
    return 2.0 * (RpRs + 1.0) * sqrt(1.0 - b^2) / tT
end

function impact_parameter(dF::Float64, tF::Float64, tT::Float64, P::Float64)
    # Calculates impact parameter in terms of observables
    sinF2 = sin(pi * tF / P)^2;
    sinT2 = sin(pi * tT / P)^2;
    return (((1.0 - sqrt(dF))^2 - (sinF2/sinT2) * (1.0 + sqrt(dF))^2)/(1.0 - sinF2/sinT2))^0.5
end

function aRs(dF::Float64, tT::Float64, P::Float64, b::Float64)
    # Calculates a/Rs in terms of observables
    sinT2 = sin(pi * tT / P)^2;
    return (((1.0 + sqrt(dF))^2 - b^2*(1.0 - sinT2))/(sinT2))^0.5
end

function stellar_density(dF::Float64, tT::Float64, P::Float64, b::Float64)
    # Calculates the stellar density (units of solar density) in terms of observables
    Gterm = 365.25^2 / 215^3; # [day^2 Msun Rsun^-3]
    sinT2 = sin(pi * tT / P)^2;
    return ((Gterm)/(P^2))*(((1.0 + sqrt(dF))^2 - b^2*(1.0 - sinT2))/(sinT2))^1.5
end

function impact_separation(b::Float64,t::Float64,v::Float64,P::Float64,toff::Float64)
    return sqrt(b.^2 + (v.*(mod(t,P)-toff)).^2)
end

function transit_model_observables(time::Array{Float64,1},params::Vector)
    #
    #Observables
    #-----------
    #dF : Delta Flux = Rp/Rs
    #tT : Transit Duration
    #tF : Flat Duration 
    #P  : Period 
    #tE : Time of first transit
    #
    
    # Write params (will want to remove this for less malloc)
    dF = params[1]
    tT = params[2]
    tF = params[3]
    P  = params[4]
    tE = params[5]
    
    # Calculate physical parameters from observables
    RpRs = sqrt(dF);
    b = impact_parameter(dF, tF, tT, P);
    vRs = transit_velocity(tT, RpRs, b);
    rhos = stellar_density(dF, tT, P, b);
    
    rel_flux = zeros(length(time))
    distance = 0.0
    for i=1:length(time)
        distance = impact_separation(b, time[i], vRs, P, tE);
        rel_flux[i] = relative_flux(distance,RpRs);
    end
    
    return rel_flux
end

function guess_params()
    
    print(time[1:5])
    return -1
end