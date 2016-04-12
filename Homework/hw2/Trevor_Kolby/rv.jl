#= rv.jl
functions to solve for the best fit period given RV data.
includes functions from hw1 (Kepler's solver).
=#

push!(LOAD_PATH, "../ExoJulia/")
using ExoJulia
using LsqFit

# comment out the next line if using notebook version, if using stest then you need this:
data = readdlm(string(pwd(), "/hw2/Trevor_Kolby/mystery_planet.txt"),Float64) ;

times = data[:,1]
RVs = data[:,2]
RV_errs = data[:,3] 

#@stest fit_RV(times,RVs,RV_errs)


M(t::Array,t_p::Float64,P::Float64) = (2*pi/P).*(t.-t_p)
M(t::Float64,t_p::Float64,P::Float64) = (2*pi/P)*(t-t_p)

EtoF(E::Float64,ecc::Float64) = 2.*atan((((1.+ecc)/(1.-ecc))^(1./2.))*tan(E/2.))
EtoF(E::Array,ecc::Float64) = 2.*atan((((1.+ecc)/(1.-ecc))^(1./2.)).*tan(E./2.))

K(h::Float64,c::Float64) = sqrt(h*h+c*c)

pomega(h::Float64,c::Float64) = atand(-c/h)

gamma(v0::Float64,K::Float64,ecc::Float64,pomega::Float64) = v0-K*ecc*cosd(pomega)


function P_guesser(x_data::Array,y_data::Array)
    #minimum period to look for is distance between data points. Max is size of data. Step is 0.5
    P_array::Array{Float64} = collect(mean(diff(x_data)):0.5:maximum(x_data)-minimum(y_data))
    sq_array::Array{Float64} = [sum(diff(y_data[sortperm(mod(x_data,P))]).^2.0) for P=P_array]
    return P_array[findmin(sq_array)[2]]
end


function kepler_solve!(Ms::Array{Float64},ecc::Float64)
    results = Float64[]
    for em in Ms
        push!(results,ExoJulia.Orbit.kepler_solve!(em, ecc))
    end
    return results
end


function hcv0(f::Array{Float64},y_data::Array{Float64},y_error::Array{Float64})
    W = diagm(1.0./(y_error.^2.0))
    F = hcat(cos(f),sin(f),[1.0 for x in f])'
    epsilon = inv(F*W*F')
    return y_data'*W*F'*epsilon
end


function fit_RV(x_data::Array{Float64},y_data::Array{Float64},y_err::Array{Float64})
    
    function v_rad(t::Array{Float64},params::Array{Float64}) #params are e,t_p,P
        if 0 <= params[1] < 1
            Es::Array{Float64} = kepler_solve!(M(t,params[2],params[3]),params[1])
            f::Array{Float64} = EtoF(Es,params[1])
            h,c,v0 = hcv0(f,y_data,y_err)
            return h.*cos(f)+c*sin(f).+v0
        else
            return Inf
        end
    end
    
    function get_params(fit_obj::LsqFit.LsqFitResult{Float64}) #returns e,t_p,P,K,pomega,gamma
        params::Array = fit_obj.param
        Es::Array = kepler_solve!(M(x_data,params[2],params[3]), params[1])
        f::Array = EtoF(Es,params[1])
        h,c,v0 = hcv0(f,y_data,y_err)
        return params[1],params[2],params[3],K(h,c),pomega(h,c),gamma(v0,K(h,c),params[1], pomega(h,c))
    end
        
    return get_params(curve_fit(v_rad,x_data,y_data,(1.0./y_err.^2.0),[0.1,x_data[1],P_guesser(x_data,y_data)]))
    
end



function eff_func(t,t0::Float64,r0::Float64,sma::Float64,ecc::Float64,t_p::Float64,P::Float64)
    Ms = M(t,t_p,P)
    M0 = M(t0,t_p,P)
    E0 = kepler_solve!(M0,ecc)
    E = kepler_solve!(Ms,ecc)
    return (sma/r0).*(cos(E.-E0).+1)
end



function gee_func(t,t0::Float64,r0::Float64,sma::Float64,ecc::Float64,t_p::Float64,P::Float64)
    Ms = M(t,t_p,P)
    M0 = M(t0,t_p,P)
    E0 = kepler_solve!(M0,ecc)
    E = kepler_solve!(Ms,ecc)
    return (t.-t0).+(P/(2.0*pi)).*(sin(E.-E0)-(E.-E0))
end



function plot_RV(times::Array{Float64},rvs::Array{Float64},period::Float64)
    phase = mod(times, period)/period;
    plot(phase, rvs, "o", ms=5)
    xlim(-.05,1.05)
    ylim(minimum(rvs)-abs(.1*minimum(rvs)),maximum(rvs)+abs(.1*maximum(rvs)))
    xlabel("Phase")
    ylabel("RV")
end
