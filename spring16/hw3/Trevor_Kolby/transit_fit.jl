using PyPlot
using LsqFit

#Angle between line connecting center of planet to where edge of planet meets edge of star
#and the line connecting planet center to star center
alpha(r::Float64,k::Float64,l::Float64) = acos((r^2.0 + k^2.0 - l^2.0)/(2.0*k*r))   

#Angle between line connecting centers of planet and star and line connecting center of star to
#where edge of planet meets edge of star
beta2(r::Float64,k::Float64,l::Float64) = acos((r^2.0 + l^2.0 - k^2.0)/(2.0*l*r)) 

#projected distance between center of planet and star 
function r(tau::Array{Float64},b::Float64,k::Float64,l::Float64,T::Float64)
    v = 2.0*sqrt((k+l)^2.0-b^2.0)/T
    return sqrt(b^2.0 .+ (v.*tau).^2.0)
end

#Time before/after center of transit, modulo the period
tau(t::Array{Float64},P::Float64,t_n::Float64) = mod(t.-t_n+(P/2),P)-(P/2.0)

#Overlapping area of two circles of radius k (planet) and l (star)
function overlap(rs::Array{Float64},k::Float64,l::Float64)
    result = Float64[]
    for r in rs
        if r >= l+k
            push!(result,0.0)
        elseif r <= l-k
            push!(result,pi*k^2.0)
        else
            push!(result,(k^2.0)*(alpha(r,k,l)-cos(alpha(r,k,l))*sin(alpha(r,k,l)))+(l^2.0)*(beta2(r,k,l)-cos(beta2(r,k,l))*sin(beta2(r,k,l))))
        end
    end
    return result
end

#Computes the transit depth as a function of time
function transit(ts::Array{Float64},params::Array{Float64})
    #params is [F_s,k,P,T,t_n,b]
    taus = tau(ts,params[3],params[5])
#    rs = r(taus,params[7],params[2],params[3],params[5])
    rs = r(taus,params[6],params[2],1.0,params[4])
#    overlaps = overlap(rs,params[2],params[3])
    overlaps = overlap(rs,params[2],1.0)
#    return params[1]*(1.0.-overlaps./(pi*params[3]^2.0))
    return params[1]*(1.0.-overlaps/pi)
end

#Can be used for secondary occultation if desired (requires an extra parameter)
function occultation(ts::Array{Float64},params::Array{Float64})
    #params is [F_p,k,l,P,T,t_n,b]
    taus = tau(ts-(params[3]/2.0),params[3],params[5])
#    rs = r(taus,params[7],params[2],params[3],params[5])
    rs = r(taus,params[6],params[2],1.0,params[4])
    overlaps = overlap(rs,params[2],1.0)
    return params[1]*(1.0.-overlaps./(pi*params[2]^2.0))
end

# First guess at period, accurate to within 1 day
function transit_P_guess_init(x_data::Array{Float64},y_data::Array{Float64})
    arr = hcat(x_data,y_data)
    lowest = sortrows(arr, by=x->(x[2])) #sort all points by increasing flux
    mins = sort(unique(round(lowest[1:50,1]))) #take the dimmest 50 points, sort by time rounded
    return maximum(diff(mins)) 
end

#First guess for the time of transit, assume its the time of minimum flux
transit_tn_guess_init(x_data::Array{Float64},y_data::Array{Float64}) = x_data[findmin(y_data)[2]]

#Now let's search around our initial period guess, phase fold the data, and search for the period 
#such that we minimize the mean value of the measured flux around the center of the transit
function transit_P_guess(x_data::Array{Float64},y_data::Array{Float64},P_guess::Float64,t_n::Float64)
    #period guess will be to the nearest day. Let's make an array of 100 points between the nearest two days
    P_arr = linspace(P_guess-1,P_guess+1,100)
    mean_arr = Float64[]
    for P in P_arr
        phases = (tau(x_data,P,t_n)/P)+0.5
        push!(mean_arr,mean(y_data[0.495.<=phases.<=0.505]))
    end
    return P_arr[mean_arr.==minimum(mean_arr)][1]
end

# Now that we have a much better period, take the median phase of the 75 faintest points, find 
# offset from the center of transit, and adjust the center of transit time accordingly
function transit_tn_guess(x_data::Array{Float64},y_data::Array{Float64},P::Float64,t_n::Float64)
    phases = (tau(x_data,P,t_n)/P)+0.5
    y_low = sort(y_data)[1:75]
    return t_n+(median(phases[findin(y_low,y_data)])-0.5)*P
end

#determines a best guess for star flux
transit_Fs_guess(y_data::Array{Float64}) = median(round(y_data))

#Assume l~R_sun. Calculate a rough transit depth, which gives ratio of radii squared.
#This will give an underestimate for when b>l-k (planet disk never fully overlaps the star)
function transit_k_guess(y_data::Array{Float64},F_s::Float64)
#    l = 6.96e10
    l = 1.0
    ysort = sort(y_data)
    depth_guess = 1.0 - median(ysort[1:100])/F_s
    return sqrt(depth_guess)*l
end

#Guess the transit duration by starting at the center of transit, and moving outward til the
#mean of the three nearest points starts to stabilize around the best fit star flux
function transit_T_guess(x_data::Array{Float64},y_data::Array{Float64},t_i::Float64,F_s::Float64)
    i_init = find(x_data .== t_i)[1]
    i = i_init
    av = 0.0
    while (F_s - av) >= 0.0015*F_s
        av = mean([y_data[i-1];y_data[i];y_data[i+1]])
        i += 1
    end
    return 2.0*(x_data[i]-x_data[i_init])
end

#Fit function just wraps up all of the guesses and then fits, returns a fit object.
#to get params, simply do result = fit_transit(x_data,y_data,y_err) then result.param
function fit_transit(x_data::Array{Float64},y_data::Array{Float64},y_err::Array{Float64})
#    l_guess = 6.96e10
    l_guess = 1.0
    P_guess_init = transit_P_guess_init(x_data,y_data)
    tn_guess_init = transit_tn_guess_init(x_data,y_data)
    P_guess = transit_P_guess(x_data,y_data,P_guess_init,tn_guess_init)
    tn_guess = transit_tn_guess(x_data,y_data,P_guess,tn_guess_init)
    F_guess = transit_Fs_guess(y_data)
    k_guess = transit_k_guess(y_data,F_guess)
    T_guess = transit_T_guess(x_data,y_data,tn_guess_init,F_guess)
    offset = (tn_guess-tn_guess_init)/tn_guess
    p0 = [F_guess,k_guess,P_guess,T_guess,tn_guess_init+abs(offset),0.0]
    fit = curve_fit(transit,x_data,y_data,(1.0./y_err).^2.0,p0)
    return fit
end

#Read in data, separate into arrays
data = readdlm("mystery_planet2.txt",Float64) 
times = data[:,1]
flux0 = data[:,2]
flux = flux0/maximum(flux0)
flux_err = data[:,3] ;

result = fit_transit(times,flux0,flux_err) ;
#@stest result = fit_transit(times,flux0,flux_err) 

function plot_transit(x_data::Array{Float64},y_data::Array{Float64},fit::LsqFit.LsqFitResult{Float64})
    #Outputs a phase folded normalized plot of light curve (points) with best fit (line)
    params = fit.param
    P = params[3]
    t_n = params[5]
    T = params[4]
    F0 = params[1]
    norm_flux = y_data/F0
    y_fit = transit(x_data,params)/F0
    phases = (tau(x_data,P,t_n)/P)+0.5
    
    v = sortperm(phases)   
    scatter(phases[v], norm_flux[v], color = "k", alpha=0.2)
    plot(phases[v],y_fit[v],"-r") 
    xlim(0.5-(1.5*T/P),0.5+(1.5*T/P));
end

function density(x_data::Array{Float64},fit::LsqFit.LsqFitResult{Float64})
    Msun = 1.9891e33  # g
    Rsun = 6.955e10   # cm
    
    P = fit.param[3]+0.19
    T = fit.param[4]
    F0 = fit.param[1]
    y_fit = transit(x_data,fit.param)/F0
    b = fit.param[6]*Rsun
    dflux = maximum(1.0-y_fit)
    
    # factor to put it in units of sun
    factor = (365.25^2.0/215.^3.0)*(Msun / Rsun^3.0)
    rho = (factor/P^2.0) * ((1+sqrt(dflux))^2.0 -(b^2.0*(1-(sinpi(T/P))^2.0))/(sinpi(T/P)^2.0))^1.5
    return rho
end

println("Stellar density: ",density(times,result))

plot_transit(times,flux0,result)
