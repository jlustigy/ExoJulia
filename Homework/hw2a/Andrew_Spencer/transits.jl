# Transits Module
# ExoJulia Package
# ASTR598
# PI: Eric Agol
# Written by: Andrew Lincowski & Spencer Wallace
# University of Washington
# Time-stamp: <transits.jl on Saturday, 23 April, 2016 at 12:45:09 PDT (linc)>

# Need to change velocity of planet to duration of transit
# Need to allow for impact parameter outside of star (still partial transit)

module Transit

# Required packages
using LsqFit

# General functions
function fastsortrows(B::AbstractMatrix,cols::Array; kws...)
  """
  Solution by: abhishekmalali (gihub)
  See: https://github.com/JuliaLang/julia/issues/9832
  """
       for i = 1:length(cols)
           if i == 1
               p =sortperm(B[:,cols[i]]; kws...);
               B = B[p,:];
           else
               i0_old = 0;
               i1_old = 0;
               i0_new = 0;
               i1_new = 0;
               for j = 1:size(B,1)-1
                   if B[j,cols[1:i-1]] == B[j+1,cols[1:i-1]] && i0_old == i0_new
                       i0_new = j;
                   elseif B[j,cols[1:i-1]] != B[j+1,cols[1:i-1]] && i0_old != i0_new && i1_new == i1_old
                       i1_new = j;
                   elseif i0_old != i0_new && j == size(B,1)-1
                       i1_new = j+1;
                   end
                   if i0_new != i0_old && i1_new != i1_old
                       p = sortperm(B[i0_new:i1_new,cols[i]]; kws...);
                       B[i0_new:i1_new,:] = B[i0_new:i1_new,:][p,:];
                       i0_old = i0_new;
                       i1_old = i1_new;
                   end
               end
           end
       end
    return B
end

########## TRANSIT FUNCTIONS ###############

# 1. Compute overlap of two circles as function of their center-center separation (use law of cosines!)

# 2. Compute transit/secondary eclipse of body (no limb darkening)


function pi_crust(delta::Float64,r::Float64)
    # Calculates overlaps of two circles (the sums of their respective pie crusts)
    # r = sqrt(a*a+b*b): separation of circles center-to-center
    # delta = (Rpl/R*)^2
    # r >= 0
    # Current limitation is impact parameter larger than R_st is not allowed

    if(r >= (1 + sqrt(delta)))
        # outside of transit
        return 0.0
    elseif(r <= (1 - sqrt(delta)))
        # Fully in transit
        return 1.0
    else
        # ingress / egress
        th_st = 2*acos((1+ r*r - delta)/(2*r))
        th_pl = 2*acos((r*r+delta-1)/(2*r*sqrt(delta)))
        alpha = (0.5*th_st - 0.5*sin(th_st) + 0.5*th_pl*delta - 0.5*delta*sin(th_pl))/(pi*delta)
        return alpha
    end
end



# 3. Fit model to mystery_planet2.txt.
    #    Calculate:
    # planet period P
    # transit depth K
    # impact parameter b
    # duration of transit T
    # density of star


function transit(delta,b,t,t0,dur,period,base_flux)
    # delta = (Rpl/R*)^2
    # b = impact parameter
    # t = time
    # t0 = time of first contact
    # dur = duration of transit
    # v = velocity
    # period = period of planet orbit
    
    # prevent bad parameters
    if (delta <= 0. || delta > 1.)
        #println("bad delta")
        return typemax(Float64)
    elseif (b < 0. || b > (1.))
        #println("bad b")
        return typemax(Float64)
    end
    
    v = 2*(sqrt(delta) + 1)*sqrt(1 - b*b)/dur
    t = mod(t,period)
    a = v*(t-t0)-sqrt(1-b^2)-sqrt(delta)
    r = sqrt(a*a+b*b)

    return base_flux*(1 - pi_crust(delta,r)*delta)
end


function transit_model(x,p)
        
    # Input model function to fit transit curve
    # x: independent variable (time)
    # p: dependent variables (fit parameters)
    
    delta = p[1]
    b = p[2]
    t0 = p[3]
    dur = p[4]
    period = p[5]
    base_flux = p[6]
    
    tr_arr = Array(Float64,length(x))
    for i in 1:length(x)
        tr_arr[i] = transit(delta,b,x[i],t0,dur,period,base_flux)
    end

    return tr_arr
end
    
############ EXECUTE MODEL ###################

    # Data import
    pldata = readdlm("mystery_planet2.txt")
    time_data = pldata[:,1]
    flux_data = pldata[:,2]
    err_data = pldata[:,3]


    #Probe useful parameter space for period of planet orbit
    minp = 10.0
    maxp = 25.0
    np = 1000.0
    period = linspace(minp,maxp,np)
    sum = Array(Float64,length(period))

    phase = Array(Float64,length(time_data))

    for (j,P) in enumerate(period)
        sum[j] = 0.0
        
        #Sort by phase, given period
        phase = mod(time_data,P)
        phase_data_arr = [phase flux_data] #combine arrays
        phase_sorted = fastsortrows(phase_data_arr, [1]) #sort by phase
        for i in 2:length(time_data)
            sum[j] += (phase_sorted[i,2]-phase_sorted[i-1,2])*(phase_sorted[i,2]-phase_sorted[i-1,2])
        end
    end

    # Estimated period fit (days)
    min_per=period[indmin(sum)]

    # Initial values for curve_fit
    # Curve_fit is VERY sensitive to initial guess

    delta = (maximum(flux_data)-minimum(flux_data))/maximum(flux_data)
    b = 0.1
    t0 = 0.1
    dur = 1.0
    flux = mean(flux_data)

    p = [delta,b,t0,dur,min_per,flux] #period P from phase folding estimate
    println("Guess parameters: p = [delta,b,t0,dur,min_per,flux]")
    for (idx, val) in enumerate(p)
        println(val)
    end

    # Run fitting routine for eccentricity & time of periastron
    fit = curve_fit(transit_model,time_data,flux_data,err_data,p)
    println("Fit parameters:")
    println(fit.param)

end
