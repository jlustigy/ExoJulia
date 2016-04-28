using LsqFit
using PyPlot

"""
    circle_overlap(d, r_s, r_p)
Computer the area of overlap of two circles as a function of the distance
between them.

#Arguments
* `d`: the straight-line distance between the centers of the two bodies
* `r_s`: the radius of the star
* `r_p`: the radius of the planet, assumed that r_p < r_s
"""
function circle_overlap(d,r_s,r_p)
  if d <= r_s - r_p
    #the planet is completely within the star
    area = pi*r_p^2
  elseif d >= r_s + r_p
    #there is no overlap
    area = 0
  else
    phi_s = 2*(acos((r_s^2+d^2-r_p^2)/(2*r_s*d)))
    phi_p = 2*(acos((r_p^2+d^2-r_s^2)/(2*r_p*d)))
    area = 0.5*(phi_p*r_p^2-r_p^2*sin(phi_p)+phi_s*r_s^2-r_s^2*sin(phi_s))
  end

  return area
end


"""
    get_transit_parameters(time, flux, err)
Compute the time of first transit (t0), period (p), sky velocity (v) in units of
stellar radius/day, radius ratio (k), impact parameter (b), and unocculted flux
level (f0) that best match the provided data values.

The eccentricity is assumed zero and the planet is assumed to be fully in
transit (no grazing).

#Arguments
* `time::Array{Float64,1}`: an array containing N time points
* `flux::Array{Float64,1}`: an array containing N flux measurements
* `err::Array{Float64,1}`: an array containing N error measurements
* `trasnitTimeLimit::Float64=1.0`: the limit on transit duration

#Return
* `params:Array{Float64,1}`: an array of the best fit parameters
                             [p, t0, f0, k, b, v]
"""
function get_transit_parameters(time, flux, err; trasnitTimeLimit=1.0)


  p, t0, f0, k, b, v = initial_values(time, flux, trasnitTimeLimit)

  #params format: [p, t0, f0, k, b, v]
  params = [p, t0, f0, k, b, v] #initial guesses

  #println("initial estimate: p=$(p), t0=$(t0), f0=$(f0), k=$(k), b=$(b), v=$(v)")

  #closure function called by curve_fit
  function model_transit(time, params)
    p, t0, f0, k, b, v = params

    #Curve fit isn't that great and sometimes tries negative values for b...
    if b < 0.0
      params[5] = 0.0
      b = 0.0
    elseif b >= 1
      params[5] = 0.9
      b = 0.9
    end


    #calculate the expected flux for each time point using the given
    #parameters.
    model_flux = zeros(Float64,length(time))

    for i=1:length(time) #the signal should repeat every period (p) samples
      t = time[i]
      dist = ( b^2 + ( (1 - b^2)^0.5 - ( (t%p) - t0 )*v )^2 )^0.5 #the distance between the planet and star centers
      r_p = k #radius of planet in terms of stellar radii, k=r_p/r_s
      r_s = 1.0 #radius of star in units of stellar radii
      overlap = circle_overlap(dist, r_s, r_p) #percent of star covered


      model_flux[i] = f0*(1-overlap) #calculate the total flux
    end
    return model_flux
  end

  fit = curve_fit(model_transit, time, flux, 1./err.^2, params)
  p, t0, f0, k, b, v = fit.param

  #plot(time,flux, color="blue")
  #plot(time, model_transit(time,fit.param), color="red")
  #show()

  return [p, t0, f0, k, b, v]
end


"""
    initial_values(time, flux)
Guess values for the first transit (t0), period (p), sky velocity (v) in units of
stellar radius/day, radius ratio (k), impact parameter (b), and unocculted flux
level (f0) that best match the provided data values.

#Return
* `vals::Array{Float64,1}`: array containing [p, t0, f0, k, b, v]
"""
function initial_values(time, flux, trasnitTimeLimit)
  avg_f_val = mean(flux)
  min = minimum(flux)

  t0 = -1.0 #initialize to -1.0
  p_vals = []
  f0_avg = []
  b_vals = []
  k = sqrt( (avg_f_val-min)/avg_f_val )

  transit_lengths = []

  previous_was_transit = false

  for i=1:length(time)
    if abs(flux[i]-avg_f_val)*0.5 > abs(flux[i]-min)
      previous_was_transit = true
      #closer to the min, we're in a dip
      if -1.0 == t0
        #this is the first dip, set t0
        t0 = time[i]
      end

      #check if this is our first entry into this transit

      if length(p_vals)==0 || trasnitTimeLimit < time[i] - p_vals[end]
        push!(p_vals, time[i])

        #estimate b from the slope of the previous point to the current
        if i > 1
          temp_est = abs( ( flux[i] - flux[i-1] )/(time[i]-time[i-1]) )
          b_est = 1 - temp_est/(1 + temp_est)
          push!(b_vals, b_est)
        end
      end

    else
      #this is not a transit, push the value to our f0 avg array
      push!(f0_avg, flux[i])

      if previous_was_transit
        previous_was_transit = false
        push!(transit_lengths, time[i] - p_vals[end])
      end
    end
  end

  b = mean(b_vals)

  #get the average transit time
  T = mean(transit_lengths)
  v = 2*(1-b^2)^0.5/T

  #loop over the p vals to find the average p
  p_sum = 0.0
  for i=1:length(p_vals)-1
    p_sum += p_vals[i+1] - p_vals[i]
  end
  p = p_sum/(length(p_vals)-1)

  f0 = mean(f0_avg)

  b = mean(b_vals)

  return [p, t0, f0, k, b, v]
end


function read_data_file(filepath)
  data = readdlm(filepath)
  time = (data[:,1])
  flux = (data[:,2])
  err = (data[:,3])
  return (time, flux, err)
end



function test()
  time, flux, err = read_data_file("mystery_planet2.txt")
  p, t0, f0, k, b, v =  get_transit_parameters(time, flux, err)
  println("Got t0=$(t0), p=$(p), v=$(v), k=$(k), b=$(b), f0=$(f0)")
end


#@stest test()
