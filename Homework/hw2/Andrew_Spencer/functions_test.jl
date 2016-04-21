# Kepler's equation
# ecc = eccentricity
# E = eccentric anomaly
# M = mean motion  "$HW" == "./$1"
#@stest time_rv()

### Define necessary functions ###

### PLANET ORBITAL PERIOD ###


using LsqFit



function time_rv()

    # Required packages
    #using LsqFit

    #Data import
    pldata = readdlm("./hw2/Andrew_Spencer/mystery_planet.txt")
    time_data = pldata[:,1]
    global RV_data = pldata[:,2]
    err_data = pldata[:,3]

    #Probe useful parameter space for period of planet orbit
    period = linspace(1,3500,100000)
    sum = Array(Real,length(period))

    for (j,P) in enumerate(period)
        sum[j] = 0.0
        
        #Sort by phase, given period
        phase = Array(Real,length(time_data))
        phase = mod(time_data,P)
        phase_data_arr = [phase RV_data] #combine arrays
        phase_sorted = fastsortrows(phase_data_arr, [1]) #sort by phase
        for i in 2:length(time_data)
            sum[j] += (phase_sorted[i,2]-phase_sorted[i-1,2])*(phase_sorted[i,2]-phase_sorted[i-1,2])
        end
    end

    min_index=indmin(sum)
    P=period[min_index]

    # Calculate phase-sorted RV data
    sum = 0.0
    #Sort by phase given period
    phase = Array(Real,length(time_data))
    phase = mod(time_data,P)
    phase_data_arr = [phase RV_data] #combine arrays
    phase_sorted = fastsortrows(phase_data_arr, [1]) #sort by phase
    for i in 2:length(time_data)
        sum += (phase_sorted[i,2]-phase_sorted[i-1,2])*(phase_sorted[i,2]-phase_sorted[i-1,2])
    end

    # Add relative path to ExoJulia this doesn't work (at least in Linux, v.0.4.5)
    #push!(LOAD_PATH, "../../../ExoJulia/")

    # import
    # using ExoJulia

    # The prescribed method for importing a module doesn't work in Linux apparently. Ubuntu Julia v.0.4.5
    # So we must directly specify the include path.
    include("/home/linc/Documents/SCHOOL/598_exoplanets/ExoJulia/ExoJulia/Orbit/orbit.jl")

    # Initial values for curve_fit
    ecc = 0.0
    time_peri = 0.0 #presumably in days
    p = [ecc,time_peri,P] #period P from phase folding estimate
    global W = W_func(err_data) #Must be called W...

    # Run fitting routine for eccentricity & time of periastron
    fit = curve_fit(vrad_model,time_data,RV_data,err_data,[ecc,time_peri,P])
    println(fit.param)
end


function find_period(time_data,RV_data,showplots=false)
    ### Finds the period of an RV dataset by minimizing sq residuals
    # Showplots option allows outputting useful plots
    
    # Linear period grid using 100,000 datapoints spanning 1 sec to the window width
    period = linspace(1,time_data[end]-time_data[1],100000)
    
    # Initiate other necessary arrays
    sum = Array(Real,length(period))
    phase = Array(Real,length(time_data))

    # Loop over all period guesses
    for (j,P) in enumerate(period)
        sum[j] = 0.0
  
        # Sort by phase given period
        phase = mod(time_data,P)
        
        # Combine phase & RV data as columns
        phase_data_arr = [phase RV_data]
        
        # Sort by phase
        phase_sorted = fastsortrows(phase_data_arr,[1])
        
        # Sum sq residuals
        for i in 2:length(time_data)
            sum[j] += (phase_sorted[i,2]-phase_sorted[i-1,2])*(phase_sorted[i,2]-phase_sorted[i-1,2])
        end
        
    end
    
    min_index = indmin(sum)
    min_per = period[min_index]
    
    # For output of plots describing fit
    if(showplots)
        P = min_per
        sum0 = 0.0
        
        # Sort by phase, given period
        phase = mod(time_data,P)
        
        # Combine arrays
        phase_data_arr = [phase RV_data]
        
        # Sort by phase
        phase_sorted = sortrows(phase_data_arr, by=x->x[1])
        
        # Sum sq residuals
        for i in 2:length(time_data)
            sum0 += (phase_sorted[i,2]-phase_sorted[i-1,2])*(phase_sorted[i,2]-phase_sorted[i-1,2])
        end
        
        # Plot sorted data
        plot(phase_sorted[:,1],phase_sorted[:,2],".")

        # Plot residuals; this component is not currently very flexible
        npts = 2000
        plot(period[min_index-npts:min_index+npts],log10(sum[min_index-npts:min_index+npts]),".")

        else min_per
    end

end

### f value (value of periastron) ###

function f_val(time,ecc,period,time_peri)
    
    M = 2*pi/period*(time-time_peri)
    
    if(length(M) > 1)
        E = Array(Float64,length(M))
        f = Array(Float64,length(M))
        
        for i in 1:length(M)
            E[i] = Orbit.kepler_solve(M[i],ecc)
            f[i] = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E[i]/2))
        end
        return f 
    else 
        E = Orbit.kepler_solve(M,ecc)
        f = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2))
        return f
    end
end

### ORBITAL ELEMENTS ###

function W_func(err_data)
        
    W = zeros(Float64,length(err_data),length(err_data))
    #W[:,:] = 0.0
    
    for i in 1:length(err_data)
        W[i,i] = (err_data[i])^(-2)
    end
    W
end
        

function B_func(time_data,P,ecc,time_peri)   
    F = Array(Float64,4,length(time_data))

    for i in 1:length(time_data)
        F[1,i] = cos(f_val(time_data[i],ecc,P,time_peri))
        F[2,i] = sin(f_val(time_data[i],ecc,P,time_peri))
        F[3,i] = 1.0
        F[4,i] = time_data[i]-time_data[1]
    end

    B = RV_data'*W*F'*inv(F*W*F')
    return B   
end

### RADIAL VELOCITY FUNCTION ###

function vrad_model(x,p)
    
    # x = time data
        
    ecc = p[1]
    time_peri = p[2]
    P = p[3]
    B = B_func(x,P,ecc,time_peri)
    #B = {h,c,v0,d}
    h = B[1]
    c = B[2]
    v0 = B[3]
    #K = sqrt(h*h + c*c)
    #curlypi = atan(-c/h)
    #gam = v0 - K*ecc*cos(curlypi)
    f = f_val(x,ecc,P,time_peri)
    vrad = h*cos(f) + c*sin(f) + v0
    return vrad

end

function vrad_calc(x,p,B)
    
    # Given parameters, calculate RV
        
    ecc = p[1]
    time_peri = p[2]
    P = p[3]
    #B = {h,c,v0,d}
    h = B[1]
    c = B[2]
    v0 = B[3]
    #K = sqrt(h*h + c*c)
    #curlypi = atan(-c/h)
    #gam = v0 - K*ecc*cos(curlypi)
    f = f_val(x,ecc,P,time_peri)
    vrad = h*cos(f) + c*sin(f) + v0
    return vrad

end

# Traditional sort algorithms are inefficient in Julia :/
# Julia's sort algorithms are even worse.

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

time_rv()
