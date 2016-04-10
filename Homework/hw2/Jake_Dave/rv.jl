######################
#
#
#
#
######################

function v_semi_amp!(P::Float64, ecc::Float64, inc::Float64, Mp::Float64, Ms::Float64)
    (2.0pi * 6.67259e-8 / (P * (1.0 - e*e)^(3./2.)) * (Mp*sin(inc))^3 / (Mp + Ms)^2)^(1./3.)
end

function v_rad!(K::Float64, w::Float64, f::Float64, ecc::Float64, gamma::Float64)
    K * (cos(w+f) + e*sin(w)) + gamma
end

function v_rad_lin!(h::Float64, f::Float64, c::Float64, v0::Float64)
    h*cos(f) + c*sin(f) + v0
end

function agol_periodogram(numbers::Array{Float64,2}, periods::Array{Float64,1})

    # Create array to store phase folded residuals
    dv = zeros(length(periods));

    # Append column to hold phase folded times for sorting
    numbers = [numbers zeros(length(numbers[:,1]))];

    # Loop over periods
    for i=1:length(periods)
        # Phase Phold!
        numbers[:,4] = mod(numbers[:,1] - numbers[1,1], periods[i])
        # Sort by phase
        numbers = fastsortrows(numbers, [4])
        # Loop over RV data, summing squared differences between adjacent points in phase space
        for j=2:length(numbers[:,2])
            dv[i] += (numbers[j,2] - numbers[j-1,2])^2
        end
        # Sort by observational time to reset grid
        numbers = fastsortrows(numbers, [1])
    end

    # Return the period where the summed squared differences are minimal
    return periods[argmin(dv)]
end
