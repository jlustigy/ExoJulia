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

function plot_results_full(pbest)
    time_hires = collect(linspace(minimum(time),maximum(time),100000))

    scatter(time_hires,transit_model_observables(time_hires,pbest),color="red", alpha=0.1)

    scatter(time,mean_flux, color="blue", alpha=0.1)
    ylim(0.995,1.005)
    #xlim(0,2.5)
    
    show()
end 

function plot_results_folded(pbest)
    
    data_fold = copy(data)
    data_fold[:,1] = mod(time - time[1], pbest[4]);
    data_fold = convert(Array{Float64,2},fastsortrows(data_fold,[1]));

    scatter(data_fold[:,1],transit_model_observables(data_fold[:,1],pbest),color="red")
    scatter(data_fold[:,1],data_fold[:,2]/mean(data_fold[:,2]),color="blue",alpha=0.5)
    ylim(0.995,1.005)
    xlim(0.,1.5)
    
    show()
end