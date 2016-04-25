module Transit

include("../../hw2/Jake_Dave/utils.jl")
include("../../hw2/Jake_Dave/orbital_utils.jl")
include("../../hw2/Jake_Dave/rv.jl")
include("transit_utils.jl")
include("transit_loglike.jl")
include("transit_model.jl")

# Import necessary modules
using Optim

# Define global variables for the data
#global time = nothing
#global flux = nothing
#global err = nothing

############# Create Setter functions #############

function set_time(val)
   global time
   time = val
   println("Set Transit.time")
end

function set_flux(val)
   global flux
   flux = val
   println("Set Transit.flux")
end

function set_err(val)
   global err
   err = val
   println("Set Transit.err")
end

############# Include Fitting functions #############

function fit_transit(;MAX_ITER::Int=1000)
    
    # Fit the transit data for the observables
    optimum = optimize(transit_loglike_observables_optim, params, iterations=MAX_ITER, method=:gradient_descent)
    param_list =["dF : Delta Flux = Rp/Rs","tT : Transit Duration","tF : Flat Duration","P  : Period","tE : Time of first transit"]
    for i=1:length(param_list)
        println("$param_list[i] = $optimum.minimum[i]")
    end
    return optimum.minimum
end
    
function observable_to_physical(params::Vector)
    # Convert the fitted observables to physical system parameters (like rp/rstar, stellar density...)
    
end

end
