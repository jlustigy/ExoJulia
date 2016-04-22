# Load in our functions
push!(LOAD_PATH,"../../../ExoJulia/")
push!(LOAD_PATH,".")

include("../../hw2/Jake_Dave/utils.jl")
include("../../hw2/Jake_Dave/orbital_utils.jl")
include("../../hw2/Jake_Dave/rv.jl")
include("transit_utils.jl")
include("transit_loglike.jl")
include("transit_model.jl")

using ExoJulia
using Optim

# Load in data
data = readdlm("mystery_planet2.txt");
time = data[:,1];
flux = data[:,2];
mean_flux = flux/median(flux)
err = data[:,3];

function transit_fit_speed_test()
    # Estimate parameters
    #Returns est_params = [ dF, tT, tF, P, tE ]
    est_params = guess_params()
    best_period = est_params[4];
    best_tE = est_params[5];

    # Fit using optim
    # params [ dF, tT, tF, P, tE ]
    params = est_params
    optimum = optimize(transit_loglike_observables_optim, params, iterations=5000, method=:bfgs)

    # Returns [P, dF, b_est, tT, rho_est]
    phys = observable_to_physical(best_optim)
    print("Best Fit Period: $(phys[1]) days.\n")
    print("Best Fit Depth: $(phys[2]).\n")
    print("Best Fit Impact Parameter: $(phys[3]).\n")
    print("Best Fit Transit Time: $(phys[4]) days.\n")
    print("Best Fit Stellar Density: $(phys[5]) solar density.\n")
    return phys
end

transit_fit_speed_test()

#@stest transit_fit_speed_test()