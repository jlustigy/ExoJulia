using Optim
using LsqFit
using ExoJulia

# Functions included
include("utils.jl")
include("orbital_utils.jl")
include("rv.jl")
include("rv_fitting.jl")
include("/Users/Jake/Projects/ExoJulia/ExoJulia/Orbit/orbit.jl")

function fit_rv(data::Array{Float64, 2}; p0=[nothing, nothing, nothing], alg::AbstractString="cf")
    #= This function takes a 2D (Nx3) array of RV data, where N is the number of data points with columns
    of time, RV, and error, and returns the best fitting period, eccentricity, and time of periastron.

    p0 = [period, ecc, tp]
    =#


    # Unpack RV data (make global?)
    global time = data[:,1];
    global rv = data[:,2];
    global err = data[:,end];

    # Set initial parameters if not specified
    p = [0.0, 0.0, 0.0]
    if p0[1] == nothing
        # Use Agol Periodogram for initial period guess
        periods = collect(linspace(minimum(time[2:end] - time[1:end-1]), time[end]-time[1], 10000))
        p[1] = agol_periodogram(data, periods)
    else
        p[1] = p0[1]
    end
    if p0[2] == nothing
        # Use random ecc [0,1) for initial guess
        p[2] = rand()
    else
        p[2] = p0[2]
    end
    if p0[3] == nothing
        # Use random time drawn from observed grid
        p[3] = rand(time)
    else
        p[3] = p0[3]
    end

    # Run solver using either curve_fit() or optimize()
    if alg == "cf"
        # Use curve_fit() to fit
        fit = curve_fit(rv_curve_forward, time, rv, 1.0./err.^2, p);
        pbest = fit.param
    elseif alg == "opt"
        # Use optimize to fit
        optimum = optimize(rv_loglike, p, autodiff=true)
        pbest = optimum.minimum
    else
        print("Choose alg = 'cf' or 'opt'")
        return
    end
    print(pbest)

    return pbest
end
