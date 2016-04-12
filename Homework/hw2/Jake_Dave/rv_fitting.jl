include("utils.jl")
include("orbital_utils.jl")
include("rv.jl")

function WH_W(err::Array{Float64,1})
    # Compose Wright & Howard  matrix
    W = zeros(Float64,length(err), length(err))
    for i=1:length(err)
        W[i,i] = 1./err[i]^2
    end
    return W
end

function WH_F(P::Float64, ecc::Float64, t::Array{Float64, 1}, tp::Float64; Nplanets::Int=1)
    # Compute Wright & Howard F matrix

    # Allocate matrix
    F = zeros(Float64, (2*Nplanets+2), length(t))

    # Fill matrix
    for i=1:Nplanets
        for j=1:length(t)
            f = f_from_t(P, ecc, t[j], tp)
            F[2*i - 1,j] = cos(f)
            F[2*i, j] = sin(f)
        end
    end

    #
    for j=1:length(t)
        F[end-1, j] = 1.0
        F[end, j] = t[j] - t[1]
    end
    return F
end

function WH_eps(F::Array{Float64, 2}, W::Array{Float64, 2})
    # Compute Wright & Howard epsilon matrix
    return inv(F * W * (F'))
end

function WH_Beta(RV::Array{Float64, 2}, F::Array{Float64, 2}, W::Array{Float64, 2})
    # Compute Wright & Howard Beta vector
    return RV * W * (F') * WH_eps(F, W)
end

function rv_forward(P::Float64, ecc::Float64, tp::Float64, t::Array{Float64,1}, rv::Array{Float64,1}, err::Array{Float64,1}; Nplanets::Int=1)
    # Calculates model RV given linear params: [h, c, v0, d, tp]

    # Allocate
    rv_mod = zeros(Float64, Nplanets, length(t))

    # Calculate Wright & Howard Beta Vector, B = [h_i, c_i, ..., h_n, c_n, v0, d] for n planets
    B = WH_Beta(rv', WH_F(P, ecc, t, tp), WH_W(err))

    # Loop over planets and observations, calculating model rv points
    for i=1:Nplanets
        for j=1:length(t)
            f = f_from_t(P, ecc, t[j], tp)
            rv_mod[i,j] = v_rad_lin(B[i,1], f, B[i,2], B[i,end-1], t[j], t[1], B[i,end])
        end
    end

    return rv_mod, B
end

function rv_loglike(rho)
    #rho = [period, ecc, tp]

    # hard bounds
    if rho[1] < 0.0
        return Inf
    end
    if rho[2] >= 1.0
        return Inf
    end
    if rho[2] < 0.0
        return Inf
    end

    # call forward model
    model, B = rv_forward(rho[1], rho[2], rho[3], time, rv, err);

    # Chi^2
    return -loglike(rv, model, err);
end

function rv_curve_forward(t, p::Vector)

    # call forward model
    model, B = rv_forward(p[1], p[2], p[3], t, rv, err);

    return reshape(model, length(model))

end

function solve_rv(data::Array{Float64, 2}; p0=[nothing, nothing, nothing], alg::AbstractString="cf")
    #= This function takes a 2D (Nx3) array of RV data, where N is the number of data points with columns
    of time, RV, and error, and returns the best fitting period, eccentricity, and time of periastron.

    p0 = [period, ecc, tp]
    =#

    # Unpack RV data (make global?)
    time = data[:,1];
    rv = data[:,2];
    err = data[:,end];

    # Set initial parameters if not specified
    p = [0.0, 0.0, 0.0]
    if p0[1] == nothing
        # Use Agol Periodogram for initial period guess
        periods = collect(linspace(minimum(time[2:end] - time[1:end-1]), time[end]-time[1], 10000))
        p[1] = agol_periodogram(numbers, periods)
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

    return pbest
end
